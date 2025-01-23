#include<InterpolationOperator.h>
#include<climits>
#include<fstream>
#include<sstream>
#include<algorithm>

#include<MathTools/rbf_interp.hpp>

using std::vector;
using std::pair;
using std::sort;

extern int verbose;

//------------------------------------------------------------

InterpolationOperator::InterpolationOperator(IoData &iod_, MPI_Comm &comm_)
                     : iod_meta(iod_.interp_driver.meta_input), 
		       iod_spatial(iod_.interp_driver.spatial_interp),
		       comm(comm_)
{
  if(strcmp(iod_meta.metafile, "") == 0) {
    print_error("*** Error: Metafile with paths to existing FSI simulations and "
                "their design parameters was not specified. Aborting.\n");
    exit_mpi();
  }

  ReadMetaFile();

}

//------------------------------------------------------------

void InterpolationOperator::Destroy()
{
  //
}

//------------------------------------------------------------

void 
InterpolationOperator::ReadMetaFile()
{

  // filename
  std::string filename = /*iod_meta.prefix +*/ iod_meta.metafile;
  
  std::fstream input;
  input.open(filename.c_str(), std::fstream::in);
  if(!input.is_open()) {
    print_error("*** Error: Cannot open file %s.\n", filename.c_str());
    exit_mpi();
  } else
    print("- Reading user-specified metafile: %s.\n", filename.c_str());

  std::string line, word;

  //Line #1 --- user's comment (skip)
  input.ignore(512, '\n');

  //Line #2 --- fields in the data
  input.ignore(2, ' '); //This line should start with ##
  getline(input, line);

  std::istringstream iss(line);
  bool found_parameters = false;
  int column = 0;
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0, 4, "Directory", 0, 4) and
         word.compare(0, 4, "DIRECTORY", 0, 4) and
	 word.compare(0, 4, "directory", 0, 4))) {
      field2col[DIRECTORY] = column;
      column += 1;
    }
    else if(!(word.compare(0, 4, "Mesh", 0, 4) and
              word.compare(0, 4, "MESH", 0, 4) and
	      word.compare(0, 4, "mesh", 0, 4))) {
      field2col[MESH_FILE] = column;
      column += 1;
    }
    else if(!(word.compare(0, 4, "Solution", 0, 4) and
              word.compare(0, 4, "SOLUTION", 0, 4) and
	      word.compare(0, 4, "solution", 0, 4))) {
      field2col[SOLUTION_FILE] = column;
      column += 1;
    }
    else if(!(word.compare(0, 4, "Parameters", 0, 4) and
              word.compare(0, 4, "PARAMETERS", 0, 4) and
	      word.compare(0, 4, "parameters", 0, 4))) {
      found_parameters = true;
      field2col[PARAMETER_START] = column;
      column += 1;
    }
    else if(found_parameters) {
      print_error("*** Error: Parameters were found before the end of "
                  "line. Move them to the end, after specifying other fields "
		  "in the metafile.\n");
      exit_mpi();
    }
    else {
      print_error("*** Error: I do not understand the word '%s' in the metafile.\n", word.c_str());
      exit_mpi(); 
    }
  }

  if(field2col[DIRECTORY] != 0 or field2col[MESH_FILE] != 1 or field2col[SOLUTION_FILE] != 2) {
    print_error("*** Error: Metafile field not in expected order which is,\n"
                "Directory  Mesh  Solution  Parameters\n");
    exit_mpi();
  }

  //Line #3 (and onwards) - directory, meshfile, solutionfile, parameter#1, ..., parameter#N
  vector<std::string> surface_files;
  vector<std::string> solution_files;
  vector<vector<double>> parameters;

  int target_index = 0;
  for(int i=0; i<INT_MAX; ++i) {
    getline(input, line);
    std::istringstream is(line);
    std::string directory, mesh, solution;

    is >> directory >> mesh >> solution;

    // check if this is the current parameter
    if(!(directory.compare(0, 4, "Target", 0, 4) and
         directory.compare(0, 4, "TARGET", 0, 4) and
	 directory.compare(0, 4, "target", 0, 4) and
         mesh.compare(0, 4, "Target", 0, 4) and
         mesh.compare(0, 4, "TARGET", 0, 4) and
	 mesh.compare(0, 4, "target", 0, 4) and
         solution.compare(0, 4, "Target", 0, 4) and
         solution.compare(0, 4, "TARGET", 0, 4) and
	 solution.compare(0, 4, "target", 0, 4))) {
      target_index = i;
    }

    // read parameters
    double value;
    vector<double> data;
    while(is >> value) data.push_back(value);

    if(is.fail()) break;

    parameters.push_back(data);
    surface_files.push_back(directory + mesh);
    solution_files.push_back(directory + solution);
  }

  if(!target_index) {
    print_error("*** Error: Did not find the target's parameters.\n");
    exit_mpi();
  }

  if(iod_meta.numPoints != (int)parameters.size()-1) {
    print_warning("*** Warning: The number of interpolation points specified in "
                  "the input file is %d. However, only %d points are provided in "
		  "the metafile. Using %d points for interpolation instead.\n",
		  iod_meta.numPoints, parameters.size()-1, parameters.size()-1);
    iod_meta.numPoints = parameters.size()-1;
  }

  // check the size of each parameter
  target = parameters[target_index];
  int dim = target.size();
  for(int i=0; i<iod_meta.numPoints+1; ++i)
    if(dim != (int)parameters[i].size()) {
      print_error("*** Error: The dimension of parameter %d does not match the rest.\n",
                  i+1);
      exit_mpi();
    }

  // sort the parameters based on distance from target parameter.
  vector<pair<double, int>> dist2target;
  for(int i=0; i<iod_meta.numPoints+1; ++i) {
 
    double dist = 0;
    for(int d=0; d<dim; ++d) 
      dist += parameters[i][d]*parameters[i][d] - target[d]*target[d];
    
    dist2target.push_back(std::make_pair(dist, i));

  }

  sort(dist2target.begin(), dist2target.end());

  // store the other points
  proximals.resize(iod_meta.numPoints, vector<double>(dim, 0.0));
  proxi_surface_files.resize(iod_meta.numPoints, "");
  proxi_solution_files.resize(iod_meta.numPoints, "");

  int index = 0;
  for(int i=0; i<iod_meta.numPoints+1; ++i) {

    if(i == target_index) continue;

    proximals[index]            = parameters[dist2target[i].second];
    proxi_surface_files[index]  = surface_files[dist2target[i].second];
    proxi_solution_files[index] = solution_files[dist2target[i].second];

    index++;

  }

  if(verbose>1) {
    print("- Interpolating from parameters:\n");
    for(int i=0; i<iod_meta.numPoints; ++i) {
      print("  o Parameter %d:", i);
      for(int j=0; j<(int)proximals.size(); ++j)
        print("  %e", proximals[i][j]);
      print("\n");
    }
    print("\n");
  }

}

//------------------------------------------------------------

void
InterpolationOperator::BuildSurfacesToSurfaceMap(TriangulatedSurface &surface)
{
  //
}

//------------------------------------------------------------

void
InterpolationOperator::ComputeApproximateForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                                std::vector<Vec3D> *force_over_area, double t)
{
  //
}

//------------------------------------------------------------

