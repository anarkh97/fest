#include<InterpolationOperator.h>
#include<climits>
#include<fstream>
#include<deque>
#include<sstream>
#include<algorithm>

#include<MathTools/rbf_interp.hpp>
#include<MathTools/affine_transformation.h>

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
  
  int num_points = iod_meta.numPoints;

  // read mesh files for each point.
  vector<TriangulatedSurface> proxi_surfaces(num_points);
  for(int i=0; i<num_points; ++i) {
    std::string filename = proxi_surface_files[i];
    ReadMeshFile(filename.c_str(), proxi_surfaces[i].X, proxi_surfaces[i].elems);

    proxi_surfaces[i].BuildConnectivities();
    proxi_surfaces[i].CalculateNormalsAndAreas();

    // scale the co-ordinates to [0, 1].
    vector<Vec3D> UV(proxi_surfaces[i].X.size(), Vec3D(0.0));
    MathTools::affine_transformation(proxi_surfaces[i].X.size(), 3, 
                                     (double*)proxi_surfaces[i].X.data(), (double*)UV.data());

    // update surface co-ordinates (X0 not used by TriangulatedSurface for normal and area
    // calculations, should be safe to update co-ordinates).
    proxi_surfaces[i].X  = UV;
    proxi_surfaces[i].X0 = proxi_surfaces[i].X;
  }



}

//------------------------------------------------------------------------------------------------

//! Code copied from EmbeddedBoundaryOperator::ReadMeshFile
void
InterpolationOperator::ReadMeshFile(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{
  string fname(filename);
  auto loc = fname.find_last_of(".");
  if(loc>=fname.size()-1) {//assume the default format (top) if file extension not detected
    
    if(verbose>1)
      print_warning("*** Warning: Mesh file does not have a file extension. Assuming file is "
                    "written in top format.\n");

    ReadMeshFileInTopFormat(filename, Xs, Es); 
    return;
  }

  string ext = fname.substr(loc+1);
  if(ext == "top" || ext == "Top" || ext == "TOP")
    ReadMeshFileInTopFormat(filename, Xs, Es); 
  else {
    print_error("*** Error: Currently FEST can only read surfaces in top format.\n");
    exit_mpi();
  }
}

//------------------------------------------------------------------------------------------------

//! Code copied from EmbeddedBoundaryOperator::ReadMeshFileInTopFormat
void
InterpolationOperator::ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{

  // read data from the surface input file.
  FILE *topFile;
  topFile = fopen(filename, "r");
  if(topFile == NULL) {
    print_error("*** Error: Cannot open embedded surface mesh file (%s).\n", filename);
    exit_mpi();
  }
 
  int MAXLINE = 500;
  char line[MAXLINE], key1[MAXLINE], key2[MAXLINE]; //, copyForType[MAXLINE];

  int num0 = 0;
  int num1 = 0;
  double x1, x2, x3;
  int node1, node2, node3;
  int type_reading = 0; //1 means reading node set; 2 means reading element set
  std::deque<std::pair<int, Vec3D>> nodeList;
  std::deque<std::array<int, 4>> elemList; // element ID + three node IDs
  int maxNode = 0, maxElem = 0;
  bool found_nodes = false;
  bool found_elems = false;


  // --------------------
  // Read the file
  // --------------------
  while(fgets(line, MAXLINE, topFile) != 0) {

    sscanf(line, "%s", key1);
    string key1_string(key1);

    if(key1[0] == '#') {
      //Do nothing. This is user's comment
    }
    else if(same_strings_insensitive(key1_string,"Nodes")){
      if(found_nodes) {//already found nodes... This is a conflict
        print_error("*** Error: Found multiple sets of nodes (keyword 'Nodes') in %s.\n", filename);
        exit_mpi();
      }
      sscanf(line, "%*s %s", key2);
      type_reading = 1;
      found_nodes = true;
    }
    else if(same_strings_insensitive(key1_string, "Elements")) {

      if(found_elems) {//already found elements... This is a conflict
        print_error("*** Error: Found multiple sets of elements (keyword 'Elements') in %s.\n", filename);
        exit_mpi();
      }
      type_reading = 2;
      found_elems = true;

    }
    else if(type_reading == 1) { //reading a node (following "Nodes Blabla")
      int count = sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
      if(count != 4) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a node.\n", line, filename);
        exit_mpi();
      }
      if(num1 < 1) {
        print_error("*** Error: detected a node with index %d in embedded surface file %s.\n", num1, filename);
        exit_mpi();
      }
      if(num1 > maxNode)
        maxNode = num1;

      nodeList.push_back({num1, {x1, x2, x3}});
    }
    else if(type_reading == 2) { // we are reading an element --- HAS TO BE A TRIANGLE!
      int count = sscanf(line, "%d %d %d %d %d", &num0, &num1, &node1, &node2, &node3);
      if(count != 5) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a triangular element.\n", line, filename);
        exit_mpi();
      }
      if(num0 < 1) {
        print_error("*** Error: detected an element with index %d in embedded surface file %s.\n", num0, filename);
        exit_mpi();
      }
      if(num0 > maxElem)
        maxElem = num0;

      elemList.push_back({num0, node1, node2, node3});
    }
    else { // found something I cannot understand...
      print_error("*** Error: Unable to interpret line %s (in %s).\n", line, filename);
      exit_mpi();
    }

  }

  fclose(topFile);

  if(!found_nodes) {
    print_error("*** Error: Unable to find node set in %s.\n", filename);
    exit_mpi();
  }
  if(!found_elems) {
    print_error("*** Error: Unable to find element set in %s.\n", filename);
    exit_mpi();
  }

  // ----------------------------
  // Now, check and store nodes
  // ----------------------------
  int nNodes = nodeList.size();
  map<int,int> old2new;
  Xs.resize(nNodes);
  int id(-1);
  if(nNodes != maxNode) { // need to renumber nodes, i.e. create "old2new"
    print_warning("Warning: The node indices of an embedded surface may have a gap: "
                  "max index = %d, number of nodes = %d. Renumbering nodes. (%s)\n",
                  maxNode, nNodes, filename);
//    assert(nNodes < maxNode);

    int current_id = 0; 
    vector<bool> nodecheck(maxNode+1, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first;
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[current_id] = it1->second; 
      old2new[id] = current_id;
      current_id++;
    }
    assert(current_id==(int)Xs.size());
  } 
  else { //in good shape
    vector<bool> nodecheck(nNodes, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first - 1; 
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id+1, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[it1->first - 1] = it1->second;
    }
  }


  // ------------------------------
  // check nodes used by elements
  // ------------------------------
  for(auto it = elemList.begin(); it != elemList.end(); it++) {

    id = (*it)[0];
    node1 = (*it)[1];
    node2 = (*it)[2];
    node3 = (*it)[3];
      
    if(old2new.empty()) {//node set is original order

      if(node1<=0 || node1 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      if(node2<=0 || node2 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

      if(node3<=0 || node3 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node3, id, filename);
        exit_mpi();
      }
    }
    else {// nodes are renumbered

      auto p1 = old2new.find(node1);
      if(p1 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      auto p2 = old2new.find(node2);
      if(p2 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

      auto p3 = old2new.find(node3);
      if(p3 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node3, id, filename);
        exit_mpi();
      }
    }
  }


  // ----------------------------
  // check and store elements
  // ----------------------------
  int nElems = elemList.size();
  Es.resize(nElems);
  if(nElems != maxElem) { // need to renumber elements.
    print_warning("Warning: The element indices of an embedded surface may have a gap: "
                  "max index = %d, number of elements = %d. Renumbering elements. (%s)\n",
                  maxElem, nElems, filename);
//    assert(nElems < maxElem);
    
    int current_id = 0; 
    vector<bool> elemcheck(maxElem+1, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0];
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      node3 = (*it)[3];
      
      if(old2new.empty()) //node set is original order
        Es[current_id] = Int3(node1-1, node2-1, node3-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        auto p3 = old2new.find(node3);
        Es[current_id] = Int3(p1->second, p2->second, p3->second);
      }      
      current_id++;
    }
  } 
  else { //element numbers in good shape

    vector<bool> elemcheck(nElems, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0] - 1;
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      node3 = (*it)[3];
      
      if(old2new.empty()) //node set is original order
        Es[id] = Int3(node1-1, node2-1, node3-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        auto p3 = old2new.find(node3);
        Es[id] = Int3(p1->second, p2->second, p3->second);
      }
    }
  }

}

//------------------------------------------------------------

void
InterpolationOperator::ComputeApproximateForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                                std::vector<Vec3D> *force_over_area, double t)
{
  //
}

//------------------------------------------------------------

