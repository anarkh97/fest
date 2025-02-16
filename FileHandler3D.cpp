#include<FileHandler3D.h>
#include<climits>
#include<fstream>
#include<array>
#include<deque>
#include<memory>
#include<cstring>
#include<sstream>
#include<iterator>
#include<algorithm>

using std::string;
using std::vector;
using std::pair;
using std::sort;

extern int verbose;

//------------------------------------------------------------

FileHandler3D::FileHandler3D(MetaInputData& meta_, MPI_Comm& comm_)
                : iod_meta(meta_), comm(comm_)
{

  if(strcmp(iod_meta.metafile, "") == 0) {
    print_error("*** Error: Metafile with paths to existing FSI simulations and "
                "their design parameters was not specified. Aborting.\n");
    exit_mpi();
  }

  ReadMetaFile();

}

//------------------------------------------------------------

string& FileHandler3D::GetMeshFileForProxim(int id)
{
  int num_points = iod_meta.numPoints;
  assert(id >= 0 and id < num_points);
  assert(!surface_files.empty());

  return surface_files[id];
}

//------------------------------------------------------------

string& FileHandler3D::GetSolnFileForProxim(int id)
{
  int num_points = iod_meta.numPoints;
  assert(id >= 0 and id < num_points);
  assert(!solution_files.empty());

  return solution_files[id];
}

//------------------------------------------------------------

vector<double>&
FileHandler3D::GetParametersForProxim(int id)
{
  int num_points = iod_meta.numPoints;
  assert(id >= 0 and id < num_points);
  assert(!associates.empty());

  return associates[id];
}

//------------------------------------------------------------

void
FileHandler3D::ReadMetaFile()
{

  // filename
  string filename = /*iod_meta.prefix +*/ iod_meta.metafile;
  
  std::fstream input;
  input.open(filename.c_str(), std::fstream::in);
  if(!input.is_open()) {
    print_error("*** Error: Cannot open file %s.\n", filename.c_str());
    exit_mpi();
  } else
    print("- Reading user-specified metafile: %s.\n", filename.c_str());

  string line, word;

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
  vector<string>    surf_files;
  vector<string>    soln_files;
  vector<vector<double>> parameters;

  int target_index = 0;
  int index = 0;
  while(getline(input, line)) {
    std::istringstream is(line);
    string directory, mesh, soln;

    is >> directory >> mesh >> soln;

    // check if this is the current parameter
    if(!(directory.compare(0, 4, "Target", 0, 4) and
         directory.compare(0, 4, "TARGET", 0, 4) and
	 directory.compare(0, 4, "target", 0, 4) and
         mesh.compare(0, 4, "Target", 0, 4) and
         mesh.compare(0, 4, "TARGET", 0, 4) and
	 mesh.compare(0, 4, "target", 0, 4) and
         soln.compare(0, 4, "Target", 0, 4) and
         soln.compare(0, 4, "TARGET", 0, 4) and
	 soln.compare(0, 4, "target", 0, 4))) {
      target_index = index;
    }

    // read parameters
    double value;
    vector<double> data;
    while(is >> value) data.push_back(value);

    // contains target as well
    parameters.push_back(data);
    surf_files.push_back(directory + mesh);
    soln_files.push_back(directory + soln);

    index++;
  }

  if(!target_index) {
    print_error("*** Error: Did not find the target's parameters.\n");
    exit_mpi();
  }

  if(iod_meta.numPoints != (int)parameters.size()-1) {
    print_warning("- Warning: The number of interpolation points specified in "
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
      dist += (parameters[i][d] - target[d])*(parameters[i][d] - target[d]);
    
    dist2target.push_back(std::make_pair(dist, i));

  }

  sort(dist2target.begin(), dist2target.end());

  // store the other points
  associates.resize(iod_meta.numPoints, vector<double>(dim, 0.0));
  surface_files.resize(iod_meta.numPoints, "");
  solution_files.resize(iod_meta.numPoints, "");

  index = 0;
  std::vector<double> distances(iod_meta.numPoints, -1);
  for(int i=0; i<iod_meta.numPoints+1; ++i) {

    if(dist2target[i].first == 0) continue; // skip target

    associates[index]     = parameters[dist2target[i].second];
    surface_files[index]  = surf_files[dist2target[i].second];
    solution_files[index] = soln_files[dist2target[i].second];
    distances[index]      = std::sqrt(dist2target[i].first); // mainly for log purpose.
    index++;

  }

  if(verbose>0) {
    print("- Interpolating from parameters:\n");
    for(int i=0; i<iod_meta.numPoints; ++i) {
      print("  o Parameter %d:", i);
      for(int j=0; j<(int)associates[i].size(); ++j)
        print("  %e", associates[i][j]);
      print(" (r = %e)\n", distances[i]);
    }
    print("\n");
  }

}

//------------------------------------------------------------

//! Code copied from EmbeddedBoundaryOperator::ReadMeshFile
void
FileHandler3D::ReadMeshFile(string &filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{

  string fname(filename);
  auto loc = fname.find_last_of(".");
  if(loc>=fname.size()-1) {//assume the default format (top) if file extension not detected
    
    if(verbose>1)
      print_warning("- Warning: Mesh file does not have a file extension. Assuming file is "
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

//------------------------------------------------------------

//! Code copied from EmbeddedBoundaryOperator::ReadMeshFileInTopFormat
void
FileHandler3D::ReadMeshFileInTopFormat(string &filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{

  // read data from the surface input file.
  FILE *topFile;
  topFile = fopen(filename.c_str(), "r");
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
    print_warning("- Warning: The node indices of an embedded surface may have a gap: "
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
    print_warning("- Warning: The element indices of an embedded surface may have a gap: "
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

//! Solution files will be large and in ascii. Maybe later we can
//! implement a vtr/vtu reader as well. However, since these files
//! are large we only read them on processor#0 and broadcast the 
//! data.
void
FileHandler3D::ReadSolutionFile(string &filename, SolutionData3D &S)
{

  // solution files are assummed to be large.
  // Consequently, we read them once and only on rank 0.
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  int soln_rows        = 0; // these are rows/nodes in the solution file.
  int soln_total_size  = 0; // this is the total size of the data when flattened.
  vector<double> soln;

  if(mpi_rank == 0) {

    // read the solution file.
    ReadSolutionFileSingleProcessor(filename, S);

    // update the sizes
    soln_rows       = S.GetRows();
    soln_total_size = S.GetSize();

    soln.resize(soln_total_size);
    // this process copies the data in
    S.Flatten(soln);

  }

  // broadcast sizes to all ranks
  MPI_Bcast(      &soln_rows, 1, MPI_INT, 0, comm);
  MPI_Bcast(&soln_total_size, 1, MPI_INT, 0, comm);

  // prepare to recieve read data on ranks other than root
  if(mpi_rank != 0)
    soln.resize(soln_total_size);

  // broadcast the data vector to all ranks
  MPI_Bcast(soln.data(), soln_total_size, MPI_DOUBLE, 0, comm);

  // rebuild SolutionData3D on ranks other than root
  if(mpi_rank != 0)
    S.Rebuild(soln, soln_rows);

  /*
  if(mpi_rank == 0) {
    auto &stored = S.GetSolutionAtTime(0.0);
    for(auto &vec : stored) {
      fprintf(stdout, "%e  %e  %e\n", vec[0], vec[1], vec[2]);
    }
    exit_mpi();
  }
  */

  MPI_Barrier(comm);

}

//------------------------------------------------------------

void
FileHandler3D::ReadSolutionFileSingleProcessor(string &filename, SolutionData3D &S)
{

  std::fstream input(filename.c_str(), std::fstream::in);
  if(!input.is_open()) {
    fprintf(stderr, "\033[0;31m*** Error: Cannot open file %s.\n\033[0m", filename.c_str());
    exit(-1);
  }

  int num_nodes;

  // start reading
  std::string line, word;
  // Line#1 --- ignored
  input.ignore(512, '\n');
  // Line#2 --- number nodes
  getline(input, line);
  std::istringstream iss(line);
  iss >> num_nodes;

  while(getline(input, line)) {

    double time;
    vector<Vec3D> snapshot(num_nodes, Vec3D(0.0));

    // read line
    std::istringstream is(line);

    // read time
    if(!(is >> time)) {
      fprintf(stderr, "\033[0;31m*** Error: Could not read a snapshot from the "
                      "solution file (%s).\n\033[0m", filename.c_str());
      exit(-1);
    }

    // read data
    bool ddone = false;
    for(int i=0; i<num_nodes; ++i) {

      // read the next line
      getline(input, line);
      is.clear();
      is.str(line);

      bool done = false;
      double x, y, z;
      if(!(is >> x >> y >> z)) {
        done = true;
        break;
      }

      if(done) {
        ddone = true;	      
	break;
      }

      snapshot[i] = Vec3D(x, y, z);

    }

    if(ddone) break;

    // grows the SolutionData3D object.
    // Note: snapshot is "moved" to an internal container.
    // It will no longer be available after this.
    S.Insert(time, std::move(snapshot));

    /*
    auto &stored = S.GetSolutionAtTime(time); 
    for(auto &vec : stored)
      fprintf(stdout, "%e  %e  %e\n", vec[0], vec[1], vec[2]);
    exit(-1);
    */

  }

  assert(!S.Empty());
  input.close();

}

//------------------------------------------------------------

