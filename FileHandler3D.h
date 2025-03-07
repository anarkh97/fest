#ifndef _FILE_HANDLER_3D_
#define _FILE_HANDLER_3D_

#include<IoData.h>
#include<Vector3D.h>
#include<SolutionData3D.h>

class MetaPoint {

  std::string directory;
  std::string surface_file;
  std::string solution_file;
  std::vector<double> x;

public:

  MetaPoint() : directory(""), surface_file(""), solution_file(""), x({}) { }; 
  ~MetaPoint() { };

  // setters
  void SetPointDirectory(std::string &dir) { directory = dir; };
  void SetPointSurfaceFile(std::string &surf) { surface_file = surf; };
  void SetPointSolutionFile(std::string &soln) { solution_file = soln; };
  void SetPointParameters(double *data, int size) {
    if(!x.empty())
      x.resize(size, 0.0); // clear previous values
    else
      x.assign(size, 0.0);
  
    for(int i=0; i<size; ++i)
      x[i] = data[i];
  };

  // getters
  const std::string& GetPointDirectory() const { return directory; };
  const std::string& GetPointSurfaceFile() const { return surface_file; };
  const std::string& GetPointSolutionFile() const { return solution_file; };
  const std::vector<double>& GetPointParameters() const { return x; };

  int GetDim() const { return x.size(); };

  // operators
  double GetDistance(const MetaPoint& other);

};

inline double MetaPoint::GetDistance(const MetaPoint &other)
{

  assert((int)x.size() == other.GetDim());

  const std::vector<double> &other_x = other.GetPointParameters();
  
  double dist = 0;
  for(int i=0; i<(int)x.size(); ++i) 
    dist += std::pow(x[i] - other_x[i], 2);

  return std::sqrt(dist);

}

class FileHandler3D {

  MetaInputData& iod_meta;
  MPI_Comm& comm;

  //! meta-level parameters provided by user.
  enum Var {TARGET=0, NEIGHBOR};
  std::map<Var,std::vector<MetaPoint>> points;

  //! maps the order of fields specified in the metafile.
  //! mainly used for book-keeping. We expect a certain 
  //! order, which is checked for while reading the metafile.
  enum Fields {TYPE=0, DIRECTORY, MESH_FILE, SOLUTION_FILE, 
               PARAMETER_START, PARAMETER_END};
  std::map<Fields, int> field2col;

public:

  FileHandler3D(MetaInputData& meta_, MPI_Comm& comm_);
  ~FileHandler3D() { }

  // Methods for reading surface and solution files.
  void ReadMetaFile();
  void ReadMeshFile(std::string &filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);
  void ReadSolutionFile(std::string &filename, SolutionData3D &S);

  std::string GetMeshFileForTarget();
  std::string GetMeshFileForNeighbor(int id);

  std::string GetSolnFileForTarget();
  std::string GetSolnFileForNeighbor(int id);

  std::vector<double> GetParametersForTarget();
  std::vector<double> GetParametersForNeighbor(int i);
  std::vector<std::vector<double>> GetParametersForAllNeighbors(); 

private:

  void ReadMeshFileInTopFormat(std::string& filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);
  void ReadSolutionFileSingleProcessor(std::string &filename, SolutionData3D &S);

};

#endif
