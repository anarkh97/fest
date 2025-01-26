#ifndef _FILE_HANDLER_3D_
#define _FILE_HANDLER_3D_

#include<IoData.h>
#include<Vector3D.h>
#include<SolutionData3D.h>

class FileHandler3D {

  MetaInputData& iod_meta;
  MPI_Comm& comm;

  //! meta-level parameters provided by user.
  std::vector<double> target;
  std::vector<std::vector<double>> associates;

  //! paths of stored (or existing) surfaces
  std::vector<std::string> surface_files;
  //! paths of stored (or existing) solutions
  std::vector<std::string> solution_files;

  //! maps the order of fields specified in the metafile.
  //! mainly used for book-keeping. We expect a certain 
  //! order, which is checked for while reading the metafile.
  enum Fields {DIRECTORY=0, MESH_FILE, SOLUTION_FILE, PARAMETER_START,
               PARAMETER_END};
  std::map<Fields, int> field2col;

public:

  FileHandler3D(MetaInputData& meta_, MPI_Comm& comm_);
  ~FileHandler3D() { }

  // Methods for reading surface and solution files.
  void ReadMetaFile();
  void ReadMeshFile(std::string &filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);
  void ReadSolutionFile(std::string &filename, SolutionData3D &S);

  std::string& GetMeshFileForProxim(int id);
  std::string& GetSolnFileForProxim(int id);

  std::vector<double>& GetParametersForTarget() { return target; };
  std::vector<double>& GetParametersForProxim(int i);
  std::vector<std::vector<double>>& GetParametersForAllProxim() { return associates; };

private:

  void ReadMeshFileInTopFormat(std::string& filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);
  void ReadSolutionFileSingleProcessor(std::string &filename, SolutionData3D &S);


};

#endif
