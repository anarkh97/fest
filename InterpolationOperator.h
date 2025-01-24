#ifndef _INTERPOLATION_OPERATOR_
#define _INTERPOLATION_OPERATOR_

#include<IoData.h>
#include<Vector3D.h>
#include<TriangulatedSurface.h>
#include<string>
#include<vector>
#include<map>

class InterpolationOperator {

  MetaInputData& iod_meta;
  SpatialInterpolationData& iod_spatial;
  MPI_Comm& comm; 

  //! paths of stored (or existing) surfaces
  std::vector<std::string> proxi_surface_files;
  //! paths of stored (or existing) solutions
  std::vector<std::string> proxi_solution_files;

  //! meta-level parameters provided by user.
  std::vector<double> target;
  std::vector<std::vector<double>> proximals;

  //! maps the order of fields specified in the metafile.
  //! mainly used for book-keeping. We expect a certain 
  //! order, which is checked for while reading the metafile.
  enum Fields {DIRECTORY=0, MESH_FILE, SOLUTION_FILE, PARAMETER_START,
               PARAMETER_END};
  std::map<Fields, int> field2col;


public:

  InterpolationOperator(IoData& iod_, MPI_Comm& comm_);
  ~InterpolationOperator() { }

  void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface);
  void Destroy();

  void ComputeApproximateForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                std::vector<Vec3D> *force_over_area, double t);

private:

  void ReadMetaFile();
  void ReadMeshFile(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es);
  void ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es);

};

#endif
