#ifndef _INTERPOLATION_OPERATOR_
#define _INTERPOLATION_OPERATOR_

#include<IoData.h>
#include<KDTree.h>
#include<Vector3D.h>
#include<SolutionData3D.h>
#include<TriangulatedSurface.h>
#include<string>
#include<vector>
#include<map>

typedef KDTree<PointIn3D,3> K3DTree;

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

  //! Maps target surface nodes to nodes from other surfaces
  //! that are provided in the metafile.
  std::map<int, std::vector<Int3>> node2nodes;

  //! The target node index is mapped to a vector of barycentric
  //! weights. The vector is of size 3*N, where N is the number
  //! of proximal (nearby) FSI simulation provided in the metafile.
  std::map<int, std::vector<double>> barycenter_map;

  //! Stores the solution data of each proximal FSI simulation
  //! provided in the metafile.
  std::vector<SolutionData3D> proxi_solutions;

public:

  InterpolationOperator(IoData& iod_, MPI_Comm& comm_);
  ~InterpolationOperator() { }

  void SetupInterpolator(TriangulatedSurface& surface_);
  void Destroy();

  void ComputeApproximateForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                std::vector<Vec3D> *force_over_area, double t);

private:

  // Methods for setup
  void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface);
  void LoadExistingSolutions();

  // Methods for reading surface and solution files.
  void ReadMetaFile();
  void ReadMeshFile(const char *filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);
  void ReadSolutionFile(const char *filename, SolutionData3D &S);
  void ReadMeshFileInTopFormat(const char *filename, std::vector<Vec3D> &Xs, std::vector<Int3> &Es);

  // Auxillary methods
  void BuildKDTree(std::vector<Vec3D> &Xs, K3DTree* tree, std::vector<PointIn3D> &p);

  // Interpolation methods
  void InterpolateInMetaSpaceNoMap(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                                   std::vector<Vec3D> *force_over_area, double t);
  void InterpolateInMetaSpaceWithMap(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                                     std::vector<Vec3D> *force_over_area, double t);
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size);


};

#endif
