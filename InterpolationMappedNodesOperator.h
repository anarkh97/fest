#ifndef _INTERPOLATION_MAPPED_NODES_OPERATOR_H_
#define _INTERPOLATION_MAPPED_NODES_OPERATOR_H_

#include<DynamicLoadOperator.h>
#include<Vector3D.h>
#include<KDTree.h>

typedef KDTree<PointIn3D,3> K3DTree;

class InterpolationMappedNodesOperator : public DynamicLoadOperator {

  //! Maps target surface nodes to nodes from other surfaces
  //! that are provided in the metafile.
  std::map<int, std::vector<Int3>> node2nodes;

  //! The target node index is mapped to a vector of barycentric
  //! weights. The vector is of size 3*N, where N is the number
  //! of proximal (nearby) FSI simulation provided in the metafile.
  std::map<int, std::vector<double>> barycenter_map;

public:

  InterpolationMappedNodesOperator(IoData& iod_, MPI_Comm& comm_);
  ~InterpolationMappedNodesOperator() { }

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

protected:

  void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface) override;
  void InterpolateInMetaSpace(TriangulatedSurface &surface, std::vector<std::vector<Vec3D>> &solutions, 
                              std::vector<Vec3D> &force, std::vector<Vec3D> *force_over_area) override;
  void InterpolateInSpace(std::vector<Vec3D>& X, int active_nodes, int dim, double* output) override;
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size) override;

private:

  void BuildKDTree(std::vector<Vec3D> &Xs, K3DTree* tree, std::vector<PointIn3D> &p);

};

#endif
