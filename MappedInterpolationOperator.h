#ifndef _MAPPED_INTERPOLATION_OPERATOR_H_
#define _MAPPED_INTERPOLATION_OPERATOR_H_

#include<DynamicLoadOperator.h>
#include<Vector3D.h>
#include<KDTree.h>

typedef KDTree<PointIn3D,3> K3DTree;

class MappedInterpolationOperator : public DynamicLoadOperator {

  //! Maps target surface nodes to nodes from other surfaces
  //! that are provided in the metafile.
  std::map<int, std::vector<Int3>> node2nodes;

  //! The target node index is mapped to a vector of barycentric
  //! weights. The vector is of size 3*N, where N is the number
  //! of proximal (nearby) FSI simulation provided in the metafile.
  std::map<int, std::vector<double>> barycenter_map;

public:

  MappedInterpolationOperator(IoData& iod_, MPI_Comm& comm_);
  ~MappedInterpolationOperator() { }

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

protected:

  void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface) override;

private:

  void BuildKDTree(std::vector<Vec3D> &Xs, K3DTree* tree, std::vector<PointIn3D> &p);

};

#endif
