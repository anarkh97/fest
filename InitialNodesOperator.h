#ifndef _INITIAL_NODES_OPERATOR_H_
#define _INITIAL_NODES_OPERATOR_H_

#include<NodalProjectionOperator.h>
#include<Vector3D.h>
#include<KDTree.h>

class InitialNodesOperator : public NodalProjectionOperator {

public:

  InitialNodesOperator(SpatialInterpolationData &data, MPI_Comm &comm_);
  ~InitialNodesOperator() { };

  bool NeedSurface() override { return true; };
  void ProjectToTargetSurface(TriangulatedSurface& target_, 
                              TriangulatedSurface* other_,
                              std::vector<Vec3D>& output) override;

  void SetupProjectionMap(TriangulatedSurface& target, 
                          TriangulatedSurface* other) override;
};

#endif
