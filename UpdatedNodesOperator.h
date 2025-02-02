#ifndef _UPDATED_NODES_OPERATOR_H_
#define _UPDATED_NODES_OPERATOR_H_

#include<NodalProjectionOperator.h>

class UpdatedNodesOperator : public NodalProjectionOperator {

public:

  UpdatedNodesOperator(SpatialInterpolationData &data, MPI_Comm &comm_);
  ~UpdatedNodesOperator() { };

  bool NeedSurface() override { return true; };
  void ProjectToTargetSurface(TriangulatedSurface& target_, 
                              TriangulatedSurface& other_,
                              std::vector<Vec3D>& output) override;

  void SetupProjectionMap(TriangulatedSurface& target, 
                          TriangulatedSurface& other) override;
};

#endif
