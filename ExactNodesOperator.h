#ifndef _EXACT_NODES_OPERATOR_H_
#define _EXACT_NODES_OPERATOR_H_

#include<NodalProjectionOperator.h>

class ExactNodesOperator : public NodalProjectionOperator {

public:
  
  ExactNodesOperator(SpatialInterpolationData& data, MPI_Comm& comm_);
  ~ExactNodesOperator() { }

  bool NeedSurface() override { return false; };
  void ProjectToTargetSurface(TriangulatedSurface& target_, 
                              TriangulatedSurface& other_,
                              std::vector<Vec3D>& output) override;

  void SetupProjectionMap(TriangulatedSurface& target, 
                          TriangulatedSurface& other) override;
};

#endif
