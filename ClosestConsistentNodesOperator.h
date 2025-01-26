#ifndef _CLOSEST_CONSISTENT_NODES_OPERATOR_H_
#define _CLOSEST_CONSISTENT_NODES_OPERATOR_H_

#include<DynamicLoadOperator.h>

/**********************************************
 * class ClosestConsistentNodesOperator uses the
 * solution from the point closest to the 
 * target. It currently assumes a one-to-one
 * mapping.
 *********************************************/

//! Simple interpolation operator 
class ClosestConsistentNodesOperator : public DynamicLoadOperator {

public:

  ClosestConsistentNodesOperator(IoData& iod_, MPI_Comm& comm_);
  ~ClosestConsistentNodesOperator() { }

  void LoadExistingSurfaces() override;
  void LoadExistingSolutions() override;

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

protected:

  // Interpolation methods
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size) override;

};

#endif
