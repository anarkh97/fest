#ifndef _CLOSEST_POINT_LOAD_OPERATOR_H_
#define _CLOSEST_POINT_LOAD_OPERATOR_H_

#include<DynamicLoadOperator.h>

/**********************************************
 * class ClosestPointLoadOperator uses the
 * solution from the point closest to the 
 * target. It currently assumes a one-to-one
 * mapping.
 *********************************************/

//! Simple interpolation operator 
class ClosestPointLoadOperator : public DynamicLoadOperator {

public:

  ClosestPointLoadOperator(IoData& iod_, MPI_Comm& comm_);
  ~ClosestPointLoadOperator() { }

  void LoadExistingSurfaces() override;

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

protected:

  // Interpolation methods
  void InterpolateInMetaSpace(TriangulatedSurface &surface, std::vector<std::vector<Vec3D>> &solutions, 
                              std::vector<Vec3D> &force, std::vector<Vec3D> *force_over_area) override; 
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size) override;

};

#endif
