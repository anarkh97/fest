#ifndef _SIMPLE_INTERPOLATION_OPERATOR_
#define _SIMPLE_INTERPOLATION_OPERATOR_

#include<DynamicLoadOperator.h>

/**********************************************
 * class SimpleInterpolationOperator interpolates
 * the interface forces from existing FSI 
 * simulations to compute the forces on the 
 * current (or target) surface. The surfaces
 * are assummed to have identical nodes and
 * elements, i.e. one-to-one mapping. 
 *********************************************/

//! Simple interpolation operator 
class SimpleInterpolationOperator : public DynamicLoadOperator {

public:

  SimpleInterpolationOperator(IoData& iod_, MPI_Comm& comm_);
  ~SimpleInterpolationOperator() { }

  void LoadExistingSurfaces() override;

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

protected:

  // Interpolation methods
  void InterpolateInMetaSpace(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                              std::vector<Vec3D> *force_over_area, double t) override;
  void InterpolationInSpace() override;
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size) override;

};

#endif
