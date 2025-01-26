#ifndef _CONSTANT_LOAD_OPERATOR_H_
#define _CONSTANT_LOAD_OPERATOR_H_ 

#include<DynamicLoadOperator.h>

/**********************************************
 * class ConstantLoadOperator sends a constant
 * pressure force to aero-s. This class is
 * mainly used for debugging. But in future
 * can be expanded to user specified pressure
 * values.
 *********************************************/

//! Simple interpolation operator 
class ConstantLoadOperator : public DynamicLoadOperator {

public:

  ConstantLoadOperator(IoData& iod_, MPI_Comm& comm_);
  ~ConstantLoadOperator() { }

  void LoadExistingSurfaces() override;
  void LoadExistingSolutions() override;

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

};

#endif
