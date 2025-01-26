#ifndef _DYNAMIC_LOAD_DRIVER_H_
#define _DYNAMIC_LOAD_DRIVER_H_
  
#include<IoData.h>
#include<ConcurrentProgramsHandler.h>
#include<TriangulatedSurface.h>
#include<DynamicLoadOperator.h>
#include<ConstantLoadOperator.h>
#include<SimpleInterpolationOperator.h>
#include<LagrangianOutput.h>
#include<Vector3D.h>
#include<vector>
#include<memory>

/**********************************************
 * The DynamicLoadDriver class computes forces 
 * on the shared fluid-structure interface using 
 * a dynamic load estimator specified in the 
 * input file. These forces are communicated to 
 * the Aero-S, structural dynamics solver,
 * effectively mimicking the fluid side in a 
 * fluid-structure interaction simulation.
 *********************************************/

//! Dynamic load driver
class DynamicLoadDriver {

  MPI_Comm& comm;
  IoData& iod;
  ConcurrentProgramsHandler& concurrent;
  LagrangianOutput lagout;

  DynamicLoadOperator *dlo;

public:

  DynamicLoadDriver(IoData &iod_, MPI_Comm &comm_, 
                    ConcurrentProgramsHandler &concurrent_);

  ~DynamicLoadDriver() { }

  void Run();
  void Destroy();

private:

  void ComputeForces(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                     std::vector<Vec3D> *force_over_area, double t);

};

#endif
