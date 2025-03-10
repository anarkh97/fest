#ifndef _DYNAMIC_LOAD_DRIVER_H_
#define _DYNAMIC_LOAD_DRIVER_H_
  
#include<DynamicDriver.h>

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
class DynamicLoadDriver : public DynamicDriver {

  ConcurrentProgramsHandler& concurrent;

public:

  DynamicLoadDriver(IoData& iod_, MPI_Comm& comm_, 
                    ConcurrentProgramsHandler& concurrent_);

  ~DynamicLoadDriver() { }

  void Run() override;

private:

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force, 
                     std::vector<Vec3D>& force_over_area, double t);

};

#endif
