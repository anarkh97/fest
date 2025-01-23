#ifndef _INTERPOLATION_LOAD_DRIVER_
#define _INTERPOLATION_LOAD_DRIVER_
  
#include<IoData.h>
#include<ConcurrentProgramsHandler.h>
#include<TriangulatedSurface.h>
#include<InterpolationOperator.h>
#include<LagrangianOutput.h>
#include<Vector3D.h>
#include<vector>
#include<memory>

/**********************************************
 * class InterpolationLoadDriver is a special
 * tool that reads solution time-history from
 * user-specified files, calculates (dynamic)
 * loads on an embedded structure, and sends
 * them to a structural dynamics solver (e.g.,
 * Aero-S).
 *********************************************/

//! Interpolation load driver
class InterpolationLoadDriver {

  MPI_Comm& comm;
  IoData& iod;
  ConcurrentProgramsHandler& concurrent;
  LagrangianOutput lagout;

  InterpolationOperator *ino;

public:

  InterpolationLoadDriver(IoData &iod_, MPI_Comm &comm_, 
                          ConcurrentProgramsHandler &concurrent_);

  ~InterpolationLoadDriver();

  void Run();

private:

  void ComputeForces(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                     std::vector<Vec3D> *force_over_area, double t);

  void ConstantPressureForce(TriangulatedSurface &surface, std::vector<Vec3D> &force,
                             std::vector<Vec3D> *force_over_area, double t);
};

#endif
