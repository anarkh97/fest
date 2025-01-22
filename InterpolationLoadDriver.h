#ifndef _INTERPOLATION_LOAD_DRIVER_
#define _INTERPOLATION_LOAD_DRIVER_
  
#include<IoData.h>
#include<ConcurrentProgramsHandler.h>
#include<TriangulatedSurface.h>
#include<LagrangianOutput.h>
#include<Vector3D.h>
#include<vector>
#include<memory>

//! Interpolation load driver
class InterpolationLoadDriver {

  MPI_Comm& comm;
  IoData& iod;
  ConcurrentProgramsHandler& concurrent;
  LagrangianOutput lagout;

public:

  InterpolationLoadDriver(IoData &iod_, MPI_Comm &comm_, 
                          ConcurrentProgramsHandler &concurrent_);

  ~InterpolationLoadDriver();

  void Run();

protected:

  void ComputeForces(TriangulatedSurface &surface, std::vector<Vec3D> &force, 
                     std::vector<Vec3D> *force_over_area, double t);

};

#endif
