#ifndef _INTERPOLATION_LOAD_DRIVER_
#define _INTERPOLATION_LOAD_DRIVER_
  
#include<IoData.h>
#include<ConcurrentProgramsHandler.h>

//! Interpolation load driver
class InterpolationLoadDriver {

  MPI_Comm& comm;
  IoData& iod;
  ConcurrentProgramsHandler& concurrent;

public:

  InterpolationLoadDriver(IoData &iod_, MPI_Comm &comm_, 
                          ConcurrentProgramsHandler &concurrent_);

  ~InterpolationLoadDriver();

  void Run();

};

#endif
