#ifndef _DYNAMIC_DRIVER_H_
#define _DYNAMIC_DRIVER_H_

#include<IoData.h>
#include<ConcurrentProgramsHandler.h>
#include<TriangulatedSurface.h>
#include<DynamicLoadOperator.h>
#include<LagrangianOutput.h>
#include<Vector3D.h>
#include<vector>
#include<memory>

/*********************************************
 * Abstract class which handles both 
 * DynamicLoadDriver and DynamicErrorDriver.
 *********************************************/

//! Dynamics driver
class DynamicDriver {

protected:

  MPI_Comm& comm;
  IoData& iod;
  LagrangianOutput lagout;

  std::shared_ptr<DynamicLoadOperator> dlo;

public:

  DynamicDriver(IoData& iod_, MPI_Comm& comm_);
  ~DynamicDriver() { };

  virtual void Run() = 0;
  void Destroy();

};

#endif 
