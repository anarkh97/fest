#ifndef _DYNAMIC_ERROR_DRIVER_H_
#define _DYNAMIC_ERROR_DRIVER_H_
  
#include<DynamicDriver.h>

/**********************************************
 *
 *********************************************/

//! Dynamic load driver
class DynamicErrorDriver : public DynamicDriver {

public:

  DynamicErrorDriver(IoData &iod_, MPI_Comm &comm_); 

  ~DynamicErrorDriver() { }

  void Run() override;

private:

  void ComputePressures(TriangulatedSurface& surf,
                        std::vector<double>& pressure,
                        double t);

  void ComputeError(std::vector<double> &pressure, 
                    double t, double &error);

};

#endif
