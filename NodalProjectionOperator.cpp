#include<NodalProjectionOperator.h>

NodalProjectionOperator::NodalProjectionOperator(SpatialInterpolationData &data,
                                                 MPI_Comm &comm_)
                       : interp_data(data), comm(comm_)
{
  // 
}

void NodalProjectionOperator::Destroy()
{
  //
}
