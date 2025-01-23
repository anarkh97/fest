#include<InterpolationOperator.h>

using std::vector;

InterpolationOperator::InterpolationOperator(IoData &iod_, MPI_Comm &comm_)
                     : iod_meta(iod_.interp_driver.meta_input), 
		       iod_spatial(iod_.interp_driver.spatial_interp),
		       comm(comm_)
{
  //
}

void InterpolationOperator::Destroy()
{
  //
}

void
InterpolationOperator::BuildSurfacesToSurfaceMap(TriangulatedSurface &surface)
{
  //
}
