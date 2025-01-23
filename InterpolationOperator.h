#ifndef _INTERPOLATION_OPERATOR_
#define _INTERPOLATION_OPERATOR_

#include<IoData.h>
#include<TriangulatedSurface.h>
#include<vector>

class InterpolationOperator {

  MetaInputData& iod_meta;
  SpatialInterpolationData& iod_spatial;
  MPI_Comm& comm; 

  //! Radial basis interpolation weights (meta-level)
  std::vector<double> meta_weights;

public:

  InterpolationOperator(IoData& iod_, MPI_Comm& comm_);
  ~InterpolationOperator() { }

  void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface);
  void Destroy();

};

#endif
