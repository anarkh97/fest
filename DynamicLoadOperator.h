#ifndef _DYNAMIC_LOAD_OPERATOR_H_
#define _DYNAMIC_LOAD_OPERATOR_H_

#include<IoData.h>
#include<Vector3D.h>
#include<FileHandler3D.h>
#include<TriangulatedSurface.h>
#include<NodalProjectionOperator.h>
#include<string>
#include<vector>
#include<map>


/**********************************************
 * class DynamicLoadOperator provides an abstract
 * interface for implementing different kinds of
 * dynamic load calculators. These classes are
 * responsible for computing the forces on the
 * shared fluid-structure interface.
 *********************************************/

//! Dynamic load operator
class DynamicLoadOperator {

protected:

  MetaInputData& iod_meta;
  SpatialInterpolationData& iod_spatial;
  MPI_Comm& comm; 
  FileHandler3D file_handler;

  NodalProjectionOperator *npo;

public:

  DynamicLoadOperator(IoData& iod_, MPI_Comm& comm);
  ~DynamicLoadOperator() { }

  virtual void LoadExistingSurfaces() = 0;
  virtual void LoadExistingSolutions() = 0;

  virtual void Destroy();
  virtual void ComputeForces(TriangulatedSurface &surface, std::vector<Vec3D> &force,
                             std::vector<Vec3D> *force_over_area, double t) = 0;

  virtual void SetupProjectionMap(TriangulatedSurface &surface) = 0;

protected:

  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size);


};

#endif
