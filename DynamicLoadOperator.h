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

  //! true solution at target point. Used for error
  //! calculations.
  SolutionData3D* true_solution;

  //! Projection operator used when surface topology
  //! of different point does not match.
  NodalProjectionOperator *npo;

public:

  DynamicLoadOperator(IoData& iod_, MPI_Comm& comm);
  ~DynamicLoadOperator() { }

  //! target point's surface. Only used when FEST is
  //! called in stand-alone mode.
  void InitializeSurface(TriangulatedSurface* surface);

  virtual void LoadExistingSurfaces() = 0;
  virtual void LoadExistingSolutions() = 0;

  virtual void Destroy();

  virtual void ComputeForces(TriangulatedSurface& surface, 
                             std::vector<Vec3D>& force,
                             std::vector<Vec3D>& force_over_area, 
                             double t) = 0;

  virtual void ComputePressures(TriangulatedSurface& surface, 
                                std::vector<double>& pressure,
                                double t) = 0;

  virtual void SetupProjectionMap(TriangulatedSurface& surface) = 0;

  double ComputeError(std::vector<double>& pressure, double t);

protected:

  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size);

};

#endif
