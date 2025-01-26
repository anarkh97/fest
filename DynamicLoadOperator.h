#ifndef _DYNAMIC_LOAD_OPERATOR_H_
#define _DYNAMIC_LOAD_OPERATOR_H_

#include<IoData.h>
#include<Vector3D.h>
#include<FileHandler3D.h>
#include<TriangulatedSurface.h>
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

  //! Stores the solution data of each proximal FSI simulation
  //! provided in the metafile.
  std::vector<SolutionData3D*> proxi_solutions;

  //! Stores the triangulated interface surfaces from proximal
  //! FSI simulations provided in the metafile
  std::vector<TriangulatedSurface*> proxi_surfaces;

public:

  DynamicLoadOperator(IoData& iod_, MPI_Comm& comm);
  ~DynamicLoadOperator() { }

  // Methods for setup
  virtual void BuildSurfacesToSurfaceMap(TriangulatedSurface& surface);
  virtual void LoadExistingSurfaces();
  virtual void LoadExistingSolutions();

  virtual void Destroy();
  virtual void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                             std::vector<Vec3D> *force_over_area, double t) = 0;

protected:

  // Interpolation methods
  virtual void InterpolateInMetaSpace(TriangulatedSurface &surface, std::vector<std::vector<Vec3D>> &solutions, 
                                      std::vector<Vec3D> &force, std::vector<Vec3D> *force_over_area); 
  virtual void InterpolateInSpace(std::vector<Vec3D>& X, int active_nodes, int dim, double* output);
  virtual void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                 double t, double* output, int size);


};

#endif
