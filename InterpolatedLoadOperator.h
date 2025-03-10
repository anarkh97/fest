#ifndef _INTERPOLATED_LOAD_OPERATOR_H_
#define _INTERPOLATED_LOAD_OPERATOR_H_

#include<DynamicLoadOperator.h>

/**********************************************
 * class InterpolatedLoadOperator interpolates
 * the interface forces from existing FSI 
 * simulations to compute the forces on the 
 * current (or target) surface. The surfaces
 * are assummed to have identical nodes and
 * elements, i.e. one-to-one mapping. 
 *********************************************/

//! Interpolation load operator 
class InterpolatedLoadOperator : public DynamicLoadOperator {

  //! Stores the solution data of each proximal FSI simulation
  //! provided in the metafile.
  std::vector<SolutionData3D*> neighbor_solutions;

  //! Stores the triangulated interface surfaces from proximal
  //! FSI simulations provided in the metafile
  std::vector<TriangulatedSurface*> neighbor_surfaces;

public:

  InterpolatedLoadOperator(IoData& iod_, MPI_Comm& comm_);
  ~InterpolatedLoadOperator();

  void LoadExistingSurfaces() override;
  void LoadExistingSolutions() override;

  void ComputeForces(TriangulatedSurface& surface, 
                     std::vector<Vec3D>& force,
                     std::vector<Vec3D>& force_over_area, 
                     double t) override;

  void ComputePressures(TriangulatedSurface& surface,
                        std::vector<double>& pressure,
                        double t) override;

  void SetupProjectionMap(TriangulatedSurface &surface) override;

protected:

  // Interpolation methods
  void InterpolateInMetaSpace(TriangulatedSurface& surface, 
                              std::vector<std::vector<Vec3D>>& solutions, 
                              std::vector<double>& pressure);

  void InterpolateInMetaSpace(TriangulatedSurface& surface, 
                              std::vector<std::vector<Vec3D>>& solutions, 
                              std::vector<Vec3D>& force, 
                              std::vector<Vec3D>& force_over_area); 

};

#endif
