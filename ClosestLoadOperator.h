#ifndef _CLOSEST_LOAD_OPERATOR_H_
#define _CLOSEST_LOAD_OPERATOR_H_

#include<SolutionData3D.h>
#include<TriangulatedSurface.h>
#include<DynamicLoadOperator.h>

/**********************************************
 * class ClosestLoadOperator uses the
 * solution from the point closest to the 
 * target. It currently assumes a one-to-one
 * mapping.
 *********************************************/

//! Closest point load operator 
class ClosestLoadOperator : public DynamicLoadOperator {

  SolutionData3D* closest_solution;
  TriangulatedSurface* closest_surface;

public:

  ClosestLoadOperator(IoData& iod_, MPI_Comm& comm_);
  ~ClosestLoadOperator();

  void LoadExistingSurfaces() override;
  void LoadExistingSolutions() override;

  void ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D>& force,
                     std::vector<Vec3D>* force_over_area, double t) override;

  void SetupProjectionMap(TriangulatedSurface& surface) override;

};

#endif
