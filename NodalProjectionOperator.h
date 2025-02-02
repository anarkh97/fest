#ifndef _NODAL_PROJECTION_OPERATOR_H_
#define _NODAL_PROJECTION_OPERATOR_H_

#include<Utils.h>
#include<IoData.h>
#include<TriangulatedSurface.h>

/**********************************************
 * class NodalProjectionOperator provides an 
 * abstract interface for implementing different 
 * kinds of nodal projections from exisiting
 * fluid-structure interaction simulations to
 * target simulation wetted surface.
 *********************************************/

class NodalProjectionOperator {

protected:

  SpatialInterpolationData& interp_data;
  MPI_Comm& comm;

  //! Maps target surface nodes to nodes from other surfaces
  //! that are provided in the metafile.
  std::map<int, std::vector<Int3>> node2nodes;

  //! The target node index is mapped to a vector of barycentric
  //! weights. The vector is of size 3*N, where N is the number
  //! of proximal (nearby) FSI simulation provided in the metafile.
  std::map<int, std::vector<double>> barycenter_map;

public:

  NodalProjectionOperator(SpatialInterpolationData& data, MPI_Comm& comm_);
  ~NodalProjectionOperator() { }

  virtual void Destroy();

  virtual bool NeedSurface() = 0;
  virtual void ProjectToTargetSurface(TriangulatedSurface& target_, 
                                      TriangulatedSurface* other_,
                                      std::vector<Vec3D>& output) = 0;
  virtual void SetupProjectionMap(TriangulatedSurface& target, 
                                  TriangulatedSurface* other) = 0;
};

#endif
