#include<ExactNodesOperator.h>

ExactNodesOperator::ExactNodesOperator(SpatialInterpolationData &data, 
                                       MPI_Comm &comm_)
                  : NodalProjectionOperator(data, comm_)
{
  //
}

void
ExactNodesOperator::SetupProjectionMap(TriangulatedSurface &target,
                                       TriangulatedSurface *other)
{

  // do nothing, this is one-to-one map

}

void
ExactNodesOperator::ProjectToTargetSurface(TriangulatedSurface &target, 
                                           TriangulatedSurface *other,
                                           std::vector<Vec3D> &output)
{
 
  if(other) {

    int &active_nodes = target.active_nodes;

    if(other->active_nodes != active_nodes) {
      print_error("*** Error: The number of nodes on target surface and surface "
                  "existing simulation do match.\n");
      exit_mpi();
    }

  }

  return;

}
