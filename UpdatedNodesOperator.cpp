#include<UpdatedNodesOperator.h>

UpdatedNodesOperator::UpdatedNodesOperator(SpatialInterpolationData &data, 
                                       MPI_Comm &comm_)
                    : NodalProjectionOperator(data, comm_)
{
  //
}

void
UpdatedNodesOperator::SetupProjectionMap(TriangulatedSurface &target,
                                         TriangulatedSurface &other)
{
  print_error("*** Error: UpdatedNodesOperator::SetupProjectionMap() not "
              "implemented yet.\n");
  exit_mpi();
}

void
UpdatedNodesOperator::ProjectToTargetSurface(TriangulatedSurface &target,
                                             TriangulatedSurface &other,
                                             std::vector<Vec3D> &output)
{
  print_error("*** Error: UpdatedNodesOperator::ProjectToTargetSurface() not "
              "implemented yet.\n");
  exit_mpi();
}
