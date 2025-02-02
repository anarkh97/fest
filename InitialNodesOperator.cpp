#include<InitialNodesOperator.h>

InitialNodesOperator::InitialNodesOperator(SpatialInterpolationData &data, 
                                       MPI_Comm &comm_)
                    : NodalProjectionOperator(data, comm_)
{
  //
}

void
InitialNodesOperator::SetupProjectionMap(TriangulatedSurface &target,
                                         TriangulatedSurface &other)
{
  print_error("*** Error: InitialNodesOperator::SetupProjectionMap() not "
              "implemented yet.\n");
  exit_mpi();
}

void
InitialNodesOperator::ProjectToTargetSurface(TriangulatedSurface &target,
                                             TriangulatedSurface &other,
                                             std::vector<Vec3D> &output)
{
  print_error("*** Error: InitialNodesOperator::ProjectToTargetSurface() not "
              "implemented yet.\n");
  exit_mpi();
}
