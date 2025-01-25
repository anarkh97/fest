#include<ConstantLoadOperator.h>
#include<MathTools/rbf_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

ConstantLoadOperator::ConstantLoadOperator(IoData &iod_, MPI_Comm &comm_)
                    : DynamicLoadOperator(iod_, comm_)
{
  //
}

//------------------------------------------------------------

void ConstantLoadOperator::LoadExistingSurfaces()
{
  // we do not need other surfaces.
}

//------------------------------------------------------------

void ConstantLoadOperator::LoadExistingSolutions()
{
  // we do not need other solutions.
}

//------------------------------------------------------------

void
ConstantLoadOperator::ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                    std::vector<Vec3D> *force_over_area, double t)
{

  // pass a constant pressure * area to aero-s for testing.

  // split the job among the processors (for parallel interpolation)
  // AN: reused from M2C DynamicLoadCalculator::InterpolateInSpace
  int active_nodes = surface.active_nodes;

  assert(active_nodes == (int)surface.X.size()); // cracking not supported

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  int nodes_per_rank = active_nodes / mpi_size;
  int remainder = active_nodes - nodes_per_rank*mpi_size; // left-over nodes.   

  assert(remainder >= 0 and remainder < mpi_size); 

  vector<int> counts(mpi_size, -1);      // stores the number of nodes each rank is responsible for.
  vector<int> start_index(mpi_size, -1); // stores the starting index for each rank in the surface/force vectors.

  for(int i=0; i<mpi_size; ++i) {

    // ranks 0 -- remainder - 1 handle an extra node.
    counts[i] = (i < remainder) ? nodes_per_rank + 1 : nodes_per_rank;
    start_index[i] = (i < remainder) ? (nodes_per_rank + 1)*i : nodes_per_rank*i + remainder;

  }

  assert(start_index.back() + counts.back() == active_nodes);

  int my_node_size = counts[mpi_rank];
  int my_start_index = start_index[mpi_rank];

  // clear old force values
  for(int index=my_start_index; index<my_node_size; ++index) 
    force[index] = Vec3D(0.);

  // compute forces
  Vec3D normalz(0.0, 0.0, 1.0);
  for(int index=my_start_index; index<my_node_size; ++index) {

    // elements associated with this node
    auto elems = surface.node2elem[index];

    // calculate force
    for(int i=0; i<(int)elems.size(); ++i) {
      double area = surface.elemArea[i]/3;
      force[index] += (1e6)*area*normalz;
    }

    if(force_over_area) 
      (*force_over_area)[index] = 1e6*normalz;

  }

  // communication
  for(int i=0; i<mpi_size; i++) {
    counts[i] *= 3;
    start_index[i] *= 3;
  }
  MPI_Allgatherv(MPI_IN_PLACE, 3*my_node_size, MPI_DOUBLE, (double*)force.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);

}

//------------------------------------------------------------

