#include<ClosestConsistentNodesOperator.h>
#include<MathTools/rbf_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

ClosestConsistentNodesOperator::ClosestConsistentNodesOperator(IoData &iod_, MPI_Comm &comm_)
                              : DynamicLoadOperator(iod_, comm_)
{
  //
}

//------------------------------------------------------------

void ClosestConsistentNodesOperator::LoadExistingSurfaces()
{
  // we do not need other surfaces as one-to-one mapping is
  // assumed.
}

//------------------------------------------------------------

void ClosestConsistentNodesOperator::LoadExistingSolutions()
{

  proxi_solutions.resize(1);
  proxi_solutions[0] = new SolutionData3D();
  string &filename   = file_handler.GetSolnFileForProxim(0);
  file_handler.ReadSolutionFile(filename, *proxi_solutions[0]);

}

//------------------------------------------------------------

void
ClosestConsistentNodesOperator::ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                              std::vector<Vec3D> *force_over_area, double t)
{
  int num_points   = iod_meta.numPoints;
  int active_nodes = surface.active_nodes;

  assert((int)surface.X.size() == active_nodes); // cracking not supported.

  // clear existing data
  for(int i=0; i<active_nodes; ++i) {
    force[i] = Vec3D(0.0);
    if(force_over_area)
      (*force_over_area)[i] = Vec3D(0.0);
  }

  // find the time interval
  double tk, tkp;
  auto time_bounds = proxi_solutions[0]->GetTimeBounds();

  if(t<time_bounds[0]) {
    tk = tkp = time_bounds[0];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else if(t>time_bounds[1]) {
    tk = tkp = time_bounds[1];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else {

    auto bracket = proxi_solutions[0]->GetTimeBracket(t);
    tk  = bracket[0];
    tkp = bracket[1];

  }

  //fprintf(stdout, "Bracket (%e, %e) for time %e.\n", tk, tkp, t);

  // Get solutions at tk and tkp
  vector<Vec3D> &Sk  = proxi_solutions[0]->GetSolutionAtTime(tk);
  vector<Vec3D> &Skp = proxi_solutions[0]->GetSolutionAtTime(tkp);
  vector<Vec3D> S(active_nodes, Vec3D(0.0));

  if((int)Sk.size() != active_nodes) {
    print_error("*** Error: The number of nodes for point 1 at time (t = %e) "
                "does not match the number of nodes on the target surface.\n",
      	        tk);
    exit_mpi();
  }
  if((int)Skp.size() != active_nodes) {
    print_error("*** Error: The number of nodes for point 1 at time (t = %e) "
                "does not match the number of nodes on the target surface.\n",
      	        tkp);
    exit_mpi();
  }

  InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                    (double*)S.data(), 3*active_nodes);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
  int nodes_per_rank = active_nodes / mpi_size;
  int remainder      = active_nodes - nodes_per_rank*mpi_size;

  assert(remainder >=0 and remainder < mpi_size);

  vector<int> counts(mpi_size, -1);
  vector<int> start_index(mpi_size, -1);

  for(int i=0; i<mpi_size; ++i) {
    counts[i]      = (i < remainder) ? nodes_per_rank + 1 : nodes_per_rank;
    start_index[i] = (i < remainder) ? (nodes_per_rank + 1)*i : nodes_per_rank*i;
  }

  assert(start_index.back() + counts.back() == active_nodes);

  int my_block_size  = counts[mpi_rank];
  int my_start_index = start_index[mpi_rank];

  for(int index=my_start_index; index<my_block_size; ++index) {

    // get current area and normal
    Vec3D patch(0.0);
    auto elems  = surface.node2elem[index];
    int  nelems = elems.size();
    for(auto it=elems.begin(); it!=elems.end(); ++it)
      patch += surface.elemArea[*it]*surface.elemNorm[*it];
    
    patch /= nelems;

    double area   = patch.norm();
    Vec3D  normal = patch/area;

    // calculate forces from pressures.
    double pressure = S[index].norm();

    if(force_over_area) {
      (*force_over_area)[index][0] = pressure*normal[0]; 
      (*force_over_area)[index][1] = pressure*normal[1]; 
      (*force_over_area)[index][2] = pressure*normal[2]; 
    }

    force[index] = pressure*area*normal;

  }

  exit(-1);

  // communication
  for(int i=0; i<mpi_size; i++) {
    counts[i] *= 3;
    start_index[i] *= 3;
  }
  MPI_Allgatherv(MPI_IN_PLACE, 3*my_block_size, MPI_DOUBLE, (double*)force.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);
  if(force_over_area)
    MPI_Allgatherv(MPI_IN_PLACE, 3*my_block_size, MPI_DOUBLE, (double*)force_over_area->data(), 
                   counts.data(), start_index.data(), MPI_DOUBLE, comm);

}

//------------------------------------------------------------

//! Copied from DynamicLoadCalculator::InterpolateInTime
void
ClosestConsistentNodesOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                               double t, double* output, int size)
{
  assert(t2>t1);
  double c1 = (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//------------------------------------------------------------

