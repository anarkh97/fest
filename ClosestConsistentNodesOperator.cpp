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

  // Get solutions at tk and tkp
  vector<Vec3D> &Sk  = proxi_solutions[0]->GetSolutionAtTime(tk);
  vector<Vec3D> &Skp = proxi_solutions[0]->GetSolutionAtTime(tkp);

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
                    (double*)force.data(), active_nodes);

  // correct forces to forces*area
  for(int i=0; i<active_nodes; ++i) {

    auto elems = surface.node2elem[i];
    
    double area;
    for(auto it=elems.begin(); it!=elems.end(); ++it) 
      area   += surface.elemArea[*it]/3;

    if(force_over_area) {
      (*force_over_area)[i][0] = force[i][0];
      (*force_over_area)[i][1] = force[i][1];
      (*force_over_area)[i][2] = force[i][2];
    }
    force[i] *= area;

  }

  MPI_Barrier(comm);

}

//------------------------------------------------------------

//! Copied from DynamicLoadCalculator::InterpolateInTime
void
ClosestConsistentNodesOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                               double t, double* output, int size)
{
  assert(t2>=t1);

  double c1 = (t2 == t1) ? 1.0 : (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//------------------------------------------------------------

