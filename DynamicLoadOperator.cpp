#include<ExactNodesOperator.h>
#include<InitialNodesOperator.h>
#include<UpdatedNodesOperator.h>
#include<DynamicLoadOperator.h>
#include<vector>
#include<string>

using std::vector;
using std::string;

extern int verbose;

//------------------------------------------------------------

DynamicLoadOperator::DynamicLoadOperator(IoData &iod_, MPI_Comm &comm_)
                    : iod_meta(iod_.calculator.meta_input),
		      iod_spatial(iod_.calculator.spatial_interp),
		      comm(comm_), file_handler(iod_meta, comm),
		      true_solution(nullptr), npo(nullptr)
{

  switch(iod_spatial.type) {

    case SpatialInterpolationData::EXACT : {
      npo = new ExactNodesOperator(iod_spatial, comm);
      break;
    }
    case SpatialInterpolationData::INITIAL : {
      npo = new InitialNodesOperator(iod_spatial, comm);
      break;
    }
    case SpatialInterpolationData::UPDATED : {
      npo = new UpdatedNodesOperator(iod_spatial, comm);
      break;
    }
    default : {

      print_error("*** Error: Something went wrong while setting up the "
                  "NodalProjectionOperator.\n");
      exit_mpi();

    }
  }

  if(file_handler.GetSolnFileForTarget() != "") {

    true_solution   = new SolutionData3D();
    string filename = file_handler.GetSolnFileForTarget();
    file_handler.ReadSolutionFile(filename, *true_solution);

  } 

}

//------------------------------------------------------------

void DynamicLoadOperator::Destroy()
{
  if(npo) {
    npo->Destroy();
    delete npo;
  }
}

//------------------------------------------------------------

double
DynamicLoadOperator::ComputeError(std::vector<Vec3D> &force_over_area, double t)
{

  if(!true_solution)
    return 0.0;

  int active_nodes = (int)force_over_area.size();

  // find the time interval
  double tk, tkp;
  auto time_bounds = true_solution->GetTimeBounds();

  if(t<time_bounds[0]) {
    tk = tkp = time_bounds[0];
    print_warning("- Warning: Calculating error at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else if(t>time_bounds[1]) {
    tk = tkp = time_bounds[1];
    print_warning("- Warning: Calculating error at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else {
    auto bracket = true_solution->GetTimeBracket(t);
    tk  = bracket[0];
    tkp = bracket[1];
  }

  vector<Vec3D> &Sk  = true_solution->GetSolutionAtTime(tk);
  vector<Vec3D> &Skp = true_solution->GetSolutionAtTime(tkp);
  vector<Vec3D> S(force_over_area.size(), 0.0);

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

  double mse = 0.0;
  for(int index=my_start_index; index<my_block_size; ++index) {

    double diff  = (S[index] - force_over_area[index]).norm();
    mse         += diff*diff;

  }

  mse /= active_nodes;

  // communication
  double error = 0.0;

  // sum error from all ranks.
  MPI_Reduce(&mse, &error, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  return error;

}

//------------------------------------------------------------

//! Copied from DynamicLoadCalculator::InterpolateInTime
void
DynamicLoadOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                        double t, double* output, int size)
{
  assert(t2>t1);
  double c1 = (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//------------------------------------------------------------

