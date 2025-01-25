#include<SimpleInterpolationOperator.h>
#include<MathTools/rbf_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

SimpleInterpolationOperator::SimpleInterpolationOperator(IoData &iod_, MPI_Comm &comm_)
                           : DynamicLoadOperator(iod_, comm_)
{
  //
}

//------------------------------------------------------------

void SimpleInterpolationOperator::LoadExistingSurfaces()
{
  // we do not need other surfaces as one-to-one mapping is
  // assumed.
}

//------------------------------------------------------------

void
SimpleInterpolationOperator::ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
                                           std::vector<Vec3D> *force_over_area, double t)
{

}

//------------------------------------------------------------

void
SimpleInterpolationOperator::InterpolateInMetaSpace(TriangulatedSurface &surface, vector<Vec3D> &force, 
                                                    vector<Vec3D> *force_over_area, double t)
{

  int active_nodes = surface.active_nodes;
  int num_points   = iod_meta.numPoints;

  vector<Vec3D> &Xs = surface.X;

  assert((int)Xs.size() == active_nodes); // cracking not supported

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

  //choose a radial basis function for interpolation
  void (*phi)(int, double[], double, double[]); //a function pointer
  switch (iod_meta.basis) {
    case MetaInputData::MULTIQUADRIC :
      phi = MathTools::phi1;  break;
    case MetaInputData::INVERSE_MULTIQUADRIC :
      phi = MathTools::phi2;  break;
    case MetaInputData::THIN_PLATE_SPLINE :
      phi = MathTools::phi3;  break;
    case MetaInputData::GAUSSIAN :
      phi = MathTools::phi4;  break;
    default :
      phi = MathTools::phi4;  break; //Gaussian
  }

  // clear old values
  for(int index=my_start_index; index<my_block_size; ++index) {
    force[index]              = Vec3D(0.0);
    (*force_over_area)[index] = Vec3D(0.0);
  }
  

  // compute forces
  for(int index=my_start_index; index<my_block_size; ++index) {

     // elements attached to this node
     auto elems = surface.node2elem[index];

     // values from existing solutions at this node
     vector<Vec3D> stored_solutions(num_points, Vec3D(0.0));

     for(int i=0; i<num_points; ++i) {
       
       // get data map
       auto solution_map = proxi_solutions[i]->GetMap();  

       // find the time interval
       auto uit = solution_map.lower_bound(t);
       auto lit = std::prev(uit);

       if(uit == solution_map.begin()) {
         lit = uit = solution_map.begin();
	 print_warning("*** Warning: Calculating forces at %e for point %d by "
                       "const extrapolation (outside data interval [%e,%e]).\n",
		       t, i+1, solution_map.begin()->first, solution_map.begin()->second);
       }
       else if(lit == solution_map.end()) {
         lit = uit = --solution_map.end();
	 print_warning("*** Warning: Calculating forces at %e for point %d by "
                       "const extrapolation (outside data interval [%e,%e]).\n",
		       t, i+1, solution_map.begin()->first, solution_map.begin()->second);
       }

       assert((uit != solution_map.end()) and (lit != solution_map.end()));

       // interpolate in time
       double t0 = lit->first;
       double t1 = uit->first;

       Vec3D &s0 = lit->second[index];
       Vec3D &s1 = uit->second[index];
       Vec3D st;

       InterpolateInTime(t0, (double*)s0, t1, (double*)s1, t, (double*)st, 3);

       // update solutions
       stored_solutions[i] = st;

     }

     //InterpolateInMetaSpace(

  }	  
}

//------------------------------------------------------------

void
SimpleInterpolationOperator::InterpolationInSpace()
{

}

//------------------------------------------------------------

//! Copied from DynamicLoadCalculator::InterpolateInTime
void
SimpleInterpolationOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                               double t, double* output, int size)
{
  assert(t2>t1);
  double c1 = (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//------------------------------------------------------------

