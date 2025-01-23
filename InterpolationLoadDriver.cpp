#include<InterpolationLoadDriver.h>
#include<cstring>

using std::vector;
using std::shared_ptr;
using std::make_shared;

extern double start_time;

//------------------------------------------------------------

InterpolationLoadDriver::InterpolationLoadDriver(IoData &iod_, MPI_Comm &comm_, 
                                                 ConcurrentProgramsHandler &concurrent_)
                       : iod(iod_), comm(comm_), concurrent(concurrent_), 
		         lagout(comm_, iod_.output)
{
  if(iod.concurrent.aeros.fsi_algo == AerosCouplingData::NONE) {
    print_error("*** Error: Currently FEST operate as a stand-alone program. "
                "It must be coupled with Aero-S.\n");
    exit_mpi();
  }

  // Setup the interpolator.
  ino = new InterpolationOperator(iod, comm);

  // Run time integration. 
  // Note that FEST does not perform time integration, strictly
  // speaking. Instead, it will compute fluid-structure interface forces 
  // based on existing FSI simulation provided by the user.
  Run();
}

//------------------------------------------------------------

InterpolationLoadDriver::~InterpolationLoadDriver()
{
  // smart pointers are automatically deleted.
  if(ino)
    ino->Destroy();
}

//------------------------------------------------------------

//! Code reused from M2C's DynamicLoadCalculator::RunForAeroS
void InterpolationLoadDriver::Run()
{
  
  // objects for the embedded surface and boundary force.
  // Aero-S will update the surface with computed displacementes.
  // We will populate force and send it to Aero-S.
  TriangulatedSurface       surface;
  vector<Vec3D>             force; 
  shared_ptr<vector<Vec3D>> force_over_area = nullptr; // only used for output

  // This will populate the surface object from Aero-S.
  concurrent.InitializeMessengers(&surface, &force);

  // All meta-level setup, such as nearest neighbor weights and surface maps,
  // will be computed here, before time stepping.
  ino->BuildSurfacesToSurfaceMap(surface);
  double overhead_time = walltime();

  // Output recieved surface
  lagout.OutputTriangulatedMesh(surface.X0, surface.elems);

  // allocated space for force over area if the output is requested.
  if(iod.output.force_over_area[0] != 0) 
    force_over_area = make_shared<vector<Vec3D>>(surface.active_nodes, Vec3D(0.)); 

  // ---------------------------------------------------
  // Main time loop (mimics the time loop in Main.cpp
  // ---------------------------------------------------
  double t = 0.0, dt = 0.0, tmax = 0.0;
  int time_step = 0;
  ComputeForces(surface, force, force_over_area.get(), t);

  concurrent.CommunicateBeforeTimeStepping();
  dt   = concurrent.GetTimeStepSize();
  tmax = concurrent.GetMaxTime();
  
  while(t<tmax) {

    time_step++;

    if(t+dt >= tmax)
      dt = tmax - t;

    //---------------------------------
    // Move forward by one time-step
    //---------------------------------
    t += dt;
    print("Step %d: t = %e, dt = %e. Computation time: %.4e s.\n", 
          time_step, t, dt, walltime()-start_time);
 
    ComputeForces(surface, force, force_over_area.get(), t); 

    if(t<tmax) {
      if(time_step==1)
        concurrent.FirstExchange();
      else
        concurrent.Exchange();
    } 

    dt   = concurrent.GetTimeStepSize();
    tmax = concurrent.GetMaxTime(); //set to a small number at final time-step

    lagout.OutputResults(t, dt, time_step, surface.X0, surface.X, force, force_over_area.get(), false);
  }

  concurrent.FinalExchange();

  lagout.OutputResults(t, dt, time_step, surface.X0, surface.X, force, force_over_area.get(), true);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m            NORMAL TERMINATION            \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("Total File I/O Overhead Time: %f sec.\n", overhead_time - start_time);
  print("Total Time Integration Time : %f sec.\n", walltime()    - overhead_time);
  print("Total Computation Time      : %f sec.\n", walltime()    - start_time);
  print("\n");

}

//------------------------------------------------------------

void
InterpolationLoadDriver::ComputeForces(TriangulatedSurface &surface, vector<Vec3D> &force, 
                                       vector<Vec3D> *force_over_area, double t)
{

#ifdef DEBUG_CONNECTION
  ConstantPressureForce(surface, force, force_over_area, t);
  return;
#endif  

  // use the interpolator to estimate interface forces.
  ino->ComputeApproximateForces(surface, force, force_over_area, t);

}

void
InterpolationLoadDriver::ConstantPressureForce(TriangulatedSurface &surface, vector<Vec3D> &force, 
                                               vector<Vec3D> *force_over_area, double t)
{
  // pass a constant pressure * area to aero-s for testing.

  // split the job among the processors (for parallel interpolation)
  // AN: reused from M2C DynamicLoadCalculator::InterpolateInSpace
  int active_nodes = surface.active_nodes;
  int active_elems = surface.active_elems;

  assert(active_nodes == (int)force.size());

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
