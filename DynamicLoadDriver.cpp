#include<DynamicLoadDriver.h>
#include<cstring>

using std::vector;
using std::shared_ptr;
using std::make_shared;

extern double start_time;
extern int verbose;

//------------------------------------------------------------

DynamicLoadDriver::DynamicLoadDriver(IoData &iod_, MPI_Comm &comm_, 
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
#ifdef DEBUG_CONNECTION // will be convertde to user option
  dlo = new ConstantLoadOperator(iod, comm);
#else
  dlo = new SimpleInterpolationOperator(iod, comm);
#endif
}

//------------------------------------------------------------

void DynamicLoadDriver::Destroy()
{
  if(dlo) {
    dlo->Destroy();
    delete dlo;
  }
}

//------------------------------------------------------------

//! Code reused from M2C's DynamicLoadCalculator::RunForAeroS
void DynamicLoadDriver::Run()
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
  dlo->LoadExistingSurfaces();
  dlo->LoadExistingSolutions();
  dlo->BuildSurfacesToSurfaceMap(surface);
  double overhead_time = walltime();

  if(verbose>1)
    print("- Time taken for file I/O: %f sec.\n", overhead_time - start_time);

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
  print("Total Time For Integration  : %f sec.\n", walltime()    - overhead_time);
  print("Total Computation Time      : %f sec.\n", walltime()    - start_time);
  print("\n");

}

//------------------------------------------------------------

void
DynamicLoadDriver::ComputeForces(TriangulatedSurface &surface, vector<Vec3D> &force, 
                                 vector<Vec3D> *force_over_area, double t)
{

  assert(dlo); // cannot be null
  dlo->ComputeForces(surface, force, force_over_area, t);

}

//------------------------------------------------------------

