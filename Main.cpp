#include<ctime>
#include<ConcurrentProgramsHandler.h>
#include<DynamicLoadDriver.h>

int verbose;
double start_time;
MPI_Comm fest_comm;

int main(int argc, char* argv[])
{

  //! Initialize MPI 
  MPI_Init(NULL,NULL); //called together with all concurrent programs -> MPI_COMM_WORLD
  start_time = walltime(); //for timing purpose only (calls MPI_Wtime)

  //! Print header (global proc #0, assumed to be a FEST proc)
  fest_comm = MPI_COMM_WORLD;
  printHeader(argc, argv);

  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  verbose = iod.interp_driver.verbose;

  //! Partition MPI, if there are concurrent programs
  MPI_Comm comm; //this is going to be the FEST communicator
  ConcurrentProgramsHandler concurrent(iod, MPI_COMM_WORLD, comm);
  fest_comm = comm; //correct it

  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();

  //! Special tool for nearest neighbor pressure interpolations.
  DynamicLoadDriver load_driver(iod, comm, concurrent);
  MPI_Barrier(fest_comm); 

  //! finalize 
  load_driver.Destroy();
  concurrent.Destroy();
  MPI_Finalize();
  
  return EXIT_SUCCESS;

}
