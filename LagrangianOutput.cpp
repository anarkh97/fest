/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<LagrangianOutput.h>
#include<cstring>
#include<cassert>
#include<fstream>
#include<iomanip> //std::setw

using std::vector;
using std::string;
using std::endl;

//------------------------------------------------------------------------------

LagrangianOutput::LagrangianOutput(MPI_Comm &comm_, LagrangianMeshOutputData &iod_lag_)
                : iod_lag(iod_lag_), comm(comm_), disp_file(NULL), force_file(NULL),
                  fovera_file(NULL)
{
  iFrame = 0;
  last_snapshot_time = -1.0;
}

//------------------------------------------------------------------------------

LagrangianOutput::~LagrangianOutput()
{
}

//------------------------------------------------------------------------------

void
LagrangianOutput::OutputTriangulatedMesh(vector<Vec3D>& X0, vector<Int3>& elems)
{

  if(!strcmp(iod_lag.orig_config,""))
    return; 

  int mpi_rank = -1;
  MPI_Comm_rank(comm, &mpi_rank);

  //Only Proc #0 writes
  if(mpi_rank == 0) {

    char outname[512];
    sprintf(outname, "%s%s", iod_lag.prefix, iod_lag.orig_config);

    std::fstream out;
    out.open(outname, std::fstream::out);

    if(!out.is_open()) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
      exit(-1);
    }

    out << "Nodes MyNodes" << endl;

    for(int i=0; i<(int)X0.size(); i++) 
      out << std::setw(10) << i+1 
          << std::setw(14) << std::scientific << X0[i][0]
          << std::setw(14) << std::scientific << X0[i][1] 
          << std::setw(14) << std::scientific << X0[i][2] << "\n";

    out << "Elements MyElems using MyNodes" << endl;

    for(int i=0; i<(int)elems.size(); i++)
      out << std::setw(10) << i+1 << "  4  "  //"4" for triangle elements
          << std::setw(10) << elems[i][0]+1
          << std::setw(10) << elems[i][1]+1
          << std::setw(10) << elems[i][2]+1 << "\n";

    out.flush();
    out.close();
  }

  MPI_Barrier(comm);

}

//------------------------------------------------------------------------------

void
LagrangianOutput::OutputResults(double t, double dt, int time_step, std::vector<Vec3D>& X0, std::vector<Vec3D>& X,
                                std::vector<Vec3D>& F, std::vector<Vec3D>* F2_ptr, bool force_write)
{

  if(iod_lag.frequency_dt<=0.0 && iod_lag.frequency<=0)
    return;

  if(!isTimeToWrite(t, dt, time_step, iod_lag.frequency_dt, iod_lag.frequency, last_snapshot_time, force_write))
    return;

  int mpi_rank = -1;
  MPI_Comm_rank(comm, &mpi_rank);

  if(mpi_rank != 0)
    goto END_OF_OUTPUT;
  
  //Only Proc #0 writes
  
/* AN: Error not really needed.
  if(!(strcmp(iod_lag.disp,"") || strcmp(iod_lag.force,""))) {
    fprintf(stdout,"\033[0;31m*** Error: Missing output file names.\n\033[0m");
    exit(-1);
  }
*/

  if(strcmp(iod_lag.disp,"")) {

    assert(X0.size() == X.size());

    char outname[512];
    sprintf(outname, "%s%s", iod_lag.prefix, iod_lag.disp);

    if(disp_file == NULL) { //create new file and write header
      disp_file = fopen(outname, "w");
      if(disp_file == NULL) {//unable to open file
        fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
        exit(-1);
      }
      fprintf(disp_file, "Vector DISP under NLDynamic for MyNodes\n");
      fprintf(disp_file, "%d\n", (int)X0.size());
      fclose(disp_file);
    }

    disp_file = fopen(outname,"a");
    if(disp_file == NULL) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
      exit(-1);
    }

    fprintf(disp_file,"%e\n", t);
    for(int i=0; i<(int)X0.size(); i++)
      fprintf(disp_file,"%12.8e    %12.8e    %12.8e\n", X[i][0]-X0[i][0], X[i][1]-X0[i][1], X[i][2]-X0[i][2]);

    fclose(disp_file);
  } 


  if(strcmp(iod_lag.force,"")) {

    char outname[512];
    sprintf(outname, "%s%s", iod_lag.prefix, iod_lag.force);

    if(force_file == NULL) { //create new file and write header
      force_file = fopen(outname, "w");
      if(force_file == NULL) {//unable to open file
        fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
        exit(-1);
      }
      fprintf(force_file, "Vector SOLUTION under NLDynamic for MyNodes\n");
      fprintf(force_file, "%d\n", (int)F.size());
      fclose(force_file);
    }

    force_file = fopen(outname,"a");
    if(force_file == NULL) {//unable to open file
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
      exit(-1);
    }
    AppendResultToFile(force_file, t, F.size(), 3, (double*)F.data());

  } 

  if(strcmp(iod_lag.force_over_area,"")) {
    assert(F2_ptr); // should not be NULL

    char outname[512];
    sprintf(outname, "%s%s", iod_lag.prefix, iod_lag.force);

    if(fovera_file == NULL) { //create new file and write header
      fovera_file = fopen(outname, "w");
      if(fovera_file == NULL) {//unable to open file
        fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
        exit(-1);
      }
      fprintf(fovera_file, "Vector SOLUTION under NLDynamic for MyNodes\n");
      fprintf(fovera_file, "%d\n", (int)F2_ptr->size());
      fclose(fovera_file);
    }

    fovera_file = fopen(outname,"a");
    if(fovera_file == NULL) {//unable to open file
      fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
      exit(-1);
    }
    AppendResultToFile(fovera_file, t, F2_ptr->size(), 3, (double*)F2_ptr->data());
  }

END_OF_OUTPUT:

  MPI_Barrier(comm);

  print("- Wrote solution on a Lagrangian mesh at %e.\n", t);

  last_snapshot_time = t;
  iFrame++;

}

//------------------------------------------------------------------------------
// This function should be called only by Proc #0
void
LagrangianOutput::AppendResultToFile(FILE* file, double time, int N, int dim, double* S)
{
  fprintf(file,"%e\n", time);
  for(int i=0; i<N; i++) {
    for(int j=0; j<dim; j++)
      fprintf(file,"%12.8e    ", S[dim*i+j]);
    fprintf(file,"\n");
  }

  fclose(file);
}

//------------------------------------------------------------------------------




