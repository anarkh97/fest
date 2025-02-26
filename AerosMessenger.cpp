/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<AerosMessenger.h>
#include<cassert>
#include<climits>

#define SUGGEST_DT_TAG 444
#define WET_SURF_TAG1 555
#define WET_SURF_TAG2 666
#define WET_SURF_TAG3 888
#define WET_SURF_TAG4 999
#define SUBCYCLING_TAG 777

#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

#define FORCE_TAG 1000
#define DISP_TAG 2000
#define INFO_TAG 3000

using std::vector;

extern int verbose;

//---------------------------------------------------------------
// TODO(KW): When I or someone has time, should re-visit the different coupling algorithms (C0, A6, etc.).
//
AerosMessenger::AerosMessenger(AerosCouplingData &iod_aeros_, MPI_Comm &fest_comm_, MPI_Comm &joint_comm_, 
                               TriangulatedSurface &surf_, vector<Vec3D> &F_)
              : iod_aeros(iod_aeros_), fest_comm(fest_comm_), joint_comm(joint_comm_),
                surface(surf_), F(F_), numStrNodes(NULL)
{

  MPI_Comm_rank(fest_comm, &fest_rank);
  MPI_Comm_size(fest_comm, &fest_size);

  // verification of MPI / communicators
  int joint_rank, joint_size;
  MPI_Comm_rank(joint_comm, &joint_rank);
  MPI_Comm_size(joint_comm, &joint_size);
  assert(fest_rank == joint_rank);
  assert(fest_size == joint_size);

  MPI_Comm_remote_size(joint_comm, &numAerosProcs);
  //fprintf(stdout,"I am [%d]. Aero-S running on %d procs.\n", fest_rank, numAerosProcs);

  //----------------------
  // copied from AERO-F/DynamicNodalTransfer.cpp (written by KW in grad school...)
  //----------------------
  bufsize = 6; //6 dofs per node (disp and velo)

  bool crack;
  int nStNodes, nStElems, totalStNodes, totalStElems;

  GetEmbeddedWetSurfaceInfo(elemType, crack, nStNodes, nStElems);

  totalStNodes = nStNodes;
  totalStElems = nStElems;

  // allocate memory for the node list
  totalNodes = totalStNodes;
  nNodes     = nStNodes;
  surface.X.resize(totalNodes, Vec3D(0.0));
  surface.X0.resize(totalNodes, Vec3D(0.0));
  surface.Udot.resize(totalNodes, Vec3D(0.0));
  F.resize(totalNodes, Vec3D(0.0));
  temp_buffer.resize(2*bufsize*totalNodes, 0.0); //a long temporary buffer, can be used for anything

  // allocate memory for the element topology list
  int tmpTopo[nStElems][4];
  switch(elemType) {
    case 3: // all triangles
      totalElems = totalStElems;
      nElems     = nStElems;
      surface.elems.resize(totalElems, Int3(0));
      GetEmbeddedWetSurface(nNodes, surface.X0.data(), nElems, (int*)surface.elems.data(), elemType);
      break;
    case 4: // quadrangles include triangles represented as degenerated quadrangles.
      GetEmbeddedWetSurface(nNodes, surface.X0.data(), nStElems, (int*)tmpTopo, elemType);
      nElems = SplitQuads((int *)tmpTopo, nStElems, surface.elems); 
      totalElems = nElems;
      break;
    default:
      print_error("*** Error: Element type (%d) of the wet surface not recognized! Must be 3 or 4.\n", elemType);
      exit_mpi();
  }
  for(int i = 0; i < nNodes; i++)
    surface.X[i] = surface.X0[i];


  surface.BuildConnectivities();
  surface.CalculateNormalsAndAreas();


  surface.active_nodes = nNodes;
  surface.active_elems = nElems;

  Negotiate(); // Following the AERO-F/S function name, although misleading

  GetInfo(); // Get algorithm number, dt, and tmax

  structureSubcycling = (algNum == 22) ? GetStructSubcyclingInfo() : 0; 
  //currently, M2C/FEST does not support "structure subcycling".

  if(algNum == 6)
    print("- Coupled with Aero-S (running on %d processors) using the A6 algorithm "
          "(Assuming fixed time-step in Aero-S).\n", numAerosProcs);
  else if(algNum == 22)
    print("- Coupled with Aero-S (running on %d processors) using the C0 algorithm.\n", numAerosProcs);



//debug
/*
  if(fest_rank==0) {
    vector<Vec3D> &x(surface.X);
    for(int i=0; i<x.size(); i++)
      fprintf(stdout,"%d %e %e %e.\n", i, x[i][0], x[i][1], x[i][2]);
    vector<Int3> &e(surface.elems);
    for(int i=0; i<e.size(); i++)
      fprintf(stdout,"%d %d %d %d.\n", i, e[i][0], e[i][1], e[i][2]);
  }
  MPI_Barrier(fest_comm);
  exit_mpi();
*/
} 

//---------------------------------------------------------------

AerosMessenger::~AerosMessenger()
{
  if(numStrNodes)
    delete [] numStrNodes;
}

//---------------------------------------------------------------

void
AerosMessenger::Destroy()
{ }

//---------------------------------------------------------------

void
AerosMessenger::GetEmbeddedWetSurfaceInfo(int &eType, bool &crack, int &nStNodes, int &nStElems)
{
  int info[4];
  if(fest_rank==0)
    MPI_Recv(info, 4, MPI_INT, MPI_ANY_SOURCE, WET_SURF_TAG1, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(info, 4, MPI_INT, 0, fest_comm);

  eType    = info[0];
  crack    = false; //! AN: cracking disabled as it would result in errors in interpolations.
  nStNodes = info[2]; 
  nStElems = info[3];

  print("- Embedded Surface from Aero-S: Number of Nodes/Elements: %d/%d, Element Type = %d, Fracture = %d.\n", 
        nStNodes, nStElems, eType, (int)crack);
}

//---------------------------------------------------------------

void
AerosMessenger::GetEmbeddedWetSurface(int nNodes, Vec3D *nodes, int nElems, int *elems, int eType)
{
  if(fest_rank==0)
    MPI_Recv((double*)nodes, nNodes*3, MPI_DOUBLE, MPI_ANY_SOURCE, WET_SURF_TAG2, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast((double*)nodes, nNodes*3, MPI_DOUBLE, 0, fest_comm);

  if(fest_rank==0)
    MPI_Recv(elems, nElems*eType, MPI_INT, MPI_ANY_SOURCE, WET_SURF_TAG3, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(elems, nElems*eType, MPI_INT, 0, fest_comm);

  if(verbose>=1)
    print("- Received nodes and elements of embedded surface from Aero-S.\n");

/*
  if(fest_rank==0) {
    for(int i=0; i<nNodes; i++) {
      fprintf(stdout,"%d  %e %e %e\n", i, nodes[i][0], nodes[i][1], nodes[i][2]);
    }
  }
*/
}

//---------------------------------------------------------------

int
AerosMessenger::SplitQuads(int *quads, int nStElems, vector<Int3> &Tria)
{
  int nTrias = 0;
  for(int i = 0; i < nStElems; ++i) {
    if(quads[i * 4 + 2] == quads[i * 4 + 3]) {
      nTrias += 1;
    }
    else {
      nTrias += 2;
    }
  }
  Tria.assign(nTrias, Int3(0));
  int count = 0;
  for(int i = 0; i < nStElems; ++i) {
    Tria[count][0] = quads[i * 4];
    Tria[count][1] = quads[i * 4 + 1];
    Tria[count][2] = quads[i * 4 + 2];
    count++;
    if(quads[i * 4 + 2] == quads[i * 4 + 3]) {
      continue;
    }
    Tria[count][0] = quads[i * 4];
    Tria[count][1] = quads[i * 4 + 2];
    Tria[count][2] = quads[i * 4 + 3];
    count++;
  }

  assert(count == nTrias);

  return nTrias;
}

//---------------------------------------------------------------

void
AerosMessenger::Negotiate()
{

  int numCPUMatchedNodes = fest_rank ? 0 : totalNodes;

  vector<int> ibuffer;
  if(numCPUMatchedNodes>0) {
    ibuffer.resize(numCPUMatchedNodes);
    for(int i=0; i<numCPUMatchedNodes; i++)
      ibuffer[i] = i; //trivial for M2C/Embedded boundary. (non-trivial in AERO-F/S w/ ALE)
  }

  // send the matched node numbers of this fluid CPU to all the structure CPUs
  vector<MPI_Request> send_requests;
  for (int proc = 0; proc < numAerosProcs; proc++) {
    send_requests.push_back(MPI_Request());
    MPI_Isend(&numCPUMatchedNodes, 1, MPI_INT, proc, NEGO_NUM_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
    if(numCPUMatchedNodes>0)
      MPI_Isend(ibuffer.data(), numCPUMatchedNodes, MPI_INT, proc, NEGO_BUF_TAG, 
                joint_comm, &(send_requests[send_requests.size()-1]));
  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  // receive the list of matched nodes that each structure CPU contains
  if(numCPUMatchedNodes > 0) {

    local2pack.assign(totalNodes,-INT_MAX);

    numStrNodes = new int[numAerosProcs][2];
    for(int proc = 0; proc < numAerosProcs; ++proc) {
      MPI_Recv(&numStrNodes[proc][0], 1, MPI_INT, proc, NEGO_NUM_TAG, joint_comm, MPI_STATUS_IGNORE);
      //fprintf(stdout,"Got %d from proc %d.\n", numStrNodes[proc][0], proc);

      if(numStrNodes[proc][0] > 0) {
        MPI_Recv(ibuffer.data(), numStrNodes[proc][0], MPI_INT, proc, NEGO_BUF_TAG, joint_comm, MPI_STATUS_IGNORE);

        pack2local.reserve(pack2local.size() + numStrNodes[proc][0]);

        for(int i = 0; i < numStrNodes[proc][0]; ++i) {
          int idx = ibuffer[i];
          assert(idx>=0 && idx<totalNodes);

          pack2local.push_back(idx);
          local2pack[idx] = pack2local.size()-1;
        }
      }
    }

    if((int)pack2local.size() != numCPUMatchedNodes) {
      fprintf(stdout, "\033[0;31m*** Error (proc %d): wrong number of matched nodes (%d instead of %d).\n\033[0m",
              fest_rank, (int)pack2local.size(), numCPUMatchedNodes);
      exit(-1);
    }

    for(auto it = local2pack.begin(); it != local2pack.end(); it++) {
      if(*it < 0) {
        fprintf(stdout, "\033[0;31m*** Error (proc %d): found unmatched node (local id: %d).\n\033[0m",
              fest_rank, (int)(it - local2pack.begin()));
        exit(-1);
      }
    }

    numStrNodes[0][1] = 0;
    for(int proc = 1; proc < numAerosProcs; ++proc) {
      numStrNodes[proc][1] = numStrNodes[proc - 1][1] + numStrNodes[proc - 1][0];
    }
  }
  
}

//---------------------------------------------------------------

void
AerosMessenger::GetInfo()
{

  double info[5];
  if(fest_rank == 0) {
    MPI_Recv(info, 5, MPI_DOUBLE, MPI_ANY_SOURCE, INFO_TAG, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(info, 5, MPI_DOUBLE, 0, fest_comm);

  algNum = int(info[0]);
  dt     = info[1];
  tmax   = info[2];

  //int rstrt = int(info[3]); //not used
  //int smode = int(info[4]); //not used

  // check for consistency in algorithm number
  if(iod_aeros.fsi_algo == AerosCouplingData::A6)
    assert(algNum == 6);
  if(iod_aeros.fsi_algo == AerosCouplingData::C0)
    assert(algNum == 22);

  if(algNum != 6 && algNum != 22) {
    print_error("*** Error: Detected unsupported FSI algorithm from Aero-S (%d).\n", algNum);
    exit_mpi();
  }

  if(algNum == 6) {
    tmax -= 0.5 * dt;
  }
  if(algNum == 22) {
    tmax += 0.5 * dt;
  }

  if(verbose>1) 
    print("- Received from Aero-S: dt = %e, tmax(+/-0.5dt) = %e.\n", dt, tmax);
}

//---------------------------------------------------------------

int
AerosMessenger::GetStructSubcyclingInfo()
{
  double info;
  if(fest_rank == 0) {
    MPI_Recv(&info, 1, MPI_DOUBLE, MPI_ANY_SOURCE, SUBCYCLING_TAG, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&info, 1, MPI_DOUBLE, 0, fest_comm);

  return (int)info;
}

//---------------------------------------------------------------

void
AerosMessenger::SendForce()
{
  //IMPORTANT: Assuming that the force has been assembled on Proc 0
  if(fest_rank == 0) {

    //TODO: Need to take care of "staggering"
    //

    // prepare package
    for(int i=0; i<(int)pack2local.size(); i++) {
      for(int j=0; j<3; j++)
        temp_buffer[3*i+j] = F[pack2local[i]][j];
    //  fprintf(stdout,"temp_buffer[%d] = %e %e %e.\n", i, temp_buffer[3*i], temp_buffer[3*i+1], temp_buffer[3*i+2]);
    }

    vector<MPI_Request> send_requests;

    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        int size = 3*numStrNodes[proc][0];

        double *localBuffer = temp_buffer.data() + 3*numStrNodes[proc][1];

        send_requests.push_back(MPI_Request());
        MPI_Isend(localBuffer, size, MPI_DOUBLE, proc, FORCE_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  }

}

//---------------------------------------------------------------

void
AerosMessenger::GetDisplacementAndVelocity()
{
  if(fest_rank == 0) {

    assert(bufsize==6);

    // get disp and velo to temp_buffer
    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        int size = bufsize*numStrNodes[proc][0]; 
        int first_index =  bufsize*numStrNodes[proc][1];

        double *localBuffer = temp_buffer.data() + first_index;

        MPI_Recv(localBuffer, size, MPI_DOUBLE, proc, DISP_TAG, joint_comm, MPI_STATUS_IGNORE);
      }
    }

    // apply disp and velo to surface and Udot
    for(int i=0; i<nNodes; i++) {
      int id = local2pack[i];
      for(int j=0; j<3; j++)
        surface.X[i][j] = surface.X0[i][j] + temp_buffer[bufsize*id+j];
      for(int j=0; j<3; j++)
        surface.Udot[i][j] = temp_buffer[bufsize*id+3+j];
      if(algNum==6) {//A6
        for(int j=0; j<3; j++)
          surface.X[i][j] += 0.5*dt*surface.Udot[i][j];
      }
    }
  }

  
  //broadcast surface and Udot
  MPI_Bcast((double*)surface.X.data(), 3*nNodes, MPI_DOUBLE, 0, fest_comm);
  MPI_Bcast((double*)surface.Udot.data(), 3*nNodes, MPI_DOUBLE, 0, fest_comm);

}

//---------------------------------------------------------------

void
AerosMessenger::SendFESTSuggestedTimeStep(double dtf0)
{
  if(fest_rank == 0) {

    vector<MPI_Request> send_requests;

    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(&dtf0, 1, MPI_DOUBLE, proc, SUGGEST_DT_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  }
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeStepping()
{
  if(algNum == 6) //A6
    CommunicateBeforeTimeSteppingForA6();
  else if(algNum == 22)
    CommunicateBeforeTimeSteppingForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchange()
{
  double refer_time = walltime();

  if(algNum == 6) //A6
    FirstExchangeForA6();
  else if(algNum == 22)
    FirstExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }

  if(verbose>0)
    print("- Total waiting time: %.4e s.\n", walltime() - refer_time);
}

//---------------------------------------------------------------

void
AerosMessenger::Exchange()
{
  double refer_time = walltime();

  if(algNum == 6) //A6
    ExchangeForA6();
  else if(algNum == 22)
    ExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }

  if(verbose>0)
    print("- Total waiting time: %.4e s.\n", walltime() - refer_time);
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchange()
{
  double refer_time = walltime();

  if(algNum == 6) //A6
    FinalExchangeForA6();
  else if(algNum == 22)
    FinalExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }

  if(verbose>0)
    print("- Total waiting time: %.4e s.\n", walltime() - refer_time);
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeSteppingForA6()
{
  dt *= 0.5;
  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchangeForA6()
{
  //nothing special
  ExchangeForA6(); 
}

//---------------------------------------------------------------

void
AerosMessenger::ExchangeForA6()
{
  SendForce();  

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchangeForA6()
{
  SendForce();
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeSteppingForC0()
{
  GetDisplacementAndVelocity();

  SendForce();

  GetInfo(); //get dt, tmax

  dt *= 0.5;
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchangeForC0()
{
  GetInfo(); //get dt, tmax

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::ExchangeForC0()
{
  SendForce();  

  GetInfo();

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchangeForC0()
{
  SendForce();
  dt = 0.0;
}

//---------------------------------------------------------------















//---------------------------------------------------------------

