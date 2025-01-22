/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _AEROS_MESSENGER_H_
#define _AEROS_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>
#include<TriangulatedSurface.h>

/*************************************************************
* Class AerosMessenger is responsible for communicating with
* the AERO-S solver, e.g. for fluid-structure interaction
* simulations.
* This class is like (but not the same as) a combination of 
* DynamicNodalTransfer, EmbeddedStructure, StructExc, 
* MatchNodes, and EmbeddedMeshMotionHandler in AERO-F.
************************************************************/

class AerosMessenger {

  AerosCouplingData &iod_aeros;

  MPI_Comm &fest_comm; //!< This is the FEST communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of FEST and AERO-S

  int fest_rank, fest_size;

  int numAerosProcs;
  int (*numStrNodes)[2];  //!< numStrNodes[AEROS-proc-num][0]: the num of nodes; [1]: index of first node
  std::vector<int> pack2local; //!< [id. in AERO-S package]  -> [node number in surface] 
  std::vector<int> local2pack; //!< [node number in surface] -> [id. in AERO-S package]
  int bufsize; //!< number of DOFs per node (6)

  int nNodes, totalNodes;
  int nElems, totalElems;
  int elemType; //!< 3 or 4
  TriangulatedSurface &surface; //!< the embedded surface
  std::vector<Vec3D> &F;

  double dt, tmax; //!< dt, tmax suggested by AERO-S. May be different from those in AERO-S (Staggering)
  int algNum; //!< algo number received from AERO-S
  int structureSubcycling;

  std::vector<double> temp_buffer;

public:

  AerosMessenger(AerosCouplingData &iod_aeros_, MPI_Comm &fest_comm_, MPI_Comm &joint_comm_, 
                 TriangulatedSurface &surf_, std::vector<Vec3D> &F_);
  ~AerosMessenger();
  void Destroy();

  double GetTimeStepSize() {return dt;}
  double GetMaxTime() {return tmax;}

  //! Functions for structure/AERO-S subcycling (not supported at present)
  int StructSubcycling() {return structureSubcycling;}
  void SendFESTSuggestedTimeStep(double dtf0);

  //! Exchange data w/ AERO-S (called before the first time step)
  void CommunicateBeforeTimeStepping();

  //! Exchange data w/ AERO-S (called at the first time step)
  void FirstExchange();

  //! Exchange data w/ AERO-S (called at every time step except first and last)
  void Exchange();

  //! Exchange data w/ AERO-S (called at the last time step)
  void FinalExchange();

protected:

  //! functions called by the constructor
  void GetEmbeddedWetSurfaceInfo(int &eType, bool &crack, int &nStNodes, int &nStElems);
  void GetEmbeddedWetSurface(int nNodes, Vec3D *nodes, int nElems, int *elems, int eType);
  int  SplitQuads(int *quads, int nStElems, std::vector<Int3> &Tria);
  void Negotiate();
  void GetInfo();
  void GetInitialPhantomNodes(int newNodes, std::vector<Vec3D>& xyz, int nNodes);
  int  GetStructSubcyclingInfo();

  //! functions call by data-exchange functions
  void CommunicateBeforeTimeSteppingForA6();
  void FirstExchangeForA6();
  void ExchangeForA6();
  void FinalExchangeForA6();
  void CommunicateBeforeTimeSteppingForC0();
  void FirstExchangeForC0();
  void ExchangeForC0();
  void FinalExchangeForC0();

  //! functions call by the above exchange functions
  void GetDisplacementAndVelocity();
  void SendForce();

};





#endif
