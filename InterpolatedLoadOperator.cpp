#include<InterpolatedLoadOperator.h>
#include<MathTools/rbf_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

InterpolatedLoadOperator::InterpolatedLoadOperator(IoData &iod_, MPI_Comm &comm_)
                        : DynamicLoadOperator(iod_, comm_)
{
   
  int num_points = iod_meta.numPoints;
  // resize solution container
  proxi_solutions.resize(num_points, nullptr);
  proxi_surfaces.resize(num_points, nullptr);

}

//------------------------------------------------------------

InterpolatedLoadOperator::~InterpolatedLoadOperator() {
 
  if(!proxi_solutions.empty())
    for(int i=0; i<(int)proxi_solutions.size(); ++i) {
      if(proxi_solutions[i]) delete proxi_solutions[i];
    }

  if(!proxi_surfaces.empty())
    for(int i=0; i<(int)proxi_surfaces.size(); ++i) {
      if(proxi_surfaces[i]) delete proxi_surfaces[i];
    }

}

//------------------------------------------------------------

void InterpolatedLoadOperator::LoadExistingSurfaces()
{

  if(!npo->NeedSurface()) return;

  int num_points = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i) {
  
    // create new "empty" triangulated surface
    proxi_surfaces[i] = new TriangulatedSurface();
    string &filename  = file_handler.GetMeshFileForProxim(i);
    file_handler.ReadMeshFile(filename, proxi_surfaces[i]->X, 
                              proxi_surfaces[i]->elems);

    proxi_surfaces[i]->X0 = proxi_surfaces[i]->X;
    proxi_surfaces[i]->BuildConnectivities();
    proxi_surfaces[i]->CalculateNormalsAndAreas();

  }

}

//------------------------------------------------------------

void InterpolatedLoadOperator::LoadExistingSolutions()
{
  int num_points = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i) {
  
    // create new "empty" solution containers
    proxi_solutions[i] = new SolutionData3D();

    // first get the solution file name from file handler
    string &filename = file_handler.GetSolnFileForProxim(i);

    // read this solution file and update the solution container
    file_handler.ReadSolutionFile(filename, *proxi_solutions[i]);

  } 
}

//------------------------------------------------------------

void
InterpolatedLoadOperator::SetupProjectionMap(TriangulatedSurface &surface)
{

  assert(npo); //cannot be null
  int num_points   = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i)
    npo->SetupProjectionMap(surface, *proxi_surfaces[i]);

}

//------------------------------------------------------------

void
InterpolatedLoadOperator::ComputeForces(TriangulatedSurface& surface, std::vector<Vec3D> &force,
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

  vector<vector<Vec3D>> proxi_forces(num_points, 
                                     vector<Vec3D>(active_nodes, Vec3D(0.0)));

  for(int i=0; i<num_points; ++i) {

    // find the time interval
    double tk, tkp;
    auto time_bounds = proxi_solutions[i]->GetTimeBounds();

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

      auto bracket = proxi_solutions[i]->GetTimeBracket(t);
      tk  = bracket[0];
      tkp = bracket[1];

    }

    // Get solutions at tk and tkp
    vector<Vec3D> &Sk  = proxi_solutions[i]->GetSolutionAtTime(tk);
    vector<Vec3D> &Skp = proxi_solutions[i]->GetSolutionAtTime(tkp);
    vector<Vec3D> &S   = proxi_forces[i];

    npo->ProjectToTargetSurface(surface, *proxi_surfaces[i], Sk);
    npo->ProjectToTargetSurface(surface, *proxi_surfaces[i], Skp);
    InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                      (double*)S.data(), 3*active_nodes);

  }

  InterpolateInMetaSpace(surface, proxi_forces, force, force_over_area);


}

//------------------------------------------------------------

void 
InterpolatedLoadOperator::InterpolateInMetaSpace(TriangulatedSurface &surface, vector<vector<Vec3D>> &solutions, 
                                                 vector<Vec3D> &force, vector<Vec3D> *force_over_area) 
{

  int active_nodes = surface.active_nodes;
  int num_points   = iod_meta.numPoints;

  // Get parameters for target and proximal points --- sorted by distance.
  vector<double> &targ              = file_handler.GetParametersForTarget();
  std::vector<vector<double>> &prox = file_handler.GetParametersForAllProxim();

  int var_dim = targ.size();
  double rmin = 0, rmax = 0;
  for(int i=0; i<var_dim; ++i) {
    rmin += (prox[           0][i] - targ[i])*(prox[           0][i] - targ[i]);
    rmax += (prox[num_points-1][i] - targ[i])*(prox[num_points-1][i] - targ[i]);
  }

  rmin = std::sqrt(rmin);
  rmax = std::sqrt(rmax);

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

  //fprintf(stdout, "rank[%d] block_size %d, start_index %d.\n", mpi_rank, my_block_size, my_start_index);

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

  // compute forces
  for(int index=my_start_index; index<my_block_size; ++index) {

    // calculate nodal area and current normal.
    Vec3D patch(0.0);
    auto elems  = surface.node2elem[index];
    int  nelems = elems.size();
    for(auto it=elems.begin(); it!=elems.end(); ++it)
      patch += surface.elemArea[*it]*surface.elemNorm[*it];
    
    patch /= nelems;

    double area   = patch.norm();
    Vec3D  normal = patch/area;

    // prepare for interpolation
    double xd[var_dim*num_points];
    for(int i=0; i<num_points; ++i)
      for(int j=0; j<var_dim; ++j)
        xd[var_dim*i + j] = prox[i][j];


    double r0; //smaller than maximum separation, larger than typical separation
    r0 = rmin + rmax;

    double fd[num_points];

    //interpolate
    for(int i=0; i<num_points; ++i)
      fd[i] = solutions[i][index].norm();

    vector<double> weight(num_points, -1.0);
    vector<double> interp(num_points, -1.0);

    MathTools::rbf_weight(var_dim, num_points, xd, r0, phi, fd, weight.data());
    MathTools::rbf_interp(var_dim, num_points, xd, r0, phi, weight.data(), 1,
                          targ.data(), interp.data());

    force[index] = interp[0]*area*normal;   
    if(force_over_area)
      (*force_over_area)[index] = interp[0]*normal;

  }	  


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
