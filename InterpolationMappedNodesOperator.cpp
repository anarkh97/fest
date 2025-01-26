#include<InterpolationMappedNodesOperator.h>

using std::vector;

extern int verbose;

//------------------------------------------------------------------------------------------------

InterpolationMappedNodesOperator::InterpolationMappedNodesOperator(IoData& iod_, MPI_Comm& comm_)
                                : DynamicLoadOperator(iod_, comm_)
{

}

//------------------------------------------------------------------------------------------------

void
InterpolationMappedNodesOperator::ComputeForces(TriangulatedSurface& surface, vector<Vec3D>& force,
                                                vector<Vec3D>* force_over_area, double t)
{

}

//------------------------------------------------------------------------------------------------

void
InterpolationMappedNodesOperator::InterpolateInMetaSpace(TriangulatedSurface &surface, vector<vector<Vec3D>> &solutions, 
                                                         vector<Vec3D> &force, vector<Vec3D> *force_over_area)
{
  //
}

//------------------------------------------------------------------------------------------------

void
InterpolationMappedNodesOperator::InterpolateInSpace(std::vector<Vec3D>& X, int active_nodes, int dim, double* output)
{
  //
}

//------------------------------------------------------------------------------------------------

void
InterpolationMappedNodesOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                                    double t, double* output, int size)
{
  //
}

//------------------------------------------------------------------------------------------------

//! Code re-used from M2C's DynamicLoadCalculator::BuildKDTree
void InterpolationMappedNodesOperator::BuildKDTree(vector<Vec3D> &Xs, K3DTree* tree, vector<PointIn3D> &p)
{

  if(tree)
    delete tree;

  int N = Xs.size();
  p.resize(N);

  for(int i=0; i<N; ++i) 
    p[i] = PointIn3D(i, Xs[i]);
  
  tree = new K3DTree(N, p.data()); 

}

//------------------------------------------------------------------------------------------------

void
InterpolationMappedNodesOperator::BuildSurfacesToSurfaceMap(TriangulatedSurface &surface)
{

/* AN: Code is incomplete...	
  int active_nodes = surface.active_nodes;
  int num_points   = iod_meta.numPoints;

  assert(active_nodes == (int)surface.X.size()); //cracking not supported.

  // read mesh files for each point --- ordered based on distance from target
  vector<TriangulatedSurface>     proxi_surfaces(num_points);
  vector<vector<PointIn3D>>       proxi_tree_data(num_points); 
  vector<std::shared_ptr<K3DTree>> proxi_trees(num_points, nullptr); // automatically deleted.

  for(int i=0; i<num_points; ++i) {
    std::string filename = proxi_surface_files[i];
    ReadMeshFile(filename.c_str(), proxi_surfaces[i].X, proxi_surfaces[i].elems);

    proxi_surfaces[i].BuildConnectivities();
    proxi_surfaces[i].CalculateNormalsAndAreas();

    // scale the co-ordinates to [0, 1].
    vector<Vec3D> UV(proxi_surfaces[i].X.size(), Vec3D(0.0));
    MathTools::affine_transformation(proxi_surfaces[i].X.size(), 3, 
                                     (double*)proxi_surfaces[i].X.data(), (double*)UV.data());

    // update surface co-ordinates (X0 not used by TriangulatedSurface for normal and area
    // calculations, should be safe to update co-ordinates).
    proxi_surfaces[i].X  = UV;
    proxi_surfaces[i].X0 = proxi_surfaces[i].X;

    // build the tree based on scaled co-ordinates
    BuildKDTree(proxi_surfaces[i].X, proxi_trees[i].get(), proxi_tree_data[i]);

  }

  // prepare for mapping
  // scale the target surface to [0, 1] 
  vector<Vec3D> ScaleX(active_nodes, Vec3D(0.0));
  MathTools::affine_transformation(active_nodes, 3, (double*)surface.X.data(), (double*)ScaleX.data());

  // prepare for parallel computation
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

  // compute and store the barycentric weights for each node.
  PointIn3D candidates[num_points][30];
  for(int index=start_index; index<my_block_size; ++index) {

    // current node on target surface
    Vec3D &Xn(surface.X[index]);

    // allocate memory
    barycenter_map[index] = vector<double>(3*num_points, 0.0);

    for(int i=0; i<num_points; ++i) {

      assert(proxi_trees[i]); // should not be null
      int nFound = proxi_trees[i]->findCandidates(Xn, candidates[i], 30);

      if(nFound < 3) {
        fprintf(stderr, "\033[0;31m*** Error: Cannot find the required number "
		        "of nodes (3) on the surface %d provided in the metafile."
			"\n\033[0m", i+1);
	exit(-1);
      }

      // sort by distance
      vector<pair<double, int>> dist2node;
      for(int i=0; i<nFound; i++)
        dist2node.push_back(std::make_pair(
                            (candidates[i].x-pnode).norm(), candidates[i].id));

      if(nFound>3)
        sort(dist2node.begin(), dist2node.end());
      dist2node.resize(3);

      // prepare for barycenter weights
      
      // compute weights
      vector<double> weights(3);

      // update barycenter map
      barycenter_map[index][3*i]   = weights[0];
      barycenter_map[index][3*i+1] = weights[1];
      barycenter_map[index][3*i+2] = weights[2];

      // update node2nodes
    }

  }	  

  // communication
  // keys
  // values
*/
}

//------------------------------------------------------------

