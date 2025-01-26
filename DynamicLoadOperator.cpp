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
		      comm(comm_), file_handler(iod_meta, comm)
{

  int num_points = iod_meta.numPoints;
  // resize solution container
  proxi_solutions.resize(num_points, NULL);
  proxi_surfaces.resize(num_points, NULL);

}

//------------------------------------------------------------

void DynamicLoadOperator::Destroy()
{
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

void DynamicLoadOperator::BuildSurfacesToSurfaceMap(TriangulatedSurface &)
{
  // Default is to do nothing.
}

//------------------------------------------------------------

void DynamicLoadOperator::LoadExistingSurfaces()
{
  int num_points = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i) {
  
    // create new "empty" triangulated surfaces 
    proxi_surfaces[i] = new TriangulatedSurface();

    // first get the mesh file name from file handler
    string &filename = file_handler.GetMeshFileForProxim(i);

    // read this mesh file and update the triangulated surface
    file_handler.ReadMeshFile(filename, (*proxi_surfaces[i]).X, 
                              (*proxi_surfaces[i]).elems);

    (*proxi_surfaces[i]).X0 = (*proxi_surfaces[i]).X;
    (*proxi_surfaces[i]).BuildConnectivities();
    (*proxi_surfaces[i]).CalculateNormalsAndAreas();

  }
}

//------------------------------------------------------------

void DynamicLoadOperator::LoadExistingSolutions()
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
DynamicLoadOperator::InterpolateInMetaSpace(TriangulatedSurface &surface, vector<vector<Vec3D>> &solutions, 
                                            vector<Vec3D> &force, vector<Vec3D> *force_over_area)
{
  // default is to do nothing.
}

//------------------------------------------------------------

void
DynamicLoadOperator::InterpolateInSpace(std::vector<Vec3D>& X, int active_nodes, int dim, double* output)
{
  // default is to do nothing.
}

//------------------------------------------------------------

void
DynamicLoadOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                        double t, double* output, int size)
{
  // default is to do nothing.
}

//------------------------------------------------------------

