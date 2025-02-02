#include<ExactNodesOperator.h>
#include<InitialNodesOperator.h>
#include<UpdatedNodesOperator.h>
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
		      comm(comm_), file_handler(iod_meta, comm),
		      npo(nullptr)
{

  switch(iod_spatial.type) {

    case SpatialInterpolationData::EXACT : {
      npo = new ExactNodesOperator(iod_spatial, comm);
      break;
    }
    case SpatialInterpolationData::INITIAL : {
      npo = new InitialNodesOperator(iod_spatial, comm);
      break;
    }
    case SpatialInterpolationData::UPDATED : {
      npo = new UpdatedNodesOperator(iod_spatial, comm);
      break;
    }
    default : {

      print_error("*** Error: Something went wrong while setting up the "
                  "NodalProjectionOperator.\n");
      exit_mpi();

    }
  } 

}

//------------------------------------------------------------

void DynamicLoadOperator::Destroy()
{
  if(npo) {
    npo->Destroy();
    delete npo;
  }
}

//------------------------------------------------------------

//! Copied from DynamicLoadCalculator::InterpolateInTime
void
DynamicLoadOperator::InterpolateInTime(double t1, double* input1, double t2, double* input2,
                                        double t, double* output, int size)
{
  assert(t2>t1);
  double c1 = (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//------------------------------------------------------------

