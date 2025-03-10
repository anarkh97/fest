#include<ConstantLoadOperator.h>
#include<ClosestLoadOperator.h>
#include<InterpolatedLoadOperator.h>
#include<DynamicDriver.h>
#include<cstring>

using std::vector;
using std::shared_ptr;
using std::make_shared;

extern double start_time;
extern int verbose;

//------------------------------------------------------------

DynamicDriver::DynamicDriver(IoData &iod_, MPI_Comm &comm_) 
             : iod(iod_), comm(comm_), lagout(comm_, iod_.output)
{

  // Setup the interpolator.
  switch(iod.calculator.type) {

    case DynamicLoadCalculatorData::CONSTANT : {
      dlo = make_shared<ConstantLoadOperator>(iod, comm);
      break;
    }
    case DynamicLoadCalculatorData::CLOSEST : {
      dlo = make_shared<ClosestLoadOperator>(iod, comm);
      break;
    }
    case DynamicLoadCalculatorData::INTERP : {
      dlo = make_shared<InterpolatedLoadOperator>(iod, comm);
      break;
    }
    case DynamicLoadCalculatorData::NONE : {
      print_warning("- Warning: The Dynamic calculator type was not specified. "
                    "Assuming constant pressure calculator.\n");
      dlo = make_shared<ConstantLoadOperator>(iod, comm);
      break;
    }

  }

}

//------------------------------------------------------------

void DynamicDriver::Destroy()
{
  if(dlo.get()) {
    dlo->Destroy();
  }
}

//------------------------------------------------------------
