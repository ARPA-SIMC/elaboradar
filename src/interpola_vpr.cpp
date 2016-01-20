#include "interpola_vpr.h"
//#include <vpr_par.h>  presente in interpola_vpr.h

namespace radarelab {

InterpolaVPR::InterpolaVPR()
    : B(NODATAVPR), E(NODATAVPR), G(NODATAVPR), C(NODATAVPR), F(NODATAVPR),
      chisqfin(100), rmsefin(100)
{
}

InterpolaVPR::~InterpolaVPR()
{
}

}
