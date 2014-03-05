#include "interpola_vpr.h"
#include <vpr_par.h>

namespace cumbac {

InterpolaVPR::InterpolaVPR()
    : B(NODATAVPR), E(NODATAVPR), G(NODATAVPR), C(NODATAVPR), F(NODATAVPR),
      chisqfin(100), rmsefin(100)
{
}

InterpolaVPR::~InterpolaVPR()
{
}

}
