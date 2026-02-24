#include "common.h"
#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include "stdbool.h"
#include <cmath>

namespace mssm2to2
{

  void setMu(const double mu) { setmudim(mu * mu); }

  void setLambda2(const double lambda2) { setlambda(lambda2); }

  void setUVDiv(const double x) { setuvdiv(x); }

} // End of namespace mssm2to2
