#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "stdbool.h"
#include "common.h"

namespace scalar2to2 {

void setMu(const double mu)
{
    setmudim(mu * mu);
}

void setLambda2(const double lambda2)
{
    setlambda(lambda2);
}

void setUVDiv(const double x)
{
    setuvdiv(x);
}

} // End of namespace scalar2to2

