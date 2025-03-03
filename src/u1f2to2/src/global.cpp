#include "global.h"
#include "libdiagonalization.h"
#include "u1f2to2.h"
#include "libcomplexop.h"

namespace u1f2to2 {


void updateSpectrum(param_t &params)
{
    updateMassExpressions(params);
}

void updateMassExpressions(param_t &params)
{
    params.m_chi_dm_1 = m_chi_dm_1(params).real();
    params.m_chi_dm_2 = m_chi_dm_2(params).real();
}

} // End of namespace u1f2to2

