#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_Z_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_Z_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*g_s*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 0.666666666666667*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = pow(m_sd_R, 2);
    const ccomplex_t IT_0007 = pow(m_Z, -2);
    const ccomplex_t IT_0008 = pow(s_13, 2);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0006 + IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0015 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = pow(m_Z, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0017 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (-2)*IT_0012;
    const ccomplex_t IT_0022 = IT_0015*IT_0021;
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = 4*s_24;
    const ccomplex_t IT_0026 = conj(IT_0024)*IT_0025;
    const ccomplex_t IT_0027 = 4*IT_0006;
    const ccomplex_t IT_0028 = s_13*s_34*IT_0007;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = s_14 + IT_0029;
    const ccomplex_t IT_0031 = 4*IT_0004;
    const ccomplex_t IT_0032 = (-4)*IT_0004;
    const ccomplex_t IT_0033 = s_23*s_34*IT_0007;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = s_24 + IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0017 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0016*IT_0037;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = s_13*s_23*IT_0007;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = s_12 + IT_0041;
    const ccomplex_t IT_0043 = IT_0022*IT_0037;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = 4*s_34;
    const ccomplex_t IT_0046 = 4*s_23;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = IT_0027*IT_0043;
    const ccomplex_t IT_0049 = IT_0026*IT_0039 + conj(IT_0020)*(IT_0025
      *IT_0043 + IT_0044*IT_0045) + IT_0020*(IT_0025*conj(IT_0043) + conj
      (IT_0044)*IT_0045) + IT_0024*(IT_0032 + IT_0025*conj(IT_0039) + IT_0027
      *conj(IT_0043) + conj(IT_0044)*IT_0046) + conj(IT_0024)*(IT_0031 + IT_0047
       + IT_0048);
    const ccomplex_t IT_0050 = pow(s_23, 2);
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0006 + IT_0052;
    const ccomplex_t IT_0054 = (IT_0032 + IT_0025*conj(IT_0039))*IT_0043 +
       conj(IT_0039)*IT_0044*IT_0045 + 4*conj(IT_0044)*(IT_0017*IT_0044 + 0.25
      *IT_0039*IT_0045 + 0.25*IT_0043*IT_0046) + conj(IT_0043)*(IT_0031 +
       IT_0025*IT_0039 + IT_0047 + IT_0048);
    const ccomplex_t IT_0055 = (-12)*IT_0005 + IT_0011*(IT_0020*IT_0026 +
       IT_0024*(conj(IT_0020)*IT_0025 + conj(IT_0024)*IT_0027)) + IT_0030*(conj
      (IT_0020)*IT_0031 + IT_0020*IT_0032) + IT_0035*(IT_0032*IT_0039 + IT_0031
      *conj(IT_0039)) + IT_0042*IT_0049 + IT_0053*IT_0054;
    return create_ccomplex_return(IT_0055);
}

