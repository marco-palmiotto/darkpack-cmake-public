#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_A_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_A_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t IT_0000 = cpow(s_13 + (-0.5)*reg_prop + (0 + _Complex_I*(
      -0.5))*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0003 = 2*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cpow(s_23 + (-0.5)*reg_prop + (0 + _Complex_I*(
      -0.5))*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = IT_0002*IT_0003;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = -IT_0011;
    const ccomplex_t IT_0015 = IT_0004 + IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0012 + (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0009;
    const ccomplex_t IT_0019 = pow(e_em, 2);
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = tan(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, 2);
    const ccomplex_t IT_0024 = cpow(IT_0020 + IT_0020*IT_0023, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0019*IT_0021*IT_0024;
    const ccomplex_t IT_0026 = 2*IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0011;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0012;
    const ccomplex_t IT_0030 = cpow(IT_0002, 2);
    const ccomplex_t IT_0031 = IT_0008*IT_0030;
    const ccomplex_t IT_0032 = -IT_0030;
    const ccomplex_t IT_0033 = IT_0011 + IT_0032;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0031 + 0.5*IT_0034;
    const ccomplex_t IT_0036 = 0.25*s_12*(IT_0006 + -conj(IT_0006)) + IT_0026 
      + 0.25*s_14*(IT_0028 + -conj(IT_0028)) + 0.25*s_24*(IT_0029 + -conj
      (IT_0029)) + 0.25*s_34*(IT_0035 + -conj(IT_0035));
    const ccomplex_t IT_0037 = -s_23*IT_0026;
    const ccomplex_t IT_0038 = s_12*s_24;
    const ccomplex_t IT_0039 = s_13*s_24;
    const ccomplex_t IT_0040 = s_12*s_34;
    const ccomplex_t IT_0041 = s_13*s_34;
    const ccomplex_t IT_0042 = s_23*s_24;
    const ccomplex_t IT_0043 = s_23*s_34;
    const ccomplex_t IT_0044 = conj(IT_0010) + conj(IT_0017);
    const ccomplex_t IT_0045 = s_23*IT_0026;
    const ccomplex_t IT_0046 = pow(m_stau_2, 2);
    const ccomplex_t IT_0047 = s_23*IT_0018;
    const ccomplex_t IT_0048 = s_23*conj(IT_0018);
    const ccomplex_t IT_0049 = conj(IT_0018)*IT_0026 + s_34*(conj(IT_0010)
      *IT_0029 + IT_0010*conj(IT_0029)) + s_24*(conj(IT_0006)*IT_0028 + IT_0006
      *conj(IT_0028) + conj(IT_0018)*IT_0029 + IT_0018*conj(IT_0029)) + IT_0006*
      (s_13*conj(IT_0017) + s_12*conj(IT_0018) + conj(IT_0006)*IT_0046) +
       IT_0018*(s_12*conj(IT_0006) + -IT_0026 + conj(IT_0018)*IT_0046) + IT_0044
      *IT_0047 + IT_0010*IT_0048 + IT_0017*(s_13*conj(IT_0006) + IT_0048);
    const ccomplex_t IT_0050 = s_23*(s_12*(conj(IT_0006)*IT_0010 + IT_0006
      *conj(IT_0010)) + s_13*(conj(IT_0006)*IT_0013 + IT_0006*conj(IT_0013)) +
       s_23*(conj(IT_0010)*IT_0017 + IT_0010*conj(IT_0017) + conj(IT_0013)
      *IT_0018 + IT_0013*conj(IT_0018))) + (-4)*IT_0026*IT_0036 + (IT_0010 +
       IT_0017)*IT_0037 + conj(IT_0006)*(IT_0029*IT_0038 + IT_0035*IT_0039) +
       IT_0006*(conj(IT_0029)*IT_0038 + conj(IT_0035)*IT_0039) + conj(IT_0028)*
      (IT_0018*IT_0038 + IT_0017*IT_0039 + IT_0010*IT_0040 + IT_0013*IT_0041) +
       IT_0028*(conj(IT_0018)*IT_0038 + conj(IT_0017)*IT_0039 + conj(IT_0010)
      *IT_0040 + conj(IT_0013)*IT_0041) + (conj(IT_0017)*IT_0029 + IT_0017*conj
      (IT_0029) + conj(IT_0018)*IT_0035 + IT_0018*conj(IT_0035))*IT_0042 + (conj
      (IT_0013)*IT_0029 + IT_0013*conj(IT_0029) + conj(IT_0010)*IT_0035 +
       IT_0010*conj(IT_0035))*IT_0043 + IT_0044*IT_0045 + IT_0046*IT_0049;
    return create_ccomplex_return(IT_0050);
}

