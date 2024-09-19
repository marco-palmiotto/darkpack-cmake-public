#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_ss_R_to_s_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_ss_R_to_s_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0007;
    const ccomplex_t IT_0009 = (-0.333333333333333)*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (-2)*IT_0012 + -IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = sin(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0019*IT_0021*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0013*IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0025;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0028 + 2*IT_0031;
    const ccomplex_t IT_0033 = (-6)*IT_0000;
    const ccomplex_t IT_0034 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0027*IT_0034;
    const ccomplex_t IT_0038 = (-6)*s_13;
    const ccomplex_t IT_0039 = s_34*IT_0038;
    const ccomplex_t IT_0040 = m_s*s_14*m_N_1;
    const ccomplex_t IT_0041 = -IT_0012;
    const ccomplex_t IT_0042 = s_34*IT_0001;
    const ccomplex_t IT_0043 = 6*IT_0042;
    const ccomplex_t IT_0044 = s_13*s_14;
    const ccomplex_t IT_0045 = (-12)*IT_0044;
    const ccomplex_t IT_0046 = IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = s_14*IT_0038;
    const ccomplex_t IT_0048 = IT_0001*IT_0038;
    const ccomplex_t IT_0049 = pow(s_13, 2);
    const ccomplex_t IT_0050 = (-6)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0028;
    const ccomplex_t IT_0052 = m_s*m_N_1;
    const ccomplex_t IT_0053 = (-6)*s_14;
    const ccomplex_t IT_0054 = (-6)*s_34;
    const ccomplex_t IT_0055 = (-6)*IT_0001;
    const ccomplex_t IT_0056 = (conj(IT_0018)*IT_0032 + IT_0018*conj(IT_0032))
      *IT_0033 + IT_0038*(conj(IT_0017)*IT_0032 + IT_0017*conj(IT_0032) + conj
      (IT_0018)*IT_0051 + IT_0018*conj(IT_0051)) + (conj(IT_0017)*IT_0037 +
       IT_0017*conj(IT_0037) + conj(IT_0036)*IT_0051 + IT_0036*conj(IT_0051))
      *IT_0053 + (conj(IT_0032)*IT_0036 + IT_0032*conj(IT_0036) + conj(IT_0018)
      *IT_0037 + IT_0018*conj(IT_0037))*IT_0054 + (conj(IT_0017)*IT_0051 +
       IT_0017*conj(IT_0051))*IT_0055;
    const ccomplex_t IT_0057 = m_s*s_34*m_N_1;
    const ccomplex_t IT_0058 = (-6)*conj(IT_0018);
    const ccomplex_t IT_0059 = s_14*IT_0000;
    const ccomplex_t IT_0060 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0061 = s_14*s_34;
    const ccomplex_t IT_0062 = (-6)*IT_0032;
    const ccomplex_t IT_0063 = (-0.166666666666667)*IT_0059;
    const ccomplex_t IT_0064 = (-6)*IT_0018;
    const ccomplex_t IT_0065 = (-0.166666666666667)*IT_0057;
    const ccomplex_t IT_0066 = s_13*(IT_0018*conj(IT_0018) + IT_0032*conj
      (IT_0032))*IT_0033 + (conj(IT_0018)*IT_0036 + IT_0018*conj(IT_0036) + conj
      (IT_0032)*IT_0037 + IT_0032*conj(IT_0037))*IT_0039 + (-6)*conj(IT_0017)*
      (IT_0031*IT_0040 + (-0.166666666666667)*IT_0041*IT_0046 + (
      -0.166666666666667)*IT_0036*IT_0047 + (-0.166666666666667)*IT_0017*IT_0048
       + (-0.166666666666667)*IT_0018*IT_0050) + IT_0017*(conj(IT_0036)*IT_0047 
      + conj(IT_0018)*IT_0050) + (conj(IT_0037)*IT_0047 + conj(IT_0032)*IT_0050)
      *IT_0051 + (-6)*(IT_0040*IT_0041 + (-0.166666666666667)*IT_0031*IT_0046 + 
      (-0.166666666666667)*IT_0037*IT_0047 + (-0.166666666666667)*IT_0032
      *IT_0050 + (-0.166666666666667)*IT_0048*IT_0051)*conj(IT_0051) + IT_0052
      *IT_0056 + IT_0041*(IT_0058*IT_0059 + IT_0057*IT_0060) + IT_0031*(IT_0057
      *IT_0058 + IT_0059*IT_0060) + (-6)*conj(IT_0041)*((-0.166666666666667)
      *IT_0017*IT_0046 + IT_0040*IT_0051 + (-4)*IT_0041*IT_0061 + IT_0063
      *IT_0064 + IT_0062*IT_0065) + (-6)*conj(IT_0031)*(IT_0017*IT_0040 + (
      -0.166666666666667)*IT_0046*IT_0051 + (-4)*IT_0031*IT_0061 + IT_0062
      *IT_0063 + IT_0064*IT_0065);
    return create_ccomplex_return(IT_0066);
}

