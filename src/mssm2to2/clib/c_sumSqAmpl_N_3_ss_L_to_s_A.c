#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_ss_L_to_s_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_ss_L_to_s_A(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0008 = (-0.333333333333333)*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = pow(m_N_3, 2);
    const ccomplex_t IT_0012 = pow(m_ss_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 0.666666666666667*IT_0007;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (-2)*IT_0015 + -IT_0020;
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = conj(N_B3)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W3)*e_em;
    const ccomplex_t IT_0028 = IT_0004*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + (-3)*IT_0029);
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = IT_0008*IT_0031;
    const ccomplex_t IT_0033 = IT_0013*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0017*IT_0031;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0016*IT_0036;
    const ccomplex_t IT_0038 = 2*IT_0034 + IT_0037;
    const ccomplex_t IT_0039 = (-6)*IT_0010;
    const ccomplex_t IT_0040 = 0.333333333333333*IT_0007;
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0036*IT_0040;
    const ccomplex_t IT_0044 = (-6)*s_13;
    const ccomplex_t IT_0045 = s_34*IT_0044;
    const ccomplex_t IT_0046 = m_s*s_14*m_N_3;
    const ccomplex_t IT_0047 = -IT_0015;
    const ccomplex_t IT_0048 = s_34*IT_0011;
    const ccomplex_t IT_0049 = 6*IT_0048;
    const ccomplex_t IT_0050 = s_13*s_14;
    const ccomplex_t IT_0051 = (-12)*IT_0050;
    const ccomplex_t IT_0052 = IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = s_14*IT_0044;
    const ccomplex_t IT_0054 = IT_0011*IT_0044;
    const ccomplex_t IT_0055 = pow(s_13, 2);
    const ccomplex_t IT_0056 = (-6)*IT_0055;
    const ccomplex_t IT_0057 = -IT_0037;
    const ccomplex_t IT_0058 = m_s*m_N_3;
    const ccomplex_t IT_0059 = (-6)*s_14;
    const ccomplex_t IT_0060 = (-6)*s_34;
    const ccomplex_t IT_0061 = (-6)*IT_0011;
    const ccomplex_t IT_0062 = (conj(IT_0021)*IT_0038 + IT_0021*conj(IT_0038))
      *IT_0039 + IT_0044*(conj(IT_0020)*IT_0038 + IT_0020*conj(IT_0038) + conj
      (IT_0021)*IT_0057 + IT_0021*conj(IT_0057)) + (conj(IT_0020)*IT_0043 +
       IT_0020*conj(IT_0043) + conj(IT_0042)*IT_0057 + IT_0042*conj(IT_0057))
      *IT_0059 + (conj(IT_0038)*IT_0042 + IT_0038*conj(IT_0042) + conj(IT_0021)
      *IT_0043 + IT_0021*conj(IT_0043))*IT_0060 + (conj(IT_0020)*IT_0057 +
       IT_0020*conj(IT_0057))*IT_0061;
    const ccomplex_t IT_0063 = m_s*s_34*m_N_3;
    const ccomplex_t IT_0064 = (-6)*conj(IT_0021);
    const ccomplex_t IT_0065 = s_14*IT_0010;
    const ccomplex_t IT_0066 = (-6)*conj(IT_0038);
    const ccomplex_t IT_0067 = s_14*s_34;
    const ccomplex_t IT_0068 = (-6)*IT_0038;
    const ccomplex_t IT_0069 = (-0.166666666666667)*IT_0065;
    const ccomplex_t IT_0070 = (-6)*IT_0021;
    const ccomplex_t IT_0071 = (-0.166666666666667)*IT_0063;
    const ccomplex_t IT_0072 = s_13*(IT_0021*conj(IT_0021) + IT_0038*conj
      (IT_0038))*IT_0039 + (conj(IT_0021)*IT_0042 + IT_0021*conj(IT_0042) + conj
      (IT_0038)*IT_0043 + IT_0038*conj(IT_0043))*IT_0045 + (-6)*conj(IT_0020)*
      (IT_0034*IT_0046 + (-0.166666666666667)*IT_0047*IT_0052 + (
      -0.166666666666667)*IT_0042*IT_0053 + (-0.166666666666667)*IT_0020*IT_0054
       + (-0.166666666666667)*IT_0021*IT_0056) + IT_0020*(conj(IT_0042)*IT_0053 
      + conj(IT_0021)*IT_0056) + (conj(IT_0043)*IT_0053 + conj(IT_0038)*IT_0056)
      *IT_0057 + (-6)*(IT_0046*IT_0047 + (-0.166666666666667)*IT_0034*IT_0052 + 
      (-0.166666666666667)*IT_0043*IT_0053 + (-0.166666666666667)*IT_0038
      *IT_0056 + (-0.166666666666667)*IT_0054*IT_0057)*conj(IT_0057) + IT_0058
      *IT_0062 + IT_0047*(IT_0064*IT_0065 + IT_0063*IT_0066) + IT_0034*(IT_0063
      *IT_0064 + IT_0065*IT_0066) + (-6)*conj(IT_0047)*((-0.166666666666667)
      *IT_0020*IT_0052 + IT_0046*IT_0057 + (-4)*IT_0047*IT_0067 + IT_0069
      *IT_0070 + IT_0068*IT_0071) + (-6)*conj(IT_0034)*(IT_0020*IT_0046 + (
      -0.166666666666667)*IT_0052*IT_0057 + (-4)*IT_0034*IT_0067 + IT_0068
      *IT_0069 + IT_0070*IT_0071);
    return create_ccomplex_return(IT_0072);
}

