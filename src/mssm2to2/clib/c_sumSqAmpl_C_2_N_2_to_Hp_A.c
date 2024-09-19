#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0000*IT_0020;
    const ccomplex_t IT_0022 = pow(m_Hp, 2);
    const ccomplex_t IT_0023 = pow(m_N_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = 2*IT_0000;
    const ccomplex_t IT_0029 = pow(m_C_2, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0023 + -IT_0029 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0028*IT_0032;
    const ccomplex_t IT_0034 = 2*IT_0026 + -IT_0033;
    const ccomplex_t IT_0035 = s_24*IT_0029;
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (-2)*IT_0039;
    const ccomplex_t IT_0041 = IT_0009*IT_0036;
    const ccomplex_t IT_0042 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0046 = IT_0037*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = IT_0044 + IT_0047;
    const ccomplex_t IT_0049 = IT_0040 + IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0000*IT_0051;
    const ccomplex_t IT_0053 = IT_0024*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0030*IT_0051;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0028*IT_0056;
    const ccomplex_t IT_0058 = (-2)*IT_0054 + IT_0057;
    const ccomplex_t IT_0059 = s_14*m_C_2*m_N_2;
    const ccomplex_t IT_0060 = s_24*m_C_2*m_N_2;
    const ccomplex_t IT_0061 = s_14*s_24;
    const ccomplex_t IT_0062 = -IT_0033;
    const ccomplex_t IT_0063 = s_12*s_24;
    const ccomplex_t IT_0064 = s_14*IT_0023;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = m_C_2*m_N_2;
    const ccomplex_t IT_0068 = 2*s_12;
    const ccomplex_t IT_0069 = IT_0000*IT_0032;
    const ccomplex_t IT_0070 = IT_0000*IT_0056;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = 2*s_14;
    const ccomplex_t IT_0073 = 2*s_24;
    const ccomplex_t IT_0074 = 2*IT_0023;
    const ccomplex_t IT_0075 = IT_0029*(conj(IT_0034)*IT_0058 + IT_0034*conj
      (IT_0058)) + 0.5*(conj(IT_0034)*IT_0057 + IT_0034*conj(IT_0057) + conj
      (IT_0058)*IT_0062 + IT_0058*conj(IT_0062))*IT_0068 + 0.5*(conj(IT_0058)
      *IT_0069 + IT_0058*conj(IT_0069) + conj(IT_0034)*IT_0071 + IT_0034*conj
      (IT_0071))*IT_0072 + 0.5*(conj(IT_0057)*IT_0069 + IT_0057*conj(IT_0069) +
       conj(IT_0062)*IT_0071 + IT_0062*conj(IT_0071))*IT_0073 + 0.5*(conj
      (IT_0057)*IT_0062 + IT_0057*conj(IT_0062))*IT_0074;
    const ccomplex_t IT_0076 = (-2)*s_12;
    const ccomplex_t IT_0077 = s_14*IT_0076;
    const ccomplex_t IT_0078 = s_24*IT_0076;
    const ccomplex_t IT_0079 = IT_0023*IT_0076;
    const ccomplex_t IT_0080 = pow(s_12, 2);
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = (-2)*conj(IT_0054)*(IT_0035*IT_0058 + -IT_0034
      *IT_0059 + (-4)*IT_0054*IT_0061 + -IT_0060*IT_0062 + 2*IT_0057*IT_0066) + 
      (-2)*conj(IT_0027)*(IT_0034*IT_0035 + -IT_0058*IT_0059 + -IT_0057*IT_0060 
      + (-4)*IT_0027*IT_0061 + 2*IT_0062*IT_0066) + 2*IT_0067*IT_0075 + conj
      (IT_0071)*(IT_0058*IT_0077 + IT_0057*IT_0078) + conj(IT_0069)*(IT_0034
      *IT_0077 + IT_0062*IT_0078) + 2*conj(IT_0062)*(IT_0054*IT_0060 + (-2)
      *IT_0027*IT_0066 + 0.5*IT_0069*IT_0078 + 0.5*IT_0062*IT_0079 + 0.5*IT_0034
      *IT_0081) + (-2)*conj(IT_0058)*(IT_0035*IT_0054 + -IT_0027*IT_0059 + (-0.5
      )*IT_0029*IT_0058*IT_0076 + (-0.5)*IT_0071*IT_0077 + (-0.5)*IT_0057
      *IT_0081) + 2*conj(IT_0057)*(IT_0027*IT_0060 + (-2)*IT_0054*IT_0066 + 0.5
      *IT_0071*IT_0078 + 0.5*IT_0057*IT_0079 + 0.5*IT_0058*IT_0081) + (-2)*conj
      (IT_0034)*(IT_0027*IT_0035 + -IT_0054*IT_0059 + (-0.5)*IT_0029*IT_0034
      *IT_0076 + (-0.5)*IT_0069*IT_0077 + (-0.5)*IT_0062*IT_0081);
    return create_ccomplex_return(IT_0082);
}

