#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_2_to_anti_b_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_2_to_anti_b_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = pow(m_b, 2);
    const ccomplex_t IT_0002 = pow(m_C_2, 2);
    const ccomplex_t IT_0003 = pow(m_st_2, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_12 + IT_0001 + -IT_0002 + -IT_0003 
      + -reg_prop, -1);
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_t*V_tb*V_u2*e_em*IT_0012*conj(U_st_11);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0017 + IT_0022;
    const ccomplex_t IT_0024 = IT_0006*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (-2)*IT_0000;
    const ccomplex_t IT_0028 = IT_0017*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0027;
    const ccomplex_t IT_0031 = IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = (-2)*IT_0024 + IT_0032;
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0010*IT_0012*IT_0035*conj(U_st_01);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0025*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0027*IT_0039;
    const ccomplex_t IT_0041 = IT_0000*IT_0037;
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0040 + (-2)*IT_0043;
    const ccomplex_t IT_0045 = (-8)*IT_0001;
    const ccomplex_t IT_0046 = -IT_0000;
    const ccomplex_t IT_0047 = -IT_0039*(IT_0027 + -IT_0046);
    const ccomplex_t IT_0048 = IT_0017*IT_0046;
    const ccomplex_t IT_0049 = IT_0021*IT_0046;
    const ccomplex_t IT_0050 = IT_0026*(IT_0028 + -IT_0030 + -IT_0048 +
       IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = (-8)*s_13;
    const ccomplex_t IT_0053 = s_34*IT_0052;
    const ccomplex_t IT_0054 = IT_0026*(IT_0028 + -IT_0030);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0043;
    const ccomplex_t IT_0058 = m_b*s_14*m_C_2;
    const ccomplex_t IT_0059 = IT_0006*(IT_0017 + -IT_0021);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = s_34*IT_0002;
    const ccomplex_t IT_0062 = 8*IT_0061;
    const ccomplex_t IT_0063 = s_13*s_14;
    const ccomplex_t IT_0064 = (-16)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = s_14*IT_0052;
    const ccomplex_t IT_0067 = IT_0002*IT_0052;
    const ccomplex_t IT_0068 = pow(s_13, 2);
    const ccomplex_t IT_0069 = (-8)*IT_0068;
    const ccomplex_t IT_0070 = m_b*m_C_2;
    const ccomplex_t IT_0071 = (-8)*s_14;
    const ccomplex_t IT_0072 = (-8)*s_34;
    const ccomplex_t IT_0073 = (-8)*IT_0002;
    const ccomplex_t IT_0074 = (conj(IT_0033)*IT_0044 + IT_0033*conj(IT_0044))
      *IT_0045 + IT_0052*(conj(IT_0033)*IT_0040 + IT_0033*conj(IT_0040) + conj
      (IT_0044)*IT_0056 + IT_0044*conj(IT_0056)) + (conj(IT_0040)*IT_0051 +
       IT_0040*conj(IT_0051) + conj(IT_0047)*IT_0056 + IT_0047*conj(IT_0056))
      *IT_0071 + (conj(IT_0033)*IT_0047 + IT_0033*conj(IT_0047) + conj(IT_0044)
      *IT_0051 + IT_0044*conj(IT_0051))*IT_0072 + (conj(IT_0040)*IT_0056 +
       IT_0040*conj(IT_0056))*IT_0073;
    const ccomplex_t IT_0075 = m_b*s_34*m_C_2;
    const ccomplex_t IT_0076 = (-8)*conj(IT_0033);
    const ccomplex_t IT_0077 = s_14*IT_0001;
    const ccomplex_t IT_0078 = (-8)*conj(IT_0044);
    const ccomplex_t IT_0079 = s_14*s_34;
    const ccomplex_t IT_0080 = (-8)*IT_0044;
    const ccomplex_t IT_0081 = (-0.125)*IT_0077;
    const ccomplex_t IT_0082 = (-8)*IT_0033;
    const ccomplex_t IT_0083 = (-0.125)*IT_0075;
    const ccomplex_t IT_0084 = s_13*(IT_0033*conj(IT_0033) + IT_0044*conj
      (IT_0044))*IT_0045 + (conj(IT_0044)*IT_0047 + IT_0044*conj(IT_0047) + conj
      (IT_0033)*IT_0051 + IT_0033*conj(IT_0051))*IT_0053 + (-8)*conj(IT_0056)*
      (IT_0057*IT_0058 + (-0.125)*IT_0060*IT_0065 + (-0.125)*IT_0051*IT_0066 + (
      -0.125)*IT_0056*IT_0067 + (-0.125)*IT_0033*IT_0069) + IT_0056*(conj
      (IT_0051)*IT_0066 + conj(IT_0033)*IT_0069) + (-8)*conj(IT_0040)*(IT_0058
      *IT_0060 + (-0.125)*IT_0057*IT_0065 + (-0.125)*IT_0047*IT_0066 + (-0.125)
      *IT_0040*IT_0067 + (-0.125)*IT_0044*IT_0069) + IT_0040*(conj(IT_0047)
      *IT_0066 + conj(IT_0044)*IT_0069) + IT_0070*IT_0074 + IT_0060*(IT_0076
      *IT_0077 + IT_0075*IT_0078) + IT_0057*(IT_0075*IT_0076 + IT_0077*IT_0078) 
      + (-8)*conj(IT_0060)*(IT_0040*IT_0058 + (-0.125)*IT_0056*IT_0065 + (-4)
      *IT_0060*IT_0079 + IT_0081*IT_0082 + IT_0080*IT_0083) + (-8)*conj(IT_0057)
      *(IT_0056*IT_0058 + (-0.125)*IT_0040*IT_0065 + (-4)*IT_0057*IT_0079 +
       IT_0080*IT_0081 + IT_0082*IT_0083);
    return create_ccomplex_return(IT_0084);
}

