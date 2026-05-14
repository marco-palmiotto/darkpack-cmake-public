#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_st_1_to_anti_b_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_st_1_to_anti_b_G(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = (-2)*IT_0000;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0003*IT_0005*IT_0007*conj(U_st_00);
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (IT_0001 + -IT_0002)*IT_0014;
    const ccomplex_t IT_0016 = IT_0002*IT_0014;
    const ccomplex_t IT_0017 = (-8)*s_13;
    const ccomplex_t IT_0018 = s_14*IT_0017;
    const ccomplex_t IT_0019 = IT_0000*IT_0009;
    const ccomplex_t IT_0020 = pow(m_st_1, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0020 
      + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = -IT_0016 + (-2)*IT_0023;
    const ccomplex_t IT_0025 = s_34*IT_0017;
    const ccomplex_t IT_0026 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0027 = IT_0007*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0002*IT_0029;
    const ccomplex_t IT_0031 = sin(beta);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = m_t*V_tb*V_u1*e_em*IT_0003*conj(U_st_10);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0002*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0041 = (IT_0030 + -IT_0039)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = m_b*m_C_1;
    const ccomplex_t IT_0044 = (-8)*s_14*IT_0043;
    const ccomplex_t IT_0045 = IT_0000*IT_0021;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = -IT_0029;
    const ccomplex_t IT_0048 = IT_0038 + IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = -IT_0039;
    const ccomplex_t IT_0051 = IT_0030 + IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0049 + IT_0052;
    const ccomplex_t IT_0054 = (-8)*s_34*IT_0043;
    const ccomplex_t IT_0055 = IT_0001*IT_0029;
    const ccomplex_t IT_0056 = IT_0001*IT_0038;
    const ccomplex_t IT_0057 = IT_0040*(IT_0030 + -IT_0039 + -IT_0055 +
       IT_0056);
    const ccomplex_t IT_0058 = -IT_0023;
    const ccomplex_t IT_0059 = s_34*IT_0011;
    const ccomplex_t IT_0060 = 8*IT_0059;
    const ccomplex_t IT_0061 = s_13*s_14;
    const ccomplex_t IT_0062 = (-16)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = IT_0011*IT_0017;
    const ccomplex_t IT_0065 = IT_0017*IT_0043;
    const ccomplex_t IT_0066 = pow(s_13, 2);
    const ccomplex_t IT_0067 = (-8)*IT_0066;
    const ccomplex_t IT_0068 = (-8)*IT_0011*IT_0043;
    const ccomplex_t IT_0069 = -(IT_0029 + -IT_0038)*IT_0046;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = m_b*s_14*m_C_1;
    const ccomplex_t IT_0072 = (-8)*IT_0071;
    const ccomplex_t IT_0073 = s_14*s_34;
    const ccomplex_t IT_0074 = 32*IT_0073;
    const ccomplex_t IT_0075 = s_14*IT_0010;
    const ccomplex_t IT_0076 = (-8)*IT_0075;
    const ccomplex_t IT_0077 = m_b*s_34*m_C_1;
    const ccomplex_t IT_0078 = (-8)*IT_0077;
    const ccomplex_t IT_0079 = (-8)*IT_0010;
    const ccomplex_t IT_0080 = IT_0043*IT_0079;
    const ccomplex_t IT_0081 = conj(IT_0015)*(IT_0016*IT_0018 + IT_0024
      *IT_0025 + IT_0042*IT_0044 + IT_0053*IT_0054) + (IT_0018*IT_0042 + IT_0016
      *IT_0044 + IT_0025*IT_0053 + IT_0024*IT_0054)*conj(IT_0057) + conj(IT_0042
      )*(IT_0015*IT_0044 + IT_0018*IT_0057 + IT_0042*IT_0064 + IT_0024*IT_0065 +
       IT_0053*IT_0067 + IT_0016*IT_0068 + IT_0063*IT_0070 + IT_0058*IT_0072) +
       conj(IT_0016)*(IT_0015*IT_0018 + IT_0044*IT_0057 + IT_0058*IT_0063 +
       IT_0016*IT_0064 + IT_0053*IT_0065 + IT_0024*IT_0067 + IT_0042*IT_0068 +
       IT_0070*IT_0072) + conj(IT_0070)*(IT_0042*IT_0063 + IT_0016*IT_0072 +
       IT_0070*IT_0074 + IT_0053*IT_0076 + IT_0024*IT_0078) + conj(IT_0058)*
      (IT_0016*IT_0063 + IT_0042*IT_0072 + IT_0058*IT_0074 + IT_0024*IT_0076 +
       IT_0053*IT_0078) + conj(IT_0053)*(IT_0015*IT_0054 + IT_0025*IT_0057 +
       IT_0016*IT_0065 + IT_0042*IT_0067 + IT_0070*IT_0076 + IT_0058*IT_0078 +
       s_13*IT_0053*IT_0079 + IT_0024*IT_0080) + conj(IT_0024)*(IT_0015*IT_0025 
      + IT_0054*IT_0057 + IT_0042*IT_0065 + IT_0016*IT_0067 + IT_0058*IT_0076 +
       IT_0070*IT_0078 + s_13*IT_0024*IT_0079 + IT_0053*IT_0080);
    return create_ccomplex_return(IT_0081);
}

