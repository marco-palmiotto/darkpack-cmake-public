#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_st_2_to_s_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_st_2_to_s_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = pow(m_sc_R, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0004;
    const ccomplex_t IT_0013 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0019 = IT_0004*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0017 + IT_0022;
    const ccomplex_t IT_0024 = IT_0011*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = m_t*V_tb*V_u2*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0031 = IT_0012*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0036 = IT_0004*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0034 + IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = -IT_0024 + -IT_0041;
    const ccomplex_t IT_0043 = m_s*IT_0034;
    const ccomplex_t IT_0044 = m_s*IT_0038;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = IT_0029*IT_0046;
    const ccomplex_t IT_0048 = m_s*IT_0021;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = m_s*IT_0017;
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = IT_0011*IT_0051;
    const ccomplex_t IT_0053 = IT_0047 + IT_0052;
    const ccomplex_t IT_0054 = m_s*s_14;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0056 = cos(beta);
    const ccomplex_t IT_0057 = cpow(IT_0056, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0000*IT_0004*IT_0057*conj(U_st_01);
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0026*IT_0059;
    const ccomplex_t IT_0061 = IT_0027*IT_0055*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0004*IT_0057*conj(U_st_01);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0006*IT_0064;
    const ccomplex_t IT_0066 = IT_0009*IT_0062*IT_0065;
    const ccomplex_t IT_0067 = IT_0061 + IT_0066;
    const ccomplex_t IT_0068 = m_b*s_13;
    const ccomplex_t IT_0069 = s_13*s_14;
    const ccomplex_t IT_0070 = s_34*IT_0008;
    const ccomplex_t IT_0071 = m_b*m_s;
    const ccomplex_t IT_0072 = 18*s_34;
    const ccomplex_t IT_0073 = conj(IT_0042)*(18*IT_0053*IT_0054 + (-18)
      *IT_0067*IT_0068 + IT_0042*(36*IT_0069 + (-18)*IT_0070)) + conj(IT_0053)*
      (18*IT_0042*IT_0054 + (-18)*IT_0067*IT_0071 + IT_0053*IT_0072) + conj
      (IT_0067)*((-18)*IT_0042*IT_0068 + (-18)*IT_0053*IT_0071 + IT_0067*IT_0072);
    return create_ccomplex_return(IT_0073);
}

