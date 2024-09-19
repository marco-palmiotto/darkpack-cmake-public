#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_sb_2_to_c_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_sb_2_to_c_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_ss_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0006 = V_tb*e_em*U_Wm1*conj(U_sb_01);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = m_b*U_d1*V_tb*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0007 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0004;
    const ccomplex_t IT_0021 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = m_b*U_d2*V_tb*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0024 = IT_0010*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = IT_0019*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0018 + -IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0032 = sin(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0004*IT_0011*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = IT_0019*IT_0031*IT_0036;
    const ccomplex_t IT_0038 = m_c*IT_0018;
    const ccomplex_t IT_0039 = m_c*IT_0029;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0004*IT_0011*IT_0033;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0015*IT_0042;
    const ccomplex_t IT_0044 = IT_0002*IT_0040*IT_0043;
    const ccomplex_t IT_0045 = -IT_0037 + IT_0038 + IT_0039 + -IT_0044;
    const ccomplex_t IT_0046 = m_c*s_14;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0004*IT_0011*IT_0033*conj(U_sb_01);
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0005*IT_0048;
    const ccomplex_t IT_0050 = IT_0002*IT_0040*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0004*IT_0011*IT_0033*conj(U_sb_01);
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0035*IT_0052;
    const ccomplex_t IT_0054 = IT_0019*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_c*IT_0055;
    const ccomplex_t IT_0057 = IT_0020*IT_0052;
    const ccomplex_t IT_0058 = IT_0019*IT_0031*IT_0057;
    const ccomplex_t IT_0059 = IT_0042*IT_0048;
    const ccomplex_t IT_0060 = IT_0002*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = m_c*IT_0061;
    const ccomplex_t IT_0063 = -IT_0050 + IT_0056 + -IT_0058 + IT_0062;
    const ccomplex_t IT_0064 = m_t*s_13;
    const ccomplex_t IT_0065 = -IT_0055 + -IT_0061;
    const ccomplex_t IT_0066 = m_c*m_t*IT_0001;
    const ccomplex_t IT_0067 = s_34*IT_0001;
    const ccomplex_t IT_0068 = (-18)*IT_0067;
    const ccomplex_t IT_0069 = s_13*s_14;
    const ccomplex_t IT_0070 = 36*IT_0069;
    const ccomplex_t IT_0071 = IT_0068 + IT_0070;
    const ccomplex_t IT_0072 = 18*s_34;
    const ccomplex_t IT_0073 = m_c*m_t;
    const ccomplex_t IT_0074 = (-18)*conj(IT_0045);
    const ccomplex_t IT_0075 = (-18)*IT_0045;
    const ccomplex_t IT_0076 = (-18)*IT_0030;
    const ccomplex_t IT_0077 = conj(IT_0030)*(18*IT_0045*IT_0046 + (-18)
      *IT_0063*IT_0064 + (-18)*IT_0065*IT_0066 + IT_0030*IT_0071) + conj(IT_0045
      )*(18*IT_0030*IT_0046 + IT_0045*IT_0072) + (IT_0064*IT_0065 + IT_0063
      *IT_0073)*IT_0074 + conj(IT_0063)*(18*IT_0046*IT_0065 + IT_0063*IT_0072 +
       IT_0073*IT_0075 + IT_0064*IT_0076) + conj(IT_0065)*(18*IT_0046*IT_0063 +
       IT_0065*IT_0071 + IT_0064*IT_0075 + IT_0066*IT_0076);
    return create_ccomplex_return(IT_0077);
}

