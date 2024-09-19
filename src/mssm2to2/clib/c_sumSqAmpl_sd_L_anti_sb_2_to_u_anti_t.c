#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sb_2_to_u_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sb_2_to_u_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_u, 2);
    const ccomplex_t IT_0002 = pow(m_sd_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0005*IT_0007*IT_0009*conj(U_sb_01);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0003*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0005;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0005*IT_0007*IT_0009*conj(U_sb_01);
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = IT_0014*IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = IT_0015*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = m_u*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0011*IT_0028;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_u*IT_0031;
    const ccomplex_t IT_0033 = -IT_0013 + -IT_0020 + IT_0026 + IT_0032;
    const ccomplex_t IT_0034 = -IT_0025 + -IT_0031;
    const ccomplex_t IT_0035 = m_u*s_14;
    const ccomplex_t IT_0036 = 18*s_34;
    const ccomplex_t IT_0037 = m_t*m_u;
    const ccomplex_t IT_0038 = V_tb*e_em*U_Wm1*conj(U_sb_01);
    const ccomplex_t IT_0039 = IT_0005*IT_0038;
    const ccomplex_t IT_0040 = cos(beta);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = IT_0005*IT_0041;
    const ccomplex_t IT_0043 = m_b*U_d1*V_tb*e_em*IT_0007*conj(U_sb_11);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0039 + (-0.5)*IT_0045);
    const ccomplex_t IT_0047 = IT_0028*IT_0046;
    const ccomplex_t IT_0048 = IT_0000*IT_0003*IT_0047;
    const ccomplex_t IT_0049 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0050 = IT_0005*IT_0049;
    const ccomplex_t IT_0051 = m_b*U_d2*V_tb*e_em*IT_0007*conj(U_sb_11);
    const ccomplex_t IT_0052 = IT_0042*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0053);
    const ccomplex_t IT_0055 = IT_0022*IT_0054;
    const ccomplex_t IT_0056 = IT_0014*IT_0015*IT_0055;
    const ccomplex_t IT_0057 = IT_0006*IT_0046;
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = m_u*IT_0059;
    const ccomplex_t IT_0061 = IT_0016*IT_0054;
    const ccomplex_t IT_0062 = IT_0015*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = m_u*IT_0063;
    const ccomplex_t IT_0065 = -IT_0048 + -IT_0056 + IT_0060 + IT_0064;
    const ccomplex_t IT_0066 = (-18)*conj(IT_0065);
    const ccomplex_t IT_0067 = m_t*s_13;
    const ccomplex_t IT_0068 = -IT_0059 + -IT_0063;
    const ccomplex_t IT_0069 = (-18)*conj(IT_0068);
    const ccomplex_t IT_0070 = s_34*IT_0002;
    const ccomplex_t IT_0071 = (-18)*IT_0070;
    const ccomplex_t IT_0072 = s_13*s_14;
    const ccomplex_t IT_0073 = 36*IT_0072;
    const ccomplex_t IT_0074 = IT_0071 + IT_0073;
    const ccomplex_t IT_0075 = m_t*m_u*IT_0002;
    const ccomplex_t IT_0076 = (-18)*conj(IT_0034);
    const ccomplex_t IT_0077 = IT_0033*(18*conj(IT_0034)*IT_0035 + conj
      (IT_0033)*IT_0036 + IT_0037*IT_0066 + IT_0067*IT_0069) + IT_0034*(18*conj
      (IT_0033)*IT_0035 + IT_0066*IT_0067 + conj(IT_0034)*IT_0074 + IT_0069
      *IT_0075) + IT_0065*((-18)*conj(IT_0033)*IT_0037 + IT_0036*conj(IT_0065) +
       18*IT_0035*conj(IT_0068) + IT_0067*IT_0076) + IT_0068*(18*IT_0035*conj
      (IT_0065) + (-18)*conj(IT_0033)*IT_0067 + conj(IT_0068)*IT_0074 + IT_0075
      *IT_0076);
    return create_ccomplex_return(IT_0077);
}

