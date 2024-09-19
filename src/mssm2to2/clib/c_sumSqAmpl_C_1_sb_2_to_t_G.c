#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_2_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_2_to_t_G(
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0003*IT_0005*U_sb_01;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = pow(m_C_1, 2);
    const ccomplex_t IT_0010 = pow(m_sb_2, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-2)*IT_0000;
    const ccomplex_t IT_0015 = pow(m_t, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0009 + IT_0015 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = (-2)*IT_0013 + -IT_0019;
    const ccomplex_t IT_0021 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0022 = IT_0005*IT_0021;
    const ccomplex_t IT_0023 = cos(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = m_b*conj(U_d1)*V_tb*e_em*IT_0001*U_sb_11;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = IT_0000*IT_0029;
    const ccomplex_t IT_0031 = IT_0011*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0016*IT_0029;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0014*IT_0034;
    const ccomplex_t IT_0036 = 2*IT_0032 + IT_0035;
    const ccomplex_t IT_0037 = (-8)*IT_0015;
    const ccomplex_t IT_0038 = -IT_0000;
    const ccomplex_t IT_0039 = IT_0018*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0034*IT_0038;
    const ccomplex_t IT_0042 = (-8)*s_13;
    const ccomplex_t IT_0043 = s_34*IT_0042;
    const ccomplex_t IT_0044 = m_t*s_14*m_C_1;
    const ccomplex_t IT_0045 = -IT_0013;
    const ccomplex_t IT_0046 = s_34*IT_0009;
    const ccomplex_t IT_0047 = 8*IT_0046;
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = (-16)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = s_14*IT_0042;
    const ccomplex_t IT_0052 = IT_0009*IT_0042;
    const ccomplex_t IT_0053 = pow(s_13, 2);
    const ccomplex_t IT_0054 = (-8)*IT_0053;
    const ccomplex_t IT_0055 = -IT_0035;
    const ccomplex_t IT_0056 = m_t*m_C_1;
    const ccomplex_t IT_0057 = (-8)*s_14;
    const ccomplex_t IT_0058 = (-8)*s_34;
    const ccomplex_t IT_0059 = (-8)*IT_0009;
    const ccomplex_t IT_0060 = (conj(IT_0020)*IT_0036 + IT_0020*conj(IT_0036))
      *IT_0037 + IT_0042*(conj(IT_0019)*IT_0036 + IT_0019*conj(IT_0036) + conj
      (IT_0020)*IT_0055 + IT_0020*conj(IT_0055)) + (conj(IT_0019)*IT_0041 +
       IT_0019*conj(IT_0041) + conj(IT_0040)*IT_0055 + IT_0040*conj(IT_0055))
      *IT_0057 + (conj(IT_0036)*IT_0040 + IT_0036*conj(IT_0040) + conj(IT_0020)
      *IT_0041 + IT_0020*conj(IT_0041))*IT_0058 + (conj(IT_0019)*IT_0055 +
       IT_0019*conj(IT_0055))*IT_0059;
    const ccomplex_t IT_0061 = m_t*s_34*m_C_1;
    const ccomplex_t IT_0062 = (-8)*conj(IT_0020);
    const ccomplex_t IT_0063 = s_14*IT_0015;
    const ccomplex_t IT_0064 = (-8)*conj(IT_0036);
    const ccomplex_t IT_0065 = s_14*s_34;
    const ccomplex_t IT_0066 = (-8)*IT_0036;
    const ccomplex_t IT_0067 = (-0.125)*IT_0063;
    const ccomplex_t IT_0068 = (-8)*IT_0020;
    const ccomplex_t IT_0069 = (-0.125)*IT_0061;
    const ccomplex_t IT_0070 = s_13*(IT_0020*conj(IT_0020) + IT_0036*conj
      (IT_0036))*IT_0037 + (conj(IT_0020)*IT_0040 + IT_0020*conj(IT_0040) + conj
      (IT_0036)*IT_0041 + IT_0036*conj(IT_0041))*IT_0043 + (-8)*conj(IT_0019)*
      (IT_0032*IT_0044 + (-0.125)*IT_0045*IT_0050 + (-0.125)*IT_0040*IT_0051 + (
      -0.125)*IT_0019*IT_0052 + (-0.125)*IT_0020*IT_0054) + IT_0019*(conj
      (IT_0040)*IT_0051 + conj(IT_0020)*IT_0054) + (conj(IT_0041)*IT_0051 + conj
      (IT_0036)*IT_0054)*IT_0055 + (-8)*(IT_0044*IT_0045 + (-0.125)*IT_0032
      *IT_0050 + (-0.125)*IT_0041*IT_0051 + (-0.125)*IT_0036*IT_0054 + (-0.125)
      *IT_0052*IT_0055)*conj(IT_0055) + IT_0056*IT_0060 + IT_0045*(IT_0062
      *IT_0063 + IT_0061*IT_0064) + IT_0032*(IT_0061*IT_0062 + IT_0063*IT_0064) 
      + (-8)*conj(IT_0045)*((-0.125)*IT_0019*IT_0050 + IT_0044*IT_0055 + (-4)
      *IT_0045*IT_0065 + IT_0067*IT_0068 + IT_0066*IT_0069) + (-8)*conj(IT_0032)
      *(IT_0019*IT_0044 + (-0.125)*IT_0050*IT_0055 + (-4)*IT_0032*IT_0065 +
       IT_0066*IT_0067 + IT_0068*IT_0069);
    return create_ccomplex_return(IT_0070);
}

