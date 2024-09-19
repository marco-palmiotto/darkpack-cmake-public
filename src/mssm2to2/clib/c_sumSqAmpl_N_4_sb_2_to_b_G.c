#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sb_2_to_b_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sb_2_to_b_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_N_4, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = m_b*N_d4*e_em*IT_0015*U_sb_01;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0009 + 1.5*IT_0018);
    const ccomplex_t IT_0020 = 0.666666666666667*IT_0019;
    const ccomplex_t IT_0021 = 1.5*IT_0020;
    const ccomplex_t IT_0022 = (-0.333333333333333)*IT_0021;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (-2)*IT_0004;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (-2)*IT_0025 + -IT_0030;
    const ccomplex_t IT_0032 = m_b*conj(N_d4)*e_em*IT_0015*U_sb_11;
    const ccomplex_t IT_0033 = IT_0014*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0039 = IT_0013*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0034 + 0.333333333333333
      *IT_0037 + -IT_0040);
    const ccomplex_t IT_0042 = 3*IT_0041;
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = IT_0027*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0026*IT_0045;
    const ccomplex_t IT_0047 = IT_0004*IT_0043;
    const ccomplex_t IT_0048 = IT_0003*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0046 + 2*IT_0049;
    const ccomplex_t IT_0051 = (-8)*IT_0000;
    const ccomplex_t IT_0052 = -IT_0004;
    const ccomplex_t IT_0053 = IT_0029*IT_0052;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0045*IT_0052;
    const ccomplex_t IT_0056 = (-8)*s_13;
    const ccomplex_t IT_0057 = s_34*IT_0056;
    const ccomplex_t IT_0058 = m_b*s_14*m_N_4;
    const ccomplex_t IT_0059 = -IT_0025;
    const ccomplex_t IT_0060 = s_34*IT_0001;
    const ccomplex_t IT_0061 = 8*IT_0060;
    const ccomplex_t IT_0062 = s_13*s_14;
    const ccomplex_t IT_0063 = (-16)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = s_14*IT_0056;
    const ccomplex_t IT_0066 = IT_0001*IT_0056;
    const ccomplex_t IT_0067 = pow(s_13, 2);
    const ccomplex_t IT_0068 = (-8)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0046;
    const ccomplex_t IT_0070 = m_b*m_N_4;
    const ccomplex_t IT_0071 = (-8)*s_14;
    const ccomplex_t IT_0072 = (-8)*s_34;
    const ccomplex_t IT_0073 = (-8)*IT_0001;
    const ccomplex_t IT_0074 = (conj(IT_0031)*IT_0050 + IT_0031*conj(IT_0050))
      *IT_0051 + IT_0056*(conj(IT_0030)*IT_0050 + IT_0030*conj(IT_0050) + conj
      (IT_0031)*IT_0069 + IT_0031*conj(IT_0069)) + (conj(IT_0030)*IT_0055 +
       IT_0030*conj(IT_0055) + conj(IT_0054)*IT_0069 + IT_0054*conj(IT_0069))
      *IT_0071 + (conj(IT_0050)*IT_0054 + IT_0050*conj(IT_0054) + conj(IT_0031)
      *IT_0055 + IT_0031*conj(IT_0055))*IT_0072 + (conj(IT_0030)*IT_0069 +
       IT_0030*conj(IT_0069))*IT_0073;
    const ccomplex_t IT_0075 = m_b*s_34*m_N_4;
    const ccomplex_t IT_0076 = (-8)*conj(IT_0031);
    const ccomplex_t IT_0077 = s_14*IT_0000;
    const ccomplex_t IT_0078 = (-8)*conj(IT_0050);
    const ccomplex_t IT_0079 = s_14*s_34;
    const ccomplex_t IT_0080 = (-8)*IT_0050;
    const ccomplex_t IT_0081 = (-0.125)*IT_0077;
    const ccomplex_t IT_0082 = (-8)*IT_0031;
    const ccomplex_t IT_0083 = (-0.125)*IT_0075;
    const ccomplex_t IT_0084 = s_13*(IT_0031*conj(IT_0031) + IT_0050*conj
      (IT_0050))*IT_0051 + (conj(IT_0031)*IT_0054 + IT_0031*conj(IT_0054) + conj
      (IT_0050)*IT_0055 + IT_0050*conj(IT_0055))*IT_0057 + (-8)*conj(IT_0030)*
      (IT_0049*IT_0058 + (-0.125)*IT_0059*IT_0064 + (-0.125)*IT_0054*IT_0065 + (
      -0.125)*IT_0030*IT_0066 + (-0.125)*IT_0031*IT_0068) + IT_0030*(conj
      (IT_0054)*IT_0065 + conj(IT_0031)*IT_0068) + (conj(IT_0055)*IT_0065 + conj
      (IT_0050)*IT_0068)*IT_0069 + (-8)*(IT_0058*IT_0059 + (-0.125)*IT_0049
      *IT_0064 + (-0.125)*IT_0055*IT_0065 + (-0.125)*IT_0050*IT_0068 + (-0.125)
      *IT_0066*IT_0069)*conj(IT_0069) + IT_0070*IT_0074 + IT_0059*(IT_0076
      *IT_0077 + IT_0075*IT_0078) + IT_0049*(IT_0075*IT_0076 + IT_0077*IT_0078) 
      + (-8)*conj(IT_0059)*((-0.125)*IT_0030*IT_0064 + IT_0058*IT_0069 + (-4)
      *IT_0059*IT_0079 + IT_0081*IT_0082 + IT_0080*IT_0083) + (-8)*conj(IT_0049)
      *(IT_0030*IT_0058 + (-0.125)*IT_0064*IT_0069 + (-4)*IT_0049*IT_0079 +
       IT_0080*IT_0081 + IT_0082*IT_0083);
    return create_ccomplex_return(IT_0084);
}

