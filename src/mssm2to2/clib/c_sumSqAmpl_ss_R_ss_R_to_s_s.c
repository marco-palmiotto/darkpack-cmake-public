#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_ss_R_to_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_ss_R_to_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = pow(m_ss_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0020 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0029 = (-0.333333333333333)*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0016 + IT_0024 + IT_0032 + IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0021*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0030*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0011*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0044 + -IT_0047 + -IT_0050 + -IT_0053;
    const ccomplex_t IT_0055 = s_34*IT_0013;
    const ccomplex_t IT_0056 = (-6)*IT_0055;
    const ccomplex_t IT_0057 = s_23*s_24;
    const ccomplex_t IT_0058 = 12*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = (-18)*IT_0055;
    const ccomplex_t IT_0061 = 36*IT_0057;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = m_s*s_24;
    const ccomplex_t IT_0064 = pow(m_sG, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -IT_0064 +
       reg_prop, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0067 = cpow(IT_0066, 2);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_sG*IT_0067;
    const ccomplex_t IT_0069 = IT_0065*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -IT_0064 +
       reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = cpow(IT_0029, 2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_N_3*IT_0072;
    const ccomplex_t IT_0074 = IT_0025*IT_0073;
    const ccomplex_t IT_0075 = cpow(IT_0010, 2);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_N_4*IT_0075;
    const ccomplex_t IT_0077 = IT_0014*IT_0076;
    const ccomplex_t IT_0078 = m_s*IT_0040;
    const ccomplex_t IT_0079 = m_s*IT_0032;
    const ccomplex_t IT_0080 = cpow(IT_0020, 2);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_1*IT_0080;
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0036, 2);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_2*IT_0083;
    const ccomplex_t IT_0085 = IT_0038*IT_0084;
    const ccomplex_t IT_0086 = m_s*IT_0024;
    const ccomplex_t IT_0087 = m_s*IT_0016;
    const ccomplex_t IT_0088 = 0.5*IT_0069 + (-0.166666666666667)*IT_0071 +
       IT_0074 + IT_0077 + -IT_0078 + -IT_0079 + IT_0082 + IT_0085 + -IT_0086 + 
      -IT_0087;
    const ccomplex_t IT_0089 = cpow(IT_0034, 2);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_2*IT_0089;
    const ccomplex_t IT_0091 = IT_0038*IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0027, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_N_3*IT_0092;
    const ccomplex_t IT_0094 = IT_0025*IT_0093;
    const ccomplex_t IT_0095 = cpow(IT_0018, 2);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_N_1*IT_0095;
    const ccomplex_t IT_0097 = IT_0022*IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0006, 2);
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_4*IT_0098;
    const ccomplex_t IT_0100 = IT_0014*IT_0099;
    const ccomplex_t IT_0101 = -IT_0078 + -IT_0079 + -IT_0086 + -IT_0087 +
       IT_0091 + IT_0094 + IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0048*IT_0084;
    const ccomplex_t IT_0103 = m_s*IT_0050;
    const ccomplex_t IT_0104 = IT_0045*IT_0073;
    const ccomplex_t IT_0105 = IT_0051*IT_0076;
    const ccomplex_t IT_0106 = m_s*IT_0044;
    const ccomplex_t IT_0107 = IT_0042*IT_0081;
    const ccomplex_t IT_0108 = m_s*IT_0047;
    const ccomplex_t IT_0109 = m_s*IT_0053;
    const ccomplex_t IT_0110 = (-0.166666666666667)*IT_0069 + 0.5*IT_0071 +
       IT_0102 + -IT_0103 + IT_0104 + IT_0105 + -IT_0106 + IT_0107 + -IT_0108 + 
      -IT_0109;
    const ccomplex_t IT_0111 = IT_0041 + 3*IT_0054;
    const ccomplex_t IT_0112 = 3*conj(IT_0054);
    const ccomplex_t IT_0113 = conj(IT_0041) + IT_0112;
    const ccomplex_t IT_0114 = IT_0045*IT_0093;
    const ccomplex_t IT_0115 = IT_0042*IT_0096;
    const ccomplex_t IT_0116 = IT_0048*IT_0090;
    const ccomplex_t IT_0117 = IT_0051*IT_0099;
    const ccomplex_t IT_0118 = -IT_0103 + -IT_0106 + -IT_0108 + -IT_0109 +
       IT_0114 + IT_0115 + IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = 3*IT_0054;
    const ccomplex_t IT_0120 = IT_0041 + IT_0119;
    const ccomplex_t IT_0121 = (18*conj(IT_0041) + 6*conj(IT_0054))*IT_0088 + 
      (18*IT_0041 + 6*IT_0054)*conj(IT_0088) + (18*conj(IT_0041) + 6*conj
      (IT_0054))*IT_0101 + (18*IT_0041 + 6*IT_0054)*conj(IT_0101) + 6*conj
      (IT_0110)*IT_0111 + 6*IT_0110*IT_0113 + 6*IT_0113*IT_0118 + 6*conj(IT_0118
      )*IT_0120;
    const ccomplex_t IT_0122 = 6*s_34;
    const ccomplex_t IT_0123 = 18*s_34;
    const ccomplex_t IT_0124 = (-18)*conj(IT_0088);
    const ccomplex_t IT_0125 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0126 = IT_0124 + IT_0125;
    const ccomplex_t IT_0127 = (-18)*conj(IT_0110);
    const ccomplex_t IT_0128 = (-6)*conj(IT_0088);
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = (-18)*IT_0110;
    const ccomplex_t IT_0131 = (-6)*IT_0110;
    const ccomplex_t IT_0132 = IT_0012*IT_0013;
    const ccomplex_t IT_0133 = (-6)*conj(IT_0054);
    const ccomplex_t IT_0134 = (-18)*conj(IT_0041);
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = (-6)*IT_0054;
    const ccomplex_t IT_0137 = m_s*s_23;
    const ccomplex_t IT_0138 = IT_0113*IT_0118 + IT_0111*conj(IT_0118) + (
      -0.166666666666667)*IT_0041*IT_0126 + (-0.166666666666667)*IT_0054*IT_0129
       + (-0.166666666666667)*conj(IT_0054)*IT_0130 + (-0.166666666666667)*conj
      (IT_0041)*IT_0131 + (-0.166666666666667)*IT_0101*((-18)*conj(IT_0041) +
       IT_0133) + (-0.166666666666667)*IT_0088*IT_0135 + (-0.166666666666667)
      *conj(IT_0101)*((-18)*IT_0041 + IT_0136);
    const ccomplex_t IT_0139 = 2*(conj(IT_0041)*IT_0054 + IT_0041*conj(IT_0054
      ))*IT_0059 + 2*(IT_0041*conj(IT_0041) + IT_0054*conj(IT_0054))*IT_0062 +
       IT_0063*IT_0121 + (conj(IT_0088)*IT_0110 + IT_0088*conj(IT_0110) + conj
      (IT_0101)*IT_0118 + IT_0101*conj(IT_0118))*IT_0122 + (IT_0088*conj(IT_0088
      ) + IT_0101*conj(IT_0101) + IT_0110*conj(IT_0110) + IT_0118*conj(IT_0118))
      *IT_0123 + IT_0012*(IT_0088*((-18)*conj(IT_0101) + (-6)*conj(IT_0118)) +
       IT_0101*IT_0126 + IT_0118*IT_0129 + conj(IT_0118)*IT_0130 + conj(IT_0101)
      *IT_0131) + IT_0132*(conj(IT_0041)*((-18)*IT_0041 + (-6)*IT_0054) +
       IT_0120*IT_0133 + IT_0041*IT_0135 + IT_0113*IT_0136) + (-6)*IT_0137
      *IT_0138;
    return create_ccomplex_return(IT_0139);
}

