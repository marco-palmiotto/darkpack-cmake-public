#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_sc_R_to_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_sc_R_to_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_R = param->m_sc_R;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = 0.666666666666667*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(m_sc_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0020 = 0.666666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0028 = 0.666666666666667*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0036 = 0.666666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0021*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0011*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0044 + IT_0047 + IT_0050 + IT_0053;
    const ccomplex_t IT_0055 = IT_0012*IT_0013;
    const ccomplex_t IT_0056 = (-18)*(conj(IT_0041)*(IT_0041 + 1./3*IT_0054) +
       0.333333333333333*(IT_0041 + 3*IT_0054)*conj(IT_0054) + 0.333333333333333
      *IT_0054*(conj(IT_0041) + 3*conj(IT_0054)) + IT_0041*(conj(IT_0041) + 1./3
      *conj(IT_0054)))*IT_0055;
    const ccomplex_t IT_0057 = s_34*IT_0013;
    const ccomplex_t IT_0058 = (-6)*IT_0057;
    const ccomplex_t IT_0059 = s_23*s_24;
    const ccomplex_t IT_0060 = 12*IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = (-18)*IT_0057;
    const ccomplex_t IT_0063 = 36*IT_0059;
    const ccomplex_t IT_0064 = IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = m_c*s_23;
    const ccomplex_t IT_0066 = cpow(IT_0006, 2);
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_N_2*IT_0066;
    const ccomplex_t IT_0068 = IT_0051*IT_0067;
    const ccomplex_t IT_0069 = cpow(IT_0034, 2);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_N_4*IT_0069;
    const ccomplex_t IT_0071 = IT_0048*IT_0070;
    const ccomplex_t IT_0072 = m_c*IT_0044;
    const ccomplex_t IT_0073 = m_c*IT_0047;
    const ccomplex_t IT_0074 = cpow(IT_0026, 2);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_N_3*IT_0074;
    const ccomplex_t IT_0076 = IT_0045*IT_0075;
    const ccomplex_t IT_0077 = m_c*IT_0053;
    const ccomplex_t IT_0078 = m_c*IT_0050;
    const ccomplex_t IT_0079 = IT_0068 + IT_0071 + -IT_0072 + -IT_0073 +
       IT_0076 + -IT_0077 + -IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0018, 2);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_1*IT_0080;
    const ccomplex_t IT_0082 = IT_0042*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = m_c*IT_0024;
    const ccomplex_t IT_0085 = m_c*IT_0016;
    const ccomplex_t IT_0086 = m_c*IT_0032;
    const ccomplex_t IT_0087 = m_c*IT_0040;
    const ccomplex_t IT_0088 = IT_0022*IT_0081;
    const ccomplex_t IT_0089 = IT_0014*IT_0067;
    const ccomplex_t IT_0090 = -IT_0084 + -IT_0085 + -IT_0086 + -IT_0087 +
       IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = IT_0038*IT_0070;
    const ccomplex_t IT_0092 = IT_0030*IT_0075;
    const ccomplex_t IT_0093 = -IT_0091 + -IT_0092;
    const ccomplex_t IT_0094 = IT_0079 + -IT_0083 + 3*IT_0090 + (-3)*IT_0093;
    const ccomplex_t IT_0095 = (-3)*conj(IT_0093);
    const ccomplex_t IT_0096 = 3*conj(IT_0090);
    const ccomplex_t IT_0097 = -conj(IT_0083);
    const ccomplex_t IT_0098 = conj(IT_0079) + IT_0095 + IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = pow(m_sG, 2);
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -IT_0099 +
       reg_prop, -1);
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0102 = cpow(IT_0101, 2);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_sG*IT_0102;
    const ccomplex_t IT_0104 = IT_0100*IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0010, 2);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_N_2*IT_0105;
    const ccomplex_t IT_0107 = IT_0014*IT_0106;
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0104 + -IT_0107;
    const ccomplex_t IT_0109 = 6*IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0028, 2);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_N_3*IT_0110;
    const ccomplex_t IT_0112 = IT_0045*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0104 + -IT_0112;
    const ccomplex_t IT_0114 = 18*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -IT_0099 +
       reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0103*IT_0115;
    const ccomplex_t IT_0117 = cpow(IT_0036, 2);
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_N_4*IT_0117;
    const ccomplex_t IT_0119 = IT_0048*IT_0118;
    const ccomplex_t IT_0120 = cpow(IT_0020, 2);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_N_1*IT_0120;
    const ccomplex_t IT_0122 = IT_0042*IT_0121;
    const ccomplex_t IT_0123 = IT_0051*IT_0106;
    const ccomplex_t IT_0124 = -IT_0072 + -IT_0073 + -IT_0077 + -IT_0078 + (
      -0.166666666666667)*IT_0116 + IT_0119 + IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = (-18)*IT_0124;
    const ccomplex_t IT_0126 = IT_0030*IT_0111;
    const ccomplex_t IT_0127 = IT_0038*IT_0118;
    const ccomplex_t IT_0128 = IT_0022*IT_0121;
    const ccomplex_t IT_0129 = -IT_0084 + -IT_0085 + -IT_0086 + -IT_0087 + 0.5
      *IT_0116 + IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = (-6)*IT_0129;
    const ccomplex_t IT_0131 = IT_0109 + IT_0114 + IT_0125 + IT_0130;
    const ccomplex_t IT_0132 = 6*conj(IT_0108);
    const ccomplex_t IT_0133 = 18*conj(IT_0113);
    const ccomplex_t IT_0134 = (-18)*conj(IT_0124);
    const ccomplex_t IT_0135 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0136 = IT_0132 + IT_0133 + IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = 6*IT_0113;
    const ccomplex_t IT_0138 = 18*IT_0108;
    const ccomplex_t IT_0139 = (-18)*IT_0129;
    const ccomplex_t IT_0140 = (-6)*IT_0124;
    const ccomplex_t IT_0141 = IT_0137 + IT_0138 + IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = 6*conj(IT_0113);
    const ccomplex_t IT_0143 = 18*conj(IT_0108);
    const ccomplex_t IT_0144 = (-18)*conj(IT_0129);
    const ccomplex_t IT_0145 = (-6)*conj(IT_0124);
    const ccomplex_t IT_0146 = IT_0142 + IT_0143 + IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = conj(IT_0054)*((-18)*IT_0079 + 18*IT_0083 + (-6
      )*IT_0090 + 6*IT_0093) + IT_0054*((-18)*conj(IT_0079) + 18*conj(IT_0083) +
       (-6)*conj(IT_0090) + 6*conj(IT_0093)) + (-6)*conj(IT_0041)*IT_0094 + (-6)
      *IT_0041*IT_0098 + conj(IT_0054)*IT_0131 + IT_0054*IT_0136 + conj(IT_0041)
      *IT_0141 + IT_0041*IT_0146;
    const ccomplex_t IT_0148 = 6*IT_0124;
    const ccomplex_t IT_0149 = 18*IT_0129;
    const ccomplex_t IT_0150 = (-18)*IT_0108;
    const ccomplex_t IT_0151 = (-6)*IT_0113;
    const ccomplex_t IT_0152 = IT_0148 + IT_0149 + IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = 6*conj(IT_0124);
    const ccomplex_t IT_0154 = 18*conj(IT_0129);
    const ccomplex_t IT_0155 = (-18)*conj(IT_0108);
    const ccomplex_t IT_0156 = (-6)*conj(IT_0113);
    const ccomplex_t IT_0157 = IT_0153 + IT_0154 + IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = 6*IT_0129;
    const ccomplex_t IT_0159 = 18*IT_0124;
    const ccomplex_t IT_0160 = (-18)*IT_0113;
    const ccomplex_t IT_0161 = (-6)*IT_0108;
    const ccomplex_t IT_0162 = IT_0158 + IT_0159 + IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = 6*conj(IT_0129);
    const ccomplex_t IT_0164 = 18*conj(IT_0124);
    const ccomplex_t IT_0165 = (-18)*conj(IT_0113);
    const ccomplex_t IT_0166 = (-6)*conj(IT_0108);
    const ccomplex_t IT_0167 = IT_0163 + IT_0164 + IT_0165 + IT_0166;
    const ccomplex_t IT_0168 = m_c*s_24;
    const ccomplex_t IT_0169 = conj(IT_0054)*(18*IT_0079 + (-18)*IT_0083 + 6
      *IT_0090 + (-6)*IT_0093) + IT_0054*(18*conj(IT_0079) + (-18)*conj(IT_0083)
       + 6*conj(IT_0090) + (-6)*conj(IT_0093)) + 6*conj(IT_0041)*IT_0094 + 6
      *IT_0041*IT_0098 + conj(IT_0041)*IT_0152 + IT_0041*IT_0157 + conj(IT_0054)
      *IT_0162 + IT_0054*IT_0167;
    const ccomplex_t IT_0170 = 6*s_34;
    const ccomplex_t IT_0171 = 18*s_34;
    const ccomplex_t IT_0172 = (-18)*s_34;
    const ccomplex_t IT_0173 = (-6)*s_34;
    const ccomplex_t IT_0174 = IT_0056 + 2*(conj(IT_0041)*IT_0054 + IT_0041
      *conj(IT_0054))*IT_0061 + 2*(IT_0041*conj(IT_0041) + IT_0054*conj(IT_0054)
      )*IT_0064 + IT_0065*IT_0147 + IT_0012*(conj(IT_0079)*IT_0131 + IT_0079
      *IT_0136 + conj(IT_0090)*IT_0141 + IT_0090*IT_0146 + conj(IT_0093)*IT_0152
       + IT_0093*IT_0157 + conj(IT_0083)*IT_0162 + IT_0083*IT_0167) + IT_0168
      *IT_0169 + (conj(IT_0079)*IT_0090 + IT_0079*conj(IT_0090) + conj(IT_0083)
      *IT_0093 + IT_0083*conj(IT_0093) + conj(IT_0108)*IT_0113 + IT_0108*conj
      (IT_0113) + conj(IT_0124)*IT_0129 + IT_0124*conj(IT_0129))*IT_0170 + 
      (IT_0079*conj(IT_0079) + IT_0083*conj(IT_0083) + IT_0090*conj(IT_0090) +
       IT_0093*conj(IT_0093) + IT_0108*conj(IT_0108) + IT_0113*conj(IT_0113) +
       IT_0124*conj(IT_0124) + IT_0129*conj(IT_0129))*IT_0171 + (conj(IT_0079)
      *IT_0083 + IT_0079*conj(IT_0083) + conj(IT_0090)*IT_0093 + IT_0090*conj
      (IT_0093) + conj(IT_0113)*IT_0124 + IT_0113*conj(IT_0124) + conj(IT_0108)
      *IT_0129 + IT_0108*conj(IT_0129))*IT_0172 + (conj(IT_0083)*IT_0090 +
       IT_0083*conj(IT_0090) + conj(IT_0079)*IT_0093 + IT_0079*conj(IT_0093) +
       conj(IT_0108)*IT_0124 + IT_0108*conj(IT_0124) + conj(IT_0113)*IT_0129 +
       IT_0113*conj(IT_0129))*IT_0173;
    return create_ccomplex_return(IT_0174);
}

