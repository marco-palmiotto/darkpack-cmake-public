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
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = pow(m_ss_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0021 = (-0.333333333333333)*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_3*
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
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0037 = (-0.333333333333333)*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = s_34*IT_0013;
    const ccomplex_t IT_0043 = (-18)*IT_0042;
    const ccomplex_t IT_0044 = s_23*s_24;
    const ccomplex_t IT_0045 = 36*IT_0044;
    const ccomplex_t IT_0046 = IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = IT_0012*IT_0013;
    const ccomplex_t IT_0048 = (-36)*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0046 + IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0022*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0011*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0030*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0038*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0052 + IT_0055 + IT_0058 + IT_0061;
    const ccomplex_t IT_0063 = (-6)*IT_0042;
    const ccomplex_t IT_0064 = 2*IT_0063;
    const ccomplex_t IT_0065 = 12*IT_0044;
    const ccomplex_t IT_0066 = 2*IT_0065;
    const ccomplex_t IT_0067 = (-12)*IT_0047;
    const ccomplex_t IT_0068 = IT_0064 + IT_0066 + IT_0067;
    const ccomplex_t IT_0069 = pow(m_sG, 2);
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -IT_0069 +
       reg_prop, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0072 = cpow(IT_0071, 2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_sG*IT_0072;
    const ccomplex_t IT_0074 = IT_0070*IT_0073;
    const ccomplex_t IT_0075 = cpow(IT_0029, 2);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_N_3*IT_0075;
    const ccomplex_t IT_0077 = IT_0056*IT_0076;
    const ccomplex_t IT_0078 = m_s*IT_0052;
    const ccomplex_t IT_0079 = m_s*IT_0055;
    const ccomplex_t IT_0080 = m_s*IT_0058;
    const ccomplex_t IT_0081 = m_s*IT_0061;
    const ccomplex_t IT_0082 = 0.5*IT_0074 + IT_0077 + -IT_0078 + -IT_0079 + 
      -IT_0080 + -IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0010, 2);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_2*IT_0083;
    const ccomplex_t IT_0085 = IT_0014*IT_0084;
    const ccomplex_t IT_0086 = m_s*IT_0024;
    const ccomplex_t IT_0087 = m_s*IT_0016;
    const ccomplex_t IT_0088 = m_s*IT_0032;
    const ccomplex_t IT_0089 = m_s*IT_0040;
    const ccomplex_t IT_0090 = (-0.166666666666667)*IT_0074 + IT_0085 + 
      -IT_0086 + -IT_0087 + -IT_0088 + -IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -IT_0069 +
       reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0073*IT_0091;
    const ccomplex_t IT_0093 = IT_0053*IT_0084;
    const ccomplex_t IT_0094 = cpow(IT_0037, 2);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_N_4*IT_0094;
    const ccomplex_t IT_0096 = IT_0059*IT_0095;
    const ccomplex_t IT_0097 = cpow(IT_0021, 2);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_1*IT_0097;
    const ccomplex_t IT_0099 = IT_0050*IT_0098;
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0092 + -IT_0093 + -IT_0096
       + -IT_0099;
    const ccomplex_t IT_0101 = IT_0017*IT_0098;
    const ccomplex_t IT_0102 = IT_0033*IT_0095;
    const ccomplex_t IT_0103 = IT_0025*IT_0076;
    const ccomplex_t IT_0104 = (-0.5)*IT_0092 + -IT_0101 + -IT_0102 + -IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0027, 2);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_N_3*IT_0105;
    const ccomplex_t IT_0107 = IT_0056*IT_0106;
    const ccomplex_t IT_0108 = -IT_0078 + -IT_0079 + -IT_0080 + -IT_0081 +
       IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0019, 2);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_N_1*IT_0109;
    const ccomplex_t IT_0111 = IT_0017*IT_0110;
    const ccomplex_t IT_0112 = cpow(IT_0035, 2);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*m_N_4*IT_0112;
    const ccomplex_t IT_0114 = IT_0033*IT_0113;
    const ccomplex_t IT_0115 = -IT_0086 + -IT_0087 + -IT_0088 + -IT_0089 +
       IT_0111 + IT_0114;
    const ccomplex_t IT_0116 = IT_0059*IT_0113;
    const ccomplex_t IT_0117 = IT_0050*IT_0110;
    const ccomplex_t IT_0118 = cpow(IT_0006, 2);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_2*IT_0118;
    const ccomplex_t IT_0120 = IT_0053*IT_0119;
    const ccomplex_t IT_0121 = -IT_0116 + -IT_0117 + -IT_0120;
    const ccomplex_t IT_0122 = IT_0014*IT_0119;
    const ccomplex_t IT_0123 = IT_0025*IT_0106;
    const ccomplex_t IT_0124 = -IT_0122 + -IT_0123;
    const ccomplex_t IT_0125 = 6*s_34;
    const ccomplex_t IT_0126 = 18*s_34;
    const ccomplex_t IT_0127 = (-18)*s_34;
    const ccomplex_t IT_0128 = (-6)*s_34;
    const ccomplex_t IT_0129 = m_s*s_23;
    const ccomplex_t IT_0130 = IT_0108 + 3*IT_0115 + -IT_0121 + (-3)*IT_0124;
    const ccomplex_t IT_0131 = 6*IT_0100;
    const ccomplex_t IT_0132 = 18*IT_0104;
    const ccomplex_t IT_0133 = (-18)*IT_0090;
    const ccomplex_t IT_0134 = (-6)*IT_0082;
    const ccomplex_t IT_0135 = (-3)*conj(IT_0124);
    const ccomplex_t IT_0136 = 3*conj(IT_0115);
    const ccomplex_t IT_0137 = -conj(IT_0121);
    const ccomplex_t IT_0138 = conj(IT_0108) + IT_0135 + IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = 6*conj(IT_0100);
    const ccomplex_t IT_0140 = 18*conj(IT_0104);
    const ccomplex_t IT_0141 = (-18)*conj(IT_0090);
    const ccomplex_t IT_0142 = (-6)*conj(IT_0082);
    const ccomplex_t IT_0143 = 6*IT_0104;
    const ccomplex_t IT_0144 = 18*IT_0100;
    const ccomplex_t IT_0145 = (-18)*IT_0082;
    const ccomplex_t IT_0146 = (-6)*IT_0090;
    const ccomplex_t IT_0147 = 6*conj(IT_0104);
    const ccomplex_t IT_0148 = 18*conj(IT_0100);
    const ccomplex_t IT_0149 = (-18)*conj(IT_0082);
    const ccomplex_t IT_0150 = (-6)*conj(IT_0090);
    const ccomplex_t IT_0151 = conj(IT_0041)*(IT_0130 + (-0.166666666666667)
      *IT_0131 + (-0.166666666666667)*IT_0132 + (-0.166666666666667)*IT_0133 + (
      -0.166666666666667)*IT_0134) + IT_0041*(IT_0138 + (-0.166666666666667)
      *IT_0139 + (-0.166666666666667)*IT_0140 + (-0.166666666666667)*IT_0141 + (
      -0.166666666666667)*IT_0142) + (-0.166666666666667)*conj(IT_0062)*((-18)
      *IT_0108 + (-6)*IT_0115 + 18*IT_0121 + 6*IT_0124 + IT_0143 + IT_0144 +
       IT_0145 + IT_0146) + (-0.166666666666667)*IT_0062*((-18)*conj(IT_0108) + 
      (-6)*conj(IT_0115) + 18*conj(IT_0121) + 6*conj(IT_0124) + IT_0147 +
       IT_0148 + IT_0149 + IT_0150);
    const ccomplex_t IT_0152 = 6*IT_0082;
    const ccomplex_t IT_0153 = 18*IT_0090;
    const ccomplex_t IT_0154 = (-18)*IT_0104;
    const ccomplex_t IT_0155 = (-6)*IT_0100;
    const ccomplex_t IT_0156 = 6*conj(IT_0082);
    const ccomplex_t IT_0157 = 18*conj(IT_0090);
    const ccomplex_t IT_0158 = (-18)*conj(IT_0104);
    const ccomplex_t IT_0159 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0160 = 6*IT_0090;
    const ccomplex_t IT_0161 = 18*IT_0082;
    const ccomplex_t IT_0162 = (-18)*IT_0100;
    const ccomplex_t IT_0163 = (-6)*IT_0104;
    const ccomplex_t IT_0164 = 6*conj(IT_0090);
    const ccomplex_t IT_0165 = 18*conj(IT_0082);
    const ccomplex_t IT_0166 = (-18)*conj(IT_0100);
    const ccomplex_t IT_0167 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0168 = conj(IT_0115)*(IT_0131 + IT_0132 + IT_0133 +
       IT_0134) + IT_0115*(IT_0139 + IT_0140 + IT_0141 + IT_0142) + conj(IT_0108
      )*(IT_0143 + IT_0144 + IT_0145 + IT_0146) + IT_0108*(IT_0147 + IT_0148 +
       IT_0149 + IT_0150) + conj(IT_0124)*(IT_0152 + IT_0153 + IT_0154 + IT_0155
      ) + IT_0124*(IT_0156 + IT_0157 + IT_0158 + IT_0159) + conj(IT_0121)*
      (IT_0160 + IT_0161 + IT_0162 + IT_0163) + IT_0121*(IT_0164 + IT_0165 +
       IT_0166 + IT_0167);
    const ccomplex_t IT_0169 = m_s*s_24;
    const ccomplex_t IT_0170 = conj(IT_0041)*(IT_0130 + 0.166666666666667
      *IT_0152 + 0.166666666666667*IT_0153 + 0.166666666666667*IT_0154 +
       0.166666666666667*IT_0155) + IT_0041*(IT_0138 + 0.166666666666667*IT_0156
       + 0.166666666666667*IT_0157 + 0.166666666666667*IT_0158 +
       0.166666666666667*IT_0159) + 0.166666666666667*conj(IT_0062)*(18*IT_0108 
      + 6*IT_0115 + (-18)*IT_0121 + (-6)*IT_0124 + IT_0160 + IT_0161 + IT_0162 +
       IT_0163) + 0.166666666666667*IT_0062*(18*conj(IT_0108) + 6*conj(IT_0115) 
      + (-18)*conj(IT_0121) + (-6)*conj(IT_0124) + IT_0164 + IT_0165 + IT_0166 +
       IT_0167);
    const ccomplex_t IT_0171 = IT_0062*(IT_0049*conj(IT_0062) + conj(IT_0041)
      *IT_0068) + IT_0041*(conj(IT_0041)*IT_0049 + conj(IT_0062)*IT_0068) + 
      (conj(IT_0082)*IT_0090 + IT_0082*conj(IT_0090) + conj(IT_0100)*IT_0104 +
       IT_0100*conj(IT_0104) + conj(IT_0108)*IT_0115 + IT_0108*conj(IT_0115) +
       conj(IT_0121)*IT_0124 + IT_0121*conj(IT_0124))*IT_0125 + (IT_0082*conj
      (IT_0082) + IT_0090*conj(IT_0090) + IT_0100*conj(IT_0100) + IT_0104*conj
      (IT_0104) + IT_0108*conj(IT_0108) + IT_0115*conj(IT_0115) + IT_0121*conj
      (IT_0121) + IT_0124*conj(IT_0124))*IT_0126 + (conj(IT_0082)*IT_0100 +
       IT_0082*conj(IT_0100) + conj(IT_0090)*IT_0104 + IT_0090*conj(IT_0104) +
       conj(IT_0108)*IT_0121 + IT_0108*conj(IT_0121) + conj(IT_0115)*IT_0124 +
       IT_0115*conj(IT_0124))*IT_0127 + (conj(IT_0090)*IT_0100 + IT_0090*conj
      (IT_0100) + conj(IT_0082)*IT_0104 + IT_0082*conj(IT_0104) + conj(IT_0115)
      *IT_0121 + IT_0115*conj(IT_0121) + conj(IT_0108)*IT_0124 + IT_0108*conj
      (IT_0124))*IT_0128 + (-6)*IT_0129*IT_0151 + IT_0012*IT_0168 + 6*IT_0169
      *IT_0170;
    return create_ccomplex_return(IT_0171);
}

