#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_ss_L_to_s_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_ss_L_to_s_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_ss_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = m_s*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0013 = (-0.333333333333333)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = conj(N_B1)*e_em;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = conj(N_W1)*e_em;
    const ccomplex_t IT_0029 = IT_0018*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + (-3)*IT_0030);
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0034 = (-0.333333333333333)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = m_s*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B3)*e_em;
    const ccomplex_t IT_0041 = IT_0011*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W3)*e_em;
    const ccomplex_t IT_0044 = IT_0018*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-3)*IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0049 = (-0.333333333333333)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_s*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B4)*e_em;
    const ccomplex_t IT_0056 = IT_0011*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W4)*e_em;
    const ccomplex_t IT_0059 = IT_0018*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-3)*IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = m_s*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B2)*e_em;
    const ccomplex_t IT_0071 = IT_0011*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W2)*e_em;
    const ccomplex_t IT_0074 = IT_0018*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + (-3)*IT_0075);
    const ccomplex_t IT_0077 = 0.166666666666667*IT_0076;
    const ccomplex_t IT_0078 = IT_0013*IT_0077;
    const ccomplex_t IT_0079 = IT_0023*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = m_s*IT_0080;
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0008 + IT_0024 + -IT_0039 
      + -IT_0054 + -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0034*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0087 = IT_0036*IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0049*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0092 = IT_0051*IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0064*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0097 = IT_0066*IT_0095*IT_0096;
    const ccomplex_t IT_0098 = -IT_0087 + -IT_0092 + -IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0084*IT_0100;
    const ccomplex_t IT_0102 = IT_0036*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = IT_0023*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0089*IT_0110;
    const ccomplex_t IT_0112 = IT_0051*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0094*IT_0115;
    const ccomplex_t IT_0117 = IT_0066*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = -IT_0103 + -IT_0108 + -IT_0113 + -IT_0118;
    const ccomplex_t IT_0120 = 0.166666666666667*IT_0007 + -IT_0038 + -IT_0053
       + -IT_0068 + -IT_0080;
    const ccomplex_t IT_0121 = m_d*m_s*IT_0002;
    const ccomplex_t IT_0122 = s_34*IT_0002;
    const ccomplex_t IT_0123 = (-18)*IT_0122;
    const ccomplex_t IT_0124 = s_23*s_24;
    const ccomplex_t IT_0125 = 36*IT_0124;
    const ccomplex_t IT_0126 = IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = m_s*s_24;
    const ccomplex_t IT_0128 = IT_0032*IT_0100;
    const ccomplex_t IT_0129 = IT_0036*IT_0086*IT_0128;
    const ccomplex_t IT_0130 = IT_0062*IT_0115;
    const ccomplex_t IT_0131 = IT_0066*IT_0096*IT_0130;
    const ccomplex_t IT_0132 = m_s*IT_0103;
    const ccomplex_t IT_0133 = m_s*IT_0108;
    const ccomplex_t IT_0134 = m_s*IT_0113;
    const ccomplex_t IT_0135 = m_s*IT_0118;
    const ccomplex_t IT_0136 = IT_0129 + IT_0131 + -IT_0132 + -IT_0133 + 
      -IT_0134 + -IT_0135;
    const ccomplex_t IT_0137 = IT_0077*IT_0105;
    const ccomplex_t IT_0138 = IT_0022*IT_0023*IT_0137;
    const ccomplex_t IT_0139 = IT_0047*IT_0110;
    const ccomplex_t IT_0140 = IT_0051*IT_0091*IT_0139;
    const ccomplex_t IT_0141 = -IT_0138 + -IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0007;
    const ccomplex_t IT_0143 = 6*IT_0009;
    const ccomplex_t IT_0144 = 18*conj(IT_0098);
    const ccomplex_t IT_0145 = (-18)*conj(IT_0082);
    const ccomplex_t IT_0146 = IT_0143 + IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = 18*IT_0098;
    const ccomplex_t IT_0148 = (-18)*IT_0082;
    const ccomplex_t IT_0149 = (-6)*IT_0009;
    const ccomplex_t IT_0150 = IT_0147 + IT_0148 + IT_0149;
    const ccomplex_t IT_0151 = conj(IT_0119)*((-18)*IT_0136 + 18*IT_0141) +
       IT_0119*((-18)*conj(IT_0136) + 18*conj(IT_0141)) + (36*IT_0009 + 6
      *IT_0082 + (-6)*conj(IT_0082) + (-6)*IT_0098 + 6*conj(IT_0098))*IT_0142 +
       IT_0120*IT_0146 + conj(IT_0120)*IT_0150;
    const ccomplex_t IT_0152 = 18*s_34;
    const ccomplex_t IT_0153 = cpow(IT_0009, 2);
    const ccomplex_t IT_0154 = (-18)*s_34;
    const ccomplex_t IT_0155 = m_d*m_s;
    const ccomplex_t IT_0156 = 18*conj(IT_0082);
    const ccomplex_t IT_0157 = (-18)*conj(IT_0098);
    const ccomplex_t IT_0158 = IT_0149 + IT_0156 + IT_0157;
    const ccomplex_t IT_0159 = 18*IT_0082;
    const ccomplex_t IT_0160 = (-18)*IT_0098;
    const ccomplex_t IT_0161 = (-6)*s_34*IT_0009;
    const ccomplex_t IT_0162 = 12*IT_0124;
    const ccomplex_t IT_0163 = (-6)*IT_0122;
    const ccomplex_t IT_0164 = m_d*s_23;
    const ccomplex_t IT_0165 = -IT_0141;
    const ccomplex_t IT_0166 = IT_0120*(conj(IT_0136) + -conj(IT_0141)) +
       0.0555555555555556*IT_0119*IT_0158 + 0.0555555555555556*conj(IT_0119)*
      (IT_0143 + IT_0159 + IT_0160) + conj(IT_0120)*(IT_0136 + IT_0165) + (
      -0.333333333333333)*IT_0142*(IT_0136 + -conj(IT_0136) + conj(IT_0141) +
       IT_0165);
    const ccomplex_t IT_0167 = 6*s_34*IT_0009*(conj(IT_0082) + IT_0098) + conj
      (IT_0119)*((-18)*IT_0120*IT_0121 + IT_0119*IT_0126) + conj(IT_0120)*((-18)
      *IT_0119*IT_0121 + IT_0120*IT_0126) + IT_0127*IT_0151 + (IT_0082*conj
      (IT_0082) + IT_0098*conj(IT_0098) + IT_0136*conj(IT_0136) + IT_0141*conj
      (IT_0141))*IT_0152 + (conj(IT_0082)*IT_0098 + IT_0082*conj(IT_0098) + conj
      (IT_0136)*IT_0141 + IT_0136*conj(IT_0141) + IT_0153)*IT_0154 + IT_0155*
      (IT_0136*IT_0146 + conj(IT_0136)*IT_0150 + IT_0141*IT_0158 + conj(IT_0141)
      *(IT_0143 + IT_0159 + IT_0160)) + (IT_0082 + conj(IT_0098))*IT_0161 + 6
      *IT_0142*((IT_0119 + -conj(IT_0119))*IT_0121 + (-0.166666666666667)
      *IT_0126*IT_0142 + (-0.166666666666667)*(IT_0120 + -conj(IT_0120))*
      (IT_0162 + IT_0163)) + 18*IT_0164*IT_0166;
    return create_ccomplex_return(IT_0167);
}

