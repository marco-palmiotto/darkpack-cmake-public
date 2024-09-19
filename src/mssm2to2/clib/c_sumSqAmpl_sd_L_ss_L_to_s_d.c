#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_ss_L_to_s_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_ss_L_to_s_d(
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
    const ccomplex_t IT_0008 = m_sG*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_B1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_W1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0014 + (-3)*IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0021, 2);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_N_1*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B4)*e_em;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W4)*e_em;
    const ccomplex_t IT_0030 = IT_0016*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-3)*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0033, 2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*m_N_4*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = cos(beta);
    const ccomplex_t IT_0040 = cpow(IT_0039, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0021*IT_0042;
    const ccomplex_t IT_0044 = IT_0024*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = m_s*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B2)*e_em;
    const ccomplex_t IT_0048 = IT_0011*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = conj(N_W2)*e_em;
    const ccomplex_t IT_0051 = IT_0016*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = m_s*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B3)*e_em;
    const ccomplex_t IT_0063 = IT_0011*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W3)*e_em;
    const ccomplex_t IT_0066 = IT_0016*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + (-3)*IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = m_s*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0033*IT_0078;
    const ccomplex_t IT_0080 = IT_0036*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = m_s*IT_0081;
    const ccomplex_t IT_0083 = (-0.166666666666667)*IT_0008 + IT_0025 +
       IT_0037 + -IT_0046 + -IT_0061 + -IT_0076 + -IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0069, 2);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_N_3*IT_0084;
    const ccomplex_t IT_0086 = IT_0073*IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0054, 2);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_N_2*IT_0087;
    const ccomplex_t IT_0089 = IT_0058*IT_0088;
    const ccomplex_t IT_0090 = -IT_0086 + -IT_0089;
    const ccomplex_t IT_0091 = -IT_0045 + -IT_0060 + -IT_0075 + -IT_0081;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0021*IT_0093;
    const ccomplex_t IT_0095 = IT_0024*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0054*IT_0098;
    const ccomplex_t IT_0100 = IT_0058*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0069*IT_0103;
    const ccomplex_t IT_0105 = IT_0073*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0016*IT_0038*IT_0040;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0033*IT_0108;
    const ccomplex_t IT_0110 = IT_0036*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0096 + -IT_0101 + -IT_0106 + -IT_0111;
    const ccomplex_t IT_0113 = m_d*m_s*IT_0002;
    const ccomplex_t IT_0114 = s_34*IT_0002;
    const ccomplex_t IT_0115 = (-18)*IT_0114;
    const ccomplex_t IT_0116 = s_23*s_24;
    const ccomplex_t IT_0117 = 36*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = m_d*m_s;
    const ccomplex_t IT_0120 = IT_0056*IT_0098;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0122 = IT_0058*IT_0120*IT_0121;
    const ccomplex_t IT_0123 = m_s*IT_0096;
    const ccomplex_t IT_0124 = m_s*IT_0101;
    const ccomplex_t IT_0125 = m_s*IT_0111;
    const ccomplex_t IT_0126 = m_s*IT_0106;
    const ccomplex_t IT_0127 = IT_0122 + -IT_0123 + -IT_0124 + -IT_0125 + 
      -IT_0126;
    const ccomplex_t IT_0128 = IT_0042*IT_0093;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0130 = IT_0024*IT_0128*IT_0129;
    const ccomplex_t IT_0131 = IT_0071*IT_0103;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0133 = IT_0073*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = IT_0078*IT_0108;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0136 = IT_0036*IT_0134*IT_0135;
    const ccomplex_t IT_0137 = -IT_0130 + -IT_0133 + -IT_0136;
    const ccomplex_t IT_0138 = 18*IT_0127;
    const ccomplex_t IT_0139 = (-18)*IT_0137;
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = 18*conj(IT_0127);
    const ccomplex_t IT_0142 = (-18)*conj(IT_0137);
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = 18*IT_0137;
    const ccomplex_t IT_0145 = (-18)*IT_0127;
    const ccomplex_t IT_0146 = IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = 18*conj(IT_0137);
    const ccomplex_t IT_0148 = (-18)*conj(IT_0127);
    const ccomplex_t IT_0149 = IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = 6*IT_0119*(IT_0009*(IT_0127 + -conj(IT_0127) + 
      -IT_0137 + conj(IT_0137)) + 0.166666666666667*conj(IT_0090)*IT_0140 +
       0.166666666666667*IT_0090*IT_0143 + 0.166666666666667*conj(IT_0083)
      *IT_0146 + 0.166666666666667*IT_0083*IT_0149);
    const ccomplex_t IT_0151 = 18*s_34;
    const ccomplex_t IT_0152 = cpow(IT_0009, 2);
    const ccomplex_t IT_0153 = (-18)*s_34;
    const ccomplex_t IT_0154 = (-6)*s_34*IT_0009;
    const ccomplex_t IT_0155 = m_s*s_24;
    const ccomplex_t IT_0156 = 6*IT_0009;
    const ccomplex_t IT_0157 = (-6)*IT_0009;
    const ccomplex_t IT_0158 = m_d*s_23;
    const ccomplex_t IT_0159 = 6*s_34*IT_0009*(conj(IT_0083) + IT_0090) + conj
      (IT_0091)*((-18)*IT_0112*IT_0113 + IT_0091*IT_0118) + conj(IT_0112)*((-18)
      *IT_0091*IT_0113 + IT_0112*IT_0118) + IT_0150 + (IT_0083*conj(IT_0083) +
       IT_0090*conj(IT_0090) + IT_0127*conj(IT_0127) + IT_0137*conj(IT_0137))
      *IT_0151 + (conj(IT_0083)*IT_0090 + IT_0083*conj(IT_0090) + conj(IT_0127)
      *IT_0137 + IT_0127*conj(IT_0137) + IT_0152)*IT_0153 + (IT_0083 + conj
      (IT_0090))*IT_0154 + IT_0155*(conj(IT_0112)*IT_0146 + IT_0112*IT_0149 +
       IT_0091*((-18)*conj(IT_0083) + 18*conj(IT_0090) + IT_0156) + conj(IT_0091
      )*((-18)*IT_0083 + 18*IT_0090 + IT_0157)) + (conj(IT_0091)*IT_0140 +
       IT_0091*IT_0143 + conj(IT_0112)*(18*IT_0083 + (-18)*IT_0090 + IT_0156) +
       IT_0112*(18*conj(IT_0083) + (-18)*conj(IT_0090) + IT_0157))*IT_0158;
    return create_ccomplex_return(IT_0159);
}

