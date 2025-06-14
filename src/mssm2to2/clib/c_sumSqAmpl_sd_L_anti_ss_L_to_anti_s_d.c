#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_ss_L_to_anti_s_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_ss_L_to_anti_s_d(
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
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
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
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W1)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = N_B1*e_em;
    const ccomplex_t IT_0021 = IT_0009*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = N_W1*e_em;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-3)*IT_0025);
    const ccomplex_t IT_0027 = 0.166666666666667*IT_0026;
    const ccomplex_t IT_0028 = IT_0019*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B2)*e_em;
    const ccomplex_t IT_0033 = IT_0009*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W2)*e_em;
    const ccomplex_t IT_0036 = IT_0014*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + (-3)*IT_0037);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = N_B2*e_em;
    const ccomplex_t IT_0041 = IT_0009*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = N_W2*e_em;
    const ccomplex_t IT_0044 = IT_0014*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-3)*IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0039*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = conj(N_B3)*e_em;
    const ccomplex_t IT_0053 = IT_0009*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W3)*e_em;
    const ccomplex_t IT_0056 = IT_0014*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-3)*IT_0057);
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = N_B3*e_em;
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = N_W3*e_em;
    const ccomplex_t IT_0064 = IT_0014*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + (-3)*IT_0065);
    const ccomplex_t IT_0067 = 0.166666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0059*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B4)*e_em;
    const ccomplex_t IT_0073 = IT_0009*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W4)*e_em;
    const ccomplex_t IT_0076 = IT_0014*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + (-3)*IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = N_B4*e_em;
    const ccomplex_t IT_0081 = IT_0009*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W4*e_em;
    const ccomplex_t IT_0084 = IT_0014*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-3)*IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = IT_0079*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.166666666666667*IT_0007 + -IT_0031 + -IT_0051
       + -IT_0071 + -IT_0091;
    const ccomplex_t IT_0093 = pow(m_W, -1);
    const ccomplex_t IT_0094 = cos(beta);
    const ccomplex_t IT_0095 = cpow(IT_0094, -1);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = IT_0049*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = IT_0069*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = IT_0089*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0014*IT_0093*IT_0095;
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = IT_0029*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = -IT_0102 + -IT_0109 + -IT_0116 + -IT_0123;
    const ccomplex_t IT_0125 = pow(m_sd_L, 2);
    const ccomplex_t IT_0126 = m_d*m_s*IT_0125;
    const ccomplex_t IT_0127 = s_34*IT_0125;
    const ccomplex_t IT_0128 = (-18)*IT_0127;
    const ccomplex_t IT_0129 = s_13*s_14;
    const ccomplex_t IT_0130 = 36*IT_0129;
    const ccomplex_t IT_0131 = IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = m_d*m_s;
    const ccomplex_t IT_0133 = m_d*IT_0102;
    const ccomplex_t IT_0134 = m_d*IT_0109;
    const ccomplex_t IT_0135 = m_d*IT_0116;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0137 = IT_0019*IT_0118;
    const ccomplex_t IT_0138 = IT_0029*IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0140 = IT_0039*IT_0097;
    const ccomplex_t IT_0141 = IT_0049*IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0143 = IT_0059*IT_0104;
    const ccomplex_t IT_0144 = IT_0069*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0146 = IT_0079*IT_0111;
    const ccomplex_t IT_0147 = IT_0089*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = m_d*IT_0123;
    const ccomplex_t IT_0149 = IT_0133 + IT_0134 + IT_0135 + -IT_0138 + 
      -IT_0141 + -IT_0144 + -IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = m_d*IT_0007;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = 6*IT_0151;
    const ccomplex_t IT_0153 = IT_0027*IT_0120;
    const ccomplex_t IT_0154 = IT_0029*IT_0136*IT_0153;
    const ccomplex_t IT_0155 = IT_0047*IT_0099;
    const ccomplex_t IT_0156 = IT_0049*IT_0139*IT_0155;
    const ccomplex_t IT_0157 = IT_0067*IT_0106;
    const ccomplex_t IT_0158 = IT_0069*IT_0142*IT_0157;
    const ccomplex_t IT_0159 = IT_0087*IT_0113;
    const ccomplex_t IT_0160 = IT_0089*IT_0145*IT_0159;
    const ccomplex_t IT_0161 = m_d*IT_0031;
    const ccomplex_t IT_0162 = m_d*IT_0051;
    const ccomplex_t IT_0163 = m_d*IT_0071;
    const ccomplex_t IT_0164 = m_d*IT_0091;
    const ccomplex_t IT_0165 = (-0.166666666666667)*IT_0150 + -IT_0154 + 
      -IT_0156 + -IT_0158 + -IT_0160 + IT_0161 + IT_0162 + IT_0163 + IT_0164;
    const ccomplex_t IT_0166 = (-18)*conj(IT_0165);
    const ccomplex_t IT_0167 = IT_0152 + IT_0166;
    const ccomplex_t IT_0168 = (-18)*IT_0165;
    const ccomplex_t IT_0169 = (-6)*IT_0151;
    const ccomplex_t IT_0170 = IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = m_s*s_14;
    const ccomplex_t IT_0172 = 18*s_34;
    const ccomplex_t IT_0173 = m_d*s_13;
    const ccomplex_t IT_0174 = (-0.5)*IT_0007;
    const ccomplex_t IT_0175 = 12*IT_0129;
    const ccomplex_t IT_0176 = (-6)*IT_0127;
    const ccomplex_t IT_0177 = -conj(IT_0165);
    const ccomplex_t IT_0178 = (IT_0124 + -conj(IT_0124))*IT_0126 + (IT_0149 +
       -conj(IT_0149))*IT_0171 + (-0.166666666666667)*IT_0131*IT_0174 + (
      -0.166666666666667)*(IT_0092 + -conj(IT_0092))*(IT_0175 + IT_0176) + (-6)
      *IT_0173*(IT_0151 + 1./6*IT_0165 + 1./6*IT_0177);
    const ccomplex_t IT_0179 = conj(IT_0092)*((-18)*IT_0124*IT_0126 + IT_0092
      *IT_0131) + conj(IT_0124)*((-18)*IT_0092*IT_0126 + IT_0124*IT_0131) +
       IT_0132*(IT_0149*IT_0167 + conj(IT_0149)*IT_0170) + ((-18)*conj(IT_0092)
      *IT_0149 + (-18)*IT_0092*conj(IT_0149) + IT_0124*IT_0167 + conj(IT_0124)
      *IT_0170)*IT_0171 + (IT_0149*conj(IT_0149) + IT_0165*conj(IT_0165))
      *IT_0172 + (18*conj(IT_0124)*IT_0149 + 18*IT_0124*conj(IT_0149) + conj
      (IT_0092)*(IT_0152 + 18*IT_0165) + IT_0092*(18*conj(IT_0165) + IT_0169))
      *IT_0173 + (-18)*s_34*IT_0151*(IT_0151 + 1./3*IT_0165 + 1./3*IT_0177) + 6
      *IT_0174*IT_0178;
    return create_ccomplex_return(IT_0179);
}

