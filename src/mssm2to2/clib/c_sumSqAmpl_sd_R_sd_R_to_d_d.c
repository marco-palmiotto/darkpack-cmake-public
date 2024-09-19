#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sd_R_to_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sd_R_to_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_R = param->m_sd_R;
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
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_d, 2);
    const ccomplex_t IT_0013 = pow(m_sd_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0020 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0016 + IT_0024;
    const ccomplex_t IT_0026 = m_d*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = pow(m_sG, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -IT_0028 +
       reg_prop, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0031 = cpow(IT_0030, 2);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_sG*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = 6*s_34;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -IT_0028 +
       reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = (-18)*s_34;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = IT_0036 + IT_0041;
    const ccomplex_t IT_0043 = 12*s_34;
    const ccomplex_t IT_0044 = cpow(IT_0006, 2);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_N_1*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0018, 2);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_N_2*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = IT_0047 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0056 = (-0.333333333333333)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = m_d*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0035*IT_0062;
    const ccomplex_t IT_0064 = IT_0021*IT_0050;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = m_d*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = 18*s_34;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0057*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = m_d*IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0068*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0079 = (-0.333333333333333)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = m_d*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0068*IT_0085;
    const ccomplex_t IT_0087 = IT_0052*IT_0068;
    const ccomplex_t IT_0088 = cpow(IT_0034, 2);
    const ccomplex_t IT_0089 = cpow(IT_0039, 2);
    const ccomplex_t IT_0090 = cpow(IT_0020, 2);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_2*IT_0090;
    const ccomplex_t IT_0092 = IT_0050*IT_0091;
    const ccomplex_t IT_0093 = IT_0011*IT_0046;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = m_d*IT_0094;
    const ccomplex_t IT_0096 = (-0.166666666666667)*IT_0038 + IT_0092 + 
      -IT_0095;
    const ccomplex_t IT_0097 = cpow(IT_0079, 2);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_3*IT_0097;
    const ccomplex_t IT_0099 = IT_0081*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = cpow(IT_0077, 2);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_N_3*IT_0101;
    const ccomplex_t IT_0103 = IT_0081*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0022*IT_0091;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0098*IT_0107;
    const ccomplex_t IT_0109 = IT_0080*IT_0107;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = m_d*IT_0110;
    const ccomplex_t IT_0112 = IT_0108 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0027 + IT_0112;
    const ccomplex_t IT_0114 = IT_0014*IT_0045;
    const ccomplex_t IT_0115 = cpow(IT_0054, 2);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*m_N_4*IT_0115;
    const ccomplex_t IT_0117 = IT_0058*IT_0116;
    const ccomplex_t IT_0118 = -IT_0114 + -IT_0117;
    const ccomplex_t IT_0119 = IT_0022*IT_0049;
    const ccomplex_t IT_0120 = IT_0102*IT_0107;
    const ccomplex_t IT_0121 = -IT_0111 + IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = m_d*IT_0024;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = m_d*IT_0016;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0121 + IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = conj(IT_0121) + conj(IT_0123) + conj(IT_0125);
    const ccomplex_t IT_0128 = IT_0035*IT_0039;
    const ccomplex_t IT_0129 = IT_0027*IT_0035;
    const ccomplex_t IT_0130 = IT_0128 + IT_0129;
    const ccomplex_t IT_0131 = IT_0035*IT_0112;
    const ccomplex_t IT_0132 = IT_0068*IT_0096;
    const ccomplex_t IT_0133 = IT_0034*IT_0040;
    const ccomplex_t IT_0134 = 36*s_34;
    const ccomplex_t IT_0135 = IT_0067*IT_0134;
    const ccomplex_t IT_0136 = conj(IT_0074) + conj(IT_0085);
    const ccomplex_t IT_0137 = IT_0043*IT_0062;
    const ccomplex_t IT_0138 = IT_0074*IT_0134;
    const ccomplex_t IT_0139 = IT_0035*conj(IT_0067);
    const ccomplex_t IT_0140 = conj(IT_0027)*IT_0068;
    const ccomplex_t IT_0141 = IT_0068*conj(IT_0112);
    const ccomplex_t IT_0142 = IT_0035*conj(IT_0062);
    const ccomplex_t IT_0143 = conj(IT_0067)*IT_0068;
    const ccomplex_t IT_0144 = IT_0068*conj(IT_0074);
    const ccomplex_t IT_0145 = IT_0068*conj(IT_0085);
    const ccomplex_t IT_0146 = (-6)*s_34;
    const ccomplex_t IT_0147 = IT_0039*IT_0068;
    const ccomplex_t IT_0148 = IT_0027*IT_0068;
    const ccomplex_t IT_0149 = IT_0034*IT_0146;
    const ccomplex_t IT_0150 = IT_0035*IT_0067;
    const ccomplex_t IT_0151 = IT_0068*IT_0112;
    const ccomplex_t IT_0152 = conj(IT_0027)*IT_0035;
    const ccomplex_t IT_0153 = IT_0034*IT_0068;
    const ccomplex_t IT_0154 = IT_0039*IT_0146;
    const ccomplex_t IT_0155 = IT_0035*conj(IT_0112);
    const ccomplex_t IT_0156 = IT_0152 + IT_0155;
    const ccomplex_t IT_0157 = IT_0043*conj(IT_0062);
    const ccomplex_t IT_0158 = conj(IT_0067)*IT_0134;
    const ccomplex_t IT_0159 = IT_0040*IT_0067;
    const ccomplex_t IT_0160 = IT_0040*IT_0074;
    const ccomplex_t IT_0161 = IT_0040*IT_0085;
    const ccomplex_t IT_0162 = IT_0062*IT_0146;
    const ccomplex_t IT_0163 = IT_0040*conj(IT_0067);
    const ccomplex_t IT_0164 = IT_0040*conj(IT_0074);
    const ccomplex_t IT_0165 = IT_0040*conj(IT_0085);
    const ccomplex_t IT_0166 = conj(IT_0062)*IT_0146;
    const ccomplex_t IT_0167 = IT_0040*IT_0062;
    const ccomplex_t IT_0168 = IT_0067*IT_0146;
    const ccomplex_t IT_0169 = IT_0074*IT_0146;
    const ccomplex_t IT_0170 = IT_0085*IT_0146;
    const ccomplex_t IT_0171 = IT_0040*conj(IT_0062);
    const ccomplex_t IT_0172 = conj(IT_0067)*IT_0146;
    const ccomplex_t IT_0173 = conj(IT_0074)*IT_0146;
    const ccomplex_t IT_0174 = conj(IT_0085)*IT_0146;
    const ccomplex_t IT_0175 = cpow(IT_0010, 2);
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_N_1*IT_0175;
    const ccomplex_t IT_0177 = IT_0046*IT_0176;
    const ccomplex_t IT_0178 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0179 = IT_0035*IT_0178;
    const ccomplex_t IT_0180 = IT_0068*conj(IT_0096);
    const ccomplex_t IT_0181 = IT_0040*conj(IT_0100);
    const ccomplex_t IT_0182 = conj(IT_0106)*IT_0146;
    const ccomplex_t IT_0183 = cpow(IT_0056, 2);
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_4*IT_0183;
    const ccomplex_t IT_0185 = IT_0070*IT_0184;
    const ccomplex_t IT_0186 = conj(IT_0177) + conj(IT_0185);
    const ccomplex_t IT_0187 = IT_0058*IT_0184;
    const ccomplex_t IT_0188 = IT_0014*IT_0176;
    const ccomplex_t IT_0189 = conj(IT_0187) + conj(IT_0188);
    const ccomplex_t IT_0190 = IT_0035*IT_0074;
    const ccomplex_t IT_0191 = IT_0035*IT_0085;
    const ccomplex_t IT_0192 = IT_0035*IT_0096;
    const ccomplex_t IT_0193 = IT_0035*IT_0177;
    const ccomplex_t IT_0194 = IT_0035*IT_0185;
    const ccomplex_t IT_0195 = IT_0062*IT_0068;
    const ccomplex_t IT_0196 = IT_0068*IT_0187;
    const ccomplex_t IT_0197 = IT_0040*IT_0106;
    const ccomplex_t IT_0198 = IT_0040*IT_0178;
    const ccomplex_t IT_0199 = IT_0100*IT_0146;
    const ccomplex_t IT_0200 = IT_0035*conj(IT_0074);
    const ccomplex_t IT_0201 = IT_0035*conj(IT_0085);
    const ccomplex_t IT_0202 = IT_0035*conj(IT_0096);
    const ccomplex_t IT_0203 = IT_0035*conj(IT_0177);
    const ccomplex_t IT_0204 = IT_0035*conj(IT_0185);
    const ccomplex_t IT_0205 = conj(IT_0062)*IT_0068;
    const ccomplex_t IT_0206 = IT_0068*IT_0178;
    const ccomplex_t IT_0207 = IT_0040*conj(IT_0106);
    const ccomplex_t IT_0208 = conj(IT_0100)*IT_0146;
    const ccomplex_t IT_0209 = IT_0068*IT_0177;
    const ccomplex_t IT_0210 = IT_0040*IT_0100;
    const ccomplex_t IT_0211 = IT_0106*IT_0146;
    const ccomplex_t IT_0212 = IT_0146*IT_0178;
    const ccomplex_t IT_0213 = IT_0035*IT_0052;
    const ccomplex_t IT_0214 = IT_0068*IT_0121;
    const ccomplex_t IT_0215 = IT_0104*IT_0146;
    const ccomplex_t IT_0216 = conj(IT_0123) + conj(IT_0125);
    const ccomplex_t IT_0217 = IT_0068*IT_0123;
    const ccomplex_t IT_0218 = IT_0035*conj(IT_0052);
    const ccomplex_t IT_0219 = IT_0139 + IT_0200 + IT_0201 + IT_0218;
    const ccomplex_t IT_0220 = IT_0205 + IT_0219;
    const ccomplex_t IT_0221 = conj(IT_0104)*IT_0146;
    const ccomplex_t IT_0222 = IT_0068*conj(IT_0121);
    const ccomplex_t IT_0223 = -IT_0095;
    const ccomplex_t IT_0224 = IT_0070*IT_0116;
    const ccomplex_t IT_0225 = conj(IT_0223) + conj(IT_0224);
    const ccomplex_t IT_0226 = IT_0035*IT_0121;
    const ccomplex_t IT_0227 = IT_0035*IT_0123;
    const ccomplex_t IT_0228 = IT_0035*IT_0125;
    const ccomplex_t IT_0229 = IT_0068*IT_0223;
    const ccomplex_t IT_0230 = IT_0040*IT_0104;
    const ccomplex_t IT_0231 = IT_0118*IT_0146;
    const ccomplex_t IT_0232 = IT_0035*conj(IT_0121);
    const ccomplex_t IT_0233 = IT_0035*conj(IT_0123);
    const ccomplex_t IT_0234 = IT_0035*conj(IT_0125);
    const ccomplex_t IT_0235 = conj(IT_0052)*IT_0068;
    const ccomplex_t IT_0236 = IT_0040*conj(IT_0104);
    const ccomplex_t IT_0237 = conj(IT_0118)*IT_0146;
    const ccomplex_t IT_0238 = m_d*s_24;
    const ccomplex_t IT_0239 = IT_0065 + IT_0083;
    const ccomplex_t IT_0240 = 6*IT_0239;
    const ccomplex_t IT_0241 = 6*IT_0072;
    const ccomplex_t IT_0242 = 6*IT_0094;
    const ccomplex_t IT_0243 = -IT_0024;
    const ccomplex_t IT_0244 = 18*IT_0243;
    const ccomplex_t IT_0245 = -IT_0060;
    const ccomplex_t IT_0246 = 18*IT_0245;
    const ccomplex_t IT_0247 = -IT_0110;
    const ccomplex_t IT_0248 = 18*IT_0247;
    const ccomplex_t IT_0249 = -IT_0016;
    const ccomplex_t IT_0250 = 18*IT_0249;
    const ccomplex_t IT_0251 = IT_0187 + IT_0188;
    const ccomplex_t IT_0252 = 6*conj(IT_0239);
    const ccomplex_t IT_0253 = 6*conj(IT_0072);
    const ccomplex_t IT_0254 = 6*conj(IT_0094);
    const ccomplex_t IT_0255 = 18*conj(IT_0243);
    const ccomplex_t IT_0256 = 18*conj(IT_0245);
    const ccomplex_t IT_0257 = 18*conj(IT_0247);
    const ccomplex_t IT_0258 = 18*conj(IT_0249);
    const ccomplex_t IT_0259 = 6*IT_0243;
    const ccomplex_t IT_0260 = 6*IT_0249;
    const ccomplex_t IT_0261 = 18*IT_0094;
    const ccomplex_t IT_0262 = 6*IT_0245;
    const ccomplex_t IT_0263 = 6*IT_0247;
    const ccomplex_t IT_0264 = 18*IT_0239;
    const ccomplex_t IT_0265 = 18*IT_0072;
    const ccomplex_t IT_0266 = 6*conj(IT_0243);
    const ccomplex_t IT_0267 = 6*conj(IT_0249);
    const ccomplex_t IT_0268 = 18*conj(IT_0094);
    const ccomplex_t IT_0269 = IT_0177 + IT_0185;
    const ccomplex_t IT_0270 = 6*conj(IT_0245);
    const ccomplex_t IT_0271 = 6*conj(IT_0247);
    const ccomplex_t IT_0272 = 18*conj(IT_0239);
    const ccomplex_t IT_0273 = 18*conj(IT_0072);
    const ccomplex_t IT_0274 = IT_0072 + IT_0239;
    const ccomplex_t IT_0275 = 6*conj(IT_0027);
    const ccomplex_t IT_0276 = 6*conj(IT_0112);
    const ccomplex_t IT_0277 = 6*conj(IT_0062);
    const ccomplex_t IT_0278 = 6*IT_0178;
    const ccomplex_t IT_0279 = 18*IT_0034;
    const ccomplex_t IT_0280 = 18*conj(IT_0067);
    const ccomplex_t IT_0281 = 18*conj(IT_0074);
    const ccomplex_t IT_0282 = 18*conj(IT_0085);
    const ccomplex_t IT_0283 = 18*conj(IT_0096);
    const ccomplex_t IT_0284 = (-18)*conj(IT_0100);
    const ccomplex_t IT_0285 = (-6)*IT_0039;
    const ccomplex_t IT_0286 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0287 = conj(IT_0245) + conj(IT_0247);
    const ccomplex_t IT_0288 = 6*IT_0067;
    const ccomplex_t IT_0289 = 6*IT_0074;
    const ccomplex_t IT_0290 = 6*IT_0085;
    const ccomplex_t IT_0291 = 6*IT_0096;
    const ccomplex_t IT_0292 = 18*IT_0039;
    const ccomplex_t IT_0293 = 18*IT_0027;
    const ccomplex_t IT_0294 = 18*IT_0112;
    const ccomplex_t IT_0295 = 18*IT_0062;
    const ccomplex_t IT_0296 = (-18)*IT_0106;
    const ccomplex_t IT_0297 = (-18)*IT_0178;
    const ccomplex_t IT_0298 = (-6)*IT_0034;
    const ccomplex_t IT_0299 = (-6)*IT_0100;
    const ccomplex_t IT_0300 = IT_0245 + IT_0247;
    const ccomplex_t IT_0301 = 6*IT_0034;
    const ccomplex_t IT_0302 = 6*conj(IT_0067);
    const ccomplex_t IT_0303 = 6*conj(IT_0074);
    const ccomplex_t IT_0304 = 6*conj(IT_0085);
    const ccomplex_t IT_0305 = 6*conj(IT_0096);
    const ccomplex_t IT_0306 = 18*conj(IT_0027);
    const ccomplex_t IT_0307 = 18*conj(IT_0112);
    const ccomplex_t IT_0308 = 18*conj(IT_0062);
    const ccomplex_t IT_0309 = 18*IT_0178;
    const ccomplex_t IT_0310 = (-18)*IT_0039;
    const ccomplex_t IT_0311 = (-18)*conj(IT_0106);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0313 = conj(IT_0072) + conj(IT_0239);
    const ccomplex_t IT_0314 = 6*IT_0039;
    const ccomplex_t IT_0315 = 6*IT_0027;
    const ccomplex_t IT_0316 = 6*IT_0112;
    const ccomplex_t IT_0317 = 6*IT_0062;
    const ccomplex_t IT_0318 = 18*IT_0067;
    const ccomplex_t IT_0319 = 18*IT_0074;
    const ccomplex_t IT_0320 = 18*IT_0085;
    const ccomplex_t IT_0321 = 18*IT_0096;
    const ccomplex_t IT_0322 = (-18)*IT_0034;
    const ccomplex_t IT_0323 = (-18)*IT_0100;
    const ccomplex_t IT_0324 = (-6)*IT_0106;
    const ccomplex_t IT_0325 = (-6)*IT_0178;
    const ccomplex_t IT_0326 = 6*IT_0052;
    const ccomplex_t IT_0327 = 6*IT_0223;
    const ccomplex_t IT_0328 = 12*IT_0067;
    const ccomplex_t IT_0329 = 12*IT_0074;
    const ccomplex_t IT_0330 = 12*IT_0085;
    const ccomplex_t IT_0331 = 18*IT_0121;
    const ccomplex_t IT_0332 = 18*IT_0123;
    const ccomplex_t IT_0333 = 18*IT_0125;
    const ccomplex_t IT_0334 = 36*IT_0062;
    const ccomplex_t IT_0335 = (-18)*IT_0118;
    const ccomplex_t IT_0336 = (-6)*IT_0104;
    const ccomplex_t IT_0337 = (conj(IT_0243) + conj(IT_0249))*(IT_0291 +
       IT_0292 + IT_0293 + IT_0294 + IT_0296 + IT_0297 + IT_0298 + IT_0299 +
       IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330 + IT_0331 + IT_0332 +
       IT_0333 + IT_0334 + IT_0335 + IT_0336);
    const ccomplex_t IT_0338 = -IT_0060 + -IT_0110;
    const ccomplex_t IT_0339 = IT_0224 + 0.166666666666667*IT_0288 +
       0.166666666666667*IT_0289 + 0.166666666666667*IT_0290 + 0.166666666666667
      *IT_0295 + 0.166666666666667*IT_0326 + 0.166666666666667*IT_0327 +
       0.166666666666667*IT_0331 + 0.166666666666667*IT_0332 + 0.166666666666667
      *IT_0333 + 0.166666666666667*IT_0335 + 0.166666666666667*IT_0336;
    const ccomplex_t IT_0340 = IT_0243 + IT_0249;
    const ccomplex_t IT_0341 = 6*conj(IT_0052);
    const ccomplex_t IT_0342 = 6*conj(IT_0223);
    const ccomplex_t IT_0343 = 12*conj(IT_0067);
    const ccomplex_t IT_0344 = 12*conj(IT_0074);
    const ccomplex_t IT_0345 = 12*conj(IT_0085);
    const ccomplex_t IT_0346 = 18*conj(IT_0121);
    const ccomplex_t IT_0347 = 18*conj(IT_0123);
    const ccomplex_t IT_0348 = 18*conj(IT_0125);
    const ccomplex_t IT_0349 = 36*conj(IT_0062);
    const ccomplex_t IT_0350 = (-18)*conj(IT_0118);
    const ccomplex_t IT_0351 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0352 = conj(IT_0224) + 0.166666666666667*IT_0302 +
       0.166666666666667*IT_0303 + 0.166666666666667*IT_0304 + 0.166666666666667
      *IT_0308 + 0.166666666666667*IT_0341 + 0.166666666666667*IT_0342 +
       0.166666666666667*IT_0346 + 0.166666666666667*IT_0347 + 0.166666666666667
      *IT_0348 + 0.166666666666667*IT_0350 + 0.166666666666667*IT_0351;
    const ccomplex_t IT_0353 = IT_0065 + IT_0072 + IT_0083;
    const ccomplex_t IT_0354 = 6*IT_0121;
    const ccomplex_t IT_0355 = 6*IT_0123;
    const ccomplex_t IT_0356 = 6*IT_0125;
    const ccomplex_t IT_0357 = 18*IT_0052;
    const ccomplex_t IT_0358 = 18*IT_0223;
    const ccomplex_t IT_0359 = (-18)*IT_0104;
    const ccomplex_t IT_0360 = (-6)*IT_0118;
    const ccomplex_t IT_0361 = 12*IT_0062 + 36*IT_0067 + 36*IT_0074 + 36
      *IT_0085 + IT_0314 + IT_0315 + IT_0316 + IT_0321 + IT_0322 + IT_0323 +
       IT_0324 + IT_0325 + IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 +
       IT_0359 + IT_0360;
    const ccomplex_t IT_0362 = 6*conj(IT_0121);
    const ccomplex_t IT_0363 = 6*conj(IT_0123);
    const ccomplex_t IT_0364 = 6*conj(IT_0125);
    const ccomplex_t IT_0365 = 18*conj(IT_0052);
    const ccomplex_t IT_0366 = 18*conj(IT_0223);
    const ccomplex_t IT_0367 = (-18)*conj(IT_0104);
    const ccomplex_t IT_0368 = (-6)*conj(IT_0118);
    const ccomplex_t IT_0369 = 12*conj(IT_0062) + 36*conj(IT_0067) + 36*conj
      (IT_0074) + 36*conj(IT_0085) + IT_0275 + IT_0276 + IT_0278 + IT_0279 +
       IT_0283 + IT_0284 + IT_0285 + IT_0286 + IT_0362 + IT_0363 + IT_0364 +
       IT_0365 + IT_0366 + IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = IT_0189*(IT_0240 + IT_0241 + IT_0242 + IT_0244 
      + IT_0246 + IT_0248 + IT_0250) + IT_0251*(IT_0252 + IT_0253 + IT_0254 +
       IT_0255 + IT_0256 + IT_0257 + IT_0258) + conj(IT_0224)*(IT_0259 + IT_0260
       + IT_0261) + IT_0186*(IT_0259 + IT_0260 + IT_0261 + IT_0262 + IT_0263 +
       IT_0264 + IT_0265) + IT_0224*(IT_0266 + IT_0267 + IT_0268) + IT_0269*
      (IT_0266 + IT_0267 + IT_0268 + IT_0270 + IT_0271 + IT_0272 + IT_0273) +
       IT_0274*(IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280 +
       IT_0281 + IT_0282 + IT_0283 + IT_0284 + IT_0285 + IT_0286) + IT_0287*
      (IT_0288 + IT_0289 + IT_0290 + IT_0291 + IT_0292 + IT_0293 + IT_0294 +
       IT_0295 + IT_0296 + IT_0297 + IT_0298 + IT_0299) + IT_0300*(IT_0301 +
       IT_0302 + IT_0303 + IT_0304 + IT_0305 + IT_0306 + IT_0307 + IT_0308 +
       IT_0309 + IT_0310 + IT_0311 + IT_0312) + IT_0313*(IT_0314 + IT_0315 +
       IT_0316 + IT_0317 + IT_0318 + IT_0319 + IT_0320 + IT_0321 + IT_0322 +
       IT_0323 + IT_0324 + IT_0325) + IT_0337 + 6*conj(IT_0338)*IT_0339 +
       IT_0340*(IT_0301 + IT_0305 + IT_0306 + IT_0307 + IT_0309 + IT_0310 +
       IT_0311 + IT_0312 + IT_0341 + IT_0342 + IT_0343 + IT_0344 + IT_0345 +
       IT_0346 + IT_0347 + IT_0348 + IT_0349 + IT_0350 + IT_0351) + 6*IT_0338
      *IT_0352 + conj(IT_0353)*(18*IT_0224 + IT_0317 + IT_0318 + IT_0319 +
       IT_0320 + IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 + IT_0359 +
       IT_0360) + conj(IT_0094)*IT_0361 + IT_0353*(18*conj(IT_0224) + IT_0277 +
       IT_0280 + IT_0281 + IT_0282 + IT_0362 + IT_0363 + IT_0364 + IT_0365 +
       IT_0366 + IT_0367 + IT_0368) + IT_0094*IT_0369;
    const ccomplex_t IT_0371 = IT_0302 + IT_0303 + IT_0304 + IT_0308;
    const ccomplex_t IT_0372 = (-36)*IT_0062;
    const ccomplex_t IT_0373 = (-18)*IT_0027;
    const ccomplex_t IT_0374 = (-18)*IT_0112;
    const ccomplex_t IT_0375 = (-12)*IT_0067;
    const ccomplex_t IT_0376 = (-18)*conj(IT_0027);
    const ccomplex_t IT_0377 = (-18)*conj(IT_0112);
    const ccomplex_t IT_0378 = (-12)*conj(IT_0067);
    const ccomplex_t IT_0379 = (-36)*IT_0067;
    const ccomplex_t IT_0380 = IT_0279 + IT_0379;
    const ccomplex_t IT_0381 = (-6)*IT_0027;
    const ccomplex_t IT_0382 = (-6)*conj(IT_0027);
    const ccomplex_t IT_0383 = (-6)*IT_0067;
    const ccomplex_t IT_0384 = (-6)*conj(IT_0067);
    const ccomplex_t IT_0385 = (-36)*IT_0074;
    const ccomplex_t IT_0386 = (-12)*IT_0062;
    const ccomplex_t IT_0387 = (-6)*IT_0112;
    const ccomplex_t IT_0388 = (-36)*IT_0085;
    const ccomplex_t IT_0389 = (-36)*conj(IT_0067);
    const ccomplex_t IT_0390 = (-12)*conj(IT_0062);
    const ccomplex_t IT_0391 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0392 = (-36)*conj(IT_0074);
    const ccomplex_t IT_0393 = (-18)*IT_0067;
    const ccomplex_t IT_0394 = (-18)*IT_0074;
    const ccomplex_t IT_0395 = (-18)*IT_0085;
    const ccomplex_t IT_0396 = (-6)*IT_0062;
    const ccomplex_t IT_0397 = 6*IT_0106;
    const ccomplex_t IT_0398 = 18*IT_0100;
    const ccomplex_t IT_0399 = (-18)*IT_0096;
    const ccomplex_t IT_0400 = conj(IT_0224)*((-18)*IT_0177 + (-18)*IT_0185 + 
      (-6)*IT_0187 + (-6)*IT_0188 + IT_0278 + IT_0279 + IT_0285 + IT_0381 +
       IT_0387 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 + IT_0398 +
       IT_0399);
    const ccomplex_t IT_0401 = (-18)*conj(IT_0067);
    const ccomplex_t IT_0402 = (-18)*conj(IT_0074);
    const ccomplex_t IT_0403 = (-18)*conj(IT_0085);
    const ccomplex_t IT_0404 = (-6)*conj(IT_0062);
    const ccomplex_t IT_0405 = (-18)*IT_0062;
    const ccomplex_t IT_0406 = (-6)*IT_0074;
    const ccomplex_t IT_0407 = (-6)*IT_0085;
    const ccomplex_t IT_0408 = 6*IT_0100;
    const ccomplex_t IT_0409 = 18*IT_0106;
    const ccomplex_t IT_0410 = (-6)*IT_0096;
    const ccomplex_t IT_0411 = 6*conj(IT_0100);
    const ccomplex_t IT_0412 = 18*conj(IT_0106);
    const ccomplex_t IT_0413 = (-18)*conj(IT_0062);
    const ccomplex_t IT_0414 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0415 = (-6)*conj(IT_0085);
    const ccomplex_t IT_0416 = (-6)*conj(IT_0096);
    const ccomplex_t IT_0417 = 6*conj(IT_0106);
    const ccomplex_t IT_0418 = 18*conj(IT_0100);
    const ccomplex_t IT_0419 = (-18)*conj(IT_0096);
    const ccomplex_t IT_0420 = IT_0224*((-18)*conj(IT_0177) + (-18)*conj
      (IT_0185) + (-6)*conj(IT_0187) + (-6)*conj(IT_0188) + IT_0314 + IT_0322 +
       IT_0325 + IT_0382 + IT_0391 + IT_0401 + IT_0402 + IT_0403 + IT_0404 +
       IT_0417 + IT_0418 + IT_0419);
    const ccomplex_t IT_0421 = 6*IT_0118;
    const ccomplex_t IT_0422 = 18*IT_0104;
    const ccomplex_t IT_0423 = (-18)*IT_0052;
    const ccomplex_t IT_0424 = (-18)*IT_0223;
    const ccomplex_t IT_0425 = (-6)*IT_0121;
    const ccomplex_t IT_0426 = (-6)*IT_0123;
    const ccomplex_t IT_0427 = (-6)*IT_0125;
    const ccomplex_t IT_0428 = 6*conj(IT_0118);
    const ccomplex_t IT_0429 = 18*conj(IT_0104);
    const ccomplex_t IT_0430 = (-18)*conj(IT_0052);
    const ccomplex_t IT_0431 = (-18)*conj(IT_0223);
    const ccomplex_t IT_0432 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0433 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0434 = (-6)*conj(IT_0125);
    const ccomplex_t IT_0435 = 6*IT_0104;
    const ccomplex_t IT_0436 = 18*IT_0118;
    const ccomplex_t IT_0437 = (-18)*IT_0121;
    const ccomplex_t IT_0438 = (-18)*IT_0123;
    const ccomplex_t IT_0439 = (-18)*IT_0125;
    const ccomplex_t IT_0440 = (-6)*IT_0052;
    const ccomplex_t IT_0441 = (-6)*IT_0223;
    const ccomplex_t IT_0442 = 6*conj(IT_0104);
    const ccomplex_t IT_0443 = 18*conj(IT_0118);
    const ccomplex_t IT_0444 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0445 = (-18)*conj(IT_0123);
    const ccomplex_t IT_0446 = (-18)*conj(IT_0125);
    const ccomplex_t IT_0447 = (-6)*conj(IT_0052);
    const ccomplex_t IT_0448 = (-6)*conj(IT_0223);
    const ccomplex_t IT_0449 = IT_0100*(IT_0277 + IT_0280 + IT_0281 + IT_0282)
       + IT_0104*(IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280 +
       IT_0281 + IT_0282 + IT_0283 + IT_0284 + IT_0285 + IT_0286) + conj(IT_0106
      )*(IT_0288 + IT_0289 + IT_0290 + IT_0295) + conj(IT_0118)*(IT_0288 +
       IT_0289 + IT_0290 + IT_0291 + IT_0292 + IT_0293 + IT_0294 + IT_0295 +
       IT_0296 + IT_0297 + IT_0298 + IT_0299) + IT_0118*(IT_0301 + IT_0302 +
       IT_0303 + IT_0304 + IT_0305 + IT_0306 + IT_0307 + IT_0308 + IT_0309 +
       IT_0310 + IT_0311 + IT_0312) + conj(IT_0100)*(IT_0317 + IT_0318 + IT_0319
       + IT_0320) + conj(IT_0104)*(IT_0314 + IT_0315 + IT_0316 + IT_0317 +
       IT_0318 + IT_0319 + IT_0320 + IT_0321 + IT_0322 + IT_0323 + IT_0324 +
       IT_0325) + IT_0106*IT_0371 + conj(IT_0062)*(IT_0301 + IT_0310 + IT_0372 +
       IT_0373 + IT_0374 + IT_0375) + IT_0062*(IT_0292 + IT_0298 + IT_0376 +
       IT_0377 + IT_0378) + conj(IT_0067)*(IT_0285 + IT_0380 + IT_0381) +
       IT_0067*(IT_0314 + IT_0322 + IT_0382) + conj(IT_0112)*IT_0383 + IT_0112
      *IT_0384 + conj(IT_0074)*(IT_0285 + IT_0380 + IT_0381 + IT_0385 + IT_0386 
      + IT_0387) + conj(IT_0085)*(IT_0285 + IT_0380 + IT_0381 + IT_0385 +
       IT_0386 + IT_0387 + IT_0388) + IT_0074*(IT_0314 + IT_0322 + IT_0382 +
       IT_0389 + IT_0390 + IT_0391) + IT_0085*(IT_0314 + IT_0322 + IT_0382 +
       IT_0389 + IT_0390 + IT_0391 + IT_0392) + conj(IT_0096)*(IT_0393 + IT_0394
       + IT_0395 + IT_0396) + (conj(IT_0052) + conj(IT_0223))*(IT_0278 + IT_0279
       + IT_0285 + IT_0381 + IT_0387 + IT_0393 + IT_0394 + IT_0395 + IT_0396 +
       IT_0397 + IT_0398 + IT_0399) + IT_0400 + IT_0096*(IT_0401 + IT_0402 +
       IT_0403 + IT_0404) + IT_0178*(IT_0371 + IT_0383 + IT_0405 + IT_0406 +
       IT_0407) + IT_0127*(IT_0301 + IT_0309 + IT_0310 + IT_0373 + IT_0374 +
       IT_0383 + IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409 + IT_0410) +
       IT_0126*(IT_0292 + IT_0297 + IT_0298 + IT_0376 + IT_0377 + IT_0384 +
       IT_0411 + IT_0412 + IT_0413 + IT_0414 + IT_0415 + IT_0416) + (IT_0052 +
       IT_0223)*(IT_0314 + IT_0322 + IT_0325 + IT_0382 + IT_0391 + IT_0401 +
       IT_0402 + IT_0403 + IT_0404 + IT_0417 + IT_0418 + IT_0419) + IT_0420 +
       IT_0186*(IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0421 + IT_0422 +
       IT_0423 + IT_0424 + IT_0425 + IT_0426 + IT_0427) + IT_0269*(IT_0401 +
       IT_0402 + IT_0403 + IT_0404 + IT_0428 + IT_0429 + IT_0430 + IT_0431 +
       IT_0432 + IT_0433 + IT_0434) + IT_0189*(IT_0383 + IT_0405 + IT_0406 +
       IT_0407 + IT_0435 + IT_0436 + IT_0437 + IT_0438 + IT_0439 + IT_0440 +
       IT_0441) + IT_0251*(IT_0384 + IT_0413 + IT_0414 + IT_0415 + IT_0442 +
       IT_0443 + IT_0444 + IT_0445 + IT_0446 + IT_0447 + IT_0448);
    const ccomplex_t IT_0450 = IT_0012*IT_0013;
    const ccomplex_t IT_0451 = (-36)*IT_0243;
    const ccomplex_t IT_0452 = (-18)*IT_0243;
    const ccomplex_t IT_0453 = (-18)*conj(IT_0243);
    const ccomplex_t IT_0454 = (-6)*IT_0243;
    const ccomplex_t IT_0455 = (-6)*conj(IT_0243);
    const ccomplex_t IT_0456 = (-18)*IT_0239;
    const ccomplex_t IT_0457 = (-6)*IT_0245;
    const ccomplex_t IT_0458 = (-6)*IT_0247;
    const ccomplex_t IT_0459 = (-6)*IT_0249;
    const ccomplex_t IT_0460 = (-18)*conj(IT_0239);
    const ccomplex_t IT_0461 = (-6)*conj(IT_0245);
    const ccomplex_t IT_0462 = (-6)*conj(IT_0247);
    const ccomplex_t IT_0463 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0464 = (-18)*IT_0353;
    const ccomplex_t IT_0465 = (-18)*IT_0072;
    const ccomplex_t IT_0466 = (-6)*IT_0338;
    const ccomplex_t IT_0467 = (-18)*conj(IT_0353);
    const ccomplex_t IT_0468 = (-18)*conj(IT_0072);
    const ccomplex_t IT_0469 = (-6)*conj(IT_0338);
    const ccomplex_t IT_0470 = (-18)*IT_0245;
    const ccomplex_t IT_0471 = (-18)*IT_0247;
    const ccomplex_t IT_0472 = (-6)*IT_0239;
    const ccomplex_t IT_0473 = (-18)*IT_0249;
    const ccomplex_t IT_0474 = (-18)*conj(IT_0245);
    const ccomplex_t IT_0475 = (-18)*conj(IT_0247);
    const ccomplex_t IT_0476 = (-6)*conj(IT_0239);
    const ccomplex_t IT_0477 = (-18)*conj(IT_0249);
    const ccomplex_t IT_0478 = conj(IT_0243)*IT_0451 + IT_0287*IT_0452 +
       IT_0300*IT_0453 + conj(IT_0239)*IT_0454 + IT_0239*IT_0455 + conj(IT_0353)
      *(IT_0454 + IT_0456 + IT_0457 + IT_0458 + IT_0459) + IT_0353*(IT_0455 +
       IT_0460 + IT_0461 + IT_0462 + IT_0463) + conj(IT_0072)*(IT_0454 + IT_0459
       + IT_0464 + IT_0466) + conj(IT_0094)*((-12)*IT_0243 + (-12)*IT_0249 +
       IT_0456 + IT_0457 + IT_0458 + IT_0464 + IT_0465 + IT_0466) + IT_0072*
      (IT_0455 + IT_0463 + IT_0467 + IT_0469) + IT_0094*((-36)*conj(IT_0094) + (
      -12)*conj(IT_0243) + (-12)*conj(IT_0249) + IT_0460 + IT_0461 + IT_0462 +
       IT_0467 + IT_0468 + IT_0469) + conj(IT_0249)*(IT_0451 + IT_0470 + IT_0471
       + IT_0472) + conj(IT_0338)*(IT_0452 + IT_0470 + IT_0471 + IT_0472 +
       IT_0473) + IT_0249*((-36)*conj(IT_0243) + (-36)*conj(IT_0249) + IT_0474 +
       IT_0475 + IT_0476) + IT_0338*(IT_0453 + IT_0474 + IT_0475 + IT_0476 +
       IT_0477);
    const ccomplex_t IT_0479 = m_d*s_23;
    const ccomplex_t IT_0480 = IT_0094*((-36)*conj(IT_0085) + IT_0314 +
       IT_0322 + IT_0325 + IT_0382 + IT_0389 + IT_0390 + IT_0391 + IT_0392 +
       IT_0417 + IT_0418 + IT_0419 + IT_0428 + IT_0429 + IT_0430 + IT_0431 +
       IT_0432 + IT_0433 + IT_0434);
    const ccomplex_t IT_0481 = (-12)*IT_0074;
    const ccomplex_t IT_0482 = (-12)*IT_0085;
    const ccomplex_t IT_0483 = (conj(IT_0243) + conj(IT_0249))*(IT_0301 +
       IT_0309 + IT_0310 + IT_0372 + IT_0373 + IT_0374 + IT_0375 + IT_0408 +
       IT_0409 + IT_0410 + IT_0435 + IT_0436 + IT_0437 + IT_0438 + IT_0439 +
       IT_0440 + IT_0441 + IT_0481 + IT_0482);
    const ccomplex_t IT_0484 = (-36)*conj(IT_0062);
    const ccomplex_t IT_0485 = (-12)*conj(IT_0074);
    const ccomplex_t IT_0486 = (-12)*conj(IT_0085);
    const ccomplex_t IT_0487 = (-18)*IT_0094;
    const ccomplex_t IT_0488 = (-18)*conj(IT_0094);
    const ccomplex_t IT_0489 = (-18)*IT_0338;
    const ccomplex_t IT_0490 = (-6)*IT_0353;
    const ccomplex_t IT_0491 = (-6)*IT_0094;
    const ccomplex_t IT_0492 = (-18)*conj(IT_0338);
    const ccomplex_t IT_0493 = (-6)*conj(IT_0353);
    const ccomplex_t IT_0494 = (-6)*conj(IT_0094);
    const ccomplex_t IT_0495 = conj(IT_0353)*(IT_0278 + IT_0279 + IT_0285 +
       IT_0381 + IT_0387 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 +
       IT_0398 + IT_0399) + conj(IT_0338)*(IT_0301 + IT_0309 + IT_0310 + IT_0373
       + IT_0374 + IT_0383 + IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409 +
       IT_0410) + IT_0338*(IT_0292 + IT_0297 + IT_0298 + IT_0376 + IT_0377 +
       IT_0384 + IT_0411 + IT_0412 + IT_0413 + IT_0414 + IT_0415 + IT_0416) +
       IT_0353*(IT_0314 + IT_0322 + IT_0325 + IT_0382 + IT_0391 + IT_0401 +
       IT_0402 + IT_0403 + IT_0404 + IT_0417 + IT_0418 + IT_0419) + IT_0313*
      (IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0421 + IT_0422 + IT_0423 +
       IT_0424 + IT_0425 + IT_0426 + IT_0427) + conj(IT_0094)*(IT_0278 + IT_0279
       + IT_0285 + IT_0379 + IT_0381 + IT_0385 + IT_0386 + IT_0387 + IT_0388 +
       IT_0397 + IT_0398 + IT_0399 + IT_0421 + IT_0422 + IT_0423 + IT_0424 +
       IT_0425 + IT_0426 + IT_0427) + IT_0274*(IT_0401 + IT_0402 + IT_0403 +
       IT_0404 + IT_0428 + IT_0429 + IT_0430 + IT_0431 + IT_0432 + IT_0433 +
       IT_0434) + IT_0287*(IT_0383 + IT_0405 + IT_0406 + IT_0407 + IT_0435 +
       IT_0436 + IT_0437 + IT_0438 + IT_0439 + IT_0440 + IT_0441) + IT_0300*
      (IT_0384 + IT_0413 + IT_0414 + IT_0415 + IT_0442 + IT_0443 + IT_0444 +
       IT_0445 + IT_0446 + IT_0447 + IT_0448) + IT_0480 + IT_0483 + IT_0340*
      (IT_0292 + IT_0297 + IT_0298 + IT_0376 + IT_0377 + IT_0378 + IT_0411 +
       IT_0412 + IT_0416 + IT_0442 + IT_0443 + IT_0444 + IT_0445 + IT_0446 +
       IT_0447 + IT_0448 + IT_0484 + IT_0485 + IT_0486) + conj(IT_0224)*(IT_0454
       + IT_0456 + IT_0457 + IT_0458 + IT_0459 + IT_0465 + IT_0487) + IT_0186*
      (IT_0454 + IT_0459 + IT_0464 + IT_0466 + IT_0487) + IT_0224*(IT_0455 +
       IT_0460 + IT_0461 + IT_0462 + IT_0463 + IT_0468 + IT_0488) + IT_0269*
      (IT_0455 + IT_0463 + IT_0467 + IT_0469 + IT_0488) + IT_0189*(IT_0452 +
       IT_0473 + IT_0489 + IT_0490 + IT_0491) + IT_0251*(IT_0453 + IT_0477 +
       IT_0492 + IT_0493 + IT_0494);
    const ccomplex_t IT_0496 = s_34*IT_0013;
    const ccomplex_t IT_0497 = (-6)*IT_0496;
    const ccomplex_t IT_0498 = s_23*s_24;
    const ccomplex_t IT_0499 = 12*IT_0498;
    const ccomplex_t IT_0500 = IT_0497 + IT_0499;
    const ccomplex_t IT_0501 = 2*conj(IT_0243);
    const ccomplex_t IT_0502 = 2*IT_0249;
    const ccomplex_t IT_0503 = 0.5*IT_0338;
    const ccomplex_t IT_0504 = (conj(IT_0243) + conj(IT_0245) + conj(IT_0247) 
      + conj(IT_0249))*IT_0274 + (IT_0243 + IT_0245 + IT_0247 + IT_0249)*IT_0313
       + (conj(IT_0243) + conj(IT_0249) + conj(IT_0338))*IT_0353 + (IT_0338 +
       IT_0340)*conj(IT_0353) + 2*IT_0094*(conj(IT_0249) + 0.5*IT_0287 + 0.5
      *conj(IT_0338) + 0.5*IT_0501) + 2*conj(IT_0094)*(IT_0243 + 0.5*IT_0245 +
       0.5*IT_0247 + 0.5*IT_0502 + IT_0503);
    const ccomplex_t IT_0505 = (-18)*IT_0496;
    const ccomplex_t IT_0506 = 36*IT_0498;
    const ccomplex_t IT_0507 = IT_0505 + IT_0506;
    const ccomplex_t IT_0508 = 0.5*(conj(IT_0243) + conj(IT_0245) + conj
      (IT_0247) + conj(IT_0249))*IT_0300 + 0.5*IT_0274*(conj(IT_0094) + IT_0313)
       + IT_0243*(conj(IT_0243) + 0.5*conj(IT_0245) + 0.5*conj(IT_0247) + conj
      (IT_0249) + 0.5*conj(IT_0338)) + 0.5*IT_0353*(conj(IT_0094) + conj(IT_0353
      )) + 0.5*IT_0094*(conj(IT_0072) + 2*conj(IT_0094) + conj(IT_0239) + conj
      (IT_0353)) + 0.5*IT_0249*(IT_0287 + conj(IT_0338) + IT_0501) + 0.5*conj
      (IT_0249)*IT_0502 + (conj(IT_0243) + conj(IT_0249) + conj(IT_0338))*IT_0503;
    const ccomplex_t IT_0509 = conj(IT_0027)*IT_0040;
    const ccomplex_t IT_0510 = IT_0040*conj(IT_0112);
    const ccomplex_t IT_0511 = conj(IT_0096)*IT_0146;
    const ccomplex_t IT_0512 = IT_0068*conj(IT_0106);
    const ccomplex_t IT_0513 = (-12)*s_34*IT_0034 + IT_0035*conj(IT_0100) +
       IT_0039*IT_0134 + IT_0148 + IT_0150 + IT_0151 + IT_0171 + IT_0172 +
       IT_0173 + IT_0174 + IT_0040*IT_0178 + IT_0190 + IT_0191 + IT_0192 +
       IT_0195 + IT_0197 + IT_0199 + IT_0509 + IT_0510 + IT_0511 + IT_0512;
    const ccomplex_t IT_0514 = IT_0027*IT_0042 + IT_0034*IT_0039*IT_0043 +
       conj(IT_0052)*(IT_0063 + IT_0069 + IT_0075 + IT_0086 + IT_0087) + IT_0040
      *(IT_0088 + IT_0089 + conj(IT_0096)*IT_0100 + IT_0096*conj(IT_0100) + conj
      (IT_0052)*IT_0104 + IT_0052*conj(IT_0104) + conj(IT_0106)*IT_0113 + conj
      (IT_0118)*IT_0126 + IT_0118*IT_0127) + conj(IT_0096)*(IT_0063 + IT_0069 +
       IT_0075 + IT_0086 + IT_0130 + IT_0131 + IT_0132 + IT_0133) + conj(IT_0067
      )*(IT_0130 + IT_0133 + IT_0135) + IT_0136*(IT_0130 + IT_0131 + IT_0133 +
       IT_0135 + IT_0137 + IT_0138) + IT_0112*(IT_0036 + IT_0041 + IT_0139 +
       IT_0140) + IT_0062*(IT_0036 + IT_0041 + IT_0043*conj(IT_0067) + conj
      (IT_0062)*IT_0134 + IT_0140 + IT_0141) + IT_0052*(IT_0142 + IT_0143 +
       IT_0144 + IT_0145 + conj(IT_0118)*IT_0146) + conj(IT_0027)*(IT_0147 +
       IT_0148 + IT_0149) + conj(IT_0062)*(IT_0043*IT_0067 + IT_0147 + IT_0148 +
       IT_0149 + IT_0151) + conj(IT_0112)*(IT_0147 + IT_0148 + IT_0149 + IT_0150
       + IT_0151) + IT_0067*(IT_0152 + IT_0153 + IT_0154) + IT_0096*(IT_0142 +
       IT_0143 + IT_0144 + IT_0145 + conj(IT_0106)*IT_0146 + IT_0153 + IT_0154 +
       IT_0156) + IT_0074*(IT_0153 + IT_0154 + IT_0156 + IT_0157 + IT_0158) +
       IT_0085*(IT_0134*IT_0136 + IT_0153 + IT_0154 + IT_0156 + IT_0157 +
       IT_0158) + conj(IT_0104)*(IT_0159 + IT_0160 + IT_0161 + IT_0162) + conj
      (IT_0100)*(IT_0113*IT_0146 + IT_0153 + IT_0154 + IT_0159 + IT_0160 +
       IT_0161 + IT_0162) + IT_0104*(IT_0068*conj(IT_0104) + IT_0035*conj
      (IT_0118) + IT_0163 + IT_0164 + IT_0165 + IT_0166) + IT_0100*(IT_0068*conj
      (IT_0100) + IT_0035*conj(IT_0106) + IT_0128 + IT_0133 + (conj(IT_0027) +
       conj(IT_0112))*IT_0146 + IT_0163 + IT_0164 + IT_0165 + IT_0166) + conj
      (IT_0118)*(IT_0167 + IT_0168 + IT_0169 + IT_0170) + conj(IT_0106)*(IT_0042
       + IT_0167 + IT_0168 + IT_0169 + IT_0170) + IT_0118*(IT_0035*conj(IT_0104)
       + IT_0068*conj(IT_0118) + conj(IT_0052)*IT_0146 + IT_0171 + IT_0172 +
       IT_0173 + IT_0174) + IT_0177*(IT_0142 + IT_0143 + IT_0144 + IT_0145 +
       IT_0153 + IT_0154 + IT_0156 + IT_0179 + IT_0180 + IT_0181 + IT_0182) +
       IT_0185*(IT_0142 + IT_0143 + IT_0144 + IT_0145 + IT_0153 + IT_0154 +
       IT_0156 + IT_0179 + IT_0180 + IT_0181 + IT_0182 + IT_0068*IT_0186) +
       IT_0189*(IT_0147 + IT_0148 + IT_0149 + IT_0150 + IT_0151 + IT_0190 +
       IT_0191 + IT_0192 + IT_0193 + IT_0194 + IT_0195 + IT_0196 + IT_0197 +
       IT_0198 + IT_0199) + IT_0187*(IT_0036 + IT_0041 + IT_0139 + IT_0140 +
       IT_0141 + IT_0200 + IT_0201 + IT_0202 + IT_0203 + IT_0204 + IT_0205 +
       IT_0206 + IT_0207 + IT_0208) + IT_0188*(IT_0036 + IT_0041 + IT_0139 +
       IT_0140 + IT_0141 + IT_0068*IT_0189 + IT_0200 + IT_0201 + IT_0202 +
       IT_0203 + IT_0204 + IT_0205 + IT_0206 + IT_0207 + IT_0208) + IT_0186*
      (IT_0063 + IT_0069 + IT_0075 + IT_0086 + IT_0130 + IT_0131 + IT_0132 +
       IT_0133 + IT_0209 + IT_0210 + IT_0211 + IT_0212) + conj(IT_0121)*(IT_0150
       + IT_0190 + IT_0191 + IT_0195 + IT_0213 + IT_0214 + IT_0215) + IT_0216*
      (IT_0150 + IT_0190 + IT_0191 + IT_0195 + IT_0213 + IT_0214 + IT_0215 +
       IT_0217) + IT_0121*(IT_0220 + IT_0221) + IT_0125*(IT_0205 + IT_0068
      *IT_0216 + IT_0219 + IT_0221 + IT_0222) + IT_0123*(IT_0220 + IT_0221 +
       IT_0222) + IT_0225*(IT_0063 + IT_0069 + IT_0075 + IT_0086 + IT_0087 +
       IT_0226 + IT_0227 + IT_0228 + IT_0229 + IT_0230 + IT_0231) + IT_0223*
      (IT_0142 + IT_0143 + IT_0144 + IT_0145 + IT_0232 + IT_0233 + IT_0234 +
       IT_0235 + IT_0236 + IT_0237) + IT_0224*(IT_0142 + IT_0143 + IT_0144 +
       IT_0145 + IT_0068*IT_0225 + IT_0232 + IT_0233 + IT_0234 + IT_0235 +
       IT_0236 + IT_0237) + IT_0238*IT_0370 + IT_0012*IT_0449 + IT_0450*IT_0478 
      + IT_0479*IT_0495 + IT_0500*IT_0504 + 2*IT_0507*IT_0508 + IT_0106*(IT_0035
      *conj(IT_0100) + IT_0147 + IT_0149 + IT_0171 + IT_0172 + IT_0173 + IT_0174
       + IT_0509 + IT_0510 + IT_0511 + IT_0512) + IT_0178*IT_0513;
    return create_ccomplex_return(IT_0514);
}

