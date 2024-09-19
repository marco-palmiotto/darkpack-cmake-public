#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_ss_L_to_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_ss_L_to_s_s(
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
    const ccomplex_t IT_0001 = pow(m_ss_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*m_N_1*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_N_2*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_N_3*IT_0024;
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_N_4*IT_0030;
    const ccomplex_t IT_0032 = IT_0027*IT_0031;
    const ccomplex_t IT_0033 = -IT_0026 + -IT_0032;
    const ccomplex_t IT_0034 = cos(theta_W);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = conj(N_B3)*e_em;
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W3)*e_em;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + (-3)*IT_0041);
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = IT_0023*IT_0043;
    const ccomplex_t IT_0045 = IT_0021*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B4)*e_em;
    const ccomplex_t IT_0048 = IT_0035*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = conj(N_W4)*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0029*IT_0054;
    const ccomplex_t IT_0056 = IT_0027*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B2)*e_em;
    const ccomplex_t IT_0059 = IT_0035*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W2)*e_em;
    const ccomplex_t IT_0062 = IT_0007*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + (-3)*IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = IT_0016*IT_0065;
    const ccomplex_t IT_0067 = IT_0014*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0046 + IT_0057 + IT_0068;
    const ccomplex_t IT_0070 = m_s*IT_0069;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = conj(N_B1)*e_em;
    const ccomplex_t IT_0073 = IT_0035*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em;
    const ccomplex_t IT_0076 = IT_0007*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + (-3)*IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = IT_0009*IT_0079;
    const ccomplex_t IT_0081 = IT_0002*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = m_s*IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0080*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0066*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0044*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0087 + IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = m_s*IT_0094;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0055*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = m_s*IT_0099;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0087 + IT_0090;
    const ccomplex_t IT_0103 = m_s*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = m_s*IT_0093;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = m_s*IT_0068;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = m_s*IT_0057;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = m_s*IT_0046;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = 18*conj(IT_0108) + 18*conj(IT_0110) + 18*conj
      (IT_0112);
    const ccomplex_t IT_0114 = IT_0104 + IT_0106;
    const ccomplex_t IT_0115 = IT_0031*IT_0097;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0018*IT_0088;
    const ccomplex_t IT_0118 = IT_0025*IT_0091;
    const ccomplex_t IT_0119 = -IT_0117 + -IT_0118;
    const ccomplex_t IT_0120 = IT_0011*IT_0085;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0116 + IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = cpow(IT_0079, 2);
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_N_1*IT_0123;
    const ccomplex_t IT_0125 = IT_0002*IT_0124;
    const ccomplex_t IT_0126 = cpow(IT_0065, 2);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_N_2*IT_0126;
    const ccomplex_t IT_0128 = IT_0014*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0130 = cpow(IT_0129, 2);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_sG*IT_0130;
    const ccomplex_t IT_0132 = pow(m_sG, 2);
    const ccomplex_t IT_0133 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0132 +
       reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0131*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -IT_0132 +
       reg_prop, -1);
    const ccomplex_t IT_0136 = IT_0131*IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0043, 2);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_N_3*IT_0137;
    const ccomplex_t IT_0139 = IT_0021*IT_0138;
    const ccomplex_t IT_0140 = -IT_0125 + -IT_0128 + (-0.5)*IT_0134 +
       0.166666666666667*IT_0136 + -IT_0139;
    const ccomplex_t IT_0141 = cpow(IT_0054, 2);
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_4*IT_0141;
    const ccomplex_t IT_0143 = IT_0027*IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0140 + IT_0144;
    const ccomplex_t IT_0146 = IT_0097*IT_0142;
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = IT_0085*IT_0124;
    const ccomplex_t IT_0149 = IT_0088*IT_0127;
    const ccomplex_t IT_0150 = IT_0091*IT_0138;
    const ccomplex_t IT_0151 = 0.166666666666667*IT_0134 + (-0.5)*IT_0136 + 
      -IT_0148 + -IT_0149 + -IT_0150;
    const ccomplex_t IT_0152 = IT_0147 + IT_0151;
    const ccomplex_t IT_0153 = (IT_0013 + IT_0020 + IT_0033)*(18*conj(IT_0013)
       + 18*conj(IT_0020) + 18*conj(IT_0033)) + IT_0071*(18*conj(IT_0071) + 18
      *conj(IT_0084)) + IT_0096*(18*conj(IT_0096) + 18*conj(IT_0101)) + IT_0101*
      (18*conj(IT_0096) + 18*conj(IT_0104) + 18*conj(IT_0106)) + IT_0084*(18
      *conj(IT_0071) + IT_0113) + (IT_0108 + IT_0110 + IT_0112)*(18*conj(IT_0084
      ) + IT_0113) + (18*conj(IT_0101) + 18*conj(IT_0104) + 18*conj(IT_0106))
      *IT_0114 + (18*conj(IT_0116) + 18*conj(IT_0119) + 18*conj(IT_0121))
      *IT_0122 + (18*conj(IT_0140) + 18*conj(IT_0144))*IT_0145 + (18*conj
      (IT_0147) + 18*conj(IT_0151))*IT_0152;
    const ccomplex_t IT_0154 = s_34*IT_0001;
    const ccomplex_t IT_0155 = (-18)*IT_0154;
    const ccomplex_t IT_0156 = s_23*s_24;
    const ccomplex_t IT_0157 = 36*IT_0156;
    const ccomplex_t IT_0158 = IT_0155 + IT_0157;
    const ccomplex_t IT_0159 = -IT_0090;
    const ccomplex_t IT_0160 = -IT_0099;
    const ccomplex_t IT_0161 = -IT_0087 + -IT_0093;
    const ccomplex_t IT_0162 = IT_0159 + IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = conj(IT_0159) + conj(IT_0160) + conj(IT_0161);
    const ccomplex_t IT_0164 = IT_0046 + IT_0068;
    const ccomplex_t IT_0165 = IT_0057 + IT_0082 + IT_0164;
    const ccomplex_t IT_0166 = conj(IT_0057) + conj(IT_0082) + conj(IT_0164);
    const ccomplex_t IT_0167 = -IT_0090 + -IT_0099;
    const ccomplex_t IT_0168 = -IT_0093;
    const ccomplex_t IT_0169 = -IT_0087;
    const ccomplex_t IT_0170 = IT_0167 + IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = conj(IT_0167) + conj(IT_0168) + conj(IT_0169);
    const ccomplex_t IT_0172 = IT_0046 + IT_0057;
    const ccomplex_t IT_0173 = IT_0068 + IT_0082;
    const ccomplex_t IT_0174 = conj(IT_0172) + conj(IT_0173);
    const ccomplex_t IT_0175 = (-6)*IT_0154;
    const ccomplex_t IT_0176 = 12*IT_0156;
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = IT_0172 + IT_0173;
    const ccomplex_t IT_0179 = 6*s_34;
    const ccomplex_t IT_0180 = IT_0013*IT_0179;
    const ccomplex_t IT_0181 = IT_0020*IT_0179;
    const ccomplex_t IT_0182 = IT_0033*IT_0179;
    const ccomplex_t IT_0183 = (-18)*s_34;
    const ccomplex_t IT_0184 = IT_0104*IT_0183;
    const ccomplex_t IT_0185 = IT_0101*IT_0183;
    const ccomplex_t IT_0186 = IT_0106*IT_0183;
    const ccomplex_t IT_0187 = IT_0180 + IT_0181 + IT_0182 + IT_0184 + IT_0185
       + IT_0186;
    const ccomplex_t IT_0188 = (-6)*s_34;
    const ccomplex_t IT_0189 = IT_0108*IT_0188;
    const ccomplex_t IT_0190 = IT_0110*IT_0188;
    const ccomplex_t IT_0191 = IT_0187 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = conj(IT_0013)*IT_0179;
    const ccomplex_t IT_0193 = conj(IT_0020)*IT_0179;
    const ccomplex_t IT_0194 = conj(IT_0033)*IT_0179;
    const ccomplex_t IT_0195 = conj(IT_0104)*IT_0183;
    const ccomplex_t IT_0196 = conj(IT_0101)*IT_0183;
    const ccomplex_t IT_0197 = conj(IT_0106)*IT_0183;
    const ccomplex_t IT_0198 = IT_0192 + IT_0193 + IT_0194 + IT_0195 + IT_0196
       + IT_0197;
    const ccomplex_t IT_0199 = conj(IT_0108)*IT_0188;
    const ccomplex_t IT_0200 = conj(IT_0110)*IT_0188;
    const ccomplex_t IT_0201 = IT_0198 + IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = conj(IT_0147) + conj(IT_0151);
    const ccomplex_t IT_0203 = IT_0101*IT_0179;
    const ccomplex_t IT_0204 = conj(IT_0101)*IT_0179;
    const ccomplex_t IT_0205 = IT_0104*IT_0179;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = IT_0013*IT_0183;
    const ccomplex_t IT_0208 = IT_0020*IT_0183;
    const ccomplex_t IT_0209 = conj(IT_0104)*IT_0179;
    const ccomplex_t IT_0210 = IT_0204 + IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0013)*IT_0183;
    const ccomplex_t IT_0212 = conj(IT_0020)*IT_0183;
    const ccomplex_t IT_0213 = IT_0106*IT_0179;
    const ccomplex_t IT_0214 = IT_0033*IT_0183;
    const ccomplex_t IT_0215 = conj(IT_0106)*IT_0179;
    const ccomplex_t IT_0216 = conj(IT_0033)*IT_0183;
    const ccomplex_t IT_0217 = conj(IT_0140) + conj(IT_0144);
    const ccomplex_t IT_0218 = IT_0013*IT_0188;
    const ccomplex_t IT_0219 = conj(IT_0013)*IT_0188;
    const ccomplex_t IT_0220 = IT_0104*IT_0188;
    const ccomplex_t IT_0221 = conj(IT_0104)*IT_0188;
    const ccomplex_t IT_0222 = IT_0020*IT_0188;
    const ccomplex_t IT_0223 = 36*s_34;
    const ccomplex_t IT_0224 = conj(IT_0020)*IT_0188;
    const ccomplex_t IT_0225 = IT_0084*IT_0183;
    const ccomplex_t IT_0226 = IT_0071*IT_0183;
    const ccomplex_t IT_0227 = IT_0101*IT_0188;
    const ccomplex_t IT_0228 = conj(IT_0084)*IT_0183;
    const ccomplex_t IT_0229 = conj(IT_0071)*IT_0183;
    const ccomplex_t IT_0230 = conj(IT_0101)*IT_0188;
    const ccomplex_t IT_0231 = IT_0116*IT_0188;
    const ccomplex_t IT_0232 = 12*s_34;
    const ccomplex_t IT_0233 = conj(IT_0116)*IT_0188;
    const ccomplex_t IT_0234 = IT_0112*IT_0188;
    const ccomplex_t IT_0235 = conj(IT_0112)*IT_0188;
    const ccomplex_t IT_0236 = IT_0084*IT_0188;
    const ccomplex_t IT_0237 = conj(IT_0084)*IT_0188;
    const ccomplex_t IT_0238 = IT_0140*IT_0179;
    const ccomplex_t IT_0239 = IT_0071*IT_0188;
    const ccomplex_t IT_0240 = conj(IT_0140)*IT_0179;
    const ccomplex_t IT_0241 = conj(IT_0071)*IT_0188;
    const ccomplex_t IT_0242 = 6*IT_0101;
    const ccomplex_t IT_0243 = 6*conj(IT_0101);
    const ccomplex_t IT_0244 = 6*IT_0084;
    const ccomplex_t IT_0245 = 6*conj(IT_0084);
    const ccomplex_t IT_0246 = 18*IT_0101;
    const ccomplex_t IT_0247 = 18*conj(IT_0101);
    const ccomplex_t IT_0248 = 6*IT_0013;
    const ccomplex_t IT_0249 = 6*IT_0020;
    const ccomplex_t IT_0250 = (-18)*IT_0104;
    const ccomplex_t IT_0251 = 6*conj(IT_0013);
    const ccomplex_t IT_0252 = 6*conj(IT_0020);
    const ccomplex_t IT_0253 = (-18)*conj(IT_0104);
    const ccomplex_t IT_0254 = (-18)*IT_0101;
    const ccomplex_t IT_0255 = (-18)*conj(IT_0101);
    const ccomplex_t IT_0256 = (-6)*IT_0101;
    const ccomplex_t IT_0257 = (-6)*conj(IT_0101);
    const ccomplex_t IT_0258 = 6*IT_0116;
    const ccomplex_t IT_0259 = 6*IT_0119;
    const ccomplex_t IT_0260 = 18*IT_0013;
    const ccomplex_t IT_0261 = 18*IT_0020;
    const ccomplex_t IT_0262 = 18*IT_0033;
    const ccomplex_t IT_0263 = (-18)*IT_0108;
    const ccomplex_t IT_0264 = (-18)*IT_0110;
    const ccomplex_t IT_0265 = (-18)*IT_0112;
    const ccomplex_t IT_0266 = (-6)*IT_0104;
    const ccomplex_t IT_0267 = (-6)*IT_0106;
    const ccomplex_t IT_0268 = 6*IT_0121;
    const ccomplex_t IT_0269 = (-18)*IT_0084;
    const ccomplex_t IT_0270 = 6*conj(IT_0116);
    const ccomplex_t IT_0271 = 6*conj(IT_0119);
    const ccomplex_t IT_0272 = 18*conj(IT_0013);
    const ccomplex_t IT_0273 = 18*conj(IT_0020);
    const ccomplex_t IT_0274 = 18*conj(IT_0033);
    const ccomplex_t IT_0275 = (-18)*conj(IT_0108);
    const ccomplex_t IT_0276 = (-18)*conj(IT_0110);
    const ccomplex_t IT_0277 = (-18)*conj(IT_0112);
    const ccomplex_t IT_0278 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0279 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0280 = 6*conj(IT_0121);
    const ccomplex_t IT_0281 = (-18)*conj(IT_0084);
    const ccomplex_t IT_0282 = 6*IT_0108;
    const ccomplex_t IT_0283 = 6*IT_0110;
    const ccomplex_t IT_0284 = 6*IT_0112;
    const ccomplex_t IT_0285 = 18*IT_0104;
    const ccomplex_t IT_0286 = 18*IT_0106;
    const ccomplex_t IT_0287 = (-18)*IT_0116;
    const ccomplex_t IT_0288 = (-18)*IT_0119;
    const ccomplex_t IT_0289 = (-18)*IT_0121;
    const ccomplex_t IT_0290 = (-6)*IT_0013;
    const ccomplex_t IT_0291 = (-6)*IT_0020;
    const ccomplex_t IT_0292 = (-6)*IT_0033;
    const ccomplex_t IT_0293 = 6*conj(IT_0108);
    const ccomplex_t IT_0294 = 6*conj(IT_0110);
    const ccomplex_t IT_0295 = 6*conj(IT_0112);
    const ccomplex_t IT_0296 = 18*conj(IT_0104);
    const ccomplex_t IT_0297 = 18*conj(IT_0106);
    const ccomplex_t IT_0298 = (-18)*conj(IT_0116);
    const ccomplex_t IT_0299 = (-18)*conj(IT_0119);
    const ccomplex_t IT_0300 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0301 = (-6)*conj(IT_0013);
    const ccomplex_t IT_0302 = (-6)*conj(IT_0020);
    const ccomplex_t IT_0303 = (-6)*conj(IT_0033);
    const ccomplex_t IT_0304 = 6*IT_0033;
    const ccomplex_t IT_0305 = 18*IT_0116;
    const ccomplex_t IT_0306 = 18*IT_0119;
    const ccomplex_t IT_0307 = 18*IT_0121;
    const ccomplex_t IT_0308 = (-18)*IT_0106;
    const ccomplex_t IT_0309 = (-6)*IT_0108;
    const ccomplex_t IT_0310 = (-6)*IT_0110;
    const ccomplex_t IT_0311 = (-6)*IT_0112;
    const ccomplex_t IT_0312 = (-6)*IT_0084;
    const ccomplex_t IT_0313 = 6*conj(IT_0033);
    const ccomplex_t IT_0314 = 18*conj(IT_0116);
    const ccomplex_t IT_0315 = 18*conj(IT_0119);
    const ccomplex_t IT_0316 = 18*conj(IT_0121);
    const ccomplex_t IT_0317 = (-18)*conj(IT_0106);
    const ccomplex_t IT_0318 = (-6)*conj(IT_0108);
    const ccomplex_t IT_0319 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0320 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0321 = (-6)*conj(IT_0084);
    const ccomplex_t IT_0322 = 6*IT_0104;
    const ccomplex_t IT_0323 = 6*IT_0106;
    const ccomplex_t IT_0324 = 18*IT_0108;
    const ccomplex_t IT_0325 = 18*IT_0110;
    const ccomplex_t IT_0326 = 18*IT_0112;
    const ccomplex_t IT_0327 = 18*IT_0084;
    const ccomplex_t IT_0328 = (-18)*IT_0013;
    const ccomplex_t IT_0329 = (-18)*IT_0020;
    const ccomplex_t IT_0330 = (-18)*IT_0033;
    const ccomplex_t IT_0331 = (-6)*IT_0116;
    const ccomplex_t IT_0332 = (-6)*IT_0119;
    const ccomplex_t IT_0333 = (-6)*IT_0121;
    const ccomplex_t IT_0334 = 6*conj(IT_0104);
    const ccomplex_t IT_0335 = 6*conj(IT_0106);
    const ccomplex_t IT_0336 = 18*conj(IT_0108);
    const ccomplex_t IT_0337 = 18*conj(IT_0110);
    const ccomplex_t IT_0338 = 18*conj(IT_0112);
    const ccomplex_t IT_0339 = 18*conj(IT_0084);
    const ccomplex_t IT_0340 = (-18)*conj(IT_0013);
    const ccomplex_t IT_0341 = (-18)*conj(IT_0020);
    const ccomplex_t IT_0342 = (-18)*conj(IT_0033);
    const ccomplex_t IT_0343 = (-6)*conj(IT_0116);
    const ccomplex_t IT_0344 = (-6)*conj(IT_0119);
    const ccomplex_t IT_0345 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0346 = conj(IT_0033)*IT_0242 + IT_0033*IT_0243 + conj
      (IT_0121)*IT_0244 + IT_0121*IT_0245 + (conj(IT_0116) + conj(IT_0119) +
       conj(IT_0121))*IT_0246 + IT_0122*IT_0247 + conj(IT_0101)*(IT_0248 +
       IT_0249 + IT_0250) + IT_0101*((-36)*conj(IT_0101) + IT_0251 + IT_0252 +
       IT_0253) + conj(IT_0106)*IT_0254 + IT_0106*IT_0255 + (conj(IT_0108) +
       conj(IT_0110) + conj(IT_0112))*IT_0256 + (IT_0108 + IT_0110 + IT_0112)
      *IT_0257 + conj(IT_0084)*((-12)*IT_0101 + IT_0258 + IT_0259 + IT_0260 +
       IT_0261 + IT_0262 + IT_0263 + IT_0264 + IT_0265 + IT_0266 + IT_0267) +
       conj(IT_0071)*(IT_0256 + IT_0258 + IT_0259 + IT_0260 + IT_0261 + IT_0262 
      + IT_0263 + IT_0264 + IT_0265 + IT_0266 + IT_0267 + IT_0268 + IT_0269) +
       IT_0084*((-36)*conj(IT_0084) + (-12)*conj(IT_0101) + IT_0270 + IT_0271 +
       IT_0272 + IT_0273 + IT_0274 + IT_0275 + IT_0276 + IT_0277 + IT_0278 +
       IT_0279) + IT_0071*(IT_0257 + IT_0270 + IT_0271 + IT_0272 + IT_0273 +
       IT_0274 + IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280 +
       IT_0281) + IT_0202*(IT_0244 + IT_0246 + IT_0282 + IT_0283 + IT_0284 +
       IT_0285 + IT_0286 + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291 +
       IT_0292) + IT_0152*(IT_0245 + IT_0247 + IT_0293 + IT_0294 + IT_0295 +
       IT_0296 + IT_0297 + IT_0298 + IT_0299 + IT_0300 + IT_0301 + IT_0302 +
       IT_0303) + conj(IT_0096)*(IT_0248 + IT_0249 + IT_0250 + IT_0254 + IT_0304
       + IT_0305 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311 +
       IT_0312) + IT_0096*(IT_0251 + IT_0252 + IT_0253 + IT_0255 + IT_0313 +
       IT_0314 + IT_0315 + IT_0316 + IT_0317 + IT_0318 + IT_0319 + IT_0320 +
       IT_0321) + IT_0217*(IT_0242 + IT_0322 + IT_0323 + IT_0324 + IT_0325 +
       IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330 + IT_0331 + IT_0332 +
       IT_0333) + IT_0145*(IT_0243 + IT_0334 + IT_0335 + IT_0336 + IT_0337 +
       IT_0338 + IT_0339 + IT_0340 + IT_0341 + IT_0342 + IT_0343 + IT_0344 +
       IT_0345);
    const ccomplex_t IT_0347 = m_s*s_23;
    const ccomplex_t IT_0348 = 6*IT_0140;
    const ccomplex_t IT_0349 = 6*IT_0144;
    const ccomplex_t IT_0350 = 18*IT_0147;
    const ccomplex_t IT_0351 = 18*IT_0151;
    const ccomplex_t IT_0352 = (-18)*IT_0096;
    const ccomplex_t IT_0353 = (-6)*IT_0071;
    const ccomplex_t IT_0354 = 6*conj(IT_0140);
    const ccomplex_t IT_0355 = 6*conj(IT_0144);
    const ccomplex_t IT_0356 = 18*conj(IT_0147);
    const ccomplex_t IT_0357 = 18*conj(IT_0151);
    const ccomplex_t IT_0358 = (-18)*conj(IT_0096);
    const ccomplex_t IT_0359 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0360 = 6*IT_0147;
    const ccomplex_t IT_0361 = 6*IT_0151;
    const ccomplex_t IT_0362 = 18*IT_0140;
    const ccomplex_t IT_0363 = 18*IT_0144;
    const ccomplex_t IT_0364 = (-18)*IT_0071;
    const ccomplex_t IT_0365 = (-6)*IT_0096;
    const ccomplex_t IT_0366 = 6*conj(IT_0147);
    const ccomplex_t IT_0367 = 6*conj(IT_0151);
    const ccomplex_t IT_0368 = 18*conj(IT_0140);
    const ccomplex_t IT_0369 = 18*conj(IT_0144);
    const ccomplex_t IT_0370 = (-18)*conj(IT_0071);
    const ccomplex_t IT_0371 = (-6)*conj(IT_0096);
    const ccomplex_t IT_0372 = IT_0174*(IT_0256 + IT_0258 + IT_0259 + IT_0260 
      + IT_0261 + IT_0262 + IT_0263 + IT_0264 + IT_0265 + IT_0266 + IT_0267 +
       IT_0268 + IT_0269) + IT_0178*(IT_0257 + IT_0270 + IT_0271 + IT_0272 +
       IT_0273 + IT_0274 + IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279 +
       IT_0280 + IT_0281) + IT_0171*(IT_0248 + IT_0249 + IT_0250 + IT_0254 +
       IT_0304 + IT_0305 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 +
       IT_0311 + IT_0312) + IT_0170*(IT_0251 + IT_0252 + IT_0253 + IT_0255 +
       IT_0313 + IT_0314 + IT_0315 + IT_0316 + IT_0317 + IT_0318 + IT_0319 +
       IT_0320 + IT_0321) + IT_0163*(IT_0254 + IT_0312 + IT_0348 + IT_0349 +
       IT_0350 + IT_0351 + IT_0352 + IT_0353) + IT_0162*(IT_0255 + IT_0321 +
       IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 + IT_0359) + IT_0166*
      (IT_0256 + IT_0269 + IT_0360 + IT_0361 + IT_0362 + IT_0363 + IT_0364 +
       IT_0365) + IT_0165*(IT_0257 + IT_0281 + IT_0366 + IT_0367 + IT_0368 +
       IT_0369 + IT_0370 + IT_0371);
    const ccomplex_t IT_0373 = m_s*s_24;
    const ccomplex_t IT_0374 = 6*IT_0071;
    const ccomplex_t IT_0375 = 18*IT_0096;
    const ccomplex_t IT_0376 = (-18)*IT_0147;
    const ccomplex_t IT_0377 = (-18)*IT_0151;
    const ccomplex_t IT_0378 = (-6)*IT_0140;
    const ccomplex_t IT_0379 = (-6)*IT_0144;
    const ccomplex_t IT_0380 = 6*conj(IT_0071);
    const ccomplex_t IT_0381 = 18*conj(IT_0096);
    const ccomplex_t IT_0382 = (-18)*conj(IT_0147);
    const ccomplex_t IT_0383 = (-18)*conj(IT_0151);
    const ccomplex_t IT_0384 = (-6)*conj(IT_0140);
    const ccomplex_t IT_0385 = (-6)*conj(IT_0144);
    const ccomplex_t IT_0386 = 6*IT_0096;
    const ccomplex_t IT_0387 = 18*IT_0071;
    const ccomplex_t IT_0388 = (-18)*IT_0140;
    const ccomplex_t IT_0389 = (-18)*IT_0144;
    const ccomplex_t IT_0390 = (-6)*IT_0147;
    const ccomplex_t IT_0391 = (-6)*IT_0151;
    const ccomplex_t IT_0392 = 6*conj(IT_0096);
    const ccomplex_t IT_0393 = 18*conj(IT_0071);
    const ccomplex_t IT_0394 = (-18)*conj(IT_0140);
    const ccomplex_t IT_0395 = (-18)*conj(IT_0144);
    const ccomplex_t IT_0396 = (-6)*conj(IT_0147);
    const ccomplex_t IT_0397 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0398 = IT_0163*(IT_0244 + IT_0246 + IT_0282 + IT_0283 
      + IT_0284 + IT_0285 + IT_0286 + IT_0287 + IT_0288 + IT_0289 + IT_0290 +
       IT_0291 + IT_0292) + IT_0162*(IT_0245 + IT_0247 + IT_0293 + IT_0294 +
       IT_0295 + IT_0296 + IT_0297 + IT_0298 + IT_0299 + IT_0300 + IT_0301 +
       IT_0302 + IT_0303) + IT_0166*(IT_0242 + IT_0322 + IT_0323 + IT_0324 +
       IT_0325 + IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330 + IT_0331 +
       IT_0332 + IT_0333) + IT_0165*(IT_0243 + IT_0334 + IT_0335 + IT_0336 +
       IT_0337 + IT_0338 + IT_0339 + IT_0340 + IT_0341 + IT_0342 + IT_0343 +
       IT_0344 + IT_0345) + IT_0171*(IT_0244 + IT_0246 + IT_0374 + IT_0375 +
       IT_0376 + IT_0377 + IT_0378 + IT_0379) + IT_0170*(IT_0245 + IT_0247 +
       IT_0380 + IT_0381 + IT_0382 + IT_0383 + IT_0384 + IT_0385) + IT_0174*
      (IT_0242 + IT_0327 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 +
       IT_0391) + IT_0178*(IT_0243 + IT_0339 + IT_0392 + IT_0393 + IT_0394 +
       IT_0395 + IT_0396 + IT_0397);
    const ccomplex_t IT_0399 = IT_0000*IT_0001;
    const ccomplex_t IT_0400 = (-18)*IT_0159;
    const ccomplex_t IT_0401 = (-18)*IT_0160;
    const ccomplex_t IT_0402 = (-18)*IT_0161;
    const ccomplex_t IT_0403 = (-6)*IT_0164;
    const ccomplex_t IT_0404 = (-18)*conj(IT_0159);
    const ccomplex_t IT_0405 = (-18)*conj(IT_0160);
    const ccomplex_t IT_0406 = (-18)*conj(IT_0161);
    const ccomplex_t IT_0407 = (-6)*conj(IT_0164);
    const ccomplex_t IT_0408 = (-18)*IT_0164;
    const ccomplex_t IT_0409 = (-6)*IT_0159;
    const ccomplex_t IT_0410 = (-6)*IT_0160;
    const ccomplex_t IT_0411 = (-6)*IT_0161;
    const ccomplex_t IT_0412 = (-18)*conj(IT_0164);
    const ccomplex_t IT_0413 = (-6)*conj(IT_0159);
    const ccomplex_t IT_0414 = (-6)*conj(IT_0160);
    const ccomplex_t IT_0415 = (-6)*conj(IT_0161);
    const ccomplex_t IT_0416 = (-18)*IT_0172;
    const ccomplex_t IT_0417 = (-18)*IT_0173;
    const ccomplex_t IT_0418 = (-6)*IT_0167;
    const ccomplex_t IT_0419 = (-6)*IT_0168;
    const ccomplex_t IT_0420 = (-6)*IT_0169;
    const ccomplex_t IT_0421 = (-18)*conj(IT_0172);
    const ccomplex_t IT_0422 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0423 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0424 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0425 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0426 = IT_0171*(IT_0400 + IT_0401 + IT_0402 + IT_0403)
       + IT_0170*(IT_0404 + IT_0405 + IT_0406 + IT_0407) + IT_0174*(IT_0408 +
       IT_0409 + IT_0410 + IT_0411) + IT_0178*(IT_0412 + IT_0413 + IT_0414 +
       IT_0415) + (conj(IT_0057) + conj(IT_0082))*(IT_0416 + IT_0417 + IT_0418 +
       IT_0419 + IT_0420) + (IT_0057 + IT_0082)*(IT_0421 + IT_0422 + IT_0423 +
       IT_0424 + IT_0425);
    const ccomplex_t IT_0427 = s_34*IT_0153 + IT_0158*(IT_0162*IT_0163 +
       IT_0165*IT_0166 + IT_0170*IT_0171 + (IT_0172 + IT_0173)*IT_0174) +
       IT_0177*(IT_0163*IT_0165 + IT_0162*IT_0166 + IT_0170*IT_0174 + IT_0171
      *IT_0178) + conj(IT_0116)*IT_0191 + IT_0116*IT_0201 + IT_0183*(conj
      (IT_0033)*IT_0108 + IT_0033*conj(IT_0108) + conj(IT_0096)*IT_0152 +
       IT_0096*IT_0202) + IT_0179*((conj(IT_0071) + conj(IT_0084))*IT_0096 + 
      (IT_0071 + IT_0084)*conj(IT_0096) + conj(IT_0106)*IT_0108 + IT_0106*conj
      (IT_0108) + conj(IT_0144)*IT_0152 + IT_0144*IT_0202) + conj(IT_0071)
      *IT_0203 + IT_0071*IT_0204 + conj(IT_0108)*(IT_0206 + IT_0207 + IT_0208) +
       IT_0108*(IT_0210 + IT_0211 + IT_0212) + conj(IT_0110)*(IT_0206 + IT_0207 
      + IT_0208 + IT_0213 + IT_0214) + IT_0110*(IT_0210 + IT_0211 + IT_0212 +
       IT_0215 + IT_0216) + IT_0188*(conj(IT_0033)*IT_0106 + IT_0033*conj
      (IT_0106) + conj(IT_0096)*IT_0145 + IT_0096*IT_0217) + (conj(IT_0101) +
       conj(IT_0104) + conj(IT_0106))*IT_0218 + IT_0114*IT_0219 + (conj(IT_0020)
       + conj(IT_0033))*IT_0220 + (IT_0020 + IT_0033)*IT_0221 + (conj(IT_0101) +
       conj(IT_0106))*IT_0222 + IT_0101*(IT_0219 + conj(IT_0101)*IT_0223 +
       IT_0224) + IT_0106*IT_0224 + IT_0217*(IT_0225 + IT_0226 + IT_0227) + conj
      (IT_0033)*IT_0227 + IT_0145*(IT_0228 + IT_0229 + IT_0230) + IT_0033
      *IT_0230 + conj(IT_0112)*(IT_0206 + IT_0207 + IT_0208 + IT_0213 + IT_0214 
      + IT_0231) + conj(IT_0084)*(IT_0205 + IT_0207 + IT_0208 + IT_0213 +
       IT_0214 + IT_0084*IT_0223 + IT_0231 + IT_0101*IT_0232) + IT_0112*(IT_0210
       + IT_0211 + IT_0212 + IT_0215 + IT_0216 + IT_0233) + IT_0084*(IT_0209 +
       IT_0211 + IT_0212 + IT_0215 + IT_0216 + conj(IT_0101)*IT_0232 + IT_0233) 
      + conj(IT_0119)*(IT_0187 + IT_0084*IT_0188 + IT_0189 + IT_0190 + IT_0234) 
      + IT_0119*(conj(IT_0084)*IT_0188 + IT_0198 + IT_0199 + IT_0200 + IT_0235) 
      + conj(IT_0121)*(IT_0191 + IT_0234 + IT_0236) + IT_0121*(IT_0201 + IT_0235
       + IT_0237) + IT_0202*(IT_0185 + IT_0236 + IT_0238 + IT_0239) + IT_0152*
      (IT_0196 + IT_0237 + IT_0240 + IT_0241) + IT_0000*IT_0346 + IT_0347
      *IT_0372 + IT_0373*IT_0398 + IT_0399*IT_0426;
    return create_ccomplex_return(IT_0427);
}

