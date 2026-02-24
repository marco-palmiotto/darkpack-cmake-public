#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_h_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = m_W*e_em;
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -4);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = IT_0003*IT_0005*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = IT_0005*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0007, -3);
    const ccomplex_t IT_0013 = IT_0003*IT_0006*IT_0012;
    const ccomplex_t IT_0014 = IT_0009*IT_0010*IT_0012;
    const ccomplex_t IT_0015 = cpow(IT_0004, -2);
    const ccomplex_t IT_0016 = cpow(IT_0007, -1);
    const ccomplex_t IT_0017 = IT_0003*IT_0006*IT_0015*IT_0016;
    const ccomplex_t IT_0018 = IT_0009*IT_0010*IT_0015*IT_0016;
    const ccomplex_t IT_0019 = cpow(IT_0007, -2);
    const ccomplex_t IT_0020 = cpow(IT_0015 + IT_0019, -1);
    const ccomplex_t IT_0021 = (IT_0008 + -IT_0011 + IT_0013 + -IT_0014 + 2
      *IT_0017 + (-2)*IT_0018)*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0004, -1);
    const ccomplex_t IT_0027 = IT_0007*IT_0026;
    const ccomplex_t IT_0028 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0004*IT_0016;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0029 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = pow(m_h, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0038 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0042 = IT_0003*IT_0026;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0045 = IT_0010*IT_0026;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0048 = IT_0003*IT_0016;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0051 = IT_0010*IT_0016;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0043 + IT_0046 + 
      -IT_0049 + -IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0036*IT_0054;
    const ccomplex_t IT_0056 = m_N_1*IT_0055;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0035;
    const ccomplex_t IT_0059 = IT_0054*IT_0058;
    const ccomplex_t IT_0060 = m_N_1*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0056 + IT_0061;
    const ccomplex_t IT_0063 = IT_0040*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0038 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0067 = IT_0042*IT_0066;
    const ccomplex_t IT_0068 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0069 = IT_0045*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0071 = IT_0048*IT_0070;
    const ccomplex_t IT_0072 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + IT_0069 + 
      -IT_0071 + -IT_0073);
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0058*IT_0075;
    const ccomplex_t IT_0077 = m_N_1*IT_0076;
    const ccomplex_t IT_0078 = IT_0036*IT_0075;
    const ccomplex_t IT_0079 = m_N_1*IT_0078;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = IT_0077 + IT_0080;
    const ccomplex_t IT_0082 = IT_0065*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0042*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0045*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0088 = IT_0042*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0090 = IT_0045*IT_0089;
    const ccomplex_t IT_0091 = IT_0084 + IT_0086 + IT_0088 + IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0093 = IT_0048*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0095 = IT_0051*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0097 = IT_0048*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0099 = IT_0051*IT_0098;
    const ccomplex_t IT_0100 = -IT_0093 + -IT_0095 + -IT_0097 + -IT_0099;
    const ccomplex_t IT_0101 = IT_0091 + IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0105 = IT_0027*IT_0104;
    const ccomplex_t IT_0106 = IT_0030*IT_0104;
    const ccomplex_t IT_0107 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0108 = IT_0027*IT_0107;
    const ccomplex_t IT_0109 = IT_0030*IT_0107;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0105 + IT_0106 + 
      -IT_0108 + -IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0103*IT_0111;
    const ccomplex_t IT_0113 = m_N_2*IT_0112;
    const ccomplex_t IT_0114 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0115 = IT_0027*IT_0114;
    const ccomplex_t IT_0116 = IT_0030*IT_0114;
    const ccomplex_t IT_0117 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0118 = IT_0027*IT_0117;
    const ccomplex_t IT_0119 = IT_0030*IT_0117;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0115 + IT_0116 + 
      -IT_0118 + -IT_0119);
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0103*IT_0121;
    const ccomplex_t IT_0123 = m_N_1*IT_0122;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0113 + IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0132 = IT_0042*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0134 = IT_0045*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0136 = IT_0042*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0138 = IT_0045*IT_0137;
    const ccomplex_t IT_0139 = IT_0132 + IT_0134 + IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0141 = IT_0048*IT_0140;
    const ccomplex_t IT_0142 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0143 = IT_0051*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0145 = IT_0048*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0147 = IT_0051*IT_0146;
    const ccomplex_t IT_0148 = -IT_0141 + -IT_0143 + -IT_0145 + -IT_0147;
    const ccomplex_t IT_0149 = IT_0139 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0153 = IT_0027*IT_0152;
    const ccomplex_t IT_0154 = IT_0030*IT_0152;
    const ccomplex_t IT_0155 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0156 = IT_0027*IT_0155;
    const ccomplex_t IT_0157 = IT_0030*IT_0155;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0153 + IT_0154 + 
      -IT_0156 + -IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0151*IT_0159;
    const ccomplex_t IT_0161 = m_N_3*IT_0160;
    const ccomplex_t IT_0162 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0163 = IT_0027*IT_0162;
    const ccomplex_t IT_0164 = IT_0030*IT_0162;
    const ccomplex_t IT_0165 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0166 = IT_0027*IT_0165;
    const ccomplex_t IT_0167 = IT_0030*IT_0165;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0163 + IT_0164 + 
      -IT_0166 + -IT_0167);
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0151*IT_0169;
    const ccomplex_t IT_0171 = m_N_1*IT_0170;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0161 + IT_0172;
    const ccomplex_t IT_0174 = IT_0130*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0178 = IT_0042*IT_0177;
    const ccomplex_t IT_0179 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0180 = IT_0045*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0182 = IT_0042*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0184 = IT_0045*IT_0183;
    const ccomplex_t IT_0185 = IT_0178 + IT_0180 + IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0187 = IT_0048*IT_0186;
    const ccomplex_t IT_0188 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0189 = IT_0051*IT_0188;
    const ccomplex_t IT_0190 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0191 = IT_0048*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0193 = IT_0051*IT_0192;
    const ccomplex_t IT_0194 = -IT_0187 + -IT_0189 + -IT_0191 + -IT_0193;
    const ccomplex_t IT_0195 = IT_0185 + IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0199 = IT_0027*IT_0198;
    const ccomplex_t IT_0200 = IT_0030*IT_0198;
    const ccomplex_t IT_0201 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0202 = IT_0027*IT_0201;
    const ccomplex_t IT_0203 = IT_0030*IT_0201;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0199 + IT_0200 + 
      -IT_0202 + -IT_0203);
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0197*IT_0205;
    const ccomplex_t IT_0207 = m_N_4*IT_0206;
    const ccomplex_t IT_0208 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0209 = IT_0027*IT_0208;
    const ccomplex_t IT_0210 = IT_0030*IT_0208;
    const ccomplex_t IT_0211 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0212 = IT_0027*IT_0211;
    const ccomplex_t IT_0213 = IT_0030*IT_0211;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*(IT_0209 + IT_0210 + 
      -IT_0212 + -IT_0213);
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = IT_0197*IT_0215;
    const ccomplex_t IT_0217 = m_N_1*IT_0216;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = IT_0207 + IT_0218;
    const ccomplex_t IT_0220 = IT_0176*IT_0219;
    const ccomplex_t IT_0221 = IT_0037 + 0.5*IT_0057 + 0.5*IT_0063 + (-0.5)
      *IT_0082 + IT_0128 + IT_0174 + IT_0220;
    const ccomplex_t IT_0222 = IT_0040*IT_0060;
    const ccomplex_t IT_0223 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0226 = IT_0042*IT_0225;
    const ccomplex_t IT_0227 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0228 = IT_0045*IT_0227;
    const ccomplex_t IT_0229 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0230 = IT_0042*IT_0229;
    const ccomplex_t IT_0231 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0232 = IT_0045*IT_0231;
    const ccomplex_t IT_0233 = IT_0226 + IT_0228 + IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0235 = IT_0048*IT_0234;
    const ccomplex_t IT_0236 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0237 = IT_0051*IT_0236;
    const ccomplex_t IT_0238 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0239 = IT_0048*IT_0238;
    const ccomplex_t IT_0240 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0241 = IT_0051*IT_0240;
    const ccomplex_t IT_0242 = -IT_0235 + -IT_0237 + -IT_0239 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0233 + IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0121*IT_0245;
    const ccomplex_t IT_0247 = m_N_2*IT_0246;
    const ccomplex_t IT_0248 = IT_0111*IT_0245;
    const ccomplex_t IT_0249 = m_N_1*IT_0248;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0247 + IT_0250;
    const ccomplex_t IT_0252 = IT_0224*IT_0251;
    const ccomplex_t IT_0253 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0256 = IT_0042*IT_0255;
    const ccomplex_t IT_0257 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0258 = IT_0045*IT_0257;
    const ccomplex_t IT_0259 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0260 = IT_0042*IT_0259;
    const ccomplex_t IT_0261 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0262 = IT_0045*IT_0261;
    const ccomplex_t IT_0263 = IT_0256 + IT_0258 + IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0265 = IT_0048*IT_0264;
    const ccomplex_t IT_0266 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0267 = IT_0051*IT_0266;
    const ccomplex_t IT_0268 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0269 = IT_0048*IT_0268;
    const ccomplex_t IT_0270 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0271 = IT_0051*IT_0270;
    const ccomplex_t IT_0272 = -IT_0265 + -IT_0267 + -IT_0269 + -IT_0271;
    const ccomplex_t IT_0273 = IT_0263 + IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = (-0.5)*IT_0274;
    const ccomplex_t IT_0276 = IT_0169*IT_0275;
    const ccomplex_t IT_0277 = m_N_3*IT_0276;
    const ccomplex_t IT_0278 = IT_0159*IT_0275;
    const ccomplex_t IT_0279 = m_N_1*IT_0278;
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = IT_0277 + IT_0280;
    const ccomplex_t IT_0282 = IT_0254*IT_0281;
    const ccomplex_t IT_0283 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0286 = IT_0042*IT_0285;
    const ccomplex_t IT_0287 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0288 = IT_0045*IT_0287;
    const ccomplex_t IT_0289 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0290 = IT_0042*IT_0289;
    const ccomplex_t IT_0291 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0292 = IT_0045*IT_0291;
    const ccomplex_t IT_0293 = IT_0286 + IT_0288 + IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0295 = IT_0048*IT_0294;
    const ccomplex_t IT_0296 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0297 = IT_0051*IT_0296;
    const ccomplex_t IT_0298 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0299 = IT_0048*IT_0298;
    const ccomplex_t IT_0300 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0301 = IT_0051*IT_0300;
    const ccomplex_t IT_0302 = -IT_0295 + -IT_0297 + -IT_0299 + -IT_0301;
    const ccomplex_t IT_0303 = IT_0293 + IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = (-0.5)*IT_0304;
    const ccomplex_t IT_0306 = IT_0205*IT_0305;
    const ccomplex_t IT_0307 = m_N_1*IT_0306;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0215*IT_0305;
    const ccomplex_t IT_0310 = m_N_4*IT_0309;
    const ccomplex_t IT_0311 = IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = IT_0284*IT_0311;
    const ccomplex_t IT_0313 = (-0.5)*IT_0082 + (-0.5)*IT_0222 + -IT_0252 + 
      -IT_0282 + -IT_0312;
    const ccomplex_t IT_0314 = conj(IT_0221) + conj(IT_0313);
    const ccomplex_t IT_0315 = IT_0040*IT_0055;
    const ccomplex_t IT_0316 = IT_0040*IT_0059;
    const ccomplex_t IT_0317 = IT_0009*IT_0026;
    const ccomplex_t IT_0318 = IT_0041*IT_0317;
    const ccomplex_t IT_0319 = IT_0006*IT_0026;
    const ccomplex_t IT_0320 = IT_0044*IT_0319;
    const ccomplex_t IT_0321 = IT_0009*IT_0016;
    const ccomplex_t IT_0322 = IT_0047*IT_0321;
    const ccomplex_t IT_0323 = IT_0006*IT_0016;
    const ccomplex_t IT_0324 = IT_0050*IT_0323;
    const ccomplex_t IT_0325 = -IT_0318 + IT_0320 + IT_0322 + -IT_0324;
    const ccomplex_t IT_0326 = e_em*(IT_0003*IT_0009 + IT_0006*IT_0010)
      *IT_0016*IT_0026;
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = IT_0325*IT_0327;
    const ccomplex_t IT_0329 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*IT_0329;
    const ccomplex_t IT_0331 = IT_0328*IT_0330;
    const ccomplex_t IT_0332 = IT_0044*IT_0317;
    const ccomplex_t IT_0333 = IT_0041*IT_0319;
    const ccomplex_t IT_0334 = IT_0050*IT_0321;
    const ccomplex_t IT_0335 = IT_0047*IT_0323;
    const ccomplex_t IT_0336 = -IT_0332 + -IT_0333 + IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = e_em*(IT_0003*IT_0006 + -IT_0009*IT_0010)
      *IT_0016*IT_0026;
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = IT_0336*IT_0338;
    const ccomplex_t IT_0340 = pow(m_Z, 2);
    const ccomplex_t IT_0341 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0340 + 
      -reg_prop, -1);
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*IT_0341;
    const ccomplex_t IT_0343 = IT_0339*IT_0342;
    const ccomplex_t IT_0344 = IT_0315 + -IT_0316 + 0.5*IT_0331 + 0.5*IT_0343;
    const ccomplex_t IT_0345 = pow(m_Z, -2);
    const ccomplex_t IT_0346 = s_14*s_24*IT_0345;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = s_12 + IT_0347;
    const ccomplex_t IT_0349 = m_N_1*IT_0348;
    const ccomplex_t IT_0350 = 2*IT_0349;
    const ccomplex_t IT_0351 = IT_0344*IT_0350;
    const ccomplex_t IT_0352 = IT_0059*IT_0065;
    const ccomplex_t IT_0353 = (-0.5)*IT_0352;
    const ccomplex_t IT_0354 = -IT_0059;
    const ccomplex_t IT_0355 = IT_0055 + IT_0354;
    const ccomplex_t IT_0356 = IT_0040*IT_0355;
    const ccomplex_t IT_0357 = (-0.5)*IT_0356;
    const ccomplex_t IT_0358 = s_12*s_14;
    const ccomplex_t IT_0359 = s_24*IT_0000;
    const ccomplex_t IT_0360 = -IT_0359;
    const ccomplex_t IT_0361 = IT_0358 + IT_0360;
    const ccomplex_t IT_0362 = 2*IT_0361;
    const ccomplex_t IT_0363 = IT_0344*IT_0362;
    const ccomplex_t IT_0364 = (-0.5)*IT_0331 + (-0.5)*IT_0343 + IT_0352;
    const ccomplex_t IT_0365 = s_12*s_24;
    const ccomplex_t IT_0366 = s_14*IT_0000;
    const ccomplex_t IT_0367 = -IT_0366;
    const ccomplex_t IT_0368 = IT_0365 + IT_0367;
    const ccomplex_t IT_0369 = 2*IT_0368;
    const ccomplex_t IT_0370 = IT_0364*IT_0369;
    const ccomplex_t IT_0371 = IT_0224*IT_0248;
    const ccomplex_t IT_0372 = IT_0254*IT_0278;
    const ccomplex_t IT_0373 = IT_0284*IT_0306;
    const ccomplex_t IT_0374 = (-0.5)*IT_0331 + (-0.5)*IT_0343 + (-2)*IT_0371 
      + (-2)*IT_0372 + (-2)*IT_0373;
    const ccomplex_t IT_0375 = (-2)*IT_0361;
    const ccomplex_t IT_0376 = IT_0374*IT_0375;
    const ccomplex_t IT_0377 = IT_0371 + IT_0372 + IT_0373;
    const ccomplex_t IT_0378 = IT_0362*IT_0374;
    const ccomplex_t IT_0379 = IT_0055*IT_0065;
    const ccomplex_t IT_0380 = IT_0176*IT_0306;
    const ccomplex_t IT_0381 = IT_0127*IT_0248;
    const ccomplex_t IT_0382 = IT_0130*IT_0278;
    const ccomplex_t IT_0383 = 0.5*IT_0331 + 0.5*IT_0343 + IT_0379 + 2*IT_0380
       + 2*IT_0381 + 2*IT_0382;
    const ccomplex_t IT_0384 = IT_0369*IT_0383;
    const ccomplex_t IT_0385 = IT_0378 + IT_0384;
    const ccomplex_t IT_0386 = IT_0344*IT_0375;
    const ccomplex_t IT_0387 = (-2)*IT_0368;
    const ccomplex_t IT_0388 = IT_0364*IT_0387;
    const ccomplex_t IT_0389 = IT_0065*IT_0078;
    const ccomplex_t IT_0390 = (-0.5)*IT_0389;
    const ccomplex_t IT_0391 = -IT_0078;
    const ccomplex_t IT_0392 = IT_0076 + IT_0391;
    const ccomplex_t IT_0393 = IT_0040*IT_0392;
    const ccomplex_t IT_0394 = (-0.5)*IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0390) + conj(IT_0394);
    const ccomplex_t IT_0396 = IT_0066*IT_0317;
    const ccomplex_t IT_0397 = IT_0068*IT_0319;
    const ccomplex_t IT_0398 = IT_0070*IT_0321;
    const ccomplex_t IT_0399 = IT_0072*IT_0323;
    const ccomplex_t IT_0400 = IT_0396 + -IT_0397 + -IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = IT_0327*IT_0400;
    const ccomplex_t IT_0402 = IT_0330*IT_0401;
    const ccomplex_t IT_0403 = IT_0068*IT_0317;
    const ccomplex_t IT_0404 = IT_0066*IT_0319;
    const ccomplex_t IT_0405 = IT_0072*IT_0321;
    const ccomplex_t IT_0406 = IT_0070*IT_0323;
    const ccomplex_t IT_0407 = IT_0403 + IT_0404 + -IT_0405 + -IT_0406;
    const ccomplex_t IT_0408 = IT_0338*IT_0407;
    const ccomplex_t IT_0409 = IT_0342*IT_0408;
    const ccomplex_t IT_0410 = IT_0389 + (-0.5)*IT_0402 + (-0.5)*IT_0409;
    const ccomplex_t IT_0411 = IT_0369*IT_0410;
    const ccomplex_t IT_0412 = IT_0216*IT_0284;
    const ccomplex_t IT_0413 = IT_0122*IT_0224;
    const ccomplex_t IT_0414 = IT_0170*IT_0254;
    const ccomplex_t IT_0415 = (-0.5)*IT_0402 + (-0.5)*IT_0409 + (-2)*IT_0412 
      + (-2)*IT_0413 + (-2)*IT_0414;
    const ccomplex_t IT_0416 = IT_0375*IT_0415;
    const ccomplex_t IT_0417 = IT_0065*IT_0076;
    const ccomplex_t IT_0418 = IT_0122*IT_0127;
    const ccomplex_t IT_0419 = IT_0130*IT_0170;
    const ccomplex_t IT_0420 = IT_0176*IT_0216;
    const ccomplex_t IT_0421 = 0.5*IT_0402 + 0.5*IT_0409 + IT_0417 + 2*IT_0418
       + 2*IT_0419 + 2*IT_0420;
    const ccomplex_t IT_0422 = IT_0387*IT_0421;
    const ccomplex_t IT_0423 = IT_0412 + IT_0413 + IT_0414;
    const ccomplex_t IT_0424 = IT_0362*IT_0415;
    const ccomplex_t IT_0425 = IT_0369*IT_0421;
    const ccomplex_t IT_0426 = IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = IT_0040*IT_0076;
    const ccomplex_t IT_0428 = IT_0040*IT_0078;
    const ccomplex_t IT_0429 = 0.5*IT_0402 + 0.5*IT_0409 + IT_0427 + -IT_0428;
    const ccomplex_t IT_0430 = IT_0375*IT_0429;
    const ccomplex_t IT_0431 = IT_0387*IT_0410;
    const ccomplex_t IT_0432 = IT_0025*IT_0058;
    const ccomplex_t IT_0433 = IT_0040*IT_0077;
    const ccomplex_t IT_0434 = IT_0062*IT_0065;
    const ccomplex_t IT_0435 = IT_0040*IT_0081;
    const ccomplex_t IT_0436 = IT_0127*IT_0251;
    const ccomplex_t IT_0437 = IT_0130*IT_0281;
    const ccomplex_t IT_0438 = IT_0176*IT_0311;
    const ccomplex_t IT_0439 = IT_0432 + 0.5*IT_0433 + (-0.5)*IT_0434 + 0.5
      *IT_0435 + IT_0436 + IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = IT_0125*IT_0224;
    const ccomplex_t IT_0441 = IT_0173*IT_0254;
    const ccomplex_t IT_0442 = IT_0040*IT_0079;
    const ccomplex_t IT_0443 = IT_0219*IT_0284;
    const ccomplex_t IT_0444 = (-0.5)*IT_0434 + -IT_0440 + -IT_0441 + (-0.5)
      *IT_0442 + -IT_0443;
    const ccomplex_t IT_0445 = conj(IT_0439) + conj(IT_0444);
    const ccomplex_t IT_0446 = IT_0350*IT_0364;
    const ccomplex_t IT_0447 = pow(s_24, 2);
    const ccomplex_t IT_0448 = IT_0345*IT_0447;
    const ccomplex_t IT_0449 = -IT_0448;
    const ccomplex_t IT_0450 = IT_0000 + IT_0449;
    const ccomplex_t IT_0451 = m_N_1*IT_0450;
    const ccomplex_t IT_0452 = -IT_0451;
    const ccomplex_t IT_0453 = 2*IT_0452;
    const ccomplex_t IT_0454 = IT_0410*IT_0453;
    const ccomplex_t IT_0455 = s_12*IT_0340;
    const ccomplex_t IT_0456 = s_14*s_24;
    const ccomplex_t IT_0457 = (-2)*IT_0456;
    const ccomplex_t IT_0458 = IT_0455 + IT_0457;
    const ccomplex_t IT_0459 = IT_0345*IT_0458;
    const ccomplex_t IT_0460 = (-0.5)*IT_0459;
    const ccomplex_t IT_0461 = s_12 + IT_0460;
    const ccomplex_t IT_0462 = 4*IT_0461;
    const ccomplex_t IT_0463 = IT_0439*IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0349;
    const ccomplex_t IT_0465 = IT_0383*IT_0464;
    const ccomplex_t IT_0466 = IT_0415*IT_0464;
    const ccomplex_t IT_0467 = (-2)*IT_0452;
    const ccomplex_t IT_0468 = IT_0421*IT_0467;
    const ccomplex_t IT_0469 = pow(s_14, 2);
    const ccomplex_t IT_0470 = IT_0345*IT_0469;
    const ccomplex_t IT_0471 = -IT_0000;
    const ccomplex_t IT_0472 = IT_0470 + IT_0471;
    const ccomplex_t IT_0473 = s_12*IT_0472;
    const ccomplex_t IT_0474 = 2*IT_0473;
    const ccomplex_t IT_0475 = s_12*IT_0348;
    const ccomplex_t IT_0476 = 2*IT_0475;
    const ccomplex_t IT_0477 = IT_0000*IT_0472;
    const ccomplex_t IT_0478 = 2*IT_0477;
    const ccomplex_t IT_0479 = IT_0000*IT_0348;
    const ccomplex_t IT_0480 = 2*IT_0479;
    const ccomplex_t IT_0481 = m_N_1*IT_0472;
    const ccomplex_t IT_0482 = 2*IT_0481;
    const ccomplex_t IT_0483 = (-2)*IT_0473;
    const ccomplex_t IT_0484 = (-2)*IT_0475;
    const ccomplex_t IT_0485 = (-2)*IT_0477;
    const ccomplex_t IT_0486 = (-2)*IT_0479;
    const ccomplex_t IT_0487 = IT_0362*IT_0395 + IT_0350*IT_0445 + conj
      (IT_0429)*IT_0474 + conj(IT_0410)*IT_0476 + conj(IT_0374)*IT_0478 + conj
      (IT_0383)*IT_0480 + IT_0314*IT_0482 + conj(IT_0415)*IT_0483 + conj(IT_0421
      )*IT_0484 + conj(IT_0344)*IT_0485 + conj(IT_0364)*IT_0486;
    const ccomplex_t IT_0488 = (-2)*IT_0481;
    const ccomplex_t IT_0489 = s_12*IT_0450;
    const ccomplex_t IT_0490 = -IT_0489;
    const ccomplex_t IT_0491 = 2*IT_0490;
    const ccomplex_t IT_0492 = IT_0000*IT_0450;
    const ccomplex_t IT_0493 = -IT_0492;
    const ccomplex_t IT_0494 = 2*IT_0493;
    const ccomplex_t IT_0495 = (-2)*IT_0490;
    const ccomplex_t IT_0496 = (-2)*IT_0493;
    const ccomplex_t IT_0497 = (conj(IT_0353) + conj(IT_0357))*IT_0387 +
       IT_0314*IT_0467 + conj(IT_0374)*IT_0476 + conj(IT_0429)*IT_0480 + conj
      (IT_0344)*IT_0484 + conj(IT_0415)*IT_0486 + conj(IT_0383)*IT_0491 + conj
      (IT_0410)*IT_0494 + conj(IT_0364)*IT_0495 + conj(IT_0421)*IT_0496;
    const ccomplex_t IT_0498 = (-0.5)*IT_0379 + -IT_0381;
    const ccomplex_t IT_0499 = -IT_0380 + -IT_0382;
    const ccomplex_t IT_0500 = conj(IT_0498) + conj(IT_0499);
    const ccomplex_t IT_0501 = s_14*IT_0340;
    const ccomplex_t IT_0502 = IT_0345*IT_0501;
    const ccomplex_t IT_0503 = (-0.25)*IT_0502;
    const ccomplex_t IT_0504 = s_14 + IT_0503;
    const ccomplex_t IT_0505 = m_N_1*IT_0504;
    const ccomplex_t IT_0506 = (-4)*IT_0505;
    const ccomplex_t IT_0507 = 2*IT_0506;
    const ccomplex_t IT_0508 = IT_0444*IT_0507;
    const ccomplex_t IT_0509 = s_24*IT_0340;
    const ccomplex_t IT_0510 = IT_0345*IT_0509;
    const ccomplex_t IT_0511 = 0.5*IT_0510;
    const ccomplex_t IT_0512 = s_24 + IT_0511;
    const ccomplex_t IT_0513 = m_N_1*IT_0512;
    const ccomplex_t IT_0514 = 2*IT_0513;
    const ccomplex_t IT_0515 = (-2)*IT_0514;
    const ccomplex_t IT_0516 = IT_0313*IT_0515;
    const ccomplex_t IT_0517 = IT_0362*conj(IT_0429);
    const ccomplex_t IT_0518 = IT_0369*conj(IT_0410);
    const ccomplex_t IT_0519 = IT_0000*IT_0340;
    const ccomplex_t IT_0520 = 3*IT_0519;
    const ccomplex_t IT_0521 = 2*IT_0520;
    const ccomplex_t IT_0522 = 2*IT_0514;
    const ccomplex_t IT_0523 = (-3)*IT_0519;
    const ccomplex_t IT_0524 = 2*IT_0523;
    const ccomplex_t IT_0525 = (-0.5)*IT_0417 + -IT_0418;
    const ccomplex_t IT_0526 = -IT_0419 + -IT_0420;
    const ccomplex_t IT_0527 = conj(IT_0525) + conj(IT_0526);
    const ccomplex_t IT_0528 = pow(m_Z, 4);
    const ccomplex_t IT_0529 = s_12*IT_0528;
    const ccomplex_t IT_0530 = IT_0345*IT_0529;
    const ccomplex_t IT_0531 = (-0.5)*IT_0530;
    const ccomplex_t IT_0532 = IT_0455 + IT_0457 + IT_0531;
    const ccomplex_t IT_0533 = 4*IT_0532;
    const ccomplex_t IT_0534 = (-0.25)*IT_0530;
    const ccomplex_t IT_0535 = IT_0456 + IT_0534;
    const ccomplex_t IT_0536 = 8*IT_0535;
    const ccomplex_t IT_0537 = (-8)*IT_0535;
    const ccomplex_t IT_0538 = (-4)*IT_0532;
    const ccomplex_t IT_0539 = IT_0375*conj(IT_0415);
    const ccomplex_t IT_0540 = IT_0387*conj(IT_0421);
    const ccomplex_t IT_0541 = (-2)*IT_0520;
    const ccomplex_t IT_0542 = (-2)*IT_0506;
    const ccomplex_t IT_0543 = (-2)*IT_0523;
    const ccomplex_t IT_0544 = IT_0517 + IT_0518 + conj(IT_0357)*IT_0521 +
       IT_0445*IT_0522 + IT_0500*IT_0524 + IT_0527*IT_0533 + conj(IT_0394)
      *IT_0536 + conj(IT_0423)*IT_0537 + conj(IT_0390)*IT_0538 + IT_0539 +
       IT_0540 + conj(IT_0377)*IT_0541 + IT_0314*IT_0542 + conj(IT_0353)*IT_0543;
    const ccomplex_t IT_0545 = conj(IT_0344)*IT_0362;
    const ccomplex_t IT_0546 = conj(IT_0364)*IT_0369;
    const ccomplex_t IT_0547 = conj(IT_0374)*IT_0375;
    const ccomplex_t IT_0548 = conj(IT_0383)*IT_0387;
    const ccomplex_t IT_0549 = conj(IT_0394)*IT_0521 + IT_0314*IT_0522 +
       IT_0524*IT_0527 + IT_0500*IT_0533 + conj(IT_0357)*IT_0536 + conj(IT_0377)
      *IT_0537 + conj(IT_0353)*IT_0538 + conj(IT_0423)*IT_0541 + IT_0445*IT_0542
       + conj(IT_0390)*IT_0543 + IT_0545 + IT_0546 + IT_0547 + IT_0548;
    const ccomplex_t IT_0550 = 0.5*IT_0502;
    const ccomplex_t IT_0551 = s_14 + IT_0550;
    const ccomplex_t IT_0552 = m_N_1*IT_0551;
    const ccomplex_t IT_0553 = 2*IT_0552;
    const ccomplex_t IT_0554 = 2*IT_0553;
    const ccomplex_t IT_0555 = (-0.25)*IT_0510;
    const ccomplex_t IT_0556 = s_24 + IT_0555;
    const ccomplex_t IT_0557 = m_N_1*IT_0556;
    const ccomplex_t IT_0558 = (-4)*IT_0557;
    const ccomplex_t IT_0559 = (-2)*IT_0558;
    const ccomplex_t IT_0560 = IT_0517 + IT_0518 + conj(IT_0353)*IT_0521 +
       conj(IT_0377)*IT_0524 + conj(IT_0423)*IT_0533 + conj(IT_0390)*IT_0536 +
       IT_0527*IT_0537 + conj(IT_0394)*IT_0538 + IT_0539 + IT_0540 + IT_0500
      *IT_0541 + conj(IT_0357)*IT_0543 + IT_0314*IT_0554 + IT_0445*IT_0559;
    const ccomplex_t IT_0561 = conj(IT_0390)*IT_0521 + conj(IT_0423)*IT_0524 +
       conj(IT_0377)*IT_0533 + conj(IT_0353)*IT_0536 + IT_0500*IT_0537 + conj
      (IT_0357)*IT_0538 + conj(IT_0394)*IT_0543 + IT_0545 + IT_0546 + IT_0547 +
       IT_0548 + IT_0445*IT_0554 + IT_0314*IT_0559;
    const ccomplex_t IT_0562 = IT_0362*conj(IT_0415);
    const ccomplex_t IT_0563 = IT_0369*conj(IT_0421);
    const ccomplex_t IT_0564 = 2*IT_0558;
    const ccomplex_t IT_0565 = IT_0375*conj(IT_0429);
    const ccomplex_t IT_0566 = IT_0387*conj(IT_0410);
    const ccomplex_t IT_0567 = (-2)*IT_0553;
    const ccomplex_t IT_0568 = IT_0500*IT_0521 + conj(IT_0357)*IT_0524 + conj
      (IT_0394)*IT_0533 + IT_0527*IT_0536 + conj(IT_0390)*IT_0537 + conj(IT_0423
      )*IT_0538 + conj(IT_0353)*IT_0541 + conj(IT_0377)*IT_0543 + IT_0562 +
       IT_0563 + IT_0445*IT_0564 + IT_0565 + IT_0566 + IT_0314*IT_0567;
    const ccomplex_t IT_0569 = IT_0362*conj(IT_0374);
    const ccomplex_t IT_0570 = IT_0369*conj(IT_0383);
    const ccomplex_t IT_0571 = conj(IT_0344)*IT_0375;
    const ccomplex_t IT_0572 = conj(IT_0364)*IT_0387;
    const ccomplex_t IT_0573 = conj(IT_0394)*IT_0524 + IT_0521*IT_0527 + conj
      (IT_0357)*IT_0533 + IT_0500*IT_0536 + conj(IT_0353)*IT_0537 + conj(IT_0377
      )*IT_0538 + conj(IT_0390)*IT_0541 + conj(IT_0423)*IT_0543 + IT_0314
      *IT_0564 + IT_0445*IT_0567 + IT_0569 + IT_0570 + IT_0571 + IT_0572;
    const ccomplex_t IT_0574 = IT_0313*IT_0507;
    const ccomplex_t IT_0575 = IT_0525*IT_0538;
    const ccomplex_t IT_0576 = IT_0357*IT_0541;
    const ccomplex_t IT_0577 = IT_0439*IT_0515;
    const ccomplex_t IT_0578 = IT_0444*IT_0515;
    const ccomplex_t IT_0579 = IT_0499*IT_0543;
    const ccomplex_t IT_0580 = IT_0314*IT_0464;
    const ccomplex_t IT_0581 = conj(IT_0439)*IT_0507;
    const ccomplex_t IT_0582 = conj(IT_0444)*IT_0507;
    const ccomplex_t IT_0583 = conj(IT_0221)*IT_0515;
    const ccomplex_t IT_0584 = conj(IT_0313)*IT_0515;
    const ccomplex_t IT_0585 = conj(IT_0353)*IT_0533;
    const ccomplex_t IT_0586 = conj(IT_0390)*IT_0524;
    const ccomplex_t IT_0587 = conj(IT_0357)*IT_0537;
    const ccomplex_t IT_0588 = IT_0500*IT_0538;
    const ccomplex_t IT_0589 = conj(IT_0394)*IT_0541;
    const ccomplex_t IT_0590 = conj(IT_0377)*IT_0536;
    const ccomplex_t IT_0591 = IT_0350*conj(IT_0429);
    const ccomplex_t IT_0592 = IT_0350*conj(IT_0364);
    const ccomplex_t IT_0593 = conj(IT_0410)*IT_0453;
    const ccomplex_t IT_0594 = conj(IT_0344)*IT_0482;
    const ccomplex_t IT_0595 = 6*IT_0000;
    const ccomplex_t IT_0596 = conj(IT_0221)*IT_0595;
    const ccomplex_t IT_0597 = conj(IT_0383)*IT_0464;
    const ccomplex_t IT_0598 = conj(IT_0415)*IT_0464;
    const ccomplex_t IT_0599 = conj(IT_0421)*IT_0467;
    const ccomplex_t IT_0600 = conj(IT_0374)*IT_0488;
    const ccomplex_t IT_0601 = conj(IT_0353)*IT_0542;
    const ccomplex_t IT_0602 = conj(IT_0313)*IT_0595;
    const ccomplex_t IT_0603 = conj(IT_0394)*IT_0559;
    const ccomplex_t IT_0604 = conj(IT_0377)*IT_0567;
    const ccomplex_t IT_0605 = conj(IT_0423)*IT_0564;
    const ccomplex_t IT_0606 = IT_0439*(IT_0500*IT_0507 + conj(IT_0390)
      *IT_0522 + conj(IT_0357)*IT_0554 + IT_0591 + IT_0592 + IT_0593 + IT_0594 +
       IT_0596 + IT_0597 + IT_0598 + IT_0599 + IT_0600 + IT_0601 + IT_0602 +
       IT_0603 + IT_0604 + IT_0605);
    const ccomplex_t IT_0607 = IT_0444*(IT_0445*IT_0462 + conj(IT_0390)
      *IT_0522 + conj(IT_0357)*IT_0554 + IT_0591 + IT_0592 + IT_0593 + IT_0594 +
       IT_0596 + IT_0597 + IT_0598 + IT_0599 + IT_0600 + IT_0601 + IT_0602 +
       IT_0603 + IT_0604 + IT_0605);
    const ccomplex_t IT_0608 = conj(IT_0221)*IT_0507;
    const ccomplex_t IT_0609 = conj(IT_0313)*IT_0507;
    const ccomplex_t IT_0610 = conj(IT_0357)*IT_0541;
    const ccomplex_t IT_0611 = conj(IT_0439)*IT_0515;
    const ccomplex_t IT_0612 = conj(IT_0444)*IT_0515;
    const ccomplex_t IT_0613 = conj(IT_0498)*IT_0543;
    const ccomplex_t IT_0614 = conj(IT_0499)*IT_0543;
    const ccomplex_t IT_0615 = conj(IT_0353)*IT_0524;
    const ccomplex_t IT_0616 = conj(IT_0390)*IT_0533;
    const ccomplex_t IT_0617 = conj(IT_0394)*IT_0537;
    const ccomplex_t IT_0618 = conj(IT_0423)*IT_0536;
    const ccomplex_t IT_0619 = conj(IT_0344)*IT_0350;
    const ccomplex_t IT_0620 = conj(IT_0364)*IT_0453;
    const ccomplex_t IT_0621 = conj(IT_0374)*IT_0464;
    const ccomplex_t IT_0622 = conj(IT_0383)*IT_0467;
    const ccomplex_t IT_0623 = conj(IT_0429)*IT_0482;
    const ccomplex_t IT_0624 = conj(IT_0390)*IT_0542;
    const ccomplex_t IT_0625 = conj(IT_0415)*IT_0488;
    const ccomplex_t IT_0626 = conj(IT_0421)*IT_0464;
    const ccomplex_t IT_0627 = IT_0350*conj(IT_0410);
    const ccomplex_t IT_0628 = IT_0445*IT_0595;
    const ccomplex_t IT_0629 = IT_0314*IT_0462;
    const ccomplex_t IT_0630 = conj(IT_0357)*IT_0559;
    const ccomplex_t IT_0631 = conj(IT_0423)*IT_0567;
    const ccomplex_t IT_0632 = IT_0500*IT_0515 + conj(IT_0353)*IT_0522 +
       IT_0507*IT_0527 + conj(IT_0394)*IT_0554 + conj(IT_0377)*IT_0564 + IT_0619
       + IT_0620 + IT_0621 + IT_0622 + IT_0623 + IT_0624 + IT_0625 + IT_0626 +
       IT_0627 + IT_0628 + IT_0629 + IT_0630 + IT_0631;
    const ccomplex_t IT_0633 = IT_0314*IT_0351 + (conj(IT_0353) + conj(IT_0357
      ))*(IT_0363 + IT_0370 + IT_0376) + conj(IT_0377)*(IT_0385 + IT_0386 +
       IT_0388) + IT_0395*(IT_0411 + IT_0416 + IT_0422) + conj(IT_0423)*(IT_0426
       + IT_0430 + IT_0431) + IT_0445*(IT_0446 + IT_0454 + IT_0463 + IT_0465 +
       IT_0466 + IT_0468) + IT_0344*(conj(IT_0344)*IT_0474 + conj(IT_0364)
      *IT_0476 + conj(IT_0415)*IT_0478 + conj(IT_0421)*IT_0480 + IT_0445*IT_0482
       + conj(IT_0374)*IT_0483 + conj(IT_0383)*IT_0484 + conj(IT_0429)*IT_0485 +
       conj(IT_0410)*IT_0486) + IT_0429*IT_0487 + IT_0415*(conj(IT_0415)*IT_0474
       + conj(IT_0421)*IT_0476 + conj(IT_0344)*IT_0478 + conj(IT_0364)*IT_0480 +
       conj(IT_0429)*IT_0483 + conj(IT_0410)*IT_0484 + conj(IT_0374)*IT_0485 +
       conj(IT_0383)*IT_0486 + IT_0314*IT_0488) + IT_0410*(IT_0314*IT_0350 +
       conj(IT_0429)*IT_0476 + conj(IT_0374)*IT_0480 + conj(IT_0415)*IT_0484 +
       conj(IT_0344)*IT_0486 + conj(IT_0410)*IT_0491 + conj(IT_0383)*IT_0494 +
       conj(IT_0421)*IT_0495 + conj(IT_0364)*IT_0496) + IT_0364*(IT_0314*IT_0453
       + conj(IT_0344)*IT_0476 + conj(IT_0415)*IT_0480 + conj(IT_0374)*IT_0484 +
       conj(IT_0429)*IT_0486 + conj(IT_0364)*IT_0491 + conj(IT_0421)*IT_0494 +
       conj(IT_0383)*IT_0495 + conj(IT_0410)*IT_0496) + IT_0383*IT_0497 +
       IT_0500*(IT_0385 + IT_0386 + IT_0388 + IT_0508 + IT_0516) + IT_0390
      *IT_0544 + IT_0353*IT_0549 + IT_0394*IT_0560 + IT_0357*IT_0561 + IT_0423
      *IT_0568 + IT_0377*IT_0573 + IT_0527*(IT_0426 + IT_0430 + IT_0431 +
       IT_0574 + IT_0575 + IT_0576 + IT_0577 + IT_0578 + IT_0579) + IT_0374*
      (conj(IT_0374)*IT_0474 + conj(IT_0383)*IT_0476 + conj(IT_0429)*IT_0478 +
       conj(IT_0410)*IT_0480 + conj(IT_0344)*IT_0483 + conj(IT_0364)*IT_0484 +
       conj(IT_0415)*IT_0485 + conj(IT_0421)*IT_0486 + IT_0445*IT_0488 + IT_0580
      ) + IT_0421*(conj(IT_0415)*IT_0476 + conj(IT_0344)*IT_0480 + conj(IT_0429)
      *IT_0484 + conj(IT_0374)*IT_0486 + conj(IT_0421)*IT_0491 + conj(IT_0364)
      *IT_0494 + conj(IT_0410)*IT_0495 + conj(IT_0383)*IT_0496 + IT_0580) +
       IT_0499*(conj(IT_0423)*IT_0521 + IT_0569 + IT_0570 + IT_0571 + IT_0572 +
       IT_0581 + IT_0582 + IT_0583 + IT_0584 + IT_0585 + IT_0586 + IT_0587 +
       IT_0588 + IT_0589 + IT_0590) + IT_0498*(conj(IT_0423)*IT_0521 + IT_0527
      *IT_0543 + IT_0569 + IT_0570 + IT_0571 + IT_0572 + IT_0581 + IT_0582 +
       IT_0583 + IT_0584 + IT_0585 + IT_0586 + IT_0587 + IT_0588 + IT_0589 +
       IT_0590) + IT_0606 + IT_0607 + IT_0525*(conj(IT_0377)*IT_0521 + IT_0562 +
       IT_0563 + IT_0565 + IT_0566 + IT_0608 + IT_0609 + IT_0610 + IT_0611 +
       IT_0612 + IT_0613 + IT_0614 + IT_0615 + IT_0616 + IT_0617 + IT_0618) +
       IT_0526*(conj(IT_0377)*IT_0521 + IT_0527*IT_0538 + IT_0562 + IT_0563 +
       IT_0565 + IT_0566 + IT_0608 + IT_0609 + IT_0610 + IT_0611 + IT_0612 +
       IT_0613 + IT_0614 + IT_0615 + IT_0616 + IT_0617 + IT_0618) + IT_0313*
      (conj(IT_0353)*IT_0522 + conj(IT_0394)*IT_0554 + conj(IT_0377)*IT_0564 +
       IT_0619 + IT_0620 + IT_0621 + IT_0622 + IT_0623 + IT_0624 + IT_0625 +
       IT_0626 + IT_0627 + IT_0628 + IT_0629 + IT_0630 + IT_0631) + IT_0221
      *IT_0632;
    return create_ccomplex_return(IT_0633);
}

