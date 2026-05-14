#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_se_L_to_anti_e_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_se_L_to_anti_e_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_e, 2);
    const ccomplex_t IT_0001 = pow(m_se_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B1*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_W1*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = IT_0004*IT_0015;
    const ccomplex_t IT_0017 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0015;
    const ccomplex_t IT_0021 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0019 + IT_0023;
    const ccomplex_t IT_0025 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = (-2)*IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = N_B2*e_em;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_W2*e_em;
    const ccomplex_t IT_0038 = IT_0009*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0043 = IT_0020*IT_0042;
    const ccomplex_t IT_0044 = (-2)*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0046 = IT_0016*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0049 = IT_0020*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0050;
    const ccomplex_t IT_0052 = IT_0044 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0041*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = N_B3*e_em;
    const ccomplex_t IT_0060 = IT_0004*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = N_W3*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0068 = IT_0020*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0071 = IT_0016*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0075 = IT_0020*IT_0074;
    const ccomplex_t IT_0076 = (-2)*IT_0075;
    const ccomplex_t IT_0077 = IT_0073 + IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0066*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = N_B4*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = N_W4*e_em;
    const ccomplex_t IT_0088 = IT_0009*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + IT_0089);
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0093 = IT_0020*IT_0092;
    const ccomplex_t IT_0094 = (-2)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0096 = IT_0016*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0099 = IT_0020*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0033 + -IT_0058 + -IT_0083 + -IT_0108;
    const ccomplex_t IT_0110 = m_e*m_C_1*IT_0001;
    const ccomplex_t IT_0111 = pow(m_W, -1);
    const ccomplex_t IT_0112 = cos(beta);
    const ccomplex_t IT_0113 = cpow(IT_0112, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0009*IT_0112;
    const ccomplex_t IT_0117 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = 2*IT_0118;
    const ccomplex_t IT_0120 = IT_0004*IT_0112;
    const ccomplex_t IT_0121 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0125 = IT_0116*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = IT_0123 + IT_0126;
    const ccomplex_t IT_0128 = IT_0119 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0115*IT_0130;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0137 = IT_0120*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0140 = IT_0116*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = IT_0138 + IT_0141;
    const ccomplex_t IT_0143 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0144 = IT_0116*IT_0143;
    const ccomplex_t IT_0145 = 2*IT_0144;
    const ccomplex_t IT_0146 = IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0135*IT_0148;
    const ccomplex_t IT_0150 = IT_0056*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0155 = IT_0116*IT_0154;
    const ccomplex_t IT_0156 = 2*IT_0155;
    const ccomplex_t IT_0157 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0158 = IT_0120*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0161 = IT_0116*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = IT_0159 + IT_0162;
    const ccomplex_t IT_0164 = IT_0156 + IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0153*IT_0166;
    const ccomplex_t IT_0168 = IT_0081*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0173 = IT_0116*IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0173;
    const ccomplex_t IT_0175 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0176 = IT_0120*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0179 = IT_0116*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = IT_0177 + IT_0180;
    const ccomplex_t IT_0182 = IT_0174 + IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0171*IT_0184;
    const ccomplex_t IT_0186 = IT_0106*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = -IT_0133 + -IT_0151 + -IT_0169 + -IT_0187;
    const ccomplex_t IT_0189 = s_23*m_C_1;
    const ccomplex_t IT_0190 = (-2)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0015*(m_W*IT_0112 + (-0.5)*IT_0000*IT_0111*IT_0113);
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0009;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0192*IT_0194;
    const ccomplex_t IT_0196 = pow(m_C_1, 2);
    const ccomplex_t IT_0197 = cpow((-2)*s_13 + IT_0000 + IT_0196 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0198 = IT_0195*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0201 = IT_0014*IT_0130;
    const ccomplex_t IT_0202 = IT_0002*IT_0200*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0204 = IT_0041*IT_0148;
    const ccomplex_t IT_0205 = IT_0056*IT_0203*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0207 = IT_0066*IT_0166;
    const ccomplex_t IT_0208 = IT_0081*IT_0206*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0210 = IT_0091*IT_0184;
    const ccomplex_t IT_0211 = IT_0106*IT_0209*IT_0210;
    const ccomplex_t IT_0212 = m_e*IT_0133;
    const ccomplex_t IT_0213 = m_e*IT_0151;
    const ccomplex_t IT_0214 = m_e*IT_0169;
    const ccomplex_t IT_0215 = m_e*IT_0187;
    const ccomplex_t IT_0216 = IT_0199 + IT_0202 + IT_0205 + IT_0208 + IT_0211
       + -IT_0212 + -IT_0213 + -IT_0214 + -IT_0215;
    const ccomplex_t IT_0217 = 0.5*conj(IT_0216);
    const ccomplex_t IT_0218 = m_e*s_12;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = IT_0192*IT_0221;
    const ccomplex_t IT_0223 = IT_0197*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0030*IT_0115;
    const ccomplex_t IT_0226 = IT_0002*IT_0200*IT_0225;
    const ccomplex_t IT_0227 = IT_0054*IT_0135;
    const ccomplex_t IT_0228 = IT_0056*IT_0203*IT_0227;
    const ccomplex_t IT_0229 = IT_0079*IT_0153;
    const ccomplex_t IT_0230 = IT_0081*IT_0206*IT_0229;
    const ccomplex_t IT_0231 = IT_0104*IT_0171;
    const ccomplex_t IT_0232 = IT_0106*IT_0209*IT_0231;
    const ccomplex_t IT_0233 = m_e*IT_0033;
    const ccomplex_t IT_0234 = m_e*IT_0058;
    const ccomplex_t IT_0235 = m_e*IT_0083;
    const ccomplex_t IT_0236 = m_e*IT_0108;
    const ccomplex_t IT_0237 = -IT_0224 + IT_0226 + IT_0228 + IT_0230 +
       IT_0232 + -IT_0233 + -IT_0234 + -IT_0235 + -IT_0236;
    const ccomplex_t IT_0238 = 0.5*conj(IT_0237);
    const ccomplex_t IT_0239 = s_12*s_23;
    const ccomplex_t IT_0240 = s_13*IT_0001;
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = 4*IT_0242;
    const ccomplex_t IT_0244 = 0.5*conj(IT_0109);
    const ccomplex_t IT_0245 = 0.5*conj(IT_0188);
    const ccomplex_t IT_0246 = m_e*m_C_1;
    const ccomplex_t IT_0247 = 2*IT_0246;
    const ccomplex_t IT_0248 = 2*IT_0109*(IT_0110*conj(IT_0188) + IT_0190
      *IT_0217 + IT_0219*IT_0238 + IT_0243*IT_0244) + 2*IT_0188*(conj(IT_0109)
      *IT_0110 + IT_0217*IT_0219 + IT_0190*IT_0238 + IT_0243*IT_0245) + 2
      *IT_0237*(s_13*conj(IT_0237) + IT_0219*IT_0244 + IT_0190*IT_0245 + IT_0217
      *IT_0247) + 2*IT_0216*(s_13*conj(IT_0216) + IT_0190*IT_0244 + IT_0219
      *IT_0245 + IT_0238*IT_0247);
    return create_ccomplex_return(IT_0248);
}

