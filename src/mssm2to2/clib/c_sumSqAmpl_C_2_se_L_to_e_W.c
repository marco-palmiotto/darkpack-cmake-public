#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_se_L_to_e_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_se_L_to_e_W(
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
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_e, 2);
    const ccomplex_t IT_0001 = pow(m_se_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0036 = IT_0007*IT_0035;
    const ccomplex_t IT_0037 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0038 = IT_0007*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-0.5)*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0034*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = 2*IT_0019 + 2*IT_0032 + 2*IT_0045 + 2*IT_0058;
    const ccomplex_t IT_0060 = pow(m_W, -2);
    const ccomplex_t IT_0061 = pow(s_14, 2);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = pow(m_C_2, 2);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = s_13*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0068 = cos(theta_W);
    const ccomplex_t IT_0069 = cpow(IT_0068, -1);
    const ccomplex_t IT_0070 = conj(N_B2)*e_em;
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W2)*e_em;
    const ccomplex_t IT_0074 = IT_0007*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0029*IT_0077;
    const ccomplex_t IT_0079 = IT_0020*IT_0067*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0081 = conj(N_B3)*e_em;
    const ccomplex_t IT_0082 = IT_0069*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W3)*e_em;
    const ccomplex_t IT_0085 = IT_0007*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + IT_0086);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0041*IT_0088;
    const ccomplex_t IT_0090 = IT_0043*IT_0080*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0092 = conj(N_B4)*e_em;
    const ccomplex_t IT_0093 = IT_0069*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_W4)*e_em;
    const ccomplex_t IT_0096 = IT_0007*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0054*IT_0099;
    const ccomplex_t IT_0101 = IT_0056*IT_0091*IT_0100;
    const ccomplex_t IT_0102 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0103 = IT_0007*IT_0102;
    const ccomplex_t IT_0104 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0105 = IT_0007*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0103 + 0.5*IT_0106);
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = IT_0099*IT_0108;
    const ccomplex_t IT_0110 = IT_0056*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = m_C_2*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0114 = conj(N_B1)*e_em;
    const ccomplex_t IT_0115 = IT_0069*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W1)*e_em;
    const ccomplex_t IT_0118 = IT_0007*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + IT_0119);
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0016*IT_0121;
    const ccomplex_t IT_0123 = IT_0002*IT_0113*IT_0122;
    const ccomplex_t IT_0124 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0125 = IT_0007*IT_0124;
    const ccomplex_t IT_0126 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0127 = IT_0007*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + 0.5*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0121*IT_0130;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = m_C_2*IT_0133;
    const ccomplex_t IT_0135 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0136 = IT_0007*IT_0135;
    const ccomplex_t IT_0137 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0138 = IT_0007*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + 0.5*IT_0139);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0077*IT_0141;
    const ccomplex_t IT_0143 = IT_0020*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = m_C_2*IT_0144;
    const ccomplex_t IT_0146 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0147 = IT_0007*IT_0146;
    const ccomplex_t IT_0148 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0149 = IT_0007*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 0.5*IT_0150);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = IT_0088*IT_0152;
    const ccomplex_t IT_0154 = IT_0043*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_C_2*IT_0155;
    const ccomplex_t IT_0157 = -IT_0079 + -IT_0090 + -IT_0101 + IT_0112 + 
      -IT_0123 + IT_0134 + IT_0145 + IT_0156;
    const ccomplex_t IT_0158 = m_e*IT_0065;
    const ccomplex_t IT_0159 = (-2)*IT_0111 + (-2)*IT_0133 + (-2)*IT_0144 + (
      -2)*IT_0155;
    const ccomplex_t IT_0160 = m_e*m_C_2;
    const ccomplex_t IT_0161 = IT_0065*IT_0160;
    const ccomplex_t IT_0162 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = cpow(s_12 + 0.5*IT_0001 + 0.5*IT_0063 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = pow(m_W, 2);
    const ccomplex_t IT_0172 = IT_0160*IT_0171;
    const ccomplex_t IT_0173 = 3*IT_0172;
    const ccomplex_t IT_0174 = s_14*IT_0171;
    const ccomplex_t IT_0175 = IT_0060*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = s_14 + IT_0176;
    const ccomplex_t IT_0178 = m_e*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = s_13*IT_0171;
    const ccomplex_t IT_0181 = pow(m_W, 4);
    const ccomplex_t IT_0182 = s_13*IT_0181;
    const ccomplex_t IT_0183 = IT_0060*IT_0182;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = s_14*s_34;
    const ccomplex_t IT_0186 = (-2)*IT_0185;
    const ccomplex_t IT_0187 = IT_0180 + IT_0184 + IT_0186;
    const ccomplex_t IT_0188 = (-3)*IT_0172;
    const ccomplex_t IT_0189 = IT_0111 + IT_0133 + IT_0144;
    const ccomplex_t IT_0190 = IT_0155 + IT_0189;
    const ccomplex_t IT_0191 = (-0.25)*IT_0183;
    const ccomplex_t IT_0192 = IT_0185 + IT_0191;
    const ccomplex_t IT_0193 = (-8)*conj(IT_0170);
    const ccomplex_t IT_0194 = conj(IT_0155) + conj(IT_0189);
    const ccomplex_t IT_0195 = s_34*IT_0171;
    const ccomplex_t IT_0196 = IT_0060*IT_0195;
    const ccomplex_t IT_0197 = (-0.25)*IT_0196;
    const ccomplex_t IT_0198 = s_34 + IT_0197;
    const ccomplex_t IT_0199 = m_C_2*IT_0198;
    const ccomplex_t IT_0200 = (-4)*IT_0199;
    const ccomplex_t IT_0201 = -IT_0169;
    const ccomplex_t IT_0202 = s_14*IT_0000;
    const ccomplex_t IT_0203 = s_13*s_34;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = (-4)*IT_0189;
    const ccomplex_t IT_0207 = s_13*s_14;
    const ccomplex_t IT_0208 = s_34*IT_0063;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = 2*IT_0159;
    const ccomplex_t IT_0212 = (-2)*IT_0159;
    const ccomplex_t IT_0213 = (-2)*conj(IT_0159);
    const ccomplex_t IT_0214 = s_14*s_34*IT_0060;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = s_13 + IT_0215;
    const ccomplex_t IT_0217 = m_C_2*IT_0216;
    const ccomplex_t IT_0218 = (-2)*IT_0157;
    const ccomplex_t IT_0219 = (-2)*conj(IT_0157);
    const ccomplex_t IT_0220 = IT_0160*IT_0216;
    const ccomplex_t IT_0221 = (-2)*conj(IT_0059);
    const ccomplex_t IT_0222 = s_13*IT_0216;
    const ccomplex_t IT_0223 = pow(s_34, 2);
    const ccomplex_t IT_0224 = IT_0060*IT_0223;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0000 + IT_0225;
    const ccomplex_t IT_0227 = m_C_2*IT_0226;
    const ccomplex_t IT_0228 = s_13*IT_0226;
    const ccomplex_t IT_0229 = (-2)*IT_0201;
    const ccomplex_t IT_0230 = 0.5*IT_0196;
    const ccomplex_t IT_0231 = s_34 + IT_0230;
    const ccomplex_t IT_0232 = m_C_2*IT_0231;
    const ccomplex_t IT_0233 = 2*IT_0232;
    const ccomplex_t IT_0234 = (-2)*conj(IT_0201);
    const ccomplex_t IT_0235 = IT_0180 + IT_0186;
    const ccomplex_t IT_0236 = IT_0060*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0022*IT_0141;
    const ccomplex_t IT_0239 = IT_0020*IT_0067*IT_0238;
    const ccomplex_t IT_0240 = IT_0034*IT_0152;
    const ccomplex_t IT_0241 = IT_0043*IT_0080*IT_0240;
    const ccomplex_t IT_0242 = IT_0047*IT_0108;
    const ccomplex_t IT_0243 = IT_0056*IT_0091*IT_0242;
    const ccomplex_t IT_0244 = m_C_2*IT_0058;
    const ccomplex_t IT_0245 = m_e*IT_0169;
    const ccomplex_t IT_0246 = IT_0009*IT_0130;
    const ccomplex_t IT_0247 = IT_0002*IT_0113*IT_0246;
    const ccomplex_t IT_0248 = m_C_2*IT_0019;
    const ccomplex_t IT_0249 = m_C_2*IT_0032;
    const ccomplex_t IT_0250 = m_C_2*IT_0045;
    const ccomplex_t IT_0251 = IT_0239 + IT_0241 + IT_0243 + -IT_0244 + 0.5
      *IT_0245 + IT_0247 + -IT_0248 + -IT_0249 + -IT_0250;
    const ccomplex_t IT_0252 = (-0.25)*IT_0175;
    const ccomplex_t IT_0253 = s_14 + IT_0252;
    const ccomplex_t IT_0254 = m_e*IT_0253;
    const ccomplex_t IT_0255 = (-4)*IT_0254;
    const ccomplex_t IT_0256 = m_e*IT_0216;
    const ccomplex_t IT_0257 = (-0.166666666666667)*IT_0256;
    const ccomplex_t IT_0258 = (-6)*conj(IT_0251)*(IT_0157*IT_0160 + (
      -0.333333333333333)*IT_0179*IT_0190 + (-0.333333333333333)*IT_0162*IT_0200
       + (-0.166666666666667)*IT_0158*IT_0211 + 0.333333333333333*IT_0059
      *IT_0217 + (-0.333333333333333)*IT_0170*IT_0255 + IT_0229*IT_0257);
    const ccomplex_t IT_0259 = (-0.333333333333333)*IT_0158*conj(IT_0159) +
       conj(IT_0157)*IT_0160 + (-0.333333333333333)*IT_0179*IT_0194 + (
      -0.333333333333333)*conj(IT_0162)*IT_0200 + (-0.166666666666667)*IT_0217
      *IT_0221 + (-0.666666666666667)*(s_13 + IT_0237)*conj(IT_0251) + (
      -0.333333333333333)*conj(IT_0170)*IT_0255 + IT_0234*IT_0257;
    const ccomplex_t IT_0260 = 2*conj(IT_0059)*(IT_0059*IT_0066 + IT_0157
      *IT_0158 + IT_0159*IT_0161) + 2*conj(IT_0162)*(IT_0170*IT_0173 + IT_0157
      *IT_0179 + (-2)*IT_0162*IT_0187 + IT_0188*IT_0190) + IT_0190*IT_0192
      *IT_0193 + 2*IT_0194*(IT_0162*IT_0188 + (-4)*IT_0170*IT_0192 + IT_0157
      *IT_0200 + IT_0201*IT_0205 + (-2)*IT_0187*(IT_0155 + (-0.25)*IT_0206)) + 2
      *IT_0210*(conj(IT_0059)*IT_0162 + IT_0059*conj(IT_0162) + conj(IT_0159)
      *IT_0190 + 0.5*IT_0194*IT_0211 + 0.5*conj(IT_0170)*IT_0212 + 0.5*IT_0170
      *IT_0213) + 2*conj(IT_0159)*(IT_0059*IT_0161 + 0.5*IT_0066*IT_0211 + 0.5
      *IT_0217*IT_0218) + IT_0159*IT_0217*IT_0219 + IT_0201*(IT_0220*IT_0221 +
       IT_0213*IT_0222 + IT_0219*IT_0227) + 2*conj(IT_0201)*(IT_0190*IT_0205 + 
      -IT_0059*IT_0220 + 0.5*IT_0212*IT_0222 + 0.5*IT_0218*IT_0227 + 0.5*IT_0228
      *IT_0229) + 2*conj(IT_0170)*(IT_0162*IT_0173 + (-2)*IT_0170*IT_0187 + 0.5
      *IT_0205*IT_0229 + IT_0157*IT_0233) + IT_0170*IT_0205*IT_0234 + 2*conj
      (IT_0157)*(IT_0059*IT_0158 + IT_0162*IT_0179 + IT_0190*IT_0200 + IT_0170
      *IT_0233 + 2*IT_0157*(s_13 + IT_0237)) + IT_0258 + (-6)*IT_0251*IT_0259;
    return create_ccomplex_return(IT_0260);
}

