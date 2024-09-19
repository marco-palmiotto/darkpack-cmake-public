#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_1_to_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_1_to_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_Hp, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = IT_0009 + IT_0015;
    const ccomplex_t IT_0017 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = 2*IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0002*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0005*IT_0010;
    const ccomplex_t IT_0026 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = IT_0004*IT_0011;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0030;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0027 + IT_0029 + 
      -IT_0031 + -IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = m_N_2*IT_0034;
    const ccomplex_t IT_0036 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = IT_0028*IT_0036;
    const ccomplex_t IT_0039 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0040 = IT_0025*IT_0039;
    const ccomplex_t IT_0041 = IT_0028*IT_0039;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0037 + IT_0038 + 
      -IT_0040 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = m_N_1*IT_0043;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0035 + IT_0045;
    const ccomplex_t IT_0047 = IT_0024*IT_0046;
    const ccomplex_t IT_0048 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0049 = IT_0012*IT_0048;
    const ccomplex_t IT_0050 = 2*IT_0049;
    const ccomplex_t IT_0051 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0055 = IT_0012*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = IT_0053 + IT_0056;
    const ccomplex_t IT_0058 = IT_0050 + IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0065 = IT_0025*IT_0064;
    const ccomplex_t IT_0066 = IT_0028*IT_0064;
    const ccomplex_t IT_0067 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0068 = IT_0025*IT_0067;
    const ccomplex_t IT_0069 = IT_0028*IT_0067;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0065 + IT_0066 + 
      -IT_0068 + -IT_0069);
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = m_N_1*IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0070;
    const ccomplex_t IT_0074 = m_N_1*IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0072 + IT_0075;
    const ccomplex_t IT_0077 = IT_0063*IT_0076;
    const ccomplex_t IT_0078 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0079 = IT_0012*IT_0078;
    const ccomplex_t IT_0080 = 2*IT_0079;
    const ccomplex_t IT_0081 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0082 = IT_0006*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0085 = IT_0012*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = IT_0083 + IT_0086;
    const ccomplex_t IT_0088 = IT_0080 + IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0095 = IT_0025*IT_0094;
    const ccomplex_t IT_0096 = IT_0028*IT_0094;
    const ccomplex_t IT_0097 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0098 = IT_0025*IT_0097;
    const ccomplex_t IT_0099 = IT_0028*IT_0097;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0095 + IT_0096 + 
      -IT_0098 + -IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = m_N_3*IT_0101;
    const ccomplex_t IT_0103 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0104 = IT_0025*IT_0103;
    const ccomplex_t IT_0105 = IT_0028*IT_0103;
    const ccomplex_t IT_0106 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0107 = IT_0025*IT_0106;
    const ccomplex_t IT_0108 = IT_0028*IT_0106;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0104 + IT_0105 + 
      -IT_0107 + -IT_0108);
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = m_N_1*IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0102 + IT_0112;
    const ccomplex_t IT_0114 = IT_0093*IT_0113;
    const ccomplex_t IT_0115 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0116 = IT_0025*IT_0115;
    const ccomplex_t IT_0117 = IT_0028*IT_0115;
    const ccomplex_t IT_0118 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0119 = IT_0025*IT_0118;
    const ccomplex_t IT_0120 = IT_0028*IT_0118;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0116 + IT_0117 + 
      -IT_0119 + -IT_0120);
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = m_N_4*IT_0122;
    const ccomplex_t IT_0124 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0125 = IT_0025*IT_0124;
    const ccomplex_t IT_0126 = IT_0028*IT_0124;
    const ccomplex_t IT_0127 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0128 = IT_0025*IT_0127;
    const ccomplex_t IT_0129 = IT_0028*IT_0127;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0125 + IT_0126 + 
      -IT_0128 + -IT_0129);
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = m_N_1*IT_0131;
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = IT_0123 + IT_0133;
    const ccomplex_t IT_0135 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0136 = IT_0012*IT_0135;
    const ccomplex_t IT_0137 = 2*IT_0136;
    const ccomplex_t IT_0138 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0139 = IT_0006*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0142 = IT_0012*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = IT_0137 + IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = IT_0134*IT_0150;
    const ccomplex_t IT_0152 = -IT_0047 + (-0.5)*IT_0077 + -IT_0114 + -IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0154 = sin(beta);
    const ccomplex_t IT_0155 = IT_0011*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = IT_0005*IT_0154;
    const ccomplex_t IT_0160 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0161 = IT_0159*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = IT_0158 + IT_0162;
    const ccomplex_t IT_0164 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0165 = IT_0155*IT_0164;
    const ccomplex_t IT_0166 = (-2)*IT_0165;
    const ccomplex_t IT_0167 = IT_0163 + IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0171 = IT_0025*IT_0170;
    const ccomplex_t IT_0172 = IT_0028*IT_0170;
    const ccomplex_t IT_0173 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0174 = IT_0028*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + -IT_0172 + (-2)
      *IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0169*IT_0176;
    const ccomplex_t IT_0178 = pow(m_N_1, 2);
    const ccomplex_t IT_0179 = cpow((-2)*s_23 + IT_0000 + IT_0178 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0180 = IT_0153*IT_0177*IT_0179;
    const ccomplex_t IT_0181 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0182 = IT_0025*IT_0181;
    const ccomplex_t IT_0183 = IT_0028*IT_0181;
    const ccomplex_t IT_0184 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0185 = IT_0028*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + -IT_0183 + (-2)
      *IT_0185);
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = IT_0169*IT_0187;
    const ccomplex_t IT_0189 = IT_0179*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = m_C_1*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0193 = cpow((-2)*s_23 + IT_0000 + IT_0178 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0194 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0195 = IT_0025*IT_0194;
    const ccomplex_t IT_0196 = IT_0028*IT_0194;
    const ccomplex_t IT_0197 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0198 = IT_0028*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + -IT_0196 + (-2)
      *IT_0198);
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0202 = IT_0159*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0205 = IT_0155*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = IT_0203 + IT_0206;
    const ccomplex_t IT_0208 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0209 = IT_0155*IT_0208;
    const ccomplex_t IT_0210 = (-2)*IT_0209;
    const ccomplex_t IT_0211 = IT_0207 + IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = IT_0200*IT_0213;
    const ccomplex_t IT_0215 = IT_0192*IT_0193*IT_0214;
    const ccomplex_t IT_0216 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0217 = IT_0025*IT_0216;
    const ccomplex_t IT_0218 = IT_0028*IT_0216;
    const ccomplex_t IT_0219 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0220 = IT_0028*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0217 + -IT_0218 + (-2)
      *IT_0220);
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0213*IT_0222;
    const ccomplex_t IT_0224 = IT_0193*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = m_C_1*IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0077 + -IT_0180 + -IT_0191 + -IT_0215
       + -IT_0226;
    const ccomplex_t IT_0228 = conj(IT_0152) + conj(IT_0227);
    const ccomplex_t IT_0229 = pow(m_Z, -2);
    const ccomplex_t IT_0230 = pow(m_Z, 2);
    const ccomplex_t IT_0231 = s_12*IT_0230;
    const ccomplex_t IT_0232 = s_14*s_24;
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = IT_0231 + IT_0233;
    const ccomplex_t IT_0235 = IT_0229*IT_0234;
    const ccomplex_t IT_0236 = (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = s_12 + IT_0236;
    const ccomplex_t IT_0238 = 4*IT_0237;
    const ccomplex_t IT_0239 = IT_0152*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*e_em*(IT_0005*IT_0010 + 
      -IT_0004*IT_0011);
    const ccomplex_t IT_0241 = cpow((-2)*s_12 + -IT_0001 + -IT_0178 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0242 = IT_0213*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = IT_0240*IT_0243;
    const ccomplex_t IT_0245 = 2*IT_0190 + 2*IT_0225 + -IT_0244;
    const ccomplex_t IT_0246 = s_14*s_24*IT_0229;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = s_12 + IT_0247;
    const ccomplex_t IT_0249 = m_C_1*IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0249;
    const ccomplex_t IT_0251 = IT_0245*IT_0250;
    const ccomplex_t IT_0252 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0253 = IT_0155*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0256 = IT_0159*IT_0255;
    const ccomplex_t IT_0257 = 1.4142135623731*IT_0256;
    const ccomplex_t IT_0258 = IT_0254 + IT_0257;
    const ccomplex_t IT_0259 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0260 = IT_0155*IT_0259;
    const ccomplex_t IT_0261 = (-2)*IT_0260;
    const ccomplex_t IT_0262 = IT_0258 + IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = 0.5*IT_0263;
    const ccomplex_t IT_0265 = IT_0091*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0101*IT_0266;
    const ccomplex_t IT_0268 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0269 = IT_0155*IT_0268;
    const ccomplex_t IT_0270 = (-2)*IT_0269;
    const ccomplex_t IT_0271 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0272 = IT_0159*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0275 = IT_0155*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = IT_0273 + IT_0276;
    const ccomplex_t IT_0278 = IT_0270 + IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = 0.5*IT_0279;
    const ccomplex_t IT_0281 = IT_0002*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0034*IT_0282;
    const ccomplex_t IT_0284 = IT_0061*IT_0213;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = IT_0071*IT_0285;
    const ccomplex_t IT_0287 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0288 = IT_0155*IT_0287;
    const ccomplex_t IT_0289 = (-2)*IT_0288;
    const ccomplex_t IT_0290 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0291 = IT_0159*IT_0290;
    const ccomplex_t IT_0292 = 1.4142135623731*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0294 = IT_0155*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = IT_0292 + IT_0295;
    const ccomplex_t IT_0297 = IT_0289 + IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = 0.5*IT_0298;
    const ccomplex_t IT_0300 = IT_0148*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = IT_0122*IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0267 + 2*IT_0283 + IT_0286 + 2*IT_0302;
    const ccomplex_t IT_0304 = pow(s_24, 2);
    const ccomplex_t IT_0305 = IT_0229*IT_0304;
    const ccomplex_t IT_0306 = -IT_0178;
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = m_C_1*IT_0307;
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = IT_0303*IT_0309;
    const ccomplex_t IT_0311 = m_C_1*m_N_1;
    const ccomplex_t IT_0312 = IT_0248*IT_0311;
    const ccomplex_t IT_0313 = IT_0060*IT_0241;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = IT_0240*IT_0314;
    const ccomplex_t IT_0316 = IT_0063*IT_0071;
    const ccomplex_t IT_0317 = IT_0315 + -IT_0316;
    const ccomplex_t IT_0318 = pow(s_14, 2);
    const ccomplex_t IT_0319 = IT_0229*IT_0318;
    const ccomplex_t IT_0320 = -IT_0001;
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = s_12*IT_0321;
    const ccomplex_t IT_0323 = s_12*s_14;
    const ccomplex_t IT_0324 = s_24*IT_0001;
    const ccomplex_t IT_0325 = -IT_0324;
    const ccomplex_t IT_0326 = IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = -IT_0267 + -IT_0283 + (-0.5)*IT_0286 + -IT_0302;
    const ccomplex_t IT_0328 = -IT_0190 + -IT_0225;
    const ccomplex_t IT_0329 = conj(IT_0327) + conj(IT_0328);
    const ccomplex_t IT_0330 = s_12*IT_0248;
    const ccomplex_t IT_0331 = 2*IT_0330;
    const ccomplex_t IT_0332 = IT_0073*IT_0285;
    const ccomplex_t IT_0333 = (-0.5)*IT_0332;
    const ccomplex_t IT_0334 = (-2)*IT_0326;
    const ccomplex_t IT_0335 = -IT_0244 + IT_0332;
    const ccomplex_t IT_0336 = (-2)*IT_0330;
    const ccomplex_t IT_0337 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0338 = IT_0012*IT_0337;
    const ccomplex_t IT_0339 = 2*IT_0338;
    const ccomplex_t IT_0340 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0341 = IT_0006*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0344 = IT_0012*IT_0343;
    const ccomplex_t IT_0345 = 1.4142135623731*IT_0344;
    const ccomplex_t IT_0346 = IT_0342 + IT_0345;
    const ccomplex_t IT_0347 = IT_0339 + IT_0346;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0348;
    const ccomplex_t IT_0350 = IT_0176*IT_0349;
    const ccomplex_t IT_0351 = IT_0179*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*IT_0351;
    const ccomplex_t IT_0353 = IT_0060*IT_0200;
    const ccomplex_t IT_0354 = IT_0193*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*IT_0354;
    const ccomplex_t IT_0356 = IT_0315 + 2*IT_0352 + 2*IT_0355;
    const ccomplex_t IT_0357 = IT_0311*IT_0321;
    const ccomplex_t IT_0358 = (-2)*IT_0357;
    const ccomplex_t IT_0359 = IT_0063*IT_0073;
    const ccomplex_t IT_0360 = IT_0131*IT_0150;
    const ccomplex_t IT_0361 = IT_0024*IT_0043;
    const ccomplex_t IT_0362 = IT_0093*IT_0110;
    const ccomplex_t IT_0363 = -IT_0359 + (-2)*IT_0360 + (-2)*IT_0361 + (-2)
      *IT_0362;
    const ccomplex_t IT_0364 = (-2)*IT_0312;
    const ccomplex_t IT_0365 = IT_0076*IT_0285;
    const ccomplex_t IT_0366 = m_N_2*IT_0043;
    const ccomplex_t IT_0367 = m_N_1*IT_0034;
    const ccomplex_t IT_0368 = -IT_0367;
    const ccomplex_t IT_0369 = IT_0366 + IT_0368;
    const ccomplex_t IT_0370 = IT_0282*IT_0369;
    const ccomplex_t IT_0371 = m_N_3*IT_0110;
    const ccomplex_t IT_0372 = m_N_1*IT_0101;
    const ccomplex_t IT_0373 = -IT_0372;
    const ccomplex_t IT_0374 = IT_0371 + IT_0373;
    const ccomplex_t IT_0375 = IT_0266*IT_0374;
    const ccomplex_t IT_0376 = m_N_4*IT_0131;
    const ccomplex_t IT_0377 = m_N_1*IT_0122;
    const ccomplex_t IT_0378 = -IT_0377;
    const ccomplex_t IT_0379 = IT_0376 + IT_0378;
    const ccomplex_t IT_0380 = IT_0301*IT_0379;
    const ccomplex_t IT_0381 = (-0.5)*IT_0365 + IT_0370 + IT_0375 + IT_0380;
    const ccomplex_t IT_0382 = IT_0187*IT_0349;
    const ccomplex_t IT_0383 = IT_0153*IT_0179*IT_0382;
    const ccomplex_t IT_0384 = IT_0060*IT_0222;
    const ccomplex_t IT_0385 = IT_0192*IT_0193*IT_0384;
    const ccomplex_t IT_0386 = m_C_1*IT_0352;
    const ccomplex_t IT_0387 = m_C_1*IT_0355;
    const ccomplex_t IT_0388 = (-0.5)*IT_0365 + -IT_0383 + -IT_0385 + -IT_0386
       + -IT_0387;
    const ccomplex_t IT_0389 = conj(IT_0381) + conj(IT_0388);
    const ccomplex_t IT_0390 = m_N_1*IT_0321;
    const ccomplex_t IT_0391 = (-2)*IT_0390;
    const ccomplex_t IT_0392 = IT_0312*conj(IT_0317) + conj(IT_0245)*IT_0322 +
       IT_0326*IT_0329 + 0.5*conj(IT_0303)*IT_0331 + 0.5*conj(IT_0333)*IT_0334 +
       0.5*conj(IT_0335)*IT_0336 + 0.5*conj(IT_0356)*IT_0358 + 0.5*conj(IT_0363)
      *IT_0364 + 0.5*IT_0389*IT_0391;
    const ccomplex_t IT_0393 = IT_0307*IT_0311;
    const ccomplex_t IT_0394 = s_12*IT_0307;
    const ccomplex_t IT_0395 = s_12*s_24;
    const ccomplex_t IT_0396 = s_14*IT_0178;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = IT_0395 + IT_0397;
    const ccomplex_t IT_0399 = m_N_1*IT_0248;
    const ccomplex_t IT_0400 = (-2)*IT_0398;
    const ccomplex_t IT_0401 = (-2)*IT_0394;
    const ccomplex_t IT_0402 = (-2)*IT_0393;
    const ccomplex_t IT_0403 = IT_0228*IT_0308 + 0.5*conj(IT_0245)*IT_0336 +
       IT_0312*conj(IT_0356) + conj(IT_0363)*IT_0393 + conj(IT_0335)*IT_0394 +
       conj(IT_0333)*IT_0398 + IT_0389*IT_0399 + 0.5*IT_0329*IT_0400 + 0.5*conj
      (IT_0303)*IT_0401 + 0.5*conj(IT_0317)*IT_0402;
    const ccomplex_t IT_0404 = 0.5*IT_0316;
    const ccomplex_t IT_0405 = 0.5*IT_0359 + IT_0360 + IT_0361 + IT_0362;
    const ccomplex_t IT_0406 = -IT_0352;
    const ccomplex_t IT_0407 = -IT_0355;
    const ccomplex_t IT_0408 = conj(IT_0405) + conj(IT_0406) + conj(IT_0407);
    const ccomplex_t IT_0409 = conj(IT_0245)*IT_0312 + 0.5*IT_0336*conj
      (IT_0356) + IT_0308*IT_0389 + conj(IT_0303)*IT_0393 + conj(IT_0317)
      *IT_0394 + IT_0228*IT_0399 + 0.5*conj(IT_0363)*IT_0401 + 0.5*conj(IT_0335)
      *IT_0402 + IT_0398*conj(IT_0404) + 0.5*IT_0400*IT_0408;
    const ccomplex_t IT_0410 = (-2)*IT_0399;
    const ccomplex_t IT_0411 = 0.5*conj(IT_0245)*IT_0331 + 0.5*conj(IT_0356)
      *IT_0364 + conj(IT_0317)*IT_0393 + conj(IT_0303)*IT_0394 + IT_0329*IT_0398
       + 0.5*conj(IT_0333)*IT_0400 + 0.5*conj(IT_0335)*IT_0401 + 0.5*conj
      (IT_0363)*IT_0402 + 0.5*IT_0389*IT_0410;
    const ccomplex_t IT_0412 = conj(IT_0406) + conj(IT_0407);
    const ccomplex_t IT_0413 = pow(m_Z, 4);
    const ccomplex_t IT_0414 = s_12*IT_0413;
    const ccomplex_t IT_0415 = IT_0229*IT_0414;
    const ccomplex_t IT_0416 = (-0.5)*IT_0415;
    const ccomplex_t IT_0417 = IT_0231 + IT_0233 + IT_0416;
    const ccomplex_t IT_0418 = (-4)*IT_0417;
    const ccomplex_t IT_0419 = IT_0406*IT_0418;
    const ccomplex_t IT_0420 = s_14*IT_0230;
    const ccomplex_t IT_0421 = IT_0229*IT_0420;
    const ccomplex_t IT_0422 = 0.5*IT_0421;
    const ccomplex_t IT_0423 = s_14 + IT_0422;
    const ccomplex_t IT_0424 = m_N_1*IT_0423;
    const ccomplex_t IT_0425 = 2*IT_0424;
    const ccomplex_t IT_0426 = (-2)*IT_0425;
    const ccomplex_t IT_0427 = IT_0227*IT_0426;
    const ccomplex_t IT_0428 = (-0.25)*IT_0421;
    const ccomplex_t IT_0429 = s_14 + IT_0428;
    const ccomplex_t IT_0430 = m_N_1*IT_0429;
    const ccomplex_t IT_0431 = (-4)*IT_0430;
    const ccomplex_t IT_0432 = s_24*IT_0230;
    const ccomplex_t IT_0433 = IT_0229*IT_0432;
    const ccomplex_t IT_0434 = 0.5*IT_0433;
    const ccomplex_t IT_0435 = s_24 + IT_0434;
    const ccomplex_t IT_0436 = m_C_1*IT_0435;
    const ccomplex_t IT_0437 = 2*IT_0436;
    const ccomplex_t IT_0438 = (-0.25)*IT_0433;
    const ccomplex_t IT_0439 = s_24 + IT_0438;
    const ccomplex_t IT_0440 = m_C_1*IT_0439;
    const ccomplex_t IT_0441 = (-4)*IT_0440;
    const ccomplex_t IT_0442 = IT_0250*conj(IT_0356);
    const ccomplex_t IT_0443 = IT_0309*conj(IT_0363);
    const ccomplex_t IT_0444 = (-2)*IT_0437;
    const ccomplex_t IT_0445 = (-2)*IT_0431;
    const ccomplex_t IT_0446 = 6*IT_0311;
    const ccomplex_t IT_0447 = 0.5*IT_0228;
    const ccomplex_t IT_0448 = IT_0308*conj(IT_0317) + 0.5*IT_0238*IT_0389 +
       0.5*conj(IT_0245)*IT_0391 + conj(IT_0335)*IT_0399 + 0.5*conj(IT_0303)
      *IT_0410 + 0.5*conj(IT_0328)*IT_0426 + conj(IT_0327)*IT_0431 + conj
      (IT_0404)*IT_0437 + IT_0412*IT_0441 + 0.5*IT_0442 + 0.5*IT_0443 + 0.5*conj
      (IT_0405)*IT_0444 + 0.5*conj(IT_0333)*IT_0445 + IT_0446*IT_0447;
    const ccomplex_t IT_0449 = IT_0312*conj(IT_0335) + 0.5*conj(IT_0317)
      *IT_0336 + IT_0322*conj(IT_0356) + 0.5*conj(IT_0245)*IT_0358 + 0.5*IT_0331
      *conj(IT_0363) + 0.5*conj(IT_0303)*IT_0364 + 0.5*IT_0250*IT_0389 + 0.5
      *IT_0334*conj(IT_0404) + IT_0326*IT_0408 + IT_0391*IT_0447;
    const ccomplex_t IT_0450 = 0.5*IT_0331*conj(IT_0356) + 0.5*conj(IT_0245)
      *IT_0364 + 0.5*IT_0309*IT_0389 + conj(IT_0335)*IT_0393 + conj(IT_0363)
      *IT_0394 + 0.5*conj(IT_0317)*IT_0401 + 0.5*conj(IT_0303)*IT_0402 + 0.5
      *IT_0400*conj(IT_0404) + IT_0398*IT_0408 + IT_0410*IT_0447;
    const ccomplex_t IT_0451 = IT_0230*IT_0311;
    const ccomplex_t IT_0452 = (-3)*IT_0451;
    const ccomplex_t IT_0453 = 4*IT_0417;
    const ccomplex_t IT_0454 = (-0.25)*IT_0415;
    const ccomplex_t IT_0455 = IT_0232 + IT_0454;
    const ccomplex_t IT_0456 = (-8)*IT_0455;
    const ccomplex_t IT_0457 = (-2)*IT_0452;
    const ccomplex_t IT_0458 = 3*IT_0451;
    const ccomplex_t IT_0459 = (-2)*IT_0458;
    const ccomplex_t IT_0460 = 0.5*IT_0334*conj(IT_0356) + conj(IT_0317)
      *IT_0398 + 0.5*conj(IT_0363)*IT_0400 + 0.5*conj(IT_0404)*IT_0418 + IT_0389
      *IT_0437 + IT_0445*IT_0447 + conj(IT_0327)*IT_0452 + 0.5*conj(IT_0405)
      *IT_0453 + 0.5*IT_0412*IT_0456 + 0.5*conj(IT_0333)*IT_0457 + 0.5*conj
      (IT_0328)*IT_0459;
    const ccomplex_t IT_0461 = 0.5*conj(IT_0327);
    const ccomplex_t IT_0462 = 0.5*conj(IT_0245)*IT_0334 + conj(IT_0335)
      *IT_0398 + 0.5*conj(IT_0303)*IT_0400 + 0.5*conj(IT_0333)*IT_0418 + IT_0228
      *IT_0437 + 0.5*IT_0389*IT_0445 + conj(IT_0405)*IT_0452 + 0.5*conj(IT_0328)
      *IT_0456 + 0.5*conj(IT_0404)*IT_0457 + 0.5*IT_0412*IT_0459 + IT_0453
      *IT_0461;
    const ccomplex_t IT_0463 = conj(IT_0381)*IT_0446;
    const ccomplex_t IT_0464 = conj(IT_0245)*IT_0250;
    const ccomplex_t IT_0465 = conj(IT_0303)*IT_0309;
    const ccomplex_t IT_0466 = conj(IT_0356)*IT_0391;
    const ccomplex_t IT_0467 = conj(IT_0363)*IT_0410;
    const ccomplex_t IT_0468 = conj(IT_0333)*IT_0437;
    const ccomplex_t IT_0469 = conj(IT_0405)*IT_0431;
    const ccomplex_t IT_0470 = conj(IT_0328)*IT_0441;
    const ccomplex_t IT_0471 = IT_0308*conj(IT_0335) + conj(IT_0317)*IT_0399 +
       0.5*IT_0412*IT_0426 + 0.5*conj(IT_0404)*IT_0445 + 0.5*conj(IT_0388)
      *IT_0446 + IT_0444*IT_0461 + 0.5*IT_0463 + 0.5*IT_0464 + 0.5*IT_0465 + 0.5
      *IT_0466 + 0.5*IT_0467 + IT_0468 + IT_0469 + IT_0470;
    const ccomplex_t IT_0472 = IT_0308*conj(IT_0335) + conj(IT_0317)*IT_0399 +
       0.5*conj(IT_0404)*IT_0445 + 0.5*conj(IT_0388)*IT_0446 + IT_0238*IT_0447 +
       IT_0444*IT_0461 + 0.5*IT_0463 + 0.5*IT_0464 + 0.5*IT_0465 + 0.5*IT_0466 +
       0.5*IT_0467 + IT_0468 + IT_0469 + IT_0470;
    const ccomplex_t IT_0473 = 8*IT_0455;
    const ccomplex_t IT_0474 = conj(IT_0405)*IT_0473;
    const ccomplex_t IT_0475 = conj(IT_0317)*IT_0400;
    const ccomplex_t IT_0476 = conj(IT_0152)*IT_0426;
    const ccomplex_t IT_0477 = conj(IT_0227)*IT_0426;
    const ccomplex_t IT_0478 = IT_0326*conj(IT_0356);
    const ccomplex_t IT_0479 = conj(IT_0363)*IT_0398;
    const ccomplex_t IT_0480 = IT_0389*IT_0441 + 0.5*conj(IT_0404)*IT_0456 +
       0.5*conj(IT_0328)*IT_0457 + conj(IT_0327)*IT_0458 + 0.5*conj(IT_0333)
      *IT_0459 + 0.5*IT_0474 + 0.5*IT_0475 + 0.5*IT_0476 + 0.5*IT_0477 + IT_0478
       + IT_0479;
    const ccomplex_t IT_0481 = 0.5*IT_0412*IT_0418 + IT_0389*IT_0441 + 0.5
      *conj(IT_0404)*IT_0456 + 0.5*conj(IT_0328)*IT_0457 + conj(IT_0327)*IT_0458
       + 0.5*conj(IT_0333)*IT_0459 + 0.5*IT_0474 + 0.5*IT_0475 + 0.5*IT_0476 +
       0.5*IT_0477 + IT_0478 + IT_0479;
    const ccomplex_t IT_0482 = 0.5*conj(IT_0405)*IT_0418 + IT_0228*IT_0431 +
       0.5*IT_0389*IT_0444 + conj(IT_0333)*IT_0452 + 0.5*conj(IT_0404)*IT_0453 +
       conj(IT_0328)*IT_0458 + IT_0457*IT_0461 + 0.5*IT_0412*IT_0473 + 0.5
      *IT_0475 + IT_0478 + IT_0479;
    const ccomplex_t IT_0483 = conj(IT_0245)*IT_0326;
    const ccomplex_t IT_0484 = conj(IT_0303)*IT_0398;
    const ccomplex_t IT_0485 = 0.5*conj(IT_0335)*IT_0400 + 0.5*conj(IT_0328)
      *IT_0418 + 0.5*IT_0389*IT_0426 + IT_0228*IT_0441 + 0.5*conj(IT_0333)
      *IT_0456 + 0.5*IT_0412*IT_0457 + conj(IT_0405)*IT_0458 + 0.5*conj(IT_0404)
      *IT_0459 + IT_0461*IT_0473 + IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = 0.5*conj(IT_0335)*IT_0400 + IT_0389*IT_0431 +
       IT_0444*IT_0447 + conj(IT_0404)*IT_0452 + 0.5*conj(IT_0333)*IT_0453 + 0.5
      *conj(IT_0405)*IT_0457 + IT_0412*IT_0458 + IT_0418*IT_0461 + 0.5*conj
      (IT_0328)*IT_0473 + IT_0483 + IT_0484;
    const ccomplex_t IT_0487 = IT_0228*(IT_0239 + IT_0251 + IT_0310) + 2
      *IT_0245*IT_0392 + 2*IT_0335*IT_0403 + 2*IT_0317*IT_0409 + 2*IT_0303
      *IT_0411 + IT_0412*(IT_0419 + IT_0427) + 2*(IT_0381 + IT_0388)*IT_0448 + 2
      *IT_0356*IT_0449 + 2*IT_0363*IT_0450 + 2*IT_0404*IT_0460 + 2*IT_0333
      *IT_0462 + 2*IT_0152*IT_0471 + 2*IT_0227*IT_0472 + 2*IT_0406*IT_0480 + 2
      *IT_0407*IT_0481 + 2*IT_0405*IT_0482 + 2*IT_0328*IT_0485 + 2*IT_0327
      *IT_0486;
    return create_ccomplex_return(IT_0487);
}

