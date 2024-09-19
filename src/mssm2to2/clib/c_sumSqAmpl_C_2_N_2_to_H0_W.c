#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_H0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_H0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_W = param->Gamma_W;
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
    const ccomplex_t IT_0000 = m_C_2*m_N_2;
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_24, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(IT_0008*IT_0009 +
       IT_0010*IT_0011)*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0013;
    const ccomplex_t IT_0017 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (-2)*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0027 = IT_0016*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = IT_0025 + IT_0028;
    const ccomplex_t IT_0030 = IT_0019 + IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, 2);
    const ccomplex_t IT_0034 = pow(m_C_2, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + -IT_0004 + IT_0033 + -IT_0034 
      + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0015*IT_0037;
    const ccomplex_t IT_0039 = pow(m_H0, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0034 + IT_0039 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0041 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0042 = IT_0009*IT_0041;
    const ccomplex_t IT_0043 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = IT_0013*(IT_0042 + IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0050 = IT_0013*IT_0049;
    const ccomplex_t IT_0051 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0052 = IT_0013*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0060 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0061 = IT_0013*IT_0060;
    const ccomplex_t IT_0062 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0063 = IT_0013*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0070 = IT_0011*IT_0069;
    const ccomplex_t IT_0071 = IT_0013*(IT_0068 + IT_0070);
    const ccomplex_t IT_0072 = 1.4142135623731*e_em*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0059*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*e_em*(IT_0009*IT_0010 + 
      -IT_0008*IT_0011)*IT_0013;
    const ccomplex_t IT_0079 = IT_0010*IT_0013;
    const ccomplex_t IT_0080 = IT_0017*IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0010*IT_0021;
    const ccomplex_t IT_0083 = IT_0023*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = IT_0026*IT_0079;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = IT_0084 + IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + -IT_0004 + -IT_0034 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0078*IT_0093;
    const ccomplex_t IT_0095 = -IT_0038 + (-2)*IT_0058 + (-2)*IT_0077 + 
      -IT_0094;
    const ccomplex_t IT_0096 = U_d2*V_Wp2;
    const ccomplex_t IT_0097 = IT_0009*IT_0096;
    const ccomplex_t IT_0098 = V_u2*U_Wm2;
    const ccomplex_t IT_0099 = IT_0011*IT_0098;
    const ccomplex_t IT_0100 = IT_0013*(IT_0097 + IT_0099);
    const ccomplex_t IT_0101 = 1.4142135623731*e_em*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0105 = IT_0013*IT_0104;
    const ccomplex_t IT_0106 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0107 = IT_0013*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 0.5*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0103*IT_0110;
    const ccomplex_t IT_0112 = IT_0059*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = 2*IT_0113;
    const ccomplex_t IT_0115 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0116 = IT_0016*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0116;
    const ccomplex_t IT_0118 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0119 = IT_0022*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0122 = IT_0016*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = IT_0120 + IT_0123;
    const ccomplex_t IT_0125 = IT_0117 + IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0091*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0078*IT_0129;
    const ccomplex_t IT_0131 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0132 = IT_0013*IT_0131;
    const ccomplex_t IT_0133 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0134 = IT_0013*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = U_d1*V_Wp2;
    const ccomplex_t IT_0139 = IT_0009*IT_0138;
    const ccomplex_t IT_0140 = V_u2*U_Wm1;
    const ccomplex_t IT_0141 = IT_0011*IT_0140;
    const ccomplex_t IT_0142 = IT_0013*(IT_0139 + IT_0141);
    const ccomplex_t IT_0143 = 1.4142135623731*e_em*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0137*IT_0145;
    const ccomplex_t IT_0147 = IT_0040*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = IT_0079*IT_0115;
    const ccomplex_t IT_0150 = 2*IT_0149;
    const ccomplex_t IT_0151 = IT_0082*IT_0118;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = IT_0079*IT_0121;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = IT_0152 + IT_0154;
    const ccomplex_t IT_0156 = IT_0150 + IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = IT_0035*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0015*IT_0160;
    const ccomplex_t IT_0162 = IT_0130 + (-2)*IT_0148 + IT_0161;
    const ccomplex_t IT_0163 = (-2)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0114 + IT_0164;
    const ccomplex_t IT_0166 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = conj(IT_0114) + IT_0167;
    const ccomplex_t IT_0169 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = s_12 + IT_0170;
    const ccomplex_t IT_0172 = IT_0000*IT_0171;
    const ccomplex_t IT_0173 = IT_0009*IT_0021;
    const ccomplex_t IT_0174 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = IT_0011*IT_0013;
    const ccomplex_t IT_0177 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0178 = IT_0176*IT_0177;
    const ccomplex_t IT_0179 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0180 = IT_0173*IT_0179;
    const ccomplex_t IT_0181 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0182 = IT_0176*IT_0181;
    const ccomplex_t IT_0183 = IT_0175 + IT_0178 + IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = IT_0011*IT_0021;
    const ccomplex_t IT_0185 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0186 = IT_0184*IT_0185;
    const ccomplex_t IT_0187 = IT_0009*IT_0013;
    const ccomplex_t IT_0188 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0191 = IT_0184*IT_0190;
    const ccomplex_t IT_0192 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0193 = IT_0187*IT_0192;
    const ccomplex_t IT_0194 = -IT_0186 + -IT_0189 + -IT_0191 + -IT_0193;
    const ccomplex_t IT_0195 = IT_0183 + IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0199 = IT_0013*IT_0198;
    const ccomplex_t IT_0200 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0201 = IT_0013*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0199 + (-0.5)*IT_0202);
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = cpow((-2)*s_23 + IT_0004 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0206 = IT_0204*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0197*IT_0207;
    const ccomplex_t IT_0209 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0210 = IT_0173*IT_0209;
    const ccomplex_t IT_0211 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0212 = IT_0176*IT_0211;
    const ccomplex_t IT_0213 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0214 = IT_0173*IT_0213;
    const ccomplex_t IT_0215 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0216 = IT_0176*IT_0215;
    const ccomplex_t IT_0217 = IT_0210 + IT_0212 + IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0219 = IT_0184*IT_0218;
    const ccomplex_t IT_0220 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0221 = IT_0187*IT_0220;
    const ccomplex_t IT_0222 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0223 = IT_0184*IT_0222;
    const ccomplex_t IT_0224 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0225 = IT_0187*IT_0224;
    const ccomplex_t IT_0226 = -IT_0219 + -IT_0221 + -IT_0223 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0217 + IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = 0.5*IT_0228;
    const ccomplex_t IT_0230 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0231 = IT_0013*IT_0230;
    const ccomplex_t IT_0232 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0233 = IT_0013*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*(IT_0231 + (-0.5)*IT_0234);
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = cpow((-2)*s_23 + IT_0004 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0238 = IT_0236*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = IT_0229*IT_0239;
    const ccomplex_t IT_0241 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0242 = IT_0173*IT_0241;
    const ccomplex_t IT_0243 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0244 = IT_0176*IT_0243;
    const ccomplex_t IT_0245 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0246 = IT_0173*IT_0245;
    const ccomplex_t IT_0247 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0248 = IT_0176*IT_0247;
    const ccomplex_t IT_0249 = IT_0242 + IT_0244 + IT_0246 + IT_0248;
    const ccomplex_t IT_0250 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0251 = IT_0184*IT_0250;
    const ccomplex_t IT_0252 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0253 = IT_0187*IT_0252;
    const ccomplex_t IT_0254 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0255 = IT_0184*IT_0254;
    const ccomplex_t IT_0256 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0257 = IT_0187*IT_0256;
    const ccomplex_t IT_0258 = -IT_0251 + -IT_0253 + -IT_0255 + -IT_0257;
    const ccomplex_t IT_0259 = IT_0249 + IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0263 = IT_0013*IT_0262;
    const ccomplex_t IT_0264 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0265 = IT_0013*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0263 + (-0.5)*IT_0266);
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = cpow((-2)*s_23 + IT_0004 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0270 = IT_0268*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0261*IT_0271;
    const ccomplex_t IT_0273 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0274 = IT_0173*IT_0273;
    const ccomplex_t IT_0275 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0276 = IT_0184*IT_0275;
    const ccomplex_t IT_0277 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0278 = IT_0187*IT_0277;
    const ccomplex_t IT_0279 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0280 = IT_0176*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*(IT_0274 + -IT_0276 + 
      -IT_0278 + IT_0280);
    const ccomplex_t IT_0282 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0283 = IT_0066*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0281*IT_0284;
    const ccomplex_t IT_0286 = IT_0130 + IT_0161 + (-2)*IT_0208 + (-2)*IT_0240
       + (-2)*IT_0272 + (-2)*IT_0285;
    const ccomplex_t IT_0287 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0288 = IT_0173*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0290 = IT_0176*IT_0289;
    const ccomplex_t IT_0291 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0292 = IT_0173*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0294 = IT_0176*IT_0293;
    const ccomplex_t IT_0295 = IT_0288 + IT_0290 + IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0297 = IT_0184*IT_0296;
    const ccomplex_t IT_0298 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0299 = IT_0187*IT_0298;
    const ccomplex_t IT_0300 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0301 = IT_0184*IT_0300;
    const ccomplex_t IT_0302 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0303 = IT_0187*IT_0302;
    const ccomplex_t IT_0304 = -IT_0297 + -IT_0299 + -IT_0301 + -IT_0303;
    const ccomplex_t IT_0305 = IT_0295 + IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = 0.5*IT_0306;
    const ccomplex_t IT_0308 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0309 = IT_0013*IT_0308;
    const ccomplex_t IT_0310 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0311 = IT_0013*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0309 + 0.5*IT_0312);
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = IT_0269*IT_0314;
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*IT_0315;
    const ccomplex_t IT_0317 = IT_0307*IT_0316;
    const ccomplex_t IT_0318 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0319 = IT_0173*IT_0318;
    const ccomplex_t IT_0320 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0321 = IT_0184*IT_0320;
    const ccomplex_t IT_0322 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0323 = IT_0187*IT_0322;
    const ccomplex_t IT_0324 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0325 = IT_0176*IT_0324;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*(IT_0319 + -IT_0321 + 
      -IT_0323 + IT_0325);
    const ccomplex_t IT_0327 = IT_0110*IT_0282;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = IT_0326*IT_0328;
    const ccomplex_t IT_0330 = 2*IT_0317 + 2*IT_0329;
    const ccomplex_t IT_0331 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0332 = IT_0173*IT_0331;
    const ccomplex_t IT_0333 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0334 = IT_0176*IT_0333;
    const ccomplex_t IT_0335 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0336 = IT_0173*IT_0335;
    const ccomplex_t IT_0337 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0338 = IT_0176*IT_0337;
    const ccomplex_t IT_0339 = IT_0332 + IT_0334 + IT_0336 + IT_0338;
    const ccomplex_t IT_0340 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0341 = IT_0184*IT_0340;
    const ccomplex_t IT_0342 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0343 = IT_0187*IT_0342;
    const ccomplex_t IT_0344 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0345 = IT_0184*IT_0344;
    const ccomplex_t IT_0346 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0347 = IT_0187*IT_0346;
    const ccomplex_t IT_0348 = -IT_0341 + -IT_0343 + -IT_0345 + -IT_0347;
    const ccomplex_t IT_0349 = IT_0339 + IT_0348;
    const ccomplex_t IT_0350 = (0 + _Complex_I*1)*IT_0349;
    const ccomplex_t IT_0351 = 0.5*IT_0350;
    const ccomplex_t IT_0352 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0353 = IT_0013*IT_0352;
    const ccomplex_t IT_0354 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0355 = IT_0013*IT_0354;
    const ccomplex_t IT_0356 = 1.4142135623731*IT_0355;
    const ccomplex_t IT_0357 = (0 + _Complex_I*1)*(IT_0353 + 0.5*IT_0356);
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = IT_0237*IT_0358;
    const ccomplex_t IT_0360 = (0 + _Complex_I*1)*IT_0359;
    const ccomplex_t IT_0361 = IT_0351*IT_0360;
    const ccomplex_t IT_0362 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0363 = IT_0173*IT_0362;
    const ccomplex_t IT_0364 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0365 = IT_0176*IT_0364;
    const ccomplex_t IT_0366 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0367 = IT_0173*IT_0366;
    const ccomplex_t IT_0368 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0369 = IT_0176*IT_0368;
    const ccomplex_t IT_0370 = IT_0363 + IT_0365 + IT_0367 + IT_0369;
    const ccomplex_t IT_0371 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0372 = IT_0184*IT_0371;
    const ccomplex_t IT_0373 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0374 = IT_0187*IT_0373;
    const ccomplex_t IT_0375 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0376 = IT_0184*IT_0375;
    const ccomplex_t IT_0377 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0378 = IT_0187*IT_0377;
    const ccomplex_t IT_0379 = -IT_0372 + -IT_0374 + -IT_0376 + -IT_0378;
    const ccomplex_t IT_0380 = IT_0370 + IT_0379;
    const ccomplex_t IT_0381 = (0 + _Complex_I*1)*IT_0380;
    const ccomplex_t IT_0382 = 0.5*IT_0381;
    const ccomplex_t IT_0383 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0384 = IT_0013*IT_0383;
    const ccomplex_t IT_0385 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0386 = IT_0013*IT_0385;
    const ccomplex_t IT_0387 = 1.4142135623731*IT_0386;
    const ccomplex_t IT_0388 = (0 + _Complex_I*1)*(IT_0384 + 0.5*IT_0387);
    const ccomplex_t IT_0389 = -IT_0388;
    const ccomplex_t IT_0390 = IT_0205*IT_0389;
    const ccomplex_t IT_0391 = (0 + _Complex_I*1)*IT_0390;
    const ccomplex_t IT_0392 = IT_0382*IT_0391;
    const ccomplex_t IT_0393 = -IT_0038 + -IT_0094 + (-2)*IT_0361 + (-2)
      *IT_0392;
    const ccomplex_t IT_0394 = (-2)*IT_0114;
    const ccomplex_t IT_0395 = 0.5*IT_0394;
    const ccomplex_t IT_0396 = IT_0162 + IT_0395;
    const ccomplex_t IT_0397 = (-2)*conj(IT_0114);
    const ccomplex_t IT_0398 = 0.5*IT_0397;
    const ccomplex_t IT_0399 = conj(IT_0162) + IT_0398;
    const ccomplex_t IT_0400 = s_12*s_24;
    const ccomplex_t IT_0401 = s_14*IT_0004;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = -IT_0113;
    const ccomplex_t IT_0405 = IT_0208 + IT_0240 + 0.5*IT_0285;
    const ccomplex_t IT_0406 = IT_0272 + 0.5*IT_0285;
    const ccomplex_t IT_0407 = -IT_0317 + -IT_0329;
    const ccomplex_t IT_0408 = IT_0361 + IT_0392;
    const ccomplex_t IT_0409 = -IT_0408;
    const ccomplex_t IT_0410 = IT_0407 + IT_0409;
    const ccomplex_t IT_0411 = -conj(IT_0408);
    const ccomplex_t IT_0412 = conj(IT_0407) + IT_0411;
    const ccomplex_t IT_0413 = 2*IT_0403*(conj(IT_0148)*IT_0396 + IT_0148
      *IT_0399 + IT_0168*(IT_0404 + IT_0405 + IT_0406) + IT_0165*(conj(IT_0404) 
      + conj(IT_0405) + conj(IT_0406)) + conj(IT_0095)*(IT_0058 + IT_0077 +
       IT_0410) + IT_0095*(conj(IT_0058) + conj(IT_0077) + IT_0412));
    const ccomplex_t IT_0414 = s_24*IT_0033;
    const ccomplex_t IT_0415 = IT_0001*IT_0414;
    const ccomplex_t IT_0416 = 0.5*IT_0415;
    const ccomplex_t IT_0417 = s_24 + IT_0416;
    const ccomplex_t IT_0418 = m_C_2*IT_0417;
    const ccomplex_t IT_0419 = 2*IT_0418;
    const ccomplex_t IT_0420 = IT_0074*IT_0110;
    const ccomplex_t IT_0421 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0422 = IT_0059*IT_0420*IT_0421;
    const ccomplex_t IT_0423 = (0 + _Complex_I*1)*m_W*e_em*(IT_0008*IT_0009 +
       IT_0010*IT_0011)*IT_0013;
    const ccomplex_t IT_0424 = IT_0110*IT_0423;
    const ccomplex_t IT_0425 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0004 + -IT_0034 + -reg_prop, -1);
    const ccomplex_t IT_0426 = IT_0424*IT_0425;
    const ccomplex_t IT_0427 = (0 + _Complex_I*1)*IT_0426;
    const ccomplex_t IT_0428 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0429 = IT_0048*IT_0137;
    const ccomplex_t IT_0430 = IT_0040*IT_0428*IT_0429;
    const ccomplex_t IT_0431 = m_N_1*IT_0261;
    const ccomplex_t IT_0432 = IT_0316*IT_0431;
    const ccomplex_t IT_0433 = m_N_2*IT_0281;
    const ccomplex_t IT_0434 = IT_0328*IT_0433;
    const ccomplex_t IT_0435 = IT_0422 + -IT_0427 + IT_0430 + -IT_0432 + 
      -IT_0434;
    const ccomplex_t IT_0436 = m_C_2*IT_0229;
    const ccomplex_t IT_0437 = IT_0239*IT_0436;
    const ccomplex_t IT_0438 = m_C_2*IT_0261;
    const ccomplex_t IT_0439 = IT_0271*IT_0438;
    const ccomplex_t IT_0440 = m_N_2*IT_0058;
    const ccomplex_t IT_0441 = m_C_2*IT_0281;
    const ccomplex_t IT_0442 = IT_0284*IT_0441;
    const ccomplex_t IT_0443 = m_N_2*IT_0077;
    const ccomplex_t IT_0444 = m_N_3*IT_0197;
    const ccomplex_t IT_0445 = IT_0391*IT_0444;
    const ccomplex_t IT_0446 = m_C_2*IT_0197;
    const ccomplex_t IT_0447 = IT_0207*IT_0446;
    const ccomplex_t IT_0448 = m_N_4*IT_0229;
    const ccomplex_t IT_0449 = IT_0360*IT_0448;
    const ccomplex_t IT_0450 = IT_0437 + IT_0439 + -IT_0440 + IT_0442 + 
      -IT_0443 + -IT_0445 + IT_0447 + -IT_0449;
    const ccomplex_t IT_0451 = IT_0435 + IT_0450;
    const ccomplex_t IT_0452 = conj(IT_0435) + conj(IT_0450);
    const ccomplex_t IT_0453 = m_N_2*IT_0113;
    const ccomplex_t IT_0454 = m_C_2*IT_0307;
    const ccomplex_t IT_0455 = IT_0316*IT_0454;
    const ccomplex_t IT_0456 = m_C_2*IT_0326;
    const ccomplex_t IT_0457 = IT_0328*IT_0456;
    const ccomplex_t IT_0458 = -IT_0453 + IT_0455 + IT_0457;
    const ccomplex_t IT_0459 = IT_0066*IT_0423;
    const ccomplex_t IT_0460 = IT_0425*IT_0459;
    const ccomplex_t IT_0461 = (0 + _Complex_I*1)*IT_0460;
    const ccomplex_t IT_0462 = m_N_1*IT_0307;
    const ccomplex_t IT_0463 = IT_0271*IT_0462;
    const ccomplex_t IT_0464 = IT_0066*IT_0103;
    const ccomplex_t IT_0465 = IT_0059*IT_0421*IT_0464;
    const ccomplex_t IT_0466 = IT_0055*IT_0145;
    const ccomplex_t IT_0467 = IT_0040*IT_0428*IT_0466;
    const ccomplex_t IT_0468 = m_N_2*IT_0148;
    const ccomplex_t IT_0469 = m_N_4*IT_0351;
    const ccomplex_t IT_0470 = IT_0239*IT_0469;
    const ccomplex_t IT_0471 = m_N_2*IT_0326;
    const ccomplex_t IT_0472 = IT_0284*IT_0471;
    const ccomplex_t IT_0473 = m_C_2*IT_0382;
    const ccomplex_t IT_0474 = IT_0391*IT_0473;
    const ccomplex_t IT_0475 = m_N_3*IT_0382;
    const ccomplex_t IT_0476 = IT_0207*IT_0475;
    const ccomplex_t IT_0477 = m_C_2*IT_0351;
    const ccomplex_t IT_0478 = IT_0360*IT_0477;
    const ccomplex_t IT_0479 = -IT_0461 + -IT_0463 + IT_0465 + IT_0467 + 
      -IT_0468 + -IT_0470 + -IT_0472 + IT_0474 + -IT_0476 + IT_0478;
    const ccomplex_t IT_0480 = conj(IT_0458) + conj(IT_0479);
    const ccomplex_t IT_0481 = IT_0058 + IT_0077;
    const ccomplex_t IT_0482 = IT_0458 + IT_0479;
    const ccomplex_t IT_0483 = conj(IT_0058) + conj(IT_0077);
    const ccomplex_t IT_0484 = (-2)*IT_0435;
    const ccomplex_t IT_0485 = (-2)*IT_0450;
    const ccomplex_t IT_0486 = IT_0484 + IT_0485;
    const ccomplex_t IT_0487 = (-2)*conj(IT_0435);
    const ccomplex_t IT_0488 = (-2)*conj(IT_0450);
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = (-0.25)*IT_0415;
    const ccomplex_t IT_0491 = s_24 + IT_0490;
    const ccomplex_t IT_0492 = m_C_2*IT_0491;
    const ccomplex_t IT_0493 = (-4)*IT_0492;
    const ccomplex_t IT_0494 = IT_0405 + IT_0406;
    const ccomplex_t IT_0495 = conj(IT_0405) + conj(IT_0406);
    const ccomplex_t IT_0496 = (-2)*IT_0458;
    const ccomplex_t IT_0497 = (-2)*IT_0479;
    const ccomplex_t IT_0498 = IT_0496 + IT_0497;
    const ccomplex_t IT_0499 = (-2)*conj(IT_0458);
    const ccomplex_t IT_0500 = (-2)*conj(IT_0479);
    const ccomplex_t IT_0501 = IT_0499 + IT_0500;
    const ccomplex_t IT_0502 = IT_0000*IT_0033;
    const ccomplex_t IT_0503 = 3*IT_0502;
    const ccomplex_t IT_0504 = -conj(IT_0148);
    const ccomplex_t IT_0505 = conj(IT_0404) + IT_0504;
    const ccomplex_t IT_0506 = -IT_0148;
    const ccomplex_t IT_0507 = IT_0404 + IT_0506;
    const ccomplex_t IT_0508 = (-2)*IT_0405;
    const ccomplex_t IT_0509 = (-2)*IT_0406;
    const ccomplex_t IT_0510 = IT_0508 + IT_0509;
    const ccomplex_t IT_0511 = (-2)*conj(IT_0405);
    const ccomplex_t IT_0512 = (-2)*conj(IT_0406);
    const ccomplex_t IT_0513 = IT_0511 + IT_0512;
    const ccomplex_t IT_0514 = s_12*IT_0033;
    const ccomplex_t IT_0515 = s_14*s_24;
    const ccomplex_t IT_0516 = (-2)*IT_0515;
    const ccomplex_t IT_0517 = IT_0514 + IT_0516;
    const ccomplex_t IT_0518 = IT_0001*IT_0517;
    const ccomplex_t IT_0519 = (-0.5)*IT_0518;
    const ccomplex_t IT_0520 = s_12 + IT_0519;
    const ccomplex_t IT_0521 = 4*IT_0458;
    const ccomplex_t IT_0522 = 6*conj(IT_0435);
    const ccomplex_t IT_0523 = pow(m_W, 4);
    const ccomplex_t IT_0524 = s_12*IT_0523;
    const ccomplex_t IT_0525 = IT_0001*IT_0524;
    const ccomplex_t IT_0526 = (-0.25)*IT_0525;
    const ccomplex_t IT_0527 = IT_0515 + IT_0526;
    const ccomplex_t IT_0528 = 8*conj(IT_0404);
    const ccomplex_t IT_0529 = 8*IT_0407;
    const ccomplex_t IT_0530 = (-8)*IT_0408;
    const ccomplex_t IT_0531 = 8*conj(IT_0407);
    const ccomplex_t IT_0532 = (-8)*conj(IT_0408);
    const ccomplex_t IT_0533 = (-0.5)*IT_0525;
    const ccomplex_t IT_0534 = IT_0514 + IT_0516 + IT_0533;
    const ccomplex_t IT_0535 = (-4)*IT_0405;
    const ccomplex_t IT_0536 = (-4)*IT_0077;
    const ccomplex_t IT_0537 = pow(s_14, 2);
    const ccomplex_t IT_0538 = IT_0001*IT_0537;
    const ccomplex_t IT_0539 = -IT_0034;
    const ccomplex_t IT_0540 = IT_0538 + IT_0539;
    const ccomplex_t IT_0541 = IT_0000*IT_0540;
    const ccomplex_t IT_0542 = (-2)*conj(IT_0286);
    const ccomplex_t IT_0543 = (-2)*IT_0286;
    const ccomplex_t IT_0544 = 0.5*conj(IT_0393);
    const ccomplex_t IT_0545 = s_12*IT_0540;
    const ccomplex_t IT_0546 = (-2)*conj(IT_0330);
    const ccomplex_t IT_0547 = 0.5*IT_0546;
    const ccomplex_t IT_0548 = conj(IT_0393) + IT_0547;
    const ccomplex_t IT_0549 = (-2)*IT_0330;
    const ccomplex_t IT_0550 = s_12*IT_0006;
    const ccomplex_t IT_0551 = 0.5*conj(IT_0162);
    const ccomplex_t IT_0552 = s_12*IT_0171;
    const ccomplex_t IT_0553 = (-2)*IT_0393;
    const ccomplex_t IT_0554 = 0.5*IT_0553;
    const ccomplex_t IT_0555 = IT_0330 + IT_0554;
    const ccomplex_t IT_0556 = (-2)*conj(IT_0393);
    const ccomplex_t IT_0557 = 0.5*IT_0556;
    const ccomplex_t IT_0558 = conj(IT_0330) + IT_0557;
    const ccomplex_t IT_0559 = m_N_2*IT_0171;
    const ccomplex_t IT_0560 = 0.5*IT_0482;
    const ccomplex_t IT_0561 = m_N_2*IT_0540;
    const ccomplex_t IT_0562 = 2*(conj(IT_0330)*IT_0451 + IT_0330*IT_0452 +
       0.5*IT_0393*IT_0487 + 0.5*IT_0480*IT_0543 + IT_0484*IT_0544 + 0.5*conj
      (IT_0450)*IT_0553 + 0.5*IT_0450*IT_0556 + IT_0542*IT_0560)*IT_0561;
    const ccomplex_t IT_0563 = m_C_2*IT_0006;
    const ccomplex_t IT_0564 = 0.5*IT_0451;
    const ccomplex_t IT_0565 = m_C_2*IT_0171;
    const ccomplex_t IT_0566 = s_14*IT_0033;
    const ccomplex_t IT_0567 = IT_0001*IT_0566;
    const ccomplex_t IT_0568 = 0.5*IT_0567;
    const ccomplex_t IT_0569 = s_14 + IT_0568;
    const ccomplex_t IT_0570 = m_N_2*IT_0569;
    const ccomplex_t IT_0571 = 2*IT_0570;
    const ccomplex_t IT_0572 = 0.5*IT_0408;
    const ccomplex_t IT_0573 = 0.5*conj(IT_0408);
    const ccomplex_t IT_0574 = (-3)*IT_0502;
    const ccomplex_t IT_0575 = (-2)*IT_0148;
    const ccomplex_t IT_0576 = (-2)*conj(IT_0148);
    const ccomplex_t IT_0577 = (-0.25)*IT_0567;
    const ccomplex_t IT_0578 = s_14 + IT_0577;
    const ccomplex_t IT_0579 = m_N_2*IT_0578;
    const ccomplex_t IT_0580 = (-4)*IT_0579;
    const ccomplex_t IT_0581 = 0.5*IT_0148;
    const ccomplex_t IT_0582 = 0.5*conj(IT_0148);
    const ccomplex_t IT_0583 = s_12*s_14;
    const ccomplex_t IT_0584 = s_24*IT_0034;
    const ccomplex_t IT_0585 = -IT_0584;
    const ccomplex_t IT_0586 = IT_0583 + IT_0585;
    const ccomplex_t IT_0587 = conj(IT_0286)*(IT_0404 + IT_0405 + IT_0406) +
       IT_0286*(conj(IT_0404) + conj(IT_0405) + conj(IT_0406)) + IT_0408*IT_0548
       + conj(IT_0408)*(IT_0393 + 0.5*IT_0549) + (conj(IT_0058) + conj(IT_0077) 
      + conj(IT_0407))*IT_0555 + (IT_0058 + IT_0077 + IT_0407)*IT_0558 + IT_0542
      *IT_0581 + IT_0543*IT_0582;
    const ccomplex_t IT_0588 = 2*IT_0007*(conj(IT_0095)*IT_0165 + IT_0095
      *IT_0168) + 2*IT_0172*(conj(IT_0095)*IT_0286 + IT_0095*conj(IT_0286) +
       IT_0168*IT_0330 + IT_0165*conj(IT_0330) + conj(IT_0393)*IT_0396 + IT_0393
      *IT_0399) + IT_0413 + 2*IT_0419*(conj(IT_0148)*IT_0451 + IT_0148*IT_0452 +
       IT_0480*IT_0481 + IT_0482*IT_0483 + 0.5*conj(IT_0404)*IT_0486 + 0.5
      *IT_0404*IT_0489) + 2*IT_0493*(IT_0408*IT_0480 + conj(IT_0408)*IT_0482 +
       IT_0452*IT_0494 + IT_0451*IT_0495 + 0.5*conj(IT_0407)*IT_0498 + 0.5
      *IT_0407*IT_0501) + (-2)*IT_0503*(IT_0410*IT_0505 + IT_0412*IT_0507 + (
      -0.5)*IT_0483*IT_0510 + (-0.5)*IT_0481*IT_0513) + 4*IT_0520*(IT_0451
      *IT_0452 + IT_0480*(IT_0479 + 0.25*IT_0521)) + 6*IT_0000*(IT_0451*IT_0480 
      + IT_0482*(conj(IT_0450) + 0.166666666666667*IT_0522)) + 8*IT_0527*
      (IT_0495*IT_0507 + -IT_0494*(conj(IT_0148) + (-0.125)*IT_0528) + 0.125
      *IT_0483*(IT_0529 + IT_0530) + 0.125*IT_0481*(IT_0531 + IT_0532)) + (-4)
      *IT_0534*(IT_0410*IT_0412 + IT_0505*IT_0507 + IT_0495*(IT_0406 + (-0.25)
      *IT_0535) + IT_0483*(IT_0058 + (-0.25)*IT_0536)) + 2*IT_0541*(conj(IT_0286
      )*IT_0330 + IT_0286*conj(IT_0330) + 0.5*IT_0393*IT_0542 + IT_0543*IT_0544)
       + 2*IT_0545*(IT_0286*conj(IT_0286) + IT_0330*conj(IT_0330) + IT_0393
      *IT_0548 + IT_0544*IT_0549) + 2*IT_0550*(IT_0095*conj(IT_0095) + IT_0114
      *conj(IT_0114) + IT_0162*IT_0399 + IT_0394*IT_0551) + 2*IT_0552*(conj
      (IT_0114)*IT_0286 + IT_0114*conj(IT_0286) + 0.5*IT_0162*IT_0542 + IT_0543
      *IT_0551 + conj(IT_0095)*IT_0555 + IT_0095*IT_0558) + 2*IT_0559*(conj
      (IT_0095)*IT_0451 + IT_0095*IT_0452 + IT_0396*IT_0480 + conj(IT_0162)
      *IT_0482 + IT_0397*IT_0560) + IT_0562 + 2*IT_0563*(conj(IT_0162)*IT_0451 +
       IT_0396*IT_0452 + IT_0095*IT_0480 + conj(IT_0095)*IT_0482 + IT_0397
      *IT_0564) + 2*(conj(IT_0330)*IT_0482 + 0.5*IT_0452*IT_0543 + IT_0480
      *IT_0555 + IT_0556*IT_0560 + IT_0542*IT_0564)*IT_0565 + 2*IT_0571*(conj
      (IT_0407)*IT_0451 + IT_0407*IT_0452 + 0.5*IT_0495*IT_0498 + 0.5*IT_0494
      *IT_0501 + IT_0489*IT_0572 + IT_0486*IT_0573) + 2*IT_0574*(conj(IT_0407)
      *IT_0494 + IT_0407*IT_0495 + IT_0513*IT_0572 + IT_0510*IT_0573 + IT_0483*
      (IT_0404 + 0.5*IT_0575) + IT_0481*(conj(IT_0404) + 0.5*IT_0576)) + 2
      *IT_0580*(IT_0404*IT_0480 + conj(IT_0404)*IT_0482 + 0.5*IT_0483*IT_0486 +
       0.5*IT_0481*IT_0489 + IT_0501*IT_0581 + IT_0498*IT_0582) + 2*IT_0586
      *IT_0587;
    return create_ccomplex_return(IT_0588);
}

