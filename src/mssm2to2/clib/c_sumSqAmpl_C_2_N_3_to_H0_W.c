#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_H0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_H0_W(
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = IT_0002*IT_0005;
    const ccomplex_t IT_0008 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 2*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0010 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = pow(m_C_2, 2);
    const ccomplex_t IT_0025 = pow(m_N_3, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + -IT_0024 + -IT_0025 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0006*IT_0028;
    const ccomplex_t IT_0030 = IT_0000*IT_0012;
    const ccomplex_t IT_0031 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = IT_0003*IT_0005;
    const ccomplex_t IT_0034 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0037 = IT_0030*IT_0036;
    const ccomplex_t IT_0038 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = IT_0032 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = IT_0003*IT_0012;
    const ccomplex_t IT_0042 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0000*IT_0005;
    const ccomplex_t IT_0045 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0048 = IT_0041*IT_0047;
    const ccomplex_t IT_0049 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = -IT_0043 + -IT_0046 + -IT_0048 + -IT_0050;
    const ccomplex_t IT_0052 = IT_0040 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0056 = IT_0005*IT_0055;
    const ccomplex_t IT_0057 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0058 = IT_0005*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + (-0.5)*IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = pow(m_H0, 2);
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0025 + IT_0062 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0068 = IT_0030*IT_0067;
    const ccomplex_t IT_0069 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0070 = IT_0041*IT_0069;
    const ccomplex_t IT_0071 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0072 = IT_0044*IT_0071;
    const ccomplex_t IT_0073 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0074 = IT_0033*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + -IT_0070 + 
      -IT_0072 + IT_0074);
    const ccomplex_t IT_0076 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0077 = IT_0005*IT_0076;
    const ccomplex_t IT_0078 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0079 = IT_0005*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + (-0.5)*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0062 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0075*IT_0085;
    const ccomplex_t IT_0087 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0088 = IT_0030*IT_0087;
    const ccomplex_t IT_0089 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0090 = IT_0033*IT_0089;
    const ccomplex_t IT_0091 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0092 = IT_0030*IT_0091;
    const ccomplex_t IT_0093 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0094 = IT_0033*IT_0093;
    const ccomplex_t IT_0095 = IT_0088 + IT_0090 + IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0097 = IT_0041*IT_0096;
    const ccomplex_t IT_0098 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0099 = IT_0044*IT_0098;
    const ccomplex_t IT_0100 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0101 = IT_0041*IT_0100;
    const ccomplex_t IT_0102 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0103 = IT_0044*IT_0102;
    const ccomplex_t IT_0104 = -IT_0097 + -IT_0099 + -IT_0101 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0095 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0109 = IT_0005*IT_0108;
    const ccomplex_t IT_0110 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0111 = IT_0005*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + (-0.5)*IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_23 + IT_0025 + IT_0062 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0107*IT_0117;
    const ccomplex_t IT_0119 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0120 = IT_0030*IT_0119;
    const ccomplex_t IT_0121 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0122 = IT_0033*IT_0121;
    const ccomplex_t IT_0123 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0124 = IT_0030*IT_0123;
    const ccomplex_t IT_0125 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0126 = IT_0033*IT_0125;
    const ccomplex_t IT_0127 = IT_0120 + IT_0122 + IT_0124 + IT_0126;
    const ccomplex_t IT_0128 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0129 = IT_0041*IT_0128;
    const ccomplex_t IT_0130 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0131 = IT_0044*IT_0130;
    const ccomplex_t IT_0132 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0133 = IT_0041*IT_0132;
    const ccomplex_t IT_0134 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0135 = IT_0044*IT_0134;
    const ccomplex_t IT_0136 = -IT_0129 + -IT_0131 + -IT_0133 + -IT_0135;
    const ccomplex_t IT_0137 = IT_0127 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0141 = IT_0005*IT_0140;
    const ccomplex_t IT_0142 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0143 = IT_0005*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + (-0.5)*IT_0144);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_23 + IT_0025 + IT_0062 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0139*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = IT_0001*IT_0005;
    const ccomplex_t IT_0154 = IT_0008*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = IT_0001*IT_0012;
    const ccomplex_t IT_0157 = IT_0014*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = IT_0017*IT_0153;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = IT_0158 + IT_0160;
    const ccomplex_t IT_0162 = IT_0155 + IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = pow(m_W, 2);
    const ccomplex_t IT_0166 = cpow((-2)*s_12 + -IT_0024 + -IT_0025 + IT_0165 
      + -reg_prop, -1);
    const ccomplex_t IT_0167 = IT_0164*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0152*IT_0168;
    const ccomplex_t IT_0170 = IT_0029 + (-2)*IT_0066 + (-2)*IT_0086 + (-2)
      *IT_0118 + (-2)*IT_0150 + IT_0169;
    const ccomplex_t IT_0171 = pow(m_W, -2);
    const ccomplex_t IT_0172 = pow(s_14, 2);
    const ccomplex_t IT_0173 = IT_0171*IT_0172;
    const ccomplex_t IT_0174 = -IT_0024;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = s_12*IT_0175;
    const ccomplex_t IT_0177 = U_d1*V_Wp2;
    const ccomplex_t IT_0178 = IT_0000*IT_0177;
    const ccomplex_t IT_0179 = V_u2*U_Wm1;
    const ccomplex_t IT_0180 = IT_0003*IT_0179;
    const ccomplex_t IT_0181 = IT_0005*(IT_0178 + IT_0180);
    const ccomplex_t IT_0182 = 1.4142135623731*e_em*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0186 = IT_0005*IT_0185;
    const ccomplex_t IT_0187 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0188 = IT_0005*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0186 + 0.5*IT_0189);
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = IT_0184*IT_0191;
    const ccomplex_t IT_0193 = cpow((-2)*s_13 + IT_0024 + IT_0062 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0194 = IT_0192*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = cpow((-2)*s_13 + IT_0062 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0197 = U_d2*V_Wp2;
    const ccomplex_t IT_0198 = IT_0000*IT_0197;
    const ccomplex_t IT_0199 = V_u2*U_Wm2;
    const ccomplex_t IT_0200 = IT_0003*IT_0199;
    const ccomplex_t IT_0201 = IT_0005*(IT_0198 + IT_0200);
    const ccomplex_t IT_0202 = 1.4142135623731*e_em*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0206 = IT_0005*IT_0205;
    const ccomplex_t IT_0207 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0208 = IT_0005*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0206 + 0.5*IT_0209);
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0204*IT_0211;
    const ccomplex_t IT_0213 = IT_0196*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = 2*IT_0195 + 2*IT_0214;
    const ccomplex_t IT_0216 = s_14*s_24*IT_0171;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = s_12 + IT_0217;
    const ccomplex_t IT_0219 = s_12*IT_0218;
    const ccomplex_t IT_0220 = pow(s_24, 2);
    const ccomplex_t IT_0221 = IT_0171*IT_0220;
    const ccomplex_t IT_0222 = -IT_0025;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = s_12*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0226 = IT_0030*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0228 = IT_0033*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0230 = IT_0030*IT_0229;
    const ccomplex_t IT_0231 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0232 = IT_0033*IT_0231;
    const ccomplex_t IT_0233 = IT_0226 + IT_0228 + IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0235 = IT_0041*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0237 = IT_0044*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0239 = IT_0041*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0241 = IT_0044*IT_0240;
    const ccomplex_t IT_0242 = -IT_0235 + -IT_0237 + -IT_0239 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0233 + IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0247 = IT_0005*IT_0246;
    const ccomplex_t IT_0248 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0249 = IT_0005*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*(IT_0247 + 0.5*IT_0250);
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0147*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0245*IT_0254;
    const ccomplex_t IT_0256 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0257 = IT_0030*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0259 = IT_0041*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0261 = IT_0044*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0263 = IT_0033*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0257 + -IT_0259 + 
      -IT_0261 + IT_0263);
    const ccomplex_t IT_0265 = IT_0083*IT_0211;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0264*IT_0266;
    const ccomplex_t IT_0268 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0269 = IT_0030*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0271 = IT_0033*IT_0270;
    const ccomplex_t IT_0272 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0273 = IT_0030*IT_0272;
    const ccomplex_t IT_0274 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0275 = IT_0033*IT_0274;
    const ccomplex_t IT_0276 = IT_0269 + IT_0271 + IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0278 = IT_0041*IT_0277;
    const ccomplex_t IT_0279 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0280 = IT_0044*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0282 = IT_0041*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0284 = IT_0044*IT_0283;
    const ccomplex_t IT_0285 = -IT_0278 + -IT_0280 + -IT_0282 + -IT_0284;
    const ccomplex_t IT_0286 = IT_0276 + IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = 0.5*IT_0287;
    const ccomplex_t IT_0289 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0290 = IT_0005*IT_0289;
    const ccomplex_t IT_0291 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0292 = IT_0005*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*(IT_0290 + 0.5*IT_0293);
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0115*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = IT_0288*IT_0297;
    const ccomplex_t IT_0299 = 2*IT_0255 + 2*IT_0267 + 2*IT_0298;
    const ccomplex_t IT_0300 = m_C_2*m_N_3;
    const ccomplex_t IT_0301 = IT_0175*IT_0300;
    const ccomplex_t IT_0302 = IT_0218*IT_0300;
    const ccomplex_t IT_0303 = IT_0029 + IT_0169;
    const ccomplex_t IT_0304 = (-2)*conj(IT_0303);
    const ccomplex_t IT_0305 = 0.5*IT_0304;
    const ccomplex_t IT_0306 = conj(IT_0215) + IT_0305;
    const ccomplex_t IT_0307 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0308 = IT_0153*IT_0307;
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0311 = IT_0156*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0314 = IT_0153*IT_0313;
    const ccomplex_t IT_0315 = 1.4142135623731*IT_0314;
    const ccomplex_t IT_0316 = IT_0312 + IT_0315;
    const ccomplex_t IT_0317 = IT_0309 + IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = 0.5*IT_0318;
    const ccomplex_t IT_0320 = IT_0026*IT_0319;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = IT_0006*IT_0321;
    const ccomplex_t IT_0323 = IT_0007*IT_0307;
    const ccomplex_t IT_0324 = (-2)*IT_0323;
    const ccomplex_t IT_0325 = IT_0013*IT_0310;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = IT_0007*IT_0313;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = IT_0326 + IT_0328;
    const ccomplex_t IT_0330 = IT_0324 + IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = (-0.5)*IT_0331;
    const ccomplex_t IT_0333 = IT_0166*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*IT_0333;
    const ccomplex_t IT_0335 = IT_0152*IT_0334;
    const ccomplex_t IT_0336 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0337 = IT_0005*IT_0336;
    const ccomplex_t IT_0338 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0339 = IT_0005*IT_0338;
    const ccomplex_t IT_0340 = 1.4142135623731*IT_0339;
    const ccomplex_t IT_0341 = (0 + _Complex_I*1)*(IT_0337 + (-0.5)*IT_0340);
    const ccomplex_t IT_0342 = -IT_0341;
    const ccomplex_t IT_0343 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0344 = IT_0000*IT_0343;
    const ccomplex_t IT_0345 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0346 = IT_0003*IT_0345;
    const ccomplex_t IT_0347 = IT_0005*(IT_0344 + IT_0346);
    const ccomplex_t IT_0348 = 1.4142135623731*e_em*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = 0.5*IT_0349;
    const ccomplex_t IT_0351 = IT_0342*IT_0350;
    const ccomplex_t IT_0352 = IT_0193*IT_0351;
    const ccomplex_t IT_0353 = (0 + _Complex_I*1)*IT_0352;
    const ccomplex_t IT_0354 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0355 = IT_0000*IT_0354;
    const ccomplex_t IT_0356 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0357 = IT_0003*IT_0356;
    const ccomplex_t IT_0358 = IT_0005*(IT_0355 + IT_0357);
    const ccomplex_t IT_0359 = 1.4142135623731*e_em*IT_0358;
    const ccomplex_t IT_0360 = (0 + _Complex_I*1)*IT_0359;
    const ccomplex_t IT_0361 = 0.5*IT_0360;
    const ccomplex_t IT_0362 = IT_0082*IT_0361;
    const ccomplex_t IT_0363 = IT_0196*IT_0362;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*IT_0363;
    const ccomplex_t IT_0365 = -IT_0322 + -IT_0335 + (-2)*IT_0353 + (-2)
      *IT_0364;
    const ccomplex_t IT_0366 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0367 = IT_0030*IT_0366;
    const ccomplex_t IT_0368 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0369 = IT_0033*IT_0368;
    const ccomplex_t IT_0370 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0371 = IT_0030*IT_0370;
    const ccomplex_t IT_0372 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0373 = IT_0033*IT_0372;
    const ccomplex_t IT_0374 = IT_0367 + IT_0369 + IT_0371 + IT_0373;
    const ccomplex_t IT_0375 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0376 = IT_0041*IT_0375;
    const ccomplex_t IT_0377 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0378 = IT_0044*IT_0377;
    const ccomplex_t IT_0379 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0380 = IT_0041*IT_0379;
    const ccomplex_t IT_0381 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0382 = IT_0044*IT_0381;
    const ccomplex_t IT_0383 = -IT_0376 + -IT_0378 + -IT_0380 + -IT_0382;
    const ccomplex_t IT_0384 = IT_0374 + IT_0383;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*IT_0384;
    const ccomplex_t IT_0386 = 0.5*IT_0385;
    const ccomplex_t IT_0387 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0388 = IT_0005*IT_0387;
    const ccomplex_t IT_0389 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0390 = IT_0005*IT_0389;
    const ccomplex_t IT_0391 = 1.4142135623731*IT_0390;
    const ccomplex_t IT_0392 = (0 + _Complex_I*1)*(IT_0388 + 0.5*IT_0391);
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = IT_0063*IT_0393;
    const ccomplex_t IT_0395 = (0 + _Complex_I*1)*IT_0394;
    const ccomplex_t IT_0396 = IT_0386*IT_0395;
    const ccomplex_t IT_0397 = -IT_0322 + -IT_0335 + (-2)*IT_0396;
    const ccomplex_t IT_0398 = (-2)*IT_0397;
    const ccomplex_t IT_0399 = 0.5*IT_0398;
    const ccomplex_t IT_0400 = IT_0299 + IT_0399;
    const ccomplex_t IT_0401 = IT_0223*IT_0300;
    const ccomplex_t IT_0402 = (-2)*IT_0303;
    const ccomplex_t IT_0403 = 0.5*IT_0402;
    const ccomplex_t IT_0404 = IT_0215 + IT_0403;
    const ccomplex_t IT_0405 = m_N_3*IT_0218;
    const ccomplex_t IT_0406 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0407 = IT_0211*IT_0406;
    const ccomplex_t IT_0408 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0024 + -IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0409 = IT_0407*IT_0408;
    const ccomplex_t IT_0410 = (0 + _Complex_I*1)*IT_0409;
    const ccomplex_t IT_0411 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0412 = IT_0191*IT_0350;
    const ccomplex_t IT_0413 = IT_0193*IT_0411*IT_0412;
    const ccomplex_t IT_0414 = IT_0211*IT_0361;
    const ccomplex_t IT_0415 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0416 = IT_0196*IT_0414*IT_0415;
    const ccomplex_t IT_0417 = m_N_1*IT_0107;
    const ccomplex_t IT_0418 = IT_0297*IT_0417;
    const ccomplex_t IT_0419 = m_N_2*IT_0139;
    const ccomplex_t IT_0420 = IT_0254*IT_0419;
    const ccomplex_t IT_0421 = m_N_3*IT_0075;
    const ccomplex_t IT_0422 = IT_0266*IT_0421;
    const ccomplex_t IT_0423 = -IT_0410 + IT_0413 + IT_0416 + -IT_0418 + 
      -IT_0420 + -IT_0422;
    const ccomplex_t IT_0424 = m_N_4*IT_0054;
    const ccomplex_t IT_0425 = IT_0395*IT_0424;
    const ccomplex_t IT_0426 = m_C_2*IT_0075;
    const ccomplex_t IT_0427 = IT_0085*IT_0426;
    const ccomplex_t IT_0428 = m_C_2*IT_0107;
    const ccomplex_t IT_0429 = IT_0117*IT_0428;
    const ccomplex_t IT_0430 = m_C_2*IT_0139;
    const ccomplex_t IT_0431 = IT_0149*IT_0430;
    const ccomplex_t IT_0432 = m_N_3*IT_0353;
    const ccomplex_t IT_0433 = m_N_3*IT_0364;
    const ccomplex_t IT_0434 = m_C_2*IT_0054;
    const ccomplex_t IT_0435 = IT_0065*IT_0434;
    const ccomplex_t IT_0436 = -IT_0425 + IT_0427 + IT_0429 + IT_0431 + 
      -IT_0432 + -IT_0433 + IT_0435;
    const ccomplex_t IT_0437 = IT_0423 + IT_0436;
    const ccomplex_t IT_0438 = m_C_2*IT_0223;
    const ccomplex_t IT_0439 = m_C_2*IT_0288;
    const ccomplex_t IT_0440 = IT_0297*IT_0439;
    const ccomplex_t IT_0441 = m_C_2*IT_0245;
    const ccomplex_t IT_0442 = IT_0254*IT_0441;
    const ccomplex_t IT_0443 = m_C_2*IT_0264;
    const ccomplex_t IT_0444 = IT_0266*IT_0443;
    const ccomplex_t IT_0445 = m_N_3*IT_0214;
    const ccomplex_t IT_0446 = m_N_3*IT_0195;
    const ccomplex_t IT_0447 = IT_0440 + IT_0442 + IT_0444 + -IT_0445 + 
      -IT_0446;
    const ccomplex_t IT_0448 = IT_0082*IT_0204;
    const ccomplex_t IT_0449 = IT_0196*IT_0415*IT_0448;
    const ccomplex_t IT_0450 = IT_0082*IT_0406;
    const ccomplex_t IT_0451 = IT_0408*IT_0450;
    const ccomplex_t IT_0452 = (0 + _Complex_I*1)*IT_0451;
    const ccomplex_t IT_0453 = IT_0184*IT_0342;
    const ccomplex_t IT_0454 = IT_0193*IT_0411*IT_0453;
    const ccomplex_t IT_0455 = m_C_2*IT_0386;
    const ccomplex_t IT_0456 = IT_0395*IT_0455;
    const ccomplex_t IT_0457 = m_N_3*IT_0264;
    const ccomplex_t IT_0458 = IT_0085*IT_0457;
    const ccomplex_t IT_0459 = m_N_1*IT_0288;
    const ccomplex_t IT_0460 = IT_0117*IT_0459;
    const ccomplex_t IT_0461 = m_N_2*IT_0245;
    const ccomplex_t IT_0462 = IT_0149*IT_0461;
    const ccomplex_t IT_0463 = m_N_4*IT_0386;
    const ccomplex_t IT_0464 = IT_0065*IT_0463;
    const ccomplex_t IT_0465 = IT_0449 + -IT_0452 + IT_0454 + IT_0456 + 
      -IT_0458 + -IT_0460 + -IT_0462 + -IT_0464;
    const ccomplex_t IT_0466 = IT_0447 + IT_0465;
    const ccomplex_t IT_0467 = s_12*IT_0165;
    const ccomplex_t IT_0468 = pow(m_W, 4);
    const ccomplex_t IT_0469 = s_12*IT_0468;
    const ccomplex_t IT_0470 = IT_0171*IT_0469;
    const ccomplex_t IT_0471 = (-0.5)*IT_0470;
    const ccomplex_t IT_0472 = s_14*s_24;
    const ccomplex_t IT_0473 = (-2)*IT_0472;
    const ccomplex_t IT_0474 = IT_0467 + IT_0471 + IT_0473;
    const ccomplex_t IT_0475 = IT_0066 + IT_0086 + IT_0118 + IT_0150;
    const ccomplex_t IT_0476 = -IT_0255 + -IT_0267 + -IT_0298;
    const ccomplex_t IT_0477 = -IT_0396;
    const ccomplex_t IT_0478 = IT_0476 + IT_0477;
    const ccomplex_t IT_0479 = -conj(IT_0396);
    const ccomplex_t IT_0480 = conj(IT_0476) + IT_0479;
    const ccomplex_t IT_0481 = s_24*IT_0165;
    const ccomplex_t IT_0482 = IT_0171*IT_0481;
    const ccomplex_t IT_0483 = (-0.25)*IT_0482;
    const ccomplex_t IT_0484 = s_24 + IT_0483;
    const ccomplex_t IT_0485 = m_C_2*IT_0484;
    const ccomplex_t IT_0486 = (-4)*IT_0485;
    const ccomplex_t IT_0487 = (-2)*IT_0447;
    const ccomplex_t IT_0488 = (-2)*IT_0465;
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = (-2)*conj(IT_0447);
    const ccomplex_t IT_0491 = (-2)*conj(IT_0465);
    const ccomplex_t IT_0492 = IT_0490 + IT_0491;
    const ccomplex_t IT_0493 = 6*conj(IT_0423);
    const ccomplex_t IT_0494 = m_N_3*IT_0175;
    const ccomplex_t IT_0495 = (-2)*IT_0215;
    const ccomplex_t IT_0496 = 0.5*IT_0495;
    const ccomplex_t IT_0497 = IT_0303 + IT_0496;
    const ccomplex_t IT_0498 = IT_0467 + IT_0473;
    const ccomplex_t IT_0499 = IT_0171*IT_0498;
    const ccomplex_t IT_0500 = (-0.5)*IT_0499;
    const ccomplex_t IT_0501 = s_12 + IT_0500;
    const ccomplex_t IT_0502 = m_C_2*IT_0218;
    const ccomplex_t IT_0503 = (-2)*IT_0170;
    const ccomplex_t IT_0504 = s_14*IT_0165;
    const ccomplex_t IT_0505 = IT_0171*IT_0504;
    const ccomplex_t IT_0506 = 0.5*IT_0505;
    const ccomplex_t IT_0507 = s_14 + IT_0506;
    const ccomplex_t IT_0508 = m_N_3*IT_0507;
    const ccomplex_t IT_0509 = 2*IT_0508;
    const ccomplex_t IT_0510 = 2*conj(IT_0436)*(IT_0365*IT_0405 + 3*IT_0300
      *IT_0466 + IT_0475*IT_0486 + IT_0400*IT_0494 + IT_0438*IT_0497 + 2*IT_0437
      *IT_0501 + 0.5*IT_0502*IT_0503 + IT_0476*IT_0509);
    const ccomplex_t IT_0511 = IT_0353 + IT_0364;
    const ccomplex_t IT_0512 = 0.5*IT_0482;
    const ccomplex_t IT_0513 = s_24 + IT_0512;
    const ccomplex_t IT_0514 = m_C_2*IT_0513;
    const ccomplex_t IT_0515 = 2*IT_0514;
    const ccomplex_t IT_0516 = -IT_0195 + -IT_0214;
    const ccomplex_t IT_0517 = (-0.25)*IT_0505;
    const ccomplex_t IT_0518 = s_14 + IT_0517;
    const ccomplex_t IT_0519 = m_N_3*IT_0518;
    const ccomplex_t IT_0520 = (-4)*IT_0519;
    const ccomplex_t IT_0521 = 4*IT_0447;
    const ccomplex_t IT_0522 = 3*IT_0300*IT_0437 + IT_0365*IT_0438 + IT_0396
      *IT_0486 + IT_0405*IT_0497 + IT_0400*IT_0502 + 0.5*IT_0494*IT_0503 +
       IT_0511*IT_0515 + IT_0516*IT_0520 + 2*IT_0501*(IT_0465 + 0.25*IT_0521);
    const ccomplex_t IT_0523 = (-2)*conj(IT_0170);
    const ccomplex_t IT_0524 = s_12*s_14;
    const ccomplex_t IT_0525 = s_24*IT_0024;
    const ccomplex_t IT_0526 = -IT_0525;
    const ccomplex_t IT_0527 = IT_0524 + IT_0526;
    const ccomplex_t IT_0528 = (-2)*IT_0299;
    const ccomplex_t IT_0529 = 0.5*IT_0528;
    const ccomplex_t IT_0530 = IT_0397 + IT_0529;
    const ccomplex_t IT_0531 = IT_0475 + IT_0516;
    const ccomplex_t IT_0532 = conj(IT_0475) + conj(IT_0516);
    const ccomplex_t IT_0533 = (-2)*conj(IT_0299);
    const ccomplex_t IT_0534 = IT_0476 + IT_0511;
    const ccomplex_t IT_0535 = (-2)*conj(IT_0397);
    const ccomplex_t IT_0536 = s_12*s_24;
    const ccomplex_t IT_0537 = s_14*IT_0025;
    const ccomplex_t IT_0538 = -IT_0537;
    const ccomplex_t IT_0539 = IT_0536 + IT_0538;
    const ccomplex_t IT_0540 = (-2)*IT_0423;
    const ccomplex_t IT_0541 = (-2)*conj(IT_0423);
    const ccomplex_t IT_0542 = IT_0165*IT_0300;
    const ccomplex_t IT_0543 = (-3)*IT_0542;
    const ccomplex_t IT_0544 = (-0.25)*IT_0470;
    const ccomplex_t IT_0545 = IT_0472 + IT_0544;
    const ccomplex_t IT_0546 = (-2)*IT_0436;
    const ccomplex_t IT_0547 = IT_0540 + IT_0546;
    const ccomplex_t IT_0548 = 3*IT_0542;
    const ccomplex_t IT_0549 = (-2)*IT_0475;
    const ccomplex_t IT_0550 = (-2)*conj(IT_0436);
    const ccomplex_t IT_0551 = IT_0541 + IT_0550;
    const ccomplex_t IT_0552 = (-2)*conj(IT_0475);
    const ccomplex_t IT_0553 = 0.5*IT_0396;
    const ccomplex_t IT_0554 = 0.5*conj(IT_0396);
    const ccomplex_t IT_0555 = 0.5*IT_0219;
    const ccomplex_t IT_0556 = (-0.5)*IT_0515;
    const ccomplex_t IT_0557 = IT_0466*IT_0502;
    const ccomplex_t IT_0558 = (-2)*conj(IT_0215);
    const ccomplex_t IT_0559 = IT_0302*IT_0397;
    const ccomplex_t IT_0560 = IT_0405*IT_0466;
    const ccomplex_t IT_0561 = IT_0437*IT_0438;
    const ccomplex_t IT_0562 = 2*IT_0170*(conj(IT_0170)*IT_0176 + conj(IT_0215
      )*IT_0219) + 2*IT_0215*(conj(IT_0170)*IT_0219 + conj(IT_0215)*IT_0224) + 2
      *IT_0299*(conj(IT_0170)*IT_0301 + IT_0302*IT_0306) + 2*conj(IT_0365)*
      (IT_0170*IT_0302 + IT_0224*IT_0365 + IT_0219*IT_0400 + IT_0401*IT_0404 +
       IT_0405*IT_0437 + IT_0438*IT_0466) + (-4)*IT_0474*(IT_0475*conj(IT_0475) 
      + IT_0478*IT_0480) + 2*IT_0486*(conj(IT_0396)*IT_0466 + conj(IT_0423)
      *IT_0475 + IT_0437*conj(IT_0475) + 0.5*conj(IT_0476)*IT_0489 + 0.5*IT_0476
      *IT_0492) + IT_0300*IT_0466*IT_0493 + 2*conj(IT_0423)*(IT_0365*IT_0405 +
       IT_0299*IT_0494 + IT_0438*IT_0497 + 2*IT_0437*IT_0501 + 0.5*IT_0502
      *IT_0503) + IT_0510 + 2*(conj(IT_0447) + conj(IT_0465))*IT_0522 + (IT_0219
      *IT_0303 + IT_0301*IT_0397 + IT_0466*IT_0494 + IT_0437*IT_0502)*IT_0523 +
       2*IT_0527*(IT_0400*(conj(IT_0476) + conj(IT_0511)) + conj(IT_0396)
      *IT_0530 + conj(IT_0170)*IT_0531 + IT_0170*IT_0532 + IT_0396*(conj(IT_0397
      ) + 0.5*IT_0533) + IT_0534*(conj(IT_0299) + 0.5*IT_0535)) + 2*(conj
      (IT_0365)*(IT_0478 + IT_0511) + IT_0365*(IT_0480 + conj(IT_0511)) + (conj
      (IT_0215) + 0.5*IT_0304)*IT_0531 + IT_0404*IT_0532)*IT_0539 + 2*conj
      (IT_0397)*(IT_0302*IT_0497 + 0.5*IT_0301*IT_0503 + IT_0176*IT_0530 + 0.5
      *IT_0494*IT_0540) + IT_0397*(IT_0176*IT_0533 + IT_0494*IT_0541) + (-4)
      *conj(IT_0511)*(IT_0474*IT_0511 + -1./2*IT_0466*IT_0515 + -1./2*IT_0516
      *IT_0543 + (-2)*IT_0478*IT_0545 + (-0.25)*IT_0520*IT_0547 + (-0.25)
      *IT_0548*IT_0549) + 8*IT_0511*(IT_0480*IT_0545 + 0.125*IT_0520*IT_0551 +
       0.125*IT_0548*IT_0552) + 2*IT_0509*(conj(IT_0423)*IT_0476 + IT_0437*conj
      (IT_0476) + 0.5*conj(IT_0475)*IT_0489 + 0.5*IT_0475*IT_0492 + IT_0551
      *IT_0553 + IT_0547*IT_0554) + 2*IT_0543*(conj(IT_0475)*IT_0476 + IT_0475
      *conj(IT_0476) + IT_0552*IT_0553 + IT_0549*IT_0554) + 2*IT_0365*(conj
      (IT_0170)*IT_0302 + IT_0306*IT_0401 + IT_0535*IT_0555) + (-4)*conj(IT_0516
      )*(IT_0474*IT_0516 + -1./2*IT_0466*IT_0520 + -1./2*IT_0511*IT_0543 + (-2)
      *IT_0475*IT_0545 + 1./2*IT_0478*IT_0548 + 1./2*IT_0547*IT_0556) + 8
      *IT_0516*(conj(IT_0475)*IT_0545 + (-0.25)*IT_0480*IT_0548 + (-0.25)
      *IT_0551*IT_0556) + IT_0535*(IT_0436*IT_0494 + IT_0557) + 2*conj(IT_0299)*
      (IT_0176*IT_0299 + IT_0170*IT_0301 + IT_0219*IT_0365 + IT_0302*IT_0404 +
       IT_0437*IT_0494 + IT_0557) + IT_0558*(IT_0224*IT_0303 + IT_0559 + IT_0560
       + IT_0561) + 2*conj(IT_0303)*(IT_0224*IT_0497 + IT_0503*IT_0555 + IT_0559
       + IT_0560 + IT_0561);
    return create_ccomplex_return(IT_0562);
}

