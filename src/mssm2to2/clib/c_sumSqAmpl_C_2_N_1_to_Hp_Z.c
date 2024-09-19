#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_1_to_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_1_to_Hp_Z(
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 2*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0000*IT_0008;
    const ccomplex_t IT_0010 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = IT_0012 + IT_0015;
    const ccomplex_t IT_0017 = IT_0006 + IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0008;
    const ccomplex_t IT_0021 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0002*IT_0007;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0022 + -IT_0024 + (-2)
      *IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0019*IT_0028;
    const ccomplex_t IT_0030 = pow(m_Hp, 2);
    const ccomplex_t IT_0031 = pow(m_N_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0030 + IT_0031 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0029*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0036;
    const ccomplex_t IT_0038 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0042 = IT_0003*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = IT_0040 + IT_0043;
    const ccomplex_t IT_0045 = IT_0037 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0049 = IT_0020*IT_0048;
    const ccomplex_t IT_0050 = IT_0023*IT_0048;
    const ccomplex_t IT_0051 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0052 = IT_0023*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + -IT_0050 + (-2)
      *IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0030 + IT_0031 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*-1)*e_em*(IT_0002*IT_0007 + 
      -IT_0001*IT_0008);
    const ccomplex_t IT_0060 = pow(m_C_2, 2);
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + -IT_0031 + -IT_0060 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0019*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0059*IT_0063;
    const ccomplex_t IT_0065 = 2*IT_0034 + 2*IT_0058 + IT_0064;
    const ccomplex_t IT_0066 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0067 = IT_0020*IT_0066;
    const ccomplex_t IT_0068 = IT_0023*IT_0066;
    const ccomplex_t IT_0069 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0070 = IT_0023*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + -IT_0068 + (-2)
      *IT_0070);
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = sin(beta);
    const ccomplex_t IT_0074 = IT_0008*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = IT_0002*IT_0073;
    const ccomplex_t IT_0079 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0077 + IT_0081;
    const ccomplex_t IT_0083 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0084 = IT_0078*IT_0083;
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = IT_0082 + IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0072*IT_0088;
    const ccomplex_t IT_0090 = IT_0056*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0093 = IT_0078*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0096 = IT_0074*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = IT_0094 + IT_0097;
    const ccomplex_t IT_0099 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0100 = IT_0078*IT_0099;
    const ccomplex_t IT_0101 = (-2)*IT_0100;
    const ccomplex_t IT_0102 = IT_0098 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0061*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0059*IT_0106;
    const ccomplex_t IT_0108 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0109 = IT_0020*IT_0108;
    const ccomplex_t IT_0110 = IT_0023*IT_0108;
    const ccomplex_t IT_0111 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0112 = IT_0023*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + -IT_0110 + (-2)
      *IT_0112);
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0104*IT_0114;
    const ccomplex_t IT_0116 = IT_0032*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 2*IT_0091 + -IT_0107 + 2*IT_0117;
    const ccomplex_t IT_0119 = pow(m_Z, -2);
    const ccomplex_t IT_0120 = pow(s_14, 2);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = -IT_0060;
    const ccomplex_t IT_0123 = IT_0121 + IT_0122;
    const ccomplex_t IT_0124 = s_12*IT_0123;
    const ccomplex_t IT_0125 = m_C_2*m_N_1;
    const ccomplex_t IT_0126 = pow(s_24, 2);
    const ccomplex_t IT_0127 = IT_0119*IT_0126;
    const ccomplex_t IT_0128 = -IT_0031;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0125*IT_0129;
    const ccomplex_t IT_0131 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0132 = IT_0020*IT_0131;
    const ccomplex_t IT_0133 = IT_0023*IT_0131;
    const ccomplex_t IT_0134 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0135 = IT_0020*IT_0134;
    const ccomplex_t IT_0136 = IT_0023*IT_0134;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0132 + IT_0133 + 
      -IT_0135 + -IT_0136);
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = cpow((-2)*s_13 + IT_0030 + IT_0060 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0140 = IT_0104*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = -IT_0107 + IT_0142;
    const ccomplex_t IT_0144 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0145 = IT_0020*IT_0144;
    const ccomplex_t IT_0146 = IT_0023*IT_0144;
    const ccomplex_t IT_0147 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = IT_0023*IT_0147;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0145 + IT_0146 + 
      -IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0153 = IT_0003*IT_0152;
    const ccomplex_t IT_0154 = 2*IT_0153;
    const ccomplex_t IT_0155 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0156 = IT_0009*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = IT_0157 + IT_0160;
    const ccomplex_t IT_0162 = IT_0154 + IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_13 + IT_0030 + IT_0060 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0151*IT_0167;
    const ccomplex_t IT_0169 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0170 = IT_0020*IT_0169;
    const ccomplex_t IT_0171 = IT_0023*IT_0169;
    const ccomplex_t IT_0172 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0173 = IT_0020*IT_0172;
    const ccomplex_t IT_0174 = IT_0023*IT_0172;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0170 + IT_0171 + 
      -IT_0173 + -IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0178 = IT_0003*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0181 = IT_0009*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0184 = IT_0003*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = IT_0182 + IT_0185;
    const ccomplex_t IT_0187 = IT_0179 + IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = cpow((-2)*s_13 + IT_0030 + IT_0060 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0191 = IT_0189*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0176*IT_0192;
    const ccomplex_t IT_0194 = IT_0019*IT_0139;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0138*IT_0195;
    const ccomplex_t IT_0197 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0198 = IT_0020*IT_0197;
    const ccomplex_t IT_0199 = IT_0023*IT_0197;
    const ccomplex_t IT_0200 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0201 = IT_0020*IT_0200;
    const ccomplex_t IT_0202 = IT_0023*IT_0200;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0198 + IT_0199 + 
      -IT_0201 + -IT_0202);
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0206 = IT_0003*IT_0205;
    const ccomplex_t IT_0207 = 2*IT_0206;
    const ccomplex_t IT_0208 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0209 = IT_0009*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0212 = IT_0003*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = IT_0210 + IT_0213;
    const ccomplex_t IT_0215 = IT_0207 + IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_13 + IT_0030 + IT_0060 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = IT_0204*IT_0220;
    const ccomplex_t IT_0222 = (-2)*IT_0168 + (-2)*IT_0193 + -IT_0196 + (-2)
      *IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0137;
    const ccomplex_t IT_0224 = IT_0195*IT_0223;
    const ccomplex_t IT_0225 = IT_0064 + -IT_0224;
    const ccomplex_t IT_0226 = (-2)*IT_0225;
    const ccomplex_t IT_0227 = 0.5*IT_0226;
    const ccomplex_t IT_0228 = IT_0222 + IT_0227;
    const ccomplex_t IT_0229 = (-2)*conj(IT_0225);
    const ccomplex_t IT_0230 = 0.5*IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0222) + IT_0230;
    const ccomplex_t IT_0232 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0233 = IT_0020*IT_0232;
    const ccomplex_t IT_0234 = IT_0023*IT_0232;
    const ccomplex_t IT_0235 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0236 = IT_0020*IT_0235;
    const ccomplex_t IT_0237 = IT_0023*IT_0235;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*(IT_0233 + IT_0234 + 
      -IT_0236 + -IT_0237);
    const ccomplex_t IT_0239 = (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0241 = IT_0078*IT_0240;
    const ccomplex_t IT_0242 = (-2)*IT_0241;
    const ccomplex_t IT_0243 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0244 = IT_0074*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0247 = IT_0078*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = IT_0245 + IT_0248;
    const ccomplex_t IT_0250 = IT_0242 + IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0251;
    const ccomplex_t IT_0253 = IT_0165*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0239*IT_0254;
    const ccomplex_t IT_0256 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0257 = IT_0020*IT_0256;
    const ccomplex_t IT_0258 = IT_0023*IT_0256;
    const ccomplex_t IT_0259 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0260 = IT_0020*IT_0259;
    const ccomplex_t IT_0261 = IT_0023*IT_0259;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*(IT_0257 + IT_0258 + 
      -IT_0260 + -IT_0261);
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0265 = IT_0078*IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0268 = IT_0074*IT_0267;
    const ccomplex_t IT_0269 = 1.4142135623731*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0271 = IT_0078*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = IT_0269 + IT_0272;
    const ccomplex_t IT_0274 = IT_0266 + IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = IT_0190*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = IT_0263*IT_0278;
    const ccomplex_t IT_0280 = IT_0141*IT_0223;
    const ccomplex_t IT_0281 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0282 = IT_0020*IT_0281;
    const ccomplex_t IT_0283 = IT_0023*IT_0281;
    const ccomplex_t IT_0284 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0285 = IT_0020*IT_0284;
    const ccomplex_t IT_0286 = IT_0023*IT_0284;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*(IT_0282 + IT_0283 + 
      -IT_0285 + -IT_0286);
    const ccomplex_t IT_0288 = (-0.5)*IT_0287;
    const ccomplex_t IT_0289 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0290 = IT_0078*IT_0289;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0293 = IT_0074*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0296 = IT_0078*IT_0295;
    const ccomplex_t IT_0297 = 1.4142135623731*IT_0296;
    const ccomplex_t IT_0298 = IT_0294 + IT_0297;
    const ccomplex_t IT_0299 = IT_0291 + IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = 0.5*IT_0300;
    const ccomplex_t IT_0302 = IT_0218*IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0302;
    const ccomplex_t IT_0304 = IT_0288*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0255 + 2*IT_0279 + IT_0280 + 2*IT_0304;
    const ccomplex_t IT_0306 = (-2)*IT_0222;
    const ccomplex_t IT_0307 = 0.5*IT_0306;
    const ccomplex_t IT_0308 = IT_0225 + IT_0307;
    const ccomplex_t IT_0309 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0310 = 0.5*IT_0309;
    const ccomplex_t IT_0311 = conj(IT_0225) + IT_0310;
    const ccomplex_t IT_0312 = pow(m_Z, 2);
    const ccomplex_t IT_0313 = s_24*IT_0312;
    const ccomplex_t IT_0314 = IT_0119*IT_0313;
    const ccomplex_t IT_0315 = (-0.25)*IT_0314;
    const ccomplex_t IT_0316 = s_24 + IT_0315;
    const ccomplex_t IT_0317 = m_C_2*IT_0316;
    const ccomplex_t IT_0318 = (-4)*IT_0317;
    const ccomplex_t IT_0319 = -IT_0058;
    const ccomplex_t IT_0320 = -IT_0034;
    const ccomplex_t IT_0321 = conj(IT_0319) + conj(IT_0320);
    const ccomplex_t IT_0322 = m_N_1*IT_0223;
    const ccomplex_t IT_0323 = m_N_1*IT_0138;
    const ccomplex_t IT_0324 = -IT_0323;
    const ccomplex_t IT_0325 = IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = IT_0141*IT_0325;
    const ccomplex_t IT_0327 = m_N_2*IT_0151;
    const ccomplex_t IT_0328 = m_N_1*IT_0239;
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = IT_0327 + IT_0329;
    const ccomplex_t IT_0331 = IT_0254*IT_0330;
    const ccomplex_t IT_0332 = m_N_3*IT_0176;
    const ccomplex_t IT_0333 = m_N_1*IT_0263;
    const ccomplex_t IT_0334 = -IT_0333;
    const ccomplex_t IT_0335 = IT_0332 + IT_0334;
    const ccomplex_t IT_0336 = IT_0278*IT_0335;
    const ccomplex_t IT_0337 = m_N_4*IT_0204;
    const ccomplex_t IT_0338 = m_N_1*IT_0288;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0303*IT_0340;
    const ccomplex_t IT_0342 = (-0.5)*IT_0326 + IT_0331 + IT_0336 + IT_0341;
    const ccomplex_t IT_0343 = m_C_2*IT_0034;
    const ccomplex_t IT_0344 = m_C_2*IT_0058;
    const ccomplex_t IT_0345 = IT_0047*IT_0072;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0347 = IT_0056*IT_0345*IT_0346;
    const ccomplex_t IT_0348 = IT_0019*IT_0114;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0350 = IT_0032*IT_0348*IT_0349;
    const ccomplex_t IT_0351 = (-0.5)*IT_0326 + -IT_0343 + -IT_0344 + -IT_0347
       + -IT_0350;
    const ccomplex_t IT_0352 = IT_0342 + IT_0351;
    const ccomplex_t IT_0353 = IT_0319 + IT_0320;
    const ccomplex_t IT_0354 = conj(IT_0342) + conj(IT_0351);
    const ccomplex_t IT_0355 = IT_0195*IT_0325;
    const ccomplex_t IT_0356 = m_N_2*IT_0239;
    const ccomplex_t IT_0357 = m_N_1*IT_0151;
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = IT_0356 + IT_0358;
    const ccomplex_t IT_0360 = IT_0167*IT_0359;
    const ccomplex_t IT_0361 = m_N_3*IT_0263;
    const ccomplex_t IT_0362 = m_N_1*IT_0176;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = IT_0361 + IT_0363;
    const ccomplex_t IT_0365 = IT_0192*IT_0364;
    const ccomplex_t IT_0366 = m_N_4*IT_0288;
    const ccomplex_t IT_0367 = m_N_1*IT_0204;
    const ccomplex_t IT_0368 = -IT_0367;
    const ccomplex_t IT_0369 = IT_0366 + IT_0368;
    const ccomplex_t IT_0370 = IT_0220*IT_0369;
    const ccomplex_t IT_0371 = (-0.5)*IT_0355 + -IT_0360 + -IT_0365 + -IT_0370;
    const ccomplex_t IT_0372 = IT_0028*IT_0104;
    const ccomplex_t IT_0373 = IT_0032*IT_0349*IT_0372;
    const ccomplex_t IT_0374 = IT_0054*IT_0088;
    const ccomplex_t IT_0375 = IT_0056*IT_0346*IT_0374;
    const ccomplex_t IT_0376 = m_C_2*IT_0091;
    const ccomplex_t IT_0377 = m_C_2*IT_0117;
    const ccomplex_t IT_0378 = (-0.5)*IT_0355 + -IT_0373 + -IT_0375 + -IT_0376
       + -IT_0377;
    const ccomplex_t IT_0379 = conj(IT_0371) + conj(IT_0378);
    const ccomplex_t IT_0380 = -IT_0117;
    const ccomplex_t IT_0381 = -IT_0091;
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = IT_0371 + IT_0378;
    const ccomplex_t IT_0384 = conj(IT_0380) + conj(IT_0381);
    const ccomplex_t IT_0385 = 0.5*IT_0314;
    const ccomplex_t IT_0386 = s_24 + IT_0385;
    const ccomplex_t IT_0387 = m_C_2*IT_0386;
    const ccomplex_t IT_0388 = 2*IT_0387;
    const ccomplex_t IT_0389 = 0.5*IT_0224;
    const ccomplex_t IT_0390 = (-0.5)*IT_0142;
    const ccomplex_t IT_0391 = IT_0168 + IT_0193 + 0.5*IT_0196 + IT_0221;
    const ccomplex_t IT_0392 = (-2)*IT_0342;
    const ccomplex_t IT_0393 = (-2)*IT_0351;
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = (-2)*conj(IT_0342);
    const ccomplex_t IT_0396 = (-2)*conj(IT_0351);
    const ccomplex_t IT_0397 = IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = -IT_0255 + -IT_0279 + (-0.5)*IT_0280 + -IT_0304;
    const ccomplex_t IT_0399 = (-2)*IT_0371;
    const ccomplex_t IT_0400 = (-2)*IT_0378;
    const ccomplex_t IT_0401 = IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = (-2)*conj(IT_0371);
    const ccomplex_t IT_0403 = (-2)*conj(IT_0378);
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 2*IT_0388*(IT_0354*IT_0389 + IT_0352*conj
      (IT_0389) + IT_0379*IT_0390 + IT_0383*conj(IT_0390) + 0.5*conj(IT_0391)
      *IT_0394 + 0.5*IT_0391*IT_0397 + 0.5*conj(IT_0398)*IT_0401 + 0.5*IT_0398
      *IT_0404);
    const ccomplex_t IT_0406 = s_14*IT_0312;
    const ccomplex_t IT_0407 = IT_0119*IT_0406;
    const ccomplex_t IT_0408 = 0.5*IT_0407;
    const ccomplex_t IT_0409 = s_14 + IT_0408;
    const ccomplex_t IT_0410 = m_N_1*IT_0409;
    const ccomplex_t IT_0411 = 2*IT_0410;
    const ccomplex_t IT_0412 = s_12*IT_0312;
    const ccomplex_t IT_0413 = s_14*s_24;
    const ccomplex_t IT_0414 = (-2)*IT_0413;
    const ccomplex_t IT_0415 = IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = IT_0119*IT_0415;
    const ccomplex_t IT_0417 = (-0.5)*IT_0416;
    const ccomplex_t IT_0418 = s_12 + IT_0417;
    const ccomplex_t IT_0419 = 4*IT_0371;
    const ccomplex_t IT_0420 = 6*conj(IT_0342);
    const ccomplex_t IT_0421 = pow(m_Z, 4);
    const ccomplex_t IT_0422 = s_12*IT_0421;
    const ccomplex_t IT_0423 = IT_0119*IT_0422;
    const ccomplex_t IT_0424 = (-0.25)*IT_0423;
    const ccomplex_t IT_0425 = IT_0413 + IT_0424;
    const ccomplex_t IT_0426 = -IT_0389;
    const ccomplex_t IT_0427 = IT_0391 + IT_0426;
    const ccomplex_t IT_0428 = 8*conj(IT_0391);
    const ccomplex_t IT_0429 = 8*IT_0398;
    const ccomplex_t IT_0430 = 8*conj(IT_0398);
    const ccomplex_t IT_0431 = IT_0123*IT_0125;
    const ccomplex_t IT_0432 = (-2)*IT_0065;
    const ccomplex_t IT_0433 = (-2)*conj(IT_0065);
    const ccomplex_t IT_0434 = s_14*s_24*IT_0119;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = s_12 + IT_0435;
    const ccomplex_t IT_0437 = m_C_2*IT_0436;
    const ccomplex_t IT_0438 = (-2)*IT_0118;
    const ccomplex_t IT_0439 = (-2)*conj(IT_0118);
    const ccomplex_t IT_0440 = s_12*s_24;
    const ccomplex_t IT_0441 = s_14*IT_0031;
    const ccomplex_t IT_0442 = -IT_0441;
    const ccomplex_t IT_0443 = IT_0440 + IT_0442;
    const ccomplex_t IT_0444 = (-2)*IT_0305;
    const ccomplex_t IT_0445 = 0.5*IT_0444;
    const ccomplex_t IT_0446 = IT_0143 + IT_0445;
    const ccomplex_t IT_0447 = (-2)*conj(IT_0305);
    const ccomplex_t IT_0448 = 0.5*IT_0447;
    const ccomplex_t IT_0449 = conj(IT_0143) + IT_0448;
    const ccomplex_t IT_0450 = conj(IT_0380) + conj(IT_0381) + conj(IT_0398);
    const ccomplex_t IT_0451 = (-2)*IT_0143;
    const ccomplex_t IT_0452 = IT_0380 + IT_0381 + IT_0398;
    const ccomplex_t IT_0453 = (-2)*conj(IT_0143);
    const ccomplex_t IT_0454 = IT_0311*IT_0389 + IT_0308*conj(IT_0389) +
       IT_0231*(IT_0319 + IT_0320 + IT_0391) + IT_0228*(conj(IT_0319) + conj
      (IT_0320) + conj(IT_0391)) + conj(IT_0390)*IT_0446 + IT_0390*IT_0449 +
       IT_0450*(IT_0305 + 0.5*IT_0451) + IT_0452*(conj(IT_0305) + 0.5*IT_0453);
    const ccomplex_t IT_0455 = m_N_1*IT_0123;
    const ccomplex_t IT_0456 = IT_0125*IT_0312;
    const ccomplex_t IT_0457 = (-3)*IT_0456;
    const ccomplex_t IT_0458 = (-2)*IT_0389;
    const ccomplex_t IT_0459 = 0.5*IT_0458;
    const ccomplex_t IT_0460 = IT_0391 + IT_0459;
    const ccomplex_t IT_0461 = (-2)*conj(IT_0389);
    const ccomplex_t IT_0462 = 0.5*IT_0461;
    const ccomplex_t IT_0463 = conj(IT_0391) + IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0319;
    const ccomplex_t IT_0465 = (-2)*IT_0320;
    const ccomplex_t IT_0466 = IT_0464 + IT_0465;
    const ccomplex_t IT_0467 = (-2)*conj(IT_0319);
    const ccomplex_t IT_0468 = (-2)*conj(IT_0320);
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = -conj(IT_0389);
    const ccomplex_t IT_0471 = (-0.5)*IT_0423;
    const ccomplex_t IT_0472 = IT_0412 + IT_0414 + IT_0471;
    const ccomplex_t IT_0473 = (-4)*IT_0380;
    const ccomplex_t IT_0474 = 3*IT_0456;
    const ccomplex_t IT_0475 = 0.5*IT_0390;
    const ccomplex_t IT_0476 = 0.5*conj(IT_0390);
    const ccomplex_t IT_0477 = IT_0125*IT_0436;
    const ccomplex_t IT_0478 = 0.5*IT_0305;
    const ccomplex_t IT_0479 = 0.5*conj(IT_0305);
    const ccomplex_t IT_0480 = s_12*IT_0436;
    const ccomplex_t IT_0481 = 0.5*conj(IT_0225);
    const ccomplex_t IT_0482 = 0.5*conj(IT_0143);
    const ccomplex_t IT_0483 = s_12*IT_0129;
    const ccomplex_t IT_0484 = m_N_1*IT_0436;
    const ccomplex_t IT_0485 = 0.5*IT_0383;
    const ccomplex_t IT_0486 = 2*IT_0484*(conj(IT_0143)*IT_0352 + IT_0143
      *IT_0354 + IT_0308*IT_0379 + conj(IT_0225)*IT_0383 + 0.5*conj(IT_0351)
      *IT_0444 + 0.5*IT_0351*IT_0447 + IT_0395*IT_0478 + IT_0392*IT_0479 +
       IT_0309*IT_0485);
    const ccomplex_t IT_0487 = m_C_2*IT_0129;
    const ccomplex_t IT_0488 = (-0.25)*IT_0407;
    const ccomplex_t IT_0489 = s_14 + IT_0488;
    const ccomplex_t IT_0490 = m_N_1*IT_0489;
    const ccomplex_t IT_0491 = (-4)*IT_0490;
    const ccomplex_t IT_0492 = 0.5*IT_0389;
    const ccomplex_t IT_0493 = 0.5*conj(IT_0389);
    const ccomplex_t IT_0494 = s_12*s_14;
    const ccomplex_t IT_0495 = s_24*IT_0060;
    const ccomplex_t IT_0496 = -IT_0495;
    const ccomplex_t IT_0497 = IT_0494 + IT_0496;
    const ccomplex_t IT_0498 = 2*(conj(IT_0065)*(IT_0319 + IT_0320 + IT_0391) 
      + IT_0065*(conj(IT_0319) + conj(IT_0320) + conj(IT_0391)) + IT_0118
      *IT_0450 + conj(IT_0118)*IT_0452 + IT_0439*IT_0475 + IT_0438*IT_0476 +
       IT_0433*IT_0492 + IT_0432*IT_0493)*IT_0497;
    const ccomplex_t IT_0499 = 2*(IT_0065*conj(IT_0065) + IT_0118*conj(IT_0118
      ))*IT_0124 + 2*IT_0130*(conj(IT_0143)*IT_0228 + IT_0143*IT_0231 + conj
      (IT_0305)*IT_0308 + IT_0305*IT_0311) + 2*IT_0318*(IT_0321*IT_0352 +
       IT_0353*IT_0354 + IT_0379*IT_0382 + IT_0383*IT_0384) + IT_0405 + (IT_0384
      *IT_0394 + IT_0382*IT_0397 + IT_0321*IT_0401 + IT_0353*IT_0404)*IT_0411 +
       4*IT_0418*(IT_0352*IT_0354 + IT_0379*(IT_0378 + 0.25*IT_0419)) + 6
      *IT_0125*(IT_0352*IT_0379 + IT_0383*(conj(IT_0351) + 0.166666666666667
      *IT_0420)) + 8*IT_0425*(IT_0321*IT_0427 + -IT_0353*(conj(IT_0389) + (
      -0.125)*IT_0428) + -IT_0384*(IT_0390 + (-0.125)*IT_0429) + -IT_0382*(conj
      (IT_0390) + (-0.125)*IT_0430)) + IT_0431*(conj(IT_0118)*IT_0432 + IT_0118
      *IT_0433) + IT_0437*(IT_0354*IT_0432 + IT_0352*IT_0433 + IT_0379*IT_0438 +
       IT_0383*IT_0439) + 2*IT_0443*IT_0454 + (conj(IT_0118)*IT_0392 + IT_0118
      *IT_0395 + IT_0379*IT_0432 + IT_0383*IT_0433 + conj(IT_0351)*IT_0438 +
       IT_0351*IT_0439)*IT_0455 + (-2)*IT_0457*(conj(IT_0398)*IT_0427 + -conj
      (IT_0390)*IT_0460 + -IT_0390*IT_0463 + (-0.5)*IT_0384*IT_0466 + (-0.5)
      *IT_0382*IT_0469 + IT_0398*(conj(IT_0391) + IT_0470)) + (-4)*IT_0472*
      (IT_0321*IT_0353 + (IT_0390 + -IT_0398)*(conj(IT_0390) + -conj(IT_0398)) +
       IT_0427*(conj(IT_0391) + IT_0470) + IT_0384*(IT_0381 + (-0.25)*IT_0473)) 
      + 2*IT_0474*(IT_0321*IT_0398 + IT_0353*conj(IT_0398) + IT_0384*IT_0460 +
       IT_0382*IT_0463 + IT_0469*IT_0475 + IT_0466*IT_0476) + 2*IT_0477*(conj
      (IT_0065)*IT_0143 + IT_0065*conj(IT_0143) + conj(IT_0118)*IT_0308 +
       IT_0118*IT_0311 + IT_0433*IT_0478 + IT_0432*IT_0479) + 2*IT_0480*(conj
      (IT_0065)*IT_0222 + IT_0065*conj(IT_0222) + conj(IT_0118)*IT_0305 +
       IT_0118*conj(IT_0305) + 0.5*IT_0225*IT_0433 + 0.5*IT_0143*IT_0439 +
       IT_0432*IT_0481 + IT_0438*IT_0482) + 2*(IT_0222*conj(IT_0222) + IT_0305
      *conj(IT_0305) + IT_0225*IT_0311 + IT_0143*IT_0449 + IT_0306*IT_0481 +
       IT_0444*IT_0482)*IT_0483 + IT_0486 + 2*(IT_0311*IT_0352 + IT_0308*IT_0354
       + conj(IT_0143)*IT_0383 + IT_0379*IT_0446 + IT_0447*IT_0485)*IT_0487 + 2
      *IT_0491*(IT_0379*IT_0391 + IT_0383*conj(IT_0391) + IT_0354*IT_0398 +
       IT_0352*conj(IT_0398) + IT_0397*IT_0475 + IT_0394*IT_0476 + IT_0404
      *IT_0492 + IT_0401*IT_0493) + IT_0498;
    return create_ccomplex_return(IT_0499);
}

