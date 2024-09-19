#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_1_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_1_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
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
    const ccomplex_t IT_0001 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0012 + 0.5*IT_0015);
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = pow(m_A0, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0019 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0022 = IT_0018*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0024 = IT_0000*IT_0023;
    const ccomplex_t IT_0025 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*(IT_0024 + IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0031 = IT_0007*IT_0030;
    const ccomplex_t IT_0032 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 0.5*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0039 = pow(m_C_1, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0019 + IT_0039 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = IT_0000*IT_0043;
    const ccomplex_t IT_0045 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = IT_0003*IT_0043;
    const ccomplex_t IT_0048 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = IT_0000*IT_0007;
    const ccomplex_t IT_0051 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = IT_0003*IT_0007;
    const ccomplex_t IT_0054 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = -IT_0046 + IT_0049 + IT_0052 + -IT_0055;
    const ccomplex_t IT_0057 = m_N_1*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0019 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0017*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0057*IT_0060;
    const ccomplex_t IT_0062 = IT_0022 + IT_0041 + (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0064 = IT_0044*IT_0063;
    const ccomplex_t IT_0065 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0066 = IT_0047*IT_0065;
    const ccomplex_t IT_0067 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0068 = IT_0050*IT_0067;
    const ccomplex_t IT_0069 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0070 = IT_0053*IT_0069;
    const ccomplex_t IT_0071 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0072 = IT_0044*IT_0071;
    const ccomplex_t IT_0073 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0074 = IT_0047*IT_0073;
    const ccomplex_t IT_0075 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0076 = IT_0050*IT_0075;
    const ccomplex_t IT_0077 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0078 = IT_0053*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0064 + 0.5*IT_0066 + 0.5*IT_0068 + (
      -0.5)*IT_0070 + (-0.5)*IT_0072 + 0.5*IT_0074 + 0.5*IT_0076 + (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = m_N_3*IT_0079;
    const ccomplex_t IT_0081 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0082 = IT_0007*IT_0081;
    const ccomplex_t IT_0083 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0084 = IT_0007*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + 0.5*IT_0085);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = pow(m_N_1, 2);
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0019 + IT_0088 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0080*IT_0091;
    const ccomplex_t IT_0093 = m_C_1*IT_0056;
    const ccomplex_t IT_0094 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0095 = IT_0007*IT_0094;
    const ccomplex_t IT_0096 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0097 = IT_0007*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + (-0.5)*IT_0098);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0058*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0093*IT_0102;
    const ccomplex_t IT_0104 = IT_0010*IT_0100;
    const ccomplex_t IT_0105 = IT_0020*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = m_N_1*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0109 = IT_0007*IT_0108;
    const ccomplex_t IT_0110 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + (-0.5)*IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0029*IT_0114;
    const ccomplex_t IT_0116 = IT_0040*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = m_N_1*IT_0117;
    const ccomplex_t IT_0119 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0120 = IT_0044*IT_0119;
    const ccomplex_t IT_0121 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0122 = IT_0047*IT_0121;
    const ccomplex_t IT_0123 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0124 = IT_0050*IT_0123;
    const ccomplex_t IT_0125 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0126 = IT_0053*IT_0125;
    const ccomplex_t IT_0127 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0128 = IT_0044*IT_0127;
    const ccomplex_t IT_0129 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0130 = IT_0047*IT_0129;
    const ccomplex_t IT_0131 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0132 = IT_0050*IT_0131;
    const ccomplex_t IT_0133 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0134 = IT_0053*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0120 + 0.5*IT_0122 + 0.5*IT_0124 + (
      -0.5)*IT_0126 + (-0.5)*IT_0128 + 0.5*IT_0130 + 0.5*IT_0132 + (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = m_N_2*IT_0135;
    const ccomplex_t IT_0137 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0138 = IT_0007*IT_0137;
    const ccomplex_t IT_0139 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0140 = IT_0007*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 0.5*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_23 + IT_0019 + IT_0088 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0136*IT_0146;
    const ccomplex_t IT_0148 = m_C_1*IT_0135;
    const ccomplex_t IT_0149 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0150 = IT_0007*IT_0149;
    const ccomplex_t IT_0151 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0152 = IT_0007*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + (-0.5)*IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0144*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = IT_0148*IT_0157;
    const ccomplex_t IT_0159 = m_C_1*IT_0079;
    const ccomplex_t IT_0160 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0161 = IT_0007*IT_0160;
    const ccomplex_t IT_0162 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0163 = IT_0007*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + (-0.5)*IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0089*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0159*IT_0168;
    const ccomplex_t IT_0170 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0171 = IT_0044*IT_0170;
    const ccomplex_t IT_0172 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0173 = IT_0047*IT_0172;
    const ccomplex_t IT_0174 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0175 = IT_0050*IT_0174;
    const ccomplex_t IT_0176 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0177 = IT_0053*IT_0176;
    const ccomplex_t IT_0178 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0179 = IT_0044*IT_0178;
    const ccomplex_t IT_0180 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0181 = IT_0047*IT_0180;
    const ccomplex_t IT_0182 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0183 = IT_0050*IT_0182;
    const ccomplex_t IT_0184 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0185 = IT_0053*IT_0184;
    const ccomplex_t IT_0186 = (-0.5)*IT_0171 + 0.5*IT_0173 + 0.5*IT_0175 + (
      -0.5)*IT_0177 + (-0.5)*IT_0179 + 0.5*IT_0181 + 0.5*IT_0183 + (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = m_N_4*IT_0186;
    const ccomplex_t IT_0188 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0189 = IT_0007*IT_0188;
    const ccomplex_t IT_0190 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0191 = IT_0007*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0189 + 0.5*IT_0192);
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = cpow((-2)*s_23 + IT_0019 + IT_0088 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0196 = IT_0194*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = IT_0187*IT_0197;
    const ccomplex_t IT_0199 = m_C_1*IT_0186;
    const ccomplex_t IT_0200 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0201 = IT_0007*IT_0200;
    const ccomplex_t IT_0202 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0203 = IT_0007*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0201 + (-0.5)*IT_0204);
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = IT_0195*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0199*IT_0208;
    const ccomplex_t IT_0210 = (-0.5)*IT_0061 + -IT_0092 + IT_0103 + -IT_0107 
      + -IT_0118 + -IT_0147 + IT_0158 + IT_0169 + -IT_0198 + IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0062) + conj(IT_0210);
    const ccomplex_t IT_0212 = pow(m_W, -2);
    const ccomplex_t IT_0213 = pow(m_W, 2);
    const ccomplex_t IT_0214 = s_12*IT_0213;
    const ccomplex_t IT_0215 = s_14*s_24;
    const ccomplex_t IT_0216 = (-2)*IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = IT_0212*IT_0217;
    const ccomplex_t IT_0219 = (-0.5)*IT_0218;
    const ccomplex_t IT_0220 = s_12 + IT_0219;
    const ccomplex_t IT_0221 = 4*IT_0220;
    const ccomplex_t IT_0222 = IT_0062*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0224 = IT_0044*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0226 = IT_0047*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0228 = IT_0050*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0230 = IT_0053*IT_0229;
    const ccomplex_t IT_0231 = IT_0224 + -IT_0226 + -IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0060*IT_0231;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = IT_0106 + IT_0117;
    const ccomplex_t IT_0235 = conj(IT_0233) + conj(IT_0234);
    const ccomplex_t IT_0236 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0237 = IT_0044*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0239 = IT_0047*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0241 = IT_0050*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0243 = IT_0053*IT_0242;
    const ccomplex_t IT_0244 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0245 = IT_0044*IT_0244;
    const ccomplex_t IT_0246 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0247 = IT_0047*IT_0246;
    const ccomplex_t IT_0248 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0249 = IT_0050*IT_0248;
    const ccomplex_t IT_0250 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0251 = IT_0053*IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0237 + (-0.5)*IT_0239 + (-0.5)*IT_0241 +
       0.5*IT_0243 + 0.5*IT_0245 + (-0.5)*IT_0247 + (-0.5)*IT_0249 + 0.5*IT_0251;
    const ccomplex_t IT_0253 = IT_0091*IT_0252;
    const ccomplex_t IT_0254 = e_em*IT_0007;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = cpow((-2)*s_12 + -IT_0039 + -IT_0088 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0257 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0258 = IT_0047*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0261 = IT_0053*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = IT_0259 + IT_0262;
    const ccomplex_t IT_0264 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0265 = IT_0053*IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = IT_0263 + IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = IT_0256*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0255*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0274 = IT_0044*IT_0273;
    const ccomplex_t IT_0275 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0276 = IT_0047*IT_0275;
    const ccomplex_t IT_0277 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0278 = IT_0050*IT_0277;
    const ccomplex_t IT_0279 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0280 = IT_0053*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0282 = IT_0044*IT_0281;
    const ccomplex_t IT_0283 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0284 = IT_0047*IT_0283;
    const ccomplex_t IT_0285 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0286 = IT_0050*IT_0285;
    const ccomplex_t IT_0287 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0288 = IT_0053*IT_0287;
    const ccomplex_t IT_0289 = 0.5*IT_0274 + (-0.5)*IT_0276 + (-0.5)*IT_0278 +
       0.5*IT_0280 + 0.5*IT_0282 + (-0.5)*IT_0284 + (-0.5)*IT_0286 + 0.5*IT_0288;
    const ccomplex_t IT_0290 = IT_0146*IT_0289;
    const ccomplex_t IT_0291 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0292 = IT_0044*IT_0291;
    const ccomplex_t IT_0293 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0294 = IT_0047*IT_0293;
    const ccomplex_t IT_0295 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0296 = IT_0050*IT_0295;
    const ccomplex_t IT_0297 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0298 = IT_0053*IT_0297;
    const ccomplex_t IT_0299 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0300 = IT_0044*IT_0299;
    const ccomplex_t IT_0301 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0302 = IT_0047*IT_0301;
    const ccomplex_t IT_0303 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0304 = IT_0050*IT_0303;
    const ccomplex_t IT_0305 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0306 = IT_0053*IT_0305;
    const ccomplex_t IT_0307 = 0.5*IT_0292 + (-0.5)*IT_0294 + (-0.5)*IT_0296 +
       0.5*IT_0298 + 0.5*IT_0300 + (-0.5)*IT_0302 + (-0.5)*IT_0304 + 0.5*IT_0306;
    const ccomplex_t IT_0308 = IT_0197*IT_0307;
    const ccomplex_t IT_0309 = -IT_0232 + (-2)*IT_0253 + -IT_0272 + (-2)
      *IT_0290 + (-2)*IT_0308;
    const ccomplex_t IT_0310 = s_12*s_14;
    const ccomplex_t IT_0311 = s_24*IT_0039;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0313;
    const ccomplex_t IT_0315 = IT_0309*IT_0314;
    const ccomplex_t IT_0316 = m_C_1*IT_0231;
    const ccomplex_t IT_0317 = IT_0060*IT_0316;
    const ccomplex_t IT_0318 = V_u1*U_Wm1;
    const ccomplex_t IT_0319 = IT_0000*IT_0318;
    const ccomplex_t IT_0320 = U_d1*V_Wp1;
    const ccomplex_t IT_0321 = IT_0003*IT_0320;
    const ccomplex_t IT_0322 = IT_0007*(IT_0319 + IT_0321);
    const ccomplex_t IT_0323 = 1.4142135623731*e_em*IT_0322;
    const ccomplex_t IT_0324 = 0.5*IT_0323;
    const ccomplex_t IT_0325 = IT_0017*IT_0324;
    const ccomplex_t IT_0326 = IT_0020*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*IT_0326;
    const ccomplex_t IT_0328 = m_N_1*IT_0327;
    const ccomplex_t IT_0329 = V_u1*U_Wm2;
    const ccomplex_t IT_0330 = IT_0000*IT_0329;
    const ccomplex_t IT_0331 = U_d2*V_Wp1;
    const ccomplex_t IT_0332 = IT_0003*IT_0331;
    const ccomplex_t IT_0333 = IT_0007*(IT_0330 + IT_0332);
    const ccomplex_t IT_0334 = 1.4142135623731*e_em*IT_0333;
    const ccomplex_t IT_0335 = 0.5*IT_0334;
    const ccomplex_t IT_0336 = IT_0036*IT_0335;
    const ccomplex_t IT_0337 = IT_0040*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = m_N_1*IT_0338;
    const ccomplex_t IT_0340 = 0.5*IT_0317 + -IT_0328 + -IT_0339;
    const ccomplex_t IT_0341 = IT_0100*IT_0324;
    const ccomplex_t IT_0342 = IT_0020*IT_0021*IT_0341;
    const ccomplex_t IT_0343 = IT_0114*IT_0335;
    const ccomplex_t IT_0344 = IT_0038*IT_0040*IT_0343;
    const ccomplex_t IT_0345 = m_N_3*IT_0252;
    const ccomplex_t IT_0346 = IT_0168*IT_0345;
    const ccomplex_t IT_0347 = m_N_1*IT_0231;
    const ccomplex_t IT_0348 = IT_0102*IT_0347;
    const ccomplex_t IT_0349 = m_C_1*IT_0289;
    const ccomplex_t IT_0350 = IT_0146*IT_0349;
    const ccomplex_t IT_0351 = m_N_2*IT_0289;
    const ccomplex_t IT_0352 = IT_0157*IT_0351;
    const ccomplex_t IT_0353 = m_C_1*IT_0252;
    const ccomplex_t IT_0354 = IT_0091*IT_0353;
    const ccomplex_t IT_0355 = m_C_1*IT_0307;
    const ccomplex_t IT_0356 = IT_0197*IT_0355;
    const ccomplex_t IT_0357 = m_N_4*IT_0307;
    const ccomplex_t IT_0358 = IT_0208*IT_0357;
    const ccomplex_t IT_0359 = 0.5*IT_0317 + IT_0342 + IT_0344 + -IT_0346 + 
      -IT_0348 + IT_0350 + -IT_0352 + IT_0354 + IT_0356 + -IT_0358;
    const ccomplex_t IT_0360 = conj(IT_0340) + conj(IT_0359);
    const ccomplex_t IT_0361 = s_14*s_24*IT_0212;
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = s_12 + IT_0362;
    const ccomplex_t IT_0364 = m_C_1*IT_0363;
    const ccomplex_t IT_0365 = (-2)*IT_0364;
    const ccomplex_t IT_0366 = IT_0309*IT_0365;
    const ccomplex_t IT_0367 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0368 = IT_0050*IT_0367;
    const ccomplex_t IT_0369 = 2*IT_0368;
    const ccomplex_t IT_0370 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0371 = IT_0044*IT_0370;
    const ccomplex_t IT_0372 = 1.4142135623731*IT_0371;
    const ccomplex_t IT_0373 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0374 = IT_0050*IT_0373;
    const ccomplex_t IT_0375 = 1.4142135623731*IT_0374;
    const ccomplex_t IT_0376 = IT_0372 + IT_0375;
    const ccomplex_t IT_0377 = IT_0369 + IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*IT_0377;
    const ccomplex_t IT_0379 = 0.5*IT_0378;
    const ccomplex_t IT_0380 = IT_0256*IT_0379;
    const ccomplex_t IT_0381 = (0 + _Complex_I*1)*IT_0380;
    const ccomplex_t IT_0382 = IT_0255*IT_0381;
    const ccomplex_t IT_0383 = pow(s_24, 2);
    const ccomplex_t IT_0384 = IT_0212*IT_0383;
    const ccomplex_t IT_0385 = -IT_0088;
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = s_12*IT_0386;
    const ccomplex_t IT_0388 = m_C_1*m_N_1;
    const ccomplex_t IT_0389 = IT_0363*IT_0388;
    const ccomplex_t IT_0390 = m_N_1*IT_0363;
    const ccomplex_t IT_0391 = m_C_1*IT_0386;
    const ccomplex_t IT_0392 = -IT_0327 + -IT_0338;
    const ccomplex_t IT_0393 = IT_0079*IT_0168;
    const ccomplex_t IT_0394 = IT_0186*IT_0208;
    const ccomplex_t IT_0395 = IT_0056*IT_0102;
    const ccomplex_t IT_0396 = IT_0135*IT_0157;
    const ccomplex_t IT_0397 = IT_0393 + IT_0394 + IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = conj(IT_0392) + conj(IT_0397);
    const ccomplex_t IT_0399 = s_12*s_24;
    const ccomplex_t IT_0400 = s_14*IT_0088;
    const ccomplex_t IT_0401 = -IT_0400;
    const ccomplex_t IT_0402 = IT_0399 + IT_0401;
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = (-2)*IT_0389;
    const ccomplex_t IT_0405 = 0.5*conj(IT_0232);
    const ccomplex_t IT_0406 = s_12*IT_0363;
    const ccomplex_t IT_0407 = (-2)*IT_0406;
    const ccomplex_t IT_0408 = IT_0382 + (-2)*IT_0393 + (-2)*IT_0394 + (-2)
      *IT_0395 + (-2)*IT_0396;
    const ccomplex_t IT_0409 = 0.5*conj(IT_0408);
    const ccomplex_t IT_0410 = (-2)*IT_0387;
    const ccomplex_t IT_0411 = 2*IT_0327 + 2*IT_0338;
    const ccomplex_t IT_0412 = 0.5*conj(IT_0411);
    const ccomplex_t IT_0413 = IT_0386*IT_0388;
    const ccomplex_t IT_0414 = (-2)*IT_0413;
    const ccomplex_t IT_0415 = (-2)*IT_0106 + (-2)*IT_0117 + -IT_0272;
    const ccomplex_t IT_0416 = 0.5*conj(IT_0415);
    const ccomplex_t IT_0417 = 2*IT_0382*(conj(IT_0382)*IT_0387 + conj(IT_0309
      )*IT_0389 + IT_0360*IT_0390 + IT_0211*IT_0391 + 0.5*IT_0398*IT_0403 +
       IT_0404*IT_0405 + IT_0407*IT_0409 + IT_0410*IT_0412 + IT_0414*IT_0416);
    const ccomplex_t IT_0418 = pow(s_14, 2);
    const ccomplex_t IT_0419 = IT_0212*IT_0418;
    const ccomplex_t IT_0420 = -IT_0039;
    const ccomplex_t IT_0421 = IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = IT_0388*IT_0421;
    const ccomplex_t IT_0423 = s_12*IT_0421;
    const ccomplex_t IT_0424 = 0.5*conj(IT_0382);
    const ccomplex_t IT_0425 = 2*IT_0406;
    const ccomplex_t IT_0426 = (-2)*IT_0422;
    const ccomplex_t IT_0427 = 0.5*conj(IT_0309);
    const ccomplex_t IT_0428 = m_N_1*IT_0421;
    const ccomplex_t IT_0429 = (-2)*IT_0428;
    const ccomplex_t IT_0430 = 0.5*IT_0360;
    const ccomplex_t IT_0431 = 0.5*IT_0211;
    const ccomplex_t IT_0432 = (-2)*IT_0390;
    const ccomplex_t IT_0433 = (-2)*IT_0391;
    const ccomplex_t IT_0434 = 0.5*IT_0232 + IT_0253 + IT_0290 + IT_0308;
    const ccomplex_t IT_0435 = -conj(IT_0309);
    const ccomplex_t IT_0436 = conj(IT_0232) + IT_0435;
    const ccomplex_t IT_0437 = s_14*IT_0213;
    const ccomplex_t IT_0438 = IT_0212*IT_0437;
    const ccomplex_t IT_0439 = (-0.25)*IT_0438;
    const ccomplex_t IT_0440 = s_14 + IT_0439;
    const ccomplex_t IT_0441 = m_N_1*IT_0440;
    const ccomplex_t IT_0442 = (-4)*IT_0441;
    const ccomplex_t IT_0443 = s_24*IT_0213;
    const ccomplex_t IT_0444 = IT_0212*IT_0443;
    const ccomplex_t IT_0445 = (-0.25)*IT_0444;
    const ccomplex_t IT_0446 = s_24 + IT_0445;
    const ccomplex_t IT_0447 = m_C_1*IT_0446;
    const ccomplex_t IT_0448 = (-4)*IT_0447;
    const ccomplex_t IT_0449 = 0.5*IT_0444;
    const ccomplex_t IT_0450 = s_24 + IT_0449;
    const ccomplex_t IT_0451 = m_C_1*IT_0450;
    const ccomplex_t IT_0452 = 2*IT_0451;
    const ccomplex_t IT_0453 = 6*IT_0388;
    const ccomplex_t IT_0454 = conj(IT_0062)*IT_0453;
    const ccomplex_t IT_0455 = conj(IT_0309)*IT_0365;
    const ccomplex_t IT_0456 = conj(IT_0408)*IT_0429;
    const ccomplex_t IT_0457 = conj(IT_0411)*IT_0432;
    const ccomplex_t IT_0458 = (-2)*IT_0448;
    const ccomplex_t IT_0459 = 0.5*conj(IT_0233);
    const ccomplex_t IT_0460 = 0.5*IT_0438;
    const ccomplex_t IT_0461 = s_14 + IT_0460;
    const ccomplex_t IT_0462 = m_N_1*IT_0461;
    const ccomplex_t IT_0463 = 2*IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0463;
    const ccomplex_t IT_0465 = 0.5*conj(IT_0397);
    const ccomplex_t IT_0466 = conj(IT_0232)*IT_0364 + conj(IT_0382)*IT_0390 +
       IT_0391*conj(IT_0415) + IT_0221*IT_0430 + conj(IT_0392)*IT_0442 + conj
      (IT_0434)*IT_0448 + conj(IT_0234)*IT_0452 + 0.5*conj(IT_0210)*IT_0453 +
       0.5*IT_0454 + 0.5*IT_0455 + 0.5*IT_0456 + 0.5*IT_0457 + IT_0458*IT_0459 +
       IT_0464*IT_0465;
    const ccomplex_t IT_0467 = 0.5*conj(IT_0434);
    const ccomplex_t IT_0468 = 2*IT_0415*(IT_0211*IT_0390 + IT_0360*IT_0391 +
       IT_0235*IT_0402 + IT_0389*conj(IT_0408) + conj(IT_0411)*IT_0413 + IT_0387
      *conj(IT_0415) + IT_0414*IT_0424 + IT_0405*IT_0425 + IT_0407*IT_0427 +
       IT_0403*IT_0467);
    const ccomplex_t IT_0469 = IT_0213*IT_0388;
    const ccomplex_t IT_0470 = 3*IT_0469;
    const ccomplex_t IT_0471 = pow(m_W, 4);
    const ccomplex_t IT_0472 = s_12*IT_0471;
    const ccomplex_t IT_0473 = IT_0212*IT_0472;
    const ccomplex_t IT_0474 = (-0.5)*IT_0473;
    const ccomplex_t IT_0475 = IT_0214 + IT_0216 + IT_0474;
    const ccomplex_t IT_0476 = (-3)*IT_0469;
    const ccomplex_t IT_0477 = (-2)*IT_0476;
    const ccomplex_t IT_0478 = (-4)*IT_0475;
    const ccomplex_t IT_0479 = (-0.25)*IT_0473;
    const ccomplex_t IT_0480 = IT_0215 + IT_0479;
    const ccomplex_t IT_0481 = (-8)*IT_0480;
    const ccomplex_t IT_0482 = 0.5*conj(IT_0234);
    const ccomplex_t IT_0483 = conj(IT_0309)*IT_0313 + IT_0314*IT_0405 +
       IT_0403*IT_0416 + IT_0360*IT_0448 + IT_0431*IT_0464 + conj(IT_0392)
      *IT_0470 + 2*conj(IT_0233)*IT_0475 + IT_0465*IT_0477 + IT_0467*IT_0478 +
       IT_0481*IT_0482;
    const ccomplex_t IT_0484 = IT_0427*IT_0429;
    const ccomplex_t IT_0485 = (-2)*IT_0452;
    const ccomplex_t IT_0486 = 0.5*conj(IT_0392);
    const ccomplex_t IT_0487 = IT_0485*IT_0486;
    const ccomplex_t IT_0488 = IT_0464*IT_0467;
    const ccomplex_t IT_0489 = (-2)*IT_0442;
    const ccomplex_t IT_0490 = IT_0482*IT_0489;
    const ccomplex_t IT_0491 = IT_0365*conj(IT_0408);
    const ccomplex_t IT_0492 = 0.5*IT_0491;
    const ccomplex_t IT_0493 = conj(IT_0411)*IT_0433;
    const ccomplex_t IT_0494 = 0.5*IT_0493;
    const ccomplex_t IT_0495 = conj(IT_0233)*IT_0463;
    const ccomplex_t IT_0496 = conj(IT_0397)*IT_0448;
    const ccomplex_t IT_0497 = (-2)*IT_0470;
    const ccomplex_t IT_0498 = 8*IT_0480;
    const ccomplex_t IT_0499 = conj(IT_0382)*IT_0403;
    const ccomplex_t IT_0500 = 0.5*IT_0499;
    const ccomplex_t IT_0501 = IT_0313*conj(IT_0408);
    const ccomplex_t IT_0502 = IT_0402*conj(IT_0411);
    const ccomplex_t IT_0503 = conj(IT_0232)*IT_0313;
    const ccomplex_t IT_0504 = IT_0402*conj(IT_0415);
    const ccomplex_t IT_0505 = 2*IT_0233*(IT_0314*IT_0427 + IT_0430*IT_0458 +
       IT_0211*IT_0463 + 2*conj(IT_0434)*IT_0475 + conj(IT_0397)*IT_0476 +
       IT_0459*IT_0478 + IT_0486*IT_0497 + IT_0482*IT_0498 + IT_0503 + IT_0504);
    const ccomplex_t IT_0506 = IT_0211*IT_0222 + IT_0235*IT_0315 + IT_0360
      *IT_0366 + IT_0417 + 2*IT_0408*(IT_0313*IT_0398 + IT_0389*conj(IT_0415) +
       conj(IT_0232)*IT_0422 + conj(IT_0408)*IT_0423 + IT_0407*IT_0424 + IT_0412
      *IT_0425 + IT_0426*IT_0427 + IT_0429*IT_0430 + IT_0365*IT_0431) + 2
      *IT_0411*(conj(IT_0232)*IT_0389 + IT_0398*IT_0402 + IT_0387*conj(IT_0411) 
      + IT_0413*conj(IT_0415) + IT_0410*IT_0424 + IT_0409*IT_0425 + IT_0404
      *IT_0427 + IT_0430*IT_0432 + IT_0431*IT_0433) + 2*IT_0309*(conj(IT_0382)
      *IT_0389 + IT_0404*IT_0412 + IT_0407*IT_0416 + IT_0409*IT_0426 + IT_0429
      *IT_0431 + IT_0313*conj(IT_0434) + -IT_0423*IT_0436) + 2*(IT_0340 +
       IT_0359)*IT_0466 + 2*IT_0232*(IT_0235*IT_0313 + IT_0360*IT_0364 + IT_0389
      *conj(IT_0411) + conj(IT_0408)*IT_0422 + IT_0404*IT_0424 + IT_0416*IT_0425
       + IT_0211*IT_0428 + IT_0423*IT_0436 + IT_0314*IT_0467) + IT_0468 + 2
      *IT_0434*IT_0483 + 2*IT_0062*(conj(IT_0382)*IT_0391 + IT_0390*conj(IT_0415
      ) + conj(IT_0232)*IT_0428 + IT_0430*IT_0453 + IT_0484 + IT_0487 + IT_0488 
      + IT_0490 + IT_0492 + IT_0494 + IT_0495 + IT_0496) + 2*IT_0210*(conj
      (IT_0382)*IT_0391 + IT_0390*conj(IT_0415) + conj(IT_0232)*IT_0428 +
       IT_0221*IT_0431 + IT_0430*IT_0453 + IT_0484 + IT_0487 + IT_0488 + IT_0490
       + IT_0492 + IT_0494 + IT_0495 + IT_0496) + 2*IT_0392*(IT_0360*IT_0442 +
       conj(IT_0434)*IT_0470 + conj(IT_0234)*IT_0476 + IT_0431*IT_0485 + IT_0478
      *IT_0486 + IT_0459*IT_0497 + IT_0465*IT_0498 + IT_0500 + IT_0501 + IT_0502
      ) + 2*IT_0397*(IT_0211*IT_0448 + IT_0430*IT_0464 + conj(IT_0233)*IT_0476 +
       IT_0467*IT_0477 + IT_0465*IT_0478 + IT_0482*IT_0497 + IT_0486*IT_0498 +
       IT_0500 + IT_0501 + IT_0502) + 2*IT_0234*(IT_0314*IT_0427 + IT_0360
      *IT_0452 + conj(IT_0392)*IT_0476 + IT_0467*IT_0481 + IT_0478*IT_0482 +
       IT_0431*IT_0489 + IT_0465*IT_0497 + IT_0459*IT_0498 + IT_0503 + IT_0504) 
      + IT_0505;
    return create_ccomplex_return(IT_0506);
}

