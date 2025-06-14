#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_se_L_to_anti_e_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_se_L_to_anti_e_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0009 + (-0.5)*IT_0012);
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_e, 2);
    const ccomplex_t IT_0017 = pow(m_se_L, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0000*IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = IT_0020*IT_0021*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0038 = IT_0005*IT_0037;
    const ccomplex_t IT_0039 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0040 = IT_0005*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + (-0.5)*IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0036*IT_0043;
    const ccomplex_t IT_0045 = IT_0033*IT_0034*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0050 = IT_0005*IT_0049;
    const ccomplex_t IT_0051 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0052 = IT_0005*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0046*IT_0056*IT_0057;
    const ccomplex_t IT_0059 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0060 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0061 = IT_0005*IT_0060;
    const ccomplex_t IT_0062 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0063 = IT_0005*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0023*IT_0066;
    const ccomplex_t IT_0068 = IT_0021*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0071 = IT_0005*IT_0070;
    const ccomplex_t IT_0072 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + 0.5*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0036*IT_0076;
    const ccomplex_t IT_0078 = IT_0034*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0081 = IT_0005*IT_0080;
    const ccomplex_t IT_0082 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0083 = IT_0005*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 0.5*IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0048*IT_0086;
    const ccomplex_t IT_0088 = IT_0057*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0091 = IT_0005*IT_0090;
    const ccomplex_t IT_0092 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0093 = IT_0005*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + 0.5*IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0007*IT_0096;
    const ccomplex_t IT_0098 = IT_0018*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0069 + -IT_0079 + -IT_0089 + -IT_0099;
    const ccomplex_t IT_0101 = pow(m_W, -2);
    const ccomplex_t IT_0102 = pow(m_W, 2);
    const ccomplex_t IT_0103 = s_34*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = (-0.25)*IT_0104;
    const ccomplex_t IT_0106 = s_34 + IT_0105;
    const ccomplex_t IT_0107 = m_C_2*IT_0106;
    const ccomplex_t IT_0108 = (-4)*IT_0107;
    const ccomplex_t IT_0109 = s_13*IT_0102;
    const ccomplex_t IT_0110 = s_14*s_34;
    const ccomplex_t IT_0111 = (-2)*IT_0110;
    const ccomplex_t IT_0112 = IT_0109 + IT_0111;
    const ccomplex_t IT_0113 = IT_0101*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = s_13 + IT_0114;
    const ccomplex_t IT_0116 = IT_0069 + IT_0079 + IT_0089 + IT_0099;
    const ccomplex_t IT_0117 = m_C_2*IT_0116;
    const ccomplex_t IT_0118 = IT_0059 + IT_0117;
    const ccomplex_t IT_0119 = m_e*m_C_2;
    const ccomplex_t IT_0120 = pow(s_14, 2);
    const ccomplex_t IT_0121 = IT_0101*IT_0120;
    const ccomplex_t IT_0122 = pow(m_C_2, 2);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0121 + IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_13 + IT_0016 + IT_0122 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0131 = IT_0129*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = 2*IT_0069 + 2*IT_0079 + 2*IT_0089 + 2*IT_0099 +
       IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0130*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0127*IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0140;
    const ccomplex_t IT_0142 = cos(theta_W);
    const ccomplex_t IT_0143 = cpow(IT_0142, -1);
    const ccomplex_t IT_0144 = N_B1*e_em;
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = N_W1*e_em;
    const ccomplex_t IT_0148 = IT_0005*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0030*IT_0151;
    const ccomplex_t IT_0153 = IT_0021*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = N_B2*e_em;
    const ccomplex_t IT_0156 = IT_0143*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = N_W2*e_em;
    const ccomplex_t IT_0159 = IT_0005*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + IT_0160);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0043*IT_0162;
    const ccomplex_t IT_0164 = IT_0034*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = N_B3*e_em;
    const ccomplex_t IT_0167 = IT_0143*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = N_W3*e_em;
    const ccomplex_t IT_0170 = IT_0005*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + IT_0171);
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0055*IT_0173;
    const ccomplex_t IT_0175 = IT_0057*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = N_B4*e_em;
    const ccomplex_t IT_0178 = IT_0143*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = N_W4*e_em;
    const ccomplex_t IT_0181 = IT_0005*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0179 + IT_0182);
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0014*IT_0184;
    const ccomplex_t IT_0186 = IT_0018*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0154 + 2*IT_0165 + 2*IT_0176 + 2*IT_0187;
    const ccomplex_t IT_0189 = (-2)*IT_0188;
    const ccomplex_t IT_0190 = IT_0141 + IT_0189;
    const ccomplex_t IT_0191 = 2*conj(IT_0140);
    const ccomplex_t IT_0192 = (-2)*conj(IT_0188);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = s_14*IT_0016;
    const ccomplex_t IT_0195 = s_13*s_34;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = -IT_0133;
    const ccomplex_t IT_0199 = 2*conj(IT_0198);
    const ccomplex_t IT_0200 = s_13*s_14;
    const ccomplex_t IT_0201 = s_34*IT_0122;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = 2*conj(IT_0134);
    const ccomplex_t IT_0205 = s_14*s_34*IT_0101;
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = s_13 + IT_0206;
    const ccomplex_t IT_0208 = m_e*IT_0207;
    const ccomplex_t IT_0209 = IT_0154 + IT_0165 + IT_0176 + IT_0187;
    const ccomplex_t IT_0210 = m_C_2*IT_0209;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0076*IT_0162;
    const ccomplex_t IT_0213 = IT_0033*IT_0034*IT_0212;
    const ccomplex_t IT_0214 = IT_0086*IT_0173;
    const ccomplex_t IT_0215 = IT_0046*IT_0057*IT_0214;
    const ccomplex_t IT_0216 = IT_0096*IT_0184;
    const ccomplex_t IT_0217 = IT_0000*IT_0018*IT_0216;
    const ccomplex_t IT_0218 = IT_0066*IT_0151;
    const ccomplex_t IT_0219 = IT_0020*IT_0021*IT_0218;
    const ccomplex_t IT_0220 = IT_0213 + IT_0215 + IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = IT_0211 + IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0139;
    const ccomplex_t IT_0223 = 2*conj(IT_0139);
    const ccomplex_t IT_0224 = pow(s_34, 2);
    const ccomplex_t IT_0225 = IT_0101*IT_0224;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = IT_0016 + IT_0226;
    const ccomplex_t IT_0228 = m_C_2*IT_0227;
    const ccomplex_t IT_0229 = m_C_2*IT_0207;
    const ccomplex_t IT_0230 = conj(IT_0211) + conj(IT_0220);
    const ccomplex_t IT_0231 = 2*IT_0134;
    const ccomplex_t IT_0232 = 2*IT_0059;
    const ccomplex_t IT_0233 = 2*conj(IT_0059);
    const ccomplex_t IT_0234 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0235 = s_13*IT_0124;
    const ccomplex_t IT_0236 = 2*IT_0188;
    const ccomplex_t IT_0237 = m_e*IT_0124;
    const ccomplex_t IT_0238 = (-2)*IT_0140;
    const ccomplex_t IT_0239 = IT_0236 + IT_0238;
    const ccomplex_t IT_0240 = 2*conj(IT_0188);
    const ccomplex_t IT_0241 = (-2)*conj(IT_0140);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = (-2)*IT_0134;
    const ccomplex_t IT_0244 = (-2)*conj(IT_0134);
    const ccomplex_t IT_0245 = s_13*IT_0227;
    const ccomplex_t IT_0246 = (-2)*IT_0139;
    const ccomplex_t IT_0247 = IT_0119*IT_0207;
    const ccomplex_t IT_0248 = s_13*IT_0207;
    const ccomplex_t IT_0249 = IT_0119*IT_0227;
    const ccomplex_t IT_0250 = (-2)*conj(IT_0139);
    const ccomplex_t IT_0251 = -IT_0154 + -IT_0165 + -IT_0176 + -IT_0187;
    const ccomplex_t IT_0252 = s_14*IT_0102;
    const ccomplex_t IT_0253 = IT_0101*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = s_14 + IT_0254;
    const ccomplex_t IT_0256 = m_e*IT_0255;
    const ccomplex_t IT_0257 = 2*IT_0256;
    const ccomplex_t IT_0258 = IT_0102*IT_0119;
    const ccomplex_t IT_0259 = (-3)*IT_0258;
    const ccomplex_t IT_0260 = pow(m_W, 4);
    const ccomplex_t IT_0261 = s_13*IT_0260;
    const ccomplex_t IT_0262 = IT_0101*IT_0261;
    const ccomplex_t IT_0263 = IT_0109 + IT_0111 + (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0198;
    const ccomplex_t IT_0265 = (-0.5)*IT_0197;
    const ccomplex_t IT_0266 = 0.333333333333333*IT_0100;
    const ccomplex_t IT_0267 = (-0.333333333333333)*IT_0251;
    const ccomplex_t IT_0268 = (-0.166666666666667)*IT_0228;
    const ccomplex_t IT_0269 = 4*IT_0211;
    const ccomplex_t IT_0270 = (-2)*conj(IT_0059)*(IT_0100*IT_0108 + (-2)
      *IT_0115*IT_0118) + IT_0125*(conj(IT_0134)*IT_0190 + IT_0134*IT_0193) +
       IT_0100*(IT_0197*IT_0199 + IT_0203*IT_0204) + IT_0208*(IT_0199*IT_0221 +
       conj(IT_0059)*IT_0222 + IT_0118*IT_0223) + (IT_0117*IT_0199 + IT_0221
      *IT_0223)*IT_0228 + IT_0229*(IT_0117*IT_0204 + IT_0193*IT_0221 + IT_0190
      *IT_0230 + conj(IT_0117)*IT_0231 + conj(IT_0134)*IT_0232 + IT_0134*IT_0233
      ) + IT_0119*IT_0221*IT_0234 + IT_0235*(conj(IT_0140)*IT_0190 + IT_0140
      *IT_0192 + conj(IT_0134)*IT_0231 + conj(IT_0188)*IT_0236) + IT_0237*((conj
      (IT_0059) + conj(IT_0117))*IT_0239 + IT_0118*IT_0242 + IT_0230*IT_0243 +
       IT_0221*IT_0244) + conj(IT_0139)*IT_0245*IT_0246 + conj(IT_0198)*(IT_0228
      *IT_0232 + IT_0239*IT_0247 + IT_0243*IT_0248 + IT_0246*IT_0249) + IT_0247*
      (conj(IT_0134)*IT_0246 + IT_0134*IT_0250) + IT_0248*(conj(IT_0139)*IT_0236
       + IT_0139*IT_0240 + conj(IT_0140)*IT_0246 + IT_0140*IT_0250) + (-2)
      *IT_0198*((-0.5)*IT_0228*IT_0233 + conj(IT_0198)*IT_0245 + (-0.5)*IT_0242
      *IT_0247 + (-0.5)*IT_0244*IT_0248 + (-0.5)*IT_0249*IT_0250) + IT_0251*
      (IT_0203*IT_0242 + IT_0197*IT_0250 + IT_0233*IT_0257) + 2*conj(IT_0251)*
      (IT_0108*IT_0221 + 0.5*IT_0203*IT_0239 + (IT_0117 + 0.5*IT_0232)*IT_0257 +
       -IT_0100*IT_0259 + (-2)*IT_0251*IT_0263 + -IT_0246*IT_0265) + (-2)*conj
      (IT_0100)*(IT_0108*IT_0118 + (-0.5)*IT_0203*IT_0231 + IT_0221*IT_0257 +
       IT_0251*IT_0259 + 2*IT_0100*IT_0263 + IT_0264*IT_0265) + 4*conj(IT_0117)*
      (IT_0115*IT_0118 + (-1.5)*IT_0119*IT_0221 + 0.25*IT_0208*IT_0222 + (-1.5)
      *IT_0108*IT_0266 + (-1.5)*IT_0257*IT_0267 + (-1.5)*IT_0264*IT_0268) + (-6)
      *IT_0230*(IT_0118*IT_0119 + (-0.166666666666667)*IT_0208*IT_0264 + IT_0257
      *IT_0266 + IT_0108*IT_0267 + IT_0222*IT_0268 + (-0.666666666666667)
      *IT_0115*(IT_0220 + 0.25*IT_0269));
    return create_ccomplex_return(IT_0270);
}

