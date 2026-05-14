#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_snu_tau_to_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_snu_tau_to_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_01);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = N_d2*e_em*m_tau*IT_0012*U_stau_01;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0008 + 0.5*IT_0015);
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = pow(m_N_2, 2);
    const ccomplex_t IT_0019 = pow(m_tau, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0019 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0025 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0024 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0001*IT_0010*IT_0012;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0029 + 0.5*IT_0032);
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = IT_0027*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0024 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0001*IT_0010*IT_0012;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0042 = IT_0001*IT_0041;
    const ccomplex_t IT_0043 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0044 = IT_0001*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.5*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0040*IT_0047;
    const ccomplex_t IT_0049 = IT_0038*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_00);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0054 = IT_0005*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_d2*e_em*m_tau*IT_0012*U_stau_00;
    const ccomplex_t IT_0057 = IT_0011*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 0.5*IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0019 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0023 + IT_0037 + IT_0050 + (-0.5)
      *IT_0064;
    const ccomplex_t IT_0066 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0067 = IT_0005*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0070 = IT_0001*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d2)*e_em*m_tau*IT_0012*U_stau_10;
    const ccomplex_t IT_0073 = IT_0011*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + IT_0071 + 
      -IT_0074);
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0061*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0052*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0081 = IT_0005*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_d2)*e_em*m_tau*IT_0012*U_stau_11;
    const ccomplex_t IT_0087 = IT_0011*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0082 + IT_0085 + 
      -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0020*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0003*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0079 + 0.5*IT_0093;
    const ccomplex_t IT_0095 = m_N_2*m_tau;
    const ccomplex_t IT_0096 = pow(m_W, -2);
    const ccomplex_t IT_0097 = pow(s_14, 2);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = -IT_0018;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = IT_0095*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0037 + (-0.5)*IT_0050;
    const ccomplex_t IT_0103 = s_13*s_14;
    const ccomplex_t IT_0104 = s_34*IT_0018;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B2)*e_em;
    const ccomplex_t IT_0110 = IT_0005*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W2)*e_em;
    const ccomplex_t IT_0113 = IT_0001*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + -IT_0114);
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0108*IT_0116;
    const ccomplex_t IT_0118 = cpow(s_12 + 0.5*IT_0018 + 0.5*IT_0024 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = pow(m_W, 2);
    const ccomplex_t IT_0123 = IT_0095*IT_0122;
    const ccomplex_t IT_0124 = 3*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0079 + (-0.5)*IT_0093 + -IT_0120;
    const ccomplex_t IT_0126 = s_13*s_34;
    const ccomplex_t IT_0127 = s_14*IT_0019;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = s_13*IT_0122;
    const ccomplex_t IT_0131 = pow(m_W, 4);
    const ccomplex_t IT_0132 = s_13*IT_0131;
    const ccomplex_t IT_0133 = IT_0096*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = s_14*s_34;
    const ccomplex_t IT_0136 = (-2)*IT_0135;
    const ccomplex_t IT_0137 = IT_0130 + IT_0134 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0139 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0140 = IT_0001*IT_0139;
    const ccomplex_t IT_0141 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + (-0.5)*IT_0143);
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = IT_0138*IT_0145;
    const ccomplex_t IT_0147 = IT_0025*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0150 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0151 = IT_0001*IT_0150;
    const ccomplex_t IT_0152 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0153 = IT_0001*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + (-0.5)*IT_0154);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0149*IT_0156;
    const ccomplex_t IT_0158 = IT_0038*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (-0.5)*IT_0148 + (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = (-0.25)*IT_0133;
    const ccomplex_t IT_0162 = IT_0135 + IT_0161;
    const ccomplex_t IT_0163 = IT_0148 + IT_0159;
    const ccomplex_t IT_0164 = (-2)*IT_0094;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0163 + IT_0165;
    const ccomplex_t IT_0167 = s_14*IT_0122;
    const ccomplex_t IT_0168 = IT_0096*IT_0167;
    const ccomplex_t IT_0169 = (-0.25)*IT_0168;
    const ccomplex_t IT_0170 = s_14 + IT_0169;
    const ccomplex_t IT_0171 = m_tau*IT_0170;
    const ccomplex_t IT_0172 = (-4)*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0174 = IT_0040*IT_0156;
    const ccomplex_t IT_0175 = IT_0038*IT_0173*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0177 = IT_0027*IT_0145;
    const ccomplex_t IT_0178 = IT_0025*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0175 + 0.5*IT_0178;
    const ccomplex_t IT_0180 = m_N_2*IT_0037;
    const ccomplex_t IT_0181 = m_N_2*IT_0050;
    const ccomplex_t IT_0182 = m_tau*IT_0120;
    const ccomplex_t IT_0183 = (-0.5)*IT_0180 + (-0.5)*IT_0181 + 0.5*IT_0182;
    const ccomplex_t IT_0184 = IT_0179 + IT_0183;
    const ccomplex_t IT_0185 = s_34*IT_0122;
    const ccomplex_t IT_0186 = IT_0096*IT_0185;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = s_34 + IT_0187;
    const ccomplex_t IT_0189 = m_N_2*IT_0188;
    const ccomplex_t IT_0190 = 2*IT_0189;
    const ccomplex_t IT_0191 = m_N_2*IT_0163;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0034*IT_0138;
    const ccomplex_t IT_0194 = IT_0025*IT_0176*IT_0193;
    const ccomplex_t IT_0195 = IT_0047*IT_0149;
    const ccomplex_t IT_0196 = IT_0038*IT_0173*IT_0195;
    const ccomplex_t IT_0197 = (-0.5)*IT_0194 + (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0192 + IT_0197;
    const ccomplex_t IT_0199 = (-3)*IT_0123;
    const ccomplex_t IT_0200 = (-0.25)*IT_0186;
    const ccomplex_t IT_0201 = s_34 + IT_0200;
    const ccomplex_t IT_0202 = m_N_2*IT_0201;
    const ccomplex_t IT_0203 = (-4)*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0168;
    const ccomplex_t IT_0205 = s_14 + IT_0204;
    const ccomplex_t IT_0206 = m_tau*IT_0205;
    const ccomplex_t IT_0207 = 2*IT_0206;
    const ccomplex_t IT_0208 = (-2)*IT_0179;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = IT_0183 + IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0023 + 0.5*IT_0064;
    const ccomplex_t IT_0212 = s_14*s_34*IT_0096;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = s_13 + IT_0213;
    const ccomplex_t IT_0215 = IT_0095*IT_0214;
    const ccomplex_t IT_0216 = pow(s_34, 2);
    const ccomplex_t IT_0217 = IT_0096*IT_0216;
    const ccomplex_t IT_0218 = -IT_0019;
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = s_13*IT_0219;
    const ccomplex_t IT_0221 = IT_0095*IT_0219;
    const ccomplex_t IT_0222 = m_N_2*IT_0219;
    const ccomplex_t IT_0223 = s_13*IT_0214;
    const ccomplex_t IT_0224 = (-2)*IT_0065;
    const ccomplex_t IT_0225 = IT_0130 + IT_0136;
    const ccomplex_t IT_0226 = IT_0096*IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = s_13 + IT_0227;
    const ccomplex_t IT_0229 = m_N_2*IT_0214;
    const ccomplex_t IT_0230 = (-2)*conj(IT_0065);
    const ccomplex_t IT_0231 = s_13*IT_0100;
    const ccomplex_t IT_0232 = (-2)*IT_0163;
    const ccomplex_t IT_0233 = (-2)*conj(IT_0163);
    const ccomplex_t IT_0234 = (-2)*IT_0211;
    const ccomplex_t IT_0235 = 2*conj(IT_0102)*(IT_0065*IT_0106 + IT_0121
      *IT_0124 + (-2)*IT_0102*IT_0137 + -IT_0160*IT_0199 + IT_0184*IT_0203 +
       IT_0198*IT_0207 + 0.5*IT_0129*IT_0234);
    const ccomplex_t IT_0236 = m_tau*IT_0214;
    const ccomplex_t IT_0237 = (-2)*conj(IT_0211);
    const ccomplex_t IT_0238 = (-2)*conj(IT_0179);
    const ccomplex_t IT_0239 = m_tau*IT_0100;
    const ccomplex_t IT_0240 = IT_0121 + IT_0160;
    const ccomplex_t IT_0241 = 0.5*IT_0215;
    const ccomplex_t IT_0242 = 4*IT_0192;
    const ccomplex_t IT_0243 = 0.333333333333333*IT_0160;
    const ccomplex_t IT_0244 = (-0.333333333333333)*IT_0102;
    const ccomplex_t IT_0245 = (-0.333333333333333)*IT_0121;
    const ccomplex_t IT_0246 = (-3)*IT_0095*IT_0184 + IT_0125*IT_0222 +
       IT_0166*IT_0229 + 0.5*IT_0234*IT_0236 + IT_0065*IT_0239 + 2*IT_0228*
      (IT_0197 + 0.25*IT_0242) + (-3)*IT_0203*IT_0243 + (-3)*IT_0207*IT_0244 + (
      -3)*IT_0190*IT_0245;
    const ccomplex_t IT_0247 = (-6)*conj(IT_0179);
    const ccomplex_t IT_0248 = IT_0095*IT_0198;
    const ccomplex_t IT_0249 = IT_0211*IT_0222 + 2*IT_0184*IT_0228 + 0.5
      *IT_0224*IT_0229 + -IT_0125*IT_0236 + (IT_0094 + 0.5*IT_0232)*IT_0239 + (
      -3)*IT_0207*IT_0243 + (-3)*IT_0203*IT_0244 + (-3)*IT_0172*IT_0245 + (-3)
      *IT_0248;
    const ccomplex_t IT_0250 = (-2)*conj(IT_0094);
    const ccomplex_t IT_0251 = IT_0198*IT_0229;
    const ccomplex_t IT_0252 = IT_0106*IT_0240;
    const ccomplex_t IT_0253 = 2*conj(IT_0065)*(IT_0094*IT_0101 + IT_0102
      *IT_0106) + 2*conj(IT_0121)*(IT_0102*IT_0124 + IT_0125*IT_0129 + (-2)
      *IT_0121*IT_0137 + 4*IT_0160*IT_0162 + IT_0106*IT_0166 + IT_0172*IT_0184 +
       IT_0190*IT_0198) + 2*conj(IT_0160)*(IT_0125*IT_0129 + (-2)*IT_0137
      *IT_0160 + 4*IT_0121*IT_0162 + IT_0106*IT_0166 + -IT_0102*IT_0199 + 
      -IT_0198*IT_0203 + -IT_0207*IT_0210) + 2*conj(IT_0211)*(IT_0163*IT_0215 +
       IT_0211*IT_0220 + IT_0125*IT_0221 + IT_0184*IT_0222 + 0.5*IT_0223*IT_0224
      ) + 2*conj(IT_0179)*(IT_0121*IT_0172 + IT_0102*IT_0203 + IT_0211*IT_0222 +
       2*IT_0184*IT_0228 + 0.5*IT_0224*IT_0229) + (IT_0101*IT_0163 + IT_0125
      *IT_0215 + IT_0211*IT_0223 + IT_0184*IT_0229)*IT_0230 + 2*IT_0231*(IT_0065
      *conj(IT_0065) + 0.5*conj(IT_0094)*IT_0232 + IT_0094*(conj(IT_0094) + 0.5
      *IT_0233)) + IT_0235 + (IT_0102*IT_0129 + IT_0094*IT_0215 + IT_0198
      *IT_0236)*IT_0237 + (IT_0160*IT_0207 + IT_0125*IT_0236)*IT_0238 + 2*
      (IT_0094*conj(IT_0179) + conj(IT_0094)*IT_0184 + conj(IT_0065)*IT_0198 +
       0.5*IT_0183*IT_0233 + 0.5*IT_0163*IT_0238)*IT_0239 + 2*conj(IT_0125)*
      (IT_0125*IT_0220 + IT_0211*IT_0221 + IT_0198*IT_0222 + IT_0166*IT_0223 + 
      -IT_0210*IT_0236 + IT_0129*IT_0240 + IT_0224*IT_0241) + 2*conj(IT_0094)*
      (IT_0065*IT_0101 + IT_0234*IT_0241) + 2*(conj(IT_0192) + conj(IT_0197))
      *IT_0246 + IT_0247*IT_0248 + 2*conj(IT_0183)*IT_0249 + IT_0250*(IT_0125
      *IT_0223 + IT_0251 + IT_0252) + 2*conj(IT_0163)*(IT_0211*IT_0215 + IT_0125
      *IT_0223 + 0.5*IT_0101*IT_0224 + IT_0163*IT_0231 + 0.5*IT_0208*IT_0239 +
       IT_0251 + IT_0252);
    return create_ccomplex_return(IT_0253);
}

