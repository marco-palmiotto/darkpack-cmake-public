#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_mu_mu(
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
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = e_em*m_mu*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*(IT_0009 + IT_0011);
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_C_1, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0010;
    const ccomplex_t IT_0021 = IT_0003*IT_0008;
    const ccomplex_t IT_0022 = IT_0005*(IT_0020 + -IT_0021);
    const ccomplex_t IT_0023 = 1.4142135623731*e_em*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = e_em*m_mu*IT_0000*IT_0005;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + (-2)*IT_0016 + IT_0028 + 
      -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0019 + -IT_0031;
    const ccomplex_t IT_0033 = pow(m_mu, 2);
    const ccomplex_t IT_0034 = s_34 + IT_0033;
    const ccomplex_t IT_0035 = s_12*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em*m_mu
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0016 + IT_0033 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = m_mu*m_C_1;
    const ccomplex_t IT_0046 = s_13 + s_14;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = V_u1*U_Wm1;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = U_d1*V_Wp1;
    const ccomplex_t IT_0051 = IT_0003*IT_0050;
    const ccomplex_t IT_0052 = IT_0005*(IT_0049 + IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0007*IT_0054;
    const ccomplex_t IT_0056 = IT_0017*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0001*IT_0050;
    const ccomplex_t IT_0059 = IT_0003*IT_0048;
    const ccomplex_t IT_0060 = IT_0005*(IT_0058 + -IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0026*IT_0062;
    const ccomplex_t IT_0064 = IT_0029*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0057 + -IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0039*IT_0068;
    const ccomplex_t IT_0070 = IT_0041*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0016*IT_0034;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0075 = IT_0068*IT_0074;
    const ccomplex_t IT_0076 = IT_0041*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = s_23 + s_24;
    const ccomplex_t IT_0080 = IT_0045*IT_0079;
    const ccomplex_t IT_0081 = sin(alpha);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0005*IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = cos(alpha);
    const ccomplex_t IT_0085 = IT_0008*IT_0084;
    const ccomplex_t IT_0086 = IT_0010*IT_0081;
    const ccomplex_t IT_0087 = IT_0005*(IT_0085 + -IT_0086);
    const ccomplex_t IT_0088 = 1.4142135623731*e_em*IT_0087;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0083*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0005*IT_0084;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0010*IT_0084;
    const ccomplex_t IT_0098 = IT_0008*IT_0081;
    const ccomplex_t IT_0099 = IT_0005*(IT_0097 + IT_0098);
    const ccomplex_t IT_0100 = 1.4142135623731*e_em*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0096*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = -IT_0094 + -IT_0106;
    const ccomplex_t IT_0108 = -s_34;
    const ccomplex_t IT_0109 = IT_0033 + IT_0108;
    const ccomplex_t IT_0110 = s_12*IT_0109;
    const ccomplex_t IT_0111 = IT_0048*IT_0084;
    const ccomplex_t IT_0112 = IT_0050*IT_0081;
    const ccomplex_t IT_0113 = IT_0005*(IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = 1.4142135623731*e_em*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0083*IT_0116;
    const ccomplex_t IT_0118 = IT_0092*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0050*IT_0084;
    const ccomplex_t IT_0121 = IT_0048*IT_0081;
    const ccomplex_t IT_0122 = IT_0005*(IT_0120 + IT_0121);
    const ccomplex_t IT_0123 = 1.4142135623731*e_em*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0096*IT_0125;
    const ccomplex_t IT_0127 = IT_0104*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = -IT_0119 + -IT_0128;
    const ccomplex_t IT_0130 = IT_0016*IT_0109;
    const ccomplex_t IT_0131 = s_12*IT_0033;
    const ccomplex_t IT_0132 = s_14*s_23;
    const ccomplex_t IT_0133 = s_12*s_34;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = s_13*s_24;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0131 + IT_0132 + IT_0134 + IT_0136;
    const ccomplex_t IT_0138 = cos(theta_W);
    const ccomplex_t IT_0139 = cpow(IT_0138, -1);
    const ccomplex_t IT_0140 = IT_0004*IT_0139;
    const ccomplex_t IT_0141 = e_em*IT_0140;
    const ccomplex_t IT_0142 = IT_0005*IT_0138;
    const ccomplex_t IT_0143 = e_em*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0141 + -IT_0143);
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0147 = IT_0140*IT_0146;
    const ccomplex_t IT_0148 = IT_0142*IT_0146;
    const ccomplex_t IT_0149 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0150 = IT_0142*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + -IT_0148 + (-2)
      *IT_0150);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0145*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = -s_14;
    const ccomplex_t IT_0159 = s_13 + IT_0158;
    const ccomplex_t IT_0160 = IT_0045*IT_0159;
    const ccomplex_t IT_0161 = (-4)*IT_0160;
    const ccomplex_t IT_0162 = -s_24;
    const ccomplex_t IT_0163 = s_23 + IT_0162;
    const ccomplex_t IT_0164 = IT_0045*IT_0163;
    const ccomplex_t IT_0165 = 4*IT_0164;
    const ccomplex_t IT_0166 = m_mu*s_13*m_C_1;
    const ccomplex_t IT_0167 = 8*IT_0072;
    const ccomplex_t IT_0168 = m_mu*s_24*m_C_1;
    const ccomplex_t IT_0169 = IT_0037*IT_0074;
    const ccomplex_t IT_0170 = IT_0041*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = 8*IT_0172;
    const ccomplex_t IT_0174 = s_34*IT_0016;
    const ccomplex_t IT_0175 = IT_0016*IT_0033;
    const ccomplex_t IT_0176 = 2*IT_0175;
    const ccomplex_t IT_0177 = 2*IT_0135;
    const ccomplex_t IT_0178 = IT_0131 + IT_0174 + IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0183 = IT_0181*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = 8*IT_0185;
    const ccomplex_t IT_0187 = (-8)*IT_0078;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0139;
    const ccomplex_t IT_0189 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0190 = IT_0140*IT_0189;
    const ccomplex_t IT_0191 = IT_0142*IT_0189;
    const ccomplex_t IT_0192 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0193 = IT_0142*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0190 + -IT_0191 + (-2)
      *IT_0193);
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = IT_0188*IT_0195;
    const ccomplex_t IT_0197 = IT_0154*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = (-8)*IT_0044;
    const ccomplex_t IT_0200 = (-4)*IT_0032;
    const ccomplex_t IT_0201 = IT_0047*IT_0066 + 0.25*IT_0129*IT_0161 + 0.25
      *IT_0107*IT_0165 + 0.25*IT_0166*IT_0167 + 0.25*IT_0168*IT_0173 + 0.25
      *IT_0178*IT_0186 + 0.25*IT_0135*(16*IT_0157 + IT_0187) + 0.25*IT_0175*(16
      *IT_0198 + IT_0199) + 0.25*IT_0080*IT_0200;
    const ccomplex_t IT_0202 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0203 = s_13 + s_24;
    const ccomplex_t IT_0204 = IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = IT_0045*IT_0204;
    const ccomplex_t IT_0206 = 8*conj(IT_0072);
    const ccomplex_t IT_0207 = 8*conj(IT_0172);
    const ccomplex_t IT_0208 = 4*IT_0078;
    const ccomplex_t IT_0209 = 4*IT_0044;
    const ccomplex_t IT_0210 = (-4)*conj(IT_0078);
    const ccomplex_t IT_0211 = (-4)*conj(IT_0044);
    const ccomplex_t IT_0212 = -s_14 + -s_23;
    const ccomplex_t IT_0213 = IT_0203 + IT_0212;
    const ccomplex_t IT_0214 = IT_0045*IT_0213;
    const ccomplex_t IT_0215 = (-8)*IT_0214;
    const ccomplex_t IT_0216 = -IT_0129;
    const ccomplex_t IT_0217 = -IT_0107;
    const ccomplex_t IT_0218 = (32*IT_0131 + 32*IT_0132 + 32*IT_0135 + 32
      *IT_0174 + 64*IT_0175)*IT_0185 + 8*IT_0205*(IT_0072 + IT_0172 + (-0.125)
      *IT_0206 + (-0.125)*IT_0207) + -IT_0178*(IT_0208 + IT_0209 + IT_0210 +
       IT_0211) + -IT_0215*(conj(IT_0107) + conj(IT_0129) + IT_0216 + IT_0217);
    const ccomplex_t IT_0219 = (-4)*conj(IT_0072);
    const ccomplex_t IT_0220 = (-2)*conj(IT_0066);
    const ccomplex_t IT_0221 = 2*IT_0164;
    const ccomplex_t IT_0222 = -conj(IT_0129);
    const ccomplex_t IT_0223 = (-2)*IT_0160;
    const ccomplex_t IT_0224 = -conj(IT_0107);
    const ccomplex_t IT_0225 = -IT_0132;
    const ccomplex_t IT_0226 = IT_0131 + IT_0133 + IT_0135 + IT_0225;
    const ccomplex_t IT_0227 = -IT_0032;
    const ccomplex_t IT_0228 = IT_0145*IT_0195;
    const ccomplex_t IT_0229 = IT_0154*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = (-4)*IT_0044;
    const ccomplex_t IT_0232 = 0.125*IT_0231;
    const ccomplex_t IT_0233 = IT_0198 + IT_0232;
    const ccomplex_t IT_0234 = 8*IT_0157;
    const ccomplex_t IT_0235 = (-4)*IT_0078;
    const ccomplex_t IT_0236 = IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = m_mu*s_23*m_C_1;
    const ccomplex_t IT_0238 = 4*IT_0072;
    const ccomplex_t IT_0239 = m_mu*s_14*m_C_1;
    const ccomplex_t IT_0240 = 4*IT_0172;
    const ccomplex_t IT_0241 = 2*IT_0132;
    const ccomplex_t IT_0242 = IT_0131 + IT_0174 + IT_0176 + IT_0241;
    const ccomplex_t IT_0243 = IT_0152*IT_0188;
    const ccomplex_t IT_0244 = IT_0154*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = 16*IT_0246;
    const ccomplex_t IT_0248 = 0.125*IT_0107*IT_0161;
    const ccomplex_t IT_0249 = 4*IT_0032;
    const ccomplex_t IT_0250 = 0.125*IT_0047;
    const ccomplex_t IT_0251 = (-0.25)*IT_0066*IT_0080 + 0.0625*IT_0129
      *IT_0165 + IT_0132*IT_0230 + 0.5*IT_0131*IT_0233 + 0.0625*IT_0174*IT_0236 
      + 0.0625*IT_0237*IT_0238 + 0.0625*IT_0239*IT_0240 + 0.0625*IT_0186*IT_0242
       + 0.0625*IT_0175*IT_0247 + 0.5*IT_0248 + 0.5*IT_0249*IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0066*IT_0080 + 0.125*IT_0129*IT_0165 +
       0.125*IT_0167*IT_0168 + 0.125*IT_0166*IT_0173 + 0.125*IT_0178*IT_0186 + 2
      *IT_0175*(IT_0157 + 0.0625*IT_0187) + 2*IT_0135*(IT_0198 + 0.0625*IT_0199)
       + IT_0131*IT_0230 + IT_0174*IT_0246 + IT_0248 + IT_0200*IT_0250;
    const ccomplex_t IT_0253 = 0.5*IT_0135;
    const ccomplex_t IT_0254 = 0.5*IT_0166;
    const ccomplex_t IT_0255 = 0.5*IT_0168;
    const ccomplex_t IT_0256 = 0.5*IT_0175;
    const ccomplex_t IT_0257 = 4*conj(IT_0072);
    const ccomplex_t IT_0258 = (-4)*IT_0072;
    const ccomplex_t IT_0259 = (-0.5)*IT_0135;
    const ccomplex_t IT_0260 = IT_0047*IT_0066;
    const ccomplex_t IT_0261 = 4*conj(IT_0078)*(IT_0044*IT_0175 + 0.25*IT_0217
      *IT_0221 + 0.25*IT_0216*IT_0223 + (-0.5)*IT_0080*IT_0227 + 0.25*IT_0166
      *IT_0258 + (-0.5)*IT_0208*IT_0259 + (-0.5)*IT_0260);
    const ccomplex_t IT_0262 = 0.0625*IT_0129*IT_0161 + 0.0625*IT_0107*IT_0165
       + IT_0175*IT_0230 + 0.5*IT_0174*IT_0233 + 0.0625*IT_0131*IT_0236 + 0.0625
      *IT_0238*IT_0239 + 0.0625*IT_0237*IT_0240 + 0.0625*IT_0186*IT_0242 +
       0.0625*IT_0132*IT_0247 + 0.0625*IT_0080*IT_0249 + (-0.25)*IT_0260;
    const ccomplex_t IT_0263 = conj(IT_0129)*IT_0165;
    const ccomplex_t IT_0264 = conj(IT_0107)*IT_0161 + IT_0263;
    const ccomplex_t IT_0265 = 8*conj(IT_0157);
    const ccomplex_t IT_0266 = IT_0210 + IT_0265;
    const ccomplex_t IT_0267 = 4*conj(IT_0032);
    const ccomplex_t IT_0268 = 4*conj(IT_0172);
    const ccomplex_t IT_0269 = (-4)*conj(IT_0066);
    const ccomplex_t IT_0270 = (-8)*IT_0185;
    const ccomplex_t IT_0271 = IT_0242*IT_0270;
    const ccomplex_t IT_0272 = conj(IT_0107)*IT_0165;
    const ccomplex_t IT_0273 = conj(IT_0129)*IT_0161 + IT_0272;
    const ccomplex_t IT_0274 = 4*conj(IT_0066);
    const ccomplex_t IT_0275 = (-8)*conj(IT_0078);
    const ccomplex_t IT_0276 = (-8)*conj(IT_0044);
    const ccomplex_t IT_0277 = (-4)*conj(IT_0032);
    const ccomplex_t IT_0278 = IT_0178*IT_0270;
    const ccomplex_t IT_0279 = 8*conj(IT_0032)*(IT_0032*IT_0035 + 0.25*IT_0044
      *IT_0047 + -(IT_0066 + 0.25*IT_0072)*IT_0073 + 0.25*IT_0078*IT_0080) + (-8
      )*conj(IT_0107)*(IT_0107*IT_0110 + 0.25*(IT_0072 + (-4)*IT_0129)*IT_0130) 
      + (-8)*conj(IT_0129)*(IT_0110*IT_0129 + -IT_0107*IT_0130 + (-0.25)*IT_0072
      *IT_0137) + 4*conj(IT_0157)*IT_0201 + -IT_0185*IT_0218 + IT_0078*(IT_0166
      *IT_0219 + IT_0047*IT_0220 + IT_0222*IT_0223 + IT_0221*IT_0224) + IT_0044*
      (IT_0168*IT_0219 + IT_0080*IT_0220 + IT_0221*IT_0222 + IT_0223*IT_0224) +
       8*conj(IT_0066)*(IT_0035*IT_0066 + 0.25*IT_0072*IT_0226 + IT_0073*IT_0227
      ) + 16*conj(IT_0230)*IT_0251 + 8*conj(IT_0198)*IT_0252 + 2*conj(IT_0172)*
      (IT_0066*IT_0073 + IT_0107*IT_0137 + IT_0130*IT_0216 + IT_0226*IT_0227 +
       IT_0240*IT_0253 + IT_0231*IT_0254 + IT_0235*IT_0255 + IT_0238*IT_0256) +
       2*IT_0172*(conj(IT_0066)*IT_0073 + -conj(IT_0129)*IT_0130 + conj(IT_0107)
      *IT_0137 + -conj(IT_0032)*IT_0226 + IT_0211*IT_0254 + IT_0210*IT_0255 +
       IT_0256*IT_0257) + 2*conj(IT_0044)*(IT_0032*IT_0047 + -IT_0066*IT_0080 +
       0.5*IT_0216*IT_0221 + 0.5*IT_0217*IT_0223 + IT_0209*IT_0253 + IT_0208
      *IT_0256 + IT_0255*IT_0258) + (-2)*conj(IT_0072)*(IT_0032*IT_0073 +
       IT_0107*IT_0130 + IT_0137*IT_0216 + -IT_0066*IT_0226 + IT_0238*IT_0259) +
       IT_0261 + 16*conj(IT_0246)*IT_0262 + IT_0230*(IT_0131*IT_0211 + IT_0237
      *IT_0257 + IT_0264 + IT_0174*IT_0266 + IT_0047*IT_0267 + IT_0239*IT_0268 +
       IT_0080*IT_0269 + IT_0271) + IT_0246*(IT_0174*IT_0211 + IT_0239*IT_0257 +
       IT_0131*IT_0266 + IT_0080*IT_0267 + IT_0237*IT_0268 + IT_0047*IT_0269 +
       IT_0271 + IT_0273) + IT_0198*(IT_0168*IT_0206 + IT_0166*IT_0207 + IT_0264
       + IT_0080*IT_0274 + IT_0175*IT_0275 + IT_0135*IT_0276 + IT_0047*IT_0277 +
       IT_0278) + IT_0157*(IT_0166*IT_0206 + IT_0168*IT_0207 + IT_0273 + IT_0047
      *IT_0274 + IT_0135*IT_0275 + IT_0175*IT_0276 + IT_0080*IT_0277 + IT_0278);
    return create_ccomplex_return(IT_0279);
}

