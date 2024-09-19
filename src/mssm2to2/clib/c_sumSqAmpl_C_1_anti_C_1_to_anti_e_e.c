#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
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
    const creal_t Gamma_nuel = param->Gamma_nuel;
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
    const ccomplex_t IT_0006 = m_e*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
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
    const ccomplex_t IT_0025 = m_e*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + (-2)*IT_0016 + IT_0028 + 
      -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0019 + -IT_0031;
    const ccomplex_t IT_0033 = pow(m_e, 2);
    const ccomplex_t IT_0034 = s_34 + IT_0033;
    const ccomplex_t IT_0035 = s_12*IT_0034;
    const ccomplex_t IT_0036 = V_u1*U_Wm1;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = U_d1*V_Wp1;
    const ccomplex_t IT_0039 = IT_0003*IT_0038;
    const ccomplex_t IT_0040 = IT_0005*(IT_0037 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0001*IT_0038;
    const ccomplex_t IT_0047 = IT_0003*IT_0036;
    const ccomplex_t IT_0048 = IT_0005*(IT_0046 + -IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0026*IT_0050;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0045 + -IT_0053;
    const ccomplex_t IT_0055 = IT_0016*IT_0034;
    const ccomplex_t IT_0056 = cos(theta_W);
    const ccomplex_t IT_0057 = cpow(IT_0056, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0057;
    const ccomplex_t IT_0059 = IT_0004*IT_0057;
    const ccomplex_t IT_0060 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = IT_0005*IT_0056;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0061 + -IT_0063 + (-2)
      *IT_0065);
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0058*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0016*IT_0033;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0016 + IT_0033 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = (-8)*IT_0080;
    const ccomplex_t IT_0082 = m_e*m_C_1;
    const ccomplex_t IT_0083 = s_13 + s_14;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (-4)*IT_0032;
    const ccomplex_t IT_0086 = e_em*IT_0059;
    const ccomplex_t IT_0087 = e_em*IT_0062;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0091 = IT_0059*IT_0090;
    const ccomplex_t IT_0092 = IT_0062*IT_0090;
    const ccomplex_t IT_0093 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0094 = IT_0062*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + -IT_0092 + (-2)
      *IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0089*IT_0096;
    const ccomplex_t IT_0098 = IT_0069*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = s_13*s_24;
    const ccomplex_t IT_0102 = s_23 + s_24;
    const ccomplex_t IT_0103 = IT_0082*IT_0102;
    const ccomplex_t IT_0104 = IT_0058*IT_0096;
    const ccomplex_t IT_0105 = IT_0069*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = sin(alpha);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = cos(alpha);
    const ccomplex_t IT_0112 = IT_0036*IT_0111;
    const ccomplex_t IT_0113 = IT_0038*IT_0108;
    const ccomplex_t IT_0114 = IT_0005*(IT_0112 + -IT_0113);
    const ccomplex_t IT_0115 = 1.4142135623731*e_em*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0110*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0111;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0038*IT_0111;
    const ccomplex_t IT_0125 = IT_0036*IT_0108;
    const ccomplex_t IT_0126 = IT_0005*(IT_0124 + IT_0125);
    const ccomplex_t IT_0127 = 1.4142135623731*e_em*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = IT_0123*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0121 + -IT_0133;
    const ccomplex_t IT_0135 = -s_14;
    const ccomplex_t IT_0136 = s_13 + IT_0135;
    const ccomplex_t IT_0137 = IT_0082*IT_0136;
    const ccomplex_t IT_0138 = (-4)*IT_0137;
    const ccomplex_t IT_0139 = IT_0010*IT_0111;
    const ccomplex_t IT_0140 = IT_0008*IT_0108;
    const ccomplex_t IT_0141 = IT_0005*(IT_0139 + IT_0140);
    const ccomplex_t IT_0142 = 1.4142135623731*e_em*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = IT_0123*IT_0144;
    const ccomplex_t IT_0146 = IT_0131*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0008*IT_0111;
    const ccomplex_t IT_0149 = IT_0010*IT_0108;
    const ccomplex_t IT_0150 = IT_0005*(IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = 1.4142135623731*e_em*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0110*IT_0153;
    const ccomplex_t IT_0155 = IT_0119*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0147 + -IT_0156;
    const ccomplex_t IT_0158 = -s_24;
    const ccomplex_t IT_0159 = s_23 + IT_0158;
    const ccomplex_t IT_0160 = IT_0082*IT_0159;
    const ccomplex_t IT_0161 = 4*IT_0160;
    const ccomplex_t IT_0162 = m_e*s_14*m_C_1;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0074*IT_0164;
    const ccomplex_t IT_0166 = IT_0077*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = 4*IT_0168;
    const ccomplex_t IT_0170 = 4*IT_0032;
    const ccomplex_t IT_0171 = m_e*s_23*m_C_1;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0075*IT_0173;
    const ccomplex_t IT_0175 = IT_0077*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = 4*IT_0177;
    const ccomplex_t IT_0179 = s_14*s_23;
    const ccomplex_t IT_0180 = 16*IT_0107;
    const ccomplex_t IT_0181 = s_12*IT_0033;
    const ccomplex_t IT_0182 = (-4)*IT_0080;
    const ccomplex_t IT_0183 = IT_0054*IT_0084 + (-0.25)*IT_0134*IT_0138 + (
      -0.25)*IT_0157*IT_0161 + (-0.25)*IT_0162*IT_0169 + (-0.25)*IT_0103*IT_0170
       + (-0.25)*IT_0171*IT_0178 + (-0.25)*IT_0179*IT_0180 + (-0.25)*IT_0181
      *IT_0182;
    const ccomplex_t IT_0184 = IT_0067*IT_0089;
    const ccomplex_t IT_0185 = IT_0069*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0164*IT_0173;
    const ccomplex_t IT_0188 = IT_0077*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = 4*IT_0189;
    const ccomplex_t IT_0191 = s_34*IT_0016;
    const ccomplex_t IT_0192 = 2*IT_0054;
    const ccomplex_t IT_0193 = m_e*s_24*m_C_1;
    const ccomplex_t IT_0194 = m_e*s_13*m_C_1;
    const ccomplex_t IT_0195 = (-2)*IT_0032;
    const ccomplex_t IT_0196 = (-2)*IT_0160;
    const ccomplex_t IT_0197 = -IT_0134;
    const ccomplex_t IT_0198 = 2*IT_0137;
    const ccomplex_t IT_0199 = -IT_0157;
    const ccomplex_t IT_0200 = 4*IT_0080;
    const ccomplex_t IT_0201 = (-2)*IT_0137;
    const ccomplex_t IT_0202 = 2*IT_0160;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0203*IT_0204;
    const ccomplex_t IT_0206 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = s_13 + s_24;
    const ccomplex_t IT_0211 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = IT_0082*IT_0212;
    const ccomplex_t IT_0214 = 2*IT_0072;
    const ccomplex_t IT_0215 = 2*IT_0179;
    const ccomplex_t IT_0216 = IT_0181 + IT_0191 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = -s_14 + -s_23;
    const ccomplex_t IT_0218 = IT_0210 + IT_0217;
    const ccomplex_t IT_0219 = IT_0082*IT_0218;
    const ccomplex_t IT_0220 = (-8)*IT_0219;
    const ccomplex_t IT_0221 = 2*IT_0101;
    const ccomplex_t IT_0222 = IT_0181 + IT_0191 + IT_0214 + IT_0221;
    const ccomplex_t IT_0223 = (-4)*conj(IT_0080);
    const ccomplex_t IT_0224 = 0.125*IT_0223;
    const ccomplex_t IT_0225 = (64*IT_0072 + 32*IT_0101 + 32*IT_0179 + 32
      *IT_0181 + 32*IT_0191)*IT_0209 + 8*(IT_0168 + -conj(IT_0168) + IT_0177 + 
      -conj(IT_0177))*IT_0213 + (-8)*(conj(IT_0107) + conj(IT_0186))*IT_0216 + -
      (conj(IT_0134) + conj(IT_0157) + IT_0197 + IT_0199)*IT_0220 + (-8)*IT_0222
      *(conj(IT_0071) + conj(IT_0100) + (-0.5)*IT_0189 + 0.5*conj(IT_0189) +
       0.125*IT_0200 + IT_0224);
    const ccomplex_t IT_0226 = (-2)*conj(IT_0100);
    const ccomplex_t IT_0227 = -conj(IT_0032);
    const ccomplex_t IT_0228 = -conj(IT_0134);
    const ccomplex_t IT_0229 = -conj(IT_0157);
    const ccomplex_t IT_0230 = conj(IT_0100) + IT_0224;
    const ccomplex_t IT_0231 = 4*conj(IT_0168);
    const ccomplex_t IT_0232 = 4*conj(IT_0177);
    const ccomplex_t IT_0233 = (-2)*conj(IT_0032);
    const ccomplex_t IT_0234 = 0.125*IT_0084;
    const ccomplex_t IT_0235 = 2*conj(IT_0054);
    const ccomplex_t IT_0236 = 0.125*IT_0103;
    const ccomplex_t IT_0237 = conj(IT_0071)*IT_0101 + 0.5*conj(IT_0107)
      *IT_0191 + 0.125*IT_0196*IT_0228 + 0.125*IT_0198*IT_0229 + IT_0072*IT_0230
       + 0.125*IT_0193*IT_0231 + 0.125*IT_0194*IT_0232 + IT_0233*IT_0234 +
       IT_0235*IT_0236;
    const ccomplex_t IT_0238 = 2*conj(IT_0071);
    const ccomplex_t IT_0239 = conj(IT_0189) + IT_0238;
    const ccomplex_t IT_0240 = (-8)*IT_0209;
    const ccomplex_t IT_0241 = 0.125*conj(IT_0134)*IT_0138;
    const ccomplex_t IT_0242 = 0.125*conj(IT_0157);
    const ccomplex_t IT_0243 = IT_0161*IT_0242;
    const ccomplex_t IT_0244 = (-4)*conj(IT_0054);
    const ccomplex_t IT_0245 = 4*conj(IT_0032);
    const ccomplex_t IT_0246 = 8*IT_0107*(IT_0181*IT_0230 + 0.125*IT_0162
      *IT_0231 + 0.125*IT_0171*IT_0232 + 0.5*IT_0191*IT_0239 + 0.125*IT_0216
      *IT_0240 + IT_0241 + IT_0243 + IT_0234*IT_0244 + IT_0236*IT_0245);
    const ccomplex_t IT_0247 = (-8)*conj(IT_0080);
    const ccomplex_t IT_0248 = 0.0625*IT_0247;
    const ccomplex_t IT_0249 = conj(IT_0100) + IT_0248;
    const ccomplex_t IT_0250 = (-4)*conj(IT_0032);
    const ccomplex_t IT_0251 = 4*conj(IT_0054);
    const ccomplex_t IT_0252 = 0.125*IT_0222*IT_0240;
    const ccomplex_t IT_0253 = 0.125*conj(IT_0134)*IT_0161 + IT_0181*conj
      (IT_0186) + conj(IT_0107)*IT_0191 + conj(IT_0168)*IT_0193 + conj(IT_0177)
      *IT_0194 + IT_0101*IT_0239 + IT_0138*IT_0242 + 2*IT_0072*IT_0249 + IT_0234
      *IT_0250 + IT_0236*IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 8*IT_0100*(conj(IT_0107)*IT_0181 + conj(IT_0186
      )*IT_0191 + conj(IT_0177)*IT_0193 + conj(IT_0168)*IT_0194 + IT_0072
      *IT_0239 + IT_0241 + IT_0243 + 2*IT_0101*IT_0249 + IT_0236*IT_0250 +
       IT_0234*IT_0251 + IT_0252);
    const ccomplex_t IT_0255 = 0.25*conj(IT_0157);
    const ccomplex_t IT_0256 = IT_0072*conj(IT_0107) + 0.0625*conj(IT_0134)
      *IT_0161 + IT_0179*conj(IT_0186) + 0.5*IT_0191*IT_0230 + 0.0625*IT_0171
      *IT_0231 + 0.0625*IT_0162*IT_0232 + 0.25*IT_0181*IT_0239 + 0.0625*IT_0216
      *IT_0240 + 0.0625*IT_0103*IT_0244 + 0.0625*IT_0084*IT_0245 + 0.25*IT_0138
      *IT_0255;
    const ccomplex_t IT_0257 = -s_34;
    const ccomplex_t IT_0258 = IT_0033 + IT_0257;
    const ccomplex_t IT_0259 = IT_0016*IT_0258;
    const ccomplex_t IT_0260 = s_12*s_34;
    const ccomplex_t IT_0261 = -IT_0179;
    const ccomplex_t IT_0262 = IT_0101 + IT_0181 + IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = -IT_0260;
    const ccomplex_t IT_0264 = -IT_0101;
    const ccomplex_t IT_0265 = IT_0179 + IT_0181 + IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = (-0.5)*IT_0262;
    const ccomplex_t IT_0267 = (-2)*conj(IT_0177)*((-0.5)*IT_0072*IT_0169 + (
      -0.5)*IT_0101*IT_0178 + (-0.5)*IT_0055*IT_0192 + (-0.5)*IT_0182*IT_0193 +
       IT_0134*IT_0259 + -IT_0157*IT_0265 + IT_0195*IT_0266);
    const ccomplex_t IT_0268 = s_12*IT_0258;
    const ccomplex_t IT_0269 = (-0.125)*conj(IT_0071);
    const ccomplex_t IT_0270 = 8*IT_0032*(conj(IT_0032)*IT_0035 + -conj
      (IT_0054)*IT_0055) + conj(IT_0071)*(IT_0072*IT_0081 + IT_0084*IT_0085) +
       conj(IT_0100)*(IT_0081*IT_0101 + IT_0085*IT_0103) + (-4)*conj(IT_0107)
      *IT_0183 + (-4)*conj(IT_0186)*(IT_0054*IT_0103 + -1./4*IT_0138*IT_0157 + 
      -1./4*IT_0134*IT_0161 + -1./4*IT_0084*IT_0170 + -1./4*IT_0169*IT_0171 + 
      -1./4*IT_0162*IT_0178 + -1./4*IT_0072*IT_0180 + -1./4*IT_0181*IT_0190 + 
      -1./4*IT_0182*IT_0191) + conj(IT_0189)*(IT_0072*IT_0182 + IT_0101*IT_0190 
      + IT_0103*IT_0192 + IT_0169*IT_0193 + IT_0178*IT_0194 + IT_0084*IT_0195 +
       IT_0196*IT_0197 + IT_0198*IT_0199) + 2*conj(IT_0080)*(IT_0032*IT_0103 +
       0.5*IT_0101*IT_0200 + 0.5*IT_0197*IT_0201 + 0.5*IT_0199*IT_0202) + 
      -IT_0209*IT_0225 + 8*IT_0054*(IT_0035*conj(IT_0054) + 0.5*conj(IT_0071)
      *IT_0103 + (-0.25)*IT_0084*(conj(IT_0080) + IT_0226) + IT_0055*IT_0227) + 
      (-2)*IT_0080*(conj(IT_0054)*IT_0084 + IT_0103*IT_0227 + (-0.5)*IT_0201
      *IT_0228 + (-0.5)*IT_0202*IT_0229) + 8*IT_0189*IT_0237 + IT_0246 + 8
      *IT_0071*IT_0253 + IT_0254 + 16*IT_0186*IT_0256 + 8*IT_0177*(conj(IT_0071)
      *IT_0194 + IT_0193*IT_0230 + 0.125*IT_0072*IT_0231 + 0.125*IT_0055*IT_0235
       + (-0.25)*conj(IT_0134)*IT_0259 + 0.125*IT_0233*IT_0262 + IT_0255*IT_0265
      ) + (-2)*conj(IT_0168)*(IT_0032*IT_0055 + (-0.5)*IT_0101*IT_0169 + 2
      *IT_0080*IT_0194 + IT_0157*IT_0259 + -IT_0134*IT_0265 + IT_0192*IT_0266) +
       (-2)*IT_0168*(conj(IT_0032)*IT_0055 + (-4)*conj(IT_0071)*IT_0193 + 2
      *IT_0194*(conj(IT_0080) + IT_0226) + conj(IT_0157)*IT_0259 + -conj(IT_0134
      )*IT_0265 + IT_0235*IT_0266) + IT_0267 + (-8)*IT_0157*((-0.125)*conj
      (IT_0100)*IT_0161 + -conj(IT_0134)*IT_0259 + conj(IT_0157)*IT_0268 +
       IT_0138*IT_0269) + (-8)*IT_0134*((-0.125)*conj(IT_0100)*IT_0138 + -conj
      (IT_0157)*IT_0259 + conj(IT_0134)*IT_0268 + IT_0161*IT_0269);
    return create_ccomplex_return(IT_0270);
}

