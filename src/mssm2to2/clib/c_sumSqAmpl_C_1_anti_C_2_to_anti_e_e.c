#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_e_e(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_e*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = V_u1*U_Wm2;
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = U_d2*V_Wp1;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*(IT_0009 + IT_0011);
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_C_1, 2);
    const ccomplex_t IT_0017 = pow(m_C_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_Z, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + IT_0021 
      + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0001*IT_0010;
    const ccomplex_t IT_0024 = IT_0003*IT_0008;
    const ccomplex_t IT_0025 = IT_0005*(IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = m_e*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0022*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0020 + -IT_0032;
    const ccomplex_t IT_0034 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = IT_0005*(IT_0035 + IT_0037);
    const ccomplex_t IT_0039 = 1.4142135623731*e_em*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = IT_0018*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0036;
    const ccomplex_t IT_0045 = IT_0003*IT_0034;
    const ccomplex_t IT_0046 = IT_0005*(IT_0044 + -IT_0045);
    const ccomplex_t IT_0047 = 1.4142135623731*e_em*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0029*IT_0048;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0051;
    const ccomplex_t IT_0053 = pow(m_e, 2);
    const ccomplex_t IT_0054 = s_34 + IT_0053;
    const ccomplex_t IT_0055 = s_12*IT_0054;
    const ccomplex_t IT_0056 = sin(alpha);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cos(alpha);
    const ccomplex_t IT_0060 = IT_0008*IT_0059;
    const ccomplex_t IT_0061 = IT_0010*IT_0056;
    const ccomplex_t IT_0062 = IT_0005*(IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = 1.4142135623731*e_em*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0058*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0059;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0010*IT_0059;
    const ccomplex_t IT_0073 = IT_0008*IT_0056;
    const ccomplex_t IT_0074 = IT_0005*(IT_0072 + IT_0073);
    const ccomplex_t IT_0075 = 1.4142135623731*e_em*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = cos(theta_W);
    const ccomplex_t IT_0084 = cpow(IT_0083, -1);
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = e_em*IT_0085;
    const ccomplex_t IT_0087 = IT_0005*IT_0083;
    const ccomplex_t IT_0088 = e_em*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0086 + -IT_0088);
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0092 = IT_0085*IT_0091;
    const ccomplex_t IT_0093 = IT_0087*IT_0091;
    const ccomplex_t IT_0094 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0095 = IT_0087*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + -IT_0093 + (-2)
      *IT_0095);
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0090*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = m_e*m_C_2;
    const ccomplex_t IT_0104 = -s_14;
    const ccomplex_t IT_0105 = s_13 + IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = (-4)*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0084;
    const ccomplex_t IT_0109 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0110 = IT_0085*IT_0109;
    const ccomplex_t IT_0111 = IT_0087*IT_0109;
    const ccomplex_t IT_0112 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0113 = IT_0087*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + -IT_0111 + (-2)
      *IT_0113);
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0108*IT_0115;
    const ccomplex_t IT_0117 = IT_0099*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = m_e*m_C_1;
    const ccomplex_t IT_0120 = -s_24;
    const ccomplex_t IT_0121 = s_23 + IT_0120;
    const ccomplex_t IT_0122 = IT_0119*IT_0121;
    const ccomplex_t IT_0123 = 4*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_13 + IT_0016 + IT_0053 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = s_12*IT_0053;
    const ccomplex_t IT_0134 = s_14*s_23;
    const ccomplex_t IT_0135 = s_12*s_34;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_13*s_24;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0133 + IT_0134 + IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = IT_0034*IT_0059;
    const ccomplex_t IT_0141 = IT_0036*IT_0056;
    const ccomplex_t IT_0142 = IT_0005*(IT_0140 + -IT_0141);
    const ccomplex_t IT_0143 = 1.4142135623731*e_em*IT_0142;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = IT_0058*IT_0145;
    const ccomplex_t IT_0147 = IT_0067*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = IT_0036*IT_0059;
    const ccomplex_t IT_0150 = IT_0034*IT_0056;
    const ccomplex_t IT_0151 = IT_0005*(IT_0149 + IT_0150);
    const ccomplex_t IT_0152 = 1.4142135623731*e_em*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = IT_0071*IT_0154;
    const ccomplex_t IT_0156 = IT_0079*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = -IT_0148 + -IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = IT_0129*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0090*IT_0115;
    const ccomplex_t IT_0167 = IT_0099*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0170 = IT_0127*IT_0161;
    const ccomplex_t IT_0171 = IT_0129*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = (-4)*conj(IT_0173);
    const ccomplex_t IT_0175 = 0.125*IT_0174;
    const ccomplex_t IT_0176 = conj(IT_0102) + IT_0175;
    const ccomplex_t IT_0177 = IT_0125*IT_0160;
    const ccomplex_t IT_0178 = IT_0129*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = (-4)*conj(IT_0180);
    const ccomplex_t IT_0182 = 0.125*IT_0181;
    const ccomplex_t IT_0183 = conj(IT_0118) + IT_0182;
    const ccomplex_t IT_0184 = m_e*s_23*m_C_1;
    const ccomplex_t IT_0185 = 4*conj(IT_0132);
    const ccomplex_t IT_0186 = m_e*s_14*m_C_2;
    const ccomplex_t IT_0187 = 4*conj(IT_0165);
    const ccomplex_t IT_0188 = m_e*s_13*m_C_2;
    const ccomplex_t IT_0189 = 8*IT_0132;
    const ccomplex_t IT_0190 = 8*conj(IT_0132);
    const ccomplex_t IT_0191 = 8*IT_0165;
    const ccomplex_t IT_0192 = 8*conj(IT_0165);
    const ccomplex_t IT_0193 = (-4)*IT_0132;
    const ccomplex_t IT_0194 = (-4)*conj(IT_0132);
    const ccomplex_t IT_0195 = (-4)*IT_0180;
    const ccomplex_t IT_0196 = m_e*s_24*m_C_1;
    const ccomplex_t IT_0197 = (-4)*IT_0173;
    const ccomplex_t IT_0198 = s_13 + s_14;
    const ccomplex_t IT_0199 = IT_0103*IT_0198;
    const ccomplex_t IT_0200 = IT_0097*IT_0108;
    const ccomplex_t IT_0201 = IT_0099*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = 4*conj(IT_0033);
    const ccomplex_t IT_0205 = 4*IT_0052;
    const ccomplex_t IT_0206 = 4*conj(IT_0052);
    const ccomplex_t IT_0207 = (-4)*conj(IT_0033);
    const ccomplex_t IT_0208 = (-4)*IT_0052;
    const ccomplex_t IT_0209 = (-4)*conj(IT_0052);
    const ccomplex_t IT_0210 = (-2)*conj(IT_0033);
    const ccomplex_t IT_0211 = conj(IT_0052)*IT_0180 + IT_0052*conj(IT_0180) +
       2*IT_0033*(conj(IT_0102) + -1./2*conj(IT_0173) + -conj(IT_0203)) + 0.5
      *IT_0102*IT_0204 + 0.5*conj(IT_0168)*IT_0205 + 0.5*IT_0168*IT_0206 + 0.5
      *IT_0203*IT_0207 + 0.5*conj(IT_0118)*IT_0208 + 0.5*IT_0118*IT_0209 + 0.5
      *IT_0173*IT_0210;
    const ccomplex_t IT_0212 = m_C_1*m_C_2;
    const ccomplex_t IT_0213 = -IT_0053;
    const ccomplex_t IT_0214 = s_34 + IT_0213;
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = (-2)*IT_0215;
    const ccomplex_t IT_0217 = -IT_0082;
    const ccomplex_t IT_0218 = -conj(IT_0082);
    const ccomplex_t IT_0219 = -IT_0158;
    const ccomplex_t IT_0220 = -conj(IT_0158);
    const ccomplex_t IT_0221 = (-2)*IT_0106;
    const ccomplex_t IT_0222 = 2*IT_0122;
    const ccomplex_t IT_0223 = s_23 + s_24;
    const ccomplex_t IT_0224 = IT_0119*IT_0223;
    const ccomplex_t IT_0225 = -conj(IT_0052);
    const ccomplex_t IT_0226 = IT_0052*conj(IT_0173) + 2*IT_0033*(conj(IT_0118
      ) + -conj(IT_0168) + -1./2*conj(IT_0180)) + 0.5*IT_0118*IT_0204 + 0.5*conj
      (IT_0203)*IT_0205 + 0.5*IT_0203*IT_0206 + 0.5*IT_0168*IT_0207 + 0.5*conj
      (IT_0102)*IT_0208 + 0.5*IT_0102*IT_0209 + 0.5*IT_0180*IT_0210 + -IT_0173
      *IT_0225;
    const ccomplex_t IT_0227 = -IT_0134;
    const ccomplex_t IT_0228 = IT_0133 + IT_0135 + IT_0137 + IT_0227;
    const ccomplex_t IT_0229 = 2*IT_0033;
    const ccomplex_t IT_0230 = 0.125*IT_0195;
    const ccomplex_t IT_0231 = IT_0118 + IT_0230;
    const ccomplex_t IT_0232 = 8*IT_0102;
    const ccomplex_t IT_0233 = IT_0197 + IT_0232;
    const ccomplex_t IT_0234 = 0.125*IT_0158;
    const ccomplex_t IT_0235 = 4*IT_0165;
    const ccomplex_t IT_0236 = 0.125*IT_0186;
    const ccomplex_t IT_0237 = 4*IT_0132;
    const ccomplex_t IT_0238 = 0.125*IT_0184;
    const ccomplex_t IT_0239 = 16*conj(IT_0168)*(0.0625*IT_0082*IT_0123 +
       IT_0134*IT_0168 + 0.5*IT_0133*IT_0231 + 0.0625*IT_0169*IT_0233 + 0.5
      *IT_0107*IT_0234 + 0.5*IT_0235*IT_0236 + 0.5*IT_0237*IT_0238);
    const ccomplex_t IT_0240 = 16*IT_0203;
    const ccomplex_t IT_0241 = IT_0054*IT_0212;
    const ccomplex_t IT_0242 = 0.25*conj(IT_0132);
    const ccomplex_t IT_0243 = (-8)*conj(IT_0180);
    const ccomplex_t IT_0244 = 0.0625*IT_0243;
    const ccomplex_t IT_0245 = conj(IT_0118) + IT_0244;
    const ccomplex_t IT_0246 = 4*IT_0180;
    const ccomplex_t IT_0247 = (-8)*IT_0173;
    const ccomplex_t IT_0248 = (-8)*conj(IT_0173);
    const ccomplex_t IT_0249 = 0.25*conj(IT_0165);
    const ccomplex_t IT_0250 = (-8)*IT_0180;
    const ccomplex_t IT_0251 = 0.25*conj(IT_0118);
    const ccomplex_t IT_0252 = IT_0173*conj(IT_0173) + IT_0237*IT_0242 + 4
      *IT_0118*IT_0245 + 0.25*conj(IT_0180)*IT_0246 + 0.25*conj(IT_0102)*IT_0247
       + 4*IT_0102*(conj(IT_0102) + 0.0625*IT_0248) + IT_0235*IT_0249 + IT_0250
      *IT_0251;
    const ccomplex_t IT_0253 = m_C_1*m_C_2*IT_0053;
    const ccomplex_t IT_0254 = IT_0173*conj(IT_0180) + 0.25*IT_0165*IT_0185 +
       4*IT_0168*conj(IT_0203) + 0.25*conj(IT_0168)*IT_0240 + 4*IT_0102*IT_0245 
      + 0.25*conj(IT_0173)*IT_0246 + 0.25*IT_0118*IT_0248 + IT_0237*IT_0249 + 4
      *conj(IT_0102)*(IT_0118 + 0.0625*IT_0250) + IT_0247*IT_0251;
    const ccomplex_t IT_0255 = -s_34;
    const ccomplex_t IT_0256 = IT_0053 + IT_0255;
    const ccomplex_t IT_0257 = s_12*IT_0256;
    const ccomplex_t IT_0258 = (-8)*IT_0215;
    const ccomplex_t IT_0259 = (-0.125)*conj(IT_0102);
    const ccomplex_t IT_0260 = (-0.125)*conj(IT_0118);
    const ccomplex_t IT_0261 = 8*(IT_0033*conj(IT_0033) + IT_0052*conj(IT_0052
      ))*IT_0055 + conj(IT_0082)*(IT_0102*IT_0107 + IT_0118*IT_0123 + 2*IT_0132
      *IT_0139) + conj(IT_0158)*(IT_0107*IT_0118 + IT_0102*IT_0123 + 2*IT_0139
      *IT_0165) + IT_0168*(conj(IT_0082)*IT_0123 + IT_0107*conj(IT_0158) + 8
      *IT_0169*IT_0176 + 8*IT_0133*IT_0183 + IT_0184*IT_0185 + IT_0186*IT_0187) 
      + IT_0188*(IT_0165*IT_0181 + conj(IT_0102)*IT_0189 + IT_0102*IT_0190 +
       conj(IT_0118)*IT_0191 + IT_0118*IT_0192 + conj(IT_0173)*IT_0193 + IT_0173
      *IT_0194 + conj(IT_0165)*IT_0195) + IT_0196*(IT_0165*IT_0174 + conj
      (IT_0118)*IT_0189 + IT_0118*IT_0190 + conj(IT_0102)*IT_0191 + IT_0102
      *IT_0192 + conj(IT_0180)*IT_0193 + IT_0180*IT_0194 + conj(IT_0165)*IT_0197
      ) + 2*IT_0199*IT_0211 + IT_0216*(conj(IT_0165)*IT_0217 + IT_0165*IT_0218 +
       conj(IT_0132)*IT_0219 + IT_0132*IT_0220) + (conj(IT_0173)*IT_0217 +
       IT_0173*IT_0218 + conj(IT_0180)*IT_0219 + IT_0180*IT_0220)*IT_0221 + 
      (conj(IT_0180)*IT_0217 + IT_0180*IT_0218 + conj(IT_0173)*IT_0219 + IT_0173
      *IT_0220)*IT_0222 + 2*IT_0224*IT_0226 + 2*IT_0228*(conj(IT_0033)*IT_0132 +
       -IT_0052*conj(IT_0165) + IT_0165*IT_0225 + 0.5*conj(IT_0132)*IT_0229) +
       IT_0203*(conj(IT_0082)*IT_0107 + IT_0123*conj(IT_0158) + 8*IT_0133
      *IT_0176 + 8*IT_0169*IT_0183 + 8*IT_0185*IT_0236 + 8*IT_0187*IT_0238) +
       IT_0239 + conj(IT_0203)*(IT_0082*IT_0107 + 8*IT_0169*IT_0231 + IT_0133
      *IT_0233 + 8*IT_0123*IT_0234 + 8*IT_0236*IT_0237 + 8*IT_0235*IT_0238 +
       IT_0134*IT_0240) + (-8)*IT_0241*(conj(IT_0052)*(IT_0033 + 0.25*IT_0132) +
       conj(IT_0033)*(IT_0052 + (-0.25)*IT_0165) + (-0.125)*conj(IT_0165)
      *IT_0229 + IT_0052*IT_0242) + 4*IT_0137*IT_0252 + 4*IT_0253*IT_0254 + (-8)
      *IT_0158*((-0.25)*IT_0139*conj(IT_0165) + conj(IT_0158)*IT_0257 + (-0.125)
      *conj(IT_0082)*IT_0258 + IT_0123*IT_0259 + IT_0107*IT_0260) + (-8)*IT_0082
      *((-0.25)*conj(IT_0132)*IT_0139 + conj(IT_0082)*IT_0257 + (-0.125)*conj
      (IT_0158)*IT_0258 + IT_0107*IT_0259 + IT_0123*IT_0260);
    return create_ccomplex_return(IT_0261);
}

