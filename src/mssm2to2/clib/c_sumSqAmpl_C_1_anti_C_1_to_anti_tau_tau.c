#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_tau_tau(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*(IT_0009 + IT_0012);
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = IT_0003*IT_0010;
    const ccomplex_t IT_0025 = IT_0008*IT_0011;
    const ccomplex_t IT_0026 = IT_0005*(IT_0024 + -IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0023*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0021 + -IT_0033;
    const ccomplex_t IT_0035 = pow(m_tau, 2);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = s_34 + IT_0036;
    const ccomplex_t IT_0038 = s_12*IT_0037;
    const ccomplex_t IT_0039 = U_d1*V_Wp1;
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = V_u1*U_Wm1;
    const ccomplex_t IT_0042 = IT_0011*IT_0041;
    const ccomplex_t IT_0043 = IT_0005*(IT_0040 + IT_0042);
    const ccomplex_t IT_0044 = 1.4142135623731*e_em*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = IT_0019*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0003*IT_0041;
    const ccomplex_t IT_0051 = IT_0011*IT_0039;
    const ccomplex_t IT_0052 = IT_0005*(IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0023*IT_0055;
    const ccomplex_t IT_0057 = IT_0031*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0049 + -IT_0058;
    const ccomplex_t IT_0060 = IT_0018*IT_0037;
    const ccomplex_t IT_0061 = cos(theta_W);
    const ccomplex_t IT_0062 = cpow(IT_0061, -1);
    const ccomplex_t IT_0063 = IT_0004*IT_0062;
    const ccomplex_t IT_0064 = e_em*IT_0063;
    const ccomplex_t IT_0065 = IT_0005*IT_0061;
    const ccomplex_t IT_0066 = e_em*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0064 + -IT_0066);
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0070 = IT_0063*IT_0069;
    const ccomplex_t IT_0071 = IT_0065*IT_0069;
    const ccomplex_t IT_0072 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0073 = IT_0065*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + -IT_0071 + (-2)
      *IT_0073);
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0068*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0062;
    const ccomplex_t IT_0082 = IT_0075*IT_0081;
    const ccomplex_t IT_0083 = IT_0077*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = m_C_1*m_tau;
    const ccomplex_t IT_0087 = -s_14;
    const ccomplex_t IT_0088 = s_13 + IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = (-4)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0035 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0100 = IT_0063*IT_0099;
    const ccomplex_t IT_0101 = IT_0065*IT_0099;
    const ccomplex_t IT_0102 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0103 = IT_0065*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + -IT_0101 + (-2)
      *IT_0103);
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0068*IT_0105;
    const ccomplex_t IT_0107 = IT_0077*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0081*IT_0105;
    const ccomplex_t IT_0110 = IT_0077*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -s_24;
    const ccomplex_t IT_0113 = s_23 + IT_0112;
    const ccomplex_t IT_0114 = IT_0086*IT_0113;
    const ccomplex_t IT_0115 = 4*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = IT_0095*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = s_12*IT_0035;
    const ccomplex_t IT_0125 = s_14*s_23;
    const ccomplex_t IT_0126 = s_12*s_34;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = s_13*s_24;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = IT_0125 + IT_0127 + IT_0129;
    const ccomplex_t IT_0131 = IT_0124 + IT_0130;
    const ccomplex_t IT_0132 = s_23 + s_24;
    const ccomplex_t IT_0133 = IT_0086*IT_0132;
    const ccomplex_t IT_0134 = sin(beta);
    const ccomplex_t IT_0135 = e_em*m_tau*IT_0000*IT_0002*IT_0005*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = IT_0001*IT_0041;
    const ccomplex_t IT_0138 = IT_0039*IT_0134;
    const ccomplex_t IT_0139 = IT_0005*(IT_0137 + IT_0138);
    const ccomplex_t IT_0140 = 1.4142135623731*e_em*IT_0139;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0136*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0144 = IT_0142*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = pow(m_Z, 2);
    const ccomplex_t IT_0147 = cpow((-2)*s_12 + (-2)*IT_0018 + IT_0146 + 
      -reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0001*IT_0039;
    const ccomplex_t IT_0149 = IT_0041*IT_0134;
    const ccomplex_t IT_0150 = IT_0005*(IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = 1.4142135623731*e_em*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = e_em*m_tau*IT_0000*IT_0005;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = IT_0147*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = -IT_0145 + -IT_0157;
    const ccomplex_t IT_0159 = 4*IT_0158;
    const ccomplex_t IT_0160 = IT_0018*IT_0035;
    const ccomplex_t IT_0161 = IT_0093*IT_0119;
    const ccomplex_t IT_0162 = IT_0095*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = 8*IT_0164;
    const ccomplex_t IT_0166 = IT_0092*IT_0117;
    const ccomplex_t IT_0167 = IT_0095*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = 8*IT_0169;
    const ccomplex_t IT_0171 = s_24*m_C_1*m_tau;
    const ccomplex_t IT_0172 = (-8)*IT_0123;
    const ccomplex_t IT_0173 = s_13*m_C_1*m_tau;
    const ccomplex_t IT_0174 = (-8)*IT_0098;
    const ccomplex_t IT_0175 = s_13 + s_14;
    const ccomplex_t IT_0176 = IT_0086*IT_0175;
    const ccomplex_t IT_0177 = IT_0001*IT_0010;
    const ccomplex_t IT_0178 = IT_0008*IT_0134;
    const ccomplex_t IT_0179 = IT_0005*(IT_0177 + IT_0178);
    const ccomplex_t IT_0180 = 1.4142135623731*e_em*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0136*IT_0181;
    const ccomplex_t IT_0183 = IT_0143*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = IT_0001*IT_0008;
    const ccomplex_t IT_0186 = IT_0010*IT_0134;
    const ccomplex_t IT_0187 = IT_0005*(IT_0185 + -IT_0186);
    const ccomplex_t IT_0188 = 1.4142135623731*e_em*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = IT_0154*IT_0189;
    const ccomplex_t IT_0191 = IT_0147*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = -IT_0184 + -IT_0192;
    const ccomplex_t IT_0194 = (-4)*IT_0193;
    const ccomplex_t IT_0195 = 4*IT_0193;
    const ccomplex_t IT_0196 = s_34*IT_0018;
    const ccomplex_t IT_0197 = 8*IT_0080;
    const ccomplex_t IT_0198 = s_23*m_C_1*m_tau;
    const ccomplex_t IT_0199 = (-4)*IT_0123;
    const ccomplex_t IT_0200 = s_14*m_C_1*m_tau;
    const ccomplex_t IT_0201 = (-4)*IT_0098;
    const ccomplex_t IT_0202 = (-4)*IT_0158;
    const ccomplex_t IT_0203 = (-2)*IT_0034;
    const ccomplex_t IT_0204 = (-4)*IT_0169;
    const ccomplex_t IT_0205 = s_34 + IT_0035;
    const ccomplex_t IT_0206 = IT_0018*IT_0205;
    const ccomplex_t IT_0207 = -IT_0124;
    const ccomplex_t IT_0208 = IT_0130 + IT_0207;
    const ccomplex_t IT_0209 = (-2)*IT_0098;
    const ccomplex_t IT_0210 = s_12*IT_0205;
    const ccomplex_t IT_0211 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = IT_0212*IT_0213;
    const ccomplex_t IT_0215 = cpow(s_12 + IT_0018 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = 0.5*IT_0217;
    const ccomplex_t IT_0219 = 8*IT_0196;
    const ccomplex_t IT_0220 = 16*IT_0160;
    const ccomplex_t IT_0221 = 8*IT_0124;
    const ccomplex_t IT_0222 = 16*IT_0125;
    const ccomplex_t IT_0223 = IT_0219 + IT_0220 + IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = 16*IT_0128;
    const ccomplex_t IT_0225 = IT_0219 + IT_0220 + IT_0221 + IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0196;
    const ccomplex_t IT_0227 = 0.5*IT_0124 + IT_0128 + IT_0160 + IT_0226;
    const ccomplex_t IT_0228 = 8*conj(IT_0164);
    const ccomplex_t IT_0229 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0230 = s_13 + s_24;
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = IT_0086*IT_0231;
    const ccomplex_t IT_0233 = (-8)*conj(IT_0123);
    const ccomplex_t IT_0234 = (-8)*conj(IT_0098);
    const ccomplex_t IT_0235 = 2*IT_0128;
    const ccomplex_t IT_0236 = 2*IT_0160;
    const ccomplex_t IT_0237 = IT_0124 + IT_0196 + IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = -s_14 + -s_23;
    const ccomplex_t IT_0239 = IT_0230 + IT_0238;
    const ccomplex_t IT_0240 = IT_0086*IT_0239;
    const ccomplex_t IT_0241 = (-8)*IT_0240;
    const ccomplex_t IT_0242 = -IT_0059;
    const ccomplex_t IT_0243 = -IT_0034;
    const ccomplex_t IT_0244 = (32*IT_0124 + 32*IT_0125 + 32*IT_0128 + 64
      *IT_0160 + 32*IT_0196)*IT_0218 + -(conj(IT_0085) + conj(IT_0108))*IT_0223 
      + -(conj(IT_0080) + conj(IT_0111))*IT_0225 + 8*IT_0227*(IT_0164 + (-0.125)
      *IT_0228) + (-8)*IT_0232*(IT_0123 + 0.125*IT_0233) + (-8)*IT_0232*(IT_0098
       + 0.125*IT_0234) + (-4)*(conj(IT_0169) + 0.25*IT_0204)*IT_0237 + -IT_0241
      *(conj(IT_0034) + conj(IT_0059) + IT_0242 + IT_0243);
    const ccomplex_t IT_0245 = (-2)*IT_0114;
    const ccomplex_t IT_0246 = 2*IT_0089;
    const ccomplex_t IT_0247 = IT_0034*IT_0115;
    const ccomplex_t IT_0248 = (-2)*conj(IT_0059);
    const ccomplex_t IT_0249 = (-2)*conj(IT_0034);
    const ccomplex_t IT_0250 = -conj(IT_0098);
    const ccomplex_t IT_0251 = (-4)*conj(IT_0123);
    const ccomplex_t IT_0252 = -conj(IT_0059);
    const ccomplex_t IT_0253 = 0.25*IT_0246;
    const ccomplex_t IT_0254 = -conj(IT_0034);
    const ccomplex_t IT_0255 = 0.25*IT_0245;
    const ccomplex_t IT_0256 = 0.5*conj(IT_0158);
    const ccomplex_t IT_0257 = conj(IT_0085)*IT_0124 + IT_0128*conj(IT_0164) +
       IT_0160*conj(IT_0169) + (-0.5)*IT_0133*conj(IT_0193) + conj(IT_0108)
      *IT_0196 + IT_0171*IT_0250 + 0.25*IT_0173*IT_0251 + IT_0252*IT_0253 +
       IT_0254*IT_0255 + IT_0176*IT_0256;
    const ccomplex_t IT_0258 = 4*conj(IT_0158);
    const ccomplex_t IT_0259 = 8*conj(IT_0169);
    const ccomplex_t IT_0260 = (-4)*conj(IT_0193);
    const ccomplex_t IT_0261 = -IT_0218;
    const ccomplex_t IT_0262 = conj(IT_0034)*IT_0115;
    const ccomplex_t IT_0263 = conj(IT_0059)*IT_0090 + 16*conj(IT_0080)
      *IT_0128 + IT_0128*IT_0228 + IT_0173*IT_0233 + IT_0171*IT_0234 + IT_0176
      *IT_0258 + IT_0160*(16*conj(IT_0111) + IT_0259) + IT_0133*IT_0260 +
       IT_0225*IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 2*conj(IT_0111);
    const ccomplex_t IT_0265 = conj(IT_0085)*IT_0196;
    const ccomplex_t IT_0266 = 0.125*conj(IT_0034)*IT_0090 + 0.125*conj
      (IT_0059)*IT_0115 + conj(IT_0108)*IT_0124 + 2*IT_0160*(conj(IT_0080) +
       0.0625*IT_0228) + 0.125*IT_0171*IT_0233 + 0.125*IT_0173*IT_0234 + 0.125
      *IT_0133*IT_0258 + 0.125*IT_0128*IT_0259 + 0.125*IT_0176*IT_0260 + 0.125
      *IT_0225*IT_0261 + IT_0128*IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = 4*IT_0169*(conj(IT_0108)*IT_0124 + IT_0160*conj
      (IT_0164) + IT_0128*conj(IT_0169) + (-0.5)*IT_0176*conj(IT_0193) + 0.25
      *IT_0171*IT_0251 + IT_0253*IT_0254 + IT_0252*IT_0255 + IT_0133*IT_0256 +
       IT_0265);
    const ccomplex_t IT_0268 = 8*conj(IT_0080);
    const ccomplex_t IT_0269 = 0.25*IT_0268;
    const ccomplex_t IT_0270 = conj(IT_0164) + IT_0269;
    const ccomplex_t IT_0271 = 4*conj(IT_0193);
    const ccomplex_t IT_0272 = (-4)*conj(IT_0098);
    const ccomplex_t IT_0273 = (-4)*conj(IT_0158);
    const ccomplex_t IT_0274 = 8*conj(IT_0111) + 4*conj(IT_0169);
    const ccomplex_t IT_0275 = 16*conj(IT_0108);
    const ccomplex_t IT_0276 = IT_0223*IT_0261;
    const ccomplex_t IT_0277 = conj(IT_0034)*IT_0090 + conj(IT_0059)*IT_0115 +
       16*conj(IT_0085)*IT_0160 + IT_0198*IT_0251 + 4*IT_0196*IT_0270 + IT_0176
      *IT_0271 + IT_0200*IT_0272 + IT_0133*IT_0273 + IT_0124*IT_0274 + IT_0125
      *IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = conj(IT_0059)*IT_0090 + 16*conj(IT_0085)
      *IT_0125 + IT_0200*IT_0251 + IT_0262 + 4*IT_0124*IT_0270 + IT_0133*IT_0271
       + IT_0198*IT_0272 + IT_0176*IT_0273 + IT_0196*IT_0274 + IT_0160*IT_0275 +
       IT_0276;
    const ccomplex_t IT_0279 = 8*IT_0034*(conj(IT_0034)*IT_0038 + -conj
      (IT_0059)*IT_0060) + 8*IT_0059*(IT_0038*conj(IT_0059) + 1./8*(conj(IT_0080
      ) + conj(IT_0085))*IT_0090 + -IT_0060*(conj(IT_0034) + 0.25*conj(IT_0098))
       + 1./8*(conj(IT_0108) + conj(IT_0111))*IT_0115 + -1./4*conj(IT_0123)
      *IT_0131) + conj(IT_0111)*(IT_0034*IT_0090 + IT_0133*IT_0159 + IT_0160
      *IT_0165 + IT_0128*IT_0170 + IT_0171*IT_0172 + IT_0173*IT_0174 + IT_0176
      *IT_0194) + conj(IT_0108)*(IT_0034*IT_0090 + IT_0176*IT_0195 + IT_0196
      *IT_0197 + IT_0198*IT_0199 + IT_0200*IT_0201 + IT_0133*IT_0202) + IT_0060
      *conj(IT_0123)*IT_0203 + conj(IT_0098)*(IT_0131*IT_0203 + IT_0173*IT_0204)
       + (conj(IT_0158)*IT_0206 + conj(IT_0193)*IT_0208)*IT_0209 + (-2)*IT_0193*
      (IT_0133*conj(IT_0164) + conj(IT_0169)*IT_0176 + -(conj(IT_0123) + (-4)
      *conj(IT_0158))*IT_0206 + (-4)*conj(IT_0193)*IT_0210 + (-0.5)*IT_0208
      *IT_0211) + 2*IT_0158*(IT_0133*conj(IT_0169) + conj(IT_0164)*IT_0176 +
       conj(IT_0123)*IT_0208 + 4*conj(IT_0158)*IT_0210 + (-4)*IT_0206*(conj
      (IT_0193) + (-0.125)*IT_0211)) + -IT_0218*IT_0244 + conj(IT_0164)*(IT_0173
      *IT_0199 + IT_0243*IT_0245 + IT_0242*IT_0246) + conj(IT_0169)*(IT_0171
      *IT_0199 + IT_0242*IT_0245 + IT_0243*IT_0246) + conj(IT_0080)*(IT_0128
      *IT_0165 + IT_0160*IT_0170 + IT_0172*IT_0173 + IT_0171*IT_0174 + IT_0159
      *IT_0176 + IT_0133*IT_0194 + IT_0247) + conj(IT_0085)*(IT_0133*IT_0195 +
       IT_0124*IT_0197 + IT_0199*IT_0200 + IT_0198*IT_0201 + IT_0176*IT_0202 +
       IT_0247) + 4*IT_0098*(conj(IT_0123)*IT_0160 + -conj(IT_0164)*IT_0171 + 
      -conj(IT_0169)*IT_0173 + 0.25*IT_0060*IT_0248 + 0.25*IT_0131*IT_0249 + 
      -IT_0128*IT_0250) + 4*IT_0123*(conj(IT_0123)*IT_0128 + conj(IT_0098)
      *IT_0160 + 0.5*conj(IT_0193)*IT_0206 + 0.25*IT_0131*IT_0248 + 0.25*IT_0060
      *IT_0249 + IT_0208*IT_0256) + 4*IT_0164*IT_0257 + IT_0080*IT_0263 + 8
      *IT_0111*IT_0266 + IT_0267 + IT_0108*IT_0277 + IT_0085*IT_0278;
    return create_ccomplex_return(IT_0279);
}

