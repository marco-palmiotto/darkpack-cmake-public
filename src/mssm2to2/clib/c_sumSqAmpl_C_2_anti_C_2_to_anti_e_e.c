#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_e_e(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = pow(m_e, 2);
    const ccomplex_t IT_0002 = pow(m_C_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = m_e*e_em*IT_0004*IT_0006*IT_0008*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0020 = IT_0005*IT_0019;
    const ccomplex_t IT_0021 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = IT_0008*(IT_0020 + IT_0022);
    const ccomplex_t IT_0024 = 1.4142135623731*e_em*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0018*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + (-2)*IT_0002 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0021;
    const ccomplex_t IT_0031 = IT_0016*IT_0019;
    const ccomplex_t IT_0032 = IT_0008*(IT_0030 + -IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = m_e*e_em*IT_0004*IT_0008;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = pow(m_Z, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + (-2)*IT_0002 + IT_0038 + 
      -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0029 + -IT_0041;
    const ccomplex_t IT_0043 = s_12*IT_0001;
    const ccomplex_t IT_0044 = s_12*s_34;
    const ccomplex_t IT_0045 = s_14*s_23;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0000 + IT_0043 + IT_0044 + IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0008;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = IT_0003*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 2*IT_0053;
    const ccomplex_t IT_0055 = s_34 + IT_0001;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0015;
    const ccomplex_t IT_0058 = cos(theta_W);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = IT_0007*IT_0059;
    const ccomplex_t IT_0061 = e_em*IT_0060;
    const ccomplex_t IT_0062 = IT_0008*IT_0058;
    const ccomplex_t IT_0063 = e_em*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0061 + -IT_0063);
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0067 = IT_0060*IT_0066;
    const ccomplex_t IT_0068 = IT_0062*IT_0066;
    const ccomplex_t IT_0069 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0070 = IT_0062*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + -IT_0068 + (-2)
      *IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0065*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = s_34*IT_0002;
    const ccomplex_t IT_0079 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0080 = IT_0060*IT_0079;
    const ccomplex_t IT_0081 = IT_0062*IT_0079;
    const ccomplex_t IT_0082 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + -IT_0081 + (-2)
      *IT_0083);
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0065*IT_0085;
    const ccomplex_t IT_0087 = IT_0074*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = 8*IT_0088;
    const ccomplex_t IT_0090 = IT_0012*IT_0048;
    const ccomplex_t IT_0091 = IT_0003*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = 8*IT_0092;
    const ccomplex_t IT_0094 = IT_0001*IT_0002;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*e_em*IT_0007*IT_0059;
    const ccomplex_t IT_0096 = IT_0085*IT_0095;
    const ccomplex_t IT_0097 = IT_0074*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = 16*IT_0098;
    const ccomplex_t IT_0100 = m_e*s_13*m_C_2;
    const ccomplex_t IT_0101 = (-8)*IT_0015;
    const ccomplex_t IT_0102 = m_e*s_24*m_C_2;
    const ccomplex_t IT_0103 = m_e*s_14*m_C_2;
    const ccomplex_t IT_0104 = m_e*s_23*m_C_2;
    const ccomplex_t IT_0105 = IT_0072*IT_0095;
    const ccomplex_t IT_0106 = IT_0074*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = (-4)*IT_0053;
    const ccomplex_t IT_0110 = IT_0010*IT_0050;
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (-4)*IT_0015;
    const ccomplex_t IT_0114 = m_e*m_C_2;
    const ccomplex_t IT_0115 = s_13 + s_14;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = V_u2*U_Wm2;
    const ccomplex_t IT_0118 = IT_0005*IT_0117;
    const ccomplex_t IT_0119 = U_d2*V_Wp2;
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = IT_0008*(IT_0118 + IT_0120);
    const ccomplex_t IT_0122 = 1.4142135623731*e_em*IT_0121;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0018*IT_0123;
    const ccomplex_t IT_0125 = IT_0027*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0005*IT_0119;
    const ccomplex_t IT_0128 = IT_0016*IT_0117;
    const ccomplex_t IT_0129 = IT_0008*(IT_0127 + -IT_0128);
    const ccomplex_t IT_0130 = 1.4142135623731*e_em*IT_0129;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0036*IT_0131;
    const ccomplex_t IT_0133 = IT_0039*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = -IT_0126 + -IT_0134;
    const ccomplex_t IT_0136 = s_23 + s_24;
    const ccomplex_t IT_0137 = IT_0114*IT_0136;
    const ccomplex_t IT_0138 = (-4)*conj(IT_0015);
    const ccomplex_t IT_0139 = (-2)*IT_0015;
    const ccomplex_t IT_0140 = cos(alpha);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_e*e_em*IT_0004*IT_0006
      *IT_0008*IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = IT_0021*IT_0140;
    const ccomplex_t IT_0144 = sin(alpha);
    const ccomplex_t IT_0145 = IT_0019*IT_0144;
    const ccomplex_t IT_0146 = IT_0008*(IT_0143 + IT_0145);
    const ccomplex_t IT_0147 = 1.4142135623731*e_em*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0142*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_12 + (-2)*IT_0002 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*m_e*e_em*IT_0004*IT_0006
      *IT_0008*IT_0144;
    const ccomplex_t IT_0155 = 0.5*IT_0154;
    const ccomplex_t IT_0156 = IT_0019*IT_0140;
    const ccomplex_t IT_0157 = IT_0021*IT_0144;
    const ccomplex_t IT_0158 = IT_0008*(IT_0156 + -IT_0157);
    const ccomplex_t IT_0159 = 1.4142135623731*e_em*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0155*IT_0161;
    const ccomplex_t IT_0163 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0153 + -IT_0165;
    const ccomplex_t IT_0167 = -s_34;
    const ccomplex_t IT_0168 = IT_0001 + IT_0167;
    const ccomplex_t IT_0169 = IT_0002*IT_0168;
    const ccomplex_t IT_0170 = -s_24;
    const ccomplex_t IT_0171 = s_23 + IT_0170;
    const ccomplex_t IT_0172 = IT_0114*IT_0171;
    const ccomplex_t IT_0173 = (-2)*IT_0172;
    const ccomplex_t IT_0174 = -IT_0092;
    const ccomplex_t IT_0175 = IT_0119*IT_0140;
    const ccomplex_t IT_0176 = IT_0117*IT_0144;
    const ccomplex_t IT_0177 = IT_0008*(IT_0175 + IT_0176);
    const ccomplex_t IT_0178 = 1.4142135623731*e_em*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = 0.5*IT_0179;
    const ccomplex_t IT_0181 = IT_0142*IT_0180;
    const ccomplex_t IT_0182 = IT_0151*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = IT_0117*IT_0140;
    const ccomplex_t IT_0185 = IT_0119*IT_0144;
    const ccomplex_t IT_0186 = IT_0008*(IT_0184 + -IT_0185);
    const ccomplex_t IT_0187 = 1.4142135623731*e_em*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0155*IT_0189;
    const ccomplex_t IT_0191 = IT_0163*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = -IT_0183 + -IT_0192;
    const ccomplex_t IT_0194 = -IT_0044;
    const ccomplex_t IT_0195 = -IT_0000;
    const ccomplex_t IT_0196 = IT_0043 + IT_0045 + IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = -s_14;
    const ccomplex_t IT_0198 = s_13 + IT_0197;
    const ccomplex_t IT_0199 = IT_0114*IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = (-4)*IT_0199;
    const ccomplex_t IT_0202 = conj(IT_0088) + conj(IT_0098);
    const ccomplex_t IT_0203 = 4*IT_0172;
    const ccomplex_t IT_0204 = conj(IT_0077) + conj(IT_0108);
    const ccomplex_t IT_0205 = 2*conj(IT_0015);
    const ccomplex_t IT_0206 = 8*IT_0169;
    const ccomplex_t IT_0207 = s_12*IT_0168;
    const ccomplex_t IT_0208 = (-8)*IT_0207;
    const ccomplex_t IT_0209 = (-2)*conj(IT_0053);
    const ccomplex_t IT_0210 = -conj(IT_0112);
    const ccomplex_t IT_0211 = -conj(IT_0092);
    const ccomplex_t IT_0212 = (-2)*conj(IT_0015);
    const ccomplex_t IT_0213 = conj(IT_0053)*IT_0169 + 0.5*IT_0202*IT_0203 +
       0.5*IT_0201*IT_0204 + 0.5*conj(IT_0166)*IT_0206 + 0.5*conj(IT_0193)
      *IT_0208 + 0.5*IT_0173*IT_0210 + 0.5*IT_0200*IT_0211 + 0.5*IT_0196*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = cpow(s_12 + IT_0002 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0216*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = s_13 + s_24;
    const ccomplex_t IT_0222 = -s_14 + -s_23;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = IT_0114*IT_0223;
    const ccomplex_t IT_0225 = (-8)*IT_0224;
    const ccomplex_t IT_0226 = 2*IT_0094;
    const ccomplex_t IT_0227 = 2*IT_0045;
    const ccomplex_t IT_0228 = IT_0043 + IT_0078 + IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = 8*conj(IT_0088);
    const ccomplex_t IT_0230 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0231 = IT_0221 + IT_0230;
    const ccomplex_t IT_0232 = IT_0114*IT_0231;
    const ccomplex_t IT_0233 = (-8)*conj(IT_0053);
    const ccomplex_t IT_0234 = (-8)*conj(IT_0015);
    const ccomplex_t IT_0235 = 2*IT_0000;
    const ccomplex_t IT_0236 = IT_0043 + IT_0078 + IT_0226 + IT_0235;
    const ccomplex_t IT_0237 = 4*conj(IT_0112);
    const ccomplex_t IT_0238 = 8*conj(IT_0098);
    const ccomplex_t IT_0239 = (-4)*IT_0092;
    const ccomplex_t IT_0240 = -IT_0225;
    const ccomplex_t IT_0241 = -(32*IT_0000 + 32*IT_0043 + 32*IT_0045 + 32
      *IT_0078 + 64*IT_0094)*IT_0220 + (conj(IT_0166) + conj(IT_0193))*IT_0225 +
       (-8)*IT_0228*(IT_0088 + IT_0108 + -conj(IT_0108) + (-0.125)*IT_0229) + 8
      *IT_0232*(IT_0015 + IT_0053 + 0.125*IT_0233 + 0.125*IT_0234) + (-8)
      *IT_0236*(IT_0077 + -conj(IT_0077) + -1./2*conj(IT_0092) + IT_0098 + 1./2
      *IT_0112 + (-0.125)*IT_0237 + (-0.125)*IT_0238 + (-0.125)*IT_0239) + 
      (IT_0166 + IT_0193)*IT_0240;
    const ccomplex_t IT_0242 = 4*IT_0116;
    const ccomplex_t IT_0243 = 8*conj(IT_0112);
    const ccomplex_t IT_0244 = 8*IT_0043;
    const ccomplex_t IT_0245 = 8*conj(IT_0092);
    const ccomplex_t IT_0246 = (-4)*IT_0137;
    const ccomplex_t IT_0247 = conj(IT_0166)*IT_0203;
    const ccomplex_t IT_0248 = conj(IT_0193)*IT_0201;
    const ccomplex_t IT_0249 = IT_0078*IT_0229 + IT_0102*IT_0233 + IT_0100
      *IT_0234 + conj(IT_0135)*IT_0242 + IT_0094*(16*conj(IT_0098) + IT_0243) +
       conj(IT_0108)*IT_0244 + IT_0000*(16*conj(IT_0077) + IT_0245) + conj
      (IT_0042)*IT_0246 + IT_0247 + IT_0248;
    const ccomplex_t IT_0250 = 4*IT_0137;
    const ccomplex_t IT_0251 = (-4)*conj(IT_0053);
    const ccomplex_t IT_0252 = (-4)*IT_0116;
    const ccomplex_t IT_0253 = 4*IT_0043*conj(IT_0092) + 16*conj(IT_0088)
      *IT_0094 + 16*IT_0045*conj(IT_0108) + IT_0103*IT_0138 + IT_0078*(IT_0237 +
       IT_0238) + conj(IT_0077)*IT_0244 + IT_0247 + IT_0248 + conj(IT_0042)
      *IT_0250 + IT_0104*IT_0251 + conj(IT_0135)*IT_0252;
    const ccomplex_t IT_0254 = 2*conj(IT_0077);
    const ccomplex_t IT_0255 = conj(IT_0092) + IT_0254;
    const ccomplex_t IT_0256 = 2*conj(IT_0098);
    const ccomplex_t IT_0257 = s_12*IT_0055;
    const ccomplex_t IT_0258 = 8*IT_0257;
    const ccomplex_t IT_0259 = (-8)*IT_0056;
    const ccomplex_t IT_0260 = -IT_0047*conj(IT_0053) + conj(IT_0092)*IT_0137 
      + (-0.5)*IT_0056*IT_0205 + (-0.5)*conj(IT_0077)*IT_0246 + (-0.5)*conj
      (IT_0108)*IT_0250 + (-2)*IT_0116*(conj(IT_0088) + -1./2*conj(IT_0112) + 
      -1./2*IT_0256) + (-0.5)*conj(IT_0042)*IT_0258 + (-0.5)*conj(IT_0135)
      *IT_0259;
    const ccomplex_t IT_0261 = conj(IT_0092)*IT_0116 + 0.5*IT_0056*IT_0209 +
       0.5*IT_0047*IT_0212 + 0.5*conj(IT_0077)*IT_0242 + 0.5*conj(IT_0108)
      *IT_0252 + (-2)*IT_0137*(conj(IT_0088) + -1./2*conj(IT_0112) + -1./2
      *IT_0256) + 0.5*conj(IT_0135)*IT_0258 + 0.5*conj(IT_0042)*IT_0259;
    const ccomplex_t IT_0262 = 0.25*IT_0102;
    const ccomplex_t IT_0263 = IT_0043*conj(IT_0088) + -conj(IT_0053)*IT_0100 
      + IT_0078*conj(IT_0108) + (-0.5)*conj(IT_0042)*IT_0116 + 0.5*conj(IT_0135)
      *IT_0137 + (-0.25)*IT_0173*conj(IT_0193) + (-0.25)*conj(IT_0166)*IT_0200 +
       IT_0094*IT_0255 + IT_0000*(conj(IT_0112) + IT_0256) + IT_0138*IT_0262;
    const ccomplex_t IT_0264 = conj(IT_0193)*IT_0203;
    const ccomplex_t IT_0265 = conj(IT_0166)*IT_0201;
    const ccomplex_t IT_0266 = 16*IT_0045*conj(IT_0088) + 4*IT_0078*conj
      (IT_0092) + 16*IT_0094*conj(IT_0108) + 4*conj(IT_0042)*IT_0116 + (-4)*conj
      (IT_0135)*IT_0137 + IT_0104*IT_0138 + IT_0043*IT_0237 + IT_0103*IT_0251 +
       IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = IT_0078*conj(IT_0108) + (-0.5)*conj(IT_0042)
      *IT_0116 + 0.5*conj(IT_0135)*IT_0137 + 0.125*IT_0043*IT_0229 + 0.125
      *IT_0100*IT_0233 + 0.125*IT_0102*IT_0234 + 0.125*IT_0000*IT_0243 + 0.125
      *IT_0094*IT_0245 + 0.125*IT_0264 + 0.125*IT_0265;
    const ccomplex_t IT_0268 = 4*IT_0000*IT_0015*conj(IT_0015) + conj(IT_0042)
      *(IT_0047*IT_0054 + IT_0056*IT_0057) + conj(IT_0077)*(IT_0078*IT_0089 +
       IT_0000*IT_0093 + IT_0094*IT_0099 + IT_0100*IT_0101) + conj(IT_0098)*
      (IT_0043*IT_0089 + IT_0093*IT_0094 + IT_0000*IT_0099 + IT_0101*IT_0102) + 
      (conj(IT_0088)*IT_0103 + IT_0104*conj(IT_0108))*IT_0109 + (conj(IT_0092)
      *IT_0100 + conj(IT_0088)*IT_0104 + IT_0103*conj(IT_0108) + IT_0102*conj
      (IT_0112))*IT_0113 + 4*IT_0092*(IT_0078*conj(IT_0088) + IT_0000*conj
      (IT_0092) + IT_0043*conj(IT_0108) + IT_0094*conj(IT_0112) + 0.5*IT_0116
      *conj(IT_0135) + (-0.5)*conj(IT_0042)*IT_0137 + 0.25*IT_0100*IT_0138) +
       IT_0047*conj(IT_0135)*IT_0139 + conj(IT_0166)*(IT_0057*IT_0169 + IT_0173
      *IT_0174) + conj(IT_0193)*(IT_0054*IT_0169 + IT_0139*IT_0196 + IT_0174
      *IT_0200) + IT_0166*(IT_0201*IT_0202 + IT_0203*IT_0204 + IT_0169*IT_0205 +
       conj(IT_0193)*IT_0206 + conj(IT_0166)*IT_0208 + IT_0196*IT_0209 + IT_0200
      *IT_0210 + IT_0173*IT_0211) + 2*IT_0193*IT_0213 + IT_0220*IT_0241 +
       IT_0077*IT_0249 + IT_0108*IT_0253 + 4*IT_0053*(IT_0000*conj(IT_0053) +
       conj(IT_0015)*IT_0094 + (-0.5)*IT_0056*conj(IT_0135) + (-0.5)*conj
      (IT_0166)*IT_0196 + -IT_0102*IT_0255 + -IT_0100*(conj(IT_0112) + IT_0256))
       + (-2)*IT_0042*IT_0260 + 2*IT_0135*IT_0261 + 4*conj(IT_0053)*(IT_0015
      *IT_0094 + IT_0239*IT_0262) + 4*IT_0112*IT_0263 + IT_0088*IT_0266 + 8
      *IT_0098*IT_0267;
    return create_ccomplex_return(IT_0268);
}

