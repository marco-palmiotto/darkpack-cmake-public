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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0002;
    const ccomplex_t IT_0005 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0002, -1);
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009 + (-2)
      *IT_0011);
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_e, 2);
    const ccomplex_t IT_0020 = s_12*IT_0019;
    const ccomplex_t IT_0021 = e_em*IT_0004;
    const ccomplex_t IT_0022 = e_em*IT_0008;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0021 + -IT_0022);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = 8*IT_0027;
    const ccomplex_t IT_0029 = s_13*s_24;
    const ccomplex_t IT_0030 = 16*IT_0018;
    const ccomplex_t IT_0031 = IT_0015*IT_0019;
    const ccomplex_t IT_0032 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = IT_0008*IT_0032;
    const ccomplex_t IT_0035 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0036 = IT_0008*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + -IT_0034 + (-2)
      *IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0024*IT_0038;
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0007;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0007;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0015 + IT_0019 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = (-8)*IT_0018;
    const ccomplex_t IT_0052 = m_e*m_C_2;
    const ccomplex_t IT_0053 = s_23 + s_24;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = pow(m_W, -1);
    const ccomplex_t IT_0056 = cos(beta);
    const ccomplex_t IT_0057 = cpow(IT_0056, -1);
    const ccomplex_t IT_0058 = sin(beta);
    const ccomplex_t IT_0059 = m_e*e_em*IT_0007*IT_0055*IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = V_u2*U_Wm2;
    const ccomplex_t IT_0062 = IT_0056*IT_0061;
    const ccomplex_t IT_0063 = U_d2*V_Wp2;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = IT_0007*(IT_0062 + IT_0064);
    const ccomplex_t IT_0066 = 1.4142135623731*e_em*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0060*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + (-2)*IT_0015 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = pow(m_Z, 2);
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + (-2)*IT_0015 + IT_0072 + 
      -reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0056*IT_0063;
    const ccomplex_t IT_0075 = IT_0058*IT_0061;
    const ccomplex_t IT_0076 = IT_0007*(IT_0074 + -IT_0075);
    const ccomplex_t IT_0077 = 1.4142135623731*e_em*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = m_e*e_em*IT_0007*IT_0055;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0073*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0071 + -IT_0083;
    const ccomplex_t IT_0085 = (-4)*IT_0027;
    const ccomplex_t IT_0086 = s_34*IT_0015;
    const ccomplex_t IT_0087 = IT_0003*IT_0038;
    const ccomplex_t IT_0088 = IT_0016*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = 8*IT_0018;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0007*IT_0055*IT_0057;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = IT_0044*IT_0093;
    const ccomplex_t IT_0095 = IT_0047*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = s_12*s_34;
    const ccomplex_t IT_0099 = s_14*s_23;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0020 + IT_0029 + IT_0098 + IT_0100;
    const ccomplex_t IT_0102 = 2*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0007*IT_0055*IT_0057;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0045*IT_0104;
    const ccomplex_t IT_0106 = IT_0047*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_34 + IT_0019;
    const ccomplex_t IT_0110 = IT_0015*IT_0109;
    const ccomplex_t IT_0111 = 2*IT_0110;
    const ccomplex_t IT_0112 = s_13 + s_14;
    const ccomplex_t IT_0113 = IT_0052*IT_0112;
    const ccomplex_t IT_0114 = 4*IT_0113;
    const ccomplex_t IT_0115 = s_12*IT_0109;
    const ccomplex_t IT_0116 = 8*IT_0115;
    const ccomplex_t IT_0117 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0118 = IT_0056*IT_0117;
    const ccomplex_t IT_0119 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0120 = IT_0058*IT_0119;
    const ccomplex_t IT_0121 = IT_0007*(IT_0118 + IT_0120);
    const ccomplex_t IT_0122 = 1.4142135623731*e_em*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0060*IT_0123;
    const ccomplex_t IT_0125 = IT_0069*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0056*IT_0119;
    const ccomplex_t IT_0128 = IT_0058*IT_0117;
    const ccomplex_t IT_0129 = IT_0007*(IT_0127 + -IT_0128);
    const ccomplex_t IT_0130 = 1.4142135623731*e_em*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0080*IT_0131;
    const ccomplex_t IT_0133 = IT_0073*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = -IT_0126 + -IT_0134;
    const ccomplex_t IT_0136 = (-8)*IT_0110;
    const ccomplex_t IT_0137 = IT_0093*IT_0104;
    const ccomplex_t IT_0138 = IT_0047*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = 2*conj(IT_0018);
    const ccomplex_t IT_0141 = conj(IT_0139) + IT_0140;
    const ccomplex_t IT_0142 = (-4)*conj(IT_0027);
    const ccomplex_t IT_0143 = (-4)*IT_0113;
    const ccomplex_t IT_0144 = (-2)*IT_0113;
    const ccomplex_t IT_0145 = sin(alpha);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_e*e_em*IT_0007*IT_0055
      *IT_0057*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = cos(alpha);
    const ccomplex_t IT_0149 = IT_0061*IT_0148;
    const ccomplex_t IT_0150 = IT_0063*IT_0145;
    const ccomplex_t IT_0151 = IT_0007*(IT_0149 + -IT_0150);
    const ccomplex_t IT_0152 = 1.4142135623731*e_em*IT_0151;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0147*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*m_e*e_em*IT_0007*IT_0055
      *IT_0057*IT_0148;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0063*IT_0148;
    const ccomplex_t IT_0162 = IT_0061*IT_0145;
    const ccomplex_t IT_0163 = IT_0007*(IT_0161 + IT_0162);
    const ccomplex_t IT_0164 = 1.4142135623731*e_em*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0160*IT_0166;
    const ccomplex_t IT_0168 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = -IT_0158 + -IT_0170;
    const ccomplex_t IT_0172 = -s_24;
    const ccomplex_t IT_0173 = s_23 + IT_0172;
    const ccomplex_t IT_0174 = IT_0052*IT_0173;
    const ccomplex_t IT_0175 = 4*IT_0174;
    const ccomplex_t IT_0176 = conj(IT_0018) + conj(IT_0027);
    const ccomplex_t IT_0177 = -s_14;
    const ccomplex_t IT_0178 = s_13 + IT_0177;
    const ccomplex_t IT_0179 = IT_0052*IT_0178;
    const ccomplex_t IT_0180 = (-4)*IT_0179;
    const ccomplex_t IT_0181 = conj(IT_0042) + conj(IT_0090);
    const ccomplex_t IT_0182 = -IT_0098;
    const ccomplex_t IT_0183 = -IT_0029;
    const ccomplex_t IT_0184 = IT_0020 + IT_0099 + IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = 2*IT_0184;
    const ccomplex_t IT_0186 = IT_0119*IT_0148;
    const ccomplex_t IT_0187 = IT_0117*IT_0145;
    const ccomplex_t IT_0188 = IT_0007*(IT_0186 + IT_0187);
    const ccomplex_t IT_0189 = 1.4142135623731*e_em*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = IT_0160*IT_0191;
    const ccomplex_t IT_0193 = IT_0168*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0117*IT_0148;
    const ccomplex_t IT_0196 = IT_0119*IT_0145;
    const ccomplex_t IT_0197 = IT_0007*(IT_0195 + -IT_0196);
    const ccomplex_t IT_0198 = 1.4142135623731*e_em*IT_0197;
    const ccomplex_t IT_0199 = 0.5*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0147*IT_0200;
    const ccomplex_t IT_0202 = IT_0156*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = -IT_0194 + -IT_0203;
    const ccomplex_t IT_0205 = -s_34;
    const ccomplex_t IT_0206 = IT_0019 + IT_0205;
    const ccomplex_t IT_0207 = IT_0015*IT_0206;
    const ccomplex_t IT_0208 = 8*IT_0207;
    const ccomplex_t IT_0209 = s_12*IT_0206;
    const ccomplex_t IT_0210 = (-8)*IT_0209;
    const ccomplex_t IT_0211 = (-2)*IT_0207;
    const ccomplex_t IT_0212 = (-2)*IT_0174;
    const ccomplex_t IT_0213 = -conj(IT_0139);
    const ccomplex_t IT_0214 = (-2)*IT_0179;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = 2*IT_0054;
    const ccomplex_t IT_0217 = (-8)*IT_0029;
    const ccomplex_t IT_0218 = (-4)*IT_0020;
    const ccomplex_t IT_0219 = m_e*s_13*m_C_2;
    const ccomplex_t IT_0220 = (-4)*IT_0219;
    const ccomplex_t IT_0221 = m_e*s_24*m_C_2;
    const ccomplex_t IT_0222 = (-4)*IT_0221;
    const ccomplex_t IT_0223 = 2*IT_0174;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0029*conj(IT_0050) + -IT_0031*IT_0141 + 0.25
      *IT_0086*IT_0142 + 0.25*conj(IT_0084)*IT_0144 + 0.25*conj(IT_0171)*IT_0215
       + 0.25*conj(IT_0135)*IT_0216 + 0.25*conj(IT_0042)*IT_0217 + 0.25*conj
      (IT_0090)*IT_0218 + 0.25*conj(IT_0097)*IT_0220 + 0.25*conj(IT_0108)
      *IT_0222 + 0.25*conj(IT_0204)*IT_0224;
    const ccomplex_t IT_0226 = 2*IT_0179;
    const ccomplex_t IT_0227 = 8*conj(IT_0139);
    const ccomplex_t IT_0228 = 8*conj(IT_0027);
    const ccomplex_t IT_0229 = 8*IT_0020;
    const ccomplex_t IT_0230 = 8*IT_0219;
    const ccomplex_t IT_0231 = 8*IT_0221;
    const ccomplex_t IT_0232 = (-4)*IT_0054;
    const ccomplex_t IT_0233 = IT_0175*conj(IT_0204);
    const ccomplex_t IT_0234 = conj(IT_0171)*IT_0180;
    const ccomplex_t IT_0235 = 16*IT_0029*conj(IT_0042) + conj(IT_0084)
      *IT_0114 + conj(IT_0050)*IT_0217 + IT_0031*(16*conj(IT_0018) + IT_0227) +
       IT_0086*IT_0228 + conj(IT_0090)*IT_0229 + conj(IT_0097)*IT_0230 + conj
      (IT_0108)*IT_0231 + conj(IT_0135)*IT_0232 + IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = 4*conj(IT_0139);
    const ccomplex_t IT_0237 = 8*conj(IT_0018);
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = 4*IT_0054;
    const ccomplex_t IT_0240 = m_e*s_14*m_C_2;
    const ccomplex_t IT_0241 = 4*IT_0240;
    const ccomplex_t IT_0242 = m_e*s_23*m_C_2;
    const ccomplex_t IT_0243 = 4*IT_0242;
    const ccomplex_t IT_0244 = 16*conj(IT_0027)*IT_0031 + 16*conj(IT_0090)
      *IT_0099 + conj(IT_0084)*IT_0143 + conj(IT_0050)*IT_0218 + conj(IT_0042)
      *IT_0229 + IT_0233 + IT_0234 + IT_0086*IT_0238 + conj(IT_0135)*IT_0239 +
       conj(IT_0097)*IT_0241 + conj(IT_0108)*IT_0243;
    const ccomplex_t IT_0245 = (-2)*IT_0101;
    const ccomplex_t IT_0246 = (-2)*IT_0110;
    const ccomplex_t IT_0247 = IT_0135*(IT_0116*conj(IT_0135) + conj(IT_0084)
      *IT_0136 + 4*IT_0113*(conj(IT_0027) + -1./2*conj(IT_0139) + -1./2*IT_0140)
       + conj(IT_0050)*IT_0216 + conj(IT_0042)*IT_0232 + conj(IT_0090)*IT_0239 +
       conj(IT_0108)*IT_0245 + conj(IT_0097)*IT_0246);
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0248*IT_0249;
    const ccomplex_t IT_0251 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = -s_14 + -s_23;
    const ccomplex_t IT_0256 = s_13 + s_24;
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = IT_0052*IT_0257;
    const ccomplex_t IT_0259 = (-8)*IT_0258;
    const ccomplex_t IT_0260 = 2*IT_0031;
    const ccomplex_t IT_0261 = 2*IT_0099;
    const ccomplex_t IT_0262 = IT_0020 + IT_0086 + IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0264 = IT_0256 + IT_0263;
    const ccomplex_t IT_0265 = IT_0052*IT_0264;
    const ccomplex_t IT_0266 = 8*IT_0265;
    const ccomplex_t IT_0267 = (-8)*IT_0265;
    const ccomplex_t IT_0268 = -IT_0259;
    const ccomplex_t IT_0269 = 2*IT_0029;
    const ccomplex_t IT_0270 = IT_0020 + IT_0086 + IT_0260 + IT_0269;
    const ccomplex_t IT_0271 = (-2)*conj(IT_0042);
    const ccomplex_t IT_0272 = -(32*IT_0020 + 32*IT_0029 + 64*IT_0031 + 32
      *IT_0086 + 32*IT_0099)*IT_0254 + (conj(IT_0171) + conj(IT_0204))*IT_0259 +
       (-8)*(IT_0027 + IT_0090 + -conj(IT_0090) + (-0.125)*IT_0228)*IT_0262 + 
      (conj(IT_0097) + conj(IT_0108))*IT_0266 + (IT_0097 + IT_0108)*IT_0267 + 
      (IT_0171 + IT_0204)*IT_0268 + (-8)*IT_0270*(IT_0042 + -1./2*IT_0050 + 1./2
      *conj(IT_0050) + (-0.125)*IT_0051 + 1./2*IT_0139 + (-0.125)*IT_0236 + (
      -0.125)*IT_0237 + 1./2*IT_0271);
    const ccomplex_t IT_0273 = 0.25*conj(IT_0171);
    const ccomplex_t IT_0274 = 0.25*conj(IT_0090);
    const ccomplex_t IT_0275 = 0.25*conj(IT_0204);
    const ccomplex_t IT_0276 = 0.25*conj(IT_0084);
    const ccomplex_t IT_0277 = 0.25*conj(IT_0050);
    const ccomplex_t IT_0278 = 0.25*conj(IT_0135);
    const ccomplex_t IT_0279 = 0.25*conj(IT_0042);
    const ccomplex_t IT_0280 = IT_0029*conj(IT_0097) + IT_0031*conj(IT_0108) +
       0.25*IT_0221*IT_0238 + conj(IT_0027)*IT_0242 + IT_0185*IT_0273 + IT_0241
      *IT_0274 + IT_0211*IT_0275 + IT_0102*IT_0276 + IT_0220*IT_0277 + IT_0246
      *IT_0278 + IT_0230*IT_0279;
    const ccomplex_t IT_0281 = IT_0031*conj(IT_0097) + IT_0029*conj(IT_0108) +
       0.25*IT_0219*IT_0238 + conj(IT_0027)*IT_0240 + IT_0211*IT_0273 + IT_0243
      *IT_0274 + IT_0185*IT_0275 + IT_0111*IT_0276 + IT_0222*IT_0277 + IT_0245
      *IT_0278 + IT_0231*IT_0279;
    const ccomplex_t IT_0282 = conj(IT_0108)*IT_0219;
    const ccomplex_t IT_0283 = conj(IT_0097)*IT_0221;
    const ccomplex_t IT_0284 = IT_0020*conj(IT_0027) + 0.5*IT_0054*conj
      (IT_0084) + IT_0086*conj(IT_0090) + (-0.5)*IT_0113*conj(IT_0135) + IT_0029
      *IT_0141 + (-0.25)*conj(IT_0171)*IT_0212 + (-0.25)*conj(IT_0204)*IT_0226 +
       -IT_0031*(conj(IT_0050) + IT_0271) + IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = conj(IT_0171)*IT_0175;
    const ccomplex_t IT_0286 = IT_0180*conj(IT_0204);
    const ccomplex_t IT_0287 = conj(IT_0018)*(IT_0020*IT_0028 + IT_0029
      *IT_0030) + IT_0031*(IT_0030*conj(IT_0042) + conj(IT_0050)*IT_0051) +
       IT_0054*conj(IT_0084)*IT_0085 + IT_0086*(IT_0028*conj(IT_0042) + conj
      (IT_0050)*IT_0085 + conj(IT_0090)*IT_0091) + IT_0084*(conj(IT_0097)
      *IT_0102 + conj(IT_0108)*IT_0111 + conj(IT_0042)*IT_0114 + conj(IT_0084)
      *IT_0116 + conj(IT_0135)*IT_0136 + 2*IT_0054*(IT_0141 + 0.5*IT_0142) +
       conj(IT_0090)*IT_0143 + conj(IT_0050)*IT_0144) + IT_0171*(IT_0175*IT_0176
       + IT_0180*IT_0181 + conj(IT_0097)*IT_0185 + conj(IT_0204)*IT_0208 + conj
      (IT_0171)*IT_0210 + conj(IT_0108)*IT_0211 + IT_0212*IT_0213 + conj(IT_0050
      )*IT_0215) + 4*IT_0050*IT_0225 + IT_0204*(IT_0176*IT_0180 + IT_0175
      *IT_0181 + conj(IT_0108)*IT_0185 + conj(IT_0171)*IT_0208 + conj(IT_0204)
      *IT_0210 + conj(IT_0097)*IT_0211 + conj(IT_0050)*IT_0224 + IT_0213*IT_0226
      ) + IT_0042*IT_0235 + IT_0090*IT_0244 + IT_0247 + IT_0254*IT_0272 + 4
      *IT_0097*IT_0280 + 4*IT_0108*IT_0281 + IT_0091*(IT_0282 + IT_0283) + 4
      *IT_0139*IT_0284 + IT_0027*(16*IT_0031*conj(IT_0090) + 16*conj(IT_0027)
      *IT_0099 + 4*IT_0113*conj(IT_0135) + IT_0020*IT_0236 + 4*conj(IT_0108)
      *IT_0240 + 4*conj(IT_0097)*IT_0242 + IT_0285 + IT_0286) + 4*IT_0018*
      (IT_0054*conj(IT_0084) + -IT_0113*conj(IT_0135) + 0.25*IT_0029*IT_0227 +
       0.25*IT_0020*IT_0228 + 0.25*IT_0285 + 0.25*IT_0286);
    return create_ccomplex_return(IT_0287);
}

