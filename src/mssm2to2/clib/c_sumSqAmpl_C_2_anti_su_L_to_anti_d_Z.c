#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_d_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_d_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = s_14*s_34;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = (-6)*IT_0006;
    const ccomplex_t IT_0008 = 12*s_13;
    const ccomplex_t IT_0009 = IT_0007 + IT_0008;
    const ccomplex_t IT_0010 = pow(m_d, 2);
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = pow(m_su_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0015*IT_0018;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + 3*IT_0023);
    const ccomplex_t IT_0025 = (-0.166666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0017*IT_0025;
    const ccomplex_t IT_0027 = IT_0013*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = m_d*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0019;
    const ccomplex_t IT_0032 = 0.333333333333333*IT_0031;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0013*IT_0030*IT_0033;
    const ccomplex_t IT_0035 = IT_0029 + -IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = cos(beta);
    const ccomplex_t IT_0040 = cpow(IT_0039, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0015*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = IT_0022*IT_0043;
    const ccomplex_t IT_0046 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0047 = IT_0022*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + -IT_0045 + (-2)
      *IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = IT_0042*IT_0049;
    const ccomplex_t IT_0051 = IT_0036*IT_0037*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0015*IT_0038*IT_0040;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = IT_0022*IT_0056;
    const ccomplex_t IT_0059 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0060 = IT_0022*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + -IT_0058 + (-2)
      *IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0055*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0053*IT_0063;
    const ccomplex_t IT_0065 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0066 = IT_0020*IT_0065;
    const ccomplex_t IT_0067 = IT_0022*IT_0065;
    const ccomplex_t IT_0068 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0069 = IT_0022*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + -IT_0067 + (-2)
      *IT_0069);
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0055*IT_0071;
    const ccomplex_t IT_0073 = IT_0053*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = m_C_2*IT_0074;
    const ccomplex_t IT_0076 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0077 = IT_0020*IT_0076;
    const ccomplex_t IT_0078 = IT_0022*IT_0076;
    const ccomplex_t IT_0079 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0080 = IT_0022*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + -IT_0078 + (-2)
      *IT_0080);
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0042*IT_0082;
    const ccomplex_t IT_0084 = IT_0037*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = m_C_2*IT_0085;
    const ccomplex_t IT_0087 = IT_0051 + IT_0064 + IT_0075 + IT_0086;
    const ccomplex_t IT_0088 = IT_0035 + IT_0087;
    const ccomplex_t IT_0089 = conj(IT_0035) + conj(IT_0087);
    const ccomplex_t IT_0090 = IT_0017*IT_0082;
    const ccomplex_t IT_0091 = IT_0036*IT_0037*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0015;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0071*IT_0093;
    const ccomplex_t IT_0095 = IT_0052*IT_0053*IT_0094;
    const ccomplex_t IT_0096 = IT_0017*IT_0049;
    const ccomplex_t IT_0097 = IT_0037*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = m_C_2*IT_0098;
    const ccomplex_t IT_0100 = IT_0062*IT_0093;
    const ccomplex_t IT_0101 = IT_0053*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = m_C_2*IT_0102;
    const ccomplex_t IT_0104 = IT_0091 + IT_0095 + IT_0099 + IT_0103;
    const ccomplex_t IT_0105 = IT_0032*IT_0042;
    const ccomplex_t IT_0106 = IT_0013*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = m_d*IT_0107;
    const ccomplex_t IT_0109 = IT_0025*IT_0042;
    const ccomplex_t IT_0110 = IT_0013*IT_0030*IT_0109;
    const ccomplex_t IT_0111 = -IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = IT_0104 + IT_0111;
    const ccomplex_t IT_0113 = conj(IT_0104) + conj(IT_0111);
    const ccomplex_t IT_0114 = m_d*m_C_2;
    const ccomplex_t IT_0115 = pow(s_14, 2);
    const ccomplex_t IT_0116 = IT_0000*IT_0115;
    const ccomplex_t IT_0117 = -IT_0011;
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = IT_0114*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*-3)*e_em*(IT_0015*IT_0018 + -1.
      /3*IT_0014*IT_0019);
    const ccomplex_t IT_0121 = 0.333333333333333*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0123 = IT_0042*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0121*IT_0124;
    const ccomplex_t IT_0126 = IT_0017*IT_0122;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0121*IT_0127;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = 6*IT_0129;
    const ccomplex_t IT_0131 = (-2)*IT_0098 + (-2)*IT_0102;
    const ccomplex_t IT_0132 = (-6)*IT_0131;
    const ccomplex_t IT_0133 = IT_0130 + IT_0132;
    const ccomplex_t IT_0134 = 6*conj(IT_0129);
    const ccomplex_t IT_0135 = (-6)*conj(IT_0131);
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = (-2)*IT_0074 + (-2)*IT_0085;
    const ccomplex_t IT_0138 = 6*IT_0131;
    const ccomplex_t IT_0139 = (-6)*IT_0129;
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = 6*conj(IT_0131);
    const ccomplex_t IT_0142 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = s_13 + IT_0145;
    const ccomplex_t IT_0147 = s_13*IT_0146;
    const ccomplex_t IT_0148 = (-2)*IT_0107 + -IT_0125;
    const ccomplex_t IT_0149 = 6*IT_0137;
    const ccomplex_t IT_0150 = (-6)*IT_0125;
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = 6*conj(IT_0137);
    const ccomplex_t IT_0153 = (-6)*conj(IT_0125);
    const ccomplex_t IT_0154 = IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0028 + IT_0128;
    const ccomplex_t IT_0156 = s_13*s_14;
    const ccomplex_t IT_0157 = s_34*IT_0011;
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = IT_0074 + IT_0085;
    const ccomplex_t IT_0161 = -IT_0107;
    const ccomplex_t IT_0162 = IT_0098 + IT_0102;
    const ccomplex_t IT_0163 = m_d*IT_0118;
    const ccomplex_t IT_0164 = m_d*IT_0146;
    const ccomplex_t IT_0165 = 6*IT_0155;
    const ccomplex_t IT_0166 = 6*conj(IT_0155);
    const ccomplex_t IT_0167 = 6*IT_0148;
    const ccomplex_t IT_0168 = 6*conj(IT_0148);
    const ccomplex_t IT_0169 = pow(s_34, 2);
    const ccomplex_t IT_0170 = IT_0000*IT_0169;
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = IT_0010 + IT_0171;
    const ccomplex_t IT_0173 = m_C_2*IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0087;
    const ccomplex_t IT_0175 = 6*conj(IT_0087);
    const ccomplex_t IT_0176 = (-18)*conj(IT_0087);
    const ccomplex_t IT_0177 = s_13*IT_0118;
    const ccomplex_t IT_0178 = 6*IT_0125;
    const ccomplex_t IT_0179 = (-6)*IT_0137;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0182 = IT_0114*IT_0146;
    const ccomplex_t IT_0183 = (-6)*IT_0155;
    const ccomplex_t IT_0184 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0185 = IT_0114*IT_0172;
    const ccomplex_t IT_0186 = s_14*IT_0001;
    const ccomplex_t IT_0187 = IT_0000*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = s_14 + IT_0188;
    const ccomplex_t IT_0190 = m_d*IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0190;
    const ccomplex_t IT_0192 = 6*conj(IT_0035);
    const ccomplex_t IT_0193 = IT_0175 + IT_0192;
    const ccomplex_t IT_0194 = s_34*IT_0001;
    const ccomplex_t IT_0195 = IT_0000*IT_0194;
    const ccomplex_t IT_0196 = (-0.25)*IT_0195;
    const ccomplex_t IT_0197 = s_34 + IT_0196;
    const ccomplex_t IT_0198 = m_C_2*IT_0197;
    const ccomplex_t IT_0199 = (-4)*IT_0198;
    const ccomplex_t IT_0200 = 6*conj(IT_0104);
    const ccomplex_t IT_0201 = 6*conj(IT_0111);
    const ccomplex_t IT_0202 = IT_0200 + IT_0201;
    const ccomplex_t IT_0203 = s_14*IT_0010;
    const ccomplex_t IT_0204 = s_13*s_34;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = (-0.25)*IT_0187;
    const ccomplex_t IT_0208 = s_14 + IT_0207;
    const ccomplex_t IT_0209 = m_d*IT_0208;
    const ccomplex_t IT_0210 = (-4)*IT_0209;
    const ccomplex_t IT_0211 = (-6)*conj(IT_0035);
    const ccomplex_t IT_0212 = (-6)*conj(IT_0087);
    const ccomplex_t IT_0213 = IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0195;
    const ccomplex_t IT_0215 = s_34 + IT_0214;
    const ccomplex_t IT_0216 = m_C_2*IT_0215;
    const ccomplex_t IT_0217 = 2*IT_0216;
    const ccomplex_t IT_0218 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0219 = (-6)*conj(IT_0111);
    const ccomplex_t IT_0220 = IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = s_13*IT_0172;
    const ccomplex_t IT_0222 = (-6)*IT_0148;
    const ccomplex_t IT_0223 = m_C_2*IT_0146;
    const ccomplex_t IT_0224 = (-6)*IT_0087;
    const ccomplex_t IT_0225 = IT_0001*IT_0114;
    const ccomplex_t IT_0226 = (-3)*IT_0225;
    const ccomplex_t IT_0227 = 3*IT_0225;
    const ccomplex_t IT_0228 = pow(m_Z, 4);
    const ccomplex_t IT_0229 = s_13*IT_0228;
    const ccomplex_t IT_0230 = IT_0000*IT_0229;
    const ccomplex_t IT_0231 = 6*IT_0230;
    const ccomplex_t IT_0232 = 24*IT_0003;
    const ccomplex_t IT_0233 = (-12)*IT_0002;
    const ccomplex_t IT_0234 = IT_0231 + IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = (-6)*IT_0230;
    const ccomplex_t IT_0236 = IT_0232 + IT_0235;
    const ccomplex_t IT_0237 = 0.166666666666667*IT_0028;
    const ccomplex_t IT_0238 = 6*IT_0035;
    const ccomplex_t IT_0239 = IT_0174 + IT_0238;
    const ccomplex_t IT_0240 = 0.166666666666667*IT_0191;
    const ccomplex_t IT_0241 = 6*IT_0104;
    const ccomplex_t IT_0242 = 6*IT_0111;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 0.166666666666667*IT_0199;
    const ccomplex_t IT_0245 = 0.166666666666667*IT_0206;
    const ccomplex_t IT_0246 = (-6)*IT_0035;
    const ccomplex_t IT_0247 = IT_0224 + IT_0246;
    const ccomplex_t IT_0248 = 0.166666666666667*IT_0210;
    const ccomplex_t IT_0249 = (-6)*IT_0104;
    const ccomplex_t IT_0250 = (-6)*IT_0111;
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = 0.166666666666667*IT_0217;
    const ccomplex_t IT_0253 = 0.166666666666667*IT_0160;
    const ccomplex_t IT_0254 = 0.166666666666667*IT_0161;
    const ccomplex_t IT_0255 = IT_0222*IT_0245;
    const ccomplex_t IT_0256 = (-6)*conj(IT_0148);
    const ccomplex_t IT_0257 = IT_0206*IT_0256;
    const ccomplex_t IT_0258 = IT_0009*(IT_0088*IT_0089 + IT_0112*IT_0113) +
       IT_0119*(conj(IT_0125)*IT_0133 + IT_0125*IT_0136 + conj(IT_0137)*IT_0140 
      + IT_0137*IT_0143) + IT_0147*(conj(IT_0148)*IT_0151 + IT_0148*IT_0154 +
       IT_0143*IT_0155 + IT_0140*conj(IT_0155)) + IT_0159*(IT_0154*(IT_0160 +
       IT_0161) + IT_0151*(conj(IT_0160) + conj(IT_0161)) + IT_0143*(IT_0028 +
       IT_0162) + IT_0140*(conj(IT_0028) + conj(IT_0162))) + (IT_0089*IT_0140 +
       IT_0088*IT_0143 + IT_0113*IT_0151 + IT_0112*IT_0154)*IT_0163 + IT_0164*
      (IT_0089*IT_0165 + IT_0088*IT_0166 + IT_0113*IT_0167 + IT_0112*IT_0168) +
       IT_0173*(IT_0113*IT_0165 + IT_0112*IT_0166 + conj(IT_0035)*IT_0167 +
       IT_0035*IT_0168 + conj(IT_0148)*IT_0174 + IT_0148*IT_0175) + IT_0114*
      (conj(IT_0035)*((-18)*IT_0104 + (-18)*IT_0111) + IT_0088*((-18)*conj
      (IT_0104) + (-18)*conj(IT_0111)) + IT_0112*IT_0176) + IT_0177*(conj
      (IT_0129)*IT_0133 + IT_0129*IT_0135 + conj(IT_0131)*IT_0138 + conj(IT_0137
      )*IT_0149 + conj(IT_0125)*IT_0180 + IT_0125*IT_0181) + IT_0182*(IT_0143
      *IT_0148 + IT_0140*conj(IT_0148) + conj(IT_0137)*IT_0165 + IT_0137*IT_0166
       + conj(IT_0125)*IT_0183 + IT_0125*IT_0184) + (conj(IT_0148)*IT_0183 +
       IT_0148*IT_0184)*IT_0185 + IT_0162*(IT_0191*IT_0193 + IT_0199*IT_0202 +
       IT_0184*IT_0206) + IT_0028*(IT_0184*IT_0206 + IT_0210*IT_0213 + IT_0217
      *IT_0220) + IT_0221*(conj(IT_0155)*IT_0183 + conj(IT_0148)*IT_0222) +
       IT_0223*(IT_0113*IT_0133 + IT_0112*IT_0136 + conj(IT_0125)*IT_0174 +
       IT_0125*IT_0175 + conj(IT_0035)*IT_0180 + 6*IT_0035*(conj(IT_0125) +
       0.166666666666667*IT_0181) + IT_0137*IT_0212 + conj(IT_0137)*IT_0224) + 6
      *conj(IT_0162)*(IT_0160*IT_0226 + -IT_0161*IT_0227 + 0.166666666666667
      *IT_0162*IT_0234 + IT_0236*IT_0237 + IT_0239*IT_0240 + IT_0243*IT_0244 +
       IT_0183*IT_0245) + 6*conj(IT_0028)*(IT_0161*IT_0226 + -IT_0160*IT_0227 +
       0.166666666666667*IT_0162*IT_0236 + IT_0234*IT_0237 + IT_0183*IT_0245 +
       IT_0247*IT_0248 + IT_0251*IT_0252) + 6*conj(IT_0161)*(IT_0028*IT_0226 + 
      -IT_0162*IT_0227 + IT_0248*IT_0251 + IT_0247*IT_0252 + IT_0236*IT_0253 +
       IT_0234*IT_0254 + IT_0255) + 6*conj(IT_0160)*(IT_0162*IT_0226 + -IT_0028
      *IT_0227 + IT_0240*IT_0243 + IT_0239*IT_0244 + IT_0234*IT_0253 + IT_0236
      *IT_0254 + IT_0255) + IT_0160*(IT_0193*IT_0199 + IT_0191*IT_0202 + IT_0257
      ) + IT_0161*(IT_0213*IT_0217 + IT_0210*IT_0220 + IT_0257);
    return create_ccomplex_return(IT_0258);
}

