#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = m_s*m_N_1;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0006;
    const ccomplex_t IT_0008 = (-0.333333333333333)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0006;
    const ccomplex_t IT_0010 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = pow(m_N_1, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0017*IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0017*IT_0019*IT_0021;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0016 + -IT_0029;
    const ccomplex_t IT_0031 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0032 = cos(alpha);
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0036 = sin(alpha);
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0040 = IT_0021*IT_0032;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0043 = IT_0021*IT_0036;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0034 + IT_0038 + 
      -IT_0041 + -IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_s*e_em*IT_0017*IT_0019
      *IT_0021*IT_0036;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = IT_0033*IT_0035;
    const ccomplex_t IT_0054 = IT_0031*IT_0037;
    const ccomplex_t IT_0055 = IT_0040*IT_0042;
    const ccomplex_t IT_0056 = IT_0039*IT_0043;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + -IT_0054 + 
      -IT_0055 + IT_0056);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_s*e_em*IT_0017*IT_0019
      *IT_0021*IT_0032;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + (-2)*IT_0013 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0052 + (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = 3*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0052 + 0.5*IT_0063;
    const ccomplex_t IT_0067 = (-3)*IT_0066;
    const ccomplex_t IT_0068 = IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = 3*conj(IT_0064);
    const ccomplex_t IT_0070 = (-3)*conj(IT_0066);
    const ccomplex_t IT_0071 = IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = conj(N_B1)*e_em;
    const ccomplex_t IT_0073 = IT_0006*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em;
    const ccomplex_t IT_0076 = IT_0021*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + (-3)*IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = N_B1*e_em;
    const ccomplex_t IT_0081 = IT_0006*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W1*e_em;
    const ccomplex_t IT_0084 = IT_0021*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-3)*IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = IT_0079*IT_0087;
    const ccomplex_t IT_0089 = IT_0027*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0022;
    const ccomplex_t IT_0092 = 0.5*IT_0024;
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = IT_0014*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = -IT_0090 + -IT_0095;
    const ccomplex_t IT_0097 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0098 = IT_0033*IT_0097;
    const ccomplex_t IT_0099 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0100 = IT_0037*IT_0099;
    const ccomplex_t IT_0101 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0102 = IT_0040*IT_0101;
    const ccomplex_t IT_0103 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0104 = IT_0043*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0098 + IT_0100 + 
      -IT_0102 + -IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0051*IT_0106;
    const ccomplex_t IT_0108 = IT_0033*IT_0099;
    const ccomplex_t IT_0109 = IT_0037*IT_0097;
    const ccomplex_t IT_0110 = IT_0040*IT_0103;
    const ccomplex_t IT_0111 = IT_0043*IT_0101;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + -IT_0109 + 
      -IT_0110 + IT_0111);
    const ccomplex_t IT_0113 = IT_0062*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0107 + (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = 3*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0107 + 0.5*IT_0113;
    const ccomplex_t IT_0117 = (-3)*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = 3*conj(IT_0114);
    const ccomplex_t IT_0120 = (-3)*conj(IT_0116);
    const ccomplex_t IT_0121 = IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = -s_24;
    const ccomplex_t IT_0123 = s_23 + IT_0122;
    const ccomplex_t IT_0124 = IT_0000*IT_0123;
    const ccomplex_t IT_0125 = (-2)*IT_0124;
    const ccomplex_t IT_0126 = (-4)*IT_0124;
    const ccomplex_t IT_0127 = IT_0006*IT_0020;
    const ccomplex_t IT_0128 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = IT_0005*IT_0021;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = IT_0130*IT_0132;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0129 + IT_0131 + 
      -IT_0133 + -IT_0134);
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = e_em*IT_0127;
    const ccomplex_t IT_0138 = e_em*IT_0130;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0137 + 3*IT_0138);
    const ccomplex_t IT_0140 = (-0.166666666666667)*IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0136*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0020;
    const ccomplex_t IT_0147 = 0.333333333333333*IT_0146;
    const ccomplex_t IT_0148 = IT_0141*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0136*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = conj(IT_0145) + conj(IT_0151);
    const ccomplex_t IT_0153 = IT_0145 + IT_0151;
    const ccomplex_t IT_0154 = (-0.5)*IT_0135;
    const ccomplex_t IT_0155 = IT_0143*IT_0154;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = IT_0149*IT_0154;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = conj(IT_0156) + conj(IT_0158);
    const ccomplex_t IT_0160 = IT_0156 + IT_0158;
    const ccomplex_t IT_0161 = 4*IT_0003;
    const ccomplex_t IT_0162 = 4*IT_0124;
    const ccomplex_t IT_0163 = (-0.5)*IT_0144;
    const ccomplex_t IT_0164 = (-0.5)*IT_0150;
    const ccomplex_t IT_0165 = 3*IT_0116;
    const ccomplex_t IT_0166 = (-3)*IT_0114;
    const ccomplex_t IT_0167 = IT_0165 + IT_0166;
    const ccomplex_t IT_0168 = 3*conj(IT_0116);
    const ccomplex_t IT_0169 = (-3)*conj(IT_0114);
    const ccomplex_t IT_0170 = IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = (-0.5)*IT_0155;
    const ccomplex_t IT_0172 = (-0.5)*IT_0157;
    const ccomplex_t IT_0173 = 3*IT_0066;
    const ccomplex_t IT_0174 = (-3)*IT_0064;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = 3*conj(IT_0066);
    const ccomplex_t IT_0177 = (-3)*conj(IT_0064);
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = (-4)*IT_0003;
    const ccomplex_t IT_0180 = 2*IT_0124;
    const ccomplex_t IT_0181 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0182 = IT_0088*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0185 = IT_0093*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0183 + IT_0186;
    const ccomplex_t IT_0188 = IT_0011*IT_0184;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0026*IT_0181;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = (-2)*IT_0003;
    const ccomplex_t IT_0194 = s_12*IT_0012;
    const ccomplex_t IT_0195 = s_14*s_23;
    const ccomplex_t IT_0196 = s_12*s_34;
    const ccomplex_t IT_0197 = s_13*s_24;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0194 + IT_0195 + IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0025*IT_0087;
    const ccomplex_t IT_0201 = IT_0181*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0010*IT_0092;
    const ccomplex_t IT_0204 = IT_0184*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = -IT_0202 + -IT_0205;
    const ccomplex_t IT_0207 = IT_0006*IT_0018;
    const ccomplex_t IT_0208 = IT_0099*IT_0207;
    const ccomplex_t IT_0209 = sin(beta);
    const ccomplex_t IT_0210 = IT_0006*IT_0209;
    const ccomplex_t IT_0211 = IT_0097*IT_0210;
    const ccomplex_t IT_0212 = IT_0018*IT_0021;
    const ccomplex_t IT_0213 = IT_0103*IT_0212;
    const ccomplex_t IT_0214 = IT_0021*IT_0209;
    const ccomplex_t IT_0215 = IT_0101*IT_0214;
    const ccomplex_t IT_0216 = -IT_0208 + -IT_0211 + IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = m_s*e_em*IT_0017*IT_0021;
    const ccomplex_t IT_0218 = 0.5*IT_0217;
    const ccomplex_t IT_0219 = pow(m_Z, 2);
    const ccomplex_t IT_0220 = cpow((-2)*s_12 + (-2)*IT_0013 + IT_0219 + 
      -reg_prop, -1);
    const ccomplex_t IT_0221 = IT_0218*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0216*IT_0222;
    const ccomplex_t IT_0224 = IT_0097*IT_0207;
    const ccomplex_t IT_0225 = IT_0099*IT_0210;
    const ccomplex_t IT_0226 = IT_0101*IT_0212;
    const ccomplex_t IT_0227 = IT_0103*IT_0214;
    const ccomplex_t IT_0228 = -IT_0224 + IT_0225 + IT_0226 + -IT_0227;
    const ccomplex_t IT_0229 = m_s*e_em*IT_0017*IT_0019*IT_0021*IT_0209;
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + (-2)*IT_0013 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = IT_0228*IT_0233;
    const ccomplex_t IT_0235 = (-0.5)*IT_0223 + (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = 6*IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0223 + 0.5*IT_0234;
    const ccomplex_t IT_0238 = (-6)*IT_0237;
    const ccomplex_t IT_0239 = IT_0236 + IT_0238;
    const ccomplex_t IT_0240 = 6*conj(IT_0235);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0237);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = IT_0023*IT_0079;
    const ccomplex_t IT_0244 = IT_0181*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = IT_0008*IT_0091;
    const ccomplex_t IT_0247 = IT_0184*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = -IT_0245 + -IT_0248;
    const ccomplex_t IT_0250 = IT_0035*IT_0207;
    const ccomplex_t IT_0251 = IT_0031*IT_0210;
    const ccomplex_t IT_0252 = IT_0042*IT_0212;
    const ccomplex_t IT_0253 = IT_0039*IT_0214;
    const ccomplex_t IT_0254 = IT_0250 + IT_0251 + -IT_0252 + -IT_0253;
    const ccomplex_t IT_0255 = IT_0222*IT_0254;
    const ccomplex_t IT_0256 = IT_0031*IT_0207;
    const ccomplex_t IT_0257 = IT_0035*IT_0210;
    const ccomplex_t IT_0258 = IT_0039*IT_0212;
    const ccomplex_t IT_0259 = IT_0042*IT_0214;
    const ccomplex_t IT_0260 = IT_0256 + -IT_0257 + -IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = IT_0233*IT_0260;
    const ccomplex_t IT_0262 = 0.5*IT_0255 + 0.5*IT_0261;
    const ccomplex_t IT_0263 = 6*IT_0262;
    const ccomplex_t IT_0264 = (-0.5)*IT_0255 + (-0.5)*IT_0261;
    const ccomplex_t IT_0265 = (-6)*IT_0264;
    const ccomplex_t IT_0266 = IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = 6*conj(IT_0262);
    const ccomplex_t IT_0268 = (-6)*conj(IT_0264);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = -IT_0195;
    const ccomplex_t IT_0271 = IT_0194 + IT_0196 + IT_0197 + IT_0270;
    const ccomplex_t IT_0272 = IT_0027*IT_0200;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = IT_0014*IT_0203;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0237;
    const ccomplex_t IT_0278 = (-6)*IT_0235;
    const ccomplex_t IT_0279 = IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = 6*conj(IT_0237);
    const ccomplex_t IT_0281 = (-6)*conj(IT_0235);
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = IT_0027*IT_0243;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0014*IT_0246;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = 6*IT_0264;
    const ccomplex_t IT_0289 = (-6)*IT_0262;
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = 6*conj(IT_0264);
    const ccomplex_t IT_0292 = (-6)*conj(IT_0262);
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = s_13 + s_14;
    const ccomplex_t IT_0295 = IT_0000*IT_0294;
    const ccomplex_t IT_0296 = conj(IT_0145) + conj(IT_0171);
    const ccomplex_t IT_0297 = 12*IT_0237;
    const ccomplex_t IT_0298 = (-12)*IT_0235;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = IT_0145 + IT_0171;
    const ccomplex_t IT_0301 = 12*conj(IT_0237);
    const ccomplex_t IT_0302 = (-12)*conj(IT_0235);
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = conj(IT_0151) + conj(IT_0172);
    const ccomplex_t IT_0305 = 12*IT_0235;
    const ccomplex_t IT_0306 = (-12)*IT_0237;
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = IT_0151 + IT_0172;
    const ccomplex_t IT_0309 = 12*conj(IT_0235);
    const ccomplex_t IT_0310 = (-12)*conj(IT_0237);
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = conj(IT_0156) + conj(IT_0163);
    const ccomplex_t IT_0313 = 12*IT_0262;
    const ccomplex_t IT_0314 = (-12)*IT_0264;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0156 + IT_0163;
    const ccomplex_t IT_0317 = 12*conj(IT_0262);
    const ccomplex_t IT_0318 = (-12)*conj(IT_0264);
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0158) + conj(IT_0164);
    const ccomplex_t IT_0321 = 12*IT_0264;
    const ccomplex_t IT_0322 = (-12)*IT_0262;
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = IT_0158 + IT_0164;
    const ccomplex_t IT_0325 = 12*conj(IT_0264);
    const ccomplex_t IT_0326 = (-12)*conj(IT_0262);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = conj(IT_0192)*IT_0239 + IT_0192*IT_0242 + conj
      (IT_0187)*IT_0266 + IT_0187*IT_0269 + conj(IT_0096)*IT_0279 + IT_0096
      *IT_0282 + conj(IT_0030)*IT_0290 + IT_0030*IT_0293 + IT_0296*IT_0299 +
       IT_0300*IT_0303 + IT_0304*IT_0307 + IT_0308*IT_0311 + IT_0312*IT_0315 +
       IT_0316*IT_0319 + IT_0320*IT_0323 + IT_0324*IT_0327;
    const ccomplex_t IT_0329 = -IT_0196;
    const ccomplex_t IT_0330 = IT_0195 + IT_0197 + IT_0329;
    const ccomplex_t IT_0331 = 6*IT_0206;
    const ccomplex_t IT_0332 = 6*conj(IT_0206);
    const ccomplex_t IT_0333 = s_23 + s_24;
    const ccomplex_t IT_0334 = IT_0000*IT_0333;
    const ccomplex_t IT_0335 = 6*IT_0192;
    const ccomplex_t IT_0336 = 6*conj(IT_0192);
    const ccomplex_t IT_0337 = 6*IT_0096;
    const ccomplex_t IT_0338 = 6*conj(IT_0096);
    const ccomplex_t IT_0339 = conj(IT_0187)*IT_0239 + IT_0187*IT_0242 + conj
      (IT_0030)*IT_0279 + IT_0030*IT_0282 + IT_0307*IT_0312 + IT_0304*IT_0315 +
       IT_0311*IT_0316 + IT_0308*IT_0319 + IT_0299*IT_0320 + IT_0296*IT_0323 +
       IT_0303*IT_0324 + IT_0300*IT_0327 + (-6)*conj(IT_0262)*(IT_0096 + (
      -0.166666666666667)*IT_0335) + (-6)*IT_0262*(conj(IT_0096) + (
      -0.166666666666667)*IT_0336) + (-6)*conj(IT_0264)*(IT_0192 + (
      -0.166666666666667)*IT_0337) + (-6)*IT_0264*(conj(IT_0192) + (
      -0.166666666666667)*IT_0338);
    const ccomplex_t IT_0340 = m_s*s_24*m_N_1;
    const ccomplex_t IT_0341 = 12*IT_0096;
    const ccomplex_t IT_0342 = 24*IT_0171;
    const ccomplex_t IT_0343 = 24*IT_0145;
    const ccomplex_t IT_0344 = IT_0341 + IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 12*conj(IT_0096);
    const ccomplex_t IT_0346 = 24*conj(IT_0171);
    const ccomplex_t IT_0347 = 24*conj(IT_0145);
    const ccomplex_t IT_0348 = IT_0345 + IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = 12*IT_0206;
    const ccomplex_t IT_0350 = 24*IT_0276;
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 12*conj(IT_0206);
    const ccomplex_t IT_0353 = 24*conj(IT_0276);
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = 12*IT_0276;
    const ccomplex_t IT_0356 = IT_0331 + IT_0355;
    const ccomplex_t IT_0357 = 12*conj(IT_0276);
    const ccomplex_t IT_0358 = IT_0332 + IT_0357;
    const ccomplex_t IT_0359 = 12*IT_0249;
    const ccomplex_t IT_0360 = 12*conj(IT_0249);
    const ccomplex_t IT_0361 = IT_0340*(conj(IT_0249)*IT_0337 + IT_0249
      *IT_0338 + conj(IT_0287)*IT_0344 + IT_0287*IT_0348 + IT_0320*IT_0351 +
       IT_0324*IT_0354 + conj(IT_0030)*IT_0356 + IT_0030*IT_0358 + IT_0296
      *IT_0359 + IT_0300*IT_0360);
    const ccomplex_t IT_0362 = s_34*IT_0013;
    const ccomplex_t IT_0363 = 12*IT_0030;
    const ccomplex_t IT_0364 = 24*IT_0158;
    const ccomplex_t IT_0365 = 24*IT_0164;
    const ccomplex_t IT_0366 = IT_0363 + IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = 12*conj(IT_0030);
    const ccomplex_t IT_0368 = 24*conj(IT_0158);
    const ccomplex_t IT_0369 = 24*conj(IT_0164);
    const ccomplex_t IT_0370 = IT_0367 + IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 12*IT_0192;
    const ccomplex_t IT_0372 = 12*conj(IT_0192);
    const ccomplex_t IT_0373 = 12*IT_0187;
    const ccomplex_t IT_0374 = 12*conj(IT_0187);
    const ccomplex_t IT_0375 = 24*IT_0163;
    const ccomplex_t IT_0376 = 24*conj(IT_0163);
    const ccomplex_t IT_0377 = conj(IT_0030)*IT_0335 + IT_0030*IT_0336 + conj
      (IT_0187)*IT_0337 + IT_0187*IT_0338 + conj(IT_0163)*(IT_0341 + IT_0342) +
       conj(IT_0156)*IT_0344 + IT_0163*(IT_0345 + IT_0346) + IT_0156*IT_0348 +
       IT_0304*IT_0366 + IT_0308*IT_0370 + IT_0320*IT_0371 + IT_0324*IT_0372 +
       IT_0296*IT_0373 + IT_0300*IT_0374 + conj(IT_0145)*IT_0375 + IT_0145
      *IT_0376;
    const ccomplex_t IT_0378 = s_34 + IT_0012;
    const ccomplex_t IT_0379 = s_12*IT_0378;
    const ccomplex_t IT_0380 = (-24)*IT_0237;
    const ccomplex_t IT_0381 = (-24)*conj(IT_0237);
    const ccomplex_t IT_0382 = 6*conj(IT_0116);
    const ccomplex_t IT_0383 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0384 = IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = 6*IT_0066;
    const ccomplex_t IT_0386 = (-6)*IT_0064;
    const ccomplex_t IT_0387 = IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = (-6)*conj(IT_0064);
    const ccomplex_t IT_0389 = (-6)*IT_0114;
    const ccomplex_t IT_0390 = m_s*s_13*m_N_1;
    const ccomplex_t IT_0391 = 6*conj(IT_0249);
    const ccomplex_t IT_0392 = conj(IT_0030)*IT_0249 + 0.166666666666667
      *IT_0296*IT_0351 + 0.166666666666667*IT_0300*IT_0354 + 0.166666666666667
      *conj(IT_0096)*IT_0356 + 0.166666666666667*IT_0096*IT_0358 +
       0.166666666666667*IT_0320*IT_0359 + 0.166666666666667*IT_0324*IT_0360 +
       0.166666666666667*conj(IT_0287)*IT_0366 + 0.166666666666667*IT_0287
      *IT_0370 + 0.166666666666667*IT_0030*IT_0391;
    const ccomplex_t IT_0393 = -s_34;
    const ccomplex_t IT_0394 = IT_0012 + IT_0393;
    const ccomplex_t IT_0395 = IT_0013*IT_0394;
    const ccomplex_t IT_0396 = 6*IT_0064;
    const ccomplex_t IT_0397 = (-6)*IT_0066;
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = 6*conj(IT_0064);
    const ccomplex_t IT_0400 = (-6)*conj(IT_0066);
    const ccomplex_t IT_0401 = IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = 6*conj(IT_0114);
    const ccomplex_t IT_0403 = (-6)*conj(IT_0116);
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 6*IT_0114;
    const ccomplex_t IT_0406 = (-24)*IT_0114;
    const ccomplex_t IT_0407 = 24*conj(IT_0066);
    const ccomplex_t IT_0408 = (-24)*conj(IT_0064);
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = 24*conj(IT_0064);
    const ccomplex_t IT_0411 = (24*IT_0064 + (-24)*IT_0066)*conj(IT_0114) + ((
      -24)*IT_0064 + 24*IT_0066)*conj(IT_0116) + IT_0249*IT_0384 + conj(IT_0206)
      *IT_0387 + IT_0206*IT_0388 + conj(IT_0249)*IT_0389 + conj(IT_0276)*IT_0398
       + IT_0276*IT_0401 + IT_0287*IT_0404 + conj(IT_0287)*IT_0405 + 6*conj
      (IT_0066)*(IT_0206 + 0.166666666666667*IT_0406) + (-6)*IT_0116*(conj
      (IT_0287) + (-0.166666666666667)*IT_0391 + (-0.166666666666667)*IT_0409) +
       IT_0114*IT_0410;
    const ccomplex_t IT_0412 = s_12*IT_0394;
    const ccomplex_t IT_0413 = m_s*s_14*m_N_1;
    const ccomplex_t IT_0414 = 24*IT_0156;
    const ccomplex_t IT_0415 = IT_0373 + IT_0375 + IT_0414;
    const ccomplex_t IT_0416 = 24*IT_0206;
    const ccomplex_t IT_0417 = IT_0355 + IT_0416;
    const ccomplex_t IT_0418 = 24*conj(IT_0206);
    const ccomplex_t IT_0419 = IT_0357 + IT_0418;
    const ccomplex_t IT_0420 = 6*IT_0276;
    const ccomplex_t IT_0421 = IT_0349 + IT_0420;
    const ccomplex_t IT_0422 = 6*conj(IT_0276);
    const ccomplex_t IT_0423 = IT_0352 + IT_0422;
    const ccomplex_t IT_0424 = 6*IT_0187;
    const ccomplex_t IT_0425 = 6*conj(IT_0187);
    const ccomplex_t IT_0426 = 12*IT_0287;
    const ccomplex_t IT_0427 = 12*conj(IT_0287);
    const ccomplex_t IT_0428 = 24*conj(IT_0156);
    const ccomplex_t IT_0429 = conj(IT_0249)*IT_0415 + IT_0304*IT_0417 +
       IT_0308*IT_0419 + conj(IT_0192)*IT_0421 + IT_0192*IT_0423 + conj(IT_0287)
      *IT_0424 + IT_0287*IT_0425 + IT_0312*IT_0426 + IT_0316*IT_0427 + IT_0249*
      (24*conj(IT_0163) + IT_0374 + IT_0428);
    const ccomplex_t IT_0430 = conj(IT_0096)*IT_0335 + IT_0096*IT_0336 +
       IT_0304*IT_0344 + IT_0308*IT_0348 + IT_0312*IT_0363 + IT_0316*IT_0367 +
       IT_0296*IT_0371 + IT_0300*IT_0372 + IT_0320*(IT_0373 + IT_0375 + IT_0414)
       + conj(IT_0030)*IT_0424 + IT_0030*IT_0425 + IT_0324*(IT_0374 + IT_0376 +
       IT_0428);
    const ccomplex_t IT_0431 = m_s*s_23*m_N_1;
    const ccomplex_t IT_0432 = 24*conj(IT_0151);
    const ccomplex_t IT_0433 = 24*conj(IT_0172);
    const ccomplex_t IT_0434 = IT_0372 + IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = 24*IT_0172;
    const ccomplex_t IT_0436 = 24*IT_0151;
    const ccomplex_t IT_0437 = conj(IT_0287)*IT_0335 + IT_0287*IT_0336 +
       IT_0312*IT_0417 + IT_0316*IT_0419 + conj(IT_0187)*IT_0421 + IT_0187
      *IT_0423 + IT_0304*IT_0426 + IT_0308*IT_0427 + IT_0249*IT_0434 + conj
      (IT_0249)*(IT_0371 + IT_0435 + IT_0436);
    const ccomplex_t IT_0438 = 48*IT_0172;
    const ccomplex_t IT_0439 = 24*IT_0192;
    const ccomplex_t IT_0440 = IT_0438 + IT_0439;
    const ccomplex_t IT_0441 = 48*IT_0163;
    const ccomplex_t IT_0442 = 24*IT_0187;
    const ccomplex_t IT_0443 = IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = 48*conj(IT_0163);
    const ccomplex_t IT_0445 = 48*IT_0156;
    const ccomplex_t IT_0446 = 48*conj(IT_0172);
    const ccomplex_t IT_0447 = 24*conj(IT_0192);
    const ccomplex_t IT_0448 = conj(IT_0206)*IT_0349 + conj(IT_0249)*IT_0359 +
       conj(IT_0187)*IT_0415 + conj(IT_0192)*(IT_0371 + IT_0435) + IT_0304
      *IT_0440 + conj(IT_0163)*IT_0443 + IT_0156*IT_0444 + conj(IT_0156)*
      (IT_0443 + IT_0445) + IT_0151*(48*conj(IT_0151) + IT_0446 + IT_0447);
    const ccomplex_t IT_0449 = IT_0012*IT_0013;
    const ccomplex_t IT_0450 = 48*IT_0171;
    const ccomplex_t IT_0451 = 24*IT_0096;
    const ccomplex_t IT_0452 = IT_0450 + IT_0451;
    const ccomplex_t IT_0453 = 48*conj(IT_0164);
    const ccomplex_t IT_0454 = 48*IT_0164;
    const ccomplex_t IT_0455 = 48*IT_0158;
    const ccomplex_t IT_0456 = 48*conj(IT_0171);
    const ccomplex_t IT_0457 = 24*conj(IT_0096);
    const ccomplex_t IT_0458 = conj(IT_0030)*IT_0344 + IT_0030*IT_0348 + conj
      (IT_0287)*IT_0356 + IT_0287*IT_0358 + conj(IT_0249)*IT_0421 + IT_0249
      *IT_0423 + IT_0187*IT_0434 + conj(IT_0187)*(IT_0371 + IT_0435 + IT_0436) +
       conj(IT_0156)*(48*IT_0151 + IT_0438 + IT_0439) + conj(IT_0163)*IT_0440 +
       IT_0151*IT_0444 + conj(IT_0151)*(IT_0441 + IT_0445) + IT_0316*(IT_0446 +
       IT_0447) + conj(IT_0158)*(48*IT_0145 + IT_0450 + IT_0451) + conj(IT_0164)
      *IT_0452 + IT_0145*IT_0453 + conj(IT_0145)*(IT_0454 + IT_0455) + IT_0324*
      (IT_0456 + IT_0457);
    const ccomplex_t IT_0459 = 24*IT_0030;
    const ccomplex_t IT_0460 = IT_0454 + IT_0459;
    const ccomplex_t IT_0461 = conj(IT_0096)*(IT_0341 + IT_0342) + conj
      (IT_0276)*IT_0355 + conj(IT_0030)*IT_0366 + conj(IT_0287)*IT_0426 +
       IT_0296*IT_0452 + IT_0158*IT_0453 + IT_0145*(48*conj(IT_0145) + IT_0456 +
       IT_0457) + conj(IT_0164)*IT_0460 + conj(IT_0158)*(IT_0455 + IT_0460);
    const ccomplex_t IT_0462 = (-0.166666666666667)*IT_0249;
    const ccomplex_t IT_0463 = (-0.166666666666667)*conj(IT_0249);
    const ccomplex_t IT_0464 = IT_0013*IT_0378;
    const ccomplex_t IT_0465 = IT_0282*IT_0287 + IT_0279*conj(IT_0287) + 24
      *conj(IT_0262)*(IT_0235 + -1./4*IT_0276 + 0.0416666666666667*IT_0331 +
       0.0416666666666667*IT_0380) + 24*IT_0262*(conj(IT_0235) + -1./4*conj
      (IT_0276) + 0.0416666666666667*IT_0332 + 0.0416666666666667*IT_0381) + (-6
      )*conj(IT_0264)*(IT_0206 + 4*IT_0235 + (-4)*IT_0237 + (-0.166666666666667)
      *IT_0420) + (-6)*IT_0264*(conj(IT_0206) + 4*conj(IT_0235) + (-4)*conj
      (IT_0237) + (-0.166666666666667)*IT_0422) + (-6)*IT_0242*IT_0462 + (-6)
      *IT_0239*IT_0463;
    const ccomplex_t IT_0466 = IT_0004*(conj(IT_0030)*IT_0068 + IT_0030
      *IT_0071 + conj(IT_0096)*IT_0118 + IT_0096*IT_0121) + (IT_0071*IT_0096 +
       IT_0068*conj(IT_0096) + conj(IT_0030)*IT_0118 + IT_0030*IT_0121)*IT_0125 
      + IT_0126*(IT_0068*IT_0152 + IT_0071*IT_0153 + IT_0118*IT_0159 + IT_0121
      *IT_0160) + (IT_0118*IT_0152 + IT_0121*IT_0153 + IT_0068*IT_0159 + IT_0071
      *IT_0160)*IT_0161 + IT_0162*((conj(IT_0163) + conj(IT_0164))*IT_0167 + 
      (IT_0163 + IT_0164)*IT_0170 + (conj(IT_0171) + conj(IT_0172))*IT_0175 + 
      (IT_0171 + IT_0172)*IT_0178) + (IT_0170*(IT_0171 + IT_0172) + IT_0167*
      (conj(IT_0171) + conj(IT_0172)) + (conj(IT_0163) + conj(IT_0164))*IT_0175 
      + (IT_0163 + IT_0164)*IT_0178)*IT_0179 + IT_0180*(IT_0170*IT_0187 +
       IT_0167*conj(IT_0187) + IT_0178*IT_0192 + IT_0175*conj(IT_0192)) + 
      (IT_0178*IT_0187 + IT_0175*conj(IT_0187) + IT_0170*IT_0192 + IT_0167*conj
      (IT_0192))*IT_0193 + IT_0199*(conj(IT_0206)*IT_0239 + IT_0206*IT_0242 +
       conj(IT_0249)*IT_0266 + IT_0249*IT_0269) + IT_0271*(conj(IT_0276)*IT_0279
       + IT_0276*IT_0282 + conj(IT_0287)*IT_0290 + IT_0287*IT_0293) + IT_0295
      *IT_0328 + 6*IT_0330*(conj(IT_0249)*IT_0287 + IT_0249*conj(IT_0287) +
       0.166666666666667*conj(IT_0276)*IT_0331 + 0.166666666666667*IT_0276
      *IT_0332) + IT_0334*IT_0339 + IT_0361 + IT_0362*IT_0377 + IT_0379*(24
      *IT_0237*conj(IT_0237) + IT_0264*((-24)*conj(IT_0262) + 24*conj(IT_0264)) 
      + IT_0262*(24*conj(IT_0262) + (-24)*conj(IT_0264)) + conj(IT_0235)*IT_0380
       + IT_0235*(24*conj(IT_0235) + IT_0381)) + (IT_0194 + IT_0195 + IT_0198 +
       IT_0329)*(IT_0276*IT_0384 + conj(IT_0287)*IT_0387 + 6*IT_0287*(conj
      (IT_0066) + 0.166666666666667*IT_0388) + 6*conj(IT_0276)*(IT_0116 +
       0.166666666666667*IT_0389)) + 6*IT_0390*IT_0392 + IT_0395*IT_0411 + (24
      *IT_0114*conj(IT_0116) + IT_0116*(24*conj(IT_0114) + (-24)*conj(IT_0116)) 
      + conj(IT_0114)*IT_0406 + IT_0064*IT_0409 + IT_0066*((-24)*conj(IT_0066) +
       IT_0410))*IT_0412 + IT_0413*IT_0429 + IT_0194*IT_0430 + IT_0431*IT_0437 +
       IT_0195*IT_0448 + IT_0449*IT_0458 + IT_0197*IT_0461 + (IT_0194 + IT_0197 
      + IT_0270 + IT_0329)*(IT_0206*IT_0404 + (-6)*conj(IT_0206)*(IT_0116 + (
      -0.166666666666667)*IT_0405) + (-6)*IT_0401*IT_0462 + (-6)*IT_0398*IT_0463
      ) + IT_0464*IT_0465;
    return create_ccomplex_return(IT_0466);
}

