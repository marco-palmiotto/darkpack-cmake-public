#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = m_d*m_N_4;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B4)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_W4)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0009 + (-3)*IT_0014);
    const ccomplex_t IT_0016 = 0.166666666666667*IT_0015;
    const ccomplex_t IT_0017 = N_B4*e_em;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W4*e_em;
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + (-3)*IT_0022);
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = pow(m_d, 2);
    const ccomplex_t IT_0027 = pow(m_N_4, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0030 + -IT_0041;
    const ccomplex_t IT_0043 = cos(alpha);
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = sin(alpha);
    const ccomplex_t IT_0048 = IT_0006*IT_0047;
    const ccomplex_t IT_0049 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = IT_0011*IT_0043;
    const ccomplex_t IT_0052 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = IT_0011*IT_0047;
    const ccomplex_t IT_0055 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0046 + -IT_0050 + 
      -IT_0053 + IT_0056);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_d*e_em*IT_0011*IT_0031
      *IT_0033*IT_0043;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + (-2)*IT_0027 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = IT_0044*IT_0049;
    const ccomplex_t IT_0065 = IT_0045*IT_0048;
    const ccomplex_t IT_0066 = IT_0051*IT_0055;
    const ccomplex_t IT_0067 = IT_0052*IT_0054;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + IT_0065 + 
      -IT_0066 + -IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_d*e_em*IT_0011*IT_0031
      *IT_0033*IT_0047;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0063 + (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = 3*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0063 + 0.5*IT_0075;
    const ccomplex_t IT_0079 = (-3)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077 + IT_0079;
    const ccomplex_t IT_0081 = 3*conj(IT_0076);
    const ccomplex_t IT_0082 = (-3)*conj(IT_0078);
    const ccomplex_t IT_0083 = IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0006;
    const ccomplex_t IT_0085 = (-0.333333333333333)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0006;
    const ccomplex_t IT_0087 = (-0.333333333333333)*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0039*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0034;
    const ccomplex_t IT_0092 = (-0.5)*IT_0036;
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = IT_0028*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = -IT_0090 + -IT_0095;
    const ccomplex_t IT_0097 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0098 = IT_0044*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0100 = IT_0048*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0102 = IT_0051*IT_0101;
    const ccomplex_t IT_0103 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0104 = IT_0054*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0098 + -IT_0100 + 
      -IT_0102 + IT_0104);
    const ccomplex_t IT_0106 = IT_0062*IT_0105;
    const ccomplex_t IT_0107 = IT_0044*IT_0099;
    const ccomplex_t IT_0108 = IT_0048*IT_0097;
    const ccomplex_t IT_0109 = IT_0051*IT_0103;
    const ccomplex_t IT_0110 = IT_0054*IT_0101;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + IT_0108 + 
      -IT_0109 + -IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0074*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0106 + (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = 3*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0106 + 0.5*IT_0113;
    const ccomplex_t IT_0117 = (-3)*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = 3*conj(IT_0114);
    const ccomplex_t IT_0120 = (-3)*conj(IT_0116);
    const ccomplex_t IT_0121 = IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = -s_24;
    const ccomplex_t IT_0123 = s_23 + IT_0122;
    const ccomplex_t IT_0124 = IT_0000*IT_0123;
    const ccomplex_t IT_0125 = (-2)*IT_0124;
    const ccomplex_t IT_0126 = 4*IT_0003;
    const ccomplex_t IT_0127 = IT_0006*IT_0010;
    const ccomplex_t IT_0128 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = IT_0005*IT_0011;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = N_d4*conj(N_d4)*e_em;
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
      *Gamma_Z) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0136*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0010;
    const ccomplex_t IT_0147 = 0.333333333333333*IT_0146;
    const ccomplex_t IT_0148 = IT_0141*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0136*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0135;
    const ccomplex_t IT_0153 = IT_0143*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = IT_0149*IT_0152;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = (-4)*IT_0124;
    const ccomplex_t IT_0158 = 2*IT_0124;
    const ccomplex_t IT_0159 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0160 = IT_0088*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0093*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0161 + IT_0164;
    const ccomplex_t IT_0166 = 3*IT_0116;
    const ccomplex_t IT_0167 = (-3)*IT_0114;
    const ccomplex_t IT_0168 = IT_0166 + IT_0167;
    const ccomplex_t IT_0169 = 3*conj(IT_0116);
    const ccomplex_t IT_0170 = (-3)*conj(IT_0114);
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = IT_0025*IT_0162;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0038*IT_0159;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = 3*IT_0078;
    const ccomplex_t IT_0178 = (-3)*IT_0076;
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = 3*conj(IT_0078);
    const ccomplex_t IT_0181 = (-3)*conj(IT_0076);
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = (-2)*IT_0003;
    const ccomplex_t IT_0184 = (-4)*IT_0003;
    const ccomplex_t IT_0185 = (-0.5)*IT_0144;
    const ccomplex_t IT_0186 = (-0.5)*IT_0150;
    const ccomplex_t IT_0187 = conj(IT_0185) + conj(IT_0186);
    const ccomplex_t IT_0188 = IT_0185 + IT_0186;
    const ccomplex_t IT_0189 = (-0.5)*IT_0153;
    const ccomplex_t IT_0190 = (-0.5)*IT_0155;
    const ccomplex_t IT_0191 = conj(IT_0189) + conj(IT_0190);
    const ccomplex_t IT_0192 = IT_0189 + IT_0190;
    const ccomplex_t IT_0193 = 4*IT_0124;
    const ccomplex_t IT_0194 = s_12*IT_0026;
    const ccomplex_t IT_0195 = s_13*s_24;
    const ccomplex_t IT_0196 = s_12*s_34;
    const ccomplex_t IT_0197 = s_14*s_23;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0194 + IT_0195 + IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0024*IT_0092;
    const ccomplex_t IT_0201 = IT_0028*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0037*IT_0087;
    const ccomplex_t IT_0204 = IT_0039*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0202 + IT_0205;
    const ccomplex_t IT_0207 = IT_0006*IT_0032;
    const ccomplex_t IT_0208 = IT_0045*IT_0207;
    const ccomplex_t IT_0209 = sin(beta);
    const ccomplex_t IT_0210 = IT_0006*IT_0209;
    const ccomplex_t IT_0211 = IT_0049*IT_0210;
    const ccomplex_t IT_0212 = IT_0011*IT_0032;
    const ccomplex_t IT_0213 = IT_0052*IT_0212;
    const ccomplex_t IT_0214 = IT_0011*IT_0209;
    const ccomplex_t IT_0215 = IT_0055*IT_0214;
    const ccomplex_t IT_0216 = -IT_0208 + -IT_0211 + IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = m_d*e_em*IT_0011*IT_0031;
    const ccomplex_t IT_0218 = 0.5*IT_0217;
    const ccomplex_t IT_0219 = pow(m_Z, 2);
    const ccomplex_t IT_0220 = cpow((-2)*s_12 + (-2)*IT_0027 + IT_0219 + 
      -reg_prop, -1);
    const ccomplex_t IT_0221 = IT_0218*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0216*IT_0222;
    const ccomplex_t IT_0224 = IT_0049*IT_0207;
    const ccomplex_t IT_0225 = IT_0045*IT_0210;
    const ccomplex_t IT_0226 = IT_0055*IT_0212;
    const ccomplex_t IT_0227 = IT_0052*IT_0214;
    const ccomplex_t IT_0228 = -IT_0224 + IT_0225 + IT_0226 + -IT_0227;
    const ccomplex_t IT_0229 = m_d*e_em*IT_0011*IT_0031*IT_0033*IT_0209;
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + (-2)*IT_0027 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = IT_0228*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0223 + 0.5*IT_0234;
    const ccomplex_t IT_0236 = 6*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0223 + (-0.5)*IT_0234;
    const ccomplex_t IT_0238 = (-6)*IT_0237;
    const ccomplex_t IT_0239 = IT_0236 + IT_0238;
    const ccomplex_t IT_0240 = 6*conj(IT_0235);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0237);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = IT_0016*IT_0091;
    const ccomplex_t IT_0244 = IT_0028*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = IT_0035*IT_0085;
    const ccomplex_t IT_0247 = IT_0039*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0245 + IT_0248;
    const ccomplex_t IT_0250 = IT_0097*IT_0207;
    const ccomplex_t IT_0251 = IT_0099*IT_0210;
    const ccomplex_t IT_0252 = IT_0101*IT_0212;
    const ccomplex_t IT_0253 = IT_0103*IT_0214;
    const ccomplex_t IT_0254 = IT_0250 + IT_0251 + -IT_0252 + -IT_0253;
    const ccomplex_t IT_0255 = IT_0222*IT_0254;
    const ccomplex_t IT_0256 = IT_0099*IT_0207;
    const ccomplex_t IT_0257 = IT_0097*IT_0210;
    const ccomplex_t IT_0258 = IT_0103*IT_0212;
    const ccomplex_t IT_0259 = IT_0101*IT_0214;
    const ccomplex_t IT_0260 = IT_0256 + -IT_0257 + -IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = IT_0233*IT_0260;
    const ccomplex_t IT_0262 = (-0.5)*IT_0255 + (-0.5)*IT_0261;
    const ccomplex_t IT_0263 = 6*IT_0262;
    const ccomplex_t IT_0264 = 0.5*IT_0255 + 0.5*IT_0261;
    const ccomplex_t IT_0265 = (-6)*IT_0264;
    const ccomplex_t IT_0266 = IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = 6*conj(IT_0262);
    const ccomplex_t IT_0268 = (-6)*conj(IT_0264);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = -IT_0195;
    const ccomplex_t IT_0271 = IT_0194 + IT_0196 + IT_0197 + IT_0270;
    const ccomplex_t IT_0272 = IT_0162*IT_0200;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = IT_0159*IT_0203;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = -IT_0273 + -IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0237;
    const ccomplex_t IT_0278 = (-6)*IT_0235;
    const ccomplex_t IT_0279 = IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = 6*conj(IT_0237);
    const ccomplex_t IT_0281 = (-6)*conj(IT_0235);
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = IT_0162*IT_0243;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0159*IT_0246;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = -IT_0284 + -IT_0286;
    const ccomplex_t IT_0288 = 6*IT_0264;
    const ccomplex_t IT_0289 = (-6)*IT_0262;
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = 6*conj(IT_0264);
    const ccomplex_t IT_0292 = (-6)*conj(IT_0262);
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = s_23 + s_24;
    const ccomplex_t IT_0295 = IT_0000*IT_0294;
    const ccomplex_t IT_0296 = conj(IT_0154) + conj(IT_0185);
    const ccomplex_t IT_0297 = 12*IT_0237;
    const ccomplex_t IT_0298 = (-12)*IT_0235;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = IT_0154 + IT_0185;
    const ccomplex_t IT_0301 = 12*conj(IT_0237);
    const ccomplex_t IT_0302 = (-12)*conj(IT_0235);
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = conj(IT_0156) + conj(IT_0186);
    const ccomplex_t IT_0305 = 12*IT_0235;
    const ccomplex_t IT_0306 = (-12)*IT_0237;
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = IT_0156 + IT_0186;
    const ccomplex_t IT_0309 = 12*conj(IT_0235);
    const ccomplex_t IT_0310 = (-12)*conj(IT_0237);
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = conj(IT_0145) + conj(IT_0189);
    const ccomplex_t IT_0313 = 12*IT_0262;
    const ccomplex_t IT_0314 = (-12)*IT_0264;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0145 + IT_0189;
    const ccomplex_t IT_0317 = 12*conj(IT_0262);
    const ccomplex_t IT_0318 = (-12)*conj(IT_0264);
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0151) + conj(IT_0190);
    const ccomplex_t IT_0321 = 12*IT_0264;
    const ccomplex_t IT_0322 = (-12)*IT_0262;
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = IT_0151 + IT_0190;
    const ccomplex_t IT_0325 = 12*conj(IT_0264);
    const ccomplex_t IT_0326 = (-12)*conj(IT_0262);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = conj(IT_0096)*IT_0239 + IT_0096*IT_0242 + conj
      (IT_0042)*IT_0266 + IT_0042*IT_0269 + conj(IT_0176)*IT_0279 + IT_0176
      *IT_0282 + conj(IT_0165)*IT_0290 + IT_0165*IT_0293 + IT_0296*IT_0299 +
       IT_0300*IT_0303 + IT_0304*IT_0307 + IT_0308*IT_0311 + IT_0312*IT_0315 +
       IT_0316*IT_0319 + IT_0320*IT_0323 + IT_0324*IT_0327;
    const ccomplex_t IT_0329 = s_13 + s_14;
    const ccomplex_t IT_0330 = IT_0000*IT_0329;
    const ccomplex_t IT_0331 = 6*IT_0096;
    const ccomplex_t IT_0332 = 6*conj(IT_0096);
    const ccomplex_t IT_0333 = 6*IT_0176;
    const ccomplex_t IT_0334 = 6*conj(IT_0176);
    const ccomplex_t IT_0335 = conj(IT_0042)*IT_0239 + IT_0042*IT_0242 + conj
      (IT_0165)*IT_0279 + IT_0165*IT_0282 + IT_0307*IT_0312 + IT_0304*IT_0315 +
       IT_0311*IT_0316 + IT_0308*IT_0319 + IT_0299*IT_0320 + IT_0296*IT_0323 +
       IT_0303*IT_0324 + IT_0300*IT_0327 + (-6)*conj(IT_0262)*(IT_0176 + (
      -0.166666666666667)*IT_0331) + (-6)*IT_0262*(conj(IT_0176) + (
      -0.166666666666667)*IT_0332) + (-6)*conj(IT_0264)*(IT_0096 + (
      -0.166666666666667)*IT_0333) + (-6)*IT_0264*(conj(IT_0096) + (
      -0.166666666666667)*IT_0334);
    const ccomplex_t IT_0336 = -IT_0196;
    const ccomplex_t IT_0337 = IT_0195 + IT_0197 + IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0206;
    const ccomplex_t IT_0339 = 6*conj(IT_0206);
    const ccomplex_t IT_0340 = m_d*s_13*m_N_4;
    const ccomplex_t IT_0341 = 12*IT_0096;
    const ccomplex_t IT_0342 = 24*IT_0186;
    const ccomplex_t IT_0343 = 24*IT_0156;
    const ccomplex_t IT_0344 = IT_0341 + IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 12*conj(IT_0096);
    const ccomplex_t IT_0346 = 24*conj(IT_0186);
    const ccomplex_t IT_0347 = 24*conj(IT_0156);
    const ccomplex_t IT_0348 = IT_0345 + IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = 12*IT_0276;
    const ccomplex_t IT_0350 = 24*IT_0206;
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 12*conj(IT_0276);
    const ccomplex_t IT_0353 = 24*conj(IT_0206);
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = 6*IT_0276;
    const ccomplex_t IT_0356 = 12*IT_0206;
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 6*conj(IT_0276);
    const ccomplex_t IT_0359 = 12*conj(IT_0206);
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = 12*IT_0287;
    const ccomplex_t IT_0362 = 12*conj(IT_0287);
    const ccomplex_t IT_0363 = IT_0340*(conj(IT_0287)*IT_0331 + IT_0287
      *IT_0332 + conj(IT_0249)*IT_0344 + IT_0249*IT_0348 + IT_0312*IT_0351 +
       IT_0316*IT_0354 + conj(IT_0042)*IT_0357 + IT_0042*IT_0360 + IT_0304
      *IT_0361 + IT_0308*IT_0362);
    const ccomplex_t IT_0364 = m_d*s_14*m_N_4;
    const ccomplex_t IT_0365 = 12*IT_0176;
    const ccomplex_t IT_0366 = 24*IT_0154;
    const ccomplex_t IT_0367 = 24*IT_0185;
    const ccomplex_t IT_0368 = IT_0365 + IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = 12*conj(IT_0176);
    const ccomplex_t IT_0370 = 24*conj(IT_0154);
    const ccomplex_t IT_0371 = 24*conj(IT_0185);
    const ccomplex_t IT_0372 = IT_0369 + IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 24*IT_0276;
    const ccomplex_t IT_0374 = IT_0356 + IT_0373;
    const ccomplex_t IT_0375 = 24*conj(IT_0276);
    const ccomplex_t IT_0376 = IT_0359 + IT_0375;
    const ccomplex_t IT_0377 = IT_0338 + IT_0349;
    const ccomplex_t IT_0378 = IT_0339 + IT_0352;
    const ccomplex_t IT_0379 = 12*IT_0249;
    const ccomplex_t IT_0380 = 12*conj(IT_0249);
    const ccomplex_t IT_0381 = IT_0364*(conj(IT_0249)*IT_0333 + IT_0249
      *IT_0334 + conj(IT_0287)*IT_0368 + IT_0287*IT_0372 + IT_0320*IT_0374 +
       IT_0324*IT_0376 + conj(IT_0165)*IT_0377 + IT_0165*IT_0378 + IT_0296
      *IT_0379 + IT_0300*IT_0380);
    const ccomplex_t IT_0382 = m_d*s_23*m_N_4;
    const ccomplex_t IT_0383 = 12*IT_0165;
    const ccomplex_t IT_0384 = 24*IT_0190;
    const ccomplex_t IT_0385 = 24*IT_0151;
    const ccomplex_t IT_0386 = IT_0383 + IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = 12*conj(IT_0165);
    const ccomplex_t IT_0388 = 24*conj(IT_0190);
    const ccomplex_t IT_0389 = 24*conj(IT_0151);
    const ccomplex_t IT_0390 = IT_0387 + IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = conj(IT_0165)*IT_0249 + IT_0165*conj(IT_0249) +
       0.166666666666667*IT_0296*IT_0374 + 0.166666666666667*IT_0300*IT_0376 +
       0.166666666666667*conj(IT_0176)*IT_0377 + 0.166666666666667*IT_0176
      *IT_0378 + 0.166666666666667*IT_0320*IT_0379 + 0.166666666666667*IT_0324
      *IT_0380 + 0.166666666666667*conj(IT_0287)*IT_0386 + 0.166666666666667
      *IT_0287*IT_0390;
    const ccomplex_t IT_0392 = s_34*IT_0027;
    const ccomplex_t IT_0393 = 12*IT_0042;
    const ccomplex_t IT_0394 = 12*conj(IT_0042);
    const ccomplex_t IT_0395 = 24*IT_0145;
    const ccomplex_t IT_0396 = 24*conj(IT_0145);
    const ccomplex_t IT_0397 = conj(IT_0165)*IT_0331 + IT_0165*IT_0332 + conj
      (IT_0042)*IT_0333 + IT_0042*IT_0334 + IT_0320*IT_0341 + IT_0324*IT_0345 +
       conj(IT_0145)*(IT_0365 + IT_0366) + conj(IT_0189)*IT_0368 + IT_0145*
      (IT_0369 + IT_0370) + IT_0189*IT_0372 + IT_0304*IT_0386 + IT_0308*IT_0390 
      + IT_0296*IT_0393 + IT_0300*IT_0394 + conj(IT_0185)*IT_0395 + IT_0185
      *IT_0396;
    const ccomplex_t IT_0398 = s_34 + IT_0026;
    const ccomplex_t IT_0399 = s_12*IT_0398;
    const ccomplex_t IT_0400 = (-24)*IT_0237;
    const ccomplex_t IT_0401 = (-24)*conj(IT_0237);
    const ccomplex_t IT_0402 = 6*conj(IT_0076);
    const ccomplex_t IT_0403 = (-6)*conj(IT_0078);
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 6*IT_0114;
    const ccomplex_t IT_0406 = (-6)*IT_0116;
    const ccomplex_t IT_0407 = IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = (-6)*conj(IT_0116);
    const ccomplex_t IT_0409 = (-6)*IT_0078;
    const ccomplex_t IT_0410 = 48*IT_0154;
    const ccomplex_t IT_0411 = 24*IT_0176;
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = 48*IT_0185;
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = 48*conj(IT_0154);
    const ccomplex_t IT_0416 = 48*conj(IT_0151);
    const ccomplex_t IT_0417 = 48*IT_0151;
    const ccomplex_t IT_0418 = 24*IT_0165;
    const ccomplex_t IT_0419 = conj(IT_0276)*IT_0349 + conj(IT_0287)*IT_0361 +
       conj(IT_0176)*IT_0368 + conj(IT_0165)*IT_0386 + conj(IT_0154)*IT_0412 +
       conj(IT_0185)*IT_0414 + IT_0185*IT_0415 + IT_0190*IT_0416 + conj(IT_0151)
      *(IT_0417 + IT_0418) + conj(IT_0190)*(48*IT_0190 + IT_0417 + IT_0418);
    const ccomplex_t IT_0420 = -s_34;
    const ccomplex_t IT_0421 = IT_0026 + IT_0420;
    const ccomplex_t IT_0422 = IT_0027*IT_0421;
    const ccomplex_t IT_0423 = (-24)*conj(IT_0116);
    const ccomplex_t IT_0424 = 24*conj(IT_0114);
    const ccomplex_t IT_0425 = IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = 6*IT_0116;
    const ccomplex_t IT_0427 = (-6)*IT_0114;
    const ccomplex_t IT_0428 = IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = 6*conj(IT_0116);
    const ccomplex_t IT_0430 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0431 = IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 6*conj(IT_0078);
    const ccomplex_t IT_0433 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0434 = IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = 6*IT_0078;
    const ccomplex_t IT_0436 = (-24)*IT_0078;
    const ccomplex_t IT_0437 = 24*conj(IT_0116);
    const ccomplex_t IT_0438 = conj(IT_0078)*((-24)*IT_0114 + 24*IT_0116) +
       conj(IT_0076)*(24*IT_0114 + (-24)*IT_0116) + IT_0249*IT_0404 + conj
      (IT_0206)*IT_0407 + IT_0206*IT_0408 + conj(IT_0249)*IT_0409 + 6*IT_0076*
      (conj(IT_0249) + -conj(IT_0287) + 0.166666666666667*IT_0425) + conj
      (IT_0276)*IT_0428 + IT_0276*IT_0431 + IT_0287*IT_0434 + conj(IT_0287)
      *IT_0435 + 6*conj(IT_0114)*(IT_0206 + 0.166666666666667*IT_0436) + IT_0078
      *IT_0437;
    const ccomplex_t IT_0439 = s_12*IT_0421;
    const ccomplex_t IT_0440 = m_d*s_24*m_N_4;
    const ccomplex_t IT_0441 = 24*IT_0189;
    const ccomplex_t IT_0442 = IT_0395 + IT_0441;
    const ccomplex_t IT_0443 = IT_0393 + IT_0442;
    const ccomplex_t IT_0444 = 6*IT_0042;
    const ccomplex_t IT_0445 = 6*conj(IT_0042);
    const ccomplex_t IT_0446 = 24*conj(IT_0189);
    const ccomplex_t IT_0447 = IT_0304*IT_0351 + IT_0308*IT_0354 + conj
      (IT_0096)*IT_0357 + IT_0096*IT_0360 + IT_0312*IT_0361 + IT_0316*IT_0362 +
       conj(IT_0249)*IT_0443 + conj(IT_0287)*IT_0444 + IT_0287*IT_0445 + IT_0249
      *(24*conj(IT_0145) + IT_0394 + IT_0446);
    const ccomplex_t IT_0448 = conj(IT_0176)*IT_0331 + IT_0176*IT_0332 +
       IT_0296*IT_0341 + IT_0300*IT_0345 + IT_0304*IT_0368 + IT_0308*IT_0372 +
       IT_0312*IT_0383 + IT_0316*IT_0387 + IT_0320*(IT_0393 + IT_0395 + IT_0441)
       + conj(IT_0165)*IT_0444 + IT_0165*IT_0445 + IT_0324*(IT_0394 + IT_0396 +
       IT_0446);
    const ccomplex_t IT_0449 = IT_0026*IT_0027;
    const ccomplex_t IT_0450 = 48*conj(IT_0156);
    const ccomplex_t IT_0451 = 48*conj(IT_0145);
    const ccomplex_t IT_0452 = 48*IT_0145;
    const ccomplex_t IT_0453 = 48*IT_0189;
    const ccomplex_t IT_0454 = 48*IT_0156;
    const ccomplex_t IT_0455 = 24*IT_0096;
    const ccomplex_t IT_0456 = 24*conj(IT_0176);
    const ccomplex_t IT_0457 = 48*IT_0186;
    const ccomplex_t IT_0458 = conj(IT_0042)*IT_0344 + IT_0042*IT_0348 + conj
      (IT_0249)*IT_0357 + IT_0249*IT_0360 + conj(IT_0165)*IT_0368 + IT_0165
      *IT_0372 + conj(IT_0287)*IT_0377 + IT_0287*IT_0378 + conj(IT_0151)*IT_0412
       + conj(IT_0190)*IT_0414 + IT_0185*IT_0416 + conj(IT_0185)*IT_0417 + conj
      (IT_0096)*IT_0442 + IT_0316*IT_0450 + IT_0186*IT_0451 + conj(IT_0186)*
      (IT_0452 + IT_0453) + conj(IT_0145)*(IT_0454 + IT_0455) + IT_0151*(IT_0415
       + IT_0456) + IT_0190*(48*conj(IT_0185) + IT_0415 + IT_0456) + conj
      (IT_0189)*(IT_0454 + IT_0455 + IT_0457);
    const ccomplex_t IT_0459 = 24*IT_0042;
    const ccomplex_t IT_0460 = IT_0452 + IT_0459;
    const ccomplex_t IT_0461 = conj(IT_0096)*IT_0344 + conj(IT_0206)*IT_0356 +
       conj(IT_0249)*IT_0379 + conj(IT_0042)*IT_0443 + IT_0186*IT_0450 + IT_0189
      *IT_0451 + conj(IT_0156)*(IT_0454 + IT_0455) + conj(IT_0186)*(IT_0454 +
       IT_0455 + IT_0457) + conj(IT_0145)*IT_0460 + conj(IT_0189)*(IT_0453 +
       IT_0460);
    const ccomplex_t IT_0462 = (-0.166666666666667)*IT_0249;
    const ccomplex_t IT_0463 = (-0.166666666666667)*conj(IT_0249);
    const ccomplex_t IT_0464 = IT_0027*IT_0398;
    const ccomplex_t IT_0465 = IT_0282*IT_0287 + IT_0279*conj(IT_0287) + (-6)
      *conj(IT_0264)*(IT_0206 + 4*IT_0235 + (-4)*IT_0237 + (-0.166666666666667)
      *IT_0355) + (-6)*IT_0264*(conj(IT_0206) + 4*conj(IT_0235) + (-4)*conj
      (IT_0237) + (-0.166666666666667)*IT_0358) + 24*conj(IT_0262)*(IT_0235 + 
      -1./4*IT_0276 + 0.0416666666666667*IT_0338 + 0.0416666666666667*IT_0400) +
       24*IT_0262*(conj(IT_0235) + -1./4*conj(IT_0276) + 0.0416666666666667
      *IT_0339 + 0.0416666666666667*IT_0401) + (-6)*IT_0242*IT_0462 + (-6)
      *IT_0239*IT_0463;
    const ccomplex_t IT_0466 = IT_0004*(conj(IT_0042)*IT_0080 + IT_0042
      *IT_0083 + conj(IT_0096)*IT_0118 + IT_0096*IT_0121) + (IT_0083*IT_0096 +
       IT_0080*conj(IT_0096) + conj(IT_0042)*IT_0118 + IT_0042*IT_0121)*IT_0125 
      + IT_0126*(IT_0083*(IT_0145 + IT_0151) + IT_0080*(conj(IT_0145) + conj
      (IT_0151)) + IT_0121*(IT_0154 + IT_0156) + IT_0118*(conj(IT_0154) + conj
      (IT_0156))) + (IT_0121*(IT_0145 + IT_0151) + IT_0118*(conj(IT_0145) + conj
      (IT_0151)) + IT_0083*(IT_0154 + IT_0156) + IT_0080*(conj(IT_0154) + conj
      (IT_0156)))*IT_0157 + IT_0158*(conj(IT_0165)*IT_0168 + IT_0165*IT_0171 +
       conj(IT_0176)*IT_0179 + IT_0176*IT_0182) + (IT_0171*IT_0176 + IT_0168
      *conj(IT_0176) + conj(IT_0165)*IT_0179 + IT_0165*IT_0182)*IT_0183 +
       IT_0184*(IT_0168*IT_0187 + IT_0171*IT_0188 + IT_0179*IT_0191 + IT_0182
      *IT_0192) + (IT_0179*IT_0187 + IT_0182*IT_0188 + IT_0168*IT_0191 + IT_0171
      *IT_0192)*IT_0193 + IT_0199*(conj(IT_0206)*IT_0239 + IT_0206*IT_0242 +
       conj(IT_0249)*IT_0266 + IT_0249*IT_0269) + IT_0271*(conj(IT_0276)*IT_0279
       + IT_0276*IT_0282 + conj(IT_0287)*IT_0290 + IT_0287*IT_0293) + IT_0295
      *IT_0328 + IT_0330*IT_0335 + 6*IT_0337*(conj(IT_0249)*IT_0287 + IT_0249
      *conj(IT_0287) + 0.166666666666667*conj(IT_0276)*IT_0338 +
       0.166666666666667*IT_0276*IT_0339) + IT_0363 + IT_0381 + 6*IT_0382
      *IT_0391 + IT_0392*IT_0397 + IT_0399*(24*IT_0237*conj(IT_0237) + IT_0264*(
      (-24)*conj(IT_0262) + 24*conj(IT_0264)) + IT_0262*(24*conj(IT_0262) + (-24
      )*conj(IT_0264)) + conj(IT_0235)*IT_0400 + IT_0235*(24*conj(IT_0235) +
       IT_0401)) + (IT_0194 + IT_0195 + IT_0198 + IT_0336)*(IT_0276*IT_0404 +
       conj(IT_0287)*IT_0407 + 6*IT_0287*(conj(IT_0114) + 0.166666666666667
      *IT_0408) + 6*conj(IT_0276)*(IT_0076 + 0.166666666666667*IT_0409)) +
       IT_0197*IT_0419 + IT_0422*IT_0438 + (24*conj(IT_0076)*IT_0078 + IT_0076*(
      (-24)*conj(IT_0076) + 24*conj(IT_0078)) + IT_0116*IT_0425 + conj(IT_0078)
      *IT_0436 + IT_0114*((-24)*conj(IT_0114) + IT_0437))*IT_0439 + IT_0440
      *IT_0447 + IT_0194*IT_0448 + IT_0449*IT_0458 + IT_0195*IT_0461 + (IT_0194 
      + IT_0197 + IT_0270 + IT_0336)*(IT_0206*IT_0434 + (-6)*conj(IT_0206)*
      (IT_0076 + (-0.166666666666667)*IT_0435) + (-6)*IT_0431*IT_0462 + (-6)
      *IT_0428*IT_0463) + IT_0464*IT_0465;
    return create_ccomplex_return(IT_0466);
}

