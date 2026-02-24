#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0004;
    const ccomplex_t IT_0010 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0004, -1);
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0015;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0011 + IT_0014 + 
      -IT_0016 + -IT_0017);
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -1);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0012*IT_0022*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0012*IT_0022*IT_0024;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = pow(m_c, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0000 + IT_0030 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0000 + IT_0030 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0035 = conj(N_B1)*e_em;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W1)*e_em;
    const ccomplex_t IT_0039 = IT_0012*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + 3*IT_0040);
    const ccomplex_t IT_0042 = 0.166666666666667*IT_0041;
    const ccomplex_t IT_0043 = N_B1*e_em;
    const ccomplex_t IT_0044 = IT_0003*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_W1*e_em;
    const ccomplex_t IT_0047 = IT_0012*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 3*IT_0048);
    const ccomplex_t IT_0050 = 0.166666666666667*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = IT_0034*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0033 + -IT_0053;
    const ccomplex_t IT_0055 = s_12*IT_0030;
    const ccomplex_t IT_0056 = 12*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = pow(m_Z, 2);
    const ccomplex_t IT_0059 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0058 + 
      -reg_prop, -1);
    const ccomplex_t IT_0060 = m_c*e_em*IT_0012*IT_0022;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cos(beta);
    const ccomplex_t IT_0065 = IT_0003*IT_0064;
    const ccomplex_t IT_0066 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0069 = IT_0003*IT_0023;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = IT_0012*IT_0064;
    const ccomplex_t IT_0072 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0075 = IT_0012*IT_0023;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = -IT_0067 + -IT_0070 + IT_0073 + IT_0076;
    const ccomplex_t IT_0078 = IT_0063*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = m_c*e_em*IT_0012*IT_0022*IT_0024*IT_0064;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0065*IT_0068;
    const ccomplex_t IT_0085 = IT_0066*IT_0069;
    const ccomplex_t IT_0086 = IT_0071*IT_0074;
    const ccomplex_t IT_0087 = IT_0072*IT_0075;
    const ccomplex_t IT_0088 = -IT_0084 + IT_0085 + IT_0086 + -IT_0087;
    const ccomplex_t IT_0089 = IT_0083*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0078 + (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = m_c*m_N_1;
    const ccomplex_t IT_0092 = s_13 + s_14;
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = 12*IT_0093;
    const ccomplex_t IT_0095 = IT_0090*IT_0094;
    const ccomplex_t IT_0096 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0097 = IT_0065*IT_0096;
    const ccomplex_t IT_0098 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0099 = IT_0069*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0101 = IT_0071*IT_0100;
    const ccomplex_t IT_0102 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0103 = IT_0075*IT_0102;
    const ccomplex_t IT_0104 = IT_0097 + IT_0099 + -IT_0101 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0063*IT_0104;
    const ccomplex_t IT_0106 = IT_0065*IT_0098;
    const ccomplex_t IT_0107 = IT_0069*IT_0096;
    const ccomplex_t IT_0108 = IT_0071*IT_0102;
    const ccomplex_t IT_0109 = IT_0075*IT_0100;
    const ccomplex_t IT_0110 = IT_0106 + -IT_0107 + -IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = IT_0083*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0105 + 0.5*IT_0111;
    const ccomplex_t IT_0113 = s_23 + s_24;
    const ccomplex_t IT_0114 = IT_0091*IT_0113;
    const ccomplex_t IT_0115 = 12*IT_0114;
    const ccomplex_t IT_0116 = IT_0112*IT_0115;
    const ccomplex_t IT_0117 = (-0.5)*IT_0027;
    const ccomplex_t IT_0118 = IT_0050*IT_0117;
    const ccomplex_t IT_0119 = IT_0034*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0003;
    const ccomplex_t IT_0122 = 0.666666666666667*IT_0121;
    const ccomplex_t IT_0123 = IT_0028*IT_0122;
    const ccomplex_t IT_0124 = IT_0031*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = IT_0120 + IT_0125;
    const ccomplex_t IT_0127 = m_c*s_14*m_N_1;
    const ccomplex_t IT_0128 = 12*IT_0127;
    const ccomplex_t IT_0129 = IT_0126*IT_0128;
    const ccomplex_t IT_0130 = (-0.5)*IT_0025;
    const ccomplex_t IT_0131 = IT_0042*IT_0130;
    const ccomplex_t IT_0132 = IT_0034*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0003;
    const ccomplex_t IT_0135 = 0.666666666666667*IT_0134;
    const ccomplex_t IT_0136 = IT_0026*IT_0135;
    const ccomplex_t IT_0137 = IT_0031*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0133 + IT_0138;
    const ccomplex_t IT_0140 = m_c*s_23*m_N_1;
    const ccomplex_t IT_0141 = 12*IT_0140;
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = IT_0122*IT_0135;
    const ccomplex_t IT_0144 = IT_0031*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = IT_0117*IT_0130;
    const ccomplex_t IT_0147 = IT_0034*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = -IT_0145 + -IT_0148;
    const ccomplex_t IT_0150 = s_34*IT_0000;
    const ccomplex_t IT_0151 = 12*IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = 0.5*IT_0078 + 0.5*IT_0089;
    const ccomplex_t IT_0154 = (-12)*IT_0093;
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = e_em*IT_0009;
    const ccomplex_t IT_0157 = e_em*IT_0013;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0156 + (-3)*IT_0157);
    const ccomplex_t IT_0159 = (-0.166666666666667)*IT_0158;
    const ccomplex_t IT_0160 = IT_0001*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0019*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = IT_0056*IT_0149;
    const ccomplex_t IT_0165 = IT_0094*IT_0112;
    const ccomplex_t IT_0166 = IT_0090*IT_0115;
    const ccomplex_t IT_0167 = IT_0126*IT_0141;
    const ccomplex_t IT_0168 = IT_0054*IT_0151;
    const ccomplex_t IT_0169 = IT_0000*IT_0030;
    const ccomplex_t IT_0170 = 48*IT_0169;
    const ccomplex_t IT_0171 = IT_0021*IT_0170;
    const ccomplex_t IT_0172 = (-12)*IT_0114;
    const ccomplex_t IT_0173 = IT_0153*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0020;
    const ccomplex_t IT_0175 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0176 = IT_0051*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0029*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0177 + IT_0180;
    const ccomplex_t IT_0182 = IT_0056*IT_0181;
    const ccomplex_t IT_0183 = IT_0115*IT_0153;
    const ccomplex_t IT_0184 = IT_0131*IT_0175;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0136*IT_0178;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = -IT_0185 + -IT_0187;
    const ccomplex_t IT_0189 = m_c*s_13*m_N_1;
    const ccomplex_t IT_0190 = 12*IT_0189;
    const ccomplex_t IT_0191 = IT_0188*IT_0190;
    const ccomplex_t IT_0192 = IT_0118*IT_0175;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = IT_0123*IT_0178;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = -IT_0193 + -IT_0195;
    const ccomplex_t IT_0197 = m_c*s_24*m_N_1;
    const ccomplex_t IT_0198 = 12*IT_0197;
    const ccomplex_t IT_0199 = IT_0196*IT_0198;
    const ccomplex_t IT_0200 = 0.5*IT_0018;
    const ccomplex_t IT_0201 = IT_0008*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = s_13*s_24;
    const ccomplex_t IT_0204 = 48*IT_0203;
    const ccomplex_t IT_0205 = IT_0202*IT_0204;
    const ccomplex_t IT_0206 = IT_0112*IT_0154;
    const ccomplex_t IT_0207 = IT_0090*IT_0172;
    const ccomplex_t IT_0208 = (-0.5)*IT_0162;
    const ccomplex_t IT_0209 = IT_0094*IT_0153;
    const ccomplex_t IT_0210 = IT_0190*IT_0196;
    const ccomplex_t IT_0211 = IT_0151*IT_0181;
    const ccomplex_t IT_0212 = IT_0090*IT_0154;
    const ccomplex_t IT_0213 = IT_0112*IT_0172;
    const ccomplex_t IT_0214 = 0.5*IT_0201;
    const ccomplex_t IT_0215 = 24*conj(IT_0214);
    const ccomplex_t IT_0216 = 24*conj(IT_0021);
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = 24*conj(IT_0163);
    const ccomplex_t IT_0219 = IT_0161*IT_0200;
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = 24*conj(IT_0220);
    const ccomplex_t IT_0222 = IT_0218 + IT_0221;
    const ccomplex_t IT_0223 = s_12*s_34;
    const ccomplex_t IT_0224 = s_14*s_23;
    const ccomplex_t IT_0225 = -IT_0203;
    const ccomplex_t IT_0226 = IT_0055 + IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0223 + IT_0226;
    const ccomplex_t IT_0228 = 6*IT_0227;
    const ccomplex_t IT_0229 = -IT_0223;
    const ccomplex_t IT_0230 = IT_0203 + IT_0224 + IT_0229;
    const ccomplex_t IT_0231 = 6*IT_0230;
    const ccomplex_t IT_0232 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0233 = sin(alpha);
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*m_c*e_em*IT_0012*IT_0022
      *IT_0024*IT_0233;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = IT_0232*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = cos(alpha);
    const ccomplex_t IT_0239 = IT_0003*IT_0238;
    const ccomplex_t IT_0240 = IT_0066*IT_0239;
    const ccomplex_t IT_0241 = IT_0003*IT_0233;
    const ccomplex_t IT_0242 = IT_0068*IT_0241;
    const ccomplex_t IT_0243 = IT_0012*IT_0238;
    const ccomplex_t IT_0244 = IT_0072*IT_0243;
    const ccomplex_t IT_0245 = IT_0012*IT_0233;
    const ccomplex_t IT_0246 = IT_0074*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*(IT_0240 + -IT_0242 + 
      -IT_0244 + IT_0246);
    const ccomplex_t IT_0248 = IT_0237*IT_0247;
    const ccomplex_t IT_0249 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*m_c*e_em*IT_0012*IT_0022
      *IT_0024*IT_0238;
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0249*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = IT_0068*IT_0239;
    const ccomplex_t IT_0255 = IT_0066*IT_0241;
    const ccomplex_t IT_0256 = IT_0074*IT_0243;
    const ccomplex_t IT_0257 = IT_0072*IT_0245;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*(IT_0254 + IT_0255 + 
      -IT_0256 + -IT_0257);
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = IT_0253*IT_0259;
    const ccomplex_t IT_0261 = (-0.5)*IT_0248 + (-0.5)*IT_0260;
    const ccomplex_t IT_0262 = -IT_0224;
    const ccomplex_t IT_0263 = IT_0055 + IT_0203 + IT_0262;
    const ccomplex_t IT_0264 = IT_0229 + IT_0263;
    const ccomplex_t IT_0265 = 6*IT_0264;
    const ccomplex_t IT_0266 = 6*IT_0189;
    const ccomplex_t IT_0267 = 6*IT_0197;
    const ccomplex_t IT_0268 = 6*IT_0169;
    const ccomplex_t IT_0269 = IT_0098*IT_0239;
    const ccomplex_t IT_0270 = IT_0096*IT_0241;
    const ccomplex_t IT_0271 = IT_0102*IT_0243;
    const ccomplex_t IT_0272 = IT_0100*IT_0245;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*(IT_0269 + IT_0270 + 
      -IT_0271 + -IT_0272);
    const ccomplex_t IT_0274 = -IT_0273;
    const ccomplex_t IT_0275 = IT_0253*IT_0274;
    const ccomplex_t IT_0276 = IT_0096*IT_0239;
    const ccomplex_t IT_0277 = IT_0098*IT_0241;
    const ccomplex_t IT_0278 = IT_0100*IT_0243;
    const ccomplex_t IT_0279 = IT_0102*IT_0245;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(IT_0276 + -IT_0277 + 
      -IT_0278 + IT_0279);
    const ccomplex_t IT_0281 = IT_0237*IT_0280;
    const ccomplex_t IT_0282 = 0.5*IT_0275 + 0.5*IT_0281;
    const ccomplex_t IT_0283 = -s_34;
    const ccomplex_t IT_0284 = IT_0030 + IT_0283;
    const ccomplex_t IT_0285 = IT_0000*IT_0284;
    const ccomplex_t IT_0286 = 6*IT_0285;
    const ccomplex_t IT_0287 = s_34 + IT_0030;
    const ccomplex_t IT_0288 = IT_0000*IT_0287;
    const ccomplex_t IT_0289 = 6*IT_0288;
    const ccomplex_t IT_0290 = 12*IT_0224;
    const ccomplex_t IT_0291 = IT_0146*IT_0175;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = IT_0143*IT_0178;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = 12*IT_0169;
    const ccomplex_t IT_0297 = (-6)*IT_0227;
    const ccomplex_t IT_0298 = 0.5*IT_0248 + 0.5*IT_0260;
    const ccomplex_t IT_0299 = (-6)*IT_0264;
    const ccomplex_t IT_0300 = (-0.5)*IT_0275 + (-0.5)*IT_0281;
    const ccomplex_t IT_0301 = (-6)*IT_0285;
    const ccomplex_t IT_0302 = (-0.5)*IT_0105 + (-0.5)*IT_0111;
    const ccomplex_t IT_0303 = (-6)*IT_0288;
    const ccomplex_t IT_0304 = IT_0141*conj(IT_0181) + IT_0127*IT_0217 +
       IT_0140*IT_0222 + conj(IT_0090)*IT_0228 + conj(IT_0126)*IT_0231 + conj
      (IT_0261)*IT_0265 + conj(IT_0054)*IT_0266 + conj(IT_0149)*IT_0267 + conj
      (IT_0139)*IT_0268 + conj(IT_0282)*IT_0286 + conj(IT_0112)*IT_0289 + conj
      (IT_0196)*IT_0290 + IT_0128*conj(IT_0295) + conj(IT_0188)*IT_0296 + conj
      (IT_0153)*IT_0297 + conj(IT_0298)*IT_0299 + conj(IT_0300)*IT_0301 + conj
      (IT_0302)*IT_0303;
    const ccomplex_t IT_0305 = 6*IT_0223 + 6*IT_0263;
    const ccomplex_t IT_0306 = 6*IT_0093;
    const ccomplex_t IT_0307 = 6*IT_0114;
    const ccomplex_t IT_0308 = s_12*IT_0287;
    const ccomplex_t IT_0309 = 24*IT_0308;
    const ccomplex_t IT_0310 = 24*IT_0288;
    const ccomplex_t IT_0311 = (-24)*IT_0308;
    const ccomplex_t IT_0312 = (-24)*IT_0288;
    const ccomplex_t IT_0313 = (-6)*IT_0093;
    const ccomplex_t IT_0314 = (-6)*IT_0114;
    const ccomplex_t IT_0315 = conj(IT_0139)*IT_0289 + conj(IT_0196)*IT_0297 +
       conj(IT_0188)*IT_0303 + conj(IT_0126)*IT_0305 + conj(IT_0054)*IT_0306 +
       conj(IT_0149)*IT_0307 + conj(IT_0153)*IT_0309 + conj(IT_0302)*IT_0310 +
       conj(IT_0090)*IT_0311 + conj(IT_0112)*IT_0312 + conj(IT_0295)*IT_0313 +
       conj(IT_0181)*IT_0314;
    const ccomplex_t IT_0316 = (-6)*IT_0223 + (-6)*IT_0263;
    const ccomplex_t IT_0317 = conj(IT_0196)*IT_0228 + conj(IT_0188)*IT_0289 +
       conj(IT_0139)*IT_0303 + conj(IT_0295)*IT_0306 + conj(IT_0181)*IT_0307 +
       conj(IT_0090)*IT_0309 + conj(IT_0112)*IT_0310 + conj(IT_0153)*IT_0311 +
       conj(IT_0302)*IT_0312 + conj(IT_0054)*IT_0313 + conj(IT_0149)*IT_0314 +
       conj(IT_0126)*IT_0316;
    const ccomplex_t IT_0318 = 0.5*IT_0219;
    const ccomplex_t IT_0319 = conj(IT_0208) + conj(IT_0318);
    const ccomplex_t IT_0320 = IT_0128*conj(IT_0181) + IT_0140*IT_0217 +
       IT_0127*IT_0222 + conj(IT_0112)*IT_0228 + conj(IT_0139)*IT_0231 + conj
      (IT_0149)*IT_0266 + conj(IT_0054)*IT_0267 + conj(IT_0126)*IT_0268 + conj
      (IT_0090)*IT_0289 + conj(IT_0188)*IT_0290 + IT_0141*conj(IT_0295) + conj
      (IT_0196)*IT_0296 + IT_0286*conj(IT_0298) + conj(IT_0282)*IT_0299 +
       IT_0265*conj(IT_0300) + conj(IT_0261)*IT_0301 + IT_0297*conj(IT_0302) +
       conj(IT_0153)*IT_0303 + IT_0198*IT_0319;
    const ccomplex_t IT_0321 = 24*conj(IT_0318);
    const ccomplex_t IT_0322 = 24*conj(IT_0208);
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = 24*conj(IT_0174);
    const ccomplex_t IT_0325 = 24*conj(IT_0202);
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0226 + IT_0229;
    const ccomplex_t IT_0328 = 6*IT_0327;
    const ccomplex_t IT_0329 = 6*IT_0127;
    const ccomplex_t IT_0330 = 6*IT_0140;
    const ccomplex_t IT_0331 = 12*IT_0203;
    const ccomplex_t IT_0332 = conj(IT_0163) + conj(IT_0220);
    const ccomplex_t IT_0333 = (-6)*IT_0327;
    const ccomplex_t IT_0334 = conj(IT_0149)*IT_0190 + conj(IT_0054)*IT_0198 +
       conj(IT_0188)*IT_0231 + conj(IT_0196)*IT_0268 + conj(IT_0261)*IT_0286 +
       conj(IT_0153)*IT_0289 + conj(IT_0126)*IT_0296 + conj(IT_0298)*IT_0301 +
       conj(IT_0090)*IT_0303 + conj(IT_0302)*IT_0305 + conj(IT_0112)*IT_0316 +
       IT_0197*IT_0323 + IT_0189*IT_0326 + conj(IT_0282)*IT_0328 + conj(IT_0181)
      *IT_0329 + conj(IT_0295)*IT_0330 + conj(IT_0139)*IT_0331 + IT_0128*IT_0332
       + conj(IT_0300)*IT_0333;
    const ccomplex_t IT_0335 = conj(IT_0174) + conj(IT_0202);
    const ccomplex_t IT_0336 = conj(IT_0021) + conj(IT_0214);
    const ccomplex_t IT_0337 = conj(IT_0126)*IT_0289 + conj(IT_0188)*IT_0297 +
       conj(IT_0196)*IT_0303 + conj(IT_0139)*IT_0305 + conj(IT_0149)*IT_0306 +
       conj(IT_0054)*IT_0307 + conj(IT_0302)*IT_0309 + conj(IT_0153)*IT_0310 +
       conj(IT_0112)*IT_0311 + conj(IT_0090)*IT_0312 + conj(IT_0181)*IT_0313 +
       conj(IT_0295)*IT_0314 + IT_0115*IT_0319 + IT_0154*IT_0332 + IT_0094
      *IT_0335 + IT_0172*IT_0336;
    const ccomplex_t IT_0338 = conj(IT_0188)*IT_0228 + conj(IT_0196)*IT_0289 +
       conj(IT_0126)*IT_0303 + conj(IT_0181)*IT_0306 + conj(IT_0295)*IT_0307 +
       conj(IT_0112)*IT_0309 + conj(IT_0090)*IT_0310 + conj(IT_0302)*IT_0311 +
       conj(IT_0153)*IT_0312 + conj(IT_0149)*IT_0313 + conj(IT_0054)*IT_0314 +
       conj(IT_0139)*IT_0316;
    const ccomplex_t IT_0339 = conj(IT_0054)*IT_0190 + conj(IT_0149)*IT_0198 +
       conj(IT_0196)*IT_0231 + conj(IT_0188)*IT_0268 + conj(IT_0139)*IT_0296 +
       IT_0286*conj(IT_0300) + conj(IT_0282)*IT_0301 + IT_0289*conj(IT_0302) +
       conj(IT_0112)*IT_0303 + conj(IT_0153)*IT_0305 + conj(IT_0090)*IT_0316 +
       IT_0189*IT_0323 + IT_0197*IT_0326 + conj(IT_0298)*IT_0328 + conj(IT_0295)
      *IT_0329 + conj(IT_0181)*IT_0330 + conj(IT_0126)*IT_0331 + conj(IT_0261)
      *IT_0333;
    const ccomplex_t IT_0340 = -s_14;
    const ccomplex_t IT_0341 = s_13 + IT_0340;
    const ccomplex_t IT_0342 = IT_0091*IT_0341;
    const ccomplex_t IT_0343 = (-2)*IT_0342;
    const ccomplex_t IT_0344 = 3*IT_0343;
    const ccomplex_t IT_0345 = -s_24;
    const ccomplex_t IT_0346 = s_23 + IT_0345;
    const ccomplex_t IT_0347 = IT_0091*IT_0346;
    const ccomplex_t IT_0348 = 2*IT_0347;
    const ccomplex_t IT_0349 = 3*IT_0348;
    const ccomplex_t IT_0350 = 6*IT_0055;
    const ccomplex_t IT_0351 = 6*IT_0150;
    const ccomplex_t IT_0352 = (-3)*IT_0343;
    const ccomplex_t IT_0353 = (-3)*IT_0348;
    const ccomplex_t IT_0354 = IT_0141*conj(IT_0188) + IT_0128*conj(IT_0196) +
       IT_0169*IT_0222 + IT_0217*IT_0224 + IT_0290*conj(IT_0295) + conj(IT_0181)
      *IT_0296 + conj(IT_0090)*IT_0306 + conj(IT_0112)*IT_0307 + conj(IT_0153)
      *IT_0313 + conj(IT_0302)*IT_0314 + IT_0056*IT_0319 + conj(IT_0126)*IT_0329
       + conj(IT_0139)*IT_0330 + IT_0151*IT_0335 + conj(IT_0298)*IT_0344 + conj
      (IT_0282)*IT_0349 + conj(IT_0054)*IT_0350 + conj(IT_0149)*IT_0351 + conj
      (IT_0261)*IT_0352 + conj(IT_0300)*IT_0353;
    const ccomplex_t IT_0355 = IT_0128*conj(IT_0188) + IT_0141*conj(IT_0196) +
       IT_0169*IT_0217 + IT_0222*IT_0224 + conj(IT_0181)*IT_0290 + conj(IT_0295)
      *IT_0296 + conj(IT_0112)*IT_0306 + conj(IT_0090)*IT_0307 + conj(IT_0302)
      *IT_0313 + conj(IT_0153)*IT_0314 + conj(IT_0139)*IT_0329 + conj(IT_0126)
      *IT_0330 + conj(IT_0282)*IT_0344 + conj(IT_0298)*IT_0349 + conj(IT_0149)
      *IT_0350 + conj(IT_0054)*IT_0351 + conj(IT_0300)*IT_0352 + conj(IT_0261)
      *IT_0353;
    const ccomplex_t IT_0356 = 2*IT_0342;
    const ccomplex_t IT_0357 = 3*IT_0356;
    const ccomplex_t IT_0358 = (-2)*IT_0347;
    const ccomplex_t IT_0359 = 3*IT_0358;
    const ccomplex_t IT_0360 = (-3)*IT_0356;
    const ccomplex_t IT_0361 = (-3)*IT_0358;
    const ccomplex_t IT_0362 = conj(IT_0139)*IT_0190 + conj(IT_0126)*IT_0198 +
       conj(IT_0188)*IT_0266 + conj(IT_0196)*IT_0267 + conj(IT_0054)*IT_0296 +
       conj(IT_0302)*IT_0306 + conj(IT_0153)*IT_0307 + conj(IT_0112)*IT_0313 +
       conj(IT_0090)*IT_0314 + IT_0169*IT_0323 + IT_0203*IT_0326 + conj(IT_0149)
      *IT_0331 + conj(IT_0181)*IT_0350 + conj(IT_0295)*IT_0351 + conj(IT_0300)
      *IT_0357 + conj(IT_0261)*IT_0359 + conj(IT_0282)*IT_0360 + conj(IT_0298)
      *IT_0361;
    const ccomplex_t IT_0363 = conj(IT_0126)*IT_0190 + conj(IT_0139)*IT_0198 +
       conj(IT_0196)*IT_0266 + conj(IT_0188)*IT_0267 + conj(IT_0149)*IT_0296 +
       conj(IT_0153)*IT_0306 + conj(IT_0302)*IT_0307 + conj(IT_0090)*IT_0313 +
       conj(IT_0112)*IT_0314 + IT_0203*IT_0323 + IT_0169*IT_0326 + conj(IT_0054)
      *IT_0331 + conj(IT_0295)*IT_0350 + conj(IT_0181)*IT_0351 + conj(IT_0261)
      *IT_0357 + conj(IT_0300)*IT_0359 + conj(IT_0298)*IT_0360 + conj(IT_0282)
      *IT_0361;
    const ccomplex_t IT_0364 = (-4)*IT_0342;
    const ccomplex_t IT_0365 = 3*IT_0364;
    const ccomplex_t IT_0366 = IT_0282*IT_0365;
    const ccomplex_t IT_0367 = 4*IT_0347;
    const ccomplex_t IT_0368 = (-3)*IT_0367;
    const ccomplex_t IT_0369 = IT_0261*IT_0368;
    const ccomplex_t IT_0370 = 4*IT_0342;
    const ccomplex_t IT_0371 = 3*IT_0370;
    const ccomplex_t IT_0372 = IT_0261*IT_0371;
    const ccomplex_t IT_0373 = (-4)*IT_0347;
    const ccomplex_t IT_0374 = (-3)*IT_0373;
    const ccomplex_t IT_0375 = IT_0282*IT_0374;
    const ccomplex_t IT_0376 = conj(IT_0163) + conj(IT_0174);
    const ccomplex_t IT_0377 = 3*IT_0373;
    const ccomplex_t IT_0378 = s_12*IT_0284;
    const ccomplex_t IT_0379 = (-24)*IT_0378;
    const ccomplex_t IT_0380 = (-24)*IT_0285;
    const ccomplex_t IT_0381 = conj(IT_0021) + conj(IT_0208);
    const ccomplex_t IT_0382 = (-3)*IT_0364;
    const ccomplex_t IT_0383 = 24*IT_0378;
    const ccomplex_t IT_0384 = 24*IT_0285;
    const ccomplex_t IT_0385 = conj(IT_0196)*IT_0265 + conj(IT_0139)*IT_0286 +
       conj(IT_0188)*IT_0301 + conj(IT_0126)*IT_0333 + conj(IT_0295)*IT_0352 +
       conj(IT_0181)*IT_0353 + conj(IT_0054)*IT_0357 + conj(IT_0149)*IT_0359 +
       IT_0376*IT_0377 + conj(IT_0261)*IT_0379 + conj(IT_0282)*IT_0380 + IT_0381
      *IT_0382 + conj(IT_0298)*IT_0383 + conj(IT_0300)*IT_0384;
    const ccomplex_t IT_0386 = 3*IT_0367;
    const ccomplex_t IT_0387 = (-3)*IT_0370;
    const ccomplex_t IT_0388 = conj(IT_0196)*IT_0286 + conj(IT_0188)*IT_0299 +
       conj(IT_0126)*IT_0301 + conj(IT_0139)*IT_0328 + conj(IT_0181)*IT_0344 +
       conj(IT_0295)*IT_0349 + conj(IT_0149)*IT_0360 + conj(IT_0054)*IT_0361 +
       conj(IT_0282)*IT_0379 + conj(IT_0261)*IT_0380 + conj(IT_0300)*IT_0383 +
       conj(IT_0298)*IT_0384 + IT_0381*IT_0386 + IT_0376*IT_0387;
    const ccomplex_t IT_0389 = conj(IT_0214) + conj(IT_0318);
    const ccomplex_t IT_0390 = conj(IT_0188)*IT_0265 + conj(IT_0126)*IT_0286 +
       conj(IT_0196)*IT_0301 + conj(IT_0139)*IT_0333 + conj(IT_0181)*IT_0352 +
       conj(IT_0295)*IT_0353 + conj(IT_0149)*IT_0357 + conj(IT_0054)*IT_0359 +
       IT_0371*IT_0376 + conj(IT_0300)*IT_0379 + conj(IT_0298)*IT_0380 + IT_0368
      *IT_0381 + (conj(IT_0202) + conj(IT_0220))*IT_0382 + conj(IT_0282)*IT_0383
       + conj(IT_0261)*IT_0384 + IT_0377*IT_0389;
    const ccomplex_t IT_0391 = conj(IT_0188)*IT_0286 + conj(IT_0196)*IT_0299 +
       conj(IT_0139)*IT_0301 + conj(IT_0126)*IT_0328 + conj(IT_0295)*IT_0344 +
       conj(IT_0181)*IT_0349 + conj(IT_0054)*IT_0360 + conj(IT_0149)*IT_0361 +
       IT_0374*IT_0376 + conj(IT_0298)*IT_0379 + conj(IT_0300)*IT_0380 + IT_0365
      *IT_0381 + conj(IT_0261)*IT_0383 + conj(IT_0282)*IT_0384 + (conj(IT_0202) 
      + conj(IT_0220))*IT_0386 + IT_0387*IT_0389;
    const ccomplex_t IT_0392 = conj(IT_0298)*IT_0365;
    const ccomplex_t IT_0393 = conj(IT_0282)*IT_0386;
    const ccomplex_t IT_0394 = conj(IT_0054)*IT_0056;
    const ccomplex_t IT_0395 = conj(IT_0090)*IT_0094;
    const ccomplex_t IT_0396 = conj(IT_0112)*IT_0115;
    const ccomplex_t IT_0397 = conj(IT_0126)*IT_0128;
    const ccomplex_t IT_0398 = conj(IT_0139)*IT_0141;
    const ccomplex_t IT_0399 = conj(IT_0149)*IT_0151;
    const ccomplex_t IT_0400 = 48*IT_0224;
    const ccomplex_t IT_0401 = conj(IT_0153)*IT_0154;
    const ccomplex_t IT_0402 = IT_0172*conj(IT_0302);
    const ccomplex_t IT_0403 = conj(IT_0261)*IT_0382;
    const ccomplex_t IT_0404 = conj(IT_0300)*IT_0368;
    const ccomplex_t IT_0405 = 24*IT_0224;
    const ccomplex_t IT_0406 = 24*IT_0127;
    const ccomplex_t IT_0407 = conj(IT_0196)*IT_0406;
    const ccomplex_t IT_0408 = 24*IT_0140;
    const ccomplex_t IT_0409 = 24*IT_0169;
    const ccomplex_t IT_0410 = conj(IT_0181)*IT_0409;
    const ccomplex_t IT_0411 = IT_0055*IT_0323;
    const ccomplex_t IT_0412 = IT_0150*IT_0326;
    const ccomplex_t IT_0413 = IT_0392 + IT_0393 + IT_0394 + IT_0395 + IT_0396
       + IT_0397 + IT_0398 + IT_0399 + IT_0336*IT_0400 + IT_0401 + IT_0402 +
       IT_0403 + IT_0404 + conj(IT_0295)*IT_0405 + IT_0407 + conj(IT_0188)
      *IT_0408 + IT_0410 + IT_0411 + IT_0412;
    const ccomplex_t IT_0414 = conj(IT_0300)*IT_0377;
    const ccomplex_t IT_0415 = conj(IT_0261)*IT_0371;
    const ccomplex_t IT_0416 = conj(IT_0282)*IT_0374;
    const ccomplex_t IT_0417 = conj(IT_0298)*IT_0387;
    const ccomplex_t IT_0418 = IT_0170*IT_0332 + IT_0394 + IT_0395 + IT_0396 +
       IT_0397 + IT_0398 + IT_0399 + IT_0336*IT_0400 + IT_0401 + IT_0402 + conj
      (IT_0295)*IT_0405 + IT_0407 + conj(IT_0188)*IT_0408 + IT_0410 + IT_0411 +
       IT_0412 + IT_0414 + IT_0415 + IT_0416 + IT_0417;
    const ccomplex_t IT_0419 = conj(IT_0282)*IT_0365;
    const ccomplex_t IT_0420 = conj(IT_0298)*IT_0386;
    const ccomplex_t IT_0421 = IT_0056*conj(IT_0181);
    const ccomplex_t IT_0422 = IT_0094*conj(IT_0302);
    const ccomplex_t IT_0423 = IT_0115*conj(IT_0153);
    const ccomplex_t IT_0424 = conj(IT_0188)*IT_0190;
    const ccomplex_t IT_0425 = conj(IT_0196)*IT_0198;
    const ccomplex_t IT_0426 = IT_0151*conj(IT_0295);
    const ccomplex_t IT_0427 = IT_0170*conj(IT_0208);
    const ccomplex_t IT_0428 = conj(IT_0112)*IT_0154;
    const ccomplex_t IT_0429 = conj(IT_0090)*IT_0172;
    const ccomplex_t IT_0430 = conj(IT_0300)*IT_0382;
    const ccomplex_t IT_0431 = conj(IT_0261)*IT_0368;
    const ccomplex_t IT_0432 = 24*IT_0203;
    const ccomplex_t IT_0433 = 24*IT_0189;
    const ccomplex_t IT_0434 = 24*IT_0197;
    const ccomplex_t IT_0435 = conj(IT_0126)*IT_0434;
    const ccomplex_t IT_0436 = conj(IT_0054)*IT_0409;
    const ccomplex_t IT_0437 = IT_0150*IT_0217;
    const ccomplex_t IT_0438 = IT_0055*IT_0222 + IT_0170*conj(IT_0318) +
       IT_0419 + IT_0420 + IT_0421 + IT_0422 + IT_0423 + IT_0424 + IT_0425 +
       IT_0426 + IT_0427 + IT_0428 + IT_0429 + IT_0430 + IT_0431 + conj(IT_0149)
      *IT_0432 + conj(IT_0139)*IT_0433 + IT_0435 + IT_0436 + IT_0437;
    const ccomplex_t IT_0439 = conj(IT_0261)*IT_0377;
    const ccomplex_t IT_0440 = conj(IT_0300)*IT_0371;
    const ccomplex_t IT_0441 = conj(IT_0298)*IT_0374;
    const ccomplex_t IT_0442 = conj(IT_0282)*IT_0387;
    const ccomplex_t IT_0443 = IT_0055*IT_0222 + IT_0170*conj(IT_0318) +
       IT_0204*IT_0335 + IT_0421 + IT_0422 + IT_0423 + IT_0424 + IT_0425 +
       IT_0426 + IT_0427 + IT_0428 + IT_0429 + conj(IT_0149)*IT_0432 + conj
      (IT_0139)*IT_0433 + IT_0435 + IT_0436 + IT_0437 + IT_0439 + IT_0440 +
       IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = IT_0056*conj(IT_0149);
    const ccomplex_t IT_0445 = IT_0094*conj(IT_0112);
    const ccomplex_t IT_0446 = conj(IT_0090)*IT_0115;
    const ccomplex_t IT_0447 = IT_0128*conj(IT_0139);
    const ccomplex_t IT_0448 = conj(IT_0126)*IT_0141;
    const ccomplex_t IT_0449 = conj(IT_0054)*IT_0151;
    const ccomplex_t IT_0450 = conj(IT_0021)*IT_0170;
    const ccomplex_t IT_0451 = IT_0154*conj(IT_0302);
    const ccomplex_t IT_0452 = conj(IT_0153)*IT_0172;
    const ccomplex_t IT_0453 = conj(IT_0181)*IT_0405;
    const ccomplex_t IT_0454 = conj(IT_0196)*IT_0408;
    const ccomplex_t IT_0455 = IT_0055*IT_0326;
    const ccomplex_t IT_0456 = IT_0150*IT_0323;
    const ccomplex_t IT_0457 = IT_0332*IT_0400;
    const ccomplex_t IT_0458 = IT_0170*conj(IT_0214) + conj(IT_0188)*IT_0406 +
       conj(IT_0295)*IT_0409 + IT_0419 + IT_0420 + IT_0430 + IT_0431 + IT_0444 +
       IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 + IT_0451 +
       IT_0452 + IT_0453 + IT_0454 + IT_0455 + IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = IT_0170*conj(IT_0214) + conj(IT_0188)*IT_0406 +
       conj(IT_0295)*IT_0409 + IT_0439 + IT_0440 + IT_0441 + IT_0442 + IT_0444 +
       IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 + IT_0451 +
       IT_0452 + IT_0453 + IT_0454 + IT_0455 + IT_0456 + IT_0457;
    const ccomplex_t IT_0460 = IT_0094*conj(IT_0153);
    const ccomplex_t IT_0461 = IT_0115*conj(IT_0302);
    const ccomplex_t IT_0462 = IT_0190*conj(IT_0196);
    const ccomplex_t IT_0463 = conj(IT_0188)*IT_0198;
    const ccomplex_t IT_0464 = IT_0151*conj(IT_0181);
    const ccomplex_t IT_0465 = conj(IT_0112)*IT_0172;
    const ccomplex_t IT_0466 = conj(IT_0054)*IT_0432;
    const ccomplex_t IT_0467 = conj(IT_0126)*IT_0433;
    const ccomplex_t IT_0468 = conj(IT_0149)*IT_0409;
    const ccomplex_t IT_0469 = IT_0055*IT_0217;
    const ccomplex_t IT_0470 = IT_0056*conj(IT_0295);
    const ccomplex_t IT_0471 = IT_0170*IT_0335;
    const ccomplex_t IT_0472 = IT_0204*IT_0319;
    const ccomplex_t IT_0473 = conj(IT_0090)*IT_0154;
    const ccomplex_t IT_0474 = IT_0150*IT_0222 + IT_0392 + IT_0393 + IT_0403 +
       IT_0404 + conj(IT_0139)*IT_0434 + IT_0460 + IT_0461 + IT_0462 + IT_0463 +
       IT_0464 + IT_0465 + IT_0466 + IT_0467 + IT_0468 + IT_0469 + IT_0470 +
       IT_0471 + IT_0472 + IT_0473;
    const ccomplex_t IT_0475 = IT_0150*IT_0222 + IT_0414 + IT_0415 + IT_0416 +
       IT_0417 + conj(IT_0139)*IT_0434 + IT_0460 + IT_0461 + IT_0462 + IT_0463 +
       IT_0464 + IT_0465 + IT_0466 + IT_0467 + IT_0468 + IT_0469 + IT_0470 +
       IT_0471 + IT_0472 + IT_0473;
    const ccomplex_t IT_0476 = conj(IT_0021)*(IT_0057 + IT_0095 + IT_0116 +
       IT_0129 + IT_0142 + IT_0152 + IT_0155) + conj(IT_0163)*(IT_0164 + IT_0165
       + IT_0166 + IT_0167 + IT_0168 + IT_0171 + IT_0173) + conj(IT_0174)*
      (IT_0182 + IT_0183 + IT_0191 + IT_0199 + IT_0205 + IT_0206 + IT_0207) +
       conj(IT_0208)*(IT_0209 + IT_0210 + IT_0211 + IT_0212 + IT_0213) + IT_0196
      *IT_0304 + IT_0153*IT_0315 + IT_0090*IT_0317 + IT_0188*IT_0320 + IT_0139
      *IT_0334 + IT_0302*IT_0337 + IT_0112*IT_0338 + IT_0126*IT_0339 + IT_0295
      *IT_0354 + IT_0181*IT_0355 + IT_0149*IT_0362 + IT_0054*IT_0363 + conj
      (IT_0220)*(IT_0164 + IT_0165 + IT_0166 + IT_0167 + IT_0168 + IT_0171 +
       IT_0173 + IT_0366 + IT_0369) + conj(IT_0202)*(IT_0182 + IT_0183 + IT_0191
       + IT_0199 + IT_0205 + IT_0206 + IT_0207 + IT_0366 + IT_0369) + conj
      (IT_0214)*(IT_0057 + IT_0095 + IT_0116 + IT_0129 + IT_0142 + IT_0152 +
       IT_0155 + IT_0372 + IT_0375) + conj(IT_0318)*(IT_0209 + IT_0210 + IT_0211
       + IT_0212 + IT_0213 + IT_0372 + IT_0375) + IT_0261*IT_0385 + IT_0282
      *IT_0388 + IT_0300*IT_0390 + IT_0298*IT_0391 + IT_0021*IT_0413 + IT_0214
      *IT_0418 + IT_0202*IT_0438 + IT_0174*IT_0443 + IT_0220*IT_0458 + IT_0163
      *IT_0459 + IT_0208*IT_0474 + IT_0318*IT_0475;
    return create_ccomplex_return(IT_0476);
}

