#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_c_c(
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
    const creal_t m_N_4 = param->m_N_4;
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = m_c*m_N_4;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 4*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0007, -1);
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0014;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0010 + IT_0013 + 
      -IT_0015 + -IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = e_em*IT_0008;
    const ccomplex_t IT_0023 = e_em*IT_0012;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0022 + (-3)*IT_0023);
    const ccomplex_t IT_0025 = (-0.166666666666667)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_N_4, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0007;
    const ccomplex_t IT_0033 = (-0.666666666666667)*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = cos(alpha);
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = sin(alpha);
    const ccomplex_t IT_0043 = IT_0006*IT_0042;
    const ccomplex_t IT_0044 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = IT_0011*IT_0038;
    const ccomplex_t IT_0047 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = IT_0011*IT_0042;
    const ccomplex_t IT_0050 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0041 + -IT_0045 + 
      -IT_0048 + IT_0051);
    const ccomplex_t IT_0053 = pow(m_W, -1);
    const ccomplex_t IT_0054 = sin(beta);
    const ccomplex_t IT_0055 = cpow(IT_0054, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_c*e_em*IT_0011*IT_0042
      *IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_12 + (-2)*IT_0026 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0052*IT_0060;
    const ccomplex_t IT_0062 = IT_0039*IT_0044;
    const ccomplex_t IT_0063 = IT_0040*IT_0043;
    const ccomplex_t IT_0064 = IT_0046*IT_0050;
    const ccomplex_t IT_0065 = IT_0047*IT_0049;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + IT_0063 + 
      -IT_0064 + -IT_0065);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_c*e_em*IT_0011*IT_0038
      *IT_0053*IT_0055;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0067*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0061 + (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = 3*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0061 + 0.5*IT_0073;
    const ccomplex_t IT_0077 = (-3)*IT_0076;
    const ccomplex_t IT_0078 = IT_0075 + IT_0077;
    const ccomplex_t IT_0079 = 3*conj(IT_0074);
    const ccomplex_t IT_0080 = (-3)*conj(IT_0076);
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0020;
    const ccomplex_t IT_0083 = IT_0029*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0035*IT_0082;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0088 = IT_0039*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0090 = IT_0043*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0092 = IT_0046*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0094 = IT_0049*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0088 + -IT_0090 + 
      -IT_0092 + IT_0094);
    const ccomplex_t IT_0096 = IT_0060*IT_0095;
    const ccomplex_t IT_0097 = IT_0039*IT_0089;
    const ccomplex_t IT_0098 = IT_0043*IT_0087;
    const ccomplex_t IT_0099 = IT_0046*IT_0093;
    const ccomplex_t IT_0100 = IT_0049*IT_0091;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + IT_0098 + 
      -IT_0099 + -IT_0100);
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = IT_0072*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0096 + (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = 3*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0096 + 0.5*IT_0103;
    const ccomplex_t IT_0107 = (-3)*IT_0106;
    const ccomplex_t IT_0108 = IT_0105 + IT_0107;
    const ccomplex_t IT_0109 = 3*conj(IT_0104);
    const ccomplex_t IT_0110 = (-3)*conj(IT_0106);
    const ccomplex_t IT_0111 = IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = -s_24;
    const ccomplex_t IT_0113 = s_23 + IT_0112;
    const ccomplex_t IT_0114 = IT_0000*IT_0113;
    const ccomplex_t IT_0115 = (-4)*IT_0114;
    const ccomplex_t IT_0116 = 2*IT_0003;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0011*IT_0053*IT_0055;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0011*IT_0053*IT_0055;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = pow(m_c, 2);
    const ccomplex_t IT_0123 = cpow((-2)*s_13 + IT_0026 + IT_0122 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0026 + IT_0122 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0127 = N_B4*e_em;
    const ccomplex_t IT_0128 = IT_0006*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = N_W4*e_em;
    const ccomplex_t IT_0131 = IT_0011*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + 3*IT_0132);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B4)*e_em;
    const ccomplex_t IT_0136 = IT_0006*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W4)*e_em;
    const ccomplex_t IT_0139 = IT_0011*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + 3*IT_0140);
    const ccomplex_t IT_0142 = 0.166666666666667*IT_0141;
    const ccomplex_t IT_0143 = IT_0134*IT_0142;
    const ccomplex_t IT_0144 = IT_0126*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = -IT_0125 + -IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0006;
    const ccomplex_t IT_0148 = 0.666666666666667*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0006;
    const ccomplex_t IT_0150 = 0.666666666666667*IT_0149;
    const ccomplex_t IT_0151 = IT_0148*IT_0150;
    const ccomplex_t IT_0152 = IT_0123*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = (-0.5)*IT_0117;
    const ccomplex_t IT_0155 = (-0.5)*IT_0119;
    const ccomplex_t IT_0156 = IT_0154*IT_0155;
    const ccomplex_t IT_0157 = IT_0126*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = -IT_0153 + -IT_0158;
    const ccomplex_t IT_0160 = (-2)*IT_0114;
    const ccomplex_t IT_0161 = (-2)*IT_0003;
    const ccomplex_t IT_0162 = cpow((-2)*s_23 + IT_0026 + IT_0122 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0151*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_23 + IT_0026 + IT_0122 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0156*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0164 + IT_0167;
    const ccomplex_t IT_0169 = 3*IT_0076;
    const ccomplex_t IT_0170 = (-3)*IT_0074;
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = 3*conj(IT_0076);
    const ccomplex_t IT_0173 = (-3)*conj(IT_0074);
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = IT_0121*IT_0162;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = IT_0143*IT_0165;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = 3*IT_0106;
    const ccomplex_t IT_0181 = (-3)*IT_0104;
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = 3*conj(IT_0106);
    const ccomplex_t IT_0184 = (-3)*conj(IT_0104);
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = 2*IT_0114;
    const ccomplex_t IT_0187 = 4*IT_0114;
    const ccomplex_t IT_0188 = (-0.5)*IT_0030;
    const ccomplex_t IT_0189 = (-0.5)*IT_0036;
    const ccomplex_t IT_0190 = conj(IT_0188) + conj(IT_0189);
    const ccomplex_t IT_0191 = IT_0188 + IT_0189;
    const ccomplex_t IT_0192 = (-0.5)*IT_0083;
    const ccomplex_t IT_0193 = (-0.5)*IT_0085;
    const ccomplex_t IT_0194 = conj(IT_0192) + conj(IT_0193);
    const ccomplex_t IT_0195 = IT_0192 + IT_0193;
    const ccomplex_t IT_0196 = (-4)*IT_0003;
    const ccomplex_t IT_0197 = s_12*IT_0122;
    const ccomplex_t IT_0198 = s_13*s_24;
    const ccomplex_t IT_0199 = s_12*s_34;
    const ccomplex_t IT_0200 = s_14*s_23;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0197 + IT_0198 + IT_0199 + IT_0201;
    const ccomplex_t IT_0203 = IT_0134*IT_0155;
    const ccomplex_t IT_0204 = IT_0126*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0120*IT_0150;
    const ccomplex_t IT_0207 = IT_0123*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0205 + IT_0208;
    const ccomplex_t IT_0210 = cos(beta);
    const ccomplex_t IT_0211 = IT_0006*IT_0210;
    const ccomplex_t IT_0212 = IT_0040*IT_0211;
    const ccomplex_t IT_0213 = IT_0006*IT_0054;
    const ccomplex_t IT_0214 = IT_0044*IT_0213;
    const ccomplex_t IT_0215 = IT_0011*IT_0210;
    const ccomplex_t IT_0216 = IT_0047*IT_0215;
    const ccomplex_t IT_0217 = IT_0011*IT_0054;
    const ccomplex_t IT_0218 = IT_0050*IT_0217;
    const ccomplex_t IT_0219 = -IT_0212 + -IT_0214 + IT_0216 + IT_0218;
    const ccomplex_t IT_0220 = m_c*e_em*IT_0011*IT_0053;
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = pow(m_Z, 2);
    const ccomplex_t IT_0223 = cpow((-2)*s_12 + (-2)*IT_0026 + IT_0222 + 
      -reg_prop, -1);
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = IT_0219*IT_0225;
    const ccomplex_t IT_0227 = IT_0044*IT_0211;
    const ccomplex_t IT_0228 = IT_0040*IT_0213;
    const ccomplex_t IT_0229 = IT_0050*IT_0215;
    const ccomplex_t IT_0230 = IT_0047*IT_0217;
    const ccomplex_t IT_0231 = -IT_0227 + IT_0228 + IT_0229 + -IT_0230;
    const ccomplex_t IT_0232 = m_c*e_em*IT_0011*IT_0053*IT_0055*IT_0210;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = cpow((-2)*s_12 + (-2)*IT_0026 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0235 = IT_0233*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0231*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0226 + 0.5*IT_0237;
    const ccomplex_t IT_0239 = 6*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0226 + (-0.5)*IT_0237;
    const ccomplex_t IT_0241 = (-6)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = 6*conj(IT_0238);
    const ccomplex_t IT_0244 = (-6)*conj(IT_0240);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0118*IT_0148;
    const ccomplex_t IT_0247 = IT_0123*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0142*IT_0154;
    const ccomplex_t IT_0250 = IT_0126*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0248 + IT_0251;
    const ccomplex_t IT_0253 = IT_0087*IT_0211;
    const ccomplex_t IT_0254 = IT_0089*IT_0213;
    const ccomplex_t IT_0255 = IT_0091*IT_0215;
    const ccomplex_t IT_0256 = IT_0093*IT_0217;
    const ccomplex_t IT_0257 = IT_0253 + IT_0254 + -IT_0255 + -IT_0256;
    const ccomplex_t IT_0258 = IT_0225*IT_0257;
    const ccomplex_t IT_0259 = IT_0089*IT_0211;
    const ccomplex_t IT_0260 = IT_0087*IT_0213;
    const ccomplex_t IT_0261 = IT_0093*IT_0215;
    const ccomplex_t IT_0262 = IT_0091*IT_0217;
    const ccomplex_t IT_0263 = IT_0259 + -IT_0260 + -IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = IT_0236*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0258 + (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = 6*IT_0265;
    const ccomplex_t IT_0267 = 0.5*IT_0258 + 0.5*IT_0264;
    const ccomplex_t IT_0268 = (-6)*IT_0267;
    const ccomplex_t IT_0269 = IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = 6*conj(IT_0265);
    const ccomplex_t IT_0271 = (-6)*conj(IT_0267);
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = -IT_0198;
    const ccomplex_t IT_0274 = IT_0197 + IT_0199 + IT_0200 + IT_0273;
    const ccomplex_t IT_0275 = IT_0165*IT_0203;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0162*IT_0206;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = -IT_0276 + -IT_0278;
    const ccomplex_t IT_0280 = 6*IT_0240;
    const ccomplex_t IT_0281 = (-6)*IT_0238;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 6*conj(IT_0240);
    const ccomplex_t IT_0284 = (-6)*conj(IT_0238);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0162*IT_0246;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = IT_0165*IT_0249;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0287 + -IT_0289;
    const ccomplex_t IT_0291 = 6*IT_0267;
    const ccomplex_t IT_0292 = (-6)*IT_0265;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = 6*conj(IT_0267);
    const ccomplex_t IT_0295 = (-6)*conj(IT_0265);
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_23 + s_24;
    const ccomplex_t IT_0298 = IT_0000*IT_0297;
    const ccomplex_t IT_0299 = conj(IT_0084) + conj(IT_0188);
    const ccomplex_t IT_0300 = 12*IT_0240;
    const ccomplex_t IT_0301 = (-12)*IT_0238;
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0084 + IT_0188;
    const ccomplex_t IT_0304 = 12*conj(IT_0240);
    const ccomplex_t IT_0305 = (-12)*conj(IT_0238);
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = conj(IT_0086) + conj(IT_0189);
    const ccomplex_t IT_0308 = 12*IT_0238;
    const ccomplex_t IT_0309 = (-12)*IT_0240;
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = IT_0086 + IT_0189;
    const ccomplex_t IT_0312 = 12*conj(IT_0238);
    const ccomplex_t IT_0313 = (-12)*conj(IT_0240);
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = conj(IT_0031) + conj(IT_0192);
    const ccomplex_t IT_0316 = 12*IT_0265;
    const ccomplex_t IT_0317 = (-12)*IT_0267;
    const ccomplex_t IT_0318 = IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = IT_0031 + IT_0192;
    const ccomplex_t IT_0320 = 12*conj(IT_0265);
    const ccomplex_t IT_0321 = (-12)*conj(IT_0267);
    const ccomplex_t IT_0322 = IT_0320 + IT_0321;
    const ccomplex_t IT_0323 = conj(IT_0037) + conj(IT_0193);
    const ccomplex_t IT_0324 = 12*IT_0267;
    const ccomplex_t IT_0325 = (-12)*IT_0265;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0037 + IT_0193;
    const ccomplex_t IT_0328 = 12*conj(IT_0267);
    const ccomplex_t IT_0329 = (-12)*conj(IT_0265);
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = conj(IT_0159)*IT_0242 + IT_0159*IT_0245 + conj
      (IT_0146)*IT_0269 + IT_0146*IT_0272 + conj(IT_0179)*IT_0282 + IT_0179
      *IT_0285 + conj(IT_0168)*IT_0293 + IT_0168*IT_0296 + IT_0299*IT_0302 +
       IT_0303*IT_0306 + IT_0307*IT_0310 + IT_0311*IT_0314 + IT_0315*IT_0318 +
       IT_0319*IT_0322 + IT_0323*IT_0326 + IT_0327*IT_0330;
    const ccomplex_t IT_0332 = -IT_0199;
    const ccomplex_t IT_0333 = IT_0198 + IT_0200 + IT_0332;
    const ccomplex_t IT_0334 = 6*IT_0209;
    const ccomplex_t IT_0335 = 6*conj(IT_0209);
    const ccomplex_t IT_0336 = s_13 + s_14;
    const ccomplex_t IT_0337 = IT_0000*IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0159;
    const ccomplex_t IT_0339 = 6*conj(IT_0159);
    const ccomplex_t IT_0340 = 6*IT_0179;
    const ccomplex_t IT_0341 = 6*conj(IT_0179);
    const ccomplex_t IT_0342 = conj(IT_0146)*IT_0242 + IT_0146*IT_0245 + conj
      (IT_0168)*IT_0282 + IT_0168*IT_0285 + IT_0310*IT_0315 + IT_0307*IT_0318 +
       IT_0314*IT_0319 + IT_0311*IT_0322 + IT_0302*IT_0323 + IT_0299*IT_0326 +
       IT_0306*IT_0327 + IT_0303*IT_0330 + (-6)*conj(IT_0265)*(IT_0179 + (
      -0.166666666666667)*IT_0338) + (-6)*IT_0265*(conj(IT_0179) + (
      -0.166666666666667)*IT_0339) + (-6)*conj(IT_0267)*(IT_0159 + (
      -0.166666666666667)*IT_0340) + (-6)*IT_0267*(conj(IT_0159) + (
      -0.166666666666667)*IT_0341);
    const ccomplex_t IT_0343 = s_34*IT_0026;
    const ccomplex_t IT_0344 = 12*IT_0179;
    const ccomplex_t IT_0345 = 24*IT_0084;
    const ccomplex_t IT_0346 = 24*IT_0188;
    const ccomplex_t IT_0347 = IT_0344 + IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = 12*IT_0168;
    const ccomplex_t IT_0349 = 24*IT_0193;
    const ccomplex_t IT_0350 = 24*IT_0037;
    const ccomplex_t IT_0351 = IT_0348 + IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 12*conj(IT_0168);
    const ccomplex_t IT_0353 = 24*conj(IT_0193);
    const ccomplex_t IT_0354 = 24*conj(IT_0037);
    const ccomplex_t IT_0355 = IT_0352 + IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = 12*conj(IT_0179);
    const ccomplex_t IT_0357 = 24*conj(IT_0084);
    const ccomplex_t IT_0358 = 24*conj(IT_0188);
    const ccomplex_t IT_0359 = IT_0356 + IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = 12*IT_0159;
    const ccomplex_t IT_0361 = 12*conj(IT_0159);
    const ccomplex_t IT_0362 = 12*IT_0146;
    const ccomplex_t IT_0363 = 12*conj(IT_0146);
    const ccomplex_t IT_0364 = 24*IT_0031;
    const ccomplex_t IT_0365 = 24*conj(IT_0031);
    const ccomplex_t IT_0366 = conj(IT_0168)*IT_0338 + IT_0168*IT_0339 + conj
      (IT_0146)*IT_0340 + IT_0146*IT_0341 + conj(IT_0031)*(IT_0344 + IT_0345) +
       conj(IT_0192)*IT_0347 + IT_0307*IT_0351 + IT_0311*IT_0355 + IT_0031*
      (IT_0356 + IT_0357) + IT_0192*IT_0359 + IT_0323*IT_0360 + IT_0327*IT_0361 
      + IT_0299*IT_0362 + IT_0303*IT_0363 + conj(IT_0188)*IT_0364 + IT_0188
      *IT_0365;
    const ccomplex_t IT_0367 = -s_34;
    const ccomplex_t IT_0368 = IT_0122 + IT_0367;
    const ccomplex_t IT_0369 = s_12*IT_0368;
    const ccomplex_t IT_0370 = (-24)*conj(IT_0106);
    const ccomplex_t IT_0371 = 24*conj(IT_0104);
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = (-24)*IT_0076;
    const ccomplex_t IT_0374 = s_34 + IT_0122;
    const ccomplex_t IT_0375 = s_12*IT_0374;
    const ccomplex_t IT_0376 = (-24)*IT_0240;
    const ccomplex_t IT_0377 = (-24)*conj(IT_0240);
    const ccomplex_t IT_0378 = 6*conj(IT_0074);
    const ccomplex_t IT_0379 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = 6*IT_0104;
    const ccomplex_t IT_0382 = (-6)*IT_0106;
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = (-6)*IT_0076;
    const ccomplex_t IT_0385 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0386 = m_c*s_14*m_N_4;
    const ccomplex_t IT_0387 = 12*IT_0209;
    const ccomplex_t IT_0388 = 24*IT_0279;
    const ccomplex_t IT_0389 = IT_0387 + IT_0388;
    const ccomplex_t IT_0390 = 12*conj(IT_0209);
    const ccomplex_t IT_0391 = 24*conj(IT_0279);
    const ccomplex_t IT_0392 = IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = 12*IT_0279;
    const ccomplex_t IT_0394 = IT_0334 + IT_0393;
    const ccomplex_t IT_0395 = 12*conj(IT_0279);
    const ccomplex_t IT_0396 = IT_0335 + IT_0395;
    const ccomplex_t IT_0397 = 12*IT_0252;
    const ccomplex_t IT_0398 = 12*conj(IT_0252);
    const ccomplex_t IT_0399 = conj(IT_0252)*IT_0340 + IT_0252*IT_0341 + conj
      (IT_0290)*IT_0347 + IT_0290*(IT_0356 + IT_0357 + IT_0358) + IT_0323
      *IT_0389 + IT_0327*IT_0392 + conj(IT_0168)*IT_0394 + IT_0168*IT_0396 +
       IT_0299*IT_0397 + IT_0303*IT_0398;
    const ccomplex_t IT_0400 = m_c*s_23*m_N_4;
    const ccomplex_t IT_0401 = 6*conj(IT_0252);
    const ccomplex_t IT_0402 = conj(IT_0168)*IT_0252 + 0.166666666666667*conj
      (IT_0290)*IT_0351 + 0.166666666666667*IT_0290*IT_0355 + 0.166666666666667
      *IT_0299*IT_0389 + 0.166666666666667*IT_0303*IT_0392 + 0.166666666666667
      *conj(IT_0179)*IT_0394 + 0.166666666666667*IT_0179*IT_0396 +
       0.166666666666667*IT_0323*IT_0397 + 0.166666666666667*IT_0327*IT_0398 +
       0.166666666666667*IT_0168*IT_0401;
    const ccomplex_t IT_0403 = IT_0026*IT_0368;
    const ccomplex_t IT_0404 = 6*conj(IT_0076);
    const ccomplex_t IT_0405 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = 6*IT_0106;
    const ccomplex_t IT_0408 = (-6)*IT_0104;
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = 6*conj(IT_0106);
    const ccomplex_t IT_0411 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = 6*IT_0076;
    const ccomplex_t IT_0414 = conj(IT_0076)*((-24)*IT_0104 + 24*IT_0106) +
       conj(IT_0074)*(24*IT_0104 + (-24)*IT_0106) + 24*IT_0076*conj(IT_0106) + 6
      *conj(IT_0104)*(IT_0209 + 0.166666666666667*IT_0373) + IT_0252*IT_0380 +
       conj(IT_0209)*IT_0383 + conj(IT_0252)*IT_0384 + IT_0209*IT_0385 + (-6)
      *IT_0074*(conj(IT_0290) + (-0.166666666666667)*IT_0372 + (
      -0.166666666666667)*IT_0401) + IT_0290*IT_0406 + conj(IT_0279)*IT_0409 +
       IT_0279*IT_0412 + conj(IT_0290)*IT_0413;
    const ccomplex_t IT_0415 = m_c*s_24*m_N_4;
    const ccomplex_t IT_0416 = 24*IT_0192;
    const ccomplex_t IT_0417 = IT_0362 + IT_0364 + IT_0416;
    const ccomplex_t IT_0418 = 24*IT_0209;
    const ccomplex_t IT_0419 = IT_0393 + IT_0418;
    const ccomplex_t IT_0420 = 24*conj(IT_0209);
    const ccomplex_t IT_0421 = IT_0395 + IT_0420;
    const ccomplex_t IT_0422 = 6*IT_0279;
    const ccomplex_t IT_0423 = IT_0387 + IT_0422;
    const ccomplex_t IT_0424 = 6*conj(IT_0279);
    const ccomplex_t IT_0425 = IT_0390 + IT_0424;
    const ccomplex_t IT_0426 = 6*IT_0146;
    const ccomplex_t IT_0427 = 6*conj(IT_0146);
    const ccomplex_t IT_0428 = 12*IT_0290;
    const ccomplex_t IT_0429 = 12*conj(IT_0290);
    const ccomplex_t IT_0430 = 24*conj(IT_0192);
    const ccomplex_t IT_0431 = conj(IT_0252)*IT_0417 + IT_0307*IT_0419 +
       IT_0311*IT_0421 + conj(IT_0159)*IT_0423 + IT_0159*IT_0425 + conj(IT_0290)
      *IT_0426 + IT_0290*IT_0427 + IT_0315*IT_0428 + IT_0319*IT_0429 + IT_0252*
      (24*conj(IT_0031) + IT_0363 + IT_0430);
    const ccomplex_t IT_0432 = conj(IT_0179)*IT_0338 + IT_0179*IT_0339 +
       IT_0307*IT_0347 + IT_0315*IT_0348 + IT_0319*IT_0352 + IT_0311*IT_0359 +
       IT_0299*IT_0360 + IT_0303*IT_0361 + IT_0323*(IT_0362 + IT_0364 + IT_0416)
       + conj(IT_0168)*IT_0426 + IT_0168*IT_0427 + IT_0327*(IT_0363 + IT_0365 +
       IT_0430);
    const ccomplex_t IT_0433 = m_c*s_13*m_N_4;
    const ccomplex_t IT_0434 = 24*conj(IT_0189);
    const ccomplex_t IT_0435 = 24*conj(IT_0086);
    const ccomplex_t IT_0436 = IT_0361 + IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = 24*IT_0086;
    const ccomplex_t IT_0438 = 24*IT_0189;
    const ccomplex_t IT_0439 = conj(IT_0290)*IT_0338 + IT_0290*IT_0339 +
       IT_0315*IT_0419 + IT_0319*IT_0421 + conj(IT_0146)*IT_0423 + IT_0146
      *IT_0425 + IT_0307*IT_0428 + IT_0311*IT_0429 + IT_0252*IT_0436 + conj
      (IT_0252)*(IT_0360 + IT_0437 + IT_0438);
    const ccomplex_t IT_0440 = 48*IT_0086;
    const ccomplex_t IT_0441 = 24*IT_0159;
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = 48*IT_0031;
    const ccomplex_t IT_0444 = 24*IT_0146;
    const ccomplex_t IT_0445 = IT_0443 + IT_0444;
    const ccomplex_t IT_0446 = 48*conj(IT_0031);
    const ccomplex_t IT_0447 = 48*IT_0192;
    const ccomplex_t IT_0448 = 48*conj(IT_0086);
    const ccomplex_t IT_0449 = 24*conj(IT_0159);
    const ccomplex_t IT_0450 = conj(IT_0209)*IT_0387 + conj(IT_0252)*IT_0397 +
       conj(IT_0146)*IT_0417 + conj(IT_0159)*(IT_0360 + IT_0437) + IT_0307
      *IT_0442 + conj(IT_0031)*IT_0445 + IT_0192*IT_0446 + conj(IT_0192)*
      (IT_0445 + IT_0447) + IT_0189*(48*conj(IT_0189) + IT_0448 + IT_0449);
    const ccomplex_t IT_0451 = 48*IT_0084;
    const ccomplex_t IT_0452 = 24*IT_0179;
    const ccomplex_t IT_0453 = IT_0451 + IT_0452;
    const ccomplex_t IT_0454 = 48*IT_0188;
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = 48*conj(IT_0084);
    const ccomplex_t IT_0457 = 48*conj(IT_0037);
    const ccomplex_t IT_0458 = 48*IT_0037;
    const ccomplex_t IT_0459 = 24*IT_0168;
    const ccomplex_t IT_0460 = conj(IT_0179)*IT_0347 + conj(IT_0168)*(IT_0348 
      + IT_0350) + conj(IT_0279)*IT_0393 + conj(IT_0290)*IT_0428 + conj(IT_0084)
      *IT_0453 + conj(IT_0188)*IT_0455 + IT_0188*IT_0456 + IT_0193*(24*conj
      (IT_0168) + 48*conj(IT_0193) + IT_0457) + IT_0323*(IT_0458 + IT_0459);
    const ccomplex_t IT_0461 = IT_0026*IT_0122;
    const ccomplex_t IT_0462 = conj(IT_0168)*IT_0347 + IT_0168*(IT_0356 +
       IT_0357) + conj(IT_0290)*IT_0394 + IT_0290*IT_0396 + conj(IT_0252)
      *IT_0423 + IT_0252*IT_0425 + IT_0146*IT_0436 + conj(IT_0146)*(IT_0360 +
       IT_0437 + IT_0438) + conj(IT_0192)*(48*IT_0189 + IT_0440 + IT_0441) +
       conj(IT_0031)*IT_0442 + IT_0189*IT_0446 + conj(IT_0189)*(IT_0443 +
       IT_0447) + IT_0319*(IT_0448 + IT_0449) + conj(IT_0037)*IT_0453 + conj
      (IT_0193)*IT_0455 + IT_0327*(24*conj(IT_0179) + IT_0456) + IT_0188*IT_0457
       + conj(IT_0188)*(48*IT_0193 + IT_0458 + IT_0459);
    const ccomplex_t IT_0463 = (-0.166666666666667)*IT_0252;
    const ccomplex_t IT_0464 = (-0.166666666666667)*conj(IT_0252);
    const ccomplex_t IT_0465 = IT_0026*IT_0374;
    const ccomplex_t IT_0466 = IT_0285*IT_0290 + IT_0282*conj(IT_0290) + 24
      *conj(IT_0265)*(IT_0238 + -1./4*IT_0279 + 0.0416666666666667*IT_0334 +
       0.0416666666666667*IT_0376) + 24*IT_0265*(conj(IT_0238) + -1./4*conj
      (IT_0279) + 0.0416666666666667*IT_0335 + 0.0416666666666667*IT_0377) + (-6
      )*conj(IT_0267)*(IT_0209 + 4*IT_0238 + (-4)*IT_0240 + (-0.166666666666667)
      *IT_0422) + (-6)*IT_0267*(conj(IT_0209) + 4*conj(IT_0238) + (-4)*conj
      (IT_0240) + (-0.166666666666667)*IT_0424) + (-6)*IT_0245*IT_0463 + (-6)
      *IT_0242*IT_0464;
    const ccomplex_t IT_0467 = IT_0004*((conj(IT_0031) + conj(IT_0037))
      *IT_0078 + (IT_0031 + IT_0037)*IT_0081 + (conj(IT_0084) + conj(IT_0086))
      *IT_0108 + (IT_0084 + IT_0086)*IT_0111) + (IT_0081*(IT_0084 + IT_0086) +
       IT_0078*(conj(IT_0084) + conj(IT_0086)) + (conj(IT_0031) + conj(IT_0037))
      *IT_0108 + (IT_0031 + IT_0037)*IT_0111)*IT_0115 + IT_0116*(IT_0081*IT_0146
       + IT_0078*conj(IT_0146) + IT_0111*IT_0159 + IT_0108*conj(IT_0159)) + 
      (IT_0111*IT_0146 + IT_0108*conj(IT_0146) + IT_0081*IT_0159 + IT_0078*conj
      (IT_0159))*IT_0160 + IT_0161*(conj(IT_0168)*IT_0171 + IT_0168*IT_0174 +
       conj(IT_0179)*IT_0182 + IT_0179*IT_0185) + (IT_0174*IT_0179 + IT_0171
      *conj(IT_0179) + conj(IT_0168)*IT_0182 + IT_0168*IT_0185)*IT_0186 +
       IT_0187*(IT_0171*IT_0190 + IT_0174*IT_0191 + IT_0182*IT_0194 + IT_0185
      *IT_0195) + (IT_0182*IT_0190 + IT_0185*IT_0191 + IT_0171*IT_0194 + IT_0174
      *IT_0195)*IT_0196 + IT_0202*(conj(IT_0209)*IT_0242 + IT_0209*IT_0245 +
       conj(IT_0252)*IT_0269 + IT_0252*IT_0272) + IT_0274*(conj(IT_0279)*IT_0282
       + IT_0279*IT_0285 + conj(IT_0290)*IT_0293 + IT_0290*IT_0296) + IT_0298
      *IT_0331 + 6*IT_0333*(conj(IT_0252)*IT_0290 + IT_0252*conj(IT_0290) +
       0.166666666666667*conj(IT_0279)*IT_0334 + 0.166666666666667*IT_0279
      *IT_0335) + IT_0337*IT_0342 + IT_0343*IT_0366 + IT_0369*(24*conj(IT_0074)
      *IT_0076 + IT_0074*((-24)*conj(IT_0074) + 24*conj(IT_0076)) + IT_0104*((
      -24)*conj(IT_0104) + 24*conj(IT_0106)) + IT_0106*IT_0372 + conj(IT_0076)
      *IT_0373) + IT_0375*(24*IT_0240*conj(IT_0240) + IT_0267*((-24)*conj
      (IT_0265) + 24*conj(IT_0267)) + IT_0265*(24*conj(IT_0265) + (-24)*conj
      (IT_0267)) + conj(IT_0238)*IT_0376 + IT_0238*(24*conj(IT_0238) + IT_0377))
       + (IT_0197 + IT_0198 + IT_0201 + IT_0332)*(IT_0279*IT_0380 + conj(IT_0290
      )*IT_0383 + 6*conj(IT_0279)*(IT_0074 + 0.166666666666667*IT_0384) + 6
      *IT_0290*(conj(IT_0104) + 0.166666666666667*IT_0385)) + IT_0386*IT_0399 +
       6*IT_0400*IT_0402 + IT_0403*IT_0414 + IT_0415*IT_0431 + IT_0197*IT_0432 +
       IT_0433*IT_0439 + IT_0198*IT_0450 + IT_0200*IT_0460 + IT_0461*IT_0462 + 
      (IT_0197 + IT_0200 + IT_0273 + IT_0332)*(IT_0209*IT_0406 + (-6)*conj
      (IT_0209)*(IT_0074 + (-0.166666666666667)*IT_0413) + (-6)*IT_0412*IT_0463 
      + (-6)*IT_0409*IT_0464) + IT_0465*IT_0466;
    return create_ccomplex_return(IT_0467);
}

