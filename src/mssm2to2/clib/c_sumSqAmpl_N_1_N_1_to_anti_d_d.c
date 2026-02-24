#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_d_d(
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
    const creal_t m_N_1 = param->m_N_1;
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
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0004, -1);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 3*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0015 = IT_0005*IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0014;
    const ccomplex_t IT_0017 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0018 = IT_0005*IT_0017;
    const ccomplex_t IT_0019 = IT_0008*IT_0017;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + IT_0016 + 
      -IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0003;
    const ccomplex_t IT_0025 = (-0.333333333333333)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0003;
    const ccomplex_t IT_0027 = (-0.333333333333333)*IT_0026;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = pow(m_d, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0000 + IT_0029 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0007*IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0007*IT_0033*IT_0035;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0000 + IT_0029 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0032 + -IT_0043;
    const ccomplex_t IT_0045 = s_12*IT_0029;
    const ccomplex_t IT_0046 = 12*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = sin(beta);
    const ccomplex_t IT_0050 = m_d*e_em*IT_0007*IT_0033*IT_0035*IT_0049;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0048*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0003*IT_0034;
    const ccomplex_t IT_0055 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = IT_0003*IT_0049;
    const ccomplex_t IT_0058 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = IT_0007*IT_0034;
    const ccomplex_t IT_0061 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = IT_0007*IT_0049;
    const ccomplex_t IT_0064 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = IT_0056 + -IT_0059 + -IT_0062 + IT_0065;
    const ccomplex_t IT_0067 = IT_0053*IT_0066;
    const ccomplex_t IT_0068 = pow(m_Z, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0068 + 
      -reg_prop, -1);
    const ccomplex_t IT_0070 = m_d*e_em*IT_0007*IT_0033;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0054*IT_0058;
    const ccomplex_t IT_0075 = IT_0055*IT_0057;
    const ccomplex_t IT_0076 = IT_0060*IT_0064;
    const ccomplex_t IT_0077 = IT_0061*IT_0063;
    const ccomplex_t IT_0078 = IT_0074 + IT_0075 + -IT_0076 + -IT_0077;
    const ccomplex_t IT_0079 = IT_0073*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0067 + 0.5*IT_0079;
    const ccomplex_t IT_0081 = m_d*m_N_1;
    const ccomplex_t IT_0082 = s_13 + s_14;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 12*IT_0083;
    const ccomplex_t IT_0085 = IT_0080*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B1)*e_em;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em;
    const ccomplex_t IT_0090 = IT_0007*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-3)*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0037*IT_0093;
    const ccomplex_t IT_0095 = IT_0041*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0036;
    const ccomplex_t IT_0098 = IT_0025*IT_0097;
    const ccomplex_t IT_0099 = IT_0030*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0096 + IT_0100;
    const ccomplex_t IT_0102 = m_d*s_14*m_N_1;
    const ccomplex_t IT_0103 = 12*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = N_B1*e_em;
    const ccomplex_t IT_0106 = IT_0003*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = N_W1*e_em;
    const ccomplex_t IT_0109 = IT_0007*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + (-3)*IT_0110);
    const ccomplex_t IT_0112 = 0.166666666666667*IT_0111;
    const ccomplex_t IT_0113 = IT_0039*IT_0112;
    const ccomplex_t IT_0114 = IT_0041*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0038;
    const ccomplex_t IT_0117 = IT_0027*IT_0116;
    const ccomplex_t IT_0118 = IT_0030*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0115 + IT_0119;
    const ccomplex_t IT_0121 = m_d*s_23*m_N_1;
    const ccomplex_t IT_0122 = 12*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = IT_0093*IT_0112;
    const ccomplex_t IT_0125 = IT_0041*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0097*IT_0116;
    const ccomplex_t IT_0128 = IT_0030*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = -IT_0126 + -IT_0129;
    const ccomplex_t IT_0131 = s_34*IT_0000;
    const ccomplex_t IT_0132 = 12*IT_0131;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0135 = IT_0054*IT_0134;
    const ccomplex_t IT_0136 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0137 = IT_0057*IT_0136;
    const ccomplex_t IT_0138 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0139 = IT_0060*IT_0138;
    const ccomplex_t IT_0140 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0141 = IT_0063*IT_0140;
    const ccomplex_t IT_0142 = -IT_0135 + IT_0137 + IT_0139 + -IT_0141;
    const ccomplex_t IT_0143 = IT_0053*IT_0142;
    const ccomplex_t IT_0144 = IT_0054*IT_0136;
    const ccomplex_t IT_0145 = IT_0057*IT_0134;
    const ccomplex_t IT_0146 = IT_0060*IT_0140;
    const ccomplex_t IT_0147 = IT_0063*IT_0138;
    const ccomplex_t IT_0148 = -IT_0144 + -IT_0145 + IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = IT_0073*IT_0148;
    const ccomplex_t IT_0150 = (-0.5)*IT_0143 + (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = s_23 + s_24;
    const ccomplex_t IT_0152 = IT_0081*IT_0151;
    const ccomplex_t IT_0153 = 12*IT_0152;
    const ccomplex_t IT_0154 = IT_0150*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0156 = 0.333333333333333*IT_0155;
    const ccomplex_t IT_0157 = IT_0001*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = IT_0021*IT_0158;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0000*IT_0029;
    const ccomplex_t IT_0162 = 48*IT_0161;
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0067 + (-0.5)*IT_0079;
    const ccomplex_t IT_0165 = (-12)*IT_0083;
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = 0.5*IT_0159;
    const ccomplex_t IT_0168 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0169 = IT_0124*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0172 = IT_0127*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0170 + IT_0173;
    const ccomplex_t IT_0175 = IT_0046*IT_0174;
    const ccomplex_t IT_0176 = IT_0084*IT_0164;
    const ccomplex_t IT_0177 = IT_0094*IT_0168;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0098*IT_0171;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = -IT_0178 + -IT_0180;
    const ccomplex_t IT_0182 = m_d*s_13*m_N_1;
    const ccomplex_t IT_0183 = 12*IT_0182;
    const ccomplex_t IT_0184 = IT_0181*IT_0183;
    const ccomplex_t IT_0185 = 0.5*IT_0020;
    const ccomplex_t IT_0186 = IT_0158*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = s_13*s_24;
    const ccomplex_t IT_0189 = 48*IT_0188;
    const ccomplex_t IT_0190 = IT_0187*IT_0189;
    const ccomplex_t IT_0191 = IT_0080*IT_0165;
    const ccomplex_t IT_0192 = (-12)*IT_0152;
    const ccomplex_t IT_0193 = IT_0150*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0022;
    const ccomplex_t IT_0195 = m_d*s_24*m_N_1;
    const ccomplex_t IT_0196 = 12*IT_0195;
    const ccomplex_t IT_0197 = IT_0181*IT_0196;
    const ccomplex_t IT_0198 = IT_0132*IT_0174;
    const ccomplex_t IT_0199 = IT_0153*IT_0164;
    const ccomplex_t IT_0200 = IT_0150*IT_0165;
    const ccomplex_t IT_0201 = IT_0080*IT_0192;
    const ccomplex_t IT_0202 = IT_0046*IT_0130;
    const ccomplex_t IT_0203 = IT_0084*IT_0150;
    const ccomplex_t IT_0204 = IT_0103*IT_0120;
    const ccomplex_t IT_0205 = IT_0101*IT_0122;
    const ccomplex_t IT_0206 = IT_0044*IT_0132;
    const ccomplex_t IT_0207 = IT_0080*IT_0153;
    const ccomplex_t IT_0208 = IT_0164*IT_0192;
    const ccomplex_t IT_0209 = 0.5*IT_0186;
    const ccomplex_t IT_0210 = 24*conj(IT_0209);
    const ccomplex_t IT_0211 = 24*conj(IT_0160);
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = 24*conj(IT_0023);
    const ccomplex_t IT_0214 = IT_0013*IT_0185;
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = 24*conj(IT_0215);
    const ccomplex_t IT_0217 = IT_0213 + IT_0216;
    const ccomplex_t IT_0218 = s_12*s_34;
    const ccomplex_t IT_0219 = s_14*s_23;
    const ccomplex_t IT_0220 = -IT_0188;
    const ccomplex_t IT_0221 = IT_0045 + IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0218 + IT_0221;
    const ccomplex_t IT_0223 = 6*IT_0222;
    const ccomplex_t IT_0224 = -IT_0218;
    const ccomplex_t IT_0225 = IT_0188 + IT_0219 + IT_0224;
    const ccomplex_t IT_0226 = 6*IT_0225;
    const ccomplex_t IT_0227 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0228 = sin(alpha);
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*m_d*e_em*IT_0007*IT_0033
      *IT_0035*IT_0228;
    const ccomplex_t IT_0230 = 0.5*IT_0229;
    const ccomplex_t IT_0231 = IT_0227*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = cos(alpha);
    const ccomplex_t IT_0234 = IT_0003*IT_0233;
    const ccomplex_t IT_0235 = IT_0055*IT_0234;
    const ccomplex_t IT_0236 = IT_0003*IT_0228;
    const ccomplex_t IT_0237 = IT_0058*IT_0236;
    const ccomplex_t IT_0238 = IT_0007*IT_0233;
    const ccomplex_t IT_0239 = IT_0061*IT_0238;
    const ccomplex_t IT_0240 = IT_0007*IT_0228;
    const ccomplex_t IT_0241 = IT_0064*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*(IT_0235 + IT_0237 + 
      -IT_0239 + -IT_0241);
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = IT_0232*IT_0243;
    const ccomplex_t IT_0245 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*m_d*e_em*IT_0007*IT_0033
      *IT_0035*IT_0233;
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0245*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0058*IT_0234;
    const ccomplex_t IT_0251 = IT_0055*IT_0236;
    const ccomplex_t IT_0252 = IT_0064*IT_0238;
    const ccomplex_t IT_0253 = IT_0061*IT_0240;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*(IT_0250 + -IT_0251 + 
      -IT_0252 + IT_0253);
    const ccomplex_t IT_0255 = IT_0249*IT_0254;
    const ccomplex_t IT_0256 = (-0.5)*IT_0244 + (-0.5)*IT_0255;
    const ccomplex_t IT_0257 = -IT_0219;
    const ccomplex_t IT_0258 = IT_0045 + IT_0188 + IT_0257;
    const ccomplex_t IT_0259 = IT_0224 + IT_0258;
    const ccomplex_t IT_0260 = 6*IT_0259;
    const ccomplex_t IT_0261 = 6*IT_0182;
    const ccomplex_t IT_0262 = 6*IT_0195;
    const ccomplex_t IT_0263 = 6*IT_0161;
    const ccomplex_t IT_0264 = IT_0136*IT_0234;
    const ccomplex_t IT_0265 = IT_0134*IT_0236;
    const ccomplex_t IT_0266 = IT_0140*IT_0238;
    const ccomplex_t IT_0267 = IT_0138*IT_0240;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0264 + -IT_0265 + 
      -IT_0266 + IT_0267);
    const ccomplex_t IT_0269 = IT_0249*IT_0268;
    const ccomplex_t IT_0270 = IT_0134*IT_0234;
    const ccomplex_t IT_0271 = IT_0136*IT_0236;
    const ccomplex_t IT_0272 = IT_0138*IT_0238;
    const ccomplex_t IT_0273 = IT_0140*IT_0240;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*(IT_0270 + IT_0271 + 
      -IT_0272 + -IT_0273);
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = IT_0232*IT_0275;
    const ccomplex_t IT_0277 = 0.5*IT_0269 + 0.5*IT_0276;
    const ccomplex_t IT_0278 = -s_34;
    const ccomplex_t IT_0279 = IT_0029 + IT_0278;
    const ccomplex_t IT_0280 = IT_0000*IT_0279;
    const ccomplex_t IT_0281 = 6*IT_0280;
    const ccomplex_t IT_0282 = s_34 + IT_0029;
    const ccomplex_t IT_0283 = IT_0000*IT_0282;
    const ccomplex_t IT_0284 = 6*IT_0283;
    const ccomplex_t IT_0285 = 12*IT_0219;
    const ccomplex_t IT_0286 = IT_0028*IT_0171;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = IT_0040*IT_0168;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = IT_0117*IT_0171;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = IT_0113*IT_0168;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = -IT_0292 + -IT_0294;
    const ccomplex_t IT_0296 = 12*IT_0161;
    const ccomplex_t IT_0297 = (-6)*IT_0222;
    const ccomplex_t IT_0298 = 0.5*IT_0244 + 0.5*IT_0255;
    const ccomplex_t IT_0299 = (-6)*IT_0259;
    const ccomplex_t IT_0300 = (-0.5)*IT_0269 + (-0.5)*IT_0276;
    const ccomplex_t IT_0301 = (-6)*IT_0280;
    const ccomplex_t IT_0302 = 0.5*IT_0143 + 0.5*IT_0149;
    const ccomplex_t IT_0303 = (-6)*IT_0283;
    const ccomplex_t IT_0304 = IT_0103*conj(IT_0174) + IT_0121*IT_0212 +
       IT_0102*IT_0217 + conj(IT_0080)*IT_0223 + conj(IT_0101)*IT_0226 + conj
      (IT_0256)*IT_0260 + conj(IT_0044)*IT_0261 + conj(IT_0130)*IT_0262 + conj
      (IT_0120)*IT_0263 + conj(IT_0277)*IT_0281 + conj(IT_0150)*IT_0284 + conj
      (IT_0181)*IT_0285 + IT_0122*conj(IT_0290) + conj(IT_0295)*IT_0296 + conj
      (IT_0164)*IT_0297 + conj(IT_0298)*IT_0299 + conj(IT_0300)*IT_0301 + conj
      (IT_0302)*IT_0303;
    const ccomplex_t IT_0305 = 0.5*IT_0214;
    const ccomplex_t IT_0306 = 24*conj(IT_0305);
    const ccomplex_t IT_0307 = 24*conj(IT_0194);
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = 24*conj(IT_0167);
    const ccomplex_t IT_0310 = 24*conj(IT_0187);
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = 6*IT_0218 + 6*IT_0258;
    const ccomplex_t IT_0313 = IT_0221 + IT_0224;
    const ccomplex_t IT_0314 = 6*IT_0313;
    const ccomplex_t IT_0315 = 6*IT_0102;
    const ccomplex_t IT_0316 = 6*IT_0121;
    const ccomplex_t IT_0317 = 12*IT_0188;
    const ccomplex_t IT_0318 = (-6)*IT_0218 + (-6)*IT_0258;
    const ccomplex_t IT_0319 = (-6)*IT_0313;
    const ccomplex_t IT_0320 = conj(IT_0044)*IT_0183 + conj(IT_0130)*IT_0196 +
       conj(IT_0181)*IT_0226 + IT_0263*conj(IT_0295) + conj(IT_0120)*IT_0296 +
       IT_0281*conj(IT_0300) + conj(IT_0277)*IT_0301 + IT_0284*conj(IT_0302) +
       conj(IT_0150)*IT_0303 + IT_0195*IT_0308 + IT_0182*IT_0311 + conj(IT_0164)
      *IT_0312 + conj(IT_0298)*IT_0314 + conj(IT_0174)*IT_0315 + conj(IT_0290)
      *IT_0316 + conj(IT_0101)*IT_0317 + conj(IT_0080)*IT_0318 + conj(IT_0256)
      *IT_0319;
    const ccomplex_t IT_0321 = conj(IT_0130)*IT_0183 + conj(IT_0044)*IT_0196 +
       conj(IT_0181)*IT_0263 + conj(IT_0256)*IT_0281 + conj(IT_0164)*IT_0284 +
       IT_0226*conj(IT_0295) + conj(IT_0101)*IT_0296 + conj(IT_0298)*IT_0301 +
       conj(IT_0080)*IT_0303 + IT_0182*IT_0308 + IT_0195*IT_0311 + conj(IT_0302)
      *IT_0312 + conj(IT_0277)*IT_0314 + conj(IT_0290)*IT_0315 + conj(IT_0174)
      *IT_0316 + conj(IT_0120)*IT_0317 + conj(IT_0150)*IT_0318 + conj(IT_0300)
      *IT_0319;
    const ccomplex_t IT_0322 = conj(IT_0194) + conj(IT_0305);
    const ccomplex_t IT_0323 = conj(IT_0167) + conj(IT_0187);
    const ccomplex_t IT_0324 = IT_0122*conj(IT_0174) + IT_0102*IT_0212 +
       IT_0121*IT_0217 + conj(IT_0150)*IT_0223 + conj(IT_0120)*IT_0226 + conj
      (IT_0130)*IT_0261 + conj(IT_0044)*IT_0262 + conj(IT_0101)*IT_0263 + conj
      (IT_0080)*IT_0284 + IT_0103*conj(IT_0290) + IT_0285*conj(IT_0295) + conj
      (IT_0181)*IT_0296 + IT_0281*conj(IT_0298) + conj(IT_0277)*IT_0299 +
       IT_0260*conj(IT_0300) + conj(IT_0256)*IT_0301 + IT_0297*conj(IT_0302) +
       conj(IT_0164)*IT_0303 + IT_0183*IT_0322 + IT_0196*IT_0323;
    const ccomplex_t IT_0325 = -s_14;
    const ccomplex_t IT_0326 = s_13 + IT_0325;
    const ccomplex_t IT_0327 = IT_0081*IT_0326;
    const ccomplex_t IT_0328 = (-2)*IT_0327;
    const ccomplex_t IT_0329 = 3*IT_0328;
    const ccomplex_t IT_0330 = -s_24;
    const ccomplex_t IT_0331 = s_23 + IT_0330;
    const ccomplex_t IT_0332 = IT_0081*IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0332;
    const ccomplex_t IT_0334 = 3*IT_0333;
    const ccomplex_t IT_0335 = 6*IT_0045;
    const ccomplex_t IT_0336 = 6*IT_0083;
    const ccomplex_t IT_0337 = 6*IT_0131;
    const ccomplex_t IT_0338 = 6*IT_0152;
    const ccomplex_t IT_0339 = (-6)*IT_0083;
    const ccomplex_t IT_0340 = (-6)*IT_0152;
    const ccomplex_t IT_0341 = (-3)*IT_0328;
    const ccomplex_t IT_0342 = (-3)*IT_0333;
    const ccomplex_t IT_0343 = IT_0103*conj(IT_0181) + IT_0161*IT_0212 +
       IT_0217*IT_0219 + conj(IT_0174)*IT_0285 + IT_0122*conj(IT_0295) + conj
      (IT_0290)*IT_0296 + conj(IT_0101)*IT_0315 + conj(IT_0120)*IT_0316 + conj
      (IT_0298)*IT_0329 + conj(IT_0277)*IT_0334 + conj(IT_0044)*IT_0335 + conj
      (IT_0080)*IT_0336 + conj(IT_0130)*IT_0337 + conj(IT_0150)*IT_0338 + conj
      (IT_0164)*IT_0339 + conj(IT_0302)*IT_0340 + conj(IT_0256)*IT_0341 + conj
      (IT_0300)*IT_0342;
    const ccomplex_t IT_0344 = IT_0122*conj(IT_0181) + IT_0161*IT_0217 +
       IT_0212*IT_0219 + IT_0285*conj(IT_0290) + IT_0103*conj(IT_0295) + conj
      (IT_0174)*IT_0296 + conj(IT_0120)*IT_0315 + conj(IT_0101)*IT_0316 +
       IT_0046*IT_0322 + IT_0132*IT_0323 + conj(IT_0277)*IT_0329 + conj(IT_0298)
      *IT_0334 + conj(IT_0130)*IT_0335 + conj(IT_0150)*IT_0336 + conj(IT_0044)
      *IT_0337 + conj(IT_0080)*IT_0338 + conj(IT_0302)*IT_0339 + conj(IT_0164)
      *IT_0340 + conj(IT_0300)*IT_0341 + conj(IT_0256)*IT_0342;
    const ccomplex_t IT_0345 = 2*IT_0327;
    const ccomplex_t IT_0346 = 3*IT_0345;
    const ccomplex_t IT_0347 = (-2)*IT_0332;
    const ccomplex_t IT_0348 = 3*IT_0347;
    const ccomplex_t IT_0349 = (-3)*IT_0345;
    const ccomplex_t IT_0350 = (-3)*IT_0347;
    const ccomplex_t IT_0351 = conj(IT_0101)*IT_0183 + conj(IT_0120)*IT_0196 +
       conj(IT_0181)*IT_0261 + IT_0262*conj(IT_0295) + conj(IT_0130)*IT_0296 +
       IT_0161*IT_0308 + IT_0188*IT_0311 + conj(IT_0044)*IT_0317 + conj(IT_0174)
      *IT_0335 + conj(IT_0164)*IT_0336 + conj(IT_0290)*IT_0337 + conj(IT_0302)
      *IT_0338 + conj(IT_0080)*IT_0339 + conj(IT_0150)*IT_0340 + conj(IT_0256)
      *IT_0346 + conj(IT_0300)*IT_0348 + conj(IT_0298)*IT_0349 + conj(IT_0277)
      *IT_0350;
    const ccomplex_t IT_0352 = conj(IT_0120)*IT_0183 + conj(IT_0101)*IT_0196 +
       conj(IT_0181)*IT_0262 + IT_0261*conj(IT_0295) + conj(IT_0044)*IT_0296 +
       IT_0188*IT_0308 + IT_0161*IT_0311 + conj(IT_0130)*IT_0317 + conj(IT_0290)
      *IT_0335 + conj(IT_0302)*IT_0336 + conj(IT_0174)*IT_0337 + conj(IT_0164)
      *IT_0338 + conj(IT_0150)*IT_0339 + conj(IT_0080)*IT_0340 + conj(IT_0300)
      *IT_0346 + conj(IT_0256)*IT_0348 + conj(IT_0277)*IT_0349 + conj(IT_0298)
      *IT_0350;
    const ccomplex_t IT_0353 = (-4)*IT_0327;
    const ccomplex_t IT_0354 = 3*IT_0353;
    const ccomplex_t IT_0355 = IT_0298*IT_0354;
    const ccomplex_t IT_0356 = (-3)*IT_0353;
    const ccomplex_t IT_0357 = IT_0256*IT_0356;
    const ccomplex_t IT_0358 = (-4)*IT_0332;
    const ccomplex_t IT_0359 = 3*IT_0358;
    const ccomplex_t IT_0360 = IT_0256*IT_0359;
    const ccomplex_t IT_0361 = (-3)*IT_0358;
    const ccomplex_t IT_0362 = IT_0298*IT_0361;
    const ccomplex_t IT_0363 = 24*IT_0283;
    const ccomplex_t IT_0364 = s_12*IT_0282;
    const ccomplex_t IT_0365 = (-24)*IT_0364;
    const ccomplex_t IT_0366 = (-24)*IT_0283;
    const ccomplex_t IT_0367 = conj(IT_0160) + conj(IT_0209);
    const ccomplex_t IT_0368 = conj(IT_0023) + conj(IT_0215);
    const ccomplex_t IT_0369 = 24*IT_0364;
    const ccomplex_t IT_0370 = conj(IT_0101)*IT_0284 + conj(IT_0295)*IT_0297 +
       conj(IT_0181)*IT_0303 + conj(IT_0120)*IT_0312 + IT_0084*IT_0322 + IT_0153
      *IT_0323 + conj(IT_0130)*IT_0336 + conj(IT_0044)*IT_0338 + conj(IT_0290)
      *IT_0339 + conj(IT_0174)*IT_0340 + conj(IT_0164)*IT_0363 + conj(IT_0150)
      *IT_0365 + conj(IT_0080)*IT_0366 + IT_0165*IT_0367 + IT_0192*IT_0368 +
       conj(IT_0302)*IT_0369;
    const ccomplex_t IT_0371 = conj(IT_0181)*IT_0284 + IT_0223*conj(IT_0295) +
       conj(IT_0101)*IT_0303 + conj(IT_0120)*IT_0318 + conj(IT_0290)*IT_0336 +
       conj(IT_0174)*IT_0338 + conj(IT_0130)*IT_0339 + conj(IT_0044)*IT_0340 +
       conj(IT_0080)*IT_0363 + conj(IT_0302)*IT_0365 + conj(IT_0164)*IT_0366 +
       conj(IT_0150)*IT_0369;
    const ccomplex_t IT_0372 = conj(IT_0181)*IT_0223 + IT_0284*conj(IT_0295) +
       conj(IT_0120)*IT_0303 + conj(IT_0101)*IT_0318 + conj(IT_0174)*IT_0336 +
       conj(IT_0290)*IT_0338 + conj(IT_0044)*IT_0339 + conj(IT_0130)*IT_0340 +
       conj(IT_0150)*IT_0363 + conj(IT_0164)*IT_0365 + conj(IT_0302)*IT_0366 +
       conj(IT_0080)*IT_0369;
    const ccomplex_t IT_0373 = conj(IT_0120)*IT_0284 + conj(IT_0181)*IT_0297 +
       conj(IT_0295)*IT_0303 + conj(IT_0101)*IT_0312 + conj(IT_0044)*IT_0336 +
       conj(IT_0130)*IT_0338 + conj(IT_0174)*IT_0339 + conj(IT_0290)*IT_0340 +
       conj(IT_0302)*IT_0363 + conj(IT_0080)*IT_0365 + conj(IT_0150)*IT_0366 +
       conj(IT_0164)*IT_0369;
    const ccomplex_t IT_0374 = conj(IT_0160) + conj(IT_0194);
    const ccomplex_t IT_0375 = 4*IT_0332;
    const ccomplex_t IT_0376 = 3*IT_0375;
    const ccomplex_t IT_0377 = s_12*IT_0279;
    const ccomplex_t IT_0378 = (-24)*IT_0377;
    const ccomplex_t IT_0379 = (-24)*IT_0280;
    const ccomplex_t IT_0380 = conj(IT_0023) + conj(IT_0167);
    const ccomplex_t IT_0381 = 4*IT_0327;
    const ccomplex_t IT_0382 = (-3)*IT_0381;
    const ccomplex_t IT_0383 = 24*IT_0377;
    const ccomplex_t IT_0384 = 24*IT_0280;
    const ccomplex_t IT_0385 = IT_0281*conj(IT_0295) + conj(IT_0181)*IT_0299 +
       conj(IT_0120)*IT_0301 + conj(IT_0101)*IT_0314 + conj(IT_0174)*IT_0329 +
       conj(IT_0290)*IT_0334 + conj(IT_0044)*IT_0349 + conj(IT_0130)*IT_0350 +
       IT_0374*IT_0376 + conj(IT_0298)*IT_0378 + conj(IT_0300)*IT_0379 + IT_0380
      *IT_0382 + conj(IT_0256)*IT_0383 + conj(IT_0277)*IT_0384;
    const ccomplex_t IT_0386 = 3*IT_0381;
    const ccomplex_t IT_0387 = (-3)*IT_0375;
    const ccomplex_t IT_0388 = conj(IT_0181)*IT_0260 + conj(IT_0120)*IT_0281 +
       conj(IT_0295)*IT_0301 + conj(IT_0101)*IT_0319 + conj(IT_0174)*IT_0341 +
       conj(IT_0290)*IT_0342 + conj(IT_0044)*IT_0346 + conj(IT_0130)*IT_0348 +
       conj(IT_0256)*IT_0378 + conj(IT_0277)*IT_0379 + conj(IT_0298)*IT_0383 +
       conj(IT_0300)*IT_0384 + IT_0380*IT_0386 + IT_0374*IT_0387;
    const ccomplex_t IT_0389 = conj(IT_0209) + conj(IT_0305);
    const ccomplex_t IT_0390 = conj(IT_0181)*IT_0281 + conj(IT_0295)*IT_0299 +
       conj(IT_0101)*IT_0301 + conj(IT_0120)*IT_0314 + conj(IT_0290)*IT_0329 +
       conj(IT_0174)*IT_0334 + conj(IT_0130)*IT_0349 + conj(IT_0044)*IT_0350 +
       IT_0354*IT_0374 + (conj(IT_0187) + conj(IT_0215))*IT_0376 + conj(IT_0277)
      *IT_0378 + conj(IT_0256)*IT_0379 + IT_0361*IT_0380 + conj(IT_0300)*IT_0383
       + conj(IT_0298)*IT_0384 + IT_0382*IT_0389;
    const ccomplex_t IT_0391 = conj(IT_0101)*IT_0281 + IT_0260*conj(IT_0295) +
       conj(IT_0181)*IT_0301 + conj(IT_0120)*IT_0319 + conj(IT_0290)*IT_0341 +
       conj(IT_0174)*IT_0342 + conj(IT_0130)*IT_0346 + conj(IT_0044)*IT_0348 +
       IT_0356*IT_0374 + conj(IT_0300)*IT_0378 + conj(IT_0298)*IT_0379 + IT_0359
      *IT_0380 + conj(IT_0277)*IT_0383 + conj(IT_0256)*IT_0384 + (conj(IT_0187) 
      + conj(IT_0215))*IT_0387 + IT_0386*IT_0389;
    const ccomplex_t IT_0392 = conj(IT_0277)*IT_0354;
    const ccomplex_t IT_0393 = conj(IT_0298)*IT_0376;
    const ccomplex_t IT_0394 = IT_0046*conj(IT_0130);
    const ccomplex_t IT_0395 = IT_0084*conj(IT_0150);
    const ccomplex_t IT_0396 = IT_0103*conj(IT_0120);
    const ccomplex_t IT_0397 = conj(IT_0101)*IT_0122;
    const ccomplex_t IT_0398 = conj(IT_0044)*IT_0132;
    const ccomplex_t IT_0399 = conj(IT_0080)*IT_0153;
    const ccomplex_t IT_0400 = 48*IT_0219;
    const ccomplex_t IT_0401 = IT_0165*conj(IT_0302);
    const ccomplex_t IT_0402 = conj(IT_0164)*IT_0192;
    const ccomplex_t IT_0403 = conj(IT_0300)*IT_0356;
    const ccomplex_t IT_0404 = conj(IT_0256)*IT_0387;
    const ccomplex_t IT_0405 = 24*IT_0219;
    const ccomplex_t IT_0406 = conj(IT_0290)*IT_0405;
    const ccomplex_t IT_0407 = 24*IT_0102;
    const ccomplex_t IT_0408 = 24*IT_0121;
    const ccomplex_t IT_0409 = conj(IT_0181)*IT_0408;
    const ccomplex_t IT_0410 = 24*IT_0161;
    const ccomplex_t IT_0411 = conj(IT_0174)*IT_0410;
    const ccomplex_t IT_0412 = IT_0045*IT_0308;
    const ccomplex_t IT_0413 = IT_0131*IT_0311;
    const ccomplex_t IT_0414 = IT_0392 + IT_0393 + IT_0394 + IT_0395 + IT_0396
       + IT_0397 + IT_0398 + IT_0399 + IT_0367*IT_0400 + IT_0401 + IT_0402 +
       IT_0403 + IT_0404 + IT_0406 + conj(IT_0295)*IT_0407 + IT_0409 + IT_0411 +
       IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = conj(IT_0256)*IT_0359;
    const ccomplex_t IT_0416 = conj(IT_0300)*IT_0386;
    const ccomplex_t IT_0417 = conj(IT_0298)*IT_0361;
    const ccomplex_t IT_0418 = conj(IT_0277)*IT_0382;
    const ccomplex_t IT_0419 = IT_0162*IT_0368 + IT_0394 + IT_0395 + IT_0396 +
       IT_0397 + IT_0398 + IT_0399 + IT_0367*IT_0400 + IT_0401 + IT_0402 +
       IT_0406 + conj(IT_0295)*IT_0407 + IT_0409 + IT_0411 + IT_0412 + IT_0413 +
       IT_0415 + IT_0416 + IT_0417 + IT_0418;
    const ccomplex_t IT_0420 = conj(IT_0298)*IT_0354;
    const ccomplex_t IT_0421 = conj(IT_0277)*IT_0376;
    const ccomplex_t IT_0422 = conj(IT_0044)*IT_0046;
    const ccomplex_t IT_0423 = conj(IT_0080)*IT_0084;
    const ccomplex_t IT_0424 = conj(IT_0101)*IT_0103;
    const ccomplex_t IT_0425 = conj(IT_0120)*IT_0122;
    const ccomplex_t IT_0426 = conj(IT_0130)*IT_0132;
    const ccomplex_t IT_0427 = conj(IT_0150)*IT_0153;
    const ccomplex_t IT_0428 = conj(IT_0160)*IT_0162;
    const ccomplex_t IT_0429 = conj(IT_0164)*IT_0165;
    const ccomplex_t IT_0430 = IT_0192*conj(IT_0302);
    const ccomplex_t IT_0431 = conj(IT_0256)*IT_0356;
    const ccomplex_t IT_0432 = conj(IT_0300)*IT_0387;
    const ccomplex_t IT_0433 = conj(IT_0181)*IT_0407;
    const ccomplex_t IT_0434 = IT_0045*IT_0311;
    const ccomplex_t IT_0435 = IT_0131*IT_0308;
    const ccomplex_t IT_0436 = IT_0162*conj(IT_0209) + IT_0368*IT_0400 + conj
      (IT_0174)*IT_0405 + conj(IT_0295)*IT_0408 + conj(IT_0290)*IT_0410 +
       IT_0420 + IT_0421 + IT_0422 + IT_0423 + IT_0424 + IT_0425 + IT_0426 +
       IT_0427 + IT_0428 + IT_0429 + IT_0430 + IT_0431 + IT_0432 + IT_0433 +
       IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = conj(IT_0300)*IT_0359;
    const ccomplex_t IT_0438 = conj(IT_0256)*IT_0386;
    const ccomplex_t IT_0439 = conj(IT_0277)*IT_0361;
    const ccomplex_t IT_0440 = conj(IT_0298)*IT_0382;
    const ccomplex_t IT_0441 = IT_0162*conj(IT_0209) + IT_0368*IT_0400 + conj
      (IT_0174)*IT_0405 + conj(IT_0295)*IT_0408 + conj(IT_0290)*IT_0410 +
       IT_0422 + IT_0423 + IT_0424 + IT_0425 + IT_0426 + IT_0427 + IT_0428 +
       IT_0429 + IT_0430 + IT_0433 + IT_0434 + IT_0435 + IT_0437 + IT_0438 +
       IT_0439 + IT_0440;
    const ccomplex_t IT_0442 = IT_0046*conj(IT_0174);
    const ccomplex_t IT_0443 = IT_0084*conj(IT_0164);
    const ccomplex_t IT_0444 = conj(IT_0181)*IT_0183;
    const ccomplex_t IT_0445 = IT_0196*conj(IT_0295);
    const ccomplex_t IT_0446 = IT_0132*conj(IT_0290);
    const ccomplex_t IT_0447 = IT_0153*conj(IT_0302);
    const ccomplex_t IT_0448 = IT_0162*conj(IT_0194);
    const ccomplex_t IT_0449 = conj(IT_0080)*IT_0165;
    const ccomplex_t IT_0450 = conj(IT_0150)*IT_0192;
    const ccomplex_t IT_0451 = 24*IT_0188;
    const ccomplex_t IT_0452 = conj(IT_0044)*IT_0451;
    const ccomplex_t IT_0453 = 24*IT_0182;
    const ccomplex_t IT_0454 = conj(IT_0101)*IT_0453;
    const ccomplex_t IT_0455 = 24*IT_0195;
    const ccomplex_t IT_0456 = conj(IT_0130)*IT_0410;
    const ccomplex_t IT_0457 = IT_0131*IT_0212;
    const ccomplex_t IT_0458 = IT_0045*IT_0217 + IT_0162*conj(IT_0305) +
       IT_0420 + IT_0421 + IT_0431 + IT_0432 + IT_0442 + IT_0443 + IT_0444 +
       IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 + IT_0452 +
       IT_0454 + conj(IT_0120)*IT_0455 + IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = IT_0045*IT_0217 + IT_0162*conj(IT_0305) +
       IT_0189*IT_0323 + IT_0437 + IT_0438 + IT_0439 + IT_0440 + IT_0442 +
       IT_0443 + IT_0444 + IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 +
       IT_0450 + IT_0452 + IT_0454 + conj(IT_0120)*IT_0455 + IT_0456 + IT_0457;
    const ccomplex_t IT_0460 = IT_0084*conj(IT_0302);
    const ccomplex_t IT_0461 = IT_0183*conj(IT_0295);
    const ccomplex_t IT_0462 = conj(IT_0181)*IT_0196;
    const ccomplex_t IT_0463 = IT_0132*conj(IT_0174);
    const ccomplex_t IT_0464 = IT_0153*conj(IT_0164);
    const ccomplex_t IT_0465 = conj(IT_0080)*IT_0192;
    const ccomplex_t IT_0466 = conj(IT_0101)*IT_0455;
    const ccomplex_t IT_0467 = conj(IT_0044)*IT_0410;
    const ccomplex_t IT_0468 = IT_0045*IT_0212;
    const ccomplex_t IT_0469 = conj(IT_0150)*IT_0165;
    const ccomplex_t IT_0470 = IT_0046*conj(IT_0290);
    const ccomplex_t IT_0471 = IT_0162*IT_0323;
    const ccomplex_t IT_0472 = IT_0189*IT_0322;
    const ccomplex_t IT_0473 = IT_0131*IT_0217 + IT_0392 + IT_0393 + IT_0403 +
       IT_0404 + conj(IT_0130)*IT_0451 + conj(IT_0120)*IT_0453 + IT_0460 +
       IT_0461 + IT_0462 + IT_0463 + IT_0464 + IT_0465 + IT_0466 + IT_0467 +
       IT_0468 + IT_0469 + IT_0470 + IT_0471 + IT_0472;
    const ccomplex_t IT_0474 = IT_0131*IT_0217 + IT_0415 + IT_0416 + IT_0417 +
       IT_0418 + conj(IT_0130)*IT_0451 + conj(IT_0120)*IT_0453 + IT_0460 +
       IT_0461 + IT_0462 + IT_0463 + IT_0464 + IT_0465 + IT_0466 + IT_0467 +
       IT_0468 + IT_0469 + IT_0470 + IT_0471 + IT_0472;
    const ccomplex_t IT_0475 = conj(IT_0023)*(IT_0047 + IT_0085 + IT_0104 +
       IT_0123 + IT_0133 + IT_0154 + IT_0163 + IT_0166) + conj(IT_0167)*(IT_0175
       + IT_0176 + IT_0184 + IT_0190 + IT_0191 + IT_0193) + conj(IT_0194)*
      (IT_0197 + IT_0198 + IT_0199 + IT_0200 + IT_0201) + conj(IT_0160)*(IT_0202
       + IT_0203 + IT_0204 + IT_0205 + IT_0206 + IT_0207 + IT_0208) + IT_0181
      *IT_0304 + IT_0101*IT_0320 + IT_0120*IT_0321 + IT_0295*IT_0324 + IT_0174
      *IT_0343 + IT_0290*IT_0344 + IT_0044*IT_0351 + IT_0130*IT_0352 + conj
      (IT_0215)*(IT_0047 + IT_0085 + IT_0104 + IT_0123 + IT_0133 + IT_0154 +
       IT_0163 + IT_0166 + IT_0355 + IT_0357) + conj(IT_0187)*(IT_0175 + IT_0176
       + IT_0184 + IT_0190 + IT_0191 + IT_0193 + IT_0355 + IT_0357) + conj
      (IT_0305)*(IT_0197 + IT_0198 + IT_0199 + IT_0200 + IT_0201 + IT_0360 +
       IT_0362) + conj(IT_0209)*(IT_0202 + IT_0203 + IT_0204 + IT_0205 + IT_0206
       + IT_0207 + IT_0208 + IT_0360 + IT_0362) + IT_0302*IT_0370 + IT_0150
      *IT_0371 + IT_0080*IT_0372 + IT_0164*IT_0373 + IT_0298*IT_0385 + IT_0256
      *IT_0388 + IT_0277*IT_0390 + IT_0300*IT_0391 + IT_0160*IT_0414 + IT_0209
      *IT_0419 + IT_0215*IT_0436 + IT_0023*IT_0441 + IT_0187*IT_0458 + IT_0167
      *IT_0459 + IT_0194*IT_0473 + IT_0305*IT_0474;
    return create_ccomplex_return(IT_0475);
}

