#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0015 = IT_0003*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0014;
    const ccomplex_t IT_0017 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = IT_0007*IT_0017;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + IT_0016 + 
      -IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = pow(m_Z, 2);
    const ccomplex_t IT_0024 = pow(m_N_1, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0022*IT_0026;
    const ccomplex_t IT_0028 = (-2)*IT_0027;
    const ccomplex_t IT_0029 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = IT_0007*IT_0029;
    const ccomplex_t IT_0032 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = IT_0007*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0030 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0037;
    const ccomplex_t IT_0040 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0041 = IT_0003*IT_0040;
    const ccomplex_t IT_0042 = IT_0007*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0038 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0048;
    const ccomplex_t IT_0050 = pow(m_Z, -2);
    const ccomplex_t IT_0051 = s_14*s_24*IT_0050;
    const ccomplex_t IT_0052 = s_12*IT_0023;
    const ccomplex_t IT_0053 = s_13*s_23;
    const ccomplex_t IT_0054 = (-2)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = (s_12 + -IT_0051)*(s_12 + (-0.5)*IT_0056);
    const ccomplex_t IT_0058 = (-2)*IT_0057;
    const ccomplex_t IT_0059 = 2*IT_0058;
    const ccomplex_t IT_0060 = conj(IT_0049)*IT_0059;
    const ccomplex_t IT_0061 = IT_0028*IT_0060;
    const ccomplex_t IT_0062 = -IT_0048;
    const ccomplex_t IT_0063 = s_14*IT_0024;
    const ccomplex_t IT_0064 = pow(s_23, 2);
    const ccomplex_t IT_0065 = s_14*IT_0064;
    const ccomplex_t IT_0066 = s_12*s_24*IT_0023;
    const ccomplex_t IT_0067 = s_12*s_23*s_34;
    const ccomplex_t IT_0068 = s_14*IT_0023*IT_0024;
    const ccomplex_t IT_0069 = 0.5*IT_0066 + -IT_0067 + (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0065 + IT_0069;
    const ccomplex_t IT_0071 = IT_0050*IT_0070;
    const ccomplex_t IT_0072 = s_12*s_24;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0063 + IT_0071 + IT_0073;
    const ccomplex_t IT_0075 = 4*IT_0074;
    const ccomplex_t IT_0076 = IT_0028*IT_0075;
    const ccomplex_t IT_0077 = conj(IT_0062)*IT_0076;
    const ccomplex_t IT_0078 = s_13*s_23*s_24;
    const ccomplex_t IT_0079 = s_13*s_34*IT_0024;
    const ccomplex_t IT_0080 = -IT_0068 + (-2)*IT_0078 + 2*IT_0079;
    const ccomplex_t IT_0081 = IT_0066 + IT_0080;
    const ccomplex_t IT_0082 = IT_0050*IT_0081;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0063;
    const ccomplex_t IT_0085 = IT_0072 + IT_0083 + IT_0084;
    const ccomplex_t IT_0086 = 4*IT_0085;
    const ccomplex_t IT_0087 = IT_0028*IT_0086;
    const ccomplex_t IT_0088 = conj(IT_0027)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0045*IT_0090;
    const ccomplex_t IT_0092 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0093 = IT_0003*IT_0092;
    const ccomplex_t IT_0094 = IT_0007*IT_0092;
    const ccomplex_t IT_0095 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = IT_0007*IT_0095;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0093 + IT_0094 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0101 = IT_0003*IT_0100;
    const ccomplex_t IT_0102 = IT_0007*IT_0100;
    const ccomplex_t IT_0103 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0104 = IT_0003*IT_0103;
    const ccomplex_t IT_0105 = IT_0007*IT_0103;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0101 + IT_0102 + 
      -IT_0104 + -IT_0105);
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = IT_0091 + IT_0111;
    const ccomplex_t IT_0113 = IT_0087*conj(IT_0112);
    const ccomplex_t IT_0114 = conj(IT_0028)*IT_0086;
    const ccomplex_t IT_0115 = IT_0027*IT_0114;
    const ccomplex_t IT_0116 = IT_0112*IT_0114;
    const ccomplex_t IT_0117 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0118 = IT_0003*IT_0117;
    const ccomplex_t IT_0119 = IT_0007*IT_0117;
    const ccomplex_t IT_0120 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0121 = IT_0003*IT_0120;
    const ccomplex_t IT_0122 = IT_0007*IT_0120;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0118 + IT_0119 + 
      -IT_0121 + -IT_0122);
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = cpow(IT_0124, 2);
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0023 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = (-2)*IT_0091 + (-2)*IT_0111 + 2*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_12 + (-2)*IT_0024 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0131 = cos(beta);
    const ccomplex_t IT_0132 = cos(alpha);
    const ccomplex_t IT_0133 = cpow(IT_0000, -4);
    const ccomplex_t IT_0134 = IT_0002*IT_0131*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = sin(beta);
    const ccomplex_t IT_0136 = sin(alpha);
    const ccomplex_t IT_0137 = IT_0002*IT_0133*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = cpow(IT_0002, -3);
    const ccomplex_t IT_0139 = IT_0131*IT_0132*IT_0138;
    const ccomplex_t IT_0140 = IT_0135*IT_0136*IT_0138;
    const ccomplex_t IT_0141 = cpow(IT_0000, -2);
    const ccomplex_t IT_0142 = IT_0006*IT_0131*IT_0132*IT_0141;
    const ccomplex_t IT_0143 = IT_0006*IT_0135*IT_0136*IT_0141;
    const ccomplex_t IT_0144 = cpow(IT_0002, -2);
    const ccomplex_t IT_0145 = cpow(IT_0141 + IT_0144, -1);
    const ccomplex_t IT_0146 = (IT_0134 + IT_0137 + IT_0139 + IT_0140 + 2
      *IT_0142 + 2*IT_0143)*IT_0145;
    const ccomplex_t IT_0147 = m_W*e_em;
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0130*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0153 = IT_0001*IT_0132;
    const ccomplex_t IT_0154 = IT_0152*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0156 = IT_0001*IT_0136;
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0159 = IT_0006*IT_0132;
    const ccomplex_t IT_0160 = IT_0158*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0162 = IT_0006*IT_0136;
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0154 + -IT_0157 + 
      -IT_0160 + IT_0163);
    const ccomplex_t IT_0165 = IT_0151*IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0167 = IT_0002*IT_0132*IT_0133*IT_0135;
    const ccomplex_t IT_0168 = IT_0002*IT_0131*IT_0133*IT_0136;
    const ccomplex_t IT_0169 = IT_0132*IT_0135*IT_0138;
    const ccomplex_t IT_0170 = IT_0131*IT_0136*IT_0138;
    const ccomplex_t IT_0171 = IT_0006*IT_0132*IT_0135*IT_0141;
    const ccomplex_t IT_0172 = IT_0006*IT_0131*IT_0136*IT_0141;
    const ccomplex_t IT_0173 = IT_0145*(IT_0167 + -IT_0168 + IT_0169 + 
      -IT_0170 + 2*IT_0171 + (-2)*IT_0172);
    const ccomplex_t IT_0174 = IT_0147*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0166*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0153*IT_0155;
    const ccomplex_t IT_0179 = IT_0152*IT_0156;
    const ccomplex_t IT_0180 = IT_0159*IT_0161;
    const ccomplex_t IT_0181 = IT_0158*IT_0162;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0178 + IT_0179 + 
      -IT_0180 + -IT_0181);
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = IT_0177*IT_0183;
    const ccomplex_t IT_0185 = -IT_0165 + -IT_0184;
    const ccomplex_t IT_0186 = pow(m_Z, -4);
    const ccomplex_t IT_0187 = s_23*s_24*s_34*IT_0186;
    const ccomplex_t IT_0188 = IT_0024 + IT_0187;
    const ccomplex_t IT_0189 = IT_0050*IT_0064;
    const ccomplex_t IT_0190 = pow(s_24, 2);
    const ccomplex_t IT_0191 = IT_0050*IT_0190;
    const ccomplex_t IT_0192 = -IT_0189 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0188 + IT_0192;
    const ccomplex_t IT_0194 = m_N_1*IT_0193;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = 2*IT_0195;
    const ccomplex_t IT_0197 = IT_0185*IT_0196;
    const ccomplex_t IT_0198 = IT_0049*IT_0059;
    const ccomplex_t IT_0199 = -IT_0191;
    const ccomplex_t IT_0200 = IT_0024 + IT_0199;
    const ccomplex_t IT_0201 = IT_0024*IT_0200;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = 6*IT_0202;
    const ccomplex_t IT_0204 = IT_0028*IT_0203;
    const ccomplex_t IT_0205 = -IT_0051;
    const ccomplex_t IT_0206 = s_12 + IT_0205;
    const ccomplex_t IT_0207 = IT_0024*IT_0206;
    const ccomplex_t IT_0208 = (-6)*IT_0207;
    const ccomplex_t IT_0209 = IT_0049*IT_0208;
    const ccomplex_t IT_0210 = (-0.5)*IT_0123;
    const ccomplex_t IT_0211 = cpow(IT_0210, 2);
    const ccomplex_t IT_0212 = IT_0127*IT_0211;
    const ccomplex_t IT_0213 = (-2)*IT_0091 + (-2)*IT_0111 + 2*IT_0212;
    const ccomplex_t IT_0214 = cpow(IT_0013, 2);
    const ccomplex_t IT_0215 = m_N_4*IT_0214;
    const ccomplex_t IT_0216 = m_N_1*IT_0022;
    const ccomplex_t IT_0217 = IT_0026*(IT_0215 + -IT_0216);
    const ccomplex_t IT_0218 = IT_0023*IT_0024;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0190 + IT_0219;
    const ccomplex_t IT_0221 = IT_0050*IT_0220;
    const ccomplex_t IT_0222 = -IT_0024;
    const ccomplex_t IT_0223 = IT_0191 + IT_0222;
    const ccomplex_t IT_0224 = (-4)*IT_0223;
    const ccomplex_t IT_0225 = IT_0221 + IT_0224;
    const ccomplex_t IT_0226 = m_N_1*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0226;
    const ccomplex_t IT_0228 = IT_0217*IT_0227;
    const ccomplex_t IT_0229 = cpow(IT_0036, 2);
    const ccomplex_t IT_0230 = m_N_3*IT_0229;
    const ccomplex_t IT_0231 = m_N_1*IT_0045;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = IT_0047*IT_0233;
    const ccomplex_t IT_0235 = cpow(IT_0099, 2);
    const ccomplex_t IT_0236 = m_N_2*IT_0235;
    const ccomplex_t IT_0237 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = IT_0236*IT_0238;
    const ccomplex_t IT_0240 = IT_0234 + IT_0239;
    const ccomplex_t IT_0241 = s_23*s_24*s_34*IT_0050;
    const ccomplex_t IT_0242 = 0.5*IT_0190 + (-0.5)*IT_0218 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0064 + IT_0242;
    const ccomplex_t IT_0244 = IT_0050*IT_0243;
    const ccomplex_t IT_0245 = IT_0024 + IT_0199 + IT_0244;
    const ccomplex_t IT_0246 = m_N_1*IT_0245;
    const ccomplex_t IT_0247 = (-4)*IT_0246;
    const ccomplex_t IT_0248 = IT_0240*IT_0247;
    const ccomplex_t IT_0249 = s_12*s_14;
    const ccomplex_t IT_0250 = s_12*s_14*IT_0023;
    const ccomplex_t IT_0251 = s_24*IT_0023*IT_0024;
    const ccomplex_t IT_0252 = pow(s_13, 2);
    const ccomplex_t IT_0253 = s_24*IT_0252;
    const ccomplex_t IT_0254 = s_12*s_13*s_34;
    const ccomplex_t IT_0255 = -IT_0251 + 2*IT_0253 + (-2)*IT_0254;
    const ccomplex_t IT_0256 = IT_0250 + IT_0255;
    const ccomplex_t IT_0257 = IT_0050*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = s_24*IT_0024;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = IT_0249 + IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = (-4)*IT_0261;
    const ccomplex_t IT_0263 = IT_0049*IT_0262;
    const ccomplex_t IT_0264 = conj(IT_0112)*IT_0263;
    const ccomplex_t IT_0265 = conj(IT_0049)*IT_0262;
    const ccomplex_t IT_0266 = IT_0027*IT_0265;
    const ccomplex_t IT_0267 = IT_0112*IT_0265;
    const ccomplex_t IT_0268 = s_13*s_14*s_23;
    const ccomplex_t IT_0269 = s_23*s_34*IT_0024;
    const ccomplex_t IT_0270 = -IT_0251 + (-2)*IT_0268 + 2*IT_0269;
    const ccomplex_t IT_0271 = IT_0250 + IT_0270;
    const ccomplex_t IT_0272 = IT_0050*IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0259;
    const ccomplex_t IT_0274 = (-2)*IT_0249;
    const ccomplex_t IT_0275 = IT_0272 + IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = (-2)*IT_0275;
    const ccomplex_t IT_0277 = conj(IT_0049)*IT_0276;
    const ccomplex_t IT_0278 = IT_0062*IT_0277;
    const ccomplex_t IT_0279 = IT_0108*IT_0238;
    const ccomplex_t IT_0280 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0022*IT_0281;
    const ccomplex_t IT_0283 = cpow((-2)*s_23 + IT_0023 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0125*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0279 + 2*IT_0282 + (-2)*IT_0285;
    const ccomplex_t IT_0287 = IT_0028*IT_0059;
    const ccomplex_t IT_0288 = IT_0028*IT_0208;
    const ccomplex_t IT_0289 = s_14*s_23*s_34*IT_0186;
    const ccomplex_t IT_0290 = s_13*s_23*IT_0050;
    const ccomplex_t IT_0291 = -IT_0051 + -IT_0290;
    const ccomplex_t IT_0292 = s_12 + IT_0289 + IT_0291;
    const ccomplex_t IT_0293 = m_N_1*IT_0292;
    const ccomplex_t IT_0294 = (-2)*IT_0293;
    const ccomplex_t IT_0295 = IT_0185*IT_0294;
    const ccomplex_t IT_0296 = pow(s_14, 2);
    const ccomplex_t IT_0297 = IT_0050*IT_0296;
    const ccomplex_t IT_0298 = (s_12 + (-0.5)*IT_0056)*(IT_0024 + -IT_0297);
    const ccomplex_t IT_0299 = 2*IT_0298;
    const ccomplex_t IT_0300 = (-2)*IT_0299;
    const ccomplex_t IT_0301 = IT_0049*conj(IT_0049)*IT_0300;
    const ccomplex_t IT_0302 = (s_12 + (-0.5)*IT_0056)*(IT_0024 + -IT_0191);
    const ccomplex_t IT_0303 = 2*IT_0302;
    const ccomplex_t IT_0304 = (-2)*IT_0303;
    const ccomplex_t IT_0305 = IT_0028*conj(IT_0028)*IT_0304;
    const ccomplex_t IT_0306 = cpow(IT_0044, 2);
    const ccomplex_t IT_0307 = m_N_3*IT_0306;
    const ccomplex_t IT_0308 = IT_0232 + IT_0307;
    const ccomplex_t IT_0309 = IT_0047*IT_0308;
    const ccomplex_t IT_0310 = cpow(IT_0107, 2);
    const ccomplex_t IT_0311 = m_N_2*IT_0310;
    const ccomplex_t IT_0312 = IT_0238*IT_0311;
    const ccomplex_t IT_0313 = IT_0309 + IT_0312;
    const ccomplex_t IT_0314 = s_14*s_24;
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = IT_0052 + IT_0315;
    const ccomplex_t IT_0317 = IT_0050*IT_0316;
    const ccomplex_t IT_0318 = -IT_0051 + (-0.25)*IT_0317;
    const ccomplex_t IT_0319 = s_12 + IT_0318;
    const ccomplex_t IT_0320 = m_N_1*IT_0319;
    const ccomplex_t IT_0321 = (-4)*IT_0320;
    const ccomplex_t IT_0322 = 2*IT_0321;
    const ccomplex_t IT_0323 = conj(IT_0028)*IT_0322;
    const ccomplex_t IT_0324 = conj(IT_0049) + conj(IT_0286);
    const ccomplex_t IT_0325 = s_14*s_23*s_34*IT_0050;
    const ccomplex_t IT_0326 = (-2)*IT_0053 + -IT_0314 + 2*IT_0325;
    const ccomplex_t IT_0327 = IT_0052 + IT_0326;
    const ccomplex_t IT_0328 = IT_0050*IT_0327;
    const ccomplex_t IT_0329 = -IT_0051 + (-0.5)*IT_0328;
    const ccomplex_t IT_0330 = s_12 + IT_0329;
    const ccomplex_t IT_0331 = m_N_1*IT_0330;
    const ccomplex_t IT_0332 = (-2)*IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0332;
    const ccomplex_t IT_0334 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0335 = IT_0153*IT_0334;
    const ccomplex_t IT_0336 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0337 = IT_0156*IT_0336;
    const ccomplex_t IT_0338 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0339 = IT_0159*IT_0338;
    const ccomplex_t IT_0340 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0341 = IT_0162*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*(IT_0335 + -IT_0337 + 
      -IT_0339 + IT_0341);
    const ccomplex_t IT_0343 = IT_0151*IT_0342;
    const ccomplex_t IT_0344 = IT_0153*IT_0336;
    const ccomplex_t IT_0345 = IT_0156*IT_0334;
    const ccomplex_t IT_0346 = IT_0159*IT_0340;
    const ccomplex_t IT_0347 = IT_0162*IT_0338;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*(IT_0344 + IT_0345 + 
      -IT_0346 + -IT_0347);
    const ccomplex_t IT_0349 = -IT_0348;
    const ccomplex_t IT_0350 = IT_0177*IT_0349;
    const ccomplex_t IT_0351 = -IT_0343 + -IT_0350;
    const ccomplex_t IT_0352 = pow(s_34, 2);
    const ccomplex_t IT_0353 = IT_0186*IT_0352;
    const ccomplex_t IT_0354 = 0.5*IT_0353;
    const ccomplex_t IT_0355 = 1 + IT_0354;
    const ccomplex_t IT_0356 = IT_0024*IT_0355;
    const ccomplex_t IT_0357 = 4*IT_0356;
    const ccomplex_t IT_0358 = (-0.5)*IT_0056;
    const ccomplex_t IT_0359 = s_12 + IT_0358;
    const ccomplex_t IT_0360 = (-2)*IT_0314;
    const ccomplex_t IT_0361 = IT_0052 + IT_0360;
    const ccomplex_t IT_0362 = IT_0050*IT_0361;
    const ccomplex_t IT_0363 = (-0.5)*IT_0362;
    const ccomplex_t IT_0364 = pow(m_Z, 4);
    const ccomplex_t IT_0365 = s_12*IT_0364;
    const ccomplex_t IT_0366 = s_13*s_23*IT_0023;
    const ccomplex_t IT_0367 = s_14*s_24*IT_0023;
    const ccomplex_t IT_0368 = s_14*s_23*s_34;
    const ccomplex_t IT_0369 = (-2)*IT_0366 + (-2)*IT_0367 + 4*IT_0368;
    const ccomplex_t IT_0370 = IT_0365 + IT_0369;
    const ccomplex_t IT_0371 = IT_0186*IT_0370;
    const ccomplex_t IT_0372 = 0.25*IT_0371;
    const ccomplex_t IT_0373 = IT_0359 + IT_0363 + IT_0372;
    const ccomplex_t IT_0374 = 8*IT_0373;
    const ccomplex_t IT_0375 = cpow(IT_0021, 2);
    const ccomplex_t IT_0376 = m_N_4*IT_0375;
    const ccomplex_t IT_0377 = -IT_0026*(IT_0216 + -IT_0376);
    const ccomplex_t IT_0378 = s_13*s_24*s_34;
    const ccomplex_t IT_0379 = s_12*IT_0352;
    const ccomplex_t IT_0380 = (-2)*IT_0366 + (-2)*IT_0367 + 2*IT_0368 + 2
      *IT_0378 + (-2)*IT_0379;
    const ccomplex_t IT_0381 = IT_0365 + IT_0380;
    const ccomplex_t IT_0382 = IT_0186*IT_0381;
    const ccomplex_t IT_0383 = (-0.25)*IT_0382;
    const ccomplex_t IT_0384 = -IT_0290;
    const ccomplex_t IT_0385 = s_12 + IT_0205 + IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = 8*IT_0385;
    const ccomplex_t IT_0387 = (-2)*IT_0352;
    const ccomplex_t IT_0388 = IT_0364 + IT_0387;
    const ccomplex_t IT_0389 = IT_0186*IT_0388;
    const ccomplex_t IT_0390 = 0.25*IT_0389;
    const ccomplex_t IT_0391 = 1 + IT_0390;
    const ccomplex_t IT_0392 = IT_0024*IT_0391;
    const ccomplex_t IT_0393 = 8*IT_0392;
    const ccomplex_t IT_0394 = IT_0124*IT_0210;
    const ccomplex_t IT_0395 = m_N_1*IT_0394;
    const ccomplex_t IT_0396 = m_N_1*IT_0211;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = IT_0395 + IT_0397;
    const ccomplex_t IT_0399 = IT_0284*IT_0398;
    const ccomplex_t IT_0400 = m_N_1*IT_0108;
    const ccomplex_t IT_0401 = IT_0238*IT_0400;
    const ccomplex_t IT_0402 = -IT_0216;
    const ccomplex_t IT_0403 = IT_0215 + IT_0402;
    const ccomplex_t IT_0404 = IT_0281*IT_0403;
    const ccomplex_t IT_0405 = IT_0399 + IT_0401 + -IT_0404;
    const ccomplex_t IT_0406 = 1.125*IT_0024;
    const ccomplex_t IT_0407 = 16*IT_0406;
    const ccomplex_t IT_0408 = (-16)*IT_0406;
    const ccomplex_t IT_0409 = IT_0050*IT_0052;
    const ccomplex_t IT_0410 = (-0.5)*IT_0409;
    const ccomplex_t IT_0411 = s_34*IT_0186;
    const ccomplex_t IT_0412 = s_14*s_23;
    const ccomplex_t IT_0413 = s_12*s_34;
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = s_13*s_24;
    const ccomplex_t IT_0416 = -IT_0415;
    const ccomplex_t IT_0417 = IT_0414 + IT_0416;
    const ccomplex_t IT_0418 = IT_0411*IT_0417;
    const ccomplex_t IT_0419 = 0.25*IT_0418;
    const ccomplex_t IT_0420 = s_12 + IT_0410 + IT_0419;
    const ccomplex_t IT_0421 = (-8)*IT_0420;
    const ccomplex_t IT_0422 = m_N_1*IT_0125;
    const ccomplex_t IT_0423 = -IT_0422;
    const ccomplex_t IT_0424 = IT_0395 + IT_0423;
    const ccomplex_t IT_0425 = IT_0284*IT_0424;
    const ccomplex_t IT_0426 = IT_0376 + IT_0402;
    const ccomplex_t IT_0427 = IT_0281*IT_0426;
    const ccomplex_t IT_0428 = IT_0401 + IT_0425 + -IT_0427;
    const ccomplex_t IT_0429 = (-8)*IT_0373;
    const ccomplex_t IT_0430 = IT_0127*IT_0424;
    const ccomplex_t IT_0431 = -IT_0400;
    const ccomplex_t IT_0432 = IT_0236 + IT_0431;
    const ccomplex_t IT_0433 = IT_0110*IT_0432;
    const ccomplex_t IT_0434 = IT_0090*IT_0233;
    const ccomplex_t IT_0435 = IT_0430 + -IT_0433 + -IT_0434;
    const ccomplex_t IT_0436 = (-8)*IT_0392;
    const ccomplex_t IT_0437 = conj(IT_0028) + conj(IT_0129);
    const ccomplex_t IT_0438 = IT_0211*IT_0284;
    const ccomplex_t IT_0439 = 2*IT_0279 + 2*IT_0282 + (-2)*IT_0438;
    const ccomplex_t IT_0440 = IT_0219 + IT_0296;
    const ccomplex_t IT_0441 = IT_0050*IT_0440;
    const ccomplex_t IT_0442 = IT_0222 + IT_0297;
    const ccomplex_t IT_0443 = (-4)*IT_0442;
    const ccomplex_t IT_0444 = IT_0441 + IT_0443;
    const ccomplex_t IT_0445 = m_N_1*IT_0444;
    const ccomplex_t IT_0446 = (-2)*IT_0445;
    const ccomplex_t IT_0447 = conj(IT_0049)*IT_0446;
    const ccomplex_t IT_0448 = conj(IT_0213)*IT_0322 + IT_0323 + IT_0324
      *IT_0333 + conj(IT_0351)*IT_0357 + conj(IT_0313)*IT_0374 + conj(IT_0377)
      *IT_0386 + conj(IT_0217)*IT_0393 + conj(IT_0405)*IT_0407 + conj(IT_0240)
      *IT_0408 + conj(IT_0185)*IT_0421 + conj(IT_0428)*IT_0429 + conj(IT_0435)
      *IT_0436 + IT_0247*IT_0437 + conj(IT_0439)*IT_0446 + IT_0447;
    const ccomplex_t IT_0449 = conj(IT_0185)*IT_0357;
    const ccomplex_t IT_0450 = (-8)*IT_0385;
    const ccomplex_t IT_0451 = conj(IT_0129)*IT_0322 + IT_0323 + conj(IT_0240)
      *IT_0374 + conj(IT_0217)*IT_0386 + conj(IT_0377)*IT_0393 + conj(IT_0313)
      *IT_0408 + conj(IT_0351)*IT_0421 + IT_0407*conj(IT_0428) + conj(IT_0405)
      *IT_0429 + conj(IT_0286)*IT_0446 + IT_0447 + IT_0449 + conj(IT_0435)
      *IT_0450;
    const ccomplex_t IT_0452 = conj(IT_0405) + conj(IT_0428);
    const ccomplex_t IT_0453 = 4*IT_0246;
    const ccomplex_t IT_0454 = IT_0028*IT_0453;
    const ccomplex_t IT_0455 = (-2)*IT_0321;
    const ccomplex_t IT_0456 = IT_0028*IT_0455;
    const ccomplex_t IT_0457 = IT_0240*IT_0333;
    const ccomplex_t IT_0458 = s_13*s_14*s_34*IT_0186;
    const ccomplex_t IT_0459 = IT_0024 + IT_0458;
    const ccomplex_t IT_0460 = IT_0050*IT_0252;
    const ccomplex_t IT_0461 = -IT_0297 + -IT_0460;
    const ccomplex_t IT_0462 = IT_0459 + IT_0461;
    const ccomplex_t IT_0463 = m_N_1*IT_0462;
    const ccomplex_t IT_0464 = -IT_0463;
    const ccomplex_t IT_0465 = (-2)*IT_0464;
    const ccomplex_t IT_0466 = IT_0185*IT_0465;
    const ccomplex_t IT_0467 = IT_0217*IT_0455;
    const ccomplex_t IT_0468 = s_13*s_24*s_34*IT_0186;
    const ccomplex_t IT_0469 = s_12 + IT_0291 + IT_0468;
    const ccomplex_t IT_0470 = m_N_1*IT_0469;
    const ccomplex_t IT_0471 = 2*IT_0470;
    const ccomplex_t IT_0472 = conj(IT_0185)*IT_0471;
    const ccomplex_t IT_0473 = IT_0196*conj(IT_0351);
    const ccomplex_t IT_0474 = -IT_0090*(IT_0231 + -IT_0307);
    const ccomplex_t IT_0475 = -IT_0474;
    const ccomplex_t IT_0476 = IT_0311 + IT_0431;
    const ccomplex_t IT_0477 = IT_0110*IT_0476;
    const ccomplex_t IT_0478 = IT_0127*IT_0398;
    const ccomplex_t IT_0479 = -IT_0477 + IT_0478;
    const ccomplex_t IT_0480 = conj(IT_0475) + conj(IT_0479);
    const ccomplex_t IT_0481 = s_13*s_24*s_34*IT_0050;
    const ccomplex_t IT_0482 = (-2)*IT_0053 + -IT_0314 + 2*IT_0481;
    const ccomplex_t IT_0483 = IT_0052 + IT_0482;
    const ccomplex_t IT_0484 = IT_0050*IT_0483;
    const ccomplex_t IT_0485 = -IT_0051 + (-0.5)*IT_0484;
    const ccomplex_t IT_0486 = s_12 + IT_0485;
    const ccomplex_t IT_0487 = m_N_1*IT_0486;
    const ccomplex_t IT_0488 = (-2)*IT_0487;
    const ccomplex_t IT_0489 = 2*IT_0488;
    const ccomplex_t IT_0490 = conj(IT_0217)*IT_0227;
    const ccomplex_t IT_0491 = conj(IT_0028)*IT_0203;
    const ccomplex_t IT_0492 = conj(IT_0049)*IT_0208;
    const ccomplex_t IT_0493 = conj(IT_0240)*IT_0247;
    const ccomplex_t IT_0494 = conj(IT_0028)*IT_0304;
    const ccomplex_t IT_0495 = (-2)*IT_0488;
    const ccomplex_t IT_0496 = (-2)*IT_0226;
    const ccomplex_t IT_0497 = IT_0060 + conj(IT_0129)*IT_0203 + IT_0208*conj
      (IT_0286) + conj(IT_0213)*IT_0304 + conj(IT_0313)*IT_0322 + IT_0059*conj
      (IT_0439) + conj(IT_0405)*IT_0453 + conj(IT_0428)*IT_0455 + IT_0472 +
       IT_0473 + IT_0480*IT_0489 + IT_0490 + IT_0491 + IT_0492 + IT_0493 +
       IT_0494 + conj(IT_0377)*IT_0495 + conj(IT_0435)*IT_0496;
    const ccomplex_t IT_0498 = IT_0049*IT_0322;
    const ccomplex_t IT_0499 = IT_0028*IT_0489;
    const ccomplex_t IT_0500 = s_13*s_14*s_34*IT_0050;
    const ccomplex_t IT_0501 = (-0.5)*IT_0218 + 0.5*IT_0296 + -IT_0500;
    const ccomplex_t IT_0502 = IT_0252 + IT_0501;
    const ccomplex_t IT_0503 = IT_0050*IT_0502;
    const ccomplex_t IT_0504 = -IT_0297;
    const ccomplex_t IT_0505 = IT_0024 + IT_0503 + IT_0504;
    const ccomplex_t IT_0506 = m_N_1*IT_0505;
    const ccomplex_t IT_0507 = (-4)*IT_0506;
    const ccomplex_t IT_0508 = IT_0049*IT_0507;
    const ccomplex_t IT_0509 = IT_0028*IT_0496;
    const ccomplex_t IT_0510 = IT_0286*IT_0322;
    const ccomplex_t IT_0511 = -IT_0413 + -IT_0415;
    const ccomplex_t IT_0512 = IT_0412 + IT_0511;
    const ccomplex_t IT_0513 = IT_0411*IT_0512;
    const ccomplex_t IT_0514 = (-0.25)*IT_0513;
    const ccomplex_t IT_0515 = s_12 + IT_0410 + IT_0514;
    const ccomplex_t IT_0516 = 8*IT_0515;
    const ccomplex_t IT_0517 = IT_0185*IT_0516;
    const ccomplex_t IT_0518 = (-2)*IT_0366 + (-2)*IT_0367 + 4*IT_0378;
    const ccomplex_t IT_0519 = IT_0365 + IT_0518;
    const ccomplex_t IT_0520 = IT_0186*IT_0519;
    const ccomplex_t IT_0521 = 0.25*IT_0520;
    const ccomplex_t IT_0522 = IT_0359 + IT_0363 + IT_0521;
    const ccomplex_t IT_0523 = 8*IT_0522;
    const ccomplex_t IT_0524 = IT_0475*IT_0523;
    const ccomplex_t IT_0525 = IT_0393*IT_0405;
    const ccomplex_t IT_0526 = IT_0217*IT_0407;
    const ccomplex_t IT_0527 = IT_0313*IT_0450;
    const ccomplex_t IT_0528 = IT_0240*IT_0436;
    const ccomplex_t IT_0529 = IT_0439*IT_0507;
    const ccomplex_t IT_0530 = IT_0129*IT_0496;
    const ccomplex_t IT_0531 = s_13*IT_0023;
    const ccomplex_t IT_0532 = IT_0411*IT_0531;
    const ccomplex_t IT_0533 = s_14*IT_0023;
    const ccomplex_t IT_0534 = IT_0050*IT_0533;
    const ccomplex_t IT_0535 = s_13*s_34;
    const ccomplex_t IT_0536 = (-2)*IT_0535;
    const ccomplex_t IT_0537 = IT_0533 + IT_0536;
    const ccomplex_t IT_0538 = IT_0050*IT_0537;
    const ccomplex_t IT_0539 = (-0.5)*IT_0532 + 0.5*IT_0534 + (-0.5)*IT_0538;
    const ccomplex_t IT_0540 = s_14 + IT_0539;
    const ccomplex_t IT_0541 = m_N_1*IT_0540;
    const ccomplex_t IT_0542 = (-2)*IT_0541;
    const ccomplex_t IT_0543 = 2*IT_0542;
    const ccomplex_t IT_0544 = conj(IT_0351)*IT_0543;
    const ccomplex_t IT_0545 = conj(IT_0185)*IT_0543;
    const ccomplex_t IT_0546 = s_24*IT_0023;
    const ccomplex_t IT_0547 = s_23*s_34;
    const ccomplex_t IT_0548 = (-2)*IT_0547;
    const ccomplex_t IT_0549 = IT_0546 + IT_0548;
    const ccomplex_t IT_0550 = IT_0050*IT_0549;
    const ccomplex_t IT_0551 = s_24*IT_0364;
    const ccomplex_t IT_0552 = s_23*s_34*IT_0023;
    const ccomplex_t IT_0553 = (-2)*IT_0552;
    const ccomplex_t IT_0554 = IT_0551 + IT_0553;
    const ccomplex_t IT_0555 = IT_0186*IT_0554;
    const ccomplex_t IT_0556 = (-0.666666666666667)*IT_0550 +
       0.166666666666667*IT_0555;
    const ccomplex_t IT_0557 = s_24 + IT_0556;
    const ccomplex_t IT_0558 = m_N_1*IT_0557;
    const ccomplex_t IT_0559 = 6*IT_0558;
    const ccomplex_t IT_0560 = 2*IT_0559;
    const ccomplex_t IT_0561 = conj(IT_0405)*IT_0560;
    const ccomplex_t IT_0562 = conj(IT_0428)*IT_0560;
    const ccomplex_t IT_0563 = s_13*s_34*IT_0050;
    const ccomplex_t IT_0564 = s_14*IT_0364;
    const ccomplex_t IT_0565 = s_13*s_34*IT_0023;
    const ccomplex_t IT_0566 = (-2)*IT_0565;
    const ccomplex_t IT_0567 = IT_0564 + IT_0566;
    const ccomplex_t IT_0568 = IT_0186*IT_0567;
    const ccomplex_t IT_0569 = (-0.666666666666667)*IT_0563 + (
      -0.166666666666667)*IT_0568;
    const ccomplex_t IT_0570 = s_14 + IT_0569;
    const ccomplex_t IT_0571 = m_N_1*IT_0570;
    const ccomplex_t IT_0572 = (-6)*IT_0571;
    const ccomplex_t IT_0573 = 2*IT_0572;
    const ccomplex_t IT_0574 = conj(IT_0377)*IT_0573;
    const ccomplex_t IT_0575 = conj(IT_0217)*IT_0573;
    const ccomplex_t IT_0576 = m_N_1*IT_0050;
    const ccomplex_t IT_0577 = IT_0050*IT_0551;
    const ccomplex_t IT_0578 = (-2)*IT_0547 + (-0.5)*IT_0577;
    const ccomplex_t IT_0579 = IT_0546 + IT_0578;
    const ccomplex_t IT_0580 = IT_0576*IT_0579;
    const ccomplex_t IT_0581 = (-2)*IT_0580;
    const ccomplex_t IT_0582 = 2*IT_0581;
    const ccomplex_t IT_0583 = conj(IT_0377)*IT_0582;
    const ccomplex_t IT_0584 = conj(IT_0217)*IT_0582;
    const ccomplex_t IT_0585 = IT_0186*IT_0564;
    const ccomplex_t IT_0586 = (-10)*s_14;
    const ccomplex_t IT_0587 = IT_0585 + IT_0586;
    const ccomplex_t IT_0588 = m_N_1*IT_0587;
    const ccomplex_t IT_0589 = 2*IT_0588;
    const ccomplex_t IT_0590 = conj(IT_0405)*IT_0589;
    const ccomplex_t IT_0591 = conj(IT_0428)*IT_0589;
    const ccomplex_t IT_0592 = IT_0075*conj(IT_0213);
    const ccomplex_t IT_0593 = IT_0075*conj(IT_0129);
    const ccomplex_t IT_0594 = conj(IT_0028)*IT_0075;
    const ccomplex_t IT_0595 = pow(m_Z, 6);
    const ccomplex_t IT_0596 = s_12*IT_0595;
    const ccomplex_t IT_0597 = s_13*s_23*IT_0364;
    const ccomplex_t IT_0598 = (-2)*IT_0597;
    const ccomplex_t IT_0599 = IT_0596 + IT_0598;
    const ccomplex_t IT_0600 = IT_0186*IT_0599;
    const ccomplex_t IT_0601 = 0.5*IT_0600;
    const ccomplex_t IT_0602 = (-4)*IT_0314;
    const ccomplex_t IT_0603 = IT_0052 + IT_0601 + IT_0602;
    const ccomplex_t IT_0604 = IT_0050*IT_0370;
    const ccomplex_t IT_0605 = -IT_0604;
    const ccomplex_t IT_0606 = IT_0603 + IT_0605;
    const ccomplex_t IT_0607 = 9*IT_0218;
    const ccomplex_t IT_0608 = (-4)*IT_0607;
    const ccomplex_t IT_0609 = conj(IT_0027) + conj(IT_0112);
    const ccomplex_t IT_0610 = 0.25*IT_0600;
    const ccomplex_t IT_0611 = IT_0052 + IT_0610;
    const ccomplex_t IT_0612 = (-16)*IT_0611;
    const ccomplex_t IT_0613 = IT_0024*IT_0050;
    const ccomplex_t IT_0614 = (-0.25)*IT_0364;
    const ccomplex_t IT_0615 = IT_0352 + IT_0614;
    const ccomplex_t IT_0616 = IT_0613*IT_0615;
    const ccomplex_t IT_0617 = (-4)*IT_0616;
    const ccomplex_t IT_0618 = (-4)*IT_0617;
    const ccomplex_t IT_0619 = IT_0276*conj(IT_0439);
    const ccomplex_t IT_0620 = IT_0276*conj(IT_0286);
    const ccomplex_t IT_0621 = s_23*IT_0023;
    const ccomplex_t IT_0622 = IT_0411*IT_0621;
    const ccomplex_t IT_0623 = IT_0050*IT_0546;
    const ccomplex_t IT_0624 = (-0.5)*IT_0550 + (-0.5)*IT_0622 + 0.5*IT_0623;
    const ccomplex_t IT_0625 = s_24 + IT_0624;
    const ccomplex_t IT_0626 = m_N_1*IT_0625;
    const ccomplex_t IT_0627 = (-2)*IT_0626;
    const ccomplex_t IT_0628 = (-2)*IT_0627;
    const ccomplex_t IT_0629 = conj(IT_0351)*IT_0628;
    const ccomplex_t IT_0630 = conj(IT_0185)*IT_0628;
    const ccomplex_t IT_0631 = conj(IT_0240) + conj(IT_0313);
    const ccomplex_t IT_0632 = 0.666666666666667*IT_0534 + (-0.166666666666667
      )*IT_0585;
    const ccomplex_t IT_0633 = s_14 + IT_0632;
    const ccomplex_t IT_0634 = m_N_1*IT_0633;
    const ccomplex_t IT_0635 = (-6)*IT_0634;
    const ccomplex_t IT_0636 = (-2)*IT_0635;
    const ccomplex_t IT_0637 = (-2)*IT_0559;
    const ccomplex_t IT_0638 = conj(IT_0240)*IT_0637;
    const ccomplex_t IT_0639 = conj(IT_0313)*IT_0637;
    const ccomplex_t IT_0640 = (-2)*IT_0572;
    const ccomplex_t IT_0641 = conj(IT_0435)*IT_0640;
    const ccomplex_t IT_0642 = conj(IT_0475)*IT_0640;
    const ccomplex_t IT_0643 = conj(IT_0479)*IT_0640;
    const ccomplex_t IT_0644 = (-2)*IT_0581;
    const ccomplex_t IT_0645 = conj(IT_0435)*IT_0644;
    const ccomplex_t IT_0646 = conj(IT_0475)*IT_0644;
    const ccomplex_t IT_0647 = conj(IT_0479)*IT_0644;
    const ccomplex_t IT_0648 = IT_0544 + IT_0545 + IT_0561 + IT_0562 + IT_0574
       + IT_0575 + IT_0583 + IT_0584 + IT_0590 + IT_0591 + IT_0592 + IT_0593 + 2
      *IT_0594 + (-8)*conj(IT_0062)*(IT_0606 + (-0.125)*IT_0608) + IT_0609*
      (IT_0612 + IT_0618) + IT_0619 + IT_0620 + IT_0629 + IT_0630 + IT_0631
      *IT_0636 + IT_0638 + IT_0639 + IT_0641 + IT_0642 + IT_0643 + IT_0645 +
       IT_0646 + IT_0647;
    const ccomplex_t IT_0649 = IT_0351*IT_0543;
    const ccomplex_t IT_0650 = IT_0185*IT_0543;
    const ccomplex_t IT_0651 = IT_0050*IT_0564;
    const ccomplex_t IT_0652 = (-2)*IT_0535 + (-0.5)*IT_0651;
    const ccomplex_t IT_0653 = IT_0533 + IT_0652;
    const ccomplex_t IT_0654 = IT_0576*IT_0653;
    const ccomplex_t IT_0655 = (-2)*IT_0654;
    const ccomplex_t IT_0656 = 2*IT_0655;
    const ccomplex_t IT_0657 = IT_0405*IT_0656;
    const ccomplex_t IT_0658 = IT_0428*IT_0656;
    const ccomplex_t IT_0659 = s_23*s_34*IT_0050;
    const ccomplex_t IT_0660 = (-0.166666666666667)*IT_0555 + (
      -0.666666666666667)*IT_0659;
    const ccomplex_t IT_0661 = s_24 + IT_0660;
    const ccomplex_t IT_0662 = m_N_1*IT_0661;
    const ccomplex_t IT_0663 = (-6)*IT_0662;
    const ccomplex_t IT_0664 = 2*IT_0663;
    const ccomplex_t IT_0665 = IT_0405*IT_0664;
    const ccomplex_t IT_0666 = IT_0428*IT_0664;
    const ccomplex_t IT_0667 = (-0.666666666666667)*IT_0538 +
       0.166666666666667*IT_0568;
    const ccomplex_t IT_0668 = s_14 + IT_0667;
    const ccomplex_t IT_0669 = m_N_1*IT_0668;
    const ccomplex_t IT_0670 = 6*IT_0669;
    const ccomplex_t IT_0671 = 2*IT_0670;
    const ccomplex_t IT_0672 = IT_0377*IT_0671;
    const ccomplex_t IT_0673 = IT_0217*IT_0671;
    const ccomplex_t IT_0674 = IT_0186*IT_0551;
    const ccomplex_t IT_0675 = (-10)*s_24;
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = m_N_1*IT_0676;
    const ccomplex_t IT_0678 = 2*IT_0677;
    const ccomplex_t IT_0679 = IT_0377*IT_0678;
    const ccomplex_t IT_0680 = IT_0217*IT_0678;
    const ccomplex_t IT_0681 = IT_0086*IT_0213;
    const ccomplex_t IT_0682 = IT_0086*IT_0129;
    const ccomplex_t IT_0683 = IT_0262*IT_0439;
    const ccomplex_t IT_0684 = IT_0262*IT_0286;
    const ccomplex_t IT_0685 = IT_0351*IT_0628;
    const ccomplex_t IT_0686 = IT_0185*IT_0628;
    const ccomplex_t IT_0687 = (-2)*IT_0655;
    const ccomplex_t IT_0688 = IT_0240*IT_0687;
    const ccomplex_t IT_0689 = IT_0313*IT_0687;
    const ccomplex_t IT_0690 = (-2)*IT_0663;
    const ccomplex_t IT_0691 = IT_0240*IT_0690;
    const ccomplex_t IT_0692 = IT_0313*IT_0690;
    const ccomplex_t IT_0693 = (-2)*IT_0670;
    const ccomplex_t IT_0694 = IT_0435*IT_0693;
    const ccomplex_t IT_0695 = IT_0475*IT_0693;
    const ccomplex_t IT_0696 = IT_0479*IT_0693;
    const ccomplex_t IT_0697 = (-2)*IT_0677;
    const ccomplex_t IT_0698 = 2*IT_0263 + IT_0649 + IT_0650 + IT_0657 +
       IT_0658 + IT_0665 + IT_0666 + IT_0672 + IT_0673 + IT_0679 + IT_0680 +
       IT_0681 + IT_0682 + IT_0683 + IT_0684 + IT_0685 + IT_0686 + IT_0688 +
       IT_0689 + IT_0691 + IT_0692 + IT_0694 + IT_0695 + IT_0696 + (IT_0435 +
       IT_0475 + IT_0479)*IT_0697;
    const ccomplex_t IT_0699 = IT_0050*IT_0519;
    const ccomplex_t IT_0700 = -IT_0699;
    const ccomplex_t IT_0701 = IT_0603 + IT_0700;
    const ccomplex_t IT_0702 = (-8)*IT_0701;
    const ccomplex_t IT_0703 = IT_0608 + IT_0702;
    const ccomplex_t IT_0704 = -IT_0212;
    const ccomplex_t IT_0705 = -IT_0128;
    const ccomplex_t IT_0706 = (-4)*IT_0701;
    const ccomplex_t IT_0707 = (-2)*IT_0607;
    const ccomplex_t IT_0708 = IT_0706 + IT_0707;
    const ccomplex_t IT_0709 = conj(IT_0405)*IT_0656;
    const ccomplex_t IT_0710 = conj(IT_0428)*IT_0656;
    const ccomplex_t IT_0711 = conj(IT_0405)*IT_0664;
    const ccomplex_t IT_0712 = conj(IT_0428)*IT_0664;
    const ccomplex_t IT_0713 = conj(IT_0377)*IT_0671;
    const ccomplex_t IT_0714 = conj(IT_0217)*IT_0671;
    const ccomplex_t IT_0715 = conj(IT_0377)*IT_0678;
    const ccomplex_t IT_0716 = conj(IT_0217)*IT_0678;
    const ccomplex_t IT_0717 = IT_0086*conj(IT_0213);
    const ccomplex_t IT_0718 = IT_0086*conj(IT_0129);
    const ccomplex_t IT_0719 = conj(IT_0062)*IT_0612;
    const ccomplex_t IT_0720 = IT_0262*conj(IT_0439);
    const ccomplex_t IT_0721 = IT_0262*conj(IT_0286);
    const ccomplex_t IT_0722 = conj(IT_0062)*IT_0618;
    const ccomplex_t IT_0723 = conj(IT_0240)*IT_0687;
    const ccomplex_t IT_0724 = conj(IT_0313)*IT_0687;
    const ccomplex_t IT_0725 = conj(IT_0240)*IT_0690;
    const ccomplex_t IT_0726 = conj(IT_0313)*IT_0690;
    const ccomplex_t IT_0727 = conj(IT_0435)*IT_0693;
    const ccomplex_t IT_0728 = conj(IT_0475)*IT_0693;
    const ccomplex_t IT_0729 = conj(IT_0479)*IT_0693;
    const ccomplex_t IT_0730 = IT_0544 + IT_0545 + IT_0629 + IT_0630 + (conj
      (IT_0435) + conj(IT_0475) + conj(IT_0479))*IT_0697 + IT_0609*IT_0703 + 
      (conj(IT_0704) + conj(IT_0705))*IT_0708 + IT_0709 + IT_0710 + IT_0711 +
       IT_0712 + IT_0713 + IT_0714 + IT_0715 + IT_0716 + IT_0717 + IT_0718 +
       IT_0719 + IT_0720 + IT_0721 + IT_0722 + IT_0723 + IT_0724 + IT_0725 +
       IT_0726 + IT_0727 + IT_0728 + IT_0729;
    const ccomplex_t IT_0731 = -IT_0279 + IT_0438;
    const ccomplex_t IT_0732 = (-8)*IT_0611;
    const ccomplex_t IT_0733 = conj(IT_0027)*IT_0732;
    const ccomplex_t IT_0734 = (-4)*IT_0606;
    const ccomplex_t IT_0735 = conj(IT_0062)*IT_0734;
    const ccomplex_t IT_0736 = (-2)*IT_0588;
    const ccomplex_t IT_0737 = conj(IT_0313)*IT_0736;
    const ccomplex_t IT_0738 = conj(IT_0062)*IT_0707;
    const ccomplex_t IT_0739 = (-2)*IT_0617;
    const ccomplex_t IT_0740 = conj(IT_0027)*IT_0739;
    const ccomplex_t IT_0741 = -IT_0279 + -IT_0282 + IT_0285;
    const ccomplex_t IT_0742 = conj(IT_0112)*IT_0732;
    const ccomplex_t IT_0743 = conj(IT_0112)*IT_0739;
    const ccomplex_t IT_0744 = IT_0277 + IT_0544 + IT_0562 + IT_0575 + IT_0583
       + IT_0590 + IT_0593 + IT_0594 + IT_0620 + IT_0630 + IT_0639 + IT_0641 +
       IT_0646 + IT_0647 + IT_0733 + IT_0735 + conj(IT_0240)*IT_0736 + IT_0738 +
       IT_0740 + IT_0742 + IT_0743;
    const ccomplex_t IT_0745 = -IT_0282;
    const ccomplex_t IT_0746 = IT_0277 + IT_0545 + IT_0561 + IT_0574 + IT_0584
       + IT_0591 + IT_0592 + IT_0594 + IT_0619 + IT_0629 + IT_0638 + IT_0642 +
       IT_0643 + IT_0645 + IT_0733 + IT_0735 + IT_0737 + IT_0738 + IT_0740 +
       IT_0742 + IT_0743;
    const ccomplex_t IT_0747 = conj(IT_0731)*IT_0732;
    const ccomplex_t IT_0748 = conj(IT_0062)*IT_0732;
    const ccomplex_t IT_0749 = conj(IT_0027)*IT_0706;
    const ccomplex_t IT_0750 = conj(IT_0112)*IT_0706;
    const ccomplex_t IT_0751 = conj(IT_0027)*IT_0707;
    const ccomplex_t IT_0752 = conj(IT_0112)*IT_0707;
    const ccomplex_t IT_0753 = conj(IT_0062)*IT_0739;
    const ccomplex_t IT_0754 = 0.666666666666667*IT_0623 + (-0.166666666666667
      )*IT_0674;
    const ccomplex_t IT_0755 = s_24 + IT_0754;
    const ccomplex_t IT_0756 = m_N_1*IT_0755;
    const ccomplex_t IT_0757 = (-6)*IT_0756;
    const ccomplex_t IT_0758 = (-2)*IT_0757;
    const ccomplex_t IT_0759 = conj(IT_0435)*IT_0758;
    const ccomplex_t IT_0760 = IT_0114 + IT_0265 + IT_0545 + IT_0629 + IT_0710
       + IT_0711 + IT_0713 + IT_0716 + IT_0717 + IT_0720 + IT_0724 + IT_0725 +
       IT_0728 + IT_0729 + IT_0747 + IT_0748 + IT_0749 + IT_0750 + IT_0751 +
       IT_0752 + IT_0753 + IT_0759;
    const ccomplex_t IT_0761 = IT_0062*IT_0732;
    const ccomplex_t IT_0762 = IT_0112*IT_0706;
    const ccomplex_t IT_0763 = IT_0112*IT_0707;
    const ccomplex_t IT_0764 = IT_0731*IT_0739;
    const ccomplex_t IT_0765 = IT_0062*IT_0739;
    const ccomplex_t IT_0766 = IT_0475*IT_0758;
    const ccomplex_t IT_0767 = IT_0479*IT_0758;
    const ccomplex_t IT_0768 = IT_0087 + IT_0263 + IT_0649 + IT_0657 + IT_0666
       + IT_0673 + IT_0679 + IT_0682 + IT_0684 + IT_0686 + IT_0688 + IT_0692 +
       IT_0694 + IT_0705*IT_0706 + IT_0704*IT_0707 + IT_0732*IT_0741 + IT_0739
      *IT_0745 + IT_0761 + IT_0762 + IT_0763 + IT_0764 + IT_0765 + IT_0766 +
       IT_0767;
    const ccomplex_t IT_0769 = conj(IT_0731)*IT_0739;
    const ccomplex_t IT_0770 = conj(IT_0475)*IT_0758;
    const ccomplex_t IT_0771 = conj(IT_0479)*IT_0758;
    const ccomplex_t IT_0772 = IT_0114 + IT_0265 + IT_0544 + IT_0630 + IT_0709
       + IT_0712 + IT_0714 + IT_0715 + IT_0718 + IT_0721 + IT_0723 + IT_0726 +
       IT_0727 + IT_0748 + IT_0749 + IT_0750 + IT_0751 + IT_0752 + IT_0753 +
       IT_0769 + IT_0770 + IT_0771;
    const ccomplex_t IT_0773 = IT_0544 + IT_0545 + IT_0629 + IT_0630 + IT_0609
      *IT_0703 + IT_0709 + IT_0710 + IT_0711 + IT_0712 + IT_0713 + IT_0714 +
       IT_0715 + IT_0716 + IT_0717 + IT_0718 + IT_0719 + IT_0720 + IT_0721 +
       IT_0722 + IT_0723 + IT_0724 + IT_0725 + IT_0726 + IT_0727 + IT_0728 +
       IT_0729 + IT_0747 + IT_0759 + IT_0769 + IT_0770 + IT_0771;
    const ccomplex_t IT_0774 = s_12*IT_0355;
    const ccomplex_t IT_0775 = (-4)*IT_0356;
    const ccomplex_t IT_0776 = 0.25*conj(IT_0439);
    const ccomplex_t IT_0777 = 0.25*conj(IT_0286);
    const ccomplex_t IT_0778 = 0.25*conj(IT_0435);
    const ccomplex_t IT_0779 = conj(IT_0028) + conj(IT_0213);
    const ccomplex_t IT_0780 = 0.25*IT_0779;
    const ccomplex_t IT_0781 = 0.25*IT_0437;
    const ccomplex_t IT_0782 = 0.25*conj(IT_0240);
    const ccomplex_t IT_0783 = (-8)*IT_0515;
    const ccomplex_t IT_0784 = 0.25*conj(IT_0217);
    const ccomplex_t IT_0785 = 8*IT_0420;
    const ccomplex_t IT_0786 = 0.25*conj(IT_0405);
    const ccomplex_t IT_0787 = 0.25*IT_0357*(conj(IT_0313) + conj(IT_0377)) +
       0.25*conj(IT_0049)*(IT_0294 + IT_0465) + conj(IT_0351)*IT_0774 + 0.25*
      (conj(IT_0185) + conj(IT_0428) + conj(IT_0475) + conj(IT_0479))*IT_0775 +
       IT_0294*IT_0776 + IT_0465*IT_0777 + IT_0516*IT_0778 + IT_0196*IT_0780 +
       IT_0471*IT_0781 + IT_0421*IT_0782 + IT_0783*IT_0784 + IT_0785*IT_0786;
    const ccomplex_t IT_0788 = 0.25*conj(IT_0377);
    const ccomplex_t IT_0789 = 0.25*conj(IT_0313);
    const ccomplex_t IT_0790 = 0.25*conj(IT_0428);
    const ccomplex_t IT_0791 = 4*IT_0185*(0.25*(conj(IT_0217) + conj(IT_0240))
      *IT_0357 + conj(IT_0185)*IT_0774 + 0.25*(conj(IT_0351) + conj(IT_0405) +
       conj(IT_0435))*IT_0775 + IT_0471*IT_0780 + IT_0783*IT_0788 + IT_0421
      *IT_0789 + IT_0785*IT_0790);
    const ccomplex_t IT_0792 = IT_0024 + IT_0504;
    const ccomplex_t IT_0793 = IT_0024*IT_0792;
    const ccomplex_t IT_0794 = -IT_0793;
    const ccomplex_t IT_0795 = 6*IT_0794;
    const ccomplex_t IT_0796 = 0.25*conj(IT_0351);
    const ccomplex_t IT_0797 = IT_0465*IT_0796;
    const ccomplex_t IT_0798 = IT_0455*IT_0788;
    const ccomplex_t IT_0799 = 0.25*conj(IT_0213);
    const ccomplex_t IT_0800 = 0.25*conj(IT_0129);
    const ccomplex_t IT_0801 = IT_0446*IT_0782;
    const ccomplex_t IT_0802 = 2*IT_0445;
    const ccomplex_t IT_0803 = (-2)*IT_0332;
    const ccomplex_t IT_0804 = conj(IT_0028)*IT_0059;
    const ccomplex_t IT_0805 = 0.25*IT_0804;
    const ccomplex_t IT_0806 = conj(IT_0028)*IT_0208;
    const ccomplex_t IT_0807 = 0.25*IT_0806;
    const ccomplex_t IT_0808 = conj(IT_0185)*IT_0294;
    const ccomplex_t IT_0809 = 0.25*IT_0808;
    const ccomplex_t IT_0810 = conj(IT_0049)*IT_0300;
    const ccomplex_t IT_0811 = 0.25*IT_0810;
    const ccomplex_t IT_0812 = conj(IT_0217)*IT_0506 + IT_0300*IT_0777 +
       IT_0507*IT_0778 + IT_0333*IT_0789 + 0.25*(conj(IT_0049) + conj(IT_0439))
      *IT_0795 + IT_0797 + IT_0798 + IT_0208*IT_0799 + IT_0059*IT_0800 + IT_0801
       + IT_0786*IT_0802 + IT_0790*IT_0803 + IT_0805 + IT_0807 + IT_0809 +
       IT_0811;
    const ccomplex_t IT_0813 = conj(IT_0240)*IT_0333;
    const ccomplex_t IT_0814 = 0.25*IT_0813;
    const ccomplex_t IT_0815 = conj(IT_0185)*IT_0465;
    const ccomplex_t IT_0816 = 0.25*IT_0815;
    const ccomplex_t IT_0817 = IT_0294*conj(IT_0351);
    const ccomplex_t IT_0818 = 0.25*IT_0817;
    const ccomplex_t IT_0819 = conj(IT_0217)*IT_0455;
    const ccomplex_t IT_0820 = 0.25*IT_0819;
    const ccomplex_t IT_0821 = 0.25*IT_0300*(conj(IT_0286) + conj(IT_0439)) + 
      (conj(IT_0217) + conj(IT_0377))*IT_0506 + (IT_0333 + IT_0446)*IT_0789 +
       0.5*(conj(IT_0049) + 0.5*conj(IT_0286) + 0.5*conj(IT_0439))*IT_0795 +
       IT_0797 + IT_0798 + IT_0801 + 0.25*IT_0452*(IT_0802 + IT_0803) + IT_0807 
      + IT_0809 + IT_0814 + IT_0816 + IT_0818 + IT_0820;
    const ccomplex_t IT_0822 = 0.25*IT_0324;
    const ccomplex_t IT_0823 = conj(IT_0377)*IT_0506 + IT_0300*IT_0776 +
       IT_0322*IT_0778 + IT_0446*IT_0789 + IT_0059*IT_0799 + IT_0208*IT_0800 +
       IT_0790*IT_0802 + IT_0786*IT_0803 + IT_0805 + IT_0807 + IT_0811 + IT_0814
       + IT_0816 + IT_0818 + IT_0820 + IT_0795*IT_0822;
    const ccomplex_t IT_0824 = (-8)*IT_0522;
    const ccomplex_t IT_0825 = conj(IT_0028)*IT_0495;
    const ccomplex_t IT_0826 = 0.25*IT_0825;
    const ccomplex_t IT_0827 = (conj(IT_0049) + conj(IT_0439))*IT_0506 +
       IT_0407*IT_0778 + IT_0227*IT_0781 + IT_0393*IT_0782 + 0.25*conj(IT_0185)
      *IT_0783 + IT_0408*IT_0784 + IT_0436*IT_0786 + IT_0523*IT_0788 + IT_0386
      *IT_0789 + IT_0450*IT_0790 + IT_0357*IT_0796 + IT_0495*IT_0799 + IT_0455
      *IT_0822 + 0.25*IT_0480*IT_0824 + IT_0826;
    const ccomplex_t IT_0828 = 0.25*IT_0449 + IT_0324*IT_0506 + IT_0386
      *IT_0782 + IT_0523*IT_0784 + IT_0450*IT_0786 + IT_0408*IT_0788 + IT_0393
      *IT_0789 + IT_0436*IT_0790 + IT_0783*IT_0796 + IT_0495*IT_0800 + IT_0778
      *IT_0824 + IT_0826;
    const ccomplex_t IT_0829 = conj(IT_0049)*IT_0802;
    const ccomplex_t IT_0830 = conj(IT_0028)*IT_0455;
    const ccomplex_t IT_0831 = conj(IT_0049)*IT_0803;
    const ccomplex_t IT_0832 = conj(IT_0351)*IT_0775;
    const ccomplex_t IT_0833 = conj(IT_0240)*IT_0407 + conj(IT_0405)*IT_0408 +
       IT_0374*conj(IT_0428) + conj(IT_0313)*IT_0429 + IT_0393*conj(IT_0435) +
       conj(IT_0217)*IT_0436 + conj(IT_0377)*IT_0450 + IT_0437*IT_0453 + conj
      (IT_0213)*IT_0455 + IT_0386*IT_0480 + conj(IT_0185)*IT_0785 + conj(IT_0439
      )*IT_0802 + conj(IT_0286)*IT_0803 + IT_0829 + IT_0830 + IT_0831 + IT_0832;
    const ccomplex_t IT_0834 = conj(IT_0185)*IT_0775;
    const ccomplex_t IT_0835 = IT_0374*conj(IT_0405) + conj(IT_0313)*IT_0407 +
       IT_0408*conj(IT_0428) + conj(IT_0240)*IT_0429 + IT_0386*conj(IT_0435) +
       conj(IT_0377)*IT_0436 + conj(IT_0217)*IT_0450 + conj(IT_0129)*IT_0455 +
       IT_0453*IT_0779 + conj(IT_0351)*IT_0785 + conj(IT_0286)*IT_0802 + conj
      (IT_0439)*IT_0803 + IT_0829 + IT_0830 + IT_0831 + IT_0834;
    const ccomplex_t IT_0836 = conj(IT_0049)*IT_0322;
    const ccomplex_t IT_0837 = conj(IT_0028)*IT_0489;
    const ccomplex_t IT_0838 = conj(IT_0049)*IT_0507;
    const ccomplex_t IT_0839 = conj(IT_0028)*IT_0496;
    const ccomplex_t IT_0840 = IT_0386*conj(IT_0405) + conj(IT_0377)*IT_0407 +
       IT_0393*conj(IT_0428) + conj(IT_0313)*IT_0436 + IT_0322*conj(IT_0439) +
       conj(IT_0240)*IT_0450 + IT_0408*IT_0480 + conj(IT_0129)*IT_0489 + conj
      (IT_0213)*IT_0496 + conj(IT_0286)*IT_0507 + conj(IT_0351)*IT_0516 + conj
      (IT_0435)*IT_0523 + conj(IT_0217)*IT_0824 + IT_0834 + IT_0836 + IT_0837 +
       IT_0838 + IT_0839;
    const ccomplex_t IT_0841 = IT_0731*IT_0732;
    const ccomplex_t IT_0842 = IT_0435*IT_0758;
    const ccomplex_t IT_0843 = IT_0087 + IT_0263 + IT_0650 + IT_0658 + IT_0665
       + IT_0672 + IT_0680 + IT_0681 + IT_0683 + IT_0685 + IT_0689 + IT_0691 +
       IT_0695 + IT_0696 + IT_0704*IT_0706 + IT_0705*IT_0707 + IT_0739*IT_0741 +
       IT_0732*IT_0745 + IT_0761 + IT_0762 + IT_0763 + IT_0765 + IT_0841 +
       IT_0842;
    const ccomplex_t IT_0844 = IT_0649 + IT_0650 + IT_0657 + IT_0658 + IT_0665
       + IT_0666 + IT_0672 + IT_0673 + IT_0679 + IT_0680 + IT_0681 + IT_0682 +
       IT_0683 + IT_0684 + IT_0685 + IT_0686 + IT_0688 + IT_0689 + IT_0691 +
       IT_0692 + IT_0694 + IT_0695 + IT_0696 + IT_0764 + IT_0766 + IT_0767 +
       IT_0841 + IT_0842;
    const ccomplex_t IT_0845 = conj(IT_0185)*IT_0196;
    const ccomplex_t IT_0846 = conj(IT_0217)*IT_0495;
    const ccomplex_t IT_0847 = IT_0247*conj(IT_0313);
    const ccomplex_t IT_0848 = IT_0060 + IT_0203*conj(IT_0213) + IT_0059*conj
      (IT_0286) + conj(IT_0129)*IT_0304 + conj(IT_0240)*IT_0322 + IT_0227*conj
      (IT_0377) + IT_0208*conj(IT_0439) + conj(IT_0428)*IT_0453 + conj(IT_0405)
      *IT_0455 + conj(IT_0351)*IT_0471 + conj(IT_0435)*IT_0489 + IT_0491 +
       IT_0492 + IT_0494 + IT_0845 + IT_0846 + IT_0847;
    const ccomplex_t IT_0849 = IT_0405*IT_0560;
    const ccomplex_t IT_0850 = IT_0428*IT_0560;
    const ccomplex_t IT_0851 = IT_0377*IT_0573;
    const ccomplex_t IT_0852 = IT_0217*IT_0573;
    const ccomplex_t IT_0853 = IT_0217*IT_0582;
    const ccomplex_t IT_0854 = IT_0428*IT_0589;
    const ccomplex_t IT_0855 = IT_0075*IT_0213;
    const ccomplex_t IT_0856 = IT_0075*IT_0129;
    const ccomplex_t IT_0857 = IT_0049*IT_0276;
    const ccomplex_t IT_0858 = IT_0276*IT_0439;
    const ccomplex_t IT_0859 = IT_0276*IT_0286;
    const ccomplex_t IT_0860 = IT_0240*IT_0637;
    const ccomplex_t IT_0861 = IT_0313*IT_0637;
    const ccomplex_t IT_0862 = IT_0475*IT_0640;
    const ccomplex_t IT_0863 = IT_0479*IT_0640;
    const ccomplex_t IT_0864 = IT_0435*IT_0644;
    const ccomplex_t IT_0865 = IT_0475*IT_0644;
    const ccomplex_t IT_0866 = IT_0479*IT_0644;
    const ccomplex_t IT_0867 = IT_0377*IT_0582;
    const ccomplex_t IT_0868 = IT_0435*IT_0640;
    const ccomplex_t IT_0869 = IT_0405*IT_0589;
    const ccomplex_t IT_0870 = (IT_0240 + IT_0313)*IT_0636 + IT_0649 + IT_0650
       + IT_0685 + IT_0686 + IT_0849 + IT_0850 + IT_0851 + IT_0852 + IT_0853 +
       IT_0854 + IT_0855 + IT_0856 + 2*IT_0857 + IT_0858 + IT_0859 + IT_0860 +
       IT_0861 + IT_0862 + IT_0863 + IT_0864 + IT_0865 + IT_0866 + IT_0867 +
       IT_0868 + IT_0869;
    const ccomplex_t IT_0871 = IT_0027*IT_0732;
    const ccomplex_t IT_0872 = IT_0112*IT_0732;
    const ccomplex_t IT_0873 = IT_0062*IT_0734;
    const ccomplex_t IT_0874 = IT_0731 + IT_0745;
    const ccomplex_t IT_0875 = IT_0062*IT_0707;
    const ccomplex_t IT_0876 = IT_0027*IT_0739;
    const ccomplex_t IT_0877 = IT_0112*IT_0739;
    const ccomplex_t IT_0878 = IT_0076 + IT_0649 + IT_0686 + IT_0705*IT_0732 +
       IT_0240*IT_0736 + IT_0704*IT_0739 + IT_0734*IT_0741 + IT_0850 + IT_0852 +
       IT_0856 + IT_0857 + IT_0859 + IT_0861 + IT_0865 + IT_0866 + IT_0867 +
       IT_0868 + IT_0869 + IT_0871 + IT_0872 + IT_0873 + IT_0707*IT_0874 +
       IT_0875 + IT_0876 + IT_0877;
    const ccomplex_t IT_0879 = conj(IT_0286)*IT_0322;
    const ccomplex_t IT_0880 = conj(IT_0213)*IT_0489;
    const ccomplex_t IT_0881 = conj(IT_0185)*IT_0516;
    const ccomplex_t IT_0882 = IT_0393*conj(IT_0405);
    const ccomplex_t IT_0883 = conj(IT_0217)*IT_0407;
    const ccomplex_t IT_0884 = IT_0408*conj(IT_0435);
    const ccomplex_t IT_0885 = conj(IT_0313)*IT_0450;
    const ccomplex_t IT_0886 = conj(IT_0377)*IT_0824;
    const ccomplex_t IT_0887 = conj(IT_0240)*IT_0436;
    const ccomplex_t IT_0888 = conj(IT_0439)*IT_0507;
    const ccomplex_t IT_0889 = conj(IT_0129)*IT_0496;
    const ccomplex_t IT_0890 = IT_0386*conj(IT_0428);
    const ccomplex_t IT_0891 = IT_0313*IT_0736;
    const ccomplex_t IT_0892 = IT_0707*IT_0741;
    const ccomplex_t IT_0893 = IT_0076 + IT_0650 + IT_0685 + IT_0849 + IT_0851
       + IT_0853 + IT_0854 + IT_0855 + IT_0857 + IT_0858 + IT_0860 + IT_0862 +
       IT_0863 + IT_0864 + IT_0871 + IT_0873 + IT_0734*IT_0874 + IT_0875 +
       IT_0876 + IT_0891 + IT_0892;
    const ccomplex_t IT_0894 = IT_0076 + IT_0650 + IT_0685 + IT_0704*IT_0732 +
       IT_0705*IT_0739 + IT_0849 + IT_0851 + IT_0853 + IT_0854 + IT_0855 +
       IT_0857 + IT_0858 + IT_0860 + IT_0862 + IT_0863 + IT_0864 + IT_0871 +
       IT_0872 + IT_0873 + IT_0734*IT_0874 + IT_0875 + IT_0876 + IT_0877 +
       IT_0891 + IT_0892;
    const ccomplex_t IT_0895 = 2*IT_0061 + 2*IT_0077 + 2*IT_0088 + 2*IT_0113 +
       2*IT_0115 + 2*IT_0116 + conj(IT_0129)*(IT_0197 + IT_0198 + IT_0204 +
       IT_0209) + conj(IT_0213)*(IT_0198 + IT_0204 + IT_0209 + IT_0228 + IT_0248
      ) + conj(IT_0028)*(IT_0197 + 2*IT_0198 + IT_0204 + IT_0209 + IT_0228 +
       IT_0248) + 2*IT_0264 + 2*IT_0266 + 2*IT_0267 + 2*IT_0278 + conj(IT_0286)*
      (IT_0287 + IT_0288 + IT_0295) + 2*IT_0301 + 2*IT_0305 + IT_0313*IT_0448 +
       IT_0240*IT_0451 + IT_0452*(IT_0454 + IT_0456) + conj(IT_0439)*(IT_0287 +
       IT_0288 + IT_0457 + IT_0466 + IT_0467) + conj(IT_0049)*(IT_0288 + IT_0295
       + IT_0457 + IT_0466 + IT_0467) + IT_0213*IT_0497 + conj(IT_0435)*(IT_0498
       + IT_0499 + IT_0508 + IT_0509) + IT_0480*(IT_0498 + IT_0499 + IT_0508 +
       IT_0509 + IT_0510 + IT_0517 + IT_0524 + IT_0525 + IT_0526 + IT_0527 +
       IT_0528 + IT_0529 + IT_0530) + IT_0062*IT_0648 + conj(IT_0027)*IT_0698 +
       IT_0027*IT_0730 + IT_0731*(IT_0277 + IT_0545 + IT_0561 + IT_0574 +
       IT_0584 + IT_0591 + IT_0592 + IT_0594 + IT_0619 + IT_0629 + IT_0638 +
       IT_0642 + IT_0643 + IT_0645 + IT_0733 + IT_0735 + IT_0737 + IT_0738 +
       IT_0740) + IT_0741*IT_0744 + IT_0745*IT_0746 + IT_0704*IT_0760 + conj
      (IT_0705)*IT_0768 + IT_0705*IT_0772 + IT_0112*IT_0773 + 4*IT_0351*IT_0787 
      + IT_0791 + 4*IT_0286*IT_0812 + 4*IT_0049*IT_0821 + 4*IT_0439*IT_0823 + 4
      *IT_0377*IT_0827 + 4*IT_0217*IT_0828 + IT_0428*IT_0833 + IT_0405*IT_0835 +
       IT_0435*IT_0840 + conj(IT_0704)*IT_0843 + conj(IT_0112)*IT_0844 + IT_0028
      *((conj(IT_0129) + conj(IT_0213))*IT_0304 + conj(IT_0351)*IT_0471 +
       IT_0472 + IT_0473 + IT_0490 + IT_0491 + IT_0492 + IT_0493 + conj(IT_0377)
      *(IT_0227 + IT_0495) + IT_0322*IT_0631 + IT_0845 + IT_0846 + IT_0847) +
       IT_0129*IT_0848 + conj(IT_0062)*IT_0870 + conj(IT_0741)*IT_0878 + IT_0475
      *(IT_0832 + IT_0836 + IT_0837 + IT_0838 + IT_0839 + IT_0879 + IT_0880 +
       IT_0881 + IT_0882 + IT_0883 + IT_0884 + IT_0885 + IT_0886 + IT_0887 +
       IT_0888 + IT_0889 + IT_0890) + IT_0479*(IT_0480*IT_0523 + IT_0832 +
       IT_0836 + IT_0837 + IT_0838 + IT_0839 + IT_0879 + IT_0880 + IT_0881 +
       IT_0882 + IT_0883 + IT_0884 + IT_0885 + IT_0886 + IT_0887 + IT_0888 +
       IT_0889 + IT_0890) + conj(IT_0731)*IT_0893 + conj(IT_0745)*IT_0894;
    return create_ccomplex_return(IT_0895);
}

