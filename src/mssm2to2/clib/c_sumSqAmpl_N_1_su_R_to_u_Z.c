#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_su_R_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_su_R_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
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
    const ccomplex_t IT_0004 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = pow(m_u, 2);
    const ccomplex_t IT_0017 = pow(m_su_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = IT_0007*IT_0022;
    const ccomplex_t IT_0025 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*IT_0025;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0023 + IT_0024 + 
      -IT_0026 + -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0031 = 0.666666666666667*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0029*IT_0034;
    const ccomplex_t IT_0036 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = IT_0007*IT_0036;
    const ccomplex_t IT_0039 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0039;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0037 + IT_0038 + 
      -IT_0040 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0045 = 0.666666666666667*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0043*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0051 = 1.33333333333333*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0053 = 0.666666666666667*IT_0052;
    const ccomplex_t IT_0054 = pow(m_N_1, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0016 + IT_0054 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0051*IT_0057;
    const ccomplex_t IT_0059 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = IT_0007*IT_0059;
    const ccomplex_t IT_0062 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = IT_0007*IT_0062;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + IT_0061 + 
      -IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0053*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0066*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0065;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = 2*IT_0021 + 2*IT_0035 + 2*IT_0049 + IT_0058 +
       IT_0070 + -IT_0072;
    const ccomplex_t IT_0074 = pow(m_Z, -2);
    const ccomplex_t IT_0075 = pow(s_14, 2);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = -IT_0054;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = s_13*IT_0078;
    const ccomplex_t IT_0080 = pow(m_W, -1);
    const ccomplex_t IT_0081 = sin(beta);
    const ccomplex_t IT_0082 = cpow(IT_0081, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0006*IT_0080*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0055*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0051*IT_0086;
    const ccomplex_t IT_0088 = IT_0067*IT_0084;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0071*IT_0089;
    const ccomplex_t IT_0091 = IT_0066*IT_0089;
    const ccomplex_t IT_0092 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0093 = IT_0003*IT_0092;
    const ccomplex_t IT_0094 = IT_0007*IT_0092;
    const ccomplex_t IT_0095 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = IT_0007*IT_0095;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0093 + IT_0094 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0006*IT_0080*IT_0082;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0018*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0099*IT_0103;
    const ccomplex_t IT_0105 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0106 = IT_0003*IT_0105;
    const ccomplex_t IT_0107 = IT_0007*IT_0105;
    const ccomplex_t IT_0108 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0109 = IT_0003*IT_0108;
    const ccomplex_t IT_0110 = IT_0007*IT_0108;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0106 + IT_0107 + 
      -IT_0109 + -IT_0110);
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0006*IT_0080*IT_0082;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0032*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0112*IT_0116;
    const ccomplex_t IT_0118 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0119 = IT_0003*IT_0118;
    const ccomplex_t IT_0120 = IT_0007*IT_0118;
    const ccomplex_t IT_0121 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0122 = IT_0003*IT_0121;
    const ccomplex_t IT_0123 = IT_0007*IT_0121;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0119 + IT_0120 + 
      -IT_0122 + -IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0006*IT_0080*IT_0082;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0046*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0125*IT_0129;
    const ccomplex_t IT_0131 = -IT_0087 + -IT_0090 + IT_0091 + (-2)*IT_0104 + 
      (-2)*IT_0117 + (-2)*IT_0130;
    const ccomplex_t IT_0132 = m_u*m_N_1;
    const ccomplex_t IT_0133 = IT_0078*IT_0132;
    const ccomplex_t IT_0134 = -IT_0021 + -IT_0035 + -IT_0049 + (-0.5)*IT_0070
       + 0.5*IT_0072;
    const ccomplex_t IT_0135 = s_13*s_14;
    const ccomplex_t IT_0136 = s_34*IT_0054;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = m_u*IT_0078;
    const ccomplex_t IT_0140 = m_N_1*IT_0071;
    const ccomplex_t IT_0141 = IT_0089*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = e_em*IT_0003;
    const ccomplex_t IT_0144 = e_em*IT_0007;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0143 + (-3)*IT_0144);
    const ccomplex_t IT_0146 = (-0.166666666666667)*IT_0145;
    const ccomplex_t IT_0147 = IT_0053*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0054 
      + -reg_prop, -1);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0150 = IT_0147*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = m_N_1*IT_0066;
    const ccomplex_t IT_0152 = IT_0089*IT_0151;
    const ccomplex_t IT_0153 = m_N_2*IT_0013;
    const ccomplex_t IT_0154 = m_N_1*IT_0099;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = IT_0103*IT_0156;
    const ccomplex_t IT_0158 = -IT_0140;
    const ccomplex_t IT_0159 = IT_0151 + IT_0158;
    const ccomplex_t IT_0160 = IT_0089*IT_0159;
    const ccomplex_t IT_0161 = m_N_3*IT_0029;
    const ccomplex_t IT_0162 = m_N_1*IT_0112;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = IT_0161 + IT_0163;
    const ccomplex_t IT_0165 = IT_0116*IT_0164;
    const ccomplex_t IT_0166 = m_N_4*IT_0043;
    const ccomplex_t IT_0167 = m_N_1*IT_0125;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = IT_0129*IT_0169;
    const ccomplex_t IT_0171 = (-0.666666666666667)*IT_0050;
    const ccomplex_t IT_0172 = IT_0053*IT_0171;
    const ccomplex_t IT_0173 = IT_0148*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = m_u*IT_0174;
    const ccomplex_t IT_0176 = -IT_0150 + (-0.5)*IT_0152 + -IT_0157 + (-0.5)
      *IT_0160 + -IT_0165 + -IT_0170 + IT_0175;
    const ccomplex_t IT_0177 = IT_0142 + IT_0176;
    const ccomplex_t IT_0178 = m_N_4*IT_0125;
    const ccomplex_t IT_0179 = m_N_1*IT_0043;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = IT_0048*IT_0181;
    const ccomplex_t IT_0183 = IT_0084*IT_0171;
    const ccomplex_t IT_0184 = IT_0148*IT_0149*IT_0183;
    const ccomplex_t IT_0185 = IT_0084*IT_0146;
    const ccomplex_t IT_0186 = IT_0148*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = m_u*IT_0187;
    const ccomplex_t IT_0189 = IT_0069*IT_0140;
    const ccomplex_t IT_0190 = IT_0069*IT_0159;
    const ccomplex_t IT_0191 = m_N_2*IT_0099;
    const ccomplex_t IT_0192 = m_N_1*IT_0013;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0020*IT_0194;
    const ccomplex_t IT_0196 = m_N_3*IT_0112;
    const ccomplex_t IT_0197 = m_N_1*IT_0029;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0034*IT_0199;
    const ccomplex_t IT_0201 = IT_0182 + IT_0184 + -IT_0188 + 0.5*IT_0189 + (
      -0.5)*IT_0190 + IT_0195 + IT_0200;
    const ccomplex_t IT_0202 = IT_0069*IT_0151;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = (-18)*conj(IT_0176);
    const ccomplex_t IT_0206 = -IT_0058 + (-2)*IT_0174;
    const ccomplex_t IT_0207 = s_14*s_34*IT_0074;
    const ccomplex_t IT_0208 = -IT_0207;
    const ccomplex_t IT_0209 = s_13 + IT_0208;
    const ccomplex_t IT_0210 = s_13*IT_0209;
    const ccomplex_t IT_0211 = IT_0087 + 2*IT_0187;
    const ccomplex_t IT_0212 = IT_0132*IT_0209;
    const ccomplex_t IT_0213 = -IT_0174;
    const ccomplex_t IT_0214 = m_N_1*IT_0209;
    const ccomplex_t IT_0215 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0216 = (-6)*conj(IT_0131);
    const ccomplex_t IT_0217 = pow(m_Z, 2);
    const ccomplex_t IT_0218 = s_14*IT_0217;
    const ccomplex_t IT_0219 = IT_0074*IT_0218;
    const ccomplex_t IT_0220 = (-0.25)*IT_0219;
    const ccomplex_t IT_0221 = s_14 + IT_0220;
    const ccomplex_t IT_0222 = m_u*IT_0221;
    const ccomplex_t IT_0223 = (-4)*IT_0222;
    const ccomplex_t IT_0224 = 0.5*IT_0219;
    const ccomplex_t IT_0225 = s_14 + IT_0224;
    const ccomplex_t IT_0226 = m_u*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0226;
    const ccomplex_t IT_0228 = 0.5*IT_0090 + (-0.5)*IT_0091 + IT_0104 +
       IT_0117 + IT_0130;
    const ccomplex_t IT_0229 = s_34*IT_0217;
    const ccomplex_t IT_0230 = IT_0074*IT_0229;
    const ccomplex_t IT_0231 = (-0.25)*IT_0230;
    const ccomplex_t IT_0232 = s_34 + IT_0231;
    const ccomplex_t IT_0233 = m_N_1*IT_0232;
    const ccomplex_t IT_0234 = (-4)*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0230;
    const ccomplex_t IT_0236 = s_34 + IT_0235;
    const ccomplex_t IT_0237 = m_N_1*IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0237;
    const ccomplex_t IT_0239 = s_13*IT_0217;
    const ccomplex_t IT_0240 = s_14*s_34;
    const ccomplex_t IT_0241 = (-2)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = IT_0074*IT_0242;
    const ccomplex_t IT_0244 = (-6)*IT_0243;
    const ccomplex_t IT_0245 = 12*s_13;
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = m_u*IT_0209;
    const ccomplex_t IT_0248 = (-6)*IT_0206;
    const ccomplex_t IT_0249 = pow(s_34, 2);
    const ccomplex_t IT_0250 = IT_0074*IT_0249;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0016 + IT_0251;
    const ccomplex_t IT_0253 = IT_0132*IT_0252;
    const ccomplex_t IT_0254 = s_14*IT_0016;
    const ccomplex_t IT_0255 = s_13*s_34;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = IT_0254 + IT_0256;
    const ccomplex_t IT_0258 = m_N_1*IT_0252;
    const ccomplex_t IT_0259 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0260 = (-6)*conj(IT_0211);
    const ccomplex_t IT_0261 = (-6)*conj(IT_0176);
    const ccomplex_t IT_0262 = (-6)*IT_0073;
    const ccomplex_t IT_0263 = (-6)*IT_0211;
    const ccomplex_t IT_0264 = 6*IT_0213;
    const ccomplex_t IT_0265 = 6*IT_0134;
    const ccomplex_t IT_0266 = 6*IT_0228;
    const ccomplex_t IT_0267 = 6*IT_0187;
    const ccomplex_t IT_0268 = 6*IT_0131*IT_0139 + IT_0132*((-18)*IT_0142 + (
      -18)*IT_0176) + IT_0204*IT_0246 + IT_0248*IT_0258 + IT_0214*IT_0262 +
       IT_0247*IT_0263 + IT_0238*IT_0264 + IT_0234*IT_0265 + IT_0227*IT_0266 +
       IT_0223*IT_0267;
    const ccomplex_t IT_0269 = (-6)*IT_0131;
    const ccomplex_t IT_0270 = 6*IT_0073*IT_0139 + IT_0132*((-18)*IT_0201 + (
      -18)*IT_0203) + IT_0177*IT_0246 + IT_0247*IT_0248 + IT_0258*IT_0263 +
       IT_0223*IT_0264 + IT_0227*IT_0265 + IT_0234*IT_0266 + IT_0238*IT_0267 +
       IT_0214*IT_0269;
    const ccomplex_t IT_0271 = IT_0132*IT_0217;
    const ccomplex_t IT_0272 = 3*IT_0271;
    const ccomplex_t IT_0273 = (-3)*IT_0271;
    const ccomplex_t IT_0274 = pow(m_Z, 4);
    const ccomplex_t IT_0275 = s_13*IT_0274;
    const ccomplex_t IT_0276 = IT_0074*IT_0275;
    const ccomplex_t IT_0277 = (-24)*IT_0240;
    const ccomplex_t IT_0278 = 6*IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = 0.166666666666667*IT_0213;
    const ccomplex_t IT_0280 = 24*IT_0240;
    const ccomplex_t IT_0281 = (-12)*IT_0239;
    const ccomplex_t IT_0282 = 6*IT_0276 + IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 0.166666666666667*IT_0134;
    const ccomplex_t IT_0284 = 0.166666666666667*IT_0228;
    const ccomplex_t IT_0285 = 0.166666666666667*IT_0187;
    const ccomplex_t IT_0286 = (-6)*IT_0176;
    const ccomplex_t IT_0287 = 0.166666666666667*IT_0212;
    const ccomplex_t IT_0288 = 0.166666666666667*IT_0210;
    const ccomplex_t IT_0289 = s_13*IT_0252;
    const ccomplex_t IT_0290 = 0.166666666666667*IT_0289;
    const ccomplex_t IT_0291 = 0.166666666666667*IT_0138;
    const ccomplex_t IT_0292 = 0.166666666666667*IT_0257;
    const ccomplex_t IT_0293 = 6*conj(IT_0073)*(IT_0073*IT_0079 + IT_0131
      *IT_0133 + IT_0134*IT_0138 + IT_0139*IT_0177) + IT_0132*IT_0204*IT_0205 + 
      (IT_0206*IT_0210 + IT_0211*IT_0212 + IT_0138*IT_0213 + IT_0204*IT_0214)
      *IT_0215 + (IT_0138*IT_0187 + IT_0210*IT_0211 + IT_0142*IT_0214)*IT_0216 +
       6*conj(IT_0176)*(IT_0073*IT_0139 + IT_0213*IT_0223 + IT_0134*IT_0227 +
       IT_0228*IT_0234 + IT_0187*IT_0238 + 0.166666666666667*IT_0177*IT_0246 +
       0.166666666666667*IT_0247*IT_0248) + (IT_0131*IT_0212 + IT_0177*IT_0247 +
       IT_0211*IT_0253 + IT_0213*IT_0257 + IT_0204*IT_0258)*IT_0259 + (IT_0204
      *IT_0247 + IT_0187*IT_0257 + IT_0142*IT_0258)*IT_0260 + (IT_0131*IT_0214 +
       IT_0211*IT_0258)*IT_0261 + (conj(IT_0201) + conj(IT_0203))*IT_0268 + conj
      (IT_0142)*IT_0270 + 6*conj(IT_0134)*(IT_0073*IT_0138 + IT_0177*IT_0227 +
       IT_0204*IT_0234 + IT_0206*IT_0257 + IT_0187*IT_0272 + IT_0228*IT_0273 +
       IT_0278*IT_0279 + IT_0282*IT_0283) + 6*conj(IT_0228)*(IT_0131*IT_0138 +
       IT_0204*IT_0227 + IT_0177*IT_0234 + IT_0211*IT_0257 + IT_0213*IT_0272 +
       IT_0134*IT_0273 + IT_0282*IT_0284 + IT_0278*IT_0285) + 6*conj(IT_0131)*
      (IT_0079*IT_0131 + IT_0073*IT_0133 + IT_0139*IT_0204 + IT_0138*IT_0228 +
       0.166666666666667*IT_0214*IT_0286 + IT_0248*IT_0287) + 6*conj(IT_0206)*
      (IT_0134*IT_0257 + IT_0262*IT_0288 + IT_0248*IT_0290) + 6*conj(IT_0211)*
      (0.166666666666667*IT_0248*IT_0253 + IT_0228*IT_0257 + 0.166666666666667
      *IT_0258*IT_0286 + IT_0262*IT_0287 + IT_0269*IT_0288 + IT_0263*IT_0290) +
       6*conj(IT_0213)*(IT_0177*IT_0223 + IT_0204*IT_0238 + IT_0228*IT_0272 +
       IT_0187*IT_0273 + IT_0279*IT_0282 + IT_0278*IT_0283 + IT_0262*IT_0291 +
       IT_0248*IT_0292) + 6*conj(IT_0187)*(IT_0204*IT_0223 + IT_0177*IT_0238 +
       IT_0134*IT_0272 + IT_0213*IT_0273 + IT_0278*IT_0284 + IT_0282*IT_0285 +
       IT_0269*IT_0291 + IT_0263*IT_0292);
    return create_ccomplex_return(IT_0293);
}

