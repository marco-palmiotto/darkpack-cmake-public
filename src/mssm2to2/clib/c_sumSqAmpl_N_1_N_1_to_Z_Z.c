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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = s_12*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_Z, -4);
    const ccomplex_t IT_0006 = s_34*IT_0005;
    const ccomplex_t IT_0007 = s_14*s_23;
    const ccomplex_t IT_0008 = s_12*s_34;
    const ccomplex_t IT_0009 = IT_0007 + IT_0008;
    const ccomplex_t IT_0010 = s_13*s_24;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0009 + IT_0011;
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = 0.25*IT_0013;
    const ccomplex_t IT_0015 = s_12 + IT_0004 + IT_0014;
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0020, -1);
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0027;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + IT_0026 + 
      -IT_0028 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0030;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = m_N_1*IT_0033;
    const ccomplex_t IT_0035 = cpow(IT_0031, 2);
    const ccomplex_t IT_0036 = m_N_1*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0034 + IT_0037;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = pow(m_N_1, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0000 + IT_0040 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0044 = IT_0021*IT_0043;
    const ccomplex_t IT_0045 = IT_0025*IT_0043;
    const ccomplex_t IT_0046 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0047 = IT_0021*IT_0046;
    const ccomplex_t IT_0048 = IT_0025*IT_0046;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0044 + IT_0045 + 
      -IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = m_N_3*IT_0051;
    const ccomplex_t IT_0053 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0054 = IT_0021*IT_0053;
    const ccomplex_t IT_0055 = IT_0025*IT_0053;
    const ccomplex_t IT_0056 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0057 = IT_0021*IT_0056;
    const ccomplex_t IT_0058 = IT_0025*IT_0056;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0054 + IT_0055 + 
      -IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0050*IT_0060;
    const ccomplex_t IT_0062 = m_N_1*IT_0061;
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = IT_0052 + IT_0063;
    const ccomplex_t IT_0065 = IT_0042*IT_0064;
    const ccomplex_t IT_0066 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0067 = IT_0021*IT_0066;
    const ccomplex_t IT_0068 = IT_0025*IT_0066;
    const ccomplex_t IT_0069 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0070 = IT_0021*IT_0069;
    const ccomplex_t IT_0071 = IT_0025*IT_0069;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0067 + IT_0068 + 
      -IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = cpow(IT_0073, 2);
    const ccomplex_t IT_0075 = m_N_2*IT_0074;
    const ccomplex_t IT_0076 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0077 = IT_0021*IT_0076;
    const ccomplex_t IT_0078 = IT_0025*IT_0076;
    const ccomplex_t IT_0079 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0080 = IT_0021*IT_0079;
    const ccomplex_t IT_0081 = IT_0025*IT_0079;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0077 + IT_0078 + 
      -IT_0080 + -IT_0081);
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0073*IT_0083;
    const ccomplex_t IT_0085 = m_N_1*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0075 + IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_23 + IT_0000 + IT_0040 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0092 = IT_0021*IT_0091;
    const ccomplex_t IT_0093 = IT_0025*IT_0091;
    const ccomplex_t IT_0094 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0095 = IT_0021*IT_0094;
    const ccomplex_t IT_0096 = IT_0025*IT_0094;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0092 + IT_0093 + 
      -IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = cpow(IT_0098, 2);
    const ccomplex_t IT_0100 = m_N_4*IT_0099;
    const ccomplex_t IT_0101 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0102 = IT_0021*IT_0101;
    const ccomplex_t IT_0103 = IT_0025*IT_0101;
    const ccomplex_t IT_0104 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0105 = IT_0021*IT_0104;
    const ccomplex_t IT_0106 = IT_0025*IT_0104;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0102 + IT_0103 + 
      -IT_0105 + -IT_0106);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0098*IT_0108;
    const ccomplex_t IT_0110 = m_N_1*IT_0109;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0100 + IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0000 + IT_0040 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = IT_0039 + -IT_0065 + -IT_0090 + -IT_0115;
    const ccomplex_t IT_0117 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0118 = cos(alpha);
    const ccomplex_t IT_0119 = IT_0019*IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0122 = sin(alpha);
    const ccomplex_t IT_0123 = IT_0019*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0126 = IT_0024*IT_0118;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0129 = IT_0024*IT_0122;
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0120 + -IT_0124 + 
      -IT_0127 + IT_0130);
    const ccomplex_t IT_0132 = cos(beta);
    const ccomplex_t IT_0133 = cpow(IT_0018, -4);
    const ccomplex_t IT_0134 = IT_0020*IT_0118*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = sin(beta);
    const ccomplex_t IT_0136 = IT_0020*IT_0122*IT_0133*IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0020, -3);
    const ccomplex_t IT_0138 = IT_0118*IT_0132*IT_0137;
    const ccomplex_t IT_0139 = IT_0122*IT_0135*IT_0137;
    const ccomplex_t IT_0140 = cpow(IT_0018, -2);
    const ccomplex_t IT_0141 = IT_0024*IT_0118*IT_0132*IT_0140;
    const ccomplex_t IT_0142 = IT_0024*IT_0122*IT_0135*IT_0140;
    const ccomplex_t IT_0143 = cpow(IT_0020, -2);
    const ccomplex_t IT_0144 = cpow(IT_0140 + IT_0143, -1);
    const ccomplex_t IT_0145 = (IT_0134 + IT_0136 + IT_0138 + IT_0139 + 2
      *IT_0141 + 2*IT_0142)*IT_0144;
    const ccomplex_t IT_0146 = m_W*e_em;
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_12 + (-2)*IT_0040 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0131*IT_0151;
    const ccomplex_t IT_0153 = IT_0119*IT_0121;
    const ccomplex_t IT_0154 = IT_0117*IT_0123;
    const ccomplex_t IT_0155 = IT_0126*IT_0128;
    const ccomplex_t IT_0156 = IT_0125*IT_0129;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + IT_0154 + 
      -IT_0155 + -IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0020*IT_0118*IT_0133*IT_0135;
    const ccomplex_t IT_0160 = IT_0020*IT_0122*IT_0132*IT_0133;
    const ccomplex_t IT_0161 = IT_0118*IT_0135*IT_0137;
    const ccomplex_t IT_0162 = IT_0122*IT_0132*IT_0137;
    const ccomplex_t IT_0163 = IT_0024*IT_0118*IT_0135*IT_0140;
    const ccomplex_t IT_0164 = IT_0024*IT_0122*IT_0132*IT_0140;
    const ccomplex_t IT_0165 = IT_0144*(IT_0159 + -IT_0160 + IT_0161 + 
      -IT_0162 + 2*IT_0163 + (-2)*IT_0164);
    const ccomplex_t IT_0166 = IT_0146*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0040 + -reg_prop, -1);
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0158*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0152 + -IT_0171;
    const ccomplex_t IT_0173 = 0.5*IT_0152;
    const ccomplex_t IT_0174 = conj(IT_0172) + -conj(IT_0173);
    const ccomplex_t IT_0175 = cpow(IT_0032, 2);
    const ccomplex_t IT_0176 = m_N_1*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0034 + IT_0177;
    const ccomplex_t IT_0179 = IT_0017*IT_0178;
    const ccomplex_t IT_0180 = cpow(IT_0083, 2);
    const ccomplex_t IT_0181 = m_N_2*IT_0180;
    const ccomplex_t IT_0182 = IT_0086 + IT_0181;
    const ccomplex_t IT_0183 = IT_0089*IT_0182;
    const ccomplex_t IT_0184 = cpow(IT_0060, 2);
    const ccomplex_t IT_0185 = m_N_3*IT_0184;
    const ccomplex_t IT_0186 = IT_0063 + IT_0185;
    const ccomplex_t IT_0187 = IT_0042*IT_0186;
    const ccomplex_t IT_0188 = cpow(IT_0108, 2);
    const ccomplex_t IT_0189 = m_N_4*IT_0188;
    const ccomplex_t IT_0190 = IT_0111 + IT_0189;
    const ccomplex_t IT_0191 = IT_0114*IT_0190;
    const ccomplex_t IT_0192 = IT_0179 + -IT_0183 + -IT_0187 + -IT_0191;
    const ccomplex_t IT_0193 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0194 = IT_0119*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0196 = IT_0123*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0198 = IT_0126*IT_0197;
    const ccomplex_t IT_0199 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0200 = IT_0129*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0194 + -IT_0196 + 
      -IT_0198 + IT_0200);
    const ccomplex_t IT_0202 = IT_0151*IT_0201;
    const ccomplex_t IT_0203 = IT_0119*IT_0195;
    const ccomplex_t IT_0204 = IT_0123*IT_0193;
    const ccomplex_t IT_0205 = IT_0126*IT_0199;
    const ccomplex_t IT_0206 = IT_0129*IT_0197;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0203 + IT_0204 + 
      -IT_0205 + -IT_0206);
    const ccomplex_t IT_0208 = -IT_0207;
    const ccomplex_t IT_0209 = IT_0170*IT_0208;
    const ccomplex_t IT_0210 = (-0.5)*IT_0202 + -IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0202;
    const ccomplex_t IT_0212 = conj(IT_0210) + -conj(IT_0211);
    const ccomplex_t IT_0213 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = IT_0033*IT_0214;
    const ccomplex_t IT_0216 = IT_0035*IT_0214;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0000 + IT_0040 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0084*IT_0218;
    const ccomplex_t IT_0220 = cpow((-2)*s_13 + IT_0000 + IT_0040 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = IT_0061*IT_0221;
    const ccomplex_t IT_0223 = cpow((-2)*s_13 + IT_0000 + IT_0040 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0109*IT_0224;
    const ccomplex_t IT_0226 = -IT_0215 + IT_0216 + (-2)*IT_0219 + (-2)
      *IT_0222 + (-2)*IT_0225;
    const ccomplex_t IT_0227 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0228 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0229 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0230 = -IT_0228 + -IT_0229;
    const ccomplex_t IT_0231 = s_12 + IT_0227 + IT_0230;
    const ccomplex_t IT_0232 = m_N_1*IT_0231;
    const ccomplex_t IT_0233 = 2*conj(IT_0172);
    const ccomplex_t IT_0234 = (-2)*conj(IT_0173);
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0237 = IT_0040 + IT_0236;
    const ccomplex_t IT_0238 = pow(s_23, 2);
    const ccomplex_t IT_0239 = IT_0002*IT_0238;
    const ccomplex_t IT_0240 = pow(s_24, 2);
    const ccomplex_t IT_0241 = IT_0002*IT_0240;
    const ccomplex_t IT_0242 = -IT_0239 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0237 + IT_0242;
    const ccomplex_t IT_0244 = m_N_1*IT_0243;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = 2*conj(IT_0210);
    const ccomplex_t IT_0247 = (-2)*conj(IT_0211);
    const ccomplex_t IT_0248 = IT_0246 + IT_0247;
    const ccomplex_t IT_0249 = s_13*IT_0000;
    const ccomplex_t IT_0250 = IT_0006*IT_0249;
    const ccomplex_t IT_0251 = s_14*IT_0000;
    const ccomplex_t IT_0252 = IT_0002*IT_0251;
    const ccomplex_t IT_0253 = s_13*s_34;
    const ccomplex_t IT_0254 = (-2)*IT_0253;
    const ccomplex_t IT_0255 = IT_0251 + IT_0254;
    const ccomplex_t IT_0256 = IT_0002*IT_0255;
    const ccomplex_t IT_0257 = (-0.5)*IT_0250 + 0.5*IT_0252 + (-0.5)*IT_0256;
    const ccomplex_t IT_0258 = s_14 + IT_0257;
    const ccomplex_t IT_0259 = m_N_1*IT_0258;
    const ccomplex_t IT_0260 = (-2)*IT_0259;
    const ccomplex_t IT_0261 = IT_0017*IT_0033;
    const ccomplex_t IT_0262 = IT_0109*IT_0114;
    const ccomplex_t IT_0263 = IT_0042*IT_0061;
    const ccomplex_t IT_0264 = IT_0017*IT_0175;
    const ccomplex_t IT_0265 = IT_0084*IT_0089;
    const ccomplex_t IT_0266 = (-0.5)*IT_0261 + -IT_0262 + -IT_0263 + 0.5
      *IT_0264 + -IT_0265;
    const ccomplex_t IT_0267 = 0.5*IT_0215 + (-0.5)*IT_0216 + IT_0219 +
       IT_0222 + IT_0225;
    const ccomplex_t IT_0268 = conj(IT_0266) + conj(IT_0267);
    const ccomplex_t IT_0269 = 2*IT_0172;
    const ccomplex_t IT_0270 = (-2)*IT_0173;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = IT_0266 + IT_0267;
    const ccomplex_t IT_0273 = IT_0017*IT_0035;
    const ccomplex_t IT_0274 = (-0.5)*IT_0261 + -IT_0262 + -IT_0263 + -IT_0265
       + 0.5*IT_0273;
    const ccomplex_t IT_0275 = IT_0175*IT_0214;
    const ccomplex_t IT_0276 = 0.5*IT_0215 + IT_0219 + IT_0222 + IT_0225 + (
      -0.5)*IT_0275;
    const ccomplex_t IT_0277 = conj(IT_0274) + conj(IT_0276);
    const ccomplex_t IT_0278 = 2*IT_0210;
    const ccomplex_t IT_0279 = (-2)*IT_0211;
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = IT_0274 + IT_0276;
    const ccomplex_t IT_0282 = IT_0261 + 2*IT_0262 + 2*IT_0263 + 2*IT_0265 + 
      -IT_0273;
    const ccomplex_t IT_0283 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0284 = s_12 + IT_0230 + IT_0283;
    const ccomplex_t IT_0285 = m_N_1*IT_0284;
    const ccomplex_t IT_0286 = 2*conj(IT_0173);
    const ccomplex_t IT_0287 = (-2)*conj(IT_0172);
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0290 = IT_0040 + IT_0289;
    const ccomplex_t IT_0291 = pow(s_13, 2);
    const ccomplex_t IT_0292 = IT_0002*IT_0291;
    const ccomplex_t IT_0293 = pow(s_14, 2);
    const ccomplex_t IT_0294 = IT_0002*IT_0293;
    const ccomplex_t IT_0295 = -IT_0292 + -IT_0294;
    const ccomplex_t IT_0296 = IT_0290 + IT_0295;
    const ccomplex_t IT_0297 = m_N_1*IT_0296;
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = 2*conj(IT_0211);
    const ccomplex_t IT_0300 = (-2)*conj(IT_0210);
    const ccomplex_t IT_0301 = IT_0299 + IT_0300;
    const ccomplex_t IT_0302 = s_23*IT_0000;
    const ccomplex_t IT_0303 = IT_0006*IT_0302;
    const ccomplex_t IT_0304 = s_24*IT_0000;
    const ccomplex_t IT_0305 = IT_0002*IT_0304;
    const ccomplex_t IT_0306 = s_23*s_34;
    const ccomplex_t IT_0307 = (-2)*IT_0306;
    const ccomplex_t IT_0308 = IT_0304 + IT_0307;
    const ccomplex_t IT_0309 = IT_0002*IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0303 + 0.5*IT_0305 + (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = s_24 + IT_0310;
    const ccomplex_t IT_0312 = m_N_1*IT_0311;
    const ccomplex_t IT_0313 = (-2)*IT_0312;
    const ccomplex_t IT_0314 = 2*IT_0173;
    const ccomplex_t IT_0315 = (-2)*IT_0172;
    const ccomplex_t IT_0316 = IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = 2*IT_0211;
    const ccomplex_t IT_0318 = (-2)*IT_0210;
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = IT_0261 + 2*IT_0262 + 2*IT_0263 + -IT_0264 + 2
      *IT_0265;
    const ccomplex_t IT_0321 = IT_0000*IT_0040;
    const ccomplex_t IT_0322 = -IT_0321;
    const ccomplex_t IT_0323 = IT_0293 + IT_0322;
    const ccomplex_t IT_0324 = IT_0002*IT_0323;
    const ccomplex_t IT_0325 = -IT_0040;
    const ccomplex_t IT_0326 = IT_0294 + IT_0325;
    const ccomplex_t IT_0327 = (-4)*IT_0326;
    const ccomplex_t IT_0328 = IT_0324 + IT_0327;
    const ccomplex_t IT_0329 = m_N_1*IT_0328;
    const ccomplex_t IT_0330 = 2*conj(IT_0116);
    const ccomplex_t IT_0331 = -IT_0215 + (-2)*IT_0219 + (-2)*IT_0222 + (-2)
      *IT_0225 + IT_0275;
    const ccomplex_t IT_0332 = s_13*s_23;
    const ccomplex_t IT_0333 = (-2)*IT_0332;
    const ccomplex_t IT_0334 = IT_0001 + IT_0333;
    const ccomplex_t IT_0335 = IT_0002*IT_0334;
    const ccomplex_t IT_0336 = (s_12 + -IT_0229)*(s_12 + (-0.5)*IT_0335);
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = 2*conj(IT_0282);
    const ccomplex_t IT_0339 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0340 = pow(m_Z, 4);
    const ccomplex_t IT_0341 = s_24*IT_0340;
    const ccomplex_t IT_0342 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0343 = (-2)*IT_0342;
    const ccomplex_t IT_0344 = IT_0341 + IT_0343;
    const ccomplex_t IT_0345 = IT_0005*IT_0344;
    const ccomplex_t IT_0346 = (-0.666666666666667)*IT_0339 + (
      -0.166666666666667)*IT_0345;
    const ccomplex_t IT_0347 = s_24 + IT_0346;
    const ccomplex_t IT_0348 = m_N_1*IT_0347;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = m_N_1*IT_0002;
    const ccomplex_t IT_0351 = s_14*IT_0340;
    const ccomplex_t IT_0352 = IT_0002*IT_0351;
    const ccomplex_t IT_0353 = (-2)*IT_0253 + (-0.5)*IT_0352;
    const ccomplex_t IT_0354 = IT_0251 + IT_0353;
    const ccomplex_t IT_0355 = IT_0350*IT_0354;
    const ccomplex_t IT_0356 = (-2)*IT_0355;
    const ccomplex_t IT_0357 = 2*conj(IT_0192);
    const ccomplex_t IT_0358 = IT_0005*IT_0341;
    const ccomplex_t IT_0359 = (-10)*s_24;
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = m_N_1*IT_0360;
    const ccomplex_t IT_0362 = IT_0112*IT_0224;
    const ccomplex_t IT_0363 = IT_0038*IT_0214;
    const ccomplex_t IT_0364 = IT_0087*IT_0218;
    const ccomplex_t IT_0365 = IT_0064*IT_0221;
    const ccomplex_t IT_0366 = IT_0362 + -IT_0363 + IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = 2*conj(IT_0366);
    const ccomplex_t IT_0368 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0369 = (-2)*IT_0368;
    const ccomplex_t IT_0370 = IT_0351 + IT_0369;
    const ccomplex_t IT_0371 = IT_0005*IT_0370;
    const ccomplex_t IT_0372 = (-0.666666666666667)*IT_0256 +
       0.166666666666667*IT_0371;
    const ccomplex_t IT_0373 = s_14 + IT_0372;
    const ccomplex_t IT_0374 = m_N_1*IT_0373;
    const ccomplex_t IT_0375 = 6*IT_0374;
    const ccomplex_t IT_0376 = IT_0178*IT_0214;
    const ccomplex_t IT_0377 = IT_0182*IT_0218;
    const ccomplex_t IT_0378 = IT_0186*IT_0221;
    const ccomplex_t IT_0379 = IT_0190*IT_0224;
    const ccomplex_t IT_0380 = -IT_0376 + IT_0377 + IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = 2*conj(IT_0380);
    const ccomplex_t IT_0382 = s_12*s_24;
    const ccomplex_t IT_0383 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0384 = s_14*IT_0000*IT_0040;
    const ccomplex_t IT_0385 = s_13*s_23*s_24;
    const ccomplex_t IT_0386 = s_13*s_34*IT_0040;
    const ccomplex_t IT_0387 = -IT_0384 + (-2)*IT_0385 + 2*IT_0386;
    const ccomplex_t IT_0388 = IT_0383 + IT_0387;
    const ccomplex_t IT_0389 = IT_0002*IT_0388;
    const ccomplex_t IT_0390 = (-0.5)*IT_0389;
    const ccomplex_t IT_0391 = s_14*IT_0040;
    const ccomplex_t IT_0392 = -IT_0391;
    const ccomplex_t IT_0393 = IT_0382 + IT_0390 + IT_0392;
    const ccomplex_t IT_0394 = 4*conj(IT_0331);
    const ccomplex_t IT_0395 = pow(s_34, 2);
    const ccomplex_t IT_0396 = IT_0005*IT_0395;
    const ccomplex_t IT_0397 = 0.5*IT_0396;
    const ccomplex_t IT_0398 = 1 + IT_0397;
    const ccomplex_t IT_0399 = s_12*IT_0398;
    const ccomplex_t IT_0400 = 4*IT_0173;
    const ccomplex_t IT_0401 = (-4)*IT_0172;
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = 4*IT_0211;
    const ccomplex_t IT_0404 = (-4)*IT_0210;
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = 4*IT_0210;
    const ccomplex_t IT_0407 = (-4)*conj(IT_0172);
    const ccomplex_t IT_0408 = (-4)*conj(IT_0210);
    const ccomplex_t IT_0409 = 1.125*IT_0040;
    const ccomplex_t IT_0410 = s_14*s_24;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = IT_0001 + IT_0411;
    const ccomplex_t IT_0413 = IT_0002*IT_0412;
    const ccomplex_t IT_0414 = -IT_0229 + (-0.25)*IT_0413;
    const ccomplex_t IT_0415 = s_12 + IT_0414;
    const ccomplex_t IT_0416 = m_N_1*IT_0415;
    const ccomplex_t IT_0417 = (-4)*IT_0416;
    const ccomplex_t IT_0418 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0419 = 0.5*IT_0240 + (-0.5)*IT_0321 + -IT_0418;
    const ccomplex_t IT_0420 = IT_0238 + IT_0419;
    const ccomplex_t IT_0421 = IT_0002*IT_0420;
    const ccomplex_t IT_0422 = -IT_0241;
    const ccomplex_t IT_0423 = IT_0040 + IT_0421 + IT_0422;
    const ccomplex_t IT_0424 = m_N_1*IT_0423;
    const ccomplex_t IT_0425 = (-2)*IT_0410;
    const ccomplex_t IT_0426 = IT_0001 + IT_0425;
    const ccomplex_t IT_0427 = IT_0002*IT_0426;
    const ccomplex_t IT_0428 = (-0.5)*IT_0427;
    const ccomplex_t IT_0429 = s_12 + IT_0428;
    const ccomplex_t IT_0430 = s_12*IT_0340;
    const ccomplex_t IT_0431 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0432 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0433 = s_14*s_23*s_34;
    const ccomplex_t IT_0434 = (-2)*IT_0431 + (-2)*IT_0432 + 4*IT_0433;
    const ccomplex_t IT_0435 = IT_0430 + IT_0434;
    const ccomplex_t IT_0436 = IT_0005*IT_0435;
    const ccomplex_t IT_0437 = 0.25*IT_0436;
    const ccomplex_t IT_0438 = (-0.5)*IT_0335;
    const ccomplex_t IT_0439 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0440 = (-2)*IT_0332 + -IT_0410 + 2*IT_0439;
    const ccomplex_t IT_0441 = IT_0001 + IT_0440;
    const ccomplex_t IT_0442 = IT_0002*IT_0441;
    const ccomplex_t IT_0443 = -IT_0229 + (-0.5)*IT_0442;
    const ccomplex_t IT_0444 = s_12 + IT_0443;
    const ccomplex_t IT_0445 = m_N_1*IT_0444;
    const ccomplex_t IT_0446 = (-2)*IT_0445;
    const ccomplex_t IT_0447 = (-2)*IT_0282;
    const ccomplex_t IT_0448 = 8*IT_0015*(IT_0172 + -IT_0173) + (-16)*IT_0192
      *IT_0409 + (-2)*IT_0226*IT_0417 + 4*IT_0331*IT_0424 + 8*IT_0116*(IT_0429 +
       IT_0437 + IT_0438) + IT_0446*IT_0447;
    const ccomplex_t IT_0449 = pow(m_Z, 6);
    const ccomplex_t IT_0450 = s_12*IT_0449;
    const ccomplex_t IT_0451 = s_13*s_23*IT_0340;
    const ccomplex_t IT_0452 = (-2)*IT_0451;
    const ccomplex_t IT_0453 = IT_0450 + IT_0452;
    const ccomplex_t IT_0454 = IT_0005*IT_0453;
    const ccomplex_t IT_0455 = 0.25*IT_0454;
    const ccomplex_t IT_0456 = IT_0001 + IT_0455;
    const ccomplex_t IT_0457 = 0.5*IT_0454;
    const ccomplex_t IT_0458 = (-4)*IT_0410;
    const ccomplex_t IT_0459 = IT_0001 + IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = IT_0002*IT_0435;
    const ccomplex_t IT_0461 = -IT_0460;
    const ccomplex_t IT_0462 = IT_0459 + IT_0461;
    const ccomplex_t IT_0463 = (-0.666666666666667)*IT_0309 +
       0.166666666666667*IT_0345;
    const ccomplex_t IT_0464 = s_24 + IT_0463;
    const ccomplex_t IT_0465 = m_N_1*IT_0464;
    const ccomplex_t IT_0466 = 6*IT_0465;
    const ccomplex_t IT_0467 = 2*IT_0116;
    const ccomplex_t IT_0468 = IT_0005*IT_0351;
    const ccomplex_t IT_0469 = (-10)*s_14;
    const ccomplex_t IT_0470 = IT_0468 + IT_0469;
    const ccomplex_t IT_0471 = m_N_1*IT_0470;
    const ccomplex_t IT_0472 = 2*IT_0192;
    const ccomplex_t IT_0473 = IT_0002*IT_0341;
    const ccomplex_t IT_0474 = (-2)*IT_0306 + (-0.5)*IT_0473;
    const ccomplex_t IT_0475 = IT_0304 + IT_0474;
    const ccomplex_t IT_0476 = IT_0350*IT_0475;
    const ccomplex_t IT_0477 = (-2)*IT_0476;
    const ccomplex_t IT_0478 = 2*IT_0366;
    const ccomplex_t IT_0479 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0480 = (-0.166666666666667)*IT_0371 + (
      -0.666666666666667)*IT_0479;
    const ccomplex_t IT_0481 = s_14 + IT_0480;
    const ccomplex_t IT_0482 = m_N_1*IT_0481;
    const ccomplex_t IT_0483 = (-6)*IT_0482;
    const ccomplex_t IT_0484 = 2*IT_0380;
    const ccomplex_t IT_0485 = s_14*IT_0238;
    const ccomplex_t IT_0486 = s_12*s_23*s_34;
    const ccomplex_t IT_0487 = 0.5*IT_0383 + (-0.5)*IT_0384 + -IT_0486;
    const ccomplex_t IT_0488 = IT_0485 + IT_0487;
    const ccomplex_t IT_0489 = IT_0002*IT_0488;
    const ccomplex_t IT_0490 = -IT_0382;
    const ccomplex_t IT_0491 = IT_0391 + IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = 4*IT_0331;
    const ccomplex_t IT_0493 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0494 = s_13*s_14*s_23;
    const ccomplex_t IT_0495 = s_24*IT_0000*IT_0040;
    const ccomplex_t IT_0496 = s_23*s_34*IT_0040;
    const ccomplex_t IT_0497 = (-2)*IT_0494 + -IT_0495 + 2*IT_0496;
    const ccomplex_t IT_0498 = IT_0493 + IT_0497;
    const ccomplex_t IT_0499 = IT_0002*IT_0498;
    const ccomplex_t IT_0500 = s_24*IT_0040;
    const ccomplex_t IT_0501 = 2*IT_0500;
    const ccomplex_t IT_0502 = s_12*s_14;
    const ccomplex_t IT_0503 = (-2)*IT_0502;
    const ccomplex_t IT_0504 = IT_0499 + IT_0501 + IT_0503;
    const ccomplex_t IT_0505 = IT_0276*IT_0456 + 0.5*IT_0274*IT_0462 + (-0.125
      )*IT_0466*IT_0467 + (-0.125)*IT_0471*IT_0472 + (-0.125)*IT_0477*IT_0478 + 
      (-0.125)*IT_0483*IT_0484 + (-0.125)*IT_0491*IT_0492 + (-0.125)*IT_0447
      *IT_0504;
    const ccomplex_t IT_0506 = (-2)*conj(IT_0282);
    const ccomplex_t IT_0507 = (-2)*conj(IT_0320);
    const ccomplex_t IT_0508 = -IT_0294;
    const ccomplex_t IT_0509 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_0510 = 0.5*IT_0293 + (-0.5)*IT_0321 + -IT_0509;
    const ccomplex_t IT_0511 = IT_0291 + IT_0510;
    const ccomplex_t IT_0512 = IT_0002*IT_0511;
    const ccomplex_t IT_0513 = IT_0040 + IT_0508 + IT_0512;
    const ccomplex_t IT_0514 = m_N_1*IT_0513;
    const ccomplex_t IT_0515 = IT_0240 + IT_0322;
    const ccomplex_t IT_0516 = IT_0002*IT_0515;
    const ccomplex_t IT_0517 = IT_0241 + IT_0325;
    const ccomplex_t IT_0518 = (-4)*IT_0517;
    const ccomplex_t IT_0519 = IT_0516 + IT_0518;
    const ccomplex_t IT_0520 = m_N_1*IT_0519;
    const ccomplex_t IT_0521 = -IT_0008 + -IT_0010;
    const ccomplex_t IT_0522 = IT_0007 + IT_0521;
    const ccomplex_t IT_0523 = IT_0006*IT_0522;
    const ccomplex_t IT_0524 = (-0.25)*IT_0523;
    const ccomplex_t IT_0525 = s_12 + IT_0004 + IT_0524;
    const ccomplex_t IT_0526 = s_13*s_24*s_34;
    const ccomplex_t IT_0527 = (-2)*IT_0431 + (-2)*IT_0432 + 4*IT_0526;
    const ccomplex_t IT_0528 = IT_0430 + IT_0527;
    const ccomplex_t IT_0529 = IT_0005*IT_0528;
    const ccomplex_t IT_0530 = 0.25*IT_0529;
    const ccomplex_t IT_0531 = IT_0429 + IT_0438 + IT_0530;
    const ccomplex_t IT_0532 = s_12*IT_0395;
    const ccomplex_t IT_0533 = (-2)*IT_0431 + (-2)*IT_0432 + 2*IT_0433 + 2
      *IT_0526 + (-2)*IT_0532;
    const ccomplex_t IT_0534 = IT_0430 + IT_0533;
    const ccomplex_t IT_0535 = IT_0005*IT_0534;
    const ccomplex_t IT_0536 = (-0.25)*IT_0535;
    const ccomplex_t IT_0537 = -IT_0229;
    const ccomplex_t IT_0538 = -IT_0228;
    const ccomplex_t IT_0539 = s_12 + IT_0536 + IT_0537 + IT_0538;
    const ccomplex_t IT_0540 = (-8)*IT_0116;
    const ccomplex_t IT_0541 = (-2)*IT_0395;
    const ccomplex_t IT_0542 = IT_0340 + IT_0541;
    const ccomplex_t IT_0543 = IT_0005*IT_0542;
    const ccomplex_t IT_0544 = 0.25*IT_0543;
    const ccomplex_t IT_0545 = 1 + IT_0544;
    const ccomplex_t IT_0546 = IT_0040*IT_0545;
    const ccomplex_t IT_0547 = (-8)*IT_0192;
    const ccomplex_t IT_0548 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_0549 = (-2)*IT_0332 + -IT_0410 + 2*IT_0548;
    const ccomplex_t IT_0550 = IT_0001 + IT_0549;
    const ccomplex_t IT_0551 = IT_0002*IT_0550;
    const ccomplex_t IT_0552 = -IT_0229 + (-0.5)*IT_0551;
    const ccomplex_t IT_0553 = s_12 + IT_0552;
    const ccomplex_t IT_0554 = m_N_1*IT_0553;
    const ccomplex_t IT_0555 = (-2)*IT_0554;
    const ccomplex_t IT_0556 = (-2)*IT_0226;
    const ccomplex_t IT_0557 = (-16)*IT_0380*IT_0409 + IT_0417*IT_0447 + 4
      *IT_0320*IT_0514 + 2*IT_0331*IT_0520 + (-8)*(IT_0172 + -IT_0173)*IT_0525 +
       8*IT_0366*IT_0531 + IT_0539*IT_0540 + IT_0546*IT_0547 + IT_0555*IT_0556;
    const ccomplex_t IT_0558 = IT_0002*IT_0040;
    const ccomplex_t IT_0559 = (-0.25)*IT_0340;
    const ccomplex_t IT_0560 = IT_0395 + IT_0559;
    const ccomplex_t IT_0561 = IT_0558*IT_0560;
    const ccomplex_t IT_0562 = (-4)*IT_0561;
    const ccomplex_t IT_0563 = 4*IT_0226;
    const ccomplex_t IT_0564 = (-2)*IT_0320;
    const ccomplex_t IT_0565 = 9*IT_0321;
    const ccomplex_t IT_0566 = (-2)*IT_0274;
    const ccomplex_t IT_0567 = 4*IT_0267*IT_0456 + 2*IT_0266*IT_0462 + (-0.5)
      *IT_0467*IT_0471 + (-0.5)*IT_0466*IT_0472 + (-0.5)*IT_0478*IT_0483 + (-0.5
      )*IT_0477*IT_0484 + IT_0276*IT_0562 + (-0.5)*IT_0491*IT_0563 + (-0.5)
      *IT_0504*IT_0564 + (-0.5)*IT_0565*IT_0566;
    const ccomplex_t IT_0568 = 4*conj(IT_0226);
    const ccomplex_t IT_0569 = (-2)*conj(IT_0274);
    const ccomplex_t IT_0570 = (-2)*IT_0331;
    const ccomplex_t IT_0571 = -IT_0211;
    const ccomplex_t IT_0572 = (-16)*IT_0366*IT_0409 + 4*IT_0282*IT_0514 + 2
      *IT_0226*IT_0520 + 8*IT_0380*IT_0531 + IT_0540*IT_0546 + IT_0539*IT_0547 +
       IT_0417*IT_0564 + IT_0555*IT_0570 + (-8)*IT_0525*(IT_0210 + IT_0571);
    const ccomplex_t IT_0573 = 2*IT_0282;
    const ccomplex_t IT_0574 = (-16)*IT_0116*IT_0409 + (-2)*IT_0331*IT_0417 +
       4*IT_0226*IT_0424 + 8*IT_0192*(IT_0429 + IT_0437 + IT_0438) + IT_0446
      *IT_0564 + 8*IT_0015*(IT_0210 + IT_0571) + IT_0329*IT_0573;
    const ccomplex_t IT_0575 = IT_0040 + IT_0508;
    const ccomplex_t IT_0576 = IT_0040*IT_0575;
    const ccomplex_t IT_0577 = -IT_0576;
    const ccomplex_t IT_0578 = s_24*IT_0291;
    const ccomplex_t IT_0579 = s_12*s_13*s_34;
    const ccomplex_t IT_0580 = -IT_0495 + 2*IT_0578 + (-2)*IT_0579;
    const ccomplex_t IT_0581 = IT_0493 + IT_0580;
    const ccomplex_t IT_0582 = IT_0002*IT_0581;
    const ccomplex_t IT_0583 = (-0.5)*IT_0582;
    const ccomplex_t IT_0584 = -IT_0500;
    const ccomplex_t IT_0585 = IT_0502 + IT_0583 + IT_0584;
    const ccomplex_t IT_0586 = s_12 + IT_0537;
    const ccomplex_t IT_0587 = IT_0040*IT_0586;
    const ccomplex_t IT_0588 = -IT_0226;
    const ccomplex_t IT_0589 = 0.166666666666667*IT_0298;
    const ccomplex_t IT_0590 = (IT_0040 + -IT_0294)*(s_12 + (-0.5)*IT_0335);
    const ccomplex_t IT_0591 = 2*IT_0590;
    const ccomplex_t IT_0592 = 0.166666666666667*IT_0591;
    const ccomplex_t IT_0593 = -IT_0331;
    const ccomplex_t IT_0594 = 0.166666666666667*IT_0285*IT_0319 +
       0.333333333333333*IT_0226*IT_0337 + 0.166666666666667*IT_0329*IT_0467 +
       0.666666666666667*IT_0366*IT_0514 + IT_0282*IT_0577 + (-0.666666666666667
      )*IT_0267*IT_0585 + IT_0316*IT_0589 + IT_0564*IT_0592 + IT_0587*IT_0593;
    const ccomplex_t IT_0595 = IT_0040*IT_0398;
    const ccomplex_t IT_0596 = 4*conj(IT_0173);
    const ccomplex_t IT_0597 = 4*conj(IT_0172);
    const ccomplex_t IT_0598 = (-4)*conj(IT_0173);
    const ccomplex_t IT_0599 = (-4)*IT_0173;
    const ccomplex_t IT_0600 = 0.25*IT_0599;
    const ccomplex_t IT_0601 = 4*conj(IT_0210)*(IT_0366 + 0.25*IT_0400 + 0.25
      *IT_0401) + conj(IT_0192)*IT_0402 + conj(IT_0116)*IT_0405 + conj(IT_0366)
      *IT_0406 + IT_0116*IT_0408 + (IT_0192 + IT_0210)*(IT_0407 + IT_0596) +
       IT_0380*(IT_0597 + IT_0598) + (-4)*IT_0211*(conj(IT_0366) + (-0.25)
      *IT_0597 + (-0.25)*IT_0598) + 4*conj(IT_0380)*(IT_0172 + IT_0600) + 4*conj
      (IT_0211)*(IT_0116 + IT_0172 + -IT_0366 + IT_0600);
    const ccomplex_t IT_0602 = IT_0002*IT_0528;
    const ccomplex_t IT_0603 = -IT_0602;
    const ccomplex_t IT_0604 = IT_0459 + IT_0603;
    const ccomplex_t IT_0605 = (-0.5)*IT_0349;
    const ccomplex_t IT_0606 = (-0.5)*IT_0375;
    const ccomplex_t IT_0607 = (-0.5)*IT_0361;
    const ccomplex_t IT_0608 = (-0.5)*IT_0393;
    const ccomplex_t IT_0609 = 4*IT_0266*IT_0456 + (-0.5)*IT_0356*IT_0467 +
       IT_0276*IT_0565 + (-0.5)*IT_0562*IT_0566 + 2*IT_0320*IT_0585 + 2*IT_0267
      *IT_0604 + IT_0472*IT_0605 + IT_0478*IT_0606 + IT_0484*IT_0607 + IT_0563
      *IT_0608;
    const ccomplex_t IT_0610 = 4*IT_0274*IT_0456 + (-0.5)*IT_0356*IT_0472 +
       IT_0266*IT_0562 + IT_0267*IT_0565 + 2*IT_0282*IT_0585 + 2*IT_0276*IT_0604
       + IT_0467*IT_0605 + IT_0484*IT_0606 + IT_0478*IT_0607 + IT_0492*IT_0608;
    const ccomplex_t IT_0611 = IT_0040 + IT_0422;
    const ccomplex_t IT_0612 = IT_0040*IT_0611;
    const ccomplex_t IT_0613 = -IT_0612;
    const ccomplex_t IT_0614 = (-0.166666666666667)*IT_0245;
    const ccomplex_t IT_0615 = (IT_0040 + -IT_0241)*(s_12 + (-0.5)*IT_0335);
    const ccomplex_t IT_0616 = 2*IT_0615;
    const ccomplex_t IT_0617 = (-0.166666666666667)*IT_0616;
    const ccomplex_t IT_0618 = (-0.166666666666667)*IT_0232*IT_0271 + (
      -0.333333333333333)*IT_0320*IT_0337 + 0.333333333333333*IT_0116*IT_0417 + 
      (-0.666666666666667)*IT_0192*IT_0424 + (-0.333333333333333)*IT_0380
      *IT_0520 + 0.333333333333333*IT_0366*IT_0555 + IT_0282*IT_0587 + IT_0593
      *IT_0613 + IT_0280*IT_0614 + IT_0556*IT_0617;
    const ccomplex_t IT_0619 = (-0.166666666666667)*IT_0232*IT_0280 +
       0.333333333333333*IT_0192*IT_0417 + (-0.666666666666667)*IT_0116*IT_0424 
      + (-0.333333333333333)*IT_0366*IT_0520 + 0.333333333333333*IT_0380*IT_0555
       + (-0.166666666666667)*IT_0337*IT_0573 + IT_0320*IT_0587 + IT_0588
      *IT_0613 + IT_0271*IT_0614 + IT_0570*IT_0617;
    const ccomplex_t IT_0620 = (-8)*conj(IT_0192);
    const ccomplex_t IT_0621 = (-0.125)*IT_0546;
    const ccomplex_t IT_0622 = (-8)*conj(IT_0116);
    const ccomplex_t IT_0623 = (-0.125)*IT_0539;
    const ccomplex_t IT_0624 = 8*IT_0015*(IT_0116*IT_0174 + IT_0192*IT_0212) +
       IT_0226*(IT_0232*IT_0235 + IT_0245*IT_0248) + IT_0260*(IT_0268*IT_0271 +
       IT_0235*IT_0272 + IT_0277*IT_0280 + IT_0248*IT_0281) + IT_0282*(IT_0285
      *IT_0288 + IT_0298*IT_0301) + IT_0313*(IT_0281*IT_0288 + IT_0272*IT_0301 +
       IT_0277*IT_0316 + IT_0268*IT_0319) + IT_0320*(IT_0288*IT_0298 + IT_0285
      *IT_0301 + IT_0329*IT_0330) + IT_0331*(IT_0235*IT_0245 + IT_0232*IT_0248 +
       IT_0337*IT_0338) + IT_0276*(IT_0330*IT_0349 + IT_0356*IT_0357 + IT_0361
      *IT_0367 + IT_0375*IT_0381 + IT_0393*IT_0394) + 4*IT_0399*(IT_0172*conj
      (IT_0172) + 0.25*conj(IT_0173)*IT_0402 + 0.25*conj(IT_0211)*IT_0405 + 0.25
      *conj(IT_0210)*IT_0406 + 0.25*IT_0173*IT_0407 + 0.25*IT_0211*IT_0408) +
       conj(IT_0116)*IT_0448 + (-8)*conj(IT_0274)*IT_0505 + IT_0116*IT_0446
      *IT_0506 + IT_0274*(IT_0330*IT_0466 + IT_0357*IT_0471 + IT_0367*IT_0477 +
       IT_0381*IT_0483 + IT_0394*IT_0491 + IT_0504*IT_0506) + IT_0192*(IT_0329
      *IT_0338 + IT_0446*IT_0507) + conj(IT_0366)*IT_0557 + (-2)*conj(IT_0266)
      *IT_0567 + IT_0267*(IT_0330*IT_0356 + IT_0349*IT_0357 + IT_0367*IT_0375 +
       IT_0361*IT_0381 + IT_0393*IT_0568 + IT_0562*IT_0569) + IT_0266*(IT_0357
      *IT_0466 + IT_0330*IT_0471 + IT_0381*IT_0477 + IT_0367*IT_0483 + IT_0504
      *IT_0507 + IT_0491*IT_0568 + IT_0565*IT_0569) + conj(IT_0380)*IT_0572 +
       conj(IT_0192)*IT_0574 + 6*conj(IT_0282)*(0.166666666666667*IT_0285
      *IT_0316 + 0.666666666666667*IT_0380*IT_0514 + IT_0320*IT_0577 + (
      -0.666666666666667)*IT_0276*IT_0585 + IT_0587*IT_0588 + IT_0319*IT_0589 +
       IT_0447*IT_0592) + 6*conj(IT_0320)*IT_0594 + IT_0595*IT_0601 + (-2)*conj
      (IT_0267)*IT_0609 + (-2)*conj(IT_0276)*IT_0610 + (-6)*conj(IT_0226)
      *IT_0618 + (-6)*conj(IT_0331)*IT_0619 + (-8)*IT_0380*((-0.125)*IT_0417
      *IT_0507 + IT_0212*IT_0525 + IT_0621*IT_0622 + IT_0620*IT_0623) + (-8)
      *IT_0366*((-0.125)*IT_0417*IT_0506 + IT_0174*IT_0525 + IT_0620*IT_0621 +
       IT_0622*IT_0623);
    return create_ccomplex_return(IT_0624);
}

