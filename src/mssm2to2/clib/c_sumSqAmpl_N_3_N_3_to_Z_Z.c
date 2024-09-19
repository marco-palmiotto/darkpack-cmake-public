#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_Z_Z(
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
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_Z, -4);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, 2);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = pow(m_Z, -2);
    const ccomplex_t IT_0006 = s_14*IT_0002;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = s_13*s_34;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = IT_0006 + IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0004 + 0.5*IT_0007 + (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = s_14 + IT_0012;
    const ccomplex_t IT_0014 = m_N_3*IT_0013;
    const ccomplex_t IT_0015 = (-2)*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0018, -1);
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0025;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + IT_0024 + 
      -IT_0026 + -IT_0027);
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0002 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = IT_0023*IT_0035;
    const ccomplex_t IT_0038 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0039 = IT_0019*IT_0038;
    const ccomplex_t IT_0040 = IT_0023*IT_0038;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0036 + IT_0037 + 
      -IT_0039 + -IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0044 = IT_0019*IT_0043;
    const ccomplex_t IT_0045 = IT_0023*IT_0043;
    const ccomplex_t IT_0046 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0047 = IT_0019*IT_0046;
    const ccomplex_t IT_0048 = IT_0023*IT_0046;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0044 + IT_0045 + 
      -IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = pow(m_N_3, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0002 + IT_0052 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0051*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0028;
    const ccomplex_t IT_0058 = IT_0029*IT_0057;
    const ccomplex_t IT_0059 = IT_0032*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0002 + IT_0052 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0051*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0059 + IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0002 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0058*IT_0065;
    const ccomplex_t IT_0067 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0068 = IT_0019*IT_0067;
    const ccomplex_t IT_0069 = IT_0023*IT_0067;
    const ccomplex_t IT_0070 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0071 = IT_0019*IT_0070;
    const ccomplex_t IT_0072 = IT_0023*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0076 = IT_0019*IT_0075;
    const ccomplex_t IT_0077 = IT_0023*IT_0075;
    const ccomplex_t IT_0078 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0079 = IT_0019*IT_0078;
    const ccomplex_t IT_0080 = IT_0023*IT_0078;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0076 + IT_0077 + 
      -IT_0079 + -IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_23 + IT_0002 + IT_0052 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0066 + -IT_0086;
    const ccomplex_t IT_0088 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0089 = IT_0019*IT_0088;
    const ccomplex_t IT_0090 = IT_0023*IT_0088;
    const ccomplex_t IT_0091 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0092 = IT_0019*IT_0091;
    const ccomplex_t IT_0093 = IT_0023*IT_0091;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0089 + IT_0090 + 
      -IT_0092 + -IT_0093);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0097 = IT_0019*IT_0096;
    const ccomplex_t IT_0098 = IT_0023*IT_0096;
    const ccomplex_t IT_0099 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0100 = IT_0019*IT_0099;
    const ccomplex_t IT_0101 = IT_0023*IT_0099;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0097 + IT_0098 + 
      -IT_0100 + -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0095*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_13 + IT_0002 + IT_0052 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0002 + IT_0052 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0104*IT_0109;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0002 + IT_0052 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0083*IT_0113;
    const ccomplex_t IT_0115 = cpow(IT_0057, 2);
    const ccomplex_t IT_0116 = IT_0065*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = cos(alpha);
    const ccomplex_t IT_0119 = IT_0017*IT_0118;
    const ccomplex_t IT_0120 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = sin(alpha);
    const ccomplex_t IT_0123 = IT_0017*IT_0122;
    const ccomplex_t IT_0124 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = IT_0022*IT_0118;
    const ccomplex_t IT_0127 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = IT_0022*IT_0122;
    const ccomplex_t IT_0130 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0131 = IT_0129*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0121 + -IT_0125 + 
      -IT_0128 + IT_0131);
    const ccomplex_t IT_0133 = cos(beta);
    const ccomplex_t IT_0134 = cpow(IT_0016, -4);
    const ccomplex_t IT_0135 = IT_0018*IT_0118*IT_0133*IT_0134;
    const ccomplex_t IT_0136 = sin(beta);
    const ccomplex_t IT_0137 = IT_0018*IT_0122*IT_0134*IT_0136;
    const ccomplex_t IT_0138 = cpow(IT_0018, -3);
    const ccomplex_t IT_0139 = IT_0118*IT_0133*IT_0138;
    const ccomplex_t IT_0140 = IT_0122*IT_0136*IT_0138;
    const ccomplex_t IT_0141 = cpow(IT_0016, -2);
    const ccomplex_t IT_0142 = IT_0022*IT_0118*IT_0133*IT_0141;
    const ccomplex_t IT_0143 = IT_0022*IT_0122*IT_0136*IT_0141;
    const ccomplex_t IT_0144 = cpow(IT_0018, -2);
    const ccomplex_t IT_0145 = cpow(IT_0141 + IT_0144, -1);
    const ccomplex_t IT_0146 = (IT_0135 + IT_0137 + IT_0139 + IT_0140 + 2
      *IT_0142 + 2*IT_0143)*IT_0145;
    const ccomplex_t IT_0147 = m_W*e_em;
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = cpow((-2)*s_12 + (-2)*IT_0052 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0132*IT_0152;
    const ccomplex_t IT_0154 = IT_0119*IT_0124;
    const ccomplex_t IT_0155 = IT_0120*IT_0123;
    const ccomplex_t IT_0156 = IT_0126*IT_0130;
    const ccomplex_t IT_0157 = IT_0127*IT_0129;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + IT_0155 + 
      -IT_0156 + -IT_0157);
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0018*IT_0118*IT_0134*IT_0136;
    const ccomplex_t IT_0161 = IT_0018*IT_0122*IT_0133*IT_0134;
    const ccomplex_t IT_0162 = IT_0118*IT_0136*IT_0138;
    const ccomplex_t IT_0163 = IT_0122*IT_0133*IT_0138;
    const ccomplex_t IT_0164 = IT_0022*IT_0118*IT_0136*IT_0141;
    const ccomplex_t IT_0165 = IT_0022*IT_0122*IT_0133*IT_0141;
    const ccomplex_t IT_0166 = IT_0145*(IT_0160 + -IT_0161 + IT_0162 + 
      -IT_0163 + 2*IT_0164 + (-2)*IT_0165);
    const ccomplex_t IT_0167 = IT_0147*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0052 + -reg_prop, -1);
    const ccomplex_t IT_0170 = IT_0168*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = IT_0159*IT_0171;
    const ccomplex_t IT_0173 = -IT_0153 + (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0173;
    const ccomplex_t IT_0175 = 0.5*IT_0172;
    const ccomplex_t IT_0176 = (-2)*IT_0175;
    const ccomplex_t IT_0177 = IT_0174 + IT_0176;
    const ccomplex_t IT_0178 = 2*conj(IT_0173);
    const ccomplex_t IT_0179 = (-2)*conj(IT_0175);
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = -IT_0086 + -IT_0110;
    const ccomplex_t IT_0182 = -IT_0030;
    const ccomplex_t IT_0183 = IT_0058 + IT_0182;
    const ccomplex_t IT_0184 = IT_0065*IT_0183;
    const ccomplex_t IT_0185 = (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = 0.5*IT_0059;
    const ccomplex_t IT_0187 = IT_0032*IT_0115;
    const ccomplex_t IT_0188 = IT_0062 + (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = IT_0107 + IT_0114;
    const ccomplex_t IT_0190 = conj(IT_0188) + conj(IT_0189);
    const ccomplex_t IT_0191 = conj(IT_0186) + IT_0190;
    const ccomplex_t IT_0192 = conj(IT_0056) + conj(IT_0181) + conj(IT_0185) +
       IT_0191;
    const ccomplex_t IT_0193 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0194 = IT_0119*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0196 = IT_0123*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0198 = IT_0126*IT_0197;
    const ccomplex_t IT_0199 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0200 = IT_0129*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0194 + -IT_0196 + 
      -IT_0198 + IT_0200);
    const ccomplex_t IT_0202 = IT_0152*IT_0201;
    const ccomplex_t IT_0203 = IT_0119*IT_0195;
    const ccomplex_t IT_0204 = IT_0123*IT_0193;
    const ccomplex_t IT_0205 = IT_0126*IT_0199;
    const ccomplex_t IT_0206 = IT_0129*IT_0197;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0203 + IT_0204 + 
      -IT_0205 + -IT_0206);
    const ccomplex_t IT_0208 = -IT_0207;
    const ccomplex_t IT_0209 = IT_0171*IT_0208;
    const ccomplex_t IT_0210 = -IT_0202 + (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = 2*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0209;
    const ccomplex_t IT_0213 = (-2)*IT_0212;
    const ccomplex_t IT_0214 = IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = IT_0188 + IT_0189;
    const ccomplex_t IT_0216 = IT_0186 + IT_0215;
    const ccomplex_t IT_0217 = 2*conj(IT_0210);
    const ccomplex_t IT_0218 = (-2)*conj(IT_0212);
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = (conj(IT_0034) + conj(IT_0056) + conj(IT_0063) 
      + conj(IT_0087) + conj(IT_0107) + conj(IT_0111) + conj(IT_0114) + conj
      (IT_0117))*IT_0177 + (IT_0034 + IT_0056 + IT_0063 + IT_0087 + IT_0107 +
       IT_0111 + IT_0114 + IT_0117)*IT_0180 + IT_0192*IT_0214 + (IT_0056 +
       IT_0181 + IT_0185 + IT_0216)*IT_0219;
    const ccomplex_t IT_0221 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0222 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0223 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0224 = -IT_0222 + -IT_0223;
    const ccomplex_t IT_0225 = s_12 + IT_0221 + IT_0224;
    const ccomplex_t IT_0226 = m_N_3*IT_0225;
    const ccomplex_t IT_0227 = IT_0033 + -IT_0059 + (-2)*IT_0062;
    const ccomplex_t IT_0228 = (-2)*IT_0107 + (-2)*IT_0114;
    const ccomplex_t IT_0229 = conj(IT_0227) + conj(IT_0228);
    const ccomplex_t IT_0230 = IT_0227 + IT_0228;
    const ccomplex_t IT_0231 = -IT_0059 + (-2)*IT_0062 + IT_0187;
    const ccomplex_t IT_0232 = conj(IT_0228) + conj(IT_0231);
    const ccomplex_t IT_0233 = IT_0228 + IT_0231;
    const ccomplex_t IT_0234 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0235 = IT_0052 + IT_0234;
    const ccomplex_t IT_0236 = pow(s_23, 2);
    const ccomplex_t IT_0237 = IT_0005*IT_0236;
    const ccomplex_t IT_0238 = pow(s_24, 2);
    const ccomplex_t IT_0239 = IT_0005*IT_0238;
    const ccomplex_t IT_0240 = -IT_0237 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0235 + IT_0240;
    const ccomplex_t IT_0242 = m_N_3*IT_0241;
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = s_23*IT_0002;
    const ccomplex_t IT_0245 = IT_0001*IT_0244;
    const ccomplex_t IT_0246 = s_24*IT_0002;
    const ccomplex_t IT_0247 = IT_0005*IT_0246;
    const ccomplex_t IT_0248 = s_23*s_34;
    const ccomplex_t IT_0249 = (-2)*IT_0248;
    const ccomplex_t IT_0250 = IT_0246 + IT_0249;
    const ccomplex_t IT_0251 = IT_0005*IT_0250;
    const ccomplex_t IT_0252 = (-0.5)*IT_0245 + 0.5*IT_0247 + (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = s_24 + IT_0252;
    const ccomplex_t IT_0254 = m_N_3*IT_0253;
    const ccomplex_t IT_0255 = (-2)*IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0175;
    const ccomplex_t IT_0257 = (-2)*IT_0173;
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = 2*conj(IT_0175);
    const ccomplex_t IT_0260 = (-2)*conj(IT_0173);
    const ccomplex_t IT_0261 = IT_0259 + IT_0260;
    const ccomplex_t IT_0262 = 2*IT_0212;
    const ccomplex_t IT_0263 = (-2)*IT_0210;
    const ccomplex_t IT_0264 = IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = 2*conj(IT_0212);
    const ccomplex_t IT_0266 = (-2)*conj(IT_0210);
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = IT_0192*IT_0258 + (IT_0056 + IT_0181 + IT_0185 
      + IT_0216)*IT_0261 + (conj(IT_0034) + conj(IT_0056) + conj(IT_0063) + conj
      (IT_0087) + conj(IT_0107) + conj(IT_0111) + conj(IT_0114) + conj(IT_0117))
      *IT_0264 + (IT_0034 + IT_0056 + IT_0063 + IT_0087 + IT_0107 + IT_0111 +
       IT_0114 + IT_0117)*IT_0267;
    const ccomplex_t IT_0269 = pow(m_Z, 4);
    const ccomplex_t IT_0270 = s_14*IT_0269;
    const ccomplex_t IT_0271 = IT_0000*IT_0270;
    const ccomplex_t IT_0272 = 0.666666666666667*IT_0007 + (-0.166666666666667
      )*IT_0271;
    const ccomplex_t IT_0273 = s_14 + IT_0272;
    const ccomplex_t IT_0274 = m_N_3*IT_0273;
    const ccomplex_t IT_0275 = (-6)*IT_0274;
    const ccomplex_t IT_0276 = conj(IT_0181) + conj(IT_0185);
    const ccomplex_t IT_0277 = cpow(IT_0074, 2);
    const ccomplex_t IT_0278 = m_N_1*IT_0277;
    const ccomplex_t IT_0279 = m_N_3*IT_0083;
    const ccomplex_t IT_0280 = IT_0085*(IT_0278 + -IT_0279);
    const ccomplex_t IT_0281 = (-2)*IT_0280;
    const ccomplex_t IT_0282 = m_N_3*IT_0058;
    const ccomplex_t IT_0283 = m_N_3*IT_0030;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = IT_0282 + IT_0284;
    const ccomplex_t IT_0286 = IT_0065*IT_0285;
    const ccomplex_t IT_0287 = cpow(IT_0095, 2);
    const ccomplex_t IT_0288 = m_N_2*IT_0287;
    const ccomplex_t IT_0289 = m_N_3*IT_0104;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = IT_0288 + IT_0290;
    const ccomplex_t IT_0292 = IT_0109*IT_0291;
    const ccomplex_t IT_0293 = cpow(IT_0050, 2);
    const ccomplex_t IT_0294 = m_N_4*IT_0293;
    const ccomplex_t IT_0295 = IT_0054*IT_0294;
    const ccomplex_t IT_0296 = (-0.5)*IT_0286 + IT_0292 + IT_0295;
    const ccomplex_t IT_0297 = (-2)*IT_0296;
    const ccomplex_t IT_0298 = m_N_3*IT_0115;
    const ccomplex_t IT_0299 = -IT_0065*(IT_0282 + -IT_0298);
    const ccomplex_t IT_0300 = 0.5*IT_0299;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0281 + IT_0297 + IT_0301;
    const ccomplex_t IT_0303 = IT_0181 + IT_0185;
    const ccomplex_t IT_0304 = (-2)*conj(IT_0280);
    const ccomplex_t IT_0305 = (-2)*conj(IT_0296);
    const ccomplex_t IT_0306 = (-2)*conj(IT_0300);
    const ccomplex_t IT_0307 = IT_0304 + IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0087) + conj(IT_0111) + conj(IT_0117);
    const ccomplex_t IT_0309 = -IT_0298;
    const ccomplex_t IT_0310 = IT_0282 + IT_0309;
    const ccomplex_t IT_0311 = IT_0065*IT_0310;
    const ccomplex_t IT_0312 = cpow(IT_0082, 2);
    const ccomplex_t IT_0313 = m_N_1*IT_0312;
    const ccomplex_t IT_0314 = -IT_0279;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0085*IT_0315;
    const ccomplex_t IT_0317 = cpow(IT_0103, 2);
    const ccomplex_t IT_0318 = m_N_2*IT_0317;
    const ccomplex_t IT_0319 = IT_0290 + IT_0318;
    const ccomplex_t IT_0320 = IT_0109*IT_0319;
    const ccomplex_t IT_0321 = cpow(IT_0042, 2);
    const ccomplex_t IT_0322 = m_N_4*IT_0321;
    const ccomplex_t IT_0323 = IT_0054*IT_0322;
    const ccomplex_t IT_0324 = (-0.5)*IT_0311 + IT_0316 + IT_0320 + IT_0323;
    const ccomplex_t IT_0325 = (-2)*IT_0324;
    const ccomplex_t IT_0326 = -IT_0065*(IT_0282 + -IT_0283);
    const ccomplex_t IT_0327 = 0.5*IT_0326;
    const ccomplex_t IT_0328 = (-2)*IT_0327;
    const ccomplex_t IT_0329 = IT_0325 + IT_0328;
    const ccomplex_t IT_0330 = IT_0087 + IT_0111 + IT_0117;
    const ccomplex_t IT_0331 = (-2)*conj(IT_0324);
    const ccomplex_t IT_0332 = (-2)*conj(IT_0327);
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = s_12*s_24*IT_0002;
    const ccomplex_t IT_0335 = s_13*s_23*s_24;
    const ccomplex_t IT_0336 = s_14*IT_0002*IT_0052;
    const ccomplex_t IT_0337 = s_13*s_34*IT_0052;
    const ccomplex_t IT_0338 = (-2)*IT_0335 + -IT_0336 + 2*IT_0337;
    const ccomplex_t IT_0339 = IT_0334 + IT_0338;
    const ccomplex_t IT_0340 = IT_0005*IT_0339;
    const ccomplex_t IT_0341 = s_14*IT_0052;
    const ccomplex_t IT_0342 = 2*IT_0341;
    const ccomplex_t IT_0343 = s_12*s_24;
    const ccomplex_t IT_0344 = (-2)*IT_0343;
    const ccomplex_t IT_0345 = IT_0340 + IT_0342 + IT_0344;
    const ccomplex_t IT_0346 = (-2)*IT_0231;
    const ccomplex_t IT_0347 = (-2)*IT_0228;
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = (-2)*conj(IT_0231);
    const ccomplex_t IT_0350 = (-2)*conj(IT_0228);
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = conj(IT_0034) + conj(IT_0063);
    const ccomplex_t IT_0353 = conj(IT_0107) + conj(IT_0114) + IT_0352;
    const ccomplex_t IT_0354 = (-2)*IT_0227;
    const ccomplex_t IT_0355 = IT_0347 + IT_0354;
    const ccomplex_t IT_0356 = IT_0034 + IT_0063;
    const ccomplex_t IT_0357 = IT_0107 + IT_0114 + IT_0356;
    const ccomplex_t IT_0358 = (-2)*conj(IT_0227);
    const ccomplex_t IT_0359 = IT_0350 + IT_0358;
    const ccomplex_t IT_0360 = s_24*IT_0269;
    const ccomplex_t IT_0361 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0362 = (-2)*IT_0361;
    const ccomplex_t IT_0363 = IT_0360 + IT_0362;
    const ccomplex_t IT_0364 = IT_0000*IT_0363;
    const ccomplex_t IT_0365 = (-0.666666666666667)*IT_0251 +
       0.166666666666667*IT_0364;
    const ccomplex_t IT_0366 = s_24 + IT_0365;
    const ccomplex_t IT_0367 = m_N_3*IT_0366;
    const ccomplex_t IT_0368 = 6*IT_0367;
    const ccomplex_t IT_0369 = m_N_3*IT_0051;
    const ccomplex_t IT_0370 = IT_0054*IT_0369;
    const ccomplex_t IT_0371 = 2*IT_0370;
    const ccomplex_t IT_0372 = conj(IT_0056)*IT_0371;
    const ccomplex_t IT_0373 = IT_0368*IT_0372;
    const ccomplex_t IT_0374 = (-10)*s_14;
    const ccomplex_t IT_0375 = IT_0271 + IT_0374;
    const ccomplex_t IT_0376 = m_N_3*IT_0375;
    const ccomplex_t IT_0377 = IT_0372*IT_0376;
    const ccomplex_t IT_0378 = 2*conj(IT_0370);
    const ccomplex_t IT_0379 = IT_0056*IT_0378;
    const ccomplex_t IT_0380 = IT_0368*IT_0379;
    const ccomplex_t IT_0381 = IT_0376*IT_0379;
    const ccomplex_t IT_0382 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0383 = (-0.166666666666667)*IT_0364 + (
      -0.666666666666667)*IT_0382;
    const ccomplex_t IT_0384 = s_24 + IT_0383;
    const ccomplex_t IT_0385 = m_N_3*IT_0384;
    const ccomplex_t IT_0386 = (-6)*IT_0385;
    const ccomplex_t IT_0387 = conj(IT_0188)*IT_0371;
    const ccomplex_t IT_0388 = conj(IT_0034)*IT_0371;
    const ccomplex_t IT_0389 = conj(IT_0189)*IT_0371;
    const ccomplex_t IT_0390 = conj(IT_0186)*IT_0371;
    const ccomplex_t IT_0391 = conj(IT_0063)*IT_0371;
    const ccomplex_t IT_0392 = conj(IT_0107)*IT_0371;
    const ccomplex_t IT_0393 = conj(IT_0114)*IT_0371;
    const ccomplex_t IT_0394 = IT_0188*IT_0378;
    const ccomplex_t IT_0395 = IT_0034*IT_0378;
    const ccomplex_t IT_0396 = IT_0189*IT_0378;
    const ccomplex_t IT_0397 = IT_0186*IT_0378;
    const ccomplex_t IT_0398 = IT_0063*IT_0378;
    const ccomplex_t IT_0399 = IT_0107*IT_0378;
    const ccomplex_t IT_0400 = IT_0114*IT_0378;
    const ccomplex_t IT_0401 = IT_0191*IT_0329 + IT_0216*IT_0333 + IT_0302
      *IT_0353 + IT_0307*IT_0357 + IT_0387 + IT_0388 + IT_0389 + IT_0390 +
       IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 +
       IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = m_N_3*IT_0005;
    const ccomplex_t IT_0403 = IT_0005*IT_0270;
    const ccomplex_t IT_0404 = (-2)*IT_0008 + (-0.5)*IT_0403;
    const ccomplex_t IT_0405 = IT_0006 + IT_0404;
    const ccomplex_t IT_0406 = IT_0402*IT_0405;
    const ccomplex_t IT_0407 = (-2)*IT_0406;
    const ccomplex_t IT_0408 = IT_0191*IT_0302 + IT_0216*IT_0307 + IT_0329
      *IT_0353 + IT_0333*IT_0357 + IT_0387 + IT_0388 + IT_0389 + IT_0390 +
       IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 +
       IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0409 = s_12*s_14*IT_0002;
    const ccomplex_t IT_0410 = pow(s_13, 2);
    const ccomplex_t IT_0411 = s_24*IT_0410;
    const ccomplex_t IT_0412 = s_12*s_13*s_34;
    const ccomplex_t IT_0413 = s_24*IT_0002*IT_0052;
    const ccomplex_t IT_0414 = 2*IT_0411 + (-2)*IT_0412 + -IT_0413;
    const ccomplex_t IT_0415 = IT_0409 + IT_0414;
    const ccomplex_t IT_0416 = IT_0005*IT_0415;
    const ccomplex_t IT_0417 = s_24*IT_0052;
    const ccomplex_t IT_0418 = 2*IT_0417;
    const ccomplex_t IT_0419 = s_12*s_14;
    const ccomplex_t IT_0420 = (-2)*IT_0419;
    const ccomplex_t IT_0421 = IT_0416 + IT_0418 + IT_0420;
    const ccomplex_t IT_0422 = 2*IT_0055;
    const ccomplex_t IT_0423 = 2*IT_0422;
    const ccomplex_t IT_0424 = 2*conj(IT_0422);
    const ccomplex_t IT_0425 = IT_0066 + 2*IT_0086 + 2*IT_0110 + -IT_0116;
    const ccomplex_t IT_0426 = 2*IT_0425;
    const ccomplex_t IT_0427 = 2*conj(IT_0425);
    const ccomplex_t IT_0428 = IT_0030*IT_0065;
    const ccomplex_t IT_0429 = IT_0066 + 2*IT_0086 + 2*IT_0110 + -IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0429;
    const ccomplex_t IT_0431 = 2*conj(IT_0429);
    const ccomplex_t IT_0432 = (conj(IT_0034) + conj(IT_0063) + conj(IT_0107) 
      + conj(IT_0114) + conj(IT_0186) + conj(IT_0188) + conj(IT_0189))*IT_0423 +
       (IT_0034 + IT_0063 + IT_0107 + IT_0114 + IT_0186 + IT_0188 + IT_0189)
      *IT_0424 + IT_0353*IT_0426 + IT_0357*IT_0427 + IT_0191*IT_0430 + IT_0216
      *IT_0431;
    const ccomplex_t IT_0433 = s_12*IT_0002;
    const ccomplex_t IT_0434 = s_13*s_23;
    const ccomplex_t IT_0435 = s_14*s_24;
    const ccomplex_t IT_0436 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0437 = (-2)*IT_0434 + -IT_0435 + 2*IT_0436;
    const ccomplex_t IT_0438 = IT_0433 + IT_0437;
    const ccomplex_t IT_0439 = IT_0005*IT_0438;
    const ccomplex_t IT_0440 = -IT_0223 + (-0.5)*IT_0439;
    const ccomplex_t IT_0441 = s_12 + IT_0440;
    const ccomplex_t IT_0442 = m_N_3*IT_0441;
    const ccomplex_t IT_0443 = (-2)*IT_0442;
    const ccomplex_t IT_0444 = IT_0106*(IT_0288 + -IT_0289);
    const ccomplex_t IT_0445 = IT_0113*(IT_0278 + -IT_0279);
    const ccomplex_t IT_0446 = conj(IT_0444) + conj(IT_0445);
    const ccomplex_t IT_0447 = IT_0444 + IT_0445;
    const ccomplex_t IT_0448 = -IT_0113*(IT_0279 + -IT_0313);
    const ccomplex_t IT_0449 = -IT_0106*(IT_0289 + -IT_0318);
    const ccomplex_t IT_0450 = conj(IT_0448) + conj(IT_0449);
    const ccomplex_t IT_0451 = IT_0448 + IT_0449;
    const ccomplex_t IT_0452 = conj(IT_0228) + IT_0229;
    const ccomplex_t IT_0453 = (-0.5)*IT_0283 + (-0.5)*IT_0298;
    const ccomplex_t IT_0454 = IT_0282 + IT_0453;
    const ccomplex_t IT_0455 = IT_0032*IT_0454;
    const ccomplex_t IT_0456 = 2*IT_0455;
    const ccomplex_t IT_0457 = 0.5*IT_0456;
    const ccomplex_t IT_0458 = 2*IT_0457;
    const ccomplex_t IT_0459 = IT_0228 + IT_0230;
    const ccomplex_t IT_0460 = 2*conj(IT_0457);
    const ccomplex_t IT_0461 = IT_0061*(IT_0322 + -IT_0369);
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = 2*IT_0462;
    const ccomplex_t IT_0464 = 2*conj(IT_0462);
    const ccomplex_t IT_0465 = IT_0061*(IT_0294 + -IT_0369);
    const ccomplex_t IT_0466 = -IT_0465;
    const ccomplex_t IT_0467 = 2*IT_0466;
    const ccomplex_t IT_0468 = 2*conj(IT_0466);
    const ccomplex_t IT_0469 = 2*IT_0231;
    const ccomplex_t IT_0470 = 2*conj(IT_0231);
    const ccomplex_t IT_0471 = IT_0348*IT_0446 + IT_0351*IT_0447 + IT_0355
      *IT_0450 + IT_0359*IT_0451 + IT_0452*IT_0458 + IT_0459*IT_0460 + IT_0229
      *IT_0463 + IT_0230*IT_0464 + conj(IT_0228)*IT_0467 + IT_0228*IT_0468 + 
      (conj(IT_0457) + conj(IT_0466))*IT_0469 + (IT_0457 + IT_0466)*IT_0470;
    const ccomplex_t IT_0472 = (-2)*IT_0434;
    const ccomplex_t IT_0473 = IT_0433 + IT_0472;
    const ccomplex_t IT_0474 = IT_0005*IT_0473;
    const ccomplex_t IT_0475 = (s_12 + -IT_0223)*(s_12 + (-0.5)*IT_0474);
    const ccomplex_t IT_0476 = (-2)*IT_0475;
    const ccomplex_t IT_0477 = conj(IT_0227) + conj(IT_0231);
    const ccomplex_t IT_0478 = IT_0227 + IT_0231;
    const ccomplex_t IT_0479 = 2*IT_0227;
    const ccomplex_t IT_0480 = 2*conj(IT_0227);
    const ccomplex_t IT_0481 = conj(IT_0228)*IT_0426 + IT_0228*IT_0427 + conj
      (IT_0228)*IT_0430 + IT_0228*IT_0431 + conj(IT_0429)*IT_0469 + IT_0429
      *IT_0470 + 2*IT_0423*(conj(IT_0228) + 1./2*IT_0477) + 2*IT_0424*(IT_0228 +
       1./2*IT_0478) + conj(IT_0425)*IT_0479 + IT_0425*IT_0480;
    const ccomplex_t IT_0482 = s_13*s_14*s_23;
    const ccomplex_t IT_0483 = s_23*s_34*IT_0052;
    const ccomplex_t IT_0484 = -IT_0413 + (-2)*IT_0482 + 2*IT_0483;
    const ccomplex_t IT_0485 = IT_0409 + IT_0484;
    const ccomplex_t IT_0486 = IT_0005*IT_0485;
    const ccomplex_t IT_0487 = (-0.5)*IT_0486;
    const ccomplex_t IT_0488 = -IT_0417;
    const ccomplex_t IT_0489 = IT_0419 + IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = 4*IT_0422;
    const ccomplex_t IT_0491 = conj(IT_0056)*IT_0489*IT_0490;
    const ccomplex_t IT_0492 = 4*conj(IT_0422);
    const ccomplex_t IT_0493 = IT_0056*IT_0489*IT_0492;
    const ccomplex_t IT_0494 = IT_0002*IT_0052;
    const ccomplex_t IT_0495 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0496 = 0.5*IT_0238 + (-0.5)*IT_0494 + -IT_0495;
    const ccomplex_t IT_0497 = IT_0236 + IT_0496;
    const ccomplex_t IT_0498 = IT_0005*IT_0497;
    const ccomplex_t IT_0499 = -IT_0239;
    const ccomplex_t IT_0500 = IT_0052 + IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = m_N_3*IT_0500;
    const ccomplex_t IT_0502 = conj(IT_0324) + conj(IT_0327);
    const ccomplex_t IT_0503 = (-4)*IT_0231;
    const ccomplex_t IT_0504 = (-4)*IT_0228;
    const ccomplex_t IT_0505 = IT_0503 + IT_0504;
    const ccomplex_t IT_0506 = IT_0324 + IT_0327;
    const ccomplex_t IT_0507 = (-4)*conj(IT_0231);
    const ccomplex_t IT_0508 = (-4)*conj(IT_0228);
    const ccomplex_t IT_0509 = IT_0507 + IT_0508;
    const ccomplex_t IT_0510 = conj(IT_0280) + conj(IT_0296) + conj(IT_0300);
    const ccomplex_t IT_0511 = (-4)*IT_0227;
    const ccomplex_t IT_0512 = IT_0504 + IT_0511;
    const ccomplex_t IT_0513 = IT_0280 + IT_0296 + IT_0300;
    const ccomplex_t IT_0514 = (-4)*conj(IT_0227);
    const ccomplex_t IT_0515 = IT_0508 + IT_0514;
    const ccomplex_t IT_0516 = 4*IT_0370;
    const ccomplex_t IT_0517 = 4*conj(IT_0370);
    const ccomplex_t IT_0518 = IT_0370*IT_0477 + conj(IT_0370)*IT_0478 + 0.25
      *IT_0502*IT_0505 + 0.25*IT_0506*IT_0509 + 0.25*IT_0510*IT_0512 + 0.25
      *IT_0513*IT_0515 + 0.5*conj(IT_0228)*IT_0516 + 0.5*IT_0228*IT_0517;
    const ccomplex_t IT_0519 = conj(IT_0181) + conj(IT_0185) + IT_0308;
    const ccomplex_t IT_0520 = IT_0181 + IT_0185 + IT_0330;
    const ccomplex_t IT_0521 = conj(IT_0056) + conj(IT_0087) + conj(IT_0111) +
       conj(IT_0117);
    const ccomplex_t IT_0522 = 4*IT_0425;
    const ccomplex_t IT_0523 = IT_0056 + IT_0087 + IT_0111 + IT_0117;
    const ccomplex_t IT_0524 = 4*conj(IT_0425);
    const ccomplex_t IT_0525 = 4*IT_0429;
    const ccomplex_t IT_0526 = 4*conj(IT_0429);
    const ccomplex_t IT_0527 = pow(s_14, 2);
    const ccomplex_t IT_0528 = IT_0005*IT_0527;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = IT_0052 + IT_0529;
    const ccomplex_t IT_0531 = IT_0052*IT_0530;
    const ccomplex_t IT_0532 = -IT_0531;
    const ccomplex_t IT_0533 = conj(IT_0422) + conj(IT_0425);
    const ccomplex_t IT_0534 = conj(IT_0429) + IT_0533;
    const ccomplex_t IT_0535 = 6*IT_0422;
    const ccomplex_t IT_0536 = IT_0422 + IT_0425;
    const ccomplex_t IT_0537 = IT_0429 + IT_0536;
    const ccomplex_t IT_0538 = 6*conj(IT_0422);
    const ccomplex_t IT_0539 = IT_0052 + IT_0499;
    const ccomplex_t IT_0540 = IT_0052*IT_0539;
    const ccomplex_t IT_0541 = -IT_0540;
    const ccomplex_t IT_0542 = 6*IT_0231;
    const ccomplex_t IT_0543 = 6*conj(IT_0231);
    const ccomplex_t IT_0544 = (-2)*IT_0435;
    const ccomplex_t IT_0545 = IT_0433 + IT_0544;
    const ccomplex_t IT_0546 = IT_0005*IT_0545;
    const ccomplex_t IT_0547 = s_12*IT_0269;
    const ccomplex_t IT_0548 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0549 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0550 = s_14*s_23*s_34;
    const ccomplex_t IT_0551 = (-2)*IT_0548 + (-2)*IT_0549 + 4*IT_0550;
    const ccomplex_t IT_0552 = IT_0547 + IT_0551;
    const ccomplex_t IT_0553 = IT_0000*IT_0552;
    const ccomplex_t IT_0554 = 8*IT_0370;
    const ccomplex_t IT_0555 = conj(IT_0370)*(s_12 + (-0.5)*IT_0474 + (-0.5)
      *IT_0546 + 0.25*IT_0553)*IT_0554;
    const ccomplex_t IT_0556 = 1.125*IT_0052;
    const ccomplex_t IT_0557 = conj(IT_0444) + conj(IT_0445) + conj(IT_0448) +
       conj(IT_0449);
    const ccomplex_t IT_0558 = 16*IT_0457;
    const ccomplex_t IT_0559 = IT_0444 + IT_0445 + IT_0448 + IT_0449;
    const ccomplex_t IT_0560 = 16*conj(IT_0457);
    const ccomplex_t IT_0561 = conj(IT_0280) + conj(IT_0296) + conj(IT_0300) +
       conj(IT_0324) + conj(IT_0327);
    const ccomplex_t IT_0562 = 16*IT_0370;
    const ccomplex_t IT_0563 = IT_0280 + IT_0296 + IT_0300 + IT_0324 + IT_0327;
    const ccomplex_t IT_0564 = 16*conj(IT_0370);
    const ccomplex_t IT_0565 = 16*IT_0466;
    const ccomplex_t IT_0566 = 16*conj(IT_0466);
    const ccomplex_t IT_0567 = (-16)*IT_0462;
    const ccomplex_t IT_0568 = (-16)*conj(IT_0462);
    const ccomplex_t IT_0569 = 16*IT_0462;
    const ccomplex_t IT_0570 = (-16)*IT_0448;
    const ccomplex_t IT_0571 = (-16)*IT_0449;
    const ccomplex_t IT_0572 = 16*conj(IT_0462);
    const ccomplex_t IT_0573 = (-16)*conj(IT_0448);
    const ccomplex_t IT_0574 = (-16)*conj(IT_0449);
    const ccomplex_t IT_0575 = (-16)*IT_0324;
    const ccomplex_t IT_0576 = (-16)*IT_0327;
    const ccomplex_t IT_0577 = (-16)*conj(IT_0324);
    const ccomplex_t IT_0578 = (-16)*conj(IT_0327);
    const ccomplex_t IT_0579 = (-32)*IT_0370*conj(IT_0370) + IT_0557*IT_0558 +
       IT_0559*IT_0560 + IT_0561*IT_0562 + IT_0563*IT_0564 + IT_0450*IT_0565 +
       IT_0451*IT_0566 + conj(IT_0466)*IT_0567 + conj(IT_0457)*((-16)*IT_0466 +
       IT_0567) + IT_0457*((-32)*conj(IT_0457) + (-16)*conj(IT_0466) + IT_0568) 
      + IT_0466*IT_0568 + IT_0446*(IT_0569 + IT_0570 + IT_0571) + IT_0447*
      (IT_0572 + IT_0573 + IT_0574) + IT_0510*(IT_0575 + IT_0576) + IT_0513*
      (IT_0577 + IT_0578);
    const ccomplex_t IT_0580 = IT_0005*IT_0433;
    const ccomplex_t IT_0581 = (-0.5)*IT_0580;
    const ccomplex_t IT_0582 = s_14*s_23;
    const ccomplex_t IT_0583 = s_13*s_24;
    const ccomplex_t IT_0584 = s_12*s_34;
    const ccomplex_t IT_0585 = -IT_0583 + -IT_0584;
    const ccomplex_t IT_0586 = IT_0582 + IT_0585;
    const ccomplex_t IT_0587 = IT_0001*IT_0586;
    const ccomplex_t IT_0588 = (-0.25)*IT_0587;
    const ccomplex_t IT_0589 = s_12 + IT_0581 + IT_0588;
    const ccomplex_t IT_0590 = 8*IT_0175;
    const ccomplex_t IT_0591 = (-8)*IT_0173;
    const ccomplex_t IT_0592 = IT_0590 + IT_0591;
    const ccomplex_t IT_0593 = 8*conj(IT_0175);
    const ccomplex_t IT_0594 = (-8)*conj(IT_0173);
    const ccomplex_t IT_0595 = IT_0593 + IT_0594;
    const ccomplex_t IT_0596 = 8*IT_0212;
    const ccomplex_t IT_0597 = (-8)*IT_0210;
    const ccomplex_t IT_0598 = IT_0596 + IT_0597;
    const ccomplex_t IT_0599 = 8*conj(IT_0212);
    const ccomplex_t IT_0600 = (-8)*conj(IT_0210);
    const ccomplex_t IT_0601 = IT_0599 + IT_0600;
    const ccomplex_t IT_0602 = 8*IT_0173;
    const ccomplex_t IT_0603 = (-8)*IT_0175;
    const ccomplex_t IT_0604 = 8*conj(IT_0173);
    const ccomplex_t IT_0605 = (-8)*conj(IT_0175);
    const ccomplex_t IT_0606 = 8*IT_0210;
    const ccomplex_t IT_0607 = (-8)*IT_0212;
    const ccomplex_t IT_0608 = 8*conj(IT_0210);
    const ccomplex_t IT_0609 = (-8)*conj(IT_0212);
    const ccomplex_t IT_0610 = IT_0450*IT_0592 + IT_0451*IT_0595 + IT_0446
      *IT_0598 + IT_0447*IT_0601 + conj(IT_0462)*(IT_0602 + IT_0603) + IT_0462*
      (IT_0604 + IT_0605) + conj(IT_0466)*(IT_0606 + IT_0607) + conj(IT_0457)*
      (IT_0602 + IT_0603 + IT_0606 + IT_0607) + IT_0466*(IT_0608 + IT_0609) +
       IT_0457*(IT_0604 + IT_0605 + IT_0608 + IT_0609);
    const ccomplex_t IT_0611 = IT_0582 + IT_0584;
    const ccomplex_t IT_0612 = -IT_0583;
    const ccomplex_t IT_0613 = IT_0611 + IT_0612;
    const ccomplex_t IT_0614 = IT_0001*IT_0613;
    const ccomplex_t IT_0615 = 0.25*IT_0614;
    const ccomplex_t IT_0616 = s_12 + IT_0581 + IT_0615;
    const ccomplex_t IT_0617 = (-0.5)*IT_0546;
    const ccomplex_t IT_0618 = (-0.5)*IT_0474;
    const ccomplex_t IT_0619 = 0.25*IT_0553;
    const ccomplex_t IT_0620 = 8*IT_0280;
    const ccomplex_t IT_0621 = 8*conj(IT_0280);
    const ccomplex_t IT_0622 = 8*IT_0324;
    const ccomplex_t IT_0623 = 8*conj(IT_0324);
    const ccomplex_t IT_0624 = 8*IT_0327;
    const ccomplex_t IT_0625 = 8*IT_0296;
    const ccomplex_t IT_0626 = 8*conj(IT_0296);
    const ccomplex_t IT_0627 = (-8)*IT_0370;
    const ccomplex_t IT_0628 = (-8)*conj(IT_0370);
    const ccomplex_t IT_0629 = IT_0510*IT_0620 + IT_0296*IT_0621 + IT_0502
      *IT_0622 + IT_0327*IT_0623 + conj(IT_0327)*IT_0624 + (conj(IT_0296) + conj
      (IT_0300))*IT_0625 + 8*IT_0300*(conj(IT_0300) + 0.125*IT_0621 + 0.125
      *IT_0626) + IT_0561*IT_0627 + IT_0563*IT_0628;
    const ccomplex_t IT_0630 = s_13*s_24*s_34;
    const ccomplex_t IT_0631 = pow(s_34, 2);
    const ccomplex_t IT_0632 = s_12*IT_0631;
    const ccomplex_t IT_0633 = (-2)*IT_0548 + (-2)*IT_0549 + 2*IT_0550 + 2
      *IT_0630 + (-2)*IT_0632;
    const ccomplex_t IT_0634 = IT_0547 + IT_0633;
    const ccomplex_t IT_0635 = IT_0000*IT_0634;
    const ccomplex_t IT_0636 = (-0.25)*IT_0635;
    const ccomplex_t IT_0637 = -IT_0223;
    const ccomplex_t IT_0638 = -IT_0222;
    const ccomplex_t IT_0639 = 8*IT_0457;
    const ccomplex_t IT_0640 = 8*IT_0466;
    const ccomplex_t IT_0641 = IT_0639 + IT_0640;
    const ccomplex_t IT_0642 = 8*conj(IT_0448);
    const ccomplex_t IT_0643 = (-8)*conj(IT_0462);
    const ccomplex_t IT_0644 = IT_0642 + IT_0643;
    const ccomplex_t IT_0645 = 8*conj(IT_0457);
    const ccomplex_t IT_0646 = 8*conj(IT_0444);
    const ccomplex_t IT_0647 = 8*IT_0462;
    const ccomplex_t IT_0648 = conj(IT_0370)*IT_0647;
    const ccomplex_t IT_0649 = 8*conj(IT_0462);
    const ccomplex_t IT_0650 = IT_0370*IT_0649;
    const ccomplex_t IT_0651 = conj(IT_0457)*IT_0554;
    const ccomplex_t IT_0652 = 8*conj(IT_0370);
    const ccomplex_t IT_0653 = IT_0457*IT_0652;
    const ccomplex_t IT_0654 = (-8)*IT_0457;
    const ccomplex_t IT_0655 = conj(IT_0280)*IT_0654;
    const ccomplex_t IT_0656 = conj(IT_0324)*IT_0654;
    const ccomplex_t IT_0657 = conj(IT_0327)*IT_0654;
    const ccomplex_t IT_0658 = conj(IT_0296)*IT_0654;
    const ccomplex_t IT_0659 = conj(IT_0300)*IT_0654;
    const ccomplex_t IT_0660 = (-8)*conj(IT_0457);
    const ccomplex_t IT_0661 = IT_0280*IT_0660;
    const ccomplex_t IT_0662 = IT_0324*IT_0660;
    const ccomplex_t IT_0663 = IT_0327*IT_0660;
    const ccomplex_t IT_0664 = IT_0296*IT_0660;
    const ccomplex_t IT_0665 = IT_0300*IT_0660;
    const ccomplex_t IT_0666 = 8*IT_0448;
    const ccomplex_t IT_0667 = 8*IT_0449;
    const ccomplex_t IT_0668 = IT_0666 + IT_0667;
    const ccomplex_t IT_0669 = (-8)*IT_0462;
    const ccomplex_t IT_0670 = 8*conj(IT_0449);
    const ccomplex_t IT_0671 = conj(IT_0444)*IT_0627;
    const ccomplex_t IT_0672 = conj(IT_0449)*IT_0627;
    const ccomplex_t IT_0673 = conj(IT_0445)*IT_0627;
    const ccomplex_t IT_0674 = IT_0444*IT_0628;
    const ccomplex_t IT_0675 = IT_0449*IT_0628;
    const ccomplex_t IT_0676 = IT_0445*IT_0628;
    const ccomplex_t IT_0677 = 8*IT_0444;
    const ccomplex_t IT_0678 = (-8)*IT_0466;
    const ccomplex_t IT_0679 = 8*IT_0445;
    const ccomplex_t IT_0680 = (-8)*conj(IT_0466);
    const ccomplex_t IT_0681 = conj(IT_0444)*IT_0620 + IT_0447*IT_0621 + conj
      (IT_0449)*IT_0622 + IT_0449*IT_0623 + 8*conj(IT_0445)*(IT_0300 + 0.125
      *IT_0620 + 0.125*IT_0625) + IT_0445*IT_0626 + conj(IT_0448)*IT_0627 +
       IT_0448*IT_0628 + conj(IT_0370)*IT_0641 + IT_0324*IT_0644 + 8*IT_0370*
      (conj(IT_0466) + 0.125*IT_0645) + (IT_0296 + IT_0300)*IT_0646 + IT_0648 +
       IT_0650 + IT_0651 + IT_0653 + IT_0655 + IT_0656 + IT_0657 + IT_0658 +
       IT_0659 + IT_0661 + IT_0662 + IT_0663 + IT_0664 + IT_0665 + conj(IT_0324)
      *(IT_0666 + IT_0669) + conj(IT_0327)*(IT_0668 + IT_0669) + IT_0327*
      (IT_0642 + IT_0643 + IT_0670) + IT_0671 + IT_0672 + IT_0673 + IT_0674 +
       IT_0675 + IT_0676 + conj(IT_0280)*IT_0678 + conj(IT_0296)*(IT_0677 +
       IT_0678) + conj(IT_0300)*(IT_0677 + IT_0678 + IT_0679) + IT_0513*IT_0680;
    const ccomplex_t IT_0682 = (-2)*IT_0631;
    const ccomplex_t IT_0683 = IT_0269 + IT_0682;
    const ccomplex_t IT_0684 = IT_0000*IT_0683;
    const ccomplex_t IT_0685 = 0.25*IT_0684;
    const ccomplex_t IT_0686 = 1 + IT_0685;
    const ccomplex_t IT_0687 = IT_0052*IT_0686;
    const ccomplex_t IT_0688 = 8*conj(IT_0327);
    const ccomplex_t IT_0689 = conj(IT_0466)*IT_0554 + IT_0446*(IT_0622 +
       IT_0624) + conj(IT_0448)*IT_0627 + IT_0448*IT_0628 + conj(IT_0370)
      *IT_0639 + IT_0370*IT_0645 + IT_0648 + IT_0650 + IT_0651 + IT_0466*IT_0652
       + IT_0653 + IT_0655 + IT_0656 + IT_0657 + IT_0658 + IT_0659 + IT_0661 +
       IT_0662 + IT_0663 + IT_0664 + IT_0665 + IT_0510*(IT_0668 + IT_0669) +
       IT_0513*(IT_0642 + IT_0643 + IT_0670) + IT_0671 + IT_0672 + IT_0673 +
       IT_0674 + IT_0675 + IT_0676 + IT_0502*IT_0678 + IT_0506*IT_0680 + IT_0447
      *(IT_0623 + IT_0688);
    const ccomplex_t IT_0690 = (-2)*IT_0548 + (-2)*IT_0549 + 4*IT_0630;
    const ccomplex_t IT_0691 = IT_0547 + IT_0690;
    const ccomplex_t IT_0692 = IT_0000*IT_0691;
    const ccomplex_t IT_0693 = 0.25*IT_0692;
    const ccomplex_t IT_0694 = conj(IT_0466)*IT_0641 + IT_0448*IT_0643 +
       IT_0449*IT_0644 + IT_0466*IT_0645 + conj(IT_0462)*IT_0647 + conj(IT_0457)
      *(IT_0558 + IT_0647) + IT_0457*IT_0649 + IT_0559*IT_0660 + conj(IT_0448)*
      (IT_0654 + IT_0666 + IT_0669) + conj(IT_0449)*(IT_0654 + IT_0668 + IT_0669
      ) + conj(IT_0444)*(IT_0654 + IT_0677 + IT_0678) + conj(IT_0445)*(IT_0654 +
       IT_0677 + IT_0678 + IT_0679) + IT_0444*IT_0680 + IT_0445*(IT_0646 +
       IT_0680);
    const ccomplex_t IT_0695 = pow(m_Z, 6);
    const ccomplex_t IT_0696 = s_12*IT_0695;
    const ccomplex_t IT_0697 = s_13*s_23*IT_0269;
    const ccomplex_t IT_0698 = (-2)*IT_0697;
    const ccomplex_t IT_0699 = IT_0696 + IT_0698;
    const ccomplex_t IT_0700 = IT_0000*IT_0699;
    const ccomplex_t IT_0701 = 0.25*IT_0700;
    const ccomplex_t IT_0702 = IT_0433 + IT_0701;
    const ccomplex_t IT_0703 = (-8)*IT_0188;
    const ccomplex_t IT_0704 = (-8)*IT_0186;
    const ccomplex_t IT_0705 = (-8)*conj(IT_0188);
    const ccomplex_t IT_0706 = (-8)*conj(IT_0186);
    const ccomplex_t IT_0707 = (-8)*IT_0056;
    const ccomplex_t IT_0708 = (-8)*conj(IT_0056);
    const ccomplex_t IT_0709 = (-8)*IT_0034;
    const ccomplex_t IT_0710 = (-8)*IT_0063;
    const ccomplex_t IT_0711 = (-8)*conj(IT_0034);
    const ccomplex_t IT_0712 = (-8)*conj(IT_0063);
    const ccomplex_t IT_0713 = conj(IT_0107) + conj(IT_0114);
    const ccomplex_t IT_0714 = (-8)*IT_0087;
    const ccomplex_t IT_0715 = (-8)*conj(IT_0087);
    const ccomplex_t IT_0716 = conj(IT_0111) + conj(IT_0117);
    const ccomplex_t IT_0717 = (-8)*IT_0107;
    const ccomplex_t IT_0718 = IT_0111 + IT_0117;
    const ccomplex_t IT_0719 = (-8)*conj(IT_0107);
    const ccomplex_t IT_0720 = IT_0276*(IT_0189 + (-0.125)*IT_0703 + (-0.125)
      *IT_0704) + IT_0303*(conj(IT_0189) + (-0.125)*IT_0705 + (-0.125)*IT_0706) 
      + (-0.125)*(conj(IT_0034) + conj(IT_0107) + conj(IT_0114) + conj(IT_0188) 
      + conj(IT_0189))*IT_0707 + (-0.125)*(IT_0034 + IT_0107 + IT_0114 + IT_0188
       + IT_0189)*IT_0708 + (-0.125)*conj(IT_0056)*(IT_0704 + IT_0710) + (-0.125
      )*conj(IT_0087)*(IT_0709 + IT_0710) + (-0.125)*IT_0056*(IT_0706 + IT_0712)
       + (-0.125)*IT_0087*(IT_0711 + IT_0712) + (-0.125)*IT_0713*IT_0714 + (
      -0.125)*IT_0189*IT_0715 + IT_0716*(IT_0114 + (-0.125)*IT_0709 + (-0.125)
      *IT_0710 + (-0.125)*IT_0717) + IT_0718*(conj(IT_0114) + (-0.125)*IT_0711 +
       (-0.125)*IT_0712 + (-0.125)*IT_0719);
    const ccomplex_t IT_0721 = s_12 + IT_0637;
    const ccomplex_t IT_0722 = IT_0052*IT_0721;
    const ccomplex_t IT_0723 = (-6)*IT_0422;
    const ccomplex_t IT_0724 = (-6)*conj(IT_0422);
    const ccomplex_t IT_0725 = (-6)*IT_0231;
    const ccomplex_t IT_0726 = (-6)*conj(IT_0231);
    const ccomplex_t IT_0727 = IT_0229*IT_0429 + IT_0230*conj(IT_0429) + (
      -0.166666666666667)*IT_0452*IT_0723 + (-0.166666666666667)*IT_0459*IT_0724
       + (-0.166666666666667)*conj(IT_0422)*IT_0725 + conj(IT_0425)*(IT_0228 + (
      -0.166666666666667)*IT_0725) + (-0.166666666666667)*IT_0422*IT_0726 +
       IT_0425*(conj(IT_0228) + (-0.166666666666667)*IT_0726);
    const ccomplex_t IT_0728 = IT_0000*IT_0631;
    const ccomplex_t IT_0729 = 0.5*IT_0728;
    const ccomplex_t IT_0730 = 1 + IT_0729;
    const ccomplex_t IT_0731 = s_12*IT_0730;
    const ccomplex_t IT_0732 = 4*IT_0175;
    const ccomplex_t IT_0733 = (-4)*IT_0173;
    const ccomplex_t IT_0734 = IT_0732 + IT_0733;
    const ccomplex_t IT_0735 = 4*IT_0212;
    const ccomplex_t IT_0736 = (-4)*IT_0210;
    const ccomplex_t IT_0737 = IT_0735 + IT_0736;
    const ccomplex_t IT_0738 = 4*IT_0173;
    const ccomplex_t IT_0739 = 4*IT_0210;
    const ccomplex_t IT_0740 = (-4)*conj(IT_0173);
    const ccomplex_t IT_0741 = (-4)*conj(IT_0210);
    const ccomplex_t IT_0742 = IT_0052*IT_0730;
    const ccomplex_t IT_0743 = 4*conj(IT_0175);
    const ccomplex_t IT_0744 = (-4)*IT_0175;
    const ccomplex_t IT_0745 = 4*conj(IT_0173);
    const ccomplex_t IT_0746 = (-4)*conj(IT_0175);
    const ccomplex_t IT_0747 = 4*conj(IT_0212);
    const ccomplex_t IT_0748 = (-4)*IT_0212;
    const ccomplex_t IT_0749 = 4*conj(IT_0210);
    const ccomplex_t IT_0750 = (-4)*conj(IT_0212);
    const ccomplex_t IT_0751 = (conj(IT_0210) + conj(IT_0466))*IT_0734 + (conj
      (IT_0370) + conj(IT_0457))*(IT_0732 + IT_0733 + IT_0735 + IT_0736) + conj
      (IT_0462)*IT_0737 + (IT_0210 + IT_0466)*(IT_0740 + IT_0743) + (conj
      (IT_0212) + conj(IT_0280) + conj(IT_0296) + conj(IT_0300) + conj(IT_0444) 
      + conj(IT_0445))*(IT_0738 + IT_0744) + (IT_0212 + IT_0280 + IT_0296 +
       IT_0300 + IT_0444 + IT_0445)*(IT_0745 + IT_0746) + IT_0462*(IT_0741 +
       IT_0747) + (IT_0370 + IT_0457)*(IT_0740 + IT_0741 + IT_0743 + IT_0747) + 
      (conj(IT_0324) + conj(IT_0327) + conj(IT_0448) + conj(IT_0449))*(IT_0739 +
       IT_0748) + (IT_0324 + IT_0327 + IT_0448 + IT_0449)*(IT_0749 + IT_0750);
    const ccomplex_t IT_0752 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0753 = (-0.5)*IT_0494 + 0.5*IT_0527 + -IT_0752;
    const ccomplex_t IT_0754 = IT_0410 + IT_0753;
    const ccomplex_t IT_0755 = IT_0005*IT_0754;
    const ccomplex_t IT_0756 = IT_0052 + IT_0529 + IT_0755;
    const ccomplex_t IT_0757 = m_N_3*IT_0756;
    const ccomplex_t IT_0758 = (-4)*IT_0422;
    const ccomplex_t IT_0759 = (-4)*conj(IT_0422);
    const ccomplex_t IT_0760 = (-4)*IT_0457;
    const ccomplex_t IT_0761 = (-4)*conj(IT_0457);
    const ccomplex_t IT_0762 = (-4)*IT_0462;
    const ccomplex_t IT_0763 = (-4)*conj(IT_0462);
    const ccomplex_t IT_0764 = IT_0450*IT_0522 + IT_0451*IT_0524 + IT_0446
      *IT_0525 + IT_0447*IT_0526 + IT_0490*IT_0557 + IT_0492*IT_0559 + conj
      (IT_0457)*IT_0758 + (-4)*conj(IT_0466)*(IT_0429 + (-0.25)*IT_0758) +
       IT_0457*IT_0759 + (-4)*IT_0466*(conj(IT_0429) + (-0.25)*IT_0759) +
       IT_0534*IT_0760 + IT_0537*IT_0761 + IT_0533*IT_0762 + IT_0536*IT_0763;
    const ccomplex_t IT_0765 = s_14*IT_0236;
    const ccomplex_t IT_0766 = s_12*s_23*s_34;
    const ccomplex_t IT_0767 = -IT_0336 + 2*IT_0765 + (-2)*IT_0766;
    const ccomplex_t IT_0768 = IT_0334 + IT_0767;
    const ccomplex_t IT_0769 = IT_0005*IT_0768;
    const ccomplex_t IT_0770 = 0.5*IT_0769;
    const ccomplex_t IT_0771 = -IT_0770;
    const ccomplex_t IT_0772 = -IT_0341;
    const ccomplex_t IT_0773 = IT_0343 + IT_0771 + IT_0772;
    const ccomplex_t IT_0774 = 0.5*IT_0700;
    const ccomplex_t IT_0775 = (-4)*IT_0435;
    const ccomplex_t IT_0776 = IT_0433 + IT_0774 + IT_0775;
    const ccomplex_t IT_0777 = IT_0005*IT_0552;
    const ccomplex_t IT_0778 = -IT_0777;
    const ccomplex_t IT_0779 = IT_0776 + IT_0778;
    const ccomplex_t IT_0780 = (-4)*IT_0056;
    const ccomplex_t IT_0781 = conj(IT_0056)*IT_0779*IT_0780;
    const ccomplex_t IT_0782 = IT_0005*IT_0691;
    const ccomplex_t IT_0783 = -IT_0782;
    const ccomplex_t IT_0784 = IT_0776 + IT_0783;
    const ccomplex_t IT_0785 = (-4)*conj(IT_0034);
    const ccomplex_t IT_0786 = (-4)*IT_0186;
    const ccomplex_t IT_0787 = (-4)*IT_0188;
    const ccomplex_t IT_0788 = (-4)*conj(IT_0186);
    const ccomplex_t IT_0789 = (-4)*IT_0063;
    const ccomplex_t IT_0790 = (-4)*IT_0034;
    const ccomplex_t IT_0791 = (-4)*conj(IT_0063);
    const ccomplex_t IT_0792 = (-4)*IT_0107;
    const ccomplex_t IT_0793 = IT_0189*IT_0785 + conj(IT_0186)*IT_0786 + (-4)
      *IT_0190*(IT_0189 + (-0.25)*IT_0786 + (-0.25)*IT_0787) + IT_0215*IT_0788 +
       conj(IT_0063)*IT_0789 + conj(IT_0034)*(IT_0789 + IT_0790) + (IT_0034 +
       IT_0107 + IT_0114)*IT_0791 + (-4)*IT_0713*(IT_0114 + (-0.25)*IT_0789 + (
      -0.25)*IT_0790 + (-0.25)*IT_0792);
    const ccomplex_t IT_0794 = (-4)*conj(IT_0056);
    const ccomplex_t IT_0795 = (-4)*IT_0087;
    const ccomplex_t IT_0796 = (-4)*conj(IT_0087);
    const ccomplex_t IT_0797 = (-4)*IT_0111;
    const ccomplex_t IT_0798 = (-4)*IT_0181;
    const ccomplex_t IT_0799 = IT_0779*(IT_0519*IT_0780 + IT_0520*IT_0794 +
       conj(IT_0087)*IT_0795 + IT_0718*IT_0796 + (-4)*IT_0716*(IT_0117 + (-0.25)
      *IT_0795 + (-0.25)*IT_0797) + (-4)*IT_0276*(IT_0185 + (-0.25)*IT_0798));
    const ccomplex_t IT_0800 = (s_12 + (-0.5)*IT_0474)*(IT_0052 + -IT_0528);
    const ccomplex_t IT_0801 = 2*IT_0800;
    const ccomplex_t IT_0802 = (-2)*IT_0422;
    const ccomplex_t IT_0803 = conj(IT_0422)*IT_0801*IT_0802;
    const ccomplex_t IT_0804 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0805 = (-2)*IT_0434 + -IT_0435 + 2*IT_0804;
    const ccomplex_t IT_0806 = IT_0433 + IT_0805;
    const ccomplex_t IT_0807 = IT_0005*IT_0806;
    const ccomplex_t IT_0808 = -IT_0223 + (-0.5)*IT_0807;
    const ccomplex_t IT_0809 = s_12 + IT_0808;
    const ccomplex_t IT_0810 = m_N_3*IT_0809;
    const ccomplex_t IT_0811 = (-2)*IT_0810;
    const ccomplex_t IT_0812 = conj(IT_0370)*IT_0802*IT_0811;
    const ccomplex_t IT_0813 = (-2)*conj(IT_0422);
    const ccomplex_t IT_0814 = IT_0370*IT_0811*IT_0813;
    const ccomplex_t IT_0815 = IT_0005*IT_0360;
    const ccomplex_t IT_0816 = (-2)*IT_0248 + (-0.5)*IT_0815;
    const ccomplex_t IT_0817 = IT_0246 + IT_0816;
    const ccomplex_t IT_0818 = IT_0402*IT_0817;
    const ccomplex_t IT_0819 = (-2)*IT_0818;
    const ccomplex_t IT_0820 = (-2)*IT_0457;
    const ccomplex_t IT_0821 = conj(IT_0056)*IT_0820;
    const ccomplex_t IT_0822 = IT_0819*IT_0821;
    const ccomplex_t IT_0823 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0824 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0825 = (-2)*IT_0824;
    const ccomplex_t IT_0826 = IT_0270 + IT_0825;
    const ccomplex_t IT_0827 = IT_0000*IT_0826;
    const ccomplex_t IT_0828 = (-0.666666666666667)*IT_0823 + (
      -0.166666666666667)*IT_0827;
    const ccomplex_t IT_0829 = s_14 + IT_0828;
    const ccomplex_t IT_0830 = m_N_3*IT_0829;
    const ccomplex_t IT_0831 = (-6)*IT_0830;
    const ccomplex_t IT_0832 = IT_0821*IT_0831;
    const ccomplex_t IT_0833 = (-2)*conj(IT_0457);
    const ccomplex_t IT_0834 = IT_0056*IT_0833;
    const ccomplex_t IT_0835 = IT_0819*IT_0834;
    const ccomplex_t IT_0836 = IT_0831*IT_0834;
    const ccomplex_t IT_0837 = (-2)*IT_0370;
    const ccomplex_t IT_0838 = IT_0425 + IT_0429;
    const ccomplex_t IT_0839 = (-2)*conj(IT_0370);
    const ccomplex_t IT_0840 = (-0.666666666666667)*IT_0011 +
       0.166666666666667*IT_0827;
    const ccomplex_t IT_0841 = s_14 + IT_0840;
    const ccomplex_t IT_0842 = m_N_3*IT_0841;
    const ccomplex_t IT_0843 = 6*IT_0842;
    const ccomplex_t IT_0844 = 2*IT_0448;
    const ccomplex_t IT_0845 = 2*IT_0449;
    const ccomplex_t IT_0846 = IT_0844 + IT_0845;
    const ccomplex_t IT_0847 = (-2)*IT_0462;
    const ccomplex_t IT_0848 = IT_0846 + IT_0847;
    const ccomplex_t IT_0849 = 2*conj(IT_0448);
    const ccomplex_t IT_0850 = 2*conj(IT_0449);
    const ccomplex_t IT_0851 = IT_0849 + IT_0850;
    const ccomplex_t IT_0852 = (-2)*conj(IT_0462);
    const ccomplex_t IT_0853 = IT_0851 + IT_0852;
    const ccomplex_t IT_0854 = 2*IT_0444;
    const ccomplex_t IT_0855 = 2*IT_0445;
    const ccomplex_t IT_0856 = IT_0854 + IT_0855;
    const ccomplex_t IT_0857 = (-2)*IT_0466;
    const ccomplex_t IT_0858 = IT_0856 + IT_0857;
    const ccomplex_t IT_0859 = 2*conj(IT_0444);
    const ccomplex_t IT_0860 = 2*conj(IT_0445);
    const ccomplex_t IT_0861 = IT_0859 + IT_0860;
    const ccomplex_t IT_0862 = (-2)*conj(IT_0466);
    const ccomplex_t IT_0863 = IT_0861 + IT_0862;
    const ccomplex_t IT_0864 = conj(IT_0034)*IT_0820;
    const ccomplex_t IT_0865 = conj(IT_0189)*IT_0820;
    const ccomplex_t IT_0866 = conj(IT_0186)*IT_0820;
    const ccomplex_t IT_0867 = conj(IT_0063)*IT_0820;
    const ccomplex_t IT_0868 = conj(IT_0107)*IT_0820;
    const ccomplex_t IT_0869 = conj(IT_0114)*IT_0820;
    const ccomplex_t IT_0870 = IT_0034*IT_0833;
    const ccomplex_t IT_0871 = IT_0189*IT_0833;
    const ccomplex_t IT_0872 = IT_0186*IT_0833;
    const ccomplex_t IT_0873 = IT_0063*IT_0833;
    const ccomplex_t IT_0874 = IT_0107*IT_0833;
    const ccomplex_t IT_0875 = IT_0114*IT_0833;
    const ccomplex_t IT_0876 = IT_0353*IT_0848 + IT_0357*IT_0853 + conj
      (IT_0188)*(IT_0820 + IT_0856 + IT_0857) + (conj(IT_0186) + conj(IT_0189))
      *IT_0858 + IT_0188*(IT_0833 + IT_0861 + IT_0862) + (IT_0186 + IT_0189)
      *IT_0863 + IT_0864 + IT_0865 + IT_0866 + IT_0867 + IT_0868 + IT_0869 +
       IT_0870 + IT_0871 + IT_0872 + IT_0873 + IT_0874 + IT_0875;
    const ccomplex_t IT_0877 = IT_0000*IT_0360;
    const ccomplex_t IT_0878 = 0.666666666666667*IT_0247 + (-0.166666666666667
      )*IT_0877;
    const ccomplex_t IT_0879 = s_24 + IT_0878;
    const ccomplex_t IT_0880 = m_N_3*IT_0879;
    const ccomplex_t IT_0881 = (-6)*IT_0880;
    const ccomplex_t IT_0882 = (conj(IT_0186)*IT_0847 + conj(IT_0188)*(IT_0820
       + IT_0847) + IT_0186*IT_0852 + IT_0188*(IT_0833 + IT_0852) + IT_0352
      *IT_0857 + IT_0356*IT_0862 + IT_0864 + IT_0866 + IT_0867 + IT_0870 +
       IT_0872 + IT_0873)*IT_0881;
    const ccomplex_t IT_0883 = (-10)*s_24;
    const ccomplex_t IT_0884 = IT_0877 + IT_0883;
    const ccomplex_t IT_0885 = m_N_3*IT_0884;
    const ccomplex_t IT_0886 = IT_0191*IT_0846 + conj(IT_0189)*IT_0847 +
       IT_0216*IT_0851 + IT_0189*IT_0852 + IT_0353*IT_0856 + IT_0713*IT_0857 +
       IT_0357*IT_0861 + IT_0189*IT_0862 + IT_0865 + IT_0868 + IT_0869 + IT_0871
       + IT_0874 + IT_0875;
    const ccomplex_t IT_0887 = conj(IT_0056)*IT_0844;
    const ccomplex_t IT_0888 = IT_0056*IT_0849;
    const ccomplex_t IT_0889 = conj(IT_0056)*IT_0854;
    const ccomplex_t IT_0890 = IT_0056*IT_0859;
    const ccomplex_t IT_0891 = conj(IT_0056)*IT_0845;
    const ccomplex_t IT_0892 = IT_0056*IT_0850;
    const ccomplex_t IT_0893 = conj(IT_0056)*IT_0855;
    const ccomplex_t IT_0894 = IT_0056*IT_0860;
    const ccomplex_t IT_0895 = conj(IT_0087)*IT_0820;
    const ccomplex_t IT_0896 = conj(IT_0111)*IT_0820;
    const ccomplex_t IT_0897 = conj(IT_0117)*IT_0820;
    const ccomplex_t IT_0898 = conj(IT_0181)*IT_0820;
    const ccomplex_t IT_0899 = conj(IT_0185)*IT_0820;
    const ccomplex_t IT_0900 = IT_0087*IT_0833;
    const ccomplex_t IT_0901 = IT_0111*IT_0833;
    const ccomplex_t IT_0902 = IT_0117*IT_0833;
    const ccomplex_t IT_0903 = IT_0181*IT_0833;
    const ccomplex_t IT_0904 = IT_0185*IT_0833;
    const ccomplex_t IT_0905 = conj(IT_0056)*IT_0847;
    const ccomplex_t IT_0906 = IT_0056*IT_0852;
    const ccomplex_t IT_0907 = conj(IT_0056)*IT_0857;
    const ccomplex_t IT_0908 = IT_0056*IT_0862;
    const ccomplex_t IT_0909 = IT_0276*IT_0848 + IT_0303*IT_0853 + IT_0308
      *IT_0858 + IT_0330*IT_0863 + IT_0887 + IT_0888 + IT_0889 + IT_0890 +
       IT_0891 + IT_0892 + IT_0893 + IT_0894 + IT_0895 + IT_0896 + IT_0897 +
       IT_0898 + IT_0899 + IT_0900 + IT_0901 + IT_0902 + IT_0903 + IT_0904 +
       IT_0905 + IT_0906 + IT_0907 + IT_0908;
    const ccomplex_t IT_0910 = IT_0308*IT_0848 + IT_0330*IT_0853 + IT_0276
      *IT_0858 + IT_0303*IT_0863 + IT_0887 + IT_0888 + IT_0889 + IT_0890 +
       IT_0891 + IT_0892 + IT_0893 + IT_0894 + IT_0895 + IT_0896 + IT_0897 +
       IT_0898 + IT_0899 + IT_0900 + IT_0901 + IT_0902 + IT_0903 + IT_0904 +
       IT_0905 + IT_0906 + IT_0907 + IT_0908;
    const ccomplex_t IT_0911 = conj(IT_0087)*IT_0371;
    const ccomplex_t IT_0912 = conj(IT_0111)*IT_0371;
    const ccomplex_t IT_0913 = conj(IT_0117)*IT_0371;
    const ccomplex_t IT_0914 = conj(IT_0181)*IT_0371;
    const ccomplex_t IT_0915 = conj(IT_0185)*IT_0371;
    const ccomplex_t IT_0916 = IT_0087*IT_0378;
    const ccomplex_t IT_0917 = IT_0111*IT_0378;
    const ccomplex_t IT_0918 = IT_0117*IT_0378;
    const ccomplex_t IT_0919 = IT_0181*IT_0378;
    const ccomplex_t IT_0920 = IT_0185*IT_0378;
    const ccomplex_t IT_0921 = conj(IT_0056)*IT_0281;
    const ccomplex_t IT_0922 = IT_0056*IT_0304;
    const ccomplex_t IT_0923 = conj(IT_0056)*IT_0325;
    const ccomplex_t IT_0924 = IT_0056*IT_0331;
    const ccomplex_t IT_0925 = conj(IT_0056)*IT_0328;
    const ccomplex_t IT_0926 = IT_0056*IT_0332;
    const ccomplex_t IT_0927 = conj(IT_0056)*IT_0297;
    const ccomplex_t IT_0928 = IT_0056*IT_0305;
    const ccomplex_t IT_0929 = conj(IT_0056)*IT_0301;
    const ccomplex_t IT_0930 = IT_0056*IT_0306;
    const ccomplex_t IT_0931 = IT_0376*(IT_0911 + IT_0912 + IT_0913 + IT_0914 
      + IT_0915 + IT_0916 + IT_0917 + IT_0918 + IT_0919 + IT_0920 + IT_0921 +
       IT_0922 + IT_0923 + IT_0924 + IT_0925 + IT_0926 + IT_0927 + IT_0928 +
       IT_0929 + IT_0930);
    const ccomplex_t IT_0932 = IT_0302*IT_0308 + IT_0276*IT_0329 + IT_0307
      *IT_0330 + IT_0303*IT_0333 + IT_0911 + IT_0912 + IT_0913 + IT_0914 +
       IT_0915 + IT_0916 + IT_0917 + IT_0918 + IT_0919 + IT_0920 + IT_0921 +
       IT_0922 + IT_0923 + IT_0924 + IT_0925 + IT_0926 + IT_0927 + IT_0928 +
       IT_0929 + IT_0930;
    const ccomplex_t IT_0933 = (-2)*IT_0425;
    const ccomplex_t IT_0934 = (-2)*IT_0429;
    const ccomplex_t IT_0935 = -IT_0494;
    const ccomplex_t IT_0936 = IT_0527 + IT_0935;
    const ccomplex_t IT_0937 = IT_0005*IT_0936;
    const ccomplex_t IT_0938 = -IT_0052;
    const ccomplex_t IT_0939 = IT_0528 + IT_0938;
    const ccomplex_t IT_0940 = (-4)*IT_0939;
    const ccomplex_t IT_0941 = IT_0937 + IT_0940;
    const ccomplex_t IT_0942 = m_N_3*IT_0941;
    const ccomplex_t IT_0943 = (-2)*conj(IT_0425);
    const ccomplex_t IT_0944 = (-2)*conj(IT_0429);
    const ccomplex_t IT_0945 = IT_0942*(conj(IT_0370)*IT_0423 + IT_0370
      *IT_0424 + IT_0371*IT_0534 + IT_0378*IT_0537 + IT_0561*IT_0802 + IT_0563
      *IT_0813 + IT_0502*IT_0933 + IT_0510*IT_0934 + IT_0506*IT_0943 + IT_0513
      *IT_0944);
    const ccomplex_t IT_0946 = -IT_0435;
    const ccomplex_t IT_0947 = IT_0433 + IT_0946;
    const ccomplex_t IT_0948 = IT_0005*IT_0947;
    const ccomplex_t IT_0949 = -IT_0223 + (-0.25)*IT_0948;
    const ccomplex_t IT_0950 = s_12 + IT_0949;
    const ccomplex_t IT_0951 = m_N_3*IT_0950;
    const ccomplex_t IT_0952 = (-4)*IT_0951;
    const ccomplex_t IT_0953 = 2*IT_0228;
    const ccomplex_t IT_0954 = IT_0469 + IT_0953;
    const ccomplex_t IT_0955 = 2*conj(IT_0228);
    const ccomplex_t IT_0956 = IT_0470 + IT_0955;
    const ccomplex_t IT_0957 = IT_0479 + IT_0953;
    const ccomplex_t IT_0958 = IT_0480 + IT_0955;
    const ccomplex_t IT_0959 = IT_0349*IT_0370 + IT_0346*conj(IT_0370) + conj
      (IT_0429)*IT_0463 + IT_0429*IT_0464 + IT_0424*(IT_0457 + IT_0462 + IT_0466
      ) + IT_0423*(conj(IT_0457) + conj(IT_0462) + conj(IT_0466)) + conj(IT_0425
      )*IT_0467 + IT_0425*IT_0468 + IT_0458*IT_0534 + IT_0460*IT_0537 + IT_0557
      *IT_0802 + IT_0559*IT_0813 + IT_0452*IT_0837 + IT_0459*IT_0839 + IT_0446
      *IT_0933 + IT_0450*IT_0934 + IT_0447*IT_0943 + IT_0451*IT_0944 + IT_0510
      *IT_0954 + IT_0513*IT_0956 + IT_0502*IT_0957 + IT_0506*IT_0958;
    const ccomplex_t IT_0960 = IT_0238 + IT_0935;
    const ccomplex_t IT_0961 = IT_0005*IT_0960;
    const ccomplex_t IT_0962 = IT_0239 + IT_0938;
    const ccomplex_t IT_0963 = (-4)*IT_0962;
    const ccomplex_t IT_0964 = IT_0961 + IT_0963;
    const ccomplex_t IT_0965 = m_N_3*IT_0964;
    const ccomplex_t IT_0966 = IT_0349*(IT_0457 + IT_0462) + IT_0346*(conj
      (IT_0457) + conj(IT_0462)) + IT_0358*IT_0466 + IT_0354*conj(IT_0466) +
       IT_0452*IT_0820 + IT_0459*IT_0833 + conj(IT_0228)*IT_0847 + IT_0228
      *IT_0852 + conj(IT_0228)*IT_0857 + IT_0228*IT_0862 + IT_0450*IT_0954 +
       IT_0451*IT_0956 + IT_0446*IT_0957 + IT_0447*IT_0958;
    const ccomplex_t IT_0967 = (IT_0052 + -IT_0239)*(s_12 + (-0.5)*IT_0474);
    const ccomplex_t IT_0968 = 2*IT_0967;
    const ccomplex_t IT_0969 = IT_0005*IT_0052;
    const ccomplex_t IT_0970 = (-0.25)*IT_0269;
    const ccomplex_t IT_0971 = IT_0631 + IT_0970;
    const ccomplex_t IT_0972 = IT_0969*IT_0971;
    const ccomplex_t IT_0973 = (-4)*IT_0972;
    const ccomplex_t IT_0974 = (-2)*IT_0188;
    const ccomplex_t IT_0975 = (-2)*IT_0186;
    const ccomplex_t IT_0976 = IT_0974 + IT_0975;
    const ccomplex_t IT_0977 = (-2)*conj(IT_0188);
    const ccomplex_t IT_0978 = (-2)*conj(IT_0186);
    const ccomplex_t IT_0979 = IT_0977 + IT_0978;
    const ccomplex_t IT_0980 = (-2)*IT_0034;
    const ccomplex_t IT_0981 = (-2)*IT_0063;
    const ccomplex_t IT_0982 = (-2)*IT_0107;
    const ccomplex_t IT_0983 = (-2)*IT_0114;
    const ccomplex_t IT_0984 = IT_0980 + IT_0981 + IT_0982 + IT_0983;
    const ccomplex_t IT_0985 = (-2)*conj(IT_0034);
    const ccomplex_t IT_0986 = (-2)*conj(IT_0063);
    const ccomplex_t IT_0987 = (-2)*conj(IT_0107);
    const ccomplex_t IT_0988 = (-2)*conj(IT_0114);
    const ccomplex_t IT_0989 = IT_0985 + IT_0986 + IT_0987 + IT_0988;
    const ccomplex_t IT_0990 = (-2)*IT_0087;
    const ccomplex_t IT_0991 = (-2)*IT_0111;
    const ccomplex_t IT_0992 = (-2)*IT_0117;
    const ccomplex_t IT_0993 = IT_0990 + IT_0991 + IT_0992;
    const ccomplex_t IT_0994 = (-2)*conj(IT_0087);
    const ccomplex_t IT_0995 = (-2)*conj(IT_0111);
    const ccomplex_t IT_0996 = (-2)*conj(IT_0117);
    const ccomplex_t IT_0997 = IT_0994 + IT_0995 + IT_0996;
    const ccomplex_t IT_0998 = (-2)*IT_0056;
    const ccomplex_t IT_0999 = (-2)*conj(IT_0056);
    const ccomplex_t IT_1000 = IT_0521*IT_0976 + IT_0523*IT_0979 + IT_0276
      *IT_0984 + IT_0303*IT_0989 + conj(IT_0189)*IT_0993 + IT_0189*IT_0997 + 
      (conj(IT_0034) + conj(IT_0063) + conj(IT_0107) + conj(IT_0114) + conj
      (IT_0189))*IT_0998 + (IT_0034 + IT_0063 + IT_0107 + IT_0114 + IT_0189)
      *IT_0999;
    const ccomplex_t IT_1001 = 9*IT_0494;
    const ccomplex_t IT_1002 = (IT_0353*IT_0976 + IT_0357*IT_0979 + conj
      (IT_0189)*IT_0984 + IT_0189*IT_0989 + IT_0276*IT_0993 + IT_0303*IT_0997 + 
      (conj(IT_0181) + conj(IT_0185) + IT_0521)*IT_0998 + (IT_0181 + IT_0185 +
       IT_0523)*IT_0999)*IT_1001;
    const ccomplex_t IT_1003 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_1004 = s_12 + IT_0224 + IT_1003;
    const ccomplex_t IT_1005 = m_N_3*IT_1004;
    const ccomplex_t IT_1006 = IT_0256*conj(IT_0422);
    const ccomplex_t IT_1007 = (IT_0259 + IT_0260 + IT_0265 + IT_0266)*IT_0422;
    const ccomplex_t IT_1008 = (IT_0257 + IT_0262 + IT_0263)*conj(IT_0422);
    const ccomplex_t IT_1009 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_1010 = IT_0052 + IT_1009;
    const ccomplex_t IT_1011 = IT_0005*IT_0410;
    const ccomplex_t IT_1012 = -IT_0528 + -IT_1011;
    const ccomplex_t IT_1013 = IT_1010 + IT_1012;
    const ccomplex_t IT_1014 = m_N_3*IT_1013;
    const ccomplex_t IT_1015 = -IT_1014;
    const ccomplex_t IT_1016 = IT_0015*IT_0220 + IT_0226*(IT_0177*IT_0229 +
       IT_0180*IT_0230 + IT_0214*IT_0232 + IT_0219*IT_0233) + (IT_0214*IT_0229 +
       IT_0219*IT_0230 + IT_0177*IT_0232 + IT_0180*IT_0233)*IT_0243 + IT_0255
      *IT_0268 + IT_0275*(IT_0276*IT_0302 + IT_0303*IT_0307 + IT_0308*IT_0329 +
       IT_0330*IT_0333) + IT_0345*(IT_0191*IT_0348 + IT_0216*IT_0351 + IT_0353
      *IT_0355 + IT_0357*IT_0359) + 2*IT_0373 + 2*IT_0377 + 2*IT_0380 + 2
      *IT_0381 + IT_0386*IT_0401 + IT_0407*IT_0408 + IT_0421*IT_0432 + IT_0443
      *IT_0471 + IT_0476*IT_0481 + 2*IT_0491 + 2*IT_0493 + 4*IT_0501*IT_0518 +
       IT_0489*(IT_0490*IT_0519 + IT_0492*IT_0520 + IT_0521*IT_0522 + IT_0523
      *IT_0524 + (conj(IT_0056) + conj(IT_0181) + conj(IT_0185))*IT_0525 + 
      (IT_0056 + IT_0181 + IT_0185)*IT_0526) + 6*IT_0532*(conj(IT_0425)*IT_0429 
      + IT_0425*conj(IT_0429) + 0.166666666666667*IT_0534*IT_0535 +
       0.166666666666667*IT_0537*IT_0538) + 6*IT_0541*(IT_0229*(IT_0228 +
       0.166666666666667*IT_0542) + IT_0230*(conj(IT_0228) + 0.166666666666667
      *IT_0543)) + 2*IT_0555 + IT_0556*IT_0579 + IT_0589*IT_0610 + (IT_0502
      *IT_0592 + IT_0506*IT_0595 + IT_0510*IT_0598 + IT_0513*IT_0601 + conj
      (IT_0370)*(IT_0602 + IT_0603 + IT_0606 + IT_0607) + IT_0370*(IT_0604 +
       IT_0605 + IT_0608 + IT_0609))*IT_0616 + (s_12 + IT_0617 + IT_0618 +
       IT_0619)*IT_0629 + (s_12 + IT_0636 + IT_0637 + IT_0638)*IT_0681 + IT_0687
      *IT_0689 + (s_12 + IT_0617 + IT_0618 + IT_0693)*IT_0694 + (-8)*IT_0702
      *IT_0720 + (-6)*IT_0722*IT_0727 + IT_0731*(conj(IT_0175)*IT_0734 + conj
      (IT_0212)*IT_0737 + conj(IT_0173)*IT_0738 + conj(IT_0210)*IT_0739 +
       IT_0175*IT_0740 + IT_0212*IT_0741) + IT_0742*IT_0751 + IT_0757*IT_0764 + 
      (IT_0276*IT_0505 + IT_0303*IT_0509 + conj(IT_0056)*(IT_0503 + 2*IT_0504 +
       IT_0511) + IT_0308*IT_0512 + IT_0056*(IT_0507 + 2*IT_0508 + IT_0514) +
       IT_0330*IT_0515)*IT_0773 + 2*IT_0781 + IT_0784*IT_0793 + IT_0799 + 2
      *IT_0803 + 2*IT_0812 + 2*IT_0814 + 2*IT_0822 + 2*IT_0832 + 2*IT_0835 + 2
      *IT_0836 + IT_0811*(IT_0430*IT_0502 + IT_0431*IT_0506 + IT_0426*IT_0510 +
       IT_0427*IT_0513 + IT_0423*IT_0561 + IT_0424*IT_0563 + (conj(IT_0425) +
       conj(IT_0429))*IT_0837 + IT_0838*IT_0839) + IT_0843*IT_0876 + IT_0882 +
       IT_0885*IT_0886 + IT_0819*IT_0909 + IT_0831*IT_0910 + IT_0931 + IT_0368
      *IT_0932 + IT_0801*(IT_0813*IT_0838 + conj(IT_0425)*(IT_0802 + IT_0933) +
       conj(IT_0429)*(IT_0802 + IT_0934)) + IT_0945 + IT_0952*IT_0959 + IT_0965
      *IT_0966 + (IT_0232*IT_0346 + 2*conj(IT_0228)*IT_0347 + IT_0228*IT_0349 +
       IT_0229*IT_0354 + IT_0228*IT_0358)*IT_0968 + IT_0973*IT_1000 + IT_1002 +
       IT_1005*(IT_0267*IT_0425 + IT_0264*conj(IT_0425) + IT_0261*IT_0429 +
       IT_0258*conj(IT_0429) + IT_1006 + IT_1007 + IT_1008) + (IT_0261*IT_0425 +
       IT_0258*conj(IT_0425) + IT_0267*IT_0429 + IT_0264*conj(IT_0429) + IT_1006
       + IT_1007 + IT_1008)*IT_1015;
    return create_ccomplex_return(IT_1016);
}

