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
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0002 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0020, -1);
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0027;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + IT_0026 + 
      -IT_0028 + -IT_0029);
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0031, 2);
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = pow(m_N_3, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0002 + IT_0035 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = IT_0025*IT_0038;
    const ccomplex_t IT_0041 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0042 = IT_0021*IT_0041;
    const ccomplex_t IT_0043 = IT_0025*IT_0041;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + IT_0040 + 
      -IT_0042 + -IT_0043);
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0047 = IT_0021*IT_0046;
    const ccomplex_t IT_0048 = IT_0025*IT_0046;
    const ccomplex_t IT_0049 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0050 = IT_0021*IT_0049;
    const ccomplex_t IT_0051 = IT_0025*IT_0049;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0047 + IT_0048 + 
      -IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0045*IT_0053;
    const ccomplex_t IT_0055 = IT_0037*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0002 + IT_0035 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0054*IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0002 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0030;
    const ccomplex_t IT_0063 = IT_0031*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0017*IT_0063;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0002 + IT_0035 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0071 = IT_0021*IT_0070;
    const ccomplex_t IT_0072 = IT_0025*IT_0070;
    const ccomplex_t IT_0073 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0074 = IT_0021*IT_0073;
    const ccomplex_t IT_0075 = IT_0025*IT_0073;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0071 + IT_0072 + 
      -IT_0074 + -IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = IT_0025*IT_0078;
    const ccomplex_t IT_0081 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0082 = IT_0021*IT_0081;
    const ccomplex_t IT_0083 = IT_0025*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0077*IT_0085;
    const ccomplex_t IT_0087 = IT_0069*IT_0086;
    const ccomplex_t IT_0088 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0089 = IT_0025*IT_0088;
    const ccomplex_t IT_0090 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0091 = IT_0021*IT_0090;
    const ccomplex_t IT_0092 = IT_0025*IT_0090;
    const ccomplex_t IT_0093 = IT_0021*IT_0088;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0089 + -IT_0091 + 
      -IT_0092 + IT_0093);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0097 = IT_0021*IT_0096;
    const ccomplex_t IT_0098 = IT_0025*IT_0096;
    const ccomplex_t IT_0099 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0100 = IT_0021*IT_0099;
    const ccomplex_t IT_0101 = IT_0025*IT_0099;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0097 + IT_0098 + 
      -IT_0100 + -IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0095*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0002 + IT_0035 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_13 + IT_0002 + IT_0035 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0105*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0002 + IT_0035 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0086*IT_0114;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = cpow(IT_0062, 2);
    const ccomplex_t IT_0118 = IT_0061*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = cos(alpha);
    const ccomplex_t IT_0121 = IT_0019*IT_0120;
    const ccomplex_t IT_0122 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = sin(alpha);
    const ccomplex_t IT_0125 = IT_0019*IT_0124;
    const ccomplex_t IT_0126 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = IT_0024*IT_0120;
    const ccomplex_t IT_0129 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = IT_0024*IT_0124;
    const ccomplex_t IT_0132 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0123 + -IT_0127 + 
      -IT_0130 + IT_0133);
    const ccomplex_t IT_0135 = cos(beta);
    const ccomplex_t IT_0136 = cpow(IT_0018, -4);
    const ccomplex_t IT_0137 = IT_0020*IT_0120*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = sin(beta);
    const ccomplex_t IT_0139 = IT_0020*IT_0124*IT_0136*IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0020, -3);
    const ccomplex_t IT_0141 = IT_0120*IT_0135*IT_0140;
    const ccomplex_t IT_0142 = IT_0124*IT_0138*IT_0140;
    const ccomplex_t IT_0143 = cpow(IT_0018, -2);
    const ccomplex_t IT_0144 = IT_0024*IT_0120*IT_0135*IT_0143;
    const ccomplex_t IT_0145 = IT_0024*IT_0124*IT_0138*IT_0143;
    const ccomplex_t IT_0146 = cpow(IT_0020, -2);
    const ccomplex_t IT_0147 = cpow(IT_0143 + IT_0146, -1);
    const ccomplex_t IT_0148 = (IT_0137 + IT_0139 + IT_0141 + IT_0142 + 2
      *IT_0144 + 2*IT_0145)*IT_0147;
    const ccomplex_t IT_0149 = m_W*e_em;
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_12 + (-2)*IT_0035 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0134*IT_0154;
    const ccomplex_t IT_0156 = IT_0121*IT_0126;
    const ccomplex_t IT_0157 = IT_0122*IT_0125;
    const ccomplex_t IT_0158 = IT_0128*IT_0132;
    const ccomplex_t IT_0159 = IT_0129*IT_0131;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0156 + IT_0157 + 
      -IT_0158 + -IT_0159);
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0020*IT_0120*IT_0136*IT_0138;
    const ccomplex_t IT_0163 = IT_0020*IT_0124*IT_0135*IT_0136;
    const ccomplex_t IT_0164 = IT_0120*IT_0138*IT_0140;
    const ccomplex_t IT_0165 = IT_0124*IT_0135*IT_0140;
    const ccomplex_t IT_0166 = IT_0024*IT_0120*IT_0138*IT_0143;
    const ccomplex_t IT_0167 = IT_0024*IT_0124*IT_0135*IT_0143;
    const ccomplex_t IT_0168 = IT_0147*(IT_0162 + -IT_0163 + IT_0164 + 
      -IT_0165 + 2*IT_0166 + (-2)*IT_0167);
    const ccomplex_t IT_0169 = IT_0149*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0035 + -reg_prop, -1);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0161*IT_0173;
    const ccomplex_t IT_0175 = -IT_0155 + (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = 2*IT_0175;
    const ccomplex_t IT_0177 = 0.5*IT_0174;
    const ccomplex_t IT_0178 = (-2)*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = 2*conj(IT_0175);
    const ccomplex_t IT_0181 = (-2)*conj(IT_0177);
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = IT_0055 + 0.5*IT_0064;
    const ccomplex_t IT_0184 = IT_0017*IT_0117;
    const ccomplex_t IT_0185 = -IT_0058 + 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0087 + IT_0112;
    const ccomplex_t IT_0187 = (-0.5)*IT_0066 + -IT_0108 + -IT_0115;
    const ccomplex_t IT_0188 = IT_0032*IT_0061;
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0191 = IT_0121*IT_0190;
    const ccomplex_t IT_0192 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0193 = IT_0125*IT_0192;
    const ccomplex_t IT_0194 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0195 = IT_0128*IT_0194;
    const ccomplex_t IT_0196 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0197 = IT_0131*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + -IT_0193 + 
      -IT_0195 + IT_0197);
    const ccomplex_t IT_0199 = IT_0154*IT_0198;
    const ccomplex_t IT_0200 = IT_0121*IT_0192;
    const ccomplex_t IT_0201 = IT_0125*IT_0190;
    const ccomplex_t IT_0202 = IT_0128*IT_0196;
    const ccomplex_t IT_0203 = IT_0131*IT_0194;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0200 + IT_0201 + 
      -IT_0202 + -IT_0203);
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0173*IT_0205;
    const ccomplex_t IT_0207 = -IT_0199 + (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = 2*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0206;
    const ccomplex_t IT_0210 = (-2)*IT_0209;
    const ccomplex_t IT_0211 = IT_0208 + IT_0210;
    const ccomplex_t IT_0212 = 2*conj(IT_0207);
    const ccomplex_t IT_0213 = (-2)*conj(IT_0209);
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = (conj(IT_0034) + conj(IT_0055) + conj(IT_0059) 
      + conj(IT_0065) + conj(IT_0067) + conj(IT_0087) + conj(IT_0109) + conj
      (IT_0112) + conj(IT_0116) + conj(IT_0119))*IT_0179 + (IT_0034 + IT_0055 +
       IT_0059 + IT_0065 + IT_0067 + IT_0087 + IT_0109 + IT_0112 + IT_0116 +
       IT_0119)*IT_0182 + (conj(IT_0183) + conj(IT_0185) + conj(IT_0186) + conj
      (IT_0187) + conj(IT_0189))*IT_0211 + (IT_0183 + IT_0185 + IT_0186 +
       IT_0187 + IT_0189)*IT_0214;
    const ccomplex_t IT_0216 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0217 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0218 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0219 = -IT_0217 + -IT_0218;
    const ccomplex_t IT_0220 = s_12 + IT_0216 + IT_0219;
    const ccomplex_t IT_0221 = m_N_3*IT_0220;
    const ccomplex_t IT_0222 = (-2)*IT_0055 + -IT_0064;
    const ccomplex_t IT_0223 = (-2)*IT_0087 + (-2)*IT_0112 + IT_0118;
    const ccomplex_t IT_0224 = conj(IT_0222) + conj(IT_0223);
    const ccomplex_t IT_0225 = IT_0222 + IT_0223;
    const ccomplex_t IT_0226 = (-2)*IT_0087 + (-2)*IT_0112 + IT_0188;
    const ccomplex_t IT_0227 = conj(IT_0222) + conj(IT_0226);
    const ccomplex_t IT_0228 = IT_0222 + IT_0226;
    const ccomplex_t IT_0229 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0230 = IT_0035 + IT_0229;
    const ccomplex_t IT_0231 = pow(s_23, 2);
    const ccomplex_t IT_0232 = IT_0005*IT_0231;
    const ccomplex_t IT_0233 = pow(s_24, 2);
    const ccomplex_t IT_0234 = IT_0005*IT_0233;
    const ccomplex_t IT_0235 = -IT_0232 + -IT_0234;
    const ccomplex_t IT_0236 = IT_0230 + IT_0235;
    const ccomplex_t IT_0237 = m_N_3*IT_0236;
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = s_23*IT_0002;
    const ccomplex_t IT_0240 = IT_0001*IT_0239;
    const ccomplex_t IT_0241 = s_24*IT_0002;
    const ccomplex_t IT_0242 = IT_0005*IT_0241;
    const ccomplex_t IT_0243 = s_23*s_34;
    const ccomplex_t IT_0244 = (-2)*IT_0243;
    const ccomplex_t IT_0245 = IT_0241 + IT_0244;
    const ccomplex_t IT_0246 = IT_0005*IT_0245;
    const ccomplex_t IT_0247 = (-0.5)*IT_0240 + 0.5*IT_0242 + (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = s_24 + IT_0247;
    const ccomplex_t IT_0249 = m_N_3*IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0177;
    const ccomplex_t IT_0252 = (-2)*IT_0175;
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 2*conj(IT_0177);
    const ccomplex_t IT_0255 = (-2)*conj(IT_0175);
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 2*IT_0209;
    const ccomplex_t IT_0258 = (-2)*IT_0207;
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = 2*conj(IT_0209);
    const ccomplex_t IT_0261 = (-2)*conj(IT_0207);
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = (conj(IT_0183) + conj(IT_0185) + conj(IT_0186) 
      + conj(IT_0187) + conj(IT_0189))*IT_0253 + (IT_0183 + IT_0185 + IT_0186 +
       IT_0187 + IT_0189)*IT_0256 + (conj(IT_0034) + conj(IT_0055) + conj
      (IT_0059) + conj(IT_0065) + conj(IT_0067) + conj(IT_0087) + conj(IT_0109) 
      + conj(IT_0112) + conj(IT_0116) + conj(IT_0119))*IT_0259 + (IT_0034 +
       IT_0055 + IT_0059 + IT_0065 + IT_0067 + IT_0087 + IT_0109 + IT_0112 +
       IT_0116 + IT_0119)*IT_0262;
    const ccomplex_t IT_0264 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0265 = IT_0035 + IT_0264;
    const ccomplex_t IT_0266 = pow(s_13, 2);
    const ccomplex_t IT_0267 = IT_0005*IT_0266;
    const ccomplex_t IT_0268 = pow(s_14, 2);
    const ccomplex_t IT_0269 = IT_0005*IT_0268;
    const ccomplex_t IT_0270 = -IT_0267 + -IT_0269;
    const ccomplex_t IT_0271 = IT_0265 + IT_0270;
    const ccomplex_t IT_0272 = m_N_3*IT_0271;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = -IT_0033 + 2*IT_0058;
    const ccomplex_t IT_0275 = IT_0066 + 2*IT_0108 + 2*IT_0115;
    const ccomplex_t IT_0276 = conj(IT_0274) + conj(IT_0275);
    const ccomplex_t IT_0277 = IT_0274 + IT_0275;
    const ccomplex_t IT_0278 = 2*IT_0058 + -IT_0184;
    const ccomplex_t IT_0279 = conj(IT_0275) + conj(IT_0278);
    const ccomplex_t IT_0280 = IT_0275 + IT_0278;
    const ccomplex_t IT_0281 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0282 = s_12 + IT_0219 + IT_0281;
    const ccomplex_t IT_0283 = m_N_3*IT_0282;
    const ccomplex_t IT_0284 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0285 = pow(m_Z, 4);
    const ccomplex_t IT_0286 = s_24*IT_0285;
    const ccomplex_t IT_0287 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = IT_0000*IT_0289;
    const ccomplex_t IT_0291 = (-0.666666666666667)*IT_0284 + (
      -0.166666666666667)*IT_0290;
    const ccomplex_t IT_0292 = s_24 + IT_0291;
    const ccomplex_t IT_0293 = m_N_3*IT_0292;
    const ccomplex_t IT_0294 = (-6)*IT_0293;
    const ccomplex_t IT_0295 = conj(IT_0055) + conj(IT_0065);
    const ccomplex_t IT_0296 = conj(IT_0087) + IT_0295;
    const ccomplex_t IT_0297 = conj(IT_0112) + conj(IT_0119) + IT_0296;
    const ccomplex_t IT_0298 = cpow(IT_0053, 2);
    const ccomplex_t IT_0299 = m_N_4*IT_0298;
    const ccomplex_t IT_0300 = m_N_3*IT_0054;
    const ccomplex_t IT_0301 = IT_0057*(IT_0299 + -IT_0300);
    const ccomplex_t IT_0302 = -IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0302;
    const ccomplex_t IT_0304 = m_N_3*IT_0032;
    const ccomplex_t IT_0305 = m_N_3*IT_0063;
    const ccomplex_t IT_0306 = IT_0017*(IT_0304 + -IT_0305);
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = 0.5*IT_0307;
    const ccomplex_t IT_0309 = 2*IT_0308;
    const ccomplex_t IT_0310 = IT_0303 + IT_0309;
    const ccomplex_t IT_0311 = m_N_3*IT_0105;
    const ccomplex_t IT_0312 = cpow(IT_0104, 2);
    const ccomplex_t IT_0313 = m_N_2*IT_0312;
    const ccomplex_t IT_0314 = IT_0107*(IT_0311 + -IT_0313);
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = (-2)*IT_0315;
    const ccomplex_t IT_0317 = -IT_0304;
    const ccomplex_t IT_0318 = IT_0305 + IT_0317;
    const ccomplex_t IT_0319 = IT_0017*IT_0318;
    const ccomplex_t IT_0320 = cpow(IT_0085, 2);
    const ccomplex_t IT_0321 = m_N_1*IT_0320;
    const ccomplex_t IT_0322 = m_N_3*IT_0086;
    const ccomplex_t IT_0323 = -IT_0322;
    const ccomplex_t IT_0324 = IT_0321 + IT_0323;
    const ccomplex_t IT_0325 = IT_0114*IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0319 + IT_0325;
    const ccomplex_t IT_0327 = (-2)*IT_0326;
    const ccomplex_t IT_0328 = IT_0310 + IT_0316 + IT_0327;
    const ccomplex_t IT_0329 = IT_0055 + IT_0065;
    const ccomplex_t IT_0330 = IT_0087 + IT_0112 + IT_0119 + IT_0329;
    const ccomplex_t IT_0331 = 2*conj(IT_0302);
    const ccomplex_t IT_0332 = 2*conj(IT_0308);
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = (-2)*conj(IT_0315);
    const ccomplex_t IT_0335 = (-2)*conj(IT_0326);
    const ccomplex_t IT_0336 = IT_0333 + IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = conj(IT_0183) + conj(IT_0186) + conj(IT_0189);
    const ccomplex_t IT_0338 = cpow(IT_0045, 2);
    const ccomplex_t IT_0339 = m_N_4*IT_0338;
    const ccomplex_t IT_0340 = -IT_0057*(IT_0300 + -IT_0339);
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = 2*IT_0341;
    const ccomplex_t IT_0343 = m_N_3*IT_0117;
    const ccomplex_t IT_0344 = IT_0017*(IT_0305 + -IT_0343);
    const ccomplex_t IT_0345 = 0.5*IT_0344;
    const ccomplex_t IT_0346 = 2*IT_0345;
    const ccomplex_t IT_0347 = IT_0342 + IT_0346;
    const ccomplex_t IT_0348 = cpow(IT_0077, 2);
    const ccomplex_t IT_0349 = m_N_1*IT_0348;
    const ccomplex_t IT_0350 = IT_0114*(IT_0322 + -IT_0349);
    const ccomplex_t IT_0351 = -IT_0350;
    const ccomplex_t IT_0352 = (-2)*IT_0351;
    const ccomplex_t IT_0353 = -IT_0343;
    const ccomplex_t IT_0354 = IT_0305 + IT_0353;
    const ccomplex_t IT_0355 = IT_0017*IT_0354;
    const ccomplex_t IT_0356 = cpow(IT_0095, 2);
    const ccomplex_t IT_0357 = m_N_2*IT_0356;
    const ccomplex_t IT_0358 = -IT_0311;
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = IT_0107*IT_0359;
    const ccomplex_t IT_0361 = (-0.5)*IT_0355 + IT_0360;
    const ccomplex_t IT_0362 = (-2)*IT_0361;
    const ccomplex_t IT_0363 = IT_0347 + IT_0352 + IT_0362;
    const ccomplex_t IT_0364 = IT_0183 + IT_0186 + IT_0189;
    const ccomplex_t IT_0365 = 2*conj(IT_0341);
    const ccomplex_t IT_0366 = 2*conj(IT_0345);
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = (-2)*conj(IT_0351);
    const ccomplex_t IT_0369 = (-2)*conj(IT_0361);
    const ccomplex_t IT_0370 = IT_0367 + IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = (-0.666666666666667)*IT_0246 +
       0.166666666666667*IT_0290;
    const ccomplex_t IT_0372 = s_24 + IT_0371;
    const ccomplex_t IT_0373 = m_N_3*IT_0372;
    const ccomplex_t IT_0374 = 6*IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0034) + conj(IT_0059);
    const ccomplex_t IT_0376 = conj(IT_0067) + conj(IT_0109) + conj(IT_0116) +
       IT_0375;
    const ccomplex_t IT_0377 = IT_0034 + IT_0059;
    const ccomplex_t IT_0378 = IT_0067 + IT_0109 + IT_0116 + IT_0377;
    const ccomplex_t IT_0379 = conj(IT_0185) + conj(IT_0187);
    const ccomplex_t IT_0380 = IT_0185 + IT_0187;
    const ccomplex_t IT_0381 = m_N_3*IT_0005;
    const ccomplex_t IT_0382 = s_14*IT_0285;
    const ccomplex_t IT_0383 = IT_0005*IT_0382;
    const ccomplex_t IT_0384 = (-2)*IT_0008 + (-0.5)*IT_0383;
    const ccomplex_t IT_0385 = IT_0006 + IT_0384;
    const ccomplex_t IT_0386 = IT_0381*IT_0385;
    const ccomplex_t IT_0387 = (-2)*IT_0386;
    const ccomplex_t IT_0388 = IT_0005*IT_0286;
    const ccomplex_t IT_0389 = (-2)*IT_0243 + (-0.5)*IT_0388;
    const ccomplex_t IT_0390 = IT_0241 + IT_0389;
    const ccomplex_t IT_0391 = IT_0381*IT_0390;
    const ccomplex_t IT_0392 = (-2)*IT_0391;
    const ccomplex_t IT_0393 = IT_0323 + IT_0349;
    const ccomplex_t IT_0394 = IT_0069*IT_0393;
    const ccomplex_t IT_0395 = IT_0111*IT_0359;
    const ccomplex_t IT_0396 = IT_0037*IT_0339;
    const ccomplex_t IT_0397 = IT_0394 + IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = 2*IT_0397;
    const ccomplex_t IT_0399 = -IT_0061*(IT_0305 + -IT_0343);
    const ccomplex_t IT_0400 = -IT_0399;
    const ccomplex_t IT_0401 = (-0.5)*IT_0400;
    const ccomplex_t IT_0402 = 2*IT_0401;
    const ccomplex_t IT_0403 = IT_0398 + IT_0402;
    const ccomplex_t IT_0404 = IT_0037*IT_0300;
    const ccomplex_t IT_0405 = IT_0061*IT_0354;
    const ccomplex_t IT_0406 = IT_0404 + 0.5*IT_0405;
    const ccomplex_t IT_0407 = (-2)*IT_0406;
    const ccomplex_t IT_0408 = IT_0403 + IT_0407;
    const ccomplex_t IT_0409 = 2*conj(IT_0397);
    const ccomplex_t IT_0410 = 2*conj(IT_0401);
    const ccomplex_t IT_0411 = IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = (-2)*conj(IT_0406);
    const ccomplex_t IT_0413 = IT_0411 + IT_0412;
    const ccomplex_t IT_0414 = IT_0069*IT_0324;
    const ccomplex_t IT_0415 = IT_0037*IT_0299;
    const ccomplex_t IT_0416 = IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = 2*IT_0416;
    const ccomplex_t IT_0418 = -IT_0111*(IT_0311 + -IT_0313);
    const ccomplex_t IT_0419 = -IT_0418;
    const ccomplex_t IT_0420 = -IT_0419;
    const ccomplex_t IT_0421 = 2*IT_0420;
    const ccomplex_t IT_0422 = -IT_0061*(IT_0304 + -IT_0305);
    const ccomplex_t IT_0423 = -IT_0422;
    const ccomplex_t IT_0424 = 0.5*IT_0423;
    const ccomplex_t IT_0425 = 2*IT_0424;
    const ccomplex_t IT_0426 = IT_0417 + IT_0421 + IT_0425;
    const ccomplex_t IT_0427 = IT_0061*IT_0318;
    const ccomplex_t IT_0428 = IT_0404 + 0.5*IT_0427;
    const ccomplex_t IT_0429 = (-2)*IT_0428;
    const ccomplex_t IT_0430 = IT_0426 + IT_0429;
    const ccomplex_t IT_0431 = 2*conj(IT_0416);
    const ccomplex_t IT_0432 = 2*conj(IT_0420);
    const ccomplex_t IT_0433 = 2*conj(IT_0424);
    const ccomplex_t IT_0434 = IT_0431 + IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = (-2)*conj(IT_0428);
    const ccomplex_t IT_0436 = IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0438 = (-2)*IT_0437;
    const ccomplex_t IT_0439 = IT_0382 + IT_0438;
    const ccomplex_t IT_0440 = IT_0000*IT_0439;
    const ccomplex_t IT_0441 = (-0.666666666666667)*IT_0011 +
       0.166666666666667*IT_0440;
    const ccomplex_t IT_0442 = s_14 + IT_0441;
    const ccomplex_t IT_0443 = m_N_3*IT_0442;
    const ccomplex_t IT_0444 = 6*IT_0443;
    const ccomplex_t IT_0445 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0446 = (-0.166666666666667)*IT_0440 + (
      -0.666666666666667)*IT_0445;
    const ccomplex_t IT_0447 = s_14 + IT_0446;
    const ccomplex_t IT_0448 = m_N_3*IT_0447;
    const ccomplex_t IT_0449 = (-6)*IT_0448;
    const ccomplex_t IT_0450 = s_12*s_14*IT_0002;
    const ccomplex_t IT_0451 = s_24*IT_0266;
    const ccomplex_t IT_0452 = s_12*s_13*s_34;
    const ccomplex_t IT_0453 = s_24*IT_0002*IT_0035;
    const ccomplex_t IT_0454 = 2*IT_0451 + (-2)*IT_0452 + -IT_0453;
    const ccomplex_t IT_0455 = IT_0450 + IT_0454;
    const ccomplex_t IT_0456 = IT_0005*IT_0455;
    const ccomplex_t IT_0457 = s_24*IT_0035;
    const ccomplex_t IT_0458 = 2*IT_0457;
    const ccomplex_t IT_0459 = s_12*s_14;
    const ccomplex_t IT_0460 = (-2)*IT_0459;
    const ccomplex_t IT_0461 = IT_0456 + IT_0458 + IT_0460;
    const ccomplex_t IT_0462 = 2*IT_0278;
    const ccomplex_t IT_0463 = 2*IT_0275;
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = 2*conj(IT_0278);
    const ccomplex_t IT_0466 = 2*conj(IT_0275);
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = 2*IT_0274;
    const ccomplex_t IT_0469 = IT_0463 + IT_0468;
    const ccomplex_t IT_0470 = 2*conj(IT_0274);
    const ccomplex_t IT_0471 = IT_0466 + IT_0470;
    const ccomplex_t IT_0472 = s_13*s_14*s_23;
    const ccomplex_t IT_0473 = s_23*s_34*IT_0035;
    const ccomplex_t IT_0474 = -IT_0453 + (-2)*IT_0472 + 2*IT_0473;
    const ccomplex_t IT_0475 = IT_0450 + IT_0474;
    const ccomplex_t IT_0476 = IT_0005*IT_0475;
    const ccomplex_t IT_0477 = (-0.5)*IT_0476;
    const ccomplex_t IT_0478 = -IT_0457;
    const ccomplex_t IT_0479 = IT_0459 + IT_0477 + IT_0478;
    const ccomplex_t IT_0480 = 4*IT_0278;
    const ccomplex_t IT_0481 = 4*IT_0275;
    const ccomplex_t IT_0482 = IT_0480 + IT_0481;
    const ccomplex_t IT_0483 = 4*conj(IT_0278);
    const ccomplex_t IT_0484 = 4*conj(IT_0275);
    const ccomplex_t IT_0485 = IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = 4*IT_0274;
    const ccomplex_t IT_0487 = IT_0481 + IT_0486;
    const ccomplex_t IT_0488 = 4*conj(IT_0274);
    const ccomplex_t IT_0489 = IT_0484 + IT_0488;
    const ccomplex_t IT_0490 = s_14*s_23;
    const ccomplex_t IT_0491 = s_13*s_24;
    const ccomplex_t IT_0492 = s_12*s_34;
    const ccomplex_t IT_0493 = -IT_0491 + -IT_0492;
    const ccomplex_t IT_0494 = IT_0490 + IT_0493;
    const ccomplex_t IT_0495 = IT_0001*IT_0494;
    const ccomplex_t IT_0496 = (-0.25)*IT_0495;
    const ccomplex_t IT_0497 = s_12*IT_0002;
    const ccomplex_t IT_0498 = IT_0005*IT_0497;
    const ccomplex_t IT_0499 = (-0.5)*IT_0498;
    const ccomplex_t IT_0500 = s_12 + IT_0496 + IT_0499;
    const ccomplex_t IT_0501 = 8*IT_0175;
    const ccomplex_t IT_0502 = (-8)*IT_0177;
    const ccomplex_t IT_0503 = IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = 8*conj(IT_0175);
    const ccomplex_t IT_0505 = (-8)*conj(IT_0177);
    const ccomplex_t IT_0506 = IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = conj(IT_0397) + conj(IT_0401);
    const ccomplex_t IT_0508 = 8*IT_0177;
    const ccomplex_t IT_0509 = (-8)*IT_0175;
    const ccomplex_t IT_0510 = IT_0508 + IT_0509;
    const ccomplex_t IT_0511 = IT_0397 + IT_0401;
    const ccomplex_t IT_0512 = 8*conj(IT_0177);
    const ccomplex_t IT_0513 = (-8)*conj(IT_0175);
    const ccomplex_t IT_0514 = IT_0512 + IT_0513;
    const ccomplex_t IT_0515 = 8*IT_0207;
    const ccomplex_t IT_0516 = (-8)*IT_0209;
    const ccomplex_t IT_0517 = IT_0515 + IT_0516;
    const ccomplex_t IT_0518 = 8*conj(IT_0207);
    const ccomplex_t IT_0519 = (-8)*conj(IT_0209);
    const ccomplex_t IT_0520 = IT_0518 + IT_0519;
    const ccomplex_t IT_0521 = conj(IT_0416) + conj(IT_0420) + conj(IT_0424);
    const ccomplex_t IT_0522 = 8*IT_0209;
    const ccomplex_t IT_0523 = (-8)*IT_0207;
    const ccomplex_t IT_0524 = IT_0522 + IT_0523;
    const ccomplex_t IT_0525 = IT_0416 + IT_0420 + IT_0424;
    const ccomplex_t IT_0526 = 8*conj(IT_0209);
    const ccomplex_t IT_0527 = (-8)*conj(IT_0207);
    const ccomplex_t IT_0528 = IT_0526 + IT_0527;
    const ccomplex_t IT_0529 = IT_0490 + IT_0492;
    const ccomplex_t IT_0530 = -IT_0491;
    const ccomplex_t IT_0531 = IT_0529 + IT_0530;
    const ccomplex_t IT_0532 = IT_0001*IT_0531;
    const ccomplex_t IT_0533 = 0.25*IT_0532;
    const ccomplex_t IT_0534 = s_12 + IT_0499 + IT_0533;
    const ccomplex_t IT_0535 = conj(IT_0341) + conj(IT_0345);
    const ccomplex_t IT_0536 = IT_0341 + IT_0345;
    const ccomplex_t IT_0537 = conj(IT_0351) + conj(IT_0361);
    const ccomplex_t IT_0538 = IT_0351 + IT_0361;
    const ccomplex_t IT_0539 = conj(IT_0302) + conj(IT_0308);
    const ccomplex_t IT_0540 = IT_0302 + IT_0308;
    const ccomplex_t IT_0541 = conj(IT_0315) + conj(IT_0326);
    const ccomplex_t IT_0542 = IT_0315 + IT_0326;
    const ccomplex_t IT_0543 = s_12*s_24;
    const ccomplex_t IT_0544 = s_12*s_24*IT_0002;
    const ccomplex_t IT_0545 = s_14*IT_0002*IT_0035;
    const ccomplex_t IT_0546 = s_14*IT_0231;
    const ccomplex_t IT_0547 = s_12*s_23*s_34;
    const ccomplex_t IT_0548 = -IT_0545 + 2*IT_0546 + (-2)*IT_0547;
    const ccomplex_t IT_0549 = IT_0544 + IT_0548;
    const ccomplex_t IT_0550 = IT_0005*IT_0549;
    const ccomplex_t IT_0551 = 0.5*IT_0550;
    const ccomplex_t IT_0552 = -IT_0551;
    const ccomplex_t IT_0553 = s_14*IT_0035;
    const ccomplex_t IT_0554 = -IT_0553;
    const ccomplex_t IT_0555 = IT_0543 + IT_0552 + IT_0554;
    const ccomplex_t IT_0556 = (-4)*IT_0222;
    const ccomplex_t IT_0557 = (-4)*IT_0226;
    const ccomplex_t IT_0558 = IT_0556 + IT_0557;
    const ccomplex_t IT_0559 = (-4)*conj(IT_0222);
    const ccomplex_t IT_0560 = (-4)*conj(IT_0226);
    const ccomplex_t IT_0561 = IT_0559 + IT_0560;
    const ccomplex_t IT_0562 = (-4)*IT_0223;
    const ccomplex_t IT_0563 = IT_0556 + IT_0562;
    const ccomplex_t IT_0564 = (-4)*conj(IT_0223);
    const ccomplex_t IT_0565 = IT_0559 + IT_0564;
    const ccomplex_t IT_0566 = IT_0000*IT_0382;
    const ccomplex_t IT_0567 = (-10)*s_14;
    const ccomplex_t IT_0568 = IT_0566 + IT_0567;
    const ccomplex_t IT_0569 = m_N_3*IT_0568;
    const ccomplex_t IT_0570 = IT_0352 + IT_0362;
    const ccomplex_t IT_0571 = IT_0368 + IT_0369;
    const ccomplex_t IT_0572 = s_13*s_23*s_24;
    const ccomplex_t IT_0573 = s_13*s_34*IT_0035;
    const ccomplex_t IT_0574 = -IT_0545 + (-2)*IT_0572 + 2*IT_0573;
    const ccomplex_t IT_0575 = IT_0544 + IT_0574;
    const ccomplex_t IT_0576 = IT_0005*IT_0575;
    const ccomplex_t IT_0577 = 2*IT_0553;
    const ccomplex_t IT_0578 = (-2)*IT_0543;
    const ccomplex_t IT_0579 = IT_0576 + IT_0577 + IT_0578;
    const ccomplex_t IT_0580 = (-2)*IT_0222;
    const ccomplex_t IT_0581 = (-2)*IT_0226;
    const ccomplex_t IT_0582 = IT_0580 + IT_0581;
    const ccomplex_t IT_0583 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0584 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0585 = IT_0583 + IT_0584;
    const ccomplex_t IT_0586 = (-2)*IT_0223;
    const ccomplex_t IT_0587 = IT_0580 + IT_0586;
    const ccomplex_t IT_0588 = (-2)*conj(IT_0223);
    const ccomplex_t IT_0589 = IT_0583 + IT_0588;
    const ccomplex_t IT_0590 = IT_0002*IT_0035;
    const ccomplex_t IT_0591 = -IT_0590;
    const ccomplex_t IT_0592 = IT_0268 + IT_0591;
    const ccomplex_t IT_0593 = IT_0005*IT_0592;
    const ccomplex_t IT_0594 = -IT_0035;
    const ccomplex_t IT_0595 = IT_0269 + IT_0594;
    const ccomplex_t IT_0596 = (-4)*IT_0595;
    const ccomplex_t IT_0597 = IT_0593 + IT_0596;
    const ccomplex_t IT_0598 = m_N_3*IT_0597;
    const ccomplex_t IT_0599 = (-2)*IT_0278;
    const ccomplex_t IT_0600 = (-2)*IT_0275;
    const ccomplex_t IT_0601 = IT_0599 + IT_0600;
    const ccomplex_t IT_0602 = (-2)*conj(IT_0278);
    const ccomplex_t IT_0603 = (-2)*conj(IT_0275);
    const ccomplex_t IT_0604 = IT_0602 + IT_0603;
    const ccomplex_t IT_0605 = (-2)*IT_0274;
    const ccomplex_t IT_0606 = IT_0600 + IT_0605;
    const ccomplex_t IT_0607 = (-2)*conj(IT_0274);
    const ccomplex_t IT_0608 = IT_0603 + IT_0607;
    const ccomplex_t IT_0609 = IT_0598*(conj(IT_0275)*IT_0310 + IT_0275
      *IT_0333 + IT_0276*IT_0347 + IT_0277*IT_0367 + IT_0462*IT_0539 + IT_0465
      *IT_0540 + IT_0541*IT_0601 + IT_0542*IT_0604 + IT_0537*IT_0606 + IT_0538
      *IT_0608);
    const ccomplex_t IT_0610 = s_13*s_23;
    const ccomplex_t IT_0611 = s_14*s_24;
    const ccomplex_t IT_0612 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0613 = (-2)*IT_0610 + -IT_0611 + 2*IT_0612;
    const ccomplex_t IT_0614 = IT_0497 + IT_0613;
    const ccomplex_t IT_0615 = IT_0005*IT_0614;
    const ccomplex_t IT_0616 = -IT_0218 + (-0.5)*IT_0615;
    const ccomplex_t IT_0617 = s_12 + IT_0616;
    const ccomplex_t IT_0618 = m_N_3*IT_0617;
    const ccomplex_t IT_0619 = (-2)*IT_0618;
    const ccomplex_t IT_0620 = 2*IT_0406;
    const ccomplex_t IT_0621 = 2*conj(IT_0406);
    const ccomplex_t IT_0622 = 2*IT_0428;
    const ccomplex_t IT_0623 = 2*conj(IT_0428);
    const ccomplex_t IT_0624 = 2*IT_0222;
    const ccomplex_t IT_0625 = 2*conj(IT_0222);
    const ccomplex_t IT_0626 = IT_0619*(IT_0521*IT_0582 + IT_0525*IT_0585 +
       IT_0507*IT_0587 + IT_0511*IT_0589 + IT_0224*IT_0620 + IT_0225*IT_0621 +
       conj(IT_0226)*IT_0622 + IT_0226*IT_0623 + conj(IT_0428)*IT_0624 + IT_0428
      *IT_0625);
    const ccomplex_t IT_0627 = (-2)*IT_0610;
    const ccomplex_t IT_0628 = IT_0497 + IT_0627;
    const ccomplex_t IT_0629 = IT_0005*IT_0628;
    const ccomplex_t IT_0630 = (s_12 + -IT_0218)*(s_12 + (-0.5)*IT_0629);
    const ccomplex_t IT_0631 = (-2)*IT_0630;
    const ccomplex_t IT_0632 = IT_0227*IT_0462 + conj(IT_0223)*IT_0463 + conj
      (IT_0226)*IT_0463 + IT_0228*IT_0465 + IT_0223*IT_0466 + IT_0226*IT_0466 +
       IT_0224*IT_0468 + IT_0225*IT_0470 + 2*conj(IT_0275)*IT_0624 + 2*IT_0275
      *IT_0625;
    const ccomplex_t IT_0633 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0634 = 0.5*IT_0233 + (-0.5)*IT_0590 + -IT_0633;
    const ccomplex_t IT_0635 = IT_0231 + IT_0634;
    const ccomplex_t IT_0636 = IT_0005*IT_0635;
    const ccomplex_t IT_0637 = -IT_0234;
    const ccomplex_t IT_0638 = IT_0035 + IT_0636 + IT_0637;
    const ccomplex_t IT_0639 = m_N_3*IT_0638;
    const ccomplex_t IT_0640 = 4*IT_0222;
    const ccomplex_t IT_0641 = 4*conj(IT_0222);
    const ccomplex_t IT_0642 = IT_0537*IT_0558 + IT_0538*IT_0561 + IT_0541
      *IT_0563 + IT_0542*IT_0565 + 4*IT_0539*(IT_0223 + 0.25*IT_0640) + 4
      *IT_0535*(IT_0226 + 0.25*IT_0640) + 4*IT_0540*(conj(IT_0223) + 0.25
      *IT_0641) + 4*IT_0536*(conj(IT_0226) + 0.25*IT_0641);
    const ccomplex_t IT_0643 = -IT_0269;
    const ccomplex_t IT_0644 = IT_0035 + IT_0643;
    const ccomplex_t IT_0645 = IT_0035*IT_0644;
    const ccomplex_t IT_0646 = -IT_0645;
    const ccomplex_t IT_0647 = 6*conj(IT_0278);
    const ccomplex_t IT_0648 = IT_0035 + IT_0637;
    const ccomplex_t IT_0649 = IT_0035*IT_0648;
    const ccomplex_t IT_0650 = -IT_0649;
    const ccomplex_t IT_0651 = 6*IT_0222;
    const ccomplex_t IT_0652 = 6*conj(IT_0222);
    const ccomplex_t IT_0653 = 1.125*IT_0035;
    const ccomplex_t IT_0654 = 16*IT_0428;
    const ccomplex_t IT_0655 = 16*conj(IT_0428);
    const ccomplex_t IT_0656 = 16*IT_0308;
    const ccomplex_t IT_0657 = 16*conj(IT_0302);
    const ccomplex_t IT_0658 = 16*conj(IT_0308);
    const ccomplex_t IT_0659 = (-16)*IT_0345;
    const ccomplex_t IT_0660 = (-16)*conj(IT_0341);
    const ccomplex_t IT_0661 = (-16)*conj(IT_0345);
    const ccomplex_t IT_0662 = 16*IT_0341;
    const ccomplex_t IT_0663 = 16*IT_0345;
    const ccomplex_t IT_0664 = (-16)*IT_0351;
    const ccomplex_t IT_0665 = (-16)*IT_0361;
    const ccomplex_t IT_0666 = 16*conj(IT_0341);
    const ccomplex_t IT_0667 = 16*conj(IT_0345);
    const ccomplex_t IT_0668 = (-16)*conj(IT_0351);
    const ccomplex_t IT_0669 = (-16)*conj(IT_0361);
    const ccomplex_t IT_0670 = 16*IT_0406;
    const ccomplex_t IT_0671 = (-16)*IT_0397;
    const ccomplex_t IT_0672 = (-16)*IT_0401;
    const ccomplex_t IT_0673 = 16*conj(IT_0406);
    const ccomplex_t IT_0674 = (-16)*conj(IT_0397);
    const ccomplex_t IT_0675 = (-16)*conj(IT_0401);
    const ccomplex_t IT_0676 = (-16)*conj(IT_0406)*IT_0428 + (-16)*IT_0406
      *conj(IT_0428) + IT_0507*IT_0654 + IT_0511*IT_0655 + IT_0537*IT_0656 +
       IT_0538*(IT_0657 + IT_0658) + IT_0539*((-16)*IT_0341 + IT_0659) + IT_0308
      *(IT_0660 + IT_0661) + IT_0302*(16*conj(IT_0351) + 16*conj(IT_0361) +
       IT_0660 + IT_0661) + IT_0541*(IT_0662 + IT_0663 + IT_0664 + IT_0665) +
       IT_0542*(IT_0666 + IT_0667 + IT_0668 + IT_0669) + IT_0521*(IT_0670 +
       IT_0671 + IT_0672) + IT_0525*(IT_0673 + IT_0674 + IT_0675);
    const ccomplex_t IT_0677 = pow(s_34, 2);
    const ccomplex_t IT_0678 = (-2)*IT_0677;
    const ccomplex_t IT_0679 = IT_0285 + IT_0678;
    const ccomplex_t IT_0680 = IT_0000*IT_0679;
    const ccomplex_t IT_0681 = 0.25*IT_0680;
    const ccomplex_t IT_0682 = 1 + IT_0681;
    const ccomplex_t IT_0683 = IT_0035*IT_0682;
    const ccomplex_t IT_0684 = 8*IT_0351;
    const ccomplex_t IT_0685 = 8*IT_0361;
    const ccomplex_t IT_0686 = (-8)*IT_0341;
    const ccomplex_t IT_0687 = (-8)*IT_0345;
    const ccomplex_t IT_0688 = IT_0684 + IT_0685 + IT_0686 + IT_0687;
    const ccomplex_t IT_0689 = (-8)*IT_0302;
    const ccomplex_t IT_0690 = (-8)*IT_0308;
    const ccomplex_t IT_0691 = IT_0689 + IT_0690;
    const ccomplex_t IT_0692 = (-8)*conj(IT_0302);
    const ccomplex_t IT_0693 = (-8)*conj(IT_0308);
    const ccomplex_t IT_0694 = IT_0692 + IT_0693;
    const ccomplex_t IT_0695 = 8*conj(IT_0341);
    const ccomplex_t IT_0696 = 8*IT_0406;
    const ccomplex_t IT_0697 = 8*conj(IT_0406);
    const ccomplex_t IT_0698 = 8*IT_0341;
    const ccomplex_t IT_0699 = 8*IT_0345;
    const ccomplex_t IT_0700 = 8*IT_0397;
    const ccomplex_t IT_0701 = 8*IT_0401;
    const ccomplex_t IT_0702 = (-8)*IT_0406;
    const ccomplex_t IT_0703 = 8*conj(IT_0397);
    const ccomplex_t IT_0704 = 8*conj(IT_0401);
    const ccomplex_t IT_0705 = (-8)*conj(IT_0406);
    const ccomplex_t IT_0706 = 8*conj(IT_0351);
    const ccomplex_t IT_0707 = 8*conj(IT_0361);
    const ccomplex_t IT_0708 = (-8)*conj(IT_0341);
    const ccomplex_t IT_0709 = (-8)*conj(IT_0345);
    const ccomplex_t IT_0710 = (-8)*IT_0428;
    const ccomplex_t IT_0711 = (-8)*conj(IT_0428);
    const ccomplex_t IT_0712 = IT_0521*IT_0688 + IT_0507*IT_0691 + IT_0511
      *IT_0694 + 8*IT_0428*(conj(IT_0345) + 0.125*IT_0695) + IT_0539*IT_0696 +
       IT_0540*IT_0697 + conj(IT_0428)*(IT_0698 + IT_0699) + IT_0541*(IT_0700 +
       IT_0701 + IT_0702) + IT_0542*(IT_0703 + IT_0704 + IT_0705) + IT_0525*
      (IT_0706 + IT_0707 + IT_0708 + IT_0709) + IT_0537*IT_0710 + IT_0538*IT_0711;
    const ccomplex_t IT_0713 = (-2)*IT_0611;
    const ccomplex_t IT_0714 = IT_0497 + IT_0713;
    const ccomplex_t IT_0715 = IT_0005*IT_0714;
    const ccomplex_t IT_0716 = (-0.5)*IT_0715;
    const ccomplex_t IT_0717 = (-0.5)*IT_0629;
    const ccomplex_t IT_0718 = s_12*IT_0285;
    const ccomplex_t IT_0719 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0720 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0721 = s_13*s_24*s_34;
    const ccomplex_t IT_0722 = (-2)*IT_0719 + (-2)*IT_0720 + 4*IT_0721;
    const ccomplex_t IT_0723 = IT_0718 + IT_0722;
    const ccomplex_t IT_0724 = IT_0000*IT_0723;
    const ccomplex_t IT_0725 = 0.25*IT_0724;
    const ccomplex_t IT_0726 = IT_0700 + IT_0702;
    const ccomplex_t IT_0727 = IT_0703 + IT_0705;
    const ccomplex_t IT_0728 = 8*IT_0416;
    const ccomplex_t IT_0729 = IT_0710 + IT_0728;
    const ccomplex_t IT_0730 = 8*conj(IT_0416);
    const ccomplex_t IT_0731 = IT_0711 + IT_0730;
    const ccomplex_t IT_0732 = 8*IT_0428;
    const ccomplex_t IT_0733 = 8*IT_0420;
    const ccomplex_t IT_0734 = 8*conj(IT_0420);
    const ccomplex_t IT_0735 = conj(IT_0406)*IT_0696 + conj(IT_0401)*(IT_0700 
      + IT_0701 + IT_0702) + IT_0397*IT_0705 + IT_0416*IT_0711 + conj(IT_0397)
      *IT_0726 + IT_0401*IT_0727 + IT_0521*IT_0729 + IT_0420*IT_0731 + conj
      (IT_0428)*IT_0732 + (conj(IT_0420) + conj(IT_0424))*IT_0733 + 8*IT_0424*
      (conj(IT_0424) + 0.125*IT_0711 + 0.125*IT_0730 + 0.125*IT_0734);
    const ccomplex_t IT_0736 = s_14*s_23*s_34;
    const ccomplex_t IT_0737 = (-2)*IT_0719 + (-2)*IT_0720 + 4*IT_0736;
    const ccomplex_t IT_0738 = IT_0718 + IT_0737;
    const ccomplex_t IT_0739 = IT_0000*IT_0738;
    const ccomplex_t IT_0740 = 0.25*IT_0739;
    const ccomplex_t IT_0741 = IT_0684 + IT_0686 + IT_0687;
    const ccomplex_t IT_0742 = IT_0706 + IT_0708 + IT_0709;
    const ccomplex_t IT_0743 = 8*IT_0315;
    const ccomplex_t IT_0744 = IT_0689 + IT_0690 + IT_0743;
    const ccomplex_t IT_0745 = 8*conj(IT_0315);
    const ccomplex_t IT_0746 = IT_0692 + IT_0693 + IT_0745;
    const ccomplex_t IT_0747 = 8*IT_0302;
    const ccomplex_t IT_0748 = 8*conj(IT_0302);
    const ccomplex_t IT_0749 = 8*IT_0326;
    const ccomplex_t IT_0750 = conj(IT_0361)*IT_0688 + IT_0315*IT_0694 +
       IT_0345*IT_0695 + IT_0535*IT_0698 + conj(IT_0345)*IT_0699 + IT_0351*
      (IT_0708 + IT_0709) + conj(IT_0351)*IT_0741 + IT_0361*IT_0742 + conj
      (IT_0315)*IT_0744 + IT_0326*IT_0746 + IT_0539*IT_0747 + 8*IT_0308*(conj
      (IT_0308) + 0.125*IT_0748) + conj(IT_0326)*(IT_0689 + IT_0690 + IT_0743 +
       IT_0749);
    const ccomplex_t IT_0751 = s_12*IT_0677;
    const ccomplex_t IT_0752 = (-2)*IT_0719 + (-2)*IT_0720 + 2*IT_0721 + 2
      *IT_0736 + (-2)*IT_0751;
    const ccomplex_t IT_0753 = IT_0718 + IT_0752;
    const ccomplex_t IT_0754 = IT_0000*IT_0753;
    const ccomplex_t IT_0755 = (-0.25)*IT_0754;
    const ccomplex_t IT_0756 = -IT_0218;
    const ccomplex_t IT_0757 = -IT_0217;
    const ccomplex_t IT_0758 = conj(IT_0401)*IT_0688 + conj(IT_0416)*IT_0691 +
       IT_0416*IT_0694 + IT_0406*IT_0695 + conj(IT_0345)*IT_0696 + IT_0345
      *IT_0697 + conj(IT_0406)*IT_0698 + conj(IT_0351)*IT_0702 + IT_0351*IT_0705
       + IT_0401*(IT_0706 + IT_0707 + IT_0708 + IT_0709) + conj(IT_0361)*IT_0726
       + IT_0361*IT_0727 + conj(IT_0315)*IT_0729 + IT_0315*IT_0731 + conj
      (IT_0308)*IT_0732 + 8*conj(IT_0326)*(IT_0424 + 0.125*IT_0710 + 0.125
      *IT_0728 + 0.125*IT_0733) + IT_0326*(IT_0711 + IT_0730 + IT_0734) + conj
      (IT_0397)*IT_0741 + IT_0397*IT_0742 + conj(IT_0420)*IT_0744 + (IT_0420 +
       IT_0424)*IT_0746 + 8*conj(IT_0428)*(IT_0308 + 0.125*IT_0747) + IT_0428
      *IT_0748 + conj(IT_0424)*(IT_0689 + IT_0690 + IT_0743 + IT_0749);
    const ccomplex_t IT_0759 = pow(m_Z, 6);
    const ccomplex_t IT_0760 = s_12*IT_0759;
    const ccomplex_t IT_0761 = s_13*s_23*IT_0285;
    const ccomplex_t IT_0762 = (-2)*IT_0761;
    const ccomplex_t IT_0763 = IT_0760 + IT_0762;
    const ccomplex_t IT_0764 = IT_0000*IT_0763;
    const ccomplex_t IT_0765 = 0.25*IT_0764;
    const ccomplex_t IT_0766 = IT_0497 + IT_0765;
    const ccomplex_t IT_0767 = (-8)*IT_0183;
    const ccomplex_t IT_0768 = IT_0186 + IT_0189;
    const ccomplex_t IT_0769 = (-8)*conj(IT_0183);
    const ccomplex_t IT_0770 = conj(IT_0186) + conj(IT_0189);
    const ccomplex_t IT_0771 = (-8)*IT_0185;
    const ccomplex_t IT_0772 = (-8)*conj(IT_0185);
    const ccomplex_t IT_0773 = (-8)*IT_0055;
    const ccomplex_t IT_0774 = (-8)*conj(IT_0055);
    const ccomplex_t IT_0775 = (-8)*IT_0059;
    const ccomplex_t IT_0776 = (-8)*conj(IT_0059);
    const ccomplex_t IT_0777 = (-8)*IT_0065;
    const ccomplex_t IT_0778 = IT_0034 + IT_0067 + IT_0109 + IT_0116;
    const ccomplex_t IT_0779 = (-8)*conj(IT_0065);
    const ccomplex_t IT_0780 = (-8)*IT_0034;
    const ccomplex_t IT_0781 = (-8)*IT_0067;
    const ccomplex_t IT_0782 = (-8)*conj(IT_0034);
    const ccomplex_t IT_0783 = (-8)*conj(IT_0067);
    const ccomplex_t IT_0784 = conj(IT_0109) + conj(IT_0116);
    const ccomplex_t IT_0785 = (-8)*IT_0087;
    const ccomplex_t IT_0786 = IT_0109 + IT_0116;
    const ccomplex_t IT_0787 = (-8)*conj(IT_0087);
    const ccomplex_t IT_0788 = conj(IT_0112) + conj(IT_0119);
    const ccomplex_t IT_0789 = (-8)*IT_0109;
    const ccomplex_t IT_0790 = IT_0112 + IT_0119;
    const ccomplex_t IT_0791 = (-8)*conj(IT_0109);
    const ccomplex_t IT_0792 = conj(IT_0185)*IT_0767 + (-8)*conj(IT_0187)*((
      -0.125)*IT_0767 + IT_0768) + IT_0185*IT_0769 + (-8)*IT_0187*((-0.125)
      *IT_0769 + IT_0770) + IT_0770*IT_0771 + IT_0768*IT_0772 + IT_0376*IT_0773 
      + IT_0378*IT_0774 + conj(IT_0065)*IT_0775 + IT_0065*IT_0776 + (conj
      (IT_0034) + conj(IT_0067) + conj(IT_0109) + conj(IT_0116))*IT_0777 +
       IT_0778*IT_0779 + conj(IT_0087)*(IT_0775 + IT_0780 + IT_0781) + IT_0087*
      (IT_0776 + IT_0782 + IT_0783) + IT_0784*IT_0785 + IT_0786*IT_0787 + (-8)
      *IT_0788*(IT_0116 + (-0.125)*IT_0775 + (-0.125)*IT_0780 + (-0.125)*IT_0781
       + (-0.125)*IT_0789) + (-8)*IT_0790*(conj(IT_0116) + (-0.125)*IT_0776 + (
      -0.125)*IT_0782 + (-0.125)*IT_0783 + (-0.125)*IT_0791);
    const ccomplex_t IT_0793 = s_12 + IT_0756;
    const ccomplex_t IT_0794 = IT_0035*IT_0793;
    const ccomplex_t IT_0795 = (-6)*conj(IT_0278);
    const ccomplex_t IT_0796 = (-6)*IT_0222;
    const ccomplex_t IT_0797 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0798 = IT_0000*IT_0677;
    const ccomplex_t IT_0799 = 0.5*IT_0798;
    const ccomplex_t IT_0800 = 1 + IT_0799;
    const ccomplex_t IT_0801 = s_12*IT_0800;
    const ccomplex_t IT_0802 = 4*IT_0177;
    const ccomplex_t IT_0803 = (-4)*IT_0175;
    const ccomplex_t IT_0804 = IT_0802 + IT_0803;
    const ccomplex_t IT_0805 = 4*IT_0209;
    const ccomplex_t IT_0806 = (-4)*IT_0207;
    const ccomplex_t IT_0807 = IT_0805 + IT_0806;
    const ccomplex_t IT_0808 = 4*IT_0175;
    const ccomplex_t IT_0809 = 4*IT_0207;
    const ccomplex_t IT_0810 = (-4)*conj(IT_0175);
    const ccomplex_t IT_0811 = (-4)*conj(IT_0207);
    const ccomplex_t IT_0812 = IT_0035*IT_0800;
    const ccomplex_t IT_0813 = 4*conj(IT_0177);
    const ccomplex_t IT_0814 = (-4)*IT_0177;
    const ccomplex_t IT_0815 = 4*conj(IT_0175);
    const ccomplex_t IT_0816 = (-4)*conj(IT_0177);
    const ccomplex_t IT_0817 = 4*conj(IT_0209);
    const ccomplex_t IT_0818 = (-4)*IT_0209;
    const ccomplex_t IT_0819 = 4*conj(IT_0207);
    const ccomplex_t IT_0820 = (-4)*conj(IT_0209);
    const ccomplex_t IT_0821 = (conj(IT_0207) + conj(IT_0302) + conj(IT_0308) 
      + conj(IT_0428))*IT_0804 + (conj(IT_0341) + conj(IT_0345) + conj(IT_0406))
      *IT_0807 + (IT_0207 + IT_0302 + IT_0308 + IT_0428)*(IT_0810 + IT_0813) + 
      (conj(IT_0209) + conj(IT_0315) + conj(IT_0326) + conj(IT_0416) + conj
      (IT_0420) + conj(IT_0424))*(IT_0808 + IT_0814) + (IT_0209 + IT_0315 +
       IT_0326 + IT_0416 + IT_0420 + IT_0424)*(IT_0815 + IT_0816) + (IT_0341 +
       IT_0345 + IT_0406)*(IT_0811 + IT_0817) + (conj(IT_0351) + conj(IT_0361) +
       conj(IT_0397) + conj(IT_0401))*(IT_0809 + IT_0818) + (IT_0351 + IT_0361 +
       IT_0397 + IT_0401)*(IT_0819 + IT_0820);
    const ccomplex_t IT_0822 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0823 = 0.5*IT_0268 + (-0.5)*IT_0590 + -IT_0822;
    const ccomplex_t IT_0824 = IT_0266 + IT_0823;
    const ccomplex_t IT_0825 = IT_0005*IT_0824;
    const ccomplex_t IT_0826 = IT_0035 + IT_0643 + IT_0825;
    const ccomplex_t IT_0827 = m_N_3*IT_0826;
    const ccomplex_t IT_0828 = (-4)*IT_0406;
    const ccomplex_t IT_0829 = (-4)*conj(IT_0406);
    const ccomplex_t IT_0830 = IT_0279*IT_0428 + IT_0280*conj(IT_0428) + (
      -0.25)*IT_0487*IT_0507 + (-0.25)*IT_0489*IT_0511 + (-0.25)*IT_0482*IT_0521
       + (-0.25)*IT_0485*IT_0525 + (-0.25)*IT_0276*IT_0828 + (-0.25)*IT_0277
      *IT_0829;
    const ccomplex_t IT_0831 = 0.5*IT_0764;
    const ccomplex_t IT_0832 = (-4)*IT_0611;
    const ccomplex_t IT_0833 = IT_0497 + IT_0831 + IT_0832;
    const ccomplex_t IT_0834 = IT_0005*IT_0738;
    const ccomplex_t IT_0835 = -IT_0834;
    const ccomplex_t IT_0836 = IT_0833 + IT_0835;
    const ccomplex_t IT_0837 = (-4)*IT_0034;
    const ccomplex_t IT_0838 = (-4)*IT_0059;
    const ccomplex_t IT_0839 = IT_0837 + IT_0838;
    const ccomplex_t IT_0840 = IT_0067 + IT_0109 + IT_0116;
    const ccomplex_t IT_0841 = (-4)*conj(IT_0034);
    const ccomplex_t IT_0842 = (-4)*IT_0185;
    const ccomplex_t IT_0843 = (-4)*conj(IT_0059);
    const ccomplex_t IT_0844 = (-4)*IT_0067;
    const ccomplex_t IT_0845 = (-4)*conj(IT_0067);
    const ccomplex_t IT_0846 = (-4)*IT_0109;
    const ccomplex_t IT_0847 = conj(IT_0059)*IT_0838 + conj(IT_0034)*IT_0839 +
       IT_0840*IT_0841 + (-4)*IT_0379*(IT_0187 + (-0.25)*IT_0842) + IT_0778
      *IT_0843 + conj(IT_0067)*(IT_0839 + IT_0844) + IT_0786*IT_0845 + (-4)
      *IT_0784*(IT_0116 + (-0.25)*IT_0837 + (-0.25)*IT_0838 + (-0.25)*IT_0844 + 
      (-0.25)*IT_0846);
    const ccomplex_t IT_0848 = IT_0005*IT_0723;
    const ccomplex_t IT_0849 = -IT_0848;
    const ccomplex_t IT_0850 = IT_0833 + IT_0849;
    const ccomplex_t IT_0851 = (-4)*IT_0183;
    const ccomplex_t IT_0852 = (-4)*conj(IT_0183);
    const ccomplex_t IT_0853 = (-4)*IT_0055;
    const ccomplex_t IT_0854 = (-4)*conj(IT_0055);
    const ccomplex_t IT_0855 = (-4)*IT_0065;
    const ccomplex_t IT_0856 = IT_0087 + IT_0112 + IT_0119;
    const ccomplex_t IT_0857 = (-4)*conj(IT_0065);
    const ccomplex_t IT_0858 = (-4)*IT_0087;
    const ccomplex_t IT_0859 = (-4)*conj(IT_0087);
    const ccomplex_t IT_0860 = (-4)*IT_0112;
    const ccomplex_t IT_0861 = (-4)*IT_0186;
    const ccomplex_t IT_0862 = conj(IT_0183)*IT_0851 + IT_0768*IT_0852 +
       IT_0296*IT_0853 + (IT_0065 + IT_0087 + IT_0112 + IT_0119)*IT_0854 + (conj
      (IT_0065) + conj(IT_0087))*IT_0855 + IT_0856*IT_0857 + conj(IT_0087)
      *IT_0858 + IT_0790*IT_0859 + (-4)*IT_0788*(IT_0119 + (-0.25)*IT_0853 + (
      -0.25)*IT_0855 + (-0.25)*IT_0858 + (-0.25)*IT_0860) + (-4)*IT_0770*
      (IT_0189 + (-0.25)*IT_0851 + (-0.25)*IT_0861);
    const ccomplex_t IT_0863 = IT_0000*IT_0286;
    const ccomplex_t IT_0864 = 0.666666666666667*IT_0242 + (-0.166666666666667
      )*IT_0863;
    const ccomplex_t IT_0865 = s_24 + IT_0864;
    const ccomplex_t IT_0866 = m_N_3*IT_0865;
    const ccomplex_t IT_0867 = (-6)*IT_0866;
    const ccomplex_t IT_0868 = (-10)*s_24;
    const ccomplex_t IT_0869 = IT_0863 + IT_0868;
    const ccomplex_t IT_0870 = m_N_3*IT_0869;
    const ccomplex_t IT_0871 = 0.666666666666667*IT_0007 + (-0.166666666666667
      )*IT_0566;
    const ccomplex_t IT_0872 = s_14 + IT_0871;
    const ccomplex_t IT_0873 = m_N_3*IT_0872;
    const ccomplex_t IT_0874 = (-6)*IT_0873;
    const ccomplex_t IT_0875 = (IT_0035 + -IT_0269)*(s_12 + (-0.5)*IT_0629);
    const ccomplex_t IT_0876 = 2*IT_0875;
    const ccomplex_t IT_0877 = IT_0233 + IT_0591;
    const ccomplex_t IT_0878 = IT_0005*IT_0877;
    const ccomplex_t IT_0879 = IT_0234 + IT_0594;
    const ccomplex_t IT_0880 = (-4)*IT_0879;
    const ccomplex_t IT_0881 = IT_0878 + IT_0880;
    const ccomplex_t IT_0882 = m_N_3*IT_0881;
    const ccomplex_t IT_0883 = 2*IT_0226;
    const ccomplex_t IT_0884 = IT_0624 + IT_0883;
    const ccomplex_t IT_0885 = 2*conj(IT_0226);
    const ccomplex_t IT_0886 = IT_0625 + IT_0885;
    const ccomplex_t IT_0887 = 2*IT_0223;
    const ccomplex_t IT_0888 = IT_0624 + IT_0887;
    const ccomplex_t IT_0889 = 2*conj(IT_0223);
    const ccomplex_t IT_0890 = IT_0625 + IT_0889;
    const ccomplex_t IT_0891 = conj(IT_0226)*IT_0407 + IT_0226*IT_0412 + conj
      (IT_0223)*IT_0429 + IT_0223*IT_0435 + conj(IT_0406)*IT_0580 + conj(IT_0428
      )*IT_0580 + IT_0406*IT_0583 + IT_0428*IT_0583 + IT_0507*IT_0884 + IT_0511
      *IT_0886 + IT_0521*IT_0888 + IT_0525*IT_0890;
    const ccomplex_t IT_0892 = (IT_0035 + -IT_0234)*(s_12 + (-0.5)*IT_0629);
    const ccomplex_t IT_0893 = 2*IT_0892;
    const ccomplex_t IT_0894 = 9*IT_0590;
    const ccomplex_t IT_0895 = (-2)*IT_0034;
    const ccomplex_t IT_0896 = (-2)*IT_0059;
    const ccomplex_t IT_0897 = (-2)*IT_0067;
    const ccomplex_t IT_0898 = (-2)*IT_0109;
    const ccomplex_t IT_0899 = (-2)*IT_0116;
    const ccomplex_t IT_0900 = IT_0895 + IT_0896 + IT_0897 + IT_0898 + IT_0899;
    const ccomplex_t IT_0901 = (-2)*conj(IT_0034);
    const ccomplex_t IT_0902 = (-2)*conj(IT_0059);
    const ccomplex_t IT_0903 = (-2)*conj(IT_0067);
    const ccomplex_t IT_0904 = (-2)*conj(IT_0109);
    const ccomplex_t IT_0905 = (-2)*conj(IT_0116);
    const ccomplex_t IT_0906 = IT_0901 + IT_0902 + IT_0903 + IT_0904 + IT_0905;
    const ccomplex_t IT_0907 = (-2)*IT_0055;
    const ccomplex_t IT_0908 = (-2)*IT_0065;
    const ccomplex_t IT_0909 = (-2)*IT_0087;
    const ccomplex_t IT_0910 = (-2)*IT_0112;
    const ccomplex_t IT_0911 = (-2)*IT_0119;
    const ccomplex_t IT_0912 = IT_0907 + IT_0908 + IT_0909 + IT_0910 + IT_0911;
    const ccomplex_t IT_0913 = (-2)*conj(IT_0055);
    const ccomplex_t IT_0914 = (-2)*conj(IT_0065);
    const ccomplex_t IT_0915 = (-2)*conj(IT_0087);
    const ccomplex_t IT_0916 = (-2)*conj(IT_0112);
    const ccomplex_t IT_0917 = (-2)*conj(IT_0119);
    const ccomplex_t IT_0918 = IT_0913 + IT_0914 + IT_0915 + IT_0916 + IT_0917;
    const ccomplex_t IT_0919 = (-2)*IT_0183;
    const ccomplex_t IT_0920 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0921 = (-2)*IT_0185;
    const ccomplex_t IT_0922 = (-2)*conj(IT_0185);
    const ccomplex_t IT_0923 = IT_0005*IT_0035;
    const ccomplex_t IT_0924 = (-0.25)*IT_0285;
    const ccomplex_t IT_0925 = IT_0677 + IT_0924;
    const ccomplex_t IT_0926 = IT_0923*IT_0925;
    const ccomplex_t IT_0927 = (-4)*IT_0926;
    const ccomplex_t IT_0928 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0929 = (-2)*IT_0610 + -IT_0611 + 2*IT_0928;
    const ccomplex_t IT_0930 = IT_0497 + IT_0929;
    const ccomplex_t IT_0931 = IT_0005*IT_0930;
    const ccomplex_t IT_0932 = -IT_0218 + (-0.5)*IT_0931;
    const ccomplex_t IT_0933 = s_12 + IT_0932;
    const ccomplex_t IT_0934 = m_N_3*IT_0933;
    const ccomplex_t IT_0935 = (-2)*IT_0934;
    const ccomplex_t IT_0936 = (-2)*conj(IT_0302);
    const ccomplex_t IT_0937 = (-2)*conj(IT_0308);
    const ccomplex_t IT_0938 = IT_0936 + IT_0937;
    const ccomplex_t IT_0939 = (-2)*IT_0341;
    const ccomplex_t IT_0940 = (-2)*IT_0345;
    const ccomplex_t IT_0941 = IT_0939 + IT_0940;
    const ccomplex_t IT_0942 = (-2)*conj(IT_0341);
    const ccomplex_t IT_0943 = (-2)*conj(IT_0345);
    const ccomplex_t IT_0944 = IT_0942 + IT_0943;
    const ccomplex_t IT_0945 = (-2)*IT_0308;
    const ccomplex_t IT_0946 = (-2)*IT_0302;
    const ccomplex_t IT_0947 = IT_0464*IT_0537 + IT_0467*IT_0538 + IT_0469
      *IT_0541 + IT_0471*IT_0542 + IT_0535*IT_0599 + IT_0536*IT_0602 + IT_0539
      *IT_0605 + IT_0308*IT_0607 + IT_0275*IT_0938 + conj(IT_0275)*IT_0941 +
       IT_0275*IT_0944 + conj(IT_0275)*IT_0945 + (conj(IT_0275) + (-0.5)*IT_0607
      )*IT_0946;
    const ccomplex_t IT_0948 = -IT_0611;
    const ccomplex_t IT_0949 = IT_0497 + IT_0948;
    const ccomplex_t IT_0950 = IT_0005*IT_0949;
    const ccomplex_t IT_0951 = -IT_0218 + (-0.25)*IT_0950;
    const ccomplex_t IT_0952 = s_12 + IT_0951;
    const ccomplex_t IT_0953 = m_N_3*IT_0952;
    const ccomplex_t IT_0954 = (-4)*IT_0953;
    const ccomplex_t IT_0955 = conj(IT_0406)*IT_0462 + IT_0406*IT_0465 + conj
      (IT_0428)*IT_0468 + IT_0428*IT_0470 + IT_0539*IT_0580 + IT_0308*IT_0583 +
       IT_0507*IT_0601 + IT_0511*IT_0604 + IT_0521*IT_0606 + IT_0525*IT_0608 +
       conj(IT_0275)*IT_0620 + IT_0275*IT_0621 + conj(IT_0275)*IT_0622 + IT_0275
      *IT_0623 + IT_0541*IT_0884 + IT_0542*IT_0886 + IT_0537*IT_0888 + IT_0538
      *IT_0890 + IT_0226*IT_0938 + IT_0224*IT_0941 + IT_0225*IT_0944 + conj
      (IT_0226)*IT_0945 + (conj(IT_0226) + (-0.5)*IT_0583)*IT_0946;
    const ccomplex_t IT_0956 = IT_0015*IT_0215 + IT_0221*(IT_0179*IT_0224 +
       IT_0182*IT_0225 + IT_0211*IT_0227 + IT_0214*IT_0228) + (IT_0211*IT_0224 +
       IT_0214*IT_0225 + IT_0179*IT_0227 + IT_0182*IT_0228)*IT_0238 + IT_0250
      *IT_0263 + IT_0273*(IT_0253*IT_0276 + IT_0256*IT_0277 + IT_0259*IT_0279 +
       IT_0262*IT_0280) + (IT_0259*IT_0276 + IT_0262*IT_0277 + IT_0253*IT_0279 +
       IT_0256*IT_0280)*IT_0283 + IT_0294*(IT_0297*IT_0328 + IT_0330*IT_0336 +
       IT_0337*IT_0363 + IT_0364*IT_0370) + IT_0374*(IT_0328*IT_0376 + IT_0336
      *IT_0378 + IT_0363*IT_0379 + IT_0370*IT_0380) + (IT_0328*IT_0337 + IT_0297
      *IT_0363 + IT_0336*IT_0364 + IT_0330*IT_0370)*IT_0387 + IT_0392*(IT_0379
      *IT_0408 + IT_0380*IT_0413 + IT_0376*IT_0430 + IT_0378*IT_0436) + (IT_0297
      *IT_0408 + IT_0330*IT_0413 + IT_0337*IT_0430 + IT_0364*IT_0436)*IT_0444 + 
      (IT_0376*IT_0408 + IT_0378*IT_0413 + IT_0379*IT_0430 + IT_0380*IT_0436)
      *IT_0449 + IT_0461*(IT_0337*IT_0464 + IT_0364*IT_0467 + IT_0297*IT_0469 +
       IT_0330*IT_0471) + IT_0479*(IT_0379*IT_0482 + IT_0380*IT_0485 + IT_0376
      *IT_0487 + IT_0378*IT_0489) + IT_0500*(conj(IT_0406)*IT_0503 + IT_0406
      *IT_0506 + IT_0507*IT_0510 + IT_0511*IT_0514 + conj(IT_0428)*IT_0517 +
       IT_0428*IT_0520 + IT_0521*IT_0524 + IT_0525*IT_0528) + IT_0534*(IT_0503
      *IT_0535 + IT_0506*IT_0536 + IT_0510*IT_0537 + IT_0514*IT_0538 + IT_0517
      *IT_0539 + IT_0520*IT_0540 + IT_0524*IT_0541 + IT_0528*IT_0542) + IT_0555*
      (IT_0379*IT_0558 + IT_0380*IT_0561 + IT_0376*IT_0563 + IT_0378*IT_0565) +
       IT_0569*(conj(IT_0187)*IT_0310 + conj(IT_0185)*IT_0328 + IT_0187*IT_0333 
      + IT_0185*IT_0336 + IT_0347*IT_0376 + IT_0367*IT_0378 + IT_0375*IT_0570 +
       IT_0377*IT_0571) + IT_0579*(IT_0337*IT_0582 + IT_0364*IT_0585 + IT_0297
      *IT_0587 + IT_0330*IT_0589) + IT_0609 + IT_0626 + IT_0631*IT_0632 +
       IT_0639*IT_0642 + 6*IT_0646*(IT_0276*IT_0280 + IT_0277*(conj(IT_0275) +
       0.166666666666667*IT_0647)) + 6*IT_0650*(IT_0224*(IT_0226 +
       0.166666666666667*IT_0651) + IT_0225*(conj(IT_0226) + 0.166666666666667
      *IT_0652)) + IT_0653*IT_0676 + IT_0683*IT_0712 + (s_12 + IT_0716 + IT_0717
       + IT_0725)*IT_0735 + (s_12 + IT_0716 + IT_0717 + IT_0740)*IT_0750 + (s_12
       + IT_0755 + IT_0756 + IT_0757)*IT_0758 + IT_0766*IT_0792 + (-6)*IT_0794*
      (IT_0224*IT_0280 + IT_0225*(conj(IT_0275) + (-0.166666666666667)*IT_0795) 
      + IT_0276*(IT_0226 + (-0.166666666666667)*IT_0796) + IT_0277*(conj(IT_0226
      ) + (-0.166666666666667)*IT_0797)) + IT_0801*(conj(IT_0177)*IT_0804 + conj
      (IT_0209)*IT_0807 + conj(IT_0175)*IT_0808 + conj(IT_0207)*IT_0809 +
       IT_0177*IT_0810 + IT_0209*IT_0811) + IT_0812*IT_0821 + (-4)*IT_0827
      *IT_0830 + IT_0836*IT_0847 + IT_0850*IT_0862 + (conj(IT_0183)*IT_0407 +
       IT_0183*IT_0412 + IT_0295*IT_0429 + IT_0329*IT_0435)*IT_0867 + (IT_0337
      *IT_0403 + IT_0364*IT_0411 + IT_0297*IT_0426 + (conj(IT_0087) + conj
      (IT_0112) + conj(IT_0119))*IT_0429 + IT_0330*IT_0434 + IT_0412*IT_0768 +
       IT_0407*IT_0770 + IT_0435*IT_0856)*IT_0870 + (conj(IT_0187)*(IT_0316 +
       IT_0327) + IT_0187*(IT_0334 + IT_0335) + (conj(IT_0067) + conj(IT_0109) +
       conj(IT_0116))*IT_0570 + IT_0571*IT_0840)*IT_0874 + (IT_0279*IT_0599 + 2
      *conj(IT_0275)*IT_0600 + IT_0275*IT_0602 + IT_0276*IT_0605 + IT_0275
      *IT_0607)*IT_0876 + IT_0882*IT_0891 + (IT_0224*IT_0580 + IT_0227*IT_0580 +
       conj(IT_0226)*IT_0581 + IT_0223*IT_0583 + IT_0226*IT_0583 + conj(IT_0223)
      *IT_0586)*IT_0893 + IT_0894*(conj(IT_0187)*IT_0900 + IT_0187*IT_0906 +
       IT_0770*IT_0912 + IT_0768*IT_0918 + IT_0297*IT_0919 + IT_0330*IT_0920 +
       IT_0376*IT_0921 + IT_0378*IT_0922) + (IT_0770*IT_0900 + IT_0768*IT_0906 +
       conj(IT_0187)*IT_0912 + IT_0187*IT_0918 + IT_0376*IT_0919 + IT_0378
      *IT_0920 + IT_0297*IT_0921 + IT_0330*IT_0922)*IT_0927 + IT_0935*IT_0947 +
       IT_0954*IT_0955;
    return create_ccomplex_return(IT_0956);
}

