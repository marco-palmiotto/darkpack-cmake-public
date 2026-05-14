#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sd_L_to_d_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sd_L_to_d_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_d*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0004 + 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0013 = (-0.333333333333333)*IT_0012;
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = pow(m_d, 2);
    const ccomplex_t IT_0016 = pow(m_N_4, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0015 + IT_0016 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0001*IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_sd_L, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0015 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0011*IT_0020;
    const ccomplex_t IT_0029 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = IT_0005*IT_0011;
    const ccomplex_t IT_0032 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*IT_0020;
    const ccomplex_t IT_0035 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0001*IT_0005;
    const ccomplex_t IT_0038 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0041 = IT_0028*IT_0040;
    const ccomplex_t IT_0042 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0043 = IT_0031*IT_0042;
    const ccomplex_t IT_0044 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0045 = IT_0034*IT_0044;
    const ccomplex_t IT_0046 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0047 = IT_0037*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0030 + 0.5*IT_0033 + 0.5*IT_0036 + (
      -0.5)*IT_0039 + (-0.5)*IT_0041 + 0.5*IT_0043 + 0.5*IT_0045 + (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = m_N_2*IT_0048;
    const ccomplex_t IT_0050 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0051 = IT_0028*IT_0050;
    const ccomplex_t IT_0052 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0053 = IT_0031*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0055 = IT_0034*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0057 = IT_0037*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0059 = IT_0028*IT_0058;
    const ccomplex_t IT_0060 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0061 = IT_0031*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0063 = IT_0034*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0065 = IT_0037*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0051 + (-0.5)*IT_0053 + (-0.5)*IT_0055 +
       0.5*IT_0057 + 0.5*IT_0059 + (-0.5)*IT_0061 + (-0.5)*IT_0063 + 0.5*IT_0065;
    const ccomplex_t IT_0067 = m_N_4*IT_0066;
    const ccomplex_t IT_0068 = IT_0049 + IT_0067;
    const ccomplex_t IT_0069 = IT_0027*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0001*IT_0002*IT_0021;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0015 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0076 = IT_0028*IT_0075;
    const ccomplex_t IT_0077 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0078 = IT_0031*IT_0077;
    const ccomplex_t IT_0079 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0080 = IT_0034*IT_0079;
    const ccomplex_t IT_0081 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0082 = IT_0037*IT_0081;
    const ccomplex_t IT_0083 = -IT_0076 + IT_0078 + IT_0080 + -IT_0082;
    const ccomplex_t IT_0084 = m_N_4*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0086 = IT_0028*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0088 = IT_0031*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0090 = IT_0034*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0092 = IT_0037*IT_0091;
    const ccomplex_t IT_0093 = IT_0086 + -IT_0088 + -IT_0090 + IT_0092;
    const ccomplex_t IT_0094 = m_N_4*IT_0093;
    const ccomplex_t IT_0095 = IT_0084 + IT_0094;
    const ccomplex_t IT_0096 = IT_0074*IT_0095;
    const ccomplex_t IT_0097 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0098 = IT_0028*IT_0097;
    const ccomplex_t IT_0099 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0100 = IT_0031*IT_0099;
    const ccomplex_t IT_0101 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0102 = IT_0034*IT_0101;
    const ccomplex_t IT_0103 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0104 = IT_0037*IT_0103;
    const ccomplex_t IT_0105 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0106 = IT_0028*IT_0105;
    const ccomplex_t IT_0107 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0108 = IT_0031*IT_0107;
    const ccomplex_t IT_0109 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0110 = IT_0034*IT_0109;
    const ccomplex_t IT_0111 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0112 = IT_0037*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0098 + 0.5*IT_0100 + 0.5*IT_0102 + (
      -0.5)*IT_0104 + (-0.5)*IT_0106 + 0.5*IT_0108 + 0.5*IT_0110 + (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = m_N_1*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0116 = IT_0028*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0118 = IT_0031*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0120 = IT_0034*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0122 = IT_0037*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0124 = IT_0028*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0126 = IT_0031*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0128 = IT_0034*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0130 = IT_0037*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0116 + (-0.5)*IT_0118 + (-0.5)*IT_0120 +
       0.5*IT_0122 + 0.5*IT_0124 + (-0.5)*IT_0126 + (-0.5)*IT_0128 + 0.5*IT_0130;
    const ccomplex_t IT_0132 = m_N_4*IT_0131;
    const ccomplex_t IT_0133 = IT_0114 + IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0001*IT_0002*IT_0021;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_23 + IT_0015 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0133*IT_0138;
    const ccomplex_t IT_0140 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0141 = IT_0028*IT_0140;
    const ccomplex_t IT_0142 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0143 = IT_0031*IT_0142;
    const ccomplex_t IT_0144 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0145 = IT_0034*IT_0144;
    const ccomplex_t IT_0146 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0147 = IT_0037*IT_0146;
    const ccomplex_t IT_0148 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0149 = IT_0028*IT_0148;
    const ccomplex_t IT_0150 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0151 = IT_0031*IT_0150;
    const ccomplex_t IT_0152 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0153 = IT_0034*IT_0152;
    const ccomplex_t IT_0154 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0155 = IT_0037*IT_0154;
    const ccomplex_t IT_0156 = (-0.5)*IT_0141 + 0.5*IT_0143 + 0.5*IT_0145 + (
      -0.5)*IT_0147 + (-0.5)*IT_0149 + 0.5*IT_0151 + 0.5*IT_0153 + (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = m_N_3*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0159 = IT_0028*IT_0158;
    const ccomplex_t IT_0160 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0161 = IT_0031*IT_0160;
    const ccomplex_t IT_0162 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0163 = IT_0034*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0165 = IT_0037*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0167 = IT_0028*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0169 = IT_0031*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0171 = IT_0034*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0173 = IT_0037*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0159 + (-0.5)*IT_0161 + (-0.5)*IT_0163 +
       0.5*IT_0165 + 0.5*IT_0167 + (-0.5)*IT_0169 + (-0.5)*IT_0171 + 0.5*IT_0173;
    const ccomplex_t IT_0175 = m_N_4*IT_0174;
    const ccomplex_t IT_0176 = IT_0157 + IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0001*IT_0002*IT_0021;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = cpow((-2)*s_23 + IT_0015 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0180 = IT_0178*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = IT_0176*IT_0181;
    const ccomplex_t IT_0183 = IT_0019 + IT_0069 + IT_0096 + IT_0139 + IT_0182;
    const ccomplex_t IT_0184 = conj(N_B3)*e_em;
    const ccomplex_t IT_0185 = IT_0011*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W3)*e_em;
    const ccomplex_t IT_0188 = IT_0001*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0186 + (-3)*IT_0189);
    const ccomplex_t IT_0191 = 0.166666666666667*IT_0190;
    const ccomplex_t IT_0192 = IT_0179*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_N_4*IT_0156;
    const ccomplex_t IT_0195 = m_N_3*IT_0174;
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = IT_0193*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0001*IT_0002*IT_0021;
    const ccomplex_t IT_0199 = 0.5*IT_0198;
    const ccomplex_t IT_0200 = IT_0009*IT_0199;
    const ccomplex_t IT_0201 = IT_0017*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = m_N_4*IT_0113;
    const ccomplex_t IT_0204 = m_N_1*IT_0131;
    const ccomplex_t IT_0205 = IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = conj(N_B1)*e_em;
    const ccomplex_t IT_0207 = IT_0011*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W1)*e_em;
    const ccomplex_t IT_0210 = IT_0001*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*(IT_0208 + (-3)*IT_0211);
    const ccomplex_t IT_0213 = 0.166666666666667*IT_0212;
    const ccomplex_t IT_0214 = IT_0136*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0205*IT_0215;
    const ccomplex_t IT_0217 = m_N_4*IT_0048;
    const ccomplex_t IT_0218 = m_N_2*IT_0066;
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = conj(N_B2)*e_em;
    const ccomplex_t IT_0221 = IT_0011*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W2)*e_em;
    const ccomplex_t IT_0224 = IT_0001*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + (-3)*IT_0225);
    const ccomplex_t IT_0227 = 0.166666666666667*IT_0226;
    const ccomplex_t IT_0228 = IT_0025*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = IT_0219*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B4)*e_em;
    const ccomplex_t IT_0232 = IT_0011*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W4)*e_em;
    const ccomplex_t IT_0235 = IT_0001*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0233 + (-3)*IT_0236);
    const ccomplex_t IT_0238 = 0.166666666666667*IT_0237;
    const ccomplex_t IT_0239 = IT_0072*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = IT_0095*IT_0240;
    const ccomplex_t IT_0242 = -IT_0197 + -IT_0202 + -IT_0216 + -IT_0230 + 
      -IT_0241;
    const ccomplex_t IT_0243 = 6*s_13;
    const ccomplex_t IT_0244 = m_d*m_N_4;
    const ccomplex_t IT_0245 = 6*IT_0183;
    const ccomplex_t IT_0246 = m_d*e_em*IT_0001*IT_0002*IT_0005*IT_0021;
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0071*IT_0247;
    const ccomplex_t IT_0249 = cpow((-2)*s_12 + IT_0015 + -IT_0016 + -IT_0024 
      + -reg_prop, -1);
    const ccomplex_t IT_0250 = IT_0248*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0113*IT_0215;
    const ccomplex_t IT_0253 = IT_0156*IT_0193;
    const ccomplex_t IT_0254 = IT_0083*IT_0240;
    const ccomplex_t IT_0255 = IT_0048*IT_0229;
    const ccomplex_t IT_0256 = -IT_0251 + IT_0252 + IT_0253 + IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = m_d*s_14;
    const ccomplex_t IT_0258 = IT_0238*IT_0247;
    const ccomplex_t IT_0259 = IT_0249*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = IT_0074*IT_0093;
    const ccomplex_t IT_0262 = IT_0131*IT_0138;
    const ccomplex_t IT_0263 = IT_0174*IT_0181;
    const ccomplex_t IT_0264 = IT_0027*IT_0066;
    const ccomplex_t IT_0265 = -IT_0260 + -IT_0261 + -IT_0262 + -IT_0263 + 
      -IT_0264;
    const ccomplex_t IT_0266 = s_34*m_N_4;
    const ccomplex_t IT_0267 = 6*conj(IT_0183);
    const ccomplex_t IT_0268 = 6*conj(IT_0242);
    const ccomplex_t IT_0269 = pow(m_A0, 2);
    const ccomplex_t IT_0270 = m_d*m_N_4*IT_0269;
    const ccomplex_t IT_0271 = s_13*IT_0269;
    const ccomplex_t IT_0272 = s_14*s_34;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = 6*IT_0242;
    const ccomplex_t IT_0276 = 0.166666666666667*IT_0266;
    const ccomplex_t IT_0277 = (IT_0183*conj(IT_0183) + IT_0242*conj(IT_0242))
      *IT_0243 + conj(IT_0242)*IT_0244*IT_0245 + (IT_0242*IT_0244 + IT_0256
      *IT_0257 + IT_0265*IT_0266)*IT_0267 + (IT_0257*IT_0265 + IT_0256*IT_0266)
      *IT_0268 + 6*conj(IT_0265)*(IT_0256*IT_0270 + -IT_0265*IT_0274 +
       0.166666666666667*IT_0257*IT_0275 + IT_0245*IT_0276) + 6*conj(IT_0256)*
      (0.166666666666667*IT_0245*IT_0257 + IT_0265*IT_0270 + -IT_0256*IT_0274 +
       IT_0275*IT_0276);
    return create_ccomplex_return(IT_0277);
}

