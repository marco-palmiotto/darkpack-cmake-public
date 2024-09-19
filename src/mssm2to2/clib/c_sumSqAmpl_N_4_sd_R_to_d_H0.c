#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sd_R_to_d_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sd_R_to_d_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0005 + -IT_0009 + 
      -IT_0014 + IT_0017);
    const ccomplex_t IT_0019 = m_N_4*IT_0018;
    const ccomplex_t IT_0020 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0023 = IT_0007*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023 + 
      -IT_0025 + IT_0027);
    const ccomplex_t IT_0029 = m_N_4*IT_0028;
    const ccomplex_t IT_0030 = IT_0019 + IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0002;
    const ccomplex_t IT_0032 = (-0.333333333333333)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_d, 2);
    const ccomplex_t IT_0034 = pow(m_sd_R, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = pow(m_W, -1);
    const ccomplex_t IT_0040 = cos(beta);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0000 + -m_d*e_em*mu_h*IT_0006*IT_0011*IT_0039*IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0011*IT_0039*IT_0041;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = pow(m_N_4, 2);
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0033 + IT_0047 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0002;
    const ccomplex_t IT_0052 = (-0.333333333333333)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0063 = IT_0015*IT_0062;
    const ccomplex_t IT_0064 = IT_0057 + IT_0059 + IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0066 = IT_0007*IT_0065;
    const ccomplex_t IT_0067 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0068 = IT_0012*IT_0067;
    const ccomplex_t IT_0069 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0070 = IT_0007*IT_0069;
    const ccomplex_t IT_0071 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0072 = IT_0012*IT_0071;
    const ccomplex_t IT_0073 = -IT_0066 + -IT_0068 + -IT_0070 + -IT_0072;
    const ccomplex_t IT_0074 = IT_0064 + IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = m_N_2*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0079 = IT_0003*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0081 = IT_0015*IT_0080;
    const ccomplex_t IT_0082 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0085 = IT_0015*IT_0084;
    const ccomplex_t IT_0086 = IT_0079 + IT_0081 + IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0088 = IT_0007*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0090 = IT_0012*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0092 = IT_0007*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0094 = IT_0012*IT_0093;
    const ccomplex_t IT_0095 = -IT_0088 + -IT_0090 + -IT_0092 + -IT_0094;
    const ccomplex_t IT_0096 = IT_0086 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = m_N_4*IT_0098;
    const ccomplex_t IT_0100 = IT_0077 + IT_0099;
    const ccomplex_t IT_0101 = IT_0055*IT_0100;
    const ccomplex_t IT_0102 = cpow(IT_0001, -2);
    const ccomplex_t IT_0103 = sin(beta);
    const ccomplex_t IT_0104 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0011
      *IT_0033*IT_0039*IT_0041 + -1./3*m_W*IT_0010*IT_0102*(IT_0000*IT_0040 + 
      -IT_0006*IT_0103));
    const ccomplex_t IT_0105 = 0.333333333333333*IT_0104;
    const ccomplex_t IT_0106 = IT_0032*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0033 + IT_0047 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0011
      *IT_0039*IT_0041;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0032*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_12 + IT_0033 + -IT_0034 + -IT_0047 
      + -reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = m_d*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0002;
    const ccomplex_t IT_0118 = (-0.333333333333333)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0123 = IT_0003*IT_0122;
    const ccomplex_t IT_0124 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0125 = IT_0015*IT_0124;
    const ccomplex_t IT_0126 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0127 = IT_0003*IT_0126;
    const ccomplex_t IT_0128 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0129 = IT_0015*IT_0128;
    const ccomplex_t IT_0130 = IT_0123 + IT_0125 + IT_0127 + IT_0129;
    const ccomplex_t IT_0131 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0132 = IT_0007*IT_0131;
    const ccomplex_t IT_0133 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0134 = IT_0012*IT_0133;
    const ccomplex_t IT_0135 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0136 = IT_0007*IT_0135;
    const ccomplex_t IT_0137 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0138 = IT_0012*IT_0137;
    const ccomplex_t IT_0139 = -IT_0132 + -IT_0134 + -IT_0136 + -IT_0138;
    const ccomplex_t IT_0140 = IT_0130 + IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = m_N_1*IT_0142;
    const ccomplex_t IT_0144 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0145 = IT_0003*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0147 = IT_0015*IT_0146;
    const ccomplex_t IT_0148 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0149 = IT_0003*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0151 = IT_0015*IT_0150;
    const ccomplex_t IT_0152 = IT_0145 + IT_0147 + IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0154 = IT_0007*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0156 = IT_0012*IT_0155;
    const ccomplex_t IT_0157 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0158 = IT_0007*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0160 = IT_0012*IT_0159;
    const ccomplex_t IT_0161 = -IT_0154 + -IT_0156 + -IT_0158 + -IT_0160;
    const ccomplex_t IT_0162 = IT_0152 + IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = m_N_4*IT_0164;
    const ccomplex_t IT_0166 = IT_0143 + IT_0165;
    const ccomplex_t IT_0167 = IT_0121*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0002;
    const ccomplex_t IT_0169 = (-0.333333333333333)*IT_0168;
    const ccomplex_t IT_0170 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0174 = IT_0003*IT_0173;
    const ccomplex_t IT_0175 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0176 = IT_0015*IT_0175;
    const ccomplex_t IT_0177 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0178 = IT_0003*IT_0177;
    const ccomplex_t IT_0179 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0180 = IT_0015*IT_0179;
    const ccomplex_t IT_0181 = IT_0174 + IT_0176 + IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0183 = IT_0007*IT_0182;
    const ccomplex_t IT_0184 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0185 = IT_0012*IT_0184;
    const ccomplex_t IT_0186 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0187 = IT_0007*IT_0186;
    const ccomplex_t IT_0188 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0189 = IT_0012*IT_0188;
    const ccomplex_t IT_0190 = -IT_0183 + -IT_0185 + -IT_0187 + -IT_0189;
    const ccomplex_t IT_0191 = IT_0181 + IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 0.5*IT_0192;
    const ccomplex_t IT_0194 = m_N_3*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0196 = IT_0003*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0198 = IT_0015*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0200 = IT_0003*IT_0199;
    const ccomplex_t IT_0201 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0202 = IT_0015*IT_0201;
    const ccomplex_t IT_0203 = IT_0196 + IT_0198 + IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0205 = IT_0007*IT_0204;
    const ccomplex_t IT_0206 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0207 = IT_0012*IT_0206;
    const ccomplex_t IT_0208 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0209 = IT_0007*IT_0208;
    const ccomplex_t IT_0210 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0211 = IT_0012*IT_0210;
    const ccomplex_t IT_0212 = -IT_0205 + -IT_0207 + -IT_0209 + -IT_0211;
    const ccomplex_t IT_0213 = IT_0203 + IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = m_N_4*IT_0215;
    const ccomplex_t IT_0217 = IT_0194 + IT_0216;
    const ccomplex_t IT_0218 = IT_0172*IT_0217;
    const ccomplex_t IT_0219 = IT_0038 + IT_0050 + IT_0101 + IT_0109 + (-2)
      *IT_0116 + IT_0167 + IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0011*IT_0039*IT_0041;
    const ccomplex_t IT_0221 = 0.5*IT_0220;
    const ccomplex_t IT_0222 = IT_0111*IT_0221;
    const ccomplex_t IT_0223 = IT_0113*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = m_d*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0011*IT_0039*IT_0041;
    const ccomplex_t IT_0227 = 0.5*IT_0226;
    const ccomplex_t IT_0228 = IT_0119*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = m_N_4*IT_0142;
    const ccomplex_t IT_0231 = m_N_1*IT_0164;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0229*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0011*IT_0039*IT_0041;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = IT_0053*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = m_N_4*IT_0076;
    const ccomplex_t IT_0239 = m_N_2*IT_0098;
    const ccomplex_t IT_0240 = IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = IT_0237*IT_0240;
    const ccomplex_t IT_0242 = IT_0035*IT_0221;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = IT_0030*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0011*IT_0039*IT_0041;
    const ccomplex_t IT_0246 = 0.5*IT_0245;
    const ccomplex_t IT_0247 = IT_0170*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = m_N_4*IT_0193;
    const ccomplex_t IT_0250 = m_N_3*IT_0215;
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = IT_0248*IT_0251;
    const ccomplex_t IT_0253 = conj(N_B4)*e_em;
    const ccomplex_t IT_0254 = IT_0002*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = conj(N_W4)*e_em;
    const ccomplex_t IT_0257 = IT_0011*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*(IT_0255 + (-3)*IT_0258);
    const ccomplex_t IT_0260 = 0.166666666666667*IT_0259;
    const ccomplex_t IT_0261 = IT_0043*IT_0260;
    const ccomplex_t IT_0262 = IT_0048*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = IT_0105*IT_0221;
    const ccomplex_t IT_0265 = IT_0107*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = 2*IT_0225 + -IT_0233 + -IT_0241 + -IT_0244 + 
      -IT_0252 + -IT_0263 + -IT_0266;
    const ccomplex_t IT_0268 = 6*s_13;
    const ccomplex_t IT_0269 = m_d*m_N_4;
    const ccomplex_t IT_0270 = 6*IT_0219;
    const ccomplex_t IT_0271 = IT_0018*IT_0243;
    const ccomplex_t IT_0272 = IT_0193*IT_0248;
    const ccomplex_t IT_0273 = IT_0142*IT_0229;
    const ccomplex_t IT_0274 = IT_0076*IT_0237;
    const ccomplex_t IT_0275 = -IT_0115 + IT_0271 + IT_0272 + IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = m_d*s_14;
    const ccomplex_t IT_0277 = IT_0121*IT_0164;
    const ccomplex_t IT_0278 = IT_0055*IT_0098;
    const ccomplex_t IT_0279 = IT_0172*IT_0215;
    const ccomplex_t IT_0280 = IT_0028*IT_0037;
    const ccomplex_t IT_0281 = IT_0224 + -IT_0277 + -IT_0278 + -IT_0279 + 
      -IT_0280;
    const ccomplex_t IT_0282 = s_34*m_N_4;
    const ccomplex_t IT_0283 = 6*conj(IT_0219);
    const ccomplex_t IT_0284 = 6*conj(IT_0267);
    const ccomplex_t IT_0285 = pow(m_H0, 2);
    const ccomplex_t IT_0286 = m_d*m_N_4*IT_0285;
    const ccomplex_t IT_0287 = s_13*IT_0285;
    const ccomplex_t IT_0288 = s_14*s_34;
    const ccomplex_t IT_0289 = (-2)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = 6*IT_0267;
    const ccomplex_t IT_0292 = 0.166666666666667*IT_0282;
    const ccomplex_t IT_0293 = (IT_0219*conj(IT_0219) + IT_0267*conj(IT_0267))
      *IT_0268 + conj(IT_0267)*IT_0269*IT_0270 + (IT_0267*IT_0269 + IT_0275
      *IT_0276 + IT_0281*IT_0282)*IT_0283 + (IT_0276*IT_0281 + IT_0275*IT_0282)
      *IT_0284 + 6*conj(IT_0281)*(IT_0275*IT_0286 + -IT_0281*IT_0290 +
       0.166666666666667*IT_0276*IT_0291 + IT_0270*IT_0292) + 6*conj(IT_0275)*
      (0.166666666666667*IT_0270*IT_0276 + IT_0281*IT_0286 + -IT_0275*IT_0290 +
       IT_0291*IT_0292);
    return create_ccomplex_return(IT_0293);
}

