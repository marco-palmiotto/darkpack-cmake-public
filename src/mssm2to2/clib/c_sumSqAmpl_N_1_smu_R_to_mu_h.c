#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_R_to_mu_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_R_to_mu_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
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
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
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
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0003
      *IT_0005*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -2);
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = sin(alpha);
    const ccomplex_t IT_0021 = (0 + _Complex_I*-1)*e_em*(IT_0000*IT_0003
      *IT_0005*IT_0008*IT_0020 + -m_W*IT_0007*IT_0018*(IT_0006*IT_0019 + IT_0004
      *IT_0020));
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0017, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = IT_0016*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_smu_R, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0034 = IT_0006*IT_0023;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0037 = IT_0020*IT_0023;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0040 = IT_0006*IT_0008;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0043 = IT_0008*IT_0020;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0035 + IT_0038 + 
      -IT_0041 + -IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = m_N_1*IT_0046;
    const ccomplex_t IT_0048 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0049 = IT_0034*IT_0048;
    const ccomplex_t IT_0050 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0051 = IT_0037*IT_0050;
    const ccomplex_t IT_0052 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0053 = IT_0040*IT_0052;
    const ccomplex_t IT_0054 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0055 = IT_0043*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0049 + IT_0051 + 
      -IT_0053 + -IT_0055);
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = m_N_1*IT_0057;
    const ccomplex_t IT_0059 = IT_0047 + IT_0058;
    const ccomplex_t IT_0060 = IT_0032*IT_0059;
    const ccomplex_t IT_0061 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0062 = IT_0034*IT_0061;
    const ccomplex_t IT_0063 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0064 = IT_0037*IT_0063;
    const ccomplex_t IT_0065 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0066 = IT_0034*IT_0065;
    const ccomplex_t IT_0067 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0068 = IT_0037*IT_0067;
    const ccomplex_t IT_0069 = IT_0062 + IT_0064 + IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0071 = IT_0040*IT_0070;
    const ccomplex_t IT_0072 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0073 = IT_0043*IT_0072;
    const ccomplex_t IT_0074 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0075 = IT_0040*IT_0074;
    const ccomplex_t IT_0076 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0077 = IT_0043*IT_0076;
    const ccomplex_t IT_0078 = -IT_0071 + -IT_0073 + -IT_0075 + -IT_0077;
    const ccomplex_t IT_0079 = IT_0069 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = m_N_2*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0034*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0037*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0088 = IT_0034*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0090 = IT_0037*IT_0089;
    const ccomplex_t IT_0091 = IT_0084 + IT_0086 + IT_0088 + IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0093 = IT_0040*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0095 = IT_0043*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0097 = IT_0040*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0099 = IT_0043*IT_0098;
    const ccomplex_t IT_0100 = -IT_0093 + -IT_0095 + -IT_0097 + -IT_0099;
    const ccomplex_t IT_0101 = IT_0091 + IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = m_N_1*IT_0103;
    const ccomplex_t IT_0105 = IT_0082 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0023;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0105*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0023;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0118 = IT_0034*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0120 = IT_0037*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0122 = IT_0034*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0124 = IT_0037*IT_0123;
    const ccomplex_t IT_0125 = IT_0118 + IT_0120 + IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0127 = IT_0040*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0129 = IT_0043*IT_0128;
    const ccomplex_t IT_0130 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0131 = IT_0040*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0133 = IT_0043*IT_0132;
    const ccomplex_t IT_0134 = -IT_0127 + -IT_0129 + -IT_0131 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0125 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = m_N_1*IT_0137;
    const ccomplex_t IT_0139 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0140 = IT_0034*IT_0139;
    const ccomplex_t IT_0141 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0142 = IT_0037*IT_0141;
    const ccomplex_t IT_0143 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0144 = IT_0034*IT_0143;
    const ccomplex_t IT_0145 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0146 = IT_0037*IT_0145;
    const ccomplex_t IT_0147 = IT_0140 + IT_0142 + IT_0144 + IT_0146;
    const ccomplex_t IT_0148 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0149 = IT_0040*IT_0148;
    const ccomplex_t IT_0150 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0151 = IT_0043*IT_0150;
    const ccomplex_t IT_0152 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0153 = IT_0040*IT_0152;
    const ccomplex_t IT_0154 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0155 = IT_0043*IT_0154;
    const ccomplex_t IT_0156 = -IT_0149 + -IT_0151 + -IT_0153 + -IT_0155;
    const ccomplex_t IT_0157 = IT_0147 + IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = m_N_3*IT_0159;
    const ccomplex_t IT_0161 = IT_0138 + IT_0160;
    const ccomplex_t IT_0162 = IT_0116*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0023;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_23 + IT_0000 + IT_0029 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0169 = IT_0034*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0171 = IT_0037*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0173 = IT_0034*IT_0172;
    const ccomplex_t IT_0174 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0175 = IT_0037*IT_0174;
    const ccomplex_t IT_0176 = IT_0169 + IT_0171 + IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0178 = IT_0040*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0180 = IT_0043*IT_0179;
    const ccomplex_t IT_0181 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0182 = IT_0040*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0184 = IT_0043*IT_0183;
    const ccomplex_t IT_0185 = -IT_0178 + -IT_0180 + -IT_0182 + -IT_0184;
    const ccomplex_t IT_0186 = IT_0176 + IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = m_N_1*IT_0188;
    const ccomplex_t IT_0190 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0191 = IT_0034*IT_0190;
    const ccomplex_t IT_0192 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0193 = IT_0037*IT_0192;
    const ccomplex_t IT_0194 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0195 = IT_0034*IT_0194;
    const ccomplex_t IT_0196 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0197 = IT_0037*IT_0196;
    const ccomplex_t IT_0198 = IT_0191 + IT_0193 + IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0200 = IT_0040*IT_0199;
    const ccomplex_t IT_0201 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0202 = IT_0043*IT_0201;
    const ccomplex_t IT_0203 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0204 = IT_0040*IT_0203;
    const ccomplex_t IT_0205 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0206 = IT_0043*IT_0205;
    const ccomplex_t IT_0207 = -IT_0200 + -IT_0202 + -IT_0204 + -IT_0206;
    const ccomplex_t IT_0208 = IT_0198 + IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = m_N_4*IT_0210;
    const ccomplex_t IT_0212 = IT_0189 + IT_0211;
    const ccomplex_t IT_0213 = IT_0167*IT_0212;
    const ccomplex_t IT_0214 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0029 
      + -reg_prop, -1);
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*e_em*m_mu*IT_0003*IT_0005
      *IT_0008*IT_0020;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = IT_0025*IT_0216;
    const ccomplex_t IT_0218 = IT_0214*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_mu*IT_0219;
    const ccomplex_t IT_0221 = IT_0015 + IT_0028 + IT_0060 + IT_0111 + IT_0162
       + IT_0213 + (-2)*IT_0220;
    const ccomplex_t IT_0222 = conj(N_B1)*e_em;
    const ccomplex_t IT_0223 = IT_0023*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = conj(N_W1)*e_em;
    const ccomplex_t IT_0226 = IT_0008*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0224 + IT_0227);
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = IT_0010*IT_0229;
    const ccomplex_t IT_0231 = IT_0002*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0234 = 0.5*IT_0233;
    const ccomplex_t IT_0235 = IT_0022*IT_0234;
    const ccomplex_t IT_0236 = IT_0016*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = m_N_1*IT_0159;
    const ccomplex_t IT_0239 = m_N_3*IT_0137;
    const ccomplex_t IT_0240 = IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0242 = 0.5*IT_0241;
    const ccomplex_t IT_0243 = IT_0114*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = IT_0240*IT_0244;
    const ccomplex_t IT_0246 = IT_0216*IT_0234;
    const ccomplex_t IT_0247 = IT_0214*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = m_mu*IT_0248;
    const ccomplex_t IT_0250 = IT_0030*IT_0234;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0059*IT_0251;
    const ccomplex_t IT_0253 = m_N_1*IT_0081;
    const ccomplex_t IT_0254 = m_N_2*IT_0103;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0257 = 0.5*IT_0256;
    const ccomplex_t IT_0258 = IT_0108*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = IT_0255*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0262 = 0.5*IT_0261;
    const ccomplex_t IT_0263 = IT_0165*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = m_N_1*IT_0210;
    const ccomplex_t IT_0266 = m_N_4*IT_0188;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = IT_0264*IT_0267;
    const ccomplex_t IT_0269 = -IT_0232 + -IT_0237 + -IT_0245 + 2*IT_0249 + 
      -IT_0252 + -IT_0260 + -IT_0268;
    const ccomplex_t IT_0270 = 2*s_13;
    const ccomplex_t IT_0271 = m_mu*m_N_1;
    const ccomplex_t IT_0272 = 2*IT_0221;
    const ccomplex_t IT_0273 = IT_0081*IT_0259;
    const ccomplex_t IT_0274 = IT_0057*IT_0251;
    const ccomplex_t IT_0275 = IT_0210*IT_0264;
    const ccomplex_t IT_0276 = IT_0159*IT_0244;
    const ccomplex_t IT_0277 = -IT_0219 + IT_0273 + IT_0274 + IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = m_mu*s_14;
    const ccomplex_t IT_0279 = IT_0032*IT_0046;
    const ccomplex_t IT_0280 = IT_0103*IT_0110;
    const ccomplex_t IT_0281 = IT_0116*IT_0137;
    const ccomplex_t IT_0282 = IT_0167*IT_0188;
    const ccomplex_t IT_0283 = IT_0248 + -IT_0279 + -IT_0280 + -IT_0281 + 
      -IT_0282;
    const ccomplex_t IT_0284 = s_34*m_N_1;
    const ccomplex_t IT_0285 = 2*conj(IT_0221);
    const ccomplex_t IT_0286 = 2*conj(IT_0269);
    const ccomplex_t IT_0287 = pow(m_h, 2);
    const ccomplex_t IT_0288 = m_mu*m_N_1*IT_0287;
    const ccomplex_t IT_0289 = s_13*IT_0287;
    const ccomplex_t IT_0290 = s_14*s_34;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = IT_0289 + IT_0291;
    const ccomplex_t IT_0293 = 2*IT_0269;
    const ccomplex_t IT_0294 = 0.5*IT_0284;
    const ccomplex_t IT_0295 = (IT_0221*conj(IT_0221) + IT_0269*conj(IT_0269))
      *IT_0270 + conj(IT_0269)*IT_0271*IT_0272 + (IT_0269*IT_0271 + IT_0277
      *IT_0278 + IT_0283*IT_0284)*IT_0285 + (IT_0278*IT_0283 + IT_0277*IT_0284)
      *IT_0286 + 2*conj(IT_0283)*(IT_0277*IT_0288 + -IT_0283*IT_0292 + 0.5
      *IT_0278*IT_0293 + IT_0272*IT_0294) + 2*conj(IT_0277)*(0.5*IT_0272*IT_0278
       + IT_0283*IT_0288 + -IT_0277*IT_0292 + IT_0293*IT_0294);
    return create_ccomplex_return(IT_0295);
}

