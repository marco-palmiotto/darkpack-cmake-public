#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_smu_L_to_mu_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_smu_L_to_mu_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0000
      *IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = pow(m_mu, 2);
    const ccomplex_t IT_0014 = pow(m_N_2, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_smu_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0013 + IT_0020 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cos(alpha);
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = IT_0003*IT_0005;
    const ccomplex_t IT_0029 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0032 = IT_0025*IT_0031;
    const ccomplex_t IT_0033 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = IT_0027 + IT_0030 + IT_0032 + IT_0034;
    const ccomplex_t IT_0036 = IT_0003*IT_0009;
    const ccomplex_t IT_0037 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0024;
    const ccomplex_t IT_0040 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    const ccomplex_t IT_0044 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = -IT_0038 + -IT_0041 + -IT_0043 + -IT_0045;
    const ccomplex_t IT_0047 = IT_0035 + IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = m_N_3*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0052 = IT_0025*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0054 = IT_0028*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0056 = IT_0025*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0058 = IT_0028*IT_0057;
    const ccomplex_t IT_0059 = IT_0052 + IT_0054 + IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0061 = IT_0036*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0063 = IT_0039*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0065 = IT_0036*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0067 = IT_0039*IT_0066;
    const ccomplex_t IT_0068 = -IT_0061 + -IT_0063 + -IT_0065 + -IT_0067;
    const ccomplex_t IT_0069 = IT_0059 + IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = m_N_2*IT_0071;
    const ccomplex_t IT_0073 = IT_0050 + IT_0072;
    const ccomplex_t IT_0074 = IT_0023*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0013 + IT_0020 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0081 = IT_0025*IT_0080;
    const ccomplex_t IT_0082 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0083 = IT_0028*IT_0082;
    const ccomplex_t IT_0084 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0085 = IT_0025*IT_0084;
    const ccomplex_t IT_0086 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0087 = IT_0028*IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0083 + IT_0085 + IT_0087;
    const ccomplex_t IT_0089 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0090 = IT_0036*IT_0089;
    const ccomplex_t IT_0091 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0092 = IT_0039*IT_0091;
    const ccomplex_t IT_0093 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0094 = IT_0036*IT_0093;
    const ccomplex_t IT_0095 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0096 = IT_0039*IT_0095;
    const ccomplex_t IT_0097 = -IT_0090 + -IT_0092 + -IT_0094 + -IT_0096;
    const ccomplex_t IT_0098 = IT_0088 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = m_N_4*IT_0100;
    const ccomplex_t IT_0102 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0103 = IT_0025*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0105 = IT_0028*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0107 = IT_0025*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0109 = IT_0028*IT_0108;
    const ccomplex_t IT_0110 = IT_0103 + IT_0105 + IT_0107 + IT_0109;
    const ccomplex_t IT_0111 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0112 = IT_0036*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0114 = IT_0039*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0116 = IT_0036*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0118 = IT_0039*IT_0117;
    const ccomplex_t IT_0119 = -IT_0112 + -IT_0114 + -IT_0116 + -IT_0118;
    const ccomplex_t IT_0120 = IT_0110 + IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = m_N_2*IT_0122;
    const ccomplex_t IT_0124 = IT_0101 + IT_0123;
    const ccomplex_t IT_0125 = IT_0079*IT_0124;
    const ccomplex_t IT_0126 = sin(beta);
    const ccomplex_t IT_0127 = cpow(IT_0008, -2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0013*IT_0024 + (-0.5)*m_W*(IT_0001*IT_0024 + -IT_0003*IT_0126)
      *(IT_0005 + -IT_0004*IT_0127));
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0132*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0137 = IT_0025*IT_0136;
    const ccomplex_t IT_0138 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0139 = IT_0036*IT_0138;
    const ccomplex_t IT_0140 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0141 = IT_0039*IT_0140;
    const ccomplex_t IT_0142 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0143 = IT_0028*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0137 + -IT_0139 + 
      -IT_0141 + IT_0143);
    const ccomplex_t IT_0145 = m_N_2*IT_0144;
    const ccomplex_t IT_0146 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0147 = IT_0025*IT_0146;
    const ccomplex_t IT_0148 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0149 = IT_0036*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0151 = IT_0039*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0153 = IT_0028*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + -IT_0149 + 
      -IT_0151 + IT_0153);
    const ccomplex_t IT_0155 = m_N_2*IT_0154;
    const ccomplex_t IT_0156 = IT_0145 + IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0013 + IT_0020 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0158 = IT_0131*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0156*IT_0159;
    const ccomplex_t IT_0161 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0162 = IT_0025*IT_0161;
    const ccomplex_t IT_0163 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0164 = IT_0028*IT_0163;
    const ccomplex_t IT_0165 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0166 = IT_0025*IT_0165;
    const ccomplex_t IT_0167 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0168 = IT_0028*IT_0167;
    const ccomplex_t IT_0169 = IT_0162 + IT_0164 + IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0171 = IT_0036*IT_0170;
    const ccomplex_t IT_0172 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0173 = IT_0039*IT_0172;
    const ccomplex_t IT_0174 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0175 = IT_0036*IT_0174;
    const ccomplex_t IT_0176 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0177 = IT_0039*IT_0176;
    const ccomplex_t IT_0178 = -IT_0171 + -IT_0173 + -IT_0175 + -IT_0177;
    const ccomplex_t IT_0179 = IT_0169 + IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0180;
    const ccomplex_t IT_0182 = m_N_1*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0184 = IT_0025*IT_0183;
    const ccomplex_t IT_0185 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0186 = IT_0028*IT_0185;
    const ccomplex_t IT_0187 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0188 = IT_0025*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0190 = IT_0028*IT_0189;
    const ccomplex_t IT_0191 = IT_0184 + IT_0186 + IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0193 = IT_0036*IT_0192;
    const ccomplex_t IT_0194 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0195 = IT_0039*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0197 = IT_0036*IT_0196;
    const ccomplex_t IT_0198 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0199 = IT_0039*IT_0198;
    const ccomplex_t IT_0200 = -IT_0193 + -IT_0195 + -IT_0197 + -IT_0199;
    const ccomplex_t IT_0201 = IT_0191 + IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = m_N_2*IT_0203;
    const ccomplex_t IT_0205 = IT_0182 + IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0207 = (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = cpow((-2)*s_23 + IT_0013 + IT_0020 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0209 = IT_0207*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = IT_0205*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0005*IT_0024;
    const ccomplex_t IT_0213 = (-0.5)*IT_0212;
    const ccomplex_t IT_0214 = IT_0131*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_12 + IT_0013 + -IT_0014 + -IT_0020 
      + -reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = m_mu*IT_0217;
    const ccomplex_t IT_0219 = IT_0017 + IT_0074 + IT_0125 + IT_0135 + IT_0160
       + IT_0211 + (-2)*IT_0218;
    const ccomplex_t IT_0220 = conj(N_B1)*e_em;
    const ccomplex_t IT_0221 = IT_0009*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W1)*e_em;
    const ccomplex_t IT_0224 = IT_0005*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + IT_0225);
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = IT_0208*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = m_N_2*IT_0181;
    const ccomplex_t IT_0231 = m_N_1*IT_0203;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0229*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B3)*e_em;
    const ccomplex_t IT_0235 = IT_0009*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W3)*e_em;
    const ccomplex_t IT_0238 = IT_0005*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + IT_0239);
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = IT_0021*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = m_N_2*IT_0049;
    const ccomplex_t IT_0245 = m_N_3*IT_0071;
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = IT_0243*IT_0246;
    const ccomplex_t IT_0248 = conj(N_B2)*e_em;
    const ccomplex_t IT_0249 = IT_0009*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = conj(N_W2)*e_em;
    const ccomplex_t IT_0252 = IT_0005*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*(IT_0250 + IT_0253);
    const ccomplex_t IT_0255 = (-0.5)*IT_0254;
    const ccomplex_t IT_0256 = IT_0213*IT_0255;
    const ccomplex_t IT_0257 = IT_0215*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = m_mu*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0261 = 0.5*IT_0260;
    const ccomplex_t IT_0262 = IT_0007*IT_0261;
    const ccomplex_t IT_0263 = IT_0015*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = conj(N_B4)*e_em;
    const ccomplex_t IT_0266 = IT_0009*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = conj(N_W4)*e_em;
    const ccomplex_t IT_0269 = IT_0005*IT_0268;
    const ccomplex_t IT_0270 = 1.4142135623731*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*(IT_0267 + IT_0270);
    const ccomplex_t IT_0272 = (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = IT_0077*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = m_N_2*IT_0100;
    const ccomplex_t IT_0276 = m_N_4*IT_0122;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = IT_0274*IT_0277;
    const ccomplex_t IT_0279 = IT_0129*IT_0255;
    const ccomplex_t IT_0280 = IT_0133*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0157*IT_0255;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = IT_0156*IT_0283;
    const ccomplex_t IT_0285 = -IT_0233 + -IT_0247 + 2*IT_0259 + -IT_0264 + 
      -IT_0278 + -IT_0281 + -IT_0284;
    const ccomplex_t IT_0286 = 2*s_13;
    const ccomplex_t IT_0287 = m_mu*m_N_2;
    const ccomplex_t IT_0288 = 2*IT_0219;
    const ccomplex_t IT_0289 = IT_0049*IT_0243;
    const ccomplex_t IT_0290 = IT_0100*IT_0274;
    const ccomplex_t IT_0291 = IT_0144*IT_0283;
    const ccomplex_t IT_0292 = IT_0181*IT_0229;
    const ccomplex_t IT_0293 = -IT_0217 + IT_0289 + IT_0290 + IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = m_mu*s_14;
    const ccomplex_t IT_0295 = IT_0023*IT_0071;
    const ccomplex_t IT_0296 = IT_0079*IT_0122;
    const ccomplex_t IT_0297 = IT_0154*IT_0159;
    const ccomplex_t IT_0298 = IT_0203*IT_0210;
    const ccomplex_t IT_0299 = IT_0258 + -IT_0295 + -IT_0296 + -IT_0297 + 
      -IT_0298;
    const ccomplex_t IT_0300 = s_34*m_N_2;
    const ccomplex_t IT_0301 = 2*conj(IT_0219);
    const ccomplex_t IT_0302 = 2*conj(IT_0285);
    const ccomplex_t IT_0303 = pow(m_H0, 2);
    const ccomplex_t IT_0304 = m_mu*m_N_2*IT_0303;
    const ccomplex_t IT_0305 = s_13*IT_0303;
    const ccomplex_t IT_0306 = s_14*s_34;
    const ccomplex_t IT_0307 = (-2)*IT_0306;
    const ccomplex_t IT_0308 = IT_0305 + IT_0307;
    const ccomplex_t IT_0309 = 2*IT_0285;
    const ccomplex_t IT_0310 = 0.5*IT_0300;
    const ccomplex_t IT_0311 = (IT_0219*conj(IT_0219) + IT_0285*conj(IT_0285))
      *IT_0286 + conj(IT_0285)*IT_0287*IT_0288 + (IT_0285*IT_0287 + IT_0293
      *IT_0294 + IT_0299*IT_0300)*IT_0301 + (IT_0294*IT_0299 + IT_0293*IT_0300)
      *IT_0302 + 2*conj(IT_0299)*(IT_0293*IT_0304 + -IT_0299*IT_0308 + 0.5
      *IT_0294*IT_0309 + IT_0288*IT_0310) + 2*conj(IT_0293)*(0.5*IT_0288*IT_0294
       + IT_0299*IT_0304 + -IT_0293*IT_0308 + IT_0309*IT_0310);
    return create_ccomplex_return(IT_0311);
}

