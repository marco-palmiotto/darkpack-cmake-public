#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_L_to_mu_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_L_to_mu_h(
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
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0003
      *IT_0005*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = IT_0002*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = sin(alpha);
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = cpow(IT_0011, -2);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0003
      *IT_0005*IT_0008*IT_0018 + (-0.5)*m_W*(IT_0004*IT_0018 + IT_0006*IT_0019)*
      (IT_0008 + -IT_0007*IT_0020));
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = pow(m_smu_L, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0023*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0033 = IT_0006*IT_0012;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0036 = IT_0012*IT_0018;
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0039 = IT_0006*IT_0008;
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0042 = IT_0008*IT_0018;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0034 + IT_0037 + 
      -IT_0040 + -IT_0043);
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = m_N_1*IT_0045;
    const ccomplex_t IT_0047 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0048 = IT_0033*IT_0047;
    const ccomplex_t IT_0049 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0050 = IT_0036*IT_0049;
    const ccomplex_t IT_0051 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0052 = IT_0039*IT_0051;
    const ccomplex_t IT_0053 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0054 = IT_0042*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0048 + IT_0050 + 
      -IT_0052 + -IT_0054);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = m_N_1*IT_0056;
    const ccomplex_t IT_0058 = IT_0046 + IT_0057;
    const ccomplex_t IT_0059 = IT_0031*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0066 = IT_0033*IT_0065;
    const ccomplex_t IT_0067 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0068 = IT_0036*IT_0067;
    const ccomplex_t IT_0069 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0070 = IT_0033*IT_0069;
    const ccomplex_t IT_0071 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0072 = IT_0036*IT_0071;
    const ccomplex_t IT_0073 = IT_0066 + IT_0068 + IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0075 = IT_0039*IT_0074;
    const ccomplex_t IT_0076 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0077 = IT_0042*IT_0076;
    const ccomplex_t IT_0078 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0079 = IT_0039*IT_0078;
    const ccomplex_t IT_0080 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0081 = IT_0042*IT_0080;
    const ccomplex_t IT_0082 = -IT_0075 + -IT_0077 + -IT_0079 + -IT_0081;
    const ccomplex_t IT_0083 = IT_0073 + IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = m_N_4*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0088 = IT_0033*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0090 = IT_0036*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0092 = IT_0033*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0094 = IT_0036*IT_0093;
    const ccomplex_t IT_0095 = IT_0088 + IT_0090 + IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0097 = IT_0039*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0099 = IT_0042*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0101 = IT_0039*IT_0100;
    const ccomplex_t IT_0102 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0103 = IT_0042*IT_0102;
    const ccomplex_t IT_0104 = -IT_0097 + -IT_0099 + -IT_0101 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0095 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = m_N_1*IT_0107;
    const ccomplex_t IT_0109 = IT_0086 + IT_0108;
    const ccomplex_t IT_0110 = IT_0064*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0117 = IT_0033*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0119 = IT_0036*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0121 = IT_0033*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0123 = IT_0036*IT_0122;
    const ccomplex_t IT_0124 = IT_0117 + IT_0119 + IT_0121 + IT_0123;
    const ccomplex_t IT_0125 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0126 = IT_0039*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0128 = IT_0042*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0130 = IT_0039*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0132 = IT_0042*IT_0131;
    const ccomplex_t IT_0133 = -IT_0126 + -IT_0128 + -IT_0130 + -IT_0132;
    const ccomplex_t IT_0134 = IT_0124 + IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = m_N_1*IT_0136;
    const ccomplex_t IT_0138 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0139 = IT_0033*IT_0138;
    const ccomplex_t IT_0140 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0141 = IT_0036*IT_0140;
    const ccomplex_t IT_0142 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0143 = IT_0033*IT_0142;
    const ccomplex_t IT_0144 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0145 = IT_0036*IT_0144;
    const ccomplex_t IT_0146 = IT_0139 + IT_0141 + IT_0143 + IT_0145;
    const ccomplex_t IT_0147 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0148 = IT_0039*IT_0147;
    const ccomplex_t IT_0149 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0150 = IT_0042*IT_0149;
    const ccomplex_t IT_0151 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0152 = IT_0039*IT_0151;
    const ccomplex_t IT_0153 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0154 = IT_0042*IT_0153;
    const ccomplex_t IT_0155 = -IT_0148 + -IT_0150 + -IT_0152 + -IT_0154;
    const ccomplex_t IT_0156 = IT_0146 + IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = m_N_2*IT_0158;
    const ccomplex_t IT_0160 = IT_0137 + IT_0159;
    const ccomplex_t IT_0161 = IT_0115*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0168 = IT_0033*IT_0167;
    const ccomplex_t IT_0169 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0170 = IT_0036*IT_0169;
    const ccomplex_t IT_0171 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0172 = IT_0033*IT_0171;
    const ccomplex_t IT_0173 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0174 = IT_0036*IT_0173;
    const ccomplex_t IT_0175 = IT_0168 + IT_0170 + IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0177 = IT_0039*IT_0176;
    const ccomplex_t IT_0178 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0179 = IT_0042*IT_0178;
    const ccomplex_t IT_0180 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0181 = IT_0039*IT_0180;
    const ccomplex_t IT_0182 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0183 = IT_0042*IT_0182;
    const ccomplex_t IT_0184 = -IT_0177 + -IT_0179 + -IT_0181 + -IT_0183;
    const ccomplex_t IT_0185 = IT_0175 + IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = m_N_3*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0190 = IT_0033*IT_0189;
    const ccomplex_t IT_0191 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0192 = IT_0036*IT_0191;
    const ccomplex_t IT_0193 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0194 = IT_0033*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0196 = IT_0036*IT_0195;
    const ccomplex_t IT_0197 = IT_0190 + IT_0192 + IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0199 = IT_0039*IT_0198;
    const ccomplex_t IT_0200 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0201 = IT_0042*IT_0200;
    const ccomplex_t IT_0202 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0203 = IT_0039*IT_0202;
    const ccomplex_t IT_0204 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0205 = IT_0042*IT_0204;
    const ccomplex_t IT_0206 = -IT_0199 + -IT_0201 + -IT_0203 + -IT_0205;
    const ccomplex_t IT_0207 = IT_0197 + IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = m_N_1*IT_0209;
    const ccomplex_t IT_0211 = IT_0188 + IT_0210;
    const ccomplex_t IT_0212 = IT_0166*IT_0211;
    const ccomplex_t IT_0213 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0028 
      + -reg_prop, -1);
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*e_em*m_mu*IT_0003*IT_0005
      *IT_0008*IT_0018;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = IT_0023*IT_0215;
    const ccomplex_t IT_0217 = IT_0213*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = m_mu*IT_0218;
    const ccomplex_t IT_0220 = IT_0017 + IT_0027 + IT_0059 + IT_0110 + IT_0161
       + IT_0212 + (-2)*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B1)*e_em;
    const ccomplex_t IT_0222 = IT_0012*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W1)*e_em;
    const ccomplex_t IT_0225 = IT_0008*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + IT_0226);
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0021*IT_0228;
    const ccomplex_t IT_0230 = IT_0025*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = IT_0215*IT_0228;
    const ccomplex_t IT_0233 = IT_0213*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = m_mu*IT_0234;
    const ccomplex_t IT_0236 = m_N_4*IT_0107;
    const ccomplex_t IT_0237 = m_N_1*IT_0085;
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = conj(N_B4)*e_em;
    const ccomplex_t IT_0240 = IT_0012*IT_0239;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W4)*e_em;
    const ccomplex_t IT_0243 = IT_0008*IT_0242;
    const ccomplex_t IT_0244 = 1.4142135623731*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*(IT_0241 + IT_0244);
    const ccomplex_t IT_0246 = (-0.5)*IT_0245;
    const ccomplex_t IT_0247 = IT_0062*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0238*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = IT_0010*IT_0251;
    const ccomplex_t IT_0253 = IT_0002*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0029*IT_0228;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0058*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B2)*e_em;
    const ccomplex_t IT_0259 = IT_0012*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W2)*e_em;
    const ccomplex_t IT_0262 = IT_0008*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0260 + IT_0263);
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = IT_0113*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = m_N_1*IT_0158;
    const ccomplex_t IT_0269 = m_N_2*IT_0136;
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = IT_0267*IT_0270;
    const ccomplex_t IT_0272 = conj(N_B3)*e_em;
    const ccomplex_t IT_0273 = IT_0012*IT_0272;
    const ccomplex_t IT_0274 = 1.4142135623731*IT_0273;
    const ccomplex_t IT_0275 = conj(N_W3)*e_em;
    const ccomplex_t IT_0276 = IT_0008*IT_0275;
    const ccomplex_t IT_0277 = 1.4142135623731*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*(IT_0274 + IT_0277);
    const ccomplex_t IT_0279 = (-0.5)*IT_0278;
    const ccomplex_t IT_0280 = IT_0164*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = m_N_3*IT_0209;
    const ccomplex_t IT_0283 = m_N_1*IT_0187;
    const ccomplex_t IT_0284 = IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = IT_0281*IT_0284;
    const ccomplex_t IT_0286 = -IT_0231 + 2*IT_0235 + -IT_0249 + -IT_0254 + 
      -IT_0257 + -IT_0271 + -IT_0285;
    const ccomplex_t IT_0287 = 2*s_13;
    const ccomplex_t IT_0288 = m_mu*m_N_1;
    const ccomplex_t IT_0289 = 2*IT_0220;
    const ccomplex_t IT_0290 = IT_0187*IT_0281;
    const ccomplex_t IT_0291 = IT_0056*IT_0256;
    const ccomplex_t IT_0292 = IT_0158*IT_0267;
    const ccomplex_t IT_0293 = IT_0085*IT_0248;
    const ccomplex_t IT_0294 = -IT_0218 + IT_0290 + IT_0291 + IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = m_mu*s_14;
    const ccomplex_t IT_0296 = IT_0031*IT_0045;
    const ccomplex_t IT_0297 = IT_0115*IT_0136;
    const ccomplex_t IT_0298 = IT_0166*IT_0209;
    const ccomplex_t IT_0299 = IT_0064*IT_0107;
    const ccomplex_t IT_0300 = IT_0234 + -IT_0296 + -IT_0297 + -IT_0298 + 
      -IT_0299;
    const ccomplex_t IT_0301 = s_34*m_N_1;
    const ccomplex_t IT_0302 = 2*conj(IT_0220);
    const ccomplex_t IT_0303 = 2*conj(IT_0286);
    const ccomplex_t IT_0304 = pow(m_h, 2);
    const ccomplex_t IT_0305 = m_mu*m_N_1*IT_0304;
    const ccomplex_t IT_0306 = s_13*IT_0304;
    const ccomplex_t IT_0307 = s_14*s_34;
    const ccomplex_t IT_0308 = (-2)*IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0286;
    const ccomplex_t IT_0311 = 0.5*IT_0301;
    const ccomplex_t IT_0312 = (IT_0220*conj(IT_0220) + IT_0286*conj(IT_0286))
      *IT_0287 + conj(IT_0286)*IT_0288*IT_0289 + (IT_0286*IT_0288 + IT_0294
      *IT_0295 + IT_0300*IT_0301)*IT_0302 + (IT_0295*IT_0300 + IT_0294*IT_0301)
      *IT_0303 + 2*conj(IT_0300)*(IT_0294*IT_0305 + -IT_0300*IT_0309 + 0.5
      *IT_0295*IT_0310 + IT_0289*IT_0311) + 2*conj(IT_0294)*(0.5*IT_0289*IT_0295
       + IT_0300*IT_0305 + -IT_0294*IT_0309 + IT_0310*IT_0311);
    return create_ccomplex_return(IT_0312);
}

