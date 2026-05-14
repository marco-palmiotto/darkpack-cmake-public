#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_smu_L_to_mu_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_smu_L_to_mu_h(
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = m_N_1*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0039 = IT_0007*IT_0038;
    const ccomplex_t IT_0040 = IT_0033 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = -IT_0042 + -IT_0044 + -IT_0046 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0040 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = m_N_3*IT_0052;
    const ccomplex_t IT_0054 = IT_0031 + IT_0053;
    const ccomplex_t IT_0055 = pow(m_W, -1);
    const ccomplex_t IT_0056 = cos(beta);
    const ccomplex_t IT_0057 = cpow(IT_0056, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = pow(m_mu, 2);
    const ccomplex_t IT_0061 = pow(m_smu_L, 2);
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0060 + IT_0061 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0059*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0054*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0000
      *IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0002;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = pow(m_N_3, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0060 + IT_0071 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = sin(beta);
    const ccomplex_t IT_0076 = cpow(IT_0001, -2);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*e_em*(IT_0006*IT_0016
      *IT_0055*IT_0057*IT_0060 + (-0.5)*m_W*(IT_0006*IT_0056 + IT_0000*IT_0075)*
      (IT_0016 + -IT_0015*IT_0076));
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_13 + IT_0060 + IT_0071 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0085 = IT_0003*IT_0084;
    const ccomplex_t IT_0086 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0087 = IT_0007*IT_0086;
    const ccomplex_t IT_0088 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0089 = IT_0017*IT_0088;
    const ccomplex_t IT_0090 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0091 = IT_0020*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + IT_0087 + 
      -IT_0089 + -IT_0091);
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = m_N_3*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0098 = IT_0007*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0100 = IT_0017*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0102 = IT_0020*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + IT_0098 + 
      -IT_0100 + -IT_0102);
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = m_N_3*IT_0104;
    const ccomplex_t IT_0106 = IT_0094 + IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0060 + IT_0061 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0079*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0106*IT_0109;
    const ccomplex_t IT_0111 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0112 = IT_0003*IT_0111;
    const ccomplex_t IT_0113 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0114 = IT_0007*IT_0113;
    const ccomplex_t IT_0115 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0116 = IT_0003*IT_0115;
    const ccomplex_t IT_0117 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0118 = IT_0007*IT_0117;
    const ccomplex_t IT_0119 = IT_0112 + IT_0114 + IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0123 = IT_0020*IT_0122;
    const ccomplex_t IT_0124 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0125 = IT_0017*IT_0124;
    const ccomplex_t IT_0126 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0127 = IT_0020*IT_0126;
    const ccomplex_t IT_0128 = -IT_0121 + -IT_0123 + -IT_0125 + -IT_0127;
    const ccomplex_t IT_0129 = IT_0119 + IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = m_N_2*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0134 = IT_0003*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0136 = IT_0007*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0138 = IT_0003*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0140 = IT_0007*IT_0139;
    const ccomplex_t IT_0141 = IT_0134 + IT_0136 + IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0143 = IT_0017*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0145 = IT_0020*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0147 = IT_0017*IT_0146;
    const ccomplex_t IT_0148 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0149 = IT_0020*IT_0148;
    const ccomplex_t IT_0150 = -IT_0143 + -IT_0145 + -IT_0147 + -IT_0149;
    const ccomplex_t IT_0151 = IT_0141 + IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = m_N_3*IT_0153;
    const ccomplex_t IT_0155 = IT_0132 + IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_23 + IT_0060 + IT_0061 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0155*IT_0160;
    const ccomplex_t IT_0162 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0163 = IT_0003*IT_0162;
    const ccomplex_t IT_0164 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0165 = IT_0007*IT_0164;
    const ccomplex_t IT_0166 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0167 = IT_0003*IT_0166;
    const ccomplex_t IT_0168 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0169 = IT_0007*IT_0168;
    const ccomplex_t IT_0170 = IT_0163 + IT_0165 + IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0172 = IT_0017*IT_0171;
    const ccomplex_t IT_0173 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0174 = IT_0020*IT_0173;
    const ccomplex_t IT_0175 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0176 = IT_0017*IT_0175;
    const ccomplex_t IT_0177 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0178 = IT_0020*IT_0177;
    const ccomplex_t IT_0179 = -IT_0172 + -IT_0174 + -IT_0176 + -IT_0178;
    const ccomplex_t IT_0180 = IT_0170 + IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = m_N_4*IT_0182;
    const ccomplex_t IT_0184 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0185 = IT_0003*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0187 = IT_0007*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0189 = IT_0003*IT_0188;
    const ccomplex_t IT_0190 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0191 = IT_0007*IT_0190;
    const ccomplex_t IT_0192 = IT_0185 + IT_0187 + IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0194 = IT_0017*IT_0193;
    const ccomplex_t IT_0195 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0196 = IT_0020*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0198 = IT_0017*IT_0197;
    const ccomplex_t IT_0199 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0200 = IT_0020*IT_0199;
    const ccomplex_t IT_0201 = -IT_0194 + -IT_0196 + -IT_0198 + -IT_0200;
    const ccomplex_t IT_0202 = IT_0192 + IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = m_N_3*IT_0204;
    const ccomplex_t IT_0206 = IT_0183 + IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = cpow((-2)*s_23 + IT_0060 + IT_0061 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = IT_0206*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*e_em*m_mu*IT_0006*IT_0016
      *IT_0055*IT_0057;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0079*IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_12 + IT_0060 + -IT_0061 + -IT_0071 
      + -reg_prop, -1);
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = m_mu*IT_0218;
    const ccomplex_t IT_0220 = IT_0065 + IT_0074 + IT_0083 + IT_0110 + IT_0161
       + IT_0212 + (-2)*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B3)*e_em;
    const ccomplex_t IT_0222 = IT_0002*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W3)*e_em;
    const ccomplex_t IT_0225 = IT_0016*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + IT_0226);
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0077*IT_0228;
    const ccomplex_t IT_0230 = IT_0081*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0016*IT_0055*IT_0057;
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = IT_0067*IT_0233;
    const ccomplex_t IT_0235 = IT_0072*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0107*IT_0228;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = IT_0106*IT_0238;
    const ccomplex_t IT_0240 = m_N_3*IT_0030;
    const ccomplex_t IT_0241 = m_N_1*IT_0052;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = conj(N_B1)*e_em;
    const ccomplex_t IT_0244 = IT_0002*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W1)*e_em;
    const ccomplex_t IT_0247 = IT_0016*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*(IT_0245 + IT_0248);
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = IT_0062*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = IT_0242*IT_0252;
    const ccomplex_t IT_0254 = m_N_3*IT_0131;
    const ccomplex_t IT_0255 = m_N_2*IT_0153;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = conj(N_B2)*e_em;
    const ccomplex_t IT_0258 = IT_0002*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W2)*e_em;
    const ccomplex_t IT_0261 = IT_0016*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*(IT_0259 + IT_0262);
    const ccomplex_t IT_0264 = (-0.5)*IT_0263;
    const ccomplex_t IT_0265 = IT_0158*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0256*IT_0266;
    const ccomplex_t IT_0268 = m_N_3*IT_0182;
    const ccomplex_t IT_0269 = m_N_4*IT_0204;
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = conj(N_B4)*e_em;
    const ccomplex_t IT_0272 = IT_0002*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = conj(N_W4)*e_em;
    const ccomplex_t IT_0275 = IT_0016*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*(IT_0273 + IT_0276);
    const ccomplex_t IT_0278 = (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = IT_0209*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = IT_0270*IT_0280;
    const ccomplex_t IT_0282 = IT_0214*IT_0228;
    const ccomplex_t IT_0283 = IT_0216*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = m_mu*IT_0284;
    const ccomplex_t IT_0286 = -IT_0231 + -IT_0236 + -IT_0239 + -IT_0253 + 
      -IT_0267 + -IT_0281 + 2*IT_0285;
    const ccomplex_t IT_0287 = 2*s_13;
    const ccomplex_t IT_0288 = m_mu*m_N_3;
    const ccomplex_t IT_0289 = 2*IT_0220;
    const ccomplex_t IT_0290 = IT_0093*IT_0238;
    const ccomplex_t IT_0291 = IT_0030*IT_0252;
    const ccomplex_t IT_0292 = IT_0131*IT_0266;
    const ccomplex_t IT_0293 = IT_0182*IT_0280;
    const ccomplex_t IT_0294 = -IT_0218 + IT_0290 + IT_0291 + IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = m_mu*s_14;
    const ccomplex_t IT_0296 = IT_0104*IT_0109;
    const ccomplex_t IT_0297 = IT_0052*IT_0064;
    const ccomplex_t IT_0298 = IT_0153*IT_0160;
    const ccomplex_t IT_0299 = IT_0204*IT_0211;
    const ccomplex_t IT_0300 = IT_0284 + -IT_0296 + -IT_0297 + -IT_0298 + 
      -IT_0299;
    const ccomplex_t IT_0301 = s_34*m_N_3;
    const ccomplex_t IT_0302 = 2*conj(IT_0220);
    const ccomplex_t IT_0303 = 2*conj(IT_0286);
    const ccomplex_t IT_0304 = pow(m_h, 2);
    const ccomplex_t IT_0305 = m_mu*m_N_3*IT_0304;
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

