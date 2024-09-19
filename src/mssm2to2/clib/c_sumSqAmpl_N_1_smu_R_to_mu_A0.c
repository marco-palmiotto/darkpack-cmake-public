#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_R_to_mu_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_R_to_mu_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_mu = param->m_mu;
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
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_mu, 2);
    const ccomplex_t IT_0005 = pow(m_smu_R, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0004 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = IT_0001*IT_0013;
    const ccomplex_t IT_0015 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0009*IT_0018;
    const ccomplex_t IT_0020 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0026 = IT_0010*IT_0025;
    const ccomplex_t IT_0027 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0032 = IT_0022*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0012 + 0.5*IT_0016 + 0.5*IT_0021 + (
      -0.5)*IT_0024 + (-0.5)*IT_0026 + 0.5*IT_0028 + 0.5*IT_0030 + (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = m_N_3*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0038 = IT_0014*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0040 = IT_0019*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0042 = IT_0022*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0044 = IT_0010*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0046 = IT_0014*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0048 = IT_0019*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0036 + (-0.5)*IT_0038 + (-0.5)*IT_0040 +
       0.5*IT_0042 + 0.5*IT_0044 + (-0.5)*IT_0046 + (-0.5)*IT_0048 + 0.5*IT_0050;
    const ccomplex_t IT_0052 = m_N_1*IT_0051;
    const ccomplex_t IT_0053 = IT_0034 + IT_0052;
    const ccomplex_t IT_0054 = IT_0008*IT_0053;
    const ccomplex_t IT_0055 = pow(m_N_1, 2);
    const ccomplex_t IT_0056 = cpow((-2)*s_13 + IT_0004 + IT_0055 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0057 = pow(m_W, -1);
    const ccomplex_t IT_0058 = e_em*m_mu*mu_h*IT_0018*IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0009, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0059*IT_0062;
    const ccomplex_t IT_0064 = IT_0056*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0067 = IT_0010*IT_0066;
    const ccomplex_t IT_0068 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0069 = IT_0014*IT_0068;
    const ccomplex_t IT_0070 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0071 = IT_0019*IT_0070;
    const ccomplex_t IT_0072 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0073 = IT_0022*IT_0072;
    const ccomplex_t IT_0074 = -IT_0067 + IT_0069 + IT_0071 + -IT_0073;
    const ccomplex_t IT_0075 = m_N_1*IT_0074;
    const ccomplex_t IT_0076 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0079 = IT_0014*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0081 = IT_0019*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0083 = IT_0022*IT_0082;
    const ccomplex_t IT_0084 = IT_0077 + -IT_0079 + -IT_0081 + IT_0083;
    const ccomplex_t IT_0085 = m_N_1*IT_0084;
    const ccomplex_t IT_0086 = IT_0075 + IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0004 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0086*IT_0091;
    const ccomplex_t IT_0093 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0094 = IT_0010*IT_0093;
    const ccomplex_t IT_0095 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0096 = IT_0014*IT_0095;
    const ccomplex_t IT_0097 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0098 = IT_0019*IT_0097;
    const ccomplex_t IT_0099 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0100 = IT_0022*IT_0099;
    const ccomplex_t IT_0101 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0102 = IT_0010*IT_0101;
    const ccomplex_t IT_0103 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0104 = IT_0014*IT_0103;
    const ccomplex_t IT_0105 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0106 = IT_0019*IT_0105;
    const ccomplex_t IT_0107 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0108 = IT_0022*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0094 + 0.5*IT_0096 + 0.5*IT_0098 + (
      -0.5)*IT_0100 + (-0.5)*IT_0102 + 0.5*IT_0104 + 0.5*IT_0106 + (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = m_N_2*IT_0109;
    const ccomplex_t IT_0111 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0112 = IT_0010*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0114 = IT_0014*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0116 = IT_0019*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0118 = IT_0022*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0120 = IT_0010*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0122 = IT_0014*IT_0121;
    const ccomplex_t IT_0123 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0124 = IT_0019*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0126 = IT_0022*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0112 + (-0.5)*IT_0114 + (-0.5)*IT_0116 +
       0.5*IT_0118 + 0.5*IT_0120 + (-0.5)*IT_0122 + (-0.5)*IT_0124 + 0.5*IT_0126;
    const ccomplex_t IT_0128 = m_N_1*IT_0127;
    const ccomplex_t IT_0129 = IT_0110 + IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = cpow((-2)*s_23 + IT_0004 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0129*IT_0134;
    const ccomplex_t IT_0136 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0137 = IT_0010*IT_0136;
    const ccomplex_t IT_0138 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0139 = IT_0014*IT_0138;
    const ccomplex_t IT_0140 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0141 = IT_0019*IT_0140;
    const ccomplex_t IT_0142 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0143 = IT_0022*IT_0142;
    const ccomplex_t IT_0144 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0145 = IT_0010*IT_0144;
    const ccomplex_t IT_0146 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0147 = IT_0014*IT_0146;
    const ccomplex_t IT_0148 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0149 = IT_0019*IT_0148;
    const ccomplex_t IT_0150 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0151 = IT_0022*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0137 + 0.5*IT_0139 + 0.5*IT_0141 + (
      -0.5)*IT_0143 + (-0.5)*IT_0145 + 0.5*IT_0147 + 0.5*IT_0149 + (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = m_N_4*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0155 = IT_0010*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0157 = IT_0014*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0159 = IT_0019*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0161 = IT_0022*IT_0160;
    const ccomplex_t IT_0162 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0163 = IT_0010*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0165 = IT_0014*IT_0164;
    const ccomplex_t IT_0166 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0167 = IT_0019*IT_0166;
    const ccomplex_t IT_0168 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0169 = IT_0022*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0155 + (-0.5)*IT_0157 + (-0.5)*IT_0159 +
       0.5*IT_0161 + 0.5*IT_0163 + (-0.5)*IT_0165 + (-0.5)*IT_0167 + 0.5*IT_0169;
    const ccomplex_t IT_0171 = m_N_1*IT_0170;
    const ccomplex_t IT_0172 = IT_0153 + IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_23 + IT_0004 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0176 = IT_0174*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0172*IT_0177;
    const ccomplex_t IT_0179 = IT_0054 + IT_0065 + IT_0092 + IT_0135 + IT_0178;
    const ccomplex_t IT_0180 = conj(N_B1)*e_em;
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W1)*e_em;
    const ccomplex_t IT_0184 = IT_0018*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + IT_0185);
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = IT_0059*IT_0187;
    const ccomplex_t IT_0189 = IT_0056*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0089*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0086*IT_0194;
    const ccomplex_t IT_0196 = m_N_1*IT_0109;
    const ccomplex_t IT_0197 = m_N_2*IT_0127;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0200 = 0.5*IT_0199;
    const ccomplex_t IT_0201 = IT_0132*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0198*IT_0202;
    const ccomplex_t IT_0204 = m_N_1*IT_0033;
    const ccomplex_t IT_0205 = m_N_3*IT_0051;
    const ccomplex_t IT_0206 = IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0208 = 0.5*IT_0207;
    const ccomplex_t IT_0209 = IT_0006*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = IT_0206*IT_0210;
    const ccomplex_t IT_0212 = m_N_1*IT_0152;
    const ccomplex_t IT_0213 = m_N_4*IT_0170;
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = IT_0175*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0214*IT_0218;
    const ccomplex_t IT_0220 = -IT_0190 + -IT_0195 + -IT_0203 + -IT_0211 + 
      -IT_0219;
    const ccomplex_t IT_0221 = 2*s_13;
    const ccomplex_t IT_0222 = m_mu*m_N_1;
    const ccomplex_t IT_0223 = 2*IT_0179;
    const ccomplex_t IT_0224 = e_em*m_mu*IT_0013*IT_0018*IT_0057*IT_0060;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = IT_0088*IT_0225;
    const ccomplex_t IT_0227 = cpow((-2)*s_12 + IT_0004 + -IT_0005 + -IT_0055 
      + -reg_prop, -1);
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = IT_0109*IT_0202;
    const ccomplex_t IT_0231 = IT_0033*IT_0210;
    const ccomplex_t IT_0232 = IT_0074*IT_0194;
    const ccomplex_t IT_0233 = IT_0152*IT_0218;
    const ccomplex_t IT_0234 = -IT_0229 + IT_0230 + IT_0231 + IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = m_mu*s_14;
    const ccomplex_t IT_0236 = IT_0192*IT_0225;
    const ccomplex_t IT_0237 = IT_0227*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = IT_0084*IT_0091;
    const ccomplex_t IT_0240 = IT_0127*IT_0134;
    const ccomplex_t IT_0241 = IT_0008*IT_0051;
    const ccomplex_t IT_0242 = IT_0170*IT_0177;
    const ccomplex_t IT_0243 = -IT_0238 + -IT_0239 + -IT_0240 + -IT_0241 + 
      -IT_0242;
    const ccomplex_t IT_0244 = s_34*m_N_1;
    const ccomplex_t IT_0245 = 2*conj(IT_0179);
    const ccomplex_t IT_0246 = 2*conj(IT_0220);
    const ccomplex_t IT_0247 = pow(m_A0, 2);
    const ccomplex_t IT_0248 = m_mu*m_N_1*IT_0247;
    const ccomplex_t IT_0249 = s_13*IT_0247;
    const ccomplex_t IT_0250 = s_14*s_34;
    const ccomplex_t IT_0251 = (-2)*IT_0250;
    const ccomplex_t IT_0252 = IT_0249 + IT_0251;
    const ccomplex_t IT_0253 = 2*IT_0220;
    const ccomplex_t IT_0254 = 0.5*IT_0244;
    const ccomplex_t IT_0255 = (IT_0179*conj(IT_0179) + IT_0220*conj(IT_0220))
      *IT_0221 + conj(IT_0220)*IT_0222*IT_0223 + (IT_0220*IT_0222 + IT_0234
      *IT_0235 + IT_0243*IT_0244)*IT_0245 + (IT_0235*IT_0243 + IT_0234*IT_0244)
      *IT_0246 + 2*conj(IT_0243)*(IT_0234*IT_0248 + -IT_0243*IT_0252 + 0.5
      *IT_0235*IT_0253 + IT_0223*IT_0254) + 2*conj(IT_0234)*(0.5*IT_0223*IT_0235
       + IT_0243*IT_0248 + -IT_0234*IT_0252 + IT_0253*IT_0254);
    return create_ccomplex_return(IT_0255);
}

