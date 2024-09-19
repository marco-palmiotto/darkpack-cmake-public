#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_L_to_mu_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_L_to_mu_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
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
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0009, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = IT_0014*IT_0016;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0012 + IT_0015 + 
      -IT_0017 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = conj(N_B1)*e_em;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W1)*e_em;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_mu, 2);
    const ccomplex_t IT_0030 = pow(m_smu_L, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0020*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0037 = IT_0010*IT_0036;
    const ccomplex_t IT_0038 = IT_0014*IT_0036;
    const ccomplex_t IT_0039 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = IT_0014*IT_0039;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0037 + IT_0038 + 
      -IT_0040 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B4)*e_em;
    const ccomplex_t IT_0045 = IT_0008*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W4)*e_em;
    const ccomplex_t IT_0048 = IT_0013*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0043*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0019;
    const ccomplex_t IT_0057 = IT_0033*IT_0056;
    const ccomplex_t IT_0058 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = IT_0014*IT_0058;
    const ccomplex_t IT_0061 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0062 = IT_0010*IT_0061;
    const ccomplex_t IT_0063 = IT_0014*IT_0061;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0059 + IT_0060 + 
      -IT_0062 + -IT_0063);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = conj(N_B2)*e_em;
    const ccomplex_t IT_0067 = IT_0008*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W2)*e_em;
    const ccomplex_t IT_0070 = IT_0013*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0065*IT_0076;
    const ccomplex_t IT_0078 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0079 = IT_0010*IT_0078;
    const ccomplex_t IT_0080 = IT_0014*IT_0078;
    const ccomplex_t IT_0081 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0082 = IT_0010*IT_0081;
    const ccomplex_t IT_0083 = IT_0014*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B3)*e_em;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W3)*e_em;
    const ccomplex_t IT_0090 = IT_0013*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + IT_0091);
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0085*IT_0096;
    const ccomplex_t IT_0098 = IT_0055 + 0.5*IT_0057 + IT_0077 + IT_0097;
    const ccomplex_t IT_0099 = pow(m_N_1, 2);
    const ccomplex_t IT_0100 = cpow((-2)*s_12 + IT_0029 + -IT_0030 + -IT_0099 
      + -reg_prop, -1);
    const ccomplex_t IT_0101 = e_em*IT_0010;
    const ccomplex_t IT_0102 = e_em*IT_0014;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0101 + -IT_0102);
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0028*IT_0104;
    const ccomplex_t IT_0106 = IT_0100*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = conj(IT_0035) + -conj(IT_0098);
    const ccomplex_t IT_0109 = pow(m_W, -1);
    const ccomplex_t IT_0110 = cos(beta);
    const ccomplex_t IT_0111 = cpow(IT_0110, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0013*IT_0109*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0009;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = IT_0100*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0120 = IT_0010*IT_0119;
    const ccomplex_t IT_0121 = IT_0014*IT_0119;
    const ccomplex_t IT_0122 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0123 = IT_0010*IT_0122;
    const ccomplex_t IT_0124 = IT_0014*IT_0122;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0120 + IT_0121 + 
      -IT_0123 + -IT_0124);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0013*IT_0109*IT_0111;
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = IT_0094*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0126*IT_0130;
    const ccomplex_t IT_0132 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0133 = IT_0010*IT_0132;
    const ccomplex_t IT_0134 = IT_0014*IT_0132;
    const ccomplex_t IT_0135 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0136 = IT_0010*IT_0135;
    const ccomplex_t IT_0137 = IT_0014*IT_0135;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0133 + IT_0134 + 
      -IT_0136 + -IT_0137);
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0013*IT_0109*IT_0111;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0052*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0139*IT_0143;
    const ccomplex_t IT_0145 = IT_0031*IT_0113;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0020*IT_0146;
    const ccomplex_t IT_0148 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0149 = IT_0010*IT_0148;
    const ccomplex_t IT_0150 = IT_0014*IT_0148;
    const ccomplex_t IT_0151 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0152 = IT_0010*IT_0151;
    const ccomplex_t IT_0153 = IT_0014*IT_0151;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0149 + IT_0150 + 
      -IT_0152 + -IT_0153);
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0013*IT_0109*IT_0111;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = IT_0074*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0155*IT_0159;
    const ccomplex_t IT_0161 = -IT_0131 + -IT_0144 + (-0.5)*IT_0147 + -IT_0160;
    const ccomplex_t IT_0162 = IT_0056*IT_0146;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0161 + -IT_0163;
    const ccomplex_t IT_0165 = conj(IT_0161) + -conj(IT_0163);
    const ccomplex_t IT_0166 = pow(m_Z, 2);
    const ccomplex_t IT_0167 = s_13*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0004;
    const ccomplex_t IT_0169 = (-2)*IT_0000;
    const ccomplex_t IT_0170 = IT_0167 + IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = m_mu*m_N_1;
    const ccomplex_t IT_0172 = pow(s_14, 2);
    const ccomplex_t IT_0173 = IT_0001*IT_0172;
    const ccomplex_t IT_0174 = -IT_0099;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = IT_0171*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*e_em*(IT_0008*IT_0009 + 
      -IT_0007*IT_0013);
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = cpow((-2)*s_13 + IT_0029 + IT_0099 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0180 = IT_0028*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = IT_0178*IT_0181;
    const ccomplex_t IT_0183 = (-2)*IT_0055 + -IT_0057 + (-2)*IT_0077 + (-2)
      *IT_0097 + -IT_0182;
    const ccomplex_t IT_0184 = IT_0113*IT_0179;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0178*IT_0185;
    const ccomplex_t IT_0187 = 2*IT_0131 + 2*IT_0144 + IT_0147 + 2*IT_0160 +
       IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = (-2)*IT_0162;
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = 2*conj(IT_0187);
    const ccomplex_t IT_0192 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = -IT_0034;
    const ccomplex_t IT_0195 = 2*IT_0162;
    const ccomplex_t IT_0196 = (-2)*IT_0187;
    const ccomplex_t IT_0197 = IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = 2*conj(IT_0162);
    const ccomplex_t IT_0199 = (-2)*conj(IT_0187);
    const ccomplex_t IT_0200 = IT_0198 + IT_0199;
    const ccomplex_t IT_0201 = m_mu*IT_0175;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0203 = IT_0104*IT_0113;
    const ccomplex_t IT_0204 = IT_0100*IT_0202*IT_0203;
    const ccomplex_t IT_0205 = m_N_1*IT_0020;
    const ccomplex_t IT_0206 = m_N_1*IT_0056;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = IT_0033*IT_0208;
    const ccomplex_t IT_0210 = m_N_2*IT_0155;
    const ccomplex_t IT_0211 = m_N_1*IT_0065;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = IT_0076*IT_0213;
    const ccomplex_t IT_0215 = m_N_3*IT_0126;
    const ccomplex_t IT_0216 = m_N_1*IT_0085;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0215 + IT_0217;
    const ccomplex_t IT_0219 = IT_0096*IT_0218;
    const ccomplex_t IT_0220 = m_N_4*IT_0139;
    const ccomplex_t IT_0221 = m_N_1*IT_0043;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = IT_0220 + IT_0222;
    const ccomplex_t IT_0224 = IT_0054*IT_0223;
    const ccomplex_t IT_0225 = m_mu*IT_0117;
    const ccomplex_t IT_0226 = -IT_0204 + (-0.5)*IT_0209 + -IT_0214 + -IT_0219
       + -IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0033*(IT_0205 + -IT_0206);
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = conj(IT_0226) + conj(IT_0228);
    const ccomplex_t IT_0230 = IT_0226 + IT_0228;
    const ccomplex_t IT_0231 = IT_0028*IT_0114;
    const ccomplex_t IT_0232 = IT_0100*IT_0202*IT_0231;
    const ccomplex_t IT_0233 = m_mu*IT_0107;
    const ccomplex_t IT_0234 = IT_0146*IT_0208;
    const ccomplex_t IT_0235 = m_N_2*IT_0065;
    const ccomplex_t IT_0236 = m_N_1*IT_0155;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0159*IT_0238;
    const ccomplex_t IT_0240 = m_N_3*IT_0085;
    const ccomplex_t IT_0241 = m_N_1*IT_0126;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = IT_0130*IT_0243;
    const ccomplex_t IT_0245 = m_N_4*IT_0043;
    const ccomplex_t IT_0246 = m_N_1*IT_0139;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = IT_0143*IT_0248;
    const ccomplex_t IT_0250 = IT_0232 + -IT_0233 + (-0.5)*IT_0234 + IT_0239 +
       IT_0244 + IT_0249;
    const ccomplex_t IT_0251 = IT_0146*(IT_0205 + -IT_0206);
    const ccomplex_t IT_0252 = (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = conj(IT_0250) + conj(IT_0252);
    const ccomplex_t IT_0254 = 2*IT_0183;
    const ccomplex_t IT_0255 = (-2)*IT_0194;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = IT_0250 + IT_0252;
    const ccomplex_t IT_0258 = 2*conj(IT_0183);
    const ccomplex_t IT_0259 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0260 = IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = s_13 + IT_0262;
    const ccomplex_t IT_0264 = s_13*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0117 + -IT_0186;
    const ccomplex_t IT_0266 = 2*IT_0107 + IT_0182;
    const ccomplex_t IT_0267 = 2*IT_0194;
    const ccomplex_t IT_0268 = (-2)*IT_0183;
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = 2*conj(IT_0194);
    const ccomplex_t IT_0271 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = s_13*s_14;
    const ccomplex_t IT_0274 = s_34*IT_0099;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = conj(IT_0118) + conj(IT_0163);
    const ccomplex_t IT_0278 = IT_0118 + IT_0163;
    const ccomplex_t IT_0279 = conj(IT_0035) + conj(IT_0107);
    const ccomplex_t IT_0280 = IT_0035 + IT_0107;
    const ccomplex_t IT_0281 = s_34*IT_0166;
    const ccomplex_t IT_0282 = IT_0001*IT_0281;
    const ccomplex_t IT_0283 = 0.5*IT_0282;
    const ccomplex_t IT_0284 = s_34 + IT_0283;
    const ccomplex_t IT_0285 = m_N_1*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0226;
    const ccomplex_t IT_0288 = 2*IT_0228;
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = 2*conj(IT_0226);
    const ccomplex_t IT_0291 = 2*conj(IT_0228);
    const ccomplex_t IT_0292 = IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = 2*IT_0250;
    const ccomplex_t IT_0294 = 2*IT_0252;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = 2*conj(IT_0250);
    const ccomplex_t IT_0297 = 2*conj(IT_0252);
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = s_14*IT_0166;
    const ccomplex_t IT_0300 = IT_0001*IT_0299;
    const ccomplex_t IT_0301 = (-0.25)*IT_0300;
    const ccomplex_t IT_0302 = s_14 + IT_0301;
    const ccomplex_t IT_0303 = m_mu*IT_0302;
    const ccomplex_t IT_0304 = (-4)*IT_0303;
    const ccomplex_t IT_0305 = 0.5*IT_0300;
    const ccomplex_t IT_0306 = s_14 + IT_0305;
    const ccomplex_t IT_0307 = m_mu*IT_0306;
    const ccomplex_t IT_0308 = 2*IT_0307;
    const ccomplex_t IT_0309 = (-2)*IT_0226;
    const ccomplex_t IT_0310 = (-2)*IT_0228;
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0313 = (-2)*conj(IT_0228);
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = (-2)*IT_0250;
    const ccomplex_t IT_0316 = (-2)*IT_0252;
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = (-2)*conj(IT_0250);
    const ccomplex_t IT_0319 = (-2)*conj(IT_0252);
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = (-0.25)*IT_0282;
    const ccomplex_t IT_0322 = s_34 + IT_0321;
    const ccomplex_t IT_0323 = m_N_1*IT_0322;
    const ccomplex_t IT_0324 = (-4)*IT_0323;
    const ccomplex_t IT_0325 = IT_0166*IT_0171;
    const ccomplex_t IT_0326 = 3*IT_0325;
    const ccomplex_t IT_0327 = 2*IT_0161;
    const ccomplex_t IT_0328 = (-0.5)*IT_0327;
    const ccomplex_t IT_0329 = IT_0163 + IT_0328;
    const ccomplex_t IT_0330 = 2*conj(IT_0161);
    const ccomplex_t IT_0331 = (-0.5)*IT_0330;
    const ccomplex_t IT_0332 = conj(IT_0163) + IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0118;
    const ccomplex_t IT_0334 = 2*conj(IT_0118);
    const ccomplex_t IT_0335 = (-3)*IT_0325;
    const ccomplex_t IT_0336 = IT_0167 + IT_0169;
    const ccomplex_t IT_0337 = IT_0001*IT_0336;
    const ccomplex_t IT_0338 = (-0.5)*IT_0337;
    const ccomplex_t IT_0339 = s_13 + IT_0338;
    const ccomplex_t IT_0340 = 4*IT_0250;
    const ccomplex_t IT_0341 = (-6)*conj(IT_0226);
    const ccomplex_t IT_0342 = s_13*IT_0175;
    const ccomplex_t IT_0343 = IT_0171*IT_0263;
    const ccomplex_t IT_0344 = 2*IT_0265;
    const ccomplex_t IT_0345 = 2*conj(IT_0265);
    const ccomplex_t IT_0346 = (-2)*IT_0265;
    const ccomplex_t IT_0347 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0348 = pow(s_34, 2);
    const ccomplex_t IT_0349 = IT_0001*IT_0348;
    const ccomplex_t IT_0350 = -IT_0349;
    const ccomplex_t IT_0351 = IT_0029 + IT_0350;
    const ccomplex_t IT_0352 = IT_0171*IT_0351;
    const ccomplex_t IT_0353 = s_13*IT_0351;
    const ccomplex_t IT_0354 = (-2)*IT_0266;
    const ccomplex_t IT_0355 = s_14*IT_0029;
    const ccomplex_t IT_0356 = s_13*s_34;
    const ccomplex_t IT_0357 = -IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = (-2)*conj(IT_0266);
    const ccomplex_t IT_0360 = conj(IT_0098)*IT_0266 + IT_0098*conj(IT_0266) +
       0.5*conj(IT_0161)*IT_0344 + 0.5*IT_0161*IT_0345 + 0.5*IT_0277*IT_0346 +
       0.5*IT_0278*IT_0347 + 0.5*IT_0279*IT_0354 + 0.5*IT_0280*IT_0359;
    const ccomplex_t IT_0361 = m_mu*IT_0263;
    const ccomplex_t IT_0362 = m_N_1*IT_0263;
    const ccomplex_t IT_0363 = m_N_1*IT_0351;
    const ccomplex_t IT_0364 = 8*IT_0006*((IT_0035 + -IT_0098)*conj(IT_0107) +
       IT_0107*IT_0108 + -conj(IT_0118)*IT_0164 + -IT_0118*IT_0165) + (-4)*
      (IT_0107*conj(IT_0107) + (IT_0035 + -IT_0098)*IT_0108 + IT_0118*conj
      (IT_0118) + IT_0164*IT_0165)*IT_0170 + IT_0176*(conj(IT_0183)*IT_0190 +
       IT_0183*IT_0193 + conj(IT_0194)*IT_0197 + IT_0194*IT_0200) + IT_0201*
      (IT_0190*IT_0229 + IT_0193*IT_0230 + IT_0253*IT_0256 + IT_0257*IT_0260) +
       IT_0264*(IT_0200*IT_0265 + IT_0197*conj(IT_0265) + conj(IT_0266)*IT_0269 
      + IT_0266*IT_0272) + IT_0276*(conj(IT_0161)*IT_0190 + IT_0161*IT_0193 +
       conj(IT_0098)*IT_0256 + IT_0098*IT_0260 + IT_0197*IT_0277 + IT_0200
      *IT_0278 + IT_0269*IT_0279 + IT_0272*IT_0280) + IT_0286*(conj(IT_0107)
      *IT_0289 + IT_0107*IT_0292 + conj(IT_0118)*IT_0295 + IT_0118*IT_0298) + 
      (conj(IT_0118)*IT_0289 + IT_0118*IT_0292 + conj(IT_0107)*IT_0295 + IT_0107
      *IT_0298)*IT_0304 + IT_0308*(conj(IT_0161)*IT_0289 + IT_0161*IT_0292 +
       conj(IT_0098)*IT_0295 + IT_0098*IT_0298 + conj(IT_0163)*IT_0311 + IT_0163
      *IT_0314 + conj(IT_0035)*IT_0317 + IT_0035*IT_0320) + (conj(IT_0098)
      *IT_0289 + IT_0098*IT_0292 + conj(IT_0161)*IT_0295 + IT_0161*IT_0298 +
       conj(IT_0035)*IT_0311 + IT_0035*IT_0314 + conj(IT_0163)*IT_0317 + IT_0163
      *IT_0320)*IT_0324 + (-2)*IT_0326*(conj(IT_0035)*IT_0118 + IT_0035*conj
      (IT_0118) + conj(IT_0107)*IT_0329 + IT_0107*IT_0332 + (-0.5)*conj(IT_0098)
      *IT_0333 + (-0.5)*IT_0098*IT_0334) + (-2)*(conj(IT_0035)*IT_0164 + IT_0035
      *IT_0165 + conj(IT_0098)*IT_0329 + IT_0098*IT_0332 + (-0.5)*conj(IT_0107)
      *IT_0333 + (-0.5)*IT_0107*IT_0334)*IT_0335 + 4*IT_0339*(IT_0229*IT_0230 +
       IT_0253*(IT_0252 + 0.25*IT_0340)) + (-6)*IT_0171*(IT_0230*IT_0253 +
       IT_0257*(conj(IT_0228) + (-0.166666666666667)*IT_0341)) + (conj(IT_0187)
      *IT_0188 + conj(IT_0162)*IT_0197 + IT_0162*IT_0199 + conj(IT_0183)*IT_0254
       + conj(IT_0194)*IT_0269 + IT_0194*IT_0271)*IT_0342 + IT_0343*(IT_0200
      *IT_0266 + IT_0197*conj(IT_0266) + conj(IT_0194)*IT_0344 + IT_0194*IT_0345
       + conj(IT_0183)*IT_0346 + IT_0183*IT_0347) + (conj(IT_0266)*IT_0346 +
       IT_0266*IT_0347)*IT_0352 + IT_0353*(conj(IT_0265)*IT_0346 + conj(IT_0266)
      *IT_0354) + 2*IT_0358*IT_0360 + (IT_0229*IT_0346 + IT_0230*IT_0347 +
       IT_0253*IT_0354 + IT_0257*IT_0359)*IT_0361 + (IT_0197*IT_0253 + IT_0200
      *IT_0257 + conj(IT_0228)*IT_0269 + IT_0228*IT_0272 + conj(IT_0194)*IT_0287
       + IT_0194*IT_0290 + conj(IT_0183)*IT_0309 + IT_0183*IT_0312)*IT_0362 + 
      (conj(IT_0266)*IT_0309 + IT_0266*IT_0312 + IT_0253*IT_0346 + IT_0257
      *IT_0347 + conj(IT_0228)*IT_0354 + IT_0228*IT_0359)*IT_0363;
    return create_ccomplex_return(IT_0364);
}

