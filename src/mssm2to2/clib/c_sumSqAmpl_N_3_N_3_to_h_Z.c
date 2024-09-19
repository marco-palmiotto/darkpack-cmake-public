#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_h_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t IT_0000 = pow(m_N_3, 2);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0000;
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(alpha);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + IT_0016 + 
      -IT_0021 + -IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0009*IT_0017;
    const ccomplex_t IT_0028 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0008*IT_0018;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0029 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0026*IT_0036;
    const ccomplex_t IT_0038 = pow(m_h, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0038 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0051 = IT_0014*IT_0050;
    const ccomplex_t IT_0052 = IT_0045 + IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0054 = IT_0019*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0056 = IT_0022*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0058 = IT_0019*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0060 = IT_0022*IT_0059;
    const ccomplex_t IT_0061 = -IT_0054 + -IT_0056 + -IT_0058 + -IT_0060;
    const ccomplex_t IT_0062 = IT_0052 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0066 = IT_0027*IT_0065;
    const ccomplex_t IT_0067 = IT_0030*IT_0065;
    const ccomplex_t IT_0068 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0069 = IT_0027*IT_0068;
    const ccomplex_t IT_0070 = IT_0030*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = IT_0043*IT_0073;
    const ccomplex_t IT_0075 = cos(beta);
    const ccomplex_t IT_0076 = IT_0009*IT_0075;
    const ccomplex_t IT_0077 = IT_0011*IT_0076;
    const ccomplex_t IT_0078 = sin(beta);
    const ccomplex_t IT_0079 = IT_0009*IT_0078;
    const ccomplex_t IT_0080 = IT_0015*IT_0079;
    const ccomplex_t IT_0081 = IT_0018*IT_0075;
    const ccomplex_t IT_0082 = IT_0020*IT_0081;
    const ccomplex_t IT_0083 = IT_0018*IT_0078;
    const ccomplex_t IT_0084 = IT_0023*IT_0083;
    const ccomplex_t IT_0085 = IT_0077 + -IT_0080 + -IT_0082 + IT_0084;
    const ccomplex_t IT_0086 = e_em*IT_0009*IT_0018*(IT_0007*IT_0075 + IT_0013
      *IT_0078);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0015*IT_0076;
    const ccomplex_t IT_0093 = IT_0011*IT_0079;
    const ccomplex_t IT_0094 = IT_0023*IT_0081;
    const ccomplex_t IT_0095 = IT_0020*IT_0083;
    const ccomplex_t IT_0096 = IT_0092 + IT_0093 + -IT_0094 + -IT_0095;
    const ccomplex_t IT_0097 = -e_em*IT_0009*IT_0018*(IT_0013*IT_0075 + 
      -IT_0007*IT_0078);
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0096*IT_0098;
    const ccomplex_t IT_0100 = pow(m_Z, 2);
    const ccomplex_t IT_0101 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0100 + 
      -reg_prop, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0099*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0107 = IT_0010*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0109 = IT_0014*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0111 = IT_0010*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0113 = IT_0014*IT_0112;
    const ccomplex_t IT_0114 = IT_0107 + IT_0109 + IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0116 = IT_0019*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0118 = IT_0022*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0120 = IT_0019*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0122 = IT_0022*IT_0121;
    const ccomplex_t IT_0123 = -IT_0116 + -IT_0118 + -IT_0120 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0114 + IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0128 = IT_0027*IT_0127;
    const ccomplex_t IT_0129 = IT_0030*IT_0127;
    const ccomplex_t IT_0130 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0131 = IT_0027*IT_0130;
    const ccomplex_t IT_0132 = IT_0030*IT_0130;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0128 + IT_0129 + 
      -IT_0131 + -IT_0132);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = IT_0126*IT_0135;
    const ccomplex_t IT_0137 = IT_0105*IT_0136;
    const ccomplex_t IT_0138 = -IT_0041 + (-2)*IT_0074 + 0.5*IT_0091 + 0.5
      *IT_0103 + (-2)*IT_0137;
    const ccomplex_t IT_0139 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0140 = IT_0010*IT_0139;
    const ccomplex_t IT_0141 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0142 = IT_0014*IT_0141;
    const ccomplex_t IT_0143 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0144 = IT_0019*IT_0143;
    const ccomplex_t IT_0145 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0146 = IT_0022*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0140 + IT_0142 + 
      -IT_0144 + -IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = IT_0036*IT_0148;
    const ccomplex_t IT_0150 = IT_0040*IT_0149;
    const ccomplex_t IT_0151 = IT_0076*IT_0139;
    const ccomplex_t IT_0152 = IT_0079*IT_0141;
    const ccomplex_t IT_0153 = IT_0081*IT_0143;
    const ccomplex_t IT_0154 = IT_0083*IT_0145;
    const ccomplex_t IT_0155 = -IT_0151 + IT_0152 + IT_0153 + -IT_0154;
    const ccomplex_t IT_0156 = IT_0087*IT_0155;
    const ccomplex_t IT_0157 = IT_0090*IT_0156;
    const ccomplex_t IT_0158 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0159 = IT_0010*IT_0158;
    const ccomplex_t IT_0160 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0161 = IT_0014*IT_0160;
    const ccomplex_t IT_0162 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0163 = IT_0010*IT_0162;
    const ccomplex_t IT_0164 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0165 = IT_0014*IT_0164;
    const ccomplex_t IT_0166 = IT_0159 + IT_0161 + IT_0163 + IT_0165;
    const ccomplex_t IT_0167 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0168 = IT_0019*IT_0167;
    const ccomplex_t IT_0169 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0170 = IT_0022*IT_0169;
    const ccomplex_t IT_0171 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0172 = IT_0019*IT_0171;
    const ccomplex_t IT_0173 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0174 = IT_0022*IT_0173;
    const ccomplex_t IT_0175 = -IT_0168 + -IT_0170 + -IT_0172 + -IT_0174;
    const ccomplex_t IT_0176 = IT_0166 + IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0180 = IT_0027*IT_0179;
    const ccomplex_t IT_0181 = IT_0030*IT_0179;
    const ccomplex_t IT_0182 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0183 = IT_0027*IT_0182;
    const ccomplex_t IT_0184 = IT_0030*IT_0182;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0180 + IT_0181 + 
      -IT_0183 + -IT_0184);
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0178*IT_0186;
    const ccomplex_t IT_0188 = cpow((-2)*s_23 + IT_0000 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0187*IT_0189;
    const ccomplex_t IT_0191 = IT_0076*IT_0141;
    const ccomplex_t IT_0192 = IT_0079*IT_0139;
    const ccomplex_t IT_0193 = IT_0081*IT_0145;
    const ccomplex_t IT_0194 = IT_0083*IT_0143;
    const ccomplex_t IT_0195 = -IT_0191 + -IT_0192 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = IT_0098*IT_0195;
    const ccomplex_t IT_0197 = IT_0102*IT_0196;
    const ccomplex_t IT_0198 = -IT_0150 + (-0.5)*IT_0157 + (-2)*IT_0190 + (
      -0.5)*IT_0197;
    const ccomplex_t IT_0199 = 2*IT_0198;
    const ccomplex_t IT_0200 = 0.5*IT_0035;
    const ccomplex_t IT_0201 = IT_0148*IT_0200;
    const ccomplex_t IT_0202 = IT_0040*IT_0201;
    const ccomplex_t IT_0203 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0204 = IT_0010*IT_0203;
    const ccomplex_t IT_0205 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0206 = IT_0014*IT_0205;
    const ccomplex_t IT_0207 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0208 = IT_0010*IT_0207;
    const ccomplex_t IT_0209 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0210 = IT_0014*IT_0209;
    const ccomplex_t IT_0211 = IT_0204 + IT_0206 + IT_0208 + IT_0210;
    const ccomplex_t IT_0212 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0213 = IT_0019*IT_0212;
    const ccomplex_t IT_0214 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0215 = IT_0022*IT_0214;
    const ccomplex_t IT_0216 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0217 = IT_0019*IT_0216;
    const ccomplex_t IT_0218 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0219 = IT_0022*IT_0218;
    const ccomplex_t IT_0220 = -IT_0213 + -IT_0215 + -IT_0217 + -IT_0219;
    const ccomplex_t IT_0221 = IT_0211 + IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0225 = IT_0030*IT_0224;
    const ccomplex_t IT_0226 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0227 = IT_0027*IT_0226;
    const ccomplex_t IT_0228 = IT_0030*IT_0226;
    const ccomplex_t IT_0229 = IT_0027*IT_0224;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*(IT_0225 + -IT_0227 + 
      -IT_0228 + IT_0229);
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = IT_0223*IT_0231;
    const ccomplex_t IT_0233 = IT_0105*IT_0232;
    const ccomplex_t IT_0234 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0235 = IT_0010*IT_0234;
    const ccomplex_t IT_0236 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0237 = IT_0014*IT_0236;
    const ccomplex_t IT_0238 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0239 = IT_0010*IT_0238;
    const ccomplex_t IT_0240 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0241 = IT_0014*IT_0240;
    const ccomplex_t IT_0242 = IT_0235 + IT_0237 + IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0244 = IT_0019*IT_0243;
    const ccomplex_t IT_0245 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0246 = IT_0022*IT_0245;
    const ccomplex_t IT_0247 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0248 = IT_0019*IT_0247;
    const ccomplex_t IT_0249 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0250 = IT_0022*IT_0249;
    const ccomplex_t IT_0251 = -IT_0244 + -IT_0246 + -IT_0248 + -IT_0250;
    const ccomplex_t IT_0252 = IT_0242 + IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = (-0.5)*IT_0253;
    const ccomplex_t IT_0255 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0256 = IT_0027*IT_0255;
    const ccomplex_t IT_0257 = IT_0030*IT_0255;
    const ccomplex_t IT_0258 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0259 = IT_0027*IT_0258;
    const ccomplex_t IT_0260 = IT_0030*IT_0258;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*(IT_0256 + IT_0257 + 
      -IT_0259 + -IT_0260);
    const ccomplex_t IT_0262 = (-0.5)*IT_0261;
    const ccomplex_t IT_0263 = IT_0254*IT_0262;
    const ccomplex_t IT_0264 = IT_0043*IT_0263;
    const ccomplex_t IT_0265 = 0.5*IT_0157 + 0.5*IT_0197 + -IT_0202 + (-2)
      *IT_0233 + (-2)*IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = IT_0199 + IT_0266;
    const ccomplex_t IT_0268 = 2*conj(IT_0198);
    const ccomplex_t IT_0269 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = IT_0026*IT_0200;
    const ccomplex_t IT_0272 = IT_0040*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0274 = IT_0010*IT_0273;
    const ccomplex_t IT_0275 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0276 = IT_0014*IT_0275;
    const ccomplex_t IT_0277 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0278 = IT_0010*IT_0277;
    const ccomplex_t IT_0279 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0280 = IT_0014*IT_0279;
    const ccomplex_t IT_0281 = IT_0274 + IT_0276 + IT_0278 + IT_0280;
    const ccomplex_t IT_0282 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0283 = IT_0019*IT_0282;
    const ccomplex_t IT_0284 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0285 = IT_0022*IT_0284;
    const ccomplex_t IT_0286 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0287 = IT_0019*IT_0286;
    const ccomplex_t IT_0288 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0289 = IT_0022*IT_0288;
    const ccomplex_t IT_0290 = -IT_0283 + -IT_0285 + -IT_0287 + -IT_0289;
    const ccomplex_t IT_0291 = IT_0281 + IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = (-0.5)*IT_0292;
    const ccomplex_t IT_0294 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0295 = IT_0027*IT_0294;
    const ccomplex_t IT_0296 = IT_0030*IT_0294;
    const ccomplex_t IT_0297 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0298 = IT_0027*IT_0297;
    const ccomplex_t IT_0299 = IT_0030*IT_0297;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*(IT_0295 + IT_0296 + 
      -IT_0298 + -IT_0299);
    const ccomplex_t IT_0301 = (-0.5)*IT_0300;
    const ccomplex_t IT_0302 = IT_0293*IT_0301;
    const ccomplex_t IT_0303 = IT_0189*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0091 + (-0.5)*IT_0103 + -IT_0272 + (
      -2)*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0265;
    const ccomplex_t IT_0306 = (-2)*IT_0198;
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = 2*conj(IT_0265);
    const ccomplex_t IT_0309 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = s_12 + IT_0312;
    const ccomplex_t IT_0314 = s_12*IT_0313;
    const ccomplex_t IT_0315 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*IT_0315;
    const ccomplex_t IT_0317 = IT_0187*IT_0316;
    const ccomplex_t IT_0318 = cpow((-2)*s_13 + IT_0038 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = IT_0149*IT_0319;
    const ccomplex_t IT_0321 = 0.5*IT_0157 + 0.5*IT_0197 + 2*IT_0317 + IT_0320;
    const ccomplex_t IT_0322 = IT_0201*IT_0319;
    const ccomplex_t IT_0323 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*IT_0323;
    const ccomplex_t IT_0325 = IT_0263*IT_0324;
    const ccomplex_t IT_0326 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*IT_0326;
    const ccomplex_t IT_0328 = IT_0232*IT_0327;
    const ccomplex_t IT_0329 = (-0.5)*IT_0157 + (-0.5)*IT_0197 + IT_0322 + 2
      *IT_0325 + 2*IT_0328;
    const ccomplex_t IT_0330 = IT_0271*IT_0319;
    const ccomplex_t IT_0331 = IT_0302*IT_0316;
    const ccomplex_t IT_0332 = 0.5*IT_0091 + 0.5*IT_0103 + IT_0330 + 2*IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0304;
    const ccomplex_t IT_0334 = (-2)*IT_0138;
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = 2*conj(IT_0304);
    const ccomplex_t IT_0337 = (-2)*conj(IT_0138);
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = IT_0037*IT_0319;
    const ccomplex_t IT_0340 = IT_0073*IT_0324;
    const ccomplex_t IT_0341 = IT_0136*IT_0327;
    const ccomplex_t IT_0342 = (-0.5)*IT_0091 + (-0.5)*IT_0103 + IT_0339 + 2
      *IT_0340 + 2*IT_0341;
    const ccomplex_t IT_0343 = 2*IT_0138;
    const ccomplex_t IT_0344 = (-2)*IT_0304;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = 2*conj(IT_0138);
    const ccomplex_t IT_0347 = (-2)*conj(IT_0304);
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = s_12*s_14;
    const ccomplex_t IT_0350 = s_24*IT_0000;
    const ccomplex_t IT_0351 = -IT_0350;
    const ccomplex_t IT_0352 = IT_0349 + IT_0351;
    const ccomplex_t IT_0353 = 0.5*IT_0150 + IT_0190;
    const ccomplex_t IT_0354 = -IT_0317 + (-0.5)*IT_0320;
    const ccomplex_t IT_0355 = (-0.5)*IT_0322 + -IT_0325 + -IT_0328;
    const ccomplex_t IT_0356 = 0.5*IT_0202 + IT_0233 + IT_0264;
    const ccomplex_t IT_0357 = conj(IT_0355) + conj(IT_0356);
    const ccomplex_t IT_0358 = IT_0355 + IT_0356;
    const ccomplex_t IT_0359 = 0.5*IT_0272 + IT_0303;
    const ccomplex_t IT_0360 = -IT_0331;
    const ccomplex_t IT_0361 = (-0.5)*IT_0330;
    const ccomplex_t IT_0362 = conj(IT_0360) + conj(IT_0361);
    const ccomplex_t IT_0363 = IT_0360 + IT_0361;
    const ccomplex_t IT_0364 = IT_0359 + IT_0363;
    const ccomplex_t IT_0365 = 0.5*IT_0041;
    const ccomplex_t IT_0366 = (-0.5)*IT_0339 + -IT_0340 + -IT_0341;
    const ccomplex_t IT_0367 = IT_0074 + IT_0137;
    const ccomplex_t IT_0368 = IT_0270*(IT_0353 + IT_0354) + IT_0267*(conj
      (IT_0353) + conj(IT_0354)) + IT_0307*IT_0357 + IT_0310*IT_0358 + IT_0335*
      (conj(IT_0359) + IT_0362) + IT_0338*IT_0364 + IT_0348*(IT_0365 + IT_0366 +
       IT_0367) + IT_0345*(conj(IT_0365) + conj(IT_0366) + conj(IT_0367));
    const ccomplex_t IT_0369 = m_N_3*IT_0005;
    const ccomplex_t IT_0370 = m_W*e_em;
    const ccomplex_t IT_0371 = cpow(IT_0008, -4);
    const ccomplex_t IT_0372 = IT_0007*IT_0017*IT_0078*IT_0371;
    const ccomplex_t IT_0373 = IT_0013*IT_0017*IT_0075*IT_0371;
    const ccomplex_t IT_0374 = cpow(IT_0017, -3);
    const ccomplex_t IT_0375 = IT_0007*IT_0078*IT_0374;
    const ccomplex_t IT_0376 = IT_0013*IT_0075*IT_0374;
    const ccomplex_t IT_0377 = cpow(IT_0008, -2);
    const ccomplex_t IT_0378 = IT_0007*IT_0018*IT_0078*IT_0377;
    const ccomplex_t IT_0379 = IT_0013*IT_0018*IT_0075*IT_0377;
    const ccomplex_t IT_0380 = cpow(IT_0017, -2);
    const ccomplex_t IT_0381 = cpow(IT_0377 + IT_0380, -1);
    const ccomplex_t IT_0382 = (IT_0372 + -IT_0373 + IT_0375 + -IT_0376 + 2
      *IT_0378 + (-2)*IT_0379)*IT_0381;
    const ccomplex_t IT_0383 = IT_0370*IT_0382;
    const ccomplex_t IT_0384 = (0 + _Complex_I*1)*IT_0383;
    const ccomplex_t IT_0385 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0386 = IT_0384*IT_0385;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*IT_0386;
    const ccomplex_t IT_0388 = IT_0036*IT_0387;
    const ccomplex_t IT_0389 = m_N_3*IT_0271;
    const ccomplex_t IT_0390 = m_N_3*IT_0037;
    const ccomplex_t IT_0391 = -IT_0390;
    const ccomplex_t IT_0392 = IT_0389 + IT_0391;
    const ccomplex_t IT_0393 = IT_0040*IT_0392;
    const ccomplex_t IT_0394 = m_N_3*IT_0149;
    const ccomplex_t IT_0395 = m_N_3*IT_0201;
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = IT_0394 + IT_0396;
    const ccomplex_t IT_0398 = IT_0319*IT_0397;
    const ccomplex_t IT_0399 = IT_0072*IT_0254;
    const ccomplex_t IT_0400 = m_N_1*IT_0399;
    const ccomplex_t IT_0401 = m_N_3*IT_0263;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = IT_0324*IT_0403;
    const ccomplex_t IT_0405 = IT_0135*IT_0223;
    const ccomplex_t IT_0406 = m_N_2*IT_0405;
    const ccomplex_t IT_0407 = m_N_3*IT_0232;
    const ccomplex_t IT_0408 = -IT_0407;
    const ccomplex_t IT_0409 = IT_0406 + IT_0408;
    const ccomplex_t IT_0410 = IT_0327*IT_0409;
    const ccomplex_t IT_0411 = IT_0186*IT_0293;
    const ccomplex_t IT_0412 = m_N_4*IT_0411;
    const ccomplex_t IT_0413 = m_N_3*IT_0302;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = IT_0189*IT_0415;
    const ccomplex_t IT_0417 = (-0.5)*IT_0388 + 0.5*IT_0393 + (-0.5)*IT_0398 +
       -IT_0404 + -IT_0410 + -IT_0416;
    const ccomplex_t IT_0418 = IT_0200*IT_0387;
    const ccomplex_t IT_0419 = IT_0064*IT_0262;
    const ccomplex_t IT_0420 = m_N_1*IT_0419;
    const ccomplex_t IT_0421 = m_N_3*IT_0073;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0420 + IT_0422;
    const ccomplex_t IT_0424 = IT_0043*IT_0423;
    const ccomplex_t IT_0425 = IT_0126*IT_0231;
    const ccomplex_t IT_0426 = m_N_2*IT_0425;
    const ccomplex_t IT_0427 = m_N_3*IT_0136;
    const ccomplex_t IT_0428 = -IT_0427;
    const ccomplex_t IT_0429 = IT_0426 + IT_0428;
    const ccomplex_t IT_0430 = IT_0105*IT_0429;
    const ccomplex_t IT_0431 = IT_0178*IT_0301;
    const ccomplex_t IT_0432 = m_N_4*IT_0431;
    const ccomplex_t IT_0433 = m_N_3*IT_0187;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = IT_0316*IT_0435;
    const ccomplex_t IT_0437 = 0.5*IT_0393 + (-0.5)*IT_0398 + 0.5*IT_0418 +
       IT_0424 + IT_0430 + IT_0436;
    const ccomplex_t IT_0438 = conj(IT_0417) + conj(IT_0437);
    const ccomplex_t IT_0439 = IT_0417 + IT_0437;
    const ccomplex_t IT_0440 = IT_0040*IT_0397;
    const ccomplex_t IT_0441 = IT_0319*IT_0392;
    const ccomplex_t IT_0442 = IT_0324*IT_0423;
    const ccomplex_t IT_0443 = IT_0327*IT_0429;
    const ccomplex_t IT_0444 = IT_0189*IT_0435;
    const ccomplex_t IT_0445 = (-0.5)*IT_0418 + 0.5*IT_0440 + (-0.5)*IT_0441 +
       -IT_0442 + -IT_0443 + -IT_0444;
    const ccomplex_t IT_0446 = IT_0043*IT_0403;
    const ccomplex_t IT_0447 = IT_0105*IT_0409;
    const ccomplex_t IT_0448 = IT_0316*IT_0415;
    const ccomplex_t IT_0449 = 0.5*IT_0388 + 0.5*IT_0440 + (-0.5)*IT_0441 +
       IT_0446 + IT_0447 + IT_0448;
    const ccomplex_t IT_0450 = conj(IT_0445) + conj(IT_0449);
    const ccomplex_t IT_0451 = IT_0445 + IT_0449;
    const ccomplex_t IT_0452 = IT_0000*IT_0313;
    const ccomplex_t IT_0453 = 2*IT_0329;
    const ccomplex_t IT_0454 = (-2)*IT_0321;
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = 2*conj(IT_0329);
    const ccomplex_t IT_0457 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0458 = IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = 2*IT_0321;
    const ccomplex_t IT_0460 = (-2)*IT_0329;
    const ccomplex_t IT_0461 = IT_0459 + IT_0460;
    const ccomplex_t IT_0462 = 2*conj(IT_0321);
    const ccomplex_t IT_0463 = (-2)*conj(IT_0329);
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = pow(s_24, 2);
    const ccomplex_t IT_0466 = IT_0001*IT_0465;
    const ccomplex_t IT_0467 = -IT_0466;
    const ccomplex_t IT_0468 = IT_0000 + IT_0467;
    const ccomplex_t IT_0469 = IT_0000*IT_0468;
    const ccomplex_t IT_0470 = -IT_0469;
    const ccomplex_t IT_0471 = s_14*IT_0100;
    const ccomplex_t IT_0472 = IT_0001*IT_0471;
    const ccomplex_t IT_0473 = (-0.25)*IT_0472;
    const ccomplex_t IT_0474 = s_14 + IT_0473;
    const ccomplex_t IT_0475 = m_N_3*IT_0474;
    const ccomplex_t IT_0476 = (-4)*IT_0475;
    const ccomplex_t IT_0477 = 2*IT_0417;
    const ccomplex_t IT_0478 = 2*IT_0437;
    const ccomplex_t IT_0479 = IT_0477 + IT_0478;
    const ccomplex_t IT_0480 = 2*conj(IT_0417);
    const ccomplex_t IT_0481 = 2*conj(IT_0437);
    const ccomplex_t IT_0482 = IT_0480 + IT_0481;
    const ccomplex_t IT_0483 = 2*IT_0445;
    const ccomplex_t IT_0484 = 2*IT_0449;
    const ccomplex_t IT_0485 = IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = 2*conj(IT_0445);
    const ccomplex_t IT_0487 = 2*conj(IT_0449);
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = (-2)*IT_0417;
    const ccomplex_t IT_0490 = (-2)*IT_0437;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = (-2)*conj(IT_0417);
    const ccomplex_t IT_0493 = (-2)*conj(IT_0437);
    const ccomplex_t IT_0494 = IT_0492 + IT_0493;
    const ccomplex_t IT_0495 = (-2)*IT_0445;
    const ccomplex_t IT_0496 = (-2)*IT_0449;
    const ccomplex_t IT_0497 = IT_0495 + IT_0496;
    const ccomplex_t IT_0498 = (-2)*conj(IT_0445);
    const ccomplex_t IT_0499 = (-2)*conj(IT_0449);
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = 0.5*IT_0472;
    const ccomplex_t IT_0502 = s_14 + IT_0501;
    const ccomplex_t IT_0503 = m_N_3*IT_0502;
    const ccomplex_t IT_0504 = 2*IT_0503;
    const ccomplex_t IT_0505 = conj(IT_0365) + conj(IT_0367);
    const ccomplex_t IT_0506 = IT_0365 + IT_0367;
    const ccomplex_t IT_0507 = s_24*IT_0100;
    const ccomplex_t IT_0508 = IT_0001*IT_0507;
    const ccomplex_t IT_0509 = 0.5*IT_0508;
    const ccomplex_t IT_0510 = s_24 + IT_0509;
    const ccomplex_t IT_0511 = m_N_3*IT_0510;
    const ccomplex_t IT_0512 = 2*IT_0511;
    const ccomplex_t IT_0513 = (-0.25)*IT_0508;
    const ccomplex_t IT_0514 = s_24 + IT_0513;
    const ccomplex_t IT_0515 = m_N_3*IT_0514;
    const ccomplex_t IT_0516 = (-4)*IT_0515;
    const ccomplex_t IT_0517 = s_12*IT_0100;
    const ccomplex_t IT_0518 = s_14*s_24;
    const ccomplex_t IT_0519 = (-2)*IT_0518;
    const ccomplex_t IT_0520 = IT_0517 + IT_0519;
    const ccomplex_t IT_0521 = IT_0001*IT_0520;
    const ccomplex_t IT_0522 = (-0.5)*IT_0521;
    const ccomplex_t IT_0523 = s_12 + IT_0522;
    const ccomplex_t IT_0524 = 4*IT_0445;
    const ccomplex_t IT_0525 = 6*conj(IT_0417);
    const ccomplex_t IT_0526 = pow(m_Z, 4);
    const ccomplex_t IT_0527 = s_12*IT_0526;
    const ccomplex_t IT_0528 = IT_0001*IT_0527;
    const ccomplex_t IT_0529 = (-0.25)*IT_0528;
    const ccomplex_t IT_0530 = IT_0518 + IT_0529;
    const ccomplex_t IT_0531 = -conj(IT_0354);
    const ccomplex_t IT_0532 = conj(IT_0355) + IT_0531;
    const ccomplex_t IT_0533 = -IT_0356;
    const ccomplex_t IT_0534 = IT_0353 + IT_0533;
    const ccomplex_t IT_0535 = -IT_0354 + IT_0355;
    const ccomplex_t IT_0536 = -conj(IT_0356);
    const ccomplex_t IT_0537 = conj(IT_0353) + IT_0536;
    const ccomplex_t IT_0538 = 8*IT_0366;
    const ccomplex_t IT_0539 = 8*IT_0359;
    const ccomplex_t IT_0540 = 8*conj(IT_0359);
    const ccomplex_t IT_0541 = 8*conj(IT_0366);
    const ccomplex_t IT_0542 = (-8)*conj(IT_0360);
    const ccomplex_t IT_0543 = conj(IT_0359)*IT_0366 + IT_0359*conj(IT_0366) +
       IT_0532*IT_0534 + IT_0535*IT_0537 + (-0.125)*IT_0505*IT_0538 + (-0.125)
      *conj(IT_0360)*IT_0539 + conj(IT_0361)*(IT_0506 + (-0.125)*IT_0539) +
       IT_0363*(IT_0505 + (-0.125)*IT_0540) + (-0.125)*IT_0506*(IT_0541 +
       IT_0542);
    const ccomplex_t IT_0544 = (-0.5)*IT_0528;
    const ccomplex_t IT_0545 = IT_0517 + IT_0519 + IT_0544;
    const ccomplex_t IT_0546 = 4*IT_0366;
    const ccomplex_t IT_0547 = 4*conj(IT_0366);
    const ccomplex_t IT_0548 = 4*conj(IT_0359);
    const ccomplex_t IT_0549 = 4*IT_0359;
    const ccomplex_t IT_0550 = (-4)*IT_0367;
    const ccomplex_t IT_0551 = IT_0359*conj(IT_0359) + IT_0366*conj(IT_0366) +
       IT_0532*IT_0535 + IT_0534*IT_0537 + (-0.25)*IT_0362*IT_0546 + IT_0363*
      (IT_0362 + (-0.25)*IT_0547) + (-0.25)*IT_0367*IT_0548 + IT_0365*(IT_0505 +
       (-0.25)*IT_0548) + (-0.25)*IT_0505*(IT_0549 + IT_0550);
    const ccomplex_t IT_0552 = m_N_3*IT_0468;
    const ccomplex_t IT_0553 = -IT_0552;
    const ccomplex_t IT_0554 = 2*IT_0342;
    const ccomplex_t IT_0555 = (-2)*IT_0332;
    const ccomplex_t IT_0556 = IT_0554 + IT_0555;
    const ccomplex_t IT_0557 = 2*conj(IT_0342);
    const ccomplex_t IT_0558 = (-2)*conj(IT_0332);
    const ccomplex_t IT_0559 = IT_0557 + IT_0558;
    const ccomplex_t IT_0560 = s_12*IT_0005;
    const ccomplex_t IT_0561 = s_12*IT_0468;
    const ccomplex_t IT_0562 = -IT_0561;
    const ccomplex_t IT_0563 = 2*IT_0332;
    const ccomplex_t IT_0564 = (-2)*IT_0342;
    const ccomplex_t IT_0565 = IT_0563 + IT_0564;
    const ccomplex_t IT_0566 = (-2)*conj(IT_0342);
    const ccomplex_t IT_0567 = s_12*s_24;
    const ccomplex_t IT_0568 = s_14*IT_0000;
    const ccomplex_t IT_0569 = -IT_0568;
    const ccomplex_t IT_0570 = IT_0567 + IT_0569;
    const ccomplex_t IT_0571 = 2*conj(IT_0332);
    const ccomplex_t IT_0572 = IT_0357*IT_0455 + IT_0358*IT_0458 + (conj
      (IT_0353) + conj(IT_0354))*IT_0461 + (IT_0353 + IT_0354)*IT_0464 + (conj
      (IT_0365) + conj(IT_0366) + conj(IT_0367))*IT_0556 + (IT_0365 + IT_0366 +
       IT_0367)*IT_0559 + (conj(IT_0359) + IT_0362)*IT_0565 + IT_0364*(IT_0566 +
       IT_0571);
    const ccomplex_t IT_0573 = m_N_3*IT_0313;
    const ccomplex_t IT_0574 = conj(IT_0437)*(IT_0343 + IT_0344 + IT_0453 +
       IT_0454) + conj(IT_0417)*IT_0455 + IT_0437*(IT_0346 + IT_0347 + IT_0456 +
       IT_0457) + IT_0417*IT_0458 + conj(IT_0138)*IT_0477 + IT_0138*IT_0480 +
       conj(IT_0304)*IT_0489 + IT_0304*IT_0492 + IT_0450*(IT_0305 + IT_0306 +
       IT_0554 + IT_0555) + IT_0451*(IT_0308 + IT_0309 + IT_0557 + IT_0558);
    const ccomplex_t IT_0575 = IT_0000*IT_0100;
    const ccomplex_t IT_0576 = 3*IT_0575;
    const ccomplex_t IT_0577 = 2*IT_0355;
    const ccomplex_t IT_0578 = (-2)*IT_0354;
    const ccomplex_t IT_0579 = IT_0577 + IT_0578;
    const ccomplex_t IT_0580 = 2*conj(IT_0355);
    const ccomplex_t IT_0581 = (-2)*conj(IT_0354);
    const ccomplex_t IT_0582 = IT_0580 + IT_0581;
    const ccomplex_t IT_0583 = 2*IT_0353;
    const ccomplex_t IT_0584 = (-2)*IT_0356;
    const ccomplex_t IT_0585 = IT_0583 + IT_0584;
    const ccomplex_t IT_0586 = 2*conj(IT_0353);
    const ccomplex_t IT_0587 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0588 = IT_0586 + IT_0587;
    const ccomplex_t IT_0589 = 2*IT_0366;
    const ccomplex_t IT_0590 = 2*conj(IT_0366);
    const ccomplex_t IT_0591 = 2*IT_0359;
    const ccomplex_t IT_0592 = 2*conj(IT_0359);
    const ccomplex_t IT_0593 = (-2)*IT_0355;
    const ccomplex_t IT_0594 = (-2)*conj(IT_0355);
    const ccomplex_t IT_0595 = (-2)*IT_0353;
    const ccomplex_t IT_0596 = (-2)*conj(IT_0353);
    const ccomplex_t IT_0597 = IT_0505*IT_0579 + IT_0506*IT_0582 + IT_0362
      *IT_0585 + IT_0363*IT_0588 + conj(IT_0356)*IT_0589 + IT_0356*IT_0590 +
       conj(IT_0354)*IT_0591 + IT_0354*IT_0592 + conj(IT_0359)*IT_0593 + IT_0359
      *IT_0594 + conj(IT_0366)*IT_0595 + IT_0366*IT_0596;
    const ccomplex_t IT_0598 = (-3)*IT_0575;
    const ccomplex_t IT_0599 = IT_0362*IT_0579 + IT_0363*IT_0582 + IT_0505
      *IT_0585 + IT_0506*IT_0588 + conj(IT_0354)*IT_0589 + IT_0354*IT_0590 +
       conj(IT_0356)*IT_0591 + IT_0356*IT_0592 + conj(IT_0366)*IT_0593 + IT_0366
      *IT_0594 + conj(IT_0359)*IT_0595 + IT_0359*IT_0596;
    const ccomplex_t IT_0600 = IT_0006*(conj(IT_0138)*IT_0267 + IT_0138
      *IT_0270 + conj(IT_0304)*IT_0307 + IT_0304*IT_0310) + IT_0314*(IT_0270
      *IT_0321 + IT_0267*conj(IT_0321) + IT_0310*IT_0329 + IT_0307*conj(IT_0329)
       + conj(IT_0332)*IT_0335 + IT_0332*IT_0338 + conj(IT_0342)*IT_0345 +
       IT_0342*IT_0348) + IT_0352*IT_0368 + IT_0369*(IT_0307*IT_0438 + IT_0310
      *IT_0439 + IT_0345*IT_0450 + IT_0348*IT_0451) + IT_0452*(IT_0310*IT_0332 +
       IT_0307*conj(IT_0332) + IT_0270*IT_0342 + IT_0267*conj(IT_0342) + conj
      (IT_0304)*IT_0455 + IT_0304*IT_0458 + conj(IT_0138)*IT_0461 + IT_0138
      *IT_0464) + (conj(IT_0332)*IT_0455 + IT_0332*IT_0458 + conj(IT_0342)
      *IT_0461 + IT_0342*IT_0464)*IT_0470 + IT_0476*(conj(IT_0354)*IT_0479 +
       IT_0354*IT_0482 + IT_0362*IT_0485 + IT_0363*IT_0488 + conj(IT_0355)
      *IT_0491 + IT_0355*IT_0494 + conj(IT_0366)*IT_0497 + IT_0366*IT_0500) +
       IT_0504*(conj(IT_0356)*IT_0479 + IT_0356*IT_0482 + conj(IT_0353)*IT_0491 
      + IT_0353*IT_0494 + conj(IT_0359)*IT_0497 + IT_0359*IT_0500 + IT_0485
      *IT_0505 + IT_0488*IT_0506) + (conj(IT_0366)*IT_0479 + IT_0366*IT_0482 +
       conj(IT_0355)*IT_0485 + IT_0355*IT_0488 + IT_0362*IT_0491 + IT_0363
      *IT_0494 + conj(IT_0354)*IT_0497 + IT_0354*IT_0500)*IT_0512 + (conj
      (IT_0359)*IT_0479 + IT_0359*IT_0482 + conj(IT_0353)*IT_0485 + IT_0353
      *IT_0488 + conj(IT_0356)*IT_0497 + IT_0356*IT_0500 + IT_0491*IT_0505 +
       IT_0494*IT_0506)*IT_0516 + 4*IT_0523*(IT_0438*IT_0439 + IT_0450*(IT_0449 
      + 0.25*IT_0524)) + 6*IT_0000*(IT_0439*IT_0450 + IT_0451*(conj(IT_0437) +
       0.166666666666667*IT_0525)) + (-8)*IT_0530*IT_0543 + (-4)*IT_0545*IT_0551
       + IT_0553*(IT_0450*IT_0455 + IT_0451*IT_0458 + conj(IT_0342)*IT_0477 +
       IT_0342*IT_0480 + conj(IT_0332)*IT_0489 + IT_0332*IT_0492 + conj(IT_0437)
      *IT_0556 + IT_0437*IT_0559) + (conj(IT_0198)*IT_0199 + conj(IT_0265)
      *IT_0307 + IT_0265*IT_0309 + conj(IT_0304)*IT_0333 + conj(IT_0138)*IT_0345
       + IT_0138*IT_0347)*IT_0560 + IT_0562*(conj(IT_0329)*IT_0453 + conj
      (IT_0321)*IT_0461 + IT_0321*IT_0463 + conj(IT_0342)*IT_0554 + conj(IT_0332
      )*IT_0565 + IT_0332*IT_0566) + IT_0570*IT_0572 + IT_0573*IT_0574 + IT_0576
      *IT_0597 + IT_0598*IT_0599;
    return create_ccomplex_return(IT_0600);
}

