#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_4_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_4_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
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
    const ccomplex_t IT_0000 = m_N_3*m_N_4;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_3, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_H0, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0008 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(alpha);
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = IT_0017*IT_0022;
    const ccomplex_t IT_0027 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0016 + -IT_0020 + 
      -IT_0025 + IT_0028);
    const ccomplex_t IT_0030 = IT_0013*IT_0021;
    const ccomplex_t IT_0031 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = IT_0012*IT_0022;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = IT_0033*IT_0035;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0032 + IT_0034 + 
      -IT_0036 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0029*IT_0039;
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = (-2)*IT_0041;
    const ccomplex_t IT_0043 = pow(m_N_4, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0008 + IT_0043 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0049 = IT_0026*IT_0048;
    const ccomplex_t IT_0050 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0051 = IT_0014*IT_0050;
    const ccomplex_t IT_0052 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0053 = IT_0026*IT_0052;
    const ccomplex_t IT_0054 = IT_0047 + IT_0049 + IT_0051 + IT_0053;
    const ccomplex_t IT_0055 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0056 = IT_0018*IT_0055;
    const ccomplex_t IT_0057 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0058 = IT_0023*IT_0057;
    const ccomplex_t IT_0059 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0060 = IT_0018*IT_0059;
    const ccomplex_t IT_0061 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0062 = IT_0023*IT_0061;
    const ccomplex_t IT_0063 = -IT_0056 + -IT_0058 + -IT_0060 + -IT_0062;
    const ccomplex_t IT_0064 = IT_0054 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0068 = IT_0030*IT_0067;
    const ccomplex_t IT_0069 = IT_0033*IT_0067;
    const ccomplex_t IT_0070 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0071 = IT_0030*IT_0070;
    const ccomplex_t IT_0072 = IT_0033*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0045*IT_0075;
    const ccomplex_t IT_0077 = cos(beta);
    const ccomplex_t IT_0078 = IT_0013*IT_0077;
    const ccomplex_t IT_0079 = IT_0055*IT_0078;
    const ccomplex_t IT_0080 = sin(beta);
    const ccomplex_t IT_0081 = IT_0013*IT_0080;
    const ccomplex_t IT_0082 = IT_0046*IT_0081;
    const ccomplex_t IT_0083 = IT_0022*IT_0077;
    const ccomplex_t IT_0084 = IT_0048*IT_0083;
    const ccomplex_t IT_0085 = IT_0022*IT_0080;
    const ccomplex_t IT_0086 = IT_0057*IT_0085;
    const ccomplex_t IT_0087 = IT_0059*IT_0078;
    const ccomplex_t IT_0088 = IT_0050*IT_0081;
    const ccomplex_t IT_0089 = IT_0052*IT_0083;
    const ccomplex_t IT_0090 = IT_0061*IT_0085;
    const ccomplex_t IT_0091 = (-0.5)*IT_0079 + 0.5*IT_0082 + 0.5*IT_0084 + (
      -0.5)*IT_0086 + (-0.5)*IT_0087 + 0.5*IT_0088 + 0.5*IT_0089 + (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = -e_em*IT_0013*IT_0022*(IT_0017*IT_0077 + 
      -IT_0011*IT_0080);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_12 + -IT_0004 + -IT_0043 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0073;
    const ccomplex_t IT_0098 = IT_0066*IT_0097;
    const ccomplex_t IT_0099 = IT_0045*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_23 + IT_0008 + IT_0043 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0103 = IT_0014*IT_0102;
    const ccomplex_t IT_0104 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0105 = IT_0026*IT_0104;
    const ccomplex_t IT_0106 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0107 = IT_0014*IT_0106;
    const ccomplex_t IT_0108 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0109 = IT_0026*IT_0108;
    const ccomplex_t IT_0110 = IT_0103 + IT_0105 + IT_0107 + IT_0109;
    const ccomplex_t IT_0111 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0112 = IT_0018*IT_0111;
    const ccomplex_t IT_0113 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0114 = IT_0023*IT_0113;
    const ccomplex_t IT_0115 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0116 = IT_0018*IT_0115;
    const ccomplex_t IT_0117 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0118 = IT_0023*IT_0117;
    const ccomplex_t IT_0119 = -IT_0112 + -IT_0114 + -IT_0116 + -IT_0118;
    const ccomplex_t IT_0120 = IT_0110 + IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0124 = IT_0030*IT_0123;
    const ccomplex_t IT_0125 = IT_0033*IT_0123;
    const ccomplex_t IT_0126 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0127 = IT_0030*IT_0126;
    const ccomplex_t IT_0128 = IT_0033*IT_0126;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0124 + IT_0125 + 
      -IT_0127 + -IT_0128);
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = IT_0122*IT_0130;
    const ccomplex_t IT_0132 = IT_0101*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_23 + IT_0008 + IT_0043 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0136 = IT_0014*IT_0135;
    const ccomplex_t IT_0137 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0138 = IT_0026*IT_0137;
    const ccomplex_t IT_0139 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0140 = IT_0014*IT_0139;
    const ccomplex_t IT_0141 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0142 = IT_0026*IT_0141;
    const ccomplex_t IT_0143 = IT_0136 + IT_0138 + IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0145 = IT_0018*IT_0144;
    const ccomplex_t IT_0146 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0147 = IT_0023*IT_0146;
    const ccomplex_t IT_0148 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0149 = IT_0018*IT_0148;
    const ccomplex_t IT_0150 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0151 = IT_0023*IT_0150;
    const ccomplex_t IT_0152 = -IT_0145 + -IT_0147 + -IT_0149 + -IT_0151;
    const ccomplex_t IT_0153 = IT_0143 + IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = 0.5*IT_0154;
    const ccomplex_t IT_0156 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0157 = IT_0033*IT_0156;
    const ccomplex_t IT_0158 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0159 = IT_0030*IT_0158;
    const ccomplex_t IT_0160 = IT_0033*IT_0158;
    const ccomplex_t IT_0161 = IT_0030*IT_0156;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0157 + -IT_0159 + 
      -IT_0160 + IT_0161);
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0155*IT_0163;
    const ccomplex_t IT_0165 = IT_0134*IT_0164;
    const ccomplex_t IT_0166 = IT_0046*IT_0078;
    const ccomplex_t IT_0167 = IT_0055*IT_0081;
    const ccomplex_t IT_0168 = IT_0057*IT_0083;
    const ccomplex_t IT_0169 = IT_0048*IT_0085;
    const ccomplex_t IT_0170 = IT_0050*IT_0078;
    const ccomplex_t IT_0171 = IT_0059*IT_0081;
    const ccomplex_t IT_0172 = IT_0061*IT_0083;
    const ccomplex_t IT_0173 = IT_0052*IT_0085;
    const ccomplex_t IT_0174 = (-0.5)*IT_0166 + (-0.5)*IT_0167 + 0.5*IT_0168 +
       0.5*IT_0169 + (-0.5)*IT_0170 + (-0.5)*IT_0171 + 0.5*IT_0172 + 0.5*IT_0173;
    const ccomplex_t IT_0175 = e_em*IT_0013*IT_0022*(IT_0011*IT_0077 + IT_0017
      *IT_0080);
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0174*IT_0176;
    const ccomplex_t IT_0178 = pow(m_Z, 2);
    const ccomplex_t IT_0179 = cpow((-2)*s_12 + -IT_0004 + -IT_0043 + IT_0178 
      + -reg_prop, -1);
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0177*IT_0180;
    const ccomplex_t IT_0182 = -IT_0076 + IT_0096 + IT_0099 + (-2)*IT_0132 + (
      -2)*IT_0165 + IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0182;
    const ccomplex_t IT_0184 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0185 = IT_0014*IT_0184;
    const ccomplex_t IT_0186 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0187 = IT_0018*IT_0186;
    const ccomplex_t IT_0188 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0189 = IT_0023*IT_0188;
    const ccomplex_t IT_0190 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0191 = IT_0026*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0185 + -IT_0187 + 
      -IT_0189 + IT_0191);
    const ccomplex_t IT_0193 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0194 = IT_0030*IT_0193;
    const ccomplex_t IT_0195 = IT_0033*IT_0193;
    const ccomplex_t IT_0196 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0197 = IT_0030*IT_0196;
    const ccomplex_t IT_0198 = IT_0033*IT_0196;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0194 + IT_0195 + 
      -IT_0197 + -IT_0198);
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = IT_0192*IT_0200;
    const ccomplex_t IT_0202 = IT_0010*IT_0201;
    const ccomplex_t IT_0203 = (-2)*IT_0202;
    const ccomplex_t IT_0204 = (-2)*IT_0203;
    const ccomplex_t IT_0205 = IT_0183 + IT_0204;
    const ccomplex_t IT_0206 = 2*conj(IT_0182);
    const ccomplex_t IT_0207 = (-2)*conj(IT_0203);
    const ccomplex_t IT_0208 = IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0210 = IT_0078*IT_0209;
    const ccomplex_t IT_0211 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0212 = IT_0081*IT_0211;
    const ccomplex_t IT_0213 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0214 = IT_0083*IT_0213;
    const ccomplex_t IT_0215 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0216 = IT_0085*IT_0215;
    const ccomplex_t IT_0217 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0218 = IT_0078*IT_0217;
    const ccomplex_t IT_0219 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0220 = IT_0081*IT_0219;
    const ccomplex_t IT_0221 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0222 = IT_0083*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0224 = IT_0085*IT_0223;
    const ccomplex_t IT_0225 = 0.5*IT_0210 + (-0.5)*IT_0212 + (-0.5)*IT_0214 +
       0.5*IT_0216 + 0.5*IT_0218 + (-0.5)*IT_0220 + (-0.5)*IT_0222 + 0.5*IT_0224;
    const ccomplex_t IT_0226 = IT_0092*IT_0225;
    const ccomplex_t IT_0227 = IT_0095*IT_0226;
    const ccomplex_t IT_0228 = IT_0078*IT_0211;
    const ccomplex_t IT_0229 = IT_0081*IT_0209;
    const ccomplex_t IT_0230 = IT_0083*IT_0215;
    const ccomplex_t IT_0231 = IT_0085*IT_0213;
    const ccomplex_t IT_0232 = IT_0078*IT_0219;
    const ccomplex_t IT_0233 = IT_0081*IT_0217;
    const ccomplex_t IT_0234 = IT_0083*IT_0223;
    const ccomplex_t IT_0235 = IT_0085*IT_0221;
    const ccomplex_t IT_0236 = 0.5*IT_0228 + 0.5*IT_0229 + (-0.5)*IT_0230 + (
      -0.5)*IT_0231 + 0.5*IT_0232 + 0.5*IT_0233 + (-0.5)*IT_0234 + (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = IT_0176*IT_0236;
    const ccomplex_t IT_0238 = IT_0180*IT_0237;
    const ccomplex_t IT_0239 = IT_0014*IT_0211;
    const ccomplex_t IT_0240 = IT_0026*IT_0213;
    const ccomplex_t IT_0241 = IT_0014*IT_0219;
    const ccomplex_t IT_0242 = IT_0026*IT_0221;
    const ccomplex_t IT_0243 = IT_0239 + IT_0240 + IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = IT_0018*IT_0209;
    const ccomplex_t IT_0245 = IT_0023*IT_0215;
    const ccomplex_t IT_0246 = IT_0018*IT_0217;
    const ccomplex_t IT_0247 = IT_0023*IT_0223;
    const ccomplex_t IT_0248 = -IT_0244 + -IT_0245 + -IT_0246 + -IT_0247;
    const ccomplex_t IT_0249 = IT_0243 + IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = IT_0074*IT_0251;
    const ccomplex_t IT_0253 = IT_0045*IT_0252;
    const ccomplex_t IT_0254 = IT_0097*IT_0251;
    const ccomplex_t IT_0255 = IT_0045*IT_0254;
    const ccomplex_t IT_0256 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0257 = IT_0030*IT_0256;
    const ccomplex_t IT_0258 = IT_0033*IT_0256;
    const ccomplex_t IT_0259 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0260 = IT_0030*IT_0259;
    const ccomplex_t IT_0261 = IT_0033*IT_0259;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*(IT_0257 + IT_0258 + 
      -IT_0260 + -IT_0261);
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0265 = IT_0014*IT_0264;
    const ccomplex_t IT_0266 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0267 = IT_0026*IT_0266;
    const ccomplex_t IT_0268 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0269 = IT_0014*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0271 = IT_0026*IT_0270;
    const ccomplex_t IT_0272 = IT_0265 + IT_0267 + IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0274 = IT_0018*IT_0273;
    const ccomplex_t IT_0275 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0276 = IT_0023*IT_0275;
    const ccomplex_t IT_0277 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0278 = IT_0018*IT_0277;
    const ccomplex_t IT_0279 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0280 = IT_0023*IT_0279;
    const ccomplex_t IT_0281 = -IT_0274 + -IT_0276 + -IT_0278 + -IT_0280;
    const ccomplex_t IT_0282 = IT_0272 + IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = 0.5*IT_0283;
    const ccomplex_t IT_0285 = IT_0263*IT_0284;
    const ccomplex_t IT_0286 = IT_0101*IT_0285;
    const ccomplex_t IT_0287 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0288 = IT_0014*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0290 = IT_0026*IT_0289;
    const ccomplex_t IT_0291 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0292 = IT_0014*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0294 = IT_0026*IT_0293;
    const ccomplex_t IT_0295 = IT_0288 + IT_0290 + IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0297 = IT_0018*IT_0296;
    const ccomplex_t IT_0298 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0299 = IT_0023*IT_0298;
    const ccomplex_t IT_0300 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0301 = IT_0018*IT_0300;
    const ccomplex_t IT_0302 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0303 = IT_0023*IT_0302;
    const ccomplex_t IT_0304 = -IT_0297 + -IT_0299 + -IT_0301 + -IT_0303;
    const ccomplex_t IT_0305 = IT_0295 + IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = 0.5*IT_0306;
    const ccomplex_t IT_0308 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0309 = IT_0030*IT_0308;
    const ccomplex_t IT_0310 = IT_0033*IT_0308;
    const ccomplex_t IT_0311 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0312 = IT_0030*IT_0311;
    const ccomplex_t IT_0313 = IT_0033*IT_0311;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*(IT_0309 + IT_0310 + 
      -IT_0312 + -IT_0313);
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = (-0.5)*IT_0315;
    const ccomplex_t IT_0317 = IT_0307*IT_0316;
    const ccomplex_t IT_0318 = IT_0134*IT_0317;
    const ccomplex_t IT_0319 = IT_0227 + IT_0238 + IT_0253 + -IT_0255 + (-2)
      *IT_0286 + (-2)*IT_0318;
    const ccomplex_t IT_0320 = 2*IT_0203;
    const ccomplex_t IT_0321 = (-2)*IT_0182;
    const ccomplex_t IT_0322 = IT_0320 + IT_0321;
    const ccomplex_t IT_0323 = 2*conj(IT_0203);
    const ccomplex_t IT_0324 = (-2)*conj(IT_0182);
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = s_12 + IT_0327;
    const ccomplex_t IT_0329 = m_N_3*IT_0328;
    const ccomplex_t IT_0330 = m_N_3*IT_0098;
    const ccomplex_t IT_0331 = IT_0045*IT_0330;
    const ccomplex_t IT_0332 = m_N_3*IT_0075;
    const ccomplex_t IT_0333 = -IT_0332;
    const ccomplex_t IT_0334 = IT_0330 + IT_0333;
    const ccomplex_t IT_0335 = IT_0045*IT_0334;
    const ccomplex_t IT_0336 = m_N_3*IT_0131;
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = IT_0122*IT_0263;
    const ccomplex_t IT_0339 = m_N_1*IT_0338;
    const ccomplex_t IT_0340 = IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0101*IT_0340;
    const ccomplex_t IT_0342 = IT_0155*IT_0316;
    const ccomplex_t IT_0343 = m_N_2*IT_0342;
    const ccomplex_t IT_0344 = m_N_3*IT_0164;
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = IT_0343 + IT_0345;
    const ccomplex_t IT_0347 = IT_0134*IT_0346;
    const ccomplex_t IT_0348 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0351 = IT_0030*IT_0350;
    const ccomplex_t IT_0352 = IT_0033*IT_0350;
    const ccomplex_t IT_0353 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0354 = IT_0030*IT_0353;
    const ccomplex_t IT_0355 = IT_0033*IT_0353;
    const ccomplex_t IT_0356 = (0 + _Complex_I*1)*(IT_0351 + IT_0352 + 
      -IT_0354 + -IT_0355);
    const ccomplex_t IT_0357 = 0.5*IT_0356;
    const ccomplex_t IT_0358 = IT_0251*IT_0357;
    const ccomplex_t IT_0359 = m_N_4*IT_0358;
    const ccomplex_t IT_0360 = (-0.5)*IT_0356;
    const ccomplex_t IT_0361 = IT_0251*IT_0360;
    const ccomplex_t IT_0362 = m_N_4*IT_0361;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = IT_0359 + IT_0363;
    const ccomplex_t IT_0365 = IT_0349*IT_0364;
    const ccomplex_t IT_0366 = 0.5*IT_0331 + 0.5*IT_0335 + IT_0341 + IT_0347 +
       (-0.5)*IT_0365;
    const ccomplex_t IT_0367 = IT_0045*IT_0332;
    const ccomplex_t IT_0368 = cpow(IT_0012, -4);
    const ccomplex_t IT_0369 = IT_0011*IT_0021*IT_0077*IT_0368;
    const ccomplex_t IT_0370 = IT_0017*IT_0021*IT_0080*IT_0368;
    const ccomplex_t IT_0371 = cpow(IT_0021, -3);
    const ccomplex_t IT_0372 = IT_0011*IT_0077*IT_0371;
    const ccomplex_t IT_0373 = IT_0017*IT_0080*IT_0371;
    const ccomplex_t IT_0374 = cpow(IT_0012, -2);
    const ccomplex_t IT_0375 = IT_0011*IT_0022*IT_0077*IT_0374;
    const ccomplex_t IT_0376 = IT_0017*IT_0022*IT_0080*IT_0374;
    const ccomplex_t IT_0377 = cpow(IT_0021, -2);
    const ccomplex_t IT_0378 = cpow(IT_0374 + IT_0377, -1);
    const ccomplex_t IT_0379 = (IT_0369 + IT_0370 + IT_0372 + IT_0373 + 2
      *IT_0375 + 2*IT_0376)*IT_0378;
    const ccomplex_t IT_0380 = m_W*e_em;
    const ccomplex_t IT_0381 = IT_0379*IT_0380;
    const ccomplex_t IT_0382 = (0 + _Complex_I*1)*IT_0381;
    const ccomplex_t IT_0383 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0004 + -IT_0043 + -reg_prop, -1);
    const ccomplex_t IT_0384 = IT_0382*IT_0383;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*IT_0384;
    const ccomplex_t IT_0386 = IT_0039*IT_0385;
    const ccomplex_t IT_0387 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0388 = (0 + _Complex_I*1)*IT_0387;
    const ccomplex_t IT_0389 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0390 = IT_0014*IT_0389;
    const ccomplex_t IT_0391 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0392 = IT_0026*IT_0391;
    const ccomplex_t IT_0393 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0394 = IT_0014*IT_0393;
    const ccomplex_t IT_0395 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0396 = IT_0026*IT_0395;
    const ccomplex_t IT_0397 = IT_0390 + IT_0392 + IT_0394 + IT_0396;
    const ccomplex_t IT_0398 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0399 = IT_0018*IT_0398;
    const ccomplex_t IT_0400 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0401 = IT_0023*IT_0400;
    const ccomplex_t IT_0402 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0403 = IT_0018*IT_0402;
    const ccomplex_t IT_0404 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0405 = IT_0023*IT_0404;
    const ccomplex_t IT_0406 = -IT_0399 + -IT_0401 + -IT_0403 + -IT_0405;
    const ccomplex_t IT_0407 = IT_0397 + IT_0406;
    const ccomplex_t IT_0408 = (0 + _Complex_I*1)*IT_0407;
    const ccomplex_t IT_0409 = 0.5*IT_0408;
    const ccomplex_t IT_0410 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0411 = IT_0030*IT_0410;
    const ccomplex_t IT_0412 = IT_0033*IT_0410;
    const ccomplex_t IT_0413 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0414 = IT_0030*IT_0413;
    const ccomplex_t IT_0415 = IT_0033*IT_0413;
    const ccomplex_t IT_0416 = (0 + _Complex_I*1)*(IT_0411 + IT_0412 + 
      -IT_0414 + -IT_0415);
    const ccomplex_t IT_0417 = (-0.5)*IT_0416;
    const ccomplex_t IT_0418 = IT_0409*IT_0417;
    const ccomplex_t IT_0419 = m_N_1*IT_0418;
    const ccomplex_t IT_0420 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0421 = IT_0030*IT_0420;
    const ccomplex_t IT_0422 = IT_0033*IT_0420;
    const ccomplex_t IT_0423 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0424 = IT_0030*IT_0423;
    const ccomplex_t IT_0425 = IT_0033*IT_0423;
    const ccomplex_t IT_0426 = (0 + _Complex_I*1)*(IT_0421 + IT_0422 + 
      -IT_0424 + -IT_0425);
    const ccomplex_t IT_0427 = (-0.5)*IT_0426;
    const ccomplex_t IT_0428 = IT_0409*IT_0427;
    const ccomplex_t IT_0429 = m_N_4*IT_0428;
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = IT_0419 + IT_0430;
    const ccomplex_t IT_0432 = IT_0388*IT_0431;
    const ccomplex_t IT_0433 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0434 = (0 + _Complex_I*1)*IT_0433;
    const ccomplex_t IT_0435 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0436 = IT_0014*IT_0435;
    const ccomplex_t IT_0437 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0438 = IT_0026*IT_0437;
    const ccomplex_t IT_0439 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0440 = IT_0014*IT_0439;
    const ccomplex_t IT_0441 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0442 = IT_0026*IT_0441;
    const ccomplex_t IT_0443 = IT_0436 + IT_0438 + IT_0440 + IT_0442;
    const ccomplex_t IT_0444 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0445 = IT_0018*IT_0444;
    const ccomplex_t IT_0446 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0447 = IT_0023*IT_0446;
    const ccomplex_t IT_0448 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0449 = IT_0018*IT_0448;
    const ccomplex_t IT_0450 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0451 = IT_0023*IT_0450;
    const ccomplex_t IT_0452 = -IT_0445 + -IT_0447 + -IT_0449 + -IT_0451;
    const ccomplex_t IT_0453 = IT_0443 + IT_0452;
    const ccomplex_t IT_0454 = (0 + _Complex_I*1)*IT_0453;
    const ccomplex_t IT_0455 = 0.5*IT_0454;
    const ccomplex_t IT_0456 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0457 = IT_0030*IT_0456;
    const ccomplex_t IT_0458 = IT_0033*IT_0456;
    const ccomplex_t IT_0459 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0460 = IT_0033*IT_0459;
    const ccomplex_t IT_0461 = IT_0030*IT_0459;
    const ccomplex_t IT_0462 = (0 + _Complex_I*1)*(IT_0457 + IT_0458 + 
      -IT_0460 + -IT_0461);
    const ccomplex_t IT_0463 = -IT_0462;
    const ccomplex_t IT_0464 = (-0.5)*IT_0463;
    const ccomplex_t IT_0465 = IT_0455*IT_0464;
    const ccomplex_t IT_0466 = m_N_4*IT_0465;
    const ccomplex_t IT_0467 = -IT_0466;
    const ccomplex_t IT_0468 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0469 = IT_0030*IT_0468;
    const ccomplex_t IT_0470 = IT_0033*IT_0468;
    const ccomplex_t IT_0471 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0472 = IT_0033*IT_0471;
    const ccomplex_t IT_0473 = IT_0030*IT_0471;
    const ccomplex_t IT_0474 = (0 + _Complex_I*1)*(IT_0469 + IT_0470 + 
      -IT_0472 + -IT_0473);
    const ccomplex_t IT_0475 = -IT_0474;
    const ccomplex_t IT_0476 = -IT_0475;
    const ccomplex_t IT_0477 = -IT_0476;
    const ccomplex_t IT_0478 = (-0.5)*IT_0477;
    const ccomplex_t IT_0479 = IT_0455*IT_0478;
    const ccomplex_t IT_0480 = m_N_2*IT_0479;
    const ccomplex_t IT_0481 = IT_0467 + IT_0480;
    const ccomplex_t IT_0482 = IT_0434*IT_0481;
    const ccomplex_t IT_0483 = cpow((-2)*s_13 + IT_0008 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0484 = (0 + _Complex_I*1)*IT_0483;
    const ccomplex_t IT_0485 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0486 = IT_0014*IT_0485;
    const ccomplex_t IT_0487 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0488 = IT_0018*IT_0487;
    const ccomplex_t IT_0489 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0490 = IT_0023*IT_0489;
    const ccomplex_t IT_0491 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0492 = IT_0026*IT_0491;
    const ccomplex_t IT_0493 = (0 + _Complex_I*1)*(IT_0486 + -IT_0488 + 
      -IT_0490 + IT_0492);
    const ccomplex_t IT_0494 = IT_0200*IT_0493;
    const ccomplex_t IT_0495 = m_N_4*IT_0494;
    const ccomplex_t IT_0496 = -IT_0495;
    const ccomplex_t IT_0497 = IT_0039*IT_0493;
    const ccomplex_t IT_0498 = m_N_3*IT_0497;
    const ccomplex_t IT_0499 = IT_0496 + IT_0498;
    const ccomplex_t IT_0500 = IT_0484*IT_0499;
    const ccomplex_t IT_0501 = IT_0039*IT_0192;
    const ccomplex_t IT_0502 = m_N_4*IT_0501;
    const ccomplex_t IT_0503 = m_N_3*IT_0201;
    const ccomplex_t IT_0504 = -IT_0503;
    const ccomplex_t IT_0505 = IT_0502 + IT_0504;
    const ccomplex_t IT_0506 = IT_0010*IT_0505;
    const ccomplex_t IT_0507 = (-0.5)*IT_0365 + (-0.5)*IT_0367 + -IT_0386 + 
      -IT_0432 + -IT_0482 + -IT_0500 + -IT_0506;
    const ccomplex_t IT_0508 = conj(IT_0366) + conj(IT_0507);
    const ccomplex_t IT_0509 = IT_0366 + IT_0507;
    const ccomplex_t IT_0510 = IT_0130*IT_0284;
    const ccomplex_t IT_0511 = m_N_1*IT_0510;
    const ccomplex_t IT_0512 = m_N_3*IT_0285;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = IT_0511 + IT_0513;
    const ccomplex_t IT_0515 = IT_0101*IT_0514;
    const ccomplex_t IT_0516 = m_N_3*IT_0252;
    const ccomplex_t IT_0517 = IT_0045*IT_0516;
    const ccomplex_t IT_0518 = m_N_3*IT_0254;
    const ccomplex_t IT_0519 = -IT_0518;
    const ccomplex_t IT_0520 = IT_0516 + IT_0519;
    const ccomplex_t IT_0521 = IT_0045*IT_0520;
    const ccomplex_t IT_0522 = IT_0163*IT_0307;
    const ccomplex_t IT_0523 = m_N_2*IT_0522;
    const ccomplex_t IT_0524 = m_N_3*IT_0317;
    const ccomplex_t IT_0525 = -IT_0524;
    const ccomplex_t IT_0526 = IT_0523 + IT_0525;
    const ccomplex_t IT_0527 = IT_0134*IT_0526;
    const ccomplex_t IT_0528 = IT_0066*IT_0360;
    const ccomplex_t IT_0529 = m_N_4*IT_0528;
    const ccomplex_t IT_0530 = IT_0066*IT_0357;
    const ccomplex_t IT_0531 = m_N_4*IT_0530;
    const ccomplex_t IT_0532 = -IT_0531;
    const ccomplex_t IT_0533 = IT_0529 + IT_0532;
    const ccomplex_t IT_0534 = IT_0349*IT_0533;
    const ccomplex_t IT_0535 = IT_0515 + 0.5*IT_0517 + 0.5*IT_0521 + IT_0527 +
       (-0.5)*IT_0534;
    const ccomplex_t IT_0536 = IT_0200*IT_0385;
    const ccomplex_t IT_0537 = IT_0045*IT_0518;
    const ccomplex_t IT_0538 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0539 = IT_0014*IT_0538;
    const ccomplex_t IT_0540 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0541 = IT_0026*IT_0540;
    const ccomplex_t IT_0542 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0543 = IT_0014*IT_0542;
    const ccomplex_t IT_0544 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0545 = IT_0026*IT_0544;
    const ccomplex_t IT_0546 = IT_0539 + IT_0541 + IT_0543 + IT_0545;
    const ccomplex_t IT_0547 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0548 = IT_0018*IT_0547;
    const ccomplex_t IT_0549 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0550 = IT_0023*IT_0549;
    const ccomplex_t IT_0551 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0552 = IT_0018*IT_0551;
    const ccomplex_t IT_0553 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0554 = IT_0023*IT_0553;
    const ccomplex_t IT_0555 = -IT_0548 + -IT_0550 + -IT_0552 + -IT_0554;
    const ccomplex_t IT_0556 = IT_0546 + IT_0555;
    const ccomplex_t IT_0557 = (0 + _Complex_I*1)*IT_0556;
    const ccomplex_t IT_0558 = 0.5*IT_0557;
    const ccomplex_t IT_0559 = IT_0427*IT_0558;
    const ccomplex_t IT_0560 = m_N_1*IT_0559;
    const ccomplex_t IT_0561 = IT_0417*IT_0558;
    const ccomplex_t IT_0562 = m_N_4*IT_0561;
    const ccomplex_t IT_0563 = -IT_0562;
    const ccomplex_t IT_0564 = IT_0560 + IT_0563;
    const ccomplex_t IT_0565 = IT_0388*IT_0564;
    const ccomplex_t IT_0566 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0567 = IT_0014*IT_0566;
    const ccomplex_t IT_0568 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0569 = IT_0026*IT_0568;
    const ccomplex_t IT_0570 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0571 = IT_0014*IT_0570;
    const ccomplex_t IT_0572 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0573 = IT_0026*IT_0572;
    const ccomplex_t IT_0574 = IT_0567 + IT_0569 + IT_0571 + IT_0573;
    const ccomplex_t IT_0575 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0576 = IT_0018*IT_0575;
    const ccomplex_t IT_0577 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0578 = IT_0023*IT_0577;
    const ccomplex_t IT_0579 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0580 = IT_0018*IT_0579;
    const ccomplex_t IT_0581 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0582 = IT_0023*IT_0581;
    const ccomplex_t IT_0583 = -IT_0576 + -IT_0578 + -IT_0580 + -IT_0582;
    const ccomplex_t IT_0584 = IT_0574 + IT_0583;
    const ccomplex_t IT_0585 = (0 + _Complex_I*1)*IT_0584;
    const ccomplex_t IT_0586 = 0.5*IT_0585;
    const ccomplex_t IT_0587 = IT_0464*IT_0586;
    const ccomplex_t IT_0588 = m_N_2*IT_0587;
    const ccomplex_t IT_0589 = IT_0478*IT_0586;
    const ccomplex_t IT_0590 = m_N_4*IT_0589;
    const ccomplex_t IT_0591 = -IT_0590;
    const ccomplex_t IT_0592 = IT_0588 + IT_0591;
    const ccomplex_t IT_0593 = IT_0434*IT_0592;
    const ccomplex_t IT_0594 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0595 = IT_0014*IT_0594;
    const ccomplex_t IT_0596 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0597 = IT_0018*IT_0596;
    const ccomplex_t IT_0598 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0599 = IT_0023*IT_0598;
    const ccomplex_t IT_0600 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0601 = IT_0026*IT_0600;
    const ccomplex_t IT_0602 = (0 + _Complex_I*1)*(IT_0595 + -IT_0597 + 
      -IT_0599 + IT_0601);
    const ccomplex_t IT_0603 = IT_0200*IT_0602;
    const ccomplex_t IT_0604 = m_N_3*IT_0603;
    const ccomplex_t IT_0605 = IT_0039*IT_0602;
    const ccomplex_t IT_0606 = m_N_4*IT_0605;
    const ccomplex_t IT_0607 = -IT_0606;
    const ccomplex_t IT_0608 = IT_0604 + IT_0607;
    const ccomplex_t IT_0609 = IT_0484*IT_0608;
    const ccomplex_t IT_0610 = IT_0029*IT_0200;
    const ccomplex_t IT_0611 = m_N_4*IT_0610;
    const ccomplex_t IT_0612 = m_N_3*IT_0040;
    const ccomplex_t IT_0613 = -IT_0612;
    const ccomplex_t IT_0614 = IT_0611 + IT_0613;
    const ccomplex_t IT_0615 = IT_0010*IT_0614;
    const ccomplex_t IT_0616 = (-0.5)*IT_0534 + -IT_0536 + (-0.5)*IT_0537 + 
      -IT_0565 + -IT_0593 + -IT_0609 + -IT_0615;
    const ccomplex_t IT_0617 = conj(IT_0535) + conj(IT_0616);
    const ccomplex_t IT_0618 = 2*IT_0319;
    const ccomplex_t IT_0619 = (-2)*IT_0042;
    const ccomplex_t IT_0620 = IT_0618 + IT_0619;
    const ccomplex_t IT_0621 = IT_0535 + IT_0616;
    const ccomplex_t IT_0622 = 2*conj(IT_0319);
    const ccomplex_t IT_0623 = (-2)*conj(IT_0042);
    const ccomplex_t IT_0624 = IT_0622 + IT_0623;
    const ccomplex_t IT_0625 = s_12*IT_0328;
    const ccomplex_t IT_0626 = IT_0484*IT_0605;
    const ccomplex_t IT_0627 = IT_0434*IT_0589;
    const ccomplex_t IT_0628 = IT_0349*IT_0530;
    const ccomplex_t IT_0629 = IT_0388*IT_0561;
    const ccomplex_t IT_0630 = 2*IT_0626 + 2*IT_0627 + IT_0628 + 2*IT_0629;
    const ccomplex_t IT_0631 = IT_0349*IT_0528;
    const ccomplex_t IT_0632 = IT_0096 + IT_0181 + IT_0631;
    const ccomplex_t IT_0633 = IT_0484*IT_0494;
    const ccomplex_t IT_0634 = IT_0349*IT_0361;
    const ccomplex_t IT_0635 = IT_0388*IT_0428;
    const ccomplex_t IT_0636 = IT_0434*IT_0465;
    const ccomplex_t IT_0637 = 2*IT_0633 + IT_0634 + 2*IT_0635 + 2*IT_0636;
    const ccomplex_t IT_0638 = IT_0349*IT_0358;
    const ccomplex_t IT_0639 = IT_0227 + IT_0238 + IT_0638;
    const ccomplex_t IT_0640 = 2*IT_0042;
    const ccomplex_t IT_0641 = (-2)*IT_0319;
    const ccomplex_t IT_0642 = IT_0640 + IT_0641;
    const ccomplex_t IT_0643 = 2*conj(IT_0042);
    const ccomplex_t IT_0644 = (-2)*conj(IT_0319);
    const ccomplex_t IT_0645 = IT_0643 + IT_0644;
    const ccomplex_t IT_0646 = s_12*s_14;
    const ccomplex_t IT_0647 = s_24*IT_0004;
    const ccomplex_t IT_0648 = -IT_0647;
    const ccomplex_t IT_0649 = IT_0646 + IT_0648;
    const ccomplex_t IT_0650 = (-0.5)*IT_0628;
    const ccomplex_t IT_0651 = 0.5*IT_0076 + (-0.5)*IT_0099 + IT_0132 + IT_0165;
    const ccomplex_t IT_0652 = -IT_0626 + -IT_0627 + -IT_0629;
    const ccomplex_t IT_0653 = conj(IT_0650) + conj(IT_0651) + conj(IT_0652);
    const ccomplex_t IT_0654 = (-0.5)*IT_0631;
    const ccomplex_t IT_0655 = conj(IT_0202) + conj(IT_0654);
    const ccomplex_t IT_0656 = IT_0202 + IT_0654;
    const ccomplex_t IT_0657 = (-0.5)*IT_0253 + 0.5*IT_0255 + IT_0286 + IT_0318;
    const ccomplex_t IT_0658 = -IT_0636;
    const ccomplex_t IT_0659 = -IT_0633 + (-0.5)*IT_0634 + -IT_0635;
    const ccomplex_t IT_0660 = conj(IT_0657) + conj(IT_0658) + conj(IT_0659);
    const ccomplex_t IT_0661 = IT_0657 + IT_0658 + IT_0659;
    const ccomplex_t IT_0662 = (-0.5)*IT_0638;
    const ccomplex_t IT_0663 = IT_0649*(IT_0208*(IT_0650 + IT_0651 + IT_0652) 
      + IT_0205*IT_0653 + IT_0322*IT_0655 + IT_0325*IT_0656 + IT_0620*IT_0660 +
       IT_0624*IT_0661 + IT_0645*(IT_0041 + IT_0662) + IT_0642*(conj(IT_0041) +
       conj(IT_0662)));
    const ccomplex_t IT_0664 = IT_0000*IT_0328;
    const ccomplex_t IT_0665 = 2*IT_0632;
    const ccomplex_t IT_0666 = (-2)*IT_0630;
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = 2*conj(IT_0632);
    const ccomplex_t IT_0669 = (-2)*conj(IT_0630);
    const ccomplex_t IT_0670 = IT_0668 + IT_0669;
    const ccomplex_t IT_0671 = 2*IT_0630;
    const ccomplex_t IT_0672 = (-2)*IT_0632;
    const ccomplex_t IT_0673 = IT_0671 + IT_0672;
    const ccomplex_t IT_0674 = 2*conj(IT_0630);
    const ccomplex_t IT_0675 = (-2)*conj(IT_0632);
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = pow(s_24, 2);
    const ccomplex_t IT_0678 = IT_0001*IT_0677;
    const ccomplex_t IT_0679 = -IT_0043;
    const ccomplex_t IT_0680 = IT_0678 + IT_0679;
    const ccomplex_t IT_0681 = IT_0000*IT_0680;
    const ccomplex_t IT_0682 = m_N_3*IT_0680;
    const ccomplex_t IT_0683 = 2*IT_0637;
    const ccomplex_t IT_0684 = (-2)*IT_0639;
    const ccomplex_t IT_0685 = IT_0683 + IT_0684;
    const ccomplex_t IT_0686 = 2*conj(IT_0637);
    const ccomplex_t IT_0687 = (-2)*conj(IT_0639);
    const ccomplex_t IT_0688 = IT_0686 + IT_0687;
    const ccomplex_t IT_0689 = s_24*IT_0178;
    const ccomplex_t IT_0690 = IT_0001*IT_0689;
    const ccomplex_t IT_0691 = 0.5*IT_0690;
    const ccomplex_t IT_0692 = s_24 + IT_0691;
    const ccomplex_t IT_0693 = m_N_3*IT_0692;
    const ccomplex_t IT_0694 = 2*IT_0693;
    const ccomplex_t IT_0695 = conj(IT_0650) + conj(IT_0652);
    const ccomplex_t IT_0696 = 2*IT_0366;
    const ccomplex_t IT_0697 = 2*IT_0507;
    const ccomplex_t IT_0698 = IT_0696 + IT_0697;
    const ccomplex_t IT_0699 = IT_0650 + IT_0652;
    const ccomplex_t IT_0700 = 2*conj(IT_0366);
    const ccomplex_t IT_0701 = 2*conj(IT_0507);
    const ccomplex_t IT_0702 = IT_0700 + IT_0701;
    const ccomplex_t IT_0703 = conj(IT_0658) + conj(IT_0659);
    const ccomplex_t IT_0704 = 2*IT_0535;
    const ccomplex_t IT_0705 = 2*IT_0616;
    const ccomplex_t IT_0706 = IT_0704 + IT_0705;
    const ccomplex_t IT_0707 = IT_0658 + IT_0659;
    const ccomplex_t IT_0708 = 2*conj(IT_0535);
    const ccomplex_t IT_0709 = 2*conj(IT_0616);
    const ccomplex_t IT_0710 = IT_0708 + IT_0709;
    const ccomplex_t IT_0711 = (-2)*IT_0366;
    const ccomplex_t IT_0712 = (-2)*IT_0507;
    const ccomplex_t IT_0713 = IT_0711 + IT_0712;
    const ccomplex_t IT_0714 = (-2)*conj(IT_0366);
    const ccomplex_t IT_0715 = (-2)*conj(IT_0507);
    const ccomplex_t IT_0716 = IT_0714 + IT_0715;
    const ccomplex_t IT_0717 = (-2)*IT_0535;
    const ccomplex_t IT_0718 = (-2)*IT_0616;
    const ccomplex_t IT_0719 = IT_0717 + IT_0718;
    const ccomplex_t IT_0720 = (-2)*conj(IT_0535);
    const ccomplex_t IT_0721 = (-2)*conj(IT_0616);
    const ccomplex_t IT_0722 = IT_0720 + IT_0721;
    const ccomplex_t IT_0723 = (-0.25)*IT_0690;
    const ccomplex_t IT_0724 = s_24 + IT_0723;
    const ccomplex_t IT_0725 = m_N_3*IT_0724;
    const ccomplex_t IT_0726 = (-4)*IT_0725;
    const ccomplex_t IT_0727 = s_14*IT_0178;
    const ccomplex_t IT_0728 = IT_0001*IT_0727;
    const ccomplex_t IT_0729 = 0.5*IT_0728;
    const ccomplex_t IT_0730 = s_14 + IT_0729;
    const ccomplex_t IT_0731 = m_N_4*IT_0730;
    const ccomplex_t IT_0732 = 2*IT_0731;
    const ccomplex_t IT_0733 = (-0.25)*IT_0728;
    const ccomplex_t IT_0734 = s_14 + IT_0733;
    const ccomplex_t IT_0735 = m_N_4*IT_0734;
    const ccomplex_t IT_0736 = (-4)*IT_0735;
    const ccomplex_t IT_0737 = s_12*IT_0178;
    const ccomplex_t IT_0738 = s_14*s_24;
    const ccomplex_t IT_0739 = (-2)*IT_0738;
    const ccomplex_t IT_0740 = IT_0737 + IT_0739;
    const ccomplex_t IT_0741 = IT_0001*IT_0740;
    const ccomplex_t IT_0742 = (-0.5)*IT_0741;
    const ccomplex_t IT_0743 = s_12 + IT_0742;
    const ccomplex_t IT_0744 = 4*IT_0535;
    const ccomplex_t IT_0745 = 6*conj(IT_0366);
    const ccomplex_t IT_0746 = pow(m_Z, 4);
    const ccomplex_t IT_0747 = s_12*IT_0746;
    const ccomplex_t IT_0748 = IT_0001*IT_0747;
    const ccomplex_t IT_0749 = (-0.25)*IT_0748;
    const ccomplex_t IT_0750 = IT_0738 + IT_0749;
    const ccomplex_t IT_0751 = 8*IT_0651;
    const ccomplex_t IT_0752 = 8*conj(IT_0651);
    const ccomplex_t IT_0753 = 8*IT_0657;
    const ccomplex_t IT_0754 = 8*conj(IT_0657);
    const ccomplex_t IT_0755 = conj(IT_0651)*IT_0654 + conj(IT_0202)*(IT_0650 
      + IT_0652 + -IT_0654) + IT_0651*conj(IT_0654) + IT_0202*(conj(IT_0650) +
       conj(IT_0652) + -conj(IT_0654)) + conj(IT_0657)*IT_0662 + conj(IT_0041)*
      (IT_0658 + IT_0659 + -IT_0662) + IT_0657*conj(IT_0662) + IT_0041*(conj
      (IT_0658) + conj(IT_0659) + -conj(IT_0662)) + (-0.125)*IT_0695*IT_0751 + (
      -0.125)*IT_0699*IT_0752 + (-0.125)*IT_0703*IT_0753 + (-0.125)*IT_0707
      *IT_0754;
    const ccomplex_t IT_0756 = (-0.5)*IT_0748;
    const ccomplex_t IT_0757 = IT_0737 + IT_0739 + IT_0756;
    const ccomplex_t IT_0758 = 4*conj(IT_0654);
    const ccomplex_t IT_0759 = 4*conj(IT_0662);
    const ccomplex_t IT_0760 = 4*IT_0654;
    const ccomplex_t IT_0761 = (-4)*IT_0652;
    const ccomplex_t IT_0762 = 4*IT_0662;
    const ccomplex_t IT_0763 = (-4)*IT_0658;
    const ccomplex_t IT_0764 = (IT_0202 + -IT_0651)*(conj(IT_0202) + -conj
      (IT_0651)) + IT_0654*conj(IT_0654) + (IT_0041 + -IT_0657)*(conj(IT_0041) +
       -conj(IT_0657)) + IT_0662*conj(IT_0662) + (-0.25)*IT_0652*IT_0758 +
       IT_0650*(IT_0695 + (-0.25)*IT_0758) + (-0.25)*IT_0707*IT_0759 + (-0.25)
      *IT_0695*(IT_0760 + IT_0761) + IT_0703*(IT_0659 + (-0.25)*IT_0762 + (-0.25
      )*IT_0763);
    const ccomplex_t IT_0765 = s_12*IT_0006;
    const ccomplex_t IT_0766 = s_12*IT_0680;
    const ccomplex_t IT_0767 = 2*IT_0639;
    const ccomplex_t IT_0768 = s_12*s_24;
    const ccomplex_t IT_0769 = s_14*IT_0043;
    const ccomplex_t IT_0770 = -IT_0769;
    const ccomplex_t IT_0771 = IT_0768 + IT_0770;
    const ccomplex_t IT_0772 = (-2)*IT_0637;
    const ccomplex_t IT_0773 = 2*conj(IT_0639);
    const ccomplex_t IT_0774 = (-2)*conj(IT_0637);
    const ccomplex_t IT_0775 = IT_0655*IT_0667 + IT_0656*IT_0670 + IT_0653
      *IT_0673 + (IT_0650 + IT_0651 + IT_0652)*IT_0676 + IT_0660*IT_0685 +
       IT_0661*IT_0688 + (conj(IT_0041) + conj(IT_0662))*(IT_0767 + IT_0772) + 
      (IT_0041 + IT_0662)*(IT_0773 + IT_0774);
    const ccomplex_t IT_0776 = m_N_4*IT_0006;
    const ccomplex_t IT_0777 = m_N_4*IT_0328;
    const ccomplex_t IT_0778 = IT_0000*IT_0178;
    const ccomplex_t IT_0779 = 3*IT_0778;
    const ccomplex_t IT_0780 = 2*IT_0654;
    const ccomplex_t IT_0781 = (-2)*IT_0650;
    const ccomplex_t IT_0782 = (-2)*IT_0652;
    const ccomplex_t IT_0783 = IT_0780 + IT_0781 + IT_0782;
    const ccomplex_t IT_0784 = 2*conj(IT_0654);
    const ccomplex_t IT_0785 = (-2)*conj(IT_0650);
    const ccomplex_t IT_0786 = (-2)*conj(IT_0652);
    const ccomplex_t IT_0787 = IT_0784 + IT_0785 + IT_0786;
    const ccomplex_t IT_0788 = 2*IT_0651;
    const ccomplex_t IT_0789 = IT_0203 + IT_0788;
    const ccomplex_t IT_0790 = 2*conj(IT_0651);
    const ccomplex_t IT_0791 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0792 = IT_0790 + IT_0791;
    const ccomplex_t IT_0793 = 2*IT_0662;
    const ccomplex_t IT_0794 = 2*conj(IT_0662);
    const ccomplex_t IT_0795 = 2*IT_0657;
    const ccomplex_t IT_0796 = 2*conj(IT_0657);
    const ccomplex_t IT_0797 = (-2)*IT_0654;
    const ccomplex_t IT_0798 = (-2)*conj(IT_0654);
    const ccomplex_t IT_0799 = (-2)*IT_0651;
    const ccomplex_t IT_0800 = (-2)*conj(IT_0651);
    const ccomplex_t IT_0801 = conj(IT_0041)*IT_0783 + IT_0041*IT_0787 +
       IT_0703*IT_0789 + IT_0707*IT_0792 + conj(IT_0202)*IT_0793 + IT_0202
      *IT_0794 + IT_0695*IT_0795 + IT_0699*IT_0796 + conj(IT_0657)*IT_0797 +
       IT_0657*IT_0798 + conj(IT_0662)*IT_0799 + IT_0662*IT_0800;
    const ccomplex_t IT_0802 = (-3)*IT_0778;
    const ccomplex_t IT_0803 = IT_0703*IT_0783 + IT_0707*IT_0787 + conj
      (IT_0041)*IT_0789 + IT_0041*IT_0792 + IT_0695*IT_0793 + IT_0699*IT_0794 +
       conj(IT_0202)*IT_0795 + IT_0202*IT_0796 + conj(IT_0662)*IT_0797 + IT_0662
      *IT_0798 + conj(IT_0657)*IT_0799 + IT_0657*IT_0800;
    const ccomplex_t IT_0804 = IT_0007*(conj(IT_0042)*IT_0205 + IT_0042
      *IT_0208 + conj(IT_0319)*IT_0322 + IT_0319*IT_0325) + IT_0329*(IT_0205
      *IT_0508 + IT_0208*IT_0509 + IT_0617*IT_0620 + IT_0621*IT_0624) + IT_0625*
      (IT_0208*IT_0630 + IT_0205*conj(IT_0630) + IT_0325*IT_0632 + IT_0322*conj
      (IT_0632) + IT_0624*IT_0637 + IT_0620*conj(IT_0637) + conj(IT_0639)
      *IT_0642 + IT_0639*IT_0645) + IT_0663 + IT_0664*(IT_0325*IT_0637 + IT_0322
      *conj(IT_0637) + IT_0208*IT_0639 + IT_0205*conj(IT_0639) + conj(IT_0319)
      *IT_0667 + IT_0319*IT_0670 + conj(IT_0042)*IT_0673 + IT_0042*IT_0676) + 
      (conj(IT_0637)*IT_0667 + IT_0637*IT_0670 + conj(IT_0639)*IT_0673 + IT_0639
      *IT_0676)*IT_0681 + IT_0682*(IT_0508*IT_0673 + IT_0509*IT_0676 + IT_0617
      *IT_0685 + IT_0621*IT_0688) + IT_0694*(IT_0695*IT_0698 + IT_0699*IT_0702 +
       IT_0703*IT_0706 + IT_0707*IT_0710 + conj(IT_0654)*IT_0713 + IT_0654
      *IT_0716 + conj(IT_0662)*IT_0719 + IT_0662*IT_0722) + (conj(IT_0202)
      *IT_0698 + IT_0202*IT_0702 + conj(IT_0041)*IT_0706 + IT_0041*IT_0710 +
       conj(IT_0651)*IT_0713 + IT_0651*IT_0716 + conj(IT_0657)*IT_0719 + IT_0657
      *IT_0722)*IT_0726 + (conj(IT_0657)*IT_0698 + IT_0657*IT_0702 + conj
      (IT_0651)*IT_0706 + IT_0651*IT_0710 + conj(IT_0041)*IT_0713 + IT_0041
      *IT_0716 + conj(IT_0202)*IT_0719 + IT_0202*IT_0722)*IT_0732 + (conj
      (IT_0662)*IT_0698 + IT_0662*IT_0702 + conj(IT_0654)*IT_0706 + IT_0654
      *IT_0710 + IT_0703*IT_0713 + IT_0707*IT_0716 + IT_0695*IT_0719 + IT_0699
      *IT_0722)*IT_0736 + 4*IT_0743*(IT_0508*IT_0509 + IT_0617*(IT_0616 + 0.25
      *IT_0744)) + 6*IT_0000*(IT_0509*IT_0617 + IT_0621*(conj(IT_0507) +
       0.166666666666667*IT_0745)) + (-8)*IT_0750*IT_0755 + (-4)*IT_0757*IT_0764
       + (conj(IT_0182)*IT_0183 + conj(IT_0203)*IT_0322 + IT_0203*IT_0324 + conj
      (IT_0319)*IT_0618 + conj(IT_0042)*IT_0642 + IT_0042*IT_0644)*IT_0765 +
       IT_0766*(conj(IT_0632)*IT_0665 + conj(IT_0630)*IT_0673 + IT_0630*IT_0675 
      + conj(IT_0637)*IT_0685 + IT_0637*IT_0687 + conj(IT_0639)*IT_0767) +
       IT_0771*IT_0775 + (IT_0205*IT_0617 + conj(IT_0507)*IT_0620 + IT_0208
      *IT_0621 + IT_0507*IT_0624 + conj(IT_0319)*IT_0696 + IT_0319*IT_0700 +
       conj(IT_0042)*IT_0711 + IT_0042*IT_0714)*IT_0776 + (IT_0617*IT_0673 +
       IT_0621*IT_0676 + conj(IT_0507)*IT_0685 + IT_0507*IT_0688 + conj(IT_0637)
      *IT_0696 + IT_0637*IT_0700 + conj(IT_0639)*IT_0711 + IT_0639*IT_0714)
      *IT_0777 + IT_0779*IT_0801 + IT_0802*IT_0803;
    return create_ccomplex_return(IT_0804);
}

