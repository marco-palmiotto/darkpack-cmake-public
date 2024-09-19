#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_4_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_4_to_h_Z(
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
    const ccomplex_t IT_0000 = m_N_2*m_N_4;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_h, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0008 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(alpha);
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = IT_0017*IT_0022;
    const ccomplex_t IT_0027 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0016 + IT_0020 + 
      -IT_0025 + -IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0013*IT_0021;
    const ccomplex_t IT_0032 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = IT_0012*IT_0022;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = IT_0031*IT_0036;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0033 + IT_0035 + 
      -IT_0037 + -IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0030*IT_0043;
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = pow(m_N_4, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0008 + IT_0046 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0050 = IT_0014*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0052 = IT_0018*IT_0051;
    const ccomplex_t IT_0053 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0054 = IT_0014*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0056 = IT_0018*IT_0055;
    const ccomplex_t IT_0057 = IT_0050 + IT_0052 + IT_0054 + IT_0056;
    const ccomplex_t IT_0058 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0059 = IT_0023*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0061 = IT_0026*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0063 = IT_0023*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0065 = IT_0026*IT_0064;
    const ccomplex_t IT_0066 = -IT_0059 + -IT_0061 + -IT_0063 + -IT_0065;
    const ccomplex_t IT_0067 = IT_0057 + IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0071 = IT_0034*IT_0070;
    const ccomplex_t IT_0072 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0073 = IT_0031*IT_0072;
    const ccomplex_t IT_0074 = IT_0034*IT_0072;
    const ccomplex_t IT_0075 = IT_0031*IT_0070;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0071 + -IT_0073 + 
      -IT_0074 + IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0069*IT_0077;
    const ccomplex_t IT_0079 = IT_0048*IT_0078;
    const ccomplex_t IT_0080 = (-2)*IT_0045 + (-2)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0008 + IT_0046 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0034*IT_0083;
    const ccomplex_t IT_0085 = IT_0031*IT_0083;
    const ccomplex_t IT_0086 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0031*IT_0086;
    const ccomplex_t IT_0088 = IT_0034*IT_0086;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0084 + IT_0085 + 
      -IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0092 = IT_0014*IT_0091;
    const ccomplex_t IT_0093 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0094 = IT_0018*IT_0093;
    const ccomplex_t IT_0095 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0096 = IT_0014*IT_0095;
    const ccomplex_t IT_0097 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0098 = IT_0018*IT_0097;
    const ccomplex_t IT_0099 = IT_0092 + IT_0094 + IT_0096 + IT_0098;
    const ccomplex_t IT_0100 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0101 = IT_0023*IT_0100;
    const ccomplex_t IT_0102 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0103 = IT_0026*IT_0102;
    const ccomplex_t IT_0104 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0105 = IT_0023*IT_0104;
    const ccomplex_t IT_0106 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0107 = IT_0026*IT_0106;
    const ccomplex_t IT_0108 = -IT_0101 + -IT_0103 + -IT_0105 + -IT_0107;
    const ccomplex_t IT_0109 = IT_0099 + IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0090*IT_0111;
    const ccomplex_t IT_0113 = IT_0082*IT_0112;
    const ccomplex_t IT_0114 = cos(beta);
    const ccomplex_t IT_0115 = IT_0013*IT_0114;
    const ccomplex_t IT_0116 = IT_0093*IT_0115;
    const ccomplex_t IT_0117 = sin(beta);
    const ccomplex_t IT_0118 = IT_0013*IT_0117;
    const ccomplex_t IT_0119 = IT_0091*IT_0118;
    const ccomplex_t IT_0120 = IT_0022*IT_0114;
    const ccomplex_t IT_0121 = IT_0102*IT_0120;
    const ccomplex_t IT_0122 = IT_0022*IT_0117;
    const ccomplex_t IT_0123 = IT_0100*IT_0122;
    const ccomplex_t IT_0124 = IT_0097*IT_0115;
    const ccomplex_t IT_0125 = IT_0095*IT_0118;
    const ccomplex_t IT_0126 = IT_0106*IT_0120;
    const ccomplex_t IT_0127 = IT_0104*IT_0122;
    const ccomplex_t IT_0128 = (-0.5)*IT_0116 + (-0.5)*IT_0119 + 0.5*IT_0121 +
       0.5*IT_0123 + (-0.5)*IT_0124 + (-0.5)*IT_0125 + 0.5*IT_0126 + 0.5*IT_0127;
    const ccomplex_t IT_0129 = -e_em*IT_0013*IT_0022*(IT_0017*IT_0114 + 
      -IT_0011*IT_0117);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = pow(m_Z, 2);
    const ccomplex_t IT_0133 = cpow((-2)*s_12 + -IT_0004 + -IT_0046 + IT_0132 
      + -reg_prop, -1);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0131*IT_0134;
    const ccomplex_t IT_0136 = IT_0091*IT_0115;
    const ccomplex_t IT_0137 = IT_0093*IT_0118;
    const ccomplex_t IT_0138 = IT_0100*IT_0120;
    const ccomplex_t IT_0139 = IT_0102*IT_0122;
    const ccomplex_t IT_0140 = IT_0095*IT_0115;
    const ccomplex_t IT_0141 = IT_0097*IT_0118;
    const ccomplex_t IT_0142 = IT_0104*IT_0120;
    const ccomplex_t IT_0143 = IT_0106*IT_0122;
    const ccomplex_t IT_0144 = (-0.5)*IT_0136 + 0.5*IT_0137 + 0.5*IT_0138 + (
      -0.5)*IT_0139 + (-0.5)*IT_0140 + 0.5*IT_0141 + 0.5*IT_0142 + (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = e_em*IT_0013*IT_0022*(IT_0011*IT_0114 + IT_0017
      *IT_0117);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_12 + -IT_0004 + -IT_0046 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147*IT_0149;
    const ccomplex_t IT_0151 = (-0.5)*IT_0089;
    const ccomplex_t IT_0152 = IT_0111*IT_0151;
    const ccomplex_t IT_0153 = IT_0082*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_23 + IT_0008 + IT_0046 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0157 = IT_0014*IT_0156;
    const ccomplex_t IT_0158 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0159 = IT_0018*IT_0158;
    const ccomplex_t IT_0160 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0161 = IT_0014*IT_0160;
    const ccomplex_t IT_0162 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0163 = IT_0018*IT_0162;
    const ccomplex_t IT_0164 = IT_0157 + IT_0159 + IT_0161 + IT_0163;
    const ccomplex_t IT_0165 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0166 = IT_0023*IT_0165;
    const ccomplex_t IT_0167 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0168 = IT_0026*IT_0167;
    const ccomplex_t IT_0169 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0170 = IT_0023*IT_0169;
    const ccomplex_t IT_0171 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0172 = IT_0026*IT_0171;
    const ccomplex_t IT_0173 = -IT_0166 + -IT_0168 + -IT_0170 + -IT_0172;
    const ccomplex_t IT_0174 = IT_0164 + IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0178 = IT_0031*IT_0177;
    const ccomplex_t IT_0179 = IT_0034*IT_0177;
    const ccomplex_t IT_0180 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0181 = IT_0031*IT_0180;
    const ccomplex_t IT_0182 = IT_0034*IT_0180;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0178 + IT_0179 + 
      -IT_0181 + -IT_0182);
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0176*IT_0184;
    const ccomplex_t IT_0186 = IT_0155*IT_0185;
    const ccomplex_t IT_0187 = -IT_0113 + IT_0135 + IT_0150 + IT_0153 + (-2)
      *IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0190 = IT_0014*IT_0189;
    const ccomplex_t IT_0191 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0192 = IT_0018*IT_0191;
    const ccomplex_t IT_0193 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0194 = IT_0014*IT_0193;
    const ccomplex_t IT_0195 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0196 = IT_0018*IT_0195;
    const ccomplex_t IT_0197 = IT_0190 + IT_0192 + IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0199 = IT_0023*IT_0198;
    const ccomplex_t IT_0200 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0201 = IT_0026*IT_0200;
    const ccomplex_t IT_0202 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0203 = IT_0023*IT_0202;
    const ccomplex_t IT_0204 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0205 = IT_0026*IT_0204;
    const ccomplex_t IT_0206 = -IT_0199 + -IT_0201 + -IT_0203 + -IT_0205;
    const ccomplex_t IT_0207 = IT_0197 + IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0211 = IT_0031*IT_0210;
    const ccomplex_t IT_0212 = IT_0034*IT_0210;
    const ccomplex_t IT_0213 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0214 = IT_0031*IT_0213;
    const ccomplex_t IT_0215 = IT_0034*IT_0213;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0211 + IT_0212 + 
      -IT_0214 + -IT_0215);
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = (-0.5)*IT_0217;
    const ccomplex_t IT_0219 = IT_0209*IT_0218;
    const ccomplex_t IT_0220 = IT_0048*IT_0219;
    const ccomplex_t IT_0221 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0222 = IT_0014*IT_0221;
    const ccomplex_t IT_0223 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0224 = IT_0018*IT_0223;
    const ccomplex_t IT_0225 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0226 = IT_0023*IT_0225;
    const ccomplex_t IT_0227 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0228 = IT_0026*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0222 + IT_0224 + 
      -IT_0226 + -IT_0228);
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0232 = IT_0031*IT_0231;
    const ccomplex_t IT_0233 = IT_0034*IT_0231;
    const ccomplex_t IT_0234 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0235 = IT_0034*IT_0234;
    const ccomplex_t IT_0236 = IT_0031*IT_0234;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0232 + IT_0233 + 
      -IT_0235 + -IT_0236);
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = IT_0230*IT_0239;
    const ccomplex_t IT_0241 = IT_0010*IT_0240;
    const ccomplex_t IT_0242 = (-2)*IT_0220 + (-2)*IT_0241;
    const ccomplex_t IT_0243 = (-2)*IT_0242;
    const ccomplex_t IT_0244 = IT_0188 + IT_0243;
    const ccomplex_t IT_0245 = 2*conj(IT_0187);
    const ccomplex_t IT_0246 = (-2)*conj(IT_0242);
    const ccomplex_t IT_0247 = IT_0245 + IT_0246;
    const ccomplex_t IT_0248 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0249 = IT_0014*IT_0248;
    const ccomplex_t IT_0250 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0251 = IT_0018*IT_0250;
    const ccomplex_t IT_0252 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0253 = IT_0014*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0255 = IT_0018*IT_0254;
    const ccomplex_t IT_0256 = IT_0249 + IT_0251 + IT_0253 + IT_0255;
    const ccomplex_t IT_0257 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0258 = IT_0023*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0260 = IT_0026*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0262 = IT_0023*IT_0261;
    const ccomplex_t IT_0263 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0264 = IT_0026*IT_0263;
    const ccomplex_t IT_0265 = -IT_0258 + -IT_0260 + -IT_0262 + -IT_0264;
    const ccomplex_t IT_0266 = IT_0256 + IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = (-0.5)*IT_0267;
    const ccomplex_t IT_0269 = IT_0090*IT_0268;
    const ccomplex_t IT_0270 = IT_0082*IT_0269;
    const ccomplex_t IT_0271 = IT_0151*IT_0268;
    const ccomplex_t IT_0272 = IT_0082*IT_0271;
    const ccomplex_t IT_0273 = IT_0115*IT_0250;
    const ccomplex_t IT_0274 = IT_0118*IT_0248;
    const ccomplex_t IT_0275 = IT_0120*IT_0259;
    const ccomplex_t IT_0276 = IT_0122*IT_0257;
    const ccomplex_t IT_0277 = IT_0115*IT_0254;
    const ccomplex_t IT_0278 = IT_0118*IT_0252;
    const ccomplex_t IT_0279 = IT_0120*IT_0263;
    const ccomplex_t IT_0280 = IT_0122*IT_0261;
    const ccomplex_t IT_0281 = 0.5*IT_0273 + 0.5*IT_0274 + (-0.5)*IT_0275 + (
      -0.5)*IT_0276 + 0.5*IT_0277 + 0.5*IT_0278 + (-0.5)*IT_0279 + (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = IT_0130*IT_0281;
    const ccomplex_t IT_0283 = IT_0134*IT_0282;
    const ccomplex_t IT_0284 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0285 = IT_0014*IT_0284;
    const ccomplex_t IT_0286 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0287 = IT_0018*IT_0286;
    const ccomplex_t IT_0288 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0289 = IT_0014*IT_0288;
    const ccomplex_t IT_0290 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0291 = IT_0018*IT_0290;
    const ccomplex_t IT_0292 = IT_0285 + IT_0287 + IT_0289 + IT_0291;
    const ccomplex_t IT_0293 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0294 = IT_0023*IT_0293;
    const ccomplex_t IT_0295 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0296 = IT_0026*IT_0295;
    const ccomplex_t IT_0297 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0298 = IT_0023*IT_0297;
    const ccomplex_t IT_0299 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0300 = IT_0026*IT_0299;
    const ccomplex_t IT_0301 = -IT_0294 + -IT_0296 + -IT_0298 + -IT_0300;
    const ccomplex_t IT_0302 = IT_0292 + IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0306 = IT_0034*IT_0305;
    const ccomplex_t IT_0307 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0308 = IT_0034*IT_0307;
    const ccomplex_t IT_0309 = IT_0031*IT_0305;
    const ccomplex_t IT_0310 = IT_0031*IT_0307;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*(IT_0306 + -IT_0308 +
       IT_0309 + -IT_0310);
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = (-0.5)*IT_0312;
    const ccomplex_t IT_0314 = IT_0304*IT_0313;
    const ccomplex_t IT_0315 = IT_0155*IT_0314;
    const ccomplex_t IT_0316 = IT_0115*IT_0248;
    const ccomplex_t IT_0317 = IT_0118*IT_0250;
    const ccomplex_t IT_0318 = IT_0120*IT_0257;
    const ccomplex_t IT_0319 = IT_0122*IT_0259;
    const ccomplex_t IT_0320 = IT_0115*IT_0252;
    const ccomplex_t IT_0321 = IT_0118*IT_0254;
    const ccomplex_t IT_0322 = IT_0120*IT_0261;
    const ccomplex_t IT_0323 = IT_0122*IT_0263;
    const ccomplex_t IT_0324 = 0.5*IT_0316 + (-0.5)*IT_0317 + (-0.5)*IT_0318 +
       0.5*IT_0319 + 0.5*IT_0320 + (-0.5)*IT_0321 + (-0.5)*IT_0322 + 0.5*IT_0323;
    const ccomplex_t IT_0325 = IT_0146*IT_0324;
    const ccomplex_t IT_0326 = IT_0149*IT_0325;
    const ccomplex_t IT_0327 = IT_0270 + -IT_0272 + IT_0283 + (-2)*IT_0315 +
       IT_0326;
    const ccomplex_t IT_0328 = 2*IT_0242;
    const ccomplex_t IT_0329 = (-2)*IT_0187;
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = 2*conj(IT_0242);
    const ccomplex_t IT_0332 = (-2)*conj(IT_0187);
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = s_12 + IT_0335;
    const ccomplex_t IT_0337 = m_N_2*IT_0336;
    const ccomplex_t IT_0338 = m_N_2*IT_0152;
    const ccomplex_t IT_0339 = IT_0082*IT_0338;
    const ccomplex_t IT_0340 = m_N_2*IT_0112;
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = IT_0338 + IT_0341;
    const ccomplex_t IT_0343 = IT_0082*IT_0342;
    const ccomplex_t IT_0344 = IT_0176*IT_0313;
    const ccomplex_t IT_0345 = m_N_1*IT_0344;
    const ccomplex_t IT_0346 = m_N_2*IT_0185;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = IT_0345 + IT_0347;
    const ccomplex_t IT_0349 = IT_0155*IT_0348;
    const ccomplex_t IT_0350 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*IT_0350;
    const ccomplex_t IT_0352 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0353 = IT_0031*IT_0352;
    const ccomplex_t IT_0354 = IT_0034*IT_0352;
    const ccomplex_t IT_0355 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0356 = IT_0031*IT_0355;
    const ccomplex_t IT_0357 = IT_0034*IT_0355;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*(IT_0353 + IT_0354 + 
      -IT_0356 + -IT_0357);
    const ccomplex_t IT_0359 = 0.5*IT_0358;
    const ccomplex_t IT_0360 = IT_0268*IT_0359;
    const ccomplex_t IT_0361 = m_N_4*IT_0360;
    const ccomplex_t IT_0362 = (-0.5)*IT_0358;
    const ccomplex_t IT_0363 = IT_0268*IT_0362;
    const ccomplex_t IT_0364 = m_N_4*IT_0363;
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = IT_0361 + IT_0365;
    const ccomplex_t IT_0367 = IT_0351*IT_0366;
    const ccomplex_t IT_0368 = 0.5*IT_0339 + 0.5*IT_0343 + IT_0349 + (-0.5)
      *IT_0367;
    const ccomplex_t IT_0369 = m_W*e_em;
    const ccomplex_t IT_0370 = cpow(IT_0012, -4);
    const ccomplex_t IT_0371 = IT_0011*IT_0021*IT_0117*IT_0370;
    const ccomplex_t IT_0372 = IT_0017*IT_0021*IT_0114*IT_0370;
    const ccomplex_t IT_0373 = cpow(IT_0021, -3);
    const ccomplex_t IT_0374 = IT_0011*IT_0117*IT_0373;
    const ccomplex_t IT_0375 = IT_0017*IT_0114*IT_0373;
    const ccomplex_t IT_0376 = cpow(IT_0012, -2);
    const ccomplex_t IT_0377 = IT_0011*IT_0022*IT_0117*IT_0376;
    const ccomplex_t IT_0378 = IT_0017*IT_0022*IT_0114*IT_0376;
    const ccomplex_t IT_0379 = cpow(IT_0021, -2);
    const ccomplex_t IT_0380 = cpow(IT_0376 + IT_0379, -1);
    const ccomplex_t IT_0381 = (IT_0371 + -IT_0372 + IT_0374 + -IT_0375 + 2
      *IT_0377 + (-2)*IT_0378)*IT_0380;
    const ccomplex_t IT_0382 = IT_0369*IT_0381;
    const ccomplex_t IT_0383 = (0 + _Complex_I*1)*IT_0382;
    const ccomplex_t IT_0384 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0004 + -IT_0046 + -reg_prop, -1);
    const ccomplex_t IT_0385 = IT_0383*IT_0384;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*IT_0385;
    const ccomplex_t IT_0387 = IT_0043*IT_0386;
    const ccomplex_t IT_0388 = IT_0082*IT_0340;
    const ccomplex_t IT_0389 = IT_0077*IT_0209;
    const ccomplex_t IT_0390 = m_N_3*IT_0389;
    const ccomplex_t IT_0391 = m_N_2*IT_0219;
    const ccomplex_t IT_0392 = -IT_0391;
    const ccomplex_t IT_0393 = IT_0390 + IT_0392;
    const ccomplex_t IT_0394 = IT_0048*IT_0393;
    const ccomplex_t IT_0395 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*IT_0395;
    const ccomplex_t IT_0397 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0398 = IT_0014*IT_0397;
    const ccomplex_t IT_0399 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0400 = IT_0018*IT_0399;
    const ccomplex_t IT_0401 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0402 = IT_0014*IT_0401;
    const ccomplex_t IT_0403 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0404 = IT_0018*IT_0403;
    const ccomplex_t IT_0405 = IT_0398 + IT_0400 + IT_0402 + IT_0404;
    const ccomplex_t IT_0406 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0407 = IT_0023*IT_0406;
    const ccomplex_t IT_0408 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0409 = IT_0026*IT_0408;
    const ccomplex_t IT_0410 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0411 = IT_0023*IT_0410;
    const ccomplex_t IT_0412 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0413 = IT_0026*IT_0412;
    const ccomplex_t IT_0414 = -IT_0407 + -IT_0409 + -IT_0411 + -IT_0413;
    const ccomplex_t IT_0415 = IT_0405 + IT_0414;
    const ccomplex_t IT_0416 = (0 + _Complex_I*1)*IT_0415;
    const ccomplex_t IT_0417 = (-0.5)*IT_0416;
    const ccomplex_t IT_0418 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0419 = IT_0031*IT_0418;
    const ccomplex_t IT_0420 = IT_0034*IT_0418;
    const ccomplex_t IT_0421 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0422 = IT_0031*IT_0421;
    const ccomplex_t IT_0423 = IT_0034*IT_0421;
    const ccomplex_t IT_0424 = (0 + _Complex_I*1)*(IT_0419 + IT_0420 + 
      -IT_0422 + -IT_0423);
    const ccomplex_t IT_0425 = (-0.5)*IT_0424;
    const ccomplex_t IT_0426 = IT_0417*IT_0425;
    const ccomplex_t IT_0427 = m_N_1*IT_0426;
    const ccomplex_t IT_0428 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0429 = IT_0031*IT_0428;
    const ccomplex_t IT_0430 = IT_0034*IT_0428;
    const ccomplex_t IT_0431 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0432 = IT_0031*IT_0431;
    const ccomplex_t IT_0433 = IT_0034*IT_0431;
    const ccomplex_t IT_0434 = (0 + _Complex_I*1)*(IT_0429 + IT_0430 + 
      -IT_0432 + -IT_0433);
    const ccomplex_t IT_0435 = (-0.5)*IT_0434;
    const ccomplex_t IT_0436 = IT_0417*IT_0435;
    const ccomplex_t IT_0437 = m_N_4*IT_0436;
    const ccomplex_t IT_0438 = -IT_0437;
    const ccomplex_t IT_0439 = IT_0427 + IT_0438;
    const ccomplex_t IT_0440 = IT_0396*IT_0439;
    const ccomplex_t IT_0441 = cpow((-2)*s_13 + IT_0008 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0442 = (0 + _Complex_I*1)*IT_0441;
    const ccomplex_t IT_0443 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0444 = IT_0014*IT_0443;
    const ccomplex_t IT_0445 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0446 = IT_0018*IT_0445;
    const ccomplex_t IT_0447 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0448 = IT_0023*IT_0447;
    const ccomplex_t IT_0449 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0450 = IT_0026*IT_0449;
    const ccomplex_t IT_0451 = (0 + _Complex_I*1)*(IT_0444 + IT_0446 + 
      -IT_0448 + -IT_0450);
    const ccomplex_t IT_0452 = -IT_0451;
    const ccomplex_t IT_0453 = IT_0043*IT_0452;
    const ccomplex_t IT_0454 = m_N_2*IT_0453;
    const ccomplex_t IT_0455 = IT_0239*IT_0452;
    const ccomplex_t IT_0456 = m_N_4*IT_0455;
    const ccomplex_t IT_0457 = -IT_0456;
    const ccomplex_t IT_0458 = IT_0454 + IT_0457;
    const ccomplex_t IT_0459 = IT_0442*IT_0458;
    const ccomplex_t IT_0460 = IT_0043*IT_0230;
    const ccomplex_t IT_0461 = m_N_4*IT_0460;
    const ccomplex_t IT_0462 = m_N_2*IT_0240;
    const ccomplex_t IT_0463 = -IT_0462;
    const ccomplex_t IT_0464 = IT_0461 + IT_0463;
    const ccomplex_t IT_0465 = IT_0010*IT_0464;
    const ccomplex_t IT_0466 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0467 = (0 + _Complex_I*1)*IT_0466;
    const ccomplex_t IT_0468 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0469 = IT_0014*IT_0468;
    const ccomplex_t IT_0470 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0471 = IT_0018*IT_0470;
    const ccomplex_t IT_0472 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0473 = IT_0014*IT_0472;
    const ccomplex_t IT_0474 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0475 = IT_0018*IT_0474;
    const ccomplex_t IT_0476 = IT_0469 + IT_0471 + IT_0473 + IT_0475;
    const ccomplex_t IT_0477 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0478 = IT_0023*IT_0477;
    const ccomplex_t IT_0479 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0480 = IT_0026*IT_0479;
    const ccomplex_t IT_0481 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0482 = IT_0023*IT_0481;
    const ccomplex_t IT_0483 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0484 = IT_0026*IT_0483;
    const ccomplex_t IT_0485 = -IT_0478 + -IT_0480 + -IT_0482 + -IT_0484;
    const ccomplex_t IT_0486 = IT_0476 + IT_0485;
    const ccomplex_t IT_0487 = (0 + _Complex_I*1)*IT_0486;
    const ccomplex_t IT_0488 = (-0.5)*IT_0487;
    const ccomplex_t IT_0489 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0490 = IT_0031*IT_0489;
    const ccomplex_t IT_0491 = IT_0034*IT_0489;
    const ccomplex_t IT_0492 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0493 = IT_0031*IT_0492;
    const ccomplex_t IT_0494 = IT_0034*IT_0492;
    const ccomplex_t IT_0495 = (0 + _Complex_I*1)*(IT_0490 + IT_0491 + 
      -IT_0493 + -IT_0494);
    const ccomplex_t IT_0496 = (-0.5)*IT_0495;
    const ccomplex_t IT_0497 = IT_0488*IT_0496;
    const ccomplex_t IT_0498 = m_N_3*IT_0497;
    const ccomplex_t IT_0499 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0500 = IT_0031*IT_0499;
    const ccomplex_t IT_0501 = IT_0034*IT_0499;
    const ccomplex_t IT_0502 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0503 = IT_0031*IT_0502;
    const ccomplex_t IT_0504 = IT_0034*IT_0502;
    const ccomplex_t IT_0505 = (0 + _Complex_I*1)*(IT_0500 + IT_0501 + 
      -IT_0503 + -IT_0504);
    const ccomplex_t IT_0506 = (-0.5)*IT_0505;
    const ccomplex_t IT_0507 = IT_0488*IT_0506;
    const ccomplex_t IT_0508 = m_N_4*IT_0507;
    const ccomplex_t IT_0509 = -IT_0508;
    const ccomplex_t IT_0510 = IT_0498 + IT_0509;
    const ccomplex_t IT_0511 = IT_0467*IT_0510;
    const ccomplex_t IT_0512 = (-0.5)*IT_0367 + -IT_0387 + (-0.5)*IT_0388 + 
      -IT_0394 + -IT_0440 + -IT_0459 + -IT_0465 + -IT_0511;
    const ccomplex_t IT_0513 = conj(IT_0368) + conj(IT_0512);
    const ccomplex_t IT_0514 = IT_0368 + IT_0512;
    const ccomplex_t IT_0515 = m_N_2*IT_0269;
    const ccomplex_t IT_0516 = IT_0082*IT_0515;
    const ccomplex_t IT_0517 = IT_0111*IT_0362;
    const ccomplex_t IT_0518 = m_N_4*IT_0517;
    const ccomplex_t IT_0519 = IT_0111*IT_0359;
    const ccomplex_t IT_0520 = m_N_4*IT_0519;
    const ccomplex_t IT_0521 = -IT_0520;
    const ccomplex_t IT_0522 = IT_0518 + IT_0521;
    const ccomplex_t IT_0523 = IT_0351*IT_0522;
    const ccomplex_t IT_0524 = m_N_2*IT_0271;
    const ccomplex_t IT_0525 = -IT_0524;
    const ccomplex_t IT_0526 = IT_0515 + IT_0525;
    const ccomplex_t IT_0527 = IT_0082*IT_0526;
    const ccomplex_t IT_0528 = IT_0184*IT_0304;
    const ccomplex_t IT_0529 = m_N_1*IT_0528;
    const ccomplex_t IT_0530 = m_N_2*IT_0314;
    const ccomplex_t IT_0531 = -IT_0530;
    const ccomplex_t IT_0532 = IT_0529 + IT_0531;
    const ccomplex_t IT_0533 = IT_0155*IT_0532;
    const ccomplex_t IT_0534 = 0.5*IT_0516 + (-0.5)*IT_0523 + 0.5*IT_0527 +
       IT_0533;
    const ccomplex_t IT_0535 = IT_0069*IT_0218;
    const ccomplex_t IT_0536 = m_N_3*IT_0535;
    const ccomplex_t IT_0537 = m_N_2*IT_0078;
    const ccomplex_t IT_0538 = -IT_0537;
    const ccomplex_t IT_0539 = IT_0536 + IT_0538;
    const ccomplex_t IT_0540 = IT_0048*IT_0539;
    const ccomplex_t IT_0541 = IT_0239*IT_0386;
    const ccomplex_t IT_0542 = IT_0082*IT_0524;
    const ccomplex_t IT_0543 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0544 = IT_0014*IT_0543;
    const ccomplex_t IT_0545 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0546 = IT_0018*IT_0545;
    const ccomplex_t IT_0547 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0548 = IT_0014*IT_0547;
    const ccomplex_t IT_0549 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0550 = IT_0018*IT_0549;
    const ccomplex_t IT_0551 = IT_0544 + IT_0546 + IT_0548 + IT_0550;
    const ccomplex_t IT_0552 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0553 = IT_0023*IT_0552;
    const ccomplex_t IT_0554 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0555 = IT_0026*IT_0554;
    const ccomplex_t IT_0556 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0557 = IT_0023*IT_0556;
    const ccomplex_t IT_0558 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0559 = IT_0026*IT_0558;
    const ccomplex_t IT_0560 = -IT_0553 + -IT_0555 + -IT_0557 + -IT_0559;
    const ccomplex_t IT_0561 = IT_0551 + IT_0560;
    const ccomplex_t IT_0562 = (0 + _Complex_I*1)*IT_0561;
    const ccomplex_t IT_0563 = (-0.5)*IT_0562;
    const ccomplex_t IT_0564 = IT_0435*IT_0563;
    const ccomplex_t IT_0565 = m_N_1*IT_0564;
    const ccomplex_t IT_0566 = IT_0425*IT_0563;
    const ccomplex_t IT_0567 = m_N_4*IT_0566;
    const ccomplex_t IT_0568 = -IT_0567;
    const ccomplex_t IT_0569 = IT_0565 + IT_0568;
    const ccomplex_t IT_0570 = IT_0396*IT_0569;
    const ccomplex_t IT_0571 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0572 = IT_0014*IT_0571;
    const ccomplex_t IT_0573 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0574 = IT_0018*IT_0573;
    const ccomplex_t IT_0575 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0576 = IT_0023*IT_0575;
    const ccomplex_t IT_0577 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0578 = IT_0026*IT_0577;
    const ccomplex_t IT_0579 = (0 + _Complex_I*1)*(IT_0572 + IT_0574 + 
      -IT_0576 + -IT_0578);
    const ccomplex_t IT_0580 = -IT_0579;
    const ccomplex_t IT_0581 = IT_0239*IT_0580;
    const ccomplex_t IT_0582 = m_N_2*IT_0581;
    const ccomplex_t IT_0583 = IT_0043*IT_0580;
    const ccomplex_t IT_0584 = m_N_4*IT_0583;
    const ccomplex_t IT_0585 = -IT_0584;
    const ccomplex_t IT_0586 = IT_0582 + IT_0585;
    const ccomplex_t IT_0587 = IT_0442*IT_0586;
    const ccomplex_t IT_0588 = IT_0030*IT_0239;
    const ccomplex_t IT_0589 = m_N_4*IT_0588;
    const ccomplex_t IT_0590 = m_N_2*IT_0044;
    const ccomplex_t IT_0591 = -IT_0590;
    const ccomplex_t IT_0592 = IT_0589 + IT_0591;
    const ccomplex_t IT_0593 = IT_0010*IT_0592;
    const ccomplex_t IT_0594 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0595 = IT_0014*IT_0594;
    const ccomplex_t IT_0596 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0597 = IT_0018*IT_0596;
    const ccomplex_t IT_0598 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0599 = IT_0014*IT_0598;
    const ccomplex_t IT_0600 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0601 = IT_0018*IT_0600;
    const ccomplex_t IT_0602 = IT_0595 + IT_0597 + IT_0599 + IT_0601;
    const ccomplex_t IT_0603 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0604 = IT_0023*IT_0603;
    const ccomplex_t IT_0605 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0606 = IT_0026*IT_0605;
    const ccomplex_t IT_0607 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0608 = IT_0023*IT_0607;
    const ccomplex_t IT_0609 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0610 = IT_0026*IT_0609;
    const ccomplex_t IT_0611 = -IT_0604 + -IT_0606 + -IT_0608 + -IT_0610;
    const ccomplex_t IT_0612 = IT_0602 + IT_0611;
    const ccomplex_t IT_0613 = (0 + _Complex_I*1)*IT_0612;
    const ccomplex_t IT_0614 = (-0.5)*IT_0613;
    const ccomplex_t IT_0615 = IT_0496*IT_0614;
    const ccomplex_t IT_0616 = m_N_4*IT_0615;
    const ccomplex_t IT_0617 = -IT_0616;
    const ccomplex_t IT_0618 = IT_0506*IT_0614;
    const ccomplex_t IT_0619 = m_N_3*IT_0618;
    const ccomplex_t IT_0620 = IT_0617 + IT_0619;
    const ccomplex_t IT_0621 = IT_0467*IT_0620;
    const ccomplex_t IT_0622 = (-0.5)*IT_0523 + -IT_0540 + -IT_0541 + (-0.5)
      *IT_0542 + -IT_0570 + -IT_0587 + -IT_0593 + -IT_0621;
    const ccomplex_t IT_0623 = conj(IT_0534) + conj(IT_0622);
    const ccomplex_t IT_0624 = 2*IT_0327;
    const ccomplex_t IT_0625 = (-2)*IT_0080;
    const ccomplex_t IT_0626 = IT_0624 + IT_0625;
    const ccomplex_t IT_0627 = IT_0534 + IT_0622;
    const ccomplex_t IT_0628 = 2*conj(IT_0327);
    const ccomplex_t IT_0629 = (-2)*conj(IT_0080);
    const ccomplex_t IT_0630 = IT_0628 + IT_0629;
    const ccomplex_t IT_0631 = s_12*IT_0336;
    const ccomplex_t IT_0632 = IT_0467*IT_0615;
    const ccomplex_t IT_0633 = IT_0351*IT_0519;
    const ccomplex_t IT_0634 = IT_0396*IT_0566;
    const ccomplex_t IT_0635 = IT_0442*IT_0583;
    const ccomplex_t IT_0636 = 2*IT_0632 + IT_0633 + 2*IT_0634 + 2*IT_0635;
    const ccomplex_t IT_0637 = IT_0351*IT_0517;
    const ccomplex_t IT_0638 = IT_0135 + IT_0150 + IT_0637;
    const ccomplex_t IT_0639 = IT_0442*IT_0455;
    const ccomplex_t IT_0640 = IT_0351*IT_0363;
    const ccomplex_t IT_0641 = IT_0396*IT_0436;
    const ccomplex_t IT_0642 = IT_0467*IT_0507;
    const ccomplex_t IT_0643 = 2*IT_0639 + IT_0640 + 2*IT_0641 + 2*IT_0642;
    const ccomplex_t IT_0644 = IT_0351*IT_0360;
    const ccomplex_t IT_0645 = IT_0283 + IT_0326 + IT_0644;
    const ccomplex_t IT_0646 = 2*IT_0080;
    const ccomplex_t IT_0647 = (-2)*IT_0327;
    const ccomplex_t IT_0648 = IT_0646 + IT_0647;
    const ccomplex_t IT_0649 = 2*conj(IT_0080);
    const ccomplex_t IT_0650 = (-2)*conj(IT_0327);
    const ccomplex_t IT_0651 = IT_0649 + IT_0650;
    const ccomplex_t IT_0652 = s_12*s_14;
    const ccomplex_t IT_0653 = s_24*IT_0004;
    const ccomplex_t IT_0654 = -IT_0653;
    const ccomplex_t IT_0655 = IT_0652 + IT_0654;
    const ccomplex_t IT_0656 = 0.5*IT_0113 + (-0.5)*IT_0153 + IT_0186;
    const ccomplex_t IT_0657 = -IT_0632 + (-0.5)*IT_0633 + -IT_0634 + -IT_0635;
    const ccomplex_t IT_0658 = conj(IT_0656) + conj(IT_0657);
    const ccomplex_t IT_0659 = IT_0656 + IT_0657;
    const ccomplex_t IT_0660 = (-0.5)*IT_0637;
    const ccomplex_t IT_0661 = (-0.5)*IT_0270 + 0.5*IT_0272 + IT_0315;
    const ccomplex_t IT_0662 = -IT_0642;
    const ccomplex_t IT_0663 = -IT_0639 + (-0.5)*IT_0640 + -IT_0641;
    const ccomplex_t IT_0664 = conj(IT_0661) + conj(IT_0662) + conj(IT_0663);
    const ccomplex_t IT_0665 = IT_0661 + IT_0662 + IT_0663;
    const ccomplex_t IT_0666 = (-0.5)*IT_0644;
    const ccomplex_t IT_0667 = conj(IT_0045) + conj(IT_0079) + conj(IT_0666);
    const ccomplex_t IT_0668 = IT_0045 + IT_0079 + IT_0666;
    const ccomplex_t IT_0669 = IT_0655*(IT_0244*IT_0658 + IT_0247*IT_0659 +
       IT_0333*(IT_0220 + IT_0241 + IT_0660) + IT_0330*(conj(IT_0220) + conj
      (IT_0241) + conj(IT_0660)) + IT_0626*IT_0664 + IT_0630*IT_0665 + IT_0648
      *IT_0667 + IT_0651*IT_0668);
    const ccomplex_t IT_0670 = IT_0000*IT_0336;
    const ccomplex_t IT_0671 = 2*IT_0638;
    const ccomplex_t IT_0672 = (-2)*IT_0636;
    const ccomplex_t IT_0673 = IT_0671 + IT_0672;
    const ccomplex_t IT_0674 = 2*conj(IT_0638);
    const ccomplex_t IT_0675 = (-2)*conj(IT_0636);
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = 2*IT_0636;
    const ccomplex_t IT_0678 = (-2)*IT_0638;
    const ccomplex_t IT_0679 = IT_0677 + IT_0678;
    const ccomplex_t IT_0680 = 2*conj(IT_0636);
    const ccomplex_t IT_0681 = (-2)*conj(IT_0638);
    const ccomplex_t IT_0682 = IT_0680 + IT_0681;
    const ccomplex_t IT_0683 = pow(s_24, 2);
    const ccomplex_t IT_0684 = IT_0001*IT_0683;
    const ccomplex_t IT_0685 = -IT_0046;
    const ccomplex_t IT_0686 = IT_0684 + IT_0685;
    const ccomplex_t IT_0687 = IT_0000*IT_0686;
    const ccomplex_t IT_0688 = m_N_2*IT_0686;
    const ccomplex_t IT_0689 = 2*IT_0643;
    const ccomplex_t IT_0690 = (-2)*IT_0645;
    const ccomplex_t IT_0691 = IT_0689 + IT_0690;
    const ccomplex_t IT_0692 = 2*conj(IT_0643);
    const ccomplex_t IT_0693 = (-2)*conj(IT_0645);
    const ccomplex_t IT_0694 = IT_0692 + IT_0693;
    const ccomplex_t IT_0695 = s_24*IT_0132;
    const ccomplex_t IT_0696 = IT_0001*IT_0695;
    const ccomplex_t IT_0697 = 0.5*IT_0696;
    const ccomplex_t IT_0698 = s_24 + IT_0697;
    const ccomplex_t IT_0699 = m_N_2*IT_0698;
    const ccomplex_t IT_0700 = 2*IT_0699;
    const ccomplex_t IT_0701 = 2*IT_0368;
    const ccomplex_t IT_0702 = 2*IT_0512;
    const ccomplex_t IT_0703 = IT_0701 + IT_0702;
    const ccomplex_t IT_0704 = 2*conj(IT_0368);
    const ccomplex_t IT_0705 = 2*conj(IT_0512);
    const ccomplex_t IT_0706 = IT_0704 + IT_0705;
    const ccomplex_t IT_0707 = conj(IT_0662) + conj(IT_0663);
    const ccomplex_t IT_0708 = 2*IT_0534;
    const ccomplex_t IT_0709 = 2*IT_0622;
    const ccomplex_t IT_0710 = IT_0708 + IT_0709;
    const ccomplex_t IT_0711 = IT_0662 + IT_0663;
    const ccomplex_t IT_0712 = 2*conj(IT_0534);
    const ccomplex_t IT_0713 = 2*conj(IT_0622);
    const ccomplex_t IT_0714 = IT_0712 + IT_0713;
    const ccomplex_t IT_0715 = (-2)*IT_0368;
    const ccomplex_t IT_0716 = (-2)*IT_0512;
    const ccomplex_t IT_0717 = IT_0715 + IT_0716;
    const ccomplex_t IT_0718 = (-2)*conj(IT_0368);
    const ccomplex_t IT_0719 = (-2)*conj(IT_0512);
    const ccomplex_t IT_0720 = IT_0718 + IT_0719;
    const ccomplex_t IT_0721 = (-2)*IT_0534;
    const ccomplex_t IT_0722 = (-2)*IT_0622;
    const ccomplex_t IT_0723 = IT_0721 + IT_0722;
    const ccomplex_t IT_0724 = (-2)*conj(IT_0534);
    const ccomplex_t IT_0725 = (-2)*conj(IT_0622);
    const ccomplex_t IT_0726 = IT_0724 + IT_0725;
    const ccomplex_t IT_0727 = (-0.25)*IT_0696;
    const ccomplex_t IT_0728 = s_24 + IT_0727;
    const ccomplex_t IT_0729 = m_N_2*IT_0728;
    const ccomplex_t IT_0730 = (-4)*IT_0729;
    const ccomplex_t IT_0731 = conj(IT_0220) + conj(IT_0241);
    const ccomplex_t IT_0732 = IT_0220 + IT_0241;
    const ccomplex_t IT_0733 = conj(IT_0045) + conj(IT_0079);
    const ccomplex_t IT_0734 = IT_0045 + IT_0079;
    const ccomplex_t IT_0735 = s_14*IT_0132;
    const ccomplex_t IT_0736 = IT_0001*IT_0735;
    const ccomplex_t IT_0737 = (-0.25)*IT_0736;
    const ccomplex_t IT_0738 = s_14 + IT_0737;
    const ccomplex_t IT_0739 = m_N_4*IT_0738;
    const ccomplex_t IT_0740 = (-4)*IT_0739;
    const ccomplex_t IT_0741 = 0.5*IT_0736;
    const ccomplex_t IT_0742 = s_14 + IT_0741;
    const ccomplex_t IT_0743 = m_N_4*IT_0742;
    const ccomplex_t IT_0744 = 2*IT_0743;
    const ccomplex_t IT_0745 = s_12*IT_0132;
    const ccomplex_t IT_0746 = s_14*s_24;
    const ccomplex_t IT_0747 = (-2)*IT_0746;
    const ccomplex_t IT_0748 = IT_0745 + IT_0747;
    const ccomplex_t IT_0749 = IT_0001*IT_0748;
    const ccomplex_t IT_0750 = (-0.5)*IT_0749;
    const ccomplex_t IT_0751 = s_12 + IT_0750;
    const ccomplex_t IT_0752 = 4*IT_0534;
    const ccomplex_t IT_0753 = 6*conj(IT_0368);
    const ccomplex_t IT_0754 = s_12*IT_0006;
    const ccomplex_t IT_0755 = s_12*IT_0686;
    const ccomplex_t IT_0756 = 2*IT_0645;
    const ccomplex_t IT_0757 = s_12*s_24;
    const ccomplex_t IT_0758 = s_14*IT_0046;
    const ccomplex_t IT_0759 = -IT_0758;
    const ccomplex_t IT_0760 = IT_0757 + IT_0759;
    const ccomplex_t IT_0761 = (-2)*IT_0643;
    const ccomplex_t IT_0762 = 2*conj(IT_0645);
    const ccomplex_t IT_0763 = (-2)*conj(IT_0643);
    const ccomplex_t IT_0764 = (conj(IT_0220) + conj(IT_0241) + conj(IT_0660))
      *IT_0673 + (IT_0220 + IT_0241 + IT_0660)*IT_0676 + IT_0658*IT_0679 +
       IT_0659*IT_0682 + IT_0664*IT_0691 + IT_0665*IT_0694 + IT_0667*(IT_0756 +
       IT_0761) + IT_0668*(IT_0762 + IT_0763);
    const ccomplex_t IT_0765 = m_N_4*IT_0006;
    const ccomplex_t IT_0766 = m_N_4*IT_0336;
    const ccomplex_t IT_0767 = pow(m_Z, 4);
    const ccomplex_t IT_0768 = s_12*IT_0767;
    const ccomplex_t IT_0769 = IT_0001*IT_0768;
    const ccomplex_t IT_0770 = (-0.25)*IT_0769;
    const ccomplex_t IT_0771 = IT_0746 + IT_0770;
    const ccomplex_t IT_0772 = 8*conj(IT_0660);
    const ccomplex_t IT_0773 = 8*IT_0661;
    const ccomplex_t IT_0774 = 8*conj(IT_0661);
    const ccomplex_t IT_0775 = 8*IT_0666;
    const ccomplex_t IT_0776 = (-8)*IT_0662;
    const ccomplex_t IT_0777 = 8*conj(IT_0666);
    const ccomplex_t IT_0778 = (-8)*conj(IT_0662);
    const ccomplex_t IT_0779 = -IT_0657;
    const ccomplex_t IT_0780 = -conj(IT_0657);
    const ccomplex_t IT_0781 = IT_0656*conj(IT_0660) + conj(IT_0661)*IT_0666 +
       IT_0661*conj(IT_0666) + (-0.125)*IT_0732*IT_0772 + (-0.125)*IT_0707
      *IT_0773 + (-0.125)*IT_0711*IT_0774 + IT_0733*(IT_0663 + (-0.125)*IT_0775 
      + (-0.125)*IT_0776) + IT_0734*(conj(IT_0663) + (-0.125)*IT_0777 + (-0.125)
      *IT_0778) + -(conj(IT_0220) + conj(IT_0241) + -conj(IT_0656))*(IT_0660 +
       IT_0779) + -(IT_0220 + IT_0241 + -IT_0656)*IT_0780;
    const ccomplex_t IT_0782 = (-0.5)*IT_0769;
    const ccomplex_t IT_0783 = IT_0745 + IT_0747 + IT_0782;
    const ccomplex_t IT_0784 = 4*conj(IT_0666);
    const ccomplex_t IT_0785 = 4*conj(IT_0656);
    const ccomplex_t IT_0786 = 4*conj(IT_0661);
    const ccomplex_t IT_0787 = 4*IT_0656;
    const ccomplex_t IT_0788 = (-4)*IT_0241;
    const ccomplex_t IT_0789 = 4*IT_0666;
    const ccomplex_t IT_0790 = (-4)*IT_0662;
    const ccomplex_t IT_0791 = 4*IT_0661;
    const ccomplex_t IT_0792 = (-4)*IT_0045;
    const ccomplex_t IT_0793 = IT_0656*conj(IT_0656) + IT_0661*conj(IT_0661) +
       IT_0666*conj(IT_0666) + (IT_0660 + IT_0779)*(conj(IT_0660) + IT_0780) + (
      -0.25)*IT_0711*IT_0784 + (-0.25)*IT_0732*IT_0785 + (-0.25)*IT_0734*IT_0786
       + IT_0731*(IT_0220 + (-0.25)*IT_0787 + (-0.25)*IT_0788) + IT_0707*
      (IT_0663 + (-0.25)*IT_0789 + (-0.25)*IT_0790) + IT_0733*(IT_0079 + (-0.25)
      *IT_0791 + (-0.25)*IT_0792);
    const ccomplex_t IT_0794 = IT_0000*IT_0132;
    const ccomplex_t IT_0795 = 3*IT_0794;
    const ccomplex_t IT_0796 = (-2)*IT_0656;
    const ccomplex_t IT_0797 = (-2)*conj(IT_0656);
    const ccomplex_t IT_0798 = (-2)*conj(IT_0660);
    const ccomplex_t IT_0799 = (-0.5)*IT_0661;
    const ccomplex_t IT_0800 = 2*conj(IT_0661);
    const ccomplex_t IT_0801 = (-0.5)*IT_0657;
    const ccomplex_t IT_0802 = 2*IT_0661;
    const ccomplex_t IT_0803 = (-0.5)*conj(IT_0657);
    const ccomplex_t IT_0804 = 2*conj(IT_0666);
    const ccomplex_t IT_0805 = (-0.5)*IT_0732;
    const ccomplex_t IT_0806 = 2*IT_0666;
    const ccomplex_t IT_0807 = (-0.5)*IT_0731;
    const ccomplex_t IT_0808 = 2*conj(IT_0656);
    const ccomplex_t IT_0809 = (-2)*conj(IT_0241);
    const ccomplex_t IT_0810 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0811 = IT_0808 + IT_0809 + IT_0810;
    const ccomplex_t IT_0812 = (-0.5)*IT_0711;
    const ccomplex_t IT_0813 = 2*IT_0656;
    const ccomplex_t IT_0814 = (-2)*IT_0241;
    const ccomplex_t IT_0815 = (-2)*IT_0220;
    const ccomplex_t IT_0816 = IT_0813 + IT_0814 + IT_0815;
    const ccomplex_t IT_0817 = (-0.5)*IT_0707;
    const ccomplex_t IT_0818 = 2*conj(IT_0660);
    const ccomplex_t IT_0819 = (-2)*conj(IT_0657);
    const ccomplex_t IT_0820 = IT_0818 + IT_0819;
    const ccomplex_t IT_0821 = (-0.5)*IT_0734;
    const ccomplex_t IT_0822 = 2*IT_0660;
    const ccomplex_t IT_0823 = (-2)*IT_0657;
    const ccomplex_t IT_0824 = IT_0822 + IT_0823;
    const ccomplex_t IT_0825 = (-0.5)*IT_0733;
    const ccomplex_t IT_0826 = IT_0660*conj(IT_0661) + (-0.5)*conj(IT_0666)
      *IT_0796 + (-0.5)*IT_0666*IT_0797 + IT_0798*IT_0799 + IT_0800*IT_0801 +
       IT_0802*IT_0803 + IT_0804*IT_0805 + IT_0806*IT_0807 + IT_0811*IT_0812 +
       IT_0816*IT_0817 + IT_0820*IT_0821 + IT_0824*IT_0825;
    const ccomplex_t IT_0827 = (-3)*IT_0794;
    const ccomplex_t IT_0828 = IT_0660*conj(IT_0666) + (-0.5)*conj(IT_0661)
      *IT_0796 + (-0.5)*IT_0666*IT_0798 + IT_0797*IT_0799 + IT_0801*IT_0804 +
       IT_0800*IT_0805 + IT_0803*IT_0806 + IT_0802*IT_0807 + IT_0812*IT_0820 +
       IT_0811*IT_0821 + IT_0817*IT_0824 + IT_0816*IT_0825;
    const ccomplex_t IT_0829 = IT_0007*(conj(IT_0080)*IT_0244 + IT_0080
      *IT_0247 + conj(IT_0327)*IT_0330 + IT_0327*IT_0333) + IT_0337*(IT_0244
      *IT_0513 + IT_0247*IT_0514 + IT_0623*IT_0626 + IT_0627*IT_0630) + IT_0631*
      (IT_0247*IT_0636 + IT_0244*conj(IT_0636) + IT_0333*IT_0638 + IT_0330*conj
      (IT_0638) + IT_0630*IT_0643 + IT_0626*conj(IT_0643) + conj(IT_0645)
      *IT_0648 + IT_0645*IT_0651) + IT_0669 + IT_0670*(IT_0333*IT_0643 + IT_0330
      *conj(IT_0643) + IT_0247*IT_0645 + IT_0244*conj(IT_0645) + conj(IT_0327)
      *IT_0673 + IT_0327*IT_0676 + conj(IT_0080)*IT_0679 + IT_0080*IT_0682) + 
      (conj(IT_0643)*IT_0673 + IT_0643*IT_0676 + conj(IT_0645)*IT_0679 + IT_0645
      *IT_0682)*IT_0687 + IT_0688*(IT_0513*IT_0679 + IT_0514*IT_0682 + IT_0623
      *IT_0691 + IT_0627*IT_0694) + IT_0700*(conj(IT_0657)*IT_0703 + IT_0657
      *IT_0706 + IT_0707*IT_0710 + IT_0711*IT_0714 + conj(IT_0660)*IT_0717 +
       IT_0660*IT_0720 + conj(IT_0666)*IT_0723 + IT_0666*IT_0726) + IT_0730*
      (conj(IT_0656)*IT_0717 + IT_0656*IT_0720 + conj(IT_0661)*IT_0723 + IT_0661
      *IT_0726 + IT_0703*IT_0731 + IT_0706*IT_0732 + IT_0710*IT_0733 + IT_0714
      *IT_0734) + (conj(IT_0666)*IT_0703 + IT_0666*IT_0706 + conj(IT_0660)
      *IT_0710 + IT_0660*IT_0714 + IT_0707*IT_0717 + IT_0711*IT_0720 + conj
      (IT_0657)*IT_0723 + IT_0657*IT_0726)*IT_0740 + (conj(IT_0661)*IT_0703 +
       IT_0661*IT_0706 + conj(IT_0656)*IT_0710 + IT_0656*IT_0714 + IT_0723
      *IT_0731 + IT_0726*IT_0732 + IT_0717*IT_0733 + IT_0720*IT_0734)*IT_0744 +
       4*IT_0751*(IT_0513*IT_0514 + IT_0623*(IT_0622 + 0.25*IT_0752)) + 6
      *IT_0000*(IT_0514*IT_0623 + IT_0627*(conj(IT_0512) + 0.166666666666667
      *IT_0753)) + (conj(IT_0187)*IT_0188 + conj(IT_0242)*IT_0330 + IT_0242
      *IT_0332 + conj(IT_0327)*IT_0624 + conj(IT_0080)*IT_0648 + IT_0080*IT_0650
      )*IT_0754 + IT_0755*(conj(IT_0638)*IT_0671 + conj(IT_0636)*IT_0679 +
       IT_0636*IT_0681 + conj(IT_0643)*IT_0691 + IT_0643*IT_0693 + conj(IT_0645)
      *IT_0756) + IT_0760*IT_0764 + (IT_0244*IT_0623 + conj(IT_0512)*IT_0626 +
       IT_0247*IT_0627 + IT_0512*IT_0630 + conj(IT_0327)*IT_0701 + IT_0327
      *IT_0704 + conj(IT_0080)*IT_0715 + IT_0080*IT_0718)*IT_0765 + (IT_0623
      *IT_0679 + IT_0627*IT_0682 + conj(IT_0512)*IT_0691 + IT_0512*IT_0694 +
       conj(IT_0643)*IT_0701 + IT_0643*IT_0704 + conj(IT_0645)*IT_0715 + IT_0645
      *IT_0718)*IT_0766 + (-8)*IT_0771*IT_0781 + (-4)*IT_0783*IT_0793 + (-2)
      *IT_0795*IT_0826 + (-2)*IT_0827*IT_0828;
    return create_ccomplex_return(IT_0829);
}

