#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_h_Z(
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
    const ccomplex_t IT_0000 = pow(m_N_2, 2);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0000;
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0012 + IT_0016 + -IT_0021 + -IT_0024;
    const ccomplex_t IT_0026 = cos(alpha);
    const ccomplex_t IT_0027 = sin(alpha);
    const ccomplex_t IT_0028 = e_em*IT_0009*IT_0018*(IT_0013*IT_0026 + 
      -IT_0007*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0025*IT_0029;
    const ccomplex_t IT_0031 = pow(m_Z, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0031 + 
      -reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = IT_0010*IT_0015;
    const ccomplex_t IT_0036 = IT_0011*IT_0014;
    const ccomplex_t IT_0037 = IT_0019*IT_0023;
    const ccomplex_t IT_0038 = IT_0020*IT_0022;
    const ccomplex_t IT_0039 = IT_0035 + -IT_0036 + -IT_0037 + IT_0038;
    const ccomplex_t IT_0040 = e_em*IT_0009*IT_0018*(IT_0007*IT_0026 + IT_0013
      *IT_0027);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = IT_0009*IT_0017;
    const ccomplex_t IT_0047 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = IT_0008*IT_0018;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = IT_0049*IT_0051;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0048 + IT_0050 + 
      -IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0009*IT_0026;
    const ccomplex_t IT_0057 = IT_0015*IT_0056;
    const ccomplex_t IT_0058 = IT_0009*IT_0027;
    const ccomplex_t IT_0059 = IT_0011*IT_0058;
    const ccomplex_t IT_0060 = IT_0018*IT_0026;
    const ccomplex_t IT_0061 = IT_0023*IT_0060;
    const ccomplex_t IT_0062 = IT_0018*IT_0027;
    const ccomplex_t IT_0063 = IT_0020*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0057 + IT_0059 + 
      -IT_0061 + -IT_0063);
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0055*IT_0065;
    const ccomplex_t IT_0067 = pow(m_h, 2);
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0067 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0066*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0072 = IT_0056*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0074 = IT_0058*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0056*IT_0075;
    const ccomplex_t IT_0077 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0058*IT_0077;
    const ccomplex_t IT_0079 = IT_0072 + IT_0074 + IT_0076 + IT_0078;
    const ccomplex_t IT_0080 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0081 = IT_0060*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0085 = IT_0060*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0062*IT_0086;
    const ccomplex_t IT_0088 = -IT_0081 + -IT_0083 + -IT_0085 + -IT_0087;
    const ccomplex_t IT_0089 = IT_0079 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0093 = IT_0046*IT_0092;
    const ccomplex_t IT_0094 = IT_0049*IT_0092;
    const ccomplex_t IT_0095 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0096 = IT_0046*IT_0095;
    const ccomplex_t IT_0097 = IT_0049*IT_0095;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0093 + IT_0094 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = cpow((-2)*s_23 + IT_0000 + IT_0067 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0105 = IT_0056*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0107 = IT_0058*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0109 = IT_0056*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0111 = IT_0058*IT_0110;
    const ccomplex_t IT_0112 = IT_0105 + IT_0107 + IT_0109 + IT_0111;
    const ccomplex_t IT_0113 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0114 = IT_0060*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0116 = IT_0062*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0118 = IT_0060*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0120 = IT_0062*IT_0119;
    const ccomplex_t IT_0121 = -IT_0114 + -IT_0116 + -IT_0118 + -IT_0120;
    const ccomplex_t IT_0122 = IT_0112 + IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0126 = IT_0046*IT_0125;
    const ccomplex_t IT_0127 = IT_0049*IT_0125;
    const ccomplex_t IT_0128 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0129 = IT_0046*IT_0128;
    const ccomplex_t IT_0130 = IT_0049*IT_0128;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0126 + IT_0127 + 
      -IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0124*IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_23 + IT_0000 + IT_0067 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = (-0.5)*IT_0034 + (-0.5)*IT_0045 + -IT_0070 + (
      -2)*IT_0103 + (-2)*IT_0136;
    const ccomplex_t IT_0138 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0139 = IT_0010*IT_0138;
    const ccomplex_t IT_0140 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0141 = IT_0014*IT_0140;
    const ccomplex_t IT_0142 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0143 = IT_0019*IT_0142;
    const ccomplex_t IT_0144 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0145 = IT_0022*IT_0144;
    const ccomplex_t IT_0146 = -IT_0139 + IT_0141 + IT_0143 + -IT_0145;
    const ccomplex_t IT_0147 = IT_0041*IT_0146;
    const ccomplex_t IT_0148 = IT_0044*IT_0147;
    const ccomplex_t IT_0149 = IT_0056*IT_0138;
    const ccomplex_t IT_0150 = IT_0058*IT_0140;
    const ccomplex_t IT_0151 = IT_0060*IT_0142;
    const ccomplex_t IT_0152 = IT_0062*IT_0144;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + IT_0150 + 
      -IT_0151 + -IT_0152);
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0055*IT_0154;
    const ccomplex_t IT_0156 = IT_0069*IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0000 + IT_0067 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0160 = IT_0056*IT_0159;
    const ccomplex_t IT_0161 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0162 = IT_0058*IT_0161;
    const ccomplex_t IT_0163 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0164 = IT_0056*IT_0163;
    const ccomplex_t IT_0165 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0166 = IT_0058*IT_0165;
    const ccomplex_t IT_0167 = IT_0160 + IT_0162 + IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0169 = IT_0060*IT_0168;
    const ccomplex_t IT_0170 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0171 = IT_0062*IT_0170;
    const ccomplex_t IT_0172 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0173 = IT_0060*IT_0172;
    const ccomplex_t IT_0174 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0175 = IT_0062*IT_0174;
    const ccomplex_t IT_0176 = -IT_0169 + -IT_0171 + -IT_0173 + -IT_0175;
    const ccomplex_t IT_0177 = IT_0167 + IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0181 = IT_0046*IT_0180;
    const ccomplex_t IT_0182 = IT_0049*IT_0180;
    const ccomplex_t IT_0183 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0184 = IT_0046*IT_0183;
    const ccomplex_t IT_0185 = IT_0049*IT_0183;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0181 + IT_0182 + 
      -IT_0184 + -IT_0185);
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = IT_0179*IT_0187;
    const ccomplex_t IT_0189 = IT_0158*IT_0188;
    const ccomplex_t IT_0190 = IT_0010*IT_0140;
    const ccomplex_t IT_0191 = IT_0014*IT_0138;
    const ccomplex_t IT_0192 = IT_0019*IT_0144;
    const ccomplex_t IT_0193 = IT_0022*IT_0142;
    const ccomplex_t IT_0194 = -IT_0190 + -IT_0191 + IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = IT_0029*IT_0194;
    const ccomplex_t IT_0196 = IT_0033*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0148 + -IT_0156 + (-2)*IT_0189 + 0.5
      *IT_0196;
    const ccomplex_t IT_0198 = 2*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0054;
    const ccomplex_t IT_0200 = IT_0154*IT_0199;
    const ccomplex_t IT_0201 = IT_0069*IT_0200;
    const ccomplex_t IT_0202 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0203 = IT_0056*IT_0202;
    const ccomplex_t IT_0204 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0205 = IT_0058*IT_0204;
    const ccomplex_t IT_0206 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0207 = IT_0056*IT_0206;
    const ccomplex_t IT_0208 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0209 = IT_0058*IT_0208;
    const ccomplex_t IT_0210 = IT_0203 + IT_0205 + IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0212 = IT_0060*IT_0211;
    const ccomplex_t IT_0213 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0214 = IT_0062*IT_0213;
    const ccomplex_t IT_0215 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0216 = IT_0060*IT_0215;
    const ccomplex_t IT_0217 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0218 = IT_0062*IT_0217;
    const ccomplex_t IT_0219 = -IT_0212 + -IT_0214 + -IT_0216 + -IT_0218;
    const ccomplex_t IT_0220 = IT_0210 + IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0224 = IT_0046*IT_0223;
    const ccomplex_t IT_0225 = IT_0049*IT_0223;
    const ccomplex_t IT_0226 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0227 = IT_0046*IT_0226;
    const ccomplex_t IT_0228 = IT_0049*IT_0226;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0224 + IT_0225 + 
      -IT_0227 + -IT_0228);
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0222*IT_0230;
    const ccomplex_t IT_0232 = IT_0102*IT_0231;
    const ccomplex_t IT_0233 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0234 = IT_0056*IT_0233;
    const ccomplex_t IT_0235 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0236 = IT_0058*IT_0235;
    const ccomplex_t IT_0237 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0238 = IT_0056*IT_0237;
    const ccomplex_t IT_0239 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0240 = IT_0058*IT_0239;
    const ccomplex_t IT_0241 = IT_0234 + IT_0236 + IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0243 = IT_0060*IT_0242;
    const ccomplex_t IT_0244 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0245 = IT_0062*IT_0244;
    const ccomplex_t IT_0246 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0247 = IT_0060*IT_0246;
    const ccomplex_t IT_0248 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0249 = IT_0062*IT_0248;
    const ccomplex_t IT_0250 = -IT_0243 + -IT_0245 + -IT_0247 + -IT_0249;
    const ccomplex_t IT_0251 = IT_0241 + IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = (-0.5)*IT_0252;
    const ccomplex_t IT_0254 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0255 = IT_0046*IT_0254;
    const ccomplex_t IT_0256 = IT_0049*IT_0254;
    const ccomplex_t IT_0257 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0258 = IT_0046*IT_0257;
    const ccomplex_t IT_0259 = IT_0049*IT_0257;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*(IT_0255 + IT_0256 + 
      -IT_0258 + -IT_0259);
    const ccomplex_t IT_0261 = (-0.5)*IT_0260;
    const ccomplex_t IT_0262 = IT_0253*IT_0261;
    const ccomplex_t IT_0263 = IT_0135*IT_0262;
    const ccomplex_t IT_0264 = (-0.5)*IT_0148 + (-0.5)*IT_0196 + -IT_0201 + (
      -2)*IT_0232 + (-2)*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0264;
    const ccomplex_t IT_0266 = IT_0198 + IT_0265;
    const ccomplex_t IT_0267 = 2*conj(IT_0197);
    const ccomplex_t IT_0268 = (-2)*conj(IT_0264);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = IT_0065*IT_0199;
    const ccomplex_t IT_0271 = IT_0069*IT_0270;
    const ccomplex_t IT_0272 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0273 = IT_0056*IT_0272;
    const ccomplex_t IT_0274 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0275 = IT_0058*IT_0274;
    const ccomplex_t IT_0276 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0277 = IT_0056*IT_0276;
    const ccomplex_t IT_0278 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0279 = IT_0058*IT_0278;
    const ccomplex_t IT_0280 = IT_0273 + IT_0275 + IT_0277 + IT_0279;
    const ccomplex_t IT_0281 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0282 = IT_0060*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0284 = IT_0062*IT_0283;
    const ccomplex_t IT_0285 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0286 = IT_0060*IT_0285;
    const ccomplex_t IT_0287 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0288 = IT_0062*IT_0287;
    const ccomplex_t IT_0289 = -IT_0282 + -IT_0284 + -IT_0286 + -IT_0288;
    const ccomplex_t IT_0290 = IT_0280 + IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = (-0.5)*IT_0291;
    const ccomplex_t IT_0293 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0294 = IT_0046*IT_0293;
    const ccomplex_t IT_0295 = IT_0049*IT_0293;
    const ccomplex_t IT_0296 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0297 = IT_0046*IT_0296;
    const ccomplex_t IT_0298 = IT_0049*IT_0296;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*(IT_0294 + IT_0295 + 
      -IT_0297 + -IT_0298);
    const ccomplex_t IT_0300 = (-0.5)*IT_0299;
    const ccomplex_t IT_0301 = IT_0292*IT_0300;
    const ccomplex_t IT_0302 = IT_0158*IT_0301;
    const ccomplex_t IT_0303 = 0.5*IT_0034 + 0.5*IT_0045 + -IT_0271 + (-2)
      *IT_0302;
    const ccomplex_t IT_0304 = 2*IT_0264;
    const ccomplex_t IT_0305 = (-2)*IT_0197;
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = 2*conj(IT_0264);
    const ccomplex_t IT_0308 = (-2)*conj(IT_0197);
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = s_12 + IT_0311;
    const ccomplex_t IT_0313 = s_12*IT_0312;
    const ccomplex_t IT_0314 = cpow((-2)*s_13 + IT_0000 + IT_0067 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0100*IT_0315;
    const ccomplex_t IT_0317 = cpow((-2)*s_13 + IT_0067 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = IT_0066*IT_0318;
    const ccomplex_t IT_0320 = cpow((-2)*s_13 + IT_0000 + IT_0067 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = IT_0133*IT_0321;
    const ccomplex_t IT_0323 = 0.5*IT_0034 + 0.5*IT_0045 + 2*IT_0316 + IT_0319
       + 2*IT_0322;
    const ccomplex_t IT_0324 = 2*IT_0137;
    const ccomplex_t IT_0325 = (-2)*IT_0303;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = 2*conj(IT_0137);
    const ccomplex_t IT_0328 = (-2)*conj(IT_0303);
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = cpow((-2)*s_13 + IT_0000 + IT_0067 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = IT_0188*IT_0331;
    const ccomplex_t IT_0333 = IT_0155*IT_0318;
    const ccomplex_t IT_0334 = (-0.5)*IT_0148 + (-0.5)*IT_0196 + 2*IT_0332 +
       IT_0333;
    const ccomplex_t IT_0335 = IT_0200*IT_0318;
    const ccomplex_t IT_0336 = IT_0231*IT_0315;
    const ccomplex_t IT_0337 = IT_0262*IT_0321;
    const ccomplex_t IT_0338 = 0.5*IT_0148 + 0.5*IT_0196 + IT_0335 + 2*IT_0336
       + 2*IT_0337;
    const ccomplex_t IT_0339 = IT_0270*IT_0318;
    const ccomplex_t IT_0340 = IT_0301*IT_0331;
    const ccomplex_t IT_0341 = (-0.5)*IT_0034 + (-0.5)*IT_0045 + IT_0339 + 2
      *IT_0340;
    const ccomplex_t IT_0342 = 2*IT_0303;
    const ccomplex_t IT_0343 = (-2)*IT_0137;
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 2*conj(IT_0303);
    const ccomplex_t IT_0346 = (-2)*conj(IT_0137);
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = s_12*s_14;
    const ccomplex_t IT_0349 = s_24*IT_0000;
    const ccomplex_t IT_0350 = -IT_0349;
    const ccomplex_t IT_0351 = IT_0348 + IT_0350;
    const ccomplex_t IT_0352 = -IT_0316 + (-0.5)*IT_0319 + -IT_0322;
    const ccomplex_t IT_0353 = 0.5*IT_0070 + IT_0103 + IT_0136;
    const ccomplex_t IT_0354 = conj(IT_0352) + conj(IT_0353);
    const ccomplex_t IT_0355 = IT_0352 + IT_0353;
    const ccomplex_t IT_0356 = -IT_0332 + (-0.5)*IT_0333;
    const ccomplex_t IT_0357 = 0.5*IT_0156;
    const ccomplex_t IT_0358 = conj(IT_0189) + conj(IT_0357);
    const ccomplex_t IT_0359 = IT_0189 + IT_0357;
    const ccomplex_t IT_0360 = (-0.5)*IT_0335 + -IT_0336 + -IT_0337;
    const ccomplex_t IT_0361 = 0.5*IT_0201 + IT_0232 + IT_0263;
    const ccomplex_t IT_0362 = conj(IT_0360) + conj(IT_0361);
    const ccomplex_t IT_0363 = IT_0360 + IT_0361;
    const ccomplex_t IT_0364 = (-0.5)*IT_0339 + -IT_0340;
    const ccomplex_t IT_0365 = 0.5*IT_0271;
    const ccomplex_t IT_0366 = conj(IT_0302) + conj(IT_0365);
    const ccomplex_t IT_0367 = IT_0302 + IT_0365;
    const ccomplex_t IT_0368 = IT_0351*(IT_0326*IT_0354 + IT_0329*IT_0355 +
       IT_0266*(conj(IT_0356) + IT_0358) + IT_0269*(IT_0356 + IT_0359) + IT_0306
      *IT_0362 + IT_0309*IT_0363 + IT_0344*(conj(IT_0364) + IT_0366) + IT_0347*
      (IT_0364 + IT_0367));
    const ccomplex_t IT_0369 = IT_0000*IT_0312;
    const ccomplex_t IT_0370 = 2*IT_0338;
    const ccomplex_t IT_0371 = (-2)*IT_0334;
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 2*conj(IT_0338);
    const ccomplex_t IT_0374 = (-2)*conj(IT_0334);
    const ccomplex_t IT_0375 = IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = 2*IT_0334;
    const ccomplex_t IT_0377 = (-2)*IT_0338;
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = 2*conj(IT_0334);
    const ccomplex_t IT_0380 = (-2)*conj(IT_0338);
    const ccomplex_t IT_0381 = IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = pow(s_24, 2);
    const ccomplex_t IT_0383 = IT_0001*IT_0382;
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = IT_0000 + IT_0384;
    const ccomplex_t IT_0386 = IT_0000*IT_0385;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = m_N_2*IT_0385;
    const ccomplex_t IT_0389 = -IT_0388;
    const ccomplex_t IT_0390 = m_W*e_em;
    const ccomplex_t IT_0391 = cpow(IT_0008, -4);
    const ccomplex_t IT_0392 = IT_0013*IT_0017*IT_0026*IT_0391;
    const ccomplex_t IT_0393 = IT_0007*IT_0017*IT_0027*IT_0391;
    const ccomplex_t IT_0394 = cpow(IT_0017, -3);
    const ccomplex_t IT_0395 = IT_0013*IT_0026*IT_0394;
    const ccomplex_t IT_0396 = IT_0007*IT_0027*IT_0394;
    const ccomplex_t IT_0397 = cpow(IT_0008, -2);
    const ccomplex_t IT_0398 = IT_0013*IT_0018*IT_0026*IT_0397;
    const ccomplex_t IT_0399 = IT_0007*IT_0018*IT_0027*IT_0397;
    const ccomplex_t IT_0400 = cpow(IT_0017, -2);
    const ccomplex_t IT_0401 = cpow(IT_0397 + IT_0400, -1);
    const ccomplex_t IT_0402 = (IT_0392 + -IT_0393 + IT_0395 + -IT_0396 + 2
      *IT_0398 + (-2)*IT_0399)*IT_0401;
    const ccomplex_t IT_0403 = IT_0390*IT_0402;
    const ccomplex_t IT_0404 = (0 + _Complex_I*1)*IT_0403;
    const ccomplex_t IT_0405 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0406 = IT_0404*IT_0405;
    const ccomplex_t IT_0407 = (0 + _Complex_I*1)*IT_0406;
    const ccomplex_t IT_0408 = IT_0199*IT_0407;
    const ccomplex_t IT_0409 = m_N_2*IT_0200;
    const ccomplex_t IT_0410 = m_N_2*IT_0155;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = IT_0409 + IT_0411;
    const ccomplex_t IT_0413 = IT_0069*IT_0412;
    const ccomplex_t IT_0414 = IT_0179*IT_0300;
    const ccomplex_t IT_0415 = m_N_1*IT_0414;
    const ccomplex_t IT_0416 = m_N_2*IT_0188;
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = IT_0415 + IT_0417;
    const ccomplex_t IT_0419 = IT_0158*IT_0418;
    const ccomplex_t IT_0420 = m_N_2*IT_0066;
    const ccomplex_t IT_0421 = m_N_2*IT_0270;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0420 + IT_0422;
    const ccomplex_t IT_0424 = IT_0318*IT_0423;
    const ccomplex_t IT_0425 = IT_0091*IT_0230;
    const ccomplex_t IT_0426 = m_N_3*IT_0425;
    const ccomplex_t IT_0427 = m_N_2*IT_0100;
    const ccomplex_t IT_0428 = -IT_0427;
    const ccomplex_t IT_0429 = IT_0426 + IT_0428;
    const ccomplex_t IT_0430 = IT_0315*IT_0429;
    const ccomplex_t IT_0431 = IT_0124*IT_0261;
    const ccomplex_t IT_0432 = m_N_4*IT_0431;
    const ccomplex_t IT_0433 = m_N_2*IT_0133;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = IT_0321*IT_0435;
    const ccomplex_t IT_0437 = 0.5*IT_0408 + 0.5*IT_0413 + IT_0419 + (-0.5)
      *IT_0424 + IT_0430 + IT_0436;
    const ccomplex_t IT_0438 = IT_0055*IT_0407;
    const ccomplex_t IT_0439 = IT_0187*IT_0292;
    const ccomplex_t IT_0440 = m_N_1*IT_0439;
    const ccomplex_t IT_0441 = m_N_2*IT_0301;
    const ccomplex_t IT_0442 = -IT_0441;
    const ccomplex_t IT_0443 = IT_0440 + IT_0442;
    const ccomplex_t IT_0444 = IT_0331*IT_0443;
    const ccomplex_t IT_0445 = IT_0099*IT_0222;
    const ccomplex_t IT_0446 = m_N_3*IT_0445;
    const ccomplex_t IT_0447 = m_N_2*IT_0231;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = IT_0446 + IT_0448;
    const ccomplex_t IT_0450 = IT_0102*IT_0449;
    const ccomplex_t IT_0451 = IT_0132*IT_0253;
    const ccomplex_t IT_0452 = m_N_4*IT_0451;
    const ccomplex_t IT_0453 = m_N_2*IT_0262;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = IT_0452 + IT_0454;
    const ccomplex_t IT_0456 = IT_0135*IT_0455;
    const ccomplex_t IT_0457 = 0.5*IT_0413 + (-0.5)*IT_0424 + (-0.5)*IT_0438 +
       -IT_0444 + -IT_0450 + -IT_0456;
    const ccomplex_t IT_0458 = conj(IT_0437) + conj(IT_0457);
    const ccomplex_t IT_0459 = IT_0437 + IT_0457;
    const ccomplex_t IT_0460 = IT_0069*IT_0423;
    const ccomplex_t IT_0461 = IT_0318*IT_0412;
    const ccomplex_t IT_0462 = IT_0158*IT_0443;
    const ccomplex_t IT_0463 = IT_0315*IT_0449;
    const ccomplex_t IT_0464 = IT_0321*IT_0455;
    const ccomplex_t IT_0465 = 0.5*IT_0438 + 0.5*IT_0460 + (-0.5)*IT_0461 +
       IT_0462 + IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = IT_0331*IT_0418;
    const ccomplex_t IT_0467 = IT_0102*IT_0429;
    const ccomplex_t IT_0468 = IT_0135*IT_0435;
    const ccomplex_t IT_0469 = (-0.5)*IT_0408 + 0.5*IT_0460 + (-0.5)*IT_0461 +
       -IT_0466 + -IT_0467 + -IT_0468;
    const ccomplex_t IT_0470 = conj(IT_0465) + conj(IT_0469);
    const ccomplex_t IT_0471 = 2*IT_0341;
    const ccomplex_t IT_0472 = (-2)*IT_0323;
    const ccomplex_t IT_0473 = IT_0471 + IT_0472;
    const ccomplex_t IT_0474 = IT_0465 + IT_0469;
    const ccomplex_t IT_0475 = 2*conj(IT_0341);
    const ccomplex_t IT_0476 = (-2)*conj(IT_0323);
    const ccomplex_t IT_0477 = IT_0475 + IT_0476;
    const ccomplex_t IT_0478 = s_24*IT_0031;
    const ccomplex_t IT_0479 = IT_0001*IT_0478;
    const ccomplex_t IT_0480 = 0.5*IT_0479;
    const ccomplex_t IT_0481 = s_24 + IT_0480;
    const ccomplex_t IT_0482 = m_N_2*IT_0481;
    const ccomplex_t IT_0483 = 2*IT_0482;
    const ccomplex_t IT_0484 = 2*IT_0437;
    const ccomplex_t IT_0485 = 2*IT_0457;
    const ccomplex_t IT_0486 = IT_0484 + IT_0485;
    const ccomplex_t IT_0487 = 2*conj(IT_0437);
    const ccomplex_t IT_0488 = 2*conj(IT_0457);
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = 2*IT_0465;
    const ccomplex_t IT_0491 = 2*IT_0469;
    const ccomplex_t IT_0492 = IT_0490 + IT_0491;
    const ccomplex_t IT_0493 = 2*conj(IT_0465);
    const ccomplex_t IT_0494 = 2*conj(IT_0469);
    const ccomplex_t IT_0495 = IT_0493 + IT_0494;
    const ccomplex_t IT_0496 = (-2)*IT_0437;
    const ccomplex_t IT_0497 = (-2)*IT_0457;
    const ccomplex_t IT_0498 = IT_0496 + IT_0497;
    const ccomplex_t IT_0499 = (-2)*conj(IT_0437);
    const ccomplex_t IT_0500 = (-2)*conj(IT_0457);
    const ccomplex_t IT_0501 = IT_0499 + IT_0500;
    const ccomplex_t IT_0502 = (-2)*IT_0465;
    const ccomplex_t IT_0503 = (-2)*IT_0469;
    const ccomplex_t IT_0504 = IT_0502 + IT_0503;
    const ccomplex_t IT_0505 = (-2)*conj(IT_0465);
    const ccomplex_t IT_0506 = (-2)*conj(IT_0469);
    const ccomplex_t IT_0507 = IT_0505 + IT_0506;
    const ccomplex_t IT_0508 = s_14*IT_0031;
    const ccomplex_t IT_0509 = IT_0001*IT_0508;
    const ccomplex_t IT_0510 = (-0.25)*IT_0509;
    const ccomplex_t IT_0511 = s_14 + IT_0510;
    const ccomplex_t IT_0512 = m_N_2*IT_0511;
    const ccomplex_t IT_0513 = (-4)*IT_0512;
    const ccomplex_t IT_0514 = 0.5*IT_0509;
    const ccomplex_t IT_0515 = s_14 + IT_0514;
    const ccomplex_t IT_0516 = m_N_2*IT_0515;
    const ccomplex_t IT_0517 = 2*IT_0516;
    const ccomplex_t IT_0518 = (-0.25)*IT_0479;
    const ccomplex_t IT_0519 = s_24 + IT_0518;
    const ccomplex_t IT_0520 = m_N_2*IT_0519;
    const ccomplex_t IT_0521 = (-4)*IT_0520;
    const ccomplex_t IT_0522 = s_12*IT_0031;
    const ccomplex_t IT_0523 = s_14*s_24;
    const ccomplex_t IT_0524 = (-2)*IT_0523;
    const ccomplex_t IT_0525 = IT_0522 + IT_0524;
    const ccomplex_t IT_0526 = IT_0001*IT_0525;
    const ccomplex_t IT_0527 = (-0.5)*IT_0526;
    const ccomplex_t IT_0528 = s_12 + IT_0527;
    const ccomplex_t IT_0529 = 4*IT_0465;
    const ccomplex_t IT_0530 = 6*conj(IT_0437);
    const ccomplex_t IT_0531 = m_N_2*IT_0005;
    const ccomplex_t IT_0532 = s_12*IT_0005;
    const ccomplex_t IT_0533 = m_N_2*IT_0312;
    const ccomplex_t IT_0534 = IT_0269*IT_0437 + IT_0266*conj(IT_0437) + 
      (IT_0342 + IT_0343 + IT_0376 + IT_0377)*IT_0470 + conj(IT_0457)*(IT_0198 +
       IT_0265 + IT_0471 + IT_0472) + (IT_0345 + IT_0346 + IT_0379 + IT_0380)
      *IT_0474 + IT_0457*(IT_0267 + IT_0268 + IT_0475 + IT_0476) + conj(IT_0341)
      *IT_0484 + IT_0341*IT_0487 + conj(IT_0323)*IT_0496 + IT_0323*IT_0499;
    const ccomplex_t IT_0535 = s_12*IT_0385;
    const ccomplex_t IT_0536 = -IT_0535;
    const ccomplex_t IT_0537 = 2*IT_0323;
    const ccomplex_t IT_0538 = s_12*s_24;
    const ccomplex_t IT_0539 = s_14*IT_0000;
    const ccomplex_t IT_0540 = -IT_0539;
    const ccomplex_t IT_0541 = IT_0538 + IT_0540;
    const ccomplex_t IT_0542 = (-2)*IT_0341;
    const ccomplex_t IT_0543 = 2*conj(IT_0323);
    const ccomplex_t IT_0544 = (-2)*conj(IT_0341);
    const ccomplex_t IT_0545 = IT_0362*IT_0372 + IT_0363*IT_0375 + (conj
      (IT_0356) + IT_0358)*IT_0378 + (IT_0356 + IT_0359)*IT_0381 + (conj(IT_0364
      ) + IT_0366)*IT_0473 + (IT_0364 + IT_0367)*IT_0477 + IT_0354*(IT_0537 +
       IT_0542) + IT_0355*(IT_0543 + IT_0544);
    const ccomplex_t IT_0546 = IT_0000*IT_0031;
    const ccomplex_t IT_0547 = (-3)*IT_0546;
    const ccomplex_t IT_0548 = 2*IT_0357;
    const ccomplex_t IT_0549 = 2*IT_0189;
    const ccomplex_t IT_0550 = (-2)*IT_0361;
    const ccomplex_t IT_0551 = IT_0548 + IT_0549 + IT_0550;
    const ccomplex_t IT_0552 = 2*conj(IT_0357);
    const ccomplex_t IT_0553 = 2*conj(IT_0189);
    const ccomplex_t IT_0554 = (-2)*conj(IT_0361);
    const ccomplex_t IT_0555 = IT_0552 + IT_0553 + IT_0554;
    const ccomplex_t IT_0556 = 2*IT_0302;
    const ccomplex_t IT_0557 = 2*IT_0365;
    const ccomplex_t IT_0558 = IT_0556 + IT_0557;
    const ccomplex_t IT_0559 = 2*conj(IT_0302);
    const ccomplex_t IT_0560 = 2*conj(IT_0365);
    const ccomplex_t IT_0561 = IT_0559 + IT_0560;
    const ccomplex_t IT_0562 = 2*IT_0360;
    const ccomplex_t IT_0563 = (-2)*IT_0356;
    const ccomplex_t IT_0564 = IT_0562 + IT_0563;
    const ccomplex_t IT_0565 = 2*conj(IT_0360);
    const ccomplex_t IT_0566 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0567 = IT_0565 + IT_0566;
    const ccomplex_t IT_0568 = (-2)*IT_0357;
    const ccomplex_t IT_0569 = (-2)*IT_0189;
    const ccomplex_t IT_0570 = IT_0568 + IT_0569;
    const ccomplex_t IT_0571 = (-2)*conj(IT_0357);
    const ccomplex_t IT_0572 = (-2)*conj(IT_0189);
    const ccomplex_t IT_0573 = IT_0571 + IT_0572;
    const ccomplex_t IT_0574 = 2*IT_0352;
    const ccomplex_t IT_0575 = 2*conj(IT_0352);
    const ccomplex_t IT_0576 = (-2)*IT_0360;
    const ccomplex_t IT_0577 = (-2)*conj(IT_0360);
    const ccomplex_t IT_0578 = conj(IT_0353)*IT_0551 + IT_0353*IT_0555 + conj
      (IT_0361)*IT_0558 + IT_0361*IT_0561 + conj(IT_0364)*IT_0564 + IT_0364
      *IT_0567 + IT_0366*IT_0570 + IT_0367*IT_0573 + conj(IT_0356)*IT_0574 +
       IT_0356*IT_0575 + conj(IT_0352)*IT_0576 + IT_0352*IT_0577;
    const ccomplex_t IT_0579 = 3*IT_0546;
    const ccomplex_t IT_0580 = conj(IT_0364)*IT_0551 + IT_0364*IT_0555 + conj
      (IT_0356)*IT_0558 + IT_0356*IT_0561 + conj(IT_0353)*IT_0564 + IT_0353
      *IT_0567 + conj(IT_0352)*IT_0570 + IT_0352*IT_0573 + conj(IT_0361)*IT_0574
       + IT_0361*IT_0575 + IT_0366*IT_0576 + IT_0367*IT_0577;
    const ccomplex_t IT_0581 = pow(m_Z, 4);
    const ccomplex_t IT_0582 = s_12*IT_0581;
    const ccomplex_t IT_0583 = IT_0001*IT_0582;
    const ccomplex_t IT_0584 = (-0.5)*IT_0583;
    const ccomplex_t IT_0585 = IT_0522 + IT_0524 + IT_0584;
    const ccomplex_t IT_0586 = -IT_0353;
    const ccomplex_t IT_0587 = IT_0367 + IT_0586;
    const ccomplex_t IT_0588 = -IT_0356 + IT_0360;
    const ccomplex_t IT_0589 = -conj(IT_0356);
    const ccomplex_t IT_0590 = conj(IT_0360) + IT_0589;
    const ccomplex_t IT_0591 = IT_0352 + -IT_0364;
    const ccomplex_t IT_0592 = conj(IT_0352) + -conj(IT_0364);
    const ccomplex_t IT_0593 = (-4)*IT_0302;
    const ccomplex_t IT_0594 = IT_0358 + -conj(IT_0361);
    const ccomplex_t IT_0595 = -IT_0361;
    const ccomplex_t IT_0596 = (-0.25)*IT_0583;
    const ccomplex_t IT_0597 = IT_0523 + IT_0596;
    const ccomplex_t IT_0598 = IT_0006*(conj(IT_0137)*IT_0266 + IT_0137
      *IT_0269 + conj(IT_0303)*IT_0306 + IT_0303*IT_0309) + IT_0313*(conj
      (IT_0323)*IT_0326 + IT_0323*IT_0329 + IT_0269*IT_0334 + IT_0266*conj
      (IT_0334) + IT_0309*IT_0338 + IT_0306*conj(IT_0338) + conj(IT_0341)
      *IT_0344 + IT_0341*IT_0347) + IT_0368 + IT_0369*(IT_0269*IT_0323 + IT_0266
      *conj(IT_0323) + IT_0309*IT_0341 + IT_0306*conj(IT_0341) + conj(IT_0303)
      *IT_0372 + IT_0303*IT_0375 + conj(IT_0137)*IT_0378 + IT_0137*IT_0381) + 
      (conj(IT_0341)*IT_0372 + IT_0341*IT_0375 + conj(IT_0323)*IT_0378 + IT_0323
      *IT_0381)*IT_0387 + IT_0389*(IT_0378*IT_0458 + IT_0381*IT_0459 + IT_0470
      *IT_0473 + IT_0474*IT_0477) + IT_0483*(conj(IT_0356)*IT_0486 + IT_0356
      *IT_0489 + conj(IT_0364)*IT_0492 + IT_0364*IT_0495 + conj(IT_0360)*IT_0498
       + IT_0360*IT_0501 + conj(IT_0352)*IT_0504 + IT_0352*IT_0507) + (conj
      (IT_0352)*IT_0486 + IT_0352*IT_0489 + conj(IT_0360)*IT_0492 + IT_0360
      *IT_0495 + conj(IT_0364)*IT_0498 + IT_0364*IT_0501 + conj(IT_0356)*IT_0504
       + IT_0356*IT_0507)*IT_0513 + (IT_0366*IT_0486 + IT_0367*IT_0489 + IT_0358
      *IT_0492 + IT_0359*IT_0495 + conj(IT_0353)*IT_0498 + IT_0353*IT_0501 +
       conj(IT_0361)*IT_0504 + IT_0361*IT_0507)*IT_0517 + (conj(IT_0361)*IT_0486
       + IT_0361*IT_0489 + conj(IT_0353)*IT_0492 + IT_0353*IT_0495 + IT_0358
      *IT_0498 + IT_0359*IT_0501 + IT_0366*IT_0504 + IT_0367*IT_0507)*IT_0521 +
       4*IT_0528*(IT_0458*IT_0459 + IT_0470*(IT_0469 + 0.25*IT_0529)) + 6
      *IT_0000*(IT_0459*IT_0470 + IT_0474*(conj(IT_0457) + 0.166666666666667
      *IT_0530)) + (IT_0347*IT_0457 + IT_0344*conj(IT_0457) + IT_0266*IT_0470 +
       IT_0269*IT_0474 + conj(IT_0303)*IT_0484 + IT_0303*IT_0487 + conj(IT_0137)
      *IT_0496 + IT_0137*IT_0499)*IT_0531 + (conj(IT_0197)*IT_0198 + conj
      (IT_0264)*IT_0306 + IT_0264*IT_0308 + conj(IT_0137)*IT_0326 + IT_0137
      *IT_0328 + conj(IT_0303)*IT_0342)*IT_0532 + IT_0533*IT_0534 + IT_0536*
      (conj(IT_0338)*IT_0370 + conj(IT_0334)*IT_0378 + IT_0334*IT_0380 + conj
      (IT_0341)*IT_0473 + IT_0341*IT_0476 + conj(IT_0323)*IT_0537) + IT_0541
      *IT_0545 + IT_0547*IT_0578 + IT_0579*IT_0580 + 4*IT_0585*(conj(IT_0353)
      *IT_0587 + -IT_0588*IT_0590 + -IT_0591*IT_0592 + -IT_0366*(IT_0365 +
       IT_0586 + (-0.25)*IT_0593) + -IT_0594*(IT_0359 + IT_0595)) + 8*((conj
      (IT_0353) + -IT_0366)*IT_0591 + -IT_0587*IT_0592 + -IT_0588*IT_0594 + 
      -IT_0590*(IT_0359 + IT_0595))*IT_0597;
    return create_ccomplex_return(IT_0598);
}

