#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = pow(m_C_2, 2);
    const ccomplex_t IT_0003 = cpow(s_12 + IT_0002 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0004 = IT_0000*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_W*e_em*(IT_0008*IT_0009 + 
      -IT_0010*IT_0011)*IT_0013;
    const ccomplex_t IT_0015 = V_u2*U_Wm2;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = U_d2*V_Wp2;
    const ccomplex_t IT_0018 = IT_0011*IT_0017;
    const ccomplex_t IT_0019 = IT_0013*(IT_0016 + -IT_0018);
    const ccomplex_t IT_0020 = 1.4142135623731*e_em*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_W*e_em*(IT_0008*IT_0010 +
       IT_0009*IT_0011)*IT_0013;
    const ccomplex_t IT_0028 = IT_0008*IT_0017;
    const ccomplex_t IT_0029 = IT_0011*IT_0015;
    const ccomplex_t IT_0030 = IT_0013*(IT_0028 + IT_0029);
    const ccomplex_t IT_0031 = 1.4142135623731*e_em*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + (-2)*IT_0002 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cos(theta_W);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*IT_0013*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = m_C_2*IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0038, -1);
    const ccomplex_t IT_0043 = IT_0012*IT_0042;
    const ccomplex_t IT_0044 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = IT_0013*IT_0038;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0045 + -IT_0047 + (-2)
      *IT_0049);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0041*IT_0054;
    const ccomplex_t IT_0056 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0057 = IT_0043*IT_0056;
    const ccomplex_t IT_0058 = IT_0046*IT_0056;
    const ccomplex_t IT_0059 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0060 = IT_0046*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + -IT_0058 + (-2)
      *IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0052*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0041*IT_0064;
    const ccomplex_t IT_0066 = IT_0026 + IT_0037 + -IT_0055 + -IT_0065;
    const ccomplex_t IT_0067 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0068 = IT_0008*IT_0067;
    const ccomplex_t IT_0069 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0070 = IT_0011*IT_0069;
    const ccomplex_t IT_0071 = IT_0013*(IT_0068 + -IT_0070);
    const ccomplex_t IT_0072 = 1.4142135623731*e_em*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0014*IT_0074;
    const ccomplex_t IT_0076 = IT_0024*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0008*IT_0069;
    const ccomplex_t IT_0079 = IT_0011*IT_0067;
    const ccomplex_t IT_0080 = IT_0013*(IT_0078 + IT_0079);
    const ccomplex_t IT_0081 = 1.4142135623731*e_em*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0027*IT_0083;
    const ccomplex_t IT_0085 = IT_0035*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0055 + IT_0065 + IT_0077 + IT_0086;
    const ccomplex_t IT_0088 = pow(m_W, -4);
    const ccomplex_t IT_0089 = pow(s_34, 2);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (s_13 + -s_23)*(1 + 0.5*IT_0090);
    const ccomplex_t IT_0092 = m_C_2*IT_0091;
    const ccomplex_t IT_0093 = 4*IT_0092;
    const ccomplex_t IT_0094 = pow(m_W, 2);
    const ccomplex_t IT_0095 = s_12*IT_0094;
    const ccomplex_t IT_0096 = s_13*s_23;
    const ccomplex_t IT_0097 = IT_0002*IT_0094;
    const ccomplex_t IT_0098 = (1 + 0.5*IT_0090)*(IT_0095 + (-2)*IT_0096 +
       IT_0097);
    const ccomplex_t IT_0099 = IT_0000*IT_0005;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = pow(m_W, -2);
    const ccomplex_t IT_0102 = s_14*s_23*s_34*IT_0101;
    const ccomplex_t IT_0103 = s_13*s_24*s_34*IT_0101;
    const ccomplex_t IT_0104 = s_13*s_23*IT_0088*IT_0089;
    const ccomplex_t IT_0105 = (-2)*IT_0104;
    const ccomplex_t IT_0106 = IT_0102 + IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = -IT_0000;
    const ccomplex_t IT_0108 = IT_0005*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = s_14*s_24*IT_0101;
    const ccomplex_t IT_0111 = (-2)*IT_0096;
    const ccomplex_t IT_0112 = IT_0095 + IT_0111;
    const ccomplex_t IT_0113 = IT_0101*IT_0112;
    const ccomplex_t IT_0114 = (s_12 + -IT_0110)*(s_12 + 1.5*IT_0002 + (-0.5)
      *IT_0113);
    const ccomplex_t IT_0115 = s_14*s_34*IT_0101;
    const ccomplex_t IT_0116 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0113)*(s_13 + 
      -IT_0115);
    const ccomplex_t IT_0117 = s_24*s_34*IT_0101;
    const ccomplex_t IT_0118 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0113)*(s_23 + 
      -IT_0117);
    const ccomplex_t IT_0119 = pow(s_23, 2);
    const ccomplex_t IT_0120 = IT_0101*IT_0119;
    const ccomplex_t IT_0121 = pow(s_24, 2);
    const ccomplex_t IT_0122 = IT_0101*IT_0121;
    const ccomplex_t IT_0123 = s_13*s_24*s_34*IT_0088;
    const ccomplex_t IT_0124 = IT_0120 + IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = s_13*s_23*IT_0101;
    const ccomplex_t IT_0126 = s_23*s_24*s_34*IT_0088;
    const ccomplex_t IT_0127 = -IT_0002 + -IT_0110 + -IT_0125 + -IT_0126;
    const ccomplex_t IT_0128 = s_12 + IT_0124 + IT_0127;
    const ccomplex_t IT_0129 = m_C_2*IT_0128;
    const ccomplex_t IT_0130 = 2*IT_0066;
    const ccomplex_t IT_0131 = 2*IT_0087;
    const ccomplex_t IT_0132 = (-2)*conj(IT_0066);
    const ccomplex_t IT_0133 = (-2)*conj(IT_0087);
    const ccomplex_t IT_0134 = IT_0130 + IT_0131 + IT_0132 + IT_0133;
    const ccomplex_t IT_0135 = pow(s_13, 2);
    const ccomplex_t IT_0136 = IT_0101*IT_0135;
    const ccomplex_t IT_0137 = pow(s_14, 2);
    const ccomplex_t IT_0138 = IT_0101*IT_0137;
    const ccomplex_t IT_0139 = s_14*s_23*s_34*IT_0088;
    const ccomplex_t IT_0140 = IT_0136 + IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = s_13*s_14*s_34*IT_0088;
    const ccomplex_t IT_0142 = -IT_0002 + -IT_0110 + -IT_0125 + -IT_0141;
    const ccomplex_t IT_0143 = s_12 + IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = m_C_2*IT_0143;
    const ccomplex_t IT_0145 = 2*conj(IT_0066);
    const ccomplex_t IT_0146 = 2*conj(IT_0087);
    const ccomplex_t IT_0147 = (-2)*IT_0066;
    const ccomplex_t IT_0148 = (-2)*IT_0087;
    const ccomplex_t IT_0149 = IT_0145 + IT_0146 + IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = IT_0089*IT_0101;
    const ccomplex_t IT_0151 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0113)*(IT_0094 
      + -IT_0150);
    const ccomplex_t IT_0152 = s_34*m_C_2*IT_0101;
    const ccomplex_t IT_0153 = s_23*s_34*IT_0101;
    const ccomplex_t IT_0154 = s_13*s_34*IT_0101;
    const ccomplex_t IT_0155 = -s_24 + -IT_0154;
    const ccomplex_t IT_0156 = s_14 + IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = IT_0152*IT_0156;
    const ccomplex_t IT_0158 = 2*IT_0000;
    const ccomplex_t IT_0159 = IT_0005*IT_0158;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = s_14*s_24;
    const ccomplex_t IT_0162 = (-2)*IT_0161;
    const ccomplex_t IT_0163 = IT_0095 + IT_0162;
    const ccomplex_t IT_0164 = IT_0101*IT_0163;
    const ccomplex_t IT_0165 = (s_12 + -IT_0125)*(s_12 + 1.5*IT_0002 + (-0.5)
      *IT_0164);
    const ccomplex_t IT_0166 = IT_0120 + IT_0122 + IT_0139;
    const ccomplex_t IT_0167 = s_12 + IT_0127 + IT_0166;
    const ccomplex_t IT_0168 = m_C_2*IT_0167;
    const ccomplex_t IT_0169 = IT_0123 + IT_0136 + IT_0138;
    const ccomplex_t IT_0170 = s_12 + IT_0142 + IT_0169;
    const ccomplex_t IT_0171 = m_C_2*IT_0170;
    const ccomplex_t IT_0172 = (s_12 + (-0.5)*IT_0113)*(s_13 + -IT_0115);
    const ccomplex_t IT_0173 = (-2)*IT_0172;
    const ccomplex_t IT_0174 = (s_12 + (-0.5)*IT_0113)*(s_23 + -IT_0117);
    const ccomplex_t IT_0175 = (-2)*IT_0174;
    const ccomplex_t IT_0176 = IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = IT_0039*IT_0054;
    const ccomplex_t IT_0178 = IT_0040*IT_0054;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = 2*IT_0179;
    const ccomplex_t IT_0181 = conj(IT_0177)*IT_0180;
    const ccomplex_t IT_0182 = 2*conj(IT_0179);
    const ccomplex_t IT_0183 = IT_0177*IT_0182;
    const ccomplex_t IT_0184 = s_13*s_23*IT_0088*IT_0137;
    const ccomplex_t IT_0185 = s_14*s_24*IT_0088*IT_0135;
    const ccomplex_t IT_0186 = s_12*IT_0002;
    const ccomplex_t IT_0187 = s_14*s_24*IT_0002*IT_0101;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = s_13*s_23*IT_0002*IT_0101;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = s_12*s_13*s_14*s_34*IT_0088;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = IT_0184 + IT_0185 + IT_0186 + IT_0188 + IT_0190
       + IT_0192;
    const ccomplex_t IT_0194 = s_13*s_14*s_23*s_24*IT_0088;
    const ccomplex_t IT_0195 = pow(m_C_2, 4);
    const ccomplex_t IT_0196 = IT_0088*IT_0121*IT_0135;
    const ccomplex_t IT_0197 = IT_0002*IT_0101*IT_0121;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0002*IT_0101*IT_0135;
    const ccomplex_t IT_0200 = -IT_0199;
    const ccomplex_t IT_0201 = s_12*s_13*s_24*s_34*IT_0088;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0194 + IT_0195 + IT_0196 + IT_0198 + IT_0200
       + IT_0202;
    const ccomplex_t IT_0204 = IT_0193 + IT_0203;
    const ccomplex_t IT_0205 = 2*IT_0039;
    const ccomplex_t IT_0206 = IT_0064*IT_0205;
    const ccomplex_t IT_0207 = IT_0040*IT_0064;
    const ccomplex_t IT_0208 = 2*IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0206)*IT_0208;
    const ccomplex_t IT_0210 = 2*conj(IT_0207);
    const ccomplex_t IT_0211 = IT_0206*IT_0210;
    const ccomplex_t IT_0212 = s_14*s_24*IT_0088*IT_0119;
    const ccomplex_t IT_0213 = s_13*s_23*IT_0088*IT_0121;
    const ccomplex_t IT_0214 = s_12*s_23*s_24*s_34*IT_0088;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0186 + IT_0188 + IT_0190 + IT_0212 + IT_0213
       + IT_0215;
    const ccomplex_t IT_0217 = IT_0088*IT_0119*IT_0137;
    const ccomplex_t IT_0218 = IT_0002*IT_0101*IT_0119;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0002*IT_0101*IT_0137;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = s_12*s_14*s_23*s_34*IT_0088;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = IT_0194 + IT_0195 + IT_0217 + IT_0219 + IT_0221
       + IT_0223;
    const ccomplex_t IT_0225 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0226 = IT_0013*IT_0225;
    const ccomplex_t IT_0227 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0228 = IT_0013*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*(IT_0226 + 0.5*IT_0229);
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0233 = IT_0013*IT_0232;
    const ccomplex_t IT_0234 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0235 = IT_0013*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0233 + 0.5*IT_0236);
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = IT_0231*IT_0238;
    const ccomplex_t IT_0240 = cpow((-2)*s_23 + IT_0002 + IT_0094 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0241 = IT_0239*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0244 = IT_0013*IT_0243;
    const ccomplex_t IT_0245 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0246 = IT_0013*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0244 + 0.5*IT_0247);
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0251 = IT_0013*IT_0250;
    const ccomplex_t IT_0252 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0253 = IT_0013*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*(IT_0251 + 0.5*IT_0254);
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = IT_0249*IT_0256;
    const ccomplex_t IT_0258 = cpow((-2)*s_23 + IT_0002 + IT_0094 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0259 = IT_0257*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = cpow((-2)*s_23 + IT_0002 + IT_0094 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0262 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0263 = IT_0013*IT_0262;
    const ccomplex_t IT_0264 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0265 = IT_0013*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0263 + 0.5*IT_0266);
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0270 = IT_0013*IT_0269;
    const ccomplex_t IT_0271 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0272 = IT_0013*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*(IT_0270 + 0.5*IT_0273);
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = IT_0268*IT_0275;
    const ccomplex_t IT_0277 = IT_0261*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0280 = IT_0013*IT_0279;
    const ccomplex_t IT_0281 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0282 = IT_0013*IT_0281;
    const ccomplex_t IT_0283 = 1.4142135623731*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*(IT_0280 + 0.5*IT_0283);
    const ccomplex_t IT_0285 = -IT_0284;
    const ccomplex_t IT_0286 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0287 = IT_0013*IT_0286;
    const ccomplex_t IT_0288 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0289 = IT_0013*IT_0288;
    const ccomplex_t IT_0290 = 1.4142135623731*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*(IT_0287 + 0.5*IT_0290);
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = IT_0285*IT_0292;
    const ccomplex_t IT_0294 = cpow((-2)*s_23 + IT_0002 + IT_0094 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0295 = IT_0293*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0054*IT_0205;
    const ccomplex_t IT_0298 = (-2)*IT_0242 + (-2)*IT_0260 + (-2)*IT_0278 + (
      -2)*IT_0296 + -IT_0297;
    const ccomplex_t IT_0299 = IT_0180*conj(IT_0298);
    const ccomplex_t IT_0300 = IT_0182*IT_0298;
    const ccomplex_t IT_0301 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0302 = IT_0013*IT_0301;
    const ccomplex_t IT_0303 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0304 = IT_0013*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*(IT_0302 + (-0.5)*IT_0305);
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0309 = IT_0013*IT_0308;
    const ccomplex_t IT_0310 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0311 = IT_0013*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0309 + (-0.5)*IT_0312);
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = IT_0307*IT_0314;
    const ccomplex_t IT_0316 = IT_0240*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0319 = IT_0013*IT_0318;
    const ccomplex_t IT_0320 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0321 = IT_0013*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*(IT_0319 + (-0.5)*IT_0322);
    const ccomplex_t IT_0324 = -IT_0323;
    const ccomplex_t IT_0325 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0326 = IT_0013*IT_0325;
    const ccomplex_t IT_0327 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0328 = IT_0013*IT_0327;
    const ccomplex_t IT_0329 = 1.4142135623731*IT_0328;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*(IT_0326 + (-0.5)*IT_0329);
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0324*IT_0331;
    const ccomplex_t IT_0333 = IT_0261*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*IT_0333;
    const ccomplex_t IT_0335 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0336 = IT_0013*IT_0335;
    const ccomplex_t IT_0337 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0338 = IT_0013*IT_0337;
    const ccomplex_t IT_0339 = 1.4142135623731*IT_0338;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*(IT_0336 + (-0.5)*IT_0339);
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0343 = IT_0013*IT_0342;
    const ccomplex_t IT_0344 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0345 = IT_0013*IT_0344;
    const ccomplex_t IT_0346 = 1.4142135623731*IT_0345;
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*(IT_0343 + (-0.5)*IT_0346);
    const ccomplex_t IT_0348 = -IT_0347;
    const ccomplex_t IT_0349 = IT_0341*IT_0348;
    const ccomplex_t IT_0350 = IT_0294*IT_0349;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*IT_0350;
    const ccomplex_t IT_0352 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0353 = IT_0013*IT_0352;
    const ccomplex_t IT_0354 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0355 = IT_0013*IT_0354;
    const ccomplex_t IT_0356 = 1.4142135623731*IT_0355;
    const ccomplex_t IT_0357 = (0 + _Complex_I*1)*(IT_0353 + (-0.5)*IT_0356);
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0360 = IT_0013*IT_0359;
    const ccomplex_t IT_0361 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0362 = IT_0013*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*(IT_0360 + (-0.5)*IT_0363);
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = IT_0358*IT_0365;
    const ccomplex_t IT_0367 = IT_0258*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = IT_0206 + (-2)*IT_0317 + (-2)*IT_0334 + (-2)
      *IT_0351 + (-2)*IT_0368;
    const ccomplex_t IT_0370 = IT_0208*conj(IT_0369);
    const ccomplex_t IT_0371 = IT_0210*IT_0369;
    const ccomplex_t IT_0372 = s_13*IT_0094;
    const ccomplex_t IT_0373 = s_14*s_34;
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = IT_0101*IT_0375;
    const ccomplex_t IT_0377 = -IT_0115 + (-0.25)*IT_0376;
    const ccomplex_t IT_0378 = s_13 + IT_0377;
    const ccomplex_t IT_0379 = IT_0002*IT_0378;
    const ccomplex_t IT_0380 = (-4)*IT_0379;
    const ccomplex_t IT_0381 = s_23*IT_0094;
    const ccomplex_t IT_0382 = s_24*s_34;
    const ccomplex_t IT_0383 = -IT_0382;
    const ccomplex_t IT_0384 = IT_0381 + IT_0383;
    const ccomplex_t IT_0385 = IT_0101*IT_0384;
    const ccomplex_t IT_0386 = -IT_0117 + (-0.25)*IT_0385;
    const ccomplex_t IT_0387 = s_23 + IT_0386;
    const ccomplex_t IT_0388 = IT_0002*IT_0387;
    const ccomplex_t IT_0389 = (-4)*IT_0388;
    const ccomplex_t IT_0390 = IT_0317 + IT_0334 + IT_0351 + IT_0368;
    const ccomplex_t IT_0391 = IT_0180*conj(IT_0390);
    const ccomplex_t IT_0392 = IT_0182*IT_0390;
    const ccomplex_t IT_0393 = IT_0242 + IT_0260 + IT_0278 + IT_0296;
    const ccomplex_t IT_0394 = IT_0208*conj(IT_0393);
    const ccomplex_t IT_0395 = IT_0210*IT_0393;
    const ccomplex_t IT_0396 = 0.5*IT_0090;
    const ccomplex_t IT_0397 = 1 + IT_0396;
    const ccomplex_t IT_0398 = s_12*IT_0397;
    const ccomplex_t IT_0399 = 4*IT_0087;
    const ccomplex_t IT_0400 = -IT_0138;
    const ccomplex_t IT_0401 = IT_0002 + IT_0400;
    const ccomplex_t IT_0402 = IT_0002*IT_0401;
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = -IT_0122;
    const ccomplex_t IT_0405 = IT_0002 + IT_0404;
    const ccomplex_t IT_0406 = IT_0002*IT_0405;
    const ccomplex_t IT_0407 = -IT_0406;
    const ccomplex_t IT_0408 = IT_0179*conj(IT_0207);
    const ccomplex_t IT_0409 = 6*IT_0408;
    const ccomplex_t IT_0410 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0113)*(IT_0002 
      + -IT_0138);
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0113)*(IT_0002 
      + -IT_0122);
    const ccomplex_t IT_0413 = -IT_0412;
    const ccomplex_t IT_0414 = cpow(IT_0109, 2);
    const ccomplex_t IT_0415 = (-8)*IT_0414;
    const ccomplex_t IT_0416 = (IT_0002 + -IT_0136)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0164);
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = (IT_0002 + -IT_0120)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0164);
    const ccomplex_t IT_0419 = -IT_0418;
    const ccomplex_t IT_0420 = cpow(IT_0160, 2);
    const ccomplex_t IT_0421 = (-8)*IT_0420;
    const ccomplex_t IT_0422 = -IT_0115;
    const ccomplex_t IT_0423 = s_13 + IT_0422;
    const ccomplex_t IT_0424 = IT_0002*IT_0423;
    const ccomplex_t IT_0425 = -IT_0117;
    const ccomplex_t IT_0426 = s_23 + IT_0425;
    const ccomplex_t IT_0427 = IT_0002*IT_0426;
    const ccomplex_t IT_0428 = IT_0424 + IT_0427;
    const ccomplex_t IT_0429 = IT_0039*IT_0064;
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = (-6)*conj(IT_0430);
    const ccomplex_t IT_0432 = IT_0179*IT_0431;
    const ccomplex_t IT_0433 = conj(IT_0177)*IT_0207;
    const ccomplex_t IT_0434 = (-6)*IT_0433;
    const ccomplex_t IT_0435 = IT_0177*conj(IT_0207);
    const ccomplex_t IT_0436 = (-6)*IT_0435;
    const ccomplex_t IT_0437 = IT_0002*IT_0397;
    const ccomplex_t IT_0438 = (-4)*IT_0066;
    const ccomplex_t IT_0439 = (-4)*conj(IT_0066);
    const ccomplex_t IT_0440 = -IT_0002;
    const ccomplex_t IT_0441 = IT_0138 + IT_0440;
    const ccomplex_t IT_0442 = (-2)*IT_0096 + 2*IT_0102 + -IT_0161;
    const ccomplex_t IT_0443 = IT_0095 + IT_0442;
    const ccomplex_t IT_0444 = IT_0101*IT_0443;
    const ccomplex_t IT_0445 = -IT_0097;
    const ccomplex_t IT_0446 = IT_0137 + IT_0445;
    const ccomplex_t IT_0447 = IT_0101*IT_0446;
    const ccomplex_t IT_0448 = -IT_0110 + (-2)*IT_0441 + (-0.5)*IT_0444 + 0.5
      *IT_0447;
    const ccomplex_t IT_0449 = s_12 + IT_0448;
    const ccomplex_t IT_0450 = m_C_2*IT_0449;
    const ccomplex_t IT_0451 = (-2)*IT_0450;
    const ccomplex_t IT_0452 = s_23*s_24*s_34*IT_0101;
    const ccomplex_t IT_0453 = (-0.5)*IT_0097 + 0.5*IT_0121 + -IT_0452;
    const ccomplex_t IT_0454 = IT_0119 + IT_0453;
    const ccomplex_t IT_0455 = IT_0101*IT_0454;
    const ccomplex_t IT_0456 = IT_0002 + IT_0455;
    const ccomplex_t IT_0457 = -IT_0161;
    const ccomplex_t IT_0458 = IT_0095 + IT_0457;
    const ccomplex_t IT_0459 = IT_0101*IT_0458;
    const ccomplex_t IT_0460 = -IT_0110;
    const ccomplex_t IT_0461 = s_12 + IT_0460;
    const ccomplex_t IT_0462 = -IT_0122 + (-0.5)*IT_0459 + 2*IT_0461;
    const ccomplex_t IT_0463 = IT_0456 + IT_0462;
    const ccomplex_t IT_0464 = m_C_2*IT_0463;
    const ccomplex_t IT_0465 = (-2)*IT_0464;
    const ccomplex_t IT_0466 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0467 = IT_0268*IT_0324;
    const ccomplex_t IT_0468 = IT_0261*IT_0466*IT_0467;
    const ccomplex_t IT_0469 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0470 = IT_0292*IT_0341;
    const ccomplex_t IT_0471 = IT_0294*IT_0469*IT_0470;
    const ccomplex_t IT_0472 = m_C_2*IT_0334;
    const ccomplex_t IT_0473 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0474 = IT_0238*IT_0307;
    const ccomplex_t IT_0475 = IT_0240*IT_0473*IT_0474;
    const ccomplex_t IT_0476 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0477 = IT_0256*IT_0365;
    const ccomplex_t IT_0478 = IT_0258*IT_0476*IT_0477;
    const ccomplex_t IT_0479 = m_C_2*IT_0351;
    const ccomplex_t IT_0480 = m_C_2*IT_0317;
    const ccomplex_t IT_0481 = m_C_2*IT_0368;
    const ccomplex_t IT_0482 = IT_0468 + IT_0471 + -IT_0472 + IT_0475 +
       IT_0478 + -IT_0479 + -IT_0480 + -IT_0481;
    const ccomplex_t IT_0483 = 2*IT_0109;
    const ccomplex_t IT_0484 = conj(IT_0482)*IT_0483;
    const ccomplex_t IT_0485 = IT_0275*IT_0331;
    const ccomplex_t IT_0486 = IT_0261*IT_0466*IT_0485;
    const ccomplex_t IT_0487 = IT_0285*IT_0348;
    const ccomplex_t IT_0488 = IT_0294*IT_0469*IT_0487;
    const ccomplex_t IT_0489 = IT_0231*IT_0314;
    const ccomplex_t IT_0490 = IT_0240*IT_0473*IT_0489;
    const ccomplex_t IT_0491 = IT_0249*IT_0358;
    const ccomplex_t IT_0492 = IT_0258*IT_0476*IT_0491;
    const ccomplex_t IT_0493 = m_C_2*IT_0278;
    const ccomplex_t IT_0494 = m_C_2*IT_0296;
    const ccomplex_t IT_0495 = m_C_2*IT_0242;
    const ccomplex_t IT_0496 = m_C_2*IT_0260;
    const ccomplex_t IT_0497 = IT_0486 + IT_0488 + IT_0490 + IT_0492 + 
      -IT_0493 + -IT_0494 + -IT_0495 + -IT_0496;
    const ccomplex_t IT_0498 = IT_0483*conj(IT_0497);
    const ccomplex_t IT_0499 = (-2)*IT_0109;
    const ccomplex_t IT_0500 = IT_0482*IT_0499;
    const ccomplex_t IT_0501 = IT_0497*IT_0499;
    const ccomplex_t IT_0502 = 2*IT_0410;
    const ccomplex_t IT_0503 = 2*IT_0412;
    const ccomplex_t IT_0504 = IT_0179*IT_0483;
    const ccomplex_t IT_0505 = IT_0207*IT_0483;
    const ccomplex_t IT_0506 = conj(IT_0179)*IT_0499;
    const ccomplex_t IT_0507 = conj(IT_0207)*IT_0499;
    const ccomplex_t IT_0508 = (-2)*IT_0116;
    const ccomplex_t IT_0509 = (-2)*IT_0118;
    const ccomplex_t IT_0510 = IT_0508 + IT_0509;
    const ccomplex_t IT_0511 = conj(IT_0177)*IT_0483;
    const ccomplex_t IT_0512 = 2*IT_0100;
    const ccomplex_t IT_0513 = conj(IT_0207)*IT_0512;
    const ccomplex_t IT_0514 = IT_0177*IT_0499;
    const ccomplex_t IT_0515 = (-2)*IT_0100;
    const ccomplex_t IT_0516 = IT_0207*IT_0515;
    const ccomplex_t IT_0517 = -IT_0096;
    const ccomplex_t IT_0518 = IT_0095 + IT_0517;
    const ccomplex_t IT_0519 = IT_0101*IT_0518;
    const ccomplex_t IT_0520 = s_13*s_14*s_34*IT_0101;
    const ccomplex_t IT_0521 = -IT_0097 + 2*IT_0137 + (-2)*IT_0520;
    const ccomplex_t IT_0522 = IT_0135 + IT_0521;
    const ccomplex_t IT_0523 = IT_0101*IT_0522;
    const ccomplex_t IT_0524 = -IT_0125;
    const ccomplex_t IT_0525 = s_12 + IT_0524;
    const ccomplex_t IT_0526 = -IT_0136 + (-0.5)*IT_0519 + 0.5*IT_0523 + 2
      *IT_0525;
    const ccomplex_t IT_0527 = IT_0002 + IT_0526;
    const ccomplex_t IT_0528 = m_C_2*IT_0527;
    const ccomplex_t IT_0529 = (-2)*IT_0528;
    const ccomplex_t IT_0530 = -IT_0096 + 2*IT_0102 + (-2)*IT_0161;
    const ccomplex_t IT_0531 = IT_0095 + IT_0530;
    const ccomplex_t IT_0532 = IT_0101*IT_0531;
    const ccomplex_t IT_0533 = IT_0119 + IT_0445;
    const ccomplex_t IT_0534 = IT_0101*IT_0533;
    const ccomplex_t IT_0535 = IT_0120 + IT_0440;
    const ccomplex_t IT_0536 = -IT_0125 + (-0.5)*IT_0532 + 0.5*IT_0534 + (-2)
      *IT_0535;
    const ccomplex_t IT_0537 = s_12 + IT_0536;
    const ccomplex_t IT_0538 = m_C_2*IT_0537;
    const ccomplex_t IT_0539 = (-2)*IT_0538;
    const ccomplex_t IT_0540 = 2*IT_0160;
    const ccomplex_t IT_0541 = IT_0482*IT_0540;
    const ccomplex_t IT_0542 = IT_0497*IT_0540;
    const ccomplex_t IT_0543 = (-2)*IT_0160;
    const ccomplex_t IT_0544 = conj(IT_0482)*IT_0543;
    const ccomplex_t IT_0545 = conj(IT_0497)*IT_0543;
    const ccomplex_t IT_0546 = IT_0002 + IT_0141;
    const ccomplex_t IT_0547 = -IT_0136 + -IT_0138;
    const ccomplex_t IT_0548 = IT_0546 + IT_0547;
    const ccomplex_t IT_0549 = IT_0002*IT_0548;
    const ccomplex_t IT_0550 = -IT_0549;
    const ccomplex_t IT_0551 = IT_0193 + IT_0550;
    const ccomplex_t IT_0552 = -IT_0110 + -IT_0125;
    const ccomplex_t IT_0553 = s_12 + IT_0123 + IT_0552;
    const ccomplex_t IT_0554 = IT_0002*IT_0553;
    const ccomplex_t IT_0555 = -IT_0554;
    const ccomplex_t IT_0556 = IT_0203 + IT_0555;
    const ccomplex_t IT_0557 = IT_0551 + IT_0556;
    const ccomplex_t IT_0558 = -IT_0297;
    const ccomplex_t IT_0559 = IT_0483*conj(IT_0558);
    const ccomplex_t IT_0560 = conj(IT_0179)*IT_0540;
    const ccomplex_t IT_0561 = conj(IT_0207)*IT_0540;
    const ccomplex_t IT_0562 = (-8)*IT_0109;
    const ccomplex_t IT_0563 = IT_0160*IT_0562;
    const ccomplex_t IT_0564 = IT_0499*IT_0558;
    const ccomplex_t IT_0565 = IT_0179*IT_0543;
    const ccomplex_t IT_0566 = IT_0207*IT_0543;
    const ccomplex_t IT_0567 = IT_0002 + IT_0126;
    const ccomplex_t IT_0568 = -IT_0120 + -IT_0122;
    const ccomplex_t IT_0569 = IT_0567 + IT_0568;
    const ccomplex_t IT_0570 = IT_0002*IT_0569;
    const ccomplex_t IT_0571 = -IT_0570;
    const ccomplex_t IT_0572 = IT_0216 + IT_0571;
    const ccomplex_t IT_0573 = s_12 + IT_0139 + IT_0552;
    const ccomplex_t IT_0574 = IT_0002*IT_0573;
    const ccomplex_t IT_0575 = -IT_0574;
    const ccomplex_t IT_0576 = IT_0224 + IT_0575;
    const ccomplex_t IT_0577 = conj(IT_0369)*IT_0483;
    const ccomplex_t IT_0578 = conj(IT_0298)*IT_0483;
    const ccomplex_t IT_0579 = IT_0369*IT_0499;
    const ccomplex_t IT_0580 = IT_0298*IT_0499;
    const ccomplex_t IT_0581 = s_23*IT_0088*IT_0089;
    const ccomplex_t IT_0582 = -IT_0581;
    const ccomplex_t IT_0583 = IT_0117 + IT_0582;
    const ccomplex_t IT_0584 = IT_0002*IT_0583;
    const ccomplex_t IT_0585 = s_13*IT_0088*IT_0089;
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = IT_0115 + IT_0586;
    const ccomplex_t IT_0588 = IT_0002*IT_0587;
    const ccomplex_t IT_0589 = s_13*s_14*s_23*s_34*IT_0088;
    const ccomplex_t IT_0590 = s_23*IT_0002;
    const ccomplex_t IT_0591 = s_24*s_34*IT_0088*IT_0135;
    const ccomplex_t IT_0592 = s_12*s_13;
    const ccomplex_t IT_0593 = s_23*IT_0101*IT_0135;
    const ccomplex_t IT_0594 = (-2)*IT_0593;
    const ccomplex_t IT_0595 = IT_0592 + IT_0594;
    const ccomplex_t IT_0596 = s_12*s_13*IT_0088*IT_0089;
    const ccomplex_t IT_0597 = -IT_0596;
    const ccomplex_t IT_0598 = s_24*s_34*IT_0002*IT_0101;
    const ccomplex_t IT_0599 = -IT_0598;
    const ccomplex_t IT_0600 = IT_0589 + IT_0590 + IT_0591 + IT_0595 + IT_0597
       + IT_0599;
    const ccomplex_t IT_0601 = IT_0588 + IT_0600;
    const ccomplex_t IT_0602 = s_14*s_34*IT_0088*IT_0119;
    const ccomplex_t IT_0603 = s_12*s_23;
    const ccomplex_t IT_0604 = s_13*s_23*s_24*s_34*IT_0088;
    const ccomplex_t IT_0605 = s_13*IT_0002;
    const ccomplex_t IT_0606 = IT_0603 + IT_0604 + IT_0605;
    const ccomplex_t IT_0607 = IT_0602 + IT_0606;
    const ccomplex_t IT_0608 = s_13*IT_0101*IT_0119;
    const ccomplex_t IT_0609 = (-2)*IT_0608;
    const ccomplex_t IT_0610 = s_12*s_23*IT_0088*IT_0089;
    const ccomplex_t IT_0611 = -IT_0610;
    const ccomplex_t IT_0612 = s_14*s_34*IT_0002*IT_0101;
    const ccomplex_t IT_0613 = -IT_0612;
    const ccomplex_t IT_0614 = IT_0584 + IT_0601 + IT_0607 + IT_0609 + IT_0611
       + IT_0613;
    const ccomplex_t IT_0615 = conj(IT_0177)*IT_0540;
    const ccomplex_t IT_0616 = IT_0100*IT_0160;
    const ccomplex_t IT_0617 = (-8)*IT_0616;
    const ccomplex_t IT_0618 = IT_0177*IT_0543;
    const ccomplex_t IT_0619 = s_12*s_14*s_34*IT_0101;
    const ccomplex_t IT_0620 = s_13*s_14*s_24*IT_0101;
    const ccomplex_t IT_0621 = -IT_0620;
    const ccomplex_t IT_0622 = s_23*IT_0101*IT_0137;
    const ccomplex_t IT_0623 = -IT_0622;
    const ccomplex_t IT_0624 = IT_0595 + IT_0621 + IT_0623;
    const ccomplex_t IT_0625 = IT_0589 + IT_0590 + IT_0591 + IT_0597 + IT_0619
       + IT_0624;
    const ccomplex_t IT_0626 = s_12*s_24*s_34*IT_0101;
    const ccomplex_t IT_0627 = IT_0602 + IT_0606 + IT_0626;
    const ccomplex_t IT_0628 = s_13*IT_0101*IT_0121;
    const ccomplex_t IT_0629 = -IT_0628;
    const ccomplex_t IT_0630 = s_14*s_23*s_24*IT_0101;
    const ccomplex_t IT_0631 = -IT_0630;
    const ccomplex_t IT_0632 = (-2)*IT_0039;
    const ccomplex_t IT_0633 = IT_0064*IT_0632;
    const ccomplex_t IT_0634 = -IT_0633;
    const ccomplex_t IT_0635 = IT_0540*conj(IT_0634);
    const ccomplex_t IT_0636 = IT_0054*IT_0632;
    const ccomplex_t IT_0637 = IT_0540*conj(IT_0636);
    const ccomplex_t IT_0638 = IT_0007*IT_0160;
    const ccomplex_t IT_0639 = (-8)*IT_0638;
    const ccomplex_t IT_0640 = IT_0543*IT_0634;
    const ccomplex_t IT_0641 = IT_0543*IT_0636;
    const ccomplex_t IT_0642 = s_34*IT_0002*IT_0101;
    const ccomplex_t IT_0643 = -IT_0154;
    const ccomplex_t IT_0644 = s_14 + IT_0643;
    const ccomplex_t IT_0645 = IT_0642*IT_0644;
    const ccomplex_t IT_0646 = 2*IT_0558;
    const ccomplex_t IT_0647 = 2*IT_0206;
    const ccomplex_t IT_0648 = IT_0609 + IT_0611 + IT_0627 + IT_0629 + IT_0631;
    const ccomplex_t IT_0649 = 2*IT_0369;
    const ccomplex_t IT_0650 = -IT_0153;
    const ccomplex_t IT_0651 = s_24 + IT_0650;
    const ccomplex_t IT_0652 = IT_0642*IT_0651;
    const ccomplex_t IT_0653 = 2*IT_0298;
    const ccomplex_t IT_0654 = s_23*m_C_2;
    const ccomplex_t IT_0655 = IT_0397*IT_0654;
    const ccomplex_t IT_0656 = (1 + 0.5*IT_0090)*(IT_0095 + (-2)*IT_0096);
    const ccomplex_t IT_0657 = (-4)*IT_0634;
    const ccomplex_t IT_0658 = s_13*m_C_2;
    const ccomplex_t IT_0659 = IT_0397*IT_0658;
    const ccomplex_t IT_0660 = (-4)*IT_0007;
    const ccomplex_t IT_0661 = s_34*IT_0101;
    const ccomplex_t IT_0662 = s_14*s_23;
    const ccomplex_t IT_0663 = s_13*s_24;
    const ccomplex_t IT_0664 = IT_0662 + IT_0663;
    const ccomplex_t IT_0665 = s_13*s_23*s_34*IT_0101;
    const ccomplex_t IT_0666 = (-2)*IT_0665;
    const ccomplex_t IT_0667 = IT_0664 + IT_0666;
    const ccomplex_t IT_0668 = IT_0661*IT_0667;
    const ccomplex_t IT_0669 = (-2)*IT_0430;
    const ccomplex_t IT_0670 = conj(IT_0634)*(IT_0106*IT_0515 + IT_0645
      *IT_0646 + IT_0625*IT_0647 + IT_0648*IT_0649 + IT_0652*IT_0653 + IT_0399
      *IT_0655 + IT_0656*IT_0657 + IT_0438*IT_0659 + IT_0098*IT_0660 + IT_0668
      *IT_0669);
    const ccomplex_t IT_0671 = m_C_2*IT_0548;
    const ccomplex_t IT_0672 = -IT_0671;
    const ccomplex_t IT_0673 = m_C_2*IT_0573;
    const ccomplex_t IT_0674 = (-2)*IT_0114;
    const ccomplex_t IT_0675 = 4*IT_0109;
    const ccomplex_t IT_0676 = (s_12 + -IT_0110)*(s_12 + (-0.5)*IT_0113);
    const ccomplex_t IT_0677 = (-2)*IT_0676;
    const ccomplex_t IT_0678 = 4*IT_0207;
    const ccomplex_t IT_0679 = m_C_2*IT_0553;
    const ccomplex_t IT_0680 = m_C_2*IT_0569;
    const ccomplex_t IT_0681 = -IT_0680;
    const ccomplex_t IT_0682 = 2*IT_0589;
    const ccomplex_t IT_0683 = IT_0590 + IT_0624 + IT_0682;
    const ccomplex_t IT_0684 = 2*IT_0604;
    const ccomplex_t IT_0685 = IT_0603 + IT_0605 + IT_0609 + IT_0629 + IT_0631
       + IT_0684;
    const ccomplex_t IT_0686 = IT_0683 + IT_0685;
    const ccomplex_t IT_0687 = IT_0483*conj(IT_0634);
    const ccomplex_t IT_0688 = IT_0483*conj(IT_0636);
    const ccomplex_t IT_0689 = IT_0180*conj(IT_0636);
    const ccomplex_t IT_0690 = IT_0182*IT_0636;
    const ccomplex_t IT_0691 = 2*IT_0007;
    const ccomplex_t IT_0692 = conj(IT_0207)*IT_0691;
    const ccomplex_t IT_0693 = IT_0208*conj(IT_0634);
    const ccomplex_t IT_0694 = IT_0210*IT_0634;
    const ccomplex_t IT_0695 = IT_0007*IT_0562;
    const ccomplex_t IT_0696 = IT_0499*IT_0634;
    const ccomplex_t IT_0697 = IT_0499*IT_0636;
    const ccomplex_t IT_0698 = (-2)*IT_0007;
    const ccomplex_t IT_0699 = IT_0207*IT_0698;
    const ccomplex_t IT_0700 = (s_12 + (-0.5)*IT_0113)*(IT_0002 + -IT_0138);
    const ccomplex_t IT_0701 = 2*IT_0700;
    const ccomplex_t IT_0702 = (s_12 + (-0.5)*IT_0113)*(IT_0002 + -IT_0122);
    const ccomplex_t IT_0703 = 2*IT_0702;
    const ccomplex_t IT_0704 = IT_0701 + IT_0703;
    const ccomplex_t IT_0705 = (-2)*IT_0207;
    const ccomplex_t IT_0706 = conj(IT_0207)*IT_0705;
    const ccomplex_t IT_0707 = (-2)*IT_0165;
    const ccomplex_t IT_0708 = 2*IT_0430;
    const ccomplex_t IT_0709 = IT_0588 + IT_0625;
    const ccomplex_t IT_0710 = 2*IT_0416;
    const ccomplex_t IT_0711 = (IT_0002 + -IT_0136)*(s_12 + (-0.5)*IT_0164);
    const ccomplex_t IT_0712 = 2*IT_0711;
    const ccomplex_t IT_0713 = (-2)*IT_0206;
    const ccomplex_t IT_0714 = IT_0097*IT_0397;
    const ccomplex_t IT_0715 = 2*conj(IT_0558);
    const ccomplex_t IT_0716 = 2*conj(IT_0206);
    const ccomplex_t IT_0717 = 2*conj(IT_0369);
    const ccomplex_t IT_0718 = 2*conj(IT_0298);
    const ccomplex_t IT_0719 = 4*IT_0007;
    const ccomplex_t IT_0720 = (-4)*conj(IT_0087);
    const ccomplex_t IT_0721 = (-2)*conj(IT_0177);
    const ccomplex_t IT_0722 = (-0.25)*IT_0106*IT_0512 + -conj(IT_0066)
      *IT_0655 + conj(IT_0636)*IT_0656 + conj(IT_0634)*IT_0714 + (-0.25)*IT_0625
      *IT_0715 + (-0.25)*IT_0645*IT_0716 + (-0.25)*IT_0652*IT_0717 + (-0.25)
      *IT_0648*IT_0718 + (-0.25)*IT_0098*IT_0719 + (-0.25)*IT_0659*IT_0720 + (
      -0.25)*IT_0668*IT_0721;
    const ccomplex_t IT_0723 = IT_0607 + IT_0609 + IT_0611 + IT_0613;
    const ccomplex_t IT_0724 = IT_0584 + IT_0609 + IT_0611 + IT_0627 + IT_0629
       + IT_0631;
    const ccomplex_t IT_0725 = (s_12 + -IT_0125)*(s_12 + (-0.5)*IT_0164);
    const ccomplex_t IT_0726 = (-2)*IT_0725;
    const ccomplex_t IT_0727 = 2*IT_0418;
    const ccomplex_t IT_0728 = (IT_0002 + -IT_0120)*(s_12 + (-0.5)*IT_0164);
    const ccomplex_t IT_0729 = 2*IT_0728;
    const ccomplex_t IT_0730 = (-2)*IT_0369;
    const ccomplex_t IT_0731 = IT_0512*(IT_0584 + IT_0607 + IT_0609 + IT_0611 
      + IT_0613) + IT_0147*IT_0673 + IT_0148*IT_0681 + IT_0540*IT_0707 + IT_0708
      *IT_0723 + IT_0691*IT_0724 + IT_0647*IT_0726 + IT_0543*IT_0727 + IT_0729
      *IT_0730;
    const ccomplex_t IT_0732 = IT_0095*IT_0101;
    const ccomplex_t IT_0733 = (-0.5)*IT_0732;
    const ccomplex_t IT_0734 = s_34*IT_0088;
    const ccomplex_t IT_0735 = s_12*s_34;
    const ccomplex_t IT_0736 = IT_0662 + IT_0735;
    const ccomplex_t IT_0737 = -IT_0663;
    const ccomplex_t IT_0738 = IT_0736 + IT_0737;
    const ccomplex_t IT_0739 = IT_0734*IT_0738;
    const ccomplex_t IT_0740 = 0.25*IT_0739;
    const ccomplex_t IT_0741 = s_12 + IT_0733 + IT_0740;
    const ccomplex_t IT_0742 = m_C_2*IT_0378;
    const ccomplex_t IT_0743 = (-4)*IT_0742;
    const ccomplex_t IT_0744 = 2*conj(IT_0430);
    const ccomplex_t IT_0745 = -IT_0110 + (-0.5)*IT_0444;
    const ccomplex_t IT_0746 = s_12 + IT_0745;
    const ccomplex_t IT_0747 = m_C_2*IT_0746;
    const ccomplex_t IT_0748 = (-2)*IT_0747;
    const ccomplex_t IT_0749 = IT_0101*IT_0372;
    const ccomplex_t IT_0750 = s_14*IT_0094;
    const ccomplex_t IT_0751 = IT_0734*IT_0750;
    const ccomplex_t IT_0752 = IT_0101*IT_0381;
    const ccomplex_t IT_0753 = s_24*IT_0094;
    const ccomplex_t IT_0754 = s_23*s_34;
    const ccomplex_t IT_0755 = (-2)*IT_0754;
    const ccomplex_t IT_0756 = IT_0753 + IT_0755;
    const ccomplex_t IT_0757 = IT_0734*IT_0756;
    const ccomplex_t IT_0758 = -s_23 + (-0.5)*IT_0749 + 0.25*IT_0751 + 0.5
      *IT_0752 + 0.25*IT_0757;
    const ccomplex_t IT_0759 = s_13 + IT_0758;
    const ccomplex_t IT_0760 = m_C_2*IT_0759;
    const ccomplex_t IT_0761 = 4*IT_0760;
    const ccomplex_t IT_0762 = -IT_0110 + (-0.25)*IT_0459;
    const ccomplex_t IT_0763 = s_12 + IT_0762;
    const ccomplex_t IT_0764 = m_C_2*IT_0763;
    const ccomplex_t IT_0765 = (-4)*IT_0764;
    const ccomplex_t IT_0766 = -IT_0136 + 0.5*IT_0523;
    const ccomplex_t IT_0767 = IT_0002 + IT_0766;
    const ccomplex_t IT_0768 = m_C_2*IT_0767;
    const ccomplex_t IT_0769 = 2*IT_0768;
    const ccomplex_t IT_0770 = (-4)*IT_0535;
    const ccomplex_t IT_0771 = IT_0534 + IT_0770;
    const ccomplex_t IT_0772 = m_C_2*IT_0771;
    const ccomplex_t IT_0773 = IT_0381 + IT_0382;
    const ccomplex_t IT_0774 = s_23*IT_0089*IT_0101;
    const ccomplex_t IT_0775 = (-2)*IT_0774;
    const ccomplex_t IT_0776 = IT_0773 + IT_0775;
    const ccomplex_t IT_0777 = IT_0101*IT_0776;
    const ccomplex_t IT_0778 = -IT_0117 + (-0.5)*IT_0376 + 2*IT_0423 + 0.5
      *IT_0777;
    const ccomplex_t IT_0779 = s_23 + IT_0778;
    const ccomplex_t IT_0780 = m_C_2*IT_0779;
    const ccomplex_t IT_0781 = (-2)*IT_0780;
    const ccomplex_t IT_0782 = (-0.5)*IT_0749 + 0.25*IT_0751;
    const ccomplex_t IT_0783 = s_13 + IT_0782;
    const ccomplex_t IT_0784 = m_C_2*IT_0783;
    const ccomplex_t IT_0785 = 4*IT_0784;
    const ccomplex_t IT_0786 = (-2)*conj(IT_0634);
    const ccomplex_t IT_0787 = (-4)*IT_0441;
    const ccomplex_t IT_0788 = IT_0447 + IT_0787;
    const ccomplex_t IT_0789 = m_C_2*IT_0788;
    const ccomplex_t IT_0790 = (-2)*conj(IT_0207);
    const ccomplex_t IT_0791 = -IT_0117 + 0.5*IT_0777;
    const ccomplex_t IT_0792 = s_23 + IT_0791;
    const ccomplex_t IT_0793 = m_C_2*IT_0792;
    const ccomplex_t IT_0794 = 2*IT_0793;
    const ccomplex_t IT_0795 = -IT_0125 + (-0.5)*IT_0532;
    const ccomplex_t IT_0796 = s_12 + IT_0795;
    const ccomplex_t IT_0797 = m_C_2*IT_0796;
    const ccomplex_t IT_0798 = (-2)*IT_0797;
    const ccomplex_t IT_0799 = (-2)*conj(IT_0369);
    const ccomplex_t IT_0800 = (-0.5)*IT_0752 + (-0.25)*IT_0757;
    const ccomplex_t IT_0801 = s_23 + IT_0800;
    const ccomplex_t IT_0802 = m_C_2*IT_0801;
    const ccomplex_t IT_0803 = (-4)*IT_0802;
    const ccomplex_t IT_0804 = (-2)*conj(IT_0636);
    const ccomplex_t IT_0805 = 0.125*IT_0437*IT_0720 + conj(IT_0066)*IT_0741 +
       0.125*IT_0743*IT_0744 + 0.125*IT_0182*IT_0748 + 0.125*IT_0691*IT_0761 +
       0.125*IT_0210*IT_0765 + 0.125*IT_0716*IT_0769 + 0.125*IT_0718*IT_0772 +
       0.125*IT_0515*IT_0781 + 0.125*IT_0785*IT_0786 + 0.125*IT_0789*IT_0790 +
       0.125*IT_0721*IT_0794 + 0.125*IT_0798*IT_0799 + 0.125*IT_0803*IT_0804;
    const ccomplex_t IT_0806 = IT_0002*IT_0461;
    const ccomplex_t IT_0807 = IT_0002 + IT_0404 + IT_0455;
    const ccomplex_t IT_0808 = m_C_2*IT_0807;
    const ccomplex_t IT_0809 = 4*IT_0179;
    const ccomplex_t IT_0810 = (-6)*IT_0430;
    const ccomplex_t IT_0811 = (-2)*IT_0179;
    const ccomplex_t IT_0812 = -IT_0590;
    const ccomplex_t IT_0813 = s_23*IT_0137;
    const ccomplex_t IT_0814 = s_13*s_14*s_24;
    const ccomplex_t IT_0815 = IT_0813 + IT_0814;
    const ccomplex_t IT_0816 = s_12*s_14*s_34;
    const ccomplex_t IT_0817 = s_23*IT_0002*IT_0094;
    const ccomplex_t IT_0818 = -IT_0816 + -IT_0817;
    const ccomplex_t IT_0819 = IT_0815 + IT_0818;
    const ccomplex_t IT_0820 = IT_0101*IT_0819;
    const ccomplex_t IT_0821 = 0.5*IT_0820;
    const ccomplex_t IT_0822 = IT_0621 + IT_0622 + IT_0812 + IT_0821;
    const ccomplex_t IT_0823 = IT_0619 + IT_0822;
    const ccomplex_t IT_0824 = s_14*s_23*s_24;
    const ccomplex_t IT_0825 = s_13*IT_0121;
    const ccomplex_t IT_0826 = IT_0824 + IT_0825;
    const ccomplex_t IT_0827 = s_12*s_24*s_34;
    const ccomplex_t IT_0828 = s_13*IT_0002*IT_0094;
    const ccomplex_t IT_0829 = -IT_0827 + -IT_0828;
    const ccomplex_t IT_0830 = IT_0826 + IT_0829;
    const ccomplex_t IT_0831 = IT_0101*IT_0830;
    const ccomplex_t IT_0832 = (-0.25)*IT_0831;
    const ccomplex_t IT_0833 = (-0.5)*IT_0605;
    const ccomplex_t IT_0834 = 0.5*IT_0628;
    const ccomplex_t IT_0835 = (-0.5)*IT_0630;
    const ccomplex_t IT_0836 = IT_0603 + IT_0832 + IT_0833 + IT_0834 + IT_0835;
    const ccomplex_t IT_0837 = (-0.5)*IT_0626;
    const ccomplex_t IT_0838 = IT_0836 + IT_0837;
    const ccomplex_t IT_0839 = (-2)*IT_0838;
    const ccomplex_t IT_0840 = IT_0823 + IT_0839;
    const ccomplex_t IT_0841 = (-2)*IT_0393;
    const ccomplex_t IT_0842 = (-0.5)*IT_0510*IT_0512 + (-0.5)*IT_0550*IT_0647
       + (-0.5)*IT_0571*IT_0649 + (-0.5)*IT_0131*IT_0672 + (-0.5)*IT_0130
      *IT_0673 + (-0.5)*IT_0674*IT_0675 + (-0.5)*IT_0148*IT_0679 + (-0.5)
      *IT_0147*IT_0681 + (-0.5)*IT_0686*IT_0691 + (-0.5)*IT_0554*IT_0713 + (-0.5
      )*IT_0574*IT_0730 + IT_0497*IT_0789 + (-3)*IT_0207*(IT_0403 + IT_0407 + (
      -2)*IT_0806) + 2*IT_0482*IT_0808 + (-0.5)*IT_0677*IT_0809 + (-0.5)*IT_0428
      *IT_0810 + (-0.5)*IT_0704*IT_0811 + IT_0840*IT_0841;
    const ccomplex_t IT_0843 = -IT_0136;
    const ccomplex_t IT_0844 = IT_0002 + IT_0843;
    const ccomplex_t IT_0845 = IT_0002*IT_0844;
    const ccomplex_t IT_0846 = -IT_0845;
    const ccomplex_t IT_0847 = IT_0002*IT_0525;
    const ccomplex_t IT_0848 = -conj(IT_0369);
    const ccomplex_t IT_0849 = 0.166666666666667*IT_0679;
    const ccomplex_t IT_0850 = 0.166666666666667*IT_0550;
    const ccomplex_t IT_0851 = 0.166666666666667*IT_0554;
    const ccomplex_t IT_0852 = 0.166666666666667*IT_0710;
    const ccomplex_t IT_0853 = 0.166666666666667*IT_0707;
    const ccomplex_t IT_0854 = 0.166666666666667*IT_0645;
    const ccomplex_t IT_0855 = 0.166666666666667*IT_0204;
    const ccomplex_t IT_0856 = 0.166666666666667*IT_0709;
    const ccomplex_t IT_0857 = 2*conj(IT_0177);
    const ccomplex_t IT_0858 = 0.166666666666667*IT_0600;
    const ccomplex_t IT_0859 = 0.166666666666667*IT_0601;
    const ccomplex_t IT_0860 = 0.166666666666667*IT_0146*IT_0672 + conj
      (IT_0206)*IT_0846 + IT_0847*IT_0848 + IT_0145*IT_0849 + IT_0210*IT_0850 +
       IT_0790*IT_0851 + IT_0540*IT_0852 + IT_0543*IT_0853 + IT_0744*IT_0854 +
       IT_0182*IT_0855 + IT_0698*IT_0856 + IT_0857*IT_0858 + IT_0515*IT_0859;
    const ccomplex_t IT_0861 = 0.25*IT_0659;
    const ccomplex_t IT_0862 = 0.25*IT_0098;
    const ccomplex_t IT_0863 = 0.25*IT_0645;
    const ccomplex_t IT_0864 = 0.25*IT_0652;
    const ccomplex_t IT_0865 = (-2)*conj(IT_0430);
    const ccomplex_t IT_0866 = 0.25*IT_0668;
    const ccomplex_t IT_0867 = 0.25*IT_0625;
    const ccomplex_t IT_0868 = 0.25*IT_0106;
    const ccomplex_t IT_0869 = 0.25*IT_0648;
    const ccomplex_t IT_0870 = (-4)*IT_0087;
    const ccomplex_t IT_0871 = (-2)*IT_0177;
    const ccomplex_t IT_0872 = IT_0066*IT_0655 + 0.25*IT_0657*IT_0714 +
       IT_0660*IT_0862 + IT_0647*IT_0863 + IT_0649*IT_0864 + IT_0646*IT_0867 +
       IT_0515*IT_0868 + IT_0653*IT_0869 + IT_0861*IT_0870 + IT_0866*IT_0871;
    const ccomplex_t IT_0873 = s_34*IT_0094;
    const ccomplex_t IT_0874 = IT_0734*IT_0873;
    const ccomplex_t IT_0875 = pow(m_W, 4);
    const ccomplex_t IT_0876 = (-2)*IT_0089;
    const ccomplex_t IT_0877 = IT_0875 + IT_0876;
    const ccomplex_t IT_0878 = IT_0101*IT_0877;
    const ccomplex_t IT_0879 = (-0.333333333333333)*IT_0874 + (
      -0.333333333333333)*IT_0878;
    const ccomplex_t IT_0880 = IT_0094 + IT_0879;
    const ccomplex_t IT_0881 = IT_0002*IT_0880;
    const ccomplex_t IT_0882 = 3*IT_0881;
    const ccomplex_t IT_0883 = s_12*IT_0875;
    const ccomplex_t IT_0884 = s_13*s_23*IT_0094;
    const ccomplex_t IT_0885 = (-2)*IT_0884;
    const ccomplex_t IT_0886 = IT_0883 + IT_0885;
    const ccomplex_t IT_0887 = IT_0101*IT_0886;
    const ccomplex_t IT_0888 = 0.5*IT_0887;
    const ccomplex_t IT_0889 = s_14*s_23*s_34;
    const ccomplex_t IT_0890 = s_13*s_24*s_34;
    const ccomplex_t IT_0891 = s_12*IT_0089;
    const ccomplex_t IT_0892 = (-2)*IT_0884 + 2*IT_0889 + 2*IT_0890 + (-2)
      *IT_0891;
    const ccomplex_t IT_0893 = IT_0883 + IT_0892;
    const ccomplex_t IT_0894 = IT_0101*IT_0893;
    const ccomplex_t IT_0895 = (-0.5)*IT_0894;
    const ccomplex_t IT_0896 = s_14*s_23*IT_0094;
    const ccomplex_t IT_0897 = s_13*s_24*IT_0094;
    const ccomplex_t IT_0898 = IT_0896 + IT_0897;
    const ccomplex_t IT_0899 = s_12*s_34*IT_0094;
    const ccomplex_t IT_0900 = -IT_0899;
    const ccomplex_t IT_0901 = IT_0898 + IT_0900;
    const ccomplex_t IT_0902 = IT_0734*IT_0901;
    const ccomplex_t IT_0903 = 0.5*IT_0902;
    const ccomplex_t IT_0904 = IT_0112 + IT_0888 + IT_0895 + IT_0903;
    const ccomplex_t IT_0905 = IT_0883 + IT_0889 + IT_0890;
    const ccomplex_t IT_0906 = (-2)*IT_0884 + -IT_0891;
    const ccomplex_t IT_0907 = IT_0905 + IT_0906;
    const ccomplex_t IT_0908 = IT_0101*IT_0907;
    const ccomplex_t IT_0909 = (-0.5)*IT_0908;
    const ccomplex_t IT_0910 = -IT_0102;
    const ccomplex_t IT_0911 = IT_0095 + IT_0103 + IT_0909 + IT_0910;
    const ccomplex_t IT_0912 = s_12*IT_0089*IT_0101;
    const ccomplex_t IT_0913 = -IT_0912;
    const ccomplex_t IT_0914 = IT_0911 + IT_0913;
    const ccomplex_t IT_0915 = IT_0002*IT_0101;
    const ccomplex_t IT_0916 = IT_0384*IT_0915;
    const ccomplex_t IT_0917 = (-0.25)*IT_0916;
    const ccomplex_t IT_0918 = IT_0427 + IT_0836 + IT_0837 + IT_0917;
    const ccomplex_t IT_0919 = s_12*s_13*IT_0094;
    const ccomplex_t IT_0920 = -IT_0816 + (-0.5)*IT_0817 + 0.5*IT_0919;
    const ccomplex_t IT_0921 = IT_0813 + IT_0920;
    const ccomplex_t IT_0922 = IT_0101*IT_0921;
    const ccomplex_t IT_0923 = 2*IT_0922;
    const ccomplex_t IT_0924 = (-0.5)*IT_0923;
    const ccomplex_t IT_0925 = IT_0592 + IT_0812 + IT_0924;
    const ccomplex_t IT_0926 = 4*IT_0160;
    const ccomplex_t IT_0927 = s_12*s_23*IT_0094;
    const ccomplex_t IT_0928 = s_14*s_34*IT_0002;
    const ccomplex_t IT_0929 = 0.5*IT_0828 + (-0.5)*IT_0927 + -IT_0928;
    const ccomplex_t IT_0930 = IT_0824 + IT_0929;
    const ccomplex_t IT_0931 = IT_0101*IT_0930;
    const ccomplex_t IT_0932 = (-2)*IT_0931;
    const ccomplex_t IT_0933 = (-0.5)*IT_0932;
    const ccomplex_t IT_0934 = -IT_0605;
    const ccomplex_t IT_0935 = IT_0603 + IT_0933 + IT_0934;
    const ccomplex_t IT_0936 = (-4)*IT_0160;
    const ccomplex_t IT_0937 = (-2)*IT_0373;
    const ccomplex_t IT_0938 = IT_0372 + IT_0937;
    const ccomplex_t IT_0939 = IT_0101*IT_0938;
    const ccomplex_t IT_0940 = 0.5*IT_0749 + (-0.5)*IT_0751 + (-0.5)*IT_0939;
    const ccomplex_t IT_0941 = s_13 + IT_0940;
    const ccomplex_t IT_0942 = m_C_2*IT_0941;
    const ccomplex_t IT_0943 = 2*IT_0942;
    const ccomplex_t IT_0944 = 0.5*IT_0943;
    const ccomplex_t IT_0945 = IT_0734*IT_0753;
    const ccomplex_t IT_0946 = (-2)*IT_0382;
    const ccomplex_t IT_0947 = IT_0381 + IT_0946;
    const ccomplex_t IT_0948 = IT_0101*IT_0947;
    const ccomplex_t IT_0949 = 0.5*IT_0752 + (-0.5)*IT_0945 + (-0.5)*IT_0948;
    const ccomplex_t IT_0950 = s_23 + IT_0949;
    const ccomplex_t IT_0951 = m_C_2*IT_0950;
    const ccomplex_t IT_0952 = 2*IT_0951;
    const ccomplex_t IT_0953 = 0.5*IT_0952;
    const ccomplex_t IT_0954 = -IT_0089;
    const ccomplex_t IT_0955 = IT_0875 + IT_0954;
    const ccomplex_t IT_0956 = IT_0101*IT_0955;
    const ccomplex_t IT_0957 = -IT_0094 + 0.25*IT_0956;
    const ccomplex_t IT_0958 = IT_0150 + IT_0957;
    const ccomplex_t IT_0959 = IT_0002*IT_0958;
    const ccomplex_t IT_0960 = 4*IT_0959;
    const ccomplex_t IT_0961 = 0.5*IT_0960;
    const ccomplex_t IT_0962 = 2*conj(IT_0497);
    const ccomplex_t IT_0963 = s_13*IT_0875;
    const ccomplex_t IT_0964 = s_14*s_34*IT_0094;
    const ccomplex_t IT_0965 = (-2)*IT_0964;
    const ccomplex_t IT_0966 = IT_0963 + IT_0965;
    const ccomplex_t IT_0967 = IT_0088*IT_0966;
    const ccomplex_t IT_0968 = (-0.666666666666667)*IT_0939 +
       0.166666666666667*IT_0967;
    const ccomplex_t IT_0969 = s_13 + IT_0968;
    const ccomplex_t IT_0970 = m_C_2*IT_0969;
    const ccomplex_t IT_0971 = 6*IT_0970;
    const ccomplex_t IT_0972 = 0.5*IT_0971;
    const ccomplex_t IT_0973 = 2*conj(IT_0482);
    const ccomplex_t IT_0974 = s_23*IT_0875;
    const ccomplex_t IT_0975 = IT_0088*IT_0974;
    const ccomplex_t IT_0976 = (-10)*s_23;
    const ccomplex_t IT_0977 = IT_0975 + IT_0976;
    const ccomplex_t IT_0978 = m_C_2*IT_0977;
    const ccomplex_t IT_0979 = 0.5*IT_0978;
    const ccomplex_t IT_0980 = IT_0877*IT_0915;
    const ccomplex_t IT_0981 = (-0.5)*IT_0980;
    const ccomplex_t IT_0982 = s_34*IT_0002*IT_0088;
    const ccomplex_t IT_0983 = IT_0873*IT_0982;
    const ccomplex_t IT_0984 = (-0.5)*IT_0983;
    const ccomplex_t IT_0985 = 1.5*IT_0097;
    const ccomplex_t IT_0986 = IT_0112 + IT_0888 + IT_0895 + IT_0903 + IT_0981
       + IT_0984 + IT_0985;
    const ccomplex_t IT_0987 = 0.5*IT_0986;
    const ccomplex_t IT_0988 = IT_0660*IT_0987;
    const ccomplex_t IT_0989 = 4*IT_0100;
    const ccomplex_t IT_0990 = -IT_0150;
    const ccomplex_t IT_0991 = IT_0094 + IT_0990;
    const ccomplex_t IT_0992 = IT_0002*IT_0991;
    const ccomplex_t IT_0993 = 2*IT_0992;
    const ccomplex_t IT_0994 = IT_0915*IT_0955;
    const ccomplex_t IT_0995 = (-0.5)*IT_0994;
    const ccomplex_t IT_0996 = IT_0911 + IT_0913 + IT_0993 + IT_0995;
    const ccomplex_t IT_0997 = 0.5*IT_0996;
    const ccomplex_t IT_0998 = (-4)*IT_0109;
    const ccomplex_t IT_0999 = (-2)*IT_0424;
    const ccomplex_t IT_1000 = IT_0375*IT_0915;
    const ccomplex_t IT_1001 = 0.5*IT_1000;
    const ccomplex_t IT_1002 = IT_0619 + IT_0822 + IT_0999 + IT_1001;
    const ccomplex_t IT_1003 = 0.5*IT_1002;
    const ccomplex_t IT_1004 = IT_0998*IT_1003;
    const ccomplex_t IT_1005 = conj(IT_0634)*IT_0882 + 2*conj(IT_0636)*IT_0904
       + (-2)*conj(IT_0177)*IT_0914 + 4*IT_0109*IT_0918 + 0.5*IT_0925*IT_0926 +
       0.5*IT_0935*IT_0936 + IT_0146*IT_0944 + IT_0132*IT_0953 + IT_0744*IT_0961
       + IT_0962*IT_0972 + IT_0973*IT_0979 + IT_0988 + IT_0989*IT_0997 + IT_1004;
    const ccomplex_t IT_1006 = 9*IT_0097;
    const ccomplex_t IT_1007 = (-4)*IT_0096;
    const ccomplex_t IT_1008 = s_14*s_24*IT_0094;
    const ccomplex_t IT_1009 = (-2)*IT_0884 + 4*IT_0889 + (-2)*IT_1008;
    const ccomplex_t IT_1010 = IT_0883 + IT_1009;
    const ccomplex_t IT_1011 = IT_0101*IT_1010;
    const ccomplex_t IT_1012 = -IT_1011;
    const ccomplex_t IT_1013 = pow(m_W, 6);
    const ccomplex_t IT_1014 = s_12*IT_1013;
    const ccomplex_t IT_1015 = s_14*s_24*IT_0875;
    const ccomplex_t IT_1016 = (-2)*IT_1015;
    const ccomplex_t IT_1017 = IT_1014 + IT_1016;
    const ccomplex_t IT_1018 = IT_0088*IT_1017;
    const ccomplex_t IT_1019 = 0.5*IT_1018;
    const ccomplex_t IT_1020 = IT_0095 + IT_1007 + IT_1012 + IT_1019;
    const ccomplex_t IT_1021 = 2*IT_0935;
    const ccomplex_t IT_1022 = (-2)*IT_0925;
    const ccomplex_t IT_1023 = 2*conj(IT_0207)*IT_0840 + conj(IT_0636)*IT_0882
       + 2*conj(IT_0634)*IT_0904 + (-2)*conj(IT_0430)*IT_0914 + 4*IT_0109
      *IT_0918 + IT_0145*IT_0944 + IT_0133*IT_0953 + IT_0857*IT_0961 + IT_0972
      *IT_0973 + IT_0962*IT_0979 + IT_0988 + IT_0989*IT_0997 + IT_1004 + -conj
      (IT_0393)*IT_1006 + (-2)*conj(IT_0390)*IT_1020 + (conj(IT_0369) + 0.25
      *IT_0936)*IT_1021 + (conj(IT_0206) + (-0.25)*IT_0926)*IT_1022;
    const ccomplex_t IT_1024 = -IT_0125 + (-0.25)*IT_0519;
    const ccomplex_t IT_1025 = s_12 + IT_1024;
    const ccomplex_t IT_1026 = m_C_2*IT_1025;
    const ccomplex_t IT_1027 = (-4)*IT_1026;
    const ccomplex_t IT_1028 = (-2)*IT_0558;
    const ccomplex_t IT_1029 = -IT_0369;
    const ccomplex_t IT_1030 = 2*IT_0177;
    const ccomplex_t IT_1031 = (-0.666666666666667)*IT_0393;
    const ccomplex_t IT_1032 = 0.166666666666667*IT_0131*IT_0672 + IT_0206
      *IT_0846 + IT_0130*IT_0849 + IT_0208*IT_0850 + IT_0705*IT_0851 + IT_0543
      *IT_0852 + IT_0540*IT_0853 + IT_0708*IT_0854 + IT_0180*IT_0855 + IT_0691
      *IT_0856 + IT_0512*IT_0859 + (-0.333333333333333)*IT_0482*IT_1027 +
       0.166666666666667*IT_0712*IT_1028 + IT_0847*IT_1029 + IT_0858*IT_1030 +
       IT_0925*IT_1031;
    const ccomplex_t IT_1033 = 1.125*IT_0002;
    const ccomplex_t IT_1034 = IT_0088*IT_1010;
    const ccomplex_t IT_1035 = 0.25*IT_1034;
    const ccomplex_t IT_1036 = (-0.5)*IT_0164;
    const ccomplex_t IT_1037 = (-0.5)*IT_0113;
    const ccomplex_t IT_1038 = (-0.25)*IT_0437;
    const ccomplex_t IT_1039 = (-0.25)*IT_0781;
    const ccomplex_t IT_1040 = IT_0512*IT_1039;
    const ccomplex_t IT_1041 = (-0.25)*IT_0761;
    const ccomplex_t IT_1042 = IT_0698*IT_1041;
    const ccomplex_t IT_1043 = (-0.25)*IT_0789;
    const ccomplex_t IT_1044 = (-0.25)*IT_0748;
    const ccomplex_t IT_1045 = (-0.25)*IT_0765;
    const ccomplex_t IT_1046 = (-0.25)*IT_0743;
    const ccomplex_t IT_1047 = (-0.25)*IT_0794;
    const ccomplex_t IT_1048 = (-0.25)*IT_0769;
    const ccomplex_t IT_1049 = (-0.25)*IT_1027;
    const ccomplex_t IT_1050 = (-2)*IT_0298;
    const ccomplex_t IT_1051 = (-0.25)*IT_0798;
    const ccomplex_t IT_1052 = (-0.25)*IT_0772;
    const ccomplex_t IT_1053 = (-2)*IT_0636;
    const ccomplex_t IT_1054 = (-0.25)*IT_0785;
    const ccomplex_t IT_1055 = (-2)*IT_0634;
    const ccomplex_t IT_1056 = (-0.25)*IT_0803;
    const ccomplex_t IT_1057 = (-2)*IT_0087*IT_0741 + IT_0207*IT_0808 + 4
      *IT_0482*IT_1033 + (-2)*IT_0497*(s_12 + IT_1035 + IT_1036 + IT_1037) +
       IT_0438*IT_1038 + IT_1040 + IT_1042 + IT_0811*IT_1043 + IT_0208*IT_1044 +
       IT_0180*IT_1045 + IT_1030*IT_1046 + IT_0669*IT_1047 + IT_0646*IT_1048 +
       IT_0713*IT_1049 + IT_1050*IT_1051 + IT_0649*IT_1052 + IT_1053*IT_1054 +
       IT_1055*IT_1056;
    const ccomplex_t IT_1058 = (-2)*conj(IT_0206);
    const ccomplex_t IT_1059 = (-2)*conj(IT_0087)*IT_0741 + conj(IT_0207)
      *IT_0808 + IT_0439*IT_1038 + IT_0515*IT_1039 + IT_0691*IT_1041 + IT_0210
      *IT_1044 + IT_0182*IT_1045 + IT_0857*IT_1046 + IT_0865*IT_1047 + IT_0715
      *IT_1048 + IT_0717*IT_1052 + IT_0804*IT_1054 + IT_0786*IT_1056 + IT_1049
      *IT_1058;
    const ccomplex_t IT_1060 = (-2)*IT_0066*IT_0741 + IT_0179*IT_0808 + 4
      *IT_0497*IT_1033 + (-2)*IT_0482*(s_12 + IT_1035 + IT_1036 + IT_1037) +
       IT_0870*IT_1038 + IT_1040 + IT_1042 + IT_0705*IT_1043 + IT_0180*IT_1044 +
       IT_0208*IT_1045 + IT_0708*IT_1046 + IT_0871*IT_1047 + IT_0647*IT_1048 +
       IT_1028*IT_1049 + IT_0730*IT_1051 + IT_0653*IT_1052 + IT_1054*IT_1055 +
       IT_1053*IT_1056;
    const ccomplex_t IT_1061 = 0.25*IT_0926;
    const ccomplex_t IT_1062 = 2*IT_0497;
    const ccomplex_t IT_1063 = 2*IT_0482;
    const ccomplex_t IT_1064 = (-4)*IT_0100;
    const ccomplex_t IT_1065 = IT_0997*IT_1064;
    const ccomplex_t IT_1066 = 0.5*IT_0823;
    const ccomplex_t IT_1067 = IT_0562*IT_0918;
    const ccomplex_t IT_1068 = 0.5*IT_1067;
    const ccomplex_t IT_1069 = (-0.25)*IT_0936;
    const ccomplex_t IT_1070 = (-4)*IT_0179*IT_0838 + IT_0634*IT_0882 + 2
      *IT_0636*IT_0904 + (-2)*IT_0177*IT_0914 + IT_0131*IT_0944 + IT_0147
      *IT_0953 + IT_0708*IT_0961 + IT_0719*IT_0987 + IT_0675*IT_1003 + IT_0841
      *IT_1020 + IT_1021*(IT_0298 + IT_1061) + IT_0972*IT_1062 + IT_0979*IT_1063
       + IT_1065 + IT_0809*IT_1066 + IT_1068 + IT_1022*(IT_0558 + IT_1069);
    const ccomplex_t IT_1071 = (-4)*IT_0207*IT_0838 + IT_0636*IT_0882 + 2
      *IT_0634*IT_0904 + (-2)*IT_0430*IT_0914 + IT_0130*IT_0944 + IT_0148
      *IT_0953 + IT_0719*IT_0987 + IT_0675*IT_1003 + -IT_0393*IT_1006 + IT_0961
      *IT_1030 + IT_1021*(IT_0369 + IT_1061) + IT_0979*IT_1062 + IT_0972*IT_1063
       + IT_1065 + IT_0678*IT_1066 + IT_1068 + IT_1022*(IT_0206 + IT_1069);
    const ccomplex_t IT_1072 = -IT_0120;
    const ccomplex_t IT_1073 = IT_0002 + IT_1072;
    const ccomplex_t IT_1074 = IT_0002*IT_1073;
    const ccomplex_t IT_1075 = -IT_1074;
    const ccomplex_t IT_1076 = (-0.166666666666667)*IT_0681;
    const ccomplex_t IT_1077 = (-0.166666666666667)*IT_0574;
    const ccomplex_t IT_1078 = (-0.166666666666667)*IT_0571;
    const ccomplex_t IT_1079 = (-0.166666666666667)*IT_0707;
    const ccomplex_t IT_1080 = (-0.166666666666667)*IT_0727;
    const ccomplex_t IT_1081 = (-0.166666666666667)*IT_0652;
    const ccomplex_t IT_1082 = (-0.166666666666667)*IT_0726;
    const ccomplex_t IT_1083 = (-0.166666666666667)*IT_0724;
    const ccomplex_t IT_1084 = (-0.166666666666667)*IT_0723;
    const ccomplex_t IT_1085 = (-0.166666666666667)*IT_0584 + (
      -0.166666666666667)*IT_0607 + (-0.166666666666667)*IT_0609 + (
      -0.166666666666667)*IT_0611 + (-0.166666666666667)*IT_0613;
    const ccomplex_t IT_1086 = (-0.166666666666667)*IT_0148*IT_0673 +
       0.333333333333333*IT_0497*IT_0798 + IT_0206*IT_0847 + IT_0935*IT_1031 + (
      -0.166666666666667)*IT_0729*IT_1050 + IT_1029*IT_1075 + IT_0147*IT_1076 +
       IT_0705*IT_1077 + IT_0208*IT_1078 + IT_0540*IT_1079 + IT_0543*IT_1080 +
       IT_0708*IT_1081 + IT_0646*IT_1082 + IT_0691*IT_1083 + IT_1030*IT_1084 +
       IT_0512*IT_1085;
    const ccomplex_t IT_1087 = (-0.166666666666667)*IT_0133*IT_0673 + conj
      (IT_0206)*IT_0847 + IT_0848*IT_1075 + IT_0132*IT_1076 + IT_0790*IT_1077 +
       IT_0210*IT_1078 + IT_0543*IT_1079 + IT_0540*IT_1080 + IT_0744*IT_1081 +
       IT_0715*IT_1082 + IT_0698*IT_1083 + IT_0857*IT_1084 + IT_0515*IT_1085;
    const ccomplex_t IT_1088 = IT_0674*IT_0998;
    const ccomplex_t IT_1089 = IT_0510*IT_0515 + IT_0146*IT_0672 + IT_0145
      *IT_0673 + IT_0133*IT_0679 + IT_0132*IT_0681 + IT_0686*IT_0698 + IT_0550
      *IT_0716 + IT_0571*IT_0717 + IT_0574*IT_0799 + (-12)*conj(IT_0207)*IT_0806
       + IT_0554*IT_1058 + IT_1088;
    const ccomplex_t IT_1090 = IT_0145*IT_0672;
    const ccomplex_t IT_1091 = IT_0133*IT_0681;
    const ccomplex_t IT_1092 = IT_0152*IT_0644;
    const ccomplex_t IT_1093 = IT_0152*IT_0651;
    const ccomplex_t IT_1094 = (-2)*IT_0151;
    const ccomplex_t IT_1095 = IT_0515*IT_1094;
    const ccomplex_t IT_1096 = IT_0543*IT_0707;
    const ccomplex_t IT_1097 = (s_12 + (-0.5)*IT_0113)*(IT_0094 + -IT_0150);
    const ccomplex_t IT_1098 = (-2)*IT_1097;
    const ccomplex_t IT_1099 = IT_0106*IT_0698;
    const ccomplex_t IT_1100 = IT_0007*((IT_0066 + -conj(IT_0066) + IT_0087 + 
      -conj(IT_0087))*IT_0093 + 8*IT_0007*IT_0098 + 8*IT_0100*IT_0106) + IT_0109
      *(16*IT_0109*IT_0114 + IT_0100*(16*IT_0116 + 16*IT_0118) + IT_0129*IT_0134
       + IT_0144*IT_0149) + 8*IT_0100*(IT_0100*IT_0151 + 0.125*IT_0149*IT_0157) 
      + IT_0160*(16*IT_0160*IT_0165 + IT_0134*IT_0168 + IT_0149*IT_0171) +
       IT_0176*(IT_0181 + IT_0183) + IT_0204*(IT_0209 + IT_0211) + (IT_0216 +
       IT_0224)*(IT_0299 + IT_0300 + IT_0370 + IT_0371) + (IT_0380 + IT_0389)*
      (IT_0391 + IT_0392 + IT_0394 + IT_0395) + 4*IT_0398*(IT_0066*conj(IT_0066)
       + 0.25*conj(IT_0087)*IT_0399) + (IT_0403 + IT_0407)*IT_0409 + (IT_0411 +
       IT_0413)*IT_0415 + (IT_0417 + IT_0419)*IT_0421 + IT_0428*(IT_0432 +
       IT_0434 + IT_0436) + IT_0437*(conj(IT_0087)*IT_0438 + IT_0087*IT_0439) + 
      (IT_0451 + IT_0465)*(IT_0484 + IT_0498 + IT_0500 + IT_0501) + (IT_0502 +
       IT_0503)*(IT_0504 + IT_0505 + IT_0506 + IT_0507) + IT_0510*(IT_0511 +
       IT_0513 + IT_0514 + IT_0516) + (IT_0529 + IT_0539)*(IT_0541 + IT_0542 +
       IT_0544 + IT_0545) + IT_0557*(IT_0559 + IT_0560 + IT_0561 + IT_0563 +
       IT_0564 + IT_0565 + IT_0566) + (IT_0572 + IT_0576)*(IT_0560 + IT_0561 +
       IT_0563 + IT_0565 + IT_0566 + IT_0577 + IT_0578 + IT_0579 + IT_0580) +
       IT_0614*(IT_0615 + IT_0617 + IT_0618) + (IT_0584 + IT_0588 + IT_0609 +
       IT_0611 + IT_0625 + IT_0627 + IT_0629 + IT_0631)*(IT_0635 + IT_0637 +
       IT_0639 + IT_0640 + IT_0641) + IT_0670 + conj(IT_0207)*(IT_0130*IT_0672 +
       IT_0131*IT_0673 + IT_0674*IT_0675 + IT_0677*IT_0678 + IT_0147*IT_0679 +
       IT_0148*IT_0681) + IT_0686*(IT_0687 + IT_0688 + IT_0689 + IT_0690 +
       IT_0692 + IT_0693 + IT_0694 + IT_0695 + IT_0696 + IT_0697 + IT_0699) +
       IT_0704*IT_0706 + conj(IT_0206)*(IT_0483*IT_0557 + IT_0512*IT_0601 +
       IT_0130*IT_0672 + IT_0131*IT_0679 + IT_0540*IT_0707 + IT_0600*IT_0708 +
       IT_0691*IT_0709 + IT_0543*IT_0710 + IT_0712*IT_0713) + (-4)*IT_0636
      *IT_0722 + conj(IT_0369)*IT_0731 + 8*IT_0482*IT_0805 + (-2)*conj(IT_0179)
      *IT_0842 + 6*IT_0558*IT_0860 + 4*IT_0634*(conj(IT_0087)*IT_0655 + IT_0439
      *IT_0861 + IT_0719*IT_0862 + IT_0715*IT_0863 + IT_0718*IT_0864 + IT_0865
      *IT_0866 + IT_0716*IT_0867 + IT_0512*IT_0868 + IT_0717*IT_0869) + 4*conj
      (IT_0636)*IT_0872 + 2*IT_0393*IT_1005 + 2*IT_0390*IT_1023 + 6*conj(IT_0558
      )*IT_1032 + (-4)*conj(IT_0497)*IT_1057 + (-4)*IT_0497*IT_1059 + (-4)*conj
      (IT_0482)*IT_1060 + 2*conj(IT_0393)*IT_1070 + 2*conj(IT_0390)*IT_1071 + (
      -6)*conj(IT_0298)*IT_1086 + (-6)*IT_0298*IT_1087 + IT_0179*IT_1089 +
       IT_0207*(IT_0146*IT_0673 + IT_0132*IT_0679 + IT_1088 + IT_1090 + IT_1091)
       + IT_0177*(IT_0106*IT_0691 + IT_0645*IT_0716 + IT_0652*IT_0717 + IT_0431
      *IT_0992 + IT_0146*IT_1092 + IT_0132*IT_1093 + IT_1095) + IT_0430*(IT_0176
      *IT_0210 + IT_0499*IT_0510 + IT_0543*IT_0614 + IT_0106*IT_0691 + IT_0145
      *IT_1092 + IT_0133*IT_1093 + IT_1095) + IT_0206*(IT_0499*IT_0557 + IT_0515
      *IT_0601 + IT_0146*IT_0679 + IT_0698*IT_0709 + IT_0540*IT_0710 + IT_0600
      *IT_0744 + IT_1090 + IT_1096) + IT_0369*(IT_0515*(IT_0584 + IT_0607 +
       IT_0609 + IT_0611 + IT_0613) + IT_0132*IT_0673 + IT_0698*IT_0724 +
       IT_0716*IT_0726 + IT_0540*IT_0727 + IT_0723*IT_0744 + IT_1091 + IT_1096) 
      + conj(IT_0430)*(IT_0176*IT_0208 + IT_0483*IT_0510 + IT_0540*IT_0614 +
       IT_0130*IT_1092 + IT_0148*IT_1093 + IT_0512*IT_1094 + IT_0708*IT_1098 +
       IT_1099) + conj(IT_0177)*(IT_0645*IT_0647 + IT_0649*IT_0652 + IT_0810
      *IT_0992 + IT_0131*IT_1092 + IT_0147*IT_1093 + IT_0512*IT_1094 + IT_1030
      *IT_1098 + IT_1099);
    return create_ccomplex_return(IT_1100);
}

