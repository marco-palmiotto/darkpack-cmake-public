#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_anti_W_W(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_3, 2);
    const ccomplex_t IT_0001 = pow(m_W, -4);
    const ccomplex_t IT_0002 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0003 = IT_0000 + IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -2);
    const ccomplex_t IT_0005 = pow(s_23, 2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = pow(s_24, 2);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = -IT_0006 + -IT_0008;
    const ccomplex_t IT_0010 = IT_0003 + IT_0009;
    const ccomplex_t IT_0011 = m_N_3*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0013, -1);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0015;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0024;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0021 + IT_0023 + 
      -IT_0025 + -IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0017*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0037 = IT_0015*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-0.5)*IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0042 = IT_0015*IT_0041;
    const ccomplex_t IT_0043 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-0.5)*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0040*IT_0047;
    const ccomplex_t IT_0049 = pow(m_W, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = 2*IT_0016;
    const ccomplex_t IT_0054 = IT_0028*IT_0053;
    const ccomplex_t IT_0055 = IT_0031*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0057 = IT_0015*IT_0056;
    const ccomplex_t IT_0058 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-0.5)*IT_0060);
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0064 = IT_0015*IT_0063;
    const ccomplex_t IT_0065 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0066 = IT_0015*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + (-0.5)*IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0062*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = 2*IT_0052 + 0.5*IT_0055 + 2*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0027;
    const ccomplex_t IT_0076 = IT_0017*IT_0075;
    const ccomplex_t IT_0077 = IT_0031*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0053*IT_0075;
    const ccomplex_t IT_0080 = IT_0031*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = cos(alpha);
    const ccomplex_t IT_0083 = IT_0018*IT_0082;
    const ccomplex_t IT_0084 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = sin(alpha);
    const ccomplex_t IT_0087 = IT_0018*IT_0086;
    const ccomplex_t IT_0088 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = IT_0015*IT_0082;
    const ccomplex_t IT_0091 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = IT_0015*IT_0086;
    const ccomplex_t IT_0094 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0085 + -IT_0089 + 
      -IT_0092 + IT_0095);
    const ccomplex_t IT_0097 = cos(beta);
    const ccomplex_t IT_0098 = sin(beta);
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_W*e_em*IT_0015*(IT_0082
      *IT_0097 + IT_0086*IT_0098);
    const ccomplex_t IT_0100 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0096*IT_0102;
    const ccomplex_t IT_0104 = m_N_3*IT_0029;
    const ccomplex_t IT_0105 = m_N_3*IT_0076;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = IT_0031*IT_0107;
    const ccomplex_t IT_0109 = IT_0083*IT_0088;
    const ccomplex_t IT_0110 = IT_0084*IT_0087;
    const ccomplex_t IT_0111 = IT_0090*IT_0094;
    const ccomplex_t IT_0112 = IT_0091*IT_0093;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + IT_0110 + 
      -IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*-1)*m_W*e_em*IT_0015*(IT_0086
      *IT_0097 + -IT_0082*IT_0098);
    const ccomplex_t IT_0116 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0114*IT_0118;
    const ccomplex_t IT_0120 = -IT_0103 + (-0.5)*IT_0108 + -IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0120;
    const ccomplex_t IT_0122 = IT_0031*(IT_0104 + -IT_0105);
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = (-2)*IT_0123;
    const ccomplex_t IT_0125 = IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = 2*conj(IT_0120);
    const ccomplex_t IT_0127 = (-2)*conj(IT_0123);
    const ccomplex_t IT_0128 = IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0077;
    const ccomplex_t IT_0130 = 0.5*IT_0080;
    const ccomplex_t IT_0131 = (-0.5)*IT_0032;
    const ccomplex_t IT_0132 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0133 = IT_0015*IT_0132;
    const ccomplex_t IT_0134 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0135 = IT_0015*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + 0.5*IT_0136);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0140 = IT_0015*IT_0139;
    const ccomplex_t IT_0141 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0142 = IT_0015*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 0.5*IT_0143);
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = IT_0138*IT_0145;
    const ccomplex_t IT_0147 = IT_0050*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0150 = IT_0015*IT_0149;
    const ccomplex_t IT_0151 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0152 = IT_0015*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + 0.5*IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0157 = IT_0015*IT_0156;
    const ccomplex_t IT_0158 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0159 = IT_0015*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + 0.5*IT_0160);
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0155*IT_0162;
    const ccomplex_t IT_0164 = IT_0071*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = (-0.5)*IT_0055 + 2*IT_0148 + 2*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0168 = IT_0083*IT_0167;
    const ccomplex_t IT_0169 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0170 = IT_0087*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0172 = IT_0090*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0174 = IT_0093*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0168 + IT_0170 + 
      -IT_0172 + -IT_0174);
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0118*IT_0176;
    const ccomplex_t IT_0178 = IT_0083*IT_0169;
    const ccomplex_t IT_0179 = IT_0087*IT_0167;
    const ccomplex_t IT_0180 = IT_0090*IT_0173;
    const ccomplex_t IT_0181 = IT_0093*IT_0171;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0178 + -IT_0179 + 
      -IT_0180 + IT_0181);
    const ccomplex_t IT_0183 = IT_0102*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0108 + -IT_0177 + -IT_0183;
    const ccomplex_t IT_0185 = 2*IT_0184;
    const ccomplex_t IT_0186 = -IT_0031*(IT_0104 + -IT_0105);
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = (-2)*IT_0188;
    const ccomplex_t IT_0190 = IT_0185 + IT_0189;
    const ccomplex_t IT_0191 = 2*conj(IT_0184);
    const ccomplex_t IT_0192 = (-2)*conj(IT_0188);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = IT_0012*((conj(IT_0033) + conj(IT_0074) + conj
      (IT_0078) + conj(IT_0081))*IT_0125 + (IT_0033 + IT_0074 + IT_0078 +
       IT_0081)*IT_0128 + (conj(IT_0129) + conj(IT_0130) + conj(IT_0131) + conj
      (IT_0166))*IT_0190 + (IT_0129 + IT_0130 + IT_0131 + IT_0166)*IT_0193);
    const ccomplex_t IT_0195 = s_34*m_N_3*IT_0004;
    const ccomplex_t IT_0196 = s_23*s_34*IT_0004;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = s_24 + IT_0197;
    const ccomplex_t IT_0199 = IT_0195*IT_0198;
    const ccomplex_t IT_0200 = IT_0016*IT_0028;
    const ccomplex_t IT_0201 = IT_0031*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0016*IT_0075;
    const ccomplex_t IT_0204 = IT_0031*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = conj(IT_0202) + conj(IT_0205);
    const ccomplex_t IT_0207 = IT_0202 + IT_0205;
    const ccomplex_t IT_0208 = (-0.5)*IT_0204;
    const ccomplex_t IT_0209 = 0.5*IT_0201;
    const ccomplex_t IT_0210 = conj(IT_0208) + conj(IT_0209);
    const ccomplex_t IT_0211 = IT_0208 + IT_0209;
    const ccomplex_t IT_0212 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0213 = s_13*s_23*IT_0004;
    const ccomplex_t IT_0214 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0215 = -IT_0213 + -IT_0214;
    const ccomplex_t IT_0216 = s_12 + IT_0212 + IT_0215;
    const ccomplex_t IT_0217 = m_N_3*IT_0216;
    const ccomplex_t IT_0218 = conj(IT_0129) + conj(IT_0131);
    const ccomplex_t IT_0219 = IT_0129 + IT_0131;
    const ccomplex_t IT_0220 = cpow((-2)*s_13 + IT_0000 + IT_0049 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0221 = IT_0163*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = cpow((-2)*s_13 + IT_0000 + IT_0049 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0224 = IT_0146*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0055 + (-2)*IT_0222 + (-2)*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0123;
    const ccomplex_t IT_0228 = (-2)*IT_0120;
    const ccomplex_t IT_0229 = IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = 2*conj(IT_0123);
    const ccomplex_t IT_0231 = (-2)*conj(IT_0120);
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = conj(IT_0033) + conj(IT_0078);
    const ccomplex_t IT_0234 = IT_0033 + IT_0078;
    const ccomplex_t IT_0235 = IT_0070*IT_0220;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0048*IT_0223;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = (-0.5)*IT_0055 + (-2)*IT_0236 + (-2)*IT_0238;
    const ccomplex_t IT_0240 = 2*IT_0188;
    const ccomplex_t IT_0241 = (-2)*IT_0184;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = 2*conj(IT_0188);
    const ccomplex_t IT_0244 = (-2)*conj(IT_0184);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0217*(IT_0125*IT_0218 + IT_0128*IT_0219 + 
      (conj(IT_0081) + conj(IT_0226))*IT_0229 + (IT_0081 + IT_0226)*IT_0232 +
       IT_0190*IT_0233 + IT_0193*IT_0234 + (conj(IT_0130) + conj(IT_0239))
      *IT_0242 + (IT_0130 + IT_0239)*IT_0245);
    const ccomplex_t IT_0247 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0248 = IT_0000 + IT_0247;
    const ccomplex_t IT_0249 = pow(s_13, 2);
    const ccomplex_t IT_0250 = IT_0004*IT_0249;
    const ccomplex_t IT_0251 = pow(s_14, 2);
    const ccomplex_t IT_0252 = IT_0004*IT_0251;
    const ccomplex_t IT_0253 = -IT_0250 + -IT_0252;
    const ccomplex_t IT_0254 = IT_0248 + IT_0253;
    const ccomplex_t IT_0255 = m_N_3*IT_0254;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0258 = s_12 + IT_0215 + IT_0257;
    const ccomplex_t IT_0259 = m_N_3*IT_0258;
    const ccomplex_t IT_0260 = conj(IT_0130) + conj(IT_0166);
    const ccomplex_t IT_0261 = IT_0130 + IT_0166;
    const ccomplex_t IT_0262 = s_13*s_34*IT_0004;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = s_14 + IT_0263;
    const ccomplex_t IT_0265 = IT_0195*IT_0264;
    const ccomplex_t IT_0266 = s_13*m_N_3;
    const ccomplex_t IT_0267 = pow(s_34, 2);
    const ccomplex_t IT_0268 = IT_0001*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = 1 + IT_0269;
    const ccomplex_t IT_0271 = IT_0266*IT_0270;
    const ccomplex_t IT_0272 = (-2)*IT_0016;
    const ccomplex_t IT_0273 = IT_0075*IT_0272;
    const ccomplex_t IT_0274 = IT_0031*IT_0273;
    const ccomplex_t IT_0275 = (-0.5)*IT_0274;
    const ccomplex_t IT_0276 = IT_0028*IT_0272;
    const ccomplex_t IT_0277 = IT_0031*IT_0276;
    const ccomplex_t IT_0278 = 0.5*IT_0277;
    const ccomplex_t IT_0279 = conj(IT_0275) + conj(IT_0278);
    const ccomplex_t IT_0280 = 4*IT_0120;
    const ccomplex_t IT_0281 = (-4)*IT_0123;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = IT_0275 + IT_0278;
    const ccomplex_t IT_0284 = 4*conj(IT_0120);
    const ccomplex_t IT_0285 = (-4)*conj(IT_0123);
    const ccomplex_t IT_0286 = IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = (-0.5)*IT_0277;
    const ccomplex_t IT_0288 = 0.5*IT_0274;
    const ccomplex_t IT_0289 = conj(IT_0287) + conj(IT_0288);
    const ccomplex_t IT_0290 = 4*IT_0184;
    const ccomplex_t IT_0291 = (-4)*IT_0188;
    const ccomplex_t IT_0292 = IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = IT_0287 + IT_0288;
    const ccomplex_t IT_0294 = 4*conj(IT_0184);
    const ccomplex_t IT_0295 = (-4)*conj(IT_0188);
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_14*s_34*IT_0004;
    const ccomplex_t IT_0298 = s_13*IT_0049;
    const ccomplex_t IT_0299 = s_14*s_34;
    const ccomplex_t IT_0300 = IT_0298 + IT_0299;
    const ccomplex_t IT_0301 = s_13*IT_0004*IT_0267;
    const ccomplex_t IT_0302 = (-2)*IT_0301;
    const ccomplex_t IT_0303 = IT_0300 + IT_0302;
    const ccomplex_t IT_0304 = IT_0004*IT_0303;
    const ccomplex_t IT_0305 = -IT_0297 + 0.5*IT_0304;
    const ccomplex_t IT_0306 = s_13 + IT_0305;
    const ccomplex_t IT_0307 = m_N_3*IT_0306;
    const ccomplex_t IT_0308 = 2*IT_0307;
    const ccomplex_t IT_0309 = m_N_3*IT_0225;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = (-2)*IT_0310;
    const ccomplex_t IT_0312 = m_N_3*IT_0222;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0313;
    const ccomplex_t IT_0315 = IT_0311 + IT_0314;
    const ccomplex_t IT_0316 = (-2)*conj(IT_0310);
    const ccomplex_t IT_0317 = (-2)*conj(IT_0313);
    const ccomplex_t IT_0318 = IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = IT_0004*IT_0298;
    const ccomplex_t IT_0320 = s_34*IT_0001;
    const ccomplex_t IT_0321 = s_14*IT_0049;
    const ccomplex_t IT_0322 = s_13*s_34;
    const ccomplex_t IT_0323 = (-2)*IT_0322;
    const ccomplex_t IT_0324 = IT_0321 + IT_0323;
    const ccomplex_t IT_0325 = IT_0320*IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0319 + (-0.25)*IT_0325;
    const ccomplex_t IT_0327 = s_13 + IT_0326;
    const ccomplex_t IT_0328 = m_N_3*IT_0327;
    const ccomplex_t IT_0329 = (-4)*IT_0328;
    const ccomplex_t IT_0330 = s_24*s_34*IT_0004;
    const ccomplex_t IT_0331 = s_23*IT_0049;
    const ccomplex_t IT_0332 = s_24*s_34;
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = s_23*IT_0004*IT_0267;
    const ccomplex_t IT_0335 = (-2)*IT_0334;
    const ccomplex_t IT_0336 = IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = IT_0004*IT_0336;
    const ccomplex_t IT_0338 = -IT_0330 + 0.5*IT_0337;
    const ccomplex_t IT_0339 = s_23 + IT_0338;
    const ccomplex_t IT_0340 = m_N_3*IT_0339;
    const ccomplex_t IT_0341 = 2*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0343 = IT_0040*IT_0145;
    const ccomplex_t IT_0344 = IT_0050*IT_0342*IT_0343;
    const ccomplex_t IT_0345 = m_N_3*IT_0073;
    const ccomplex_t IT_0346 = m_N_3*IT_0052;
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0348 = IT_0062*IT_0162;
    const ccomplex_t IT_0349 = IT_0071*IT_0347*IT_0348;
    const ccomplex_t IT_0350 = IT_0344 + -IT_0345 + -IT_0346 + IT_0349;
    const ccomplex_t IT_0351 = 2*IT_0350;
    const ccomplex_t IT_0352 = 2*conj(IT_0350);
    const ccomplex_t IT_0353 = IT_0000*IT_0254;
    const ccomplex_t IT_0354 = -IT_0353;
    const ccomplex_t IT_0355 = 2*conj(IT_0239);
    const ccomplex_t IT_0356 = 2*IT_0033;
    const ccomplex_t IT_0357 = conj(IT_0130)*IT_0356;
    const ccomplex_t IT_0358 = 2*conj(IT_0033);
    const ccomplex_t IT_0359 = IT_0130*IT_0358;
    const ccomplex_t IT_0360 = 2*IT_0129;
    const ccomplex_t IT_0361 = conj(IT_0081)*IT_0360;
    const ccomplex_t IT_0362 = 2*conj(IT_0129);
    const ccomplex_t IT_0363 = IT_0081*IT_0362;
    const ccomplex_t IT_0364 = 2*IT_0130;
    const ccomplex_t IT_0365 = conj(IT_0078)*IT_0364;
    const ccomplex_t IT_0366 = 2*conj(IT_0130);
    const ccomplex_t IT_0367 = IT_0078*IT_0366;
    const ccomplex_t IT_0368 = 2*IT_0131;
    const ccomplex_t IT_0369 = conj(IT_0081)*IT_0368;
    const ccomplex_t IT_0370 = 2*conj(IT_0131);
    const ccomplex_t IT_0371 = IT_0081*IT_0370;
    const ccomplex_t IT_0372 = IT_0000*IT_0010;
    const ccomplex_t IT_0373 = -IT_0372;
    const ccomplex_t IT_0374 = 2*conj(IT_0166);
    const ccomplex_t IT_0375 = s_12*IT_0000;
    const ccomplex_t IT_0376 = s_14*s_24*IT_0001*IT_0005;
    const ccomplex_t IT_0377 = s_13*s_23*IT_0001*IT_0007;
    const ccomplex_t IT_0378 = s_14*s_24*IT_0000*IT_0004;
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = s_13*s_23*IT_0000*IT_0004;
    const ccomplex_t IT_0381 = -IT_0380;
    const ccomplex_t IT_0382 = s_12*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0383 = -IT_0382;
    const ccomplex_t IT_0384 = IT_0375 + IT_0376 + IT_0377 + IT_0379 + IT_0381
       + IT_0383;
    const ccomplex_t IT_0385 = IT_0001*IT_0005*IT_0251;
    const ccomplex_t IT_0386 = s_13*s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0387 = pow(m_N_3, 4);
    const ccomplex_t IT_0388 = IT_0000*IT_0004*IT_0005;
    const ccomplex_t IT_0389 = -IT_0388;
    const ccomplex_t IT_0390 = IT_0000*IT_0004*IT_0251;
    const ccomplex_t IT_0391 = -IT_0390;
    const ccomplex_t IT_0392 = s_12*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = IT_0385 + IT_0386 + IT_0387 + IT_0389 + IT_0391
       + IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0074)*IT_0356;
    const ccomplex_t IT_0396 = conj(IT_0081)*IT_0356;
    const ccomplex_t IT_0397 = IT_0074*IT_0358;
    const ccomplex_t IT_0398 = IT_0081*IT_0358;
    const ccomplex_t IT_0399 = 2*IT_0074;
    const ccomplex_t IT_0400 = conj(IT_0078)*IT_0399;
    const ccomplex_t IT_0401 = 2*conj(IT_0074);
    const ccomplex_t IT_0402 = IT_0078*IT_0401;
    const ccomplex_t IT_0403 = conj(IT_0130)*IT_0360;
    const ccomplex_t IT_0404 = conj(IT_0166)*IT_0360;
    const ccomplex_t IT_0405 = IT_0130*IT_0362;
    const ccomplex_t IT_0406 = IT_0166*IT_0362;
    const ccomplex_t IT_0407 = conj(IT_0131)*IT_0364;
    const ccomplex_t IT_0408 = IT_0131*IT_0366;
    const ccomplex_t IT_0409 = conj(IT_0166)*IT_0368;
    const ccomplex_t IT_0410 = IT_0166*IT_0370;
    const ccomplex_t IT_0411 = 2*IT_0078;
    const ccomplex_t IT_0412 = conj(IT_0081)*IT_0411;
    const ccomplex_t IT_0413 = 2*conj(IT_0078);
    const ccomplex_t IT_0414 = IT_0081*IT_0413;
    const ccomplex_t IT_0415 = -IT_0299;
    const ccomplex_t IT_0416 = IT_0298 + IT_0415;
    const ccomplex_t IT_0417 = IT_0004*IT_0416;
    const ccomplex_t IT_0418 = -IT_0297 + (-0.25)*IT_0417;
    const ccomplex_t IT_0419 = s_13 + IT_0418;
    const ccomplex_t IT_0420 = IT_0000*IT_0419;
    const ccomplex_t IT_0421 = (-4)*IT_0420;
    const ccomplex_t IT_0422 = -IT_0332;
    const ccomplex_t IT_0423 = IT_0331 + IT_0422;
    const ccomplex_t IT_0424 = IT_0004*IT_0423;
    const ccomplex_t IT_0425 = -IT_0330 + (-0.25)*IT_0424;
    const ccomplex_t IT_0426 = s_23 + IT_0425;
    const ccomplex_t IT_0427 = IT_0000*IT_0426;
    const ccomplex_t IT_0428 = (-4)*IT_0427;
    const ccomplex_t IT_0429 = IT_0236 + IT_0238;
    const ccomplex_t IT_0430 = IT_0356*conj(IT_0429);
    const ccomplex_t IT_0431 = -IT_0148 + -IT_0165;
    const ccomplex_t IT_0432 = IT_0356*conj(IT_0431);
    const ccomplex_t IT_0433 = IT_0358*IT_0429;
    const ccomplex_t IT_0434 = IT_0358*IT_0431;
    const ccomplex_t IT_0435 = -IT_0073;
    const ccomplex_t IT_0436 = IT_0360*conj(IT_0435);
    const ccomplex_t IT_0437 = -IT_0052;
    const ccomplex_t IT_0438 = IT_0360*conj(IT_0437);
    const ccomplex_t IT_0439 = conj(IT_0225)*IT_0360;
    const ccomplex_t IT_0440 = conj(IT_0222)*IT_0360;
    const ccomplex_t IT_0441 = IT_0362*IT_0435;
    const ccomplex_t IT_0442 = IT_0362*IT_0437;
    const ccomplex_t IT_0443 = IT_0225*IT_0362;
    const ccomplex_t IT_0444 = IT_0222*IT_0362;
    const ccomplex_t IT_0445 = IT_0368*conj(IT_0435);
    const ccomplex_t IT_0446 = IT_0368*conj(IT_0437);
    const ccomplex_t IT_0447 = conj(IT_0225)*IT_0368;
    const ccomplex_t IT_0448 = conj(IT_0222)*IT_0368;
    const ccomplex_t IT_0449 = IT_0370*IT_0435;
    const ccomplex_t IT_0450 = IT_0370*IT_0437;
    const ccomplex_t IT_0451 = IT_0225*IT_0370;
    const ccomplex_t IT_0452 = IT_0222*IT_0370;
    const ccomplex_t IT_0453 = IT_0411*conj(IT_0429);
    const ccomplex_t IT_0454 = IT_0411*conj(IT_0431);
    const ccomplex_t IT_0455 = IT_0413*IT_0429;
    const ccomplex_t IT_0456 = IT_0413*IT_0431;
    const ccomplex_t IT_0457 = IT_0430 + IT_0432 + IT_0433 + IT_0434 + IT_0436
       + IT_0438 + IT_0439 + IT_0440 + IT_0441 + IT_0442 + IT_0443 + IT_0444 +
       IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 + IT_0451 +
       IT_0452 + IT_0453 + IT_0454 + IT_0455 + IT_0456;
    const ccomplex_t IT_0458 = s_12*IT_0049;
    const ccomplex_t IT_0459 = s_13*s_23;
    const ccomplex_t IT_0460 = (-2)*IT_0459;
    const ccomplex_t IT_0461 = IT_0458 + IT_0460;
    const ccomplex_t IT_0462 = IT_0004*IT_0461;
    const ccomplex_t IT_0463 = (s_13 + -IT_0297)*(s_12 + (-0.5)*IT_0462);
    const ccomplex_t IT_0464 = (-2)*IT_0463;
    const ccomplex_t IT_0465 = (s_23 + -IT_0330)*(s_12 + (-0.5)*IT_0462);
    const ccomplex_t IT_0466 = (-2)*IT_0465;
    const ccomplex_t IT_0467 = 2*IT_0208;
    const ccomplex_t IT_0468 = conj(IT_0131)*IT_0467;
    const ccomplex_t IT_0469 = 2*conj(IT_0208);
    const ccomplex_t IT_0470 = IT_0131*IT_0469;
    const ccomplex_t IT_0471 = conj(IT_0202)*IT_0356;
    const ccomplex_t IT_0472 = conj(IT_0205)*IT_0356;
    const ccomplex_t IT_0473 = IT_0202*IT_0358;
    const ccomplex_t IT_0474 = IT_0205*IT_0358;
    const ccomplex_t IT_0475 = conj(IT_0209)*IT_0360;
    const ccomplex_t IT_0476 = IT_0209*IT_0362;
    const ccomplex_t IT_0477 = conj(IT_0209)*IT_0368;
    const ccomplex_t IT_0478 = IT_0209*IT_0370;
    const ccomplex_t IT_0479 = conj(IT_0205)*IT_0411;
    const ccomplex_t IT_0480 = IT_0205*IT_0413;
    const ccomplex_t IT_0481 = 2*IT_0202;
    const ccomplex_t IT_0482 = conj(IT_0078)*IT_0481;
    const ccomplex_t IT_0483 = 2*conj(IT_0202);
    const ccomplex_t IT_0484 = IT_0078*IT_0483;
    const ccomplex_t IT_0485 = IT_0004*IT_0267;
    const ccomplex_t IT_0486 = (s_12 + (-0.5)*IT_0462)*(IT_0049 + -IT_0485);
    const ccomplex_t IT_0487 = (-2)*IT_0486;
    const ccomplex_t IT_0488 = 2*IT_0209;
    const ccomplex_t IT_0489 = s_14*s_34*IT_0001*IT_0005;
    const ccomplex_t IT_0490 = s_12*s_23;
    const ccomplex_t IT_0491 = s_13*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0492 = s_13*IT_0000;
    const ccomplex_t IT_0493 = IT_0490 + IT_0491 + IT_0492;
    const ccomplex_t IT_0494 = s_13*IT_0004*IT_0005;
    const ccomplex_t IT_0495 = (-2)*IT_0494;
    const ccomplex_t IT_0496 = s_12*s_23*IT_0001*IT_0267;
    const ccomplex_t IT_0497 = -IT_0496;
    const ccomplex_t IT_0498 = s_14*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0499 = -IT_0498;
    const ccomplex_t IT_0500 = IT_0467 + IT_0488;
    const ccomplex_t IT_0501 = 2*conj(IT_0209);
    const ccomplex_t IT_0502 = IT_0469 + IT_0501;
    const ccomplex_t IT_0503 = conj(IT_0130)*IT_0467;
    const ccomplex_t IT_0504 = IT_0130*IT_0469;
    const ccomplex_t IT_0505 = conj(IT_0209)*IT_0364;
    const ccomplex_t IT_0506 = IT_0209*IT_0366;
    const ccomplex_t IT_0507 = conj(IT_0081)*IT_0481;
    const ccomplex_t IT_0508 = IT_0081*IT_0483;
    const ccomplex_t IT_0509 = 2*conj(IT_0205);
    const ccomplex_t IT_0510 = IT_0081*IT_0509;
    const ccomplex_t IT_0511 = (IT_0489 + IT_0493 + IT_0495 + IT_0497 +
       IT_0499)*(conj(IT_0205)*IT_0399 + IT_0205*IT_0401 + conj(IT_0074)*IT_0481
       + IT_0074*IT_0483 + conj(IT_0166)*IT_0500 + IT_0166*IT_0502 + IT_0503 +
       IT_0504 + IT_0505 + IT_0506 + IT_0507 + IT_0508 + IT_0510);
    const ccomplex_t IT_0512 = -IT_0459;
    const ccomplex_t IT_0513 = IT_0458 + IT_0512;
    const ccomplex_t IT_0514 = IT_0004*IT_0513;
    const ccomplex_t IT_0515 = -IT_0213 + (-0.25)*IT_0514;
    const ccomplex_t IT_0516 = s_12 + IT_0515;
    const ccomplex_t IT_0517 = m_N_3*IT_0516;
    const ccomplex_t IT_0518 = (-4)*IT_0517;
    const ccomplex_t IT_0519 = IT_0220*IT_0347*IT_0348;
    const ccomplex_t IT_0520 = m_N_3*IT_0236;
    const ccomplex_t IT_0521 = m_N_3*IT_0238;
    const ccomplex_t IT_0522 = IT_0223*IT_0342*IT_0343;
    const ccomplex_t IT_0523 = -IT_0519 + IT_0520 + IT_0521 + -IT_0522;
    const ccomplex_t IT_0524 = 2*IT_0081;
    const ccomplex_t IT_0525 = IT_0399 + IT_0524;
    const ccomplex_t IT_0526 = 2*conj(IT_0081);
    const ccomplex_t IT_0527 = IT_0401 + IT_0526;
    const ccomplex_t IT_0528 = m_N_3*IT_0148;
    const ccomplex_t IT_0529 = m_N_3*IT_0165;
    const ccomplex_t IT_0530 = conj(IT_0528) + conj(IT_0529);
    const ccomplex_t IT_0531 = (-2)*IT_0081;
    const ccomplex_t IT_0532 = (-2)*IT_0226;
    const ccomplex_t IT_0533 = IT_0531 + IT_0532;
    const ccomplex_t IT_0534 = IT_0528 + IT_0529;
    const ccomplex_t IT_0535 = (-2)*conj(IT_0081);
    const ccomplex_t IT_0536 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0537 = IT_0535 + IT_0536;
    const ccomplex_t IT_0538 = IT_0069*IT_0155;
    const ccomplex_t IT_0539 = IT_0071*IT_0347*IT_0538;
    const ccomplex_t IT_0540 = IT_0047*IT_0138;
    const ccomplex_t IT_0541 = IT_0050*IT_0342*IT_0540;
    const ccomplex_t IT_0542 = IT_0539 + IT_0541;
    const ccomplex_t IT_0543 = 2*IT_0542;
    const ccomplex_t IT_0544 = 2*conj(IT_0542);
    const ccomplex_t IT_0545 = IT_0220*IT_0347*IT_0538;
    const ccomplex_t IT_0546 = -IT_0545;
    const ccomplex_t IT_0547 = IT_0223*IT_0342*IT_0540;
    const ccomplex_t IT_0548 = -IT_0547;
    const ccomplex_t IT_0549 = conj(IT_0350) + conj(IT_0546) + conj(IT_0548);
    const ccomplex_t IT_0550 = IT_0350 + IT_0546 + IT_0548;
    const ccomplex_t IT_0551 = conj(IT_0546) + conj(IT_0548);
    const ccomplex_t IT_0552 = 2*IT_0166;
    const ccomplex_t IT_0553 = IT_0546 + IT_0548;
    const ccomplex_t IT_0554 = 2*IT_0226;
    const ccomplex_t IT_0555 = 2*conj(IT_0226);
    const ccomplex_t IT_0556 = IT_0260*IT_0315 + IT_0261*IT_0318 + conj
      (IT_0239)*IT_0351 + IT_0239*IT_0352 + conj(IT_0523)*IT_0525 + IT_0523
      *IT_0527 + IT_0530*IT_0533 + IT_0534*IT_0537 + conj(IT_0081)*IT_0543 +
       IT_0081*IT_0544 + IT_0364*IT_0549 + IT_0366*IT_0550 + IT_0551*IT_0552 +
       IT_0374*IT_0553 + conj(IT_0542)*IT_0554 + IT_0542*IT_0555;
    const ccomplex_t IT_0557 = s_14*s_24;
    const ccomplex_t IT_0558 = (-2)*IT_0557;
    const ccomplex_t IT_0559 = IT_0458 + IT_0558;
    const ccomplex_t IT_0560 = IT_0004*IT_0559;
    const ccomplex_t IT_0561 = (s_12 + -IT_0213)*(s_12 + (-0.5)*IT_0560);
    const ccomplex_t IT_0562 = (-2)*IT_0561;
    const ccomplex_t IT_0563 = s_24*s_34*IT_0001*IT_0249;
    const ccomplex_t IT_0564 = s_12*s_13;
    const ccomplex_t IT_0565 = s_13*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0566 = s_23*IT_0000;
    const ccomplex_t IT_0567 = IT_0564 + IT_0565 + IT_0566;
    const ccomplex_t IT_0568 = s_23*IT_0004*IT_0249;
    const ccomplex_t IT_0569 = (-2)*IT_0568;
    const ccomplex_t IT_0570 = s_12*s_13*IT_0001*IT_0267;
    const ccomplex_t IT_0571 = -IT_0570;
    const ccomplex_t IT_0572 = s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0573 = -IT_0572;
    const ccomplex_t IT_0574 = IT_0563 + IT_0567 + IT_0569 + IT_0571 + IT_0573;
    const ccomplex_t IT_0575 = s_13*s_23*IT_0001*IT_0251;
    const ccomplex_t IT_0576 = s_14*s_24*IT_0001*IT_0249;
    const ccomplex_t IT_0577 = s_12*s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0578 = -IT_0577;
    const ccomplex_t IT_0579 = IT_0375 + IT_0379 + IT_0381 + IT_0575 + IT_0576
       + IT_0578;
    const ccomplex_t IT_0580 = IT_0001*IT_0007*IT_0249;
    const ccomplex_t IT_0581 = IT_0000*IT_0004*IT_0007;
    const ccomplex_t IT_0582 = -IT_0581;
    const ccomplex_t IT_0583 = IT_0000*IT_0004*IT_0249;
    const ccomplex_t IT_0584 = -IT_0583;
    const ccomplex_t IT_0585 = s_12*s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = IT_0386 + IT_0387 + IT_0580 + IT_0582 + IT_0584
       + IT_0586;
    const ccomplex_t IT_0588 = IT_0579 + IT_0587;
    const ccomplex_t IT_0589 = conj(IT_0226)*IT_0356;
    const ccomplex_t IT_0590 = IT_0226*IT_0358;
    const ccomplex_t IT_0591 = conj(IT_0239)*IT_0360;
    const ccomplex_t IT_0592 = IT_0239*IT_0362;
    const ccomplex_t IT_0593 = conj(IT_0078)*IT_0554;
    const ccomplex_t IT_0594 = IT_0078*IT_0555;
    const ccomplex_t IT_0595 = s_12*s_14*s_34*IT_0004;
    const ccomplex_t IT_0596 = s_13*s_14*s_24*IT_0004;
    const ccomplex_t IT_0597 = -IT_0596;
    const ccomplex_t IT_0598 = s_23*IT_0004*IT_0251;
    const ccomplex_t IT_0599 = -IT_0598;
    const ccomplex_t IT_0600 = 2*IT_0239;
    const ccomplex_t IT_0601 = conj(IT_0275)*IT_0364;
    const ccomplex_t IT_0602 = conj(IT_0278)*IT_0364;
    const ccomplex_t IT_0603 = IT_0275*IT_0366;
    const ccomplex_t IT_0604 = IT_0278*IT_0366;
    const ccomplex_t IT_0605 = conj(IT_0288)*IT_0524;
    const ccomplex_t IT_0606 = IT_0288*IT_0526;
    const ccomplex_t IT_0607 = 2*IT_0287;
    const ccomplex_t IT_0608 = conj(IT_0081)*IT_0607;
    const ccomplex_t IT_0609 = 2*conj(IT_0287);
    const ccomplex_t IT_0610 = IT_0081*IT_0609;
    const ccomplex_t IT_0611 = 2*conj(IT_0275);
    const ccomplex_t IT_0612 = (IT_0563 + IT_0567 + IT_0569 + IT_0571 +
       IT_0595 + IT_0597 + IT_0599)*(IT_0278*IT_0355 + IT_0289*IT_0554 + IT_0293
      *IT_0555 + conj(IT_0278)*IT_0600 + IT_0601 + IT_0602 + IT_0603 + IT_0604 +
       IT_0605 + IT_0606 + IT_0608 + IT_0610 + IT_0239*IT_0611);
    const ccomplex_t IT_0613 = s_12*s_24*s_34*IT_0004;
    const ccomplex_t IT_0614 = s_13*IT_0004*IT_0007;
    const ccomplex_t IT_0615 = -IT_0614;
    const ccomplex_t IT_0616 = s_14*s_23*s_24*IT_0004;
    const ccomplex_t IT_0617 = -IT_0616;
    const ccomplex_t IT_0618 = 2*IT_0275;
    const ccomplex_t IT_0619 = (IT_0489 + IT_0493 + IT_0495 + IT_0497 +
       IT_0613 + IT_0615 + IT_0617)*(IT_0278*IT_0374 + IT_0289*IT_0399 + IT_0293
      *IT_0401 + conj(IT_0278)*IT_0552 + IT_0601 + IT_0602 + IT_0603 + IT_0604 +
       IT_0605 + IT_0606 + IT_0608 + IT_0610 + IT_0166*IT_0611 + conj(IT_0166)
      *IT_0618);
    const ccomplex_t IT_0620 = 2*IT_0565;
    const ccomplex_t IT_0621 = IT_0564 + IT_0566 + IT_0569 + IT_0597 + IT_0599
       + IT_0620;
    const ccomplex_t IT_0622 = 2*IT_0491;
    const ccomplex_t IT_0623 = IT_0490 + IT_0492 + IT_0495 + IT_0615 + IT_0617
       + IT_0622;
    const ccomplex_t IT_0624 = conj(IT_0287)*IT_0356;
    const ccomplex_t IT_0625 = conj(IT_0288)*IT_0356;
    const ccomplex_t IT_0626 = IT_0287*IT_0358;
    const ccomplex_t IT_0627 = IT_0288*IT_0358;
    const ccomplex_t IT_0628 = conj(IT_0275)*IT_0360;
    const ccomplex_t IT_0629 = conj(IT_0278)*IT_0360;
    const ccomplex_t IT_0630 = IT_0275*IT_0362;
    const ccomplex_t IT_0631 = IT_0278*IT_0362;
    const ccomplex_t IT_0632 = conj(IT_0278)*IT_0368;
    const ccomplex_t IT_0633 = IT_0278*IT_0370;
    const ccomplex_t IT_0634 = conj(IT_0288)*IT_0411;
    const ccomplex_t IT_0635 = IT_0288*IT_0413;
    const ccomplex_t IT_0636 = conj(IT_0078)*IT_0607;
    const ccomplex_t IT_0637 = IT_0078*IT_0609;
    const ccomplex_t IT_0638 = conj(IT_0131)*IT_0618;
    const ccomplex_t IT_0639 = IT_0131*IT_0611;
    const ccomplex_t IT_0640 = (s_12 + -IT_0214)*(s_12 + (-0.5)*IT_0462);
    const ccomplex_t IT_0641 = (-2)*IT_0640;
    const ccomplex_t IT_0642 = 4*IT_0033;
    const ccomplex_t IT_0643 = -IT_0250;
    const ccomplex_t IT_0644 = IT_0000 + IT_0643;
    const ccomplex_t IT_0645 = IT_0000*IT_0644;
    const ccomplex_t IT_0646 = -IT_0645;
    const ccomplex_t IT_0647 = -IT_0006;
    const ccomplex_t IT_0648 = IT_0000 + IT_0647;
    const ccomplex_t IT_0649 = IT_0000*IT_0648;
    const ccomplex_t IT_0650 = -IT_0649;
    const ccomplex_t IT_0651 = conj(IT_0081)*IT_0130;
    const ccomplex_t IT_0652 = 6*IT_0651;
    const ccomplex_t IT_0653 = IT_0081*conj(IT_0130);
    const ccomplex_t IT_0654 = 6*IT_0653;
    const ccomplex_t IT_0655 = -IT_0252;
    const ccomplex_t IT_0656 = IT_0000 + IT_0655;
    const ccomplex_t IT_0657 = IT_0000*IT_0656;
    const ccomplex_t IT_0658 = -IT_0657;
    const ccomplex_t IT_0659 = -IT_0008;
    const ccomplex_t IT_0660 = IT_0000 + IT_0659;
    const ccomplex_t IT_0661 = IT_0000*IT_0660;
    const ccomplex_t IT_0662 = -IT_0661;
    const ccomplex_t IT_0663 = IT_0658 + IT_0662;
    const ccomplex_t IT_0664 = 6*IT_0033;
    const ccomplex_t IT_0665 = conj(IT_0131)*IT_0664;
    const ccomplex_t IT_0666 = 6*conj(IT_0033);
    const ccomplex_t IT_0667 = IT_0131*IT_0666;
    const ccomplex_t IT_0668 = conj(IT_0078)*IT_0129;
    const ccomplex_t IT_0669 = 6*IT_0668;
    const ccomplex_t IT_0670 = IT_0078*conj(IT_0129);
    const ccomplex_t IT_0671 = 6*IT_0670;
    const ccomplex_t IT_0672 = conj(IT_0078)*IT_0131;
    const ccomplex_t IT_0673 = 6*IT_0672;
    const ccomplex_t IT_0674 = IT_0078*conj(IT_0131);
    const ccomplex_t IT_0675 = 6*IT_0674;
    const ccomplex_t IT_0676 = IT_0000*IT_0049;
    const ccomplex_t IT_0677 = s_13*s_14*s_34*IT_0004;
    const ccomplex_t IT_0678 = 0.5*IT_0251 + (-0.5)*IT_0676 + -IT_0677;
    const ccomplex_t IT_0679 = IT_0249 + IT_0678;
    const ccomplex_t IT_0680 = IT_0004*IT_0679;
    const ccomplex_t IT_0681 = IT_0000 + IT_0655 + IT_0680;
    const ccomplex_t IT_0682 = m_N_3*IT_0681;
    const ccomplex_t IT_0683 = 4*conj(IT_0129);
    const ccomplex_t IT_0684 = conj(IT_0131) + 0.25*IT_0683;
    const ccomplex_t IT_0685 = -IT_0676;
    const ccomplex_t IT_0686 = IT_0007 + IT_0685;
    const ccomplex_t IT_0687 = IT_0004*IT_0686;
    const ccomplex_t IT_0688 = -IT_0000;
    const ccomplex_t IT_0689 = IT_0008 + IT_0688;
    const ccomplex_t IT_0690 = (-4)*IT_0689;
    const ccomplex_t IT_0691 = IT_0687 + IT_0690;
    const ccomplex_t IT_0692 = m_N_3*IT_0691;
    const ccomplex_t IT_0693 = IT_0358 + IT_0413;
    const ccomplex_t IT_0694 = s_13*s_24*s_34*IT_0004;
    const ccomplex_t IT_0695 = -IT_0459 + (-2)*IT_0557 + 2*IT_0694;
    const ccomplex_t IT_0696 = IT_0458 + IT_0695;
    const ccomplex_t IT_0697 = IT_0004*IT_0696;
    const ccomplex_t IT_0698 = -IT_0213 + (-0.5)*IT_0697;
    const ccomplex_t IT_0699 = s_12 + IT_0698;
    const ccomplex_t IT_0700 = m_N_3*IT_0699;
    const ccomplex_t IT_0701 = (-2)*IT_0700;
    const ccomplex_t IT_0702 = IT_0526 + IT_0555;
    const ccomplex_t IT_0703 = IT_0004*IT_0331;
    const ccomplex_t IT_0704 = s_24*IT_0049;
    const ccomplex_t IT_0705 = IT_0320*IT_0704;
    const ccomplex_t IT_0706 = (-0.5)*IT_0703 + 0.25*IT_0705;
    const ccomplex_t IT_0707 = s_23 + IT_0706;
    const ccomplex_t IT_0708 = m_N_3*IT_0707;
    const ccomplex_t IT_0709 = 4*IT_0708;
    const ccomplex_t IT_0710 = 2*conj(IT_0288);
    const ccomplex_t IT_0711 = IT_0609 + IT_0710;
    const ccomplex_t IT_0712 = 2*conj(IT_0278);
    const ccomplex_t IT_0713 = IT_0611 + IT_0712;
    const ccomplex_t IT_0714 = pow(m_W, 4);
    const ccomplex_t IT_0715 = s_12*IT_0714;
    const ccomplex_t IT_0716 = s_13*s_23*IT_0049;
    const ccomplex_t IT_0717 = s_14*s_23*s_34;
    const ccomplex_t IT_0718 = s_13*s_24*s_34;
    const ccomplex_t IT_0719 = s_12*IT_0267;
    const ccomplex_t IT_0720 = s_14*s_24*IT_0049;
    const ccomplex_t IT_0721 = (-2)*IT_0716 + 2*IT_0717 + 2*IT_0718 + (-2)
      *IT_0719 + (-2)*IT_0720;
    const ccomplex_t IT_0722 = IT_0715 + IT_0721;
    const ccomplex_t IT_0723 = IT_0001*IT_0722;
    const ccomplex_t IT_0724 = (-0.25)*IT_0723;
    const ccomplex_t IT_0725 = -IT_0214;
    const ccomplex_t IT_0726 = -IT_0213;
    const ccomplex_t IT_0727 = s_12 + IT_0724 + IT_0725 + IT_0726;
    const ccomplex_t IT_0728 = 8*conj(IT_0528);
    const ccomplex_t IT_0729 = (-8)*conj(IT_0542);
    const ccomplex_t IT_0730 = IT_0728 + IT_0729;
    const ccomplex_t IT_0731 = IT_0004*IT_0458;
    const ccomplex_t IT_0732 = (-0.5)*IT_0731;
    const ccomplex_t IT_0733 = s_14*s_23;
    const ccomplex_t IT_0734 = s_13*s_24;
    const ccomplex_t IT_0735 = s_12*s_34;
    const ccomplex_t IT_0736 = -IT_0734 + -IT_0735;
    const ccomplex_t IT_0737 = IT_0733 + IT_0736;
    const ccomplex_t IT_0738 = IT_0320*IT_0737;
    const ccomplex_t IT_0739 = (-0.25)*IT_0738;
    const ccomplex_t IT_0740 = s_12 + IT_0732 + IT_0739;
    const ccomplex_t IT_0741 = 8*conj(IT_0123);
    const ccomplex_t IT_0742 = (-8)*conj(IT_0120);
    const ccomplex_t IT_0743 = IT_0741 + IT_0742;
    const ccomplex_t IT_0744 = (-2)*IT_0716 + 4*IT_0718 + (-2)*IT_0720;
    const ccomplex_t IT_0745 = IT_0715 + IT_0744;
    const ccomplex_t IT_0746 = IT_0001*IT_0745;
    const ccomplex_t IT_0747 = 0.25*IT_0746;
    const ccomplex_t IT_0748 = (-0.5)*IT_0560;
    const ccomplex_t IT_0749 = (-0.5)*IT_0462;
    const ccomplex_t IT_0750 = s_12 + IT_0747 + IT_0748 + IT_0749;
    const ccomplex_t IT_0751 = 8*conj(IT_0546);
    const ccomplex_t IT_0752 = (-8)*conj(IT_0310);
    const ccomplex_t IT_0753 = (-8)*conj(IT_0313);
    const ccomplex_t IT_0754 = IT_0751 + IT_0752 + IT_0753;
    const ccomplex_t IT_0755 = IT_0249 + IT_0685;
    const ccomplex_t IT_0756 = IT_0004*IT_0755;
    const ccomplex_t IT_0757 = IT_0250 + IT_0688;
    const ccomplex_t IT_0758 = (-4)*IT_0757;
    const ccomplex_t IT_0759 = IT_0756 + IT_0758;
    const ccomplex_t IT_0760 = m_N_3*IT_0759;
    const ccomplex_t IT_0761 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0762 = (-2)*conj(IT_0130);
    const ccomplex_t IT_0763 = IT_0761 + IT_0762;
    const ccomplex_t IT_0764 = s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0765 = 2*IT_0007 + -IT_0676 + (-2)*IT_0764;
    const ccomplex_t IT_0766 = IT_0005 + IT_0765;
    const ccomplex_t IT_0767 = IT_0004*IT_0766;
    const ccomplex_t IT_0768 = -IT_0006 + 0.5*IT_0767;
    const ccomplex_t IT_0769 = IT_0000 + IT_0768;
    const ccomplex_t IT_0770 = m_N_3*IT_0769;
    const ccomplex_t IT_0771 = 2*IT_0770;
    const ccomplex_t IT_0772 = (-2)*conj(IT_0074);
    const ccomplex_t IT_0773 = IT_0535 + IT_0772;
    const ccomplex_t IT_0774 = (-2)*IT_0459 + -IT_0557 + 2*IT_0694;
    const ccomplex_t IT_0775 = IT_0458 + IT_0774;
    const ccomplex_t IT_0776 = IT_0004*IT_0775;
    const ccomplex_t IT_0777 = -IT_0214 + (-0.5)*IT_0776;
    const ccomplex_t IT_0778 = s_12 + IT_0777;
    const ccomplex_t IT_0779 = m_N_3*IT_0778;
    const ccomplex_t IT_0780 = (-2)*IT_0779;
    const ccomplex_t IT_0781 = (-2)*conj(IT_0129);
    const ccomplex_t IT_0782 = (-2)*conj(IT_0131);
    const ccomplex_t IT_0783 = IT_0781 + IT_0782;
    const ccomplex_t IT_0784 = m_N_3*IT_0426;
    const ccomplex_t IT_0785 = (-4)*IT_0784;
    const ccomplex_t IT_0786 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0787 = (-2)*conj(IT_0205);
    const ccomplex_t IT_0788 = IT_0786 + IT_0787;
    const ccomplex_t IT_0789 = (-2)*IT_0267;
    const ccomplex_t IT_0790 = IT_0714 + IT_0789;
    const ccomplex_t IT_0791 = IT_0001*IT_0790;
    const ccomplex_t IT_0792 = 0.25*IT_0791;
    const ccomplex_t IT_0793 = 1 + IT_0792;
    const ccomplex_t IT_0794 = IT_0000*IT_0793;
    const ccomplex_t IT_0795 = (-8)*conj(IT_0350);
    const ccomplex_t IT_0796 = 0.25*IT_0308*IT_0502 + IT_0682*IT_0684 + 0.25
      *IT_0692*IT_0693 + 0.25*IT_0701*IT_0702 + 0.25*IT_0709*IT_0711 + 0.25
      *IT_0329*IT_0713 + 0.25*IT_0727*IT_0730 + 0.25*IT_0740*IT_0743 + 0.25
      *IT_0750*IT_0754 + 0.25*IT_0760*IT_0763 + 0.25*IT_0771*IT_0773 + 0.25
      *IT_0780*IT_0783 + 0.25*IT_0785*IT_0788 + 0.25*IT_0794*IT_0795;
    const ccomplex_t IT_0797 = 0.25*IT_0308*IT_0502 + IT_0682*IT_0684 + 0.25
      *IT_0692*IT_0693 + 0.25*IT_0701*IT_0702 + 0.25*IT_0709*IT_0711 + 0.25
      *IT_0329*IT_0713 + 0.25*IT_0727*IT_0729 + 0.25*IT_0740*IT_0743 + 0.25
      *IT_0750*(IT_0752 + IT_0753) + 0.25*IT_0760*IT_0763 + 0.25*IT_0771*IT_0773
       + 0.25*IT_0780*IT_0783 + 0.25*IT_0785*IT_0788 + 0.25*IT_0794*IT_0795;
    const ccomplex_t IT_0798 = -IT_0485;
    const ccomplex_t IT_0799 = IT_0049 + IT_0798;
    const ccomplex_t IT_0800 = IT_0000*IT_0799;
    const ccomplex_t IT_0801 = (-6)*IT_0202;
    const ccomplex_t IT_0802 = (-6)*conj(IT_0202);
    const ccomplex_t IT_0803 = -IT_0297;
    const ccomplex_t IT_0804 = s_13 + IT_0803;
    const ccomplex_t IT_0805 = IT_0000*IT_0804;
    const ccomplex_t IT_0806 = -IT_0330;
    const ccomplex_t IT_0807 = s_23 + IT_0806;
    const ccomplex_t IT_0808 = IT_0000*IT_0807;
    const ccomplex_t IT_0809 = IT_0805 + IT_0808;
    const ccomplex_t IT_0810 = (-6)*conj(IT_0208);
    const ccomplex_t IT_0811 = IT_0078*IT_0810;
    const ccomplex_t IT_0812 = (-6)*IT_0033;
    const ccomplex_t IT_0813 = conj(IT_0209)*IT_0812;
    const ccomplex_t IT_0814 = (-6)*conj(IT_0033);
    const ccomplex_t IT_0815 = IT_0209*IT_0814;
    const ccomplex_t IT_0816 = IT_0129*conj(IT_0205);
    const ccomplex_t IT_0817 = (-6)*IT_0816;
    const ccomplex_t IT_0818 = conj(IT_0129)*IT_0205;
    const ccomplex_t IT_0819 = (-6)*IT_0818;
    const ccomplex_t IT_0820 = IT_0131*conj(IT_0205);
    const ccomplex_t IT_0821 = (-6)*IT_0820;
    const ccomplex_t IT_0822 = conj(IT_0131)*IT_0205;
    const ccomplex_t IT_0823 = (-6)*IT_0822;
    const ccomplex_t IT_0824 = conj(IT_0131)*IT_0801;
    const ccomplex_t IT_0825 = IT_0131*IT_0802;
    const ccomplex_t IT_0826 = (-6)*IT_0209;
    const ccomplex_t IT_0827 = conj(IT_0078)*IT_0826;
    const ccomplex_t IT_0828 = (-6)*conj(IT_0209);
    const ccomplex_t IT_0829 = IT_0078*IT_0828;
    const ccomplex_t IT_0830 = s_12*IT_0270;
    const ccomplex_t IT_0831 = (-4)*conj(IT_0120);
    const ccomplex_t IT_0832 = conj(IT_0123) + 0.25*IT_0831;
    const ccomplex_t IT_0833 = (-4)*conj(IT_0184);
    const ccomplex_t IT_0834 = conj(IT_0188) + 0.25*IT_0833;
    const ccomplex_t IT_0835 = (-4)*IT_0120;
    const ccomplex_t IT_0836 = (-4)*IT_0184;
    const ccomplex_t IT_0837 = s_23*m_N_3;
    const ccomplex_t IT_0838 = IT_0270*IT_0837;
    const ccomplex_t IT_0839 = IT_0123 + 0.25*IT_0835;
    const ccomplex_t IT_0840 = IT_0000*IT_0270;
    const ccomplex_t IT_0841 = IT_0184 + IT_0350;
    const ccomplex_t IT_0842 = conj(IT_0310) + conj(IT_0313) + conj(IT_0542);
    const ccomplex_t IT_0843 = 0.25*IT_0286*IT_0523 + 0.25*IT_0282*conj
      (IT_0523) + conj(IT_0188)*(0.25*IT_0280 + 0.25*IT_0281 + IT_0310 + IT_0313
       + IT_0542) + 0.25*IT_0296*(IT_0534 + IT_0546 + IT_0548) + 0.25*IT_0292*
      (IT_0530 + conj(IT_0546) + conj(IT_0548)) + 0.25*(IT_0310 + IT_0313 +
       IT_0542)*IT_0833 + (conj(IT_0184) + conj(IT_0350))*IT_0839 + (conj
      (IT_0123) + 0.25*IT_0831)*IT_0841 + IT_0188*(0.25*IT_0284 + 0.25*IT_0285 +
       IT_0842) + 0.25*IT_0836*IT_0842;
    const ccomplex_t IT_0844 = (1 + 0.5*IT_0268)*(IT_0458 + (-2)*IT_0459);
    const ccomplex_t IT_0845 = (-4)*IT_0287;
    const ccomplex_t IT_0846 = IT_0270*IT_0676;
    const ccomplex_t IT_0847 = (-4)*IT_0275;
    const ccomplex_t IT_0848 = (-4)*IT_0278;
    const ccomplex_t IT_0849 = IT_0847 + IT_0848;
    const ccomplex_t IT_0850 = s_23*s_34;
    const ccomplex_t IT_0851 = (-2)*IT_0850;
    const ccomplex_t IT_0852 = IT_0704 + IT_0851;
    const ccomplex_t IT_0853 = IT_0320*IT_0852;
    const ccomplex_t IT_0854 = (-0.5)*IT_0703 + (-0.25)*IT_0853;
    const ccomplex_t IT_0855 = s_23 + IT_0854;
    const ccomplex_t IT_0856 = m_N_3*IT_0855;
    const ccomplex_t IT_0857 = (-4)*IT_0856;
    const ccomplex_t IT_0858 = IT_0320*IT_0321;
    const ccomplex_t IT_0859 = (-0.5)*IT_0319 + 0.25*IT_0858;
    const ccomplex_t IT_0860 = s_13 + IT_0859;
    const ccomplex_t IT_0861 = m_N_3*IT_0860;
    const ccomplex_t IT_0862 = 4*IT_0861;
    const ccomplex_t IT_0863 = (-4)*conj(IT_0287);
    const ccomplex_t IT_0864 = (-4)*conj(IT_0275);
    const ccomplex_t IT_0865 = 2*IT_0310;
    const ccomplex_t IT_0866 = 2*IT_0313;
    const ccomplex_t IT_0867 = IT_0865 + IT_0866;
    const ccomplex_t IT_0868 = s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0869 = -IT_0459 + (-2)*IT_0557 + 2*IT_0868;
    const ccomplex_t IT_0870 = IT_0458 + IT_0869;
    const ccomplex_t IT_0871 = IT_0004*IT_0870;
    const ccomplex_t IT_0872 = -IT_0213 + (-0.5)*IT_0871;
    const ccomplex_t IT_0873 = s_12 + IT_0872;
    const ccomplex_t IT_0874 = m_N_3*IT_0873;
    const ccomplex_t IT_0875 = (-2)*IT_0874;
    const ccomplex_t IT_0876 = 4*IT_0130;
    const ccomplex_t IT_0877 = 6*IT_0074;
    const ccomplex_t IT_0878 = 6*IT_0226;
    const ccomplex_t IT_0879 = s_12 + IT_0726;
    const ccomplex_t IT_0880 = IT_0000*IT_0879;
    const ccomplex_t IT_0881 = (-6)*IT_0074;
    const ccomplex_t IT_0882 = (-6)*IT_0226;
    const ccomplex_t IT_0883 = 2*IT_0251 + -IT_0676 + (-2)*IT_0677;
    const ccomplex_t IT_0884 = IT_0249 + IT_0883;
    const ccomplex_t IT_0885 = IT_0004*IT_0884;
    const ccomplex_t IT_0886 = -IT_0250 + 0.5*IT_0885;
    const ccomplex_t IT_0887 = IT_0000 + IT_0886;
    const ccomplex_t IT_0888 = m_N_3*IT_0887;
    const ccomplex_t IT_0889 = 2*IT_0888;
    const ccomplex_t IT_0890 = (-2)*IT_0542;
    const ccomplex_t IT_0891 = IT_0251 + IT_0685;
    const ccomplex_t IT_0892 = IT_0004*IT_0891;
    const ccomplex_t IT_0893 = IT_0252 + IT_0688;
    const ccomplex_t IT_0894 = (-4)*IT_0893;
    const ccomplex_t IT_0895 = IT_0892 + IT_0894;
    const ccomplex_t IT_0896 = m_N_3*IT_0895;
    const ccomplex_t IT_0897 = 4*conj(IT_0033);
    const ccomplex_t IT_0898 = (-2)*IT_0459 + -IT_0557 + 2*IT_0868;
    const ccomplex_t IT_0899 = IT_0458 + IT_0898;
    const ccomplex_t IT_0900 = IT_0004*IT_0899;
    const ccomplex_t IT_0901 = -IT_0214 + (-0.5)*IT_0900;
    const ccomplex_t IT_0902 = s_12 + IT_0901;
    const ccomplex_t IT_0903 = m_N_3*IT_0902;
    const ccomplex_t IT_0904 = (-2)*IT_0903;
    const ccomplex_t IT_0905 = (-2)*conj(IT_0542);
    const ccomplex_t IT_0906 = m_N_3*IT_0419;
    const ccomplex_t IT_0907 = (-4)*IT_0906;
    const ccomplex_t IT_0908 = 2*IT_0205;
    const ccomplex_t IT_0909 = IT_0481 + IT_0908;
    const ccomplex_t IT_0910 = (-6)*IT_0208;
    const ccomplex_t IT_0911 = (-2)*IT_0350;
    const ccomplex_t IT_0912 = IT_0364 + IT_0600;
    const ccomplex_t IT_0913 = (IT_0000 + -IT_0006)*(s_12 + (-0.5)*IT_0560);
    const ccomplex_t IT_0914 = 2*IT_0913;
    const ccomplex_t IT_0915 = (-2)*IT_0130;
    const ccomplex_t IT_0916 = (-2)*IT_0166;
    const ccomplex_t IT_0917 = IT_0915 + IT_0916;
    const ccomplex_t IT_0918 = IT_0005 + IT_0685;
    const ccomplex_t IT_0919 = IT_0004*IT_0918;
    const ccomplex_t IT_0920 = IT_0006 + IT_0688;
    const ccomplex_t IT_0921 = (-4)*IT_0920;
    const ccomplex_t IT_0922 = IT_0919 + IT_0921;
    const ccomplex_t IT_0923 = m_N_3*IT_0922;
    const ccomplex_t IT_0924 = 2*conj(IT_0310);
    const ccomplex_t IT_0925 = 2*conj(IT_0313);
    const ccomplex_t IT_0926 = IT_0924 + IT_0925;
    const ccomplex_t IT_0927 = (-2)*conj(IT_0350);
    const ccomplex_t IT_0928 = (-2)*IT_0716 + 4*IT_0717 + (-2)*IT_0720;
    const ccomplex_t IT_0929 = IT_0715 + IT_0928;
    const ccomplex_t IT_0930 = IT_0001*IT_0929;
    const ccomplex_t IT_0931 = 0.25*IT_0930;
    const ccomplex_t IT_0932 = s_12 + IT_0748 + IT_0749 + IT_0931;
    const ccomplex_t IT_0933 = IT_0733 + IT_0735;
    const ccomplex_t IT_0934 = -IT_0734;
    const ccomplex_t IT_0935 = IT_0933 + IT_0934;
    const ccomplex_t IT_0936 = IT_0320*IT_0935;
    const ccomplex_t IT_0937 = 0.25*IT_0936;
    const ccomplex_t IT_0938 = s_12 + IT_0732 + IT_0937;
    const ccomplex_t IT_0939 = 8*IT_0120;
    const ccomplex_t IT_0940 = (-8)*IT_0123;
    const ccomplex_t IT_0941 = IT_0939 + IT_0940;
    const ccomplex_t IT_0942 = 0.5*IT_0007 + (-0.5)*IT_0676 + -IT_0764;
    const ccomplex_t IT_0943 = IT_0005 + IT_0942;
    const ccomplex_t IT_0944 = IT_0004*IT_0943;
    const ccomplex_t IT_0945 = IT_0000 + IT_0659 + IT_0944;
    const ccomplex_t IT_0946 = m_N_3*IT_0945;
    const ccomplex_t IT_0947 = 8*IT_0310;
    const ccomplex_t IT_0948 = (-2)*IT_0033;
    const ccomplex_t IT_0949 = (-2)*IT_0074;
    const ccomplex_t IT_0950 = 0.125*IT_0341*IT_0481 + 0.125*IT_0607*IT_0857 +
       IT_0542*IT_0932 + 0.125*IT_0938*IT_0941 + 0.5*(IT_0078 + 0.25*IT_0642)
      *IT_0946 + IT_0727*(IT_0313 + 0.125*IT_0947) + 0.125*IT_0904*IT_0948 +
       0.125*IT_0923*IT_0949;
    const ccomplex_t IT_0951 = 8*conj(IT_0120);
    const ccomplex_t IT_0952 = (-8)*conj(IT_0123);
    const ccomplex_t IT_0953 = IT_0951 + IT_0952;
    const ccomplex_t IT_0954 = 8*conj(IT_0310);
    const ccomplex_t IT_0955 = (-2)*conj(IT_0033);
    const ccomplex_t IT_0956 = IT_0000*IT_0216;
    const ccomplex_t IT_0957 = conj(IT_0130)*IT_0948;
    const ccomplex_t IT_0958 = IT_0130*IT_0955;
    const ccomplex_t IT_0959 = (-2)*IT_0129;
    const ccomplex_t IT_0960 = conj(IT_0081)*IT_0959;
    const ccomplex_t IT_0961 = IT_0081*IT_0781;
    const ccomplex_t IT_0962 = conj(IT_0078)*IT_0915;
    const ccomplex_t IT_0963 = IT_0078*IT_0762;
    const ccomplex_t IT_0964 = (-2)*IT_0131;
    const ccomplex_t IT_0965 = conj(IT_0081)*IT_0964;
    const ccomplex_t IT_0966 = IT_0081*IT_0782;
    const ccomplex_t IT_0967 = s_12 + IT_0725;
    const ccomplex_t IT_0968 = IT_0000*IT_0967;
    const ccomplex_t IT_0969 = 4*IT_0078;
    const ccomplex_t IT_0970 = IT_0642 + IT_0969;
    const ccomplex_t IT_0971 = (-2)*IT_0239;
    const ccomplex_t IT_0972 = IT_0000*IT_0258;
    const ccomplex_t IT_0973 = IT_0373*IT_0552 + IT_0354*IT_0600 + IT_0315
      *IT_0692 + IT_0351*IT_0896 + IT_0890*IT_0904 + IT_0809*IT_0910 + 4*IT_0542
      *IT_0946 + (-12)*IT_0131*IT_0968 + IT_0641*IT_0970 + IT_0956*IT_0971 +
       IT_0916*IT_0972;
    const ccomplex_t IT_0974 = (-2)*conj(IT_0166);
    const ccomplex_t IT_0975 = (IT_0000 + -IT_0252)*(s_12 + (-0.5)*IT_0462);
    const ccomplex_t IT_0976 = 2*IT_0975;
    const ccomplex_t IT_0977 = (IT_0000 + -IT_0008)*(s_12 + (-0.5)*IT_0462);
    const ccomplex_t IT_0978 = 2*IT_0977;
    const ccomplex_t IT_0979 = conj(IT_0033)*IT_0948;
    const ccomplex_t IT_0980 = conj(IT_0078)*IT_0948;
    const ccomplex_t IT_0981 = IT_0078*IT_0955;
    const ccomplex_t IT_0982 = conj(IT_0129)*IT_0959;
    const ccomplex_t IT_0983 = conj(IT_0131)*IT_0959;
    const ccomplex_t IT_0984 = IT_0131*IT_0781;
    const ccomplex_t IT_0985 = conj(IT_0131)*IT_0964;
    const ccomplex_t IT_0986 = (-2)*IT_0078;
    const ccomplex_t IT_0987 = conj(IT_0078)*IT_0986;
    const ccomplex_t IT_0988 = -IT_0557;
    const ccomplex_t IT_0989 = IT_0458 + IT_0988;
    const ccomplex_t IT_0990 = IT_0004*IT_0989;
    const ccomplex_t IT_0991 = -IT_0214 + (-0.25)*IT_0990;
    const ccomplex_t IT_0992 = s_12 + IT_0991;
    const ccomplex_t IT_0993 = m_N_3*IT_0992;
    const ccomplex_t IT_0994 = (-4)*IT_0993;
    const ccomplex_t IT_0995 = IT_0360 + IT_0368;
    const ccomplex_t IT_0996 = IT_0362 + IT_0370;
    const ccomplex_t IT_0997 = IT_0959 + IT_0964;
    const ccomplex_t IT_0998 = (-2)*conj(IT_0078);
    const ccomplex_t IT_0999 = (conj(IT_0310) + conj(IT_0313))*IT_0356 + 
      (IT_0310 + IT_0313)*IT_0358 + IT_0523*IT_0783 + IT_0218*IT_0890 + IT_0219
      *IT_0905 + conj(IT_0078)*(IT_0867 + IT_0911) + IT_0078*(IT_0926 + IT_0927)
       + IT_0549*IT_0948 + IT_0550*IT_0955 + IT_0551*IT_0986 + IT_0530*IT_0995 +
       IT_0534*IT_0996 + conj(IT_0523)*IT_0997 + IT_0553*IT_0998;
    const ccomplex_t IT_1000 = conj(IT_0130)*IT_0915;
    const ccomplex_t IT_1001 = conj(IT_0081)*IT_0531;
    const ccomplex_t IT_1002 = (IT_0000 + -IT_0250)*(s_12 + (-0.5)*IT_0560);
    const ccomplex_t IT_1003 = 2*IT_1002;
    const ccomplex_t IT_1004 = 1.125*IT_0000;
    const ccomplex_t IT_1005 = 8*IT_0313;
    const ccomplex_t IT_1006 = IT_0947 + IT_1005;
    const ccomplex_t IT_1007 = 4*IT_0129;
    const ccomplex_t IT_1008 = (-2)*IT_0202;
    const ccomplex_t IT_1009 = IT_0341*IT_0467 + IT_0618*IT_0857 + IT_0607
      *IT_0862 + IT_0399*IT_0875 + IT_0532*IT_0889 + IT_0356*IT_0896 + IT_0915
      *IT_0923 + 8*(IT_0184 + -IT_0188)*IT_0938 + IT_0904*IT_0959 + (-16)
      *IT_0542*IT_1004 + IT_0794*IT_1006 + 4*IT_0946*(IT_0131 + 0.25*IT_1007) +
       IT_0907*IT_1008;
    const ccomplex_t IT_1010 = (IT_0401 + IT_0555)*IT_0562 + 6*conj(IT_0166)
      *IT_0650 + IT_0318*IT_0701 + IT_0352*IT_0875 + ((-12)*conj(IT_0130) + (-6)
      *conj(IT_0166))*IT_0880 + IT_0772*IT_0914 + IT_0905*IT_0923 + IT_0771
      *IT_0926 + IT_0889*IT_0927 + IT_0536*IT_1003;
    const ccomplex_t IT_1011 = (-2)*IT_0287;
    const ccomplex_t IT_1012 = 0.125*IT_0532*IT_0701 + IT_0542*IT_0727 + 0.125
      *IT_0399*IT_0771 + 0.125*IT_0481*IT_0785 + 0.125*IT_0740*IT_0941 + 0.125
      *IT_0692*IT_0948 + 0.125*IT_0750*IT_1006 + 0.125*IT_0709*IT_1011;
    const ccomplex_t IT_1013 = s_34*IT_0004;
    const ccomplex_t IT_1014 = IT_0733 + IT_0734;
    const ccomplex_t IT_1015 = s_13*s_23*s_34*IT_0004;
    const ccomplex_t IT_1016 = (-2)*IT_1015;
    const ccomplex_t IT_1017 = IT_1014 + IT_1016;
    const ccomplex_t IT_1018 = IT_1013*IT_1017;
    const ccomplex_t IT_1019 = (-2)*IT_0208;
    const ccomplex_t IT_1020 = (-2)*conj(IT_0208);
    const ccomplex_t IT_1021 = (-2)*IT_0209;
    const ccomplex_t IT_1022 = (-2)*conj(IT_0209);
    const ccomplex_t IT_1023 = (-2)*IT_0205;
    const ccomplex_t IT_1024 = (-2)*conj(IT_0287);
    const ccomplex_t IT_1025 = (-2)*IT_0275;
    const ccomplex_t IT_1026 = (-2)*conj(IT_0275);
    const ccomplex_t IT_1027 = IT_0293*IT_0786 + IT_0288*IT_0787 + IT_0289
      *IT_1008 + conj(IT_0205)*IT_1011 + IT_0279*IT_1019 + IT_0283*IT_1020 +
       conj(IT_0278)*IT_1021 + IT_0278*IT_1022 + conj(IT_0288)*IT_1023 + IT_0205
      *IT_1024 + conj(IT_0209)*IT_1025 + IT_0209*IT_1026;
    const ccomplex_t IT_1028 = 2*IT_0288;
    const ccomplex_t IT_1029 = IT_0607 + IT_1028;
    const ccomplex_t IT_1030 = 2*IT_0278;
    const ccomplex_t IT_1031 = IT_0618 + IT_1030;
    const ccomplex_t IT_1032 = 8*IT_0528;
    const ccomplex_t IT_1033 = 8*IT_0529;
    const ccomplex_t IT_1034 = (-8)*IT_0542;
    const ccomplex_t IT_1035 = (-8)*IT_0350;
    const ccomplex_t IT_1036 = 8*IT_0523;
    const ccomplex_t IT_1037 = (-8)*(IT_0184 + -IT_0188)*IT_0740 + IT_0533
      *IT_0760 + IT_0308*IT_0909 + IT_0701*IT_0912 + IT_0771*IT_0917 + IT_0682
      *IT_0970 + IT_0780*(IT_0948 + IT_0986) + IT_0692*IT_0995 + (16*IT_0310 +
       16*IT_0313 + (-16)*IT_0546 + (-16)*IT_0548)*IT_1004 + IT_0785*(IT_1019 +
       IT_1021) + IT_0329*IT_1029 + IT_0709*IT_1031 + IT_0794*(IT_1032 + IT_1033
       + IT_1034) + IT_0727*IT_1035 + IT_0750*IT_1036;
    const ccomplex_t IT_1038 = -conj(IT_0188);
    const ccomplex_t IT_1039 = IT_0341*IT_0469 + IT_0611*IT_0857 + IT_0609
      *IT_0862 + IT_0401*IT_0875 + IT_0536*IT_0889 + IT_0358*IT_0896 + IT_0781
      *IT_0904 + IT_0786*IT_0907 + IT_0762*IT_0923 + 8*conj(IT_0350)*IT_0932 + 4
      *IT_0684*IT_0946 + 8*IT_0794*(conj(IT_0313) + 0.125*IT_0954) + (-16)*conj
      (IT_0542)*IT_1004 + 8*IT_0938*(conj(IT_0184) + IT_1038);
    const ccomplex_t IT_1040 = IT_0355 + IT_0366;
    const ccomplex_t IT_1041 = IT_0483 + IT_0509;
    const ccomplex_t IT_1042 = 4*conj(IT_0078);
    const ccomplex_t IT_1043 = IT_0897 + IT_1042;
    const ccomplex_t IT_1044 = IT_0762 + IT_0974;
    const ccomplex_t IT_1045 = IT_0329*IT_0711 + IT_0709*IT_0713 + IT_0537
      *IT_0760 + IT_0730*IT_0794 + IT_0727*IT_0795 + IT_0692*IT_0996 + IT_0780*
      (IT_0955 + IT_0998) + (16*conj(IT_0310) + 16*conj(IT_0313) + (-16)*conj
      (IT_0546))*IT_1004 + IT_0785*(IT_1020 + IT_1022) + (-8)*IT_0740*(conj
      (IT_0184) + IT_1038) + IT_0701*IT_1040 + IT_0308*IT_1041 + IT_0682*IT_1043
       + IT_0771*IT_1044;
    const ccomplex_t IT_1046 = IT_0915 + IT_0971;
    const ccomplex_t IT_1047 = 0.166666666666667*IT_1003;
    const ccomplex_t IT_1048 = 0.166666666666667*IT_0562;
    const ccomplex_t IT_1049 = 0.166666666666667*IT_0889;
    const ccomplex_t IT_1050 = 0.166666666666667*IT_0354*IT_0356 +
       0.166666666666667*IT_0500*IT_0574 + 0.166666666666667*IT_0368*IT_0588 +
       0.166666666666667*(IT_0563 + IT_0567 + IT_0569 + IT_0571 + IT_0595 +
       IT_0597 + IT_0599)*IT_0618 + 0.166666666666667*IT_0760*IT_0867 +
       0.166666666666667*IT_0646*IT_0878 + IT_0081*(IT_0646 + -IT_0880) +
       0.166666666666667*IT_0880*IT_0881 + 0.166666666666667*IT_0948*IT_0956 +
       IT_1046*IT_1047 + IT_0364*IT_1048 + IT_0890*IT_1049;
    const ccomplex_t IT_1051 = 4*IT_0081;
    const ccomplex_t IT_1052 = IT_0239*IT_0646 + IT_0166*IT_0650 +
       0.166666666666667*IT_0315*IT_0701 + 0.166666666666667*IT_0771*IT_0867 +
       0.166666666666667*IT_0351*IT_0875 + (-2)*(IT_0130 + 1./2*IT_0166 + 1./2
      *IT_0239)*IT_0880 + 0.166666666666667*(IT_0489 + IT_0490 + IT_0491 +
       IT_0492 + IT_0495 + IT_0497 + IT_0499 + IT_0563 + IT_0567 + IT_0569 +
       IT_0571 + IT_0573)*IT_0908 + 0.166666666666667*IT_0890*IT_0923 +
       0.166666666666667*IT_0914*IT_0949 + IT_0532*IT_1047 + IT_0911*IT_1049 +
       IT_1048*(IT_0399 + IT_0554 + IT_1051);
    const ccomplex_t IT_1053 = pow(m_W, 6);
    const ccomplex_t IT_1054 = s_12*IT_1053;
    const ccomplex_t IT_1055 = s_14*s_24*IT_0714;
    const ccomplex_t IT_1056 = (-2)*IT_1055;
    const ccomplex_t IT_1057 = IT_1054 + IT_1056;
    const ccomplex_t IT_1058 = IT_0001*IT_1057;
    const ccomplex_t IT_1059 = 0.5*IT_1058;
    const ccomplex_t IT_1060 = (-4)*IT_0459;
    const ccomplex_t IT_1061 = IT_0458 + IT_1059 + IT_1060;
    const ccomplex_t IT_1062 = IT_0004*IT_0745;
    const ccomplex_t IT_1063 = -IT_1062;
    const ccomplex_t IT_1064 = IT_1061 + IT_1063;
    const ccomplex_t IT_1065 = conj(IT_0222) + conj(IT_0225);
    const ccomplex_t IT_1066 = -IT_0566;
    const ccomplex_t IT_1067 = IT_0564 + IT_1066;
    const ccomplex_t IT_1068 = s_13*s_14*s_24;
    const ccomplex_t IT_1069 = s_23*IT_0000*IT_0049;
    const ccomplex_t IT_1070 = s_12*s_13*IT_0049;
    const ccomplex_t IT_1071 = s_24*s_34*IT_0000;
    const ccomplex_t IT_1072 = 0.5*IT_1069 + (-0.5)*IT_1070 + -IT_1071;
    const ccomplex_t IT_1073 = IT_1068 + IT_1072;
    const ccomplex_t IT_1074 = IT_0004*IT_1073;
    const ccomplex_t IT_1075 = (-2)*IT_1074;
    const ccomplex_t IT_1076 = (-0.5)*IT_1075;
    const ccomplex_t IT_1077 = IT_1067 + IT_1076;
    const ccomplex_t IT_1078 = 4*conj(IT_0081);
    const ccomplex_t IT_1079 = 4*conj(IT_0226);
    const ccomplex_t IT_1080 = IT_1078 + IT_1079;
    const ccomplex_t IT_1081 = 0.25*IT_1058;
    const ccomplex_t IT_1082 = IT_0458 + IT_1081;
    const ccomplex_t IT_1083 = (-8)*conj(IT_0435);
    const ccomplex_t IT_1084 = (-8)*conj(IT_0437);
    const ccomplex_t IT_1085 = IT_1083 + IT_1084;
    const ccomplex_t IT_1086 = -IT_0613;
    const ccomplex_t IT_1087 = s_14*s_23*s_24;
    const ccomplex_t IT_1088 = s_13*IT_0007;
    const ccomplex_t IT_1089 = IT_1087 + IT_1088;
    const ccomplex_t IT_1090 = s_12*s_24*s_34;
    const ccomplex_t IT_1091 = s_13*IT_0000*IT_0049;
    const ccomplex_t IT_1092 = -IT_1090 + -IT_1091;
    const ccomplex_t IT_1093 = IT_1089 + IT_1092;
    const ccomplex_t IT_1094 = IT_0004*IT_1093;
    const ccomplex_t IT_1095 = (-0.5)*IT_1094;
    const ccomplex_t IT_1096 = IT_0492 + IT_0615 + IT_0616 + IT_1086 + IT_1095;
    const ccomplex_t IT_1097 = (-4)*conj(IT_0033);
    const ccomplex_t IT_1098 = (-4)*conj(IT_0078);
    const ccomplex_t IT_1099 = IT_1097 + IT_1098;
    const ccomplex_t IT_1100 = -IT_0492;
    const ccomplex_t IT_1101 = IT_0490 + IT_1100;
    const ccomplex_t IT_1102 = s_12*s_23*IT_0049;
    const ccomplex_t IT_1103 = -IT_1090 + (-0.5)*IT_1091 + 0.5*IT_1102;
    const ccomplex_t IT_1104 = IT_1088 + IT_1103;
    const ccomplex_t IT_1105 = IT_0004*IT_1104;
    const ccomplex_t IT_1106 = 2*IT_1105;
    const ccomplex_t IT_1107 = (-0.5)*IT_1106;
    const ccomplex_t IT_1108 = IT_1101 + IT_1107;
    const ccomplex_t IT_1109 = (-4)*conj(IT_0074);
    const ccomplex_t IT_1110 = (-4)*conj(IT_0081);
    const ccomplex_t IT_1111 = IT_1109 + IT_1110;
    const ccomplex_t IT_1112 = (-2)*IT_0332;
    const ccomplex_t IT_1113 = IT_0331 + IT_1112;
    const ccomplex_t IT_1114 = IT_0004*IT_1113;
    const ccomplex_t IT_1115 = 0.5*IT_0703 + (-0.5)*IT_0705 + (-0.5)*IT_1114;
    const ccomplex_t IT_1116 = s_23 + IT_1115;
    const ccomplex_t IT_1117 = m_N_3*IT_1116;
    const ccomplex_t IT_1118 = 2*IT_1117;
    const ccomplex_t IT_1119 = (-0.25)*IT_1118;
    const ccomplex_t IT_1120 = (-2)*IT_0299;
    const ccomplex_t IT_1121 = IT_0298 + IT_1120;
    const ccomplex_t IT_1122 = IT_0004*IT_1121;
    const ccomplex_t IT_1123 = 0.5*IT_0319 + (-0.5)*IT_0858 + (-0.5)*IT_1122;
    const ccomplex_t IT_1124 = s_13 + IT_1123;
    const ccomplex_t IT_1125 = m_N_3*IT_1124;
    const ccomplex_t IT_1126 = 2*IT_1125;
    const ccomplex_t IT_1127 = (-0.25)*IT_1126;
    const ccomplex_t IT_1128 = -IT_0267;
    const ccomplex_t IT_1129 = IT_0714 + IT_1128;
    const ccomplex_t IT_1130 = IT_0004*IT_1129;
    const ccomplex_t IT_1131 = -IT_0049 + 0.25*IT_1130;
    const ccomplex_t IT_1132 = IT_0485 + IT_1131;
    const ccomplex_t IT_1133 = IT_0000*IT_1132;
    const ccomplex_t IT_1134 = 4*IT_1133;
    const ccomplex_t IT_1135 = (-0.25)*IT_1134;
    const ccomplex_t IT_1136 = s_34*IT_0049;
    const ccomplex_t IT_1137 = IT_0320*IT_1136;
    const ccomplex_t IT_1138 = IT_0004*IT_0790;
    const ccomplex_t IT_1139 = (-0.333333333333333)*IT_1137 + (
      -0.333333333333333)*IT_1138;
    const ccomplex_t IT_1140 = IT_0049 + IT_1139;
    const ccomplex_t IT_1141 = IT_0000*IT_1140;
    const ccomplex_t IT_1142 = 3*IT_1141;
    const ccomplex_t IT_1143 = (-0.25)*IT_1142;
    const ccomplex_t IT_1144 = (-2)*conj(IT_0546);
    const ccomplex_t IT_1145 = (-2)*conj(IT_0548);
    const ccomplex_t IT_1146 = IT_0926 + IT_1144 + IT_1145;
    const ccomplex_t IT_1147 = s_23*IT_0714;
    const ccomplex_t IT_1148 = s_24*s_34*IT_0049;
    const ccomplex_t IT_1149 = (-2)*IT_1148;
    const ccomplex_t IT_1150 = IT_1147 + IT_1149;
    const ccomplex_t IT_1151 = IT_0001*IT_1150;
    const ccomplex_t IT_1152 = (-0.666666666666667)*IT_1114 +
       0.166666666666667*IT_1151;
    const ccomplex_t IT_1153 = s_23 + IT_1152;
    const ccomplex_t IT_1154 = m_N_3*IT_1153;
    const ccomplex_t IT_1155 = 6*IT_1154;
    const ccomplex_t IT_1156 = (-0.25)*IT_1155;
    const ccomplex_t IT_1157 = (-2)*conj(IT_0523);
    const ccomplex_t IT_1158 = s_13*IT_0714;
    const ccomplex_t IT_1159 = IT_0001*IT_1158;
    const ccomplex_t IT_1160 = (-10)*s_13;
    const ccomplex_t IT_1161 = IT_1159 + IT_1160;
    const ccomplex_t IT_1162 = m_N_3*IT_1161;
    const ccomplex_t IT_1163 = (-0.25)*IT_1162;
    const ccomplex_t IT_1164 = 2*conj(IT_0528);
    const ccomplex_t IT_1165 = 2*conj(IT_0529);
    const ccomplex_t IT_1166 = IT_1164 + IT_1165;
    const ccomplex_t IT_1167 = IT_0905 + IT_1166;
    const ccomplex_t IT_1168 = m_N_3*IT_0004;
    const ccomplex_t IT_1169 = IT_0004*IT_1147;
    const ccomplex_t IT_1170 = (-2)*IT_0332 + (-0.5)*IT_1169;
    const ccomplex_t IT_1171 = IT_0331 + IT_1170;
    const ccomplex_t IT_1172 = IT_1168*IT_1171;
    const ccomplex_t IT_1173 = (-2)*IT_1172;
    const ccomplex_t IT_1174 = (-0.25)*IT_1173;
    const ccomplex_t IT_1175 = s_14*s_34*IT_0049;
    const ccomplex_t IT_1176 = (-2)*IT_1175;
    const ccomplex_t IT_1177 = IT_1158 + IT_1176;
    const ccomplex_t IT_1178 = IT_0001*IT_1177;
    const ccomplex_t IT_1179 = (-0.666666666666667)*IT_0297 + (
      -0.166666666666667)*IT_1178;
    const ccomplex_t IT_1180 = s_13 + IT_1179;
    const ccomplex_t IT_1181 = m_N_3*IT_1180;
    const ccomplex_t IT_1182 = (-6)*IT_1181;
    const ccomplex_t IT_1183 = (-0.25)*IT_1182;
    const ccomplex_t IT_1184 = (-2)*conj(IT_0431);
    const ccomplex_t IT_1185 = IT_0000*IT_0004;
    const ccomplex_t IT_1186 = (-0.25)*IT_0714;
    const ccomplex_t IT_1187 = IT_0267 + IT_1186;
    const ccomplex_t IT_1188 = IT_1185*IT_1187;
    const ccomplex_t IT_1189 = (-4)*IT_1188;
    const ccomplex_t IT_1190 = (-0.25)*IT_1189;
    const ccomplex_t IT_1191 = (-2)*conj(IT_0429);
    const ccomplex_t IT_1192 = 9*IT_0676;
    const ccomplex_t IT_1193 = (-0.25)*IT_1192;
    const ccomplex_t IT_1194 = 4*conj(IT_0287);
    const ccomplex_t IT_1195 = 4*conj(IT_0288);
    const ccomplex_t IT_1196 = IT_1194 + IT_1195;
    const ccomplex_t IT_1197 = (-2)*IT_0716;
    const ccomplex_t IT_1198 = IT_0715 + IT_1197;
    const ccomplex_t IT_1199 = IT_0004*IT_1198;
    const ccomplex_t IT_1200 = 0.5*IT_1199;
    const ccomplex_t IT_1201 = (-2)*IT_0716 + 2*IT_0717 + 2*IT_0718 + (-2)
      *IT_0719;
    const ccomplex_t IT_1202 = IT_0715 + IT_1201;
    const ccomplex_t IT_1203 = IT_0004*IT_1202;
    const ccomplex_t IT_1204 = (-0.5)*IT_1203;
    const ccomplex_t IT_1205 = s_14*s_23*IT_0049;
    const ccomplex_t IT_1206 = s_13*s_24*IT_0049;
    const ccomplex_t IT_1207 = IT_1205 + IT_1206;
    const ccomplex_t IT_1208 = s_12*s_34*IT_0049;
    const ccomplex_t IT_1209 = -IT_1208;
    const ccomplex_t IT_1210 = IT_1207 + IT_1209;
    const ccomplex_t IT_1211 = IT_0320*IT_1210;
    const ccomplex_t IT_1212 = 0.5*IT_1211;
    const ccomplex_t IT_1213 = IT_0458 + IT_0460 + IT_1200 + IT_1204 + IT_1212;
    const ccomplex_t IT_1214 = (-0.25)*IT_1213;
    const ccomplex_t IT_1215 = (-4)*conj(IT_0202);
    const ccomplex_t IT_1216 = (-4)*conj(IT_0205);
    const ccomplex_t IT_1217 = IT_1215 + IT_1216;
    const ccomplex_t IT_1218 = s_12*IT_0004*IT_0267;
    const ccomplex_t IT_1219 = -IT_1218;
    const ccomplex_t IT_1220 = IT_0715 + IT_0717 + IT_0718;
    const ccomplex_t IT_1221 = (-2)*IT_0716 + -IT_0719;
    const ccomplex_t IT_1222 = IT_1220 + IT_1221;
    const ccomplex_t IT_1223 = IT_0004*IT_1222;
    const ccomplex_t IT_1224 = (-0.5)*IT_1223;
    const ccomplex_t IT_1225 = -IT_0694;
    const ccomplex_t IT_1226 = IT_0458 + IT_0868 + IT_1219 + IT_1224 + IT_1225;
    const ccomplex_t IT_1227 = (-0.25)*IT_1226;
    const ccomplex_t IT_1228 = (-8)*conj(IT_0033);
    const ccomplex_t IT_1229 = (-8)*conj(IT_0078);
    const ccomplex_t IT_1230 = IT_1228 + IT_1229;
    const ccomplex_t IT_1231 = (-0.5)*IT_0595;
    const ccomplex_t IT_1232 = s_23*IT_0251;
    const ccomplex_t IT_1233 = IT_1068 + IT_1232;
    const ccomplex_t IT_1234 = s_12*s_14*s_34;
    const ccomplex_t IT_1235 = -IT_1069 + -IT_1234;
    const ccomplex_t IT_1236 = IT_1233 + IT_1235;
    const ccomplex_t IT_1237 = IT_0004*IT_1236;
    const ccomplex_t IT_1238 = (-0.25)*IT_1237;
    const ccomplex_t IT_1239 = (-0.5)*IT_0566;
    const ccomplex_t IT_1240 = (-0.5)*IT_0596;
    const ccomplex_t IT_1241 = 0.5*IT_0598;
    const ccomplex_t IT_1242 = IT_0564 + IT_1231 + IT_1238 + IT_1239 + IT_1240
       + IT_1241;
    const ccomplex_t IT_1243 = (-0.25)*IT_1242;
    const ccomplex_t IT_1244 = IT_1064*IT_1065 + (-0.25)*IT_1077*IT_1080 + (
      -0.25)*IT_1082*IT_1085 + (-0.25)*IT_1096*IT_1099 + (-0.25)*IT_1108*IT_1111
       + IT_0128*IT_1119 + IT_0245*IT_1127 + IT_0502*IT_1135 + IT_0713*IT_1143 +
       IT_1146*IT_1156 + IT_1157*IT_1163 + IT_1167*IT_1174 + IT_0927*IT_1183 +
       IT_1184*IT_1190 + IT_1191*IT_1193 + IT_1196*IT_1214 + IT_1217*IT_1227 +
       IT_1230*IT_1243;
    const ccomplex_t IT_1245 = IT_0004*IT_0929;
    const ccomplex_t IT_1246 = -IT_1245;
    const ccomplex_t IT_1247 = IT_1061 + IT_1246;
    const ccomplex_t IT_1248 = (-4)*IT_0435;
    const ccomplex_t IT_1249 = (-0.25)*IT_1248;
    const ccomplex_t IT_1250 = IT_0437 + IT_1249;
    const ccomplex_t IT_1251 = 0.5*IT_1237;
    const ccomplex_t IT_1252 = IT_0595 + IT_0597 + IT_0598 + IT_1066 + IT_1251;
    const ccomplex_t IT_1253 = s_14*s_34*IT_0000;
    const ccomplex_t IT_1254 = 0.5*IT_1091 + (-0.5)*IT_1102 + -IT_1253;
    const ccomplex_t IT_1255 = IT_1087 + IT_1254;
    const ccomplex_t IT_1256 = IT_0004*IT_1255;
    const ccomplex_t IT_1257 = (-2)*IT_1256;
    const ccomplex_t IT_1258 = (-0.5)*IT_1257;
    const ccomplex_t IT_1259 = IT_1101 + IT_1258;
    const ccomplex_t IT_1260 = 4*IT_0074;
    const ccomplex_t IT_1261 = IT_1051 + IT_1260;
    const ccomplex_t IT_1262 = (-0.5)*IT_1069 + 0.5*IT_1070 + -IT_1234;
    const ccomplex_t IT_1263 = IT_1232 + IT_1262;
    const ccomplex_t IT_1264 = IT_0004*IT_1263;
    const ccomplex_t IT_1265 = 2*IT_1264;
    const ccomplex_t IT_1266 = (-0.5)*IT_1265;
    const ccomplex_t IT_1267 = IT_1067 + IT_1266;
    const ccomplex_t IT_1268 = (-4)*IT_0081;
    const ccomplex_t IT_1269 = (-4)*IT_0226;
    const ccomplex_t IT_1270 = IT_1268 + IT_1269;
    const ccomplex_t IT_1271 = IT_0001*IT_1147;
    const ccomplex_t IT_1272 = 0.666666666666667*IT_0703 + (-0.166666666666667
      )*IT_1271;
    const ccomplex_t IT_1273 = s_23 + IT_1272;
    const ccomplex_t IT_1274 = m_N_3*IT_1273;
    const ccomplex_t IT_1275 = (-6)*IT_1274;
    const ccomplex_t IT_1276 = (-2)*IT_0546;
    const ccomplex_t IT_1277 = (-2)*IT_0548;
    const ccomplex_t IT_1278 = IT_0867 + IT_1276 + IT_1277;
    const ccomplex_t IT_1279 = (-0.666666666666667)*IT_0330 + (
      -0.166666666666667)*IT_1151;
    const ccomplex_t IT_1280 = s_23 + IT_1279;
    const ccomplex_t IT_1281 = m_N_3*IT_1280;
    const ccomplex_t IT_1282 = (-6)*IT_1281;
    const ccomplex_t IT_1283 = (-0.25)*IT_1282;
    const ccomplex_t IT_1284 = (-2)*IT_0523;
    const ccomplex_t IT_1285 = IT_0004*IT_1158;
    const ccomplex_t IT_1286 = (-2)*IT_0299 + (-0.5)*IT_1285;
    const ccomplex_t IT_1287 = IT_0298 + IT_1286;
    const ccomplex_t IT_1288 = IT_1168*IT_1287;
    const ccomplex_t IT_1289 = (-2)*IT_1288;
    const ccomplex_t IT_1290 = (-0.25)*IT_1289;
    const ccomplex_t IT_1291 = (-0.666666666666667)*IT_1122 +
       0.166666666666667*IT_1178;
    const ccomplex_t IT_1292 = s_13 + IT_1291;
    const ccomplex_t IT_1293 = m_N_3*IT_1292;
    const ccomplex_t IT_1294 = 6*IT_1293;
    const ccomplex_t IT_1295 = (-0.25)*IT_1294;
    const ccomplex_t IT_1296 = 2*IT_0528;
    const ccomplex_t IT_1297 = 2*IT_0529;
    const ccomplex_t IT_1298 = IT_1296 + IT_1297;
    const ccomplex_t IT_1299 = (-10)*s_23;
    const ccomplex_t IT_1300 = IT_1271 + IT_1299;
    const ccomplex_t IT_1301 = m_N_3*IT_1300;
    const ccomplex_t IT_1302 = (-0.25)*IT_1301;
    const ccomplex_t IT_1303 = 4*IT_0287;
    const ccomplex_t IT_1304 = 4*IT_0288;
    const ccomplex_t IT_1305 = IT_1303 + IT_1304;
    const ccomplex_t IT_1306 = (-4)*IT_0202;
    const ccomplex_t IT_1307 = (-4)*IT_0205;
    const ccomplex_t IT_1308 = IT_1306 + IT_1307;
    const ccomplex_t IT_1309 = -IT_0868;
    const ccomplex_t IT_1310 = IT_0458 + IT_0694 + IT_1219 + IT_1224 + IT_1309;
    const ccomplex_t IT_1311 = (-0.25)*IT_1310;
    const ccomplex_t IT_1312 = (-8)*IT_0033;
    const ccomplex_t IT_1313 = (-8)*IT_0078;
    const ccomplex_t IT_1314 = IT_1312 + IT_1313;
    const ccomplex_t IT_1315 = (-0.5)*IT_0613;
    const ccomplex_t IT_1316 = (-0.25)*IT_1094;
    const ccomplex_t IT_1317 = (-0.5)*IT_0492;
    const ccomplex_t IT_1318 = 0.5*IT_0614;
    const ccomplex_t IT_1319 = (-0.5)*IT_0616;
    const ccomplex_t IT_1320 = IT_0490 + IT_1315 + IT_1316 + IT_1317 + IT_1318
       + IT_1319;
    const ccomplex_t IT_1321 = (-0.25)*IT_1320;
    const ccomplex_t IT_1322 = IT_0125*IT_1119 + IT_0242*IT_1127 + IT_0500
      *IT_1135 + IT_1031*IT_1143 + IT_1247*IT_1250 + (-0.25)*IT_0970*IT_1252 + (
      -0.25)*IT_1259*IT_1261 + (-0.25)*IT_1267*IT_1270 + (-0.25)*IT_0890*IT_1275
       + IT_1278*IT_1283 + IT_1284*IT_1290 + IT_0911*IT_1295 + IT_1298*IT_1302 +
       IT_1214*IT_1305 + IT_1308*IT_1311 + IT_1314*IT_1321;
    const ccomplex_t IT_1323 = IT_0890 + IT_1298;
    const ccomplex_t IT_1324 = 4*IT_0131;
    const ccomplex_t IT_1325 = IT_1007 + IT_1324;
    const ccomplex_t IT_1326 = 4*IT_0275;
    const ccomplex_t IT_1327 = 4*IT_0278;
    const ccomplex_t IT_1328 = IT_1326 + IT_1327;
    const ccomplex_t IT_1329 = (-8)*IT_0129;
    const ccomplex_t IT_1330 = (-8)*IT_0131;
    const ccomplex_t IT_1331 = IT_1329 + IT_1330;
    const ccomplex_t IT_1332 = (-2)*IT_0435;
    const ccomplex_t IT_1333 = (-2)*IT_0437;
    const ccomplex_t IT_1334 = IT_1332 + IT_1333;
    const ccomplex_t IT_1335 = (-4)*IT_0209;
    const ccomplex_t IT_1336 = (-0.25)*IT_1335;
    const ccomplex_t IT_1337 = 0.25*IT_0876;
    const ccomplex_t IT_1338 = (-4)*IT_0130;
    const ccomplex_t IT_1339 = (-0.25)*IT_1338;
    const ccomplex_t IT_1340 = IT_0190*IT_1118 + IT_0229*IT_1126 + IT_0909
      *IT_1134 + IT_1029*IT_1142 + IT_1282*IT_1284 + IT_1278*IT_1289 + IT_0911
      *IT_1301 + IT_1294*IT_1323 + IT_1252*IT_1325 + IT_1213*IT_1328 + IT_1320
      *IT_1331 + IT_1192*IT_1334 + (-4)*IT_1310*(IT_0208 + IT_1336) + 4*IT_1259*
      (IT_0166 + IT_1337) + (-4)*IT_1267*(IT_0239 + IT_1339);
    const ccomplex_t IT_1341 = -IT_1077;
    const ccomplex_t IT_1342 = 0.666666666666667*IT_0319 + (-0.166666666666667
      )*IT_1159;
    const ccomplex_t IT_1343 = s_13 + IT_1342;
    const ccomplex_t IT_1344 = m_N_3*IT_1343;
    const ccomplex_t IT_1345 = (-6)*IT_1344;
    const ccomplex_t IT_1346 = (-0.25)*IT_1345;
    const ccomplex_t IT_1347 = IT_0219*IT_1096 + IT_0190*IT_1119 + IT_0229
      *IT_1127 + IT_0909*IT_1135 + IT_1029*IT_1143 + IT_0867*IT_1163 + IT_0911
      *IT_1174 + IT_1156*IT_1284 + IT_1183*IT_1323 + IT_1214*IT_1328 + IT_1243
      *IT_1331 + IT_1190*IT_1334 + IT_1226*(IT_0208 + IT_1336) + IT_1108*
      (IT_0166 + IT_1339) + (IT_0239 + IT_1337)*IT_1341 + (IT_1276 + IT_1277)
      *IT_1346;
    const ccomplex_t IT_1348 = 4*conj(IT_0130);
    const ccomplex_t IT_1349 = IT_0193*IT_1119;
    const ccomplex_t IT_1350 = IT_0232*IT_1127;
    const ccomplex_t IT_1351 = IT_1041*IT_1135;
    const ccomplex_t IT_1352 = IT_0711*IT_1143;
    const ccomplex_t IT_1353 = (-2)*conj(IT_0435);
    const ccomplex_t IT_1354 = (-2)*conj(IT_0437);
    const ccomplex_t IT_1355 = IT_1353 + IT_1354;
    const ccomplex_t IT_1356 = 4*conj(IT_0275);
    const ccomplex_t IT_1357 = 4*conj(IT_0278);
    const ccomplex_t IT_1358 = IT_1356 + IT_1357;
    const ccomplex_t IT_1359 = IT_1214*IT_1358;
    const ccomplex_t IT_1360 = (-4)*conj(IT_0208);
    const ccomplex_t IT_1361 = (-4)*conj(IT_0209);
    const ccomplex_t IT_1362 = IT_1360 + IT_1361;
    const ccomplex_t IT_1363 = (-8)*conj(IT_0129);
    const ccomplex_t IT_1364 = (-8)*conj(IT_0131);
    const ccomplex_t IT_1365 = IT_1363 + IT_1364;
    const ccomplex_t IT_1366 = (-4)*conj(IT_0130);
    const ccomplex_t IT_1367 = (-0.25)*IT_1366;
    const ccomplex_t IT_1368 = 2*conj(IT_0429)*IT_1082 + conj(IT_0431)*IT_1247
       + -IT_0684*IT_1252 + IT_1157*IT_1283 + IT_1146*IT_1290 + IT_1167*IT_1295 
      + IT_0927*IT_1302 + -IT_1259*(conj(IT_0166) + 0.25*IT_1348) + IT_1349 +
       IT_1350 + IT_1351 + IT_1352 + IT_1193*IT_1355 + IT_1359 + IT_1311*IT_1362
       + IT_1321*IT_1365 + IT_1267*(conj(IT_0239) + IT_1367);
    const ccomplex_t IT_1369 = conj(IT_0429)*IT_1064 + 2*conj(IT_0431)*IT_1082
       + IT_0218*IT_1096 + IT_1156*IT_1157 + IT_0926*IT_1163 + IT_0927*IT_1174 +
       IT_1167*IT_1183 + (IT_1144 + IT_1145)*IT_1346 + IT_1341*(conj(IT_0239) +
       0.25*IT_1348) + IT_1349 + IT_1350 + IT_1351 + IT_1352 + IT_1190*IT_1355 +
       IT_1359 + IT_1227*IT_1362 + IT_1243*IT_1365 + IT_1108*(conj(IT_0166) +
       IT_1367);
    const ccomplex_t IT_1370 = s_34*IT_0000*IT_0004;
    const ccomplex_t IT_1371 = IT_0198*IT_1370;
    const ccomplex_t IT_1372 = IT_0481 + IT_0607;
    const ccomplex_t IT_1373 = IT_0483 + IT_0609;
    const ccomplex_t IT_1374 = conj(IT_0081)*IT_0467;
    const ccomplex_t IT_1375 = IT_0081*IT_0469;
    const ccomplex_t IT_1376 = conj(IT_0205)*IT_0364;
    const ccomplex_t IT_1377 = conj(IT_0288)*IT_0364;
    const ccomplex_t IT_1378 = IT_0205*IT_0366;
    const ccomplex_t IT_1379 = IT_0288*IT_0366;
    const ccomplex_t IT_1380 = conj(IT_0205) + conj(IT_0288);
    const ccomplex_t IT_1381 = IT_0205 + IT_0288;
    const ccomplex_t IT_1382 = conj(IT_0130)*IT_0481;
    const ccomplex_t IT_1383 = IT_0130*IT_0483;
    const ccomplex_t IT_1384 = conj(IT_0081)*IT_0488;
    const ccomplex_t IT_1385 = IT_0081*IT_0501;
    const ccomplex_t IT_1386 = conj(IT_0130)*IT_0607;
    const ccomplex_t IT_1387 = IT_0130*IT_0609;
    const ccomplex_t IT_1388 = conj(IT_0081)*IT_0618;
    const ccomplex_t IT_1389 = IT_0081*IT_0611;
    const ccomplex_t IT_1390 = IT_0278*IT_0526;
    const ccomplex_t IT_1391 = conj(IT_0278)*IT_0524;
    const ccomplex_t IT_1392 = (IT_0210 + conj(IT_0275) + conj(IT_0278))
      *IT_0399 + (IT_0211 + IT_0275 + IT_0278)*IT_0401 + conj(IT_0166)*IT_1372 +
       IT_0166*IT_1373 + IT_1374 + IT_1375 + IT_1376 + IT_1377 + IT_1378 +
       IT_1379 + IT_0552*IT_1380 + IT_0374*IT_1381 + IT_1382 + IT_1383 + IT_1384
       + IT_1385 + IT_1386 + IT_1387 + IT_1388 + IT_1389 + IT_1390 + IT_1391;
    const ccomplex_t IT_1393 = IT_0264*IT_1370;
    const ccomplex_t IT_1394 = (IT_0210 + conj(IT_0275) + conj(IT_0278))
      *IT_0554 + (IT_0211 + IT_0275 + IT_0278)*IT_0555 + conj(IT_0239)*IT_1372 +
       IT_0239*IT_1373 + IT_1374 + IT_1375 + IT_1376 + IT_1377 + IT_1378 +
       IT_1379 + IT_0600*IT_1380 + IT_0355*IT_1381 + IT_1382 + IT_1383 + IT_1384
       + IT_1385 + IT_1386 + IT_1387 + IT_1388 + IT_1389 + IT_1390 + IT_1391;
    const ccomplex_t IT_1395 = (-12)*IT_0078;
    const ccomplex_t IT_1396 = (-12)*IT_0033 + IT_1395;
    const ccomplex_t IT_1397 = IT_0968*IT_1396;
    const ccomplex_t IT_1398 = IT_0354*IT_0554;
    const ccomplex_t IT_1399 = (-4)*IT_0310;
    const ccomplex_t IT_1400 = (-4)*IT_0313;
    const ccomplex_t IT_1401 = IT_1399 + IT_1400;
    const ccomplex_t IT_1402 = IT_0682*IT_1401;
    const ccomplex_t IT_1403 = IT_0532*IT_0956;
    const ccomplex_t IT_1404 = IT_0354*IT_0555;
    const ccomplex_t IT_1405 = IT_0373*IT_0401;
    const ccomplex_t IT_1406 = IT_0780*IT_0926;
    const ccomplex_t IT_1407 = IT_0544*IT_0896;
    const ccomplex_t IT_1408 = IT_1404 + IT_1405 + IT_1406 + IT_1407;
    const ccomplex_t IT_1409 = (-12)*conj(IT_0033);
    const ccomplex_t IT_1410 = IT_0772*IT_0972;
    const ccomplex_t IT_1411 = (-4)*conj(IT_0310);
    const ccomplex_t IT_1412 = (-4)*conj(IT_0313);
    const ccomplex_t IT_1413 = IT_1411 + IT_1412;
    const ccomplex_t IT_1414 = IT_0682*IT_1413;
    const ccomplex_t IT_1415 = IT_0536*IT_0956;
    const ccomplex_t IT_1416 = (-6)*conj(IT_0074);
    const ccomplex_t IT_1417 = 6*conj(IT_0226);
    const ccomplex_t IT_1418 = IT_0646*IT_1417;
    const ccomplex_t IT_1419 = IT_0889*IT_0905;
    const ccomplex_t IT_1420 = IT_0760*IT_0926;
    const ccomplex_t IT_1421 = (-6)*conj(IT_0226);
    const ccomplex_t IT_1422 = IT_0544*IT_0875;
    const ccomplex_t IT_1423 = 6*conj(IT_0074);
    const ccomplex_t IT_1424 = IT_0650*IT_1423;
    const ccomplex_t IT_1425 = IT_0544*IT_0862;
    const ccomplex_t IT_1426 = IT_0846*IT_0863;
    const ccomplex_t IT_1427 = IT_0356 + IT_0411;
    const ccomplex_t IT_1428 = 8*IT_0123;
    const ccomplex_t IT_1429 = (-8)*IT_0120;
    const ccomplex_t IT_1430 = IT_1428 + IT_1429;
    const ccomplex_t IT_1431 = 8*IT_0546;
    const ccomplex_t IT_1432 = (-8)*IT_0310;
    const ccomplex_t IT_1433 = (-8)*IT_0313;
    const ccomplex_t IT_1434 = IT_1431 + IT_1432 + IT_1433;
    const ccomplex_t IT_1435 = IT_0350*IT_1004;
    const ccomplex_t IT_1436 = 16*IT_1435;
    const ccomplex_t IT_1437 = IT_0896*IT_0997;
    const ccomplex_t IT_1438 = (-4)*IT_0033;
    const ccomplex_t IT_1439 = (-4)*IT_0078;
    const ccomplex_t IT_1440 = IT_1438 + IT_1439;
    const ccomplex_t IT_1441 = IT_0946*IT_1440;
    const ccomplex_t IT_1442 = IT_1008 + IT_1023;
    const ccomplex_t IT_1443 = IT_0341*IT_1442;
    const ccomplex_t IT_1444 = IT_0875*IT_0917;
    const ccomplex_t IT_1445 = (-2)*IT_0288;
    const ccomplex_t IT_1446 = IT_1011 + IT_1445;
    const ccomplex_t IT_1447 = IT_0857*IT_1446;
    const ccomplex_t IT_1448 = (-2)*IT_0278;
    const ccomplex_t IT_1449 = IT_1025 + IT_1448;
    const ccomplex_t IT_1450 = IT_0862*IT_1449;
    const ccomplex_t IT_1451 = IT_0500*IT_0907 + IT_0889*IT_0912 + IT_0525
      *IT_0923 + IT_0932*(IT_1032 + IT_1034) + IT_0904*IT_1427 + IT_0938*IT_1430
       + IT_0727*IT_1434 + IT_1436 + IT_1437 + IT_1441 + IT_1443 + IT_1444 +
       IT_1447 + IT_1450;
    const ccomplex_t IT_1452 = 8*IT_0548;
    const ccomplex_t IT_1453 = IT_0500*IT_0907 + IT_0889*IT_0912 + IT_0525
      *IT_0923 + IT_0932*(IT_1032 + IT_1033 + IT_1034) + IT_0794*IT_1036 +
       IT_0904*IT_1427 + IT_0938*IT_1430 + IT_1436 + IT_1437 + IT_1441 + IT_1443
       + IT_1444 + IT_1447 + IT_1450 + IT_0727*(IT_1431 + IT_1432 + IT_1433 +
       IT_1452);
    const ccomplex_t IT_1454 = IT_0682*IT_1325;
    const ccomplex_t IT_1455 = IT_0780*IT_0997;
    const ccomplex_t IT_1456 = IT_0785*IT_1442;
    const ccomplex_t IT_1457 = IT_0524 + IT_0554;
    const ccomplex_t IT_1458 = IT_0701*IT_1457;
    const ccomplex_t IT_1459 = IT_0760*IT_1046;
    const ccomplex_t IT_1460 = IT_0531 + IT_0949;
    const ccomplex_t IT_1461 = IT_0771*IT_1460;
    const ccomplex_t IT_1462 = IT_0794*IT_1035;
    const ccomplex_t IT_1463 = IT_0740*IT_1430;
    const ccomplex_t IT_1464 = conj(IT_0546)*(IT_0308*IT_0500 + IT_0709
      *IT_1029 + IT_0329*IT_1031 + IT_0727*IT_1034 + IT_0692*IT_1427 + IT_0750
      *IT_1434 + IT_1454 + IT_1455 + IT_1456 + IT_1458 + IT_1459 + IT_1461 +
       IT_1462 + IT_1463);
    const ccomplex_t IT_1465 = IT_0308*IT_0500 + (-16)*IT_0523*IT_1004 +
       IT_0709*IT_1029 + IT_0329*IT_1031 + 8*IT_0727*(IT_0529 + 0.125*IT_1032 +
       0.125*IT_1034) + IT_0692*IT_1427 + IT_0750*(IT_1431 + IT_1432 + IT_1433 +
       IT_1452) + IT_1454 + IT_1455 + IT_1456 + IT_1458 + IT_1459 + IT_1461 +
       IT_1462 + IT_1463;
    const ccomplex_t IT_1466 = IT_0483*IT_0785;
    const ccomplex_t IT_1467 = IT_0401*IT_0771;
    const ccomplex_t IT_1468 = IT_1466 + IT_1467;
    const ccomplex_t IT_1469 = IT_0692*IT_0955;
    const ccomplex_t IT_1470 = IT_0536*IT_0701;
    const ccomplex_t IT_1471 = IT_0709*IT_1024;
    const ccomplex_t IT_1472 = IT_0740*IT_0953;
    const ccomplex_t IT_1473 = 4*IT_0226;
    const ccomplex_t IT_1474 = (-8)*IT_0435;
    const ccomplex_t IT_1475 = (-8)*IT_0437;
    const ccomplex_t IT_1476 = (-4)*IT_0074;
    const ccomplex_t IT_1477 = (-4)*IT_0222;
    const ccomplex_t IT_1478 = (-2)*IT_0429;
    const ccomplex_t IT_1479 = (-2)*IT_0431;
    const ccomplex_t IT_1480 = IT_0242*IT_1126;
    const ccomplex_t IT_1481 = IT_1213*IT_1305;
    const ccomplex_t IT_1482 = IT_0125*IT_1118 + IT_0500*IT_1134 + IT_1031
      *IT_1142 + IT_0911*IT_1182 + IT_1155*IT_1278 + IT_1162*IT_1284 + IT_1226
      *IT_1308 + IT_1242*IT_1314 + IT_1173*IT_1323 + IT_1096*IT_1440 + IT_1077*
      (IT_1051 + IT_1473) + IT_1082*(IT_1474 + IT_1475) + IT_1108*(IT_1268 +
       IT_1476) + IT_1064*IT_1477 + IT_1192*IT_1478 + IT_1189*IT_1479 + IT_1480 
      + IT_1481;
    const ccomplex_t IT_1483 = (-0.25)*IT_0125*IT_1118 + (-0.25)*IT_0500
      *IT_1134 + (-0.25)*IT_1031*IT_1142 + IT_1247*IT_1250 + (-0.25)*IT_0970
      *IT_1252 + (-0.25)*IT_1259*IT_1261 + (-0.25)*IT_1267*IT_1270 + (-0.25)
      *IT_0890*IT_1275 + (-0.25)*IT_1278*IT_1282 + (-0.25)*IT_1284*IT_1289 + (
      -0.25)*IT_0911*IT_1294 + (-0.25)*IT_1298*IT_1301 + (-0.25)*IT_1308*IT_1310
       + (-0.25)*IT_1314*IT_1320 + (-0.25)*IT_1480 + (-0.25)*IT_1481;
    const ccomplex_t IT_1484 = IT_0128*IT_1118;
    const ccomplex_t IT_1485 = IT_0245*IT_1126;
    const ccomplex_t IT_1486 = IT_0502*IT_1134;
    const ccomplex_t IT_1487 = IT_0713*IT_1142;
    const ccomplex_t IT_1488 = IT_1196*IT_1213;
    const ccomplex_t IT_1489 = IT_1077*IT_1080 + IT_1082*IT_1085 + IT_1096
      *IT_1099 + IT_1108*IT_1111 + IT_1146*IT_1155 + IT_1157*IT_1162 + IT_1167
      *IT_1173 + IT_0927*IT_1182 + IT_1184*IT_1189 + IT_1191*IT_1192 + IT_1217
      *IT_1226 + IT_1230*IT_1242 + IT_1484 + IT_1485 + IT_1486 + IT_1487 +
       IT_1488;
    const ccomplex_t IT_1490 = IT_1146*IT_1282;
    const ccomplex_t IT_1491 = IT_1157*IT_1289;
    const ccomplex_t IT_1492 = IT_0905*IT_1275;
    const ccomplex_t IT_1493 = IT_0927*IT_1294;
    const ccomplex_t IT_1494 = IT_1166*IT_1301;
    const ccomplex_t IT_1495 = IT_1217*IT_1310;
    const ccomplex_t IT_1496 = (-4)*conj(IT_0226);
    const ccomplex_t IT_1497 = IT_1110 + IT_1496;
    const ccomplex_t IT_1498 = IT_1267*IT_1497;
    const ccomplex_t IT_1499 = IT_1043*IT_1252;
    const ccomplex_t IT_1500 = IT_1230*IT_1320;
    const ccomplex_t IT_1501 = 4*conj(IT_0074);
    const ccomplex_t IT_1502 = IT_1078 + IT_1501;
    const ccomplex_t IT_1503 = IT_1259*IT_1502;
    const ccomplex_t IT_1504 = 16*conj(IT_0350);
    const ccomplex_t IT_1505 = IT_1004*IT_1504;
    const ccomplex_t IT_1506 = IT_0693*IT_0904;
    const ccomplex_t IT_1507 = IT_0783*IT_0896;
    const ccomplex_t IT_1508 = IT_0946*IT_1099;
    const ccomplex_t IT_1509 = IT_0341*IT_0788;
    const ccomplex_t IT_1510 = IT_0502*IT_0907;
    const ccomplex_t IT_1511 = IT_0889*IT_1040;
    const ccomplex_t IT_1512 = IT_0527*IT_0923;
    const ccomplex_t IT_1513 = IT_0875*IT_1044;
    const ccomplex_t IT_1514 = (-2)*conj(IT_0288);
    const ccomplex_t IT_1515 = IT_1024 + IT_1514;
    const ccomplex_t IT_1516 = IT_0857*IT_1515;
    const ccomplex_t IT_1517 = (-2)*conj(IT_0278);
    const ccomplex_t IT_1518 = IT_1026 + IT_1517;
    const ccomplex_t IT_1519 = IT_0862*IT_1518;
    const ccomplex_t IT_1520 = IT_0743*IT_0938;
    const ccomplex_t IT_1521 = IT_0727*IT_0754;
    const ccomplex_t IT_1522 = IT_0194 + IT_0199*(IT_0125*IT_0206 + IT_0128
      *IT_0207 + IT_0190*IT_0210 + IT_0193*IT_0211) + IT_0246 + (IT_0232*
      (IT_0129 + IT_0130 + IT_0131 + IT_0239) + IT_0229*(conj(IT_0129) + conj
      (IT_0130) + conj(IT_0131) + conj(IT_0239)) + (conj(IT_0081) + conj(IT_0226
      ) + IT_0233)*IT_0242 + (IT_0081 + IT_0226 + IT_0234)*IT_0245)*IT_0256 +
       IT_0259*((conj(IT_0074) + conj(IT_0081))*IT_0190 + (IT_0074 + IT_0081)
      *IT_0193 + IT_0229*IT_0233 + IT_0232*IT_0234 + IT_0218*IT_0242 + IT_0219
      *IT_0245 + IT_0125*IT_0260 + IT_0128*IT_0261) + (IT_0210*IT_0229 + IT_0211
      *IT_0232 + IT_0206*IT_0242 + IT_0207*IT_0245)*IT_0265 + IT_0271*(IT_0279
      *IT_0282 + IT_0283*IT_0286 + IT_0289*IT_0292 + IT_0293*IT_0296) + IT_0308*
      (IT_0210*IT_0315 + IT_0211*IT_0318) + (IT_0279*IT_0315 + IT_0283*IT_0318)
      *IT_0329 + IT_0341*(conj(IT_0209)*IT_0351 + IT_0209*IT_0352) + IT_0354*
      (IT_0078*IT_0355 + IT_0357 + IT_0239*IT_0358 + IT_0359 + IT_0361 + IT_0363
       + IT_0365 + IT_0367 + IT_0369 + IT_0371) + IT_0373*(conj(IT_0166)*IT_0356
       + IT_0357 + IT_0166*IT_0358 + IT_0359 + IT_0361 + IT_0363 + IT_0365 +
       IT_0367 + IT_0369 + IT_0371 + IT_0078*IT_0374) + (IT_0384 + IT_0394)*
      (IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0400 + IT_0402 + IT_0403 +
       IT_0404 + IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409 + IT_0410 +
       IT_0412 + IT_0414) + (IT_0421 + IT_0428)*IT_0457 + (IT_0464 + IT_0466)*
      (conj(IT_0208)*IT_0360 + IT_0208*IT_0362 + IT_0468 + IT_0470 + IT_0471 +
       IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 + IT_0478 +
       IT_0479 + IT_0480 + IT_0482 + IT_0484) + IT_0487*(IT_0210*IT_0467 +
       IT_0209*IT_0469 + conj(IT_0202)*IT_0481 + conj(IT_0209)*IT_0488) +
       IT_0511 + IT_0518*IT_0556 + (conj(IT_0074)*IT_0554 + IT_0074*IT_0555)
      *IT_0562 + (conj(IT_0226)*IT_0481 + IT_0226*IT_0483 + IT_0239*IT_0502 +
       IT_0503 + IT_0504 + IT_0505 + IT_0506 + IT_0507 + IT_0508 + IT_0510 +
       conj(IT_0205)*IT_0554 + IT_0205*IT_0555)*IT_0574 + IT_0588*(IT_0239
      *IT_0370 + IT_0396 + IT_0398 + IT_0403 + IT_0405 + IT_0407 + IT_0408 +
       IT_0412 + IT_0414 + IT_0589 + IT_0590 + IT_0591 + IT_0592 + IT_0593 +
       IT_0594) + IT_0612 + IT_0619 + (IT_0621 + IT_0623)*(IT_0624 + IT_0625 +
       IT_0626 + IT_0627 + IT_0628 + IT_0629 + IT_0630 + IT_0631 + IT_0632 +
       IT_0633 + IT_0634 + IT_0635 + IT_0636 + IT_0637 + IT_0638 + IT_0639) +
       conj(IT_0033)*IT_0641*IT_0642 + (IT_0646 + IT_0650)*(IT_0652 + IT_0654) +
       IT_0663*(IT_0665 + IT_0667 + IT_0669 + IT_0671 + IT_0673 + IT_0675) + 4
      *IT_0548*IT_0796 + 4*IT_0546*IT_0797 + IT_0800*(IT_0210*IT_0801 + IT_0211
      *IT_0802) + IT_0809*(IT_0811 + conj(IT_0208)*IT_0812 + IT_0813 + IT_0208
      *IT_0814 + IT_0815 + IT_0817 + IT_0819 + IT_0821 + IT_0823 + IT_0824 +
       IT_0825 + IT_0827 + IT_0829) + 4*IT_0830*(0.25*conj(IT_0120)*IT_0280 +
       0.25*conj(IT_0184)*IT_0290 + IT_0123*IT_0832 + IT_0188*IT_0834 + 0.25
      *conj(IT_0123)*IT_0835 + 0.25*conj(IT_0188)*IT_0836) + 4*IT_0838*(IT_0293
      *IT_0832 + IT_0283*IT_0834 + IT_0279*(IT_0188 + 0.25*IT_0836) + IT_0289
      *IT_0839) + 4*IT_0840*IT_0843 + conj(IT_0287)*IT_0844*IT_0845 + conj
      (IT_0288)*(IT_0315*IT_0709 + IT_0844*IT_0845 + IT_0846*IT_0849 + IT_0543
      *IT_0857 + IT_0351*IT_0862) + conj(IT_0275)*(IT_0845*IT_0846 + IT_0844
      *IT_0847 + IT_0543*IT_0862) + conj(IT_0278)*(IT_0845*IT_0846 + IT_0844
      *IT_0849 + IT_0351*IT_0857 + IT_0543*IT_0862) + IT_0288*(IT_0318*IT_0709 +
       IT_0544*IT_0857 + IT_0352*IT_0862 + (-4)*IT_0844*(conj(IT_0288) + (-0.25)
      *IT_0863) + (-4)*IT_0846*(conj(IT_0278) + (-0.25)*IT_0864)) + conj(IT_0130
      )*(IT_0760*IT_0867 + IT_0543*IT_0875 + IT_0562*IT_0876 + IT_0650*IT_0877 +
       IT_0646*IT_0878 + IT_0880*(IT_0881 + IT_0882) + IT_0889*IT_0890) +
       IT_0078*(IT_0318*IT_0692 + IT_0352*IT_0896 + IT_0641*IT_0897 + IT_0904
      *IT_0905) + (IT_0210*IT_0890 + IT_0211*IT_0905)*IT_0907 + conj(IT_0205)*
      (IT_0341*IT_0543 + IT_0785*IT_0867 + IT_0487*IT_0909 + IT_0800*(IT_0826 +
       IT_0910) + IT_0907*IT_0911) + conj(IT_0166)*(IT_0543*IT_0875 + IT_0650
      *IT_0877 + IT_0880*IT_0882 + IT_0562*IT_0912 + IT_0914*IT_0917 + IT_0911
      *IT_0923) + IT_0205*(IT_0483*IT_0487 + IT_0341*IT_0544 + IT_0800*(IT_0810 
      + IT_0828) + IT_0785*IT_0926 + IT_0907*IT_0927) + 8*conj(IT_0542)*IT_0950 
      + IT_0542*(IT_0341*IT_0483 + IT_0609*IT_0857 + IT_0772*IT_0923 + IT_0897
      *IT_0946 + IT_0938*IT_0953 + IT_0727*IT_0954 + IT_0904*IT_0955) + IT_0956*
      (IT_0078*IT_0761 + IT_0239*IT_0955 + IT_0957 + IT_0958 + IT_0960 + IT_0961
       + IT_0962 + IT_0963 + IT_0965 + IT_0966) + conj(IT_0078)*IT_0973 +
       IT_0972*(conj(IT_0166)*IT_0948 + IT_0166*IT_0955 + IT_0957 + IT_0958 +
       IT_0960 + IT_0961 + IT_0962 + IT_0963 + IT_0965 + IT_0966 + IT_0078
      *IT_0974) + (IT_0976 + IT_0978)*(IT_0979 + IT_0980 + IT_0981 + IT_0982 +
       IT_0983 + IT_0984 + IT_0985 + IT_0987) + IT_0994*IT_0999 + IT_0914*(conj
      (IT_0074)*IT_0949 + IT_1000 + IT_1001) + (conj(IT_0226)*IT_0532 + IT_1000 
      + IT_1001)*IT_1003 + conj(IT_0350)*IT_1009 + IT_0081*IT_1010 + conj
      (IT_0310)*(IT_0532*IT_0701 + IT_0399*IT_0771 + IT_0481*IT_0785 + IT_0740
      *IT_0941 + IT_0750*IT_0947 + IT_0692*IT_0948 + IT_0709*IT_1011) + 8*conj
      (IT_0313)*IT_1012 + IT_1018*IT_1027 + conj(IT_0523)*IT_1037 + IT_0350
      *IT_1039 + IT_0523*IT_1045 + 6*conj(IT_0239)*IT_1050 + 6*conj(IT_0081)
      *IT_1052 + (-4)*IT_0225*IT_1244 + (-4)*conj(IT_0435)*IT_1322 + conj
      (IT_0431)*IT_1340 + (-4)*conj(IT_0429)*IT_1347 + (-4)*IT_0431*IT_1368 + (
      -4)*IT_0429*IT_1369 + IT_1371*IT_1392 + IT_1393*IT_1394 + conj(IT_0129)*
      (IT_0373*IT_0399 + IT_0663*IT_0664 + IT_0801*IT_0809 + IT_0780*IT_0867 +
       IT_0543*IT_0896 + IT_0949*IT_0972 + IT_0641*IT_1007 + IT_1397 + IT_1398 +
       IT_1402 + IT_1403) + conj(IT_0131)*(IT_0373*IT_0399 + IT_0780*IT_0867 +
       IT_0543*IT_0896 + IT_0904*IT_0911 + IT_0949*IT_0972 + IT_0641*IT_1325 +
       IT_1397 + IT_1398 + IT_1402 + IT_1403) + IT_0131*(IT_0641*IT_0683 +
       IT_0904*IT_0927 + IT_1408 + IT_0968*IT_1409 + IT_1410 + IT_1414 + IT_1415
      ) + IT_0129*(IT_0663*IT_0666 + IT_0802*IT_0809 + IT_1408 + IT_0968*((-12)
      *conj(IT_0078) + IT_1409) + IT_1410 + IT_1414 + IT_1415) + IT_0239*
      (IT_0366*IT_0562 + IT_0762*IT_1003 + IT_0880*IT_1416 + IT_1418 + IT_1419 +
       IT_1420) + IT_0166*(IT_0762*IT_0914 + IT_0923*IT_0927 + IT_0562*IT_1040 +
       IT_0880*IT_1421 + IT_1422 + IT_1424) + IT_0130*(IT_1418 + IT_1419 +
       IT_1420 + IT_0880*(IT_1416 + IT_1421) + IT_1422 + IT_1424) + IT_0275*
      (IT_1425 + IT_1426) + IT_0278*(IT_0352*IT_0857 + IT_0844*IT_0864 + IT_1425
       + IT_1426) + conj(IT_0528)*IT_1451 + conj(IT_0529)*IT_1453 + IT_1464 +
       conj(IT_0548)*IT_1465 + IT_0310*(IT_1468 + IT_1469 + IT_1470 + IT_1471 +
       IT_1472) + IT_0313*(IT_0750*IT_0954 + IT_1468 + IT_1469 + IT_1470 +
       IT_1471 + IT_1472) + IT_1065*IT_1482 + (-4)*conj(IT_0437)*IT_1483 +
       IT_0222*IT_1489 + (IT_0435 + IT_0437)*(IT_1484 + IT_1485 + IT_1486 +
       IT_1487 + IT_1488 + IT_1490 + IT_1491 + IT_1492 + IT_1493 + IT_1494 +
       IT_1495 + IT_1498 + IT_1499 + IT_1500 + IT_1503) + IT_0528*(IT_0729
      *IT_0932 + IT_1505 + IT_1506 + IT_1507 + IT_1508 + IT_1509 + IT_1510 +
       IT_1511 + IT_1512 + IT_1513 + IT_1516 + IT_1519 + IT_1520 + IT_1521) +
       IT_0529*(IT_0730*IT_0932 + IT_1505 + IT_1506 + IT_1507 + IT_1508 +
       IT_1509 + IT_1510 + IT_1511 + IT_1512 + IT_1513 + IT_1516 + IT_1519 +
       IT_1520 + IT_1521);
    return create_ccomplex_return(IT_1522);
}

