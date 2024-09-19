#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_W_W(
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0000, -1);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = IT_0000*IT_0002;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0008 + -IT_0010 + (-2)
      *IT_0012);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_1, 2);
    const ccomplex_t IT_0016 = pow(m_C_2, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0015 + -IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0004*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = sin(alpha);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0022
      *IT_0023 + -IT_0024*IT_0025);
    const ccomplex_t IT_0027 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0030 = IT_0025*IT_0029;
    const ccomplex_t IT_0031 = IT_0002*(IT_0028 + -IT_0030);
    const ccomplex_t IT_0032 = 1.4142135623731*e_em*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0026*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0015 + -IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0022
      *IT_0024 + IT_0023*IT_0025);
    const ccomplex_t IT_0041 = IT_0022*IT_0029;
    const ccomplex_t IT_0042 = IT_0025*IT_0027;
    const ccomplex_t IT_0043 = IT_0002*(IT_0041 + IT_0042);
    const ccomplex_t IT_0044 = 1.4142135623731*e_em*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0040*IT_0046;
    const ccomplex_t IT_0048 = IT_0039*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0003;
    const ccomplex_t IT_0051 = m_C_2*IT_0050;
    const ccomplex_t IT_0052 = IT_0019*IT_0051;
    const ccomplex_t IT_0053 = m_C_1*IT_0050;
    const ccomplex_t IT_0054 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = IT_0009*IT_0054;
    const ccomplex_t IT_0057 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0058 = IT_0009*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + -IT_0056 + (-2)
      *IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0053*IT_0062;
    const ccomplex_t IT_0064 = IT_0038 + IT_0049 + IT_0052 + IT_0063;
    const ccomplex_t IT_0065 = pow(m_W, -2);
    const ccomplex_t IT_0066 = pow(s_13, 2);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = pow(s_14, 2);
    const ccomplex_t IT_0069 = IT_0065*IT_0068;
    const ccomplex_t IT_0070 = IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = pow(m_W, -4);
    const ccomplex_t IT_0072 = s_13*s_14*s_34*IT_0071;
    const ccomplex_t IT_0073 = -IT_0015 + -IT_0072;
    const ccomplex_t IT_0074 = IT_0070 + IT_0073;
    const ccomplex_t IT_0075 = m_C_2*IT_0074;
    const ccomplex_t IT_0076 = 2*IT_0075;
    const ccomplex_t IT_0077 = IT_0064*IT_0076;
    const ccomplex_t IT_0078 = IT_0003*IT_0019;
    const ccomplex_t IT_0079 = s_13*s_14*s_23*s_34*IT_0071;
    const ccomplex_t IT_0080 = s_24*s_34*IT_0066*IT_0071;
    const ccomplex_t IT_0081 = s_23*IT_0015;
    const ccomplex_t IT_0082 = s_12*s_13;
    const ccomplex_t IT_0083 = s_23*IT_0065*IT_0066;
    const ccomplex_t IT_0084 = (-2)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082 + IT_0084;
    const ccomplex_t IT_0086 = pow(s_34, 2);
    const ccomplex_t IT_0087 = s_12*s_13*IT_0071*IT_0086;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0079 + IT_0080 + IT_0081 + IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = s_24*s_34*IT_0015*IT_0065;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = 2*IT_0092;
    const ccomplex_t IT_0094 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0095 = IT_0002*IT_0094;
    const ccomplex_t IT_0096 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0097 = IT_0002*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 0.5*IT_0098);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0102 = IT_0002*IT_0101;
    const ccomplex_t IT_0103 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0104 = IT_0002*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0102 + 0.5*IT_0105);
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0100*IT_0107;
    const ccomplex_t IT_0109 = pow(m_W, 2);
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0016 + IT_0109 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0114 = IT_0002*IT_0113;
    const ccomplex_t IT_0115 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0116 = IT_0002*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + 0.5*IT_0117);
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0121 = IT_0002*IT_0120;
    const ccomplex_t IT_0122 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0123 = IT_0002*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0119*IT_0126;
    const ccomplex_t IT_0128 = cpow((-2)*s_23 + IT_0016 + IT_0109 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0134 = IT_0002*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0139 = IT_0002*IT_0138;
    const ccomplex_t IT_0140 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0141 = IT_0002*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 0.5*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0137*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_23 + IT_0016 + IT_0109 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0150 = IT_0002*IT_0149;
    const ccomplex_t IT_0151 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0152 = IT_0002*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + 0.5*IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0157 = IT_0002*IT_0156;
    const ccomplex_t IT_0158 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0159 = IT_0002*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + 0.5*IT_0160);
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0155*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_23 + IT_0016 + IT_0109 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = -IT_0020 + (-2)*IT_0112 + (-2)*IT_0130 + (-2)
      *IT_0148 + (-2)*IT_0166;
    const ccomplex_t IT_0168 = pow(s_23, 2);
    const ccomplex_t IT_0169 = s_14*s_34*IT_0071*IT_0168;
    const ccomplex_t IT_0170 = s_13*s_23*s_24*s_34*IT_0071;
    const ccomplex_t IT_0171 = s_13*IT_0016;
    const ccomplex_t IT_0172 = s_12*s_23;
    const ccomplex_t IT_0173 = s_13*IT_0065*IT_0168;
    const ccomplex_t IT_0174 = (-2)*IT_0173;
    const ccomplex_t IT_0175 = IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = s_12*s_23*IT_0071*IT_0086;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0169 + IT_0170 + IT_0171 + IT_0175 + IT_0177;
    const ccomplex_t IT_0179 = s_14*s_34*IT_0016*IT_0065;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0181;
    const ccomplex_t IT_0183 = s_34*m_C_2*IT_0065;
    const ccomplex_t IT_0184 = s_13*s_34*IT_0065;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = s_14 + IT_0185;
    const ccomplex_t IT_0187 = IT_0183*IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = IT_0065*IT_0086;
    const ccomplex_t IT_0190 = s_12*IT_0109;
    const ccomplex_t IT_0191 = s_13*s_23;
    const ccomplex_t IT_0192 = (-2)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = IT_0065*IT_0193;
    const ccomplex_t IT_0195 = (IT_0109 + -IT_0189)*(s_12 + (-0.5)*IT_0194);
    const ccomplex_t IT_0196 = (-2)*IT_0195;
    const ccomplex_t IT_0197 = 2*IT_0196;
    const ccomplex_t IT_0198 = IT_0004*IT_0062;
    const ccomplex_t IT_0199 = s_34*m_C_1*m_C_2*IT_0065;
    const ccomplex_t IT_0200 = IT_0186*IT_0199;
    const ccomplex_t IT_0201 = 2*IT_0200;
    const ccomplex_t IT_0202 = conj(IT_0198)*IT_0201;
    const ccomplex_t IT_0203 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0204 = IT_0002*IT_0203;
    const ccomplex_t IT_0205 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0206 = IT_0002*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (-0.5)*IT_0207);
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0211 = IT_0002*IT_0210;
    const ccomplex_t IT_0212 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0213 = IT_0002*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0211 + (-0.5)*IT_0214);
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0209*IT_0216;
    const ccomplex_t IT_0218 = IT_0110*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0221 = IT_0002*IT_0220;
    const ccomplex_t IT_0222 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0223 = IT_0002*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + (-0.5)*IT_0224);
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0228 = IT_0002*IT_0227;
    const ccomplex_t IT_0229 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0230 = IT_0002*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + (-0.5)*IT_0231);
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0226*IT_0233;
    const ccomplex_t IT_0235 = IT_0128*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0238 = IT_0002*IT_0237;
    const ccomplex_t IT_0239 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0240 = IT_0002*IT_0239;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*(IT_0238 + (-0.5)*IT_0241);
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0245 = IT_0002*IT_0244;
    const ccomplex_t IT_0246 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0247 = IT_0002*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*(IT_0245 + (-0.5)*IT_0248);
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0243*IT_0250;
    const ccomplex_t IT_0252 = IT_0146*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0255 = IT_0002*IT_0254;
    const ccomplex_t IT_0256 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0257 = IT_0002*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*(IT_0255 + (-0.5)*IT_0258);
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0262 = IT_0002*IT_0261;
    const ccomplex_t IT_0263 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0264 = IT_0002*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*(IT_0262 + (-0.5)*IT_0265);
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = IT_0260*IT_0267;
    const ccomplex_t IT_0269 = IT_0164*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = IT_0198 + (-2)*IT_0219 + (-2)*IT_0236 + (-2)
      *IT_0253 + (-2)*IT_0270;
    const ccomplex_t IT_0272 = s_23*s_34*IT_0065;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = s_24 + IT_0273;
    const ccomplex_t IT_0275 = IT_0199*IT_0274;
    const ccomplex_t IT_0276 = 2*IT_0275;
    const ccomplex_t IT_0277 = conj(IT_0271)*IT_0276;
    const ccomplex_t IT_0278 = IT_0155*IT_0267;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0280 = IT_0164*IT_0278*IT_0279;
    const ccomplex_t IT_0281 = IT_0119*IT_0233;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0283 = IT_0128*IT_0281*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0285 = IT_0100*IT_0216;
    const ccomplex_t IT_0286 = IT_0110*IT_0284*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0288 = IT_0137*IT_0250;
    const ccomplex_t IT_0289 = IT_0146*IT_0287*IT_0288;
    const ccomplex_t IT_0290 = m_C_2*IT_0166;
    const ccomplex_t IT_0291 = m_C_2*IT_0112;
    const ccomplex_t IT_0292 = m_C_2*IT_0130;
    const ccomplex_t IT_0293 = m_C_2*IT_0148;
    const ccomplex_t IT_0294 = IT_0280 + IT_0283 + IT_0286 + IT_0289 + 
      -IT_0290 + -IT_0291 + -IT_0292 + -IT_0293;
    const ccomplex_t IT_0295 = s_14*s_34*IT_0065;
    const ccomplex_t IT_0296 = s_13*IT_0109;
    const ccomplex_t IT_0297 = s_14*s_34;
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = IT_0296 + IT_0298;
    const ccomplex_t IT_0300 = IT_0065*IT_0299;
    const ccomplex_t IT_0301 = -IT_0295 + (-0.25)*IT_0300;
    const ccomplex_t IT_0302 = s_13 + IT_0301;
    const ccomplex_t IT_0303 = m_C_2*IT_0302;
    const ccomplex_t IT_0304 = (-4)*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = IT_0219 + IT_0236 + IT_0253 + IT_0270;
    const ccomplex_t IT_0307 = m_C_1*m_C_2;
    const ccomplex_t IT_0308 = -IT_0086;
    const ccomplex_t IT_0309 = pow(m_W, 4);
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = IT_0065*IT_0310;
    const ccomplex_t IT_0312 = -IT_0109 + 0.25*IT_0311;
    const ccomplex_t IT_0313 = IT_0189 + IT_0312;
    const ccomplex_t IT_0314 = IT_0307*IT_0313;
    const ccomplex_t IT_0315 = 4*IT_0314;
    const ccomplex_t IT_0316 = 2*IT_0315;
    const ccomplex_t IT_0317 = IT_0003*IT_0062;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = -IT_0189;
    const ccomplex_t IT_0320 = IT_0109 + IT_0319;
    const ccomplex_t IT_0321 = IT_0307*IT_0320;
    const ccomplex_t IT_0322 = (-6)*IT_0321;
    const ccomplex_t IT_0323 = IT_0112 + IT_0130 + IT_0148 + IT_0166;
    const ccomplex_t IT_0324 = s_13*s_24*s_34*IT_0065;
    const ccomplex_t IT_0325 = s_12*IT_0065*IT_0086;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = s_12*IT_0309;
    const ccomplex_t IT_0328 = s_14*s_23*s_34;
    const ccomplex_t IT_0329 = s_13*s_24*s_34;
    const ccomplex_t IT_0330 = IT_0327 + IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = s_13*s_23*IT_0109;
    const ccomplex_t IT_0332 = s_12*IT_0086;
    const ccomplex_t IT_0333 = (-2)*IT_0331 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0330 + IT_0333;
    const ccomplex_t IT_0335 = IT_0065*IT_0334;
    const ccomplex_t IT_0336 = (-0.5)*IT_0335;
    const ccomplex_t IT_0337 = s_14*s_23*s_34*IT_0065;
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = IT_0190 + IT_0324 + IT_0326 + IT_0336 + IT_0338;
    const ccomplex_t IT_0340 = (-4)*IT_0339;
    const ccomplex_t IT_0341 = U_d2*V_Wp1;
    const ccomplex_t IT_0342 = IT_0022*IT_0341;
    const ccomplex_t IT_0343 = V_u1*U_Wm2;
    const ccomplex_t IT_0344 = IT_0025*IT_0343;
    const ccomplex_t IT_0345 = IT_0002*(IT_0342 + IT_0344);
    const ccomplex_t IT_0346 = 1.4142135623731*e_em*IT_0345;
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = IT_0040*IT_0348;
    const ccomplex_t IT_0350 = IT_0039*IT_0349;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*IT_0350;
    const ccomplex_t IT_0352 = IT_0019*IT_0053;
    const ccomplex_t IT_0353 = IT_0051*IT_0062;
    const ccomplex_t IT_0354 = IT_0022*IT_0343;
    const ccomplex_t IT_0355 = IT_0025*IT_0341;
    const ccomplex_t IT_0356 = IT_0002*(IT_0354 + -IT_0355);
    const ccomplex_t IT_0357 = 1.4142135623731*e_em*IT_0356;
    const ccomplex_t IT_0358 = 0.5*IT_0357;
    const ccomplex_t IT_0359 = (0 + _Complex_I*1)*IT_0358;
    const ccomplex_t IT_0360 = IT_0026*IT_0359;
    const ccomplex_t IT_0361 = IT_0036*IT_0360;
    const ccomplex_t IT_0362 = (0 + _Complex_I*1)*IT_0361;
    const ccomplex_t IT_0363 = IT_0351 + -IT_0352 + -IT_0353 + IT_0362;
    const ccomplex_t IT_0364 = s_34*m_C_1*IT_0065;
    const ccomplex_t IT_0365 = IT_0274*IT_0364;
    const ccomplex_t IT_0366 = (-2)*IT_0365;
    const ccomplex_t IT_0367 = (-2)*IT_0003;
    const ccomplex_t IT_0368 = IT_0019*IT_0367;
    const ccomplex_t IT_0369 = s_34*IT_0065;
    const ccomplex_t IT_0370 = s_14*s_23;
    const ccomplex_t IT_0371 = s_13*s_24;
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = s_13*s_23*s_34*IT_0065;
    const ccomplex_t IT_0374 = (-2)*IT_0373;
    const ccomplex_t IT_0375 = IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = IT_0369*IT_0375;
    const ccomplex_t IT_0377 = (-2)*IT_0376;
    const ccomplex_t IT_0378 = IT_0144*IT_0243;
    const ccomplex_t IT_0379 = IT_0146*IT_0287*IT_0378;
    const ccomplex_t IT_0380 = IT_0162*IT_0260;
    const ccomplex_t IT_0381 = IT_0164*IT_0279*IT_0380;
    const ccomplex_t IT_0382 = IT_0107*IT_0209;
    const ccomplex_t IT_0383 = IT_0110*IT_0284*IT_0382;
    const ccomplex_t IT_0384 = IT_0126*IT_0226;
    const ccomplex_t IT_0385 = IT_0128*IT_0282*IT_0384;
    const ccomplex_t IT_0386 = m_C_2*IT_0270;
    const ccomplex_t IT_0387 = m_C_2*IT_0219;
    const ccomplex_t IT_0388 = m_C_2*IT_0236;
    const ccomplex_t IT_0389 = m_C_2*IT_0253;
    const ccomplex_t IT_0390 = IT_0379 + IT_0381 + IT_0383 + IT_0385 + 
      -IT_0386 + -IT_0387 + -IT_0388 + -IT_0389;
    const ccomplex_t IT_0391 = s_24*s_34*IT_0065;
    const ccomplex_t IT_0392 = s_23*IT_0109;
    const ccomplex_t IT_0393 = s_24*s_34;
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = s_23*IT_0065*IT_0086;
    const ccomplex_t IT_0396 = (-2)*IT_0395;
    const ccomplex_t IT_0397 = IT_0394 + IT_0396;
    const ccomplex_t IT_0398 = IT_0065*IT_0397;
    const ccomplex_t IT_0399 = -IT_0391 + 0.5*IT_0398;
    const ccomplex_t IT_0400 = s_23 + IT_0399;
    const ccomplex_t IT_0401 = m_C_1*IT_0400;
    const ccomplex_t IT_0402 = 2*IT_0401;
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = conj(IT_0021)*IT_0093 + conj(IT_0167)*IT_0182 +
       conj(IT_0064)*IT_0188 + conj(IT_0078)*IT_0197 + IT_0202 + IT_0277 + conj
      (IT_0294)*IT_0305 + conj(IT_0306)*IT_0316 + conj(IT_0318)*IT_0322 + conj
      (IT_0323)*IT_0340 + conj(IT_0363)*IT_0366 + conj(IT_0368)*IT_0377 + conj
      (IT_0390)*IT_0403;
    const ccomplex_t IT_0405 = IT_0062*IT_0367;
    const ccomplex_t IT_0406 = -IT_0405;
    const ccomplex_t IT_0407 = IT_0093*conj(IT_0198) + conj(IT_0021)*IT_0201 +
       IT_0182*conj(IT_0271) + conj(IT_0167)*IT_0276 + IT_0197*conj(IT_0318) +
       conj(IT_0078)*IT_0322 + IT_0316*conj(IT_0323) + conj(IT_0306)*IT_0340 +
       IT_0188*conj(IT_0363) + conj(IT_0064)*IT_0366 + IT_0305*conj(IT_0390) +
       conj(IT_0294)*IT_0403 + IT_0377*conj(IT_0406);
    const ccomplex_t IT_0408 = IT_0071*IT_0086;
    const ccomplex_t IT_0409 = 0.5*IT_0408;
    const ccomplex_t IT_0410 = 1 + IT_0409;
    const ccomplex_t IT_0411 = IT_0307*IT_0410;
    const ccomplex_t IT_0412 = (-4)*IT_0411;
    const ccomplex_t IT_0413 = s_13*s_14*s_34*IT_0065;
    const ccomplex_t IT_0414 = IT_0015*IT_0109;
    const ccomplex_t IT_0415 = 2*IT_0068 + (-2)*IT_0413 + -IT_0414;
    const ccomplex_t IT_0416 = IT_0066 + IT_0415;
    const ccomplex_t IT_0417 = IT_0065*IT_0416;
    const ccomplex_t IT_0418 = -IT_0015 + (-0.5)*IT_0417;
    const ccomplex_t IT_0419 = IT_0067 + IT_0418;
    const ccomplex_t IT_0420 = m_C_2*IT_0419;
    const ccomplex_t IT_0421 = (-2)*IT_0420;
    const ccomplex_t IT_0422 = 2*IT_0421;
    const ccomplex_t IT_0423 = IT_0016*IT_0109;
    const ccomplex_t IT_0424 = -IT_0423;
    const ccomplex_t IT_0425 = IT_0168 + IT_0424;
    const ccomplex_t IT_0426 = IT_0065*IT_0425;
    const ccomplex_t IT_0427 = IT_0065*IT_0168;
    const ccomplex_t IT_0428 = -IT_0016;
    const ccomplex_t IT_0429 = IT_0427 + IT_0428;
    const ccomplex_t IT_0430 = (-4)*IT_0429;
    const ccomplex_t IT_0431 = IT_0426 + IT_0430;
    const ccomplex_t IT_0432 = m_C_1*IT_0431;
    const ccomplex_t IT_0433 = 2*IT_0432;
    const ccomplex_t IT_0434 = (-2)*IT_0297;
    const ccomplex_t IT_0435 = IT_0296 + IT_0434;
    const ccomplex_t IT_0436 = IT_0065*IT_0435;
    const ccomplex_t IT_0437 = s_13*IT_0309;
    const ccomplex_t IT_0438 = s_14*s_34*IT_0109;
    const ccomplex_t IT_0439 = (-2)*IT_0438;
    const ccomplex_t IT_0440 = IT_0437 + IT_0439;
    const ccomplex_t IT_0441 = IT_0071*IT_0440;
    const ccomplex_t IT_0442 = (-0.666666666666667)*IT_0436 +
       0.166666666666667*IT_0441;
    const ccomplex_t IT_0443 = s_13 + IT_0442;
    const ccomplex_t IT_0444 = m_C_2*IT_0443;
    const ccomplex_t IT_0445 = 6*IT_0444;
    const ccomplex_t IT_0446 = 2*IT_0445;
    const ccomplex_t IT_0447 = s_23*IT_0309;
    const ccomplex_t IT_0448 = IT_0071*IT_0447;
    const ccomplex_t IT_0449 = (-10)*s_23;
    const ccomplex_t IT_0450 = IT_0448 + IT_0449;
    const ccomplex_t IT_0451 = m_C_1*IT_0450;
    const ccomplex_t IT_0452 = 2*IT_0451;
    const ccomplex_t IT_0453 = IT_0065*IT_0190;
    const ccomplex_t IT_0454 = (-0.5)*IT_0453;
    const ccomplex_t IT_0455 = s_34*IT_0071;
    const ccomplex_t IT_0456 = s_12*s_34;
    const ccomplex_t IT_0457 = IT_0370 + IT_0456;
    const ccomplex_t IT_0458 = -IT_0371;
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = IT_0455*IT_0459;
    const ccomplex_t IT_0461 = 0.25*IT_0460;
    const ccomplex_t IT_0462 = s_12 + IT_0454 + IT_0461;
    const ccomplex_t IT_0463 = 8*IT_0462;
    const ccomplex_t IT_0464 = s_14*s_24;
    const ccomplex_t IT_0465 = (-2)*IT_0464;
    const ccomplex_t IT_0466 = IT_0190 + IT_0465;
    const ccomplex_t IT_0467 = IT_0065*IT_0466;
    const ccomplex_t IT_0468 = s_14*s_24*IT_0109;
    const ccomplex_t IT_0469 = 4*IT_0328 + (-2)*IT_0331 + (-2)*IT_0468;
    const ccomplex_t IT_0470 = IT_0327 + IT_0469;
    const ccomplex_t IT_0471 = IT_0071*IT_0470;
    const ccomplex_t IT_0472 = 8*s_12 + (-4)*IT_0194 + (-4)*IT_0467 + 2*IT_0471;
    const ccomplex_t IT_0473 = 1.125*IT_0307;
    const ccomplex_t IT_0474 = (-16)*IT_0473;
    const ccomplex_t IT_0475 = s_13*s_23*IT_0065;
    const ccomplex_t IT_0476 = -IT_0191 + 2*IT_0337 + (-2)*IT_0464;
    const ccomplex_t IT_0477 = IT_0190 + IT_0476;
    const ccomplex_t IT_0478 = IT_0065*IT_0477;
    const ccomplex_t IT_0479 = -IT_0475 + (-0.5)*IT_0478;
    const ccomplex_t IT_0480 = s_12 + IT_0479;
    const ccomplex_t IT_0481 = m_C_2*IT_0480;
    const ccomplex_t IT_0482 = (-2)*IT_0481;
    const ccomplex_t IT_0483 = (-2)*IT_0482;
    const ccomplex_t IT_0484 = -IT_0191;
    const ccomplex_t IT_0485 = IT_0190 + IT_0484;
    const ccomplex_t IT_0486 = IT_0065*IT_0485;
    const ccomplex_t IT_0487 = -IT_0475 + (-0.25)*IT_0486;
    const ccomplex_t IT_0488 = s_12 + IT_0487;
    const ccomplex_t IT_0489 = m_C_1*IT_0488;
    const ccomplex_t IT_0490 = (-4)*IT_0489;
    const ccomplex_t IT_0491 = (-2)*IT_0490;
    const ccomplex_t IT_0492 = IT_0065*IT_0296;
    const ccomplex_t IT_0493 = s_14*IT_0109;
    const ccomplex_t IT_0494 = IT_0455*IT_0493;
    const ccomplex_t IT_0495 = (-0.5)*IT_0492 + 0.25*IT_0494;
    const ccomplex_t IT_0496 = s_13 + IT_0495;
    const ccomplex_t IT_0497 = m_C_2*IT_0496;
    const ccomplex_t IT_0498 = 4*IT_0497;
    const ccomplex_t IT_0499 = (-2)*IT_0498;
    const ccomplex_t IT_0500 = IT_0065*IT_0392;
    const ccomplex_t IT_0501 = s_24*IT_0109;
    const ccomplex_t IT_0502 = s_23*s_34;
    const ccomplex_t IT_0503 = (-2)*IT_0502;
    const ccomplex_t IT_0504 = IT_0501 + IT_0503;
    const ccomplex_t IT_0505 = IT_0455*IT_0504;
    const ccomplex_t IT_0506 = (-0.5)*IT_0500 + (-0.25)*IT_0505;
    const ccomplex_t IT_0507 = s_23 + IT_0506;
    const ccomplex_t IT_0508 = m_C_1*IT_0507;
    const ccomplex_t IT_0509 = (-4)*IT_0508;
    const ccomplex_t IT_0510 = (-2)*IT_0509;
    const ccomplex_t IT_0511 = conj(IT_0078)*IT_0305 + conj(IT_0318)*IT_0403 +
       conj(IT_0363)*IT_0412 + conj(IT_0021)*IT_0422 + conj(IT_0271)*IT_0433 +
       conj(IT_0323)*IT_0446 + conj(IT_0306)*IT_0452 + conj(IT_0064)*IT_0463 +
       conj(IT_0294)*IT_0472 + conj(IT_0390)*IT_0474 + conj(IT_0167)*IT_0483 +
       conj(IT_0198)*IT_0491 + conj(IT_0368)*IT_0499 + conj(IT_0406)*IT_0510;
    const ccomplex_t IT_0512 = IT_0305*conj(IT_0318) + conj(IT_0078)*IT_0403 +
       conj(IT_0064)*IT_0412 + conj(IT_0198)*IT_0422 + conj(IT_0167)*IT_0433 +
       conj(IT_0306)*IT_0446 + conj(IT_0323)*IT_0452 + conj(IT_0363)*IT_0463 +
       conj(IT_0390)*IT_0472 + conj(IT_0294)*IT_0474 + conj(IT_0271)*IT_0483 +
       conj(IT_0021)*IT_0491 + conj(IT_0406)*IT_0499 + conj(IT_0368)*IT_0510;
    const ccomplex_t IT_0513 = (-0.5)*IT_0436 + 0.5*IT_0492 + (-0.5)*IT_0494;
    const ccomplex_t IT_0514 = s_13 + IT_0513;
    const ccomplex_t IT_0515 = m_C_2*IT_0514;
    const ccomplex_t IT_0516 = 2*IT_0515;
    const ccomplex_t IT_0517 = 2*IT_0516;
    const ccomplex_t IT_0518 = s_34*IT_0109;
    const ccomplex_t IT_0519 = IT_0455*IT_0518;
    const ccomplex_t IT_0520 = (-2)*IT_0086;
    const ccomplex_t IT_0521 = IT_0309 + IT_0520;
    const ccomplex_t IT_0522 = IT_0065*IT_0521;
    const ccomplex_t IT_0523 = (-0.333333333333333)*IT_0519 + (
      -0.333333333333333)*IT_0522;
    const ccomplex_t IT_0524 = IT_0109 + IT_0523;
    const ccomplex_t IT_0525 = IT_0307*IT_0524;
    const ccomplex_t IT_0526 = 3*IT_0525;
    const ccomplex_t IT_0527 = 2*IT_0526;
    const ccomplex_t IT_0528 = s_12*s_23*IT_0109;
    const ccomplex_t IT_0529 = s_14*s_23*s_24;
    const ccomplex_t IT_0530 = s_13*IT_0016*IT_0109;
    const ccomplex_t IT_0531 = s_14*s_34*IT_0016;
    const ccomplex_t IT_0532 = (-2)*IT_0529 + -IT_0530 + 2*IT_0531;
    const ccomplex_t IT_0533 = IT_0528 + IT_0532;
    const ccomplex_t IT_0534 = IT_0065*IT_0533;
    const ccomplex_t IT_0535 = (-0.5)*IT_0534;
    const ccomplex_t IT_0536 = -IT_0171;
    const ccomplex_t IT_0537 = IT_0172 + IT_0535 + IT_0536;
    const ccomplex_t IT_0538 = 4*IT_0537;
    const ccomplex_t IT_0539 = (-2)*IT_0331;
    const ccomplex_t IT_0540 = IT_0327 + IT_0539;
    const ccomplex_t IT_0541 = IT_0065*IT_0540;
    const ccomplex_t IT_0542 = 0.5*IT_0541;
    const ccomplex_t IT_0543 = 2*IT_0328 + 2*IT_0329 + (-2)*IT_0331 + (-2)
      *IT_0332;
    const ccomplex_t IT_0544 = IT_0327 + IT_0543;
    const ccomplex_t IT_0545 = IT_0065*IT_0544;
    const ccomplex_t IT_0546 = (-0.5)*IT_0545;
    const ccomplex_t IT_0547 = s_14*s_23*IT_0109;
    const ccomplex_t IT_0548 = s_13*s_24*IT_0109;
    const ccomplex_t IT_0549 = IT_0547 + IT_0548;
    const ccomplex_t IT_0550 = s_12*s_34*IT_0109;
    const ccomplex_t IT_0551 = -IT_0550;
    const ccomplex_t IT_0552 = IT_0549 + IT_0551;
    const ccomplex_t IT_0553 = IT_0455*IT_0552;
    const ccomplex_t IT_0554 = 0.5*IT_0553;
    const ccomplex_t IT_0555 = IT_0190 + IT_0192 + IT_0542 + IT_0546 + IT_0554;
    const ccomplex_t IT_0556 = 4*IT_0555;
    const ccomplex_t IT_0557 = s_12*s_13*IT_0109;
    const ccomplex_t IT_0558 = s_23*IT_0068;
    const ccomplex_t IT_0559 = s_12*s_14*s_34;
    const ccomplex_t IT_0560 = s_23*IT_0015*IT_0109;
    const ccomplex_t IT_0561 = 2*IT_0558 + (-2)*IT_0559 + -IT_0560;
    const ccomplex_t IT_0562 = IT_0557 + IT_0561;
    const ccomplex_t IT_0563 = IT_0065*IT_0562;
    const ccomplex_t IT_0564 = (-0.5)*IT_0563;
    const ccomplex_t IT_0565 = -IT_0081;
    const ccomplex_t IT_0566 = IT_0082 + IT_0564 + IT_0565;
    const ccomplex_t IT_0567 = (-4)*IT_0566;
    const ccomplex_t IT_0568 = pow(m_W, 6);
    const ccomplex_t IT_0569 = s_12*IT_0568;
    const ccomplex_t IT_0570 = s_14*s_24*IT_0309;
    const ccomplex_t IT_0571 = (-2)*IT_0570;
    const ccomplex_t IT_0572 = IT_0569 + IT_0571;
    const ccomplex_t IT_0573 = IT_0071*IT_0572;
    const ccomplex_t IT_0574 = IT_0065*IT_0470;
    const ccomplex_t IT_0575 = (-4)*IT_0190 + 16*IT_0191 + (-2)*IT_0573 + 4
      *IT_0574;
    const ccomplex_t IT_0576 = IT_0455*IT_0501;
    const ccomplex_t IT_0577 = (-2)*IT_0393;
    const ccomplex_t IT_0578 = IT_0392 + IT_0577;
    const ccomplex_t IT_0579 = IT_0065*IT_0578;
    const ccomplex_t IT_0580 = 0.5*IT_0500 + (-0.5)*IT_0576 + (-0.5)*IT_0579;
    const ccomplex_t IT_0581 = s_23 + IT_0580;
    const ccomplex_t IT_0582 = m_C_1*IT_0581;
    const ccomplex_t IT_0583 = 2*IT_0582;
    const ccomplex_t IT_0584 = (-2)*IT_0583;
    const ccomplex_t IT_0585 = IT_0109*IT_0307;
    const ccomplex_t IT_0586 = 9*IT_0585;
    const ccomplex_t IT_0587 = (-2)*IT_0586;
    const ccomplex_t IT_0588 = IT_0316*conj(IT_0318) + conj(IT_0078)*IT_0340 +
       conj(IT_0294)*IT_0446 + conj(IT_0390)*IT_0452 + conj(IT_0064)*IT_0517 +
       conj(IT_0406)*IT_0527 + conj(IT_0167)*IT_0538 + conj(IT_0368)*IT_0556 +
       conj(IT_0021)*IT_0567 + conj(IT_0323)*IT_0575 + conj(IT_0363)*IT_0584 +
       conj(IT_0306)*IT_0587;
    const ccomplex_t IT_0589 = conj(IT_0078)*IT_0316 + conj(IT_0318)*IT_0340 +
       conj(IT_0390)*IT_0446 + conj(IT_0294)*IT_0452 + conj(IT_0363)*IT_0517 +
       conj(IT_0368)*IT_0527 + conj(IT_0271)*IT_0538 + conj(IT_0406)*IT_0556 +
       conj(IT_0198)*IT_0567 + conj(IT_0306)*IT_0575 + conj(IT_0064)*IT_0584 +
       conj(IT_0323)*IT_0587;
    const ccomplex_t IT_0590 = IT_0019*IT_0050;
    const ccomplex_t IT_0591 = -IT_0590;
    const ccomplex_t IT_0592 = s_14*s_24*IT_0065;
    const ccomplex_t IT_0593 = (s_12 + (-0.5)*IT_0194)*(s_12 + -IT_0592);
    const ccomplex_t IT_0594 = (-2)*IT_0593;
    const ccomplex_t IT_0595 = pow(s_24, 2);
    const ccomplex_t IT_0596 = IT_0065*IT_0595;
    const ccomplex_t IT_0597 = -(s_12 + (-0.5)*IT_0194)*(IT_0016 + -IT_0596);
    const ccomplex_t IT_0598 = (-2)*IT_0597;
    const ccomplex_t IT_0599 = (-2)*IT_0598;
    const ccomplex_t IT_0600 = 0.25*IT_0599;
    const ccomplex_t IT_0601 = -(IT_0015 + -IT_0069)*(s_12 + (-0.5)*IT_0194);
    const ccomplex_t IT_0602 = (-2)*IT_0601;
    const ccomplex_t IT_0603 = (-2)*IT_0602;
    const ccomplex_t IT_0604 = 0.25*IT_0603;
    const ccomplex_t IT_0605 = IT_0594 + IT_0600 + IT_0604;
    const ccomplex_t IT_0606 = s_13*IT_0595;
    const ccomplex_t IT_0607 = IT_0529 + IT_0606;
    const ccomplex_t IT_0608 = s_12*s_24*s_34;
    const ccomplex_t IT_0609 = -IT_0530 + -IT_0608;
    const ccomplex_t IT_0610 = IT_0607 + IT_0609;
    const ccomplex_t IT_0611 = IT_0065*IT_0610;
    const ccomplex_t IT_0612 = (-0.25)*IT_0611;
    const ccomplex_t IT_0613 = s_12*s_24*s_34*IT_0065;
    const ccomplex_t IT_0614 = (-0.5)*IT_0613;
    const ccomplex_t IT_0615 = s_13*IT_0065*IT_0595;
    const ccomplex_t IT_0616 = 0.5*IT_0615;
    const ccomplex_t IT_0617 = s_14*s_23*s_24*IT_0065;
    const ccomplex_t IT_0618 = (-0.5)*IT_0617;
    const ccomplex_t IT_0619 = (-0.5)*IT_0171;
    const ccomplex_t IT_0620 = IT_0172 + IT_0612 + IT_0614 + IT_0616 + IT_0618
       + IT_0619;
    const ccomplex_t IT_0621 = s_13*s_14*s_24*IT_0065;
    const ccomplex_t IT_0622 = s_23*IT_0065*IT_0068;
    const ccomplex_t IT_0623 = -IT_0622;
    const ccomplex_t IT_0624 = s_13*s_14*s_24;
    const ccomplex_t IT_0625 = IT_0558 + IT_0624;
    const ccomplex_t IT_0626 = -IT_0559 + -IT_0560;
    const ccomplex_t IT_0627 = IT_0625 + IT_0626;
    const ccomplex_t IT_0628 = IT_0065*IT_0627;
    const ccomplex_t IT_0629 = (-0.5)*IT_0628;
    const ccomplex_t IT_0630 = s_12*s_14*s_34*IT_0065;
    const ccomplex_t IT_0631 = -IT_0630;
    const ccomplex_t IT_0632 = IT_0081 + IT_0621 + IT_0623 + IT_0629 + IT_0631;
    const ccomplex_t IT_0633 = (-4)*IT_0632;
    const ccomplex_t IT_0634 = (-0.125)*IT_0633;
    const ccomplex_t IT_0635 = IT_0620 + IT_0634;
    const ccomplex_t IT_0636 = IT_0050*IT_0062;
    const ccomplex_t IT_0637 = -IT_0015;
    const ccomplex_t IT_0638 = IT_0069 + IT_0637;
    const ccomplex_t IT_0639 = IT_0307*IT_0638;
    const ccomplex_t IT_0640 = 6*IT_0639;
    const ccomplex_t IT_0641 = IT_0428 + IT_0596;
    const ccomplex_t IT_0642 = IT_0307*IT_0641;
    const ccomplex_t IT_0643 = 6*IT_0642;
    const ccomplex_t IT_0644 = -IT_0592;
    const ccomplex_t IT_0645 = s_12 + IT_0644;
    const ccomplex_t IT_0646 = IT_0307*IT_0645;
    const ccomplex_t IT_0647 = (-12)*IT_0646;
    const ccomplex_t IT_0648 = IT_0640 + IT_0643 + IT_0647;
    const ccomplex_t IT_0649 = conj(IT_0064)*IT_0076;
    const ccomplex_t IT_0650 = IT_0427 + IT_0596;
    const ccomplex_t IT_0651 = s_23*s_24*s_34*IT_0071;
    const ccomplex_t IT_0652 = -IT_0016 + -IT_0651;
    const ccomplex_t IT_0653 = IT_0650 + IT_0652;
    const ccomplex_t IT_0654 = m_C_1*IT_0653;
    const ccomplex_t IT_0655 = (-2)*IT_0654;
    const ccomplex_t IT_0656 = conj(IT_0363)*IT_0655;
    const ccomplex_t IT_0657 = s_14*s_23*s_34*IT_0071;
    const ccomplex_t IT_0658 = -IT_0475 + -IT_0592;
    const ccomplex_t IT_0659 = s_12 + IT_0657 + IT_0658;
    const ccomplex_t IT_0660 = m_C_1*IT_0659;
    const ccomplex_t IT_0661 = 2*IT_0660;
    const ccomplex_t IT_0662 = 0.25*conj(IT_0363);
    const ccomplex_t IT_0663 = s_13*s_24*s_34*IT_0071;
    const ccomplex_t IT_0664 = s_12 + IT_0658 + IT_0663;
    const ccomplex_t IT_0665 = m_C_2*IT_0664;
    const ccomplex_t IT_0666 = (-2)*IT_0665;
    const ccomplex_t IT_0667 = 0.25*conj(IT_0064);
    const ccomplex_t IT_0668 = IT_0074*IT_0307;
    const ccomplex_t IT_0669 = 2*IT_0668;
    const ccomplex_t IT_0670 = IT_0307*IT_0664;
    const ccomplex_t IT_0671 = (-2)*IT_0670;
    const ccomplex_t IT_0672 = IT_0669 + IT_0671;
    const ccomplex_t IT_0673 = 0.25*conj(IT_0198);
    const ccomplex_t IT_0674 = IT_0307*IT_0653;
    const ccomplex_t IT_0675 = 2*IT_0674;
    const ccomplex_t IT_0676 = IT_0307*IT_0659;
    const ccomplex_t IT_0677 = (-2)*IT_0676;
    const ccomplex_t IT_0678 = IT_0675 + IT_0677;
    const ccomplex_t IT_0679 = 0.25*conj(IT_0271);
    const ccomplex_t IT_0680 = -IT_0295;
    const ccomplex_t IT_0681 = s_13 + IT_0680;
    const ccomplex_t IT_0682 = IT_0307*IT_0681;
    const ccomplex_t IT_0683 = (-6)*IT_0682;
    const ccomplex_t IT_0684 = -IT_0391;
    const ccomplex_t IT_0685 = s_23 + IT_0684;
    const ccomplex_t IT_0686 = IT_0307*IT_0685;
    const ccomplex_t IT_0687 = (-6)*IT_0686;
    const ccomplex_t IT_0688 = IT_0683 + IT_0687;
    const ccomplex_t IT_0689 = 0.25*conj(IT_0318);
    const ccomplex_t IT_0690 = s_13*s_23*IT_0068*IT_0071;
    const ccomplex_t IT_0691 = s_14*s_24*IT_0066*IT_0071;
    const ccomplex_t IT_0692 = s_12*IT_0015;
    const ccomplex_t IT_0693 = s_14*s_24*IT_0015*IT_0065;
    const ccomplex_t IT_0694 = -IT_0693;
    const ccomplex_t IT_0695 = s_13*s_23*IT_0015*IT_0065;
    const ccomplex_t IT_0696 = -IT_0695;
    const ccomplex_t IT_0697 = s_12*s_13*s_14*s_34*IT_0071;
    const ccomplex_t IT_0698 = -IT_0697;
    const ccomplex_t IT_0699 = IT_0690 + IT_0691 + IT_0692 + IT_0694 + IT_0696
       + IT_0698;
    const ccomplex_t IT_0700 = 2*IT_0699;
    const ccomplex_t IT_0701 = s_13*s_14*s_23*s_24*IT_0071;
    const ccomplex_t IT_0702 = IT_0066*IT_0071*IT_0595;
    const ccomplex_t IT_0703 = IT_0015*IT_0016;
    const ccomplex_t IT_0704 = IT_0016*IT_0065*IT_0066;
    const ccomplex_t IT_0705 = -IT_0704;
    const ccomplex_t IT_0706 = IT_0015*IT_0065*IT_0595;
    const ccomplex_t IT_0707 = -IT_0706;
    const ccomplex_t IT_0708 = s_12*s_13*s_24*s_34*IT_0071;
    const ccomplex_t IT_0709 = -IT_0708;
    const ccomplex_t IT_0710 = IT_0701 + IT_0702 + IT_0703 + IT_0705 + IT_0707
       + IT_0709;
    const ccomplex_t IT_0711 = 2*IT_0710;
    const ccomplex_t IT_0712 = IT_0700 + IT_0711;
    const ccomplex_t IT_0713 = 0.25*conj(IT_0021);
    const ccomplex_t IT_0714 = IT_0068*IT_0071*IT_0168;
    const ccomplex_t IT_0715 = IT_0016*IT_0065*IT_0068;
    const ccomplex_t IT_0716 = -IT_0715;
    const ccomplex_t IT_0717 = IT_0015*IT_0065*IT_0168;
    const ccomplex_t IT_0718 = -IT_0717;
    const ccomplex_t IT_0719 = s_12*s_14*s_23*s_34*IT_0071;
    const ccomplex_t IT_0720 = -IT_0719;
    const ccomplex_t IT_0721 = IT_0701 + IT_0703 + IT_0714 + IT_0716 + IT_0718
       + IT_0720;
    const ccomplex_t IT_0722 = 2*IT_0721;
    const ccomplex_t IT_0723 = s_14*s_24*IT_0071*IT_0168;
    const ccomplex_t IT_0724 = s_13*s_23*IT_0071*IT_0595;
    const ccomplex_t IT_0725 = s_12*IT_0016;
    const ccomplex_t IT_0726 = s_14*s_24*IT_0016*IT_0065;
    const ccomplex_t IT_0727 = -IT_0726;
    const ccomplex_t IT_0728 = s_13*s_23*IT_0016*IT_0065;
    const ccomplex_t IT_0729 = -IT_0728;
    const ccomplex_t IT_0730 = s_12*s_23*s_24*s_34*IT_0071;
    const ccomplex_t IT_0731 = -IT_0730;
    const ccomplex_t IT_0732 = IT_0723 + IT_0724 + IT_0725 + IT_0727 + IT_0729
       + IT_0731;
    const ccomplex_t IT_0733 = 2*IT_0732;
    const ccomplex_t IT_0734 = IT_0722 + IT_0733;
    const ccomplex_t IT_0735 = 0.25*conj(IT_0167);
    const ccomplex_t IT_0736 = (s_12 + (-0.5)*IT_0194)*(s_13 + -IT_0295);
    const ccomplex_t IT_0737 = (-2)*IT_0736;
    const ccomplex_t IT_0738 = 2*IT_0737;
    const ccomplex_t IT_0739 = (s_12 + (-0.5)*IT_0194)*(s_23 + -IT_0391);
    const ccomplex_t IT_0740 = (-2)*IT_0739;
    const ccomplex_t IT_0741 = 2*IT_0740;
    const ccomplex_t IT_0742 = IT_0738 + IT_0741;
    const ccomplex_t IT_0743 = 0.25*conj(IT_0078);
    const ccomplex_t IT_0744 = -IT_0621;
    const ccomplex_t IT_0745 = 2*IT_0079;
    const ccomplex_t IT_0746 = IT_0081 + IT_0085 + IT_0623 + IT_0744 + IT_0745;
    const ccomplex_t IT_0747 = 2*IT_0746;
    const ccomplex_t IT_0748 = -IT_0615;
    const ccomplex_t IT_0749 = -IT_0617;
    const ccomplex_t IT_0750 = 2*IT_0170;
    const ccomplex_t IT_0751 = IT_0171 + IT_0175 + IT_0748 + IT_0749 + IT_0750;
    const ccomplex_t IT_0752 = 2*IT_0751;
    const ccomplex_t IT_0753 = IT_0747 + IT_0752;
    const ccomplex_t IT_0754 = 0.25*conj(IT_0368);
    const ccomplex_t IT_0755 = (-2)*IT_0191 + 2*IT_0337 + -IT_0464;
    const ccomplex_t IT_0756 = IT_0190 + IT_0755;
    const ccomplex_t IT_0757 = IT_0065*IT_0756;
    const ccomplex_t IT_0758 = -IT_0592 + (-0.5)*IT_0757;
    const ccomplex_t IT_0759 = s_12 + IT_0758;
    const ccomplex_t IT_0760 = m_C_1*IT_0759;
    const ccomplex_t IT_0761 = (-2)*IT_0760;
    const ccomplex_t IT_0762 = 2*IT_0761;
    const ccomplex_t IT_0763 = s_23*s_24*s_34*IT_0065;
    const ccomplex_t IT_0764 = (-0.5)*IT_0423 + 0.5*IT_0595 + -IT_0763;
    const ccomplex_t IT_0765 = IT_0168 + IT_0764;
    const ccomplex_t IT_0766 = IT_0065*IT_0765;
    const ccomplex_t IT_0767 = -IT_0596;
    const ccomplex_t IT_0768 = IT_0016 + IT_0766 + IT_0767;
    const ccomplex_t IT_0769 = m_C_1*IT_0768;
    const ccomplex_t IT_0770 = (-4)*IT_0769;
    const ccomplex_t IT_0771 = IT_0762 + IT_0770;
    const ccomplex_t IT_0772 = 0.25*conj(IT_0390);
    const ccomplex_t IT_0773 = -IT_0464;
    const ccomplex_t IT_0774 = IT_0190 + IT_0773;
    const ccomplex_t IT_0775 = IT_0065*IT_0774;
    const ccomplex_t IT_0776 = -IT_0592 + (-0.25)*IT_0775;
    const ccomplex_t IT_0777 = s_12 + IT_0776;
    const ccomplex_t IT_0778 = m_C_2*IT_0777;
    const ccomplex_t IT_0779 = (-4)*IT_0778;
    const ccomplex_t IT_0780 = 2*IT_0779;
    const ccomplex_t IT_0781 = -IT_0414;
    const ccomplex_t IT_0782 = IT_0068 + IT_0781;
    const ccomplex_t IT_0783 = IT_0065*IT_0782;
    const ccomplex_t IT_0784 = (-4)*IT_0638;
    const ccomplex_t IT_0785 = IT_0783 + IT_0784;
    const ccomplex_t IT_0786 = m_C_2*IT_0785;
    const ccomplex_t IT_0787 = (-2)*IT_0786;
    const ccomplex_t IT_0788 = IT_0780 + IT_0787;
    const ccomplex_t IT_0789 = 0.25*conj(IT_0294);
    const ccomplex_t IT_0790 = IT_0302*IT_0307;
    const ccomplex_t IT_0791 = (-4)*IT_0790;
    const ccomplex_t IT_0792 = 2*IT_0791;
    const ccomplex_t IT_0793 = -IT_0393;
    const ccomplex_t IT_0794 = IT_0392 + IT_0793;
    const ccomplex_t IT_0795 = IT_0065*IT_0794;
    const ccomplex_t IT_0796 = -IT_0391 + (-0.25)*IT_0795;
    const ccomplex_t IT_0797 = s_23 + IT_0796;
    const ccomplex_t IT_0798 = IT_0307*IT_0797;
    const ccomplex_t IT_0799 = (-4)*IT_0798;
    const ccomplex_t IT_0800 = 2*IT_0799;
    const ccomplex_t IT_0801 = IT_0792 + IT_0800;
    const ccomplex_t IT_0802 = 0.25*conj(IT_0306);
    const ccomplex_t IT_0803 = conj(IT_0591)*IT_0605 + (-2)*conj(IT_0323)
      *IT_0635 + 0.25*conj(IT_0636)*IT_0648 + 0.25*IT_0649 + 0.25*IT_0656 +
       IT_0661*IT_0662 + IT_0666*IT_0667 + IT_0672*IT_0673 + IT_0678*IT_0679 +
       IT_0688*IT_0689 + IT_0712*IT_0713 + IT_0734*IT_0735 + IT_0742*IT_0743 +
       IT_0753*IT_0754 + IT_0771*IT_0772 + IT_0788*IT_0789 + IT_0801*IT_0802;
    const ccomplex_t IT_0804 = conj(IT_0064)*IT_0655;
    const ccomplex_t IT_0805 = IT_0076 + IT_0666;
    const ccomplex_t IT_0806 = 0.25*conj(IT_0406);
    const ccomplex_t IT_0807 = 0.25*conj(IT_0323);
    const ccomplex_t IT_0808 = (-2)*conj(IT_0306)*IT_0635 + IT_0605*conj
      (IT_0636) + 0.25*conj(IT_0591)*IT_0648 + IT_0661*IT_0667 + IT_0673*IT_0712
       + IT_0672*IT_0713 + IT_0679*IT_0734 + IT_0678*IT_0735 + IT_0689*IT_0742 +
       IT_0688*IT_0743 + IT_0772*IT_0788 + IT_0771*IT_0789 + 0.25*IT_0804 +
       IT_0662*IT_0805 + IT_0753*IT_0806 + IT_0801*IT_0807;
    const ccomplex_t IT_0809 = s_12*IT_0410;
    const ccomplex_t IT_0810 = m_C_1*IT_0664;
    const ccomplex_t IT_0811 = 2*IT_0810;
    const ccomplex_t IT_0812 = m_C_2*IT_0659;
    const ccomplex_t IT_0813 = (-2)*IT_0812;
    const ccomplex_t IT_0814 = s_23*m_C_1;
    const ccomplex_t IT_0815 = IT_0410*IT_0814;
    const ccomplex_t IT_0816 = 4*IT_0815;
    const ccomplex_t IT_0817 = s_13*m_C_2;
    const ccomplex_t IT_0818 = IT_0410*IT_0817;
    const ccomplex_t IT_0819 = (-4)*IT_0818;
    const ccomplex_t IT_0820 = 0.25*conj(IT_0390)*IT_0412 + IT_0412*IT_0662 +
       IT_0655*IT_0679 + IT_0366*IT_0689 + IT_0188*IT_0743 + IT_0463*IT_0789 +
       IT_0584*IT_0802 + IT_0517*IT_0807 + conj(IT_0064)*IT_0809 + IT_0673
      *IT_0811 + IT_0735*IT_0813 + IT_0806*IT_0816 + IT_0754*IT_0819;
    const ccomplex_t IT_0821 = 0.25*conj(IT_0294)*IT_0412 + IT_0412*IT_0667 +
       IT_0076*IT_0673 + IT_0188*IT_0689 + IT_0655*IT_0735 + IT_0366*IT_0743 +
       IT_0463*IT_0772 + IT_0517*IT_0802 + IT_0584*IT_0807 + conj(IT_0363)
      *IT_0809 + IT_0713*IT_0811 + IT_0679*IT_0813 + IT_0754*IT_0816 + IT_0806
      *IT_0819;
    const ccomplex_t IT_0822 = IT_0178 + IT_0613 + IT_0748 + IT_0749;
    const ccomplex_t IT_0823 = conj(IT_0078) + conj(IT_0368);
    const ccomplex_t IT_0824 = 0.5*conj(IT_0363);
    const ccomplex_t IT_0825 = (s_12 + (-0.5)*IT_0467)*(s_12 + -IT_0475);
    const ccomplex_t IT_0826 = (-2)*IT_0825;
    const ccomplex_t IT_0827 = 2*IT_0826;
    const ccomplex_t IT_0828 = 0.5*conj(IT_0198);
    const ccomplex_t IT_0829 = -(IT_0016 + -IT_0427)*(s_12 + (-0.5)*IT_0467);
    const ccomplex_t IT_0830 = (-2)*IT_0829;
    const ccomplex_t IT_0831 = (-2)*IT_0830;
    const ccomplex_t IT_0832 = 0.5*conj(IT_0271);
    const ccomplex_t IT_0833 = 0.5*conj(IT_0318);
    const ccomplex_t IT_0834 = -IT_0475;
    const ccomplex_t IT_0835 = s_12 + IT_0834;
    const ccomplex_t IT_0836 = IT_0307*IT_0835;
    const ccomplex_t IT_0837 = (-6)*IT_0836;
    const ccomplex_t IT_0838 = 0.5*conj(IT_0021);
    const ccomplex_t IT_0839 = IT_0307*IT_0429;
    const ccomplex_t IT_0840 = 6*IT_0839;
    const ccomplex_t IT_0841 = 0.5*conj(IT_0167);
    const ccomplex_t IT_0842 = 0.5*conj(IT_0390);
    const ccomplex_t IT_0843 = 0.5*conj(IT_0294);
    const ccomplex_t IT_0844 = 0.5*conj(IT_0306);
    const ccomplex_t IT_0845 = 0.5*IT_0804 + conj(IT_0406)*IT_0822 + 0.5
      *IT_0276*IT_0823 + IT_0813*IT_0824 + IT_0827*IT_0828 + IT_0831*IT_0832 +
       IT_0182*IT_0833 + IT_0837*IT_0838 + IT_0840*IT_0841 + IT_0483*IT_0842 +
       IT_0433*IT_0843 + IT_0538*IT_0844;
    const ccomplex_t IT_0846 = IT_0089 + IT_0623 + IT_0630 + IT_0744;
    const ccomplex_t IT_0847 = 0.5*conj(IT_0064);
    const ccomplex_t IT_0848 = -(IT_0015 + -IT_0067)*(s_12 + (-0.5)*IT_0467);
    const ccomplex_t IT_0849 = (-2)*IT_0848;
    const ccomplex_t IT_0850 = (-2)*IT_0849;
    const ccomplex_t IT_0851 = IT_0067 + IT_0637;
    const ccomplex_t IT_0852 = IT_0307*IT_0851;
    const ccomplex_t IT_0853 = 6*IT_0852;
    const ccomplex_t IT_0854 = 0.5*IT_0201*IT_0823 + IT_0076*IT_0824 + IT_0827
      *IT_0832 + IT_0093*IT_0833 + IT_0837*IT_0841 + IT_0422*IT_0842 + IT_0491
      *IT_0843 + IT_0567*IT_0844 + conj(IT_0406)*IT_0846 + IT_0811*IT_0847 +
       IT_0828*IT_0850 + IT_0838*IT_0853;
    const ccomplex_t IT_0855 = (IT_0190 + (-2)*IT_0191)*(1 + 0.5*IT_0408);
    const ccomplex_t IT_0856 = (-4)*IT_0855;
    const ccomplex_t IT_0857 = 0.5*conj(IT_0406);
    const ccomplex_t IT_0858 = m_C_1*m_C_2*IT_0109;
    const ccomplex_t IT_0859 = IT_0410*IT_0858;
    const ccomplex_t IT_0860 = (-4)*IT_0859;
    const ccomplex_t IT_0861 = 0.5*conj(IT_0368);
    const ccomplex_t IT_0862 = 0.5*conj(IT_0323);
    const ccomplex_t IT_0863 = conj(IT_0271)*IT_0822 + IT_0819*IT_0824 +
       IT_0377*IT_0833 + IT_0201*IT_0838 + IT_0276*IT_0841 + IT_0499*IT_0842 +
       IT_0510*IT_0843 + IT_0556*IT_0844 + conj(IT_0198)*IT_0846 + IT_0816
      *IT_0847 + IT_0856*IT_0857 + IT_0860*IT_0861 + IT_0527*IT_0862;
    const ccomplex_t IT_0864 = 0.5*conj(IT_0078);
    const ccomplex_t IT_0865 = conj(IT_0318) + conj(IT_0406);
    const ccomplex_t IT_0866 = 0.5*IT_0865;
    const ccomplex_t IT_0867 = 0.5*IT_0656 + conj(IT_0368)*IT_0822 + IT_0828
      *IT_0837 + IT_0827*IT_0838 + IT_0832*IT_0840 + IT_0831*IT_0841 + IT_0433
      *IT_0842 + IT_0483*IT_0843 + IT_0813*IT_0847 + IT_0538*IT_0862 + IT_0182
      *IT_0864 + IT_0276*IT_0866;
    const ccomplex_t IT_0868 = 0.5*IT_0202 + 0.5*IT_0277 + conj(IT_0167)
      *IT_0822 + IT_0816*IT_0824 + IT_0510*IT_0842 + IT_0499*IT_0843 + IT_0527
      *IT_0844 + conj(IT_0021)*IT_0846 + IT_0819*IT_0847 + IT_0857*IT_0860 +
       IT_0856*IT_0861 + IT_0556*IT_0862 + IT_0377*IT_0864;
    const ccomplex_t IT_0869 = 0.5*IT_0649 + IT_0811*IT_0824 + IT_0832*IT_0837
       + IT_0827*IT_0841 + IT_0491*IT_0842 + IT_0422*IT_0843 + conj(IT_0368)
      *IT_0846 + IT_0838*IT_0850 + IT_0828*IT_0853 + IT_0567*IT_0862 + IT_0093
      *IT_0864 + IT_0201*IT_0866;
    const ccomplex_t IT_0870 = IT_0655 + IT_0661;
    const ccomplex_t IT_0871 = (-0.125)*IT_0064;
    const ccomplex_t IT_0872 = (-0.125)*IT_0198;
    const ccomplex_t IT_0873 = (-0.125)*IT_0271;
    const ccomplex_t IT_0874 = (-0.125)*IT_0318;
    const ccomplex_t IT_0875 = (-0.125)*IT_0021;
    const ccomplex_t IT_0876 = (-0.125)*IT_0167;
    const ccomplex_t IT_0877 = (-0.125)*IT_0078;
    const ccomplex_t IT_0878 = (-0.125)*IT_0390;
    const ccomplex_t IT_0879 = (-0.125)*IT_0294;
    const ccomplex_t IT_0880 = (-0.125)*IT_0077 + IT_0323*IT_0635 + (-0.125)
      *IT_0368*IT_0753 + (-0.125)*IT_0306*IT_0801 + (-0.125)*IT_0363*IT_0870 +
       IT_0666*IT_0871 + IT_0672*IT_0872 + IT_0678*IT_0873 + IT_0688*IT_0874 +
       IT_0712*IT_0875 + IT_0734*IT_0876 + IT_0742*IT_0877 + IT_0771*IT_0878 +
       IT_0788*IT_0879;
    const ccomplex_t IT_0881 = IT_0306*IT_0635 + (-0.125)*IT_0406*IT_0753 + (
      -0.125)*IT_0323*IT_0801 + (-0.125)*IT_0363*IT_0805 + IT_0870*IT_0871 +
       IT_0712*IT_0872 + IT_0734*IT_0873 + IT_0742*IT_0874 + IT_0672*IT_0875 +
       IT_0678*IT_0876 + IT_0688*IT_0877 + IT_0788*IT_0878 + IT_0771*IT_0879;
    const ccomplex_t IT_0882 = conj(IT_0021)*IT_0077 + IT_0078*IT_0404 +
       IT_0318*IT_0407 + IT_0294*IT_0511 + IT_0390*IT_0512 + IT_0323*IT_0588 +
       IT_0306*IT_0589 + 4*IT_0591*IT_0803 + 4*IT_0636*IT_0808 + 4*IT_0064
      *IT_0820 + 4*IT_0363*IT_0821 + 2*IT_0271*IT_0845 + 2*IT_0198*IT_0854 + 2
      *IT_0406*IT_0863 + 2*IT_0167*IT_0867 + 2*IT_0368*IT_0868 + 2*IT_0021
      *IT_0869 + (-8)*conj(IT_0591)*IT_0880 + (-8)*conj(IT_0636)*IT_0881;
    return create_ccomplex_return(IT_0882);
}

