#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_C_1, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cpow(IT_0006, -1);
    const ccomplex_t IT_0016 = m_c*e_em*IT_0003*IT_0010*IT_0014*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = IT_0002*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_Z, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + IT_0021 
      + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0003*IT_0007;
    const ccomplex_t IT_0024 = IT_0004*IT_0006;
    const ccomplex_t IT_0025 = IT_0010*(IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = m_c*e_em*IT_0010*IT_0014;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0022*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0020 + -IT_0032;
    const ccomplex_t IT_0034 = pow(m_c, 2);
    const ccomplex_t IT_0035 = s_34 + IT_0034;
    const ccomplex_t IT_0036 = s_12*IT_0035;
    const ccomplex_t IT_0037 = V_u1*U_Wm2;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = U_d2*V_Wp1;
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = IT_0010*(IT_0038 + IT_0040);
    const ccomplex_t IT_0042 = 1.4142135623731*e_em*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = IT_0002*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0003*IT_0039;
    const ccomplex_t IT_0048 = IT_0006*IT_0037;
    const ccomplex_t IT_0049 = IT_0010*(IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = 1.4142135623731*e_em*IT_0049;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = IT_0022*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = -IT_0046 + -IT_0054;
    const ccomplex_t IT_0056 = m_C_1*m_C_2;
    const ccomplex_t IT_0057 = IT_0035*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0000 + -IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0059 = cos(alpha);
    const ccomplex_t IT_0060 = IT_0037*IT_0059;
    const ccomplex_t IT_0061 = sin(alpha);
    const ccomplex_t IT_0062 = IT_0039*IT_0061;
    const ccomplex_t IT_0063 = IT_0010*(IT_0060 + -IT_0062);
    const ccomplex_t IT_0064 = 1.4142135623731*e_em*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_c*e_em*IT_0010*IT_0014
      *IT_0015*IT_0059;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0058*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0039*IT_0059;
    const ccomplex_t IT_0074 = IT_0037*IT_0061;
    const ccomplex_t IT_0075 = IT_0010*(IT_0073 + IT_0074);
    const ccomplex_t IT_0076 = 1.4142135623731*e_em*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_c*e_em*IT_0010*IT_0014
      *IT_0015*IT_0061;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0072*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0071 + -IT_0083;
    const ccomplex_t IT_0085 = -s_34;
    const ccomplex_t IT_0086 = IT_0034 + IT_0085;
    const ccomplex_t IT_0087 = s_12*IT_0086;
    const ccomplex_t IT_0088 = IT_0007*IT_0059;
    const ccomplex_t IT_0089 = IT_0004*IT_0061;
    const ccomplex_t IT_0090 = IT_0010*(IT_0088 + IT_0089);
    const ccomplex_t IT_0091 = 1.4142135623731*e_em*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = IT_0080*IT_0093;
    const ccomplex_t IT_0095 = IT_0072*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0004*IT_0059;
    const ccomplex_t IT_0098 = IT_0007*IT_0061;
    const ccomplex_t IT_0099 = IT_0010*(IT_0097 + -IT_0098);
    const ccomplex_t IT_0100 = 1.4142135623731*e_em*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0068*IT_0102;
    const ccomplex_t IT_0104 = IT_0058*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = -IT_0096 + -IT_0105;
    const ccomplex_t IT_0107 = -IT_0034;
    const ccomplex_t IT_0108 = s_34 + IT_0107;
    const ccomplex_t IT_0109 = IT_0056*IT_0108;
    const ccomplex_t IT_0110 = (-8)*IT_0109;
    const ccomplex_t IT_0111 = 3*IT_0084;
    const ccomplex_t IT_0112 = 3*conj(IT_0084);
    const ccomplex_t IT_0113 = cos(theta_W);
    const ccomplex_t IT_0114 = cpow(IT_0113, -1);
    const ccomplex_t IT_0115 = IT_0009*IT_0114;
    const ccomplex_t IT_0116 = e_em*IT_0115;
    const ccomplex_t IT_0117 = IT_0010*IT_0113;
    const ccomplex_t IT_0118 = e_em*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0116 + (-3)*IT_0118);
    const ccomplex_t IT_0120 = (-0.166666666666667)*IT_0119;
    const ccomplex_t IT_0121 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0122 = IT_0115*IT_0121;
    const ccomplex_t IT_0123 = IT_0117*IT_0121;
    const ccomplex_t IT_0124 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0125 = IT_0117*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0122 + -IT_0123 + (-2)
      *IT_0125);
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0120*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0000 + -IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = m_c*m_C_1;
    const ccomplex_t IT_0133 = -s_24;
    const ccomplex_t IT_0134 = s_23 + IT_0133;
    const ccomplex_t IT_0135 = IT_0132*IT_0134;
    const ccomplex_t IT_0136 = 4*IT_0135;
    const ccomplex_t IT_0137 = IT_0112*IT_0136;
    const ccomplex_t IT_0138 = m_c*m_C_2;
    const ccomplex_t IT_0139 = -s_14;
    const ccomplex_t IT_0140 = s_13 + IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = (-4)*IT_0141;
    const ccomplex_t IT_0143 = 3*conj(IT_0106);
    const ccomplex_t IT_0144 = IT_0142*IT_0143;
    const ccomplex_t IT_0145 = IT_0137 + IT_0144;
    const ccomplex_t IT_0146 = s_13 + s_14;
    const ccomplex_t IT_0147 = IT_0138*IT_0146;
    const ccomplex_t IT_0148 = 12*conj(IT_0033);
    const ccomplex_t IT_0149 = s_14*s_23;
    const ccomplex_t IT_0150 = cpow(IT_0003, -1);
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_cd*e_em*IT_0010*IT_0014*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_cd)*e_em*IT_0010*IT_0014*IT_0150;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0010*IT_0014*IT_0150;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0010*IT_0014*IT_0150;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*V_cd*e_em*conj(U_Wm1)*IT_0010;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*conj(V_cd)*e_em*U_Wm2*IT_0010;
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0010;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*conj(V_cs)*e_em*U_Wm2*IT_0010;
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0177 = IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = V_cb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0180 = IT_0010*IT_0179;
    const ccomplex_t IT_0181 = IT_0010*IT_0150;
    const ccomplex_t IT_0182 = m_b*U_d2*V_cb*e_em*IT_0014*conj(U_sb_10);
    const ccomplex_t IT_0183 = IT_0181*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0180 + (-0.5)*IT_0184);
    const ccomplex_t IT_0186 = V_cb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0187 = IT_0010*IT_0186;
    const ccomplex_t IT_0188 = m_b*conj(U_d1)*V_cb*e_em*IT_0014*U_sb_10;
    const ccomplex_t IT_0189 = IT_0181*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0187 + (-0.5)*IT_0190);
    const ccomplex_t IT_0192 = IT_0185*IT_0191;
    const ccomplex_t IT_0193 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0194 = IT_0192*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = V_cb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0197 = IT_0010*IT_0196;
    const ccomplex_t IT_0198 = m_b*U_d2*V_cb*e_em*IT_0014*conj(U_sb_11);
    const ccomplex_t IT_0199 = IT_0181*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + (-0.5)*IT_0200);
    const ccomplex_t IT_0202 = V_cb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0203 = IT_0010*IT_0202;
    const ccomplex_t IT_0204 = m_b*conj(U_d1)*V_cb*e_em*IT_0014*U_sb_11;
    const ccomplex_t IT_0205 = IT_0181*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0203 + (-0.5)*IT_0206);
    const ccomplex_t IT_0208 = IT_0201*IT_0207;
    const ccomplex_t IT_0209 = cpow((-2)*s_23 + IT_0001 + IT_0034 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = IT_0158 + IT_0166 + IT_0172 + IT_0178 + IT_0195
       + IT_0211;
    const ccomplex_t IT_0213 = 24*conj(IT_0212);
    const ccomplex_t IT_0214 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0215 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0216 = IT_0115*IT_0215;
    const ccomplex_t IT_0217 = IT_0117*IT_0215;
    const ccomplex_t IT_0218 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0219 = IT_0117*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*(IT_0216 + -IT_0217 + (-2)
      *IT_0219);
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = IT_0120*IT_0221;
    const ccomplex_t IT_0223 = IT_0129*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = 24*conj(IT_0225);
    const ccomplex_t IT_0227 = m_C_1*m_C_2*IT_0034;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u1
      *e_em*IT_0010*IT_0014*IT_0015;
    const ccomplex_t IT_0229 = 0.5*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0010*IT_0014*IT_0015;
    const ccomplex_t IT_0231 = 0.5*IT_0230;
    const ccomplex_t IT_0232 = IT_0229*IT_0231;
    const ccomplex_t IT_0233 = IT_0170*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u2)*e_em*IT_0010*IT_0014*IT_0015*conj(U_sb_00);
    const ccomplex_t IT_0236 = 0.5*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0010*IT_0014*IT_0015*U_sb_00;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = IT_0236*IT_0238;
    const ccomplex_t IT_0240 = IT_0193*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u2)*e_em*IT_0010*IT_0014*IT_0015*conj(U_sb_01);
    const ccomplex_t IT_0243 = 0.5*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0010*IT_0014*IT_0015*U_sb_01;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = IT_0243*IT_0245;
    const ccomplex_t IT_0247 = IT_0209*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0010*IT_0014*IT_0015;
    const ccomplex_t IT_0250 = 0.5*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0010*IT_0014*IT_0015;
    const ccomplex_t IT_0252 = 0.5*IT_0251;
    const ccomplex_t IT_0253 = IT_0250*IT_0252;
    const ccomplex_t IT_0254 = IT_0176*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0234 + IT_0241 + IT_0248 + IT_0255;
    const ccomplex_t IT_0257 = 24*conj(IT_0256);
    const ccomplex_t IT_0258 = m_c*s_23*m_C_1;
    const ccomplex_t IT_0259 = IT_0168*IT_0229;
    const ccomplex_t IT_0260 = IT_0170*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0174*IT_0250;
    const ccomplex_t IT_0263 = IT_0176*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = IT_0185*IT_0238;
    const ccomplex_t IT_0266 = IT_0193*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = IT_0201*IT_0245;
    const ccomplex_t IT_0269 = IT_0209*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = -IT_0261 + -IT_0264 + -IT_0267 + -IT_0270;
    const ccomplex_t IT_0272 = 24*conj(IT_0271);
    const ccomplex_t IT_0273 = m_c*s_14*m_C_2;
    const ccomplex_t IT_0274 = IT_0173*IT_0252;
    const ccomplex_t IT_0275 = IT_0176*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0191*IT_0236;
    const ccomplex_t IT_0278 = IT_0193*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0167*IT_0231;
    const ccomplex_t IT_0281 = IT_0170*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0207*IT_0243;
    const ccomplex_t IT_0284 = IT_0209*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = -IT_0276 + -IT_0279 + -IT_0282 + -IT_0285;
    const ccomplex_t IT_0287 = 24*conj(IT_0286);
    const ccomplex_t IT_0288 = s_23 + s_24;
    const ccomplex_t IT_0289 = IT_0132*IT_0288;
    const ccomplex_t IT_0290 = (-12)*conj(IT_0055);
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0114;
    const ccomplex_t IT_0292 = (-0.666666666666667)*IT_0291;
    const ccomplex_t IT_0293 = IT_0221*IT_0292;
    const ccomplex_t IT_0294 = IT_0129*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = -IT_0295;
    const ccomplex_t IT_0297 = IT_0112*IT_0142;
    const ccomplex_t IT_0298 = IT_0136*IT_0143;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = s_12*IT_0034;
    const ccomplex_t IT_0301 = 12*conj(IT_0212);
    const ccomplex_t IT_0302 = 12*conj(IT_0256);
    const ccomplex_t IT_0303 = m_c*s_13*m_C_2;
    const ccomplex_t IT_0304 = 12*conj(IT_0271);
    const ccomplex_t IT_0305 = m_c*s_24*m_C_1;
    const ccomplex_t IT_0306 = 12*conj(IT_0286);
    const ccomplex_t IT_0307 = 12*conj(IT_0055);
    const ccomplex_t IT_0308 = (-12)*conj(IT_0033);
    const ccomplex_t IT_0309 = IT_0127*IT_0292;
    const ccomplex_t IT_0310 = IT_0129*IT_0309;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*IT_0310;
    const ccomplex_t IT_0312 = s_13*s_24;
    const ccomplex_t IT_0313 = s_12*s_34;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = -IT_0149;
    const ccomplex_t IT_0316 = IT_0300 + IT_0312 + IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = (-2)*IT_0109;
    const ccomplex_t IT_0318 = -IT_0312;
    const ccomplex_t IT_0319 = IT_0149 + IT_0300 + IT_0313 + IT_0318;
    const ccomplex_t IT_0320 = 6*conj(IT_0033);
    const ccomplex_t IT_0321 = (-6)*conj(IT_0055);
    const ccomplex_t IT_0322 = (-6)*IT_0286*((-0.166666666666667)*IT_0273
      *IT_0301 + (-0.166666666666667)*IT_0258*IT_0302 + (-0.166666666666667)
      *IT_0227*IT_0304 + conj(IT_0106)*IT_0316 + (-0.166666666666667)*IT_0112
      *IT_0317 + (-0.166666666666667)*IT_0319*IT_0320 + (-0.166666666666667)
      *IT_0057*IT_0321);
    const ccomplex_t IT_0323 = 2*IT_0135;
    const ccomplex_t IT_0324 = (-2)*IT_0141;
    const ccomplex_t IT_0325 = 12*IT_0212;
    const ccomplex_t IT_0326 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0327 = (-6)*IT_0055;
    const ccomplex_t IT_0328 = 0.166666666666667*IT_0289;
    const ccomplex_t IT_0329 = 12*IT_0271;
    const ccomplex_t IT_0330 = 0.166666666666667*IT_0258;
    const ccomplex_t IT_0331 = 0.166666666666667*IT_0323;
    const ccomplex_t IT_0332 = 3*IT_0106;
    const ccomplex_t IT_0333 = 0.166666666666667*IT_0324;
    const ccomplex_t IT_0334 = 12*IT_0256;
    const ccomplex_t IT_0335 = 0.166666666666667*IT_0273;
    const ccomplex_t IT_0336 = 0.166666666666667*IT_0227;
    const ccomplex_t IT_0337 = 12*IT_0286;
    const ccomplex_t IT_0338 = 0.166666666666667*IT_0317;
    const ccomplex_t IT_0339 = 0.166666666666667*IT_0319;
    const ccomplex_t IT_0340 = 12*IT_0055;
    const ccomplex_t IT_0341 = (-12)*IT_0033;
    const ccomplex_t IT_0342 = 48*IT_0311;
    const ccomplex_t IT_0343 = IT_0111*IT_0142;
    const ccomplex_t IT_0344 = conj(IT_0225)*(48*IT_0225*IT_0312 + IT_0214
      *IT_0325 + IT_0303*IT_0329 + IT_0136*IT_0332 + IT_0300*IT_0334 + IT_0305
      *IT_0337 + IT_0147*IT_0340 + IT_0289*IT_0341 + IT_0227*IT_0342 + IT_0343);
    const ccomplex_t IT_0345 = 24*IT_0256;
    const ccomplex_t IT_0346 = 48*IT_0296;
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = 12*IT_0033;
    const ccomplex_t IT_0349 = 24*IT_0225;
    const ccomplex_t IT_0350 = 24*IT_0271;
    const ccomplex_t IT_0351 = 24*IT_0286;
    const ccomplex_t IT_0352 = 24*IT_0311;
    const ccomplex_t IT_0353 = 48*IT_0131;
    const ccomplex_t IT_0354 = (-12)*IT_0055;
    const ccomplex_t IT_0355 = IT_0136*IT_0332 + IT_0343 + IT_0149*IT_0347 +
       IT_0289*IT_0348 + IT_0300*IT_0349 + IT_0273*IT_0350 + IT_0258*IT_0351 +
       IT_0214*IT_0352 + IT_0227*(24*IT_0212 + IT_0353) + IT_0147*IT_0354;
    const ccomplex_t IT_0356 = IT_0111*IT_0136;
    const ccomplex_t IT_0357 = 48*IT_0225*IT_0227 + IT_0300*(24*IT_0131 +
       IT_0325) + IT_0305*IT_0329 + IT_0142*IT_0332 + IT_0214*(24*IT_0296 +
       IT_0334) + IT_0303*IT_0337 + IT_0289*IT_0340 + IT_0147*IT_0341 + IT_0312
      *IT_0342 + IT_0356;
    const ccomplex_t IT_0358 = IT_0142*IT_0332 + IT_0227*IT_0347 + IT_0147
      *IT_0348 + IT_0214*IT_0349 + IT_0258*IT_0350 + IT_0273*IT_0351 + IT_0300
      *IT_0352 + IT_0149*(24*IT_0212 + IT_0353) + IT_0289*IT_0354 + IT_0356;
    const ccomplex_t IT_0359 = IT_0055*(24*IT_0036*conj(IT_0055) + (-24)*conj
      (IT_0033)*IT_0057) + IT_0033*(24*conj(IT_0033)*IT_0036 + (-24)*conj
      (IT_0055)*IT_0057) + (-24)*IT_0084*conj(IT_0084)*IT_0087 + conj(IT_0106)
      *IT_0110*IT_0111 + IT_0106*((-24)*IT_0087*conj(IT_0106) + IT_0110*IT_0112)
       + IT_0131*(IT_0145 + IT_0147*IT_0148 + IT_0149*IT_0213 + IT_0214*IT_0226 
      + IT_0227*IT_0257 + IT_0258*IT_0272 + IT_0273*IT_0287 + IT_0289*IT_0290) +
       IT_0296*(IT_0213*IT_0227 + IT_0149*IT_0257 + IT_0272*IT_0273 + IT_0258
      *IT_0287 + IT_0148*IT_0289 + IT_0147*IT_0290 + IT_0299 + IT_0226*IT_0300) 
      + IT_0225*(IT_0299 + IT_0214*IT_0301 + IT_0300*IT_0302 + IT_0303*IT_0304 +
       IT_0305*IT_0306 + IT_0147*IT_0307 + IT_0289*IT_0308) + (IT_0145 + IT_0300
      *IT_0301 + IT_0214*IT_0302 + IT_0304*IT_0305 + IT_0303*IT_0306 + IT_0289
      *IT_0307 + IT_0147*IT_0308)*IT_0311 + IT_0322 + IT_0256*(IT_0227*IT_0301 +
       IT_0273*IT_0304 + IT_0289*IT_0320 + IT_0147*IT_0321 + IT_0143*IT_0323 +
       IT_0112*IT_0324) + 6*IT_0212*(conj(IT_0033)*IT_0147 + IT_0321*IT_0328 +
       IT_0304*IT_0330 + IT_0112*IT_0331 + IT_0143*IT_0333) + 6*conj(IT_0212)*
      (IT_0033*IT_0147 + IT_0325*IT_0326 + IT_0327*IT_0328 + IT_0329*IT_0330 +
       IT_0111*IT_0331 + IT_0332*IT_0333) + 6*conj(IT_0256)*(IT_0033*IT_0289 +
       0.166666666666667*IT_0147*IT_0327 + IT_0331*IT_0332 + IT_0111*IT_0333 +
       IT_0326*IT_0334 + IT_0329*IT_0335 + IT_0325*IT_0336) + (-6)*conj(IT_0286)
      *(IT_0106*IT_0316 + -IT_0033*IT_0319 + (-0.166666666666667)*IT_0057
      *IT_0327 + -IT_0330*IT_0334 + -IT_0325*IT_0335 + -IT_0329*IT_0336 + 
      -IT_0326*IT_0337 + -IT_0111*IT_0338) + 6*IT_0271*(conj(IT_0033)*IT_0057 + 
      -conj(IT_0084)*IT_0316 + IT_0143*IT_0338 + IT_0321*IT_0339) + 6*conj
      (IT_0271)*(IT_0033*IT_0057 + -IT_0084*IT_0316 + IT_0326*IT_0329 + IT_0332
      *IT_0338 + IT_0327*IT_0339) + IT_0344 + conj(IT_0296)*IT_0355 + conj
      (IT_0311)*IT_0357 + conj(IT_0131)*IT_0358;
    return create_ccomplex_return(IT_0359);
}

