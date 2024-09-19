#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_b*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = V_u1*U_Wm2;
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = U_d2*V_Wp1;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*(IT_0009 + IT_0011);
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_C_1, 2);
    const ccomplex_t IT_0017 = pow(m_C_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_Z, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + IT_0021 
      + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0001*IT_0010;
    const ccomplex_t IT_0024 = IT_0003*IT_0008;
    const ccomplex_t IT_0025 = IT_0005*(IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0022*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0020 + -IT_0032;
    const ccomplex_t IT_0034 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = IT_0005*(IT_0035 + IT_0037);
    const ccomplex_t IT_0039 = 1.4142135623731*e_em*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = IT_0018*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0036;
    const ccomplex_t IT_0045 = IT_0003*IT_0034;
    const ccomplex_t IT_0046 = IT_0005*(IT_0044 + -IT_0045);
    const ccomplex_t IT_0047 = 1.4142135623731*e_em*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0029*IT_0048;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0051;
    const ccomplex_t IT_0053 = pow(m_b, 2);
    const ccomplex_t IT_0054 = s_34 + IT_0053;
    const ccomplex_t IT_0055 = s_12*IT_0054;
    const ccomplex_t IT_0056 = cos(alpha);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0010*IT_0056;
    const ccomplex_t IT_0060 = sin(alpha);
    const ccomplex_t IT_0061 = IT_0008*IT_0060;
    const ccomplex_t IT_0062 = IT_0005*(IT_0059 + IT_0061);
    const ccomplex_t IT_0063 = 1.4142135623731*e_em*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0058*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0060;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0008*IT_0056;
    const ccomplex_t IT_0073 = IT_0010*IT_0060;
    const ccomplex_t IT_0074 = IT_0005*(IT_0072 + -IT_0073);
    const ccomplex_t IT_0075 = 1.4142135623731*e_em*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = IT_0036*IT_0056;
    const ccomplex_t IT_0084 = IT_0034*IT_0060;
    const ccomplex_t IT_0085 = IT_0005*(IT_0083 + IT_0084);
    const ccomplex_t IT_0086 = 1.4142135623731*e_em*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0058*IT_0088;
    const ccomplex_t IT_0090 = IT_0067*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0034*IT_0056;
    const ccomplex_t IT_0093 = IT_0036*IT_0060;
    const ccomplex_t IT_0094 = IT_0005*(IT_0092 + -IT_0093);
    const ccomplex_t IT_0095 = 1.4142135623731*e_em*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0071*IT_0097;
    const ccomplex_t IT_0099 = IT_0079*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = -IT_0091 + -IT_0100;
    const ccomplex_t IT_0102 = -s_34;
    const ccomplex_t IT_0103 = IT_0053 + IT_0102;
    const ccomplex_t IT_0104 = s_12*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0108 = IT_0005*IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0003, -1);
    const ccomplex_t IT_0110 = IT_0005*IT_0109;
    const ccomplex_t IT_0111 = m_t*V_tb*conj(V_u2)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0108 + (-0.5)*IT_0113);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0106*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_cb*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*V_cb*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0126 = IT_0124*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_01);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = V_tb*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0131 = IT_0005*IT_0130;
    const ccomplex_t IT_0132 = m_t*V_tb*conj(V_u2)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0133 = IT_0110*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + (-0.5)*IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0129*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0005*IT_0141*V_ub_mod;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005
      *IT_0144*V_ub_mod;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = -IT_0127 + IT_0140 + -IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*e_em
      *IT_0000*IT_0002*IT_0005*IT_0144*V_ub_mod;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005*IT_0141
      *V_ub_mod;
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = IT_0147*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0161 = IT_0110*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0166 = IT_0005*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = IT_0164 + IT_0169;
    const ccomplex_t IT_0171 = IT_0159*IT_0170;
    const ccomplex_t IT_0172 = IT_0117*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0000*IT_0002*IT_0005*U_st_01;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0177 = IT_0110*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0182 = IT_0005*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0180 + IT_0185;
    const ccomplex_t IT_0187 = IT_0175*IT_0186;
    const ccomplex_t IT_0188 = IT_0138*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = -IT_0157 + IT_0173 + IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_cb
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*V_cb*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0192*IT_0194;
    const ccomplex_t IT_0196 = IT_0125*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = s_12*IT_0053;
    const ccomplex_t IT_0199 = s_14*s_23;
    const ccomplex_t IT_0200 = s_12*s_34;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = s_13*s_24;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = IT_0198 + IT_0199 + IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = m_C_1*m_C_2;
    const ccomplex_t IT_0206 = -IT_0053;
    const ccomplex_t IT_0207 = s_34 + IT_0206;
    const ccomplex_t IT_0208 = IT_0205*IT_0207;
    const ccomplex_t IT_0209 = (-8)*IT_0208;
    const ccomplex_t IT_0210 = 3*IT_0082;
    const ccomplex_t IT_0211 = 3*conj(IT_0082);
    const ccomplex_t IT_0212 = m_b*m_C_1;
    const ccomplex_t IT_0213 = -s_24;
    const ccomplex_t IT_0214 = s_23 + IT_0213;
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = 4*IT_0215;
    const ccomplex_t IT_0217 = cos(theta_W);
    const ccomplex_t IT_0218 = cpow(IT_0217, -1);
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0218;
    const ccomplex_t IT_0220 = 0.333333333333333*IT_0219;
    const ccomplex_t IT_0221 = IT_0004*IT_0218;
    const ccomplex_t IT_0222 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0223 = IT_0221*IT_0222;
    const ccomplex_t IT_0224 = IT_0005*IT_0217;
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0227 = IT_0224*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0223 + -IT_0225 + (-2)
      *IT_0227);
    const ccomplex_t IT_0229 = 0.5*IT_0228;
    const ccomplex_t IT_0230 = IT_0220*IT_0229;
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = e_em*IT_0221;
    const ccomplex_t IT_0235 = e_em*IT_0224;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*(IT_0234 + 3*IT_0235);
    const ccomplex_t IT_0237 = (-0.166666666666667)*IT_0236;
    const ccomplex_t IT_0238 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0239 = IT_0221*IT_0238;
    const ccomplex_t IT_0240 = IT_0224*IT_0238;
    const ccomplex_t IT_0241 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0242 = IT_0224*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*(IT_0239 + -IT_0240 + (-2)
      *IT_0242);
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = IT_0237*IT_0244;
    const ccomplex_t IT_0246 = IT_0231*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = 3*IT_0101;
    const ccomplex_t IT_0250 = 3*conj(IT_0101);
    const ccomplex_t IT_0251 = m_b*m_C_2;
    const ccomplex_t IT_0252 = -s_14;
    const ccomplex_t IT_0253 = s_13 + IT_0252;
    const ccomplex_t IT_0254 = IT_0251*IT_0253;
    const ccomplex_t IT_0255 = (-4)*IT_0254;
    const ccomplex_t IT_0256 = m_b*s_13*m_C_2;
    const ccomplex_t IT_0257 = IT_0122*IT_0192;
    const ccomplex_t IT_0258 = IT_0125*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = IT_0106*IT_0159;
    const ccomplex_t IT_0261 = IT_0117*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = IT_0129*IT_0175;
    const ccomplex_t IT_0264 = IT_0138*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0259 + IT_0262 + IT_0265;
    const ccomplex_t IT_0267 = 12*IT_0266;
    const ccomplex_t IT_0268 = IT_0143*IT_0152;
    const ccomplex_t IT_0269 = IT_0147*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = (-12)*IT_0271;
    const ccomplex_t IT_0273 = IT_0267 + IT_0272;
    const ccomplex_t IT_0274 = 12*conj(IT_0266);
    const ccomplex_t IT_0275 = (-12)*conj(IT_0271);
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = 12*IT_0271;
    const ccomplex_t IT_0278 = (-12)*IT_0266;
    const ccomplex_t IT_0279 = IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = 12*conj(IT_0271);
    const ccomplex_t IT_0281 = (-12)*conj(IT_0266);
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 24*IT_0150;
    const ccomplex_t IT_0284 = (-24)*IT_0120;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = 24*conj(IT_0150);
    const ccomplex_t IT_0287 = (-24)*conj(IT_0120);
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0000*IT_0005*IT_0109*IT_0141*V_ub_mod;
    const ccomplex_t IT_0290 = 0.5*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0000*IT_0005*IT_0109*IT_0144*V_ub_mod;
    const ccomplex_t IT_0292 = (-0.5)*IT_0291;
    const ccomplex_t IT_0293 = IT_0290*IT_0292;
    const ccomplex_t IT_0294 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0295 = IT_0293*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0145*IT_0154;
    const ccomplex_t IT_0298 = IT_0147*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = IT_0115*IT_0170;
    const ccomplex_t IT_0301 = IT_0117*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = IT_0136*IT_0186;
    const ccomplex_t IT_0304 = IT_0138*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = IT_0296 + IT_0299 + IT_0302 + IT_0305;
    const ccomplex_t IT_0307 = 12*IT_0190;
    const ccomplex_t IT_0308 = 12*conj(IT_0190);
    const ccomplex_t IT_0309 = 24*conj(IT_0190);
    const ccomplex_t IT_0310 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0000*IT_0005*IT_0109;
    const ccomplex_t IT_0311 = 0.5*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u2)*e_em*IT_0000*IT_0005*IT_0109;
    const ccomplex_t IT_0313 = (-0.5)*IT_0312;
    const ccomplex_t IT_0314 = IT_0311*IT_0313;
    const ccomplex_t IT_0315 = cpow((-2)*s_13 + IT_0016 + IT_0053 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0316 = IT_0314*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = IT_0123*IT_0194;
    const ccomplex_t IT_0319 = IT_0125*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = -IT_0317 + -IT_0320;
    const ccomplex_t IT_0322 = (-12)*IT_0190;
    const ccomplex_t IT_0323 = (-12)*conj(IT_0190);
    const ccomplex_t IT_0324 = 12*IT_0321;
    const ccomplex_t IT_0325 = (-12)*IT_0306;
    const ccomplex_t IT_0326 = 12*conj(IT_0321);
    const ccomplex_t IT_0327 = (-12)*conj(IT_0306);
    const ccomplex_t IT_0328 = 24*IT_0190*conj(IT_0248) + conj(IT_0150)
      *IT_0273 + IT_0150*IT_0276 + conj(IT_0120)*IT_0279 + IT_0120*IT_0282 +
       conj(IT_0233)*IT_0285 + IT_0233*IT_0288 + conj(IT_0306)*IT_0307 + IT_0306
      *IT_0308 + IT_0248*IT_0309 + conj(IT_0321)*IT_0322 + IT_0321*IT_0323 +
       conj(IT_0197)*((-24)*IT_0248 + IT_0324 + IT_0325) + IT_0197*((-24)*conj
      (IT_0248) + IT_0326 + IT_0327);
    const ccomplex_t IT_0329 = IT_0220*IT_0244;
    const ccomplex_t IT_0330 = IT_0231*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = -IT_0331;
    const ccomplex_t IT_0333 = m_b*s_14*m_C_2;
    const ccomplex_t IT_0334 = (-12)*IT_0197;
    const ccomplex_t IT_0335 = IT_0307 + IT_0334;
    const ccomplex_t IT_0336 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0337 = 24*IT_0233 + IT_0267 + IT_0272;
    const ccomplex_t IT_0338 = m_b*s_23*m_C_1;
    const ccomplex_t IT_0339 = 12*IT_0150;
    const ccomplex_t IT_0340 = (-12)*IT_0120;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 48*IT_0332;
    const ccomplex_t IT_0343 = 12*IT_0306;
    const ccomplex_t IT_0344 = 24*IT_0248;
    const ccomplex_t IT_0345 = (-12)*IT_0321;
    const ccomplex_t IT_0346 = IT_0229*IT_0237;
    const ccomplex_t IT_0347 = IT_0231*IT_0346;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*IT_0347;
    const ccomplex_t IT_0349 = m_b*s_24*m_C_1;
    const ccomplex_t IT_0350 = IT_0324 + IT_0325;
    const ccomplex_t IT_0351 = 12*conj(IT_0306);
    const ccomplex_t IT_0352 = (-12)*conj(IT_0321);
    const ccomplex_t IT_0353 = 24*IT_0190*conj(IT_0233) + conj(IT_0197)*((-24)
      *IT_0233 + IT_0277 + IT_0278) + IT_0197*((-24)*conj(IT_0233) + IT_0280 +
       IT_0281) + conj(IT_0248)*IT_0285 + IT_0248*IT_0288 + conj(IT_0266)
      *IT_0307 + IT_0266*IT_0308 + IT_0233*IT_0309 + conj(IT_0271)*IT_0322 +
       IT_0271*IT_0323 + IT_0120*(IT_0326 + IT_0327) + conj(IT_0150)*(IT_0343 +
       IT_0345) + conj(IT_0120)*IT_0350 + IT_0150*(IT_0351 + IT_0352);
    const ccomplex_t IT_0354 = (-12)*conj(IT_0197);
    const ccomplex_t IT_0355 = IT_0308 + IT_0354;
    const ccomplex_t IT_0356 = 24*conj(IT_0233) + IT_0274 + IT_0275;
    const ccomplex_t IT_0357 = 12*conj(IT_0150);
    const ccomplex_t IT_0358 = (-12)*conj(IT_0120);
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = 24*conj(IT_0248);
    const ccomplex_t IT_0361 = 24*conj(IT_0266);
    const ccomplex_t IT_0362 = (-24)*conj(IT_0271);
    const ccomplex_t IT_0363 = 48*conj(IT_0233);
    const ccomplex_t IT_0364 = IT_0361 + IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = 12*IT_0120;
    const ccomplex_t IT_0366 = (-12)*IT_0150;
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = 24*IT_0306;
    const ccomplex_t IT_0369 = (-24)*IT_0321;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 24*conj(IT_0306);
    const ccomplex_t IT_0372 = (-24)*conj(IT_0321);
    const ccomplex_t IT_0373 = 24*IT_0266;
    const ccomplex_t IT_0374 = (-24)*IT_0271;
    const ccomplex_t IT_0375 = (-12)*conj(IT_0150);
    const ccomplex_t IT_0376 = conj(IT_0266)*IT_0267 + conj(IT_0271)*IT_0279 +
       IT_0271*IT_0281 + conj(IT_0190)*IT_0307 + conj(IT_0197)*IT_0322 + IT_0197
      *(12*conj(IT_0197) + IT_0323) + IT_0321*IT_0327 + conj(IT_0150)*IT_0339 +
       conj(IT_0306)*IT_0343 + conj(IT_0321)*IT_0350 + IT_0233*IT_0364 + conj
      (IT_0120)*IT_0367 + conj(IT_0248)*IT_0370 + IT_0248*(48*conj(IT_0248) +
       IT_0371 + IT_0372) + conj(IT_0233)*(IT_0373 + IT_0374) + IT_0120*IT_0375;
    const ccomplex_t IT_0377 = m_C_1*m_C_2*IT_0053;
    const ccomplex_t IT_0378 = IT_0276*IT_0306 + IT_0273*conj(IT_0306) + conj
      (IT_0150)*IT_0307 + IT_0150*IT_0308 + IT_0282*IT_0321 + IT_0279*conj
      (IT_0321) + conj(IT_0120)*(12*IT_0197 + IT_0322) + IT_0120*IT_0323 + 48
      *conj(IT_0332)*IT_0348 + IT_0342*conj(IT_0348) + IT_0248*IT_0364 + conj
      (IT_0197)*IT_0367 + conj(IT_0233)*IT_0370 + IT_0233*(IT_0371 + IT_0372) +
       conj(IT_0248)*(48*IT_0233 + IT_0373 + IT_0374) + IT_0197*IT_0375;
    const ccomplex_t IT_0379 = s_13 + s_14;
    const ccomplex_t IT_0380 = IT_0251*IT_0379;
    const ccomplex_t IT_0381 = 6*IT_0033;
    const ccomplex_t IT_0382 = 6*conj(IT_0033);
    const ccomplex_t IT_0383 = 12*IT_0033;
    const ccomplex_t IT_0384 = 12*conj(IT_0033);
    const ccomplex_t IT_0385 = 12*IT_0052;
    const ccomplex_t IT_0386 = 12*conj(IT_0052);
    const ccomplex_t IT_0387 = (-12)*IT_0033;
    const ccomplex_t IT_0388 = (-12)*conj(IT_0033);
    const ccomplex_t IT_0389 = (-12)*IT_0052;
    const ccomplex_t IT_0390 = (-12)*conj(IT_0052);
    const ccomplex_t IT_0391 = (-6)*IT_0033;
    const ccomplex_t IT_0392 = (-6)*conj(IT_0033);
    const ccomplex_t IT_0393 = conj(IT_0052)*(IT_0266 + -IT_0271) + IT_0052*
      (conj(IT_0266) + -conj(IT_0271)) + (-0.166666666666667)*conj(IT_0306)
      *IT_0381 + (-0.166666666666667)*IT_0306*IT_0382 + (-0.166666666666667)
      *conj(IT_0248)*IT_0383 + (-0.166666666666667)*IT_0248*IT_0384 + (
      -0.166666666666667)*conj(IT_0348)*IT_0385 + (-0.166666666666667)*IT_0348
      *IT_0386 + (-0.166666666666667)*conj(IT_0332)*IT_0387 + (
      -0.166666666666667)*IT_0332*IT_0388 + (-0.166666666666667)*conj(IT_0233)
      *IT_0389 + (-0.166666666666667)*IT_0233*IT_0390 + (-0.166666666666667)
      *conj(IT_0321)*IT_0391 + (-0.166666666666667)*IT_0321*IT_0392;
    const ccomplex_t IT_0394 = IT_0054*IT_0205;
    const ccomplex_t IT_0395 = 6*IT_0052;
    const ccomplex_t IT_0396 = 6*conj(IT_0052);
    const ccomplex_t IT_0397 = 4*conj(IT_0052)*(IT_0033 + 1./4*IT_0190) + 4
      *IT_0052*(conj(IT_0033) + 1./4*conj(IT_0190)) + (-0.166666666666667)*conj
      (IT_0150)*IT_0381 + (-0.166666666666667)*IT_0150*IT_0382 + (
      -0.166666666666667)*conj(IT_0120)*IT_0391 + (-0.166666666666667)*IT_0120
      *IT_0392 + (-0.166666666666667)*conj(IT_0197)*IT_0395 + (
      -0.166666666666667)*IT_0197*IT_0396;
    const ccomplex_t IT_0398 = s_23 + s_24;
    const ccomplex_t IT_0399 = IT_0212*IT_0398;
    const ccomplex_t IT_0400 = (-6)*IT_0052;
    const ccomplex_t IT_0401 = (-6)*conj(IT_0052);
    const ccomplex_t IT_0402 = conj(IT_0266)*IT_0381 + IT_0266*IT_0382 + conj
      (IT_0233)*IT_0383 + IT_0233*IT_0384 + conj(IT_0332)*IT_0385 + IT_0332
      *IT_0386 + conj(IT_0348)*IT_0387 + IT_0348*IT_0388 + conj(IT_0248)*IT_0389
       + IT_0248*IT_0390 + conj(IT_0271)*IT_0391 + IT_0271*IT_0392 + conj
      (IT_0321)*IT_0395 + IT_0321*IT_0396 + conj(IT_0306)*IT_0400 + IT_0306
      *IT_0401;
    const ccomplex_t IT_0403 = -IT_0199;
    const ccomplex_t IT_0404 = IT_0198 + IT_0200 + IT_0202 + IT_0403;
    const ccomplex_t IT_0405 = (-2)*IT_0215;
    const ccomplex_t IT_0406 = 2*IT_0215;
    const ccomplex_t IT_0407 = 2*IT_0254;
    const ccomplex_t IT_0408 = (-2)*IT_0254;
    const ccomplex_t IT_0409 = (-2)*IT_0208;
    const ccomplex_t IT_0410 = 2*IT_0208;
    const ccomplex_t IT_0411 = (-3)*IT_0082;
    const ccomplex_t IT_0412 = (-3)*conj(IT_0082);
    const ccomplex_t IT_0413 = (-3)*IT_0101;
    const ccomplex_t IT_0414 = (-3)*conj(IT_0101);
    const ccomplex_t IT_0415 = (24*IT_0033*conj(IT_0033) + 24*IT_0052*conj
      (IT_0052))*IT_0055 + ((-24)*IT_0082*conj(IT_0082) + (-24)*IT_0101*conj
      (IT_0101))*IT_0104 + (-6)*(conj(IT_0101)*(IT_0120 + -IT_0150) + IT_0101*
      (conj(IT_0120) + -conj(IT_0150)) + -conj(IT_0082)*(IT_0190 + -IT_0197) + 
      -IT_0082*(conj(IT_0190) + -conj(IT_0197)))*IT_0204 + IT_0209*(conj(IT_0101
      )*IT_0210 + IT_0101*IT_0211) + IT_0216*(IT_0211*IT_0233 + IT_0210*conj
      (IT_0233) + conj(IT_0248)*IT_0249 + IT_0248*IT_0250) + (IT_0211*IT_0248 +
       IT_0210*conj(IT_0248) + conj(IT_0233)*IT_0249 + IT_0233*IT_0250)*IT_0255 
      + IT_0256*IT_0328 + conj(IT_0332)*(IT_0216*IT_0249 + IT_0210*IT_0255 +
       IT_0333*IT_0335 + IT_0336*IT_0337 + IT_0338*IT_0341 + IT_0199*IT_0342 +
       IT_0198*(IT_0343 + IT_0344 + IT_0345)) + (IT_0210*IT_0216 + IT_0249
      *IT_0255 + IT_0198*IT_0337 + IT_0335*IT_0338 + IT_0333*IT_0341 + IT_0336*
      (IT_0343 + IT_0344 + IT_0345) + 48*IT_0199*IT_0348)*conj(IT_0348) +
       IT_0349*IT_0353 + IT_0332*(IT_0216*IT_0250 + IT_0211*IT_0255 + IT_0333
      *IT_0355 + IT_0336*IT_0356 + IT_0338*IT_0359 + IT_0198*(IT_0351 + IT_0352 
      + IT_0360)) + IT_0348*(IT_0211*IT_0216 + IT_0250*IT_0255 + IT_0338*IT_0355
       + IT_0198*IT_0356 + IT_0333*IT_0359 + IT_0336*(IT_0351 + IT_0352 +
       IT_0360)) + IT_0202*IT_0376 + IT_0377*IT_0378 + (-6)*IT_0380*IT_0393 + (
      -6)*IT_0394*IT_0397 + IT_0399*IT_0402 + (conj(IT_0190)*IT_0381 + IT_0190
      *IT_0382 + conj(IT_0197)*IT_0391 + IT_0197*IT_0392 + conj(IT_0120)*IT_0395
       + IT_0120*IT_0396 + conj(IT_0150)*IT_0400 + IT_0150*IT_0401)*IT_0404 + 
      (conj(IT_0266)*IT_0405 + conj(IT_0271)*IT_0406 + conj(IT_0306)*IT_0407 +
       conj(IT_0321)*IT_0408 + conj(IT_0150)*IT_0409 + conj(IT_0120)*IT_0410)
      *IT_0411 + (IT_0266*IT_0405 + IT_0271*IT_0406 + IT_0306*IT_0407 + IT_0321
      *IT_0408 + IT_0150*IT_0409 + IT_0120*IT_0410)*IT_0412 + (conj(IT_0306)
      *IT_0405 + conj(IT_0321)*IT_0406 + conj(IT_0266)*IT_0407 + conj(IT_0271)
      *IT_0408 + conj(IT_0190)*IT_0409 + conj(IT_0197)*IT_0410)*IT_0413 + 
      (IT_0306*IT_0405 + IT_0321*IT_0406 + IT_0266*IT_0407 + IT_0271*IT_0408 +
       IT_0190*IT_0409 + IT_0197*IT_0410)*IT_0414;
    return create_ccomplex_return(IT_0415);
}

