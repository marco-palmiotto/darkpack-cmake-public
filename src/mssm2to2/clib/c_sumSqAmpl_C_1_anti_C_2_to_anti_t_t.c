#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0007);
    const ccomplex_t IT_0009 = (-0.166666666666667)*IT_0008;
    const ccomplex_t IT_0010 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = pow(m_C_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0025 = (-0.666666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = cos(alpha);
    const ccomplex_t IT_0032 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = sin(alpha);
    const ccomplex_t IT_0035 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0005*(IT_0033 + IT_0036);
    const ccomplex_t IT_0038 = 1.4142135623731*e_em*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, -1);
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_t*e_em*IT_0005*IT_0034
      *IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0040*IT_0045;
    const ccomplex_t IT_0047 = IT_0030*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0031*IT_0035;
    const ccomplex_t IT_0051 = IT_0032*IT_0034;
    const ccomplex_t IT_0052 = IT_0005*(IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_t*e_em*IT_0005*IT_0031
      *IT_0041*IT_0043;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0049*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0048 + -IT_0060;
    const ccomplex_t IT_0062 = m_t*m_C_1;
    const ccomplex_t IT_0063 = -s_24;
    const ccomplex_t IT_0064 = s_23 + IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = 4*IT_0065;
    const ccomplex_t IT_0067 = 3*IT_0066;
    const ccomplex_t IT_0068 = conj(IT_0061)*IT_0067;
    const ccomplex_t IT_0069 = U_d2*V_Wp1;
    const ccomplex_t IT_0070 = IT_0031*IT_0069;
    const ccomplex_t IT_0071 = V_u1*U_Wm2;
    const ccomplex_t IT_0072 = IT_0034*IT_0071;
    const ccomplex_t IT_0073 = IT_0005*(IT_0070 + IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0045*IT_0076;
    const ccomplex_t IT_0078 = IT_0030*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0031*IT_0071;
    const ccomplex_t IT_0081 = IT_0034*IT_0069;
    const ccomplex_t IT_0082 = IT_0005*(IT_0080 + -IT_0081);
    const ccomplex_t IT_0083 = 1.4142135623731*e_em*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0057*IT_0085;
    const ccomplex_t IT_0087 = IT_0049*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = -IT_0079 + -IT_0088;
    const ccomplex_t IT_0090 = m_t*m_C_2;
    const ccomplex_t IT_0091 = -s_14;
    const ccomplex_t IT_0092 = s_13 + IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (-4)*IT_0093;
    const ccomplex_t IT_0095 = 3*IT_0094;
    const ccomplex_t IT_0096 = conj(IT_0089)*IT_0095;
    const ccomplex_t IT_0097 = IT_0068 + IT_0096;
    const ccomplex_t IT_0098 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0099 = IT_0003*IT_0098;
    const ccomplex_t IT_0100 = IT_0006*IT_0098;
    const ccomplex_t IT_0101 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0102 = IT_0006*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + -IT_0100 + (-2)
      *IT_0102);
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0025*IT_0104;
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 3*conj(IT_0107);
    const ccomplex_t IT_0109 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0110 = pow(m_t, 2);
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u1
      *e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *conj(V_u2)*e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = IT_0111*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u2)*e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0005*IT_0041*IT_0043*conj(U_sb_00);
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0041*IT_0043*U_sb_00;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0005*IT_0041*IT_0043*conj(U_sb_01);
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0041*IT_0043*U_sb_01;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = IT_0137*IT_0139;
    const ccomplex_t IT_0141 = IT_0135*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0118 + IT_0126 + IT_0134 + IT_0142;
    const ccomplex_t IT_0144 = s_23 + s_24;
    const ccomplex_t IT_0145 = IT_0062*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0147 = cos(beta);
    const ccomplex_t IT_0148 = IT_0071*IT_0147;
    const ccomplex_t IT_0149 = IT_0042*IT_0069;
    const ccomplex_t IT_0150 = IT_0005*(IT_0148 + IT_0149);
    const ccomplex_t IT_0151 = 1.4142135623731*e_em*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = m_t*e_em*IT_0005*IT_0041*IT_0043*IT_0147;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = IT_0146*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = pow(m_Z, 2);
    const ccomplex_t IT_0159 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + IT_0158 
      + -reg_prop, -1);
    const ccomplex_t IT_0160 = IT_0069*IT_0147;
    const ccomplex_t IT_0161 = IT_0042*IT_0071;
    const ccomplex_t IT_0162 = IT_0005*(IT_0160 + -IT_0161);
    const ccomplex_t IT_0163 = 1.4142135623731*e_em*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = m_t*e_em*IT_0005*IT_0041;
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = IT_0164*IT_0166;
    const ccomplex_t IT_0168 = IT_0159*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = -IT_0157 + -IT_0169;
    const ccomplex_t IT_0171 = m_t*s_24*m_C_1;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm2*IT_0005;
    const ccomplex_t IT_0173 = IT_0121*IT_0172;
    const ccomplex_t IT_0174 = IT_0119*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = V_tb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0177 = IT_0005*IT_0176;
    const ccomplex_t IT_0178 = cpow(IT_0147, -1);
    const ccomplex_t IT_0179 = IT_0005*IT_0178;
    const ccomplex_t IT_0180 = m_b*U_d2*V_tb*e_em*IT_0041*conj(U_sb_10);
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0177 + (-0.5)*IT_0182);
    const ccomplex_t IT_0184 = IT_0131*IT_0183;
    const ccomplex_t IT_0185 = IT_0127*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*conj(V_td)*e_em*U_Wm2*IT_0005;
    const ccomplex_t IT_0188 = IT_0113*IT_0187;
    const ccomplex_t IT_0189 = IT_0111*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0192 = IT_0005*IT_0191;
    const ccomplex_t IT_0193 = m_b*U_d2*V_tb*e_em*IT_0041*conj(U_sb_11);
    const ccomplex_t IT_0194 = IT_0179*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0192 + (-0.5)*IT_0195);
    const ccomplex_t IT_0197 = IT_0139*IT_0196;
    const ccomplex_t IT_0198 = IT_0135*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = -IT_0175 + -IT_0186 + -IT_0190 + -IT_0199;
    const ccomplex_t IT_0201 = s_12*IT_0110;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0005*IT_0041*IT_0178;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0005*IT_0041*IT_0178;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0203*IT_0205;
    const ccomplex_t IT_0207 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0208 = IT_0206*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0005*IT_0041*IT_0178;
    const ccomplex_t IT_0211 = (-0.5)*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0005*IT_0041*IT_0178;
    const ccomplex_t IT_0213 = (-0.5)*IT_0212;
    const ccomplex_t IT_0214 = IT_0211*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0219 = IT_0005*IT_0218;
    const ccomplex_t IT_0220 = m_b*conj(U_d1)*V_tb*e_em*IT_0041*U_sb_10;
    const ccomplex_t IT_0221 = IT_0179*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0219 + (-0.5)*IT_0222);
    const ccomplex_t IT_0224 = IT_0183*IT_0223;
    const ccomplex_t IT_0225 = IT_0127*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0228 = IT_0187*IT_0227;
    const ccomplex_t IT_0229 = IT_0111*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0232 = IT_0172*IT_0231;
    const ccomplex_t IT_0233 = IT_0119*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0236 = IT_0005*IT_0235;
    const ccomplex_t IT_0237 = m_b*conj(U_d1)*V_tb*e_em*IT_0041*U_sb_11;
    const ccomplex_t IT_0238 = IT_0179*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + (-0.5)*IT_0239);
    const ccomplex_t IT_0241 = IT_0196*IT_0240;
    const ccomplex_t IT_0242 = IT_0135*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = IT_0209 + IT_0217 + IT_0226 + IT_0230 + IT_0234
       + IT_0243;
    const ccomplex_t IT_0245 = m_t*s_13*m_C_2;
    const ccomplex_t IT_0246 = IT_0115*IT_0227;
    const ccomplex_t IT_0247 = IT_0111*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0137*IT_0240;
    const ccomplex_t IT_0250 = IT_0135*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0123*IT_0231;
    const ccomplex_t IT_0253 = IT_0119*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0129*IT_0223;
    const ccomplex_t IT_0256 = IT_0127*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = -IT_0248 + -IT_0251 + -IT_0254 + -IT_0257;
    const ccomplex_t IT_0259 = 12*conj(IT_0107);
    const ccomplex_t IT_0260 = m_C_1*m_C_2*IT_0110;
    const ccomplex_t IT_0261 = m_t*s_23*m_C_1;
    const ccomplex_t IT_0262 = s_14*s_23;
    const ccomplex_t IT_0263 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0264 = IT_0009*IT_0104;
    const ccomplex_t IT_0265 = IT_0020*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = 24*conj(IT_0266);
    const ccomplex_t IT_0268 = IT_0035*IT_0147;
    const ccomplex_t IT_0269 = IT_0032*IT_0042;
    const ccomplex_t IT_0270 = IT_0005*(IT_0268 + IT_0269);
    const ccomplex_t IT_0271 = 1.4142135623731*e_em*IT_0270;
    const ccomplex_t IT_0272 = (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = IT_0154*IT_0272;
    const ccomplex_t IT_0274 = IT_0146*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = IT_0032*IT_0147;
    const ccomplex_t IT_0277 = IT_0035*IT_0042;
    const ccomplex_t IT_0278 = IT_0005*(IT_0276 + -IT_0277);
    const ccomplex_t IT_0279 = 1.4142135623731*e_em*IT_0278;
    const ccomplex_t IT_0280 = 0.5*IT_0279;
    const ccomplex_t IT_0281 = IT_0166*IT_0280;
    const ccomplex_t IT_0282 = IT_0159*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = -IT_0275 + -IT_0283;
    const ccomplex_t IT_0285 = s_13 + s_14;
    const ccomplex_t IT_0286 = IT_0090*IT_0285;
    const ccomplex_t IT_0287 = m_C_1*m_C_2;
    const ccomplex_t IT_0288 = s_34 + IT_0110;
    const ccomplex_t IT_0289 = IT_0287*IT_0288;
    const ccomplex_t IT_0290 = s_12*s_34;
    const ccomplex_t IT_0291 = s_13*s_24;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = IT_0201 + IT_0262 + IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = 6*IT_0286;
    const ccomplex_t IT_0295 = 12*IT_0145;
    const ccomplex_t IT_0296 = s_12*IT_0288;
    const ccomplex_t IT_0297 = 24*IT_0296;
    const ccomplex_t IT_0298 = (-12)*IT_0145;
    const ccomplex_t IT_0299 = 6*IT_0143*IT_0145 + ((-12)*IT_0107 + 12*IT_0266
      )*IT_0286 + ((-24)*IT_0170 + 6*IT_0200)*IT_0289 + 6*IT_0258*IT_0293 +
       IT_0244*IT_0294 + IT_0029*IT_0295 + IT_0284*IT_0297 + IT_0023*IT_0298;
    const ccomplex_t IT_0300 = 2*IT_0065;
    const ccomplex_t IT_0301 = 3*IT_0300;
    const ccomplex_t IT_0302 = (-2)*IT_0093;
    const ccomplex_t IT_0303 = 3*IT_0302;
    const ccomplex_t IT_0304 = 12*IT_0262;
    const ccomplex_t IT_0305 = 12*IT_0261;
    const ccomplex_t IT_0306 = 12*IT_0263;
    const ccomplex_t IT_0307 = 12*IT_0260;
    const ccomplex_t IT_0308 = 12*IT_0109;
    const ccomplex_t IT_0309 = 12*IT_0107;
    const ccomplex_t IT_0310 = 24*IT_0266;
    const ccomplex_t IT_0311 = 24*IT_0260;
    const ccomplex_t IT_0312 = (-6)*IT_0145;
    const ccomplex_t IT_0313 = IT_0284*IT_0294 + IT_0089*IT_0301 + IT_0061
      *IT_0303 + IT_0244*IT_0304 + IT_0200*IT_0305 + IT_0258*IT_0306 + IT_0143
      *IT_0307 + IT_0023*IT_0308 + IT_0201*IT_0309 + IT_0262*IT_0310 + IT_0029
      *IT_0311 + IT_0170*IT_0312;
    const ccomplex_t IT_0314 = 12*IT_0286;
    const ccomplex_t IT_0315 = (-12)*IT_0286;
    const ccomplex_t IT_0316 = (-6)*IT_0289;
    const ccomplex_t IT_0317 = (-6)*IT_0286;
    const ccomplex_t IT_0318 = (-24)*IT_0284*IT_0289 + (-6)*IT_0200*IT_0293 +
       IT_0170*IT_0297 + IT_0145*((-12)*IT_0266 + IT_0309) + IT_0244*IT_0312 +
       IT_0023*IT_0314 + IT_0029*IT_0315 + IT_0258*IT_0316 + IT_0143*IT_0317;
    const ccomplex_t IT_0319 = 48*IT_0107;
    const ccomplex_t IT_0320 = 24*IT_0109;
    const ccomplex_t IT_0321 = 48*IT_0260;
    const ccomplex_t IT_0322 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0323 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0324 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0325 = 0.166666666666667*IT_0170;
    const ccomplex_t IT_0326 = 0.166666666666667*IT_0200;
    const ccomplex_t IT_0327 = 0.166666666666667*IT_0244;
    const ccomplex_t IT_0328 = 0.166666666666667*IT_0258;
    const ccomplex_t IT_0329 = 12*IT_0201;
    const ccomplex_t IT_0330 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0331 = 24*IT_0262;
    const ccomplex_t IT_0332 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0333 = IT_0145*IT_0284 + 0.166666666666667*IT_0109
      *IT_0309 + 0.166666666666667*IT_0260*IT_0310 + IT_0304*IT_0322 + IT_0303
      *IT_0323 + IT_0301*IT_0324 + IT_0317*IT_0325 + IT_0306*IT_0326 + IT_0307
      *IT_0327 + IT_0305*IT_0328 + IT_0329*IT_0330 + IT_0331*IT_0332;
    const ccomplex_t IT_0334 = -IT_0110;
    const ccomplex_t IT_0335 = s_34 + IT_0334;
    const ccomplex_t IT_0336 = IT_0287*IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = 3*IT_0337;
    const ccomplex_t IT_0339 = -IT_0290;
    const ccomplex_t IT_0340 = -IT_0262;
    const ccomplex_t IT_0341 = IT_0201 + IT_0291 + IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = (-6)*IT_0341;
    const ccomplex_t IT_0343 = 12*IT_0171;
    const ccomplex_t IT_0344 = 24*IT_0261;
    const ccomplex_t IT_0345 = IT_0284*IT_0293 + 0.166666666666667*IT_0245
      *IT_0309 + 0.166666666666667*IT_0263*IT_0310 + IT_0305*IT_0322 + IT_0316
      *IT_0325 + IT_0307*IT_0326 + IT_0306*IT_0327 + IT_0304*IT_0328 + IT_0323
      *IT_0338 + IT_0324*IT_0342 + IT_0330*IT_0343 + IT_0332*IT_0344;
    const ccomplex_t IT_0346 = 12*IT_0245;
    const ccomplex_t IT_0347 = 24*IT_0263;
    const ccomplex_t IT_0348 = IT_0284*IT_0289 + -IT_0170*IT_0293 +
       0.166666666666667*IT_0171*IT_0309 + 0.166666666666667*IT_0261*IT_0310 +
       IT_0306*IT_0322 + IT_0304*IT_0326 + IT_0305*IT_0327 + IT_0307*IT_0328 +
       IT_0324*IT_0338 + IT_0323*IT_0342 + IT_0330*IT_0346 + IT_0332*IT_0347;
    const ccomplex_t IT_0349 = -s_34;
    const ccomplex_t IT_0350 = IT_0110 + IT_0349;
    const ccomplex_t IT_0351 = s_12*IT_0350;
    const ccomplex_t IT_0352 = (-24)*IT_0351;
    const ccomplex_t IT_0353 = 0.333333333333333*IT_0089;
    const ccomplex_t IT_0354 = (-8)*IT_0336;
    const ccomplex_t IT_0355 = 3*IT_0354;
    const ccomplex_t IT_0356 = 0.333333333333333*IT_0061;
    const ccomplex_t IT_0357 = 0.333333333333333*IT_0200;
    const ccomplex_t IT_0358 = 0.333333333333333*IT_0244;
    const ccomplex_t IT_0359 = 0.333333333333333*IT_0258;
    const ccomplex_t IT_0360 = 3*IT_0107;
    const ccomplex_t IT_0361 = 0.333333333333333*IT_0360;
    const ccomplex_t IT_0362 = IT_0143*IT_0303 + 3*IT_0352*IT_0353 + 3*IT_0355
      *IT_0356 + 3*IT_0342*IT_0357 + 3*IT_0301*IT_0358 + 3*IT_0338*IT_0359 + 3
      *IT_0066*(IT_0266 + IT_0361);
    const ccomplex_t IT_0363 = IT_0143*IT_0301 + 3*IT_0353*IT_0355 + 3*IT_0352
      *IT_0356 + 3*IT_0338*IT_0357 + 3*IT_0303*IT_0358 + 3*IT_0342*IT_0359 + 3
      *IT_0094*(IT_0266 + IT_0361);
    const ccomplex_t IT_0364 = 24*IT_0201;
    const ccomplex_t IT_0365 = IT_0089*IT_0095;
    const ccomplex_t IT_0366 = IT_0061*IT_0067;
    const ccomplex_t IT_0367 = 48*IT_0023*IT_0291 + IT_0284*IT_0298 + IT_0244
      *IT_0308 + IT_0109*IT_0310 + IT_0170*IT_0314 + IT_0260*IT_0319 + IT_0143
      *IT_0329 + IT_0258*IT_0343 + IT_0200*IT_0346 + IT_0029*IT_0364 + IT_0365 +
       IT_0366;
    const ccomplex_t IT_0368 = 48*IT_0029*IT_0262 + 48*IT_0260*IT_0266 +
       IT_0284*IT_0295 + IT_0244*IT_0311 + IT_0170*IT_0315 + IT_0107*IT_0320 +
       IT_0143*IT_0331 + IT_0258*IT_0344 + IT_0200*IT_0347 + IT_0023*IT_0364 +
       IT_0365 + IT_0366;
    const ccomplex_t IT_0369 = (IT_0023 + IT_0029)*IT_0097 + (IT_0066*IT_0089 
      + IT_0061*IT_0094)*IT_0108 + (IT_0109*IT_0143 + IT_0145*IT_0170 + IT_0171
      *IT_0200 + IT_0201*IT_0244 + IT_0245*IT_0258)*IT_0259 + (IT_0023*IT_0109 +
       IT_0107*IT_0201 + IT_0143*IT_0260 + IT_0200*IT_0261 + IT_0244*IT_0262 +
       IT_0258*IT_0263)*IT_0267 + conj(IT_0284)*IT_0299 + conj(IT_0244)*IT_0313 
      + conj(IT_0170)*IT_0318 + conj(IT_0107)*((-12)*IT_0284*IT_0286 + IT_0201
      *IT_0310 + IT_0291*IT_0319 + IT_0029*IT_0320 + IT_0023*IT_0321) + 3*conj
      (IT_0266)*(IT_0066*IT_0089 + IT_0061*IT_0094 + (-4)*IT_0145*IT_0170 + 16
      *IT_0262*IT_0266 + 4*IT_0284*IT_0286 + 1./3*IT_0029*IT_0321) + 6*conj
      (IT_0143)*IT_0333 + 6*conj(IT_0258)*IT_0345 + 6*conj(IT_0200)*IT_0348 +
       conj(IT_0089)*IT_0362 + conj(IT_0061)*IT_0363 + conj(IT_0023)*IT_0367 +
       conj(IT_0029)*IT_0368;
    return create_ccomplex_return(IT_0369);
}

