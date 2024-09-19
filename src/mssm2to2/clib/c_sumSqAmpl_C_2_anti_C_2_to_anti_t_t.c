#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_t_t(
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (IT_0004 + IT_0007)*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = cpow(IT_0005, -1);
    const ccomplex_t IT_0015 = m_t*e_em*IT_0002*IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0002*IT_0006;
    const ccomplex_t IT_0021 = IT_0003*IT_0005;
    const ccomplex_t IT_0022 = IT_0009*(IT_0020 + -IT_0021);
    const ccomplex_t IT_0023 = 1.4142135623731*e_em*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = m_t*e_em*IT_0009*IT_0013;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0028 + 
      -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0019 + -IT_0031;
    const ccomplex_t IT_0033 = pow(m_t, 2);
    const ccomplex_t IT_0034 = s_34 + IT_0033;
    const ccomplex_t IT_0035 = s_12*IT_0034;
    const ccomplex_t IT_0036 = V_u2*U_Wm2;
    const ccomplex_t IT_0037 = IT_0002*IT_0036;
    const ccomplex_t IT_0038 = U_d2*V_Wp2;
    const ccomplex_t IT_0039 = IT_0005*IT_0038;
    const ccomplex_t IT_0040 = IT_0009*(IT_0037 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0016*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0002*IT_0038;
    const ccomplex_t IT_0047 = IT_0005*IT_0036;
    const ccomplex_t IT_0048 = IT_0009*(IT_0046 + -IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0026*IT_0050;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0045 + -IT_0053;
    const ccomplex_t IT_0055 = IT_0000*IT_0034;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0057 = cos(alpha);
    const ccomplex_t IT_0058 = IT_0006*IT_0057;
    const ccomplex_t IT_0059 = sin(alpha);
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = IT_0009*(IT_0058 + IT_0060);
    const ccomplex_t IT_0062 = 1.4142135623731*e_em*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_t*e_em*IT_0009*IT_0013
      *IT_0014*IT_0059;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = IT_0056*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0003*IT_0057;
    const ccomplex_t IT_0072 = IT_0006*IT_0059;
    const ccomplex_t IT_0073 = IT_0009*(IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_t*e_em*IT_0009*IT_0013
      *IT_0014*IT_0057;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = IT_0070*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = -s_34;
    const ccomplex_t IT_0084 = IT_0033 + IT_0083;
    const ccomplex_t IT_0085 = s_12*IT_0084;
    const ccomplex_t IT_0086 = IT_0038*IT_0057;
    const ccomplex_t IT_0087 = IT_0036*IT_0059;
    const ccomplex_t IT_0088 = IT_0009*(IT_0086 + IT_0087);
    const ccomplex_t IT_0089 = 1.4142135623731*e_em*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0066*IT_0091;
    const ccomplex_t IT_0093 = IT_0056*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0036*IT_0057;
    const ccomplex_t IT_0096 = IT_0038*IT_0059;
    const ccomplex_t IT_0097 = IT_0009*(IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = 1.4142135623731*e_em*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0078*IT_0100;
    const ccomplex_t IT_0102 = IT_0070*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = -IT_0094 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0000*IT_0084;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0107 = 0.666666666666667*IT_0106;
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = cpow(s_12 + IT_0000 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = s_12*IT_0033;
    const ccomplex_t IT_0114 = s_14*s_23;
    const ccomplex_t IT_0115 = IT_0000*IT_0033;
    const ccomplex_t IT_0116 = s_34*IT_0000;
    const ccomplex_t IT_0117 = s_13*s_24;
    const ccomplex_t IT_0118 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0009*IT_0013*IT_0014*conj(U_sb_00);
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0009*IT_0013*IT_0014*U_sb_00;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = IT_0118*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0009*IT_0013*IT_0014*conj(U_sb_01);
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0009*IT_0013*IT_0014*U_sb_01;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = IT_0126*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0135 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u2
      *e_em*IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *conj(V_u2)*e_em*IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = IT_0136*IT_0138;
    const ccomplex_t IT_0140 = IT_0134*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u2)*e_em*IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0125 + IT_0133 + IT_0141 + IT_0149;
    const ccomplex_t IT_0151 = cpow(IT_0002, -1);
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0009*IT_0013*IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0009*IT_0013*IT_0151;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = IT_0153*IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0158 = IT_0156*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0009*IT_0013*IT_0151;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0009*IT_0013*IT_0151;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_23 + IT_0000 + IT_0033 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = V_tb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0169 = IT_0009*IT_0168;
    const ccomplex_t IT_0170 = IT_0009*IT_0151;
    const ccomplex_t IT_0171 = m_b*U_d2*V_tb*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0169 + (-0.5)*IT_0173);
    const ccomplex_t IT_0175 = m_b*conj(U_d2)*V_tb*e_em*IT_0013*U_sb_10;
    const ccomplex_t IT_0176 = IT_0170*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0179 = IT_0009*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0177 + (-2)*IT_0179);
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0174*IT_0181;
    const ccomplex_t IT_0183 = IT_0118*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0186 = IT_0009*IT_0185;
    const ccomplex_t IT_0187 = m_b*U_d2*V_tb*e_em*IT_0013*conj(U_sb_11);
    const ccomplex_t IT_0188 = IT_0170*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0186 + (-0.5)*IT_0189);
    const ccomplex_t IT_0191 = m_b*conj(U_d2)*V_tb*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0192 = IT_0170*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0195 = IT_0009*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0193 + (-2)*IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0190*IT_0197;
    const ccomplex_t IT_0199 = IT_0126*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm2)*IT_0009;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*conj(V_td)*e_em*U_Wm2*IT_0009;
    const ccomplex_t IT_0203 = IT_0201*IT_0202;
    const ccomplex_t IT_0204 = IT_0134*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0009;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm2*IT_0009;
    const ccomplex_t IT_0208 = IT_0206*IT_0207;
    const ccomplex_t IT_0209 = IT_0147*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = IT_0159 + IT_0167 + IT_0184 + IT_0200 + IT_0205
       + IT_0210;
    const ccomplex_t IT_0212 = 12*IT_0116;
    const ccomplex_t IT_0213 = 24*IT_0115;
    const ccomplex_t IT_0214 = 24*IT_0114;
    const ccomplex_t IT_0215 = 12*IT_0113;
    const ccomplex_t IT_0216 = IT_0212 + IT_0213 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = cos(theta_W);
    const ccomplex_t IT_0218 = cpow(IT_0217, -1);
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0218;
    const ccomplex_t IT_0220 = (-0.666666666666667)*IT_0219;
    const ccomplex_t IT_0221 = IT_0008*IT_0218;
    const ccomplex_t IT_0222 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0223 = IT_0221*IT_0222;
    const ccomplex_t IT_0224 = IT_0009*IT_0217;
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0227 = IT_0224*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0223 + -IT_0225 + (-2)
      *IT_0227);
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = IT_0220*IT_0229;
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = e_em*IT_0221;
    const ccomplex_t IT_0236 = e_em*IT_0224;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0235 + (-3)*IT_0236);
    const ccomplex_t IT_0238 = (-0.166666666666667)*IT_0237;
    const ccomplex_t IT_0239 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0240 = IT_0221*IT_0239;
    const ccomplex_t IT_0241 = IT_0224*IT_0239;
    const ccomplex_t IT_0242 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0243 = IT_0224*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0240 + -IT_0241 + (-2)
      *IT_0243);
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = IT_0238*IT_0245;
    const ccomplex_t IT_0247 = IT_0231*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = 24*IT_0116;
    const ccomplex_t IT_0250 = 48*IT_0115;
    const ccomplex_t IT_0251 = 48*IT_0114;
    const ccomplex_t IT_0252 = 24*IT_0113;
    const ccomplex_t IT_0253 = IT_0249 + IT_0250 + IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = IT_0229*IT_0238;
    const ccomplex_t IT_0255 = IT_0231*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0220*IT_0245;
    const ccomplex_t IT_0259 = IT_0231*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = 48*IT_0117;
    const ccomplex_t IT_0262 = IT_0249 + IT_0250 + IT_0252 + IT_0261;
    const ccomplex_t IT_0263 = m_t*m_C_2;
    const ccomplex_t IT_0264 = -s_14 + -s_23;
    const ccomplex_t IT_0265 = s_13 + s_24;
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = IT_0263*IT_0266;
    const ccomplex_t IT_0268 = (-8)*IT_0267;
    const ccomplex_t IT_0269 = 3*conj(IT_0104);
    const ccomplex_t IT_0270 = 3*conj(IT_0082);
    const ccomplex_t IT_0271 = 2*s_14 + 2*s_23;
    const ccomplex_t IT_0272 = IT_0265 + IT_0271;
    const ccomplex_t IT_0273 = IT_0263*IT_0272;
    const ccomplex_t IT_0274 = IT_0120*IT_0181;
    const ccomplex_t IT_0275 = IT_0118*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0128*IT_0197;
    const ccomplex_t IT_0278 = IT_0126*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0138*IT_0201;
    const ccomplex_t IT_0281 = IT_0134*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0145*IT_0206;
    const ccomplex_t IT_0284 = IT_0147*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = -IT_0276 + -IT_0279 + -IT_0282 + -IT_0285;
    const ccomplex_t IT_0287 = 12*conj(IT_0286);
    const ccomplex_t IT_0288 = IT_0122*IT_0174;
    const ccomplex_t IT_0289 = IT_0118*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = IT_0130*IT_0190;
    const ccomplex_t IT_0292 = IT_0126*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*IT_0292;
    const ccomplex_t IT_0294 = IT_0143*IT_0207;
    const ccomplex_t IT_0295 = IT_0147*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0136*IT_0202;
    const ccomplex_t IT_0298 = IT_0134*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = -IT_0290 + -IT_0293 + -IT_0296 + -IT_0299;
    const ccomplex_t IT_0301 = 12*conj(IT_0300);
    const ccomplex_t IT_0302 = IT_0112*(96*IT_0113 + 96*IT_0114 + 192*IT_0115 
      + 96*IT_0116 + 96*IT_0117) + (IT_0150 + -conj(IT_0150) + IT_0211 + -conj
      (IT_0211))*IT_0216 + -(conj(IT_0234) + conj(IT_0248))*IT_0253 + -(conj
      (IT_0257) + conj(IT_0260))*IT_0262 + 3*IT_0268*(IT_0082 + IT_0104 + (
      -0.333333333333333)*IT_0269 + (-0.333333333333333)*IT_0270) + -IT_0273*((
      -12)*IT_0286 + IT_0287) + -IT_0273*((-12)*IT_0300 + IT_0301);
    const ccomplex_t IT_0303 = s_12*s_34;
    const ccomplex_t IT_0304 = -IT_0117;
    const ccomplex_t IT_0305 = IT_0113 + IT_0114 + IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = -IT_0303;
    const ccomplex_t IT_0307 = -IT_0114;
    const ccomplex_t IT_0308 = IT_0113 + IT_0117 + IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = 12*IT_0286;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0311 = (-6)*IT_0054;
    const ccomplex_t IT_0312 = 0.166666666666667*IT_0055;
    const ccomplex_t IT_0313 = 12*IT_0211;
    const ccomplex_t IT_0314 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0314;
    const ccomplex_t IT_0316 = 12*IT_0150;
    const ccomplex_t IT_0317 = m_t*s_23*m_C_2;
    const ccomplex_t IT_0318 = 0.166666666666667*IT_0317;
    const ccomplex_t IT_0319 = 12*IT_0300;
    const ccomplex_t IT_0320 = 0.166666666666667*IT_0115;
    const ccomplex_t IT_0321 = (-6)*conj(IT_0054);
    const ccomplex_t IT_0322 = 12*conj(IT_0211);
    const ccomplex_t IT_0323 = 12*conj(IT_0150);
    const ccomplex_t IT_0324 = s_13 + s_14;
    const ccomplex_t IT_0325 = IT_0263*IT_0324;
    const ccomplex_t IT_0326 = s_23 + s_24;
    const ccomplex_t IT_0327 = IT_0263*IT_0326;
    const ccomplex_t IT_0328 = 0.166666666666667*IT_0327;
    const ccomplex_t IT_0329 = 3*IT_0082;
    const ccomplex_t IT_0330 = -s_14;
    const ccomplex_t IT_0331 = s_13 + IT_0330;
    const ccomplex_t IT_0332 = IT_0263*IT_0331;
    const ccomplex_t IT_0333 = (-2)*IT_0332;
    const ccomplex_t IT_0334 = 0.166666666666667*IT_0333;
    const ccomplex_t IT_0335 = 3*IT_0104;
    const ccomplex_t IT_0336 = -s_24;
    const ccomplex_t IT_0337 = s_23 + IT_0336;
    const ccomplex_t IT_0338 = IT_0263*IT_0337;
    const ccomplex_t IT_0339 = 2*IT_0338;
    const ccomplex_t IT_0340 = 0.166666666666667*IT_0339;
    const ccomplex_t IT_0341 = 0.166666666666667*IT_0325;
    const ccomplex_t IT_0342 = 0.166666666666667*IT_0305;
    const ccomplex_t IT_0343 = 48*IT_0234;
    const ccomplex_t IT_0344 = 24*IT_0150 + IT_0343;
    const ccomplex_t IT_0345 = 4*IT_0338;
    const ccomplex_t IT_0346 = 12*IT_0032;
    const ccomplex_t IT_0347 = (-4)*IT_0332;
    const ccomplex_t IT_0348 = IT_0329*IT_0347;
    const ccomplex_t IT_0349 = m_t*s_13*m_C_2;
    const ccomplex_t IT_0350 = m_t*s_24*m_C_2;
    const ccomplex_t IT_0351 = (-12)*IT_0032;
    const ccomplex_t IT_0352 = 12*IT_0054;
    const ccomplex_t IT_0353 = IT_0329*IT_0345;
    const ccomplex_t IT_0354 = 12*conj(IT_0032);
    const ccomplex_t IT_0355 = (-12)*conj(IT_0054);
    const ccomplex_t IT_0356 = IT_0270*IT_0347;
    const ccomplex_t IT_0357 = IT_0269*IT_0345;
    const ccomplex_t IT_0358 = -IT_0112;
    const ccomplex_t IT_0359 = IT_0253*IT_0358;
    const ccomplex_t IT_0360 = IT_0115*(24*conj(IT_0150) + 48*conj(IT_0234)) +
       IT_0114*(24*conj(IT_0211) + 48*conj(IT_0248)) + 24*IT_0116*conj(IT_0257) 
      + 24*IT_0113*conj(IT_0260) + 24*conj(IT_0286)*IT_0314 + 24*conj(IT_0300)
      *IT_0317 + IT_0325*IT_0354 + IT_0327*IT_0355 + IT_0356 + IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = IT_0269*IT_0347;
    const ccomplex_t IT_0362 = IT_0270*IT_0345;
    const ccomplex_t IT_0363 = 24*IT_0114*conj(IT_0150) + 24*IT_0115*conj
      (IT_0211) + 24*IT_0113*conj(IT_0257) + 24*IT_0116*conj(IT_0260) + 24*conj
      (IT_0300)*IT_0314 + 24*conj(IT_0286)*IT_0317 + IT_0327*IT_0354 + IT_0325
      *IT_0355 + IT_0359 + IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 24*conj(IT_0234);
    const ccomplex_t IT_0365 = IT_0323 + IT_0364;
    const ccomplex_t IT_0366 = 24*conj(IT_0248);
    const ccomplex_t IT_0367 = IT_0322 + IT_0366;
    const ccomplex_t IT_0368 = 12*conj(IT_0054);
    const ccomplex_t IT_0369 = (-12)*conj(IT_0032);
    const ccomplex_t IT_0370 = 48*conj(IT_0260);
    const ccomplex_t IT_0371 = IT_0262*IT_0358;
    const ccomplex_t IT_0372 = IT_0260*(48*IT_0115*conj(IT_0257) + IT_0287
      *IT_0349 + IT_0301*IT_0350 + IT_0356 + IT_0357 + IT_0116*IT_0365 + IT_0113
      *IT_0367 + IT_0327*IT_0368 + IT_0325*IT_0369 + IT_0117*IT_0370 + IT_0371);
    const ccomplex_t IT_0373 = IT_0257*(48*IT_0117*conj(IT_0257) + IT_0301
      *IT_0349 + IT_0287*IT_0350 + IT_0361 + IT_0362 + IT_0113*IT_0365 + IT_0116
      *IT_0367 + IT_0325*IT_0368 + IT_0327*IT_0369 + IT_0115*IT_0370 + IT_0371);
    const ccomplex_t IT_0374 = IT_0054*(24*IT_0035*conj(IT_0054) + (-24)*conj
      (IT_0032)*IT_0055) + IT_0032*(24*conj(IT_0032)*IT_0035 + (-24)*conj
      (IT_0054)*IT_0055) + IT_0104*((-24)*IT_0085*conj(IT_0104) + 24*conj
      (IT_0082)*IT_0105) + IT_0082*((-24)*conj(IT_0082)*IT_0085 + 24*conj
      (IT_0104)*IT_0105) + -IT_0112*IT_0302 + 6*conj(IT_0286)*(IT_0104*IT_0105 +
       IT_0032*IT_0305 + -IT_0082*IT_0308 + IT_0309*IT_0310 + IT_0311*IT_0312 +
       IT_0313*IT_0315 + IT_0316*IT_0318 + IT_0319*IT_0320) + 6*IT_0286*(conj
      (IT_0104)*IT_0105 + conj(IT_0032)*IT_0305 + -conj(IT_0082)*IT_0308 +
       IT_0301*IT_0320 + IT_0312*IT_0321 + IT_0315*IT_0322 + IT_0318*IT_0323) +
       6*IT_0211*(IT_0301*IT_0318 + IT_0320*IT_0323 + conj(IT_0032)*IT_0325 +
       IT_0321*IT_0328 + IT_0270*IT_0334 + IT_0269*IT_0340) + 6*conj(IT_0211)*
      (IT_0310*IT_0313 + IT_0318*IT_0319 + IT_0316*IT_0320 + IT_0032*IT_0325 +
       IT_0311*IT_0328 + IT_0329*IT_0334 + IT_0335*IT_0340) + 6*conj(IT_0150)*
      (IT_0310*IT_0316 + IT_0315*IT_0319 + IT_0032*IT_0327 + IT_0334*IT_0335 +
       IT_0329*IT_0340 + IT_0311*IT_0341) + 6*IT_0150*(IT_0301*IT_0315 + conj
      (IT_0032)*IT_0327 + IT_0269*IT_0334 + IT_0270*IT_0340 + IT_0321*IT_0341) +
       6*conj(IT_0300)*(IT_0032*IT_0055 + IT_0082*IT_0105 + -IT_0104*IT_0308 +
       IT_0310*IT_0319 + IT_0311*IT_0342) + 6*IT_0300*(conj(IT_0032)*IT_0055 +
       conj(IT_0082)*IT_0105 + -conj(IT_0104)*IT_0308 + IT_0321*IT_0342) + conj
      (IT_0248)*(24*IT_0114*IT_0211 + 24*IT_0286*IT_0314 + 24*IT_0300*IT_0317 + 
      (-12)*IT_0054*IT_0327 + IT_0115*IT_0344 + IT_0335*IT_0345 + IT_0325
      *IT_0346 + IT_0348) + conj(IT_0260)*(IT_0113*IT_0313 + IT_0116*IT_0316 +
       IT_0335*IT_0345 + IT_0348 + IT_0309*IT_0349 + IT_0319*IT_0350 + IT_0325
      *IT_0351 + IT_0327*IT_0352) + conj(IT_0234)*(24*IT_0115*IT_0211 + 24
      *IT_0300*IT_0314 + 24*IT_0286*IT_0317 + (-12)*IT_0054*IT_0325 + IT_0114
      *IT_0344 + IT_0327*IT_0346 + IT_0335*IT_0347 + IT_0353) + conj(IT_0257)*
      (IT_0116*IT_0313 + IT_0113*IT_0316 + IT_0335*IT_0347 + IT_0319*IT_0349 +
       IT_0309*IT_0350 + IT_0327*IT_0351 + IT_0325*IT_0352 + IT_0353) + IT_0248
      *IT_0360 + IT_0234*IT_0363 + IT_0372 + IT_0373;
    return create_ccomplex_return(IT_0374);
}

