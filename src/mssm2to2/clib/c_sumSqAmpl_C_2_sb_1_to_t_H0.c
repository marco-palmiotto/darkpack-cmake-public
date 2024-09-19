#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_1_to_t_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_1_to_t_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = U_d2*V_Wp2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = V_u2*U_Wm2;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*conj(U_d2)*V_tb*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0023 = IT_0010*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + (-2)*IT_0023);
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0030 = U_d1*V_Wp2;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = V_u2*U_Wm1;
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = IT_0010*(IT_0031 + IT_0033);
    const ccomplex_t IT_0035 = 1.4142135623731*e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0039 = IT_0010*IT_0038;
    const ccomplex_t IT_0040 = m_b*conj(U_d1)*V_tb*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + (-0.5)*IT_0042);
    const ccomplex_t IT_0044 = IT_0037*IT_0043;
    const ccomplex_t IT_0045 = IT_0029*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = pow(m_C_2, 2);
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + -IT_0001 + -IT_0047 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0049 = sin(beta);
    const ccomplex_t IT_0050 = cpow(IT_0049, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0010*IT_0018*IT_0050*U_sb_00;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0010
      *IT_0018*IT_0050;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0028 + -IT_0046 + -IT_0057;
    const ccomplex_t IT_0059 = pow(m_H0, 2);
    const ccomplex_t IT_0060 = m_t*m_C_2*IT_0059;
    const ccomplex_t IT_0061 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0062 = IT_0003*IT_0061;
    const ccomplex_t IT_0063 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0064 = IT_0006*IT_0063;
    const ccomplex_t IT_0065 = IT_0010*(IT_0062 + IT_0064);
    const ccomplex_t IT_0066 = 1.4142135623731*e_em*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0052*IT_0068;
    const ccomplex_t IT_0070 = IT_0002*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0025*IT_0054;
    const ccomplex_t IT_0073 = IT_0048*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0010*IT_0018*IT_0050*U_sb_00;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0080 = IT_0006*IT_0079;
    const ccomplex_t IT_0081 = IT_0010*(IT_0078 + IT_0080);
    const ccomplex_t IT_0082 = 1.4142135623731*e_em*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0076*IT_0084;
    const ccomplex_t IT_0086 = IT_0029*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0071 + IT_0074 + IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0010*IT_0018*IT_0050*U_sb_01;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = cpow(V_tb, 2);
    const ccomplex_t IT_0093 = m_b*mu_h*IT_0018;
    const ccomplex_t IT_0094 = -IT_0006*IT_0049;
    const ccomplex_t IT_0095 = IT_0003*IT_0015;
    const ccomplex_t IT_0096 = IT_0094 + IT_0095;
    const ccomplex_t IT_0097 = cos(theta_W);
    const ccomplex_t IT_0098 = cpow(IT_0097, -2);
    const ccomplex_t IT_0099 = pow(m_b, 2);
    const ccomplex_t IT_0100 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0096*(IT_0010 + 1.
      /3*IT_0009*IT_0098) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)*conj
      (U_sb_01) + -1./2*IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0009*IT_0096*IT_0098 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)
      *conj(U_sb_11));
    const ccomplex_t IT_0101 = (0 + _Complex_I*-1)*(e_em*IT_0100 + (
      -1.4142135623731)*A_b*IT_0003*IT_0092*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0091*IT_0102;
    const ccomplex_t IT_0104 = IT_0089*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0107 = IT_0014*IT_0052;
    const ccomplex_t IT_0108 = IT_0002*IT_0106*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0010*IT_0018*IT_0050;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_10;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = IT_0109*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0118 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0096*(IT_0010 + 1.
      /3*IT_0009*IT_0098) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)*conj
      (U_sb_00) + -1./2*IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0009*IT_0096*IT_0098 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)
      *conj(U_sb_10));
    const ccomplex_t IT_0119 = (0 + _Complex_I*-1)*(e_em*IT_0118 + (
      -1.4142135623731)*A_b*IT_0003*IT_0092*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0052*IT_0120;
    const ccomplex_t IT_0122 = IT_0117*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0125 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u2
      *e_em*IT_0010*IT_0018*IT_0050;
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_10;
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = IT_0126*IT_0128;
    const ccomplex_t IT_0130 = IT_0124*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0133 = IT_0037*IT_0076;
    const ccomplex_t IT_0134 = IT_0029*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = m_t*IT_0057;
    const ccomplex_t IT_0136 = m_C_2*IT_0087;
    const ccomplex_t IT_0137 = m_C_2*IT_0071;
    const ccomplex_t IT_0138 = IT_0105 + -IT_0108 + IT_0116 + IT_0123 +
       IT_0131 + -IT_0134 + (-2)*IT_0135 + -IT_0136 + -IT_0137;
    const ccomplex_t IT_0139 = m_t*s_14;
    const ccomplex_t IT_0140 = 6*IT_0139;
    const ccomplex_t IT_0141 = s_34*m_C_2;
    const ccomplex_t IT_0142 = 6*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0144 = IT_0128*IT_0143;
    const ccomplex_t IT_0145 = IT_0124*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0010*IT_0016*IT_0018;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_00;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = IT_0147*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0156 = IT_0113*IT_0155;
    const ccomplex_t IT_0157 = IT_0109*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0010*IT_0016*IT_0018;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_00;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0096*(IT_0010 + 1.
      /3*IT_0009*IT_0098) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)*conj
      (U_sb_00) + -1./2*IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0009*IT_0096*IT_0098 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)
      *conj(U_sb_10));
    const ccomplex_t IT_0168 = (0 + _Complex_I*-1)*(e_em*IT_0167 + (
      -1.4142135623731)*A_b*IT_0003*IT_0092*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0025*IT_0169;
    const ccomplex_t IT_0171 = IT_0117*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0174 = IT_0010*IT_0173;
    const ccomplex_t IT_0175 = m_b*conj(U_d2)*V_tb*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0176 = IT_0017*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + (-0.5)*IT_0177);
    const ccomplex_t IT_0179 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0096*(IT_0010 + 1.
      /3*IT_0009*IT_0098) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)*conj
      (U_sb_01) + -1./2*IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0010*IT_0016*IT_0093*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0009*IT_0096*IT_0098 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0099)
      *conj(U_sb_11));
    const ccomplex_t IT_0180 = (0 + _Complex_I*-1)*(e_em*IT_0179 + (
      -1.4142135623731)*A_b*IT_0003*IT_0092*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0178*IT_0181;
    const ccomplex_t IT_0183 = IT_0089*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = IT_0043*IT_0084;
    const ccomplex_t IT_0186 = IT_0029*IT_0132*IT_0185;
    const ccomplex_t IT_0187 = IT_0025*IT_0068;
    const ccomplex_t IT_0188 = IT_0002*IT_0106*IT_0187;
    const ccomplex_t IT_0189 = m_t*IT_0074;
    const ccomplex_t IT_0190 = m_C_2*IT_0046;
    const ccomplex_t IT_0191 = m_C_2*IT_0028;
    const ccomplex_t IT_0192 = -IT_0146 + -IT_0154 + -IT_0158 + -IT_0166 + 
      -IT_0172 + -IT_0184 + IT_0186 + IT_0188 + 2*IT_0189 + IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = 0.166666666666667*conj(IT_0192);
    const ccomplex_t IT_0194 = s_13*IT_0059;
    const ccomplex_t IT_0195 = s_14*s_34;
    const ccomplex_t IT_0196 = (-2)*IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = (-6)*IT_0197;
    const ccomplex_t IT_0199 = 0.166666666666667*conj(IT_0058);
    const ccomplex_t IT_0200 = m_t*m_C_2;
    const ccomplex_t IT_0201 = 6*IT_0200;
    const ccomplex_t IT_0202 = 0.166666666666667*conj(IT_0088);
    const ccomplex_t IT_0203 = 6*IT_0058*(IT_0060*conj(IT_0088) +
       0.166666666666667*conj(IT_0138)*IT_0140 + IT_0142*IT_0193 + IT_0198
      *IT_0199) + 6*IT_0192*(s_13*conj(IT_0192) + IT_0142*IT_0199 +
       0.166666666666667*conj(IT_0138)*IT_0201 + IT_0140*IT_0202) + 6*IT_0138*
      (s_13*conj(IT_0138) + IT_0140*IT_0199 + IT_0193*IT_0201 + IT_0142*IT_0202)
       + 6*IT_0088*(conj(IT_0058)*IT_0060 + 0.166666666666667*conj(IT_0138)
      *IT_0142 + IT_0140*IT_0193 + IT_0198*IT_0202);
    return create_ccomplex_return(IT_0203);
}

