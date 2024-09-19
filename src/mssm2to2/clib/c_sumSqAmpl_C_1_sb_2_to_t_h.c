#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_2_to_t_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_2_to_t_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = V_u1*U_Wm1;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = U_d1*V_Wp1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (IT_0003 + -IT_0006)*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = sin(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0008*IT_0013*IT_0015*U_sb_01;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_t, 2);
    const ccomplex_t IT_0020 = pow(m_sb_2, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0000*IT_0018*IT_0021;
    const ccomplex_t IT_0023 = pow(m_C_1, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0004*U_sb_11;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0008*IT_0013*IT_0015*U_sb_00;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = cpow(V_tb, 2);
    const ccomplex_t IT_0036 = m_b*mu_h*IT_0013;
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = IT_0004*IT_0037;
    const ccomplex_t IT_0040 = IT_0001*IT_0014;
    const ccomplex_t IT_0041 = IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -2);
    const ccomplex_t IT_0044 = pow(m_b, 2);
    const ccomplex_t IT_0045 = 2*conj(U_sb_00)*(((-0.5)*m_W*IT_0041*(IT_0008 +
       1./3*IT_0007*IT_0043) + IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*U_sb_01 
      + 1./2*IT_0001*IT_0008*IT_0036*IT_0038*U_sb_11) + conj(U_sb_10)*(IT_0001
      *IT_0008*IT_0036*IT_0038*U_sb_01 + (-0.666666666666667)*(m_W*IT_0007
      *IT_0041*IT_0043 + (-3)*IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*U_sb_11);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(e_em*IT_0045 +
       1.4142135623731*A_b*IT_0004*IT_0035*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0034*IT_0047;
    const ccomplex_t IT_0049 = IT_0032*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0053 = V_u1*U_Wm2;
    const ccomplex_t IT_0054 = IT_0001*IT_0053;
    const ccomplex_t IT_0055 = U_d2*V_Wp1;
    const ccomplex_t IT_0056 = IT_0004*IT_0055;
    const ccomplex_t IT_0057 = IT_0008*(IT_0054 + -IT_0056);
    const ccomplex_t IT_0058 = 1.4142135623731*e_em*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0008*IT_0013*IT_0015*U_sb_01;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0051*IT_0052*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0066 = 2*U_sb_01*(((-0.5)*m_W*IT_0041*(IT_0008 + 1./3
      *IT_0007*IT_0043) + IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*conj(U_sb_01)
       + 1./2*IT_0001*IT_0008*IT_0036*IT_0038*conj(U_sb_11)) + U_sb_11*(IT_0001
      *IT_0008*IT_0036*IT_0038*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0007
      *IT_0041*IT_0043 + (-3)*IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*conj
      (U_sb_11));
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(e_em*IT_0066 +
       1.4142135623731*A_b*IT_0004*IT_0035*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0017*IT_0068;
    const ccomplex_t IT_0070 = IT_0065*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u1
      *e_em*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0004*U_sb_11;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0072*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_12 + -IT_0020 + -IT_0023 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_t*e_em*IT_0001*IT_0008
      *IT_0013*IT_0015;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = m_t*IT_0085;
    const ccomplex_t IT_0087 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0090 = IT_0004*IT_0089;
    const ccomplex_t IT_0091 = IT_0008*(IT_0088 + -IT_0090);
    const ccomplex_t IT_0092 = 1.4142135623731*e_em*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0017*IT_0094;
    const ccomplex_t IT_0096 = IT_0021*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0097;
    const ccomplex_t IT_0099 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0100 = IT_0001*IT_0099;
    const ccomplex_t IT_0101 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0102 = IT_0004*IT_0101;
    const ccomplex_t IT_0103 = IT_0008*(IT_0100 + -IT_0102);
    const ccomplex_t IT_0104 = 1.4142135623731*e_em*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0062*IT_0106;
    const ccomplex_t IT_0108 = IT_0052*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = m_C_1*IT_0109;
    const ccomplex_t IT_0111 = -IT_0022 + IT_0031 + IT_0050 + -IT_0064 +
       IT_0071 + IT_0079 + (-2)*IT_0086 + -IT_0098 + -IT_0110;
    const ccomplex_t IT_0112 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0113 = IT_0008*IT_0112;
    const ccomplex_t IT_0114 = IT_0008*IT_0038;
    const ccomplex_t IT_0115 = m_b*conj(U_d1)*V_tb*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0113 + (-0.5)*IT_0117);
    const ccomplex_t IT_0119 = IT_0012*IT_0118;
    const ccomplex_t IT_0120 = IT_0021*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = m_C_1*IT_0121;
    const ccomplex_t IT_0123 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0124 = IT_0008*IT_0123;
    const ccomplex_t IT_0125 = m_b*conj(U_d2)*V_tb*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0126 = IT_0114*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + (-0.5)*IT_0127);
    const ccomplex_t IT_0129 = IT_0060*IT_0128;
    const ccomplex_t IT_0130 = IT_0052*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = m_C_1*IT_0131;
    const ccomplex_t IT_0133 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0134 = IT_0008*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(U_d1)*V_tb*e_em*IT_0013*U_sb_10;
    const ccomplex_t IT_0136 = IT_0114*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + (-0.5)*IT_0137);
    const ccomplex_t IT_0139 = 2*conj(U_sb_00)*(((-0.5)*m_W*IT_0041*(IT_0008 +
       1./3*IT_0007*IT_0043) + IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*U_sb_01 
      + 1./2*IT_0001*IT_0008*IT_0036*IT_0038*U_sb_11) + conj(U_sb_10)*(IT_0001
      *IT_0008*IT_0036*IT_0038*U_sb_01 + (-0.666666666666667)*(m_W*IT_0007
      *IT_0041*IT_0043 + (-3)*IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*U_sb_11);
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(e_em*IT_0139 +
       1.4142135623731*A_b*IT_0004*IT_0035*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = IT_0032*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 2*U_sb_01*(((-0.5)*m_W*IT_0041*(IT_0008 + 1./3
      *IT_0007*IT_0043) + IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*conj(U_sb_01)
       + 1./2*IT_0001*IT_0008*IT_0036*IT_0038*conj(U_sb_11)) + U_sb_11*(IT_0001
      *IT_0008*IT_0036*IT_0038*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0007
      *IT_0041*IT_0043 + (-3)*IT_0004*IT_0008*IT_0013*IT_0038*IT_0044)*conj
      (U_sb_11));
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(e_em*IT_0145 +
       1.4142135623731*A_b*IT_0004*IT_0035*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = IT_0118*IT_0147;
    const ccomplex_t IT_0149 = IT_0065*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0008*IT_0013*IT_0038;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0004*U_sb_01;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm1)*IT_0008;
    const ccomplex_t IT_0160 = IT_0076*IT_0159;
    const ccomplex_t IT_0161 = IT_0072*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0008*IT_0013*IT_0038;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0004*U_sb_01;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = IT_0165*IT_0167;
    const ccomplex_t IT_0169 = IT_0163*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0008;
    const ccomplex_t IT_0172 = IT_0028*IT_0171;
    const ccomplex_t IT_0173 = IT_0024*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = IT_0094*IT_0118;
    const ccomplex_t IT_0176 = IT_0000*IT_0021*IT_0175;
    const ccomplex_t IT_0177 = IT_0106*IT_0128;
    const ccomplex_t IT_0178 = IT_0051*IT_0052*IT_0177;
    const ccomplex_t IT_0179 = IT_0082*IT_0118;
    const ccomplex_t IT_0180 = IT_0080*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = m_t*IT_0181;
    const ccomplex_t IT_0183 = IT_0122 + IT_0132 + -IT_0144 + -IT_0150 + 
      -IT_0158 + -IT_0162 + -IT_0170 + -IT_0174 + IT_0176 + IT_0178 + 2*IT_0182;
    const ccomplex_t IT_0184 = 6*s_13;
    const ccomplex_t IT_0185 = m_t*m_C_1;
    const ccomplex_t IT_0186 = 6*IT_0111;
    const ccomplex_t IT_0187 = -IT_0085 + -IT_0121 + -IT_0131;
    const ccomplex_t IT_0188 = m_t*s_14;
    const ccomplex_t IT_0189 = IT_0097 + IT_0109 + IT_0181;
    const ccomplex_t IT_0190 = s_34*m_C_1;
    const ccomplex_t IT_0191 = 6*conj(IT_0111);
    const ccomplex_t IT_0192 = 6*conj(IT_0183);
    const ccomplex_t IT_0193 = pow(m_h, 2);
    const ccomplex_t IT_0194 = m_t*m_C_1*IT_0193;
    const ccomplex_t IT_0195 = s_13*IT_0193;
    const ccomplex_t IT_0196 = s_14*s_34;
    const ccomplex_t IT_0197 = (-2)*IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = 6*IT_0183;
    const ccomplex_t IT_0200 = 0.166666666666667*IT_0190;
    const ccomplex_t IT_0201 = (IT_0111*conj(IT_0111) + IT_0183*conj(IT_0183))
      *IT_0184 + conj(IT_0183)*IT_0185*IT_0186 + (IT_0183*IT_0185 + IT_0187
      *IT_0188 + IT_0189*IT_0190)*IT_0191 + (IT_0188*IT_0189 + IT_0187*IT_0190)
      *IT_0192 + 6*conj(IT_0189)*(IT_0187*IT_0194 + -IT_0189*IT_0198 +
       0.166666666666667*IT_0188*IT_0199 + IT_0186*IT_0200) + 6*conj(IT_0187)*
      (0.166666666666667*IT_0186*IT_0188 + IT_0189*IT_0194 + -IT_0187*IT_0198 +
       IT_0199*IT_0200);
    return create_ccomplex_return(IT_0201);
}

