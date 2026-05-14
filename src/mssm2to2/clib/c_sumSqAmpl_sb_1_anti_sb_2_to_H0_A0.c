#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_H0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(V_tb, 2);
    const ccomplex_t IT_0009 = A_b*IT_0008*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = A_b*IT_0008*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0004 + 0.5*IT_0006 + (-0.5)*IT_0011 +
       0.5*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -2);
    const ccomplex_t IT_0018 = IT_0000*IT_0017;
    const ccomplex_t IT_0019 = IT_0001 + IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0007, 3);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = cos(alpha);
    const ccomplex_t IT_0025 = cpow(IT_0007, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_W*e_em*IT_0019*(IT_0020
      *IT_0021 + -IT_0023*(IT_0023*(IT_0007*IT_0021 + -IT_0022*IT_0024) + 
      -IT_0024*IT_0025));
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0015*IT_0027;
    const ccomplex_t IT_0029 = pow(m_sb_1, 2);
    const ccomplex_t IT_0030 = pow(m_sb_2, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + -IT_0029 + -IT_0030 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0007
      *U_sb_00;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0024*conj(U_sb_01);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = pow(m_H0, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0030 + IT_0039 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0007
      *U_sb_00;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0024*conj(U_sb_01);
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0030 + IT_0039 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0007
      *U_sb_10;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0024*conj(U_sb_11);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0030 + IT_0039 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0007
      *U_sb_10;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0024*conj(U_sb_11);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0030 + IT_0039 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0024*U_sb_10;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*A_b*V_tb*V_ts*IT_0007*conj
      (U_sb_11);
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0029 + IT_0039 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*-1)*m_W*e_em*IT_0007*IT_0019
      *IT_0022*(IT_0007*IT_0021 + -IT_0022*IT_0024);
    const ccomplex_t IT_0076 = cpow(IT_0022, -1);
    const ccomplex_t IT_0077 = IT_0001*IT_0007*IT_0076;
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = IT_0005*IT_0077;
    const ccomplex_t IT_0080 = IT_0009*IT_0022;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0012*IT_0022;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0078 + 0.5*IT_0079 + 0.5*IT_0081 + (
      -0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0075*IT_0084;
    const ccomplex_t IT_0086 = pow(m_Z, 2);
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + -IT_0029 + -IT_0030 + IT_0086 
      + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0029 + -IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0091 = cpow(IT_0016, -1);
    const ccomplex_t IT_0092 = -e_em*IT_0001*(IT_0021*IT_0022 + -IT_0007
      *IT_0024)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0001
      *IT_0016 + 1./3*IT_0000*IT_0091)*U_sb_00*conj(U_sb_01) + IT_0000*IT_0091
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0094 = (-0.333333333333333)*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = s_13 + s_23;
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0092;
    const ccomplex_t IT_0099 = (-0.666666666666667)*IT_0093;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = s_12 + IT_0030;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0102;
    const ccomplex_t IT_0104 = IT_0092*IT_0099;
    const ccomplex_t IT_0105 = s_23*IT_0104;
    const ccomplex_t IT_0106 = IT_0094*IT_0098;
    const ccomplex_t IT_0107 = 0.5*IT_0029 + 0.5*IT_0030;
    const ccomplex_t IT_0108 = s_12 + IT_0107;
    const ccomplex_t IT_0109 = IT_0106*IT_0108;
    const ccomplex_t IT_0110 = -IT_0105 + (-2)*IT_0109;
    const ccomplex_t IT_0111 = IT_0103 + IT_0110;
    const ccomplex_t IT_0112 = IT_0090*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_23 + IT_0030 + IT_0039 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0115 = m_b*e_em*mu_h*IT_0002*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0116 = IT_0001*IT_0115;
    const ccomplex_t IT_0117 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = A_b*IT_0008*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0120 = IT_0007*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = A_b*IT_0008*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0123 = IT_0007*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0116 + 0.5*IT_0118 + (-0.5)*IT_0121 +
       0.5*IT_0124;
    const ccomplex_t IT_0126 = m_b*mu_h*IT_0002;
    const ccomplex_t IT_0127 = -IT_0007*IT_0021;
    const ccomplex_t IT_0128 = IT_0022*IT_0024;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = pow(m_b, 2);
    const ccomplex_t IT_0131 = (-2)*U_sb_00*(((-0.5)*m_W*(IT_0001 + 1./3
      *IT_0000*IT_0017)*IT_0129 + IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)*conj
      (U_sb_01) + -1./2*IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_11)) + U_sb_10
      *(IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0000*IT_0017*IT_0129 + (-3)*IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)
      *conj(U_sb_11));
    const ccomplex_t IT_0132 = (0 + _Complex_I*-1)*(e_em*IT_0131 + (
      -1.4142135623731)*A_b*IT_0008*IT_0024*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = IT_0125*IT_0133;
    const ccomplex_t IT_0135 = IT_0114*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_13 + IT_0029 + IT_0039 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0138 = 1.4142135623731*A_b*V_tb*V_td*IT_0007*conj
      (U_sb_11);
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0024*U_sb_10;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = IT_0137*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0029 + IT_0039 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0146 = 1.4142135623731*A_b*V_tb*V_td*IT_0007*conj
      (U_sb_01);
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0024*U_sb_00;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147*IT_0149;
    const ccomplex_t IT_0151 = IT_0145*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_13 + IT_0029 + IT_0039 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0154 = 1.4142135623731*A_b*V_tb*V_ts*IT_0007*conj
      (U_sb_01);
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0024*U_sb_00;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = IT_0155*IT_0157;
    const ccomplex_t IT_0159 = IT_0153*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_13 + IT_0029 + IT_0039 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0162 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0163 = IT_0001*IT_0162;
    const ccomplex_t IT_0164 = m_b*e_em*mu_h*IT_0002*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0165 = IT_0001*IT_0164;
    const ccomplex_t IT_0166 = A_b*IT_0008*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0167 = IT_0007*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = A_b*IT_0008*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0170 = IT_0007*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0163 + 0.5*IT_0165 + (-0.5)*IT_0168 +
       0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0133*IT_0172;
    const ccomplex_t IT_0174 = IT_0161*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0177 = (-2)*U_sb_00*(((-0.5)*m_W*(IT_0001 + 1./3
      *IT_0000*IT_0017)*IT_0129 + IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)*conj
      (U_sb_00) + -1./2*IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_10)) + U_sb_10
      *(IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0000*IT_0017*IT_0129 + (-3)*IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)
      *conj(U_sb_10));
    const ccomplex_t IT_0178 = (0 + _Complex_I*-1)*(e_em*IT_0177 + (
      -1.4142135623731)*A_b*IT_0008*IT_0024*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0015*IT_0179;
    const ccomplex_t IT_0181 = IT_0176*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0184 = (-2)*U_sb_01*(((-0.5)*m_W*(IT_0001 + 1./3
      *IT_0000*IT_0017)*IT_0129 + IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)*conj
      (U_sb_01) + -1./2*IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_11)) + U_sb_11
      *(IT_0001*IT_0021*IT_0076*IT_0126*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0000*IT_0017*IT_0129 + (-3)*IT_0001*IT_0002*IT_0024*IT_0076*IT_0130)
      *conj(U_sb_11));
    const ccomplex_t IT_0185 = (0 + _Complex_I*-1)*(e_em*IT_0184 + (
      -1.4142135623731)*A_b*IT_0008*IT_0024*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0015*IT_0186;
    const ccomplex_t IT_0188 = IT_0183*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = -IT_0033 + IT_0042 + IT_0050 + IT_0058 +
       IT_0066 + IT_0074 + -IT_0089 + IT_0113 + IT_0136 + IT_0144 + IT_0152 +
       IT_0160 + IT_0175 + IT_0182 + IT_0189;
    return create_ccomplex_return(3*IT_0190*conj(IT_0190));
}

