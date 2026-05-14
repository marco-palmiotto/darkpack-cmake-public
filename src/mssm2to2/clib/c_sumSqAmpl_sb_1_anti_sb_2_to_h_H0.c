#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_h_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(alpha);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = pow(m_W, -2);
    const ccomplex_t IT_0006 = pow(m_b, 2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0002*(0.5
      *(IT_0004 + 3*(1 + -IT_0005*IT_0006*IT_0008)*IT_0010)*U_sb_00*conj(U_sb_01
      ) + (IT_0004 + (-1.5)*IT_0005*IT_0006*IT_0008*IT_0010)*U_sb_10*conj
      (U_sb_11));
    const ccomplex_t IT_0012 = (-0.333333333333333)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_h, 2);
    const ccomplex_t IT_0014 = pow(m_sb_1, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_00;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = IT_0015*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = pow(m_sb_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_01);
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_00;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0002*conj(U_sb_01);
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0001*U_sb_00;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0013 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0041 = cpow(V_tb, 2);
    const ccomplex_t IT_0042 = pow(m_W, -1);
    const ccomplex_t IT_0043 = m_b*mu_h*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0007, -1);
    const ccomplex_t IT_0045 = cpow(IT_0009, -1);
    const ccomplex_t IT_0046 = sin(beta);
    const ccomplex_t IT_0047 = -IT_0002*IT_0046;
    const ccomplex_t IT_0048 = IT_0001*IT_0007;
    const ccomplex_t IT_0049 = IT_0047 + IT_0048;
    const ccomplex_t IT_0050 = (-2)*U_sb_00*((IT_0001*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0049)
      *conj(U_sb_01) + -1./2*IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + (-2)*(IT_0001
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0049)*conj
      (U_sb_11));
    const ccomplex_t IT_0051 = (0 + _Complex_I*-1)*(e_em*IT_0050 + (
      -1.4142135623731)*A_b*IT_0001*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0002*IT_0007;
    const ccomplex_t IT_0054 = IT_0001*IT_0046;
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = 2*U_sb_01*((IT_0002*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0055)
      *conj(U_sb_01) + 1./2*IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_11*(IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + 2*(IT_0002
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0055)*conj
      (U_sb_11));
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(e_em*IT_0056 +
       1.4142135623731*A_b*IT_0002*IT_0041*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0052*IT_0058;
    const ccomplex_t IT_0060 = IT_0040*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0002*U_sb_10;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0001*conj(U_sb_11);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_11);
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_10;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = IT_0070*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0002*conj(U_sb_11);
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0001*U_sb_10;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0080*IT_0082;
    const ccomplex_t IT_0084 = IT_0078*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_12 + -IT_0014 + -IT_0023 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0004*IT_0009;
    const ccomplex_t IT_0088 = IT_0045 + IT_0087;
    const ccomplex_t IT_0089 = cpow(IT_0002, 3);
    const ccomplex_t IT_0090 = cpow(IT_0002, 2);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_W*e_em*IT_0088*(IT_0007
      *IT_0089 + 5*IT_0001*((-1)*IT_0001*(IT_0002*IT_0007 + 1./5*IT_0001*IT_0046
      ) + IT_0046*IT_0090));
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = (-2)*U_sb_00*((IT_0001*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0049)
      *conj(U_sb_01) + -1./2*IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + (-2)*(IT_0001
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0049)*conj
      (U_sb_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_b*IT_0001*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0092*IT_0095;
    const ccomplex_t IT_0097 = IT_0086*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0014 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0100 = cpow(IT_0001, 2);
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_W*e_em*IT_0088*(IT_0046
      *IT_0089 + (-5)*IT_0001*(IT_0001*IT_0002*IT_0046 + 1*IT_0007*(IT_0090 + 
      -1./5*IT_0100)));
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = 2*U_sb_00*((IT_0002*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0055)
      *conj(U_sb_01) + 1./2*IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + 2*(IT_0002
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0055)*conj
      (U_sb_11));
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(e_em*IT_0103 +
       1.4142135623731*A_b*IT_0002*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0102*IT_0105;
    const ccomplex_t IT_0107 = IT_0099*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0013 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0110 = (-2)*U_sb_00*((IT_0001*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0049)
      *conj(U_sb_01) + -1./2*IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + (-2)*(IT_0001
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0049)*conj
      (U_sb_11));
    const ccomplex_t IT_0111 = (0 + _Complex_I*-1)*(e_em*IT_0110 + (
      -1.4142135623731)*A_b*IT_0001*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = 2*U_sb_00*((IT_0002*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0055)
      *conj(U_sb_00) + 1./2*IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_10)) +
       U_sb_10*(IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_00) + 2*(IT_0002
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0055)*conj
      (U_sb_10));
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(e_em*IT_0113 +
       1.4142135623731*A_b*IT_0002*IT_0041*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0112*IT_0115;
    const ccomplex_t IT_0117 = IT_0109*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_10;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_11);
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0002*U_sb_00;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0136 = (-2)*U_sb_00*((IT_0001*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0049)
      *conj(U_sb_00) + -1./2*IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_10)) +
       U_sb_10*(IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_00) + (-2)*(IT_0001
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0049)*conj
      (U_sb_10));
    const ccomplex_t IT_0137 = (0 + _Complex_I*-1)*(e_em*IT_0136 + (
      -1.4142135623731)*A_b*IT_0001*IT_0041*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = 2*U_sb_00*((IT_0002*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0055)
      *conj(U_sb_01) + 1./2*IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + 2*(IT_0002
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0055)*conj
      (U_sb_11));
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(e_em*IT_0139 +
       1.4142135623731*A_b*IT_0002*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = IT_0135*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0146 = (-2)*U_sb_01*((IT_0001*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0049)
      *conj(U_sb_01) + -1./2*IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_11*(IT_0002*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + (-2)*(IT_0001
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0049)*conj
      (U_sb_11));
    const ccomplex_t IT_0147 = (0 + _Complex_I*-1)*(e_em*IT_0146 + (
      -1.4142135623731)*A_b*IT_0001*IT_0041*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = 2*U_sb_00*((IT_0002*IT_0006*IT_0042*IT_0044
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0004*IT_0009 + 3*IT_0045)*IT_0055)
      *conj(U_sb_01) + 1./2*IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_11)) +
       U_sb_10*(IT_0001*IT_0043*IT_0044*IT_0045*conj(U_sb_01) + 2*(IT_0002
      *IT_0006*IT_0042*IT_0044*IT_0045 + -1./3*m_W*IT_0004*IT_0009*IT_0055)*conj
      (U_sb_11));
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(e_em*IT_0149 +
       1.4142135623731*A_b*IT_0002*IT_0041*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0148*IT_0151;
    const ccomplex_t IT_0153 = IT_0145*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0012 + IT_0022 + IT_0031 + IT_0039 + IT_0061
       + IT_0069 + IT_0077 + IT_0085 + -IT_0098 + -IT_0108 + IT_0118 + IT_0126 +
       IT_0134 + IT_0144 + IT_0154;
    return create_ccomplex_return(3*IT_0155*conj(IT_0155));
}

