#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_A0_A0(
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
    const ccomplex_t IT_0000 = sin(beta);
    const ccomplex_t IT_0001 = 1.4142135623731*A_b*V_tb*V_ts*IT_0000*conj
      (U_sb_10);
    const ccomplex_t IT_0002 = 0.5*IT_0001;
    const ccomplex_t IT_0003 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0000
      *U_sb_10;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = pow(m_A0, 2);
    const ccomplex_t IT_0007 = pow(m_sb_1, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*A_b*V_tb*V_td*IT_0000*conj
      (U_sb_00);
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0000
      *U_sb_00;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*A_b*V_tb*V_ts*IT_0000*conj
      (U_sb_00);
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0000
      *U_sb_00;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0006 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0034 = sin(theta_W);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = pow(m_W, -1);
    const ccomplex_t IT_0037 = m_b*e_em*mu_h*IT_0036*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = m_b*e_em*mu_h*IT_0036*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0040 = IT_0035*IT_0039;
    const ccomplex_t IT_0041 = cpow(V_tb, 2);
    const ccomplex_t IT_0042 = A_b*IT_0041*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0043 = IT_0000*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = A_b*IT_0041*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0046 = IT_0000*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0038 + 0.5*IT_0040 + (-0.5)*IT_0044 +
       0.5*IT_0047;
    const ccomplex_t IT_0049 = cpow(IT_0048, 2);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0033*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0006 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = m_b*e_em*mu_h*IT_0036*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0055 = IT_0035*IT_0054;
    const ccomplex_t IT_0056 = m_b*e_em*mu_h*IT_0036*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0057 = IT_0035*IT_0056;
    const ccomplex_t IT_0058 = A_b*IT_0041*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0059 = IT_0000*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = A_b*IT_0041*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0062 = IT_0000*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0055 + (-0.5)*IT_0057 + 0.5*IT_0060 + (
      -0.5)*IT_0063;
    const ccomplex_t IT_0065 = m_b*e_em*mu_h*IT_0036*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0066 = IT_0035*IT_0065;
    const ccomplex_t IT_0067 = m_b*e_em*mu_h*IT_0036*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0068 = IT_0035*IT_0067;
    const ccomplex_t IT_0069 = A_b*IT_0041*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0070 = IT_0000*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = A_b*IT_0041*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0073 = IT_0000*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (-0.5)*IT_0066 + 0.5*IT_0068 + (-0.5)*IT_0071 +
       0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0076*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = cos(theta_W);
    const ccomplex_t IT_0084 = cpow(IT_0083, -2);
    const ccomplex_t IT_0085 = IT_0034*IT_0084;
    const ccomplex_t IT_0086 = IT_0035 + IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0000, 3);
    const ccomplex_t IT_0088 = sin(alpha);
    const ccomplex_t IT_0089 = cos(beta);
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = cos(alpha);
    const ccomplex_t IT_0092 = cpow(IT_0000, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_W*e_em*IT_0086*(IT_0087
      *IT_0088 + -IT_0090*(IT_0090*(IT_0000*IT_0088 + -IT_0089*IT_0091) + 
      -IT_0091*IT_0092));
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = m_b*mu_h*IT_0036;
    const ccomplex_t IT_0096 = cpow(IT_0089, -1);
    const ccomplex_t IT_0097 = -IT_0000*IT_0088;
    const ccomplex_t IT_0098 = IT_0089*IT_0091;
    const ccomplex_t IT_0099 = IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = pow(m_b, 2);
    const ccomplex_t IT_0101 = (-2)*U_sb_00*(((-0.5)*m_W*(IT_0035 + 1./3
      *IT_0034*IT_0084)*IT_0099 + IT_0035*IT_0036*IT_0091*IT_0096*IT_0100)*conj
      (U_sb_00) + -1./2*IT_0035*IT_0088*IT_0095*IT_0096*conj(U_sb_10)) + U_sb_10
      *(IT_0035*IT_0088*IT_0095*IT_0096*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0034*IT_0084*IT_0099 + (-3)*IT_0035*IT_0036*IT_0091*IT_0096*IT_0100)
      *conj(U_sb_10));
    const ccomplex_t IT_0102 = (0 + _Complex_I*-1)*(e_em*IT_0101 + (
      -1.4142135623731)*A_b*IT_0041*IT_0091*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0094*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + (-2)*IT_0007 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*A_b*V_tb*V_td*IT_0000*conj
      (U_sb_10);
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0000
      *U_sb_10;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0112*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0005*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_W*e_em*IT_0086*(IT_0087
      *IT_0091 + -IT_0089*(IT_0089*(IT_0088*IT_0089 + IT_0000*IT_0091) + 
      -IT_0088*IT_0092));
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0088*IT_0089;
    const ccomplex_t IT_0126 = IT_0000*IT_0091;
    const ccomplex_t IT_0127 = IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = 2*U_sb_00*((IT_0035*IT_0036*IT_0088*IT_0096
      *IT_0100 + (-0.5)*m_W*(IT_0035 + 1./3*IT_0034*IT_0084)*IT_0127)*conj
      (U_sb_00) + 1./2*IT_0035*IT_0091*IT_0095*IT_0096*conj(U_sb_10)) + U_sb_10*
      (IT_0035*IT_0091*IT_0095*IT_0096*conj(U_sb_00) + 2*(IT_0035*IT_0036
      *IT_0088*IT_0096*IT_0100 + -1./3*m_W*IT_0034*IT_0084*IT_0127)*conj(U_sb_10
      ));
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(e_em*IT_0128 +
       1.4142135623731*A_b*IT_0041*IT_0088*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0124*IT_0130;
    const ccomplex_t IT_0132 = IT_0122*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = pow(e_em, 2);
    const ccomplex_t IT_0135 = -IT_0092;
    const ccomplex_t IT_0136 = cpow(IT_0089, 2);
    const ccomplex_t IT_0137 = pow(m_W, -2);
    const ccomplex_t IT_0138 = cpow(IT_0089, -2);
    const ccomplex_t IT_0139 = cpow(IT_0034, -2);
    const ccomplex_t IT_0140 = 0.5*(IT_0136*(IT_0084 + 3*IT_0139) + IT_0135*
      (IT_0084 + 3*(1 + (-2)*IT_0100*IT_0137*IT_0138)*IT_0139))*U_sb_00*conj
      (U_sb_00) + (IT_0084*(IT_0135 + IT_0136) + 3*IT_0092*IT_0100*IT_0137
      *IT_0138*IT_0139)*U_sb_10*conj(U_sb_10);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0134*IT_0140;
    const ccomplex_t IT_0142 = (-0.166666666666667)*IT_0141;
    const ccomplex_t IT_0143 = IT_0010 + IT_0018 + IT_0021 + IT_0029 + IT_0032
       + IT_0051 + IT_0053 + IT_0079 + IT_0082 + -IT_0107 + IT_0115 + IT_0118 +
       IT_0121 + -IT_0133 + IT_0142;
    return create_ccomplex_return(3*IT_0143*conj(IT_0143));
}

