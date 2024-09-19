#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_H0_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_H0_H0(
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
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -2);
    const ccomplex_t IT_0009 = pow(m_b, 2);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = 0.5*(IT_0004*(IT_0002 + 3*IT_0013) + IT_0007*
      (IT_0002 + 3*(1 + (-2)*IT_0008*IT_0009*IT_0011)*IT_0013))*U_sb_01*conj
      (U_sb_01) + (IT_0002*IT_0004 + IT_0007*(IT_0002 + (-3)*IT_0008*IT_0009
      *IT_0011*IT_0013))*U_sb_11*conj(U_sb_11);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0000*IT_0014;
    const ccomplex_t IT_0016 = (-0.166666666666667)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_sb_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0002*IT_0012;
    const ccomplex_t IT_0020 = cpow(IT_0012, -1);
    const ccomplex_t IT_0021 = IT_0019 + IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0003, 3);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_W*e_em*IT_0021*(IT_0010
      *IT_0022 + 5*IT_0005*(IT_0004*IT_0023 + (-1)*IT_0005*(IT_0003*IT_0010 + 1.
      /5*IT_0005*IT_0023)));
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = cpow(V_tb, 2);
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = m_b*mu_h*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0010, -1);
    const ccomplex_t IT_0030 = IT_0003*IT_0010;
    const ccomplex_t IT_0031 = IT_0005*IT_0023;
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = 2*U_sb_01*((IT_0003*IT_0009*IT_0020*IT_0027
      *IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + 3*IT_0020)*IT_0032)
      *conj(U_sb_01) + 1./2*IT_0005*IT_0020*IT_0028*IT_0029*conj(U_sb_11)) +
       U_sb_11*(IT_0005*IT_0020*IT_0028*IT_0029*conj(U_sb_01) + 2*(IT_0003
      *IT_0009*IT_0020*IT_0027*IT_0029 + -1./3*m_W*IT_0002*IT_0012*IT_0032)*conj
      (U_sb_11));
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(e_em*IT_0033 +
       1.4142135623731*A_b*IT_0003*IT_0026*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = IT_0018*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0005*conj(U_sb_11);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0005*U_sb_11;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = pow(m_H0, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0017 + IT_0044 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0017 + IT_0044 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0043*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0005*conj(U_sb_11);
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0005*U_sb_11;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_13 + IT_0017 + IT_0044 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0017 + IT_0044 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0055*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0005*conj(U_sb_01);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0005*U_sb_01;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0017 + IT_0044 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0017 + IT_0044 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0066*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0005*conj(U_sb_01);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0005*U_sb_01;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0017 + IT_0044 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0017 + IT_0044 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0077*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0044 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0085 = -IT_0003*IT_0023;
    const ccomplex_t IT_0086 = IT_0005*IT_0010;
    const ccomplex_t IT_0087 = IT_0085 + IT_0086;
    const ccomplex_t IT_0088 = (-2)*U_sb_01*((IT_0005*IT_0009*IT_0020*IT_0027
      *IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + 3*IT_0020)*IT_0087)
      *conj(U_sb_01) + -1./2*IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_11)) +
       U_sb_11*(IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_01) + (-2)*(IT_0005
      *IT_0009*IT_0020*IT_0027*IT_0029 + -1./3*m_W*IT_0002*IT_0012*IT_0087)*conj
      (U_sb_11));
    const ccomplex_t IT_0089 = (0 + _Complex_I*-1)*(e_em*IT_0088 + (
      -1.4142135623731)*A_b*IT_0005*IT_0026*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = cpow(IT_0090, 2);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0084*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0044 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0017 + IT_0044 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0097 = (-2)*conj(U_sb_00)*((IT_0005*IT_0009*IT_0020
      *IT_0027*IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + 3*IT_0020)
      *IT_0087)*U_sb_01 + -1./2*IT_0003*IT_0020*IT_0028*IT_0029*U_sb_11) + conj
      (U_sb_10)*(IT_0003*IT_0020*IT_0028*IT_0029*U_sb_01 + (-2)*(IT_0005*IT_0009
      *IT_0020*IT_0027*IT_0029 + -1./3*m_W*IT_0002*IT_0012*IT_0087)*U_sb_11);
    const ccomplex_t IT_0098 = (0 + _Complex_I*-1)*(e_em*IT_0097 + (
      -1.4142135623731)*A_b*IT_0005*IT_0026*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = (-2)*U_sb_00*((IT_0005*IT_0009*IT_0020*IT_0027
      *IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + 3*IT_0020)*IT_0087)
      *conj(U_sb_01) + -1./2*IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_11)) +
       U_sb_10*(IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_01) + (-2)*(IT_0005
      *IT_0009*IT_0020*IT_0027*IT_0029 + -1./3*m_W*IT_0002*IT_0012*IT_0087)*conj
      (U_sb_11));
    const ccomplex_t IT_0101 = (0 + _Complex_I*-1)*(e_em*IT_0100 + (
      -1.4142135623731)*A_b*IT_0005*IT_0026*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0099*IT_0102;
    const ccomplex_t IT_0104 = IT_0096*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_W*e_em*IT_0021*(IT_0022
      *IT_0023 + -IT_0005*((IT_0004 + -IT_0006)*IT_0010 + IT_0003*IT_0005
      *IT_0023));
    const ccomplex_t IT_0108 = (-1.5)*IT_0107;
    const ccomplex_t IT_0109 = (-2)*U_sb_01*((IT_0005*IT_0009*IT_0020*IT_0027
      *IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + 3*IT_0020)*IT_0087)
      *conj(U_sb_01) + -1./2*IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_11)) +
       U_sb_11*(IT_0003*IT_0020*IT_0028*IT_0029*conj(U_sb_01) + (-2)*(IT_0005
      *IT_0009*IT_0020*IT_0027*IT_0029 + -1./3*m_W*IT_0002*IT_0012*IT_0087)*conj
      (U_sb_11));
    const ccomplex_t IT_0110 = (0 + _Complex_I*-1)*(e_em*IT_0109 + (
      -1.4142135623731)*A_b*IT_0005*IT_0026*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = IT_0106*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0017 + IT_0044 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0103*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0016 + -IT_0038 + IT_0047 + IT_0050 +
       IT_0058 + IT_0061 + IT_0069 + IT_0072 + IT_0080 + IT_0083 + IT_0093 +
       IT_0095 + IT_0105 + -IT_0114 + IT_0117;
    return create_ccomplex_return(3*IT_0118*conj(IT_0118));
}

