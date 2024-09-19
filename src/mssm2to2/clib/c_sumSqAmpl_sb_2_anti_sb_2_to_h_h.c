#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_h_h(
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
    const ccomplex_t IT_0001 = sin(alpha);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -2);
    const ccomplex_t IT_0008 = pow(m_W, -2);
    const ccomplex_t IT_0009 = pow(m_b, 2);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = 0.5*(IT_0005*(IT_0007 + 3*IT_0013) + IT_0003*
      (IT_0007 + 3*(1 + (-2)*IT_0008*IT_0009*IT_0011)*IT_0013))*U_sb_01*conj
      (U_sb_01) + ((IT_0003 + IT_0005)*IT_0007 + 3*IT_0002*IT_0008*IT_0009
      *IT_0011*IT_0013)*U_sb_11*conj(U_sb_11);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0000*IT_0014;
    const ccomplex_t IT_0016 = (-0.166666666666667)*IT_0015;
    const ccomplex_t IT_0017 = cpow(V_tb, 2);
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*mu_h*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0010, -1);
    const ccomplex_t IT_0021 = cpow(IT_0012, -1);
    const ccomplex_t IT_0022 = IT_0001*IT_0010;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = IT_0022 + IT_0024;
    const ccomplex_t IT_0026 = 2*conj(U_sb_00)*((IT_0001*IT_0009*IT_0018
      *IT_0020*IT_0021 + (-0.166666666666667)*m_W*(IT_0007*IT_0012 + 3*IT_0021)
      *IT_0025)*U_sb_01 + 1./2*IT_0004*IT_0019*IT_0020*IT_0021*U_sb_11) + conj
      (U_sb_10)*(IT_0004*IT_0019*IT_0020*IT_0021*U_sb_01 + 2*(IT_0001*IT_0009
      *IT_0018*IT_0020*IT_0021 + -1./3*m_W*IT_0007*IT_0012*IT_0025)*U_sb_11);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(e_em*IT_0026 +
       1.4142135623731*A_b*IT_0001*IT_0017*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = 2*U_sb_00*((IT_0001*IT_0009*IT_0018*IT_0020
      *IT_0021 + (-0.166666666666667)*m_W*(IT_0007*IT_0012 + 3*IT_0021)*IT_0025)
      *conj(U_sb_01) + 1./2*IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_11)) +
       U_sb_10*(IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_01) + 2*(IT_0001
      *IT_0009*IT_0018*IT_0020*IT_0021 + -1./3*m_W*IT_0007*IT_0012*IT_0025)*conj
      (U_sb_11));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(e_em*IT_0029 +
       1.4142135623731*A_b*IT_0001*IT_0017*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = pow(m_h, 2);
    const ccomplex_t IT_0034 = pow(m_sb_2, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0033 + IT_0034 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0039 = 2*U_sb_00*((IT_0001*IT_0009*IT_0018*IT_0020
      *IT_0021 + (-0.166666666666667)*m_W*(IT_0007*IT_0012 + 3*IT_0021)*IT_0025)
      *conj(U_sb_01) + 1./2*IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_11)) +
       U_sb_10*(IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_01) + 2*(IT_0001
      *IT_0009*IT_0018*IT_0020*IT_0021 + -1./3*m_W*IT_0007*IT_0012*IT_0025)*conj
      (U_sb_11));
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(e_em*IT_0039 +
       1.4142135623731*A_b*IT_0001*IT_0017*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0028*IT_0041;
    const ccomplex_t IT_0043 = IT_0038*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + (-2)*IT_0034 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0007*IT_0012;
    const ccomplex_t IT_0047 = IT_0021 + IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0001, 3);
    const ccomplex_t IT_0049 = (0 + _Complex_I*-5)*m_W*e_em*IT_0047*(IT_0004*
      (1*(IT_0002 + -1./5*IT_0005)*IT_0010 + IT_0001*IT_0004*IT_0023) + -1./5
      *IT_0023*IT_0048);
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = -IT_0001*IT_0023;
    const ccomplex_t IT_0052 = IT_0004*IT_0010;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = (-2)*U_sb_01*((IT_0004*IT_0009*IT_0018*IT_0020
      *IT_0021 + (-0.166666666666667)*m_W*(IT_0007*IT_0012 + 3*IT_0021)*IT_0053)
      *conj(U_sb_01) + -1./2*IT_0001*IT_0019*IT_0020*IT_0021*conj(U_sb_11)) +
       U_sb_11*(IT_0001*IT_0019*IT_0020*IT_0021*conj(U_sb_01) + (-2)*(IT_0004
      *IT_0009*IT_0018*IT_0020*IT_0021 + -1./3*m_W*IT_0007*IT_0012*IT_0053)*conj
      (U_sb_11));
    const ccomplex_t IT_0055 = (0 + _Complex_I*-1)*(e_em*IT_0054 + (
      -1.4142135623731)*A_b*IT_0004*IT_0017*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0050*IT_0056;
    const ccomplex_t IT_0058 = IT_0045*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0034 + -reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_W*e_em*IT_0047*(IT_0004*
      (IT_0002*IT_0023 + -IT_0004*(IT_0001*IT_0010 + IT_0004*IT_0023)) + IT_0010
      *IT_0048);
    const ccomplex_t IT_0062 = 1.5*IT_0061;
    const ccomplex_t IT_0063 = 2*U_sb_01*((IT_0001*IT_0009*IT_0018*IT_0020
      *IT_0021 + (-0.166666666666667)*m_W*(IT_0007*IT_0012 + 3*IT_0021)*IT_0025)
      *conj(U_sb_01) + 1./2*IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_11)) +
       U_sb_11*(IT_0004*IT_0019*IT_0020*IT_0021*conj(U_sb_01) + 2*(IT_0001
      *IT_0009*IT_0018*IT_0020*IT_0021 + -1./3*m_W*IT_0007*IT_0012*IT_0025)*conj
      (U_sb_11));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(e_em*IT_0063 +
       1.4142135623731*A_b*IT_0001*IT_0017*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0062*IT_0065;
    const ccomplex_t IT_0067 = IT_0060*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_11);
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_11;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0033 + IT_0034 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0073*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0001*conj(U_sb_11);
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0001*U_sb_11;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0033 + IT_0034 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0084*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_01;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_13 + IT_0033 + IT_0034 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0095*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0001*U_sb_01;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0033 + IT_0034 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0106*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0033 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0114 = cpow(IT_0065, 2);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_13 + IT_0033 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = IT_0016 + IT_0037 + IT_0044 + -IT_0059 + 
      -IT_0068 + IT_0076 + IT_0079 + IT_0087 + IT_0090 + IT_0098 + IT_0101 +
       IT_0109 + IT_0112 + IT_0116 + IT_0118;
    return create_ccomplex_return(3*IT_0119*conj(IT_0119));
}

