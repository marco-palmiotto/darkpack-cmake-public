#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_R_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_R_to_h_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
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
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0002 = (-0.5)*IT_0001;
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_01;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = pow(m_h, 2);
    const ccomplex_t IT_0008 = pow(m_sd_R, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sd_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0008 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0003 + m_d*e_em*mu_h*IT_0000*IT_0014*IT_0016*IT_0018);
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = cos(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -2);
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = IT_0018 + IT_0023;
    const ccomplex_t IT_0025 = sin(beta);
    const ccomplex_t IT_0026 = cpow(IT_0003, 3);
    const ccomplex_t IT_0027 = cpow(IT_0000, 2);
    const ccomplex_t IT_0028 = cpow(IT_0003, 2);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_W*e_em*IT_0024*(IT_0025
      *IT_0026 + (-5)*IT_0000*(IT_0000*IT_0003*IT_0025 + (-0.2)*IT_0015*(IT_0027
       + (-5)*IT_0028)));
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = IT_0013*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_d, 2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0018 +
       1./3*IT_0017*IT_0022)*(IT_0000*IT_0015 + -IT_0003*IT_0025) + IT_0000
      *IT_0014*IT_0016*IT_0018*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0020*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*(m_W*IT_0017*IT_0022*
      (IT_0000*IT_0015 + -IT_0003*IT_0025) + (-3)*IT_0000*IT_0014*IT_0016
      *IT_0018*IT_0034);
    const ccomplex_t IT_0042 = 0.333333333333333*IT_0041;
    const ccomplex_t IT_0043 = IT_0020*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0007 + IT_0012 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0000 + -m_d*e_em*mu_h*IT_0003*IT_0014*IT_0016*IT_0018);
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0018 +
       1./3*IT_0017*IT_0022)*(IT_0003*IT_0015 + IT_0000*IT_0025) + IT_0003
      *IT_0014*IT_0016*IT_0018*IT_0034);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0007 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*(m_W*IT_0017*IT_0022*
      (IT_0003*IT_0015 + IT_0000*IT_0025) + (-3)*IT_0003*IT_0014*IT_0016*IT_0018
      *IT_0034);
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0007 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0007 + IT_0012 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_00;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = IT_0060*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0003*conj(U_sb_11);
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_01;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0007 + IT_0012 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + -IT_0008 + -IT_0012 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_W*e_em*IT_0024*(IT_0015
      *IT_0026 + 5*IT_0000*((-1)*IT_0000*(IT_0003*IT_0015 + 1./5*IT_0000*IT_0025
      ) + IT_0025*IT_0028));
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0048*IT_0078;
    const ccomplex_t IT_0080 = IT_0076*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_00;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = IT_0082*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0011 + -IT_0033 + IT_0040 + IT_0046 +
       IT_0053 + IT_0059 + IT_0067 + IT_0075 + -IT_0081 + IT_0089;
    return create_ccomplex_return(3*IT_0090*conj(IT_0090));
}

