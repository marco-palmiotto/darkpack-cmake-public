#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_h_H0(
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
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(alpha);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = pow(m_W, -2);
    const ccomplex_t IT_0006 = pow(m_d, 2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0002*
      (IT_0004 + (-1.5)*IT_0005*IT_0006*IT_0008*IT_0010);
    const ccomplex_t IT_0012 = (-0.333333333333333)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_sd_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cpow(IT_0007, -1);
    const ccomplex_t IT_0018 = cpow(IT_0009, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0002*IT_0007 + IT_0001*IT_0015) + (-3)*IT_0002*IT_0006*IT_0016*IT_0017
      *IT_0018);
    const ccomplex_t IT_0020 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0021 = IT_0004*IT_0009;
    const ccomplex_t IT_0022 = IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0002, 3);
    const ccomplex_t IT_0024 = cpow(IT_0001, 2);
    const ccomplex_t IT_0025 = cpow(IT_0002, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0015
      *IT_0023 + (-5)*IT_0001*(IT_0001*IT_0002*IT_0015 + (-0.2)*IT_0007*(IT_0024
       + (-5)*IT_0025)));
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0020*IT_0027;
    const ccomplex_t IT_0029 = IT_0014*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0001 + -m_d*e_em*mu_h*IT_0002*IT_0016*IT_0017*IT_0018);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0002 + m_d*e_em*mu_h*IT_0001*IT_0016*IT_0017*IT_0018);
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = pow(m_h, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0013 + IT_0036 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0001;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0002;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0013 + IT_0036 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0013 + IT_0036 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0035*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0001;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0002;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_13 + IT_0013 + IT_0036 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0007 + -IT_0002*IT_0015) + (-3)*IT_0001*IT_0006*IT_0016
      *IT_0017*IT_0018);
    const ccomplex_t IT_0060 = 0.333333333333333*IT_0059;
    const ccomplex_t IT_0061 = IT_0020*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0036 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0036 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0061*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_00);
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_00;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0013 + IT_0036 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_01);
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0001*U_sb_01;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_13 + IT_0013 + IT_0036 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_01;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0013 + IT_0036 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + (-2)*IT_0013 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0007
      *IT_0023 + 5*IT_0001*((-1)*IT_0001*(IT_0002*IT_0007 + 1./5*IT_0001*IT_0015
      ) + IT_0015*IT_0025));
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0060*IT_0094;
    const ccomplex_t IT_0096 = IT_0092*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0001*conj(U_sb_00);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_00;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0013 + IT_0036 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0012 + -IT_0030 + IT_0039 + IT_0047 +
       IT_0050 + IT_0058 + IT_0064 + IT_0067 + IT_0075 + IT_0083 + IT_0091 + 
      -IT_0097 + IT_0105;
    return create_ccomplex_return(3*IT_0106*conj(IT_0106));
}

