#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_A0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
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
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = cpow(IT_0001, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 +
       3*IT_0006) + -IT_0008*(IT_0004 + 3*IT_0006*(1 + (-2)*IT_0009*IT_0010
      *IT_0011)));
    const ccomplex_t IT_0013 = (-0.0833333333333333)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_ss_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cpow(IT_0001, -1);
    const ccomplex_t IT_0018 = sin(alpha);
    const ccomplex_t IT_0019 = cpow(IT_0005, -1);
    const ccomplex_t IT_0020 = cos(alpha);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0010*IT_0016
      *IT_0017*IT_0018*IT_0019 + (-0.166666666666667)*m_W*(IT_0004*IT_0005 + 3
      *IT_0019)*(IT_0001*IT_0018 + IT_0007*IT_0020));
    const ccomplex_t IT_0022 = IT_0004*IT_0005;
    const ccomplex_t IT_0023 = IT_0019 + IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0007, 3);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_W*e_em*IT_0023*(IT_0001*
      (IT_0008*IT_0018 + -IT_0001*(IT_0001*IT_0018 + IT_0007*IT_0020)) + IT_0020
      *IT_0024);
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = IT_0015*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*(IT_0010*IT_0016
      *IT_0017*IT_0019*IT_0020 + 0.166666666666667*m_W*(IT_0004*IT_0005 + 3
      *IT_0019)*(IT_0007*IT_0018 + -IT_0001*IT_0020));
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = -IT_0001;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_W*e_em*IT_0023*(IT_0018
      *IT_0024 + IT_0032*(IT_0008*IT_0020 + -(IT_0007*IT_0018 + -IT_0001*IT_0020
      )*IT_0032));
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0031*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + (-2)*IT_0014 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0007;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0007;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = pow(m_A0, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0014 + IT_0044 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0014 + IT_0044 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0043*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = m_s*e_em*mu_h*IT_0016;
    const ccomplex_t IT_0052 = IT_0019*IT_0051;
    const ccomplex_t IT_0053 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0054 = IT_0007*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0052 + 0.5*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0052 + (-0.5)*IT_0055;
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0014 + IT_0044 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0014 + IT_0044 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0058*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*A_b*V_tb*V_ts*IT_0007*conj
      (U_sb_10);
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0007
      *U_sb_10;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0014 + IT_0044 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0014 + IT_0044 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0069*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0014 + IT_0044 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0077 = 1.4142135623731*A_b*V_tb*V_ts*IT_0007*conj
      (U_sb_11);
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0007
      *U_sb_11;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0076*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0014 + IT_0044 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0081*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0013 + -IT_0029 + -IT_0038 + IT_0047 +
       IT_0050 + IT_0061 + IT_0064 + IT_0072 + IT_0075 + IT_0083 + IT_0086;
    return create_ccomplex_return(3*IT_0087*conj(IT_0087));
}

