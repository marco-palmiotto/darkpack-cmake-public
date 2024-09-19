#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_R_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_R_to_h_h(
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
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(alpha);
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0000 + m_s*e_em*mu_h*IT_0001*IT_0003*IT_0004*IT_0006);
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = pow(m_s, 2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0006*IT_0009 + (-0.5)*m_W*(IT_0000*IT_0002 + IT_0004*IT_0010)*
      (IT_0006 + 1./3*IT_0005*IT_0012));
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0015 + (0 + _Complex_I*1)
      *m_ss_L*Gamma_sl + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_ss_R, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0015 + IT_0019 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0014*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0006*IT_0009 + -1./3*m_W*IT_0005*(IT_0000*IT_0002 + IT_0004
      *IT_0010)*IT_0012);
    const ccomplex_t IT_0024 = (-0.333333333333333)*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0015 + (0 + _Complex_I*1)
      *m_ss_R*Gamma_sr + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_ss_L, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0015 + IT_0029 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0015 + IT_0019 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0000*U_sb_01;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0015 + IT_0029 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0004 + -m_s*e_em*mu_h*IT_0000*IT_0001*IT_0003*IT_0006);
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0005*IT_0012;
    const ccomplex_t IT_0047 = IT_0006 + IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0000, 3);
    const ccomplex_t IT_0049 = cpow(IT_0004, 2);
    const ccomplex_t IT_0050 = cpow(IT_0000, 2);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_W*e_em*IT_0047*(IT_0010
      *IT_0048 + (-5)*IT_0004*(IT_0000*IT_0004*IT_0010 + (-0.2)*IT_0002*(IT_0049
       + (-5)*IT_0050)));
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0045*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + -IT_0019 + -IT_0029 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_W*e_em*IT_0047*(IT_0002
      *IT_0048 + -IT_0004*(IT_0004*(IT_0000*IT_0002 + IT_0004*IT_0010) + 
      -IT_0010*IT_0050));
    const ccomplex_t IT_0058 = 1.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0008*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0019 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0015 + IT_0019 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0000*U_sb_00;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0015 + IT_0029 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0018 + IT_0022 + IT_0028 + IT_0032 + IT_0040
       + IT_0043 + -IT_0056 + -IT_0062 + IT_0070 + IT_0073;
    return create_ccomplex_return(3*IT_0074*conj(IT_0074));
}

