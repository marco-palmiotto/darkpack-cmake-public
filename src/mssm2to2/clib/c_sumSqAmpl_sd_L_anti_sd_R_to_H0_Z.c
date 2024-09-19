#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_R_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_R_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + IT_0004*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = m_d*e_em*mu_h*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0002, -1);
    const ccomplex_t IT_0013 = IT_0004*IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0016 = IT_0002*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0014 + (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_Z, 2);
    const ccomplex_t IT_0020 = pow(m_sd_L, 2);
    const ccomplex_t IT_0021 = pow(m_sd_R, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + IT_0019 + -IT_0020 + -IT_0021 
      + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0006 + 3
      *IT_0000*IT_0007);
    const ccomplex_t IT_0027 = 0.333333333333333*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0003 + -m_d*e_em*mu_h*IT_0005*IT_0007*IT_0010*IT_0012);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = pow(m_H0, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0021 + IT_0030 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = e_em*IT_0001*(IT_0003*IT_0004 + -IT_0002
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0036 = IT_0007*IT_0011;
    const ccomplex_t IT_0037 = IT_0004*IT_0015;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0036 + 0.5*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = IT_0025 + -IT_0034 + IT_0043;
    const ccomplex_t IT_0045 = pow(m_Z, -2);
    const ccomplex_t IT_0046 = pow(s_34, 2);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0049 = (-0.666666666666667)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0020 + IT_0030 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0029*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0049*IT_0052;
    const ccomplex_t IT_0054 = IT_0034 + IT_0053;
    const ccomplex_t IT_0055 = s_24*s_34*IT_0045;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = s_23 + IT_0056;
    const ccomplex_t IT_0058 = pow(s_24, 2);
    const ccomplex_t IT_0059 = IT_0045*IT_0058;
    return create_ccomplex_return((-3)*IT_0044*(conj(IT_0044)*(IT_0030 + 
      -IT_0047) + conj(IT_0054)*IT_0057) + (-3)*IT_0054*(conj(IT_0044)*IT_0057 +
       conj(IT_0054)*(IT_0021 + -IT_0059)));
}

