#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_sc_R_to_H0_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_sc_R_to_H0_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_c, 2);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + 1.5*IT_0006*IT_0007*IT_0008*IT_0010*IT_0012);
    const ccomplex_t IT_0014 = (-0.333333333333333)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = cpow(IT_0009, -1);
    const ccomplex_t IT_0017 = cpow(IT_0011, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0003
      *IT_0015*IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_H0, 2);
    const ccomplex_t IT_0023 = pow(m_sc_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0022 + (0 + _Complex_I*1)
      *m_sc_R*Gamma_cr + reg_prop, -1);
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*(1.5*IT_0005*IT_0008
      *IT_0015*IT_0016*IT_0017 + -m_W*IT_0002*IT_0011*(IT_0005*IT_0009 + 
      -IT_0003*IT_0027));
    const ccomplex_t IT_0029 = (-0.666666666666667)*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0022 + (0 + _Complex_I*1)
      *m_sc_R*Gamma_cr + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*e_em*((-1.5)*IT_0003*IT_0008
      *IT_0015*IT_0016*IT_0017 + m_W*IT_0002*IT_0011*(IT_0003*IT_0009 + IT_0005
      *IT_0027));
    const ccomplex_t IT_0038 = 0.666666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0002*IT_0011;
    const ccomplex_t IT_0040 = IT_0017 + IT_0039;
    const ccomplex_t IT_0041 = cpow(IT_0005, 3);
    const ccomplex_t IT_0042 = (0 + _Complex_I*5)*m_W*e_em*IT_0040*(IT_0003*
      (IT_0006*IT_0009 + (-0.2)*IT_0003*(IT_0003*IT_0009 + 5*IT_0005*IT_0027)) +
       1./5*IT_0027*IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*-1)*m_W*e_em*IT_0040*(IT_0003*
      (IT_0003*IT_0005*IT_0009 + -(IT_0004 + -IT_0006)*IT_0027) + -IT_0009
      *IT_0041);
    const ccomplex_t IT_0049 = (-1.5)*IT_0048;
    const ccomplex_t IT_0050 = IT_0029*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0014 + IT_0025 + IT_0032 + IT_0034 + IT_0036
       + -IT_0047 + -IT_0053;
    return create_ccomplex_return(3*IT_0054*conj(IT_0054));
}

