#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_e, 2);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + IT_0006*IT_0007*IT_0008*IT_0010*IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = cpow(IT_0009, -1);
    const ccomplex_t IT_0017 = cpow(IT_0011, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0003
      *IT_0015*IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_h, 2);
    const ccomplex_t IT_0023 = pow(m_se_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0022 + (0 + _Complex_I*1)
      *m_se_R*Gamma_er + reg_prop, -1);
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = (0 + _Complex_I*-1)*e_em*(IT_0005*IT_0008
      *IT_0015*IT_0016*IT_0017 + -m_W*IT_0002*IT_0011*(IT_0005*IT_0009 + IT_0003
      *IT_0027));
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0022 + (0 + _Complex_I*1)
      *m_se_R*Gamma_er + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0002*IT_0011;
    const ccomplex_t IT_0037 = IT_0017 + IT_0036;
    const ccomplex_t IT_0038 = cpow(IT_0005, 3);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_W*e_em*IT_0037*(IT_0003*
      (IT_0006*IT_0027 + -IT_0003*(IT_0005*IT_0009 + IT_0003*IT_0027)) + IT_0009
      *IT_0038);
    const ccomplex_t IT_0040 = 1.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*-1)*e_em*(IT_0003*IT_0008
      *IT_0015*IT_0016*IT_0017 + -m_W*IT_0002*IT_0011*(IT_0003*IT_0009 + 
      -IT_0005*IT_0027));
    const ccomplex_t IT_0046 = (0 + _Complex_I*-5)*m_W*e_em*IT_0037*(IT_0003*(
      (-0.2)*(IT_0004 + (-5)*IT_0006)*IT_0009 + IT_0003*IT_0005*IT_0027) + -1./5
      *IT_0027*IT_0038);
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0044*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0021*IT_0051;
    const ccomplex_t IT_0053 = IT_0014 + IT_0025 + IT_0032 + IT_0034 + 
      -IT_0043 + -IT_0050 + IT_0052;
    return create_ccomplex_return(IT_0053*conj(IT_0053));
}

