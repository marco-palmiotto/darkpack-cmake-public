#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_sc_R_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_sc_R_to_A0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_c, 2);
    const ccomplex_t IT_0009 = cpow(IT_0003, -2);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*IT_0004 + 
      -IT_0006*(IT_0002 + (-1.5)*IT_0007*IT_0008*IT_0009*IT_0011));
    const ccomplex_t IT_0013 = (-0.333333333333333)*IT_0012;
    const ccomplex_t IT_0014 = cos(alpha);
    const ccomplex_t IT_0015 = sin(alpha);
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cpow(IT_0003, -1);
    const ccomplex_t IT_0018 = cpow(IT_0010, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0010*
      (IT_0005*IT_0014 + -IT_0003*IT_0015) + 1.5*IT_0008*IT_0015*IT_0016*IT_0017
      *IT_0018);
    const ccomplex_t IT_0020 = (-0.666666666666667)*IT_0019;
    const ccomplex_t IT_0021 = IT_0002*IT_0010;
    const ccomplex_t IT_0022 = IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0003, 3);
    const ccomplex_t IT_0024 = -IT_0005;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0015
      *IT_0023 + IT_0024*(IT_0004*IT_0014 + (IT_0005*IT_0014 + -IT_0003*IT_0015)
      *IT_0024));
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = pow(m_sc_R, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + (-2)*IT_0028 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_c*e_em*mu_h*IT_0016*IT_0018;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0032;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = pow(m_A0, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0028 + IT_0036 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0010*
      (IT_0003*IT_0014 + IT_0005*IT_0015) + (-1.5)*IT_0008*IT_0014*IT_0016
      *IT_0017*IT_0018);
    const ccomplex_t IT_0041 = 0.666666666666667*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0005*
      (IT_0004*IT_0015 + -IT_0005*(IT_0003*IT_0014 + IT_0005*IT_0015)) + IT_0014
      *IT_0023);
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0028 + IT_0036 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0035*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0013 + -IT_0031 + IT_0039 + -IT_0047 +
       IT_0050;
    return create_ccomplex_return(3*IT_0051*conj(IT_0051));
}

