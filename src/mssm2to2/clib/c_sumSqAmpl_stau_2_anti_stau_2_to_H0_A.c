#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_H0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_H0_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = mu_h*m_tau*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = -IT_0007*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0005;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = (-2)*U_stau_01*(((-0.5)*m_W*IT_0013*(IT_0009 + 
      -IT_0008*IT_0015) + IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj
      (U_stau_01) + -1./2*IT_0004*IT_0006*IT_0007*IT_0009*conj(U_stau_11)) +
       U_stau_11*(IT_0004*IT_0006*IT_0007*IT_0009*conj(U_stau_01) + 2*(m_W
      *IT_0008*IT_0013*IT_0015 + -IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj
      (U_stau_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0002*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_H0, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0020 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0020 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0023 + IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0027;
    const ccomplex_t IT_0029 = pow(m_stau_2, 2);
    const ccomplex_t IT_0030 = IT_0001*IT_0023;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = -s_23;
    const ccomplex_t IT_0033 = IT_0000*IT_0027;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = -s_24;
    const ccomplex_t IT_0036 = -s_34;
    const ccomplex_t IT_0037 = -IT_0028*(conj(IT_0028)*IT_0029 + -conj(IT_0031
      )*IT_0032 + -conj(IT_0034)*IT_0035) + IT_0034*(conj(IT_0028)*IT_0035 +
       conj(IT_0031)*IT_0036) + -IT_0031*(IT_0020*conj(IT_0031) + -conj(IT_0028)
      *IT_0032 + -conj(IT_0034)*IT_0036);
    return create_ccomplex_return(IT_0037);
}

