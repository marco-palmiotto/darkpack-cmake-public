#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sd_R_to_u_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sd_R_to_u_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = pow(m_sd_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0005 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0005*IT_0018;
    const ccomplex_t IT_0020 = 2*IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = s_14*IT_0000;
    const ccomplex_t IT_0022 = IT_0004*IT_0012;
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = s_14*s_34;
    const ccomplex_t IT_0027 = -IT_0019 + (-2)*IT_0025;
    const ccomplex_t IT_0028 = s_34*IT_0001;
    const ccomplex_t IT_0029 = s_13*s_14;
    const ccomplex_t IT_0030 = (-12)*IT_0029;
    const ccomplex_t IT_0031 = 6*IT_0028 + IT_0030;
    const ccomplex_t IT_0032 = s_13*s_34;
    const ccomplex_t IT_0033 = (-2)*IT_0032;
    const ccomplex_t IT_0034 = IT_0021 + IT_0033;
    const ccomplex_t IT_0035 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0036 = IT_0018*IT_0035;
    const ccomplex_t IT_0037 = (-6)*s_13;
    const ccomplex_t IT_0038 = s_34*IT_0037;
    const ccomplex_t IT_0039 = pow(s_13, 2);
    const ccomplex_t IT_0040 = (-6)*IT_0039;
    const ccomplex_t IT_0041 = (-6)*conj(IT_0015)*(IT_0020*IT_0021 + (
      -0.166666666666667)*(24*IT_0015 + (-24)*IT_0025)*IT_0026 + (
      -0.166666666666667)*IT_0027*IT_0031) + (-24)*conj(IT_0025)*((IT_0015 + 
      -IT_0025)*IT_0026 + 1./4*IT_0027*IT_0028 + -1./4*IT_0020*IT_0034) + conj
      (IT_0036)*(s_14*IT_0027*IT_0037 + IT_0020*IT_0038) + (-6)*conj(IT_0027)*
      (IT_0025*IT_0028 + (-0.166666666666667)*IT_0015*IT_0031 + (
      -0.166666666666667)*(IT_0001*IT_0027 + s_14*IT_0036)*IT_0037 + (
      -0.166666666666667)*IT_0020*IT_0040) + (-6)*conj(IT_0020)*(s_13*IT_0000
      *IT_0020 + IT_0015*IT_0021 + -IT_0025*IT_0034 + (-0.166666666666667)
      *IT_0036*IT_0038 + (-0.166666666666667)*IT_0027*IT_0040);
    return create_ccomplex_return(IT_0041);
}

