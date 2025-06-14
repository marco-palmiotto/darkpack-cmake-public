#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_H0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_H0_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0002 = cpow(V_tb, 2);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = -IT_0008*IT_0011;
    const ccomplex_t IT_0013 = IT_0003*IT_0006;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0010 + 1.
      /3*IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj
      (U_sb_01) + -1./2*IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_11)) + U_sb_10
      *(IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0009*IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)
      *conj(U_sb_11));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0002*IT_0003*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_H0, 2);
    const ccomplex_t IT_0022 = pow(m_sb_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_sb_1, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0021 + IT_0026 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0020*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0025 + IT_0029;
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = IT_0001*IT_0025;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (-3)*s_23;
    const ccomplex_t IT_0035 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = (-3)*s_24;
    const ccomplex_t IT_0039 = (-3)*s_34;
    const ccomplex_t IT_0040 = (-3)*IT_0031*(IT_0022*conj(IT_0031) + (
      -0.333333333333333)*conj(IT_0033)*IT_0034 + (-0.333333333333333)*conj
      (IT_0037)*IT_0038) + IT_0037*(conj(IT_0031)*IT_0038 + conj(IT_0033)
      *IT_0039) + (-3)*IT_0033*(IT_0021*conj(IT_0033) + (-0.333333333333333)
      *conj(IT_0031)*IT_0034 + (-0.333333333333333)*conj(IT_0037)*IT_0039);
    return create_ccomplex_return(IT_0040);
}

