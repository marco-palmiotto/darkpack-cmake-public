#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_h_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_h_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0000 = pow(m_h, 2);
    const ccomplex_t IT_0001 = pow(m_sb_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cpow(V_tb, 2);
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = m_b*mu_h*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0004*IT_0007;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = IT_0012 + IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -2);
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = 2*U_sb_00*(((-0.5)*m_W*IT_0015*(IT_0011 + 1./3
      *IT_0010*IT_0017) + IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*conj(U_sb_01)
       + 1./2*IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_11)) + U_sb_10*(IT_0006
      *IT_0008*IT_0009*IT_0011*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0010
      *IT_0015*IT_0017 + (-3)*IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*conj
      (U_sb_11));
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(e_em*IT_0019 +
       1.4142135623731*A_b*IT_0003*IT_0004*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_sb_1, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0023 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0030 = 0.333333333333333*IT_0029;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0029;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = (-3)*s_24;
    const ccomplex_t IT_0036 = IT_0023*IT_0033;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = (-3)*s_34;
    const ccomplex_t IT_0039 = (-0.333333333333333)*conj(IT_0032);
    const ccomplex_t IT_0040 = IT_0032*(conj(IT_0034)*IT_0035 + conj(IT_0037)
      *IT_0038) + (-3)*IT_0034*(IT_0001*conj(IT_0034) + s_23*conj(IT_0037) +
       IT_0035*IT_0039) + (-3)*IT_0037*(s_23*conj(IT_0034) + IT_0000*conj
      (IT_0037) + IT_0038*IT_0039);
    return create_ccomplex_return(IT_0040);
}

