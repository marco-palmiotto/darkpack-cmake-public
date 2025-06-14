#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_anti_nu_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_anti_nu_e_nu_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_se_R, 2);
    const ccomplex_t IT_0011 = cpow(s_23 + (-0.5)*IT_0010 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = e_em*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0014;
    const ccomplex_t IT_0019 = e_em*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0017 + IT_0019);
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0015;
    const ccomplex_t IT_0023 = (-2)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_23 + (-0.5)*IT_0010 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0013 + -IT_0027 + 0.5*IT_0035;
    const ccomplex_t IT_0037 = s_13*s_14;
    const ccomplex_t IT_0038 = s_34*IT_0010;
    return create_ccomplex_return(4*IT_0036*conj(IT_0036)*(IT_0037 + (-0.5)
      *IT_0038));
}

