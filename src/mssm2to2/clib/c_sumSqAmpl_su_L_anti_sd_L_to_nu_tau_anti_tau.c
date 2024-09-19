#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_L_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_L_to_nu_tau_anti_tau(
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
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_tau = param->m_tau;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = pow(m_sd_L, 2);
    const ccomplex_t IT_0005 = pow(m_su_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0004 + -IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = m_tau*IT_0010;
    const ccomplex_t IT_0012 = -IT_0009;
    const ccomplex_t IT_0013 = m_tau*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0011 + IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0001*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_d, 2);
    const ccomplex_t IT_0021 = pow(m_u, 2);
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, 2);
    const ccomplex_t IT_0024 = sin(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001*(IT_0017*(IT_0020 + -IT_0021) + -m_W*(IT_0023 + -IT_0025));
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0019*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + IT_0029 
      + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0034 = cpow(IT_0022, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0001*IT_0017*IT_0024*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0024, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001*(IT_0017*IT_0020*IT_0024*IT_0034 + (-2)*IT_0022*(m_W*IT_0024 + (
      -0.5)*IT_0017*IT_0021*IT_0037));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = IT_0033*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0016 + -IT_0032 + -IT_0042;
    const ccomplex_t IT_0044 = IT_0008*IT_0012;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = s_13*m_tau;
    const ccomplex_t IT_0047 = s_13*s_14;
    const ccomplex_t IT_0048 = s_34*IT_0005;
    return create_ccomplex_return(6*IT_0043*(s_34*conj(IT_0043) + -conj
      (IT_0045)*IT_0046) + (-6)*IT_0045*(conj(IT_0043)*IT_0046 + (
      -0.166666666666667)*conj(IT_0045)*(12*IT_0047 + (-6)*IT_0048)));
}

