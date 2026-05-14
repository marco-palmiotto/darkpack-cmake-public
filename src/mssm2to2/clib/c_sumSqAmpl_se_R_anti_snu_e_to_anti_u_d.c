#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_snu_e_to_anti_u_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_snu_e_to_anti_u_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_34 = param->s_34;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0000*IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, 2);
    const ccomplex_t IT_0012 = pow(m_se_R, 2);
    const ccomplex_t IT_0013 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0013 
      + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0000*IT_0002;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0000*IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0016 + -IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0000*IT_0002;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0027;
    const ccomplex_t IT_0029 = IT_0014*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0007, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0000*IT_0002*IT_0005*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = IT_0022*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0030 + -IT_0036;
    const ccomplex_t IT_0038 = m_d*m_u;
    const ccomplex_t IT_0039 = 6*s_34;
    return create_ccomplex_return((-6)*conj(IT_0025)*(IT_0037*IT_0038 + (
      -0.166666666666667)*IT_0025*IT_0039) + (-6)*conj(IT_0037)*(IT_0025*IT_0038
       + (-0.166666666666667)*IT_0037*IT_0039));
}

