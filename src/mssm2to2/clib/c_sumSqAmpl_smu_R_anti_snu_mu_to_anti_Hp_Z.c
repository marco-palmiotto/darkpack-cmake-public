#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_snu_mu_to_anti_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_snu_mu_to_anti_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_smu_R, 2);
    const ccomplex_t IT_0009 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = pow(m_Z, -2);
    const ccomplex_t IT_0016 = pow(s_14, 2);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0019 = (-2)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_Hp, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0009 + IT_0020 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0019*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = s_14*s_34*IT_0015;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = s_13 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0003;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0008 + IT_0020 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0032 = IT_0007*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = -IT_0013 + IT_0024 + IT_0034;
    const ccomplex_t IT_0036 = s_14*s_24*IT_0015;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = s_12 + IT_0037;
    const ccomplex_t IT_0039 = -IT_0014;
    const ccomplex_t IT_0040 = -conj(IT_0014);
    const ccomplex_t IT_0041 = pow(s_34, 2);
    const ccomplex_t IT_0042 = IT_0015*IT_0041;
    const ccomplex_t IT_0043 = s_24*s_34*IT_0015;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = s_23 + IT_0044;
    const ccomplex_t IT_0046 = pow(s_24, 2);
    const ccomplex_t IT_0047 = IT_0015*IT_0046;
    const ccomplex_t IT_0048 = -IT_0035;
    const ccomplex_t IT_0049 = (conj(IT_0014)*(IT_0008 + -IT_0017) + conj
      (IT_0025)*IT_0028 + conj(IT_0035)*IT_0038)*IT_0039 + IT_0035*IT_0038
      *IT_0040 + IT_0025*(IT_0028*IT_0040 + -conj(IT_0025)*(IT_0020 + -IT_0042) 
      + -conj(IT_0035)*IT_0045) + (conj(IT_0025)*IT_0045 + conj(IT_0035)*
      (IT_0009 + -IT_0047))*IT_0048;
    return create_ccomplex_return(IT_0049);
}

