#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_snu_e_to_anti_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_snu_e_to_anti_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = pow(m_e, 2);
    const ccomplex_t IT_0008 = cpow(IT_0005, -1);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0003*
      (m_W*IT_0005 + (-0.5)*IT_0006*IT_0007*IT_0008)*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_se_L, 2);
    const ccomplex_t IT_0013 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0004*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = pow(m_Z, -2);
    const ccomplex_t IT_0020 = pow(s_14, 2);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = -IT_0004;
    const ccomplex_t IT_0023 = pow(m_Hp, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0013 + IT_0023 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0022*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = s_14*s_34*IT_0019;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = s_13 + IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0003;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0012 + IT_0023 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = -IT_0017 + IT_0027 + IT_0037;
    const ccomplex_t IT_0039 = s_14*s_24*IT_0019;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = s_12 + IT_0040;
    const ccomplex_t IT_0042 = -IT_0018;
    const ccomplex_t IT_0043 = -conj(IT_0018);
    const ccomplex_t IT_0044 = pow(s_34, 2);
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = s_24*s_34*IT_0019;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = s_23 + IT_0047;
    const ccomplex_t IT_0049 = pow(s_24, 2);
    const ccomplex_t IT_0050 = IT_0019*IT_0049;
    const ccomplex_t IT_0051 = -IT_0038;
    const ccomplex_t IT_0052 = (conj(IT_0018)*(IT_0012 + -IT_0021) + conj
      (IT_0028)*IT_0031 + conj(IT_0038)*IT_0041)*IT_0042 + IT_0038*IT_0041
      *IT_0043 + IT_0028*(IT_0031*IT_0043 + -conj(IT_0028)*(IT_0023 + -IT_0045) 
      + -conj(IT_0038)*IT_0048) + (conj(IT_0028)*IT_0048 + conj(IT_0038)*
      (IT_0013 + -IT_0050))*IT_0051;
    return create_ccomplex_return(IT_0052);
}

