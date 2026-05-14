#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_mu, 2);
    const ccomplex_t IT_0011 = cpow(IT_0001, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + -IT_0008*(IT_0004 + -IT_0006*(1 + (-2)*IT_0009*IT_0010*IT_0011
      )));
    const ccomplex_t IT_0013 = 0.25*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0005, -1);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0004*IT_0005;
    const ccomplex_t IT_0017 = IT_0015 + IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0007, 3);
    const ccomplex_t IT_0019 = sin(alpha);
    const ccomplex_t IT_0020 = -IT_0001;
    const ccomplex_t IT_0021 = 3*IT_0014 + IT_0016;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*(IT_0017*IT_0018
      *IT_0019 + IT_0020*(IT_0008*IT_0021*IT_0022 + -IT_0020*(IT_0007*IT_0019
      *IT_0021 + -IT_0001*IT_0017*IT_0022)));
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0014 + IT_0016;
    const ccomplex_t IT_0026 = (0 + _Complex_I*-1)*m_W*e_em*(IT_0007*IT_0019 +
       -IT_0001*IT_0022)*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + (-2)*IT_0029 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cpow(IT_0001, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007
      *IT_0014*(m_W*IT_0001 + (-0.5)*IT_0010*IT_0033*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0036, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = pow(m_Hp, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0029 + IT_0039 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0014*IT_0033;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0043, 2);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0029 + IT_0039 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0003, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*-1)*e_em*(IT_0003*IT_0014 + 
      -IT_0005*IT_0048);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0048;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = s_14 + s_24;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = s_12 + IT_0029;
    const ccomplex_t IT_0056 = 0.5*IT_0049;
    const ccomplex_t IT_0057 = -IT_0050;
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = IT_0055*IT_0058;
    const ccomplex_t IT_0060 = IT_0054 + IT_0059;
    const ccomplex_t IT_0061 = IT_0049*IT_0057;
    const ccomplex_t IT_0062 = s_24*IT_0061;
    const ccomplex_t IT_0063 = IT_0051*IT_0056;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = -IT_0062 + (-2)*IT_0064;
    const ccomplex_t IT_0066 = IT_0060 + IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = -IT_0014 + IT_0016;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_W*e_em*(IT_0018*IT_0022
      *IT_0070 + IT_0001*(IT_0008*IT_0019*IT_0021 + -IT_0001*(IT_0007*IT_0021
      *IT_0022 + IT_0001*IT_0019*IT_0070)));
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_W*e_em*(IT_0001*IT_0019 +
       IT_0007*IT_0022)*IT_0025;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0013 + -IT_0032 + IT_0041 + IT_0047 +
       IT_0069 + -IT_0078;
    return create_ccomplex_return(IT_0079*conj(IT_0079));
}

