#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_snu_e_to_anti_t_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_snu_e_to_anti_t_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_e = param->m_e;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_se_L = param->m_se_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = pow(m_se_L, 2);
    const ccomplex_t IT_0007 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = pow(m_e, 2);
    const ccomplex_t IT_0016 = cpow(IT_0013, -1);
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001*
      (m_W*IT_0013 + (-0.5)*IT_0014*IT_0015*IT_0016)*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0014*IT_0016*IT_0017;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0014;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0013, 2);
    const ccomplex_t IT_0029 = cpow(IT_0017, 2);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001*
      (IT_0014*IT_0015 + -m_W*(IT_0028 + -IT_0029));
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0027*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, 2);
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + IT_0033 
      + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0002;
    const ccomplex_t IT_0038 = m_b*IT_0037;
    const ccomplex_t IT_0039 = m_b*IT_0003;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0038 + IT_0040;
    const ccomplex_t IT_0042 = IT_0010*IT_0041;
    const ccomplex_t IT_0043 = -IT_0025 + -IT_0036 + -IT_0042;
    const ccomplex_t IT_0044 = m_b*s_13;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0014;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0031*IT_0046;
    const ccomplex_t IT_0048 = IT_0034*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0017, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0013*IT_0014*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0019*IT_0052;
    const ccomplex_t IT_0054 = IT_0023*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_t*IT_0037;
    const ccomplex_t IT_0057 = m_t*IT_0003;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = IT_0010*IT_0059;
    const ccomplex_t IT_0061 = -IT_0049 + -IT_0055 + IT_0060;
    const ccomplex_t IT_0062 = m_t*s_14;
    const ccomplex_t IT_0063 = s_13*s_14;
    const ccomplex_t IT_0064 = s_34*IT_0006;
    const ccomplex_t IT_0065 = m_b*m_t;
    const ccomplex_t IT_0066 = 6*s_34;
    const ccomplex_t IT_0067 = 6*conj(IT_0012)*(IT_0043*IT_0044 + -IT_0061
      *IT_0062 + 0.166666666666667*IT_0012*(12*IT_0063 + (-6)*IT_0064)) + 6*conj
      (IT_0043)*(IT_0012*IT_0044 + -IT_0061*IT_0065 + 0.166666666666667*IT_0043
      *IT_0066) + (-6)*conj(IT_0061)*(IT_0012*IT_0062 + IT_0043*IT_0065 + (
      -0.166666666666667)*IT_0061*IT_0066);
    return create_ccomplex_return(IT_0067);
}

