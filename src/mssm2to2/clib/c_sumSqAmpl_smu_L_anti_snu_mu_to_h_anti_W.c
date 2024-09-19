#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_snu_mu_to_h_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_snu_mu_to_h_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = pow(m_mu, 2);
    const ccomplex_t IT_0010 = cpow(IT_0002, 2);
    const ccomplex_t IT_0011 = cpow(IT_0001, 2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005*
      (IT_0008*IT_0009 + -m_W*(IT_0010 + -IT_0011));
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, 2);
    const ccomplex_t IT_0015 = pow(m_smu_L, 2);
    const ccomplex_t IT_0016 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0016 
      + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0002, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0005*(m_W*IT_0002 + (-0.5)*IT_0008*IT_0009*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0015 + -IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0031*IT_0035;
    const ccomplex_t IT_0037 = -IT_0020 + -IT_0029 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -2);
    const ccomplex_t IT_0039 = pow(s_14, 2);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = -IT_0030;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -2);
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = IT_0005 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = pow(m_h, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0048 + reg_prop + (0 +
       _Complex_I*1)*m_snu_mu*Gamma_numl, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0041*IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = s_14*s_34*IT_0038;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = s_13 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0005
      *IT_0008*IT_0009*IT_0023 + (-0.5)*m_W*(IT_0000*IT_0001 + IT_0002*IT_0003)*
      (IT_0005 + -IT_0004*IT_0043));
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0048 + (0 + _Complex_I*1)
      *m_smu_L*Gamma_ml + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0041*IT_0060;
    const ccomplex_t IT_0062 = -IT_0041;
    const ccomplex_t IT_0063 = IT_0031 + IT_0062;
    const ccomplex_t IT_0064 = IT_0035*IT_0063;
    const ccomplex_t IT_0065 = -IT_0020 + -IT_0029 + IT_0052 + IT_0061 + 
      -IT_0064;
    const ccomplex_t IT_0066 = s_14*s_24*IT_0038;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = s_12 + IT_0067;
    const ccomplex_t IT_0069 = -IT_0037;
    const ccomplex_t IT_0070 = -conj(IT_0037);
    const ccomplex_t IT_0071 = pow(s_34, 2);
    const ccomplex_t IT_0072 = IT_0038*IT_0071;
    const ccomplex_t IT_0073 = s_24*s_34*IT_0038;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = s_23 + IT_0074;
    const ccomplex_t IT_0076 = pow(s_24, 2);
    const ccomplex_t IT_0077 = IT_0038*IT_0076;
    const ccomplex_t IT_0078 = -IT_0065;
    const ccomplex_t IT_0079 = (conj(IT_0037)*(IT_0015 + -IT_0040) + conj
      (IT_0053)*IT_0056 + conj(IT_0065)*IT_0068)*IT_0069 + IT_0065*IT_0068
      *IT_0070 + IT_0053*(IT_0056*IT_0070 + -conj(IT_0053)*(IT_0048 + -IT_0072) 
      + -conj(IT_0065)*IT_0075) + (conj(IT_0053)*IT_0075 + conj(IT_0065)*
      (IT_0016 + -IT_0077))*IT_0078;
    return create_ccomplex_return(IT_0079);
}

