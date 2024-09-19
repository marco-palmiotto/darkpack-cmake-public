#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = IT_0005 + IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = sin(alpha);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0001
      *IT_0003 + -IT_0012*IT_0013);
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0013;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0003
      *IT_0012 + IT_0001*IT_0013);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = cpow(IT_0008, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0004*IT_0031;
    const ccomplex_t IT_0035 = e_em*IT_0034;
    const ccomplex_t IT_0036 = IT_0005*IT_0008;
    const ccomplex_t IT_0037 = e_em*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0035 + 3*IT_0037);
    const ccomplex_t IT_0039 = (-0.166666666666667)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0033*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0031;
    const ccomplex_t IT_0046 = 0.333333333333333*IT_0045;
    const ccomplex_t IT_0047 = IT_0040*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0033*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -s_14;
    const ccomplex_t IT_0052 = s_13 + IT_0051;
    const ccomplex_t IT_0053 = m_b*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0032;
    const ccomplex_t IT_0055 = m_b*IT_0054;
    const ccomplex_t IT_0056 = m_b*IT_0033;
    const ccomplex_t IT_0057 = -(IT_0042 + -IT_0048)*(IT_0055 + -IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0030;
    const ccomplex_t IT_0061 = s_34 + IT_0060;
    const ccomplex_t IT_0062 = m_b*s_13;
    const ccomplex_t IT_0063 = IT_0017*IT_0030;
    const ccomplex_t IT_0064 = m_b*s_14;
    const ccomplex_t IT_0065 = s_34*IT_0017;
    const ccomplex_t IT_0066 = (-6)*IT_0065;
    const ccomplex_t IT_0067 = s_13*s_14;
    const ccomplex_t IT_0068 = 12*IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = conj(IT_0029)*IT_0053;
    const ccomplex_t IT_0071 = (-6)*conj(IT_0044);
    const ccomplex_t IT_0072 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0073 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0074 = 6*s_34;
    const ccomplex_t IT_0075 = conj(IT_0029)*IT_0061;
    const ccomplex_t IT_0076 = IT_0029*(conj(IT_0029)*(12*s_34 + (-12)*IT_0030
      ) + 6*(conj(IT_0044) + conj(IT_0050))*IT_0053 + 6*(conj(IT_0058) + conj
      (IT_0059))*IT_0061) + 6*IT_0044*(conj(IT_0059)*IT_0062 + -conj(IT_0050)
      *IT_0063 + -conj(IT_0058)*IT_0064 + 0.166666666666667*conj(IT_0044)
      *IT_0069 + IT_0070) + 6*IT_0050*(conj(IT_0058)*IT_0062 + 0.166666666666667
      *conj(IT_0050)*IT_0069 + IT_0070 + 0.166666666666667*IT_0063*IT_0071 +
       IT_0072*IT_0073) + 6*IT_0058*(conj(IT_0050)*IT_0062 + 0.166666666666667
      *IT_0030*IT_0072 + IT_0071*IT_0073 + 0.166666666666667*conj(IT_0058)
      *IT_0074 + IT_0075) + (-6)*IT_0059*(IT_0030*conj(IT_0058) + -conj(IT_0044)
      *IT_0062 + conj(IT_0050)*IT_0064 + (-0.166666666666667)*conj(IT_0059)
      *IT_0074 + -IT_0075);
    return create_ccomplex_return(IT_0076);
}

