#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_anti_Hp_Hp(
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
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_mu, 2);
    const ccomplex_t IT_0009 = cpow(IT_0003, -2);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + IT_0006*IT_0007*IT_0008*IT_0009*IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0015 = 2*IT_0014;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = s_13 + s_23;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = pow(m_smu_R, 2);
    const ccomplex_t IT_0020 = s_12 + IT_0019;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0015, 2);
    const ccomplex_t IT_0024 = s_23*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0014, 2);
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0024 + (-2)*IT_0026;
    const ccomplex_t IT_0028 = IT_0022 + IT_0027;
    const ccomplex_t IT_0029 = cpow(s_12 + IT_0019 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0001, -1);
    const ccomplex_t IT_0033 = cpow(IT_0010, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*(IT_0010*IT_0032 + 
      -IT_0001*IT_0033);
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0032;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0020*IT_0038;
    const ccomplex_t IT_0040 = (-2)*IT_0036;
    const ccomplex_t IT_0041 = IT_0035*IT_0040;
    const ccomplex_t IT_0042 = IT_0020*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0040;
    const ccomplex_t IT_0044 = -s_23;
    const ccomplex_t IT_0045 = IT_0020 + IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0034*IT_0037;
    const ccomplex_t IT_0048 = (-0.5)*s_13 + (-0.5)*s_23;
    const ccomplex_t IT_0049 = IT_0020 + IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (-0.5)*IT_0042 + 0.5*IT_0046 + -IT_0050;
    const ccomplex_t IT_0052 = IT_0039 + IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cos(alpha);
    const ccomplex_t IT_0057 = sin(alpha);
    const ccomplex_t IT_0058 = pow(m_W, -1);
    const ccomplex_t IT_0059 = cpow(IT_0003, -1);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0010*
      (IT_0005*IT_0056 + IT_0003*IT_0057) + -IT_0008*IT_0033*IT_0057*IT_0058
      *IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0002*IT_0010;
    const ccomplex_t IT_0063 = -IT_0033 + IT_0062;
    const ccomplex_t IT_0064 = cpow(IT_0005, 3);
    const ccomplex_t IT_0065 = 3*IT_0033 + IT_0062;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_W*e_em*(IT_0056*IT_0063
      *IT_0064 + IT_0003*(IT_0006*IT_0057*IT_0065 + -IT_0003*(IT_0003*IT_0057
      *IT_0063 + IT_0005*IT_0056*IT_0065)));
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0061*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0010*
      (IT_0003*IT_0056 + -IT_0005*IT_0057) + -IT_0008*IT_0033*IT_0056*IT_0058
      *IT_0059);
    const ccomplex_t IT_0073 = -IT_0033;
    const ccomplex_t IT_0074 = IT_0062 + IT_0073;
    const ccomplex_t IT_0075 = -IT_0003;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_W*e_em*(IT_0057*IT_0064
      *IT_0074 + IT_0075*(IT_0006*IT_0056*IT_0065 + -(IT_0005*IT_0057*IT_0065 + 
      -IT_0003*IT_0056*IT_0074)*IT_0075));
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0072*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + (-2)*IT_0019 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0033*IT_0058;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0083, 2);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = pow(m_Hp, 2);
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0019 + IT_0086 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0013 + (-0.5)*IT_0031 + (-2)*IT_0055 + 
      -IT_0071 + -IT_0081 + IT_0088;
    return create_ccomplex_return(IT_0089*conj(IT_0089));
}

