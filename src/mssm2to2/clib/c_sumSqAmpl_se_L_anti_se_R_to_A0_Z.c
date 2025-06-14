#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_R_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_R_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = cpow(IT_0004, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0005
      *IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = pow(m_se_L, 2);
    const ccomplex_t IT_0014 = pow(m_se_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + -IT_0013 + -IT_0014 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0008*IT_0017;
    const ccomplex_t IT_0019 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0002
      *IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0013 + -IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, -2);
    const ccomplex_t IT_0029 = pow(s_14, 2);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0006 + 
      -IT_0000*IT_0007);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = m_e*e_em*mu_h*IT_0007*IT_0009;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = pow(m_A0, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0014 + IT_0035 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0032*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = s_14*s_34*IT_0028;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_13 + IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0045 = (-2)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0013 + IT_0035 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0034*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0045*IT_0048;
    const ccomplex_t IT_0050 = -IT_0018 + -IT_0026 + IT_0039 + IT_0049;
    const ccomplex_t IT_0051 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = s_12 + IT_0052;
    const ccomplex_t IT_0054 = -IT_0027;
    const ccomplex_t IT_0055 = -conj(IT_0027);
    const ccomplex_t IT_0056 = pow(s_34, 2);
    const ccomplex_t IT_0057 = IT_0028*IT_0056;
    const ccomplex_t IT_0058 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = s_23 + IT_0059;
    const ccomplex_t IT_0061 = pow(s_24, 2);
    const ccomplex_t IT_0062 = IT_0028*IT_0061;
    const ccomplex_t IT_0063 = -IT_0050;
    const ccomplex_t IT_0064 = (conj(IT_0027)*(IT_0013 + -IT_0030) + conj
      (IT_0040)*IT_0043 + conj(IT_0050)*IT_0053)*IT_0054 + IT_0050*IT_0053
      *IT_0055 + IT_0040*(IT_0043*IT_0055 + -conj(IT_0040)*(IT_0035 + -IT_0057) 
      + -conj(IT_0050)*IT_0060) + (conj(IT_0040)*IT_0060 + conj(IT_0050)*
      (IT_0014 + -IT_0062))*IT_0063;
    return create_ccomplex_return(IT_0064);
}

