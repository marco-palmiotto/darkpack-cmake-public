#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_h_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(alpha);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = pow(m_W, -2);
    const ccomplex_t IT_0006 = pow(m_u, 2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0002*
      (IT_0004 + (-3)*(1 + -IT_0005*IT_0006*IT_0008)*IT_0010);
    const ccomplex_t IT_0012 = (-0.166666666666667)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = cpow(IT_0007, -1);
    const ccomplex_t IT_0015 = cpow(IT_0009, -1);
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0006
      *IT_0013*IT_0014*IT_0015 + 0.166666666666667*m_W*(IT_0004*IT_0009 + (-3)
      *IT_0015)*(IT_0002*IT_0007 + -IT_0001*IT_0016));
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0004*IT_0009;
    const ccomplex_t IT_0020 = IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0002, 3);
    const ccomplex_t IT_0022 = cpow(IT_0002, 2);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*IT_0020*(IT_0016
      *IT_0021 + 5*IT_0001*((-0.2)*IT_0001*(IT_0001*IT_0007 + 5*IT_0002*IT_0016)
       + IT_0007*IT_0022));
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = pow(m_su_L, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + (-2)*IT_0026 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0006
      *IT_0013*IT_0014*IT_0015 + 0.166666666666667*m_W*(IT_0004*IT_0009 + (-3)
      *IT_0015)*(IT_0001*IT_0007 + IT_0002*IT_0016));
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0001, 2);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_W*e_em*IT_0020*(IT_0007
      *IT_0021 + (-5)*IT_0001*(IT_0001*IT_0002*IT_0007 + 1*IT_0016*(IT_0022 + 
      -1./5*IT_0032)));
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0031*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0001
      *IT_0013*IT_0014*IT_0015;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0002
      *IT_0013*IT_0014*IT_0015;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = pow(m_h, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0026 + IT_0044 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0018*IT_0031;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0044 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0044 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0048*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0026 + IT_0044 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0043*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0012 + -IT_0029 + -IT_0038 + IT_0047 +
       IT_0051 + IT_0054 + IT_0057;
    return create_ccomplex_return(3*IT_0058*conj(IT_0058));
}

