#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_1_to_t_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_1_to_t_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0001 = sin(alpha);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_t*mu_h*IT_0002;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = -IT_0001*IT_0005 + IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = pow(m_t, 2);
    const ccomplex_t IT_0015 = (-2)*U_st_00*((0.5*m_W*IT_0011*(IT_0008 + -1./3
      *IT_0007*IT_0013) + IT_0001*IT_0002*IT_0006*IT_0008*IT_0014)*conj(U_st_00)
       + -1./2*IT_0003*IT_0004*IT_0006*IT_0008*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0004*IT_0006*IT_0008*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0007
      *IT_0011*IT_0013 + 1.5*IT_0001*IT_0002*IT_0006*IT_0008*IT_0014)*conj
      (U_st_10));
    const ccomplex_t IT_0016 = (0 + _Complex_I*-1)*(e_em*IT_0015 + (
      -1.4142135623731)*A_t*IT_0001*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0000*IT_0017;
    const ccomplex_t IT_0019 = pow(m_sG, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0014 + IT_0019 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0024 = (-2)*U_st_00*((0.5*m_W*IT_0011*(IT_0008 + -1./3
      *IT_0007*IT_0013) + IT_0001*IT_0002*IT_0006*IT_0008*IT_0014)*conj(U_st_01)
       + -1./2*IT_0003*IT_0004*IT_0006*IT_0008*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0004*IT_0006*IT_0008*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0007
      *IT_0011*IT_0013 + 1.5*IT_0001*IT_0002*IT_0006*IT_0008*IT_0014)*conj
      (U_st_11));
    const ccomplex_t IT_0025 = (0 + _Complex_I*-1)*(e_em*IT_0024 + (
      -1.4142135623731)*A_t*IT_0001*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0014 + IT_0019 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_t*e_em*IT_0001*IT_0002
      *IT_0006*IT_0008;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0000*IT_0032;
    const ccomplex_t IT_0034 = pow(m_st_1, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + -IT_0019 + -IT_0034 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = m_t*IT_0037;
    const ccomplex_t IT_0039 = IT_0022 + IT_0030 + (-2)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = IT_0020*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0045 = IT_0026*IT_0044;
    const ccomplex_t IT_0046 = IT_0028*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0032*IT_0040;
    const ccomplex_t IT_0049 = IT_0035*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = m_t*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0047 + 2*IT_0051;
    const ccomplex_t IT_0053 = 8*s_13;
    const ccomplex_t IT_0054 = m_t*m_sG;
    const ccomplex_t IT_0055 = 8*IT_0039;
    const ccomplex_t IT_0056 = -IT_0037;
    const ccomplex_t IT_0057 = m_t*s_14;
    const ccomplex_t IT_0058 = m_sG*s_34;
    const ccomplex_t IT_0059 = 8*conj(IT_0039);
    const ccomplex_t IT_0060 = 8*conj(IT_0052);
    const ccomplex_t IT_0061 = pow(m_H0, 2);
    const ccomplex_t IT_0062 = m_t*m_sG*IT_0061;
    const ccomplex_t IT_0063 = s_13*IT_0061;
    const ccomplex_t IT_0064 = s_14*s_34;
    const ccomplex_t IT_0065 = (-2)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = 8*IT_0052;
    const ccomplex_t IT_0068 = 0.125*IT_0058;
    const ccomplex_t IT_0069 = (IT_0039*conj(IT_0039) + IT_0052*conj(IT_0052))
      *IT_0053 + conj(IT_0052)*IT_0054*IT_0055 + (IT_0052*IT_0054 + IT_0056
      *IT_0057 + IT_0050*IT_0058)*IT_0059 + (IT_0050*IT_0057 + IT_0056*IT_0058)
      *IT_0060 + 8*conj(IT_0050)*(IT_0056*IT_0062 + -IT_0050*IT_0066 + 0.125
      *IT_0057*IT_0067 + IT_0055*IT_0068) + 8*conj(IT_0056)*(0.125*IT_0055
      *IT_0057 + IT_0050*IT_0062 + -IT_0056*IT_0066 + IT_0067*IT_0068);
    return create_ccomplex_return(IT_0069);
}

