#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_R_to_anti_d_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_R_to_anti_d_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = pow(m_d, 2);
    const ccomplex_t IT_0010 = pow(m_su_R, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-0.333333333333333)*IT_0000;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_C_1, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + IT_0009 + -IT_0010 + -IT_0016 
      + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = s_14*s_34;
    const ccomplex_t IT_0021 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0022 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0009 + IT_0016 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (IT_0021 + -IT_0022)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0022*IT_0025;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-6)*s_13;
    const ccomplex_t IT_0031 = conj(IT_0029)*IT_0030;
    const ccomplex_t IT_0032 = (-2)*IT_0013;
    const ccomplex_t IT_0033 = 6*s_13;
    const ccomplex_t IT_0034 = conj(IT_0032)*IT_0033;
    const ccomplex_t IT_0035 = IT_0031 + IT_0034;
    const ccomplex_t IT_0036 = 2*IT_0019 + IT_0028;
    const ccomplex_t IT_0037 = s_34*IT_0030;
    const ccomplex_t IT_0038 = s_34*IT_0016;
    const ccomplex_t IT_0039 = s_13*s_14;
    const ccomplex_t IT_0040 = (-2)*IT_0039;
    const ccomplex_t IT_0041 = IT_0038 + IT_0040;
    const ccomplex_t IT_0042 = pow(s_13, 2);
    const ccomplex_t IT_0043 = (-6)*IT_0042;
    const ccomplex_t IT_0044 = s_14*IT_0009;
    const ccomplex_t IT_0045 = s_13*s_34;
    const ccomplex_t IT_0046 = (-2)*IT_0045;
    const ccomplex_t IT_0047 = IT_0044 + IT_0046;
    const ccomplex_t IT_0048 = (-6)*IT_0013;
    const ccomplex_t IT_0049 = (-6)*IT_0019;
    const ccomplex_t IT_0050 = 6*IT_0042;
    const ccomplex_t IT_0051 = (-6)*conj(IT_0013);
    const ccomplex_t IT_0052 = (-6)*conj(IT_0019);
    const ccomplex_t IT_0053 = 0.166666666666667*conj(IT_0032);
    const ccomplex_t IT_0054 = (IT_0013 + IT_0019)*(24*conj(IT_0013) + 24*conj
      (IT_0019))*IT_0020 + IT_0027*(s_14*IT_0035 + conj(IT_0036)*IT_0037) + 6
      *IT_0029*(0.166666666666667*s_14*conj(IT_0027)*IT_0030 + 0.166666666666667
      *IT_0016*IT_0035 + conj(IT_0013)*IT_0038 + conj(IT_0019)*IT_0041 +
       0.166666666666667*conj(IT_0036)*IT_0043) + conj(IT_0036)*(IT_0047*IT_0048
       + IT_0044*IT_0049) + IT_0032*(s_14*conj(IT_0027)*IT_0033 + IT_0016*
      (IT_0030*conj(IT_0032) + conj(IT_0029)*IT_0033) + conj(IT_0036)*IT_0050 +
       IT_0038*IT_0051 + IT_0041*IT_0052) + (-6)*IT_0036*(s_13*IT_0009*conj
      (IT_0036) + (-0.166666666666667)*conj(IT_0027)*IT_0037 + (
      -0.166666666666667)*conj(IT_0029)*IT_0043 + (-0.166666666666667)*conj
      (IT_0032)*IT_0050 + (-0.166666666666667)*IT_0047*IT_0051 + (
      -0.166666666666667)*IT_0044*IT_0052) + 6*IT_0038*(IT_0013*conj(IT_0029) +
       IT_0048*IT_0053) + 6*IT_0041*(IT_0019*conj(IT_0029) + IT_0049*IT_0053);
    return create_ccomplex_return(IT_0054);
}

