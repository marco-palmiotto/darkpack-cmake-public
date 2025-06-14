#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_R_to_anti_d_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_R_to_anti_d_A(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0002 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_d, 2);
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (IT_0001 + -IT_0002)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0000*IT_0009;
    const ccomplex_t IT_0018 = pow(m_su_R, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0010 + IT_0018 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (-2)*IT_0021;
    const ccomplex_t IT_0023 = 6*s_13*s_14;
    const ccomplex_t IT_0024 = IT_0002*IT_0014;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (-6)*s_13;
    const ccomplex_t IT_0027 = (-0.333333333333333)*IT_0000;
    const ccomplex_t IT_0028 = IT_0009*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0024 + 2*IT_0031;
    const ccomplex_t IT_0033 = s_34*IT_0026;
    const ccomplex_t IT_0034 = pow(s_13, 2);
    const ccomplex_t IT_0035 = 6*IT_0034;
    const ccomplex_t IT_0036 = (-6)*IT_0034;
    const ccomplex_t IT_0037 = s_34*IT_0011;
    const ccomplex_t IT_0038 = s_13*s_14;
    const ccomplex_t IT_0039 = (-2)*IT_0038;
    const ccomplex_t IT_0040 = IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = 6*IT_0025;
    const ccomplex_t IT_0042 = (-6)*IT_0022;
    const ccomplex_t IT_0043 = IT_0041 + IT_0042;
    const ccomplex_t IT_0044 = s_14*IT_0010;
    const ccomplex_t IT_0045 = (-6)*IT_0032;
    const ccomplex_t IT_0046 = s_13*s_34;
    const ccomplex_t IT_0047 = (-2)*IT_0046;
    const ccomplex_t IT_0048 = IT_0044 + IT_0047;
    const ccomplex_t IT_0049 = 6*conj(IT_0025);
    const ccomplex_t IT_0050 = (-6)*conj(IT_0022);
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0053 = s_14*s_34;
    const ccomplex_t IT_0054 = 24*conj(IT_0031);
    const ccomplex_t IT_0055 = 24*conj(IT_0021) + IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = IT_0011*IT_0026;
    const ccomplex_t IT_0058 = s_13*IT_0011;
    const ccomplex_t IT_0059 = conj(IT_0016)*(IT_0022*IT_0023 + s_14*IT_0025
      *IT_0026 + IT_0032*IT_0033) + (-6)*conj(IT_0032)*(s_13*IT_0010*IT_0032 + (
      -0.166666666666667)*IT_0016*IT_0033 + (-0.166666666666667)*IT_0022*IT_0035
       + (-0.166666666666667)*IT_0025*IT_0036) + conj(IT_0031)*(IT_0040*IT_0043 
      + IT_0044*IT_0045) + conj(IT_0021)*(IT_0037*IT_0043 + IT_0045*IT_0048) +
       IT_0031*(IT_0040*IT_0051 + IT_0044*IT_0052 + IT_0056) + IT_0021*(IT_0037
      *IT_0051 + IT_0048*IT_0052 + IT_0056) + conj(IT_0025)*(s_14*IT_0016
      *IT_0026 + IT_0032*IT_0036 + IT_0025*IT_0057 + 6*IT_0022*IT_0058) + conj
      (IT_0022)*(IT_0016*IT_0023 + IT_0032*IT_0035 + IT_0022*IT_0057 + 6*IT_0025
      *IT_0058);
    return create_ccomplex_return(IT_0059);
}

