#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sd_L_to_u_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sd_L_to_u_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_u, 2);
    const ccomplex_t IT_0006 = pow(m_C_2, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = IT_0000*IT_0004;
    const ccomplex_t IT_0012 = pow(m_sd_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0005 + -IT_0006 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0010 + 2*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0003*IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = IT_0000*IT_0021;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0024 + (-2)*IT_0027;
    const ccomplex_t IT_0029 = (-8)*IT_0005;
    const ccomplex_t IT_0030 = -IT_0000;
    const ccomplex_t IT_0031 = IT_0009*IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0030;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (-8)*s_13;
    const ccomplex_t IT_0035 = s_34*IT_0034;
    const ccomplex_t IT_0036 = -IT_0010;
    const ccomplex_t IT_0037 = -IT_0027;
    const ccomplex_t IT_0038 = m_u*s_14*m_C_2;
    const ccomplex_t IT_0039 = s_34*IT_0006;
    const ccomplex_t IT_0040 = 8*IT_0039;
    const ccomplex_t IT_0041 = s_13*s_14;
    const ccomplex_t IT_0042 = (-16)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040 + IT_0042;
    const ccomplex_t IT_0044 = s_14*IT_0034;
    const ccomplex_t IT_0045 = IT_0006*IT_0034;
    const ccomplex_t IT_0046 = pow(s_13, 2);
    const ccomplex_t IT_0047 = (-8)*IT_0046;
    const ccomplex_t IT_0048 = m_u*m_C_2;
    const ccomplex_t IT_0049 = (-8)*s_14;
    const ccomplex_t IT_0050 = (-8)*s_34;
    const ccomplex_t IT_0051 = (-8)*IT_0006;
    const ccomplex_t IT_0052 = (conj(IT_0016)*IT_0028 + IT_0016*conj(IT_0028))
      *IT_0029 + IT_0034*(conj(IT_0016)*IT_0024 + IT_0016*conj(IT_0024) + conj
      (IT_0028)*IT_0036 + IT_0028*conj(IT_0036)) + (conj(IT_0024)*IT_0031 +
       IT_0024*conj(IT_0031) + conj(IT_0033)*IT_0036 + IT_0033*conj(IT_0036))
      *IT_0049 + (conj(IT_0028)*IT_0031 + IT_0028*conj(IT_0031) + conj(IT_0016)
      *IT_0033 + IT_0016*conj(IT_0033))*IT_0050 + (conj(IT_0024)*IT_0036 +
       IT_0024*conj(IT_0036))*IT_0051;
    const ccomplex_t IT_0053 = m_u*s_34*m_C_2;
    const ccomplex_t IT_0054 = (-8)*conj(IT_0016);
    const ccomplex_t IT_0055 = s_14*IT_0005;
    const ccomplex_t IT_0056 = (-8)*conj(IT_0028);
    const ccomplex_t IT_0057 = s_14*s_34;
    const ccomplex_t IT_0058 = (-8)*IT_0028;
    const ccomplex_t IT_0059 = (-0.125)*IT_0055;
    const ccomplex_t IT_0060 = (-8)*IT_0016;
    const ccomplex_t IT_0061 = (-0.125)*IT_0053;
    const ccomplex_t IT_0062 = s_13*(IT_0016*conj(IT_0016) + IT_0028*conj
      (IT_0028))*IT_0029 + (conj(IT_0016)*IT_0031 + IT_0016*conj(IT_0031) + conj
      (IT_0028)*IT_0033 + IT_0028*conj(IT_0033))*IT_0035 + (-8)*conj(IT_0036)*
      (IT_0037*IT_0038 + (-0.125)*IT_0015*IT_0043 + (-0.125)*IT_0031*IT_0044 + (
      -0.125)*IT_0036*IT_0045 + (-0.125)*IT_0016*IT_0047) + IT_0036*(conj
      (IT_0031)*IT_0044 + conj(IT_0016)*IT_0047) + (-8)*conj(IT_0024)*(IT_0015
      *IT_0038 + (-0.125)*IT_0037*IT_0043 + (-0.125)*IT_0033*IT_0044 + (-0.125)
      *IT_0024*IT_0045 + (-0.125)*IT_0028*IT_0047) + IT_0024*(conj(IT_0033)
      *IT_0044 + conj(IT_0028)*IT_0047) + IT_0048*IT_0052 + IT_0015*(IT_0054
      *IT_0055 + IT_0053*IT_0056) + IT_0037*(IT_0053*IT_0054 + IT_0055*IT_0056) 
      + (-8)*conj(IT_0015)*(IT_0024*IT_0038 + (-0.125)*IT_0036*IT_0043 + (-4)
      *IT_0015*IT_0057 + IT_0059*IT_0060 + IT_0058*IT_0061) + (-8)*conj(IT_0037)
      *(IT_0036*IT_0038 + (-0.125)*IT_0024*IT_0043 + (-4)*IT_0037*IT_0057 +
       IT_0058*IT_0059 + IT_0060*IT_0061);
    return create_ccomplex_return(IT_0062);
}

