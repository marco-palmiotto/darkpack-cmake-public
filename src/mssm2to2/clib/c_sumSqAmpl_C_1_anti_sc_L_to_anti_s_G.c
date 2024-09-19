#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_sc_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0006;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = s_14*s_34;
    const ccomplex_t IT_0013 = (-2)*IT_0004;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = s_34*IT_0001;
    const ccomplex_t IT_0020 = 8*IT_0019;
    const ccomplex_t IT_0021 = s_13*s_14;
    const ccomplex_t IT_0022 = (-16)*IT_0021;
    const ccomplex_t IT_0023 = IT_0020 + IT_0022;
    const ccomplex_t IT_0024 = -IT_0004;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0006*IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0014*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -(IT_0013 + -IT_0024)*IT_0031;
    const ccomplex_t IT_0033 = IT_0013*IT_0031;
    const ccomplex_t IT_0034 = (-8)*s_13;
    const ccomplex_t IT_0035 = s_14*IT_0034;
    const ccomplex_t IT_0036 = IT_0004*IT_0029;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0033 + (-2)*IT_0038;
    const ccomplex_t IT_0040 = s_34*IT_0034;
    const ccomplex_t IT_0041 = m_s*m_C_1;
    const ccomplex_t IT_0042 = (-8)*s_14*IT_0041;
    const ccomplex_t IT_0043 = 2*IT_0011 + IT_0017;
    const ccomplex_t IT_0044 = (-8)*s_34*IT_0041;
    const ccomplex_t IT_0045 = -IT_0016*(IT_0013 + -IT_0024);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = s_14*IT_0000;
    const ccomplex_t IT_0048 = m_s*s_14*m_C_1;
    const ccomplex_t IT_0049 = m_s*s_34*m_C_1;
    const ccomplex_t IT_0050 = (-8)*IT_0011;
    const ccomplex_t IT_0051 = -IT_0038;
    const ccomplex_t IT_0052 = IT_0001*IT_0034;
    const ccomplex_t IT_0053 = IT_0034*IT_0041;
    const ccomplex_t IT_0054 = pow(s_13, 2);
    const ccomplex_t IT_0055 = (-8)*IT_0054;
    const ccomplex_t IT_0056 = (-8)*IT_0001*IT_0041;
    const ccomplex_t IT_0057 = (-8)*conj(IT_0011);
    const ccomplex_t IT_0058 = (-8)*IT_0048;
    const ccomplex_t IT_0059 = (-8)*IT_0000;
    const ccomplex_t IT_0060 = IT_0041*IT_0059;
    const ccomplex_t IT_0061 = (-8)*IT_0047;
    const ccomplex_t IT_0062 = (-8)*IT_0049;
    const ccomplex_t IT_0063 = IT_0011*(32*conj(IT_0011)*IT_0012 + conj
      (IT_0018)*IT_0023) + IT_0032*(conj(IT_0033)*IT_0035 + conj(IT_0039)
      *IT_0040 + conj(IT_0018)*IT_0042 + conj(IT_0043)*IT_0044) + (conj(IT_0018)
      *IT_0035 + conj(IT_0033)*IT_0042 + IT_0040*conj(IT_0043) + conj(IT_0039)
      *IT_0044)*IT_0046 + (conj(IT_0043)*IT_0047 + conj(IT_0033)*IT_0048 + conj
      (IT_0039)*IT_0049)*IT_0050 + IT_0033*(conj(IT_0032)*IT_0035 + IT_0042*conj
      (IT_0046) + IT_0023*conj(IT_0051) + conj(IT_0033)*IT_0052 + conj(IT_0043)
      *IT_0053 + conj(IT_0039)*IT_0055 + conj(IT_0018)*IT_0056 + IT_0048*IT_0057
      ) + IT_0018*(conj(IT_0011)*IT_0023 + conj(IT_0032)*IT_0042 + IT_0035*conj
      (IT_0046) + conj(IT_0018)*IT_0052 + conj(IT_0039)*IT_0053 + conj(IT_0043)
      *IT_0055 + conj(IT_0033)*IT_0056 + conj(IT_0051)*IT_0058) + IT_0039*(conj
      (IT_0032)*IT_0040 + IT_0044*conj(IT_0046) + conj(IT_0018)*IT_0053 + conj
      (IT_0033)*IT_0055 + IT_0049*IT_0057 + s_13*conj(IT_0039)*IT_0059 + conj
      (IT_0043)*IT_0060 + conj(IT_0051)*IT_0061) + IT_0051*(IT_0023*conj(IT_0033
      ) + 32*IT_0012*conj(IT_0051) + conj(IT_0018)*IT_0058 + conj(IT_0039)
      *IT_0061 + conj(IT_0043)*IT_0062) + IT_0043*(conj(IT_0032)*IT_0044 +
       IT_0040*conj(IT_0046) + conj(IT_0033)*IT_0053 + conj(IT_0018)*IT_0055 +
       IT_0047*IT_0057 + s_13*conj(IT_0043)*IT_0059 + conj(IT_0039)*IT_0060 +
       conj(IT_0051)*IT_0062);
    return create_ccomplex_return(IT_0063);
}

