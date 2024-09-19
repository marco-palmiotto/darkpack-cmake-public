#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_ss_L_to_c_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_ss_L_to_c_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0002;
    const ccomplex_t IT_0004 = IT_0000*IT_0003;
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = pow(m_C_1, 2);
    const ccomplex_t IT_0007 = pow(m_ss_L, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + IT_0005 + -IT_0006 + -IT_0007 
      + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = s_14*s_34;
    const ccomplex_t IT_0012 = (-2)*IT_0000;
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = s_34*IT_0006;
    const ccomplex_t IT_0019 = 8*IT_0018;
    const ccomplex_t IT_0020 = s_13*s_14;
    const ccomplex_t IT_0021 = (-16)*IT_0020;
    const ccomplex_t IT_0022 = IT_0019 + IT_0021;
    const ccomplex_t IT_0023 = -IT_0000;
    const ccomplex_t IT_0024 = IT_0015*IT_0023;
    const ccomplex_t IT_0025 = (-8)*s_13;
    const ccomplex_t IT_0026 = s_14*IT_0025;
    const ccomplex_t IT_0027 = 2*IT_0010 + IT_0016;
    const ccomplex_t IT_0028 = s_34*IT_0025;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = sin(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0002*IT_0029*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0013*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0012*IT_0035;
    const ccomplex_t IT_0037 = m_c*m_C_1;
    const ccomplex_t IT_0038 = (-8)*s_14*IT_0037;
    const ccomplex_t IT_0039 = IT_0000*IT_0033;
    const ccomplex_t IT_0040 = IT_0008*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0036 + (-2)*IT_0041;
    const ccomplex_t IT_0043 = (-8)*s_34*IT_0037;
    const ccomplex_t IT_0044 = IT_0023*IT_0035;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = m_c*s_14*m_C_1;
    const ccomplex_t IT_0047 = m_c*s_34*m_C_1;
    const ccomplex_t IT_0048 = s_14*IT_0005;
    const ccomplex_t IT_0049 = (-8)*IT_0010;
    const ccomplex_t IT_0050 = -IT_0041;
    const ccomplex_t IT_0051 = IT_0006*IT_0025;
    const ccomplex_t IT_0052 = IT_0025*IT_0037;
    const ccomplex_t IT_0053 = pow(s_13, 2);
    const ccomplex_t IT_0054 = (-8)*IT_0053;
    const ccomplex_t IT_0055 = (-8)*IT_0006*IT_0037;
    const ccomplex_t IT_0056 = (-8)*conj(IT_0010);
    const ccomplex_t IT_0057 = (-8)*IT_0046;
    const ccomplex_t IT_0058 = (-8)*IT_0005;
    const ccomplex_t IT_0059 = IT_0037*IT_0058;
    const ccomplex_t IT_0060 = (-8)*IT_0048;
    const ccomplex_t IT_0061 = (-8)*IT_0047;
    const ccomplex_t IT_0062 = IT_0010*(32*conj(IT_0010)*IT_0011 + conj
      (IT_0017)*IT_0022) + IT_0024*(conj(IT_0017)*IT_0026 + conj(IT_0027)
      *IT_0028 + conj(IT_0036)*IT_0038 + conj(IT_0042)*IT_0043) + (IT_0026*conj
      (IT_0036) + conj(IT_0017)*IT_0038 + IT_0028*conj(IT_0042) + conj(IT_0027)
      *IT_0043)*IT_0045 + (conj(IT_0036)*IT_0046 + conj(IT_0042)*IT_0047 + conj
      (IT_0027)*IT_0048)*IT_0049 + IT_0036*(conj(IT_0024)*IT_0038 + IT_0026*conj
      (IT_0045) + IT_0022*conj(IT_0050) + conj(IT_0036)*IT_0051 + conj(IT_0027)
      *IT_0052 + conj(IT_0042)*IT_0054 + conj(IT_0017)*IT_0055 + IT_0046*IT_0056
      ) + IT_0017*(conj(IT_0010)*IT_0022 + conj(IT_0024)*IT_0026 + IT_0038*conj
      (IT_0045) + conj(IT_0017)*IT_0051 + conj(IT_0042)*IT_0052 + conj(IT_0027)
      *IT_0054 + conj(IT_0036)*IT_0055 + conj(IT_0050)*IT_0057) + IT_0042*(conj
      (IT_0024)*IT_0043 + IT_0028*conj(IT_0045) + conj(IT_0017)*IT_0052 + conj
      (IT_0036)*IT_0054 + IT_0047*IT_0056 + s_13*conj(IT_0042)*IT_0058 + conj
      (IT_0027)*IT_0059 + conj(IT_0050)*IT_0060) + IT_0050*(IT_0022*conj(IT_0036
      ) + 32*IT_0011*conj(IT_0050) + conj(IT_0017)*IT_0057 + conj(IT_0042)
      *IT_0060 + conj(IT_0027)*IT_0061) + IT_0027*(conj(IT_0024)*IT_0028 +
       IT_0043*conj(IT_0045) + conj(IT_0036)*IT_0052 + conj(IT_0017)*IT_0054 +
       IT_0048*IT_0056 + s_13*conj(IT_0027)*IT_0058 + conj(IT_0042)*IT_0059 +
       conj(IT_0050)*IT_0061);
    return create_ccomplex_return(IT_0062);
}

