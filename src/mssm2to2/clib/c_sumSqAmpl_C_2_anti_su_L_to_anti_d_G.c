#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_d_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_d_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_d, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = pow(m_su_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0006;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = (-2)*IT_0004;
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = 2*IT_0011 + IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = cos(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0006*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0013*IT_0022;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0012*IT_0027;
    const ccomplex_t IT_0029 = (-2)*IT_0025 + -IT_0028;
    const ccomplex_t IT_0030 = (-8)*IT_0000;
    const ccomplex_t IT_0031 = -IT_0004;
    const ccomplex_t IT_0032 = -IT_0015*(IT_0012 + -IT_0031);
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = -IT_0027*(IT_0012 + -IT_0031);
    const ccomplex_t IT_0035 = (-8)*s_13;
    const ccomplex_t IT_0036 = s_34*IT_0035;
    const ccomplex_t IT_0037 = -IT_0016;
    const ccomplex_t IT_0038 = -IT_0025;
    const ccomplex_t IT_0039 = m_d*s_14*m_C_2;
    const ccomplex_t IT_0040 = s_34*IT_0001;
    const ccomplex_t IT_0041 = 8*IT_0040;
    const ccomplex_t IT_0042 = s_13*s_14;
    const ccomplex_t IT_0043 = (-16)*IT_0042;
    const ccomplex_t IT_0044 = IT_0041 + IT_0043;
    const ccomplex_t IT_0045 = s_14*IT_0035;
    const ccomplex_t IT_0046 = IT_0001*IT_0035;
    const ccomplex_t IT_0047 = pow(s_13, 2);
    const ccomplex_t IT_0048 = (-8)*IT_0047;
    const ccomplex_t IT_0049 = m_d*m_C_2;
    const ccomplex_t IT_0050 = (-8)*s_14;
    const ccomplex_t IT_0051 = (-8)*s_34;
    const ccomplex_t IT_0052 = (-8)*IT_0001;
    const ccomplex_t IT_0053 = (conj(IT_0017)*IT_0029 + IT_0017*conj(IT_0029))
      *IT_0030 + IT_0035*(conj(IT_0017)*IT_0028 + IT_0017*conj(IT_0028) + conj
      (IT_0029)*IT_0037 + IT_0029*conj(IT_0037)) + (conj(IT_0028)*IT_0033 +
       IT_0028*conj(IT_0033) + conj(IT_0034)*IT_0037 + IT_0034*conj(IT_0037))
      *IT_0050 + (conj(IT_0029)*IT_0033 + IT_0029*conj(IT_0033) + conj(IT_0017)
      *IT_0034 + IT_0017*conj(IT_0034))*IT_0051 + (conj(IT_0028)*IT_0037 +
       IT_0028*conj(IT_0037))*IT_0052;
    const ccomplex_t IT_0054 = m_d*s_34*m_C_2;
    const ccomplex_t IT_0055 = (-8)*conj(IT_0017);
    const ccomplex_t IT_0056 = s_14*IT_0000;
    const ccomplex_t IT_0057 = (-8)*conj(IT_0029);
    const ccomplex_t IT_0058 = (-0.125)*IT_0028;
    const ccomplex_t IT_0059 = s_14*s_34;
    const ccomplex_t IT_0060 = (-8)*IT_0029;
    const ccomplex_t IT_0061 = (-0.125)*IT_0056;
    const ccomplex_t IT_0062 = (-8)*IT_0017;
    const ccomplex_t IT_0063 = (-0.125)*IT_0054;
    const ccomplex_t IT_0064 = s_13*(IT_0017*conj(IT_0017) + IT_0029*conj
      (IT_0029))*IT_0030 + (conj(IT_0017)*IT_0033 + IT_0017*conj(IT_0033) + conj
      (IT_0029)*IT_0034 + IT_0029*conj(IT_0034))*IT_0036 + (-8)*conj(IT_0037)*
      (IT_0038*IT_0039 + (-0.125)*IT_0011*IT_0044 + (-0.125)*IT_0033*IT_0045 + (
      -0.125)*IT_0037*IT_0046 + (-0.125)*IT_0017*IT_0048) + IT_0037*(conj
      (IT_0033)*IT_0045 + conj(IT_0017)*IT_0048) + IT_0028*(conj(IT_0034)
      *IT_0045 + conj(IT_0029)*IT_0048) + IT_0049*IT_0053 + IT_0011*(IT_0055
      *IT_0056 + IT_0054*IT_0057) + IT_0038*(IT_0054*IT_0055 + IT_0056*IT_0057) 
      + (-8)*conj(IT_0028)*(IT_0011*IT_0039 + (-0.125)*IT_0038*IT_0044 + (-0.125
      )*IT_0034*IT_0045 + (-0.125)*IT_0029*IT_0048 + IT_0046*IT_0058) + (-8)
      *conj(IT_0011)*(IT_0028*IT_0039 + (-0.125)*IT_0037*IT_0044 + (-4)*IT_0011
      *IT_0059 + IT_0061*IT_0062 + IT_0060*IT_0063) + (-8)*conj(IT_0038)*
      (IT_0037*IT_0039 + IT_0044*IT_0058 + (-4)*IT_0038*IT_0059 + IT_0060
      *IT_0061 + IT_0062*IT_0063);
    return create_ccomplex_return(IT_0064);
}

