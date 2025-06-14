#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sc_R_to_c_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sc_R_to_c_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_c = param->m_c;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = -g_s;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sc_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0012 +
       reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = pow(s_13, 2);
    const ccomplex_t IT_0018 = 1.33333333333333*IT_0017;
    const ccomplex_t IT_0019 = -IT_0000;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0023 = s_34*IT_0022;
    const ccomplex_t IT_0024 = IT_0008 + (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0025 = (-10.6666666666667)*IT_0017;
    const ccomplex_t IT_0026 = 1.33333333333333*s_13;
    const ccomplex_t IT_0027 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0028 = (-1.33333333333333)*IT_0017;
    const ccomplex_t IT_0029 = IT_0000*IT_0002;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (-2)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0032;
    const ccomplex_t IT_0035 = s_34*IT_0004;
    const ccomplex_t IT_0036 = s_13*s_14;
    const ccomplex_t IT_0037 = -IT_0016;
    const ccomplex_t IT_0038 = s_14*s_34;
    const ccomplex_t IT_0039 = (-2)*IT_0036;
    const ccomplex_t IT_0040 = IT_0035 + IT_0039;
    const ccomplex_t IT_0041 = s_14*IT_0003;
    const ccomplex_t IT_0042 = 1.33333333333333*conj(IT_0009);
    const ccomplex_t IT_0043 = (-1.33333333333333)*IT_0009;
    const ccomplex_t IT_0044 = IT_0016*((-10.6666666666667)*IT_0035 +
       21.3333333333333*IT_0036) + (10.6666666666667*IT_0035 + (
      -21.3333333333333)*IT_0036)*IT_0037 + 21.3333333333333*IT_0034*IT_0038 +
       1.33333333333333*((-0.5)*IT_0024 + 0.5*conj(IT_0024) + IT_0037)*IT_0040 +
       (-0.5)*IT_0041*(21.3333333333333*IT_0033 + IT_0042 + IT_0043);
    const ccomplex_t IT_0045 = (-0.5)*IT_0016;
    const ccomplex_t IT_0046 = s_13*s_34;
    const ccomplex_t IT_0047 = (-21.3333333333333)*IT_0046;
    const ccomplex_t IT_0048 = 10.6666666666667*IT_0041;
    const ccomplex_t IT_0049 = IT_0047 + IT_0048;
    const ccomplex_t IT_0050 = -IT_0009;
    const ccomplex_t IT_0051 = conj(IT_0009) + IT_0050;
    const ccomplex_t IT_0052 = 21.3333333333333*IT_0046;
    const ccomplex_t IT_0053 = (-10.6666666666667)*IT_0041;
    const ccomplex_t IT_0054 = IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = (-2)*IT_0033;
    const ccomplex_t IT_0056 = (-2)*IT_0046;
    const ccomplex_t IT_0057 = IT_0041 + IT_0056;
    const ccomplex_t IT_0058 = IT_0035*(24*IT_0016 + 12*IT_0024 + (-12)*conj
      (IT_0024) + (-48)*IT_0037) + IT_0038*((-96)*IT_0034 + (-96)*IT_0045) +
       IT_0049*IT_0051 + IT_0054*IT_0055 + 2.66666666666667*(IT_0033 + 0.375
      *IT_0042 + 0.375*IT_0043)*IT_0057;
    const ccomplex_t IT_0059 = (-24)*IT_0016;
    const ccomplex_t IT_0060 = IT_0037*(s_14*((-12)*IT_0021 + 12*conj(IT_0021)
      ) + s_13*((-12)*IT_0009 + 12*conj(IT_0009) + 24*IT_0033) + IT_0004*((-12)
      *IT_0024 + 12*conj(IT_0024) + 24*IT_0037 + IT_0059));
    const ccomplex_t IT_0061 = (-10.6666666666667)*IT_0003*IT_0009*conj
      (IT_0009) + IT_0004*((-10.6666666666667)*IT_0024*conj(IT_0024) + IT_0016*
      (10.6666666666667*IT_0016 + 1.33333333333333*IT_0024 + (-1.33333333333333)
      *conj(IT_0024))) + IT_0033*((-21.3333333333333)*s_13*IT_0016 + IT_0003*((
      -1.33333333333333)*IT_0009 + 1.33333333333333*conj(IT_0009) +
       10.6666666666667*IT_0033)) + IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0016;
    const ccomplex_t IT_0063 = IT_0051*IT_0054 + IT_0049*IT_0055 + (
      -2.66666666666667)*((-0.5)*IT_0009 + 0.5*conj(IT_0009) + IT_0033)*IT_0057 
      + IT_0035*((-12)*IT_0024 + 12*conj(IT_0024) + 48*IT_0037 + IT_0059) +
       IT_0038*(96*IT_0034 + 192*IT_0045 + (-96)*IT_0062);
    const ccomplex_t IT_0064 = IT_0009*(IT_0016*IT_0018 + conj(IT_0021)
      *IT_0023 + conj(IT_0024)*IT_0025) + s_14*(IT_0022*(conj(IT_0021)*IT_0024 +
       IT_0021*conj(IT_0024)) + IT_0016*(IT_0021*IT_0026 + conj(IT_0021)*IT_0027
      )) + conj(IT_0009)*(IT_0021*IT_0023 + IT_0024*IT_0025 + IT_0016*IT_0028) +
       (IT_0018*conj(IT_0024) + s_34*(conj(IT_0021)*IT_0026 + IT_0021*IT_0027) +
       IT_0024*IT_0028)*IT_0033 + (-2)*IT_0034*IT_0044 + IT_0045*IT_0058 + s_13
      *IT_0061 + IT_0062*IT_0063;
    return create_ccomplex_return(IT_0064);
}

