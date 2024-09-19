#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sc_L_to_c_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sc_L_to_c_G(
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = -g_s;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = pow(m_sc_L, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0011 +
       reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 12*s_13*s_14;
    const ccomplex_t IT_0017 = -IT_0015;
    const ccomplex_t IT_0018 = (-12)*s_13*s_14;
    const ccomplex_t IT_0019 = (-2)*IT_0000;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0023 = s_34*IT_0022;
    const ccomplex_t IT_0024 = pow(s_13, 2);
    const ccomplex_t IT_0025 = (-12)*IT_0024;
    const ccomplex_t IT_0026 = (-10.6666666666667)*IT_0024;
    const ccomplex_t IT_0027 = IT_0000*IT_0002;
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0011 
      + -reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = s_13*s_14;
    const ccomplex_t IT_0032 = s_34*IT_0004;
    const ccomplex_t IT_0033 = s_14*s_34;
    const ccomplex_t IT_0034 = (-0.5)*IT_0032;
    const ccomplex_t IT_0035 = IT_0031 + IT_0034;
    const ccomplex_t IT_0036 = s_14*IT_0003;
    const ccomplex_t IT_0037 = 2*IT_0030;
    const ccomplex_t IT_0038 = 1.33333333333333*conj(IT_0020);
    const ccomplex_t IT_0039 = (-1.33333333333333)*IT_0020;
    const ccomplex_t IT_0040 = IT_0015*((-21.3333333333333)*IT_0031 +
       10.6666666666667*IT_0032) + IT_0017*(21.3333333333333*IT_0031 + (
      -10.6666666666667)*IT_0032) + 21.3333333333333*IT_0030*IT_0033 + (
      -2.66666666666667)*(IT_0015 + -IT_0017 + (-0.5)*IT_0021 + 0.5*conj(IT_0021
      ))*IT_0035 + (-0.5)*IT_0036*(21.3333333333333*IT_0037 + IT_0038 + IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0015;
    const ccomplex_t IT_0042 = s_13*s_34;
    const ccomplex_t IT_0043 = (-2)*IT_0042;
    const ccomplex_t IT_0044 = IT_0036 + IT_0043;
    const ccomplex_t IT_0045 = IT_0020 + -conj(IT_0020);
    const ccomplex_t IT_0046 = 21.3333333333333*IT_0042;
    const ccomplex_t IT_0047 = (-10.6666666666667)*IT_0036;
    const ccomplex_t IT_0048 = IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = conj(IT_0020) + IT_0021;
    const ccomplex_t IT_0050 = (-21.3333333333333)*IT_0042;
    const ccomplex_t IT_0051 = 10.6666666666667*IT_0036;
    const ccomplex_t IT_0052 = IT_0050 + IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0037;
    const ccomplex_t IT_0054 = 12*conj(IT_0021);
    const ccomplex_t IT_0055 = IT_0033*(96*IT_0030 + (-96)*IT_0041) + (
      -2.66666666666667)*IT_0044*(IT_0037 + (-0.5)*IT_0045) + IT_0048*IT_0049 +
       IT_0052*IT_0053 + IT_0032*(48*IT_0015 + (-48)*IT_0017 + (-12)*IT_0021 +
       IT_0054);
    const ccomplex_t IT_0056 = 24*IT_0015;
    const ccomplex_t IT_0057 = 12*IT_0021;
    const ccomplex_t IT_0058 = (-48)*IT_0015;
    const ccomplex_t IT_0059 = (-12)*conj(IT_0021);
    const ccomplex_t IT_0060 = (-10.6666666666667)*IT_0003*IT_0020*conj
      (IT_0020) + s_13*(12*IT_0017*IT_0020 + 12*IT_0015*conj(IT_0020)) + (
      -1.33333333333333)*IT_0037*(s_34*(IT_0008 + -conj(IT_0008)) + (-8)*IT_0003
      *(IT_0037 + (-0.125)*IT_0045) + (-0.75)*s_13*((-24)*IT_0017 + (
      -1.33333333333333)*IT_0021 + 1.33333333333333*conj(IT_0021) + IT_0056)) +
       IT_0004*((-10.6666666666667)*IT_0021*conj(IT_0021) + IT_0015*((-12)
      *IT_0021 + IT_0054 + IT_0056) + IT_0017*(24*IT_0017 + IT_0057 + IT_0058 +
       IT_0059));
    const ccomplex_t IT_0061 = 0.5*IT_0015;
    const ccomplex_t IT_0062 = 2.66666666666667*(IT_0037 + 0.375*IT_0038 +
       0.375*IT_0039)*IT_0044 + IT_0049*IT_0052 + IT_0048*IT_0053 + IT_0032*(48
      *IT_0017 + IT_0057 + IT_0058 + IT_0059) + IT_0033*((-96)*IT_0030 + 192
      *IT_0041 + (-96)*IT_0061);
    const ccomplex_t IT_0063 = conj(IT_0008)*(IT_0015*IT_0016 + IT_0017
      *IT_0018 + s_14*IT_0021*IT_0022 + IT_0020*IT_0023) + IT_0008*(IT_0016
      *IT_0017 + IT_0015*IT_0018 + s_14*conj(IT_0021)*IT_0022 + conj(IT_0020)
      *IT_0023) + conj(IT_0020)*(IT_0017*IT_0025 + IT_0021*IT_0026) + IT_0020*
      (IT_0015*IT_0025 + conj(IT_0021)*IT_0026) + (-2)*IT_0030*IT_0040 + IT_0041
      *IT_0055 + s_13*IT_0060 + IT_0061*IT_0062;
    return create_ccomplex_return(IT_0063);
}

