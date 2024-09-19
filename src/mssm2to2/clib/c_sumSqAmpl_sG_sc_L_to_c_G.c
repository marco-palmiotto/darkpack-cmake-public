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
    const ccomplex_t IT_0000 = m_sG*s_13;
    const ccomplex_t IT_0001 = -g_s;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_c, 2);
    const ccomplex_t IT_0005 = pow(m_sG, 2);
    const ccomplex_t IT_0006 = pow(m_sc_L, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0004 + -IT_0005 + IT_0006 +
       reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*m_sG*IT_0009;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0012;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0016 = (-2)*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0004 + IT_0005 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*-1)*IT_0009 + -IT_0020;
    const ccomplex_t IT_0022 = 9.33333333333333*IT_0013;
    const ccomplex_t IT_0023 = 9.33333333333333*IT_0014;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = (-9.33333333333333)*IT_0013;
    const ccomplex_t IT_0026 = (-9.33333333333333)*IT_0014;
    const ccomplex_t IT_0027 = IT_0025 + IT_0026;
    const ccomplex_t IT_0028 = m_sG*IT_0004;
    const ccomplex_t IT_0029 = IT_0002*IT_0015;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + IT_0004 + -IT_0005 + -IT_0006 
      + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = 2*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0013, 2);
    const ccomplex_t IT_0035 = cpow(IT_0014, 2);
    const ccomplex_t IT_0036 = (-37.3333333333333)*s_13;
    const ccomplex_t IT_0037 = -IT_0015;
    const ccomplex_t IT_0038 = IT_0019*IT_0037;
    const ccomplex_t IT_0039 = m_sG*s_34;
    const ccomplex_t IT_0040 = IT_0011 + IT_0021;
    const ccomplex_t IT_0041 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0042 = s_14*IT_0041;
    const ccomplex_t IT_0043 = s_34*IT_0041;
    const ccomplex_t IT_0044 = pow(s_13, 2);
    const ccomplex_t IT_0045 = (-10.6666666666667)*IT_0044;
    const ccomplex_t IT_0046 = (-10.6666666666667)*s_13*IT_0005;
    const ccomplex_t IT_0047 = (-2.66666666666667)*IT_0044;
    const ccomplex_t IT_0048 = IT_0038 + -conj(IT_0038);
    const ccomplex_t IT_0049 = (-5.33333333333333)*IT_0011;
    const ccomplex_t IT_0050 = (-2.66666666666667)*IT_0021;
    const ccomplex_t IT_0051 = 1.33333333333333*conj(IT_0020);
    const ccomplex_t IT_0052 = (-1.33333333333333)*IT_0020;
    const ccomplex_t IT_0053 = (-74.6666666666667)*IT_0013*IT_0014 + (
      -10.6666666666667)*IT_0004*IT_0020*conj(IT_0020) + IT_0011*
      (10.6666666666667*s_13*IT_0020 + IT_0005*(10.6666666666667*IT_0011 +
       10.6666666666667*IT_0021) + 10.6666666666667*s_14*IT_0038) + IT_0010*
      (s_13*(2.66666666666667*conj(IT_0020) + 42.6666666666667*IT_0033) + (
      -2.66666666666667)*s_14*IT_0048 + IT_0005*(42.6666666666667*IT_0010 +
       2.66666666666667*conj(IT_0021) + IT_0049 + IT_0050)) + (-1.33333333333333
      )*IT_0033*(s_13*(IT_0021 + -conj(IT_0021)) + s_34*IT_0048 + (-0.75)
      *IT_0004*(10.6666666666667*IT_0033 + IT_0051 + IT_0052));
    const ccomplex_t IT_0054 = s_34*IT_0005;
    const ccomplex_t IT_0055 = s_13*s_14;
    const ccomplex_t IT_0056 = s_14*s_34;
    const ccomplex_t IT_0057 = (-2)*IT_0055;
    const ccomplex_t IT_0058 = IT_0054 + IT_0057;
    const ccomplex_t IT_0059 = s_14*IT_0004;
    const ccomplex_t IT_0060 = ((-37.3333333333333)*IT_0013 + (
      -37.3333333333333)*IT_0014)*IT_0039 + (-4)*IT_0010*(10.6666666666667
      *IT_0054 + (-21.3333333333333)*IT_0055) + (-42.6666666666667)*IT_0032
      *IT_0056 + 2.66666666666667*(IT_0011 + 0.5*IT_0021 + (-0.5)*conj(IT_0021))
      *IT_0058 + (21.3333333333333*IT_0033 + IT_0051 + IT_0052)*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0010;
    const ccomplex_t IT_0062 = s_13*s_34;
    const ccomplex_t IT_0063 = (-2)*IT_0062;
    const ccomplex_t IT_0064 = IT_0059 + IT_0063;
    const ccomplex_t IT_0065 = (-21.3333333333333)*IT_0062;
    const ccomplex_t IT_0066 = 10.6666666666667*IT_0059;
    const ccomplex_t IT_0067 = IT_0065 + IT_0066;
    const ccomplex_t IT_0068 = IT_0020 + -conj(IT_0020);
    const ccomplex_t IT_0069 = 149.333333333333*IT_0014;
    const ccomplex_t IT_0070 = IT_0039*(149.333333333333*IT_0013 + IT_0069);
    const ccomplex_t IT_0071 = ((-10.6666666666667)*IT_0010 + 42.6666666666667
      *IT_0011 + 21.3333333333333*IT_0021 + (-21.3333333333333)*conj(IT_0021))
      *IT_0054 + IT_0056*((-21.3333333333333)*IT_0032 + (-170.666666666667)
      *IT_0061) + 5.33333333333333*IT_0033*IT_0064 + (-2)*IT_0067*IT_0068 +
       IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0010;
    const ccomplex_t IT_0073 = (85.3333333333333*IT_0010 + 2.66666666666667
      *conj(IT_0021) + IT_0049 + IT_0050)*IT_0054 + (-4)*IT_0033*IT_0067 +
       2.66666666666667*IT_0064*IT_0068 + IT_0070 + IT_0056*(170.666666666667
      *IT_0032 + 42.6666666666667*IT_0061 + (-170.666666666667)*IT_0072);
    const ccomplex_t IT_0074 = IT_0000*((37.3333333333333*IT_0010 +
       18.6666666666667*IT_0011)*(IT_0013 + IT_0014) + IT_0021*IT_0024 + conj
      (IT_0021)*IT_0027) + IT_0028*(IT_0020*IT_0024 + conj(IT_0020)*IT_0027 + 
      (18.6666666666667*IT_0013 + 18.6666666666667*IT_0014)*IT_0033) + (IT_0034 
      + IT_0035)*IT_0036 + conj(IT_0038)*(IT_0027*IT_0039 + IT_0040*IT_0042 +
       IT_0020*IT_0043) + IT_0038*(IT_0024*IT_0039 + conj(IT_0021)*IT_0042 +
       conj(IT_0020)*IT_0043) + conj(IT_0020)*IT_0040*IT_0045 + conj(IT_0021)*
      (IT_0020*IT_0045 + IT_0040*IT_0046) + (IT_0010*IT_0020 + IT_0011*IT_0033)
      *IT_0047 + s_13*IT_0053 + IT_0032*IT_0060 + IT_0061*IT_0071 + IT_0072
      *IT_0073;
    return create_ccomplex_return(IT_0074);
}

