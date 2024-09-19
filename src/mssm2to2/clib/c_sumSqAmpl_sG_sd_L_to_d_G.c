#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sd_L_to_d_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sd_L_to_d_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_d = param->m_d;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_d, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sd_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*m_sG*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = 0.5*IT_0009;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = (-37.3333333333333)*s_13;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = m_sG*s_34;
    const ccomplex_t IT_0022 = (-9.33333333333333)*IT_0010;
    const ccomplex_t IT_0023 = (-9.33333333333333)*IT_0012;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = (-2)*IT_0015;
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*-2)*IT_0008 + -IT_0026;
    const ccomplex_t IT_0028 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0029 = s_34*IT_0028;
    const ccomplex_t IT_0030 = pow(s_13, 2);
    const ccomplex_t IT_0031 = (-10.6666666666667)*IT_0030;
    const ccomplex_t IT_0032 = s_34*IT_0004;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0034 = (-10.6666666666667)*IT_0033;
    const ccomplex_t IT_0035 = 9.33333333333333*IT_0010;
    const ccomplex_t IT_0036 = 9.33333333333333*IT_0012;
    const ccomplex_t IT_0037 = IT_0035 + IT_0036;
    const ccomplex_t IT_0038 = IT_0001*IT_0015;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0005 
      + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0041;
    const ccomplex_t IT_0043 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0044 = 2*IT_0033;
    const ccomplex_t IT_0045 = (-1.33333333333333)*IT_0030;
    const ccomplex_t IT_0046 = 10.6666666666667*IT_0033;
    const ccomplex_t IT_0047 = s_14*s_34;
    const ccomplex_t IT_0048 = s_14*IT_0003;
    const ccomplex_t IT_0049 = s_13*s_34;
    const ccomplex_t IT_0050 = (-2)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = (-21.3333333333333)*IT_0049;
    const ccomplex_t IT_0053 = 10.6666666666667*IT_0048;
    const ccomplex_t IT_0054 = IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = IT_0026 + -conj(IT_0026);
    const ccomplex_t IT_0056 = (74.6666666666667*IT_0010 + 74.6666666666667
      *IT_0012)*IT_0021;
    const ccomplex_t IT_0057 = -IT_0033;
    const ccomplex_t IT_0058 = 1.33333333333333*conj(IT_0027);
    const ccomplex_t IT_0059 = IT_0042*IT_0054 + (-0.666666666666667)*IT_0051
      *IT_0055 + (-0.5)*IT_0056 + (-0.5)*IT_0047*(85.3333333333333*IT_0041 +
       IT_0046 + (-42.6666666666667)*IT_0057) + (-0.5)*IT_0032*((
      -1.33333333333333)*IT_0027 + 21.3333333333333*IT_0044 + IT_0058);
    const ccomplex_t IT_0060 = m_sG*s_13;
    const ccomplex_t IT_0061 = 18.6666666666667*IT_0012;
    const ccomplex_t IT_0062 = m_sG*IT_0003;
    const ccomplex_t IT_0063 = 21.3333333333333*IT_0042;
    const ccomplex_t IT_0064 = (-1.33333333333333)*IT_0026;
    const ccomplex_t IT_0065 = (-74.6666666666667)*IT_0010*IT_0012 + (
      -10.6666666666667)*IT_0003*IT_0026*conj(IT_0026) + (-10.6666666666667)
      *IT_0004*IT_0027*conj(IT_0027) + 1.33333333333333*IT_0044*(s_14*conj
      (IT_0020) + 0.75*IT_0004*((-1.33333333333333)*IT_0027 + 10.6666666666667
      *IT_0044 + IT_0058) + s_13*(conj(IT_0026) + 0.75*IT_0063)) +
       1.33333333333333*IT_0042*(s_34*conj(IT_0020) + s_13*conj(IT_0027) + 0.75
      *IT_0003*(1.33333333333333*conj(IT_0026) + 10.6666666666667*IT_0042 +
       IT_0064));
    const ccomplex_t IT_0066 = s_13*s_14;
    const ccomplex_t IT_0067 = IT_0032 + (-2)*IT_0066;
    const ccomplex_t IT_0068 = ((-37.3333333333333)*IT_0010 + (
      -37.3333333333333)*IT_0012)*IT_0021 + (IT_0034 + (-42.6666666666667)
      *IT_0041)*IT_0047 + 1.33333333333333*IT_0048*(conj(IT_0026) + 0.75*IT_0063
       + 0.75*IT_0064) + (-2)*IT_0044*(10.6666666666667*IT_0032 + (
      -21.3333333333333)*IT_0066) + 1.33333333333333*(IT_0027 + -conj(IT_0027))
      *IT_0067;
    const ccomplex_t IT_0069 = (IT_0011 + IT_0013)*IT_0014 + conj(IT_0020)*
      (IT_0021*IT_0024 + s_14*IT_0027*IT_0028 + IT_0026*IT_0029) + conj(IT_0027)
      *(IT_0026*IT_0031 + IT_0032*IT_0034) + IT_0020*(conj(IT_0026)*IT_0029 +
       IT_0021*IT_0037 + s_34*IT_0042*IT_0043 + s_14*(conj(IT_0027)*IT_0028 +
       IT_0043*IT_0044)) + IT_0026*IT_0044*IT_0045 + IT_0027*(conj(IT_0026)
      *IT_0031 + IT_0042*IT_0045 + IT_0032*IT_0046) + (-2.66666666666667)
      *IT_0033*(IT_0032*IT_0044 + 16*IT_0033*IT_0047 + -IT_0042*IT_0051 + 0.375
      *IT_0054*IT_0055 + (-0.375)*IT_0056) + (-2)*IT_0057*IT_0059 + IT_0060*
      (IT_0024*conj(IT_0027) + IT_0027*IT_0037 + IT_0044*(18.6666666666667
      *IT_0010 + IT_0061)) + (IT_0024*conj(IT_0026) + IT_0026*IT_0037 + IT_0042*
      (18.6666666666667*IT_0010 + IT_0061))*IT_0062 + s_13*IT_0065 + IT_0041
      *IT_0068;
    return create_ccomplex_return(IT_0069);
}

