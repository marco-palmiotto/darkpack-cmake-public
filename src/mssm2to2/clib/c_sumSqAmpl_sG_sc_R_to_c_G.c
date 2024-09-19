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
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sc_R, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = (-2)*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0013 = (-2)*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*2)*IT_0008 + IT_0017;
    const ccomplex_t IT_0019 = s_34*IT_0004;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_sG*IT_0008;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0020;
    const ccomplex_t IT_0023 = m_sG*s_34;
    const ccomplex_t IT_0024 = IT_0001*IT_0012;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0005 
      + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = s_14*s_34;
    const ccomplex_t IT_0030 = (-2)*IT_0027;
    const ccomplex_t IT_0031 = s_14*IT_0003;
    const ccomplex_t IT_0032 = s_13*s_34;
    const ccomplex_t IT_0033 = (-2)*IT_0032;
    const ccomplex_t IT_0034 = IT_0031 + IT_0033;
    const ccomplex_t IT_0035 = -IT_0017;
    const ccomplex_t IT_0036 = (-21.3333333333333)*IT_0032;
    const ccomplex_t IT_0037 = 10.6666666666667*IT_0031;
    const ccomplex_t IT_0038 = IT_0036 + IT_0037;
    const ccomplex_t IT_0039 = (IT_0011 + (-4)*IT_0018 + 4*conj(IT_0018))
      *IT_0019 + (-0.375)*(74.6666666666667*IT_0021 + 74.6666666666667*IT_0022)
      *IT_0023 + (-0.375)*(10.6666666666667*IT_0009 + (-42.6666666666667)
      *IT_0010 + (-10.6666666666667)*IT_0028)*IT_0029 + -IT_0030*IT_0034 + 0.375
      *(IT_0035 + -conj(IT_0035))*IT_0038;
    const ccomplex_t IT_0040 = m_sG*s_13;
    const ccomplex_t IT_0041 = (-9.33333333333333)*IT_0021;
    const ccomplex_t IT_0042 = (-9.33333333333333)*IT_0022;
    const ccomplex_t IT_0043 = IT_0041 + IT_0042;
    const ccomplex_t IT_0044 = IT_0021 + IT_0022;
    const ccomplex_t IT_0045 = 18.6666666666667*IT_0011;
    const ccomplex_t IT_0046 = m_sG*IT_0003;
    const ccomplex_t IT_0047 = 18.6666666666667*IT_0030;
    const ccomplex_t IT_0048 = -IT_0012;
    const ccomplex_t IT_0049 = IT_0016*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = 74.6666666666667*IT_0009;
    const ccomplex_t IT_0052 = cpow(IT_0021, 2);
    const ccomplex_t IT_0053 = cpow(IT_0022, 2);
    const ccomplex_t IT_0054 = (-37.3333333333333)*s_13;
    const ccomplex_t IT_0055 = 9.33333333333333*IT_0021;
    const ccomplex_t IT_0056 = 9.33333333333333*IT_0022;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0059 = s_34*IT_0058;
    const ccomplex_t IT_0060 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0061 = pow(s_13, 2);
    const ccomplex_t IT_0062 = (-10.6666666666667)*IT_0061;
    const ccomplex_t IT_0063 = (-1.33333333333333)*IT_0061;
    const ccomplex_t IT_0064 = (-1.33333333333333)*IT_0009;
    const ccomplex_t IT_0065 = 1.33333333333333*conj(IT_0018);
    const ccomplex_t IT_0066 = 21.3333333333333*IT_0030;
    const ccomplex_t IT_0067 = (-1.33333333333333)*IT_0035;
    const ccomplex_t IT_0068 = (-10.6666666666667)*IT_0004*IT_0018*conj
      (IT_0018) + (-74.6666666666667)*IT_0021*IT_0022 + (-10.6666666666667)
      *IT_0003*IT_0035*conj(IT_0035) + 1.33333333333333*IT_0011*(s_14*conj
      (IT_0050) + 0.75*IT_0004*(10.6666666666667*IT_0011 + (-1.33333333333333)
      *IT_0018 + IT_0065) + s_13*(conj(IT_0035) + 0.75*IT_0066)) +
       1.33333333333333*IT_0030*(s_13*conj(IT_0018) + s_34*conj(IT_0050) + 0.75
      *IT_0003*(10.6666666666667*IT_0030 + 1.33333333333333*conj(IT_0035) +
       IT_0067));
    const ccomplex_t IT_0069 = s_13*s_14;
    const ccomplex_t IT_0070 = IT_0019 + (-2)*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*((-37.3333333333333)*IT_0021 + (
      -37.3333333333333)*IT_0022)*IT_0023 + (-0.5)*(85.3333333333333*IT_0009 + (
      -42.6666666666667)*IT_0028)*IT_0029 + (-0.666666666666667)*IT_0031*(conj
      (IT_0035) + 0.75*IT_0066 + 0.75*IT_0067) + IT_0011*(10.6666666666667
      *IT_0019 + (-21.3333333333333)*IT_0069) + (-0.666666666666667)*(IT_0018 + 
      -conj(IT_0018))*IT_0070;
    const ccomplex_t IT_0072 = (-2.66666666666667)*IT_0010*IT_0039 + IT_0040*
      (conj(IT_0018)*IT_0043 + IT_0044*IT_0045) + IT_0046*(conj(IT_0035)*IT_0043
       + IT_0044*IT_0047) + IT_0023*(IT_0043*conj(IT_0050) + IT_0044*IT_0051) + 
      (IT_0052 + IT_0053)*IT_0054 + IT_0050*(IT_0023*IT_0057 + conj(IT_0035)
      *IT_0059 + s_34*IT_0030*IT_0060 + s_14*(conj(IT_0018)*IT_0058 + IT_0011
      *IT_0060)) + 1.33333333333333*IT_0035*(IT_0009*IT_0034 + 0.75*IT_0046
      *IT_0057 + 0.75*conj(IT_0050)*IT_0059 + 0.75*conj(IT_0018)*IT_0062 + 0.75
      *IT_0011*IT_0063) + IT_0034*conj(IT_0035)*IT_0064 + IT_0018*(IT_0040
      *IT_0057 + s_14*conj(IT_0050)*IT_0058 + conj(IT_0035)*IT_0062 + IT_0030
      *IT_0063 + IT_0019*IT_0064) + IT_0009*((-42.6666666666667)*IT_0009*IT_0029
       + (-2)*IT_0030*IT_0038 + IT_0019*(21.3333333333333*IT_0011 + IT_0065)) +
       s_13*IT_0068 + (-2)*IT_0028*IT_0071;
    return create_ccomplex_return(IT_0072);
}

