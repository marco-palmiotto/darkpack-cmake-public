#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sG_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sG_to_anti_c_c(
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
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0004;
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = s_12*IT_0003;
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0002*IT_0015;
    const ccomplex_t IT_0017 = s_34*IT_0004;
    const ccomplex_t IT_0018 = s_13*s_24;
    const ccomplex_t IT_0019 = s_14*s_23;
    const ccomplex_t IT_0020 = 0.5*IT_0014;
    const ccomplex_t IT_0021 = IT_0019 + IT_0020;
    const ccomplex_t IT_0022 = -g_s;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cpow(s_12 + IT_0004 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.25*IT_0028;
    const ccomplex_t IT_0030 = 10.6666666666667*IT_0029;
    const ccomplex_t IT_0031 = IT_0016*IT_0030;
    const ccomplex_t IT_0032 = (-10.6666666666667)*IT_0029;
    const ccomplex_t IT_0033 = conj(IT_0016)*IT_0032;
    const ccomplex_t IT_0034 = 0.5*s_34;
    const ccomplex_t IT_0035 = IT_0003 + IT_0034;
    const ccomplex_t IT_0036 = IT_0004*IT_0035;
    const ccomplex_t IT_0037 = (-0.25)*IT_0028;
    const ccomplex_t IT_0038 = 85.3333333333333*IT_0037;
    const ccomplex_t IT_0039 = (-85.3333333333333)*IT_0037;
    const ccomplex_t IT_0040 = ((-85.3333333333333)*IT_0007 + 85.3333333333333
      *conj(IT_0007) + (-85.3333333333333)*IT_0010 + 85.3333333333333*conj
      (IT_0010))*IT_0029 + IT_0031 + IT_0033 + (10.6666666666667*IT_0007 + (
      -10.6666666666667)*conj(IT_0007) + 10.6666666666667*IT_0010 + (
      -10.6666666666667)*conj(IT_0010) + 341.333333333333*IT_0029)*IT_0037 +
       conj(IT_0016)*IT_0038 + conj(IT_0013)*(IT_0032 + IT_0038) + IT_0016
      *IT_0039 + IT_0013*(IT_0030 + IT_0039);
    const ccomplex_t IT_0041 = 21.3333333333333*IT_0013;
    const ccomplex_t IT_0042 = (-1.33333333333333)*IT_0010;
    const ccomplex_t IT_0043 = 21.3333333333333*IT_0014;
    const ccomplex_t IT_0044 = 42.6666666666667*IT_0019;
    const ccomplex_t IT_0045 = IT_0043 + IT_0044;
    const ccomplex_t IT_0046 = (-16)*s_34 + (-32)*IT_0003;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = 85.3333333333333*IT_0018;
    const ccomplex_t IT_0049 = 85.3333333333333*IT_0014;
    const ccomplex_t IT_0050 = 85.3333333333333*IT_0019;
    const ccomplex_t IT_0051 = IT_0048 + IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = (-8)*IT_0051;
    const ccomplex_t IT_0053 = 16*s_34 + 32*IT_0003;
    const ccomplex_t IT_0054 = IT_0004*IT_0053;
    const ccomplex_t IT_0055 = (-21.3333333333333)*IT_0054;
    const ccomplex_t IT_0056 = 21.3333333333333*IT_0047 + IT_0052 + IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0018;
    const ccomplex_t IT_0058 = IT_0014 + IT_0057;
    const ccomplex_t IT_0059 = -conj(IT_0010);
    const ccomplex_t IT_0060 = (-16)*((-10.6666666666667)*IT_0014 + (
      -10.6666666666667)*IT_0018 + (-10.6666666666667)*IT_0019)*IT_0029 + 2*
      (conj(IT_0013) + conj(IT_0016))*IT_0045 + IT_0037*IT_0056 +
       5.33333333333333*IT_0058*(IT_0007 + -conj(IT_0007) + IT_0010 + IT_0059);
    const ccomplex_t IT_0061 = 42.6666666666667*IT_0018;
    const ccomplex_t IT_0062 = IT_0043 + IT_0061;
    const ccomplex_t IT_0063 = (-1.33333333333333)*conj(IT_0007);
    const ccomplex_t IT_0064 = (-1.33333333333333)*conj(IT_0010);
    const ccomplex_t IT_0065 = (-2)*IT_0037;
    const ccomplex_t IT_0066 = IT_0045*IT_0065;
    const ccomplex_t IT_0067 = IT_0007*(21.3333333333333*conj(IT_0010)*IT_0011
       + (-1.33333333333333)*conj(IT_0013)*IT_0014 + (-1.33333333333333)*conj
      (IT_0016)*IT_0017 + 21.3333333333333*conj(IT_0007)*IT_0018) + IT_0010*
      (21.3333333333333*conj(IT_0007)*IT_0011 + 21.3333333333333*conj(IT_0010)
      *IT_0018) + IT_0021*(IT_0031 + conj(IT_0013)*IT_0032 + IT_0033) + IT_0036
      *IT_0040 + conj(IT_0016)*(IT_0011*IT_0041 + IT_0014*IT_0042) + conj
      (IT_0013)*(IT_0019*IT_0041 + IT_0017*IT_0042) + IT_0037*IT_0060 + IT_0029*
      (IT_0029*IT_0056 + (-2)*(IT_0007 + -conj(IT_0007) + IT_0010 + IT_0059)
      *IT_0062) + IT_0016*(21.3333333333333*IT_0011*conj(IT_0013) +
       21.3333333333333*conj(IT_0016)*IT_0019 + IT_0017*IT_0063 + IT_0014
      *IT_0064 + IT_0066) + IT_0013*(IT_0021*IT_0030 + IT_0014*IT_0063 + IT_0017
      *IT_0064 + IT_0066);
    return create_ccomplex_return(IT_0067);
}

