#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sG_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sG_to_anti_d_d(
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
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_d, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = s_13*s_24;
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0003*IT_0004;
    const ccomplex_t IT_0013 = -g_s;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = cpow(s_12 + IT_0004 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.25*IT_0019;
    const ccomplex_t IT_0021 = (-0.25)*IT_0019;
    const ccomplex_t IT_0022 = s_12*IT_0003;
    const ccomplex_t IT_0023 = s_14*s_23;
    const ccomplex_t IT_0024 = 21.3333333333333*IT_0022;
    const ccomplex_t IT_0025 = 42.6666666666667*IT_0008;
    const ccomplex_t IT_0026 = IT_0024 + IT_0025;
    const ccomplex_t IT_0027 = (-16)*s_34 + (-32)*IT_0003;
    const ccomplex_t IT_0028 = IT_0004*IT_0027;
    const ccomplex_t IT_0029 = 85.3333333333333*IT_0023;
    const ccomplex_t IT_0030 = 85.3333333333333*IT_0022;
    const ccomplex_t IT_0031 = 85.3333333333333*IT_0008;
    const ccomplex_t IT_0032 = IT_0029 + IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = (-8)*IT_0032;
    const ccomplex_t IT_0034 = 16*s_34 + 32*IT_0003;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = (-21.3333333333333)*IT_0035;
    const ccomplex_t IT_0037 = 21.3333333333333*IT_0028 + IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0022;
    const ccomplex_t IT_0039 = IT_0008 + IT_0038;
    const ccomplex_t IT_0040 = 10.6666666666667*IT_0011;
    const ccomplex_t IT_0041 = (-10.6666666666667)*conj(IT_0011) + IT_0040;
    const ccomplex_t IT_0042 = 10.6666666666667*IT_0007*IT_0021;
    const ccomplex_t IT_0043 = (-10.6666666666667)*conj(IT_0007)*IT_0021;
    const ccomplex_t IT_0044 = 0.5*s_34;
    const ccomplex_t IT_0045 = IT_0003 + IT_0044;
    const ccomplex_t IT_0046 = IT_0004*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0002*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0002*IT_0049;
    const ccomplex_t IT_0051 = 10.6666666666667*IT_0020;
    const ccomplex_t IT_0052 = (-85.3333333333333)*IT_0021;
    const ccomplex_t IT_0053 = 85.3333333333333*IT_0021;
    const ccomplex_t IT_0054 = (-10.6666666666667)*IT_0020;
    const ccomplex_t IT_0055 = IT_0020*((-85.3333333333333)*IT_0007 +
       85.3333333333333*conj(IT_0007) + (-85.3333333333333)*IT_0011 +
       85.3333333333333*conj(IT_0011) + 341.333333333333*IT_0021) + IT_0021
      *IT_0041 + IT_0042 + IT_0043 + (IT_0048 + IT_0050)*(IT_0051 + IT_0052) + 
      (conj(IT_0048) + conj(IT_0050))*(IT_0053 + IT_0054);
    const ccomplex_t IT_0056 = s_34*IT_0004;
    const ccomplex_t IT_0057 = 42.6666666666667*IT_0023;
    const ccomplex_t IT_0058 = IT_0024 + IT_0057;
    const ccomplex_t IT_0059 = 2*IT_0021;
    const ccomplex_t IT_0060 = 21.3333333333333*IT_0048;
    const ccomplex_t IT_0061 = 2*IT_0023;
    const ccomplex_t IT_0062 = IT_0022 + IT_0061;
    const ccomplex_t IT_0063 = (-5.33333333333333)*IT_0020;
    const ccomplex_t IT_0064 = (-1.33333333333333)*IT_0007;
    const ccomplex_t IT_0065 = 5.33333333333333*IT_0020;
    const ccomplex_t IT_0066 = (-2)*IT_0021;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = IT_0062*IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = (-1.33333333333333)*conj(IT_0011);
    const ccomplex_t IT_0070 = (-1.33333333333333)*conj(IT_0007);
    const ccomplex_t IT_0071 = IT_0011*(21.3333333333333*IT_0008*conj(IT_0011)
       + 21.3333333333333*conj(IT_0007)*IT_0012) + IT_0007*(21.3333333333333
      *conj(IT_0007)*IT_0008 + 21.3333333333333*conj(IT_0011)*IT_0012) + IT_0020
      *((-16)*IT_0021*((-10.6666666666667)*IT_0008 + (-10.6666666666667)*IT_0022
       + (-10.6666666666667)*IT_0023) + (-2)*(IT_0007 + -conj(IT_0007) + IT_0011
       + -conj(IT_0011))*IT_0026 + IT_0020*IT_0037) + IT_0021*(IT_0021*IT_0037 +
       IT_0039*IT_0041) + IT_0039*(IT_0042 + IT_0043) + IT_0046*IT_0055 + conj
      (IT_0048)*(21.3333333333333*IT_0012*IT_0050 + (-1.33333333333333)*IT_0011
      *IT_0056 + IT_0058*IT_0059 + IT_0023*IT_0060 + IT_0062*IT_0063 + IT_0022
      *IT_0064) + (-1.33333333333333)*conj(IT_0050)*(IT_0011*IT_0022 + (-16)
      *IT_0023*IT_0050 + (-0.75)*IT_0058*IT_0059 + (-0.75)*IT_0012*IT_0060 + (
      -0.75)*IT_0062*IT_0063 + (-0.75)*IT_0056*IT_0064) + IT_0048*(IT_0068 +
       IT_0056*IT_0069 + IT_0022*IT_0070) + IT_0050*(IT_0068 + IT_0022*IT_0069 +
       IT_0056*IT_0070);
    return create_ccomplex_return(IT_0071);
}

