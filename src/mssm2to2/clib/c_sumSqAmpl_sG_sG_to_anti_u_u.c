#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sG_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sG_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_u = param->m_u;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0004;
    const ccomplex_t IT_0012 = s_13*s_24;
    const ccomplex_t IT_0013 = -g_s;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = cpow(s_12 + IT_0004 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (-0.25)*IT_0019;
    const ccomplex_t IT_0021 = 0.25*IT_0019;
    const ccomplex_t IT_0022 = s_14*s_23;
    const ccomplex_t IT_0023 = s_12*IT_0003;
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0003 + IT_0004 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = 0.5*s_34;
    const ccomplex_t IT_0029 = IT_0003 + IT_0028;
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = 21.3333333333333*IT_0023;
    const ccomplex_t IT_0032 = 42.6666666666667*IT_0022;
    const ccomplex_t IT_0033 = IT_0031 + IT_0032;
    const ccomplex_t IT_0034 = 2*IT_0012;
    const ccomplex_t IT_0035 = IT_0023 + IT_0034;
    const ccomplex_t IT_0036 = -conj(IT_0010);
    const ccomplex_t IT_0037 = IT_0007 + -conj(IT_0007) + IT_0010 + IT_0036;
    const ccomplex_t IT_0038 = 85.3333333333333*IT_0012;
    const ccomplex_t IT_0039 = 85.3333333333333*IT_0023;
    const ccomplex_t IT_0040 = 85.3333333333333*IT_0022;
    const ccomplex_t IT_0041 = IT_0038 + IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = (-8)*IT_0041;
    const ccomplex_t IT_0043 = (-16)*s_34 + (-32)*IT_0003;
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = 21.3333333333333*IT_0044;
    const ccomplex_t IT_0046 = 16*s_34 + 32*IT_0003;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = (-21.3333333333333)*IT_0047;
    const ccomplex_t IT_0049 = IT_0042 + IT_0045 + IT_0048;
    const ccomplex_t IT_0050 = (-16)*IT_0021*((-10.6666666666667)*IT_0012 + (
      -10.6666666666667)*IT_0022 + (-10.6666666666667)*IT_0023) + 
      (10.6666666666667*IT_0007 + (-10.6666666666667)*conj(IT_0007) +
       10.6666666666667*IT_0010 + (-10.6666666666667)*conj(IT_0010) +
       341.333333333333*IT_0021 + (-85.3333333333333)*IT_0025 + 85.3333333333333
      *conj(IT_0025) + (-85.3333333333333)*IT_0027 + 85.3333333333333*conj
      (IT_0027))*IT_0030 + 2*(conj(IT_0025) + conj(IT_0027))*IT_0033 +
       5.33333333333333*IT_0035*IT_0037 + IT_0020*IT_0049;
    const ccomplex_t IT_0051 = 42.6666666666667*IT_0012;
    const ccomplex_t IT_0052 = IT_0031 + IT_0051;
    const ccomplex_t IT_0053 = s_34*IT_0004;
    const ccomplex_t IT_0054 = 21.3333333333333*IT_0025;
    const ccomplex_t IT_0055 = (-1.33333333333333)*IT_0010;
    const ccomplex_t IT_0056 = 0.5*IT_0023;
    const ccomplex_t IT_0057 = IT_0022 + IT_0056;
    const ccomplex_t IT_0058 = IT_0030 + IT_0057;
    const ccomplex_t IT_0059 = (-10.6666666666667)*IT_0021;
    const ccomplex_t IT_0060 = (-0.75)*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (-1.33333333333333)*conj(IT_0007);
    const ccomplex_t IT_0062 = (-1.33333333333333)*conj(IT_0010);
    const ccomplex_t IT_0063 = 10.6666666666667*IT_0021;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = (-2)*IT_0020;
    const ccomplex_t IT_0066 = IT_0033*IT_0065;
    const ccomplex_t IT_0067 = conj(IT_0007)*(21.3333333333333*IT_0010*IT_0011
       + 21.3333333333333*IT_0007*IT_0012) + conj(IT_0010)*(21.3333333333333
      *IT_0007*IT_0011 + 21.3333333333333*IT_0010*IT_0012) + IT_0020*IT_0050 +
       IT_0021*(((-85.3333333333333)*IT_0007 + 85.3333333333333*conj(IT_0007) + 
      (-85.3333333333333)*IT_0010 + 85.3333333333333*conj(IT_0010))*IT_0030 +
       IT_0021*IT_0049 + (-2)*IT_0037*IT_0052) + (-1.33333333333333)*conj
      (IT_0027)*(IT_0007*IT_0053 + (-0.75)*IT_0011*IT_0054 + (-0.75)*IT_0023
      *IT_0055 + IT_0060) + (-1.33333333333333)*conj(IT_0025)*(IT_0007*IT_0023 +
       (-0.75)*IT_0022*IT_0054 + (-0.75)*IT_0053*IT_0055 + IT_0060) + IT_0027*
      (21.3333333333333*IT_0011*conj(IT_0025) + 21.3333333333333*IT_0022*conj
      (IT_0027) + IT_0053*IT_0061 + IT_0023*IT_0062 + IT_0064 + IT_0066) +
       IT_0025*(IT_0023*IT_0061 + IT_0053*IT_0062 + IT_0064 + IT_0066);
    return create_ccomplex_return(IT_0067);
}

