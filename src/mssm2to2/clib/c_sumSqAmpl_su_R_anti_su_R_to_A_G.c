#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_A_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_A_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s*e_em;
    const ccomplex_t IT_0001 = 1.33333333333333*IT_0000;
    const ccomplex_t IT_0002 = pow(m_su_R, 2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0006 = (-1.33333333333333)*IT_0005;
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_su_R
      *Gamma_ur + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_su_R
      *Gamma_ur + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-2)*IT_0003;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = IT_0013*(IT_0007 + -IT_0015);
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = conj(IT_0011) + conj(IT_0017);
    const ccomplex_t IT_0019 = IT_0007*IT_0013;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0009*(IT_0007 + -IT_0015);
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(IT_0020) + conj(IT_0022);
    const ccomplex_t IT_0024 = 4*s_12;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0027 = IT_0004*IT_0026;
    const ccomplex_t IT_0028 = IT_0013*IT_0027;
    const ccomplex_t IT_0029 = -IT_0027;
    const ccomplex_t IT_0030 = IT_0007 + IT_0029;
    const ccomplex_t IT_0031 = IT_0009*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0028 + 0.5*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0010;
    const ccomplex_t IT_0034 = conj(IT_0032) + conj(IT_0033);
    const ccomplex_t IT_0035 = 4*s_13;
    const ccomplex_t IT_0036 = IT_0001*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0016;
    const ccomplex_t IT_0038 = IT_0014*IT_0026;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0027 + IT_0039;
    const ccomplex_t IT_0041 = IT_0013*IT_0040;
    const ccomplex_t IT_0042 = IT_0007 + IT_0038;
    const ccomplex_t IT_0043 = -IT_0015 + -IT_0027;
    const ccomplex_t IT_0044 = IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0041 + 0.5*IT_0045;
    const ccomplex_t IT_0047 = conj(IT_0037) + conj(IT_0046);
    const ccomplex_t IT_0048 = 4*s_23;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0031 + 0.5*IT_0041;
    const ccomplex_t IT_0051 = s_23*IT_0035;
    const ccomplex_t IT_0052 = pow(s_13, 2);
    const ccomplex_t IT_0053 = 4*IT_0052;
    const ccomplex_t IT_0054 = pow(s_23, 2);
    const ccomplex_t IT_0055 = 4*IT_0054;
    const ccomplex_t IT_0056 = s_13*IT_0024;
    const ccomplex_t IT_0057 = s_23*IT_0024;
    const ccomplex_t IT_0058 = IT_0002*IT_0035;
    const ccomplex_t IT_0059 = IT_0002*IT_0048;
    const ccomplex_t IT_0060 = (-4)*s_13*IT_0001;
    const ccomplex_t IT_0061 = IT_0002*IT_0024;
    const ccomplex_t IT_0062 = pow(s_12, 2);
    const ccomplex_t IT_0063 = 4*IT_0062;
    const ccomplex_t IT_0064 = pow(m_su_R, 4);
    const ccomplex_t IT_0065 = 4*IT_0064;
    const ccomplex_t IT_0066 = (-4)*IT_0001*IT_0002;
    const ccomplex_t IT_0067 = s_12*IT_0001;
    const ccomplex_t IT_0068 = (-0.25)*conj(IT_0032);
    const ccomplex_t IT_0069 = (-0.25)*conj(IT_0033);
    const ccomplex_t IT_0070 = (-0.25)*conj(IT_0037);
    const ccomplex_t IT_0071 = (-0.25)*conj(IT_0050)*IT_0051;
    const ccomplex_t IT_0072 = (-0.25)*conj(IT_0022);
    const ccomplex_t IT_0073 = (-0.25)*IT_0018*IT_0061;
    const ccomplex_t IT_0074 = (-0.25)*conj(IT_0046);
    const ccomplex_t IT_0075 = conj(IT_0020)*IT_0065 + (-4)*IT_0067 + (-4)
      *IT_0059*IT_0068 + (-4)*IT_0056*IT_0069 + (-4)*IT_0058*IT_0070 + (-4)
      *IT_0071 + (-4)*IT_0063*IT_0072 + (-4)*IT_0073 + (-4)*IT_0057*IT_0074;
    const ccomplex_t IT_0076 = conj(IT_0020)*IT_0063 + (-4)*IT_0067 + (-4)
      *IT_0056*IT_0068 + (-4)*IT_0059*IT_0069 + (-4)*IT_0057*IT_0070 + (-4)
      *IT_0071 + (-4)*IT_0065*IT_0072 + (-4)*IT_0073 + (-4)*IT_0058*IT_0074;
    const ccomplex_t IT_0077 = conj(IT_0011)*IT_0056;
    const ccomplex_t IT_0078 = (-4)*s_23*IT_0001;
    const ccomplex_t IT_0079 = conj(IT_0017)*IT_0059;
    const ccomplex_t IT_0080 = IT_0001*((-16)*IT_0001 + 4*IT_0002*IT_0018) +
       IT_0023*IT_0025 + IT_0034*IT_0036 + IT_0047*IT_0049 + IT_0050*(IT_0023
      *IT_0051 + conj(IT_0011)*IT_0053 + conj(IT_0017)*IT_0055) + IT_0032*(conj
      (IT_0037)*IT_0051 + conj(IT_0033)*IT_0053 + conj(IT_0022)*IT_0056 + conj
      (IT_0017)*IT_0057 + conj(IT_0011)*IT_0058 + conj(IT_0020)*IT_0059 +
       IT_0060) + IT_0033*(conj(IT_0046)*IT_0051 + conj(IT_0032)*IT_0053 + conj
      (IT_0020)*IT_0056 + conj(IT_0017)*IT_0057 + conj(IT_0011)*IT_0058 + conj
      (IT_0022)*IT_0059 + IT_0060) + IT_0011*(conj(IT_0050)*IT_0053 + IT_0047
      *IT_0056 + IT_0034*IT_0058 + IT_0023*IT_0061 + conj(IT_0017)*IT_0063 +
       conj(IT_0011)*IT_0065 + IT_0066) + IT_0017*(conj(IT_0050)*IT_0055 +
       IT_0034*IT_0057 + IT_0047*IT_0059 + IT_0023*IT_0061 + conj(IT_0011)
      *IT_0063 + conj(IT_0017)*IT_0065 + IT_0066) + IT_0020*IT_0075 + IT_0022
      *IT_0076 + IT_0037*(conj(IT_0032)*IT_0051 + conj(IT_0046)*IT_0055 + conj
      (IT_0022)*IT_0057 + conj(IT_0020)*IT_0058 + IT_0077 + IT_0078 + IT_0079) +
       IT_0046*(conj(IT_0033)*IT_0051 + conj(IT_0037)*IT_0055 + conj(IT_0020)
      *IT_0057 + conj(IT_0022)*IT_0058 + IT_0077 + IT_0078 + IT_0079);
    return create_ccomplex_return(IT_0080);
}

