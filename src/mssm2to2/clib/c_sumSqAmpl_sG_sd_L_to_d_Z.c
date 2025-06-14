#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sd_L_to_d_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sd_L_to_d_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 3
      *IT_0000*IT_0003);
    const ccomplex_t IT_0005 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0007 = pow(m_d, 2);
    const ccomplex_t IT_0008 = pow(m_sG, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_Z, -2);
    const ccomplex_t IT_0015 = pow(s_14, 2);
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = IT_0008 + IT_0017;
    const ccomplex_t IT_0019 = s_13*IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + 3*IT_0023);
    const ccomplex_t IT_0025 = (-0.166666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0006*IT_0025;
    const ccomplex_t IT_0027 = pow(m_sd_L, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0012 + 2*IT_0030;
    const ccomplex_t IT_0032 = s_14*s_34*IT_0014;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = s_13 + IT_0033;
    const ccomplex_t IT_0035 = s_13*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0037 = 0.333333333333333*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0040 = IT_0028*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = m_d*IT_0030;
    const ccomplex_t IT_0042 = IT_0040 + -IT_0041;
    const ccomplex_t IT_0043 = pow(m_Z, 2);
    const ccomplex_t IT_0044 = s_13*IT_0043;
    const ccomplex_t IT_0045 = s_14*s_34;
    const ccomplex_t IT_0046 = (-2)*IT_0045;
    const ccomplex_t IT_0047 = IT_0044 + IT_0046;
    const ccomplex_t IT_0048 = IT_0014*IT_0047;
    const ccomplex_t IT_0049 = m_d*IT_0018;
    const ccomplex_t IT_0050 = (-8)*conj(IT_0013);
    const ccomplex_t IT_0051 = 0.125*IT_0050;
    const ccomplex_t IT_0052 = IT_0013 + IT_0051;
    const ccomplex_t IT_0053 = m_d*IT_0034;
    const ccomplex_t IT_0054 = 8*IT_0031;
    const ccomplex_t IT_0055 = (-8)*conj(IT_0031);
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = pow(m_Z, 4);
    const ccomplex_t IT_0058 = s_13*IT_0057;
    const ccomplex_t IT_0059 = IT_0014*IT_0058;
    const ccomplex_t IT_0060 = s_14*IT_0043;
    const ccomplex_t IT_0061 = IT_0014*IT_0060;
    const ccomplex_t IT_0062 = (-0.25)*IT_0061;
    const ccomplex_t IT_0063 = s_14 + IT_0062;
    const ccomplex_t IT_0064 = m_d*IT_0063;
    const ccomplex_t IT_0065 = (-4)*IT_0064;
    const ccomplex_t IT_0066 = s_13*s_14;
    const ccomplex_t IT_0067 = s_34*IT_0008;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = s_14*IT_0007;
    const ccomplex_t IT_0071 = s_13*s_34;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = pow(s_34, 2);
    const ccomplex_t IT_0075 = IT_0014*IT_0074;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0007 + IT_0076;
    const ccomplex_t IT_0078 = s_13*IT_0077;
    const ccomplex_t IT_0079 = (-8)*IT_0013*(conj(IT_0013)*IT_0019 + conj
      (IT_0031)*IT_0035) + -IT_0042*(IT_0042*(16*s_13 + (-8)*IT_0048) + (-8)
      *IT_0049*IT_0052 + -IT_0053*IT_0056) + -IT_0030*(IT_0030*((-16)*IT_0044 +
       32*IT_0045 + 8*IT_0059) + 16*IT_0042*IT_0065 + (-8)*IT_0052*IT_0069 + 
      -IT_0056*IT_0073) + (-8)*IT_0031*((-0.125)*IT_0035*IT_0050 + conj(IT_0031)
      *IT_0078);
    return create_ccomplex_return(IT_0079);
}

