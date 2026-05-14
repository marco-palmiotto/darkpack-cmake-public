#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_su_L_to_d_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_su_L_to_d_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = pow(m_d, 2);
    const ccomplex_t IT_0006 = pow(m_sG, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -2);
    const ccomplex_t IT_0013 = pow(s_14, 2);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0006 + IT_0015;
    const ccomplex_t IT_0017 = s_13*IT_0016;
    const ccomplex_t IT_0018 = pow(m_u, 2);
    const ccomplex_t IT_0019 = pow(m_su_L, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + -IT_0006 + IT_0018 + -IT_0019 
      + -reg_prop, -1);
    const ccomplex_t IT_0021 = 0.5*IT_0002;
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0010 + 2*IT_0024;
    const ccomplex_t IT_0026 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = s_13 + IT_0027;
    const ccomplex_t IT_0029 = s_13*IT_0028;
    const ccomplex_t IT_0030 = m_d*IT_0024;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, 2);
    const ccomplex_t IT_0033 = s_13*IT_0032;
    const ccomplex_t IT_0034 = s_14*s_34;
    const ccomplex_t IT_0035 = (-2)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033 + IT_0035;
    const ccomplex_t IT_0037 = IT_0012*IT_0036;
    const ccomplex_t IT_0038 = s_14*IT_0032;
    const ccomplex_t IT_0039 = IT_0012*IT_0038;
    const ccomplex_t IT_0040 = (-0.25)*IT_0039;
    const ccomplex_t IT_0041 = s_14 + IT_0040;
    const ccomplex_t IT_0042 = m_d*IT_0041;
    const ccomplex_t IT_0043 = (-4)*IT_0042;
    const ccomplex_t IT_0044 = m_d*IT_0016;
    const ccomplex_t IT_0045 = (-8)*conj(IT_0011);
    const ccomplex_t IT_0046 = 0.125*IT_0045;
    const ccomplex_t IT_0047 = IT_0011 + IT_0046;
    const ccomplex_t IT_0048 = m_d*IT_0028;
    const ccomplex_t IT_0049 = 8*IT_0025;
    const ccomplex_t IT_0050 = (-8)*conj(IT_0025);
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = pow(m_W, 4);
    const ccomplex_t IT_0053 = s_13*IT_0052;
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = s_13*s_14;
    const ccomplex_t IT_0056 = s_34*IT_0006;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0055 + IT_0057;
    const ccomplex_t IT_0059 = s_14*IT_0005;
    const ccomplex_t IT_0060 = s_13*s_34;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0059 + IT_0061;
    const ccomplex_t IT_0063 = pow(s_34, 2);
    const ccomplex_t IT_0064 = IT_0012*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0005 + IT_0065;
    const ccomplex_t IT_0067 = s_13*IT_0066;
    const ccomplex_t IT_0068 = (-8)*IT_0011*(conj(IT_0011)*IT_0017 + conj
      (IT_0025)*IT_0029) + -IT_0031*(IT_0031*(16*s_13 + (-8)*IT_0037) + 16
      *IT_0024*IT_0043 + (-8)*IT_0044*IT_0047 + -IT_0048*IT_0051) + -IT_0024*
      (IT_0024*((-16)*IT_0033 + 32*IT_0034 + 8*IT_0054) + (-8)*IT_0047*IT_0058 +
       -IT_0051*IT_0062) + (-8)*IT_0025*((-0.125)*IT_0029*IT_0045 + conj(IT_0025
      )*IT_0067);
    return create_ccomplex_return(IT_0068);
}

