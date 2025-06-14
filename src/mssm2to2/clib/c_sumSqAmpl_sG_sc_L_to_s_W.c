#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sc_L_to_s_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sc_L_to_s_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = pow(m_s, 2);
    const ccomplex_t IT_0006 = pow(m_sG, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
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
    const ccomplex_t IT_0018 = pow(m_c, 2);
    const ccomplex_t IT_0019 = pow(m_sc_L, 2);
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
    const ccomplex_t IT_0030 = m_s*IT_0024;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = m_s*IT_0016;
    const ccomplex_t IT_0033 = s_13*s_14;
    const ccomplex_t IT_0034 = s_34*IT_0006;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0033 + IT_0035;
    const ccomplex_t IT_0037 = (-8)*IT_0011;
    const ccomplex_t IT_0038 = pow(s_34, 2);
    const ccomplex_t IT_0039 = IT_0012*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0005 + IT_0040;
    const ccomplex_t IT_0042 = s_13*IT_0041;
    const ccomplex_t IT_0043 = m_s*IT_0028;
    const ccomplex_t IT_0044 = s_14*IT_0005;
    const ccomplex_t IT_0045 = s_13*s_34;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0044 + IT_0046;
    const ccomplex_t IT_0048 = (-8)*conj(IT_0011);
    const ccomplex_t IT_0049 = pow(m_W, 2);
    const ccomplex_t IT_0050 = s_13*IT_0049;
    const ccomplex_t IT_0051 = s_14*s_34;
    const ccomplex_t IT_0052 = (-2)*IT_0051;
    const ccomplex_t IT_0053 = IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = s_14*IT_0049;
    const ccomplex_t IT_0056 = IT_0012*IT_0055;
    const ccomplex_t IT_0057 = (-0.25)*IT_0056;
    const ccomplex_t IT_0058 = s_14 + IT_0057;
    const ccomplex_t IT_0059 = m_s*IT_0058;
    const ccomplex_t IT_0060 = (-4)*IT_0059;
    const ccomplex_t IT_0061 = pow(m_W, 4);
    const ccomplex_t IT_0062 = s_13*IT_0061;
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = (conj(IT_0011)*IT_0017 + conj(IT_0025)*IT_0029 
      + conj(IT_0031)*IT_0032 + conj(IT_0024)*IT_0036)*IT_0037 + (-8)*IT_0025*
      (conj(IT_0025)*IT_0042 + conj(IT_0031)*IT_0043 + conj(IT_0024)*IT_0047 + (
      -0.125)*IT_0029*IT_0048) + (-8)*IT_0031*(conj(IT_0025)*IT_0043 + (-0.125)
      *IT_0032*IT_0048 + (-0.125)*conj(IT_0031)*(16*s_13 + (-8)*IT_0054) + -conj
      (IT_0024)*IT_0060) + (-8)*IT_0024*(conj(IT_0025)*IT_0047 + (-0.125)
      *IT_0036*IT_0048 + -conj(IT_0031)*IT_0060 + (-0.125)*conj(IT_0024)*((-16)
      *IT_0050 + 32*IT_0051 + 8*IT_0063));
    return create_ccomplex_return(IT_0064);
}

