#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_1_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_1_to_t_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = m_b*conj(U_d2)*V_tb*e_em*IT_0009*U_sb_10;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0012 + (-2)*IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = s_14*s_34;
    const ccomplex_t IT_0021 = (-2)*IT_0003;
    const ccomplex_t IT_0022 = pow(m_t, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0000 + IT_0022 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = s_34*IT_0000;
    const ccomplex_t IT_0029 = 8*IT_0028;
    const ccomplex_t IT_0030 = s_13*s_14;
    const ccomplex_t IT_0031 = (-16)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029 + IT_0031;
    const ccomplex_t IT_0033 = -IT_0003;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = (-8)*s_13;
    const ccomplex_t IT_0036 = s_14*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0019 + IT_0026;
    const ccomplex_t IT_0038 = s_34*IT_0035;
    const ccomplex_t IT_0039 = sin(beta);
    const ccomplex_t IT_0040 = cpow(IT_0039, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0007*IT_0009*IT_0040*U_sb_00;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0023*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0021*IT_0044;
    const ccomplex_t IT_0046 = m_t*m_C_2;
    const ccomplex_t IT_0047 = (-8)*s_14*IT_0046;
    const ccomplex_t IT_0048 = IT_0003*IT_0042;
    const ccomplex_t IT_0049 = IT_0002*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0045 + (-2)*IT_0050;
    const ccomplex_t IT_0052 = (-8)*s_34*IT_0046;
    const ccomplex_t IT_0053 = IT_0033*IT_0044;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0056 = m_t*s_34*m_C_2;
    const ccomplex_t IT_0057 = s_14*IT_0022;
    const ccomplex_t IT_0058 = (-8)*IT_0019;
    const ccomplex_t IT_0059 = -IT_0050;
    const ccomplex_t IT_0060 = IT_0000*IT_0035;
    const ccomplex_t IT_0061 = IT_0035*IT_0046;
    const ccomplex_t IT_0062 = pow(s_13, 2);
    const ccomplex_t IT_0063 = (-8)*IT_0062;
    const ccomplex_t IT_0064 = (-8)*IT_0000*IT_0046;
    const ccomplex_t IT_0065 = (-8)*conj(IT_0019);
    const ccomplex_t IT_0066 = (-8)*IT_0055;
    const ccomplex_t IT_0067 = (-8)*IT_0022;
    const ccomplex_t IT_0068 = IT_0046*IT_0067;
    const ccomplex_t IT_0069 = (-8)*IT_0057;
    const ccomplex_t IT_0070 = (-8)*IT_0056;
    const ccomplex_t IT_0071 = IT_0019*(32*conj(IT_0019)*IT_0020 + conj
      (IT_0027)*IT_0032) + IT_0034*(conj(IT_0027)*IT_0036 + conj(IT_0037)
      *IT_0038 + conj(IT_0045)*IT_0047 + conj(IT_0051)*IT_0052) + (IT_0036*conj
      (IT_0045) + conj(IT_0027)*IT_0047 + IT_0038*conj(IT_0051) + conj(IT_0037)
      *IT_0052)*IT_0054 + (conj(IT_0045)*IT_0055 + conj(IT_0051)*IT_0056 + conj
      (IT_0037)*IT_0057)*IT_0058 + IT_0045*(conj(IT_0034)*IT_0047 + IT_0036*conj
      (IT_0054) + IT_0032*conj(IT_0059) + conj(IT_0045)*IT_0060 + conj(IT_0037)
      *IT_0061 + conj(IT_0051)*IT_0063 + conj(IT_0027)*IT_0064 + IT_0055*IT_0065
      ) + IT_0027*(conj(IT_0019)*IT_0032 + conj(IT_0034)*IT_0036 + IT_0047*conj
      (IT_0054) + conj(IT_0027)*IT_0060 + conj(IT_0051)*IT_0061 + conj(IT_0037)
      *IT_0063 + conj(IT_0045)*IT_0064 + conj(IT_0059)*IT_0066) + IT_0051*(conj
      (IT_0034)*IT_0052 + IT_0038*conj(IT_0054) + conj(IT_0027)*IT_0061 + conj
      (IT_0045)*IT_0063 + IT_0056*IT_0065 + s_13*conj(IT_0051)*IT_0067 + conj
      (IT_0037)*IT_0068 + conj(IT_0059)*IT_0069) + IT_0059*(IT_0032*conj(IT_0045
      ) + 32*IT_0020*conj(IT_0059) + conj(IT_0027)*IT_0066 + conj(IT_0051)
      *IT_0069 + conj(IT_0037)*IT_0070) + IT_0037*(conj(IT_0034)*IT_0038 +
       IT_0052*conj(IT_0054) + conj(IT_0045)*IT_0061 + conj(IT_0027)*IT_0063 +
       IT_0057*IT_0065 + s_13*conj(IT_0037)*IT_0067 + conj(IT_0051)*IT_0068 +
       conj(IT_0059)*IT_0070);
    return create_ccomplex_return(IT_0071);
}

