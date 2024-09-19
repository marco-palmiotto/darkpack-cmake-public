#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_t_G(
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = m_b*conj(U_d2)*V_tb*e_em*IT_0008*U_sb_11;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0004 + (-0.5)*IT_0011);
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = pow(m_C_2, 2);
    const ccomplex_t IT_0015 = pow(m_sb_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + -IT_0014 + -IT_0015 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = s_14*s_34;
    const ccomplex_t IT_0020 = (-2)*IT_0000;
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0014 + IT_0021 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0020*IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = s_34*IT_0014;
    const ccomplex_t IT_0028 = 8*IT_0027;
    const ccomplex_t IT_0029 = s_13*s_14;
    const ccomplex_t IT_0030 = (-16)*IT_0029;
    const ccomplex_t IT_0031 = IT_0028 + IT_0030;
    const ccomplex_t IT_0032 = -IT_0000;
    const ccomplex_t IT_0033 = IT_0024*IT_0032;
    const ccomplex_t IT_0034 = (-8)*s_13;
    const ccomplex_t IT_0035 = s_14*IT_0034;
    const ccomplex_t IT_0036 = 2*IT_0018 + IT_0025;
    const ccomplex_t IT_0037 = s_34*IT_0034;
    const ccomplex_t IT_0038 = sin(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0002*IT_0008*IT_0039*U_sb_01;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0022*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = m_t*m_C_2;
    const ccomplex_t IT_0046 = (-8)*s_14*IT_0045;
    const ccomplex_t IT_0047 = IT_0000*IT_0041;
    const ccomplex_t IT_0048 = IT_0016*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0044 + (-2)*IT_0049;
    const ccomplex_t IT_0051 = (-8)*s_34*IT_0045;
    const ccomplex_t IT_0052 = IT_0032*IT_0043;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0014*IT_0034;
    const ccomplex_t IT_0055 = IT_0034*IT_0045;
    const ccomplex_t IT_0056 = pow(s_13, 2);
    const ccomplex_t IT_0057 = (-8)*IT_0056;
    const ccomplex_t IT_0058 = (-8)*IT_0014*IT_0045;
    const ccomplex_t IT_0059 = -IT_0049;
    const ccomplex_t IT_0060 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0061 = (-8)*IT_0060;
    const ccomplex_t IT_0062 = s_14*IT_0021;
    const ccomplex_t IT_0063 = (-8)*IT_0062;
    const ccomplex_t IT_0064 = m_t*s_34*m_C_2;
    const ccomplex_t IT_0065 = (-8)*IT_0064;
    const ccomplex_t IT_0066 = (-8)*IT_0018;
    const ccomplex_t IT_0067 = (-8)*conj(IT_0018);
    const ccomplex_t IT_0068 = (-8)*IT_0021;
    const ccomplex_t IT_0069 = IT_0045*IT_0068;
    const ccomplex_t IT_0070 = IT_0018*(32*conj(IT_0018)*IT_0019 + conj
      (IT_0026)*IT_0031) + IT_0033*(conj(IT_0026)*IT_0035 + conj(IT_0036)
      *IT_0037 + conj(IT_0044)*IT_0046 + conj(IT_0050)*IT_0051) + (IT_0035*conj
      (IT_0044) + conj(IT_0026)*IT_0046 + IT_0037*conj(IT_0050) + conj(IT_0036)
      *IT_0051)*IT_0053 + IT_0026*(conj(IT_0018)*IT_0031 + conj(IT_0033)*IT_0035
       + IT_0046*conj(IT_0053) + conj(IT_0026)*IT_0054 + conj(IT_0050)*IT_0055 +
       conj(IT_0036)*IT_0057 + conj(IT_0044)*IT_0058 + conj(IT_0059)*IT_0061) +
       IT_0059*(IT_0031*conj(IT_0044) + 32*IT_0019*conj(IT_0059) + conj(IT_0026)
      *IT_0061 + conj(IT_0050)*IT_0063 + conj(IT_0036)*IT_0065) + (conj(IT_0044)
      *IT_0060 + conj(IT_0036)*IT_0062 + conj(IT_0050)*IT_0064)*IT_0066 +
       IT_0044*(conj(IT_0033)*IT_0046 + IT_0035*conj(IT_0053) + conj(IT_0044)
      *IT_0054 + conj(IT_0036)*IT_0055 + conj(IT_0050)*IT_0057 + conj(IT_0026)
      *IT_0058 + IT_0031*conj(IT_0059) + IT_0060*IT_0067) + IT_0050*(conj
      (IT_0033)*IT_0051 + IT_0037*conj(IT_0053) + conj(IT_0026)*IT_0055 + conj
      (IT_0044)*IT_0057 + conj(IT_0059)*IT_0063 + IT_0064*IT_0067 + s_13*conj
      (IT_0050)*IT_0068 + conj(IT_0036)*IT_0069) + IT_0036*(conj(IT_0033)
      *IT_0037 + IT_0051*conj(IT_0053) + conj(IT_0044)*IT_0055 + conj(IT_0026)
      *IT_0057 + conj(IT_0059)*IT_0065 + IT_0062*IT_0067 + s_13*conj(IT_0036)
      *IT_0068 + conj(IT_0050)*IT_0069);
    return create_ccomplex_return(IT_0070);
}

