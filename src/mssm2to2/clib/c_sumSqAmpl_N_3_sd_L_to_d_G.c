#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sd_L_to_d_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sd_L_to_d_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = pow(m_d, 2);
    const ccomplex_t IT_0010 = pow(m_N_3, 2);
    const ccomplex_t IT_0011 = pow(m_sd_L, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + IT_0009 + -IT_0010 + -IT_0011 
      + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (-2)*IT_0000;
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (-2)*IT_0014 + -IT_0019;
    const ccomplex_t IT_0021 = cos(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = conj(N_B3)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W3)*e_em;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-3)*IT_0028);
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = IT_0000*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0016*IT_0030;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0015*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = (-8)*IT_0009;
    const ccomplex_t IT_0039 = -IT_0000;
    const ccomplex_t IT_0040 = IT_0018*IT_0039;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0035*IT_0039;
    const ccomplex_t IT_0043 = (-8)*s_13;
    const ccomplex_t IT_0044 = s_34*IT_0043;
    const ccomplex_t IT_0045 = m_d*s_14*m_N_3;
    const ccomplex_t IT_0046 = -IT_0014;
    const ccomplex_t IT_0047 = s_34*IT_0010;
    const ccomplex_t IT_0048 = 8*IT_0047;
    const ccomplex_t IT_0049 = s_13*s_14;
    const ccomplex_t IT_0050 = (-16)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = s_14*IT_0043;
    const ccomplex_t IT_0053 = IT_0010*IT_0043;
    const ccomplex_t IT_0054 = pow(s_13, 2);
    const ccomplex_t IT_0055 = (-8)*IT_0054;
    const ccomplex_t IT_0056 = -IT_0036;
    const ccomplex_t IT_0057 = m_d*m_N_3;
    const ccomplex_t IT_0058 = (-8)*s_14;
    const ccomplex_t IT_0059 = (-8)*s_34;
    const ccomplex_t IT_0060 = (-8)*IT_0010;
    const ccomplex_t IT_0061 = (conj(IT_0020)*IT_0037 + IT_0020*conj(IT_0037))
      *IT_0038 + IT_0043*(conj(IT_0019)*IT_0037 + IT_0019*conj(IT_0037) + conj
      (IT_0020)*IT_0056 + IT_0020*conj(IT_0056)) + (conj(IT_0019)*IT_0042 +
       IT_0019*conj(IT_0042) + conj(IT_0041)*IT_0056 + IT_0041*conj(IT_0056))
      *IT_0058 + (conj(IT_0037)*IT_0041 + IT_0037*conj(IT_0041) + conj(IT_0020)
      *IT_0042 + IT_0020*conj(IT_0042))*IT_0059 + (conj(IT_0019)*IT_0056 +
       IT_0019*conj(IT_0056))*IT_0060;
    const ccomplex_t IT_0062 = m_d*s_34*m_N_3;
    const ccomplex_t IT_0063 = (-8)*conj(IT_0020);
    const ccomplex_t IT_0064 = s_14*IT_0009;
    const ccomplex_t IT_0065 = (-8)*conj(IT_0037);
    const ccomplex_t IT_0066 = s_14*s_34;
    const ccomplex_t IT_0067 = (-8)*IT_0037;
    const ccomplex_t IT_0068 = (-0.125)*IT_0064;
    const ccomplex_t IT_0069 = (-8)*IT_0020;
    const ccomplex_t IT_0070 = (-0.125)*IT_0062;
    const ccomplex_t IT_0071 = s_13*(IT_0020*conj(IT_0020) + IT_0037*conj
      (IT_0037))*IT_0038 + (conj(IT_0020)*IT_0041 + IT_0020*conj(IT_0041) + conj
      (IT_0037)*IT_0042 + IT_0037*conj(IT_0042))*IT_0044 + (-8)*conj(IT_0019)*
      (IT_0033*IT_0045 + (-0.125)*IT_0046*IT_0051 + (-0.125)*IT_0041*IT_0052 + (
      -0.125)*IT_0019*IT_0053 + (-0.125)*IT_0020*IT_0055) + IT_0019*(conj
      (IT_0041)*IT_0052 + conj(IT_0020)*IT_0055) + (conj(IT_0042)*IT_0052 + conj
      (IT_0037)*IT_0055)*IT_0056 + (-8)*(IT_0045*IT_0046 + (-0.125)*IT_0033
      *IT_0051 + (-0.125)*IT_0042*IT_0052 + (-0.125)*IT_0037*IT_0055 + (-0.125)
      *IT_0053*IT_0056)*conj(IT_0056) + IT_0057*IT_0061 + IT_0046*(IT_0063
      *IT_0064 + IT_0062*IT_0065) + IT_0033*(IT_0062*IT_0063 + IT_0064*IT_0065) 
      + (-8)*conj(IT_0046)*((-0.125)*IT_0019*IT_0051 + IT_0045*IT_0056 + (-4)
      *IT_0046*IT_0066 + IT_0068*IT_0069 + IT_0067*IT_0070) + (-8)*conj(IT_0033)
      *(IT_0019*IT_0045 + (-0.125)*IT_0051*IT_0056 + (-4)*IT_0033*IT_0066 +
       IT_0067*IT_0068 + IT_0069*IT_0070);
    return create_ccomplex_return(IT_0071);
}

