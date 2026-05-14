#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_su_L_to_u_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_su_L_to_u_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = conj(N_B1)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_W1)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + 3*IT_0011);
    const ccomplex_t IT_0013 = 0.166666666666667*IT_0012;
    const ccomplex_t IT_0014 = pow(m_u, 2);
    const ccomplex_t IT_0015 = pow(m_N_1, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0014 + IT_0015 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = pow(m_su_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0020 
      + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0000*IT_0013;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0019 + 2*IT_0024;
    const ccomplex_t IT_0026 = pow(m_W, -1);
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0008*IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = IT_0016*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = IT_0000*IT_0030;
    const ccomplex_t IT_0035 = IT_0021*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0033 + (-2)*IT_0036;
    const ccomplex_t IT_0038 = (-8)*IT_0014;
    const ccomplex_t IT_0039 = -IT_0000;
    const ccomplex_t IT_0040 = IT_0018*IT_0039;
    const ccomplex_t IT_0041 = IT_0032*IT_0039;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (-8)*s_13;
    const ccomplex_t IT_0044 = s_34*IT_0043;
    const ccomplex_t IT_0045 = -IT_0019;
    const ccomplex_t IT_0046 = -IT_0036;
    const ccomplex_t IT_0047 = m_u*s_14*m_N_1;
    const ccomplex_t IT_0048 = s_34*IT_0015;
    const ccomplex_t IT_0049 = 8*IT_0048;
    const ccomplex_t IT_0050 = s_13*s_14;
    const ccomplex_t IT_0051 = (-16)*IT_0050;
    const ccomplex_t IT_0052 = IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = s_14*IT_0043;
    const ccomplex_t IT_0054 = IT_0015*IT_0043;
    const ccomplex_t IT_0055 = pow(s_13, 2);
    const ccomplex_t IT_0056 = (-8)*IT_0055;
    const ccomplex_t IT_0057 = m_u*m_N_1;
    const ccomplex_t IT_0058 = (-8)*s_14;
    const ccomplex_t IT_0059 = (-8)*s_34;
    const ccomplex_t IT_0060 = (-8)*IT_0015;
    const ccomplex_t IT_0061 = (conj(IT_0025)*IT_0037 + IT_0025*conj(IT_0037))
      *IT_0038 + IT_0043*(conj(IT_0025)*IT_0033 + IT_0025*conj(IT_0033) + conj
      (IT_0037)*IT_0045 + IT_0037*conj(IT_0045)) + (conj(IT_0033)*IT_0040 +
       IT_0033*conj(IT_0040) + conj(IT_0042)*IT_0045 + IT_0042*conj(IT_0045))
      *IT_0058 + (conj(IT_0037)*IT_0040 + IT_0037*conj(IT_0040) + conj(IT_0025)
      *IT_0042 + IT_0025*conj(IT_0042))*IT_0059 + (conj(IT_0033)*IT_0045 +
       IT_0033*conj(IT_0045))*IT_0060;
    const ccomplex_t IT_0062 = m_u*s_34*m_N_1;
    const ccomplex_t IT_0063 = (-8)*conj(IT_0025);
    const ccomplex_t IT_0064 = s_14*IT_0014;
    const ccomplex_t IT_0065 = (-8)*conj(IT_0037);
    const ccomplex_t IT_0066 = s_14*s_34;
    const ccomplex_t IT_0067 = (-8)*IT_0037;
    const ccomplex_t IT_0068 = (-0.125)*IT_0064;
    const ccomplex_t IT_0069 = (-8)*IT_0025;
    const ccomplex_t IT_0070 = (-0.125)*IT_0062;
    const ccomplex_t IT_0071 = s_13*(IT_0025*conj(IT_0025) + IT_0037*conj
      (IT_0037))*IT_0038 + (conj(IT_0025)*IT_0040 + IT_0025*conj(IT_0040) + conj
      (IT_0037)*IT_0042 + IT_0037*conj(IT_0042))*IT_0044 + (-8)*conj(IT_0045)*
      (IT_0046*IT_0047 + (-0.125)*IT_0024*IT_0052 + (-0.125)*IT_0040*IT_0053 + (
      -0.125)*IT_0045*IT_0054 + (-0.125)*IT_0025*IT_0056) + IT_0045*(conj
      (IT_0040)*IT_0053 + conj(IT_0025)*IT_0056) + (-8)*conj(IT_0033)*(IT_0024
      *IT_0047 + (-0.125)*IT_0046*IT_0052 + (-0.125)*IT_0042*IT_0053 + (-0.125)
      *IT_0033*IT_0054 + (-0.125)*IT_0037*IT_0056) + IT_0033*(conj(IT_0042)
      *IT_0053 + conj(IT_0037)*IT_0056) + IT_0057*IT_0061 + IT_0024*(IT_0063
      *IT_0064 + IT_0062*IT_0065) + IT_0046*(IT_0062*IT_0063 + IT_0064*IT_0065) 
      + (-8)*conj(IT_0024)*(IT_0033*IT_0047 + (-0.125)*IT_0045*IT_0052 + (-4)
      *IT_0024*IT_0066 + IT_0068*IT_0069 + IT_0067*IT_0070) + (-8)*conj(IT_0046)
      *(IT_0045*IT_0047 + (-0.125)*IT_0033*IT_0052 + (-4)*IT_0046*IT_0066 +
       IT_0067*IT_0068 + IT_0069*IT_0070);
    return create_ccomplex_return(IT_0071);
}

