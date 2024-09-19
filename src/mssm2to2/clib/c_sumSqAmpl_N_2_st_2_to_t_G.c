#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_2_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_2_to_t_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 1.4142135623731*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = m_t*N_u2*e_em*IT_0011*U_st_01;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0005 + (-0.75)*IT_0014);
    const ccomplex_t IT_0016 = 0.666666666666667*IT_0015;
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = pow(m_N_2, 2);
    const ccomplex_t IT_0019 = pow(m_st_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (-2)*IT_0000;
    const ccomplex_t IT_0024 = pow(m_t, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0018 + IT_0024 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0023*IT_0027;
    const ccomplex_t IT_0029 = (-2)*IT_0022 + -IT_0028;
    const ccomplex_t IT_0030 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0031 = IT_0002*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = m_t*conj(N_u2)*e_em*IT_0011*U_st_11;
    const ccomplex_t IT_0037 = IT_0010*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0032 + 3*IT_0035 + 3
      *IT_0038);
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = IT_0000*IT_0040;
    const ccomplex_t IT_0042 = IT_0020*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0025*IT_0040;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0023*IT_0045;
    const ccomplex_t IT_0047 = 2*IT_0043 + IT_0046;
    const ccomplex_t IT_0048 = (-8)*IT_0024;
    const ccomplex_t IT_0049 = -IT_0000;
    const ccomplex_t IT_0050 = IT_0027*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0045*IT_0049;
    const ccomplex_t IT_0053 = (-8)*s_13;
    const ccomplex_t IT_0054 = s_34*IT_0053;
    const ccomplex_t IT_0055 = m_t*s_14*m_N_2;
    const ccomplex_t IT_0056 = -IT_0022;
    const ccomplex_t IT_0057 = s_34*IT_0018;
    const ccomplex_t IT_0058 = 8*IT_0057;
    const ccomplex_t IT_0059 = s_13*s_14;
    const ccomplex_t IT_0060 = (-16)*IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = s_14*IT_0053;
    const ccomplex_t IT_0063 = IT_0018*IT_0053;
    const ccomplex_t IT_0064 = pow(s_13, 2);
    const ccomplex_t IT_0065 = (-8)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0046;
    const ccomplex_t IT_0067 = m_t*m_N_2;
    const ccomplex_t IT_0068 = (-8)*s_14;
    const ccomplex_t IT_0069 = (-8)*s_34;
    const ccomplex_t IT_0070 = (-8)*IT_0018;
    const ccomplex_t IT_0071 = (conj(IT_0029)*IT_0047 + IT_0029*conj(IT_0047))
      *IT_0048 + IT_0053*(conj(IT_0028)*IT_0047 + IT_0028*conj(IT_0047) + conj
      (IT_0029)*IT_0066 + IT_0029*conj(IT_0066)) + (conj(IT_0028)*IT_0052 +
       IT_0028*conj(IT_0052) + conj(IT_0051)*IT_0066 + IT_0051*conj(IT_0066))
      *IT_0068 + (conj(IT_0047)*IT_0051 + IT_0047*conj(IT_0051) + conj(IT_0029)
      *IT_0052 + IT_0029*conj(IT_0052))*IT_0069 + (conj(IT_0028)*IT_0066 +
       IT_0028*conj(IT_0066))*IT_0070;
    const ccomplex_t IT_0072 = m_t*s_34*m_N_2;
    const ccomplex_t IT_0073 = (-8)*conj(IT_0029);
    const ccomplex_t IT_0074 = s_14*IT_0024;
    const ccomplex_t IT_0075 = (-8)*conj(IT_0047);
    const ccomplex_t IT_0076 = s_14*s_34;
    const ccomplex_t IT_0077 = (-8)*IT_0047;
    const ccomplex_t IT_0078 = (-0.125)*IT_0074;
    const ccomplex_t IT_0079 = (-8)*IT_0029;
    const ccomplex_t IT_0080 = (-0.125)*IT_0072;
    const ccomplex_t IT_0081 = s_13*(IT_0029*conj(IT_0029) + IT_0047*conj
      (IT_0047))*IT_0048 + (conj(IT_0029)*IT_0051 + IT_0029*conj(IT_0051) + conj
      (IT_0047)*IT_0052 + IT_0047*conj(IT_0052))*IT_0054 + (-8)*conj(IT_0028)*
      (IT_0043*IT_0055 + (-0.125)*IT_0056*IT_0061 + (-0.125)*IT_0051*IT_0062 + (
      -0.125)*IT_0028*IT_0063 + (-0.125)*IT_0029*IT_0065) + IT_0028*(conj
      (IT_0051)*IT_0062 + conj(IT_0029)*IT_0065) + (conj(IT_0052)*IT_0062 + conj
      (IT_0047)*IT_0065)*IT_0066 + (-8)*(IT_0055*IT_0056 + (-0.125)*IT_0043
      *IT_0061 + (-0.125)*IT_0052*IT_0062 + (-0.125)*IT_0047*IT_0065 + (-0.125)
      *IT_0063*IT_0066)*conj(IT_0066) + IT_0067*IT_0071 + IT_0056*(IT_0073
      *IT_0074 + IT_0072*IT_0075) + IT_0043*(IT_0072*IT_0073 + IT_0074*IT_0075) 
      + (-8)*conj(IT_0056)*((-0.125)*IT_0028*IT_0061 + IT_0055*IT_0066 + (-4)
      *IT_0056*IT_0076 + IT_0078*IT_0079 + IT_0077*IT_0080) + (-8)*conj(IT_0043)
      *(IT_0028*IT_0055 + (-0.125)*IT_0061*IT_0066 + (-4)*IT_0043*IT_0076 +
       IT_0077*IT_0078 + IT_0079*IT_0080);
    return create_ccomplex_return(IT_0081);
}

