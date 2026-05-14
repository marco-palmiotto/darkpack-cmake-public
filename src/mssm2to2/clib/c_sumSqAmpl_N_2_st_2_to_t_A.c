#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_2_to_t_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_2_to_t_A(
        cparam_t const *param
        )
{
    clearcache();
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_t*N_u2*e_em*IT_0012*U_st_01;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0006 + (-0.75)*IT_0015);
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0016;
    const ccomplex_t IT_0018 = pow(m_t, 2);
    const ccomplex_t IT_0019 = pow(m_N_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0025 = IT_0017*IT_0024;
    const ccomplex_t IT_0026 = pow(m_st_2, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + -IT_0019 + -IT_0026 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0023 + (-2)*IT_0029;
    const ccomplex_t IT_0031 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0032 = IT_0003*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = m_t*conj(N_u2)*e_em*IT_0012*U_st_11;
    const ccomplex_t IT_0038 = IT_0011*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + 3*IT_0036 + 3
      *IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0024*IT_0041;
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0020*IT_0041;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 2*IT_0044 + IT_0047;
    const ccomplex_t IT_0049 = (-6)*IT_0018;
    const ccomplex_t IT_0050 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0051 = IT_0022*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0046*IT_0050;
    const ccomplex_t IT_0054 = (-6)*s_13;
    const ccomplex_t IT_0055 = s_34*IT_0054;
    const ccomplex_t IT_0056 = m_t*s_14*m_N_2;
    const ccomplex_t IT_0057 = -IT_0029;
    const ccomplex_t IT_0058 = s_34*IT_0019;
    const ccomplex_t IT_0059 = 6*IT_0058;
    const ccomplex_t IT_0060 = s_13*s_14;
    const ccomplex_t IT_0061 = (-12)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059 + IT_0061;
    const ccomplex_t IT_0063 = s_14*IT_0054;
    const ccomplex_t IT_0064 = IT_0019*IT_0054;
    const ccomplex_t IT_0065 = pow(s_13, 2);
    const ccomplex_t IT_0066 = (-6)*IT_0065;
    const ccomplex_t IT_0067 = -IT_0047;
    const ccomplex_t IT_0068 = m_t*m_N_2;
    const ccomplex_t IT_0069 = (-6)*s_14;
    const ccomplex_t IT_0070 = (-6)*s_34;
    const ccomplex_t IT_0071 = (-6)*IT_0019;
    const ccomplex_t IT_0072 = (conj(IT_0030)*IT_0048 + IT_0030*conj(IT_0048))
      *IT_0049 + IT_0054*(conj(IT_0023)*IT_0048 + IT_0023*conj(IT_0048) + conj
      (IT_0030)*IT_0067 + IT_0030*conj(IT_0067)) + (conj(IT_0023)*IT_0053 +
       IT_0023*conj(IT_0053) + conj(IT_0052)*IT_0067 + IT_0052*conj(IT_0067))
      *IT_0069 + (conj(IT_0048)*IT_0052 + IT_0048*conj(IT_0052) + conj(IT_0030)
      *IT_0053 + IT_0030*conj(IT_0053))*IT_0070 + (conj(IT_0023)*IT_0067 +
       IT_0023*conj(IT_0067))*IT_0071;
    const ccomplex_t IT_0073 = m_t*s_34*m_N_2;
    const ccomplex_t IT_0074 = (-6)*conj(IT_0030);
    const ccomplex_t IT_0075 = s_14*IT_0018;
    const ccomplex_t IT_0076 = (-6)*conj(IT_0048);
    const ccomplex_t IT_0077 = s_14*s_34;
    const ccomplex_t IT_0078 = (-6)*IT_0048;
    const ccomplex_t IT_0079 = (-0.166666666666667)*IT_0075;
    const ccomplex_t IT_0080 = (-6)*IT_0030;
    const ccomplex_t IT_0081 = (-0.166666666666667)*IT_0073;
    const ccomplex_t IT_0082 = s_13*(IT_0030*conj(IT_0030) + IT_0048*conj
      (IT_0048))*IT_0049 + (conj(IT_0030)*IT_0052 + IT_0030*conj(IT_0052) + conj
      (IT_0048)*IT_0053 + IT_0048*conj(IT_0053))*IT_0055 + (-6)*conj(IT_0023)*
      (IT_0044*IT_0056 + (-0.166666666666667)*IT_0057*IT_0062 + (
      -0.166666666666667)*IT_0052*IT_0063 + (-0.166666666666667)*IT_0023*IT_0064
       + (-0.166666666666667)*IT_0030*IT_0066) + IT_0023*(conj(IT_0052)*IT_0063 
      + conj(IT_0030)*IT_0066) + (conj(IT_0053)*IT_0063 + conj(IT_0048)*IT_0066)
      *IT_0067 + (-6)*(IT_0056*IT_0057 + (-0.166666666666667)*IT_0044*IT_0062 + 
      (-0.166666666666667)*IT_0053*IT_0063 + (-0.166666666666667)*IT_0048
      *IT_0066 + (-0.166666666666667)*IT_0064*IT_0067)*conj(IT_0067) + IT_0068
      *IT_0072 + IT_0057*(IT_0074*IT_0075 + IT_0073*IT_0076) + IT_0044*(IT_0073
      *IT_0074 + IT_0075*IT_0076) + (-6)*conj(IT_0057)*((-0.166666666666667)
      *IT_0023*IT_0062 + IT_0056*IT_0067 + (-4)*IT_0057*IT_0077 + IT_0079
      *IT_0080 + IT_0078*IT_0081) + (-6)*conj(IT_0044)*(IT_0023*IT_0056 + (
      -0.166666666666667)*IT_0062*IT_0067 + (-4)*IT_0044*IT_0077 + IT_0078
      *IT_0079 + IT_0080*IT_0081);
    return create_ccomplex_return(IT_0082);
}

