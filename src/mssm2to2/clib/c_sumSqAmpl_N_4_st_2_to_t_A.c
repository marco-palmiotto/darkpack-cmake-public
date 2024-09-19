#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_st_2_to_t_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_st_2_to_t_A(
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = m_t*N_u4*e_em*IT_0007*U_st_01;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0010 + (
      -1.33333333333333)*IT_0015);
    const ccomplex_t IT_0017 = (-0.75)*IT_0016;
    const ccomplex_t IT_0018 = 0.666666666666667*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = pow(m_N_4, 2);
    const ccomplex_t IT_0021 = pow(m_st_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0026 = pow(m_t, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0020 + IT_0026 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0025*IT_0029;
    const ccomplex_t IT_0031 = (-2)*IT_0024 + -IT_0030;
    const ccomplex_t IT_0032 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0033 = IT_0012*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0036 = IT_0005*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = m_t*conj(N_u4)*e_em*IT_0007*U_st_11;
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0034 + 3*IT_0037 + 3
      *IT_0040);
    const ccomplex_t IT_0042 = 0.166666666666667*IT_0041;
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0025*IT_0044;
    const ccomplex_t IT_0046 = IT_0001*IT_0042;
    const ccomplex_t IT_0047 = IT_0022*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0045 + 2*IT_0048;
    const ccomplex_t IT_0050 = (-6)*IT_0026;
    const ccomplex_t IT_0051 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0044*IT_0051;
    const ccomplex_t IT_0055 = (-6)*s_13;
    const ccomplex_t IT_0056 = s_34*IT_0055;
    const ccomplex_t IT_0057 = m_t*s_14*m_N_4;
    const ccomplex_t IT_0058 = -IT_0024;
    const ccomplex_t IT_0059 = s_34*IT_0020;
    const ccomplex_t IT_0060 = 6*IT_0059;
    const ccomplex_t IT_0061 = s_13*s_14;
    const ccomplex_t IT_0062 = (-12)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = s_14*IT_0055;
    const ccomplex_t IT_0065 = IT_0020*IT_0055;
    const ccomplex_t IT_0066 = pow(s_13, 2);
    const ccomplex_t IT_0067 = (-6)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0045;
    const ccomplex_t IT_0069 = m_t*m_N_4;
    const ccomplex_t IT_0070 = (-6)*s_14;
    const ccomplex_t IT_0071 = (-6)*s_34;
    const ccomplex_t IT_0072 = (-6)*IT_0020;
    const ccomplex_t IT_0073 = (conj(IT_0031)*IT_0049 + IT_0031*conj(IT_0049))
      *IT_0050 + IT_0055*(conj(IT_0030)*IT_0049 + IT_0030*conj(IT_0049) + conj
      (IT_0031)*IT_0068 + IT_0031*conj(IT_0068)) + (conj(IT_0030)*IT_0054 +
       IT_0030*conj(IT_0054) + conj(IT_0053)*IT_0068 + IT_0053*conj(IT_0068))
      *IT_0070 + (conj(IT_0049)*IT_0053 + IT_0049*conj(IT_0053) + conj(IT_0031)
      *IT_0054 + IT_0031*conj(IT_0054))*IT_0071 + (conj(IT_0030)*IT_0068 +
       IT_0030*conj(IT_0068))*IT_0072;
    const ccomplex_t IT_0074 = m_t*s_34*m_N_4;
    const ccomplex_t IT_0075 = (-6)*conj(IT_0031);
    const ccomplex_t IT_0076 = s_14*IT_0026;
    const ccomplex_t IT_0077 = (-6)*conj(IT_0049);
    const ccomplex_t IT_0078 = s_14*s_34;
    const ccomplex_t IT_0079 = (-6)*IT_0049;
    const ccomplex_t IT_0080 = (-0.166666666666667)*IT_0076;
    const ccomplex_t IT_0081 = (-6)*IT_0031;
    const ccomplex_t IT_0082 = (-0.166666666666667)*IT_0074;
    const ccomplex_t IT_0083 = s_13*(IT_0031*conj(IT_0031) + IT_0049*conj
      (IT_0049))*IT_0050 + (conj(IT_0031)*IT_0053 + IT_0031*conj(IT_0053) + conj
      (IT_0049)*IT_0054 + IT_0049*conj(IT_0054))*IT_0056 + (-6)*conj(IT_0030)*
      (IT_0048*IT_0057 + (-0.166666666666667)*IT_0058*IT_0063 + (
      -0.166666666666667)*IT_0053*IT_0064 + (-0.166666666666667)*IT_0030*IT_0065
       + (-0.166666666666667)*IT_0031*IT_0067) + IT_0030*(conj(IT_0053)*IT_0064 
      + conj(IT_0031)*IT_0067) + (conj(IT_0054)*IT_0064 + conj(IT_0049)*IT_0067)
      *IT_0068 + (-6)*(IT_0057*IT_0058 + (-0.166666666666667)*IT_0048*IT_0063 + 
      (-0.166666666666667)*IT_0054*IT_0064 + (-0.166666666666667)*IT_0049
      *IT_0067 + (-0.166666666666667)*IT_0065*IT_0068)*conj(IT_0068) + IT_0069
      *IT_0073 + IT_0058*(IT_0075*IT_0076 + IT_0074*IT_0077) + IT_0048*(IT_0074
      *IT_0075 + IT_0076*IT_0077) + (-6)*conj(IT_0058)*((-0.166666666666667)
      *IT_0030*IT_0063 + IT_0057*IT_0068 + (-4)*IT_0058*IT_0078 + IT_0080
      *IT_0081 + IT_0079*IT_0082) + (-6)*conj(IT_0048)*(IT_0030*IT_0057 + (
      -0.166666666666667)*IT_0063*IT_0068 + (-4)*IT_0048*IT_0078 + IT_0079
      *IT_0080 + IT_0081*IT_0082);
    return create_ccomplex_return(IT_0083);
}

