#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_1_to_t_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_1_to_t_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_sb_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0005 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0006*IT_0008*IT_0010*U_sb_00;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = pow(m_t, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0002 + IT_0017 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0004*IT_0012;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0024 = IT_0010*IT_0023;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = m_b*conj(U_d1)*V_tb*e_em*IT_0006*U_sb_10;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0030);
    const ccomplex_t IT_0032 = IT_0005*IT_0031;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0004*IT_0031;
    const ccomplex_t IT_0036 = IT_0018*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0000*(IT_0022*((-24)*conj(IT_0016) + 24*conj
      (IT_0022)) + IT_0016*(24*conj(IT_0016) + (-24)*conj(IT_0022)) + IT_0037*((
      -24)*conj(IT_0034) + 24*conj(IT_0037)) + IT_0034*(24*conj(IT_0034) + (-24)
      *conj(IT_0037)));
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0001 + IT_0017 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0005*IT_0041;
    const ccomplex_t IT_0043 = 2*IT_0021 + IT_0042;
    const ccomplex_t IT_0044 = s_34*IT_0001;
    const ccomplex_t IT_0045 = 6*IT_0044;
    const ccomplex_t IT_0046 = s_13*s_14;
    const ccomplex_t IT_0047 = (-12)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045 + IT_0047;
    const ccomplex_t IT_0049 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0050 = IT_0041*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = (-6)*s_13;
    const ccomplex_t IT_0053 = s_14*IT_0052;
    const ccomplex_t IT_0054 = IT_0031*IT_0039;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0005*IT_0055;
    const ccomplex_t IT_0057 = (-2)*IT_0037 + -IT_0056;
    const ccomplex_t IT_0058 = IT_0049*IT_0055;
    const ccomplex_t IT_0059 = (-2)*IT_0015 + -IT_0042;
    const ccomplex_t IT_0060 = s_14*IT_0017;
    const ccomplex_t IT_0061 = s_13*s_34;
    const ccomplex_t IT_0062 = (-2)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = s_34*IT_0052;
    const ccomplex_t IT_0065 = 2*IT_0034 + IT_0056;
    const ccomplex_t IT_0066 = (-6)*IT_0017;
    const ccomplex_t IT_0067 = pow(s_13, 2);
    const ccomplex_t IT_0068 = (-6)*IT_0067;
    const ccomplex_t IT_0069 = IT_0001*IT_0052;
    const ccomplex_t IT_0070 = m_t*m_C_1;
    const ccomplex_t IT_0071 = (-6)*s_14;
    const ccomplex_t IT_0072 = (-6)*s_34;
    const ccomplex_t IT_0073 = (-6)*IT_0001;
    const ccomplex_t IT_0074 = IT_0052*(conj(IT_0057)*IT_0059 + IT_0057*conj
      (IT_0059) + conj(IT_0043)*IT_0065 + IT_0043*conj(IT_0065)) + (conj(IT_0059
      )*IT_0065 + IT_0059*conj(IT_0065))*IT_0066 + (conj(IT_0051)*IT_0057 +
       IT_0051*conj(IT_0057) + conj(IT_0043)*IT_0058 + IT_0043*conj(IT_0058))
      *IT_0071 + (conj(IT_0058)*IT_0059 + IT_0058*conj(IT_0059) + conj(IT_0051)
      *IT_0065 + IT_0051*conj(IT_0065))*IT_0072 + (conj(IT_0043)*IT_0057 +
       IT_0043*conj(IT_0057))*IT_0073;
    const ccomplex_t IT_0075 = (-6)*IT_0043;
    const ccomplex_t IT_0076 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0077 = (-6)*IT_0057;
    const ccomplex_t IT_0078 = (-6)*conj(IT_0057);
    const ccomplex_t IT_0079 = m_t*s_14*m_C_1;
    const ccomplex_t IT_0080 = conj(IT_0034) + conj(IT_0037);
    const ccomplex_t IT_0081 = IT_0034 + IT_0037;
    const ccomplex_t IT_0082 = conj(IT_0016) + conj(IT_0022);
    const ccomplex_t IT_0083 = IT_0016 + IT_0022;
    const ccomplex_t IT_0084 = (-6)*IT_0059;
    const ccomplex_t IT_0085 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0086 = (-6)*IT_0065;
    const ccomplex_t IT_0087 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0088 = m_t*s_34*m_C_1;
    const ccomplex_t IT_0089 = IT_0038 + IT_0043*(conj(IT_0016)*IT_0048 + conj
      (IT_0051)*IT_0053) + IT_0057*(conj(IT_0034)*IT_0048 + IT_0053*conj(IT_0058
      )) + 6*IT_0059*(conj(IT_0022)*IT_0063 + 0.166666666666667*conj(IT_0051)
      *IT_0064) + 6*(conj(IT_0037)*IT_0063 + 0.166666666666667*conj(IT_0058)
      *IT_0064)*IT_0065 + 6*conj(IT_0059)*(IT_0022*IT_0063 + 0.166666666666667
      *IT_0051*IT_0064 + 0.166666666666667*s_13*IT_0059*IT_0066 +
       0.166666666666667*IT_0043*IT_0068) + 6*conj(IT_0065)*(IT_0037*IT_0063 +
       0.166666666666667*IT_0058*IT_0064 + 0.166666666666667*s_13*IT_0065
      *IT_0066 + 0.166666666666667*IT_0057*IT_0068) + conj(IT_0043)*(IT_0016
      *IT_0048 + IT_0051*IT_0053 + IT_0059*IT_0068 + IT_0043*IT_0069) + conj
      (IT_0057)*(IT_0034*IT_0048 + IT_0053*IT_0058 + IT_0065*IT_0068 + IT_0057
      *IT_0069) + IT_0070*IT_0074 + IT_0044*(conj(IT_0022)*IT_0075 + IT_0022
      *IT_0076 + conj(IT_0037)*IT_0077 + IT_0037*IT_0078) + IT_0079*(IT_0075
      *IT_0080 + IT_0076*IT_0081 + IT_0077*IT_0082 + IT_0078*IT_0083) + IT_0060*
      (conj(IT_0016)*IT_0084 + IT_0016*IT_0085 + conj(IT_0034)*IT_0086 + IT_0034
      *IT_0087) + (IT_0080*IT_0084 + IT_0081*IT_0085 + IT_0082*IT_0086 + IT_0083
      *IT_0087)*IT_0088;
    return create_ccomplex_return(IT_0089);
}

