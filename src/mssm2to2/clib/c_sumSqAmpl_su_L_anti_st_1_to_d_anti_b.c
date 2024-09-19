#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_st_1_to_d_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_st_1_to_d_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(beta);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = m_t*V_tb*V_u1*e_em*IT_0005*conj(U_st_10);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0010 + IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0003;
    const ccomplex_t IT_0018 = pow(m_d, 2);
    const ccomplex_t IT_0019 = pow(m_su_L, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = m_t*V_tb*V_u2*e_em*IT_0005*conj(U_st_10);
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0028 + IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0003;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0034*IT_0038;
    const ccomplex_t IT_0040 = -IT_0023 + -IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0042 = cos(beta);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0003*IT_0005*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0016*IT_0045;
    const ccomplex_t IT_0047 = IT_0020*IT_0041*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0003*IT_0005*IT_0043;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0034*IT_0050;
    const ccomplex_t IT_0052 = IT_0036*IT_0048*IT_0051;
    const ccomplex_t IT_0053 = m_d*IT_0028;
    const ccomplex_t IT_0054 = m_d*IT_0032;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = IT_0038*IT_0056;
    const ccomplex_t IT_0058 = m_d*IT_0010;
    const ccomplex_t IT_0059 = m_d*IT_0014;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = IT_0022*IT_0061;
    const ccomplex_t IT_0063 = -IT_0047 + -IT_0052 + IT_0057 + IT_0062;
    const ccomplex_t IT_0064 = m_d*s_14;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0003*IT_0005*IT_0043*conj(U_st_00);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0017*IT_0066;
    const ccomplex_t IT_0068 = IT_0020*IT_0041*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0003*IT_0005*IT_0043*conj(U_st_00);
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0035*IT_0070;
    const ccomplex_t IT_0072 = IT_0036*IT_0048*IT_0071;
    const ccomplex_t IT_0073 = IT_0045*IT_0066;
    const ccomplex_t IT_0074 = IT_0020*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = m_d*IT_0075;
    const ccomplex_t IT_0077 = IT_0050*IT_0070;
    const ccomplex_t IT_0078 = IT_0036*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = m_d*IT_0079;
    const ccomplex_t IT_0081 = IT_0068 + IT_0072 + -IT_0076 + -IT_0080;
    const ccomplex_t IT_0082 = m_b*s_13;
    const ccomplex_t IT_0083 = IT_0075 + IT_0079;
    const ccomplex_t IT_0084 = m_b*m_d*IT_0019;
    const ccomplex_t IT_0085 = s_34*IT_0019;
    const ccomplex_t IT_0086 = (-18)*IT_0085;
    const ccomplex_t IT_0087 = s_13*s_14;
    const ccomplex_t IT_0088 = 36*IT_0087;
    const ccomplex_t IT_0089 = IT_0086 + IT_0088;
    const ccomplex_t IT_0090 = 18*s_34;
    const ccomplex_t IT_0091 = m_b*m_d;
    const ccomplex_t IT_0092 = (-18)*IT_0063;
    const ccomplex_t IT_0093 = (-18)*conj(IT_0063);
    const ccomplex_t IT_0094 = (-18)*conj(IT_0040);
    const ccomplex_t IT_0095 = IT_0040*(18*conj(IT_0063)*IT_0064 + (-18)*conj
      (IT_0081)*IT_0082 + (-18)*conj(IT_0083)*IT_0084 + conj(IT_0040)*IT_0089) +
       IT_0063*(18*conj(IT_0040)*IT_0064 + conj(IT_0063)*IT_0090) + (IT_0082
      *conj(IT_0083) + conj(IT_0081)*IT_0091)*IT_0092 + IT_0081*(18*IT_0064*conj
      (IT_0083) + conj(IT_0081)*IT_0090 + IT_0091*IT_0093 + IT_0082*IT_0094) +
       IT_0083*(18*IT_0064*conj(IT_0081) + conj(IT_0083)*IT_0089 + IT_0082
      *IT_0093 + IT_0084*IT_0094);
    return create_ccomplex_return(IT_0095);
}

