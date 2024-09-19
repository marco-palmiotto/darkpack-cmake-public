#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_4_to_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_4_to_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = pow(m_C_1, 2);
    const ccomplex_t IT_0022 = pow(m_N_4, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0000*IT_0025;
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (-2)*IT_0030;
    const ccomplex_t IT_0032 = IT_0009*IT_0027;
    const ccomplex_t IT_0033 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = IT_0035 + IT_0038;
    const ccomplex_t IT_0040 = IT_0031 + IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0000*IT_0042;
    const ccomplex_t IT_0044 = pow(m_Hp, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0022 + IT_0044 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 2*IT_0000;
    const ccomplex_t IT_0049 = IT_0023*IT_0042;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (-2)*IT_0047 + IT_0051;
    const ccomplex_t IT_0053 = m_C_1*m_N_4;
    const ccomplex_t IT_0054 = 2*s_14*IT_0053;
    const ccomplex_t IT_0055 = 2*s_24*IT_0053;
    const ccomplex_t IT_0056 = IT_0000*IT_0020;
    const ccomplex_t IT_0057 = IT_0045*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0025*IT_0048;
    const ccomplex_t IT_0060 = 2*IT_0058 + -IT_0059;
    const ccomplex_t IT_0061 = (-2)*s_12;
    const ccomplex_t IT_0062 = s_14*IT_0061;
    const ccomplex_t IT_0063 = -IT_0059;
    const ccomplex_t IT_0064 = s_24*IT_0061;
    const ccomplex_t IT_0065 = IT_0000*IT_0050;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = s_12*s_24;
    const ccomplex_t IT_0068 = s_14*IT_0022;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = 2*s_12*IT_0053;
    const ccomplex_t IT_0072 = 2*IT_0022*IT_0053;
    const ccomplex_t IT_0073 = -IT_0058;
    const ccomplex_t IT_0074 = s_24*m_C_1*m_N_4;
    const ccomplex_t IT_0075 = 2*IT_0074;
    const ccomplex_t IT_0076 = IT_0022*IT_0061;
    const ccomplex_t IT_0077 = pow(s_12, 2);
    const ccomplex_t IT_0078 = (-2)*IT_0077;
    const ccomplex_t IT_0079 = (-0.25)*conj(IT_0026)*IT_0055 + (-0.25)*IT_0064
      *conj(IT_0066) + conj(IT_0047)*IT_0070 + (-0.25)*conj(IT_0060)*IT_0071 + (
      -0.25)*conj(IT_0063)*IT_0072 + (-0.25)*conj(IT_0073)*IT_0075 + (-0.25)
      *conj(IT_0051)*IT_0076 + (-0.25)*conj(IT_0052)*IT_0078;
    const ccomplex_t IT_0080 = s_14*s_24;
    const ccomplex_t IT_0081 = s_14*m_C_1*m_N_4;
    const ccomplex_t IT_0082 = 2*IT_0081;
    const ccomplex_t IT_0083 = (-4)*IT_0070;
    const ccomplex_t IT_0084 = s_24*IT_0021;
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*conj(IT_0063);
    const ccomplex_t IT_0087 = 2*IT_0021*IT_0053;
    const ccomplex_t IT_0088 = 0.5*conj(IT_0052);
    const ccomplex_t IT_0089 = 0.5*conj(IT_0051);
    const ccomplex_t IT_0090 = 0.5*conj(IT_0066);
    const ccomplex_t IT_0091 = 0.5*conj(IT_0073);
    const ccomplex_t IT_0092 = 0.5*conj(IT_0026);
    const ccomplex_t IT_0093 = (-0.5)*conj(IT_0051);
    const ccomplex_t IT_0094 = (-0.5)*conj(IT_0026)*IT_0054 + (-0.5)*IT_0021
      *conj(IT_0052)*IT_0061 + (-0.5)*IT_0062*conj(IT_0066) + (-0.5)*conj
      (IT_0063)*IT_0071 + (-0.5)*conj(IT_0073)*IT_0082 + conj(IT_0047)*IT_0084 +
       (-0.5)*conj(IT_0060)*IT_0087 + IT_0078*IT_0093;
    const ccomplex_t IT_0095 = IT_0026*(conj(IT_0052)*IT_0054 + conj(IT_0051)
      *IT_0055 + conj(IT_0060)*IT_0062 + conj(IT_0063)*IT_0064) + (IT_0054*conj
      (IT_0060) + conj(IT_0052)*IT_0062 + IT_0055*conj(IT_0063) + conj(IT_0051)
      *IT_0064)*IT_0066 + (-4)*IT_0051*IT_0079 + 8*IT_0073*(0.125*conj(IT_0051)
      *IT_0075 + conj(IT_0073)*IT_0080 + 0.125*conj(IT_0052)*IT_0082 + 0.125
      *conj(IT_0063)*IT_0083 + 0.125*conj(IT_0060)*IT_0085) + 2*IT_0060*(0.5
      *IT_0021*conj(IT_0060)*IT_0061 + conj(IT_0047)*IT_0081 + IT_0078*IT_0086 +
       IT_0087*IT_0088 + IT_0071*IT_0089 + IT_0054*IT_0090 + IT_0085*IT_0091 +
       IT_0062*IT_0092) + 2*IT_0063*(conj(IT_0047)*IT_0074 + 0.5*conj(IT_0060)
      *IT_0078 + IT_0076*IT_0086 + IT_0071*IT_0088 + IT_0072*IT_0089 + IT_0055
      *IT_0090 + IT_0083*IT_0091 + IT_0064*IT_0092) + 8*IT_0047*(0.25*conj
      (IT_0063)*IT_0074 + conj(IT_0047)*IT_0080 + 0.25*conj(IT_0060)*IT_0081 + (
      -0.25)*conj(IT_0052)*IT_0084 + IT_0070*IT_0093) + (-2)*IT_0052*IT_0094;
    return create_ccomplex_return(IT_0095);
}

