#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sG_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sG_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = s_13*s_24;
    const ccomplex_t IT_0002 = s_12*s_34;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0000 + IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = pow(m_s, 2);
    const ccomplex_t IT_0015 = pow(m_sG, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0011;
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = pow(m_N_4, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0014 + IT_0022 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 8*conj(IT_0019);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0014 + IT_0022 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0027;
    const ccomplex_t IT_0034 = IT_0005*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = 8*IT_0038;
    const ccomplex_t IT_0040 = 8*conj(IT_0038);
    const ccomplex_t IT_0041 = m_sG*s_34*m_N_4;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = conj(N_B4)*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W4)*e_em;
    const ccomplex_t IT_0048 = IT_0010*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0005*IT_0051;
    const ccomplex_t IT_0053 = IT_0035*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = N_B4*e_em;
    const ccomplex_t IT_0056 = IT_0043*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_W4*e_em;
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-3)*IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = IT_0005*IT_0062;
    const ccomplex_t IT_0064 = IT_0023*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = 8*IT_0066;
    const ccomplex_t IT_0068 = 8*conj(IT_0066);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0043;
    const ccomplex_t IT_0070 = (-0.333333333333333)*IT_0069;
    const ccomplex_t IT_0071 = IT_0005*IT_0070;
    const ccomplex_t IT_0072 = IT_0016*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0043;
    const ccomplex_t IT_0075 = (-0.333333333333333)*IT_0074;
    const ccomplex_t IT_0076 = IT_0005*IT_0075;
    const ccomplex_t IT_0077 = IT_0030*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = 8*IT_0079;
    const ccomplex_t IT_0081 = 8*conj(IT_0079);
    const ccomplex_t IT_0082 = s_12*IT_0014;
    const ccomplex_t IT_0083 = m_s*s_23*m_N_4;
    const ccomplex_t IT_0084 = 16*conj(IT_0019);
    const ccomplex_t IT_0085 = conj(IT_0025) + 0.125*IT_0084;
    const ccomplex_t IT_0086 = 8*IT_0032;
    const ccomplex_t IT_0087 = 16*IT_0038;
    const ccomplex_t IT_0088 = IT_0086 + IT_0087;
    const ccomplex_t IT_0089 = 8*conj(IT_0032);
    const ccomplex_t IT_0090 = 16*conj(IT_0038);
    const ccomplex_t IT_0091 = IT_0089 + IT_0090;
    const ccomplex_t IT_0092 = 16*IT_0019;
    const ccomplex_t IT_0093 = m_s*m_sG*s_14;
    const ccomplex_t IT_0094 = m_s*s_24*m_N_4;
    const ccomplex_t IT_0095 = 16*IT_0066;
    const ccomplex_t IT_0096 = 16*conj(IT_0066);
    const ccomplex_t IT_0097 = 16*IT_0079;
    const ccomplex_t IT_0098 = 16*conj(IT_0079);
    const ccomplex_t IT_0099 = m_s*m_sG*s_13;
    const ccomplex_t IT_0100 = 16*IT_0073;
    const ccomplex_t IT_0101 = 16*IT_0032;
    const ccomplex_t IT_0102 = m_sG*m_N_4*IT_0014;
    const ccomplex_t IT_0103 = IT_0026*IT_0032 + IT_0025*(16*conj(IT_0032) +
       IT_0040) + 16*IT_0054*conj(IT_0073) + conj(IT_0019)*IT_0087 + 8*IT_0019*
      (conj(IT_0032) + 0.125*IT_0090) + conj(IT_0079)*IT_0095 + IT_0079*IT_0096 
      + conj(IT_0054)*IT_0100 + conj(IT_0025)*(IT_0039 + IT_0101);
    const ccomplex_t IT_0104 = 8*IT_0004*(IT_0019*conj(IT_0025) + 0.125
      *IT_0025*IT_0026 + 0.125*conj(IT_0032)*IT_0039 + 0.125*IT_0032*IT_0040) +
       IT_0041*(conj(IT_0054)*IT_0067 + IT_0054*IT_0068 + conj(IT_0073)*IT_0080 
      + IT_0073*IT_0081) + (IT_0068*IT_0073 + IT_0067*conj(IT_0073) + conj
      (IT_0054)*IT_0080 + IT_0054*IT_0081)*IT_0082 + 8*IT_0083*(IT_0073*IT_0085 
      + 0.125*conj(IT_0054)*IT_0088 + 0.125*IT_0054*IT_0091 + conj(IT_0073)*
      (IT_0025 + 0.125*IT_0092)) + 8*(IT_0054*IT_0085 + 0.125*conj(IT_0073)
      *IT_0088 + 0.125*IT_0073*IT_0091 + conj(IT_0054)*(IT_0025 + 0.125*IT_0092)
      )*IT_0093 + IT_0094*(conj(IT_0019)*IT_0067 + IT_0019*IT_0068 + IT_0040
      *IT_0079 + IT_0039*conj(IT_0079) + conj(IT_0025)*IT_0095 + IT_0025*IT_0096
       + conj(IT_0032)*IT_0097 + IT_0032*IT_0098) + (IT_0040*IT_0066 + IT_0039
      *conj(IT_0066) + conj(IT_0019)*IT_0080 + IT_0019*IT_0081 + conj(IT_0032)
      *IT_0095 + IT_0032*IT_0096 + conj(IT_0025)*IT_0097 + IT_0025*IT_0098)
      *IT_0099 + IT_0000*(16*IT_0054*conj(IT_0054) + conj(IT_0038)*IT_0087 +
       conj(IT_0019)*IT_0092 + conj(IT_0073)*IT_0100) + IT_0001*(16*IT_0025*conj
      (IT_0025) + conj(IT_0066)*IT_0095 + conj(IT_0079)*IT_0097 + conj(IT_0032)
      *IT_0101) + IT_0102*IT_0103;
    return create_ccomplex_return(IT_0104);
}

