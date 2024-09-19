#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sG_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sG_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = m_sG*s_34*m_N_1;
    const ccomplex_t IT_0001 = pow(m_c, 2);
    const ccomplex_t IT_0002 = pow(m_sG, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_W1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0009 + 3*IT_0014);
    const ccomplex_t IT_0016 = 0.166666666666667*IT_0015;
    const ccomplex_t IT_0017 = IT_0004*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_N_1, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0001 + IT_0020 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0022 = N_B1*e_em;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W1*e_em;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + 3*IT_0027);
    const ccomplex_t IT_0029 = 0.166666666666667*IT_0028;
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = IT_0021*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = 8*IT_0033;
    const ccomplex_t IT_0035 = 8*conj(IT_0033);
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0006;
    const ccomplex_t IT_0038 = 0.666666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0001 + IT_0020 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0006;
    const ccomplex_t IT_0044 = 0.666666666666667*IT_0043;
    const ccomplex_t IT_0045 = IT_0004*IT_0044;
    const ccomplex_t IT_0046 = IT_0042*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = 8*IT_0048;
    const ccomplex_t IT_0050 = 8*conj(IT_0048);
    const ccomplex_t IT_0051 = s_12*IT_0001;
    const ccomplex_t IT_0052 = s_14*s_23;
    const ccomplex_t IT_0053 = s_13*s_24;
    const ccomplex_t IT_0054 = s_12*s_34;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0052 + IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = pow(m_W, -1);
    const ccomplex_t IT_0058 = sin(beta);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0011*IT_0057*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0004*IT_0061;
    const ccomplex_t IT_0063 = IT_0042*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0060;
    const ccomplex_t IT_0066 = IT_0004*IT_0065;
    const ccomplex_t IT_0067 = IT_0003*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = 8*IT_0069;
    const ccomplex_t IT_0071 = 8*conj(IT_0069);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0011*IT_0057*IT_0059;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0004*IT_0073;
    const ccomplex_t IT_0075 = IT_0021*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0072;
    const ccomplex_t IT_0078 = IT_0004*IT_0077;
    const ccomplex_t IT_0079 = IT_0036*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = 8*IT_0081;
    const ccomplex_t IT_0083 = 8*conj(IT_0081);
    const ccomplex_t IT_0084 = m_c*s_24*m_N_1;
    const ccomplex_t IT_0085 = 16*IT_0033;
    const ccomplex_t IT_0086 = 16*conj(IT_0033);
    const ccomplex_t IT_0087 = 16*IT_0048;
    const ccomplex_t IT_0088 = 16*conj(IT_0048);
    const ccomplex_t IT_0089 = m_c*m_sG*s_13;
    const ccomplex_t IT_0090 = 16*IT_0069;
    const ccomplex_t IT_0091 = 16*IT_0081;
    const ccomplex_t IT_0092 = 16*IT_0041;
    const ccomplex_t IT_0093 = 16*IT_0064;
    const ccomplex_t IT_0094 = 16*IT_0076;
    const ccomplex_t IT_0095 = m_sG*m_N_1*IT_0001;
    const ccomplex_t IT_0096 = 16*conj(IT_0069);
    const ccomplex_t IT_0097 = 16*IT_0019*conj(IT_0041) + IT_0071*IT_0076 +
       IT_0064*IT_0083 + conj(IT_0048)*IT_0085 + IT_0048*IT_0086 + conj(IT_0081)
      *IT_0090 + conj(IT_0019)*IT_0092 + conj(IT_0076)*(IT_0070 + IT_0093) +
       conj(IT_0064)*(IT_0082 + IT_0094) + IT_0081*IT_0096;
    const ccomplex_t IT_0098 = m_c*s_23*m_N_1;
    const ccomplex_t IT_0099 = 16*conj(IT_0081);
    const ccomplex_t IT_0100 = 0.125*IT_0099;
    const ccomplex_t IT_0101 = conj(IT_0076) + IT_0100;
    const ccomplex_t IT_0102 = 8*IT_0064;
    const ccomplex_t IT_0103 = IT_0090 + IT_0102;
    const ccomplex_t IT_0104 = 8*conj(IT_0064);
    const ccomplex_t IT_0105 = IT_0096 + IT_0104;
    const ccomplex_t IT_0106 = 0.125*IT_0091;
    const ccomplex_t IT_0107 = m_c*m_sG*s_14;
    const ccomplex_t IT_0108 = IT_0000*(conj(IT_0019)*IT_0034 + IT_0019
      *IT_0035 + conj(IT_0041)*IT_0049 + IT_0041*IT_0050) + (IT_0035*IT_0041 +
       IT_0034*conj(IT_0041) + conj(IT_0019)*IT_0049 + IT_0019*IT_0050)*IT_0051 
      + IT_0056*(conj(IT_0064)*IT_0070 + IT_0064*IT_0071 + conj(IT_0076)*IT_0082
       + IT_0076*IT_0083) + IT_0084*(conj(IT_0048)*IT_0070 + IT_0048*IT_0071 +
       IT_0035*IT_0081 + IT_0034*conj(IT_0081) + conj(IT_0076)*IT_0085 + IT_0076
      *IT_0086 + conj(IT_0064)*IT_0087 + IT_0064*IT_0088) + (conj(IT_0033)
      *IT_0070 + IT_0033*IT_0071 + IT_0050*IT_0081 + IT_0049*conj(IT_0081) +
       conj(IT_0064)*IT_0085 + IT_0064*IT_0086 + conj(IT_0076)*IT_0087 + IT_0076
      *IT_0088)*IT_0089 + IT_0052*(16*IT_0019*conj(IT_0019) + conj(IT_0069)
      *IT_0090 + conj(IT_0081)*IT_0091 + conj(IT_0041)*IT_0092) + IT_0053*(conj
      (IT_0033)*IT_0085 + conj(IT_0048)*IT_0087 + conj(IT_0064)*IT_0093 + conj
      (IT_0076)*IT_0094) + IT_0095*IT_0097 + 8*IT_0098*(IT_0041*IT_0101 + 0.125
      *conj(IT_0019)*IT_0103 + 0.125*IT_0019*IT_0105 + conj(IT_0041)*(IT_0076 +
       IT_0106)) + 8*(IT_0019*IT_0101 + 0.125*conj(IT_0041)*IT_0103 + 0.125
      *IT_0041*IT_0105 + conj(IT_0019)*(IT_0076 + IT_0106))*IT_0107;
    return create_ccomplex_return(IT_0108);
}

