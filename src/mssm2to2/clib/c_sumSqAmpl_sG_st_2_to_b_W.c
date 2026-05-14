#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_2_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_2_to_b_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_sG, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = -IT_0010 + -IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -2);
    const ccomplex_t IT_0020 = pow(s_14, 2);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0006 + IT_0022;
    const ccomplex_t IT_0024 = s_13*IT_0023;
    const ccomplex_t IT_0025 = pow(m_st_2, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + -IT_0006 + -IT_0025 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0010 + IT_0017 + 2*IT_0032;
    const ccomplex_t IT_0034 = s_14*s_34*IT_0019;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = s_13 + IT_0035;
    const ccomplex_t IT_0037 = s_13*IT_0036;
    const ccomplex_t IT_0038 = m_b*IT_0032;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = m_b*IT_0023;
    const ccomplex_t IT_0041 = s_13*s_14;
    const ccomplex_t IT_0042 = s_34*IT_0006;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0041 + IT_0043;
    const ccomplex_t IT_0045 = (-8)*IT_0018;
    const ccomplex_t IT_0046 = pow(s_34, 2);
    const ccomplex_t IT_0047 = IT_0019*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0005 + IT_0048;
    const ccomplex_t IT_0050 = s_13*IT_0049;
    const ccomplex_t IT_0051 = m_b*IT_0036;
    const ccomplex_t IT_0052 = s_14*IT_0005;
    const ccomplex_t IT_0053 = s_13*s_34;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = (-8)*IT_0033;
    const ccomplex_t IT_0057 = pow(m_W, 2);
    const ccomplex_t IT_0058 = s_13*IT_0057;
    const ccomplex_t IT_0059 = s_14*s_34;
    const ccomplex_t IT_0060 = pow(m_W, 4);
    const ccomplex_t IT_0061 = s_13*IT_0060;
    const ccomplex_t IT_0062 = IT_0019*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0065 = IT_0029*IT_0064;
    const ccomplex_t IT_0066 = IT_0026*IT_0063*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = s_34*IT_0057;
    const ccomplex_t IT_0069 = IT_0019*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = s_34 + IT_0070;
    const ccomplex_t IT_0072 = m_sG*IT_0071;
    const ccomplex_t IT_0073 = 2*IT_0072;
    const ccomplex_t IT_0074 = s_14*IT_0057;
    const ccomplex_t IT_0075 = IT_0019*IT_0074;
    const ccomplex_t IT_0076 = (-0.25)*IT_0075;
    const ccomplex_t IT_0077 = s_14 + IT_0076;
    const ccomplex_t IT_0078 = m_b*IT_0077;
    const ccomplex_t IT_0079 = (-4)*IT_0078;
    const ccomplex_t IT_0080 = (-8)*conj(IT_0018);
    const ccomplex_t IT_0081 = (-8)*conj(IT_0033);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0083 = IT_0014*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0012*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0087 = IT_0007*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = IT_0085 + IT_0089;
    const ccomplex_t IT_0091 = -IT_0085 + -IT_0089;
    const ccomplex_t IT_0092 = m_b*m_sG;
    const ccomplex_t IT_0093 = IT_0023*IT_0092;
    const ccomplex_t IT_0094 = IT_0036*IT_0092;
    const ccomplex_t IT_0095 = m_sG*IT_0036;
    const ccomplex_t IT_0096 = (-8)*IT_0090;
    const ccomplex_t IT_0097 = (-8)*conj(IT_0090);
    const ccomplex_t IT_0098 = IT_0049*IT_0092;
    const ccomplex_t IT_0099 = m_sG*IT_0049;
    const ccomplex_t IT_0100 = (-8)*IT_0091;
    const ccomplex_t IT_0101 = (-8)*conj(IT_0067);
    const ccomplex_t IT_0102 = (-8)*conj(IT_0091);
    const ccomplex_t IT_0103 = (-2)*IT_0059;
    const ccomplex_t IT_0104 = IT_0058 + IT_0103;
    const ccomplex_t IT_0105 = IT_0019*IT_0104;
    const ccomplex_t IT_0106 = (-8)*IT_0105;
    const ccomplex_t IT_0107 = 16*s_13 + IT_0106;
    const ccomplex_t IT_0108 = 8*conj(IT_0032);
    const ccomplex_t IT_0109 = (conj(IT_0018)*IT_0024 + conj(IT_0033)*IT_0037 
      + conj(IT_0039)*IT_0040 + conj(IT_0032)*IT_0044)*IT_0045 + (conj(IT_0033)
      *IT_0050 + conj(IT_0039)*IT_0051 + conj(IT_0032)*IT_0055)*IT_0056 +
       IT_0032*(conj(IT_0032)*((-16)*IT_0058 + 32*IT_0059 + 8*IT_0062) + 8*conj
      (IT_0067)*IT_0073 + 8*conj(IT_0039)*IT_0079 + IT_0044*IT_0080 + IT_0055
      *IT_0081) + (IT_0040*conj(IT_0067) + IT_0024*conj(IT_0090) + IT_0037*conj
      (IT_0091) + conj(IT_0018)*IT_0093 + conj(IT_0033)*IT_0094 + conj(IT_0039)
      *IT_0095)*IT_0096 + IT_0037*IT_0091*IT_0097 + (IT_0051*conj(IT_0067) +
       IT_0050*conj(IT_0091) + conj(IT_0018)*IT_0094 + conj(IT_0033)*IT_0098 +
       conj(IT_0039)*IT_0099)*IT_0100 + IT_0018*(IT_0093*IT_0097 + IT_0095
      *IT_0101 + IT_0094*IT_0102) + IT_0033*(IT_0037*IT_0080 + IT_0094*IT_0097 +
       IT_0099*IT_0101 + IT_0098*IT_0102) + IT_0067*((-24)*conj(IT_0039)*IT_0092
       + (-8)*conj(IT_0018)*IT_0095 + IT_0040*IT_0097 + (-8)*conj(IT_0033)
      *IT_0099 + IT_0051*IT_0102 + conj(IT_0067)*IT_0107 + IT_0073*IT_0108) +
       IT_0039*(IT_0040*IT_0080 + IT_0051*IT_0081 + (-24)*conj(IT_0067)*IT_0092 
      + IT_0095*IT_0097 + IT_0099*IT_0102 + conj(IT_0039)*IT_0107 + IT_0079
      *IT_0108);
    return create_ccomplex_return(IT_0109);
}

