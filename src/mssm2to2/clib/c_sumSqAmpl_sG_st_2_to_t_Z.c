#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_2_to_t_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_2_to_t_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t m_t = param->m_t;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = e_em*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0003, -1);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = e_em*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + (-3)*IT_0008);
    const ccomplex_t IT_0010 = (-0.166666666666667)*IT_0009;
    const ccomplex_t IT_0011 = IT_0000*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sG, 2);
    const ccomplex_t IT_0013 = pow(m_st_2, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_Z, -2);
    const ccomplex_t IT_0018 = pow(m_Z, 2);
    const ccomplex_t IT_0019 = s_14*IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (-0.25)*IT_0020;
    const ccomplex_t IT_0022 = s_14 + IT_0021;
    const ccomplex_t IT_0023 = m_t*IT_0022;
    const ccomplex_t IT_0024 = (-4)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0003;
    const ccomplex_t IT_0027 = (-0.666666666666667)*IT_0026;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = IT_0014*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = s_34*IT_0018;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = s_34 + IT_0034;
    const ccomplex_t IT_0036 = m_sG*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0039 = IT_0000*IT_0027;
    const ccomplex_t IT_0040 = IT_0014*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = m_t*IT_0016;
    const ccomplex_t IT_0042 = IT_0040 + -IT_0041;
    const ccomplex_t IT_0043 = 8*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0002*IT_0003 
      + (-3)*IT_0001*IT_0006)*conj(U_st_00)*U_st_01 + IT_0002*IT_0003*conj
      (U_st_10)*U_st_11);
    const ccomplex_t IT_0045 = 1.33333333333333*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0047 = pow(m_t, 2);
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0012 + IT_0047 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0002*IT_0003 
      + (-3)*IT_0001*IT_0006)*U_st_01*conj(U_st_01) + IT_0002*IT_0003*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0053 = 1.33333333333333*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0012 + IT_0047 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0025*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = IT_0051 + IT_0057;
    const ccomplex_t IT_0059 = pow(s_14, 2);
    const ccomplex_t IT_0060 = IT_0017*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0012 + IT_0061;
    const ccomplex_t IT_0063 = s_13*IT_0062;
    const ccomplex_t IT_0064 = (-2)*IT_0030 + -IT_0051 + -IT_0057;
    const ccomplex_t IT_0065 = s_14*s_34*IT_0017;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = s_13 + IT_0066;
    const ccomplex_t IT_0068 = s_13*IT_0067;
    const ccomplex_t IT_0069 = IT_0010*IT_0025;
    const ccomplex_t IT_0070 = IT_0014*IT_0038*IT_0069;
    const ccomplex_t IT_0071 = m_t*IT_0030;
    const ccomplex_t IT_0072 = -IT_0070 + IT_0071;
    const ccomplex_t IT_0073 = m_t*IT_0062;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0075 = IT_0048*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0045*IT_0076;
    const ccomplex_t IT_0078 = IT_0000*IT_0054;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0053*IT_0079;
    const ccomplex_t IT_0081 = -IT_0077 + -IT_0080;
    const ccomplex_t IT_0082 = m_t*m_sG;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = 2*IT_0016 + IT_0077 + IT_0080;
    const ccomplex_t IT_0085 = IT_0067*IT_0082;
    const ccomplex_t IT_0086 = m_sG*IT_0067;
    const ccomplex_t IT_0087 = s_13*s_14;
    const ccomplex_t IT_0088 = s_34*IT_0012;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = (-8)*IT_0058;
    const ccomplex_t IT_0092 = (-8)*conj(IT_0058);
    const ccomplex_t IT_0093 = pow(s_34, 2);
    const ccomplex_t IT_0094 = IT_0017*IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0047 + IT_0095;
    const ccomplex_t IT_0097 = s_13*IT_0096;
    const ccomplex_t IT_0098 = m_t*IT_0067;
    const ccomplex_t IT_0099 = IT_0082*IT_0096;
    const ccomplex_t IT_0100 = m_sG*IT_0096;
    const ccomplex_t IT_0101 = s_14*IT_0047;
    const ccomplex_t IT_0102 = s_13*s_34;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = (-8)*IT_0064;
    const ccomplex_t IT_0106 = IT_0018*IT_0082;
    const ccomplex_t IT_0107 = (-3)*IT_0106;
    const ccomplex_t IT_0108 = pow(m_Z, 4);
    const ccomplex_t IT_0109 = s_13*IT_0108;
    const ccomplex_t IT_0110 = IT_0017*IT_0109;
    const ccomplex_t IT_0111 = 8*IT_0110;
    const ccomplex_t IT_0112 = s_14*s_34;
    const ccomplex_t IT_0113 = 32*IT_0112;
    const ccomplex_t IT_0114 = s_13*IT_0018;
    const ccomplex_t IT_0115 = (-16)*IT_0114;
    const ccomplex_t IT_0116 = IT_0111 + IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = 8*conj(IT_0042);
    const ccomplex_t IT_0118 = (-8)*conj(IT_0064);
    const ccomplex_t IT_0119 = (-8)*IT_0072;
    const ccomplex_t IT_0120 = (-8)*conj(IT_0072);
    const ccomplex_t IT_0121 = (-8)*IT_0081;
    const ccomplex_t IT_0122 = (-8)*conj(IT_0081);
    const ccomplex_t IT_0123 = (-8)*IT_0084;
    const ccomplex_t IT_0124 = (-8)*conj(IT_0084);
    const ccomplex_t IT_0125 = (-2)*IT_0112;
    const ccomplex_t IT_0126 = IT_0114 + IT_0125;
    const ccomplex_t IT_0127 = IT_0017*IT_0126;
    const ccomplex_t IT_0128 = (-8)*IT_0127;
    const ccomplex_t IT_0129 = 16*s_13;
    const ccomplex_t IT_0130 = (conj(IT_0016)*IT_0024 + conj(IT_0031)*IT_0037)
      *IT_0043 + (conj(IT_0058)*IT_0063 + conj(IT_0064)*IT_0068 + conj(IT_0072)
      *IT_0073 + conj(IT_0081)*IT_0083 + conj(IT_0084)*IT_0085 + conj(IT_0042)
      *IT_0086 + conj(IT_0031)*IT_0090)*IT_0091 + IT_0064*IT_0068*IT_0092 + 
      (conj(IT_0081)*IT_0085 + conj(IT_0064)*IT_0097 + conj(IT_0072)*IT_0098 +
       conj(IT_0084)*IT_0099 + conj(IT_0042)*IT_0100 + conj(IT_0031)*IT_0104)
      *IT_0105 + 8*IT_0031*(IT_0024*conj(IT_0072) + 0.125*IT_0090*IT_0092 + conj
      (IT_0016)*IT_0107 + 0.125*conj(IT_0031)*IT_0116 + 0.125*IT_0037*IT_0117 +
       0.125*IT_0104*IT_0118) + (conj(IT_0081)*IT_0086 + conj(IT_0084)*IT_0100)
      *IT_0119 + IT_0081*(IT_0083*IT_0092 + IT_0085*IT_0118 + IT_0086*IT_0120) +
       (conj(IT_0042)*IT_0073 + IT_0063*conj(IT_0081) + IT_0068*conj(IT_0084) +
       conj(IT_0016)*IT_0090)*IT_0121 + IT_0084*(IT_0085*IT_0092 + IT_0099
      *IT_0118 + IT_0100*IT_0120 + IT_0068*IT_0122) + (conj(IT_0084)*IT_0097 +
       conj(IT_0042)*IT_0098 + conj(IT_0016)*IT_0104)*IT_0123 + 8*IT_0016*
      (IT_0037*conj(IT_0072) + conj(IT_0031)*IT_0107 + 0.125*conj(IT_0016)
      *IT_0116 + 0.125*IT_0024*IT_0117 + 0.125*IT_0090*IT_0122 + 0.125*IT_0104
      *IT_0124) + IT_0042*((-24)*conj(IT_0072)*IT_0082 + IT_0086*IT_0092 +
       IT_0100*IT_0118 + IT_0073*IT_0122 + IT_0098*IT_0124 + conj(IT_0042)*
      (IT_0128 + IT_0129)) + 8*IT_0072*(IT_0024*conj(IT_0031) + conj(IT_0016)
      *IT_0037 + (-3)*conj(IT_0042)*IT_0082 + 1./8*IT_0073*IT_0092 + 1./8
      *IT_0098*IT_0118 + 1./8*conj(IT_0072)*(IT_0128 + IT_0129));
    return create_ccomplex_return(IT_0130);
}

