#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_2_to_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_2_to_b_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0007, -1);
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0009 + 3*IT_0012);
    const ccomplex_t IT_0014 = (-0.166666666666667)*IT_0013;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_Z, -2);
    const ccomplex_t IT_0019 = pow(m_Z, 2);
    const ccomplex_t IT_0020 = s_14*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (-0.25)*IT_0021;
    const ccomplex_t IT_0023 = s_14 + IT_0022;
    const ccomplex_t IT_0024 = m_b*IT_0023;
    const ccomplex_t IT_0025 = (-4)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0007;
    const ccomplex_t IT_0028 = 0.333333333333333*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = s_34*IT_0019;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = s_34 + IT_0035;
    const ccomplex_t IT_0037 = m_sG*IT_0036;
    const ccomplex_t IT_0038 = 2*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0040 = IT_0004*IT_0028;
    const ccomplex_t IT_0041 = IT_0003*IT_0039*IT_0040;
    const ccomplex_t IT_0042 = m_b*IT_0017;
    const ccomplex_t IT_0043 = IT_0041 + -IT_0042;
    const ccomplex_t IT_0044 = 8*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0006
      *IT_0007 + 3*IT_0005*IT_0010)*conj(U_sb_00)*U_sb_01 + IT_0006*IT_0007*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0046 = (-0.666666666666667)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0046*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0006
      *IT_0007 + 3*IT_0005*IT_0010)*U_sb_01*conj(U_sb_01) + IT_0006*IT_0007
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0053 = (-0.666666666666667)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0026*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = IT_0051 + IT_0057;
    const ccomplex_t IT_0059 = pow(s_14, 2);
    const ccomplex_t IT_0060 = IT_0018*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0001 + IT_0061;
    const ccomplex_t IT_0063 = s_13*IT_0062;
    const ccomplex_t IT_0064 = (-2)*IT_0031 + -IT_0051 + -IT_0057;
    const ccomplex_t IT_0065 = s_14*s_34*IT_0018;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = s_13 + IT_0066;
    const ccomplex_t IT_0068 = s_13*IT_0067;
    const ccomplex_t IT_0069 = IT_0014*IT_0026;
    const ccomplex_t IT_0070 = IT_0003*IT_0039*IT_0069;
    const ccomplex_t IT_0071 = m_b*IT_0031;
    const ccomplex_t IT_0072 = -IT_0070 + IT_0071;
    const ccomplex_t IT_0073 = m_b*IT_0062;
    const ccomplex_t IT_0074 = IT_0004*IT_0054;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0053*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0078 = IT_0048*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0046*IT_0079;
    const ccomplex_t IT_0081 = -IT_0076 + -IT_0080;
    const ccomplex_t IT_0082 = m_b*m_sG;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = 2*IT_0017 + IT_0076 + IT_0080;
    const ccomplex_t IT_0085 = IT_0067*IT_0082;
    const ccomplex_t IT_0086 = m_sG*IT_0067;
    const ccomplex_t IT_0087 = s_13*s_14;
    const ccomplex_t IT_0088 = s_34*IT_0001;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = (-8)*IT_0058;
    const ccomplex_t IT_0092 = (-8)*conj(IT_0058);
    const ccomplex_t IT_0093 = pow(s_34, 2);
    const ccomplex_t IT_0094 = IT_0018*IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0000 + IT_0095;
    const ccomplex_t IT_0097 = s_13*IT_0096;
    const ccomplex_t IT_0098 = m_b*IT_0067;
    const ccomplex_t IT_0099 = IT_0082*IT_0096;
    const ccomplex_t IT_0100 = m_sG*IT_0096;
    const ccomplex_t IT_0101 = s_14*IT_0000;
    const ccomplex_t IT_0102 = s_13*s_34;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = (-8)*IT_0064;
    const ccomplex_t IT_0106 = (-8)*IT_0081;
    const ccomplex_t IT_0107 = (-8)*IT_0084;
    const ccomplex_t IT_0108 = (-8)*IT_0072;
    const ccomplex_t IT_0109 = (-8)*conj(IT_0064);
    const ccomplex_t IT_0110 = (-8)*conj(IT_0072);
    const ccomplex_t IT_0111 = (-8)*conj(IT_0081);
    const ccomplex_t IT_0112 = s_13*IT_0019;
    const ccomplex_t IT_0113 = s_14*s_34;
    const ccomplex_t IT_0114 = (-2)*IT_0113;
    const ccomplex_t IT_0115 = IT_0112 + IT_0114;
    const ccomplex_t IT_0116 = IT_0018*IT_0115;
    const ccomplex_t IT_0117 = (-8)*IT_0116;
    const ccomplex_t IT_0118 = 16*s_13;
    const ccomplex_t IT_0119 = (-8)*conj(IT_0084);
    const ccomplex_t IT_0120 = IT_0019*IT_0082;
    const ccomplex_t IT_0121 = (-3)*IT_0120;
    const ccomplex_t IT_0122 = pow(m_Z, 4);
    const ccomplex_t IT_0123 = s_13*IT_0122;
    const ccomplex_t IT_0124 = IT_0018*IT_0123;
    const ccomplex_t IT_0125 = 8*IT_0124;
    const ccomplex_t IT_0126 = 32*IT_0113;
    const ccomplex_t IT_0127 = (-16)*IT_0112;
    const ccomplex_t IT_0128 = IT_0125 + IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = 8*conj(IT_0043);
    const ccomplex_t IT_0130 = 0.125*IT_0090;
    const ccomplex_t IT_0131 = 0.125*IT_0104;
    const ccomplex_t IT_0132 = (conj(IT_0017)*IT_0025 + conj(IT_0032)*IT_0038)
      *IT_0044 + (conj(IT_0058)*IT_0063 + conj(IT_0064)*IT_0068 + conj(IT_0072)
      *IT_0073 + conj(IT_0081)*IT_0083 + conj(IT_0084)*IT_0085 + conj(IT_0043)
      *IT_0086 + conj(IT_0032)*IT_0090)*IT_0091 + IT_0064*IT_0068*IT_0092 + 
      (conj(IT_0081)*IT_0085 + conj(IT_0064)*IT_0097 + conj(IT_0072)*IT_0098 +
       conj(IT_0084)*IT_0099 + conj(IT_0043)*IT_0100 + conj(IT_0032)*IT_0104)
      *IT_0105 + (conj(IT_0043)*IT_0073 + IT_0063*conj(IT_0081) + IT_0068*conj
      (IT_0084) + conj(IT_0017)*IT_0090)*IT_0106 + (conj(IT_0084)*IT_0097 + conj
      (IT_0043)*IT_0098 + conj(IT_0017)*IT_0104)*IT_0107 + (conj(IT_0081)
      *IT_0086 + conj(IT_0084)*IT_0100)*IT_0108 + IT_0081*(IT_0083*IT_0092 +
       IT_0085*IT_0109 + IT_0086*IT_0110) + IT_0084*(IT_0085*IT_0092 + IT_0099
      *IT_0109 + IT_0100*IT_0110 + IT_0068*IT_0111) + 8*IT_0072*(IT_0025*conj
      (IT_0032) + conj(IT_0017)*IT_0038 + (-3)*conj(IT_0043)*IT_0082 + 1./8
      *IT_0073*IT_0092 + 1./8*IT_0098*IT_0109 + 1./8*conj(IT_0072)*(IT_0117 +
       IT_0118)) + IT_0043*((-24)*conj(IT_0072)*IT_0082 + IT_0086*IT_0092 +
       IT_0100*IT_0109 + IT_0073*IT_0111 + conj(IT_0043)*(IT_0117 + IT_0118) +
       IT_0098*IT_0119) + 8*IT_0032*(IT_0025*conj(IT_0072) + conj(IT_0017)
      *IT_0121 + 0.125*conj(IT_0032)*IT_0128 + 0.125*IT_0038*IT_0129 + IT_0092
      *IT_0130 + IT_0109*IT_0131) + 8*IT_0017*(IT_0038*conj(IT_0072) + conj
      (IT_0032)*IT_0121 + 0.125*conj(IT_0017)*IT_0128 + 0.125*IT_0025*IT_0129 +
       IT_0111*IT_0130 + IT_0119*IT_0131);
    return create_ccomplex_return(IT_0132);
}

