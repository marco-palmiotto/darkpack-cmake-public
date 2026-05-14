#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_A_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0003
      *IT_0004 + 3*IT_0002*IT_0005)*U_sb_00*conj(U_sb_01) + IT_0003*IT_0004
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0007 = (-0.333333333333333)*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_sb_1
      *Gamma_b1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = (-0.666666666666667)*IT_0006;
    const ccomplex_t IT_0014 = IT_0001*IT_0013;
    const ccomplex_t IT_0015 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_sb_2
      *Gamma_b2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (IT_0008 + -IT_0014)*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = pow(m_Z, -2);
    const ccomplex_t IT_0020 = s_14*s_24*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = s_12 + IT_0021;
    const ccomplex_t IT_0023 = 3*s_12;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0017;
    const ccomplex_t IT_0026 = s_14*s_34*IT_0019;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = s_13 + IT_0027;
    const ccomplex_t IT_0029 = IT_0023*IT_0028;
    const ccomplex_t IT_0030 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0031 = IT_0007*IT_0030;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = -IT_0031;
    const ccomplex_t IT_0034 = IT_0008 + IT_0033;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0032 + 0.5*IT_0035;
    const ccomplex_t IT_0037 = pow(m_sb_1, 2);
    const ccomplex_t IT_0038 = pow(s_14, 2);
    const ccomplex_t IT_0039 = IT_0019*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0037 + IT_0040;
    const ccomplex_t IT_0042 = 3*s_13;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0013*IT_0030;
    const ccomplex_t IT_0045 = IT_0008 + IT_0044;
    const ccomplex_t IT_0046 = -IT_0014 + -IT_0031;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = IT_0010*IT_0047;
    const ccomplex_t IT_0049 = -IT_0044;
    const ccomplex_t IT_0050 = IT_0031 + IT_0049;
    const ccomplex_t IT_0051 = IT_0016*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0048 + (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0022*IT_0042;
    const ccomplex_t IT_0054 = (-0.5)*IT_0035 + 0.5*IT_0051;
    const ccomplex_t IT_0055 = IT_0028*IT_0042;
    const ccomplex_t IT_0056 = pow(e_em, 2);
    const ccomplex_t IT_0057 = cpow(IT_0002, -2);
    const ccomplex_t IT_0058 = cpow(IT_0004, -2);
    const ccomplex_t IT_0059 = tan(theta_W);
    const ccomplex_t IT_0060 = cpow(IT_0059, 2);
    const ccomplex_t IT_0061 = cpow(1 + IT_0060, -1);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0003*IT_0004*IT_0056
      *IT_0061*((-0.5)*(IT_0057 + 3*IT_0058)*U_sb_00*conj(U_sb_01) + IT_0057
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0063 = (-0.222222222222222)*IT_0062;
    const ccomplex_t IT_0064 = IT_0008*IT_0016;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = 3*IT_0037;
    const ccomplex_t IT_0067 = IT_0010*(IT_0008 + -IT_0014);
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0022*IT_0066;
    const ccomplex_t IT_0070 = 0.5*IT_0011;
    const ccomplex_t IT_0071 = IT_0028*IT_0066;
    const ccomplex_t IT_0072 = conj(IT_0018)*IT_0024 + conj(IT_0025)*IT_0029 +
       conj(IT_0036)*IT_0043 + conj(IT_0052)*IT_0053 + conj(IT_0054)*IT_0055 + 3
      *IT_0041*(conj(IT_0063) + 0.333333333333333*IT_0023*conj(IT_0065) +
       0.333333333333333*conj(IT_0012)*IT_0066) + conj(IT_0068)*IT_0069 + conj
      (IT_0070)*IT_0071;
    const ccomplex_t IT_0073 = pow(m_sb_2, 2);
    const ccomplex_t IT_0074 = pow(s_24, 2);
    const ccomplex_t IT_0075 = IT_0019*IT_0074;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0073 + IT_0076;
    const ccomplex_t IT_0078 = IT_0023*IT_0077;
    const ccomplex_t IT_0079 = s_24*s_34*IT_0019;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = s_23 + IT_0080;
    const ccomplex_t IT_0082 = IT_0023*IT_0081;
    const ccomplex_t IT_0083 = IT_0042*IT_0077;
    const ccomplex_t IT_0084 = IT_0042*IT_0081;
    const ccomplex_t IT_0085 = IT_0066*IT_0081;
    const ccomplex_t IT_0086 = 3*IT_0022;
    const ccomplex_t IT_0087 = conj(IT_0063)*IT_0086;
    const ccomplex_t IT_0088 = IT_0068*(conj(IT_0036)*IT_0053 + IT_0024*conj
      (IT_0065) + conj(IT_0012)*IT_0069 + IT_0066*conj(IT_0068)*IT_0077 + conj
      (IT_0018)*IT_0078 + conj(IT_0025)*IT_0082 + conj(IT_0052)*IT_0083 + conj
      (IT_0054)*IT_0084 + conj(IT_0070)*IT_0085 + IT_0087);
    const ccomplex_t IT_0089 = 3*s_23;
    const ccomplex_t IT_0090 = IT_0041*IT_0089;
    const ccomplex_t IT_0091 = IT_0022*IT_0089;
    const ccomplex_t IT_0092 = IT_0028*IT_0089;
    const ccomplex_t IT_0093 = 3*IT_0073;
    const ccomplex_t IT_0094 = IT_0022*IT_0093;
    const ccomplex_t IT_0095 = IT_0028*IT_0093;
    const ccomplex_t IT_0096 = IT_0065*(IT_0024*conj(IT_0068) + IT_0029*conj
      (IT_0070) + IT_0087 + conj(IT_0036)*IT_0090 + conj(IT_0052)*IT_0091 + conj
      (IT_0054)*IT_0092 + IT_0041*(conj(IT_0012)*IT_0023 + conj(IT_0065)*IT_0093
      ) + conj(IT_0018)*IT_0094 + conj(IT_0025)*IT_0095);
    const ccomplex_t IT_0097 = 3*IT_0028;
    const ccomplex_t IT_0098 = conj(IT_0063)*IT_0097;
    const ccomplex_t IT_0099 = IT_0077*IT_0089;
    const ccomplex_t IT_0100 = IT_0081*IT_0089;
    const ccomplex_t IT_0101 = IT_0081*IT_0093;
    const ccomplex_t IT_0102 = 3*IT_0077;
    const ccomplex_t IT_0103 = conj(IT_0012)*IT_0024 + conj(IT_0068)*IT_0078 +
       conj(IT_0070)*IT_0082 + conj(IT_0036)*IT_0091 + conj(IT_0018)*IT_0077
      *IT_0093 + conj(IT_0065)*IT_0094 + conj(IT_0052)*IT_0099 + conj(IT_0054)
      *IT_0100 + conj(IT_0025)*IT_0101 + conj(IT_0063)*IT_0102;
    const ccomplex_t IT_0104 = 3*IT_0081;
    const ccomplex_t IT_0105 = conj(IT_0063)*IT_0104;
    const ccomplex_t IT_0106 = pow(s_34, 2);
    const ccomplex_t IT_0107 = IT_0019*IT_0106;
    const ccomplex_t IT_0108 = (-3)*s_13*IT_0107;
    const ccomplex_t IT_0109 = (-3)*s_23*IT_0107;
    const ccomplex_t IT_0110 = (-3)*IT_0107;
    const ccomplex_t IT_0111 = 0.333333333333333*conj(IT_0012)*IT_0041 + conj
      (IT_0063) + 0.111111111111111*(conj(IT_0065) + conj(IT_0068))*IT_0086 +
       0.111111111111111*(conj(IT_0036) + conj(IT_0070))*IT_0097 +
       0.111111111111111*conj(IT_0018)*IT_0102 + 0.111111111111111*(conj(IT_0025
      ) + conj(IT_0052))*IT_0104 + 0.111111111111111*conj(IT_0054)*IT_0110;
    const ccomplex_t IT_0112 = (-0.333333333333333)*conj(IT_0065);
    const ccomplex_t IT_0113 = (-0.333333333333333)*conj(IT_0036);
    const ccomplex_t IT_0114 = (-0.333333333333333)*conj(IT_0068);
    const ccomplex_t IT_0115 = (-0.333333333333333)*conj(IT_0018);
    const ccomplex_t IT_0116 = (-0.333333333333333)*conj(IT_0052);
    const ccomplex_t IT_0117 = (-0.333333333333333)*conj(IT_0054);
    const ccomplex_t IT_0118 = (-0.333333333333333)*conj(IT_0012)*IT_0071 + (
      -0.333333333333333)*IT_0098 + (s_12*conj(IT_0025) + IT_0037*conj(IT_0070))
      *IT_0107 + IT_0029*IT_0112 + IT_0055*IT_0113 + IT_0085*IT_0114 + IT_0082
      *IT_0115 + IT_0084*IT_0116 + IT_0108*IT_0117;
    const ccomplex_t IT_0119 = (-0.333333333333333)*conj(IT_0012)*IT_0029 + (
      -0.333333333333333)*IT_0105 + (s_12*conj(IT_0070) + conj(IT_0025)*IT_0073)
      *IT_0107 + IT_0095*IT_0112 + IT_0092*IT_0113 + IT_0082*IT_0114 + IT_0101
      *IT_0115 + IT_0100*IT_0116 + IT_0109*IT_0117;
    const ccomplex_t IT_0120 = IT_0012*IT_0072 + IT_0088 + IT_0096 + IT_0036*
      (conj(IT_0012)*IT_0043 + IT_0053*conj(IT_0068) + IT_0055*conj(IT_0070) +
       conj(IT_0065)*IT_0090 + conj(IT_0018)*IT_0091 + conj(IT_0025)*IT_0092 +
       IT_0098) + IT_0018*IT_0103 + IT_0052*(conj(IT_0012)*IT_0053 + conj
      (IT_0068)*IT_0083 + conj(IT_0070)*IT_0084 + conj(IT_0065)*IT_0091 + conj
      (IT_0018)*IT_0099 + conj(IT_0025)*IT_0100 + IT_0105) + IT_0054*(conj
      (IT_0012)*IT_0055 + conj(IT_0068)*IT_0084 + conj(IT_0065)*IT_0092 + conj
      (IT_0018)*IT_0100 + conj(IT_0070)*IT_0108 + conj(IT_0025)*IT_0109 + conj
      (IT_0063)*IT_0110) + 9*IT_0063*IT_0111 + (-3)*IT_0070*IT_0118 + (-3)
      *IT_0025*IT_0119;
    return create_ccomplex_return(IT_0120);
}

