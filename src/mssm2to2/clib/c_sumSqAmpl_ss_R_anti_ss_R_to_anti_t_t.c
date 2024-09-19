#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = cpow(IT_0002, -1);
    const ccomplex_t IT_0010 = cpow(IT_0006, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0003 + -IT_0004*IT_0005)*IT_0006 + (-3)*IT_0003*IT_0007*IT_0008
      *IT_0009*IT_0010);
    const ccomplex_t IT_0012 = 0.333333333333333*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0004, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_t*e_em*IT_0005*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = pow(m_ss_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0003
      *IT_0004 + IT_0002*IT_0005)*IT_0006 + (-3)*IT_0005*IT_0007*IT_0008*IT_0009
      *IT_0010);
    const ccomplex_t IT_0022 = (-0.333333333333333)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_t*e_em*IT_0003*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_t, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = -IT_0038 + -IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0000, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0048;
    const ccomplex_t IT_0050 = (-0.666666666666667)*IT_0049;
    const ccomplex_t IT_0051 = IT_0006*IT_0048;
    const ccomplex_t IT_0052 = e_em*IT_0051;
    const ccomplex_t IT_0053 = IT_0000*IT_0010;
    const ccomplex_t IT_0054 = e_em*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0052 + (-3)*IT_0054);
    const ccomplex_t IT_0056 = (-0.166666666666667)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0050*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0050*IT_0057;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0050*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0017*IT_0030;
    const ccomplex_t IT_0067 = m_t*IT_0050;
    const ccomplex_t IT_0068 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0069 = m_t*IT_0068;
    const ccomplex_t IT_0070 = (IT_0059 + -IT_0063)*(IT_0067 + -IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = (-18)*IT_0030;
    const ccomplex_t IT_0073 = 18*s_34;
    const ccomplex_t IT_0074 = IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = s_34*IT_0017;
    const ccomplex_t IT_0076 = (-6)*IT_0075;
    const ccomplex_t IT_0077 = s_13*s_14;
    const ccomplex_t IT_0078 = 12*IT_0077;
    const ccomplex_t IT_0079 = IT_0076 + IT_0078;
    const ccomplex_t IT_0080 = (-18)*IT_0075;
    const ccomplex_t IT_0081 = 36*IT_0077;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = (-18)*IT_0071;
    const ccomplex_t IT_0084 = IT_0038 + IT_0046;
    const ccomplex_t IT_0085 = m_t*IT_0084;
    const ccomplex_t IT_0086 = (-6)*IT_0071;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0088 = (-2)*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0088*IT_0091;
    const ccomplex_t IT_0093 = 0.25*IT_0092;
    const ccomplex_t IT_0094 = (-12)*IT_0075;
    const ccomplex_t IT_0095 = (-12)*IT_0066;
    const ccomplex_t IT_0096 = 24*IT_0077;
    const ccomplex_t IT_0097 = IT_0094 + IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0099 = 0.666666666666667*IT_0098;
    const ccomplex_t IT_0100 = cpow(IT_0099, 2);
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0089*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (-0.0833333333333333)*IT_0092;
    const ccomplex_t IT_0105 = (-2)*IT_0104;
    const ccomplex_t IT_0106 = (-36)*IT_0075;
    const ccomplex_t IT_0107 = (-36)*IT_0066;
    const ccomplex_t IT_0108 = 72*IT_0077;
    const ccomplex_t IT_0109 = IT_0106 + IT_0107 + IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0104, 2);
    const ccomplex_t IT_0111 = cpow(IT_0093, 2);
    const ccomplex_t IT_0112 = (-18)*IT_0066;
    const ccomplex_t IT_0113 = IT_0080 + IT_0081 + IT_0112;
    const ccomplex_t IT_0114 = IT_0103 + IT_0104;
    const ccomplex_t IT_0115 = conj(IT_0061) + conj(IT_0065);
    const ccomplex_t IT_0116 = IT_0061 + IT_0065;
    const ccomplex_t IT_0117 = -IT_0104;
    const ccomplex_t IT_0118 = -IT_0103;
    const ccomplex_t IT_0119 = IT_0117 + IT_0118;
    const ccomplex_t IT_0120 = -IT_0093;
    const ccomplex_t IT_0121 = (-6)*IT_0066;
    const ccomplex_t IT_0122 = IT_0076 + IT_0078 + IT_0121;
    const ccomplex_t IT_0123 = 6*IT_0029;
    const ccomplex_t IT_0124 = -IT_0030;
    const ccomplex_t IT_0125 = m_t*s_13;
    const ccomplex_t IT_0126 = 18*conj(IT_0071);
    const ccomplex_t IT_0127 = 6*IT_0085;
    const ccomplex_t IT_0128 = 6*conj(IT_0085);
    const ccomplex_t IT_0129 = 18*conj(IT_0070);
    const ccomplex_t IT_0130 = 18*conj(IT_0065)*IT_0071 + 6*conj(IT_0047)*
      (IT_0070 + 3*IT_0085) + 6*IT_0047*(conj(IT_0070) + 3*conj(IT_0085)) +
       IT_0065*IT_0126 + conj(IT_0061)*(18*IT_0070 + IT_0127) + IT_0061*(IT_0128
       + IT_0129);
    const ccomplex_t IT_0131 = m_t*s_14;
    const ccomplex_t IT_0132 = (-18)*conj(IT_0071);
    const ccomplex_t IT_0133 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0134 = (-6)*IT_0085;
    const ccomplex_t IT_0135 = -s_14;
    const ccomplex_t IT_0136 = s_13 + IT_0135;
    const ccomplex_t IT_0137 = m_t*IT_0136;
    const ccomplex_t IT_0138 = 6*conj(IT_0029);
    const ccomplex_t IT_0139 = 18*IT_0029;
    const ccomplex_t IT_0140 = 18*conj(IT_0029);
    const ccomplex_t IT_0141 = 36*conj(IT_0029);
    const ccomplex_t IT_0142 = ((-12)*IT_0029 + 12*conj(IT_0029) + (-6)
      *IT_0070 + 6*conj(IT_0070) + 6*conj(IT_0071) + (-18)*IT_0085 + 18*conj
      (IT_0085) + IT_0086)*IT_0093 + conj(IT_0047)*IT_0123 + IT_0047*IT_0138 +
       IT_0115*IT_0139 + IT_0116*IT_0140 + IT_0114*((-36)*IT_0029 + (-18)
      *IT_0070 + IT_0083 + IT_0126 + IT_0128 + IT_0129 + IT_0134 + IT_0141);
    const ccomplex_t IT_0143 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0144 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + (-6)*((conj(IT_0047) + 3*conj(IT_0061))*IT_0065 + (IT_0047 + 3*IT_0061)
      *conj(IT_0065))*IT_0066 + IT_0071*conj(IT_0071)*IT_0073 + (conj(IT_0029)*
      (IT_0070 + IT_0071) + IT_0029*(conj(IT_0070) + conj(IT_0071)))*IT_0074 + 
      (conj(IT_0047)*IT_0061 + IT_0047*conj(IT_0061))*IT_0079 + (IT_0047*conj
      (IT_0047) + IT_0061*conj(IT_0061) + IT_0065*conj(IT_0065))*IT_0082 +
       IT_0030*(conj(IT_0070)*IT_0083 + conj(IT_0085)*IT_0086) + (-2)*IT_0093
      *IT_0097*(IT_0103 + (-0.5)*IT_0105) + -IT_0109*(IT_0103*(IT_0103 + 
      -IT_0105) + IT_0110 + IT_0111) + IT_0113*(conj(IT_0047)*IT_0093 + IT_0114
      *IT_0115 + IT_0116*IT_0119 + IT_0047*IT_0120) + (conj(IT_0047)*IT_0114 +
       IT_0093*IT_0115 + IT_0047*IT_0119 + IT_0116*IT_0120)*IT_0122 + conj
      (IT_0085)*IT_0123*(s_34 + IT_0124) + IT_0125*IT_0130 + IT_0131*(conj
      (IT_0061)*IT_0083 + IT_0065*((-18)*conj(IT_0070) + (-6)*conj(IT_0085)) +
       conj(IT_0047)*IT_0086 + IT_0061*IT_0132 + IT_0047*IT_0133 + conj(IT_0065)
      *((-18)*IT_0070 + IT_0134)) + IT_0137*IT_0142 + 6*IT_0070*
      (0.166666666666667*conj(IT_0070)*IT_0073 + s_34*conj(IT_0085) + IT_0132
      *IT_0143) + IT_0127*(s_34*conj(IT_0070) + 0.166666666666667*IT_0073*conj
      (IT_0085) + 0.166666666666667*(s_34 + IT_0124)*IT_0138 + IT_0133*IT_0143);
    return create_ccomplex_return(IT_0144);
}

