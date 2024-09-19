#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_u, 2);
    const ccomplex_t IT_0013 = cpow(IT_0010, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + 1.5*IT_0000*IT_0005*IT_0011*IT_0012
      *IT_0013);
    const ccomplex_t IT_0015 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = pow(m_su_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + (-1.5)*IT_0000*IT_0003*IT_0005
      *IT_0012*IT_0013);
    const ccomplex_t IT_0024 = 0.666666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = cpow(IT_0008, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0031;
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = m_b*IT_0033;
    const ccomplex_t IT_0035 = 1.33333333333333*IT_0032;
    const ccomplex_t IT_0036 = m_b*IT_0035;
    const ccomplex_t IT_0037 = IT_0004*IT_0031;
    const ccomplex_t IT_0038 = e_em*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0008;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0040);
    const ccomplex_t IT_0042 = (-0.166666666666667)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.333333333333333*IT_0032;
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (IT_0034 + -IT_0036)*(IT_0045 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = (-18)*IT_0030;
    const ccomplex_t IT_0052 = 18*s_34;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = IT_0035*IT_0045;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0000*IT_0005*IT_0013*IT_0056*V_ub_mod;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0000*IT_0005*IT_0013*IT_0059*V_ub_mod;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0058*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0000*IT_0005*IT_0013*IT_0056*V_ub_mod;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0000*IT_0005*IT_0013*IT_0059*V_ub_mod;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0065 + IT_0073;
    const ccomplex_t IT_0075 = s_34*IT_0017;
    const ccomplex_t IT_0076 = (-6)*IT_0075;
    const ccomplex_t IT_0077 = s_13*s_14;
    const ccomplex_t IT_0078 = 12*IT_0077;
    const ccomplex_t IT_0079 = IT_0076 + IT_0078;
    const ccomplex_t IT_0080 = IT_0035*IT_0048;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = (-18)*IT_0075;
    const ccomplex_t IT_0083 = 36*IT_0077;
    const ccomplex_t IT_0084 = IT_0082 + IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0086 = (-0.333333333333333)*IT_0085;
    const ccomplex_t IT_0087 = (-1.33333333333333)*IT_0085;
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0093 = (-2)*IT_0092;
    const ccomplex_t IT_0094 = IT_0089*IT_0092;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0091 + (-0.0833333333333333)*IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0097, 2);
    const ccomplex_t IT_0099 = (-36)*IT_0075;
    const ccomplex_t IT_0100 = IT_0017*IT_0030;
    const ccomplex_t IT_0101 = (-36)*IT_0100;
    const ccomplex_t IT_0102 = 72*IT_0077;
    const ccomplex_t IT_0103 = IT_0099 + IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = 0.25*IT_0096;
    const ccomplex_t IT_0105 = (-18)*IT_0050;
    const ccomplex_t IT_0106 = m_b*IT_0074;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = (-6)*IT_0107;
    const ccomplex_t IT_0109 = IT_0105 + IT_0108;
    const ccomplex_t IT_0110 = (-18)*conj(IT_0050);
    const ccomplex_t IT_0111 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = (-18)*IT_0055;
    const ccomplex_t IT_0114 = (-6)*IT_0074;
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = (-18)*conj(IT_0055);
    const ccomplex_t IT_0117 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = m_b*s_14;
    const ccomplex_t IT_0120 = m_b*s_13;
    const ccomplex_t IT_0121 = 6*conj(IT_0050);
    const ccomplex_t IT_0122 = 18*conj(IT_0107);
    const ccomplex_t IT_0123 = 6*conj(IT_0107);
    const ccomplex_t IT_0124 = 18*conj(IT_0050);
    const ccomplex_t IT_0125 = 3*IT_0050*(conj(IT_0055) + 0.333333333333333
      *conj(IT_0074)) + 3*conj(IT_0049)*IT_0081 + 3*IT_0049*conj(IT_0081) + 
      (conj(IT_0055) + 3*conj(IT_0074))*IT_0107 + 0.166666666666667*IT_0074*
      (IT_0121 + IT_0122) + 0.166666666666667*IT_0055*(IT_0123 + IT_0124);
    const ccomplex_t IT_0126 = -s_14;
    const ccomplex_t IT_0127 = s_13 + IT_0126;
    const ccomplex_t IT_0128 = m_b*IT_0127;
    const ccomplex_t IT_0129 = 6*IT_0029;
    const ccomplex_t IT_0130 = 6*conj(IT_0029);
    const ccomplex_t IT_0131 = conj(IT_0055) + conj(IT_0081);
    const ccomplex_t IT_0132 = 18*IT_0029;
    const ccomplex_t IT_0133 = IT_0055 + IT_0081;
    const ccomplex_t IT_0134 = 18*conj(IT_0029);
    const ccomplex_t IT_0135 = IT_0104*((-12)*IT_0029 + 12*conj(IT_0029) + (-6
      )*IT_0049 + 6*conj(IT_0049) + (-6)*IT_0050 + (-18)*IT_0107 + IT_0121 +
       IT_0122) + IT_0097*((-36)*IT_0029 + 36*conj(IT_0029) + (-18)*IT_0049 + 18
      *conj(IT_0049) + IT_0105 + IT_0108 + IT_0123 + IT_0124) + conj(IT_0074)
      *IT_0129 + IT_0074*IT_0130 + IT_0131*IT_0132 + IT_0133*IT_0134;
    const ccomplex_t IT_0136 = (-6)*IT_0100;
    const ccomplex_t IT_0137 = IT_0076 + IT_0078 + IT_0136;
    const ccomplex_t IT_0138 = -IT_0104;
    const ccomplex_t IT_0139 = -IT_0097;
    const ccomplex_t IT_0140 = (-18)*IT_0100;
    const ccomplex_t IT_0141 = IT_0082 + IT_0083 + IT_0140;
    const ccomplex_t IT_0142 = -IT_0030;
    const ccomplex_t IT_0143 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + (IT_0049*conj(IT_0049) + IT_0050*conj(IT_0050))*IT_0052 + (conj(IT_0029
      )*(IT_0049 + IT_0050) + IT_0029*(conj(IT_0049) + conj(IT_0050)))*IT_0053 +
       (conj(IT_0055)*IT_0074 + IT_0055*conj(IT_0074))*IT_0079 + (IT_0055*conj
      (IT_0055) + IT_0074*conj(IT_0074) + IT_0081*conj(IT_0081))*IT_0084 + 
      -IT_0098*IT_0103 + (-2)*IT_0104*(IT_0097*((-12)*IT_0075 + 24*IT_0077 + (
      -12)*IT_0100) + 0.5*IT_0103*IT_0104) + IT_0030*(conj(IT_0049)*IT_0109 +
       IT_0049*IT_0112) + IT_0100*(conj(IT_0081)*IT_0115 + IT_0081*IT_0118) + 
      (conj(IT_0081)*IT_0109 + IT_0081*IT_0112 + conj(IT_0049)*IT_0115 + IT_0049
      *IT_0118)*IT_0119 + 6*IT_0120*IT_0125 + IT_0128*IT_0135 + IT_0137*(conj
      (IT_0074)*IT_0097 + IT_0104*IT_0131 + IT_0133*IT_0138 + IT_0074*IT_0139) +
       (conj(IT_0074)*IT_0104 + IT_0097*IT_0131 + IT_0074*IT_0138 + IT_0133
      *IT_0139)*IT_0141 + 6*conj(IT_0107)*(s_34*IT_0050 + 0.166666666666667
      *IT_0052*IT_0107 + 0.166666666666667*IT_0129*(s_34 + IT_0142)) + 6*IT_0107
      *(s_34*conj(IT_0050) + 0.166666666666667*IT_0130*(s_34 + IT_0142));
    return create_ccomplex_return(IT_0143);
}

