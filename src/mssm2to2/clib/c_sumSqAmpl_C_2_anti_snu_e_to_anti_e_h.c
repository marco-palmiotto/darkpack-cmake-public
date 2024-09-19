#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_e_to_anti_e_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_e_to_anti_e_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_h, 2);
    const ccomplex_t IT_0001 = m_e*m_C_2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_e, 2);
    const ccomplex_t IT_0003 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_23 + IT_0002 + IT_0003 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = V_u2*U_Wm1;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = U_d1*V_Wp2;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (IT_0007 + -IT_0010)*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0012*IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0002 + IT_0003 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0026 = V_u2*U_Wm2;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = U_d2*V_Wp2;
    const ccomplex_t IT_0029 = IT_0008*IT_0028;
    const ccomplex_t IT_0030 = IT_0012*(IT_0027 + -IT_0029);
    const ccomplex_t IT_0031 = 1.4142135623731*e_em*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0012*IT_0017*IT_0019;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0024 + IT_0038;
    const ccomplex_t IT_0040 = pow(m_C_2, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + IT_0002 + -IT_0003 + -IT_0040 
      + -reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_e*e_em*IT_0008*IT_0012
      *IT_0017*IT_0019;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0012;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = conj(IT_0039) + conj(IT_0048);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0012;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0053 = IT_0005*IT_0052;
    const ccomplex_t IT_0054 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0055 = IT_0008*IT_0054;
    const ccomplex_t IT_0056 = IT_0012*(IT_0053 + -IT_0055);
    const ccomplex_t IT_0057 = 1.4142135623731*e_em*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0051*IT_0059;
    const ccomplex_t IT_0061 = IT_0004*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0064 = IT_0005*IT_0063;
    const ccomplex_t IT_0065 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0066 = IT_0008*IT_0065;
    const ccomplex_t IT_0067 = IT_0012*(IT_0064 + -IT_0066);
    const ccomplex_t IT_0068 = 1.4142135623731*e_em*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0045*IT_0070;
    const ccomplex_t IT_0072 = IT_0025*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = -IT_0062 + -IT_0073;
    const ccomplex_t IT_0075 = IT_0035*IT_0043;
    const ccomplex_t IT_0076 = IT_0041*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0074 + IT_0078;
    const ccomplex_t IT_0080 = IT_0039 + IT_0048;
    const ccomplex_t IT_0081 = conj(IT_0074) + conj(IT_0078);
    const ccomplex_t IT_0082 = m_C_2*IT_0039;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0002 + IT_0040 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0085 = cos(theta_W);
    const ccomplex_t IT_0086 = cpow(IT_0085, -2);
    const ccomplex_t IT_0087 = IT_0011*IT_0086;
    const ccomplex_t IT_0088 = IT_0012 + IT_0087;
    const ccomplex_t IT_0089 = sin(beta);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_W*e_em*IT_0088*(IT_0008
      *IT_0018 + IT_0005*IT_0089);
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0035*IT_0091;
    const ccomplex_t IT_0093 = IT_0084*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0096 = IT_0021*IT_0059;
    const ccomplex_t IT_0097 = IT_0004*IT_0095*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0099 = IT_0035*IT_0070;
    const ccomplex_t IT_0100 = IT_0025*IT_0098*IT_0099;
    const ccomplex_t IT_0101 = m_e*IT_0077;
    const ccomplex_t IT_0102 = -IT_0094 + IT_0097 + IT_0100 + 2*IT_0101;
    const ccomplex_t IT_0103 = IT_0033*IT_0045;
    const ccomplex_t IT_0104 = IT_0025*IT_0098*IT_0103;
    const ccomplex_t IT_0105 = IT_0045*IT_0091;
    const ccomplex_t IT_0106 = IT_0084*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = m_e*IT_0048;
    const ccomplex_t IT_0109 = IT_0016*IT_0051;
    const ccomplex_t IT_0110 = IT_0004*IT_0095*IT_0109;
    const ccomplex_t IT_0111 = -IT_0104 + IT_0107 + (-2)*IT_0108 + -IT_0110;
    const ccomplex_t IT_0112 = IT_0062 + IT_0073;
    const ccomplex_t IT_0113 = m_C_2*IT_0112;
    const ccomplex_t IT_0114 = 2*s_13;
    const ccomplex_t IT_0115 = (-2)*s_13;
    const ccomplex_t IT_0116 = m_e*m_C_2;
    const ccomplex_t IT_0117 = (-2)*IT_0111;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0113 + IT_0118;
    const ccomplex_t IT_0120 = (-2)*conj(IT_0111);
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = conj(IT_0113) + IT_0121;
    const ccomplex_t IT_0123 = (-2)*IT_0102;
    const ccomplex_t IT_0124 = (-2)*conj(IT_0102);
    const ccomplex_t IT_0125 = s_13*IT_0000;
    const ccomplex_t IT_0126 = s_14*s_34;
    const ccomplex_t IT_0127 = (-2)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = (-2)*IT_0039;
    const ccomplex_t IT_0130 = (-2)*IT_0074;
    const ccomplex_t IT_0131 = s_34*m_C_2;
    const ccomplex_t IT_0132 = 0.5*IT_0124;
    const ccomplex_t IT_0133 = conj(IT_0083) + IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0123;
    const ccomplex_t IT_0135 = m_e*s_14;
    const ccomplex_t IT_0136 = 2*IT_0001*(IT_0049*IT_0079 + IT_0080*IT_0081) +
       (IT_0083*conj(IT_0083) + IT_0102*conj(IT_0102) + IT_0111*conj(IT_0111) +
       IT_0113*conj(IT_0113))*IT_0114 + (conj(IT_0083)*IT_0102 + IT_0083*conj
      (IT_0102) + conj(IT_0111)*IT_0113 + IT_0111*conj(IT_0113))*IT_0115 + 2
      *IT_0116*(conj(IT_0102)*IT_0111 + IT_0102*conj(IT_0111) + conj(IT_0083)
      *IT_0119 + IT_0083*IT_0122 + 0.5*conj(IT_0113)*IT_0123 + 0.5*IT_0113
      *IT_0124) + (-2)*IT_0128*(IT_0049*(IT_0048 + (-0.5)*IT_0129) + IT_0081*
      (IT_0078 + (-0.5)*IT_0130)) + 2*IT_0131*(IT_0081*IT_0119 + IT_0079*IT_0122
       + IT_0080*IT_0133 + IT_0049*(IT_0083 + IT_0134)) + 2*(IT_0049*IT_0119 +
       IT_0080*IT_0122 + IT_0079*IT_0133 + IT_0081*(IT_0083 + IT_0134))*IT_0135;
    return create_ccomplex_return(IT_0136);
}

