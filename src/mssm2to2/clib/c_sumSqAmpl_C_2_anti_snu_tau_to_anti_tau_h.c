#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = m_C_2*m_tau;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0002;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = IT_0002*(IT_0007 + -IT_0010);
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0018 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0022 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0005*IT_0022;
    const ccomplex_t IT_0024 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = IT_0002*(IT_0023 + -IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0002;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = IT_0021*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0020 + IT_0034;
    const ccomplex_t IT_0036 = m_C_2*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = V_u2*U_Wm1;
    const ccomplex_t IT_0039 = IT_0005*IT_0038;
    const ccomplex_t IT_0040 = U_d1*V_Wp2;
    const ccomplex_t IT_0041 = IT_0008*IT_0040;
    const ccomplex_t IT_0042 = IT_0002*(IT_0039 + -IT_0041);
    const ccomplex_t IT_0043 = 1.4142135623731*e_em*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0004*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0048 = IT_0018*IT_0046*IT_0047;
    const ccomplex_t IT_0049 = pow(m_C_2, 2);
    const ccomplex_t IT_0050 = cpow(s_12 + (-0.5)*IT_0016 + 0.5*IT_0017 + 0.5
      *IT_0049 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0051 = pow(m_W, -1);
    const ccomplex_t IT_0052 = cos(beta);
    const ccomplex_t IT_0053 = cpow(IT_0052, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*m_tau*IT_0002*IT_0008
      *IT_0051*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0031*IT_0055;
    const ccomplex_t IT_0057 = IT_0050*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = m_tau*IT_0058;
    const ccomplex_t IT_0060 = cpow(s_13 + (-0.5)*IT_0016 + (-0.5)*IT_0049 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0061 = cos(theta_W);
    const ccomplex_t IT_0062 = cpow(IT_0061, -2);
    const ccomplex_t IT_0063 = IT_0001*IT_0062;
    const ccomplex_t IT_0064 = IT_0002 + IT_0063;
    const ccomplex_t IT_0065 = sin(beta);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_W*e_em*IT_0064*(IT_0008
      *IT_0052 + IT_0005*IT_0065);
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0031*IT_0067;
    const ccomplex_t IT_0069 = IT_0060*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0072 = V_u2*U_Wm2;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = U_d2*V_Wp2;
    const ccomplex_t IT_0075 = IT_0008*IT_0074;
    const ccomplex_t IT_0076 = IT_0002*(IT_0073 + -IT_0075);
    const ccomplex_t IT_0077 = 1.4142135623731*e_em*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0031*IT_0079;
    const ccomplex_t IT_0081 = IT_0021*IT_0071*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0048 + IT_0059 + (-0.5)*IT_0070 + 0.5
      *IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0002*IT_0051*IT_0053;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0014*IT_0084;
    const ccomplex_t IT_0086 = IT_0018*IT_0047*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0002*IT_0051*IT_0053;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0029*IT_0088;
    const ccomplex_t IT_0090 = IT_0021*IT_0071*IT_0089;
    const ccomplex_t IT_0091 = IT_0055*IT_0088;
    const ccomplex_t IT_0092 = IT_0050*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = m_tau*IT_0093;
    const ccomplex_t IT_0095 = IT_0067*IT_0088;
    const ccomplex_t IT_0096 = IT_0060*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0086 + (-0.5)*IT_0090 + -IT_0094 +
       0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0045*IT_0084;
    const ccomplex_t IT_0100 = IT_0018*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0079*IT_0088;
    const ccomplex_t IT_0103 = IT_0021*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0101 + IT_0104;
    const ccomplex_t IT_0106 = m_C_2*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = 2*IT_0037;
    const ccomplex_t IT_0109 = (-2)*IT_0082;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = 2*conj(IT_0037);
    const ccomplex_t IT_0112 = (-2)*conj(IT_0082);
    const ccomplex_t IT_0113 = IT_0111 + IT_0112;
    const ccomplex_t IT_0114 = s_14*m_tau;
    const ccomplex_t IT_0115 = (-0.5)*IT_0058;
    const ccomplex_t IT_0116 = (-0.5)*IT_0101 + (-0.5)*IT_0104;
    const ccomplex_t IT_0117 = conj(IT_0115) + conj(IT_0116);
    const ccomplex_t IT_0118 = IT_0115 + IT_0116;
    const ccomplex_t IT_0119 = 0.5*IT_0093;
    const ccomplex_t IT_0120 = 0.5*IT_0020 + 0.5*IT_0034;
    const ccomplex_t IT_0121 = conj(IT_0119) + conj(IT_0120);
    const ccomplex_t IT_0122 = 2*IT_0107;
    const ccomplex_t IT_0123 = (-2)*IT_0098;
    const ccomplex_t IT_0124 = IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = IT_0119 + IT_0120;
    const ccomplex_t IT_0126 = 2*conj(IT_0107);
    const ccomplex_t IT_0127 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0128 = IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = s_34*m_C_2;
    const ccomplex_t IT_0130 = 2*s_13;
    const ccomplex_t IT_0131 = pow(m_h, 2);
    const ccomplex_t IT_0132 = m_C_2*m_tau*IT_0131;
    const ccomplex_t IT_0133 = 2*IT_0119;
    const ccomplex_t IT_0134 = 2*conj(IT_0119);
    const ccomplex_t IT_0135 = (-2)*s_13;
    const ccomplex_t IT_0136 = s_13*IT_0131;
    const ccomplex_t IT_0137 = s_14*s_34;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = (-2)*IT_0119;
    const ccomplex_t IT_0141 = (-2)*IT_0115;
    const ccomplex_t IT_0142 = (-2)*IT_0000*((conj(IT_0037) + -conj(IT_0082))
      *IT_0098 + (IT_0037 + -IT_0082)*conj(IT_0098) + (-0.5)*conj(IT_0107)
      *IT_0110 + (-0.5)*IT_0107*IT_0113) + IT_0114*(IT_0110*IT_0117 + IT_0113
      *IT_0118 + IT_0121*IT_0124 + IT_0125*IT_0128) + (IT_0110*IT_0121 + IT_0117
      *IT_0124 + IT_0113*IT_0125 + IT_0118*IT_0128)*IT_0129 + (IT_0037*conj
      (IT_0037) + IT_0082*conj(IT_0082) + IT_0098*conj(IT_0098) + IT_0107*conj
      (IT_0107))*IT_0130 + 2*IT_0132*(IT_0117*(IT_0120 + 0.5*IT_0133) + IT_0118*
      (conj(IT_0120) + 0.5*IT_0134)) + (conj(IT_0037)*IT_0082 + IT_0037*conj
      (IT_0082) + conj(IT_0098)*IT_0107 + IT_0098*conj(IT_0107))*IT_0135 + (-2)
      *IT_0139*(IT_0121*(IT_0120 + (-0.5)*IT_0140) + IT_0117*(IT_0116 + (-0.5)
      *IT_0141));
    return create_ccomplex_return(IT_0142);
}

