#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_se_R_to_nu_e_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_se_R_to_nu_e_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0005 + -IT_0008 + (-2)
      *IT_0010);
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0006*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_se_R, 2);
    const ccomplex_t IT_0020 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = IT_0007*IT_0023;
    const ccomplex_t IT_0026 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0027 = IT_0007*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0025 + (-2)
      *IT_0027);
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0006*IT_0013*IT_0015;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0037 = (-2)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_C_2, 2);
    const ccomplex_t IT_0039 = cpow(s_13 + (-0.5)*IT_0038 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0037*IT_0041;
    const ccomplex_t IT_0043 = -IT_0022 + -IT_0035 + 0.5*IT_0042;
    const ccomplex_t IT_0044 = pow(m_Z, -2);
    const ccomplex_t IT_0045 = pow(s_14, 2);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = -IT_0038;
    const ccomplex_t IT_0048 = IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = s_13*IT_0048;
    const ccomplex_t IT_0050 = cpow(s_12 + 0.5*IT_0019 + 0.5*IT_0038 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0051 = e_em*IT_0003;
    const ccomplex_t IT_0052 = e_em*IT_0007;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0051 + IT_0052);
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0031*IT_0054;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0042 + -IT_0057;
    const ccomplex_t IT_0059 = pow(s_34, 2);
    const ccomplex_t IT_0060 = s_13*IT_0044*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0062 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = IT_0007*IT_0062;
    const ccomplex_t IT_0065 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0066 = IT_0007*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + -IT_0064 + (-2)
      *IT_0066);
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0017*IT_0068;
    const ccomplex_t IT_0070 = IT_0020*IT_0061*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0072 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0073 = IT_0003*IT_0072;
    const ccomplex_t IT_0074 = IT_0007*IT_0072;
    const ccomplex_t IT_0075 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0076 = IT_0007*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + -IT_0074 + (-2)
      *IT_0076);
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0031*IT_0078;
    const ccomplex_t IT_0080 = IT_0033*IT_0071*IT_0079;
    const ccomplex_t IT_0081 = m_C_2*IT_0035;
    const ccomplex_t IT_0082 = m_C_2*IT_0022;
    const ccomplex_t IT_0083 = 0.5*IT_0070 + 0.5*IT_0080 + 0.5*IT_0081 + 0.5
      *IT_0082;
    const ccomplex_t IT_0084 = m_C_2*IT_0044*IT_0059;
    const ccomplex_t IT_0085 = (-0.5)*IT_0057;
    const ccomplex_t IT_0086 = s_13*s_34;
    const ccomplex_t IT_0087 = 2*IT_0058;
    const ccomplex_t IT_0088 = pow(m_Z, 2);
    const ccomplex_t IT_0089 = s_34*IT_0088;
    const ccomplex_t IT_0090 = IT_0044*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = s_34 + IT_0091;
    const ccomplex_t IT_0093 = m_C_2*IT_0092;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = 2*IT_0083;
    const ccomplex_t IT_0096 = s_14*s_34*IT_0044;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = s_13 + IT_0097;
    const ccomplex_t IT_0099 = s_13*IT_0098;
    const ccomplex_t IT_0100 = m_C_2*IT_0098;
    const ccomplex_t IT_0101 = s_13*s_14;
    const ccomplex_t IT_0102 = s_34*IT_0038;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = (-2)*IT_0043;
    const ccomplex_t IT_0106 = (-2)*conj(IT_0043);
    const ccomplex_t IT_0107 = s_13*IT_0088;
    const ccomplex_t IT_0108 = s_14*s_34;
    const ccomplex_t IT_0109 = (-2)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107 + IT_0109;
    const ccomplex_t IT_0111 = IT_0044*IT_0110;
    const ccomplex_t IT_0112 = 2*conj(IT_0058);
    const ccomplex_t IT_0113 = pow(m_Z, 4);
    const ccomplex_t IT_0114 = s_13*IT_0113;
    const ccomplex_t IT_0115 = IT_0044*IT_0114;
    const ccomplex_t IT_0116 = IT_0107 + IT_0109 + (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = 2*conj(IT_0083);
    const ccomplex_t IT_0118 = (-0.25)*IT_0115;
    const ccomplex_t IT_0119 = IT_0108 + IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0022 + 0.5*IT_0035;
    const ccomplex_t IT_0121 = 2*conj(IT_0120);
    const ccomplex_t IT_0122 = (-0.25)*IT_0090;
    const ccomplex_t IT_0123 = s_34 + IT_0122;
    const ccomplex_t IT_0124 = m_C_2*IT_0123;
    const ccomplex_t IT_0125 = (-4)*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = 2*IT_0043*conj(IT_0043)*IT_0049 + (conj(IT_0058
      )*IT_0060 + conj(IT_0083)*IT_0084 + conj(IT_0085)*IT_0086)*IT_0087 + conj
      (IT_0085)*IT_0094*IT_0095 + (conj(IT_0058)*IT_0099 + conj(IT_0083)*IT_0100
       + conj(IT_0085)*IT_0104)*IT_0105 + IT_0058*IT_0099*IT_0106 + 4*IT_0083*
      (0.25*IT_0100*IT_0106 + conj(IT_0083)*(s_13 + (-0.5)*IT_0111) + 0.25
      *IT_0084*IT_0112) + (-4)*IT_0085*((-0.25)*IT_0104*IT_0106 + (-0.25)
      *IT_0086*IT_0112 + conj(IT_0085)*IT_0116 + (-0.25)*IT_0094*IT_0117 +
       IT_0119*IT_0121) + (-2)*conj(IT_0120)*(IT_0058*IT_0086 + -IT_0043*IT_0104
       + IT_0095*IT_0126) + (-2)*IT_0120*(conj(IT_0058)*IT_0086 + -conj(IT_0043)
      *IT_0104 + 4*conj(IT_0085)*IT_0119 + IT_0116*IT_0121 + IT_0117*IT_0126);
    return create_ccomplex_return(IT_0127);
}

