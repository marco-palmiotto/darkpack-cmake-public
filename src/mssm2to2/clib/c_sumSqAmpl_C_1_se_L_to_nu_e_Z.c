#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_se_L_to_nu_e_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_se_L_to_nu_e_Z(
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
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0003;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0007 + -IT_0009 + (-2)
      *IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0016 = pow(m_se_L, 2);
    const ccomplex_t IT_0017 = cpow(s_23 + (-0.5)*IT_0016 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0020 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0020;
    const ccomplex_t IT_0023 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + -IT_0022 + (-2)
      *IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0029 = cpow(s_23 + (-0.5)*IT_0016 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0028*IT_0029;
    const ccomplex_t IT_0031 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0032 = IT_0005*IT_0031;
    const ccomplex_t IT_0033 = IT_0008*IT_0031;
    const ccomplex_t IT_0034 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0035 = IT_0008*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + -IT_0033 + (-2)
      *IT_0035);
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0019*IT_0037;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = m_C_1*IT_0040;
    const ccomplex_t IT_0042 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0043 = IT_0005*IT_0042;
    const ccomplex_t IT_0044 = IT_0008*IT_0042;
    const ccomplex_t IT_0045 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0046 = IT_0008*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + -IT_0044 + (-2)
      *IT_0046);
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0002*IT_0048;
    const ccomplex_t IT_0050 = IT_0017*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = m_C_1*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0018 + 0.5*IT_0030 + 0.5*IT_0041 + 0.5
      *IT_0052;
    const ccomplex_t IT_0054 = pow(m_Z, -2);
    const ccomplex_t IT_0055 = pow(m_Z, 2);
    const ccomplex_t IT_0056 = s_13*IT_0055;
    const ccomplex_t IT_0057 = s_14*s_34;
    const ccomplex_t IT_0058 = (-2)*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = IT_0054*IT_0059;
    const ccomplex_t IT_0061 = e_em*IT_0005;
    const ccomplex_t IT_0062 = e_em*IT_0008;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0061 + IT_0062);
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0002*IT_0064;
    const ccomplex_t IT_0066 = pow(m_C_1, 2);
    const ccomplex_t IT_0067 = cpow(s_12 + 0.5*IT_0016 + 0.5*IT_0066 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*-1)*e_em*(IT_0001*IT_0003 + 
      -IT_0000*IT_0004);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = cpow(s_13 + (-0.5)*IT_0066 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0002*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0071*IT_0074;
    const ccomplex_t IT_0076 = -IT_0069 + (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = pow(s_34, 2);
    const ccomplex_t IT_0078 = m_C_1*IT_0054*IT_0077;
    const ccomplex_t IT_0079 = 2*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0040 + 0.5*IT_0051;
    const ccomplex_t IT_0081 = s_34*IT_0055;
    const ccomplex_t IT_0082 = IT_0054*IT_0081;
    const ccomplex_t IT_0083 = (-0.25)*IT_0082;
    const ccomplex_t IT_0084 = s_34 + IT_0083;
    const ccomplex_t IT_0085 = m_C_1*IT_0084;
    const ccomplex_t IT_0086 = (-4)*IT_0085;
    const ccomplex_t IT_0087 = 2*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0069;
    const ccomplex_t IT_0089 = 0.5*IT_0082;
    const ccomplex_t IT_0090 = s_34 + IT_0089;
    const ccomplex_t IT_0091 = m_C_1*IT_0090;
    const ccomplex_t IT_0092 = 2*IT_0091;
    const ccomplex_t IT_0093 = 2*IT_0092;
    const ccomplex_t IT_0094 = -IT_0040 + -IT_0051 + 0.5*IT_0075;
    const ccomplex_t IT_0095 = s_14*s_34*IT_0054;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = s_13 + IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0097;
    const ccomplex_t IT_0099 = (-2)*IT_0098;
    const ccomplex_t IT_0100 = pow(s_14, 2);
    const ccomplex_t IT_0101 = IT_0054*IT_0100;
    const ccomplex_t IT_0102 = -IT_0066;
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = s_13*IT_0103;
    const ccomplex_t IT_0105 = s_13*s_14;
    const ccomplex_t IT_0106 = s_34*IT_0066;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0105 + IT_0107;
    const ccomplex_t IT_0109 = s_13*IT_0097;
    const ccomplex_t IT_0110 = (-2)*IT_0109;
    const ccomplex_t IT_0111 = s_13*s_34;
    const ccomplex_t IT_0112 = s_13*IT_0054*IT_0077;
    const ccomplex_t IT_0113 = 2*IT_0111;
    const ccomplex_t IT_0114 = (-0.5)*conj(IT_0088);
    const ccomplex_t IT_0115 = (-0.5)*conj(IT_0053);
    const ccomplex_t IT_0116 = pow(m_Z, 4);
    const ccomplex_t IT_0117 = s_13*IT_0116;
    const ccomplex_t IT_0118 = IT_0054*IT_0117;
    const ccomplex_t IT_0119 = IT_0056 + IT_0058 + (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = (-4)*IT_0119;
    const ccomplex_t IT_0121 = (-0.25)*IT_0118;
    const ccomplex_t IT_0122 = IT_0057 + IT_0121;
    const ccomplex_t IT_0123 = (-8)*IT_0122;
    const ccomplex_t IT_0124 = 4*IT_0053*(conj(IT_0053)*(s_13 + (-0.5)*IT_0060
      ) + 0.25*conj(IT_0076)*IT_0079 + 0.25*conj(IT_0080)*IT_0087 + 0.25*conj
      (IT_0088)*IT_0093 + 0.25*conj(IT_0094)*IT_0099) + 2*IT_0094*(0.5*conj
      (IT_0053)*IT_0099 + conj(IT_0094)*IT_0104 + (conj(IT_0080) + -conj(IT_0088
      ))*IT_0108 + 0.5*conj(IT_0076)*IT_0110) + (-2)*IT_0076*((-0.5)*conj
      (IT_0094)*IT_0110 + conj(IT_0080)*IT_0111 + -conj(IT_0076)*IT_0112 +
       IT_0113*IT_0114 + IT_0079*IT_0115) + (-2)*IT_0088*(conj(IT_0094)*IT_0108 
      + (-0.5)*conj(IT_0076)*IT_0113 + IT_0093*IT_0115 + IT_0114*IT_0120 + (-0.5
      )*conj(IT_0080)*IT_0123) + 2*IT_0080*(conj(IT_0094)*IT_0108 + -conj
      (IT_0076)*IT_0111 + -IT_0087*IT_0115 + 0.5*conj(IT_0080)*IT_0120 + 
      -IT_0114*IT_0123);
    return create_ccomplex_return(IT_0124);
}

