#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_L_to_anti_d_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_L_to_anti_d_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ur = param->Gamma_ur;
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
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*(IT_0006 + IT_0009);
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = pow(m_d, 2);
    const ccomplex_t IT_0015 = pow(m_su_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = IT_0001*(IT_0022 + IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0020*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0018 + IT_0031;
    const ccomplex_t IT_0033 = m_C_1*IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = cpow(IT_0004, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0001*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = V_u1*U_Wm1;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = U_d1*V_Wp1;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = IT_0001*(IT_0039 + IT_0041);
    const ccomplex_t IT_0043 = 1.4142135623731*e_em*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = IT_0016*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0001*IT_0034*IT_0035;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = V_u1*U_Wm2;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = U_d2*V_Wp1;
    const ccomplex_t IT_0053 = IT_0007*IT_0052;
    const ccomplex_t IT_0054 = IT_0001*(IT_0051 + IT_0053);
    const ccomplex_t IT_0055 = 1.4142135623731*e_em*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0049*IT_0056;
    const ccomplex_t IT_0058 = IT_0029*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0047 + IT_0059;
    const ccomplex_t IT_0061 = m_C_1*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = pow(m_A0, 2);
    const ccomplex_t IT_0064 = m_d*m_C_1*IT_0063;
    const ccomplex_t IT_0065 = m_d*e_em*IT_0001*IT_0007*IT_0034*IT_0035;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0037*IT_0066;
    const ccomplex_t IT_0068 = pow(m_C_1, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0068 
      + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0018 + -IT_0031;
    const ccomplex_t IT_0073 = conj(IT_0071) + conj(IT_0072);
    const ccomplex_t IT_0074 = IT_0003*IT_0066;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0060 + IT_0076;
    const ccomplex_t IT_0078 = IT_0071 + IT_0072;
    const ccomplex_t IT_0079 = conj(IT_0060) + conj(IT_0076);
    const ccomplex_t IT_0080 = m_d*m_C_1;
    const ccomplex_t IT_0081 = m_u*e_em*mu_h*IT_0001*IT_0034;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0007, -1);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0034*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0082*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0014 + IT_0068 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0091 = IT_0003*IT_0044;
    const ccomplex_t IT_0092 = IT_0016*IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0094 = IT_0020*IT_0056;
    const ccomplex_t IT_0095 = IT_0029*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0089 + -IT_0092 + -IT_0095;
    const ccomplex_t IT_0097 = IT_0012*IT_0037;
    const ccomplex_t IT_0098 = IT_0016*IT_0090*IT_0097;
    const ccomplex_t IT_0099 = IT_0027*IT_0049;
    const ccomplex_t IT_0100 = IT_0029*IT_0093*IT_0099;
    const ccomplex_t IT_0101 = IT_0098 + IT_0100;
    const ccomplex_t IT_0102 = 6*IT_0033;
    const ccomplex_t IT_0103 = 6*conj(IT_0033);
    const ccomplex_t IT_0104 = s_13*IT_0063;
    const ccomplex_t IT_0105 = s_14*s_34;
    const ccomplex_t IT_0106 = (-2)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = (-6)*IT_0071;
    const ccomplex_t IT_0109 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0110 = (-6)*IT_0076;
    const ccomplex_t IT_0111 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0112 = (-6)*s_13;
    const ccomplex_t IT_0113 = m_d*s_14;
    const ccomplex_t IT_0114 = (-6)*IT_0113;
    const ccomplex_t IT_0115 = s_34*m_C_1;
    const ccomplex_t IT_0116 = (-6)*IT_0115;
    const ccomplex_t IT_0117 = (-6)*IT_0080;
    const ccomplex_t IT_0118 = (-6)*IT_0107;
    const ccomplex_t IT_0119 = 0.166666666666667*conj(IT_0096);
    const ccomplex_t IT_0120 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0121 = 0.166666666666667*conj(IT_0101);
    const ccomplex_t IT_0122 = 6*IT_0113*(IT_0062*IT_0073 + conj(IT_0062)
      *IT_0078 + 0.166666666666667*IT_0079*IT_0102 + 0.166666666666667*IT_0077
      *IT_0103 + 0.166666666666667*IT_0096*IT_0111 + IT_0110*IT_0119 + IT_0109
      *IT_0120 + IT_0108*IT_0121);
    const ccomplex_t IT_0123 = 6*IT_0115*(conj(IT_0062)*IT_0077 + IT_0062
      *IT_0079 + 0.166666666666667*IT_0073*IT_0102 + 0.166666666666667*IT_0078
      *IT_0103 + 0.166666666666667*IT_0096*IT_0109 + IT_0108*IT_0119 + IT_0111
      *IT_0120 + IT_0110*IT_0121);
    const ccomplex_t IT_0124 = 0.166666666666667*conj(IT_0062);
    const ccomplex_t IT_0125 = 0.166666666666667*conj(IT_0060);
    const ccomplex_t IT_0126 = 0.166666666666667*conj(IT_0072);
    const ccomplex_t IT_0127 = 6*s_13*(IT_0033*conj(IT_0033) + IT_0062*conj
      (IT_0062)) + 6*IT_0064*(IT_0073*IT_0077 + IT_0078*IT_0079) + (-6)*IT_0080*
      ((conj(IT_0033) + -conj(IT_0096))*IT_0101 + (IT_0033 + -IT_0096)*conj
      (IT_0101) + (-0.166666666666667)*conj(IT_0062)*IT_0102 + (
      -0.166666666666667)*IT_0062*IT_0103) + IT_0107*(IT_0073*IT_0108 + IT_0072
      *IT_0109 + IT_0079*IT_0110 + IT_0060*IT_0111) + conj(IT_0101)*(IT_0062
      *IT_0112 + IT_0072*IT_0114 + IT_0060*IT_0116) + conj(IT_0096)*(IT_0033
      *IT_0112 + IT_0060*IT_0114 + IT_0072*IT_0116 + IT_0062*IT_0117) + (IT_0060
      *conj(IT_0060) + IT_0072*conj(IT_0072))*IT_0118 + IT_0122 + IT_0123 + 6
      *IT_0101*(s_13*conj(IT_0101) + IT_0112*IT_0124 + IT_0116*IT_0125 + IT_0114
      *IT_0126) + 6*IT_0096*(s_13*conj(IT_0096) + 0.166666666666667*conj(IT_0033
      )*IT_0112 + IT_0117*IT_0124 + IT_0114*IT_0125 + IT_0116*IT_0126);
    return create_ccomplex_return(IT_0127);
}

