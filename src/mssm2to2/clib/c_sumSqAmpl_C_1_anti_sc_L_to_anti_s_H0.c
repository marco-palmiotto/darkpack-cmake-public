#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*(IT_0006 + IT_0009);
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_s, 2);
    const ccomplex_t IT_0016 = pow(m_sc_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*(IT_0023 + IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0019 + IT_0033;
    const ccomplex_t IT_0035 = m_C_1*IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, -1);
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0001*IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = U_d1*V_Wp1;
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = V_u1*U_Wm1;
    const ccomplex_t IT_0044 = IT_0007*IT_0043;
    const ccomplex_t IT_0045 = IT_0001*(IT_0042 + IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0040*IT_0048;
    const ccomplex_t IT_0050 = IT_0031*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0001*IT_0036*IT_0038;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = U_d2*V_Wp1;
    const ccomplex_t IT_0055 = IT_0004*IT_0054;
    const ccomplex_t IT_0056 = V_u1*U_Wm2;
    const ccomplex_t IT_0057 = IT_0007*IT_0056;
    const ccomplex_t IT_0058 = IT_0001*(IT_0055 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0051 + IT_0064;
    const ccomplex_t IT_0066 = m_C_1*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = pow(m_c, 2);
    const ccomplex_t IT_0069 = sin(beta);
    const ccomplex_t IT_0070 = cpow(IT_0069, -1);
    const ccomplex_t IT_0071 = cos(theta_W);
    const ccomplex_t IT_0072 = cpow(IT_0071, -2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0007
      *IT_0036*IT_0068*IT_0070 + 0.5*m_W*(IT_0004*IT_0037 + -IT_0007*IT_0069)*
      (IT_0001 + -1./3*IT_0000*IT_0072));
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0040*IT_0074;
    const ccomplex_t IT_0076 = pow(m_C_1, 2);
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0015 + IT_0076 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0029*IT_0040;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0082 = IT_0031*IT_0080*IT_0081;
    const ccomplex_t IT_0083 = IT_0013*IT_0053;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0085 = IT_0017*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_s*e_em*IT_0001*IT_0004
      *IT_0036*IT_0038;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0040*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + IT_0015 + -IT_0016 + -IT_0076 
      + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = m_s*IT_0091;
    const ccomplex_t IT_0093 = -IT_0079 + IT_0082 + IT_0085 + 2*IT_0092;
    const ccomplex_t IT_0094 = m_s*m_C_1;
    const ccomplex_t IT_0095 = IT_0021*IT_0074;
    const ccomplex_t IT_0096 = IT_0077*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0003*IT_0061;
    const ccomplex_t IT_0099 = IT_0017*IT_0084*IT_0098;
    const ccomplex_t IT_0100 = IT_0021*IT_0087;
    const ccomplex_t IT_0101 = IT_0089*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = m_s*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0001
      *IT_0004*IT_0036*IT_0070;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0036*IT_0070;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0015 + IT_0076 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0021*IT_0048;
    const ccomplex_t IT_0113 = IT_0031*IT_0081*IT_0112;
    const ccomplex_t IT_0114 = IT_0097 + -IT_0099 + (-2)*IT_0103 + IT_0111 + 
      -IT_0113;
    const ccomplex_t IT_0115 = pow(m_H0, 2);
    const ccomplex_t IT_0116 = m_s*m_C_1*IT_0115;
    const ccomplex_t IT_0117 = conj(IT_0065) + conj(IT_0102);
    const ccomplex_t IT_0118 = -IT_0019 + -IT_0033;
    const ccomplex_t IT_0119 = -IT_0091;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = IT_0065 + IT_0102;
    const ccomplex_t IT_0122 = conj(IT_0118) + conj(IT_0119);
    const ccomplex_t IT_0123 = m_s*s_14;
    const ccomplex_t IT_0124 = (-6)*IT_0123;
    const ccomplex_t IT_0125 = conj(IT_0093)*IT_0124;
    const ccomplex_t IT_0126 = s_34*m_C_1;
    const ccomplex_t IT_0127 = (-6)*IT_0126;
    const ccomplex_t IT_0128 = conj(IT_0114)*IT_0127;
    const ccomplex_t IT_0129 = IT_0125 + IT_0128;
    const ccomplex_t IT_0130 = (-6)*IT_0102;
    const ccomplex_t IT_0131 = (-6)*conj(IT_0102);
    const ccomplex_t IT_0132 = s_13*IT_0115;
    const ccomplex_t IT_0133 = s_14*s_34;
    const ccomplex_t IT_0134 = (-2)*IT_0133;
    const ccomplex_t IT_0135 = IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = (-6)*s_13;
    const ccomplex_t IT_0137 = (-6)*IT_0094;
    const ccomplex_t IT_0138 = IT_0093*IT_0124;
    const ccomplex_t IT_0139 = (-6)*IT_0135;
    const ccomplex_t IT_0140 = 6*s_13*(IT_0035*conj(IT_0035) + IT_0067*conj
      (IT_0067) + IT_0093*conj(IT_0093)) + 6*IT_0094*(conj(IT_0035)*(IT_0067 + 
      -IT_0093) + IT_0035*(conj(IT_0067) + -conj(IT_0093)) + IT_0093*conj
      (IT_0114)) + 6*IT_0116*(IT_0117*IT_0120 + IT_0121*IT_0122) + (IT_0118 +
       IT_0119)*IT_0129 + 6*IT_0126*(IT_0067*IT_0117 + conj(IT_0035)*IT_0120 +
       conj(IT_0067)*IT_0121 + IT_0035*IT_0122 + 0.166666666666667*conj(IT_0093)
      *IT_0130 + 0.166666666666667*IT_0093*IT_0131) + 6*IT_0123*(IT_0035*IT_0117
       + conj(IT_0067)*IT_0120 + conj(IT_0035)*IT_0121 + IT_0067*IT_0122 +
       0.166666666666667*conj(IT_0114)*IT_0130 + 0.166666666666667*IT_0114
      *IT_0131) + (IT_0117*IT_0130 + IT_0065*IT_0131)*IT_0135 + conj(IT_0093)*
      (IT_0065*IT_0127 + IT_0067*IT_0136) + IT_0093*(conj(IT_0065)*IT_0127 +
       conj(IT_0067)*IT_0136) + conj(IT_0114)*(IT_0065*IT_0124 + IT_0035*IT_0136
       + IT_0067*IT_0137) + 6*IT_0114*(conj(IT_0093)*IT_0094 + s_13*conj(IT_0114
      ) + 0.166666666666667*conj(IT_0065)*IT_0124 + 0.166666666666667*IT_0122
      *IT_0127 + 0.166666666666667*conj(IT_0035)*IT_0136 + 0.166666666666667
      *conj(IT_0067)*IT_0137) + IT_0122*IT_0138 + (IT_0065*conj(IT_0065) +
       IT_0120*IT_0122)*IT_0139;
    return create_ccomplex_return(IT_0140);
}

