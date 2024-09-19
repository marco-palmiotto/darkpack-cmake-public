#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_d_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_d_A0(
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
    const ccomplex_t IT_0000 = pow(m_A0, 2);
    const ccomplex_t IT_0001 = m_d*m_C_2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_d, 2);
    const ccomplex_t IT_0003 = pow(m_su_L, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_23 + IT_0002 + IT_0003 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = V_u2*U_Wm1;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = U_d1*V_Wp2;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (IT_0007 + IT_0010)*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cpow(IT_0005, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0012*IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = IT_0004*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0002 + IT_0003 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0024 = V_u2*U_Wm2;
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = U_d2*V_Wp2;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = IT_0012*(IT_0025 + IT_0027);
    const ccomplex_t IT_0029 = 1.4142135623731*e_em*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0012*IT_0016*IT_0017;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = IT_0023*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0022 + IT_0035;
    const ccomplex_t IT_0037 = pow(m_C_2, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + IT_0002 + -IT_0003 + -IT_0037 
      + -reg_prop, -1);
    const ccomplex_t IT_0039 = m_d*e_em*IT_0008*IT_0012*IT_0016*IT_0017;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0012;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = conj(IT_0036) + conj(IT_0045);
    const ccomplex_t IT_0047 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0048 = IT_0005*IT_0047;
    const ccomplex_t IT_0049 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = IT_0012*(IT_0048 + IT_0050);
    const ccomplex_t IT_0052 = 1.4142135623731*e_em*IT_0051;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0042*IT_0053;
    const ccomplex_t IT_0055 = IT_0023*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0012;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0060 = IT_0005*IT_0059;
    const ccomplex_t IT_0061 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0062 = IT_0008*IT_0061;
    const ccomplex_t IT_0063 = IT_0012*(IT_0060 + IT_0062);
    const ccomplex_t IT_0064 = 1.4142135623731*e_em*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0058*IT_0065;
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0056 + -IT_0068;
    const ccomplex_t IT_0070 = IT_0032*IT_0040;
    const ccomplex_t IT_0071 = IT_0038*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = IT_0036 + IT_0045;
    const ccomplex_t IT_0075 = conj(IT_0069) + conj(IT_0072);
    const ccomplex_t IT_0076 = m_u*e_em*mu_h*IT_0012*IT_0016;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = cpow(IT_0008, -1);
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0012*IT_0016*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0077*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0002 + IT_0037 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0086 = IT_0030*IT_0042;
    const ccomplex_t IT_0087 = IT_0023*IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0089 = IT_0015*IT_0058;
    const ccomplex_t IT_0090 = IT_0004*IT_0088*IT_0089;
    const ccomplex_t IT_0091 = IT_0084 + -IT_0087 + -IT_0090;
    const ccomplex_t IT_0092 = m_C_2*IT_0036;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0032*IT_0053;
    const ccomplex_t IT_0095 = IT_0023*IT_0085*IT_0094;
    const ccomplex_t IT_0096 = IT_0019*IT_0065;
    const ccomplex_t IT_0097 = IT_0004*IT_0088*IT_0096;
    const ccomplex_t IT_0098 = IT_0095 + IT_0097;
    const ccomplex_t IT_0099 = IT_0056 + IT_0068;
    const ccomplex_t IT_0100 = m_C_2*IT_0099;
    const ccomplex_t IT_0101 = 6*s_13;
    const ccomplex_t IT_0102 = (-6)*s_13;
    const ccomplex_t IT_0103 = m_d*m_C_2;
    const ccomplex_t IT_0104 = (-6)*IT_0091;
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0100 + IT_0105;
    const ccomplex_t IT_0107 = (-6)*conj(IT_0091);
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = conj(IT_0100) + IT_0108;
    const ccomplex_t IT_0110 = (-6)*IT_0098;
    const ccomplex_t IT_0111 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0112 = s_34*m_C_2;
    const ccomplex_t IT_0113 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0114 = IT_0093 + IT_0113;
    const ccomplex_t IT_0115 = m_d*s_14;
    const ccomplex_t IT_0116 = s_13*IT_0000;
    const ccomplex_t IT_0117 = s_14*s_34;
    const ccomplex_t IT_0118 = (-2)*IT_0117;
    const ccomplex_t IT_0119 = IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = (-6)*IT_0036;
    const ccomplex_t IT_0121 = (-6)*IT_0069;
    const ccomplex_t IT_0122 = 6*IT_0001*(IT_0046*IT_0073 + IT_0074*IT_0075) +
       (IT_0091*conj(IT_0091) + IT_0093*conj(IT_0093) + IT_0098*conj(IT_0098) +
       IT_0100*conj(IT_0100))*IT_0101 + (conj(IT_0093)*IT_0098 + IT_0093*conj
      (IT_0098) + conj(IT_0091)*IT_0100 + IT_0091*conj(IT_0100))*IT_0102 + 6
      *IT_0103*(conj(IT_0091)*IT_0098 + IT_0091*conj(IT_0098) + conj(IT_0093)
      *IT_0106 + IT_0093*IT_0109 + 0.166666666666667*conj(IT_0100)*IT_0110 +
       0.166666666666667*IT_0100*IT_0111) + 6*IT_0112*(IT_0075*IT_0106 + IT_0073
      *IT_0109 + IT_0074*(conj(IT_0093) + 0.166666666666667*IT_0111) + IT_0046
      *IT_0114) + 6*(IT_0046*IT_0106 + IT_0074*IT_0109 + IT_0073*(conj(IT_0093) 
      + 0.166666666666667*IT_0111) + IT_0075*IT_0114)*IT_0115 + (-6)*IT_0119*
      (IT_0046*(IT_0045 + (-0.166666666666667)*IT_0120) + IT_0075*(IT_0072 + (
      -0.166666666666667)*IT_0121));
    return create_ccomplex_return(IT_0122);
}

