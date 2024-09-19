#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_sc_R_to_anti_s_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_sc_R_to_anti_s_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*(IT_0009 + IT_0011);
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = pow(m_s, 2);
    const ccomplex_t IT_0017 = pow(m_sc_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = IT_0004*(IT_0024 + IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0020 + -IT_0033;
    const ccomplex_t IT_0035 = cpow(IT_0007, -1);
    const ccomplex_t IT_0036 = m_s*e_em*IT_0000*IT_0001*IT_0004*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = pow(m_C_2, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0039 
      + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = pow(m_A0, 2);
    const ccomplex_t IT_0044 = m_s*m_C_2*IT_0043;
    const ccomplex_t IT_0045 = s_13*IT_0043;
    const ccomplex_t IT_0046 = s_14*s_34;
    const ccomplex_t IT_0047 = (-2)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045 + IT_0047;
    const ccomplex_t IT_0049 = V_u2*U_Wm1;
    const ccomplex_t IT_0050 = IT_0007*IT_0049;
    const ccomplex_t IT_0051 = U_d1*V_Wp2;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = IT_0004*(IT_0050 + IT_0052);
    const ccomplex_t IT_0054 = 1.4142135623731*e_em*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0022*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0058 = IT_0031*IT_0056*IT_0057;
    const ccomplex_t IT_0059 = m_c*e_em*mu_h*IT_0000*IT_0004;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0004;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0016 + IT_0039 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = V_u2*U_Wm2;
    const ccomplex_t IT_0068 = IT_0007*IT_0067;
    const ccomplex_t IT_0069 = U_d2*V_Wp2;
    const ccomplex_t IT_0070 = IT_0001*IT_0069;
    const ccomplex_t IT_0071 = IT_0004*(IT_0068 + IT_0070);
    const ccomplex_t IT_0072 = 1.4142135623731*e_em*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0006*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0076 = IT_0018*IT_0074*IT_0075;
    const ccomplex_t IT_0077 = -IT_0058 + IT_0066 + -IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0000*IT_0004*IT_0035;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0060*IT_0079;
    const ccomplex_t IT_0081 = IT_0064*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = IT_0020 + IT_0033;
    const ccomplex_t IT_0085 = m_C_2*IT_0084;
    const ccomplex_t IT_0086 = 6*s_13;
    const ccomplex_t IT_0087 = (-6)*s_13;
    const ccomplex_t IT_0088 = m_s*s_14;
    const ccomplex_t IT_0089 = (-6)*IT_0077;
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = IT_0085 + IT_0090;
    const ccomplex_t IT_0092 = (-6)*conj(IT_0077);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = conj(IT_0085) + IT_0093;
    const ccomplex_t IT_0095 = (-6)*IT_0083;
    const ccomplex_t IT_0096 = (-6)*conj(IT_0083);
    const ccomplex_t IT_0097 = s_34*m_C_2;
    const ccomplex_t IT_0098 = m_s*m_C_2;
    const ccomplex_t IT_0099 = 6*conj(IT_0034)*(IT_0042*IT_0044 + -IT_0034
      *IT_0048) + 6*conj(IT_0042)*(IT_0034*IT_0044 + -IT_0042*IT_0048) + 
      (IT_0077*conj(IT_0077) + IT_0083*conj(IT_0083) + IT_0085*conj(IT_0085))
      *IT_0086 + (conj(IT_0077)*IT_0085 + IT_0077*conj(IT_0085))*IT_0087 + 6
      *IT_0088*(conj(IT_0042)*IT_0091 + IT_0042*IT_0094 + 0.166666666666667*conj
      (IT_0034)*IT_0095 + 0.166666666666667*IT_0034*IT_0096) + 6*(conj(IT_0034)
      *IT_0091 + IT_0034*IT_0094 + 0.166666666666667*conj(IT_0042)*IT_0095 +
       0.166666666666667*IT_0042*IT_0096)*IT_0097 + 6*(conj(IT_0077)*IT_0083 +
       IT_0077*conj(IT_0083) + 0.166666666666667*conj(IT_0085)*IT_0095 +
       0.166666666666667*IT_0085*IT_0096)*IT_0098;
    return create_ccomplex_return(IT_0099);
}

