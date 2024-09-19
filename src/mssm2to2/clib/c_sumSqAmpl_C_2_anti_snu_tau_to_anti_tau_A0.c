#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = m_C_2*m_tau;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0002 = pow(m_tau, 2);
    const ccomplex_t IT_0003 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0004 = cpow(s_23 + (-0.5)*IT_0002 + (-0.5)*IT_0003 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0006;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = V_u2*U_Wm1;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = U_d1*V_Wp2;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*(IT_0011 + IT_0014);
    const ccomplex_t IT_0016 = 1.4142135623731*e_em*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0008*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0004*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0021 = cpow(s_23 + (-0.5)*IT_0002 + (-0.5)*IT_0003 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0006;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = V_u2*U_Wm2;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = U_d2*V_Wp2;
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = IT_0006*(IT_0025 + IT_0027);
    const ccomplex_t IT_0029 = 1.4142135623731*e_em*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = IT_0020*IT_0021*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0019 + 0.5*IT_0032;
    const ccomplex_t IT_0034 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0037 = IT_0012*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*(IT_0035 + IT_0037);
    const ccomplex_t IT_0039 = 1.4142135623731*e_em*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0008*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0047 = IT_0012*IT_0046;
    const ccomplex_t IT_0048 = IT_0006*(IT_0045 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0023*IT_0050;
    const ccomplex_t IT_0052 = IT_0021*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0043 + IT_0053;
    const ccomplex_t IT_0055 = m_C_2*IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = pow(m_W, -1);
    const ccomplex_t IT_0058 = cpow(IT_0009, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0006*IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0040*IT_0060;
    const ccomplex_t IT_0062 = IT_0001*IT_0004*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0006*IT_0057*IT_0058;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0050*IT_0064;
    const ccomplex_t IT_0066 = IT_0020*IT_0021*IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0062 + (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0017*IT_0060;
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0030*IT_0064;
    const ccomplex_t IT_0072 = IT_0021*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0070 + IT_0073;
    const ccomplex_t IT_0075 = m_C_2*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = 2*IT_0056;
    const ccomplex_t IT_0078 = (-2)*IT_0033;
    const ccomplex_t IT_0079 = IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = 2*conj(IT_0056);
    const ccomplex_t IT_0081 = (-2)*conj(IT_0033);
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = s_14*m_tau;
    const ccomplex_t IT_0084 = e_em*m_tau*IT_0006*IT_0012*IT_0057*IT_0058;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0023*IT_0085;
    const ccomplex_t IT_0087 = pow(m_C_2, 2);
    const ccomplex_t IT_0088 = cpow(s_12 + (-0.5)*IT_0002 + 0.5*IT_0003 + 0.5
      *IT_0087 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0070 + (-0.5)*IT_0073;
    const ccomplex_t IT_0093 = conj(IT_0091) + conj(IT_0092);
    const ccomplex_t IT_0094 = IT_0091 + IT_0092;
    const ccomplex_t IT_0095 = IT_0064*IT_0085;
    const ccomplex_t IT_0096 = IT_0088*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0043 + 0.5*IT_0053;
    const ccomplex_t IT_0100 = conj(IT_0098) + conj(IT_0099);
    const ccomplex_t IT_0101 = 2*IT_0076;
    const ccomplex_t IT_0102 = (-2)*IT_0067;
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = IT_0098 + IT_0099;
    const ccomplex_t IT_0105 = 2*conj(IT_0076);
    const ccomplex_t IT_0106 = (-2)*conj(IT_0067);
    const ccomplex_t IT_0107 = IT_0105 + IT_0106;
    const ccomplex_t IT_0108 = s_34*m_C_2;
    const ccomplex_t IT_0109 = 2*s_13;
    const ccomplex_t IT_0110 = pow(m_A0, 2);
    const ccomplex_t IT_0111 = m_C_2*m_tau*IT_0110;
    const ccomplex_t IT_0112 = 2*IT_0098;
    const ccomplex_t IT_0113 = 2*conj(IT_0098);
    const ccomplex_t IT_0114 = (-2)*s_13;
    const ccomplex_t IT_0115 = s_13*IT_0110;
    const ccomplex_t IT_0116 = s_14*s_34;
    const ccomplex_t IT_0117 = (-2)*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = (-2)*IT_0098;
    const ccomplex_t IT_0120 = (-2)*IT_0091;
    const ccomplex_t IT_0121 = 2*IT_0000*((conj(IT_0033) + -conj(IT_0056))
      *IT_0067 + (IT_0033 + -IT_0056)*conj(IT_0067) + 0.5*conj(IT_0076)*IT_0079 
      + 0.5*IT_0076*IT_0082) + IT_0083*(IT_0079*IT_0093 + IT_0082*IT_0094 +
       IT_0100*IT_0103 + IT_0104*IT_0107) + (IT_0079*IT_0100 + IT_0093*IT_0103 +
       IT_0082*IT_0104 + IT_0094*IT_0107)*IT_0108 + (IT_0033*conj(IT_0033) +
       IT_0056*conj(IT_0056) + IT_0067*conj(IT_0067) + IT_0076*conj(IT_0076))
      *IT_0109 + 2*IT_0111*(IT_0093*(IT_0099 + 0.5*IT_0112) + IT_0094*(conj
      (IT_0099) + 0.5*IT_0113)) + (conj(IT_0033)*IT_0056 + IT_0033*conj(IT_0056)
       + conj(IT_0067)*IT_0076 + IT_0067*conj(IT_0076))*IT_0114 + (-2)*IT_0118*
      (IT_0100*(IT_0099 + (-0.5)*IT_0119) + IT_0093*(IT_0092 + (-0.5)*IT_0120));
    return create_ccomplex_return(IT_0121);
}

