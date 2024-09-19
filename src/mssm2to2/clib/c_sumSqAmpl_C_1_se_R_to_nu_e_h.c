#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_se_R_to_nu_e_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_se_R_to_nu_e_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_el = param->Gamma_el;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = V_u1*U_Wm1;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = U_d1*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*(IT_0009 + -IT_0012);
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = pow(m_se_R, 2);
    const ccomplex_t IT_0019 = cpow(s_23 + (-0.5)*IT_0018 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = V_u1*U_Wm2;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = U_d2*V_Wp1;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = IT_0004*(IT_0025 + -IT_0027);
    const ccomplex_t IT_0029 = 1.4142135623731*e_em*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_23 + (-0.5)*IT_0018 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0021 + 0.5*IT_0035;
    const ccomplex_t IT_0037 = pow(m_h, 2);
    const ccomplex_t IT_0038 = s_13*IT_0037;
    const ccomplex_t IT_0039 = s_14*s_34;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0000
      *IT_0002*IT_0004*IT_0007;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = pow(m_C_1, 2);
    const ccomplex_t IT_0045 = cpow(s_13 + (-0.5)*IT_0044 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cos(theta_W);
    const ccomplex_t IT_0049 = cpow(IT_0048, -2);
    const ccomplex_t IT_0050 = sin(beta);
    const ccomplex_t IT_0051 = pow(m_e, 2);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0049*
      (IT_0001*IT_0010 + IT_0007*IT_0050) + -IT_0000*IT_0002*IT_0004*IT_0010
      *IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0006*IT_0053;
    const ccomplex_t IT_0055 = cpow(s_13 + (-0.5)*IT_0044 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0059 = IT_0007*IT_0058;
    const ccomplex_t IT_0060 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0061 = IT_0010*IT_0060;
    const ccomplex_t IT_0062 = IT_0004*(IT_0059 + -IT_0061);
    const ccomplex_t IT_0063 = 1.4142135623731*e_em*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0006*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0068 = IT_0019*IT_0066*IT_0067;
    const ccomplex_t IT_0069 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0070 = IT_0007*IT_0069;
    const ccomplex_t IT_0071 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0072 = IT_0010*IT_0071;
    const ccomplex_t IT_0073 = IT_0004*(IT_0070 + -IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0023*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0079 = IT_0033*IT_0077*IT_0078;
    const ccomplex_t IT_0080 = m_C_1*IT_0021;
    const ccomplex_t IT_0081 = m_C_1*IT_0035;
    const ccomplex_t IT_0082 = 0.5*IT_0047 + 0.5*IT_0057 + (-0.5)*IT_0068 + (
      -0.5)*IT_0079 + (-0.5)*IT_0080 + (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = s_34*m_C_1;
    const ccomplex_t IT_0084 = 2*IT_0083;
    return create_ccomplex_return(2*IT_0082*(s_13*conj(IT_0082) + 0.5*conj
      (IT_0036)*IT_0084) + (-2)*IT_0036*(conj(IT_0036)*(IT_0038 + (-2)*IT_0039) 
      + (-0.5)*conj(IT_0082)*IT_0084));
}

