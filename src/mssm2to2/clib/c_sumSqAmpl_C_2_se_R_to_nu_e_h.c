#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_se_R_to_nu_e_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_se_R_to_nu_e_h(
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (IT_0003 + -IT_0006)*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_se_R, 2);
    const ccomplex_t IT_0020 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0000*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0023 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*(IT_0024 + -IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0022*IT_0033*IT_0034;
    const ccomplex_t IT_0036 = V_u2*U_Wm1;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = U_d1*V_Wp2;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = IT_0008*(IT_0037 + -IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = IT_0020*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_C_2*IT_0046;
    const ccomplex_t IT_0048 = V_u2*U_Wm2;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = U_d2*V_Wp2;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = IT_0008*(IT_0049 + -IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0032*IT_0055;
    const ccomplex_t IT_0057 = IT_0034*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = m_C_2*IT_0058;
    const ccomplex_t IT_0060 = cos(theta_W);
    const ccomplex_t IT_0061 = cpow(IT_0060, -2);
    const ccomplex_t IT_0062 = sin(beta);
    const ccomplex_t IT_0063 = pow(m_e, 2);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*e_em*(m_W*IT_0007*IT_0061*
      (IT_0004*IT_0014 + IT_0001*IT_0062) + -IT_0004*IT_0008*IT_0013*IT_0015
      *IT_0063);
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0032*IT_0065;
    const ccomplex_t IT_0067 = pow(m_C_2, 2);
    const ccomplex_t IT_0068 = cpow(s_13 + (-0.5)*IT_0067 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0001
      *IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0008;
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = cpow(s_13 + (-0.5)*IT_0067 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0021 + (-0.5)*IT_0035 + (-0.5)
      *IT_0047 + (-0.5)*IT_0059 + 0.5*IT_0070 + 0.5*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0046 + 0.5*IT_0058;
    const ccomplex_t IT_0080 = s_34*m_C_2;
    const ccomplex_t IT_0081 = pow(m_h, 2);
    const ccomplex_t IT_0082 = s_13*IT_0081;
    const ccomplex_t IT_0083 = s_14*s_34;
    return create_ccomplex_return(2*IT_0078*(s_13*conj(IT_0078) + conj(IT_0079
      )*IT_0080) + 2*IT_0079*(conj(IT_0078)*IT_0080 + -conj(IT_0079)*(IT_0082 + 
      (-2)*IT_0083)));
}

