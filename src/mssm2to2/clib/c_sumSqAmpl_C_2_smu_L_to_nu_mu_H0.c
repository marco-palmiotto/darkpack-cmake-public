#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_smu_L_to_nu_mu_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_smu_L_to_nu_mu_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
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
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0000
      *IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = cpow(s_13 + (-0.5)*IT_0011 + 0.5*m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0016 = cos(alpha);
    const ccomplex_t IT_0017 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = IT_0005*(IT_0018 + IT_0020);
    const ccomplex_t IT_0022 = 1.4142135623731*e_em*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0027 = pow(m_smu_L, 2);
    const ccomplex_t IT_0028 = cpow(s_23 + (-0.5)*IT_0027 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0005;
    const ccomplex_t IT_0031 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0034 = IT_0003*IT_0033;
    const ccomplex_t IT_0035 = IT_0005*(IT_0032 + IT_0034);
    const ccomplex_t IT_0036 = 1.4142135623731*e_em*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0030*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0041 = cpow(s_23 + (-0.5)*IT_0027 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0040*IT_0041;
    const ccomplex_t IT_0043 = U_d1*V_Wp2;
    const ccomplex_t IT_0044 = IT_0016*IT_0043;
    const ccomplex_t IT_0045 = V_u2*U_Wm1;
    const ccomplex_t IT_0046 = IT_0003*IT_0045;
    const ccomplex_t IT_0047 = IT_0005*(IT_0044 + IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0015*IT_0050;
    const ccomplex_t IT_0052 = IT_0028*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_C_2*IT_0053;
    const ccomplex_t IT_0055 = U_d2*V_Wp2;
    const ccomplex_t IT_0056 = IT_0016*IT_0055;
    const ccomplex_t IT_0057 = V_u2*U_Wm2;
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = IT_0005*(IT_0056 + IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0030*IT_0062;
    const ccomplex_t IT_0064 = IT_0041*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = m_C_2*IT_0065;
    const ccomplex_t IT_0067 = pow(m_mu, 2);
    const ccomplex_t IT_0068 = sin(beta);
    const ccomplex_t IT_0069 = cos(theta_W);
    const ccomplex_t IT_0070 = cpow(IT_0069, -2);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0016*IT_0067 + (-0.5)*m_W*(IT_0001*IT_0016 + -IT_0003*IT_0068)
      *(IT_0005 + -IT_0004*IT_0070));
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0030*IT_0072;
    const ccomplex_t IT_0074 = cpow(s_13 + (-0.5)*IT_0011 + 0.5*m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0014 + (-0.5)*IT_0029 + (-0.5)*IT_0042 +
       (-0.5)*IT_0054 + (-0.5)*IT_0066 + 0.5*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0053 + 0.5*IT_0065;
    const ccomplex_t IT_0079 = s_34*m_C_2;
    const ccomplex_t IT_0080 = pow(m_H0, 2);
    const ccomplex_t IT_0081 = s_13*IT_0080;
    const ccomplex_t IT_0082 = s_14*s_34;
    return create_ccomplex_return(2*IT_0077*(s_13*conj(IT_0077) + conj(IT_0078
      )*IT_0079) + 2*IT_0078*(conj(IT_0077)*IT_0079 + -conj(IT_0078)*(IT_0081 + 
      (-2)*IT_0082)));
}

