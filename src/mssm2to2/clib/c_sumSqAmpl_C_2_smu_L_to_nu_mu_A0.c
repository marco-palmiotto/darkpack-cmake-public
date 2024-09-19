#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_smu_L_to_nu_mu_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_smu_L_to_nu_mu_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = e_em*m_mu*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = pow(m_C_2, 2);
    const ccomplex_t IT_0011 = cpow(s_13 + (-0.5)*IT_0010 + 0.5*m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0002;
    const ccomplex_t IT_0015 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0002*(IT_0016 + IT_0019);
    const ccomplex_t IT_0021 = 1.4142135623731*e_em*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0025 = pow(m_smu_L, 2);
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0025 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0002;
    const ccomplex_t IT_0029 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0032 = IT_0017*IT_0031;
    const ccomplex_t IT_0033 = IT_0002*(IT_0030 + IT_0032);
    const ccomplex_t IT_0034 = 1.4142135623731*e_em*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0038 = cpow(s_23 + (-0.5)*IT_0025 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0037*IT_0038;
    const ccomplex_t IT_0040 = V_u2*U_Wm1;
    const ccomplex_t IT_0041 = IT_0005*IT_0040;
    const ccomplex_t IT_0042 = U_d1*V_Wp2;
    const ccomplex_t IT_0043 = IT_0017*IT_0042;
    const ccomplex_t IT_0044 = IT_0002*(IT_0041 + IT_0043);
    const ccomplex_t IT_0045 = 1.4142135623731*e_em*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = IT_0026*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = m_C_2*IT_0049;
    const ccomplex_t IT_0051 = V_u2*U_Wm2;
    const ccomplex_t IT_0052 = IT_0005*IT_0051;
    const ccomplex_t IT_0053 = U_d2*V_Wp2;
    const ccomplex_t IT_0054 = IT_0017*IT_0053;
    const ccomplex_t IT_0055 = IT_0002*(IT_0052 + IT_0054);
    const ccomplex_t IT_0056 = 1.4142135623731*e_em*IT_0055;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0028*IT_0057;
    const ccomplex_t IT_0059 = IT_0038*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = m_C_2*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0013 + (-0.5)*IT_0027 + (-0.5)*IT_0039 +
       (-0.5)*IT_0050 + (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0049 + 0.5*IT_0060;
    const ccomplex_t IT_0064 = s_34*m_C_2;
    const ccomplex_t IT_0065 = pow(m_A0, 2);
    const ccomplex_t IT_0066 = s_13*IT_0065;
    const ccomplex_t IT_0067 = s_14*s_34;
    return create_ccomplex_return(2*IT_0062*(s_13*conj(IT_0062) + conj(IT_0063
      )*IT_0064) + 2*IT_0063*(conj(IT_0062)*IT_0064 + -conj(IT_0063)*(IT_0066 + 
      (-2)*IT_0067)));
}

