#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_se_L_to_nu_e_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_se_L_to_nu_e_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = V_u1*U_Wm1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = U_d1*V_Wp1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*(IT_0005 + IT_0008);
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = pow(m_se_L, 2);
    const ccomplex_t IT_0014 = cpow(s_23 + (-0.5)*IT_0013 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0018 = V_u1*U_Wm2;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = U_d2*V_Wp1;
    const ccomplex_t IT_0021 = IT_0006*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*(IT_0019 + IT_0021);
    const ccomplex_t IT_0023 = 1.4142135623731*e_em*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0017*IT_0024;
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0013 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0016 + 0.5*IT_0028;
    const ccomplex_t IT_0030 = pow(m_A0, 2);
    const ccomplex_t IT_0031 = s_13*IT_0030;
    const ccomplex_t IT_0032 = s_14*s_34;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_e*e_em*mu_h*IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0003, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0001*IT_0033*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = pow(m_C_1, 2);
    const ccomplex_t IT_0041 = cpow(s_13 + (-0.5)*IT_0040 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0045 = IT_0003*IT_0044;
    const ccomplex_t IT_0046 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0047 = IT_0006*IT_0046;
    const ccomplex_t IT_0048 = IT_0001*(IT_0045 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0002*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0053 = IT_0014*IT_0051*IT_0052;
    const ccomplex_t IT_0054 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0055 = IT_0003*IT_0054;
    const ccomplex_t IT_0056 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0057 = IT_0006*IT_0056;
    const ccomplex_t IT_0058 = IT_0001*(IT_0055 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0063 = IT_0026*IT_0061*IT_0062;
    const ccomplex_t IT_0064 = m_C_1*IT_0016;
    const ccomplex_t IT_0065 = m_C_1*IT_0028;
    const ccomplex_t IT_0066 = 0.5*IT_0043 + (-0.5)*IT_0053 + (-0.5)*IT_0063 +
       (-0.5)*IT_0064 + (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = s_34*m_C_1;
    const ccomplex_t IT_0068 = 2*IT_0067;
    return create_ccomplex_return(2*IT_0066*(s_13*conj(IT_0066) + 0.5*conj
      (IT_0029)*IT_0068) + (-2)*IT_0029*(conj(IT_0029)*(IT_0031 + (-2)*IT_0032) 
      + (-0.5)*conj(IT_0066)*IT_0068));
}

