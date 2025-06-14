#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_st_1_to_nu_e_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_st_1_to_nu_e_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_e = param->m_e;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_R = param->m_se_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = m_t*V_tb*conj(V_u1)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0019 = pow(m_se_R, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = m_t*V_tb*conj(V_u2)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0033 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0031*IT_0032*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0021 + 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0000*IT_0002*IT_0004*U_st_00;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = IT_0020*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0000*IT_0002*IT_0004*U_st_00;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0023*IT_0042;
    const ccomplex_t IT_0044 = IT_0033*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0040 + IT_0045;
    const ccomplex_t IT_0047 = m_b*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = (-0.5)*IT_0040 + (-0.5)*IT_0045;
    const ccomplex_t IT_0050 = m_b*s_23;
    const ccomplex_t IT_0051 = s_23*s_24;
    const ccomplex_t IT_0052 = pow(m_st_1, 2);
    const ccomplex_t IT_0053 = s_34*IT_0052;
    const ccomplex_t IT_0054 = (-6)*s_34;
    const ccomplex_t IT_0055 = 6*s_34*(IT_0035*conj(IT_0035) + IT_0048*conj
      (IT_0048)) + 6*((conj(IT_0035) + -conj(IT_0048))*IT_0049 + (IT_0035 + 
      -IT_0048)*conj(IT_0049))*IT_0050 + IT_0049*conj(IT_0049)*(12*IT_0051 + (-6
      )*IT_0053) + (conj(IT_0035)*IT_0048 + IT_0035*conj(IT_0048))*IT_0054;
    return create_ccomplex_return(IT_0055);
}

