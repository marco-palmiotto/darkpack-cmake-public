#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_st_1_to_nu_e_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_st_1_to_nu_e_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
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
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = m_t*V_tb*conj(V_u1)*e_em*IT_0008*U_st_10;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0004 + (-0.5)*IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0016 = pow(m_se_L, 2);
    const ccomplex_t IT_0017 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0020 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = m_t*V_tb*conj(V_u2)*e_em*IT_0008*U_st_10;
    const ccomplex_t IT_0023 = IT_0007*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + (-0.5)*IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0029 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0018 + 0.5*IT_0030;
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0001*IT_0008*IT_0033*U_st_00;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0035;
    const ccomplex_t IT_0037 = IT_0017*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0001*IT_0008*IT_0033*U_st_00;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = IT_0029*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0038 + IT_0043;
    const ccomplex_t IT_0045 = m_b*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0038 + (-0.5)*IT_0043;
    const ccomplex_t IT_0048 = m_b*s_23;
    const ccomplex_t IT_0049 = s_23*s_24;
    const ccomplex_t IT_0050 = pow(m_st_1, 2);
    const ccomplex_t IT_0051 = s_34*IT_0050;
    const ccomplex_t IT_0052 = (-6)*s_34;
    const ccomplex_t IT_0053 = 6*s_34*(IT_0031*conj(IT_0031) + IT_0046*conj
      (IT_0046)) + 6*((conj(IT_0031) + -conj(IT_0046))*IT_0047 + (IT_0031 + 
      -IT_0046)*conj(IT_0047))*IT_0048 + IT_0047*conj(IT_0047)*(12*IT_0049 + (-6
      )*IT_0051) + (conj(IT_0031)*IT_0046 + IT_0031*conj(IT_0046))*IT_0052;
    return create_ccomplex_return(IT_0053);
}

