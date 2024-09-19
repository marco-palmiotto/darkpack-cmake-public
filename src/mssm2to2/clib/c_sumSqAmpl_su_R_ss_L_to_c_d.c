#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_ss_L_to_c_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_ss_L_to_c_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_c, 2);
    const ccomplex_t IT_0011 = pow(m_ss_L, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0010 + IT_0011 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0010 + IT_0011 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014 + IT_0022;
    const ccomplex_t IT_0024 = m_d*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0029 = IT_0012*IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0004;
    const ccomplex_t IT_0031 = IT_0018*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0033 = IT_0020*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = IT_0029 + IT_0033;
    const ccomplex_t IT_0035 = m_c*m_d;
    const ccomplex_t IT_0036 = -IT_0014 + -IT_0022;
    const ccomplex_t IT_0037 = m_d*s_13;
    const ccomplex_t IT_0038 = 18*s_34;
    const ccomplex_t IT_0039 = m_c*s_14;
    const ccomplex_t IT_0040 = (-18)*conj(IT_0025);
    const ccomplex_t IT_0041 = s_13*s_14;
    const ccomplex_t IT_0042 = pow(m_su_R, 2);
    const ccomplex_t IT_0043 = s_34*IT_0042;
    const ccomplex_t IT_0044 = IT_0025*((-18)*conj(IT_0034)*IT_0035 + (-18)
      *conj(IT_0036)*IT_0037 + conj(IT_0025)*IT_0038) + IT_0034*(conj(IT_0034)
      *IT_0038 + 18*conj(IT_0036)*IT_0039 + IT_0035*IT_0040) + IT_0036*(18*conj
      (IT_0034)*IT_0039 + IT_0037*IT_0040 + conj(IT_0036)*(36*IT_0041 + (-18)
      *IT_0043));
    return create_ccomplex_return(IT_0044);
}

