#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_sb_2_to_d_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_sb_2_to_d_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_R = param->m_su_R;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0000*IT_0002*IT_0004*U_sb_01;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_d, 2);
    const ccomplex_t IT_0011 = pow(m_su_R, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0000*IT_0002*IT_0004*U_sb_01;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014 + IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0004*IT_0026;
    const ccomplex_t IT_0028 = m_b*conj(U_d2)*V_tb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0030 + (-2)*IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = IT_0020*IT_0024*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = m_t*s_13;
    const ccomplex_t IT_0039 = s_13*s_14;
    const ccomplex_t IT_0040 = s_34*IT_0011;
    const ccomplex_t IT_0041 = m_d*IT_0023;
    const ccomplex_t IT_0042 = m_d*m_t;
    const ccomplex_t IT_0043 = 18*s_34;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0045 = m_b*conj(U_d1)*V_tb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0046 = IT_0027*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0047 + (-2)*IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = IT_0012*IT_0044*IT_0052;
    const ccomplex_t IT_0054 = (-18)*s_34;
    const ccomplex_t IT_0055 = m_d*s_14;
    const ccomplex_t IT_0056 = (-18)*IT_0023;
    const ccomplex_t IT_0057 = (-18)*conj(IT_0023);
    const ccomplex_t IT_0058 = IT_0023*(18*conj(IT_0037)*IT_0038 + conj
      (IT_0023)*(36*IT_0039 + (-18)*IT_0040)) + IT_0037*(18*conj(IT_0023)
      *IT_0038 + (-18)*conj(IT_0041)*IT_0042 + conj(IT_0037)*IT_0043 + conj
      (IT_0053)*IT_0054) + (IT_0038*conj(IT_0053) + conj(IT_0041)*IT_0055)
      *IT_0056 + IT_0053*(18*conj(IT_0041)*IT_0042 + IT_0043*conj(IT_0053) +
       conj(IT_0037)*IT_0054 + IT_0038*IT_0057) + IT_0041*(conj(IT_0041)*IT_0043
       + IT_0042*((-18)*conj(IT_0037) + 18*conj(IT_0053)) + IT_0055*IT_0057);
    return create_ccomplex_return(IT_0058);
}

