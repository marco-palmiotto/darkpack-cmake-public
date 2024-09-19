#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_sb_2_to_c_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_sb_2_to_c_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_c, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0004*IT_0006*IT_0008*conj(U_sb_01);
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0004*IT_0008*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = IT_0000*IT_0003*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0004*IT_0008*IT_0012;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0004*IT_0006*IT_0008*conj(U_sb_01);
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = -IT_0016 + -IT_0024;
    const ccomplex_t IT_0026 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = IT_0008*IT_0012;
    const ccomplex_t IT_0029 = m_b*U_d2*V_tb*e_em*IT_0004*conj(U_sb_11);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0027 + (-0.5)*IT_0031);
    const ccomplex_t IT_0033 = IT_0018*IT_0032;
    const ccomplex_t IT_0034 = IT_0023*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = V_tb*e_em*U_Wm1*conj(U_sb_01);
    const ccomplex_t IT_0037 = IT_0008*IT_0036;
    const ccomplex_t IT_0038 = m_b*U_d1*V_tb*e_em*IT_0004*conj(U_sb_11);
    const ccomplex_t IT_0039 = IT_0028*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = IT_0014*IT_0041;
    const ccomplex_t IT_0043 = IT_0003*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = -IT_0035 + -IT_0044;
    const ccomplex_t IT_0046 = m_t*s_13;
    const ccomplex_t IT_0047 = m_c*m_t;
    const ccomplex_t IT_0048 = IT_0035 + IT_0044;
    const ccomplex_t IT_0049 = m_c*IT_0048;
    const ccomplex_t IT_0050 = 18*s_34;
    const ccomplex_t IT_0051 = m_c*s_14;
    const ccomplex_t IT_0052 = (-18)*conj(IT_0025);
    const ccomplex_t IT_0053 = s_13*s_14;
    const ccomplex_t IT_0054 = s_34*IT_0002;
    const ccomplex_t IT_0055 = IT_0025*((-18)*conj(IT_0045)*IT_0046 + (-18)
      *IT_0047*conj(IT_0049) + conj(IT_0025)*IT_0050) + IT_0049*(conj(IT_0049)
      *IT_0050 + 18*conj(IT_0045)*IT_0051 + IT_0047*IT_0052) + IT_0045*(18*conj
      (IT_0049)*IT_0051 + IT_0046*IT_0052 + conj(IT_0045)*(36*IT_0053 + (-18)
      *IT_0054));
    return create_ccomplex_return(IT_0055);
}

