#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_nu_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_nu_e_nu_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = e_em*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0003, -1);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = e_em*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + IT_0008);
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = IT_0007*IT_0011;
    const ccomplex_t IT_0014 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0012 + -IT_0013 + (-2)
      *IT_0015);
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = pow(m_C_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*U_Wm2*IT_0006;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = cos(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0006*IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0006*IT_0029*IT_0031;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = cpow(s_23 + (-0.5)*IT_0019 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0028 + (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = IT_0007*IT_0041;
    const ccomplex_t IT_0044 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + -IT_0043 + (-2)
      *IT_0045);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0010*IT_0047;
    const ccomplex_t IT_0049 = IT_0020*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = s_34*IT_0019;
    const ccomplex_t IT_0053 = 4*IT_0052;
    const ccomplex_t IT_0054 = s_13*s_24;
    const ccomplex_t IT_0055 = IT_0000*(IT_0022*(16*conj(IT_0022) + 8*conj
      (IT_0040)) + 8*IT_0040*(conj(IT_0022) + 0.5*conj(IT_0040))) + 8*conj
      (IT_0051)*(IT_0022*IT_0052 + 0.125*IT_0040*IT_0053) + 8*IT_0051*(conj
      (IT_0022)*IT_0052 + 0.125*conj(IT_0040)*IT_0053 + 2*conj(IT_0051)*IT_0054);
    return create_ccomplex_return(IT_0055);
}

