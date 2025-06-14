#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_ss_L_to_h_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_ss_L_to_h_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0005*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = pow(m_sc_R, 2);
    const ccomplex_t IT_0012 = pow(m_ss_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0002, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0005*IT_0008*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0022 
      + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0018*IT_0025;
    const ccomplex_t IT_0027 = IT_0016 + IT_0026;
    const ccomplex_t IT_0028 = pow(m_h, 2);
    const ccomplex_t IT_0029 = pow(m_W, -2);
    const ccomplex_t IT_0030 = pow(s_34, 2);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0005;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0003
      *IT_0005*IT_0008*IT_0019;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0011 + IT_0028 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = s_24*s_34*IT_0029;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = s_23 + IT_0041;
    const ccomplex_t IT_0043 = pow(s_24, 2);
    const ccomplex_t IT_0044 = IT_0029*IT_0043;
    return create_ccomplex_return((-3)*IT_0027*(conj(IT_0027)*(IT_0028 + 
      -IT_0031) + conj(IT_0039)*IT_0042) + (-3)*IT_0039*(conj(IT_0027)*IT_0042 +
       conj(IT_0039)*(IT_0012 + -IT_0044)));
}

