#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = cpow(IT_0006, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*(IT_0002*IT_0003*IT_0005*IT_0006 + (-2)*IT_0004*(m_W*IT_0006 + (-0.5
      )*IT_0002*IT_0007*IT_0008))*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = pow(m_Hp, 2);
    const ccomplex_t IT_0014 = pow(m_sc_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = pow(m_ss_L, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + -IT_0014 + -IT_0019 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0012*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0000*IT_0022;
    const ccomplex_t IT_0024 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0013 + IT_0019 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = -IT_0018 + IT_0023 + -IT_0028;
    const ccomplex_t IT_0030 = 2*IT_0000;
    const ccomplex_t IT_0031 = IT_0022*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (-3)*s_14;
    const ccomplex_t IT_0034 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0035 = IT_0017*IT_0034;
    const ccomplex_t IT_0036 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = -IT_0031 + IT_0035 + IT_0037;
    const ccomplex_t IT_0039 = (-3)*s_24;
    const ccomplex_t IT_0040 = -IT_0037;
    const ccomplex_t IT_0041 = (-3)*s_34;
    const ccomplex_t IT_0042 = (-3)*s_13;
    const ccomplex_t IT_0043 = (-3)*s_23;
    const ccomplex_t IT_0044 = (-0.333333333333333)*conj(IT_0029);
    const ccomplex_t IT_0045 = (-0.333333333333333)*conj(IT_0040);
    const ccomplex_t IT_0046 = IT_0029*(conj(IT_0032)*IT_0033 + conj(IT_0038)
      *IT_0039 + conj(IT_0040)*IT_0041) + (-3)*IT_0040*(IT_0013*conj(IT_0040) + 
      (-0.333333333333333)*conj(IT_0032)*IT_0042 + (-0.333333333333333)*conj
      (IT_0038)*IT_0043 + IT_0041*IT_0044) + (-3)*IT_0032*(IT_0014*conj(IT_0032)
       + s_12*conj(IT_0038) + IT_0033*IT_0044 + IT_0042*IT_0045) + (-3)*IT_0038*
      (s_12*conj(IT_0032) + IT_0019*conj(IT_0038) + IT_0039*IT_0044 + IT_0043
      *IT_0045);
    return create_ccomplex_return(IT_0046);
}

