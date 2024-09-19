#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sc_L_to_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sc_L_to_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = pow(m_sG, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(m_sc_L, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0008 + IT_0012 + -IT_0013 
      + -reg_prop, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0000*IT_0002*IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0005*IT_0019;
    const ccomplex_t IT_0021 = IT_0014*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = m_sG*IT_0022;
    const ccomplex_t IT_0024 = IT_0011 + IT_0023;
    const ccomplex_t IT_0025 = 8*s_13;
    const ccomplex_t IT_0026 = cpow(IT_0017, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0002*(IT_0000*IT_0007*IT_0016*IT_0017 + (-2)*IT_0015*(m_W*IT_0017
       + (-0.5)*IT_0000*IT_0012*IT_0026));
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0000*IT_0002*IT_0015*IT_0026;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0005*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0037 = IT_0014*IT_0035*IT_0036;
    const ccomplex_t IT_0038 = -IT_0032 + IT_0037;
    const ccomplex_t IT_0039 = m_s*s_12;
    const ccomplex_t IT_0040 = m_s*m_sG;
    const ccomplex_t IT_0041 = 8*IT_0024;
    const ccomplex_t IT_0042 = s_12*s_23;
    const ccomplex_t IT_0043 = s_13*IT_0013;
    const ccomplex_t IT_0044 = m_sG*s_23;
    const ccomplex_t IT_0045 = 8*conj(IT_0024);
    const ccomplex_t IT_0046 = IT_0024*conj(IT_0024)*IT_0025 + 8*conj(IT_0038)
      *(0.125*IT_0025*IT_0038 + IT_0022*IT_0039 + 0.125*IT_0040*IT_0041) + 8
      *conj(IT_0022)*(IT_0038*IT_0039 + 0.125*IT_0022*(16*IT_0042 + (-8)*IT_0043
      ) + 0.125*IT_0041*IT_0044) + (IT_0038*IT_0040 + IT_0022*IT_0044)*IT_0045;
    return create_ccomplex_return(IT_0046);
}

