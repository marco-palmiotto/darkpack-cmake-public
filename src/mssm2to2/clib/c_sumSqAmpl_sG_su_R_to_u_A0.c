#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_su_R_to_u_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_su_R_to_u_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = m_u*e_em*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = pow(m_u, 2);
    const ccomplex_t IT_0008 = pow(m_sG, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = pow(m_su_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0013 
      + -reg_prop, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = m_u*e_em*IT_0000*IT_0002*IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0005*IT_0019;
    const ccomplex_t IT_0021 = IT_0014*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = m_sG*s_34;
    const ccomplex_t IT_0025 = pow(m_A0, 2);
    const ccomplex_t IT_0026 = s_13*IT_0025;
    const ccomplex_t IT_0027 = s_14*s_34;
    return create_ccomplex_return(8*s_13*IT_0012*conj(IT_0012) + 8*IT_0023*(
      (IT_0012 + conj(IT_0012))*IT_0024 + -IT_0023*(IT_0026 + (-2)*IT_0027)));
}

