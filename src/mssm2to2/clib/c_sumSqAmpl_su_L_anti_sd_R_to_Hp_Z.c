#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_R_to_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_R_to_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_sd_R, 2);
    const ccomplex_t IT_0009 = pow(m_su_L, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + (-3)
      *IT_0000*IT_0003);
    const ccomplex_t IT_0015 = 0.333333333333333*IT_0014;
    const ccomplex_t IT_0016 = pow(m_Hp, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + -IT_0020;
    const ccomplex_t IT_0022 = pow(m_Z, -2);
    const ccomplex_t IT_0023 = pow(s_34, 2);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0026 = (-0.666666666666667)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0009 + IT_0016 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0007*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = IT_0020 + IT_0030;
    const ccomplex_t IT_0032 = s_24*s_34*IT_0022;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = s_23 + IT_0033;
    const ccomplex_t IT_0035 = pow(s_24, 2);
    const ccomplex_t IT_0036 = IT_0022*IT_0035;
    return create_ccomplex_return((-3)*IT_0021*(conj(IT_0021)*(IT_0016 + 
      -IT_0024) + conj(IT_0031)*IT_0034) + (-3)*IT_0031*(conj(IT_0021)*IT_0034 +
       conj(IT_0031)*(IT_0008 + -IT_0036)));
}

