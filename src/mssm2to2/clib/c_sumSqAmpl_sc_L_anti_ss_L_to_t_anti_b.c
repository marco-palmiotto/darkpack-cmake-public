#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_t_anti_b(
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
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = pow(m_sc_L, 2);
    const ccomplex_t IT_0007 = pow(m_ss_L, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0002;
    const ccomplex_t IT_0014 = m_t*IT_0013;
    const ccomplex_t IT_0015 = m_t*IT_0003;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0014 + IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = sin(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0019*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = pow(m_s, 2);
    const ccomplex_t IT_0026 = cpow(IT_0020, -1);
    const ccomplex_t IT_0027 = pow(m_c, 2);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001*(IT_0019*IT_0021*IT_0025*IT_0026 + (-2)*IT_0020*(m_W*IT_0021
       + (-0.5)*IT_0019*IT_0022*IT_0027));
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0019;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0020, 2);
    const ccomplex_t IT_0037 = cpow(IT_0021, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*(-1.4142135623731))*conj(V_cs)
      *e_em*IT_0001*(IT_0019*(IT_0025 + -IT_0027) + -m_W*(IT_0036 + -IT_0037));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0035*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, 2);
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + IT_0041 
      + -reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = -IT_0018 + -IT_0033 + -IT_0044;
    const ccomplex_t IT_0046 = m_t*s_14;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0019*IT_0021*IT_0026;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0029*IT_0048;
    const ccomplex_t IT_0050 = IT_0031*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = m_b*IT_0013;
    const ccomplex_t IT_0053 = m_b*IT_0003;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = IT_0010*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0019;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0039*IT_0058;
    const ccomplex_t IT_0060 = IT_0042*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0051 + IT_0056 + -IT_0061;
    const ccomplex_t IT_0063 = m_b*s_13;
    const ccomplex_t IT_0064 = s_13*s_14;
    const ccomplex_t IT_0065 = s_34*IT_0006;
    const ccomplex_t IT_0066 = m_b*m_t;
    const ccomplex_t IT_0067 = 18*s_34;
    const ccomplex_t IT_0068 = conj(IT_0012)*(18*IT_0045*IT_0046 + (-18)
      *IT_0062*IT_0063 + IT_0012*(36*IT_0064 + (-18)*IT_0065)) + conj(IT_0045)*
      (18*IT_0012*IT_0046 + (-18)*IT_0062*IT_0066 + IT_0045*IT_0067) + conj
      (IT_0062)*((-18)*IT_0012*IT_0063 + (-18)*IT_0045*IT_0066 + IT_0062*IT_0067);
    return create_ccomplex_return(IT_0068);
}

