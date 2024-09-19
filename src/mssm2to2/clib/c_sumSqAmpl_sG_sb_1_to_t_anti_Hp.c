#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_1_to_t_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_1_to_t_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0000*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = pow(m_sG, 2);
    const ccomplex_t IT_0012 = pow(m_sb_1, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0015 = IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0017 = cpow(IT_0001, -1);
    const ccomplex_t IT_0018 = m_t*U_sb_00*(m_t*IT_0001*IT_0003*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0002*IT_0017*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0003*IT_0017*conj(U_st_11)));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0002*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0000*IT_0018 +
       m_W*IT_0001*IT_0002*U_sb_00*conj(U_st_01)) + A_t*IT_0001*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0011 + IT_0021 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = IT_0013*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = m_sG*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0031 = m_t*U_sb_00*(m_t*IT_0001*IT_0003*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0002*IT_0017*U_sb_00*conj(U_st_00) +
       U_sb_10*(mu_h*conj(U_st_00) + m_t*IT_0003*IT_0017*conj(U_st_10)));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0002*U_sb_10
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0000*IT_0031 +
       m_W*IT_0001*IT_0002*U_sb_00*conj(U_st_00)) + A_t*IT_0001*U_sb_00*conj
      (U_st_10));
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0011 + IT_0021 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0015 + IT_0024 + IT_0029 + IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0039 = IT_0032*IT_0038;
    const ccomplex_t IT_0040 = IT_0034*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0043 = IT_0019*IT_0042;
    const ccomplex_t IT_0044 = IT_0022*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0000*IT_0002*IT_0005*IT_0017;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0025*IT_0047;
    const ccomplex_t IT_0049 = IT_0013*IT_0014*IT_0048;
    const ccomplex_t IT_0050 = IT_0008*IT_0047;
    const ccomplex_t IT_0051 = IT_0013*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = m_sG*IT_0052;
    const ccomplex_t IT_0054 = -IT_0041 + -IT_0045 + IT_0049 + -IT_0053;
    const ccomplex_t IT_0055 = 8*s_13;
    const ccomplex_t IT_0056 = m_t*m_sG;
    const ccomplex_t IT_0057 = 8*IT_0037;
    const ccomplex_t IT_0058 = -IT_0052;
    const ccomplex_t IT_0059 = m_t*s_12;
    const ccomplex_t IT_0060 = m_sG*s_23;
    const ccomplex_t IT_0061 = 8*conj(IT_0037);
    const ccomplex_t IT_0062 = 8*conj(IT_0054);
    const ccomplex_t IT_0063 = m_t*m_sG*IT_0012;
    const ccomplex_t IT_0064 = s_13*IT_0012;
    const ccomplex_t IT_0065 = (-8)*IT_0064;
    const ccomplex_t IT_0066 = s_12*s_23;
    const ccomplex_t IT_0067 = 16*IT_0066;
    const ccomplex_t IT_0068 = IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = 8*IT_0054;
    const ccomplex_t IT_0070 = 0.125*IT_0059;
    const ccomplex_t IT_0071 = 0.125*IT_0060;
    const ccomplex_t IT_0072 = (IT_0037*conj(IT_0037) + IT_0054*conj(IT_0054))
      *IT_0055 + conj(IT_0054)*IT_0056*IT_0057 + (IT_0054*IT_0056 + IT_0058
      *IT_0059 + IT_0028*IT_0060)*IT_0061 + (IT_0028*IT_0059 + IT_0058*IT_0060)
      *IT_0062 + 8*conj(IT_0028)*(IT_0058*IT_0063 + 0.125*IT_0028*IT_0068 +
       IT_0069*IT_0070 + IT_0057*IT_0071) + 8*conj(IT_0058)*(IT_0028*IT_0063 +
       0.125*IT_0058*IT_0068 + IT_0057*IT_0070 + IT_0069*IT_0071);
    return create_ccomplex_return(IT_0072);
}

