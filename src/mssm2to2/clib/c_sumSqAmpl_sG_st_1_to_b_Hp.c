#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_1_to_b_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_1_to_b_Hp(
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = pow(m_sG, 2);
    const ccomplex_t IT_0012 = pow(m_st_1, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0010*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0005, -1);
    const ccomplex_t IT_0016 = m_t*conj(U_sb_00)*(m_t*IT_0003*IT_0015*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0004*IT_0005*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0004*IT_0015*U_st_10));
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002
      *IT_0016 + m_W*IT_0003*IT_0005*conj(U_sb_00)*U_st_00) + A_t*IT_0003*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = pow(m_b, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0011 + IT_0020 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = m_sG*IT_0027;
    const ccomplex_t IT_0029 = m_t*conj(U_sb_01)*(m_t*IT_0003*IT_0015*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0004*IT_0005*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0004*IT_0015*U_st_10));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002
      *IT_0029 + m_W*IT_0003*IT_0005*conj(U_sb_01)*U_st_00) + A_t*IT_0003*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0011 + IT_0020 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0014 + IT_0023 + IT_0028 + IT_0035;
    const ccomplex_t IT_0037 = 8*s_13;
    const ccomplex_t IT_0038 = m_b*s_12;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0002*IT_0003*IT_0007*IT_0015;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = m_sG*s_23;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0047 = IT_0030*IT_0046;
    const ccomplex_t IT_0048 = IT_0033*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0024*IT_0040;
    const ccomplex_t IT_0051 = IT_0000*IT_0013*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0053 = IT_0017*IT_0052;
    const ccomplex_t IT_0054 = IT_0021*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_sG*IT_0043;
    const ccomplex_t IT_0057 = -IT_0049 + IT_0051 + -IT_0055 + -IT_0056;
    const ccomplex_t IT_0058 = 8*conj(IT_0057);
    const ccomplex_t IT_0059 = m_b*m_sG;
    const ccomplex_t IT_0060 = 8*IT_0036;
    const ccomplex_t IT_0061 = 8*conj(IT_0036);
    const ccomplex_t IT_0062 = m_b*m_sG*IT_0012;
    const ccomplex_t IT_0063 = s_13*IT_0012;
    const ccomplex_t IT_0064 = (-8)*IT_0063;
    const ccomplex_t IT_0065 = s_12*s_23;
    const ccomplex_t IT_0066 = 16*IT_0065;
    const ccomplex_t IT_0067 = IT_0064 + IT_0066;
    const ccomplex_t IT_0068 = 0.125*IT_0038;
    const ccomplex_t IT_0069 = 8*IT_0057;
    const ccomplex_t IT_0070 = 0.125*IT_0045;
    const ccomplex_t IT_0071 = IT_0036*conj(IT_0036)*IT_0037 + (IT_0027
      *IT_0038 + IT_0044*IT_0045)*IT_0058 + conj(IT_0057)*(IT_0037*IT_0057 +
       IT_0059*IT_0060) + (IT_0038*IT_0044 + IT_0027*IT_0045 + IT_0057*IT_0059)
      *IT_0061 + 8*conj(IT_0027)*(IT_0044*IT_0062 + 0.125*IT_0027*IT_0067 +
       IT_0068*IT_0069 + IT_0060*IT_0070) + 8*conj(IT_0044)*(IT_0027*IT_0062 +
       0.125*IT_0044*IT_0067 + IT_0060*IT_0068 + IT_0069*IT_0070);
    return create_ccomplex_return(IT_0071);
}

