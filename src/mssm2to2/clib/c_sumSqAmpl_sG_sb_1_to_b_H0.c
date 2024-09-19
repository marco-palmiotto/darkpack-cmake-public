#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_1_to_b_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_1_to_b_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0001 = cpow(V_tb, 2);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = m_b*mu_h*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = -IT_0007*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0005;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_b, 2);
    const ccomplex_t IT_0017 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0013*(IT_0009 + 1.
      /3*IT_0008*IT_0015) + IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj
      (U_sb_01) + -1./2*IT_0004*IT_0006*IT_0007*IT_0009*conj(U_sb_11)) + U_sb_10
      *(IT_0004*IT_0006*IT_0007*IT_0009*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0008*IT_0013*IT_0015 + (-3)*IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)
      *conj(U_sb_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_b*IT_0001*IT_0002*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0019;
    const ccomplex_t IT_0021 = pow(m_sG, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0016 + IT_0021 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_b*e_em*IT_0002*IT_0003
      *IT_0006*IT_0009;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = pow(m_sb_1, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + IT_0016 + -IT_0021 + -IT_0029 
      + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = m_b*IT_0032;
    const ccomplex_t IT_0034 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0013*(IT_0009 + 1.
      /3*IT_0008*IT_0015) + IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj
      (U_sb_00) + -1./2*IT_0004*IT_0006*IT_0007*IT_0009*conj(U_sb_10)) + U_sb_10
      *(IT_0004*IT_0006*IT_0007*IT_0009*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0008*IT_0013*IT_0015 + (-3)*IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)
      *conj(U_sb_10));
    const ccomplex_t IT_0035 = (0 + _Complex_I*-1)*(e_em*IT_0034 + (
      -1.4142135623731)*A_b*IT_0001*IT_0002*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0016 + IT_0021 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0024 + (-2)*IT_0033 + IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0043 = IT_0026*IT_0042;
    const ccomplex_t IT_0044 = IT_0030*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = m_b*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0048 = IT_0019*IT_0047;
    const ccomplex_t IT_0049 = IT_0022*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0036*IT_0042;
    const ccomplex_t IT_0052 = IT_0038*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 2*IT_0046 + -IT_0050 + -IT_0053;
    const ccomplex_t IT_0055 = 8*s_13;
    const ccomplex_t IT_0056 = m_b*m_sG;
    const ccomplex_t IT_0057 = 8*IT_0041;
    const ccomplex_t IT_0058 = -IT_0032;
    const ccomplex_t IT_0059 = m_b*s_14;
    const ccomplex_t IT_0060 = m_sG*s_34;
    const ccomplex_t IT_0061 = 8*conj(IT_0041);
    const ccomplex_t IT_0062 = 8*conj(IT_0054);
    const ccomplex_t IT_0063 = pow(m_H0, 2);
    const ccomplex_t IT_0064 = m_b*m_sG*IT_0063;
    const ccomplex_t IT_0065 = s_13*IT_0063;
    const ccomplex_t IT_0066 = s_14*s_34;
    const ccomplex_t IT_0067 = (-2)*IT_0066;
    const ccomplex_t IT_0068 = IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = 8*IT_0054;
    const ccomplex_t IT_0070 = 0.125*IT_0060;
    const ccomplex_t IT_0071 = (IT_0041*conj(IT_0041) + IT_0054*conj(IT_0054))
      *IT_0055 + conj(IT_0054)*IT_0056*IT_0057 + (IT_0054*IT_0056 + IT_0058
      *IT_0059 + IT_0045*IT_0060)*IT_0061 + (IT_0045*IT_0059 + IT_0058*IT_0060)
      *IT_0062 + 8*conj(IT_0045)*(IT_0058*IT_0064 + -IT_0045*IT_0068 + 0.125
      *IT_0059*IT_0069 + IT_0057*IT_0070) + 8*conj(IT_0058)*(0.125*IT_0057
      *IT_0059 + IT_0045*IT_0064 + -IT_0058*IT_0068 + IT_0069*IT_0070);
    return create_ccomplex_return(IT_0071);
}

