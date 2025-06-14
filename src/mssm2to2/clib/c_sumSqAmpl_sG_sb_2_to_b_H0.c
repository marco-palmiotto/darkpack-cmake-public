#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_2_to_b_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_2_to_b_H0(
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = pow(m_sG, 2);
    const ccomplex_t IT_0012 = pow(m_sb_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = m_b*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0018 = cpow(V_tb, 2);
    const ccomplex_t IT_0019 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0020 = sin(alpha);
    const ccomplex_t IT_0021 = sin(beta);
    const ccomplex_t IT_0022 = -IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0003;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = cos(theta_W);
    const ccomplex_t IT_0026 = cpow(IT_0025, -2);
    const ccomplex_t IT_0027 = (-2)*conj(U_sb_00)*((IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0010 + (-0.5)*m_W*IT_0024*(IT_0005 + 1./3*IT_0004*IT_0026))
      *U_sb_01 + -1./2*IT_0002*IT_0005*IT_0019*IT_0020*U_sb_11) + conj(U_sb_10)*
      (IT_0002*IT_0005*IT_0019*IT_0020*U_sb_01 + (-2)*(IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0010 + -1./3*m_W*IT_0004*IT_0024*IT_0026)*U_sb_11);
    const ccomplex_t IT_0028 = (0 + _Complex_I*-1)*(e_em*IT_0027 + (
      -1.4142135623731)*A_b*IT_0003*IT_0018*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0017*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0035 = (-2)*U_sb_01*((IT_0000*IT_0002*IT_0003*IT_0005
      *IT_0010 + (-0.5)*m_W*IT_0024*(IT_0005 + 1./3*IT_0004*IT_0026))*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0019*IT_0020*conj(U_sb_11)) + U_sb_11
      *(IT_0002*IT_0005*IT_0019*IT_0020*conj(U_sb_01) + (-2)*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0010 + -1./3*m_W*IT_0004*IT_0024*IT_0026)*conj(U_sb_11
      ));
    const ccomplex_t IT_0036 = (0 + _Complex_I*-1)*(e_em*IT_0035 + (
      -1.4142135623731)*A_b*IT_0003*IT_0018*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0008*IT_0037;
    const ccomplex_t IT_0039 = IT_0034*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (-2)*IT_0016 + IT_0033 + IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = IT_0013*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = m_b*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0048 = (-2)*conj(U_sb_00)*((IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0010 + (-0.5)*m_W*IT_0024*(IT_0005 + 1./3*IT_0004*IT_0026))
      *U_sb_01 + -1./2*IT_0002*IT_0005*IT_0019*IT_0020*U_sb_11) + conj(U_sb_10)*
      (IT_0002*IT_0005*IT_0019*IT_0020*U_sb_01 + (-2)*(IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0010 + -1./3*m_W*IT_0004*IT_0024*IT_0026)*U_sb_11);
    const ccomplex_t IT_0049 = (0 + _Complex_I*-1)*(e_em*IT_0048 + (
      -1.4142135623731)*A_b*IT_0003*IT_0018*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = IT_0031*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = (-2)*U_sb_01*((IT_0000*IT_0002*IT_0003*IT_0005
      *IT_0010 + (-0.5)*m_W*IT_0024*(IT_0005 + 1./3*IT_0004*IT_0026))*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0019*IT_0020*conj(U_sb_11)) + U_sb_11
      *(IT_0002*IT_0005*IT_0019*IT_0020*conj(U_sb_01) + (-2)*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0010 + -1./3*m_W*IT_0004*IT_0024*IT_0026)*conj(U_sb_11
      ));
    const ccomplex_t IT_0055 = (0 + _Complex_I*-1)*(e_em*IT_0054 + (
      -1.4142135623731)*A_b*IT_0003*IT_0018*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0042*IT_0056;
    const ccomplex_t IT_0058 = IT_0034*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = 2*IT_0046 + -IT_0053 + -IT_0059;
    const ccomplex_t IT_0061 = 8*s_13;
    const ccomplex_t IT_0062 = m_b*m_sG;
    const ccomplex_t IT_0063 = 8*IT_0041;
    const ccomplex_t IT_0064 = -IT_0015;
    const ccomplex_t IT_0065 = m_b*s_14;
    const ccomplex_t IT_0066 = m_sG*s_34;
    const ccomplex_t IT_0067 = 8*conj(IT_0041);
    const ccomplex_t IT_0068 = 8*conj(IT_0060);
    const ccomplex_t IT_0069 = pow(m_H0, 2);
    const ccomplex_t IT_0070 = m_b*m_sG*IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0069;
    const ccomplex_t IT_0072 = s_14*s_34;
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = IT_0071 + IT_0073;
    const ccomplex_t IT_0075 = 8*IT_0060;
    const ccomplex_t IT_0076 = 0.125*IT_0066;
    const ccomplex_t IT_0077 = (IT_0041*conj(IT_0041) + IT_0060*conj(IT_0060))
      *IT_0061 + conj(IT_0060)*IT_0062*IT_0063 + (IT_0060*IT_0062 + IT_0064
      *IT_0065 + IT_0045*IT_0066)*IT_0067 + (IT_0045*IT_0065 + IT_0064*IT_0066)
      *IT_0068 + 8*conj(IT_0045)*(IT_0064*IT_0070 + -IT_0045*IT_0074 + 0.125
      *IT_0065*IT_0075 + IT_0063*IT_0076) + 8*conj(IT_0064)*(0.125*IT_0063
      *IT_0065 + IT_0045*IT_0070 + -IT_0064*IT_0074 + IT_0075*IT_0076);
    return create_ccomplex_return(IT_0077);
}

