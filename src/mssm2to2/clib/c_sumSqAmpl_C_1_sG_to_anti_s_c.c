#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sG_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sG_to_anti_s_c(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = m_c*m_s*m_sG*m_C_1;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0005*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = pow(m_sG, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0010 + IT_0020 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 8*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = 16*IT_0027;
    const ccomplex_t IT_0029 = IT_0024 + IT_0028;
    const ccomplex_t IT_0030 = 8*conj(IT_0023);
    const ccomplex_t IT_0031 = 16*conj(IT_0027);
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0010 + IT_0020 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 8*IT_0027;
    const ccomplex_t IT_0037 = 16*IT_0023;
    const ccomplex_t IT_0038 = IT_0036 + IT_0037;
    const ccomplex_t IT_0039 = 8*conj(IT_0027);
    const ccomplex_t IT_0040 = 16*conj(IT_0023);
    const ccomplex_t IT_0041 = IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = s_14*s_23;
    const ccomplex_t IT_0043 = s_13*s_24;
    const ccomplex_t IT_0044 = s_12*s_34;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0042 + IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = 8*IT_0035;
    const ccomplex_t IT_0048 = 8*conj(IT_0035);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = IT_0021*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = m_s*m_sG*s_14;
    const ccomplex_t IT_0055 = m_c*s_23*m_C_1;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0008*IT_0057;
    const ccomplex_t IT_0059 = IT_0012*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = m_sG*s_34*m_C_1;
    const ccomplex_t IT_0063 = 8*IT_0062;
    const ccomplex_t IT_0064 = m_c*m_sG*s_13;
    const ccomplex_t IT_0065 = m_s*s_24*m_C_1;
    const ccomplex_t IT_0066 = 16*IT_0035;
    const ccomplex_t IT_0067 = 16*IT_0014;
    const ccomplex_t IT_0068 = IT_0000*(conj(IT_0014)*IT_0029 + IT_0014
      *IT_0032 + conj(IT_0035)*IT_0038 + IT_0035*IT_0041) + IT_0046*(IT_0024
      *conj(IT_0027) + IT_0027*IT_0030 + conj(IT_0014)*IT_0047 + IT_0014*IT_0048
      ) + IT_0053*(16*IT_0042*conj(IT_0053) + 8*(conj(IT_0014) + 2*conj(IT_0035)
      )*IT_0054 + IT_0041*IT_0055 + conj(IT_0061)*IT_0063) + IT_0061*(16*IT_0043
      *conj(IT_0061) + conj(IT_0053)*IT_0063 + IT_0032*IT_0064 + (16*conj
      (IT_0014) + IT_0048)*IT_0065) + IT_0042*(conj(IT_0023)*IT_0037 + conj
      (IT_0035)*IT_0066) + conj(IT_0053)*(IT_0038*IT_0055 + 8*IT_0054*(IT_0014 +
       0.125*IT_0066)) + IT_0043*(conj(IT_0027)*IT_0028 + conj(IT_0014)*IT_0067)
       + conj(IT_0061)*(IT_0029*IT_0064 + IT_0065*(IT_0047 + IT_0067));
    return create_ccomplex_return(IT_0068);
}

