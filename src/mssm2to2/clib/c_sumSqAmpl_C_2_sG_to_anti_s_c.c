#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sG_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sG_to_anti_s_c(
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
    const creal_t m_C_2 = param->m_C_2;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = pow(m_s, 2);
    const ccomplex_t IT_0006 = pow(m_sG, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = m_sG*s_34*m_C_2;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0005 + IT_0015 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = 8*conj(IT_0019);
    const ccomplex_t IT_0021 = s_13*s_24;
    const ccomplex_t IT_0022 = 16*IT_0019;
    const ccomplex_t IT_0023 = pow(m_W, -1);
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0001*IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = IT_0016*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = m_s*m_sG*s_14;
    const ccomplex_t IT_0032 = 8*IT_0010;
    const ccomplex_t IT_0033 = m_s*s_24*m_C_2;
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0001*IT_0023*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0005 + IT_0015 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = m_c*s_23*m_C_2;
    const ccomplex_t IT_0043 = 8*conj(IT_0010);
    const ccomplex_t IT_0044 = m_c*m_sG*s_13;
    const ccomplex_t IT_0045 = 16*conj(IT_0019);
    const ccomplex_t IT_0046 = 8*IT_0019;
    const ccomplex_t IT_0047 = s_14*s_23;
    const ccomplex_t IT_0048 = 16*IT_0010;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0028*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0007*IT_0038;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 16*conj(IT_0010);
    const ccomplex_t IT_0055 = s_12*s_34;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0021 + IT_0047 + IT_0056;
    const ccomplex_t IT_0058 = 16*IT_0053;
    const ccomplex_t IT_0059 = 8*IT_0053;
    const ccomplex_t IT_0060 = 16*IT_0041;
    const ccomplex_t IT_0061 = 8*conj(IT_0051);
    const ccomplex_t IT_0062 = m_c*m_s*m_sG*m_C_2;
    const ccomplex_t IT_0063 = 16*(conj(IT_0041)*(IT_0030 + 1./2*IT_0051) + 1.
      /2*(IT_0030 + 2*IT_0051)*conj(IT_0053) + 0.0625*conj(IT_0030)*(IT_0059 +
       IT_0060) + 0.0625*(IT_0041 + 0.125*IT_0058)*IT_0061)*IT_0062;
    const ccomplex_t IT_0064 = IT_0010*IT_0011*IT_0020 + conj(IT_0019)*IT_0021
      *IT_0022 + conj(IT_0030)*(IT_0031*IT_0032 + IT_0022*IT_0033) + IT_0041*
      (IT_0042*IT_0043 + IT_0044*IT_0045) + conj(IT_0010)*(IT_0011*IT_0046 +
       IT_0047*IT_0048) + (IT_0033*IT_0046 + IT_0031*IT_0048)*conj(IT_0051) +
       IT_0053*(IT_0020*IT_0044 + IT_0042*IT_0054) + IT_0051*(IT_0020*IT_0033 +
       16*IT_0047*conj(IT_0051) + IT_0031*IT_0054 + 8*conj(IT_0030)*IT_0057) + 8
      *conj(IT_0053)*(0.125*IT_0044*IT_0046 + 0.125*IT_0042*IT_0048 + IT_0041
      *IT_0057 + 0.125*IT_0047*IT_0058) + conj(IT_0041)*(IT_0032*IT_0042 +
       IT_0022*IT_0044 + IT_0057*IT_0059 + IT_0021*IT_0060) + IT_0030*(16
      *IT_0021*conj(IT_0030) + IT_0031*IT_0043 + IT_0033*IT_0045 + IT_0057
      *IT_0061) + IT_0063;
    return create_ccomplex_return(IT_0064);
}

