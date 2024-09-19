#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sG_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sG_to_u_anti_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = pow(m_u, 2);
    const ccomplex_t IT_0006 = pow(m_C_2, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = m_sG*s_34*m_C_2;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_sG, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0005 + IT_0015 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = 8*conj(IT_0019);
    const ccomplex_t IT_0021 = s_14*s_23;
    const ccomplex_t IT_0022 = 16*IT_0019;
    const ccomplex_t IT_0023 = m_d*s_23*m_C_2;
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0001*IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = 8*IT_0019;
    const ccomplex_t IT_0034 = IT_0016*IT_0029;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0024*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0005 + IT_0015 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = m_u*s_24*m_C_2;
    const ccomplex_t IT_0045 = 8*conj(IT_0010);
    const ccomplex_t IT_0046 = m_u*m_sG*s_14;
    const ccomplex_t IT_0047 = 16*conj(IT_0019);
    const ccomplex_t IT_0048 = s_13*s_24;
    const ccomplex_t IT_0049 = s_12*s_34;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0021 + IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = m_d*m_sG*s_13;
    const ccomplex_t IT_0053 = 16*IT_0010;
    const ccomplex_t IT_0054 = 8*IT_0010;
    const ccomplex_t IT_0055 = IT_0007*IT_0040;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = 16*conj(IT_0010);
    const ccomplex_t IT_0058 = 16*IT_0056;
    const ccomplex_t IT_0059 = 8*IT_0056;
    const ccomplex_t IT_0060 = 16*IT_0043;
    const ccomplex_t IT_0061 = m_d*m_u*m_sG*m_C_2;
    const ccomplex_t IT_0062 = IT_0010*IT_0011*IT_0020 + conj(IT_0019)*IT_0021
      *IT_0022 + IT_0023*(conj(IT_0032)*IT_0033 + IT_0022*conj(IT_0035)) +
       IT_0043*(IT_0044*IT_0045 + IT_0046*IT_0047) + IT_0035*(16*IT_0021*conj
      (IT_0035) + IT_0023*IT_0047 + 8*conj(IT_0032)*IT_0051 + IT_0045*IT_0052) +
       conj(IT_0010)*(IT_0011*IT_0033 + IT_0048*IT_0053) + IT_0052*(conj(IT_0032
      )*IT_0053 + conj(IT_0035)*IT_0054) + IT_0056*(IT_0020*IT_0046 + IT_0044
      *IT_0057) + IT_0032*(IT_0020*IT_0023 + 16*conj(IT_0032)*IT_0048 + 8*conj
      (IT_0035)*IT_0051 + IT_0052*IT_0057) + 8*conj(IT_0056)*(0.125*IT_0033
      *IT_0046 + IT_0043*IT_0051 + 0.125*IT_0044*IT_0053 + 0.125*IT_0048*IT_0058
      ) + conj(IT_0043)*(IT_0022*IT_0046 + IT_0044*IT_0054 + IT_0051*IT_0059 +
       IT_0021*IT_0060) + 16*((0.5*IT_0032 + IT_0035)*conj(IT_0043) + (IT_0032 +
       1./2*IT_0035)*conj(IT_0056) + 1./2*conj(IT_0032)*(IT_0043 + 0.125*IT_0058
      ) + 0.0625*conj(IT_0035)*(IT_0059 + IT_0060))*IT_0061;
    return create_ccomplex_return(IT_0062);
}

