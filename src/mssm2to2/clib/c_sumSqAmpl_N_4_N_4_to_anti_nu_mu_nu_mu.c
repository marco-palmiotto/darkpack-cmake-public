#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_nu_mu_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_nu_mu_nu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = conj(N_B4)*e_em;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 1.4142135623731*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_W4)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0005 + -IT_0010);
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = N_B4*e_em;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = N_W4*e_em;
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0015 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0012*IT_0020;
    const ccomplex_t IT_0022 = pow(m_N_4, 2);
    const ccomplex_t IT_0023 = cpow(s_13 + (-0.5)*IT_0022 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0006;
    const ccomplex_t IT_0028 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0001*IT_0007;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0029 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = e_em*IT_0027;
    const ccomplex_t IT_0038 = e_em*IT_0030;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0037 + IT_0038);
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0036*IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0035;
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = conj(IT_0045) + conj(IT_0048);
    const ccomplex_t IT_0050 = 4*IT_0026;
    const ccomplex_t IT_0051 = IT_0045 + IT_0048;
    const ccomplex_t IT_0052 = 8*conj(IT_0026);
    const ccomplex_t IT_0053 = s_34*IT_0022;
    const ccomplex_t IT_0054 = cpow(s_23 + (-0.5)*IT_0022 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0055 = IT_0021*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = 4*IT_0057;
    const ccomplex_t IT_0059 = 4*conj(IT_0057);
    const ccomplex_t IT_0060 = 0.5*IT_0044;
    const ccomplex_t IT_0061 = (-0.5)*IT_0047;
    const ccomplex_t IT_0062 = conj(IT_0060) + conj(IT_0061);
    const ccomplex_t IT_0063 = 8*IT_0045;
    const ccomplex_t IT_0064 = IT_0060 + IT_0061;
    const ccomplex_t IT_0065 = 4*conj(IT_0026);
    const ccomplex_t IT_0066 = 8*conj(IT_0045);
    const ccomplex_t IT_0067 = conj(IT_0026)*IT_0057 + IT_0026*conj(IT_0057) +
       0.5*IT_0049*IT_0058 + 0.5*IT_0051*IT_0059 + 4*IT_0062*(IT_0048 + 0.125
      *IT_0050 + 0.125*IT_0063) + 4*IT_0064*(conj(IT_0048) + 0.125*IT_0065 +
       0.125*IT_0066);
    const ccomplex_t IT_0068 = s_14*s_23;
    const ccomplex_t IT_0069 = 8*conj(IT_0057);
    const ccomplex_t IT_0070 = 8*IT_0057;
    const ccomplex_t IT_0071 = 16*IT_0060;
    const ccomplex_t IT_0072 = 8*IT_0000*((IT_0026 + 2*IT_0045 + 2*IT_0048)
      *IT_0049 + 0.125*conj(IT_0026)*IT_0050 + 0.125*IT_0051*IT_0052) + 2
      *IT_0053*IT_0067 + IT_0068*(conj(IT_0057)*IT_0058 + IT_0064*IT_0069 +
       IT_0062*(16*IT_0061 + IT_0070 + IT_0071));
    return create_ccomplex_return(IT_0072);
}

