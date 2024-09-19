#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_sd_L_to_e_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_sd_L_to_e_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_e, 2);
    const ccomplex_t IT_0010 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + IT_0020;
    const ccomplex_t IT_0022 = sin(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0003*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0003*IT_0023;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = IT_0011*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0029 + IT_0035;
    const ccomplex_t IT_0037 = s_34*IT_0010;
    const ccomplex_t IT_0038 = (-6)*IT_0037;
    const ccomplex_t IT_0039 = s_13*s_14;
    const ccomplex_t IT_0040 = 12*IT_0039;
    const ccomplex_t IT_0041 = IT_0038 + IT_0040;
    const ccomplex_t IT_0042 = m_e*m_u;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0044 = IT_0007*IT_0031;
    const ccomplex_t IT_0045 = IT_0011*IT_0043*IT_0044;
    const ccomplex_t IT_0046 = m_e*IT_0029;
    const ccomplex_t IT_0047 = m_e*IT_0035;
    const ccomplex_t IT_0048 = -IT_0045 + IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0050 = IT_0016*IT_0025;
    const ccomplex_t IT_0051 = IT_0018*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = IT_0014*IT_0026;
    const ccomplex_t IT_0053 = IT_0018*IT_0049*IT_0052;
    const ccomplex_t IT_0054 = m_e*IT_0013;
    const ccomplex_t IT_0055 = IT_0002*IT_0032;
    const ccomplex_t IT_0056 = IT_0011*IT_0043*IT_0055;
    const ccomplex_t IT_0057 = m_e*IT_0020;
    const ccomplex_t IT_0058 = IT_0054 + -IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = (-6)*IT_0048;
    const ccomplex_t IT_0060 = IT_0051 + 0.166666666666667*IT_0059;
    const ccomplex_t IT_0061 = (-6)*conj(IT_0048);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = conj(IT_0051) + IT_0062;
    const ccomplex_t IT_0064 = m_e*s_14;
    const ccomplex_t IT_0065 = -IT_0053 + IT_0058;
    const ccomplex_t IT_0066 = -conj(IT_0053) + conj(IT_0058);
    const ccomplex_t IT_0067 = 6*s_34;
    const ccomplex_t IT_0068 = (-6)*s_34;
    const ccomplex_t IT_0069 = m_e*m_u*IT_0010;
    const ccomplex_t IT_0070 = (-6)*IT_0021;
    const ccomplex_t IT_0071 = (-6)*conj(IT_0021);
    const ccomplex_t IT_0072 = m_u*s_13;
    const ccomplex_t IT_0073 = (IT_0021*conj(IT_0021) + IT_0036*conj(IT_0036))
      *IT_0041 + 6*IT_0042*((conj(IT_0048) + -conj(IT_0051))*IT_0053 + (IT_0048 
      + -IT_0051)*conj(IT_0053) + conj(IT_0058)*IT_0060 + IT_0058*IT_0063) + 6
      *IT_0064*(conj(IT_0036)*IT_0060 + IT_0036*IT_0063 + -conj(IT_0021)*IT_0065
       + -IT_0021*IT_0066) + (IT_0048*conj(IT_0048) + IT_0051*conj(IT_0051) +
       IT_0053*conj(IT_0053) + IT_0058*conj(IT_0058))*IT_0067 + (conj(IT_0048)
      *IT_0051 + IT_0048*conj(IT_0051) + conj(IT_0053)*IT_0058 + IT_0053*conj
      (IT_0058))*IT_0068 + IT_0069*(conj(IT_0036)*IT_0070 + IT_0036*IT_0071) + 6
      *(conj(IT_0021)*IT_0048 + IT_0021*conj(IT_0048) + conj(IT_0036)*IT_0065 +
       IT_0036*IT_0066 + 0.166666666666667*conj(IT_0051)*IT_0070 +
       0.166666666666667*IT_0051*IT_0071)*IT_0072;
    return create_ccomplex_return(IT_0073);
}

