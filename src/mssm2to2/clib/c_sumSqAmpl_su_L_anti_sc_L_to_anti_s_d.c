#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sc_L_to_anti_s_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sc_L_to_anti_s_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_s, 2);
    const ccomplex_t IT_0012 = pow(m_sc_L, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0015*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0024 = IT_0009*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = m_d*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0029 = IT_0019*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_d*IT_0031;
    const ccomplex_t IT_0033 = IT_0014 + IT_0022 + -IT_0027 + -IT_0032;
    const ccomplex_t IT_0034 = IT_0026 + IT_0031;
    const ccomplex_t IT_0035 = m_d*s_13;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0017*IT_0037;
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = IT_0013*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0040 + IT_0045;
    const ccomplex_t IT_0047 = m_s*s_14;
    const ccomplex_t IT_0048 = 18*s_34;
    const ccomplex_t IT_0049 = m_d*m_s;
    const ccomplex_t IT_0050 = IT_0023*IT_0042;
    const ccomplex_t IT_0051 = IT_0000*IT_0013*IT_0050;
    const ccomplex_t IT_0052 = IT_0028*IT_0037;
    const ccomplex_t IT_0053 = IT_0015*IT_0021*IT_0052;
    const ccomplex_t IT_0054 = m_d*IT_0045;
    const ccomplex_t IT_0055 = m_d*IT_0040;
    const ccomplex_t IT_0056 = IT_0051 + IT_0053 + -IT_0054 + -IT_0055;
    const ccomplex_t IT_0057 = (-18)*conj(IT_0056);
    const ccomplex_t IT_0058 = pow(m_su_L, 2);
    const ccomplex_t IT_0059 = m_d*m_s*IT_0058;
    const ccomplex_t IT_0060 = s_34*IT_0058;
    const ccomplex_t IT_0061 = (-18)*IT_0060;
    const ccomplex_t IT_0062 = s_13*s_14;
    const ccomplex_t IT_0063 = 36*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = (-18)*conj(IT_0033);
    const ccomplex_t IT_0066 = (-18)*conj(IT_0034);
    const ccomplex_t IT_0067 = IT_0033*(18*conj(IT_0034)*IT_0035 + (-18)*conj
      (IT_0046)*IT_0047 + conj(IT_0033)*IT_0048 + IT_0049*IT_0057) + IT_0034*(18
      *conj(IT_0033)*IT_0035 + IT_0047*IT_0057 + (-18)*conj(IT_0046)*IT_0059 +
       conj(IT_0034)*IT_0064) + IT_0056*(18*IT_0035*conj(IT_0046) + IT_0048*conj
      (IT_0056) + IT_0049*IT_0065 + IT_0047*IT_0066) + IT_0046*(18*IT_0035*conj
      (IT_0056) + conj(IT_0046)*IT_0064 + IT_0047*IT_0065 + IT_0059*IT_0066);
    return create_ccomplex_return(IT_0067);
}

