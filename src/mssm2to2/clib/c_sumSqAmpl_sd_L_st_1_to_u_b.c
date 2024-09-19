#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_st_1_to_u_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_st_1_to_u_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0001*IT_0003*IT_0005*U_st_00;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_u, 2);
    const ccomplex_t IT_0010 = pow(m_sd_L, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0001*IT_0003*IT_0005*U_st_00;
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
    const ccomplex_t IT_0026 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0023;
    const ccomplex_t IT_0029 = m_t*V_tb*conj(V_u1)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0027 + (-0.5)*IT_0031);
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = IT_0018*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0003*IT_0023;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = m_t*V_tb*conj(V_u2)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0042 = IT_0028*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + (-0.5)*IT_0043);
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0038*IT_0045;
    const ccomplex_t IT_0047 = IT_0011*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = -IT_0036 + -IT_0048;
    const ccomplex_t IT_0050 = pow(m_st_1, 2);
    const ccomplex_t IT_0051 = s_34*IT_0050;
    const ccomplex_t IT_0052 = (-18)*IT_0051;
    const ccomplex_t IT_0053 = s_23*s_24;
    const ccomplex_t IT_0054 = 36*IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = m_b*s_23;
    const ccomplex_t IT_0057 = IT_0002*IT_0045;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0059 = IT_0011*IT_0057*IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = IT_0014*IT_0033;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0063 = IT_0018*IT_0061*IT_0062;
    const ccomplex_t IT_0064 = m_b*IT_0020;
    const ccomplex_t IT_0065 = m_b*IT_0013;
    const ccomplex_t IT_0066 = IT_0063 + IT_0064 + IT_0065;
    const ccomplex_t IT_0067 = IT_0016*IT_0025;
    const ccomplex_t IT_0068 = IT_0018*IT_0062*IT_0067;
    const ccomplex_t IT_0069 = m_b*IT_0036;
    const ccomplex_t IT_0070 = m_b*IT_0048;
    const ccomplex_t IT_0071 = -IT_0068 + -IT_0069 + -IT_0070;
    const ccomplex_t IT_0072 = (-18)*IT_0071;
    const ccomplex_t IT_0073 = IT_0007*IT_0038;
    const ccomplex_t IT_0074 = IT_0011*IT_0058*IT_0073;
    const ccomplex_t IT_0075 = 18*IT_0074;
    const ccomplex_t IT_0076 = IT_0072 + IT_0075;
    const ccomplex_t IT_0077 = (-18)*conj(IT_0071);
    const ccomplex_t IT_0078 = 18*conj(IT_0074);
    const ccomplex_t IT_0079 = IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = 18*s_34;
    const ccomplex_t IT_0081 = (-18)*s_34;
    const ccomplex_t IT_0082 = m_b*m_u;
    const ccomplex_t IT_0083 = (-18)*IT_0060;
    const ccomplex_t IT_0084 = (-18)*conj(IT_0060);
    const ccomplex_t IT_0085 = m_b*m_u*IT_0050;
    const ccomplex_t IT_0086 = (-18)*IT_0021;
    const ccomplex_t IT_0087 = (-18)*conj(IT_0021);
    const ccomplex_t IT_0088 = m_u*s_24;
    const ccomplex_t IT_0089 = (IT_0021*conj(IT_0021) + IT_0049*conj(IT_0049))
      *IT_0055 + IT_0056*(conj(IT_0021)*(18*IT_0060 + (-18)*IT_0066) + IT_0021*
      (18*conj(IT_0060) + (-18)*conj(IT_0066)) + conj(IT_0049)*IT_0076 + IT_0049
      *IT_0079) + (IT_0060*conj(IT_0060) + IT_0066*conj(IT_0066) + IT_0071*conj
      (IT_0071) + IT_0074*conj(IT_0074))*IT_0080 + (conj(IT_0060)*IT_0066 +
       IT_0060*conj(IT_0066) + conj(IT_0071)*IT_0074 + IT_0071*conj(IT_0074))
      *IT_0081 + IT_0082*(18*conj(IT_0060)*IT_0071 + 18*IT_0060*conj(IT_0071) +
       conj(IT_0066)*IT_0076 + IT_0066*IT_0079 + conj(IT_0074)*IT_0083 + IT_0074
      *IT_0084) + IT_0085*(conj(IT_0049)*IT_0086 + IT_0049*IT_0087) + (18*conj
      (IT_0021)*IT_0071 + 18*IT_0021*conj(IT_0071) + conj(IT_0049)*(18*IT_0066 +
       IT_0083) + IT_0049*(18*conj(IT_0066) + IT_0084) + conj(IT_0074)*IT_0086 +
       IT_0074*IT_0087)*IT_0088;
    return create_ccomplex_return(IT_0089);
}

