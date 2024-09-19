#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_sb_2_to_tau_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_sb_2_to_tau_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0001*IT_0003*IT_0005*U_sb_01;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_tau, 2);
    const ccomplex_t IT_0010 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0011 = cpow(s_13 + (-0.5)*IT_0009 + (-0.5)*IT_0010 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cpow(s_13 + (-0.5)*IT_0009 + (-0.5)*IT_0010 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0003*IT_0005*U_sb_01;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0013 + 0.5*IT_0020;
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0001*IT_0003*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0023;
    const ccomplex_t IT_0027 = m_b*conj(U_d2)*V_tb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0029 + (-2)*IT_0031);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = m_b*conj(U_d1)*V_tb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0038 = IT_0026*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0039 + (-2)*IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0001*IT_0003*IT_0023;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0036 + 0.5*IT_0048;
    const ccomplex_t IT_0050 = m_t*m_tau*IT_0010;
    const ccomplex_t IT_0051 = s_34*IT_0010;
    const ccomplex_t IT_0052 = (-6)*IT_0051;
    const ccomplex_t IT_0053 = s_13*s_14;
    const ccomplex_t IT_0054 = 12*IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = m_t*s_13;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0058 = IT_0017*IT_0045;
    const ccomplex_t IT_0059 = IT_0014*IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0061 = IT_0007*IT_0025;
    const ccomplex_t IT_0062 = IT_0011*IT_0060*IT_0061;
    const ccomplex_t IT_0063 = m_tau*IT_0020;
    const ccomplex_t IT_0064 = m_tau*IT_0013;
    const ccomplex_t IT_0065 = (-0.5)*IT_0059 + (-0.5)*IT_0062 + 0.5*IT_0063 +
       0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0015*IT_0043;
    const ccomplex_t IT_0067 = IT_0014*IT_0057*IT_0066;
    const ccomplex_t IT_0068 = m_tau*IT_0048;
    const ccomplex_t IT_0069 = m_tau*IT_0036;
    const ccomplex_t IT_0070 = (-0.5)*IT_0067 + 0.5*IT_0068 + 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0002*IT_0033;
    const ccomplex_t IT_0072 = IT_0011*IT_0060*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0070 + IT_0074;
    const ccomplex_t IT_0076 = conj(IT_0070) + -conj(IT_0073);
    const ccomplex_t IT_0077 = 6*s_34;
    const ccomplex_t IT_0078 = (-6)*s_34;
    const ccomplex_t IT_0079 = s_14*m_tau;
    const ccomplex_t IT_0080 = (-6)*IT_0065;
    const ccomplex_t IT_0081 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0082 = m_t*m_tau;
    const ccomplex_t IT_0083 = (-6)*conj(IT_0021)*(IT_0049*IT_0050 + (
      -0.166666666666667)*IT_0021*IT_0055) + (-6)*conj(IT_0049)*(IT_0021*IT_0050
       + (-0.166666666666667)*IT_0049*IT_0055) + 6*IT_0056*(conj(IT_0049)
      *IT_0065 + IT_0049*conj(IT_0065) + conj(IT_0021)*IT_0075 + IT_0021*IT_0076
      ) + (IT_0065*conj(IT_0065) + IT_0070*conj(IT_0070) + IT_0073*conj(IT_0073)
      )*IT_0077 + (conj(IT_0070)*IT_0073 + IT_0070*conj(IT_0073))*IT_0078 + (-6)
      *IT_0079*(conj(IT_0049)*IT_0075 + IT_0049*IT_0076 + (-0.166666666666667)
      *conj(IT_0021)*IT_0080 + (-0.166666666666667)*IT_0021*IT_0081) + 6*(conj
      (IT_0065)*IT_0073 + IT_0065*conj(IT_0073) + 0.166666666666667*conj(IT_0070
      )*IT_0080 + 0.166666666666667*IT_0070*IT_0081)*IT_0082;
    return create_ccomplex_return(IT_0083);
}

