#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_L_to_nu_mu_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_L_to_nu_mu_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = pow(m_mu, 2);
    const ccomplex_t IT_0003 = cpow(IT_0000, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*e_em*(m_W
      *IT_0000 + (-0.5)*IT_0001*IT_0002*IT_0003)*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = conj(N_B1)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W1)*e_em;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + -IT_0016);
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = pow(m_N_1, 2);
    const ccomplex_t IT_0021 = cpow(s_13 + (-0.5)*IT_0020 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0025 = IT_0000*IT_0006;
    const ccomplex_t IT_0026 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 2*IT_0027;
    const ccomplex_t IT_0029 = IT_0000*IT_0010;
    const ccomplex_t IT_0030 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = IT_0032 + IT_0035;
    const ccomplex_t IT_0037 = IT_0028 + IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0024*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0042 = pow(m_smu_L, 2);
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0042 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0040*IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0046 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0047 = IT_0025*IT_0046;
    const ccomplex_t IT_0048 = 2*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0050 = IT_0029*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0053 = IT_0025*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = IT_0051 + IT_0054;
    const ccomplex_t IT_0056 = IT_0048 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0045*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0061 = cpow(s_23 + (-0.5)*IT_0042 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0059*IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0001*IT_0003*IT_0004*IT_0006;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + IT_0002 + -IT_0020 + -IT_0042 
      + -reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = m_N_1*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0013 + IT_0016);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0064*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0076 = IT_0068*IT_0074*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0023 + (-0.5)*IT_0044 + (-0.5)*IT_0062 +
       -IT_0071 + IT_0076;
    const ccomplex_t IT_0078 = IT_0004*IT_0006;
    const ccomplex_t IT_0079 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0004*IT_0010;
    const ccomplex_t IT_0083 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0087 = IT_0078*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = IT_0081 + IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0045*IT_0092;
    const ccomplex_t IT_0094 = IT_0061*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0097 = IT_0078*IT_0096;
    const ccomplex_t IT_0098 = (-2)*IT_0097;
    const ccomplex_t IT_0099 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0100 = IT_0082*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0103 = IT_0078*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = IT_0101 + IT_0104;
    const ccomplex_t IT_0106 = IT_0098 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0024*IT_0108;
    const ccomplex_t IT_0110 = IT_0043*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0070 + (-0.5)*IT_0095 + (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = s_23*m_N_1;
    const ccomplex_t IT_0114 = s_12*s_23;
    const ccomplex_t IT_0115 = s_13*IT_0042;
    return create_ccomplex_return(2*IT_0077*(s_13*conj(IT_0077) + conj(IT_0112
      )*IT_0113) + 2*IT_0112*(conj(IT_0077)*IT_0113 + 2*conj(IT_0112)*(IT_0114 +
       (-0.5)*IT_0115)));
}

