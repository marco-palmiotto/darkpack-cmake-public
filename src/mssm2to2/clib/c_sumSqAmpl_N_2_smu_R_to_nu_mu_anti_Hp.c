#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_smu_R_to_nu_mu_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_smu_R_to_nu_mu_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = conj(N_W2)*e_em;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0009 + -IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = pow(m_N_2, 2);
    const ccomplex_t IT_0017 = cpow(s_13 + (-0.5)*IT_0016 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0002*IT_0020;
    const ccomplex_t IT_0025 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = 2*IT_0026;
    const ccomplex_t IT_0028 = IT_0006*IT_0020;
    const ccomplex_t IT_0029 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0033 = IT_0024*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = IT_0031 + IT_0034;
    const ccomplex_t IT_0036 = IT_0027 + IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0023*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0041 = pow(m_smu_R, 2);
    const ccomplex_t IT_0042 = cpow(s_23 + (-0.5)*IT_0041 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0039*IT_0040*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0021;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0047 = IT_0024*IT_0046;
    const ccomplex_t IT_0048 = 2*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0050 = IT_0028*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0053 = IT_0024*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = IT_0051 + IT_0054;
    const ccomplex_t IT_0056 = IT_0048 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0045*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0061 = cpow(s_23 + (-0.5)*IT_0041 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0059*IT_0060*IT_0061;
    const ccomplex_t IT_0063 = sin(beta);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0002*IT_0021*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0006;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = pow(m_mu, 2);
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + -IT_0016 + -IT_0041 + IT_0069 
      + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = m_N_2*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0021;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0065*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0078 = IT_0070*IT_0076*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0019 + (-0.5)*IT_0043 + (-0.5)*IT_0062 +
       -IT_0073 + IT_0078;
    const ccomplex_t IT_0080 = IT_0002*IT_0063;
    const ccomplex_t IT_0081 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (-2)*IT_0082;
    const ccomplex_t IT_0084 = IT_0006*IT_0063;
    const ccomplex_t IT_0085 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0089 = IT_0080*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = IT_0087 + IT_0090;
    const ccomplex_t IT_0092 = IT_0083 + IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0045*IT_0094;
    const ccomplex_t IT_0096 = IT_0061*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0099 = IT_0080*IT_0098;
    const ccomplex_t IT_0100 = (-2)*IT_0099;
    const ccomplex_t IT_0101 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0102 = IT_0084*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0105 = IT_0080*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = IT_0103 + IT_0106;
    const ccomplex_t IT_0108 = IT_0100 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0023*IT_0110;
    const ccomplex_t IT_0112 = IT_0042*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = -IT_0072 + (-0.5)*IT_0097 + (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = s_23*m_N_2;
    const ccomplex_t IT_0116 = s_12*s_23;
    const ccomplex_t IT_0117 = s_13*IT_0041;
    return create_ccomplex_return(2*IT_0079*(s_13*conj(IT_0079) + conj(IT_0114
      )*IT_0115) + 2*IT_0114*(conj(IT_0079)*IT_0115 + 2*conj(IT_0114)*(IT_0116 +
       (-0.5)*IT_0117)));
}

