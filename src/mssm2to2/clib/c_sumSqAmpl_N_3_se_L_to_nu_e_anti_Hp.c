#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_se_L_to_nu_e_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_se_L_to_nu_e_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
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
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = conj(N_u3)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = conj(N_B3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0002*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0023 = pow(m_se_L, 2);
    const ccomplex_t IT_0024 = cpow(s_23 + (-0.5)*IT_0023 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0022*IT_0024;
    const ccomplex_t IT_0026 = pow(m_e, 2);
    const ccomplex_t IT_0027 = pow(m_N_3, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + -IT_0023 + IT_0026 + -IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0030 = pow(m_W, -1);
    const ccomplex_t IT_0031 = cpow(IT_0003, -1);
    const ccomplex_t IT_0032 = sin(beta);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0001*IT_0030*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B3)*e_em;
    const ccomplex_t IT_0036 = IT_0009*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W3)*e_em;
    const ccomplex_t IT_0039 = IT_0001*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = IT_0028*IT_0029*IT_0043;
    const ccomplex_t IT_0045 = cpow(s_13 + (-0.5)*IT_0027 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001*
      (m_W*IT_0003 + (-0.5)*IT_0026*IT_0030*IT_0031)*IT_0032;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0037 + -IT_0040);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0001*IT_0030*IT_0031;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0034*IT_0054;
    const ccomplex_t IT_0056 = IT_0028*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_N_3*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0060 = conj(N_u3)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0061 = IT_0004*IT_0060;
    const ccomplex_t IT_0062 = 2*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0064 = IT_0010*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = IT_0065 + IT_0068;
    const ccomplex_t IT_0070 = IT_0062 + IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0059*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0075 = cpow(s_23 + (-0.5)*IT_0023 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0073*IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0025 + IT_0044 + 0.5*IT_0052 + 
      -IT_0058 + (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0001*IT_0032;
    const ccomplex_t IT_0079 = N_d3*e_em*U_Wm1;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0009*IT_0032;
    const ccomplex_t IT_0083 = N_B3*U_d1*e_em;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W3*U_d1*e_em;
    const ccomplex_t IT_0087 = IT_0078*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = IT_0081 + IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0059*IT_0092;
    const ccomplex_t IT_0094 = IT_0075*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = N_W3*U_d2*e_em;
    const ccomplex_t IT_0097 = IT_0078*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = N_B3*U_d2*e_em;
    const ccomplex_t IT_0100 = IT_0082*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = IT_0098 + IT_0101;
    const ccomplex_t IT_0103 = N_d3*e_em*U_Wm2;
    const ccomplex_t IT_0104 = IT_0078*IT_0103;
    const ccomplex_t IT_0105 = (-2)*IT_0104;
    const ccomplex_t IT_0106 = IT_0102 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0002*IT_0108;
    const ccomplex_t IT_0110 = IT_0024*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0057 + (-0.5)*IT_0095 + (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = s_23*m_N_3;
    const ccomplex_t IT_0114 = s_12*s_23;
    const ccomplex_t IT_0115 = s_13*IT_0023;
    return create_ccomplex_return(2*IT_0077*(s_13*conj(IT_0077) + conj(IT_0112
      )*IT_0113) + 2*IT_0112*(conj(IT_0077)*IT_0113 + 2*conj(IT_0112)*(IT_0114 +
       (-0.5)*IT_0115)));
}

