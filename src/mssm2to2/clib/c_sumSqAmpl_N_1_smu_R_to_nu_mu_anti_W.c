#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_R_to_nu_mu_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_R_to_nu_mu_anti_W(
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
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0000 = pow(m_smu_R, 2);
    const ccomplex_t IT_0001 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0010 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0025 = IT_0006*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0018 + IT_0031;
    const ccomplex_t IT_0033 = m_N_1*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0036 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0008*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0035*IT_0043;
    const ccomplex_t IT_0045 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0046 = IT_0006*IT_0045;
    const ccomplex_t IT_0047 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0048 = IT_0006*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-0.5)*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0021*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0054 = IT_0019*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = pow(m_mu, 2);
    const ccomplex_t IT_0056 = pow(m_N_1, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + -IT_0000 + IT_0055 + -IT_0056 
      + -reg_prop, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0006;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cos(theta_W);
    const ccomplex_t IT_0062 = cpow(IT_0061, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0060*IT_0064;
    const ccomplex_t IT_0066 = IT_0057*IT_0058*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0044 + 0.5*IT_0054 + -IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = IT_0060*IT_0069;
    const ccomplex_t IT_0071 = IT_0057*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = 2*IT_0072;
    const ccomplex_t IT_0074 = pow(m_W, -2);
    const ccomplex_t IT_0075 = pow(s_34, 2);
    const ccomplex_t IT_0076 = m_N_1*IT_0074*IT_0075;
    const ccomplex_t IT_0077 = pow(m_W, 2);
    const ccomplex_t IT_0078 = s_13*IT_0077;
    const ccomplex_t IT_0079 = s_14*s_34;
    const ccomplex_t IT_0080 = (-2)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = IT_0074*IT_0081;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = s_13 + IT_0083;
    const ccomplex_t IT_0085 = s_14*s_34*IT_0074;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = s_13 + IT_0086;
    const ccomplex_t IT_0088 = m_N_1*IT_0087;
    const ccomplex_t IT_0089 = -IT_0018 + -IT_0031;
    const ccomplex_t IT_0090 = 2*conj(IT_0089);
    const ccomplex_t IT_0091 = 0.5*IT_0088*IT_0090;
    const ccomplex_t IT_0092 = s_34*IT_0077;
    const ccomplex_t IT_0093 = IT_0074*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = s_34 + IT_0094;
    const ccomplex_t IT_0096 = m_N_1*IT_0095;
    const ccomplex_t IT_0097 = 2*IT_0096;
    const ccomplex_t IT_0098 = conj(IT_0072)*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0018 + 0.5*IT_0031;
    const ccomplex_t IT_0100 = (-0.25)*IT_0093;
    const ccomplex_t IT_0101 = s_34 + IT_0100;
    const ccomplex_t IT_0102 = m_N_1*IT_0101;
    const ccomplex_t IT_0103 = (-4)*IT_0102;
    const ccomplex_t IT_0104 = -conj(IT_0099)*IT_0103;
    const ccomplex_t IT_0105 = conj(IT_0073)*IT_0076 + 2*(conj(IT_0034) + conj
      (IT_0067))*IT_0084 + IT_0091 + IT_0098 + IT_0104;
    const ccomplex_t IT_0106 = pow(s_14, 2);
    const ccomplex_t IT_0107 = IT_0074*IT_0106;
    const ccomplex_t IT_0108 = -IT_0056;
    const ccomplex_t IT_0109 = IT_0107 + IT_0108;
    const ccomplex_t IT_0110 = s_13*IT_0109;
    const ccomplex_t IT_0111 = s_13*IT_0087;
    const ccomplex_t IT_0112 = conj(IT_0034) + conj(IT_0067);
    const ccomplex_t IT_0113 = s_13*s_14;
    const ccomplex_t IT_0114 = s_34*IT_0056;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = conj(IT_0072) + conj(IT_0099);
    const ccomplex_t IT_0118 = 2*IT_0089;
    const ccomplex_t IT_0119 = s_13*IT_0074*IT_0075;
    const ccomplex_t IT_0120 = s_13*s_34;
    const ccomplex_t IT_0121 = pow(m_W, 4);
    const ccomplex_t IT_0122 = s_13*IT_0121;
    const ccomplex_t IT_0123 = IT_0074*IT_0122;
    const ccomplex_t IT_0124 = IT_0078 + IT_0080 + (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = (-0.25)*IT_0123;
    const ccomplex_t IT_0126 = IT_0079 + IT_0125;
    const ccomplex_t IT_0127 = conj(IT_0073)*IT_0120;
    const ccomplex_t IT_0128 = 0.5*IT_0090*IT_0116;
    const ccomplex_t IT_0129 = 2*(IT_0034 + IT_0067)*IT_0105 + (conj(IT_0089)
      *IT_0110 + conj(IT_0073)*IT_0111 + IT_0088*IT_0112 + IT_0116*IT_0117)
      *IT_0118 + 2*IT_0073*(0.5*IT_0090*IT_0111 + IT_0076*IT_0112 + conj(IT_0073
      )*IT_0119 + IT_0117*IT_0120) + (-4)*IT_0099*(0.5*IT_0103*IT_0112 + conj
      (IT_0099)*IT_0124 + (-2)*conj(IT_0072)*IT_0126 + (-0.5)*IT_0127 + (-0.5)
      *IT_0128) + (-4)*IT_0072*((-0.5)*IT_0097*IT_0112 + conj(IT_0072)*IT_0124 +
       (-2)*conj(IT_0099)*IT_0126 + (-0.5)*IT_0127 + (-0.5)*IT_0128);
    return create_ccomplex_return(IT_0129);
}

