#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_sc_L_to_nu_tau_anti_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_sc_L_to_nu_tau_anti_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B1)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0019 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B2)*e_em;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W2)*e_em;
    const ccomplex_t IT_0026 = IT_0006*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0034 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0032*IT_0033*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*e_em;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W3)*e_em;
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0048 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0046*IT_0047*IT_0048;
    const ccomplex_t IT_0050 = conj(N_B4)*e_em;
    const ccomplex_t IT_0051 = IT_0001*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W4)*e_em;
    const ccomplex_t IT_0054 = IT_0006*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + -IT_0055);
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0062 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0060*IT_0061*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0021 + 0.5*IT_0035 + 0.5*IT_0049 + 0.5
      *IT_0063;
    const ccomplex_t IT_0065 = N_B2*e_em;
    const ccomplex_t IT_0066 = IT_0001*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = N_W2*e_em;
    const ccomplex_t IT_0069 = IT_0006*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + 3*IT_0070);
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = IT_0029*IT_0072;
    const ccomplex_t IT_0074 = IT_0034*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = N_B3*e_em;
    const ccomplex_t IT_0077 = IT_0001*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = N_W3*e_em;
    const ccomplex_t IT_0080 = IT_0006*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + 3*IT_0081);
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0043*IT_0083;
    const ccomplex_t IT_0085 = IT_0048*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = N_B4*e_em;
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = N_W4*e_em;
    const ccomplex_t IT_0091 = IT_0006*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + 3*IT_0092);
    const ccomplex_t IT_0094 = 0.166666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0057*IT_0094;
    const ccomplex_t IT_0096 = IT_0062*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = N_B1*e_em;
    const ccomplex_t IT_0099 = IT_0001*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = N_W1*e_em;
    const ccomplex_t IT_0102 = IT_0006*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0011*IT_0105;
    const ccomplex_t IT_0107 = IT_0020*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0075 + 0.5*IT_0086 + 0.5*IT_0097 + 0.5
      *IT_0108;
    const ccomplex_t IT_0110 = m_c*s_13;
    const ccomplex_t IT_0111 = s_13*s_14;
    const ccomplex_t IT_0112 = s_34*IT_0019;
    return create_ccomplex_return(6*IT_0064*(s_34*conj(IT_0064) + -conj
      (IT_0109)*IT_0110) + (-6)*IT_0109*(conj(IT_0064)*IT_0110 + (
      -0.166666666666667)*conj(IT_0109)*(12*IT_0111 + (-6)*IT_0112)));
}

