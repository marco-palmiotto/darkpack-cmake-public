#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_sc_R_to_nu_mu_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_sc_R_to_nu_mu_c(
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
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0019 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B1)*e_em;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W1)*e_em;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + -IT_0028);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0035 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B3)*e_em;
    const ccomplex_t IT_0038 = IT_0001*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W3)*e_em;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0049 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0048*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B4)*e_em;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W4)*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + -IT_0056);
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0063 = cpow(s_13 + (-0.5)*IT_0019 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0061*IT_0062*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0036 + 0.5*IT_0050 + 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0067 = 0.666666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0030*IT_0067;
    const ccomplex_t IT_0069 = IT_0035*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = IT_0044*IT_0072;
    const ccomplex_t IT_0074 = IT_0049*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0077 = 0.666666666666667*IT_0076;
    const ccomplex_t IT_0078 = IT_0058*IT_0077;
    const ccomplex_t IT_0079 = IT_0063*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0082 = 0.666666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = IT_0020*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0070 + 0.5*IT_0075 + 0.5*IT_0080 + 0.5
      *IT_0085;
    const ccomplex_t IT_0087 = m_c*s_13;
    const ccomplex_t IT_0088 = s_13*s_14;
    const ccomplex_t IT_0089 = s_34*IT_0019;
    const ccomplex_t IT_0090 = (-6)*s_34;
    const ccomplex_t IT_0091 = 6*s_34*(IT_0022*conj(IT_0022) + IT_0065*conj
      (IT_0065)) + 6*((conj(IT_0022) + -conj(IT_0065))*IT_0086 + (IT_0022 + 
      -IT_0065)*conj(IT_0086))*IT_0087 + IT_0086*conj(IT_0086)*(12*IT_0088 + (-6
      )*IT_0089) + (conj(IT_0022)*IT_0065 + IT_0022*conj(IT_0065))*IT_0090;
    return create_ccomplex_return(IT_0091);
}

