#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_ss_R_to_nu_tau_anti_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_ss_R_to_nu_tau_anti_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0001;
    const ccomplex_t IT_0003 = (-0.333333333333333)*IT_0002;
    const ccomplex_t IT_0004 = conj(N_B1)*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_W1)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0016 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0017 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0001;
    const ccomplex_t IT_0020 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0021 = conj(N_B2)*e_em;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W2)*e_em;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + -IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0020*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0001;
    const ccomplex_t IT_0034 = (-0.333333333333333)*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B3)*e_em;
    const ccomplex_t IT_0036 = IT_0001*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W3)*e_em;
    const ccomplex_t IT_0039 = IT_0008*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0045 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0001;
    const ccomplex_t IT_0048 = (-0.333333333333333)*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B4)*e_em;
    const ccomplex_t IT_0050 = IT_0001*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W4)*e_em;
    const ccomplex_t IT_0053 = IT_0008*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + -IT_0054);
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0048*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0059 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0018 + 0.5*IT_0032 + 0.5*IT_0046 + 0.5
      *IT_0060;
    const ccomplex_t IT_0062 = pow(m_W, -1);
    const ccomplex_t IT_0063 = cos(beta);
    const ccomplex_t IT_0064 = cpow(IT_0063, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0008*IT_0062*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0028*IT_0066;
    const ccomplex_t IT_0068 = IT_0031*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0008*IT_0062*IT_0064;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0042*IT_0071;
    const ccomplex_t IT_0073 = IT_0045*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0008*IT_0062*IT_0064;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0056*IT_0076;
    const ccomplex_t IT_0078 = IT_0059*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0008*IT_0062*IT_0064;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0013*IT_0081;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0069 + 0.5*IT_0074 + 0.5*IT_0079 + 0.5
      *IT_0084;
    const ccomplex_t IT_0086 = m_s*s_13;
    const ccomplex_t IT_0087 = s_13*s_14;
    const ccomplex_t IT_0088 = s_34*IT_0016;
    return create_ccomplex_return(6*IT_0061*(s_34*conj(IT_0061) + -conj
      (IT_0085)*IT_0086) + (-6)*IT_0085*(conj(IT_0061)*IT_0086 + (
      -0.166666666666667)*conj(IT_0085)*(12*IT_0087 + (-6)*IT_0088)));
}

