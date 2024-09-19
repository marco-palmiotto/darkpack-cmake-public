#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_smu_R_to_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_smu_R_to_mu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_mu, 2);
    const ccomplex_t IT_0013 = pow(m_smu_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0011*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0022*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0046*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0016 + -IT_0024 + -IT_0032 + IT_0035 + 
      -IT_0038 + IT_0041 + -IT_0049 + IT_0052;
    const ccomplex_t IT_0054 = IT_0012*IT_0013;
    const ccomplex_t IT_0055 = s_23*s_24;
    const ccomplex_t IT_0056 = s_34*IT_0013;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055 + IT_0057;
    const ccomplex_t IT_0059 = m_mu*s_23;
    const ccomplex_t IT_0060 = cpow(IT_0026, 2);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_N_1*IT_0060;
    const ccomplex_t IT_0062 = IT_0033*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0006, 2);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_N_2*IT_0063;
    const ccomplex_t IT_0065 = IT_0014*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0043, 2);
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_N_4*IT_0066;
    const ccomplex_t IT_0068 = IT_0050*IT_0067;
    const ccomplex_t IT_0069 = -IT_0062 + -IT_0065 + -IT_0068;
    const ccomplex_t IT_0070 = IT_0030*IT_0061;
    const ccomplex_t IT_0071 = m_mu*IT_0024;
    const ccomplex_t IT_0072 = m_mu*IT_0038;
    const ccomplex_t IT_0073 = m_mu*IT_0041;
    const ccomplex_t IT_0074 = m_mu*IT_0049;
    const ccomplex_t IT_0075 = m_mu*IT_0032;
    const ccomplex_t IT_0076 = IT_0036*IT_0064;
    const ccomplex_t IT_0077 = cpow(IT_0021, 2);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_N_3*IT_0077;
    const ccomplex_t IT_0079 = IT_0017*IT_0078;
    const ccomplex_t IT_0080 = IT_0039*IT_0078;
    const ccomplex_t IT_0081 = IT_0047*IT_0067;
    const ccomplex_t IT_0082 = m_mu*IT_0035;
    const ccomplex_t IT_0083 = m_mu*IT_0016;
    const ccomplex_t IT_0084 = m_mu*IT_0052;
    const ccomplex_t IT_0085 = IT_0070 + -IT_0071 + -IT_0072 + -IT_0073 + 
      -IT_0074 + -IT_0075 + IT_0076 + IT_0079 + IT_0080 + IT_0081 + -IT_0082 + 
      -IT_0083 + -IT_0084;
    const ccomplex_t IT_0086 = IT_0069 + -IT_0085;
    const ccomplex_t IT_0087 = -conj(IT_0085);
    const ccomplex_t IT_0088 = conj(IT_0069) + IT_0087;
    const ccomplex_t IT_0089 = cpow(IT_0028, 2);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_1*IT_0089;
    const ccomplex_t IT_0091 = IT_0030*IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0045, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_N_4*IT_0092;
    const ccomplex_t IT_0094 = IT_0047*IT_0093;
    const ccomplex_t IT_0095 = cpow(IT_0010, 2);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_N_2*IT_0095;
    const ccomplex_t IT_0097 = IT_0014*IT_0096;
    const ccomplex_t IT_0098 = -IT_0091 + -IT_0094 + -IT_0097;
    const ccomplex_t IT_0099 = 2*IT_0098;
    const ccomplex_t IT_0100 = IT_0033*IT_0090;
    const ccomplex_t IT_0101 = IT_0036*IT_0096;
    const ccomplex_t IT_0102 = cpow(IT_0019, 2);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_N_3*IT_0102;
    const ccomplex_t IT_0104 = IT_0017*IT_0103;
    const ccomplex_t IT_0105 = IT_0039*IT_0103;
    const ccomplex_t IT_0106 = IT_0050*IT_0093;
    const ccomplex_t IT_0107 = -IT_0071 + -IT_0072 + -IT_0073 + -IT_0074 + 
      -IT_0075 + -IT_0082 + -IT_0083 + -IT_0084 + IT_0100 + IT_0101 + IT_0104 +
       IT_0105 + IT_0106;
    const ccomplex_t IT_0108 = (-2)*IT_0107;
    const ccomplex_t IT_0109 = IT_0099 + IT_0108;
    const ccomplex_t IT_0110 = 2*conj(IT_0098);
    const ccomplex_t IT_0111 = (-2)*conj(IT_0107);
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0107;
    const ccomplex_t IT_0114 = (-2)*IT_0098;
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = 2*conj(IT_0107);
    const ccomplex_t IT_0117 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = m_mu*s_24;
    const ccomplex_t IT_0120 = 2*s_34;
    const ccomplex_t IT_0121 = (-2)*s_34;
    const ccomplex_t IT_0122 = (-4)*conj(IT_0053)*(IT_0053*IT_0054 + (-2)
      *IT_0053*IT_0058) + 2*IT_0059*(conj(IT_0053)*IT_0086 + IT_0053*IT_0088 +
       0.5*conj(IT_0053)*IT_0109 + 0.5*IT_0053*IT_0112) + IT_0012*(conj(IT_0085)
      *IT_0109 + IT_0085*IT_0112 + conj(IT_0069)*IT_0115 + IT_0069*IT_0118) + (
      -2)*(conj(IT_0053)*IT_0086 + IT_0053*IT_0088 + (-0.5)*conj(IT_0053)
      *IT_0115 + (-0.5)*IT_0053*IT_0118)*IT_0119 + (IT_0069*conj(IT_0069) +
       IT_0085*conj(IT_0085) + IT_0098*conj(IT_0098) + IT_0107*conj(IT_0107))
      *IT_0120 + (conj(IT_0069)*IT_0085 + IT_0069*conj(IT_0085) + conj(IT_0098)
      *IT_0107 + IT_0098*conj(IT_0107))*IT_0121;
    return create_ccomplex_return(IT_0122);
}

