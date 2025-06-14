#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_smu_L_to_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_smu_L_to_mu_mu(
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
    const creal_t m_smu_L = param->m_smu_L;
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
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_N_1*IT_0007;
    const ccomplex_t IT_0009 = pow(m_mu, 2);
    const ccomplex_t IT_0010 = pow(m_smu_L, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_N_2*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_N_3*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = cpow(IT_0032, 2);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_N_4*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0034*IT_0037;
    const ccomplex_t IT_0039 = -IT_0012 + -IT_0014 + -IT_0020 + -IT_0022 + 
      -IT_0028 + -IT_0030 + -IT_0036 + -IT_0038;
    const ccomplex_t IT_0040 = cos(theta_W);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = conj(N_B1)*e_em;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W1)*e_em;
    const ccomplex_t IT_0046 = IT_0004*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0049, 2);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_N_1*IT_0050;
    const ccomplex_t IT_0052 = IT_0013*IT_0051;
    const ccomplex_t IT_0053 = conj(N_B2)*e_em;
    const ccomplex_t IT_0054 = IT_0041*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W2)*e_em;
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = cpow(IT_0060, 2);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_N_2*IT_0061;
    const ccomplex_t IT_0063 = IT_0021*IT_0062;
    const ccomplex_t IT_0064 = IT_0019*IT_0062;
    const ccomplex_t IT_0065 = conj(N_B3)*e_em;
    const ccomplex_t IT_0066 = IT_0041*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W3)*e_em;
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = cpow(IT_0072, 2);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*m_N_3*IT_0073;
    const ccomplex_t IT_0075 = IT_0027*IT_0074;
    const ccomplex_t IT_0076 = IT_0029*IT_0074;
    const ccomplex_t IT_0077 = conj(N_B4)*e_em;
    const ccomplex_t IT_0078 = IT_0041*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W4)*e_em;
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + IT_0082);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = cpow(IT_0084, 2);
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_N_4*IT_0085;
    const ccomplex_t IT_0087 = IT_0035*IT_0086;
    const ccomplex_t IT_0088 = IT_0037*IT_0086;
    const ccomplex_t IT_0089 = IT_0011*IT_0051;
    const ccomplex_t IT_0090 = -IT_0052 + -IT_0063 + -IT_0064 + -IT_0075 + 
      -IT_0076 + -IT_0087 + -IT_0088 + -IT_0089;
    const ccomplex_t IT_0091 = IT_0006*IT_0049;
    const ccomplex_t IT_0092 = IT_0013*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0016*IT_0060;
    const ccomplex_t IT_0095 = IT_0019*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0021*IT_0094;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0032*IT_0084;
    const ccomplex_t IT_0100 = IT_0035*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0037*IT_0099;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0011*IT_0091;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0024*IT_0072;
    const ccomplex_t IT_0107 = IT_0027*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0029*IT_0106;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0093 + -IT_0096 + IT_0098 + -IT_0101 +
       IT_0103 + -IT_0105 + -IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = IT_0009*IT_0010;
    const ccomplex_t IT_0113 = s_23*s_24;
    const ccomplex_t IT_0114 = s_34*IT_0010;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = IT_0093 + IT_0096 + IT_0098 + IT_0101 + IT_0103
       + IT_0105 + IT_0108 + IT_0110;
    const ccomplex_t IT_0118 = m_mu*IT_0117;
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = conj(IT_0039) + conj(IT_0090);
    const ccomplex_t IT_0121 = (-2)*s_34;
    const ccomplex_t IT_0122 = conj(IT_0119)*IT_0121;
    const ccomplex_t IT_0123 = 2*IT_0119;
    const ccomplex_t IT_0124 = 2*conj(IT_0119);
    const ccomplex_t IT_0125 = (-2)*IT_0039;
    const ccomplex_t IT_0126 = (-2)*conj(IT_0039);
    const ccomplex_t IT_0127 = m_mu*s_24;
    const ccomplex_t IT_0128 = m_mu*s_23;
    const ccomplex_t IT_0129 = (-2)*IT_0119;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = (-2)*conj(IT_0119);
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = 2*s_34*(IT_0039*conj(IT_0039) + IT_0090*conj
      (IT_0090)) + (-4)*conj(IT_0111)*(IT_0111*IT_0112 + (-2)*IT_0111*IT_0116) +
       4*IT_0119*(s_34*conj(IT_0119) + 0.25*IT_0120*IT_0121) + (IT_0039 +
       IT_0090)*IT_0122 + (-4)*IT_0009*(IT_0119*conj(IT_0119) + (-0.25)*IT_0120
      *IT_0123 + (-0.25)*(IT_0039 + IT_0090)*IT_0124 + (-0.25)*conj(IT_0090)
      *IT_0125 + (-0.25)*IT_0090*IT_0126) + (-2)*(conj(IT_0111)*(IT_0090 + (-0.5
      )*IT_0123) + IT_0111*(conj(IT_0090) + (-0.5)*IT_0124) + (-0.5)*conj
      (IT_0111)*(IT_0123 + IT_0125) + (-0.5)*IT_0111*(IT_0124 + IT_0126))
      *IT_0127 + 2*IT_0128*(conj(IT_0111)*(IT_0039 + IT_0130) + conj(IT_0111)*
      (IT_0090 + IT_0130) + IT_0111*(conj(IT_0039) + IT_0132) + IT_0111*(conj
      (IT_0090) + IT_0132));
    return create_ccomplex_return(IT_0133);
}

