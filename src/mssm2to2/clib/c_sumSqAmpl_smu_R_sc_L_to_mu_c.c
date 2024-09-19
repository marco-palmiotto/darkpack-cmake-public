#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_sc_L_to_mu_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_sc_L_to_mu_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
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
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
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
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = pow(m_mu, 2);
    const ccomplex_t IT_0016 = pow(m_smu_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0015 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = conj(N_B2)*e_em;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W2)*e_em;
    const ccomplex_t IT_0024 = IT_0006*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + 3*IT_0025);
    const ccomplex_t IT_0027 = 0.166666666666667*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0015 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B3)*e_em;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W3)*e_em;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 3*IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0015 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B4)*e_em;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W4)*e_em;
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0053);
    const ccomplex_t IT_0055 = 0.166666666666667*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0015 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0019 + -IT_0033 + -IT_0047 + -IT_0061;
    const ccomplex_t IT_0063 = pow(m_W, -1);
    const ccomplex_t IT_0064 = cos(beta);
    const ccomplex_t IT_0065 = cpow(IT_0064, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0006*IT_0063*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = sin(beta);
    const ccomplex_t IT_0069 = cpow(IT_0068, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0006*IT_0063*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0067*IT_0071;
    const ccomplex_t IT_0073 = IT_0017*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0006*IT_0063*IT_0065;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0006*IT_0063*IT_0069;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = IT_0045*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0006*IT_0063*IT_0065;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0006*IT_0063*IT_0069;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0059*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0006*IT_0063*IT_0065;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0006*IT_0063*IT_0069;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = IT_0031*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = -IT_0074 + -IT_0081 + -IT_0088 + -IT_0095;
    const ccomplex_t IT_0097 = m_c*m_mu*IT_0016;
    const ccomplex_t IT_0098 = s_34*IT_0016;
    const ccomplex_t IT_0099 = (-6)*IT_0098;
    const ccomplex_t IT_0100 = s_13*s_14;
    const ccomplex_t IT_0101 = 12*IT_0100;
    const ccomplex_t IT_0102 = IT_0099 + IT_0101;
    const ccomplex_t IT_0103 = m_c*s_13;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0105 = IT_0011*IT_0067;
    const ccomplex_t IT_0106 = IT_0017*IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0108 = IT_0055*IT_0083;
    const ccomplex_t IT_0109 = IT_0059*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = m_mu*IT_0019;
    const ccomplex_t IT_0111 = m_mu*IT_0033;
    const ccomplex_t IT_0112 = m_mu*IT_0047;
    const ccomplex_t IT_0113 = m_mu*IT_0061;
    const ccomplex_t IT_0114 = IT_0106 + IT_0109 + -IT_0110 + -IT_0111 + 
      -IT_0112 + -IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0116 = IT_0027*IT_0090;
    const ccomplex_t IT_0117 = IT_0031*IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0119 = IT_0041*IT_0076;
    const ccomplex_t IT_0120 = IT_0045*IT_0118*IT_0119;
    const ccomplex_t IT_0121 = -IT_0117 + -IT_0120;
    const ccomplex_t IT_0122 = IT_0114 + -IT_0121;
    const ccomplex_t IT_0123 = -conj(IT_0121);
    const ccomplex_t IT_0124 = conj(IT_0114) + IT_0123;
    const ccomplex_t IT_0125 = IT_0029*IT_0092;
    const ccomplex_t IT_0126 = IT_0031*IT_0115*IT_0125;
    const ccomplex_t IT_0127 = m_mu*IT_0081;
    const ccomplex_t IT_0128 = m_mu*IT_0088;
    const ccomplex_t IT_0129 = m_mu*IT_0074;
    const ccomplex_t IT_0130 = m_mu*IT_0095;
    const ccomplex_t IT_0131 = IT_0126 + -IT_0127 + -IT_0128 + -IT_0129 + 
      -IT_0130;
    const ccomplex_t IT_0132 = 6*IT_0131;
    const ccomplex_t IT_0133 = IT_0013*IT_0071;
    const ccomplex_t IT_0134 = IT_0017*IT_0104*IT_0133;
    const ccomplex_t IT_0135 = IT_0043*IT_0078;
    const ccomplex_t IT_0136 = IT_0045*IT_0118*IT_0135;
    const ccomplex_t IT_0137 = IT_0057*IT_0085;
    const ccomplex_t IT_0138 = IT_0059*IT_0107*IT_0137;
    const ccomplex_t IT_0139 = -IT_0134 + -IT_0136 + -IT_0138;
    const ccomplex_t IT_0140 = (-6)*IT_0139;
    const ccomplex_t IT_0141 = IT_0132 + IT_0140;
    const ccomplex_t IT_0142 = 6*conj(IT_0131);
    const ccomplex_t IT_0143 = (-6)*conj(IT_0139);
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = m_c*m_mu;
    const ccomplex_t IT_0146 = 6*IT_0139;
    const ccomplex_t IT_0147 = (-6)*IT_0131;
    const ccomplex_t IT_0148 = IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = 6*conj(IT_0139);
    const ccomplex_t IT_0150 = (-6)*conj(IT_0131);
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = 6*s_34;
    const ccomplex_t IT_0153 = (-6)*s_34;
    const ccomplex_t IT_0154 = (-0.166666666666667)*IT_0096;
    const ccomplex_t IT_0155 = m_mu*s_14;
    const ccomplex_t IT_0156 = (-6)*conj(IT_0062)*(IT_0096*IT_0097 + (
      -0.166666666666667)*IT_0062*IT_0102) + 6*IT_0103*(conj(IT_0096)*IT_0122 +
       IT_0096*IT_0124 + 0.166666666666667*conj(IT_0062)*IT_0141 +
       0.166666666666667*IT_0062*IT_0144) + IT_0145*(conj(IT_0121)*IT_0141 +
       IT_0121*IT_0144 + conj(IT_0114)*IT_0148 + IT_0114*IT_0151) + (IT_0114
      *conj(IT_0114) + IT_0121*conj(IT_0121) + IT_0131*conj(IT_0131) + IT_0139
      *conj(IT_0139))*IT_0152 + (conj(IT_0114)*IT_0121 + IT_0114*conj(IT_0121) +
       conj(IT_0131)*IT_0139 + IT_0131*conj(IT_0139))*IT_0153 + (-6)*conj
      (IT_0096)*(IT_0062*IT_0097 + IT_0102*IT_0154) + (-6)*(conj(IT_0062)
      *IT_0122 + IT_0062*IT_0124 + (-0.166666666666667)*conj(IT_0096)*IT_0148 +
       IT_0151*IT_0154)*IT_0155;
    return create_ccomplex_return(IT_0156);
}

