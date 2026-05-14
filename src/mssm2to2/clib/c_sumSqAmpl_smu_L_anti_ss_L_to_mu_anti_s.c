#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_ss_L_to_mu_anti_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_ss_L_to_mu_anti_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = N_B1*e_em;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = N_W1*e_em;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_mu, 2);
    const ccomplex_t IT_0022 = pow(m_smu_L, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B2)*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W2)*e_em;
    const ccomplex_t IT_0030 = IT_0006*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + IT_0031);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = N_B2*e_em;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_W2*e_em;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-3)*IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B4)*e_em;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W4)*e_em;
    const ccomplex_t IT_0050 = IT_0006*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = N_B4*e_em;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W4*e_em;
    const ccomplex_t IT_0058 = IT_0006*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + (-3)*IT_0059);
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0067 = N_B3*e_em;
    const ccomplex_t IT_0068 = IT_0001*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_W3*e_em;
    const ccomplex_t IT_0071 = IT_0006*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-3)*IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B3)*e_em;
    const ccomplex_t IT_0076 = IT_0001*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W3)*e_em;
    const ccomplex_t IT_0079 = IT_0006*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = IT_0066*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0025 + -IT_0045 + -IT_0065 + -IT_0085;
    const ccomplex_t IT_0087 = pow(m_W, -1);
    const ccomplex_t IT_0088 = cos(beta);
    const ccomplex_t IT_0089 = cpow(IT_0088, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0041*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0094 = IT_0043*IT_0092*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0074*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0099 = IT_0066*IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0061*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0104 = IT_0063*IT_0102*IT_0103;
    const ccomplex_t IT_0105 = m_mu*IT_0025;
    const ccomplex_t IT_0106 = m_mu*IT_0045;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0019*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0111 = IT_0023*IT_0109*IT_0110;
    const ccomplex_t IT_0112 = m_mu*IT_0085;
    const ccomplex_t IT_0113 = m_mu*IT_0065;
    const ccomplex_t IT_0114 = -IT_0094 + -IT_0099 + -IT_0104 + IT_0105 +
       IT_0106 + -IT_0111 + IT_0112 + IT_0113;
    const ccomplex_t IT_0115 = m_mu*s_14;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = IT_0091*IT_0117;
    const ccomplex_t IT_0119 = IT_0043*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = m_mu*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0101*IT_0123;
    const ccomplex_t IT_0125 = IT_0063*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = m_mu*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = IT_0011*IT_0129;
    const ccomplex_t IT_0131 = IT_0023*IT_0110*IT_0130;
    const ccomplex_t IT_0132 = IT_0033*IT_0117;
    const ccomplex_t IT_0133 = IT_0043*IT_0093*IT_0132;
    const ccomplex_t IT_0134 = IT_0053*IT_0123;
    const ccomplex_t IT_0135 = IT_0063*IT_0103*IT_0134;
    const ccomplex_t IT_0136 = IT_0108*IT_0129;
    const ccomplex_t IT_0137 = IT_0023*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = m_mu*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0006*IT_0087*IT_0089;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0096*IT_0141;
    const ccomplex_t IT_0143 = IT_0066*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = m_mu*IT_0144;
    const ccomplex_t IT_0146 = IT_0082*IT_0141;
    const ccomplex_t IT_0147 = IT_0066*IT_0098*IT_0146;
    const ccomplex_t IT_0148 = IT_0121 + IT_0127 + -IT_0131 + -IT_0133 + 
      -IT_0135 + IT_0139 + IT_0145 + -IT_0147;
    const ccomplex_t IT_0149 = m_s*s_13;
    const ccomplex_t IT_0150 = -IT_0120 + -IT_0126 + -IT_0138 + -IT_0144;
    const ccomplex_t IT_0151 = m_s*m_mu*IT_0022;
    const ccomplex_t IT_0152 = s_34*IT_0022;
    const ccomplex_t IT_0153 = (-6)*IT_0152;
    const ccomplex_t IT_0154 = s_13*s_14;
    const ccomplex_t IT_0155 = 12*IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = m_s*m_mu;
    const ccomplex_t IT_0158 = 6*s_34;
    const ccomplex_t IT_0159 = -conj(IT_0086);
    const ccomplex_t IT_0160 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0161 = 6*IT_0086*(conj(IT_0114)*IT_0115 + -conj
      (IT_0148)*IT_0149 + -conj(IT_0150)*IT_0151 + 0.166666666666667*conj
      (IT_0086)*IT_0156) + (-6)*IT_0114*(IT_0149*conj(IT_0150) + conj(IT_0148)
      *IT_0157 + (-0.166666666666667)*conj(IT_0114)*IT_0158 + IT_0115*IT_0159) +
       6*IT_0150*(IT_0115*conj(IT_0148) + 0.166666666666667*conj(IT_0150)
      *IT_0156 + IT_0151*IT_0159 + 0.166666666666667*IT_0149*IT_0160) + 6
      *IT_0148*(IT_0115*conj(IT_0150) + 0.166666666666667*conj(IT_0148)*IT_0158 
      + IT_0149*IT_0159 + 0.166666666666667*IT_0157*IT_0160);
    return create_ccomplex_return(IT_0161);
}

