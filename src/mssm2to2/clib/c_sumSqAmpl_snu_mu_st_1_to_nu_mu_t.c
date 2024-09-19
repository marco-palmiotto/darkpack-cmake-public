#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_st_1_to_nu_mu_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_st_1_to_nu_mu_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
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
    const ccomplex_t IT_0012 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_t*conj(N_u1)*e_em*IT_0021*U_st_10;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0014 + 3*IT_0017 + 3
      *IT_0024);
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0029 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0030 = cpow(s_13 + (-0.5)*IT_0029 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0027*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B2)*e_em;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W2)*e_em;
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = m_t*conj(N_u2)*e_em*IT_0021*U_st_10;
    const ccomplex_t IT_0047 = IT_0020*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0042 + 3*IT_0045 + 3
      *IT_0048);
    const ccomplex_t IT_0050 = 0.166666666666667*IT_0049;
    const ccomplex_t IT_0051 = IT_0039*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0053 = cpow(s_13 + (-0.5)*IT_0029 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B3)*e_em;
    const ccomplex_t IT_0056 = IT_0001*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W3)*e_em;
    const ccomplex_t IT_0059 = IT_0006*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + -IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0064 = IT_0001*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0067 = IT_0006*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = m_t*conj(N_u3)*e_em*IT_0021*U_st_10;
    const ccomplex_t IT_0070 = IT_0020*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + 3*IT_0068 + 3
      *IT_0071);
    const ccomplex_t IT_0073 = 0.166666666666667*IT_0072;
    const ccomplex_t IT_0074 = IT_0062*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0076 = cpow(s_13 + (-0.5)*IT_0029 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0074*IT_0075*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B4)*e_em;
    const ccomplex_t IT_0079 = IT_0001*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W4)*e_em;
    const ccomplex_t IT_0082 = IT_0006*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0090 = IT_0006*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = m_t*conj(N_u4)*e_em*IT_0021*U_st_10;
    const ccomplex_t IT_0093 = IT_0020*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091 + 3
      *IT_0094);
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0097 = IT_0085*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0099 = cpow(s_13 + (-0.5)*IT_0029 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0097*IT_0098*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0031 + 0.5*IT_0054 + 0.5*IT_0077 + 0.5
      *IT_0100;
    const ccomplex_t IT_0102 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = m_t*N_u1*e_em*IT_0021*U_st_00;
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + (-0.75)*IT_0107);
    const ccomplex_t IT_0109 = 0.666666666666667*IT_0108;
    const ccomplex_t IT_0110 = IT_0011*IT_0109;
    const ccomplex_t IT_0111 = IT_0030*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = m_t*N_u2*e_em*IT_0021*U_st_00;
    const ccomplex_t IT_0117 = IT_0020*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.75)*IT_0118);
    const ccomplex_t IT_0120 = 0.666666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0039*IT_0120;
    const ccomplex_t IT_0122 = IT_0053*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0125 = IT_0001*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_t*N_u3*e_em*IT_0021*U_st_00;
    const ccomplex_t IT_0128 = IT_0020*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.75)*IT_0129);
    const ccomplex_t IT_0131 = 0.666666666666667*IT_0130;
    const ccomplex_t IT_0132 = IT_0062*IT_0131;
    const ccomplex_t IT_0133 = IT_0076*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0136 = IT_0001*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = m_t*N_u4*e_em*IT_0021*U_st_00;
    const ccomplex_t IT_0139 = IT_0020*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + (-0.75)*IT_0140);
    const ccomplex_t IT_0142 = 0.666666666666667*IT_0141;
    const ccomplex_t IT_0143 = IT_0085*IT_0142;
    const ccomplex_t IT_0144 = IT_0099*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = IT_0112 + IT_0123 + IT_0134 + IT_0145;
    const ccomplex_t IT_0147 = m_t*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0112 + 0.5*IT_0123 + 0.5*IT_0134 + 0.5
      *IT_0145;
    const ccomplex_t IT_0150 = m_t*s_23;
    const ccomplex_t IT_0151 = s_23*s_24;
    const ccomplex_t IT_0152 = pow(m_st_1, 2);
    const ccomplex_t IT_0153 = s_34*IT_0152;
    const ccomplex_t IT_0154 = (-6)*s_34;
    const ccomplex_t IT_0155 = 6*s_34*(IT_0101*conj(IT_0101) + IT_0148*conj
      (IT_0148)) + 6*((conj(IT_0101) + -conj(IT_0148))*IT_0149 + (IT_0101 + 
      -IT_0148)*conj(IT_0149))*IT_0150 + IT_0149*conj(IT_0149)*(12*IT_0151 + (-6
      )*IT_0153) + (conj(IT_0101)*IT_0148 + IT_0101*conj(IT_0148))*IT_0154;
    return create_ccomplex_return(IT_0155);
}

