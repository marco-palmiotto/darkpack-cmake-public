#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_sb_1_to_nu_e_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_sb_1_to_nu_e_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0001 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W2)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + -IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*conj(N_d2)*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0023 = IT_0009*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + -IT_0024 +
       0.333333333333333*IT_0027);
    const ccomplex_t IT_0029 = 3*IT_0028;
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = IT_0000*IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0034 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0035 = conj(N_B4)*e_em;
    const ccomplex_t IT_0036 = IT_0004*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W4)*e_em;
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = m_b*conj(N_d4)*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0050 = IT_0009*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0045 + 0.333333333333333
      *IT_0048 + -IT_0051);
    const ccomplex_t IT_0053 = 3*IT_0052;
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0042*IT_0054;
    const ccomplex_t IT_0056 = IT_0033*IT_0034*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0058 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B1)*e_em;
    const ccomplex_t IT_0060 = IT_0004*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W1)*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = m_b*conj(N_d1)*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0068 = IT_0017*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0074 = IT_0004*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0069 + -IT_0072 +
       0.333333333333333*IT_0075);
    const ccomplex_t IT_0077 = 3*IT_0076;
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0066*IT_0078;
    const ccomplex_t IT_0080 = IT_0057*IT_0058*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0082 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0083 = m_b*conj(N_d3)*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0084 = IT_0017*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0090 = IT_0009*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + 0.333333333333333
      *IT_0088 + -IT_0091);
    const ccomplex_t IT_0093 = 3*IT_0092;
    const ccomplex_t IT_0094 = 0.166666666666667*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B3)*e_em;
    const ccomplex_t IT_0096 = IT_0004*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W3)*e_em;
    const ccomplex_t IT_0099 = IT_0009*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + -IT_0100);
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0094*IT_0102;
    const ccomplex_t IT_0104 = IT_0081*IT_0082*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0032 + (-0.5)*IT_0056 + (-0.5)
      *IT_0080 + (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0107 = IT_0004*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = m_b*N_d1*e_em*IT_0018*U_sb_00;
    const ccomplex_t IT_0110 = IT_0017*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + 1.5*IT_0111);
    const ccomplex_t IT_0113 = (-0.333333333333333)*IT_0112;
    const ccomplex_t IT_0114 = IT_0066*IT_0113;
    const ccomplex_t IT_0115 = IT_0058*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0118 = IT_0004*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = m_b*N_d2*e_em*IT_0018*U_sb_00;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + 1.5*IT_0122);
    const ccomplex_t IT_0124 = (-0.333333333333333)*IT_0123;
    const ccomplex_t IT_0125 = IT_0014*IT_0124;
    const ccomplex_t IT_0126 = IT_0002*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0129 = IT_0004*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = m_b*N_d3*e_em*IT_0018*U_sb_00;
    const ccomplex_t IT_0132 = IT_0017*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + 1.5*IT_0133);
    const ccomplex_t IT_0135 = (-0.333333333333333)*IT_0134;
    const ccomplex_t IT_0136 = IT_0102*IT_0135;
    const ccomplex_t IT_0137 = IT_0082*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0140 = IT_0004*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = m_b*N_d4*e_em*IT_0018*U_sb_00;
    const ccomplex_t IT_0143 = IT_0017*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + 1.5*IT_0144);
    const ccomplex_t IT_0146 = (-0.333333333333333)*IT_0145;
    const ccomplex_t IT_0147 = IT_0042*IT_0146;
    const ccomplex_t IT_0148 = IT_0034*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0116 + 0.5*IT_0127 + 0.5*IT_0138 + 0.5
      *IT_0149;
    const ccomplex_t IT_0151 = m_b*s_13;
    const ccomplex_t IT_0152 = s_13*s_14;
    const ccomplex_t IT_0153 = s_34*IT_0001;
    return create_ccomplex_return(6*IT_0105*(s_34*conj(IT_0105) + conj(IT_0150
      )*IT_0151) + 6*IT_0150*(conj(IT_0105)*IT_0151 + 0.166666666666667*conj
      (IT_0150)*(12*IT_0152 + (-6)*IT_0153)));
}

