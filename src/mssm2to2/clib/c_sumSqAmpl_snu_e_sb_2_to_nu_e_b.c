#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_sb_2_to_nu_e_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_sb_2_to_nu_e_b(
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
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
    const ccomplex_t IT_0019 = m_b*conj(N_d2)*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0023 = IT_0009*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + -IT_0024 +
       0.333333333333333*IT_0027);
    const ccomplex_t IT_0029 = 3*IT_0028;
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = IT_0000*IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B1)*e_em;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W1)*e_em;
    const ccomplex_t IT_0038 = IT_0009*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0043 = IT_0004*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0046 = IT_0009*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = m_b*conj(N_d1)*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0044 + (-3)*IT_0047 + 3
      *IT_0050);
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0041*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0055 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0053*IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0058 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B4)*e_em;
    const ccomplex_t IT_0060 = IT_0004*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W4)*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = m_b*conj(N_d4)*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0068 = IT_0017*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B4)*e_em*U_sb_01;
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
    const ccomplex_t IT_0083 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0084 = IT_0004*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0087 = IT_0009*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_b*conj(N_d3)*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0090 = IT_0017*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + (-3)*IT_0088 + 3
      *IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B3)*e_em;
    const ccomplex_t IT_0095 = IT_0004*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W3)*e_em;
    const ccomplex_t IT_0098 = IT_0009*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + -IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0093*IT_0101;
    const ccomplex_t IT_0103 = IT_0081*IT_0082*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0056 + 0.5*IT_0080 + 0.5*IT_0103;
    const ccomplex_t IT_0105 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0106 = IT_0004*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = m_b*N_d1*e_em*IT_0018*U_sb_01;
    const ccomplex_t IT_0109 = IT_0017*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + 1.5*IT_0110);
    const ccomplex_t IT_0112 = (-0.333333333333333)*IT_0111;
    const ccomplex_t IT_0113 = IT_0041*IT_0112;
    const ccomplex_t IT_0114 = IT_0055*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0117 = IT_0004*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = m_b*N_d2*e_em*IT_0018*U_sb_01;
    const ccomplex_t IT_0120 = IT_0017*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + 1.5*IT_0121);
    const ccomplex_t IT_0123 = (-0.333333333333333)*IT_0122;
    const ccomplex_t IT_0124 = IT_0014*IT_0123;
    const ccomplex_t IT_0125 = IT_0002*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0128 = IT_0004*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_b*N_d3*e_em*IT_0018*U_sb_01;
    const ccomplex_t IT_0131 = IT_0017*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + 1.5*IT_0132);
    const ccomplex_t IT_0134 = (-0.333333333333333)*IT_0133;
    const ccomplex_t IT_0135 = IT_0101*IT_0134;
    const ccomplex_t IT_0136 = IT_0082*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0139 = IT_0004*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = m_b*N_d4*e_em*IT_0018*U_sb_01;
    const ccomplex_t IT_0142 = IT_0017*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 1.5*IT_0143);
    const ccomplex_t IT_0145 = (-0.333333333333333)*IT_0144;
    const ccomplex_t IT_0146 = IT_0066*IT_0145;
    const ccomplex_t IT_0147 = IT_0058*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0115 + 0.5*IT_0126 + 0.5*IT_0137 + 0.5
      *IT_0148;
    const ccomplex_t IT_0150 = m_b*s_13;
    const ccomplex_t IT_0151 = s_13*s_14;
    const ccomplex_t IT_0152 = s_34*IT_0001;
    const ccomplex_t IT_0153 = (-6)*s_34;
    const ccomplex_t IT_0154 = 6*s_34*(IT_0033*conj(IT_0033) + IT_0104*conj
      (IT_0104)) + 6*((conj(IT_0033) + -conj(IT_0104))*IT_0149 + (IT_0033 + 
      -IT_0104)*conj(IT_0149))*IT_0150 + IT_0149*conj(IT_0149)*(12*IT_0151 + (-6
      )*IT_0152) + (conj(IT_0033)*IT_0104 + IT_0033*conj(IT_0104))*IT_0153;
    return create_ccomplex_return(IT_0154);
}

