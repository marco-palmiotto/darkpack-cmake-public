#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_sb_1_to_nu_mu_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_sb_1_to_nu_mu_anti_b(
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
    const ccomplex_t IT_0015 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0016 = IT_0004*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0009*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_b*conj(N_d2)*e_em*IT_0021*conj(U_sb_00);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0017 + 1.5*IT_0024);
    const ccomplex_t IT_0026 = (-0.333333333333333)*IT_0025;
    const ccomplex_t IT_0027 = IT_0014*IT_0026;
    const ccomplex_t IT_0028 = IT_0000*IT_0002*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0030 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = conj(N_B4)*e_em;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W4)*e_em;
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0040 = IT_0004*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = m_b*conj(N_d4)*e_em*IT_0021*conj(U_sb_00);
    const ccomplex_t IT_0043 = IT_0020*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + 1.5*IT_0044);
    const ccomplex_t IT_0046 = (-0.333333333333333)*IT_0045;
    const ccomplex_t IT_0047 = IT_0038*IT_0046;
    const ccomplex_t IT_0048 = IT_0029*IT_0030*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0050 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0051 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_b*conj(N_d3)*e_em*IT_0021*conj(U_sb_00);
    const ccomplex_t IT_0055 = IT_0020*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 1.5*IT_0056);
    const ccomplex_t IT_0058 = (-0.333333333333333)*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B3)*e_em;
    const ccomplex_t IT_0060 = IT_0004*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = IT_0049*IT_0050*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0070 = cpow(s_13 + (-0.5)*IT_0001 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0071 = conj(N_B1)*e_em;
    const ccomplex_t IT_0072 = IT_0004*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W1)*e_em;
    const ccomplex_t IT_0075 = IT_0009*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + -IT_0076);
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0080 = IT_0004*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = m_b*conj(N_d1)*e_em*IT_0021*conj(U_sb_00);
    const ccomplex_t IT_0083 = IT_0020*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 1.5*IT_0084);
    const ccomplex_t IT_0086 = (-0.333333333333333)*IT_0085;
    const ccomplex_t IT_0087 = IT_0078*IT_0086;
    const ccomplex_t IT_0088 = IT_0069*IT_0070*IT_0087;
    const ccomplex_t IT_0089 = 0.5*IT_0028 + 0.5*IT_0048 + 0.5*IT_0068 + 0.5
      *IT_0088;
    const ccomplex_t IT_0090 = m_b*N_d1*e_em*IT_0021*conj(U_sb_10);
    const ccomplex_t IT_0091 = IT_0020*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0094 = IT_0009*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0097 = IT_0004*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + -IT_0095 +
       0.333333333333333*IT_0098);
    const ccomplex_t IT_0100 = 3*IT_0099;
    const ccomplex_t IT_0101 = 0.166666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0078*IT_0101;
    const ccomplex_t IT_0103 = IT_0070*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_b*N_d2*e_em*IT_0021*conj(U_sb_10);
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0109 = IT_0009*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0112 = IT_0004*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0107 + -IT_0110 +
       0.333333333333333*IT_0113);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = (-3)*IT_0115;
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0118 = IT_0014*IT_0117;
    const ccomplex_t IT_0119 = IT_0002*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = m_b*N_d3*e_em*IT_0021*conj(U_sb_10);
    const ccomplex_t IT_0122 = IT_0020*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0125 = IT_0004*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0128 = IT_0009*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0123 + 0.333333333333333
      *IT_0126 + -IT_0129);
    const ccomplex_t IT_0131 = 3*IT_0130;
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0131;
    const ccomplex_t IT_0133 = IT_0066*IT_0132;
    const ccomplex_t IT_0134 = IT_0050*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0137 = IT_0004*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0140 = IT_0009*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = m_b*N_d4*e_em*IT_0021*conj(U_sb_10);
    const ccomplex_t IT_0143 = IT_0020*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0138 + (-3)*IT_0141 + 3
      *IT_0144);
    const ccomplex_t IT_0146 = 0.166666666666667*IT_0145;
    const ccomplex_t IT_0147 = IT_0038*IT_0146;
    const ccomplex_t IT_0148 = IT_0030*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0104 + 0.5*IT_0120 + 0.5*IT_0135 + 0.5
      *IT_0149;
    const ccomplex_t IT_0151 = m_b*s_13;
    const ccomplex_t IT_0152 = s_13*s_14;
    const ccomplex_t IT_0153 = s_34*IT_0001;
    return create_ccomplex_return(6*IT_0089*(s_34*conj(IT_0089) + -conj
      (IT_0150)*IT_0151) + (-6)*IT_0150*(conj(IT_0089)*IT_0151 + (
      -0.166666666666667)*conj(IT_0150)*(12*IT_0152 + (-6)*IT_0153)));
}

