#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_snu_mu_to_mu_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_snu_mu_to_mu_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0001 = m_mu*m_N_1*IT_0000;
    const ccomplex_t IT_0002 = pow(m_mu, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 2*IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = IT_0017 + IT_0020;
    const ccomplex_t IT_0022 = IT_0011 + IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0006*IT_0024;
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0030 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0031 = IT_0008*IT_0030;
    const ccomplex_t IT_0032 = 2*IT_0031;
    const ccomplex_t IT_0033 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0034 = IT_0014*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0037 = IT_0008*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = IT_0035 + IT_0038;
    const ccomplex_t IT_0040 = IT_0032 + IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0029*IT_0042;
    const ccomplex_t IT_0044 = IT_0028*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = -IT_0027 + -IT_0045;
    const ccomplex_t IT_0047 = pow(m_N_1, 2);
    const ccomplex_t IT_0048 = cpow(s_12 + 0.5*IT_0000 + 0.5*IT_0047 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0049 = pow(m_W, -1);
    const ccomplex_t IT_0050 = cpow(IT_0007, -1);
    const ccomplex_t IT_0051 = sin(beta);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0005*IT_0049*IT_0050*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B1)*e_em;
    const ccomplex_t IT_0055 = IT_0013*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W1)*e_em;
    const ccomplex_t IT_0058 = IT_0005*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + -IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = IT_0048*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0005*IT_0049*IT_0050;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0005*IT_0051;
    const ccomplex_t IT_0069 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = IT_0013*IT_0051;
    const ccomplex_t IT_0073 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = IT_0071 + IT_0075;
    const ccomplex_t IT_0077 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0078 = IT_0068*IT_0077;
    const ccomplex_t IT_0079 = (-2)*IT_0078;
    const ccomplex_t IT_0080 = IT_0076 + IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0067*IT_0082;
    const ccomplex_t IT_0084 = IT_0028*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em*m_mu
      *IT_0005*IT_0049*IT_0050;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0089 = IT_0072*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0092 = IT_0068*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0096 = IT_0068*IT_0095;
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = IT_0094 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0087*IT_0100;
    const ccomplex_t IT_0102 = IT_0003*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = -IT_0085 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0065 + IT_0104;
    const ccomplex_t IT_0106 = conj(IT_0065) + conj(IT_0104);
    const ccomplex_t IT_0107 = s_12*s_23;
    const ccomplex_t IT_0108 = s_13*IT_0000;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107 + IT_0109;
    const ccomplex_t IT_0111 = m_mu*m_N_1;
    const ccomplex_t IT_0112 = IT_0027 + IT_0045;
    const ccomplex_t IT_0113 = m_mu*IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = m_N_1*IT_0064;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0117 = IT_0024*IT_0087;
    const ccomplex_t IT_0118 = IT_0003*IT_0116*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_13 + IT_0002 + IT_0047 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005*
      (m_W*IT_0007 + (-0.5)*IT_0002*IT_0049*IT_0050)*IT_0051;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0005*IT_0049*IT_0050;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_13 + IT_0002 + IT_0047 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0005*IT_0049;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0013;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0136 = IT_0042*IT_0067;
    const ccomplex_t IT_0137 = IT_0028*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0115 + -IT_0118 + IT_0126 + IT_0134 +
       -IT_0137;
    const ccomplex_t IT_0139 = IT_0085 + IT_0103;
    const ccomplex_t IT_0140 = m_mu*IT_0139;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0029*IT_0082;
    const ccomplex_t IT_0143 = IT_0028*IT_0135*IT_0142;
    const ccomplex_t IT_0144 = IT_0006*IT_0100;
    const ccomplex_t IT_0145 = IT_0003*IT_0116*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0005*IT_0049*IT_0050;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = IT_0129*IT_0147;
    const ccomplex_t IT_0149 = IT_0127*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0056 + IT_0059);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0121*IT_0152;
    const ccomplex_t IT_0154 = IT_0119*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = -IT_0143 + -IT_0145 + -IT_0150 + -IT_0155;
    const ccomplex_t IT_0157 = 2*IT_0138;
    const ccomplex_t IT_0158 = (-2)*IT_0114;
    const ccomplex_t IT_0159 = IT_0157 + IT_0158;
    const ccomplex_t IT_0160 = 2*conj(IT_0138);
    const ccomplex_t IT_0161 = (-2)*conj(IT_0114);
    const ccomplex_t IT_0162 = IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = s_23*m_N_1;
    const ccomplex_t IT_0164 = 2*IT_0156;
    const ccomplex_t IT_0165 = (-2)*IT_0141;
    const ccomplex_t IT_0166 = IT_0164 + IT_0165;
    const ccomplex_t IT_0167 = 2*conj(IT_0156);
    const ccomplex_t IT_0168 = (-2)*conj(IT_0141);
    const ccomplex_t IT_0169 = IT_0167 + IT_0168;
    const ccomplex_t IT_0170 = m_mu*s_12;
    const ccomplex_t IT_0171 = 2*s_13;
    const ccomplex_t IT_0172 = (-2)*s_13;
    const ccomplex_t IT_0173 = 2*IT_0001*(conj(IT_0046)*IT_0105 + IT_0046
      *IT_0106) + 4*(IT_0046*conj(IT_0046) + IT_0105*IT_0106)*IT_0110 + 2
      *IT_0111*((conj(IT_0114) + -conj(IT_0138))*IT_0141 + (IT_0114 + -IT_0138)
      *conj(IT_0141) + 0.5*conj(IT_0156)*IT_0159 + 0.5*IT_0156*IT_0162) +
       IT_0163*(IT_0106*IT_0159 + IT_0105*IT_0162 + conj(IT_0046)*IT_0166 +
       IT_0046*IT_0169) + (conj(IT_0046)*IT_0159 + IT_0046*IT_0162 + IT_0106
      *IT_0166 + IT_0105*IT_0169)*IT_0170 + (IT_0114*conj(IT_0114) + IT_0138
      *conj(IT_0138) + IT_0141*conj(IT_0141) + IT_0156*conj(IT_0156))*IT_0171 + 
      (conj(IT_0114)*IT_0138 + IT_0114*conj(IT_0138) + conj(IT_0141)*IT_0156 +
       IT_0141*conj(IT_0156))*IT_0172;
    return create_ccomplex_return(IT_0173);
}

