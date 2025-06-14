#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_ss_R_to_c_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_ss_R_to_c_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0006;
    const ccomplex_t IT_0008 = (-0.333333333333333)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0009*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0008*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0004*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0010, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0009*IT_0014*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = m_N_1*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *mu_h*IT_0009*IT_0014;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0009*IT_0012*IT_0014;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = pow(m_c, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*V_cs
      *e_em*IT_0009*IT_0012*IT_0014*IT_0019;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = 0.666666666666667*IT_0007;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0002 + IT_0031 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0009*IT_0014*IT_0019;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0011*IT_0014;
    const ccomplex_t IT_0046 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (-2)*IT_0047;
    const ccomplex_t IT_0049 = IT_0006*IT_0011;
    const ccomplex_t IT_0050 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0054 = IT_0045*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = IT_0052 + IT_0055;
    const ccomplex_t IT_0057 = IT_0048 + IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = IT_0044*IT_0059;
    const ccomplex_t IT_0061 = IT_0042*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = m_c*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0002 + IT_0031 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0065 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0066 = IT_0045*IT_0065;
    const ccomplex_t IT_0067 = (-2)*IT_0066;
    const ccomplex_t IT_0068 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0069 = IT_0049*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0072 = IT_0045*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = IT_0070 + IT_0073;
    const ccomplex_t IT_0075 = IT_0067 + IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0009*IT_0014*IT_0019;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = IT_0064*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = m_c*IT_0082;
    const ccomplex_t IT_0084 = -IT_0018 + IT_0025 + IT_0034 + IT_0041 + 
      -IT_0063 + -IT_0083;
    const ccomplex_t IT_0085 = conj(N_B1)*e_em;
    const ccomplex_t IT_0086 = IT_0006*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W1)*e_em;
    const ccomplex_t IT_0089 = IT_0014*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + 3*IT_0090);
    const ccomplex_t IT_0092 = 0.166666666666667*IT_0091;
    const ccomplex_t IT_0093 = IT_0027*IT_0092;
    const ccomplex_t IT_0094 = IT_0032*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0097 = IT_0010*IT_0014;
    const ccomplex_t IT_0098 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = 2*IT_0099;
    const ccomplex_t IT_0101 = IT_0006*IT_0010;
    const ccomplex_t IT_0102 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0106 = IT_0097*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = IT_0104 + IT_0107;
    const ccomplex_t IT_0109 = IT_0100 + IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = IT_0044*IT_0111;
    const ccomplex_t IT_0113 = IT_0042*IT_0096*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0009*IT_0011*IT_0014*IT_0019;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0021*IT_0115;
    const ccomplex_t IT_0117 = IT_0003*IT_0004*IT_0116;
    const ccomplex_t IT_0118 = IT_0008*IT_0115;
    const ccomplex_t IT_0119 = IT_0003*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = m_N_1*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0123 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0124 = IT_0097*IT_0123;
    const ccomplex_t IT_0125 = 2*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0127 = IT_0101*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0130 = IT_0097*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = IT_0125 + IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0079*IT_0135;
    const ccomplex_t IT_0137 = IT_0064*IT_0122*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0009*IT_0012*IT_0014;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = IT_0036*IT_0139;
    const ccomplex_t IT_0141 = IT_0039*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0095 + IT_0113 + IT_0117 + -IT_0121 +
       IT_0137 + -IT_0142;
    const ccomplex_t IT_0144 = 6*s_13;
    const ccomplex_t IT_0145 = m_c*m_N_1;
    const ccomplex_t IT_0146 = 6*IT_0084;
    const ccomplex_t IT_0147 = IT_0062 + IT_0082 + -IT_0120;
    const ccomplex_t IT_0148 = m_c*s_12;
    const ccomplex_t IT_0149 = s_23*m_N_1;
    const ccomplex_t IT_0150 = 6*conj(IT_0084);
    const ccomplex_t IT_0151 = 6*conj(IT_0143);
    const ccomplex_t IT_0152 = m_c*m_N_1*IT_0002;
    const ccomplex_t IT_0153 = s_13*IT_0002;
    const ccomplex_t IT_0154 = (-6)*IT_0153;
    const ccomplex_t IT_0155 = s_12*s_23;
    const ccomplex_t IT_0156 = 12*IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = 6*IT_0143;
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0148;
    const ccomplex_t IT_0160 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0161 = (IT_0084*conj(IT_0084) + IT_0143*conj(IT_0143))
      *IT_0144 + conj(IT_0143)*IT_0145*IT_0146 + (IT_0143*IT_0145 + IT_0147
      *IT_0148 + IT_0024*IT_0149)*IT_0150 + (IT_0024*IT_0148 + IT_0147*IT_0149)
      *IT_0151 + 6*conj(IT_0024)*(IT_0147*IT_0152 + 0.166666666666667*IT_0024
      *IT_0157 + IT_0158*IT_0159 + IT_0146*IT_0160) + 6*conj(IT_0147)*(IT_0024
      *IT_0152 + 0.166666666666667*IT_0147*IT_0157 + IT_0146*IT_0159 + IT_0158
      *IT_0160);
    return create_ccomplex_return(IT_0161);
}

