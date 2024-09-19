#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_2_to_nu_tau_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_2_to_nu_tau_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
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
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = conj(U_d1)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0010 + (-0.5)*IT_0017);
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0021 = pow(m_stau_2, 2);
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = conj(U_d2)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0034 = IT_0013*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-0.5)*IT_0035);
    const ccomplex_t IT_0037 = IT_0030*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0039 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0037*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = pow(m_N_1, 2);
    const ccomplex_t IT_0042 = pow(m_tau, 2);
    const ccomplex_t IT_0043 = cpow(s_12 + 0.5*IT_0021 + 0.5*IT_0041 + (-0.5)
      *IT_0042 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = cos(theta_W);
    const ccomplex_t IT_0048 = cpow(IT_0047, -1);
    const ccomplex_t IT_0049 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = N_d1*e_em*m_tau*IT_0014*U_stau_01;
    const ccomplex_t IT_0053 = IT_0013*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + 0.5*IT_0054);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0046*IT_0056;
    const ccomplex_t IT_0058 = IT_0043*IT_0044*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0023 + 0.5*IT_0040 + 0.5*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0063 = IT_0001*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0036*IT_0066;
    const ccomplex_t IT_0068 = IT_0039*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0071 = IT_0001*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0073 = IT_0001*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + 0.5*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0018*IT_0076;
    const ccomplex_t IT_0078 = IT_0022*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0069 + IT_0079;
    const ccomplex_t IT_0081 = m_N_1*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = pow(m_W, -2);
    const ccomplex_t IT_0084 = pow(m_W, 2);
    const ccomplex_t IT_0085 = s_13*IT_0084;
    const ccomplex_t IT_0086 = s_14*s_34;
    const ccomplex_t IT_0087 = (-2)*IT_0086;
    const ccomplex_t IT_0088 = IT_0085 + IT_0087;
    const ccomplex_t IT_0089 = IT_0083*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = s_13 + IT_0090;
    const ccomplex_t IT_0092 = pow(s_34, 2);
    const ccomplex_t IT_0093 = m_N_1*IT_0083*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0095 = IT_0048*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0098 = IT_0001*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = conj(N_d1)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0101 = IT_0013*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + IT_0099 + 
      -IT_0102);
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0046*IT_0104;
    const ccomplex_t IT_0106 = IT_0043*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *U_stau_01;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = conj(N_B1)*e_em;
    const ccomplex_t IT_0111 = IT_0048*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W1)*e_em;
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = cpow(s_13 + (-0.5)*IT_0041 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0109*IT_0120;
    const ccomplex_t IT_0122 = -IT_0107 + (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = 2*IT_0122;
    const ccomplex_t IT_0124 = (-0.5)*IT_0093*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0107;
    const ccomplex_t IT_0126 = s_34*IT_0084;
    const ccomplex_t IT_0127 = IT_0083*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = s_34 + IT_0128;
    const ccomplex_t IT_0130 = m_N_1*IT_0129;
    const ccomplex_t IT_0131 = 2*IT_0130;
    const ccomplex_t IT_0132 = -IT_0125*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0069 + 0.5*IT_0079;
    const ccomplex_t IT_0134 = (-0.25)*IT_0127;
    const ccomplex_t IT_0135 = s_34 + IT_0134;
    const ccomplex_t IT_0136 = m_N_1*IT_0135;
    const ccomplex_t IT_0137 = (-4)*IT_0136;
    const ccomplex_t IT_0138 = IT_0133*IT_0137;
    const ccomplex_t IT_0139 = (-2)*(IT_0059 + IT_0082)*IT_0091 + IT_0124 +
       IT_0132 + IT_0138;
    const ccomplex_t IT_0140 = s_14*s_34*IT_0083;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = s_13 + IT_0141;
    const ccomplex_t IT_0143 = s_13*IT_0142;
    const ccomplex_t IT_0144 = -IT_0069 + -IT_0079;
    const ccomplex_t IT_0145 = 2*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0121;
    const ccomplex_t IT_0147 = (-2)*IT_0146;
    const ccomplex_t IT_0148 = IT_0145 + IT_0147;
    const ccomplex_t IT_0149 = 2*conj(IT_0144);
    const ccomplex_t IT_0150 = (-2)*conj(IT_0146);
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = m_N_1*IT_0142;
    const ccomplex_t IT_0153 = IT_0059 + IT_0082;
    const ccomplex_t IT_0154 = s_13*s_14;
    const ccomplex_t IT_0155 = s_34*IT_0041;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = IT_0125 + IT_0133;
    const ccomplex_t IT_0159 = s_13*IT_0083*IT_0092;
    const ccomplex_t IT_0160 = pow(m_W, 4);
    const ccomplex_t IT_0161 = s_13*IT_0160;
    const ccomplex_t IT_0162 = IT_0083*IT_0161;
    const ccomplex_t IT_0163 = IT_0085 + IT_0087 + (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = (-0.25)*IT_0162;
    const ccomplex_t IT_0165 = IT_0086 + IT_0164;
    const ccomplex_t IT_0166 = s_13*s_34;
    const ccomplex_t IT_0167 = 2*conj(IT_0122);
    const ccomplex_t IT_0168 = pow(s_14, 2);
    const ccomplex_t IT_0169 = IT_0083*IT_0168;
    const ccomplex_t IT_0170 = -IT_0041;
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = s_13*IT_0171;
    const ccomplex_t IT_0173 = (-2)*(conj(IT_0059) + conj(IT_0082))*IT_0139 +
       IT_0143*(conj(IT_0122)*IT_0148 + IT_0122*IT_0151) + IT_0152*((conj
      (IT_0059) + conj(IT_0082))*IT_0148 + IT_0151*IT_0153) + IT_0157*((conj
      (IT_0125) + conj(IT_0133))*IT_0148 + IT_0151*IT_0158) + conj(IT_0122)
      *IT_0123*IT_0159 + (-4)*conj(IT_0133)*(0.5*IT_0137*IT_0153 + IT_0133
      *IT_0163 + (-2)*IT_0125*IT_0165 + (-0.25)*IT_0123*IT_0166) + (-4)*conj
      (IT_0125)*((-0.5)*IT_0131*IT_0153 + IT_0125*IT_0163 + (-2)*IT_0133*IT_0165
       + (-0.25)*IT_0123*IT_0166) + (IT_0093*IT_0153 + IT_0158*IT_0166)*IT_0167 
      + (-2)*(conj(IT_0144)*((-0.5)*IT_0145 + IT_0146) + (IT_0144 + -IT_0146)
      *conj(IT_0146))*IT_0172;
    return create_ccomplex_return(IT_0173);
}

