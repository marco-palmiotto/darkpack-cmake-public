#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_stau_1_to_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_stau_1_to_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_1, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = N_d1*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0007 + 0.5*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0020 = IT_0004*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_d1)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + IT_0024 + 
      -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0018*IT_0029;
    const ccomplex_t IT_0031 = IT_0002*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0038 = IT_0011*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_d3)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0041 = IT_0012*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0036 + IT_0039 + 
      -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0046 = IT_0004*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_d3*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0049 = IT_0012*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + 0.5*IT_0050);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = IT_0033*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0030*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0063 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0064 = IT_0004*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0067 = IT_0011*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_d2)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0070 = IT_0012*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + IT_0068 + 
      -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0075 = IT_0004*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = N_d2*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0078 = IT_0012*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + 0.5*IT_0079);
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = IT_0073*IT_0081;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0082*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0089 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0090 = IT_0004*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0093 = IT_0011*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_d4)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0096 = IT_0012*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0091 + IT_0094 + 
      -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0101 = IT_0004*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = N_d4*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0104 = IT_0012*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0102 + 0.5*IT_0105);
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0107;
    const ccomplex_t IT_0109 = IT_0088*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = 0.5*IT_0032 + 0.5*IT_0055 + (-0.5)*IT_0058 + (
      -0.5)*IT_0061 + 0.5*IT_0084 + (-0.5)*IT_0087 + 0.5*IT_0110 + (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = IT_0000*IT_0001;
    const ccomplex_t IT_0116 = s_23*s_24;
    const ccomplex_t IT_0117 = s_34*IT_0001;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = s_24*m_tau;
    const ccomplex_t IT_0121 = m_tau*IT_0032;
    const ccomplex_t IT_0122 = m_tau*IT_0084;
    const ccomplex_t IT_0123 = m_tau*IT_0087;
    const ccomplex_t IT_0124 = m_tau*IT_0055;
    const ccomplex_t IT_0125 = m_tau*IT_0058;
    const ccomplex_t IT_0126 = m_tau*IT_0110;
    const ccomplex_t IT_0127 = m_tau*IT_0113;
    const ccomplex_t IT_0128 = cpow(IT_0073, 2);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_2*IT_0128;
    const ccomplex_t IT_0130 = IT_0085*IT_0129;
    const ccomplex_t IT_0131 = cpow(IT_0044, 2);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_3*IT_0131;
    const ccomplex_t IT_0133 = IT_0056*IT_0132;
    const ccomplex_t IT_0134 = IT_0033*IT_0132;
    const ccomplex_t IT_0135 = cpow(IT_0099, 2);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_4*IT_0135;
    const ccomplex_t IT_0137 = IT_0111*IT_0136;
    const ccomplex_t IT_0138 = m_tau*IT_0061;
    const ccomplex_t IT_0139 = cpow(IT_0029, 2);
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*m_N_1*IT_0139;
    const ccomplex_t IT_0141 = IT_0059*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0121 + 0.5*IT_0122 + 0.5*IT_0123 + 0.5
      *IT_0124 + 0.5*IT_0125 + 0.5*IT_0126 + 0.5*IT_0127 + (-0.5)*IT_0130 + (
      -0.5)*IT_0133 + (-0.5)*IT_0134 + (-0.5)*IT_0137 + 0.5*IT_0138 + (-0.5)
      *IT_0141;
    const ccomplex_t IT_0143 = IT_0002*IT_0140;
    const ccomplex_t IT_0144 = IT_0062*IT_0129;
    const ccomplex_t IT_0145 = IT_0088*IT_0136;
    const ccomplex_t IT_0146 = 0.5*IT_0143 + 0.5*IT_0144 + 0.5*IT_0145;
    const ccomplex_t IT_0147 = IT_0142 + -IT_0146;
    const ccomplex_t IT_0148 = -conj(IT_0146);
    const ccomplex_t IT_0149 = conj(IT_0142) + IT_0148;
    const ccomplex_t IT_0150 = cpow(IT_0081, 2);
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*m_N_2*IT_0150;
    const ccomplex_t IT_0152 = IT_0062*IT_0151;
    const ccomplex_t IT_0153 = cpow(IT_0052, 2);
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*m_N_3*IT_0153;
    const ccomplex_t IT_0155 = IT_0033*IT_0154;
    const ccomplex_t IT_0156 = IT_0056*IT_0154;
    const ccomplex_t IT_0157 = cpow(IT_0107, 2);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_N_4*IT_0157;
    const ccomplex_t IT_0159 = IT_0111*IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0018, 2);
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_N_1*IT_0160;
    const ccomplex_t IT_0162 = IT_0059*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0121 + 0.5*IT_0122 + 0.5*IT_0123 + 0.5
      *IT_0124 + 0.5*IT_0125 + 0.5*IT_0126 + 0.5*IT_0127 + 0.5*IT_0138 + (-0.5)
      *IT_0152 + (-0.5)*IT_0155 + (-0.5)*IT_0156 + (-0.5)*IT_0159 + (-0.5)
      *IT_0162;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = IT_0085*IT_0151;
    const ccomplex_t IT_0166 = IT_0088*IT_0158;
    const ccomplex_t IT_0167 = IT_0002*IT_0161;
    const ccomplex_t IT_0168 = 0.5*IT_0165 + 0.5*IT_0166 + 0.5*IT_0167;
    const ccomplex_t IT_0169 = (-2)*IT_0168;
    const ccomplex_t IT_0170 = IT_0164 + IT_0169;
    const ccomplex_t IT_0171 = 2*conj(IT_0163);
    const ccomplex_t IT_0172 = (-2)*conj(IT_0168);
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0168;
    const ccomplex_t IT_0175 = (-2)*IT_0163;
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = 2*conj(IT_0168);
    const ccomplex_t IT_0178 = (-2)*conj(IT_0163);
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = s_23*m_tau;
    const ccomplex_t IT_0181 = 2*s_34;
    const ccomplex_t IT_0182 = (-2)*s_34;
    const ccomplex_t IT_0183 = (-4)*conj(IT_0114)*(IT_0114*IT_0115 + (-2)
      *IT_0114*IT_0119) + 2*IT_0120*(conj(IT_0114)*IT_0147 + IT_0114*IT_0149 +
       0.5*conj(IT_0114)*IT_0170 + 0.5*IT_0114*IT_0173) + IT_0000*(conj(IT_0146)
      *IT_0170 + IT_0146*IT_0173 + conj(IT_0142)*IT_0176 + IT_0142*IT_0179) + (
      -2)*(conj(IT_0114)*IT_0147 + IT_0114*IT_0149 + (-0.5)*conj(IT_0114)
      *IT_0176 + (-0.5)*IT_0114*IT_0179)*IT_0180 + (IT_0142*conj(IT_0142) +
       IT_0146*conj(IT_0146) + IT_0163*conj(IT_0163) + IT_0168*conj(IT_0168))
      *IT_0181 + (conj(IT_0142)*IT_0146 + IT_0142*conj(IT_0146) + conj(IT_0163)
      *IT_0168 + IT_0163*conj(IT_0168))*IT_0182;
    return create_ccomplex_return(IT_0183);
}

