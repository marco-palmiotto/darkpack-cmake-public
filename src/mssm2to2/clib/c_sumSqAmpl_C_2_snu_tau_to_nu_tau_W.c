#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_snu_tau_to_nu_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_snu_tau_to_nu_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
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
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0001 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = conj(N_B1)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_W1)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0017 = IT_0008*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0015 + 0.5*IT_0018);
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B4)*e_em;
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W4)*e_em;
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0033 = IT_0008*IT_0032;
    const ccomplex_t IT_0034 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0035 = IT_0008*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = conj(N_B2)*e_em;
    const ccomplex_t IT_0045 = IT_0003*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W2)*e_em;
    const ccomplex_t IT_0048 = IT_0008*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + -IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0053 = IT_0008*IT_0052;
    const ccomplex_t IT_0054 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0055 = IT_0008*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 0.5*IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0051*IT_0058;
    const ccomplex_t IT_0060 = IT_0043*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B3)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W3)*e_em;
    const ccomplex_t IT_0066 = IT_0008*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0073 = IT_0008*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + 0.5*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0069*IT_0076;
    const ccomplex_t IT_0078 = cpow(s_23 + (-0.5)*IT_0000 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0023 + (-0.5)*IT_0042 + (-0.5)
      *IT_0061 + (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008
      *conj(U_stau_01);
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0085 = IT_0008*IT_0084;
    const ccomplex_t IT_0086 = cos(beta);
    const ccomplex_t IT_0087 = cpow(IT_0086, -1);
    const ccomplex_t IT_0088 = IT_0008*IT_0087;
    const ccomplex_t IT_0089 = pow(m_W, -1);
    const ccomplex_t IT_0090 = conj(U_d2)*e_em*m_tau*IT_0089*U_stau_11;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0085 + (-0.5)*IT_0092);
    const ccomplex_t IT_0094 = pow(m_C_2, 2);
    const ccomplex_t IT_0095 = cpow(s_13 + (-0.5)*IT_0094 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0083*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008
      *conj(U_stau_00);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0102 = IT_0008*IT_0101;
    const ccomplex_t IT_0103 = conj(U_d2)*e_em*m_tau*IT_0089*U_stau_10;
    const ccomplex_t IT_0104 = IT_0088*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0102 + (-0.5)*IT_0105);
    const ccomplex_t IT_0107 = cpow(s_13 + (-0.5)*IT_0094 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0100*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0098 + (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = s_13*s_34;
    const ccomplex_t IT_0113 = m_C_2*IT_0080;
    const ccomplex_t IT_0114 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0115 = IT_0008*IT_0114;
    const ccomplex_t IT_0116 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0117 = IT_0008*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.5)*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0013*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0123 = IT_0001*IT_0121*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0125 = IT_0008*IT_0124;
    const ccomplex_t IT_0126 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0127 = IT_0008*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0051*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0133 = IT_0043*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0135 = IT_0008*IT_0134;
    const ccomplex_t IT_0136 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0137 = IT_0008*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + (-0.5)*IT_0138);
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0069*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0143 = IT_0078*IT_0141*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0145 = IT_0008*IT_0144;
    const ccomplex_t IT_0146 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0147 = IT_0008*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + (-0.5)*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = IT_0031*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0153 = IT_0040*IT_0151*IT_0152;
    const ccomplex_t IT_0154 = m_C_2*IT_0023;
    const ccomplex_t IT_0155 = m_C_2*IT_0061;
    const ccomplex_t IT_0156 = m_C_2*IT_0042;
    const ccomplex_t IT_0157 = (-0.5)*IT_0113 + 0.5*IT_0123 + 0.5*IT_0133 +
       0.5*IT_0143 + 0.5*IT_0153 + (-0.5)*IT_0154 + (-0.5)*IT_0155 + (-0.5)
      *IT_0156;
    const ccomplex_t IT_0158 = pow(m_W, -2);
    const ccomplex_t IT_0159 = pow(m_W, 2);
    const ccomplex_t IT_0160 = s_34*IT_0159;
    const ccomplex_t IT_0161 = IT_0158*IT_0160;
    const ccomplex_t IT_0162 = (-0.25)*IT_0161;
    const ccomplex_t IT_0163 = s_34 + IT_0162;
    const ccomplex_t IT_0164 = m_C_2*IT_0163;
    const ccomplex_t IT_0165 = (-4)*IT_0164;
    const ccomplex_t IT_0166 = s_13*IT_0159;
    const ccomplex_t IT_0167 = s_14*s_34;
    const ccomplex_t IT_0168 = pow(m_W, 4);
    const ccomplex_t IT_0169 = s_13*IT_0168;
    const ccomplex_t IT_0170 = IT_0158*IT_0169;
    const ccomplex_t IT_0171 = IT_0023 + IT_0042 + IT_0061 + IT_0080 + 0.5
      *IT_0098 + 0.5*IT_0110;
    const ccomplex_t IT_0172 = s_13*s_14;
    const ccomplex_t IT_0173 = s_34*IT_0094;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = s_14*s_34*IT_0158;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = s_13 + IT_0177;
    const ccomplex_t IT_0179 = m_C_2*IT_0178;
    const ccomplex_t IT_0180 = pow(s_14, 2);
    const ccomplex_t IT_0181 = IT_0158*IT_0180;
    const ccomplex_t IT_0182 = -IT_0094;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = s_13*IT_0183;
    const ccomplex_t IT_0185 = s_13*IT_0178;
    const ccomplex_t IT_0186 = pow(s_34, 2);
    const ccomplex_t IT_0187 = s_13*IT_0158*IT_0186;
    const ccomplex_t IT_0188 = m_C_2*IT_0158*IT_0186;
    const ccomplex_t IT_0189 = 2*IT_0111;
    const ccomplex_t IT_0190 = (-2)*IT_0167;
    const ccomplex_t IT_0191 = IT_0166 + IT_0190;
    const ccomplex_t IT_0192 = IT_0158*IT_0191;
    const ccomplex_t IT_0193 = (-2)*conj(IT_0171);
    const ccomplex_t IT_0194 = (-2)*IT_0081*(conj(IT_0111)*IT_0112 + -conj
      (IT_0157)*IT_0165 + 2*conj(IT_0081)*(IT_0166 + (-2)*IT_0167 + (-0.5)
      *IT_0170) + -conj(IT_0171)*IT_0175) + 2*IT_0171*(conj(IT_0081)*IT_0175 + 
      -conj(IT_0157)*IT_0179 + conj(IT_0171)*IT_0184 + -conj(IT_0111)*IT_0185) +
       (conj(IT_0111)*IT_0187 + conj(IT_0157)*IT_0188)*IT_0189 + 2*IT_0157*(conj
      (IT_0081)*IT_0165 + conj(IT_0111)*IT_0188 + 2*conj(IT_0157)*(s_13 + (-0.5)
      *IT_0192) + 0.5*IT_0179*IT_0193) + (-2)*IT_0111*(conj(IT_0081)*IT_0112 + (
      -0.5)*IT_0185*IT_0193);
    return create_ccomplex_return(IT_0194);
}

