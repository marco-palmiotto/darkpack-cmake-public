#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_snu_tau_to_nu_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_snu_tau_to_nu_tau_W(
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
    const creal_t m_C_1 = param->m_C_1;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + (-0.5)*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0021 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B3)*e_em;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W3)*e_em;
    const ccomplex_t IT_0028 = IT_0006*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0035 = IT_0006*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + (-0.5)*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0041 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0040*IT_0041;
    const ccomplex_t IT_0043 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0046 = IT_0006*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + 0.5*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0011*IT_0049;
    const ccomplex_t IT_0051 = IT_0022*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = m_C_1*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B1)*e_em;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W1)*e_em;
    const ccomplex_t IT_0058 = IT_0006*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + -IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0063 = IT_0006*IT_0062;
    const ccomplex_t IT_0064 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0065 = IT_0006*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + (-0.5)*IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0061*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0071 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0069*IT_0070*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B4)*e_em;
    const ccomplex_t IT_0074 = IT_0001*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = conj(N_W4)*e_em;
    const ccomplex_t IT_0077 = IT_0006*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + -IT_0078);
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0082 = IT_0006*IT_0081;
    const ccomplex_t IT_0083 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0084 = IT_0006*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-0.5)*IT_0085);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0090 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0088*IT_0089*IT_0090;
    const ccomplex_t IT_0092 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0093 = IT_0006*IT_0092;
    const ccomplex_t IT_0094 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0095 = IT_0006*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + 0.5*IT_0096);
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0061*IT_0098;
    const ccomplex_t IT_0100 = IT_0071*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_C_1*IT_0101;
    const ccomplex_t IT_0103 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0104 = IT_0006*IT_0103;
    const ccomplex_t IT_0105 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0106 = IT_0006*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + 0.5*IT_0107);
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = IT_0031*IT_0109;
    const ccomplex_t IT_0111 = IT_0041*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = m_C_1*IT_0112;
    const ccomplex_t IT_0114 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0115 = IT_0006*IT_0114;
    const ccomplex_t IT_0116 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0117 = IT_0006*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + 0.5*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0080*IT_0120;
    const ccomplex_t IT_0122 = IT_0090*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_C_1*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0023 + 0.5*IT_0042 + (-0.5)*IT_0053 +
       0.5*IT_0072 + 0.5*IT_0091 + (-0.5)*IT_0102 + (-0.5)*IT_0113 + (-0.5)
      *IT_0124;
    const ccomplex_t IT_0126 = pow(m_W, -2);
    const ccomplex_t IT_0127 = s_14*s_34*IT_0126;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_13 + IT_0128;
    const ccomplex_t IT_0130 = m_C_1*IT_0129;
    const ccomplex_t IT_0131 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0132 = IT_0006*IT_0131;
    const ccomplex_t IT_0133 = cos(beta);
    const ccomplex_t IT_0134 = cpow(IT_0133, -1);
    const ccomplex_t IT_0135 = IT_0006*IT_0134;
    const ccomplex_t IT_0136 = pow(m_W, -1);
    const ccomplex_t IT_0137 = conj(U_d1)*e_em*m_tau*IT_0136*U_stau_11;
    const ccomplex_t IT_0138 = IT_0135*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0132 + (-0.5)*IT_0139);
    const ccomplex_t IT_0141 = pow(m_C_1, 2);
    const ccomplex_t IT_0142 = cpow(s_13 + (-0.5)*IT_0141 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0143 = IT_0140*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0006
      *conj(U_stau_01);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0149 = IT_0006*IT_0148;
    const ccomplex_t IT_0150 = conj(U_d1)*e_em*m_tau*IT_0136*U_stau_10;
    const ccomplex_t IT_0151 = IT_0135*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + (-0.5)*IT_0152);
    const ccomplex_t IT_0154 = cpow(s_13 + (-0.5)*IT_0141 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0006
      *conj(U_stau_00);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0156*IT_0158;
    const ccomplex_t IT_0160 = IT_0052 + IT_0101 + IT_0112 + IT_0123 + 0.5
      *IT_0147 + 0.5*IT_0159;
    const ccomplex_t IT_0161 = pow(m_W, 2);
    const ccomplex_t IT_0162 = s_13*IT_0161;
    const ccomplex_t IT_0163 = s_14*s_34;
    const ccomplex_t IT_0164 = (-2)*IT_0163;
    const ccomplex_t IT_0165 = IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = IT_0126*IT_0165;
    const ccomplex_t IT_0167 = (-0.5)*IT_0147 + (-0.5)*IT_0159;
    const ccomplex_t IT_0168 = pow(s_34, 2);
    const ccomplex_t IT_0169 = m_C_1*IT_0126*IT_0168;
    const ccomplex_t IT_0170 = 2*IT_0169;
    const ccomplex_t IT_0171 = (-0.5)*IT_0052 + (-0.5)*IT_0101 + (-0.5)
      *IT_0112 + (-0.5)*IT_0123;
    const ccomplex_t IT_0172 = s_34*IT_0161;
    const ccomplex_t IT_0173 = IT_0126*IT_0172;
    const ccomplex_t IT_0174 = (-0.25)*IT_0173;
    const ccomplex_t IT_0175 = s_34 + IT_0174;
    const ccomplex_t IT_0176 = m_C_1*IT_0175;
    const ccomplex_t IT_0177 = (-4)*IT_0176;
    const ccomplex_t IT_0178 = 2*IT_0177;
    const ccomplex_t IT_0179 = pow(m_W, 4);
    const ccomplex_t IT_0180 = s_13*IT_0179;
    const ccomplex_t IT_0181 = IT_0126*IT_0180;
    const ccomplex_t IT_0182 = s_13*s_14;
    const ccomplex_t IT_0183 = s_34*IT_0141;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = s_13*s_34;
    const ccomplex_t IT_0187 = (-2)*IT_0186;
    const ccomplex_t IT_0188 = s_13*IT_0126*IT_0168;
    const ccomplex_t IT_0189 = s_13*IT_0129;
    const ccomplex_t IT_0190 = (-2)*IT_0189;
    const ccomplex_t IT_0191 = pow(s_14, 2);
    const ccomplex_t IT_0192 = IT_0126*IT_0191;
    const ccomplex_t IT_0193 = -IT_0141;
    const ccomplex_t IT_0194 = IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = s_13*IT_0194;
    const ccomplex_t IT_0196 = (-2)*IT_0125*(IT_0130*conj(IT_0160) + (-2)*conj
      (IT_0125)*(s_13 + (-0.5)*IT_0166) + (-0.5)*conj(IT_0167)*IT_0170 + (-0.5)
      *conj(IT_0171)*IT_0178) + (-4)*IT_0171*((-0.25)*conj(IT_0125)*IT_0178 +
       conj(IT_0171)*(IT_0162 + (-2)*IT_0163 + (-0.5)*IT_0181) + (-0.5)*conj
      (IT_0160)*IT_0185 + (-0.25)*conj(IT_0167)*IT_0187) + 2*IT_0167*(0.5*conj
      (IT_0125)*IT_0170 + 0.5*conj(IT_0171)*IT_0187 + conj(IT_0167)*IT_0188 +
       0.5*conj(IT_0160)*IT_0190) + (-2)*IT_0160*(conj(IT_0125)*IT_0130 + -conj
      (IT_0171)*IT_0185 + (-0.5)*conj(IT_0167)*IT_0190 + -conj(IT_0160)*IT_0195);
    return create_ccomplex_return(IT_0196);
}

