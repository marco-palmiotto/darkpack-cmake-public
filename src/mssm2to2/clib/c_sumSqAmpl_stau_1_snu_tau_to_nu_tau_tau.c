#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_snu_tau_to_nu_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_snu_tau_to_nu_tau_tau(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
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
    const ccomplex_t IT_0012 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = conj(N_d1)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0014 + IT_0017 + 
      -IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0029 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0030 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0027*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B4)*e_em;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W4)*e_em;
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_d4)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0047 = IT_0020*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0042 + IT_0045 + 
      -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0039*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0053 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0006;
    const ccomplex_t IT_0056 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0057 = IT_0006*IT_0056;
    const ccomplex_t IT_0058 = conj(U_d1)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0059 = IT_0020*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-0.5)*IT_0060);
    const ccomplex_t IT_0062 = IT_0055*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0064 = pow(m_stau_1, 2);
    const ccomplex_t IT_0065 = cpow(s_13 + (-0.5)*IT_0064 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0062*IT_0063*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0006;
    const ccomplex_t IT_0068 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0069 = IT_0006*IT_0068;
    const ccomplex_t IT_0070 = conj(U_d2)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0071 = IT_0020*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.5)*IT_0072);
    const ccomplex_t IT_0074 = IT_0067*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0076 = cpow(s_13 + (-0.5)*IT_0064 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0074*IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0079 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0080 = conj(N_B2)*e_em;
    const ccomplex_t IT_0081 = IT_0001*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W2)*e_em;
    const ccomplex_t IT_0084 = IT_0006*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + -IT_0085);
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0089 = IT_0001*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0092 = IT_0006*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_d2)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0095 = IT_0020*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0090 + IT_0093 + 
      -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0087*IT_0098;
    const ccomplex_t IT_0100 = IT_0078*IT_0079*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0102 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0103 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0104 = IT_0001*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0107 = IT_0006*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = conj(N_d3)*e_em*m_tau*IT_0021*U_stau_10;
    const ccomplex_t IT_0110 = IT_0020*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0105 + IT_0108 + 
      -IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = conj(N_B3)*e_em;
    const ccomplex_t IT_0115 = IT_0001*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W3)*e_em;
    const ccomplex_t IT_0118 = IT_0006*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + -IT_0119);
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0113*IT_0121;
    const ccomplex_t IT_0123 = IT_0101*IT_0102*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0031 + 0.5*IT_0054 + 0.5*IT_0066 + 0.5
      *IT_0077 + 0.5*IT_0100 + 0.5*IT_0123;
    const ccomplex_t IT_0125 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = N_d4*e_em*m_tau*IT_0021*U_stau_00;
    const ccomplex_t IT_0129 = IT_0020*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 0.5*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0039*IT_0132;
    const ccomplex_t IT_0134 = IT_0053*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0137 = IT_0001*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_d2*e_em*m_tau*IT_0021*U_stau_00;
    const ccomplex_t IT_0140 = IT_0020*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 0.5*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0087*IT_0143;
    const ccomplex_t IT_0145 = IT_0079*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0148 = IT_0001*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = N_d3*e_em*m_tau*IT_0021*U_stau_00;
    const ccomplex_t IT_0151 = IT_0020*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + 0.5*IT_0152);
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0121*IT_0154;
    const ccomplex_t IT_0156 = IT_0102*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0159 = IT_0001*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = N_d1*e_em*m_tau*IT_0021*U_stau_00;
    const ccomplex_t IT_0162 = IT_0020*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + 0.5*IT_0163);
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0011*IT_0165;
    const ccomplex_t IT_0167 = IT_0030*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0135 + IT_0146 + IT_0157 + IT_0168;
    const ccomplex_t IT_0170 = m_tau*IT_0169;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0006*IT_0019*IT_0021;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = IT_0061*IT_0173;
    const ccomplex_t IT_0175 = IT_0065*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0006*IT_0019*IT_0021;
    const ccomplex_t IT_0178 = 0.5*IT_0177;
    const ccomplex_t IT_0179 = IT_0073*IT_0178;
    const ccomplex_t IT_0180 = IT_0076*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = 0.5*IT_0135 + 0.5*IT_0146 + 0.5*IT_0157 + 0.5
      *IT_0168 + (-0.5)*IT_0176 + (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = s_13*m_tau;
    const ccomplex_t IT_0184 = s_13*s_14;
    const ccomplex_t IT_0185 = s_34*IT_0064;
    const ccomplex_t IT_0186 = (-2)*s_34;
    const ccomplex_t IT_0187 = 2*s_34*(IT_0124*conj(IT_0124) + IT_0171*conj
      (IT_0171)) + 2*((conj(IT_0124) + -conj(IT_0171))*IT_0182 + (IT_0124 + 
      -IT_0171)*conj(IT_0182))*IT_0183 + 4*IT_0182*conj(IT_0182)*(IT_0184 + (
      -0.5)*IT_0185) + (conj(IT_0124)*IT_0171 + IT_0124*conj(IT_0171))*IT_0186;
    return create_ccomplex_return(IT_0187);
}

