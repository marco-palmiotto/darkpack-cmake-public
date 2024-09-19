#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_snu_tau_to_nu_tau_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_snu_tau_to_nu_tau_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t A_tau = param->A_tau;
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
    const ccomplex_t IT_0000 = sin(beta);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cpow(IT_0001, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*(A_tau*IT_0000*conj
      (U_stau_10) + (-1.4142135623731)*e_em*IT_0005*(m_W*IT_0000*IT_0001*conj
      (U_stau_00) + (-0.5)*m_tau*IT_0002*(m_tau*IT_0000*IT_0003*conj(U_stau_00) 
      + mu_h*conj(U_stau_10))));
    const ccomplex_t IT_0007 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0005;
    const ccomplex_t IT_0010 = conj(U_d2)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0012);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow(s_13 + (-0.5)*IT_0015 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(A_tau*IT_0000*conj
      (U_stau_11) + (-1.4142135623731)*e_em*IT_0005*(m_W*IT_0000*IT_0001*conj
      (U_stau_01) + (-0.5)*m_tau*IT_0002*(m_tau*IT_0000*IT_0003*conj(U_stau_01) 
      + mu_h*conj(U_stau_11))));
    const ccomplex_t IT_0020 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = conj(U_d2)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0023 = IT_0009*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + (-0.5)*IT_0024);
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = cpow(s_13 + (-0.5)*IT_0015 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cos(theta_W);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = conj(N_B1)*e_em;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W1)*e_em;
    const ccomplex_t IT_0036 = IT_0005*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0000*IT_0005;
    const ccomplex_t IT_0041 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = IT_0000*IT_0031;
    const ccomplex_t IT_0045 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = IT_0043 + IT_0047;
    const ccomplex_t IT_0049 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0050 = IT_0040*IT_0049;
    const ccomplex_t IT_0051 = (-2)*IT_0050;
    const ccomplex_t IT_0052 = IT_0048 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0039*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0057 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0058 = cpow(s_23 + (-0.5)*IT_0057 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0055*IT_0056*IT_0058;
    const ccomplex_t IT_0060 = conj(N_B3)*e_em;
    const ccomplex_t IT_0061 = IT_0031*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = conj(N_W3)*e_em;
    const ccomplex_t IT_0064 = IT_0005*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + -IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0069 = IT_0040*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0072 = IT_0044*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0075 = IT_0040*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = IT_0073 + IT_0076;
    const ccomplex_t IT_0078 = IT_0070 + IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0067*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0057 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0082*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B4)*e_em;
    const ccomplex_t IT_0086 = IT_0031*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W4)*e_em;
    const ccomplex_t IT_0089 = IT_0005*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + -IT_0090);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0094 = IT_0040*IT_0093;
    const ccomplex_t IT_0095 = (-2)*IT_0094;
    const ccomplex_t IT_0096 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0097 = IT_0044*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0100 = IT_0040*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = IT_0098 + IT_0101;
    const ccomplex_t IT_0103 = IT_0095 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0092*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0108 = cpow(s_23 + (-0.5)*IT_0057 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0106*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B2)*e_em;
    const ccomplex_t IT_0111 = IT_0031*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W2)*e_em;
    const ccomplex_t IT_0114 = IT_0005*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0119 = IT_0040*IT_0118;
    const ccomplex_t IT_0120 = (-2)*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0122 = IT_0044*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0125 = IT_0040*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = IT_0123 + IT_0126;
    const ccomplex_t IT_0128 = IT_0120 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0117*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0133 = cpow(s_23 + (-0.5)*IT_0057 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0131*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = 0.5*IT_0018 + 0.5*IT_0029 + 0.5*IT_0059 + 0.5
      *IT_0084 + 0.5*IT_0109 + 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0001*IT_0005;
    const ccomplex_t IT_0137 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0138;
    const ccomplex_t IT_0140 = IT_0001*IT_0031;
    const ccomplex_t IT_0141 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = IT_0143 + IT_0146;
    const ccomplex_t IT_0148 = IT_0139 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0117*IT_0150;
    const ccomplex_t IT_0152 = IT_0133*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0155 = IT_0136*IT_0154;
    const ccomplex_t IT_0156 = 2*IT_0155;
    const ccomplex_t IT_0157 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0158 = IT_0140*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0161 = IT_0136*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = IT_0159 + IT_0162;
    const ccomplex_t IT_0164 = IT_0156 + IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0067*IT_0166;
    const ccomplex_t IT_0168 = IT_0083*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0171 = IT_0136*IT_0170;
    const ccomplex_t IT_0172 = 2*IT_0171;
    const ccomplex_t IT_0173 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0174 = IT_0140*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0177 = IT_0136*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = IT_0175 + IT_0178;
    const ccomplex_t IT_0180 = IT_0172 + IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0092*IT_0182;
    const ccomplex_t IT_0184 = IT_0108*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0187 = IT_0136*IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0190 = IT_0140*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0193 = IT_0136*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = IT_0191 + IT_0194;
    const ccomplex_t IT_0196 = IT_0188 + IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = IT_0039*IT_0198;
    const ccomplex_t IT_0200 = IT_0058*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0153 + 0.5*IT_0169 + 0.5*IT_0185 + 0.5
      *IT_0201;
    const ccomplex_t IT_0203 = s_23*m_C_2;
    const ccomplex_t IT_0204 = s_12*s_23;
    const ccomplex_t IT_0205 = s_13*IT_0057;
    return create_ccomplex_return(2*IT_0135*(s_13*conj(IT_0135) + conj(IT_0202
      )*IT_0203) + 2*IT_0202*(conj(IT_0135)*IT_0203 + 2*conj(IT_0202)*(IT_0204 +
       (-0.5)*IT_0205)));
}

