#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_snu_tau_to_nu_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_snu_tau_to_nu_tau_tau(
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
    const creal_t m_stau_2 = param->m_stau_2;
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
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0001 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0002 = cpow(s_23 + (-0.5)*IT_0001 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = conj(N_d3)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + IT_0012 + 
      -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B3)*e_em;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*e_em;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = IT_0000*IT_0002*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0009;
    const ccomplex_t IT_0033 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = conj(U_d2)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0036 = IT_0015*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + (-0.5)*IT_0037);
    const ccomplex_t IT_0039 = IT_0032*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0041 = pow(m_stau_2, 2);
    const ccomplex_t IT_0042 = cpow(s_13 + (-0.5)*IT_0041 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0039*IT_0040*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0031 + 0.5*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0009;
    const ccomplex_t IT_0046 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0047 = IT_0009*IT_0046;
    const ccomplex_t IT_0048 = conj(U_d1)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0049 = IT_0015*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + (-0.5)*IT_0050);
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0054 = cpow(s_13 + (-0.5)*IT_0041 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B1)*e_em;
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W1)*e_em;
    const ccomplex_t IT_0060 = IT_0009*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_d1)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0071 = IT_0015*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + IT_0069 + 
      -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0063*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0077 = cpow(s_23 + (-0.5)*IT_0001 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0075*IT_0076*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B2)*e_em;
    const ccomplex_t IT_0080 = IT_0004*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W2)*e_em;
    const ccomplex_t IT_0083 = IT_0009*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + -IT_0084);
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0088 = IT_0004*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0091 = IT_0009*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_d2)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0094 = IT_0015*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + IT_0092 + 
      -IT_0095);
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0086*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0100 = cpow(s_23 + (-0.5)*IT_0001 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0098*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = conj(N_B4)*e_em;
    const ccomplex_t IT_0103 = IT_0004*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W4)*e_em;
    const ccomplex_t IT_0106 = IT_0009*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0111 = IT_0004*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0114 = IT_0009*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_d4)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0117 = IT_0015*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + IT_0115 + 
      -IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0109*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0123 = cpow(s_23 + (-0.5)*IT_0001 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0121*IT_0122*IT_0123;
    const ccomplex_t IT_0125 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0126 = IT_0004*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = N_d1*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0129 = IT_0015*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 0.5*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0063*IT_0132;
    const ccomplex_t IT_0134 = IT_0077*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_tau*IT_0135;
    const ccomplex_t IT_0137 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0138 = IT_0004*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = N_d2*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0141 = IT_0015*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 0.5*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0086*IT_0144;
    const ccomplex_t IT_0146 = IT_0100*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = m_tau*IT_0147;
    const ccomplex_t IT_0149 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0150 = IT_0004*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = N_d3*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0153 = IT_0015*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 0.5*IT_0154);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0029*IT_0156;
    const ccomplex_t IT_0158 = IT_0002*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = m_tau*IT_0159;
    const ccomplex_t IT_0161 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0162 = IT_0004*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = N_d4*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0165 = IT_0015*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + 0.5*IT_0166);
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = IT_0109*IT_0168;
    const ccomplex_t IT_0170 = IT_0123*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_tau*IT_0171;
    const ccomplex_t IT_0173 = (-0.5)*IT_0055 + (-0.5)*IT_0078 + (-0.5)
      *IT_0101 + (-0.5)*IT_0124 + 0.5*IT_0136 + 0.5*IT_0148 + 0.5*IT_0160 + 0.5
      *IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0009*IT_0014*IT_0016;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = IT_0038*IT_0175;
    const ccomplex_t IT_0177 = IT_0042*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0009*IT_0014*IT_0016;
    const ccomplex_t IT_0180 = 0.5*IT_0179;
    const ccomplex_t IT_0181 = IT_0051*IT_0180;
    const ccomplex_t IT_0182 = IT_0054*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0135 + 0.5*IT_0147 + 0.5*IT_0159 + 0.5
      *IT_0171 + (-0.5)*IT_0178 + (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = s_13*m_tau;
    const ccomplex_t IT_0186 = s_13*s_14;
    const ccomplex_t IT_0187 = s_34*IT_0041;
    const ccomplex_t IT_0188 = (-2)*s_34;
    const ccomplex_t IT_0189 = 2*s_34*(IT_0044*conj(IT_0044) + IT_0173*conj
      (IT_0173)) + 2*((conj(IT_0044) + -conj(IT_0173))*IT_0184 + (IT_0044 + 
      -IT_0173)*conj(IT_0184))*IT_0185 + 4*IT_0184*conj(IT_0184)*(IT_0186 + (
      -0.5)*IT_0187) + (conj(IT_0044)*IT_0173 + IT_0044*conj(IT_0173))*IT_0188;
    return create_ccomplex_return(IT_0189);
}

