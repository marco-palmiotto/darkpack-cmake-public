#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_snu_tau_to_tau_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_snu_tau_to_tau_Hp(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0001 = m_N_3*m_tau*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 2*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = IT_0015 + IT_0018;
    const ccomplex_t IT_0020 = IT_0009 + IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = pow(m_tau, 2);
    const ccomplex_t IT_0025 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0024 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0003;
    const ccomplex_t IT_0029 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0030 = IT_0006*IT_0029;
    const ccomplex_t IT_0031 = 2*IT_0030;
    const ccomplex_t IT_0032 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0033 = IT_0012*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = IT_0034 + IT_0037;
    const ccomplex_t IT_0039 = IT_0031 + IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0028*IT_0041;
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0024 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0027 + 0.5*IT_0045;
    const ccomplex_t IT_0047 = pow(m_W, -1);
    const ccomplex_t IT_0048 = cpow(IT_0005, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0003*IT_0047*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = sin(beta);
    const ccomplex_t IT_0052 = IT_0003*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = IT_0011*IT_0051;
    const ccomplex_t IT_0057 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = IT_0055 + IT_0059;
    const ccomplex_t IT_0061 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0062 = IT_0052*IT_0061;
    const ccomplex_t IT_0063 = (-2)*IT_0062;
    const ccomplex_t IT_0064 = IT_0060 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0050*IT_0066;
    const ccomplex_t IT_0068 = IT_0043*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0003*IT_0047*IT_0048;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0073 = IT_0052*IT_0072;
    const ccomplex_t IT_0074 = (-2)*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0076 = IT_0056*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0079 = IT_0052*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = IT_0077 + IT_0080;
    const ccomplex_t IT_0082 = IT_0074 + IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0071*IT_0084;
    const ccomplex_t IT_0086 = IT_0025*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0069 + 0.5*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0003*IT_0047*IT_0048*IT_0051;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B3)*e_em;
    const ccomplex_t IT_0092 = IT_0011*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W3)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0090*IT_0098;
    const ccomplex_t IT_0100 = pow(m_N_3, 2);
    const ccomplex_t IT_0101 = cpow(s_12 + 0.5*IT_0000 + 0.5*IT_0100 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0088 + IT_0104;
    const ccomplex_t IT_0106 = conj(IT_0088) + conj(IT_0104);
    const ccomplex_t IT_0107 = m_N_3*m_tau;
    const ccomplex_t IT_0108 = m_tau*IT_0045;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0110 = IT_0041*IT_0050;
    const ccomplex_t IT_0111 = IT_0043*IT_0109*IT_0110;
    const ccomplex_t IT_0112 = m_tau*IT_0027;
    const ccomplex_t IT_0113 = 0.5*IT_0108 + (-0.5)*IT_0111 + 0.5*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0115 = IT_0022*IT_0071;
    const ccomplex_t IT_0116 = IT_0025*IT_0114*IT_0115;
    const ccomplex_t IT_0117 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0100 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(A_tau*IT_0051*conj
      (U_stau_10) + (-1.4142135623731)*e_em*IT_0003*(m_W*IT_0005*IT_0051*conj
      (U_stau_00) + (-0.5)*m_tau*IT_0047*(m_tau*IT_0048*IT_0051*conj(U_stau_00) 
      + mu_h*conj(U_stau_10))));
    const ccomplex_t IT_0119 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0120 = IT_0011*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = IT_0003*IT_0048;
    const ccomplex_t IT_0123 = N_d3*e_em*m_tau*IT_0047*U_stau_00;
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0125);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0118*IT_0127;
    const ccomplex_t IT_0129 = IT_0117*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0100 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(A_tau*IT_0051*conj
      (U_stau_11) + (-1.4142135623731)*e_em*IT_0003*(m_W*IT_0005*IT_0051*conj
      (U_stau_01) + (-0.5)*m_tau*IT_0047*(m_tau*IT_0048*IT_0051*conj(U_stau_01) 
      + mu_h*conj(U_stau_11))));
    const ccomplex_t IT_0133 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0134 = IT_0011*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = N_d3*e_em*m_tau*IT_0047*U_stau_01;
    const ccomplex_t IT_0137 = IT_0122*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + 0.5*IT_0138);
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0132*IT_0140;
    const ccomplex_t IT_0142 = IT_0131*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = m_N_3*IT_0103;
    const ccomplex_t IT_0145 = 0.5*IT_0116 + (-0.5)*IT_0130 + (-0.5)*IT_0143 +
       (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = IT_0004*IT_0084;
    const ccomplex_t IT_0147 = IT_0025*IT_0114*IT_0146;
    const ccomplex_t IT_0148 = m_tau*IT_0069;
    const ccomplex_t IT_0149 = m_tau*IT_0087;
    const ccomplex_t IT_0150 = (-0.5)*IT_0147 + 0.5*IT_0148 + 0.5*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0152 = IT_0011*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0155 = IT_0003*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = conj(N_d3)*e_em*m_tau*IT_0047*U_stau_10;
    const ccomplex_t IT_0158 = IT_0122*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0153 + IT_0156 + 
      -IT_0159);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0118*IT_0161;
    const ccomplex_t IT_0163 = IT_0117*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0028*IT_0066;
    const ccomplex_t IT_0166 = IT_0043*IT_0109*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0168 = IT_0011*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0171 = IT_0003*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_d3)*e_em*m_tau*IT_0047*U_stau_11;
    const ccomplex_t IT_0174 = IT_0122*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + IT_0172 + 
      -IT_0175);
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0132*IT_0177;
    const ccomplex_t IT_0179 = IT_0131*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0164 + 0.5*IT_0166 + 0.5*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0145;
    const ccomplex_t IT_0183 = (-2)*IT_0113;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = 2*conj(IT_0145);
    const ccomplex_t IT_0186 = (-2)*conj(IT_0113);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = s_23*m_N_3;
    const ccomplex_t IT_0189 = 2*IT_0181;
    const ccomplex_t IT_0190 = (-2)*IT_0150;
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 2*conj(IT_0181);
    const ccomplex_t IT_0193 = (-2)*conj(IT_0150);
    const ccomplex_t IT_0194 = IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = s_12*m_tau;
    const ccomplex_t IT_0196 = 2*s_13;
    const ccomplex_t IT_0197 = s_12*s_23;
    const ccomplex_t IT_0198 = s_13*IT_0000;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = 4*IT_0088;
    const ccomplex_t IT_0202 = (-2)*s_13;
    const ccomplex_t IT_0203 = 2*IT_0001*(conj(IT_0046)*IT_0105 + IT_0046
      *IT_0106) + 2*IT_0107*((conj(IT_0113) + -conj(IT_0145))*IT_0150 + (IT_0113
       + -IT_0145)*conj(IT_0150) + 0.5*conj(IT_0181)*IT_0184 + 0.5*IT_0181
      *IT_0187) + IT_0188*(IT_0106*IT_0184 + IT_0105*IT_0187 + conj(IT_0046)
      *IT_0191 + IT_0046*IT_0194) + (conj(IT_0046)*IT_0184 + IT_0046*IT_0187 +
       IT_0106*IT_0191 + IT_0105*IT_0194)*IT_0195 + (IT_0113*conj(IT_0113) +
       IT_0145*conj(IT_0145) + IT_0150*conj(IT_0150) + IT_0181*conj(IT_0181))
      *IT_0196 + 4*IT_0200*(IT_0046*conj(IT_0046) + IT_0106*(IT_0104 + 0.25
      *IT_0201)) + (conj(IT_0113)*IT_0145 + IT_0113*conj(IT_0145) + conj(IT_0150
      )*IT_0181 + IT_0150*conj(IT_0181))*IT_0202;
    return create_ccomplex_return(IT_0203);
}

