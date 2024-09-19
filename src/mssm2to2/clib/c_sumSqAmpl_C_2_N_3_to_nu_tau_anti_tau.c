#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_nu_tau_anti_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 1.4142135623731*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = N_d3*e_em*m_tau*IT_0014*conj(U_stau_11);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0005 + IT_0010 + 
      -IT_0017);
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = conj(U_d2)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + (-0.5)*IT_0024);
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = pow(m_C_2, 2);
    const ccomplex_t IT_0028 = cpow(s_13 + (-0.5)*IT_0027 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_d3*e_em*m_tau*IT_0014*conj(U_stau_10);
    const ccomplex_t IT_0038 = IT_0013*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + IT_0036 + 
      -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = conj(U_d2)*e_em*m_tau*IT_0014*U_stau_10;
    const ccomplex_t IT_0045 = IT_0013*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + (-0.5)*IT_0046);
    const ccomplex_t IT_0048 = IT_0041*IT_0047;
    const ccomplex_t IT_0049 = cpow(s_13 + (-0.5)*IT_0027 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0030 + 0.5*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0056 = IT_0007*IT_0055;
    const ccomplex_t IT_0057 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0058 = IT_0007*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + 0.5*IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0054*IT_0061;
    const ccomplex_t IT_0063 = pow(m_N_3, 2);
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0027 + -IT_0063 + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 4*IT_0052;
    const ccomplex_t IT_0068 = 0.125*IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0071 = IT_0002*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_d3)*e_em*m_tau*IT_0014*conj(U_stau_00);
    const ccomplex_t IT_0074 = IT_0013*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + 0.5*IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0047*IT_0077;
    const ccomplex_t IT_0079 = IT_0049*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0082 = IT_0002*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_d3)*e_em*m_tau*IT_0014*conj(U_stau_01);
    const ccomplex_t IT_0085 = IT_0013*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + 0.5*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0025*IT_0088;
    const ccomplex_t IT_0090 = IT_0028*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0080 + (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = 4*IT_0092;
    const ccomplex_t IT_0094 = s_14*s_23;
    const ccomplex_t IT_0095 = cpow(s_23 + (-0.5)*IT_0063 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0096 = conj(N_B3)*e_em;
    const ccomplex_t IT_0097 = IT_0002*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W3)*e_em;
    const ccomplex_t IT_0100 = IT_0007*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0007;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0113 = IT_0007*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + (-0.5)*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0054*IT_0116;
    const ccomplex_t IT_0118 = IT_0064*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = 4*IT_0109;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0007*IT_0012*IT_0014;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0103*IT_0122;
    const ccomplex_t IT_0124 = IT_0095*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = 4*IT_0126;
    const ccomplex_t IT_0128 = s_12*s_34;
    const ccomplex_t IT_0129 = pow(m_W, 2);
    const ccomplex_t IT_0130 = cpow((-2)*s_12 + -IT_0027 + -IT_0063 + IT_0129 
      + -reg_prop, -1);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0007*IT_0014;
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = sin(beta);
    const ccomplex_t IT_0134 = IT_0007*IT_0133;
    const ccomplex_t IT_0135 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = 2*IT_0136;
    const ccomplex_t IT_0138 = IT_0002*IT_0133;
    const ccomplex_t IT_0139 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0143 = IT_0134*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0141 + IT_0144;
    const ccomplex_t IT_0146 = IT_0137 + IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0132*IT_0148;
    const ccomplex_t IT_0150 = IT_0130*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0007*IT_0012*IT_0014*IT_0133;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = IT_0007*IT_0011;
    const ccomplex_t IT_0155 = IT_0135*IT_0154;
    const ccomplex_t IT_0156 = 2*IT_0155;
    const ccomplex_t IT_0157 = IT_0002*IT_0011;
    const ccomplex_t IT_0158 = IT_0139*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = IT_0142*IT_0154;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = IT_0156 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0153*IT_0165;
    const ccomplex_t IT_0167 = cpow((-2)*s_12 + -IT_0027 + -IT_0063 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = IT_0151 + IT_0169;
    const ccomplex_t IT_0171 = 4*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0173 = IT_0134*IT_0172;
    const ccomplex_t IT_0174 = (-2)*IT_0173;
    const ccomplex_t IT_0175 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0176 = IT_0138*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0179 = IT_0134*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = IT_0177 + IT_0180;
    const ccomplex_t IT_0182 = IT_0174 + IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0153*IT_0184;
    const ccomplex_t IT_0186 = IT_0167*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0154*IT_0172;
    const ccomplex_t IT_0189 = (-2)*IT_0188;
    const ccomplex_t IT_0190 = IT_0157*IT_0175;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = IT_0154*IT_0178;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0189 + IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0132*IT_0197;
    const ccomplex_t IT_0199 = IT_0130*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = -IT_0187 + -IT_0200;
    const ccomplex_t IT_0202 = 4*IT_0201;
    const ccomplex_t IT_0203 = -IT_0000;
    const ccomplex_t IT_0204 = (-2)*IT_0201;
    const ccomplex_t IT_0205 = (-2)*conj(IT_0201);
    const ccomplex_t IT_0206 = s_13*m_N_3*m_tau;
    const ccomplex_t IT_0207 = 2*IT_0170;
    const ccomplex_t IT_0208 = 2*conj(IT_0170);
    const ccomplex_t IT_0209 = 2*IT_0201;
    const ccomplex_t IT_0210 = 2*conj(IT_0201);
    const ccomplex_t IT_0211 = 4*conj(IT_0052);
    const ccomplex_t IT_0212 = 0.5*conj(IT_0126);
    const ccomplex_t IT_0213 = 0.25*IT_0067*conj(IT_0092) + IT_0119*conj
      (IT_0170) + 0.25*conj(IT_0119)*IT_0171 + 2*conj(IT_0066)*(IT_0092 + 0.125
      *IT_0127 + 0.125*IT_0202) + 0.25*conj(IT_0109)*IT_0207 + 0.25*IT_0109
      *IT_0208 + 0.5*conj(IT_0052)*(IT_0126 + 0.5*IT_0209) + 0.5*IT_0052*(conj
      (IT_0126) + 0.5*IT_0210) + 0.25*IT_0092*IT_0211 + 2*IT_0066*(conj(IT_0092)
       + 0.5*conj(IT_0201) + IT_0212);
    const ccomplex_t IT_0214 = s_23*m_C_2*m_tau;
    const ccomplex_t IT_0215 = 4*conj(IT_0109);
    const ccomplex_t IT_0216 = (-4)*conj(IT_0170);
    const ccomplex_t IT_0217 = (-2)*conj(IT_0170);
    const ccomplex_t IT_0218 = (-0.5)*IT_0201;
    const ccomplex_t IT_0219 = (-0.5)*(conj(IT_0052) + 2*conj(IT_0066))
      *IT_0170 + IT_0119*(conj(IT_0092) + 2*conj(IT_0126) + -conj(IT_0201)) +
       0.5*conj(IT_0109)*(IT_0092 + 0.5*IT_0204) + 0.5*IT_0109*(conj(IT_0092) +
       0.5*IT_0205) + 0.5*IT_0120*IT_0212 + 0.25*IT_0126*IT_0215 + 0.25*IT_0066
      *IT_0216 + 0.25*IT_0052*IT_0217 + 2*conj(IT_0119)*(0.125*IT_0093 + IT_0126
       + IT_0218);
    const ccomplex_t IT_0220 = s_34*m_C_2*m_N_3;
    const ccomplex_t IT_0221 = -conj(IT_0052)*IT_0109 + -IT_0052*conj(IT_0109)
       + (-4)*IT_0069*conj(IT_0119) + (-0.5)*conj(IT_0066)*IT_0120 + IT_0170*
      (conj(IT_0092) + 2*conj(IT_0201)) + (-0.5)*conj(IT_0126)*IT_0207 + (-0.5)
      *IT_0126*IT_0208 + (-4)*IT_0119*(conj(IT_0066) + 0.125*IT_0211) + (-0.5)
      *IT_0066*IT_0215 + (-0.5)*IT_0092*IT_0217 + IT_0216*IT_0218;
    const ccomplex_t IT_0222 = -IT_0094;
    const ccomplex_t IT_0223 = IT_0128 + IT_0222;
    const ccomplex_t IT_0224 = IT_0203 + IT_0223;
    const ccomplex_t IT_0225 = IT_0000 + IT_0223;
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = 8*IT_0000*((IT_0052 + 2*IT_0066)*conj(IT_0066) 
      + conj(IT_0052)*IT_0069 + 0.125*conj(IT_0092)*IT_0093) + 8*IT_0094*(
      (IT_0109 + 2*IT_0119)*conj(IT_0119) + conj(IT_0109)*(IT_0119 + 0.125
      *IT_0120) + 0.125*conj(IT_0126)*IT_0127) + IT_0128*(conj(IT_0170)*IT_0171 
      + conj(IT_0201)*IT_0202) + (IT_0094 + IT_0128 + IT_0203)*(conj(IT_0126)
      *IT_0204 + IT_0126*IT_0205) + 4*IT_0206*IT_0213 + 4*IT_0214*IT_0219 + (-2)
      *IT_0220*IT_0221 + (-2)*conj(IT_0092)*(IT_0126*IT_0224 + IT_0209*IT_0226) 
      + (-2)*IT_0092*(conj(IT_0126)*IT_0224 + IT_0210*IT_0226);
    return create_ccomplex_return(IT_0227);
}

