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
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0001 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = conj(U_d2)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0006 + (-0.5)*IT_0013);
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = conj(N_d3)*e_em*m_tau*IT_0010*conj(U_stau_01);
    const ccomplex_t IT_0021 = IT_0009*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + 0.5*IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = IT_0014*IT_0024;
    const ccomplex_t IT_0026 = IT_0002*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow(s_13 + (-0.5)*IT_0001 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0029 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = conj(U_d2)*e_em*m_tau*IT_0010*U_stau_10;
    const ccomplex_t IT_0032 = IT_0009*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-0.5)*IT_0033);
    const ccomplex_t IT_0035 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0036 = IT_0016*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_d3)*e_em*m_tau*IT_0010*conj(U_stau_00);
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + 0.5*IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = IT_0028*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0027 + (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0048 = IT_0016*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = N_d3*e_em*m_tau*IT_0010*conj(U_stau_11);
    const ccomplex_t IT_0054 = IT_0009*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0049 + IT_0052 + 
      -IT_0055);
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0014*IT_0057;
    const ccomplex_t IT_0059 = IT_0002*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0062 = IT_0016*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_d3*e_em*m_tau*IT_0010*conj(U_stau_10);
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0063 + IT_0066 + 
      -IT_0069);
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0034*IT_0071;
    const ccomplex_t IT_0073 = IT_0028*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0060 + 0.5*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0079 = IT_0004*IT_0078;
    const ccomplex_t IT_0080 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + 0.5*IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = IT_0077*IT_0084;
    const ccomplex_t IT_0086 = pow(m_N_3, 2);
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0086 + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = 4*IT_0075;
    const ccomplex_t IT_0091 = 0.125*IT_0090;
    const ccomplex_t IT_0092 = IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = s_12*s_34;
    const ccomplex_t IT_0094 = pow(m_W, 2);
    const ccomplex_t IT_0095 = cpow((-2)*s_12 + -IT_0001 + -IT_0086 + IT_0094 
      + -reg_prop, -1);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0004*IT_0010;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = sin(beta);
    const ccomplex_t IT_0099 = IT_0004*IT_0098;
    const ccomplex_t IT_0100 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = 2*IT_0101;
    const ccomplex_t IT_0103 = IT_0016*IT_0098;
    const ccomplex_t IT_0104 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0108 = IT_0099*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = IT_0106 + IT_0109;
    const ccomplex_t IT_0111 = IT_0102 + IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0097*IT_0113;
    const ccomplex_t IT_0115 = IT_0095*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0004*IT_0008*IT_0010*IT_0098;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0004*IT_0007;
    const ccomplex_t IT_0120 = IT_0100*IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0120;
    const ccomplex_t IT_0122 = IT_0007*IT_0016;
    const ccomplex_t IT_0123 = IT_0104*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = IT_0107*IT_0119;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = IT_0124 + IT_0126;
    const ccomplex_t IT_0128 = IT_0121 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0118*IT_0130;
    const ccomplex_t IT_0132 = cpow((-2)*s_12 + -IT_0001 + -IT_0086 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0116 + IT_0134;
    const ccomplex_t IT_0136 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0137 = IT_0119*IT_0136;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0140 = IT_0122*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0143 = IT_0119*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0141 + IT_0144;
    const ccomplex_t IT_0146 = IT_0138 + IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0097*IT_0148;
    const ccomplex_t IT_0150 = IT_0095*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0099*IT_0136;
    const ccomplex_t IT_0153 = (-2)*IT_0152;
    const ccomplex_t IT_0154 = IT_0103*IT_0139;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = IT_0099*IT_0142;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = IT_0155 + IT_0157;
    const ccomplex_t IT_0159 = IT_0153 + IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = IT_0118*IT_0161;
    const ccomplex_t IT_0163 = IT_0132*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0151 + -IT_0164;
    const ccomplex_t IT_0166 = 4*IT_0165;
    const ccomplex_t IT_0167 = s_14*s_23;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0004;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = conj(N_B3)*e_em;
    const ccomplex_t IT_0171 = IT_0016*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W3)*e_em;
    const ccomplex_t IT_0174 = IT_0004*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0172 + -IT_0175);
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0169*IT_0177;
    const ccomplex_t IT_0179 = cpow(s_23 + (-0.5)*IT_0086 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0180 = IT_0178*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0184 = IT_0004*IT_0183;
    const ccomplex_t IT_0185 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0186 = IT_0004*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + (-0.5)*IT_0187);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0077*IT_0189;
    const ccomplex_t IT_0191 = IT_0087*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 4*IT_0182;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0004*IT_0008*IT_0010;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = IT_0177*IT_0195;
    const ccomplex_t IT_0197 = IT_0179*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = 4*IT_0199;
    const ccomplex_t IT_0201 = -IT_0000;
    const ccomplex_t IT_0202 = (-2)*IT_0165;
    const ccomplex_t IT_0203 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0204 = s_13*m_N_3*m_tau;
    const ccomplex_t IT_0205 = 2*IT_0165;
    const ccomplex_t IT_0206 = 2*conj(IT_0165);
    const ccomplex_t IT_0207 = 4*conj(IT_0075);
    const ccomplex_t IT_0208 = 0.5*conj(IT_0046);
    const ccomplex_t IT_0209 = 2*conj(IT_0135);
    const ccomplex_t IT_0210 = 0.5*IT_0182;
    const ccomplex_t IT_0211 = 0.5*conj(IT_0199);
    const ccomplex_t IT_0212 = conj(IT_0135)*IT_0192 + 0.5*IT_0135*(conj
      (IT_0182) + 2*conj(IT_0192)) + 2*conj(IT_0089)*(IT_0046 + 0.125*IT_0166 +
       0.125*IT_0200) + 0.5*conj(IT_0075)*(IT_0199 + 0.5*IT_0205) + 0.5*IT_0075*
      (conj(IT_0199) + 0.5*IT_0206) + 0.25*IT_0046*IT_0207 + 0.5*IT_0090*IT_0208
       + 0.5*IT_0209*IT_0210 + 2*IT_0089*(conj(IT_0046) + 0.5*conj(IT_0165) +
       IT_0211);
    const ccomplex_t IT_0213 = s_34*m_C_2*m_N_3;
    const ccomplex_t IT_0214 = 2*IT_0182;
    const ccomplex_t IT_0215 = (-4)*conj(IT_0135);
    const ccomplex_t IT_0216 = (-2)*conj(IT_0135);
    const ccomplex_t IT_0217 = 0.5*IT_0199;
    const ccomplex_t IT_0218 = 4*conj(IT_0182);
    const ccomplex_t IT_0219 = 0.5*IT_0089;
    const ccomplex_t IT_0220 = IT_0075*conj(IT_0182) + 4*IT_0092*conj(IT_0192)
       + 0.5*conj(IT_0089)*IT_0193 + -IT_0135*(conj(IT_0046) + 2*conj(IT_0165) +
       -conj(IT_0199)) + 4*IT_0192*(conj(IT_0089) + 0.125*IT_0207) + 0.5*conj
      (IT_0075)*IT_0214 + 0.5*IT_0165*IT_0215 + 0.5*IT_0046*IT_0216 + IT_0209
      *IT_0217 + IT_0218*IT_0219;
    const ccomplex_t IT_0221 = s_23*m_C_2*m_tau;
    const ccomplex_t IT_0222 = (-0.5)*(conj(IT_0075) + 2*conj(IT_0089))
      *IT_0135 + conj(IT_0192)*(IT_0046 + -IT_0165 + 2*IT_0199) + IT_0192*(conj
      (IT_0046) + -conj(IT_0165) + 2*conj(IT_0199)) + 0.5*conj(IT_0182)*(IT_0046
       + 0.5*IT_0202) + 0.5*IT_0203*IT_0210 + 0.5*IT_0193*IT_0211 + 0.5*IT_0208
      *IT_0214 + 0.25*IT_0075*IT_0216 + 0.5*IT_0217*IT_0218 + 0.5*IT_0215*IT_0219;
    const ccomplex_t IT_0223 = -IT_0167;
    const ccomplex_t IT_0224 = IT_0093 + IT_0223;
    const ccomplex_t IT_0225 = IT_0201 + IT_0224;
    const ccomplex_t IT_0226 = IT_0000 + IT_0224;
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = 4*IT_0000*(IT_0046*conj(IT_0046) + 2*(IT_0075 +
       2*IT_0089)*conj(IT_0089) + 2*conj(IT_0075)*IT_0092) + 4*IT_0093*(IT_0135
      *conj(IT_0135) + 0.25*conj(IT_0165)*IT_0166) + 8*IT_0167*((IT_0182 + 2
      *IT_0192)*conj(IT_0192) + conj(IT_0182)*(IT_0192 + 0.125*IT_0193) + 0.125
      *conj(IT_0199)*IT_0200) + (IT_0093 + IT_0167 + IT_0201)*(conj(IT_0199)
      *IT_0202 + IT_0199*IT_0203) + 4*IT_0204*IT_0212 + 2*IT_0213*IT_0220 + 4
      *IT_0221*IT_0222 + (-2)*conj(IT_0046)*(IT_0199*IT_0225 + IT_0205*IT_0227) 
      + (-2)*IT_0046*(conj(IT_0199)*IT_0225 + IT_0206*IT_0227);
    return create_ccomplex_return(IT_0228);
}

