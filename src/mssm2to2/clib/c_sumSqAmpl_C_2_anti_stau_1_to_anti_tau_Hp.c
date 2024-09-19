#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_stau_1_to_anti_tau_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_stau_1_to_anti_tau_Hp(
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
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
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
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0022 = IT_0009*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0001, -1);
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = N_d3*e_em*m_tau*IT_0029*conj(U_stau_10);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0023 + IT_0026 + 
      -IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0020*IT_0034;
    const ccomplex_t IT_0036 = pow(m_tau, 2);
    const ccomplex_t IT_0037 = pow(m_stau_1, 2);
    const ccomplex_t IT_0038 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0000*IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0003;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(A_tau*IT_0042*conj
      (U_stau_10) + (-1.4142135623731)*e_em*IT_0003*(m_W*IT_0001*IT_0042*conj
      (U_stau_00) + (-0.5)*m_tau*IT_0029*(m_tau*IT_0027*IT_0042*conj(U_stau_00) 
      + mu_h*conj(U_stau_10))));
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = pow(m_C_2, 2);
    const ccomplex_t IT_0046 = cpow(s_13 + (-0.5)*IT_0036 + (-0.5)*IT_0045 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0050 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 2*IT_0051;
    const ccomplex_t IT_0053 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0054 = IT_0010*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = IT_0055 + IT_0058;
    const ccomplex_t IT_0060 = IT_0052 + IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0064 = IT_0009*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0067 = IT_0003*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = N_d1*e_em*m_tau*IT_0029*conj(U_stau_10);
    const ccomplex_t IT_0070 = IT_0028*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + IT_0068 + 
      -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0062*IT_0073;
    const ccomplex_t IT_0075 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0049*IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0078 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0079 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0080 = IT_0009*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = N_d2*e_em*m_tau*IT_0029*conj(U_stau_10);
    const ccomplex_t IT_0086 = IT_0028*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0081 + IT_0084 + 
      -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = 2*IT_0091;
    const ccomplex_t IT_0093 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0094 = IT_0010*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0097 = IT_0004*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = IT_0095 + IT_0098;
    const ccomplex_t IT_0100 = IT_0092 + IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0089*IT_0102;
    const ccomplex_t IT_0104 = IT_0077*IT_0078*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0106 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0107 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0108 = IT_0009*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = N_d4*e_em*m_tau*IT_0029*conj(U_stau_10);
    const ccomplex_t IT_0114 = IT_0028*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + IT_0112 + 
      -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0119 = IT_0004*IT_0118;
    const ccomplex_t IT_0120 = 2*IT_0119;
    const ccomplex_t IT_0121 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0122 = IT_0010*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0125 = IT_0004*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = IT_0123 + IT_0126;
    const ccomplex_t IT_0128 = IT_0120 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0117*IT_0130;
    const ccomplex_t IT_0132 = IT_0105*IT_0106*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0134 = IT_0009*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d1)*e_em*m_tau*IT_0029*conj(U_stau_00);
    const ccomplex_t IT_0137 = IT_0028*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + 0.5*IT_0138);
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0062*IT_0140;
    const ccomplex_t IT_0142 = IT_0075*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = m_tau*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0146 = IT_0009*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_d2)*e_em*m_tau*IT_0029*conj(U_stau_00);
    const ccomplex_t IT_0149 = IT_0028*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 0.5*IT_0150);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = IT_0102*IT_0152;
    const ccomplex_t IT_0154 = IT_0078*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_tau*IT_0155;
    const ccomplex_t IT_0157 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0158 = IT_0009*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = conj(N_d3)*e_em*m_tau*IT_0029*conj(U_stau_00);
    const ccomplex_t IT_0161 = IT_0028*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + 0.5*IT_0162);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0020*IT_0164;
    const ccomplex_t IT_0166 = IT_0038*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = m_tau*IT_0167;
    const ccomplex_t IT_0169 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0170 = IT_0009*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d4)*e_em*m_tau*IT_0029*conj(U_stau_00);
    const ccomplex_t IT_0173 = IT_0028*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + 0.5*IT_0174);
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0130*IT_0176;
    const ccomplex_t IT_0178 = IT_0106*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = m_tau*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0039 + (-0.5)*IT_0048 + (-0.5)
      *IT_0076 + (-0.5)*IT_0104 + (-0.5)*IT_0132 + 0.5*IT_0144 + 0.5*IT_0156 +
       0.5*IT_0168 + 0.5*IT_0180;
    const ccomplex_t IT_0182 = IT_0003*IT_0042;
    const ccomplex_t IT_0183 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = (-2)*IT_0184;
    const ccomplex_t IT_0186 = IT_0009*IT_0042;
    const ccomplex_t IT_0187 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0188 = IT_0186*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0191 = IT_0182*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = IT_0189 + IT_0192;
    const ccomplex_t IT_0194 = IT_0185 + IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = IT_0152*IT_0196;
    const ccomplex_t IT_0198 = IT_0077*IT_0078*IT_0197;
    const ccomplex_t IT_0199 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0200 = IT_0182*IT_0199;
    const ccomplex_t IT_0201 = (-2)*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0203 = IT_0186*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0206 = IT_0182*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = IT_0204 + IT_0207;
    const ccomplex_t IT_0209 = IT_0201 + IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = IT_0164*IT_0211;
    const ccomplex_t IT_0213 = IT_0000*IT_0038*IT_0212;
    const ccomplex_t IT_0214 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0215 = IT_0182*IT_0214;
    const ccomplex_t IT_0216 = (-2)*IT_0215;
    const ccomplex_t IT_0217 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0218 = IT_0186*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0221 = IT_0182*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = IT_0219 + IT_0222;
    const ccomplex_t IT_0224 = IT_0216 + IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0225;
    const ccomplex_t IT_0227 = IT_0176*IT_0226;
    const ccomplex_t IT_0228 = IT_0105*IT_0106*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0230 = IT_0182*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0233 = IT_0186*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = IT_0231 + IT_0234;
    const ccomplex_t IT_0236 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0237 = IT_0182*IT_0236;
    const ccomplex_t IT_0238 = (-2)*IT_0237;
    const ccomplex_t IT_0239 = IT_0235 + IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = 0.5*IT_0240;
    const ccomplex_t IT_0242 = IT_0073*IT_0241;
    const ccomplex_t IT_0243 = IT_0075*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = m_tau*IT_0244;
    const ccomplex_t IT_0246 = IT_0140*IT_0241;
    const ccomplex_t IT_0247 = IT_0049*IT_0075*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0003*IT_0027*IT_0029;
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = IT_0043*IT_0249;
    const ccomplex_t IT_0251 = IT_0046*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = IT_0089*IT_0196;
    const ccomplex_t IT_0254 = IT_0078*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = m_tau*IT_0255;
    const ccomplex_t IT_0257 = IT_0034*IT_0211;
    const ccomplex_t IT_0258 = IT_0038*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = m_tau*IT_0259;
    const ccomplex_t IT_0261 = IT_0117*IT_0226;
    const ccomplex_t IT_0262 = IT_0106*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = m_tau*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0198 + (-0.5)*IT_0213 + (-0.5)
      *IT_0228 + 0.5*IT_0245 + (-0.5)*IT_0247 + 0.5*IT_0252 + 0.5*IT_0256 + 0.5
      *IT_0260 + 0.5*IT_0264;
    const ccomplex_t IT_0266 = m_C_2*m_tau;
    const ccomplex_t IT_0267 = 0.5*IT_0244 + 0.5*IT_0255 + 0.5*IT_0259 + 0.5
      *IT_0263;
    const ccomplex_t IT_0268 = s_23*m_C_2;
    const ccomplex_t IT_0269 = 0.5*IT_0143 + 0.5*IT_0155 + 0.5*IT_0167 + 0.5
      *IT_0179;
    const ccomplex_t IT_0270 = s_12*m_tau;
    const ccomplex_t IT_0271 = 2*s_13;
    const ccomplex_t IT_0272 = (-2)*IT_0265;
    const ccomplex_t IT_0273 = m_C_2*m_tau*IT_0037;
    const ccomplex_t IT_0274 = s_12*s_23;
    const ccomplex_t IT_0275 = s_13*IT_0037;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0279 = 0.5*IT_0268;
    const ccomplex_t IT_0280 = (-2)*conj(IT_0181);
    const ccomplex_t IT_0281 = 0.5*IT_0270;
    const ccomplex_t IT_0282 = 2*IT_0181*(conj(IT_0265)*IT_0266 + -conj
      (IT_0267)*IT_0268 + -conj(IT_0269)*IT_0270 + 0.5*conj(IT_0181)*IT_0271) +
       2*IT_0265*(conj(IT_0181)*IT_0266 + 0.5*conj(IT_0265)*IT_0271) + (IT_0268
      *conj(IT_0269) + conj(IT_0267)*IT_0270)*IT_0272 + 2*IT_0267*(conj(IT_0269)
      *IT_0273 + 2*conj(IT_0267)*IT_0277 + IT_0279*IT_0280 + IT_0278*IT_0281) +
       2*IT_0269*(conj(IT_0267)*IT_0273 + 2*conj(IT_0269)*IT_0277 + IT_0278
      *IT_0279 + IT_0280*IT_0281);
    return create_ccomplex_return(IT_0282);
}

