#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_stau_1_to_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_stau_1_to_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = N_d3*e_em*m_tau*IT_0017*U_stau_00;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0013 + 0.5*IT_0020);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = pow(m_tau, 2);
    const ccomplex_t IT_0025 = pow(m_stau_1, 2);
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-0.5)*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0037 = IT_0010*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = N_d4*e_em*m_tau*IT_0017*U_stau_00;
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 0.5*IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0035*IT_0043;
    const ccomplex_t IT_0045 = cpow(s_23 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow(s_23 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0049 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0050 = IT_0001*IT_0049;
    const ccomplex_t IT_0051 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0057 = IT_0010*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = N_d2*e_em*m_tau*IT_0017*U_stau_00;
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = IT_0048*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow(s_23 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0069 = IT_0001*IT_0068;
    const ccomplex_t IT_0070 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0071 = IT_0001*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.5)*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0076 = IT_0010*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_d1*e_em*m_tau*IT_0017*U_stau_00;
    const ccomplex_t IT_0079 = IT_0016*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = IT_0067*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0028 + -IT_0047 + -IT_0066 + -IT_0085;
    const ccomplex_t IT_0087 = pow(m_W, -2);
    const ccomplex_t IT_0088 = pow(s_14, 2);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = pow(m_C_2, 2);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = s_13*IT_0092;
    const ccomplex_t IT_0094 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0095 = IT_0001*IT_0094;
    const ccomplex_t IT_0096 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0097 = IT_0001*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 0.5*IT_0098);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0102 = IT_0010*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0105 = IT_0001*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_d4)*e_em*m_tau*IT_0017*U_stau_10;
    const ccomplex_t IT_0108 = IT_0016*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + IT_0106 + 
      -IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0111;
    const ccomplex_t IT_0113 = IT_0045*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0116 = IT_0001*IT_0115;
    const ccomplex_t IT_0117 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + 0.5*IT_0119);
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0123 = IT_0010*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_d1)*e_em*m_tau*IT_0017*U_stau_10;
    const ccomplex_t IT_0129 = IT_0016*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0124 + IT_0127 + 
      -IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0121*IT_0132;
    const ccomplex_t IT_0134 = IT_0067*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0137 = IT_0001*IT_0136;
    const ccomplex_t IT_0138 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0139 = IT_0001*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + 0.5*IT_0140);
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0144 = IT_0010*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0147 = IT_0001*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = conj(N_d2)*e_em*m_tau*IT_0017*U_stau_10;
    const ccomplex_t IT_0150 = IT_0016*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + IT_0148 + 
      -IT_0151);
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = IT_0142*IT_0153;
    const ccomplex_t IT_0155 = IT_0048*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0158 = IT_0001*IT_0157;
    const ccomplex_t IT_0159 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0160 = IT_0001*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0158 + 0.5*IT_0161);
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0165 = IT_0010*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0168 = IT_0001*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_d3)*e_em*m_tau*IT_0017*U_stau_10;
    const ccomplex_t IT_0171 = IT_0016*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0166 + IT_0169 + 
      -IT_0172);
    const ccomplex_t IT_0174 = (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = IT_0163*IT_0174;
    const ccomplex_t IT_0176 = IT_0026*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0114 + IT_0135 + IT_0156 + IT_0177;
    const ccomplex_t IT_0179 = m_C_2*m_tau;
    const ccomplex_t IT_0180 = IT_0092*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0182 = IT_0008*IT_0174;
    const ccomplex_t IT_0183 = IT_0026*IT_0181*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0185 = IT_0035*IT_0111;
    const ccomplex_t IT_0186 = IT_0045*IT_0184*IT_0185;
    const ccomplex_t IT_0187 = m_C_2*IT_0156;
    const ccomplex_t IT_0188 = m_C_2*IT_0177;
    const ccomplex_t IT_0189 = m_C_2*IT_0114;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0191 = IT_0074*IT_0132;
    const ccomplex_t IT_0192 = IT_0067*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = m_C_2*IT_0135;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0195 = IT_0055*IT_0153;
    const ccomplex_t IT_0196 = IT_0048*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0183 + 0.5*IT_0186 + (-0.5)*IT_0187 + (
      -0.5)*IT_0188 + (-0.5)*IT_0189 + 0.5*IT_0192 + (-0.5)*IT_0193 + 0.5*IT_0196;
    const ccomplex_t IT_0198 = m_tau*IT_0092;
    const ccomplex_t IT_0199 = 0.5*IT_0028 + 0.5*IT_0047 + 0.5*IT_0066 + 0.5
      *IT_0085;
    const ccomplex_t IT_0200 = s_13*s_14;
    const ccomplex_t IT_0201 = s_34*IT_0090;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0207 = IT_0001*IT_0206;
    const ccomplex_t IT_0208 = conj(U_d2)*e_em*m_tau*IT_0017*U_stau_10;
    const ccomplex_t IT_0209 = IT_0016*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + (-0.5)*IT_0210);
    const ccomplex_t IT_0212 = IT_0205*IT_0211;
    const ccomplex_t IT_0213 = cpow(s_12 + 0.5*IT_0025 + 0.5*IT_0090 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0214 = IT_0212*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = pow(m_W, 2);
    const ccomplex_t IT_0218 = IT_0179*IT_0217;
    const ccomplex_t IT_0219 = 3*IT_0218;
    const ccomplex_t IT_0220 = s_14*IT_0217;
    const ccomplex_t IT_0221 = IT_0087*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = s_14 + IT_0222;
    const ccomplex_t IT_0224 = m_tau*IT_0223;
    const ccomplex_t IT_0225 = 2*IT_0224;
    const ccomplex_t IT_0226 = IT_0043*IT_0100;
    const ccomplex_t IT_0227 = IT_0045*IT_0184*IT_0226;
    const ccomplex_t IT_0228 = m_C_2*IT_0028;
    const ccomplex_t IT_0229 = m_C_2*IT_0047;
    const ccomplex_t IT_0230 = m_tau*IT_0215;
    const ccomplex_t IT_0231 = IT_0022*IT_0163;
    const ccomplex_t IT_0232 = IT_0026*IT_0181*IT_0231;
    const ccomplex_t IT_0233 = m_C_2*IT_0066;
    const ccomplex_t IT_0234 = IT_0082*IT_0121;
    const ccomplex_t IT_0235 = IT_0067*IT_0190*IT_0234;
    const ccomplex_t IT_0236 = IT_0063*IT_0142;
    const ccomplex_t IT_0237 = IT_0048*IT_0194*IT_0236;
    const ccomplex_t IT_0238 = m_C_2*IT_0085;
    const ccomplex_t IT_0239 = (-0.5)*IT_0227 + 0.5*IT_0228 + 0.5*IT_0229 +
       0.5*IT_0230 + (-0.5)*IT_0232 + 0.5*IT_0233 + (-0.5)*IT_0235 + (-0.5)
      *IT_0237 + 0.5*IT_0238;
    const ccomplex_t IT_0240 = s_34*IT_0217;
    const ccomplex_t IT_0241 = IT_0087*IT_0240;
    const ccomplex_t IT_0242 = (-0.25)*IT_0241;
    const ccomplex_t IT_0243 = s_34 + IT_0242;
    const ccomplex_t IT_0244 = m_C_2*IT_0243;
    const ccomplex_t IT_0245 = (-4)*IT_0244;
    const ccomplex_t IT_0246 = s_13*IT_0217;
    const ccomplex_t IT_0247 = pow(m_W, 4);
    const ccomplex_t IT_0248 = s_13*IT_0247;
    const ccomplex_t IT_0249 = IT_0087*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = s_14*s_34;
    const ccomplex_t IT_0252 = (-2)*IT_0251;
    const ccomplex_t IT_0253 = IT_0246 + IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = (-3)*IT_0218;
    const ccomplex_t IT_0255 = (-0.5)*IT_0177;
    const ccomplex_t IT_0256 = (-0.5)*IT_0114 + (-0.5)*IT_0135 + (-0.5)*IT_0156;
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = (-0.25)*IT_0249;
    const ccomplex_t IT_0259 = IT_0251 + IT_0258;
    const ccomplex_t IT_0260 = (-4)*IT_0255;
    const ccomplex_t IT_0261 = s_13*s_34;
    const ccomplex_t IT_0262 = s_14*IT_0024;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = IT_0261 + IT_0263;
    const ccomplex_t IT_0265 = -IT_0215;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = 2*(conj(IT_0255) + conj(IT_0256))*(IT_0178
      *IT_0203 + IT_0225*IT_0239 + IT_0197*IT_0245 + IT_0199*IT_0254 + (-4)
      *IT_0216*IT_0259 + (-2)*IT_0253*(IT_0256 + (-0.25)*IT_0260) + 0.5*IT_0264
      *IT_0266);
    const ccomplex_t IT_0268 = s_14*s_34*IT_0087;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = s_13 + IT_0269;
    const ccomplex_t IT_0271 = m_tau*IT_0270;
    const ccomplex_t IT_0272 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0273 = IT_0179*IT_0270;
    const ccomplex_t IT_0274 = m_C_2*IT_0270;
    const ccomplex_t IT_0275 = (-2)*conj(IT_0086);
    const ccomplex_t IT_0276 = pow(s_34, 2);
    const ccomplex_t IT_0277 = IT_0087*IT_0276;
    const ccomplex_t IT_0278 = -IT_0024;
    const ccomplex_t IT_0279 = IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = s_13*IT_0279;
    const ccomplex_t IT_0281 = m_C_2*IT_0279;
    const ccomplex_t IT_0282 = s_13*IT_0270;
    const ccomplex_t IT_0283 = (-2)*IT_0178;
    const ccomplex_t IT_0284 = 0.5*IT_0241;
    const ccomplex_t IT_0285 = s_34 + IT_0284;
    const ccomplex_t IT_0286 = m_C_2*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0286;
    const ccomplex_t IT_0288 = (-0.25)*IT_0221;
    const ccomplex_t IT_0289 = s_14 + IT_0288;
    const ccomplex_t IT_0290 = m_tau*IT_0289;
    const ccomplex_t IT_0291 = (-4)*IT_0290;
    const ccomplex_t IT_0292 = (-2)*conj(IT_0178);
    const ccomplex_t IT_0293 = IT_0246 + IT_0252;
    const ccomplex_t IT_0294 = IT_0087*IT_0293;
    const ccomplex_t IT_0295 = (-0.5)*IT_0294;
    const ccomplex_t IT_0296 = (-0.333333333333333)*IT_0216;
    const ccomplex_t IT_0297 = (-0.333333333333333)*IT_0199;
    const ccomplex_t IT_0298 = (-0.333333333333333)*IT_0257;
    const ccomplex_t IT_0299 = IT_0179*IT_0197 + (-0.333333333333333)*IT_0178
      *IT_0198 + (-0.166666666666667)*IT_0266*IT_0271 + 0.333333333333333
      *IT_0086*IT_0274 + (-0.666666666666667)*IT_0239*(s_13 + IT_0295) + IT_0291
      *IT_0296 + IT_0245*IT_0297 + IT_0225*IT_0298;
    const ccomplex_t IT_0300 = (-0.333333333333333)*IT_0086*IT_0198 + IT_0179
      *IT_0239 + (-0.333333333333333)*IT_0265*IT_0281 + (-0.166666666666667)
      *IT_0274*IT_0283 + (-0.666666666666667)*IT_0197*(s_13 + IT_0295) + IT_0287
      *IT_0296 + IT_0225*IT_0297 + IT_0245*IT_0298;
    const ccomplex_t IT_0301 = 2*conj(IT_0086)*(IT_0086*IT_0093 + IT_0178
      *IT_0180 + IT_0197*IT_0198 + IT_0199*IT_0203) + 2*conj(IT_0178)*(IT_0093
      *IT_0178 + IT_0086*IT_0180 + IT_0198*IT_0239 + IT_0203*IT_0257) + 2*conj
      (IT_0199)*(IT_0086*IT_0203 + IT_0216*IT_0219 + IT_0197*IT_0225 + IT_0239
      *IT_0245 + (-2)*IT_0199*IT_0253 + IT_0254*IT_0257) + IT_0267 + (IT_0257
      *IT_0264 + IT_0239*IT_0271)*IT_0272 + (IT_0265*IT_0273 + IT_0239*IT_0274)
      *IT_0275 + 2*conj(IT_0265)*(IT_0216*IT_0264 + -IT_0086*IT_0273 + IT_0265
      *IT_0280 + IT_0197*IT_0281 + 0.5*IT_0282*IT_0283) + 2*conj(IT_0216)*
      (IT_0199*IT_0219 + (-2)*IT_0216*IT_0253 + (-4)*IT_0257*IT_0259 + IT_0264
      *IT_0265 + 0.5*IT_0203*IT_0283 + IT_0197*IT_0287 + IT_0239*IT_0291) + 
      (IT_0203*IT_0216 + IT_0197*IT_0274 + IT_0265*IT_0282)*IT_0292 + (-6)*conj
      (IT_0239)*IT_0299 + (-6)*conj(IT_0197)*IT_0300;
    return create_ccomplex_return(IT_0301);
}

