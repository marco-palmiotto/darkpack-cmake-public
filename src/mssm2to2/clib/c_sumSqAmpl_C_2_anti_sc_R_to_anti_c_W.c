#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_sc_R_to_anti_c_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_sc_R_to_anti_c_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
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
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sc_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0004;
    const ccomplex_t IT_0006 = 0.666666666666667*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0010 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0004;
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0004;
    const ccomplex_t IT_0034 = 0.666666666666667*IT_0033;
    const ccomplex_t IT_0035 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0036 = IT_0008*IT_0035;
    const ccomplex_t IT_0037 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0038 = IT_0008*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0034*IT_0041;
    const ccomplex_t IT_0043 = IT_0032*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0004;
    const ccomplex_t IT_0046 = 0.666666666666667*IT_0045;
    const ccomplex_t IT_0047 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0048 = IT_0008*IT_0047;
    const ccomplex_t IT_0049 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + 0.5*IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0018 + -IT_0031 + -IT_0044 + -IT_0057;
    const ccomplex_t IT_0059 = pow(m_W, -2);
    const ccomplex_t IT_0060 = pow(m_W, 2);
    const ccomplex_t IT_0061 = s_34*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = (-0.25)*IT_0062;
    const ccomplex_t IT_0064 = s_34 + IT_0063;
    const ccomplex_t IT_0065 = m_C_2*IT_0064;
    const ccomplex_t IT_0066 = (-4)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0069 = IT_0008*IT_0068;
    const ccomplex_t IT_0070 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.5)*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0034*IT_0074;
    const ccomplex_t IT_0076 = IT_0032*IT_0067*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0078 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0079 = IT_0008*IT_0078;
    const ccomplex_t IT_0080 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0081 = IT_0008*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + (-0.5)*IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = IT_0021*IT_0084;
    const ccomplex_t IT_0086 = IT_0019*IT_0077*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0088 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0089 = IT_0008*IT_0088;
    const ccomplex_t IT_0090 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0091 = IT_0008*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + (-0.5)*IT_0092);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = IT_0046*IT_0094;
    const ccomplex_t IT_0096 = IT_0055*IT_0087*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0098 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0099 = IT_0008*IT_0098;
    const ccomplex_t IT_0100 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0101 = IT_0008*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + (-0.5)*IT_0102);
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0006*IT_0104;
    const ccomplex_t IT_0106 = IT_0002*IT_0097*IT_0105;
    const ccomplex_t IT_0107 = -IT_0076 + -IT_0086 + -IT_0096 + -IT_0106;
    const ccomplex_t IT_0108 = (-6)*IT_0107;
    const ccomplex_t IT_0109 = m_C_2*IT_0057;
    const ccomplex_t IT_0110 = m_C_2*IT_0044;
    const ccomplex_t IT_0111 = m_C_2*IT_0018;
    const ccomplex_t IT_0112 = m_C_2*IT_0031;
    const ccomplex_t IT_0113 = pow(m_b, 2);
    const ccomplex_t IT_0114 = pow(m_C_2, 2);
    const ccomplex_t IT_0115 = cpow((-2)*s_12 + -IT_0001 + IT_0113 + -IT_0114 
      + -reg_prop, -1);
    const ccomplex_t IT_0116 = pow(m_W, -1);
    const ccomplex_t IT_0117 = sin(beta);
    const ccomplex_t IT_0118 = cpow(IT_0117, -1);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = IT_0115*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = m_c*IT_0125;
    const ccomplex_t IT_0127 = pow(m_d, 2);
    const ccomplex_t IT_0128 = cpow((-2)*s_12 + -IT_0001 + -IT_0114 + IT_0127 
      + -reg_prop, -1);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*IT_0008;
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = IT_0128*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_c*IT_0135;
    const ccomplex_t IT_0137 = pow(m_s, 2);
    const ccomplex_t IT_0138 = cpow((-2)*s_12 + -IT_0001 + -IT_0114 + IT_0137 
      + -reg_prop, -1);
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0008;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = IT_0140*IT_0142;
    const ccomplex_t IT_0144 = IT_0138*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = m_c*IT_0145;
    const ccomplex_t IT_0147 = IT_0109 + IT_0110 + IT_0111 + IT_0112 + IT_0126
       + IT_0136 + IT_0146;
    const ccomplex_t IT_0148 = (-6)*IT_0147;
    const ccomplex_t IT_0149 = IT_0108 + IT_0148;
    const ccomplex_t IT_0150 = s_14*IT_0060;
    const ccomplex_t IT_0151 = IT_0059*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = s_14 + IT_0152;
    const ccomplex_t IT_0154 = m_c*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = IT_0084*IT_0157;
    const ccomplex_t IT_0159 = IT_0019*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0162 = 0.5*IT_0161;
    const ccomplex_t IT_0163 = IT_0094*IT_0162;
    const ccomplex_t IT_0164 = IT_0055*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = IT_0104*IT_0167;
    const ccomplex_t IT_0169 = IT_0002*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0008*IT_0116*IT_0118;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0074*IT_0172;
    const ccomplex_t IT_0174 = IT_0032*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0160 + IT_0165 + IT_0170 + IT_0175;
    const ccomplex_t IT_0177 = m_C_2*IT_0176;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = (-6)*IT_0178;
    const ccomplex_t IT_0180 = IT_0041*IT_0172;
    const ccomplex_t IT_0181 = IT_0032*IT_0067*IT_0180;
    const ccomplex_t IT_0182 = IT_0028*IT_0157;
    const ccomplex_t IT_0183 = IT_0019*IT_0077*IT_0182;
    const ccomplex_t IT_0184 = IT_0053*IT_0162;
    const ccomplex_t IT_0185 = IT_0055*IT_0087*IT_0184;
    const ccomplex_t IT_0186 = IT_0015*IT_0167;
    const ccomplex_t IT_0187 = IT_0002*IT_0097*IT_0186;
    const ccomplex_t IT_0188 = IT_0181 + IT_0183 + IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = (-6)*IT_0188;
    const ccomplex_t IT_0190 = IT_0179 + IT_0189;
    const ccomplex_t IT_0191 = -IT_0160 + -IT_0165 + -IT_0170 + -IT_0175;
    const ccomplex_t IT_0192 = m_c*m_C_2;
    const ccomplex_t IT_0193 = IT_0060*IT_0192;
    const ccomplex_t IT_0194 = (-3)*IT_0193;
    const ccomplex_t IT_0195 = IT_0125 + IT_0135 + IT_0145;
    const ccomplex_t IT_0196 = 3*IT_0193;
    const ccomplex_t IT_0197 = 2*IT_0018 + 2*IT_0031 + 2*IT_0044 + 2*IT_0057;
    const ccomplex_t IT_0198 = s_13*s_14;
    const ccomplex_t IT_0199 = s_34*IT_0114;
    const ccomplex_t IT_0200 = -IT_0199;
    const ccomplex_t IT_0201 = IT_0198 + IT_0200;
    const ccomplex_t IT_0202 = pow(m_W, 4);
    const ccomplex_t IT_0203 = s_13*IT_0202;
    const ccomplex_t IT_0204 = IT_0059*IT_0203;
    const ccomplex_t IT_0205 = s_14*s_34;
    const ccomplex_t IT_0206 = 24*IT_0205;
    const ccomplex_t IT_0207 = s_13*IT_0060;
    const ccomplex_t IT_0208 = (-12)*IT_0207;
    const ccomplex_t IT_0209 = 6*IT_0204 + IT_0206 + IT_0208;
    const ccomplex_t IT_0210 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0211 = (-6)*conj(IT_0147);
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0214 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = 2*IT_0160 + 2*IT_0165 + 2*IT_0170 + 2*IT_0175;
    const ccomplex_t IT_0217 = pow(s_14, 2);
    const ccomplex_t IT_0218 = IT_0059*IT_0217;
    const ccomplex_t IT_0219 = -IT_0114;
    const ccomplex_t IT_0220 = IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = s_13*IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0125 + 2*IT_0135 + 2*IT_0145;
    const ccomplex_t IT_0223 = s_14*s_34*IT_0059;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = s_13 + IT_0224;
    const ccomplex_t IT_0226 = s_13*IT_0225;
    const ccomplex_t IT_0227 = m_c*IT_0220;
    const ccomplex_t IT_0228 = conj(IT_0107) + conj(IT_0147);
    const ccomplex_t IT_0229 = 6*IT_0216;
    const ccomplex_t IT_0230 = m_c*IT_0225;
    const ccomplex_t IT_0231 = pow(s_34, 2);
    const ccomplex_t IT_0232 = IT_0059*IT_0231;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0000 + IT_0233;
    const ccomplex_t IT_0235 = m_C_2*IT_0234;
    const ccomplex_t IT_0236 = conj(IT_0178) + conj(IT_0188);
    const ccomplex_t IT_0237 = 6*conj(IT_0216);
    const ccomplex_t IT_0238 = IT_0192*IT_0220;
    const ccomplex_t IT_0239 = m_C_2*IT_0225;
    const ccomplex_t IT_0240 = (-6)*IT_0216;
    const ccomplex_t IT_0241 = (-2)*IT_0205;
    const ccomplex_t IT_0242 = IT_0207 + IT_0241;
    const ccomplex_t IT_0243 = IT_0059*IT_0242;
    const ccomplex_t IT_0244 = (-6)*IT_0243;
    const ccomplex_t IT_0245 = 12*s_13;
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = (-18)*conj(IT_0107);
    const ccomplex_t IT_0248 = (-18)*conj(IT_0147);
    const ccomplex_t IT_0249 = IT_0247 + IT_0248;
    const ccomplex_t IT_0250 = (-6)*conj(IT_0216);
    const ccomplex_t IT_0251 = (-6)*conj(IT_0197);
    const ccomplex_t IT_0252 = s_13*IT_0234;
    const ccomplex_t IT_0253 = IT_0192*IT_0225;
    const ccomplex_t IT_0254 = s_14*IT_0000;
    const ccomplex_t IT_0255 = s_13*s_34;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = IT_0254 + IT_0256;
    const ccomplex_t IT_0258 = (-6)*IT_0222;
    const ccomplex_t IT_0259 = (-0.25)*IT_0151;
    const ccomplex_t IT_0260 = s_14 + IT_0259;
    const ccomplex_t IT_0261 = m_c*IT_0260;
    const ccomplex_t IT_0262 = (-4)*IT_0261;
    const ccomplex_t IT_0263 = 0.5*IT_0062;
    const ccomplex_t IT_0264 = s_34 + IT_0263;
    const ccomplex_t IT_0265 = m_C_2*IT_0264;
    const ccomplex_t IT_0266 = 2*IT_0265;
    const ccomplex_t IT_0267 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0268 = 6*conj(IT_0222);
    const ccomplex_t IT_0269 = 6*conj(IT_0191);
    const ccomplex_t IT_0270 = (-6)*IT_0204;
    const ccomplex_t IT_0271 = IT_0206 + IT_0270;
    const ccomplex_t IT_0272 = 0.166666666666667*conj(IT_0191);
    const ccomplex_t IT_0273 = 0.166666666666667*conj(IT_0195);
    const ccomplex_t IT_0274 = 0.166666666666667*IT_0201*IT_0237;
    const ccomplex_t IT_0275 = 0.166666666666667*IT_0257*IT_0267;
    const ccomplex_t IT_0276 = conj(IT_0058)*(IT_0066*IT_0149 + IT_0155
      *IT_0190) + (-6)*IT_0058*(conj(IT_0191)*IT_0194 + -conj(IT_0195)*IT_0196 +
       -conj(IT_0197)*IT_0201 + (-0.166666666666667)*conj(IT_0058)*IT_0209 + (
      -0.166666666666667)*IT_0066*IT_0212 + (-0.166666666666667)*IT_0155*IT_0215
      ) + ((conj(IT_0191) + conj(IT_0195))*IT_0201 + conj(IT_0216)*IT_0221 +
       conj(IT_0222)*IT_0226 + IT_0227*IT_0228)*IT_0229 + 6*IT_0222*(IT_0228
      *IT_0230 + IT_0235*IT_0236 + 0.166666666666667*IT_0226*IT_0237) + (conj
      (IT_0197)*IT_0238 + IT_0236*IT_0239)*IT_0240 + 6*IT_0188*(IT_0066*conj
      (IT_0191) + conj(IT_0222)*IT_0235 + 0.166666666666667*IT_0236*IT_0246 +
       0.166666666666667*IT_0192*IT_0249 + 0.166666666666667*IT_0239*IT_0250 +
       0.166666666666667*IT_0227*IT_0251) + (conj(IT_0222)*IT_0252 + conj
      (IT_0197)*IT_0253 + conj(IT_0191)*IT_0257)*IT_0258 + conj(IT_0195)*
      (IT_0257*IT_0258 + IT_0149*IT_0262 + IT_0190*IT_0266) + 6*IT_0197*(conj
      (IT_0058)*IT_0201 + conj(IT_0197)*IT_0221 + -IT_0227*IT_0236 + IT_0228
      *IT_0239 + 0.166666666666667*IT_0238*IT_0250 + 0.166666666666667*IT_0253
      *IT_0267) + IT_0178*(IT_0236*IT_0246 + IT_0192*IT_0249 + IT_0239*IT_0250 +
       IT_0227*IT_0251 + IT_0235*IT_0268 + IT_0066*IT_0269) + (IT_0107 + IT_0147
      )*(((-18)*conj(IT_0178) + (-18)*conj(IT_0188))*IT_0192 + IT_0227*IT_0237 +
       6*conj(IT_0197)*IT_0239 + IT_0228*IT_0246 + IT_0230*IT_0268 + IT_0155
      *IT_0269) + 6*IT_0195*(conj(IT_0058)*IT_0196 + 0.166666666666667*IT_0212
      *IT_0262 + 0.166666666666667*IT_0215*IT_0266 + IT_0271*IT_0272 + IT_0209
      *IT_0273 + IT_0274 + IT_0275) + (-6)*IT_0191*(conj(IT_0058)*IT_0194 + 
      -IT_0155*IT_0228 + -IT_0066*IT_0236 + -IT_0209*IT_0272 + -IT_0271*IT_0273 
      + -IT_0274 + -IT_0275);
    return create_ccomplex_return(IT_0276);
}

