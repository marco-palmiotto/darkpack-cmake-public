#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_R_to_anti_u_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_R_to_anti_u_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
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
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_u, 2);
    const ccomplex_t IT_0016 = pow(m_su_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = IT_0033*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = (-2)*IT_0019 + (-2)*IT_0032 + (-2)*IT_0045 + (
      -2)*IT_0058;
    const ccomplex_t IT_0060 = pow(m_W, -2);
    const ccomplex_t IT_0061 = pow(s_14, 2);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = pow(m_C_2, 2);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = s_13*IT_0065;
    const ccomplex_t IT_0067 = cos(theta_W);
    const ccomplex_t IT_0068 = cpow(IT_0067, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0068;
    const ccomplex_t IT_0070 = 0.666666666666667*IT_0069;
    const ccomplex_t IT_0071 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0072 = IT_0004*IT_0071;
    const ccomplex_t IT_0073 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0074 = IT_0004*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + 0.5*IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0070*IT_0077;
    const ccomplex_t IT_0079 = IT_0056*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0068;
    const ccomplex_t IT_0082 = 0.666666666666667*IT_0081;
    const ccomplex_t IT_0083 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0084 = IT_0004*IT_0083;
    const ccomplex_t IT_0085 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0004*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + 0.5*IT_0087);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0082*IT_0089;
    const ccomplex_t IT_0091 = IT_0017*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0068;
    const ccomplex_t IT_0094 = 0.666666666666667*IT_0093;
    const ccomplex_t IT_0095 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0096 = IT_0004*IT_0095;
    const ccomplex_t IT_0097 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0098 = IT_0004*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + 0.5*IT_0099);
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0094*IT_0101;
    const ccomplex_t IT_0103 = IT_0020*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0068;
    const ccomplex_t IT_0106 = 0.666666666666667*IT_0105;
    const ccomplex_t IT_0107 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0108 = IT_0004*IT_0107;
    const ccomplex_t IT_0109 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0110 = IT_0004*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + 0.5*IT_0111);
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = IT_0033*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0080 + 2*IT_0092 + 2*IT_0104 + 2*IT_0116;
    const ccomplex_t IT_0118 = m_u*m_C_2;
    const ccomplex_t IT_0119 = IT_0065*IT_0118;
    const ccomplex_t IT_0120 = IT_0047*IT_0077;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0122 = IT_0056*IT_0120*IT_0121;
    const ccomplex_t IT_0123 = IT_0006*IT_0089;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0125 = IT_0017*IT_0123*IT_0124;
    const ccomplex_t IT_0126 = m_C_2*IT_0032;
    const ccomplex_t IT_0127 = m_C_2*IT_0058;
    const ccomplex_t IT_0128 = m_C_2*IT_0019;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0130 = IT_0022*IT_0101;
    const ccomplex_t IT_0131 = IT_0020*IT_0129*IT_0130;
    const ccomplex_t IT_0132 = IT_0035*IT_0113;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0134 = IT_0033*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = m_C_2*IT_0045;
    const ccomplex_t IT_0136 = IT_0122 + IT_0125 + -IT_0126 + -IT_0127 + 
      -IT_0128 + IT_0131 + IT_0134 + -IT_0135;
    const ccomplex_t IT_0137 = s_14*s_34*IT_0060;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = s_13 + IT_0138;
    const ccomplex_t IT_0140 = m_C_2*IT_0139;
    const ccomplex_t IT_0141 = IT_0029*IT_0094;
    const ccomplex_t IT_0142 = IT_0020*IT_0129*IT_0141;
    const ccomplex_t IT_0143 = IT_0042*IT_0106;
    const ccomplex_t IT_0144 = IT_0033*IT_0133*IT_0143;
    const ccomplex_t IT_0145 = IT_0013*IT_0082;
    const ccomplex_t IT_0146 = IT_0017*IT_0124*IT_0145;
    const ccomplex_t IT_0147 = IT_0054*IT_0070;
    const ccomplex_t IT_0148 = IT_0056*IT_0121*IT_0147;
    const ccomplex_t IT_0149 = -IT_0142 + -IT_0144 + -IT_0146 + -IT_0148;
    const ccomplex_t IT_0150 = m_C_2*IT_0104;
    const ccomplex_t IT_0151 = m_C_2*IT_0116;
    const ccomplex_t IT_0152 = m_C_2*IT_0080;
    const ccomplex_t IT_0153 = m_C_2*IT_0092;
    const ccomplex_t IT_0154 = pow(m_b, 2);
    const ccomplex_t IT_0155 = cpow((-2)*s_12 + -IT_0016 + -IT_0063 + IT_0154 
      + -reg_prop, -1);
    const ccomplex_t IT_0156 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0000*IT_0002*IT_0004*IT_0156*V_ub_mod;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0159*V_ub_mod;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = IT_0158*IT_0161;
    const ccomplex_t IT_0163 = IT_0155*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = m_u*IT_0164;
    const ccomplex_t IT_0166 = pow(m_d, 2);
    const ccomplex_t IT_0167 = cpow((-2)*s_12 + -IT_0016 + -IT_0063 + IT_0166 
      + -reg_prop, -1);
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0004;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = IT_0167*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = m_u*IT_0174;
    const ccomplex_t IT_0176 = pow(m_s, 2);
    const ccomplex_t IT_0177 = cpow((-2)*s_12 + -IT_0016 + -IT_0063 + IT_0176 
      + -reg_prop, -1);
    const ccomplex_t IT_0178 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0004;
    const ccomplex_t IT_0181 = 0.5*IT_0180;
    const ccomplex_t IT_0182 = IT_0179*IT_0181;
    const ccomplex_t IT_0183 = IT_0177*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = m_u*IT_0184;
    const ccomplex_t IT_0186 = IT_0150 + IT_0151 + IT_0152 + IT_0153 + IT_0165
       + IT_0175 + IT_0185;
    const ccomplex_t IT_0187 = m_u*IT_0065;
    const ccomplex_t IT_0188 = (-6)*IT_0059;
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = pow(m_W, 2);
    const ccomplex_t IT_0191 = s_13*IT_0190;
    const ccomplex_t IT_0192 = s_14*s_34;
    const ccomplex_t IT_0193 = (-2)*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0060*IT_0194;
    const ccomplex_t IT_0196 = (-6)*IT_0195;
    const ccomplex_t IT_0197 = 12*s_13;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = IT_0149 + IT_0186;
    const ccomplex_t IT_0200 = IT_0198*IT_0199;
    const ccomplex_t IT_0201 = m_u*IT_0139;
    const ccomplex_t IT_0202 = 2*IT_0164 + 2*IT_0174 + 2*IT_0184;
    const ccomplex_t IT_0203 = 6*IT_0202;
    const ccomplex_t IT_0204 = 6*IT_0117*IT_0140;
    const ccomplex_t IT_0205 = (-18)*IT_0118*IT_0136;
    const ccomplex_t IT_0206 = IT_0189 + IT_0200 + IT_0201*IT_0203 + IT_0204 +
       IT_0205;
    const ccomplex_t IT_0207 = s_13*IT_0139;
    const ccomplex_t IT_0208 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0209 = pow(s_34, 2);
    const ccomplex_t IT_0210 = IT_0060*IT_0209;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0015 + IT_0211;
    const ccomplex_t IT_0213 = m_C_2*IT_0212;
    const ccomplex_t IT_0214 = IT_0019 + IT_0032 + IT_0045 + IT_0058;
    const ccomplex_t IT_0215 = s_14*IT_0015;
    const ccomplex_t IT_0216 = s_13*s_34;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0215 + IT_0217;
    const ccomplex_t IT_0219 = s_13*IT_0212;
    const ccomplex_t IT_0220 = (-6)*IT_0202;
    const ccomplex_t IT_0221 = -IT_0080 + -IT_0092 + -IT_0104 + -IT_0116;
    const ccomplex_t IT_0222 = s_13*s_14;
    const ccomplex_t IT_0223 = s_34*IT_0063;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = IT_0118*IT_0139;
    const ccomplex_t IT_0227 = (-6)*conj(IT_0202);
    const ccomplex_t IT_0228 = s_34*IT_0190;
    const ccomplex_t IT_0229 = IT_0060*IT_0228;
    const ccomplex_t IT_0230 = (-0.25)*IT_0229;
    const ccomplex_t IT_0231 = s_34 + IT_0230;
    const ccomplex_t IT_0232 = m_C_2*IT_0231;
    const ccomplex_t IT_0233 = (-4)*IT_0232;
    const ccomplex_t IT_0234 = (-6)*conj(IT_0149);
    const ccomplex_t IT_0235 = (-6)*conj(IT_0186);
    const ccomplex_t IT_0236 = IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = s_14*IT_0190;
    const ccomplex_t IT_0238 = IT_0060*IT_0237;
    const ccomplex_t IT_0239 = 0.5*IT_0238;
    const ccomplex_t IT_0240 = s_14 + IT_0239;
    const ccomplex_t IT_0241 = m_u*IT_0240;
    const ccomplex_t IT_0242 = 2*IT_0241;
    const ccomplex_t IT_0243 = (-6)*conj(IT_0136);
    const ccomplex_t IT_0244 = IT_0164 + IT_0174 + IT_0184;
    const ccomplex_t IT_0245 = (-0.25)*IT_0238;
    const ccomplex_t IT_0246 = s_14 + IT_0245;
    const ccomplex_t IT_0247 = m_u*IT_0246;
    const ccomplex_t IT_0248 = (-4)*IT_0247;
    const ccomplex_t IT_0249 = 0.5*IT_0229;
    const ccomplex_t IT_0250 = s_34 + IT_0249;
    const ccomplex_t IT_0251 = m_C_2*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0251;
    const ccomplex_t IT_0253 = IT_0118*IT_0190;
    const ccomplex_t IT_0254 = 3*IT_0253;
    const ccomplex_t IT_0255 = pow(m_W, 4);
    const ccomplex_t IT_0256 = s_13*IT_0255;
    const ccomplex_t IT_0257 = IT_0060*IT_0256;
    const ccomplex_t IT_0258 = 24*IT_0192;
    const ccomplex_t IT_0259 = (-12)*IT_0191;
    const ccomplex_t IT_0260 = 6*IT_0257 + IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = (-6)*IT_0149;
    const ccomplex_t IT_0262 = (-6)*IT_0186;
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = (-6)*IT_0136;
    const ccomplex_t IT_0265 = (-24)*IT_0192;
    const ccomplex_t IT_0266 = 6*IT_0257 + IT_0265;
    const ccomplex_t IT_0267 = 0.166666666666667*IT_0214;
    const ccomplex_t IT_0268 = 6*conj(IT_0244)*(0.166666666666667*IT_0218
      *IT_0220 + 0.166666666666667*IT_0188*IT_0225 + IT_0221*IT_0254 +
       0.166666666666667*IT_0244*IT_0260 + 0.166666666666667*IT_0248*IT_0263 +
       0.166666666666667*IT_0252*IT_0264 + IT_0266*IT_0267);
    const ccomplex_t IT_0269 = (-3)*IT_0253;
    const ccomplex_t IT_0270 = 0.166666666666667*IT_0242;
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0233;
    const ccomplex_t IT_0272 = 6*conj(IT_0059)*(IT_0059*IT_0066 + IT_0117
      *IT_0119 + IT_0136*IT_0140) + (conj(IT_0149) + conj(IT_0186))*IT_0206 + 
      (IT_0187*IT_0199 + IT_0202*IT_0207)*IT_0208 + 6*conj(IT_0136)*(IT_0059
      *IT_0140 + 1./6*IT_0118*((-18)*IT_0149 + (-18)*IT_0186) + -IT_0117*IT_0187
       + 1./6*IT_0136*IT_0198 + 1./6*IT_0203*IT_0213) + 6*conj(IT_0202)*(IT_0199
      *IT_0201 + 0.166666666666667*IT_0188*IT_0207 + IT_0136*IT_0213 + IT_0214
      *IT_0218 + 0.166666666666667*IT_0219*IT_0220) + 6*conj(IT_0117)*(IT_0066
      *IT_0117 + IT_0059*IT_0119 + -IT_0136*IT_0187 + IT_0140*IT_0199 + IT_0221
      *IT_0225 + 0.166666666666667*IT_0220*IT_0226) + IT_0117*IT_0226*IT_0227 +
       IT_0221*(IT_0233*IT_0236 + IT_0242*IT_0243) + IT_0244*(IT_0208*IT_0225 +
       IT_0218*IT_0227 + IT_0236*IT_0248 + IT_0243*IT_0252) + IT_0268 + 6
      *IT_0214*(conj(IT_0059)*IT_0225 + IT_0236*IT_0270 + IT_0243*IT_0271) + 6
      *conj(IT_0221)*(IT_0117*IT_0225 + IT_0244*IT_0254 + 0.166666666666667
      *IT_0221*IT_0260 + IT_0214*IT_0269 + IT_0264*IT_0270 + IT_0263*IT_0271) +
       6*conj(IT_0214)*(IT_0202*IT_0218 + IT_0059*IT_0225 + 0.166666666666667
      *IT_0244*IT_0266 + IT_0260*IT_0267 + IT_0221*IT_0269 + IT_0263*IT_0270 +
       IT_0264*IT_0271);
    return create_ccomplex_return(IT_0272);
}

