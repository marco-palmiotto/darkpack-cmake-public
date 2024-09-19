#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_u_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_u_W(
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
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B1*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = N_W1*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + (-0.5)*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = pow(m_u, 2);
    const ccomplex_t IT_0021 = pow(m_su_L, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0020 + IT_0021 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = N_B4*e_em;
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = N_W4*e_em;
    const ccomplex_t IT_0029 = IT_0006*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + 3*IT_0030);
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0034 = IT_0006*IT_0033;
    const ccomplex_t IT_0035 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + (-0.5)*IT_0037);
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0032*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0020 + IT_0021 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = N_B3*e_em;
    const ccomplex_t IT_0045 = IT_0001*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_W3*e_em;
    const ccomplex_t IT_0048 = IT_0006*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 3*IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0053 = IT_0006*IT_0052;
    const ccomplex_t IT_0054 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + (-0.5)*IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0051*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0020 + IT_0021 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = N_B2*e_em;
    const ccomplex_t IT_0064 = IT_0001*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = N_W2*e_em;
    const ccomplex_t IT_0067 = IT_0006*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + 3*IT_0068);
    const ccomplex_t IT_0070 = 0.166666666666667*IT_0069;
    const ccomplex_t IT_0071 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0072 = IT_0006*IT_0071;
    const ccomplex_t IT_0073 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0074 = IT_0006*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + (-0.5)*IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0070*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0020 + IT_0021 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (-2)*IT_0024 + (-2)*IT_0043 + (-2)*IT_0062 + (
      -2)*IT_0081;
    const ccomplex_t IT_0083 = pow(m_W, -2);
    const ccomplex_t IT_0084 = pow(s_14, 2);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = pow(m_C_2, 2);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0085 + IT_0087;
    const ccomplex_t IT_0089 = s_13*IT_0088;
    const ccomplex_t IT_0090 = pow(m_W, -1);
    const ccomplex_t IT_0091 = sin(beta);
    const ccomplex_t IT_0092 = cpow(IT_0091, -1);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0006*IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0096 = IT_0006*IT_0095;
    const ccomplex_t IT_0097 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0098 = IT_0006*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + 0.5*IT_0099);
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0094*IT_0101;
    const ccomplex_t IT_0103 = IT_0060*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0006*IT_0090*IT_0092;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0108 = IT_0006*IT_0107;
    const ccomplex_t IT_0109 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0110 = IT_0006*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + 0.5*IT_0111);
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = IT_0041*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0006*IT_0090*IT_0092;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0120 = IT_0006*IT_0119;
    const ccomplex_t IT_0121 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0122 = IT_0006*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 0.5*IT_0123);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0118*IT_0125;
    const ccomplex_t IT_0127 = IT_0022*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0006*IT_0090*IT_0092;
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0132 = IT_0006*IT_0131;
    const ccomplex_t IT_0133 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0134 = IT_0006*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0130*IT_0137;
    const ccomplex_t IT_0139 = IT_0079*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0104 + 2*IT_0116 + 2*IT_0128 + 2*IT_0140;
    const ccomplex_t IT_0142 = m_u*m_C_2;
    const ccomplex_t IT_0143 = IT_0088*IT_0142;
    const ccomplex_t IT_0144 = pow(m_b, 2);
    const ccomplex_t IT_0145 = cpow((-2)*s_12 + -IT_0021 + -IT_0086 + IT_0144 
      + -reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0147 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0006
      *IT_0147*V_ub_mod;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0151 = cos(beta);
    const ccomplex_t IT_0152 = cpow(IT_0151, -1);
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *e_em*IT_0006*IT_0090*IT_0150*IT_0152*V_ub_mod;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0149*IT_0154;
    const ccomplex_t IT_0156 = IT_0145*IT_0146*IT_0155;
    const ccomplex_t IT_0157 = pow(m_d, 2);
    const ccomplex_t IT_0158 = cpow((-2)*s_12 + -IT_0021 + -IT_0086 + IT_0157 
      + -reg_prop, -1);
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0006;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0006*IT_0090*IT_0152;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = IT_0158*IT_0159*IT_0164;
    const ccomplex_t IT_0166 = m_C_2*IT_0081;
    const ccomplex_t IT_0167 = IT_0011*IT_0125;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0169 = IT_0022*IT_0167*IT_0168;
    const ccomplex_t IT_0170 = IT_0070*IT_0137;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0172 = IT_0079*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = IT_0051*IT_0101;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0175 = IT_0060*IT_0173*IT_0174;
    const ccomplex_t IT_0176 = IT_0032*IT_0113;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0178 = IT_0041*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = pow(m_s, 2);
    const ccomplex_t IT_0180 = cpow((-2)*s_12 + -IT_0021 + -IT_0086 + IT_0179 
      + -reg_prop, -1);
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0006;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0006*IT_0090*IT_0152;
    const ccomplex_t IT_0185 = (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = IT_0183*IT_0185;
    const ccomplex_t IT_0187 = IT_0180*IT_0181*IT_0186;
    const ccomplex_t IT_0188 = m_C_2*IT_0062;
    const ccomplex_t IT_0189 = m_C_2*IT_0024;
    const ccomplex_t IT_0190 = m_C_2*IT_0043;
    const ccomplex_t IT_0191 = IT_0156 + IT_0165 + -IT_0166 + IT_0169 +
       IT_0172 + IT_0175 + IT_0178 + IT_0187 + -IT_0188 + -IT_0189 + -IT_0190;
    const ccomplex_t IT_0192 = s_14*s_34*IT_0083;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = s_13 + IT_0193;
    const ccomplex_t IT_0195 = m_C_2*IT_0194;
    const ccomplex_t IT_0196 = IT_0018*IT_0118;
    const ccomplex_t IT_0197 = IT_0022*IT_0168*IT_0196;
    const ccomplex_t IT_0198 = IT_0058*IT_0094;
    const ccomplex_t IT_0199 = IT_0060*IT_0174*IT_0198;
    const ccomplex_t IT_0200 = IT_0039*IT_0106;
    const ccomplex_t IT_0201 = IT_0041*IT_0177*IT_0200;
    const ccomplex_t IT_0202 = IT_0077*IT_0130;
    const ccomplex_t IT_0203 = IT_0079*IT_0171*IT_0202;
    const ccomplex_t IT_0204 = -IT_0197 + -IT_0199 + -IT_0201 + -IT_0203;
    const ccomplex_t IT_0205 = m_C_2*IT_0104;
    const ccomplex_t IT_0206 = m_C_2*IT_0128;
    const ccomplex_t IT_0207 = m_C_2*IT_0140;
    const ccomplex_t IT_0208 = m_C_2*IT_0116;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0006*IT_0150
      *V_ub_mod;
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = IT_0149*IT_0210;
    const ccomplex_t IT_0212 = IT_0145*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = m_u*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0006;
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0161*IT_0216;
    const ccomplex_t IT_0218 = IT_0158*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_u*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*V_us*e_em*V_Wp2*IT_0006;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = IT_0183*IT_0222;
    const ccomplex_t IT_0224 = IT_0180*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = m_u*IT_0225;
    const ccomplex_t IT_0227 = IT_0205 + IT_0206 + IT_0207 + IT_0208 + IT_0214
       + IT_0220 + IT_0226;
    const ccomplex_t IT_0228 = m_u*IT_0088;
    const ccomplex_t IT_0229 = (-6)*IT_0082;
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = pow(m_W, 2);
    const ccomplex_t IT_0232 = s_13*IT_0231;
    const ccomplex_t IT_0233 = s_14*s_34;
    const ccomplex_t IT_0234 = (-2)*IT_0233;
    const ccomplex_t IT_0235 = IT_0232 + IT_0234;
    const ccomplex_t IT_0236 = IT_0083*IT_0235;
    const ccomplex_t IT_0237 = (-6)*IT_0236;
    const ccomplex_t IT_0238 = 12*s_13;
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = IT_0204 + IT_0227;
    const ccomplex_t IT_0241 = IT_0239*IT_0240;
    const ccomplex_t IT_0242 = m_u*IT_0194;
    const ccomplex_t IT_0243 = 2*IT_0213 + 2*IT_0219 + 2*IT_0225;
    const ccomplex_t IT_0244 = 6*IT_0243;
    const ccomplex_t IT_0245 = 6*IT_0141*IT_0195;
    const ccomplex_t IT_0246 = (-18)*IT_0142*IT_0191;
    const ccomplex_t IT_0247 = IT_0230 + IT_0241 + IT_0242*IT_0244 + IT_0245 +
       IT_0246;
    const ccomplex_t IT_0248 = s_13*IT_0194;
    const ccomplex_t IT_0249 = (-6)*conj(IT_0082);
    const ccomplex_t IT_0250 = pow(s_34, 2);
    const ccomplex_t IT_0251 = IT_0083*IT_0250;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0020 + IT_0252;
    const ccomplex_t IT_0254 = m_C_2*IT_0253;
    const ccomplex_t IT_0255 = IT_0024 + IT_0043 + IT_0062 + IT_0081;
    const ccomplex_t IT_0256 = s_14*IT_0020;
    const ccomplex_t IT_0257 = s_13*s_34;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0256 + IT_0258;
    const ccomplex_t IT_0260 = s_13*IT_0253;
    const ccomplex_t IT_0261 = (-6)*IT_0243;
    const ccomplex_t IT_0262 = -IT_0104 + -IT_0116 + -IT_0128 + -IT_0140;
    const ccomplex_t IT_0263 = s_13*s_14;
    const ccomplex_t IT_0264 = s_34*IT_0086;
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = IT_0142*IT_0194;
    const ccomplex_t IT_0268 = (-6)*conj(IT_0243);
    const ccomplex_t IT_0269 = s_34*IT_0231;
    const ccomplex_t IT_0270 = IT_0083*IT_0269;
    const ccomplex_t IT_0271 = (-0.25)*IT_0270;
    const ccomplex_t IT_0272 = s_34 + IT_0271;
    const ccomplex_t IT_0273 = m_C_2*IT_0272;
    const ccomplex_t IT_0274 = (-4)*IT_0273;
    const ccomplex_t IT_0275 = (-6)*conj(IT_0204);
    const ccomplex_t IT_0276 = (-6)*conj(IT_0227);
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = s_14*IT_0231;
    const ccomplex_t IT_0279 = IT_0083*IT_0278;
    const ccomplex_t IT_0280 = 0.5*IT_0279;
    const ccomplex_t IT_0281 = s_14 + IT_0280;
    const ccomplex_t IT_0282 = m_u*IT_0281;
    const ccomplex_t IT_0283 = 2*IT_0282;
    const ccomplex_t IT_0284 = (-6)*conj(IT_0191);
    const ccomplex_t IT_0285 = IT_0213 + IT_0219 + IT_0225;
    const ccomplex_t IT_0286 = (-0.25)*IT_0279;
    const ccomplex_t IT_0287 = s_14 + IT_0286;
    const ccomplex_t IT_0288 = m_u*IT_0287;
    const ccomplex_t IT_0289 = (-4)*IT_0288;
    const ccomplex_t IT_0290 = 0.5*IT_0270;
    const ccomplex_t IT_0291 = s_34 + IT_0290;
    const ccomplex_t IT_0292 = m_C_2*IT_0291;
    const ccomplex_t IT_0293 = 2*IT_0292;
    const ccomplex_t IT_0294 = IT_0142*IT_0231;
    const ccomplex_t IT_0295 = 3*IT_0294;
    const ccomplex_t IT_0296 = pow(m_W, 4);
    const ccomplex_t IT_0297 = s_13*IT_0296;
    const ccomplex_t IT_0298 = IT_0083*IT_0297;
    const ccomplex_t IT_0299 = 24*IT_0233;
    const ccomplex_t IT_0300 = (-12)*IT_0232;
    const ccomplex_t IT_0301 = 6*IT_0298 + IT_0299 + IT_0300;
    const ccomplex_t IT_0302 = (-6)*IT_0204;
    const ccomplex_t IT_0303 = (-6)*IT_0227;
    const ccomplex_t IT_0304 = IT_0302 + IT_0303;
    const ccomplex_t IT_0305 = (-6)*IT_0191;
    const ccomplex_t IT_0306 = (-24)*IT_0233;
    const ccomplex_t IT_0307 = 6*IT_0298 + IT_0306;
    const ccomplex_t IT_0308 = 0.166666666666667*IT_0255;
    const ccomplex_t IT_0309 = 6*conj(IT_0285)*(0.166666666666667*IT_0259
      *IT_0261 + 0.166666666666667*IT_0229*IT_0266 + IT_0262*IT_0295 +
       0.166666666666667*IT_0285*IT_0301 + 0.166666666666667*IT_0289*IT_0304 +
       0.166666666666667*IT_0293*IT_0305 + IT_0307*IT_0308);
    const ccomplex_t IT_0310 = (-3)*IT_0294;
    const ccomplex_t IT_0311 = 0.166666666666667*IT_0283;
    const ccomplex_t IT_0312 = 0.166666666666667*IT_0274;
    const ccomplex_t IT_0313 = 6*conj(IT_0082)*(IT_0082*IT_0089 + IT_0141
      *IT_0143 + IT_0191*IT_0195) + (conj(IT_0204) + conj(IT_0227))*IT_0247 + 
      (IT_0228*IT_0240 + IT_0243*IT_0248)*IT_0249 + 6*conj(IT_0191)*(IT_0082
      *IT_0195 + 1./6*IT_0142*((-18)*IT_0204 + (-18)*IT_0227) + -IT_0141*IT_0228
       + 1./6*IT_0191*IT_0239 + 1./6*IT_0244*IT_0254) + 6*conj(IT_0243)*(IT_0240
      *IT_0242 + 0.166666666666667*IT_0229*IT_0248 + IT_0191*IT_0254 + IT_0255
      *IT_0259 + 0.166666666666667*IT_0260*IT_0261) + 6*conj(IT_0141)*(IT_0089
      *IT_0141 + IT_0082*IT_0143 + -IT_0191*IT_0228 + IT_0195*IT_0240 + IT_0262
      *IT_0266 + 0.166666666666667*IT_0261*IT_0267) + IT_0141*IT_0267*IT_0268 +
       IT_0262*(IT_0274*IT_0277 + IT_0283*IT_0284) + IT_0285*(IT_0249*IT_0266 +
       IT_0259*IT_0268 + IT_0277*IT_0289 + IT_0284*IT_0293) + IT_0309 + 6
      *IT_0255*(conj(IT_0082)*IT_0266 + IT_0277*IT_0311 + IT_0284*IT_0312) + 6
      *conj(IT_0262)*(IT_0141*IT_0266 + IT_0285*IT_0295 + 0.166666666666667
      *IT_0262*IT_0301 + IT_0255*IT_0310 + IT_0305*IT_0311 + IT_0304*IT_0312) +
       6*conj(IT_0255)*(IT_0243*IT_0259 + IT_0082*IT_0266 + 0.166666666666667
      *IT_0285*IT_0307 + IT_0301*IT_0308 + IT_0262*IT_0310 + IT_0304*IT_0311 +
       IT_0305*IT_0312);
    return create_ccomplex_return(IT_0313);
}

