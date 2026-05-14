#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_1_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_1_to_b_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = pow(m_sb_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = m_b*conj(U_d2)*V_cb*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = V_cb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0014 + (-2)*IT_0016);
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -2);
    const ccomplex_t IT_0023 = pow(m_W, 4);
    const ccomplex_t IT_0024 = s_13*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = 6*IT_0025;
    const ccomplex_t IT_0027 = s_14*s_34;
    const ccomplex_t IT_0028 = (-24)*IT_0027;
    const ccomplex_t IT_0029 = IT_0026 + IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0002 + IT_0030 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0032 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0033 = IT_0005*IT_0032;
    const ccomplex_t IT_0034 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0035 = IT_0005*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = m_b*conj(N_d3)*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0043 = IT_0005*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = cos(theta_W);
    const ccomplex_t IT_0046 = cpow(IT_0045, -1);
    const ccomplex_t IT_0047 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0041 + -IT_0044 +
       0.333333333333333*IT_0049);
    const ccomplex_t IT_0051 = 3*IT_0050;
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0038*IT_0052;
    const ccomplex_t IT_0054 = IT_0031*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0002 + IT_0030 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0057 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0058 = IT_0005*IT_0057;
    const ccomplex_t IT_0059 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0060 = IT_0005*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = m_b*conj(N_d1)*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0065 = IT_0010*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0068 = IT_0005*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0071 = IT_0046*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + -IT_0069 +
       0.333333333333333*IT_0072);
    const ccomplex_t IT_0074 = 3*IT_0073;
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = IT_0063*IT_0075;
    const ccomplex_t IT_0077 = IT_0056*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0002 + IT_0030 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0080 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0081 = IT_0005*IT_0080;
    const ccomplex_t IT_0082 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0083 = IT_0005*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 0.5*IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = m_b*conj(N_d2)*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0088 = IT_0010*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0091 = IT_0005*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0094 = IT_0046*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + -IT_0092 +
       0.333333333333333*IT_0095);
    const ccomplex_t IT_0097 = 3*IT_0096;
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = IT_0086*IT_0098;
    const ccomplex_t IT_0100 = IT_0079*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0002 + IT_0030 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0103 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0104 = IT_0005*IT_0103;
    const ccomplex_t IT_0105 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0106 = IT_0005*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + 0.5*IT_0107);
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = m_b*conj(N_d4)*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0111 = IT_0010*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0114 = IT_0046*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0117 = IT_0005*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + 0.333333333333333
      *IT_0115 + -IT_0118);
    const ccomplex_t IT_0120 = 3*IT_0119;
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0122 = IT_0109*IT_0121;
    const ccomplex_t IT_0123 = IT_0102*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0078 + IT_0101 + IT_0124;
    const ccomplex_t IT_0126 = IT_0055 + IT_0125;
    const ccomplex_t IT_0127 = 24*IT_0027;
    const ccomplex_t IT_0128 = pow(m_W, 2);
    const ccomplex_t IT_0129 = s_13*IT_0128;
    const ccomplex_t IT_0130 = (-12)*IT_0129;
    const ccomplex_t IT_0131 = IT_0026 + IT_0127 + IT_0130;
    const ccomplex_t IT_0132 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0005;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(U_d2)*V_tb*e_em*IT_0011*U_sb_10;
    const ccomplex_t IT_0136 = IT_0010*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0139 = IT_0005*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0137 + (-2)*IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0134*IT_0141;
    const ccomplex_t IT_0143 = IT_0132*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = pow(m_u, 2);
    const ccomplex_t IT_0146 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + IT_0145 
      + -reg_prop, -1);
    const ccomplex_t IT_0147 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0148 = IT_0005*IT_0147;
    const ccomplex_t IT_0149 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = IT_0005*IT_0009*IT_0147;
    const ccomplex_t IT_0152 = m_b*conj(U_d2)*e_em*IT_0011*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0150 + (-0.5)*IT_0154);
    const ccomplex_t IT_0156 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *IT_0156*V_ub_mod;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = IT_0155*IT_0158;
    const ccomplex_t IT_0160 = IT_0146*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0144 + IT_0161;
    const ccomplex_t IT_0163 = IT_0021 + IT_0162;
    const ccomplex_t IT_0164 = IT_0029*IT_0163;
    const ccomplex_t IT_0165 = (IT_0055 + IT_0125)*IT_0131 + IT_0164;
    const ccomplex_t IT_0166 = s_34*IT_0128;
    const ccomplex_t IT_0167 = IT_0022*IT_0166;
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = s_34 + IT_0168;
    const ccomplex_t IT_0170 = m_C_2*IT_0169;
    const ccomplex_t IT_0171 = 2*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0173 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0174 = IT_0005*IT_0173;
    const ccomplex_t IT_0175 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0176 = IT_0005*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + (-0.5)*IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0075*IT_0179;
    const ccomplex_t IT_0181 = IT_0056*IT_0172*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0183 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0184 = IT_0005*IT_0183;
    const ccomplex_t IT_0185 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0186 = IT_0005*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + (-0.5)*IT_0187);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0098*IT_0189;
    const ccomplex_t IT_0191 = IT_0079*IT_0182*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0193 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0194 = IT_0005*IT_0193;
    const ccomplex_t IT_0195 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0196 = IT_0005*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0194 + (-0.5)*IT_0197);
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = IT_0052*IT_0199;
    const ccomplex_t IT_0201 = IT_0031*IT_0192*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0203 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0204 = IT_0005*IT_0203;
    const ccomplex_t IT_0205 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0206 = IT_0005*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (-0.5)*IT_0207);
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0121*IT_0209;
    const ccomplex_t IT_0211 = IT_0102*IT_0202*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0213 = sin(beta);
    const ccomplex_t IT_0214 = cpow(IT_0213, -1);
    const ccomplex_t IT_0215 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0005*IT_0011*IT_0214*U_sb_00;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = IT_0007*IT_0216;
    const ccomplex_t IT_0218 = IT_0003*IT_0212*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0005*IT_0011*IT_0214*U_sb_00;
    const ccomplex_t IT_0221 = 0.5*IT_0220;
    const ccomplex_t IT_0222 = IT_0134*IT_0221;
    const ccomplex_t IT_0223 = IT_0132*IT_0219*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0005*IT_0011*IT_0147*IT_0214*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0226 = 0.5*IT_0225;
    const ccomplex_t IT_0227 = IT_0158*IT_0226;
    const ccomplex_t IT_0228 = IT_0146*IT_0224*IT_0227;
    const ccomplex_t IT_0229 = m_C_2*IT_0078;
    const ccomplex_t IT_0230 = m_C_2*IT_0101;
    const ccomplex_t IT_0231 = m_C_2*IT_0055;
    const ccomplex_t IT_0232 = m_C_2*IT_0124;
    const ccomplex_t IT_0233 = -IT_0181 + -IT_0191 + -IT_0201 + -IT_0211 + 
      -IT_0218 + -IT_0223 + -IT_0228 + IT_0229 + IT_0230 + IT_0231 + IT_0232;
    const ccomplex_t IT_0234 = s_14*IT_0128;
    const ccomplex_t IT_0235 = IT_0022*IT_0234;
    const ccomplex_t IT_0236 = (-0.25)*IT_0235;
    const ccomplex_t IT_0237 = s_14 + IT_0236;
    const ccomplex_t IT_0238 = m_b*IT_0237;
    const ccomplex_t IT_0239 = (-4)*IT_0238;
    const ccomplex_t IT_0240 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0241 = IT_0046*IT_0240;
    const ccomplex_t IT_0242 = 1.4142135623731*IT_0241;
    const ccomplex_t IT_0243 = m_b*N_d1*e_em*IT_0011*U_sb_00;
    const ccomplex_t IT_0244 = IT_0010*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + 1.5*IT_0245);
    const ccomplex_t IT_0247 = (-0.333333333333333)*IT_0246;
    const ccomplex_t IT_0248 = IT_0063*IT_0247;
    const ccomplex_t IT_0249 = IT_0056*IT_0172*IT_0248;
    const ccomplex_t IT_0250 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0251 = IT_0046*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = m_b*N_d2*e_em*IT_0011*U_sb_00;
    const ccomplex_t IT_0254 = IT_0010*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*(IT_0252 + 1.5*IT_0255);
    const ccomplex_t IT_0257 = (-0.333333333333333)*IT_0256;
    const ccomplex_t IT_0258 = IT_0086*IT_0257;
    const ccomplex_t IT_0259 = IT_0079*IT_0182*IT_0258;
    const ccomplex_t IT_0260 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0261 = IT_0046*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = m_b*N_d3*e_em*IT_0011*U_sb_00;
    const ccomplex_t IT_0264 = IT_0010*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*(IT_0262 + 1.5*IT_0265);
    const ccomplex_t IT_0267 = (-0.333333333333333)*IT_0266;
    const ccomplex_t IT_0268 = IT_0038*IT_0267;
    const ccomplex_t IT_0269 = IT_0031*IT_0192*IT_0268;
    const ccomplex_t IT_0270 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0271 = IT_0046*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = m_b*N_d4*e_em*IT_0011*U_sb_00;
    const ccomplex_t IT_0274 = IT_0010*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*(IT_0272 + 1.5*IT_0275);
    const ccomplex_t IT_0277 = (-0.333333333333333)*IT_0276;
    const ccomplex_t IT_0278 = IT_0109*IT_0277;
    const ccomplex_t IT_0279 = IT_0102*IT_0202*IT_0278;
    const ccomplex_t IT_0280 = IT_0179*IT_0247;
    const ccomplex_t IT_0281 = IT_0056*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = m_C_2*IT_0282;
    const ccomplex_t IT_0284 = IT_0189*IT_0257;
    const ccomplex_t IT_0285 = IT_0079*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = m_C_2*IT_0286;
    const ccomplex_t IT_0288 = IT_0199*IT_0267;
    const ccomplex_t IT_0289 = IT_0031*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = m_C_2*IT_0290;
    const ccomplex_t IT_0292 = IT_0209*IT_0277;
    const ccomplex_t IT_0293 = IT_0102*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = m_C_2*IT_0294;
    const ccomplex_t IT_0296 = m_b*IT_0021;
    const ccomplex_t IT_0297 = m_b*IT_0144;
    const ccomplex_t IT_0298 = m_b*IT_0161;
    const ccomplex_t IT_0299 = IT_0249 + IT_0259 + IT_0269 + IT_0279 + 
      -IT_0283 + -IT_0287 + -IT_0291 + -IT_0295 + -IT_0296 + -IT_0297 + -IT_0298;
    const ccomplex_t IT_0300 = m_b*m_C_2;
    const ccomplex_t IT_0301 = IT_0128*IT_0300;
    const ccomplex_t IT_0302 = 3*IT_0301;
    const ccomplex_t IT_0303 = -IT_0282 + -IT_0286 + -IT_0290 + -IT_0294;
    const ccomplex_t IT_0304 = 6*conj(IT_0021);
    const ccomplex_t IT_0305 = (-2)*IT_0055 + (-2)*IT_0078 + (-2)*IT_0101 + (
      -2)*IT_0124;
    const ccomplex_t IT_0306 = s_13*s_14;
    const ccomplex_t IT_0307 = s_34*IT_0001;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = 6*IT_0309;
    const ccomplex_t IT_0311 = 2*IT_0021 + 2*IT_0144 + 2*IT_0161;
    const ccomplex_t IT_0312 = s_14*IT_0030;
    const ccomplex_t IT_0313 = s_13*s_34;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = IT_0312 + IT_0314;
    const ccomplex_t IT_0316 = 6*IT_0315;
    const ccomplex_t IT_0317 = (-3)*IT_0301;
    const ccomplex_t IT_0318 = 6*IT_0317;
    const ccomplex_t IT_0319 = (-0.25)*IT_0167;
    const ccomplex_t IT_0320 = s_34 + IT_0319;
    const ccomplex_t IT_0321 = m_C_2*IT_0320;
    const ccomplex_t IT_0322 = (-4)*IT_0321;
    const ccomplex_t IT_0323 = 6*IT_0322;
    const ccomplex_t IT_0324 = 0.5*IT_0235;
    const ccomplex_t IT_0325 = s_14 + IT_0324;
    const ccomplex_t IT_0326 = m_b*IT_0325;
    const ccomplex_t IT_0327 = 2*IT_0326;
    const ccomplex_t IT_0328 = 6*IT_0327;
    const ccomplex_t IT_0329 = 2*IT_0282 + 2*IT_0286 + 2*IT_0290 + 2*IT_0294;
    const ccomplex_t IT_0330 = 6*IT_0055;
    const ccomplex_t IT_0331 = 6*IT_0021;
    const ccomplex_t IT_0332 = 6*IT_0302;
    const ccomplex_t IT_0333 = pow(s_14, 2);
    const ccomplex_t IT_0334 = IT_0022*IT_0333;
    const ccomplex_t IT_0335 = -IT_0001;
    const ccomplex_t IT_0336 = IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = s_13*IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0337;
    const ccomplex_t IT_0339 = m_b*IT_0336;
    const ccomplex_t IT_0340 = 6*IT_0339;
    const ccomplex_t IT_0341 = IT_0300*IT_0336;
    const ccomplex_t IT_0342 = 6*IT_0341;
    const ccomplex_t IT_0343 = s_14*s_34*IT_0022;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = s_13 + IT_0344;
    const ccomplex_t IT_0346 = m_C_2*IT_0345;
    const ccomplex_t IT_0347 = (-6)*IT_0346;
    const ccomplex_t IT_0348 = IT_0300*IT_0345;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = 12*s_13;
    const ccomplex_t IT_0351 = (-2)*IT_0027;
    const ccomplex_t IT_0352 = IT_0129 + IT_0351;
    const ccomplex_t IT_0353 = IT_0022*IT_0352;
    const ccomplex_t IT_0354 = (-6)*IT_0353;
    const ccomplex_t IT_0355 = IT_0350 + IT_0354;
    const ccomplex_t IT_0356 = 6*IT_0171;
    const ccomplex_t IT_0357 = (-18)*IT_0300;
    const ccomplex_t IT_0358 = pow(s_34, 2);
    const ccomplex_t IT_0359 = IT_0022*IT_0358;
    const ccomplex_t IT_0360 = -IT_0359;
    const ccomplex_t IT_0361 = IT_0030 + IT_0360;
    const ccomplex_t IT_0362 = m_C_2*IT_0361;
    const ccomplex_t IT_0363 = (-6)*IT_0362;
    const ccomplex_t IT_0364 = conj(IT_0233)*(IT_0125*IT_0323 + IT_0303
      *IT_0328 + IT_0322*IT_0330 + IT_0171*IT_0331 + IT_0329*IT_0340 + IT_0305
      *IT_0347 + IT_0233*IT_0355 + IT_0162*IT_0356 + IT_0299*IT_0357 + IT_0311
      *IT_0363);
    const ccomplex_t IT_0365 = 6*IT_0239;
    const ccomplex_t IT_0366 = m_b*IT_0345;
    const ccomplex_t IT_0367 = (-6)*IT_0366;
    const ccomplex_t IT_0368 = conj(IT_0299)*(IT_0303*IT_0323 + IT_0125
      *IT_0328 + IT_0327*IT_0330 + IT_0239*IT_0331 + IT_0305*IT_0340 + IT_0329
      *IT_0347 + IT_0299*IT_0355 + IT_0233*IT_0357 + IT_0162*IT_0365 + IT_0311
      *IT_0367);
    const ccomplex_t IT_0369 = (-6)*conj(IT_0021);
    const ccomplex_t IT_0370 = IT_0305*IT_0309;
    const ccomplex_t IT_0371 = 6*conj(IT_0055);
    const ccomplex_t IT_0372 = (-6)*IT_0315;
    const ccomplex_t IT_0373 = (-0.166666666666667)*IT_0233;
    const ccomplex_t IT_0374 = (-0.166666666666667)*IT_0299;
    const ccomplex_t IT_0375 = s_13*IT_0361;
    const ccomplex_t IT_0376 = s_13*IT_0345;
    const ccomplex_t IT_0377 = (-6)*IT_0376;
    const ccomplex_t IT_0378 = (-0.166666666666667)*IT_0305;
    const ccomplex_t IT_0379 = (-0.166666666666667)*IT_0125;
    const ccomplex_t IT_0380 = (-0.166666666666667)*IT_0330;
    const ccomplex_t IT_0381 = (-6)*conj(IT_0311)*((-0.166666666666667)
      *IT_0329*IT_0349 + (-0.166666666666667)*IT_0162*IT_0372 + IT_0363*IT_0373 
      + IT_0367*IT_0374 + IT_0311*IT_0375 + IT_0377*IT_0378 + IT_0316*IT_0379 +
       IT_0315*(IT_0021 + IT_0380));
    const ccomplex_t IT_0382 = conj(IT_0305)*(IT_0329*IT_0342 + (-6)*IT_0347
      *IT_0373 + (-6)*IT_0340*IT_0374 + IT_0311*IT_0377 + (-6)*IT_0338*IT_0378 +
       (-6)*IT_0310*IT_0379 + (-6)*IT_0309*(IT_0163 + IT_0380));
    const ccomplex_t IT_0383 = conj(IT_0021)*(IT_0029*IT_0126 + IT_0131
      *IT_0163) + conj(IT_0055)*IT_0165 + (IT_0171*IT_0233 + IT_0239*IT_0299 +
       IT_0302*IT_0303)*IT_0304 + conj(IT_0125)*(IT_0165 + IT_0305*IT_0310 +
       IT_0311*IT_0316 + IT_0303*IT_0318 + IT_0233*IT_0323 + IT_0299*IT_0328) +
       conj(IT_0303)*(IT_0131*IT_0303 + IT_0125*IT_0318 + IT_0299*IT_0323 +
       IT_0233*IT_0328 + IT_0310*IT_0329 + IT_0317*IT_0330 + IT_0302*IT_0331 +
       IT_0162*IT_0332) + conj(IT_0329)*(IT_0303*IT_0310 + IT_0329*IT_0338 +
       IT_0233*IT_0340 + IT_0305*IT_0342 + IT_0299*IT_0347 + IT_0311*IT_0349) +
       IT_0364 + IT_0368 + IT_0369*(IT_0311*IT_0315 + IT_0370) + (IT_0311
      *IT_0315 + IT_0303*IT_0317 + IT_0233*IT_0322 + IT_0299*IT_0327 + IT_0370)
      *IT_0371 + conj(IT_0162)*(IT_0029*IT_0126 + IT_0131*IT_0163 + IT_0303
      *IT_0332 + (-6)*IT_0370 + IT_0311*IT_0372 + (-6)*IT_0356*IT_0373 + (-6)
      *IT_0365*IT_0374) + IT_0381 + IT_0382;
    return create_ccomplex_return(IT_0383);
}

