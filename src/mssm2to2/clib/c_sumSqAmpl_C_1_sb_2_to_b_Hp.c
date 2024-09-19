#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_2_to_b_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_2_to_b_Hp(
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
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_2 = param->m_sb_2;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0022 = IT_0009*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0001, -1);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = pow(m_W, -1);
    const ccomplex_t IT_0027 = m_b*N_d3*e_em*IT_0026*U_sb_01;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + 1.5*IT_0029);
    const ccomplex_t IT_0031 = (-0.333333333333333)*IT_0030;
    const ccomplex_t IT_0032 = IT_0020*IT_0031;
    const ccomplex_t IT_0033 = pow(m_b, 2);
    const ccomplex_t IT_0034 = pow(m_sb_2, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0000*IT_0032*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0038 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = 2*IT_0039;
    const ccomplex_t IT_0041 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0042 = IT_0010*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0045 = IT_0004*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = IT_0043 + IT_0046;
    const ccomplex_t IT_0048 = IT_0040 + IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0052 = IT_0009*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_b*N_d4*e_em*IT_0026*U_sb_01;
    const ccomplex_t IT_0055 = IT_0025*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 1.5*IT_0056);
    const ccomplex_t IT_0058 = (-0.333333333333333)*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0037*IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0063 = pow(m_C_1, 2);
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + -IT_0034 + -IT_0063 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0065 = sin(beta);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0003*IT_0024*IT_0026*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = cpow(IT_0065, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0026*IT_0068*U_sb_01;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0067*IT_0070;
    const ccomplex_t IT_0072 = IT_0062*IT_0064*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0075 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0076 = IT_0009*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_b*N_d2*e_em*IT_0026*U_sb_01;
    const ccomplex_t IT_0079 = IT_0025*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 1.5*IT_0080);
    const ccomplex_t IT_0082 = (-0.333333333333333)*IT_0081;
    const ccomplex_t IT_0083 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0084 = IT_0010*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = 2*IT_0091;
    const ccomplex_t IT_0093 = IT_0089 + IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0082*IT_0095;
    const ccomplex_t IT_0097 = IT_0073*IT_0074*IT_0096;
    const ccomplex_t IT_0098 = pow(m_u, 2);
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + -IT_0034 + -IT_0063 + IT_0098 
      + -reg_prop, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0101 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0003*IT_0024*IT_0026*IT_0065*IT_0101*V_ub_mod;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0003*IT_0026*IT_0068*IT_0104*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0103*IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0100*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0033 + IT_0034 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0111 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0112 = IT_0004*IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0112;
    const ccomplex_t IT_0114 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0115 = IT_0010*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0118 = IT_0004*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = IT_0116 + IT_0119;
    const ccomplex_t IT_0121 = IT_0113 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0125 = IT_0009*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_b*N_d1*e_em*IT_0026*U_sb_01;
    const ccomplex_t IT_0128 = IT_0025*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + 1.5*IT_0129);
    const ccomplex_t IT_0131 = (-0.333333333333333)*IT_0130;
    const ccomplex_t IT_0132 = IT_0123*IT_0131;
    const ccomplex_t IT_0133 = IT_0109*IT_0110*IT_0132;
    const ccomplex_t IT_0134 = pow(m_c, 2);
    const ccomplex_t IT_0135 = cpow((-2)*s_12 + -IT_0034 + -IT_0063 + IT_0134 
      + -reg_prop, -1);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0003*IT_0024*IT_0026*IT_0065;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0003*IT_0026*IT_0068*U_sb_01;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = IT_0135*IT_0136*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0144 = IT_0009*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0147 = IT_0003*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = m_b*conj(N_d1)*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0150 = IT_0025*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + (-3)*IT_0148 + 3
      *IT_0151);
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0152;
    const ccomplex_t IT_0154 = IT_0123*IT_0153;
    const ccomplex_t IT_0155 = IT_0110*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_b*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0159 = IT_0009*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0162 = IT_0003*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(N_d2)*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0165 = IT_0025*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + (-3)*IT_0163 + 3
      *IT_0166);
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0169 = IT_0095*IT_0168;
    const ccomplex_t IT_0170 = IT_0074*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_b*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0174 = IT_0009*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0177 = IT_0003*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = m_b*conj(N_d3)*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0180 = IT_0025*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0175 + (-3)*IT_0178 + 3
      *IT_0181);
    const ccomplex_t IT_0183 = 0.166666666666667*IT_0182;
    const ccomplex_t IT_0184 = IT_0020*IT_0183;
    const ccomplex_t IT_0185 = IT_0035*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = m_b*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0189 = IT_0009*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0192 = IT_0003*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*conj(N_d4)*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0195 = IT_0025*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0190 + (-3)*IT_0193 + 3
      *IT_0196);
    const ccomplex_t IT_0198 = 0.166666666666667*IT_0197;
    const ccomplex_t IT_0199 = IT_0050*IT_0198;
    const ccomplex_t IT_0200 = IT_0060*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_b*IT_0201;
    const ccomplex_t IT_0203 = V_cb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0204 = IT_0003*IT_0203;
    const ccomplex_t IT_0205 = m_b*conj(U_d1)*V_cb*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0206 = IT_0025*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (-0.5)*IT_0207);
    const ccomplex_t IT_0209 = IT_0138*IT_0208;
    const ccomplex_t IT_0210 = IT_0135*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = m_C_1*IT_0211;
    const ccomplex_t IT_0213 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0214 = IT_0003*IT_0213;
    const ccomplex_t IT_0215 = m_b*conj(U_d1)*V_tb*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0216 = IT_0025*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*(IT_0214 + (-0.5)*IT_0217);
    const ccomplex_t IT_0219 = IT_0067*IT_0218;
    const ccomplex_t IT_0220 = IT_0064*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = m_C_1*IT_0221;
    const ccomplex_t IT_0223 = IT_0003*IT_0104;
    const ccomplex_t IT_0224 = e_em*conj(U_Wm1)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0225 = IT_0223*IT_0224;
    const ccomplex_t IT_0226 = IT_0003*IT_0024*IT_0104;
    const ccomplex_t IT_0227 = m_b*conj(U_d1)*e_em*IT_0026*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*(IT_0225 + (-0.5)*IT_0229);
    const ccomplex_t IT_0231 = IT_0103*IT_0230;
    const ccomplex_t IT_0232 = IT_0099*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = m_C_1*IT_0233;
    const ccomplex_t IT_0235 = -IT_0036 + -IT_0061 + -IT_0072 + -IT_0097 + 
      -IT_0108 + -IT_0133 + -IT_0142 + IT_0157 + IT_0172 + IT_0187 + IT_0202 +
       IT_0212 + IT_0222 + IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0001*IT_0003*IT_0026*IT_0068*IT_0101*V_ub_mod;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0230*IT_0237;
    const ccomplex_t IT_0239 = IT_0099*IT_0100*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0001*IT_0003*IT_0026*IT_0068;
    const ccomplex_t IT_0241 = 0.5*IT_0240;
    const ccomplex_t IT_0242 = IT_0140*IT_0241;
    const ccomplex_t IT_0243 = IT_0135*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = m_C_1*IT_0244;
    const ccomplex_t IT_0246 = IT_0003*IT_0065;
    const ccomplex_t IT_0247 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0248 = IT_0246*IT_0247;
    const ccomplex_t IT_0249 = (-2)*IT_0248;
    const ccomplex_t IT_0250 = IT_0009*IT_0065;
    const ccomplex_t IT_0251 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0255 = IT_0246*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = IT_0253 + IT_0256;
    const ccomplex_t IT_0258 = IT_0249 + IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0259;
    const ccomplex_t IT_0261 = IT_0168*IT_0260;
    const ccomplex_t IT_0262 = IT_0073*IT_0074*IT_0261;
    const ccomplex_t IT_0263 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0264 = IT_0246*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0264;
    const ccomplex_t IT_0266 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0267 = IT_0250*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0270 = IT_0246*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = IT_0268 + IT_0271;
    const ccomplex_t IT_0273 = IT_0265 + IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = 0.5*IT_0274;
    const ccomplex_t IT_0276 = IT_0198*IT_0275;
    const ccomplex_t IT_0277 = IT_0037*IT_0060*IT_0276;
    const ccomplex_t IT_0278 = IT_0208*IT_0241;
    const ccomplex_t IT_0279 = IT_0135*IT_0136*IT_0278;
    const ccomplex_t IT_0280 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0281 = IT_0250*IT_0280;
    const ccomplex_t IT_0282 = 1.4142135623731*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0284 = IT_0246*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = IT_0282 + IT_0285;
    const ccomplex_t IT_0287 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0288 = IT_0246*IT_0287;
    const ccomplex_t IT_0289 = (-2)*IT_0288;
    const ccomplex_t IT_0290 = IT_0286 + IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = 0.5*IT_0291;
    const ccomplex_t IT_0293 = IT_0131*IT_0292;
    const ccomplex_t IT_0294 = IT_0110*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = m_b*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0003*IT_0026*IT_0068;
    const ccomplex_t IT_0298 = 0.5*IT_0297;
    const ccomplex_t IT_0299 = IT_0070*IT_0298;
    const ccomplex_t IT_0300 = IT_0064*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = m_C_1*IT_0301;
    const ccomplex_t IT_0303 = IT_0106*IT_0237;
    const ccomplex_t IT_0304 = IT_0099*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = m_C_1*IT_0305;
    const ccomplex_t IT_0307 = IT_0153*IT_0292;
    const ccomplex_t IT_0308 = IT_0109*IT_0110*IT_0307;
    const ccomplex_t IT_0309 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0310 = IT_0246*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0313 = IT_0250*IT_0312;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = IT_0311 + IT_0314;
    const ccomplex_t IT_0316 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0317 = IT_0246*IT_0316;
    const ccomplex_t IT_0318 = (-2)*IT_0317;
    const ccomplex_t IT_0319 = IT_0315 + IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = 0.5*IT_0320;
    const ccomplex_t IT_0322 = IT_0183*IT_0321;
    const ccomplex_t IT_0323 = IT_0000*IT_0035*IT_0322;
    const ccomplex_t IT_0324 = IT_0218*IT_0298;
    const ccomplex_t IT_0325 = IT_0062*IT_0064*IT_0324;
    const ccomplex_t IT_0326 = IT_0082*IT_0260;
    const ccomplex_t IT_0327 = IT_0074*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = m_b*IT_0328;
    const ccomplex_t IT_0330 = IT_0031*IT_0321;
    const ccomplex_t IT_0331 = IT_0035*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = m_b*IT_0332;
    const ccomplex_t IT_0334 = IT_0058*IT_0275;
    const ccomplex_t IT_0335 = IT_0060*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*IT_0335;
    const ccomplex_t IT_0337 = m_b*IT_0336;
    const ccomplex_t IT_0338 = IT_0239 + -IT_0245 + -IT_0262 + -IT_0277 +
       IT_0279 + IT_0296 + -IT_0302 + -IT_0306 + -IT_0308 + -IT_0323 + IT_0325 +
       IT_0329 + IT_0333 + IT_0337;
    const ccomplex_t IT_0339 = 6*s_13;
    const ccomplex_t IT_0340 = m_b*m_C_1;
    const ccomplex_t IT_0341 = 6*IT_0235;
    const ccomplex_t IT_0342 = -IT_0156 + -IT_0171 + -IT_0186 + -IT_0201 + 
      -IT_0244 + -IT_0301 + -IT_0305;
    const ccomplex_t IT_0343 = m_b*s_12;
    const ccomplex_t IT_0344 = IT_0211 + IT_0221 + IT_0233 + -IT_0295 + 
      -IT_0328 + -IT_0332 + -IT_0336;
    const ccomplex_t IT_0345 = s_23*m_C_1;
    const ccomplex_t IT_0346 = 6*conj(IT_0235);
    const ccomplex_t IT_0347 = 6*conj(IT_0338);
    const ccomplex_t IT_0348 = m_b*m_C_1*IT_0034;
    const ccomplex_t IT_0349 = s_13*IT_0034;
    const ccomplex_t IT_0350 = (-6)*IT_0349;
    const ccomplex_t IT_0351 = s_12*s_23;
    const ccomplex_t IT_0352 = 12*IT_0351;
    const ccomplex_t IT_0353 = IT_0350 + IT_0352;
    const ccomplex_t IT_0354 = 6*IT_0338;
    const ccomplex_t IT_0355 = 0.166666666666667*IT_0343;
    const ccomplex_t IT_0356 = 0.166666666666667*IT_0345;
    const ccomplex_t IT_0357 = (IT_0235*conj(IT_0235) + IT_0338*conj(IT_0338))
      *IT_0339 + conj(IT_0338)*IT_0340*IT_0341 + (IT_0338*IT_0340 + IT_0342
      *IT_0343 + IT_0344*IT_0345)*IT_0346 + (IT_0343*IT_0344 + IT_0342*IT_0345)
      *IT_0347 + 6*conj(IT_0344)*(IT_0342*IT_0348 + 0.166666666666667*IT_0344
      *IT_0353 + IT_0354*IT_0355 + IT_0341*IT_0356) + 6*conj(IT_0342)*(IT_0344
      *IT_0348 + 0.166666666666667*IT_0342*IT_0353 + IT_0341*IT_0355 + IT_0354
      *IT_0356);
    return create_ccomplex_return(IT_0357);
}

