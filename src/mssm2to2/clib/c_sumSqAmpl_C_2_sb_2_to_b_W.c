#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_b_W(
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -2);
    const ccomplex_t IT_0001 = pow(m_W, 4);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 6*IT_0003;
    const ccomplex_t IT_0005 = s_14*s_34;
    const ccomplex_t IT_0006 = (-24)*IT_0005;
    const ccomplex_t IT_0007 = IT_0004 + IT_0006;
    const ccomplex_t IT_0008 = pow(m_b, 2);
    const ccomplex_t IT_0009 = pow(m_sb_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + 0.5*IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = m_b*conj(N_d3)*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0024 + (-3)*IT_0027 + 3
      *IT_0034);
    const ccomplex_t IT_0036 = 0.166666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0019*IT_0036;
    const ccomplex_t IT_0038 = IT_0010*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0041 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0042 = IT_0012*IT_0041;
    const ccomplex_t IT_0043 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0044 = IT_0012*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.5*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0049 = IT_0021*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0052 = IT_0012*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_b*conj(N_d2)*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0055 = IT_0030*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + (-3)*IT_0053 + 3
      *IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = IT_0040*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0063 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0064 = IT_0012*IT_0063;
    const ccomplex_t IT_0065 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0066 = IT_0012*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + 0.5*IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0071 = IT_0021*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0074 = IT_0012*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = m_b*conj(N_d4)*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0077 = IT_0030*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0072 + (-3)*IT_0075 + 3
      *IT_0078);
    const ccomplex_t IT_0080 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = IT_0062*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0085 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0086 = IT_0012*IT_0085;
    const ccomplex_t IT_0087 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0088 = IT_0012*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + 0.5*IT_0089);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0093 = IT_0021*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0096 = IT_0012*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = m_b*conj(N_d1)*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0099 = IT_0030*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0094 + (-3)*IT_0097 + 3
      *IT_0100);
    const ccomplex_t IT_0102 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0103 = IT_0091*IT_0102;
    const ccomplex_t IT_0104 = IT_0084*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0061 + IT_0083 + IT_0105;
    const ccomplex_t IT_0107 = conj(IT_0039) + conj(IT_0106);
    const ccomplex_t IT_0108 = pow(m_c, 2);
    const ccomplex_t IT_0109 = pow(m_C_2, 2);
    const ccomplex_t IT_0110 = cpow((-2)*s_12 + -IT_0009 + IT_0108 + -IT_0109 
      + -reg_prop, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0012;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = m_b*conj(U_d2)*V_cb*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0114 = IT_0030*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = V_cb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0117 = IT_0012*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0115 + (-2)*IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0112*IT_0119;
    const ccomplex_t IT_0121 = IT_0110*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = cpow((-2)*s_12 + -IT_0009 + -IT_0109 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0124 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0125 = IT_0012*IT_0124;
    const ccomplex_t IT_0126 = m_b*conj(U_d2)*V_tb*e_em*IT_0031*U_sb_11;
    const ccomplex_t IT_0127 = IT_0030*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = IT_0123*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = pow(m_u, 2);
    const ccomplex_t IT_0136 = cpow((-2)*s_12 + -IT_0009 + -IT_0109 + IT_0135 
      + -reg_prop, -1);
    const ccomplex_t IT_0137 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0138 = IT_0012*IT_0137;
    const ccomplex_t IT_0139 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = IT_0012*IT_0029*IT_0137;
    const ccomplex_t IT_0142 = m_b*conj(U_d2)*e_em*IT_0031*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0143 = IT_0141*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0140 + (-0.5)*IT_0144);
    const ccomplex_t IT_0146 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012
      *IT_0146*V_ub_mod;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0145*IT_0148;
    const ccomplex_t IT_0150 = IT_0136*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0134 + IT_0151;
    const ccomplex_t IT_0153 = IT_0122 + IT_0152;
    const ccomplex_t IT_0154 = conj(IT_0122) + conj(IT_0152);
    const ccomplex_t IT_0155 = 24*IT_0005;
    const ccomplex_t IT_0156 = pow(m_W, 2);
    const ccomplex_t IT_0157 = s_13*IT_0156;
    const ccomplex_t IT_0158 = (-12)*IT_0157;
    const ccomplex_t IT_0159 = IT_0004 + IT_0155 + IT_0158;
    const ccomplex_t IT_0160 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0161 = IT_0012*IT_0160;
    const ccomplex_t IT_0162 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0163 = IT_0012*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + (-0.5)*IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0168 = IT_0021*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = m_b*N_d4*e_em*IT_0031*U_sb_01;
    const ccomplex_t IT_0171 = IT_0030*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0169 + 1.5*IT_0172);
    const ccomplex_t IT_0174 = (-0.333333333333333)*IT_0173;
    const ccomplex_t IT_0175 = IT_0166*IT_0174;
    const ccomplex_t IT_0176 = IT_0062*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0179 = IT_0012*IT_0178;
    const ccomplex_t IT_0180 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0181 = IT_0012*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0179 + (-0.5)*IT_0182);
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0186 = IT_0021*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = m_b*N_d2*e_em*IT_0031*U_sb_01;
    const ccomplex_t IT_0189 = IT_0030*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0187 + 1.5*IT_0190);
    const ccomplex_t IT_0192 = (-0.333333333333333)*IT_0191;
    const ccomplex_t IT_0193 = IT_0184*IT_0192;
    const ccomplex_t IT_0194 = IT_0040*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0197 = IT_0012*IT_0196;
    const ccomplex_t IT_0198 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0199 = IT_0012*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + (-0.5)*IT_0200);
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0204 = IT_0021*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = m_b*N_d1*e_em*IT_0031*U_sb_01;
    const ccomplex_t IT_0207 = IT_0030*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0205 + 1.5*IT_0208);
    const ccomplex_t IT_0210 = (-0.333333333333333)*IT_0209;
    const ccomplex_t IT_0211 = IT_0202*IT_0210;
    const ccomplex_t IT_0212 = IT_0084*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0215 = IT_0012*IT_0214;
    const ccomplex_t IT_0216 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0217 = IT_0012*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0215 + (-0.5)*IT_0218);
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0222 = IT_0021*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = m_b*N_d3*e_em*IT_0031*U_sb_01;
    const ccomplex_t IT_0225 = IT_0030*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + 1.5*IT_0226);
    const ccomplex_t IT_0228 = (-0.333333333333333)*IT_0227;
    const ccomplex_t IT_0229 = IT_0220*IT_0228;
    const ccomplex_t IT_0230 = IT_0010*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = -IT_0177 + -IT_0195 + -IT_0213 + -IT_0231;
    const ccomplex_t IT_0233 = s_34*IT_0156;
    const ccomplex_t IT_0234 = IT_0000*IT_0233;
    const ccomplex_t IT_0235 = (-0.25)*IT_0234;
    const ccomplex_t IT_0236 = s_34 + IT_0235;
    const ccomplex_t IT_0237 = m_C_2*IT_0236;
    const ccomplex_t IT_0238 = (-4)*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0240 = sin(beta);
    const ccomplex_t IT_0241 = cpow(IT_0240, -1);
    const ccomplex_t IT_0242 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0012*IT_0031*IT_0241*U_sb_01;
    const ccomplex_t IT_0243 = 0.5*IT_0242;
    const ccomplex_t IT_0244 = IT_0112*IT_0243;
    const ccomplex_t IT_0245 = IT_0110*IT_0239*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0012*IT_0031*IT_0241*U_sb_01;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = IT_0131*IT_0248;
    const ccomplex_t IT_0250 = IT_0123*IT_0246*IT_0249;
    const ccomplex_t IT_0251 = m_C_2*IT_0039;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0253 = IT_0080*IT_0166;
    const ccomplex_t IT_0254 = IT_0062*IT_0252*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0012*IT_0031*IT_0137*IT_0241*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0257 = 0.5*IT_0256;
    const ccomplex_t IT_0258 = IT_0148*IT_0257;
    const ccomplex_t IT_0259 = IT_0136*IT_0255*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0261 = IT_0102*IT_0202;
    const ccomplex_t IT_0262 = IT_0084*IT_0260*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0264 = IT_0058*IT_0184;
    const ccomplex_t IT_0265 = IT_0040*IT_0263*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0267 = IT_0036*IT_0220;
    const ccomplex_t IT_0268 = IT_0010*IT_0266*IT_0267;
    const ccomplex_t IT_0269 = m_C_2*IT_0105;
    const ccomplex_t IT_0270 = m_C_2*IT_0061;
    const ccomplex_t IT_0271 = m_C_2*IT_0083;
    const ccomplex_t IT_0272 = -IT_0245 + -IT_0250 + IT_0251 + -IT_0254 + 
      -IT_0259 + -IT_0262 + -IT_0265 + -IT_0268 + IT_0269 + IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = s_14*IT_0156;
    const ccomplex_t IT_0274 = IT_0000*IT_0273;
    const ccomplex_t IT_0275 = 0.5*IT_0274;
    const ccomplex_t IT_0276 = s_14 + IT_0275;
    const ccomplex_t IT_0277 = m_b*IT_0276;
    const ccomplex_t IT_0278 = 2*IT_0277;
    const ccomplex_t IT_0279 = IT_0019*IT_0228;
    const ccomplex_t IT_0280 = IT_0010*IT_0266*IT_0279;
    const ccomplex_t IT_0281 = m_C_2*IT_0231;
    const ccomplex_t IT_0282 = m_b*IT_0134;
    const ccomplex_t IT_0283 = m_b*IT_0151;
    const ccomplex_t IT_0284 = m_C_2*IT_0213;
    const ccomplex_t IT_0285 = m_C_2*IT_0195;
    const ccomplex_t IT_0286 = m_C_2*IT_0177;
    const ccomplex_t IT_0287 = IT_0091*IT_0210;
    const ccomplex_t IT_0288 = IT_0084*IT_0260*IT_0287;
    const ccomplex_t IT_0289 = IT_0047*IT_0192;
    const ccomplex_t IT_0290 = IT_0040*IT_0263*IT_0289;
    const ccomplex_t IT_0291 = IT_0069*IT_0174;
    const ccomplex_t IT_0292 = IT_0062*IT_0252*IT_0291;
    const ccomplex_t IT_0293 = m_b*IT_0122;
    const ccomplex_t IT_0294 = IT_0280 + -IT_0281 + -IT_0282 + -IT_0283 + 
      -IT_0284 + -IT_0285 + -IT_0286 + IT_0288 + IT_0290 + IT_0292 + -IT_0293;
    const ccomplex_t IT_0295 = m_b*m_C_2;
    const ccomplex_t IT_0296 = IT_0156*IT_0295;
    const ccomplex_t IT_0297 = (-3)*IT_0296;
    const ccomplex_t IT_0298 = (-2)*IT_0039 + (-2)*IT_0061 + (-2)*IT_0083 + (
      -2)*IT_0105;
    const ccomplex_t IT_0299 = s_13*s_14;
    const ccomplex_t IT_0300 = s_34*IT_0109;
    const ccomplex_t IT_0301 = -IT_0300;
    const ccomplex_t IT_0302 = IT_0299 + IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0122 + 2*IT_0134 + 2*IT_0151;
    const ccomplex_t IT_0304 = s_14*IT_0008;
    const ccomplex_t IT_0305 = s_13*s_34;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = 6*conj(IT_0039);
    const ccomplex_t IT_0309 = 0.5*IT_0234;
    const ccomplex_t IT_0310 = s_34 + IT_0309;
    const ccomplex_t IT_0311 = m_C_2*IT_0310;
    const ccomplex_t IT_0312 = 2*IT_0311;
    const ccomplex_t IT_0313 = (-0.25)*IT_0274;
    const ccomplex_t IT_0314 = s_14 + IT_0313;
    const ccomplex_t IT_0315 = m_b*IT_0314;
    const ccomplex_t IT_0316 = (-4)*IT_0315;
    const ccomplex_t IT_0317 = 3*IT_0296;
    const ccomplex_t IT_0318 = 6*conj(IT_0122);
    const ccomplex_t IT_0319 = 6*IT_0302;
    const ccomplex_t IT_0320 = 6*IT_0307;
    const ccomplex_t IT_0321 = 6*IT_0297;
    const ccomplex_t IT_0322 = 6*IT_0238;
    const ccomplex_t IT_0323 = 6*IT_0278;
    const ccomplex_t IT_0324 = 6*IT_0039;
    const ccomplex_t IT_0325 = 6*IT_0122;
    const ccomplex_t IT_0326 = 2*IT_0177 + 2*IT_0195 + 2*IT_0213 + 2*IT_0231;
    const ccomplex_t IT_0327 = 6*IT_0317;
    const ccomplex_t IT_0328 = (-6)*conj(IT_0122);
    const ccomplex_t IT_0329 = 6*IT_0312;
    const ccomplex_t IT_0330 = 6*IT_0316;
    const ccomplex_t IT_0331 = (-6)*IT_0302;
    const ccomplex_t IT_0332 = (-6)*IT_0307;
    const ccomplex_t IT_0333 = pow(s_14, 2);
    const ccomplex_t IT_0334 = IT_0000*IT_0333;
    const ccomplex_t IT_0335 = -IT_0109;
    const ccomplex_t IT_0336 = IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = s_13*IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0337;
    const ccomplex_t IT_0339 = m_b*IT_0336;
    const ccomplex_t IT_0340 = 6*IT_0339;
    const ccomplex_t IT_0341 = IT_0295*IT_0336;
    const ccomplex_t IT_0342 = 6*IT_0341;
    const ccomplex_t IT_0343 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = s_13 + IT_0344;
    const ccomplex_t IT_0346 = m_C_2*IT_0345;
    const ccomplex_t IT_0347 = (-6)*IT_0346;
    const ccomplex_t IT_0348 = IT_0295*IT_0345;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = 12*s_13;
    const ccomplex_t IT_0351 = (-2)*IT_0005;
    const ccomplex_t IT_0352 = IT_0157 + IT_0351;
    const ccomplex_t IT_0353 = IT_0000*IT_0352;
    const ccomplex_t IT_0354 = (-6)*IT_0353;
    const ccomplex_t IT_0355 = IT_0350 + IT_0354;
    const ccomplex_t IT_0356 = (-18)*IT_0295;
    const ccomplex_t IT_0357 = pow(s_34, 2);
    const ccomplex_t IT_0358 = IT_0000*IT_0357;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = IT_0008 + IT_0359;
    const ccomplex_t IT_0361 = m_C_2*IT_0360;
    const ccomplex_t IT_0362 = (-6)*IT_0361;
    const ccomplex_t IT_0363 = conj(IT_0272)*(IT_0106*IT_0322 + IT_0232
      *IT_0323 + IT_0238*IT_0324 + IT_0312*IT_0325 + IT_0152*IT_0329 + IT_0326
      *IT_0340 + IT_0298*IT_0347 + IT_0272*IT_0355 + IT_0294*IT_0356 + IT_0303
      *IT_0362);
    const ccomplex_t IT_0364 = m_b*IT_0345;
    const ccomplex_t IT_0365 = (-6)*IT_0364;
    const ccomplex_t IT_0366 = conj(IT_0294)*(IT_0232*IT_0322 + IT_0106
      *IT_0323 + IT_0278*IT_0324 + IT_0316*IT_0325 + IT_0152*IT_0330 + IT_0298
      *IT_0340 + IT_0326*IT_0347 + IT_0294*IT_0355 + IT_0272*IT_0356 + IT_0303
      *IT_0365);
    const ccomplex_t IT_0367 = IT_0122 + (-0.166666666666667)*IT_0324;
    const ccomplex_t IT_0368 = s_13*IT_0345;
    const ccomplex_t IT_0369 = (-6)*IT_0368;
    const ccomplex_t IT_0370 = (-0.166666666666667)*IT_0272;
    const ccomplex_t IT_0371 = (-0.166666666666667)*IT_0298;
    const ccomplex_t IT_0372 = (-0.166666666666667)*IT_0294;
    const ccomplex_t IT_0373 = (-0.166666666666667)*IT_0152;
    const ccomplex_t IT_0374 = (-0.166666666666667)*IT_0106;
    const ccomplex_t IT_0375 = s_13*IT_0360;
    const ccomplex_t IT_0376 = IT_0007*(IT_0107*IT_0153 + (IT_0039 + IT_0106)
      *IT_0154) + IT_0159*((IT_0039 + IT_0106)*IT_0107 + IT_0153*IT_0154 +
       IT_0232*conj(IT_0232)) + (IT_0238*IT_0272 + IT_0278*IT_0294 + IT_0232
      *IT_0297 + IT_0298*IT_0302 + IT_0303*IT_0307)*IT_0308 + (IT_0272*IT_0312 +
       IT_0294*IT_0316 + IT_0232*IT_0317)*IT_0318 + conj(IT_0106)*(IT_0298
      *IT_0319 + IT_0303*IT_0320 + IT_0232*IT_0321 + IT_0272*IT_0322 + IT_0294
      *IT_0323) + conj(IT_0232)*(IT_0106*IT_0321 + IT_0294*IT_0322 + IT_0272
      *IT_0323 + IT_0297*IT_0324 + IT_0317*IT_0325 + IT_0319*IT_0326 + IT_0152
      *IT_0327) + (IT_0298*IT_0302 + IT_0303*IT_0307)*IT_0328 + conj(IT_0152)*
      (IT_0232*IT_0327 + IT_0272*IT_0329 + IT_0294*IT_0330 + IT_0298*IT_0331 +
       IT_0303*IT_0332) + conj(IT_0326)*(IT_0232*IT_0319 + IT_0326*IT_0338 +
       IT_0272*IT_0340 + IT_0298*IT_0342 + IT_0294*IT_0347 + IT_0303*IT_0349) +
       IT_0363 + IT_0366 + (-6)*conj(IT_0298)*((-0.166666666666667)*IT_0326
      *IT_0342 + IT_0302*IT_0367 + (-0.166666666666667)*IT_0303*IT_0369 +
       IT_0347*IT_0370 + IT_0338*IT_0371 + IT_0340*IT_0372 + IT_0331*IT_0373 +
       IT_0319*IT_0374) + (-6)*conj(IT_0303)*((-0.166666666666667)*IT_0326
      *IT_0349 + IT_0307*IT_0367 + IT_0362*IT_0370 + IT_0369*IT_0371 + IT_0365
      *IT_0372 + IT_0332*IT_0373 + IT_0320*IT_0374 + IT_0303*IT_0375);
    return create_ccomplex_return(IT_0376);
}

