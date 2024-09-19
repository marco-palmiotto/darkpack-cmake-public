#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sb_2_to_t_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sb_2_to_t_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_t*N_u3*e_em*IT_0012*U_st_01;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0008 + (-0.75)*IT_0015);
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0016;
    const ccomplex_t IT_0018 = pow(m_t, 2);
    const ccomplex_t IT_0019 = pow(m_N_3, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = pow(m_sb_2, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0018 + IT_0024 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0026 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + (-0.5)*IT_0030);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0010*IT_0012*U_sb_01;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0018 + IT_0024 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0039 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0042 = IT_0001*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + (-0.5)*IT_0043);
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0001*IT_0010*IT_0012*U_sb_01;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0038*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_01*conj(U_st_00);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0054 = IT_0005*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*N_u3*e_em*IT_0012*U_st_00;
    const ccomplex_t IT_0057 = IT_0011*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.75)*IT_0058);
    const ccomplex_t IT_0060 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = IT_0023 + 2*IT_0037 + 2*IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -2);
    const ccomplex_t IT_0067 = pow(s_14, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = -IT_0019;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0073 = cos(beta);
    const ccomplex_t IT_0074 = cpow(IT_0073, -1);
    const ccomplex_t IT_0075 = IT_0001*IT_0074;
    const ccomplex_t IT_0076 = m_b*conj(U_d1)*V_tb*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0078 + (-2)*IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0032*IT_0082;
    const ccomplex_t IT_0084 = IT_0025*IT_0072*IT_0083;
    const ccomplex_t IT_0085 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0086 = IT_0001*IT_0085;
    const ccomplex_t IT_0087 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + 0.5*IT_0089);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0082*IT_0091;
    const ccomplex_t IT_0093 = IT_0025*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = m_N_3*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0097 = m_b*conj(U_d2)*V_tb*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0098 = IT_0075*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0101 = IT_0001*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0099 + (-2)*IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0045*IT_0103;
    const ccomplex_t IT_0105 = IT_0038*IT_0096*IT_0104;
    const ccomplex_t IT_0106 = pow(m_b, 2);
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + -IT_0019 + -IT_0024 + IT_0106 
      + -reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0112 = IT_0005*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = m_b*N_d3*e_em*IT_0012*U_sb_01;
    const ccomplex_t IT_0115 = IT_0075*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(IT_0113 + 1.5*IT_0116);
    const ccomplex_t IT_0118 = (-0.333333333333333)*IT_0117;
    const ccomplex_t IT_0119 = IT_0110*IT_0118;
    const ccomplex_t IT_0120 = IT_0107*IT_0108*IT_0119;
    const ccomplex_t IT_0121 = -IT_0084 + IT_0095 + -IT_0105 + -IT_0120;
    const ccomplex_t IT_0122 = m_t*IT_0070;
    const ccomplex_t IT_0123 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0124 = IT_0005*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_t*conj(N_u3)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0130 = IT_0011*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0125 + 3*IT_0128 + 3
      *IT_0131);
    const ccomplex_t IT_0133 = 0.166666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0061*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0052*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0138 = IT_0005*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = m_t*conj(N_u3)*e_em*IT_0012*U_st_11;
    const ccomplex_t IT_0144 = IT_0011*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0139 + 3*IT_0142 + 3
      *IT_0145);
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0146;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0003*IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0094 + -IT_0136 + -IT_0150;
    const ccomplex_t IT_0152 = m_t*m_N_3;
    const ccomplex_t IT_0153 = IT_0070*IT_0152;
    const ccomplex_t IT_0154 = -IT_0037 + -IT_0050;
    const ccomplex_t IT_0155 = pow(m_W, 2);
    const ccomplex_t IT_0156 = IT_0152*IT_0155;
    const ccomplex_t IT_0157 = 3*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0159 = IT_0005*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0162 = IT_0001*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(N_d3)*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0165 = IT_0075*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + (-3)*IT_0163 + 3
      *IT_0166);
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0169 = IT_0110*IT_0168;
    const ccomplex_t IT_0170 = IT_0107*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = s_14*IT_0155;
    const ccomplex_t IT_0173 = IT_0066*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = s_14 + IT_0174;
    const ccomplex_t IT_0176 = m_t*IT_0175;
    const ccomplex_t IT_0177 = 2*IT_0176;
    const ccomplex_t IT_0178 = pow(m_W, 4);
    const ccomplex_t IT_0179 = s_13*IT_0178;
    const ccomplex_t IT_0180 = IT_0066*IT_0179;
    const ccomplex_t IT_0181 = s_14*s_34;
    const ccomplex_t IT_0182 = 24*IT_0181;
    const ccomplex_t IT_0183 = s_13*IT_0155;
    const ccomplex_t IT_0184 = (-12)*IT_0183;
    const ccomplex_t IT_0185 = 6*IT_0180 + IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = (-3)*IT_0156;
    const ccomplex_t IT_0187 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0188 = IT_0001*IT_0187;
    const ccomplex_t IT_0189 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0190 = IT_0001*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0188 + 0.5*IT_0191);
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0103*IT_0193;
    const ccomplex_t IT_0195 = IT_0038*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = -conj(IT_0094);
    const ccomplex_t IT_0199 = conj(IT_0197) + IT_0198;
    const ccomplex_t IT_0200 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = s_13 + IT_0201;
    const ccomplex_t IT_0203 = s_13*IT_0202;
    const ccomplex_t IT_0204 = -IT_0023 + -IT_0064;
    const ccomplex_t IT_0205 = (-6)*IT_0065;
    const ccomplex_t IT_0206 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0207 = IT_0152*IT_0202;
    const ccomplex_t IT_0208 = pow(s_34, 2);
    const ccomplex_t IT_0209 = IT_0066*IT_0208;
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = IT_0018 + IT_0210;
    const ccomplex_t IT_0212 = s_13*IT_0211;
    const ccomplex_t IT_0213 = m_t*IT_0202;
    const ccomplex_t IT_0214 = (-6)*IT_0204;
    const ccomplex_t IT_0215 = m_N_3*IT_0202;
    const ccomplex_t IT_0216 = s_34*IT_0155;
    const ccomplex_t IT_0217 = IT_0066*IT_0216;
    const ccomplex_t IT_0218 = (-0.25)*IT_0217;
    const ccomplex_t IT_0219 = s_34 + IT_0218;
    const ccomplex_t IT_0220 = m_N_3*IT_0219;
    const ccomplex_t IT_0221 = (-4)*IT_0220;
    const ccomplex_t IT_0222 = (-6)*IT_0121;
    const ccomplex_t IT_0223 = IT_0136 + IT_0150 + 2*IT_0171;
    const ccomplex_t IT_0224 = (-6)*IT_0151;
    const ccomplex_t IT_0225 = IT_0152*IT_0211;
    const ccomplex_t IT_0226 = (-6)*IT_0223;
    const ccomplex_t IT_0227 = m_N_3*IT_0211;
    const ccomplex_t IT_0228 = (-6)*IT_0180;
    const ccomplex_t IT_0229 = (-24)*IT_0181;
    const ccomplex_t IT_0230 = 12*IT_0183;
    const ccomplex_t IT_0231 = IT_0228 + IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = IT_0182 + IT_0228;
    const ccomplex_t IT_0233 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0234 = conj(IT_0154)*IT_0186;
    const ccomplex_t IT_0235 = m_N_3*IT_0196;
    const ccomplex_t IT_0236 = m_N_3*IT_0037;
    const ccomplex_t IT_0237 = m_N_3*IT_0050;
    const ccomplex_t IT_0238 = m_t*IT_0171;
    const ccomplex_t IT_0239 = -IT_0236 + -IT_0237 + -IT_0238;
    const ccomplex_t IT_0240 = IT_0034*IT_0091;
    const ccomplex_t IT_0241 = IT_0025*IT_0072*IT_0240;
    const ccomplex_t IT_0242 = IT_0047*IT_0193;
    const ccomplex_t IT_0243 = IT_0038*IT_0096*IT_0242;
    const ccomplex_t IT_0244 = IT_0241 + IT_0243;
    const ccomplex_t IT_0245 = (-2)*IT_0181;
    const ccomplex_t IT_0246 = IT_0183 + IT_0245;
    const ccomplex_t IT_0247 = IT_0066*IT_0246;
    const ccomplex_t IT_0248 = (-6)*IT_0247;
    const ccomplex_t IT_0249 = 12*s_13;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0196;
    const ccomplex_t IT_0252 = 0.166666666666667*IT_0224;
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 6*IT_0154;
    const ccomplex_t IT_0255 = 0.5*IT_0217;
    const ccomplex_t IT_0256 = s_34 + IT_0255;
    const ccomplex_t IT_0257 = m_N_3*IT_0256;
    const ccomplex_t IT_0258 = 2*IT_0257;
    const ccomplex_t IT_0259 = 6*IT_0171;
    const ccomplex_t IT_0260 = 6*IT_0094 + (-6)*IT_0197;
    const ccomplex_t IT_0261 = 6*IT_0065*IT_0122 + IT_0213*IT_0214 + IT_0226
      *IT_0227 + IT_0152*((-18)*IT_0239 + (-18)*IT_0244) + (IT_0121 + IT_0235)
      *IT_0250 + 6*IT_0215*IT_0253 + IT_0177*IT_0254 + IT_0258*IT_0259 + IT_0221
      *IT_0260;
    const ccomplex_t IT_0262 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0263 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0264 = 0.166666666666667*IT_0263;
    const ccomplex_t IT_0265 = conj(IT_0251) + IT_0264;
    const ccomplex_t IT_0266 = (-6)*conj(IT_0204);
    const ccomplex_t IT_0267 = 0.166666666666667*IT_0207;
    const ccomplex_t IT_0268 = 0.166666666666667*IT_0215;
    const ccomplex_t IT_0269 = s_14*IT_0018;
    const ccomplex_t IT_0270 = s_13*s_34;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = conj(IT_0171) + conj(IT_0197);
    const ccomplex_t IT_0274 = (-6)*conj(IT_0223);
    const ccomplex_t IT_0275 = IT_0171 + IT_0197;
    const ccomplex_t IT_0276 = 0.166666666666667*IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0180 + IT_0229;
    const ccomplex_t IT_0278 = 0.166666666666667*conj(IT_0094);
    const ccomplex_t IT_0279 = 0.166666666666667*conj(IT_0171);
    const ccomplex_t IT_0280 = IT_0239 + IT_0244;
    const ccomplex_t IT_0281 = (-0.25)*IT_0173;
    const ccomplex_t IT_0282 = s_14 + IT_0281;
    const ccomplex_t IT_0283 = m_t*IT_0282;
    const ccomplex_t IT_0284 = (-4)*IT_0283;
    const ccomplex_t IT_0285 = (-6)*conj(IT_0251);
    const ccomplex_t IT_0286 = 0.166666666666667*IT_0285;
    const ccomplex_t IT_0287 = conj(IT_0151) + IT_0286;
    const ccomplex_t IT_0288 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0289 = -IT_0199;
    const ccomplex_t IT_0290 = 0.166666666666667*IT_0213;
    const ccomplex_t IT_0291 = 0.166666666666667*IT_0227;
    const ccomplex_t IT_0292 = s_13*s_14;
    const ccomplex_t IT_0293 = s_34*IT_0019;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = (-6)*IT_0251;
    const ccomplex_t IT_0297 = 0.166666666666667*IT_0296;
    const ccomplex_t IT_0298 = IT_0205*IT_0215 + IT_0213*IT_0226 + IT_0214
      *IT_0227 + IT_0152*((-18)*IT_0121 + (-18)*IT_0235) + IT_0221*IT_0254 +
       IT_0177*IT_0260 + IT_0250*IT_0280 + IT_0259*IT_0284 + 6*IT_0122*(IT_0151 
      + IT_0297);
    const ccomplex_t IT_0299 = conj(IT_0065)*IT_0122;
    const ccomplex_t IT_0300 = conj(IT_0154)*IT_0177;
    const ccomplex_t IT_0301 = 6*IT_0065*(conj(IT_0065)*IT_0071 + conj(IT_0121
      )*IT_0122 + conj(IT_0151)*IT_0153) + 6*IT_0154*(IT_0157*conj(IT_0171) +
       conj(IT_0121)*IT_0177 + 0.166666666666667*conj(IT_0154)*IT_0185 + 
      -IT_0186*IT_0199) + IT_0203*(conj(IT_0204)*IT_0205 + IT_0204*IT_0206) + 
      (conj(IT_0151)*IT_0207 + conj(IT_0204)*IT_0212 + conj(IT_0121)*IT_0213)
      *IT_0214 + (conj(IT_0151)*IT_0215 + conj(IT_0197)*IT_0221)*IT_0222 + conj
      (IT_0223)*(IT_0205*IT_0207 + IT_0203*IT_0224 + IT_0214*IT_0225 + IT_0212
      *IT_0226 + IT_0222*IT_0227) + IT_0197*(IT_0185*conj(IT_0197) + conj
      (IT_0094)*IT_0231 + conj(IT_0171)*IT_0232 + IT_0221*IT_0233 + (-6)*IT_0234
      ) + conj(IT_0235)*IT_0261 + 6*conj(IT_0251)*(IT_0204*IT_0207 +
       0.166666666666667*IT_0071*IT_0224 + IT_0205*IT_0262) + 6*IT_0251*(conj
      (IT_0204)*IT_0207 + conj(IT_0121)*IT_0215 + IT_0203*conj(IT_0223) +
       IT_0206*IT_0262 + IT_0071*IT_0265) + 6*IT_0151*(IT_0071*conj(IT_0151) +
       conj(IT_0065)*IT_0153 + IT_0266*IT_0267 + IT_0233*IT_0268) + 6*IT_0272*
      (conj(IT_0154)*IT_0204 + IT_0154*conj(IT_0204) + conj(IT_0094)*IT_0223 +
       IT_0094*conj(IT_0223) + 0.166666666666667*IT_0226*IT_0273 + IT_0274
      *IT_0276) + IT_0259*(conj(IT_0154)*IT_0157 + 0.166666666666667*conj
      (IT_0197)*IT_0232 + conj(IT_0121)*IT_0258 + IT_0277*IT_0278 + IT_0185
      *IT_0279) + 6*IT_0094*(conj(IT_0121)*IT_0221 + 0.166666666666667*conj
      (IT_0197)*IT_0231 + IT_0234 + IT_0185*IT_0278 + IT_0277*IT_0279) + 6
      *IT_0223*(IT_0203*IT_0265 + 0.166666666666667*IT_0225*IT_0266 + IT_0206
      *IT_0267 + IT_0233*IT_0291) + 6*IT_0280*(conj(IT_0154)*IT_0221 + IT_0206
      *IT_0268 + conj(IT_0171)*IT_0284 + IT_0122*IT_0287 + 0.166666666666667
      *IT_0152*IT_0288 + IT_0177*IT_0289 + IT_0274*IT_0290 + IT_0266*IT_0291) +
       6*IT_0295*(conj(IT_0065)*IT_0154 + IT_0065*conj(IT_0154) + IT_0253
      *IT_0273 + conj(IT_0251)*IT_0275 + IT_0263*IT_0276 + IT_0094*IT_0287 +
       conj(IT_0094)*(IT_0151 + IT_0297)) + (conj(IT_0239) + conj(IT_0244))
      *IT_0298 + 6*IT_0121*(conj(IT_0094)*IT_0221 + 0.166666666666667*conj
      (IT_0121)*IT_0250 + IT_0215*conj(IT_0251) + conj(IT_0171)*IT_0258 +
       IT_0266*IT_0290 + IT_0299 + IT_0300) + 6*IT_0235*(0.166666666666667*conj
      (IT_0121)*IT_0250 + conj(IT_0171)*IT_0258 + IT_0215*IT_0265 + IT_0221
      *IT_0289 + IT_0266*IT_0290 + IT_0274*IT_0291 + IT_0299 + IT_0300);
    return create_ccomplex_return(IT_0301);
}

