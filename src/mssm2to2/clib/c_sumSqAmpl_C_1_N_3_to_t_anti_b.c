#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_3_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_3_to_t_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + 0.5*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_C_1, 2);
    const ccomplex_t IT_0013 = pow(m_N_3, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = s_13*s_24;
    const ccomplex_t IT_0018 = m_b*m_t*s_12;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0001*IT_0019*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = cos(theta_W);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = sin(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = m_t*N_u3*e_em*IT_0019*U_st_01;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0028 + (-0.75)*IT_0034);
    const ccomplex_t IT_0036 = 0.666666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0023*IT_0036;
    const ccomplex_t IT_0038 = pow(m_t, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0013 + IT_0038 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0001*IT_0019*IT_0021*conj(U_st_00);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0045 = IT_0025*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = m_t*N_u3*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0048 = IT_0031*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-0.75)*IT_0049);
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0013 + IT_0038 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = -IT_0041 + -IT_0055;
    const ccomplex_t IT_0057 = s_34*m_C_1*m_N_3;
    const ccomplex_t IT_0058 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0059 = IT_0025*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0062 = IT_0001*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = m_t*conj(N_u3)*e_em*IT_0019*U_st_11;
    const ccomplex_t IT_0065 = IT_0031*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0060 + 3*IT_0063 + 3
      *IT_0066);
    const ccomplex_t IT_0068 = 0.166666666666667*IT_0067;
    const ccomplex_t IT_0069 = m_t*V_tb*V_u1*e_em*IT_0019*conj(U_st_11);
    const ccomplex_t IT_0070 = IT_0031*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0075 = IT_0001*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0073 + IT_0078;
    const ccomplex_t IT_0080 = IT_0068*IT_0079;
    const ccomplex_t IT_0081 = IT_0039*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0084 = IT_0025*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_t*conj(N_u3)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0090 = IT_0031*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + 3*IT_0088 + 3
      *IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = m_t*V_tb*V_u1*e_em*IT_0019*conj(U_st_10);
    const ccomplex_t IT_0095 = IT_0031*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0100 = IT_0001*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0098 + IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = IT_0053*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0082 + IT_0107;
    const ccomplex_t IT_0109 = m_t*s_24*m_C_1;
    const ccomplex_t IT_0110 = pow(m_W, 2);
    const ccomplex_t IT_0111 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0110 
      + -reg_prop, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0019;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0001*IT_0029;
    const ccomplex_t IT_0115 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0116;
    const ccomplex_t IT_0118 = IT_0025*IT_0029;
    const ccomplex_t IT_0119 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0123 = IT_0114*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = IT_0117 + IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0113*IT_0128;
    const ccomplex_t IT_0130 = IT_0111*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0001*IT_0019*IT_0020*IT_0030;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0001*IT_0020;
    const ccomplex_t IT_0135 = IT_0115*IT_0134;
    const ccomplex_t IT_0136 = 2*IT_0135;
    const ccomplex_t IT_0137 = IT_0020*IT_0025;
    const ccomplex_t IT_0138 = IT_0119*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = IT_0122*IT_0134;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = IT_0139 + IT_0141;
    const ccomplex_t IT_0143 = IT_0136 + IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0133*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0131 + IT_0149;
    const ccomplex_t IT_0151 = IT_0036*IT_0079;
    const ccomplex_t IT_0152 = IT_0039*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0051*IT_0104;
    const ccomplex_t IT_0155 = IT_0053*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0153 + -IT_0156;
    const ccomplex_t IT_0158 = IT_0150 + IT_0157;
    const ccomplex_t IT_0159 = 12*conj(IT_0016);
    const ccomplex_t IT_0160 = cpow((-2)*s_13 + IT_0012 + IT_0038 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0161 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0162 = IT_0025*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = IT_0001*IT_0021;
    const ccomplex_t IT_0165 = m_b*conj(N_d3)*e_em*IT_0019*conj(U_sb_00);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0163 + 1.5*IT_0167);
    const ccomplex_t IT_0169 = (-0.333333333333333)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0019*IT_0030*U_sb_00;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = IT_0160*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_13 + IT_0012 + IT_0038 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0176 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0177 = IT_0025*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = m_b*conj(N_d3)*e_em*IT_0019*conj(U_sb_01);
    const ccomplex_t IT_0180 = IT_0164*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0178 + 1.5*IT_0181);
    const ccomplex_t IT_0183 = (-0.333333333333333)*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0019*IT_0030*U_sb_01;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0183*IT_0185;
    const ccomplex_t IT_0187 = IT_0175*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = -IT_0174 + -IT_0188;
    const ccomplex_t IT_0190 = m_b*s_13*m_N_3;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0019;
    const ccomplex_t IT_0192 = (-0.5)*IT_0191;
    const ccomplex_t IT_0193 = IT_0128*IT_0192;
    const ccomplex_t IT_0194 = IT_0111*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0001*IT_0019*IT_0021*IT_0029;
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = IT_0145*IT_0197;
    const ccomplex_t IT_0199 = IT_0147*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0195 + IT_0200;
    const ccomplex_t IT_0202 = m_t*s_14*m_N_3;
    const ccomplex_t IT_0203 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0204 = IT_0025*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0207 = IT_0001*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = m_b*N_d3*e_em*IT_0019*conj(U_sb_10);
    const ccomplex_t IT_0210 = IT_0164*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*(IT_0205 + (-3)*IT_0208 + 3
      *IT_0211);
    const ccomplex_t IT_0213 = 0.166666666666667*IT_0212;
    const ccomplex_t IT_0214 = IT_0171*IT_0213;
    const ccomplex_t IT_0215 = IT_0160*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0218 = IT_0025*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0221 = IT_0001*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = m_b*N_d3*e_em*IT_0019*conj(U_sb_11);
    const ccomplex_t IT_0224 = IT_0164*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0219 + (-3)*IT_0222 + 3
      *IT_0225);
    const ccomplex_t IT_0227 = 0.166666666666667*IT_0226;
    const ccomplex_t IT_0228 = IT_0185*IT_0227;
    const ccomplex_t IT_0229 = IT_0175*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = IT_0216 + IT_0230;
    const ccomplex_t IT_0232 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0233 = IT_0001*IT_0232;
    const ccomplex_t IT_0234 = m_b*conj(U_d1)*V_tb*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0235 = IT_0164*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0233 + (-0.5)*IT_0236);
    const ccomplex_t IT_0238 = IT_0213*IT_0237;
    const ccomplex_t IT_0239 = IT_0160*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0242 = IT_0001*IT_0241;
    const ccomplex_t IT_0243 = m_b*conj(U_d1)*V_tb*e_em*IT_0019*U_sb_11;
    const ccomplex_t IT_0244 = IT_0164*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + (-0.5)*IT_0245);
    const ccomplex_t IT_0247 = IT_0227*IT_0246;
    const ccomplex_t IT_0248 = IT_0175*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = -IT_0240 + -IT_0249;
    const ccomplex_t IT_0251 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0252 = IT_0134*IT_0251;
    const ccomplex_t IT_0253 = (-2)*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0255 = IT_0137*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0258 = IT_0134*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = IT_0256 + IT_0259;
    const ccomplex_t IT_0261 = IT_0253 + IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = IT_0113*IT_0263;
    const ccomplex_t IT_0265 = IT_0111*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0114*IT_0257;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = IT_0118*IT_0254;
    const ccomplex_t IT_0270 = 1.4142135623731*IT_0269;
    const ccomplex_t IT_0271 = IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = IT_0114*IT_0251;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = IT_0133*IT_0276;
    const ccomplex_t IT_0278 = IT_0147*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = -IT_0266 + -IT_0279;
    const ccomplex_t IT_0281 = m_b*s_23*m_C_1;
    const ccomplex_t IT_0282 = IT_0169*IT_0237;
    const ccomplex_t IT_0283 = IT_0160*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0183*IT_0246;
    const ccomplex_t IT_0286 = IT_0175*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = IT_0284 + IT_0287;
    const ccomplex_t IT_0289 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0290 = IT_0001*IT_0289;
    const ccomplex_t IT_0291 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0292 = IT_0001*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*(IT_0290 + (-0.5)*IT_0293);
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0010*IT_0295;
    const ccomplex_t IT_0297 = IT_0014*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = 12*conj(IT_0298);
    const ccomplex_t IT_0300 = s_14*s_23;
    const ccomplex_t IT_0301 = 24*IT_0016;
    const ccomplex_t IT_0302 = m_b*m_t*m_C_1*m_N_3;
    const ccomplex_t IT_0303 = 24*conj(IT_0016);
    const ccomplex_t IT_0304 = 24*conj(IT_0298);
    const ccomplex_t IT_0305 = IT_0043*IT_0093;
    const ccomplex_t IT_0306 = IT_0053*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = IT_0023*IT_0068;
    const ccomplex_t IT_0309 = IT_0039*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*IT_0309;
    const ccomplex_t IT_0311 = IT_0307 + IT_0310;
    const ccomplex_t IT_0312 = 6*IT_0302;
    const ccomplex_t IT_0313 = conj(IT_0231)*IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0190;
    const ccomplex_t IT_0315 = conj(IT_0250)*IT_0314;
    const ccomplex_t IT_0316 = IT_0159*IT_0190;
    const ccomplex_t IT_0317 = IT_0313 + IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = (-12)*conj(IT_0016);
    const ccomplex_t IT_0319 = IT_0192*IT_0263;
    const ccomplex_t IT_0320 = IT_0111*IT_0319;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = IT_0197*IT_0276;
    const ccomplex_t IT_0323 = IT_0147*IT_0322;
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*IT_0323;
    const ccomplex_t IT_0325 = -IT_0321 + -IT_0324;
    const ccomplex_t IT_0326 = (-12)*conj(IT_0298);
    const ccomplex_t IT_0327 = s_12*s_34;
    const ccomplex_t IT_0328 = -IT_0300;
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = IT_0017 + IT_0329;
    const ccomplex_t IT_0331 = 6*IT_0330;
    const ccomplex_t IT_0332 = 6*IT_0281;
    const ccomplex_t IT_0333 = 6*IT_0202;
    const ccomplex_t IT_0334 = 12*IT_0298;
    const ccomplex_t IT_0335 = 12*IT_0017;
    const ccomplex_t IT_0336 = 12*IT_0109;
    const ccomplex_t IT_0337 = 12*IT_0302;
    const ccomplex_t IT_0338 = 12*IT_0190;
    const ccomplex_t IT_0339 = (-6)*IT_0018;
    const ccomplex_t IT_0340 = -IT_0017;
    const ccomplex_t IT_0341 = IT_0329 + IT_0340;
    const ccomplex_t IT_0342 = (-6)*IT_0341;
    const ccomplex_t IT_0343 = (-6)*IT_0057;
    const ccomplex_t IT_0344 = IT_0190*IT_0301 + IT_0158*IT_0312 + IT_0325
      *IT_0331 + IT_0108*IT_0332 + IT_0056*IT_0333 + IT_0281*IT_0334 + IT_0288
      *IT_0335 + IT_0189*IT_0336 + IT_0231*IT_0337 + IT_0250*IT_0338 + IT_0280
      *IT_0339 + IT_0311*IT_0342 + IT_0201*IT_0343;
    const ccomplex_t IT_0345 = IT_0109*IT_0301 + IT_0150*IT_0331 + IT_0056
      *IT_0332 + IT_0108*IT_0333 + IT_0202*IT_0334 + IT_0231*IT_0335 + IT_0250
      *IT_0336 + IT_0288*IT_0337 + IT_0189*IT_0338 + IT_0201*IT_0339 + IT_0157
      *IT_0342 + IT_0280*IT_0343;
    const ccomplex_t IT_0346 = 6*IT_0018;
    const ccomplex_t IT_0347 = 6*IT_0109;
    const ccomplex_t IT_0348 = 6*IT_0057;
    const ccomplex_t IT_0349 = 12*IT_0016;
    const ccomplex_t IT_0350 = 12*IT_0300;
    const ccomplex_t IT_0351 = 12*IT_0281;
    const ccomplex_t IT_0352 = 12*IT_0202;
    const ccomplex_t IT_0353 = 24*IT_0298;
    const ccomplex_t IT_0354 = (-6)*IT_0281;
    const ccomplex_t IT_0355 = (-6)*IT_0202;
    const ccomplex_t IT_0356 = IT_0201*IT_0314 + IT_0288*IT_0332 + IT_0231
      *IT_0333 + IT_0056*IT_0337 + IT_0189*IT_0346 + IT_0280*IT_0347 + IT_0250
      *IT_0348 + IT_0057*IT_0349 + IT_0108*IT_0350 + IT_0311*IT_0351 + IT_0157
      *IT_0352 + IT_0300*IT_0353 + IT_0325*IT_0354 + IT_0150*IT_0355;
    const ccomplex_t IT_0357 = IT_0301*IT_0302 + IT_0158*IT_0314 + IT_0018
      *IT_0334 + IT_0189*IT_0335 + IT_0288*IT_0336 + IT_0250*IT_0337 + IT_0231
      *IT_0338 + IT_0108*IT_0346 + (IT_0311 + IT_0325)*IT_0347 + IT_0056*IT_0348
       + IT_0280*IT_0354 + IT_0201*IT_0355;
    const ccomplex_t IT_0358 = IT_0017*IT_0301 + IT_0057*IT_0334 + IT_0250
      *IT_0335 + IT_0231*IT_0336 + IT_0189*IT_0337 + IT_0288*IT_0338 + IT_0056
      *IT_0346 + IT_0158*IT_0347 + IT_0108*IT_0348 + IT_0201*IT_0354 + IT_0280
      *IT_0355;
    const ccomplex_t IT_0359 = IT_0280*IT_0314 + IT_0231*IT_0332 + IT_0288
      *IT_0333 + IT_0108*IT_0337 + IT_0250*IT_0346 + IT_0201*IT_0347 + IT_0189
      *IT_0348 + IT_0018*IT_0349 + IT_0056*IT_0350 + IT_0157*IT_0351 + IT_0311
      *IT_0352 + IT_0302*IT_0353 + IT_0150*IT_0354 + IT_0325*IT_0355;
    const ccomplex_t IT_0360 = 12*IT_0327;
    const ccomplex_t IT_0361 = (-12)*IT_0018;
    const ccomplex_t IT_0362 = (-12)*IT_0057;
    const ccomplex_t IT_0363 = (-12)*IT_0016;
    const ccomplex_t IT_0364 = IT_0108*IT_0314 + IT_0190*IT_0334 + IT_0280
      *IT_0337 + IT_0231*IT_0339 + IT_0288*IT_0343 + IT_0157*IT_0346 + IT_0056
      *IT_0347 + IT_0311*IT_0348 + IT_0250*IT_0354 + IT_0189*IT_0355 + IT_0201
      *IT_0360 + IT_0150*IT_0361 + IT_0325*IT_0362 + IT_0281*IT_0363;
    const ccomplex_t IT_0365 = IT_0056*IT_0314 + IT_0109*IT_0334 + IT_0201
      *IT_0337 + IT_0288*IT_0339 + IT_0231*IT_0343 + IT_0311*IT_0346 + IT_0108
      *IT_0347 + IT_0157*IT_0348 + IT_0189*IT_0354 + IT_0250*IT_0355 + IT_0280
      *IT_0360 + IT_0325*IT_0361 + IT_0150*IT_0362 + IT_0202*IT_0363;
    const ccomplex_t IT_0366 = IT_0250*IT_0314;
    const ccomplex_t IT_0367 = IT_0190*IT_0349;
    const ccomplex_t IT_0368 = IT_0300 + IT_0327 + IT_0340;
    const ccomplex_t IT_0369 = (-6)*IT_0368;
    const ccomplex_t IT_0370 = (-6)*IT_0302;
    const ccomplex_t IT_0371 = IT_0189*IT_0347;
    const ccomplex_t IT_0372 = IT_0231*IT_0312 + IT_0157*IT_0337 + IT_0288
      *IT_0342 + IT_0280*IT_0346 + IT_0201*IT_0348 + IT_0311*IT_0350 + IT_0108
      *IT_0351 + IT_0056*IT_0352 + IT_0281*IT_0353 + IT_0366 + IT_0367 + IT_0325
      *IT_0369 + IT_0150*IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = (-12)*IT_0298;
    const ccomplex_t IT_0374 = IT_0231*IT_0312 + IT_0288*IT_0331 + IT_0150
      *IT_0337 + IT_0108*IT_0354 + IT_0056*IT_0355 + IT_0325*IT_0360 + IT_0280
      *IT_0361 + IT_0201*IT_0362 + IT_0366 + IT_0367 + IT_0311*IT_0369 + IT_0157
      *IT_0370 + IT_0371 + IT_0281*IT_0373;
    const ccomplex_t IT_0375 = IT_0109*IT_0349;
    const ccomplex_t IT_0376 = IT_0189*IT_0314;
    const ccomplex_t IT_0377 = IT_0250*IT_0347;
    const ccomplex_t IT_0378 = IT_0288*IT_0312;
    const ccomplex_t IT_0379 = IT_0311*IT_0337 + IT_0231*IT_0342 + IT_0201
      *IT_0346 + IT_0280*IT_0348 + IT_0157*IT_0350 + IT_0056*IT_0351 + IT_0108
      *IT_0352 + IT_0202*IT_0353 + IT_0150*IT_0369 + IT_0325*IT_0370 + IT_0375 +
       IT_0376 + IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = IT_0231*IT_0331 + IT_0325*IT_0337 + IT_0056
      *IT_0354 + IT_0108*IT_0355 + IT_0150*IT_0360 + IT_0201*IT_0361 + IT_0280
      *IT_0362 + IT_0157*IT_0369 + IT_0311*IT_0370 + IT_0202*IT_0373 + IT_0375 +
       IT_0376 + IT_0377 + IT_0378;
    const ccomplex_t IT_0381 = 48*IT_0016*conj(IT_0016)*IT_0017 + (IT_0018
      *IT_0056 + IT_0057*IT_0108 + IT_0109*IT_0158)*IT_0159 + (IT_0018*IT_0189 +
       IT_0190*IT_0201 + IT_0202*IT_0231 + IT_0057*IT_0250 + IT_0109*IT_0280 +
       IT_0281*IT_0288)*IT_0299 + conj(IT_0298)*(48*IT_0298*IT_0300 + IT_0057
      *IT_0301) + (IT_0109*IT_0231 + IT_0017*IT_0250 + IT_0190*IT_0288 + IT_0057
      *IT_0298 + IT_0189*IT_0302)*IT_0303 + (IT_0108*IT_0300 + IT_0056*IT_0302)
      *IT_0304 + IT_0311*(IT_0281*IT_0304 + IT_0317) + IT_0201*IT_0281*IT_0318 +
       IT_0202*(IT_0157*IT_0304 + IT_0280*IT_0318 + IT_0150*IT_0326) + IT_0325*
      (IT_0317 + IT_0281*IT_0326) + conj(IT_0288)*IT_0344 + conj(IT_0231)
      *IT_0345 + conj(IT_0108)*IT_0356 + conj(IT_0189)*IT_0357 + conj(IT_0250)
      *IT_0358 + conj(IT_0056)*IT_0359 + conj(IT_0201)*IT_0364 + conj(IT_0280)
      *IT_0365 + conj(IT_0311)*IT_0372 + conj(IT_0325)*IT_0374 + conj(IT_0157)
      *IT_0379 + conj(IT_0150)*IT_0380;
    return create_ccomplex_return(IT_0381);
}

