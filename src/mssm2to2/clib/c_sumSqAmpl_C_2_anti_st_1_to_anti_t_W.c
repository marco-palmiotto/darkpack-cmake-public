#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_1_to_anti_t_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_1_to_anti_t_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_t*N_u1*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0013 + 3*IT_0016 + 3
      *IT_0023);
    const ccomplex_t IT_0025 = 0.166666666666667*IT_0024;
    const ccomplex_t IT_0026 = IT_0008*IT_0025;
    const ccomplex_t IT_0027 = pow(m_t, 2);
    const ccomplex_t IT_0028 = pow(m_st_1, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0027 + IT_0028 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0027 + IT_0028 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0033 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0034 = IT_0010*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = m_t*N_u2*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0040 = IT_0019*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0035 + 3*IT_0038 + 3
      *IT_0041);
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0045 = IT_0001*IT_0044;
    const ccomplex_t IT_0046 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + (-0.5)*IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0043*IT_0050;
    const ccomplex_t IT_0052 = IT_0032*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0057 = IT_0001*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.5)*IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0062 = IT_0010*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0065 = IT_0001*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = m_t*N_u3*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0068 = IT_0019*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0063 + 3*IT_0066 + 3
      *IT_0069);
    const ccomplex_t IT_0071 = 0.166666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0060*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0027 + IT_0028 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0077 = IT_0001*IT_0076;
    const ccomplex_t IT_0078 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0079 = IT_0001*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + (-0.5)*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0084 = IT_0010*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_t*N_u4*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0090 = IT_0019*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + 3*IT_0088 + 3
      *IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0082*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0027 + IT_0028 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-2)*IT_0031 + (-2)*IT_0053 + (-2)*IT_0075 + (
      -2)*IT_0097;
    const ccomplex_t IT_0099 = pow(m_W, -2);
    const ccomplex_t IT_0100 = pow(s_14, 2);
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = pow(m_C_2, 2);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = s_13*IT_0104;
    const ccomplex_t IT_0106 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0107 = IT_0001*IT_0106;
    const ccomplex_t IT_0108 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0109 = IT_0001*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + 0.5*IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = m_t*conj(N_u4)*e_em*IT_0020*conj(U_st_00);
    const ccomplex_t IT_0117 = IT_0019*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.75)*IT_0118);
    const ccomplex_t IT_0120 = 0.666666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0112*IT_0120;
    const ccomplex_t IT_0122 = IT_0095*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0125 = IT_0001*IT_0124;
    const ccomplex_t IT_0126 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + 0.5*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0132 = IT_0010*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = m_t*conj(N_u3)*e_em*IT_0020*conj(U_st_00);
    const ccomplex_t IT_0135 = IT_0019*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-0.75)*IT_0136);
    const ccomplex_t IT_0138 = 0.666666666666667*IT_0137;
    const ccomplex_t IT_0139 = IT_0130*IT_0138;
    const ccomplex_t IT_0140 = IT_0073*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0143 = IT_0001*IT_0142;
    const ccomplex_t IT_0144 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0145 = IT_0001*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + 0.5*IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0150 = IT_0010*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = m_t*conj(N_u1)*e_em*IT_0020*conj(U_st_00);
    const ccomplex_t IT_0153 = IT_0019*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + (-0.75)*IT_0154);
    const ccomplex_t IT_0156 = 0.666666666666667*IT_0155;
    const ccomplex_t IT_0157 = IT_0148*IT_0156;
    const ccomplex_t IT_0158 = IT_0029*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0161 = IT_0001*IT_0160;
    const ccomplex_t IT_0162 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0163 = IT_0001*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + 0.5*IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0168 = IT_0010*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = m_t*conj(N_u2)*e_em*IT_0020*conj(U_st_00);
    const ccomplex_t IT_0171 = IT_0019*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0169 + (-0.75)*IT_0172);
    const ccomplex_t IT_0174 = 0.666666666666667*IT_0173;
    const ccomplex_t IT_0175 = IT_0166*IT_0174;
    const ccomplex_t IT_0176 = IT_0032*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = 2*IT_0123 + 2*IT_0141 + 2*IT_0159 + 2*IT_0177;
    const ccomplex_t IT_0179 = m_t*m_C_2;
    const ccomplex_t IT_0180 = IT_0104*IT_0179;
    const ccomplex_t IT_0181 = IT_0043*IT_0166;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0183 = IT_0032*IT_0181*IT_0182;
    const ccomplex_t IT_0184 = pow(m_b, 2);
    const ccomplex_t IT_0185 = cpow((-2)*s_12 + -IT_0028 + -IT_0102 + IT_0184 
      + -reg_prop, -1);
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0187 = cos(beta);
    const ccomplex_t IT_0188 = cpow(IT_0187, -1);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0001*IT_0020*IT_0188*conj(U_st_00);
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0190*IT_0192;
    const ccomplex_t IT_0194 = IT_0185*IT_0186*IT_0193;
    const ccomplex_t IT_0195 = pow(m_d, 2);
    const ccomplex_t IT_0196 = cpow((-2)*s_12 + -IT_0028 + -IT_0102 + IT_0195 
      + -reg_prop, -1);
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0001*IT_0020*IT_0188*conj(U_st_00);
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0001;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = IT_0199*IT_0201;
    const ccomplex_t IT_0203 = IT_0196*IT_0197*IT_0202;
    const ccomplex_t IT_0204 = m_C_2*IT_0031;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0206 = IT_0025*IT_0148;
    const ccomplex_t IT_0207 = IT_0029*IT_0205*IT_0206;
    const ccomplex_t IT_0208 = IT_0071*IT_0130;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0210 = IT_0073*IT_0208*IT_0209;
    const ccomplex_t IT_0211 = pow(m_s, 2);
    const ccomplex_t IT_0212 = cpow((-2)*s_12 + -IT_0028 + -IT_0102 + IT_0211 
      + -reg_prop, -1);
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0001*IT_0020*IT_0188*conj(U_st_00);
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0001;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = IT_0215*IT_0217;
    const ccomplex_t IT_0219 = IT_0212*IT_0213*IT_0218;
    const ccomplex_t IT_0220 = m_C_2*IT_0053;
    const ccomplex_t IT_0221 = m_C_2*IT_0075;
    const ccomplex_t IT_0222 = m_C_2*IT_0097;
    const ccomplex_t IT_0223 = IT_0093*IT_0112;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0225 = IT_0095*IT_0223*IT_0224;
    const ccomplex_t IT_0226 = IT_0183 + IT_0194 + IT_0203 + -IT_0204 +
       IT_0207 + IT_0210 + IT_0219 + -IT_0220 + -IT_0221 + -IT_0222 + IT_0225;
    const ccomplex_t IT_0227 = s_14*s_34*IT_0099;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = s_13 + IT_0228;
    const ccomplex_t IT_0230 = m_C_2*IT_0229;
    const ccomplex_t IT_0231 = IT_0008*IT_0156;
    const ccomplex_t IT_0232 = IT_0029*IT_0205*IT_0231;
    const ccomplex_t IT_0233 = IT_0082*IT_0120;
    const ccomplex_t IT_0234 = IT_0095*IT_0224*IT_0233;
    const ccomplex_t IT_0235 = IT_0050*IT_0174;
    const ccomplex_t IT_0236 = IT_0032*IT_0182*IT_0235;
    const ccomplex_t IT_0237 = IT_0060*IT_0138;
    const ccomplex_t IT_0238 = IT_0073*IT_0209*IT_0237;
    const ccomplex_t IT_0239 = -IT_0232 + -IT_0234 + -IT_0236 + -IT_0238;
    const ccomplex_t IT_0240 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0241 = IT_0001*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = m_t*IT_0243;
    const ccomplex_t IT_0245 = m_t*V_tb*V_u2*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0246 = IT_0019*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = m_t*IT_0249;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0244 + IT_0251;
    const ccomplex_t IT_0253 = IT_0185*IT_0192;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0252*IT_0254;
    const ccomplex_t IT_0256 = IT_0212*IT_0217;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = V_ts*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0259 = IT_0001*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = m_t*IT_0261;
    const ccomplex_t IT_0263 = m_t*V_ts*V_u2*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0264 = IT_0019*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = m_t*IT_0267;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = IT_0262 + IT_0269;
    const ccomplex_t IT_0271 = IT_0257*IT_0270;
    const ccomplex_t IT_0272 = m_C_2*IT_0159;
    const ccomplex_t IT_0273 = m_C_2*IT_0177;
    const ccomplex_t IT_0274 = m_C_2*IT_0141;
    const ccomplex_t IT_0275 = m_C_2*IT_0123;
    const ccomplex_t IT_0276 = V_td*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0277 = IT_0001*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = m_t*IT_0279;
    const ccomplex_t IT_0281 = m_t*V_td*V_u2*e_em*IT_0020*conj(U_st_10);
    const ccomplex_t IT_0282 = IT_0019*IT_0281;
    const ccomplex_t IT_0283 = 1.4142135623731*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = (-0.5)*IT_0284;
    const ccomplex_t IT_0286 = m_t*IT_0285;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = IT_0280 + IT_0287;
    const ccomplex_t IT_0289 = IT_0196*IT_0201;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = IT_0288*IT_0290;
    const ccomplex_t IT_0292 = IT_0255 + IT_0271 + IT_0272 + IT_0273 + IT_0274
       + IT_0275 + IT_0291;
    const ccomplex_t IT_0293 = m_t*IT_0104;
    const ccomplex_t IT_0294 = (-6)*IT_0098;
    const ccomplex_t IT_0295 = IT_0293*IT_0294;
    const ccomplex_t IT_0296 = pow(m_W, 2);
    const ccomplex_t IT_0297 = s_13*IT_0296;
    const ccomplex_t IT_0298 = s_14*s_34;
    const ccomplex_t IT_0299 = (-2)*IT_0298;
    const ccomplex_t IT_0300 = IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0099*IT_0300;
    const ccomplex_t IT_0302 = (-6)*IT_0301;
    const ccomplex_t IT_0303 = 12*s_13;
    const ccomplex_t IT_0304 = IT_0302 + IT_0303;
    const ccomplex_t IT_0305 = IT_0239 + IT_0292;
    const ccomplex_t IT_0306 = IT_0304*IT_0305;
    const ccomplex_t IT_0307 = m_t*IT_0229;
    const ccomplex_t IT_0308 = -IT_0261;
    const ccomplex_t IT_0309 = IT_0267 + IT_0308;
    const ccomplex_t IT_0310 = IT_0257*IT_0309;
    const ccomplex_t IT_0311 = -IT_0243;
    const ccomplex_t IT_0312 = IT_0249 + IT_0311;
    const ccomplex_t IT_0313 = IT_0254*IT_0312;
    const ccomplex_t IT_0314 = -IT_0279;
    const ccomplex_t IT_0315 = IT_0285 + IT_0314;
    const ccomplex_t IT_0316 = IT_0290*IT_0315;
    const ccomplex_t IT_0317 = (-2)*IT_0310 + (-2)*IT_0313 + (-2)*IT_0316;
    const ccomplex_t IT_0318 = 6*IT_0317;
    const ccomplex_t IT_0319 = 6*IT_0178*IT_0230;
    const ccomplex_t IT_0320 = (-18)*IT_0179*IT_0226;
    const ccomplex_t IT_0321 = IT_0295 + IT_0306 + IT_0307*IT_0318 + IT_0319 +
       IT_0320;
    const ccomplex_t IT_0322 = s_13*IT_0229;
    const ccomplex_t IT_0323 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0324 = pow(s_34, 2);
    const ccomplex_t IT_0325 = IT_0099*IT_0324;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0027 + IT_0326;
    const ccomplex_t IT_0328 = m_C_2*IT_0327;
    const ccomplex_t IT_0329 = IT_0031 + IT_0053 + IT_0075 + IT_0097;
    const ccomplex_t IT_0330 = s_14*IT_0027;
    const ccomplex_t IT_0331 = s_13*s_34;
    const ccomplex_t IT_0332 = -IT_0331;
    const ccomplex_t IT_0333 = IT_0330 + IT_0332;
    const ccomplex_t IT_0334 = s_13*IT_0327;
    const ccomplex_t IT_0335 = (-6)*IT_0317;
    const ccomplex_t IT_0336 = -IT_0123 + -IT_0141 + -IT_0159 + -IT_0177;
    const ccomplex_t IT_0337 = s_13*s_14;
    const ccomplex_t IT_0338 = s_34*IT_0102;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0179*IT_0229;
    const ccomplex_t IT_0342 = (-6)*conj(IT_0317);
    const ccomplex_t IT_0343 = s_34*IT_0296;
    const ccomplex_t IT_0344 = IT_0099*IT_0343;
    const ccomplex_t IT_0345 = (-0.25)*IT_0344;
    const ccomplex_t IT_0346 = s_34 + IT_0345;
    const ccomplex_t IT_0347 = m_C_2*IT_0346;
    const ccomplex_t IT_0348 = (-4)*IT_0347;
    const ccomplex_t IT_0349 = (-6)*conj(IT_0239);
    const ccomplex_t IT_0350 = (-6)*conj(IT_0292);
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = s_14*IT_0296;
    const ccomplex_t IT_0353 = IT_0099*IT_0352;
    const ccomplex_t IT_0354 = 0.5*IT_0353;
    const ccomplex_t IT_0355 = s_14 + IT_0354;
    const ccomplex_t IT_0356 = m_t*IT_0355;
    const ccomplex_t IT_0357 = 2*IT_0356;
    const ccomplex_t IT_0358 = (-6)*conj(IT_0226);
    const ccomplex_t IT_0359 = -IT_0310 + -IT_0313 + -IT_0316;
    const ccomplex_t IT_0360 = (-0.25)*IT_0353;
    const ccomplex_t IT_0361 = s_14 + IT_0360;
    const ccomplex_t IT_0362 = m_t*IT_0361;
    const ccomplex_t IT_0363 = (-4)*IT_0362;
    const ccomplex_t IT_0364 = 0.5*IT_0344;
    const ccomplex_t IT_0365 = s_34 + IT_0364;
    const ccomplex_t IT_0366 = m_C_2*IT_0365;
    const ccomplex_t IT_0367 = 2*IT_0366;
    const ccomplex_t IT_0368 = IT_0179*IT_0296;
    const ccomplex_t IT_0369 = 3*IT_0368;
    const ccomplex_t IT_0370 = pow(m_W, 4);
    const ccomplex_t IT_0371 = s_13*IT_0370;
    const ccomplex_t IT_0372 = IT_0099*IT_0371;
    const ccomplex_t IT_0373 = 24*IT_0298;
    const ccomplex_t IT_0374 = (-12)*IT_0297;
    const ccomplex_t IT_0375 = 6*IT_0372 + IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = (-6)*IT_0239;
    const ccomplex_t IT_0377 = (-6)*IT_0292;
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = (-6)*IT_0226;
    const ccomplex_t IT_0380 = (-24)*IT_0298;
    const ccomplex_t IT_0381 = 6*IT_0372 + IT_0380;
    const ccomplex_t IT_0382 = 0.166666666666667*IT_0329;
    const ccomplex_t IT_0383 = 6*conj(IT_0359)*(0.166666666666667*IT_0333
      *IT_0335 + 0.166666666666667*IT_0294*IT_0340 + IT_0336*IT_0369 +
       0.166666666666667*IT_0359*IT_0375 + 0.166666666666667*IT_0363*IT_0378 +
       0.166666666666667*IT_0367*IT_0379 + IT_0381*IT_0382);
    const ccomplex_t IT_0384 = (-3)*IT_0368;
    const ccomplex_t IT_0385 = 0.166666666666667*IT_0357;
    const ccomplex_t IT_0386 = 0.166666666666667*IT_0348;
    const ccomplex_t IT_0387 = 6*conj(IT_0098)*(IT_0098*IT_0105 + IT_0178
      *IT_0180 + IT_0226*IT_0230) + (conj(IT_0239) + conj(IT_0292))*IT_0321 + 
      (IT_0293*IT_0305 + IT_0317*IT_0322)*IT_0323 + 6*conj(IT_0226)*(IT_0098
      *IT_0230 + 1./6*IT_0179*((-18)*IT_0239 + (-18)*IT_0292) + -IT_0178*IT_0293
       + 1./6*IT_0226*IT_0304 + 1./6*IT_0318*IT_0328) + 6*conj(IT_0317)*(IT_0305
      *IT_0307 + 0.166666666666667*IT_0294*IT_0322 + IT_0226*IT_0328 + IT_0329
      *IT_0333 + 0.166666666666667*IT_0334*IT_0335) + 6*conj(IT_0178)*(IT_0105
      *IT_0178 + IT_0098*IT_0180 + -IT_0226*IT_0293 + IT_0230*IT_0305 + IT_0336
      *IT_0340 + 0.166666666666667*IT_0335*IT_0341) + IT_0178*IT_0341*IT_0342 +
       IT_0336*(IT_0348*IT_0351 + IT_0357*IT_0358) + IT_0359*(IT_0323*IT_0340 +
       IT_0333*IT_0342 + IT_0351*IT_0363 + IT_0358*IT_0367) + IT_0383 + 6
      *IT_0329*(conj(IT_0098)*IT_0340 + IT_0351*IT_0385 + IT_0358*IT_0386) + 6
      *conj(IT_0336)*(IT_0178*IT_0340 + IT_0359*IT_0369 + 0.166666666666667
      *IT_0336*IT_0375 + IT_0329*IT_0384 + IT_0379*IT_0385 + IT_0378*IT_0386) +
       6*conj(IT_0329)*(IT_0317*IT_0333 + IT_0098*IT_0340 + 0.166666666666667
      *IT_0359*IT_0381 + IT_0375*IT_0382 + IT_0336*IT_0384 + IT_0378*IT_0385 +
       IT_0379*IT_0386);
    return create_ccomplex_return(IT_0387);
}

