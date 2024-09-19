#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_ss_L_to_anti_s_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_ss_L_to_anti_s_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0008 + 0.5*IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_s, 2);
    const ccomplex_t IT_0016 = pow(m_ss_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0034 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + 0.5*IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = IT_0033*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + 0.5*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0019 + IT_0032 + IT_0045 + IT_0058;
    const ccomplex_t IT_0060 = m_C_1*IT_0059;
    const ccomplex_t IT_0061 = m_s*m_C_1;
    const ccomplex_t IT_0062 = cos(theta_W);
    const ccomplex_t IT_0063 = cpow(IT_0062, -1);
    const ccomplex_t IT_0064 = N_B1*e_em;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W1*e_em;
    const ccomplex_t IT_0068 = IT_0004*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-3)*IT_0069);
    const ccomplex_t IT_0071 = 0.166666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0040*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0074 = IT_0033*IT_0072*IT_0073;
    const ccomplex_t IT_0075 = N_B2*e_em;
    const ccomplex_t IT_0076 = IT_0063*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_W2*e_em;
    const ccomplex_t IT_0079 = IT_0004*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + (-3)*IT_0080);
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0028*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0085 = IT_0030*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0089 = IT_0004*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.5)*IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0071*IT_0092;
    const ccomplex_t IT_0094 = IT_0033*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = m_C_1*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0098 = IT_0004*IT_0097;
    const ccomplex_t IT_0099 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0100 = IT_0004*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + (-0.5)*IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0082*IT_0103;
    const ccomplex_t IT_0105 = IT_0030*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = m_C_1*IT_0106;
    const ccomplex_t IT_0108 = N_B3*e_em;
    const ccomplex_t IT_0109 = IT_0063*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = N_W3*e_em;
    const ccomplex_t IT_0112 = IT_0004*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-3)*IT_0113);
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = IT_0054*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0118 = IT_0056*IT_0116*IT_0117;
    const ccomplex_t IT_0119 = N_B4*e_em;
    const ccomplex_t IT_0120 = IT_0063*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_W4*e_em;
    const ccomplex_t IT_0123 = IT_0004*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + (-3)*IT_0124);
    const ccomplex_t IT_0126 = 0.166666666666667*IT_0125;
    const ccomplex_t IT_0127 = IT_0013*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0129 = IT_0017*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0131 = IT_0004*IT_0130;
    const ccomplex_t IT_0132 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0133 = IT_0004*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + (-0.5)*IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0126*IT_0136;
    const ccomplex_t IT_0138 = IT_0017*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = m_C_1*IT_0139;
    const ccomplex_t IT_0141 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0142 = IT_0004*IT_0141;
    const ccomplex_t IT_0143 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0144 = IT_0004*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + (-0.5)*IT_0145);
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = IT_0115*IT_0147;
    const ccomplex_t IT_0149 = IT_0056*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = m_C_1*IT_0150;
    const ccomplex_t IT_0152 = IT_0074 + IT_0085 + -IT_0096 + -IT_0107 +
       IT_0118 + IT_0129 + -IT_0140 + -IT_0151;
    const ccomplex_t IT_0153 = IT_0021*IT_0103;
    const ccomplex_t IT_0154 = IT_0030*IT_0084*IT_0153;
    const ccomplex_t IT_0155 = IT_0047*IT_0147;
    const ccomplex_t IT_0156 = IT_0056*IT_0117*IT_0155;
    const ccomplex_t IT_0157 = IT_0006*IT_0136;
    const ccomplex_t IT_0158 = IT_0017*IT_0128*IT_0157;
    const ccomplex_t IT_0159 = IT_0042*IT_0092;
    const ccomplex_t IT_0160 = IT_0033*IT_0073*IT_0159;
    const ccomplex_t IT_0161 = -IT_0154 + -IT_0156 + -IT_0158 + -IT_0160;
    const ccomplex_t IT_0162 = pow(m_W, -2);
    const ccomplex_t IT_0163 = pow(m_W, 2);
    const ccomplex_t IT_0164 = s_13*IT_0163;
    const ccomplex_t IT_0165 = s_14*s_34;
    const ccomplex_t IT_0166 = (-2)*IT_0165;
    const ccomplex_t IT_0167 = IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = IT_0162*IT_0167;
    const ccomplex_t IT_0169 = (-6)*IT_0168;
    const ccomplex_t IT_0170 = 12*s_13;
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = (conj(IT_0060) + conj(IT_0161))*IT_0171;
    const ccomplex_t IT_0173 = (-18)*IT_0061*conj(IT_0152) + IT_0172;
    const ccomplex_t IT_0174 = s_14*s_34*IT_0162;
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = s_13 + IT_0175;
    const ccomplex_t IT_0177 = m_s*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0004;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = pow(m_C_1, 2);
    const ccomplex_t IT_0183 = cpow((-2)*s_13 + IT_0015 + IT_0182 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0184 = IT_0181*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0179*IT_0185;
    const ccomplex_t IT_0187 = cpow((-2)*s_13 + IT_0015 + IT_0182 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = V_ts*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0190 = IT_0004*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0004*U_st_01;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0192*IT_0194;
    const ccomplex_t IT_0196 = sin(beta);
    const ccomplex_t IT_0197 = cpow(IT_0196, -1);
    const ccomplex_t IT_0198 = IT_0004*IT_0197;
    const ccomplex_t IT_0199 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0200 = IT_0198*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0194*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0195 + IT_0205;
    const ccomplex_t IT_0207 = IT_0188*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0004;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_13 + IT_0015 + IT_0182 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*V_us*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0210*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = IT_0209*IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_13 + IT_0015 + IT_0182 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = V_ts*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0219 = IT_0004*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0004*U_st_00;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0226 = IT_0198*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = IT_0223*IT_0229;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = IT_0224 + IT_0231;
    const ccomplex_t IT_0233 = IT_0217*IT_0232;
    const ccomplex_t IT_0234 = IT_0186 + IT_0207 + IT_0215 + IT_0233;
    const ccomplex_t IT_0235 = m_C_1*IT_0176;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0183*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = IT_0179*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0242 = (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = IT_0210*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = IT_0209*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_00);
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0216*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0223*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_01);
    const ccomplex_t IT_0252 = (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = IT_0187*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0194*IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0019 + 2*IT_0032 + 2*IT_0045 + 2*IT_0058 +
       IT_0240 + IT_0245 + IT_0250 + IT_0255;
    const ccomplex_t IT_0257 = pow(s_34, 2);
    const ccomplex_t IT_0258 = IT_0162*IT_0257;
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = IT_0015 + IT_0259;
    const ccomplex_t IT_0261 = m_C_1*IT_0260;
    const ccomplex_t IT_0262 = -IT_0240 + -IT_0245 + -IT_0250 + -IT_0255;
    const ccomplex_t IT_0263 = 6*IT_0060;
    const ccomplex_t IT_0264 = pow(s_14, 2);
    const ccomplex_t IT_0265 = IT_0162*IT_0264;
    const ccomplex_t IT_0266 = -IT_0182;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = m_s*IT_0267;
    const ccomplex_t IT_0269 = (-2)*IT_0095 + (-2)*IT_0106 + (-2)*IT_0139 + (
      -2)*IT_0150 + -IT_0186 + -IT_0207 + -IT_0215 + -IT_0233;
    const ccomplex_t IT_0270 = s_14*IT_0163;
    const ccomplex_t IT_0271 = IT_0162*IT_0270;
    const ccomplex_t IT_0272 = 0.5*IT_0271;
    const ccomplex_t IT_0273 = s_14 + IT_0272;
    const ccomplex_t IT_0274 = m_s*IT_0273;
    const ccomplex_t IT_0275 = 2*IT_0274;
    const ccomplex_t IT_0276 = IT_0095 + IT_0106 + IT_0139 + IT_0150;
    const ccomplex_t IT_0277 = s_34*IT_0163;
    const ccomplex_t IT_0278 = IT_0162*IT_0277;
    const ccomplex_t IT_0279 = (-0.25)*IT_0278;
    const ccomplex_t IT_0280 = s_34 + IT_0279;
    const ccomplex_t IT_0281 = m_C_1*IT_0280;
    const ccomplex_t IT_0282 = (-4)*IT_0281;
    const ccomplex_t IT_0283 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0284 = (-6)*IT_0060;
    const ccomplex_t IT_0285 = s_14*IT_0015;
    const ccomplex_t IT_0286 = s_13*s_34;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = IT_0285 + IT_0287;
    const ccomplex_t IT_0289 = 6*IT_0288;
    const ccomplex_t IT_0290 = 6*conj(IT_0060);
    const ccomplex_t IT_0291 = 6*IT_0177;
    const ccomplex_t IT_0292 = 6*IT_0261;
    const ccomplex_t IT_0293 = s_13*IT_0176;
    const ccomplex_t IT_0294 = (-6)*IT_0293;
    const ccomplex_t IT_0295 = s_13*IT_0260;
    const ccomplex_t IT_0296 = (-6)*IT_0295;
    const ccomplex_t IT_0297 = IT_0061*IT_0176;
    const ccomplex_t IT_0298 = (-6)*IT_0297;
    const ccomplex_t IT_0299 = IT_0061*IT_0260;
    const ccomplex_t IT_0300 = (-6)*IT_0299;
    const ccomplex_t IT_0301 = pow(m_W, 4);
    const ccomplex_t IT_0302 = s_13*IT_0301;
    const ccomplex_t IT_0303 = IT_0162*IT_0302;
    const ccomplex_t IT_0304 = 6*IT_0303;
    const ccomplex_t IT_0305 = 24*IT_0165;
    const ccomplex_t IT_0306 = (-12)*IT_0164;
    const ccomplex_t IT_0307 = IT_0304 + IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = s_13*s_14;
    const ccomplex_t IT_0309 = s_34*IT_0182;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = 6*IT_0311;
    const ccomplex_t IT_0313 = IT_0061*IT_0163;
    const ccomplex_t IT_0314 = (-3)*IT_0313;
    const ccomplex_t IT_0315 = 6*IT_0314;
    const ccomplex_t IT_0316 = (-6)*conj(IT_0060);
    const ccomplex_t IT_0317 = (-6)*IT_0275;
    const ccomplex_t IT_0318 = (-6)*IT_0282;
    const ccomplex_t IT_0319 = 6*IT_0235;
    const ccomplex_t IT_0320 = (-6)*IT_0268;
    const ccomplex_t IT_0321 = s_13*IT_0267;
    const ccomplex_t IT_0322 = IT_0061*IT_0267;
    const ccomplex_t IT_0323 = 6*IT_0322;
    const ccomplex_t IT_0324 = 0.166666666666667*conj(IT_0234);
    const ccomplex_t IT_0325 = 0.166666666666667*conj(IT_0161);
    const ccomplex_t IT_0326 = 0.166666666666667*conj(IT_0262);
    const ccomplex_t IT_0327 = 0.166666666666667*conj(IT_0152);
    const ccomplex_t IT_0328 = 6*IT_0269*(0.166666666666667*conj(IT_0276)
      *IT_0312 + 0.166666666666667*IT_0268*IT_0316 + conj(IT_0269)*IT_0321 +
       0.166666666666667*conj(IT_0256)*IT_0323 + IT_0294*IT_0324 + IT_0320
      *IT_0325 + IT_0298*IT_0326 + IT_0319*IT_0327);
    const ccomplex_t IT_0329 = 6*IT_0256*(0.166666666666667*IT_0235*IT_0290 +
       0.166666666666667*conj(IT_0283)*IT_0312 + conj(IT_0256)*IT_0321 +
       0.166666666666667*conj(IT_0269)*IT_0323 + IT_0298*IT_0324 + IT_0319
      *IT_0325 + IT_0294*IT_0326 + IT_0320*IT_0327);
    const ccomplex_t IT_0330 = IT_0060*IT_0173 + (IT_0177*conj(IT_0234) +
       IT_0235*conj(IT_0256) + IT_0261*conj(IT_0262))*IT_0263 + (IT_0268*conj
      (IT_0269) + IT_0275*conj(IT_0276) + IT_0282*conj(IT_0283))*IT_0284 +
       IT_0262*(conj(IT_0283)*IT_0289 + IT_0261*IT_0290 + conj(IT_0152)*IT_0291 
      + conj(IT_0161)*IT_0292 + conj(IT_0256)*IT_0294 + conj(IT_0262)*IT_0296 +
       conj(IT_0269)*IT_0298 + conj(IT_0234)*IT_0300) + IT_0234*(conj(IT_0276)
      *IT_0289 + IT_0177*IT_0290 + conj(IT_0161)*IT_0291 + conj(IT_0152)*IT_0292
       + conj(IT_0269)*IT_0294 + conj(IT_0234)*IT_0296 + conj(IT_0256)*IT_0298 +
       conj(IT_0262)*IT_0300) + IT_0276*(conj(IT_0234)*IT_0289 + conj(IT_0276)
      *IT_0307 + conj(IT_0269)*IT_0312 + conj(IT_0283)*IT_0315 + IT_0275*IT_0316
       + conj(IT_0161)*IT_0317 + conj(IT_0152)*IT_0318) + IT_0283*(conj(IT_0262)
      *IT_0289 + conj(IT_0283)*IT_0307 + conj(IT_0256)*IT_0312 + conj(IT_0276)
      *IT_0315 + IT_0282*IT_0316 + conj(IT_0152)*IT_0317 + conj(IT_0161)*IT_0318
      ) + IT_0152*(IT_0061*((-18)*conj(IT_0060) + (-18)*conj(IT_0161)) + conj
      (IT_0152)*IT_0171 + conj(IT_0262)*IT_0291 + conj(IT_0234)*IT_0292 + conj
      (IT_0283)*IT_0317 + conj(IT_0276)*IT_0318 + conj(IT_0269)*IT_0319 + conj
      (IT_0256)*IT_0320) + IT_0161*(IT_0173 + conj(IT_0234)*IT_0291 + conj
      (IT_0262)*IT_0292 + conj(IT_0276)*IT_0317 + conj(IT_0283)*IT_0318 + conj
      (IT_0256)*IT_0319 + conj(IT_0269)*IT_0320) + IT_0328 + IT_0329;
    return create_ccomplex_return(IT_0330);
}

