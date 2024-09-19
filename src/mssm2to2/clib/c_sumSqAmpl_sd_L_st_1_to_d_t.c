#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_st_1_to_d_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_st_1_to_d_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_st_1, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_d, 2);
    const ccomplex_t IT_0007 = pow(m_sd_L, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = sin(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = m_t*conj(N_u2)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0020 + 3*IT_0023 + 3
      *IT_0029);
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = IT_0015*IT_0031;
    const ccomplex_t IT_0033 = IT_0008*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = m_t*conj(N_u3)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0045 = IT_0026*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0040 + 3*IT_0043 + 3
      *IT_0046);
    const ccomplex_t IT_0048 = 0.166666666666667*IT_0047;
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = IT_0035*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0056 = IT_0017*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0059 = IT_0013*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = m_t*conj(N_u4)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0062 = IT_0026*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0057 + 3*IT_0060 + 3
      *IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = IT_0052*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0073 = IT_0017*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0076 = IT_0013*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_t*conj(N_u1)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0079 = IT_0026*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0074 + 3*IT_0077 + 3
      *IT_0080);
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0071*IT_0082;
    const ccomplex_t IT_0084 = IT_0069*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0034 + -IT_0051 + -IT_0068 + -IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm2)*IT_0013;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0009*IT_0011*IT_0013*U_st_00;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0087*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm1)*IT_0013;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0009*IT_0011*IT_0013*U_st_00;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = IT_0094*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0093 + IT_0100;
    const ccomplex_t IT_0102 = conj(N_B2)*e_em;
    const ccomplex_t IT_0103 = IT_0017*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W2)*e_em;
    const ccomplex_t IT_0106 = IT_0013*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + (-3)*IT_0107);
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0110 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0111 = IT_0017*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = m_t*N_u2*e_em*IT_0009*U_st_00;
    const ccomplex_t IT_0114 = IT_0026*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + (-0.75)*IT_0115);
    const ccomplex_t IT_0117 = 0.666666666666667*IT_0116;
    const ccomplex_t IT_0118 = IT_0109*IT_0117;
    const ccomplex_t IT_0119 = IT_0008*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B4)*e_em;
    const ccomplex_t IT_0122 = IT_0017*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W4)*e_em;
    const ccomplex_t IT_0125 = IT_0013*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + (-3)*IT_0126);
    const ccomplex_t IT_0128 = 0.166666666666667*IT_0127;
    const ccomplex_t IT_0129 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0130 = IT_0017*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = m_t*N_u4*e_em*IT_0009*U_st_00;
    const ccomplex_t IT_0133 = IT_0026*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + (-0.75)*IT_0134);
    const ccomplex_t IT_0136 = 0.666666666666667*IT_0135;
    const ccomplex_t IT_0137 = IT_0128*IT_0136;
    const ccomplex_t IT_0138 = IT_0052*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B3)*e_em;
    const ccomplex_t IT_0141 = IT_0017*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W3)*e_em;
    const ccomplex_t IT_0144 = IT_0013*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + (-3)*IT_0145);
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0146;
    const ccomplex_t IT_0148 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0149 = IT_0017*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = m_t*N_u3*e_em*IT_0009*U_st_00;
    const ccomplex_t IT_0152 = IT_0026*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + (-0.75)*IT_0153);
    const ccomplex_t IT_0155 = 0.666666666666667*IT_0154;
    const ccomplex_t IT_0156 = IT_0147*IT_0155;
    const ccomplex_t IT_0157 = IT_0035*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = pow(m_sG, 2);
    const ccomplex_t IT_0160 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -IT_0159 +
       reg_prop, -1);
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = IT_0160*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B1)*e_em;
    const ccomplex_t IT_0167 = IT_0017*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W1)*e_em;
    const ccomplex_t IT_0170 = IT_0013*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + (-3)*IT_0171);
    const ccomplex_t IT_0173 = 0.166666666666667*IT_0172;
    const ccomplex_t IT_0174 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0175 = IT_0017*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = m_t*N_u1*e_em*IT_0009*U_st_00;
    const ccomplex_t IT_0178 = IT_0026*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0176 + (-0.75)*IT_0179);
    const ccomplex_t IT_0181 = 0.666666666666667*IT_0180;
    const ccomplex_t IT_0182 = IT_0173*IT_0181;
    const ccomplex_t IT_0183 = IT_0069*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = -IT_0120 + -IT_0139 + -IT_0158 +
       0.166666666666667*IT_0165 + -IT_0184;
    const ccomplex_t IT_0186 = (-0.5)*IT_0165;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u1
      *e_em*IT_0009*IT_0013*IT_0025;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = conj(V_td)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0190 = IT_0013*IT_0189;
    const ccomplex_t IT_0191 = m_t*conj(V_td)*conj(V_u1)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0192 = IT_0026*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0190 + (-0.5)*IT_0193);
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = IT_0188*IT_0195;
    const ccomplex_t IT_0197 = IT_0094*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u2
      *e_em*IT_0009*IT_0013*IT_0025;
    const ccomplex_t IT_0200 = 0.5*IT_0199;
    const ccomplex_t IT_0201 = conj(V_td)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0202 = IT_0013*IT_0201;
    const ccomplex_t IT_0203 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0009*U_st_10;
    const ccomplex_t IT_0204 = IT_0026*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0202 + (-0.5)*IT_0205);
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0200*IT_0207;
    const ccomplex_t IT_0209 = IT_0087*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = -IT_0198 + -IT_0210;
    const ccomplex_t IT_0212 = IT_0186 + IT_0211;
    const ccomplex_t IT_0213 = conj(IT_0186) + conj(IT_0211);
    const ccomplex_t IT_0214 = (-18)*IT_0001;
    const ccomplex_t IT_0215 = 36*IT_0003;
    const ccomplex_t IT_0216 = IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = m_d*m_t*IT_0000;
    const ccomplex_t IT_0218 = (-18)*IT_0101;
    const ccomplex_t IT_0219 = (-6)*IT_0086;
    const ccomplex_t IT_0220 = IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = (-18)*conj(IT_0101);
    const ccomplex_t IT_0222 = (-6)*conj(IT_0086);
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = 3*(conj(IT_0086) + 0.333333333333333*conj
      (IT_0101))*IT_0185 + 3*(IT_0086 + 0.333333333333333*IT_0101)*conj(IT_0185)
       + (conj(IT_0086) + 3*conj(IT_0101))*IT_0186 + (IT_0086 + 3*IT_0101)*conj
      (IT_0186) + (-0.166666666666667)*conj(IT_0211)*IT_0220 + (
      -0.166666666666667)*IT_0211*IT_0223;
    const ccomplex_t IT_0225 = m_t*s_23;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0227 = IT_0088*IT_0207;
    const ccomplex_t IT_0228 = IT_0087*IT_0226*IT_0227;
    const ccomplex_t IT_0229 = m_t*IT_0165;
    const ccomplex_t IT_0230 = -IT_0228 + (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0233 = IT_0161*IT_0232;
    const ccomplex_t IT_0234 = IT_0160*IT_0231*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0236 = IT_0095*IT_0195;
    const ccomplex_t IT_0237 = IT_0094*IT_0235*IT_0236;
    const ccomplex_t IT_0238 = m_d*IT_0100;
    const ccomplex_t IT_0239 = m_d*IT_0093;
    const ccomplex_t IT_0240 = (-0.5)*IT_0234 + IT_0237 + IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = m_t*IT_0184;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0243 = IT_0082*IT_0173;
    const ccomplex_t IT_0244 = IT_0069*IT_0242*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0246 = IT_0031*IT_0109;
    const ccomplex_t IT_0247 = IT_0008*IT_0245*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0249 = IT_0065*IT_0128;
    const ccomplex_t IT_0250 = IT_0052*IT_0248*IT_0249;
    const ccomplex_t IT_0251 = m_t*IT_0120;
    const ccomplex_t IT_0252 = m_t*IT_0158;
    const ccomplex_t IT_0253 = m_t*IT_0139;
    const ccomplex_t IT_0254 = 0.166666666666667*IT_0229 + -IT_0241 + IT_0244 
      + IT_0247 + IT_0250 + -IT_0251 + -IT_0252 + -IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0256 = IT_0048*IT_0147;
    const ccomplex_t IT_0257 = IT_0035*IT_0255*IT_0256;
    const ccomplex_t IT_0258 = 0.166666666666667*IT_0234 + -IT_0257;
    const ccomplex_t IT_0259 = IT_0230 + -IT_0240 + 3*IT_0254 + (-3)*IT_0258;
    const ccomplex_t IT_0260 = (-3)*conj(IT_0258);
    const ccomplex_t IT_0261 = -conj(IT_0240);
    const ccomplex_t IT_0262 = conj(IT_0230) + 3*conj(IT_0254) + IT_0260 +
       IT_0261;
    const ccomplex_t IT_0263 = IT_0097*IT_0188;
    const ccomplex_t IT_0264 = IT_0094*IT_0235*IT_0263;
    const ccomplex_t IT_0265 = m_d*IT_0198;
    const ccomplex_t IT_0266 = m_d*IT_0210;
    const ccomplex_t IT_0267 = -IT_0264 + -IT_0265 + -IT_0266;
    const ccomplex_t IT_0268 = 6*IT_0267;
    const ccomplex_t IT_0269 = IT_0015*IT_0117;
    const ccomplex_t IT_0270 = IT_0008*IT_0245*IT_0269;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = 18*IT_0271;
    const ccomplex_t IT_0273 = IT_0054*IT_0136;
    const ccomplex_t IT_0274 = IT_0052*IT_0248*IT_0273;
    const ccomplex_t IT_0275 = IT_0071*IT_0181;
    const ccomplex_t IT_0276 = IT_0069*IT_0242*IT_0275;
    const ccomplex_t IT_0277 = IT_0037*IT_0155;
    const ccomplex_t IT_0278 = IT_0035*IT_0255*IT_0277;
    const ccomplex_t IT_0279 = m_t*IT_0085;
    const ccomplex_t IT_0280 = m_t*IT_0034;
    const ccomplex_t IT_0281 = m_t*IT_0051;
    const ccomplex_t IT_0282 = m_t*IT_0068;
    const ccomplex_t IT_0283 = IT_0274 + IT_0276 + IT_0278 + -IT_0279 + 
      -IT_0280 + -IT_0281 + -IT_0282;
    const ccomplex_t IT_0284 = (-18)*IT_0283;
    const ccomplex_t IT_0285 = IT_0268 + IT_0272 + IT_0284;
    const ccomplex_t IT_0286 = 6*conj(IT_0267);
    const ccomplex_t IT_0287 = 18*conj(IT_0271);
    const ccomplex_t IT_0288 = (-18)*conj(IT_0283);
    const ccomplex_t IT_0289 = IT_0286 + IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = 6*IT_0271;
    const ccomplex_t IT_0291 = 18*IT_0267;
    const ccomplex_t IT_0292 = (-6)*IT_0283;
    const ccomplex_t IT_0293 = IT_0290 + IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = 6*conj(IT_0271);
    const ccomplex_t IT_0295 = 18*conj(IT_0267);
    const ccomplex_t IT_0296 = (-6)*conj(IT_0283);
    const ccomplex_t IT_0297 = IT_0294 + IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = 6*IT_0258;
    const ccomplex_t IT_0299 = 18*IT_0240;
    const ccomplex_t IT_0300 = (-18)*IT_0230;
    const ccomplex_t IT_0301 = (-6)*IT_0254;
    const ccomplex_t IT_0302 = IT_0298 + IT_0299 + IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = 6*conj(IT_0258);
    const ccomplex_t IT_0304 = 18*conj(IT_0240);
    const ccomplex_t IT_0305 = (-18)*conj(IT_0230);
    const ccomplex_t IT_0306 = (-6)*conj(IT_0254);
    const ccomplex_t IT_0307 = IT_0303 + IT_0304 + IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = IT_0090*IT_0200;
    const ccomplex_t IT_0309 = IT_0087*IT_0226*IT_0308;
    const ccomplex_t IT_0310 = conj(IT_0185)*IT_0259 + IT_0185*IT_0262 + (
      -0.166666666666667)*conj(IT_0086)*IT_0285 + (-0.166666666666667)*IT_0086
      *IT_0289 + (-0.166666666666667)*conj(IT_0101)*IT_0293 + (
      -0.166666666666667)*IT_0101*IT_0297 + (-0.166666666666667)*IT_0213*IT_0302
       + (-0.166666666666667)*IT_0212*IT_0307 + (-0.166666666666667)*IT_0223
      *IT_0309 + (-0.166666666666667)*IT_0220*conj(IT_0309);
    const ccomplex_t IT_0311 = 6*s_34;
    const ccomplex_t IT_0312 = 18*s_34;
    const ccomplex_t IT_0313 = (-18)*s_34;
    const ccomplex_t IT_0314 = (-6)*s_34;
    const ccomplex_t IT_0315 = m_d*m_t;
    const ccomplex_t IT_0316 = 18*IT_0283;
    const ccomplex_t IT_0317 = (-18)*IT_0271;
    const ccomplex_t IT_0318 = (-6)*IT_0267;
    const ccomplex_t IT_0319 = 18*conj(IT_0283);
    const ccomplex_t IT_0320 = (-18)*conj(IT_0271);
    const ccomplex_t IT_0321 = (-6)*conj(IT_0267);
    const ccomplex_t IT_0322 = 6*IT_0283;
    const ccomplex_t IT_0323 = (-18)*IT_0267;
    const ccomplex_t IT_0324 = (-6)*IT_0271;
    const ccomplex_t IT_0325 = 6*conj(IT_0283);
    const ccomplex_t IT_0326 = (-18)*conj(IT_0267);
    const ccomplex_t IT_0327 = (-6)*conj(IT_0271);
    const ccomplex_t IT_0328 = conj(IT_0254)*IT_0285 + IT_0254*IT_0289 + conj
      (IT_0230)*IT_0293 + IT_0230*IT_0297 + IT_0307*IT_0309 + IT_0302*conj
      (IT_0309) + conj(IT_0258)*(IT_0316 + IT_0317 + IT_0318) + IT_0258*(IT_0319
       + IT_0320 + IT_0321) + conj(IT_0240)*(IT_0322 + IT_0323 + IT_0324) +
       IT_0240*(IT_0325 + IT_0326 + IT_0327);
    const ccomplex_t IT_0329 = m_d*s_24;
    const ccomplex_t IT_0330 = conj(IT_0101)*(18*IT_0230 + (-18)*IT_0240 + 6
      *IT_0254 + (-6)*IT_0258) + IT_0101*(18*conj(IT_0230) + (-18)*conj(IT_0240)
       + 6*conj(IT_0254) + (-6)*conj(IT_0258)) + 6*conj(IT_0086)*IT_0259 + 6
      *IT_0086*IT_0262 + 6*conj(IT_0185)*(IT_0309 + 0.166666666666667*IT_0316 +
       0.166666666666667*IT_0317 + 0.166666666666667*IT_0318) + 6*IT_0185*(conj
      (IT_0309) + 0.166666666666667*IT_0319 + 0.166666666666667*IT_0320 +
       0.166666666666667*IT_0321) + IT_0213*(18*IT_0309 + IT_0322 + IT_0323 +
       IT_0324) + IT_0212*(18*conj(IT_0309) + IT_0325 + IT_0326 + IT_0327);
    const ccomplex_t IT_0331 = IT_0005*(conj(IT_0086)*IT_0101 + IT_0086*conj
      (IT_0101) + conj(IT_0185)*IT_0212 + IT_0185*IT_0213) + (IT_0086*conj
      (IT_0086) + IT_0101*conj(IT_0101) + IT_0185*conj(IT_0185) + IT_0212
      *IT_0213)*IT_0216 + (-6)*IT_0217*IT_0224 + (-6)*IT_0225*IT_0310 + (conj
      (IT_0230)*IT_0254 + IT_0230*conj(IT_0254) + conj(IT_0240)*IT_0258 +
       IT_0240*conj(IT_0258) + conj(IT_0267)*IT_0271 + IT_0267*conj(IT_0271) +
       conj(IT_0283)*IT_0309 + IT_0283*conj(IT_0309))*IT_0311 + (IT_0230*conj
      (IT_0230) + IT_0240*conj(IT_0240) + IT_0254*conj(IT_0254) + IT_0258*conj
      (IT_0258) + IT_0267*conj(IT_0267) + IT_0271*conj(IT_0271) + IT_0283*conj
      (IT_0283) + IT_0309*conj(IT_0309))*IT_0312 + (conj(IT_0230)*IT_0240 +
       IT_0230*conj(IT_0240) + conj(IT_0254)*IT_0258 + IT_0254*conj(IT_0258) +
       conj(IT_0271)*IT_0283 + IT_0271*conj(IT_0283) + conj(IT_0267)*IT_0309 +
       IT_0267*conj(IT_0309))*IT_0313 + (conj(IT_0240)*IT_0254 + IT_0240*conj
      (IT_0254) + conj(IT_0230)*IT_0258 + IT_0230*conj(IT_0258) + conj(IT_0267)
      *IT_0283 + IT_0267*conj(IT_0283) + conj(IT_0271)*IT_0309 + IT_0271*conj
      (IT_0309))*IT_0314 + IT_0315*IT_0328 + IT_0329*IT_0330;
    return create_ccomplex_return(IT_0331);
}

