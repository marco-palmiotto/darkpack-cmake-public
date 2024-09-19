#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_st_1_to_s_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_st_1_to_s_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0001*IT_0003*IT_0005*U_st_00;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_s, 2);
    const ccomplex_t IT_0010 = pow(m_st_1, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0001*IT_0003*IT_0005*U_st_00;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = pow(m_sG, 2);
    const ccomplex_t IT_0026 = pow(m_ss_L, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0009 + -IT_0025 + IT_0026 +
       reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0009 + IT_0026 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0031 = cos(theta_W);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = conj(N_B1)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W1)*e_em;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-3)*IT_0038);
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = sin(beta);
    const ccomplex_t IT_0045 = cpow(IT_0044, -1);
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = m_t*N_u1*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0043 + (-0.75)*IT_0049);
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = IT_0030*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0009 + IT_0026 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0056 = conj(N_B2)*e_em;
    const ccomplex_t IT_0057 = IT_0032*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W2)*e_em;
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-3)*IT_0061);
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0065 = IT_0032*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = m_t*N_u2*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0068 = IT_0046*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.75)*IT_0069);
    const ccomplex_t IT_0071 = 0.666666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = IT_0055*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_13 + IT_0009 + IT_0026 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0076 = m_t*N_u3*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0077 = IT_0046*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0080 = IT_0032*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + (
      -1.33333333333333)*IT_0081);
    const ccomplex_t IT_0083 = (-0.75)*IT_0082;
    const ccomplex_t IT_0084 = (-1.33333333333333)*IT_0083;
    const ccomplex_t IT_0085 = (-0.75)*IT_0084;
    const ccomplex_t IT_0086 = 0.666666666666667*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B3)*e_em;
    const ccomplex_t IT_0088 = IT_0032*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W3)*e_em;
    const ccomplex_t IT_0091 = IT_0001*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + (-3)*IT_0092);
    const ccomplex_t IT_0094 = 0.166666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0086*IT_0094;
    const ccomplex_t IT_0096 = IT_0075*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_13 + IT_0009 + IT_0026 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0099 = m_t*N_u4*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0100 = IT_0046*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0103 = IT_0032*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + (
      -1.33333333333333)*IT_0104);
    const ccomplex_t IT_0106 = (-0.75)*IT_0105;
    const ccomplex_t IT_0107 = 0.666666666666667*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B4)*e_em;
    const ccomplex_t IT_0109 = IT_0032*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = conj(N_W4)*e_em;
    const ccomplex_t IT_0112 = IT_0001*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-3)*IT_0113);
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = IT_0107*IT_0115;
    const ccomplex_t IT_0117 = IT_0098*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0029 + -IT_0054 + -IT_0074
       + -IT_0097 + -IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0029;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0001*IT_0003*IT_0045;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0124 = IT_0001*IT_0123;
    const ccomplex_t IT_0125 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0126 = IT_0046*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + (-0.5)*IT_0127);
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = IT_0122*IT_0129;
    const ccomplex_t IT_0131 = IT_0011*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0001*IT_0003*IT_0045;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0136 = IT_0001*IT_0135;
    const ccomplex_t IT_0137 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0138 = IT_0046*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + (-0.5)*IT_0139);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0134*IT_0141;
    const ccomplex_t IT_0143 = IT_0018*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = -IT_0132 + -IT_0144;
    const ccomplex_t IT_0146 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0147 = IT_0032*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0150 = IT_0001*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = m_t*conj(N_u1)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0153 = IT_0046*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0148 + 3*IT_0151 + 3
      *IT_0154);
    const ccomplex_t IT_0156 = 0.166666666666667*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = IT_0156*IT_0158;
    const ccomplex_t IT_0160 = IT_0030*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0163 = IT_0032*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0166 = IT_0001*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = m_t*conj(N_u2)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0169 = IT_0046*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0164 + 3*IT_0167 + 3
      *IT_0170);
    const ccomplex_t IT_0172 = 0.166666666666667*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0174 = (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = IT_0172*IT_0174;
    const ccomplex_t IT_0176 = IT_0055*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0179 = IT_0032*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0182 = IT_0001*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = m_t*conj(N_u3)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0185 = IT_0046*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0180 + 3*IT_0183 + 3
      *IT_0186);
    const ccomplex_t IT_0188 = 0.166666666666667*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0188*IT_0190;
    const ccomplex_t IT_0192 = IT_0075*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0197 = IT_0032*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0200 = IT_0001*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = m_t*conj(N_u4)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0203 = IT_0046*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0198 + 3*IT_0201 + 3
      *IT_0204);
    const ccomplex_t IT_0206 = 0.166666666666667*IT_0205;
    const ccomplex_t IT_0207 = IT_0195*IT_0206;
    const ccomplex_t IT_0208 = IT_0098*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = -IT_0161 + -IT_0177 + -IT_0193 + -IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0021)*(IT_0119 + 3*IT_0120 + 3*IT_0145)
       + IT_0021*(conj(IT_0119) + 3*conj(IT_0120) + 3*conj(IT_0145)) + (
      -0.166666666666667)*((-18)*conj(IT_0119) + (-6)*conj(IT_0120) + (-6)*conj
      (IT_0145))*IT_0210 + (-0.166666666666667)*((-18)*IT_0119 + (-6)*IT_0120 + 
      (-6)*IT_0145)*conj(IT_0210);
    const ccomplex_t IT_0212 = m_s*m_t*IT_0010;
    const ccomplex_t IT_0213 = s_34*IT_0010;
    const ccomplex_t IT_0214 = (-6)*IT_0213;
    const ccomplex_t IT_0215 = s_23*s_24;
    const ccomplex_t IT_0216 = 12*IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = IT_0120 + IT_0145;
    const ccomplex_t IT_0219 = conj(IT_0120) + conj(IT_0145);
    const ccomplex_t IT_0220 = (-18)*IT_0213;
    const ccomplex_t IT_0221 = 36*IT_0215;
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = m_s*m_t;
    const ccomplex_t IT_0224 = m_t*IT_0029;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0227 = IT_0086*IT_0190;
    const ccomplex_t IT_0228 = IT_0075*IT_0226*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0230 = IT_0107*IT_0195;
    const ccomplex_t IT_0231 = IT_0098*IT_0229*IT_0230;
    const ccomplex_t IT_0232 = -IT_0228 + -IT_0231;
    const ccomplex_t IT_0233 = 6*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0235 = IT_0071*IT_0174;
    const ccomplex_t IT_0236 = IT_0055*IT_0234*IT_0235;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = 6*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0240 = IT_0051*IT_0158;
    const ccomplex_t IT_0241 = IT_0030*IT_0239*IT_0240;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = 6*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0245 = IT_0016*IT_0134;
    const ccomplex_t IT_0246 = IT_0018*IT_0244*IT_0245;
    const ccomplex_t IT_0247 = m_s*IT_0132;
    const ccomplex_t IT_0248 = m_s*IT_0144;
    const ccomplex_t IT_0249 = -IT_0246 + -IT_0247 + -IT_0248;
    const ccomplex_t IT_0250 = 18*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0252 = IT_0007*IT_0122;
    const ccomplex_t IT_0253 = IT_0011*IT_0251*IT_0252;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = 18*IT_0254;
    const ccomplex_t IT_0256 = IT_0161 + IT_0177 + IT_0193 + IT_0209;
    const ccomplex_t IT_0257 = m_t*IT_0256;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = (-6)*IT_0258;
    const ccomplex_t IT_0260 = IT_0233 + IT_0238 + IT_0243 + IT_0250 + IT_0255
       + IT_0259;
    const ccomplex_t IT_0261 = 6*conj(IT_0232);
    const ccomplex_t IT_0262 = 6*conj(IT_0237);
    const ccomplex_t IT_0263 = 6*conj(IT_0242);
    const ccomplex_t IT_0264 = 18*conj(IT_0249);
    const ccomplex_t IT_0265 = 18*conj(IT_0254);
    const ccomplex_t IT_0266 = (-6)*conj(IT_0258);
    const ccomplex_t IT_0267 = IT_0261 + IT_0262 + IT_0263 + IT_0264 + IT_0265
       + IT_0266;
    const ccomplex_t IT_0268 = (-6)*IT_0054 + (-6)*IT_0074 + (-6)*IT_0097 + (
      -6)*IT_0118;
    const ccomplex_t IT_0269 = IT_0029 + IT_0268;
    const ccomplex_t IT_0270 = m_t*IT_0269;
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0270;
    const ccomplex_t IT_0272 = 6*IT_0249;
    const ccomplex_t IT_0273 = 6*IT_0254;
    const ccomplex_t IT_0274 = 18*IT_0232;
    const ccomplex_t IT_0275 = 18*IT_0237;
    const ccomplex_t IT_0276 = 18*IT_0242;
    const ccomplex_t IT_0277 = (-18)*IT_0258;
    const ccomplex_t IT_0278 = IT_0272 + IT_0273 + IT_0274 + IT_0275 + IT_0276
       + IT_0277;
    const ccomplex_t IT_0279 = 6*conj(IT_0249);
    const ccomplex_t IT_0280 = 6*conj(IT_0254);
    const ccomplex_t IT_0281 = 18*conj(IT_0232);
    const ccomplex_t IT_0282 = 18*conj(IT_0237);
    const ccomplex_t IT_0283 = 18*conj(IT_0242);
    const ccomplex_t IT_0284 = (-18)*conj(IT_0258);
    const ccomplex_t IT_0285 = IT_0279 + IT_0280 + IT_0281 + IT_0282 + IT_0283
       + IT_0284;
    const ccomplex_t IT_0286 = IT_0002*IT_0129;
    const ccomplex_t IT_0287 = IT_0011*IT_0251*IT_0286;
    const ccomplex_t IT_0288 = IT_0014*IT_0141;
    const ccomplex_t IT_0289 = IT_0018*IT_0244*IT_0288;
    const ccomplex_t IT_0290 = m_s*IT_0013;
    const ccomplex_t IT_0291 = m_s*IT_0020;
    const ccomplex_t IT_0292 = IT_0287 + IT_0289 + IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0294 = IT_0022*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0296 = IT_0027*IT_0294*IT_0295;
    const ccomplex_t IT_0297 = (-0.5)*IT_0296;
    const ccomplex_t IT_0298 = conj(IT_0292) + conj(IT_0297);
    const ccomplex_t IT_0299 = 6*IT_0258;
    const ccomplex_t IT_0300 = (-18)*IT_0249;
    const ccomplex_t IT_0301 = (-18)*IT_0254;
    const ccomplex_t IT_0302 = (-6)*IT_0232;
    const ccomplex_t IT_0303 = (-6)*IT_0237;
    const ccomplex_t IT_0304 = (-6)*IT_0242;
    const ccomplex_t IT_0305 = IT_0299 + IT_0300 + IT_0301 + IT_0302 + IT_0303
       + IT_0304;
    const ccomplex_t IT_0306 = IT_0292 + IT_0297;
    const ccomplex_t IT_0307 = 6*conj(IT_0258);
    const ccomplex_t IT_0308 = (-18)*conj(IT_0249);
    const ccomplex_t IT_0309 = (-18)*conj(IT_0254);
    const ccomplex_t IT_0310 = (-6)*conj(IT_0232);
    const ccomplex_t IT_0311 = (-6)*conj(IT_0237);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0242);
    const ccomplex_t IT_0313 = IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311
       + IT_0312;
    const ccomplex_t IT_0314 = IT_0063*IT_0172;
    const ccomplex_t IT_0315 = IT_0055*IT_0234*IT_0314;
    const ccomplex_t IT_0316 = IT_0115*IT_0206;
    const ccomplex_t IT_0317 = IT_0098*IT_0229*IT_0316;
    const ccomplex_t IT_0318 = -IT_0315 + -IT_0317;
    const ccomplex_t IT_0319 = IT_0094*IT_0188;
    const ccomplex_t IT_0320 = IT_0075*IT_0226*IT_0319;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0040*IT_0156;
    const ccomplex_t IT_0323 = IT_0030*IT_0239*IT_0322;
    const ccomplex_t IT_0324 = 0.166666666666667*IT_0296 + -IT_0323;
    const ccomplex_t IT_0325 = conj(IT_0318) + conj(IT_0321) + conj(IT_0324);
    const ccomplex_t IT_0326 = 18*IT_0258;
    const ccomplex_t IT_0327 = (-18)*IT_0232;
    const ccomplex_t IT_0328 = (-18)*IT_0237;
    const ccomplex_t IT_0329 = (-18)*IT_0242;
    const ccomplex_t IT_0330 = (-6)*IT_0249;
    const ccomplex_t IT_0331 = (-6)*IT_0254;
    const ccomplex_t IT_0332 = IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330
       + IT_0331;
    const ccomplex_t IT_0333 = IT_0318 + IT_0321 + IT_0324;
    const ccomplex_t IT_0334 = 18*conj(IT_0258);
    const ccomplex_t IT_0335 = (-18)*conj(IT_0232);
    const ccomplex_t IT_0336 = (-18)*conj(IT_0237);
    const ccomplex_t IT_0337 = (-18)*conj(IT_0242);
    const ccomplex_t IT_0338 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0339 = (-6)*conj(IT_0254);
    const ccomplex_t IT_0340 = IT_0334 + IT_0335 + IT_0336 + IT_0337 + IT_0338
       + IT_0339;
    const ccomplex_t IT_0341 = IT_0232 + IT_0242;
    const ccomplex_t IT_0342 = IT_0321 + IT_0324;
    const ccomplex_t IT_0343 = conj(IT_0321) + conj(IT_0324);
    const ccomplex_t IT_0344 = 6*s_34;
    const ccomplex_t IT_0345 = conj(IT_0232) + conj(IT_0242);
    const ccomplex_t IT_0346 = conj(IT_0249) + conj(IT_0254);
    const ccomplex_t IT_0347 = 18*s_34;
    const ccomplex_t IT_0348 = IT_0237*IT_0344;
    const ccomplex_t IT_0349 = IT_0254*IT_0347;
    const ccomplex_t IT_0350 = (-18)*s_34;
    const ccomplex_t IT_0351 = IT_0237*IT_0347;
    const ccomplex_t IT_0352 = IT_0258*IT_0350;
    const ccomplex_t IT_0353 = IT_0254*IT_0344;
    const ccomplex_t IT_0354 = conj(IT_0254)*IT_0344;
    const ccomplex_t IT_0355 = conj(IT_0237)*IT_0347;
    const ccomplex_t IT_0356 = conj(IT_0258)*IT_0350;
    const ccomplex_t IT_0357 = IT_0318*IT_0347;
    const ccomplex_t IT_0358 = IT_0271*IT_0350;
    const ccomplex_t IT_0359 = IT_0292*IT_0344;
    const ccomplex_t IT_0360 = IT_0321*IT_0347;
    const ccomplex_t IT_0361 = conj(IT_0292)*IT_0344;
    const ccomplex_t IT_0362 = conj(IT_0318)*IT_0347;
    const ccomplex_t IT_0363 = conj(IT_0271)*IT_0350;
    const ccomplex_t IT_0364 = (-6)*s_34;
    const ccomplex_t IT_0365 = conj(IT_0258)*IT_0364;
    const ccomplex_t IT_0366 = IT_0318*IT_0344;
    const ccomplex_t IT_0367 = IT_0292*IT_0347;
    const ccomplex_t IT_0368 = IT_0271*IT_0364;
    const ccomplex_t IT_0369 = conj(IT_0318)*IT_0344;
    const ccomplex_t IT_0370 = conj(IT_0271)*IT_0364;
    const ccomplex_t IT_0371 = m_s*s_24;
    const ccomplex_t IT_0372 = (-3)*conj(IT_0271);
    const ccomplex_t IT_0373 = -conj(IT_0225);
    const ccomplex_t IT_0374 = (-0.166666666666667)*IT_0219*IT_0305 + (
      -0.166666666666667)*IT_0218*IT_0313 + (-0.166666666666667)*conj(IT_0021)*
      (18*IT_0225 + 6*IT_0271 + (-18)*IT_0292 + (-18)*IT_0297 + (-6)*IT_0318 + (
      -6)*IT_0321 + (-6)*IT_0324) + -conj(IT_0210)*(IT_0225 + 3*IT_0271 + 
      -IT_0292 + -IT_0297 + (-3)*IT_0318 + (-3)*IT_0321 + (-3)*IT_0324) + (
      -0.166666666666667)*IT_0021*(18*conj(IT_0225) + 6*conj(IT_0271) + (-18)
      *conj(IT_0292) + (-18)*conj(IT_0297) + (-6)*conj(IT_0318) + (-6)*conj
      (IT_0321) + (-6)*conj(IT_0324)) + (-0.166666666666667)*conj(IT_0119)
      *IT_0332 + (-0.166666666666667)*IT_0119*IT_0340 + IT_0210*(conj(IT_0292) +
       conj(IT_0297) + 3*conj(IT_0318) + 3*conj(IT_0321) + 3*conj(IT_0324) +
       IT_0372 + IT_0373);
    const ccomplex_t IT_0375 = m_t*s_23;
    const ccomplex_t IT_0376 = 6*IT_0318;
    const ccomplex_t IT_0377 = 6*IT_0321;
    const ccomplex_t IT_0378 = 6*IT_0324;
    const ccomplex_t IT_0379 = 18*IT_0292;
    const ccomplex_t IT_0380 = 18*IT_0297;
    const ccomplex_t IT_0381 = (-18)*IT_0225;
    const ccomplex_t IT_0382 = (-6)*IT_0271;
    const ccomplex_t IT_0383 = 6*conj(IT_0318);
    const ccomplex_t IT_0384 = 6*conj(IT_0321);
    const ccomplex_t IT_0385 = 6*conj(IT_0324);
    const ccomplex_t IT_0386 = 18*conj(IT_0292);
    const ccomplex_t IT_0387 = 18*conj(IT_0297);
    const ccomplex_t IT_0388 = (-18)*conj(IT_0225);
    const ccomplex_t IT_0389 = (-6)*conj(IT_0271);
    const ccomplex_t IT_0390 = 0.166666666666667*conj(IT_0021)*IT_0260 +
       0.166666666666667*IT_0021*IT_0267 + 0.166666666666667*conj(IT_0210)
      *IT_0278 + 0.166666666666667*IT_0210*IT_0285 + -conj(IT_0119)*(IT_0225 + 3
      *IT_0271 + -IT_0292 + -IT_0297 + (-3)*IT_0318 + (-3)*IT_0321 + (-3)
      *IT_0324) + IT_0119*(conj(IT_0292) + conj(IT_0297) + 3*conj(IT_0318) + 3
      *conj(IT_0321) + 3*conj(IT_0324) + IT_0372 + IT_0373) + 0.166666666666667
      *IT_0219*(IT_0376 + IT_0377 + IT_0378 + IT_0379 + IT_0380 + IT_0381 +
       IT_0382) + 0.166666666666667*IT_0218*(IT_0383 + IT_0384 + IT_0385 +
       IT_0386 + IT_0387 + IT_0388 + IT_0389);
    const ccomplex_t IT_0391 = (-6)*IT_0211*IT_0212 + IT_0217*(conj(IT_0021)
      *IT_0210 + IT_0021*conj(IT_0210) + conj(IT_0119)*IT_0218 + IT_0119*IT_0219
      ) + (IT_0021*conj(IT_0021) + IT_0119*conj(IT_0119) + IT_0210*conj(IT_0210)
       + IT_0218*IT_0219)*IT_0222 + IT_0223*(conj(IT_0225)*IT_0260 + IT_0225
      *IT_0267 + conj(IT_0271)*IT_0278 + IT_0271*IT_0285 + IT_0298*IT_0305 +
       IT_0306*IT_0313 + IT_0325*IT_0332 + IT_0333*IT_0340) + ((conj(IT_0232) +
       conj(IT_0237) + conj(IT_0242))*IT_0249 + conj(IT_0237)*IT_0254 + conj
      (IT_0225)*IT_0271 + IT_0225*conj(IT_0271) + conj(IT_0249)*IT_0341 + conj
      (IT_0297)*IT_0342 + IT_0297*IT_0343)*IT_0344 + (IT_0225*conj(IT_0225) +
       IT_0258*conj(IT_0258) + IT_0271*conj(IT_0271) + IT_0297*IT_0298 + IT_0324
      *IT_0343 + IT_0341*IT_0345 + IT_0249*IT_0346)*IT_0347 + IT_0346*(IT_0348 +
       IT_0349) + (IT_0225*IT_0298 + conj(IT_0225)*IT_0306)*IT_0350 + conj
      (IT_0237)*(IT_0351 + IT_0352) + IT_0345*(IT_0351 + IT_0352 + IT_0353) +
       IT_0341*(IT_0354 + IT_0355 + IT_0356) + IT_0237*IT_0356 + conj(IT_0318)*
      (IT_0357 + IT_0358) + IT_0343*(IT_0357 + IT_0358 + IT_0359 + IT_0360) +
       IT_0342*(IT_0361 + IT_0362 + IT_0363) + IT_0318*IT_0363 + (IT_0225
      *IT_0325 + conj(IT_0225)*IT_0333 + IT_0258*IT_0346)*IT_0364 + (IT_0249 +
       IT_0254)*IT_0365 + IT_0298*(IT_0366 + IT_0367 + IT_0368) + IT_0306*
      (IT_0369 + IT_0370) + (-6)*IT_0371*IT_0374 + 6*IT_0375*IT_0390;
    return create_ccomplex_return(IT_0391);
}

