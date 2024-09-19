#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_st_2_to_d_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_st_2_to_d_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0001*conj(U_st_01);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0001;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = pow(m_sd_L, 2);
    const ccomplex_t IT_0007 = pow(m_st_2, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = pow(m_d, 2);
    const ccomplex_t IT_0016 = pow(m_sG, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0006 + IT_0015 + -IT_0016 +
       reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = -IT_0011 + (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0001*IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = cos(theta_W);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = sin(beta);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = m_t*conj(N_u1)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0030 + (-0.75)*IT_0036);
    const ccomplex_t IT_0038 = 0.666666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0025*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0006 + IT_0015 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0001*IT_0021*IT_0023;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B2)*e_em*conj(U_st_11);
    const ccomplex_t IT_0046 = IT_0027*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = m_t*conj(N_u2)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0049 = IT_0033*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + (-0.75)*IT_0050);
    const ccomplex_t IT_0052 = 0.666666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0006 + IT_0015 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0001*IT_0021*IT_0023;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B3)*e_em*conj(U_st_11);
    const ccomplex_t IT_0060 = IT_0027*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = m_t*conj(N_u3)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0063 = IT_0033*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.75)*IT_0064);
    const ccomplex_t IT_0066 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0006 + IT_0015 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0001*IT_0021*IT_0023;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0074 = IT_0027*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = m_t*conj(N_u4)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0077 = IT_0033*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + (-0.75)*IT_0078);
    const ccomplex_t IT_0080 = 0.666666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0072*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0006 + IT_0015 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = -IT_0042 + -IT_0056 + -IT_0070 + -IT_0084;
    const ccomplex_t IT_0086 = conj(N_B1)*e_em;
    const ccomplex_t IT_0087 = IT_0027*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em;
    const ccomplex_t IT_0090 = IT_0001*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-3)*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = N_B1*e_em*conj(U_st_01);
    const ccomplex_t IT_0095 = IT_0027*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = N_W1*e_em*conj(U_st_01);
    const ccomplex_t IT_0098 = IT_0001*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = m_t*N_u1*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0101 = IT_0033*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + 3*IT_0099 + 3
      *IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = IT_0040*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B2)*e_em;
    const ccomplex_t IT_0109 = IT_0027*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = conj(N_W2)*e_em;
    const ccomplex_t IT_0112 = IT_0001*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-3)*IT_0113);
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = N_B2*e_em*conj(U_st_01);
    const ccomplex_t IT_0117 = IT_0027*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W2*e_em*conj(U_st_01);
    const ccomplex_t IT_0120 = IT_0001*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_t*N_u2*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0123 = IT_0033*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0118 + 3*IT_0121 + 3
      *IT_0124);
    const ccomplex_t IT_0126 = 0.166666666666667*IT_0125;
    const ccomplex_t IT_0127 = IT_0115*IT_0126;
    const ccomplex_t IT_0128 = IT_0054*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em;
    const ccomplex_t IT_0131 = IT_0027*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W3)*e_em;
    const ccomplex_t IT_0134 = IT_0001*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-3)*IT_0135);
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0138 = N_B3*e_em*conj(U_st_01);
    const ccomplex_t IT_0139 = IT_0027*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = N_W3*e_em*conj(U_st_01);
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = m_t*N_u3*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0145 = IT_0033*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0140 + 3*IT_0143 + 3
      *IT_0146);
    const ccomplex_t IT_0148 = 0.166666666666667*IT_0147;
    const ccomplex_t IT_0149 = IT_0137*IT_0148;
    const ccomplex_t IT_0150 = IT_0068*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*e_em;
    const ccomplex_t IT_0153 = IT_0027*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W4)*e_em;
    const ccomplex_t IT_0156 = IT_0001*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + (-3)*IT_0157);
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0161 = IT_0027*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0164 = IT_0001*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = m_t*N_u4*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0167 = IT_0033*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0162 + 3*IT_0165 + 3
      *IT_0168);
    const ccomplex_t IT_0170 = 0.166666666666667*IT_0169;
    const ccomplex_t IT_0171 = IT_0159*IT_0170;
    const ccomplex_t IT_0172 = IT_0082*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = 0.166666666666667*IT_0019 + -IT_0107 + -IT_0129
       + -IT_0151 + -IT_0173;
    const ccomplex_t IT_0175 = s_34*IT_0006;
    const ccomplex_t IT_0176 = (-18)*IT_0175;
    const ccomplex_t IT_0177 = s_13*s_14;
    const ccomplex_t IT_0178 = 36*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = m_d*m_t*IT_0006;
    const ccomplex_t IT_0181 = 3*IT_0174;
    const ccomplex_t IT_0182 = IT_0020 + IT_0181;
    const ccomplex_t IT_0183 = 3*conj(IT_0174);
    const ccomplex_t IT_0184 = conj(IT_0020) + IT_0183;
    const ccomplex_t IT_0185 = m_d*s_14;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0187 = IT_0012*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0189 = IT_0017*IT_0187*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0191 = IT_0038*IT_0093;
    const ccomplex_t IT_0192 = IT_0040*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0194 = IT_0052*IT_0115;
    const ccomplex_t IT_0195 = IT_0054*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0197 = IT_0066*IT_0137;
    const ccomplex_t IT_0198 = IT_0068*IT_0196*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0200 = IT_0080*IT_0159;
    const ccomplex_t IT_0201 = IT_0082*IT_0199*IT_0200;
    const ccomplex_t IT_0202 = m_d*IT_0042;
    const ccomplex_t IT_0203 = m_d*IT_0056;
    const ccomplex_t IT_0204 = m_d*IT_0070;
    const ccomplex_t IT_0205 = m_d*IT_0084;
    const ccomplex_t IT_0206 = 0.166666666666667*IT_0189 + -IT_0192 + -IT_0195
       + -IT_0198 + -IT_0201 + IT_0202 + IT_0203 + IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = pow(m_W, 2);
    const ccomplex_t IT_0208 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + IT_0207 
      + -reg_prop, -1);
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *e_em*IT_0001*IT_0021;
    const ccomplex_t IT_0210 = 0.5*IT_0209;
    const ccomplex_t IT_0211 = pow(m_t, 2);
    const ccomplex_t IT_0212 = cpow(IT_0022, 2);
    const ccomplex_t IT_0213 = cpow(IT_0031, 2);
    const ccomplex_t IT_0214 = A_t*IT_0031*conj(U_st_11) + (-0.707106781186548
      )*e_em*IT_0001*((IT_0021*(IT_0015 + -IT_0211) + -m_W*(IT_0212 + -IT_0213))
      *conj(U_st_01) + m_t*mu_h*IT_0021*IT_0022*IT_0032*conj(U_st_11));
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*V_td*IT_0214;
    const ccomplex_t IT_0216 = IT_0210*IT_0215;
    const ccomplex_t IT_0217 = IT_0208*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = (-0.5)*IT_0002;
    const ccomplex_t IT_0220 = m_t*IT_0219;
    const ccomplex_t IT_0221 = m_t*IT_0003;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = IT_0220 + IT_0222;
    const ccomplex_t IT_0224 = IT_0010*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *e_em*IT_0001*IT_0021*IT_0022*IT_0032;
    const ccomplex_t IT_0226 = 0.5*IT_0225;
    const ccomplex_t IT_0227 = A_t*IT_0022*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0001*(m_W*IT_0022*IT_0031*conj(U_st_01) + (-0.5)*IT_0021*(IT_0015
      *IT_0023*IT_0031*conj(U_st_01) + m_t*(m_t*IT_0022*IT_0032*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*V_td*IT_0227;
    const ccomplex_t IT_0229 = IT_0226*IT_0228;
    const ccomplex_t IT_0230 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0231 = IT_0229*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = (-0.5)*IT_0189 + -IT_0218 + IT_0224 + -IT_0232;
    const ccomplex_t IT_0234 = IT_0025*IT_0104;
    const ccomplex_t IT_0235 = IT_0040*IT_0190*IT_0234;
    const ccomplex_t IT_0236 = IT_0044*IT_0126;
    const ccomplex_t IT_0237 = IT_0054*IT_0193*IT_0236;
    const ccomplex_t IT_0238 = IT_0058*IT_0148;
    const ccomplex_t IT_0239 = IT_0068*IT_0196*IT_0238;
    const ccomplex_t IT_0240 = IT_0072*IT_0170;
    const ccomplex_t IT_0241 = IT_0082*IT_0199*IT_0240;
    const ccomplex_t IT_0242 = m_d*IT_0107;
    const ccomplex_t IT_0243 = m_d*IT_0129;
    const ccomplex_t IT_0244 = m_d*IT_0151;
    const ccomplex_t IT_0245 = m_d*IT_0173;
    const ccomplex_t IT_0246 = m_d*IT_0019;
    const ccomplex_t IT_0247 = -IT_0235 + -IT_0237 + -IT_0239 + -IT_0241 +
       IT_0242 + IT_0243 + IT_0244 + IT_0245 + (-0.166666666666667)*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_td)
      *e_em*IT_0001*IT_0021*IT_0023*IT_0031;
    const ccomplex_t IT_0249 = 0.5*IT_0248;
    const ccomplex_t IT_0250 = IT_0228*IT_0249;
    const ccomplex_t IT_0251 = IT_0230*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = m_d*IT_0219;
    const ccomplex_t IT_0254 = m_d*IT_0003;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = IT_0253 + IT_0255;
    const ccomplex_t IT_0257 = IT_0010*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_td)
      *e_em*IT_0001*IT_0021;
    const ccomplex_t IT_0259 = (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0215*IT_0259;
    const ccomplex_t IT_0261 = IT_0208*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = 0.5*IT_0246 + -IT_0252 + -IT_0257 + -IT_0262;
    const ccomplex_t IT_0264 = IT_0247 + 3*IT_0263;
    const ccomplex_t IT_0265 = 3*conj(IT_0263);
    const ccomplex_t IT_0266 = conj(IT_0247) + IT_0265;
    const ccomplex_t IT_0267 = conj(IT_0085)*(18*IT_0206 + 6*IT_0233) +
       IT_0085*(18*conj(IT_0206) + 6*conj(IT_0233)) + conj(IT_0174)*(18*IT_0247 
      + 6*IT_0263) + IT_0174*(18*conj(IT_0247) + 6*conj(IT_0263)) + 6*conj
      (IT_0020)*IT_0264 + 6*IT_0020*IT_0266;
    const ccomplex_t IT_0268 = m_d*m_t;
    const ccomplex_t IT_0269 = (-18)*IT_0247;
    const ccomplex_t IT_0270 = (-6)*IT_0263;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = (-18)*conj(IT_0247);
    const ccomplex_t IT_0273 = (-6)*conj(IT_0263);
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = m_t*s_13;
    const ccomplex_t IT_0276 = (-18)*(0.333333333333333*IT_0184*IT_0206 +
       0.333333333333333*IT_0182*conj(IT_0206) + (conj(IT_0020) +
       0.333333333333333*conj(IT_0174))*IT_0233 + (IT_0020 + 0.333333333333333
      *IT_0174)*conj(IT_0233) + (-0.0555555555555556)*conj(IT_0085)*IT_0271 + (
      -0.0555555555555556)*IT_0085*IT_0274)*IT_0275;
    const ccomplex_t IT_0277 = 6*s_34;
    const ccomplex_t IT_0278 = 18*s_34;
    const ccomplex_t IT_0279 = 12*IT_0177;
    const ccomplex_t IT_0280 = (-6)*IT_0175;
    const ccomplex_t IT_0281 = (IT_0020*conj(IT_0020) + IT_0085*conj(IT_0085) 
      + IT_0174*conj(IT_0174))*IT_0179 + (-6)*IT_0180*(conj(IT_0085)*IT_0182 +
       IT_0085*IT_0184) + IT_0185*IT_0267 + (-6)*IT_0268*(conj(IT_0233)*IT_0264 
      + IT_0233*IT_0266 + (-0.166666666666667)*conj(IT_0206)*IT_0271 + (
      -0.166666666666667)*IT_0206*IT_0274) + IT_0276 + (conj(IT_0206)*IT_0233 +
       IT_0206*conj(IT_0233) + conj(IT_0247)*IT_0263 + IT_0247*conj(IT_0263))
      *IT_0277 + (IT_0206*conj(IT_0206) + IT_0233*conj(IT_0233) + IT_0247*conj
      (IT_0247) + IT_0263*conj(IT_0263))*IT_0278 + (conj(IT_0020)*IT_0174 +
       IT_0020*conj(IT_0174))*(IT_0279 + IT_0280);
    return create_ccomplex_return(IT_0281);
}

