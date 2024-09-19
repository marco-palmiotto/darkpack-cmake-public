#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_st_1_to_s_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_st_1_to_s_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0001*conj(U_st_00);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0001;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = pow(m_ss_L, 2);
    const ccomplex_t IT_0007 = pow(m_st_1, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = pow(m_s, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0006 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_B1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = conj(N_W1)*e_em;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + (-3)*IT_0022);
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = sin(beta);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = m_t*N_u1*e_em*IT_0034*conj(U_st_10);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0027 + 3*IT_0030 + 3
      *IT_0037);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = IT_0024*IT_0039;
    const ccomplex_t IT_0041 = IT_0014*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0006 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0044 = conj(N_B4)*e_em;
    const ccomplex_t IT_0045 = IT_0016*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W4)*e_em;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0053 = IT_0016*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0056 = IT_0001*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = m_t*N_u4*e_em*IT_0034*conj(U_st_10);
    const ccomplex_t IT_0059 = IT_0033*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0054 + 3*IT_0057 + 3
      *IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = IT_0051*IT_0062;
    const ccomplex_t IT_0064 = IT_0043*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0006 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0067 = conj(N_B2)*e_em;
    const ccomplex_t IT_0068 = IT_0016*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W2)*e_em;
    const ccomplex_t IT_0071 = IT_0001*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-3)*IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0076 = IT_0016*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0079 = IT_0001*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = m_t*N_u2*e_em*IT_0034*conj(U_st_10);
    const ccomplex_t IT_0082 = IT_0033*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0077 + 3*IT_0080 + 3
      *IT_0083);
    const ccomplex_t IT_0085 = 0.166666666666667*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = IT_0066*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0006 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0090 = conj(N_B3)*e_em;
    const ccomplex_t IT_0091 = IT_0016*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W3)*e_em;
    const ccomplex_t IT_0094 = IT_0001*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + (-3)*IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0099 = IT_0016*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_t*N_u3*e_em*IT_0034*conj(U_st_10);
    const ccomplex_t IT_0105 = IT_0033*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103 + 3
      *IT_0106);
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = IT_0097*IT_0108;
    const ccomplex_t IT_0110 = IT_0089*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0042 + -IT_0065 + -IT_0088 + -IT_0111;
    const ccomplex_t IT_0113 = pow(m_sG, 2);
    const ccomplex_t IT_0114 = cpow((-2)*s_13 + IT_0006 + IT_0013 + -IT_0113 +
       reg_prop, -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = IT_0114*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0119;
    const ccomplex_t IT_0122 = cos(beta);
    const ccomplex_t IT_0123 = cpow(IT_0122, -1);
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0001*IT_0034*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0127 = IT_0016*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_t*conj(N_u2)*e_em*IT_0034*conj(U_st_00);
    const ccomplex_t IT_0130 = IT_0033*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + (-0.75)*IT_0131);
    const ccomplex_t IT_0133 = 0.666666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0125*IT_0133;
    const ccomplex_t IT_0135 = IT_0066*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0001*IT_0034*IT_0123;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0140 = IT_0016*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = m_t*conj(N_u3)*e_em*IT_0034*conj(U_st_00);
    const ccomplex_t IT_0143 = IT_0033*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + (-0.75)*IT_0144);
    const ccomplex_t IT_0146 = 0.666666666666667*IT_0145;
    const ccomplex_t IT_0147 = IT_0138*IT_0146;
    const ccomplex_t IT_0148 = IT_0089*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0001*IT_0034*IT_0123;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0153 = IT_0016*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_t*conj(N_u4)*e_em*IT_0034*conj(U_st_00);
    const ccomplex_t IT_0156 = IT_0033*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + (-0.75)*IT_0157);
    const ccomplex_t IT_0159 = 0.666666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0151*IT_0159;
    const ccomplex_t IT_0161 = IT_0043*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0001*IT_0034*IT_0123;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0166 = IT_0016*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = m_t*conj(N_u1)*e_em*IT_0034*conj(U_st_00);
    const ccomplex_t IT_0169 = IT_0033*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0167 + (-0.75)*IT_0170);
    const ccomplex_t IT_0172 = 0.666666666666667*IT_0171;
    const ccomplex_t IT_0173 = IT_0164*IT_0172;
    const ccomplex_t IT_0174 = IT_0014*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = -IT_0136 + -IT_0149 + -IT_0162 + -IT_0175;
    const ccomplex_t IT_0177 = m_s*m_t*IT_0006;
    const ccomplex_t IT_0178 = s_34*IT_0006;
    const ccomplex_t IT_0179 = (-18)*IT_0178;
    const ccomplex_t IT_0180 = s_13*s_14;
    const ccomplex_t IT_0181 = 36*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = IT_0012 + IT_0120;
    const ccomplex_t IT_0184 = conj(IT_0012) + conj(IT_0120);
    const ccomplex_t IT_0185 = conj(IT_0112) + conj(IT_0121);
    const ccomplex_t IT_0186 = (-6)*IT_0178;
    const ccomplex_t IT_0187 = 12*IT_0180;
    const ccomplex_t IT_0188 = IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = IT_0112 + IT_0121;
    const ccomplex_t IT_0190 = IT_0085*IT_0125;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0192 = IT_0066*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = IT_0108*IT_0138;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0195 = IT_0089*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = IT_0062*IT_0151;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0198 = IT_0043*IT_0196*IT_0197;
    const ccomplex_t IT_0199 = m_s*IT_0119;
    const ccomplex_t IT_0200 = m_s*IT_0065;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0202 = IT_0039*IT_0164;
    const ccomplex_t IT_0203 = IT_0014*IT_0201*IT_0202;
    const ccomplex_t IT_0204 = m_s*IT_0042;
    const ccomplex_t IT_0205 = m_s*IT_0088;
    const ccomplex_t IT_0206 = m_s*IT_0111;
    const ccomplex_t IT_0207 = -IT_0192 + -IT_0195 + -IT_0198 + (
      -0.166666666666667)*IT_0199 + IT_0200 + -IT_0203 + IT_0204 + IT_0205 +
       IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0001*IT_0034;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = pow(m_t, 2);
    const ccomplex_t IT_0211 = cpow(IT_0122, 2);
    const ccomplex_t IT_0212 = cpow(IT_0031, 2);
    const ccomplex_t IT_0213 = A_t*IT_0031*conj(U_st_10) + (-0.707106781186548
      )*e_em*IT_0001*((IT_0034*(IT_0013 + -IT_0210) + -m_W*(IT_0211 + -IT_0212))
      *conj(U_st_00) + m_t*mu_h*IT_0032*IT_0034*IT_0122*conj(U_st_10));
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*V_ts*IT_0213;
    const ccomplex_t IT_0215 = IT_0209*IT_0214;
    const ccomplex_t IT_0216 = pow(m_W, 2);
    const ccomplex_t IT_0217 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + IT_0216 
      + -reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0215*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = A_t*IT_0122*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0001*(m_W*IT_0031*IT_0122*conj(U_st_00) + (-0.5)*IT_0034*(IT_0013
      *IT_0031*IT_0123*conj(U_st_00) + m_t*(m_t*IT_0032*IT_0122*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*V_ts*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0001*IT_0031*IT_0034*IT_0123;
    const ccomplex_t IT_0223 = 0.5*IT_0222;
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0226 = IT_0224*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = (-0.5)*IT_0002;
    const ccomplex_t IT_0229 = m_s*IT_0228;
    const ccomplex_t IT_0230 = m_s*IT_0003;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = IT_0229 + IT_0231;
    const ccomplex_t IT_0233 = IT_0010*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0199 + -IT_0219 + -IT_0227 + -IT_0233;
    const ccomplex_t IT_0235 = IT_0024*IT_0172;
    const ccomplex_t IT_0236 = IT_0014*IT_0201*IT_0235;
    const ccomplex_t IT_0237 = m_s*IT_0175;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0240 = IT_0115*IT_0239;
    const ccomplex_t IT_0241 = IT_0114*IT_0238*IT_0240;
    const ccomplex_t IT_0242 = IT_0074*IT_0133;
    const ccomplex_t IT_0243 = IT_0066*IT_0191*IT_0242;
    const ccomplex_t IT_0244 = IT_0097*IT_0146;
    const ccomplex_t IT_0245 = IT_0089*IT_0194*IT_0244;
    const ccomplex_t IT_0246 = IT_0051*IT_0159;
    const ccomplex_t IT_0247 = IT_0043*IT_0197*IT_0246;
    const ccomplex_t IT_0248 = m_s*IT_0136;
    const ccomplex_t IT_0249 = m_s*IT_0149;
    const ccomplex_t IT_0250 = m_s*IT_0162;
    const ccomplex_t IT_0251 = -IT_0236 + IT_0237 + 0.166666666666667*IT_0241 
      + -IT_0243 + -IT_0245 + -IT_0247 + IT_0248 + IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = (-0.5)*IT_0241;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0001*IT_0034;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = IT_0214*IT_0254;
    const ccomplex_t IT_0256 = IT_0217*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = m_t*IT_0228;
    const ccomplex_t IT_0259 = m_t*IT_0003;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = IT_0010*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0001*IT_0032*IT_0034*IT_0122;
    const ccomplex_t IT_0264 = 0.5*IT_0263;
    const ccomplex_t IT_0265 = IT_0221*IT_0264;
    const ccomplex_t IT_0266 = IT_0225*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = -IT_0257 + IT_0262 + -IT_0267;
    const ccomplex_t IT_0269 = IT_0252 + IT_0268;
    const ccomplex_t IT_0270 = conj(IT_0252) + conj(IT_0268);
    const ccomplex_t IT_0271 = 6*s_34;
    const ccomplex_t IT_0272 = 18*s_34;
    const ccomplex_t IT_0273 = IT_0252*IT_0272;
    const ccomplex_t IT_0274 = m_s*s_14;
    const ccomplex_t IT_0275 = 6*IT_0234;
    const ccomplex_t IT_0276 = 18*IT_0207;
    const ccomplex_t IT_0277 = 6*conj(IT_0234);
    const ccomplex_t IT_0278 = 18*conj(IT_0207);
    const ccomplex_t IT_0279 = 6*IT_0207;
    const ccomplex_t IT_0280 = 18*IT_0234;
    const ccomplex_t IT_0281 = 6*conj(IT_0207);
    const ccomplex_t IT_0282 = 18*conj(IT_0234);
    const ccomplex_t IT_0283 = 3*conj(IT_0176)*(IT_0251 + 1./3*IT_0252 + 1./3
      *IT_0268) + 3*IT_0176*(conj(IT_0251) + 1./3*conj(IT_0252) + 1./3*conj
      (IT_0268)) + 0.166666666666667*IT_0185*(IT_0275 + IT_0276) +
       0.166666666666667*IT_0189*(IT_0277 + IT_0278) + 0.166666666666667*IT_0184
      *(IT_0279 + IT_0280) + 0.166666666666667*IT_0183*(IT_0281 + IT_0282);
    const ccomplex_t IT_0284 = m_s*m_t;
    const ccomplex_t IT_0285 = (-18)*conj(IT_0207);
    const ccomplex_t IT_0286 = (-6)*conj(IT_0234);
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = (-18)*IT_0234;
    const ccomplex_t IT_0289 = (-18)*conj(IT_0234);
    const ccomplex_t IT_0290 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0291 = (-6)*IT_0234;
    const ccomplex_t IT_0292 = m_t*s_13;
    const ccomplex_t IT_0293 = (-18)*IT_0252;
    const ccomplex_t IT_0294 = (-18)*IT_0268;
    const ccomplex_t IT_0295 = (-6)*IT_0251;
    const ccomplex_t IT_0296 = (-18)*conj(IT_0252);
    const ccomplex_t IT_0297 = (-18)*conj(IT_0268);
    const ccomplex_t IT_0298 = (-6)*conj(IT_0251);
    const ccomplex_t IT_0299 = (-18)*IT_0251;
    const ccomplex_t IT_0300 = (-6)*IT_0252;
    const ccomplex_t IT_0301 = (-6)*IT_0268;
    const ccomplex_t IT_0302 = (-18)*conj(IT_0251);
    const ccomplex_t IT_0303 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0304 = (-6)*conj(IT_0268);
    const ccomplex_t IT_0305 = IT_0176*IT_0287 + conj(IT_0176)*((-18)*IT_0207 
      + IT_0291) + IT_0184*(IT_0293 + IT_0294 + IT_0295) + IT_0183*(IT_0296 +
       IT_0297 + IT_0298) + IT_0185*(IT_0299 + IT_0300 + IT_0301) + IT_0189*
      (IT_0302 + IT_0303 + IT_0304);
    const ccomplex_t IT_0306 = (-6)*((conj(IT_0012) + 3*conj(IT_0112) + conj
      (IT_0120) + 3*conj(IT_0121))*IT_0176 + (IT_0012 + 3*IT_0112 + IT_0120 + 3
      *IT_0121)*conj(IT_0176))*IT_0177 + IT_0182*(IT_0176*conj(IT_0176) +
       IT_0183*IT_0184 + (IT_0112 + IT_0121)*IT_0185) + IT_0188*(IT_0183*IT_0185
       + IT_0184*IT_0189) + (conj(IT_0207)*IT_0234 + IT_0207*conj(IT_0234) +
       conj(IT_0251)*IT_0269 + IT_0251*IT_0270)*IT_0271 + (IT_0207*conj(IT_0207)
       + IT_0234*conj(IT_0234) + IT_0251*conj(IT_0251) + IT_0268*IT_0270)
      *IT_0272 + IT_0270*IT_0273 + 6*IT_0274*IT_0283 + (-18)*IT_0284*(IT_0207*
      (conj(IT_0251) + 1./3*conj(IT_0252) + 1./3*conj(IT_0268)) + (
      -0.0555555555555555)*IT_0251*IT_0287 + (-0.0555555555555555)*IT_0270
      *IT_0288 + (-0.0555555555555555)*IT_0269*(IT_0289 + IT_0290) + (
      -0.0555555555555555)*conj(IT_0251)*IT_0291) + IT_0292*IT_0305;
    return create_ccomplex_return(IT_0306);
}

