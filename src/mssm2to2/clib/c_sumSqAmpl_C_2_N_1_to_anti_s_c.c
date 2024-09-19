#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_1_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_1_to_anti_s_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = s_14*s_23;
    const ccomplex_t IT_0002 = s_13*s_24;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0000 + IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = pow(m_s, 2);
    const ccomplex_t IT_0006 = pow(m_N_1, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_B1*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = N_W1*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = sin(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, 2);
    const ccomplex_t IT_0029 = pow(m_C_2, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0006 + IT_0028 + -IT_0029 
      + -reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0014*IT_0020;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0014*IT_0021;
    const ccomplex_t IT_0034 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 2*IT_0035;
    const ccomplex_t IT_0037 = IT_0009*IT_0021;
    const ccomplex_t IT_0038 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = IT_0040 + IT_0043;
    const ccomplex_t IT_0045 = IT_0036 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0032*IT_0047;
    const ccomplex_t IT_0049 = IT_0030*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + -IT_0006 + -IT_0029 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0052 = cos(beta);
    const ccomplex_t IT_0053 = IT_0014*IT_0052;
    const ccomplex_t IT_0054 = IT_0034*IT_0053;
    const ccomplex_t IT_0055 = 2*IT_0054;
    const ccomplex_t IT_0056 = IT_0009*IT_0052;
    const ccomplex_t IT_0057 = IT_0038*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = IT_0041*IT_0053;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = IT_0055 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0014*IT_0020*IT_0022*IT_0052;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = IT_0051*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0050 + IT_0069;
    const ccomplex_t IT_0071 = 6*IT_0070;
    const ccomplex_t IT_0072 = 6*conj(IT_0070);
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0074 = cpow(IT_0052, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0014*IT_0020*IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0009;
    const ccomplex_t IT_0078 = (-0.333333333333333)*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = IT_0073*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0014;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0014*IT_0020*IT_0074;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = IT_0007*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0014*IT_0020;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0092 = IT_0053*IT_0091;
    const ccomplex_t IT_0093 = (-2)*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0095 = IT_0056*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0098 = IT_0053*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = IT_0096 + IT_0099;
    const ccomplex_t IT_0101 = IT_0093 + IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0090*IT_0103;
    const ccomplex_t IT_0105 = IT_0030*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0014*IT_0020*IT_0021*IT_0074;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0033*IT_0097;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = IT_0037*IT_0094;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = IT_0110 + IT_0112;
    const ccomplex_t IT_0114 = IT_0033*IT_0091;
    const ccomplex_t IT_0115 = (-2)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0108*IT_0118;
    const ccomplex_t IT_0120 = IT_0051*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = -IT_0106 + -IT_0121;
    const ccomplex_t IT_0123 = 6*IT_0122;
    const ccomplex_t IT_0124 = 6*conj(IT_0122);
    const ccomplex_t IT_0125 = 12*IT_0070;
    const ccomplex_t IT_0126 = IT_0064*IT_0108;
    const ccomplex_t IT_0127 = IT_0051*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0047*IT_0090;
    const ccomplex_t IT_0130 = IT_0030*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = 12*IT_0132;
    const ccomplex_t IT_0134 = IT_0032*IT_0103;
    const ccomplex_t IT_0135 = IT_0030*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0066*IT_0118;
    const ccomplex_t IT_0138 = IT_0051*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = -IT_0136 + -IT_0139;
    const ccomplex_t IT_0141 = 12*IT_0140;
    const ccomplex_t IT_0142 = 12*IT_0122;
    const ccomplex_t IT_0143 = m_s*s_14*m_N_1;
    const ccomplex_t IT_0144 = cpow((-2)*s_13 + IT_0005 + IT_0029 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0145 = conj(N_B1)*e_em;
    const ccomplex_t IT_0146 = IT_0009*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_W1)*e_em;
    const ccomplex_t IT_0149 = IT_0014*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 3*IT_0150);
    const ccomplex_t IT_0152 = 0.166666666666667*IT_0151;
    const ccomplex_t IT_0153 = IT_0076*IT_0152;
    const ccomplex_t IT_0154 = IT_0144*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0157 = IT_0014*IT_0156;
    const ccomplex_t IT_0158 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0159 = IT_0014*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + 0.5*IT_0160);
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0014;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0162*IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0167 = IT_0165*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0024*IT_0084;
    const ccomplex_t IT_0170 = IT_0007*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0009;
    const ccomplex_t IT_0174 = 0.666666666666667*IT_0173;
    const ccomplex_t IT_0175 = IT_0024*IT_0174;
    const ccomplex_t IT_0176 = cpow((-2)*s_13 + IT_0005 + IT_0029 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0177 = IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0014;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0180*IT_0182;
    const ccomplex_t IT_0184 = IT_0144*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0178 + IT_0185;
    const ccomplex_t IT_0187 = 6*IT_0186;
    const ccomplex_t IT_0188 = IT_0071 + IT_0187;
    const ccomplex_t IT_0189 = 6*conj(IT_0186);
    const ccomplex_t IT_0190 = IT_0072 + IT_0189;
    const ccomplex_t IT_0191 = IT_0152*IT_0180;
    const ccomplex_t IT_0192 = IT_0144*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = IT_0024*IT_0195;
    const ccomplex_t IT_0197 = IT_0176*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = -IT_0193 + -IT_0198;
    const ccomplex_t IT_0200 = 6*IT_0132;
    const ccomplex_t IT_0201 = 6*conj(IT_0132);
    const ccomplex_t IT_0202 = IT_0076*IT_0182;
    const ccomplex_t IT_0203 = IT_0144*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = 6*IT_0205;
    const ccomplex_t IT_0207 = 6*conj(IT_0205);
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0014*IT_0020*IT_0074;
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = IT_0076*IT_0209;
    const ccomplex_t IT_0211 = IT_0073*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0019*IT_0082;
    const ccomplex_t IT_0214 = IT_0007*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = -IT_0212 + -IT_0215;
    const ccomplex_t IT_0217 = 6*IT_0216;
    const ccomplex_t IT_0218 = 6*conj(IT_0216);
    const ccomplex_t IT_0219 = 12*IT_0172;
    const ccomplex_t IT_0220 = 12*conj(IT_0172);
    const ccomplex_t IT_0221 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0222 = IT_0014*IT_0221;
    const ccomplex_t IT_0223 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0224 = IT_0014*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + (-0.5)*IT_0225);
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0164*IT_0227;
    const ccomplex_t IT_0229 = IT_0166*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = 12*conj(IT_0132);
    const ccomplex_t IT_0232 = 12*IT_0216;
    const ccomplex_t IT_0233 = 12*conj(IT_0216);
    const ccomplex_t IT_0234 = 12*IT_0155;
    const ccomplex_t IT_0235 = (24*conj(IT_0088) + 12*conj(IT_0122) + 12*conj
      (IT_0155))*IT_0168 + conj(IT_0172)*IT_0188 + IT_0172*IT_0190 + conj
      (IT_0199)*IT_0200 + IT_0199*IT_0201 + conj(IT_0140)*IT_0206 + IT_0140
      *IT_0207 + (conj(IT_0122) + conj(IT_0155))*IT_0217 + (IT_0122 + IT_0155)
      *IT_0218 + conj(IT_0027)*IT_0219 + IT_0027*IT_0220 + IT_0133*conj(IT_0230)
       + IT_0230*IT_0231 + conj(IT_0088)*IT_0232 + IT_0088*IT_0233 + conj
      (IT_0168)*(24*IT_0088 + IT_0142 + IT_0234);
    const ccomplex_t IT_0236 = 12*IT_0199;
    const ccomplex_t IT_0237 = 24*IT_0230;
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = 12*IT_0186;
    const ccomplex_t IT_0240 = 12*IT_0205;
    const ccomplex_t IT_0241 = 12*IT_0088;
    const ccomplex_t IT_0242 = 12*IT_0027;
    const ccomplex_t IT_0243 = m_c*m_s*s_12;
    const ccomplex_t IT_0244 = 6*IT_0199;
    const ccomplex_t IT_0245 = 6*conj(IT_0199);
    const ccomplex_t IT_0246 = 6*IT_0140;
    const ccomplex_t IT_0247 = 6*conj(IT_0140);
    const ccomplex_t IT_0248 = 12*conj(IT_0205);
    const ccomplex_t IT_0249 = (-12)*IT_0070;
    const ccomplex_t IT_0250 = (-12)*conj(IT_0070);
    const ccomplex_t IT_0251 = (-12)*IT_0140;
    const ccomplex_t IT_0252 = (-12)*conj(IT_0140);
    const ccomplex_t IT_0253 = (-6)*IT_0132;
    const ccomplex_t IT_0254 = (-6)*conj(IT_0132);
    const ccomplex_t IT_0255 = (-6)*IT_0140;
    const ccomplex_t IT_0256 = (-6)*conj(IT_0140);
    const ccomplex_t IT_0257 = conj(IT_0186)*IT_0200 + IT_0186*IT_0201 +
       IT_0207*IT_0216 + IT_0206*conj(IT_0216) + IT_0220*IT_0230 + IT_0219*conj
      (IT_0230) + conj(IT_0168)*IT_0240 + conj(IT_0172)*IT_0244 + IT_0172
      *IT_0245 + conj(IT_0155)*IT_0246 + IT_0155*IT_0247 + IT_0168*IT_0248 +
       conj(IT_0132)*IT_0249 + IT_0132*IT_0250 + conj(IT_0122)*IT_0251 + IT_0122
      *IT_0252 + conj(IT_0027)*IT_0253 + IT_0027*IT_0254 + conj(IT_0088)*IT_0255
       + IT_0088*IT_0256;
    const ccomplex_t IT_0258 = s_34*m_C_2*m_N_1;
    const ccomplex_t IT_0259 = 12*conj(IT_0199);
    const ccomplex_t IT_0260 = (-12)*IT_0132;
    const ccomplex_t IT_0261 = (-12)*conj(IT_0132);
    const ccomplex_t IT_0262 = conj(IT_0155)*IT_0200 + IT_0155*IT_0201 + conj
      (IT_0172)*IT_0206 + IT_0172*IT_0207 + conj(IT_0199)*IT_0217 + IT_0199
      *IT_0218 + conj(IT_0230)*IT_0232 + IT_0230*IT_0233 + conj(IT_0168)*IT_0238
       + conj(IT_0140)*(IT_0187 + IT_0249) + IT_0140*(IT_0189 + IT_0250) + conj
      (IT_0088)*IT_0253 + IT_0088*IT_0254 + conj(IT_0027)*IT_0255 + IT_0027
      *IT_0256 + IT_0168*(24*conj(IT_0230) + IT_0259) + conj(IT_0122)*IT_0260 +
       IT_0122*IT_0261;
    const ccomplex_t IT_0263 = m_s*s_24*m_C_2;
    const ccomplex_t IT_0264 = (-6)*IT_0070;
    const ccomplex_t IT_0265 = IT_0239 + IT_0264;
    const ccomplex_t IT_0266 = 12*conj(IT_0186);
    const ccomplex_t IT_0267 = (-6)*conj(IT_0070);
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = 12*conj(IT_0088);
    const ccomplex_t IT_0270 = (-6)*IT_0122;
    const ccomplex_t IT_0271 = (-6)*conj(IT_0122);
    const ccomplex_t IT_0272 = conj(IT_0027)*IT_0206 + IT_0027*IT_0207 + conj
      (IT_0155)*IT_0236 + conj(IT_0230)*((-12)*IT_0122 + 24*IT_0155 + IT_0241) +
       conj(IT_0088)*IT_0244 + IT_0088*IT_0245 + conj(IT_0216)*IT_0253 + IT_0216
      *IT_0254 + conj(IT_0172)*IT_0255 + IT_0172*IT_0256 + IT_0155*IT_0259 +
       conj(IT_0168)*IT_0260 + IT_0168*IT_0261 + conj(IT_0205)*IT_0265 + IT_0205
      *IT_0268 + IT_0230*((-12)*conj(IT_0122) + 24*conj(IT_0155) + IT_0269) +
       conj(IT_0199)*IT_0270 + IT_0199*IT_0271;
    const ccomplex_t IT_0273 = -IT_0001;
    const ccomplex_t IT_0274 = IT_0000 + IT_0273;
    const ccomplex_t IT_0275 = IT_0002 + IT_0274;
    const ccomplex_t IT_0276 = m_c*s_23*m_C_2;
    const ccomplex_t IT_0277 = 12*conj(IT_0070);
    const ccomplex_t IT_0278 = 0.166666666666667*conj(IT_0199);
    const ccomplex_t IT_0279 = 0.166666666666667*IT_0088;
    const ccomplex_t IT_0280 = 0.166666666666667*conj(IT_0088);
    const ccomplex_t IT_0281 = 0.166666666666667*IT_0027;
    const ccomplex_t IT_0282 = 0.166666666666667*conj(IT_0027);
    const ccomplex_t IT_0283 = 0.5*(0.166666666666667*IT_0124 + conj(IT_0155))
      *IT_0172 + 0.5*(0.166666666666667*IT_0123 + IT_0155)*conj(IT_0172) +
       0.0833333333333333*IT_0201*IT_0205 + 0.0833333333333333*IT_0200*conj
      (IT_0205) + 0.0833333333333333*IT_0190*IT_0216 + 0.0833333333333333
      *IT_0188*conj(IT_0216) + conj(IT_0140)*IT_0230 + 0.0833333333333333
      *IT_0141*conj(IT_0230) + 2*conj(IT_0168)*(IT_0027 + 0.0416666666666667
      *IT_0125 + 0.0416666666666667*IT_0239) + 0.0833333333333333*IT_0199
      *IT_0247 + 2*IT_0168*(conj(IT_0027) + 0.0416666666666667*IT_0266 +
       0.0416666666666667*IT_0277) + 0.5*IT_0246*IT_0278 + 0.5*IT_0220*IT_0279 +
       0.5*IT_0219*IT_0280 + 0.5*IT_0233*IT_0281 + 0.5*IT_0232*IT_0282;
    const ccomplex_t IT_0284 = m_c*m_s*m_C_2*m_N_1;
    const ccomplex_t IT_0285 = (-6)*IT_0186;
    const ccomplex_t IT_0286 = (-6)*conj(IT_0186);
    const ccomplex_t IT_0287 = 0.0416666666666667*IT_0133*conj(IT_0140) + conj
      (IT_0205)*IT_0230 + IT_0205*conj(IT_0230) + 0.0416666666666667*IT_0140
      *IT_0231 + conj(IT_0172)*(IT_0168 + 0.0416666666666667*IT_0232) + IT_0172*
      (conj(IT_0168) + 0.0416666666666667*IT_0233) + 0.0416666666666667*IT_0199
      *IT_0248 + 0.25*conj(IT_0155)*(IT_0027 + 0.166666666666667*IT_0239 +
       0.166666666666667*IT_0264) + 0.25*IT_0155*(conj(IT_0027) +
       0.166666666666667*IT_0266 + 0.166666666666667*IT_0267) + 0.25*IT_0240
      *IT_0278 + 0.25*IT_0190*IT_0279 + 0.25*IT_0188*IT_0280 + 0.25*(IT_0124 +
       IT_0269)*IT_0281 + 0.25*(IT_0123 + IT_0241)*IT_0282 + 0.0416666666666667
      *conj(IT_0122)*(IT_0125 + IT_0285) + 0.0416666666666667*IT_0122*(IT_0277 +
       IT_0286);
    const ccomplex_t IT_0288 = m_c*s_13*m_N_1;
    const ccomplex_t IT_0289 = (-0.166666666666667)*IT_0027;
    const ccomplex_t IT_0290 = (-0.166666666666667)*conj(IT_0027);
    const ccomplex_t IT_0291 = conj(IT_0122)*IT_0205 + IT_0122*conj(IT_0205) +
       (-0.166666666666667)*conj(IT_0088)*IT_0206 + (-0.166666666666667)*IT_0088
      *IT_0207 + conj(IT_0140)*IT_0216 + IT_0140*conj(IT_0216) + (
      -0.166666666666667)*conj(IT_0155)*IT_0240 + (-0.166666666666667)*IT_0155
      *IT_0248 + (-4)*conj(IT_0230)*(IT_0186 + 0.0416666666666667*IT_0242 +
       0.0416666666666667*IT_0249) + (-4)*IT_0230*(0.5*conj(IT_0027) + conj
      (IT_0186) + 0.0416666666666667*IT_0250) + (-0.166666666666667)*conj
      (IT_0168)*IT_0251 + (-0.166666666666667)*IT_0168*IT_0252 + (
      -0.166666666666667)*conj(IT_0172)*IT_0253 + (-0.166666666666667)*IT_0172
      *IT_0254 + (-0.166666666666667)*conj(IT_0199)*IT_0265 + (
      -0.166666666666667)*IT_0199*IT_0268 + IT_0245*IT_0289 + IT_0244*IT_0290;
    const ccomplex_t IT_0292 = IT_0003 + IT_0274;
    const ccomplex_t IT_0293 = IT_0004*(conj(IT_0027)*IT_0071 + IT_0027
      *IT_0072 + conj(IT_0088)*IT_0123 + IT_0088*IT_0124) + IT_0000*(conj
      (IT_0070)*IT_0125 + conj(IT_0132)*IT_0133 + conj(IT_0140)*IT_0141 + conj
      (IT_0122)*IT_0142) + IT_0143*IT_0235 + IT_0002*((24*IT_0199 + 48*IT_0230)
      *conj(IT_0230) + conj(IT_0155)*IT_0234 + conj(IT_0199)*IT_0238 + conj
      (IT_0186)*IT_0239 + conj(IT_0205)*IT_0240) + IT_0001*(conj(IT_0168)*(48
      *IT_0168 + 24*IT_0216) + conj(IT_0172)*IT_0219 + conj(IT_0216)*(24*IT_0168
       + IT_0232) + conj(IT_0088)*IT_0241 + conj(IT_0027)*IT_0242) + IT_0243
      *IT_0257 + IT_0258*IT_0262 + IT_0263*IT_0272 + (conj(IT_0186)*IT_0264 +
       IT_0186*IT_0267 + conj(IT_0155)*IT_0270 + IT_0155*IT_0271)*IT_0275 + 12
      *IT_0276*IT_0283 + 24*IT_0284*IT_0287 + (-6)*IT_0288*IT_0291 + (-6)*(conj
      (IT_0088)*IT_0155 + IT_0088*conj(IT_0155) + IT_0286*IT_0289 + IT_0285
      *IT_0290)*IT_0292;
    return create_ccomplex_return(IT_0293);
}

