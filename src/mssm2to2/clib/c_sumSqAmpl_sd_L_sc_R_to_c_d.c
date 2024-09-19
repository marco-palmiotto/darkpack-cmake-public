#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_sc_R_to_c_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_sc_R_to_c_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_L = param->m_sd_L;
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sc_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = m_d*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0013 = 0.666666666666667*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0026 = 0.666666666666667*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = -IT_0024 + -IT_0032;
    const ccomplex_t IT_0034 = conj(N_B2)*e_em;
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W2)*e_em;
    const ccomplex_t IT_0038 = IT_0018*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-3)*IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = IT_0023*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = m_d*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B3)*e_em;
    const ccomplex_t IT_0047 = IT_0011*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W3)*e_em;
    const ccomplex_t IT_0050 = IT_0018*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + (-3)*IT_0051);
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = IT_0026*IT_0053;
    const ccomplex_t IT_0055 = IT_0031*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = m_d*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B4)*e_em;
    const ccomplex_t IT_0059 = IT_0011*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W4)*e_em;
    const ccomplex_t IT_0062 = IT_0018*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + (-3)*IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0067 = 0.666666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = m_d*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0077*IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0014*IT_0016*IT_0018;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0067*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0085 = IT_0069*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B1)*e_em;
    const ccomplex_t IT_0087 = IT_0011*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em;
    const ccomplex_t IT_0090 = IT_0018*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-3)*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0074*IT_0093;
    const ccomplex_t IT_0095 = IT_0079*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = m_d*IT_0096;
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0008 + -IT_0045 + -IT_0057
       + -IT_0072 + IT_0080 + IT_0085 + -IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*V_cd*e_em*conj(U_Wm1)*IT_0018;
    const ccomplex_t IT_0100 = sin(beta);
    const ccomplex_t IT_0101 = cpow(IT_0100, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u1)*e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0099*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0106 = pow(m_sd_L, 2);
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0000 + IT_0106 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0104*IT_0105*IT_0107;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = IT_0065*IT_0111;
    const ccomplex_t IT_0113 = IT_0069*IT_0084*IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = IT_0053*IT_0116;
    const ccomplex_t IT_0118 = IT_0030*IT_0031*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0093*IT_0120;
    const ccomplex_t IT_0122 = IT_0078*IT_0079*IT_0121;
    const ccomplex_t IT_0123 = -IT_0118 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0114 + IT_0123;
    const ccomplex_t IT_0125 = conj(IT_0114) + conj(IT_0123);
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0020*IT_0127;
    const ccomplex_t IT_0129 = IT_0023*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0076*IT_0120;
    const ccomplex_t IT_0132 = IT_0079*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0130 + IT_0133;
    const ccomplex_t IT_0135 = m_d*IT_0134;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u1
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = IT_0103*IT_0138;
    const ccomplex_t IT_0140 = IT_0107*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_13 + IT_0000 + IT_0106 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0014*IT_0018*IT_0101;
    const ccomplex_t IT_0146 = (-0.5)*IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = IT_0142*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0141 + IT_0149;
    const ccomplex_t IT_0151 = m_c*IT_0150;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*V_cd*e_em*conj(U_Wm2)*IT_0018;
    const ccomplex_t IT_0154 = IT_0146*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0156 = IT_0142*IT_0154*IT_0155;
    const ccomplex_t IT_0157 = IT_0152 + IT_0156;
    const ccomplex_t IT_0158 = conj(IT_0152) + conj(IT_0156);
    const ccomplex_t IT_0159 = 6*s_34;
    const ccomplex_t IT_0160 = IT_0041*IT_0127;
    const ccomplex_t IT_0161 = IT_0022*IT_0023*IT_0160;
    const ccomplex_t IT_0162 = IT_0028*IT_0116;
    const ccomplex_t IT_0163 = IT_0031*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = m_d*IT_0164;
    const ccomplex_t IT_0166 = IT_0082*IT_0111;
    const ccomplex_t IT_0167 = IT_0069*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = m_d*IT_0168;
    const ccomplex_t IT_0170 = IT_0161 + -IT_0165 + -IT_0169;
    const ccomplex_t IT_0171 = conj(IT_0136) + conj(IT_0170);
    const ccomplex_t IT_0172 = 18*s_34;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = cpow(IT_0009, 2);
    const ccomplex_t IT_0175 = (-18)*s_34;
    const ccomplex_t IT_0176 = IT_0114*IT_0172;
    const ccomplex_t IT_0177 = IT_0170*IT_0175;
    const ccomplex_t IT_0178 = conj(IT_0170)*IT_0175;
    const ccomplex_t IT_0179 = IT_0136 + IT_0170;
    const ccomplex_t IT_0180 = (-6)*s_34;
    const ccomplex_t IT_0181 = IT_0009*IT_0180;
    const ccomplex_t IT_0182 = IT_0159*IT_0170;
    const ccomplex_t IT_0183 = IT_0152*IT_0172;
    const ccomplex_t IT_0184 = IT_0114*IT_0180;
    const ccomplex_t IT_0185 = IT_0159*conj(IT_0170);
    const ccomplex_t IT_0186 = conj(IT_0114)*IT_0180;
    const ccomplex_t IT_0187 = m_c*m_d;
    const ccomplex_t IT_0188 = 6*IT_0033;
    const ccomplex_t IT_0189 = (-18)*IT_0009;
    const ccomplex_t IT_0190 = (-6)*IT_0098;
    const ccomplex_t IT_0191 = IT_0188 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 6*conj(IT_0033);
    const ccomplex_t IT_0193 = 18*IT_0009;
    const ccomplex_t IT_0194 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0195 = IT_0192 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = 18*IT_0033;
    const ccomplex_t IT_0197 = (-18)*IT_0098;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = 18*conj(IT_0098);
    const ccomplex_t IT_0200 = (-18)*conj(IT_0033);
    const ccomplex_t IT_0201 = IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = 6*IT_0009;
    const ccomplex_t IT_0203 = 18*IT_0098;
    const ccomplex_t IT_0204 = (-18)*IT_0033;
    const ccomplex_t IT_0205 = 18*conj(IT_0033);
    const ccomplex_t IT_0206 = (-18)*conj(IT_0098);
    const ccomplex_t IT_0207 = (-6)*IT_0009;
    const ccomplex_t IT_0208 = 6*IT_0098;
    const ccomplex_t IT_0209 = (-6)*IT_0033;
    const ccomplex_t IT_0210 = (-6)*IT_0109*(conj(IT_0033) + -conj(IT_0098) + 
      (-0.166666666666667)*IT_0189) + IT_0158*IT_0191 + IT_0157*IT_0195 +
       IT_0171*IT_0198 + IT_0124*IT_0201 + (conj(IT_0114) + conj(IT_0123) +
       IT_0136 + IT_0170)*IT_0202 + IT_0125*(IT_0203 + IT_0204) + IT_0179*
      (IT_0205 + IT_0206) + (IT_0114 + IT_0123 + conj(IT_0136) + conj(IT_0170))
      *IT_0207 + conj(IT_0109)*(IT_0193 + IT_0208 + IT_0209);
    const ccomplex_t IT_0211 = m_d*s_13;
    const ccomplex_t IT_0212 = -IT_0130 + -IT_0133 + -IT_0164 + -IT_0168;
    const ccomplex_t IT_0213 = -IT_0044 + -IT_0056;
    const ccomplex_t IT_0214 = -IT_0071 + -IT_0096;
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = 0.166666666666667*IT_0007;
    const ccomplex_t IT_0217 = conj(IT_0213) + conj(IT_0214);
    const ccomplex_t IT_0218 = 6*IT_0109;
    const ccomplex_t IT_0219 = 18*IT_0114;
    const ccomplex_t IT_0220 = 18*IT_0123;
    const ccomplex_t IT_0221 = (-18)*IT_0170;
    const ccomplex_t IT_0222 = (-18)*IT_0136;
    const ccomplex_t IT_0223 = (-6)*IT_0152;
    const ccomplex_t IT_0224 = (-0.5)*IT_0007;
    const ccomplex_t IT_0225 = (-6)*IT_0212;
    const ccomplex_t IT_0226 = (-6)*conj(IT_0212);
    const ccomplex_t IT_0227 = conj(IT_0150)*IT_0191 + IT_0150*IT_0195 + 6*
      (conj(IT_0109) + 3*conj(IT_0114) + 3*conj(IT_0123) + (-3)*conj(IT_0136) + 
      -conj(IT_0152) + (-3)*conj(IT_0170))*IT_0212 + (IT_0202 + IT_0205 +
       IT_0206)*IT_0215 + (12*IT_0009 + IT_0203 + IT_0204 + IT_0205 + IT_0206)
      *IT_0216 + (IT_0198 + IT_0207)*IT_0217 + conj(IT_0212)*(IT_0218 + IT_0219 
      + IT_0220 + IT_0221 + IT_0222 + IT_0223) + (36*IT_0009 + IT_0192 + IT_0194
       + IT_0208 + IT_0209)*IT_0224 + conj(IT_0156)*IT_0225 + IT_0156*IT_0226;
    const ccomplex_t IT_0228 = m_c*m_d*IT_0106;
    const ccomplex_t IT_0229 = (-18)*IT_0212;
    const ccomplex_t IT_0230 = IT_0213 + IT_0214 + IT_0216;
    const ccomplex_t IT_0231 = (-18)*conj(IT_0212);
    const ccomplex_t IT_0232 = m_c*s_14;
    const ccomplex_t IT_0233 = 6*conj(IT_0170);
    const ccomplex_t IT_0234 = 6*conj(IT_0136);
    const ccomplex_t IT_0235 = 18*conj(IT_0152);
    const ccomplex_t IT_0236 = (-18)*conj(IT_0109);
    const ccomplex_t IT_0237 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0238 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0239 = 6*IT_0152;
    const ccomplex_t IT_0240 = 18*IT_0170;
    const ccomplex_t IT_0241 = 18*IT_0136;
    const ccomplex_t IT_0242 = (-18)*IT_0114;
    const ccomplex_t IT_0243 = (-18)*IT_0123;
    const ccomplex_t IT_0244 = (-6)*IT_0109;
    const ccomplex_t IT_0245 = 6*conj(IT_0152);
    const ccomplex_t IT_0246 = 18*conj(IT_0170);
    const ccomplex_t IT_0247 = 18*conj(IT_0136);
    const ccomplex_t IT_0248 = (-18)*conj(IT_0114);
    const ccomplex_t IT_0249 = (-18)*conj(IT_0123);
    const ccomplex_t IT_0250 = (-6)*conj(IT_0109);
    const ccomplex_t IT_0251 = 3*IT_0156;
    const ccomplex_t IT_0252 = (-3)*IT_0109;
    const ccomplex_t IT_0253 = -IT_0123;
    const ccomplex_t IT_0254 = -IT_0114 + IT_0136 + 3*IT_0152 + (-3)*conj
      (IT_0156) + IT_0170 + (-0.166666666666667)*IT_0233 + (-0.166666666666667)
      *IT_0234 + (-0.166666666666667)*IT_0235 + (-0.166666666666667)*IT_0236 + (
      -0.166666666666667)*IT_0237 + (-0.166666666666667)*IT_0238 + IT_0251 +
       IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = -IT_0224;
    const ccomplex_t IT_0256 = IT_0156 + -conj(IT_0156) + (-0.166666666666667)
      *IT_0218 + (-0.166666666666667)*IT_0219 + (-0.166666666666667)*IT_0220 + (
      -0.166666666666667)*IT_0221 + (-0.166666666666667)*IT_0222 + (
      -0.166666666666667)*IT_0223 + (-0.166666666666667)*IT_0245 + (
      -0.166666666666667)*IT_0246 + (-0.166666666666667)*IT_0247 + (
      -0.166666666666667)*IT_0248 + (-0.166666666666667)*IT_0249 + (
      -0.166666666666667)*IT_0250;
    const ccomplex_t IT_0257 = -IT_0216;
    const ccomplex_t IT_0258 = (IT_0201 + IT_0207)*IT_0212 + (IT_0202 +
       IT_0203 + IT_0204)*conj(IT_0212) + IT_0150*(18*conj(IT_0156) + IT_0233 +
       IT_0234 + IT_0235 + IT_0236 + IT_0237 + IT_0238) + 6*IT_0217*(IT_0156 +
       0.166666666666667*IT_0239 + 0.166666666666667*IT_0240 + 0.166666666666667
      *IT_0241 + 0.166666666666667*IT_0242 + 0.166666666666667*IT_0243 +
       0.166666666666667*IT_0244) + 6*IT_0215*(conj(IT_0156) + 0.166666666666667
      *IT_0245 + 0.166666666666667*IT_0246 + 0.166666666666667*IT_0247 +
       0.166666666666667*IT_0248 + 0.166666666666667*IT_0249 + 0.166666666666667
      *IT_0250) + (-6)*conj(IT_0150)*(IT_0114 + -IT_0136 + (-3)*IT_0152 + 
      -IT_0170 + -IT_0251 + -IT_0252 + -IT_0253) + 6*IT_0254*IT_0255 + 6*IT_0256
      *IT_0257;
    const ccomplex_t IT_0259 = s_34*IT_0106;
    const ccomplex_t IT_0260 = (-18)*IT_0259;
    const ccomplex_t IT_0261 = s_13*s_14;
    const ccomplex_t IT_0262 = 36*IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = -IT_0150;
    const ccomplex_t IT_0265 = -conj(IT_0213);
    const ccomplex_t IT_0266 = -conj(IT_0214);
    const ccomplex_t IT_0267 = -IT_0214;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = (-6)*IT_0259;
    const ccomplex_t IT_0270 = 12*IT_0261;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = (IT_0009*(IT_0033 + conj(IT_0098)) + conj
      (IT_0109)*IT_0124 + IT_0109*IT_0125 + conj(IT_0136)*IT_0157 + IT_0136
      *IT_0158)*IT_0159 + (IT_0033*conj(IT_0033) + IT_0098*conj(IT_0098) +
       IT_0109*conj(IT_0109) + IT_0123*IT_0125 + IT_0156*IT_0158 + IT_0136
      *IT_0171)*IT_0172 + IT_0171*IT_0173 + (conj(IT_0033)*IT_0098 + IT_0033
      *conj(IT_0098) + IT_0125*IT_0136 + IT_0124*conj(IT_0136) + conj(IT_0109)
      *IT_0157 + IT_0109*IT_0158 + IT_0174)*IT_0175 + IT_0125*(IT_0176 + IT_0177
      ) + IT_0124*IT_0178 + (conj(IT_0123)*IT_0157 + IT_0123*IT_0158 + IT_0109
      *IT_0171 + conj(IT_0109)*IT_0179)*IT_0180 + (conj(IT_0033) + IT_0098)
      *IT_0181 + IT_0158*(IT_0182 + IT_0183 + IT_0184) + IT_0157*(IT_0185 +
       IT_0186) + IT_0187*IT_0210 + IT_0211*IT_0227 + IT_0228*(18*IT_0212
      *IT_0216 + conj(IT_0150)*IT_0225 + IT_0150*IT_0226 + 6*IT_0224*(IT_0212 +
       0.166666666666667*IT_0226) + IT_0217*IT_0229 + IT_0230*IT_0231) + IT_0232
      *IT_0258 + IT_0263*(IT_0150*conj(IT_0150) + IT_0212*conj(IT_0212) + 
      (IT_0213 + IT_0214)*IT_0217 + IT_0224*(conj(IT_0150) + IT_0255 + IT_0264) 
      + IT_0257*(IT_0213 + IT_0216 + IT_0265 + IT_0266 + IT_0268)) + (IT_0150
      *IT_0217 + conj(IT_0150)*IT_0230 + IT_0216*IT_0264 + IT_0255*(IT_0213 + 2
      *IT_0216 + IT_0265 + IT_0266 + IT_0268))*IT_0271;
    return create_ccomplex_return(IT_0272);
}

