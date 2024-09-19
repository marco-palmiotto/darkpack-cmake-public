#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_se_L_to_e_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_se_L_to_e_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
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
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0005;
    const ccomplex_t IT_0009 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 2*IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = IT_0017 + IT_0020;
    const ccomplex_t IT_0022 = IT_0011 + IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = pow(m_e, 2);
    const ccomplex_t IT_0027 = pow(m_se_L, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0000*IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0034 = IT_0008*IT_0033;
    const ccomplex_t IT_0035 = 2*IT_0034;
    const ccomplex_t IT_0036 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0037 = IT_0014*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0040 = IT_0008*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = IT_0038 + IT_0041;
    const ccomplex_t IT_0043 = IT_0035 + IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0030*IT_0046*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B3)*e_em;
    const ccomplex_t IT_0050 = IT_0013*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W3)*e_em;
    const ccomplex_t IT_0053 = IT_0005*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + IT_0054);
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0045*IT_0056;
    const ccomplex_t IT_0058 = IT_0047*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = m_e*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0065 = IT_0008*IT_0064;
    const ccomplex_t IT_0066 = 2*IT_0065;
    const ccomplex_t IT_0067 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0014*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = IT_0066 + IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0063*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0061*IT_0077*IT_0078;
    const ccomplex_t IT_0080 = sin(beta);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0001*IT_0003*IT_0005*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0005;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = pow(m_C_2, 2);
    const ccomplex_t IT_0086 = cpow(s_12 + 0.5*IT_0027 + 0.5*IT_0085 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = m_C_2*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0095 = IT_0008*IT_0094;
    const ccomplex_t IT_0096 = 2*IT_0095;
    const ccomplex_t IT_0097 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0098 = IT_0014*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0101 = IT_0008*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = IT_0099 + IT_0102;
    const ccomplex_t IT_0104 = IT_0096 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0093*IT_0106;
    const ccomplex_t IT_0108 = IT_0090*IT_0091*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B1)*e_em;
    const ccomplex_t IT_0110 = IT_0013*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W1)*e_em;
    const ccomplex_t IT_0113 = IT_0005*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + IT_0114);
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0076*IT_0116;
    const ccomplex_t IT_0118 = IT_0078*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = m_e*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B2)*e_em;
    const ccomplex_t IT_0122 = IT_0013*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W2)*e_em;
    const ccomplex_t IT_0125 = IT_0005*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + IT_0126);
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = IT_0024*IT_0128;
    const ccomplex_t IT_0130 = IT_0028*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = m_e*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B4)*e_em;
    const ccomplex_t IT_0134 = IT_0013*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_W4)*e_em;
    const ccomplex_t IT_0137 = IT_0005*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + IT_0138);
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0106*IT_0140;
    const ccomplex_t IT_0142 = IT_0091*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = m_e*IT_0143;
    const ccomplex_t IT_0145 = -IT_0029 + -IT_0048 + IT_0060 + -IT_0079 + (
      -0.5)*IT_0089 + -IT_0108 + IT_0120 + IT_0132 + IT_0144;
    const ccomplex_t IT_0146 = IT_0005*IT_0080;
    const ccomplex_t IT_0147 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = (-2)*IT_0148;
    const ccomplex_t IT_0150 = IT_0013*IT_0080;
    const ccomplex_t IT_0151 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0155 = IT_0146*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = IT_0153 + IT_0156;
    const ccomplex_t IT_0158 = IT_0149 + IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0032*IT_0160;
    const ccomplex_t IT_0162 = IT_0047*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = m_e*IT_0163;
    const ccomplex_t IT_0165 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0166 = IT_0146*IT_0165;
    const ccomplex_t IT_0167 = (-2)*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0169 = IT_0150*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0172 = IT_0146*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = IT_0170 + IT_0173;
    const ccomplex_t IT_0175 = IT_0167 + IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = IT_0093*IT_0177;
    const ccomplex_t IT_0179 = IT_0091*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = m_e*IT_0180;
    const ccomplex_t IT_0182 = IT_0056*IT_0160;
    const ccomplex_t IT_0183 = IT_0030*IT_0047*IT_0182;
    const ccomplex_t IT_0184 = IT_0140*IT_0177;
    const ccomplex_t IT_0185 = IT_0090*IT_0091*IT_0184;
    const ccomplex_t IT_0186 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0187 = IT_0146*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0190 = IT_0150*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = IT_0188 + IT_0191;
    const ccomplex_t IT_0193 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0194 = IT_0146*IT_0193;
    const ccomplex_t IT_0195 = (-2)*IT_0194;
    const ccomplex_t IT_0196 = IT_0192 + IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = IT_0116*IT_0198;
    const ccomplex_t IT_0200 = IT_0061*IT_0078*IT_0199;
    const ccomplex_t IT_0201 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0202 = IT_0146*IT_0201;
    const ccomplex_t IT_0203 = (-2)*IT_0202;
    const ccomplex_t IT_0204 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0205 = IT_0150*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0208 = IT_0146*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = IT_0206 + IT_0209;
    const ccomplex_t IT_0211 = IT_0203 + IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = IT_0007*IT_0213;
    const ccomplex_t IT_0215 = IT_0028*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = m_e*IT_0216;
    const ccomplex_t IT_0218 = IT_0128*IT_0213;
    const ccomplex_t IT_0219 = IT_0000*IT_0028*IT_0218;
    const ccomplex_t IT_0220 = IT_0063*IT_0198;
    const ccomplex_t IT_0221 = IT_0078*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = m_e*IT_0222;
    const ccomplex_t IT_0224 = IT_0164 + IT_0181 + -IT_0183 + -IT_0185 + 
      -IT_0200 + IT_0217 + -IT_0219 + IT_0223;
    const ccomplex_t IT_0225 = 2*s_13;
    const ccomplex_t IT_0226 = m_e*m_C_2;
    const ccomplex_t IT_0227 = 2*IT_0145;
    const ccomplex_t IT_0228 = -IT_0059 + -IT_0119 + -IT_0131 + -IT_0143;
    const ccomplex_t IT_0229 = m_e*s_12;
    const ccomplex_t IT_0230 = (-0.5)*IT_0088 + -IT_0163 + -IT_0180 + -IT_0216
       + -IT_0222;
    const ccomplex_t IT_0231 = s_23*m_C_2;
    const ccomplex_t IT_0232 = 2*conj(IT_0145);
    const ccomplex_t IT_0233 = 2*conj(IT_0224);
    const ccomplex_t IT_0234 = m_e*m_C_2*IT_0027;
    const ccomplex_t IT_0235 = s_12*s_23;
    const ccomplex_t IT_0236 = s_13*IT_0027;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = 2*IT_0224;
    const ccomplex_t IT_0240 = 0.5*IT_0231;
    const ccomplex_t IT_0241 = (IT_0145*conj(IT_0145) + IT_0224*conj(IT_0224))
      *IT_0225 + conj(IT_0224)*IT_0226*IT_0227 + (IT_0224*IT_0226 + IT_0228
      *IT_0229 + IT_0230*IT_0231)*IT_0232 + (IT_0229*IT_0230 + IT_0228*IT_0231)
      *IT_0233 + 2*conj(IT_0230)*(IT_0228*IT_0234 + 2*IT_0230*IT_0238 + 0.5
      *IT_0229*IT_0239 + IT_0227*IT_0240) + 2*conj(IT_0228)*(0.5*IT_0227*IT_0229
       + IT_0230*IT_0234 + 2*IT_0228*IT_0238 + IT_0239*IT_0240);
    return create_ccomplex_return(IT_0241);
}

