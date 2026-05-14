#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_st_1_to_u_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_st_1_to_u_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_su_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = conj(N_B1)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = conj(N_W1)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0013 + 3*IT_0018);
    const ccomplex_t IT_0020 = 0.166666666666667*IT_0019;
    const ccomplex_t IT_0021 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = sin(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = IT_0015*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = m_t*N_u1*e_em*IT_0027*U_st_00;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0023 + (-0.75)*IT_0030);
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0031;
    const ccomplex_t IT_0033 = IT_0020*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B3)*e_em;
    const ccomplex_t IT_0038 = IT_0010*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W3)*e_em;
    const ccomplex_t IT_0041 = IT_0015*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 3*IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = m_t*N_u3*e_em*IT_0027*U_st_00;
    const ccomplex_t IT_0049 = IT_0026*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + (-0.75)*IT_0050);
    const ccomplex_t IT_0052 = 0.666666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_B4)*e_em;
    const ccomplex_t IT_0058 = IT_0010*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W4)*e_em;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + 3*IT_0062);
    const ccomplex_t IT_0064 = 0.166666666666667*IT_0063;
    const ccomplex_t IT_0065 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0066 = IT_0010*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = m_t*N_u4*e_em*IT_0027*U_st_00;
    const ccomplex_t IT_0069 = IT_0026*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + (-0.75)*IT_0070);
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0078 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0079 = IT_0010*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = m_t*N_u2*e_em*IT_0027*U_st_00;
    const ccomplex_t IT_0082 = IT_0026*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + (-0.75)*IT_0083);
    const ccomplex_t IT_0085 = 0.666666666666667*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B2)*e_em;
    const ccomplex_t IT_0087 = IT_0010*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W2)*e_em;
    const ccomplex_t IT_0090 = IT_0015*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0085*IT_0093;
    const ccomplex_t IT_0095 = IT_0077*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0008 + -IT_0036 + -IT_0056
       + -IT_0076 + -IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0008;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0015*IT_0025*IT_0027;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0102 = IT_0010*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0105 = IT_0015*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_t*conj(N_u1)*e_em*IT_0027*U_st_10;
    const ccomplex_t IT_0108 = IT_0026*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + 3*IT_0106 + 3
      *IT_0109);
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0111;
    const ccomplex_t IT_0113 = IT_0034*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0015*IT_0025*IT_0027;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0118 = IT_0010*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0121 = IT_0015*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = m_t*conj(N_u3)*e_em*IT_0027*U_st_10;
    const ccomplex_t IT_0124 = IT_0026*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0119 + 3*IT_0122 + 3
      *IT_0125);
    const ccomplex_t IT_0127 = 0.166666666666667*IT_0126;
    const ccomplex_t IT_0128 = IT_0116*IT_0127;
    const ccomplex_t IT_0129 = IT_0054*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0015*IT_0025*IT_0027;
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0134 = IT_0010*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0137 = IT_0015*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = m_t*conj(N_u4)*e_em*IT_0027*U_st_10;
    const ccomplex_t IT_0140 = IT_0026*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0135 + 3*IT_0138 + 3
      *IT_0141);
    const ccomplex_t IT_0143 = 0.166666666666667*IT_0142;
    const ccomplex_t IT_0144 = IT_0132*IT_0143;
    const ccomplex_t IT_0145 = IT_0074*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0148 = IT_0010*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0151 = IT_0015*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = m_t*conj(N_u2)*e_em*IT_0027*U_st_10;
    const ccomplex_t IT_0154 = IT_0026*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0149 + 3*IT_0152 + 3
      *IT_0155);
    const ccomplex_t IT_0157 = 0.166666666666667*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0015*IT_0025*IT_0027;
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0157*IT_0159;
    const ccomplex_t IT_0161 = IT_0077*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = -IT_0114 + -IT_0130 + -IT_0146 + -IT_0162;
    const ccomplex_t IT_0164 = m_t*m_u*IT_0005;
    const ccomplex_t IT_0165 = s_34*IT_0005;
    const ccomplex_t IT_0166 = (-18)*IT_0165;
    const ccomplex_t IT_0167 = s_13*s_14;
    const ccomplex_t IT_0168 = 36*IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = m_t*m_u;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0172 = IT_0020*IT_0111;
    const ccomplex_t IT_0173 = IT_0034*IT_0171*IT_0172;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0176 = IT_0064*IT_0143;
    const ccomplex_t IT_0177 = IT_0074*IT_0175*IT_0176;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = conj(IT_0174) + conj(IT_0178);
    const ccomplex_t IT_0180 = m_u*IT_0008;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = 6*IT_0181;
    const ccomplex_t IT_0183 = IT_0032*IT_0100;
    const ccomplex_t IT_0184 = IT_0034*IT_0171*IT_0183;
    const ccomplex_t IT_0185 = IT_0052*IT_0116;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0187 = IT_0054*IT_0185*IT_0186;
    const ccomplex_t IT_0188 = IT_0072*IT_0132;
    const ccomplex_t IT_0189 = IT_0074*IT_0175*IT_0188;
    const ccomplex_t IT_0190 = m_u*IT_0036;
    const ccomplex_t IT_0191 = m_u*IT_0056;
    const ccomplex_t IT_0192 = m_u*IT_0076;
    const ccomplex_t IT_0193 = m_u*IT_0096;
    const ccomplex_t IT_0194 = 0.166666666666667*IT_0180 + IT_0184 + IT_0187 +
       IT_0189 + -IT_0190 + -IT_0191 + -IT_0192 + -IT_0193;
    const ccomplex_t IT_0195 = 18*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0197 = IT_0085*IT_0159;
    const ccomplex_t IT_0198 = IT_0077*IT_0196*IT_0197;
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = (-18)*IT_0199;
    const ccomplex_t IT_0201 = IT_0182 + IT_0195 + IT_0200;
    const ccomplex_t IT_0202 = IT_0174 + IT_0178;
    const ccomplex_t IT_0203 = 6*conj(IT_0181);
    const ccomplex_t IT_0204 = 18*conj(IT_0194);
    const ccomplex_t IT_0205 = (-18)*conj(IT_0199);
    const ccomplex_t IT_0206 = IT_0203 + IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0208 = IT_0000*IT_0207;
    const ccomplex_t IT_0209 = IT_0006*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = m_sG*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = 6*IT_0199;
    const ccomplex_t IT_0214 = (-18)*IT_0181;
    const ccomplex_t IT_0215 = (-6)*IT_0194;
    const ccomplex_t IT_0216 = IT_0213 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = 6*conj(IT_0199);
    const ccomplex_t IT_0218 = (-18)*conj(IT_0181);
    const ccomplex_t IT_0219 = (-6)*conj(IT_0194);
    const ccomplex_t IT_0220 = IT_0217 + IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = IT_0044*IT_0127;
    const ccomplex_t IT_0222 = IT_0054*IT_0186*IT_0221;
    const ccomplex_t IT_0223 = m_u*IT_0114;
    const ccomplex_t IT_0224 = m_u*IT_0130;
    const ccomplex_t IT_0225 = m_u*IT_0146;
    const ccomplex_t IT_0226 = IT_0093*IT_0157;
    const ccomplex_t IT_0227 = IT_0077*IT_0196*IT_0226;
    const ccomplex_t IT_0228 = m_u*IT_0162;
    const ccomplex_t IT_0229 = (-0.166666666666667)*IT_0211 + IT_0222 + 
      -IT_0223 + -IT_0224 + -IT_0225 + IT_0227 + -IT_0228;
    const ccomplex_t IT_0230 = 18*IT_0199;
    const ccomplex_t IT_0231 = (-18)*IT_0194;
    const ccomplex_t IT_0232 = (-6)*IT_0181;
    const ccomplex_t IT_0233 = IT_0230 + IT_0231 + IT_0232;
    const ccomplex_t IT_0234 = 18*conj(IT_0199);
    const ccomplex_t IT_0235 = (-18)*conj(IT_0194);
    const ccomplex_t IT_0236 = (-6)*conj(IT_0181);
    const ccomplex_t IT_0237 = IT_0234 + IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = 6*s_34;
    const ccomplex_t IT_0239 = 18*s_34;
    const ccomplex_t IT_0240 = (-18)*s_34;
    const ccomplex_t IT_0241 = IT_0178*IT_0239;
    const ccomplex_t IT_0242 = IT_0229*IT_0240;
    const ccomplex_t IT_0243 = conj(IT_0229)*IT_0240;
    const ccomplex_t IT_0244 = (-6)*s_34;
    const ccomplex_t IT_0245 = 12*IT_0167;
    const ccomplex_t IT_0246 = (-6)*IT_0165;
    const ccomplex_t IT_0247 = m_t*s_13;
    const ccomplex_t IT_0248 = (-3)*IT_0178;
    const ccomplex_t IT_0249 = 3*conj(IT_0229);
    const ccomplex_t IT_0250 = (-3)*conj(IT_0174);
    const ccomplex_t IT_0251 = (-3)*conj(IT_0178);
    const ccomplex_t IT_0252 = 0.0555555555555556*conj(IT_0163)*IT_0201 +
       0.0555555555555556*IT_0163*IT_0206 + conj(IT_0098)*((-0.333333333333333)
      *IT_0174 + (-0.333333333333333)*IT_0178 + IT_0212 + 0.333333333333333
      *IT_0229) + IT_0098*((-0.333333333333333)*conj(IT_0174) + (
      -0.333333333333333)*conj(IT_0178) + conj(IT_0212) + 0.333333333333333*conj
      (IT_0229)) + (-1)*conj(IT_0097)*(IT_0174 + -1./3*IT_0212 + -IT_0229 + -1.
      /3*IT_0248) + 0.333333333333333*IT_0097*(conj(IT_0212) + IT_0249 + IT_0250
       + IT_0251);
    const ccomplex_t IT_0253 = m_u*s_14;
    const ccomplex_t IT_0254 = conj(IT_0098)*IT_0216 + IT_0098*IT_0220 + conj
      (IT_0097)*IT_0233 + IT_0097*IT_0237 + 18*conj(IT_0163)*(IT_0174 + -1./3
      *IT_0212 + -IT_0229 + -1./3*IT_0248) + (-6)*IT_0163*(conj(IT_0212) +
       IT_0249 + IT_0250 + IT_0251);
    const ccomplex_t IT_0255 = (((-18)*conj(IT_0097) + (-6)*conj(IT_0098))
      *IT_0163 + ((-18)*IT_0097 + (-6)*IT_0098)*conj(IT_0163))*IT_0164 + 
      (IT_0097*conj(IT_0097) + IT_0098*conj(IT_0098) + IT_0163*conj(IT_0163))
      *IT_0169 + IT_0170*(IT_0179*IT_0201 + IT_0202*IT_0206 + conj(IT_0212)
      *IT_0216 + IT_0212*IT_0220 + conj(IT_0229)*IT_0233 + IT_0229*IT_0237) + 
      (conj(IT_0181)*IT_0194 + IT_0181*conj(IT_0194) + conj(IT_0212)*IT_0229 +
       IT_0212*conj(IT_0229))*IT_0238 + (IT_0174*IT_0179 + IT_0181*conj(IT_0181)
       + IT_0194*conj(IT_0194) + IT_0199*conj(IT_0199) + IT_0212*conj(IT_0212) +
       IT_0229*conj(IT_0229))*IT_0239 + (conj(IT_0194)*IT_0199 + IT_0194*conj
      (IT_0199))*IT_0240 + IT_0179*(IT_0241 + IT_0242) + IT_0202*IT_0243 + (conj
      (IT_0181)*IT_0199 + IT_0181*conj(IT_0199) + IT_0179*IT_0212 + IT_0202*conj
      (IT_0212))*IT_0244 + (conj(IT_0097)*IT_0098 + IT_0097*conj(IT_0098))*
      (IT_0245 + IT_0246) + 18*IT_0247*IT_0252 + IT_0253*IT_0254;
    return create_ccomplex_return(IT_0255);
}

