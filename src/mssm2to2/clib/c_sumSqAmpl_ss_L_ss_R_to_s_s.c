#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_ss_R_to_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_ss_R_to_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_ss_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B4)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W4)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0009;
    const ccomplex_t IT_0021 = (-0.333333333333333)*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B3)*e_em;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W3)*e_em;
    const ccomplex_t IT_0030 = IT_0014*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-3)*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0009;
    const ccomplex_t IT_0035 = (-0.333333333333333)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0041 = conj(N_B1)*e_em;
    const ccomplex_t IT_0042 = IT_0009*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W1)*e_em;
    const ccomplex_t IT_0045 = IT_0014*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + (-3)*IT_0046);
    const ccomplex_t IT_0048 = 0.166666666666667*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0009;
    const ccomplex_t IT_0050 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0009;
    const ccomplex_t IT_0056 = (-0.333333333333333)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_B2)*e_em;
    const ccomplex_t IT_0058 = IT_0009*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W2)*e_em;
    const ccomplex_t IT_0061 = IT_0014*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + (-3)*IT_0062);
    const ccomplex_t IT_0064 = 0.166666666666667*IT_0063;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = 0.166666666666667*IT_0007 + -IT_0025 + -IT_0039
       + -IT_0053 + -IT_0067;
    const ccomplex_t IT_0069 = pow(m_ss_R, 2);
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0069 +
       reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0002*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0003 + IT_0069 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0036*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0003 + IT_0069 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0022*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0003 + IT_0069 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0051*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0003 + IT_0069 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0065*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (-0.166666666666667)*IT_0071 + IT_0074 +
       IT_0077 + IT_0080 + IT_0083;
    const ccomplex_t IT_0085 = pow(m_W, -1);
    const ccomplex_t IT_0086 = cos(beta);
    const ccomplex_t IT_0087 = cpow(IT_0086, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0089*IT_0091;
    const ccomplex_t IT_0093 = IT_0023*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0096*IT_0098;
    const ccomplex_t IT_0100 = IT_0054*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = IT_0037*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0014*IT_0085*IT_0087;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = IT_0040*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0071 + -IT_0094 + -IT_0101 + -IT_0108 + 
      -IT_0115;
    const ccomplex_t IT_0117 = IT_0081*IT_0099;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0072*IT_0106;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0075*IT_0092;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0078*IT_0113;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0007 + IT_0118 + IT_0120 + IT_0122 +
       IT_0124;
    const ccomplex_t IT_0126 = IT_0003*IT_0069;
    const ccomplex_t IT_0127 = (-6)*(conj(IT_0068)*(IT_0084 + 3*IT_0116) +
       IT_0068*(conj(IT_0084) + 3*conj(IT_0116)) + (-0.166666666666667)*((-18)
      *conj(IT_0084) + (-6)*conj(IT_0116))*IT_0125 + (-0.166666666666667)*((-18)
      *IT_0084 + (-6)*IT_0116)*conj(IT_0125))*IT_0126;
    const ccomplex_t IT_0128 = s_34*IT_0069;
    const ccomplex_t IT_0129 = (-6)*IT_0128;
    const ccomplex_t IT_0130 = s_23*s_24;
    const ccomplex_t IT_0131 = 12*IT_0130;
    const ccomplex_t IT_0132 = IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = (-18)*IT_0128;
    const ccomplex_t IT_0134 = 36*IT_0130;
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0137 = IT_0019*IT_0089;
    const ccomplex_t IT_0138 = IT_0023*IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0140 = IT_0064*IT_0096;
    const ccomplex_t IT_0141 = IT_0054*IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0143 = IT_0033*IT_0103;
    const ccomplex_t IT_0144 = IT_0037*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0146 = IT_0048*IT_0110;
    const ccomplex_t IT_0147 = IT_0040*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = -IT_0138 + -IT_0141 + -IT_0144 + -IT_0147;
    const ccomplex_t IT_0149 = IT_0118 + IT_0120 + IT_0122 + IT_0124;
    const ccomplex_t IT_0150 = m_s*IT_0149;
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = 6*IT_0151;
    const ccomplex_t IT_0153 = IT_0101 + IT_0115;
    const ccomplex_t IT_0154 = m_s*IT_0153;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = 18*IT_0155;
    const ccomplex_t IT_0157 = IT_0094 + IT_0108;
    const ccomplex_t IT_0158 = m_s*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = 18*IT_0159;
    const ccomplex_t IT_0161 = IT_0056*IT_0098;
    const ccomplex_t IT_0162 = IT_0054*IT_0139*IT_0161;
    const ccomplex_t IT_0163 = IT_0035*IT_0105;
    const ccomplex_t IT_0164 = IT_0037*IT_0142*IT_0163;
    const ccomplex_t IT_0165 = IT_0021*IT_0091;
    const ccomplex_t IT_0166 = IT_0023*IT_0136*IT_0165;
    const ccomplex_t IT_0167 = -IT_0162 + -IT_0164 + -IT_0166;
    const ccomplex_t IT_0168 = (-18)*IT_0167;
    const ccomplex_t IT_0169 = IT_0050*IT_0112;
    const ccomplex_t IT_0170 = IT_0040*IT_0145*IT_0169;
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = (-18)*IT_0171;
    const ccomplex_t IT_0173 = IT_0075*IT_0136*IT_0165;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = (-6)*IT_0174;
    const ccomplex_t IT_0176 = IT_0081*IT_0139*IT_0161;
    const ccomplex_t IT_0177 = IT_0078*IT_0145*IT_0169;
    const ccomplex_t IT_0178 = -IT_0176 + -IT_0177;
    const ccomplex_t IT_0179 = (-6)*IT_0178;
    const ccomplex_t IT_0180 = IT_0072*IT_0142*IT_0163;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = (-6)*IT_0181;
    const ccomplex_t IT_0183 = IT_0152 + IT_0156 + IT_0160 + IT_0168 + IT_0172
       + IT_0175 + IT_0179 + IT_0182;
    const ccomplex_t IT_0184 = 6*conj(IT_0151);
    const ccomplex_t IT_0185 = 18*conj(IT_0155);
    const ccomplex_t IT_0186 = 18*conj(IT_0159);
    const ccomplex_t IT_0187 = (-18)*conj(IT_0167);
    const ccomplex_t IT_0188 = (-18)*conj(IT_0171);
    const ccomplex_t IT_0189 = (-6)*conj(IT_0174);
    const ccomplex_t IT_0190 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0191 = (-6)*conj(IT_0181);
    const ccomplex_t IT_0192 = IT_0184 + IT_0185 + IT_0186 + IT_0187 + IT_0188
       + IT_0189 + IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = (-3)*IT_0007 + (-6)*IT_0074 + (-6)*IT_0077 + (
      -6)*IT_0080 + (-6)*IT_0083;
    const ccomplex_t IT_0194 = IT_0071 + IT_0193;
    const ccomplex_t IT_0195 = m_s*IT_0194;
    const ccomplex_t IT_0196 = 0.166666666666667*IT_0195;
    const ccomplex_t IT_0197 = 6*IT_0167;
    const ccomplex_t IT_0198 = 6*IT_0171;
    const ccomplex_t IT_0199 = 18*IT_0174;
    const ccomplex_t IT_0200 = 18*IT_0178;
    const ccomplex_t IT_0201 = 18*IT_0181;
    const ccomplex_t IT_0202 = (-18)*IT_0151;
    const ccomplex_t IT_0203 = (-6)*IT_0155;
    const ccomplex_t IT_0204 = (-6)*IT_0159;
    const ccomplex_t IT_0205 = IT_0197 + IT_0198 + IT_0199 + IT_0200 + IT_0201
       + IT_0202 + IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = 6*conj(IT_0167);
    const ccomplex_t IT_0207 = 6*conj(IT_0171);
    const ccomplex_t IT_0208 = 18*conj(IT_0174);
    const ccomplex_t IT_0209 = 18*conj(IT_0178);
    const ccomplex_t IT_0210 = 18*conj(IT_0181);
    const ccomplex_t IT_0211 = (-18)*conj(IT_0151);
    const ccomplex_t IT_0212 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0213 = (-6)*conj(IT_0159);
    const ccomplex_t IT_0214 = IT_0206 + IT_0207 + IT_0208 + IT_0209 + IT_0210
       + IT_0211 + IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = (-6)*IT_0025;
    const ccomplex_t IT_0216 = IT_0007 + IT_0215;
    const ccomplex_t IT_0217 = m_s*IT_0216;
    const ccomplex_t IT_0218 = 0.166666666666667*IT_0217;
    const ccomplex_t IT_0219 = 2*IT_0039 + 2*IT_0053 + 2*IT_0067;
    const ccomplex_t IT_0220 = IT_0071 + IT_0219;
    const ccomplex_t IT_0221 = m_s*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = conj(IT_0218) + conj(IT_0222);
    const ccomplex_t IT_0224 = 6*IT_0174;
    const ccomplex_t IT_0225 = 6*IT_0178;
    const ccomplex_t IT_0226 = 6*IT_0181;
    const ccomplex_t IT_0227 = 18*IT_0167;
    const ccomplex_t IT_0228 = 18*IT_0171;
    const ccomplex_t IT_0229 = (-18)*IT_0155;
    const ccomplex_t IT_0230 = (-18)*IT_0159;
    const ccomplex_t IT_0231 = (-6)*IT_0151;
    const ccomplex_t IT_0232 = IT_0224 + IT_0225 + IT_0226 + IT_0227 + IT_0228
       + IT_0229 + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0218 + IT_0222;
    const ccomplex_t IT_0234 = 6*conj(IT_0174);
    const ccomplex_t IT_0235 = 6*conj(IT_0178);
    const ccomplex_t IT_0236 = 6*conj(IT_0181);
    const ccomplex_t IT_0237 = 18*conj(IT_0167);
    const ccomplex_t IT_0238 = 18*conj(IT_0171);
    const ccomplex_t IT_0239 = (-18)*conj(IT_0155);
    const ccomplex_t IT_0240 = (-18)*conj(IT_0159);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0242 = IT_0234 + IT_0235 + IT_0236 + IT_0237 + IT_0238
       + IT_0239 + IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = IT_0072*IT_0142*IT_0143;
    const ccomplex_t IT_0244 = IT_0075*IT_0136*IT_0137;
    const ccomplex_t IT_0245 = IT_0078*IT_0145*IT_0146;
    const ccomplex_t IT_0246 = IT_0081*IT_0139*IT_0140;
    const ccomplex_t IT_0247 = -IT_0243 + -IT_0244 + -IT_0245 + -IT_0246;
    const ccomplex_t IT_0248 = 6*IT_0155;
    const ccomplex_t IT_0249 = 6*IT_0159;
    const ccomplex_t IT_0250 = 18*IT_0151;
    const ccomplex_t IT_0251 = (-18)*IT_0174;
    const ccomplex_t IT_0252 = (-18)*IT_0178;
    const ccomplex_t IT_0253 = (-18)*IT_0181;
    const ccomplex_t IT_0254 = (-6)*IT_0167;
    const ccomplex_t IT_0255 = (-6)*IT_0171;
    const ccomplex_t IT_0256 = IT_0248 + IT_0249 + IT_0250 + IT_0251 + IT_0252
       + IT_0253 + IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 6*conj(IT_0155);
    const ccomplex_t IT_0258 = 6*conj(IT_0159);
    const ccomplex_t IT_0259 = 18*conj(IT_0151);
    const ccomplex_t IT_0260 = (-18)*conj(IT_0174);
    const ccomplex_t IT_0261 = (-18)*conj(IT_0178);
    const ccomplex_t IT_0262 = (-18)*conj(IT_0181);
    const ccomplex_t IT_0263 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0264 = (-6)*conj(IT_0171);
    const ccomplex_t IT_0265 = IT_0257 + IT_0258 + IT_0259 + IT_0260 + IT_0261
       + IT_0262 + IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = IT_0174 + IT_0181;
    const ccomplex_t IT_0267 = conj(IT_0174) + conj(IT_0181);
    const ccomplex_t IT_0268 = IT_0155 + IT_0159;
    const ccomplex_t IT_0269 = conj(IT_0155) + conj(IT_0159);
    const ccomplex_t IT_0270 = 6*s_34;
    const ccomplex_t IT_0271 = conj(IT_0196)*IT_0270;
    const ccomplex_t IT_0272 = IT_0167 + IT_0171;
    const ccomplex_t IT_0273 = conj(IT_0167) + conj(IT_0171);
    const ccomplex_t IT_0274 = 18*s_34;
    const ccomplex_t IT_0275 = IT_0196*IT_0270;
    const ccomplex_t IT_0276 = IT_0218*IT_0274;
    const ccomplex_t IT_0277 = (-18)*s_34;
    const ccomplex_t IT_0278 = IT_0178*IT_0274;
    const ccomplex_t IT_0279 = IT_0151*IT_0277;
    const ccomplex_t IT_0280 = conj(IT_0151)*IT_0277;
    const ccomplex_t IT_0281 = IT_0178*IT_0270;
    const ccomplex_t IT_0282 = IT_0155*IT_0277;
    const ccomplex_t IT_0283 = (-6)*s_34;
    const ccomplex_t IT_0284 = conj(IT_0178)*IT_0270;
    const ccomplex_t IT_0285 = conj(IT_0155)*IT_0277;
    const ccomplex_t IT_0286 = conj(IT_0151)*IT_0283;
    const ccomplex_t IT_0287 = IT_0155*IT_0274;
    const ccomplex_t IT_0288 = IT_0178*IT_0283;
    const ccomplex_t IT_0289 = conj(IT_0178)*IT_0283;
    const ccomplex_t IT_0290 = IT_0167*IT_0270;
    const ccomplex_t IT_0291 = IT_0174*IT_0274;
    const ccomplex_t IT_0292 = IT_0155*IT_0283;
    const ccomplex_t IT_0293 = conj(IT_0167)*IT_0270;
    const ccomplex_t IT_0294 = conj(IT_0178)*IT_0274;
    const ccomplex_t IT_0295 = conj(IT_0155)*IT_0283;
    const ccomplex_t IT_0296 = m_s*s_23;
    const ccomplex_t IT_0297 = -IT_0148;
    const ccomplex_t IT_0298 = 0.333333333333333*IT_0196 + IT_0218 + IT_0222 +
       (-0.333333333333333)*IT_0247 + IT_0297;
    const ccomplex_t IT_0299 = 3*IT_0298;
    const ccomplex_t IT_0300 = -conj(IT_0148);
    const ccomplex_t IT_0301 = (-0.333333333333333)*conj(IT_0247);
    const ccomplex_t IT_0302 = 0.333333333333333*conj(IT_0196);
    const ccomplex_t IT_0303 = conj(IT_0218) + conj(IT_0222) + IT_0300 +
       IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = 3*IT_0303;
    const ccomplex_t IT_0305 = 3*conj(IT_0196);
    const ccomplex_t IT_0306 = (-3)*conj(IT_0247);
    const ccomplex_t IT_0307 = conj(IT_0084)*IT_0205 + IT_0084*IT_0214 + conj
      (IT_0116)*IT_0232 + IT_0116*IT_0242 + (-18)*conj(IT_0125)*(IT_0196 + 1./3
      *IT_0218 + 1./3*IT_0222 + -IT_0247 + 1./3*IT_0297) + (-6)*conj(IT_0068)
      *IT_0299 + (-6)*IT_0068*IT_0304 + (-6)*IT_0125*(conj(IT_0218) + conj
      (IT_0222) + IT_0300 + IT_0305 + IT_0306);
    const ccomplex_t IT_0308 = m_s*s_24;
    const ccomplex_t IT_0309 = conj(IT_0068)*IT_0183 + IT_0068*IT_0192 + conj
      (IT_0125)*IT_0256 + IT_0125*IT_0265 + 18*conj(IT_0084)*(IT_0196 + 1./3
      *IT_0218 + 1./3*IT_0222 + -IT_0247 + 1./3*IT_0297) + 6*conj(IT_0116)
      *IT_0299 + 6*IT_0116*IT_0304 + 6*IT_0084*(conj(IT_0218) + conj(IT_0222) +
       IT_0300 + IT_0305 + IT_0306);
    const ccomplex_t IT_0310 = IT_0127 + (conj(IT_0084)*IT_0116 + IT_0084*conj
      (IT_0116) + conj(IT_0068)*IT_0125 + IT_0068*conj(IT_0125))*IT_0132 + 
      (IT_0068*conj(IT_0068) + IT_0084*conj(IT_0084) + IT_0116*conj(IT_0116) +
       IT_0125*conj(IT_0125))*IT_0135 + IT_0003*(conj(IT_0148)*IT_0183 + IT_0148
      *IT_0192 + conj(IT_0196)*IT_0205 + IT_0196*IT_0214 + IT_0223*IT_0232 +
       IT_0233*IT_0242 + conj(IT_0247)*IT_0256 + IT_0247*IT_0265) + (conj
      (IT_0148)*IT_0247 + IT_0148*conj(IT_0247) + conj(IT_0171)*IT_0266 +
       IT_0171*IT_0267 + conj(IT_0151)*IT_0268 + IT_0151*IT_0269)*IT_0270 +
       IT_0233*IT_0271 + (IT_0148*conj(IT_0148) + IT_0151*conj(IT_0151) +
       IT_0196*conj(IT_0196) + IT_0222*IT_0223 + IT_0247*conj(IT_0247) + IT_0181
      *IT_0267 + IT_0159*IT_0269 + IT_0272*IT_0273)*IT_0274 + IT_0223*(IT_0275 +
       IT_0276) + (IT_0148*IT_0223 + conj(IT_0148)*IT_0233 + conj(IT_0196)
      *IT_0247 + IT_0196*conj(IT_0247) + conj(IT_0159)*IT_0272 + IT_0159*IT_0273
      )*IT_0277 + conj(IT_0178)*(IT_0278 + IT_0279) + IT_0178*IT_0280 + IT_0273*
      (IT_0281 + IT_0282) + (conj(IT_0148)*IT_0196 + IT_0148*conj(IT_0196) +
       IT_0223*IT_0247 + IT_0233*conj(IT_0247) + conj(IT_0159)*IT_0266 + IT_0159
      *IT_0267 + IT_0151*IT_0273)*IT_0283 + IT_0272*(IT_0284 + IT_0285 + IT_0286
      ) + IT_0269*(IT_0287 + IT_0288) + IT_0268*IT_0289 + IT_0267*(IT_0278 +
       IT_0279 + IT_0290 + IT_0291 + IT_0292) + IT_0266*(IT_0280 + IT_0293 +
       IT_0294 + IT_0295) + IT_0296*IT_0307 + IT_0308*IT_0309;
    return create_ccomplex_return(IT_0310);
}

