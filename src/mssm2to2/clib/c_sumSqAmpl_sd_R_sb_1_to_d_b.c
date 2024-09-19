#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sb_1_to_d_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sb_1_to_d_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_R = param->m_sd_R;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_d, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sd_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0013 = (-0.333333333333333)*IT_0012;
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = m_b*conj(N_d1)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0024 = IT_0011*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0027 = IT_0017*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0022 + 0.333333333333333
      *IT_0025 + -IT_0028);
    const ccomplex_t IT_0030 = 3*IT_0029;
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = IT_0013*IT_0031;
    const ccomplex_t IT_0033 = IT_0009*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0037 = (-0.333333333333333)*IT_0036;
    const ccomplex_t IT_0038 = m_b*conj(N_d2)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0039 = IT_0018*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0042 = IT_0011*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0045 = IT_0017*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0040 + 0.333333333333333
      *IT_0043 + -IT_0046);
    const ccomplex_t IT_0048 = 3*IT_0047;
    const ccomplex_t IT_0049 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_0050 = IT_0037*IT_0049;
    const ccomplex_t IT_0051 = IT_0035*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0054;
    const ccomplex_t IT_0056 = m_b*conj(N_d3)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0057 = IT_0018*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0060 = IT_0011*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + 0.333333333333333
      *IT_0061 + -IT_0064);
    const ccomplex_t IT_0066 = 3*IT_0065;
    const ccomplex_t IT_0067 = 0.166666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0055*IT_0067;
    const ccomplex_t IT_0069 = IT_0053*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0072;
    const ccomplex_t IT_0074 = m_b*conj(N_d4)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0075 = IT_0018*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0078 = IT_0017*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0081 = IT_0011*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + -IT_0079 +
       0.333333333333333*IT_0082);
    const ccomplex_t IT_0084 = 3*IT_0083;
    const ccomplex_t IT_0085 = 0.166666666666667*IT_0084;
    const ccomplex_t IT_0086 = IT_0073*IT_0085;
    const ccomplex_t IT_0087 = IT_0071*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = 0.166666666666667*IT_0008 + -IT_0034 + -IT_0052
       + -IT_0070 + -IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0008;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0015*IT_0017*IT_0019;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0094 = IT_0011*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = m_b*N_d1*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0097 = IT_0018*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 1.5*IT_0098);
    const ccomplex_t IT_0100 = (-0.333333333333333)*IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = IT_0009*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0105 = IT_0011*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d2*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0108 = IT_0018*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + 1.5*IT_0109);
    const ccomplex_t IT_0111 = (-0.333333333333333)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0015*IT_0017*IT_0019;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = IT_0035*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0015*IT_0017*IT_0019;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0120 = IT_0011*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_b*N_d3*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0123 = IT_0018*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 1.5*IT_0124);
    const ccomplex_t IT_0126 = (-0.333333333333333)*IT_0125;
    const ccomplex_t IT_0127 = IT_0118*IT_0126;
    const ccomplex_t IT_0128 = IT_0053*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0015*IT_0017*IT_0019;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0133 = IT_0011*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_b*N_d4*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0136 = IT_0018*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + 1.5*IT_0137);
    const ccomplex_t IT_0139 = (-0.333333333333333)*IT_0138;
    const ccomplex_t IT_0140 = IT_0131*IT_0139;
    const ccomplex_t IT_0141 = IT_0071*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0103 + -IT_0116 + -IT_0129 + -IT_0142;
    const ccomplex_t IT_0144 = m_b*m_d*IT_0002;
    const ccomplex_t IT_0145 = s_34*IT_0002;
    const ccomplex_t IT_0146 = (-18)*IT_0145;
    const ccomplex_t IT_0147 = s_13*s_14;
    const ccomplex_t IT_0148 = 36*IT_0147;
    const ccomplex_t IT_0149 = IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = m_b*s_13;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0152 = IT_0013*IT_0100;
    const ccomplex_t IT_0153 = IT_0009*IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0155 = IT_0073*IT_0139;
    const ccomplex_t IT_0156 = IT_0071*IT_0154*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0158 = IT_0004*IT_0157;
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = m_sG*IT_0160;
    const ccomplex_t IT_0162 = m_d*IT_0103;
    const ccomplex_t IT_0163 = m_d*IT_0116;
    const ccomplex_t IT_0164 = m_d*IT_0129;
    const ccomplex_t IT_0165 = m_d*IT_0142;
    const ccomplex_t IT_0166 = IT_0153 + IT_0156 + (-0.166666666666667)
      *IT_0161 + -IT_0162 + -IT_0163 + -IT_0164 + -IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0168 = IT_0037*IT_0111;
    const ccomplex_t IT_0169 = IT_0035*IT_0167*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0171 = IT_0055*IT_0126;
    const ccomplex_t IT_0172 = IT_0053*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = -IT_0169 + -IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0161;
    const ccomplex_t IT_0175 = -conj(IT_0173);
    const ccomplex_t IT_0176 = conj(IT_0166) + 3*conj(IT_0174) + IT_0175;
    const ccomplex_t IT_0177 = m_d*IT_0008;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = IT_0067*IT_0118;
    const ccomplex_t IT_0180 = IT_0053*IT_0170*IT_0179;
    const ccomplex_t IT_0181 = m_d*IT_0034;
    const ccomplex_t IT_0182 = m_d*IT_0052;
    const ccomplex_t IT_0183 = m_d*IT_0070;
    const ccomplex_t IT_0184 = m_d*IT_0088;
    const ccomplex_t IT_0185 = 0.166666666666667*IT_0177 + IT_0180 + -IT_0181 
      + -IT_0182 + -IT_0183 + -IT_0184;
    const ccomplex_t IT_0186 = IT_0085*IT_0131;
    const ccomplex_t IT_0187 = IT_0071*IT_0154*IT_0186;
    const ccomplex_t IT_0188 = IT_0031*IT_0092;
    const ccomplex_t IT_0189 = IT_0009*IT_0151*IT_0188;
    const ccomplex_t IT_0190 = IT_0049*IT_0113;
    const ccomplex_t IT_0191 = IT_0035*IT_0167*IT_0190;
    const ccomplex_t IT_0192 = -IT_0187 + -IT_0189 + -IT_0191;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0185 + IT_0193;
    const ccomplex_t IT_0195 = -conj(IT_0192);
    const ccomplex_t IT_0196 = conj(IT_0185) + IT_0195;
    const ccomplex_t IT_0197 = 6*IT_0174;
    const ccomplex_t IT_0198 = 18*IT_0166;
    const ccomplex_t IT_0199 = (-18)*IT_0173;
    const ccomplex_t IT_0200 = IT_0197 + IT_0198 + IT_0199;
    const ccomplex_t IT_0201 = 6*conj(IT_0174);
    const ccomplex_t IT_0202 = 18*conj(IT_0166);
    const ccomplex_t IT_0203 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0204 = IT_0201 + IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = conj(IT_0090)*(IT_0166 + -IT_0173 + 3*IT_0174) 
      + IT_0090*IT_0176 + conj(IT_0143)*(IT_0178 + 3*IT_0194) + IT_0143*(conj
      (IT_0178) + 3*IT_0196) + 0.166666666666667*conj(IT_0089)*IT_0200 +
       0.166666666666667*IT_0089*IT_0204;
    const ccomplex_t IT_0206 = m_b*m_d;
    const ccomplex_t IT_0207 = 18*IT_0173;
    const ccomplex_t IT_0208 = (-18)*IT_0166;
    const ccomplex_t IT_0209 = (-6)*IT_0174;
    const ccomplex_t IT_0210 = IT_0207 + IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = 18*conj(IT_0173);
    const ccomplex_t IT_0212 = (-18)*conj(IT_0166);
    const ccomplex_t IT_0213 = (-6)*conj(IT_0174);
    const ccomplex_t IT_0214 = IT_0211 + IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = (-6)*IT_0206*(IT_0176*IT_0178 + (IT_0166 + 
      -IT_0173 + 3*IT_0174)*conj(IT_0178) + (-0.166666666666667)*conj(IT_0192)
      *IT_0200 + (-0.166666666666667)*IT_0192*IT_0204 + (-0.166666666666667)
      *conj(IT_0185)*IT_0210 + (-0.166666666666667)*IT_0185*IT_0214);
    const ccomplex_t IT_0216 = m_d*s_14;
    const ccomplex_t IT_0217 = (-6)*conj(IT_0089)*(IT_0178 + 3*IT_0185 + (-3)
      *IT_0192) + (-6)*IT_0089*(conj(IT_0178) + 3*conj(IT_0185) + (-3)*conj
      (IT_0192)) + (-18)*conj(IT_0090)*(IT_0178 + 0.333333333333333*IT_0194) + (
      -18)*IT_0090*(conj(IT_0178) + 0.333333333333333*IT_0196) + conj(IT_0143)
      *IT_0210 + IT_0143*IT_0214;
    const ccomplex_t IT_0218 = 6*s_34;
    const ccomplex_t IT_0219 = 18*s_34;
    const ccomplex_t IT_0220 = (-18)*s_34;
    const ccomplex_t IT_0221 = (-6)*s_34;
    const ccomplex_t IT_0222 = 12*IT_0147;
    const ccomplex_t IT_0223 = (-6)*IT_0145;
    const ccomplex_t IT_0224 = (((-18)*conj(IT_0089) + (-6)*conj(IT_0090))
      *IT_0143 + ((-18)*IT_0089 + (-6)*IT_0090)*conj(IT_0143))*IT_0144 + 
      (IT_0089*conj(IT_0089) + IT_0090*conj(IT_0090) + IT_0143*conj(IT_0143))
      *IT_0149 + 6*IT_0150*IT_0205 + IT_0215 + IT_0216*IT_0217 + (conj(IT_0166)
      *IT_0174 + IT_0166*conj(IT_0174) + conj(IT_0178)*IT_0185 + IT_0178*conj
      (IT_0185))*IT_0218 + (IT_0166*conj(IT_0166) + IT_0173*conj(IT_0173) +
       IT_0174*conj(IT_0174) + IT_0178*conj(IT_0178) + IT_0185*conj(IT_0185) +
       IT_0192*conj(IT_0192))*IT_0219 + (conj(IT_0166)*IT_0173 + IT_0166*conj
      (IT_0173) + conj(IT_0185)*IT_0192 + IT_0185*conj(IT_0192))*IT_0220 + (conj
      (IT_0173)*IT_0174 + IT_0173*conj(IT_0174) + conj(IT_0178)*IT_0192 +
       IT_0178*conj(IT_0192))*IT_0221 + (conj(IT_0089)*IT_0090 + IT_0089*conj
      (IT_0090))*(IT_0222 + IT_0223);
    return create_ccomplex_return(IT_0224);
}

