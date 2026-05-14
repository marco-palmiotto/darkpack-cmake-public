#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_sb_1_to_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_sb_1_to_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
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
    const creal_t m_sb_1 = param->m_sb_1;
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
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sb_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_sG*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_b*conj(N_d2)*e_em*IT_0021*U_sb_10;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017 + 3
      *IT_0024);
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0026, 2);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_N_2*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0007 + -IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0033 = cpow(IT_0032, 2);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_sG*IT_0033;
    const ccomplex_t IT_0035 = IT_0003*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0004*IT_0032;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = m_b*IT_0039;
    const ccomplex_t IT_0041 = IT_0003*IT_0037;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = m_b*IT_0042;
    const ccomplex_t IT_0044 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = m_b*N_d2*e_em*IT_0021*U_sb_00;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 1.5*IT_0049);
    const ccomplex_t IT_0051 = (-0.333333333333333)*IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0051, 2);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_N_2*IT_0052;
    const ccomplex_t IT_0054 = IT_0029*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0056 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0057 = IT_0009*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0060 = IT_0014*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = m_b*conj(N_d1)*e_em*IT_0021*U_sb_10;
    const ccomplex_t IT_0063 = IT_0020*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + (-3)*IT_0061 + 3
      *IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d1*e_em*IT_0021*U_sb_00;
    const ccomplex_t IT_0071 = IT_0020*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 1.5*IT_0072);
    const ccomplex_t IT_0074 = (-0.333333333333333)*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0055*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_b*IT_0077;
    const ccomplex_t IT_0079 = IT_0026*IT_0051;
    const ccomplex_t IT_0080 = IT_0029*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = m_b*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0084 = IT_0009*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0087 = IT_0014*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_b*conj(N_d3)*e_em*IT_0021*U_sb_10;
    const ccomplex_t IT_0090 = IT_0020*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + (-3)*IT_0088 + 3
      *IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0095 = IT_0009*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = m_b*N_d3*e_em*IT_0021*U_sb_00;
    const ccomplex_t IT_0098 = IT_0020*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + 1.5*IT_0099);
    const ccomplex_t IT_0101 = (-0.333333333333333)*IT_0100;
    const ccomplex_t IT_0102 = IT_0093*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = m_b*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0108 = IT_0009*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0111 = IT_0014*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = m_b*conj(N_d4)*e_em*IT_0021*U_sb_10;
    const ccomplex_t IT_0114 = IT_0020*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + (-3)*IT_0112 + 3
      *IT_0115);
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0118 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0119 = IT_0009*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = m_b*N_d4*e_em*IT_0021*U_sb_00;
    const ccomplex_t IT_0122 = IT_0020*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 1.5*IT_0123);
    const ccomplex_t IT_0125 = (-0.333333333333333)*IT_0124;
    const ccomplex_t IT_0126 = IT_0117*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = m_b*IT_0129;
    const ccomplex_t IT_0131 = (-0.166666666666667)*IT_0035 + (-0.5)*IT_0040 +
       0.166666666666667*IT_0043 + IT_0054 + -IT_0078 + -IT_0082 + -IT_0106 + 
      -IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0093, 2);
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*m_N_3*IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0135 = IT_0133*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0028*IT_0136;
    const ccomplex_t IT_0138 = cpow(IT_0117, 2);
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*m_N_4*IT_0138;
    const ccomplex_t IT_0140 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0143 = cpow(IT_0066, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_N_1*IT_0143;
    const ccomplex_t IT_0145 = IT_0142*IT_0144;
    const ccomplex_t IT_0146 = (-0.5)*IT_0007 + -IT_0135 + -IT_0137 + -IT_0141
       + -IT_0145;
    const ccomplex_t IT_0147 = 6*s_34;
    const ccomplex_t IT_0148 = IT_0102*IT_0134;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = m_b*IT_0149;
    const ccomplex_t IT_0151 = IT_0126*IT_0140;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_b*IT_0152;
    const ccomplex_t IT_0154 = cpow(IT_0074, 2);
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*m_N_1*IT_0154;
    const ccomplex_t IT_0156 = IT_0142*IT_0155;
    const ccomplex_t IT_0157 = cpow(IT_0125, 2);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_N_4*IT_0157;
    const ccomplex_t IT_0159 = IT_0140*IT_0158;
    const ccomplex_t IT_0160 = IT_0075*IT_0142;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = m_b*IT_0161;
    const ccomplex_t IT_0163 = IT_0079*IT_0136;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = m_b*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0035 + 0.166666666666667*IT_0040 + (-0.5
      )*IT_0043 + -IT_0150 + -IT_0153 + IT_0156 + IT_0159 + -IT_0162 + -IT_0165;
    const ccomplex_t IT_0167 = 6*IT_0000;
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0039 + 0.5*IT_0042 + 
      -IT_0149 + -IT_0152 + -IT_0161 + -IT_0164;
    const ccomplex_t IT_0169 = m_b*s_23;
    const ccomplex_t IT_0170 = 6*IT_0169;
    const ccomplex_t IT_0171 = (-0.5)*IT_0039 + (-0.166666666666667)*IT_0042 +
       IT_0077 + IT_0081 + IT_0105 + IT_0129;
    const ccomplex_t IT_0172 = 18*IT_0169;
    const ccomplex_t IT_0173 = IT_0055*IT_0144;
    const ccomplex_t IT_0174 = IT_0006*IT_0036;
    const ccomplex_t IT_0175 = IT_0127*IT_0139;
    const ccomplex_t IT_0176 = IT_0103*IT_0133;
    const ccomplex_t IT_0177 = (-0.5)*IT_0040 + 0.166666666666667*IT_0043 + 
      -IT_0078 + -IT_0082 + -IT_0106 + -IT_0130 + IT_0173 + 0.5*IT_0174 +
       IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (-18)*s_34;
    const ccomplex_t IT_0179 = IT_0034*IT_0036;
    const ccomplex_t IT_0180 = IT_0055*IT_0155;
    const ccomplex_t IT_0181 = cpow(IT_0101, 2);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_3*IT_0181;
    const ccomplex_t IT_0183 = IT_0103*IT_0182;
    const ccomplex_t IT_0184 = IT_0127*IT_0158;
    const ccomplex_t IT_0185 = (-0.5)*IT_0179 + -IT_0180 + -IT_0183 + -IT_0184;
    const ccomplex_t IT_0186 = (-18)*IT_0000;
    const ccomplex_t IT_0187 = m_b*s_24;
    const ccomplex_t IT_0188 = (-18)*IT_0187;
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0040 + (-0.5)*IT_0043 + 
      -IT_0150 + -IT_0153 + -IT_0162 + -IT_0165 + (-0.166666666666667)*IT_0174;
    const ccomplex_t IT_0190 = (-6)*s_34;
    const ccomplex_t IT_0191 = IT_0053*IT_0136;
    const ccomplex_t IT_0192 = IT_0134*IT_0182;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0179 + -IT_0191 + -IT_0192;
    const ccomplex_t IT_0194 = (-6)*IT_0000;
    const ccomplex_t IT_0195 = (-6)*IT_0187;
    const ccomplex_t IT_0196 = 18*s_34*conj(IT_0031) + 18*IT_0000*conj(IT_0131
      ) + conj(IT_0146)*IT_0147 + conj(IT_0166)*IT_0167 + conj(IT_0168)*IT_0170 
      + conj(IT_0171)*IT_0172 + conj(IT_0177)*IT_0178 + conj(IT_0185)*IT_0186 +
       conj(IT_0171)*IT_0188 + conj(IT_0189)*IT_0190 + conj(IT_0193)*IT_0194 +
       conj(IT_0168)*IT_0195;
    const ccomplex_t IT_0197 = IT_0170*conj(IT_0171) + conj(IT_0168)*IT_0172 +
       IT_0167*conj(IT_0177) + conj(IT_0166)*IT_0178 + IT_0147*conj(IT_0185) +
       conj(IT_0146)*IT_0186 + conj(IT_0168)*IT_0188 + 18*IT_0000*conj(IT_0189) 
      + conj(IT_0131)*IT_0190 + 18*s_34*conj(IT_0193) + conj(IT_0031)*IT_0194 +
       conj(IT_0171)*IT_0195;
    const ccomplex_t IT_0198 = conj(IT_0168)*IT_0170 + conj(IT_0171)*IT_0172 +
       18*IT_0000*conj(IT_0177) + conj(IT_0131)*IT_0178 + 18*s_34*conj(IT_0185) 
      + conj(IT_0031)*IT_0186 + conj(IT_0171)*IT_0188 + IT_0167*conj(IT_0189) +
       conj(IT_0166)*IT_0190 + IT_0147*conj(IT_0193) + conj(IT_0146)*IT_0194 +
       conj(IT_0168)*IT_0195;
    const ccomplex_t IT_0199 = 18*s_34*conj(IT_0146) + conj(IT_0031)*IT_0147 +
       18*IT_0000*conj(IT_0166) + conj(IT_0131)*IT_0167 + IT_0170*conj(IT_0171) 
      + conj(IT_0168)*IT_0172 + conj(IT_0168)*IT_0188 + IT_0178*conj(IT_0189) +
       conj(IT_0177)*IT_0190 + IT_0186*conj(IT_0193) + conj(IT_0185)*IT_0194 +
       conj(IT_0171)*IT_0195;
    const ccomplex_t IT_0200 = 6*IT_0187;
    const ccomplex_t IT_0201 = 18*IT_0187;
    const ccomplex_t IT_0202 = (-18)*IT_0169;
    const ccomplex_t IT_0203 = (-6)*IT_0169;
    const ccomplex_t IT_0204 = 18*IT_0000*conj(IT_0031) + 18*s_34*conj(IT_0131
      ) + IT_0147*conj(IT_0166) + conj(IT_0146)*IT_0167 + IT_0178*conj(IT_0185) 
      + conj(IT_0177)*IT_0186 + IT_0190*conj(IT_0193) + conj(IT_0189)*IT_0194 +
       conj(IT_0168)*IT_0200 + conj(IT_0171)*IT_0201 + conj(IT_0171)*IT_0202 +
       conj(IT_0168)*IT_0203;
    const ccomplex_t IT_0205 = IT_0147*conj(IT_0177) + conj(IT_0146)*IT_0178 +
       IT_0167*conj(IT_0185) + conj(IT_0166)*IT_0186 + 18*s_34*conj(IT_0189) +
       conj(IT_0031)*IT_0190 + 18*IT_0000*conj(IT_0193) + conj(IT_0131)*IT_0194 
      + conj(IT_0171)*IT_0200 + conj(IT_0168)*IT_0201 + conj(IT_0168)*IT_0202 +
       conj(IT_0171)*IT_0203;
    const ccomplex_t IT_0206 = 18*s_34*conj(IT_0177) + conj(IT_0031)*IT_0178 +
       18*IT_0000*conj(IT_0185) + conj(IT_0131)*IT_0186 + IT_0147*conj(IT_0189) 
      + conj(IT_0146)*IT_0190 + IT_0167*conj(IT_0193) + conj(IT_0166)*IT_0194 +
       conj(IT_0168)*IT_0200 + conj(IT_0171)*IT_0201 + conj(IT_0171)*IT_0202 +
       conj(IT_0168)*IT_0203;
    const ccomplex_t IT_0207 = 18*IT_0000*conj(IT_0146) + conj(IT_0131)
      *IT_0147 + 18*s_34*conj(IT_0166) + conj(IT_0031)*IT_0167 + IT_0186*conj
      (IT_0189) + conj(IT_0185)*IT_0190 + IT_0178*conj(IT_0193) + conj(IT_0177)
      *IT_0194 + conj(IT_0171)*IT_0200 + conj(IT_0168)*IT_0201 + conj(IT_0168)
      *IT_0202 + conj(IT_0171)*IT_0203;
    const ccomplex_t IT_0208 = s_34*IT_0002;
    const ccomplex_t IT_0209 = (-6)*IT_0208;
    const ccomplex_t IT_0210 = s_23*s_24;
    const ccomplex_t IT_0211 = 12*IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = (-18)*IT_0208;
    const ccomplex_t IT_0214 = 36*IT_0210;
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = IT_0000*IT_0002;
    const ccomplex_t IT_0217 = (-18)*IT_0216;
    const ccomplex_t IT_0218 = (-6)*IT_0216;
    const ccomplex_t IT_0219 = conj(IT_0146)*IT_0170 + conj(IT_0031)*IT_0172 +
       conj(IT_0185)*IT_0188 + conj(IT_0193)*IT_0195 + conj(IT_0166)*IT_0200 +
       conj(IT_0131)*IT_0201 + conj(IT_0177)*IT_0202 + conj(IT_0189)*IT_0203 +
       conj(IT_0168)*IT_0212 + conj(IT_0171)*IT_0215 + conj(IT_0171)*IT_0217 +
       conj(IT_0168)*IT_0218;
    const ccomplex_t IT_0220 = IT_0170*conj(IT_0185) + conj(IT_0146)*IT_0188 +
       IT_0172*conj(IT_0193) + conj(IT_0031)*IT_0195 + conj(IT_0177)*IT_0200 +
       conj(IT_0189)*IT_0201 + conj(IT_0166)*IT_0202 + conj(IT_0131)*IT_0203 +
       conj(IT_0171)*IT_0212 + conj(IT_0168)*IT_0215 + conj(IT_0168)*IT_0217 +
       conj(IT_0171)*IT_0218;
    const ccomplex_t IT_0221 = IT_0172*conj(IT_0185) + conj(IT_0031)*IT_0188 +
       IT_0170*conj(IT_0193) + conj(IT_0146)*IT_0195 + conj(IT_0189)*IT_0200 +
       conj(IT_0177)*IT_0201 + conj(IT_0131)*IT_0202 + conj(IT_0166)*IT_0203 +
       conj(IT_0168)*IT_0212 + conj(IT_0171)*IT_0215 + conj(IT_0171)*IT_0217 +
       conj(IT_0168)*IT_0218;
    const ccomplex_t IT_0222 = conj(IT_0031)*IT_0170 + conj(IT_0146)*IT_0172 +
       IT_0188*conj(IT_0193) + conj(IT_0185)*IT_0195 + conj(IT_0131)*IT_0200 +
       conj(IT_0166)*IT_0201 + conj(IT_0189)*IT_0202 + conj(IT_0177)*IT_0203 +
       conj(IT_0171)*IT_0212 + conj(IT_0168)*IT_0215 + conj(IT_0168)*IT_0217 +
       conj(IT_0171)*IT_0218;
    const ccomplex_t IT_0223 = IT_0031*IT_0196 + IT_0193*IT_0197 + IT_0185
      *IT_0198 + IT_0146*IT_0199 + IT_0131*IT_0204 + IT_0189*IT_0205 + IT_0177
      *IT_0206 + IT_0166*IT_0207 + IT_0171*IT_0219 + IT_0168*IT_0220 + IT_0171
      *IT_0221 + IT_0168*IT_0222;
    return create_ccomplex_return(IT_0223);
}

