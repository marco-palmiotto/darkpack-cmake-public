#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_sb_2_to_s_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_sb_2_to_s_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_s = param->m_s;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0012;
    const ccomplex_t IT_0014 = (-0.333333333333333)*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_b*conj(N_d2)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0028 = IT_0012*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + -IT_0026 +
       0.333333333333333*IT_0029);
    const ccomplex_t IT_0031 = 3*IT_0030;
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = IT_0014*IT_0032;
    const ccomplex_t IT_0034 = IT_0010*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0012;
    const ccomplex_t IT_0038 = (-0.333333333333333)*IT_0037;
    const ccomplex_t IT_0039 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0043 = IT_0018*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = m_b*conj(N_d3)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0046 = IT_0019*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0041 + (-3)*IT_0044 + 3
      *IT_0047);
    const ccomplex_t IT_0049 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_0050 = IT_0038*IT_0049;
    const ccomplex_t IT_0051 = IT_0036*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0012;
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0054;
    const ccomplex_t IT_0056 = m_b*conj(N_d4)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0057 = IT_0019*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0060 = IT_0018*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + -IT_0061 +
       0.333333333333333*IT_0064);
    const ccomplex_t IT_0066 = 3*IT_0065;
    const ccomplex_t IT_0067 = 0.166666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0055*IT_0067;
    const ccomplex_t IT_0069 = IT_0053*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0012;
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0072;
    const ccomplex_t IT_0074 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0075 = IT_0012*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0078 = IT_0018*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = m_b*conj(N_d1)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0081 = IT_0019*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + (-3)*IT_0079 + 3
      *IT_0082);
    const ccomplex_t IT_0084 = 0.166666666666667*IT_0083;
    const ccomplex_t IT_0085 = IT_0073*IT_0084;
    const ccomplex_t IT_0086 = IT_0071*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0008 + -IT_0035 + -IT_0052
       + -IT_0070 + -IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0092 = IT_0012*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = m_b*N_d4*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0095 = IT_0019*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + 1.5*IT_0096);
    const ccomplex_t IT_0098 = (-0.333333333333333)*IT_0097;
    const ccomplex_t IT_0099 = IT_0090*IT_0098;
    const ccomplex_t IT_0100 = IT_0053*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0105 = IT_0012*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d2*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0108 = IT_0019*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + 1.5*IT_0109);
    const ccomplex_t IT_0111 = (-0.333333333333333)*IT_0110;
    const ccomplex_t IT_0112 = IT_0103*IT_0111;
    const ccomplex_t IT_0113 = IT_0010*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0118 = IT_0012*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = m_b*N_d3*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0121 = IT_0019*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + 1.5*IT_0122);
    const ccomplex_t IT_0124 = (-0.333333333333333)*IT_0123;
    const ccomplex_t IT_0125 = IT_0116*IT_0124;
    const ccomplex_t IT_0126 = IT_0036*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0131 = IT_0012*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_b*N_d1*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0134 = IT_0019*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 1.5*IT_0135);
    const ccomplex_t IT_0137 = (-0.333333333333333)*IT_0136;
    const ccomplex_t IT_0138 = IT_0129*IT_0137;
    const ccomplex_t IT_0139 = IT_0071*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = -IT_0101 + -IT_0114 + -IT_0127 + -IT_0140;
    const ccomplex_t IT_0142 = m_b*m_s*IT_0002;
    const ccomplex_t IT_0143 = s_34*IT_0002;
    const ccomplex_t IT_0144 = (-18)*IT_0143;
    const ccomplex_t IT_0145 = s_13*s_14;
    const ccomplex_t IT_0146 = 36*IT_0145;
    const ccomplex_t IT_0147 = IT_0144 + IT_0146;
    const ccomplex_t IT_0148 = m_b*m_s;
    const ccomplex_t IT_0149 = m_s*IT_0008;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = m_s*IT_0140;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0153 = IT_0038*IT_0124;
    const ccomplex_t IT_0154 = IT_0036*IT_0152*IT_0153;
    const ccomplex_t IT_0155 = m_s*IT_0101;
    const ccomplex_t IT_0156 = m_s*IT_0127;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0158 = IT_0055*IT_0098;
    const ccomplex_t IT_0159 = IT_0053*IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0161 = IT_0004*IT_0160;
    const ccomplex_t IT_0162 = IT_0003*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = m_sG*IT_0163;
    const ccomplex_t IT_0165 = m_s*IT_0114;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0167 = IT_0073*IT_0137;
    const ccomplex_t IT_0168 = IT_0071*IT_0166*IT_0167;
    const ccomplex_t IT_0169 = -IT_0151 + IT_0154 + -IT_0155 + -IT_0156 +
       IT_0159 + (-0.166666666666667)*IT_0164 + -IT_0165 + IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0171 = IT_0014*IT_0111;
    const ccomplex_t IT_0172 = IT_0010*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0164;
    const ccomplex_t IT_0175 = IT_0169 + -IT_0173 + 3*IT_0174;
    const ccomplex_t IT_0176 = 3*conj(IT_0174);
    const ccomplex_t IT_0177 = -conj(IT_0173);
    const ccomplex_t IT_0178 = conj(IT_0169) + IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = IT_0067*IT_0090;
    const ccomplex_t IT_0180 = IT_0053*IT_0157*IT_0179;
    const ccomplex_t IT_0181 = IT_0084*IT_0129;
    const ccomplex_t IT_0182 = IT_0071*IT_0166*IT_0181;
    const ccomplex_t IT_0183 = -IT_0180 + -IT_0182;
    const ccomplex_t IT_0184 = 6*IT_0174;
    const ccomplex_t IT_0185 = 18*IT_0169;
    const ccomplex_t IT_0186 = (-18)*IT_0173;
    const ccomplex_t IT_0187 = IT_0184 + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 6*conj(IT_0174);
    const ccomplex_t IT_0189 = 18*conj(IT_0169);
    const ccomplex_t IT_0190 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0191 = IT_0188 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = IT_0032*IT_0103;
    const ccomplex_t IT_0193 = IT_0010*IT_0170*IT_0192;
    const ccomplex_t IT_0194 = m_s*IT_0035;
    const ccomplex_t IT_0195 = m_s*IT_0052;
    const ccomplex_t IT_0196 = IT_0049*IT_0116;
    const ccomplex_t IT_0197 = IT_0036*IT_0152*IT_0196;
    const ccomplex_t IT_0198 = m_s*IT_0087;
    const ccomplex_t IT_0199 = m_s*IT_0070;
    const ccomplex_t IT_0200 = 0.166666666666667*IT_0149 + IT_0193 + -IT_0194 
      + -IT_0195 + IT_0197 + -IT_0198 + -IT_0199;
    const ccomplex_t IT_0201 = 18*IT_0173;
    const ccomplex_t IT_0202 = (-18)*IT_0169;
    const ccomplex_t IT_0203 = (-6)*IT_0174;
    const ccomplex_t IT_0204 = IT_0201 + IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = 18*conj(IT_0173);
    const ccomplex_t IT_0206 = (-18)*conj(IT_0169);
    const ccomplex_t IT_0207 = (-6)*conj(IT_0174);
    const ccomplex_t IT_0208 = IT_0205 + IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = 6*s_34;
    const ccomplex_t IT_0210 = 18*s_34;
    const ccomplex_t IT_0211 = (-18)*s_34;
    const ccomplex_t IT_0212 = (-6)*s_34;
    const ccomplex_t IT_0213 = 12*IT_0145;
    const ccomplex_t IT_0214 = (-6)*IT_0143;
    const ccomplex_t IT_0215 = m_s*s_14;
    const ccomplex_t IT_0216 = (-3)*conj(IT_0183);
    const ccomplex_t IT_0217 = 3*conj(IT_0200);
    const ccomplex_t IT_0218 = conj(IT_0150) + IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = 3*IT_0200;
    const ccomplex_t IT_0220 = conj(IT_0009)*((-18)*IT_0150 + 6*IT_0183 + (-6)
      *IT_0200) + IT_0009*((-18)*conj(IT_0150) + 6*conj(IT_0183) + (-6)*conj
      (IT_0200)) + conj(IT_0141)*IT_0204 + IT_0141*IT_0208 + (-6)*IT_0088
      *IT_0218 + (-6)*conj(IT_0088)*(IT_0150 + (-3)*IT_0183 + IT_0219);
    const ccomplex_t IT_0221 = m_b*s_13;
    const ccomplex_t IT_0222 = (-6)*((conj(IT_0009) + 3*conj(IT_0088))*IT_0141
       + (IT_0009 + 3*IT_0088)*conj(IT_0141))*IT_0142 + (IT_0009*conj(IT_0009) +
       IT_0088*conj(IT_0088) + IT_0141*conj(IT_0141))*IT_0147 + (-6)*IT_0148*
      (conj(IT_0150)*IT_0175 + IT_0150*IT_0178 + (-0.166666666666667)*conj
      (IT_0183)*IT_0187 + (-0.166666666666667)*IT_0183*IT_0191 + (
      -0.166666666666667)*conj(IT_0200)*IT_0204 + (-0.166666666666667)*IT_0200
      *IT_0208) + (conj(IT_0169)*IT_0174 + IT_0169*conj(IT_0174) + conj(IT_0150)
      *IT_0200 + IT_0150*conj(IT_0200))*IT_0209 + (IT_0150*conj(IT_0150) +
       IT_0169*conj(IT_0169) + IT_0173*conj(IT_0173) + IT_0174*conj(IT_0174) +
       IT_0183*conj(IT_0183) + IT_0200*conj(IT_0200))*IT_0210 + (conj(IT_0169)
      *IT_0173 + IT_0169*conj(IT_0173) + conj(IT_0183)*IT_0200 + IT_0183*conj
      (IT_0200))*IT_0211 + (conj(IT_0173)*IT_0174 + IT_0173*conj(IT_0174) + conj
      (IT_0150)*IT_0183 + IT_0150*conj(IT_0183))*IT_0212 + (conj(IT_0009)
      *IT_0088 + IT_0009*conj(IT_0088))*(IT_0213 + IT_0214) + IT_0215*IT_0220 +
       6*(conj(IT_0009)*IT_0175 + IT_0009*IT_0178 + 0.166666666666667*conj
      (IT_0088)*IT_0187 + 0.166666666666667*IT_0088*IT_0191 + IT_0141*IT_0218 +
       conj(IT_0141)*(IT_0150 + (-3)*IT_0183 + IT_0219))*IT_0221;
    return create_ccomplex_return(IT_0222);
}

