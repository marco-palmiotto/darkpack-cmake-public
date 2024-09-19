#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sb_1_to_d_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sb_1_to_d_anti_b(
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
    const ccomplex_t IT_0001 = pow(m_sd_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0007;
    const ccomplex_t IT_0011 = m_b*N_d1*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0013 + 0.333333333333333
      *IT_0018 + -IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (-3)*IT_0023;
    const ccomplex_t IT_0025 = 0.166666666666667*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = m_b*N_d2*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0039 = IT_0015*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0034 + -IT_0037 + (
      -0.333333333333333)*IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = (-3)*IT_0043;
    const ccomplex_t IT_0045 = 0.166666666666667*IT_0044;
    const ccomplex_t IT_0046 = IT_0031*IT_0045;
    const ccomplex_t IT_0047 = IT_0029*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = m_b*N_d3*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0053 = IT_0010*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0056 = IT_0015*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0059 = IT_0007*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0054 + 0.333333333333333
      *IT_0057 + -IT_0060);
    const ccomplex_t IT_0062 = 3*IT_0061;
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = IT_0051*IT_0063;
    const ccomplex_t IT_0065 = IT_0049*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0071 = IT_0015*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0074 = IT_0007*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = m_b*N_d4*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0072 + (-3)*IT_0075 + 3
      *IT_0078);
    const ccomplex_t IT_0080 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = IT_0067*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0028 + -IT_0048 + -IT_0066 + -IT_0083;
    const ccomplex_t IT_0085 = pow(m_sG, 2);
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0085 +
       reg_prop, -1);
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = IT_0086*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0015;
    const ccomplex_t IT_0094 = (-0.333333333333333)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0096 = IT_0015*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = m_b*conj(N_d1)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0099 = IT_0010*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + 1.5*IT_0100);
    const ccomplex_t IT_0102 = (-0.333333333333333)*IT_0101;
    const ccomplex_t IT_0103 = IT_0094*IT_0102;
    const ccomplex_t IT_0104 = IT_0002*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0015;
    const ccomplex_t IT_0107 = (-0.333333333333333)*IT_0106;
    const ccomplex_t IT_0108 = m_b*conj(N_d2)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0109 = IT_0010*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0112 = IT_0015*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + 0.666666666666667
      *IT_0113);
    const ccomplex_t IT_0115 = 1.5*IT_0114;
    const ccomplex_t IT_0116 = (-0.333333333333333)*IT_0115;
    const ccomplex_t IT_0117 = IT_0107*IT_0116;
    const ccomplex_t IT_0118 = IT_0029*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0015;
    const ccomplex_t IT_0121 = (-0.333333333333333)*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0123 = IT_0015*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = m_b*conj(N_d3)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0126 = IT_0010*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + 1.5*IT_0127);
    const ccomplex_t IT_0129 = (-0.333333333333333)*IT_0128;
    const ccomplex_t IT_0130 = IT_0121*IT_0129;
    const ccomplex_t IT_0131 = IT_0049*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0015;
    const ccomplex_t IT_0134 = (-0.333333333333333)*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(N_d4)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0136 = IT_0010*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0139 = IT_0015*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + 0.666666666666667
      *IT_0140);
    const ccomplex_t IT_0142 = 1.5*IT_0141;
    const ccomplex_t IT_0143 = (-0.333333333333333)*IT_0142;
    const ccomplex_t IT_0144 = IT_0134*IT_0143;
    const ccomplex_t IT_0145 = IT_0067*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0091 + -IT_0105 + -IT_0119
       + -IT_0132 + -IT_0146;
    const ccomplex_t IT_0148 = m_b*m_d*IT_0001;
    const ccomplex_t IT_0149 = s_34*IT_0001;
    const ccomplex_t IT_0150 = (-18)*IT_0149;
    const ccomplex_t IT_0151 = s_13*s_14;
    const ccomplex_t IT_0152 = 36*IT_0151;
    const ccomplex_t IT_0153 = IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = m_d*s_14;
    const ccomplex_t IT_0155 = m_d*IT_0091;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0158 = IT_0051*IT_0129;
    const ccomplex_t IT_0159 = IT_0049*IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0161 = IT_0009*IT_0102;
    const ccomplex_t IT_0162 = IT_0002*IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(1.5*IT_0110 + IT_0113);
    const ccomplex_t IT_0165 = (-0.333333333333333)*IT_0164;
    const ccomplex_t IT_0166 = IT_0031*IT_0165;
    const ccomplex_t IT_0167 = IT_0029*IT_0163*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(1.5*IT_0137 + IT_0140);
    const ccomplex_t IT_0170 = (-0.333333333333333)*IT_0169;
    const ccomplex_t IT_0171 = IT_0069*IT_0170;
    const ccomplex_t IT_0172 = IT_0067*IT_0168*IT_0171;
    const ccomplex_t IT_0173 = m_d*IT_0105;
    const ccomplex_t IT_0174 = m_d*IT_0119;
    const ccomplex_t IT_0175 = m_d*IT_0132;
    const ccomplex_t IT_0176 = m_d*IT_0146;
    const ccomplex_t IT_0177 = (-0.166666666666667)*IT_0155 + -IT_0159 + 
      -IT_0162 + -IT_0167 + -IT_0172 + IT_0173 + IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0179 = IT_0087*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0181 = IT_0086*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*3)*(IT_0013 + 1./3*IT_0018 + 
      -IT_0021);
    const ccomplex_t IT_0184 = (-0.333333333333333)*IT_0183;
    const ccomplex_t IT_0185 = (-3)*IT_0184;
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0185;
    const ccomplex_t IT_0187 = IT_0094*IT_0186;
    const ccomplex_t IT_0188 = IT_0002*IT_0160*IT_0187;
    const ccomplex_t IT_0189 = (-3)*IT_0041;
    const ccomplex_t IT_0190 = 0.166666666666667*IT_0189;
    const ccomplex_t IT_0191 = IT_0107*IT_0190;
    const ccomplex_t IT_0192 = IT_0029*IT_0163*IT_0191;
    const ccomplex_t IT_0193 = IT_0063*IT_0121;
    const ccomplex_t IT_0194 = IT_0049*IT_0157*IT_0193;
    const ccomplex_t IT_0195 = IT_0080*IT_0134;
    const ccomplex_t IT_0196 = IT_0067*IT_0168*IT_0195;
    const ccomplex_t IT_0197 = m_d*IT_0028;
    const ccomplex_t IT_0198 = m_d*IT_0048;
    const ccomplex_t IT_0199 = m_d*IT_0066;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*((-0.333333333333333)
      *IT_0072 + IT_0075 + -IT_0078);
    const ccomplex_t IT_0201 = (-3)*IT_0200;
    const ccomplex_t IT_0202 = 0.166666666666667*IT_0201;
    const ccomplex_t IT_0203 = IT_0069*IT_0202;
    const ccomplex_t IT_0204 = IT_0067*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_d*IT_0205;
    const ccomplex_t IT_0207 = 0.166666666666667*IT_0181 + -IT_0188 + -IT_0192
       + -IT_0194 + -IT_0196 + IT_0197 + IT_0198 + IT_0199 + IT_0206;
    const ccomplex_t IT_0208 = 3*IT_0177;
    const ccomplex_t IT_0209 = IT_0156 + IT_0208;
    const ccomplex_t IT_0210 = 3*conj(IT_0177);
    const ccomplex_t IT_0211 = conj(IT_0156) + IT_0210;
    const ccomplex_t IT_0212 = conj(IT_0092)*(18*IT_0156 + 6*IT_0177) +
       IT_0092*(18*conj(IT_0156) + 6*conj(IT_0177)) + 6*conj(IT_0084)*(IT_0182 +
       3*IT_0207) + 6*IT_0084*(conj(IT_0182) + 3*conj(IT_0207)) + 6*conj(IT_0147
      )*IT_0209 + 6*IT_0147*IT_0211;
    const ccomplex_t IT_0213 = 6*s_34;
    const ccomplex_t IT_0214 = 18*s_34;
    const ccomplex_t IT_0215 = m_b*s_13;
    const ccomplex_t IT_0216 = (-18)*conj(IT_0182);
    const ccomplex_t IT_0217 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = (-18)*conj(IT_0207);
    const ccomplex_t IT_0220 = (-6)*conj(IT_0182);
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = (-18)*IT_0207;
    const ccomplex_t IT_0223 = (-6)*IT_0207;
    const ccomplex_t IT_0224 = m_b*m_d;
    const ccomplex_t IT_0225 = 12*IT_0151;
    const ccomplex_t IT_0226 = (-6)*IT_0149;
    const ccomplex_t IT_0227 = (-6)*(conj(IT_0084)*(IT_0092 + 3*IT_0147) +
       IT_0084*(conj(IT_0092) + 3*conj(IT_0147)))*IT_0148 + (IT_0084*conj
      (IT_0084) + IT_0092*conj(IT_0092) + IT_0147*conj(IT_0147))*IT_0153 +
       IT_0154*IT_0212 + (conj(IT_0156)*IT_0177 + IT_0156*conj(IT_0177) + conj
      (IT_0182)*IT_0207 + IT_0182*conj(IT_0207))*IT_0213 + (IT_0156*conj(IT_0156
      ) + IT_0177*conj(IT_0177) + IT_0182*conj(IT_0182) + IT_0207*conj(IT_0207))
      *IT_0214 + IT_0215*(((-18)*conj(IT_0092) + (-6)*conj(IT_0147))*IT_0182 + (
      -6)*conj(IT_0084)*IT_0209 + (-6)*IT_0084*IT_0211 + IT_0092*IT_0218 +
       IT_0147*IT_0221 + conj(IT_0147)*IT_0222 + conj(IT_0092)*IT_0223) + (((-18
      )*conj(IT_0156) + (-6)*conj(IT_0177))*IT_0182 + IT_0156*IT_0218 + IT_0177
      *IT_0221 + conj(IT_0177)*IT_0222 + conj(IT_0156)*IT_0223)*IT_0224 + (conj
      (IT_0092)*IT_0147 + IT_0092*conj(IT_0147))*(IT_0225 + IT_0226);
    return create_ccomplex_return(IT_0227);
}

