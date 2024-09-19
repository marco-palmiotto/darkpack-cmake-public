#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_st_1_to_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_st_1_to_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
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
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t IT_0000 = pow(m_st_1, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = m_t*conj(N_u1)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0010 + 3*IT_0015 + 3
      *IT_0022);
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = m_t*N_u1*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0029 = IT_0018*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + (-0.75)*IT_0030);
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0031;
    const ccomplex_t IT_0033 = IT_0024*IT_0032;
    const ccomplex_t IT_0034 = pow(m_t, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0000 + IT_0034 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = pow(m_sG, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0034 + -IT_0039 +
       reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0040*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.166666666666667*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -IT_0039 +
       reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0043*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*conj(N_u4)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0057 = IT_0018*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0052 + 3*IT_0055 + 3
      *IT_0058);
    const ccomplex_t IT_0060 = 0.166666666666667*IT_0059;
    const ccomplex_t IT_0061 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0062 = IT_0007*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = m_t*N_u4*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0065 = IT_0018*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + (-0.75)*IT_0066);
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0067;
    const ccomplex_t IT_0069 = IT_0060*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0000 + IT_0034 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0000 + IT_0034 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0074 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0075 = IT_0007*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0078 = IT_0012*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = m_t*conj(N_u2)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0081 = IT_0018*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + 3*IT_0079 + 3
      *IT_0082);
    const ccomplex_t IT_0084 = 0.166666666666667*IT_0083;
    const ccomplex_t IT_0085 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0086 = IT_0007*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = m_t*N_u2*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0089 = IT_0018*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.75)*IT_0090);
    const ccomplex_t IT_0092 = 0.666666666666667*IT_0091;
    const ccomplex_t IT_0093 = IT_0084*IT_0092;
    const ccomplex_t IT_0094 = IT_0073*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0097 = IT_0007*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0100 = IT_0012*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = m_t*conj(N_u3)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0103 = IT_0018*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0098 + 3*IT_0101 + 3
      *IT_0104);
    const ccomplex_t IT_0106 = 0.166666666666667*IT_0105;
    const ccomplex_t IT_0107 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0108 = IT_0007*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = m_t*N_u3*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0111 = IT_0018*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + (-0.75)*IT_0112);
    const ccomplex_t IT_0114 = 0.666666666666667*IT_0113;
    const ccomplex_t IT_0115 = IT_0106*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_13 + IT_0000 + IT_0034 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0049 + -IT_0072 + -IT_0095 + -IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0069*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0033*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0093*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0115*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0122 + IT_0125 + IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0045;
    const ccomplex_t IT_0134 = -IT_0095;
    const ccomplex_t IT_0135 = -IT_0037 + -IT_0118;
    const ccomplex_t IT_0136 = 0.5*IT_0049;
    const ccomplex_t IT_0137 = -IT_0072;
    const ccomplex_t IT_0138 = IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = IT_0038 + IT_0138;
    const ccomplex_t IT_0140 = (-0.166666666666667)*IT_0049;
    const ccomplex_t IT_0141 = conj(IT_0136) + conj(IT_0137);
    const ccomplex_t IT_0142 = conj(IT_0038) + IT_0141;
    const ccomplex_t IT_0143 = IT_0128 + IT_0131;
    const ccomplex_t IT_0144 = IT_0122 + IT_0125 + IT_0143;
    const ccomplex_t IT_0145 = conj(IT_0122) + conj(IT_0125) + conj(IT_0143);
    const ccomplex_t IT_0146 = (conj(IT_0038) + conj(IT_0046) + conj(IT_0119))
      *IT_0132 + (IT_0038 + IT_0046 + IT_0119)*conj(IT_0132) + conj(IT_0133)*
      (IT_0119 + IT_0134 + IT_0135 + IT_0139) + 2*(IT_0046 + 0.5*IT_0119 + 0.5
      *IT_0134 + 0.5*IT_0135 + 0.5*IT_0139)*conj(IT_0140) + IT_0133*(conj
      (IT_0119) + conj(IT_0134) + conj(IT_0135) + IT_0142) + 2*IT_0140*(conj
      (IT_0046) + 0.5*conj(IT_0119) + 0.5*conj(IT_0134) + 0.5*conj(IT_0135) +
       0.5*IT_0142) + (conj(IT_0046) + conj(IT_0134) + conj(IT_0135) + IT_0141)
      *IT_0144 + (IT_0046 + IT_0134 + IT_0135 + IT_0138)*IT_0145;
    const ccomplex_t IT_0147 = (-18)*IT_0001;
    const ccomplex_t IT_0148 = 36*IT_0003;
    const ccomplex_t IT_0149 = IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = conj(IT_0134) + conj(IT_0135) + IT_0141;
    const ccomplex_t IT_0151 = IT_0038 + IT_0119;
    const ccomplex_t IT_0152 = conj(IT_0038) + conj(IT_0119);
    const ccomplex_t IT_0153 = IT_0132*conj(IT_0132) + conj(IT_0046)*(IT_0119 
      + IT_0134 + IT_0135 + IT_0139) + IT_0046*(conj(IT_0119) + conj(IT_0134) +
       conj(IT_0135) + IT_0142) + conj(IT_0133)*(IT_0122 + IT_0125 + IT_0132 +
       IT_0143) + 2*conj(IT_0140)*(0.5*IT_0122 + 0.5*IT_0125 + 0.5*IT_0132 +
       IT_0133 + 0.5*IT_0143) + IT_0133*(conj(IT_0122) + conj(IT_0125) + conj
      (IT_0132) + conj(IT_0143)) + 2*IT_0140*(0.5*conj(IT_0122) + 0.5*conj
      (IT_0125) + 0.5*conj(IT_0132) + conj(IT_0133) + conj(IT_0140) + 0.5*conj
      (IT_0143)) + IT_0144*IT_0145 + (IT_0134 + IT_0135 + IT_0138)*IT_0150 +
       IT_0151*IT_0152;
    const ccomplex_t IT_0154 = m_t*IT_0045;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = cpow(IT_0042, 2);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_sG*IT_0156;
    const ccomplex_t IT_0158 = IT_0040*IT_0157;
    const ccomplex_t IT_0159 = IT_0047*IT_0157;
    const ccomplex_t IT_0160 = 0.166666666666667*IT_0154 + (-0.166666666666667
      )*IT_0158 + 0.5*IT_0159;
    const ccomplex_t IT_0161 = 6*s_34;
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = conj(IT_0160)*IT_0161;
    const ccomplex_t IT_0164 = 0.166666666666667*IT_0154;
    const ccomplex_t IT_0165 = m_t*IT_0049;
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = m_t*IT_0118;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = m_t*IT_0072;
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = cpow(IT_0024, 2);
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*m_N_1*IT_0171;
    const ccomplex_t IT_0173 = IT_0035*IT_0172;
    const ccomplex_t IT_0174 = cpow(IT_0041, 2);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*m_sG*IT_0174;
    const ccomplex_t IT_0176 = IT_0040*IT_0175;
    const ccomplex_t IT_0177 = IT_0047*IT_0175;
    const ccomplex_t IT_0178 = -IT_0169 + IT_0173 + (-0.166666666666667)
      *IT_0176 + 0.5*IT_0177;
    const ccomplex_t IT_0179 = cpow(IT_0060, 2);
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_4*IT_0179;
    const ccomplex_t IT_0181 = IT_0070*IT_0180;
    const ccomplex_t IT_0182 = cpow(IT_0084, 2);
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*m_N_2*IT_0182;
    const ccomplex_t IT_0184 = IT_0073*IT_0183;
    const ccomplex_t IT_0185 = cpow(IT_0106, 2);
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_N_3*IT_0185;
    const ccomplex_t IT_0187 = IT_0116*IT_0186;
    const ccomplex_t IT_0188 = IT_0181 + IT_0184 + IT_0187;
    const ccomplex_t IT_0189 = IT_0037 + IT_0095 + IT_0118;
    const ccomplex_t IT_0190 = m_t*IT_0189;
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = 2*IT_0037;
    const ccomplex_t IT_0193 = IT_0049 + IT_0192;
    const ccomplex_t IT_0194 = m_t*IT_0193;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = m_t*IT_0095;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = conj(IT_0195) + conj(IT_0197);
    const ccomplex_t IT_0199 = IT_0155*IT_0161;
    const ccomplex_t IT_0200 = IT_0195 + IT_0197;
    const ccomplex_t IT_0201 = conj(IT_0155)*IT_0161;
    const ccomplex_t IT_0202 = m_t*s_24;
    const ccomplex_t IT_0203 = 6*IT_0155;
    const ccomplex_t IT_0204 = conj(IT_0046)*IT_0202*IT_0203;
    const ccomplex_t IT_0205 = 6*conj(IT_0155);
    const ccomplex_t IT_0206 = IT_0046*IT_0202*IT_0205;
    const ccomplex_t IT_0207 = 18*s_34;
    const ccomplex_t IT_0208 = IT_0160*IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0160)*IT_0207;
    const ccomplex_t IT_0210 = cpow(IT_0032, 2);
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*m_N_1*IT_0210;
    const ccomplex_t IT_0212 = IT_0123*IT_0211;
    const ccomplex_t IT_0213 = cpow(IT_0092, 2);
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*m_N_2*IT_0213;
    const ccomplex_t IT_0215 = IT_0126*IT_0214;
    const ccomplex_t IT_0216 = IT_0212 + IT_0215;
    const ccomplex_t IT_0217 = IT_0161*conj(IT_0195);
    const ccomplex_t IT_0218 = IT_0163 + IT_0217;
    const ccomplex_t IT_0219 = conj(IT_0155)*IT_0207;
    const ccomplex_t IT_0220 = m_t*IT_0132;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = IT_0161*conj(IT_0164);
    const ccomplex_t IT_0223 = IT_0195*IT_0207;
    const ccomplex_t IT_0224 = IT_0161*conj(IT_0216);
    const ccomplex_t IT_0225 = conj(IT_0195)*IT_0207;
    const ccomplex_t IT_0226 = IT_0161*IT_0216;
    const ccomplex_t IT_0227 = IT_0197*IT_0207;
    const ccomplex_t IT_0228 = m_t*IT_0128;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = IT_0161*conj(IT_0229);
    const ccomplex_t IT_0231 = conj(IT_0197)*IT_0207;
    const ccomplex_t IT_0232 = IT_0164*IT_0207;
    const ccomplex_t IT_0233 = IT_0161*conj(IT_0221);
    const ccomplex_t IT_0234 = conj(IT_0164)*IT_0207;
    const ccomplex_t IT_0235 = IT_0161*IT_0195;
    const ccomplex_t IT_0236 = IT_0162 + IT_0235;
    const ccomplex_t IT_0237 = IT_0155*IT_0207;
    const ccomplex_t IT_0238 = IT_0207*IT_0216;
    const ccomplex_t IT_0239 = IT_0161*conj(IT_0197);
    const ccomplex_t IT_0240 = IT_0207*conj(IT_0216);
    const ccomplex_t IT_0241 = IT_0161*IT_0197;
    const ccomplex_t IT_0242 = IT_0207*IT_0229;
    const ccomplex_t IT_0243 = cpow(IT_0068, 2);
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*m_N_4*IT_0243;
    const ccomplex_t IT_0245 = IT_0120*IT_0244;
    const ccomplex_t IT_0246 = cpow(IT_0114, 2);
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*m_N_3*IT_0246;
    const ccomplex_t IT_0248 = IT_0129*IT_0247;
    const ccomplex_t IT_0249 = IT_0245 + IT_0248;
    const ccomplex_t IT_0250 = IT_0161*conj(IT_0168);
    const ccomplex_t IT_0251 = IT_0218 + IT_0239 + IT_0250;
    const ccomplex_t IT_0252 = IT_0207*conj(IT_0229);
    const ccomplex_t IT_0253 = IT_0161*IT_0229;
    const ccomplex_t IT_0254 = IT_0168*IT_0207;
    const ccomplex_t IT_0255 = IT_0161*conj(IT_0249);
    const ccomplex_t IT_0256 = conj(IT_0168)*IT_0207;
    const ccomplex_t IT_0257 = IT_0161*IT_0168;
    const ccomplex_t IT_0258 = IT_0236 + IT_0241 + IT_0257;
    const ccomplex_t IT_0259 = IT_0207*IT_0249;
    const ccomplex_t IT_0260 = m_t*IT_0122;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0161*conj(IT_0170);
    const ccomplex_t IT_0263 = IT_0207*conj(IT_0249);
    const ccomplex_t IT_0264 = IT_0161*IT_0249;
    const ccomplex_t IT_0265 = IT_0170*IT_0207;
    const ccomplex_t IT_0266 = IT_0073*IT_0214;
    const ccomplex_t IT_0267 = IT_0116*IT_0247;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = IT_0161*IT_0261;
    const ccomplex_t IT_0270 = m_t*IT_0131;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = IT_0161*IT_0271;
    const ccomplex_t IT_0273 = IT_0207*IT_0268;
    const ccomplex_t IT_0274 = IT_0161*conj(IT_0261);
    const ccomplex_t IT_0275 = IT_0161*conj(IT_0271);
    const ccomplex_t IT_0276 = conj(IT_0170)*IT_0207;
    const ccomplex_t IT_0277 = IT_0070*IT_0244;
    const ccomplex_t IT_0278 = 0.5*IT_0158 + (-0.166666666666667)*IT_0159 +
       0.166666666666667*IT_0165;
    const ccomplex_t IT_0279 = IT_0161*conj(IT_0278);
    const ccomplex_t IT_0280 = m_t*IT_0125;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = IT_0161*conj(IT_0281);
    const ccomplex_t IT_0283 = IT_0207*conj(IT_0268);
    const ccomplex_t IT_0284 = IT_0035*IT_0211;
    const ccomplex_t IT_0285 = conj(IT_0277) + conj(IT_0284);
    const ccomplex_t IT_0286 = IT_0161*IT_0170;
    const ccomplex_t IT_0287 = IT_0207*IT_0261;
    const ccomplex_t IT_0288 = IT_0207*conj(IT_0261);
    const ccomplex_t IT_0289 = IT_0207*IT_0271;
    const ccomplex_t IT_0290 = IT_0161*conj(IT_0268);
    const ccomplex_t IT_0291 = IT_0207*conj(IT_0271);
    const ccomplex_t IT_0292 = conj(IT_0278) + conj(IT_0281);
    const ccomplex_t IT_0293 = IT_0161*IT_0268;
    const ccomplex_t IT_0294 = IT_0207*IT_0281;
    const ccomplex_t IT_0295 = IT_0161*IT_0164;
    const ccomplex_t IT_0296 = IT_0207*IT_0221;
    const ccomplex_t IT_0297 = IT_0126*IT_0183;
    const ccomplex_t IT_0298 = IT_0120*IT_0180;
    const ccomplex_t IT_0299 = IT_0129*IT_0186;
    const ccomplex_t IT_0300 = IT_0297 + IT_0298 + IT_0299;
    const ccomplex_t IT_0301 = IT_0161*conj(IT_0178);
    const ccomplex_t IT_0302 = IT_0207*conj(IT_0221);
    const ccomplex_t IT_0303 = IT_0161*IT_0221;
    const ccomplex_t IT_0304 = IT_0178*IT_0207;
    const ccomplex_t IT_0305 = IT_0161*conj(IT_0300);
    const ccomplex_t IT_0306 = conj(IT_0178)*IT_0207;
    const ccomplex_t IT_0307 = IT_0161*IT_0178;
    const ccomplex_t IT_0308 = IT_0207*IT_0300;
    const ccomplex_t IT_0309 = IT_0123*IT_0172;
    const ccomplex_t IT_0310 = 0.5*IT_0176 + (-0.166666666666667)*IT_0177 +
       IT_0309;
    const ccomplex_t IT_0311 = IT_0161*conj(IT_0188);
    const ccomplex_t IT_0312 = IT_0207*conj(IT_0300);
    const ccomplex_t IT_0313 = 0.166666666666667*IT_0165;
    const ccomplex_t IT_0314 = conj(IT_0166) + conj(IT_0191);
    const ccomplex_t IT_0315 = IT_0161*IT_0300;
    const ccomplex_t IT_0316 = IT_0188*IT_0207;
    const ccomplex_t IT_0317 = IT_0161*conj(IT_0310);
    const ccomplex_t IT_0318 = conj(IT_0188)*IT_0207;
    const ccomplex_t IT_0319 = IT_0161*IT_0188;
    const ccomplex_t IT_0320 = IT_0207*IT_0310;
    const ccomplex_t IT_0321 = IT_0161*IT_0310;
    const ccomplex_t IT_0322 = IT_0191*IT_0207;
    const ccomplex_t IT_0323 = IT_0161*IT_0278;
    const ccomplex_t IT_0324 = IT_0161*IT_0281;
    const ccomplex_t IT_0325 = IT_0207*IT_0277;
    const ccomplex_t IT_0326 = 6*IT_0278;
    const ccomplex_t IT_0327 = 6*IT_0216;
    const ccomplex_t IT_0328 = 6*IT_0229;
    const ccomplex_t IT_0329 = 6*IT_0249;
    const ccomplex_t IT_0330 = 6*IT_0261;
    const ccomplex_t IT_0331 = 6*IT_0271;
    const ccomplex_t IT_0332 = 6*IT_0281;
    const ccomplex_t IT_0333 = 18*IT_0160;
    const ccomplex_t IT_0334 = 18*IT_0195;
    const ccomplex_t IT_0335 = 18*IT_0197;
    const ccomplex_t IT_0336 = 18*IT_0268;
    const ccomplex_t IT_0337 = 18*IT_0168;
    const ccomplex_t IT_0338 = 18*IT_0170;
    const ccomplex_t IT_0339 = 6*conj(IT_0278);
    const ccomplex_t IT_0340 = 6*conj(IT_0216);
    const ccomplex_t IT_0341 = 6*conj(IT_0229);
    const ccomplex_t IT_0342 = 6*conj(IT_0249);
    const ccomplex_t IT_0343 = 6*conj(IT_0261);
    const ccomplex_t IT_0344 = 6*conj(IT_0271);
    const ccomplex_t IT_0345 = 6*conj(IT_0281);
    const ccomplex_t IT_0346 = 18*conj(IT_0160);
    const ccomplex_t IT_0347 = 18*conj(IT_0195);
    const ccomplex_t IT_0348 = 18*conj(IT_0197);
    const ccomplex_t IT_0349 = 18*conj(IT_0268);
    const ccomplex_t IT_0350 = 18*conj(IT_0168);
    const ccomplex_t IT_0351 = 18*conj(IT_0170);
    const ccomplex_t IT_0352 = 6*IT_0160;
    const ccomplex_t IT_0353 = 6*IT_0195;
    const ccomplex_t IT_0354 = 6*IT_0197;
    const ccomplex_t IT_0355 = 6*IT_0268;
    const ccomplex_t IT_0356 = 6*IT_0168;
    const ccomplex_t IT_0357 = 6*IT_0170;
    const ccomplex_t IT_0358 = 18*IT_0155;
    const ccomplex_t IT_0359 = 18*IT_0278;
    const ccomplex_t IT_0360 = 18*IT_0216;
    const ccomplex_t IT_0361 = 18*IT_0229;
    const ccomplex_t IT_0362 = 18*IT_0249;
    const ccomplex_t IT_0363 = 18*IT_0261;
    const ccomplex_t IT_0364 = 18*IT_0271;
    const ccomplex_t IT_0365 = 18*IT_0281;
    const ccomplex_t IT_0366 = 6*conj(IT_0160);
    const ccomplex_t IT_0367 = 6*conj(IT_0195);
    const ccomplex_t IT_0368 = 6*conj(IT_0197);
    const ccomplex_t IT_0369 = 6*conj(IT_0268);
    const ccomplex_t IT_0370 = 6*conj(IT_0168);
    const ccomplex_t IT_0371 = 6*conj(IT_0170);
    const ccomplex_t IT_0372 = 18*conj(IT_0155);
    const ccomplex_t IT_0373 = 18*conj(IT_0278);
    const ccomplex_t IT_0374 = 18*conj(IT_0216);
    const ccomplex_t IT_0375 = 18*conj(IT_0229);
    const ccomplex_t IT_0376 = 18*conj(IT_0249);
    const ccomplex_t IT_0377 = 18*conj(IT_0261);
    const ccomplex_t IT_0378 = 18*conj(IT_0271);
    const ccomplex_t IT_0379 = 18*conj(IT_0281);
    const ccomplex_t IT_0380 = 6*IT_0164;
    const ccomplex_t IT_0381 = 6*IT_0166;
    const ccomplex_t IT_0382 = 6*IT_0178;
    const ccomplex_t IT_0383 = 6*IT_0188;
    const ccomplex_t IT_0384 = 6*IT_0191;
    const ccomplex_t IT_0385 = 18*IT_0313;
    const ccomplex_t IT_0386 = 18*IT_0221;
    const ccomplex_t IT_0387 = 18*IT_0300;
    const ccomplex_t IT_0388 = 18*IT_0310;
    const ccomplex_t IT_0389 = 36*IT_0155 + IT_0352 + IT_0353 + IT_0354 +
       IT_0355 + IT_0356 + IT_0357 + IT_0359 + IT_0360 + IT_0361 + IT_0362 +
       IT_0363 + IT_0364 + IT_0365 + IT_0380 + IT_0381 + IT_0382 + IT_0383 +
       IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388;
    const ccomplex_t IT_0390 = IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0356
       + IT_0357 + 2*IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362 + IT_0363 +
       IT_0364 + IT_0365 + IT_0380 + IT_0381 + IT_0382 + IT_0383 + IT_0384 +
       IT_0385 + IT_0386 + IT_0387 + IT_0388;
    const ccomplex_t IT_0391 = 6*conj(IT_0164);
    const ccomplex_t IT_0392 = 6*conj(IT_0166);
    const ccomplex_t IT_0393 = 6*conj(IT_0178);
    const ccomplex_t IT_0394 = 6*conj(IT_0188);
    const ccomplex_t IT_0395 = 6*conj(IT_0191);
    const ccomplex_t IT_0396 = 18*conj(IT_0313);
    const ccomplex_t IT_0397 = 18*conj(IT_0221);
    const ccomplex_t IT_0398 = 18*conj(IT_0300);
    const ccomplex_t IT_0399 = 18*conj(IT_0310);
    const ccomplex_t IT_0400 = 36*conj(IT_0155) + IT_0366 + IT_0367 + IT_0368 
      + IT_0369 + IT_0370 + IT_0371 + IT_0373 + IT_0374 + IT_0375 + IT_0376 +
       IT_0377 + IT_0378 + IT_0379 + IT_0391 + IT_0392 + IT_0393 + IT_0394 +
       IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370
       + IT_0371 + 2*IT_0372 + IT_0373 + IT_0374 + IT_0375 + IT_0376 + IT_0377 +
       IT_0378 + IT_0379 + IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0395 +
       IT_0396 + IT_0397 + IT_0398 + IT_0399;
    const ccomplex_t IT_0402 = 6*IT_0313;
    const ccomplex_t IT_0403 = 6*IT_0221;
    const ccomplex_t IT_0404 = 6*IT_0300;
    const ccomplex_t IT_0405 = 6*IT_0310;
    const ccomplex_t IT_0406 = 18*IT_0164;
    const ccomplex_t IT_0407 = 18*IT_0166;
    const ccomplex_t IT_0408 = 18*IT_0178;
    const ccomplex_t IT_0409 = 18*IT_0188;
    const ccomplex_t IT_0410 = 18*IT_0191;
    const ccomplex_t IT_0411 = IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330
       + IT_0331 + IT_0332 + IT_0333 + IT_0334 + IT_0335 + IT_0336 + IT_0337 +
       IT_0338 + IT_0402 + IT_0403 + IT_0404 + IT_0405 + IT_0406 + IT_0407 +
       IT_0408 + IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = 6*conj(IT_0313);
    const ccomplex_t IT_0413 = 6*conj(IT_0221);
    const ccomplex_t IT_0414 = 6*conj(IT_0300);
    const ccomplex_t IT_0415 = 6*conj(IT_0310);
    const ccomplex_t IT_0416 = 18*conj(IT_0164);
    const ccomplex_t IT_0417 = 18*conj(IT_0166);
    const ccomplex_t IT_0418 = 18*conj(IT_0178);
    const ccomplex_t IT_0419 = 18*conj(IT_0188);
    const ccomplex_t IT_0420 = 18*conj(IT_0191);
    const ccomplex_t IT_0421 = IT_0339 + IT_0340 + IT_0341 + IT_0342 + IT_0343
       + IT_0344 + IT_0345 + IT_0346 + IT_0347 + IT_0348 + IT_0349 + IT_0350 +
       IT_0351 + IT_0412 + IT_0413 + IT_0414 + IT_0415 + IT_0416 + IT_0417 +
       IT_0418 + IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = 6*IT_0140;
    const ccomplex_t IT_0423 = 6*IT_0133;
    const ccomplex_t IT_0424 = 6*IT_0143;
    const ccomplex_t IT_0425 = 6*IT_0122;
    const ccomplex_t IT_0426 = 6*IT_0125;
    const ccomplex_t IT_0427 = 18*IT_0136;
    const ccomplex_t IT_0428 = 18*IT_0137;
    const ccomplex_t IT_0429 = 18*IT_0046;
    const ccomplex_t IT_0430 = 18*IT_0134;
    const ccomplex_t IT_0431 = 18*IT_0135;
    const ccomplex_t IT_0432 = IT_0277 + IT_0284;
    const ccomplex_t IT_0433 = 6*conj(IT_0140);
    const ccomplex_t IT_0434 = 6*conj(IT_0133);
    const ccomplex_t IT_0435 = 6*conj(IT_0143);
    const ccomplex_t IT_0436 = 6*conj(IT_0122);
    const ccomplex_t IT_0437 = 6*conj(IT_0125);
    const ccomplex_t IT_0438 = 18*conj(IT_0136);
    const ccomplex_t IT_0439 = 18*conj(IT_0137);
    const ccomplex_t IT_0440 = 18*conj(IT_0046);
    const ccomplex_t IT_0441 = 18*conj(IT_0134);
    const ccomplex_t IT_0442 = 18*conj(IT_0135);
    const ccomplex_t IT_0443 = IT_0150*(IT_0203 + IT_0326 + IT_0327 + IT_0328 
      + IT_0329 + IT_0330 + IT_0331 + IT_0332 + IT_0333 + IT_0334 + IT_0335 +
       IT_0336 + IT_0337 + IT_0338) + (IT_0134 + IT_0135 + IT_0138)*(IT_0205 +
       IT_0339 + IT_0340 + IT_0341 + IT_0342 + IT_0343 + IT_0344 + IT_0345 +
       IT_0346 + IT_0347 + IT_0348 + IT_0349 + IT_0350 + IT_0351) + IT_0145*
      (IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 +
       IT_0359 + IT_0360 + IT_0361 + IT_0362 + IT_0363 + IT_0364 + IT_0365) +
       IT_0144*(IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 + IT_0371 +
       IT_0372 + IT_0373 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 +
       IT_0379) + conj(IT_0132)*(IT_0358 + IT_0380 + IT_0381 + IT_0382 + IT_0383
       + IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388) + conj(IT_0140)
      *IT_0389 + conj(IT_0133)*IT_0390 + IT_0132*(IT_0372 + IT_0391 + IT_0392 +
       IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0399) +
       IT_0140*IT_0400 + IT_0133*IT_0401 + IT_0152*(IT_0203 + IT_0402 + IT_0403 
      + IT_0404 + IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409 + IT_0410) +
       conj(IT_0046)*IT_0411 + IT_0151*(IT_0205 + IT_0412 + IT_0413 + IT_0414 +
       IT_0415 + IT_0416 + IT_0417 + IT_0418 + IT_0419 + IT_0420) + IT_0046
      *IT_0421 + IT_0285*(IT_0422 + IT_0423 + IT_0424 + IT_0425 + IT_0426 +
       IT_0427 + IT_0428 + IT_0429 + IT_0430 + IT_0431) + IT_0432*(IT_0433 +
       IT_0434 + IT_0435 + IT_0436 + IT_0437 + IT_0438 + IT_0439 + IT_0440 +
       IT_0441 + IT_0442);
    const ccomplex_t IT_0444 = m_t*s_23;
    const ccomplex_t IT_0445 = (-18)*IT_0155;
    const ccomplex_t IT_0446 = conj(IT_0133)*IT_0444*IT_0445;
    const ccomplex_t IT_0447 = (-18)*conj(IT_0155);
    const ccomplex_t IT_0448 = IT_0133*IT_0444*IT_0447;
    const ccomplex_t IT_0449 = (-6)*IT_0155;
    const ccomplex_t IT_0450 = conj(IT_0046)*IT_0444*IT_0449;
    const ccomplex_t IT_0451 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0452 = IT_0046*IT_0444*IT_0451;
    const ccomplex_t IT_0453 = IT_0000*IT_0034;
    const ccomplex_t IT_0454 = (-18)*IT_0140;
    const ccomplex_t IT_0455 = (-18)*conj(IT_0140);
    const ccomplex_t IT_0456 = (-6)*IT_0140;
    const ccomplex_t IT_0457 = (-6)*conj(IT_0140);
    const ccomplex_t IT_0458 = (-6)*IT_0046;
    const ccomplex_t IT_0459 = (-6)*conj(IT_0046);
    const ccomplex_t IT_0460 = (-18)*IT_0136;
    const ccomplex_t IT_0461 = (-18)*IT_0137;
    const ccomplex_t IT_0462 = (-18)*IT_0046;
    const ccomplex_t IT_0463 = (-18)*IT_0134;
    const ccomplex_t IT_0464 = (-18)*IT_0135;
    const ccomplex_t IT_0465 = (-6)*IT_0133;
    const ccomplex_t IT_0466 = (-6)*IT_0143;
    const ccomplex_t IT_0467 = (-18)*conj(IT_0136);
    const ccomplex_t IT_0468 = (-18)*conj(IT_0137);
    const ccomplex_t IT_0469 = (-18)*conj(IT_0046);
    const ccomplex_t IT_0470 = (-18)*conj(IT_0134);
    const ccomplex_t IT_0471 = (-18)*conj(IT_0135);
    const ccomplex_t IT_0472 = (-6)*conj(IT_0133);
    const ccomplex_t IT_0473 = (-6)*conj(IT_0143);
    const ccomplex_t IT_0474 = (-18)*IT_0133;
    const ccomplex_t IT_0475 = (-18)*IT_0143;
    const ccomplex_t IT_0476 = (-6)*IT_0136;
    const ccomplex_t IT_0477 = (-6)*IT_0137;
    const ccomplex_t IT_0478 = (-6)*IT_0134;
    const ccomplex_t IT_0479 = (-6)*IT_0135;
    const ccomplex_t IT_0480 = (-18)*conj(IT_0133);
    const ccomplex_t IT_0481 = (-18)*conj(IT_0143);
    const ccomplex_t IT_0482 = (-6)*conj(IT_0136);
    const ccomplex_t IT_0483 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0484 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0485 = (-6)*conj(IT_0135);
    const ccomplex_t IT_0486 = (-18)*IT_0132;
    const ccomplex_t IT_0487 = (-6)*IT_0038;
    const ccomplex_t IT_0488 = (-6)*IT_0119;
    const ccomplex_t IT_0489 = (-18)*conj(IT_0132);
    const ccomplex_t IT_0490 = (-6)*conj(IT_0038);
    const ccomplex_t IT_0491 = (-6)*conj(IT_0119);
    const ccomplex_t IT_0492 = conj(IT_0143)*IT_0454 + IT_0143*IT_0455 +
       IT_0150*IT_0456 + (IT_0134 + IT_0135 + IT_0138)*IT_0457 + conj(IT_0140)
      *IT_0458 + IT_0140*((-36)*conj(IT_0140) + IT_0459) + (conj(IT_0038) + conj
      (IT_0046) + conj(IT_0119))*(IT_0456 + IT_0460 + IT_0461 + IT_0462 +
       IT_0463 + IT_0464 + IT_0465 + IT_0466) + (IT_0038 + IT_0046 + IT_0119)*
      (IT_0457 + IT_0467 + IT_0468 + IT_0469 + IT_0470 + IT_0471 + IT_0472 +
       IT_0473) + conj(IT_0132)*(IT_0454 + IT_0458 + IT_0474 + IT_0475 + IT_0476
       + IT_0477 + IT_0478 + IT_0479) + 2*conj(IT_0133)*(IT_0454 + 0.5*IT_0458 +
       0.5*IT_0474 + 0.5*IT_0475 + 0.5*IT_0476 + 0.5*IT_0477 + 0.5*IT_0478 + 0.5
      *IT_0479) + IT_0132*(IT_0455 + IT_0459 + IT_0480 + IT_0481 + IT_0482 +
       IT_0483 + IT_0484 + IT_0485) + 2*IT_0133*(IT_0455 + 0.5*IT_0459 + 0.5
      *IT_0480 + 0.5*IT_0481 + 0.5*IT_0482 + 0.5*IT_0483 + 0.5*IT_0484 + 0.5
      *IT_0485) + (conj(IT_0122) + conj(IT_0125))*(IT_0454 + IT_0458 + IT_0474 +
       IT_0486 + IT_0487 + IT_0488) + (IT_0122 + IT_0125)*(IT_0455 + IT_0459 +
       IT_0480 + IT_0489 + IT_0490 + IT_0491);
    const ccomplex_t IT_0493 = IT_0046*conj(IT_0046)*IT_0149;
    const ccomplex_t IT_0494 = IT_0005*IT_0046*conj(IT_0133);
    const ccomplex_t IT_0495 = IT_0005*conj(IT_0046)*IT_0133;
    const ccomplex_t IT_0496 = IT_0133*conj(IT_0133)*IT_0149;
    const ccomplex_t IT_0497 = (-18)*IT_0164;
    const ccomplex_t IT_0498 = (-18)*IT_0166;
    const ccomplex_t IT_0499 = (-18)*IT_0178;
    const ccomplex_t IT_0500 = (-18)*IT_0188;
    const ccomplex_t IT_0501 = (-18)*IT_0191;
    const ccomplex_t IT_0502 = (-6)*IT_0313;
    const ccomplex_t IT_0503 = (-6)*IT_0221;
    const ccomplex_t IT_0504 = (-6)*IT_0300;
    const ccomplex_t IT_0505 = (-6)*IT_0310;
    const ccomplex_t IT_0506 = IT_0449 + IT_0497 + IT_0498 + IT_0499 + IT_0500
       + IT_0501 + IT_0502 + IT_0503 + IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = (-18)*conj(IT_0164);
    const ccomplex_t IT_0508 = (-18)*conj(IT_0166);
    const ccomplex_t IT_0509 = (-18)*conj(IT_0178);
    const ccomplex_t IT_0510 = (-18)*conj(IT_0188);
    const ccomplex_t IT_0511 = (-18)*conj(IT_0191);
    const ccomplex_t IT_0512 = (-6)*conj(IT_0313);
    const ccomplex_t IT_0513 = (-6)*conj(IT_0221);
    const ccomplex_t IT_0514 = (-6)*conj(IT_0300);
    const ccomplex_t IT_0515 = (-6)*conj(IT_0310);
    const ccomplex_t IT_0516 = IT_0451 + IT_0507 + IT_0508 + IT_0509 + IT_0510
       + IT_0511 + IT_0512 + IT_0513 + IT_0514 + IT_0515;
    const ccomplex_t IT_0517 = (-6)*IT_0160;
    const ccomplex_t IT_0518 = (-18)*IT_0278;
    const ccomplex_t IT_0519 = (-18)*IT_0216;
    const ccomplex_t IT_0520 = (-18)*IT_0229;
    const ccomplex_t IT_0521 = (-18)*IT_0249;
    const ccomplex_t IT_0522 = (-18)*IT_0261;
    const ccomplex_t IT_0523 = (-18)*IT_0271;
    const ccomplex_t IT_0524 = (-18)*IT_0281;
    const ccomplex_t IT_0525 = (-6)*IT_0195;
    const ccomplex_t IT_0526 = (-6)*IT_0197;
    const ccomplex_t IT_0527 = (-6)*IT_0268;
    const ccomplex_t IT_0528 = (-6)*IT_0168;
    const ccomplex_t IT_0529 = (-6)*IT_0170;
    const ccomplex_t IT_0530 = (-18)*conj(IT_0278);
    const ccomplex_t IT_0531 = (-18)*conj(IT_0216);
    const ccomplex_t IT_0532 = (-18)*conj(IT_0229);
    const ccomplex_t IT_0533 = (-18)*conj(IT_0249);
    const ccomplex_t IT_0534 = (-18)*conj(IT_0261);
    const ccomplex_t IT_0535 = (-18)*conj(IT_0271);
    const ccomplex_t IT_0536 = (-18)*conj(IT_0281);
    const ccomplex_t IT_0537 = (-6)*conj(IT_0160);
    const ccomplex_t IT_0538 = (-6)*conj(IT_0195);
    const ccomplex_t IT_0539 = (-6)*conj(IT_0197);
    const ccomplex_t IT_0540 = (-6)*conj(IT_0268);
    const ccomplex_t IT_0541 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0542 = (-6)*conj(IT_0170);
    const ccomplex_t IT_0543 = (-18)*IT_0160;
    const ccomplex_t IT_0544 = (-18)*IT_0195;
    const ccomplex_t IT_0545 = (-18)*IT_0197;
    const ccomplex_t IT_0546 = (-18)*IT_0268;
    const ccomplex_t IT_0547 = (-18)*IT_0168;
    const ccomplex_t IT_0548 = (-18)*IT_0170;
    const ccomplex_t IT_0549 = (-6)*IT_0278;
    const ccomplex_t IT_0550 = (-6)*IT_0216;
    const ccomplex_t IT_0551 = (-6)*IT_0229;
    const ccomplex_t IT_0552 = (-6)*IT_0249;
    const ccomplex_t IT_0553 = (-6)*IT_0261;
    const ccomplex_t IT_0554 = (-6)*IT_0271;
    const ccomplex_t IT_0555 = (-6)*IT_0281;
    const ccomplex_t IT_0556 = (-18)*conj(IT_0160);
    const ccomplex_t IT_0557 = (-18)*conj(IT_0195);
    const ccomplex_t IT_0558 = (-18)*conj(IT_0197);
    const ccomplex_t IT_0559 = (-18)*conj(IT_0268);
    const ccomplex_t IT_0560 = (-18)*conj(IT_0168);
    const ccomplex_t IT_0561 = (-18)*conj(IT_0170);
    const ccomplex_t IT_0562 = (-6)*conj(IT_0278);
    const ccomplex_t IT_0563 = (-6)*conj(IT_0216);
    const ccomplex_t IT_0564 = (-6)*conj(IT_0229);
    const ccomplex_t IT_0565 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0566 = (-6)*conj(IT_0261);
    const ccomplex_t IT_0567 = (-6)*conj(IT_0271);
    const ccomplex_t IT_0568 = (-6)*conj(IT_0281);
    const ccomplex_t IT_0569 = (-36)*conj(IT_0155);
    const ccomplex_t IT_0570 = (conj(IT_0216) + conj(IT_0229) + conj(IT_0249) 
      + conj(IT_0261) + conj(IT_0271) + conj(IT_0278) + conj(IT_0281))*IT_0445 +
       (IT_0216 + IT_0229 + IT_0249 + IT_0261 + IT_0271 + IT_0278 + IT_0281)
      *IT_0447 + (conj(IT_0168) + conj(IT_0170) + IT_0198 + conj(IT_0268))
      *IT_0449 + (IT_0168 + IT_0170 + IT_0200 + IT_0268)*IT_0451 + IT_0285
      *IT_0506 + IT_0432*IT_0516 + conj(IT_0155)*IT_0517 + (conj(IT_0221) + conj
      (IT_0300) + conj(IT_0310) + conj(IT_0313))*(IT_0445 + IT_0517 + IT_0518 +
       IT_0519 + IT_0520 + IT_0521 + IT_0522 + IT_0523 + IT_0524 + IT_0525 +
       IT_0526 + IT_0527 + IT_0528 + IT_0529) + (IT_0221 + IT_0300 + IT_0310 +
       IT_0313)*(IT_0447 + IT_0530 + IT_0531 + IT_0532 + IT_0533 + IT_0534 +
       IT_0535 + IT_0536 + IT_0537 + IT_0538 + IT_0539 + IT_0540 + IT_0541 +
       IT_0542) + (conj(IT_0164) + conj(IT_0166) + conj(IT_0178) + conj(IT_0188)
       + conj(IT_0191))*(IT_0449 + IT_0543 + IT_0544 + IT_0545 + IT_0546 +
       IT_0547 + IT_0548 + IT_0549 + IT_0550 + IT_0551 + IT_0552 + IT_0553 +
       IT_0554 + IT_0555) + (IT_0164 + IT_0166 + IT_0178 + IT_0188 + IT_0191)*
      (IT_0451 + IT_0556 + IT_0557 + IT_0558 + IT_0559 + IT_0560 + IT_0561 +
       IT_0562 + IT_0563 + IT_0564 + IT_0565 + IT_0566 + IT_0567 + IT_0568) +
       IT_0155*(IT_0537 + IT_0569);
    const ccomplex_t IT_0571 = IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501
       + IT_0502 + IT_0503 + IT_0504 + IT_0505 + IT_0543 + IT_0544 + IT_0545 +
       IT_0546 + IT_0547 + IT_0548 + IT_0549 + IT_0550 + IT_0551 + IT_0552 +
       IT_0553 + IT_0554 + IT_0555;
    const ccomplex_t IT_0572 = IT_0507 + IT_0508 + IT_0509 + IT_0510 + IT_0511
       + IT_0512 + IT_0513 + IT_0514 + IT_0515 + IT_0556 + IT_0557 + IT_0558 +
       IT_0559 + IT_0560 + IT_0561 + IT_0562 + IT_0563 + IT_0564 + IT_0565 +
       IT_0566 + IT_0567 + IT_0568;
    const ccomplex_t IT_0573 = (-18)*IT_0313;
    const ccomplex_t IT_0574 = (-18)*IT_0221;
    const ccomplex_t IT_0575 = (-18)*IT_0300;
    const ccomplex_t IT_0576 = (-18)*IT_0310;
    const ccomplex_t IT_0577 = (-6)*IT_0164;
    const ccomplex_t IT_0578 = (-6)*IT_0166;
    const ccomplex_t IT_0579 = (-6)*IT_0178;
    const ccomplex_t IT_0580 = (-6)*IT_0188;
    const ccomplex_t IT_0581 = (-6)*IT_0191;
    const ccomplex_t IT_0582 = IT_0517 + IT_0518 + IT_0519 + IT_0520 + IT_0521
       + IT_0522 + IT_0523 + IT_0524 + IT_0525 + IT_0526 + IT_0527 + IT_0528 +
       IT_0529 + IT_0573 + IT_0574 + IT_0575 + IT_0576 + IT_0577 + IT_0578 +
       IT_0579 + IT_0580 + IT_0581;
    const ccomplex_t IT_0583 = (-36)*IT_0155 + IT_0517 + IT_0518 + IT_0519 +
       IT_0520 + IT_0521 + IT_0522 + IT_0523 + IT_0524 + IT_0525 + IT_0526 +
       IT_0527 + IT_0528 + IT_0529 + IT_0573 + IT_0574 + IT_0575 + IT_0576 +
       IT_0577 + IT_0578 + IT_0579 + IT_0580 + IT_0581;
    const ccomplex_t IT_0584 = (-18)*conj(IT_0313);
    const ccomplex_t IT_0585 = (-18)*conj(IT_0221);
    const ccomplex_t IT_0586 = (-18)*conj(IT_0300);
    const ccomplex_t IT_0587 = (-18)*conj(IT_0310);
    const ccomplex_t IT_0588 = (-6)*conj(IT_0164);
    const ccomplex_t IT_0589 = (-6)*conj(IT_0166);
    const ccomplex_t IT_0590 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0591 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0592 = (-6)*conj(IT_0191);
    const ccomplex_t IT_0593 = IT_0530 + IT_0531 + IT_0532 + IT_0533 + IT_0534
       + IT_0535 + IT_0536 + IT_0537 + IT_0538 + IT_0539 + IT_0540 + IT_0541 +
       IT_0542 + IT_0584 + IT_0585 + IT_0586 + IT_0587 + IT_0588 + IT_0589 +
       IT_0590 + IT_0591 + IT_0592;
    const ccomplex_t IT_0594 = (-18)*IT_0038;
    const ccomplex_t IT_0595 = (-18)*IT_0119;
    const ccomplex_t IT_0596 = (-6)*IT_0132;
    const ccomplex_t IT_0597 = (-18)*conj(IT_0038);
    const ccomplex_t IT_0598 = (-18)*conj(IT_0119);
    const ccomplex_t IT_0599 = (-6)*conj(IT_0132);
    const ccomplex_t IT_0600 = IT_0530 + IT_0531 + IT_0532 + IT_0533 + IT_0534
       + IT_0535 + IT_0536 + IT_0537 + IT_0538 + IT_0539 + IT_0540 + IT_0541 +
       IT_0542 + IT_0569 + IT_0584 + IT_0585 + IT_0586 + IT_0587 + IT_0588 +
       IT_0589 + IT_0590 + IT_0591 + IT_0592;
    const ccomplex_t IT_0601 = IT_0150*IT_0506 + (IT_0134 + IT_0135 + IT_0138)
      *IT_0516 + conj(IT_0132)*(IT_0445 + IT_0517 + IT_0518 + IT_0519 + IT_0520 
      + IT_0521 + IT_0522 + IT_0523 + IT_0524 + IT_0525 + IT_0526 + IT_0527 +
       IT_0528 + IT_0529) + IT_0132*(IT_0447 + IT_0530 + IT_0531 + IT_0532 +
       IT_0533 + IT_0534 + IT_0535 + IT_0536 + IT_0537 + IT_0538 + IT_0539 +
       IT_0540 + IT_0541 + IT_0542) + IT_0152*(IT_0449 + IT_0543 + IT_0544 +
       IT_0545 + IT_0546 + IT_0547 + IT_0548 + IT_0549 + IT_0550 + IT_0551 +
       IT_0552 + IT_0553 + IT_0554 + IT_0555) + IT_0151*(IT_0451 + IT_0556 +
       IT_0557 + IT_0558 + IT_0559 + IT_0560 + IT_0561 + IT_0562 + IT_0563 +
       IT_0564 + IT_0565 + IT_0566 + IT_0567 + IT_0568) + conj(IT_0046)*IT_0571 
      + IT_0046*IT_0572 + IT_0145*(IT_0445 + IT_0573 + IT_0574 + IT_0575 +
       IT_0576 + IT_0577 + IT_0578 + IT_0579 + IT_0580 + IT_0581) + conj(IT_0133
      )*IT_0582 + conj(IT_0140)*IT_0583 + IT_0144*(IT_0447 + IT_0584 + IT_0585 +
       IT_0586 + IT_0587 + IT_0588 + IT_0589 + IT_0590 + IT_0591 + IT_0592) +
       IT_0133*IT_0593 + IT_0285*(IT_0456 + IT_0462 + IT_0465 + IT_0594 +
       IT_0595 + IT_0596) + IT_0432*(IT_0457 + IT_0469 + IT_0472 + IT_0597 +
       IT_0598 + IT_0599) + IT_0140*IT_0600;
    const ccomplex_t IT_0602 = IT_0005*IT_0146 + IT_0149*IT_0153 + conj
      (IT_0155)*IT_0162 + IT_0155*(36*s_34*conj(IT_0155) + IT_0163) + (conj
      (IT_0164) + conj(IT_0166) + conj(IT_0168) + conj(IT_0170) + conj(IT_0178) 
      + conj(IT_0188) + conj(IT_0191) + IT_0198)*IT_0199 + (IT_0164 + IT_0168 +
       IT_0170 + IT_0178 + IT_0188 + IT_0191 + IT_0200)*IT_0201 + 2*IT_0204 + 2
      *IT_0206 + (conj(IT_0160) + conj(IT_0195))*IT_0208 + IT_0195*IT_0209 +
       IT_0216*(IT_0218 + IT_0219) + IT_0221*(IT_0219 + IT_0222) + conj(IT_0195)
      *IT_0223 + IT_0197*(IT_0209 + IT_0224 + IT_0225) + conj(IT_0197)*(IT_0208 
      + IT_0223 + IT_0226 + IT_0227) + IT_0168*(IT_0209 + IT_0224 + IT_0225 +
       IT_0230 + IT_0231) + conj(IT_0164)*IT_0232 + IT_0178*(IT_0233 + IT_0234) 
      + conj(IT_0216)*(IT_0236 + IT_0237 + IT_0238) + IT_0229*(IT_0218 + IT_0219
       + IT_0239 + IT_0240) + conj(IT_0229)*(IT_0236 + IT_0237 + IT_0238 +
       IT_0241 + IT_0242) + IT_0249*(IT_0219 + IT_0240 + IT_0251 + IT_0252) +
       conj(IT_0168)*(IT_0208 + IT_0223 + IT_0226 + IT_0227 + IT_0253 + IT_0254)
       + IT_0170*(IT_0209 + IT_0224 + IT_0225 + IT_0230 + IT_0231 + IT_0255 +
       IT_0256) + conj(IT_0249)*(IT_0237 + IT_0238 + IT_0242 + IT_0258 + IT_0259
      ) + IT_0261*(IT_0219 + IT_0240 + IT_0251 + IT_0252 + IT_0262 + IT_0263) +
       conj(IT_0170)*(IT_0208 + IT_0223 + IT_0226 + IT_0227 + IT_0253 + IT_0254 
      + IT_0264 + IT_0265) + conj(IT_0268)*(IT_0199 + IT_0208 + IT_0223 +
       IT_0226 + IT_0227 + IT_0253 + IT_0254 + IT_0264 + IT_0265 + IT_0269 +
       IT_0272 + IT_0273) + IT_0268*(IT_0201 + IT_0209 + IT_0224 + IT_0225 +
       IT_0230 + IT_0231 + IT_0255 + IT_0256 + IT_0274 + IT_0275 + IT_0276) +
       IT_0277*(IT_0201 + IT_0209 + IT_0224 + IT_0225 + IT_0230 + IT_0231 +
       IT_0255 + IT_0256 + IT_0274 + IT_0275 + IT_0276 + IT_0279 + IT_0282 +
       IT_0283) + IT_0284*(IT_0201 + IT_0209 + IT_0224 + IT_0225 + IT_0230 +
       IT_0231 + IT_0255 + IT_0256 + IT_0274 + IT_0275 + IT_0276 + IT_0279 +
       IT_0282 + IT_0283 + IT_0207*IT_0285) + conj(IT_0261)*(IT_0237 + IT_0238 +
       IT_0242 + IT_0258 + IT_0259 + IT_0286 + IT_0287) + IT_0271*(IT_0219 +
       IT_0240 + IT_0251 + IT_0252 + IT_0262 + IT_0263 + IT_0288) + conj(IT_0271
      )*(IT_0237 + IT_0238 + IT_0242 + IT_0258 + IT_0259 + IT_0286 + IT_0287 +
       IT_0289) + IT_0281*(IT_0218 + IT_0219 + IT_0239 + IT_0240 + IT_0250 +
       IT_0252 + IT_0262 + IT_0263 + IT_0288 + IT_0290 + IT_0291) + IT_0278*
      (IT_0218 + IT_0219 + IT_0239 + IT_0240 + IT_0250 + IT_0252 + IT_0262 +
       IT_0263 + IT_0288 + IT_0290 + IT_0291 + IT_0207*IT_0292) + IT_0292*
      (IT_0236 + IT_0237 + IT_0238 + IT_0241 + IT_0242 + IT_0257 + IT_0259 +
       IT_0286 + IT_0287 + IT_0289 + IT_0293 + IT_0294) + conj(IT_0221)*(IT_0237
       + IT_0295 + IT_0296) + IT_0300*(IT_0219 + IT_0222 + IT_0301 + IT_0302) +
       conj(IT_0178)*(IT_0232 + IT_0303 + IT_0304) + IT_0188*(IT_0233 + IT_0234 
      + IT_0305 + IT_0306) + conj(IT_0300)*(IT_0237 + IT_0295 + IT_0296 +
       IT_0307 + IT_0308) + IT_0310*(IT_0219 + IT_0222 + IT_0301 + IT_0302 +
       IT_0311 + IT_0312) + IT_0313*(IT_0219 + IT_0222 + IT_0301 + IT_0302 +
       IT_0311 + IT_0312 + IT_0207*(conj(IT_0310) + conj(IT_0313)) + IT_0161
      *IT_0314) + conj(IT_0188)*(IT_0232 + IT_0303 + IT_0304 + IT_0315 + IT_0316
      ) + IT_0191*(IT_0233 + IT_0234 + IT_0305 + IT_0306 + IT_0317 + IT_0318) +
       IT_0166*(IT_0201 + IT_0233 + IT_0234 + IT_0305 + IT_0306 + IT_0161*conj
      (IT_0313) + IT_0207*IT_0314 + IT_0317 + IT_0318) + conj(IT_0310)*(IT_0237 
      + IT_0295 + IT_0296 + IT_0307 + IT_0308 + IT_0319 + IT_0320) + conj
      (IT_0313)*(IT_0161*IT_0191 + IT_0237 + IT_0295 + IT_0296 + IT_0307 +
       IT_0308 + IT_0319 + IT_0320) + IT_0314*(IT_0232 + IT_0303 + IT_0304 +
       IT_0315 + IT_0316 + IT_0321 + IT_0322) + IT_0285*(IT_0199 + IT_0208 +
       IT_0223 + IT_0226 + IT_0227 + IT_0253 + IT_0254 + IT_0264 + IT_0265 +
       IT_0269 + IT_0272 + IT_0273 + IT_0323 + IT_0324 + IT_0325) + IT_0202
      *IT_0443 + 2*IT_0446 + 2*IT_0448 + 2*IT_0450 + 2*IT_0452 + IT_0453*IT_0492
       + 2*IT_0493 + 2*IT_0494 + 2*IT_0495 + 2*IT_0496 + IT_0034*IT_0570 +
       IT_0444*IT_0601;
    return create_ccomplex_return(IT_0602);
}

