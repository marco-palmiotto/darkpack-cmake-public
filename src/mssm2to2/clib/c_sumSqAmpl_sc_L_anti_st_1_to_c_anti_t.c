#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_st_1_to_c_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_st_1_to_c_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
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
    const creal_t m_sc_L = param->m_sc_L;
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
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sc_L, 2);
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
    const ccomplex_t IT_0021 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = IT_0015*IT_0028;
    const ccomplex_t IT_0030 = pow(m_W, -1);
    const ccomplex_t IT_0031 = m_t*N_u1*e_em*IT_0030*conj(U_st_10);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0023 + 3*IT_0026 + 3
      *IT_0033);
    const ccomplex_t IT_0035 = 0.166666666666667*IT_0034;
    const ccomplex_t IT_0036 = IT_0020*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B2)*e_em;
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*e_em;
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 3*IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0052 = IT_0015*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_t*N_u2*e_em*IT_0030*conj(U_st_10);
    const ccomplex_t IT_0055 = IT_0029*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0053 + 3
      *IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B3)*e_em;
    const ccomplex_t IT_0064 = IT_0010*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*e_em;
    const ccomplex_t IT_0067 = IT_0015*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + 3*IT_0068);
    const ccomplex_t IT_0070 = 0.166666666666667*IT_0069;
    const ccomplex_t IT_0071 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0072 = IT_0010*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0075 = IT_0015*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = m_t*N_u3*e_em*IT_0030*conj(U_st_10);
    const ccomplex_t IT_0078 = IT_0029*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0073 + 3*IT_0076 + 3
      *IT_0079);
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0082 = IT_0070*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B4)*e_em;
    const ccomplex_t IT_0087 = IT_0010*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W4)*e_em;
    const ccomplex_t IT_0090 = IT_0015*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0095 = IT_0010*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0098 = IT_0015*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = m_t*N_u4*e_em*IT_0030*conj(U_st_10);
    const ccomplex_t IT_0101 = IT_0029*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + 3*IT_0099 + 3
      *IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0008 + -IT_0039 + -IT_0062
       + -IT_0085 + -IT_0108;
    const ccomplex_t IT_0110 = (-0.5)*IT_0008;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0015*IT_0028*IT_0030;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = m_t*conj(N_u1)*e_em*IT_0030*conj(U_st_00);
    const ccomplex_t IT_0117 = IT_0029*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.75)*IT_0118);
    const ccomplex_t IT_0120 = 0.666666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0112*IT_0120;
    const ccomplex_t IT_0122 = IT_0037*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0015*IT_0028*IT_0030;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0127 = IT_0010*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_t*conj(N_u2)*e_em*IT_0030*conj(U_st_00);
    const ccomplex_t IT_0130 = IT_0029*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + (-0.75)*IT_0131);
    const ccomplex_t IT_0133 = 0.666666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0125*IT_0133;
    const ccomplex_t IT_0135 = IT_0060*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0015*IT_0028*IT_0030;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0140 = IT_0010*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = m_t*conj(N_u3)*e_em*IT_0030*conj(U_st_00);
    const ccomplex_t IT_0143 = IT_0029*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + (-0.75)*IT_0144);
    const ccomplex_t IT_0146 = 0.666666666666667*IT_0145;
    const ccomplex_t IT_0147 = IT_0138*IT_0146;
    const ccomplex_t IT_0148 = IT_0083*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0015*IT_0028*IT_0030;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0153 = IT_0010*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_t*conj(N_u4)*e_em*IT_0030*conj(U_st_00);
    const ccomplex_t IT_0156 = IT_0029*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + (-0.75)*IT_0157);
    const ccomplex_t IT_0159 = 0.666666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0151*IT_0159;
    const ccomplex_t IT_0161 = IT_0106*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = -IT_0123 + -IT_0136 + -IT_0149 + -IT_0162;
    const ccomplex_t IT_0164 = s_34*IT_0005;
    const ccomplex_t IT_0165 = (-18)*IT_0164;
    const ccomplex_t IT_0166 = s_13*s_14;
    const ccomplex_t IT_0167 = 36*IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = m_c*s_14;
    const ccomplex_t IT_0170 = IT_0035*IT_0112;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0172 = IT_0037*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = IT_0058*IT_0125;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0175 = IT_0060*IT_0173*IT_0174;
    const ccomplex_t IT_0176 = IT_0081*IT_0138;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0178 = IT_0083*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = IT_0104*IT_0151;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0181 = IT_0106*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = m_c*IT_0008;
    const ccomplex_t IT_0183 = m_c*IT_0039;
    const ccomplex_t IT_0184 = m_c*IT_0062;
    const ccomplex_t IT_0185 = m_c*IT_0085;
    const ccomplex_t IT_0186 = m_c*IT_0108;
    const ccomplex_t IT_0187 = -IT_0172 + -IT_0175 + -IT_0178 + -IT_0181 + (
      -0.166666666666667)*IT_0182 + IT_0183 + IT_0184 + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0189 = IT_0000*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0191 = IT_0006*IT_0189*IT_0190;
    const ccomplex_t IT_0192 = IT_0020*IT_0120;
    const ccomplex_t IT_0193 = IT_0037*IT_0171*IT_0192;
    const ccomplex_t IT_0194 = IT_0047*IT_0133;
    const ccomplex_t IT_0195 = IT_0060*IT_0174*IT_0194;
    const ccomplex_t IT_0196 = IT_0070*IT_0146;
    const ccomplex_t IT_0197 = IT_0083*IT_0177*IT_0196;
    const ccomplex_t IT_0198 = IT_0093*IT_0159;
    const ccomplex_t IT_0199 = IT_0106*IT_0180*IT_0198;
    const ccomplex_t IT_0200 = m_c*IT_0123;
    const ccomplex_t IT_0201 = m_c*IT_0136;
    const ccomplex_t IT_0202 = m_c*IT_0149;
    const ccomplex_t IT_0203 = m_c*IT_0162;
    const ccomplex_t IT_0204 = 0.166666666666667*IT_0191 + -IT_0193 + -IT_0195
       + -IT_0197 + -IT_0199 + IT_0200 + IT_0201 + IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = (-0.5)*IT_0191;
    const ccomplex_t IT_0206 = 0.5*IT_0182;
    const ccomplex_t IT_0207 = IT_0109 + 3*IT_0110;
    const ccomplex_t IT_0208 = conj(IT_0109) + 3*conj(IT_0110);
    const ccomplex_t IT_0209 = (18*conj(IT_0109) + 6*conj(IT_0110))*IT_0187 + 
      (18*IT_0109 + 6*IT_0110)*conj(IT_0187) + conj(IT_0163)*(18*IT_0204 + 6
      *IT_0205) + IT_0163*(18*conj(IT_0204) + 6*conj(IT_0205)) + 6*conj(IT_0206)
      *IT_0207 + 6*IT_0206*IT_0208;
    const ccomplex_t IT_0210 = m_c*m_t*IT_0005;
    const ccomplex_t IT_0211 = (-6)*IT_0110;
    const ccomplex_t IT_0212 = (-18)*IT_0109;
    const ccomplex_t IT_0213 = IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0215 = (-18)*conj(IT_0109);
    const ccomplex_t IT_0216 = IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = 6*s_34;
    const ccomplex_t IT_0218 = 18*s_34;
    const ccomplex_t IT_0219 = 12*IT_0166;
    const ccomplex_t IT_0220 = (-6)*IT_0164;
    const ccomplex_t IT_0221 = m_t*s_13;
    const ccomplex_t IT_0222 = (-18)*conj(IT_0187);
    const ccomplex_t IT_0223 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0224 = IT_0222 + IT_0223;
    const ccomplex_t IT_0225 = (-6)*IT_0206;
    const ccomplex_t IT_0226 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0227 = m_c*m_t;
    const ccomplex_t IT_0228 = (IT_0109*conj(IT_0109) + IT_0110*conj(IT_0110) 
      + IT_0163*conj(IT_0163))*IT_0168 + IT_0169*IT_0209 + IT_0210*(conj(IT_0163
      )*IT_0213 + IT_0163*IT_0216) + (conj(IT_0204)*IT_0205 + IT_0204*conj
      (IT_0205) + conj(IT_0187)*IT_0206 + IT_0187*conj(IT_0206))*IT_0217 + 
      (IT_0187*conj(IT_0187) + IT_0204*conj(IT_0204) + IT_0205*conj(IT_0205) +
       IT_0206*conj(IT_0206))*IT_0218 + (conj(IT_0109)*IT_0110 + IT_0109*conj
      (IT_0110))*(IT_0219 + IT_0220) + (-6)*IT_0221*(IT_0205*IT_0208 + -1./6
      *conj(IT_0204)*IT_0213 + -1./6*IT_0204*IT_0216 + -1./6*IT_0163*IT_0224 + 
      -1./6*conj(IT_0163)*((-18)*IT_0187 + IT_0225) + -1./6*IT_0207*IT_0226) + (
      (-18)*conj(IT_0205)*IT_0206 + (-6)*IT_0205*(conj(IT_0187) + 3*conj(IT_0206
      )) + IT_0204*IT_0224 + conj(IT_0204)*IT_0225 + IT_0187*((-18)*conj(IT_0204
      ) + IT_0226))*IT_0227;
    return create_ccomplex_return(IT_0228);
}

