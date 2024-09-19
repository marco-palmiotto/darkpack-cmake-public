#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_st_1_to_u_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_st_1_to_u_anti_t(
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
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_su_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_B1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_W1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0014 + 3*IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = sin(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = IT_0016*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = m_t*N_u1*e_em*IT_0031*conj(U_st_10);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0024 + 3*IT_0027 + 3
      *IT_0034);
    const ccomplex_t IT_0036 = 0.166666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0021*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = conj(N_B2)*e_em;
    const ccomplex_t IT_0042 = IT_0011*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W2)*e_em;
    const ccomplex_t IT_0045 = IT_0016*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + 3*IT_0046);
    const ccomplex_t IT_0048 = 0.166666666666667*IT_0047;
    const ccomplex_t IT_0049 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0050 = IT_0011*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0053 = IT_0016*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = m_t*N_u2*e_em*IT_0031*conj(U_st_10);
    const ccomplex_t IT_0056 = IT_0030*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0051 + 3*IT_0054 + 3
      *IT_0057);
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = IT_0048*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0065 = conj(N_B3)*e_em;
    const ccomplex_t IT_0066 = IT_0011*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W3)*e_em;
    const ccomplex_t IT_0069 = IT_0016*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + 3*IT_0070);
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0074 = IT_0016*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = m_t*N_u3*e_em*IT_0031*conj(U_st_10);
    const ccomplex_t IT_0077 = IT_0030*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0075 + IT_0078 +
       0.333333333333333*IT_0081);
    const ccomplex_t IT_0083 = 3*IT_0082;
    const ccomplex_t IT_0084 = 0.166666666666667*IT_0083;
    const ccomplex_t IT_0085 = IT_0072*IT_0084;
    const ccomplex_t IT_0086 = IT_0064*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0089 = conj(N_B4)*e_em;
    const ccomplex_t IT_0090 = IT_0011*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W4)*e_em;
    const ccomplex_t IT_0093 = IT_0016*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + 3*IT_0094);
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0097 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0098 = IT_0016*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = m_t*N_u4*e_em*IT_0031*conj(U_st_10);
    const ccomplex_t IT_0101 = IT_0030*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0104 = IT_0011*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0099 + IT_0102 +
       0.333333333333333*IT_0105);
    const ccomplex_t IT_0107 = 3*IT_0106;
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = IT_0096*IT_0108;
    const ccomplex_t IT_0110 = IT_0088*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 0.166666666666667*IT_0008 + -IT_0040 + -IT_0063
       + -IT_0087 + -IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0016*IT_0029*IT_0031;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0116 = IT_0011*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = m_t*conj(N_u1)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0119 = IT_0030*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0117 + (-0.75)*IT_0120);
    const ccomplex_t IT_0122 = 0.666666666666667*IT_0121;
    const ccomplex_t IT_0123 = IT_0114*IT_0122;
    const ccomplex_t IT_0124 = IT_0038*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0016*IT_0029*IT_0031;
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0129 = IT_0011*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = m_t*conj(N_u2)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0132 = IT_0030*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + (-0.75)*IT_0133);
    const ccomplex_t IT_0135 = 0.666666666666667*IT_0134;
    const ccomplex_t IT_0136 = IT_0127*IT_0135;
    const ccomplex_t IT_0137 = IT_0061*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0016*IT_0029*IT_0031;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0142 = IT_0011*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = m_t*conj(N_u3)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0145 = IT_0030*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + (-0.75)*IT_0146);
    const ccomplex_t IT_0148 = 0.666666666666667*IT_0147;
    const ccomplex_t IT_0149 = IT_0140*IT_0148;
    const ccomplex_t IT_0150 = IT_0064*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0016*IT_0029*IT_0031;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0155 = IT_0011*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = m_t*conj(N_u4)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0158 = IT_0030*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0156 + (-0.75)*IT_0159);
    const ccomplex_t IT_0161 = 0.666666666666667*IT_0160;
    const ccomplex_t IT_0162 = IT_0153*IT_0161;
    const ccomplex_t IT_0163 = IT_0088*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0125 + -IT_0138 + -IT_0151 + -IT_0164;
    const ccomplex_t IT_0166 = s_34*IT_0005;
    const ccomplex_t IT_0167 = (-18)*IT_0166;
    const ccomplex_t IT_0168 = s_13*s_14;
    const ccomplex_t IT_0169 = 36*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = m_u*s_14;
    const ccomplex_t IT_0172 = m_u*IT_0008;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0175 = IT_0000*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0177 = IT_0006*IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0180 = IT_0021*IT_0122;
    const ccomplex_t IT_0181 = IT_0038*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0183 = IT_0096*IT_0161;
    const ccomplex_t IT_0184 = IT_0088*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0186 = IT_0048*IT_0135;
    const ccomplex_t IT_0187 = IT_0061*IT_0185*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0189 = IT_0072*IT_0148;
    const ccomplex_t IT_0190 = IT_0064*IT_0188*IT_0189;
    const ccomplex_t IT_0191 = m_u*IT_0125;
    const ccomplex_t IT_0192 = m_u*IT_0138;
    const ccomplex_t IT_0193 = m_u*IT_0151;
    const ccomplex_t IT_0194 = m_u*IT_0164;
    const ccomplex_t IT_0195 = 0.166666666666667*IT_0177 + -IT_0181 + -IT_0184
       + -IT_0187 + -IT_0190 + IT_0191 + IT_0192 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = IT_0036*IT_0114;
    const ccomplex_t IT_0197 = IT_0038*IT_0179*IT_0196;
    const ccomplex_t IT_0198 = IT_0059*IT_0127;
    const ccomplex_t IT_0199 = IT_0061*IT_0185*IT_0198;
    const ccomplex_t IT_0200 = IT_0084*IT_0140;
    const ccomplex_t IT_0201 = IT_0064*IT_0188*IT_0200;
    const ccomplex_t IT_0202 = m_u*IT_0040;
    const ccomplex_t IT_0203 = m_u*IT_0063;
    const ccomplex_t IT_0204 = m_u*IT_0087;
    const ccomplex_t IT_0205 = m_u*IT_0111;
    const ccomplex_t IT_0206 = IT_0108*IT_0153;
    const ccomplex_t IT_0207 = IT_0088*IT_0182*IT_0206;
    const ccomplex_t IT_0208 = (-0.166666666666667)*IT_0172 + -IT_0197 + 
      -IT_0199 + -IT_0201 + IT_0202 + IT_0203 + IT_0204 + IT_0205 + -IT_0207;
    const ccomplex_t IT_0209 = 3*IT_0112;
    const ccomplex_t IT_0210 = IT_0009 + IT_0209;
    const ccomplex_t IT_0211 = 3*conj(IT_0112);
    const ccomplex_t IT_0212 = conj(IT_0009) + IT_0211;
    const ccomplex_t IT_0213 = (18*conj(IT_0009) + 6*conj(IT_0112))*IT_0173 + 
      (18*IT_0009 + 6*IT_0112)*conj(IT_0173) + 6*conj(IT_0165)*(IT_0178 + 3
      *IT_0195) + 6*IT_0165*(conj(IT_0178) + 3*conj(IT_0195)) + 6*conj(IT_0208)
      *IT_0210 + 6*IT_0208*IT_0212;
    const ccomplex_t IT_0214 = m_t*m_u*IT_0005;
    const ccomplex_t IT_0215 = 6*s_34;
    const ccomplex_t IT_0216 = 18*s_34;
    const ccomplex_t IT_0217 = 12*IT_0168;
    const ccomplex_t IT_0218 = (-6)*IT_0166;
    const ccomplex_t IT_0219 = m_t*s_13;
    const ccomplex_t IT_0220 = (-18)*conj(IT_0208);
    const ccomplex_t IT_0221 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = (-18)*IT_0208;
    const ccomplex_t IT_0224 = (-6)*conj(IT_0195);
    const ccomplex_t IT_0225 = ((-18)*conj(IT_0009) + (-6)*conj(IT_0112))
      *IT_0178 + ((-18)*IT_0009 + (-6)*IT_0112)*conj(IT_0178) + (-6)*IT_0195
      *IT_0212 + IT_0165*IT_0222 + (-6)*conj(IT_0165)*(IT_0173 + (
      -0.166666666666667)*IT_0223) + IT_0210*IT_0224;
    const ccomplex_t IT_0226 = m_t*m_u;
    const ccomplex_t IT_0227 = (IT_0009*conj(IT_0009) + IT_0112*conj(IT_0112) 
      + IT_0165*conj(IT_0165))*IT_0170 + IT_0171*IT_0213 + (-6)*(conj(IT_0165)
      *IT_0210 + IT_0165*IT_0212)*IT_0214 + (conj(IT_0178)*IT_0195 + IT_0178
      *conj(IT_0195) + conj(IT_0173)*IT_0208 + IT_0173*conj(IT_0208))*IT_0215 + 
      (IT_0173*conj(IT_0173) + IT_0178*conj(IT_0178) + IT_0195*conj(IT_0195) +
       IT_0208*conj(IT_0208))*IT_0216 + (conj(IT_0009)*IT_0112 + IT_0009*conj
      (IT_0112))*(IT_0217 + IT_0218) + IT_0219*IT_0225 + (-6)*(conj(IT_0178)
      *IT_0208 + -1./6*IT_0178*((-18)*conj(IT_0173) + (-6)*conj(IT_0208)) + -1.
      /6*IT_0195*IT_0222 + -1./6*conj(IT_0195)*IT_0223 + -1./6*IT_0173*((-18)
      *conj(IT_0178) + IT_0224))*IT_0226;
    return create_ccomplex_return(IT_0227);
}

