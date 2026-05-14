#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_1_to_t_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_1_to_t_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0005 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0006*IT_0008*U_sb_00;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = pow(m_t, 2);
    const ccomplex_t IT_0012 = pow(m_C_1, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0018 = (-0.666666666666667)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0006*IT_0008*U_sb_01;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = pow(m_sb_1, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0011 + IT_0025 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0001*IT_0002;
    const ccomplex_t IT_0028 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0000*IT_0003;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0029 + -IT_0031 + (-2)
      *IT_0033);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0006*IT_0008*U_sb_00;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0011 + IT_0025 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0042 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = IT_0030*IT_0042;
    const ccomplex_t IT_0045 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0046 = IT_0030*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + -IT_0044 + (-2)
      *IT_0046);
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = IT_0041*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0016 + IT_0024 + 2*IT_0040 + 2*IT_0051;
    const ccomplex_t IT_0053 = pow(m_Z, -2);
    const ccomplex_t IT_0054 = pow(s_14, 2);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = -IT_0012;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = s_13*IT_0057;
    const ccomplex_t IT_0059 = m_t*m_C_1;
    const ccomplex_t IT_0060 = pow(m_Z, 2);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = 3*IT_0061;
    const ccomplex_t IT_0063 = -IT_0040 + -IT_0051;
    const ccomplex_t IT_0064 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0065 = IT_0003*IT_0064;
    const ccomplex_t IT_0066 = cos(beta);
    const ccomplex_t IT_0067 = cpow(IT_0066, -1);
    const ccomplex_t IT_0068 = IT_0003*IT_0067;
    const ccomplex_t IT_0069 = m_b*conj(U_d1)*V_tb*e_em*IT_0006*U_sb_10;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + (-0.5)*IT_0071);
    const ccomplex_t IT_0073 = e_em*IT_0027;
    const ccomplex_t IT_0074 = e_em*IT_0030;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0073 + (-3)*IT_0074);
    const ccomplex_t IT_0076 = (-0.166666666666667)*IT_0075;
    const ccomplex_t IT_0077 = IT_0072*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_12 + -IT_0012 + -IT_0025 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0082 = (-0.666666666666667)*IT_0081;
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = IT_0078*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = m_b*conj(U_d2)*V_tb*e_em*IT_0006*U_sb_10;
    const ccomplex_t IT_0090 = IT_0068*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-0.5)*IT_0091);
    const ccomplex_t IT_0093 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0094 = IT_0027*IT_0093;
    const ccomplex_t IT_0095 = IT_0030*IT_0093;
    const ccomplex_t IT_0096 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0097 = IT_0030*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + -IT_0095 + (-2)
      *IT_0097);
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = IT_0092*IT_0099;
    const ccomplex_t IT_0101 = IT_0026*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0105 = IT_0027*IT_0104;
    const ccomplex_t IT_0106 = IT_0030*IT_0104;
    const ccomplex_t IT_0107 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0108 = IT_0030*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + -IT_0106 + (-2)
      *IT_0108);
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0072*IT_0110;
    const ccomplex_t IT_0112 = IT_0041*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0103 + IT_0114;
    const ccomplex_t IT_0116 = conj(IT_0103) + conj(IT_0114);
    const ccomplex_t IT_0117 = (-3)*IT_0061;
    const ccomplex_t IT_0118 = pow(m_Z, 4);
    const ccomplex_t IT_0119 = s_13*IT_0118;
    const ccomplex_t IT_0120 = IT_0053*IT_0119;
    const ccomplex_t IT_0121 = s_14*s_34;
    const ccomplex_t IT_0122 = 24*IT_0121;
    const ccomplex_t IT_0123 = s_13*IT_0060;
    const ccomplex_t IT_0124 = (-12)*IT_0123;
    const ccomplex_t IT_0125 = 6*IT_0120 + IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = (-24)*IT_0121;
    const ccomplex_t IT_0127 = 6*IT_0120 + IT_0126;
    const ccomplex_t IT_0128 = -IT_0016 + -IT_0024 + (-2)*IT_0085;
    const ccomplex_t IT_0129 = s_14*s_34*IT_0053;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = s_13 + IT_0130;
    const ccomplex_t IT_0132 = s_13*IT_0131;
    const ccomplex_t IT_0133 = (-6)*conj(IT_0052);
    const ccomplex_t IT_0134 = (-6)*IT_0052;
    const ccomplex_t IT_0135 = pow(s_34, 2);
    const ccomplex_t IT_0136 = IT_0053*IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0011 + IT_0137;
    const ccomplex_t IT_0139 = s_13*IT_0138;
    const ccomplex_t IT_0140 = (-6)*IT_0128;
    const ccomplex_t IT_0141 = IT_0013*IT_0072;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0005*IT_0142;
    const ccomplex_t IT_0144 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0145 = IT_0003*IT_0144;
    const ccomplex_t IT_0146 = m_b*conj(U_d1)*V_tb*e_em*IT_0006*U_sb_11;
    const ccomplex_t IT_0147 = IT_0068*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + (-0.5)*IT_0148);
    const ccomplex_t IT_0150 = IT_0021*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0018*IT_0151;
    const ccomplex_t IT_0153 = 2*IT_0080 + IT_0143 + IT_0152;
    const ccomplex_t IT_0154 = IT_0059*IT_0131;
    const ccomplex_t IT_0155 = 2*IT_0102 + 2*IT_0113 + -IT_0143 + -IT_0152;
    const ccomplex_t IT_0156 = (-6)*IT_0155;
    const ccomplex_t IT_0157 = IT_0059*IT_0138;
    const ccomplex_t IT_0158 = (-6)*IT_0153;
    const ccomplex_t IT_0159 = m_C_1*IT_0138;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0161 = IT_0035*IT_0092;
    const ccomplex_t IT_0162 = IT_0026*IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0164 = IT_0010*IT_0076;
    const ccomplex_t IT_0165 = IT_0078*IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0167 = IT_0048*IT_0072;
    const ccomplex_t IT_0168 = IT_0041*IT_0166*IT_0167;
    const ccomplex_t IT_0169 = m_C_1*IT_0113;
    const ccomplex_t IT_0170 = m_C_1*IT_0102;
    const ccomplex_t IT_0171 = m_t*IT_0085;
    const ccomplex_t IT_0172 = -IT_0162 + -IT_0165 + -IT_0168 + -IT_0169 + 
      -IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = (-6)*IT_0172;
    const ccomplex_t IT_0174 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0175 = (-6)*conj(IT_0128);
    const ccomplex_t IT_0176 = (-6)*conj(IT_0172);
    const ccomplex_t IT_0177 = IT_0010*IT_0110;
    const ccomplex_t IT_0178 = IT_0041*IT_0166*IT_0177;
    const ccomplex_t IT_0179 = IT_0037*IT_0099;
    const ccomplex_t IT_0180 = IT_0026*IT_0160*IT_0179;
    const ccomplex_t IT_0181 = IT_0072*IT_0082;
    const ccomplex_t IT_0182 = IT_0078*IT_0163*IT_0181;
    const ccomplex_t IT_0183 = m_t*IT_0080;
    const ccomplex_t IT_0184 = m_C_1*IT_0051;
    const ccomplex_t IT_0185 = m_C_1*IT_0040;
    const ccomplex_t IT_0186 = -IT_0178 + -IT_0180 + IT_0182 + -IT_0183 + 
      -IT_0184 + -IT_0185;
    const ccomplex_t IT_0187 = s_14*IT_0060;
    const ccomplex_t IT_0188 = IT_0053*IT_0187;
    const ccomplex_t IT_0189 = (-0.25)*IT_0188;
    const ccomplex_t IT_0190 = s_14 + IT_0189;
    const ccomplex_t IT_0191 = m_t*IT_0190;
    const ccomplex_t IT_0192 = (-4)*IT_0191;
    const ccomplex_t IT_0193 = s_34*IT_0060;
    const ccomplex_t IT_0194 = IT_0053*IT_0193;
    const ccomplex_t IT_0195 = (-0.25)*IT_0194;
    const ccomplex_t IT_0196 = s_34 + IT_0195;
    const ccomplex_t IT_0197 = m_C_1*IT_0196;
    const ccomplex_t IT_0198 = (-4)*IT_0197;
    const ccomplex_t IT_0199 = 0.5*IT_0194;
    const ccomplex_t IT_0200 = s_34 + IT_0199;
    const ccomplex_t IT_0201 = m_C_1*IT_0200;
    const ccomplex_t IT_0202 = 2*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0188;
    const ccomplex_t IT_0204 = s_14 + IT_0203;
    const ccomplex_t IT_0205 = m_t*IT_0204;
    const ccomplex_t IT_0206 = 2*IT_0205;
    const ccomplex_t IT_0207 = m_C_1*IT_0131;
    const ccomplex_t IT_0208 = m_t*IT_0131;
    const ccomplex_t IT_0209 = m_t*IT_0057;
    const ccomplex_t IT_0210 = 6*IT_0155;
    const ccomplex_t IT_0211 = IT_0140*IT_0159 + (-18)*IT_0059*IT_0172 + 6
      *IT_0080*IT_0192 + 6*IT_0063*IT_0198 + 6*IT_0086*IT_0202 + 6*IT_0115
      *IT_0206 + IT_0134*IT_0207 + IT_0158*IT_0208 + IT_0209*IT_0210;
    const ccomplex_t IT_0212 = 12*s_13;
    const ccomplex_t IT_0213 = (-2)*IT_0121;
    const ccomplex_t IT_0214 = IT_0123 + IT_0213;
    const ccomplex_t IT_0215 = IT_0053*IT_0214;
    const ccomplex_t IT_0216 = (-6)*IT_0215;
    const ccomplex_t IT_0217 = IT_0212 + IT_0216;
    const ccomplex_t IT_0218 = (-6)*conj(IT_0153);
    const ccomplex_t IT_0219 = 6*conj(IT_0155);
    const ccomplex_t IT_0220 = (-18)*IT_0059*conj(IT_0172) + IT_0159*IT_0175 +
       6*conj(IT_0080)*IT_0192 + 6*conj(IT_0063)*IT_0198 + 6*conj(IT_0086)
      *IT_0202 + 6*IT_0116*IT_0206 + IT_0133*IT_0207 + conj(IT_0186)*IT_0217 +
       IT_0208*IT_0218 + IT_0209*IT_0219;
    const ccomplex_t IT_0221 = s_13*s_14;
    const ccomplex_t IT_0222 = s_34*IT_0012;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = IT_0221 + IT_0223;
    const ccomplex_t IT_0225 = 0.166666666666667*conj(IT_0086);
    const ccomplex_t IT_0226 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0227 = 0.166666666666667*conj(IT_0080);
    const ccomplex_t IT_0228 = s_14*IT_0011;
    const ccomplex_t IT_0229 = s_13*s_34;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = 0.166666666666667*IT_0208;
    const ccomplex_t IT_0233 = IT_0057*IT_0059;
    const ccomplex_t IT_0234 = 0.166666666666667*IT_0154;
    const ccomplex_t IT_0235 = 0.166666666666667*IT_0207;
    const ccomplex_t IT_0236 = 6*IT_0052*conj(IT_0052)*IT_0058 + 6*IT_0062*
      (conj(IT_0063)*IT_0080 + IT_0063*conj(IT_0080) + conj(IT_0086)*IT_0115 +
       IT_0086*IT_0116) + 6*(conj(IT_0080)*IT_0086 + IT_0080*conj(IT_0086) +
       conj(IT_0063)*IT_0115 + IT_0063*IT_0116)*IT_0117 + (IT_0063*conj(IT_0063)
       + IT_0080*conj(IT_0080) + IT_0086*conj(IT_0086) + (IT_0103 + IT_0114)
      *IT_0116)*IT_0125 + (conj(IT_0063)*IT_0086 + IT_0063*conj(IT_0086) + conj
      (IT_0080)*IT_0115 + IT_0080*IT_0116)*IT_0127 + IT_0128*IT_0132*IT_0133 +
       conj(IT_0128)*(IT_0132*IT_0134 + IT_0139*IT_0140) + conj(IT_0153)*
      (IT_0134*IT_0154 + IT_0132*IT_0156 + IT_0140*IT_0157 + IT_0139*IT_0158 +
       IT_0159*IT_0173) + IT_0153*(IT_0133*IT_0154 + IT_0132*IT_0174 + IT_0157
      *IT_0175 + IT_0159*IT_0176) + conj(IT_0186)*IT_0211 + IT_0186*IT_0220 + 6
      *IT_0224*(conj(IT_0052)*IT_0063 + IT_0052*conj(IT_0063) +
       0.166666666666667*IT_0086*IT_0133 + IT_0116*IT_0155 + IT_0115*conj
      (IT_0155) + IT_0134*IT_0225 + IT_0174*IT_0226 + IT_0156*IT_0227) + 6*(conj
      (IT_0063)*IT_0128 + IT_0063*conj(IT_0128) + IT_0116*IT_0153 + IT_0115*conj
      (IT_0153) + 0.166666666666667*IT_0086*IT_0175 + IT_0140*IT_0225 + IT_0218
      *IT_0226 + IT_0158*IT_0227)*IT_0231 + 6*conj(IT_0172)*(IT_0086*IT_0192 +
       IT_0115*IT_0198 + IT_0080*IT_0202 + IT_0063*IT_0206 + IT_0052*IT_0209 +
       IT_0140*IT_0232) + 6*IT_0172*(conj(IT_0086)*IT_0192 + IT_0116*IT_0198 +
       conj(IT_0080)*IT_0202 + conj(IT_0063)*IT_0206 + conj(IT_0052)*IT_0209 +
       0.166666666666667*conj(IT_0172)*IT_0217 + IT_0175*IT_0232) + IT_0219*
      (IT_0052*IT_0233 + IT_0140*IT_0234 + IT_0173*IT_0235) + IT_0210*(IT_0058
      *conj(IT_0155) + conj(IT_0052)*IT_0233 + IT_0175*IT_0234 + IT_0176*IT_0235);
    return create_ccomplex_return(IT_0236);
}

