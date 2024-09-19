#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_su_R_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_su_R_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = pow(m_u, 2);
    const ccomplex_t IT_0017 = pow(m_su_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = IT_0007*IT_0023;
    const ccomplex_t IT_0026 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0027 = IT_0003*IT_0026;
    const ccomplex_t IT_0028 = IT_0007*IT_0026;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0024 + IT_0025 + 
      -IT_0027 + -IT_0028);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0037;
    const ccomplex_t IT_0040 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0041 = IT_0003*IT_0040;
    const ccomplex_t IT_0042 = IT_0007*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0038 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0046 = 0.666666666666667*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0012;
    const ccomplex_t IT_0052 = IT_0020*IT_0051;
    const ccomplex_t IT_0053 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = IT_0007*IT_0053;
    const ccomplex_t IT_0056 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = IT_0007*IT_0056;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0054 + IT_0055 + 
      -IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0062 = 0.666666666666667*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0060*IT_0065;
    const ccomplex_t IT_0067 = -IT_0036 + -IT_0050 + (-0.5)*IT_0052 + -IT_0066;
    const ccomplex_t IT_0068 = pow(m_W, -1);
    const ccomplex_t IT_0069 = sin(beta);
    const ccomplex_t IT_0070 = cpow(IT_0069, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0006*IT_0068*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0018*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0051*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0013*IT_0074;
    const ccomplex_t IT_0078 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0079 = IT_0003*IT_0078;
    const ccomplex_t IT_0080 = IT_0007*IT_0078;
    const ccomplex_t IT_0081 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = IT_0007*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0006*IT_0068*IT_0070;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0033*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0085*IT_0089;
    const ccomplex_t IT_0091 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = IT_0007*IT_0091;
    const ccomplex_t IT_0094 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = IT_0007*IT_0094;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0092 + IT_0093 + 
      -IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0006*IT_0068*IT_0070;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0063*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0098*IT_0102;
    const ccomplex_t IT_0104 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0105 = IT_0003*IT_0104;
    const ccomplex_t IT_0106 = IT_0007*IT_0104;
    const ccomplex_t IT_0107 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0108 = IT_0003*IT_0107;
    const ccomplex_t IT_0109 = IT_0007*IT_0107;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0105 + IT_0106 + 
      -IT_0108 + -IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0006*IT_0068*IT_0070;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0047*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0111*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0077 + IT_0090 + IT_0103 + IT_0116;
    const ccomplex_t IT_0118 = pow(m_Z, -2);
    const ccomplex_t IT_0119 = pow(m_Z, 4);
    const ccomplex_t IT_0120 = s_13*IT_0119;
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = (-6)*IT_0121;
    const ccomplex_t IT_0123 = s_14*s_34;
    const ccomplex_t IT_0124 = (-24)*IT_0123;
    const ccomplex_t IT_0125 = pow(m_Z, 2);
    const ccomplex_t IT_0126 = s_13*IT_0125;
    const ccomplex_t IT_0127 = 12*IT_0126;
    const ccomplex_t IT_0128 = IT_0122 + IT_0124 + IT_0127;
    const ccomplex_t IT_0129 = 6*IT_0121;
    const ccomplex_t IT_0130 = 24*IT_0123;
    const ccomplex_t IT_0131 = (-12)*IT_0126;
    const ccomplex_t IT_0132 = IT_0129 + IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = (-2)*IT_0123;
    const ccomplex_t IT_0134 = IT_0126 + IT_0133;
    const ccomplex_t IT_0135 = IT_0118*IT_0134;
    const ccomplex_t IT_0136 = (-6)*IT_0135;
    const ccomplex_t IT_0137 = 12*s_13;
    const ccomplex_t IT_0138 = IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = e_em*IT_0003;
    const ccomplex_t IT_0140 = e_em*IT_0007;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0139 + (-3)*IT_0140);
    const ccomplex_t IT_0142 = (-0.166666666666667)*IT_0141;
    const ccomplex_t IT_0143 = IT_0015*IT_0142;
    const ccomplex_t IT_0144 = pow(m_N_4, 2);
    const ccomplex_t IT_0145 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0144 
      + -reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0147 = IT_0143*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = m_N_4*IT_0013;
    const ccomplex_t IT_0149 = m_N_4*IT_0051;
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = IT_0074*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0154 = (-0.666666666666667)*IT_0153;
    const ccomplex_t IT_0155 = IT_0015*IT_0154;
    const ccomplex_t IT_0156 = IT_0145*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = m_u*IT_0157;
    const ccomplex_t IT_0159 = -IT_0147 + (-0.5)*IT_0152 + IT_0158;
    const ccomplex_t IT_0160 = m_N_1*IT_0030;
    const ccomplex_t IT_0161 = m_N_4*IT_0085;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0160 + IT_0162;
    const ccomplex_t IT_0164 = IT_0089*IT_0163;
    const ccomplex_t IT_0165 = m_N_2*IT_0060;
    const ccomplex_t IT_0166 = m_N_4*IT_0098;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = IT_0102*IT_0168;
    const ccomplex_t IT_0170 = m_N_3*IT_0044;
    const ccomplex_t IT_0171 = m_N_4*IT_0111;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = IT_0115*IT_0173;
    const ccomplex_t IT_0175 = (-0.5)*IT_0152 + IT_0164 + IT_0169 + IT_0174;
    const ccomplex_t IT_0176 = IT_0159 + IT_0175;
    const ccomplex_t IT_0177 = conj(IT_0159) + conj(IT_0175);
    const ccomplex_t IT_0178 = IT_0072*IT_0154;
    const ccomplex_t IT_0179 = IT_0145*IT_0146*IT_0178;
    const ccomplex_t IT_0180 = IT_0020*IT_0151;
    const ccomplex_t IT_0181 = IT_0072*IT_0142;
    const ccomplex_t IT_0182 = IT_0145*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = m_u*IT_0183;
    const ccomplex_t IT_0185 = IT_0179 + (-0.5)*IT_0180 + -IT_0184;
    const ccomplex_t IT_0186 = m_N_1*IT_0085;
    const ccomplex_t IT_0187 = m_N_4*IT_0030;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0186 + IT_0188;
    const ccomplex_t IT_0190 = IT_0035*IT_0189;
    const ccomplex_t IT_0191 = m_N_2*IT_0098;
    const ccomplex_t IT_0192 = m_N_4*IT_0060;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0065*IT_0194;
    const ccomplex_t IT_0196 = m_N_3*IT_0111;
    const ccomplex_t IT_0197 = m_N_4*IT_0044;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0049*IT_0199;
    const ccomplex_t IT_0201 = (-0.5)*IT_0180 + -IT_0190 + -IT_0195 + -IT_0200;
    const ccomplex_t IT_0202 = IT_0185 + IT_0201;
    const ccomplex_t IT_0203 = conj(IT_0185) + conj(IT_0201);
    const ccomplex_t IT_0204 = m_u*m_N_4;
    const ccomplex_t IT_0205 = IT_0125*IT_0204;
    const ccomplex_t IT_0206 = (-3)*IT_0205;
    const ccomplex_t IT_0207 = 6*conj(IT_0022);
    const ccomplex_t IT_0208 = (-0.166666666666667)*IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0067) + IT_0208;
    const ccomplex_t IT_0210 = pow(s_14, 2);
    const ccomplex_t IT_0211 = IT_0118*IT_0210;
    const ccomplex_t IT_0212 = -IT_0144;
    const ccomplex_t IT_0213 = IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = IT_0204*IT_0213;
    const ccomplex_t IT_0215 = 1.33333333333333*IT_0153;
    const ccomplex_t IT_0216 = cpow((-2)*s_13 + IT_0016 + IT_0144 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0217 = IT_0072*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0215*IT_0218;
    const ccomplex_t IT_0220 = -IT_0075 + -IT_0219;
    const ccomplex_t IT_0221 = IT_0015*IT_0216;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0215*IT_0222;
    const ccomplex_t IT_0224 = IT_0021 + IT_0223;
    const ccomplex_t IT_0225 = 6*IT_0224;
    const ccomplex_t IT_0226 = 2*IT_0036 + 2*IT_0050 + IT_0052 + 2*IT_0066;
    const ccomplex_t IT_0227 = (-6)*IT_0226;
    const ccomplex_t IT_0228 = IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = 6*conj(IT_0224);
    const ccomplex_t IT_0230 = (-6)*conj(IT_0226);
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = -IT_0077 + (-2)*IT_0090 + (-2)*IT_0103 + (-2)
      *IT_0116;
    const ccomplex_t IT_0233 = 6*IT_0226;
    const ccomplex_t IT_0234 = (-6)*IT_0224;
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = 6*conj(IT_0226);
    const ccomplex_t IT_0237 = (-6)*conj(IT_0224);
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = m_u*IT_0213;
    const ccomplex_t IT_0240 = 6*IT_0220;
    const ccomplex_t IT_0241 = (-6)*IT_0232;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = 6*conj(IT_0220);
    const ccomplex_t IT_0244 = (-6)*conj(IT_0232);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = s_14*s_34*IT_0118;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = s_13 + IT_0247;
    const ccomplex_t IT_0249 = s_13*IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0157 + -IT_0223;
    const ccomplex_t IT_0251 = 2*IT_0183 + IT_0219;
    const ccomplex_t IT_0252 = 6*IT_0232;
    const ccomplex_t IT_0253 = (-6)*IT_0220;
    const ccomplex_t IT_0254 = IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = 6*conj(IT_0232);
    const ccomplex_t IT_0256 = (-6)*conj(IT_0220);
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = s_13*s_14;
    const ccomplex_t IT_0259 = s_34*IT_0144;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = -IT_0157;
    const ccomplex_t IT_0263 = conj(IT_0067) + conj(IT_0262);
    const ccomplex_t IT_0264 = IT_0067 + IT_0262;
    const ccomplex_t IT_0265 = conj(IT_0117) + conj(IT_0183);
    const ccomplex_t IT_0266 = IT_0117 + IT_0183;
    const ccomplex_t IT_0267 = s_34*IT_0125;
    const ccomplex_t IT_0268 = IT_0118*IT_0267;
    const ccomplex_t IT_0269 = (-0.25)*IT_0268;
    const ccomplex_t IT_0270 = s_34 + IT_0269;
    const ccomplex_t IT_0271 = m_N_4*IT_0270;
    const ccomplex_t IT_0272 = (-4)*IT_0271;
    const ccomplex_t IT_0273 = 6*IT_0159;
    const ccomplex_t IT_0274 = 6*IT_0175;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 6*conj(IT_0159);
    const ccomplex_t IT_0277 = 6*conj(IT_0175);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = 6*IT_0185;
    const ccomplex_t IT_0280 = 6*IT_0201;
    const ccomplex_t IT_0281 = IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = 6*conj(IT_0185);
    const ccomplex_t IT_0283 = 6*conj(IT_0201);
    const ccomplex_t IT_0284 = IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = (-6)*IT_0159;
    const ccomplex_t IT_0286 = (-6)*IT_0175;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = (-6)*conj(IT_0159);
    const ccomplex_t IT_0289 = (-6)*conj(IT_0175);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = (-6)*IT_0185;
    const ccomplex_t IT_0292 = (-6)*IT_0201;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = (-6)*conj(IT_0185);
    const ccomplex_t IT_0295 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_14*IT_0125;
    const ccomplex_t IT_0298 = IT_0118*IT_0297;
    const ccomplex_t IT_0299 = 0.5*IT_0298;
    const ccomplex_t IT_0300 = s_14 + IT_0299;
    const ccomplex_t IT_0301 = m_u*IT_0300;
    const ccomplex_t IT_0302 = 2*IT_0301;
    const ccomplex_t IT_0303 = 3*IT_0205;
    const ccomplex_t IT_0304 = IT_0122 + IT_0130;
    const ccomplex_t IT_0305 = IT_0124 + IT_0129;
    const ccomplex_t IT_0306 = 0.5*IT_0268;
    const ccomplex_t IT_0307 = s_34 + IT_0306;
    const ccomplex_t IT_0308 = m_N_4*IT_0307;
    const ccomplex_t IT_0309 = 2*IT_0308;
    const ccomplex_t IT_0310 = (-0.25)*IT_0298;
    const ccomplex_t IT_0311 = s_14 + IT_0310;
    const ccomplex_t IT_0312 = m_u*IT_0311;
    const ccomplex_t IT_0313 = (-4)*IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0262;
    const ccomplex_t IT_0315 = 6*conj(IT_0183)*((IT_0022 + -IT_0067)*IT_0303 +
       0.166666666666667*IT_0117*IT_0304 + 0.166666666666667*IT_0076*IT_0305 +
       0.166666666666667*IT_0275*IT_0309 + 0.166666666666667*IT_0281*IT_0313 +
       0.166666666666667*IT_0206*IT_0314);
    const ccomplex_t IT_0316 = 6*conj(IT_0262);
    const ccomplex_t IT_0317 = (-18)*conj(IT_0159);
    const ccomplex_t IT_0318 = s_13*IT_0213;
    const ccomplex_t IT_0319 = IT_0204*IT_0248;
    const ccomplex_t IT_0320 = 6*IT_0250;
    const ccomplex_t IT_0321 = 6*conj(IT_0250);
    const ccomplex_t IT_0322 = (-6)*IT_0250;
    const ccomplex_t IT_0323 = (-6)*conj(IT_0250);
    const ccomplex_t IT_0324 = pow(s_34, 2);
    const ccomplex_t IT_0325 = IT_0118*IT_0324;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0016 + IT_0326;
    const ccomplex_t IT_0328 = IT_0204*IT_0327;
    const ccomplex_t IT_0329 = s_13*IT_0327;
    const ccomplex_t IT_0330 = (-6)*IT_0251;
    const ccomplex_t IT_0331 = m_u*IT_0248;
    const ccomplex_t IT_0332 = (-6)*conj(IT_0251);
    const ccomplex_t IT_0333 = s_14*IT_0016;
    const ccomplex_t IT_0334 = s_13*s_34;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = conj(IT_0076)*IT_0251 + IT_0076*conj(IT_0251) +
       0.166666666666667*conj(IT_0022)*IT_0320 + 0.166666666666667*IT_0022
      *IT_0321 + 0.166666666666667*IT_0263*IT_0322 + 0.166666666666667*IT_0264
      *IT_0323 + 0.166666666666667*IT_0265*IT_0330 + 0.166666666666667*IT_0266
      *IT_0332;
    const ccomplex_t IT_0338 = m_N_4*IT_0248;
    const ccomplex_t IT_0339 = m_N_4*IT_0327;
    const ccomplex_t IT_0340 = (-0.166666666666667)*IT_0313;
    const ccomplex_t IT_0341 = (-0.166666666666667)*IT_0309;
    const ccomplex_t IT_0342 = (conj(IT_0022)*IT_0067 + IT_0022*conj(IT_0067) 
      + conj(IT_0076)*IT_0117 + IT_0076*conj(IT_0117))*IT_0128 + (IT_0022*conj
      (IT_0022) + IT_0067*conj(IT_0067) + IT_0076*conj(IT_0076) + IT_0117*conj
      (IT_0117))*IT_0132 + IT_0138*(IT_0176*IT_0177 + IT_0202*IT_0203) + (-6)
      *IT_0206*((conj(IT_0022) + -conj(IT_0067))*IT_0117 + -(IT_0022 + -IT_0067)
      *(conj(IT_0076) + -conj(IT_0117)) + IT_0076*IT_0209) + IT_0214*(conj
      (IT_0220)*IT_0228 + IT_0220*IT_0231 + conj(IT_0232)*IT_0235 + IT_0232
      *IT_0238) + IT_0239*(IT_0177*IT_0228 + IT_0176*IT_0231 + IT_0203*IT_0242 +
       IT_0202*IT_0245) + IT_0249*(IT_0238*IT_0250 + IT_0235*conj(IT_0250) +
       conj(IT_0251)*IT_0254 + IT_0251*IT_0257) + IT_0261*(conj(IT_0022)*IT_0228
       + IT_0022*IT_0231 + conj(IT_0076)*IT_0242 + IT_0076*IT_0245 + IT_0235
      *IT_0263 + IT_0238*IT_0264 + IT_0254*IT_0265 + IT_0257*IT_0266) + IT_0272*
      (conj(IT_0076)*IT_0275 + IT_0076*IT_0278 + conj(IT_0022)*IT_0281 + IT_0022
      *IT_0284 + conj(IT_0117)*IT_0287 + IT_0117*IT_0290 + conj(IT_0067)*IT_0293
       + IT_0067*IT_0296) + (conj(IT_0022)*IT_0275 + IT_0022*IT_0278 + conj
      (IT_0076)*IT_0281 + IT_0076*IT_0284 + conj(IT_0067)*IT_0287 + IT_0067
      *IT_0290 + conj(IT_0117)*IT_0293 + IT_0117*IT_0296)*IT_0302 + IT_0315 +
       IT_0303*(conj(IT_0076)*IT_0314 + IT_0076*IT_0316) + IT_0204*(IT_0176*((
      -18)*conj(IT_0185) + (-18)*conj(IT_0201)) + IT_0202*((-18)*conj(IT_0175) +
       IT_0317)) + (conj(IT_0224)*IT_0225 + conj(IT_0226)*IT_0235 + IT_0226
      *IT_0237 + conj(IT_0220)*IT_0240 + conj(IT_0232)*IT_0254 + IT_0232*IT_0256
      )*IT_0318 + IT_0319*(IT_0238*IT_0251 + IT_0235*conj(IT_0251) + conj
      (IT_0232)*IT_0320 + IT_0232*IT_0321 + conj(IT_0220)*IT_0322 + IT_0220
      *IT_0323) + (conj(IT_0251)*IT_0322 + IT_0251*IT_0323)*IT_0328 + IT_0329*
      (conj(IT_0250)*IT_0322 + conj(IT_0251)*IT_0330) + IT_0331*(IT_0177*IT_0322
       + IT_0176*IT_0323 + IT_0203*IT_0330 + IT_0202*IT_0332) + 6*IT_0336
      *IT_0337 + (IT_0203*IT_0235 + IT_0202*IT_0238 + conj(IT_0175)*IT_0254 +
       IT_0175*IT_0257 + conj(IT_0232)*IT_0273 + IT_0232*IT_0276 + conj(IT_0220)
      *IT_0285 + IT_0220*IT_0288)*IT_0338 + (conj(IT_0251)*IT_0285 + IT_0251
      *IT_0288 + IT_0203*IT_0322 + IT_0202*IT_0323 + conj(IT_0175)*IT_0330 +
       IT_0175*IT_0332)*IT_0339 + IT_0183*(IT_0132*conj(IT_0183) + (-6)*IT_0209
      *IT_0303 + conj(IT_0117)*IT_0304 + conj(IT_0076)*IT_0305 + IT_0206*IT_0316
       + (-6)*IT_0284*IT_0340 + (-6)*IT_0278*IT_0341) + (-6)*conj(IT_0262)*((
      -0.166666666666667)*IT_0132*IT_0262 + IT_0117*IT_0303 + (
      -0.166666666666667)*IT_0067*IT_0304 + (-0.166666666666667)*IT_0022*IT_0305
       + IT_0275*IT_0340 + IT_0281*IT_0341) + (-6)*IT_0262*(conj(IT_0117)
      *IT_0303 + (-0.166666666666667)*conj(IT_0067)*IT_0304 + (
      -0.166666666666667)*conj(IT_0022)*IT_0305 + IT_0278*IT_0340 + IT_0284
      *IT_0341);
    return create_ccomplex_return(IT_0342);
}

