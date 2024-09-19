#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_su_R_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_su_R_to_u_Z(
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
    const ccomplex_t IT_0004 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0016 = 0.666666666666667*IT_0015;
    const ccomplex_t IT_0017 = pow(m_u, 2);
    const ccomplex_t IT_0018 = pow(m_su_R, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0014*IT_0021;
    const ccomplex_t IT_0023 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = IT_0003*IT_0023;
    const ccomplex_t IT_0026 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0027 = IT_0003*IT_0026;
    const ccomplex_t IT_0028 = IT_0007*IT_0026;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0024 + IT_0025 + 
      -IT_0027 + -IT_0028);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0037;
    const ccomplex_t IT_0040 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = IT_0003*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0038 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = -IT_0022 + (-0.5)*IT_0036 + -IT_0051;
    const ccomplex_t IT_0053 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = IT_0007*IT_0053;
    const ccomplex_t IT_0056 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = IT_0007*IT_0056;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0054 + IT_0055 + 
      -IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0062 = 0.666666666666667*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0060*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0029;
    const ccomplex_t IT_0068 = IT_0035*IT_0067;
    const ccomplex_t IT_0069 = -IT_0066 + (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = pow(m_W, -1);
    const ccomplex_t IT_0071 = sin(beta);
    const ccomplex_t IT_0072 = cpow(IT_0071, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0006*IT_0070*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0033*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0067*IT_0076;
    const ccomplex_t IT_0078 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0079 = IT_0003*IT_0078;
    const ccomplex_t IT_0080 = IT_0007*IT_0078;
    const ccomplex_t IT_0081 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0082 = IT_0007*IT_0081;
    const ccomplex_t IT_0083 = IT_0003*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0006*IT_0070*IT_0072;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0048*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0088*IT_0092;
    const ccomplex_t IT_0094 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0095 = IT_0007*IT_0094;
    const ccomplex_t IT_0096 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0097 = IT_0003*IT_0096;
    const ccomplex_t IT_0098 = IT_0007*IT_0096;
    const ccomplex_t IT_0099 = IT_0003*IT_0094;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0095 + -IT_0097 + 
      -IT_0098 + IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0006*IT_0070*IT_0072;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0019*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0101*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0077 + IT_0093 + IT_0106;
    const ccomplex_t IT_0108 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0109 = IT_0007*IT_0108;
    const ccomplex_t IT_0110 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = IT_0003*IT_0108;
    const ccomplex_t IT_0113 = IT_0003*IT_0110;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0109 + -IT_0111 +
       IT_0112 + -IT_0113);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0006*IT_0070*IT_0072;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0063*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0116*IT_0120;
    const ccomplex_t IT_0122 = IT_0030*IT_0076;
    const ccomplex_t IT_0123 = IT_0121 + 0.5*IT_0122;
    const ccomplex_t IT_0124 = pow(m_Z, -2);
    const ccomplex_t IT_0125 = pow(m_Z, 4);
    const ccomplex_t IT_0126 = s_13*IT_0125;
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = (-6)*IT_0127;
    const ccomplex_t IT_0129 = s_14*s_34;
    const ccomplex_t IT_0130 = (-24)*IT_0129;
    const ccomplex_t IT_0131 = pow(m_Z, 2);
    const ccomplex_t IT_0132 = s_13*IT_0131;
    const ccomplex_t IT_0133 = 12*IT_0132;
    const ccomplex_t IT_0134 = IT_0128 + IT_0130 + IT_0133;
    const ccomplex_t IT_0135 = 6*IT_0127;
    const ccomplex_t IT_0136 = 24*IT_0129;
    const ccomplex_t IT_0137 = (-12)*IT_0132;
    const ccomplex_t IT_0138 = IT_0135 + IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = (-2)*IT_0129;
    const ccomplex_t IT_0140 = IT_0132 + IT_0139;
    const ccomplex_t IT_0141 = IT_0124*IT_0140;
    const ccomplex_t IT_0142 = (-6)*IT_0141;
    const ccomplex_t IT_0143 = 12*s_13;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = pow(m_N_2, 2);
    const ccomplex_t IT_0146 = cpow((-2)*s_12 + IT_0017 + -IT_0018 + -IT_0145 
      + -reg_prop, -1);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0148 = e_em*IT_0003;
    const ccomplex_t IT_0149 = e_em*IT_0007;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0148 + (-3)*IT_0149);
    const ccomplex_t IT_0151 = (-0.166666666666667)*IT_0150;
    const ccomplex_t IT_0152 = IT_0032*IT_0151;
    const ccomplex_t IT_0153 = IT_0146*IT_0147*IT_0152;
    const ccomplex_t IT_0154 = m_N_2*IT_0030;
    const ccomplex_t IT_0155 = m_N_2*IT_0067;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = IT_0076*IT_0157;
    const ccomplex_t IT_0159 = m_N_3*IT_0014;
    const ccomplex_t IT_0160 = m_N_2*IT_0101;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = IT_0105*IT_0162;
    const ccomplex_t IT_0164 = m_N_4*IT_0045;
    const ccomplex_t IT_0165 = m_N_2*IT_0088;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = IT_0092*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0170 = (-0.666666666666667)*IT_0169;
    const ccomplex_t IT_0171 = IT_0032*IT_0170;
    const ccomplex_t IT_0172 = IT_0146*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_u*IT_0173;
    const ccomplex_t IT_0175 = -IT_0153 + (-0.5)*IT_0158 + -IT_0163 + -IT_0168
       + IT_0174;
    const ccomplex_t IT_0176 = m_N_1*IT_0060;
    const ccomplex_t IT_0177 = m_N_2*IT_0116;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = IT_0120*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0158 + IT_0180;
    const ccomplex_t IT_0182 = IT_0175 + IT_0181;
    const ccomplex_t IT_0183 = conj(IT_0175) + conj(IT_0181);
    const ccomplex_t IT_0184 = IT_0074*IT_0170;
    const ccomplex_t IT_0185 = IT_0146*IT_0147*IT_0184;
    const ccomplex_t IT_0186 = m_N_3*IT_0101;
    const ccomplex_t IT_0187 = m_N_2*IT_0014;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0186 + IT_0188;
    const ccomplex_t IT_0190 = IT_0021*IT_0189;
    const ccomplex_t IT_0191 = IT_0074*IT_0151;
    const ccomplex_t IT_0192 = IT_0146*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_u*IT_0193;
    const ccomplex_t IT_0195 = IT_0035*IT_0157;
    const ccomplex_t IT_0196 = m_N_4*IT_0088;
    const ccomplex_t IT_0197 = m_N_2*IT_0045;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0050*IT_0199;
    const ccomplex_t IT_0201 = IT_0185 + IT_0190 + -IT_0194 + (-0.5)*IT_0195 +
       IT_0200;
    const ccomplex_t IT_0202 = m_N_1*IT_0116;
    const ccomplex_t IT_0203 = m_N_2*IT_0060;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0065*IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0195 + -IT_0206;
    const ccomplex_t IT_0208 = IT_0201 + IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0201) + conj(IT_0207);
    const ccomplex_t IT_0210 = m_u*m_N_2;
    const ccomplex_t IT_0211 = pow(s_14, 2);
    const ccomplex_t IT_0212 = IT_0124*IT_0211;
    const ccomplex_t IT_0213 = -IT_0145;
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = IT_0210*IT_0214;
    const ccomplex_t IT_0216 = 1.33333333333333*IT_0169;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0017 + IT_0145 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0074*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = IT_0216*IT_0219;
    const ccomplex_t IT_0221 = -IT_0077 + (-2)*IT_0093 + (-2)*IT_0106 + 
      -IT_0220;
    const ccomplex_t IT_0222 = IT_0032*IT_0217;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = IT_0216*IT_0223;
    const ccomplex_t IT_0225 = 2*IT_0022 + IT_0036 + 2*IT_0051 + IT_0224;
    const ccomplex_t IT_0226 = 6*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0066 + IT_0068;
    const ccomplex_t IT_0228 = (-6)*IT_0227;
    const ccomplex_t IT_0229 = IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = 6*conj(IT_0225);
    const ccomplex_t IT_0231 = (-6)*conj(IT_0227);
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0121 + -IT_0122;
    const ccomplex_t IT_0234 = 6*IT_0227;
    const ccomplex_t IT_0235 = (-6)*IT_0225;
    const ccomplex_t IT_0236 = IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = 6*conj(IT_0227);
    const ccomplex_t IT_0238 = (-6)*conj(IT_0225);
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = m_u*IT_0214;
    const ccomplex_t IT_0241 = 6*IT_0221;
    const ccomplex_t IT_0242 = (-6)*IT_0233;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 6*conj(IT_0221);
    const ccomplex_t IT_0245 = (-6)*conj(IT_0233);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = s_14*s_34*IT_0124;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = s_13 + IT_0248;
    const ccomplex_t IT_0250 = s_13*IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0173 + -IT_0224;
    const ccomplex_t IT_0252 = 2*IT_0193 + IT_0220;
    const ccomplex_t IT_0253 = 6*IT_0233;
    const ccomplex_t IT_0254 = (-6)*IT_0221;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = 6*conj(IT_0233);
    const ccomplex_t IT_0257 = (-6)*conj(IT_0221);
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = s_13*s_14;
    const ccomplex_t IT_0260 = s_34*IT_0145;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = -IT_0173;
    const ccomplex_t IT_0264 = conj(IT_0069) + conj(IT_0263);
    const ccomplex_t IT_0265 = IT_0069 + IT_0263;
    const ccomplex_t IT_0266 = conj(IT_0123) + conj(IT_0193);
    const ccomplex_t IT_0267 = IT_0123 + IT_0193;
    const ccomplex_t IT_0268 = s_34*IT_0131;
    const ccomplex_t IT_0269 = IT_0124*IT_0268;
    const ccomplex_t IT_0270 = (-0.25)*IT_0269;
    const ccomplex_t IT_0271 = s_34 + IT_0270;
    const ccomplex_t IT_0272 = m_N_2*IT_0271;
    const ccomplex_t IT_0273 = (-4)*IT_0272;
    const ccomplex_t IT_0274 = 6*IT_0175;
    const ccomplex_t IT_0275 = 6*IT_0181;
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = 6*conj(IT_0175);
    const ccomplex_t IT_0278 = 6*conj(IT_0181);
    const ccomplex_t IT_0279 = IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = 6*IT_0201;
    const ccomplex_t IT_0281 = 6*IT_0207;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 6*conj(IT_0201);
    const ccomplex_t IT_0284 = 6*conj(IT_0207);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = (-6)*IT_0175;
    const ccomplex_t IT_0287 = (-6)*IT_0181;
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = (-6)*conj(IT_0175);
    const ccomplex_t IT_0290 = (-6)*conj(IT_0181);
    const ccomplex_t IT_0291 = IT_0289 + IT_0290;
    const ccomplex_t IT_0292 = (-6)*IT_0201;
    const ccomplex_t IT_0293 = (-6)*IT_0207;
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0296 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = s_14*IT_0131;
    const ccomplex_t IT_0299 = IT_0124*IT_0298;
    const ccomplex_t IT_0300 = 0.5*IT_0299;
    const ccomplex_t IT_0301 = s_14 + IT_0300;
    const ccomplex_t IT_0302 = m_u*IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0302;
    const ccomplex_t IT_0304 = IT_0131*IT_0210;
    const ccomplex_t IT_0305 = 3*IT_0304;
    const ccomplex_t IT_0306 = 6*IT_0263;
    const ccomplex_t IT_0307 = 6*conj(IT_0263);
    const ccomplex_t IT_0308 = (-18)*conj(IT_0175);
    const ccomplex_t IT_0309 = s_13*IT_0214;
    const ccomplex_t IT_0310 = pow(s_34, 2);
    const ccomplex_t IT_0311 = IT_0124*IT_0310;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0017 + IT_0312;
    const ccomplex_t IT_0314 = s_13*IT_0313;
    const ccomplex_t IT_0315 = (-6)*IT_0251;
    const ccomplex_t IT_0316 = (-6)*IT_0252;
    const ccomplex_t IT_0317 = IT_0210*IT_0313;
    const ccomplex_t IT_0318 = (-6)*conj(IT_0251);
    const ccomplex_t IT_0319 = IT_0210*IT_0249;
    const ccomplex_t IT_0320 = 6*IT_0251;
    const ccomplex_t IT_0321 = 6*conj(IT_0251);
    const ccomplex_t IT_0322 = m_u*IT_0249;
    const ccomplex_t IT_0323 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0324 = s_14*IT_0017;
    const ccomplex_t IT_0325 = s_13*s_34;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0324 + IT_0326;
    const ccomplex_t IT_0328 = conj(IT_0107)*IT_0252 + IT_0107*conj(IT_0252) +
       0.166666666666667*IT_0264*IT_0315 + 0.166666666666667*IT_0266*IT_0316 +
       0.166666666666667*IT_0265*IT_0318 + 0.166666666666667*conj(IT_0052)
      *IT_0320 + 0.166666666666667*IT_0052*IT_0321 + 0.166666666666667*IT_0267
      *IT_0323;
    const ccomplex_t IT_0329 = m_N_2*IT_0313;
    const ccomplex_t IT_0330 = m_N_2*IT_0249;
    const ccomplex_t IT_0331 = (-3)*IT_0304;
    const ccomplex_t IT_0332 = 6*IT_0052;
    const ccomplex_t IT_0333 = (-0.166666666666667)*IT_0332;
    const ccomplex_t IT_0334 = IT_0069 + IT_0333;
    const ccomplex_t IT_0335 = 6*conj(IT_0052);
    const ccomplex_t IT_0336 = (-0.166666666666667)*IT_0335;
    const ccomplex_t IT_0337 = conj(IT_0069) + IT_0336;
    const ccomplex_t IT_0338 = (-0.166666666666667)*conj(IT_0193);
    const ccomplex_t IT_0339 = IT_0128 + IT_0136;
    const ccomplex_t IT_0340 = IT_0130 + IT_0135;
    const ccomplex_t IT_0341 = (-0.25)*IT_0299;
    const ccomplex_t IT_0342 = s_14 + IT_0341;
    const ccomplex_t IT_0343 = m_u*IT_0342;
    const ccomplex_t IT_0344 = (-4)*IT_0343;
    const ccomplex_t IT_0345 = (-0.166666666666667)*IT_0344;
    const ccomplex_t IT_0346 = 0.5*IT_0269;
    const ccomplex_t IT_0347 = s_34 + IT_0346;
    const ccomplex_t IT_0348 = m_N_2*IT_0347;
    const ccomplex_t IT_0349 = 2*IT_0348;
    const ccomplex_t IT_0350 = (-0.166666666666667)*IT_0349;
    const ccomplex_t IT_0351 = (conj(IT_0052)*IT_0069 + IT_0052*conj(IT_0069) 
      + conj(IT_0107)*IT_0123 + IT_0107*conj(IT_0123))*IT_0134 + (IT_0052*conj
      (IT_0052) + IT_0069*conj(IT_0069) + IT_0107*conj(IT_0107) + IT_0123*conj
      (IT_0123))*IT_0138 + IT_0144*(IT_0182*IT_0183 + IT_0208*IT_0209) + IT_0215
      *(conj(IT_0221)*IT_0229 + IT_0221*IT_0232 + conj(IT_0233)*IT_0236 +
       IT_0233*IT_0239) + IT_0240*(IT_0183*IT_0229 + IT_0182*IT_0232 + IT_0209
      *IT_0243 + IT_0208*IT_0246) + IT_0250*(IT_0239*IT_0251 + IT_0236*conj
      (IT_0251) + conj(IT_0252)*IT_0255 + IT_0252*IT_0258) + IT_0262*(conj
      (IT_0052)*IT_0229 + IT_0052*IT_0232 + conj(IT_0107)*IT_0243 + IT_0107
      *IT_0246 + IT_0236*IT_0264 + IT_0239*IT_0265 + IT_0255*IT_0266 + IT_0258
      *IT_0267) + IT_0273*(conj(IT_0107)*IT_0276 + IT_0107*IT_0279 + conj
      (IT_0052)*IT_0282 + IT_0052*IT_0285 + conj(IT_0123)*IT_0288 + IT_0123
      *IT_0291 + conj(IT_0069)*IT_0294 + IT_0069*IT_0297) + (conj(IT_0052)
      *IT_0276 + IT_0052*IT_0279 + conj(IT_0107)*IT_0282 + IT_0107*IT_0285 +
       conj(IT_0069)*IT_0288 + IT_0069*IT_0291 + conj(IT_0123)*IT_0294 + IT_0123
      *IT_0297)*IT_0303 + IT_0305*(conj(IT_0107)*IT_0306 + IT_0107*IT_0307) +
       IT_0210*(IT_0182*((-18)*conj(IT_0201) + (-18)*conj(IT_0207)) + IT_0208*((
      -18)*conj(IT_0181) + IT_0308)) + (conj(IT_0225)*IT_0226 + conj(IT_0227)
      *IT_0236 + IT_0227*IT_0238 + conj(IT_0221)*IT_0241 + conj(IT_0233)*IT_0255
       + IT_0233*IT_0257)*IT_0309 + conj(IT_0251)*IT_0314*IT_0315 + conj(IT_0252
      )*(IT_0314*IT_0316 + IT_0315*IT_0317) + IT_0252*IT_0317*IT_0318 + IT_0319*
      (IT_0239*IT_0252 + IT_0236*conj(IT_0252) + conj(IT_0221)*IT_0315 + IT_0221
      *IT_0318 + conj(IT_0233)*IT_0320 + IT_0233*IT_0321) + IT_0322*(IT_0183
      *IT_0315 + IT_0209*IT_0316 + IT_0182*IT_0318 + IT_0208*IT_0323) + 6
      *IT_0327*IT_0328 + (conj(IT_0252)*IT_0286 + IT_0252*IT_0289 + IT_0209
      *IT_0315 + conj(IT_0181)*IT_0316 + IT_0208*IT_0318 + IT_0181*IT_0323)
      *IT_0329 + (IT_0209*IT_0236 + IT_0208*IT_0239 + conj(IT_0181)*IT_0255 +
       IT_0181*IT_0258 + conj(IT_0233)*IT_0274 + IT_0233*IT_0277 + conj(IT_0221)
      *IT_0286 + IT_0221*IT_0289)*IT_0330 + (-6)*IT_0331*((conj(IT_0052) + -conj
      (IT_0069))*IT_0123 + (IT_0052 + -IT_0069)*conj(IT_0123) + (
      -0.166666666666667)*IT_0193*IT_0307 + conj(IT_0107)*IT_0334 + IT_0107
      *IT_0337 + IT_0306*IT_0338) + (-6)*conj(IT_0193)*(IT_0305*IT_0334 + -1./6
      *IT_0123*IT_0339 + -1./6*IT_0107*IT_0340 + IT_0282*IT_0345 + IT_0276
      *IT_0350) + (-6)*IT_0193*(IT_0305*IT_0337 + IT_0138*IT_0338 + -1./6*conj
      (IT_0123)*IT_0339 + -1./6*conj(IT_0107)*IT_0340 + IT_0285*IT_0345 +
       IT_0279*IT_0350) + (-6)*conj(IT_0263)*((-0.166666666666667)*IT_0138
      *IT_0263 + IT_0123*IT_0305 + (-0.166666666666667)*IT_0069*IT_0339 + (
      -0.166666666666667)*IT_0052*IT_0340 + IT_0276*IT_0345 + IT_0282*IT_0350) +
       (-6)*IT_0263*(conj(IT_0123)*IT_0305 + (-0.166666666666667)*conj(IT_0069)
      *IT_0339 + (-0.166666666666667)*conj(IT_0052)*IT_0340 + IT_0279*IT_0345 +
       IT_0285*IT_0350);
    return create_ccomplex_return(IT_0351);
}

