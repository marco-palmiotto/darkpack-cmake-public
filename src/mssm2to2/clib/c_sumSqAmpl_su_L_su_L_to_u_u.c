#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_su_L_to_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_su_L_to_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
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
    const ccomplex_t IT_0000 = pow(m_su_L, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-18)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 36*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_B2)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = conj(N_W2)*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0010 + 3*IT_0015);
    const ccomplex_t IT_0017 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = pow(m_u, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = conj(N_B4)*e_em;
    const ccomplex_t IT_0030 = IT_0007*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W4)*e_em;
    const ccomplex_t IT_0033 = IT_0012*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 3*IT_0034);
    const ccomplex_t IT_0036 = 0.166666666666667*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = conj(N_B1)*e_em;
    const ccomplex_t IT_0046 = IT_0007*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W1)*e_em;
    const ccomplex_t IT_0049 = IT_0012*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + 3*IT_0050);
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0044*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B3)*e_em;
    const ccomplex_t IT_0060 = IT_0007*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em;
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0058*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0057 + -IT_0071;
    const ccomplex_t IT_0073 = IT_0028 + IT_0043 + IT_0072;
    const ccomplex_t IT_0074 = conj(IT_0028) + conj(IT_0043) + conj(IT_0072);
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0023*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0069*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0077 + IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0055*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0039*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = conj(IT_0081) + conj(IT_0084) + conj(IT_0087);
    const ccomplex_t IT_0090 = -IT_0027 + -IT_0042;
    const ccomplex_t IT_0091 = -IT_0071;
    const ccomplex_t IT_0092 = -IT_0057;
    const ccomplex_t IT_0093 = IT_0090 + IT_0091 + IT_0092;
    const ccomplex_t IT_0094 = conj(IT_0090) + conj(IT_0091) + conj(IT_0092);
    const ccomplex_t IT_0095 = IT_0080 + IT_0087;
    const ccomplex_t IT_0096 = IT_0077 + IT_0084;
    const ccomplex_t IT_0097 = conj(IT_0095) + conj(IT_0096);
    const ccomplex_t IT_0098 = (-6)*IT_0001;
    const ccomplex_t IT_0099 = 12*IT_0003;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = IT_0095 + IT_0096;
    const ccomplex_t IT_0102 = m_u*IT_0071;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = m_u*IT_0077;
    const ccomplex_t IT_0105 = cpow(IT_0022, 2);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_N_2*IT_0105;
    const ccomplex_t IT_0107 = IT_0075*IT_0106;
    const ccomplex_t IT_0108 = -IT_0104 + IT_0107;
    const ccomplex_t IT_0109 = 6*s_34;
    const ccomplex_t IT_0110 = conj(IT_0108)*IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0054, 2);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_N_1*IT_0111;
    const ccomplex_t IT_0113 = IT_0082*IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (-6)*s_34;
    const ccomplex_t IT_0116 = conj(IT_0114)*IT_0115;
    const ccomplex_t IT_0117 = IT_0110 + IT_0116;
    const ccomplex_t IT_0118 = cpow(IT_0068, 2);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_3*IT_0118;
    const ccomplex_t IT_0120 = IT_0078*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0038, 2);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_N_4*IT_0121;
    const ccomplex_t IT_0123 = IT_0085*IT_0122;
    const ccomplex_t IT_0124 = -IT_0120 + -IT_0123;
    const ccomplex_t IT_0125 = IT_0025*IT_0106;
    const ccomplex_t IT_0126 = IT_0058*IT_0119;
    const ccomplex_t IT_0127 = -IT_0125 + -IT_0126;
    const ccomplex_t IT_0128 = cpow(IT_0017, 2);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_2*IT_0128;
    const ccomplex_t IT_0130 = IT_0025*IT_0129;
    const ccomplex_t IT_0131 = cpow(IT_0066, 2);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_3*IT_0131;
    const ccomplex_t IT_0133 = IT_0058*IT_0132;
    const ccomplex_t IT_0134 = -IT_0130 + -IT_0133;
    const ccomplex_t IT_0135 = cpow(IT_0052, 2);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_1*IT_0135;
    const ccomplex_t IT_0137 = IT_0082*IT_0136;
    const ccomplex_t IT_0138 = IT_0075*IT_0129;
    const ccomplex_t IT_0139 = -IT_0137 + -IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0036, 2);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_N_4*IT_0140;
    const ccomplex_t IT_0142 = IT_0085*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0139 + IT_0143;
    const ccomplex_t IT_0145 = conj(IT_0139) + conj(IT_0143);
    const ccomplex_t IT_0146 = IT_0077 + IT_0080 + IT_0087;
    const ccomplex_t IT_0147 = m_u*IT_0146;
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = m_u*IT_0042;
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = IT_0109*IT_0150;
    const ccomplex_t IT_0152 = IT_0109*conj(IT_0150);
    const ccomplex_t IT_0153 = pow(m_sG, 2);
    const ccomplex_t IT_0154 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -IT_0153 +
       reg_prop, -1);
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0156 = cpow(IT_0155, 2);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_sG*IT_0156;
    const ccomplex_t IT_0158 = IT_0154*IT_0157;
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0159, 2);
    const ccomplex_t IT_0161 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -IT_0153 +
       reg_prop, -1);
    const ccomplex_t IT_0162 = IT_0157*IT_0161;
    const ccomplex_t IT_0163 = (-0.166666666666667)*IT_0162;
    const ccomplex_t IT_0164 = cpow(IT_0163, 2);
    const ccomplex_t IT_0165 = (-0.5)*IT_0158;
    const ccomplex_t IT_0166 = cpow(IT_0165, 2);
    const ccomplex_t IT_0167 = IT_0027 + IT_0057;
    const ccomplex_t IT_0168 = m_u*IT_0167;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = IT_0103 + IT_0169;
    const ccomplex_t IT_0171 = conj(IT_0103) + conj(IT_0169);
    const ccomplex_t IT_0172 = IT_0044*IT_0136;
    const ccomplex_t IT_0173 = IT_0040*IT_0141;
    const ccomplex_t IT_0174 = 0.5*IT_0162 + IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = IT_0027 + IT_0057 + IT_0071;
    const ccomplex_t IT_0176 = m_u*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0174 + IT_0177;
    const ccomplex_t IT_0179 = conj(IT_0174) + conj(IT_0177);
    const ccomplex_t IT_0180 = (-18)*s_34;
    const ccomplex_t IT_0181 = IT_0109*IT_0169;
    const ccomplex_t IT_0182 = IT_0114*IT_0180;
    const ccomplex_t IT_0183 = IT_0109*conj(IT_0169);
    const ccomplex_t IT_0184 = conj(IT_0114)*IT_0180;
    const ccomplex_t IT_0185 = IT_0150*IT_0180;
    const ccomplex_t IT_0186 = conj(IT_0150)*IT_0180;
    const ccomplex_t IT_0187 = m_u*IT_0087;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = m_u*IT_0080;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = m_u*IT_0084;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = conj(IT_0188) + conj(IT_0190) + conj(IT_0192);
    const ccomplex_t IT_0194 = IT_0103*IT_0109;
    const ccomplex_t IT_0195 = IT_0124*IT_0180;
    const ccomplex_t IT_0196 = IT_0188 + IT_0190 + IT_0192;
    const ccomplex_t IT_0197 = conj(IT_0103)*IT_0109;
    const ccomplex_t IT_0198 = conj(IT_0124)*IT_0180;
    const ccomplex_t IT_0199 = IT_0108 + IT_0188 + IT_0190 + IT_0192;
    const ccomplex_t IT_0200 = IT_0114*IT_0115;
    const ccomplex_t IT_0201 = IT_0108*IT_0109;
    const ccomplex_t IT_0202 = IT_0148*IT_0180;
    const ccomplex_t IT_0203 = IT_0115*IT_0150;
    const ccomplex_t IT_0204 = IT_0109*IT_0159;
    const ccomplex_t IT_0205 = IT_0163*IT_0180;
    const ccomplex_t IT_0206 = IT_0109*IT_0174;
    const ccomplex_t IT_0207 = IT_0151 + IT_0206;
    const ccomplex_t IT_0208 = 36*s_34;
    const ccomplex_t IT_0209 = conj(IT_0148)*IT_0180;
    const ccomplex_t IT_0210 = IT_0139*IT_0180;
    const ccomplex_t IT_0211 = IT_0143*IT_0180;
    const ccomplex_t IT_0212 = IT_0115*conj(IT_0150);
    const ccomplex_t IT_0213 = IT_0115*IT_0159;
    const ccomplex_t IT_0214 = IT_0109*IT_0165;
    const ccomplex_t IT_0215 = IT_0115*IT_0163;
    const ccomplex_t IT_0216 = IT_0115*IT_0124;
    const ccomplex_t IT_0217 = IT_0040*IT_0122;
    const ccomplex_t IT_0218 = IT_0044*IT_0112;
    const ccomplex_t IT_0219 = conj(IT_0217) + conj(IT_0218);
    const ccomplex_t IT_0220 = IT_0109*IT_0188;
    const ccomplex_t IT_0221 = IT_0109*IT_0190;
    const ccomplex_t IT_0222 = IT_0109*IT_0192;
    const ccomplex_t IT_0223 = IT_0127*IT_0180;
    const ccomplex_t IT_0224 = IT_0115*conj(IT_0124);
    const ccomplex_t IT_0225 = IT_0217 + IT_0218;
    const ccomplex_t IT_0226 = IT_0109*conj(IT_0188);
    const ccomplex_t IT_0227 = IT_0109*conj(IT_0190);
    const ccomplex_t IT_0228 = IT_0109*conj(IT_0192);
    const ccomplex_t IT_0229 = conj(IT_0127)*IT_0180;
    const ccomplex_t IT_0230 = IT_0109*IT_0148;
    const ccomplex_t IT_0231 = IT_0115*conj(IT_0148);
    const ccomplex_t IT_0232 = IT_0109*IT_0163;
    const ccomplex_t IT_0233 = IT_0159*IT_0180;
    const ccomplex_t IT_0234 = IT_0115*IT_0165;
    const ccomplex_t IT_0235 = IT_0115*IT_0139;
    const ccomplex_t IT_0236 = IT_0115*IT_0143;
    const ccomplex_t IT_0237 = IT_0109*conj(IT_0148);
    const ccomplex_t IT_0238 = IT_0115*conj(IT_0139);
    const ccomplex_t IT_0239 = IT_0115*conj(IT_0143);
    const ccomplex_t IT_0240 = IT_0078*IT_0132;
    const ccomplex_t IT_0241 = IT_0109*IT_0177;
    const ccomplex_t IT_0242 = IT_0165*IT_0180;
    const ccomplex_t IT_0243 = IT_0115*IT_0134;
    const ccomplex_t IT_0244 = IT_0109*conj(IT_0174);
    const ccomplex_t IT_0245 = IT_0109*conj(IT_0177);
    const ccomplex_t IT_0246 = conj(IT_0139)*IT_0180;
    const ccomplex_t IT_0247 = conj(IT_0143)*IT_0180;
    const ccomplex_t IT_0248 = IT_0115*conj(IT_0134);
    const ccomplex_t IT_0249 = 6*conj(IT_0114);
    const ccomplex_t IT_0250 = 6*conj(IT_0124);
    const ccomplex_t IT_0251 = 18*conj(IT_0127);
    const ccomplex_t IT_0252 = (-18)*conj(IT_0103);
    const ccomplex_t IT_0253 = (-18)*conj(IT_0150);
    const ccomplex_t IT_0254 = (-18)*conj(IT_0169);
    const ccomplex_t IT_0255 = (-6)*conj(IT_0108);
    const ccomplex_t IT_0256 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0257 = (-6)*conj(IT_0190);
    const ccomplex_t IT_0258 = (-6)*conj(IT_0192);
    const ccomplex_t IT_0259 = IT_0249 + IT_0250 + IT_0251 + IT_0252 + IT_0253
       + IT_0254 + IT_0255 + IT_0256 + IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = 6*IT_0103;
    const ccomplex_t IT_0261 = 6*IT_0150;
    const ccomplex_t IT_0262 = 6*IT_0169;
    const ccomplex_t IT_0263 = 18*IT_0108;
    const ccomplex_t IT_0264 = 18*IT_0188;
    const ccomplex_t IT_0265 = 18*IT_0190;
    const ccomplex_t IT_0266 = 18*IT_0192;
    const ccomplex_t IT_0267 = (-18)*IT_0114;
    const ccomplex_t IT_0268 = (-18)*IT_0124;
    const ccomplex_t IT_0269 = (-6)*IT_0127;
    const ccomplex_t IT_0270 = IT_0260 + IT_0261 + IT_0262 + IT_0263 + IT_0264
       + IT_0265 + IT_0266 + IT_0267 + IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = 6*conj(IT_0103);
    const ccomplex_t IT_0272 = 6*conj(IT_0150);
    const ccomplex_t IT_0273 = 6*conj(IT_0169);
    const ccomplex_t IT_0274 = 18*conj(IT_0108);
    const ccomplex_t IT_0275 = 18*conj(IT_0188);
    const ccomplex_t IT_0276 = 18*conj(IT_0190);
    const ccomplex_t IT_0277 = 18*conj(IT_0192);
    const ccomplex_t IT_0278 = (-18)*conj(IT_0114);
    const ccomplex_t IT_0279 = (-18)*conj(IT_0124);
    const ccomplex_t IT_0280 = (-6)*conj(IT_0127);
    const ccomplex_t IT_0281 = IT_0271 + IT_0272 + IT_0273 + IT_0274 + IT_0275
       + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = 6*IT_0108;
    const ccomplex_t IT_0283 = 6*IT_0188;
    const ccomplex_t IT_0284 = 6*IT_0190;
    const ccomplex_t IT_0285 = 6*IT_0192;
    const ccomplex_t IT_0286 = 18*IT_0103;
    const ccomplex_t IT_0287 = 18*IT_0150;
    const ccomplex_t IT_0288 = 18*IT_0169;
    const ccomplex_t IT_0289 = (-18)*IT_0127;
    const ccomplex_t IT_0290 = (-6)*IT_0114;
    const ccomplex_t IT_0291 = (-6)*IT_0124;
    const ccomplex_t IT_0292 = IT_0282 + IT_0283 + IT_0284 + IT_0285 + IT_0286
       + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = 6*conj(IT_0108);
    const ccomplex_t IT_0294 = 6*conj(IT_0188);
    const ccomplex_t IT_0295 = 6*conj(IT_0190);
    const ccomplex_t IT_0296 = 6*conj(IT_0192);
    const ccomplex_t IT_0297 = 18*conj(IT_0103);
    const ccomplex_t IT_0298 = 18*conj(IT_0150);
    const ccomplex_t IT_0299 = 18*conj(IT_0169);
    const ccomplex_t IT_0300 = (-18)*conj(IT_0127);
    const ccomplex_t IT_0301 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0302 = (-6)*conj(IT_0124);
    const ccomplex_t IT_0303 = IT_0293 + IT_0294 + IT_0295 + IT_0296 + IT_0297
       + IT_0298 + IT_0299 + IT_0300 + IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = 6*IT_0127;
    const ccomplex_t IT_0305 = 18*IT_0114;
    const ccomplex_t IT_0306 = 18*IT_0124;
    const ccomplex_t IT_0307 = (-18)*IT_0108;
    const ccomplex_t IT_0308 = (-18)*IT_0188;
    const ccomplex_t IT_0309 = (-18)*IT_0190;
    const ccomplex_t IT_0310 = (-18)*IT_0192;
    const ccomplex_t IT_0311 = (-6)*IT_0103;
    const ccomplex_t IT_0312 = (-6)*IT_0150;
    const ccomplex_t IT_0313 = (-6)*IT_0169;
    const ccomplex_t IT_0314 = IT_0304 + IT_0305 + IT_0306 + IT_0307 + IT_0308
       + IT_0309 + IT_0310 + IT_0311 + IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = IT_0148 + IT_0192 + IT_0240;
    const ccomplex_t IT_0316 = 6*conj(IT_0127);
    const ccomplex_t IT_0317 = 18*conj(IT_0114);
    const ccomplex_t IT_0318 = 18*conj(IT_0124);
    const ccomplex_t IT_0319 = (-18)*conj(IT_0108);
    const ccomplex_t IT_0320 = (-18)*conj(IT_0188);
    const ccomplex_t IT_0321 = (-18)*conj(IT_0190);
    const ccomplex_t IT_0322 = (-18)*conj(IT_0192);
    const ccomplex_t IT_0323 = (-6)*conj(IT_0103);
    const ccomplex_t IT_0324 = (-6)*conj(IT_0150);
    const ccomplex_t IT_0325 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0326 = IT_0316 + IT_0317 + IT_0318 + IT_0319 + IT_0320
       + IT_0321 + IT_0322 + IT_0323 + IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = 6*IT_0114;
    const ccomplex_t IT_0328 = 6*IT_0124;
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = (-18)*IT_0103;
    const ccomplex_t IT_0331 = (-18)*IT_0169;
    const ccomplex_t IT_0332 = (-6)*IT_0108;
    const ccomplex_t IT_0333 = IT_0163*(IT_0260 + IT_0261 + IT_0262 + IT_0263 
      + IT_0264 + IT_0265 + IT_0266 + IT_0267 + IT_0268 + IT_0269 + IT_0316 +
       IT_0317 + IT_0318 + IT_0319 + IT_0320 + IT_0321 + IT_0322 + IT_0323 +
       IT_0324 + IT_0325);
    const ccomplex_t IT_0334 = IT_0165*(IT_0271 + IT_0272 + IT_0273 + IT_0274 
      + IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280 + IT_0304 +
       IT_0305 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311 +
       IT_0312 + IT_0313);
    const ccomplex_t IT_0335 = 18*IT_0127;
    const ccomplex_t IT_0336 = (-18)*IT_0150;
    const ccomplex_t IT_0337 = (-6)*IT_0188;
    const ccomplex_t IT_0338 = (-6)*IT_0190;
    const ccomplex_t IT_0339 = (-6)*IT_0192;
    const ccomplex_t IT_0340 = IT_0159*(IT_0293 + IT_0294 + IT_0295 + IT_0296 
      + IT_0297 + IT_0298 + IT_0299 + IT_0300 + IT_0301 + IT_0302 + IT_0327 +
       IT_0328 + IT_0330 + IT_0331 + IT_0332 + IT_0335 + IT_0336 + IT_0337 +
       IT_0338 + IT_0339);
    const ccomplex_t IT_0341 = 6*IT_0165;
    const ccomplex_t IT_0342 = 6*IT_0139;
    const ccomplex_t IT_0343 = 6*IT_0143;
    const ccomplex_t IT_0344 = 18*IT_0159;
    const ccomplex_t IT_0345 = 18*IT_0134;
    const ccomplex_t IT_0346 = (-18)*IT_0174;
    const ccomplex_t IT_0347 = (-18)*IT_0177;
    const ccomplex_t IT_0348 = (-6)*IT_0163;
    const ccomplex_t IT_0349 = (-6)*IT_0148;
    const ccomplex_t IT_0350 = IT_0240 + (-0.166666666666667)*IT_0336 + (
      -0.166666666666667)*IT_0339 + (-0.166666666666667)*IT_0341 + (
      -0.166666666666667)*IT_0342 + (-0.166666666666667)*IT_0343 + (
      -0.166666666666667)*IT_0344 + (-0.166666666666667)*IT_0345 + (
      -0.166666666666667)*IT_0346 + (-0.166666666666667)*IT_0347 + (
      -0.166666666666667)*IT_0348 + (-0.166666666666667)*IT_0349;
    const ccomplex_t IT_0351 = 6*IT_0163;
    const ccomplex_t IT_0352 = 6*conj(IT_0139);
    const ccomplex_t IT_0353 = 6*conj(IT_0143);
    const ccomplex_t IT_0354 = 18*conj(IT_0134);
    const ccomplex_t IT_0355 = (-18)*IT_0159;
    const ccomplex_t IT_0356 = (-18)*conj(IT_0174);
    const ccomplex_t IT_0357 = (-18)*conj(IT_0177);
    const ccomplex_t IT_0358 = (-6)*IT_0165;
    const ccomplex_t IT_0359 = (-6)*conj(IT_0148);
    const ccomplex_t IT_0360 = conj(IT_0240) + (-0.166666666666667)*IT_0253 + 
      (-0.166666666666667)*IT_0258 + (-0.166666666666667)*IT_0351 + (
      -0.166666666666667)*IT_0352 + (-0.166666666666667)*IT_0353 + (
      -0.166666666666667)*IT_0354 + (-0.166666666666667)*IT_0355 + (
      -0.166666666666667)*IT_0356 + (-0.166666666666667)*IT_0357 + (
      -0.166666666666667)*IT_0358 + (-0.166666666666667)*IT_0359;
    const ccomplex_t IT_0361 = IT_0150*((-36)*conj(IT_0150) + IT_0249 +
       IT_0250 + IT_0252 + IT_0254 + IT_0255) + IT_0178*IT_0259 + IT_0145
      *IT_0270 + IT_0144*IT_0281 + conj(IT_0127)*IT_0287 + conj(IT_0134)*IT_0292
       + IT_0127*IT_0298 + IT_0134*IT_0303 + IT_0193*IT_0312 + (conj(IT_0148) +
       conj(IT_0192) + conj(IT_0240))*IT_0314 + IT_0196*IT_0324 + IT_0315
      *IT_0326 + conj(IT_0150)*(IT_0329 + IT_0330 + IT_0331 + IT_0332) + IT_0333
       + IT_0334 + IT_0179*(IT_0329 + IT_0330 + IT_0331 + IT_0332 + IT_0335 +
       IT_0336 + IT_0337 + IT_0338 + IT_0339) + IT_0340 + (-6)*IT_0219*IT_0350 +
       (-6)*IT_0225*IT_0360;
    const ccomplex_t IT_0362 = m_u*s_23;
    const ccomplex_t IT_0363 = (-18)*IT_0081;
    const ccomplex_t IT_0364 = (-6)*IT_0028;
    const ccomplex_t IT_0365 = (-6)*IT_0043;
    const ccomplex_t IT_0366 = (-6)*IT_0072;
    const ccomplex_t IT_0367 = IT_0363 + IT_0364 + IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = (-18)*conj(IT_0081);
    const ccomplex_t IT_0369 = (-6)*conj(IT_0028);
    const ccomplex_t IT_0370 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0371 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0372 = IT_0368 + IT_0369 + IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 6*IT_0159;
    const ccomplex_t IT_0374 = 6*IT_0134;
    const ccomplex_t IT_0375 = 18*IT_0165;
    const ccomplex_t IT_0376 = 18*IT_0139;
    const ccomplex_t IT_0377 = 18*IT_0143;
    const ccomplex_t IT_0378 = (-18)*IT_0163;
    const ccomplex_t IT_0379 = (-18)*IT_0148;
    const ccomplex_t IT_0380 = (-6)*IT_0174;
    const ccomplex_t IT_0381 = (-6)*IT_0177;
    const ccomplex_t IT_0382 = conj(IT_0084) + conj(IT_0087);
    const ccomplex_t IT_0383 = 6*conj(IT_0134);
    const ccomplex_t IT_0384 = 18*IT_0163;
    const ccomplex_t IT_0385 = 18*conj(IT_0139);
    const ccomplex_t IT_0386 = 18*conj(IT_0143);
    const ccomplex_t IT_0387 = (-18)*IT_0165;
    const ccomplex_t IT_0388 = (-18)*conj(IT_0148);
    const ccomplex_t IT_0389 = (-6)*IT_0159;
    const ccomplex_t IT_0390 = (-6)*conj(IT_0174);
    const ccomplex_t IT_0391 = (-6)*conj(IT_0177);
    const ccomplex_t IT_0392 = IT_0084 + IT_0087;
    const ccomplex_t IT_0393 = (-18)*IT_0090;
    const ccomplex_t IT_0394 = (-18)*IT_0091;
    const ccomplex_t IT_0395 = (-18)*IT_0092;
    const ccomplex_t IT_0396 = (-6)*IT_0095;
    const ccomplex_t IT_0397 = (-6)*IT_0096;
    const ccomplex_t IT_0398 = (-18)*conj(IT_0090);
    const ccomplex_t IT_0399 = (-18)*conj(IT_0091);
    const ccomplex_t IT_0400 = (-18)*conj(IT_0092);
    const ccomplex_t IT_0401 = (-6)*conj(IT_0095);
    const ccomplex_t IT_0402 = (-6)*conj(IT_0096);
    const ccomplex_t IT_0403 = IT_0093*IT_0259 + IT_0097*IT_0314 + IT_0101
      *IT_0326 + IT_0094*(IT_0329 + IT_0330 + IT_0331 + IT_0332 + IT_0335 +
       IT_0336 + IT_0337 + IT_0338 + IT_0339) + IT_0074*(IT_0336 + IT_0339 +
       IT_0341 + IT_0342 + IT_0343 + IT_0344 + IT_0345 + IT_0346 + IT_0347 +
       IT_0348 + IT_0349) + IT_0073*(IT_0253 + IT_0258 + IT_0351 + IT_0352 +
       IT_0353 + IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 + IT_0359) +
       conj(IT_0240)*IT_0367 + IT_0240*IT_0372 + conj(IT_0081)*(IT_0310 +
       IT_0312 + IT_0373 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 +
       IT_0379 + IT_0380 + IT_0381) + ((-18)*IT_0240 + IT_0310 + IT_0312 +
       IT_0373 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 + IT_0379 +
       IT_0380 + IT_0381)*IT_0382 + IT_0081*(IT_0322 + IT_0324 + IT_0383 +
       IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 +
       IT_0391) + ((-18)*conj(IT_0240) + IT_0322 + IT_0324 + IT_0383 + IT_0384 +
       IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 + IT_0391)
      *IT_0392 + IT_0219*(IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397) +
       IT_0225*(IT_0398 + IT_0399 + IT_0400 + IT_0401 + IT_0402);
    const ccomplex_t IT_0404 = IT_0000*IT_0024;
    const ccomplex_t IT_0405 = (-18)*IT_0028;
    const ccomplex_t IT_0406 = (-18)*IT_0043;
    const ccomplex_t IT_0407 = (-18)*IT_0072;
    const ccomplex_t IT_0408 = (-6)*IT_0081;
    const ccomplex_t IT_0409 = (-18)*conj(IT_0028);
    const ccomplex_t IT_0410 = (-18)*conj(IT_0043);
    const ccomplex_t IT_0411 = (-18)*conj(IT_0072);
    const ccomplex_t IT_0412 = (-6)*conj(IT_0081);
    const ccomplex_t IT_0413 = (-18)*IT_0095;
    const ccomplex_t IT_0414 = (-18)*IT_0096;
    const ccomplex_t IT_0415 = (-6)*IT_0090;
    const ccomplex_t IT_0416 = (-6)*IT_0091;
    const ccomplex_t IT_0417 = (-6)*IT_0092;
    const ccomplex_t IT_0418 = (-18)*conj(IT_0095);
    const ccomplex_t IT_0419 = (-18)*conj(IT_0096);
    const ccomplex_t IT_0420 = (-6)*conj(IT_0090);
    const ccomplex_t IT_0421 = (-6)*conj(IT_0091);
    const ccomplex_t IT_0422 = (-6)*conj(IT_0092);
    const ccomplex_t IT_0423 = IT_0097*IT_0367 + IT_0101*IT_0372 + IT_0094*
      (IT_0405 + IT_0406 + IT_0407 + IT_0408) + IT_0093*(IT_0409 + IT_0410 +
       IT_0411 + IT_0412) + IT_0382*(IT_0413 + IT_0414 + IT_0415 + IT_0416 +
       IT_0417) + IT_0392*(IT_0418 + IT_0419 + IT_0420 + IT_0421 + IT_0422);
    const ccomplex_t IT_0424 = 18*conj(IT_0240);
    const ccomplex_t IT_0425 = (IT_0114 + IT_0124)*(18*conj(IT_0114) + 18*conj
      (IT_0124)) + 18*IT_0127*conj(IT_0127) + 18*IT_0134*conj(IT_0134) + (18
      *conj(IT_0139) + 18*conj(IT_0143))*IT_0144 + IT_0159*(18*conj(IT_0134) +
       18*IT_0174 + 18*IT_0177) + (18*conj(IT_0174) + 18*conj(IT_0177))*IT_0178 
      + (18*conj(IT_0108) + 18*conj(IT_0188) + 18*conj(IT_0190) + 18*conj
      (IT_0192))*IT_0199 + conj(IT_0150)*(18*IT_0103 + 18*IT_0169 + 18*IT_0174 +
       18*IT_0177 + 18*IT_0217 + 18*IT_0218) + (IT_0170 + IT_0217 + IT_0218)*(18
      *conj(IT_0103) + 18*conj(IT_0169) + 18*conj(IT_0217) + 18*conj(IT_0218)) +
       IT_0150*(18*conj(IT_0103) + 18*IT_0159 + 18*conj(IT_0169) + 18*conj
      (IT_0174) + 18*conj(IT_0177) + 18*conj(IT_0217) + 18*conj(IT_0218)) +
       IT_0165*(18*conj(IT_0139) + 18*conj(IT_0143) + 18*IT_0148 + (-12)*IT_0159
       + 18*IT_0192 + 18*IT_0240) + IT_0315*(18*conj(IT_0148) + 18*conj(IT_0192)
       + IT_0424) + IT_0163*(18*IT_0139 + 18*IT_0143 + 18*conj(IT_0148) + 12
      *IT_0159 + 18*conj(IT_0192) + IT_0424);
    const ccomplex_t IT_0426 = m_u*s_24;
    const ccomplex_t IT_0427 = 6*IT_0081;
    const ccomplex_t IT_0428 = 6*IT_0087;
    const ccomplex_t IT_0429 = 18*IT_0028;
    const ccomplex_t IT_0430 = 18*IT_0043;
    const ccomplex_t IT_0431 = 18*IT_0072;
    const ccomplex_t IT_0432 = 6*conj(IT_0081);
    const ccomplex_t IT_0433 = 6*conj(IT_0087);
    const ccomplex_t IT_0434 = 18*conj(IT_0028);
    const ccomplex_t IT_0435 = 18*conj(IT_0043);
    const ccomplex_t IT_0436 = 18*conj(IT_0072);
    const ccomplex_t IT_0437 = IT_0225 + 0.166666666666667*IT_0260 +
       0.166666666666667*IT_0261 + 0.166666666666667*IT_0262 + 0.166666666666667
      *IT_0263 + 0.166666666666667*IT_0264 + 0.166666666666667*IT_0265 +
       0.166666666666667*IT_0266 + 0.166666666666667*IT_0267 + 0.166666666666667
      *IT_0268 + 0.166666666666667*IT_0269;
    const ccomplex_t IT_0438 = IT_0219 + 0.166666666666667*IT_0271 +
       0.166666666666667*IT_0272 + 0.166666666666667*IT_0273 + 0.166666666666667
      *IT_0274 + 0.166666666666667*IT_0275 + 0.166666666666667*IT_0276 +
       0.166666666666667*IT_0277 + 0.166666666666667*IT_0278 + 0.166666666666667
      *IT_0279 + 0.166666666666667*IT_0280;
    const ccomplex_t IT_0439 = 6*IT_0148;
    const ccomplex_t IT_0440 = 18*IT_0174;
    const ccomplex_t IT_0441 = 18*IT_0177;
    const ccomplex_t IT_0442 = (-18)*IT_0134;
    const ccomplex_t IT_0443 = (-6)*IT_0139;
    const ccomplex_t IT_0444 = (-6)*IT_0143;
    const ccomplex_t IT_0445 = IT_0240 + 0.166666666666667*IT_0285 +
       0.166666666666667*IT_0287 + 0.166666666666667*IT_0351 + 0.166666666666667
      *IT_0355 + 0.166666666666667*IT_0358 + 0.166666666666667*IT_0439 +
       0.166666666666667*IT_0440 + 0.166666666666667*IT_0441 + 0.166666666666667
      *IT_0442 + 0.166666666666667*IT_0443 + 0.166666666666667*IT_0444;
    const ccomplex_t IT_0446 = 6*conj(IT_0148);
    const ccomplex_t IT_0447 = 18*conj(IT_0174);
    const ccomplex_t IT_0448 = 18*conj(IT_0177);
    const ccomplex_t IT_0449 = (-18)*conj(IT_0134);
    const ccomplex_t IT_0450 = (-6)*conj(IT_0139);
    const ccomplex_t IT_0451 = (-6)*conj(IT_0143);
    const ccomplex_t IT_0452 = conj(IT_0240) + 0.166666666666667*IT_0296 +
       0.166666666666667*IT_0298 + 0.166666666666667*IT_0341 + 0.166666666666667
      *IT_0344 + 0.166666666666667*IT_0348 + 0.166666666666667*IT_0446 +
       0.166666666666667*IT_0447 + 0.166666666666667*IT_0448 + 0.166666666666667
      *IT_0449 + 0.166666666666667*IT_0450 + 0.166666666666667*IT_0451;
    const ccomplex_t IT_0453 = 6*IT_0174;
    const ccomplex_t IT_0454 = 6*IT_0177;
    const ccomplex_t IT_0455 = 18*IT_0148;
    const ccomplex_t IT_0456 = (-18)*IT_0139;
    const ccomplex_t IT_0457 = (-18)*IT_0143;
    const ccomplex_t IT_0458 = (-6)*IT_0134;
    const ccomplex_t IT_0459 = 6*conj(IT_0174);
    const ccomplex_t IT_0460 = 6*conj(IT_0177);
    const ccomplex_t IT_0461 = 18*conj(IT_0148);
    const ccomplex_t IT_0462 = (-18)*conj(IT_0139);
    const ccomplex_t IT_0463 = (-18)*conj(IT_0143);
    const ccomplex_t IT_0464 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0465 = (conj(IT_0081) + conj(IT_0087))*IT_0270 + 
      (IT_0081 + IT_0087)*IT_0281 + IT_0074*IT_0292 + IT_0073*IT_0303 + IT_0219*
      (IT_0427 + IT_0428 + IT_0429 + IT_0430 + IT_0431) + IT_0225*(IT_0432 +
       IT_0433 + IT_0434 + IT_0435 + IT_0436) + 6*conj(IT_0084)*IT_0437 + 6
      *IT_0084*IT_0438 + 6*IT_0094*IT_0445 + 6*IT_0093*IT_0452 + IT_0097*(18
      *IT_0240 + IT_0261 + IT_0266 + IT_0384 + IT_0387 + IT_0389 + IT_0453 +
       IT_0454 + IT_0455 + IT_0456 + IT_0457 + IT_0458) + IT_0101*(IT_0272 +
       IT_0277 + IT_0373 + IT_0375 + IT_0378 + IT_0424 + IT_0459 + IT_0460 +
       IT_0461 + IT_0462 + IT_0463 + IT_0464);
    const ccomplex_t IT_0466 = IT_0005*(IT_0073*IT_0074 + IT_0088*IT_0089 +
       IT_0093*IT_0094 + (IT_0095 + IT_0096)*IT_0097) + IT_0100*(IT_0074*IT_0088
       + IT_0073*IT_0089 + IT_0093*IT_0097 + IT_0094*IT_0101) + IT_0103*IT_0117 
      + IT_0109*((conj(IT_0114) + conj(IT_0124))*IT_0127 + (IT_0114 + IT_0124)
      *conj(IT_0127) + conj(IT_0134)*IT_0144 + IT_0134*IT_0145) + conj(IT_0148)
      *IT_0151 + IT_0148*IT_0152 + IT_0116*IT_0169 + (conj(IT_0108)*IT_0124 +
       IT_0108*conj(IT_0124) + IT_0160 + IT_0164 + IT_0166 + conj(IT_0127)
      *IT_0170 + IT_0127*IT_0171 + conj(IT_0134)*IT_0178 + IT_0134*IT_0179)
      *IT_0180 + conj(IT_0108)*(IT_0181 + IT_0182) + IT_0108*(IT_0183 + IT_0184)
       + conj(IT_0127)*IT_0185 + IT_0127*IT_0186 + IT_0193*(IT_0151 + IT_0181 +
       IT_0182 + IT_0194 + IT_0195) + IT_0196*(IT_0152 + IT_0183 + IT_0184 +
       IT_0197 + IT_0198) + IT_0115*(conj(IT_0134)*IT_0148 + conj(IT_0124)
      *IT_0170 + IT_0124*IT_0171 + IT_0165*conj(IT_0174) + IT_0127*(conj(IT_0108
      ) + conj(IT_0188) + conj(IT_0190) + conj(IT_0192)) + conj(IT_0127)*IT_0199
      ) + conj(IT_0169)*IT_0200 + conj(IT_0103)*(IT_0200 + IT_0201) + IT_0163*
      (IT_0152 + IT_0202 + IT_0203) + IT_0145*(IT_0202 + IT_0203 + IT_0204 +
       IT_0205) + IT_0165*(IT_0207 + IT_0163*IT_0208 + IT_0209 + IT_0210 +
       IT_0211 + IT_0212) + IT_0144*(IT_0209 + IT_0212 + IT_0213) + conj(IT_0134
      )*(IT_0185 + IT_0214 + IT_0215) + conj(IT_0150)*(IT_0200 + IT_0201 +
       IT_0150*IT_0208 + IT_0216) + IT_0219*(IT_0200 + IT_0201 + IT_0216 +
       IT_0220 + IT_0221 + IT_0222 + IT_0223) + IT_0150*(IT_0117 + IT_0224) +
       IT_0225*(IT_0110 + IT_0116 + IT_0224 + IT_0226 + IT_0227 + IT_0228 +
       IT_0229) + IT_0159*(IT_0186 + IT_0230 + IT_0231) + IT_0134*(IT_0186 +
       IT_0231 + IT_0232 + IT_0233 + IT_0234) + conj(IT_0174)*(IT_0230 + IT_0232
       + IT_0233 + IT_0235 + IT_0236) + conj(IT_0177)*(IT_0230 + IT_0232 +
       IT_0233 + IT_0234 + IT_0235 + IT_0236) + IT_0174*(IT_0215 + IT_0237 +
       IT_0238 + IT_0239) + IT_0177*(IT_0214 + IT_0215 + IT_0237 + IT_0238 +
       IT_0239) + (conj(IT_0192) + conj(IT_0240))*(IT_0207 + IT_0210 + IT_0211 +
       IT_0213 + IT_0241 + IT_0242 + IT_0243) + (IT_0192 + IT_0240)*(IT_0152 +
       IT_0204 + IT_0205 + IT_0244 + IT_0245 + IT_0246 + IT_0247 + IT_0248) +
       IT_0024*IT_0361 + IT_0362*IT_0403 + IT_0404*IT_0423 + s_34*IT_0425 +
       IT_0426*IT_0465;
    return create_ccomplex_return(IT_0466);
}

