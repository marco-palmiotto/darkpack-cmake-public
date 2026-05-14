#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sd_L_to_d_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sd_L_to_d_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_d, 2);
    const ccomplex_t IT_0020 = pow(m_sd_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0013*IT_0023;
    const ccomplex_t IT_0025 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*IT_0025;
    const ccomplex_t IT_0028 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = IT_0007*IT_0028;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0026 + IT_0027 + 
      -IT_0029 + -IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0024 + -IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0012;
    const ccomplex_t IT_0041 = IT_0023*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0043 = IT_0003*IT_0042;
    const ccomplex_t IT_0044 = IT_0007*IT_0042;
    const ccomplex_t IT_0045 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0046 = IT_0003*IT_0045;
    const ccomplex_t IT_0047 = IT_0007*IT_0045;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0043 + IT_0044 + 
      -IT_0046 + -IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0049*IT_0054;
    const ccomplex_t IT_0056 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = IT_0007*IT_0056;
    const ccomplex_t IT_0059 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = IT_0007*IT_0059;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0057 + IT_0058 + 
      -IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0041 + -IT_0055 + -IT_0069;
    const ccomplex_t IT_0071 = conj(N_B3)*e_em;
    const ccomplex_t IT_0072 = IT_0001*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W3)*e_em;
    const ccomplex_t IT_0075 = IT_0006*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-3)*IT_0076);
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0040*IT_0080;
    const ccomplex_t IT_0082 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = IT_0007*IT_0082;
    const ccomplex_t IT_0085 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = IT_0007*IT_0085;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0083 + IT_0084 + 
      -IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B4)*e_em;
    const ccomplex_t IT_0091 = IT_0001*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*e_em;
    const ccomplex_t IT_0094 = IT_0006*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + (-3)*IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = IT_0035*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0089*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0081 + IT_0100;
    const ccomplex_t IT_0102 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0103 = IT_0003*IT_0102;
    const ccomplex_t IT_0104 = IT_0007*IT_0102;
    const ccomplex_t IT_0105 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0106 = IT_0003*IT_0105;
    const ccomplex_t IT_0107 = IT_0007*IT_0105;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0103 + IT_0104 + 
      -IT_0106 + -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B1)*e_em;
    const ccomplex_t IT_0111 = IT_0001*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W1)*e_em;
    const ccomplex_t IT_0114 = IT_0006*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + (-3)*IT_0115);
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0118 = IT_0052*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0109*IT_0119;
    const ccomplex_t IT_0121 = IT_0013*IT_0080;
    const ccomplex_t IT_0122 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0123 = IT_0003*IT_0122;
    const ccomplex_t IT_0124 = IT_0007*IT_0122;
    const ccomplex_t IT_0125 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0126 = IT_0003*IT_0125;
    const ccomplex_t IT_0127 = IT_0007*IT_0125;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0123 + IT_0124 + 
      -IT_0126 + -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B2)*e_em;
    const ccomplex_t IT_0131 = IT_0001*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W2)*e_em;
    const ccomplex_t IT_0134 = IT_0006*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-3)*IT_0135);
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0138 = IT_0066*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = IT_0129*IT_0139;
    const ccomplex_t IT_0141 = IT_0120 + 0.5*IT_0121 + IT_0140;
    const ccomplex_t IT_0142 = pow(m_Z, -2);
    const ccomplex_t IT_0143 = pow(m_Z, 4);
    const ccomplex_t IT_0144 = s_13*IT_0143;
    const ccomplex_t IT_0145 = IT_0142*IT_0144;
    const ccomplex_t IT_0146 = (-6)*IT_0145;
    const ccomplex_t IT_0147 = s_14*s_34;
    const ccomplex_t IT_0148 = (-24)*IT_0147;
    const ccomplex_t IT_0149 = pow(m_Z, 2);
    const ccomplex_t IT_0150 = s_13*IT_0149;
    const ccomplex_t IT_0151 = 12*IT_0150;
    const ccomplex_t IT_0152 = IT_0146 + IT_0148 + IT_0151;
    const ccomplex_t IT_0153 = 6*IT_0145;
    const ccomplex_t IT_0154 = 24*IT_0147;
    const ccomplex_t IT_0155 = (-12)*IT_0150;
    const ccomplex_t IT_0156 = IT_0153 + IT_0154 + IT_0155;
    const ccomplex_t IT_0157 = (-2)*IT_0147;
    const ccomplex_t IT_0158 = IT_0150 + IT_0157;
    const ccomplex_t IT_0159 = IT_0142*IT_0158;
    const ccomplex_t IT_0160 = (-6)*IT_0159;
    const ccomplex_t IT_0161 = 12*s_13;
    const ccomplex_t IT_0162 = IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = e_em*IT_0003;
    const ccomplex_t IT_0164 = e_em*IT_0007;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0163 + 3*IT_0164);
    const ccomplex_t IT_0166 = (-0.166666666666667)*IT_0165;
    const ccomplex_t IT_0167 = IT_0018*IT_0166;
    const ccomplex_t IT_0168 = pow(m_N_3, 2);
    const ccomplex_t IT_0169 = cpow((-2)*s_12 + IT_0019 + -IT_0020 + -IT_0168 
      + -reg_prop, -1);
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0171 = IT_0167*IT_0169*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0173 = 0.333333333333333*IT_0172;
    const ccomplex_t IT_0174 = IT_0018*IT_0173;
    const ccomplex_t IT_0175 = IT_0169*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = m_d*IT_0176;
    const ccomplex_t IT_0178 = m_N_3*IT_0013;
    const ccomplex_t IT_0179 = m_N_3*IT_0040;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = IT_0080*IT_0181;
    const ccomplex_t IT_0183 = m_N_4*IT_0032;
    const ccomplex_t IT_0184 = m_N_3*IT_0089;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = IT_0099*IT_0186;
    const ccomplex_t IT_0188 = -IT_0171 + IT_0177 + (-0.5)*IT_0182 + -IT_0187;
    const ccomplex_t IT_0189 = m_N_2*IT_0063;
    const ccomplex_t IT_0190 = m_N_3*IT_0129;
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = IT_0139*IT_0192;
    const ccomplex_t IT_0194 = m_N_1*IT_0049;
    const ccomplex_t IT_0195 = m_N_3*IT_0109;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = IT_0119*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0182 + IT_0193 + IT_0198;
    const ccomplex_t IT_0200 = IT_0188 + IT_0199;
    const ccomplex_t IT_0201 = conj(IT_0188) + conj(IT_0199);
    const ccomplex_t IT_0202 = IT_0078*IT_0173;
    const ccomplex_t IT_0203 = IT_0169*IT_0170*IT_0202;
    const ccomplex_t IT_0204 = IT_0078*IT_0166;
    const ccomplex_t IT_0205 = IT_0169*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = m_d*IT_0206;
    const ccomplex_t IT_0208 = IT_0023*IT_0181;
    const ccomplex_t IT_0209 = m_N_4*IT_0089;
    const ccomplex_t IT_0210 = m_N_3*IT_0032;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = IT_0037*IT_0212;
    const ccomplex_t IT_0214 = IT_0203 + -IT_0207 + (-0.5)*IT_0208 + IT_0213;
    const ccomplex_t IT_0215 = m_N_1*IT_0109;
    const ccomplex_t IT_0216 = m_N_3*IT_0049;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0215 + IT_0217;
    const ccomplex_t IT_0219 = IT_0054*IT_0218;
    const ccomplex_t IT_0220 = m_N_2*IT_0129;
    const ccomplex_t IT_0221 = m_N_3*IT_0063;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = IT_0220 + IT_0222;
    const ccomplex_t IT_0224 = IT_0068*IT_0223;
    const ccomplex_t IT_0225 = (-0.5)*IT_0208 + -IT_0219 + -IT_0224;
    const ccomplex_t IT_0226 = IT_0214 + IT_0225;
    const ccomplex_t IT_0227 = conj(IT_0214) + conj(IT_0225);
    const ccomplex_t IT_0228 = m_d*m_N_3;
    const ccomplex_t IT_0229 = pow(s_14, 2);
    const ccomplex_t IT_0230 = IT_0142*IT_0229;
    const ccomplex_t IT_0231 = -IT_0168;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0228*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 3
      *IT_0000*IT_0006);
    const ccomplex_t IT_0235 = 0.333333333333333*IT_0234;
    const ccomplex_t IT_0236 = cpow((-2)*s_13 + IT_0019 + IT_0168 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0237 = IT_0078*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = IT_0235*IT_0238;
    const ccomplex_t IT_0240 = -IT_0081 + (-2)*IT_0100 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0018*IT_0236;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = IT_0235*IT_0242;
    const ccomplex_t IT_0244 = IT_0024 + 2*IT_0038 + IT_0243;
    const ccomplex_t IT_0245 = 6*IT_0244;
    const ccomplex_t IT_0246 = IT_0041 + 2*IT_0055 + 2*IT_0069;
    const ccomplex_t IT_0247 = (-6)*IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = 6*conj(IT_0244);
    const ccomplex_t IT_0250 = (-6)*conj(IT_0246);
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = (-2)*IT_0120 + -IT_0121 + (-2)*IT_0140;
    const ccomplex_t IT_0253 = 6*IT_0246;
    const ccomplex_t IT_0254 = (-6)*IT_0244;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = 6*conj(IT_0246);
    const ccomplex_t IT_0257 = (-6)*conj(IT_0244);
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = m_d*IT_0232;
    const ccomplex_t IT_0260 = 6*IT_0240;
    const ccomplex_t IT_0261 = (-6)*IT_0252;
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = 6*conj(IT_0240);
    const ccomplex_t IT_0264 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0265 = IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = s_14*s_34*IT_0142;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = s_13 + IT_0267;
    const ccomplex_t IT_0269 = s_13*IT_0268;
    const ccomplex_t IT_0270 = (-2)*IT_0176 + -IT_0243;
    const ccomplex_t IT_0271 = 2*IT_0206 + IT_0239;
    const ccomplex_t IT_0272 = 6*IT_0252;
    const ccomplex_t IT_0273 = (-6)*IT_0240;
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = 6*conj(IT_0252);
    const ccomplex_t IT_0276 = (-6)*conj(IT_0240);
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = s_13*s_14;
    const ccomplex_t IT_0279 = s_34*IT_0168;
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = IT_0278 + IT_0280;
    const ccomplex_t IT_0282 = -IT_0176;
    const ccomplex_t IT_0283 = conj(IT_0070) + conj(IT_0282);
    const ccomplex_t IT_0284 = conj(IT_0141) + conj(IT_0206);
    const ccomplex_t IT_0285 = IT_0141 + IT_0206;
    const ccomplex_t IT_0286 = s_34*IT_0149;
    const ccomplex_t IT_0287 = IT_0142*IT_0286;
    const ccomplex_t IT_0288 = (-0.25)*IT_0287;
    const ccomplex_t IT_0289 = s_34 + IT_0288;
    const ccomplex_t IT_0290 = m_N_3*IT_0289;
    const ccomplex_t IT_0291 = (-4)*IT_0290;
    const ccomplex_t IT_0292 = 6*IT_0188;
    const ccomplex_t IT_0293 = 6*IT_0199;
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = 6*conj(IT_0188);
    const ccomplex_t IT_0296 = 6*conj(IT_0199);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = 6*IT_0214;
    const ccomplex_t IT_0299 = 6*IT_0225;
    const ccomplex_t IT_0300 = IT_0298 + IT_0299;
    const ccomplex_t IT_0301 = 6*conj(IT_0214);
    const ccomplex_t IT_0302 = 6*conj(IT_0225);
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = (-6)*IT_0188;
    const ccomplex_t IT_0305 = (-6)*IT_0199;
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0308 = (-6)*conj(IT_0199);
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = (-6)*IT_0214;
    const ccomplex_t IT_0311 = (-6)*IT_0225;
    const ccomplex_t IT_0312 = IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = (-6)*conj(IT_0214);
    const ccomplex_t IT_0314 = (-6)*conj(IT_0225);
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = s_14*IT_0149;
    const ccomplex_t IT_0317 = IT_0142*IT_0316;
    const ccomplex_t IT_0318 = 0.5*IT_0317;
    const ccomplex_t IT_0319 = s_14 + IT_0318;
    const ccomplex_t IT_0320 = m_d*IT_0319;
    const ccomplex_t IT_0321 = 2*IT_0320;
    const ccomplex_t IT_0322 = IT_0149*IT_0228;
    const ccomplex_t IT_0323 = 3*IT_0322;
    const ccomplex_t IT_0324 = 6*IT_0282;
    const ccomplex_t IT_0325 = 6*conj(IT_0282);
    const ccomplex_t IT_0326 = (-18)*conj(IT_0188);
    const ccomplex_t IT_0327 = s_13*IT_0232;
    const ccomplex_t IT_0328 = IT_0228*IT_0268;
    const ccomplex_t IT_0329 = 6*IT_0270;
    const ccomplex_t IT_0330 = 6*conj(IT_0270);
    const ccomplex_t IT_0331 = (-6)*IT_0270;
    const ccomplex_t IT_0332 = (-6)*conj(IT_0270);
    const ccomplex_t IT_0333 = pow(s_34, 2);
    const ccomplex_t IT_0334 = IT_0142*IT_0333;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0019 + IT_0335;
    const ccomplex_t IT_0337 = IT_0228*IT_0336;
    const ccomplex_t IT_0338 = s_13*IT_0336;
    const ccomplex_t IT_0339 = (-6)*IT_0271;
    const ccomplex_t IT_0340 = m_d*IT_0268;
    const ccomplex_t IT_0341 = (-6)*conj(IT_0271);
    const ccomplex_t IT_0342 = s_14*IT_0019;
    const ccomplex_t IT_0343 = s_13*s_34;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = IT_0342 + IT_0344;
    const ccomplex_t IT_0346 = conj(IT_0101)*IT_0271 + IT_0101*conj(IT_0271) +
       0.166666666666667*conj(IT_0039)*IT_0329 + 0.166666666666667*IT_0039
      *IT_0330 + 0.166666666666667*IT_0283*IT_0331 + 0.166666666666667*(IT_0070 
      + IT_0282)*IT_0332 + 0.166666666666667*IT_0284*IT_0339 + 0.166666666666667
      *IT_0285*IT_0341;
    const ccomplex_t IT_0347 = m_N_3*IT_0268;
    const ccomplex_t IT_0348 = m_N_3*IT_0336;
    const ccomplex_t IT_0349 = (-3)*IT_0322;
    const ccomplex_t IT_0350 = 6*IT_0039;
    const ccomplex_t IT_0351 = (-0.166666666666667)*IT_0350;
    const ccomplex_t IT_0352 = IT_0070 + IT_0351;
    const ccomplex_t IT_0353 = 6*conj(IT_0039);
    const ccomplex_t IT_0354 = (-0.166666666666667)*IT_0353;
    const ccomplex_t IT_0355 = conj(IT_0070) + IT_0354;
    const ccomplex_t IT_0356 = (-0.166666666666667)*conj(IT_0206);
    const ccomplex_t IT_0357 = IT_0146 + IT_0154;
    const ccomplex_t IT_0358 = IT_0148 + IT_0153;
    const ccomplex_t IT_0359 = (-0.25)*IT_0317;
    const ccomplex_t IT_0360 = s_14 + IT_0359;
    const ccomplex_t IT_0361 = m_d*IT_0360;
    const ccomplex_t IT_0362 = (-4)*IT_0361;
    const ccomplex_t IT_0363 = (-0.166666666666667)*IT_0362;
    const ccomplex_t IT_0364 = 0.5*IT_0287;
    const ccomplex_t IT_0365 = s_34 + IT_0364;
    const ccomplex_t IT_0366 = m_N_3*IT_0365;
    const ccomplex_t IT_0367 = 2*IT_0366;
    const ccomplex_t IT_0368 = (-0.166666666666667)*IT_0367;
    const ccomplex_t IT_0369 = (conj(IT_0039)*IT_0070 + IT_0039*conj(IT_0070) 
      + conj(IT_0101)*IT_0141 + IT_0101*conj(IT_0141))*IT_0152 + (IT_0039*conj
      (IT_0039) + IT_0070*conj(IT_0070) + IT_0101*conj(IT_0101) + IT_0141*conj
      (IT_0141))*IT_0156 + IT_0162*(IT_0200*IT_0201 + IT_0226*IT_0227) + IT_0233
      *(conj(IT_0240)*IT_0248 + IT_0240*IT_0251 + conj(IT_0252)*IT_0255 +
       IT_0252*IT_0258) + IT_0259*(IT_0201*IT_0248 + IT_0200*IT_0251 + IT_0227
      *IT_0262 + IT_0226*IT_0265) + IT_0269*(IT_0258*IT_0270 + IT_0255*conj
      (IT_0270) + conj(IT_0271)*IT_0274 + IT_0271*IT_0277) + IT_0281*(conj
      (IT_0039)*IT_0248 + IT_0039*IT_0251 + conj(IT_0101)*IT_0262 + IT_0101
      *IT_0265 + IT_0258*(IT_0070 + IT_0282) + IT_0255*IT_0283 + IT_0274*IT_0284
       + IT_0277*IT_0285) + IT_0291*(conj(IT_0101)*IT_0294 + IT_0101*IT_0297 +
       conj(IT_0039)*IT_0300 + IT_0039*IT_0303 + conj(IT_0141)*IT_0306 + IT_0141
      *IT_0309 + conj(IT_0070)*IT_0312 + IT_0070*IT_0315) + (conj(IT_0039)
      *IT_0294 + IT_0039*IT_0297 + conj(IT_0101)*IT_0300 + IT_0101*IT_0303 +
       conj(IT_0070)*IT_0306 + IT_0070*IT_0309 + conj(IT_0141)*IT_0312 + IT_0141
      *IT_0315)*IT_0321 + IT_0323*(conj(IT_0101)*IT_0324 + IT_0101*IT_0325) +
       IT_0228*(IT_0200*((-18)*conj(IT_0214) + (-18)*conj(IT_0225)) + IT_0226*((
      -18)*conj(IT_0199) + IT_0326)) + (conj(IT_0244)*IT_0245 + conj(IT_0246)
      *IT_0255 + IT_0246*IT_0257 + conj(IT_0240)*IT_0260 + conj(IT_0252)*IT_0274
       + IT_0252*IT_0276)*IT_0327 + IT_0328*(IT_0258*IT_0271 + IT_0255*conj
      (IT_0271) + conj(IT_0252)*IT_0329 + IT_0252*IT_0330 + conj(IT_0240)
      *IT_0331 + IT_0240*IT_0332) + (conj(IT_0271)*IT_0331 + IT_0271*IT_0332)
      *IT_0337 + IT_0338*(conj(IT_0270)*IT_0331 + conj(IT_0271)*IT_0339) +
       IT_0340*(IT_0201*IT_0331 + IT_0200*IT_0332 + IT_0227*IT_0339 + IT_0226
      *IT_0341) + 6*IT_0345*IT_0346 + (IT_0227*IT_0255 + IT_0226*IT_0258 + conj
      (IT_0199)*IT_0274 + IT_0199*IT_0277 + conj(IT_0252)*IT_0292 + IT_0252
      *IT_0295 + conj(IT_0240)*IT_0304 + IT_0240*IT_0307)*IT_0347 + (conj
      (IT_0271)*IT_0304 + IT_0271*IT_0307 + IT_0227*IT_0331 + IT_0226*IT_0332 +
       conj(IT_0199)*IT_0339 + IT_0199*IT_0341)*IT_0348 + (-6)*IT_0349*((conj
      (IT_0039) + -conj(IT_0070))*IT_0141 + (IT_0039 + -IT_0070)*conj(IT_0141) +
       (-0.166666666666667)*IT_0206*IT_0325 + conj(IT_0101)*IT_0352 + IT_0101
      *IT_0355 + IT_0324*IT_0356) + (-6)*conj(IT_0206)*(IT_0323*IT_0352 + -1./6
      *IT_0141*IT_0357 + -1./6*IT_0101*IT_0358 + IT_0300*IT_0363 + IT_0294
      *IT_0368) + (-6)*IT_0206*(IT_0323*IT_0355 + IT_0156*IT_0356 + -1./6*conj
      (IT_0141)*IT_0357 + -1./6*conj(IT_0101)*IT_0358 + IT_0303*IT_0363 +
       IT_0297*IT_0368) + (-6)*conj(IT_0282)*((-0.166666666666667)*IT_0156
      *IT_0282 + IT_0141*IT_0323 + (-0.166666666666667)*IT_0070*IT_0357 + (
      -0.166666666666667)*IT_0039*IT_0358 + IT_0294*IT_0363 + IT_0300*IT_0368) +
       (-6)*IT_0282*(conj(IT_0141)*IT_0323 + (-0.166666666666667)*conj(IT_0070)
      *IT_0357 + (-0.166666666666667)*conj(IT_0039)*IT_0358 + IT_0297*IT_0363 +
       IT_0303*IT_0368);
    return create_ccomplex_return(IT_0369);
}

