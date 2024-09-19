#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_3_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_3_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = IT_0005 + IT_0011 + IT_0013 + IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0006;
    const ccomplex_t IT_0018 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0008;
    const ccomplex_t IT_0021 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0016 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = cpow(IT_0001, -4);
    const ccomplex_t IT_0033 = IT_0000*IT_0007*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = IT_0006*IT_0007*IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0007, -3);
    const ccomplex_t IT_0037 = IT_0000*IT_0031*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0034*IT_0036;
    const ccomplex_t IT_0039 = cpow(IT_0001, -2);
    const ccomplex_t IT_0040 = IT_0000*IT_0008*IT_0031*IT_0039;
    const ccomplex_t IT_0041 = IT_0006*IT_0008*IT_0034*IT_0039;
    const ccomplex_t IT_0042 = cpow(IT_0007, -2);
    const ccomplex_t IT_0043 = cpow(IT_0039 + IT_0042, -1);
    const ccomplex_t IT_0044 = (IT_0033 + IT_0035 + IT_0037 + IT_0038 + 2
      *IT_0040 + 2*IT_0041)*IT_0043;
    const ccomplex_t IT_0045 = m_W*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_N_1, 2);
    const ccomplex_t IT_0049 = pow(m_N_3, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + -IT_0048 + -IT_0049 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0030*IT_0052;
    const ccomplex_t IT_0054 = IT_0003*IT_0018;
    const ccomplex_t IT_0055 = IT_0004*IT_0017;
    const ccomplex_t IT_0056 = IT_0003*IT_0023;
    const ccomplex_t IT_0057 = IT_0012*IT_0017;
    const ccomplex_t IT_0058 = IT_0054 + IT_0055 + IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = IT_0010*IT_0020;
    const ccomplex_t IT_0060 = IT_0009*IT_0021;
    const ccomplex_t IT_0061 = IT_0014*IT_0020;
    const ccomplex_t IT_0062 = IT_0009*IT_0025;
    const ccomplex_t IT_0063 = -IT_0059 + -IT_0060 + -IT_0061 + -IT_0062;
    const ccomplex_t IT_0064 = IT_0058 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0000*IT_0007*IT_0032*IT_0034;
    const ccomplex_t IT_0068 = IT_0006*IT_0007*IT_0031*IT_0032;
    const ccomplex_t IT_0069 = IT_0000*IT_0034*IT_0036;
    const ccomplex_t IT_0070 = IT_0006*IT_0031*IT_0036;
    const ccomplex_t IT_0071 = IT_0000*IT_0008*IT_0034*IT_0039;
    const ccomplex_t IT_0072 = IT_0006*IT_0008*IT_0031*IT_0039;
    const ccomplex_t IT_0073 = IT_0043*(IT_0067 + -IT_0068 + IT_0069 + 
      -IT_0070 + 2*IT_0071 + (-2)*IT_0072);
    const ccomplex_t IT_0074 = IT_0045*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0048 + -IT_0049 + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0066*IT_0078;
    const ccomplex_t IT_0080 = -IT_0053 + -IT_0079;
    const ccomplex_t IT_0081 = pow(m_Z, -4);
    const ccomplex_t IT_0082 = pow(s_34, 2);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = 1 + IT_0084;
    const ccomplex_t IT_0086 = s_12*IT_0085;
    const ccomplex_t IT_0087 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0090 = IT_0009*IT_0089;
    const ccomplex_t IT_0091 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0094 = IT_0009*IT_0093;
    const ccomplex_t IT_0095 = IT_0088 + IT_0090 + IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0097 = IT_0017*IT_0096;
    const ccomplex_t IT_0098 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0099 = IT_0020*IT_0098;
    const ccomplex_t IT_0100 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0101 = IT_0017*IT_0100;
    const ccomplex_t IT_0102 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0103 = IT_0020*IT_0102;
    const ccomplex_t IT_0104 = -IT_0097 + -IT_0099 + -IT_0101 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0095 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0052*IT_0107;
    const ccomplex_t IT_0109 = IT_0003*IT_0096;
    const ccomplex_t IT_0110 = IT_0017*IT_0087;
    const ccomplex_t IT_0111 = IT_0003*IT_0100;
    const ccomplex_t IT_0112 = IT_0017*IT_0091;
    const ccomplex_t IT_0113 = IT_0109 + IT_0110 + IT_0111 + IT_0112;
    const ccomplex_t IT_0114 = IT_0020*IT_0089;
    const ccomplex_t IT_0115 = IT_0009*IT_0098;
    const ccomplex_t IT_0116 = IT_0020*IT_0093;
    const ccomplex_t IT_0117 = IT_0009*IT_0102;
    const ccomplex_t IT_0118 = -IT_0114 + -IT_0115 + -IT_0116 + -IT_0117;
    const ccomplex_t IT_0119 = IT_0113 + IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0078*IT_0121;
    const ccomplex_t IT_0123 = -IT_0108 + -IT_0122;
    const ccomplex_t IT_0124 = m_N_1*m_N_3;
    const ccomplex_t IT_0125 = IT_0085*IT_0124;
    const ccomplex_t IT_0126 = pow(m_Z, -2);
    const ccomplex_t IT_0127 = s_13*s_34*IT_0126;
    const ccomplex_t IT_0128 = pow(m_Z, 4);
    const ccomplex_t IT_0129 = s_14*IT_0128;
    const ccomplex_t IT_0130 = pow(m_Z, 2);
    const ccomplex_t IT_0131 = s_13*s_34*IT_0130;
    const ccomplex_t IT_0132 = (-2)*IT_0131;
    const ccomplex_t IT_0133 = IT_0129 + IT_0132;
    const ccomplex_t IT_0134 = IT_0081*IT_0133;
    const ccomplex_t IT_0135 = (-0.666666666666667)*IT_0127 + (
      -0.166666666666667)*IT_0134;
    const ccomplex_t IT_0136 = s_14 + IT_0135;
    const ccomplex_t IT_0137 = m_N_3*IT_0136;
    const ccomplex_t IT_0138 = (-6)*IT_0137;
    const ccomplex_t IT_0139 = IT_0002*IT_0007;
    const ccomplex_t IT_0140 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = IT_0001*IT_0008;
    const ccomplex_t IT_0143 = IT_0140*IT_0142;
    const ccomplex_t IT_0144 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0145 = IT_0139*IT_0144;
    const ccomplex_t IT_0146 = IT_0142*IT_0144;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0141 + IT_0143 + 
      -IT_0145 + -IT_0146);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0150 = IT_0139*IT_0149;
    const ccomplex_t IT_0151 = IT_0142*IT_0149;
    const ccomplex_t IT_0152 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0153 = IT_0139*IT_0152;
    const ccomplex_t IT_0154 = IT_0142*IT_0152;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0150 + IT_0151 + 
      -IT_0153 + -IT_0154);
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0148*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_23 + IT_0049 + IT_0130 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0157*IT_0159;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0163 = IT_0139*IT_0162;
    const ccomplex_t IT_0164 = IT_0142*IT_0162;
    const ccomplex_t IT_0165 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0166 = IT_0139*IT_0165;
    const ccomplex_t IT_0167 = IT_0142*IT_0165;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0163 + IT_0164 + 
      -IT_0166 + -IT_0167);
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0171 = IT_0139*IT_0170;
    const ccomplex_t IT_0172 = IT_0142*IT_0170;
    const ccomplex_t IT_0173 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0174 = IT_0139*IT_0173;
    const ccomplex_t IT_0175 = IT_0142*IT_0173;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0171 + IT_0172 + 
      -IT_0174 + -IT_0175);
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0169*IT_0177;
    const ccomplex_t IT_0179 = cpow((-2)*s_23 + IT_0049 + IT_0130 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0178*IT_0180;
    const ccomplex_t IT_0182 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0183 = IT_0139*IT_0182;
    const ccomplex_t IT_0184 = IT_0142*IT_0182;
    const ccomplex_t IT_0185 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0186 = IT_0139*IT_0185;
    const ccomplex_t IT_0187 = IT_0142*IT_0185;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0183 + IT_0184 + 
      -IT_0186 + -IT_0187);
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0191 = IT_0139*IT_0190;
    const ccomplex_t IT_0192 = IT_0142*IT_0190;
    const ccomplex_t IT_0193 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0194 = IT_0139*IT_0193;
    const ccomplex_t IT_0195 = IT_0142*IT_0193;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0191 + IT_0192 + 
      -IT_0194 + -IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0189*IT_0197;
    const ccomplex_t IT_0199 = cpow((-2)*s_23 + IT_0049 + IT_0130 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0198*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0188;
    const ccomplex_t IT_0203 = IT_0197*IT_0202;
    const ccomplex_t IT_0204 = IT_0200*IT_0203;
    const ccomplex_t IT_0205 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0206 = IT_0139*IT_0205;
    const ccomplex_t IT_0207 = IT_0142*IT_0205;
    const ccomplex_t IT_0208 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0209 = IT_0139*IT_0208;
    const ccomplex_t IT_0210 = IT_0142*IT_0208;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0206 + IT_0207 + 
      -IT_0209 + -IT_0210);
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = IT_0197*IT_0212;
    const ccomplex_t IT_0214 = cpow((-2)*s_23 + IT_0130 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0213*IT_0215;
    const ccomplex_t IT_0217 = 0.5*IT_0201 + (-0.5)*IT_0204 + -IT_0216;
    const ccomplex_t IT_0218 = conj(IT_0161) + conj(IT_0181) + conj(IT_0217);
    const ccomplex_t IT_0219 = cpow((-2)*s_13 + IT_0048 + IT_0130 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0222 = IT_0139*IT_0221;
    const ccomplex_t IT_0223 = IT_0142*IT_0221;
    const ccomplex_t IT_0224 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0225 = IT_0139*IT_0224;
    const ccomplex_t IT_0226 = IT_0142*IT_0224;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0222 + IT_0223 + 
      -IT_0225 + -IT_0226);
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0169*IT_0228;
    const ccomplex_t IT_0230 = m_N_2*IT_0229;
    const ccomplex_t IT_0231 = m_N_3*IT_0178;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = IT_0220*IT_0233;
    const ccomplex_t IT_0235 = cpow((-2)*s_13 + IT_0130 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0238 = IT_0139*IT_0237;
    const ccomplex_t IT_0239 = IT_0142*IT_0237;
    const ccomplex_t IT_0240 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0241 = IT_0139*IT_0240;
    const ccomplex_t IT_0242 = IT_0142*IT_0240;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*(IT_0238 + IT_0239 + 
      -IT_0241 + -IT_0242);
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = IT_0189*IT_0244;
    const ccomplex_t IT_0246 = m_N_1*IT_0245;
    const ccomplex_t IT_0247 = m_N_3*IT_0198;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = IT_0246 + IT_0248;
    const ccomplex_t IT_0250 = IT_0236*IT_0249;
    const ccomplex_t IT_0251 = -IT_0234 + -IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0251;
    const ccomplex_t IT_0253 = cpow((-2)*s_13 + IT_0048 + IT_0130 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = m_N_3*IT_0157;
    const ccomplex_t IT_0256 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0257 = IT_0139*IT_0256;
    const ccomplex_t IT_0258 = IT_0142*IT_0256;
    const ccomplex_t IT_0259 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0260 = IT_0139*IT_0259;
    const ccomplex_t IT_0261 = IT_0142*IT_0259;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*(IT_0257 + IT_0258 + 
      -IT_0260 + -IT_0261);
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = IT_0148*IT_0263;
    const ccomplex_t IT_0265 = m_N_4*IT_0264;
    const ccomplex_t IT_0266 = IT_0254*(IT_0255 + -IT_0265);
    const ccomplex_t IT_0267 = (-2)*IT_0266;
    const ccomplex_t IT_0268 = (-0.5)*IT_0211;
    const ccomplex_t IT_0269 = IT_0197*IT_0268;
    const ccomplex_t IT_0270 = m_N_3*IT_0269;
    const ccomplex_t IT_0271 = m_N_3*IT_0213;
    const ccomplex_t IT_0272 = cpow((-2)*s_13 + IT_0048 + IT_0130 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = (IT_0270 + -IT_0271)*IT_0273;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = (-2)*IT_0275;
    const ccomplex_t IT_0277 = IT_0252 + IT_0267 + IT_0276;
    const ccomplex_t IT_0278 = IT_0161 + IT_0181 + IT_0217;
    const ccomplex_t IT_0279 = 2*conj(IT_0251);
    const ccomplex_t IT_0280 = (-2)*conj(IT_0266);
    const ccomplex_t IT_0281 = (-2)*conj(IT_0275);
    const ccomplex_t IT_0282 = IT_0279 + IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0284 = IT_0139*IT_0283;
    const ccomplex_t IT_0285 = IT_0142*IT_0283;
    const ccomplex_t IT_0286 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0287 = IT_0139*IT_0286;
    const ccomplex_t IT_0288 = IT_0142*IT_0286;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*(IT_0284 + IT_0285 + 
      -IT_0287 + -IT_0288);
    const ccomplex_t IT_0290 = (-0.5)*IT_0289;
    const ccomplex_t IT_0291 = IT_0263*IT_0290;
    const ccomplex_t IT_0292 = IT_0159*IT_0291;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0295 = IT_0139*IT_0294;
    const ccomplex_t IT_0296 = IT_0142*IT_0294;
    const ccomplex_t IT_0297 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0298 = IT_0139*IT_0297;
    const ccomplex_t IT_0299 = IT_0142*IT_0297;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*(IT_0295 + IT_0296 + 
      -IT_0298 + -IT_0299);
    const ccomplex_t IT_0301 = (-0.5)*IT_0300;
    const ccomplex_t IT_0302 = IT_0228*IT_0301;
    const ccomplex_t IT_0303 = IT_0180*IT_0302;
    const ccomplex_t IT_0304 = IT_0244*IT_0268;
    const ccomplex_t IT_0305 = IT_0215*IT_0304;
    const ccomplex_t IT_0306 = IT_0303 + -IT_0305;
    const ccomplex_t IT_0307 = IT_0202*IT_0244;
    const ccomplex_t IT_0308 = -IT_0245;
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = IT_0200*IT_0309;
    const ccomplex_t IT_0311 = 0.5*IT_0310;
    const ccomplex_t IT_0312 = conj(IT_0293) + conj(IT_0306) + conj(IT_0311);
    const ccomplex_t IT_0313 = IT_0177*IT_0301;
    const ccomplex_t IT_0314 = m_N_2*IT_0313;
    const ccomplex_t IT_0315 = m_N_3*IT_0302;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = IT_0220*IT_0317;
    const ccomplex_t IT_0319 = m_N_1*IT_0203;
    const ccomplex_t IT_0320 = m_N_3*IT_0307;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0319 + IT_0321;
    const ccomplex_t IT_0323 = IT_0236*IT_0322;
    const ccomplex_t IT_0324 = -IT_0318 + -IT_0323;
    const ccomplex_t IT_0325 = 2*IT_0324;
    const ccomplex_t IT_0326 = IT_0212*IT_0244;
    const ccomplex_t IT_0327 = m_N_3*IT_0326;
    const ccomplex_t IT_0328 = m_N_3*IT_0304;
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = IT_0327 + IT_0329;
    const ccomplex_t IT_0331 = IT_0273*IT_0330;
    const ccomplex_t IT_0332 = IT_0156*IT_0290;
    const ccomplex_t IT_0333 = m_N_4*IT_0332;
    const ccomplex_t IT_0334 = m_N_3*IT_0291;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = IT_0254*IT_0336;
    const ccomplex_t IT_0338 = -IT_0331 + -IT_0337;
    const ccomplex_t IT_0339 = (-2)*IT_0338;
    const ccomplex_t IT_0340 = IT_0325 + IT_0339;
    const ccomplex_t IT_0341 = IT_0293 + IT_0306 + IT_0311;
    const ccomplex_t IT_0342 = 2*conj(IT_0324);
    const ccomplex_t IT_0343 = (-2)*conj(IT_0338);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = s_14*IT_0130;
    const ccomplex_t IT_0346 = s_13*s_34;
    const ccomplex_t IT_0347 = (-2)*IT_0346;
    const ccomplex_t IT_0348 = IT_0345 + IT_0347;
    const ccomplex_t IT_0349 = IT_0126*IT_0348;
    const ccomplex_t IT_0350 = 0.166666666666667*IT_0134 + (-0.666666666666667
      )*IT_0349;
    const ccomplex_t IT_0351 = s_14 + IT_0350;
    const ccomplex_t IT_0352 = m_N_3*IT_0351;
    const ccomplex_t IT_0353 = 6*IT_0352;
    const ccomplex_t IT_0354 = IT_0178*IT_0220;
    const ccomplex_t IT_0355 = IT_0198*IT_0236;
    const ccomplex_t IT_0356 = -IT_0354 + -IT_0355;
    const ccomplex_t IT_0357 = IT_0157*IT_0254;
    const ccomplex_t IT_0358 = IT_0213*IT_0273;
    const ccomplex_t IT_0359 = conj(IT_0357) + conj(IT_0358);
    const ccomplex_t IT_0360 = conj(IT_0356) + IT_0359;
    const ccomplex_t IT_0361 = IT_0357 + IT_0358;
    const ccomplex_t IT_0362 = IT_0356 + IT_0361;
    const ccomplex_t IT_0363 = IT_0254*IT_0291;
    const ccomplex_t IT_0364 = IT_0273*IT_0304;
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = IT_0220*IT_0302;
    const ccomplex_t IT_0367 = IT_0236*IT_0307;
    const ccomplex_t IT_0368 = -IT_0366 + -IT_0367;
    const ccomplex_t IT_0369 = conj(IT_0365) + conj(IT_0368);
    const ccomplex_t IT_0370 = IT_0365 + IT_0368;
    const ccomplex_t IT_0371 = m_N_1*IT_0126;
    const ccomplex_t IT_0372 = s_24*IT_0130;
    const ccomplex_t IT_0373 = s_23*s_34;
    const ccomplex_t IT_0374 = s_24*IT_0128;
    const ccomplex_t IT_0375 = IT_0126*IT_0374;
    const ccomplex_t IT_0376 = (-2)*IT_0373 + (-0.5)*IT_0375;
    const ccomplex_t IT_0377 = IT_0372 + IT_0376;
    const ccomplex_t IT_0378 = IT_0371*IT_0377;
    const ccomplex_t IT_0379 = (-2)*IT_0378;
    const ccomplex_t IT_0380 = m_N_3*IT_0126;
    const ccomplex_t IT_0381 = IT_0126*IT_0129;
    const ccomplex_t IT_0382 = (-2)*IT_0346 + (-0.5)*IT_0381;
    const ccomplex_t IT_0383 = IT_0345 + IT_0382;
    const ccomplex_t IT_0384 = IT_0380*IT_0383;
    const ccomplex_t IT_0385 = (-2)*IT_0384;
    const ccomplex_t IT_0386 = m_N_1*IT_0291;
    const ccomplex_t IT_0387 = IT_0159*IT_0386;
    const ccomplex_t IT_0388 = 2*IT_0387;
    const ccomplex_t IT_0389 = m_N_1*IT_0302;
    const ccomplex_t IT_0390 = -IT_0389;
    const ccomplex_t IT_0391 = IT_0230 + IT_0390;
    const ccomplex_t IT_0392 = IT_0180*IT_0391;
    const ccomplex_t IT_0393 = m_N_1*IT_0304;
    const ccomplex_t IT_0394 = -IT_0393;
    const ccomplex_t IT_0395 = IT_0270 + IT_0394;
    const ccomplex_t IT_0396 = IT_0215*IT_0395;
    const ccomplex_t IT_0397 = IT_0159*IT_0265;
    const ccomplex_t IT_0398 = -IT_0392 + IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = (-2)*IT_0398;
    const ccomplex_t IT_0400 = m_N_1*IT_0307;
    const ccomplex_t IT_0401 = IT_0200*(IT_0246 + -IT_0400);
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = IT_0388 + IT_0399 + IT_0403;
    const ccomplex_t IT_0405 = 2*conj(IT_0387);
    const ccomplex_t IT_0406 = (-2)*conj(IT_0398);
    const ccomplex_t IT_0407 = (-2)*conj(IT_0402);
    const ccomplex_t IT_0408 = IT_0405 + IT_0406 + IT_0407;
    const ccomplex_t IT_0409 = m_N_1*IT_0157;
    const ccomplex_t IT_0410 = IT_0159*IT_0409;
    const ccomplex_t IT_0411 = 2*IT_0410;
    const ccomplex_t IT_0412 = m_N_1*IT_0198;
    const ccomplex_t IT_0413 = IT_0200*(IT_0319 + -IT_0412);
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = (-2)*IT_0414;
    const ccomplex_t IT_0416 = m_N_1*IT_0178;
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = IT_0314 + IT_0417;
    const ccomplex_t IT_0419 = IT_0180*IT_0418;
    const ccomplex_t IT_0420 = m_N_1*IT_0213;
    const ccomplex_t IT_0421 = -IT_0420;
    const ccomplex_t IT_0422 = IT_0327 + IT_0421;
    const ccomplex_t IT_0423 = IT_0215*IT_0422;
    const ccomplex_t IT_0424 = IT_0159*IT_0333;
    const ccomplex_t IT_0425 = -IT_0419 + IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = (-2)*IT_0425;
    const ccomplex_t IT_0427 = IT_0411 + IT_0415 + IT_0426;
    const ccomplex_t IT_0428 = 2*conj(IT_0410);
    const ccomplex_t IT_0429 = (-2)*conj(IT_0414);
    const ccomplex_t IT_0430 = (-2)*conj(IT_0425);
    const ccomplex_t IT_0431 = IT_0428 + IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = s_23*s_34*IT_0126;
    const ccomplex_t IT_0433 = s_23*s_34*IT_0130;
    const ccomplex_t IT_0434 = (-2)*IT_0433;
    const ccomplex_t IT_0435 = IT_0374 + IT_0434;
    const ccomplex_t IT_0436 = IT_0081*IT_0435;
    const ccomplex_t IT_0437 = (-0.666666666666667)*IT_0432 + (
      -0.166666666666667)*IT_0436;
    const ccomplex_t IT_0438 = s_24 + IT_0437;
    const ccomplex_t IT_0439 = m_N_1*IT_0438;
    const ccomplex_t IT_0440 = (-6)*IT_0439;
    const ccomplex_t IT_0441 = (-2)*IT_0373;
    const ccomplex_t IT_0442 = IT_0372 + IT_0441;
    const ccomplex_t IT_0443 = IT_0126*IT_0442;
    const ccomplex_t IT_0444 = 0.166666666666667*IT_0436 + (-0.666666666666667
      )*IT_0443;
    const ccomplex_t IT_0445 = s_24 + IT_0444;
    const ccomplex_t IT_0446 = m_N_1*IT_0445;
    const ccomplex_t IT_0447 = 6*IT_0446;
    const ccomplex_t IT_0448 = s_12*s_24;
    const ccomplex_t IT_0449 = s_12*s_24*IT_0130;
    const ccomplex_t IT_0450 = s_14*IT_0049*IT_0130;
    const ccomplex_t IT_0451 = s_13*s_23*s_24;
    const ccomplex_t IT_0452 = s_13*s_34*IT_0049;
    const ccomplex_t IT_0453 = -IT_0450 + (-2)*IT_0451 + 2*IT_0452;
    const ccomplex_t IT_0454 = IT_0449 + IT_0453;
    const ccomplex_t IT_0455 = IT_0126*IT_0454;
    const ccomplex_t IT_0456 = (-0.5)*IT_0455;
    const ccomplex_t IT_0457 = s_14*IT_0049;
    const ccomplex_t IT_0458 = -IT_0457;
    const ccomplex_t IT_0459 = IT_0448 + IT_0456 + IT_0458;
    const ccomplex_t IT_0460 = (-2)*IT_0363 + (-2)*IT_0364;
    const ccomplex_t IT_0461 = 4*IT_0460;
    const ccomplex_t IT_0462 = 2*IT_0366 + 2*IT_0367;
    const ccomplex_t IT_0463 = 4*IT_0462;
    const ccomplex_t IT_0464 = IT_0461 + IT_0463;
    const ccomplex_t IT_0465 = 4*conj(IT_0460);
    const ccomplex_t IT_0466 = 4*conj(IT_0462);
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = (-2)*IT_0357 + (-2)*IT_0358;
    const ccomplex_t IT_0469 = 4*IT_0468;
    const ccomplex_t IT_0470 = 2*IT_0354 + 2*IT_0355;
    const ccomplex_t IT_0471 = 4*IT_0470;
    const ccomplex_t IT_0472 = IT_0469 + IT_0471;
    const ccomplex_t IT_0473 = 4*conj(IT_0468);
    const ccomplex_t IT_0474 = 4*conj(IT_0470);
    const ccomplex_t IT_0475 = IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = pow(s_23, 2);
    const ccomplex_t IT_0477 = s_14*IT_0476;
    const ccomplex_t IT_0478 = s_12*s_23*s_34;
    const ccomplex_t IT_0479 = 0.5*IT_0449 + (-0.5)*IT_0450 + -IT_0478;
    const ccomplex_t IT_0480 = IT_0477 + IT_0479;
    const ccomplex_t IT_0481 = IT_0126*IT_0480;
    const ccomplex_t IT_0482 = -IT_0448;
    const ccomplex_t IT_0483 = IT_0457 + IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = s_12*s_14*IT_0130;
    const ccomplex_t IT_0485 = s_13*s_14*s_23;
    const ccomplex_t IT_0486 = s_24*IT_0048*IT_0130;
    const ccomplex_t IT_0487 = s_23*s_34*IT_0048;
    const ccomplex_t IT_0488 = (-2)*IT_0485 + -IT_0486 + 2*IT_0487;
    const ccomplex_t IT_0489 = IT_0484 + IT_0488;
    const ccomplex_t IT_0490 = IT_0126*IT_0489;
    const ccomplex_t IT_0491 = s_24*IT_0048;
    const ccomplex_t IT_0492 = 2*IT_0491;
    const ccomplex_t IT_0493 = s_12*s_14;
    const ccomplex_t IT_0494 = (-2)*IT_0493;
    const ccomplex_t IT_0495 = IT_0490 + IT_0492 + IT_0494;
    const ccomplex_t IT_0496 = 2*IT_0160;
    const ccomplex_t IT_0497 = (-2)*IT_0496;
    const ccomplex_t IT_0498 = (-2)*IT_0181 + -IT_0201 + IT_0204 + 2*IT_0216;
    const ccomplex_t IT_0499 = (-2)*IT_0498;
    const ccomplex_t IT_0500 = IT_0497 + IT_0499;
    const ccomplex_t IT_0501 = (-2)*conj(IT_0496);
    const ccomplex_t IT_0502 = (-2)*conj(IT_0498);
    const ccomplex_t IT_0503 = IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = 2*IT_0292;
    const ccomplex_t IT_0505 = (-2)*IT_0504;
    const ccomplex_t IT_0506 = IT_0200*IT_0307;
    const ccomplex_t IT_0507 = IT_0200*IT_0245;
    const ccomplex_t IT_0508 = (-2)*IT_0303 + 2*IT_0305 + -IT_0506 + IT_0507;
    const ccomplex_t IT_0509 = (-2)*IT_0508;
    const ccomplex_t IT_0510 = IT_0505 + IT_0509;
    const ccomplex_t IT_0511 = (-2)*conj(IT_0504);
    const ccomplex_t IT_0512 = (-2)*conj(IT_0508);
    const ccomplex_t IT_0513 = IT_0511 + IT_0512;
    const ccomplex_t IT_0514 = s_34*IT_0081;
    const ccomplex_t IT_0515 = s_13*IT_0130;
    const ccomplex_t IT_0516 = IT_0514*IT_0515;
    const ccomplex_t IT_0517 = IT_0126*IT_0345;
    const ccomplex_t IT_0518 = (-0.5)*IT_0349 + (-0.5)*IT_0516 + 0.5*IT_0517;
    const ccomplex_t IT_0519 = s_14 + IT_0518;
    const ccomplex_t IT_0520 = m_N_3*IT_0519;
    const ccomplex_t IT_0521 = (-2)*IT_0520;
    const ccomplex_t IT_0522 = 2*IT_0123;
    const ccomplex_t IT_0523 = 2*conj(IT_0123);
    const ccomplex_t IT_0524 = 2*IT_0080;
    const ccomplex_t IT_0525 = 2*conj(IT_0080);
    const ccomplex_t IT_0526 = (conj(IT_0161) + conj(IT_0181) + conj(IT_0217) 
      + conj(IT_0356) + conj(IT_0357) + conj(IT_0358))*IT_0522 + (IT_0161 +
       IT_0181 + IT_0217 + IT_0356 + IT_0357 + IT_0358)*IT_0523 + (conj(IT_0293)
       + conj(IT_0306) + conj(IT_0311) + conj(IT_0365) + conj(IT_0368))*IT_0524 
      + (IT_0293 + IT_0306 + IT_0311 + IT_0365 + IT_0368)*IT_0525;
    const ccomplex_t IT_0527 = s_13*s_24*s_34*IT_0081;
    const ccomplex_t IT_0528 = s_13*s_23*IT_0126;
    const ccomplex_t IT_0529 = s_14*s_24*IT_0126;
    const ccomplex_t IT_0530 = -IT_0528 + -IT_0529;
    const ccomplex_t IT_0531 = s_12 + IT_0527 + IT_0530;
    const ccomplex_t IT_0532 = m_N_3*IT_0531;
    const ccomplex_t IT_0533 = conj(IT_0468) + conj(IT_0470);
    const ccomplex_t IT_0534 = IT_0468 + IT_0470;
    const ccomplex_t IT_0535 = conj(IT_0460) + conj(IT_0462);
    const ccomplex_t IT_0536 = IT_0460 + IT_0462;
    const ccomplex_t IT_0537 = IT_0126*IT_0476;
    const ccomplex_t IT_0538 = pow(s_24, 2);
    const ccomplex_t IT_0539 = IT_0126*IT_0538;
    const ccomplex_t IT_0540 = IT_0537 + IT_0539;
    const ccomplex_t IT_0541 = s_23*s_24*s_34*IT_0081;
    const ccomplex_t IT_0542 = -IT_0049 + -IT_0541;
    const ccomplex_t IT_0543 = IT_0540 + IT_0542;
    const ccomplex_t IT_0544 = m_N_1*IT_0543;
    const ccomplex_t IT_0545 = s_12*IT_0130;
    const ccomplex_t IT_0546 = s_13*s_23;
    const ccomplex_t IT_0547 = (-2)*IT_0546;
    const ccomplex_t IT_0548 = IT_0545 + IT_0547;
    const ccomplex_t IT_0549 = IT_0126*IT_0548;
    const ccomplex_t IT_0550 = (s_12 + -IT_0529)*(s_12 + (-0.5)*IT_0549);
    const ccomplex_t IT_0551 = (-2)*IT_0550;
    const ccomplex_t IT_0552 = 2*IT_0496;
    const ccomplex_t IT_0553 = 2*conj(IT_0496);
    const ccomplex_t IT_0554 = 2*IT_0504;
    const ccomplex_t IT_0555 = 2*conj(IT_0504);
    const ccomplex_t IT_0556 = 2*IT_0498;
    const ccomplex_t IT_0557 = 2*conj(IT_0498);
    const ccomplex_t IT_0558 = 2*IT_0508;
    const ccomplex_t IT_0559 = 2*conj(IT_0508);
    const ccomplex_t IT_0560 = 2*IT_0460;
    const ccomplex_t IT_0561 = 2*conj(IT_0460);
    const ccomplex_t IT_0562 = 2*IT_0468;
    const ccomplex_t IT_0563 = 2*conj(IT_0468);
    const ccomplex_t IT_0564 = IT_0533*IT_0552 + IT_0534*IT_0553 + IT_0535
      *IT_0554 + IT_0536*IT_0555 + conj(IT_0470)*IT_0556 + IT_0470*IT_0557 +
       conj(IT_0462)*IT_0558 + IT_0462*IT_0559 + conj(IT_0508)*IT_0560 + IT_0508
      *IT_0561 + conj(IT_0498)*IT_0562 + IT_0498*IT_0563;
    const ccomplex_t IT_0565 = s_14*s_24;
    const ccomplex_t IT_0566 = -IT_0565;
    const ccomplex_t IT_0567 = IT_0545 + IT_0566;
    const ccomplex_t IT_0568 = IT_0126*IT_0567;
    const ccomplex_t IT_0569 = -IT_0529 + (-0.25)*IT_0568;
    const ccomplex_t IT_0570 = s_12 + IT_0569;
    const ccomplex_t IT_0571 = m_N_1*IT_0570;
    const ccomplex_t IT_0572 = (-4)*IT_0571;
    const ccomplex_t IT_0573 = 2*IT_0266;
    const ccomplex_t IT_0574 = 2*IT_0275;
    const ccomplex_t IT_0575 = IT_0573 + IT_0574;
    const ccomplex_t IT_0576 = 2*conj(IT_0266);
    const ccomplex_t IT_0577 = 2*conj(IT_0275);
    const ccomplex_t IT_0578 = IT_0576 + IT_0577;
    const ccomplex_t IT_0579 = conj(IT_0266) + conj(IT_0275);
    const ccomplex_t IT_0580 = IT_0266 + IT_0275;
    const ccomplex_t IT_0581 = 2*IT_0338;
    const ccomplex_t IT_0582 = 2*conj(IT_0338);
    const ccomplex_t IT_0583 = conj(IT_0324)*IT_0500 + IT_0324*IT_0503 + conj
      (IT_0251)*IT_0510 + IT_0251*IT_0513 + conj(IT_0338)*IT_0552 + IT_0338
      *IT_0553 + conj(IT_0508)*IT_0575 + IT_0508*IT_0578 + IT_0554*IT_0579 +
       IT_0555*IT_0580 + conj(IT_0498)*IT_0581 + IT_0498*IT_0582;
    const ccomplex_t IT_0584 = pow(s_14, 2);
    const ccomplex_t IT_0585 = IT_0048*IT_0130;
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = IT_0584 + IT_0586;
    const ccomplex_t IT_0588 = IT_0126*IT_0587;
    const ccomplex_t IT_0589 = IT_0126*IT_0584;
    const ccomplex_t IT_0590 = -IT_0048;
    const ccomplex_t IT_0591 = IT_0589 + IT_0590;
    const ccomplex_t IT_0592 = (-4)*IT_0591;
    const ccomplex_t IT_0593 = IT_0588 + IT_0592;
    const ccomplex_t IT_0594 = m_N_3*IT_0593;
    const ccomplex_t IT_0595 = conj(IT_0398) + conj(IT_0402);
    const ccomplex_t IT_0596 = IT_0398 + IT_0402;
    const ccomplex_t IT_0597 = conj(IT_0414) + conj(IT_0425);
    const ccomplex_t IT_0598 = IT_0414 + IT_0425;
    const ccomplex_t IT_0599 = conj(IT_0496) + conj(IT_0498);
    const ccomplex_t IT_0600 = IT_0496 + IT_0498;
    const ccomplex_t IT_0601 = IT_0594*(IT_0428*IT_0508 + IT_0411*conj(IT_0508
      ) + conj(IT_0410)*IT_0554 + IT_0410*IT_0555 + IT_0500*IT_0595 + IT_0503
      *IT_0596 + IT_0510*IT_0597 + IT_0513*IT_0598 + IT_0388*IT_0599 + IT_0405
      *IT_0600);
    const ccomplex_t IT_0602 = IT_0081*IT_0129;
    const ccomplex_t IT_0603 = (-10)*s_14;
    const ccomplex_t IT_0604 = IT_0602 + IT_0603;
    const ccomplex_t IT_0605 = m_N_3*IT_0604;
    const ccomplex_t IT_0606 = conj(IT_0181) + conj(IT_0217);
    const ccomplex_t IT_0607 = IT_0181 + IT_0217;
    const ccomplex_t IT_0608 = conj(IT_0306) + conj(IT_0311);
    const ccomplex_t IT_0609 = IT_0306 + IT_0311;
    const ccomplex_t IT_0610 = IT_0124*IT_0591;
    const ccomplex_t IT_0611 = IT_0504 + IT_0508;
    const ccomplex_t IT_0612 = 6*conj(IT_0504);
    const ccomplex_t IT_0613 = -IT_0049;
    const ccomplex_t IT_0614 = IT_0539 + IT_0613;
    const ccomplex_t IT_0615 = IT_0124*IT_0614;
    const ccomplex_t IT_0616 = 6*IT_0460;
    const ccomplex_t IT_0617 = 6*conj(IT_0460);
    const ccomplex_t IT_0618 = 1.125*IT_0124;
    const ccomplex_t IT_0619 = 16*IT_0410;
    const ccomplex_t IT_0620 = 16*conj(IT_0410);
    const ccomplex_t IT_0621 = 16*IT_0387;
    const ccomplex_t IT_0622 = (-16)*IT_0398;
    const ccomplex_t IT_0623 = (-16)*IT_0402;
    const ccomplex_t IT_0624 = 16*conj(IT_0387);
    const ccomplex_t IT_0625 = (-16)*conj(IT_0398);
    const ccomplex_t IT_0626 = (-16)*conj(IT_0402);
    const ccomplex_t IT_0627 = ((-16)*conj(IT_0251) + 16*conj(IT_0266) + 16
      *conj(IT_0275))*IT_0324 + ((-16)*IT_0251 + 16*IT_0266 + 16*IT_0275)*conj
      (IT_0324) + (16*conj(IT_0251) + (-16)*conj(IT_0266) + (-16)*conj(IT_0275))
      *IT_0338 + (16*IT_0251 + (-16)*IT_0266 + (-16)*IT_0275)*conj(IT_0338) + (
      -16)*conj(IT_0387)*IT_0410 + (-16)*IT_0387*conj(IT_0410) + IT_0595*IT_0619
       + IT_0596*IT_0620 + IT_0597*(IT_0621 + IT_0622 + IT_0623) + IT_0598*
      (IT_0624 + IT_0625 + IT_0626);
    const ccomplex_t IT_0628 = IT_0126*IT_0545;
    const ccomplex_t IT_0629 = (-0.5)*IT_0628;
    const ccomplex_t IT_0630 = s_14*s_23;
    const ccomplex_t IT_0631 = s_13*s_24;
    const ccomplex_t IT_0632 = s_12*s_34;
    const ccomplex_t IT_0633 = -IT_0631 + -IT_0632;
    const ccomplex_t IT_0634 = IT_0630 + IT_0633;
    const ccomplex_t IT_0635 = IT_0514*IT_0634;
    const ccomplex_t IT_0636 = (-0.25)*IT_0635;
    const ccomplex_t IT_0637 = s_12 + IT_0629 + IT_0636;
    const ccomplex_t IT_0638 = 8*IT_0123;
    const ccomplex_t IT_0639 = 8*conj(IT_0123);
    const ccomplex_t IT_0640 = 8*IT_0080;
    const ccomplex_t IT_0641 = 8*conj(IT_0080);
    const ccomplex_t IT_0642 = (-8)*IT_0123;
    const ccomplex_t IT_0643 = (-8)*conj(IT_0123);
    const ccomplex_t IT_0644 = (-8)*IT_0080;
    const ccomplex_t IT_0645 = (-8)*conj(IT_0080);
    const ccomplex_t IT_0646 = IT_0630 + IT_0632;
    const ccomplex_t IT_0647 = -IT_0631;
    const ccomplex_t IT_0648 = IT_0646 + IT_0647;
    const ccomplex_t IT_0649 = IT_0514*IT_0648;
    const ccomplex_t IT_0650 = 0.25*IT_0649;
    const ccomplex_t IT_0651 = s_12 + IT_0629 + IT_0650;
    const ccomplex_t IT_0652 = (-2)*IT_0565;
    const ccomplex_t IT_0653 = IT_0545 + IT_0652;
    const ccomplex_t IT_0654 = IT_0126*IT_0653;
    const ccomplex_t IT_0655 = (-0.5)*IT_0654;
    const ccomplex_t IT_0656 = (-0.5)*IT_0549;
    const ccomplex_t IT_0657 = s_12*IT_0128;
    const ccomplex_t IT_0658 = s_13*s_23*IT_0130;
    const ccomplex_t IT_0659 = s_14*s_24*IT_0130;
    const ccomplex_t IT_0660 = s_14*s_23*s_34;
    const ccomplex_t IT_0661 = (-2)*IT_0658 + (-2)*IT_0659 + 4*IT_0660;
    const ccomplex_t IT_0662 = IT_0657 + IT_0661;
    const ccomplex_t IT_0663 = IT_0081*IT_0662;
    const ccomplex_t IT_0664 = 0.25*IT_0663;
    const ccomplex_t IT_0665 = 8*IT_0398;
    const ccomplex_t IT_0666 = (-8)*IT_0387;
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = 8*conj(IT_0398);
    const ccomplex_t IT_0669 = (-8)*conj(IT_0387);
    const ccomplex_t IT_0670 = IT_0668 + IT_0669;
    const ccomplex_t IT_0671 = 8*IT_0387;
    const ccomplex_t IT_0672 = 8*IT_0410;
    const ccomplex_t IT_0673 = 8*IT_0402;
    const ccomplex_t IT_0674 = 8*IT_0425;
    const ccomplex_t IT_0675 = (-8)*IT_0410;
    const ccomplex_t IT_0676 = (-8)*conj(IT_0410);
    const ccomplex_t IT_0677 = 8*conj(IT_0425);
    const ccomplex_t IT_0678 = conj(IT_0398)*IT_0667 + IT_0398*IT_0669 +
       IT_0402*IT_0670 + conj(IT_0387)*IT_0671 + conj(IT_0410)*IT_0672 + conj
      (IT_0402)*(IT_0665 + IT_0666 + IT_0673) + IT_0597*(IT_0674 + IT_0675) +
       IT_0425*IT_0676 + 8*IT_0414*(conj(IT_0414) + 0.125*IT_0676 + 0.125
      *IT_0677);
    const ccomplex_t IT_0679 = pow(m_Z, 6);
    const ccomplex_t IT_0680 = s_12*IT_0679;
    const ccomplex_t IT_0681 = s_13*s_23*IT_0128;
    const ccomplex_t IT_0682 = (-2)*IT_0681;
    const ccomplex_t IT_0683 = IT_0680 + IT_0682;
    const ccomplex_t IT_0684 = IT_0081*IT_0683;
    const ccomplex_t IT_0685 = 0.25*IT_0684;
    const ccomplex_t IT_0686 = IT_0545 + IT_0685;
    const ccomplex_t IT_0687 = (-8)*IT_0293;
    const ccomplex_t IT_0688 = (-8)*conj(IT_0293);
    const ccomplex_t IT_0689 = (-8)*IT_0365;
    const ccomplex_t IT_0690 = (-8)*conj(IT_0365);
    const ccomplex_t IT_0691 = (-8)*IT_0357;
    const ccomplex_t IT_0692 = (-8)*conj(IT_0357);
    const ccomplex_t IT_0693 = conj(IT_0356) + conj(IT_0358);
    const ccomplex_t IT_0694 = (-8)*IT_0161;
    const ccomplex_t IT_0695 = IT_0356 + IT_0358;
    const ccomplex_t IT_0696 = (-8)*conj(IT_0161);
    const ccomplex_t IT_0697 = (-8)*IT_0358;
    const ccomplex_t IT_0698 = (-8)*conj(IT_0358);
    const ccomplex_t IT_0699 = conj(IT_0365)*IT_0687 + (-8)*conj(IT_0368)*
      (IT_0609 + (-0.125)*IT_0687) + IT_0365*IT_0688 + (-8)*IT_0368*(IT_0608 + (
      -0.125)*IT_0688) + IT_0608*IT_0689 + IT_0609*IT_0690 + conj(IT_0161)
      *IT_0691 + IT_0161*IT_0692 + IT_0693*IT_0694 + IT_0695*IT_0696 + (-8)
      *IT_0606*(IT_0356 + (-0.125)*IT_0691 + (-0.125)*IT_0697) + (-8)*IT_0607*
      (conj(IT_0356) + (-0.125)*IT_0692 + (-0.125)*IT_0698);
    const ccomplex_t IT_0700 = -IT_0529;
    const ccomplex_t IT_0701 = s_12 + IT_0700;
    const ccomplex_t IT_0702 = IT_0124*IT_0701;
    const ccomplex_t IT_0703 = (-6)*conj(IT_0504);
    const ccomplex_t IT_0704 = (-6)*IT_0460;
    const ccomplex_t IT_0705 = (-6)*conj(IT_0460);
    const ccomplex_t IT_0706 = pow(s_13, 2);
    const ccomplex_t IT_0707 = s_24*IT_0706;
    const ccomplex_t IT_0708 = s_12*s_13*s_34;
    const ccomplex_t IT_0709 = -IT_0486 + 2*IT_0707 + (-2)*IT_0708;
    const ccomplex_t IT_0710 = IT_0484 + IT_0709;
    const ccomplex_t IT_0711 = IT_0126*IT_0710;
    const ccomplex_t IT_0712 = (-0.5)*IT_0711;
    const ccomplex_t IT_0713 = -IT_0491;
    const ccomplex_t IT_0714 = IT_0493 + IT_0712 + IT_0713;
    const ccomplex_t IT_0715 = (-4)*IT_0496;
    const ccomplex_t IT_0716 = (-4)*IT_0498;
    const ccomplex_t IT_0717 = (-4)*conj(IT_0496);
    const ccomplex_t IT_0718 = (-4)*conj(IT_0498);
    const ccomplex_t IT_0719 = (-4)*IT_0504;
    const ccomplex_t IT_0720 = (-4)*IT_0508;
    const ccomplex_t IT_0721 = (-4)*conj(IT_0504);
    const ccomplex_t IT_0722 = (-4)*conj(IT_0508);
    const ccomplex_t IT_0723 = IT_0049*IT_0130;
    const ccomplex_t IT_0724 = s_23*s_24*s_34*IT_0126;
    const ccomplex_t IT_0725 = 0.5*IT_0538 + (-0.5)*IT_0723 + -IT_0724;
    const ccomplex_t IT_0726 = IT_0476 + IT_0725;
    const ccomplex_t IT_0727 = IT_0126*IT_0726;
    const ccomplex_t IT_0728 = -IT_0539;
    const ccomplex_t IT_0729 = IT_0049 + IT_0727 + IT_0728;
    const ccomplex_t IT_0730 = m_N_1*IT_0729;
    const ccomplex_t IT_0731 = (-4)*IT_0468;
    const ccomplex_t IT_0732 = (-4)*IT_0470;
    const ccomplex_t IT_0733 = (-4)*conj(IT_0468);
    const ccomplex_t IT_0734 = (-4)*conj(IT_0470);
    const ccomplex_t IT_0735 = (-4)*IT_0460;
    const ccomplex_t IT_0736 = (-4)*IT_0462;
    const ccomplex_t IT_0737 = (-4)*conj(IT_0460);
    const ccomplex_t IT_0738 = (-4)*conj(IT_0462);
    const ccomplex_t IT_0739 = conj(IT_0387)*(0.25*IT_0461 + IT_0462) +
       IT_0387*(conj(IT_0462) + 0.25*IT_0465) + conj(IT_0410)*(0.25*IT_0469 +
       IT_0470) + IT_0410*(conj(IT_0470) + 0.25*IT_0473) + 0.25*IT_0597*(IT_0731
       + IT_0732) + 0.25*IT_0598*(IT_0733 + IT_0734) + 0.25*IT_0595*(IT_0735 +
       IT_0736) + 0.25*IT_0596*(IT_0737 + IT_0738);
    const ccomplex_t IT_0740 = s_13*s_14*s_34*IT_0126;
    const ccomplex_t IT_0741 = 0.5*IT_0584 + (-0.5)*IT_0585 + -IT_0740;
    const ccomplex_t IT_0742 = IT_0706 + IT_0741;
    const ccomplex_t IT_0743 = IT_0126*IT_0742;
    const ccomplex_t IT_0744 = -IT_0589;
    const ccomplex_t IT_0745 = IT_0048 + IT_0743 + IT_0744;
    const ccomplex_t IT_0746 = m_N_3*IT_0745;
    const ccomplex_t IT_0747 = (-4)*IT_0266;
    const ccomplex_t IT_0748 = (-4)*IT_0275;
    const ccomplex_t IT_0749 = (-4)*conj(IT_0266);
    const ccomplex_t IT_0750 = (-4)*conj(IT_0275);
    const ccomplex_t IT_0751 = conj(IT_0324)*IT_0504 + IT_0324*conj(IT_0504) +
       (conj(IT_0324) + -conj(IT_0338))*IT_0508 + (IT_0324 + -IT_0338)*conj
      (IT_0508) + 0.25*conj(IT_0338)*IT_0719 + 0.25*IT_0338*IT_0721 + IT_0599*
      (IT_0251 + 0.25*IT_0747 + 0.25*IT_0748) + IT_0600*(conj(IT_0251) + 0.25
      *IT_0749 + 0.25*IT_0750);
    const ccomplex_t IT_0752 = 0.5*IT_0684;
    const ccomplex_t IT_0753 = (-4)*IT_0565;
    const ccomplex_t IT_0754 = IT_0545 + IT_0752 + IT_0753;
    const ccomplex_t IT_0755 = s_13*s_24*s_34;
    const ccomplex_t IT_0756 = (-2)*IT_0658 + (-2)*IT_0659 + 4*IT_0755;
    const ccomplex_t IT_0757 = IT_0657 + IT_0756;
    const ccomplex_t IT_0758 = IT_0126*IT_0757;
    const ccomplex_t IT_0759 = -IT_0758;
    const ccomplex_t IT_0760 = IT_0754 + IT_0759;
    const ccomplex_t IT_0761 = (-4)*IT_0365;
    const ccomplex_t IT_0762 = (-4)*IT_0357;
    const ccomplex_t IT_0763 = (-4)*conj(IT_0357);
    const ccomplex_t IT_0764 = (-4)*IT_0358;
    const ccomplex_t IT_0765 = IT_0126*IT_0662;
    const ccomplex_t IT_0766 = -IT_0765;
    const ccomplex_t IT_0767 = IT_0754 + IT_0766;
    const ccomplex_t IT_0768 = (-4)*IT_0293;
    const ccomplex_t IT_0769 = (-4)*conj(IT_0293);
    const ccomplex_t IT_0770 = (-4)*IT_0161;
    const ccomplex_t IT_0771 = (-4)*conj(IT_0161);
    const ccomplex_t IT_0772 = (-4)*IT_0181;
    const ccomplex_t IT_0773 = (-4)*IT_0306;
    const ccomplex_t IT_0774 = conj(IT_0293)*IT_0768 + IT_0609*IT_0769 + conj
      (IT_0161)*IT_0770 + IT_0607*IT_0771 + (-4)*IT_0606*(IT_0217 + (-0.25)
      *IT_0770 + (-0.25)*IT_0772) + (-4)*IT_0608*(IT_0311 + (-0.25)*IT_0768 + (
      -0.25)*IT_0773);
    const ccomplex_t IT_0775 = s_23*IT_0130;
    const ccomplex_t IT_0776 = IT_0514*IT_0775;
    const ccomplex_t IT_0777 = IT_0126*IT_0372;
    const ccomplex_t IT_0778 = (-0.5)*IT_0443 + (-0.5)*IT_0776 + 0.5*IT_0777;
    const ccomplex_t IT_0779 = s_24 + IT_0778;
    const ccomplex_t IT_0780 = m_N_1*IT_0779;
    const ccomplex_t IT_0781 = (-2)*IT_0780;
    const ccomplex_t IT_0782 = (-2)*IT_0123;
    const ccomplex_t IT_0783 = (-2)*conj(IT_0123);
    const ccomplex_t IT_0784 = (-2)*IT_0080;
    const ccomplex_t IT_0785 = (-2)*conj(IT_0080);
    const ccomplex_t IT_0786 = (conj(IT_0293) + conj(IT_0306) + conj(IT_0311) 
      + conj(IT_0365) + conj(IT_0368))*IT_0782 + (IT_0293 + IT_0306 + IT_0311 +
       IT_0365 + IT_0368)*IT_0783 + (conj(IT_0161) + conj(IT_0181) + conj
      (IT_0217) + conj(IT_0356) + conj(IT_0357) + conj(IT_0358))*IT_0784 + 
      (IT_0161 + IT_0181 + IT_0217 + IT_0356 + IT_0357 + IT_0358)*IT_0785;
    const ccomplex_t IT_0787 = s_14*s_23*s_34*IT_0081;
    const ccomplex_t IT_0788 = s_12 + IT_0530 + IT_0787;
    const ccomplex_t IT_0789 = m_N_1*IT_0788;
    const ccomplex_t IT_0790 = conj(IT_0504) + conj(IT_0508);
    const ccomplex_t IT_0791 = IT_0126*IT_0706;
    const ccomplex_t IT_0792 = IT_0589 + IT_0791;
    const ccomplex_t IT_0793 = s_13*s_14*s_34*IT_0081;
    const ccomplex_t IT_0794 = -IT_0048 + -IT_0793;
    const ccomplex_t IT_0795 = IT_0792 + IT_0794;
    const ccomplex_t IT_0796 = m_N_3*IT_0795;
    const ccomplex_t IT_0797 = -(s_12 + (-0.5)*IT_0549)*(IT_0048 + -IT_0589);
    const ccomplex_t IT_0798 = (-2)*IT_0797;
    const ccomplex_t IT_0799 = -(IT_0049 + -IT_0539)*(s_12 + (-0.5)*IT_0549);
    const ccomplex_t IT_0800 = (-2)*IT_0799;
    const ccomplex_t IT_0801 = (-2)*IT_0460;
    const ccomplex_t IT_0802 = (-2)*conj(IT_0460);
    const ccomplex_t IT_0803 = (-2)*IT_0468;
    const ccomplex_t IT_0804 = (-2)*conj(IT_0468);
    const ccomplex_t IT_0805 = (-2)*IT_0470;
    const ccomplex_t IT_0806 = (-2)*IT_0462;
    const ccomplex_t IT_0807 = s_13*s_24*s_34*IT_0126;
    const ccomplex_t IT_0808 = (-2)*IT_0546 + -IT_0565 + 2*IT_0807;
    const ccomplex_t IT_0809 = IT_0545 + IT_0808;
    const ccomplex_t IT_0810 = IT_0126*IT_0809;
    const ccomplex_t IT_0811 = -IT_0529 + (-0.5)*IT_0810;
    const ccomplex_t IT_0812 = s_12 + IT_0811;
    const ccomplex_t IT_0813 = m_N_3*IT_0812;
    const ccomplex_t IT_0814 = (-2)*IT_0813;
    const ccomplex_t IT_0815 = conj(IT_0338)*IT_0560 + IT_0338*IT_0561 + (-2)
      *conj(IT_0470)*(IT_0251 + (-0.5)*IT_0573 + (-0.5)*IT_0574) + conj(IT_0468)
      *IT_0575 + IT_0534*IT_0578 + (-2)*conj(IT_0462)*(IT_0324 + (-0.5)*IT_0581)
       + IT_0462*IT_0582 + IT_0324*IT_0802 + IT_0251*IT_0804 + conj(IT_0251)*
      (IT_0803 + IT_0805) + conj(IT_0324)*(IT_0801 + IT_0806);
    const ccomplex_t IT_0816 = -IT_0723;
    const ccomplex_t IT_0817 = IT_0538 + IT_0816;
    const ccomplex_t IT_0818 = IT_0126*IT_0817;
    const ccomplex_t IT_0819 = (-4)*IT_0614;
    const ccomplex_t IT_0820 = IT_0818 + IT_0819;
    const ccomplex_t IT_0821 = m_N_1*IT_0820;
    const ccomplex_t IT_0822 = 2*IT_0462;
    const ccomplex_t IT_0823 = IT_0560 + IT_0822;
    const ccomplex_t IT_0824 = 2*conj(IT_0462);
    const ccomplex_t IT_0825 = IT_0561 + IT_0824;
    const ccomplex_t IT_0826 = 2*IT_0470;
    const ccomplex_t IT_0827 = IT_0562 + IT_0826;
    const ccomplex_t IT_0828 = 2*conj(IT_0470);
    const ccomplex_t IT_0829 = IT_0563 + IT_0828;
    const ccomplex_t IT_0830 = IT_0267 + IT_0276;
    const ccomplex_t IT_0831 = IT_0280 + IT_0281;
    const ccomplex_t IT_0832 = IT_0343*IT_0470 + IT_0339*conj(IT_0470) +
       IT_0579*IT_0801 + IT_0580*IT_0802 + conj(IT_0338)*IT_0803 + IT_0338
      *IT_0804 + conj(IT_0251)*IT_0823 + IT_0251*IT_0825 + conj(IT_0324)*IT_0827
       + IT_0324*IT_0829 + conj(IT_0462)*IT_0830 + IT_0462*IT_0831;
    const ccomplex_t IT_0833 = IT_0081*IT_0374;
    const ccomplex_t IT_0834 = (-10)*s_24;
    const ccomplex_t IT_0835 = IT_0833 + IT_0834;
    const ccomplex_t IT_0836 = m_N_1*IT_0835;
    const ccomplex_t IT_0837 = 0.666666666666667*IT_0777 + (-0.166666666666667
      )*IT_0833;
    const ccomplex_t IT_0838 = s_24 + IT_0837;
    const ccomplex_t IT_0839 = m_N_1*IT_0838;
    const ccomplex_t IT_0840 = (-6)*IT_0839;
    const ccomplex_t IT_0841 = s_14*s_23*s_34*IT_0126;
    const ccomplex_t IT_0842 = (-2)*IT_0546 + -IT_0565 + 2*IT_0841;
    const ccomplex_t IT_0843 = IT_0545 + IT_0842;
    const ccomplex_t IT_0844 = IT_0126*IT_0843;
    const ccomplex_t IT_0845 = -IT_0529 + (-0.5)*IT_0844;
    const ccomplex_t IT_0846 = s_12 + IT_0845;
    const ccomplex_t IT_0847 = m_N_1*IT_0846;
    const ccomplex_t IT_0848 = (-2)*IT_0847;
    const ccomplex_t IT_0849 = (-2)*IT_0387;
    const ccomplex_t IT_0850 = (-2)*conj(IT_0387);
    const ccomplex_t IT_0851 = (-2)*IT_0410;
    const ccomplex_t IT_0852 = (-2)*conj(IT_0410);
    const ccomplex_t IT_0853 = conj(IT_0410)*IT_0497 + IT_0410*IT_0501 + conj
      (IT_0387)*IT_0505 + IT_0387*IT_0511 + (IT_0554 + IT_0558)*IT_0595 + 
      (IT_0555 + IT_0559)*IT_0596 + (IT_0552 + IT_0556)*IT_0597 + (IT_0553 +
       IT_0557)*IT_0598 + conj(IT_0508)*IT_0849 + IT_0508*IT_0850 + conj(IT_0498
      )*IT_0851 + IT_0498*IT_0852;
    const ccomplex_t IT_0854 = m_N_3*IT_0570;
    const ccomplex_t IT_0855 = (-4)*IT_0854;
    const ccomplex_t IT_0856 = (conj(IT_0410)*IT_0801 + IT_0410*IT_0802 +
       IT_0597*IT_0823 + IT_0598*IT_0825 + IT_0595*IT_0827 + IT_0596*IT_0829 +
       IT_0533*IT_0849 + IT_0534*IT_0850 + conj(IT_0462)*IT_0851 + IT_0462
      *IT_0852)*IT_0855;
    const ccomplex_t IT_0857 = 0.666666666666667*IT_0517 + (-0.166666666666667
      )*IT_0602;
    const ccomplex_t IT_0858 = s_14 + IT_0857;
    const ccomplex_t IT_0859 = m_N_3*IT_0858;
    const ccomplex_t IT_0860 = (-6)*IT_0859;
    const ccomplex_t IT_0861 = IT_0124*IT_0130;
    const ccomplex_t IT_0862 = 9*IT_0861;
    const ccomplex_t IT_0863 = (-2)*IT_0357;
    const ccomplex_t IT_0864 = (-2)*IT_0358;
    const ccomplex_t IT_0865 = (-2)*IT_0356;
    const ccomplex_t IT_0866 = IT_0863 + IT_0864 + IT_0865;
    const ccomplex_t IT_0867 = (-2)*conj(IT_0357);
    const ccomplex_t IT_0868 = (-2)*conj(IT_0358);
    const ccomplex_t IT_0869 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0870 = IT_0867 + IT_0868 + IT_0869;
    const ccomplex_t IT_0871 = (-2)*IT_0161;
    const ccomplex_t IT_0872 = (-2)*IT_0181;
    const ccomplex_t IT_0873 = (-2)*IT_0217;
    const ccomplex_t IT_0874 = IT_0871 + IT_0872 + IT_0873;
    const ccomplex_t IT_0875 = (-2)*conj(IT_0161);
    const ccomplex_t IT_0876 = (-2)*conj(IT_0181);
    const ccomplex_t IT_0877 = (-2)*conj(IT_0217);
    const ccomplex_t IT_0878 = IT_0875 + IT_0876 + IT_0877;
    const ccomplex_t IT_0879 = (-2)*IT_0293;
    const ccomplex_t IT_0880 = (-2)*conj(IT_0293);
    const ccomplex_t IT_0881 = (-2)*IT_0365;
    const ccomplex_t IT_0882 = (-2)*conj(IT_0365);
    const ccomplex_t IT_0883 = m_N_1*m_N_3*IT_0126;
    const ccomplex_t IT_0884 = (-0.25)*IT_0128;
    const ccomplex_t IT_0885 = IT_0082 + IT_0884;
    const ccomplex_t IT_0886 = IT_0883*IT_0885;
    const ccomplex_t IT_0887 = (-4)*IT_0886;
    const ccomplex_t IT_0888 = (-4)*IT_0125;
    const ccomplex_t IT_0889 = -IT_0123;
    const ccomplex_t IT_0890 = -conj(IT_0123);
    const ccomplex_t IT_0891 = -IT_0080;
    const ccomplex_t IT_0892 = -conj(IT_0080);
    const ccomplex_t IT_0893 = conj(IT_0123)*(IT_0338 + IT_0410) + IT_0123*
      (conj(IT_0338) + conj(IT_0410)) + IT_0080*(conj(IT_0387) + IT_0579) + conj
      (IT_0080)*(IT_0387 + IT_0580) + (conj(IT_0324) + conj(IT_0414) + conj
      (IT_0425))*IT_0889 + (IT_0324 + IT_0414 + IT_0425)*IT_0890 + (conj(IT_0251
      ) + conj(IT_0398) + conj(IT_0402))*IT_0891 + (IT_0251 + IT_0398 + IT_0402)
      *IT_0892;
    const ccomplex_t IT_0894 = IT_0081*IT_0757;
    const ccomplex_t IT_0895 = 0.25*IT_0894;
    const ccomplex_t IT_0896 = 8*IT_0266;
    const ccomplex_t IT_0897 = 8*conj(IT_0266);
    const ccomplex_t IT_0898 = 8*IT_0275;
    const ccomplex_t IT_0899 = 0.125*conj(IT_0275);
    const ccomplex_t IT_0900 = (-8)*conj(IT_0266);
    const ccomplex_t IT_0901 = (-8)*conj(IT_0275);
    const ccomplex_t IT_0902 = IT_0900 + IT_0901;
    const ccomplex_t IT_0903 = 0.125*IT_0251;
    const ccomplex_t IT_0904 = 8*IT_0251;
    const ccomplex_t IT_0905 = (-8)*IT_0266;
    const ccomplex_t IT_0906 = (-8)*IT_0275;
    const ccomplex_t IT_0907 = IT_0904 + IT_0905 + IT_0906;
    const ccomplex_t IT_0908 = 0.125*conj(IT_0251);
    const ccomplex_t IT_0909 = (-8)*conj(IT_0338);
    const ccomplex_t IT_0910 = 0.125*IT_0324;
    const ccomplex_t IT_0911 = 8*IT_0324;
    const ccomplex_t IT_0912 = (-8)*IT_0338;
    const ccomplex_t IT_0913 = IT_0911 + IT_0912;
    const ccomplex_t IT_0914 = 0.125*conj(IT_0324);
    const ccomplex_t IT_0915 = 8*(s_12 + IT_0655 + IT_0656 + IT_0895)*(IT_0338
      *conj(IT_0338) + 0.125*IT_0579*IT_0896 + 0.125*IT_0275*IT_0897 + IT_0898
      *IT_0899 + IT_0902*IT_0903 + IT_0907*IT_0908 + IT_0909*IT_0910 + IT_0913
      *IT_0914);
    const ccomplex_t IT_0916 = (-2)*IT_0082;
    const ccomplex_t IT_0917 = IT_0128 + IT_0916;
    const ccomplex_t IT_0918 = IT_0081*IT_0917;
    const ccomplex_t IT_0919 = 0.25*IT_0918;
    const ccomplex_t IT_0920 = 1 + IT_0919;
    const ccomplex_t IT_0921 = IT_0124*IT_0920;
    const ccomplex_t IT_0922 = 8*conj(IT_0387);
    const ccomplex_t IT_0923 = 0.125*conj(IT_0338);
    const ccomplex_t IT_0924 = 8*conj(IT_0251);
    const ccomplex_t IT_0925 = IT_0900 + IT_0901 + IT_0924;
    const ccomplex_t IT_0926 = 0.125*IT_0598;
    const ccomplex_t IT_0927 = 8*IT_0324*(conj(IT_0402) + 0.125*IT_0668 +
       0.125*IT_0669) + 8*IT_0410*(conj(IT_0275) + 0.125*IT_0897) + conj(IT_0410
      )*(IT_0896 + IT_0898) + 8*IT_0676*IT_0903 + IT_0597*IT_0907 + 8*IT_0675
      *IT_0908 + IT_0596*IT_0909 + IT_0595*IT_0912 + 8*(IT_0665 + IT_0666 +
       IT_0673)*IT_0914 + IT_0338*IT_0922 + 8*IT_0671*IT_0923 + 8*IT_0925*IT_0926;
    const ccomplex_t IT_0928 = s_12*IT_0082;
    const ccomplex_t IT_0929 = (-2)*IT_0658 + (-2)*IT_0659 + 2*IT_0660 + 2
      *IT_0755 + (-2)*IT_0928;
    const ccomplex_t IT_0930 = IT_0657 + IT_0929;
    const ccomplex_t IT_0931 = IT_0081*IT_0930;
    const ccomplex_t IT_0932 = (-0.25)*IT_0931;
    const ccomplex_t IT_0933 = -IT_0528;
    const ccomplex_t IT_0934 = IT_0338*conj(IT_0410) + conj(IT_0324)*(IT_0414 
      + 0.125*IT_0674 + 0.125*IT_0675) + 0.125*conj(IT_0387)*IT_0896 + 0.125
      *IT_0387*IT_0897 + IT_0671*IT_0899 + 0.125*IT_0398*IT_0902 + IT_0670
      *IT_0903 + 0.125*conj(IT_0398)*(IT_0905 + IT_0906) + 0.125*conj(IT_0402)
      *IT_0907 + IT_0667*IT_0908 + (IT_0676 + IT_0677)*IT_0910 + 0.125*conj
      (IT_0425)*IT_0912 + 0.125*conj(IT_0414)*IT_0913 + 0.125*IT_0275*IT_0922 +
       IT_0672*IT_0923 + 0.125*IT_0402*IT_0925 + IT_0909*IT_0926;
    const ccomplex_t IT_0935 = 4*IT_0123*(IT_0086*conj(IT_0123) + -conj
      (IT_0080)*IT_0125) + 4*IT_0080*(conj(IT_0080)*IT_0086 + -conj(IT_0123)
      *IT_0125) + IT_0138*(IT_0218*IT_0277 + IT_0278*IT_0282 + IT_0312*IT_0340 +
       IT_0341*IT_0344) + IT_0353*(IT_0277*IT_0360 + IT_0282*IT_0362 + IT_0340
      *IT_0369 + IT_0344*IT_0370) + (IT_0277*IT_0312 + IT_0218*IT_0340 + IT_0282
      *IT_0341 + IT_0278*IT_0344)*IT_0379 + IT_0385*(IT_0360*IT_0404 + IT_0362
      *IT_0408 + IT_0369*IT_0427 + IT_0370*IT_0431) + (IT_0369*IT_0404 + IT_0370
      *IT_0408 + IT_0360*IT_0427 + IT_0362*IT_0431)*IT_0440 + (IT_0312*IT_0404 +
       IT_0341*IT_0408 + IT_0218*IT_0427 + IT_0278*IT_0431)*IT_0447 + IT_0459*
      (IT_0369*IT_0464 + IT_0370*IT_0467 + IT_0360*IT_0472 + IT_0362*IT_0475) + 
      (IT_0312*IT_0464 + IT_0341*IT_0467 + IT_0218*IT_0472 + IT_0278*IT_0475)
      *IT_0483 + IT_0495*(IT_0218*IT_0500 + IT_0278*IT_0503 + IT_0312*IT_0510 +
       IT_0341*IT_0513) + IT_0521*IT_0526 + IT_0532*(IT_0522*IT_0533 + IT_0523
      *IT_0534 + IT_0524*IT_0535 + IT_0525*IT_0536) + (IT_0524*IT_0533 + IT_0525
      *IT_0534 + IT_0522*IT_0535 + IT_0523*IT_0536)*IT_0544 + IT_0551*IT_0564 +
       IT_0572*IT_0583 + IT_0601 + IT_0605*(conj(IT_0161)*IT_0404 + IT_0161
      *IT_0408 + conj(IT_0293)*IT_0427 + IT_0293*IT_0431 + IT_0388*IT_0606 +
       IT_0405*IT_0607 + IT_0411*IT_0608 + IT_0428*IT_0609) + 6*IT_0610*(IT_0599
      *IT_0611 + IT_0600*(conj(IT_0508) + 0.166666666666667*IT_0612)) + 6
      *IT_0615*(IT_0533*(IT_0462 + 0.166666666666667*IT_0616) + IT_0534*(conj
      (IT_0462) + 0.166666666666667*IT_0617)) + IT_0618*IT_0627 + IT_0637*
      (IT_0579*IT_0638 + IT_0580*IT_0639 + conj(IT_0338)*IT_0640 + IT_0338
      *IT_0641 + conj(IT_0251)*IT_0642 + IT_0251*IT_0643 + conj(IT_0324)*IT_0644
       + IT_0324*IT_0645) + (conj(IT_0387)*IT_0638 + IT_0387*IT_0639 + conj
      (IT_0410)*IT_0640 + IT_0410*IT_0641 + IT_0595*IT_0642 + IT_0596*IT_0643 +
       IT_0597*IT_0644 + IT_0598*IT_0645)*IT_0651 + (s_12 + IT_0655 + IT_0656 +
       IT_0664)*IT_0678 + IT_0686*IT_0699 + (-6)*IT_0702*(IT_0533*IT_0611 +
       IT_0534*(conj(IT_0508) + (-0.166666666666667)*IT_0703) + IT_0599*(IT_0462
       + (-0.166666666666667)*IT_0704) + IT_0600*(conj(IT_0462) + (
      -0.166666666666667)*IT_0705)) + IT_0714*(IT_0360*(IT_0715 + IT_0716) +
       IT_0362*(IT_0717 + IT_0718) + IT_0369*(IT_0719 + IT_0720) + IT_0370*
      (IT_0721 + IT_0722)) + 4*IT_0730*IT_0739 + 4*IT_0746*IT_0751 + (-4)
      *IT_0760*(IT_0369*(IT_0368 + (-0.25)*IT_0761) + (-0.25)*conj(IT_0357)
      *IT_0762 + (-0.25)*IT_0695*IT_0763 + IT_0693*(IT_0356 + (-0.25)*IT_0762 + 
      (-0.25)*IT_0764)) + IT_0767*IT_0774 + IT_0781*IT_0786 + IT_0789*(IT_0611
      *IT_0783 + IT_0599*IT_0784 + IT_0600*IT_0785 + IT_0782*IT_0790) + (IT_0599
      *IT_0782 + IT_0600*IT_0783 + IT_0611*IT_0785 + IT_0784*IT_0790)*IT_0796 + 
      (conj(IT_0498)*IT_0499 + IT_0498*IT_0501 + conj(IT_0508)*IT_0509 + IT_0508
      *IT_0511 + IT_0497*IT_0599 + IT_0505*IT_0790)*IT_0798 + IT_0800*(IT_0535
      *IT_0801 + IT_0462*IT_0802 + IT_0533*IT_0803 + IT_0470*IT_0804 + conj
      (IT_0470)*IT_0805 + conj(IT_0462)*IT_0806) + IT_0814*IT_0815 + IT_0821
      *IT_0832 + (IT_0343*IT_0356 + IT_0339*conj(IT_0356) + IT_0325*IT_0360 +
       IT_0342*IT_0362 + IT_0279*IT_0365 + IT_0252*conj(IT_0365) + IT_0282
      *IT_0368 + IT_0277*conj(IT_0368))*IT_0836 + (IT_0339*IT_0359 + IT_0343
      *IT_0361 + conj(IT_0365)*IT_0830 + IT_0365*IT_0831)*IT_0840 + IT_0848
      *IT_0853 + IT_0856 + ((IT_0399 + IT_0403)*IT_0606 + (IT_0406 + IT_0407)
      *IT_0607 + (IT_0415 + IT_0426)*IT_0608 + (IT_0429 + IT_0430)*IT_0609)
      *IT_0860 + IT_0862*(conj(IT_0368)*IT_0866 + IT_0368*IT_0870 + IT_0608
      *IT_0874 + IT_0609*IT_0878 + IT_0218*IT_0879 + IT_0278*IT_0880 + IT_0360
      *IT_0881 + IT_0362*IT_0882) + (IT_0608*IT_0866 + IT_0609*IT_0870 + conj
      (IT_0368)*IT_0874 + IT_0368*IT_0878 + IT_0360*IT_0879 + IT_0362*IT_0880 +
       IT_0218*IT_0881 + IT_0278*IT_0882)*IT_0887 + IT_0888*IT_0893 + IT_0915 +
       IT_0921*IT_0927 + 8*(s_12 + IT_0700 + IT_0932 + IT_0933)*IT_0934;
    return create_ccomplex_return(IT_0935);
}

