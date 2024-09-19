#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_3_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_3_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0005 + 0.5*IT_0009 + 0.5*IT_0014 + (
      -0.5)*IT_0017 + (-0.5)*IT_0019 + 0.5*IT_0021 + 0.5*IT_0023 + (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cpow(IT_0000, -1);
    const ccomplex_t IT_0029 = m_s*e_em*IT_0006*IT_0011*IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_N_1, 2);
    const ccomplex_t IT_0032 = pow(m_N_3, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + -IT_0031 + -IT_0032 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = IT_0003*IT_0008;
    const ccomplex_t IT_0038 = IT_0004*IT_0007;
    const ccomplex_t IT_0039 = IT_0012*IT_0016;
    const ccomplex_t IT_0040 = IT_0013*IT_0015;
    const ccomplex_t IT_0041 = IT_0003*IT_0020;
    const ccomplex_t IT_0042 = IT_0007*IT_0018;
    const ccomplex_t IT_0043 = IT_0012*IT_0024;
    const ccomplex_t IT_0044 = IT_0015*IT_0022;
    const ccomplex_t IT_0045 = (-0.5)*IT_0037 + (-0.5)*IT_0038 + 0.5*IT_0039 +
       0.5*IT_0040 + (-0.5)*IT_0041 + (-0.5)*IT_0042 + 0.5*IT_0043 + 0.5*IT_0044;
    const ccomplex_t IT_0046 = m_s*e_em*IT_0011*IT_0027;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0031 + -IT_0032 + IT_0048 
      + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = -IT_0036 + -IT_0052;
    const ccomplex_t IT_0054 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0055 = IT_0003*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0057 = IT_0007*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0059 = IT_0012*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0069 = IT_0015*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0055 + (-0.5)*IT_0057 + (-0.5)*IT_0059 +
       0.5*IT_0061 + 0.5*IT_0063 + (-0.5)*IT_0065 + (-0.5)*IT_0067 + 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0035*IT_0070;
    const ccomplex_t IT_0072 = IT_0003*IT_0056;
    const ccomplex_t IT_0073 = IT_0007*IT_0054;
    const ccomplex_t IT_0074 = IT_0012*IT_0060;
    const ccomplex_t IT_0075 = IT_0015*IT_0058;
    const ccomplex_t IT_0076 = IT_0003*IT_0064;
    const ccomplex_t IT_0077 = IT_0007*IT_0062;
    const ccomplex_t IT_0078 = IT_0012*IT_0068;
    const ccomplex_t IT_0079 = IT_0015*IT_0066;
    const ccomplex_t IT_0080 = 0.5*IT_0072 + 0.5*IT_0073 + (-0.5)*IT_0074 + (
      -0.5)*IT_0075 + 0.5*IT_0076 + 0.5*IT_0077 + (-0.5)*IT_0078 + (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0051*IT_0080;
    const ccomplex_t IT_0082 = -IT_0071 + -IT_0081;
    const ccomplex_t IT_0083 = pow(m_s, 2);
    const ccomplex_t IT_0084 = s_34 + IT_0083;
    const ccomplex_t IT_0085 = s_12*IT_0084;
    const ccomplex_t IT_0086 = cos(alpha);
    const ccomplex_t IT_0087 = IT_0002*IT_0086;
    const ccomplex_t IT_0088 = IT_0004*IT_0087;
    const ccomplex_t IT_0089 = sin(alpha);
    const ccomplex_t IT_0090 = IT_0002*IT_0089;
    const ccomplex_t IT_0091 = IT_0008*IT_0090;
    const ccomplex_t IT_0092 = IT_0018*IT_0087;
    const ccomplex_t IT_0093 = IT_0020*IT_0090;
    const ccomplex_t IT_0094 = IT_0088 + IT_0091 + IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = IT_0011*IT_0086;
    const ccomplex_t IT_0096 = IT_0013*IT_0095;
    const ccomplex_t IT_0097 = IT_0011*IT_0089;
    const ccomplex_t IT_0098 = IT_0016*IT_0097;
    const ccomplex_t IT_0099 = IT_0022*IT_0095;
    const ccomplex_t IT_0100 = IT_0024*IT_0097;
    const ccomplex_t IT_0101 = -IT_0096 + -IT_0098 + -IT_0099 + -IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_s*e_em*IT_0011*IT_0027
      *IT_0028*IT_0089;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0031 + -IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0104*IT_0109;
    const ccomplex_t IT_0111 = IT_0008*IT_0087;
    const ccomplex_t IT_0112 = IT_0013*IT_0097;
    const ccomplex_t IT_0113 = IT_0020*IT_0087;
    const ccomplex_t IT_0114 = IT_0022*IT_0097;
    const ccomplex_t IT_0115 = IT_0111 + IT_0112 + IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = IT_0004*IT_0090;
    const ccomplex_t IT_0117 = IT_0016*IT_0095;
    const ccomplex_t IT_0118 = IT_0018*IT_0090;
    const ccomplex_t IT_0119 = IT_0024*IT_0095;
    const ccomplex_t IT_0120 = -IT_0116 + -IT_0117 + -IT_0118 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0115 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_s*e_em*IT_0011*IT_0027
      *IT_0028*IT_0086;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_12 + -IT_0031 + -IT_0032 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0123*IT_0128;
    const ccomplex_t IT_0130 = -IT_0110 + -IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0031 + IT_0083 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0132 = N_B1*e_em;
    const ccomplex_t IT_0133 = IT_0002*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = N_W1*e_em;
    const ccomplex_t IT_0136 = IT_0011*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + (-3)*IT_0137);
    const ccomplex_t IT_0139 = 0.166666666666667*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0011*IT_0027*IT_0028;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = IT_0131*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0031 + IT_0083 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0011*IT_0027*IT_0028;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0002;
    const ccomplex_t IT_0149 = (-0.333333333333333)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147*IT_0149;
    const ccomplex_t IT_0151 = IT_0145*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0144 + IT_0152;
    const ccomplex_t IT_0154 = s_12*IT_0083;
    const ccomplex_t IT_0155 = s_14*s_23;
    const ccomplex_t IT_0156 = s_13*s_24;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = s_12*s_34;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0154 + IT_0155 + IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0002;
    const ccomplex_t IT_0162 = (-0.333333333333333)*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0140;
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_23 + IT_0032 + IT_0083 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0146;
    const ccomplex_t IT_0169 = N_B3*e_em;
    const ccomplex_t IT_0170 = IT_0002*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = N_W3*e_em;
    const ccomplex_t IT_0173 = IT_0011*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + (-3)*IT_0174);
    const ccomplex_t IT_0176 = 0.166666666666667*IT_0175;
    const ccomplex_t IT_0177 = IT_0168*IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_23 + IT_0032 + IT_0083 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0177*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = -IT_0167 + -IT_0180;
    const ccomplex_t IT_0182 = -IT_0155;
    const ccomplex_t IT_0183 = IT_0154 + IT_0156 + IT_0159 + IT_0182;
    const ccomplex_t IT_0184 = -s_34;
    const ccomplex_t IT_0185 = IT_0083 + IT_0184;
    const ccomplex_t IT_0186 = s_12*IT_0185;
    const ccomplex_t IT_0187 = IT_0054*IT_0087;
    const ccomplex_t IT_0188 = IT_0056*IT_0090;
    const ccomplex_t IT_0189 = IT_0062*IT_0087;
    const ccomplex_t IT_0190 = IT_0064*IT_0090;
    const ccomplex_t IT_0191 = IT_0187 + IT_0188 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = IT_0058*IT_0095;
    const ccomplex_t IT_0193 = IT_0060*IT_0097;
    const ccomplex_t IT_0194 = IT_0066*IT_0095;
    const ccomplex_t IT_0195 = IT_0068*IT_0097;
    const ccomplex_t IT_0196 = -IT_0192 + -IT_0193 + -IT_0194 + -IT_0195;
    const ccomplex_t IT_0197 = IT_0191 + IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0109*IT_0199;
    const ccomplex_t IT_0201 = IT_0056*IT_0087;
    const ccomplex_t IT_0202 = IT_0058*IT_0097;
    const ccomplex_t IT_0203 = IT_0064*IT_0087;
    const ccomplex_t IT_0204 = IT_0066*IT_0097;
    const ccomplex_t IT_0205 = IT_0201 + IT_0202 + IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = IT_0054*IT_0090;
    const ccomplex_t IT_0207 = IT_0060*IT_0095;
    const ccomplex_t IT_0208 = IT_0062*IT_0090;
    const ccomplex_t IT_0209 = IT_0068*IT_0095;
    const ccomplex_t IT_0210 = -IT_0206 + -IT_0207 + -IT_0208 + -IT_0209;
    const ccomplex_t IT_0211 = IT_0205 + IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = IT_0128*IT_0213;
    const ccomplex_t IT_0215 = -IT_0200 + -IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0011*IT_0027*IT_0028;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = conj(N_B3)*e_em;
    const ccomplex_t IT_0219 = IT_0002*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = conj(N_W3)*e_em;
    const ccomplex_t IT_0222 = IT_0011*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*(IT_0220 + (-3)*IT_0223);
    const ccomplex_t IT_0225 = 0.166666666666667*IT_0224;
    const ccomplex_t IT_0226 = IT_0217*IT_0225;
    const ccomplex_t IT_0227 = IT_0131*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0002;
    const ccomplex_t IT_0230 = (-0.333333333333333)*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0011*IT_0027*IT_0028;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = IT_0230*IT_0232;
    const ccomplex_t IT_0234 = IT_0145*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0228 + IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0216;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0002;
    const ccomplex_t IT_0239 = (-0.333333333333333)*IT_0238;
    const ccomplex_t IT_0240 = IT_0237*IT_0239;
    const ccomplex_t IT_0241 = IT_0165*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = conj(N_B1)*e_em;
    const ccomplex_t IT_0244 = IT_0002*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W1)*e_em;
    const ccomplex_t IT_0247 = IT_0011*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*(IT_0245 + (-3)*IT_0248);
    const ccomplex_t IT_0250 = 0.166666666666667*IT_0249;
    const ccomplex_t IT_0251 = (-0.5)*IT_0231;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = IT_0178*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = -IT_0242 + -IT_0254;
    const ccomplex_t IT_0256 = m_N_1*m_N_3;
    const ccomplex_t IT_0257 = -IT_0083;
    const ccomplex_t IT_0258 = s_34 + IT_0257;
    const ccomplex_t IT_0259 = IT_0256*IT_0258;
    const ccomplex_t IT_0260 = (-8)*IT_0259;
    const ccomplex_t IT_0261 = 3*IT_0130;
    const ccomplex_t IT_0262 = 3*conj(IT_0130);
    const ccomplex_t IT_0263 = m_s*m_N_3;
    const ccomplex_t IT_0264 = -s_14;
    const ccomplex_t IT_0265 = s_13 + IT_0264;
    const ccomplex_t IT_0266 = IT_0263*IT_0265;
    const ccomplex_t IT_0267 = 2*IT_0266;
    const ccomplex_t IT_0268 = IT_0139*IT_0225;
    const ccomplex_t IT_0269 = IT_0131*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = IT_0147*IT_0232;
    const ccomplex_t IT_0272 = IT_0145*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = -IT_0270 + -IT_0273;
    const ccomplex_t IT_0275 = IT_0141*IT_0217;
    const ccomplex_t IT_0276 = IT_0131*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = IT_0149*IT_0230;
    const ccomplex_t IT_0279 = IT_0145*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = -IT_0277 + -IT_0280;
    const ccomplex_t IT_0282 = 3*IT_0215;
    const ccomplex_t IT_0283 = 3*conj(IT_0215);
    const ccomplex_t IT_0284 = m_s*m_N_1;
    const ccomplex_t IT_0285 = -s_24;
    const ccomplex_t IT_0286 = s_23 + IT_0285;
    const ccomplex_t IT_0287 = IT_0284*IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0155 + IT_0156 + IT_0159;
    const ccomplex_t IT_0290 = 6*IT_0153;
    const ccomplex_t IT_0291 = 6*conj(IT_0153);
    const ccomplex_t IT_0292 = m_s*s_13*m_N_3;
    const ccomplex_t IT_0293 = 6*IT_0181;
    const ccomplex_t IT_0294 = 12*IT_0153;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = 6*conj(IT_0181);
    const ccomplex_t IT_0297 = 12*conj(IT_0153);
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = IT_0002*IT_0010;
    const ccomplex_t IT_0300 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0301 = IT_0299*IT_0300;
    const ccomplex_t IT_0302 = IT_0001*IT_0011;
    const ccomplex_t IT_0303 = IT_0300*IT_0302;
    const ccomplex_t IT_0304 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0305 = IT_0299*IT_0304;
    const ccomplex_t IT_0306 = IT_0302*IT_0304;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*(IT_0301 + IT_0303 + 
      -IT_0305 + -IT_0306);
    const ccomplex_t IT_0308 = (-0.5)*IT_0307;
    const ccomplex_t IT_0309 = e_em*IT_0299;
    const ccomplex_t IT_0310 = e_em*IT_0302;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*(IT_0309 + 3*IT_0310);
    const ccomplex_t IT_0312 = (-0.166666666666667)*IT_0311;
    const ccomplex_t IT_0313 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0031 + -IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0314 = IT_0312*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0308*IT_0315;
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = 12*IT_0181;
    const ccomplex_t IT_0319 = 24*IT_0153;
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = 12*conj(IT_0181);
    const ccomplex_t IT_0322 = 24*conj(IT_0153);
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0325 = IT_0299*IT_0324;
    const ccomplex_t IT_0326 = IT_0302*IT_0324;
    const ccomplex_t IT_0327 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0328 = IT_0299*IT_0327;
    const ccomplex_t IT_0329 = IT_0302*IT_0327;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*(IT_0325 + IT_0326 + 
      -IT_0328 + -IT_0329);
    const ccomplex_t IT_0331 = (-0.5)*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0010;
    const ccomplex_t IT_0333 = 0.333333333333333*IT_0332;
    const ccomplex_t IT_0334 = IT_0313*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0331*IT_0335;
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = 12*IT_0255;
    const ccomplex_t IT_0339 = 24*IT_0236;
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = 12*conj(IT_0255);
    const ccomplex_t IT_0342 = 24*conj(IT_0236);
    const ccomplex_t IT_0343 = IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = 6*IT_0281;
    const ccomplex_t IT_0345 = 6*conj(IT_0281);
    const ccomplex_t IT_0346 = 12*IT_0281;
    const ccomplex_t IT_0347 = 12*conj(IT_0281);
    const ccomplex_t IT_0348 = IT_0292*(conj(IT_0274)*IT_0295 + IT_0274
      *IT_0298 + conj(IT_0317)*IT_0320 + IT_0317*IT_0323 + conj(IT_0337)*IT_0340
       + IT_0337*IT_0343 + conj(IT_0255)*IT_0344 + IT_0255*IT_0345 + conj
      (IT_0236)*IT_0346 + IT_0236*IT_0347);
    const ccomplex_t IT_0349 = m_s*s_14*m_N_3;
    const ccomplex_t IT_0350 = IT_0162*IT_0239;
    const ccomplex_t IT_0351 = IT_0165*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*IT_0351;
    const ccomplex_t IT_0353 = IT_0168*IT_0251;
    const ccomplex_t IT_0354 = IT_0178*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*IT_0354;
    const ccomplex_t IT_0356 = IT_0352 + IT_0355;
    const ccomplex_t IT_0357 = IT_0290 + IT_0318;
    const ccomplex_t IT_0358 = IT_0291 + IT_0321;
    const ccomplex_t IT_0359 = IT_0308*IT_0335;
    const ccomplex_t IT_0360 = -IT_0359;
    const ccomplex_t IT_0361 = 24*IT_0181;
    const ccomplex_t IT_0362 = IT_0294 + IT_0361;
    const ccomplex_t IT_0363 = 24*conj(IT_0181);
    const ccomplex_t IT_0364 = IT_0297 + IT_0363;
    const ccomplex_t IT_0365 = IT_0315*IT_0331;
    const ccomplex_t IT_0366 = -IT_0365;
    const ccomplex_t IT_0367 = 12*IT_0236;
    const ccomplex_t IT_0368 = 24*IT_0255;
    const ccomplex_t IT_0369 = IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = 12*conj(IT_0236);
    const ccomplex_t IT_0371 = 24*conj(IT_0255);
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = IT_0163*IT_0237;
    const ccomplex_t IT_0374 = IT_0165*IT_0373;
    const ccomplex_t IT_0375 = (0 + _Complex_I*1)*IT_0374;
    const ccomplex_t IT_0376 = IT_0176*IT_0250;
    const ccomplex_t IT_0377 = IT_0178*IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*IT_0377;
    const ccomplex_t IT_0379 = IT_0375 + IT_0378;
    const ccomplex_t IT_0380 = 6*IT_0379;
    const ccomplex_t IT_0381 = 6*conj(IT_0379);
    const ccomplex_t IT_0382 = 12*IT_0379;
    const ccomplex_t IT_0383 = 12*conj(IT_0379);
    const ccomplex_t IT_0384 = IT_0349*(conj(IT_0356)*IT_0357 + IT_0356
      *IT_0358 + conj(IT_0360)*IT_0362 + IT_0360*IT_0364 + conj(IT_0366)*IT_0369
       + IT_0366*IT_0372 + conj(IT_0236)*IT_0380 + IT_0236*IT_0381 + conj
      (IT_0255)*IT_0382 + IT_0255*IT_0383);
    const ccomplex_t IT_0385 = m_s*s_24*m_N_1;
    const ccomplex_t IT_0386 = 6*IT_0274;
    const ccomplex_t IT_0387 = 6*conj(IT_0274);
    const ccomplex_t IT_0388 = 12*IT_0274;
    const ccomplex_t IT_0389 = 12*conj(IT_0274);
    const ccomplex_t IT_0390 = IT_0385*(conj(IT_0281)*IT_0295 + IT_0281
      *IT_0298 + IT_0323*IT_0337 + IT_0320*conj(IT_0337) + conj(IT_0317)*IT_0340
       + IT_0317*IT_0343 + conj(IT_0255)*IT_0386 + IT_0255*IT_0387 + conj
      (IT_0236)*IT_0388 + IT_0236*IT_0389);
    const ccomplex_t IT_0391 = 24*IT_0317;
    const ccomplex_t IT_0392 = IT_0388 + IT_0391;
    const ccomplex_t IT_0393 = 24*conj(IT_0317);
    const ccomplex_t IT_0394 = IT_0389 + IT_0393;
    const ccomplex_t IT_0395 = 12*IT_0356;
    const ccomplex_t IT_0396 = 12*conj(IT_0356);
    const ccomplex_t IT_0397 = IT_0347*IT_0366 + IT_0346*conj(IT_0366) +
       IT_0345*IT_0379 + IT_0344*conj(IT_0379) + conj(IT_0337)*(24*IT_0366 +
       IT_0382) + IT_0337*(24*conj(IT_0366) + IT_0383) + conj(IT_0356)*IT_0386 +
       IT_0356*IT_0387 + conj(IT_0360)*IT_0392 + IT_0360*IT_0394 + conj(IT_0317)
      *IT_0395 + IT_0317*IT_0396;
    const ccomplex_t IT_0398 = s_34*m_N_1*m_N_3;
    const ccomplex_t IT_0399 = IT_0345*IT_0356 + IT_0344*conj(IT_0356) + (24
      *conj(IT_0337) + IT_0347)*IT_0360 + (24*IT_0337 + IT_0346)*conj(IT_0360) +
       conj(IT_0274)*IT_0380 + IT_0274*IT_0381 + conj(IT_0317)*IT_0382 + IT_0317
      *IT_0383 + conj(IT_0366)*IT_0392 + IT_0366*IT_0394 + conj(IT_0337)*IT_0395
       + IT_0337*IT_0396;
    const ccomplex_t IT_0400 = m_s*s_23*m_N_1;
    const ccomplex_t IT_0401 = conj(IT_0236)*IT_0356 + IT_0236*conj(IT_0356) +
       0.166666666666667*IT_0364*IT_0366 + 0.166666666666667*IT_0362*conj
      (IT_0366) + 0.166666666666667*conj(IT_0360)*IT_0369 + 0.166666666666667
      *IT_0360*IT_0372 + 0.166666666666667*IT_0358*IT_0379 + 0.166666666666667
      *IT_0357*conj(IT_0379) + 0.166666666666667*conj(IT_0255)*IT_0395 +
       0.166666666666667*IT_0255*IT_0396;
    const ccomplex_t IT_0402 = 24*conj(IT_0281);
    const ccomplex_t IT_0403 = 48*conj(IT_0337);
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 24*IT_0281;
    const ccomplex_t IT_0406 = 24*IT_0274;
    const ccomplex_t IT_0407 = 24*conj(IT_0274);
    const ccomplex_t IT_0408 = 24*IT_0356;
    const ccomplex_t IT_0409 = 48*IT_0360;
    const ccomplex_t IT_0410 = IT_0408 + IT_0409;
    const ccomplex_t IT_0411 = 24*IT_0379;
    const ccomplex_t IT_0412 = 24*conj(IT_0379);
    const ccomplex_t IT_0413 = 24*conj(IT_0356);
    const ccomplex_t IT_0414 = m_N_1*m_N_3*IT_0083;
    const ccomplex_t IT_0415 = conj(IT_0236)*IT_0295 + IT_0236*IT_0298 + conj
      (IT_0274)*IT_0346 + IT_0274*IT_0347 + conj(IT_0255)*IT_0357 + IT_0255
      *IT_0358 + conj(IT_0356)*IT_0382 + IT_0356*IT_0383 + IT_0317*IT_0404 +
       conj(IT_0317)*(48*IT_0337 + IT_0405) + conj(IT_0337)*IT_0406 + IT_0337
      *IT_0407 + conj(IT_0366)*IT_0410 + conj(IT_0360)*(48*IT_0366 + IT_0411) +
       IT_0360*IT_0412 + IT_0366*IT_0413;
    const ccomplex_t IT_0416 = s_13 + s_14;
    const ccomplex_t IT_0417 = IT_0263*IT_0416;
    const ccomplex_t IT_0418 = 6*IT_0053;
    const ccomplex_t IT_0419 = 6*conj(IT_0053);
    const ccomplex_t IT_0420 = 12*conj(IT_0053);
    const ccomplex_t IT_0421 = 12*IT_0082;
    const ccomplex_t IT_0422 = 12*conj(IT_0082);
    const ccomplex_t IT_0423 = (-12)*conj(IT_0053);
    const ccomplex_t IT_0424 = (-12)*IT_0082;
    const ccomplex_t IT_0425 = (-12)*conj(IT_0082);
    const ccomplex_t IT_0426 = (-6)*IT_0053;
    const ccomplex_t IT_0427 = (-6)*conj(IT_0053);
    const ccomplex_t IT_0428 = IT_0053*((-12)*conj(IT_0317) + 12*conj(IT_0360)
      ) + (-6)*conj(IT_0082)*((-0.166666666666667)*IT_0344 + IT_0379) + (-6)
      *IT_0082*((-0.166666666666667)*IT_0345 + conj(IT_0379)) + conj(IT_0356)
      *IT_0418 + IT_0356*IT_0419 + IT_0360*IT_0420 + conj(IT_0337)*IT_0421 +
       IT_0337*IT_0422 + IT_0317*IT_0423 + conj(IT_0366)*IT_0424 + IT_0366
      *IT_0425 + conj(IT_0274)*IT_0426 + IT_0274*IT_0427;
    const ccomplex_t IT_0429 = IT_0154 + IT_0156 + IT_0158 + IT_0182;
    const ccomplex_t IT_0430 = 6*IT_0082;
    const ccomplex_t IT_0431 = 6*conj(IT_0082);
    const ccomplex_t IT_0432 = IT_0084*IT_0256;
    const ccomplex_t IT_0433 = (conj(IT_0082)*((-24)*IT_0053 + (-6)*IT_0181 +
       IT_0290) + IT_0082*((-24)*conj(IT_0053) + (-6)*conj(IT_0181) + IT_0291) +
       conj(IT_0255)*IT_0418 + IT_0255*IT_0419 + conj(IT_0236)*IT_0426 + IT_0236
      *IT_0427)*IT_0432;
    const ccomplex_t IT_0434 = s_23 + s_24;
    const ccomplex_t IT_0435 = IT_0284*IT_0434;
    const ccomplex_t IT_0436 = (-6)*IT_0082;
    const ccomplex_t IT_0437 = (-6)*conj(IT_0082);
    const ccomplex_t IT_0438 = IT_0053*((-12)*conj(IT_0337) + 12*conj(IT_0366)
      ) + conj(IT_0379)*IT_0418 + IT_0379*IT_0419 + IT_0366*IT_0420 + conj
      (IT_0317)*IT_0421 + IT_0317*IT_0422 + IT_0337*IT_0423 + conj(IT_0360)
      *IT_0424 + IT_0360*IT_0425 + conj(IT_0281)*IT_0426 + IT_0281*IT_0427 +
       conj(IT_0274)*IT_0430 + IT_0274*IT_0431 + conj(IT_0356)*IT_0436 + IT_0356
      *IT_0437;
    const ccomplex_t IT_0439 = IT_0154 + IT_0155 + IT_0157 + IT_0158;
    const ccomplex_t IT_0440 = (-2)*IT_0259;
    const ccomplex_t IT_0441 = (-3)*IT_0130;
    const ccomplex_t IT_0442 = (-3)*conj(IT_0130);
    const ccomplex_t IT_0443 = (-3)*IT_0215;
    const ccomplex_t IT_0444 = (-3)*conj(IT_0215);
    const ccomplex_t IT_0445 = (-2)*IT_0266;
    const ccomplex_t IT_0446 = 2*IT_0287;
    const ccomplex_t IT_0447 = 4*IT_0287;
    const ccomplex_t IT_0448 = conj(IT_0337) + conj(IT_0366);
    const ccomplex_t IT_0449 = IT_0337 + IT_0366;
    const ccomplex_t IT_0450 = conj(IT_0317) + conj(IT_0360);
    const ccomplex_t IT_0451 = IT_0317 + IT_0360;
    const ccomplex_t IT_0452 = (-4)*IT_0266;
    const ccomplex_t IT_0453 = (24*IT_0053*conj(IT_0053) + 24*IT_0082*conj
      (IT_0082))*IT_0085 + (-6)*conj(IT_0130)*(IT_0153*IT_0160 + -IT_0181
      *IT_0183) + (-6)*IT_0130*(conj(IT_0153)*IT_0160 + -conj(IT_0181)*IT_0183 +
       4*conj(IT_0130)*IT_0186) + (-6)*conj(IT_0215)*(IT_0160*IT_0236 + -IT_0183
      *IT_0255 + (-0.166666666666667)*IT_0260*IT_0261) + IT_0215*((-24)*IT_0186
      *conj(IT_0215) + (-6)*IT_0160*conj(IT_0236) + 6*IT_0183*conj(IT_0255) +
       IT_0260*IT_0262) + IT_0267*(IT_0262*IT_0274 + IT_0261*conj(IT_0274) +
       conj(IT_0281)*IT_0282 + IT_0281*IT_0283) + (IT_0262*IT_0281 + IT_0261
      *conj(IT_0281) + conj(IT_0274)*IT_0282 + IT_0274*IT_0283)*IT_0288 + 6
      *IT_0289*(conj(IT_0236)*IT_0255 + IT_0236*conj(IT_0255) +
       0.166666666666667*conj(IT_0181)*IT_0290 + 0.166666666666667*IT_0181
      *IT_0291) + IT_0348 + IT_0384 + IT_0390 + IT_0154*IT_0397 + IT_0398
      *IT_0399 + 6*IT_0400*IT_0401 + IT_0156*(conj(IT_0153)*IT_0294 + conj
      (IT_0281)*IT_0346 + conj(IT_0236)*IT_0367 + conj(IT_0274)*IT_0388 +
       IT_0337*IT_0404 + conj(IT_0337)*IT_0405 + conj(IT_0317)*IT_0406 + IT_0317
      *(48*conj(IT_0317) + IT_0407)) + IT_0155*(conj(IT_0181)*IT_0318 + conj
      (IT_0255)*IT_0338 + conj(IT_0379)*IT_0382 + conj(IT_0356)*IT_0395 + conj
      (IT_0360)*IT_0410 + conj(IT_0366)*IT_0411 + IT_0366*(48*conj(IT_0366) +
       IT_0412) + IT_0360*IT_0413) + IT_0414*IT_0415 + IT_0417*IT_0428 + IT_0429
      *(conj(IT_0153)*IT_0426 + IT_0153*IT_0427 + conj(IT_0236)*IT_0430 +
       IT_0236*IT_0431) + IT_0433 + IT_0435*IT_0438 + (conj(IT_0181)*IT_0418 +
       IT_0181*IT_0419 + conj(IT_0255)*IT_0436 + IT_0255*IT_0437)*IT_0439 +
       IT_0440*(conj(IT_0236)*IT_0261 + IT_0236*IT_0262 + conj(IT_0153)*IT_0282 
      + IT_0153*IT_0283 + conj(IT_0255)*IT_0441 + IT_0255*IT_0442 + conj(IT_0181
      )*IT_0443 + IT_0181*IT_0444) + (conj(IT_0356)*IT_0441 + IT_0356*IT_0442 +
       conj(IT_0379)*IT_0443 + IT_0379*IT_0444)*IT_0445 + (conj(IT_0379)*IT_0441
       + IT_0379*IT_0442 + conj(IT_0356)*IT_0443 + IT_0356*IT_0444)*IT_0446 +
       IT_0447*(IT_0441*IT_0448 + IT_0442*IT_0449 + IT_0443*IT_0450 + IT_0444
      *IT_0451) + (IT_0443*IT_0448 + IT_0444*IT_0449 + IT_0441*IT_0450 + IT_0442
      *IT_0451)*IT_0452;
    return create_ccomplex_return(IT_0453);
}

