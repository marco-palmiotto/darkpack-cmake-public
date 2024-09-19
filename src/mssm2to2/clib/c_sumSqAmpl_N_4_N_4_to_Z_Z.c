#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_Z_Z(
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_Z, -4);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, 2);
    const ccomplex_t IT_0003 = s_23*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = pow(m_Z, -2);
    const ccomplex_t IT_0006 = s_24*IT_0002;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = s_23*s_34;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = IT_0006 + IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0004 + 0.5*IT_0007 + (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = s_24 + IT_0012;
    const ccomplex_t IT_0014 = m_N_4*IT_0013;
    const ccomplex_t IT_0015 = (-2)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_N_4, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0002 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0021, -1);
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0028;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + IT_0027 + 
      -IT_0029 + -IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0034 = IT_0022*IT_0033;
    const ccomplex_t IT_0035 = IT_0026*IT_0033;
    const ccomplex_t IT_0036 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0037 = IT_0022*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0036;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0034 + IT_0035 + 
      -IT_0037 + -IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = IT_0018*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0002 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0047 = IT_0022*IT_0046;
    const ccomplex_t IT_0048 = IT_0026*IT_0046;
    const ccomplex_t IT_0049 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = IT_0026*IT_0049;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0047 + IT_0048 + 
      -IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0057 = IT_0026*IT_0056;
    const ccomplex_t IT_0058 = IT_0022*IT_0056;
    const ccomplex_t IT_0059 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0060 = IT_0026*IT_0059;
    const ccomplex_t IT_0061 = IT_0022*IT_0059;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0057 + IT_0058 + 
      -IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0055*IT_0066;
    const ccomplex_t IT_0068 = IT_0045*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0002 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0072 = IT_0022*IT_0071;
    const ccomplex_t IT_0073 = IT_0026*IT_0071;
    const ccomplex_t IT_0074 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0075 = IT_0026*IT_0074;
    const ccomplex_t IT_0076 = IT_0022*IT_0074;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0072 + IT_0073 + 
      -IT_0075 + -IT_0076);
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0083 = IT_0022*IT_0082;
    const ccomplex_t IT_0084 = IT_0026*IT_0082;
    const ccomplex_t IT_0085 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0086 = IT_0026*IT_0085;
    const ccomplex_t IT_0087 = IT_0022*IT_0085;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0083 + IT_0084 + 
      -IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0081*IT_0090;
    const ccomplex_t IT_0092 = IT_0070*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_23 + IT_0002 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0042*IT_0094;
    const ccomplex_t IT_0096 = -IT_0068 + -IT_0092 + -IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0002 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0100 = IT_0022*IT_0099;
    const ccomplex_t IT_0101 = IT_0026*IT_0099;
    const ccomplex_t IT_0102 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0103 = IT_0026*IT_0102;
    const ccomplex_t IT_0104 = IT_0022*IT_0102;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0100 + IT_0101 + 
      -IT_0103 + -IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0108;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = IT_0098*IT_0111;
    const ccomplex_t IT_0113 = cpow(IT_0110, 2);
    const ccomplex_t IT_0114 = IT_0098*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0002 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0067*IT_0116;
    const ccomplex_t IT_0118 = cpow((-2)*s_13 + IT_0002 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0091*IT_0119;
    const ccomplex_t IT_0121 = 0.5*IT_0112 + (-0.5)*IT_0114 + IT_0117 + IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_23 + IT_0002 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow(IT_0109, 2);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0111 + IT_0125;
    const ccomplex_t IT_0127 = IT_0123*IT_0126;
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = cos(alpha);
    const ccomplex_t IT_0130 = IT_0020*IT_0129;
    const ccomplex_t IT_0131 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = sin(alpha);
    const ccomplex_t IT_0134 = IT_0020*IT_0133;
    const ccomplex_t IT_0135 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = IT_0025*IT_0129;
    const ccomplex_t IT_0138 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = IT_0025*IT_0133;
    const ccomplex_t IT_0141 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0132 + -IT_0136 + 
      -IT_0139 + IT_0142);
    const ccomplex_t IT_0144 = cos(beta);
    const ccomplex_t IT_0145 = cpow(IT_0019, -4);
    const ccomplex_t IT_0146 = IT_0021*IT_0129*IT_0144*IT_0145;
    const ccomplex_t IT_0147 = sin(beta);
    const ccomplex_t IT_0148 = IT_0021*IT_0133*IT_0145*IT_0147;
    const ccomplex_t IT_0149 = cpow(IT_0021, -3);
    const ccomplex_t IT_0150 = IT_0129*IT_0144*IT_0149;
    const ccomplex_t IT_0151 = IT_0133*IT_0147*IT_0149;
    const ccomplex_t IT_0152 = cpow(IT_0019, -2);
    const ccomplex_t IT_0153 = IT_0025*IT_0129*IT_0144*IT_0152;
    const ccomplex_t IT_0154 = IT_0025*IT_0133*IT_0147*IT_0152;
    const ccomplex_t IT_0155 = cpow(IT_0021, -2);
    const ccomplex_t IT_0156 = cpow(IT_0152 + IT_0155, -1);
    const ccomplex_t IT_0157 = (IT_0146 + IT_0148 + IT_0150 + IT_0151 + 2
      *IT_0153 + 2*IT_0154)*IT_0156;
    const ccomplex_t IT_0158 = m_W*e_em;
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0143*IT_0163;
    const ccomplex_t IT_0165 = IT_0130*IT_0135;
    const ccomplex_t IT_0166 = IT_0131*IT_0134;
    const ccomplex_t IT_0167 = IT_0137*IT_0141;
    const ccomplex_t IT_0168 = IT_0138*IT_0140;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0165 + IT_0166 + 
      -IT_0167 + -IT_0168);
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = IT_0021*IT_0129*IT_0145*IT_0147;
    const ccomplex_t IT_0172 = IT_0021*IT_0133*IT_0144*IT_0145;
    const ccomplex_t IT_0173 = IT_0129*IT_0147*IT_0149;
    const ccomplex_t IT_0174 = IT_0133*IT_0144*IT_0149;
    const ccomplex_t IT_0175 = IT_0025*IT_0129*IT_0147*IT_0152;
    const ccomplex_t IT_0176 = IT_0025*IT_0133*IT_0144*IT_0152;
    const ccomplex_t IT_0177 = IT_0156*(IT_0171 + -IT_0172 + IT_0173 + 
      -IT_0174 + 2*IT_0175 + (-2)*IT_0176);
    const ccomplex_t IT_0178 = IT_0158*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = IT_0170*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0164 + 0.5*IT_0183;
    const ccomplex_t IT_0185 = 2*IT_0184;
    const ccomplex_t IT_0186 = (-0.5)*IT_0164 + (-0.5)*IT_0183;
    const ccomplex_t IT_0187 = (-2)*IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = 2*conj(IT_0184);
    const ccomplex_t IT_0190 = (-2)*conj(IT_0186);
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = IT_0043 + IT_0120;
    const ccomplex_t IT_0193 = IT_0098*IT_0124;
    const ccomplex_t IT_0194 = 0.5*IT_0112 + IT_0117 + (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = -IT_0113;
    const ccomplex_t IT_0196 = IT_0111 + IT_0195;
    const ccomplex_t IT_0197 = IT_0123*IT_0196;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0200 = IT_0130*IT_0199;
    const ccomplex_t IT_0201 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0202 = IT_0134*IT_0201;
    const ccomplex_t IT_0203 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0204 = IT_0137*IT_0203;
    const ccomplex_t IT_0205 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0206 = IT_0140*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0200 + -IT_0202 + 
      -IT_0204 + IT_0206);
    const ccomplex_t IT_0208 = IT_0163*IT_0207;
    const ccomplex_t IT_0209 = IT_0130*IT_0201;
    const ccomplex_t IT_0210 = IT_0134*IT_0199;
    const ccomplex_t IT_0211 = IT_0137*IT_0205;
    const ccomplex_t IT_0212 = IT_0140*IT_0203;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + IT_0210 + 
      -IT_0211 + -IT_0212);
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0182*IT_0214;
    const ccomplex_t IT_0216 = 0.5*IT_0208 + 0.5*IT_0215;
    const ccomplex_t IT_0217 = 2*IT_0216;
    const ccomplex_t IT_0218 = (-0.5)*IT_0208 + (-0.5)*IT_0215;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = 2*conj(IT_0216);
    const ccomplex_t IT_0222 = (-2)*conj(IT_0218);
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = IT_0015*((conj(IT_0043) + conj(IT_0096) + conj
      (IT_0121) + conj(IT_0128))*IT_0188 + (IT_0043 + IT_0096 + IT_0121 +
       IT_0128)*IT_0191 + (conj(IT_0096) + conj(IT_0192) + conj(IT_0194) + conj
      (IT_0198))*IT_0220 + (IT_0096 + IT_0192 + IT_0194 + IT_0198)*IT_0223);
    const ccomplex_t IT_0225 = 2*IT_0068 + 2*IT_0092 + 2*IT_0095;
    const ccomplex_t IT_0226 = pow(s_14, 2);
    const ccomplex_t IT_0227 = IT_0005*IT_0226;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = IT_0016 + IT_0228;
    const ccomplex_t IT_0230 = IT_0016*IT_0229;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = (-2)*IT_0043 + -IT_0112 + IT_0114 + (-2)
      *IT_0117 + (-2)*IT_0120;
    const ccomplex_t IT_0233 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = s_12 + IT_0234;
    const ccomplex_t IT_0236 = IT_0016*IT_0235;
    const ccomplex_t IT_0237 = s_12*s_14;
    const ccomplex_t IT_0238 = s_12*s_14*IT_0002;
    const ccomplex_t IT_0239 = s_24*IT_0002*IT_0016;
    const ccomplex_t IT_0240 = s_13*s_14*s_23;
    const ccomplex_t IT_0241 = s_23*s_34*IT_0016;
    const ccomplex_t IT_0242 = -IT_0239 + (-2)*IT_0240 + 2*IT_0241;
    const ccomplex_t IT_0243 = IT_0238 + IT_0242;
    const ccomplex_t IT_0244 = IT_0005*IT_0243;
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = s_24*IT_0016;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0237 + IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = conj(IT_0096) + conj(IT_0198);
    const ccomplex_t IT_0250 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0251 = IT_0016 + IT_0250;
    const ccomplex_t IT_0252 = pow(s_13, 2);
    const ccomplex_t IT_0253 = IT_0005*IT_0252;
    const ccomplex_t IT_0254 = -IT_0227 + -IT_0253;
    const ccomplex_t IT_0255 = IT_0251 + IT_0254;
    const ccomplex_t IT_0256 = m_N_4*IT_0255;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0259 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0260 = -IT_0233 + -IT_0259;
    const ccomplex_t IT_0261 = s_12 + IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = m_N_4*IT_0261;
    const ccomplex_t IT_0263 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0264 = IT_0016 + IT_0263;
    const ccomplex_t IT_0265 = pow(s_23, 2);
    const ccomplex_t IT_0266 = IT_0005*IT_0265;
    const ccomplex_t IT_0267 = pow(s_24, 2);
    const ccomplex_t IT_0268 = IT_0005*IT_0267;
    const ccomplex_t IT_0269 = -IT_0266 + -IT_0268;
    const ccomplex_t IT_0270 = IT_0264 + IT_0269;
    const ccomplex_t IT_0271 = m_N_4*IT_0270;
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0186;
    const ccomplex_t IT_0274 = (-2)*IT_0184;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 2*conj(IT_0186);
    const ccomplex_t IT_0277 = (-2)*conj(IT_0184);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = (-2)*IT_0043 + -IT_0112 + (-2)*IT_0117 + (-2)
      *IT_0120 + IT_0193;
    const ccomplex_t IT_0280 = 2*IT_0218;
    const ccomplex_t IT_0281 = (-2)*IT_0216;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 2*conj(IT_0218);
    const ccomplex_t IT_0284 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0287 = s_12 + IT_0260 + IT_0286;
    const ccomplex_t IT_0288 = m_N_4*IT_0287;
    const ccomplex_t IT_0289 = s_13*IT_0002;
    const ccomplex_t IT_0290 = IT_0001*IT_0289;
    const ccomplex_t IT_0291 = s_14*IT_0002;
    const ccomplex_t IT_0292 = IT_0005*IT_0291;
    const ccomplex_t IT_0293 = s_13*s_34;
    const ccomplex_t IT_0294 = (-2)*IT_0293;
    const ccomplex_t IT_0295 = IT_0291 + IT_0294;
    const ccomplex_t IT_0296 = IT_0005*IT_0295;
    const ccomplex_t IT_0297 = (-0.5)*IT_0290 + 0.5*IT_0292 + (-0.5)*IT_0296;
    const ccomplex_t IT_0298 = s_14 + IT_0297;
    const ccomplex_t IT_0299 = m_N_4*IT_0298;
    const ccomplex_t IT_0300 = (-2)*IT_0299;
    const ccomplex_t IT_0301 = ((conj(IT_0096) + conj(IT_0192) + conj(IT_0194)
       + conj(IT_0198))*IT_0275 + (IT_0096 + IT_0192 + IT_0194 + IT_0198)
      *IT_0278 + (conj(IT_0043) + conj(IT_0096) + conj(IT_0121) + conj(IT_0128))
      *IT_0282 + (IT_0043 + IT_0096 + IT_0121 + IT_0128)*IT_0285)*IT_0300;
    const ccomplex_t IT_0302 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0303 = pow(m_Z, 4);
    const ccomplex_t IT_0304 = s_24*IT_0303;
    const ccomplex_t IT_0305 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0306 = (-2)*IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = IT_0000*IT_0307;
    const ccomplex_t IT_0309 = (-0.666666666666667)*IT_0302 + (
      -0.166666666666667)*IT_0308;
    const ccomplex_t IT_0310 = s_24 + IT_0309;
    const ccomplex_t IT_0311 = m_N_4*IT_0310;
    const ccomplex_t IT_0312 = (-6)*IT_0311;
    const ccomplex_t IT_0313 = conj(IT_0043) + conj(IT_0121);
    const ccomplex_t IT_0314 = m_N_4*IT_0111;
    const ccomplex_t IT_0315 = m_N_4*IT_0124;
    const ccomplex_t IT_0316 = IT_0123*(IT_0314 + -IT_0315);
    const ccomplex_t IT_0317 = 2*IT_0316;
    const ccomplex_t IT_0318 = cpow(IT_0066, 2);
    const ccomplex_t IT_0319 = m_N_1*IT_0318;
    const ccomplex_t IT_0320 = m_N_4*IT_0067;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0319 + IT_0321;
    const ccomplex_t IT_0323 = IT_0045*IT_0322;
    const ccomplex_t IT_0324 = cpow(IT_0081, 2);
    const ccomplex_t IT_0325 = m_N_2*IT_0324;
    const ccomplex_t IT_0326 = m_N_4*IT_0091;
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = IT_0325 + IT_0327;
    const ccomplex_t IT_0329 = IT_0070*IT_0328;
    const ccomplex_t IT_0330 = cpow(IT_0041, 2);
    const ccomplex_t IT_0331 = m_N_3*IT_0330;
    const ccomplex_t IT_0332 = m_N_4*IT_0042;
    const ccomplex_t IT_0333 = -IT_0332;
    const ccomplex_t IT_0334 = IT_0331 + IT_0333;
    const ccomplex_t IT_0335 = IT_0094*IT_0334;
    const ccomplex_t IT_0336 = IT_0323 + IT_0329 + IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = IT_0317 + IT_0337;
    const ccomplex_t IT_0339 = IT_0043 + IT_0121;
    const ccomplex_t IT_0340 = 2*conj(IT_0316);
    const ccomplex_t IT_0341 = (-2)*conj(IT_0336);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = conj(IT_0192) + conj(IT_0194);
    const ccomplex_t IT_0344 = m_N_4*IT_0113;
    const ccomplex_t IT_0345 = IT_0123*(IT_0314 + -IT_0344);
    const ccomplex_t IT_0346 = 2*IT_0345;
    const ccomplex_t IT_0347 = cpow(IT_0032, 2);
    const ccomplex_t IT_0348 = m_N_3*IT_0347;
    const ccomplex_t IT_0349 = -IT_0094*(IT_0332 + -IT_0348);
    const ccomplex_t IT_0350 = (-2)*IT_0349;
    const ccomplex_t IT_0351 = cpow(IT_0090, 2);
    const ccomplex_t IT_0352 = m_N_2*IT_0351;
    const ccomplex_t IT_0353 = IT_0070*(IT_0326 + -IT_0352);
    const ccomplex_t IT_0354 = -IT_0353;
    const ccomplex_t IT_0355 = (-2)*IT_0354;
    const ccomplex_t IT_0356 = cpow(IT_0055, 2);
    const ccomplex_t IT_0357 = m_N_1*IT_0356;
    const ccomplex_t IT_0358 = -IT_0045*(IT_0320 + -IT_0357);
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = -IT_0359;
    const ccomplex_t IT_0361 = (-2)*IT_0360;
    const ccomplex_t IT_0362 = IT_0346 + IT_0350 + IT_0355 + IT_0361;
    const ccomplex_t IT_0363 = IT_0192 + IT_0194;
    const ccomplex_t IT_0364 = 2*conj(IT_0345);
    const ccomplex_t IT_0365 = (-2)*conj(IT_0349);
    const ccomplex_t IT_0366 = (-2)*conj(IT_0354);
    const ccomplex_t IT_0367 = (-2)*conj(IT_0360);
    const ccomplex_t IT_0368 = IT_0364 + IT_0365 + IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = (-0.666666666666667)*IT_0011 +
       0.166666666666667*IT_0308;
    const ccomplex_t IT_0370 = s_24 + IT_0369;
    const ccomplex_t IT_0371 = m_N_4*IT_0370;
    const ccomplex_t IT_0372 = 6*IT_0371;
    const ccomplex_t IT_0373 = conj(IT_0096) + conj(IT_0128);
    const ccomplex_t IT_0374 = IT_0096 + IT_0128;
    const ccomplex_t IT_0375 = IT_0096 + IT_0198;
    const ccomplex_t IT_0376 = m_N_4*IT_0005;
    const ccomplex_t IT_0377 = s_14*IT_0303;
    const ccomplex_t IT_0378 = IT_0005*IT_0377;
    const ccomplex_t IT_0379 = (-2)*IT_0293 + (-0.5)*IT_0378;
    const ccomplex_t IT_0380 = IT_0291 + IT_0379;
    const ccomplex_t IT_0381 = IT_0376*IT_0380;
    const ccomplex_t IT_0382 = (-2)*IT_0381;
    const ccomplex_t IT_0383 = IT_0002*IT_0016;
    const ccomplex_t IT_0384 = 9*IT_0383;
    const ccomplex_t IT_0385 = (-2)*IT_0096;
    const ccomplex_t IT_0386 = (-0.5)*IT_0385;
    const ccomplex_t IT_0387 = IT_0198 + IT_0386;
    const ccomplex_t IT_0388 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0389 = (-0.5)*IT_0388;
    const ccomplex_t IT_0390 = conj(IT_0198) + IT_0389;
    const ccomplex_t IT_0391 = (-2)*IT_0192;
    const ccomplex_t IT_0392 = (-2)*IT_0194;
    const ccomplex_t IT_0393 = IT_0391 + IT_0392;
    const ccomplex_t IT_0394 = (-2)*conj(IT_0192);
    const ccomplex_t IT_0395 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = IT_0005*IT_0016;
    const ccomplex_t IT_0398 = pow(s_34, 2);
    const ccomplex_t IT_0399 = (-0.25)*IT_0303;
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = IT_0397*IT_0400;
    const ccomplex_t IT_0402 = (-4)*IT_0401;
    const ccomplex_t IT_0403 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0404 = (-2)*IT_0403;
    const ccomplex_t IT_0405 = IT_0377 + IT_0404;
    const ccomplex_t IT_0406 = IT_0000*IT_0405;
    const ccomplex_t IT_0407 = (-0.666666666666667)*IT_0296 +
       0.166666666666667*IT_0406;
    const ccomplex_t IT_0408 = s_14 + IT_0407;
    const ccomplex_t IT_0409 = m_N_4*IT_0408;
    const ccomplex_t IT_0410 = 6*IT_0409;
    const ccomplex_t IT_0411 = -IT_0018*(IT_0332 + -IT_0348);
    const ccomplex_t IT_0412 = 2*IT_0411;
    const ccomplex_t IT_0413 = IT_0119*(IT_0326 + -IT_0352);
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = 2*IT_0414;
    const ccomplex_t IT_0416 = -IT_0116*(IT_0320 + -IT_0357);
    const ccomplex_t IT_0417 = 2*IT_0416;
    const ccomplex_t IT_0418 = IT_0098*(IT_0314 + -IT_0344);
    const ccomplex_t IT_0419 = -IT_0418;
    const ccomplex_t IT_0420 = 2*IT_0419;
    const ccomplex_t IT_0421 = IT_0412 + IT_0415 + IT_0417 + IT_0420;
    const ccomplex_t IT_0422 = 2*conj(IT_0411);
    const ccomplex_t IT_0423 = 2*conj(IT_0414);
    const ccomplex_t IT_0424 = 2*conj(IT_0416);
    const ccomplex_t IT_0425 = 2*conj(IT_0419);
    const ccomplex_t IT_0426 = IT_0422 + IT_0423 + IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = IT_0018*IT_0334;
    const ccomplex_t IT_0428 = -IT_0315;
    const ccomplex_t IT_0429 = IT_0314 + IT_0428;
    const ccomplex_t IT_0430 = IT_0098*IT_0429;
    const ccomplex_t IT_0431 = IT_0116*IT_0322;
    const ccomplex_t IT_0432 = IT_0119*IT_0328;
    const ccomplex_t IT_0433 = IT_0427 + -IT_0430 + IT_0431 + IT_0432;
    const ccomplex_t IT_0434 = 2*IT_0433;
    const ccomplex_t IT_0435 = 2*conj(IT_0433);
    const ccomplex_t IT_0436 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0437 = (-0.166666666666667)*IT_0406 + (
      -0.666666666666667)*IT_0436;
    const ccomplex_t IT_0438 = s_14 + IT_0437;
    const ccomplex_t IT_0439 = m_N_4*IT_0438;
    const ccomplex_t IT_0440 = (-6)*IT_0439;
    const ccomplex_t IT_0441 = IT_0000*IT_0304;
    const ccomplex_t IT_0442 = (-10)*s_24;
    const ccomplex_t IT_0443 = IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = m_N_4*IT_0443;
    const ccomplex_t IT_0445 = IT_0005*IT_0304;
    const ccomplex_t IT_0446 = (-2)*IT_0008 + (-0.5)*IT_0445;
    const ccomplex_t IT_0447 = IT_0006 + IT_0446;
    const ccomplex_t IT_0448 = IT_0376*IT_0447;
    const ccomplex_t IT_0449 = (-2)*IT_0448;
    const ccomplex_t IT_0450 = s_12*IT_0002;
    const ccomplex_t IT_0451 = s_13*s_23;
    const ccomplex_t IT_0452 = (-2)*IT_0451;
    const ccomplex_t IT_0453 = IT_0450 + IT_0452;
    const ccomplex_t IT_0454 = IT_0005*IT_0453;
    const ccomplex_t IT_0455 = (s_12 + -IT_0233)*(s_12 + (-0.5)*IT_0454);
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = 2*IT_0225;
    const ccomplex_t IT_0458 = 2*conj(IT_0225);
    const ccomplex_t IT_0459 = s_24*IT_0252;
    const ccomplex_t IT_0460 = s_12*s_13*s_34;
    const ccomplex_t IT_0461 = -IT_0239 + 2*IT_0459 + (-2)*IT_0460;
    const ccomplex_t IT_0462 = IT_0238 + IT_0461;
    const ccomplex_t IT_0463 = IT_0005*IT_0462;
    const ccomplex_t IT_0464 = 2*IT_0246;
    const ccomplex_t IT_0465 = (-2)*IT_0237;
    const ccomplex_t IT_0466 = IT_0463 + IT_0464 + IT_0465;
    const ccomplex_t IT_0467 = -IT_0383;
    const ccomplex_t IT_0468 = IT_0267 + IT_0467;
    const ccomplex_t IT_0469 = IT_0005*IT_0468;
    const ccomplex_t IT_0470 = -IT_0016;
    const ccomplex_t IT_0471 = IT_0268 + IT_0470;
    const ccomplex_t IT_0472 = (-4)*IT_0471;
    const ccomplex_t IT_0473 = IT_0469 + IT_0472;
    const ccomplex_t IT_0474 = m_N_4*IT_0473;
    const ccomplex_t IT_0475 = conj(IT_0411) + conj(IT_0414) + conj(IT_0416) +
       conj(IT_0419);
    const ccomplex_t IT_0476 = 2*IT_0279;
    const ccomplex_t IT_0477 = IT_0411 + IT_0414 + IT_0416 + IT_0419;
    const ccomplex_t IT_0478 = 2*conj(IT_0279);
    const ccomplex_t IT_0479 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0480 = 0.5*IT_0226 + (-0.5)*IT_0383 + -IT_0479;
    const ccomplex_t IT_0481 = IT_0252 + IT_0480;
    const ccomplex_t IT_0482 = IT_0005*IT_0481;
    const ccomplex_t IT_0483 = IT_0016 + IT_0228 + IT_0482;
    const ccomplex_t IT_0484 = m_N_4*IT_0483;
    const ccomplex_t IT_0485 = 4*IT_0225;
    const ccomplex_t IT_0486 = 4*conj(IT_0225);
    const ccomplex_t IT_0487 = s_14*s_24;
    const ccomplex_t IT_0488 = (-2)*IT_0487;
    const ccomplex_t IT_0489 = IT_0450 + IT_0488;
    const ccomplex_t IT_0490 = IT_0005*IT_0489;
    const ccomplex_t IT_0491 = (-0.5)*IT_0490;
    const ccomplex_t IT_0492 = (-0.5)*IT_0454;
    const ccomplex_t IT_0493 = s_12*IT_0303;
    const ccomplex_t IT_0494 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0495 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0496 = s_13*s_24*s_34;
    const ccomplex_t IT_0497 = (-2)*IT_0494 + (-2)*IT_0495 + 4*IT_0496;
    const ccomplex_t IT_0498 = IT_0493 + IT_0497;
    const ccomplex_t IT_0499 = IT_0000*IT_0498;
    const ccomplex_t IT_0500 = 0.25*IT_0499;
    const ccomplex_t IT_0501 = 8*IT_0433;
    const ccomplex_t IT_0502 = 8*IT_0411;
    const ccomplex_t IT_0503 = IT_0414 + IT_0416;
    const ccomplex_t IT_0504 = 8*conj(IT_0411);
    const ccomplex_t IT_0505 = conj(IT_0414) + conj(IT_0416) + conj(IT_0419);
    const ccomplex_t IT_0506 = 8*IT_0414;
    const ccomplex_t IT_0507 = 8*conj(IT_0414);
    const ccomplex_t IT_0508 = conj(IT_0416) + conj(IT_0419);
    const ccomplex_t IT_0509 = 8*IT_0416;
    const ccomplex_t IT_0510 = 8*conj(IT_0416);
    const ccomplex_t IT_0511 = conj(IT_0433)*IT_0501 + IT_0475*IT_0502 +
       IT_0503*IT_0504 + IT_0505*IT_0506 + IT_0416*IT_0507 + IT_0508*IT_0509 + 8
      *IT_0419*(conj(IT_0419) + 0.125*IT_0504 + 0.125*IT_0507 + 0.125*IT_0510);
    const ccomplex_t IT_0512 = s_14*s_23*s_34;
    const ccomplex_t IT_0513 = (-2)*IT_0494 + (-2)*IT_0495 + 4*IT_0512;
    const ccomplex_t IT_0514 = IT_0493 + IT_0513;
    const ccomplex_t IT_0515 = IT_0000*IT_0514;
    const ccomplex_t IT_0516 = 0.25*IT_0515;
    const ccomplex_t IT_0517 = -IT_0316;
    const ccomplex_t IT_0518 = IT_0336 + IT_0517;
    const ccomplex_t IT_0519 = IT_0349 + IT_0354 + IT_0360;
    const ccomplex_t IT_0520 = conj(IT_0349) + conj(IT_0354);
    const ccomplex_t IT_0521 = conj(IT_0360) + IT_0520;
    const ccomplex_t IT_0522 = -conj(IT_0345);
    const ccomplex_t IT_0523 = IT_0521 + IT_0522;
    const ccomplex_t IT_0524 = 8*IT_0336;
    const ccomplex_t IT_0525 = 8*IT_0349;
    const ccomplex_t IT_0526 = IT_0354 + IT_0360;
    const ccomplex_t IT_0527 = 8*conj(IT_0349);
    const ccomplex_t IT_0528 = 8*IT_0354;
    const ccomplex_t IT_0529 = 8*conj(IT_0354);
    const ccomplex_t IT_0530 = 8*IT_0360;
    const ccomplex_t IT_0531 = conj(IT_0316)*IT_0518 + conj(IT_0345)*IT_0519 +
       IT_0345*IT_0523 + conj(IT_0336)*(IT_0316 + (-0.125)*IT_0524) + (-0.125)
      *IT_0521*IT_0525 + (-0.125)*IT_0526*IT_0527 + (-0.125)*(conj(IT_0354) +
       conj(IT_0360))*IT_0528 + (-0.125)*IT_0360*IT_0529 + (-0.125)*conj(IT_0360
      )*IT_0530;
    const ccomplex_t IT_0532 = 1.125*IT_0016;
    const ccomplex_t IT_0533 = (-16)*IT_0336;
    const ccomplex_t IT_0534 = (-16)*conj(IT_0336);
    const ccomplex_t IT_0535 = (-16)*IT_0433;
    const ccomplex_t IT_0536 = (-16)*conj(IT_0433);
    const ccomplex_t IT_0537 = 16*conj(IT_0336)*IT_0345 + 16*IT_0336*conj
      (IT_0345) + conj(IT_0316)*((-16)*IT_0345 + 16*IT_0349 + 16*IT_0354 + 16
      *IT_0360) + IT_0316*((-16)*conj(IT_0345) + 16*conj(IT_0349) + 16*conj
      (IT_0354) + 16*conj(IT_0360)) + IT_0521*IT_0533 + IT_0519*IT_0534 +
       IT_0475*IT_0535 + IT_0477*IT_0536;
    const ccomplex_t IT_0538 = IT_0005*IT_0450;
    const ccomplex_t IT_0539 = (-0.5)*IT_0538;
    const ccomplex_t IT_0540 = s_14*s_23;
    const ccomplex_t IT_0541 = s_12*s_34;
    const ccomplex_t IT_0542 = IT_0540 + IT_0541;
    const ccomplex_t IT_0543 = s_13*s_24;
    const ccomplex_t IT_0544 = -IT_0543;
    const ccomplex_t IT_0545 = IT_0542 + IT_0544;
    const ccomplex_t IT_0546 = IT_0001*IT_0545;
    const ccomplex_t IT_0547 = 0.25*IT_0546;
    const ccomplex_t IT_0548 = s_12 + IT_0539 + IT_0547;
    const ccomplex_t IT_0549 = -conj(IT_0316);
    const ccomplex_t IT_0550 = conj(IT_0336) + IT_0549;
    const ccomplex_t IT_0551 = 8*conj(IT_0184);
    const ccomplex_t IT_0552 = (-8)*conj(IT_0186);
    const ccomplex_t IT_0553 = IT_0551 + IT_0552;
    const ccomplex_t IT_0554 = 8*IT_0216;
    const ccomplex_t IT_0555 = (-8)*IT_0218;
    const ccomplex_t IT_0556 = IT_0554 + IT_0555;
    const ccomplex_t IT_0557 = 8*conj(IT_0216);
    const ccomplex_t IT_0558 = (-8)*conj(IT_0218);
    const ccomplex_t IT_0559 = IT_0557 + IT_0558;
    const ccomplex_t IT_0560 = (-8)*IT_0186;
    const ccomplex_t IT_0561 = (conj(IT_0184) + -conj(IT_0186))*IT_0316 + 
      -IT_0186*conj(IT_0316) + (conj(IT_0216) + -conj(IT_0218))*IT_0345 + 
      (IT_0216 + -IT_0218)*conj(IT_0345) + -IT_0184*IT_0550 + (-0.125)*IT_0336
      *IT_0553 + (-0.125)*IT_0521*IT_0556 + (-0.125)*IT_0519*IT_0559 + (-0.125)
      *conj(IT_0336)*IT_0560;
    const ccomplex_t IT_0562 = -IT_0541 + -IT_0543;
    const ccomplex_t IT_0563 = IT_0540 + IT_0562;
    const ccomplex_t IT_0564 = IT_0001*IT_0563;
    const ccomplex_t IT_0565 = (-0.25)*IT_0564;
    const ccomplex_t IT_0566 = s_12 + IT_0539 + IT_0565;
    const ccomplex_t IT_0567 = s_12*IT_0398;
    const ccomplex_t IT_0568 = (-2)*IT_0494 + (-2)*IT_0495 + 2*IT_0496 + 2
      *IT_0512 + (-2)*IT_0567;
    const ccomplex_t IT_0569 = IT_0493 + IT_0568;
    const ccomplex_t IT_0570 = IT_0000*IT_0569;
    const ccomplex_t IT_0571 = (-0.25)*IT_0570;
    const ccomplex_t IT_0572 = -IT_0259;
    const ccomplex_t IT_0573 = (-8)*IT_0411;
    const ccomplex_t IT_0574 = (-8)*IT_0414;
    const ccomplex_t IT_0575 = (-8)*IT_0416;
    const ccomplex_t IT_0576 = (-8)*IT_0419;
    const ccomplex_t IT_0577 = IT_0573 + IT_0574 + IT_0575 + IT_0576;
    const ccomplex_t IT_0578 = (-8)*conj(IT_0411);
    const ccomplex_t IT_0579 = (-8)*conj(IT_0414);
    const ccomplex_t IT_0580 = (-8)*conj(IT_0416);
    const ccomplex_t IT_0581 = (-8)*conj(IT_0419);
    const ccomplex_t IT_0582 = IT_0578 + IT_0579 + IT_0580 + IT_0581;
    const ccomplex_t IT_0583 = 8*conj(IT_0336);
    const ccomplex_t IT_0584 = (-8)*IT_0433;
    const ccomplex_t IT_0585 = (-8)*conj(IT_0433);
    const ccomplex_t IT_0586 = conj(IT_0354)*IT_0506 + 8*conj(IT_0360)*
      (IT_0419 + 0.125*IT_0502 + 0.125*IT_0506 + 0.125*IT_0509) + IT_0360
      *IT_0510 + IT_0504*IT_0519 + IT_0502*IT_0520 + conj(IT_0433)*IT_0524 +
       IT_0505*IT_0525 + IT_0507*IT_0526 + (IT_0419 + IT_0503)*IT_0527 + IT_0508
      *IT_0528 + (IT_0416 + IT_0419)*IT_0529 + conj(IT_0419)*IT_0530 + conj
      (IT_0345)*IT_0577 + IT_0345*IT_0582 + IT_0433*IT_0583 + conj(IT_0316)
      *IT_0584 + IT_0316*IT_0585;
    const ccomplex_t IT_0587 = (-2)*IT_0398;
    const ccomplex_t IT_0588 = IT_0303 + IT_0587;
    const ccomplex_t IT_0589 = IT_0000*IT_0588;
    const ccomplex_t IT_0590 = 0.25*IT_0589;
    const ccomplex_t IT_0591 = 1 + IT_0590;
    const ccomplex_t IT_0592 = IT_0016*IT_0591;
    const ccomplex_t IT_0593 = 8*conj(IT_0433);
    const ccomplex_t IT_0594 = pow(m_Z, 6);
    const ccomplex_t IT_0595 = s_12*IT_0594;
    const ccomplex_t IT_0596 = s_13*s_23*IT_0303;
    const ccomplex_t IT_0597 = (-2)*IT_0596;
    const ccomplex_t IT_0598 = IT_0595 + IT_0597;
    const ccomplex_t IT_0599 = IT_0000*IT_0598;
    const ccomplex_t IT_0600 = 0.25*IT_0599;
    const ccomplex_t IT_0601 = IT_0450 + IT_0600;
    const ccomplex_t IT_0602 = (-8)*IT_0192;
    const ccomplex_t IT_0603 = (-8)*conj(IT_0192);
    const ccomplex_t IT_0604 = (-8)*IT_0043;
    const ccomplex_t IT_0605 = (-8)*conj(IT_0043);
    const ccomplex_t IT_0606 = 6*conj(IT_0225);
    const ccomplex_t IT_0607 = (-6)*conj(IT_0279);
    const ccomplex_t IT_0608 = IT_0000*IT_0398;
    const ccomplex_t IT_0609 = 0.5*IT_0608;
    const ccomplex_t IT_0610 = 1 + IT_0609;
    const ccomplex_t IT_0611 = s_12*IT_0610;
    const ccomplex_t IT_0612 = (-4)*conj(IT_0184);
    const ccomplex_t IT_0613 = 0.25*IT_0612;
    const ccomplex_t IT_0614 = conj(IT_0186) + IT_0613;
    const ccomplex_t IT_0615 = (-4)*IT_0184;
    const ccomplex_t IT_0616 = (-4)*IT_0216;
    const ccomplex_t IT_0617 = (-4)*conj(IT_0216);
    const ccomplex_t IT_0618 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0619 = 0.5*IT_0267 + (-0.5)*IT_0383 + -IT_0618;
    const ccomplex_t IT_0620 = IT_0265 + IT_0619;
    const ccomplex_t IT_0621 = IT_0005*IT_0620;
    const ccomplex_t IT_0622 = -IT_0268;
    const ccomplex_t IT_0623 = IT_0016 + IT_0621 + IT_0622;
    const ccomplex_t IT_0624 = m_N_4*IT_0623;
    const ccomplex_t IT_0625 = (-4)*IT_0279;
    const ccomplex_t IT_0626 = (-4)*conj(IT_0279);
    const ccomplex_t IT_0627 = s_12*s_24;
    const ccomplex_t IT_0628 = s_12*s_24*IT_0002;
    const ccomplex_t IT_0629 = s_14*IT_0002*IT_0016;
    const ccomplex_t IT_0630 = s_14*IT_0265;
    const ccomplex_t IT_0631 = s_12*s_23*s_34;
    const ccomplex_t IT_0632 = -IT_0629 + 2*IT_0630 + (-2)*IT_0631;
    const ccomplex_t IT_0633 = IT_0628 + IT_0632;
    const ccomplex_t IT_0634 = IT_0005*IT_0633;
    const ccomplex_t IT_0635 = 0.5*IT_0634;
    const ccomplex_t IT_0636 = -IT_0635;
    const ccomplex_t IT_0637 = s_14*IT_0016;
    const ccomplex_t IT_0638 = -IT_0637;
    const ccomplex_t IT_0639 = IT_0627 + IT_0636 + IT_0638;
    const ccomplex_t IT_0640 = (-4)*IT_0232;
    const ccomplex_t IT_0641 = (-4)*conj(IT_0232);
    const ccomplex_t IT_0642 = 0.5*IT_0599;
    const ccomplex_t IT_0643 = (-4)*IT_0487;
    const ccomplex_t IT_0644 = IT_0450 + IT_0642 + IT_0643;
    const ccomplex_t IT_0645 = IT_0005*IT_0498;
    const ccomplex_t IT_0646 = -IT_0645;
    const ccomplex_t IT_0647 = IT_0644 + IT_0646;
    const ccomplex_t IT_0648 = (-4)*IT_0192;
    const ccomplex_t IT_0649 = (-4)*IT_0043;
    const ccomplex_t IT_0650 = IT_0005*IT_0514;
    const ccomplex_t IT_0651 = -IT_0650;
    const ccomplex_t IT_0652 = IT_0644 + IT_0651;
    const ccomplex_t IT_0653 = (-4)*IT_0096;
    const ccomplex_t IT_0654 = IT_0000*IT_0377;
    const ccomplex_t IT_0655 = (-10)*s_14;
    const ccomplex_t IT_0656 = IT_0654 + IT_0655;
    const ccomplex_t IT_0657 = m_N_4*IT_0656;
    const ccomplex_t IT_0658 = 0.666666666666667*IT_0292 + (-0.166666666666667
      )*IT_0654;
    const ccomplex_t IT_0659 = s_14 + IT_0658;
    const ccomplex_t IT_0660 = m_N_4*IT_0659;
    const ccomplex_t IT_0661 = (-6)*IT_0660;
    const ccomplex_t IT_0662 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0663 = (-2)*IT_0451 + -IT_0487 + 2*IT_0662;
    const ccomplex_t IT_0664 = IT_0450 + IT_0663;
    const ccomplex_t IT_0665 = IT_0005*IT_0664;
    const ccomplex_t IT_0666 = -IT_0233 + (-0.5)*IT_0665;
    const ccomplex_t IT_0667 = s_12 + IT_0666;
    const ccomplex_t IT_0668 = m_N_4*IT_0667;
    const ccomplex_t IT_0669 = (-2)*IT_0668;
    const ccomplex_t IT_0670 = IT_0226 + IT_0467;
    const ccomplex_t IT_0671 = IT_0005*IT_0670;
    const ccomplex_t IT_0672 = IT_0227 + IT_0470;
    const ccomplex_t IT_0673 = (-4)*IT_0672;
    const ccomplex_t IT_0674 = IT_0671 + IT_0673;
    const ccomplex_t IT_0675 = m_N_4*IT_0674;
    const ccomplex_t IT_0676 = (IT_0016 + -IT_0227)*(s_12 + (-0.5)*IT_0454);
    const ccomplex_t IT_0677 = 2*IT_0676;
    const ccomplex_t IT_0678 = (-2)*IT_0225;
    const ccomplex_t IT_0679 = 1.5*IT_0225*IT_0231 + (-1.5)*IT_0232*IT_0236 +
       0.25*IT_0188*IT_0257 + 0.25*IT_0220*IT_0262 + IT_0248*IT_0375 + 0.5
      *IT_0363*IT_0466 + 0.25*IT_0456*IT_0476 + IT_0433*IT_0484 + (-0.5)*
      (IT_0345 + -IT_0519)*IT_0669 + 0.25*IT_0338*IT_0675 + 0.25*IT_0677*IT_0678;
    const ccomplex_t IT_0680 = (-6)*IT_0279;
    const ccomplex_t IT_0681 = 2*IT_0336;
    const ccomplex_t IT_0682 = 2*IT_0232;
    const ccomplex_t IT_0683 = (IT_0197 + IT_0225)*IT_0231 + 0.166666666666667
      *IT_0220*IT_0257 + 0.166666666666667*IT_0188*IT_0262 + 0.666666666666667
      *IT_0248*IT_0374 + 0.333333333333333*IT_0339*IT_0466 + 0.666666666666667
      *IT_0477*IT_0484 + 0.166666666666667*IT_0362*IT_0675 + 0.166666666666667
      *IT_0677*IT_0678 + 0.166666666666667*IT_0236*IT_0680 + (-0.333333333333333
      )*IT_0669*(IT_0316 + (-0.5)*IT_0681) + 0.166666666666667*IT_0456*IT_0682;
    const ccomplex_t IT_0684 = 2*conj(IT_0336);
    const ccomplex_t IT_0685 = (-2)*conj(IT_0225);
    const ccomplex_t IT_0686 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0687 = (-2)*IT_0451 + -IT_0487 + 2*IT_0686;
    const ccomplex_t IT_0688 = IT_0450 + IT_0687;
    const ccomplex_t IT_0689 = IT_0005*IT_0688;
    const ccomplex_t IT_0690 = -IT_0233 + (-0.5)*IT_0689;
    const ccomplex_t IT_0691 = s_12 + IT_0690;
    const ccomplex_t IT_0692 = m_N_4*IT_0691;
    const ccomplex_t IT_0693 = (-2)*IT_0692;
    const ccomplex_t IT_0694 = (-2)*IT_0279;
    const ccomplex_t IT_0695 = (-2)*conj(IT_0279);
    const ccomplex_t IT_0696 = (-2)*IT_0232;
    const ccomplex_t IT_0697 = (-2)*conj(IT_0232);
    const ccomplex_t IT_0698 = s_13*s_23*s_24;
    const ccomplex_t IT_0699 = s_13*s_34*IT_0016;
    const ccomplex_t IT_0700 = -IT_0629 + (-2)*IT_0698 + 2*IT_0699;
    const ccomplex_t IT_0701 = IT_0628 + IT_0700;
    const ccomplex_t IT_0702 = IT_0005*IT_0701;
    const ccomplex_t IT_0703 = 2*IT_0637;
    const ccomplex_t IT_0704 = (-2)*IT_0627;
    const ccomplex_t IT_0705 = IT_0702 + IT_0703 + IT_0704;
    const ccomplex_t IT_0706 = IT_0016*IT_0610;
    const ccomplex_t IT_0707 = (-4)*IT_0186;
    const ccomplex_t IT_0708 = (-4)*conj(IT_0186);
    const ccomplex_t IT_0709 = conj(IT_0216)*IT_0316 + IT_0216*conj(IT_0316) +
       (IT_0216 + IT_0349 + IT_0354 + IT_0360 + IT_0411 + IT_0414 + IT_0416 +
       IT_0419)*IT_0614 + (conj(IT_0216) + conj(IT_0349) + conj(IT_0354) + conj
      (IT_0360) + conj(IT_0411) + conj(IT_0414) + conj(IT_0416) + conj(IT_0419))
      *(IT_0186 + 0.25*IT_0615) + 0.25*(conj(IT_0336) + conj(IT_0433))*IT_0616 +
       0.25*(IT_0336 + IT_0433)*IT_0617 + conj(IT_0345)*(IT_0184 + 0.25*IT_0707)
       + conj(IT_0218)*(IT_0184 + IT_0336 + IT_0433 + IT_0517 + 0.25*IT_0707) +
       IT_0345*(conj(IT_0184) + 0.25*IT_0708) + IT_0218*(conj(IT_0184) + conj
      (IT_0336) + conj(IT_0433) + IT_0549 + 0.25*IT_0708);
    const ccomplex_t IT_0710 = -IT_0487;
    const ccomplex_t IT_0711 = IT_0450 + IT_0710;
    const ccomplex_t IT_0712 = IT_0005*IT_0711;
    const ccomplex_t IT_0713 = -IT_0233 + (-0.25)*IT_0712;
    const ccomplex_t IT_0714 = s_12 + IT_0713;
    const ccomplex_t IT_0715 = m_N_4*IT_0714;
    const ccomplex_t IT_0716 = (-4)*IT_0715;
    const ccomplex_t IT_0717 = (-0.5)*IT_0685;
    const ccomplex_t IT_0718 = conj(IT_0197) + IT_0717;
    const ccomplex_t IT_0719 = 2*conj(IT_0232);
    const ccomplex_t IT_0720 = (-2)*IT_0433;
    const ccomplex_t IT_0721 = (-2)*conj(IT_0433);
    const ccomplex_t IT_0722 = (-2)*IT_0475;
    const ccomplex_t IT_0723 = (-0.5)*conj(IT_0336)*IT_0476 + (-0.5)*IT_0336
      *IT_0478 + (-0.5)*IT_0521*IT_0682 + (-0.5)*conj(IT_0316)*IT_0694 + (-0.5)
      *IT_0316*IT_0695 + (-0.5)*conj(IT_0345)*IT_0696 + (-0.5)*IT_0345*IT_0697 +
       IT_0477*IT_0718 + (-0.5)*IT_0519*IT_0719 + (-0.5)*(conj(IT_0127) + conj
      (IT_0225))*IT_0720 + (-0.5)*(IT_0127 + IT_0225)*IT_0721 + (-0.5)*(IT_0197 
      + (-0.5)*IT_0678)*IT_0722;
    const ccomplex_t IT_0724 = IT_0016 + IT_0622;
    const ccomplex_t IT_0725 = IT_0016*IT_0724;
    const ccomplex_t IT_0726 = -IT_0725;
    const ccomplex_t IT_0727 = (IT_0016 + -IT_0268)*(s_12 + (-0.5)*IT_0454);
    const ccomplex_t IT_0728 = 2*IT_0727;
    const ccomplex_t IT_0729 = 0.166666666666667*IT_0728;
    const ccomplex_t IT_0730 = (-0.5)*IT_0262;
    const ccomplex_t IT_0731 = (-0.5)*IT_0456;
    const ccomplex_t IT_0732 = (-0.5)*IT_0675;
    const ccomplex_t IT_0733 = 1.5*conj(IT_0225)*IT_0231 + (-1.5)*conj(IT_0232
      )*IT_0236 + IT_0248*IT_0249 + 0.25*IT_0191*IT_0257 + 0.5*IT_0343*IT_0466 +
       conj(IT_0433)*IT_0484 + 0.5*IT_0523*IT_0669 + (-0.5)*IT_0677*IT_0718 + (
      -0.5)*IT_0223*IT_0730 + (-0.5)*IT_0478*IT_0731 + (-0.5)*IT_0342*IT_0732;
    const ccomplex_t IT_0734 = 0.5*IT_0223*IT_0257 + 2*IT_0248*IT_0373 +
       IT_0313*IT_0466 + 3*IT_0231*(conj(IT_0197) + 0.166666666666667*IT_0606) +
       0.5*IT_0236*IT_0607 + -IT_0669*(conj(IT_0316) + (-0.5)*IT_0684) + 
      -IT_0677*(conj(IT_0127) + (-0.5)*IT_0685) + -IT_0484*IT_0722 + -IT_0191
      *IT_0730 + -IT_0719*IT_0731 + -IT_0368*IT_0732;
    const ccomplex_t IT_0735 = IT_0224 + 6*IT_0225*(conj(IT_0225)*IT_0231 + 
      -conj(IT_0232)*IT_0236 + 0.666666666666667*IT_0248*IT_0249 +
       0.166666666666667*IT_0191*IT_0257 + 0.166666666666667*IT_0223*IT_0262) +
       IT_0272*(conj(IT_0232)*IT_0275 + IT_0232*IT_0278 + conj(IT_0279)*IT_0282 
      + IT_0279*IT_0285) + (IT_0278*IT_0279 + IT_0275*conj(IT_0279) + conj
      (IT_0232)*IT_0282 + IT_0232*IT_0285)*IT_0288 + IT_0301 + IT_0312*(IT_0313
      *IT_0338 + IT_0339*IT_0342 + IT_0343*IT_0362 + IT_0363*IT_0368) + IT_0372*
      (IT_0249*IT_0362 + IT_0338*IT_0373 + IT_0342*IT_0374 + IT_0368*IT_0375) + 
      (IT_0338*IT_0343 + IT_0313*IT_0362 + IT_0342*IT_0363 + IT_0339*IT_0368)
      *IT_0382 + (-2)*IT_0384*(IT_0373*IT_0387 + IT_0374*IT_0390 + (-0.5)
      *IT_0313*IT_0393 + (-0.5)*IT_0339*IT_0396) + (-2)*(IT_0313*IT_0387 +
       IT_0339*IT_0390 + (-0.5)*IT_0373*IT_0393 + (-0.5)*IT_0374*IT_0396)
      *IT_0402 + IT_0410*(IT_0313*IT_0421 + IT_0339*IT_0426 + IT_0343*IT_0434 +
       IT_0363*IT_0435) + (IT_0373*IT_0421 + IT_0374*IT_0426 + IT_0249*IT_0434 +
       IT_0375*IT_0435)*IT_0440 + (IT_0343*IT_0421 + IT_0363*IT_0426 + IT_0313
      *IT_0434 + IT_0339*IT_0435)*IT_0444 + (IT_0249*IT_0421 + IT_0375*IT_0426 +
       IT_0373*IT_0434 + IT_0374*IT_0435)*IT_0449 + IT_0456*(conj(IT_0232)
      *IT_0457 + conj(IT_0279)*IT_0457 + IT_0232*IT_0458 + IT_0279*IT_0458) + 
      (IT_0313*IT_0457 + IT_0343*IT_0457 + IT_0339*IT_0458 + IT_0363*IT_0458)
      *IT_0466 + IT_0474*(conj(IT_0232)*IT_0434 + IT_0232*IT_0435 + IT_0475
      *IT_0476 + IT_0477*IT_0478) + IT_0248*(IT_0373*IT_0485 + IT_0374*IT_0486) 
      + 4*IT_0484*(conj(IT_0225)*IT_0433 + IT_0225*conj(IT_0433) + 0.25*IT_0475
      *IT_0485 + 0.25*IT_0477*IT_0486) + (s_12 + IT_0491 + IT_0492 + IT_0500)
      *IT_0511 + (-8)*(s_12 + IT_0491 + IT_0492 + IT_0516)*IT_0531 + IT_0532
      *IT_0537 + (-8)*IT_0548*IT_0561 + (IT_0433*IT_0553 + IT_0475*IT_0556 +
       IT_0477*IT_0559 + 8*conj(IT_0433)*(IT_0184 + 0.125*IT_0560))*IT_0566 + 
      (s_12 + IT_0234 + IT_0571 + IT_0572)*IT_0586 + IT_0592*(IT_0501*IT_0521 +
       IT_0475*IT_0524 + conj(IT_0316)*IT_0577 + IT_0316*IT_0582 + IT_0477
      *IT_0583 + conj(IT_0345)*IT_0584 + IT_0345*IT_0585 + IT_0519*IT_0593) + (
      -8)*IT_0601*(IT_0249*(IT_0194 + (-0.125)*IT_0602) + IT_0375*(conj(IT_0194)
       + (-0.125)*IT_0603) + IT_0373*(IT_0121 + (-0.125)*IT_0604) + IT_0374*
      (conj(IT_0121) + (-0.125)*IT_0605)) + IT_0225*(IT_0223*IT_0257 + IT_0191
      *IT_0262 + IT_0231*IT_0606 + IT_0236*IT_0607) + 4*IT_0611*(IT_0184*conj
      (IT_0184) + IT_0216*conj(IT_0216) + IT_0186*IT_0614 + 0.25*conj(IT_0186)
      *IT_0615 + 0.25*conj(IT_0218)*IT_0616 + IT_0218*(conj(IT_0218) + 0.25
      *IT_0617)) + 4*IT_0624*(conj(IT_0279)*IT_0345 + IT_0279*conj(IT_0345) + 
      -conj(IT_0232)*IT_0518 + -IT_0232*IT_0550 + 0.25*IT_0521*IT_0625 + 0.25
      *IT_0519*IT_0626) + IT_0639*(IT_0249*IT_0625 + IT_0375*IT_0626 + IT_0373
      *IT_0640 + IT_0374*IT_0641) + (-4)*IT_0647*(IT_0343*(IT_0194 + (-0.25)
      *IT_0648) + IT_0313*(IT_0121 + (-0.25)*IT_0649)) + (-4)*IT_0652*(IT_0373*
      (IT_0128 + (-0.25)*IT_0653) + IT_0249*(IT_0198 + (-0.25)*IT_0653)) + 
      (IT_0249*IT_0317 + conj(IT_0198)*IT_0337 + IT_0198*IT_0341 + conj(IT_0096)
      *IT_0346 + conj(IT_0128)*IT_0362 + IT_0096*IT_0364 + IT_0128*IT_0368 +
       IT_0340*IT_0375)*IT_0657 + (conj(IT_0096)*IT_0337 + IT_0096*IT_0341 +
       conj(IT_0096)*(IT_0350 + IT_0355 + IT_0361) + IT_0096*(IT_0365 + IT_0366 
      + IT_0367))*IT_0661 + (-6)*conj(IT_0225)*(IT_0232*IT_0236 + (
      -0.166666666666667)*IT_0188*IT_0257 + (-0.166666666666667)*IT_0220*IT_0262
       + (-0.666666666666667)*IT_0248*IT_0375 + (-0.166666666666667)*IT_0677
      *IT_0678) + 4*conj(IT_0197)*IT_0679 + conj(IT_0225)*(IT_0220*IT_0257 +
       IT_0188*IT_0262 + IT_0677*IT_0678 + IT_0236*IT_0680) + 6*conj(IT_0127)
      *IT_0683 + IT_0669*(IT_0458*IT_0519 + IT_0457*IT_0521 + conj(IT_0316)
      *IT_0678 + conj(IT_0345)*IT_0678 + conj(IT_0225)*IT_0681 + IT_0225*IT_0684
       + IT_0316*IT_0685 + IT_0345*IT_0685) + IT_0675*(conj(IT_0316)*IT_0457 +
       conj(IT_0345)*IT_0457 + IT_0316*IT_0458 + IT_0345*IT_0458 + conj(IT_0336)
      *IT_0678 + IT_0521*IT_0678 + IT_0336*IT_0685 + IT_0519*IT_0685) + IT_0693*
      (conj(IT_0433)*IT_0694 + IT_0433*IT_0695 + IT_0475*IT_0696 + IT_0477
      *IT_0697) + (IT_0343*IT_0694 + IT_0363*IT_0695 + IT_0313*IT_0696 + IT_0339
      *IT_0697)*IT_0705 + 4*IT_0706*IT_0709 + (-2)*IT_0716*IT_0723 + 6*conj
      (IT_0279)*(IT_0232*IT_0726 + IT_0694*IT_0729) + 6*conj(IT_0232)*(IT_0279
      *IT_0726 + IT_0696*IT_0729) + 4*IT_0197*IT_0733 + 2*IT_0127*IT_0734;
    return create_ccomplex_return(IT_0735);
}

