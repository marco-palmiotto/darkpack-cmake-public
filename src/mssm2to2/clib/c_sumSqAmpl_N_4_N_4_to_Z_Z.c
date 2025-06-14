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
    const ccomplex_t IT_0001 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0004 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0005 = -IT_0003 + -IT_0004;
    const ccomplex_t IT_0006 = s_12 + IT_0001 + IT_0005;
    const ccomplex_t IT_0007 = m_N_4*IT_0006;
    const ccomplex_t IT_0008 = pow(m_N_4, 2);
    const ccomplex_t IT_0009 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0010 = IT_0008 + IT_0009;
    const ccomplex_t IT_0011 = pow(s_23, 2);
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = pow(s_24, 2);
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = -IT_0012 + -IT_0014;
    const ccomplex_t IT_0016 = IT_0010 + IT_0015;
    const ccomplex_t IT_0017 = m_N_4*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0021, -1);
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0028;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + IT_0027 + 
      -IT_0029 + -IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0034 = IT_0022*IT_0033;
    const ccomplex_t IT_0035 = IT_0026*IT_0033;
    const ccomplex_t IT_0036 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0037 = IT_0022*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0036;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0034 + IT_0035 + 
      -IT_0037 + -IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0032*IT_0040;
    const ccomplex_t IT_0042 = pow(m_Z, 2);
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0008 + IT_0042 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0047 = IT_0022*IT_0046;
    const ccomplex_t IT_0048 = IT_0026*IT_0046;
    const ccomplex_t IT_0049 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = IT_0026*IT_0049;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0047 + IT_0048 + 
      -IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0055 = IT_0022*IT_0054;
    const ccomplex_t IT_0056 = IT_0026*IT_0054;
    const ccomplex_t IT_0057 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0058 = IT_0022*IT_0057;
    const ccomplex_t IT_0059 = IT_0026*IT_0057;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0055 + IT_0056 + 
      -IT_0058 + -IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0008 + IT_0042 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0067 = IT_0022*IT_0066;
    const ccomplex_t IT_0068 = IT_0026*IT_0066;
    const ccomplex_t IT_0069 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0070 = IT_0022*IT_0069;
    const ccomplex_t IT_0071 = IT_0026*IT_0069;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0067 + IT_0068 + 
      -IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0075 = IT_0022*IT_0074;
    const ccomplex_t IT_0076 = IT_0026*IT_0074;
    const ccomplex_t IT_0077 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0078 = IT_0022*IT_0077;
    const ccomplex_t IT_0079 = IT_0026*IT_0077;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0075 + IT_0076 + 
      -IT_0078 + -IT_0079);
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0073*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0008 + IT_0042 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = (-2)*IT_0045 + (-2)*IT_0065 + (-2)*IT_0085;
    const ccomplex_t IT_0087 = cos(alpha);
    const ccomplex_t IT_0088 = IT_0020*IT_0087;
    const ccomplex_t IT_0089 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = sin(alpha);
    const ccomplex_t IT_0092 = IT_0020*IT_0091;
    const ccomplex_t IT_0093 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = IT_0025*IT_0087;
    const ccomplex_t IT_0096 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = IT_0025*IT_0091;
    const ccomplex_t IT_0099 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0090 + IT_0094 + 
      -IT_0097 + -IT_0100);
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = m_W*e_em;
    const ccomplex_t IT_0104 = cpow(IT_0019, -4);
    const ccomplex_t IT_0105 = sin(beta);
    const ccomplex_t IT_0106 = IT_0021*IT_0087*IT_0104*IT_0105;
    const ccomplex_t IT_0107 = cos(beta);
    const ccomplex_t IT_0108 = IT_0021*IT_0091*IT_0104*IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0021, -3);
    const ccomplex_t IT_0110 = IT_0087*IT_0105*IT_0109;
    const ccomplex_t IT_0111 = IT_0091*IT_0107*IT_0109;
    const ccomplex_t IT_0112 = cpow(IT_0019, -2);
    const ccomplex_t IT_0113 = IT_0025*IT_0087*IT_0105*IT_0112;
    const ccomplex_t IT_0114 = IT_0025*IT_0091*IT_0107*IT_0112;
    const ccomplex_t IT_0115 = cpow(IT_0021, -2);
    const ccomplex_t IT_0116 = cpow(IT_0112 + IT_0115, -1);
    const ccomplex_t IT_0117 = (IT_0106 + -IT_0108 + IT_0110 + -IT_0111 + 2
      *IT_0113 + (-2)*IT_0114)*IT_0116;
    const ccomplex_t IT_0118 = IT_0103*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0102*IT_0122;
    const ccomplex_t IT_0124 = IT_0088*IT_0093;
    const ccomplex_t IT_0125 = IT_0089*IT_0092;
    const ccomplex_t IT_0126 = IT_0095*IT_0099;
    const ccomplex_t IT_0127 = IT_0096*IT_0098;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + -IT_0125 + 
      -IT_0126 + IT_0127);
    const ccomplex_t IT_0129 = IT_0021*IT_0087*IT_0104*IT_0107;
    const ccomplex_t IT_0130 = IT_0021*IT_0091*IT_0104*IT_0105;
    const ccomplex_t IT_0131 = IT_0087*IT_0107*IT_0109;
    const ccomplex_t IT_0132 = IT_0091*IT_0105*IT_0109;
    const ccomplex_t IT_0133 = IT_0025*IT_0087*IT_0107*IT_0112;
    const ccomplex_t IT_0134 = IT_0025*IT_0091*IT_0105*IT_0112;
    const ccomplex_t IT_0135 = IT_0116*(IT_0129 + IT_0130 + IT_0131 + IT_0132 
      + 2*IT_0133 + 2*IT_0134);
    const ccomplex_t IT_0136 = IT_0103*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_12 + (-2)*IT_0008 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0128*IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0123 + (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = 2*IT_0142;
    const ccomplex_t IT_0144 = conj(IT_0086)*IT_0143;
    const ccomplex_t IT_0145 = 2*conj(IT_0142);
    const ccomplex_t IT_0146 = IT_0086*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0123 + 0.5*IT_0141;
    const ccomplex_t IT_0148 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0149 = IT_0086*IT_0148;
    const ccomplex_t IT_0150 = IT_0144 + IT_0146 + IT_0149;
    const ccomplex_t IT_0151 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0152 = IT_0008 + IT_0151;
    const ccomplex_t IT_0153 = pow(s_13, 2);
    const ccomplex_t IT_0154 = IT_0002*IT_0153;
    const ccomplex_t IT_0155 = pow(s_14, 2);
    const ccomplex_t IT_0156 = IT_0002*IT_0155;
    const ccomplex_t IT_0157 = -IT_0154 + -IT_0156;
    const ccomplex_t IT_0158 = IT_0152 + IT_0157;
    const ccomplex_t IT_0159 = m_N_4*IT_0158;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0162 = IT_0022*IT_0161;
    const ccomplex_t IT_0163 = IT_0026*IT_0161;
    const ccomplex_t IT_0164 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0165 = IT_0022*IT_0164;
    const ccomplex_t IT_0166 = IT_0026*IT_0164;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0162 + IT_0163 + 
      -IT_0165 + -IT_0166);
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = cpow(IT_0168, 2);
    const ccomplex_t IT_0170 = cpow((-2)*s_23 + IT_0042 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = (-2)*IT_0172;
    const ccomplex_t IT_0174 = cpow((-2)*s_23 + IT_0008 + IT_0042 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0082*IT_0175;
    const ccomplex_t IT_0177 = cpow((-2)*s_23 + IT_0008 + IT_0042 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0062*IT_0178;
    const ccomplex_t IT_0180 = cpow((-2)*s_23 + IT_0008 + IT_0042 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = IT_0041*IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0176 + 2*IT_0179 + 2*IT_0182;
    const ccomplex_t IT_0184 = conj(IT_0173) + conj(IT_0183);
    const ccomplex_t IT_0185 = 2*IT_0147;
    const ccomplex_t IT_0186 = (-2)*IT_0142;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = IT_0173 + IT_0183;
    const ccomplex_t IT_0189 = 2*conj(IT_0147);
    const ccomplex_t IT_0190 = (-2)*conj(IT_0142);
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = (-0.5)*IT_0167;
    const ccomplex_t IT_0193 = cpow(IT_0192, 2);
    const ccomplex_t IT_0194 = IT_0171*IT_0193;
    const ccomplex_t IT_0195 = (-2)*IT_0194;
    const ccomplex_t IT_0196 = conj(IT_0183) + conj(IT_0195);
    const ccomplex_t IT_0197 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0198 = IT_0088*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0200 = IT_0092*IT_0199;
    const ccomplex_t IT_0201 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0202 = IT_0095*IT_0201;
    const ccomplex_t IT_0203 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0204 = IT_0098*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0198 + -IT_0200 + 
      -IT_0202 + IT_0204);
    const ccomplex_t IT_0206 = IT_0140*IT_0205;
    const ccomplex_t IT_0207 = IT_0088*IT_0199;
    const ccomplex_t IT_0208 = IT_0092*IT_0197;
    const ccomplex_t IT_0209 = IT_0095*IT_0203;
    const ccomplex_t IT_0210 = IT_0098*IT_0201;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + IT_0208 + 
      -IT_0209 + -IT_0210);
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0122*IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0206 + 0.5*IT_0213;
    const ccomplex_t IT_0215 = 2*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0206 + (-0.5)*IT_0213;
    const ccomplex_t IT_0217 = (-2)*IT_0216;
    const ccomplex_t IT_0218 = IT_0215 + IT_0217;
    const ccomplex_t IT_0219 = 2*conj(IT_0214);
    const ccomplex_t IT_0220 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0223 = s_12 + IT_0005 + IT_0222;
    const ccomplex_t IT_0224 = m_N_4*IT_0223;
    const ccomplex_t IT_0225 = s_34*IT_0000;
    const ccomplex_t IT_0226 = s_23*IT_0042;
    const ccomplex_t IT_0227 = IT_0225*IT_0226;
    const ccomplex_t IT_0228 = s_24*IT_0042;
    const ccomplex_t IT_0229 = IT_0002*IT_0228;
    const ccomplex_t IT_0230 = s_23*s_34;
    const ccomplex_t IT_0231 = (-2)*IT_0230;
    const ccomplex_t IT_0232 = IT_0228 + IT_0231;
    const ccomplex_t IT_0233 = IT_0002*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0227 + 0.5*IT_0229 + (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = s_24 + IT_0234;
    const ccomplex_t IT_0236 = m_N_4*IT_0235;
    const ccomplex_t IT_0237 = (-2)*IT_0236;
    const ccomplex_t IT_0238 = -IT_0176 + -IT_0182;
    const ccomplex_t IT_0239 = IT_0045 + IT_0065;
    const ccomplex_t IT_0240 = cpow((-2)*s_13 + IT_0042 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = IT_0193*IT_0241;
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = -IT_0179;
    const ccomplex_t IT_0245 = conj(IT_0194) + conj(IT_0244);
    const ccomplex_t IT_0246 = conj(IT_0085) + conj(IT_0238) + conj(IT_0239) +
       conj(IT_0243) + IT_0245;
    const ccomplex_t IT_0247 = IT_0194 + IT_0244;
    const ccomplex_t IT_0248 = IT_0085 + IT_0238 + IT_0239 + IT_0243 + IT_0247;
    const ccomplex_t IT_0249 = IT_0045 + IT_0065 + IT_0085;
    const ccomplex_t IT_0250 = IT_0169*IT_0241;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = -IT_0176 + -IT_0179 + -IT_0182;
    const ccomplex_t IT_0253 = conj(IT_0172) + conj(IT_0252);
    const ccomplex_t IT_0254 = conj(IT_0249) + conj(IT_0251) + IT_0253;
    const ccomplex_t IT_0255 = IT_0172 + IT_0252;
    const ccomplex_t IT_0256 = IT_0249 + IT_0251 + IT_0255;
    const ccomplex_t IT_0257 = s_13*IT_0042;
    const ccomplex_t IT_0258 = IT_0225*IT_0257;
    const ccomplex_t IT_0259 = s_14*IT_0042;
    const ccomplex_t IT_0260 = IT_0002*IT_0259;
    const ccomplex_t IT_0261 = s_13*s_34;
    const ccomplex_t IT_0262 = (-2)*IT_0261;
    const ccomplex_t IT_0263 = IT_0259 + IT_0262;
    const ccomplex_t IT_0264 = IT_0002*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0258 + 0.5*IT_0260 + (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = s_14 + IT_0265;
    const ccomplex_t IT_0267 = m_N_4*IT_0266;
    const ccomplex_t IT_0268 = (-2)*IT_0267;
    const ccomplex_t IT_0269 = (-2)*IT_0147;
    const ccomplex_t IT_0270 = IT_0143 + IT_0269;
    const ccomplex_t IT_0271 = IT_0145 + IT_0148;
    const ccomplex_t IT_0272 = 2*IT_0216;
    const ccomplex_t IT_0273 = (-2)*IT_0214;
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = 2*conj(IT_0216);
    const ccomplex_t IT_0276 = (-2)*conj(IT_0214);
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0279 = pow(m_Z, 4);
    const ccomplex_t IT_0280 = s_24*IT_0279;
    const ccomplex_t IT_0281 = s_23*s_34*IT_0042;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = IT_0280 + IT_0282;
    const ccomplex_t IT_0284 = IT_0000*IT_0283;
    const ccomplex_t IT_0285 = (-0.666666666666667)*IT_0278 + (
      -0.166666666666667)*IT_0284;
    const ccomplex_t IT_0286 = s_24 + IT_0285;
    const ccomplex_t IT_0287 = m_N_4*IT_0286;
    const ccomplex_t IT_0288 = (-6)*IT_0287;
    const ccomplex_t IT_0289 = conj(IT_0085) + conj(IT_0239);
    const ccomplex_t IT_0290 = conj(IT_0243) + IT_0289;
    const ccomplex_t IT_0291 = IT_0168*IT_0192;
    const ccomplex_t IT_0292 = m_N_4*IT_0291;
    const ccomplex_t IT_0293 = m_N_4*IT_0193;
    const ccomplex_t IT_0294 = IT_0171*(IT_0292 + -IT_0293);
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = -IT_0295;
    const ccomplex_t IT_0297 = 2*IT_0296;
    const ccomplex_t IT_0298 = cpow(IT_0061, 2);
    const ccomplex_t IT_0299 = m_N_1*IT_0298;
    const ccomplex_t IT_0300 = m_N_4*IT_0062;
    const ccomplex_t IT_0301 = IT_0178*(IT_0299 + -IT_0300);
    const ccomplex_t IT_0302 = (-2)*IT_0301;
    const ccomplex_t IT_0303 = cpow(IT_0040, 2);
    const ccomplex_t IT_0304 = m_N_2*IT_0303;
    const ccomplex_t IT_0305 = m_N_4*IT_0041;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = IT_0181*IT_0307;
    const ccomplex_t IT_0309 = cpow(IT_0081, 2);
    const ccomplex_t IT_0310 = m_N_3*IT_0309;
    const ccomplex_t IT_0311 = m_N_4*IT_0082;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = IT_0175*IT_0313;
    const ccomplex_t IT_0315 = IT_0308 + IT_0314;
    const ccomplex_t IT_0316 = (-2)*IT_0315;
    const ccomplex_t IT_0317 = IT_0297 + IT_0302 + IT_0316;
    const ccomplex_t IT_0318 = IT_0085 + IT_0239;
    const ccomplex_t IT_0319 = IT_0243 + IT_0318;
    const ccomplex_t IT_0320 = 2*conj(IT_0296);
    const ccomplex_t IT_0321 = (-2)*conj(IT_0301);
    const ccomplex_t IT_0322 = (-2)*conj(IT_0315);
    const ccomplex_t IT_0323 = IT_0320 + IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = conj(IT_0249) + conj(IT_0251);
    const ccomplex_t IT_0325 = m_N_4*IT_0169;
    const ccomplex_t IT_0326 = IT_0171*(IT_0292 + -IT_0325);
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = 2*IT_0328;
    const ccomplex_t IT_0330 = cpow(IT_0073, 2);
    const ccomplex_t IT_0331 = m_N_3*IT_0330;
    const ccomplex_t IT_0332 = -IT_0175*(IT_0311 + -IT_0331);
    const ccomplex_t IT_0333 = (-2)*IT_0332;
    const ccomplex_t IT_0334 = cpow(IT_0053, 2);
    const ccomplex_t IT_0335 = m_N_1*IT_0334;
    const ccomplex_t IT_0336 = -IT_0300;
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = IT_0178*IT_0337;
    const ccomplex_t IT_0339 = cpow(IT_0032, 2);
    const ccomplex_t IT_0340 = m_N_2*IT_0339;
    const ccomplex_t IT_0341 = IT_0306 + IT_0340;
    const ccomplex_t IT_0342 = IT_0181*IT_0341;
    const ccomplex_t IT_0343 = IT_0338 + IT_0342;
    const ccomplex_t IT_0344 = (-2)*IT_0343;
    const ccomplex_t IT_0345 = IT_0329 + IT_0333 + IT_0344;
    const ccomplex_t IT_0346 = IT_0249 + IT_0251;
    const ccomplex_t IT_0347 = 2*conj(IT_0328);
    const ccomplex_t IT_0348 = (-2)*conj(IT_0332);
    const ccomplex_t IT_0349 = (-2)*conj(IT_0343);
    const ccomplex_t IT_0350 = IT_0347 + IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = m_N_4*IT_0002;
    const ccomplex_t IT_0352 = s_14*IT_0279;
    const ccomplex_t IT_0353 = IT_0002*IT_0352;
    const ccomplex_t IT_0354 = (-2)*IT_0261 + (-0.5)*IT_0353;
    const ccomplex_t IT_0355 = IT_0259 + IT_0354;
    const ccomplex_t IT_0356 = IT_0351*IT_0355;
    const ccomplex_t IT_0357 = (-2)*IT_0356;
    const ccomplex_t IT_0358 = s_14*s_23;
    const ccomplex_t IT_0359 = s_13*s_24;
    const ccomplex_t IT_0360 = s_12*s_34;
    const ccomplex_t IT_0361 = -IT_0359 + -IT_0360;
    const ccomplex_t IT_0362 = IT_0358 + IT_0361;
    const ccomplex_t IT_0363 = IT_0225*IT_0362;
    const ccomplex_t IT_0364 = (-0.25)*IT_0363;
    const ccomplex_t IT_0365 = s_12*IT_0042;
    const ccomplex_t IT_0366 = IT_0002*IT_0365;
    const ccomplex_t IT_0367 = (-0.5)*IT_0366;
    const ccomplex_t IT_0368 = s_12 + IT_0364 + IT_0367;
    const ccomplex_t IT_0369 = IT_0241*(IT_0292 + -IT_0325);
    const ccomplex_t IT_0370 = -IT_0369;
    const ccomplex_t IT_0371 = -IT_0370;
    const ccomplex_t IT_0372 = 8*IT_0142;
    const ccomplex_t IT_0373 = (-0.125)*IT_0372;
    const ccomplex_t IT_0374 = IT_0147 + IT_0373;
    const ccomplex_t IT_0375 = IT_0299 + IT_0336;
    const ccomplex_t IT_0376 = IT_0064*IT_0375;
    const ccomplex_t IT_0377 = IT_0044*IT_0307;
    const ccomplex_t IT_0378 = IT_0084*IT_0313;
    const ccomplex_t IT_0379 = IT_0376 + IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = 8*IT_0147;
    const ccomplex_t IT_0381 = (-8)*IT_0142;
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = 8*conj(IT_0147);
    const ccomplex_t IT_0384 = (-8)*conj(IT_0142);
    const ccomplex_t IT_0385 = IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = 8*conj(IT_0142);
    const ccomplex_t IT_0387 = (-8)*conj(IT_0147);
    const ccomplex_t IT_0388 = IT_0386 + IT_0387;
    const ccomplex_t IT_0389 = -IT_0084*(IT_0311 + -IT_0331);
    const ccomplex_t IT_0390 = IT_0064*IT_0337;
    const ccomplex_t IT_0391 = IT_0044*IT_0341;
    const ccomplex_t IT_0392 = IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = conj(IT_0389) + conj(IT_0392);
    const ccomplex_t IT_0394 = 8*IT_0214;
    const ccomplex_t IT_0395 = (-8)*IT_0216;
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = IT_0389 + IT_0392;
    const ccomplex_t IT_0398 = 8*conj(IT_0214);
    const ccomplex_t IT_0399 = (-8)*conj(IT_0216);
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = IT_0241*(IT_0292 + -IT_0293);
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = 8*IT_0216;
    const ccomplex_t IT_0405 = (-8)*IT_0214;
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = 8*conj(IT_0216);
    const ccomplex_t IT_0408 = (-8)*conj(IT_0214);
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = conj(IT_0371)*IT_0374 + (-0.125)*conj(IT_0379)
      *IT_0382 + (-0.125)*IT_0379*IT_0385 + (-0.125)*IT_0371*IT_0388 + (-0.125)
      *IT_0393*IT_0396 + (-0.125)*IT_0397*IT_0400 + (-0.125)*conj(IT_0403)
      *IT_0406 + (-0.125)*IT_0403*IT_0409;
    const ccomplex_t IT_0411 = IT_0358 + IT_0360;
    const ccomplex_t IT_0412 = -IT_0359;
    const ccomplex_t IT_0413 = IT_0411 + IT_0412;
    const ccomplex_t IT_0414 = IT_0225*IT_0413;
    const ccomplex_t IT_0415 = 0.25*IT_0414;
    const ccomplex_t IT_0416 = s_12 + IT_0367 + IT_0415;
    const ccomplex_t IT_0417 = conj(IT_0301) + conj(IT_0315);
    const ccomplex_t IT_0418 = IT_0301 + IT_0315;
    const ccomplex_t IT_0419 = conj(IT_0332) + conj(IT_0343);
    const ccomplex_t IT_0420 = IT_0332 + IT_0343;
    const ccomplex_t IT_0421 = conj(IT_0296)*IT_0374 + (-0.125)*IT_0296
      *IT_0388 + (-0.125)*conj(IT_0328)*IT_0406 + (-0.125)*IT_0328*IT_0409 + (
      -0.125)*IT_0382*IT_0417 + (-0.125)*IT_0385*IT_0418 + (-0.125)*IT_0396
      *IT_0419 + (-0.125)*IT_0400*IT_0420;
    const ccomplex_t IT_0422 = s_12*s_14;
    const ccomplex_t IT_0423 = s_12*s_14*IT_0042;
    const ccomplex_t IT_0424 = s_24*IT_0008*IT_0042;
    const ccomplex_t IT_0425 = s_24*IT_0153;
    const ccomplex_t IT_0426 = s_12*s_13*s_34;
    const ccomplex_t IT_0427 = -IT_0424 + 2*IT_0425 + (-2)*IT_0426;
    const ccomplex_t IT_0428 = IT_0423 + IT_0427;
    const ccomplex_t IT_0429 = IT_0002*IT_0428;
    const ccomplex_t IT_0430 = (-0.5)*IT_0429;
    const ccomplex_t IT_0431 = s_24*IT_0008;
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = IT_0422 + IT_0430 + IT_0432;
    const ccomplex_t IT_0434 = (-4)*IT_0195;
    const ccomplex_t IT_0435 = (-4)*IT_0183;
    const ccomplex_t IT_0436 = IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = (-4)*conj(IT_0195);
    const ccomplex_t IT_0438 = (-4)*conj(IT_0183);
    const ccomplex_t IT_0439 = IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = (-4)*IT_0173;
    const ccomplex_t IT_0441 = IT_0435 + IT_0440;
    const ccomplex_t IT_0442 = (-4)*conj(IT_0173);
    const ccomplex_t IT_0443 = IT_0438 + IT_0442;
    const ccomplex_t IT_0444 = IT_0000*IT_0352;
    const ccomplex_t IT_0445 = 0.666666666666667*IT_0260 + (-0.166666666666667
      )*IT_0444;
    const ccomplex_t IT_0446 = s_14 + IT_0445;
    const ccomplex_t IT_0447 = m_N_4*IT_0446;
    const ccomplex_t IT_0448 = (-6)*IT_0447;
    const ccomplex_t IT_0449 = conj(IT_0238) + conj(IT_0244);
    const ccomplex_t IT_0450 = IT_0333 + IT_0344;
    const ccomplex_t IT_0451 = IT_0238 + IT_0244;
    const ccomplex_t IT_0452 = IT_0348 + IT_0349;
    const ccomplex_t IT_0453 = IT_0302 + IT_0316;
    const ccomplex_t IT_0454 = IT_0321 + IT_0322;
    const ccomplex_t IT_0455 = s_13*s_23;
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = IT_0365 + IT_0456;
    const ccomplex_t IT_0458 = IT_0002*IT_0457;
    const ccomplex_t IT_0459 = (s_12 + -IT_0004)*(s_12 + (-0.5)*IT_0458);
    const ccomplex_t IT_0460 = (-2)*IT_0459;
    const ccomplex_t IT_0461 = 2*IT_0086;
    const ccomplex_t IT_0462 = IT_0173 + IT_0195;
    const ccomplex_t IT_0463 = 2*conj(IT_0086);
    const ccomplex_t IT_0464 = IT_0008*IT_0042;
    const ccomplex_t IT_0465 = -IT_0464;
    const ccomplex_t IT_0466 = IT_0013 + IT_0465;
    const ccomplex_t IT_0467 = IT_0002*IT_0466;
    const ccomplex_t IT_0468 = -IT_0008;
    const ccomplex_t IT_0469 = IT_0014 + IT_0468;
    const ccomplex_t IT_0470 = (-4)*IT_0469;
    const ccomplex_t IT_0471 = IT_0467 + IT_0470;
    const ccomplex_t IT_0472 = m_N_4*IT_0471;
    const ccomplex_t IT_0473 = 2*IT_0250;
    const ccomplex_t IT_0474 = (-2)*IT_0473;
    const ccomplex_t IT_0475 = (-2)*IT_0086;
    const ccomplex_t IT_0476 = IT_0474 + IT_0475;
    const ccomplex_t IT_0477 = (-2)*conj(IT_0473);
    const ccomplex_t IT_0478 = (-2)*conj(IT_0086);
    const ccomplex_t IT_0479 = IT_0477 + IT_0478;
    const ccomplex_t IT_0480 = 2*IT_0242;
    const ccomplex_t IT_0481 = (-2)*IT_0480;
    const ccomplex_t IT_0482 = IT_0475 + IT_0481;
    const ccomplex_t IT_0483 = (-2)*conj(IT_0480);
    const ccomplex_t IT_0484 = IT_0478 + IT_0483;
    const ccomplex_t IT_0485 = (-10)*s_14;
    const ccomplex_t IT_0486 = IT_0444 + IT_0485;
    const ccomplex_t IT_0487 = m_N_4*IT_0486;
    const ccomplex_t IT_0488 = 2*IT_0194;
    const ccomplex_t IT_0489 = 2*conj(IT_0194);
    const ccomplex_t IT_0490 = 2*IT_0172;
    const ccomplex_t IT_0491 = 2*conj(IT_0172);
    const ccomplex_t IT_0492 = conj(IT_0252)*IT_0297 + IT_0252*IT_0320 +
       IT_0329*IT_0449 + conj(IT_0194)*IT_0450 + IT_0347*IT_0451 + IT_0194
      *IT_0452 + conj(IT_0172)*IT_0453 + IT_0172*IT_0454 + conj(IT_0328)*IT_0488
       + IT_0328*IT_0489 + conj(IT_0296)*IT_0490 + IT_0296*IT_0491;
    const ccomplex_t IT_0493 = (-0.666666666666667)*IT_0233 +
       0.166666666666667*IT_0284;
    const ccomplex_t IT_0494 = s_24 + IT_0493;
    const ccomplex_t IT_0495 = m_N_4*IT_0494;
    const ccomplex_t IT_0496 = 6*IT_0495;
    const ccomplex_t IT_0497 = conj(IT_0252)*IT_0345 + IT_0252*IT_0350 +
       IT_0317*IT_0449 + conj(IT_0172)*IT_0450 + IT_0323*IT_0451 + IT_0172
      *IT_0452 + conj(IT_0194)*IT_0453 + IT_0194*IT_0454 + conj(IT_0296)*IT_0488
       + IT_0296*IT_0489 + conj(IT_0328)*IT_0490 + IT_0328*IT_0491;
    const ccomplex_t IT_0498 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_0499 = 0.5*IT_0155 + (-0.5)*IT_0464 + -IT_0498;
    const ccomplex_t IT_0500 = IT_0153 + IT_0499;
    const ccomplex_t IT_0501 = IT_0002*IT_0500;
    const ccomplex_t IT_0502 = -IT_0156;
    const ccomplex_t IT_0503 = IT_0008 + IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = m_N_4*IT_0503;
    const ccomplex_t IT_0505 = 4*IT_0183;
    const ccomplex_t IT_0506 = 4*conj(IT_0183);
    const ccomplex_t IT_0507 = IT_0184*IT_0379 + IT_0188*conj(IT_0379) + 0.25
      *conj(IT_0403)*IT_0436 + 0.25*IT_0403*IT_0439 + 0.25*conj(IT_0371)*IT_0441
       + 0.25*IT_0371*IT_0443 + IT_0393*(IT_0195 + 0.25*IT_0505) + IT_0397*(conj
      (IT_0195) + 0.25*IT_0506);
    const ccomplex_t IT_0508 = IT_0008 + IT_0502;
    const ccomplex_t IT_0509 = IT_0008*IT_0508;
    const ccomplex_t IT_0510 = -IT_0509;
    const ccomplex_t IT_0511 = 6*conj(IT_0183);
    const ccomplex_t IT_0512 = 1.125*IT_0008;
    const ccomplex_t IT_0513 = (-16)*IT_0379;
    const ccomplex_t IT_0514 = (-16)*conj(IT_0379);
    const ccomplex_t IT_0515 = (-16)*IT_0315;
    const ccomplex_t IT_0516 = (-16)*conj(IT_0301);
    const ccomplex_t IT_0517 = (-16)*conj(IT_0315);
    const ccomplex_t IT_0518 = ((-16)*conj(IT_0296) + 16*conj(IT_0301) + 16
      *conj(IT_0315))*IT_0328 + ((-16)*IT_0296 + 16*IT_0301 + 16*IT_0315)*conj
      (IT_0328) + 16*conj(IT_0379)*IT_0403 + conj(IT_0371)*(16*IT_0389 + 16
      *IT_0392 + (-16)*IT_0403) + 16*IT_0379*conj(IT_0403) + IT_0371*(16*conj
      (IT_0389) + 16*conj(IT_0392) + (-16)*conj(IT_0403)) + IT_0393*IT_0513 +
       IT_0397*IT_0514 + IT_0419*(16*IT_0296 + (-16)*IT_0301 + IT_0515) +
       IT_0420*(16*conj(IT_0296) + IT_0516 + IT_0517);
    const ccomplex_t IT_0519 = s_14*s_24;
    const ccomplex_t IT_0520 = (-2)*IT_0519;
    const ccomplex_t IT_0521 = IT_0365 + IT_0520;
    const ccomplex_t IT_0522 = IT_0002*IT_0521;
    const ccomplex_t IT_0523 = (-0.5)*IT_0522;
    const ccomplex_t IT_0524 = (-0.5)*IT_0458;
    const ccomplex_t IT_0525 = s_12*IT_0279;
    const ccomplex_t IT_0526 = s_13*s_23*IT_0042;
    const ccomplex_t IT_0527 = s_14*s_24*IT_0042;
    const ccomplex_t IT_0528 = s_13*s_24*s_34;
    const ccomplex_t IT_0529 = (-2)*IT_0526 + (-2)*IT_0527 + 4*IT_0528;
    const ccomplex_t IT_0530 = IT_0525 + IT_0529;
    const ccomplex_t IT_0531 = IT_0000*IT_0530;
    const ccomplex_t IT_0532 = 0.25*IT_0531;
    const ccomplex_t IT_0533 = 8*IT_0371;
    const ccomplex_t IT_0534 = (-8)*IT_0379;
    const ccomplex_t IT_0535 = IT_0533 + IT_0534;
    const ccomplex_t IT_0536 = 8*IT_0403;
    const ccomplex_t IT_0537 = (-8)*IT_0389;
    const ccomplex_t IT_0538 = (-8)*IT_0392;
    const ccomplex_t IT_0539 = IT_0536 + IT_0537 + IT_0538;
    const ccomplex_t IT_0540 = (-8)*conj(IT_0389);
    const ccomplex_t IT_0541 = (-8)*conj(IT_0392);
    const ccomplex_t IT_0542 = IT_0540 + IT_0541;
    const ccomplex_t IT_0543 = 8*IT_0379;
    const ccomplex_t IT_0544 = 8*IT_0389;
    const ccomplex_t IT_0545 = 8*conj(IT_0389);
    const ccomplex_t IT_0546 = (-8)*conj(IT_0379);
    const ccomplex_t IT_0547 = (s_12 + IT_0523 + IT_0524 + IT_0532)*(conj
      (IT_0371)*IT_0535 + conj(IT_0403)*IT_0539 + IT_0403*IT_0542 + conj(IT_0379
      )*IT_0543 + IT_0393*IT_0544 + 8*IT_0392*(conj(IT_0392) + 0.125*IT_0545) +
       IT_0371*IT_0546);
    const ccomplex_t IT_0548 = s_14*s_23*s_34;
    const ccomplex_t IT_0549 = pow(s_34, 2);
    const ccomplex_t IT_0550 = s_12*IT_0549;
    const ccomplex_t IT_0551 = (-2)*IT_0526 + (-2)*IT_0527 + 2*IT_0528 + 2
      *IT_0548 + (-2)*IT_0550;
    const ccomplex_t IT_0552 = IT_0525 + IT_0551;
    const ccomplex_t IT_0553 = IT_0000*IT_0552;
    const ccomplex_t IT_0554 = (-0.25)*IT_0553;
    const ccomplex_t IT_0555 = -IT_0004;
    const ccomplex_t IT_0556 = -IT_0003;
    const ccomplex_t IT_0557 = 8*conj(IT_0371);
    const ccomplex_t IT_0558 = IT_0546 + IT_0557;
    const ccomplex_t IT_0559 = (-8)*conj(IT_0332);
    const ccomplex_t IT_0560 = (-8)*conj(IT_0343);
    const ccomplex_t IT_0561 = IT_0559 + IT_0560;
    const ccomplex_t IT_0562 = (-8)*conj(IT_0301);
    const ccomplex_t IT_0563 = (-8)*conj(IT_0315);
    const ccomplex_t IT_0564 = IT_0562 + IT_0563;
    const ccomplex_t IT_0565 = 8*IT_0332;
    const ccomplex_t IT_0566 = 8*conj(IT_0332);
    const ccomplex_t IT_0567 = 8*IT_0301;
    const ccomplex_t IT_0568 = 8*conj(IT_0301);
    const ccomplex_t IT_0569 = 8*conj(IT_0379);
    const ccomplex_t IT_0570 = 8*IT_0343;
    const ccomplex_t IT_0571 = (-8)*IT_0301;
    const ccomplex_t IT_0572 = (-8)*IT_0315;
    const ccomplex_t IT_0573 = (-8)*IT_0332;
    const ccomplex_t IT_0574 = (-8)*IT_0343;
    const ccomplex_t IT_0575 = conj(IT_0296)*IT_0535 + conj(IT_0328)*IT_0539 +
       IT_0328*IT_0542 + conj(IT_0315)*IT_0543 + 8*conj(IT_0343)*(IT_0392 +
       0.125*IT_0544) + IT_0343*IT_0545 + IT_0296*IT_0558 + IT_0403*IT_0561 +
       IT_0371*IT_0564 + IT_0393*IT_0565 + IT_0397*IT_0566 + conj(IT_0379)
      *IT_0567 + IT_0379*IT_0568 + IT_0315*IT_0569 + conj(IT_0392)*IT_0570 +
       conj(IT_0371)*(IT_0571 + IT_0572) + 8*conj(IT_0403)*(IT_0328 + 0.125
      *IT_0573 + 0.125*IT_0574);
    const ccomplex_t IT_0576 = (-2)*IT_0526 + (-2)*IT_0527 + 4*IT_0548;
    const ccomplex_t IT_0577 = IT_0525 + IT_0576;
    const ccomplex_t IT_0578 = IT_0000*IT_0577;
    const ccomplex_t IT_0579 = 0.25*IT_0578;
    const ccomplex_t IT_0580 = 8*IT_0296;
    const ccomplex_t IT_0581 = IT_0571 + IT_0572 + IT_0580;
    const ccomplex_t IT_0582 = IT_0573 + IT_0574;
    const ccomplex_t IT_0583 = 8*IT_0315;
    const ccomplex_t IT_0584 = 8*IT_0328*(conj(IT_0328) + 0.125*IT_0559 +
       0.125*IT_0560) + IT_0296*IT_0564 + IT_0419*IT_0565 + IT_0343*IT_0566 +
       IT_0417*IT_0567 + IT_0315*IT_0568 + conj(IT_0343)*IT_0570 + conj(IT_0296)
      *IT_0581 + conj(IT_0328)*IT_0582 + conj(IT_0315)*IT_0583;
    const ccomplex_t IT_0585 = (-2)*IT_0549;
    const ccomplex_t IT_0586 = IT_0279 + IT_0585;
    const ccomplex_t IT_0587 = IT_0000*IT_0586;
    const ccomplex_t IT_0588 = 0.25*IT_0587;
    const ccomplex_t IT_0589 = 1 + IT_0588;
    const ccomplex_t IT_0590 = IT_0008*IT_0589;
    const ccomplex_t IT_0591 = 8*conj(IT_0315);
    const ccomplex_t IT_0592 = conj(IT_0328)*IT_0535 + 8*conj(IT_0296)*
      (IT_0403 + 0.125*IT_0537 + 0.125*IT_0538) + IT_0296*IT_0542 + IT_0419
      *IT_0543 + IT_0328*IT_0558 + IT_0371*IT_0561 + IT_0403*IT_0564 + IT_0420
      *IT_0569 + conj(IT_0403)*IT_0581 + conj(IT_0371)*IT_0582 + IT_0393*
      (IT_0567 + IT_0583) + IT_0397*(IT_0568 + IT_0591);
    const ccomplex_t IT_0593 = pow(m_Z, 6);
    const ccomplex_t IT_0594 = s_12*IT_0593;
    const ccomplex_t IT_0595 = s_13*s_23*IT_0279;
    const ccomplex_t IT_0596 = (-2)*IT_0595;
    const ccomplex_t IT_0597 = IT_0594 + IT_0596;
    const ccomplex_t IT_0598 = IT_0000*IT_0597;
    const ccomplex_t IT_0599 = 0.25*IT_0598;
    const ccomplex_t IT_0600 = IT_0365 + IT_0599;
    const ccomplex_t IT_0601 = (-8)*IT_0172;
    const ccomplex_t IT_0602 = (-8)*IT_0252;
    const ccomplex_t IT_0603 = (-8)*conj(IT_0172);
    const ccomplex_t IT_0604 = (-8)*conj(IT_0252);
    const ccomplex_t IT_0605 = (-8)*IT_0194;
    const ccomplex_t IT_0606 = (-8)*IT_0244;
    const ccomplex_t IT_0607 = (-8)*conj(IT_0194);
    const ccomplex_t IT_0608 = (-8)*conj(IT_0244);
    const ccomplex_t IT_0609 = conj(IT_0239) + conj(IT_0243);
    const ccomplex_t IT_0610 = (-8)*IT_0238;
    const ccomplex_t IT_0611 = IT_0239 + IT_0243;
    const ccomplex_t IT_0612 = (-8)*conj(IT_0238);
    const ccomplex_t IT_0613 = IT_0324*(IT_0601 + IT_0602) + IT_0346*(IT_0603 
      + IT_0604) + (-8)*conj(IT_0085)*(IT_0238 + (-0.125)*IT_0605 + (-0.125)
      *IT_0606) + (-8)*IT_0085*(conj(IT_0238) + (-0.125)*IT_0607 + (-0.125)
      *IT_0608) + IT_0609*(IT_0605 + IT_0606 + IT_0610) + IT_0611*(IT_0607 +
       IT_0608 + IT_0612);
    const ccomplex_t IT_0614 = s_12 + IT_0555;
    const ccomplex_t IT_0615 = IT_0008*IT_0614;
    const ccomplex_t IT_0616 = (-6)*IT_0086;
    const ccomplex_t IT_0617 = (-6)*conj(IT_0086);
    const ccomplex_t IT_0618 = IT_0000*IT_0549;
    const ccomplex_t IT_0619 = 0.5*IT_0618;
    const ccomplex_t IT_0620 = 1 + IT_0619;
    const ccomplex_t IT_0621 = s_12*IT_0620;
    const ccomplex_t IT_0622 = (-4)*conj(IT_0147);
    const ccomplex_t IT_0623 = 0.25*IT_0622;
    const ccomplex_t IT_0624 = conj(IT_0142) + IT_0623;
    const ccomplex_t IT_0625 = (-4)*IT_0147;
    const ccomplex_t IT_0626 = (-4)*IT_0214;
    const ccomplex_t IT_0627 = (-4)*conj(IT_0214);
    const ccomplex_t IT_0628 = -IT_0014;
    const ccomplex_t IT_0629 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0630 = 0.5*IT_0013 + (-0.5)*IT_0464 + -IT_0629;
    const ccomplex_t IT_0631 = IT_0011 + IT_0630;
    const ccomplex_t IT_0632 = IT_0002*IT_0631;
    const ccomplex_t IT_0633 = IT_0008 + IT_0628 + IT_0632;
    const ccomplex_t IT_0634 = m_N_4*IT_0633;
    const ccomplex_t IT_0635 = (-4)*IT_0086;
    const ccomplex_t IT_0636 = (-4)*conj(IT_0086);
    const ccomplex_t IT_0637 = 0.5*IT_0598;
    const ccomplex_t IT_0638 = (-4)*IT_0519;
    const ccomplex_t IT_0639 = IT_0365 + IT_0637 + IT_0638;
    const ccomplex_t IT_0640 = IT_0002*IT_0577;
    const ccomplex_t IT_0641 = -IT_0640;
    const ccomplex_t IT_0642 = IT_0639 + IT_0641;
    const ccomplex_t IT_0643 = (-4)*IT_0194;
    const ccomplex_t IT_0644 = (-4)*conj(IT_0194);
    const ccomplex_t IT_0645 = (-4)*IT_0172;
    const ccomplex_t IT_0646 = (-4)*IT_0244;
    const ccomplex_t IT_0647 = IT_0002*IT_0530;
    const ccomplex_t IT_0648 = -IT_0647;
    const ccomplex_t IT_0649 = IT_0639 + IT_0648;
    const ccomplex_t IT_0650 = (-4)*IT_0249;
    const ccomplex_t IT_0651 = (-4)*IT_0085;
    const ccomplex_t IT_0652 = (-4)*conj(IT_0085);
    const ccomplex_t IT_0653 = (-4)*IT_0239;
    const ccomplex_t IT_0654 = IT_0000*IT_0280;
    const ccomplex_t IT_0655 = (-10)*s_24;
    const ccomplex_t IT_0656 = IT_0654 + IT_0655;
    const ccomplex_t IT_0657 = m_N_4*IT_0656;
    const ccomplex_t IT_0658 = 2*IT_0389;
    const ccomplex_t IT_0659 = 2*IT_0392;
    const ccomplex_t IT_0660 = IT_0658 + IT_0659;
    const ccomplex_t IT_0661 = (-2)*IT_0403;
    const ccomplex_t IT_0662 = IT_0660 + IT_0661;
    const ccomplex_t IT_0663 = 2*conj(IT_0389);
    const ccomplex_t IT_0664 = 2*conj(IT_0392);
    const ccomplex_t IT_0665 = IT_0663 + IT_0664;
    const ccomplex_t IT_0666 = (-2)*conj(IT_0403);
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = 2*IT_0379;
    const ccomplex_t IT_0669 = 2*conj(IT_0379);
    const ccomplex_t IT_0670 = (-2)*IT_0371;
    const ccomplex_t IT_0671 = (-2)*conj(IT_0371);
    const ccomplex_t IT_0672 = s_13*s_34*IT_0042;
    const ccomplex_t IT_0673 = (-2)*IT_0672;
    const ccomplex_t IT_0674 = IT_0352 + IT_0673;
    const ccomplex_t IT_0675 = IT_0000*IT_0674;
    const ccomplex_t IT_0676 = (-0.666666666666667)*IT_0264 +
       0.166666666666667*IT_0675;
    const ccomplex_t IT_0677 = s_14 + IT_0676;
    const ccomplex_t IT_0678 = m_N_4*IT_0677;
    const ccomplex_t IT_0679 = 6*IT_0678;
    const ccomplex_t IT_0680 = 0.666666666666667*IT_0229 + (-0.166666666666667
      )*IT_0654;
    const ccomplex_t IT_0681 = s_24 + IT_0680;
    const ccomplex_t IT_0682 = m_N_4*IT_0681;
    const ccomplex_t IT_0683 = (-6)*IT_0682;
    const ccomplex_t IT_0684 = -IT_0519;
    const ccomplex_t IT_0685 = IT_0365 + IT_0684;
    const ccomplex_t IT_0686 = IT_0002*IT_0685;
    const ccomplex_t IT_0687 = -IT_0004 + (-0.25)*IT_0686;
    const ccomplex_t IT_0688 = s_12 + IT_0687;
    const ccomplex_t IT_0689 = m_N_4*IT_0688;
    const ccomplex_t IT_0690 = (-4)*IT_0689;
    const ccomplex_t IT_0691 = 2*conj(IT_0332);
    const ccomplex_t IT_0692 = 2*conj(IT_0343);
    const ccomplex_t IT_0693 = IT_0691 + IT_0692;
    const ccomplex_t IT_0694 = 2*IT_0195;
    const ccomplex_t IT_0695 = 2*IT_0183;
    const ccomplex_t IT_0696 = IT_0694 + IT_0695;
    const ccomplex_t IT_0697 = 2*conj(IT_0195);
    const ccomplex_t IT_0698 = 2*conj(IT_0183);
    const ccomplex_t IT_0699 = IT_0697 + IT_0698;
    const ccomplex_t IT_0700 = 2*IT_0173;
    const ccomplex_t IT_0701 = IT_0695 + IT_0700;
    const ccomplex_t IT_0702 = 2*conj(IT_0173);
    const ccomplex_t IT_0703 = IT_0698 + IT_0702;
    const ccomplex_t IT_0704 = 2*IT_0301;
    const ccomplex_t IT_0705 = 2*IT_0315;
    const ccomplex_t IT_0706 = IT_0704 + IT_0705;
    const ccomplex_t IT_0707 = 2*conj(IT_0301);
    const ccomplex_t IT_0708 = 2*conj(IT_0315);
    const ccomplex_t IT_0709 = IT_0707 + IT_0708;
    const ccomplex_t IT_0710 = 2*IT_0343;
    const ccomplex_t IT_0711 = (-2)*IT_0379;
    const ccomplex_t IT_0712 = (-2)*conj(IT_0379);
    const ccomplex_t IT_0713 = (-2)*IT_0389;
    const ccomplex_t IT_0714 = (-2)*IT_0392;
    const ccomplex_t IT_0715 = (-2)*conj(IT_0389);
    const ccomplex_t IT_0716 = (-2)*conj(IT_0392);
    const ccomplex_t IT_0717 = (-2)*IT_0183;
    const ccomplex_t IT_0718 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0719 = (conj(IT_0301) + conj(IT_0315) + conj(IT_0332) 
      + conj(IT_0343))*IT_0461 + (IT_0301 + IT_0315 + IT_0332 + IT_0343)*IT_0463
       + conj(IT_0296)*IT_0476 + IT_0296*IT_0479 + conj(IT_0328)*IT_0482 +
       IT_0328*IT_0484 + IT_0480*IT_0693 + conj(IT_0371)*IT_0696 + IT_0371
      *IT_0699 + conj(IT_0403)*IT_0701 + IT_0403*IT_0703 + conj(IT_0473)*IT_0706
       + IT_0473*IT_0709 + 2*conj(IT_0480)*(IT_0332 + 0.5*IT_0710) + IT_0196
      *IT_0711 + (IT_0183 + IT_0195)*IT_0712 + conj(IT_0173)*(IT_0713 + IT_0714)
       + IT_0173*(IT_0715 + IT_0716) + IT_0393*IT_0717 + IT_0397*IT_0718;
    const ccomplex_t IT_0720 = (IT_0008 + -IT_0156)*(s_12 + (-0.5)*IT_0458);
    const ccomplex_t IT_0721 = 2*IT_0720;
    const ccomplex_t IT_0722 = (-2)*IT_0195;
    const ccomplex_t IT_0723 = IT_0717 + IT_0722;
    const ccomplex_t IT_0724 = (-2)*IT_0173;
    const ccomplex_t IT_0725 = IT_0717 + IT_0724;
    const ccomplex_t IT_0726 = IT_0155 + IT_0465;
    const ccomplex_t IT_0727 = IT_0002*IT_0726;
    const ccomplex_t IT_0728 = IT_0156 + IT_0468;
    const ccomplex_t IT_0729 = (-4)*IT_0728;
    const ccomplex_t IT_0730 = IT_0727 + IT_0729;
    const ccomplex_t IT_0731 = m_N_4*IT_0730;
    const ccomplex_t IT_0732 = conj(IT_0195)*IT_0450 + IT_0195*IT_0452 + conj
      (IT_0173)*IT_0453 + IT_0173*IT_0454 + conj(IT_0328)*IT_0696 + IT_0328
      *IT_0699 + conj(IT_0296)*IT_0701 + IT_0296*IT_0703 + IT_0417*IT_0717 +
       IT_0419*IT_0717 + IT_0418*IT_0718 + IT_0420*IT_0718;
    const ccomplex_t IT_0733 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0734 = (-2)*IT_0455 + -IT_0519 + 2*IT_0733;
    const ccomplex_t IT_0735 = IT_0365 + IT_0734;
    const ccomplex_t IT_0736 = IT_0002*IT_0735;
    const ccomplex_t IT_0737 = -IT_0004 + (-0.5)*IT_0736;
    const ccomplex_t IT_0738 = s_12 + IT_0737;
    const ccomplex_t IT_0739 = m_N_4*IT_0738;
    const ccomplex_t IT_0740 = (-2)*IT_0739;
    const ccomplex_t IT_0741 = (-2)*conj(IT_0173);
    const ccomplex_t IT_0742 = IT_0718 + IT_0741;
    const ccomplex_t IT_0743 = (-2)*conj(IT_0195);
    const ccomplex_t IT_0744 = IT_0173*IT_0693 + IT_0419*IT_0695 + IT_0343
      *IT_0698 + 2*IT_0332*(conj(IT_0173) + 0.5*IT_0698) + IT_0196*IT_0706 + 
      (IT_0183 + IT_0195)*IT_0709 + conj(IT_0173)*IT_0710 + conj(IT_0296)
      *IT_0723 + conj(IT_0328)*IT_0725 + IT_0328*IT_0742 + IT_0296*(IT_0718 +
       IT_0743);
    const ccomplex_t IT_0745 = (IT_0008 + -IT_0014)*(s_12 + (-0.5)*IT_0458);
    const ccomplex_t IT_0746 = 2*IT_0745;
    const ccomplex_t IT_0747 = conj(IT_0086)*IT_0475*IT_0746;
    const ccomplex_t IT_0748 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_0749 = (-2)*IT_0455 + -IT_0519 + 2*IT_0748;
    const ccomplex_t IT_0750 = IT_0365 + IT_0749;
    const ccomplex_t IT_0751 = IT_0002*IT_0750;
    const ccomplex_t IT_0752 = -IT_0004 + (-0.5)*IT_0751;
    const ccomplex_t IT_0753 = s_12 + IT_0752;
    const ccomplex_t IT_0754 = m_N_4*IT_0753;
    const ccomplex_t IT_0755 = (-2)*IT_0754;
    const ccomplex_t IT_0756 = s_13*s_14*s_23;
    const ccomplex_t IT_0757 = s_23*s_34*IT_0008;
    const ccomplex_t IT_0758 = -IT_0424 + (-2)*IT_0756 + 2*IT_0757;
    const ccomplex_t IT_0759 = IT_0423 + IT_0758;
    const ccomplex_t IT_0760 = IT_0002*IT_0759;
    const ccomplex_t IT_0761 = 2*IT_0431;
    const ccomplex_t IT_0762 = (-2)*IT_0422;
    const ccomplex_t IT_0763 = IT_0760 + IT_0761 + IT_0762;
    const ccomplex_t IT_0764 = conj(IT_0238) + IT_0245;
    const ccomplex_t IT_0765 = IT_0238 + IT_0247;
    const ccomplex_t IT_0766 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0767 = (-2)*conj(IT_0172);
    const ccomplex_t IT_0768 = IT_0173*conj(IT_0173) + IT_0195*conj(IT_0195) +
       conj(IT_0172)*IT_0717 + IT_0172*IT_0718 + IT_0449*IT_0722 + conj(IT_0252)
      *IT_0725 + IT_0252*IT_0742 + IT_0451*IT_0743 + IT_0717*IT_0764 + IT_0718
      *IT_0765 + IT_0195*IT_0766 + IT_0173*IT_0767;
    const ccomplex_t IT_0769 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0770 = (-0.166666666666667)*IT_0675 + (
      -0.666666666666667)*IT_0769;
    const ccomplex_t IT_0771 = s_14 + IT_0770;
    const ccomplex_t IT_0772 = m_N_4*IT_0771;
    const ccomplex_t IT_0773 = (-6)*IT_0772;
    const ccomplex_t IT_0774 = IT_0173*conj(IT_0371) + IT_0195*conj(IT_0403) +
       IT_0449*IT_0661 + IT_0451*IT_0666 + IT_0253*IT_0668 + IT_0255*IT_0669 +
       conj(IT_0252)*IT_0670 + IT_0252*IT_0671 + IT_0660*IT_0764 + IT_0665
      *IT_0765 + IT_0403*IT_0766 + IT_0371*IT_0767;
    const ccomplex_t IT_0775 = IT_0002*IT_0280;
    const ccomplex_t IT_0776 = (-2)*IT_0230 + (-0.5)*IT_0775;
    const ccomplex_t IT_0777 = IT_0228 + IT_0776;
    const ccomplex_t IT_0778 = IT_0351*IT_0777;
    const ccomplex_t IT_0779 = (-2)*IT_0778;
    const ccomplex_t IT_0780 = IT_0195*conj(IT_0371) + IT_0173*conj(IT_0403) +
       conj(IT_0172)*IT_0660 + conj(IT_0252)*IT_0662 + IT_0172*IT_0665 + IT_0252
      *IT_0667 + IT_0449*IT_0670 + IT_0451*IT_0671 + IT_0668*IT_0764 + IT_0669
      *IT_0765 + IT_0371*IT_0766 + IT_0403*IT_0767;
    const ccomplex_t IT_0781 = IT_0002*IT_0008;
    const ccomplex_t IT_0782 = (-0.25)*IT_0279;
    const ccomplex_t IT_0783 = IT_0549 + IT_0782;
    const ccomplex_t IT_0784 = IT_0781*IT_0783;
    const ccomplex_t IT_0785 = (-4)*IT_0784;
    const ccomplex_t IT_0786 = (-2)*IT_0252;
    const ccomplex_t IT_0787 = IT_0173 + IT_0786;
    const ccomplex_t IT_0788 = (-2)*conj(IT_0252);
    const ccomplex_t IT_0789 = IT_0767 + IT_0788;
    const ccomplex_t IT_0790 = (-2)*IT_0249;
    const ccomplex_t IT_0791 = (-2)*conj(IT_0249);
    const ccomplex_t IT_0792 = IT_0195*conj(IT_0249) + (-2)*conj(IT_0251)*((
      -0.5)*IT_0195 + IT_0451) + IT_0249*IT_0766 + (-2)*IT_0251*(IT_0449 + (-0.5
      )*IT_0766) + IT_0290*IT_0787 + IT_0319*IT_0789 + IT_0449*IT_0790 + IT_0451
      *IT_0791;
    const ccomplex_t IT_0793 = 9*IT_0464;
    const ccomplex_t IT_0794 = IT_0008*IT_0620;
    const ccomplex_t IT_0795 = (-4)*IT_0216;
    const ccomplex_t IT_0796 = (-4)*conj(IT_0216);
    const ccomplex_t IT_0797 = (-4)*IT_0142;
    const ccomplex_t IT_0798 = 0.25*IT_0797;
    const ccomplex_t IT_0799 = (-4)*conj(IT_0142);
    const ccomplex_t IT_0800 = 0.25*IT_0799;
    const ccomplex_t IT_0801 = (IT_0214 + IT_0332 + IT_0343 + IT_0389 +
       IT_0392)*IT_0624 + (conj(IT_0214) + conj(IT_0332) + conj(IT_0343) + conj
      (IT_0389) + conj(IT_0392))*(IT_0142 + 0.25*IT_0625) + 0.25*(conj(IT_0301) 
      + conj(IT_0315) + conj(IT_0379))*IT_0626 + 0.25*(IT_0301 + IT_0315 +
       IT_0379)*IT_0627 + (conj(IT_0296) + conj(IT_0371))*(IT_0214 + 0.25
      *IT_0795) + (IT_0296 + IT_0371)*(conj(IT_0214) + 0.25*IT_0796) + (conj
      (IT_0328) + conj(IT_0403))*(IT_0147 + IT_0798) + conj(IT_0216)*(IT_0147 +
       IT_0301 + IT_0315 + IT_0379 + IT_0798) + (IT_0328 + IT_0403)*(conj
      (IT_0147) + IT_0800) + IT_0216*(conj(IT_0147) + conj(IT_0301) + conj
      (IT_0315) + conj(IT_0379) + IT_0800);
    const ccomplex_t IT_0802 = s_12*s_24;
    const ccomplex_t IT_0803 = s_12*s_24*IT_0042;
    const ccomplex_t IT_0804 = s_14*IT_0008*IT_0042;
    const ccomplex_t IT_0805 = s_13*s_23*s_24;
    const ccomplex_t IT_0806 = s_13*s_34*IT_0008;
    const ccomplex_t IT_0807 = -IT_0804 + (-2)*IT_0805 + 2*IT_0806;
    const ccomplex_t IT_0808 = IT_0803 + IT_0807;
    const ccomplex_t IT_0809 = IT_0002*IT_0808;
    const ccomplex_t IT_0810 = (-0.5)*IT_0809;
    const ccomplex_t IT_0811 = s_14*IT_0008;
    const ccomplex_t IT_0812 = -IT_0811;
    const ccomplex_t IT_0813 = IT_0802 + IT_0810 + IT_0812;
    const ccomplex_t IT_0814 = s_14*IT_0011;
    const ccomplex_t IT_0815 = s_12*s_23*s_34;
    const ccomplex_t IT_0816 = 0.5*IT_0803 + (-0.5)*IT_0804 + -IT_0815;
    const ccomplex_t IT_0817 = IT_0814 + IT_0816;
    const ccomplex_t IT_0818 = IT_0002*IT_0817;
    const ccomplex_t IT_0819 = -IT_0802;
    const ccomplex_t IT_0820 = IT_0811 + IT_0818 + IT_0819;
    const ccomplex_t IT_0821 = IT_0007 + IT_0018;
    const ccomplex_t IT_0822 = IT_0008 + IT_0628;
    const ccomplex_t IT_0823 = IT_0008*IT_0822;
    const ccomplex_t IT_0824 = -IT_0823;
    const ccomplex_t IT_0825 = 6*conj(IT_0086);
    const ccomplex_t IT_0826 = (-6)*conj(IT_0183);
    const ccomplex_t IT_0827 = (-4)*conj(IT_0301);
    const ccomplex_t IT_0828 = 0.25*IT_0827;
    const ccomplex_t IT_0829 = (-4)*conj(IT_0315);
    const ccomplex_t IT_0830 = 0.25*IT_0829;
    const ccomplex_t IT_0831 = 0.25*IT_0472*IT_0669 + 0.5*IT_0460*IT_0698 + 
      (conj(IT_0085) + conj(IT_0239) + conj(IT_0243) + conj(IT_0249) + conj
      (IT_0251))*IT_0813 + (conj(IT_0172) + conj(IT_0194) + conj(IT_0238) + conj
      (IT_0244) + conj(IT_0252))*IT_0820 + 0.25*IT_0277*IT_0821 + 0.25*IT_0824
      *IT_0825 + 0.25*IT_0615*IT_0826 + IT_0634*(conj(IT_0296) + conj(IT_0328) +
       IT_0828 + IT_0830);
    const ccomplex_t IT_0832 = (-0.333333333333333)*IT_0755;
    const ccomplex_t IT_0833 = (-0.166666666666667)*IT_0007;
    const ccomplex_t IT_0834 = (-0.166666666666667)*IT_0018;
    const ccomplex_t IT_0835 = (-0.166666666666667)*IT_0460;
    const ccomplex_t IT_0836 = (-0.166666666666667)*IT_0746;
    const ccomplex_t IT_0837 = 6*IT_0086;
    const ccomplex_t IT_0838 = (-0.166666666666667)*IT_0824;
    const ccomplex_t IT_0839 = IT_0837*IT_0838;
    const ccomplex_t IT_0840 = (-0.166666666666667)*IT_0472;
    const ccomplex_t IT_0841 = (IT_0183 + IT_0195)*IT_0615 + (
      -0.666666666666667)*(IT_0328 + -IT_0420)*IT_0634 + (-0.666666666666667)
      *IT_0346*IT_0813 + (-0.666666666666667)*IT_0255*IT_0820 + (IT_0371 + 0.5
      *IT_0711)*IT_0832 + IT_0270*IT_0833 + IT_0274*IT_0834 + IT_0701*IT_0835 +
       IT_0476*IT_0836 + IT_0839 + IT_0660*IT_0840;
    const ccomplex_t IT_0842 = IT_0478*IT_0836;
    const ccomplex_t IT_0843 = (-1.5)*IT_0196*IT_0615 + (conj(IT_0328) + 
      -IT_0419)*IT_0634 + IT_0324*IT_0813 + IT_0253*IT_0820 + (-1.5)*(conj
      (IT_0371) + 0.5*IT_0712)*IT_0832 + (-1.5)*IT_0271*IT_0833 + (-1.5)*IT_0277
      *IT_0834 + (-1.5)*IT_0703*IT_0835 + (-1.5)*IT_0825*IT_0838 + (-1.5)
      *IT_0665*IT_0840 + (-1.5)*IT_0842;
    const ccomplex_t IT_0844 = -IT_0824;
    const ccomplex_t IT_0845 = (-4)*IT_0301;
    const ccomplex_t IT_0846 = 0.25*IT_0845;
    const ccomplex_t IT_0847 = (-4)*IT_0315;
    const ccomplex_t IT_0848 = 0.25*IT_0847;
    const ccomplex_t IT_0849 = (-0.666666666666667)*IT_0634;
    const ccomplex_t IT_0850 = IT_0188*IT_0615 + (-0.666666666666667)*IT_0319
      *IT_0813 + (-0.666666666666667)*IT_0765*IT_0820 + (IT_0403 + 0.5*IT_0713 +
       0.5*IT_0714)*IT_0832 + IT_0274*IT_0833 + IT_0270*IT_0834 + IT_0696
      *IT_0835 + IT_0482*IT_0836 + IT_0668*IT_0840 + (IT_0473 +
       0.166666666666667*IT_0837)*IT_0844 + (IT_0296 + IT_0846 + IT_0848)*IT_0849;
    const ccomplex_t IT_0851 = IT_0183*IT_0615 + (-0.666666666666667)*(IT_0085
       + IT_0239 + IT_0243 + IT_0249 + IT_0251)*IT_0813 + (-0.666666666666667)*
      (IT_0172 + IT_0194 + IT_0238 + IT_0244 + IT_0252)*IT_0820 + (
      -0.166666666666667)*(IT_0269 + IT_0272 + IT_0273)*IT_0821 + 2*IT_0695
      *IT_0835 + IT_0839 + IT_0668*IT_0840 + (IT_0296 + IT_0328 + IT_0846 +
       IT_0848)*IT_0849;
    const ccomplex_t IT_0852 = IT_0290*IT_0813 + IT_0764*IT_0820 + (-1.5)
      *IT_0615*(conj(IT_0173) + (-0.166666666666667)*IT_0826) + (-1.5)*(conj
      (IT_0403) + 0.5*IT_0715 + 0.5*IT_0716)*IT_0832 + (-1.5)*IT_0277*IT_0833 + 
      (-1.5)*IT_0271*IT_0834 + (-1.5)*IT_0699*IT_0835 + (-1.5)*IT_0669*IT_0840 +
       (-1.5)*IT_0842 + (-1.5)*(conj(IT_0473) + 0.166666666666667*IT_0825)
      *IT_0844 + (-1.5)*(conj(IT_0296) + IT_0828 + IT_0830)*IT_0849;
    const ccomplex_t IT_0853 = (IT_0007 + IT_0018)*IT_0150 + IT_0160*(IT_0184
      *IT_0187 + IT_0188*IT_0191 + IT_0196*IT_0218 + (IT_0183 + IT_0195)*IT_0221
      ) + (IT_0191*(IT_0183 + IT_0195) + IT_0187*IT_0196 + IT_0184*IT_0218 +
       IT_0188*IT_0221)*IT_0224 + IT_0237*(IT_0187*IT_0246 + IT_0191*IT_0248 +
       IT_0218*IT_0254 + IT_0221*IT_0256) + IT_0268*(IT_0254*IT_0270 + IT_0256
      *IT_0271 + IT_0246*IT_0274 + IT_0248*IT_0277) + IT_0288*(IT_0290*IT_0317 +
       IT_0319*IT_0323 + IT_0324*IT_0345 + IT_0346*IT_0350) + (IT_0317*IT_0324 +
       IT_0290*IT_0345 + IT_0323*IT_0346 + IT_0319*IT_0350)*IT_0357 + (-8)
      *IT_0368*IT_0410 + (-8)*IT_0416*IT_0421 + IT_0433*(IT_0290*IT_0436 +
       IT_0319*IT_0439 + IT_0324*IT_0441 + IT_0346*IT_0443) + IT_0448*(IT_0449
      *IT_0450 + IT_0451*IT_0452 + conj(IT_0252)*IT_0453 + IT_0252*IT_0454) +
       IT_0460*((conj(IT_0173) + conj(IT_0195))*IT_0461 + IT_0462*IT_0463) +
       IT_0472*(IT_0393*IT_0461 + IT_0397*IT_0463 + conj(IT_0403)*IT_0476 +
       IT_0403*IT_0479 + conj(IT_0371)*IT_0482 + IT_0371*IT_0484) + IT_0487
      *IT_0492 + IT_0496*IT_0497 + 4*IT_0504*IT_0507 + 6*IT_0510*(IT_0188
      *IT_0196 + (IT_0183 + IT_0195)*(conj(IT_0173) + 0.166666666666667*IT_0511)
      ) + IT_0512*IT_0518 + IT_0547 + (s_12 + IT_0554 + IT_0555 + IT_0556)
      *IT_0575 + (s_12 + IT_0523 + IT_0524 + IT_0579)*IT_0584 + IT_0590*IT_0592 
      + IT_0600*IT_0613 + IT_0615*((conj(IT_0173) + conj(IT_0183) + conj(IT_0195
      ))*IT_0616 + (IT_0183 + IT_0462)*IT_0617) + 4*IT_0621*(IT_0147*conj
      (IT_0147) + IT_0214*conj(IT_0214) + IT_0142*IT_0624 + 0.25*conj(IT_0142)
      *IT_0625 + 0.25*conj(IT_0216)*IT_0626 + IT_0216*(conj(IT_0216) + 0.25
      *IT_0627)) + IT_0634*(IT_0419*IT_0635 + IT_0420*IT_0636) + IT_0642*(conj
      (IT_0194)*IT_0643 + IT_0451*IT_0644 + (-4)*IT_0253*(IT_0252 + (-0.25)
      *IT_0645) + (-4)*IT_0449*(IT_0238 + (-0.25)*IT_0643 + (-0.25)*IT_0646)) + 
      (-4)*IT_0649*(IT_0324*(IT_0251 + (-0.25)*IT_0650) + (-0.25)*conj(IT_0085)
      *IT_0651 + (-0.25)*IT_0611*IT_0652 + IT_0609*(IT_0243 + (-0.25)*IT_0651 + 
      (-0.25)*IT_0653)) + IT_0657*(conj(IT_0251)*IT_0660 + conj(IT_0249)*IT_0662
       + IT_0251*IT_0665 + IT_0249*IT_0667 + IT_0290*IT_0668 + IT_0319*IT_0669 +
       IT_0289*IT_0670 + IT_0318*IT_0671) + (IT_0290*IT_0662 + IT_0319*IT_0667 +
       IT_0324*(IT_0668 + IT_0670) + IT_0346*(IT_0669 + IT_0671))*IT_0679 + 
      (conj(IT_0251)*IT_0661 + IT_0251*IT_0666 + conj(IT_0243)*IT_0670 + IT_0243
      *IT_0671)*IT_0683 + IT_0690*IT_0719 + 2*IT_0721*(conj(IT_0183)*IT_0717 +
       1./2*IT_0173*IT_0718 + 1./2*IT_0195*IT_0718 + 1./2*conj(IT_0195)*IT_0723 
      + 1./2*conj(IT_0173)*IT_0725) + IT_0731*IT_0732 + IT_0740*IT_0744 + 2
      *IT_0747 + ((conj(IT_0371) + conj(IT_0403))*IT_0461 + (IT_0371 + IT_0403)
      *IT_0463 + (conj(IT_0379) + conj(IT_0389) + conj(IT_0392))*IT_0475 + 
      (IT_0379 + IT_0389 + IT_0392)*IT_0478)*IT_0755 + IT_0763*IT_0768 + IT_0773
      *IT_0774 + IT_0779*IT_0780 + IT_0785*IT_0792 + (IT_0195*IT_0253 + IT_0255
      *IT_0766 + IT_0449*IT_0787 + IT_0451*IT_0789 + (-2)*IT_0290*(IT_0251 + (
      -0.5)*IT_0790) + (-2)*IT_0319*(conj(IT_0251) + (-0.5)*IT_0791))*IT_0793 +
       4*IT_0794*IT_0801 + 4*IT_0086*IT_0831 + (-6)*conj(IT_0473)*IT_0841 + 4
      *IT_0473*IT_0843 + (-6)*conj(IT_0480)*IT_0850 + (-6)*conj(IT_0086)*IT_0851
       + 4*IT_0480*IT_0852;
    return create_ccomplex_return(IT_0853);
}

