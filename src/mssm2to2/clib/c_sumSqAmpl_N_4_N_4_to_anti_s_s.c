#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_s_s(
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
    const creal_t m_N_4 = param->m_N_4;
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = m_s*m_N_4;
    const ccomplex_t IT_0001 = -s_24;
    const ccomplex_t IT_0002 = s_23 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B4)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_W4)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0009 + (-3)*IT_0014);
    const ccomplex_t IT_0016 = 0.166666666666667*IT_0015;
    const ccomplex_t IT_0017 = N_B4*e_em;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W4*e_em;
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + (-3)*IT_0022);
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = pow(m_s, 2);
    const ccomplex_t IT_0027 = pow(m_N_4, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0030 + IT_0041;
    const ccomplex_t IT_0043 = cos(alpha);
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = sin(alpha);
    const ccomplex_t IT_0048 = IT_0006*IT_0047;
    const ccomplex_t IT_0049 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = IT_0011*IT_0043;
    const ccomplex_t IT_0052 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = IT_0011*IT_0047;
    const ccomplex_t IT_0055 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0046 + IT_0050 + 
      -IT_0053 + -IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_s*e_em*IT_0011*IT_0031
      *IT_0033*IT_0047;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = IT_0044*IT_0049;
    const ccomplex_t IT_0066 = IT_0045*IT_0048;
    const ccomplex_t IT_0067 = IT_0051*IT_0055;
    const ccomplex_t IT_0068 = IT_0052*IT_0054;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + -IT_0066 + 
      -IT_0067 + IT_0068);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_s*e_em*IT_0011*IT_0031
      *IT_0033*IT_0043;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + (-2)*IT_0027 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0064 + 0.5*IT_0075;
    const ccomplex_t IT_0077 = 3*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0064 + (-0.5)*IT_0075;
    const ccomplex_t IT_0079 = (-3)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077 + IT_0079;
    const ccomplex_t IT_0081 = 3*conj(IT_0076);
    const ccomplex_t IT_0082 = (-3)*conj(IT_0078);
    const ccomplex_t IT_0083 = IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0034;
    const ccomplex_t IT_0085 = (-0.5)*IT_0036;
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = IT_0028*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0006;
    const ccomplex_t IT_0090 = (-0.333333333333333)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0006;
    const ccomplex_t IT_0092 = (-0.333333333333333)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = IT_0039*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0088 + IT_0095;
    const ccomplex_t IT_0097 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0098 = IT_0044*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0100 = IT_0048*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0102 = IT_0051*IT_0101;
    const ccomplex_t IT_0103 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0104 = IT_0054*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0098 + IT_0100 + 
      -IT_0102 + -IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0063*IT_0106;
    const ccomplex_t IT_0108 = IT_0044*IT_0099;
    const ccomplex_t IT_0109 = IT_0048*IT_0097;
    const ccomplex_t IT_0110 = IT_0051*IT_0103;
    const ccomplex_t IT_0111 = IT_0054*IT_0101;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + -IT_0109 + 
      -IT_0110 + IT_0111);
    const ccomplex_t IT_0113 = IT_0074*IT_0112;
    const ccomplex_t IT_0114 = 0.5*IT_0107 + 0.5*IT_0113;
    const ccomplex_t IT_0115 = 3*IT_0114;
    const ccomplex_t IT_0116 = (-0.5)*IT_0107 + (-0.5)*IT_0113;
    const ccomplex_t IT_0117 = (-3)*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = 3*conj(IT_0114);
    const ccomplex_t IT_0120 = (-3)*conj(IT_0116);
    const ccomplex_t IT_0121 = IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = -s_14;
    const ccomplex_t IT_0123 = s_13 + IT_0122;
    const ccomplex_t IT_0124 = IT_0000*IT_0123;
    const ccomplex_t IT_0125 = (-2)*IT_0124;
    const ccomplex_t IT_0126 = 4*IT_0003;
    const ccomplex_t IT_0127 = IT_0006*IT_0010;
    const ccomplex_t IT_0128 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = IT_0005*IT_0011;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = IT_0127*IT_0132;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0129 + IT_0131 + 
      -IT_0133 + -IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = e_em*IT_0127;
    const ccomplex_t IT_0141 = e_em*IT_0130;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0140 + 3*IT_0141);
    const ccomplex_t IT_0143 = (-0.166666666666667)*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0139*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0010;
    const ccomplex_t IT_0150 = 0.333333333333333*IT_0149;
    const ccomplex_t IT_0151 = IT_0144*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0139*IT_0152;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = conj(IT_0148) + conj(IT_0154);
    const ccomplex_t IT_0156 = IT_0148 + IT_0154;
    const ccomplex_t IT_0157 = (-0.5)*IT_0138;
    const ccomplex_t IT_0158 = IT_0146*IT_0157;
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0152*IT_0157;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = (-4)*IT_0124;
    const ccomplex_t IT_0163 = 4*IT_0124;
    const ccomplex_t IT_0164 = 0.5*IT_0153;
    const ccomplex_t IT_0165 = 0.5*IT_0147;
    const ccomplex_t IT_0166 = conj(IT_0164) + conj(IT_0165);
    const ccomplex_t IT_0167 = 3*IT_0078;
    const ccomplex_t IT_0168 = (-3)*IT_0076;
    const ccomplex_t IT_0169 = IT_0167 + IT_0168;
    const ccomplex_t IT_0170 = IT_0164 + IT_0165;
    const ccomplex_t IT_0171 = 3*conj(IT_0078);
    const ccomplex_t IT_0172 = (-3)*conj(IT_0076);
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0158;
    const ccomplex_t IT_0175 = 0.5*IT_0160;
    const ccomplex_t IT_0176 = 3*IT_0116;
    const ccomplex_t IT_0177 = (-3)*IT_0114;
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = 3*conj(IT_0116);
    const ccomplex_t IT_0180 = (-3)*conj(IT_0114);
    const ccomplex_t IT_0181 = IT_0179 + IT_0180;
    const ccomplex_t IT_0182 = (-2)*IT_0003;
    const ccomplex_t IT_0183 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0184 = IT_0086*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0187 = IT_0093*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = -IT_0185 + -IT_0188;
    const ccomplex_t IT_0190 = IT_0025*IT_0183;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0038*IT_0186;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = -IT_0191 + -IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0124;
    const ccomplex_t IT_0196 = (-4)*IT_0003;
    const ccomplex_t IT_0197 = s_12*IT_0026;
    const ccomplex_t IT_0198 = s_14*s_23;
    const ccomplex_t IT_0199 = s_12*s_34;
    const ccomplex_t IT_0200 = s_13*s_24;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0197 + IT_0198 + IT_0199 + IT_0201;
    const ccomplex_t IT_0203 = IT_0037*IT_0090;
    const ccomplex_t IT_0204 = IT_0039*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0024*IT_0085;
    const ccomplex_t IT_0207 = IT_0028*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = -IT_0205 + -IT_0208;
    const ccomplex_t IT_0210 = IT_0006*IT_0032;
    const ccomplex_t IT_0211 = IT_0045*IT_0210;
    const ccomplex_t IT_0212 = sin(beta);
    const ccomplex_t IT_0213 = IT_0006*IT_0212;
    const ccomplex_t IT_0214 = IT_0049*IT_0213;
    const ccomplex_t IT_0215 = IT_0011*IT_0032;
    const ccomplex_t IT_0216 = IT_0052*IT_0215;
    const ccomplex_t IT_0217 = IT_0011*IT_0212;
    const ccomplex_t IT_0218 = IT_0055*IT_0217;
    const ccomplex_t IT_0219 = -IT_0211 + IT_0214 + IT_0216 + -IT_0218;
    const ccomplex_t IT_0220 = m_s*e_em*IT_0011*IT_0031*IT_0033*IT_0212;
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = cpow((-2)*s_12 + (-2)*IT_0027 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0223 = IT_0221*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0219*IT_0224;
    const ccomplex_t IT_0226 = IT_0049*IT_0210;
    const ccomplex_t IT_0227 = IT_0045*IT_0213;
    const ccomplex_t IT_0228 = IT_0055*IT_0215;
    const ccomplex_t IT_0229 = IT_0052*IT_0217;
    const ccomplex_t IT_0230 = -IT_0226 + -IT_0227 + IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = m_s*e_em*IT_0011*IT_0031;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = pow(m_Z, 2);
    const ccomplex_t IT_0234 = cpow((-2)*s_12 + (-2)*IT_0027 + IT_0233 + 
      -reg_prop, -1);
    const ccomplex_t IT_0235 = IT_0232*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0230*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0225 + (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = 6*IT_0238;
    const ccomplex_t IT_0240 = 0.5*IT_0225 + 0.5*IT_0237;
    const ccomplex_t IT_0241 = (-6)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = 6*conj(IT_0238);
    const ccomplex_t IT_0244 = (-6)*conj(IT_0240);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0016*IT_0084;
    const ccomplex_t IT_0247 = IT_0028*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0035*IT_0092;
    const ccomplex_t IT_0250 = IT_0039*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = -IT_0248 + -IT_0251;
    const ccomplex_t IT_0253 = IT_0097*IT_0210;
    const ccomplex_t IT_0254 = IT_0099*IT_0213;
    const ccomplex_t IT_0255 = IT_0101*IT_0215;
    const ccomplex_t IT_0256 = IT_0103*IT_0217;
    const ccomplex_t IT_0257 = IT_0253 + -IT_0254 + -IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = IT_0224*IT_0257;
    const ccomplex_t IT_0259 = IT_0099*IT_0210;
    const ccomplex_t IT_0260 = IT_0097*IT_0213;
    const ccomplex_t IT_0261 = IT_0103*IT_0215;
    const ccomplex_t IT_0262 = IT_0101*IT_0217;
    const ccomplex_t IT_0263 = IT_0259 + IT_0260 + -IT_0261 + -IT_0262;
    const ccomplex_t IT_0264 = IT_0236*IT_0263;
    const ccomplex_t IT_0265 = 0.5*IT_0258 + 0.5*IT_0264;
    const ccomplex_t IT_0266 = 6*IT_0265;
    const ccomplex_t IT_0267 = (-0.5)*IT_0258 + (-0.5)*IT_0264;
    const ccomplex_t IT_0268 = (-6)*IT_0267;
    const ccomplex_t IT_0269 = IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = 6*conj(IT_0265);
    const ccomplex_t IT_0271 = (-6)*conj(IT_0267);
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = -IT_0198;
    const ccomplex_t IT_0274 = IT_0197 + IT_0199 + IT_0200 + IT_0273;
    const ccomplex_t IT_0275 = IT_0183*IT_0206;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0186*IT_0203;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = 6*IT_0240;
    const ccomplex_t IT_0281 = (-6)*IT_0238;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 6*conj(IT_0240);
    const ccomplex_t IT_0284 = (-6)*conj(IT_0238);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0183*IT_0246;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = IT_0186*IT_0249;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = 6*IT_0267;
    const ccomplex_t IT_0292 = (-6)*IT_0265;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = 6*conj(IT_0267);
    const ccomplex_t IT_0295 = (-6)*conj(IT_0265);
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_13 + s_14;
    const ccomplex_t IT_0298 = IT_0000*IT_0297;
    const ccomplex_t IT_0299 = conj(IT_0159) + conj(IT_0165);
    const ccomplex_t IT_0300 = 12*IT_0240;
    const ccomplex_t IT_0301 = (-12)*IT_0238;
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0159 + IT_0165;
    const ccomplex_t IT_0304 = 12*conj(IT_0240);
    const ccomplex_t IT_0305 = (-12)*conj(IT_0238);
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = conj(IT_0161) + conj(IT_0164);
    const ccomplex_t IT_0308 = 12*IT_0238;
    const ccomplex_t IT_0309 = (-12)*IT_0240;
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = IT_0161 + IT_0164;
    const ccomplex_t IT_0312 = 12*conj(IT_0238);
    const ccomplex_t IT_0313 = (-12)*conj(IT_0240);
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = conj(IT_0148) + conj(IT_0174);
    const ccomplex_t IT_0316 = 12*IT_0265;
    const ccomplex_t IT_0317 = (-12)*IT_0267;
    const ccomplex_t IT_0318 = IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = IT_0148 + IT_0174;
    const ccomplex_t IT_0320 = 12*conj(IT_0265);
    const ccomplex_t IT_0321 = (-12)*conj(IT_0267);
    const ccomplex_t IT_0322 = IT_0320 + IT_0321;
    const ccomplex_t IT_0323 = conj(IT_0154) + conj(IT_0175);
    const ccomplex_t IT_0324 = 12*IT_0267;
    const ccomplex_t IT_0325 = (-12)*IT_0265;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0154 + IT_0175;
    const ccomplex_t IT_0328 = 12*conj(IT_0267);
    const ccomplex_t IT_0329 = (-12)*conj(IT_0265);
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = conj(IT_0096)*IT_0242 + IT_0096*IT_0245 + conj
      (IT_0042)*IT_0269 + IT_0042*IT_0272 + conj(IT_0194)*IT_0282 + IT_0194
      *IT_0285 + conj(IT_0189)*IT_0293 + IT_0189*IT_0296 + IT_0299*IT_0302 +
       IT_0303*IT_0306 + IT_0307*IT_0310 + IT_0311*IT_0314 + IT_0315*IT_0318 +
       IT_0319*IT_0322 + IT_0323*IT_0326 + IT_0327*IT_0330;
    const ccomplex_t IT_0332 = 6*IT_0078;
    const ccomplex_t IT_0333 = (-6)*IT_0076;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = 6*conj(IT_0078);
    const ccomplex_t IT_0336 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0116;
    const ccomplex_t IT_0339 = (-6)*IT_0114;
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = 6*conj(IT_0116);
    const ccomplex_t IT_0342 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0343 = IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = -IT_0199;
    const ccomplex_t IT_0345 = s_23 + s_24;
    const ccomplex_t IT_0346 = IT_0000*IT_0345;
    const ccomplex_t IT_0347 = 6*IT_0096;
    const ccomplex_t IT_0348 = 6*conj(IT_0096);
    const ccomplex_t IT_0349 = 6*IT_0194;
    const ccomplex_t IT_0350 = 6*conj(IT_0194);
    const ccomplex_t IT_0351 = conj(IT_0042)*IT_0242 + IT_0042*IT_0245 + conj
      (IT_0189)*IT_0282 + IT_0189*IT_0285 + IT_0310*IT_0315 + IT_0307*IT_0318 +
       IT_0314*IT_0319 + IT_0311*IT_0322 + IT_0302*IT_0323 + IT_0299*IT_0326 +
       IT_0306*IT_0327 + IT_0303*IT_0330 + (-6)*conj(IT_0265)*(IT_0194 + (
      -0.166666666666667)*IT_0347) + (-6)*IT_0265*(conj(IT_0194) + (
      -0.166666666666667)*IT_0348) + (-6)*conj(IT_0267)*(IT_0096 + (
      -0.166666666666667)*IT_0349) + (-6)*IT_0267*(conj(IT_0096) + (
      -0.166666666666667)*IT_0350);
    const ccomplex_t IT_0352 = m_s*s_13*m_N_4;
    const ccomplex_t IT_0353 = 12*IT_0189;
    const ccomplex_t IT_0354 = 24*IT_0154;
    const ccomplex_t IT_0355 = 24*IT_0175;
    const ccomplex_t IT_0356 = IT_0353 + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = 12*conj(IT_0189);
    const ccomplex_t IT_0358 = 24*conj(IT_0154);
    const ccomplex_t IT_0359 = 24*conj(IT_0175);
    const ccomplex_t IT_0360 = IT_0357 + IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = 24*IT_0159;
    const ccomplex_t IT_0362 = 24*IT_0165;
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 24*conj(IT_0159);
    const ccomplex_t IT_0365 = 24*conj(IT_0165);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = 6*IT_0209;
    const ccomplex_t IT_0368 = 12*IT_0279;
    const ccomplex_t IT_0369 = IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = 6*conj(IT_0209);
    const ccomplex_t IT_0371 = 12*conj(IT_0279);
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 12*IT_0209;
    const ccomplex_t IT_0374 = 12*conj(IT_0209);
    const ccomplex_t IT_0375 = 12*IT_0252;
    const ccomplex_t IT_0376 = 12*conj(IT_0252);
    const ccomplex_t IT_0377 = conj(IT_0189)*IT_0252 + IT_0189*conj(IT_0252) +
       0.166666666666667*conj(IT_0290)*IT_0356 + 0.166666666666667*IT_0290
      *IT_0360 + 0.166666666666667*conj(IT_0279)*IT_0363 + 0.166666666666667
      *IT_0279*IT_0366 + 0.166666666666667*conj(IT_0194)*IT_0369 +
       0.166666666666667*IT_0194*IT_0372 + 0.166666666666667*IT_0299*IT_0373 +
       0.166666666666667*IT_0303*IT_0374 + 0.166666666666667*IT_0323*IT_0375 +
       0.166666666666667*IT_0327*IT_0376;
    const ccomplex_t IT_0378 = m_s*s_24*m_N_4;
    const ccomplex_t IT_0379 = 12*IT_0194;
    const ccomplex_t IT_0380 = IT_0363 + IT_0379;
    const ccomplex_t IT_0381 = 12*conj(IT_0194);
    const ccomplex_t IT_0382 = IT_0366 + IT_0381;
    const ccomplex_t IT_0383 = conj(IT_0252)*IT_0349 + IT_0252*IT_0350 + conj
      (IT_0189)*IT_0369 + IT_0189*IT_0372 + IT_0323*(24*IT_0279 + IT_0373) +
       IT_0327*(24*conj(IT_0279) + IT_0374) + IT_0299*IT_0375 + IT_0303*IT_0376 
      + conj(IT_0290)*IT_0380 + IT_0290*IT_0382;
    const ccomplex_t IT_0384 = m_s*s_23*m_N_4;
    const ccomplex_t IT_0385 = 12*IT_0096;
    const ccomplex_t IT_0386 = 24*IT_0161;
    const ccomplex_t IT_0387 = 24*IT_0164;
    const ccomplex_t IT_0388 = IT_0385 + IT_0386 + IT_0387;
    const ccomplex_t IT_0389 = 12*conj(IT_0096);
    const ccomplex_t IT_0390 = 24*conj(IT_0161);
    const ccomplex_t IT_0391 = 24*conj(IT_0164);
    const ccomplex_t IT_0392 = IT_0389 + IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = 24*IT_0209;
    const ccomplex_t IT_0394 = IT_0368 + IT_0393;
    const ccomplex_t IT_0395 = 24*conj(IT_0209);
    const ccomplex_t IT_0396 = IT_0371 + IT_0395;
    const ccomplex_t IT_0397 = 6*IT_0279;
    const ccomplex_t IT_0398 = IT_0373 + IT_0397;
    const ccomplex_t IT_0399 = 6*conj(IT_0279);
    const ccomplex_t IT_0400 = IT_0374 + IT_0399;
    const ccomplex_t IT_0401 = 12*IT_0290;
    const ccomplex_t IT_0402 = 12*conj(IT_0290);
    const ccomplex_t IT_0403 = IT_0384*(conj(IT_0290)*IT_0347 + IT_0290
      *IT_0348 + conj(IT_0252)*IT_0388 + IT_0252*IT_0392 + IT_0315*IT_0394 +
       IT_0319*IT_0396 + conj(IT_0042)*IT_0398 + IT_0042*IT_0400 + IT_0307
      *IT_0401 + IT_0311*IT_0402);
    const ccomplex_t IT_0404 = m_s*s_14*m_N_4;
    const ccomplex_t IT_0405 = 6*IT_0042;
    const ccomplex_t IT_0406 = 6*conj(IT_0042);
    const ccomplex_t IT_0407 = 12*conj(IT_0042);
    const ccomplex_t IT_0408 = 12*IT_0042;
    const ccomplex_t IT_0409 = 24*IT_0148;
    const ccomplex_t IT_0410 = IT_0307*IT_0394 + IT_0311*IT_0396 + conj
      (IT_0096)*IT_0398 + IT_0096*IT_0400 + IT_0315*IT_0401 + IT_0319*IT_0402 +
       conj(IT_0290)*IT_0405 + IT_0290*IT_0406 + IT_0252*(24*conj(IT_0148) + 24
      *conj(IT_0174) + IT_0407) + conj(IT_0252)*(24*IT_0174 + IT_0408 + IT_0409);
    const ccomplex_t IT_0411 = s_34*IT_0027;
    const ccomplex_t IT_0412 = IT_0361 + IT_0362 + IT_0379;
    const ccomplex_t IT_0413 = IT_0364 + IT_0365 + IT_0381;
    const ccomplex_t IT_0414 = 24*conj(IT_0148);
    const ccomplex_t IT_0415 = conj(IT_0189)*IT_0347 + IT_0189*IT_0348 + conj
      (IT_0042)*IT_0349 + IT_0042*IT_0350 + IT_0307*IT_0353 + IT_0311*IT_0357 +
       conj(IT_0148)*(IT_0361 + IT_0379) + IT_0148*(IT_0364 + IT_0381) + IT_0323
      *IT_0388 + IT_0327*IT_0392 + IT_0303*IT_0407 + IT_0299*IT_0408 + conj
      (IT_0165)*IT_0409 + conj(IT_0174)*IT_0412 + IT_0174*IT_0413 + IT_0165
      *IT_0414;
    const ccomplex_t IT_0416 = conj(IT_0194)*IT_0347 + IT_0194*IT_0348 + conj
      (IT_0148)*IT_0353 + conj(IT_0174)*IT_0356 + IT_0148*IT_0357 + IT_0174
      *IT_0360 + IT_0299*IT_0385 + IT_0303*IT_0389 + conj(IT_0189)*IT_0405 +
       IT_0189*IT_0406 + IT_0323*(IT_0408 + IT_0409) + IT_0307*IT_0412 + IT_0311
      *IT_0413 + IT_0327*(IT_0407 + IT_0414);
    const ccomplex_t IT_0417 = -s_34;
    const ccomplex_t IT_0418 = IT_0026 + IT_0417;
    const ccomplex_t IT_0419 = s_12*IT_0418;
    const ccomplex_t IT_0420 = (-24)*IT_0078;
    const ccomplex_t IT_0421 = IT_0419*(24*conj(IT_0076)*IT_0078 + IT_0076*((
      -24)*conj(IT_0076) + 24*conj(IT_0078)) + IT_0114*((-24)*conj(IT_0114) + 24
      *conj(IT_0116)) + IT_0116*(24*conj(IT_0114) + (-24)*conj(IT_0116)) + conj
      (IT_0078)*IT_0420);
    const ccomplex_t IT_0422 = IT_0027*IT_0418;
    const ccomplex_t IT_0423 = 6*IT_0076;
    const ccomplex_t IT_0424 = (-6)*IT_0078;
    const ccomplex_t IT_0425 = IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = 6*conj(IT_0076);
    const ccomplex_t IT_0427 = (-6)*conj(IT_0078);
    const ccomplex_t IT_0428 = IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = 6*IT_0114;
    const ccomplex_t IT_0430 = (-6)*IT_0116;
    const ccomplex_t IT_0431 = IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = (-6)*conj(IT_0116);
    const ccomplex_t IT_0433 = conj(IT_0078)*((-24)*IT_0114 + 24*IT_0116) +
       conj(IT_0076)*(24*IT_0114 + (-24)*IT_0116) + ((-24)*IT_0076 + 24*IT_0078)
      *conj(IT_0116) + conj(IT_0290)*IT_0334 + IT_0290*IT_0337 + conj(IT_0279)
      *IT_0340 + IT_0279*IT_0343 + conj(IT_0114)*(24*IT_0076 + 6*IT_0209 +
       IT_0420) + conj(IT_0252)*IT_0425 + IT_0252*IT_0428 + conj(IT_0209)
      *IT_0431 + IT_0209*IT_0432;
    const ccomplex_t IT_0434 = s_34 + IT_0026;
    const ccomplex_t IT_0435 = s_12*IT_0434;
    const ccomplex_t IT_0436 = (-24)*conj(IT_0238);
    const ccomplex_t IT_0437 = 24*conj(IT_0240);
    const ccomplex_t IT_0438 = IT_0436 + IT_0437;
    const ccomplex_t IT_0439 = (-24)*conj(IT_0240);
    const ccomplex_t IT_0440 = 24*conj(IT_0238);
    const ccomplex_t IT_0441 = IT_0027*IT_0434;
    const ccomplex_t IT_0442 = IT_0245*IT_0252 + IT_0242*conj(IT_0252) +
       IT_0285*IT_0290 + IT_0282*conj(IT_0290) + conj(IT_0265)*(24*IT_0238 + (
      -24)*IT_0240 + (-6)*IT_0279 + IT_0367) + (-6)*conj(IT_0267)*(IT_0209 + 4
      *IT_0238 + (-4)*IT_0240 + -1./6*IT_0397) + (-6)*IT_0267*(conj(IT_0209) + (
      -0.166666666666667)*IT_0399 + (-0.166666666666667)*IT_0438) + (-6)*IT_0265
      *(conj(IT_0279) + (-0.166666666666667)*IT_0370 + (-0.166666666666667)
      *IT_0439 + (-0.166666666666667)*IT_0440);
    const ccomplex_t IT_0443 = 48*conj(IT_0161);
    const ccomplex_t IT_0444 = 24*conj(IT_0096);
    const ccomplex_t IT_0445 = 48*conj(IT_0164);
    const ccomplex_t IT_0446 = IT_0443 + IT_0444 + IT_0445;
    const ccomplex_t IT_0447 = 48*IT_0161;
    const ccomplex_t IT_0448 = 24*IT_0096;
    const ccomplex_t IT_0449 = IT_0447 + IT_0448;
    const ccomplex_t IT_0450 = 48*IT_0148;
    const ccomplex_t IT_0451 = 48*conj(IT_0148);
    const ccomplex_t IT_0452 = conj(IT_0209)*IT_0373 + conj(IT_0252)*IT_0375 +
       conj(IT_0096)*(IT_0385 + IT_0386) + conj(IT_0042)*(IT_0408 + IT_0409) +
       IT_0164*IT_0446 + IT_0307*IT_0449 + IT_0315*(24*IT_0042 + IT_0450) +
       IT_0174*(24*conj(IT_0042) + 48*conj(IT_0174) + IT_0451);
    const ccomplex_t IT_0453 = IT_0026*IT_0027;
    const ccomplex_t IT_0454 = 48*IT_0159;
    const ccomplex_t IT_0455 = 24*IT_0194;
    const ccomplex_t IT_0456 = IT_0454 + IT_0455;
    const ccomplex_t IT_0457 = 48*conj(IT_0159);
    const ccomplex_t IT_0458 = 48*IT_0175;
    const ccomplex_t IT_0459 = 48*IT_0154;
    const ccomplex_t IT_0460 = 48*conj(IT_0154);
    const ccomplex_t IT_0461 = 48*IT_0165;
    const ccomplex_t IT_0462 = conj(IT_0290)*IT_0369 + IT_0290*IT_0372 + conj
      (IT_0189)*IT_0380 + IT_0189*IT_0382 + conj(IT_0042)*IT_0388 + IT_0042
      *IT_0392 + conj(IT_0252)*IT_0398 + IT_0252*IT_0400 + IT_0148*(IT_0443 +
       IT_0444) + IT_0174*IT_0446 + conj(IT_0174)*(48*IT_0164 + IT_0447 +
       IT_0448) + conj(IT_0148)*IT_0449 + conj(IT_0164)*IT_0450 + IT_0164
      *IT_0451 + conj(IT_0154)*IT_0456 + IT_0327*(24*conj(IT_0194) + IT_0457) +
       conj(IT_0165)*(IT_0458 + IT_0459) + IT_0165*IT_0460 + conj(IT_0175)*
      (IT_0454 + IT_0455 + IT_0461);
    const ccomplex_t IT_0463 = 24*IT_0189;
    const ccomplex_t IT_0464 = conj(IT_0189)*IT_0356 + conj(IT_0279)*IT_0368 +
       conj(IT_0194)*IT_0380 + conj(IT_0290)*IT_0401 + conj(IT_0159)*IT_0456 +
       IT_0165*IT_0457 + IT_0175*IT_0460 + conj(IT_0165)*(IT_0454 + IT_0455 +
       IT_0461) + conj(IT_0154)*(IT_0459 + IT_0463) + conj(IT_0175)*(IT_0458 +
       IT_0459 + IT_0463);
    const ccomplex_t IT_0465 = 0.166666666666667*conj(IT_0279);
    const ccomplex_t IT_0466 = IT_0198 + IT_0200 + IT_0344;
    const ccomplex_t IT_0467 = IT_0004*(conj(IT_0042)*IT_0080 + IT_0042
      *IT_0083 + conj(IT_0096)*IT_0118 + IT_0096*IT_0121) + (IT_0083*IT_0096 +
       IT_0080*conj(IT_0096) + conj(IT_0042)*IT_0118 + IT_0042*IT_0121)*IT_0125 
      + IT_0126*(IT_0080*IT_0155 + IT_0083*IT_0156 + IT_0121*(IT_0159 + IT_0161)
       + IT_0118*(conj(IT_0159) + conj(IT_0161))) + (IT_0118*IT_0155 + IT_0121
      *IT_0156 + IT_0083*(IT_0159 + IT_0161) + IT_0080*(conj(IT_0159) + conj
      (IT_0161)))*IT_0162 + IT_0163*(IT_0166*IT_0169 + IT_0170*IT_0173 + (conj
      (IT_0174) + conj(IT_0175))*IT_0178 + (IT_0174 + IT_0175)*IT_0181) +
       IT_0182*(IT_0173*IT_0189 + IT_0169*conj(IT_0189) + IT_0181*IT_0194 +
       IT_0178*conj(IT_0194)) + (IT_0181*IT_0189 + IT_0178*conj(IT_0189) +
       IT_0173*IT_0194 + IT_0169*conj(IT_0194))*IT_0195 + (IT_0173*(IT_0174 +
       IT_0175) + IT_0169*(conj(IT_0174) + conj(IT_0175)) + IT_0166*IT_0178 +
       IT_0170*IT_0181)*IT_0196 + IT_0202*(conj(IT_0209)*IT_0242 + IT_0209
      *IT_0245 + conj(IT_0252)*IT_0269 + IT_0252*IT_0272) + IT_0274*(conj
      (IT_0279)*IT_0282 + IT_0279*IT_0285 + conj(IT_0290)*IT_0293 + IT_0290
      *IT_0296) + IT_0298*IT_0331 + (conj(IT_0209)*IT_0334 + IT_0209*IT_0337 +
       conj(IT_0252)*IT_0340 + IT_0252*IT_0343)*(IT_0197 + IT_0200 + IT_0273 +
       IT_0344) + IT_0346*IT_0351 + 6*IT_0352*IT_0377 + IT_0378*IT_0383 +
       IT_0403 + IT_0404*IT_0410 + IT_0411*IT_0415 + IT_0197*IT_0416 + IT_0421 +
       IT_0422*IT_0433 + IT_0435*(IT_0267*((-24)*conj(IT_0265) + 24*conj(IT_0267
      )) + IT_0265*(24*conj(IT_0265) + (-24)*conj(IT_0267)) + IT_0240*IT_0438 +
       IT_0238*(IT_0439 + IT_0440)) + IT_0441*IT_0442 + IT_0198*IT_0452 +
       IT_0453*IT_0462 + IT_0200*IT_0464 + (IT_0197 + IT_0198 + IT_0201 +
       IT_0344)*(IT_0279*IT_0428 + conj(IT_0290)*IT_0431 + 6*IT_0290*(conj
      (IT_0114) + 0.166666666666667*IT_0432) + 6*IT_0425*IT_0465) + 6*(conj
      (IT_0252)*IT_0290 + IT_0252*conj(IT_0290) + 0.166666666666667*IT_0279
      *IT_0370 + IT_0367*IT_0465)*IT_0466;
    return create_ccomplex_return(IT_0467);
}

