#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_3_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_3_to_anti_s_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_N_3, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 2*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = IT_0019 + IT_0022;
    const ccomplex_t IT_0024 = IT_0013 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = IT_0006*IT_0031;
    const ccomplex_t IT_0033 = IT_0011*IT_0032;
    const ccomplex_t IT_0034 = 2*IT_0033;
    const ccomplex_t IT_0035 = IT_0015*IT_0031;
    const ccomplex_t IT_0036 = IT_0017*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = IT_0020*IT_0032;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = IT_0034 + IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0009, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0004*IT_0006*IT_0031*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = IT_0030*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0029 + IT_0049;
    const ccomplex_t IT_0051 = s_12*s_34;
    const ccomplex_t IT_0052 = s_14*s_23;
    const ccomplex_t IT_0053 = s_13*s_24;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0051 + IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0004*IT_0006*IT_0044;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = N_B3*e_em;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = N_W3*e_em;
    const ccomplex_t IT_0062 = IT_0006*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + (-3)*IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = IT_0057*IT_0065;
    const ccomplex_t IT_0067 = pow(m_s, 2);
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0002 + IT_0067 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 6*conj(IT_0070);
    const ccomplex_t IT_0072 = m_c*s_23*m_C_1;
    const ccomplex_t IT_0073 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0074 = IT_0006*IT_0073;
    const ccomplex_t IT_0075 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0076 = IT_0006*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + 0.5*IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0006;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow(IT_0031, -1);
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0004*IT_0006*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0057*IT_0088;
    const ccomplex_t IT_0090 = IT_0068*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0004*IT_0006*IT_0086;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B3)*e_em;
    const ccomplex_t IT_0096 = IT_0015*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W3)*e_em;
    const ccomplex_t IT_0099 = IT_0006*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + 3*IT_0100);
    const ccomplex_t IT_0102 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0103 = IT_0094*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_13 + IT_0001 + IT_0067 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 6*conj(IT_0106);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0109 = IT_0065*IT_0108;
    const ccomplex_t IT_0110 = IT_0068*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0004*IT_0006*IT_0086;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0094*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_23 + IT_0002 + IT_0067 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = -IT_0111 + -IT_0117;
    const ccomplex_t IT_0119 = 12*conj(IT_0070);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0015;
    const ccomplex_t IT_0121 = 0.666666666666667*IT_0120;
    const ccomplex_t IT_0122 = IT_0057*IT_0121;
    const ccomplex_t IT_0123 = cpow((-2)*s_13 + IT_0001 + IT_0067 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0006;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0004*IT_0006*IT_0044;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = IT_0104*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0125 + IT_0132;
    const ccomplex_t IT_0134 = IT_0050 + IT_0133;
    const ccomplex_t IT_0135 = 12*conj(IT_0085);
    const ccomplex_t IT_0136 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0137 = IT_0010*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0140 = IT_0016*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = IT_0138 + IT_0141;
    const ccomplex_t IT_0143 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0144 = IT_0010*IT_0143;
    const ccomplex_t IT_0145 = (-2)*IT_0144;
    const ccomplex_t IT_0146 = IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0046*IT_0148;
    const ccomplex_t IT_0150 = IT_0030*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0032*IT_0143;
    const ccomplex_t IT_0153 = (-2)*IT_0152;
    const ccomplex_t IT_0154 = IT_0035*IT_0139;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = IT_0032*IT_0136;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = IT_0155 + IT_0157;
    const ccomplex_t IT_0159 = IT_0153 + IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0008*IT_0161;
    const ccomplex_t IT_0163 = IT_0003*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0151 + -IT_0164;
    const ccomplex_t IT_0166 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0167 = IT_0006*IT_0166;
    const ccomplex_t IT_0168 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0169 = IT_0006*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0167 + (-0.5)*IT_0170);
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0081*IT_0172;
    const ccomplex_t IT_0174 = IT_0083*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = 12*conj(IT_0175);
    const ccomplex_t IT_0177 = 12*IT_0070;
    const ccomplex_t IT_0178 = 24*conj(IT_0085);
    const ccomplex_t IT_0179 = m_s*s_14*m_N_3;
    const ccomplex_t IT_0180 = 12*IT_0106;
    const ccomplex_t IT_0181 = 12*conj(IT_0106);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0004*IT_0006*IT_0009*IT_0086;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = IT_0148*IT_0183;
    const ccomplex_t IT_0185 = IT_0030*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0004*IT_0006;
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = IT_0161*IT_0188;
    const ccomplex_t IT_0190 = IT_0003*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = -IT_0186 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0026*IT_0188;
    const ccomplex_t IT_0194 = IT_0003*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0043*IT_0183;
    const ccomplex_t IT_0197 = IT_0030*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = IT_0195 + IT_0198;
    const ccomplex_t IT_0200 = IT_0088*IT_0108;
    const ccomplex_t IT_0201 = IT_0068*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0015;
    const ccomplex_t IT_0204 = (-0.333333333333333)*IT_0203;
    const ccomplex_t IT_0205 = IT_0094*IT_0204;
    const ccomplex_t IT_0206 = IT_0115*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0202 + IT_0207;
    const ccomplex_t IT_0209 = IT_0102*IT_0127;
    const ccomplex_t IT_0210 = IT_0104*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0004*IT_0006*IT_0044;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = IT_0057*IT_0213;
    const ccomplex_t IT_0215 = IT_0123*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = -IT_0211 + -IT_0216;
    const ccomplex_t IT_0218 = 24*conj(IT_0175);
    const ccomplex_t IT_0219 = m_c*s_13*m_N_3;
    const ccomplex_t IT_0220 = IT_0094*IT_0129;
    const ccomplex_t IT_0221 = IT_0104*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = (-12)*conj(IT_0085);
    const ccomplex_t IT_0225 = (-12)*conj(IT_0175);
    const ccomplex_t IT_0226 = m_s*s_24*m_C_1;
    const ccomplex_t IT_0227 = s_34*m_C_1*m_N_3;
    const ccomplex_t IT_0228 = 24*IT_0085;
    const ccomplex_t IT_0229 = (-6)*conj(IT_0070);
    const ccomplex_t IT_0230 = m_c*m_s*s_12;
    const ccomplex_t IT_0231 = -IT_0052;
    const ccomplex_t IT_0232 = IT_0051 + IT_0231;
    const ccomplex_t IT_0233 = IT_0053 + IT_0232;
    const ccomplex_t IT_0234 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0235 = m_c*m_s*m_C_1*m_N_3;
    const ccomplex_t IT_0236 = 6*IT_0070;
    const ccomplex_t IT_0237 = 24*IT_0175;
    const ccomplex_t IT_0238 = (-6)*IT_0106;
    const ccomplex_t IT_0239 = IT_0092*IT_0178 + conj(IT_0133)*IT_0180 +
       IT_0133*IT_0181 + IT_0071*(IT_0106 + IT_0192) + IT_0119*IT_0208 + IT_0177
      *conj(IT_0208) + IT_0218*IT_0223 + conj(IT_0092)*IT_0228 + IT_0050*IT_0234
       + (conj(IT_0106) + conj(IT_0192))*IT_0236 + conj(IT_0223)*IT_0237 + conj
      (IT_0050)*IT_0238;
    const ccomplex_t IT_0240 = IT_0054 + IT_0232;
    const ccomplex_t IT_0241 = 6*IT_0055;
    const ccomplex_t IT_0242 = 6*IT_0226;
    const ccomplex_t IT_0243 = 6*IT_0235;
    const ccomplex_t IT_0244 = 6*IT_0219;
    const ccomplex_t IT_0245 = 12*IT_0175;
    const ccomplex_t IT_0246 = 12*IT_0052;
    const ccomplex_t IT_0247 = 12*IT_0072;
    const ccomplex_t IT_0248 = 12*IT_0179;
    const ccomplex_t IT_0249 = (-6)*IT_0230;
    const ccomplex_t IT_0250 = (-6)*IT_0227;
    const ccomplex_t IT_0251 = IT_0179*IT_0228 + IT_0238*IT_0240 + IT_0192
      *IT_0241 + IT_0217*IT_0242 + IT_0134*IT_0243 + IT_0223*IT_0244 + IT_0226
      *IT_0245 + IT_0208*IT_0246 + IT_0092*IT_0247 + IT_0118*IT_0248 + IT_0165
      *IT_0249 + IT_0199*IT_0250;
    const ccomplex_t IT_0252 = 6*IT_0072;
    const ccomplex_t IT_0253 = 6*IT_0179;
    const ccomplex_t IT_0254 = 12*IT_0085;
    const ccomplex_t IT_0255 = 12*IT_0051;
    const ccomplex_t IT_0256 = 12*IT_0235;
    const ccomplex_t IT_0257 = (-12)*IT_0175;
    const ccomplex_t IT_0258 = (-12)*IT_0230;
    const ccomplex_t IT_0259 = (-12)*IT_0227;
    const ccomplex_t IT_0260 = (-6)*IT_0226;
    const ccomplex_t IT_0261 = (-6)*IT_0235;
    const ccomplex_t IT_0262 = (-6)*IT_0219;
    const ccomplex_t IT_0263 = IT_0233*IT_0238 + IT_0208*IT_0241 + IT_0092
      *IT_0252 + IT_0118*IT_0253 + IT_0179*IT_0254 + IT_0192*IT_0255 + IT_0050
      *IT_0256 + IT_0226*IT_0257 + IT_0165*IT_0258 + IT_0199*IT_0259 + IT_0217
      *IT_0260 + IT_0133*IT_0261 + IT_0223*IT_0262;
    const ccomplex_t IT_0264 = 6*IT_0106;
    const ccomplex_t IT_0265 = 6*IT_0230;
    const ccomplex_t IT_0266 = 6*IT_0227;
    const ccomplex_t IT_0267 = IT_0072*IT_0177 + IT_0052*IT_0228 + IT_0227
      *IT_0245 + IT_0118*IT_0246 + IT_0208*IT_0248 + IT_0134*IT_0252 + IT_0192
      *IT_0253 + IT_0092*IT_0256 + IT_0199*IT_0260 + IT_0165*IT_0262 + IT_0179
      *IT_0264 + IT_0223*IT_0265 + IT_0217*IT_0266;
    const ccomplex_t IT_0268 = 12*IT_0053;
    const ccomplex_t IT_0269 = 12*IT_0219;
    const ccomplex_t IT_0270 = IT_0180*IT_0226 + IT_0219*IT_0236 + IT_0053
      *IT_0237 + IT_0208*IT_0242 + IT_0165*IT_0252 + IT_0199*IT_0253 + IT_0227
      *IT_0254 + IT_0223*IT_0256 + IT_0192*IT_0260 + IT_0050*IT_0262 + IT_0092
      *IT_0265 + IT_0118*IT_0266 + IT_0217*IT_0268 + IT_0133*IT_0269;
    const ccomplex_t IT_0271 = IT_0177*IT_0179 + IT_0230*IT_0245 + IT_0092
      *IT_0246 + IT_0208*IT_0247 + IT_0192*IT_0252 + IT_0134*IT_0253 + IT_0118
      *IT_0256 + IT_0165*IT_0260 + IT_0199*IT_0262 + IT_0072*IT_0264 + IT_0217
      *IT_0265 + IT_0223*IT_0266;
    const ccomplex_t IT_0272 = 12*IT_0226;
    const ccomplex_t IT_0273 = IT_0180*IT_0219 + IT_0226*IT_0236 + IT_0208
      *IT_0244 + IT_0199*IT_0252 + IT_0165*IT_0253 + IT_0230*IT_0254 + IT_0217
      *IT_0256 + IT_0050*IT_0260 + IT_0192*IT_0262 + IT_0118*IT_0265 + IT_0092
      *IT_0266 + IT_0223*IT_0268 + IT_0133*IT_0272;
    const ccomplex_t IT_0274 = (-6)*IT_0070;
    const ccomplex_t IT_0275 = (-12)*IT_0085;
    const ccomplex_t IT_0276 = IT_0179*IT_0245 + IT_0208*IT_0250 + IT_0223
      *IT_0252 + IT_0217*IT_0253 + IT_0199*IT_0255 + IT_0165*IT_0256 + IT_0050
      *IT_0258 + IT_0192*IT_0259 + IT_0118*IT_0260 + IT_0092*IT_0262 + IT_0227
      *IT_0264 + IT_0133*IT_0265 + IT_0230*IT_0274 + IT_0226*IT_0275;
    const ccomplex_t IT_0277 = IT_0072*IT_0245 + IT_0208*IT_0249 + IT_0217
      *IT_0252 + IT_0223*IT_0253 + IT_0165*IT_0255 + IT_0199*IT_0256 + IT_0192
      *IT_0258 + IT_0050*IT_0259 + IT_0092*IT_0260 + IT_0118*IT_0262 + IT_0230
      *IT_0264 + IT_0133*IT_0266 + IT_0227*IT_0274 + IT_0219*IT_0275;
    const ccomplex_t IT_0278 = IT_0072*IT_0254;
    const ccomplex_t IT_0279 = (-6)*IT_0233;
    const ccomplex_t IT_0280 = IT_0118*IT_0252;
    const ccomplex_t IT_0281 = IT_0092*IT_0253;
    const ccomplex_t IT_0282 = IT_0208*IT_0243;
    const ccomplex_t IT_0283 = IT_0219*IT_0237 + IT_0192*IT_0261 + IT_0199
      *IT_0265 + IT_0165*IT_0266 + IT_0133*IT_0268 + IT_0217*IT_0269 + IT_0223
      *IT_0272 + IT_0240*IT_0274 + IT_0278 + IT_0050*IT_0279 + IT_0280 + IT_0281
       + IT_0282;
    const ccomplex_t IT_0284 = IT_0055*IT_0236 + IT_0050*IT_0255 + IT_0192
      *IT_0256 + IT_0219*IT_0257 + IT_0199*IT_0258 + IT_0165*IT_0259 + IT_0223
      *IT_0260 + IT_0217*IT_0262 + IT_0278 + IT_0133*IT_0279 + IT_0280 + IT_0281
       + IT_0282;
    const ccomplex_t IT_0285 = IT_0050*IT_0055*IT_0071 + IT_0072*(24*conj
      (IT_0070)*IT_0085 + 24*IT_0070*conj(IT_0085) + IT_0092*IT_0107 + IT_0118
      *IT_0119 + IT_0134*IT_0135 + IT_0165*IT_0176) + IT_0052*(48*IT_0085*conj
      (IT_0085) + conj(IT_0070)*IT_0177 + IT_0118*IT_0178) + IT_0179*(IT_0107
      *IT_0118 + IT_0092*IT_0119 + conj(IT_0085)*IT_0180 + IT_0085*IT_0181 +
       IT_0135*IT_0192 + IT_0176*IT_0199 + IT_0178*IT_0208) + IT_0053*(48
      *IT_0175*conj(IT_0175) + conj(IT_0106)*IT_0180 + IT_0217*IT_0218) +
       IT_0219*(IT_0119*IT_0175 + conj(IT_0175)*IT_0177 + IT_0071*IT_0217 +
       IT_0133*IT_0218 + IT_0181*IT_0223 + IT_0165*IT_0224 + IT_0050*IT_0225) + 
      (24*conj(IT_0106)*IT_0175 + 24*IT_0106*conj(IT_0175) + IT_0176*IT_0208 +
       IT_0181*IT_0217 + IT_0071*IT_0223 + IT_0199*IT_0224 + IT_0192*IT_0225)
      *IT_0226 + IT_0227*(IT_0118*IT_0176 + IT_0175*IT_0178 + IT_0107*IT_0199 +
       IT_0135*IT_0217 + conj(IT_0175)*IT_0228 + IT_0165*IT_0229) + (IT_0107
      *IT_0165 + IT_0092*IT_0176 + IT_0135*IT_0223 + IT_0199*IT_0229)*IT_0230 +
       IT_0192*IT_0233*IT_0234 + IT_0235*IT_0239 + (IT_0133*IT_0229 + IT_0208
      *IT_0234)*IT_0240 + conj(IT_0208)*IT_0251 + conj(IT_0192)*IT_0263 + conj
      (IT_0118)*IT_0267 + conj(IT_0217)*IT_0270 + conj(IT_0092)*IT_0271 + conj
      (IT_0223)*IT_0273 + conj(IT_0199)*IT_0276 + conj(IT_0165)*IT_0277 + conj
      (IT_0133)*IT_0283 + conj(IT_0050)*IT_0284;
    return create_ccomplex_return(IT_0285);
}

