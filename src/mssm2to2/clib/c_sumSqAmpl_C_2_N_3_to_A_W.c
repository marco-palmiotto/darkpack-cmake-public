#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_A_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + 0.5*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = pow(m_N_3, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0018 + (-0.5)*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0013*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = s_13*s_23;
    const ccomplex_t IT_0028 = pow(m_W, -2);
    const ccomplex_t IT_0029 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = s_12 + IT_0030;
    const ccomplex_t IT_0032 = m_N_3*IT_0031;
    const ccomplex_t IT_0033 = (-2)*IT_0000;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = m_C_2*IT_0000;
    const ccomplex_t IT_0037 = IT_0015*IT_0036;
    const ccomplex_t IT_0038 = 2*IT_0037;
    const ccomplex_t IT_0039 = pow(m_W, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0039 
      + -reg_prop, -1);
    const ccomplex_t IT_0041 = sin(beta);
    const ccomplex_t IT_0042 = IT_0003*IT_0041;
    const ccomplex_t IT_0043 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = 2*IT_0044;
    const ccomplex_t IT_0046 = cos(theta_W);
    const ccomplex_t IT_0047 = cpow(IT_0046, -1);
    const ccomplex_t IT_0048 = IT_0041*IT_0047;
    const ccomplex_t IT_0049 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0053 = IT_0042*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = IT_0051 + IT_0054;
    const ccomplex_t IT_0056 = IT_0045 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = IT_0040*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = m_N_3*IT_0000;
    const ccomplex_t IT_0064 = IT_0025*IT_0063;
    const ccomplex_t IT_0065 = -IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = (-2)*IT_0065;
    const ccomplex_t IT_0067 = IT_0038 + IT_0066;
    const ccomplex_t IT_0068 = 2*conj(IT_0037);
    const ccomplex_t IT_0069 = (-2)*conj(IT_0065);
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = IT_0015*IT_0033;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0015*IT_0063;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = 2*IT_0074;
    const ccomplex_t IT_0076 = cos(beta);
    const ccomplex_t IT_0077 = IT_0003*IT_0076;
    const ccomplex_t IT_0078 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (-2)*IT_0079;
    const ccomplex_t IT_0081 = IT_0047*IT_0076;
    const ccomplex_t IT_0082 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0077*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0080 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0059*IT_0091;
    const ccomplex_t IT_0093 = IT_0040*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0025*IT_0036;
    const ccomplex_t IT_0096 = IT_0094 + -IT_0095;
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = IT_0075 + IT_0097;
    const ccomplex_t IT_0099 = 2*conj(IT_0074);
    const ccomplex_t IT_0100 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0101 = IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = pow(s_24, 2);
    const ccomplex_t IT_0103 = IT_0028*IT_0102;
    const ccomplex_t IT_0104 = -IT_0012;
    const ccomplex_t IT_0105 = IT_0103 + IT_0104;
    const ccomplex_t IT_0106 = m_C_2*IT_0105;
    const ccomplex_t IT_0107 = pow(s_14, 2);
    const ccomplex_t IT_0108 = IT_0028*IT_0107;
    const ccomplex_t IT_0109 = -IT_0011;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = m_N_3*IT_0110;
    const ccomplex_t IT_0112 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0000*IT_0023;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = -IT_0034 + IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0065;
    const ccomplex_t IT_0118 = (-2)*IT_0037;
    const ccomplex_t IT_0119 = IT_0117 + IT_0118;
    const ccomplex_t IT_0120 = 2*conj(IT_0065);
    const ccomplex_t IT_0121 = (-2)*conj(IT_0037);
    const ccomplex_t IT_0122 = IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = IT_0000*IT_0010;
    const ccomplex_t IT_0124 = IT_0112*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = -IT_0071 + IT_0125;
    const ccomplex_t IT_0127 = 2*IT_0096;
    const ccomplex_t IT_0128 = (-2)*IT_0074;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = 2*conj(IT_0096);
    const ccomplex_t IT_0131 = (-2)*conj(IT_0074);
    const ccomplex_t IT_0132 = IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = m_C_2*IT_0031;
    const ccomplex_t IT_0134 = s_13*m_N_3;
    const ccomplex_t IT_0135 = -IT_0096;
    const ccomplex_t IT_0136 = IT_0074 + IT_0135;
    const ccomplex_t IT_0137 = -conj(IT_0096);
    const ccomplex_t IT_0138 = conj(IT_0074) + IT_0137;
    const ccomplex_t IT_0139 = 6*IT_0037;
    const ccomplex_t IT_0140 = (-6)*IT_0065;
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = 6*conj(IT_0037);
    const ccomplex_t IT_0143 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = s_23*m_C_2;
    const ccomplex_t IT_0146 = 6*IT_0065;
    const ccomplex_t IT_0147 = (-6)*IT_0037;
    const ccomplex_t IT_0148 = IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = 6*conj(IT_0065);
    const ccomplex_t IT_0150 = (-6)*conj(IT_0037);
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = s_12*IT_0110;
    const ccomplex_t IT_0153 = s_12*IT_0105;
    const ccomplex_t IT_0154 = IT_0000*IT_0015;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = 4*IT_0155*conj(IT_0155);
    const ccomplex_t IT_0157 = IT_0000*IT_0025;
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = 4*IT_0158*conj(IT_0158);
    const ccomplex_t IT_0160 = 6*s_12;
    const ccomplex_t IT_0161 = pow(s_34, 2);
    const ccomplex_t IT_0162 = s_12*IT_0028*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = s_13*s_24*s_34*IT_0028;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = s_14*s_23*s_34*IT_0028;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = IT_0027 + IT_0163 + IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = -IT_0000;
    const ccomplex_t IT_0170 = IT_0015*IT_0169;
    const ccomplex_t IT_0171 = 8*IT_0170;
    const ccomplex_t IT_0172 = 8*conj(IT_0170);
    const ccomplex_t IT_0173 = IT_0025*IT_0169;
    const ccomplex_t IT_0174 = 8*IT_0173;
    const ccomplex_t IT_0175 = 8*conj(IT_0173);
    const ccomplex_t IT_0176 = m_C_2*m_N_3*IT_0028*IT_0161;
    const ccomplex_t IT_0177 = (-0.5)*IT_0125;
    const ccomplex_t IT_0178 = 4*IT_0170;
    const ccomplex_t IT_0179 = 4*IT_0016;
    const ccomplex_t IT_0180 = (-0.5)*IT_0115;
    const ccomplex_t IT_0181 = 4*conj(IT_0016);
    const ccomplex_t IT_0182 = 4*conj(IT_0170)*(IT_0026 + 3*IT_0173) + IT_0175
      *IT_0177 + 4*(IT_0026 + 0.25*IT_0174)*conj(IT_0177) + 4*conj(IT_0026)*
      (IT_0177 + 0.25*IT_0178) + conj(IT_0173)*(12*IT_0170 + IT_0179) + IT_0172
      *IT_0180 + (IT_0171 + IT_0179)*conj(IT_0180) + (IT_0173 + IT_0180)*IT_0181;
    const ccomplex_t IT_0183 = m_C_2*m_N_3;
    const ccomplex_t IT_0184 = IT_0110*IT_0183;
    const ccomplex_t IT_0185 = 2*IT_0155;
    const ccomplex_t IT_0186 = 2*conj(IT_0155);
    const ccomplex_t IT_0187 = 8*conj(IT_0155)*IT_0158;
    const ccomplex_t IT_0188 = IT_0105*IT_0183;
    const ccomplex_t IT_0189 = (-6)*s_12;
    const ccomplex_t IT_0190 = (-6)*s_13;
    const ccomplex_t IT_0191 = (-6)*s_23;
    const ccomplex_t IT_0192 = (-6)*IT_0011;
    const ccomplex_t IT_0193 = (-6)*IT_0012;
    const ccomplex_t IT_0194 = conj(IT_0096)*IT_0141 + IT_0096*IT_0144 + conj
      (IT_0074)*IT_0148 + IT_0074*IT_0151 + (conj(IT_0072)*IT_0116 + IT_0072
      *conj(IT_0116) + conj(IT_0035)*IT_0126 + IT_0035*conj(IT_0126))*IT_0189 + 
      (conj(IT_0116)*IT_0170 + IT_0116*conj(IT_0170) + conj(IT_0126)*IT_0173 +
       IT_0126*conj(IT_0173))*IT_0190 + (conj(IT_0035)*IT_0170 + IT_0035*conj
      (IT_0170) + conj(IT_0072)*IT_0173 + IT_0072*conj(IT_0173))*IT_0191 + (conj
      (IT_0116)*IT_0126 + IT_0116*conj(IT_0126))*IT_0192 + (conj(IT_0035)
      *IT_0072 + IT_0035*conj(IT_0072))*IT_0193;
    const ccomplex_t IT_0195 = s_12*IT_0039;
    const ccomplex_t IT_0196 = s_14*s_24;
    const ccomplex_t IT_0197 = (-2)*IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = IT_0028*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = s_12 + IT_0200;
    const ccomplex_t IT_0202 = (-4)*s_12;
    const ccomplex_t IT_0203 = (-4)*s_13;
    const ccomplex_t IT_0204 = (-4)*s_23;
    const ccomplex_t IT_0205 = (-4)*IT_0011;
    const ccomplex_t IT_0206 = (-4)*IT_0012;
    const ccomplex_t IT_0207 = (conj(IT_0035)*IT_0116 + IT_0035*conj(IT_0116) 
      + conj(IT_0072)*IT_0126 + IT_0072*conj(IT_0126))*IT_0202 + (conj(IT_0126)
      *IT_0170 + IT_0126*conj(IT_0170) + conj(IT_0116)*IT_0173 + IT_0116*conj
      (IT_0173))*IT_0203 + (conj(IT_0072)*IT_0170 + IT_0072*conj(IT_0170) + conj
      (IT_0035)*IT_0173 + IT_0035*conj(IT_0173))*IT_0204 + (IT_0116*conj(IT_0116
      ) + IT_0126*conj(IT_0126))*IT_0205 + (IT_0035*conj(IT_0035) + IT_0072*conj
      (IT_0072))*IT_0206;
    const ccomplex_t IT_0208 = s_14*s_34*IT_0028;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = s_13 + IT_0209;
    const ccomplex_t IT_0211 = s_12*IT_0210;
    const ccomplex_t IT_0212 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = s_23 + IT_0213;
    const ccomplex_t IT_0215 = s_12*IT_0214;
    const ccomplex_t IT_0216 = (-4)*IT_0170;
    const ccomplex_t IT_0217 = conj(IT_0155)*IT_0216;
    const ccomplex_t IT_0218 = (-4)*conj(IT_0170);
    const ccomplex_t IT_0219 = IT_0155*IT_0218;
    const ccomplex_t IT_0220 = (-4)*IT_0173;
    const ccomplex_t IT_0221 = conj(IT_0158)*IT_0220;
    const ccomplex_t IT_0222 = (-4)*conj(IT_0173);
    const ccomplex_t IT_0223 = IT_0158*IT_0222;
    const ccomplex_t IT_0224 = s_23*IT_0028*IT_0107;
    const ccomplex_t IT_0225 = s_23*IT_0011;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = s_12*s_14*s_34*IT_0028;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = s_12*s_13;
    const ccomplex_t IT_0230 = 2*IT_0229;
    const ccomplex_t IT_0231 = s_13*s_14*s_24*IT_0028;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = IT_0224 + IT_0226 + IT_0228 + IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = s_14*s_23*s_24*IT_0028;
    const ccomplex_t IT_0235 = s_13*IT_0012;
    const ccomplex_t IT_0236 = s_12*s_24*s_34*IT_0028;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = s_13*IT_0028*IT_0102;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = IT_0234 + IT_0235 + IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = (-4)*IT_0155*conj(IT_0177);
    const ccomplex_t IT_0242 = (-4)*conj(IT_0155)*IT_0177;
    const ccomplex_t IT_0243 = (-4)*IT_0158*conj(IT_0180);
    const ccomplex_t IT_0244 = (-4)*conj(IT_0158)*IT_0180;
    const ccomplex_t IT_0245 = IT_0011*IT_0031;
    const ccomplex_t IT_0246 = 2*IT_0126;
    const ccomplex_t IT_0247 = IT_0012*IT_0031;
    const ccomplex_t IT_0248 = 2*IT_0072;
    const ccomplex_t IT_0249 = s_12*IT_0031;
    const ccomplex_t IT_0250 = (-8)*IT_0155;
    const ccomplex_t IT_0251 = IT_0011*IT_0105;
    const ccomplex_t IT_0252 = (-2)*IT_0126;
    const ccomplex_t IT_0253 = IT_0012*IT_0110;
    const ccomplex_t IT_0254 = (-2)*IT_0072;
    const ccomplex_t IT_0255 = 2*conj(IT_0126);
    const ccomplex_t IT_0256 = 2*conj(IT_0072);
    const ccomplex_t IT_0257 = (-2)*conj(IT_0126);
    const ccomplex_t IT_0258 = (-2)*conj(IT_0072);
    const ccomplex_t IT_0259 = IT_0031*IT_0183;
    const ccomplex_t IT_0260 = IT_0183*IT_0210;
    const ccomplex_t IT_0261 = IT_0183*IT_0214;
    const ccomplex_t IT_0262 = 2*IT_0116;
    const ccomplex_t IT_0263 = 2*IT_0035;
    const ccomplex_t IT_0264 = (-8)*IT_0158;
    const ccomplex_t IT_0265 = (-2)*IT_0116;
    const ccomplex_t IT_0266 = (-2)*IT_0035;
    const ccomplex_t IT_0267 = 0.125*IT_0032*IT_0067 + 0.125*IT_0098*IT_0106 +
       0.125*IT_0111*IT_0119 + 0.125*IT_0129*IT_0133 + 0.125*IT_0184*IT_0246 +
       0.125*IT_0188*IT_0248 + 0.125*(IT_0252 + IT_0254)*IT_0259 + IT_0155*
      (IT_0184 + IT_0188 + (-2)*IT_0259) + (-1.25)*IT_0170*(IT_0260 + IT_0261) +
       0.125*IT_0245*IT_0262 + 0.125*IT_0247*IT_0263 + 0.125*IT_0249*IT_0264 +
       0.125*IT_0251*IT_0265 + 0.125*IT_0253*IT_0266;
    const ccomplex_t IT_0268 = 2*conj(IT_0116);
    const ccomplex_t IT_0269 = 2*conj(IT_0035);
    const ccomplex_t IT_0270 = (-2)*conj(IT_0116);
    const ccomplex_t IT_0271 = (-2)*conj(IT_0035);
    const ccomplex_t IT_0272 = IT_0032*IT_0070 + IT_0101*IT_0106 + IT_0111
      *IT_0122 + IT_0132*IT_0133 + IT_0184*IT_0255 + IT_0188*IT_0256 + ((-16)
      *conj(IT_0155) + IT_0257 + IT_0258)*IT_0259 + IT_0245*IT_0268 + IT_0247
      *IT_0269 + IT_0251*IT_0270 + IT_0253*IT_0271;
    const ccomplex_t IT_0273 = (-2)*IT_0155;
    const ccomplex_t IT_0274 = (-2)*conj(IT_0155);
    const ccomplex_t IT_0275 = IT_0224 + IT_0226 + IT_0228 + IT_0231;
    const ccomplex_t IT_0276 = (-2)*IT_0170;
    const ccomplex_t IT_0277 = conj(IT_0155)*IT_0276;
    const ccomplex_t IT_0278 = (-2)*conj(IT_0170);
    const ccomplex_t IT_0279 = IT_0155*IT_0278;
    const ccomplex_t IT_0280 = conj(IT_0158)*IT_0173;
    const ccomplex_t IT_0281 = (-2)*IT_0280;
    const ccomplex_t IT_0282 = IT_0158*conj(IT_0173);
    const ccomplex_t IT_0283 = (-2)*IT_0282;
    const ccomplex_t IT_0284 = conj(IT_0016)*IT_0273;
    const ccomplex_t IT_0285 = IT_0016*IT_0274;
    const ccomplex_t IT_0286 = conj(IT_0026)*IT_0158;
    const ccomplex_t IT_0287 = (-2)*IT_0286;
    const ccomplex_t IT_0288 = IT_0026*conj(IT_0158);
    const ccomplex_t IT_0289 = (-2)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0235;
    const ccomplex_t IT_0291 = IT_0234 + IT_0237 + IT_0238 + IT_0290;
    const ccomplex_t IT_0292 = m_N_3*IT_0210;
    const ccomplex_t IT_0293 = IT_0011*IT_0214;
    const ccomplex_t IT_0294 = IT_0012*IT_0210;
    const ccomplex_t IT_0295 = m_C_2*IT_0214;
    const ccomplex_t IT_0296 = -IT_0136;
    const ccomplex_t IT_0297 = (-0.25)*IT_0293;
    const ccomplex_t IT_0298 = (-0.25)*IT_0294;
    const ccomplex_t IT_0299 = s_23*IT_0011*IT_0039;
    const ccomplex_t IT_0300 = s_24*s_34*IT_0011;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0299 + IT_0301;
    const ccomplex_t IT_0303 = IT_0028*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = IT_0225 + IT_0304;
    const ccomplex_t IT_0306 = s_12*s_23*IT_0039;
    const ccomplex_t IT_0307 = s_13*IT_0102;
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = s_14*s_23*s_24;
    const ccomplex_t IT_0310 = s_12*s_24*s_34;
    const ccomplex_t IT_0311 = s_13*IT_0012*IT_0039;
    const ccomplex_t IT_0312 = -IT_0309 + -IT_0310 + (-0.5)*IT_0311;
    const ccomplex_t IT_0313 = IT_0308 + IT_0312;
    const ccomplex_t IT_0314 = IT_0028*IT_0313;
    const ccomplex_t IT_0315 = s_12*s_23;
    const ccomplex_t IT_0316 = (-2)*IT_0315;
    const ccomplex_t IT_0317 = IT_0235 + IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = s_13*s_23*IT_0039;
    const ccomplex_t IT_0319 = s_12*IT_0161;
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = s_14*s_23*s_34;
    const ccomplex_t IT_0322 = s_13*s_24*s_34;
    const ccomplex_t IT_0323 = -IT_0321 + -IT_0322;
    const ccomplex_t IT_0324 = IT_0320 + IT_0323;
    const ccomplex_t IT_0325 = IT_0028*IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = IT_0027 + IT_0326;
    const ccomplex_t IT_0328 = s_23*IT_0039;
    const ccomplex_t IT_0329 = s_24*s_34;
    const ccomplex_t IT_0330 = (-2)*IT_0329;
    const ccomplex_t IT_0331 = IT_0328 + IT_0330;
    const ccomplex_t IT_0332 = IT_0028*IT_0331;
    const ccomplex_t IT_0333 = (-0.5)*IT_0332;
    const ccomplex_t IT_0334 = s_23 + IT_0333;
    const ccomplex_t IT_0335 = m_C_2*IT_0334;
    const ccomplex_t IT_0336 = (-2)*IT_0335;
    const ccomplex_t IT_0337 = (-0.25)*IT_0336;
    const ccomplex_t IT_0338 = s_13*IT_0039;
    const ccomplex_t IT_0339 = s_14*s_34;
    const ccomplex_t IT_0340 = (-2)*IT_0339;
    const ccomplex_t IT_0341 = IT_0338 + IT_0340;
    const ccomplex_t IT_0342 = IT_0028*IT_0341;
    const ccomplex_t IT_0343 = (-0.5)*IT_0342;
    const ccomplex_t IT_0344 = s_13 + IT_0343;
    const ccomplex_t IT_0345 = m_N_3*IT_0344;
    const ccomplex_t IT_0346 = (-2)*IT_0345;
    const ccomplex_t IT_0347 = (-0.25)*IT_0346;
    const ccomplex_t IT_0348 = IT_0028*IT_0338;
    const ccomplex_t IT_0349 = (-0.25)*IT_0348;
    const ccomplex_t IT_0350 = s_13 + IT_0349;
    const ccomplex_t IT_0351 = IT_0183*IT_0350;
    const ccomplex_t IT_0352 = (-4)*IT_0351;
    const ccomplex_t IT_0353 = (-0.25)*IT_0352;
    const ccomplex_t IT_0354 = IT_0028*IT_0328;
    const ccomplex_t IT_0355 = (-0.25)*IT_0354;
    const ccomplex_t IT_0356 = s_23 + IT_0355;
    const ccomplex_t IT_0357 = IT_0183*IT_0356;
    const ccomplex_t IT_0358 = (-4)*IT_0357;
    const ccomplex_t IT_0359 = (-0.25)*IT_0358;
    const ccomplex_t IT_0360 = s_34*IT_0028;
    const ccomplex_t IT_0361 = s_14*s_23;
    const ccomplex_t IT_0362 = s_13*s_24;
    const ccomplex_t IT_0363 = s_12*s_34;
    const ccomplex_t IT_0364 = -IT_0362 + -IT_0363;
    const ccomplex_t IT_0365 = IT_0361 + IT_0364;
    const ccomplex_t IT_0366 = IT_0360*IT_0365;
    const ccomplex_t IT_0367 = (-0.25)*IT_0366;
    const ccomplex_t IT_0368 = (-2)*IT_0322;
    const ccomplex_t IT_0369 = IT_0318 + IT_0368;
    const ccomplex_t IT_0370 = IT_0028*IT_0369;
    const ccomplex_t IT_0371 = (-8)*IT_0370;
    const ccomplex_t IT_0372 = 16*IT_0027;
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = conj(IT_0126)*IT_0305 + -conj(IT_0072)*IT_0317 
      + 2*conj(IT_0016)*IT_0327 + IT_0122*IT_0337 + IT_0101*IT_0347 + IT_0268
      *IT_0353 + IT_0269*IT_0359 + IT_0218*IT_0367 + (-0.25)*conj(IT_0177)
      *IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0116)*IT_0305 + -conj(IT_0035)*IT_0317 
      + 2*conj(IT_0026)*IT_0327 + IT_0132*IT_0337 + IT_0070*IT_0347 + IT_0255
      *IT_0353 + IT_0256*IT_0359 + IT_0222*IT_0367 + (-0.25)*conj(IT_0180)
      *IT_0373;
    const ccomplex_t IT_0376 = (-2)*IT_0016;
    const ccomplex_t IT_0377 = IT_0276 + IT_0376;
    const ccomplex_t IT_0378 = (-2)*conj(IT_0016);
    const ccomplex_t IT_0379 = IT_0278 + IT_0378;
    const ccomplex_t IT_0380 = (-10)*IT_0158*conj(IT_0170);
    const ccomplex_t IT_0381 = (-10)*IT_0155*conj(IT_0173);
    const ccomplex_t IT_0382 = (-10)*conj(IT_0155)*IT_0173;
    const ccomplex_t IT_0383 = (-8)*conj(IT_0155)*IT_0180;
    const ccomplex_t IT_0384 = conj(IT_0177)*IT_0264;
    const ccomplex_t IT_0385 = (-8)*conj(IT_0158)*IT_0177;
    const ccomplex_t IT_0386 = conj(IT_0026) + conj(IT_0173);
    const ccomplex_t IT_0387 = conj(IT_0026)*IT_0273;
    const ccomplex_t IT_0388 = IT_0026*IT_0274;
    const ccomplex_t IT_0389 = conj(IT_0158)*IT_0376;
    const ccomplex_t IT_0390 = IT_0158*IT_0378;
    const ccomplex_t IT_0391 = conj(IT_0180)*IT_0250;
    const ccomplex_t IT_0392 = IT_0260*((IT_0026 + IT_0173)*IT_0257 + conj
      (IT_0116)*IT_0377 + IT_0116*IT_0379 + IT_0380 + IT_0381 + IT_0382 +
       IT_0383 + IT_0384 + IT_0385 + IT_0252*IT_0386 + IT_0387 + IT_0388 +
       IT_0389 + IT_0390 + IT_0391);
    const ccomplex_t IT_0393 = IT_0261*((IT_0026 + IT_0173)*IT_0258 + conj
      (IT_0035)*IT_0377 + IT_0035*IT_0379 + IT_0380 + IT_0381 + IT_0382 +
       IT_0383 + IT_0384 + IT_0385 + IT_0254*IT_0386 + IT_0387 + IT_0388 +
       IT_0389 + IT_0390 + IT_0391);
    const ccomplex_t IT_0394 = (12*IT_0016*conj(IT_0016) + 12*IT_0026*conj
      (IT_0026))*IT_0027 + IT_0032*(conj(IT_0035)*IT_0067 + IT_0035*IT_0070 +
       conj(IT_0072)*IT_0098 + IT_0072*IT_0101) + (IT_0070*IT_0072 + IT_0067
      *conj(IT_0072) + conj(IT_0035)*IT_0098 + IT_0035*IT_0101)*IT_0106 +
       IT_0111*(conj(IT_0116)*IT_0119 + IT_0116*IT_0122 + conj(IT_0126)*IT_0129 
      + IT_0126*IT_0132) + (IT_0122*IT_0126 + IT_0119*conj(IT_0126) + conj
      (IT_0116)*IT_0129 + IT_0116*IT_0132)*IT_0133 + 6*IT_0134*(conj(IT_0016)
      *IT_0136 + IT_0016*IT_0138 + 0.166666666666667*conj(IT_0026)*IT_0141 +
       0.166666666666667*IT_0026*IT_0144) + (-6)*IT_0145*(conj(IT_0026)*IT_0136 
      + IT_0026*IT_0138 + (-0.166666666666667)*conj(IT_0016)*IT_0148 + (
      -0.166666666666667)*IT_0016*IT_0151) + (IT_0152 + IT_0153)*(IT_0156 +
       IT_0159) + (IT_0037*conj(IT_0037) + IT_0065*conj(IT_0065) + IT_0074*conj
      (IT_0074) + IT_0096*conj(IT_0096))*IT_0160 + IT_0168*(conj(IT_0016)
      *IT_0171 + IT_0016*IT_0172 + conj(IT_0026)*IT_0174 + IT_0026*IT_0175) +
       IT_0176*IT_0182 + IT_0184*(conj(IT_0116)*IT_0185 + IT_0116*IT_0186 +
       IT_0187) + (conj(IT_0035)*IT_0185 + IT_0035*IT_0186 + IT_0187)*IT_0188 + 
      (conj(IT_0037)*IT_0065 + IT_0037*conj(IT_0065) + conj(IT_0074)*IT_0096 +
       IT_0074*conj(IT_0096))*IT_0189 + IT_0183*IT_0194 + IT_0201*IT_0207 + 
      (IT_0211 + IT_0215)*(IT_0217 + IT_0219 + IT_0221 + IT_0223) + (IT_0233 +
       IT_0240)*(IT_0241 + IT_0242 + IT_0243 + IT_0244) + conj(IT_0155)*(IT_0032
      *IT_0098 + IT_0067*IT_0106 + IT_0111*IT_0129 + IT_0119*IT_0133 + IT_0245
      *IT_0246 + IT_0247*IT_0248 + IT_0249*IT_0250 + IT_0251*IT_0252 + IT_0253
      *IT_0254) + IT_0155*(IT_0032*IT_0101 + IT_0070*IT_0106 + IT_0111*IT_0132 +
       IT_0122*IT_0133 + IT_0245*IT_0255 + IT_0247*IT_0256 + IT_0251*IT_0257 +
       IT_0253*IT_0258) + 8*conj(IT_0158)*IT_0267 + IT_0158*IT_0272 + IT_0259*(
      (conj(IT_0035) + conj(IT_0116))*IT_0273 + (IT_0035 + IT_0116)*IT_0274) +
       IT_0275*(conj(IT_0016)*IT_0252 + IT_0016*IT_0257 + conj(IT_0026)*IT_0265 
      + IT_0026*IT_0270 + IT_0277 + IT_0279 + IT_0281 + IT_0283 + IT_0284 +
       IT_0285 + IT_0287 + IT_0289) + (conj(IT_0016)*IT_0254 + IT_0016*IT_0258 +
       conj(IT_0026)*IT_0266 + IT_0026*IT_0271 + IT_0277 + IT_0279 + IT_0281 +
       IT_0283 + IT_0284 + IT_0285 + IT_0287 + IT_0289)*IT_0291 + 4*conj(IT_0173
      )*(IT_0162*IT_0173 + 0.25*IT_0168*IT_0174 + (IT_0037 + -IT_0065)*IT_0292 +
       0.25*IT_0262*IT_0293 + 0.25*IT_0263*IT_0294 + IT_0295*IT_0296) + (-4)
      *conj(IT_0170)*((-0.25)*IT_0168*IT_0171 + (-0.25)*IT_0162*IT_0178 + 
      (IT_0037 + -IT_0065)*IT_0295 + IT_0292*IT_0296 + IT_0246*IT_0297 + IT_0248
      *IT_0298) + 4*IT_0170*(IT_0138*IT_0292 + -(conj(IT_0037) + -conj(IT_0065))
      *IT_0295 + -IT_0255*IT_0297 + -IT_0256*IT_0298) + 4*IT_0173*((conj(IT_0037
      ) + -conj(IT_0065))*IT_0292 + -IT_0138*IT_0295 + -IT_0268*IT_0297 + 
      -IT_0269*IT_0298) + (-4)*conj(IT_0177)*(IT_0126*IT_0305 + -IT_0072*IT_0317
       + 2*IT_0016*IT_0327 + IT_0119*IT_0337 + IT_0098*IT_0347 + IT_0262*IT_0353
       + IT_0263*IT_0359 + IT_0216*IT_0367) + (-4)*conj(IT_0180)*(IT_0116
      *IT_0305 + -IT_0035*IT_0317 + 2*IT_0026*IT_0327 + IT_0129*IT_0337 +
       IT_0067*IT_0347 + IT_0246*IT_0353 + IT_0248*IT_0359 + IT_0220*IT_0367) + 
      (-4)*IT_0177*IT_0374 + (-4)*IT_0180*IT_0375 + IT_0392 + IT_0393;
    return create_ccomplex_return(IT_0394);
}

