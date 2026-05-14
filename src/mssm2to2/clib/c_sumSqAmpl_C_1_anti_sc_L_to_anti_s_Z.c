#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sc_L_to_anti_s_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_sc_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = e_em*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0008;
    const ccomplex_t IT_0013 = e_em*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0011 + 3*IT_0013);
    const ccomplex_t IT_0015 = (-0.166666666666667)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = m_s*m_C_1;
    const ccomplex_t IT_0020 = pow(m_Z, 2);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (-3)*IT_0021;
    const ccomplex_t IT_0023 = pow(m_W, -1);
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0005*IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0009;
    const ccomplex_t IT_0029 = 0.333333333333333*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = IT_0012*IT_0034;
    const ccomplex_t IT_0037 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0038 = IT_0012*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0036 + (-2)
      *IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = IT_0012*IT_0048;
    const ccomplex_t IT_0051 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0052 = IT_0012*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + -IT_0050 + (-2)
      *IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = IT_0045*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (-2)*IT_0044 + (-2)*IT_0057;
    const ccomplex_t IT_0059 = s_13*s_14;
    const ccomplex_t IT_0060 = s_34*IT_0001;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0059 + IT_0061;
    const ccomplex_t IT_0063 = pow(m_Z, -2);
    const ccomplex_t IT_0064 = pow(m_Z, 4);
    const ccomplex_t IT_0065 = s_13*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = 6*IT_0066;
    const ccomplex_t IT_0068 = s_14*s_34;
    const ccomplex_t IT_0069 = 24*IT_0068;
    const ccomplex_t IT_0070 = s_13*IT_0020;
    const ccomplex_t IT_0071 = (-12)*IT_0070;
    const ccomplex_t IT_0072 = IT_0067 + IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = IT_0044 + IT_0057;
    const ccomplex_t IT_0074 = (-6)*IT_0066;
    const ccomplex_t IT_0075 = IT_0069 + IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*-3)*e_em*(IT_0005*IT_0008 + -1.
      /3*IT_0004*IT_0009);
    const ccomplex_t IT_0077 = 0.333333333333333*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0007*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0077*IT_0080;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = 2*IT_0018 + IT_0081;
    const ccomplex_t IT_0084 = s_14*IT_0000;
    const ccomplex_t IT_0085 = s_13*s_34;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0084 + IT_0086;
    const ccomplex_t IT_0088 = (-6)*IT_0018;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0005*IT_0023*IT_0025;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0054*IT_0091;
    const ccomplex_t IT_0093 = IT_0045*IT_0089*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0095 = IT_0027*IT_0040;
    const ccomplex_t IT_0096 = IT_0042*IT_0094*IT_0095;
    const ccomplex_t IT_0097 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0098 = IT_0010*IT_0097;
    const ccomplex_t IT_0099 = IT_0012*IT_0097;
    const ccomplex_t IT_0100 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0101 = IT_0012*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0099 + (-2)
      *IT_0101);
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0027*IT_0103;
    const ccomplex_t IT_0105 = IT_0042*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = m_C_1*IT_0106;
    const ccomplex_t IT_0108 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0109 = IT_0010*IT_0108;
    const ccomplex_t IT_0110 = IT_0012*IT_0108;
    const ccomplex_t IT_0111 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0112 = IT_0012*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + -IT_0110 + (-2)
      *IT_0112);
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0091*IT_0114;
    const ccomplex_t IT_0116 = IT_0045*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = m_C_1*IT_0117;
    const ccomplex_t IT_0119 = IT_0093 + IT_0096 + IT_0107 + IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0121 = IT_0007*IT_0029;
    const ccomplex_t IT_0122 = IT_0003*IT_0120*IT_0121;
    const ccomplex_t IT_0123 = m_s*IT_0018;
    const ccomplex_t IT_0124 = -IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = conj(IT_0119) + conj(IT_0124);
    const ccomplex_t IT_0126 = s_14*IT_0020;
    const ccomplex_t IT_0127 = IT_0063*IT_0126;
    const ccomplex_t IT_0128 = (-0.25)*IT_0127;
    const ccomplex_t IT_0129 = s_14 + IT_0128;
    const ccomplex_t IT_0130 = m_s*IT_0129;
    const ccomplex_t IT_0131 = (-4)*IT_0130;
    const ccomplex_t IT_0132 = IT_0088*IT_0131;
    const ccomplex_t IT_0133 = (-2)*IT_0106 + (-2)*IT_0117;
    const ccomplex_t IT_0134 = 6*IT_0062;
    const ccomplex_t IT_0135 = IT_0133*IT_0134;
    const ccomplex_t IT_0136 = IT_0027*IT_0078;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0077*IT_0137;
    const ccomplex_t IT_0139 = (-6)*IT_0062;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = (-2)*IT_0032 + -IT_0138;
    const ccomplex_t IT_0142 = (-6)*IT_0087;
    const ccomplex_t IT_0143 = IT_0141*IT_0142;
    const ccomplex_t IT_0144 = IT_0106 + IT_0117;
    const ccomplex_t IT_0145 = 3*IT_0021;
    const ccomplex_t IT_0146 = conj(IT_0133)*IT_0134;
    const ccomplex_t IT_0147 = 6*IT_0022;
    const ccomplex_t IT_0148 = IT_0047*IT_0114;
    const ccomplex_t IT_0149 = IT_0045*IT_0089*IT_0148;
    const ccomplex_t IT_0150 = m_C_1*IT_0057;
    const ccomplex_t IT_0151 = IT_0007*IT_0103;
    const ccomplex_t IT_0152 = IT_0042*IT_0094*IT_0151;
    const ccomplex_t IT_0153 = m_C_1*IT_0044;
    const ccomplex_t IT_0154 = IT_0149 + IT_0150 + IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = IT_0015*IT_0027;
    const ccomplex_t IT_0156 = IT_0003*IT_0120*IT_0155;
    const ccomplex_t IT_0157 = m_s*IT_0032;
    const ccomplex_t IT_0158 = IT_0156 + -IT_0157;
    const ccomplex_t IT_0159 = conj(IT_0154) + conj(IT_0158);
    const ccomplex_t IT_0160 = 0.5*IT_0127;
    const ccomplex_t IT_0161 = s_14 + IT_0160;
    const ccomplex_t IT_0162 = m_s*IT_0161;
    const ccomplex_t IT_0163 = 2*IT_0162;
    const ccomplex_t IT_0164 = 6*IT_0163;
    const ccomplex_t IT_0165 = s_34*IT_0020;
    const ccomplex_t IT_0166 = IT_0063*IT_0165;
    const ccomplex_t IT_0167 = (-0.25)*IT_0166;
    const ccomplex_t IT_0168 = s_34 + IT_0167;
    const ccomplex_t IT_0169 = m_C_1*IT_0168;
    const ccomplex_t IT_0170 = (-4)*IT_0169;
    const ccomplex_t IT_0171 = 6*IT_0170;
    const ccomplex_t IT_0172 = (-6)*conj(IT_0018);
    const ccomplex_t IT_0173 = conj(IT_0138)*IT_0139;
    const ccomplex_t IT_0174 = conj(IT_0141)*IT_0142;
    const ccomplex_t IT_0175 = pow(s_14, 2);
    const ccomplex_t IT_0176 = IT_0063*IT_0175;
    const ccomplex_t IT_0177 = -IT_0001;
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = s_13*IT_0178;
    const ccomplex_t IT_0180 = 6*IT_0179;
    const ccomplex_t IT_0181 = IT_0019*IT_0178;
    const ccomplex_t IT_0182 = 6*IT_0181;
    const ccomplex_t IT_0183 = s_14*s_34*IT_0063;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = s_13 + IT_0184;
    const ccomplex_t IT_0186 = m_C_1*IT_0185;
    const ccomplex_t IT_0187 = 6*IT_0186;
    const ccomplex_t IT_0188 = (-6)*IT_0179;
    const ccomplex_t IT_0189 = s_13*IT_0185;
    const ccomplex_t IT_0190 = (-6)*IT_0189;
    const ccomplex_t IT_0191 = (-6)*IT_0181;
    const ccomplex_t IT_0192 = IT_0019*IT_0185;
    const ccomplex_t IT_0193 = (-6)*IT_0192;
    const ccomplex_t IT_0194 = m_s*IT_0178;
    const ccomplex_t IT_0195 = (-6)*IT_0194;
    const ccomplex_t IT_0196 = IT_0082*(conj(IT_0073)*IT_0139 + IT_0062
      *IT_0172 + conj(IT_0082)*IT_0180 + conj(IT_0138)*IT_0182 + IT_0159*IT_0187
       + conj(IT_0058)*IT_0188 + conj(IT_0083)*IT_0190 + conj(IT_0133)*IT_0191 +
       conj(IT_0141)*IT_0193 + IT_0125*IT_0195);
    const ccomplex_t IT_0197 = 6*IT_0189;
    const ccomplex_t IT_0198 = 6*IT_0192;
    const ccomplex_t IT_0199 = (-6)*IT_0186;
    const ccomplex_t IT_0200 = 6*conj(IT_0018);
    const ccomplex_t IT_0201 = 6*IT_0194;
    const ccomplex_t IT_0202 = IT_0058*(conj(IT_0073)*IT_0134 + conj(IT_0058)
      *IT_0180 + conj(IT_0133)*IT_0182 + conj(IT_0082)*IT_0188 + conj(IT_0138)
      *IT_0191 + conj(IT_0083)*IT_0197 + conj(IT_0141)*IT_0198 + IT_0159*IT_0199
       + IT_0062*IT_0200 + IT_0125*IT_0201);
    const ccomplex_t IT_0203 = pow(s_34, 2);
    const ccomplex_t IT_0204 = IT_0063*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0000 + IT_0205;
    const ccomplex_t IT_0207 = m_C_1*IT_0206;
    const ccomplex_t IT_0208 = 6*IT_0207;
    const ccomplex_t IT_0209 = s_13*IT_0206;
    const ccomplex_t IT_0210 = (-6)*IT_0209;
    const ccomplex_t IT_0211 = IT_0019*IT_0206;
    const ccomplex_t IT_0212 = (-6)*IT_0211;
    const ccomplex_t IT_0213 = m_s*IT_0185;
    const ccomplex_t IT_0214 = 6*IT_0213;
    const ccomplex_t IT_0215 = IT_0083*(conj(IT_0073)*IT_0142 + IT_0087
      *IT_0172 + conj(IT_0082)*IT_0190 + conj(IT_0138)*IT_0193 + conj(IT_0058)
      *IT_0197 + conj(IT_0133)*IT_0198 + IT_0159*IT_0208 + conj(IT_0083)*IT_0210
       + conj(IT_0141)*IT_0212 + IT_0125*IT_0214);
    const ccomplex_t IT_0216 = (-6)*IT_0131;
    const ccomplex_t IT_0217 = 0.5*IT_0166;
    const ccomplex_t IT_0218 = s_34 + IT_0217;
    const ccomplex_t IT_0219 = m_C_1*IT_0218;
    const ccomplex_t IT_0220 = 2*IT_0219;
    const ccomplex_t IT_0221 = (-6)*IT_0220;
    const ccomplex_t IT_0222 = (-6)*IT_0145;
    const ccomplex_t IT_0223 = IT_0133*IT_0201;
    const ccomplex_t IT_0224 = IT_0141*IT_0214;
    const ccomplex_t IT_0225 = IT_0088*IT_0220;
    const ccomplex_t IT_0226 = IT_0138*IT_0195;
    const ccomplex_t IT_0227 = (-2)*IT_0068;
    const ccomplex_t IT_0228 = IT_0070 + IT_0227;
    const ccomplex_t IT_0229 = IT_0063*IT_0228;
    const ccomplex_t IT_0230 = (-6)*IT_0229;
    const ccomplex_t IT_0231 = 12*s_13;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0154*IT_0232;
    const ccomplex_t IT_0234 = conj(IT_0133)*IT_0201;
    const ccomplex_t IT_0235 = conj(IT_0141)*IT_0214;
    const ccomplex_t IT_0236 = conj(IT_0082)*IT_0187;
    const ccomplex_t IT_0237 = conj(IT_0083)*IT_0208;
    const ccomplex_t IT_0238 = (-18)*IT_0019;
    const ccomplex_t IT_0239 = conj(IT_0119)*IT_0238;
    const ccomplex_t IT_0240 = IT_0172*IT_0220;
    const ccomplex_t IT_0241 = conj(IT_0138)*IT_0195;
    const ccomplex_t IT_0242 = conj(IT_0058)*IT_0199;
    const ccomplex_t IT_0243 = conj(IT_0124)*IT_0238;
    const ccomplex_t IT_0244 = conj(IT_0073)*IT_0171;
    const ccomplex_t IT_0245 = conj(IT_0144)*IT_0164;
    const ccomplex_t IT_0246 = conj(IT_0033)*IT_0216;
    const ccomplex_t IT_0247 = conj(IT_0058)*IT_0201;
    const ccomplex_t IT_0248 = conj(IT_0083)*IT_0214;
    const ccomplex_t IT_0249 = IT_0131*IT_0172;
    const ccomplex_t IT_0250 = conj(IT_0082)*IT_0195;
    const ccomplex_t IT_0251 = conj(IT_0138)*IT_0187;
    const ccomplex_t IT_0252 = conj(IT_0133)*IT_0199;
    const ccomplex_t IT_0253 = conj(IT_0141)*IT_0208;
    const ccomplex_t IT_0254 = IT_0159*IT_0238;
    const ccomplex_t IT_0255 = conj(IT_0073)*IT_0164;
    const ccomplex_t IT_0256 = conj(IT_0144)*IT_0171;
    const ccomplex_t IT_0257 = conj(IT_0033)*IT_0221;
    const ccomplex_t IT_0258 = 6*IT_0018*(IT_0022*conj(IT_0033) + conj(IT_0058
      )*IT_0062 + 0.166666666666667*conj(IT_0018)*IT_0072 + 0.166666666666667
      *conj(IT_0073)*IT_0075) + (IT_0062*conj(IT_0082) + conj(IT_0083)*IT_0087)
      *IT_0088 + IT_0125*IT_0132 + conj(IT_0033)*(IT_0135 + IT_0140 + IT_0143) +
       conj(IT_0144)*(IT_0135 + IT_0140 + IT_0143 + IT_0088*IT_0145) + IT_0144*
      (conj(IT_0033)*IT_0075 + IT_0072*conj(IT_0144) + IT_0146 + conj(IT_0073)
      *IT_0147 + IT_0159*IT_0164 + IT_0125*IT_0171 + IT_0145*IT_0172 + IT_0173 +
       IT_0174) + IT_0138*(conj(IT_0138)*IT_0180 + conj(IT_0082)*IT_0182 +
       IT_0125*IT_0187 + conj(IT_0133)*IT_0188 + conj(IT_0141)*IT_0190 + conj
      (IT_0058)*IT_0191 + conj(IT_0083)*IT_0193) + IT_0196 + IT_0133*(conj
      (IT_0133)*IT_0180 + conj(IT_0058)*IT_0182 + conj(IT_0138)*IT_0188 + conj
      (IT_0082)*IT_0191 + conj(IT_0141)*IT_0197 + conj(IT_0083)*IT_0198 +
       IT_0125*IT_0199) + IT_0202 + IT_0141*(conj(IT_0138)*IT_0190 + conj
      (IT_0082)*IT_0193 + conj(IT_0133)*IT_0197 + conj(IT_0058)*IT_0198 +
       IT_0125*IT_0208 + conj(IT_0141)*IT_0210 + conj(IT_0083)*IT_0212) +
       IT_0215 + IT_0073*(IT_0072*conj(IT_0073) + conj(IT_0018)*IT_0075 + conj
      (IT_0058)*IT_0134 + conj(IT_0082)*IT_0139 + conj(IT_0083)*IT_0142 + conj
      (IT_0144)*IT_0147 + IT_0125*IT_0164 + IT_0159*IT_0171 + conj(IT_0033)
      *IT_0222) + IT_0033*(conj(IT_0033)*IT_0072 + IT_0075*conj(IT_0144) +
       IT_0146 + IT_0173 + IT_0174 + IT_0022*IT_0200 + IT_0159*IT_0216 + IT_0125
      *IT_0221 + conj(IT_0073)*IT_0222) + IT_0159*(IT_0223 + IT_0224 + IT_0225 +
       IT_0226 + IT_0233) + IT_0154*(IT_0234 + IT_0235 + IT_0236 + IT_0237 +
       IT_0239 + IT_0240 + IT_0241 + IT_0242 + IT_0243 + IT_0244 + IT_0245 +
       IT_0246) + IT_0158*(IT_0159*IT_0232 + IT_0234 + IT_0235 + IT_0236 +
       IT_0237 + IT_0239 + IT_0240 + IT_0241 + IT_0242 + IT_0243 + IT_0244 +
       IT_0245 + IT_0246) + (IT_0119 + IT_0124)*(IT_0125*IT_0232 + IT_0247 +
       IT_0248 + IT_0249 + IT_0250 + IT_0251 + IT_0252 + IT_0253 + IT_0254 +
       IT_0255 + IT_0256 + IT_0257);
    return create_ccomplex_return(IT_0258);
}

