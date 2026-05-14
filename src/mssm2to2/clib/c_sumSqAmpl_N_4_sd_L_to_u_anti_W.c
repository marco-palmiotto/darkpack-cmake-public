#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sd_L_to_u_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sd_L_to_u_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_u, 2);
    const ccomplex_t IT_0010 = pow(m_N_4, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0018 + (-0.5)*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = pow(m_sd_L, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0009 + IT_0025 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-0.5)*IT_0035);
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0009 + IT_0025 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0014 + 2*IT_0028 + 2*IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, -2);
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = -IT_0010;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = s_13*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0050 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0051 = IT_0001*IT_0050;
    const ccomplex_t IT_0052 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0053 = IT_0001*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + 0.5*IT_0054);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0049*IT_0056;
    const ccomplex_t IT_0058 = IT_0039*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = m_N_4*IT_0059;
    const ccomplex_t IT_0061 = pow(m_d, 2);
    const ccomplex_t IT_0062 = cpow((-2)*s_12 + -IT_0010 + -IT_0025 + IT_0061 
      + -reg_prop, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0064 = cos(beta);
    const ccomplex_t IT_0065 = cpow(IT_0064, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0001*IT_0004*IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0002;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = IT_0062*IT_0063*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0072 = IT_0037*IT_0049;
    const ccomplex_t IT_0073 = IT_0039*IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0076 = IT_0023*IT_0075;
    const ccomplex_t IT_0077 = IT_0026*IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0060 + -IT_0070 + -IT_0073 + -IT_0077;
    const ccomplex_t IT_0079 = m_u*IT_0047;
    const ccomplex_t IT_0080 = cos(theta_W);
    const ccomplex_t IT_0081 = cpow(IT_0080, -1);
    const ccomplex_t IT_0082 = conj(N_B4)*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_W4)*e_em;
    const ccomplex_t IT_0086 = IT_0001*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + 3*IT_0087);
    const ccomplex_t IT_0089 = 0.166666666666667*IT_0088;
    const ccomplex_t IT_0090 = IT_0011*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = (-2)*IT_0059 + -IT_0092;
    const ccomplex_t IT_0094 = m_u*m_N_4;
    const ccomplex_t IT_0095 = IT_0047*IT_0094;
    const ccomplex_t IT_0096 = -IT_0028 + -IT_0041;
    const ccomplex_t IT_0097 = pow(m_W, 2);
    const ccomplex_t IT_0098 = IT_0094*IT_0097;
    const ccomplex_t IT_0099 = 3*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0084 + (-3)*IT_0087);
    const ccomplex_t IT_0101 = 0.166666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0068*IT_0101;
    const ccomplex_t IT_0103 = IT_0062*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = s_14*IT_0097;
    const ccomplex_t IT_0106 = IT_0043*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = s_14 + IT_0107;
    const ccomplex_t IT_0109 = m_u*IT_0108;
    const ccomplex_t IT_0110 = 2*IT_0109;
    const ccomplex_t IT_0111 = pow(m_W, 4);
    const ccomplex_t IT_0112 = s_13*IT_0111;
    const ccomplex_t IT_0113 = IT_0043*IT_0112;
    const ccomplex_t IT_0114 = s_14*s_34;
    const ccomplex_t IT_0115 = 24*IT_0114;
    const ccomplex_t IT_0116 = s_13*IT_0097;
    const ccomplex_t IT_0117 = (-12)*IT_0116;
    const ccomplex_t IT_0118 = 6*IT_0113 + IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = (-3)*IT_0098;
    const ccomplex_t IT_0120 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0123 = IT_0001*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0075*IT_0126;
    const ccomplex_t IT_0128 = IT_0026*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = -conj(IT_0059);
    const ccomplex_t IT_0132 = conj(IT_0130) + IT_0131;
    const ccomplex_t IT_0133 = s_14*s_34*IT_0043;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = s_13 + IT_0134;
    const ccomplex_t IT_0136 = IT_0094*IT_0135;
    const ccomplex_t IT_0137 = -IT_0014;
    const ccomplex_t IT_0138 = pow(s_34, 2);
    const ccomplex_t IT_0139 = IT_0043*IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0009 + IT_0140;
    const ccomplex_t IT_0142 = s_13*IT_0141;
    const ccomplex_t IT_0143 = m_u*IT_0135;
    const ccomplex_t IT_0144 = (-6)*IT_0137;
    const ccomplex_t IT_0145 = s_13*IT_0135;
    const ccomplex_t IT_0146 = (-6)*IT_0042;
    const ccomplex_t IT_0147 = (-6)*conj(IT_0042);
    const ccomplex_t IT_0148 = m_N_4*IT_0135;
    const ccomplex_t IT_0149 = s_34*IT_0097;
    const ccomplex_t IT_0150 = IT_0043*IT_0149;
    const ccomplex_t IT_0151 = (-0.25)*IT_0150;
    const ccomplex_t IT_0152 = s_34 + IT_0151;
    const ccomplex_t IT_0153 = m_N_4*IT_0152;
    const ccomplex_t IT_0154 = (-4)*IT_0153;
    const ccomplex_t IT_0155 = (-6)*IT_0078;
    const ccomplex_t IT_0156 = IT_0092 + 2*IT_0104;
    const ccomplex_t IT_0157 = IT_0094*IT_0141;
    const ccomplex_t IT_0158 = (-6)*IT_0093;
    const ccomplex_t IT_0159 = (-6)*IT_0156;
    const ccomplex_t IT_0160 = m_N_4*IT_0141;
    const ccomplex_t IT_0161 = (-6)*IT_0113;
    const ccomplex_t IT_0162 = (-24)*IT_0114;
    const ccomplex_t IT_0163 = 12*IT_0116;
    const ccomplex_t IT_0164 = IT_0161 + IT_0162 + IT_0163;
    const ccomplex_t IT_0165 = IT_0115 + IT_0161;
    const ccomplex_t IT_0166 = (-6)*conj(IT_0078);
    const ccomplex_t IT_0167 = m_N_4*IT_0129;
    const ccomplex_t IT_0168 = m_N_4*IT_0041;
    const ccomplex_t IT_0169 = m_N_4*IT_0028;
    const ccomplex_t IT_0170 = m_u*IT_0104;
    const ccomplex_t IT_0171 = -IT_0168 + -IT_0169 + -IT_0170;
    const ccomplex_t IT_0172 = IT_0030*IT_0056;
    const ccomplex_t IT_0173 = IT_0039*IT_0071*IT_0172;
    const ccomplex_t IT_0174 = IT_0016*IT_0126;
    const ccomplex_t IT_0175 = IT_0026*IT_0074*IT_0174;
    const ccomplex_t IT_0176 = IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = (-2)*IT_0114;
    const ccomplex_t IT_0178 = IT_0116 + IT_0177;
    const ccomplex_t IT_0179 = IT_0043*IT_0178;
    const ccomplex_t IT_0180 = (-6)*IT_0179;
    const ccomplex_t IT_0181 = 12*s_13;
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0129;
    const ccomplex_t IT_0184 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = 6*IT_0096;
    const ccomplex_t IT_0187 = 0.5*IT_0150;
    const ccomplex_t IT_0188 = s_34 + IT_0187;
    const ccomplex_t IT_0189 = m_N_4*IT_0188;
    const ccomplex_t IT_0190 = 2*IT_0189;
    const ccomplex_t IT_0191 = 6*IT_0104;
    const ccomplex_t IT_0192 = 6*IT_0059 + (-6)*IT_0130;
    const ccomplex_t IT_0193 = 6*IT_0042*IT_0079 + IT_0143*IT_0144 + IT_0159
      *IT_0160 + IT_0094*((-18)*IT_0171 + (-18)*IT_0176) + (IT_0078 + IT_0167)
      *IT_0182 + 6*IT_0148*IT_0185 + IT_0110*IT_0186 + IT_0190*IT_0191 + IT_0154
      *IT_0192;
    const ccomplex_t IT_0194 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0195 = (-6)*conj(IT_0093);
    const ccomplex_t IT_0196 = 0.166666666666667*IT_0195;
    const ccomplex_t IT_0197 = conj(IT_0183) + IT_0196;
    const ccomplex_t IT_0198 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0199 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0200 = 0.166666666666667*IT_0148;
    const ccomplex_t IT_0201 = s_14*IT_0009;
    const ccomplex_t IT_0202 = s_13*s_34;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = conj(IT_0104) + conj(IT_0130);
    const ccomplex_t IT_0206 = (-6)*conj(IT_0156);
    const ccomplex_t IT_0207 = IT_0104 + IT_0130;
    const ccomplex_t IT_0208 = 0.166666666666667*IT_0207;
    const ccomplex_t IT_0209 = 6*IT_0113 + IT_0162;
    const ccomplex_t IT_0210 = 0.166666666666667*conj(IT_0059);
    const ccomplex_t IT_0211 = 0.166666666666667*conj(IT_0104);
    const ccomplex_t IT_0212 = IT_0171 + IT_0176;
    const ccomplex_t IT_0213 = (-0.25)*IT_0106;
    const ccomplex_t IT_0214 = s_14 + IT_0213;
    const ccomplex_t IT_0215 = m_u*IT_0214;
    const ccomplex_t IT_0216 = (-4)*IT_0215;
    const ccomplex_t IT_0217 = (-6)*conj(IT_0183);
    const ccomplex_t IT_0218 = 0.166666666666667*IT_0217;
    const ccomplex_t IT_0219 = conj(IT_0093) + IT_0218;
    const ccomplex_t IT_0220 = (-18)*conj(IT_0078);
    const ccomplex_t IT_0221 = -IT_0132;
    const ccomplex_t IT_0222 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0223 = 0.166666666666667*IT_0160;
    const ccomplex_t IT_0224 = (-6)*IT_0183;
    const ccomplex_t IT_0225 = 0.166666666666667*IT_0224;
    const ccomplex_t IT_0226 = IT_0146*IT_0148 + IT_0143*IT_0159 + IT_0144
      *IT_0160 + IT_0094*((-18)*IT_0078 + (-18)*IT_0167) + IT_0154*IT_0186 +
       IT_0110*IT_0192 + IT_0182*IT_0212 + IT_0191*IT_0216 + 6*IT_0079*(IT_0093 
      + IT_0225);
    const ccomplex_t IT_0227 = s_13*s_14;
    const ccomplex_t IT_0228 = s_34*IT_0010;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = IT_0227 + IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0042)*IT_0079;
    const ccomplex_t IT_0232 = conj(IT_0096)*IT_0110;
    const ccomplex_t IT_0233 = 6*IT_0042*(conj(IT_0042)*IT_0048 + conj(IT_0078
      )*IT_0079 + conj(IT_0093)*IT_0095) + 6*IT_0096*(IT_0099*conj(IT_0104) +
       conj(IT_0078)*IT_0110 + 0.166666666666667*conj(IT_0096)*IT_0118 + 
      -IT_0119*IT_0132) + (conj(IT_0093)*IT_0136 + conj(IT_0137)*IT_0142 + conj
      (IT_0078)*IT_0143)*IT_0144 + IT_0145*(conj(IT_0137)*IT_0146 + IT_0137
      *IT_0147) + (conj(IT_0093)*IT_0148 + conj(IT_0130)*IT_0154)*IT_0155 + conj
      (IT_0156)*(IT_0136*IT_0146 + IT_0144*IT_0157 + IT_0145*IT_0158 + IT_0142
      *IT_0159 + IT_0155*IT_0160) + (-6)*IT_0130*(conj(IT_0096)*IT_0119 + (
      -0.166666666666667)*IT_0118*conj(IT_0130) + (-0.166666666666667)*conj
      (IT_0059)*IT_0164 + (-0.166666666666667)*conj(IT_0104)*IT_0165 + (
      -0.166666666666667)*IT_0154*IT_0166) + conj(IT_0167)*IT_0193 + 6*conj
      (IT_0183)*(IT_0136*IT_0137 + 0.166666666666667*IT_0048*IT_0158 + IT_0146
      *IT_0194) + 6*IT_0183*(IT_0136*conj(IT_0137) + conj(IT_0078)*IT_0148 +
       IT_0145*conj(IT_0156) + IT_0147*IT_0194 + IT_0048*IT_0197) + 6*IT_0093*
      (IT_0048*conj(IT_0093) + conj(IT_0042)*IT_0095 + IT_0198*IT_0199 + IT_0166
      *IT_0200) + 6*IT_0204*(conj(IT_0096)*IT_0137 + IT_0096*conj(IT_0137) +
       conj(IT_0059)*IT_0156 + IT_0059*conj(IT_0156) + 0.166666666666667*IT_0159
      *IT_0205 + IT_0206*IT_0208) + IT_0191*(conj(IT_0096)*IT_0099 +
       0.166666666666667*conj(IT_0130)*IT_0165 + conj(IT_0078)*IT_0190 + IT_0209
      *IT_0210 + IT_0118*IT_0211) + 6*IT_0059*(conj(IT_0096)*IT_0119 + conj
      (IT_0078)*IT_0154 + 0.166666666666667*conj(IT_0130)*IT_0164 + IT_0118
      *IT_0210 + IT_0209*IT_0211) + 6*IT_0156*(IT_0145*IT_0197 +
       0.166666666666667*IT_0157*IT_0198 + IT_0147*IT_0199 + IT_0166*IT_0223) +
       6*IT_0212*(conj(IT_0096)*IT_0154 + IT_0147*IT_0200 + conj(IT_0104)
      *IT_0216 + IT_0079*IT_0219 + 0.166666666666667*IT_0094*IT_0220 + IT_0110
      *IT_0221 + IT_0206*IT_0222 + IT_0198*IT_0223) + (conj(IT_0171) + conj
      (IT_0176))*IT_0226 + 6*(conj(IT_0042)*IT_0096 + IT_0042*conj(IT_0096) +
       IT_0185*IT_0205 + conj(IT_0183)*IT_0207 + IT_0195*IT_0208 + IT_0059
      *IT_0219 + conj(IT_0059)*(IT_0093 + IT_0225))*IT_0230 + 6*IT_0078*(conj
      (IT_0059)*IT_0154 + 0.166666666666667*conj(IT_0078)*IT_0182 + IT_0148*conj
      (IT_0183) + conj(IT_0104)*IT_0190 + IT_0198*IT_0222 + IT_0231 + IT_0232) +
       6*IT_0167*(0.166666666666667*conj(IT_0078)*IT_0182 + conj(IT_0104)
      *IT_0190 + IT_0148*IT_0197 + IT_0154*IT_0221 + IT_0198*IT_0222 + IT_0206
      *IT_0223 + IT_0231 + IT_0232);
    return create_ccomplex_return(IT_0233);
}

