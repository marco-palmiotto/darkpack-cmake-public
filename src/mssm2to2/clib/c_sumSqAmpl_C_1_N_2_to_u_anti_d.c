#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_2_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_2_to_u_anti_d(
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
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_N_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 2*IT_0012;
    const ccomplex_t IT_0014 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0009*IT_0018;
    const ccomplex_t IT_0020 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = IT_0016 + IT_0022;
    const ccomplex_t IT_0024 = IT_0013 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = cpow(IT_0009, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0004*IT_0006*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0018*IT_0031;
    const ccomplex_t IT_0036 = IT_0020*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0031;
    const ccomplex_t IT_0039 = IT_0014*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = IT_0037 + IT_0040;
    const ccomplex_t IT_0042 = IT_0011*IT_0038;
    const ccomplex_t IT_0043 = 2*IT_0042;
    const ccomplex_t IT_0044 = IT_0041 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0034*IT_0046;
    const ccomplex_t IT_0048 = IT_0030*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0029 + IT_0049;
    const ccomplex_t IT_0051 = s_13*s_24;
    const ccomplex_t IT_0052 = s_12*s_34;
    const ccomplex_t IT_0053 = s_14*s_23;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = IT_0051 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0004*IT_0006*IT_0032;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = N_B2*e_em;
    const ccomplex_t IT_0060 = IT_0018*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = N_W2*e_em;
    const ccomplex_t IT_0063 = IT_0006*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-3)*IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = pow(m_u, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0001 + IT_0068 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = 6*conj(IT_0071);
    const ccomplex_t IT_0073 = m_u*s_24*m_C_1;
    const ccomplex_t IT_0074 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0075 = IT_0006*IT_0074;
    const ccomplex_t IT_0076 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0077 = IT_0006*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + 0.5*IT_0078);
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0006;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0080*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0031, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0004*IT_0006*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0058*IT_0089;
    const ccomplex_t IT_0091 = IT_0069*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0004*IT_0006*IT_0087;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = conj(N_B2)*e_em;
    const ccomplex_t IT_0097 = IT_0018*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W2)*e_em;
    const ccomplex_t IT_0100 = IT_0006*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + 3*IT_0101);
    const ccomplex_t IT_0103 = 0.166666666666667*IT_0102;
    const ccomplex_t IT_0104 = IT_0095*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0002 + IT_0068 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 6*conj(IT_0107);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0110 = IT_0066*IT_0109;
    const ccomplex_t IT_0111 = IT_0069*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0004*IT_0006*IT_0087;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0095*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_13 + IT_0001 + IT_0068 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = -IT_0112 + -IT_0118;
    const ccomplex_t IT_0120 = 12*conj(IT_0071);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0018;
    const ccomplex_t IT_0122 = 0.666666666666667*IT_0121;
    const ccomplex_t IT_0123 = IT_0058*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0002 + IT_0068 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0006;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0004*IT_0006*IT_0032;
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = IT_0105*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0126 + IT_0133;
    const ccomplex_t IT_0135 = IT_0050 + IT_0134;
    const ccomplex_t IT_0136 = 12*conj(IT_0086);
    const ccomplex_t IT_0137 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0138 = IT_0010*IT_0137;
    const ccomplex_t IT_0139 = (-2)*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0141 = IT_0019*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0144 = IT_0010*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = IT_0139 + IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0034*IT_0149;
    const ccomplex_t IT_0151 = IT_0030*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0038*IT_0137;
    const ccomplex_t IT_0154 = (-2)*IT_0153;
    const ccomplex_t IT_0155 = IT_0035*IT_0140;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = IT_0038*IT_0143;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = IT_0154 + IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0008*IT_0162;
    const ccomplex_t IT_0164 = IT_0003*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0152 + -IT_0165;
    const ccomplex_t IT_0167 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0168 = IT_0006*IT_0167;
    const ccomplex_t IT_0169 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0170 = IT_0006*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + (-0.5)*IT_0171);
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = IT_0082*IT_0173;
    const ccomplex_t IT_0175 = IT_0084*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = 12*conj(IT_0176);
    const ccomplex_t IT_0178 = 12*IT_0071;
    const ccomplex_t IT_0179 = 24*conj(IT_0086);
    const ccomplex_t IT_0180 = m_d*s_13*m_N_2;
    const ccomplex_t IT_0181 = 12*IT_0107;
    const ccomplex_t IT_0182 = 12*conj(IT_0107);
    const ccomplex_t IT_0183 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0004*IT_0006*IT_0009*IT_0087;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0149*IT_0184;
    const ccomplex_t IT_0186 = IT_0030*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0004*IT_0006;
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0162*IT_0189;
    const ccomplex_t IT_0191 = IT_0003*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = -IT_0187 + -IT_0192;
    const ccomplex_t IT_0194 = IT_0046*IT_0184;
    const ccomplex_t IT_0195 = IT_0030*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = IT_0026*IT_0189;
    const ccomplex_t IT_0198 = IT_0003*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = IT_0196 + IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0018;
    const ccomplex_t IT_0202 = (-0.333333333333333)*IT_0201;
    const ccomplex_t IT_0203 = IT_0095*IT_0202;
    const ccomplex_t IT_0204 = IT_0116*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0089*IT_0109;
    const ccomplex_t IT_0207 = IT_0069*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0205 + IT_0208;
    const ccomplex_t IT_0210 = IT_0103*IT_0128;
    const ccomplex_t IT_0211 = IT_0105*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0004*IT_0006*IT_0032;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0058*IT_0214;
    const ccomplex_t IT_0216 = IT_0124*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = -IT_0212 + -IT_0217;
    const ccomplex_t IT_0219 = 24*conj(IT_0176);
    const ccomplex_t IT_0220 = m_d*s_23*m_C_1;
    const ccomplex_t IT_0221 = IT_0095*IT_0130;
    const ccomplex_t IT_0222 = IT_0105*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = (-12)*conj(IT_0086);
    const ccomplex_t IT_0226 = (-12)*conj(IT_0176);
    const ccomplex_t IT_0227 = m_u*s_14*m_N_2;
    const ccomplex_t IT_0228 = s_34*m_C_1*m_N_2;
    const ccomplex_t IT_0229 = 24*IT_0086;
    const ccomplex_t IT_0230 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0231 = m_d*m_u*s_12;
    const ccomplex_t IT_0232 = -IT_0051;
    const ccomplex_t IT_0233 = IT_0052 + IT_0053 + IT_0232;
    const ccomplex_t IT_0234 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0235 = IT_0055 + IT_0232;
    const ccomplex_t IT_0236 = m_d*m_u*m_C_1*m_N_2;
    const ccomplex_t IT_0237 = 6*IT_0071;
    const ccomplex_t IT_0238 = 24*IT_0176;
    const ccomplex_t IT_0239 = (-6)*IT_0107;
    const ccomplex_t IT_0240 = IT_0093*IT_0179 + conj(IT_0134)*IT_0181 +
       IT_0134*IT_0182 + IT_0072*(IT_0107 + IT_0193) + IT_0120*IT_0209 + IT_0178
      *conj(IT_0209) + IT_0219*IT_0224 + conj(IT_0093)*IT_0229 + IT_0050*IT_0234
       + (conj(IT_0107) + conj(IT_0193))*IT_0237 + conj(IT_0224)*IT_0238 + conj
      (IT_0050)*IT_0239;
    const ccomplex_t IT_0241 = 6*IT_0056;
    const ccomplex_t IT_0242 = 6*IT_0220;
    const ccomplex_t IT_0243 = 6*IT_0236;
    const ccomplex_t IT_0244 = 6*IT_0227;
    const ccomplex_t IT_0245 = 12*IT_0176;
    const ccomplex_t IT_0246 = 12*IT_0051;
    const ccomplex_t IT_0247 = 12*IT_0073;
    const ccomplex_t IT_0248 = 12*IT_0180;
    const ccomplex_t IT_0249 = (-6)*IT_0231;
    const ccomplex_t IT_0250 = (-6)*IT_0228;
    const ccomplex_t IT_0251 = IT_0180*IT_0229 + IT_0235*IT_0239 + IT_0193
      *IT_0241 + IT_0218*IT_0242 + IT_0135*IT_0243 + IT_0224*IT_0244 + IT_0220
      *IT_0245 + IT_0209*IT_0246 + IT_0093*IT_0247 + IT_0119*IT_0248 + IT_0166
      *IT_0249 + IT_0200*IT_0250;
    const ccomplex_t IT_0252 = 6*IT_0073;
    const ccomplex_t IT_0253 = 6*IT_0180;
    const ccomplex_t IT_0254 = 12*IT_0086;
    const ccomplex_t IT_0255 = 12*IT_0052;
    const ccomplex_t IT_0256 = 12*IT_0236;
    const ccomplex_t IT_0257 = (-12)*IT_0176;
    const ccomplex_t IT_0258 = (-12)*IT_0231;
    const ccomplex_t IT_0259 = (-12)*IT_0228;
    const ccomplex_t IT_0260 = (-6)*IT_0220;
    const ccomplex_t IT_0261 = (-6)*IT_0236;
    const ccomplex_t IT_0262 = (-6)*IT_0227;
    const ccomplex_t IT_0263 = IT_0233*IT_0239 + IT_0209*IT_0241 + IT_0093
      *IT_0252 + IT_0119*IT_0253 + IT_0180*IT_0254 + IT_0193*IT_0255 + IT_0050
      *IT_0256 + IT_0220*IT_0257 + IT_0166*IT_0258 + IT_0200*IT_0259 + IT_0218
      *IT_0260 + IT_0134*IT_0261 + IT_0224*IT_0262;
    const ccomplex_t IT_0264 = 6*IT_0107;
    const ccomplex_t IT_0265 = 6*IT_0231;
    const ccomplex_t IT_0266 = 6*IT_0228;
    const ccomplex_t IT_0267 = IT_0073*IT_0178 + IT_0051*IT_0229 + IT_0228
      *IT_0245 + IT_0119*IT_0246 + IT_0209*IT_0248 + IT_0135*IT_0252 + IT_0193
      *IT_0253 + IT_0093*IT_0256 + IT_0200*IT_0260 + IT_0166*IT_0262 + IT_0180
      *IT_0264 + IT_0224*IT_0265 + IT_0218*IT_0266;
    const ccomplex_t IT_0268 = 12*IT_0053;
    const ccomplex_t IT_0269 = 12*IT_0220;
    const ccomplex_t IT_0270 = IT_0181*IT_0227 + IT_0220*IT_0237 + IT_0209
      *IT_0244 + IT_0200*IT_0252 + IT_0166*IT_0253 + IT_0231*IT_0254 + IT_0218
      *IT_0256 + IT_0050*IT_0260 + IT_0193*IT_0262 + IT_0119*IT_0265 + IT_0093
      *IT_0266 + IT_0224*IT_0268 + IT_0134*IT_0269;
    const ccomplex_t IT_0271 = 12*IT_0227;
    const ccomplex_t IT_0272 = IT_0181*IT_0220 + IT_0227*IT_0237 + IT_0053
      *IT_0238 + IT_0209*IT_0242 + IT_0166*IT_0252 + IT_0200*IT_0253 + IT_0228
      *IT_0254 + IT_0224*IT_0256 + IT_0193*IT_0260 + IT_0050*IT_0262 + IT_0093
      *IT_0265 + IT_0119*IT_0266 + IT_0218*IT_0268 + IT_0134*IT_0271;
    const ccomplex_t IT_0273 = IT_0178*IT_0180 + IT_0231*IT_0245 + IT_0093
      *IT_0246 + IT_0209*IT_0247 + IT_0193*IT_0252 + IT_0135*IT_0253 + IT_0119
      *IT_0256 + IT_0166*IT_0260 + IT_0200*IT_0262 + IT_0073*IT_0264 + IT_0218
      *IT_0265 + IT_0224*IT_0266;
    const ccomplex_t IT_0274 = (-6)*IT_0071;
    const ccomplex_t IT_0275 = (-12)*IT_0086;
    const ccomplex_t IT_0276 = IT_0180*IT_0245 + IT_0209*IT_0250 + IT_0224
      *IT_0252 + IT_0218*IT_0253 + IT_0200*IT_0255 + IT_0166*IT_0256 + IT_0050
      *IT_0258 + IT_0193*IT_0259 + IT_0119*IT_0260 + IT_0093*IT_0262 + IT_0228
      *IT_0264 + IT_0134*IT_0265 + IT_0231*IT_0274 + IT_0220*IT_0275;
    const ccomplex_t IT_0277 = IT_0073*IT_0245 + IT_0209*IT_0249 + IT_0218
      *IT_0252 + IT_0224*IT_0253 + IT_0166*IT_0255 + IT_0200*IT_0256 + IT_0193
      *IT_0258 + IT_0050*IT_0259 + IT_0093*IT_0260 + IT_0119*IT_0262 + IT_0231
      *IT_0264 + IT_0134*IT_0266 + IT_0228*IT_0274 + IT_0227*IT_0275;
    const ccomplex_t IT_0278 = IT_0073*IT_0254;
    const ccomplex_t IT_0279 = (-6)*IT_0233;
    const ccomplex_t IT_0280 = IT_0119*IT_0252;
    const ccomplex_t IT_0281 = IT_0093*IT_0253;
    const ccomplex_t IT_0282 = IT_0209*IT_0243;
    const ccomplex_t IT_0283 = IT_0227*IT_0238 + IT_0193*IT_0261 + IT_0200
      *IT_0265 + IT_0166*IT_0266 + IT_0134*IT_0268 + IT_0224*IT_0269 + IT_0218
      *IT_0271 + IT_0235*IT_0274 + IT_0278 + IT_0050*IT_0279 + IT_0280 + IT_0281
       + IT_0282;
    const ccomplex_t IT_0284 = IT_0056*IT_0237 + IT_0050*IT_0255 + IT_0193
      *IT_0256 + IT_0227*IT_0257 + IT_0200*IT_0258 + IT_0166*IT_0259 + IT_0224
      *IT_0260 + IT_0218*IT_0262 + IT_0278 + IT_0134*IT_0279 + IT_0280 + IT_0281
       + IT_0282;
    const ccomplex_t IT_0285 = IT_0050*IT_0056*IT_0072 + IT_0073*(24*conj
      (IT_0071)*IT_0086 + 24*IT_0071*conj(IT_0086) + IT_0093*IT_0108 + IT_0119
      *IT_0120 + IT_0135*IT_0136 + IT_0166*IT_0177) + IT_0051*(48*IT_0086*conj
      (IT_0086) + conj(IT_0071)*IT_0178 + IT_0119*IT_0179) + IT_0180*(IT_0108
      *IT_0119 + IT_0093*IT_0120 + conj(IT_0086)*IT_0181 + IT_0086*IT_0182 +
       IT_0136*IT_0193 + IT_0177*IT_0200 + IT_0179*IT_0209) + IT_0053*(48
      *IT_0176*conj(IT_0176) + conj(IT_0107)*IT_0181 + IT_0218*IT_0219) +
       IT_0220*(24*conj(IT_0107)*IT_0176 + 24*IT_0107*conj(IT_0176) + IT_0177
      *IT_0209 + IT_0182*IT_0218 + IT_0072*IT_0224 + IT_0200*IT_0225 + IT_0193
      *IT_0226) + (IT_0120*IT_0176 + conj(IT_0176)*IT_0178 + IT_0072*IT_0218 +
       IT_0134*IT_0219 + IT_0182*IT_0224 + IT_0166*IT_0225 + IT_0050*IT_0226)
      *IT_0227 + IT_0228*(IT_0119*IT_0177 + IT_0176*IT_0179 + IT_0108*IT_0200 +
       IT_0136*IT_0218 + conj(IT_0176)*IT_0229 + IT_0166*IT_0230) + (IT_0108
      *IT_0166 + IT_0093*IT_0177 + IT_0136*IT_0224 + IT_0200*IT_0230)*IT_0231 +
       IT_0193*IT_0233*IT_0234 + (IT_0134*IT_0230 + IT_0209*IT_0234)*IT_0235 +
       IT_0236*IT_0240 + conj(IT_0209)*IT_0251 + conj(IT_0193)*IT_0263 + conj
      (IT_0119)*IT_0267 + conj(IT_0224)*IT_0270 + conj(IT_0218)*IT_0272 + conj
      (IT_0093)*IT_0273 + conj(IT_0200)*IT_0276 + conj(IT_0166)*IT_0277 + conj
      (IT_0134)*IT_0283 + conj(IT_0050)*IT_0284;
    return create_ccomplex_return(IT_0285);
}

