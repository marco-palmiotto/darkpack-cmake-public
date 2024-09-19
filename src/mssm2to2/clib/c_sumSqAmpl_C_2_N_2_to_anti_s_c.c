#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_anti_s_c(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = s_14*s_23;
    const ccomplex_t IT_0002 = s_13*s_24;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0000 + IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = pow(m_s, 2);
    const ccomplex_t IT_0006 = pow(m_N_2, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0016;
    const ccomplex_t IT_0018 = (-0.333333333333333)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0012;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0021 + IT_0028;
    const ccomplex_t IT_0030 = pow(m_C_2, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + -IT_0006 + -IT_0030 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0032 = sin(beta);
    const ccomplex_t IT_0033 = IT_0012*IT_0032;
    const ccomplex_t IT_0034 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (-2)*IT_0035;
    const ccomplex_t IT_0037 = IT_0016*IT_0032;
    const ccomplex_t IT_0038 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = IT_0040 + IT_0043;
    const ccomplex_t IT_0045 = IT_0036 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0008*IT_0010*IT_0012*IT_0032;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = IT_0031*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = pow(m_W, 2);
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + -IT_0006 + -IT_0030 + IT_0053 
      + -reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0008*IT_0012;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0009*IT_0012;
    const ccomplex_t IT_0058 = IT_0034*IT_0057;
    const ccomplex_t IT_0059 = (-2)*IT_0058;
    const ccomplex_t IT_0060 = IT_0009*IT_0016;
    const ccomplex_t IT_0061 = IT_0038*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = IT_0041*IT_0057;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = IT_0059 + IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0056*IT_0068;
    const ccomplex_t IT_0070 = IT_0054*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0052 + -IT_0071;
    const ccomplex_t IT_0073 = 6*IT_0072;
    const ccomplex_t IT_0074 = 6*conj(IT_0072);
    const ccomplex_t IT_0075 = cpow(IT_0032, -1);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0008*IT_0012*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = N_B2*e_em;
    const ccomplex_t IT_0079 = IT_0016*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = N_W2*e_em;
    const ccomplex_t IT_0082 = IT_0012*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + (-3)*IT_0083);
    const ccomplex_t IT_0085 = 0.166666666666667*IT_0084;
    const ccomplex_t IT_0086 = IT_0077*IT_0085;
    const ccomplex_t IT_0087 = IT_0026*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0008*IT_0012;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0092 = IT_0033*IT_0091;
    const ccomplex_t IT_0093 = 2*IT_0092;
    const ccomplex_t IT_0094 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0095 = IT_0037*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0098 = IT_0033*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = IT_0096 + IT_0099;
    const ccomplex_t IT_0101 = IT_0093 + IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0090*IT_0103;
    const ccomplex_t IT_0105 = IT_0054*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0057*IT_0091;
    const ccomplex_t IT_0108 = 2*IT_0107;
    const ccomplex_t IT_0109 = IT_0060*IT_0094;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = IT_0057*IT_0097;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = IT_0110 + IT_0112;
    const ccomplex_t IT_0114 = IT_0108 + IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0008*IT_0009*IT_0012*IT_0075;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = IT_0031*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = IT_0106 + IT_0121;
    const ccomplex_t IT_0123 = 6*IT_0122;
    const ccomplex_t IT_0124 = 6*conj(IT_0122);
    const ccomplex_t IT_0125 = IT_0047*IT_0118;
    const ccomplex_t IT_0126 = IT_0031*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0068*IT_0090;
    const ccomplex_t IT_0129 = IT_0054*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = -IT_0127 + -IT_0130;
    const ccomplex_t IT_0132 = 12*IT_0131;
    const ccomplex_t IT_0133 = 12*IT_0072;
    const ccomplex_t IT_0134 = 12*IT_0122;
    const ccomplex_t IT_0135 = IT_0049*IT_0116;
    const ccomplex_t IT_0136 = IT_0031*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0056*IT_0103;
    const ccomplex_t IT_0139 = IT_0054*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0137 + IT_0140;
    const ccomplex_t IT_0142 = 12*IT_0141;
    const ccomplex_t IT_0143 = m_c*s_23*m_C_2;
    const ccomplex_t IT_0144 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0145 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0146 = IT_0012*IT_0145;
    const ccomplex_t IT_0147 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0148 = IT_0012*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + (-0.5)*IT_0149);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0012;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = IT_0151*IT_0153;
    const ccomplex_t IT_0155 = IT_0144*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0024*IT_0077;
    const ccomplex_t IT_0158 = IT_0026*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_13 + IT_0005 + IT_0030 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0162 = conj(N_B2)*e_em;
    const ccomplex_t IT_0163 = IT_0016*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W2)*e_em;
    const ccomplex_t IT_0166 = IT_0012*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0164 + 3*IT_0167);
    const ccomplex_t IT_0169 = 0.166666666666667*IT_0168;
    const ccomplex_t IT_0170 = IT_0014*IT_0169;
    const ccomplex_t IT_0171 = IT_0161*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = 0.166666666666667*IT_0073 + IT_0172;
    const ccomplex_t IT_0174 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0175 = conj(IT_0172) + IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = IT_0014*IT_0177;
    const ccomplex_t IT_0179 = IT_0007*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0022*IT_0085;
    const ccomplex_t IT_0182 = IT_0026*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = -IT_0180 + -IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0016;
    const ccomplex_t IT_0186 = 0.666666666666667*IT_0185;
    const ccomplex_t IT_0187 = IT_0077*IT_0186;
    const ccomplex_t IT_0188 = cpow((-2)*s_13 + IT_0005 + IT_0030 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0012;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0008*IT_0012*IT_0075;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0192*IT_0194;
    const ccomplex_t IT_0196 = IT_0161*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = IT_0190 + IT_0197;
    const ccomplex_t IT_0199 = 6*IT_0198;
    const ccomplex_t IT_0200 = IT_0123 + IT_0199;
    const ccomplex_t IT_0201 = 6*conj(IT_0198);
    const ccomplex_t IT_0202 = IT_0124 + IT_0201;
    const ccomplex_t IT_0203 = IT_0169*IT_0192;
    const ccomplex_t IT_0204 = IT_0161*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0008*IT_0012*IT_0075;
    const ccomplex_t IT_0207 = 0.5*IT_0206;
    const ccomplex_t IT_0208 = IT_0077*IT_0207;
    const ccomplex_t IT_0209 = IT_0188*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = -IT_0205 + -IT_0210;
    const ccomplex_t IT_0212 = 6*IT_0131;
    const ccomplex_t IT_0213 = 6*conj(IT_0131);
    const ccomplex_t IT_0214 = IT_0014*IT_0194;
    const ccomplex_t IT_0215 = IT_0161*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = 6*IT_0141;
    const ccomplex_t IT_0219 = 6*conj(IT_0141);
    const ccomplex_t IT_0220 = 12*IT_0160;
    const ccomplex_t IT_0221 = 12*conj(IT_0160);
    const ccomplex_t IT_0222 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0223 = IT_0012*IT_0222;
    const ccomplex_t IT_0224 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0225 = IT_0012*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + 0.5*IT_0226);
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = IT_0153*IT_0228;
    const ccomplex_t IT_0230 = IT_0144*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = 12*IT_0198;
    const ccomplex_t IT_0233 = 12*conj(IT_0198);
    const ccomplex_t IT_0234 = 12*conj(IT_0122);
    const ccomplex_t IT_0235 = 12*IT_0184;
    const ccomplex_t IT_0236 = 12*conj(IT_0184);
    const ccomplex_t IT_0237 = 12*conj(IT_0131)*IT_0156 + IT_0132*conj(IT_0156
      ) + 6*conj(IT_0160)*IT_0173 + 6*IT_0160*IT_0175 + conj(IT_0184)*IT_0200 +
       IT_0184*IT_0202 + conj(IT_0211)*IT_0212 + IT_0211*IT_0213 + conj(IT_0217)
      *IT_0218 + IT_0217*IT_0219 + conj(IT_0029)*IT_0220 + IT_0029*IT_0221 +
       conj(IT_0231)*(24*IT_0088 + IT_0134 + IT_0232) + IT_0231*(24*conj(IT_0088
      ) + IT_0233 + IT_0234) + conj(IT_0088)*IT_0235 + IT_0088*IT_0236;
    const ccomplex_t IT_0238 = 12*IT_0029;
    const ccomplex_t IT_0239 = 12*IT_0088;
    const ccomplex_t IT_0240 = m_s*s_14*m_N_2;
    const ccomplex_t IT_0241 = 6*IT_0217;
    const ccomplex_t IT_0242 = 6*conj(IT_0217);
    const ccomplex_t IT_0243 = 6*IT_0184;
    const ccomplex_t IT_0244 = 6*conj(IT_0184);
    const ccomplex_t IT_0245 = 12*conj(IT_0141);
    const ccomplex_t IT_0246 = 12*IT_0172;
    const ccomplex_t IT_0247 = IT_0142*conj(IT_0156) + conj(IT_0160)*IT_0200 +
       IT_0160*IT_0202 + conj(IT_0211)*IT_0218 + IT_0211*IT_0219 + conj(IT_0088)
      *IT_0220 + IT_0088*IT_0221 + (24*conj(IT_0029) + 12*conj(IT_0072) + 12
      *conj(IT_0172))*IT_0231 + conj(IT_0029)*IT_0235 + IT_0029*IT_0236 + conj
      (IT_0131)*IT_0241 + IT_0131*IT_0242 + (conj(IT_0072) + conj(IT_0172))
      *IT_0243 + (IT_0072 + IT_0172)*IT_0244 + IT_0156*IT_0245 + conj(IT_0231)*
      (24*IT_0029 + IT_0133 + IT_0246);
    const ccomplex_t IT_0248 = 12*IT_0211;
    const ccomplex_t IT_0249 = 24*IT_0156;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 12*IT_0217;
    const ccomplex_t IT_0252 = -IT_0001;
    const ccomplex_t IT_0253 = IT_0000 + IT_0252;
    const ccomplex_t IT_0254 = IT_0002 + IT_0253;
    const ccomplex_t IT_0255 = (-6)*IT_0072;
    const ccomplex_t IT_0256 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0257 = (-6)*IT_0122;
    const ccomplex_t IT_0258 = (-6)*conj(IT_0122);
    const ccomplex_t IT_0259 = s_34*m_C_2*m_N_2;
    const ccomplex_t IT_0260 = 12*conj(IT_0211);
    const ccomplex_t IT_0261 = (-12)*IT_0122;
    const ccomplex_t IT_0262 = (-12)*conj(IT_0122);
    const ccomplex_t IT_0263 = (-12)*IT_0141;
    const ccomplex_t IT_0264 = (-12)*conj(IT_0141);
    const ccomplex_t IT_0265 = (-6)*IT_0131;
    const ccomplex_t IT_0266 = (-6)*conj(IT_0131);
    const ccomplex_t IT_0267 = (-6)*IT_0141;
    const ccomplex_t IT_0268 = (-6)*conj(IT_0141);
    const ccomplex_t IT_0269 = conj(IT_0172)*IT_0218 + IT_0172*IT_0219 + conj
      (IT_0156)*IT_0235 + IT_0156*IT_0236 + conj(IT_0160)*IT_0241 + IT_0160
      *IT_0242 + conj(IT_0211)*IT_0243 + IT_0211*IT_0244 + conj(IT_0231)*IT_0250
       + IT_0231*(24*conj(IT_0156) + IT_0260) + conj(IT_0131)*(IT_0199 + IT_0261
      ) + IT_0131*(IT_0201 + IT_0262) + conj(IT_0072)*IT_0263 + IT_0072*IT_0264 
      + conj(IT_0088)*IT_0265 + IT_0088*IT_0266 + conj(IT_0029)*IT_0267 +
       IT_0029*IT_0268;
    const ccomplex_t IT_0270 = m_s*s_24*m_C_2;
    const ccomplex_t IT_0271 = IT_0232 + IT_0257;
    const ccomplex_t IT_0272 = IT_0233 + IT_0258;
    const ccomplex_t IT_0273 = 6*IT_0211;
    const ccomplex_t IT_0274 = 6*conj(IT_0211);
    const ccomplex_t IT_0275 = 12*conj(IT_0029);
    const ccomplex_t IT_0276 = conj(IT_0156)*((-12)*IT_0072 + 24*IT_0172 +
       IT_0238) + conj(IT_0088)*IT_0241 + IT_0088*IT_0242 + conj(IT_0172)
      *IT_0248 + conj(IT_0211)*IT_0255 + IT_0211*IT_0256 + IT_0172*IT_0260 +
       conj(IT_0231)*IT_0263 + IT_0231*IT_0264 + conj(IT_0160)*IT_0265 + IT_0160
      *IT_0266 + conj(IT_0184)*IT_0267 + IT_0184*IT_0268 + conj(IT_0217)*IT_0271
       + IT_0217*IT_0272 + conj(IT_0029)*IT_0273 + IT_0029*IT_0274 + IT_0156*((
      -12)*conj(IT_0072) + 24*conj(IT_0172) + IT_0275);
    const ccomplex_t IT_0277 = m_c*m_s*s_12;
    const ccomplex_t IT_0278 = 12*conj(IT_0217);
    const ccomplex_t IT_0279 = (-12)*conj(IT_0131);
    const ccomplex_t IT_0280 = (-12)*conj(IT_0072)*IT_0131 + conj(IT_0172)
      *IT_0212 + IT_0172*IT_0213 + conj(IT_0198)*IT_0218 + IT_0198*IT_0219 +
       conj(IT_0156)*IT_0220 + IT_0156*IT_0221 + conj(IT_0184)*IT_0241 + IT_0184
      *IT_0242 + conj(IT_0231)*IT_0251 + conj(IT_0141)*IT_0261 + IT_0141*IT_0262
       + conj(IT_0029)*IT_0265 + IT_0029*IT_0266 + conj(IT_0088)*IT_0267 +
       IT_0088*IT_0268 + conj(IT_0160)*IT_0273 + IT_0160*IT_0274 + IT_0231
      *IT_0278 + IT_0072*IT_0279;
    const ccomplex_t IT_0281 = IT_0003 + IT_0253;
    const ccomplex_t IT_0282 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0283 = (-0.166666666666667)*IT_0088;
    const ccomplex_t IT_0284 = m_c*m_s*m_C_2*m_N_2;
    const ccomplex_t IT_0285 = (-0.166666666666667)*IT_0211;
    const ccomplex_t IT_0286 = (-0.166666666666667)*conj(IT_0211);
    const ccomplex_t IT_0287 = (-0.166666666666667)*IT_0029;
    const ccomplex_t IT_0288 = (-0.166666666666667)*conj(IT_0029);
    const ccomplex_t IT_0289 = (-0.166666666666667)*IT_0172;
    const ccomplex_t IT_0290 = (-0.166666666666667)*conj(IT_0172);
    const ccomplex_t IT_0291 = 0.0416666666666667*conj(IT_0131)*IT_0142 + (
      -0.25)*conj(IT_0072)*((-0.166666666666667)*IT_0134 + IT_0198) + conj
      (IT_0156)*IT_0217 + IT_0156*conj(IT_0217) + conj(IT_0160)*(IT_0231 +
       0.0416666666666667*IT_0235) + IT_0160*(conj(IT_0231) + 0.0416666666666667
      *IT_0236) + 0.25*conj(IT_0088)*(IT_0173 + 0.166666666666667*IT_0238) +
       0.0416666666666667*IT_0131*IT_0245 + 0.25*IT_0088*(IT_0175 +
       0.166666666666667*IT_0275) + 0.0416666666666667*IT_0072*(IT_0234 +
       IT_0282) + (-0.25)*IT_0278*IT_0285 + (-0.25)*IT_0251*IT_0286 + (-0.25)
      *IT_0202*IT_0287 + (-0.25)*IT_0200*IT_0288 + (-0.25)*IT_0272*IT_0289 + (
      -0.25)*IT_0271*IT_0290;
    const ccomplex_t IT_0292 = m_c*s_13*m_N_2;
    const ccomplex_t IT_0293 = conj(IT_0131)*IT_0184 + conj(IT_0072)*IT_0217 +
       IT_0072*conj(IT_0217) + IT_0131*(conj(IT_0184) + 2*conj(IT_0231)) + (-4)
      *conj(IT_0156)*(IT_0198 + 0.0416666666666667*IT_0239 + 0.0416666666666667
      *IT_0261) + (-4)*IT_0156*(0.5*conj(IT_0088) + conj(IT_0198) +
       0.0416666666666667*IT_0262) + (-0.166666666666667)*conj(IT_0160)*IT_0267 
      + (-0.166666666666667)*IT_0160*IT_0268 + (-0.166666666666667)*conj(IT_0088
      )*IT_0273 + (-0.166666666666667)*IT_0231*IT_0279 + IT_0274*IT_0283 +
       IT_0272*IT_0285 + IT_0271*IT_0286 + IT_0242*IT_0287 + IT_0241*IT_0288 +
       IT_0278*IT_0289 + IT_0251*IT_0290;
    const ccomplex_t IT_0294 = IT_0004*(conj(IT_0029)*IT_0073 + IT_0029
      *IT_0074 + conj(IT_0088)*IT_0123 + IT_0088*IT_0124) + IT_0000*(conj
      (IT_0131)*IT_0132 + conj(IT_0072)*IT_0133 + conj(IT_0122)*IT_0134 + conj
      (IT_0141)*IT_0142) + IT_0143*IT_0237 + IT_0001*(conj(IT_0160)*IT_0220 + 
      (24*IT_0184 + 48*IT_0231)*conj(IT_0231) + conj(IT_0184)*(24*IT_0231 +
       IT_0235) + conj(IT_0029)*IT_0238 + conj(IT_0088)*IT_0239) + IT_0240
      *IT_0247 + IT_0002*(conj(IT_0156)*(48*IT_0156 + 24*IT_0211) + conj(IT_0198
      )*IT_0232 + conj(IT_0172)*IT_0246 + conj(IT_0211)*IT_0250 + conj(IT_0217)
      *IT_0251) + IT_0254*(conj(IT_0172)*IT_0255 + IT_0172*IT_0256 + conj
      (IT_0198)*IT_0257 + IT_0198*IT_0258) + IT_0259*IT_0269 + IT_0270*IT_0276 +
       IT_0277*IT_0280 + (-6)*IT_0281*(conj(IT_0029)*IT_0172 + IT_0029*conj
      (IT_0172) + conj(IT_0088)*IT_0198 + IT_0282*IT_0283) + 24*IT_0284*IT_0291 
      + (-6)*IT_0292*IT_0293;
    return create_ccomplex_return(IT_0294);
}

