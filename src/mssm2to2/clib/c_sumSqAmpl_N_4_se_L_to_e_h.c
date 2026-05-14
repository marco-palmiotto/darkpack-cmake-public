#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_se_L_to_e_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_se_L_to_e_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0005 + IT_0009 + 
      -IT_0014 + -IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = m_N_4*IT_0019;
    const ccomplex_t IT_0021 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0022 = IT_0003*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0028 = IT_0015*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0022 + IT_0024 + 
      -IT_0026 + -IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = m_N_4*IT_0030;
    const ccomplex_t IT_0032 = IT_0020 + IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_e, 2);
    const ccomplex_t IT_0039 = pow(m_se_L, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0032*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0000
      *IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0002;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = pow(m_N_4, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0038 + IT_0049 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0059 = IT_0003*IT_0058;
    const ccomplex_t IT_0060 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0061 = IT_0007*IT_0060;
    const ccomplex_t IT_0062 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = IT_0059 + IT_0061 + IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0068 = IT_0012*IT_0067;
    const ccomplex_t IT_0069 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0070 = IT_0015*IT_0069;
    const ccomplex_t IT_0071 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0072 = IT_0012*IT_0071;
    const ccomplex_t IT_0073 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0074 = IT_0015*IT_0073;
    const ccomplex_t IT_0075 = -IT_0068 + -IT_0070 + -IT_0072 + -IT_0074;
    const ccomplex_t IT_0076 = IT_0066 + IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = m_N_2*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0081 = IT_0003*IT_0080;
    const ccomplex_t IT_0082 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0083 = IT_0007*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0085 = IT_0003*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0007*IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0083 + IT_0085 + IT_0087;
    const ccomplex_t IT_0089 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0090 = IT_0012*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0092 = IT_0015*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0094 = IT_0012*IT_0093;
    const ccomplex_t IT_0095 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0096 = IT_0015*IT_0095;
    const ccomplex_t IT_0097 = -IT_0090 + -IT_0092 + -IT_0094 + -IT_0096;
    const ccomplex_t IT_0098 = IT_0088 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = m_N_4*IT_0100;
    const ccomplex_t IT_0102 = IT_0079 + IT_0101;
    const ccomplex_t IT_0103 = IT_0057*IT_0102;
    const ccomplex_t IT_0104 = sin(beta);
    const ccomplex_t IT_0105 = cpow(IT_0001, -2);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*e_em*(IT_0006*IT_0011
      *IT_0033*IT_0035*IT_0038 + (-0.5)*m_W*(IT_0006*IT_0034 + IT_0000*IT_0104)*
      (IT_0011 + -IT_0010*IT_0105));
    const ccomplex_t IT_0107 = IT_0037*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_13 + IT_0038 + IT_0049 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_e*e_em*IT_0006*IT_0011
      *IT_0033*IT_0035;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = IT_0037*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_12 + IT_0038 + -IT_0039 + -IT_0049 
      + -reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = m_e*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0124 = IT_0003*IT_0123;
    const ccomplex_t IT_0125 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0126 = IT_0007*IT_0125;
    const ccomplex_t IT_0127 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0128 = IT_0003*IT_0127;
    const ccomplex_t IT_0129 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0130 = IT_0007*IT_0129;
    const ccomplex_t IT_0131 = IT_0124 + IT_0126 + IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0133 = IT_0012*IT_0132;
    const ccomplex_t IT_0134 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0135 = IT_0015*IT_0134;
    const ccomplex_t IT_0136 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0137 = IT_0012*IT_0136;
    const ccomplex_t IT_0138 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0139 = IT_0015*IT_0138;
    const ccomplex_t IT_0140 = -IT_0133 + -IT_0135 + -IT_0137 + -IT_0139;
    const ccomplex_t IT_0141 = IT_0131 + IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = m_N_1*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0146 = IT_0003*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0148 = IT_0007*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0150 = IT_0003*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0152 = IT_0007*IT_0151;
    const ccomplex_t IT_0153 = IT_0146 + IT_0148 + IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0155 = IT_0012*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0157 = IT_0015*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0159 = IT_0012*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0161 = IT_0015*IT_0160;
    const ccomplex_t IT_0162 = -IT_0155 + -IT_0157 + -IT_0159 + -IT_0161;
    const ccomplex_t IT_0163 = IT_0153 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = m_N_4*IT_0165;
    const ccomplex_t IT_0167 = IT_0144 + IT_0166;
    const ccomplex_t IT_0168 = IT_0122*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0175 = IT_0003*IT_0174;
    const ccomplex_t IT_0176 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0177 = IT_0007*IT_0176;
    const ccomplex_t IT_0178 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0179 = IT_0003*IT_0178;
    const ccomplex_t IT_0180 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0181 = IT_0007*IT_0180;
    const ccomplex_t IT_0182 = IT_0175 + IT_0177 + IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0184 = IT_0012*IT_0183;
    const ccomplex_t IT_0185 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0186 = IT_0015*IT_0185;
    const ccomplex_t IT_0187 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0188 = IT_0012*IT_0187;
    const ccomplex_t IT_0189 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0190 = IT_0015*IT_0189;
    const ccomplex_t IT_0191 = -IT_0184 + -IT_0186 + -IT_0188 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0182 + IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = m_N_3*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0197 = IT_0003*IT_0196;
    const ccomplex_t IT_0198 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0199 = IT_0007*IT_0198;
    const ccomplex_t IT_0200 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0201 = IT_0003*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0203 = IT_0007*IT_0202;
    const ccomplex_t IT_0204 = IT_0197 + IT_0199 + IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0206 = IT_0012*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0208 = IT_0015*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0210 = IT_0012*IT_0209;
    const ccomplex_t IT_0211 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0212 = IT_0015*IT_0211;
    const ccomplex_t IT_0213 = -IT_0206 + -IT_0208 + -IT_0210 + -IT_0212;
    const ccomplex_t IT_0214 = IT_0204 + IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = m_N_4*IT_0216;
    const ccomplex_t IT_0218 = IT_0195 + IT_0217;
    const ccomplex_t IT_0219 = IT_0173*IT_0218;
    const ccomplex_t IT_0220 = IT_0043 + IT_0052 + IT_0103 + IT_0110 + (-2)
      *IT_0117 + IT_0168 + IT_0219;
    const ccomplex_t IT_0221 = conj(N_B4)*e_em;
    const ccomplex_t IT_0222 = IT_0002*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W4)*e_em;
    const ccomplex_t IT_0225 = IT_0011*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + IT_0226);
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0112*IT_0228;
    const ccomplex_t IT_0230 = IT_0114*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = m_e*IT_0231;
    const ccomplex_t IT_0233 = conj(N_B1)*e_em;
    const ccomplex_t IT_0234 = IT_0002*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W1)*e_em;
    const ccomplex_t IT_0237 = IT_0011*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0235 + IT_0238);
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = IT_0120*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = m_N_4*IT_0143;
    const ccomplex_t IT_0244 = m_N_1*IT_0165;
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0242*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B2)*e_em;
    const ccomplex_t IT_0248 = IT_0002*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = conj(N_W2)*e_em;
    const ccomplex_t IT_0251 = IT_0011*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0249 + IT_0252);
    const ccomplex_t IT_0254 = (-0.5)*IT_0253;
    const ccomplex_t IT_0255 = IT_0055*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = m_N_4*IT_0078;
    const ccomplex_t IT_0258 = m_N_2*IT_0100;
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = IT_0256*IT_0259;
    const ccomplex_t IT_0261 = IT_0040*IT_0228;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = IT_0032*IT_0262;
    const ccomplex_t IT_0264 = conj(N_B3)*e_em;
    const ccomplex_t IT_0265 = IT_0002*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = conj(N_W3)*e_em;
    const ccomplex_t IT_0268 = IT_0011*IT_0267;
    const ccomplex_t IT_0269 = 1.4142135623731*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*(IT_0266 + IT_0269);
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = IT_0171*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = m_N_4*IT_0194;
    const ccomplex_t IT_0275 = m_N_3*IT_0216;
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = IT_0273*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0011*IT_0033*IT_0035;
    const ccomplex_t IT_0279 = 0.5*IT_0278;
    const ccomplex_t IT_0280 = IT_0045*IT_0279;
    const ccomplex_t IT_0281 = IT_0050*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0106*IT_0228;
    const ccomplex_t IT_0284 = IT_0108*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0232 + -IT_0246 + -IT_0260 + -IT_0263 + 
      -IT_0277 + -IT_0282 + -IT_0285;
    const ccomplex_t IT_0287 = 2*s_13;
    const ccomplex_t IT_0288 = m_e*m_N_4;
    const ccomplex_t IT_0289 = 2*IT_0220;
    const ccomplex_t IT_0290 = IT_0019*IT_0262;
    const ccomplex_t IT_0291 = IT_0194*IT_0273;
    const ccomplex_t IT_0292 = IT_0143*IT_0242;
    const ccomplex_t IT_0293 = IT_0078*IT_0256;
    const ccomplex_t IT_0294 = -IT_0116 + IT_0290 + IT_0291 + IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = m_e*s_14;
    const ccomplex_t IT_0296 = IT_0122*IT_0165;
    const ccomplex_t IT_0297 = IT_0057*IT_0100;
    const ccomplex_t IT_0298 = IT_0173*IT_0216;
    const ccomplex_t IT_0299 = IT_0030*IT_0042;
    const ccomplex_t IT_0300 = IT_0231 + -IT_0296 + -IT_0297 + -IT_0298 + 
      -IT_0299;
    const ccomplex_t IT_0301 = s_34*m_N_4;
    const ccomplex_t IT_0302 = 2*conj(IT_0220);
    const ccomplex_t IT_0303 = 2*conj(IT_0286);
    const ccomplex_t IT_0304 = pow(m_h, 2);
    const ccomplex_t IT_0305 = m_e*m_N_4*IT_0304;
    const ccomplex_t IT_0306 = s_13*IT_0304;
    const ccomplex_t IT_0307 = s_14*s_34;
    const ccomplex_t IT_0308 = (-2)*IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0286;
    const ccomplex_t IT_0311 = 0.5*IT_0301;
    const ccomplex_t IT_0312 = (IT_0220*conj(IT_0220) + IT_0286*conj(IT_0286))
      *IT_0287 + conj(IT_0286)*IT_0288*IT_0289 + (IT_0286*IT_0288 + IT_0294
      *IT_0295 + IT_0300*IT_0301)*IT_0302 + (IT_0295*IT_0300 + IT_0294*IT_0301)
      *IT_0303 + 2*conj(IT_0300)*(IT_0294*IT_0305 + -IT_0300*IT_0309 + 0.5
      *IT_0295*IT_0310 + IT_0289*IT_0311) + 2*conj(IT_0294)*(0.5*IT_0289*IT_0295
       + IT_0300*IT_0305 + -IT_0294*IT_0309 + IT_0310*IT_0311);
    return create_ccomplex_return(IT_0312);
}

