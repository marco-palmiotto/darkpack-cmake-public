#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_se_R_to_e_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_se_R_to_e_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = pow(m_e, 2);
    const ccomplex_t IT_0008 = pow(m_N_4, 2);
    const ccomplex_t IT_0009 = pow(m_se_R, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0009 
      + -reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*IT_0012*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = IT_0012*IT_0015;
    const ccomplex_t IT_0022 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0015, -1);
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0027;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + IT_0026 + 
      -IT_0028 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0030;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0031 + IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0014*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0034*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0041 = IT_0021*IT_0040;
    const ccomplex_t IT_0042 = IT_0025*IT_0040;
    const ccomplex_t IT_0043 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0044 = IT_0021*IT_0043;
    const ccomplex_t IT_0045 = IT_0025*IT_0043;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0041 + IT_0042 + 
      -IT_0044 + -IT_0045);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0012;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0047*IT_0052;
    const ccomplex_t IT_0054 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0055 = IT_0021*IT_0054;
    const ccomplex_t IT_0056 = IT_0025*IT_0054;
    const ccomplex_t IT_0057 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0058 = IT_0025*IT_0057;
    const ccomplex_t IT_0059 = IT_0021*IT_0057;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0055 + IT_0056 + 
      -IT_0058 + -IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0012;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0064*IT_0069;
    const ccomplex_t IT_0071 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0072 = IT_0021*IT_0071;
    const ccomplex_t IT_0073 = IT_0025*IT_0071;
    const ccomplex_t IT_0074 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0075 = IT_0021*IT_0074;
    const ccomplex_t IT_0076 = IT_0025*IT_0074;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0072 + IT_0073 + 
      -IT_0075 + -IT_0076);
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0012;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0078*IT_0083;
    const ccomplex_t IT_0085 = -IT_0053 + -IT_0070 + -IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0039 + IT_0086;
    const ccomplex_t IT_0088 = conj(IT_0039) + -conj(IT_0085);
    const ccomplex_t IT_0089 = pow(m_W, -1);
    const ccomplex_t IT_0090 = cos(beta);
    const ccomplex_t IT_0091 = cpow(IT_0090, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0024*IT_0089*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = e_em*IT_0021;
    const ccomplex_t IT_0095 = e_em*IT_0025;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0093*IT_0097;
    const ccomplex_t IT_0099 = IT_0010*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0035*IT_0093;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0034*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0106 = IT_0021*IT_0105;
    const ccomplex_t IT_0107 = IT_0025*IT_0105;
    const ccomplex_t IT_0108 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0109 = IT_0021*IT_0108;
    const ccomplex_t IT_0110 = IT_0025*IT_0108;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0106 + IT_0107 + 
      -IT_0109 + -IT_0110);
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0024*IT_0089*IT_0091;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0081*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0112*IT_0116;
    const ccomplex_t IT_0118 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0119 = IT_0021*IT_0118;
    const ccomplex_t IT_0120 = IT_0025*IT_0118;
    const ccomplex_t IT_0121 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0122 = IT_0021*IT_0121;
    const ccomplex_t IT_0123 = IT_0025*IT_0121;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0119 + IT_0120 + 
      -IT_0122 + -IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0024*IT_0089*IT_0091;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0050*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0125*IT_0129;
    const ccomplex_t IT_0131 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0132 = IT_0021*IT_0131;
    const ccomplex_t IT_0133 = IT_0025*IT_0131;
    const ccomplex_t IT_0134 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0135 = IT_0025*IT_0134;
    const ccomplex_t IT_0136 = IT_0021*IT_0134;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0132 + IT_0133 + 
      -IT_0135 + -IT_0136);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0024*IT_0089*IT_0091;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0067*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0139*IT_0143;
    const ccomplex_t IT_0145 = IT_0117 + IT_0130 + IT_0144;
    const ccomplex_t IT_0146 = IT_0104 + -IT_0145;
    const ccomplex_t IT_0147 = conj(IT_0104) + -conj(IT_0145);
    const ccomplex_t IT_0148 = pow(m_Z, 2);
    const ccomplex_t IT_0149 = s_13*IT_0148;
    const ccomplex_t IT_0150 = (-0.5)*IT_0004;
    const ccomplex_t IT_0151 = (-2)*IT_0000;
    const ccomplex_t IT_0152 = IT_0149 + IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = m_e*m_N_4;
    const ccomplex_t IT_0154 = pow(s_14, 2);
    const ccomplex_t IT_0155 = IT_0001*IT_0154;
    const ccomplex_t IT_0156 = -IT_0008;
    const ccomplex_t IT_0157 = IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = IT_0153*IT_0157;
    const ccomplex_t IT_0159 = (-2)*IT_0016;
    const ccomplex_t IT_0160 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0161 = IT_0093*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0159*IT_0162;
    const ccomplex_t IT_0164 = IT_0031*IT_0102;
    const ccomplex_t IT_0165 = IT_0032*IT_0102;
    const ccomplex_t IT_0166 = -IT_0163 + IT_0164 + -IT_0165;
    const ccomplex_t IT_0167 = IT_0014*IT_0160;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0159*IT_0168;
    const ccomplex_t IT_0170 = IT_0031*IT_0037;
    const ccomplex_t IT_0171 = IT_0032*IT_0037;
    const ccomplex_t IT_0172 = IT_0169 + IT_0170 + -IT_0171;
    const ccomplex_t IT_0173 = 2*IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0053 + 2*IT_0070 + 2*IT_0084;
    const ccomplex_t IT_0175 = (-2)*IT_0174;
    const ccomplex_t IT_0176 = IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = 2*conj(IT_0172);
    const ccomplex_t IT_0178 = (-2)*conj(IT_0174);
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = (-2)*IT_0117 + (-2)*IT_0130 + (-2)*IT_0144;
    const ccomplex_t IT_0181 = 2*IT_0174;
    const ccomplex_t IT_0182 = (-2)*IT_0172;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = 2*conj(IT_0174);
    const ccomplex_t IT_0185 = (-2)*conj(IT_0172);
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = m_e*IT_0157;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0189 = IT_0014*IT_0097;
    const ccomplex_t IT_0190 = IT_0010*IT_0188*IT_0189;
    const ccomplex_t IT_0191 = m_N_4*IT_0031;
    const ccomplex_t IT_0192 = IT_0102*IT_0191;
    const ccomplex_t IT_0193 = m_N_4*IT_0032;
    const ccomplex_t IT_0194 = -IT_0191;
    const ccomplex_t IT_0195 = IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = IT_0102*IT_0195;
    const ccomplex_t IT_0197 = m_e*IT_0019;
    const ccomplex_t IT_0198 = -IT_0190 + (-0.5)*IT_0192 + 0.5*IT_0196 +
       IT_0197;
    const ccomplex_t IT_0199 = IT_0102*IT_0193;
    const ccomplex_t IT_0200 = m_N_1*IT_0047;
    const ccomplex_t IT_0201 = m_N_4*IT_0125;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = IT_0129*IT_0203;
    const ccomplex_t IT_0205 = m_N_2*IT_0064;
    const ccomplex_t IT_0206 = m_N_4*IT_0139;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = IT_0143*IT_0208;
    const ccomplex_t IT_0210 = m_N_3*IT_0078;
    const ccomplex_t IT_0211 = m_N_4*IT_0112;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = IT_0116*IT_0213;
    const ccomplex_t IT_0215 = 0.5*IT_0199 + IT_0204 + IT_0209 + IT_0214;
    const ccomplex_t IT_0216 = conj(IT_0198) + conj(IT_0215);
    const ccomplex_t IT_0217 = IT_0198 + IT_0215;
    const ccomplex_t IT_0218 = IT_0016*IT_0093;
    const ccomplex_t IT_0219 = IT_0010*IT_0188*IT_0218;
    const ccomplex_t IT_0220 = m_e*IT_0100;
    const ccomplex_t IT_0221 = IT_0037*IT_0193;
    const ccomplex_t IT_0222 = IT_0037*IT_0195;
    const ccomplex_t IT_0223 = IT_0219 + -IT_0220 + 0.5*IT_0221 + 0.5*IT_0222;
    const ccomplex_t IT_0224 = IT_0037*IT_0191;
    const ccomplex_t IT_0225 = m_N_1*IT_0125;
    const ccomplex_t IT_0226 = m_N_4*IT_0047;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = IT_0052*IT_0228;
    const ccomplex_t IT_0230 = m_N_2*IT_0139;
    const ccomplex_t IT_0231 = m_N_4*IT_0064;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = IT_0069*IT_0233;
    const ccomplex_t IT_0235 = m_N_3*IT_0112;
    const ccomplex_t IT_0236 = m_N_4*IT_0078;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0083*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0224 + -IT_0229 + -IT_0234 + -IT_0239;
    const ccomplex_t IT_0241 = conj(IT_0223) + conj(IT_0240);
    const ccomplex_t IT_0242 = 2*IT_0166;
    const ccomplex_t IT_0243 = (-2)*IT_0180;
    const ccomplex_t IT_0244 = IT_0242 + IT_0243;
    const ccomplex_t IT_0245 = IT_0223 + IT_0240;
    const ccomplex_t IT_0246 = 2*conj(IT_0166);
    const ccomplex_t IT_0247 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0248 = IT_0246 + IT_0247;
    const ccomplex_t IT_0249 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = s_13 + IT_0250;
    const ccomplex_t IT_0252 = s_13*IT_0251;
    const ccomplex_t IT_0253 = (-2)*IT_0019 + -IT_0169;
    const ccomplex_t IT_0254 = 2*IT_0100 + IT_0163;
    const ccomplex_t IT_0255 = 2*IT_0180;
    const ccomplex_t IT_0256 = (-2)*IT_0166;
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = 2*conj(IT_0180);
    const ccomplex_t IT_0259 = (-2)*conj(IT_0166);
    const ccomplex_t IT_0260 = IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = s_13*s_14;
    const ccomplex_t IT_0262 = s_34*IT_0008;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = IT_0261 + IT_0263;
    const ccomplex_t IT_0265 = conj(IT_0020) + conj(IT_0085);
    const ccomplex_t IT_0266 = IT_0020 + IT_0085;
    const ccomplex_t IT_0267 = conj(IT_0100) + conj(IT_0145);
    const ccomplex_t IT_0268 = IT_0100 + IT_0145;
    const ccomplex_t IT_0269 = s_34*IT_0148;
    const ccomplex_t IT_0270 = IT_0001*IT_0269;
    const ccomplex_t IT_0271 = 0.5*IT_0270;
    const ccomplex_t IT_0272 = s_34 + IT_0271;
    const ccomplex_t IT_0273 = m_N_4*IT_0272;
    const ccomplex_t IT_0274 = 2*IT_0273;
    const ccomplex_t IT_0275 = 2*IT_0198;
    const ccomplex_t IT_0276 = 2*IT_0215;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = 2*conj(IT_0198);
    const ccomplex_t IT_0279 = 2*conj(IT_0215);
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = 2*IT_0223;
    const ccomplex_t IT_0282 = 2*IT_0240;
    const ccomplex_t IT_0283 = IT_0281 + IT_0282;
    const ccomplex_t IT_0284 = 2*conj(IT_0223);
    const ccomplex_t IT_0285 = 2*conj(IT_0240);
    const ccomplex_t IT_0286 = IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = s_14*IT_0148;
    const ccomplex_t IT_0288 = IT_0001*IT_0287;
    const ccomplex_t IT_0289 = (-0.25)*IT_0288;
    const ccomplex_t IT_0290 = s_14 + IT_0289;
    const ccomplex_t IT_0291 = m_e*IT_0290;
    const ccomplex_t IT_0292 = (-4)*IT_0291;
    const ccomplex_t IT_0293 = (-0.25)*IT_0270;
    const ccomplex_t IT_0294 = s_34 + IT_0293;
    const ccomplex_t IT_0295 = m_N_4*IT_0294;
    const ccomplex_t IT_0296 = (-4)*IT_0295;
    const ccomplex_t IT_0297 = (-2)*IT_0198;
    const ccomplex_t IT_0298 = (-2)*IT_0215;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0301 = (-2)*conj(IT_0215);
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = (-2)*IT_0223;
    const ccomplex_t IT_0304 = (-2)*IT_0240;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = (-2)*conj(IT_0223);
    const ccomplex_t IT_0307 = (-2)*conj(IT_0240);
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = 0.5*IT_0288;
    const ccomplex_t IT_0310 = s_14 + IT_0309;
    const ccomplex_t IT_0311 = m_e*IT_0310;
    const ccomplex_t IT_0312 = 2*IT_0311;
    const ccomplex_t IT_0313 = IT_0148*IT_0153;
    const ccomplex_t IT_0314 = 3*IT_0313;
    const ccomplex_t IT_0315 = 2*IT_0039;
    const ccomplex_t IT_0316 = (-0.5)*IT_0315;
    const ccomplex_t IT_0317 = IT_0085 + IT_0316;
    const ccomplex_t IT_0318 = 2*conj(IT_0039);
    const ccomplex_t IT_0319 = (-0.5)*IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0085) + IT_0319;
    const ccomplex_t IT_0321 = 2*IT_0020;
    const ccomplex_t IT_0322 = 2*conj(IT_0020);
    const ccomplex_t IT_0323 = (-3)*IT_0313;
    const ccomplex_t IT_0324 = IT_0149 + IT_0151;
    const ccomplex_t IT_0325 = IT_0001*IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = s_13 + IT_0326;
    const ccomplex_t IT_0328 = 4*IT_0223;
    const ccomplex_t IT_0329 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0330 = pow(s_34, 2);
    const ccomplex_t IT_0331 = IT_0001*IT_0330;
    const ccomplex_t IT_0332 = -IT_0331;
    const ccomplex_t IT_0333 = IT_0007 + IT_0332;
    const ccomplex_t IT_0334 = s_13*IT_0333;
    const ccomplex_t IT_0335 = (-2)*IT_0253;
    const ccomplex_t IT_0336 = IT_0153*IT_0333;
    const ccomplex_t IT_0337 = (-2)*conj(IT_0253);
    const ccomplex_t IT_0338 = IT_0153*IT_0251;
    const ccomplex_t IT_0339 = 2*IT_0253;
    const ccomplex_t IT_0340 = 2*conj(IT_0253);
    const ccomplex_t IT_0341 = s_13*IT_0157;
    const ccomplex_t IT_0342 = (-2)*IT_0254;
    const ccomplex_t IT_0343 = m_e*IT_0251;
    const ccomplex_t IT_0344 = (-2)*conj(IT_0254);
    const ccomplex_t IT_0345 = s_14*IT_0007;
    const ccomplex_t IT_0346 = s_13*s_34;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = IT_0345 + IT_0347;
    const ccomplex_t IT_0349 = conj(IT_0104)*IT_0254 + IT_0104*conj(IT_0254) +
       0.5*IT_0265*IT_0335 + 0.5*IT_0266*IT_0337 + 0.5*conj(IT_0039)*IT_0339 +
       0.5*IT_0039*IT_0340 + 0.5*IT_0267*IT_0342 + 0.5*IT_0268*IT_0344;
    const ccomplex_t IT_0350 = m_N_4*IT_0251;
    const ccomplex_t IT_0351 = m_N_4*IT_0333;
    const ccomplex_t IT_0352 = (-8)*IT_0006*(conj(IT_0020)*IT_0087 + IT_0020
      *IT_0088 + conj(IT_0100)*IT_0146 + IT_0100*IT_0147) + (-4)*(IT_0020*conj
      (IT_0020) + IT_0087*IT_0088 + IT_0100*conj(IT_0100) + IT_0146*IT_0147)
      *IT_0152 + IT_0158*(conj(IT_0166)*IT_0176 + IT_0166*IT_0179 + conj(IT_0180
      )*IT_0183 + IT_0180*IT_0186) + IT_0187*(IT_0176*IT_0216 + IT_0179*IT_0217 
      + IT_0241*IT_0244 + IT_0245*IT_0248) + IT_0252*(IT_0186*IT_0253 + IT_0183
      *conj(IT_0253) + conj(IT_0254)*IT_0257 + IT_0254*IT_0260) + IT_0264*(conj
      (IT_0039)*IT_0176 + IT_0039*IT_0179 + conj(IT_0104)*IT_0244 + IT_0104
      *IT_0248 + IT_0183*IT_0265 + IT_0186*IT_0266 + IT_0257*IT_0267 + IT_0260
      *IT_0268) + IT_0274*(conj(IT_0100)*IT_0277 + IT_0100*IT_0280 + conj
      (IT_0020)*IT_0283 + IT_0020*IT_0286) + (conj(IT_0020)*IT_0277 + IT_0020
      *IT_0280 + conj(IT_0100)*IT_0283 + IT_0100*IT_0286)*IT_0292 + IT_0296*
      (conj(IT_0104)*IT_0277 + IT_0104*IT_0280 + conj(IT_0039)*IT_0283 + IT_0039
      *IT_0286 + conj(IT_0145)*IT_0299 + IT_0145*IT_0302 + conj(IT_0085)*IT_0305
       + IT_0085*IT_0308) + (conj(IT_0039)*IT_0277 + IT_0039*IT_0280 + conj
      (IT_0104)*IT_0283 + IT_0104*IT_0286 + conj(IT_0085)*IT_0299 + IT_0085
      *IT_0302 + conj(IT_0145)*IT_0305 + IT_0145*IT_0308)*IT_0312 + (-2)*IT_0314
      *(conj(IT_0020)*IT_0145 + IT_0020*conj(IT_0145) + conj(IT_0100)*IT_0317 +
       IT_0100*IT_0320 + (-0.5)*conj(IT_0104)*IT_0321 + (-0.5)*IT_0104*IT_0322) 
      + (-2)*(IT_0088*IT_0145 + IT_0087*conj(IT_0145) + conj(IT_0104)*IT_0317 +
       IT_0104*IT_0320 + (-0.5)*conj(IT_0100)*IT_0321 + (-0.5)*IT_0100*IT_0322)
      *IT_0323 + 4*IT_0327*(IT_0216*IT_0217 + IT_0241*(IT_0240 + 0.25*IT_0328)) 
      + (-6)*IT_0153*(IT_0217*IT_0241 + IT_0245*(conj(IT_0215) + (
      -0.166666666666667)*IT_0329)) + conj(IT_0253)*IT_0334*IT_0335 + IT_0254
      *IT_0336*IT_0337 + IT_0338*(IT_0186*IT_0254 + IT_0183*conj(IT_0254) + conj
      (IT_0166)*IT_0335 + IT_0166*IT_0337 + conj(IT_0180)*IT_0339 + IT_0180
      *IT_0340) + (conj(IT_0172)*IT_0173 + conj(IT_0174)*IT_0183 + IT_0174
      *IT_0185 + conj(IT_0166)*IT_0242 + conj(IT_0180)*IT_0257 + IT_0180*IT_0259
      )*IT_0341 + conj(IT_0254)*(IT_0335*IT_0336 + IT_0334*IT_0342) + IT_0343*
      (IT_0216*IT_0335 + IT_0217*IT_0337 + IT_0241*IT_0342 + IT_0245*IT_0344) +
       2*IT_0348*IT_0349 + (IT_0183*IT_0241 + IT_0186*IT_0245 + conj(IT_0215)
      *IT_0257 + IT_0215*IT_0260 + conj(IT_0180)*IT_0275 + IT_0180*IT_0278 +
       conj(IT_0166)*IT_0297 + IT_0166*IT_0300)*IT_0350 + (conj(IT_0254)*IT_0297
       + IT_0254*IT_0300 + IT_0241*IT_0335 + IT_0245*IT_0337 + conj(IT_0215)
      *IT_0342 + IT_0215*IT_0344)*IT_0351;
    return create_ccomplex_return(IT_0352);
}

