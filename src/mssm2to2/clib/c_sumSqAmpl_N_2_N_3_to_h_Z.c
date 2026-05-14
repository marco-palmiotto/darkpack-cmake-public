#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_3_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_3_to_h_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t IT_0000 = m_N_2*m_N_3;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = sin(alpha);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0019;
    const ccomplex_t IT_0024 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0013 + IT_0017 + 
      -IT_0022 + -IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0010*IT_0018;
    const ccomplex_t IT_0029 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = IT_0009*IT_0019;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = IT_0031*IT_0033;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0030 + IT_0032 + 
      -IT_0034 + -IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0027*IT_0037;
    const ccomplex_t IT_0039 = pow(m_h, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0046 = IT_0015*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0048 = IT_0011*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0050 = IT_0015*IT_0049;
    const ccomplex_t IT_0051 = IT_0044 + IT_0046 + IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0053 = IT_0020*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0055 = IT_0023*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0059 = IT_0023*IT_0058;
    const ccomplex_t IT_0060 = -IT_0053 + -IT_0055 + -IT_0057 + -IT_0059;
    const ccomplex_t IT_0061 = IT_0051 + IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0065 = IT_0028*IT_0064;
    const ccomplex_t IT_0066 = IT_0031*IT_0064;
    const ccomplex_t IT_0067 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0068 = IT_0028*IT_0067;
    const ccomplex_t IT_0069 = IT_0031*IT_0067;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0065 + IT_0066 + 
      -IT_0068 + -IT_0069);
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = pow(m_N_3, 2);
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0039 + IT_0073 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0072*IT_0075;
    const ccomplex_t IT_0077 = (-2)*IT_0042 + (-2)*IT_0076;
    const ccomplex_t IT_0078 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0079 = IT_0011*IT_0078;
    const ccomplex_t IT_0080 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0081 = IT_0015*IT_0080;
    const ccomplex_t IT_0082 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0085 = IT_0015*IT_0084;
    const ccomplex_t IT_0086 = IT_0079 + IT_0081 + IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0088 = IT_0020*IT_0087;
    const ccomplex_t IT_0089 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0090 = IT_0023*IT_0089;
    const ccomplex_t IT_0091 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0092 = IT_0020*IT_0091;
    const ccomplex_t IT_0093 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0094 = IT_0023*IT_0093;
    const ccomplex_t IT_0095 = -IT_0088 + -IT_0090 + -IT_0092 + -IT_0094;
    const ccomplex_t IT_0096 = IT_0086 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0100 = IT_0028*IT_0099;
    const ccomplex_t IT_0101 = IT_0031*IT_0099;
    const ccomplex_t IT_0102 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0103 = IT_0028*IT_0102;
    const ccomplex_t IT_0104 = IT_0031*IT_0102;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0100 + IT_0101 + 
      -IT_0103 + -IT_0104);
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0098*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0039 + IT_0073 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = cos(beta);
    const ccomplex_t IT_0112 = IT_0010*IT_0111;
    const ccomplex_t IT_0113 = IT_0080*IT_0112;
    const ccomplex_t IT_0114 = sin(beta);
    const ccomplex_t IT_0115 = IT_0010*IT_0114;
    const ccomplex_t IT_0116 = IT_0078*IT_0115;
    const ccomplex_t IT_0117 = IT_0019*IT_0111;
    const ccomplex_t IT_0118 = IT_0089*IT_0117;
    const ccomplex_t IT_0119 = IT_0019*IT_0114;
    const ccomplex_t IT_0120 = IT_0087*IT_0119;
    const ccomplex_t IT_0121 = IT_0084*IT_0112;
    const ccomplex_t IT_0122 = IT_0082*IT_0115;
    const ccomplex_t IT_0123 = IT_0093*IT_0117;
    const ccomplex_t IT_0124 = IT_0091*IT_0119;
    const ccomplex_t IT_0125 = (-0.5)*IT_0113 + (-0.5)*IT_0116 + 0.5*IT_0118 +
       0.5*IT_0120 + (-0.5)*IT_0121 + (-0.5)*IT_0122 + 0.5*IT_0123 + 0.5*IT_0124;
    const ccomplex_t IT_0126 = -e_em*IT_0010*IT_0019*(IT_0014*IT_0111 + 
      -IT_0008*IT_0114);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = pow(m_Z, 2);
    const ccomplex_t IT_0130 = cpow((-2)*s_12 + -IT_0004 + -IT_0073 + IT_0129 
      + -reg_prop, -1);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0128*IT_0131;
    const ccomplex_t IT_0133 = IT_0078*IT_0112;
    const ccomplex_t IT_0134 = IT_0080*IT_0115;
    const ccomplex_t IT_0135 = IT_0087*IT_0117;
    const ccomplex_t IT_0136 = IT_0089*IT_0119;
    const ccomplex_t IT_0137 = IT_0082*IT_0112;
    const ccomplex_t IT_0138 = IT_0084*IT_0115;
    const ccomplex_t IT_0139 = IT_0091*IT_0117;
    const ccomplex_t IT_0140 = IT_0093*IT_0119;
    const ccomplex_t IT_0141 = (-0.5)*IT_0133 + 0.5*IT_0134 + 0.5*IT_0135 + (
      -0.5)*IT_0136 + (-0.5)*IT_0137 + 0.5*IT_0138 + 0.5*IT_0139 + (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = e_em*IT_0010*IT_0019*(IT_0008*IT_0111 + IT_0014
      *IT_0114);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_12 + -IT_0004 + -IT_0073 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0105;
    const ccomplex_t IT_0149 = IT_0098*IT_0148;
    const ccomplex_t IT_0150 = IT_0109*IT_0149;
    const ccomplex_t IT_0151 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0152 = IT_0011*IT_0151;
    const ccomplex_t IT_0153 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0154 = IT_0015*IT_0153;
    const ccomplex_t IT_0155 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0156 = IT_0011*IT_0155;
    const ccomplex_t IT_0157 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0158 = IT_0015*IT_0157;
    const ccomplex_t IT_0159 = IT_0152 + IT_0154 + IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0161 = IT_0020*IT_0160;
    const ccomplex_t IT_0162 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0163 = IT_0023*IT_0162;
    const ccomplex_t IT_0164 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0165 = IT_0020*IT_0164;
    const ccomplex_t IT_0166 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0167 = IT_0023*IT_0166;
    const ccomplex_t IT_0168 = -IT_0161 + -IT_0163 + -IT_0165 + -IT_0167;
    const ccomplex_t IT_0169 = IT_0159 + IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0173 = IT_0028*IT_0172;
    const ccomplex_t IT_0174 = IT_0031*IT_0172;
    const ccomplex_t IT_0175 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0176 = IT_0028*IT_0175;
    const ccomplex_t IT_0177 = IT_0031*IT_0175;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0173 + IT_0174 + 
      -IT_0176 + -IT_0177);
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0171*IT_0179;
    const ccomplex_t IT_0181 = cpow((-2)*s_23 + IT_0039 + IT_0073 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = IT_0180*IT_0182;
    const ccomplex_t IT_0184 = -IT_0110 + IT_0132 + IT_0147 + IT_0150 + (-2)
      *IT_0183;
    const ccomplex_t IT_0185 = 2*IT_0184;
    const ccomplex_t IT_0186 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0187 = IT_0011*IT_0186;
    const ccomplex_t IT_0188 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0189 = IT_0015*IT_0188;
    const ccomplex_t IT_0190 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0191 = IT_0020*IT_0190;
    const ccomplex_t IT_0192 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0193 = IT_0023*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0187 + IT_0189 + 
      -IT_0191 + -IT_0193);
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0197 = IT_0028*IT_0196;
    const ccomplex_t IT_0198 = IT_0031*IT_0196;
    const ccomplex_t IT_0199 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0200 = IT_0028*IT_0199;
    const ccomplex_t IT_0201 = IT_0031*IT_0199;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0197 + IT_0198 + 
      -IT_0200 + -IT_0201);
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0195*IT_0203;
    const ccomplex_t IT_0205 = IT_0041*IT_0204;
    const ccomplex_t IT_0206 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0207 = IT_0011*IT_0206;
    const ccomplex_t IT_0208 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0209 = IT_0015*IT_0208;
    const ccomplex_t IT_0210 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0211 = IT_0011*IT_0210;
    const ccomplex_t IT_0212 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0213 = IT_0015*IT_0212;
    const ccomplex_t IT_0214 = IT_0207 + IT_0209 + IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0216 = IT_0020*IT_0215;
    const ccomplex_t IT_0217 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0218 = IT_0023*IT_0217;
    const ccomplex_t IT_0219 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0220 = IT_0020*IT_0219;
    const ccomplex_t IT_0221 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0222 = IT_0023*IT_0221;
    const ccomplex_t IT_0223 = -IT_0216 + -IT_0218 + -IT_0220 + -IT_0222;
    const ccomplex_t IT_0224 = IT_0214 + IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0228 = IT_0028*IT_0227;
    const ccomplex_t IT_0229 = IT_0031*IT_0227;
    const ccomplex_t IT_0230 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0231 = IT_0028*IT_0230;
    const ccomplex_t IT_0232 = IT_0031*IT_0230;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0228 + IT_0229 + 
      -IT_0231 + -IT_0232);
    const ccomplex_t IT_0234 = (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = IT_0226*IT_0234;
    const ccomplex_t IT_0236 = IT_0075*IT_0235;
    const ccomplex_t IT_0237 = (-2)*IT_0205 + (-2)*IT_0236;
    const ccomplex_t IT_0238 = (-2)*IT_0237;
    const ccomplex_t IT_0239 = IT_0185 + IT_0238;
    const ccomplex_t IT_0240 = 2*conj(IT_0184);
    const ccomplex_t IT_0241 = (-2)*conj(IT_0237);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0244 = IT_0112*IT_0243;
    const ccomplex_t IT_0245 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0246 = IT_0115*IT_0245;
    const ccomplex_t IT_0247 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0248 = IT_0117*IT_0247;
    const ccomplex_t IT_0249 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0250 = IT_0119*IT_0249;
    const ccomplex_t IT_0251 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0252 = IT_0112*IT_0251;
    const ccomplex_t IT_0253 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0254 = IT_0115*IT_0253;
    const ccomplex_t IT_0255 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0256 = IT_0117*IT_0255;
    const ccomplex_t IT_0257 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0258 = IT_0119*IT_0257;
    const ccomplex_t IT_0259 = 0.5*IT_0244 + 0.5*IT_0246 + (-0.5)*IT_0248 + (
      -0.5)*IT_0250 + 0.5*IT_0252 + 0.5*IT_0254 + (-0.5)*IT_0256 + (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0127*IT_0259;
    const ccomplex_t IT_0261 = IT_0131*IT_0260;
    const ccomplex_t IT_0262 = IT_0011*IT_0245;
    const ccomplex_t IT_0263 = IT_0015*IT_0243;
    const ccomplex_t IT_0264 = IT_0011*IT_0253;
    const ccomplex_t IT_0265 = IT_0015*IT_0251;
    const ccomplex_t IT_0266 = IT_0262 + IT_0263 + IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = IT_0020*IT_0249;
    const ccomplex_t IT_0268 = IT_0023*IT_0247;
    const ccomplex_t IT_0269 = IT_0020*IT_0257;
    const ccomplex_t IT_0270 = IT_0023*IT_0255;
    const ccomplex_t IT_0271 = -IT_0267 + -IT_0268 + -IT_0269 + -IT_0270;
    const ccomplex_t IT_0272 = IT_0266 + IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = (-0.5)*IT_0273;
    const ccomplex_t IT_0275 = IT_0106*IT_0274;
    const ccomplex_t IT_0276 = IT_0109*IT_0275;
    const ccomplex_t IT_0277 = IT_0112*IT_0245;
    const ccomplex_t IT_0278 = IT_0115*IT_0243;
    const ccomplex_t IT_0279 = IT_0117*IT_0249;
    const ccomplex_t IT_0280 = IT_0119*IT_0247;
    const ccomplex_t IT_0281 = IT_0112*IT_0253;
    const ccomplex_t IT_0282 = IT_0115*IT_0251;
    const ccomplex_t IT_0283 = IT_0117*IT_0257;
    const ccomplex_t IT_0284 = IT_0119*IT_0255;
    const ccomplex_t IT_0285 = 0.5*IT_0277 + (-0.5)*IT_0278 + (-0.5)*IT_0279 +
       0.5*IT_0280 + 0.5*IT_0281 + (-0.5)*IT_0282 + (-0.5)*IT_0283 + 0.5*IT_0284;
    const ccomplex_t IT_0286 = IT_0143*IT_0285;
    const ccomplex_t IT_0287 = IT_0146*IT_0286;
    const ccomplex_t IT_0288 = IT_0148*IT_0274;
    const ccomplex_t IT_0289 = IT_0109*IT_0288;
    const ccomplex_t IT_0290 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0291 = IT_0011*IT_0290;
    const ccomplex_t IT_0292 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0293 = IT_0015*IT_0292;
    const ccomplex_t IT_0294 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0295 = IT_0011*IT_0294;
    const ccomplex_t IT_0296 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0297 = IT_0015*IT_0296;
    const ccomplex_t IT_0298 = IT_0291 + IT_0293 + IT_0295 + IT_0297;
    const ccomplex_t IT_0299 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0300 = IT_0020*IT_0299;
    const ccomplex_t IT_0301 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0302 = IT_0023*IT_0301;
    const ccomplex_t IT_0303 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0304 = IT_0020*IT_0303;
    const ccomplex_t IT_0305 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0306 = IT_0023*IT_0305;
    const ccomplex_t IT_0307 = -IT_0300 + -IT_0302 + -IT_0304 + -IT_0306;
    const ccomplex_t IT_0308 = IT_0298 + IT_0307;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0312 = IT_0028*IT_0311;
    const ccomplex_t IT_0313 = IT_0031*IT_0311;
    const ccomplex_t IT_0314 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0315 = IT_0028*IT_0314;
    const ccomplex_t IT_0316 = IT_0031*IT_0314;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*(IT_0312 + IT_0313 + 
      -IT_0315 + -IT_0316);
    const ccomplex_t IT_0318 = (-0.5)*IT_0317;
    const ccomplex_t IT_0319 = IT_0310*IT_0318;
    const ccomplex_t IT_0320 = IT_0182*IT_0319;
    const ccomplex_t IT_0321 = IT_0261 + IT_0276 + IT_0287 + -IT_0289 + (-2)
      *IT_0320;
    const ccomplex_t IT_0322 = 2*IT_0237;
    const ccomplex_t IT_0323 = (-2)*IT_0184;
    const ccomplex_t IT_0324 = IT_0322 + IT_0323;
    const ccomplex_t IT_0325 = 2*conj(IT_0237);
    const ccomplex_t IT_0326 = (-2)*conj(IT_0184);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = s_12 + IT_0329;
    const ccomplex_t IT_0331 = s_12*IT_0330;
    const ccomplex_t IT_0332 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0333 = IT_0028*IT_0332;
    const ccomplex_t IT_0334 = IT_0031*IT_0332;
    const ccomplex_t IT_0335 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0336 = IT_0028*IT_0335;
    const ccomplex_t IT_0337 = IT_0031*IT_0335;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*(IT_0333 + IT_0334 + 
      -IT_0336 + -IT_0337);
    const ccomplex_t IT_0339 = 0.5*IT_0338;
    const ccomplex_t IT_0340 = IT_0274*IT_0339;
    const ccomplex_t IT_0341 = cpow((-2)*s_13 + IT_0004 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*IT_0341;
    const ccomplex_t IT_0343 = IT_0340*IT_0342;
    const ccomplex_t IT_0344 = cpow((-2)*s_13 + IT_0004 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*IT_0344;
    const ccomplex_t IT_0346 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0347 = IT_0011*IT_0346;
    const ccomplex_t IT_0348 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0349 = IT_0015*IT_0348;
    const ccomplex_t IT_0350 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0351 = IT_0011*IT_0350;
    const ccomplex_t IT_0352 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0353 = IT_0015*IT_0352;
    const ccomplex_t IT_0354 = IT_0347 + IT_0349 + IT_0351 + IT_0353;
    const ccomplex_t IT_0355 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0356 = IT_0020*IT_0355;
    const ccomplex_t IT_0357 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0358 = IT_0023*IT_0357;
    const ccomplex_t IT_0359 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0360 = IT_0020*IT_0359;
    const ccomplex_t IT_0361 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0362 = IT_0023*IT_0361;
    const ccomplex_t IT_0363 = -IT_0356 + -IT_0358 + -IT_0360 + -IT_0362;
    const ccomplex_t IT_0364 = IT_0354 + IT_0363;
    const ccomplex_t IT_0365 = (0 + _Complex_I*1)*IT_0364;
    const ccomplex_t IT_0366 = (-0.5)*IT_0365;
    const ccomplex_t IT_0367 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0368 = IT_0028*IT_0367;
    const ccomplex_t IT_0369 = IT_0031*IT_0367;
    const ccomplex_t IT_0370 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0371 = IT_0028*IT_0370;
    const ccomplex_t IT_0372 = IT_0031*IT_0370;
    const ccomplex_t IT_0373 = (0 + _Complex_I*1)*(IT_0368 + IT_0369 + 
      -IT_0371 + -IT_0372);
    const ccomplex_t IT_0374 = (-0.5)*IT_0373;
    const ccomplex_t IT_0375 = IT_0366*IT_0374;
    const ccomplex_t IT_0376 = IT_0345*IT_0375;
    const ccomplex_t IT_0377 = IT_0261 + IT_0287 + IT_0343 + 2*IT_0376;
    const ccomplex_t IT_0378 = 2*IT_0077;
    const ccomplex_t IT_0379 = (-2)*IT_0321;
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = 2*conj(IT_0077);
    const ccomplex_t IT_0382 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0385 = IT_0011*IT_0384;
    const ccomplex_t IT_0386 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0387 = IT_0015*IT_0386;
    const ccomplex_t IT_0388 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0389 = IT_0020*IT_0388;
    const ccomplex_t IT_0390 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0391 = IT_0023*IT_0390;
    const ccomplex_t IT_0392 = (0 + _Complex_I*1)*(IT_0385 + IT_0387 + 
      -IT_0389 + -IT_0391);
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = IT_0037*IT_0393;
    const ccomplex_t IT_0395 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*IT_0395;
    const ccomplex_t IT_0397 = IT_0394*IT_0396;
    const ccomplex_t IT_0398 = IT_0098*IT_0339;
    const ccomplex_t IT_0399 = IT_0342*IT_0398;
    const ccomplex_t IT_0400 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0401 = IT_0011*IT_0400;
    const ccomplex_t IT_0402 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0403 = IT_0015*IT_0402;
    const ccomplex_t IT_0404 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0405 = IT_0011*IT_0404;
    const ccomplex_t IT_0406 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0407 = IT_0015*IT_0406;
    const ccomplex_t IT_0408 = IT_0401 + IT_0403 + IT_0405 + IT_0407;
    const ccomplex_t IT_0409 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0410 = IT_0020*IT_0409;
    const ccomplex_t IT_0411 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0412 = IT_0023*IT_0411;
    const ccomplex_t IT_0413 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0414 = IT_0020*IT_0413;
    const ccomplex_t IT_0415 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0416 = IT_0023*IT_0415;
    const ccomplex_t IT_0417 = -IT_0410 + -IT_0412 + -IT_0414 + -IT_0416;
    const ccomplex_t IT_0418 = IT_0408 + IT_0417;
    const ccomplex_t IT_0419 = (0 + _Complex_I*1)*IT_0418;
    const ccomplex_t IT_0420 = (-0.5)*IT_0419;
    const ccomplex_t IT_0421 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0422 = IT_0028*IT_0421;
    const ccomplex_t IT_0423 = IT_0031*IT_0421;
    const ccomplex_t IT_0424 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0425 = IT_0028*IT_0424;
    const ccomplex_t IT_0426 = IT_0031*IT_0424;
    const ccomplex_t IT_0427 = (0 + _Complex_I*1)*(IT_0422 + IT_0423 + 
      -IT_0425 + -IT_0426);
    const ccomplex_t IT_0428 = (-0.5)*IT_0427;
    const ccomplex_t IT_0429 = IT_0420*IT_0428;
    const ccomplex_t IT_0430 = cpow((-2)*s_13 + IT_0004 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0431 = (0 + _Complex_I*1)*IT_0430;
    const ccomplex_t IT_0432 = IT_0429*IT_0431;
    const ccomplex_t IT_0433 = 2*IT_0397 + IT_0399 + 2*IT_0432;
    const ccomplex_t IT_0434 = (-0.5)*IT_0338;
    const ccomplex_t IT_0435 = IT_0098*IT_0434;
    const ccomplex_t IT_0436 = IT_0342*IT_0435;
    const ccomplex_t IT_0437 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0438 = IT_0011*IT_0437;
    const ccomplex_t IT_0439 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0440 = IT_0015*IT_0439;
    const ccomplex_t IT_0441 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0442 = IT_0011*IT_0441;
    const ccomplex_t IT_0443 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0444 = IT_0015*IT_0443;
    const ccomplex_t IT_0445 = IT_0438 + IT_0440 + IT_0442 + IT_0444;
    const ccomplex_t IT_0446 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0447 = IT_0020*IT_0446;
    const ccomplex_t IT_0448 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0449 = IT_0023*IT_0448;
    const ccomplex_t IT_0450 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0451 = IT_0020*IT_0450;
    const ccomplex_t IT_0452 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0453 = IT_0023*IT_0452;
    const ccomplex_t IT_0454 = -IT_0447 + -IT_0449 + -IT_0451 + -IT_0453;
    const ccomplex_t IT_0455 = IT_0445 + IT_0454;
    const ccomplex_t IT_0456 = (0 + _Complex_I*1)*IT_0455;
    const ccomplex_t IT_0457 = (-0.5)*IT_0456;
    const ccomplex_t IT_0458 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0459 = IT_0028*IT_0458;
    const ccomplex_t IT_0460 = IT_0031*IT_0458;
    const ccomplex_t IT_0461 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0462 = IT_0028*IT_0461;
    const ccomplex_t IT_0463 = IT_0031*IT_0461;
    const ccomplex_t IT_0464 = (0 + _Complex_I*1)*(IT_0459 + IT_0460 + 
      -IT_0462 + -IT_0463);
    const ccomplex_t IT_0465 = (-0.5)*IT_0464;
    const ccomplex_t IT_0466 = IT_0457*IT_0465;
    const ccomplex_t IT_0467 = IT_0345*IT_0466;
    const ccomplex_t IT_0468 = IT_0132 + IT_0147 + IT_0436 + 2*IT_0467;
    const ccomplex_t IT_0469 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0470 = IT_0011*IT_0469;
    const ccomplex_t IT_0471 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0472 = IT_0015*IT_0471;
    const ccomplex_t IT_0473 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0474 = IT_0011*IT_0473;
    const ccomplex_t IT_0475 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0476 = IT_0015*IT_0475;
    const ccomplex_t IT_0477 = IT_0470 + IT_0472 + IT_0474 + IT_0476;
    const ccomplex_t IT_0478 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0479 = IT_0020*IT_0478;
    const ccomplex_t IT_0480 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0481 = IT_0023*IT_0480;
    const ccomplex_t IT_0482 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0483 = IT_0020*IT_0482;
    const ccomplex_t IT_0484 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0485 = IT_0023*IT_0484;
    const ccomplex_t IT_0486 = -IT_0479 + -IT_0481 + -IT_0483 + -IT_0485;
    const ccomplex_t IT_0487 = IT_0477 + IT_0486;
    const ccomplex_t IT_0488 = (0 + _Complex_I*1)*IT_0487;
    const ccomplex_t IT_0489 = (-0.5)*IT_0488;
    const ccomplex_t IT_0490 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0491 = IT_0028*IT_0490;
    const ccomplex_t IT_0492 = IT_0031*IT_0490;
    const ccomplex_t IT_0493 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0494 = IT_0028*IT_0493;
    const ccomplex_t IT_0495 = IT_0031*IT_0493;
    const ccomplex_t IT_0496 = (0 + _Complex_I*1)*(IT_0491 + IT_0492 + 
      -IT_0494 + -IT_0495);
    const ccomplex_t IT_0497 = (-0.5)*IT_0496;
    const ccomplex_t IT_0498 = IT_0489*IT_0497;
    const ccomplex_t IT_0499 = IT_0431*IT_0498;
    const ccomplex_t IT_0500 = IT_0274*IT_0434;
    const ccomplex_t IT_0501 = IT_0342*IT_0500;
    const ccomplex_t IT_0502 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0503 = IT_0011*IT_0502;
    const ccomplex_t IT_0504 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0505 = IT_0015*IT_0504;
    const ccomplex_t IT_0506 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0507 = IT_0020*IT_0506;
    const ccomplex_t IT_0508 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0509 = IT_0023*IT_0508;
    const ccomplex_t IT_0510 = (0 + _Complex_I*1)*(IT_0503 + IT_0505 + 
      -IT_0507 + -IT_0509);
    const ccomplex_t IT_0511 = -IT_0510;
    const ccomplex_t IT_0512 = IT_0203*IT_0511;
    const ccomplex_t IT_0513 = IT_0396*IT_0512;
    const ccomplex_t IT_0514 = 2*IT_0499 + IT_0501 + 2*IT_0513;
    const ccomplex_t IT_0515 = 2*IT_0321;
    const ccomplex_t IT_0516 = (-2)*IT_0077;
    const ccomplex_t IT_0517 = IT_0515 + IT_0516;
    const ccomplex_t IT_0518 = 2*conj(IT_0321);
    const ccomplex_t IT_0519 = (-2)*conj(IT_0077);
    const ccomplex_t IT_0520 = IT_0518 + IT_0519;
    const ccomplex_t IT_0521 = s_12*s_14;
    const ccomplex_t IT_0522 = s_24*IT_0004;
    const ccomplex_t IT_0523 = -IT_0522;
    const ccomplex_t IT_0524 = IT_0521 + IT_0523;
    const ccomplex_t IT_0525 = (-0.5)*IT_0343 + -IT_0376;
    const ccomplex_t IT_0526 = conj(IT_0042) + conj(IT_0076) + conj(IT_0525);
    const ccomplex_t IT_0527 = IT_0042 + IT_0076 + IT_0525;
    const ccomplex_t IT_0528 = (-0.5)*IT_0399;
    const ccomplex_t IT_0529 = 0.5*IT_0110 + (-0.5)*IT_0150 + IT_0183;
    const ccomplex_t IT_0530 = -IT_0397 + -IT_0432;
    const ccomplex_t IT_0531 = (-0.5)*IT_0436 + -IT_0467;
    const ccomplex_t IT_0532 = -IT_0513;
    const ccomplex_t IT_0533 = (-0.5)*IT_0276 + 0.5*IT_0289 + IT_0320;
    const ccomplex_t IT_0534 = -IT_0499 + (-0.5)*IT_0501;
    const ccomplex_t IT_0535 = IT_0380*IT_0526 + IT_0383*IT_0527 + IT_0242*
      (IT_0528 + IT_0529 + IT_0530) + IT_0239*(conj(IT_0528) + conj(IT_0529) +
       conj(IT_0530)) + IT_0327*(IT_0205 + IT_0236 + IT_0531) + IT_0324*(conj
      (IT_0205) + conj(IT_0236) + conj(IT_0531)) + IT_0520*(IT_0532 + IT_0533 +
       IT_0534) + IT_0517*(conj(IT_0532) + conj(IT_0533) + conj(IT_0534));
    const ccomplex_t IT_0536 = m_N_2*IT_0330;
    const ccomplex_t IT_0537 = IT_0171*IT_0318;
    const ccomplex_t IT_0538 = m_N_1*IT_0537;
    const ccomplex_t IT_0539 = m_N_2*IT_0180;
    const ccomplex_t IT_0540 = -IT_0539;
    const ccomplex_t IT_0541 = IT_0538 + IT_0540;
    const ccomplex_t IT_0542 = IT_0182*IT_0541;
    const ccomplex_t IT_0543 = m_N_3*IT_0340;
    const ccomplex_t IT_0544 = m_N_3*IT_0500;
    const ccomplex_t IT_0545 = -IT_0544;
    const ccomplex_t IT_0546 = IT_0543 + IT_0545;
    const ccomplex_t IT_0547 = IT_0342*IT_0546;
    const ccomplex_t IT_0548 = IT_0366*IT_0465;
    const ccomplex_t IT_0549 = m_N_4*IT_0548;
    const ccomplex_t IT_0550 = m_N_3*IT_0375;
    const ccomplex_t IT_0551 = -IT_0550;
    const ccomplex_t IT_0552 = IT_0549 + IT_0551;
    const ccomplex_t IT_0553 = IT_0345*IT_0552;
    const ccomplex_t IT_0554 = m_N_2*IT_0149;
    const ccomplex_t IT_0555 = IT_0109*IT_0554;
    const ccomplex_t IT_0556 = m_N_2*IT_0107;
    const ccomplex_t IT_0557 = -IT_0556;
    const ccomplex_t IT_0558 = IT_0554 + IT_0557;
    const ccomplex_t IT_0559 = IT_0109*IT_0558;
    const ccomplex_t IT_0560 = IT_0542 + (-0.5)*IT_0547 + IT_0553 + 0.5
      *IT_0555 + 0.5*IT_0559;
    const ccomplex_t IT_0561 = IT_0109*IT_0556;
    const ccomplex_t IT_0562 = IT_0037*IT_0195;
    const ccomplex_t IT_0563 = m_N_3*IT_0562;
    const ccomplex_t IT_0564 = m_N_2*IT_0204;
    const ccomplex_t IT_0565 = -IT_0564;
    const ccomplex_t IT_0566 = IT_0563 + IT_0565;
    const ccomplex_t IT_0567 = IT_0041*IT_0566;
    const ccomplex_t IT_0568 = m_N_3*IT_0512;
    const ccomplex_t IT_0569 = -IT_0568;
    const ccomplex_t IT_0570 = IT_0037*IT_0511;
    const ccomplex_t IT_0571 = m_N_2*IT_0570;
    const ccomplex_t IT_0572 = IT_0569 + IT_0571;
    const ccomplex_t IT_0573 = IT_0396*IT_0572;
    const ccomplex_t IT_0574 = m_W*e_em;
    const ccomplex_t IT_0575 = cpow(IT_0009, -4);
    const ccomplex_t IT_0576 = IT_0008*IT_0018*IT_0114*IT_0575;
    const ccomplex_t IT_0577 = IT_0014*IT_0018*IT_0111*IT_0575;
    const ccomplex_t IT_0578 = cpow(IT_0018, -3);
    const ccomplex_t IT_0579 = IT_0008*IT_0114*IT_0578;
    const ccomplex_t IT_0580 = IT_0014*IT_0111*IT_0578;
    const ccomplex_t IT_0581 = cpow(IT_0009, -2);
    const ccomplex_t IT_0582 = IT_0008*IT_0019*IT_0114*IT_0581;
    const ccomplex_t IT_0583 = IT_0014*IT_0019*IT_0111*IT_0581;
    const ccomplex_t IT_0584 = cpow(IT_0018, -2);
    const ccomplex_t IT_0585 = cpow(IT_0581 + IT_0584, -1);
    const ccomplex_t IT_0586 = (IT_0576 + -IT_0577 + IT_0579 + -IT_0580 + 2
      *IT_0582 + (-2)*IT_0583)*IT_0585;
    const ccomplex_t IT_0587 = IT_0574*IT_0586;
    const ccomplex_t IT_0588 = (0 + _Complex_I*1)*IT_0587;
    const ccomplex_t IT_0589 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0004 + -IT_0073 + -reg_prop, -1);
    const ccomplex_t IT_0590 = IT_0588*IT_0589;
    const ccomplex_t IT_0591 = (0 + _Complex_I*1)*IT_0590;
    const ccomplex_t IT_0592 = IT_0037*IT_0591;
    const ccomplex_t IT_0593 = IT_0428*IT_0489;
    const ccomplex_t IT_0594 = m_N_1*IT_0593;
    const ccomplex_t IT_0595 = m_N_3*IT_0498;
    const ccomplex_t IT_0596 = -IT_0595;
    const ccomplex_t IT_0597 = IT_0594 + IT_0596;
    const ccomplex_t IT_0598 = IT_0431*IT_0597;
    const ccomplex_t IT_0599 = IT_0071*IT_0226;
    const ccomplex_t IT_0600 = m_N_4*IT_0599;
    const ccomplex_t IT_0601 = m_N_2*IT_0235;
    const ccomplex_t IT_0602 = -IT_0601;
    const ccomplex_t IT_0603 = IT_0600 + IT_0602;
    const ccomplex_t IT_0604 = IT_0075*IT_0603;
    const ccomplex_t IT_0605 = (-0.5)*IT_0547 + (-0.5)*IT_0561 + -IT_0567 + 
      -IT_0573 + -IT_0592 + -IT_0598 + -IT_0604;
    const ccomplex_t IT_0606 = conj(IT_0560) + conj(IT_0605);
    const ccomplex_t IT_0607 = IT_0560 + IT_0605;
    const ccomplex_t IT_0608 = m_N_2*IT_0275;
    const ccomplex_t IT_0609 = IT_0109*IT_0608;
    const ccomplex_t IT_0610 = m_N_2*IT_0288;
    const ccomplex_t IT_0611 = -IT_0610;
    const ccomplex_t IT_0612 = IT_0608 + IT_0611;
    const ccomplex_t IT_0613 = IT_0109*IT_0612;
    const ccomplex_t IT_0614 = m_N_3*IT_0435;
    const ccomplex_t IT_0615 = m_N_3*IT_0398;
    const ccomplex_t IT_0616 = -IT_0615;
    const ccomplex_t IT_0617 = IT_0614 + IT_0616;
    const ccomplex_t IT_0618 = IT_0342*IT_0617;
    const ccomplex_t IT_0619 = IT_0179*IT_0310;
    const ccomplex_t IT_0620 = m_N_1*IT_0619;
    const ccomplex_t IT_0621 = m_N_2*IT_0319;
    const ccomplex_t IT_0622 = -IT_0621;
    const ccomplex_t IT_0623 = IT_0620 + IT_0622;
    const ccomplex_t IT_0624 = IT_0182*IT_0623;
    const ccomplex_t IT_0625 = IT_0374*IT_0457;
    const ccomplex_t IT_0626 = m_N_4*IT_0625;
    const ccomplex_t IT_0627 = m_N_3*IT_0466;
    const ccomplex_t IT_0628 = -IT_0627;
    const ccomplex_t IT_0629 = IT_0626 + IT_0628;
    const ccomplex_t IT_0630 = IT_0345*IT_0629;
    const ccomplex_t IT_0631 = 0.5*IT_0609 + 0.5*IT_0613 + (-0.5)*IT_0618 +
       IT_0624 + IT_0630;
    const ccomplex_t IT_0632 = IT_0109*IT_0610;
    const ccomplex_t IT_0633 = IT_0420*IT_0497;
    const ccomplex_t IT_0634 = m_N_1*IT_0633;
    const ccomplex_t IT_0635 = m_N_3*IT_0429;
    const ccomplex_t IT_0636 = -IT_0635;
    const ccomplex_t IT_0637 = IT_0634 + IT_0636;
    const ccomplex_t IT_0638 = IT_0431*IT_0637;
    const ccomplex_t IT_0639 = IT_0203*IT_0591;
    const ccomplex_t IT_0640 = IT_0027*IT_0203;
    const ccomplex_t IT_0641 = m_N_3*IT_0640;
    const ccomplex_t IT_0642 = m_N_2*IT_0038;
    const ccomplex_t IT_0643 = -IT_0642;
    const ccomplex_t IT_0644 = IT_0641 + IT_0643;
    const ccomplex_t IT_0645 = IT_0041*IT_0644;
    const ccomplex_t IT_0646 = IT_0203*IT_0393;
    const ccomplex_t IT_0647 = m_N_2*IT_0646;
    const ccomplex_t IT_0648 = m_N_3*IT_0394;
    const ccomplex_t IT_0649 = -IT_0648;
    const ccomplex_t IT_0650 = IT_0647 + IT_0649;
    const ccomplex_t IT_0651 = IT_0396*IT_0650;
    const ccomplex_t IT_0652 = IT_0063*IT_0234;
    const ccomplex_t IT_0653 = m_N_4*IT_0652;
    const ccomplex_t IT_0654 = m_N_2*IT_0072;
    const ccomplex_t IT_0655 = -IT_0654;
    const ccomplex_t IT_0656 = IT_0653 + IT_0655;
    const ccomplex_t IT_0657 = IT_0075*IT_0656;
    const ccomplex_t IT_0658 = (-0.5)*IT_0618 + (-0.5)*IT_0632 + -IT_0638 + 
      -IT_0639 + -IT_0645 + -IT_0651 + -IT_0657;
    const ccomplex_t IT_0659 = conj(IT_0631) + conj(IT_0658);
    const ccomplex_t IT_0660 = IT_0631 + IT_0658;
    const ccomplex_t IT_0661 = IT_0000*IT_0330;
    const ccomplex_t IT_0662 = 2*IT_0468;
    const ccomplex_t IT_0663 = (-2)*IT_0433;
    const ccomplex_t IT_0664 = IT_0662 + IT_0663;
    const ccomplex_t IT_0665 = 2*conj(IT_0468);
    const ccomplex_t IT_0666 = (-2)*conj(IT_0433);
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = 2*IT_0433;
    const ccomplex_t IT_0669 = (-2)*IT_0468;
    const ccomplex_t IT_0670 = IT_0668 + IT_0669;
    const ccomplex_t IT_0671 = 2*conj(IT_0433);
    const ccomplex_t IT_0672 = (-2)*conj(IT_0468);
    const ccomplex_t IT_0673 = IT_0671 + IT_0672;
    const ccomplex_t IT_0674 = pow(s_24, 2);
    const ccomplex_t IT_0675 = IT_0001*IT_0674;
    const ccomplex_t IT_0676 = -IT_0073;
    const ccomplex_t IT_0677 = IT_0675 + IT_0676;
    const ccomplex_t IT_0678 = IT_0000*IT_0677;
    const ccomplex_t IT_0679 = m_N_2*IT_0677;
    const ccomplex_t IT_0680 = 2*IT_0514;
    const ccomplex_t IT_0681 = (-2)*IT_0377;
    const ccomplex_t IT_0682 = IT_0680 + IT_0681;
    const ccomplex_t IT_0683 = 2*conj(IT_0514);
    const ccomplex_t IT_0684 = (-2)*conj(IT_0377);
    const ccomplex_t IT_0685 = IT_0683 + IT_0684;
    const ccomplex_t IT_0686 = s_24*IT_0129;
    const ccomplex_t IT_0687 = IT_0001*IT_0686;
    const ccomplex_t IT_0688 = 0.5*IT_0687;
    const ccomplex_t IT_0689 = s_24 + IT_0688;
    const ccomplex_t IT_0690 = m_N_2*IT_0689;
    const ccomplex_t IT_0691 = 2*IT_0690;
    const ccomplex_t IT_0692 = conj(IT_0528) + conj(IT_0530);
    const ccomplex_t IT_0693 = 2*IT_0560;
    const ccomplex_t IT_0694 = 2*IT_0605;
    const ccomplex_t IT_0695 = IT_0693 + IT_0694;
    const ccomplex_t IT_0696 = IT_0528 + IT_0530;
    const ccomplex_t IT_0697 = 2*conj(IT_0560);
    const ccomplex_t IT_0698 = 2*conj(IT_0605);
    const ccomplex_t IT_0699 = IT_0697 + IT_0698;
    const ccomplex_t IT_0700 = conj(IT_0532) + conj(IT_0534);
    const ccomplex_t IT_0701 = 2*IT_0631;
    const ccomplex_t IT_0702 = 2*IT_0658;
    const ccomplex_t IT_0703 = IT_0701 + IT_0702;
    const ccomplex_t IT_0704 = IT_0532 + IT_0534;
    const ccomplex_t IT_0705 = 2*conj(IT_0631);
    const ccomplex_t IT_0706 = 2*conj(IT_0658);
    const ccomplex_t IT_0707 = IT_0705 + IT_0706;
    const ccomplex_t IT_0708 = (-2)*IT_0560;
    const ccomplex_t IT_0709 = (-2)*IT_0605;
    const ccomplex_t IT_0710 = IT_0708 + IT_0709;
    const ccomplex_t IT_0711 = (-2)*conj(IT_0560);
    const ccomplex_t IT_0712 = (-2)*conj(IT_0605);
    const ccomplex_t IT_0713 = IT_0711 + IT_0712;
    const ccomplex_t IT_0714 = (-2)*IT_0631;
    const ccomplex_t IT_0715 = (-2)*IT_0658;
    const ccomplex_t IT_0716 = IT_0714 + IT_0715;
    const ccomplex_t IT_0717 = (-2)*conj(IT_0631);
    const ccomplex_t IT_0718 = (-2)*conj(IT_0658);
    const ccomplex_t IT_0719 = IT_0717 + IT_0718;
    const ccomplex_t IT_0720 = (-0.25)*IT_0687;
    const ccomplex_t IT_0721 = s_24 + IT_0720;
    const ccomplex_t IT_0722 = m_N_2*IT_0721;
    const ccomplex_t IT_0723 = (-4)*IT_0722;
    const ccomplex_t IT_0724 = conj(IT_0205) + conj(IT_0236);
    const ccomplex_t IT_0725 = IT_0205 + IT_0236;
    const ccomplex_t IT_0726 = conj(IT_0042) + conj(IT_0076);
    const ccomplex_t IT_0727 = IT_0042 + IT_0076;
    const ccomplex_t IT_0728 = s_14*IT_0129;
    const ccomplex_t IT_0729 = IT_0001*IT_0728;
    const ccomplex_t IT_0730 = 0.5*IT_0729;
    const ccomplex_t IT_0731 = s_14 + IT_0730;
    const ccomplex_t IT_0732 = m_N_3*IT_0731;
    const ccomplex_t IT_0733 = 2*IT_0732;
    const ccomplex_t IT_0734 = (-0.25)*IT_0729;
    const ccomplex_t IT_0735 = s_14 + IT_0734;
    const ccomplex_t IT_0736 = m_N_3*IT_0735;
    const ccomplex_t IT_0737 = (-4)*IT_0736;
    const ccomplex_t IT_0738 = s_12*IT_0129;
    const ccomplex_t IT_0739 = s_14*s_24;
    const ccomplex_t IT_0740 = (-2)*IT_0739;
    const ccomplex_t IT_0741 = IT_0738 + IT_0740;
    const ccomplex_t IT_0742 = IT_0001*IT_0741;
    const ccomplex_t IT_0743 = (-0.5)*IT_0742;
    const ccomplex_t IT_0744 = s_12 + IT_0743;
    const ccomplex_t IT_0745 = 4*IT_0631;
    const ccomplex_t IT_0746 = 6*conj(IT_0560);
    const ccomplex_t IT_0747 = pow(m_Z, 4);
    const ccomplex_t IT_0748 = s_12*IT_0747;
    const ccomplex_t IT_0749 = IT_0001*IT_0748;
    const ccomplex_t IT_0750 = (-0.25)*IT_0749;
    const ccomplex_t IT_0751 = IT_0739 + IT_0750;
    const ccomplex_t IT_0752 = 8*conj(IT_0531);
    const ccomplex_t IT_0753 = 8*IT_0525;
    const ccomplex_t IT_0754 = 8*conj(IT_0525);
    const ccomplex_t IT_0755 = 8*IT_0529;
    const ccomplex_t IT_0756 = (-8)*IT_0236;
    const ccomplex_t IT_0757 = 8*conj(IT_0529);
    const ccomplex_t IT_0758 = (-8)*conj(IT_0236);
    const ccomplex_t IT_0759 = 8*IT_0533;
    const ccomplex_t IT_0760 = (-8)*IT_0076;
    const ccomplex_t IT_0761 = 8*conj(IT_0533);
    const ccomplex_t IT_0762 = (-8)*conj(IT_0076);
    const ccomplex_t IT_0763 = -(conj(IT_0236) + -conj(IT_0529))*IT_0531 +
       conj(IT_0205)*(IT_0528 + IT_0530 + -IT_0531) + IT_0529*conj(IT_0531) +
       conj(IT_0525)*IT_0533 + IT_0525*conj(IT_0533) + (-0.125)*IT_0236*IT_0752 
      + IT_0205*(IT_0692 + (-0.125)*IT_0752) + (-0.125)*IT_0726*IT_0753 + (
      -0.125)*IT_0727*IT_0754 + (-0.125)*IT_0692*(IT_0755 + IT_0756) + (-0.125)
      *IT_0696*(IT_0757 + IT_0758) + IT_0700*(IT_0042 + (-0.125)*IT_0759 + (
      -0.125)*IT_0760) + IT_0704*(conj(IT_0042) + (-0.125)*IT_0761 + (-0.125)
      *IT_0762);
    const ccomplex_t IT_0764 = (-0.5)*IT_0749;
    const ccomplex_t IT_0765 = IT_0738 + IT_0740 + IT_0764;
    const ccomplex_t IT_0766 = 4*conj(IT_0531);
    const ccomplex_t IT_0767 = 4*conj(IT_0529);
    const ccomplex_t IT_0768 = 4*conj(IT_0525);
    const ccomplex_t IT_0769 = 4*conj(IT_0533);
    const ccomplex_t IT_0770 = 4*IT_0525;
    const ccomplex_t IT_0771 = (-4)*IT_0532;
    const ccomplex_t IT_0772 = 4*IT_0529;
    const ccomplex_t IT_0773 = (-4)*IT_0236;
    const ccomplex_t IT_0774 = 4*IT_0531;
    const ccomplex_t IT_0775 = (-4)*IT_0530;
    const ccomplex_t IT_0776 = 4*IT_0533;
    const ccomplex_t IT_0777 = (-4)*IT_0076;
    const ccomplex_t IT_0778 = IT_0525*conj(IT_0525) + IT_0529*conj(IT_0529) +
       IT_0531*conj(IT_0531) + IT_0533*conj(IT_0533) + (-0.25)*IT_0530*IT_0766 +
       IT_0528*(IT_0692 + (-0.25)*IT_0766) + (-0.25)*IT_0236*IT_0767 + IT_0205*
      (IT_0724 + (-0.25)*IT_0767) + (-0.25)*IT_0704*IT_0768 + (-0.25)*IT_0727
      *IT_0769 + IT_0700*(IT_0534 + (-0.25)*IT_0770 + (-0.25)*IT_0771) + (-0.25)
      *IT_0724*(IT_0772 + IT_0773) + (-0.25)*IT_0692*(IT_0774 + IT_0775) +
       IT_0726*(IT_0042 + (-0.25)*IT_0776 + (-0.25)*IT_0777);
    const ccomplex_t IT_0779 = s_12*IT_0006;
    const ccomplex_t IT_0780 = s_12*IT_0677;
    const ccomplex_t IT_0781 = 2*IT_0377;
    const ccomplex_t IT_0782 = s_12*s_24;
    const ccomplex_t IT_0783 = s_14*IT_0073;
    const ccomplex_t IT_0784 = -IT_0783;
    const ccomplex_t IT_0785 = IT_0782 + IT_0784;
    const ccomplex_t IT_0786 = (-2)*IT_0514;
    const ccomplex_t IT_0787 = 2*conj(IT_0377);
    const ccomplex_t IT_0788 = (-2)*conj(IT_0514);
    const ccomplex_t IT_0789 = (conj(IT_0205) + conj(IT_0236) + conj(IT_0531))
      *IT_0664 + (IT_0205 + IT_0236 + IT_0531)*IT_0667 + (conj(IT_0528) + conj
      (IT_0529) + conj(IT_0530))*IT_0670 + (IT_0528 + IT_0529 + IT_0530)*IT_0673
       + (conj(IT_0532) + conj(IT_0533) + conj(IT_0534))*IT_0682 + (IT_0532 +
       IT_0533 + IT_0534)*IT_0685 + IT_0526*(IT_0781 + IT_0786) + IT_0527*
      (IT_0787 + IT_0788);
    const ccomplex_t IT_0790 = m_N_3*IT_0006;
    const ccomplex_t IT_0791 = m_N_3*IT_0330;
    const ccomplex_t IT_0792 = IT_0000*IT_0129;
    const ccomplex_t IT_0793 = 3*IT_0792;
    const ccomplex_t IT_0794 = 2*IT_0531;
    const ccomplex_t IT_0795 = (-2)*IT_0528;
    const ccomplex_t IT_0796 = (-2)*IT_0530;
    const ccomplex_t IT_0797 = IT_0794 + IT_0795 + IT_0796;
    const ccomplex_t IT_0798 = 2*conj(IT_0531);
    const ccomplex_t IT_0799 = (-2)*conj(IT_0528);
    const ccomplex_t IT_0800 = (-2)*conj(IT_0530);
    const ccomplex_t IT_0801 = IT_0798 + IT_0799 + IT_0800;
    const ccomplex_t IT_0802 = 2*IT_0529;
    const ccomplex_t IT_0803 = (-2)*IT_0205;
    const ccomplex_t IT_0804 = (-2)*IT_0236;
    const ccomplex_t IT_0805 = IT_0802 + IT_0803 + IT_0804;
    const ccomplex_t IT_0806 = 2*conj(IT_0529);
    const ccomplex_t IT_0807 = (-2)*conj(IT_0205);
    const ccomplex_t IT_0808 = (-2)*conj(IT_0236);
    const ccomplex_t IT_0809 = IT_0806 + IT_0807 + IT_0808;
    const ccomplex_t IT_0810 = 2*IT_0525;
    const ccomplex_t IT_0811 = 2*conj(IT_0525);
    const ccomplex_t IT_0812 = 2*IT_0533;
    const ccomplex_t IT_0813 = 2*conj(IT_0533);
    const ccomplex_t IT_0814 = (-2)*IT_0531;
    const ccomplex_t IT_0815 = (-2)*conj(IT_0531);
    const ccomplex_t IT_0816 = (-2)*IT_0529;
    const ccomplex_t IT_0817 = (-2)*conj(IT_0529);
    const ccomplex_t IT_0818 = IT_0726*IT_0797 + IT_0727*IT_0801 + IT_0700
      *IT_0805 + IT_0704*IT_0809 + IT_0724*IT_0810 + IT_0725*IT_0811 + IT_0692
      *IT_0812 + IT_0696*IT_0813 + conj(IT_0533)*IT_0814 + IT_0533*IT_0815 +
       conj(IT_0525)*IT_0816 + IT_0525*IT_0817;
    const ccomplex_t IT_0819 = (-3)*IT_0792;
    const ccomplex_t IT_0820 = IT_0700*IT_0797 + IT_0704*IT_0801 + IT_0726
      *IT_0805 + IT_0727*IT_0809 + IT_0692*IT_0810 + IT_0696*IT_0811 + IT_0724
      *IT_0812 + IT_0725*IT_0813 + conj(IT_0525)*IT_0814 + IT_0525*IT_0815 +
       conj(IT_0533)*IT_0816 + IT_0533*IT_0817;
    const ccomplex_t IT_0821 = IT_0007*(conj(IT_0077)*IT_0239 + IT_0077
      *IT_0242 + conj(IT_0321)*IT_0324 + IT_0321*IT_0327) + IT_0331*(conj
      (IT_0377)*IT_0380 + IT_0377*IT_0383 + IT_0242*IT_0433 + IT_0239*conj
      (IT_0433) + IT_0327*IT_0468 + IT_0324*conj(IT_0468) + conj(IT_0514)
      *IT_0517 + IT_0514*IT_0520) + IT_0524*IT_0535 + IT_0536*(IT_0239*IT_0606 +
       IT_0242*IT_0607 + IT_0517*IT_0659 + IT_0520*IT_0660) + IT_0661*(IT_0242
      *IT_0377 + IT_0239*conj(IT_0377) + IT_0327*IT_0514 + IT_0324*conj(IT_0514)
       + conj(IT_0321)*IT_0664 + IT_0321*IT_0667 + conj(IT_0077)*IT_0670 +
       IT_0077*IT_0673) + (conj(IT_0514)*IT_0664 + IT_0514*IT_0667 + conj
      (IT_0377)*IT_0670 + IT_0377*IT_0673)*IT_0678 + IT_0679*(IT_0606*IT_0670 +
       IT_0607*IT_0673 + IT_0659*IT_0682 + IT_0660*IT_0685) + IT_0691*(IT_0692
      *IT_0695 + IT_0696*IT_0699 + IT_0700*IT_0703 + IT_0704*IT_0707 + conj
      (IT_0531)*IT_0710 + IT_0531*IT_0713 + conj(IT_0525)*IT_0716 + IT_0525
      *IT_0719) + IT_0723*(conj(IT_0529)*IT_0710 + IT_0529*IT_0713 + conj
      (IT_0533)*IT_0716 + IT_0533*IT_0719 + IT_0695*IT_0724 + IT_0699*IT_0725 +
       IT_0703*IT_0726 + IT_0707*IT_0727) + (conj(IT_0533)*IT_0695 + IT_0533
      *IT_0699 + conj(IT_0529)*IT_0703 + IT_0529*IT_0707 + IT_0716*IT_0724 +
       IT_0719*IT_0725 + IT_0710*IT_0726 + IT_0713*IT_0727)*IT_0733 + (conj
      (IT_0525)*IT_0695 + IT_0525*IT_0699 + conj(IT_0531)*IT_0703 + IT_0531
      *IT_0707 + IT_0700*IT_0710 + IT_0704*IT_0713 + IT_0692*IT_0716 + IT_0696
      *IT_0719)*IT_0737 + 4*IT_0744*(IT_0606*IT_0607 + IT_0659*(IT_0658 + 0.25
      *IT_0745)) + 6*IT_0000*(IT_0607*IT_0659 + IT_0660*(conj(IT_0605) +
       0.166666666666667*IT_0746)) + (-8)*IT_0751*IT_0763 + (-4)*IT_0765*IT_0778
       + (conj(IT_0184)*IT_0185 + conj(IT_0237)*IT_0324 + IT_0237*IT_0326 + conj
      (IT_0077)*IT_0380 + IT_0077*IT_0382 + conj(IT_0321)*IT_0515)*IT_0779 +
       IT_0780*(conj(IT_0468)*IT_0662 + conj(IT_0433)*IT_0670 + IT_0433*IT_0672 
      + conj(IT_0514)*IT_0682 + IT_0514*IT_0684 + conj(IT_0377)*IT_0781) +
       IT_0785*IT_0789 + (IT_0520*IT_0605 + IT_0517*conj(IT_0605) + IT_0239
      *IT_0659 + IT_0242*IT_0660 + conj(IT_0321)*IT_0693 + IT_0321*IT_0697 +
       conj(IT_0077)*IT_0708 + IT_0077*IT_0711)*IT_0790 + (IT_0659*IT_0670 +
       IT_0660*IT_0673 + conj(IT_0605)*IT_0682 + IT_0605*IT_0685 + conj(IT_0514)
      *IT_0693 + IT_0514*IT_0697 + conj(IT_0377)*IT_0708 + IT_0377*IT_0711)
      *IT_0791 + IT_0793*IT_0818 + IT_0819*IT_0820;
    return create_ccomplex_return(IT_0821);
}

