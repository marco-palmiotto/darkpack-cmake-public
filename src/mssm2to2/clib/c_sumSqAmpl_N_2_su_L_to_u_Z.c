#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_su_L_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_su_L_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
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
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_ul = param->Gamma_ul;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + (-3)
      *IT_0000*IT_0003);
    const ccomplex_t IT_0005 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_u, 2);
    const ccomplex_t IT_0012 = pow(m_N_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0002;
    const ccomplex_t IT_0018 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0003;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0022;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0019 + IT_0021 + 
      -IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_su_L, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0010*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0025;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0035 = IT_0017*IT_0034;
    const ccomplex_t IT_0036 = IT_0020*IT_0034;
    const ccomplex_t IT_0037 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0038 = IT_0017*IT_0037;
    const ccomplex_t IT_0039 = IT_0020*IT_0037;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0035 + IT_0036 + 
      -IT_0038 + -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = IT_0020*IT_0048;
    const ccomplex_t IT_0051 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0052 = IT_0017*IT_0051;
    const ccomplex_t IT_0053 = IT_0020*IT_0051;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0049 + IT_0050 + 
      -IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0055*IT_0060;
    const ccomplex_t IT_0062 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = IT_0020*IT_0062;
    const ccomplex_t IT_0065 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0066 = IT_0017*IT_0065;
    const ccomplex_t IT_0067 = IT_0020*IT_0065;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0063 + IT_0064 + 
      -IT_0066 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = IT_0016 + IT_0031 + -IT_0033 + 2*IT_0047 + (-2)
      *IT_0061 + 2*IT_0075;
    const ccomplex_t IT_0077 = pow(m_Z, -2);
    const ccomplex_t IT_0078 = pow(s_14, 2);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = -IT_0012;
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = s_13*IT_0081;
    const ccomplex_t IT_0083 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0017*IT_0083;
    const ccomplex_t IT_0085 = IT_0020*IT_0083;
    const ccomplex_t IT_0086 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0017*IT_0086;
    const ccomplex_t IT_0088 = IT_0020*IT_0086;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0084 + IT_0085 + 
      -IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B1)*e_em;
    const ccomplex_t IT_0092 = IT_0001*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W1)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + 3*IT_0096);
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = IT_0058*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0090*IT_0100;
    const ccomplex_t IT_0102 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0103 = IT_0017*IT_0102;
    const ccomplex_t IT_0104 = IT_0020*IT_0102;
    const ccomplex_t IT_0105 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0106 = IT_0017*IT_0105;
    const ccomplex_t IT_0107 = IT_0020*IT_0105;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0103 + IT_0104 + 
      -IT_0106 + -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B3)*e_em;
    const ccomplex_t IT_0111 = IT_0001*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W3)*e_em;
    const ccomplex_t IT_0114 = IT_0003*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + 3*IT_0115);
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0118 = IT_0044*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0109*IT_0119;
    const ccomplex_t IT_0121 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0122 = IT_0017*IT_0121;
    const ccomplex_t IT_0123 = IT_0020*IT_0121;
    const ccomplex_t IT_0124 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0125 = IT_0017*IT_0124;
    const ccomplex_t IT_0126 = IT_0020*IT_0124;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0122 + IT_0123 + 
      -IT_0125 + -IT_0126);
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B4)*e_em;
    const ccomplex_t IT_0130 = IT_0001*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W4)*e_em;
    const ccomplex_t IT_0133 = IT_0003*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + 3*IT_0134);
    const ccomplex_t IT_0136 = 0.166666666666667*IT_0135;
    const ccomplex_t IT_0137 = IT_0072*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0128*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B2)*e_em;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W2)*e_em;
    const ccomplex_t IT_0144 = IT_0003*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + 3*IT_0145);
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0146;
    const ccomplex_t IT_0148 = IT_0013*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0005*IT_0149;
    const ccomplex_t IT_0151 = IT_0028*IT_0147;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0026*IT_0152;
    const ccomplex_t IT_0154 = IT_0032*IT_0152;
    const ccomplex_t IT_0155 = 2*IT_0101 + (-2)*IT_0120 + (-2)*IT_0139 + 
      -IT_0150 + IT_0153 + -IT_0154;
    const ccomplex_t IT_0156 = m_u*m_N_2;
    const ccomplex_t IT_0157 = IT_0081*IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0031 + 0.5*IT_0033 + -IT_0047 +
       IT_0061 + -IT_0075;
    const ccomplex_t IT_0159 = s_13*s_14;
    const ccomplex_t IT_0160 = s_34*IT_0012;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = m_u*IT_0081;
    const ccomplex_t IT_0164 = e_em*IT_0017;
    const ccomplex_t IT_0165 = e_em*IT_0020;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0164 + (-3)*IT_0165);
    const ccomplex_t IT_0167 = (-0.166666666666667)*IT_0166;
    const ccomplex_t IT_0168 = IT_0010*IT_0167;
    const ccomplex_t IT_0169 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0171 = IT_0168*IT_0169*IT_0170;
    const ccomplex_t IT_0172 = m_N_2*IT_0090;
    const ccomplex_t IT_0173 = IT_0100*IT_0172;
    const ccomplex_t IT_0174 = m_N_4*IT_0069;
    const ccomplex_t IT_0175 = m_N_2*IT_0128;
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0174 + IT_0176;
    const ccomplex_t IT_0178 = IT_0138*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0180 = (-0.666666666666667)*IT_0179;
    const ccomplex_t IT_0181 = IT_0010*IT_0180;
    const ccomplex_t IT_0182 = IT_0169*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = m_u*IT_0183;
    const ccomplex_t IT_0185 = m_N_2*IT_0026;
    const ccomplex_t IT_0186 = IT_0152*IT_0185;
    const ccomplex_t IT_0187 = m_N_2*IT_0032;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0185 + IT_0188;
    const ccomplex_t IT_0190 = IT_0152*IT_0189;
    const ccomplex_t IT_0191 = m_N_3*IT_0041;
    const ccomplex_t IT_0192 = m_N_2*IT_0109;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0119*IT_0194;
    const ccomplex_t IT_0196 = -IT_0171 + -IT_0173 + -IT_0178 + IT_0184 + (
      -0.5)*IT_0186 + (-0.5)*IT_0190 + -IT_0195;
    const ccomplex_t IT_0197 = m_N_1*IT_0055;
    const ccomplex_t IT_0198 = IT_0100*IT_0197;
    const ccomplex_t IT_0199 = IT_0152*IT_0187;
    const ccomplex_t IT_0200 = IT_0198 + 0.5*IT_0199;
    const ccomplex_t IT_0201 = IT_0196 + IT_0200;
    const ccomplex_t IT_0202 = IT_0147*IT_0167;
    const ccomplex_t IT_0203 = IT_0169*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = m_u*IT_0204;
    const ccomplex_t IT_0206 = IT_0030*IT_0187;
    const ccomplex_t IT_0207 = IT_0030*IT_0189;
    const ccomplex_t IT_0208 = IT_0147*IT_0180;
    const ccomplex_t IT_0209 = IT_0169*IT_0170*IT_0208;
    const ccomplex_t IT_0210 = m_N_2*IT_0055;
    const ccomplex_t IT_0211 = IT_0060*IT_0210;
    const ccomplex_t IT_0212 = m_N_3*IT_0109;
    const ccomplex_t IT_0213 = m_N_2*IT_0041;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = IT_0046*IT_0215;
    const ccomplex_t IT_0217 = m_N_4*IT_0128;
    const ccomplex_t IT_0218 = m_N_2*IT_0069;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = IT_0074*IT_0220;
    const ccomplex_t IT_0222 = -IT_0205 + 0.5*IT_0206 + (-0.5)*IT_0207 +
       IT_0209 + IT_0211 + IT_0216 + IT_0221;
    const ccomplex_t IT_0223 = IT_0030*IT_0185;
    const ccomplex_t IT_0224 = m_N_1*IT_0090;
    const ccomplex_t IT_0225 = IT_0060*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*IT_0223 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0222 + IT_0226;
    const ccomplex_t IT_0228 = (-18)*conj(IT_0196);
    const ccomplex_t IT_0229 = -IT_0016 + (-2)*IT_0183;
    const ccomplex_t IT_0230 = s_14*s_34*IT_0077;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = s_13 + IT_0231;
    const ccomplex_t IT_0233 = s_13*IT_0232;
    const ccomplex_t IT_0234 = IT_0150 + 2*IT_0204;
    const ccomplex_t IT_0235 = IT_0156*IT_0232;
    const ccomplex_t IT_0236 = -IT_0183;
    const ccomplex_t IT_0237 = m_N_2*IT_0232;
    const ccomplex_t IT_0238 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0239 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0240 = pow(m_Z, 2);
    const ccomplex_t IT_0241 = s_14*IT_0240;
    const ccomplex_t IT_0242 = IT_0077*IT_0241;
    const ccomplex_t IT_0243 = (-0.25)*IT_0242;
    const ccomplex_t IT_0244 = s_14 + IT_0243;
    const ccomplex_t IT_0245 = m_u*IT_0244;
    const ccomplex_t IT_0246 = (-4)*IT_0245;
    const ccomplex_t IT_0247 = 0.5*IT_0242;
    const ccomplex_t IT_0248 = s_14 + IT_0247;
    const ccomplex_t IT_0249 = m_u*IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0249;
    const ccomplex_t IT_0251 = -IT_0101 + IT_0120 + IT_0139 + (-0.5)*IT_0153 +
       0.5*IT_0154;
    const ccomplex_t IT_0252 = s_34*IT_0240;
    const ccomplex_t IT_0253 = IT_0077*IT_0252;
    const ccomplex_t IT_0254 = (-0.25)*IT_0253;
    const ccomplex_t IT_0255 = s_34 + IT_0254;
    const ccomplex_t IT_0256 = m_N_2*IT_0255;
    const ccomplex_t IT_0257 = (-4)*IT_0256;
    const ccomplex_t IT_0258 = 0.5*IT_0253;
    const ccomplex_t IT_0259 = s_34 + IT_0258;
    const ccomplex_t IT_0260 = m_N_2*IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0260;
    const ccomplex_t IT_0262 = s_13*IT_0240;
    const ccomplex_t IT_0263 = s_14*s_34;
    const ccomplex_t IT_0264 = (-2)*IT_0263;
    const ccomplex_t IT_0265 = IT_0262 + IT_0264;
    const ccomplex_t IT_0266 = IT_0077*IT_0265;
    const ccomplex_t IT_0267 = (-6)*IT_0266;
    const ccomplex_t IT_0268 = 12*s_13;
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = m_u*IT_0232;
    const ccomplex_t IT_0271 = (-6)*IT_0229;
    const ccomplex_t IT_0272 = pow(s_34, 2);
    const ccomplex_t IT_0273 = IT_0077*IT_0272;
    const ccomplex_t IT_0274 = -IT_0273;
    const ccomplex_t IT_0275 = IT_0011 + IT_0274;
    const ccomplex_t IT_0276 = IT_0156*IT_0275;
    const ccomplex_t IT_0277 = s_14*IT_0011;
    const ccomplex_t IT_0278 = s_13*s_34;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = IT_0277 + IT_0279;
    const ccomplex_t IT_0281 = m_N_2*IT_0275;
    const ccomplex_t IT_0282 = (-6)*conj(IT_0229);
    const ccomplex_t IT_0283 = (-6)*conj(IT_0234);
    const ccomplex_t IT_0284 = (-6)*conj(IT_0196);
    const ccomplex_t IT_0285 = (-6)*IT_0076;
    const ccomplex_t IT_0286 = (-6)*IT_0234;
    const ccomplex_t IT_0287 = 6*IT_0236;
    const ccomplex_t IT_0288 = 6*IT_0158;
    const ccomplex_t IT_0289 = 6*IT_0251;
    const ccomplex_t IT_0290 = 6*IT_0204;
    const ccomplex_t IT_0291 = 6*IT_0155*IT_0163 + IT_0156*((-18)*IT_0196 + (
      -18)*IT_0200) + IT_0227*IT_0269 + IT_0271*IT_0281 + IT_0237*IT_0285 +
       IT_0270*IT_0286 + IT_0261*IT_0287 + IT_0257*IT_0288 + IT_0250*IT_0289 +
       IT_0246*IT_0290;
    const ccomplex_t IT_0292 = (-6)*IT_0155;
    const ccomplex_t IT_0293 = 6*IT_0076*IT_0163 + IT_0156*((-18)*IT_0222 + (
      -18)*IT_0226) + IT_0201*IT_0269 + IT_0270*IT_0271 + IT_0281*IT_0286 +
       IT_0246*IT_0287 + IT_0250*IT_0288 + IT_0257*IT_0289 + IT_0261*IT_0290 +
       IT_0237*IT_0292;
    const ccomplex_t IT_0294 = IT_0156*IT_0240;
    const ccomplex_t IT_0295 = (-3)*IT_0294;
    const ccomplex_t IT_0296 = 3*IT_0294;
    const ccomplex_t IT_0297 = pow(m_Z, 4);
    const ccomplex_t IT_0298 = s_13*IT_0297;
    const ccomplex_t IT_0299 = IT_0077*IT_0298;
    const ccomplex_t IT_0300 = (-24)*IT_0263;
    const ccomplex_t IT_0301 = 6*IT_0299 + IT_0300;
    const ccomplex_t IT_0302 = 0.166666666666667*IT_0236;
    const ccomplex_t IT_0303 = 24*IT_0263;
    const ccomplex_t IT_0304 = (-12)*IT_0262;
    const ccomplex_t IT_0305 = 6*IT_0299 + IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0307 = (-6)*IT_0196;
    const ccomplex_t IT_0308 = 0.166666666666667*IT_0235;
    const ccomplex_t IT_0309 = 0.166666666666667*IT_0251;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0204;
    const ccomplex_t IT_0311 = 0.166666666666667*IT_0233;
    const ccomplex_t IT_0312 = s_13*IT_0275;
    const ccomplex_t IT_0313 = 0.166666666666667*IT_0312;
    const ccomplex_t IT_0314 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0280;
    const ccomplex_t IT_0316 = 6*conj(IT_0076)*(IT_0076*IT_0082 + IT_0155
      *IT_0157 + IT_0158*IT_0162 + IT_0163*IT_0201) + IT_0156*IT_0227*IT_0228 + 
      (IT_0229*IT_0233 + IT_0234*IT_0235 + IT_0162*IT_0236 + IT_0227*IT_0237)
      *IT_0238 + (IT_0162*IT_0204 + IT_0233*IT_0234 + IT_0200*IT_0237)*IT_0239 +
       6*conj(IT_0196)*(IT_0076*IT_0163 + IT_0236*IT_0246 + IT_0158*IT_0250 +
       IT_0251*IT_0257 + IT_0204*IT_0261 + 0.166666666666667*IT_0201*IT_0269 +
       0.166666666666667*IT_0270*IT_0271) + (IT_0155*IT_0235 + IT_0201*IT_0270 +
       IT_0234*IT_0276 + IT_0236*IT_0280 + IT_0227*IT_0281)*IT_0282 + (IT_0227
      *IT_0270 + IT_0204*IT_0280 + IT_0200*IT_0281)*IT_0283 + (IT_0155*IT_0237 +
       IT_0234*IT_0281)*IT_0284 + (conj(IT_0222) + conj(IT_0226))*IT_0291 + conj
      (IT_0200)*IT_0293 + 6*conj(IT_0158)*(IT_0076*IT_0162 + IT_0201*IT_0250 +
       IT_0227*IT_0257 + IT_0229*IT_0280 + IT_0251*IT_0295 + IT_0204*IT_0296 +
       IT_0301*IT_0302 + IT_0305*IT_0306) + 6*conj(IT_0155)*(IT_0082*IT_0155 +
       IT_0076*IT_0157 + IT_0163*IT_0227 + IT_0162*IT_0251 + 0.166666666666667
      *IT_0237*IT_0307 + IT_0271*IT_0308) + 6*conj(IT_0251)*(IT_0155*IT_0162 +
       IT_0227*IT_0250 + IT_0201*IT_0257 + IT_0234*IT_0280 + IT_0158*IT_0295 +
       IT_0236*IT_0296 + IT_0305*IT_0309 + IT_0301*IT_0310) + 6*conj(IT_0229)*
      (IT_0158*IT_0280 + IT_0285*IT_0311 + IT_0271*IT_0313) + 6*conj(IT_0234)*
      (0.166666666666667*IT_0271*IT_0276 + IT_0251*IT_0280 + 0.166666666666667
      *IT_0281*IT_0307 + IT_0285*IT_0308 + IT_0292*IT_0311 + IT_0286*IT_0313) +
       6*conj(IT_0236)*(IT_0201*IT_0246 + IT_0227*IT_0261 + IT_0204*IT_0295 +
       IT_0251*IT_0296 + IT_0302*IT_0305 + IT_0301*IT_0306 + IT_0285*IT_0314 +
       IT_0271*IT_0315) + 6*conj(IT_0204)*(IT_0227*IT_0246 + IT_0201*IT_0261 +
       IT_0236*IT_0295 + IT_0158*IT_0296 + IT_0301*IT_0309 + IT_0305*IT_0310 +
       IT_0292*IT_0314 + IT_0286*IT_0315);
    return create_ccomplex_return(IT_0316);
}

