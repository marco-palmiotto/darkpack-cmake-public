#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_su_L_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_su_L_to_u_Z(
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
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_u, 2);
    const ccomplex_t IT_0020 = pow(m_su_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0013*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + (-3)
      *IT_0000*IT_0006);
    const ccomplex_t IT_0026 = 0.333333333333333*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_N_1, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0019 + IT_0029 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0032;
    const ccomplex_t IT_0034 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0035 = IT_0003*IT_0034;
    const ccomplex_t IT_0036 = IT_0007*IT_0034;
    const ccomplex_t IT_0037 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0037;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0035 + IT_0036 + 
      -IT_0038 + -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = IT_0007*IT_0048;
    const ccomplex_t IT_0051 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0052 = IT_0003*IT_0051;
    const ccomplex_t IT_0053 = IT_0007*IT_0051;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0049 + IT_0050 + 
      -IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0006*IT_0014*IT_0016;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0055*IT_0060;
    const ccomplex_t IT_0062 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = IT_0007*IT_0062;
    const ccomplex_t IT_0065 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = IT_0007*IT_0065;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0063 + IT_0064 + 
      -IT_0066 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0028*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0069*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0068;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = 2*IT_0024 + IT_0033 + 2*IT_0047 + 2*IT_0061 +
       IT_0073 + -IT_0075;
    const ccomplex_t IT_0077 = pow(m_Z, -2);
    const ccomplex_t IT_0078 = pow(s_14, 2);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = -IT_0029;
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = s_13*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B1)*e_em;
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W1)*e_em;
    const ccomplex_t IT_0087 = IT_0006*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0085 + 3*IT_0088);
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = IT_0070*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0074*IT_0092;
    const ccomplex_t IT_0094 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = IT_0007*IT_0094;
    const ccomplex_t IT_0097 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = IT_0007*IT_0097;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0095 + IT_0096 + 
      -IT_0098 + -IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = conj(N_B3)*e_em;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W3)*e_em;
    const ccomplex_t IT_0106 = IT_0006*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + 3*IT_0107);
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0110 = IT_0021*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0101*IT_0111;
    const ccomplex_t IT_0113 = IT_0030*IT_0090;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0026*IT_0114;
    const ccomplex_t IT_0116 = IT_0069*IT_0092;
    const ccomplex_t IT_0117 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0118 = IT_0003*IT_0117;
    const ccomplex_t IT_0119 = IT_0007*IT_0117;
    const ccomplex_t IT_0120 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0121 = IT_0003*IT_0120;
    const ccomplex_t IT_0122 = IT_0007*IT_0120;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0118 + IT_0119 + 
      -IT_0121 + -IT_0122);
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B2)*e_em;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W2)*e_em;
    const ccomplex_t IT_0129 = IT_0006*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 3*IT_0130);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0131;
    const ccomplex_t IT_0133 = IT_0044*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0124*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0137 = IT_0003*IT_0136;
    const ccomplex_t IT_0138 = IT_0007*IT_0136;
    const ccomplex_t IT_0139 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0140 = IT_0003*IT_0139;
    const ccomplex_t IT_0141 = IT_0007*IT_0139;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0137 + IT_0138 + 
      -IT_0140 + -IT_0141);
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = conj(N_B4)*e_em;
    const ccomplex_t IT_0145 = IT_0001*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = conj(N_W4)*e_em;
    const ccomplex_t IT_0148 = IT_0006*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 3*IT_0149);
    const ccomplex_t IT_0151 = 0.166666666666667*IT_0150;
    const ccomplex_t IT_0152 = IT_0058*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0143*IT_0153;
    const ccomplex_t IT_0155 = -IT_0093 + (-2)*IT_0112 + -IT_0115 + IT_0116 + 
      (-2)*IT_0135 + (-2)*IT_0154;
    const ccomplex_t IT_0156 = m_u*m_N_1;
    const ccomplex_t IT_0157 = IT_0081*IT_0156;
    const ccomplex_t IT_0158 = -IT_0024 + -IT_0047 + -IT_0061 + (-0.5)*IT_0073
       + 0.5*IT_0075;
    const ccomplex_t IT_0159 = s_13*s_14;
    const ccomplex_t IT_0160 = s_34*IT_0029;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = m_u*IT_0081;
    const ccomplex_t IT_0164 = m_N_1*IT_0074;
    const ccomplex_t IT_0165 = IT_0092*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = m_N_2*IT_0041;
    const ccomplex_t IT_0168 = m_N_1*IT_0124;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = IT_0134*IT_0170;
    const ccomplex_t IT_0172 = cpow((-2)*s_12 + IT_0019 + -IT_0020 + -IT_0029 
      + -reg_prop, -1);
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0174 = e_em*IT_0003;
    const ccomplex_t IT_0175 = e_em*IT_0007;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0174 + (-3)*IT_0175);
    const ccomplex_t IT_0177 = (-0.166666666666667)*IT_0176;
    const ccomplex_t IT_0178 = IT_0028*IT_0177;
    const ccomplex_t IT_0179 = IT_0172*IT_0173*IT_0178;
    const ccomplex_t IT_0180 = m_N_1*IT_0069;
    const ccomplex_t IT_0181 = IT_0092*IT_0180;
    const ccomplex_t IT_0182 = -IT_0164;
    const ccomplex_t IT_0183 = IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = IT_0092*IT_0183;
    const ccomplex_t IT_0185 = m_N_3*IT_0013;
    const ccomplex_t IT_0186 = m_N_1*IT_0101;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = IT_0111*IT_0188;
    const ccomplex_t IT_0190 = m_N_4*IT_0055;
    const ccomplex_t IT_0191 = m_N_1*IT_0143;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = IT_0153*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0196 = (-0.666666666666667)*IT_0195;
    const ccomplex_t IT_0197 = IT_0028*IT_0196;
    const ccomplex_t IT_0198 = IT_0172*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = m_u*IT_0199;
    const ccomplex_t IT_0201 = -IT_0171 + -IT_0179 + (-0.5)*IT_0181 + (-0.5)
      *IT_0184 + -IT_0189 + -IT_0194 + IT_0200;
    const ccomplex_t IT_0202 = IT_0166 + IT_0201;
    const ccomplex_t IT_0203 = IT_0072*IT_0180;
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = m_N_3*IT_0101;
    const ccomplex_t IT_0206 = m_N_1*IT_0013;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = IT_0023*IT_0208;
    const ccomplex_t IT_0210 = m_N_4*IT_0143;
    const ccomplex_t IT_0211 = m_N_1*IT_0055;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = IT_0060*IT_0213;
    const ccomplex_t IT_0215 = IT_0090*IT_0196;
    const ccomplex_t IT_0216 = IT_0172*IT_0173*IT_0215;
    const ccomplex_t IT_0217 = IT_0090*IT_0177;
    const ccomplex_t IT_0218 = IT_0172*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_u*IT_0219;
    const ccomplex_t IT_0221 = IT_0072*IT_0164;
    const ccomplex_t IT_0222 = IT_0072*IT_0183;
    const ccomplex_t IT_0223 = m_N_2*IT_0124;
    const ccomplex_t IT_0224 = m_N_1*IT_0041;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = IT_0046*IT_0226;
    const ccomplex_t IT_0228 = IT_0209 + IT_0214 + IT_0216 + -IT_0220 + 0.5
      *IT_0221 + (-0.5)*IT_0222 + IT_0227;
    const ccomplex_t IT_0229 = IT_0204 + IT_0228;
    const ccomplex_t IT_0230 = (-18)*conj(IT_0201);
    const ccomplex_t IT_0231 = -IT_0033 + (-2)*IT_0199;
    const ccomplex_t IT_0232 = s_14*s_34*IT_0077;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = s_13 + IT_0233;
    const ccomplex_t IT_0235 = s_13*IT_0234;
    const ccomplex_t IT_0236 = IT_0115 + 2*IT_0219;
    const ccomplex_t IT_0237 = IT_0156*IT_0234;
    const ccomplex_t IT_0238 = -IT_0199;
    const ccomplex_t IT_0239 = m_N_1*IT_0234;
    const ccomplex_t IT_0240 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0242 = pow(m_Z, 2);
    const ccomplex_t IT_0243 = s_14*IT_0242;
    const ccomplex_t IT_0244 = IT_0077*IT_0243;
    const ccomplex_t IT_0245 = (-0.25)*IT_0244;
    const ccomplex_t IT_0246 = s_14 + IT_0245;
    const ccomplex_t IT_0247 = m_u*IT_0246;
    const ccomplex_t IT_0248 = (-4)*IT_0247;
    const ccomplex_t IT_0249 = 0.5*IT_0244;
    const ccomplex_t IT_0250 = s_14 + IT_0249;
    const ccomplex_t IT_0251 = m_u*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0251;
    const ccomplex_t IT_0253 = 0.5*IT_0093 + IT_0112 + (-0.5)*IT_0116 +
       IT_0135 + IT_0154;
    const ccomplex_t IT_0254 = s_34*IT_0242;
    const ccomplex_t IT_0255 = IT_0077*IT_0254;
    const ccomplex_t IT_0256 = (-0.25)*IT_0255;
    const ccomplex_t IT_0257 = s_34 + IT_0256;
    const ccomplex_t IT_0258 = m_N_1*IT_0257;
    const ccomplex_t IT_0259 = (-4)*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0255;
    const ccomplex_t IT_0261 = s_34 + IT_0260;
    const ccomplex_t IT_0262 = m_N_1*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0262;
    const ccomplex_t IT_0264 = s_13*IT_0242;
    const ccomplex_t IT_0265 = s_14*s_34;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = IT_0264 + IT_0266;
    const ccomplex_t IT_0268 = IT_0077*IT_0267;
    const ccomplex_t IT_0269 = (-6)*IT_0268;
    const ccomplex_t IT_0270 = 12*s_13;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = m_u*IT_0234;
    const ccomplex_t IT_0273 = (-6)*IT_0231;
    const ccomplex_t IT_0274 = pow(s_34, 2);
    const ccomplex_t IT_0275 = IT_0077*IT_0274;
    const ccomplex_t IT_0276 = -IT_0275;
    const ccomplex_t IT_0277 = IT_0019 + IT_0276;
    const ccomplex_t IT_0278 = IT_0156*IT_0277;
    const ccomplex_t IT_0279 = s_14*IT_0019;
    const ccomplex_t IT_0280 = s_13*s_34;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = IT_0279 + IT_0281;
    const ccomplex_t IT_0283 = m_N_1*IT_0277;
    const ccomplex_t IT_0284 = (-6)*conj(IT_0231);
    const ccomplex_t IT_0285 = (-6)*conj(IT_0236);
    const ccomplex_t IT_0286 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0287 = (-6)*IT_0076;
    const ccomplex_t IT_0288 = (-6)*IT_0236;
    const ccomplex_t IT_0289 = 6*IT_0158;
    const ccomplex_t IT_0290 = 6*IT_0253;
    const ccomplex_t IT_0291 = 6*IT_0219;
    const ccomplex_t IT_0292 = 6*IT_0155*IT_0163 + IT_0156*((-18)*IT_0166 + (
      -18)*IT_0201) + 6*IT_0238*IT_0263 + IT_0229*IT_0271 + IT_0273*IT_0283 +
       IT_0239*IT_0287 + IT_0272*IT_0288 + IT_0259*IT_0289 + IT_0252*IT_0290 +
       IT_0248*IT_0291;
    const ccomplex_t IT_0293 = (-6)*IT_0155;
    const ccomplex_t IT_0294 = 6*IT_0076*IT_0163 + IT_0156*((-18)*IT_0204 + (
      -18)*IT_0228) + 6*IT_0238*IT_0248 + IT_0202*IT_0271 + IT_0272*IT_0273 +
       IT_0283*IT_0288 + IT_0252*IT_0289 + IT_0259*IT_0290 + IT_0263*IT_0291 +
       IT_0239*IT_0293;
    const ccomplex_t IT_0295 = IT_0156*IT_0242;
    const ccomplex_t IT_0296 = 3*IT_0295;
    const ccomplex_t IT_0297 = (-3)*IT_0295;
    const ccomplex_t IT_0298 = pow(m_Z, 4);
    const ccomplex_t IT_0299 = s_13*IT_0298;
    const ccomplex_t IT_0300 = IT_0077*IT_0299;
    const ccomplex_t IT_0301 = (-24)*IT_0265;
    const ccomplex_t IT_0302 = 6*IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = 24*IT_0265;
    const ccomplex_t IT_0304 = (-12)*IT_0264;
    const ccomplex_t IT_0305 = 6*IT_0300 + IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0253;
    const ccomplex_t IT_0308 = 0.166666666666667*IT_0219;
    const ccomplex_t IT_0309 = (-6)*IT_0201;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0237;
    const ccomplex_t IT_0311 = 0.166666666666667*IT_0235;
    const ccomplex_t IT_0312 = s_13*IT_0277;
    const ccomplex_t IT_0313 = 0.166666666666667*IT_0312;
    const ccomplex_t IT_0314 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0282;
    const ccomplex_t IT_0316 = 6*conj(IT_0076)*(IT_0076*IT_0082 + IT_0155
      *IT_0157 + IT_0158*IT_0162 + IT_0163*IT_0202) + IT_0156*IT_0229*IT_0230 + 
      (IT_0231*IT_0235 + IT_0236*IT_0237 + IT_0162*IT_0238 + IT_0229*IT_0239)
      *IT_0240 + (IT_0162*IT_0219 + IT_0235*IT_0236 + IT_0166*IT_0239)*IT_0241 +
       6*conj(IT_0201)*(IT_0076*IT_0163 + IT_0238*IT_0248 + IT_0158*IT_0252 +
       IT_0253*IT_0259 + IT_0219*IT_0263 + 0.166666666666667*IT_0202*IT_0271 +
       0.166666666666667*IT_0272*IT_0273) + (IT_0155*IT_0237 + IT_0202*IT_0272 +
       IT_0236*IT_0278 + IT_0238*IT_0282 + IT_0229*IT_0283)*IT_0284 + (IT_0229
      *IT_0272 + IT_0219*IT_0282 + IT_0166*IT_0283)*IT_0285 + (IT_0155*IT_0239 +
       IT_0236*IT_0283)*IT_0286 + (conj(IT_0204) + conj(IT_0228))*IT_0292 + conj
      (IT_0166)*IT_0294 + 6*conj(IT_0158)*(IT_0076*IT_0162 + IT_0202*IT_0252 +
       IT_0229*IT_0259 + IT_0231*IT_0282 + IT_0219*IT_0296 + IT_0253*IT_0297 +
       0.166666666666667*IT_0238*IT_0302 + IT_0305*IT_0306) + 6*conj(IT_0253)*
      (IT_0155*IT_0162 + IT_0229*IT_0252 + IT_0202*IT_0259 + IT_0236*IT_0282 +
       IT_0238*IT_0296 + IT_0158*IT_0297 + IT_0305*IT_0307 + IT_0302*IT_0308) +
       6*conj(IT_0155)*(IT_0082*IT_0155 + IT_0076*IT_0157 + IT_0163*IT_0229 +
       IT_0162*IT_0253 + 0.166666666666667*IT_0239*IT_0309 + IT_0273*IT_0310) +
       6*conj(IT_0231)*(IT_0158*IT_0282 + IT_0287*IT_0311 + IT_0273*IT_0313) + 6
      *conj(IT_0236)*(0.166666666666667*IT_0273*IT_0278 + IT_0253*IT_0282 +
       0.166666666666667*IT_0283*IT_0309 + IT_0287*IT_0310 + IT_0293*IT_0311 +
       IT_0288*IT_0313) + 6*conj(IT_0238)*(IT_0202*IT_0248 + IT_0229*IT_0263 +
       IT_0253*IT_0296 + IT_0219*IT_0297 + 0.166666666666667*IT_0238*IT_0305 +
       IT_0302*IT_0306 + IT_0287*IT_0314 + IT_0273*IT_0315) + 6*conj(IT_0219)*
      (IT_0229*IT_0248 + IT_0202*IT_0263 + IT_0158*IT_0296 + IT_0238*IT_0297 +
       IT_0302*IT_0307 + IT_0305*IT_0308 + IT_0293*IT_0314 + IT_0288*IT_0315);
    return create_ccomplex_return(IT_0316);
}

