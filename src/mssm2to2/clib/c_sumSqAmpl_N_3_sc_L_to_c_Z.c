#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sc_L_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sc_L_to_c_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
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
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_c, 2);
    const ccomplex_t IT_0012 = pow(m_N_3, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0002;
    const ccomplex_t IT_0018 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0003;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0022;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0019 + IT_0021 + 
      -IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_sc_L, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0010*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0020*IT_0032;
    const ccomplex_t IT_0035 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0036 = IT_0017*IT_0035;
    const ccomplex_t IT_0037 = IT_0020*IT_0035;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0033 + IT_0034 + 
      -IT_0036 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0047 = IT_0017*IT_0046;
    const ccomplex_t IT_0048 = IT_0020*IT_0046;
    const ccomplex_t IT_0049 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0050 = IT_0017*IT_0049;
    const ccomplex_t IT_0051 = IT_0020*IT_0049;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0047 + IT_0048 + 
      -IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0053*IT_0058;
    const ccomplex_t IT_0060 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = IT_0020*IT_0060;
    const ccomplex_t IT_0063 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0064 = IT_0017*IT_0063;
    const ccomplex_t IT_0065 = IT_0020*IT_0063;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0061 + IT_0062 + 
      -IT_0064 + -IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0011 + IT_0027 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0067*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0025;
    const ccomplex_t IT_0075 = IT_0030*IT_0074;
    const ccomplex_t IT_0076 = IT_0016 + IT_0031 + (-2)*IT_0045 + 2*IT_0059 + 
      (-2)*IT_0073 + -IT_0075;
    const ccomplex_t IT_0077 = pow(m_Z, -2);
    const ccomplex_t IT_0078 = pow(s_14, 2);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = -IT_0012;
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = s_13*IT_0081;
    const ccomplex_t IT_0083 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0084 = IT_0017*IT_0083;
    const ccomplex_t IT_0085 = IT_0020*IT_0083;
    const ccomplex_t IT_0086 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0087 = IT_0017*IT_0086;
    const ccomplex_t IT_0088 = IT_0020*IT_0086;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0084 + IT_0085 + 
      -IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B4)*e_em;
    const ccomplex_t IT_0092 = IT_0001*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W4)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + 3*IT_0096);
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = IT_0056*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0090*IT_0100;
    const ccomplex_t IT_0102 = conj(N_B3)*e_em;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W3)*e_em;
    const ccomplex_t IT_0106 = IT_0003*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + 3*IT_0107);
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0110 = IT_0013*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0005*IT_0111;
    const ccomplex_t IT_0113 = IT_0028*IT_0109;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0026*IT_0114;
    const ccomplex_t IT_0116 = IT_0074*IT_0114;
    const ccomplex_t IT_0117 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0118 = IT_0017*IT_0117;
    const ccomplex_t IT_0119 = IT_0020*IT_0117;
    const ccomplex_t IT_0120 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = IT_0020*IT_0120;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0118 + IT_0119 + 
      -IT_0121 + -IT_0122);
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B1)*e_em;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*e_em;
    const ccomplex_t IT_0129 = IT_0003*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 3*IT_0130);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0131;
    const ccomplex_t IT_0133 = IT_0070*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0124*IT_0134;
    const ccomplex_t IT_0136 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0137 = IT_0017*IT_0136;
    const ccomplex_t IT_0138 = IT_0020*IT_0136;
    const ccomplex_t IT_0139 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0140 = IT_0017*IT_0139;
    const ccomplex_t IT_0141 = IT_0020*IT_0139;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0137 + IT_0138 + 
      -IT_0140 + -IT_0141);
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = conj(N_B2)*e_em;
    const ccomplex_t IT_0145 = IT_0001*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = conj(N_W2)*e_em;
    const ccomplex_t IT_0148 = IT_0003*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 3*IT_0149);
    const ccomplex_t IT_0151 = 0.166666666666667*IT_0150;
    const ccomplex_t IT_0152 = IT_0042*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0143*IT_0153;
    const ccomplex_t IT_0155 = (-2)*IT_0101 + -IT_0112 + IT_0115 + -IT_0116 +
       2*IT_0135 + 2*IT_0154;
    const ccomplex_t IT_0156 = m_c*m_N_3;
    const ccomplex_t IT_0157 = IT_0081*IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0031 + IT_0045 + -IT_0059 + IT_0073 +
       0.5*IT_0075;
    const ccomplex_t IT_0159 = s_13*s_14;
    const ccomplex_t IT_0160 = s_34*IT_0012;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = m_c*IT_0081;
    const ccomplex_t IT_0164 = e_em*IT_0017;
    const ccomplex_t IT_0165 = e_em*IT_0020;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0164 + (-3)*IT_0165);
    const ccomplex_t IT_0167 = (-0.166666666666667)*IT_0166;
    const ccomplex_t IT_0168 = IT_0010*IT_0167;
    const ccomplex_t IT_0169 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0171 = IT_0168*IT_0169*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0173 = (-0.666666666666667)*IT_0172;
    const ccomplex_t IT_0174 = IT_0010*IT_0173;
    const ccomplex_t IT_0175 = IT_0169*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = m_c*IT_0176;
    const ccomplex_t IT_0178 = m_N_3*IT_0026;
    const ccomplex_t IT_0179 = IT_0114*IT_0178;
    const ccomplex_t IT_0180 = m_N_3*IT_0074;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = IT_0178 + IT_0181;
    const ccomplex_t IT_0183 = IT_0114*IT_0182;
    const ccomplex_t IT_0184 = m_N_3*IT_0124;
    const ccomplex_t IT_0185 = IT_0134*IT_0184;
    const ccomplex_t IT_0186 = m_N_3*IT_0143;
    const ccomplex_t IT_0187 = IT_0153*IT_0186;
    const ccomplex_t IT_0188 = m_N_4*IT_0053;
    const ccomplex_t IT_0189 = m_N_3*IT_0090;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = IT_0100*IT_0191;
    const ccomplex_t IT_0193 = -IT_0171 + IT_0177 + (-0.5)*IT_0179 + (-0.5)
      *IT_0183 + -IT_0185 + -IT_0187 + -IT_0192;
    const ccomplex_t IT_0194 = IT_0114*IT_0180;
    const ccomplex_t IT_0195 = m_N_1*IT_0067;
    const ccomplex_t IT_0196 = IT_0134*IT_0195;
    const ccomplex_t IT_0197 = m_N_2*IT_0039;
    const ccomplex_t IT_0198 = IT_0153*IT_0197;
    const ccomplex_t IT_0199 = 0.5*IT_0194 + IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0193 + IT_0199;
    const ccomplex_t IT_0201 = IT_0109*IT_0173;
    const ccomplex_t IT_0202 = IT_0169*IT_0170*IT_0201;
    const ccomplex_t IT_0203 = IT_0109*IT_0167;
    const ccomplex_t IT_0204 = IT_0169*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_c*IT_0205;
    const ccomplex_t IT_0207 = IT_0030*IT_0180;
    const ccomplex_t IT_0208 = IT_0030*IT_0182;
    const ccomplex_t IT_0209 = m_N_3*IT_0067;
    const ccomplex_t IT_0210 = IT_0072*IT_0209;
    const ccomplex_t IT_0211 = m_N_3*IT_0039;
    const ccomplex_t IT_0212 = IT_0044*IT_0211;
    const ccomplex_t IT_0213 = m_N_4*IT_0090;
    const ccomplex_t IT_0214 = m_N_3*IT_0053;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = IT_0058*IT_0216;
    const ccomplex_t IT_0218 = IT_0202 + -IT_0206 + 0.5*IT_0207 + (-0.5)
      *IT_0208 + IT_0210 + IT_0212 + IT_0217;
    const ccomplex_t IT_0219 = IT_0030*IT_0178;
    const ccomplex_t IT_0220 = m_N_1*IT_0124;
    const ccomplex_t IT_0221 = IT_0072*IT_0220;
    const ccomplex_t IT_0222 = m_N_2*IT_0143;
    const ccomplex_t IT_0223 = IT_0044*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0219 + -IT_0221 + -IT_0223;
    const ccomplex_t IT_0225 = IT_0218 + IT_0224;
    const ccomplex_t IT_0226 = (-18)*conj(IT_0193);
    const ccomplex_t IT_0227 = -IT_0016 + (-2)*IT_0176;
    const ccomplex_t IT_0228 = s_14*s_34*IT_0077;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = s_13 + IT_0229;
    const ccomplex_t IT_0231 = s_13*IT_0230;
    const ccomplex_t IT_0232 = IT_0112 + 2*IT_0205;
    const ccomplex_t IT_0233 = IT_0156*IT_0230;
    const ccomplex_t IT_0234 = -IT_0176;
    const ccomplex_t IT_0235 = m_N_3*IT_0230;
    const ccomplex_t IT_0236 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0237 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0238 = pow(m_Z, 2);
    const ccomplex_t IT_0239 = s_14*IT_0238;
    const ccomplex_t IT_0240 = IT_0077*IT_0239;
    const ccomplex_t IT_0241 = (-0.25)*IT_0240;
    const ccomplex_t IT_0242 = s_14 + IT_0241;
    const ccomplex_t IT_0243 = m_c*IT_0242;
    const ccomplex_t IT_0244 = (-4)*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0240;
    const ccomplex_t IT_0246 = s_14 + IT_0245;
    const ccomplex_t IT_0247 = m_c*IT_0246;
    const ccomplex_t IT_0248 = 2*IT_0247;
    const ccomplex_t IT_0249 = IT_0101 + (-0.5)*IT_0115 + 0.5*IT_0116 + 
      -IT_0135 + -IT_0154;
    const ccomplex_t IT_0250 = s_34*IT_0238;
    const ccomplex_t IT_0251 = IT_0077*IT_0250;
    const ccomplex_t IT_0252 = (-0.25)*IT_0251;
    const ccomplex_t IT_0253 = s_34 + IT_0252;
    const ccomplex_t IT_0254 = m_N_3*IT_0253;
    const ccomplex_t IT_0255 = (-4)*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0251;
    const ccomplex_t IT_0257 = s_34 + IT_0256;
    const ccomplex_t IT_0258 = m_N_3*IT_0257;
    const ccomplex_t IT_0259 = 2*IT_0258;
    const ccomplex_t IT_0260 = s_13*IT_0238;
    const ccomplex_t IT_0261 = s_14*s_34;
    const ccomplex_t IT_0262 = (-2)*IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = IT_0077*IT_0263;
    const ccomplex_t IT_0265 = (-6)*IT_0264;
    const ccomplex_t IT_0266 = 12*s_13;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = m_c*IT_0230;
    const ccomplex_t IT_0269 = (-6)*IT_0227;
    const ccomplex_t IT_0270 = pow(s_34, 2);
    const ccomplex_t IT_0271 = IT_0077*IT_0270;
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = IT_0011 + IT_0272;
    const ccomplex_t IT_0274 = IT_0156*IT_0273;
    const ccomplex_t IT_0275 = s_14*IT_0011;
    const ccomplex_t IT_0276 = s_13*s_34;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = IT_0275 + IT_0277;
    const ccomplex_t IT_0279 = m_N_3*IT_0273;
    const ccomplex_t IT_0280 = (-6)*conj(IT_0227);
    const ccomplex_t IT_0281 = (-6)*conj(IT_0232);
    const ccomplex_t IT_0282 = (-6)*conj(IT_0193);
    const ccomplex_t IT_0283 = (-6)*IT_0076;
    const ccomplex_t IT_0284 = (-6)*IT_0232;
    const ccomplex_t IT_0285 = 6*IT_0158;
    const ccomplex_t IT_0286 = 6*IT_0249;
    const ccomplex_t IT_0287 = 6*IT_0205;
    const ccomplex_t IT_0288 = 6*IT_0155*IT_0163 + IT_0156*((-18)*IT_0193 + (
      -18)*IT_0199) + 6*IT_0234*IT_0259 + IT_0225*IT_0267 + IT_0269*IT_0279 +
       IT_0235*IT_0283 + IT_0268*IT_0284 + IT_0255*IT_0285 + IT_0248*IT_0286 +
       IT_0244*IT_0287;
    const ccomplex_t IT_0289 = (-6)*IT_0155;
    const ccomplex_t IT_0290 = 6*IT_0076*IT_0163 + IT_0156*((-18)*IT_0218 + (
      -18)*IT_0224) + 6*IT_0234*IT_0244 + IT_0200*IT_0267 + IT_0268*IT_0269 +
       IT_0279*IT_0284 + IT_0248*IT_0285 + IT_0255*IT_0286 + IT_0259*IT_0287 +
       IT_0235*IT_0289;
    const ccomplex_t IT_0291 = IT_0156*IT_0238;
    const ccomplex_t IT_0292 = 3*IT_0291;
    const ccomplex_t IT_0293 = (-3)*IT_0291;
    const ccomplex_t IT_0294 = pow(m_Z, 4);
    const ccomplex_t IT_0295 = s_13*IT_0294;
    const ccomplex_t IT_0296 = IT_0077*IT_0295;
    const ccomplex_t IT_0297 = (-24)*IT_0261;
    const ccomplex_t IT_0298 = 6*IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 24*IT_0261;
    const ccomplex_t IT_0300 = (-12)*IT_0260;
    const ccomplex_t IT_0301 = 6*IT_0296 + IT_0299 + IT_0300;
    const ccomplex_t IT_0302 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0303 = (-6)*IT_0193;
    const ccomplex_t IT_0304 = 0.166666666666667*IT_0233;
    const ccomplex_t IT_0305 = 0.166666666666667*IT_0249;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0205;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0231;
    const ccomplex_t IT_0308 = s_13*IT_0273;
    const ccomplex_t IT_0309 = 0.166666666666667*IT_0308;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0311 = 0.166666666666667*IT_0278;
    const ccomplex_t IT_0312 = 6*conj(IT_0076)*(IT_0076*IT_0082 + IT_0155
      *IT_0157 + IT_0158*IT_0162 + IT_0163*IT_0200) + IT_0156*IT_0225*IT_0226 + 
      (IT_0227*IT_0231 + IT_0232*IT_0233 + IT_0162*IT_0234 + IT_0225*IT_0235)
      *IT_0236 + (IT_0162*IT_0205 + IT_0231*IT_0232 + IT_0199*IT_0235)*IT_0237 +
       6*conj(IT_0193)*(IT_0076*IT_0163 + IT_0234*IT_0244 + IT_0158*IT_0248 +
       IT_0249*IT_0255 + IT_0205*IT_0259 + 0.166666666666667*IT_0200*IT_0267 +
       0.166666666666667*IT_0268*IT_0269) + (IT_0155*IT_0233 + IT_0200*IT_0268 +
       IT_0232*IT_0274 + IT_0234*IT_0278 + IT_0225*IT_0279)*IT_0280 + (IT_0225
      *IT_0268 + IT_0205*IT_0278 + IT_0199*IT_0279)*IT_0281 + (IT_0155*IT_0235 +
       IT_0232*IT_0279)*IT_0282 + (conj(IT_0218) + conj(IT_0224))*IT_0288 + conj
      (IT_0199)*IT_0290 + 6*conj(IT_0158)*(IT_0076*IT_0162 + IT_0200*IT_0248 +
       IT_0225*IT_0255 + IT_0227*IT_0278 + IT_0205*IT_0292 + IT_0249*IT_0293 +
       0.166666666666667*IT_0234*IT_0298 + IT_0301*IT_0302) + 6*conj(IT_0155)*
      (IT_0082*IT_0155 + IT_0076*IT_0157 + IT_0163*IT_0225 + IT_0162*IT_0249 +
       0.166666666666667*IT_0235*IT_0303 + IT_0269*IT_0304) + 6*conj(IT_0249)*
      (IT_0155*IT_0162 + IT_0225*IT_0248 + IT_0200*IT_0255 + IT_0232*IT_0278 +
       IT_0234*IT_0292 + IT_0158*IT_0293 + IT_0301*IT_0305 + IT_0298*IT_0306) +
       6*conj(IT_0227)*(IT_0158*IT_0278 + IT_0283*IT_0307 + IT_0269*IT_0309) + 6
      *conj(IT_0232)*(0.166666666666667*IT_0269*IT_0274 + IT_0249*IT_0278 +
       0.166666666666667*IT_0279*IT_0303 + IT_0283*IT_0304 + IT_0289*IT_0307 +
       IT_0284*IT_0309) + 6*conj(IT_0234)*(IT_0200*IT_0244 + IT_0225*IT_0259 +
       IT_0249*IT_0292 + IT_0205*IT_0293 + 0.166666666666667*IT_0234*IT_0301 +
       IT_0298*IT_0302 + IT_0283*IT_0310 + IT_0269*IT_0311) + 6*conj(IT_0205)*
      (IT_0225*IT_0244 + IT_0200*IT_0259 + IT_0158*IT_0292 + IT_0234*IT_0293 +
       IT_0298*IT_0305 + IT_0301*IT_0306 + IT_0289*IT_0310 + IT_0284*IT_0311);
    return create_ccomplex_return(IT_0312);
}

