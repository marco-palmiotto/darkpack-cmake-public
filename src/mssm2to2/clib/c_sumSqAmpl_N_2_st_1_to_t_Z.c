#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_1_to_t_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_1_to_t_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_t = param->m_t;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0000, -1);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0002*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_t*N_u2*e_em*IT_0020*U_st_00;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + (-0.75)*IT_0023);
    const ccomplex_t IT_0025 = 0.666666666666667*IT_0024;
    const ccomplex_t IT_0026 = pow(m_t, 2);
    const ccomplex_t IT_0027 = pow(m_st_1, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0013*IT_0030;
    const ccomplex_t IT_0032 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = IT_0003*IT_0032;
    const ccomplex_t IT_0035 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0036 = IT_0007*IT_0035;
    const ccomplex_t IT_0037 = IT_0003*IT_0035;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0033 + IT_0034 + 
      -IT_0036 + -IT_0037);
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0042 = IT_0006*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = m_t*N_u3*e_em*IT_0020*U_st_00;
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + (-0.75)*IT_0046);
    const ccomplex_t IT_0048 = 0.666666666666667*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0054 = IT_0007*IT_0053;
    const ccomplex_t IT_0055 = IT_0003*IT_0053;
    const ccomplex_t IT_0056 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = IT_0007*IT_0056;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0054 + IT_0055 + 
      -IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0063 = IT_0006*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = m_t*N_u4*e_em*IT_0020*U_st_00;
    const ccomplex_t IT_0066 = IT_0019*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + (-0.75)*IT_0067);
    const ccomplex_t IT_0069 = 0.666666666666667*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0061*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0031 + -IT_0052 + -IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0012;
    const ccomplex_t IT_0076 = IT_0030*IT_0075;
    const ccomplex_t IT_0077 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0078 = IT_0007*IT_0077;
    const ccomplex_t IT_0079 = IT_0003*IT_0077;
    const ccomplex_t IT_0080 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0081 = IT_0007*IT_0080;
    const ccomplex_t IT_0082 = IT_0003*IT_0080;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0078 + IT_0079 + 
      -IT_0081 + -IT_0082);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0086 = IT_0006*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = m_t*N_u1*e_em*IT_0020*U_st_00;
    const ccomplex_t IT_0089 = IT_0019*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.75)*IT_0090);
    const ccomplex_t IT_0092 = 0.666666666666667*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0084*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0076 + -IT_0096;
    const ccomplex_t IT_0098 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0099 = IT_0006*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0102 = IT_0002*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_t*conj(N_u2)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0105 = IT_0019*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103 + 3
      *IT_0106);
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = IT_0028*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0075*IT_0110;
    const ccomplex_t IT_0112 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0113 = IT_0003*IT_0112;
    const ccomplex_t IT_0114 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0115 = IT_0007*IT_0114;
    const ccomplex_t IT_0116 = IT_0003*IT_0114;
    const ccomplex_t IT_0117 = IT_0007*IT_0112;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0113 + -IT_0115 + 
      -IT_0116 + IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0121 = IT_0006*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0124 = IT_0002*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = m_t*conj(N_u3)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0127 = IT_0019*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + 3*IT_0125 + 3
      *IT_0128);
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0129;
    const ccomplex_t IT_0131 = IT_0049*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0119*IT_0132;
    const ccomplex_t IT_0134 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0135 = IT_0007*IT_0134;
    const ccomplex_t IT_0136 = IT_0003*IT_0134;
    const ccomplex_t IT_0137 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0138 = IT_0003*IT_0137;
    const ccomplex_t IT_0139 = IT_0007*IT_0137;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0135 + IT_0136 + 
      -IT_0138 + -IT_0139);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0146 = IT_0006*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0149 = IT_0002*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = m_t*conj(N_u4)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0152 = IT_0019*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + 3*IT_0150 + 3
      *IT_0153);
    const ccomplex_t IT_0155 = 0.166666666666667*IT_0154;
    const ccomplex_t IT_0156 = IT_0070*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = IT_0144*IT_0157;
    const ccomplex_t IT_0159 = 0.5*IT_0111 + IT_0133 + IT_0158;
    const ccomplex_t IT_0160 = IT_0013*IT_0110;
    const ccomplex_t IT_0161 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0162 = IT_0003*IT_0161;
    const ccomplex_t IT_0163 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0164 = IT_0003*IT_0163;
    const ccomplex_t IT_0165 = IT_0007*IT_0161;
    const ccomplex_t IT_0166 = IT_0007*IT_0163;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0162 + -IT_0164 +
       IT_0165 + -IT_0166);
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0171 = IT_0006*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0174 = IT_0002*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = m_t*conj(N_u1)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0177 = IT_0019*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0172 + 3*IT_0175 + 3
      *IT_0178);
    const ccomplex_t IT_0180 = 0.166666666666667*IT_0179;
    const ccomplex_t IT_0181 = IT_0093*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = IT_0169*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0160 + IT_0183;
    const ccomplex_t IT_0185 = pow(m_Z, -2);
    const ccomplex_t IT_0186 = pow(m_Z, 4);
    const ccomplex_t IT_0187 = s_13*IT_0186;
    const ccomplex_t IT_0188 = IT_0185*IT_0187;
    const ccomplex_t IT_0189 = (-6)*IT_0188;
    const ccomplex_t IT_0190 = s_14*s_34;
    const ccomplex_t IT_0191 = (-24)*IT_0190;
    const ccomplex_t IT_0192 = pow(m_Z, 2);
    const ccomplex_t IT_0193 = s_13*IT_0192;
    const ccomplex_t IT_0194 = 12*IT_0193;
    const ccomplex_t IT_0195 = IT_0189 + IT_0191 + IT_0194;
    const ccomplex_t IT_0196 = 6*IT_0188;
    const ccomplex_t IT_0197 = 24*IT_0190;
    const ccomplex_t IT_0198 = (-12)*IT_0193;
    const ccomplex_t IT_0199 = IT_0196 + IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = (-2)*IT_0190;
    const ccomplex_t IT_0201 = IT_0193 + IT_0200;
    const ccomplex_t IT_0202 = IT_0185*IT_0201;
    const ccomplex_t IT_0203 = (-6)*IT_0202;
    const ccomplex_t IT_0204 = 12*s_13;
    const ccomplex_t IT_0205 = IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = e_em*IT_0007;
    const ccomplex_t IT_0207 = e_em*IT_0003;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0206 + (-3)*IT_0207);
    const ccomplex_t IT_0209 = (-0.166666666666667)*IT_0208;
    const ccomplex_t IT_0210 = IT_0025*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0212 = pow(m_N_2, 2);
    const ccomplex_t IT_0213 = cpow((-2)*s_12 + -IT_0027 + -IT_0212 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0214 = IT_0210*IT_0211*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0216 = (-0.666666666666667)*IT_0215;
    const ccomplex_t IT_0217 = IT_0025*IT_0216;
    const ccomplex_t IT_0218 = IT_0213*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_t*IT_0219;
    const ccomplex_t IT_0221 = m_N_3*IT_0040;
    const ccomplex_t IT_0222 = m_N_2*IT_0119;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = IT_0221 + IT_0223;
    const ccomplex_t IT_0225 = IT_0132*IT_0224;
    const ccomplex_t IT_0226 = m_N_2*IT_0013;
    const ccomplex_t IT_0227 = m_N_2*IT_0075;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = IT_0110*IT_0229;
    const ccomplex_t IT_0231 = m_N_4*IT_0061;
    const ccomplex_t IT_0232 = m_N_2*IT_0144;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0231 + IT_0233;
    const ccomplex_t IT_0235 = IT_0157*IT_0234;
    const ccomplex_t IT_0236 = -IT_0214 + IT_0220 + -IT_0225 + (-0.5)*IT_0230 
      + -IT_0235;
    const ccomplex_t IT_0237 = m_N_1*IT_0084;
    const ccomplex_t IT_0238 = m_N_2*IT_0169;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = IT_0182*IT_0240;
    const ccomplex_t IT_0242 = (-0.5)*IT_0230 + IT_0241;
    const ccomplex_t IT_0243 = IT_0236 + IT_0242;
    const ccomplex_t IT_0244 = conj(IT_0236) + conj(IT_0242);
    const ccomplex_t IT_0245 = IT_0108*IT_0216;
    const ccomplex_t IT_0246 = IT_0211*IT_0213*IT_0245;
    const ccomplex_t IT_0247 = IT_0108*IT_0209;
    const ccomplex_t IT_0248 = IT_0213*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = m_t*IT_0249;
    const ccomplex_t IT_0251 = IT_0030*IT_0229;
    const ccomplex_t IT_0252 = m_N_3*IT_0119;
    const ccomplex_t IT_0253 = m_N_2*IT_0040;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = IT_0252 + IT_0254;
    const ccomplex_t IT_0256 = IT_0051*IT_0255;
    const ccomplex_t IT_0257 = m_N_4*IT_0144;
    const ccomplex_t IT_0258 = m_N_2*IT_0061;
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = IT_0257 + IT_0259;
    const ccomplex_t IT_0261 = IT_0072*IT_0260;
    const ccomplex_t IT_0262 = IT_0246 + -IT_0250 + (-0.5)*IT_0251 + IT_0256 +
       IT_0261;
    const ccomplex_t IT_0263 = m_N_1*IT_0169;
    const ccomplex_t IT_0264 = m_N_2*IT_0084;
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = IT_0095*IT_0266;
    const ccomplex_t IT_0268 = (-0.5)*IT_0251 + -IT_0267;
    const ccomplex_t IT_0269 = IT_0262 + IT_0268;
    const ccomplex_t IT_0270 = conj(IT_0262) + conj(IT_0268);
    const ccomplex_t IT_0271 = m_t*m_N_2;
    const ccomplex_t IT_0272 = IT_0192*IT_0271;
    const ccomplex_t IT_0273 = (-3)*IT_0272;
    const ccomplex_t IT_0274 = 6*IT_0074;
    const ccomplex_t IT_0275 = (-0.166666666666667)*IT_0274;
    const ccomplex_t IT_0276 = IT_0097 + IT_0275;
    const ccomplex_t IT_0277 = 6*conj(IT_0074);
    const ccomplex_t IT_0278 = (-0.166666666666667)*IT_0277;
    const ccomplex_t IT_0279 = conj(IT_0097) + IT_0278;
    const ccomplex_t IT_0280 = pow(s_14, 2);
    const ccomplex_t IT_0281 = IT_0185*IT_0280;
    const ccomplex_t IT_0282 = -IT_0212;
    const ccomplex_t IT_0283 = IT_0281 + IT_0282;
    const ccomplex_t IT_0284 = IT_0271*IT_0283;
    const ccomplex_t IT_0285 = -IT_0160 + (-2)*IT_0183;
    const ccomplex_t IT_0286 = IT_0076 + 2*IT_0096;
    const ccomplex_t IT_0287 = 6*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0000
      *IT_0002 + -1./3*IT_0001*IT_0006)*U_st_00*conj(U_st_01) + IT_0001*IT_0006
      *U_st_10*conj(U_st_11));
    const ccomplex_t IT_0289 = 1.33333333333333*IT_0288;
    const ccomplex_t IT_0290 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0291 = IT_0006*IT_0290;
    const ccomplex_t IT_0292 = 1.4142135623731*IT_0291;
    const ccomplex_t IT_0293 = m_t*N_u2*e_em*IT_0020*U_st_01;
    const ccomplex_t IT_0294 = IT_0019*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*(IT_0292 + (-0.75)*IT_0295);
    const ccomplex_t IT_0297 = 0.666666666666667*IT_0296;
    const ccomplex_t IT_0298 = cpow((-2)*s_13 + IT_0026 + IT_0212 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0299 = IT_0297*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0289*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0000
      *IT_0002 + -1./3*IT_0001*IT_0006)*U_st_00*conj(U_st_00) + IT_0001*IT_0006
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0303 = 1.33333333333333*IT_0302;
    const ccomplex_t IT_0304 = cpow((-2)*s_13 + IT_0026 + IT_0212 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0305 = IT_0025*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = IT_0303*IT_0306;
    const ccomplex_t IT_0308 = IT_0031 + 2*IT_0052 + 2*IT_0073 + IT_0301 +
       IT_0307;
    const ccomplex_t IT_0309 = (-6)*IT_0308;
    const ccomplex_t IT_0310 = IT_0287 + IT_0309;
    const ccomplex_t IT_0311 = 6*conj(IT_0286);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0308);
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = IT_0108*IT_0304;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0303*IT_0315;
    const ccomplex_t IT_0317 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0318 = IT_0006*IT_0317;
    const ccomplex_t IT_0319 = 1.4142135623731*IT_0318;
    const ccomplex_t IT_0320 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0321 = IT_0002*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = m_t*conj(N_u2)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0324 = IT_0019*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*(IT_0319 + 3*IT_0322 + 3
      *IT_0325);
    const ccomplex_t IT_0327 = 0.166666666666667*IT_0326;
    const ccomplex_t IT_0328 = IT_0298*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = IT_0289*IT_0329;
    const ccomplex_t IT_0331 = -IT_0111 + (-2)*IT_0133 + (-2)*IT_0158 + 
      -IT_0316 + -IT_0330;
    const ccomplex_t IT_0332 = 6*IT_0308;
    const ccomplex_t IT_0333 = (-6)*IT_0286;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = 6*conj(IT_0308);
    const ccomplex_t IT_0336 = (-6)*conj(IT_0286);
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = m_t*IT_0283;
    const ccomplex_t IT_0339 = 6*IT_0331;
    const ccomplex_t IT_0340 = (-6)*IT_0285;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 6*conj(IT_0331);
    const ccomplex_t IT_0343 = (-6)*conj(IT_0285);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = s_14*s_34*IT_0185;
    const ccomplex_t IT_0346 = -IT_0345;
    const ccomplex_t IT_0347 = s_13 + IT_0346;
    const ccomplex_t IT_0348 = s_13*IT_0347;
    const ccomplex_t IT_0349 = (-2)*IT_0219 + -IT_0301 + -IT_0307;
    const ccomplex_t IT_0350 = 2*IT_0249 + IT_0316 + IT_0330;
    const ccomplex_t IT_0351 = 6*IT_0285;
    const ccomplex_t IT_0352 = (-6)*IT_0331;
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = 6*conj(IT_0285);
    const ccomplex_t IT_0355 = (-6)*conj(IT_0331);
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = s_13*s_14;
    const ccomplex_t IT_0358 = s_34*IT_0212;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = -IT_0219;
    const ccomplex_t IT_0362 = conj(IT_0097) + conj(IT_0361);
    const ccomplex_t IT_0363 = IT_0097 + IT_0361;
    const ccomplex_t IT_0364 = conj(IT_0184) + conj(IT_0249);
    const ccomplex_t IT_0365 = IT_0184 + IT_0249;
    const ccomplex_t IT_0366 = s_34*IT_0192;
    const ccomplex_t IT_0367 = IT_0185*IT_0366;
    const ccomplex_t IT_0368 = (-0.25)*IT_0367;
    const ccomplex_t IT_0369 = s_34 + IT_0368;
    const ccomplex_t IT_0370 = m_N_2*IT_0369;
    const ccomplex_t IT_0371 = (-4)*IT_0370;
    const ccomplex_t IT_0372 = 6*IT_0236;
    const ccomplex_t IT_0373 = 6*IT_0242;
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = 6*conj(IT_0236);
    const ccomplex_t IT_0376 = 6*conj(IT_0242);
    const ccomplex_t IT_0377 = IT_0375 + IT_0376;
    const ccomplex_t IT_0378 = 6*IT_0262;
    const ccomplex_t IT_0379 = 6*IT_0268;
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = 6*conj(IT_0262);
    const ccomplex_t IT_0382 = 6*conj(IT_0268);
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = (-6)*IT_0236;
    const ccomplex_t IT_0385 = (-6)*IT_0242;
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = (-6)*conj(IT_0236);
    const ccomplex_t IT_0388 = (-6)*conj(IT_0242);
    const ccomplex_t IT_0389 = IT_0387 + IT_0388;
    const ccomplex_t IT_0390 = (-6)*IT_0262;
    const ccomplex_t IT_0391 = (-6)*IT_0268;
    const ccomplex_t IT_0392 = IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = (-6)*conj(IT_0262);
    const ccomplex_t IT_0394 = (-6)*conj(IT_0268);
    const ccomplex_t IT_0395 = IT_0393 + IT_0394;
    const ccomplex_t IT_0396 = s_14*IT_0192;
    const ccomplex_t IT_0397 = IT_0185*IT_0396;
    const ccomplex_t IT_0398 = 0.5*IT_0397;
    const ccomplex_t IT_0399 = s_14 + IT_0398;
    const ccomplex_t IT_0400 = m_t*IT_0399;
    const ccomplex_t IT_0401 = 2*IT_0400;
    const ccomplex_t IT_0402 = 3*IT_0272;
    const ccomplex_t IT_0403 = 6*IT_0361;
    const ccomplex_t IT_0404 = 6*conj(IT_0361);
    const ccomplex_t IT_0405 = (-18)*conj(IT_0236);
    const ccomplex_t IT_0406 = s_13*IT_0283;
    const ccomplex_t IT_0407 = IT_0271*IT_0347;
    const ccomplex_t IT_0408 = 6*IT_0349;
    const ccomplex_t IT_0409 = 6*conj(IT_0349);
    const ccomplex_t IT_0410 = (-6)*IT_0349;
    const ccomplex_t IT_0411 = (-6)*conj(IT_0349);
    const ccomplex_t IT_0412 = pow(s_34, 2);
    const ccomplex_t IT_0413 = IT_0185*IT_0412;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0026 + IT_0414;
    const ccomplex_t IT_0416 = IT_0271*IT_0415;
    const ccomplex_t IT_0417 = s_13*IT_0415;
    const ccomplex_t IT_0418 = (-6)*IT_0350;
    const ccomplex_t IT_0419 = m_t*IT_0347;
    const ccomplex_t IT_0420 = (-6)*conj(IT_0350);
    const ccomplex_t IT_0421 = s_14*IT_0026;
    const ccomplex_t IT_0422 = s_13*s_34;
    const ccomplex_t IT_0423 = -IT_0422;
    const ccomplex_t IT_0424 = IT_0421 + IT_0423;
    const ccomplex_t IT_0425 = conj(IT_0159)*IT_0350 + IT_0159*conj(IT_0350) +
       0.166666666666667*conj(IT_0074)*IT_0408 + 0.166666666666667*IT_0074
      *IT_0409 + 0.166666666666667*IT_0362*IT_0410 + 0.166666666666667*IT_0363
      *IT_0411 + 0.166666666666667*IT_0364*IT_0418 + 0.166666666666667*IT_0365
      *IT_0420;
    const ccomplex_t IT_0426 = m_N_2*IT_0347;
    const ccomplex_t IT_0427 = m_N_2*IT_0415;
    const ccomplex_t IT_0428 = IT_0189 + IT_0197;
    const ccomplex_t IT_0429 = IT_0191 + IT_0196;
    const ccomplex_t IT_0430 = (-0.25)*IT_0397;
    const ccomplex_t IT_0431 = s_14 + IT_0430;
    const ccomplex_t IT_0432 = m_t*IT_0431;
    const ccomplex_t IT_0433 = (-4)*IT_0432;
    const ccomplex_t IT_0434 = (-0.166666666666667)*IT_0433;
    const ccomplex_t IT_0435 = 0.5*IT_0367;
    const ccomplex_t IT_0436 = s_34 + IT_0435;
    const ccomplex_t IT_0437 = m_N_2*IT_0436;
    const ccomplex_t IT_0438 = 2*IT_0437;
    const ccomplex_t IT_0439 = (-0.166666666666667)*IT_0438;
    const ccomplex_t IT_0440 = (-0.166666666666667)*IT_0273;
    const ccomplex_t IT_0441 = (conj(IT_0074)*IT_0097 + IT_0074*conj(IT_0097) 
      + conj(IT_0159)*IT_0184 + IT_0159*conj(IT_0184))*IT_0195 + (IT_0074*conj
      (IT_0074) + IT_0097*conj(IT_0097) + IT_0159*conj(IT_0159) + IT_0184*conj
      (IT_0184))*IT_0199 + IT_0205*(IT_0243*IT_0244 + IT_0269*IT_0270) + (-6)
      *IT_0273*((conj(IT_0074) + -conj(IT_0097))*IT_0184 + (IT_0074 + -IT_0097)
      *conj(IT_0184) + conj(IT_0159)*IT_0276 + IT_0159*IT_0279) + IT_0284*(conj
      (IT_0285)*IT_0310 + IT_0285*IT_0313 + conj(IT_0331)*IT_0334 + IT_0331
      *IT_0337) + IT_0338*(IT_0244*IT_0334 + IT_0243*IT_0337 + IT_0270*IT_0341 +
       IT_0269*IT_0344) + IT_0348*(IT_0313*IT_0349 + IT_0310*conj(IT_0349) +
       conj(IT_0350)*IT_0353 + IT_0350*IT_0356) + IT_0360*(conj(IT_0074)*IT_0334
       + IT_0074*IT_0337 + conj(IT_0159)*IT_0341 + IT_0159*IT_0344 + IT_0310
      *IT_0362 + IT_0313*IT_0363 + IT_0353*IT_0364 + IT_0356*IT_0365) + IT_0371*
      (conj(IT_0159)*IT_0374 + IT_0159*IT_0377 + conj(IT_0074)*IT_0380 + IT_0074
      *IT_0383 + conj(IT_0184)*IT_0386 + IT_0184*IT_0389 + conj(IT_0097)*IT_0392
       + IT_0097*IT_0395) + (conj(IT_0074)*IT_0374 + IT_0074*IT_0377 + conj
      (IT_0159)*IT_0380 + IT_0159*IT_0383 + conj(IT_0097)*IT_0386 + IT_0097
      *IT_0389 + conj(IT_0184)*IT_0392 + IT_0184*IT_0395)*IT_0401 + IT_0402*
      (conj(IT_0159)*IT_0403 + IT_0159*IT_0404) + IT_0271*(IT_0243*((-18)*conj
      (IT_0262) + (-18)*conj(IT_0268)) + IT_0269*((-18)*conj(IT_0242) + IT_0405)
      ) + (conj(IT_0286)*IT_0310 + IT_0286*IT_0312 + conj(IT_0308)*IT_0332 +
       conj(IT_0331)*IT_0339 + conj(IT_0285)*IT_0353 + IT_0285*IT_0355)*IT_0406 
      + IT_0407*(IT_0313*IT_0350 + IT_0310*conj(IT_0350) + conj(IT_0285)*IT_0408
       + IT_0285*IT_0409 + conj(IT_0331)*IT_0410 + IT_0331*IT_0411) + (conj
      (IT_0350)*IT_0410 + IT_0350*IT_0411)*IT_0416 + IT_0417*(conj(IT_0349)
      *IT_0410 + conj(IT_0350)*IT_0418) + IT_0419*(IT_0244*IT_0410 + IT_0243
      *IT_0411 + IT_0270*IT_0418 + IT_0269*IT_0420) + 6*IT_0424*IT_0425 + 
      (IT_0270*IT_0310 + IT_0269*IT_0313 + conj(IT_0242)*IT_0353 + IT_0242
      *IT_0356 + conj(IT_0285)*IT_0372 + IT_0285*IT_0375 + conj(IT_0331)*IT_0384
       + IT_0331*IT_0387)*IT_0426 + (conj(IT_0350)*IT_0384 + IT_0350*IT_0387 +
       IT_0270*IT_0410 + IT_0269*IT_0411 + conj(IT_0242)*IT_0418 + IT_0242
      *IT_0420)*IT_0427 + (-6)*conj(IT_0361)*((-0.166666666666667)*IT_0199
      *IT_0361 + IT_0184*IT_0402 + (-0.166666666666667)*IT_0097*IT_0428 + (
      -0.166666666666667)*IT_0074*IT_0429 + IT_0374*IT_0434 + IT_0380*IT_0439) +
       (-6)*IT_0361*(conj(IT_0184)*IT_0402 + (-0.166666666666667)*conj(IT_0097)
      *IT_0428 + (-0.166666666666667)*conj(IT_0074)*IT_0429 + IT_0377*IT_0434 +
       IT_0383*IT_0439) + (-6)*conj(IT_0249)*(IT_0276*IT_0402 + -1./6*IT_0184
      *IT_0428 + -1./6*IT_0159*IT_0429 + IT_0380*IT_0434 + IT_0374*IT_0439 +
       IT_0403*IT_0440) + IT_0249*(IT_0199*conj(IT_0249) + (-6)*IT_0279*IT_0402 
      + conj(IT_0184)*IT_0428 + conj(IT_0159)*IT_0429 + (-6)*IT_0383*IT_0434 + (
      -6)*IT_0377*IT_0439 + (-6)*IT_0404*IT_0440);
    return create_ccomplex_return(IT_0441);
}

