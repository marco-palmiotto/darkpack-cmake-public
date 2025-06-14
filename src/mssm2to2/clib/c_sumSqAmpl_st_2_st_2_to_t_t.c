#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_st_2_to_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_st_2_to_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_st_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = sin(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = m_t*conj(N_u2)*e_em*IT_0025*U_st_11;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0016 + 3*IT_0021 + 3
      *IT_0028);
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0032 = IT_0013*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = m_t*N_u2*e_em*IT_0025*U_st_01;
    const ccomplex_t IT_0035 = IT_0024*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + (-0.75)*IT_0036);
    const ccomplex_t IT_0038 = 0.666666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0030*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0044 = IT_0013*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0047 = IT_0018*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = m_t*conj(N_u1)*e_em*IT_0025*U_st_11;
    const ccomplex_t IT_0050 = IT_0024*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0045 + 3*IT_0048 + 3
      *IT_0051);
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0055 = IT_0013*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = m_t*N_u1*e_em*IT_0025*U_st_01;
    const ccomplex_t IT_0058 = IT_0024*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + (-0.75)*IT_0059);
    const ccomplex_t IT_0061 = 0.666666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0067 = IT_0013*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0070 = IT_0018*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = m_t*conj(N_u3)*e_em*IT_0025*U_st_11;
    const ccomplex_t IT_0073 = IT_0024*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + 3*IT_0071 + 3
      *IT_0074);
    const ccomplex_t IT_0076 = 0.166666666666667*IT_0075;
    const ccomplex_t IT_0077 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0078 = IT_0013*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = m_t*N_u3*e_em*IT_0025*U_st_01;
    const ccomplex_t IT_0081 = IT_0024*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + (-0.75)*IT_0082);
    const ccomplex_t IT_0084 = 0.666666666666667*IT_0083;
    const ccomplex_t IT_0085 = IT_0076*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0090 = IT_0013*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0093 = IT_0018*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = m_t*conj(N_u4)*e_em*IT_0025*U_st_11;
    const ccomplex_t IT_0096 = IT_0024*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0091 + 3*IT_0094 + 3
      *IT_0097);
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0101 = IT_0013*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = m_t*N_u4*e_em*IT_0025*U_st_01;
    const ccomplex_t IT_0104 = IT_0024*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0102 + (-0.75)*IT_0105);
    const ccomplex_t IT_0107 = 0.666666666666667*IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 0.166666666666667*IT_0008 + 0.5*IT_0011 + 
      -IT_0042 + -IT_0065 + -IT_0088 + -IT_0111;
    const ccomplex_t IT_0113 = s_34*IT_0002;
    const ccomplex_t IT_0114 = (-18)*IT_0113;
    const ccomplex_t IT_0115 = s_23*s_24;
    const ccomplex_t IT_0116 = 36*IT_0115;
    const ccomplex_t IT_0117 = IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = (-6)*IT_0113;
    const ccomplex_t IT_0119 = 12*IT_0115;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0062*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0125 = IT_0039*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0085*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0131 = IT_0108*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0008 + (-0.166666666666667)*IT_0011 +
       IT_0123 + IT_0126 + IT_0129 + IT_0132;
    const ccomplex_t IT_0134 = 2*IT_0112;
    const ccomplex_t IT_0135 = IT_0000*IT_0002;
    const ccomplex_t IT_0136 = cpow(IT_0005, 2);
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*m_sG*IT_0136;
    const ccomplex_t IT_0138 = IT_0003*IT_0137;
    const ccomplex_t IT_0139 = cpow(IT_0084, 2);
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*m_N_3*IT_0139;
    const ccomplex_t IT_0141 = IT_0127*IT_0140;
    const ccomplex_t IT_0142 = m_t*IT_0008;
    const ccomplex_t IT_0143 = m_t*IT_0011;
    const ccomplex_t IT_0144 = m_t*IT_0123;
    const ccomplex_t IT_0145 = m_t*IT_0126;
    const ccomplex_t IT_0146 = m_t*IT_0129;
    const ccomplex_t IT_0147 = m_t*IT_0132;
    const ccomplex_t IT_0148 = 0.5*IT_0138 + IT_0141 + (-0.5)*IT_0142 +
       0.166666666666667*IT_0143 + -IT_0144 + -IT_0145 + -IT_0146 + -IT_0147;
    const ccomplex_t IT_0149 = cpow(IT_0038, 2);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*m_N_2*IT_0149;
    const ccomplex_t IT_0151 = IT_0040*IT_0150;
    const ccomplex_t IT_0152 = m_t*IT_0065;
    const ccomplex_t IT_0153 = m_t*IT_0042;
    const ccomplex_t IT_0154 = m_t*IT_0088;
    const ccomplex_t IT_0155 = m_t*IT_0111;
    const ccomplex_t IT_0156 = (-0.166666666666667)*IT_0138 +
       0.166666666666667*IT_0142 + (-0.5)*IT_0143 + IT_0151 + -IT_0152 + 
      -IT_0153 + -IT_0154 + -IT_0155;
    const ccomplex_t IT_0157 = IT_0009*IT_0137;
    const ccomplex_t IT_0158 = cpow(IT_0061, 2);
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*m_N_1*IT_0158;
    const ccomplex_t IT_0160 = IT_0121*IT_0159;
    const ccomplex_t IT_0161 = IT_0124*IT_0150;
    const ccomplex_t IT_0162 = cpow(IT_0107, 2);
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_4*IT_0162;
    const ccomplex_t IT_0164 = IT_0130*IT_0163;
    const ccomplex_t IT_0165 = 0.166666666666667*IT_0157 + -IT_0160 + -IT_0161
       + -IT_0164;
    const ccomplex_t IT_0166 = IT_0086*IT_0140;
    const ccomplex_t IT_0167 = IT_0063*IT_0159;
    const ccomplex_t IT_0168 = IT_0109*IT_0163;
    const ccomplex_t IT_0169 = (-0.5)*IT_0157 + -IT_0166 + -IT_0167 + -IT_0168;
    const ccomplex_t IT_0170 = cpow(IT_0053, 2);
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_1*IT_0170;
    const ccomplex_t IT_0172 = IT_0121*IT_0171;
    const ccomplex_t IT_0173 = cpow(IT_0030, 2);
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*m_N_2*IT_0173;
    const ccomplex_t IT_0175 = IT_0124*IT_0174;
    const ccomplex_t IT_0176 = (-0.5)*IT_0142 + 0.166666666666667*IT_0143 + 
      -IT_0144 + -IT_0145 + -IT_0146 + -IT_0147 + IT_0172 + IT_0175;
    const ccomplex_t IT_0177 = cpow(IT_0076, 2);
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*m_N_3*IT_0177;
    const ccomplex_t IT_0179 = IT_0086*IT_0178;
    const ccomplex_t IT_0180 = 0.166666666666667*IT_0142 + (-0.5)*IT_0143 + 
      -IT_0152 + -IT_0153 + -IT_0154 + -IT_0155 + IT_0179;
    const ccomplex_t IT_0181 = cpow(IT_0004, 2);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_sG*IT_0181;
    const ccomplex_t IT_0183 = IT_0009*IT_0182;
    const ccomplex_t IT_0184 = IT_0003*IT_0182;
    const ccomplex_t IT_0185 = IT_0127*IT_0178;
    const ccomplex_t IT_0186 = cpow(IT_0099, 2);
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*m_N_4*IT_0186;
    const ccomplex_t IT_0188 = IT_0130*IT_0187;
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0183 + (-0.5)*IT_0184 + 
      -IT_0185 + -IT_0188;
    const ccomplex_t IT_0190 = IT_0063*IT_0171;
    const ccomplex_t IT_0191 = IT_0040*IT_0174;
    const ccomplex_t IT_0192 = IT_0109*IT_0187;
    const ccomplex_t IT_0193 = (-0.5)*IT_0183 + 0.166666666666667*IT_0184 + 
      -IT_0190 + -IT_0191 + -IT_0192;
    const ccomplex_t IT_0194 = 6*s_34;
    const ccomplex_t IT_0195 = 18*s_34;
    const ccomplex_t IT_0196 = (-18)*s_34;
    const ccomplex_t IT_0197 = (-6)*s_34;
    const ccomplex_t IT_0198 = m_t*s_23;
    const ccomplex_t IT_0199 = IT_0176 + 3*IT_0180 + -IT_0189 + (-3)*IT_0193;
    const ccomplex_t IT_0200 = (-3)*conj(IT_0193);
    const ccomplex_t IT_0201 = -conj(IT_0189);
    const ccomplex_t IT_0202 = conj(IT_0176) + 3*conj(IT_0180) + IT_0200 +
       IT_0201;
    const ccomplex_t IT_0203 = 6*IT_0165;
    const ccomplex_t IT_0204 = 18*IT_0169;
    const ccomplex_t IT_0205 = (-18)*IT_0156;
    const ccomplex_t IT_0206 = (-6)*IT_0148;
    const ccomplex_t IT_0207 = IT_0203 + IT_0204 + IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = 6*conj(IT_0165);
    const ccomplex_t IT_0209 = 18*conj(IT_0169);
    const ccomplex_t IT_0210 = (-18)*conj(IT_0156);
    const ccomplex_t IT_0211 = (-6)*conj(IT_0148);
    const ccomplex_t IT_0212 = IT_0208 + IT_0209 + IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = 6*IT_0169;
    const ccomplex_t IT_0214 = 18*IT_0165;
    const ccomplex_t IT_0215 = (-18)*IT_0148;
    const ccomplex_t IT_0216 = (-6)*IT_0156;
    const ccomplex_t IT_0217 = 6*conj(IT_0169);
    const ccomplex_t IT_0218 = 18*conj(IT_0165);
    const ccomplex_t IT_0219 = (-18)*conj(IT_0148);
    const ccomplex_t IT_0220 = (-6)*conj(IT_0156);
    const ccomplex_t IT_0221 = conj(IT_0112)*IT_0199 + IT_0112*IT_0202 + (
      -0.166666666666667)*conj(IT_0112)*IT_0207 + (-0.166666666666667)*IT_0112
      *IT_0212 + (-0.166666666666667)*conj(IT_0133)*((-18)*IT_0176 + (-6)
      *IT_0180 + 18*IT_0189 + 6*IT_0193 + IT_0213 + IT_0214 + IT_0215 + IT_0216)
       + (-0.166666666666667)*IT_0133*((-18)*conj(IT_0176) + (-6)*conj(IT_0180) 
      + 18*conj(IT_0189) + 6*conj(IT_0193) + IT_0217 + IT_0218 + IT_0219 +
       IT_0220);
    const ccomplex_t IT_0222 = 6*IT_0148;
    const ccomplex_t IT_0223 = 18*IT_0156;
    const ccomplex_t IT_0224 = (-18)*IT_0169;
    const ccomplex_t IT_0225 = (-6)*IT_0165;
    const ccomplex_t IT_0226 = IT_0222 + IT_0223 + IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = 6*conj(IT_0148);
    const ccomplex_t IT_0228 = 18*conj(IT_0156);
    const ccomplex_t IT_0229 = (-18)*conj(IT_0169);
    const ccomplex_t IT_0230 = (-6)*conj(IT_0165);
    const ccomplex_t IT_0231 = IT_0227 + IT_0228 + IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = 6*IT_0156;
    const ccomplex_t IT_0233 = 18*IT_0148;
    const ccomplex_t IT_0234 = (-18)*IT_0165;
    const ccomplex_t IT_0235 = (-6)*IT_0169;
    const ccomplex_t IT_0236 = 6*conj(IT_0156);
    const ccomplex_t IT_0237 = 18*conj(IT_0148);
    const ccomplex_t IT_0238 = (-18)*conj(IT_0165);
    const ccomplex_t IT_0239 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0240 = conj(IT_0180)*IT_0207 + IT_0180*IT_0212 + conj
      (IT_0176)*(IT_0213 + IT_0214 + IT_0215 + IT_0216) + IT_0176*(IT_0217 +
       IT_0218 + IT_0219 + IT_0220) + conj(IT_0193)*IT_0226 + IT_0193*IT_0231 +
       conj(IT_0189)*(IT_0232 + IT_0233 + IT_0234 + IT_0235) + IT_0189*(IT_0236 
      + IT_0237 + IT_0238 + IT_0239);
    const ccomplex_t IT_0241 = m_t*s_24;
    const ccomplex_t IT_0242 = conj(IT_0112)*IT_0199 + IT_0112*IT_0202 +
       0.166666666666667*conj(IT_0112)*IT_0226 + 0.166666666666667*IT_0112
      *IT_0231 + 0.166666666666667*conj(IT_0133)*(18*IT_0176 + 6*IT_0180 + (-18)
      *IT_0189 + (-6)*IT_0193 + IT_0232 + IT_0233 + IT_0234 + IT_0235) +
       0.166666666666667*IT_0133*(18*conj(IT_0176) + 6*conj(IT_0180) + (-18)
      *conj(IT_0189) + (-6)*conj(IT_0193) + IT_0236 + IT_0237 + IT_0238 +
       IT_0239);
    const ccomplex_t IT_0243 = IT_0117*IT_0133*conj(IT_0133);
    const ccomplex_t IT_0244 = 2*IT_0112*conj(IT_0112)*IT_0117 + 2*IT_0120*
      (conj(IT_0112)*IT_0133 + 1./2*conj(IT_0133)*IT_0134) + ((-36)*IT_0112*conj
      (IT_0112) + IT_0133*((-12)*conj(IT_0112) + (-36)*conj(IT_0133)) + (-6)
      *conj(IT_0133)*IT_0134)*IT_0135 + (conj(IT_0148)*IT_0156 + IT_0148*conj
      (IT_0156) + conj(IT_0165)*IT_0169 + IT_0165*conj(IT_0169) + conj(IT_0176)
      *IT_0180 + IT_0176*conj(IT_0180) + conj(IT_0189)*IT_0193 + IT_0189*conj
      (IT_0193))*IT_0194 + (IT_0148*conj(IT_0148) + IT_0156*conj(IT_0156) +
       IT_0165*conj(IT_0165) + IT_0169*conj(IT_0169) + IT_0176*conj(IT_0176) +
       IT_0180*conj(IT_0180) + IT_0189*conj(IT_0189) + IT_0193*conj(IT_0193))
      *IT_0195 + (conj(IT_0148)*IT_0165 + IT_0148*conj(IT_0165) + conj(IT_0156)
      *IT_0169 + IT_0156*conj(IT_0169) + conj(IT_0176)*IT_0189 + IT_0176*conj
      (IT_0189) + conj(IT_0180)*IT_0193 + IT_0180*conj(IT_0193))*IT_0196 + (conj
      (IT_0156)*IT_0165 + IT_0156*conj(IT_0165) + conj(IT_0148)*IT_0169 +
       IT_0148*conj(IT_0169) + conj(IT_0180)*IT_0189 + IT_0180*conj(IT_0189) +
       conj(IT_0176)*IT_0193 + IT_0176*conj(IT_0193))*IT_0197 + (-6)*IT_0198
      *IT_0221 + IT_0000*IT_0240 + 6*IT_0241*IT_0242 + 2*IT_0243;
    return create_ccomplex_return(IT_0244);
}

