#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_sc_R_to_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_sc_R_to_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sc_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W1)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0012 + 3*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0009;
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B2)*e_em;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W2)*e_em;
    const ccomplex_t IT_0030 = IT_0014*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + 3*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0009;
    const ccomplex_t IT_0035 = 0.666666666666667*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B3)*e_em;
    const ccomplex_t IT_0041 = IT_0009*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W3)*e_em;
    const ccomplex_t IT_0044 = IT_0014*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 3*IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0009;
    const ccomplex_t IT_0049 = 0.666666666666667*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B4)*e_em;
    const ccomplex_t IT_0055 = IT_0009*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W4)*e_em;
    const ccomplex_t IT_0058 = IT_0014*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + 3*IT_0059);
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0009;
    const ccomplex_t IT_0063 = 0.666666666666667*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = 0.166666666666667*IT_0007 + -IT_0025 + -IT_0039
       + -IT_0053 + -IT_0067;
    const ccomplex_t IT_0069 = pow(m_W, -1);
    const ccomplex_t IT_0070 = sin(beta);
    const ccomplex_t IT_0071 = cpow(IT_0070, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = pow(m_sc_R, 2);
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0003 + IT_0077 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0003 + IT_0077 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0003 + IT_0077 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0014*IT_0069*IT_0071;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0003 + IT_0077 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0007 + IT_0080 + IT_0088 + IT_0096 +
       IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0077 +
       reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0002*IT_0106;
    const ccomplex_t IT_0108 = IT_0023*IT_0076;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0037*IT_0085;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0051*IT_0093;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0065*IT_0101;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0107 + -IT_0109 + -IT_0111 + -IT_0113 + 
      -IT_0115;
    const ccomplex_t IT_0117 = IT_0064*IT_0102;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0050*IT_0094;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0036*IT_0086;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0022*IT_0078;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = (-0.166666666666667)*IT_0107 + IT_0118 +
       IT_0120 + IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = s_34*IT_0077;
    const ccomplex_t IT_0127 = (-6)*IT_0126;
    const ccomplex_t IT_0128 = s_23*s_24;
    const ccomplex_t IT_0129 = 12*IT_0128;
    const ccomplex_t IT_0130 = IT_0127 + IT_0129;
    const ccomplex_t IT_0131 = (-18)*IT_0126;
    const ccomplex_t IT_0132 = 36*IT_0128;
    const ccomplex_t IT_0133 = IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = m_c*IT_0007;
    const ccomplex_t IT_0135 = m_c*IT_0118;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0137 = IT_0019*IT_0073;
    const ccomplex_t IT_0138 = IT_0078*IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0134 + -IT_0135 + IT_0138;
    const ccomplex_t IT_0140 = (-0.333333333333333)*IT_0007 + 2*IT_0025 + 2
      *IT_0039 + 2*IT_0053 + 2*IT_0067;
    const ccomplex_t IT_0141 = IT_0107 + IT_0140;
    const ccomplex_t IT_0142 = m_c*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = 6*s_34;
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = 18*s_34;
    const ccomplex_t IT_0147 = IT_0139*IT_0146;
    const ccomplex_t IT_0148 = IT_0145 + IT_0147;
    const ccomplex_t IT_0149 = IT_0003*IT_0077;
    const ccomplex_t IT_0150 = (-18)*IT_0068;
    const ccomplex_t IT_0151 = (-6)*IT_0105;
    const ccomplex_t IT_0152 = IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = (-18)*conj(IT_0068);
    const ccomplex_t IT_0154 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0155 = IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = (-18)*IT_0105;
    const ccomplex_t IT_0157 = (-6)*IT_0068;
    const ccomplex_t IT_0158 = IT_0156 + IT_0157;
    const ccomplex_t IT_0159 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0160 = (-6)*conj(IT_0068);
    const ccomplex_t IT_0161 = IT_0159 + IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0163 = IT_0033*IT_0082;
    const ccomplex_t IT_0164 = IT_0086*IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0166 = IT_0047*IT_0090;
    const ccomplex_t IT_0167 = IT_0094*IT_0165*IT_0166;
    const ccomplex_t IT_0168 = -IT_0164 + -IT_0167;
    const ccomplex_t IT_0169 = IT_0023*IT_0136*IT_0137;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0171 = IT_0061*IT_0098;
    const ccomplex_t IT_0172 = IT_0065*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = IT_0051*IT_0165*IT_0166;
    const ccomplex_t IT_0174 = -IT_0169 + -IT_0172 + -IT_0173;
    const ccomplex_t IT_0175 = IT_0035*IT_0084;
    const ccomplex_t IT_0176 = IT_0086*IT_0162*IT_0175;
    const ccomplex_t IT_0177 = IT_0049*IT_0092;
    const ccomplex_t IT_0178 = IT_0094*IT_0165*IT_0177;
    const ccomplex_t IT_0179 = -IT_0176 + -IT_0178;
    const ccomplex_t IT_0180 = IT_0021*IT_0075;
    const ccomplex_t IT_0181 = IT_0023*IT_0136*IT_0180;
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = IT_0051*IT_0165*IT_0177;
    const ccomplex_t IT_0184 = IT_0063*IT_0100;
    const ccomplex_t IT_0185 = IT_0065*IT_0170*IT_0184;
    const ccomplex_t IT_0186 = -IT_0183 + -IT_0185;
    const ccomplex_t IT_0187 = IT_0182 + IT_0186;
    const ccomplex_t IT_0188 = conj(IT_0182) + conj(IT_0186);
    const ccomplex_t IT_0189 = IT_0037*IT_0162*IT_0163;
    const ccomplex_t IT_0190 = (-6)*IT_0120 + (-6)*IT_0122 + (-6)*IT_0124;
    const ccomplex_t IT_0191 = IT_0107 + IT_0190;
    const ccomplex_t IT_0192 = m_c*IT_0191;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0192;
    const ccomplex_t IT_0194 = IT_0102*IT_0170*IT_0171;
    const ccomplex_t IT_0195 = IT_0139 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = conj(IT_0139) + conj(IT_0193) + conj(IT_0194);
    const ccomplex_t IT_0197 = conj(IT_0143)*IT_0144;
    const ccomplex_t IT_0198 = IT_0037*IT_0162*IT_0175;
    const ccomplex_t IT_0199 = m_c*IT_0109;
    const ccomplex_t IT_0200 = m_c*IT_0111;
    const ccomplex_t IT_0201 = m_c*IT_0113;
    const ccomplex_t IT_0202 = m_c*IT_0115;
    const ccomplex_t IT_0203 = IT_0198 + -IT_0199 + -IT_0200 + -IT_0201 + 
      -IT_0202;
    const ccomplex_t IT_0204 = IT_0078*IT_0136*IT_0180;
    const ccomplex_t IT_0205 = IT_0102*IT_0170*IT_0184;
    const ccomplex_t IT_0206 = m_c*IT_0080;
    const ccomplex_t IT_0207 = m_c*IT_0088;
    const ccomplex_t IT_0208 = m_c*IT_0096;
    const ccomplex_t IT_0209 = m_c*IT_0104;
    const ccomplex_t IT_0210 = IT_0205 + -IT_0206 + -IT_0207 + -IT_0208 + 
      -IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0204) + conj(IT_0210);
    const ccomplex_t IT_0212 = conj(IT_0143) + conj(IT_0189);
    const ccomplex_t IT_0213 = conj(IT_0193) + conj(IT_0194);
    const ccomplex_t IT_0214 = IT_0143*IT_0146;
    const ccomplex_t IT_0215 = IT_0143 + IT_0189;
    const ccomplex_t IT_0216 = (-18)*s_34;
    const ccomplex_t IT_0217 = (-6)*s_34;
    const ccomplex_t IT_0218 = IT_0144*IT_0203;
    const ccomplex_t IT_0219 = IT_0146*IT_0210;
    const ccomplex_t IT_0220 = IT_0179*IT_0216;
    const ccomplex_t IT_0221 = IT_0186*IT_0217;
    const ccomplex_t IT_0222 = IT_0144*conj(IT_0203);
    const ccomplex_t IT_0223 = conj(IT_0179)*IT_0216;
    const ccomplex_t IT_0224 = conj(IT_0182)*IT_0217;
    const ccomplex_t IT_0225 = conj(IT_0186)*IT_0217;
    const ccomplex_t IT_0226 = IT_0146*IT_0193;
    const ccomplex_t IT_0227 = IT_0168*IT_0216;
    const ccomplex_t IT_0228 = IT_0174*IT_0217;
    const ccomplex_t IT_0229 = conj(IT_0139)*IT_0146;
    const ccomplex_t IT_0230 = conj(IT_0168)*IT_0216;
    const ccomplex_t IT_0231 = conj(IT_0174)*IT_0217;
    const ccomplex_t IT_0232 = m_c*s_24;
    const ccomplex_t IT_0233 = 6*IT_0210;
    const ccomplex_t IT_0234 = 18*IT_0203;
    const ccomplex_t IT_0235 = (-18)*IT_0182;
    const ccomplex_t IT_0236 = (-18)*IT_0186;
    const ccomplex_t IT_0237 = (-6)*IT_0179;
    const ccomplex_t IT_0238 = 6*conj(IT_0210);
    const ccomplex_t IT_0239 = 18*conj(IT_0203);
    const ccomplex_t IT_0240 = (-18)*conj(IT_0182);
    const ccomplex_t IT_0241 = (-18)*conj(IT_0186);
    const ccomplex_t IT_0242 = (-6)*conj(IT_0179);
    const ccomplex_t IT_0243 = 6*IT_0203;
    const ccomplex_t IT_0244 = 18*IT_0210;
    const ccomplex_t IT_0245 = (-18)*IT_0179;
    const ccomplex_t IT_0246 = (-6)*IT_0182;
    const ccomplex_t IT_0247 = (-6)*IT_0186;
    const ccomplex_t IT_0248 = 6*conj(IT_0203);
    const ccomplex_t IT_0249 = 18*conj(IT_0210);
    const ccomplex_t IT_0250 = (-18)*conj(IT_0179);
    const ccomplex_t IT_0251 = (-6)*conj(IT_0182);
    const ccomplex_t IT_0252 = (-6)*conj(IT_0186);
    const ccomplex_t IT_0253 = 1*conj(IT_0116)*(IT_0139 + 3*IT_0143 + -IT_0168
       + (-3)*IT_0174 + 3*IT_0189 + IT_0193 + IT_0194) + 3*conj(IT_0125)*
      (IT_0139 + 1./3*IT_0143 + -IT_0168 + -1./3*IT_0174 + 1./3*IT_0189 +
       IT_0193 + IT_0194) + 1*IT_0116*(conj(IT_0139) + 3*conj(IT_0143) + -conj
      (IT_0168) + (-3)*conj(IT_0174) + 3*conj(IT_0189) + conj(IT_0193) + conj
      (IT_0194)) + 3*IT_0125*(conj(IT_0139) + 1./3*conj(IT_0143) + -conj(IT_0168
      ) + -1./3*conj(IT_0174) + 1./3*conj(IT_0189) + conj(IT_0193) + conj
      (IT_0194)) + conj(IT_0068)*(IT_0204 + 0.166666666666667*IT_0233 +
       0.166666666666667*IT_0234 + 0.166666666666667*IT_0235 + 0.166666666666667
      *IT_0236 + 0.166666666666667*IT_0237) + IT_0068*(conj(IT_0204) +
       0.166666666666667*IT_0238 + 0.166666666666667*IT_0239 + 0.166666666666667
      *IT_0240 + 0.166666666666667*IT_0241 + 0.166666666666667*IT_0242) +
       0.166666666666667*conj(IT_0105)*(18*IT_0204 + IT_0243 + IT_0244 + IT_0245
       + IT_0246 + IT_0247) + 0.166666666666667*IT_0105*(18*conj(IT_0204) +
       IT_0248 + IT_0249 + IT_0250 + IT_0251 + IT_0252);
    const ccomplex_t IT_0254 = m_c*s_23;
    const ccomplex_t IT_0255 = 6*IT_0182;
    const ccomplex_t IT_0256 = 6*IT_0186;
    const ccomplex_t IT_0257 = 18*IT_0179;
    const ccomplex_t IT_0258 = (-18)*IT_0210;
    const ccomplex_t IT_0259 = (-6)*IT_0203;
    const ccomplex_t IT_0260 = 6*conj(IT_0182);
    const ccomplex_t IT_0261 = 6*conj(IT_0186);
    const ccomplex_t IT_0262 = 18*conj(IT_0179);
    const ccomplex_t IT_0263 = (-18)*conj(IT_0210);
    const ccomplex_t IT_0264 = (-6)*conj(IT_0203);
    const ccomplex_t IT_0265 = 6*IT_0179;
    const ccomplex_t IT_0266 = 18*IT_0182;
    const ccomplex_t IT_0267 = 18*IT_0186;
    const ccomplex_t IT_0268 = (-18)*IT_0203;
    const ccomplex_t IT_0269 = (-6)*IT_0210;
    const ccomplex_t IT_0270 = 6*conj(IT_0179);
    const ccomplex_t IT_0271 = 18*conj(IT_0182);
    const ccomplex_t IT_0272 = 18*conj(IT_0186);
    const ccomplex_t IT_0273 = (-18)*conj(IT_0203);
    const ccomplex_t IT_0274 = (-6)*conj(IT_0210);
    const ccomplex_t IT_0275 = 6*IT_0174;
    const ccomplex_t IT_0276 = 18*IT_0168;
    const ccomplex_t IT_0277 = (-18)*IT_0139;
    const ccomplex_t IT_0278 = (-18)*IT_0193;
    const ccomplex_t IT_0279 = (-6)*IT_0143;
    const ccomplex_t IT_0280 = 6*conj(IT_0174);
    const ccomplex_t IT_0281 = 18*conj(IT_0168);
    const ccomplex_t IT_0282 = (-18)*conj(IT_0139);
    const ccomplex_t IT_0283 = (-18)*conj(IT_0193);
    const ccomplex_t IT_0284 = (-6)*conj(IT_0143);
    const ccomplex_t IT_0285 = IT_0155*IT_0189 + IT_0152*conj(IT_0189) + (-6)
      *conj(IT_0068)*(IT_0139 + 3*IT_0143 + -IT_0168 + (-3)*IT_0174 + IT_0193) +
       (-6)*IT_0068*(conj(IT_0139) + 3*conj(IT_0143) + -conj(IT_0168) + (-3)
      *conj(IT_0174) + conj(IT_0193)) + IT_0161*IT_0194 + IT_0158*conj(IT_0194) 
      + conj(IT_0125)*((-18)*IT_0204 + IT_0255 + IT_0256 + IT_0257 + IT_0258 +
       IT_0259) + IT_0125*((-18)*conj(IT_0204) + IT_0260 + IT_0261 + IT_0262 +
       IT_0263 + IT_0264) + (-6)*conj(IT_0116)*(IT_0204 + (-0.166666666666667)
      *IT_0265 + (-0.166666666666667)*IT_0266 + (-0.166666666666667)*IT_0267 + (
      -0.166666666666667)*IT_0268 + (-0.166666666666667)*IT_0269) + (-6)*IT_0116
      *(conj(IT_0204) + (-0.166666666666667)*IT_0270 + (-0.166666666666667)
      *IT_0271 + (-0.166666666666667)*IT_0272 + (-0.166666666666667)*IT_0273 + (
      -0.166666666666667)*IT_0274) + conj(IT_0105)*(IT_0275 + IT_0276 + IT_0277 
      + IT_0278 + IT_0279) + IT_0105*(IT_0280 + IT_0281 + IT_0282 + IT_0283 +
       IT_0284);
    const ccomplex_t IT_0286 = conj(IT_0174)*(IT_0233 + IT_0234 + IT_0235 +
       IT_0236 + IT_0237) + IT_0174*(IT_0238 + IT_0239 + IT_0240 + IT_0241 +
       IT_0242) + conj(IT_0168)*(IT_0243 + IT_0244 + IT_0245 + IT_0246 + IT_0247
      ) + IT_0168*(IT_0248 + IT_0249 + IT_0250 + IT_0251 + IT_0252) + IT_0196*
      (IT_0255 + IT_0256 + IT_0257 + IT_0258 + IT_0259) + IT_0195*(IT_0260 +
       IT_0261 + IT_0262 + IT_0263 + IT_0264) + IT_0212*(IT_0265 + IT_0266 +
       IT_0267 + IT_0268 + IT_0269) + IT_0215*(IT_0270 + IT_0271 + IT_0272 +
       IT_0273 + IT_0274) + (-6)*conj(IT_0204)*(IT_0189 + 3*IT_0194 + -1./6
      *IT_0275 + -1./6*IT_0276 + -1./6*IT_0277 + -1./6*IT_0278 + -1./6*IT_0279) 
      + (-6)*IT_0204*(conj(IT_0189) + 3*conj(IT_0194) + -1./6*IT_0280 + -1./6
      *IT_0281 + -1./6*IT_0282 + -1./6*IT_0283 + -1./6*IT_0284);
    const ccomplex_t IT_0287 = (conj(IT_0068)*IT_0105 + IT_0068*conj(IT_0105) 
      + conj(IT_0116)*IT_0125 + IT_0116*conj(IT_0125))*IT_0130 + (IT_0068*conj
      (IT_0068) + IT_0105*conj(IT_0105) + IT_0116*conj(IT_0116) + IT_0125*conj
      (IT_0125))*IT_0133 + conj(IT_0139)*IT_0148 + IT_0149*(conj(IT_0116)
      *IT_0152 + IT_0116*IT_0155 + conj(IT_0125)*IT_0158 + IT_0125*IT_0161) +
       IT_0144*(conj(IT_0168)*IT_0174 + IT_0168*conj(IT_0174) + conj(IT_0179)
      *IT_0187 + IT_0179*IT_0188 + conj(IT_0189)*IT_0195 + IT_0189*IT_0196) +
       IT_0139*IT_0197 + IT_0146*(IT_0168*conj(IT_0168) + IT_0174*conj(IT_0174) 
      + IT_0179*conj(IT_0179) + IT_0187*IT_0188 + IT_0203*conj(IT_0203) +
       IT_0204*IT_0211 + IT_0189*IT_0212 + IT_0194*IT_0213) + IT_0212*IT_0214 + 
      (conj(IT_0139)*IT_0168 + IT_0139*conj(IT_0168) + IT_0188*IT_0203 + IT_0187
      *conj(IT_0203) + IT_0174*IT_0212 + conj(IT_0174)*IT_0215)*IT_0216 + (conj
      (IT_0139)*IT_0174 + IT_0139*conj(IT_0174) + conj(IT_0179)*IT_0203 +
       IT_0179*conj(IT_0203) + IT_0182*IT_0211 + IT_0168*IT_0212 + conj(IT_0168)
      *IT_0215)*IT_0217 + IT_0211*(IT_0218 + IT_0219 + IT_0220 + IT_0221) + 
      (IT_0204 + IT_0210)*(IT_0222 + IT_0223 + IT_0224 + IT_0225) + IT_0213*
      (IT_0148 + IT_0226 + IT_0227 + IT_0228) + (IT_0193 + IT_0194)*(IT_0197 +
       IT_0229 + IT_0230 + IT_0231) + 6*IT_0232*IT_0253 + IT_0254*IT_0285 +
       IT_0003*IT_0286;
    return create_ccomplex_return(IT_0287);
}

