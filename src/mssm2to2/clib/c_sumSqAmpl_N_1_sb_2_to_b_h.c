#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sb_2_to_b_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sb_2_to_b_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = m_b*N_d1*e_em*IT_0013*U_sb_00;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0007 + 1.5*IT_0016);
    const ccomplex_t IT_0018 = (-0.333333333333333)*IT_0017;
    const ccomplex_t IT_0019 = cpow(V_tb, 2);
    const ccomplex_t IT_0020 = sin(alpha);
    const ccomplex_t IT_0021 = m_b*mu_h*IT_0013;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = IT_0008*IT_0020;
    const ccomplex_t IT_0024 = sin(beta);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = IT_0023 + IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0003, -2);
    const ccomplex_t IT_0028 = 2*conj(U_sb_00)*((IT_0000*IT_0009*IT_0011
      *IT_0013*IT_0020 + (-0.5)*m_W*IT_0026*(IT_0011 + 1./3*IT_0010*IT_0027))
      *U_sb_01 + 1./2*IT_0009*IT_0011*IT_0021*IT_0022*U_sb_11) + conj(U_sb_10)*
      (IT_0009*IT_0011*IT_0021*IT_0022*U_sb_01 + 2*(IT_0000*IT_0009*IT_0011
      *IT_0013*IT_0020 + -1./3*m_W*IT_0010*IT_0026*IT_0027)*U_sb_11);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(e_em*IT_0028 +
       1.4142135623731*A_b*IT_0019*IT_0020*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0018*IT_0030;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_sb_2, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0034 
      + -reg_prop, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_b*e_em*IT_0009*IT_0011
      *IT_0013*IT_0020;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = m_b*N_d1*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0042 = IT_0012*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + 1.5*IT_0043);
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = IT_0037*IT_0045;
    const ccomplex_t IT_0047 = IT_0035*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = m_b*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0051 = 2*U_sb_01*((IT_0000*IT_0009*IT_0011*IT_0013
      *IT_0020 + (-0.5)*m_W*IT_0026*(IT_0011 + 1./3*IT_0010*IT_0027))*conj
      (U_sb_01) + 1./2*IT_0009*IT_0011*IT_0021*IT_0022*conj(U_sb_11)) + U_sb_11*
      (IT_0009*IT_0011*IT_0021*IT_0022*conj(U_sb_01) + 2*(IT_0000*IT_0009
      *IT_0011*IT_0013*IT_0020 + -1./3*m_W*IT_0010*IT_0026*IT_0027)*conj(U_sb_11
      ));
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(e_em*IT_0051 +
       1.4142135623731*A_b*IT_0019*IT_0020*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0045*IT_0053;
    const ccomplex_t IT_0055 = IT_0050*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0045*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0061 = IT_0004*IT_0022;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0064 = IT_0004*IT_0020;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0067 = IT_0011*IT_0022;
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0070 = IT_0011*IT_0020;
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0062 + IT_0065 + 
      -IT_0068 + -IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = m_N_1*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0076 = IT_0061*IT_0075;
    const ccomplex_t IT_0077 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0078 = IT_0064*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0080 = IT_0067*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0082 = IT_0070*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + IT_0078 + 
      -IT_0080 + -IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = m_N_1*IT_0084;
    const ccomplex_t IT_0086 = IT_0074 + IT_0085;
    const ccomplex_t IT_0087 = IT_0059*IT_0086;
    const ccomplex_t IT_0088 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0089 = IT_0004*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = m_b*N_d2*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0092 = IT_0012*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + 1.5*IT_0093);
    const ccomplex_t IT_0095 = (-0.333333333333333)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0100 = IT_0061*IT_0099;
    const ccomplex_t IT_0101 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0102 = IT_0064*IT_0101;
    const ccomplex_t IT_0103 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0104 = IT_0061*IT_0103;
    const ccomplex_t IT_0105 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0106 = IT_0064*IT_0105;
    const ccomplex_t IT_0107 = IT_0100 + IT_0102 + IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0109 = IT_0067*IT_0108;
    const ccomplex_t IT_0110 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0111 = IT_0070*IT_0110;
    const ccomplex_t IT_0112 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0113 = IT_0067*IT_0112;
    const ccomplex_t IT_0114 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0115 = IT_0070*IT_0114;
    const ccomplex_t IT_0116 = -IT_0109 + -IT_0111 + -IT_0113 + -IT_0115;
    const ccomplex_t IT_0117 = IT_0107 + IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = m_N_2*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0122 = IT_0061*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0124 = IT_0064*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0126 = IT_0061*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0128 = IT_0064*IT_0127;
    const ccomplex_t IT_0129 = IT_0122 + IT_0124 + IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0131 = IT_0067*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0133 = IT_0070*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0135 = IT_0067*IT_0134;
    const ccomplex_t IT_0136 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0137 = IT_0070*IT_0136;
    const ccomplex_t IT_0138 = -IT_0131 + -IT_0133 + -IT_0135 + -IT_0137;
    const ccomplex_t IT_0139 = IT_0129 + IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = m_N_1*IT_0141;
    const ccomplex_t IT_0143 = IT_0120 + IT_0142;
    const ccomplex_t IT_0144 = IT_0098*IT_0143;
    const ccomplex_t IT_0145 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0146 = IT_0061*IT_0145;
    const ccomplex_t IT_0147 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0148 = IT_0064*IT_0147;
    const ccomplex_t IT_0149 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0150 = IT_0061*IT_0149;
    const ccomplex_t IT_0151 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0152 = IT_0064*IT_0151;
    const ccomplex_t IT_0153 = IT_0146 + IT_0148 + IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0155 = IT_0067*IT_0154;
    const ccomplex_t IT_0156 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0157 = IT_0070*IT_0156;
    const ccomplex_t IT_0158 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0159 = IT_0067*IT_0158;
    const ccomplex_t IT_0160 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0161 = IT_0070*IT_0160;
    const ccomplex_t IT_0162 = -IT_0155 + -IT_0157 + -IT_0159 + -IT_0161;
    const ccomplex_t IT_0163 = IT_0153 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = m_N_3*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0168 = IT_0061*IT_0167;
    const ccomplex_t IT_0169 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0170 = IT_0064*IT_0169;
    const ccomplex_t IT_0171 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0172 = IT_0061*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0174 = IT_0064*IT_0173;
    const ccomplex_t IT_0175 = IT_0168 + IT_0170 + IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0177 = IT_0067*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0179 = IT_0070*IT_0178;
    const ccomplex_t IT_0180 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0181 = IT_0067*IT_0180;
    const ccomplex_t IT_0182 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0183 = IT_0070*IT_0182;
    const ccomplex_t IT_0184 = -IT_0177 + -IT_0179 + -IT_0181 + -IT_0183;
    const ccomplex_t IT_0185 = IT_0175 + IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = m_N_1*IT_0187;
    const ccomplex_t IT_0189 = IT_0166 + IT_0188;
    const ccomplex_t IT_0190 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0191 = IT_0004*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = m_b*N_d3*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0194 = IT_0012*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0192 + 1.5*IT_0195);
    const ccomplex_t IT_0197 = (-0.333333333333333)*IT_0196;
    const ccomplex_t IT_0198 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0199 = IT_0197*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0189*IT_0200;
    const ccomplex_t IT_0202 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0203 = IT_0061*IT_0202;
    const ccomplex_t IT_0204 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0205 = IT_0064*IT_0204;
    const ccomplex_t IT_0206 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0207 = IT_0061*IT_0206;
    const ccomplex_t IT_0208 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0209 = IT_0064*IT_0208;
    const ccomplex_t IT_0210 = IT_0203 + IT_0205 + IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0212 = IT_0067*IT_0211;
    const ccomplex_t IT_0213 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0214 = IT_0070*IT_0213;
    const ccomplex_t IT_0215 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0216 = IT_0067*IT_0215;
    const ccomplex_t IT_0217 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0218 = IT_0070*IT_0217;
    const ccomplex_t IT_0219 = -IT_0212 + -IT_0214 + -IT_0216 + -IT_0218;
    const ccomplex_t IT_0220 = IT_0210 + IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = m_N_4*IT_0222;
    const ccomplex_t IT_0224 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0225 = IT_0061*IT_0224;
    const ccomplex_t IT_0226 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0227 = IT_0064*IT_0226;
    const ccomplex_t IT_0228 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0229 = IT_0061*IT_0228;
    const ccomplex_t IT_0230 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0231 = IT_0064*IT_0230;
    const ccomplex_t IT_0232 = IT_0225 + IT_0227 + IT_0229 + IT_0231;
    const ccomplex_t IT_0233 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0234 = IT_0067*IT_0233;
    const ccomplex_t IT_0235 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0236 = IT_0070*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0238 = IT_0067*IT_0237;
    const ccomplex_t IT_0239 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0240 = IT_0070*IT_0239;
    const ccomplex_t IT_0241 = -IT_0234 + -IT_0236 + -IT_0238 + -IT_0240;
    const ccomplex_t IT_0242 = IT_0232 + IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = m_N_1*IT_0244;
    const ccomplex_t IT_0246 = IT_0223 + IT_0245;
    const ccomplex_t IT_0247 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0248 = IT_0004*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = m_b*N_d4*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0251 = IT_0012*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0249 + 1.5*IT_0252);
    const ccomplex_t IT_0254 = (-0.333333333333333)*IT_0253;
    const ccomplex_t IT_0255 = cpow((-2)*s_23 + IT_0000 + IT_0034 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0256 = IT_0254*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0246*IT_0257;
    const ccomplex_t IT_0259 = IT_0033 + (-2)*IT_0049 + IT_0056 + IT_0087 +
       IT_0144 + IT_0201 + IT_0258;
    const ccomplex_t IT_0260 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0261 = IT_0011*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = m_b*conj(N_d2)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0264 = IT_0012*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0267 = IT_0004*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*(IT_0262 + -IT_0265 + (
      -0.333333333333333)*IT_0268);
    const ccomplex_t IT_0270 = (-3)*IT_0269;
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0270;
    const ccomplex_t IT_0272 = IT_0096*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = m_N_2*IT_0141;
    const ccomplex_t IT_0275 = m_N_1*IT_0119;
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = IT_0273*IT_0276;
    const ccomplex_t IT_0278 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0279 = IT_0004*IT_0278;
    const ccomplex_t IT_0280 = 1.4142135623731*IT_0279;
    const ccomplex_t IT_0281 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0282 = IT_0011*IT_0281;
    const ccomplex_t IT_0283 = 1.4142135623731*IT_0282;
    const ccomplex_t IT_0284 = m_b*conj(N_d1)*e_em*IT_0013*U_sb_10;
    const ccomplex_t IT_0285 = IT_0012*IT_0284;
    const ccomplex_t IT_0286 = 1.4142135623731*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*(IT_0280 + (-3)*IT_0283 + 3
      *IT_0286);
    const ccomplex_t IT_0288 = (-0.333333333333333)*IT_0287;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = 3*IT_0289;
    const ccomplex_t IT_0291 = 0.166666666666667*IT_0290;
    const ccomplex_t IT_0292 = 2*conj(U_sb_00)*((IT_0000*IT_0009*IT_0011
      *IT_0013*IT_0020 + (-0.5)*m_W*IT_0026*(IT_0011 + 1./3*IT_0010*IT_0027))
      *U_sb_01 + 1./2*IT_0009*IT_0011*IT_0021*IT_0022*U_sb_11) + conj(U_sb_10)*
      (IT_0009*IT_0011*IT_0021*IT_0022*U_sb_01 + 2*(IT_0000*IT_0009*IT_0011
      *IT_0013*IT_0020 + -1./3*m_W*IT_0010*IT_0026*IT_0027)*U_sb_11);
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*(e_em*IT_0292 +
       1.4142135623731*A_b*IT_0019*IT_0020*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0294 = 0.5*IT_0293;
    const ccomplex_t IT_0295 = IT_0291*IT_0294;
    const ccomplex_t IT_0296 = IT_0002*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = m_b*conj(N_d1)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0299 = IT_0012*IT_0298;
    const ccomplex_t IT_0300 = 1.4142135623731*IT_0299;
    const ccomplex_t IT_0301 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0302 = IT_0011*IT_0301;
    const ccomplex_t IT_0303 = 1.4142135623731*IT_0302;
    const ccomplex_t IT_0304 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0305 = IT_0004*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*(IT_0300 + -IT_0303 +
       0.333333333333333*IT_0306);
    const ccomplex_t IT_0308 = 3*IT_0307;
    const ccomplex_t IT_0309 = 0.166666666666667*IT_0308;
    const ccomplex_t IT_0310 = IT_0037*IT_0309;
    const ccomplex_t IT_0311 = IT_0035*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = m_b*IT_0312;
    const ccomplex_t IT_0314 = 2*U_sb_01*((IT_0000*IT_0009*IT_0011*IT_0013
      *IT_0020 + (-0.5)*m_W*IT_0026*(IT_0011 + 1./3*IT_0010*IT_0027))*conj
      (U_sb_01) + 1./2*IT_0009*IT_0011*IT_0021*IT_0022*conj(U_sb_11)) + U_sb_11*
      (IT_0009*IT_0011*IT_0021*IT_0022*conj(U_sb_01) + 2*(IT_0000*IT_0009
      *IT_0011*IT_0013*IT_0020 + -1./3*m_W*IT_0010*IT_0026*IT_0027)*conj(U_sb_11
      ));
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*(e_em*IT_0314 +
       1.4142135623731*A_b*IT_0019*IT_0020*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0316 = 0.5*IT_0315;
    const ccomplex_t IT_0317 = IT_0309*IT_0316;
    const ccomplex_t IT_0318 = IT_0050*IT_0317;
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = IT_0057*IT_0309;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = IT_0086*IT_0321;
    const ccomplex_t IT_0323 = m_N_1*IT_0165;
    const ccomplex_t IT_0324 = m_N_3*IT_0187;
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = m_b*conj(N_d3)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0327 = IT_0012*IT_0326;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0330 = IT_0011*IT_0329;
    const ccomplex_t IT_0331 = 1.4142135623731*IT_0330;
    const ccomplex_t IT_0332 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0333 = IT_0004*IT_0332;
    const ccomplex_t IT_0334 = 1.4142135623731*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*(IT_0328 + -IT_0331 +
       0.333333333333333*IT_0334);
    const ccomplex_t IT_0336 = 3*IT_0335;
    const ccomplex_t IT_0337 = 0.166666666666667*IT_0336;
    const ccomplex_t IT_0338 = IT_0198*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = IT_0325*IT_0339;
    const ccomplex_t IT_0341 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0342 = IT_0004*IT_0341;
    const ccomplex_t IT_0343 = 1.4142135623731*IT_0342;
    const ccomplex_t IT_0344 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0345 = IT_0011*IT_0344;
    const ccomplex_t IT_0346 = 1.4142135623731*IT_0345;
    const ccomplex_t IT_0347 = m_b*conj(N_d4)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0348 = IT_0012*IT_0347;
    const ccomplex_t IT_0349 = 1.4142135623731*IT_0348;
    const ccomplex_t IT_0350 = (0 + _Complex_I*1)*(IT_0343 + (-3)*IT_0346 + 3
      *IT_0349);
    const ccomplex_t IT_0351 = 0.166666666666667*IT_0350;
    const ccomplex_t IT_0352 = IT_0255*IT_0351;
    const ccomplex_t IT_0353 = (0 + _Complex_I*1)*IT_0352;
    const ccomplex_t IT_0354 = m_N_1*IT_0222;
    const ccomplex_t IT_0355 = m_N_4*IT_0244;
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = IT_0353*IT_0356;
    const ccomplex_t IT_0358 = -IT_0277 + -IT_0297 + 2*IT_0313 + -IT_0319 + 
      -IT_0322 + -IT_0340 + -IT_0357;
    const ccomplex_t IT_0359 = 6*s_13;
    const ccomplex_t IT_0360 = m_b*m_N_1;
    const ccomplex_t IT_0361 = 6*IT_0259;
    const ccomplex_t IT_0362 = IT_0165*IT_0339;
    const ccomplex_t IT_0363 = IT_0073*IT_0321;
    const ccomplex_t IT_0364 = IT_0119*IT_0273;
    const ccomplex_t IT_0365 = IT_0222*IT_0353;
    const ccomplex_t IT_0366 = -IT_0048 + IT_0362 + IT_0363 + IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = m_b*s_14;
    const ccomplex_t IT_0368 = IT_0059*IT_0084;
    const ccomplex_t IT_0369 = IT_0187*IT_0200;
    const ccomplex_t IT_0370 = IT_0098*IT_0141;
    const ccomplex_t IT_0371 = IT_0244*IT_0257;
    const ccomplex_t IT_0372 = IT_0312 + -IT_0368 + -IT_0369 + -IT_0370 + 
      -IT_0371;
    const ccomplex_t IT_0373 = s_34*m_N_1;
    const ccomplex_t IT_0374 = 6*conj(IT_0259);
    const ccomplex_t IT_0375 = 6*conj(IT_0358);
    const ccomplex_t IT_0376 = pow(m_h, 2);
    const ccomplex_t IT_0377 = m_b*m_N_1*IT_0376;
    const ccomplex_t IT_0378 = s_13*IT_0376;
    const ccomplex_t IT_0379 = s_14*s_34;
    const ccomplex_t IT_0380 = (-2)*IT_0379;
    const ccomplex_t IT_0381 = IT_0378 + IT_0380;
    const ccomplex_t IT_0382 = 6*IT_0358;
    const ccomplex_t IT_0383 = 0.166666666666667*IT_0373;
    const ccomplex_t IT_0384 = (IT_0259*conj(IT_0259) + IT_0358*conj(IT_0358))
      *IT_0359 + conj(IT_0358)*IT_0360*IT_0361 + (IT_0358*IT_0360 + IT_0366
      *IT_0367 + IT_0372*IT_0373)*IT_0374 + (IT_0367*IT_0372 + IT_0366*IT_0373)
      *IT_0375 + 6*conj(IT_0372)*(IT_0366*IT_0377 + -IT_0372*IT_0381 +
       0.166666666666667*IT_0367*IT_0382 + IT_0361*IT_0383) + 6*conj(IT_0366)*
      (0.166666666666667*IT_0361*IT_0367 + IT_0372*IT_0377 + -IT_0366*IT_0381 +
       IT_0382*IT_0383);
    return create_ccomplex_return(IT_0384);
}

