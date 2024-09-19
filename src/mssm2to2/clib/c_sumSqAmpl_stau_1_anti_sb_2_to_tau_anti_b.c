#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_sb_2_to_tau_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_sb_2_to_tau_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_1, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = m_b*N_d3*e_em*IT_0016*conj(U_sb_11);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + (-3)*IT_0012 + 3
      *IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = conj(N_d3)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0029 = IT_0015*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + IT_0027 + 
      -IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0021*IT_0032;
    const ccomplex_t IT_0034 = IT_0002*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0037 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0038 = IT_0004*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0041 = IT_0009*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_d1)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0039 + IT_0042 + 
      -IT_0045);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0049 = IT_0009*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = m_b*N_d1*e_em*IT_0016*conj(U_sb_11);
    const ccomplex_t IT_0052 = IT_0015*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0055 = IT_0004*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + -IT_0053 + (
      -0.333333333333333)*IT_0056);
    const ccomplex_t IT_0058 = (-3)*IT_0057;
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = IT_0047*IT_0059;
    const ccomplex_t IT_0061 = IT_0036*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d2*e_em*IT_0016*conj(U_sb_11);
    const ccomplex_t IT_0071 = IT_0015*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + (-3)*IT_0069 + 3
      *IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0076 = IT_0004*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0079 = IT_0009*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = conj(N_d2)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0082 = IT_0015*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0077 + IT_0080 + 
      -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = IT_0063*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0090 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0094 = IT_0009*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = m_b*N_d4*e_em*IT_0016*conj(U_sb_11);
    const ccomplex_t IT_0097 = IT_0015*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + (-3)*IT_0095 + 3
      *IT_0098);
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0102 = IT_0004*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0105 = IT_0009*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_d4)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0108 = IT_0015*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + IT_0106 + 
      -IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0111;
    const ccomplex_t IT_0113 = IT_0089*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0035 + 0.5*IT_0062 + 0.5*IT_0088 + 0.5
      *IT_0114;
    const ccomplex_t IT_0116 = m_tau*IT_0062;
    const ccomplex_t IT_0117 = m_tau*IT_0088;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0119 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0120 = IT_0004*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_d2*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0123 = IT_0015*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0074*IT_0126;
    const ccomplex_t IT_0128 = IT_0063*IT_0118*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0130 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0131 = IT_0004*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = N_d4*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0134 = IT_0015*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0100*IT_0137;
    const ccomplex_t IT_0139 = IT_0089*IT_0129*IT_0138;
    const ccomplex_t IT_0140 = m_tau*IT_0035;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0142 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0143 = IT_0004*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_d1*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0146 = IT_0015*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + 0.5*IT_0147);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = IT_0059*IT_0149;
    const ccomplex_t IT_0151 = IT_0036*IT_0141*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0153 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0154 = IT_0004*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = N_d3*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0157 = IT_0015*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + 0.5*IT_0158);
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0021*IT_0160;
    const ccomplex_t IT_0162 = IT_0002*IT_0152*IT_0161;
    const ccomplex_t IT_0163 = m_tau*IT_0114;
    const ccomplex_t IT_0164 = (-0.5)*IT_0116 + (-0.5)*IT_0117 + 0.5*IT_0128 +
       0.5*IT_0139 + (-0.5)*IT_0140 + 0.5*IT_0151 + 0.5*IT_0162 + (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = s_14*m_tau;
    const ccomplex_t IT_0166 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0167 = IT_0004*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = m_b*conj(N_d1)*e_em*IT_0016*conj(U_sb_01);
    const ccomplex_t IT_0170 = IT_0015*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + 1.5*IT_0171);
    const ccomplex_t IT_0173 = (-0.333333333333333)*IT_0172;
    const ccomplex_t IT_0174 = IT_0047*IT_0173;
    const ccomplex_t IT_0175 = IT_0036*IT_0141*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0177 = IT_0004*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = m_b*conj(N_d2)*e_em*IT_0016*conj(U_sb_01);
    const ccomplex_t IT_0180 = IT_0015*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0178 + 1.5*IT_0181);
    const ccomplex_t IT_0183 = (-0.333333333333333)*IT_0182;
    const ccomplex_t IT_0184 = IT_0085*IT_0183;
    const ccomplex_t IT_0185 = IT_0063*IT_0118*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0187 = IT_0004*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = m_b*conj(N_d4)*e_em*IT_0016*conj(U_sb_01);
    const ccomplex_t IT_0190 = IT_0015*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0188 + 1.5*IT_0191);
    const ccomplex_t IT_0193 = (-0.333333333333333)*IT_0192;
    const ccomplex_t IT_0194 = IT_0111*IT_0193;
    const ccomplex_t IT_0195 = IT_0089*IT_0129*IT_0194;
    const ccomplex_t IT_0196 = IT_0149*IT_0173;
    const ccomplex_t IT_0197 = IT_0036*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = m_tau*IT_0198;
    const ccomplex_t IT_0200 = IT_0126*IT_0183;
    const ccomplex_t IT_0201 = IT_0063*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = m_tau*IT_0202;
    const ccomplex_t IT_0204 = IT_0137*IT_0193;
    const ccomplex_t IT_0205 = IT_0089*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = m_tau*IT_0206;
    const ccomplex_t IT_0208 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0209 = IT_0004*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = m_b*conj(N_d3)*e_em*IT_0016*conj(U_sb_01);
    const ccomplex_t IT_0212 = IT_0015*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*(IT_0210 + 1.5*IT_0213);
    const ccomplex_t IT_0215 = (-0.333333333333333)*IT_0214;
    const ccomplex_t IT_0216 = IT_0032*IT_0215;
    const ccomplex_t IT_0217 = IT_0002*IT_0152*IT_0216;
    const ccomplex_t IT_0218 = IT_0160*IT_0215;
    const ccomplex_t IT_0219 = IT_0002*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = m_tau*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0175 + 0.5*IT_0185 + 0.5*IT_0195 + (-0.5
      )*IT_0199 + (-0.5)*IT_0203 + (-0.5)*IT_0207 + 0.5*IT_0217 + (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = m_b*s_13;
    const ccomplex_t IT_0224 = 0.5*IT_0198 + 0.5*IT_0202 + 0.5*IT_0206 + 0.5
      *IT_0220;
    const ccomplex_t IT_0225 = m_b*m_tau*IT_0001;
    const ccomplex_t IT_0226 = s_34*IT_0001;
    const ccomplex_t IT_0227 = (-6)*IT_0226;
    const ccomplex_t IT_0228 = s_13*s_14;
    const ccomplex_t IT_0229 = 12*IT_0228;
    const ccomplex_t IT_0230 = IT_0227 + IT_0229;
    const ccomplex_t IT_0231 = m_b*m_tau;
    const ccomplex_t IT_0232 = 6*s_34;
    const ccomplex_t IT_0233 = -conj(IT_0115);
    const ccomplex_t IT_0234 = (-6)*conj(IT_0164);
    const ccomplex_t IT_0235 = 6*IT_0115*(conj(IT_0164)*IT_0165 + -conj
      (IT_0222)*IT_0223 + -conj(IT_0224)*IT_0225 + 0.166666666666667*conj
      (IT_0115)*IT_0230) + (-6)*IT_0164*(IT_0223*conj(IT_0224) + conj(IT_0222)
      *IT_0231 + (-0.166666666666667)*conj(IT_0164)*IT_0232 + IT_0165*IT_0233) +
       6*IT_0224*(IT_0165*conj(IT_0222) + 0.166666666666667*conj(IT_0224)
      *IT_0230 + IT_0225*IT_0233 + 0.166666666666667*IT_0223*IT_0234) + 6
      *IT_0222*(IT_0165*conj(IT_0224) + 0.166666666666667*conj(IT_0222)*IT_0232 
      + IT_0223*IT_0233 + 0.166666666666667*IT_0231*IT_0234);
    return create_ccomplex_return(IT_0235);
}

