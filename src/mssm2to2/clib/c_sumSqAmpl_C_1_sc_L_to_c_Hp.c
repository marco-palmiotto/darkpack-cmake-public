#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sc_L_to_c_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sc_L_to_c_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sc_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W2)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + 3*IT_0012);
    const ccomplex_t IT_0014 = 0.166666666666667*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = IT_0004*IT_0015;
    const ccomplex_t IT_0017 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0015;
    const ccomplex_t IT_0021 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0019 + IT_0023;
    const ccomplex_t IT_0025 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = 2*IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B1)*e_em;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W1)*e_em;
    const ccomplex_t IT_0038 = IT_0009*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 3*IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0043 = IT_0020*IT_0042;
    const ccomplex_t IT_0044 = 2*IT_0043;
    const ccomplex_t IT_0045 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0046 = IT_0016*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0049 = IT_0020*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0050;
    const ccomplex_t IT_0052 = IT_0044 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0041*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B3)*e_em;
    const ccomplex_t IT_0060 = IT_0004*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0068 = IT_0020*IT_0067;
    const ccomplex_t IT_0069 = 2*IT_0068;
    const ccomplex_t IT_0070 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0071 = IT_0016*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0074 = IT_0020*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = IT_0072 + IT_0075;
    const ccomplex_t IT_0077 = IT_0069 + IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0066*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B4)*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = conj(N_W4)*e_em;
    const ccomplex_t IT_0088 = IT_0009*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + 3*IT_0089);
    const ccomplex_t IT_0091 = 0.166666666666667*IT_0090;
    const ccomplex_t IT_0092 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0093 = IT_0020*IT_0092;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0096 = IT_0016*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0099 = IT_0020*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0033 + -IT_0058 + -IT_0083 + -IT_0108;
    const ccomplex_t IT_0110 = m_c*m_C_1*IT_0001;
    const ccomplex_t IT_0111 = pow(m_W, -1);
    const ccomplex_t IT_0112 = sin(beta);
    const ccomplex_t IT_0113 = cpow(IT_0112, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0004*IT_0112;
    const ccomplex_t IT_0117 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = IT_0009*IT_0112;
    const ccomplex_t IT_0121 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = IT_0119 + IT_0123;
    const ccomplex_t IT_0125 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0126 = IT_0120*IT_0125;
    const ccomplex_t IT_0127 = (-2)*IT_0126;
    const ccomplex_t IT_0128 = IT_0124 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0115*IT_0130;
    const ccomplex_t IT_0132 = IT_0056*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0137 = IT_0120*IT_0136;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0140 = IT_0116*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0143 = IT_0120*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0141 + IT_0144;
    const ccomplex_t IT_0146 = IT_0138 + IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0135*IT_0148;
    const ccomplex_t IT_0150 = IT_0106*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0155 = IT_0120*IT_0154;
    const ccomplex_t IT_0156 = (-2)*IT_0155;
    const ccomplex_t IT_0157 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0158 = IT_0116*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0161 = IT_0120*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = IT_0159 + IT_0162;
    const ccomplex_t IT_0164 = IT_0156 + IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0153*IT_0166;
    const ccomplex_t IT_0168 = IT_0002*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0173 = IT_0120*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0176 = IT_0116*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = IT_0174 + IT_0177;
    const ccomplex_t IT_0179 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0180 = IT_0120*IT_0179;
    const ccomplex_t IT_0181 = (-2)*IT_0180;
    const ccomplex_t IT_0182 = IT_0178 + IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0171*IT_0184;
    const ccomplex_t IT_0186 = IT_0081*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = -IT_0133 + -IT_0151 + -IT_0169 + -IT_0187;
    const ccomplex_t IT_0189 = m_c*s_12;
    const ccomplex_t IT_0190 = 6*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0192 = IT_0030*IT_0153;
    const ccomplex_t IT_0193 = IT_0002*IT_0191*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0195 = IT_0079*IT_0171;
    const ccomplex_t IT_0196 = IT_0081*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0198 = IT_0104*IT_0135;
    const ccomplex_t IT_0199 = IT_0106*IT_0197*IT_0198;
    const ccomplex_t IT_0200 = pow(m_C_1, 2);
    const ccomplex_t IT_0201 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0202 = cpow(IT_0015, -1);
    const ccomplex_t IT_0203 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0009*(m_W*IT_0015*IT_0112*conj(U_sb_00) + (-0.5)*IT_0111*(IT_0000
      *IT_0015*IT_0113*conj(U_sb_00) + m_b*(m_b*IT_0112*IT_0202*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0009*IT_0111*IT_0113*U_sb_00;
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = IT_0204*IT_0206;
    const ccomplex_t IT_0208 = IT_0201*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0211 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0009*(m_W*IT_0015*IT_0112*conj(U_sb_01) + (-0.5)*IT_0111*(IT_0000
      *IT_0015*IT_0113*conj(U_sb_01) + m_b*(m_b*IT_0112*IT_0202*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0009*IT_0111*IT_0113*U_sb_01;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = IT_0210*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = m_c*IT_0058;
    const ccomplex_t IT_0219 = pow(m_d, 2);
    const ccomplex_t IT_0220 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cd)
      *e_em*IT_0009*(IT_0015*(m_W*IT_0112 + (-0.5)*IT_0000*IT_0111*IT_0113) + 
      -1./2*IT_0111*IT_0112*IT_0202*IT_0219);
    const ccomplex_t IT_0221 = 0.5*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u1
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0223 = 0.5*IT_0222;
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0226 = IT_0224*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = m_c*IT_0033;
    const ccomplex_t IT_0229 = m_c*IT_0083;
    const ccomplex_t IT_0230 = m_c*IT_0108;
    const ccomplex_t IT_0231 = pow(m_s, 2);
    const ccomplex_t IT_0232 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cs)
      *e_em*IT_0009*(IT_0015*(m_W*IT_0112 + (-0.5)*IT_0000*IT_0111*IT_0113) + 
      -1./2*IT_0111*IT_0112*IT_0202*IT_0231);
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0009*IT_0111*IT_0113;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = IT_0233*IT_0235;
    const ccomplex_t IT_0237 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0238 = IT_0236*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0241 = IT_0054*IT_0115;
    const ccomplex_t IT_0242 = IT_0056*IT_0240*IT_0241;
    const ccomplex_t IT_0243 = -IT_0193 + -IT_0196 + -IT_0199 + IT_0209 +
       IT_0217 + IT_0218 + IT_0227 + IT_0228 + IT_0229 + IT_0230 + IT_0239 + 
      -IT_0242;
    const ccomplex_t IT_0244 = 0.166666666666667*conj(IT_0243);
    const ccomplex_t IT_0245 = s_23*m_C_1;
    const ccomplex_t IT_0246 = 6*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*mu_h*IT_0009*IT_0111;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_cd*e_em*IT_0009*IT_0111*IT_0202;
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = IT_0248*IT_0250;
    const ccomplex_t IT_0252 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0253 = IT_0251*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0009*IT_0111;
    const ccomplex_t IT_0256 = 0.5*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0009*IT_0111*IT_0202;
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0256*IT_0258;
    const ccomplex_t IT_0260 = cpow((-2)*s_13 + IT_0000 + IT_0200 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0261 = IT_0259*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*V_cd*e_em*conj(U_Wm1)*IT_0009;
    const ccomplex_t IT_0264 = IT_0221*IT_0263;
    const ccomplex_t IT_0265 = IT_0225*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0041*IT_0130;
    const ccomplex_t IT_0268 = IT_0056*IT_0240*IT_0267;
    const ccomplex_t IT_0269 = IT_0066*IT_0184;
    const ccomplex_t IT_0270 = IT_0081*IT_0194*IT_0269;
    const ccomplex_t IT_0271 = IT_0091*IT_0148;
    const ccomplex_t IT_0272 = IT_0106*IT_0197*IT_0271;
    const ccomplex_t IT_0273 = V_cb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0274 = IT_0009*IT_0273;
    const ccomplex_t IT_0275 = IT_0009*IT_0202;
    const ccomplex_t IT_0276 = m_b*conj(U_d1)*V_cb*e_em*IT_0111*U_sb_10;
    const ccomplex_t IT_0277 = IT_0275*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0274 + (-0.5)*IT_0278);
    const ccomplex_t IT_0280 = IT_0204*IT_0279;
    const ccomplex_t IT_0281 = IT_0201*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = m_c*IT_0133;
    const ccomplex_t IT_0284 = m_c*IT_0187;
    const ccomplex_t IT_0285 = IT_0014*IT_0166;
    const ccomplex_t IT_0286 = IT_0002*IT_0191*IT_0285;
    const ccomplex_t IT_0287 = V_cb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0288 = IT_0009*IT_0287;
    const ccomplex_t IT_0289 = m_b*conj(U_d1)*V_cb*e_em*IT_0111*U_sb_11;
    const ccomplex_t IT_0290 = IT_0275*IT_0289;
    const ccomplex_t IT_0291 = 1.4142135623731*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*(IT_0288 + (-0.5)*IT_0291);
    const ccomplex_t IT_0293 = IT_0212*IT_0292;
    const ccomplex_t IT_0294 = IT_0210*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = m_c*IT_0169;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0009;
    const ccomplex_t IT_0298 = IT_0233*IT_0297;
    const ccomplex_t IT_0299 = IT_0237*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = m_c*IT_0151;
    const ccomplex_t IT_0302 = -IT_0254 + -IT_0262 + -IT_0266 + -IT_0268 + 
      -IT_0270 + -IT_0272 + -IT_0282 + IT_0283 + IT_0284 + -IT_0286 + -IT_0295 +
       IT_0296 + -IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = 0.166666666666667*conj(IT_0302);
    const ccomplex_t IT_0304 = s_13*IT_0001;
    const ccomplex_t IT_0305 = (-6)*IT_0304;
    const ccomplex_t IT_0306 = s_12*s_23;
    const ccomplex_t IT_0307 = 12*IT_0306;
    const ccomplex_t IT_0308 = IT_0305 + IT_0307;
    const ccomplex_t IT_0309 = 0.166666666666667*conj(IT_0109);
    const ccomplex_t IT_0310 = 0.166666666666667*conj(IT_0188);
    const ccomplex_t IT_0311 = m_c*m_C_1;
    const ccomplex_t IT_0312 = 6*IT_0311;
    const ccomplex_t IT_0313 = 6*IT_0109*(IT_0110*conj(IT_0188) + IT_0190
      *IT_0244 + IT_0246*IT_0303 + IT_0308*IT_0309) + 6*IT_0188*(conj(IT_0109)
      *IT_0110 + IT_0244*IT_0246 + IT_0190*IT_0303 + IT_0308*IT_0310) + 6
      *IT_0302*(s_13*conj(IT_0302) + IT_0246*IT_0309 + IT_0190*IT_0310 + IT_0244
      *IT_0312) + 6*IT_0243*(s_13*conj(IT_0243) + IT_0190*IT_0309 + IT_0246
      *IT_0310 + IT_0303*IT_0312);
    return create_ccomplex_return(IT_0313);
}

