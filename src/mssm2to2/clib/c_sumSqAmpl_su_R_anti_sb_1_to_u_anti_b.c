#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sb_1_to_u_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sb_1_to_u_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_su_R, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_13*s_14;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_u, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0007 +
       reg_prop, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = cos(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = sin(theta_W);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = pow(m_W, -1);
    const ccomplex_t IT_0027 = m_b*conj(N_d1)*e_em*IT_0026*conj(U_sb_00);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0020 + 1.5*IT_0029);
    const ccomplex_t IT_0031 = (-0.333333333333333)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0017;
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0017;
    const ccomplex_t IT_0039 = 0.666666666666667*IT_0038;
    const ccomplex_t IT_0040 = m_b*conj(N_d3)*e_em*IT_0026*conj(U_sb_00);
    const ccomplex_t IT_0041 = IT_0025*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.666666666666667
      *IT_0045);
    const ccomplex_t IT_0047 = 1.5*IT_0046;
    const ccomplex_t IT_0048 = (-0.333333333333333)*IT_0047;
    const ccomplex_t IT_0049 = IT_0039*IT_0048;
    const ccomplex_t IT_0050 = IT_0037*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0017;
    const ccomplex_t IT_0054 = 0.666666666666667*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0056 = IT_0017*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = m_b*conj(N_d4)*e_em*IT_0026*conj(U_sb_00);
    const ccomplex_t IT_0059 = IT_0025*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0060);
    const ccomplex_t IT_0062 = (-0.333333333333333)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0067 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0068 = IT_0017*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_b*conj(N_d2)*e_em*IT_0026*conj(U_sb_00);
    const ccomplex_t IT_0071 = IT_0025*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 1.5*IT_0072);
    const ccomplex_t IT_0074 = (-0.333333333333333)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0017;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0066*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0036 + -IT_0051 + -IT_0065 + -IT_0079;
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0013;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = conj(IT_0080) + conj(IT_0081);
    const ccomplex_t IT_0084 = (-18)*IT_0001;
    const ccomplex_t IT_0085 = 36*IT_0003;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = m_b*N_d2*e_em*IT_0026*conj(U_sb_10);
    const ccomplex_t IT_0088 = IT_0025*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0091 = IT_0024*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0094 = IT_0017*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + -IT_0092 +
       0.333333333333333*IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = (-3)*IT_0097;
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = sin(beta);
    const ccomplex_t IT_0101 = cpow(IT_0100, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0024*IT_0026*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0099*IT_0103;
    const ccomplex_t IT_0105 = IT_0066*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d1*e_em*IT_0026*conj(U_sb_10);
    const ccomplex_t IT_0108 = IT_0025*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0111 = IT_0024*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0114 = IT_0017*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + -IT_0112 +
       0.333333333333333*IT_0115);
    const ccomplex_t IT_0117 = 3*IT_0116;
    const ccomplex_t IT_0118 = 0.166666666666667*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0024*IT_0026*IT_0101;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = IT_0015*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0024*IT_0026*IT_0101;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = m_b*N_d3*e_em*IT_0026*conj(U_sb_10);
    const ccomplex_t IT_0127 = IT_0025*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0130 = IT_0017*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0133 = IT_0024*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0128 + 0.333333333333333
      *IT_0131 + -IT_0134);
    const ccomplex_t IT_0136 = 3*IT_0135;
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0138 = IT_0125*IT_0137;
    const ccomplex_t IT_0139 = IT_0037*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0024*IT_0026*IT_0101;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0144 = IT_0017*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0147 = IT_0024*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = m_b*N_d4*e_em*IT_0026*conj(U_sb_10);
    const ccomplex_t IT_0150 = IT_0025*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + (-3)*IT_0148 + 3
      *IT_0151);
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0152;
    const ccomplex_t IT_0154 = IT_0142*IT_0153;
    const ccomplex_t IT_0155 = IT_0052*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0106 + -IT_0123 + -IT_0140 + -IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0159 = IT_0039*IT_0137;
    const ccomplex_t IT_0160 = IT_0037*IT_0158*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0162 = IT_0009*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0164 = IT_0008*IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0166 = IT_0033*IT_0118;
    const ccomplex_t IT_0167 = IT_0015*IT_0165*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0169 = IT_0076*IT_0099;
    const ccomplex_t IT_0170 = IT_0066*IT_0168*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0172 = IT_0054*IT_0153;
    const ccomplex_t IT_0173 = IT_0052*IT_0171*IT_0172;
    const ccomplex_t IT_0174 = m_u*IT_0123;
    const ccomplex_t IT_0175 = m_u*IT_0106;
    const ccomplex_t IT_0176 = m_u*IT_0140;
    const ccomplex_t IT_0177 = m_u*IT_0156;
    const ccomplex_t IT_0178 = -IT_0160 + 0.166666666666667*IT_0164 + -IT_0167
       + -IT_0170 + -IT_0173 + IT_0174 + IT_0175 + IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0180 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0024*IT_0026*IT_0179*V_ub_mod;
    const ccomplex_t IT_0181 = 0.5*IT_0180;
    const ccomplex_t IT_0182 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0183 = (0 + _Complex_I*1.4142135623731)*m_u*e_em*mu_h
      *IT_0021*IT_0024*IT_0026*IT_0101*IT_0182*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0181*IT_0184;
    const ccomplex_t IT_0186 = pow(m_W, 2);
    const ccomplex_t IT_0187 = pow(m_sb_1, 2);
    const ccomplex_t IT_0188 = cpow((-2)*s_12 + -IT_0000 + IT_0186 + -IT_0187 
      + -reg_prop, -1);
    const ccomplex_t IT_0189 = IT_0185*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0024*IT_0026*IT_0182*(mu_h*conj(U_sb_00) + m_b*IT_0022*IT_0101*conj
      (U_sb_10))*V_ub_mod;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0021*IT_0024*IT_0026*IT_0101*IT_0179*V_ub_mod;
    const ccomplex_t IT_0194 = 0.5*IT_0193;
    const ccomplex_t IT_0195 = IT_0192*IT_0194;
    const ccomplex_t IT_0196 = cpow((-2)*s_12 + -IT_0000 + -IT_0187 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0197 = IT_0195*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0164 + -IT_0190 + -IT_0198;
    const ccomplex_t IT_0200 = IT_0048*IT_0125;
    const ccomplex_t IT_0201 = IT_0037*IT_0158*IT_0200;
    const ccomplex_t IT_0202 = IT_0031*IT_0120;
    const ccomplex_t IT_0203 = IT_0015*IT_0165*IT_0202;
    const ccomplex_t IT_0204 = IT_0074*IT_0103;
    const ccomplex_t IT_0205 = IT_0066*IT_0168*IT_0204;
    const ccomplex_t IT_0206 = IT_0062*IT_0142;
    const ccomplex_t IT_0207 = IT_0052*IT_0171*IT_0206;
    const ccomplex_t IT_0208 = m_u*IT_0013;
    const ccomplex_t IT_0209 = m_u*IT_0036;
    const ccomplex_t IT_0210 = m_u*IT_0051;
    const ccomplex_t IT_0211 = m_u*IT_0065;
    const ccomplex_t IT_0212 = m_u*IT_0079;
    const ccomplex_t IT_0213 = -IT_0201 + -IT_0203 + -IT_0205 + -IT_0207 + (
      -0.166666666666667)*IT_0208 + IT_0209 + IT_0210 + IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0024*IT_0026*IT_0179*V_ub_mod;
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = IT_0184*IT_0215;
    const ccomplex_t IT_0217 = IT_0188*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0022*IT_0024*IT_0026*IT_0100*IT_0179*V_ub_mod;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = IT_0192*IT_0220;
    const ccomplex_t IT_0222 = IT_0196*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = 0.5*IT_0208 + -IT_0218 + -IT_0223;
    const ccomplex_t IT_0225 = 6*s_34;
    const ccomplex_t IT_0226 = 18*s_34;
    const ccomplex_t IT_0227 = m_u*s_14;
    const ccomplex_t IT_0228 = 6*IT_0224;
    const ccomplex_t IT_0229 = 18*IT_0213;
    const ccomplex_t IT_0230 = 6*conj(IT_0224);
    const ccomplex_t IT_0231 = 18*conj(IT_0213);
    const ccomplex_t IT_0232 = conj(IT_0157)*(18*IT_0178 + 6*IT_0199) +
       IT_0157*(18*conj(IT_0178) + 6*conj(IT_0199)) + 6*conj(IT_0014)*(IT_0213 +
       3*IT_0224) + 6*IT_0014*(conj(IT_0213) + 3*conj(IT_0224)) + IT_0083*
      (IT_0228 + IT_0229) + IT_0082*(IT_0230 + IT_0231);
    const ccomplex_t IT_0233 = m_b*m_u*IT_0000;
    const ccomplex_t IT_0234 = (-18)*IT_0157;
    const ccomplex_t IT_0235 = (-18)*conj(IT_0157);
    const ccomplex_t IT_0236 = m_b*m_u;
    const ccomplex_t IT_0237 = (-18)*IT_0178;
    const ccomplex_t IT_0238 = (-6)*IT_0199;
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = (-18)*conj(IT_0178);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0199);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = (-18)*conj(IT_0199);
    const ccomplex_t IT_0244 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = (-18)*IT_0199;
    const ccomplex_t IT_0247 = (-0.166666666666667)*IT_0246;
    const ccomplex_t IT_0248 = m_b*s_13;
    const ccomplex_t IT_0249 = (conj(IT_0157)*((-18)*IT_0213 + (-6)*IT_0224) +
       IT_0157*((-18)*conj(IT_0213) + (-6)*conj(IT_0224)) + IT_0083*IT_0239 +
       IT_0082*IT_0242 + IT_0014*IT_0245 + (-6)*conj(IT_0014)*(IT_0178 + IT_0247
      ))*IT_0248;
    const ccomplex_t IT_0250 = IT_0005*(conj(IT_0014)*IT_0082 + IT_0014
      *IT_0083) + IT_0086*(IT_0014*conj(IT_0014) + IT_0082*IT_0083 + IT_0157
      *conj(IT_0157)) + (conj(IT_0178)*IT_0199 + IT_0178*conj(IT_0199) + conj
      (IT_0213)*IT_0224 + IT_0213*conj(IT_0224))*IT_0225 + (IT_0178*conj(IT_0178
      ) + IT_0199*conj(IT_0199) + IT_0213*conj(IT_0213) + IT_0224*conj(IT_0224))
      *IT_0226 + IT_0227*IT_0232 + (-6)*IT_0233*(conj(IT_0014)*IT_0157 + IT_0014
      *conj(IT_0157) + (-0.166666666666667)*IT_0083*IT_0234 + (
      -0.166666666666667)*IT_0082*IT_0235) + IT_0236*(conj(IT_0213)*IT_0239 +
       IT_0213*IT_0242 + IT_0224*IT_0245 + (-6)*conj(IT_0224)*(IT_0178 + IT_0247
      )) + IT_0249;
    return create_ccomplex_return(IT_0250);
}

