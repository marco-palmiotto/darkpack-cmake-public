#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_sb_1_to_c_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_sb_1_to_c_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t V_cb = param->V_cb;
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
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sc_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = m_b*N_d1*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0012 + -IT_0018 + (
      -0.333333333333333)*IT_0023);
    const ccomplex_t IT_0025 = (-3)*IT_0024;
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = IT_0002*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0034 = IT_0020*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0037 = IT_0007*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = m_b*N_d3*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0040 = IT_0015*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0035 + (-3)*IT_0038 + 3
      *IT_0041);
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = IT_0032*IT_0043;
    const ccomplex_t IT_0045 = IT_0030*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = m_b*N_d2*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0054 = IT_0015*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0052 + -IT_0055 + (
      -0.333333333333333)*IT_0058);
    const ccomplex_t IT_0060 = (-3)*IT_0059;
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0049*IT_0061;
    const ccomplex_t IT_0063 = IT_0047*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0069 = IT_0020*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0072 = IT_0007*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = m_b*N_d4*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0075 = IT_0015*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0070 + (-3)*IT_0073 + 3
      *IT_0076);
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0067*IT_0078;
    const ccomplex_t IT_0080 = IT_0065*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0029 + -IT_0046 + -IT_0064 + -IT_0081;
    const ccomplex_t IT_0083 = pow(m_sG, 2);
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0083 +
       reg_prop, -1);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = IT_0084*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0020;
    const ccomplex_t IT_0092 = 0.666666666666667*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0094 = IT_0020*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = m_b*conj(N_d3)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0097 = IT_0015*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 1.5*IT_0098);
    const ccomplex_t IT_0100 = (-0.333333333333333)*IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = IT_0030*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0020;
    const ccomplex_t IT_0105 = 0.666666666666667*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0107 = IT_0020*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = m_b*conj(N_d2)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0110 = IT_0015*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + 1.5*IT_0111);
    const ccomplex_t IT_0113 = (-0.333333333333333)*IT_0112;
    const ccomplex_t IT_0114 = IT_0105*IT_0113;
    const ccomplex_t IT_0115 = IT_0047*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0020;
    const ccomplex_t IT_0118 = 0.666666666666667*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0120 = IT_0020*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_b*conj(N_d4)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0123 = IT_0015*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 1.5*IT_0124);
    const ccomplex_t IT_0126 = (-0.333333333333333)*IT_0125;
    const ccomplex_t IT_0127 = IT_0118*IT_0126;
    const ccomplex_t IT_0128 = IT_0065*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0020;
    const ccomplex_t IT_0131 = 0.666666666666667*IT_0130;
    const ccomplex_t IT_0132 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0133 = IT_0020*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(N_d1)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0136 = IT_0015*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + 1.5*IT_0137);
    const ccomplex_t IT_0139 = (-0.333333333333333)*IT_0138;
    const ccomplex_t IT_0140 = IT_0131*IT_0139;
    const ccomplex_t IT_0141 = IT_0002*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = 0.166666666666667*IT_0089 + -IT_0103 + -IT_0116
       + -IT_0129 + -IT_0142;
    const ccomplex_t IT_0144 = s_34*IT_0001;
    const ccomplex_t IT_0145 = (-18)*IT_0144;
    const ccomplex_t IT_0146 = s_13*s_14;
    const ccomplex_t IT_0147 = 36*IT_0146;
    const ccomplex_t IT_0148 = IT_0145 + IT_0147;
    const ccomplex_t IT_0149 = pow(m_sb_1, 2);
    const ccomplex_t IT_0150 = cpow((-2)*s_12 + -IT_0001 + -IT_0149 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0003*IT_0005*IT_0007*IT_0013;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0003*IT_0007*(mu_h*conj(U_sb_00) + m_b*IT_0005*IT_0014*conj(U_sb_10));
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = IT_0150*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0160 = IT_0085*IT_0159;
    const ccomplex_t IT_0161 = IT_0084*IT_0158*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0003*IT_0007;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *mu_h*IT_0003*IT_0005*IT_0007*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = pow(m_W, 2);
    const ccomplex_t IT_0168 = cpow((-2)*s_12 + -IT_0001 + -IT_0149 + IT_0167 
      + -reg_prop, -1);
    const ccomplex_t IT_0169 = IT_0166*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = -IT_0157 + (-0.5)*IT_0161 + -IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0173 = IT_0043*IT_0092;
    const ccomplex_t IT_0174 = IT_0030*IT_0172*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0176 = IT_0026*IT_0131;
    const ccomplex_t IT_0177 = IT_0002*IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0179 = IT_0061*IT_0105;
    const ccomplex_t IT_0180 = IT_0047*IT_0178*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0182 = IT_0078*IT_0118;
    const ccomplex_t IT_0183 = IT_0065*IT_0181*IT_0182;
    const ccomplex_t IT_0184 = m_c*IT_0029;
    const ccomplex_t IT_0185 = m_c*IT_0064;
    const ccomplex_t IT_0186 = m_c*IT_0046;
    const ccomplex_t IT_0187 = m_c*IT_0081;
    const ccomplex_t IT_0188 = 0.166666666666667*IT_0161 + -IT_0174 + -IT_0177
       + -IT_0180 + -IT_0183 + IT_0184 + IT_0185 + IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0003*IT_0007;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0165*IT_0190;
    const ccomplex_t IT_0192 = IT_0168*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0003*IT_0004*IT_0007*IT_0014;
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = IT_0154*IT_0195;
    const ccomplex_t IT_0197 = IT_0150*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = m_c*IT_0089;
    const ccomplex_t IT_0200 = -IT_0193 + -IT_0198 + 0.5*IT_0199;
    const ccomplex_t IT_0201 = IT_0009*IT_0139;
    const ccomplex_t IT_0202 = IT_0002*IT_0175*IT_0201;
    const ccomplex_t IT_0203 = m_c*IT_0142;
    const ccomplex_t IT_0204 = IT_0032*IT_0100;
    const ccomplex_t IT_0205 = IT_0030*IT_0172*IT_0204;
    const ccomplex_t IT_0206 = IT_0049*IT_0113;
    const ccomplex_t IT_0207 = IT_0047*IT_0178*IT_0206;
    const ccomplex_t IT_0208 = m_c*IT_0103;
    const ccomplex_t IT_0209 = IT_0067*IT_0126;
    const ccomplex_t IT_0210 = IT_0065*IT_0181*IT_0209;
    const ccomplex_t IT_0211 = m_c*IT_0116;
    const ccomplex_t IT_0212 = m_c*IT_0129;
    const ccomplex_t IT_0213 = (-0.166666666666667)*IT_0199 + -IT_0202 +
       IT_0203 + -IT_0205 + -IT_0207 + IT_0208 + -IT_0210 + IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = 6*s_34;
    const ccomplex_t IT_0215 = 18*s_34;
    const ccomplex_t IT_0216 = m_b*m_c;
    const ccomplex_t IT_0217 = (-18)*conj(IT_0171);
    const ccomplex_t IT_0218 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = (-18)*conj(IT_0188);
    const ccomplex_t IT_0221 = (-6)*conj(IT_0171);
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = (-18)*IT_0188;
    const ccomplex_t IT_0224 = (-6)*IT_0188;
    const ccomplex_t IT_0225 = m_b*s_13;
    const ccomplex_t IT_0226 = 0.333333333333333*conj(IT_0143);
    const ccomplex_t IT_0227 = conj(IT_0090) + IT_0226;
    const ccomplex_t IT_0228 = conj(IT_0082)*(IT_0200 + 3*IT_0213) + IT_0082*
      (conj(IT_0200) + 3*conj(IT_0213)) + (-0.166666666666667)*IT_0090*IT_0219 +
       (-0.166666666666667)*IT_0143*IT_0222 + (-0.166666666666667)*conj(IT_0143)
      *IT_0223 + (-0.166666666666667)*conj(IT_0090)*IT_0224 + 3*IT_0171*IT_0227;
    const ccomplex_t IT_0229 = 12*IT_0146;
    const ccomplex_t IT_0230 = (-6)*IT_0144;
    const ccomplex_t IT_0231 = m_c*s_14;
    const ccomplex_t IT_0232 = 3*IT_0143;
    const ccomplex_t IT_0233 = 3*conj(IT_0143);
    const ccomplex_t IT_0234 = conj(IT_0082)*(IT_0171 + 3*IT_0188) + IT_0082*
      (conj(IT_0171) + 3*conj(IT_0188)) + 3*(IT_0090 + 0.333333333333333*IT_0143
      )*conj(IT_0200) + 3*IT_0200*IT_0227 + conj(IT_0213)*(IT_0090 + IT_0232) +
       IT_0213*(conj(IT_0090) + IT_0233);
    const ccomplex_t IT_0235 = m_b*m_c*IT_0001;
    const ccomplex_t IT_0236 = (IT_0082*conj(IT_0082) + IT_0090*conj(IT_0090) 
      + IT_0143*conj(IT_0143))*IT_0148 + (conj(IT_0171)*IT_0188 + IT_0171*conj
      (IT_0188) + conj(IT_0200)*IT_0213 + IT_0200*conj(IT_0213))*IT_0214 + 
      (IT_0171*conj(IT_0171) + IT_0188*conj(IT_0188) + IT_0200*conj(IT_0200) +
       IT_0213*conj(IT_0213))*IT_0215 + IT_0216*(IT_0171*((-18)*conj(IT_0200) + 
      (-6)*conj(IT_0213)) + IT_0200*IT_0219 + IT_0213*IT_0222 + conj(IT_0213)
      *IT_0223 + conj(IT_0200)*IT_0224) + (-6)*IT_0225*IT_0228 + (conj(IT_0090)
      *IT_0143 + IT_0090*conj(IT_0143))*(IT_0229 + IT_0230) + 6*IT_0231*IT_0234 
      + (-6)*(conj(IT_0082)*(IT_0090 + IT_0232) + IT_0082*(conj(IT_0090) +
       IT_0233))*IT_0235;
    return create_ccomplex_return(IT_0236);
}

