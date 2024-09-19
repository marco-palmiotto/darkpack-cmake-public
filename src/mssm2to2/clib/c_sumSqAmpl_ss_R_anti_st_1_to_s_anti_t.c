#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_st_1_to_s_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_st_1_to_s_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = (-0.333333333333333)*IT_0002;
    const ccomplex_t IT_0004 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = m_t*N_u2*e_em*IT_0015*conj(U_st_10);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0006 + 3*IT_0011 + 3
      *IT_0018);
    const ccomplex_t IT_0020 = 0.166666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0023 = pow(m_s, 2);
    const ccomplex_t IT_0024 = pow(m_ss_R, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0021*IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0028 = (-0.333333333333333)*IT_0027;
    const ccomplex_t IT_0029 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0033 = IT_0008*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = m_t*N_u3*e_em*IT_0015*conj(U_st_10);
    const ccomplex_t IT_0036 = IT_0014*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0031 + 3*IT_0034 + 3
      *IT_0037);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0040*IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = m_t*N_u4*e_em*IT_0015*conj(U_st_10);
    const ccomplex_t IT_0053 = IT_0014*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0048 + 3*IT_0051 + 3
      *IT_0054);
    const ccomplex_t IT_0056 = 0.166666666666667*IT_0055;
    const ccomplex_t IT_0057 = IT_0045*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = pow(m_sG, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -IT_0064 +
       reg_prop, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0067 = IT_0063*IT_0065*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0069 = cos(beta);
    const ccomplex_t IT_0070 = cpow(IT_0069, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0008*IT_0015*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0074 = IT_0001*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0077 = IT_0008*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = m_t*N_u1*e_em*IT_0015*conj(U_st_10);
    const ccomplex_t IT_0080 = IT_0014*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0075 + 3*IT_0078 + 3
      *IT_0081);
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0072*IT_0083;
    const ccomplex_t IT_0085 = IT_0068*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = m_s*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0008*IT_0015*IT_0070;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0020*IT_0089;
    const ccomplex_t IT_0091 = IT_0025*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = m_s*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0008*IT_0015*IT_0070;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0056*IT_0095;
    const ccomplex_t IT_0097 = IT_0059*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = m_s*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0102 = (-0.333333333333333)*IT_0101;
    const ccomplex_t IT_0103 = IT_0083*IT_0102;
    const ccomplex_t IT_0104 = IT_0068*IT_0100*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0008*IT_0015*IT_0070;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0039*IT_0106;
    const ccomplex_t IT_0108 = IT_0042*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = m_s*IT_0109;
    const ccomplex_t IT_0111 = -IT_0026 + -IT_0043 + -IT_0060 +
       0.166666666666667*IT_0067 + IT_0087 + IT_0093 + IT_0099 + -IT_0104 +
       IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0012*conj
      (U_st_00) + 0.707106781186548*m_s*e_em*IT_0008*IT_0015*(mu_h*conj(U_st_00)
       + m_t*IT_0013*IT_0070*conj(U_st_10)));
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0008*IT_0012*IT_0015*IT_0070;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = pow(m_st_1, 2);
    const ccomplex_t IT_0117 = cpow((-2)*s_12 + -IT_0024 + -IT_0116 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0008*IT_0015;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0069 + (
      -0.707106781186548)*m_s*e_em*mu_h*IT_0008*IT_0012*IT_0015*IT_0070)*conj
      (U_st_00);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = pow(m_W, 2);
    const ccomplex_t IT_0126 = cpow((-2)*s_12 + -IT_0024 + -IT_0116 + IT_0125 
      + -reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = (-0.5)*IT_0067 + -IT_0119 + -IT_0128;
    const ccomplex_t IT_0130 = -IT_0086 + -IT_0092 + -IT_0098 + -IT_0109;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0008*IT_0015;
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = IT_0123*IT_0132;
    const ccomplex_t IT_0134 = IT_0126*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0137 = IT_0061*IT_0136;
    const ccomplex_t IT_0138 = IT_0065*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = m_s*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0008*IT_0013*IT_0015*IT_0069;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = IT_0112*IT_0142;
    const ccomplex_t IT_0144 = IT_0117*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = -IT_0135 + 0.5*IT_0140 + -IT_0145;
    const ccomplex_t IT_0147 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0148 = IT_0001*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = m_t*conj(N_u1)*e_em*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0151 = IT_0014*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + (-0.75)*IT_0152);
    const ccomplex_t IT_0154 = 0.666666666666667*IT_0153;
    const ccomplex_t IT_0155 = IT_0102*IT_0154;
    const ccomplex_t IT_0156 = IT_0068*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0159 = IT_0001*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = m_t*conj(N_u3)*e_em*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0162 = IT_0014*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + (-0.75)*IT_0163);
    const ccomplex_t IT_0165 = 0.666666666666667*IT_0164;
    const ccomplex_t IT_0166 = IT_0028*IT_0165;
    const ccomplex_t IT_0167 = IT_0042*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0170 = IT_0001*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = m_t*conj(N_u2)*e_em*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0173 = IT_0014*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + (-0.75)*IT_0174);
    const ccomplex_t IT_0176 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0177 = IT_0003*IT_0176;
    const ccomplex_t IT_0178 = IT_0025*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = m_t*conj(N_u4)*e_em*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0184 = IT_0014*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + (-0.75)*IT_0185);
    const ccomplex_t IT_0187 = 0.666666666666667*IT_0186;
    const ccomplex_t IT_0188 = IT_0045*IT_0187;
    const ccomplex_t IT_0189 = IT_0059*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = 0.166666666666667*IT_0139 + -IT_0157 + -IT_0168
       + -IT_0179 + -IT_0190;
    const ccomplex_t IT_0192 = (-0.5)*IT_0139;
    const ccomplex_t IT_0193 = IT_0072*IT_0154;
    const ccomplex_t IT_0194 = IT_0068*IT_0100*IT_0193;
    const ccomplex_t IT_0195 = IT_0089*IT_0176;
    const ccomplex_t IT_0196 = IT_0022*IT_0025*IT_0195;
    const ccomplex_t IT_0197 = IT_0106*IT_0165;
    const ccomplex_t IT_0198 = IT_0041*IT_0042*IT_0197;
    const ccomplex_t IT_0199 = IT_0095*IT_0187;
    const ccomplex_t IT_0200 = IT_0058*IT_0059*IT_0199;
    const ccomplex_t IT_0201 = m_s*IT_0179;
    const ccomplex_t IT_0202 = m_s*IT_0168;
    const ccomplex_t IT_0203 = m_s*IT_0190;
    const ccomplex_t IT_0204 = m_s*IT_0157;
    const ccomplex_t IT_0205 = (-0.166666666666667)*IT_0140 + -IT_0194 + 
      -IT_0196 + -IT_0198 + -IT_0200 + IT_0201 + IT_0202 + IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = (18*conj(IT_0111) + 6*conj(IT_0129))*IT_0130 + 
      (18*IT_0111 + 6*IT_0129)*conj(IT_0130) + 6*conj(IT_0146)*(IT_0191 + 3
      *IT_0192) + 6*IT_0146*(conj(IT_0191) + 3*conj(IT_0192)) + (18*conj(IT_0191
      ) + 6*conj(IT_0192))*IT_0205 + (18*IT_0191 + 6*IT_0192)*conj(IT_0205);
    const ccomplex_t IT_0207 = m_s*s_14;
    const ccomplex_t IT_0208 = s_34*IT_0024;
    const ccomplex_t IT_0209 = (-18)*IT_0208;
    const ccomplex_t IT_0210 = s_13*s_14;
    const ccomplex_t IT_0211 = 36*IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = m_s*m_t*IT_0024;
    const ccomplex_t IT_0214 = (-6)*conj(IT_0192);
    const ccomplex_t IT_0215 = (-18)*conj(IT_0191);
    const ccomplex_t IT_0216 = IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = 6*s_34;
    const ccomplex_t IT_0218 = 18*s_34;
    const ccomplex_t IT_0219 = m_s*m_t;
    const ccomplex_t IT_0220 = (-6)*conj(IT_0146) + (-18)*conj(IT_0205);
    const ccomplex_t IT_0221 = (-18)*conj(IT_0111);
    const ccomplex_t IT_0222 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = (-18)*conj(IT_0129);
    const ccomplex_t IT_0225 = (-6)*conj(IT_0111);
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = (-18)*IT_0129;
    const ccomplex_t IT_0228 = (-6)*IT_0129;
    const ccomplex_t IT_0229 = m_t*s_13;
    const ccomplex_t IT_0230 = 12*IT_0210;
    const ccomplex_t IT_0231 = (-6)*IT_0208;
    const ccomplex_t IT_0232 = IT_0206*IT_0207 + (IT_0130*conj(IT_0130) +
       IT_0191*conj(IT_0191) + IT_0192*conj(IT_0192))*IT_0212 + IT_0213*(conj
      (IT_0130)*((-18)*IT_0191 + (-6)*IT_0192) + IT_0130*IT_0216) + (conj
      (IT_0111)*IT_0129 + IT_0111*conj(IT_0129) + conj(IT_0146)*IT_0205 +
       IT_0146*conj(IT_0205))*IT_0217 + (IT_0111*conj(IT_0111) + IT_0129*conj
      (IT_0129) + IT_0146*conj(IT_0146) + IT_0205*conj(IT_0205))*IT_0218 +
       IT_0219*(IT_0111*IT_0220 + IT_0205*IT_0223 + IT_0146*IT_0226 + conj
      (IT_0146)*IT_0227 + conj(IT_0205)*IT_0228) + (-6)*(conj(IT_0130)*(IT_0146 
      + 3*IT_0205) + -1./6*IT_0111*IT_0216 + -1./6*IT_0130*IT_0220 + -1./6
      *IT_0191*IT_0223 + -1./6*IT_0192*IT_0226 + -1./6*conj(IT_0192)*IT_0227 + 
      -1./6*conj(IT_0191)*IT_0228)*IT_0229 + (conj(IT_0191)*IT_0192 + IT_0191
      *conj(IT_0192))*(IT_0230 + IT_0231);
    return create_ccomplex_return(IT_0232);
}

