#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_st_1_to_d_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_st_1_to_d_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_R = param->m_sd_R;
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
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = (-0.333333333333333)*IT_0002;
    const ccomplex_t IT_0004 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_t*conj(N_u3)*e_em*IT_0012*conj(U_st_00);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0006 + (-0.75)*IT_0015);
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = pow(m_d, 2);
    const ccomplex_t IT_0020 = pow(m_sd_R, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0025 = (-0.333333333333333)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = m_t*conj(N_u4)*e_em*IT_0012*conj(U_st_00);
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-0.75)*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0025*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0039 = (-0.333333333333333)*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_t*conj(N_u1)*e_em*IT_0012*conj(U_st_00);
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-0.75)*IT_0045);
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0039*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = pow(m_sG, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -IT_0052 +
       reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0062 = IT_0001*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = m_t*conj(N_u2)*e_em*IT_0012*conj(U_st_00);
    const ccomplex_t IT_0065 = IT_0011*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + (-0.75)*IT_0066);
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0067;
    const ccomplex_t IT_0069 = IT_0060*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0023 + -IT_0037 + -IT_0051 +
       0.166666666666667*IT_0058 + -IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0058;
    const ccomplex_t IT_0075 = cos(beta);
    const ccomplex_t IT_0076 = cpow(IT_0075, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0010*IT_0012*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = m_t*N_u4*e_em*IT_0012*conj(U_st_10);
    const ccomplex_t IT_0086 = IT_0011*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0081 + 3*IT_0084 + 3
      *IT_0087);
    const ccomplex_t IT_0089 = 0.166666666666667*IT_0088;
    const ccomplex_t IT_0090 = IT_0078*IT_0089;
    const ccomplex_t IT_0091 = IT_0035*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0010*IT_0012*IT_0076;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0096 = IT_0001*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0099 = IT_0010*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = m_t*N_u1*e_em*IT_0012*conj(U_st_10);
    const ccomplex_t IT_0102 = IT_0011*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0097 + 3*IT_0100 + 3
      *IT_0103);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0094*IT_0105;
    const ccomplex_t IT_0107 = IT_0049*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0010*IT_0012*IT_0076;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0112 = IT_0001*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0115 = IT_0010*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = m_t*N_u2*e_em*IT_0012*conj(U_st_10);
    const ccomplex_t IT_0118 = IT_0011*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0113 + 3*IT_0116 + 3
      *IT_0119);
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0122 = IT_0110*IT_0121;
    const ccomplex_t IT_0123 = IT_0070*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0129 = IT_0010*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = m_t*N_u3*e_em*IT_0012*conj(U_st_10);
    const ccomplex_t IT_0132 = IT_0011*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0127 + 3*IT_0130 + 3
      *IT_0133);
    const ccomplex_t IT_0135 = 0.166666666666667*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0010*IT_0012*IT_0076;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0135*IT_0137;
    const ccomplex_t IT_0139 = IT_0021*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = -IT_0092 + -IT_0108 + -IT_0124 + -IT_0140;
    const ccomplex_t IT_0142 = m_d*m_t*IT_0020;
    const ccomplex_t IT_0143 = s_34*IT_0020;
    const ccomplex_t IT_0144 = (-18)*IT_0143;
    const ccomplex_t IT_0145 = s_13*s_14;
    const ccomplex_t IT_0146 = 36*IT_0145;
    const ccomplex_t IT_0147 = IT_0144 + IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_td)
      *e_em*IT_0010*IT_0012;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*V_td*(A_b*IT_0075 + (
      -0.707106781186548)*m_d*e_em*mu_h*IT_0007*IT_0010*IT_0012*IT_0076)*conj
      (U_st_00);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = pow(m_W, 2);
    const ccomplex_t IT_0154 = pow(m_st_1, 2);
    const ccomplex_t IT_0155 = cpow((-2)*s_12 + -IT_0020 + IT_0153 + -IT_0154 
      + -reg_prop, -1);
    const ccomplex_t IT_0156 = IT_0152*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*V_td*(A_b*IT_0007*conj
      (U_st_00) + 0.707106781186548*m_d*e_em*IT_0010*IT_0012*(mu_h*conj(U_st_00)
       + m_t*IT_0008*IT_0076*conj(U_st_10)));
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_td)
      *e_em*IT_0007*IT_0010*IT_0012*IT_0076;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0158*IT_0160;
    const ccomplex_t IT_0162 = cpow((-2)*s_12 + -IT_0020 + -IT_0154 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0166 = IT_0054*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0168 = IT_0053*IT_0166*IT_0167;
    const ccomplex_t IT_0169 = -IT_0157 + -IT_0164 + (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0060*IT_0121;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0172 = IT_0070*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = IT_0003*IT_0135;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0175 = IT_0021*IT_0173*IT_0174;
    const ccomplex_t IT_0176 = IT_0025*IT_0089;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0178 = IT_0035*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = m_d*IT_0108;
    const ccomplex_t IT_0180 = m_d*IT_0124;
    const ccomplex_t IT_0181 = m_d*IT_0140;
    const ccomplex_t IT_0182 = m_d*IT_0092;
    const ccomplex_t IT_0183 = IT_0039*IT_0105;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0185 = IT_0049*IT_0183*IT_0184;
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0168 + -IT_0172 + -IT_0175
       + -IT_0178 + IT_0179 + IT_0180 + IT_0181 + IT_0182 + -IT_0185;
    const ccomplex_t IT_0187 = IT_0047*IT_0094;
    const ccomplex_t IT_0188 = IT_0049*IT_0184*IT_0187;
    const ccomplex_t IT_0189 = IT_0068*IT_0110;
    const ccomplex_t IT_0190 = IT_0070*IT_0171*IT_0189;
    const ccomplex_t IT_0191 = IT_0017*IT_0137;
    const ccomplex_t IT_0192 = IT_0021*IT_0174*IT_0191;
    const ccomplex_t IT_0193 = IT_0033*IT_0078;
    const ccomplex_t IT_0194 = IT_0035*IT_0177*IT_0193;
    const ccomplex_t IT_0195 = m_d*IT_0072;
    const ccomplex_t IT_0196 = m_d*IT_0023;
    const ccomplex_t IT_0197 = m_d*IT_0037;
    const ccomplex_t IT_0198 = m_d*IT_0051;
    const ccomplex_t IT_0199 = m_d*IT_0058;
    const ccomplex_t IT_0200 = -IT_0188 + -IT_0190 + -IT_0192 + -IT_0194 +
       IT_0195 + IT_0196 + IT_0197 + IT_0198 + (-0.166666666666667)*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *e_em*IT_0010*IT_0012;
    const ccomplex_t IT_0202 = 0.5*IT_0201;
    const ccomplex_t IT_0203 = IT_0151*IT_0202;
    const ccomplex_t IT_0204 = IT_0155*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *e_em*IT_0008*IT_0010*IT_0012*IT_0075;
    const ccomplex_t IT_0207 = 0.5*IT_0206;
    const ccomplex_t IT_0208 = IT_0158*IT_0207;
    const ccomplex_t IT_0209 = IT_0162*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0199 + -IT_0205 + -IT_0210;
    const ccomplex_t IT_0212 = 6*s_34;
    const ccomplex_t IT_0213 = 18*s_34;
    const ccomplex_t IT_0214 = m_t*s_13;
    const ccomplex_t IT_0215 = (-18)*conj(IT_0169);
    const ccomplex_t IT_0216 = (-6)*conj(IT_0186);
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = (-18)*conj(IT_0186);
    const ccomplex_t IT_0219 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0220 = IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = (-18)*IT_0186;
    const ccomplex_t IT_0222 = (-6)*IT_0186;
    const ccomplex_t IT_0223 = (conj(IT_0073) + 3*conj(IT_0074))*IT_0169 + 3
      *conj(IT_0141)*(IT_0200 + 0.333333333333333*IT_0211) + 3*IT_0141*(conj
      (IT_0200) + 0.333333333333333*conj(IT_0211)) + (-0.166666666666667)
      *IT_0074*IT_0217 + (-0.166666666666667)*IT_0073*IT_0220 + (
      -0.166666666666667)*conj(IT_0073)*IT_0221 + (-0.166666666666667)*conj
      (IT_0074)*IT_0222;
    const ccomplex_t IT_0224 = 12*IT_0145;
    const ccomplex_t IT_0225 = (-6)*IT_0143;
    const ccomplex_t IT_0226 = m_d*m_t;
    const ccomplex_t IT_0227 = 3*conj(IT_0211);
    const ccomplex_t IT_0228 = m_d*s_14;
    const ccomplex_t IT_0229 = conj(IT_0141)*(IT_0169 + 3*IT_0186) + IT_0141*
      (conj(IT_0169) + 3*conj(IT_0186)) + conj(IT_0074)*(IT_0200 + 3*IT_0211) +
       0.166666666666667*conj(IT_0073)*(18*IT_0200 + 6*IT_0211) +
       0.166666666666667*IT_0073*(18*conj(IT_0200) + 6*conj(IT_0211)) + IT_0074*
      (conj(IT_0200) + IT_0227);
    const ccomplex_t IT_0230 = (((-18)*conj(IT_0073) + (-6)*conj(IT_0074))
      *IT_0141 + ((-18)*IT_0073 + (-6)*IT_0074)*conj(IT_0141))*IT_0142 + 
      (IT_0073*conj(IT_0073) + IT_0074*conj(IT_0074) + IT_0141*conj(IT_0141))
      *IT_0147 + (conj(IT_0169)*IT_0186 + IT_0169*conj(IT_0186) + conj(IT_0200)
      *IT_0211 + IT_0200*conj(IT_0211))*IT_0212 + (IT_0169*conj(IT_0169) +
       IT_0186*conj(IT_0186) + IT_0200*conj(IT_0200) + IT_0211*conj(IT_0211))
      *IT_0213 + (-6)*IT_0214*IT_0223 + (conj(IT_0073)*IT_0074 + IT_0073*conj
      (IT_0074))*(IT_0224 + IT_0225) + IT_0226*(IT_0211*IT_0217 + IT_0200
      *IT_0220 + conj(IT_0200)*IT_0221 + conj(IT_0211)*IT_0222 + (-6)*IT_0169*
      (conj(IT_0200) + IT_0227)) + 6*IT_0228*IT_0229;
    return create_ccomplex_return(IT_0230);
}

