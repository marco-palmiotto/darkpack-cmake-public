#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_st_2_to_u_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_st_2_to_u_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_R = param->m_su_R;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_su_R, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0008 = IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0010 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = N_B1*e_em*conj(U_st_01);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = N_W1*e_em*conj(U_st_01);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = sin(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = m_t*N_u1*e_em*IT_0024*conj(U_st_11);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0015 + 3*IT_0020 + 3
      *IT_0027);
    const ccomplex_t IT_0029 = 0.166666666666667*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0012;
    const ccomplex_t IT_0031 = 0.666666666666667*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = IT_0009*IT_0010*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0012;
    const ccomplex_t IT_0037 = 0.666666666666667*IT_0036;
    const ccomplex_t IT_0038 = N_B2*e_em*conj(U_st_01);
    const ccomplex_t IT_0039 = IT_0012*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = N_W2*e_em*conj(U_st_01);
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = m_t*N_u2*e_em*IT_0024*conj(U_st_11);
    const ccomplex_t IT_0045 = IT_0023*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0040 + 3*IT_0043 + 3
      *IT_0046);
    const ccomplex_t IT_0048 = 0.166666666666667*IT_0047;
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = IT_0034*IT_0035*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0012;
    const ccomplex_t IT_0052 = 0.666666666666667*IT_0051;
    const ccomplex_t IT_0053 = N_B3*e_em*conj(U_st_01);
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_W3*e_em*conj(U_st_01);
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = m_t*N_u3*e_em*IT_0024*conj(U_st_11);
    const ccomplex_t IT_0060 = IT_0023*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0055 + 3*IT_0058 + 3
      *IT_0061);
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0064*IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0012;
    const ccomplex_t IT_0069 = 0.666666666666667*IT_0068;
    const ccomplex_t IT_0070 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0071 = IT_0012*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0074 = IT_0017*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = m_t*N_u4*e_em*IT_0024*conj(U_st_11);
    const ccomplex_t IT_0077 = IT_0023*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0072 + 3*IT_0075 + 3
      *IT_0078);
    const ccomplex_t IT_0080 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0017*IT_0022*IT_0024;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0029*IT_0086;
    const ccomplex_t IT_0088 = IT_0010*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_u*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0017*IT_0022*IT_0024;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0048*IT_0092;
    const ccomplex_t IT_0094 = IT_0035*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = m_u*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0017*IT_0022*IT_0024;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0063*IT_0098;
    const ccomplex_t IT_0100 = IT_0066*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_u*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0017*IT_0022*IT_0024;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0080*IT_0104;
    const ccomplex_t IT_0106 = IT_0083*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = m_u*IT_0107;
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0008 + -IT_0033 + -IT_0050
       + -IT_0067 + -IT_0084 + IT_0090 + IT_0096 + IT_0102 + IT_0108;
    const ccomplex_t IT_0110 = (-0.5)*IT_0008;
    const ccomplex_t IT_0111 = -IT_0089 + -IT_0095 + -IT_0101 + -IT_0107;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0113 = IT_0000*IT_0112;
    const ccomplex_t IT_0114 = IT_0006*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = m_u*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = m_t*conj(N_u1)*e_em*IT_0024*conj(U_st_01);
    const ccomplex_t IT_0119 = IT_0023*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0122 = IT_0012*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + (
      -1.33333333333333)*IT_0123);
    const ccomplex_t IT_0125 = (-0.75)*IT_0124;
    const ccomplex_t IT_0126 = 0.666666666666667*IT_0125;
    const ccomplex_t IT_0127 = IT_0031*IT_0126;
    const ccomplex_t IT_0128 = IT_0010*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em*conj(U_st_11);
    const ccomplex_t IT_0131 = IT_0012*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_t*conj(N_u3)*e_em*IT_0024*conj(U_st_01);
    const ccomplex_t IT_0134 = IT_0023*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-0.75)*IT_0135);
    const ccomplex_t IT_0137 = 0.666666666666667*IT_0136;
    const ccomplex_t IT_0138 = IT_0052*IT_0137;
    const ccomplex_t IT_0139 = IT_0066*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B2)*e_em*conj(U_st_11);
    const ccomplex_t IT_0142 = IT_0012*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = m_t*conj(N_u2)*e_em*IT_0024*conj(U_st_01);
    const ccomplex_t IT_0145 = IT_0023*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + (-0.75)*IT_0146);
    const ccomplex_t IT_0148 = 0.666666666666667*IT_0147;
    const ccomplex_t IT_0149 = IT_0037*IT_0148;
    const ccomplex_t IT_0150 = IT_0035*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0153 = IT_0012*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_t*conj(N_u4)*e_em*IT_0024*conj(U_st_01);
    const ccomplex_t IT_0156 = IT_0023*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + (-0.75)*IT_0157);
    const ccomplex_t IT_0159 = 0.666666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0069*IT_0159;
    const ccomplex_t IT_0161 = IT_0083*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = 0.166666666666667*IT_0115 + -IT_0129 + -IT_0140
       + -IT_0151 + -IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0115;
    const ccomplex_t IT_0165 = IT_0104*IT_0159;
    const ccomplex_t IT_0166 = IT_0082*IT_0083*IT_0165;
    const ccomplex_t IT_0167 = m_u*IT_0129;
    const ccomplex_t IT_0168 = IT_0098*IT_0137;
    const ccomplex_t IT_0169 = IT_0065*IT_0066*IT_0168;
    const ccomplex_t IT_0170 = IT_0086*IT_0126;
    const ccomplex_t IT_0171 = IT_0009*IT_0010*IT_0170;
    const ccomplex_t IT_0172 = IT_0092*IT_0148;
    const ccomplex_t IT_0173 = IT_0034*IT_0035*IT_0172;
    const ccomplex_t IT_0174 = m_u*IT_0151;
    const ccomplex_t IT_0175 = m_u*IT_0140;
    const ccomplex_t IT_0176 = m_u*IT_0162;
    const ccomplex_t IT_0177 = (-0.166666666666667)*IT_0116 + -IT_0166 +
       IT_0167 + -IT_0169 + -IT_0171 + -IT_0173 + IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (18*conj(IT_0109) + 6*conj(IT_0110))*IT_0111 + 
      (18*IT_0109 + 6*IT_0110)*conj(IT_0111) + 6*conj(IT_0117)*(IT_0163 + 3
      *IT_0164) + 6*IT_0117*(conj(IT_0163) + 3*conj(IT_0164)) + (18*conj(IT_0163
      ) + 6*conj(IT_0164))*IT_0177 + (18*IT_0163 + 6*IT_0164)*conj(IT_0177);
    const ccomplex_t IT_0179 = m_u*s_14;
    const ccomplex_t IT_0180 = m_t*m_u*IT_0005;
    const ccomplex_t IT_0181 = s_34*IT_0005;
    const ccomplex_t IT_0182 = (-18)*IT_0181;
    const ccomplex_t IT_0183 = s_13*s_14;
    const ccomplex_t IT_0184 = 36*IT_0183;
    const ccomplex_t IT_0185 = IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = 6*s_34;
    const ccomplex_t IT_0187 = 18*s_34;
    const ccomplex_t IT_0188 = m_t*m_u;
    const ccomplex_t IT_0189 = conj(IT_0117) + 3*conj(IT_0177);
    const ccomplex_t IT_0190 = (-18)*conj(IT_0109);
    const ccomplex_t IT_0191 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0192 = IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = (-18)*conj(IT_0110);
    const ccomplex_t IT_0194 = (-6)*conj(IT_0109);
    const ccomplex_t IT_0195 = IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = (-18)*IT_0110;
    const ccomplex_t IT_0197 = (-6)*IT_0110;
    const ccomplex_t IT_0198 = m_t*s_13;
    const ccomplex_t IT_0199 = 3*IT_0109*(conj(IT_0163) + 0.333333333333333
      *conj(IT_0164)) + conj(IT_0111)*(IT_0117 + 3*IT_0177) + IT_0111*IT_0189 + 
      (-0.166666666666667)*IT_0163*IT_0192 + (-0.166666666666667)*IT_0164
      *IT_0195 + (-0.166666666666667)*conj(IT_0164)*IT_0196 + (
      -0.166666666666667)*conj(IT_0163)*IT_0197;
    const ccomplex_t IT_0200 = 12*IT_0183;
    const ccomplex_t IT_0201 = (-6)*IT_0181;
    const ccomplex_t IT_0202 = IT_0178*IT_0179 + (conj(IT_0111)*((-18)*IT_0163
       + (-6)*IT_0164) + IT_0111*((-18)*conj(IT_0163) + (-6)*conj(IT_0164)))
      *IT_0180 + (IT_0111*conj(IT_0111) + IT_0163*conj(IT_0163) + IT_0164*conj
      (IT_0164))*IT_0185 + (conj(IT_0109)*IT_0110 + IT_0109*conj(IT_0110) + conj
      (IT_0117)*IT_0177 + IT_0117*conj(IT_0177))*IT_0186 + (IT_0109*conj(IT_0109
      ) + IT_0110*conj(IT_0110) + IT_0117*conj(IT_0117) + IT_0177*conj(IT_0177))
      *IT_0187 + (-6)*IT_0188*(IT_0109*IT_0189 + (-0.166666666666667)*IT_0177
      *IT_0192 + (-0.166666666666667)*IT_0117*IT_0195 + (-0.166666666666667)
      *conj(IT_0117)*IT_0196 + (-0.166666666666667)*conj(IT_0177)*IT_0197) + (-6
      )*IT_0198*IT_0199 + (conj(IT_0163)*IT_0164 + IT_0163*conj(IT_0164))*
      (IT_0200 + IT_0201);
    return create_ccomplex_return(IT_0202);
}

