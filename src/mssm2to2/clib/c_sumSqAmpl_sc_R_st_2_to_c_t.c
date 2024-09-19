#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_st_2_to_c_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_st_2_to_c_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_R = param->m_sc_R;
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
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sc_R, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0008 = IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0010 = IT_0000*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = m_t*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0008 + (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0017;
    const ccomplex_t IT_0019 = 0.666666666666667*IT_0018;
    const ccomplex_t IT_0020 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = sin(theta_W);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = m_t*N_u1*e_em*IT_0028*U_st_01;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0022 + (-0.75)*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0015*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0017;
    const ccomplex_t IT_0039 = 0.666666666666667*IT_0038;
    const ccomplex_t IT_0040 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_t*N_u2*e_em*IT_0028*U_st_01;
    const ccomplex_t IT_0044 = IT_0027*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-0.75)*IT_0045);
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0039*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0037*IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0017;
    const ccomplex_t IT_0052 = 0.666666666666667*IT_0051;
    const ccomplex_t IT_0053 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0054 = IT_0017*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*N_u4*e_em*IT_0028*U_st_01;
    const ccomplex_t IT_0057 = IT_0027*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.75)*IT_0058);
    const ccomplex_t IT_0060 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0061 = IT_0052*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0061*IT_0062*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0066 = IT_0017*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0069 = IT_0026*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_t*conj(N_u1)*e_em*IT_0028*U_st_11;
    const ccomplex_t IT_0072 = IT_0027*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + 3*IT_0070 + 3
      *IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = IT_0019*IT_0075;
    const ccomplex_t IT_0077 = IT_0035*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = m_t*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0081 = IT_0017*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0084 = IT_0026*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = m_t*conj(N_u2)*e_em*IT_0028*U_st_11;
    const ccomplex_t IT_0087 = IT_0027*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0082 + 3*IT_0085 + 3
      *IT_0088);
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = IT_0039*IT_0090;
    const ccomplex_t IT_0092 = IT_0049*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = m_t*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0017;
    const ccomplex_t IT_0098 = 0.666666666666667*IT_0097;
    const ccomplex_t IT_0099 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0100 = IT_0017*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = m_t*N_u3*e_em*IT_0028*U_st_01;
    const ccomplex_t IT_0103 = IT_0027*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + (-0.75)*IT_0104);
    const ccomplex_t IT_0106 = 0.666666666666667*IT_0105;
    const ccomplex_t IT_0107 = IT_0098*IT_0106;
    const ccomplex_t IT_0108 = IT_0095*IT_0096*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0110 = IT_0017*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0113 = IT_0026*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = m_t*conj(N_u3)*e_em*IT_0028*U_st_11;
    const ccomplex_t IT_0116 = IT_0027*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0111 + 3*IT_0114 + 3
      *IT_0117);
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0118;
    const ccomplex_t IT_0120 = IT_0098*IT_0119;
    const ccomplex_t IT_0121 = IT_0096*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = m_t*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0125 = IT_0017*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0128 = IT_0026*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_t*conj(N_u4)*e_em*IT_0028*U_st_11;
    const ccomplex_t IT_0131 = IT_0027*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0126 + 3*IT_0129 + 3
      *IT_0132);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = IT_0052*IT_0134;
    const ccomplex_t IT_0136 = IT_0063*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = m_t*IT_0137;
    const ccomplex_t IT_0139 = (-0.166666666666667)*IT_0008 +
       0.166666666666667*IT_0013 + IT_0036 + IT_0050 + IT_0064 + -IT_0079 + 
      -IT_0094 + IT_0108 + -IT_0123 + -IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0024*IT_0026*IT_0028;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0033*IT_0141;
    const ccomplex_t IT_0143 = IT_0035*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0024*IT_0026*IT_0028;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = IT_0047*IT_0146;
    const ccomplex_t IT_0148 = IT_0049*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0024*IT_0026*IT_0028;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0106*IT_0151;
    const ccomplex_t IT_0153 = IT_0096*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0024*IT_0026*IT_0028;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = IT_0060*IT_0156;
    const ccomplex_t IT_0158 = IT_0063*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = -IT_0144 + -IT_0149 + -IT_0154 + -IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0012;
    const ccomplex_t IT_0162 = 0.166666666666667*IT_0012 + -IT_0078 + -IT_0093
       + -IT_0122 + -IT_0137;
    const ccomplex_t IT_0163 = pow(m_st_2, 2);
    const ccomplex_t IT_0164 = s_34*IT_0163;
    const ccomplex_t IT_0165 = (-18)*IT_0164;
    const ccomplex_t IT_0166 = s_23*s_24;
    const ccomplex_t IT_0167 = 36*IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = m_c*m_t*IT_0163;
    const ccomplex_t IT_0170 = 3*IT_0162;
    const ccomplex_t IT_0171 = IT_0161 + IT_0170;
    const ccomplex_t IT_0172 = 3*conj(IT_0162);
    const ccomplex_t IT_0173 = conj(IT_0161) + IT_0172;
    const ccomplex_t IT_0174 = m_c*s_24;
    const ccomplex_t IT_0175 = IT_0075*IT_0141;
    const ccomplex_t IT_0176 = IT_0015*IT_0035*IT_0175;
    const ccomplex_t IT_0177 = IT_0090*IT_0146;
    const ccomplex_t IT_0178 = IT_0037*IT_0049*IT_0177;
    const ccomplex_t IT_0179 = m_t*IT_0144;
    const ccomplex_t IT_0180 = IT_0119*IT_0151;
    const ccomplex_t IT_0181 = IT_0095*IT_0096*IT_0180;
    const ccomplex_t IT_0182 = IT_0134*IT_0156;
    const ccomplex_t IT_0183 = IT_0062*IT_0063*IT_0182;
    const ccomplex_t IT_0184 = m_t*IT_0149;
    const ccomplex_t IT_0185 = m_t*IT_0154;
    const ccomplex_t IT_0186 = m_t*IT_0159;
    const ccomplex_t IT_0187 = IT_0176 + IT_0178 + -IT_0179 + IT_0181 +
       IT_0183 + -IT_0184 + -IT_0185 + -IT_0186;
    const ccomplex_t IT_0188 = 18*s_34;
    const ccomplex_t IT_0189 = m_t*s_23;
    const ccomplex_t IT_0190 = (-18)*conj(IT_0139);
    const ccomplex_t IT_0191 = (-6)*conj(IT_0014);
    const ccomplex_t IT_0192 = IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = (-18)*IT_0139;
    const ccomplex_t IT_0194 = (-0.166666666666667)*((-18)*conj(IT_0014) + (-6
      )*conj(IT_0139))*IT_0161 + IT_0139*conj(IT_0161) + (-0.166666666666667)
      *IT_0014*((-18)*conj(IT_0161) + (-6)*conj(IT_0162)) + 3*conj(IT_0160)
      *IT_0187 + 3*IT_0160*conj(IT_0187) + (-0.166666666666667)*IT_0162*IT_0192 
      + (-0.166666666666667)*conj(IT_0162)*IT_0193;
    const ccomplex_t IT_0195 = m_c*m_t;
    const ccomplex_t IT_0196 = 12*IT_0166;
    const ccomplex_t IT_0197 = (-6)*IT_0164;
    const ccomplex_t IT_0198 = 6*s_34*(conj(IT_0014)*IT_0139 + IT_0014*conj
      (IT_0139)) + (IT_0160*conj(IT_0160) + IT_0161*conj(IT_0161) + IT_0162*conj
      (IT_0162))*IT_0168 + (-6)*IT_0169*(conj(IT_0160)*IT_0171 + IT_0160*IT_0173
      ) + 6*IT_0174*((conj(IT_0014) + 3*conj(IT_0139))*IT_0160 + (IT_0014 + 3
      *IT_0139)*conj(IT_0160) + IT_0173*IT_0187 + IT_0171*conj(IT_0187)) + 
      (IT_0014*conj(IT_0014) + IT_0139*conj(IT_0139) + IT_0187*conj(IT_0187))
      *IT_0188 + (-6)*IT_0189*IT_0194 + (IT_0187*IT_0192 + (-6)*conj(IT_0187)*
      (IT_0014 + (-0.166666666666667)*IT_0193))*IT_0195 + (conj(IT_0161)*IT_0162
       + IT_0161*conj(IT_0162))*(IT_0196 + IT_0197);
    return create_ccomplex_return(IT_0198);
}

