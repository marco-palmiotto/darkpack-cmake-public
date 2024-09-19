#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_ss_L_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_ss_L_to_c_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0007;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0011 + -IT_0013 + (-2)
      *IT_0015);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = pow(m_c, 2);
    const ccomplex_t IT_0020 = pow(m_ss_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*3)*e_em*(IT_0004*IT_0007 + 1./3
      *IT_0003*IT_0008);
    const ccomplex_t IT_0025 = 0.333333333333333*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = pow(m_C_2, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0019 + IT_0028 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0025*IT_0031;
    const ccomplex_t IT_0033 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = IT_0012*IT_0033;
    const ccomplex_t IT_0036 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0037 = IT_0012*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + -IT_0035 + (-2)
      *IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0027*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = 2*IT_0023 + IT_0032 + 2*IT_0043;
    const ccomplex_t IT_0045 = pow(m_Z, -2);
    const ccomplex_t IT_0046 = pow(s_14, 2);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = -IT_0028;
    const ccomplex_t IT_0049 = IT_0047 + IT_0048;
    const ccomplex_t IT_0050 = s_13*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + IT_0019 + -IT_0020 + -IT_0028 
      + -reg_prop, -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0008;
    const ccomplex_t IT_0053 = (-0.666666666666667)*IT_0052;
    const ccomplex_t IT_0054 = IT_0027*IT_0053;
    const ccomplex_t IT_0055 = IT_0051*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = -IT_0032 + (-2)*IT_0056;
    const ccomplex_t IT_0058 = s_14*s_34*IT_0045;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = s_13 + IT_0059;
    const ccomplex_t IT_0061 = s_13*IT_0060;
    const ccomplex_t IT_0062 = (-6)*conj(IT_0044);
    const ccomplex_t IT_0063 = (-6)*IT_0044;
    const ccomplex_t IT_0064 = pow(s_34, 2);
    const ccomplex_t IT_0065 = IT_0045*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0019 + IT_0066;
    const ccomplex_t IT_0068 = s_13*IT_0067;
    const ccomplex_t IT_0069 = (-6)*IT_0057;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0004;
    const ccomplex_t IT_0071 = IT_0029*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0025*IT_0072;
    const ccomplex_t IT_0074 = e_em*IT_0009;
    const ccomplex_t IT_0075 = e_em*IT_0012;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0074 + (-3)*IT_0075);
    const ccomplex_t IT_0077 = (-0.166666666666667)*IT_0076;
    const ccomplex_t IT_0078 = IT_0070*IT_0077;
    const ccomplex_t IT_0079 = IT_0051*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0073 + 2*IT_0080;
    const ccomplex_t IT_0082 = (-6)*IT_0081;
    const ccomplex_t IT_0083 = m_c*m_C_2;
    const ccomplex_t IT_0084 = IT_0060*IT_0083;
    const ccomplex_t IT_0085 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0009*IT_0085;
    const ccomplex_t IT_0087 = IT_0012*IT_0085;
    const ccomplex_t IT_0088 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0089 = IT_0012*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + -IT_0087 + (-2)
      *IT_0089);
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0070*IT_0091;
    const ccomplex_t IT_0093 = IT_0041*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0096 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0097 = IT_0009*IT_0096;
    const ccomplex_t IT_0098 = IT_0012*IT_0096;
    const ccomplex_t IT_0099 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0100 = IT_0012*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + -IT_0098 + (-2)
      *IT_0100);
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0095*IT_0102;
    const ccomplex_t IT_0104 = IT_0021*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = -IT_0073 + 2*IT_0094 + 2*IT_0105;
    const ccomplex_t IT_0107 = (-6)*IT_0106;
    const ccomplex_t IT_0108 = IT_0067*IT_0083;
    const ccomplex_t IT_0109 = m_C_2*IT_0067;
    const ccomplex_t IT_0110 = IT_0027*IT_0077;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0112 = IT_0051*IT_0110*IT_0111;
    const ccomplex_t IT_0113 = m_C_2*IT_0094;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0115 = IT_0039*IT_0070;
    const ccomplex_t IT_0116 = IT_0041*IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0118 = IT_0017*IT_0095;
    const ccomplex_t IT_0119 = IT_0021*IT_0117*IT_0118;
    const ccomplex_t IT_0120 = m_C_2*IT_0105;
    const ccomplex_t IT_0121 = m_c*IT_0056;
    const ccomplex_t IT_0122 = -IT_0112 + -IT_0113 + -IT_0116 + -IT_0119 + 
      -IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = (-6)*IT_0122;
    const ccomplex_t IT_0124 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0125 = (-6)*conj(IT_0057);
    const ccomplex_t IT_0126 = (-6)*conj(IT_0122);
    const ccomplex_t IT_0127 = IT_0027*IT_0091;
    const ccomplex_t IT_0128 = IT_0041*IT_0114*IT_0127;
    const ccomplex_t IT_0129 = IT_0053*IT_0070;
    const ccomplex_t IT_0130 = IT_0051*IT_0111*IT_0129;
    const ccomplex_t IT_0131 = m_c*IT_0080;
    const ccomplex_t IT_0132 = m_C_2*IT_0043;
    const ccomplex_t IT_0133 = m_C_2*IT_0023;
    const ccomplex_t IT_0134 = IT_0006*IT_0102;
    const ccomplex_t IT_0135 = IT_0021*IT_0117*IT_0134;
    const ccomplex_t IT_0136 = -IT_0128 + IT_0130 + -IT_0131 + -IT_0132 + 
      -IT_0133 + -IT_0135;
    const ccomplex_t IT_0137 = pow(m_Z, 2);
    const ccomplex_t IT_0138 = s_14*IT_0137;
    const ccomplex_t IT_0139 = IT_0045*IT_0138;
    const ccomplex_t IT_0140 = (-0.25)*IT_0139;
    const ccomplex_t IT_0141 = s_14 + IT_0140;
    const ccomplex_t IT_0142 = m_c*IT_0141;
    const ccomplex_t IT_0143 = (-4)*IT_0142;
    const ccomplex_t IT_0144 = -IT_0094 + -IT_0105;
    const ccomplex_t IT_0145 = 0.5*IT_0139;
    const ccomplex_t IT_0146 = s_14 + IT_0145;
    const ccomplex_t IT_0147 = m_c*IT_0146;
    const ccomplex_t IT_0148 = 2*IT_0147;
    const ccomplex_t IT_0149 = -IT_0056;
    const ccomplex_t IT_0150 = s_34*IT_0137;
    const ccomplex_t IT_0151 = IT_0045*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = s_34 + IT_0152;
    const ccomplex_t IT_0154 = m_C_2*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = -IT_0023 + -IT_0043;
    const ccomplex_t IT_0157 = (-0.25)*IT_0151;
    const ccomplex_t IT_0158 = s_34 + IT_0157;
    const ccomplex_t IT_0159 = m_C_2*IT_0158;
    const ccomplex_t IT_0160 = (-4)*IT_0159;
    const ccomplex_t IT_0161 = 12*s_13;
    const ccomplex_t IT_0162 = s_13*IT_0137;
    const ccomplex_t IT_0163 = s_14*s_34;
    const ccomplex_t IT_0164 = (-2)*IT_0163;
    const ccomplex_t IT_0165 = IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = IT_0045*IT_0165;
    const ccomplex_t IT_0167 = (-6)*IT_0166;
    const ccomplex_t IT_0168 = IT_0161 + IT_0167;
    const ccomplex_t IT_0169 = m_c*IT_0060;
    const ccomplex_t IT_0170 = m_C_2*IT_0060;
    const ccomplex_t IT_0171 = m_c*IT_0049;
    const ccomplex_t IT_0172 = 6*IT_0106;
    const ccomplex_t IT_0173 = IT_0069*IT_0109 + (-18)*IT_0083*IT_0122 + 6
      *IT_0080*IT_0143 + 6*IT_0144*IT_0148 + 6*IT_0149*IT_0155 + 6*IT_0156
      *IT_0160 + IT_0136*IT_0168 + IT_0082*IT_0169 + IT_0063*IT_0170 + IT_0171
      *IT_0172;
    const ccomplex_t IT_0174 = (-6)*conj(IT_0081);
    const ccomplex_t IT_0175 = 6*conj(IT_0106);
    const ccomplex_t IT_0176 = 6*conj(IT_0156);
    const ccomplex_t IT_0177 = 6*conj(IT_0144);
    const ccomplex_t IT_0178 = 6*conj(IT_0080);
    const ccomplex_t IT_0179 = IT_0136*((-18)*IT_0083*conj(IT_0122) + IT_0109
      *IT_0125 + 6*conj(IT_0149)*IT_0155 + IT_0062*IT_0170 + IT_0169*IT_0174 +
       IT_0171*IT_0175 + IT_0160*IT_0176 + IT_0148*IT_0177 + IT_0143*IT_0178);
    const ccomplex_t IT_0180 = IT_0083*IT_0137;
    const ccomplex_t IT_0181 = (-3)*IT_0180;
    const ccomplex_t IT_0182 = 3*IT_0180;
    const ccomplex_t IT_0183 = pow(m_Z, 4);
    const ccomplex_t IT_0184 = s_13*IT_0183;
    const ccomplex_t IT_0185 = IT_0045*IT_0184;
    const ccomplex_t IT_0186 = (-24)*IT_0163;
    const ccomplex_t IT_0187 = 6*IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0189 = 24*IT_0163;
    const ccomplex_t IT_0190 = (-12)*IT_0162;
    const ccomplex_t IT_0191 = 6*IT_0185 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 0.166666666666667*IT_0156;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0144;
    const ccomplex_t IT_0194 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0195 = s_14*IT_0019;
    const ccomplex_t IT_0196 = s_13*s_34;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = 0.166666666666667*conj(IT_0149);
    const ccomplex_t IT_0200 = 0.166666666666667*conj(IT_0080);
    const ccomplex_t IT_0201 = s_13*s_14;
    const ccomplex_t IT_0202 = s_34*IT_0028;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = IT_0049*IT_0083;
    const ccomplex_t IT_0206 = 0.166666666666667*IT_0084;
    const ccomplex_t IT_0207 = 0.166666666666667*IT_0170;
    const ccomplex_t IT_0208 = 0.166666666666667*IT_0169;
    const ccomplex_t IT_0209 = 6*IT_0044*conj(IT_0044)*IT_0050 + IT_0057
      *IT_0061*IT_0062 + conj(IT_0057)*(IT_0061*IT_0063 + IT_0068*IT_0069) +
       conj(IT_0081)*(IT_0068*IT_0082 + IT_0063*IT_0084 + IT_0061*IT_0107 +
       IT_0069*IT_0108 + IT_0109*IT_0123) + IT_0081*(IT_0062*IT_0084 + IT_0061
      *IT_0124 + IT_0108*IT_0125 + IT_0109*IT_0126) + conj(IT_0136)*IT_0173 +
       IT_0179 + 6*conj(IT_0149)*(IT_0080*IT_0181 + IT_0144*IT_0182 + IT_0188
      *IT_0191 + IT_0187*IT_0192) + IT_0176*(IT_0144*IT_0181 + IT_0080*IT_0182 +
       IT_0187*IT_0188 + IT_0191*IT_0192) + IT_0177*(IT_0156*IT_0181 + IT_0149
      *IT_0182 + IT_0191*IT_0193 + IT_0187*IT_0194) + IT_0178*(IT_0149*IT_0181 +
       IT_0156*IT_0182 + IT_0187*IT_0193 + IT_0191*IT_0194) + 6*IT_0198*(conj
      (IT_0081)*IT_0144 + IT_0081*conj(IT_0144) + conj(IT_0057)*IT_0156 +
       IT_0057*conj(IT_0156) + IT_0125*IT_0188 + IT_0174*IT_0194 + IT_0069
      *IT_0199 + IT_0082*IT_0200) + 6*(conj(IT_0106)*IT_0144 + IT_0106*conj
      (IT_0144) + conj(IT_0044)*IT_0156 + IT_0044*conj(IT_0156) + IT_0062
      *IT_0188 + IT_0124*IT_0194 + IT_0063*IT_0199 + IT_0107*IT_0200)*IT_0204 +
       IT_0175*(IT_0044*IT_0205 + IT_0069*IT_0206 + IT_0123*IT_0207) + IT_0172*
      (IT_0050*conj(IT_0106) + conj(IT_0044)*IT_0205 + IT_0125*IT_0206 + IT_0126
      *IT_0207) + 6*conj(IT_0122)*(IT_0143*IT_0149 + IT_0080*IT_0155 + IT_0148
      *IT_0156 + IT_0144*IT_0160 + IT_0044*IT_0171 + IT_0069*IT_0208) + 6
      *IT_0122*(IT_0143*conj(IT_0149) + conj(IT_0080)*IT_0155 + IT_0148*conj
      (IT_0156) + conj(IT_0144)*IT_0160 + 0.166666666666667*conj(IT_0122)
      *IT_0168 + conj(IT_0044)*IT_0171 + IT_0125*IT_0208);
    return create_ccomplex_return(IT_0209);
}

