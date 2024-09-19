#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_ss_L_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_ss_L_to_c_Z(
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0003;
    const ccomplex_t IT_0008 = e_em*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0006 + (-3)*IT_0008);
    const ccomplex_t IT_0010 = (-0.166666666666667)*IT_0009;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(m_C_1, 2);
    const ccomplex_t IT_0014 = pow(m_ss_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + IT_0012 + -IT_0013 + -IT_0014 
      + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_Z, -2);
    const ccomplex_t IT_0019 = pow(m_Z, 4);
    const ccomplex_t IT_0020 = s_13*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 6*IT_0021;
    const ccomplex_t IT_0023 = s_14*s_34;
    const ccomplex_t IT_0024 = 24*IT_0023;
    const ccomplex_t IT_0025 = pow(m_Z, 2);
    const ccomplex_t IT_0026 = s_13*IT_0025;
    const ccomplex_t IT_0027 = (-12)*IT_0026;
    const ccomplex_t IT_0028 = IT_0022 + IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0030 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0031 = IT_0005*IT_0030;
    const ccomplex_t IT_0032 = IT_0007*IT_0030;
    const ccomplex_t IT_0033 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0034 = IT_0007*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + -IT_0032 + (-2)
      *IT_0034);
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0012 + IT_0014 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0042 = IT_0005*IT_0041;
    const ccomplex_t IT_0043 = IT_0007*IT_0041;
    const ccomplex_t IT_0044 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + -IT_0043 + (-2)
      *IT_0045);
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0002*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0012 + IT_0014 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0040 + -IT_0051;
    const ccomplex_t IT_0053 = (-24)*IT_0023;
    const ccomplex_t IT_0054 = IT_0022 + IT_0053;
    const ccomplex_t IT_0055 = m_c*m_C_1;
    const ccomplex_t IT_0056 = IT_0025*IT_0055;
    const ccomplex_t IT_0057 = (-3)*IT_0056;
    const ccomplex_t IT_0058 = pow(m_W, -1);
    const ccomplex_t IT_0059 = sin(beta);
    const ccomplex_t IT_0060 = cpow(IT_0059, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0058*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0004;
    const ccomplex_t IT_0064 = (-0.666666666666667)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = IT_0015*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = s_34*IT_0025;
    const ccomplex_t IT_0070 = IT_0018*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = s_34 + IT_0071;
    const ccomplex_t IT_0073 = m_C_1*IT_0072;
    const ccomplex_t IT_0074 = 2*IT_0073;
    const ccomplex_t IT_0075 = IT_0010*IT_0062;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0077 = IT_0015*IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0079 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0080 = IT_0005*IT_0079;
    const ccomplex_t IT_0081 = IT_0007*IT_0079;
    const ccomplex_t IT_0082 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0083 = IT_0007*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + -IT_0081 + (-2)
      *IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0002*IT_0085;
    const ccomplex_t IT_0087 = IT_0049*IT_0078*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0089 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0090 = IT_0005*IT_0089;
    const ccomplex_t IT_0091 = IT_0007*IT_0089;
    const ccomplex_t IT_0092 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + -IT_0091 + (-2)
      *IT_0093);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0029*IT_0095;
    const ccomplex_t IT_0097 = IT_0038*IT_0088*IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0040;
    const ccomplex_t IT_0099 = m_C_1*IT_0051;
    const ccomplex_t IT_0100 = m_c*IT_0067;
    const ccomplex_t IT_0101 = -IT_0077 + -IT_0087 + -IT_0097 + -IT_0098 + 
      -IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = s_14*IT_0025;
    const ccomplex_t IT_0103 = IT_0018*IT_0102;
    const ccomplex_t IT_0104 = (-0.25)*IT_0103;
    const ccomplex_t IT_0105 = s_14 + IT_0104;
    const ccomplex_t IT_0106 = m_c*IT_0105;
    const ccomplex_t IT_0107 = (-4)*IT_0106;
    const ccomplex_t IT_0108 = IT_0002*IT_0064;
    const ccomplex_t IT_0109 = IT_0015*IT_0076*IT_0108;
    const ccomplex_t IT_0110 = IT_0047*IT_0062;
    const ccomplex_t IT_0111 = IT_0049*IT_0078*IT_0110;
    const ccomplex_t IT_0112 = IT_0062*IT_0085;
    const ccomplex_t IT_0113 = IT_0049*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = m_C_1*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0001*IT_0058*IT_0060;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0036*IT_0117;
    const ccomplex_t IT_0119 = IT_0038*IT_0088*IT_0118;
    const ccomplex_t IT_0120 = m_c*IT_0017;
    const ccomplex_t IT_0121 = IT_0095*IT_0117;
    const ccomplex_t IT_0122 = IT_0038*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_C_1*IT_0123;
    const ccomplex_t IT_0125 = IT_0109 + -IT_0111 + -IT_0115 + -IT_0119 + 
      -IT_0120 + -IT_0124;
    const ccomplex_t IT_0126 = 3*IT_0056;
    const ccomplex_t IT_0127 = -IT_0114 + -IT_0123;
    const ccomplex_t IT_0128 = 6*conj(IT_0017);
    const ccomplex_t IT_0129 = (0 + _Complex_I*3)*e_em*(IT_0001*IT_0003 + 1./3
      *IT_0000*IT_0004);
    const ccomplex_t IT_0130 = 0.333333333333333*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0130*IT_0133;
    const ccomplex_t IT_0135 = 2*IT_0040 + 2*IT_0051 + -IT_0134;
    const ccomplex_t IT_0136 = s_13*s_14;
    const ccomplex_t IT_0137 = s_34*IT_0013;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = 6*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0017 + IT_0134;
    const ccomplex_t IT_0142 = s_14*IT_0012;
    const ccomplex_t IT_0143 = s_13*s_34;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0142 + IT_0144;
    const ccomplex_t IT_0146 = 6*IT_0145;
    const ccomplex_t IT_0147 = 6*IT_0057;
    const ccomplex_t IT_0148 = 6*IT_0126;
    const ccomplex_t IT_0149 = 0.5*IT_0103;
    const ccomplex_t IT_0150 = s_14 + IT_0149;
    const ccomplex_t IT_0151 = m_c*IT_0150;
    const ccomplex_t IT_0152 = 2*IT_0151;
    const ccomplex_t IT_0153 = 6*IT_0152;
    const ccomplex_t IT_0154 = (-0.25)*IT_0070;
    const ccomplex_t IT_0155 = s_34 + IT_0154;
    const ccomplex_t IT_0156 = m_C_1*IT_0155;
    const ccomplex_t IT_0157 = (-4)*IT_0156;
    const ccomplex_t IT_0158 = 6*IT_0157;
    const ccomplex_t IT_0159 = 6*IT_0017;
    const ccomplex_t IT_0160 = IT_0062*IT_0131;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0130*IT_0161;
    const ccomplex_t IT_0163 = 2*IT_0114 + 2*IT_0123 + IT_0162;
    const ccomplex_t IT_0164 = (-2)*IT_0067 + -IT_0162;
    const ccomplex_t IT_0165 = (-6)*conj(IT_0017);
    const ccomplex_t IT_0166 = 6*IT_0074;
    const ccomplex_t IT_0167 = 6*IT_0107;
    const ccomplex_t IT_0168 = (-6)*IT_0139;
    const ccomplex_t IT_0169 = (-6)*IT_0145;
    const ccomplex_t IT_0170 = pow(s_14, 2);
    const ccomplex_t IT_0171 = IT_0018*IT_0170;
    const ccomplex_t IT_0172 = -IT_0013;
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = s_13*IT_0173;
    const ccomplex_t IT_0175 = 6*IT_0174;
    const ccomplex_t IT_0176 = m_c*IT_0173;
    const ccomplex_t IT_0177 = 6*IT_0176;
    const ccomplex_t IT_0178 = IT_0055*IT_0173;
    const ccomplex_t IT_0179 = 6*IT_0178;
    const ccomplex_t IT_0180 = s_14*s_34*IT_0018;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = s_13 + IT_0181;
    const ccomplex_t IT_0183 = s_13*IT_0182;
    const ccomplex_t IT_0184 = (-6)*IT_0183;
    const ccomplex_t IT_0185 = IT_0055*IT_0182;
    const ccomplex_t IT_0186 = (-6)*IT_0185;
    const ccomplex_t IT_0187 = m_C_1*IT_0182;
    const ccomplex_t IT_0188 = (-6)*IT_0187;
    const ccomplex_t IT_0189 = 12*s_13;
    const ccomplex_t IT_0190 = (-2)*IT_0023;
    const ccomplex_t IT_0191 = IT_0026 + IT_0190;
    const ccomplex_t IT_0192 = IT_0018*IT_0191;
    const ccomplex_t IT_0193 = (-6)*IT_0192;
    const ccomplex_t IT_0194 = IT_0189 + IT_0193;
    const ccomplex_t IT_0195 = (-18)*IT_0055;
    const ccomplex_t IT_0196 = m_c*IT_0182;
    const ccomplex_t IT_0197 = (-6)*IT_0196;
    const ccomplex_t IT_0198 = pow(s_34, 2);
    const ccomplex_t IT_0199 = IT_0018*IT_0198;
    const ccomplex_t IT_0200 = -IT_0199;
    const ccomplex_t IT_0201 = IT_0012 + IT_0200;
    const ccomplex_t IT_0202 = m_C_1*IT_0201;
    const ccomplex_t IT_0203 = (-6)*IT_0202;
    const ccomplex_t IT_0204 = conj(IT_0125)*(IT_0052*IT_0153 + IT_0127
      *IT_0158 + IT_0107*IT_0159 + IT_0068*IT_0166 + IT_0135*IT_0177 + IT_0163
      *IT_0188 + IT_0125*IT_0194 + IT_0101*IT_0195 + IT_0141*IT_0197 + IT_0164
      *IT_0203);
    const ccomplex_t IT_0205 = conj(IT_0101)*(IT_0127*IT_0153 + IT_0052
      *IT_0158 + IT_0074*IT_0159 + IT_0068*IT_0167 + IT_0163*IT_0177 + IT_0135
      *IT_0188 + IT_0101*IT_0194 + IT_0125*IT_0195 + IT_0164*IT_0197 + IT_0141
      *IT_0203);
    const ccomplex_t IT_0206 = s_13*IT_0201;
    const ccomplex_t IT_0207 = (-6)*IT_0206;
    const ccomplex_t IT_0208 = IT_0055*IT_0201;
    const ccomplex_t IT_0209 = (-6)*IT_0208;
    const ccomplex_t IT_0210 = (-0.166666666666667)*IT_0052;
    const ccomplex_t IT_0211 = IT_0017*IT_0139 + (-0.166666666666667)*IT_0135
      *IT_0175 + (-0.166666666666667)*IT_0125*IT_0177 + (-0.166666666666667)
      *IT_0163*IT_0179 + (-0.166666666666667)*IT_0141*IT_0184 + (
      -0.166666666666667)*IT_0164*IT_0186 + (-0.166666666666667)*IT_0101*IT_0188
       + IT_0140*IT_0210;
    const ccomplex_t IT_0212 = IT_0017*IT_0145 + (-0.166666666666667)*IT_0135
      *IT_0184 + (-0.166666666666667)*IT_0163*IT_0186 + (-0.166666666666667)
      *IT_0125*IT_0197 + (-0.166666666666667)*IT_0101*IT_0203 + (
      -0.166666666666667)*IT_0141*IT_0207 + (-0.166666666666667)*IT_0164*IT_0209
       + IT_0146*IT_0210;
    const ccomplex_t IT_0213 = conj(IT_0017)*(IT_0017*IT_0028 + IT_0052
      *IT_0054) + (IT_0057*IT_0068 + IT_0074*IT_0101 + IT_0107*IT_0125 + IT_0126
      *IT_0127)*IT_0128 + conj(IT_0052)*(IT_0028*IT_0052 + IT_0017*IT_0054 +
       IT_0135*IT_0140 + IT_0141*IT_0146 + IT_0127*IT_0147 + IT_0068*IT_0148 +
       IT_0125*IT_0153 + IT_0101*IT_0158) + conj(IT_0127)*(IT_0054*IT_0068 +
       IT_0028*IT_0127 + IT_0052*IT_0147 + IT_0101*IT_0153 + IT_0125*IT_0158 +
       IT_0126*IT_0159 + IT_0140*IT_0163 + IT_0146*IT_0164) + (IT_0135*IT_0139 +
       IT_0141*IT_0145)*IT_0165 + conj(IT_0068)*(IT_0028*IT_0068 + IT_0054
      *IT_0127 + IT_0052*IT_0148 + IT_0057*IT_0159 + IT_0125*IT_0166 + IT_0101
      *IT_0167 + IT_0163*IT_0168 + IT_0164*IT_0169) + conj(IT_0163)*(IT_0127
      *IT_0140 + IT_0068*IT_0168 + IT_0163*IT_0175 + IT_0101*IT_0177 + IT_0135
      *IT_0179 + IT_0164*IT_0184 + IT_0141*IT_0186 + IT_0125*IT_0188) + IT_0204 
      + IT_0205 + conj(IT_0164)*(IT_0127*IT_0146 + IT_0068*IT_0169 + IT_0163
      *IT_0184 + IT_0135*IT_0186 + IT_0101*IT_0197 + IT_0125*IT_0203 + IT_0164
      *IT_0207 + IT_0141*IT_0209) + (-6)*conj(IT_0135)*IT_0211 + (-6)*conj
      (IT_0141)*IT_0212;
    return create_ccomplex_return(IT_0213);
}

