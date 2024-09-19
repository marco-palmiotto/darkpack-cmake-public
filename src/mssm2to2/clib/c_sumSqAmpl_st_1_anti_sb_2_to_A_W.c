#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_2_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_2_to_A_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*V_tb*IT_0000
      *IT_0002*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0004 = 0.166666666666667*IT_0003;
    const ccomplex_t IT_0005 = pow(m_sb_2, 2);
    const ccomplex_t IT_0006 = pow(m_st_1, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0005 + -IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = (s_13 + s_23)*(s_12 + 0.5*IT_0005 + 0.5*IT_0006);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0002*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0014 = 2*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = IT_0009*IT_0010*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0005 + 0.5*IT_0006;
    const ccomplex_t IT_0018 = s_12 + IT_0017;
    const ccomplex_t IT_0019 = IT_0012*IT_0013;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (s_12 + IT_0005)*(s_12 + 0.5*IT_0005 + 0.5
      *IT_0006);
    const ccomplex_t IT_0022 = -IT_0011;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = IT_0009*IT_0021*IT_0023;
    const ccomplex_t IT_0025 = IT_0016 + IT_0020 + IT_0024;
    const ccomplex_t IT_0026 = s_12 + IT_0005;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = s_13 + s_23;
    const ccomplex_t IT_0029 = IT_0015*IT_0028;
    const ccomplex_t IT_0030 = IT_0014*IT_0022;
    const ccomplex_t IT_0031 = s_23*IT_0030;
    const ccomplex_t IT_0032 = (s_13 + s_23)*(s_12 + IT_0005);
    const ccomplex_t IT_0033 = IT_0009*IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow(s_12 + 0.5*IT_0005 + 0.5*IT_0006, 2);
    const ccomplex_t IT_0035 = IT_0009*IT_0019*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0027 + (-0.5)*IT_0029 + 0.5*IT_0031 +
       (-0.5)*IT_0033 + (-2)*IT_0035;
    const ccomplex_t IT_0037 = IT_0025 + IT_0036;
    const ccomplex_t IT_0038 = IT_0008*IT_0037;
    const ccomplex_t IT_0039 = IT_0004 + (-2)*IT_0038;
    const ccomplex_t IT_0040 = pow(s_24, 2);
    const ccomplex_t IT_0041 = IT_0009*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0005 + IT_0042;
    const ccomplex_t IT_0044 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_st_1
      *Gamma_t1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (-1.33333333333333)*IT_0013;
    const ccomplex_t IT_0047 = IT_0012*IT_0046;
    const ccomplex_t IT_0048 = (-0.666666666666667)*IT_0013;
    const ccomplex_t IT_0049 = IT_0022*IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = IT_0012*IT_0048;
    const ccomplex_t IT_0052 = IT_0022*IT_0046;
    const ccomplex_t IT_0053 = -IT_0051 + -IT_0052;
    const ccomplex_t IT_0054 = IT_0050 + IT_0053;
    const ccomplex_t IT_0055 = IT_0045*IT_0054;
    const ccomplex_t IT_0056 = 0.333333333333333*IT_0013;
    const ccomplex_t IT_0057 = IT_0012*IT_0056;
    const ccomplex_t IT_0058 = IT_0022*IT_0056;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = IT_0057 + IT_0059;
    const ccomplex_t IT_0061 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_sb_2
      *Gamma_b2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = -IT_0013;
    const ccomplex_t IT_0065 = IT_0022*IT_0064;
    const ccomplex_t IT_0066 = IT_0012*IT_0064;
    const ccomplex_t IT_0067 = IT_0009*IT_0018*IT_0066;
    const ccomplex_t IT_0068 = -IT_0065 + (-2)*IT_0067;
    const ccomplex_t IT_0069 = IT_0009*IT_0026*IT_0065;
    const ccomplex_t IT_0070 = IT_0066 + IT_0069;
    const ccomplex_t IT_0071 = IT_0068 + IT_0070;
    const ccomplex_t IT_0072 = IT_0008*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0055 + (-0.5)*IT_0063 + IT_0072;
    const ccomplex_t IT_0074 = -IT_0052;
    const ccomplex_t IT_0075 = IT_0047 + IT_0074;
    const ccomplex_t IT_0076 = IT_0045*IT_0075;
    const ccomplex_t IT_0077 = IT_0009*IT_0023*IT_0026;
    const ccomplex_t IT_0078 = (-2)*IT_0013;
    const ccomplex_t IT_0079 = IT_0012*IT_0078;
    const ccomplex_t IT_0080 = IT_0022*IT_0078;
    const ccomplex_t IT_0081 = IT_0009*IT_0026*IT_0080;
    const ccomplex_t IT_0082 = IT_0019 + IT_0077 + IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = IT_0009*IT_0018*IT_0019;
    const ccomplex_t IT_0084 = IT_0009*IT_0018*IT_0079;
    const ccomplex_t IT_0085 = -IT_0080 + (-2)*IT_0083 + (-2)*IT_0084;
    const ccomplex_t IT_0086 = IT_0082 + IT_0085;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0076 + -IT_0087;
    const ccomplex_t IT_0089 = 3*s_13;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = 0.666666666666667*IT_0013;
    const ccomplex_t IT_0092 = IT_0012*IT_0091;
    const ccomplex_t IT_0093 = IT_0022*IT_0091;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = IT_0062*IT_0095;
    const ccomplex_t IT_0097 = -IT_0023 + -IT_0080 + (-2)*IT_0083 + (-2)
      *IT_0084;
    const ccomplex_t IT_0098 = IT_0082 + IT_0097;
    const ccomplex_t IT_0099 = IT_0008*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0096 + -IT_0099;
    const ccomplex_t IT_0101 = 3*s_23;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = IT_0090 + IT_0102;
    const ccomplex_t IT_0104 = conj(IT_0088)*IT_0089;
    const ccomplex_t IT_0105 = conj(IT_0100)*IT_0101;
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = 3*s_12;
    const ccomplex_t IT_0108 = 3*IT_0005;
    const ccomplex_t IT_0109 = IT_0100*IT_0108;
    const ccomplex_t IT_0110 = 3*IT_0006;
    const ccomplex_t IT_0111 = IT_0088*IT_0110;
    const ccomplex_t IT_0112 = pow(s_14, 2);
    const ccomplex_t IT_0113 = IT_0009*IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0006 + IT_0114;
    const ccomplex_t IT_0116 = (-2)*IT_0083 + (-2)*IT_0084;
    const ccomplex_t IT_0117 = IT_0082 + IT_0116;
    const ccomplex_t IT_0118 = IT_0008*IT_0117;
    const ccomplex_t IT_0119 = IT_0045*IT_0047;
    const ccomplex_t IT_0120 = -IT_0118 + (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = 3*IT_0039;
    const ccomplex_t IT_0122 = 3*conj(IT_0039);
    const ccomplex_t IT_0123 = s_24*s_34*IT_0009;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = s_23 + IT_0124;
    const ccomplex_t IT_0126 = -IT_0051;
    const ccomplex_t IT_0127 = IT_0047 + IT_0126;
    const ccomplex_t IT_0128 = IT_0045*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0063 + (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = (-2)*IT_0067;
    const ccomplex_t IT_0131 = IT_0070 + IT_0130;
    const ccomplex_t IT_0132 = IT_0008*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0119 + IT_0132;
    const ccomplex_t IT_0134 = IT_0072 + (-0.5)*IT_0096;
    const ccomplex_t IT_0135 = IT_0106*IT_0129 + IT_0103*conj(IT_0129) + conj
      (IT_0088)*IT_0110*IT_0133 + IT_0111*conj(IT_0133) + (conj(IT_0100)*IT_0108
       + IT_0122)*IT_0134 + conj(IT_0073)*(IT_0121 + IT_0089*IT_0133 + IT_0101
      *IT_0134) + (IT_0109 + IT_0121)*conj(IT_0134) + IT_0107*(conj(IT_0100)
      *IT_0133 + IT_0100*conj(IT_0133) + conj(IT_0088)*IT_0134 + IT_0088*conj
      (IT_0134)) + IT_0073*(IT_0122 + IT_0089*conj(IT_0133) + IT_0101*conj
      (IT_0134));
    const ccomplex_t IT_0136 = pow(s_34, 2);
    const ccomplex_t IT_0137 = IT_0009*IT_0136;
    const ccomplex_t IT_0138 = (-3)*s_13;
    const ccomplex_t IT_0139 = (-3)*s_23;
    const ccomplex_t IT_0140 = (-0.333333333333333)*conj(IT_0129);
    const ccomplex_t IT_0141 = IT_0039*conj(IT_0129) + IT_0129*(conj(IT_0039) 
      + (-0.333333333333333)*conj(IT_0133)*IT_0138 + (-0.333333333333333)*conj
      (IT_0134)*IT_0139) + IT_0133*(IT_0006*conj(IT_0133) + s_12*conj(IT_0134) +
       IT_0138*IT_0140) + IT_0134*(s_12*conj(IT_0133) + IT_0005*conj(IT_0134) +
       IT_0139*IT_0140);
    const ccomplex_t IT_0142 = IT_0057*IT_0062;
    const ccomplex_t IT_0143 = 0.5*IT_0128 + IT_0132 + (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0089*IT_0120;
    const ccomplex_t IT_0145 = -IT_0023 + (-2)*IT_0083 + (-2)*IT_0084;
    const ccomplex_t IT_0146 = IT_0082 + IT_0145;
    const ccomplex_t IT_0147 = IT_0008*IT_0146;
    const ccomplex_t IT_0148 = IT_0062*IT_0092;
    const ccomplex_t IT_0149 = -IT_0147 + 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0101*IT_0149;
    const ccomplex_t IT_0151 = IT_0144 + IT_0150;
    const ccomplex_t IT_0152 = IT_0089*conj(IT_0120);
    const ccomplex_t IT_0153 = IT_0101*conj(IT_0149);
    const ccomplex_t IT_0154 = IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = IT_0108*IT_0149;
    const ccomplex_t IT_0156 = IT_0107*IT_0120;
    const ccomplex_t IT_0157 = IT_0107*IT_0149;
    const ccomplex_t IT_0158 = IT_0110*IT_0120;
    const ccomplex_t IT_0159 = s_14*s_24*IT_0009;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = s_12 + IT_0160;
    const ccomplex_t IT_0162 = IT_0101*IT_0143;
    const ccomplex_t IT_0163 = IT_0155 + IT_0156 + IT_0162;
    const ccomplex_t IT_0164 = IT_0101*conj(IT_0143);
    const ccomplex_t IT_0165 = IT_0108*conj(IT_0149);
    const ccomplex_t IT_0166 = IT_0107*conj(IT_0120) + IT_0164 + IT_0165;
    const ccomplex_t IT_0167 = IT_0089*conj(IT_0143);
    const ccomplex_t IT_0168 = IT_0110*conj(IT_0120);
    const ccomplex_t IT_0169 = IT_0122 + IT_0107*conj(IT_0149) + IT_0167 +
       IT_0168;
    const ccomplex_t IT_0170 = IT_0089*IT_0143;
    const ccomplex_t IT_0171 = s_14*s_34*IT_0009;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = s_13 + IT_0172;
    const ccomplex_t IT_0174 = 9*IT_0039*conj(IT_0039) + IT_0043*(conj(IT_0073
      )*IT_0103 + IT_0073*IT_0106 + conj(IT_0100)*(IT_0088*IT_0107 + IT_0109) +
       conj(IT_0088)*(IT_0100*IT_0107 + IT_0111)) + IT_0043*(conj(IT_0100)
      *IT_0121 + IT_0100*IT_0122) + IT_0115*(conj(IT_0120)*IT_0121 + IT_0120
      *IT_0122) + IT_0125*IT_0135 + (-3)*IT_0137*IT_0141 + IT_0115*(conj(IT_0143
      )*IT_0151 + IT_0143*IT_0154 + conj(IT_0149)*(IT_0155 + IT_0156) + conj
      (IT_0120)*(IT_0157 + IT_0158)) + IT_0161*(IT_0122*IT_0149 + IT_0121*conj
      (IT_0149) + conj(IT_0073)*IT_0151 + IT_0073*IT_0154 + conj(IT_0100)
      *IT_0163 + IT_0100*IT_0166 + IT_0088*IT_0169 + conj(IT_0088)*(IT_0121 +
       IT_0157 + IT_0158 + IT_0170)) + (IT_0122*IT_0143 + IT_0121*conj(IT_0143) 
      + conj(IT_0129)*IT_0151 + IT_0129*IT_0154 + conj(IT_0134)*IT_0163 +
       IT_0134*IT_0166 + IT_0133*IT_0169 + conj(IT_0133)*(IT_0121 + IT_0157 +
       IT_0158 + IT_0170))*IT_0173;
    return create_ccomplex_return(IT_0174);
}

