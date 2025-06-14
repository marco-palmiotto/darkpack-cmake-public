#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_A_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.166666666666667*IT_0003;
    const ccomplex_t IT_0005 = pow(m_sc_L, 2);
    const ccomplex_t IT_0006 = pow(m_ss_L, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0005 + -IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0002;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0005 + 0.5*IT_0006;
    const ccomplex_t IT_0014 = s_12 + IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 2*IT_0011;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = s_13 + s_23;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = -IT_0009;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = s_23*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0020;
    const ccomplex_t IT_0024 = s_12 + IT_0006;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0019 + 0.5*IT_0022 + (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0015 + IT_0026;
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = pow(m_c, 2);
    const ccomplex_t IT_0031 = pow(m_s, 2);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, 2);
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, 2);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0002*(IT_0029*(IT_0030 + -IT_0031) + m_W*(IT_0033 + -IT_0035));
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = pow(m_W, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + -IT_0005 + -IT_0006 + IT_0040 
      + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0004 + (-2)*IT_0028 + -IT_0043;
    const ccomplex_t IT_0045 = pow(m_W, -2);
    const ccomplex_t IT_0046 = pow(s_14, 2);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0005 + IT_0048;
    const ccomplex_t IT_0050 = (-1.33333333333333)*IT_0011;
    const ccomplex_t IT_0051 = IT_0010*IT_0050;
    const ccomplex_t IT_0052 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_sc_L
      *Gamma_cl + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (-2)*IT_0011;
    const ccomplex_t IT_0056 = IT_0010*IT_0055;
    const ccomplex_t IT_0057 = IT_0012 + IT_0056;
    const ccomplex_t IT_0058 = IT_0008*IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0054 + -IT_0058;
    const ccomplex_t IT_0060 = 3*IT_0044;
    const ccomplex_t IT_0061 = 3*conj(IT_0044);
    const ccomplex_t IT_0062 = pow(s_24, 2);
    const ccomplex_t IT_0063 = IT_0045*IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0006 + IT_0064;
    const ccomplex_t IT_0066 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_ss_L
      *Gamma_sl + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0011;
    const ccomplex_t IT_0069 = IT_0010*IT_0068;
    const ccomplex_t IT_0070 = IT_0020*IT_0068;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = IT_0067*IT_0072;
    const ccomplex_t IT_0074 = IT_0020*IT_0055;
    const ccomplex_t IT_0075 = -IT_0023 + -IT_0074;
    const ccomplex_t IT_0076 = IT_0057 + IT_0075;
    const ccomplex_t IT_0077 = IT_0008*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0073 + -IT_0077;
    const ccomplex_t IT_0079 = s_14*s_24*IT_0045;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = s_12 + IT_0080;
    const ccomplex_t IT_0082 = -IT_0011;
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = IT_0020*IT_0082;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = (-0.666666666666667)*IT_0011;
    const ccomplex_t IT_0089 = IT_0020*IT_0088;
    const ccomplex_t IT_0090 = IT_0051 + IT_0089;
    const ccomplex_t IT_0091 = IT_0020*IT_0050;
    const ccomplex_t IT_0092 = IT_0010*IT_0088;
    const ccomplex_t IT_0093 = -IT_0091 + -IT_0092;
    const ccomplex_t IT_0094 = IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = IT_0053*IT_0094;
    const ccomplex_t IT_0096 = 0.333333333333333*IT_0011;
    const ccomplex_t IT_0097 = IT_0010*IT_0096;
    const ccomplex_t IT_0098 = IT_0020*IT_0096;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0097 + IT_0099;
    const ccomplex_t IT_0101 = IT_0067*IT_0100;
    const ccomplex_t IT_0102 = IT_0087 + 0.5*IT_0095 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = 3*s_13;
    const ccomplex_t IT_0104 = IT_0059*IT_0103;
    const ccomplex_t IT_0105 = -IT_0023;
    const ccomplex_t IT_0106 = IT_0057 + IT_0105;
    const ccomplex_t IT_0107 = IT_0008*IT_0106;
    const ccomplex_t IT_0108 = IT_0067*IT_0069;
    const ccomplex_t IT_0109 = -IT_0107 + 0.5*IT_0108;
    const ccomplex_t IT_0110 = 3*s_23;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = IT_0104 + IT_0111;
    const ccomplex_t IT_0113 = conj(IT_0059)*IT_0103;
    const ccomplex_t IT_0114 = conj(IT_0109)*IT_0110;
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = -IT_0074;
    const ccomplex_t IT_0117 = IT_0057 + IT_0116;
    const ccomplex_t IT_0118 = IT_0008*IT_0117;
    const ccomplex_t IT_0119 = -IT_0091;
    const ccomplex_t IT_0120 = IT_0051 + IT_0119;
    const ccomplex_t IT_0121 = IT_0053*IT_0120;
    const ccomplex_t IT_0122 = -IT_0118 + (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0067*IT_0097;
    const ccomplex_t IT_0124 = IT_0008*IT_0083;
    const ccomplex_t IT_0125 = -IT_0092;
    const ccomplex_t IT_0126 = IT_0051 + IT_0125;
    const ccomplex_t IT_0127 = IT_0053*IT_0126;
    const ccomplex_t IT_0128 = (-0.5)*IT_0123 + IT_0124 + 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0103*IT_0128;
    const ccomplex_t IT_0130 = 3*IT_0005;
    const ccomplex_t IT_0131 = IT_0059*IT_0130;
    const ccomplex_t IT_0132 = IT_0060 + IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = IT_0103*conj(IT_0128);
    const ccomplex_t IT_0134 = conj(IT_0059)*IT_0130;
    const ccomplex_t IT_0135 = IT_0061 + IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = IT_0110*IT_0128;
    const ccomplex_t IT_0137 = 3*IT_0006;
    const ccomplex_t IT_0138 = IT_0109*IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = IT_0110*conj(IT_0128);
    const ccomplex_t IT_0141 = conj(IT_0109)*IT_0137;
    const ccomplex_t IT_0142 = IT_0140 + IT_0141;
    const ccomplex_t IT_0143 = 3*s_12;
    const ccomplex_t IT_0144 = IT_0061*IT_0109 + IT_0060*conj(IT_0109) + conj
      (IT_0102)*IT_0112 + IT_0102*IT_0115 + conj(IT_0122)*IT_0132 + IT_0122
      *IT_0135 + conj(IT_0078)*IT_0139 + IT_0078*IT_0142 + (conj(IT_0059)
      *IT_0078 + IT_0059*conj(IT_0078) + conj(IT_0109)*IT_0122 + IT_0109*conj
      (IT_0122))*IT_0143;
    const ccomplex_t IT_0145 = s_14*s_34*IT_0045;
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = s_13 + IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0101 + (-0.5)*IT_0127;
    const ccomplex_t IT_0149 = 0.5*IT_0054 + IT_0124;
    const ccomplex_t IT_0150 = (-0.5)*IT_0073 + IT_0087;
    const ccomplex_t IT_0151 = IT_0061*IT_0128 + IT_0060*conj(IT_0128) +
       IT_0115*IT_0148 + IT_0112*conj(IT_0148) + IT_0135*IT_0149 + IT_0132*conj
      (IT_0149) + IT_0142*IT_0150 + IT_0139*conj(IT_0150) + IT_0143*(conj
      (IT_0109)*IT_0149 + IT_0109*conj(IT_0149) + conj(IT_0059)*IT_0150 +
       IT_0059*conj(IT_0150));
    const ccomplex_t IT_0152 = pow(s_34, 2);
    const ccomplex_t IT_0153 = IT_0045*IT_0152;
    const ccomplex_t IT_0154 = (-3)*s_13;
    const ccomplex_t IT_0155 = (-3)*s_23;
    const ccomplex_t IT_0156 = (-0.333333333333333)*conj(IT_0148);
    const ccomplex_t IT_0157 = IT_0044*conj(IT_0148) + IT_0148*(conj(IT_0044) 
      + (-0.333333333333333)*conj(IT_0149)*IT_0154 + (-0.333333333333333)*conj
      (IT_0150)*IT_0155) + IT_0149*(IT_0005*conj(IT_0149) + s_12*conj(IT_0150) +
       IT_0154*IT_0156) + IT_0150*(s_12*conj(IT_0149) + IT_0006*conj(IT_0150) +
       IT_0155*IT_0156);
    const ccomplex_t IT_0158 = IT_0103*conj(IT_0122);
    const ccomplex_t IT_0159 = conj(IT_0078)*IT_0110;
    const ccomplex_t IT_0160 = IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = IT_0103*IT_0122;
    const ccomplex_t IT_0162 = IT_0122*IT_0130;
    const ccomplex_t IT_0163 = conj(IT_0078)*IT_0137;
    const ccomplex_t IT_0164 = s_24*s_34*IT_0045;
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = s_23 + IT_0165;
    const ccomplex_t IT_0167 = conj(IT_0122)*IT_0130*IT_0149 + conj(IT_0102)*
      (IT_0060 + IT_0103*IT_0149 + IT_0110*IT_0150) + (IT_0060 + IT_0078*IT_0137
      )*conj(IT_0150) + IT_0102*(IT_0061 + IT_0103*conj(IT_0149) + IT_0110*conj
      (IT_0150)) + IT_0143*(conj(IT_0078)*IT_0149 + IT_0078*conj(IT_0149) + conj
      (IT_0122)*IT_0150 + IT_0122*conj(IT_0150)) + IT_0148*IT_0160 + conj
      (IT_0148)*(IT_0078*IT_0110 + IT_0161) + conj(IT_0149)*IT_0162 + IT_0150*
      (IT_0061 + IT_0163);
    const ccomplex_t IT_0168 = 9*IT_0044*conj(IT_0044) + IT_0049*(conj(IT_0059
      )*IT_0060 + IT_0059*IT_0061) + IT_0065*(IT_0061*IT_0078 + IT_0060*conj
      (IT_0078)) + IT_0049*(IT_0115*IT_0128 + IT_0112*conj(IT_0128) + conj
      (IT_0109)*(IT_0138 + IT_0059*IT_0143) + conj(IT_0059)*(IT_0131 + IT_0109
      *IT_0143)) + IT_0081*IT_0144 + IT_0147*IT_0151 + (-3)*IT_0153*IT_0157 +
       IT_0065*(conj(IT_0078)*IT_0122*IT_0143 + IT_0102*IT_0160 + conj(IT_0102)
      *IT_0161 + conj(IT_0122)*IT_0162 + IT_0078*(conj(IT_0102)*IT_0110 + conj
      (IT_0122)*IT_0143 + IT_0163)) + IT_0166*IT_0167;
    return create_ccomplex_return(IT_0168);
}

