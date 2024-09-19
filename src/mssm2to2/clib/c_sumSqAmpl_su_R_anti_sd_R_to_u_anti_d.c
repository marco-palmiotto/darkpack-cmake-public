#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_R_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_R_to_u_anti_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_R = param->m_su_R;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_su_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = s_34*IT_0002;
    const ccomplex_t IT_0010 = (-6)*IT_0009;
    const ccomplex_t IT_0011 = s_13*s_14;
    const ccomplex_t IT_0012 = 12*IT_0011;
    const ccomplex_t IT_0013 = IT_0010 + IT_0012;
    const ccomplex_t IT_0014 = 0.166666666666667*IT_0007;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0016;
    const ccomplex_t IT_0018 = (-0.333333333333333)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0016;
    const ccomplex_t IT_0020 = 0.666666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0016;
    const ccomplex_t IT_0026 = (-0.333333333333333)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0016;
    const ccomplex_t IT_0028 = 0.666666666666667*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0016;
    const ccomplex_t IT_0034 = (-0.333333333333333)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0016;
    const ccomplex_t IT_0036 = 0.666666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0016;
    const ccomplex_t IT_0042 = (-0.333333333333333)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0016;
    const ccomplex_t IT_0044 = 0.666666666666667*IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = -IT_0024 + -IT_0032 + -IT_0040 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0049 + -conj(IT_0049);
    const ccomplex_t IT_0051 = (-18)*IT_0009;
    const ccomplex_t IT_0052 = 36*IT_0011;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = pow(m_W, -1);
    const ccomplex_t IT_0055 = sin(beta);
    const ccomplex_t IT_0056 = cpow(IT_0055, -1);
    const ccomplex_t IT_0057 = sin(theta_W);
    const ccomplex_t IT_0058 = cpow(IT_0057, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0054*IT_0056*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cos(beta);
    const ccomplex_t IT_0062 = cpow(IT_0061, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0054*IT_0058*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0060*IT_0064;
    const ccomplex_t IT_0066 = IT_0038*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0054*IT_0056*IT_0058;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0054*IT_0058*IT_0062;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = IT_0046*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0054*IT_0056*IT_0058;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0054*IT_0058*IT_0062;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = IT_0022*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0054*IT_0056*IT_0058;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0054*IT_0058*IT_0062;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0030*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = -IT_0067 + -IT_0074 + -IT_0081 + -IT_0088;
    const ccomplex_t IT_0090 = cpow(IT_0008, 2);
    const ccomplex_t IT_0091 = pow(m_sd_R, 2);
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + -IT_0002 + -IT_0091 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0054*IT_0056*IT_0058*IT_0062;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0054*IT_0056*IT_0058*IT_0061;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0094*IT_0096;
    const ccomplex_t IT_0098 = IT_0092*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0044*IT_0071;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0103 = IT_0046*IT_0101*IT_0102;
    const ccomplex_t IT_0104 = IT_0020*IT_0078;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0106 = IT_0022*IT_0104*IT_0105;
    const ccomplex_t IT_0107 = IT_0028*IT_0085;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0109 = IT_0030*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = m_u*IT_0067;
    const ccomplex_t IT_0111 = m_u*IT_0074;
    const ccomplex_t IT_0112 = m_u*IT_0081;
    const ccomplex_t IT_0113 = IT_0036*IT_0064;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0115 = IT_0038*IT_0113*IT_0114;
    const ccomplex_t IT_0116 = m_u*IT_0088;
    const ccomplex_t IT_0117 = -IT_0103 + -IT_0106 + -IT_0109 + IT_0110 +
       IT_0111 + IT_0112 + -IT_0115 + IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0054*IT_0055*IT_0058*IT_0062;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0094*IT_0119;
    const ccomplex_t IT_0121 = IT_0092*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0018*IT_0076;
    const ccomplex_t IT_0125 = IT_0022*IT_0105*IT_0124;
    const ccomplex_t IT_0126 = IT_0026*IT_0083;
    const ccomplex_t IT_0127 = IT_0030*IT_0108*IT_0126;
    const ccomplex_t IT_0128 = m_u*IT_0007;
    const ccomplex_t IT_0129 = m_u*IT_0024;
    const ccomplex_t IT_0130 = m_u*IT_0032;
    const ccomplex_t IT_0131 = m_u*IT_0040;
    const ccomplex_t IT_0132 = m_u*IT_0048;
    const ccomplex_t IT_0133 = IT_0034*IT_0060;
    const ccomplex_t IT_0134 = IT_0038*IT_0114*IT_0133;
    const ccomplex_t IT_0135 = IT_0042*IT_0069;
    const ccomplex_t IT_0136 = IT_0046*IT_0102*IT_0135;
    const ccomplex_t IT_0137 = -IT_0125 + -IT_0127 + (-0.166666666666667)
      *IT_0128 + IT_0129 + IT_0130 + IT_0131 + IT_0132 + -IT_0134 + -IT_0136;
    const ccomplex_t IT_0138 = 6*s_34;
    const ccomplex_t IT_0139 = 18*s_34;
    const ccomplex_t IT_0140 = 0.5*IT_0128;
    const ccomplex_t IT_0141 = (-18)*s_34;
    const ccomplex_t IT_0142 = m_d*m_u*IT_0002;
    const ccomplex_t IT_0143 = (-18)*conj(IT_0089);
    const ccomplex_t IT_0144 = m_u*s_14;
    const ccomplex_t IT_0145 = 6*IT_0100;
    const ccomplex_t IT_0146 = 18*IT_0117;
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = 6*IT_0140;
    const ccomplex_t IT_0149 = 6*conj(IT_0123);
    const ccomplex_t IT_0150 = 18*conj(IT_0137);
    const ccomplex_t IT_0151 = (-6)*IT_0140;
    const ccomplex_t IT_0152 = (-18)*IT_0137;
    const ccomplex_t IT_0153 = (-6)*IT_0123;
    const ccomplex_t IT_0154 = IT_0089*(conj(IT_0100) + 3*conj(IT_0117)) +
       0.166666666666667*IT_0008*((-18)*IT_0123 + 18*conj(IT_0123) + (-6)
      *IT_0137 + 6*conj(IT_0137) + (-36)*IT_0140) + 0.166666666666667*conj
      (IT_0089)*IT_0147 + conj(IT_0049)*(IT_0123 + 3*IT_0137 + 0.166666666666667
      *IT_0148) + 0.166666666666667*IT_0049*(IT_0149 + IT_0150 + IT_0151) +
       0.166666666666667*IT_0014*((-12)*IT_0140 + IT_0149 + IT_0150 + IT_0152 +
       IT_0153);
    const ccomplex_t IT_0155 = m_d*s_13;
    const ccomplex_t IT_0156 = 6*IT_0117;
    const ccomplex_t IT_0157 = 18*IT_0100;
    const ccomplex_t IT_0158 = (-18)*conj(IT_0100);
    const ccomplex_t IT_0159 = (-6)*conj(IT_0117);
    const ccomplex_t IT_0160 = IT_0156 + IT_0157 + IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = (-18)*IT_0117;
    const ccomplex_t IT_0162 = (-6)*IT_0100;
    const ccomplex_t IT_0163 = IT_0161 + IT_0162;
    const ccomplex_t IT_0164 = (-18)*conj(IT_0117);
    const ccomplex_t IT_0165 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0166 = IT_0164 + IT_0165;
    const ccomplex_t IT_0167 = (-6)*IT_0089*(conj(IT_0123) + 3*conj(IT_0137) +
       (-0.166666666666667)*IT_0148) + conj(IT_0089)*(IT_0151 + IT_0152 +
       IT_0153) + IT_0008*IT_0160 + conj(IT_0049)*IT_0163 + IT_0014*(IT_0147 +
       IT_0164 + IT_0165) + IT_0049*IT_0166;
    const ccomplex_t IT_0168 = m_d*m_u;
    const ccomplex_t IT_0169 = (-2)*IT_0008*IT_0013*(IT_0014 + 0.5*IT_0050) +
       IT_0053*(IT_0049*conj(IT_0049) + -IT_0014*(IT_0014 + IT_0050) + IT_0089
      *conj(IT_0089) + -IT_0090) + (conj(IT_0100)*IT_0117 + IT_0100*conj(IT_0117
      ) + conj(IT_0123)*IT_0137 + IT_0123*conj(IT_0137))*IT_0138 + (IT_0100*conj
      (IT_0100) + IT_0117*conj(IT_0117) + IT_0123*conj(IT_0123) + IT_0137*conj
      (IT_0137))*IT_0139 + (-6)*IT_0140*(s_34*IT_0137 + (-0.166666666666667)
      *conj(IT_0137)*IT_0138 + (-0.166666666666667)*conj(IT_0123)*IT_0139 + (
      -0.166666666666667)*(IT_0123 + IT_0140)*IT_0141) + (-6)*IT_0142*((-1)*
      (IT_0008 + 3*IT_0014 + (-3)*conj(IT_0049))*IT_0089 + IT_0008*conj(IT_0089)
       + (-0.166666666666667)*(IT_0014 + IT_0049)*IT_0143) + 6*IT_0144*IT_0154 +
       IT_0155*IT_0167 + (((-18)*IT_0100 + (-6)*IT_0117)*conj(IT_0123) + IT_0123
      *(IT_0158 + IT_0159) + IT_0140*IT_0160 + conj(IT_0137)*IT_0163 + IT_0137
      *IT_0166)*IT_0168;
    return create_ccomplex_return(IT_0169);
}

