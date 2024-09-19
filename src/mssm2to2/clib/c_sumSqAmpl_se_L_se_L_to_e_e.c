#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_se_L_to_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_se_L_to_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_N_1*IT_0007;
    const ccomplex_t IT_0009 = pow(m_e, 2);
    const ccomplex_t IT_0010 = pow(m_se_L, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*m_N_2*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0020, 2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_N_3*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0028, 2);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_N_4*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = cos(theta_W);
    const ccomplex_t IT_0034 = cpow(IT_0033, -1);
    const ccomplex_t IT_0035 = conj(N_B1)*e_em;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W1)*e_em;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0006*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_e*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*e_em;
    const ccomplex_t IT_0049 = IT_0034*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W2)*e_em;
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0014*IT_0055;
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = m_e*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0056*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = m_e*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B3)*e_em;
    const ccomplex_t IT_0065 = IT_0034*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W3)*e_em;
    const ccomplex_t IT_0068 = IT_0004*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + IT_0069);
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0020*IT_0071;
    const ccomplex_t IT_0073 = IT_0023*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = m_e*IT_0074;
    const ccomplex_t IT_0076 = IT_0025*IT_0072;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_e*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B4)*e_em;
    const ccomplex_t IT_0080 = IT_0034*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W4)*e_em;
    const ccomplex_t IT_0083 = IT_0004*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + IT_0084);
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0028*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = m_e*IT_0090;
    const ccomplex_t IT_0092 = IT_0031*IT_0087;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = m_e*IT_0093;
    const ccomplex_t IT_0095 = IT_0011*IT_0043;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = m_e*IT_0096;
    const ccomplex_t IT_0098 = IT_0012 + IT_0018 + IT_0024 + IT_0026 + IT_0032
       + -IT_0047 + -IT_0059 + -IT_0063 + -IT_0075 + -IT_0078 + -IT_0091 + 
      -IT_0094 + -IT_0097;
    const ccomplex_t IT_0099 = IT_0008*IT_0044;
    const ccomplex_t IT_0100 = IT_0016*IT_0060;
    const ccomplex_t IT_0101 = IT_0030*IT_0088;
    const ccomplex_t IT_0102 = -IT_0099 + -IT_0100 + -IT_0101;
    const ccomplex_t IT_0103 = cpow(IT_0042, 2);
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_N_1*IT_0103;
    const ccomplex_t IT_0105 = IT_0044*IT_0104;
    const ccomplex_t IT_0106 = cpow(IT_0055, 2);
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_N_2*IT_0106;
    const ccomplex_t IT_0108 = IT_0017*IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0071, 2);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_N_3*IT_0109;
    const ccomplex_t IT_0111 = IT_0023*IT_0110;
    const ccomplex_t IT_0112 = IT_0025*IT_0110;
    const ccomplex_t IT_0113 = cpow(IT_0086, 2);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_N_4*IT_0113;
    const ccomplex_t IT_0115 = IT_0031*IT_0114;
    const ccomplex_t IT_0116 = -IT_0047 + -IT_0059 + -IT_0063 + -IT_0075 + 
      -IT_0078 + -IT_0091 + -IT_0094 + -IT_0097 + IT_0105 + IT_0108 + IT_0111 +
       IT_0112 + IT_0115;
    const ccomplex_t IT_0117 = IT_0011*IT_0104;
    const ccomplex_t IT_0118 = IT_0060*IT_0107;
    const ccomplex_t IT_0119 = IT_0088*IT_0114;
    const ccomplex_t IT_0120 = -IT_0117 + -IT_0118 + -IT_0119;
    const ccomplex_t IT_0121 = 2*s_34;
    const ccomplex_t IT_0122 = (-2)*s_34;
    const ccomplex_t IT_0123 = -IT_0046 + -IT_0058 + IT_0062 + -IT_0074 +
       IT_0077 + -IT_0090 + IT_0093 + IT_0096;
    const ccomplex_t IT_0124 = m_e*s_23;
    const ccomplex_t IT_0125 = IT_0116 + -IT_0120;
    const ccomplex_t IT_0126 = 2*IT_0102;
    const ccomplex_t IT_0127 = (-2)*IT_0098;
    const ccomplex_t IT_0128 = m_e*s_24;
    const ccomplex_t IT_0129 = 2*IT_0098;
    const ccomplex_t IT_0130 = (-2)*IT_0102;
    const ccomplex_t IT_0131 = 2*conj(IT_0102);
    const ccomplex_t IT_0132 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0133 = 2*conj(IT_0098);
    const ccomplex_t IT_0134 = (-2)*conj(IT_0102);
    const ccomplex_t IT_0135 = s_23*s_24;
    const ccomplex_t IT_0136 = s_34*IT_0010;
    const ccomplex_t IT_0137 = IT_0009*IT_0010;
    const ccomplex_t IT_0138 = -conj(IT_0120);
    const ccomplex_t IT_0139 = conj(IT_0116) + IT_0138;
    const ccomplex_t IT_0140 = 8*IT_0123*(conj(IT_0123)*(IT_0135 + (-0.5)
      *IT_0136 + (-0.5)*IT_0137) + (-0.25)*IT_0124*((-0.5)*IT_0131 + (-0.5)
      *IT_0132 + IT_0139) + 0.25*IT_0128*(0.5*IT_0133 + 0.5*IT_0134 + IT_0139));
    const ccomplex_t IT_0141 = (IT_0098*conj(IT_0098) + IT_0102*conj(IT_0102) 
      + IT_0116*conj(IT_0116) + IT_0120*conj(IT_0120))*IT_0121 + (conj(IT_0098)
      *IT_0102 + IT_0098*conj(IT_0102) + conj(IT_0116)*IT_0120 + IT_0116*conj
      (IT_0120))*IT_0122 + (-2)*conj(IT_0123)*(IT_0124*(IT_0125 + (-0.5)*IT_0126
       + (-0.5)*IT_0127) + -IT_0128*(IT_0125 + 0.5*IT_0129 + 0.5*IT_0130)) +
       IT_0009*(conj(IT_0116)*(IT_0126 + IT_0127) + conj(IT_0120)*(IT_0129 +
       IT_0130) + IT_0116*(IT_0131 + IT_0132) + IT_0120*(IT_0133 + IT_0134)) +
       IT_0140;
    return create_ccomplex_return(IT_0141);
}

