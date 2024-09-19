#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_sc_R_to_c_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_sc_R_to_c_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_R = param->m_sc_R;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sc_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0004;
    const ccomplex_t IT_0006 = 0.666666666666667*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = IT_0002*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0004;
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0004;
    const ccomplex_t IT_0029 = 0.666666666666667*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0004;
    const ccomplex_t IT_0036 = 0.666666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = IT_0002*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0036*IT_0048;
    const ccomplex_t IT_0050 = IT_0038*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0021*IT_0053;
    const ccomplex_t IT_0055 = IT_0017*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0029*IT_0058;
    const ccomplex_t IT_0060 = IT_0025*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0046 + -IT_0051 + -IT_0056 + -IT_0061;
    const ccomplex_t IT_0063 = pow(m_su_R, 2);
    const ccomplex_t IT_0064 = m_c*m_u*IT_0063;
    const ccomplex_t IT_0065 = s_34*IT_0063;
    const ccomplex_t IT_0066 = (-18)*IT_0065;
    const ccomplex_t IT_0067 = s_13*s_14;
    const ccomplex_t IT_0068 = 36*IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = m_u*s_13;
    const ccomplex_t IT_0071 = m_u*IT_0016;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0073 = IT_0027*IT_0058;
    const ccomplex_t IT_0074 = IT_0025*IT_0072*IT_0073;
    const ccomplex_t IT_0075 = IT_0034*IT_0048;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0077 = IT_0038*IT_0075*IT_0076;
    const ccomplex_t IT_0078 = m_u*IT_0032;
    const ccomplex_t IT_0079 = m_u*IT_0024;
    const ccomplex_t IT_0080 = m_u*IT_0040;
    const ccomplex_t IT_0081 = -IT_0071 + IT_0074 + IT_0077 + -IT_0078 + 
      -IT_0079 + -IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0083 = IT_0013*IT_0043;
    const ccomplex_t IT_0084 = IT_0002*IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0086 = IT_0019*IT_0053;
    const ccomplex_t IT_0087 = IT_0017*IT_0085*IT_0086;
    const ccomplex_t IT_0088 = -IT_0084 + -IT_0087;
    const ccomplex_t IT_0089 = pow(m_sG, 2);
    const ccomplex_t IT_0090 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -IT_0089 +
       reg_prop, -1);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0092 = cpow(IT_0091, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_sG*IT_0092;
    const ccomplex_t IT_0094 = IT_0090*IT_0093;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = 6*IT_0095;
    const ccomplex_t IT_0097 = cpow(IT_0029, 2);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_1*IT_0097;
    const ccomplex_t IT_0099 = IT_0025*IT_0098;
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0094 + -IT_0099;
    const ccomplex_t IT_0101 = 18*IT_0100;
    const ccomplex_t IT_0102 = cpow(IT_0021, 2);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_N_3*IT_0102;
    const ccomplex_t IT_0104 = IT_0017*IT_0103;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = 18*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0036, 2);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_4*IT_0107;
    const ccomplex_t IT_0109 = IT_0038*IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = 18*IT_0110;
    const ccomplex_t IT_0112 = m_u*IT_0061;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = (-18)*IT_0113;
    const ccomplex_t IT_0115 = cpow(IT_0006, 2);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*m_N_2*IT_0115;
    const ccomplex_t IT_0117 = IT_0002*IT_0116;
    const ccomplex_t IT_0118 = m_u*IT_0046;
    const ccomplex_t IT_0119 = m_u*IT_0051;
    const ccomplex_t IT_0120 = IT_0117 + -IT_0118 + -IT_0119;
    const ccomplex_t IT_0121 = (-18)*IT_0120;
    const ccomplex_t IT_0122 = m_u*IT_0056;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = (-18)*IT_0123;
    const ccomplex_t IT_0125 = IT_0096 + IT_0101 + IT_0106 + IT_0111 + IT_0114
       + IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = 18*conj(IT_0100);
    const ccomplex_t IT_0127 = 18*conj(IT_0105);
    const ccomplex_t IT_0128 = 18*conj(IT_0110);
    const ccomplex_t IT_0129 = (-18)*conj(IT_0113);
    const ccomplex_t IT_0130 = (-18)*conj(IT_0120);
    const ccomplex_t IT_0131 = (-18)*conj(IT_0123);
    const ccomplex_t IT_0132 = (-6)*IT_0095;
    const ccomplex_t IT_0133 = IT_0126 + IT_0127 + IT_0128 + IT_0129 + IT_0130
       + IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = conj(IT_0120) + conj(IT_0123);
    const ccomplex_t IT_0135 = conj(IT_0105) + conj(IT_0110);
    const ccomplex_t IT_0136 = 18*s_34;
    const ccomplex_t IT_0137 = IT_0105 + IT_0110;
    const ccomplex_t IT_0138 = (-18)*s_34;
    const ccomplex_t IT_0139 = IT_0100*IT_0136;
    const ccomplex_t IT_0140 = IT_0113*IT_0138;
    const ccomplex_t IT_0141 = conj(IT_0113)*IT_0138;
    const ccomplex_t IT_0142 = IT_0120*IT_0136;
    const ccomplex_t IT_0143 = IT_0100*IT_0138;
    const ccomplex_t IT_0144 = conj(IT_0113)*IT_0136;
    const ccomplex_t IT_0145 = conj(IT_0100)*IT_0138;
    const ccomplex_t IT_0146 = IT_0105*IT_0136;
    const ccomplex_t IT_0147 = IT_0120*IT_0138;
    const ccomplex_t IT_0148 = conj(IT_0100)*IT_0136;
    const ccomplex_t IT_0149 = conj(IT_0120)*IT_0138;
    const ccomplex_t IT_0150 = m_c*m_u;
    const ccomplex_t IT_0151 = 18*IT_0113;
    const ccomplex_t IT_0152 = 18*IT_0120;
    const ccomplex_t IT_0153 = 18*IT_0123;
    const ccomplex_t IT_0154 = (-18)*IT_0100;
    const ccomplex_t IT_0155 = (-18)*IT_0105;
    const ccomplex_t IT_0156 = (-18)*IT_0110;
    const ccomplex_t IT_0157 = IT_0132 + IT_0151 + IT_0152 + IT_0153 + IT_0154
       + IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = 18*conj(IT_0113);
    const ccomplex_t IT_0159 = 18*conj(IT_0120);
    const ccomplex_t IT_0160 = 18*conj(IT_0123);
    const ccomplex_t IT_0161 = (-18)*conj(IT_0100);
    const ccomplex_t IT_0162 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0163 = (-18)*conj(IT_0110);
    const ccomplex_t IT_0164 = m_c*s_14;
    const ccomplex_t IT_0165 = (conj(IT_0062)*(18*IT_0081 + (-18)*IT_0088) +
       IT_0062*(18*conj(IT_0081) + (-18)*conj(IT_0088)) + conj(IT_0041)*IT_0157 
      + IT_0041*(IT_0096 + IT_0158 + IT_0159 + IT_0160 + IT_0161 + IT_0162 +
       IT_0163))*IT_0164;
    const ccomplex_t IT_0166 = 6*s_34;
    const ccomplex_t IT_0167 = (-6)*s_34;
    const ccomplex_t IT_0168 = conj(IT_0041)*((-18)*IT_0062*IT_0064 + IT_0041
      *IT_0069) + conj(IT_0062)*((-18)*IT_0041*IT_0064 + IT_0062*IT_0069) +
       IT_0070*(conj(IT_0041)*((-18)*IT_0081 + 18*IT_0088) + IT_0041*((-18)*conj
      (IT_0081) + 18*conj(IT_0088)) + conj(IT_0062)*IT_0125 + IT_0062*IT_0133) +
       (IT_0081*conj(IT_0081) + IT_0088*conj(IT_0088) + IT_0113*(conj(IT_0113) +
       conj(IT_0120) + conj(IT_0123)) + IT_0123*IT_0134 + IT_0110*IT_0135)
      *IT_0136 + (conj(IT_0081)*IT_0088 + IT_0081*conj(IT_0088) + IT_0123
      *IT_0135 + conj(IT_0123)*IT_0137)*IT_0138 + conj(IT_0100)*(IT_0139 +
       IT_0140) + IT_0100*IT_0141 + IT_0134*(IT_0142 + IT_0143) + (IT_0120 +
       IT_0123)*(IT_0144 + IT_0145) + IT_0135*(IT_0139 + IT_0140 + IT_0146 +
       IT_0147) + IT_0137*(IT_0141 + IT_0148 + IT_0149) + IT_0150*(conj(IT_0081)
      *IT_0125 + IT_0081*IT_0133 + conj(IT_0088)*IT_0157 + IT_0088*(IT_0096 +
       IT_0158 + IT_0159 + IT_0160 + IT_0161 + IT_0162 + IT_0163)) + IT_0165 +
       IT_0095*(IT_0095*IT_0138 + (conj(IT_0100) + conj(IT_0105) + conj(IT_0110)
       + IT_0113 + IT_0120 + IT_0123)*IT_0166 + (IT_0100 + IT_0105 + IT_0110 +
       conj(IT_0113) + conj(IT_0120) + conj(IT_0123))*IT_0167);
    return create_ccomplex_return(IT_0168);
}

