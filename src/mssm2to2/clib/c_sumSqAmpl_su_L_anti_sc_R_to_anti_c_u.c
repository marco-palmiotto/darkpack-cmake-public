#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sc_R_to_anti_c_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sc_R_to_anti_c_u(
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
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sc_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_sG*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_B1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_W1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0014 + 3*IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0011;
    const ccomplex_t IT_0023 = 0.666666666666667*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0009*IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0028 = conj(N_B2)*e_em;
    const ccomplex_t IT_0029 = IT_0011*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = conj(N_W2)*e_em;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + 3*IT_0033);
    const ccomplex_t IT_0035 = 0.166666666666667*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0011;
    const ccomplex_t IT_0037 = 0.666666666666667*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0027*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0042 = conj(N_B3)*e_em;
    const ccomplex_t IT_0043 = IT_0011*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W3)*e_em;
    const ccomplex_t IT_0046 = IT_0016*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + 3*IT_0047);
    const ccomplex_t IT_0049 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0011;
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0041*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0056 = conj(N_B4)*e_em;
    const ccomplex_t IT_0057 = IT_0011*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W4)*e_em;
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 3*IT_0061);
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0011;
    const ccomplex_t IT_0065 = 0.666666666666667*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0055*IT_0066*IT_0067;
    const ccomplex_t IT_0069 = pow(m_W, -1);
    const ccomplex_t IT_0070 = sin(beta);
    const ccomplex_t IT_0071 = cpow(IT_0070, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0023*IT_0073;
    const ccomplex_t IT_0075 = IT_0025*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = m_u*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0037*IT_0079;
    const ccomplex_t IT_0081 = IT_0039*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = m_u*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0051*IT_0085;
    const ccomplex_t IT_0087 = IT_0053*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = m_u*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0065*IT_0091;
    const ccomplex_t IT_0093 = IT_0067*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = m_u*IT_0094;
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0007 + -IT_0026 + -IT_0040
       + -IT_0054 + -IT_0068 + IT_0077 + IT_0083 + IT_0089 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0073*IT_0098;
    const ccomplex_t IT_0100 = IT_0009*IT_0025*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0079*IT_0102;
    const ccomplex_t IT_0104 = IT_0027*IT_0039*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0085*IT_0106;
    const ccomplex_t IT_0108 = IT_0041*IT_0053*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0016*IT_0069*IT_0071;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0091*IT_0110;
    const ccomplex_t IT_0112 = IT_0055*IT_0067*IT_0111;
    const ccomplex_t IT_0113 = IT_0021*IT_0098;
    const ccomplex_t IT_0114 = IT_0025*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = m_u*IT_0115;
    const ccomplex_t IT_0117 = IT_0035*IT_0102;
    const ccomplex_t IT_0118 = IT_0039*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = m_u*IT_0119;
    const ccomplex_t IT_0121 = IT_0049*IT_0106;
    const ccomplex_t IT_0122 = IT_0053*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_u*IT_0123;
    const ccomplex_t IT_0125 = IT_0063*IT_0110;
    const ccomplex_t IT_0126 = IT_0067*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = m_u*IT_0127;
    const ccomplex_t IT_0129 = -IT_0100 + -IT_0104 + -IT_0108 + -IT_0112 +
       IT_0116 + IT_0120 + IT_0124 + IT_0128;
    const ccomplex_t IT_0130 = m_c*m_u;
    const ccomplex_t IT_0131 = -IT_0115 + -IT_0119 + -IT_0123 + -IT_0127;
    const ccomplex_t IT_0132 = -IT_0076 + -IT_0082 + -IT_0088 + -IT_0094;
    const ccomplex_t IT_0133 = pow(m_su_L, 2);
    const ccomplex_t IT_0134 = m_c*m_u*IT_0133;
    const ccomplex_t IT_0135 = s_34*IT_0133;
    const ccomplex_t IT_0136 = (-18)*IT_0135;
    const ccomplex_t IT_0137 = s_13*s_14;
    const ccomplex_t IT_0138 = 36*IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = 18*s_34;
    const ccomplex_t IT_0141 = (-18)*IT_0129;
    const ccomplex_t IT_0142 = (-18)*conj(IT_0129);
    const ccomplex_t IT_0143 = m_c*s_14;
    const ccomplex_t IT_0144 = 6*IT_0008;
    const ccomplex_t IT_0145 = (-6)*IT_0008;
    const ccomplex_t IT_0146 = m_u*s_13;
    const ccomplex_t IT_0147 = (-18)*IT_0008*(s_34*(IT_0008 + 1./3*IT_0096 + 
      -1./3*conj(IT_0096)) + -1./3*(IT_0129 + -conj(IT_0129))*IT_0130) + conj
      (IT_0131)*((-18)*IT_0132*IT_0134 + IT_0131*IT_0139) + conj(IT_0132)*((-18)
      *IT_0131*IT_0134 + IT_0132*IT_0139) + (IT_0096*conj(IT_0096) + IT_0129
      *conj(IT_0129))*IT_0140 + IT_0130*(conj(IT_0096)*IT_0141 + IT_0096*IT_0142
      ) + IT_0143*(conj(IT_0132)*IT_0141 + IT_0132*IT_0142 + IT_0131*((-18)*conj
      (IT_0096) + IT_0144) + conj(IT_0131)*((-18)*IT_0096 + IT_0145)) + (18*conj
      (IT_0129)*IT_0131 + 18*IT_0129*conj(IT_0131) + conj(IT_0132)*(18*IT_0096 +
       IT_0144) + IT_0132*(18*conj(IT_0096) + IT_0145))*IT_0146;
    return create_ccomplex_return(IT_0147);
}

