#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cpow(V_tb, 2);
    const ccomplex_t IT_0009 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = -IT_0010*IT_0011;
    const ccomplex_t IT_0013 = IT_0001*IT_0003;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0014*(IT_0005 + 1.
      /3*IT_0004*IT_0016) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_11)) + U_sb_11
      *(IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0004*IT_0014*IT_0016 + (-3)*IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)
      *conj(U_sb_11));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0003*IT_0008*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_sb_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0010;
    const ccomplex_t IT_0029 = IT_0003*IT_0011;
    const ccomplex_t IT_0030 = IT_0028 + IT_0029;
    const ccomplex_t IT_0031 = 2*U_sb_01*((IT_0000*IT_0002*IT_0005*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0016)*IT_0030)*conj
      (U_sb_01) + 1./2*IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_11)) + U_sb_11*
      (IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0010*IT_0017 + -1./3*m_W*IT_0004*IT_0016*IT_0030)*conj(U_sb_11
      ));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_b*IT_0008*IT_0010*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = pow(m_e, 2);
    const ccomplex_t IT_0040 = cpow(IT_0015, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0015 + 1./3*IT_0004*IT_0040)*U_sb_01*conj(U_sb_01) + IT_0004*IT_0040
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0042 = (-0.666666666666667)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0040;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0042*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0004*IT_0040;
    const ccomplex_t IT_0050 = e_em*IT_0049;
    const ccomplex_t IT_0051 = IT_0005*IT_0015;
    const ccomplex_t IT_0052 = e_em*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0050 + -IT_0052);
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0044*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0042*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -s_14;
    const ccomplex_t IT_0060 = s_13 + IT_0059;
    const ccomplex_t IT_0061 = m_e*IT_0060;
    const ccomplex_t IT_0062 = (-0.333333333333333)*IT_0041;
    const ccomplex_t IT_0063 = m_e*IT_0062;
    const ccomplex_t IT_0064 = m_e*IT_0042;
    const ccomplex_t IT_0065 = (IT_0046 + -IT_0056)*(IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = -IT_0039;
    const ccomplex_t IT_0069 = s_34 + IT_0068;
    const ccomplex_t IT_0070 = IT_0022*IT_0039;
    const ccomplex_t IT_0071 = m_e*s_14;
    const ccomplex_t IT_0072 = (-6)*IT_0058;
    const ccomplex_t IT_0073 = m_e*e_em*IT_0000*IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0076 = IT_0005*IT_0075;
    const ccomplex_t IT_0077 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0078 = IT_0005*IT_0077;
    const ccomplex_t IT_0079 = A_b*IT_0008*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = A_b*IT_0008*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (-0.5)*IT_0076 + 0.5*IT_0078 + (-0.5)*IT_0081 +
       0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + (-2)*IT_0022 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_e*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0093 = IT_0075*IT_0092;
    const ccomplex_t IT_0094 = IT_0077*IT_0092;
    const ccomplex_t IT_0095 = IT_0001*IT_0079;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = IT_0001*IT_0082;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (-0.5)*IT_0093 + 0.5*IT_0094 + 0.5*IT_0096 + (
      -0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = pow(m_Z, 2);
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + (-2)*IT_0022 + IT_0101 + 
      -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0089 + -IT_0104;
    const ccomplex_t IT_0106 = s_34 + IT_0039;
    const ccomplex_t IT_0107 = (-6)*conj(IT_0067);
    const ccomplex_t IT_0108 = s_13 + s_14;
    const ccomplex_t IT_0109 = m_e*IT_0108;
    const ccomplex_t IT_0110 = -conj(IT_0058);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = 0.666666666666667*IT_0111;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = cpow(s_12 + IT_0022 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = s_13*s_14;
    const ccomplex_t IT_0120 = s_34*IT_0022;
    const ccomplex_t IT_0121 = (-6)*IT_0067;
    const ccomplex_t IT_0122 = (-6)*IT_0120;
    const ccomplex_t IT_0123 = (-6)*IT_0070;
    const ccomplex_t IT_0124 = 12*IT_0119;
    const ccomplex_t IT_0125 = IT_0122 + IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = IT_0118*((-12)*IT_0070 + 24*IT_0119 + (-12)
      *IT_0120) + -IT_0061*((-12)*IT_0038 + 12*conj(IT_0038) + (-6)*IT_0066 + 6
      *conj(IT_0066) + 6*conj(IT_0067) + IT_0121) + (IT_0048 + -conj(IT_0048) +
       IT_0058 + IT_0110)*IT_0125;
    const ccomplex_t IT_0127 = m_e*s_13;
    const ccomplex_t IT_0128 = IT_0122 + IT_0124;
    const ccomplex_t IT_0129 = conj(IT_0038)*IT_0061;
    const ccomplex_t IT_0130 = (-6)*conj(IT_0058);
    const ccomplex_t IT_0131 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0132 = 6*s_34;
    const ccomplex_t IT_0133 = conj(IT_0038)*IT_0069;
    const ccomplex_t IT_0134 = conj(IT_0105)*IT_0106;
    const ccomplex_t IT_0135 = IT_0038*(conj(IT_0038)*(12*s_34 + (-12)*IT_0039
      ) + 6*(conj(IT_0048) + conj(IT_0058))*IT_0061 + 6*(conj(IT_0066) + conj
      (IT_0067))*IT_0069) + (conj(IT_0048)*IT_0070 + conj(IT_0066)*IT_0071)
      *IT_0072 + IT_0105*((12*s_34 + 12*IT_0039)*conj(IT_0105) + 6*IT_0106*(conj
      (IT_0066) + 0.166666666666667*IT_0107) + 6*IT_0109*(conj(IT_0048) +
       IT_0110)) + -IT_0118*IT_0126 + (-6)*IT_0058*(conj(IT_0105)*IT_0109 + 
      -conj(IT_0067)*IT_0127 + (-0.166666666666667)*conj(IT_0058)*IT_0128 + 
      -IT_0129) + 6*IT_0048*(conj(IT_0105)*IT_0109 + conj(IT_0066)*IT_0127 +
       0.166666666666667*conj(IT_0048)*IT_0128 + IT_0129 + 0.166666666666667
      *IT_0070*IT_0130 + IT_0107*IT_0131) + 6*IT_0067*(conj(IT_0058)*IT_0127 +
       0.166666666666667*conj(IT_0067)*IT_0132 + IT_0133) + IT_0121*(IT_0039
      *conj(IT_0066) + conj(IT_0048)*IT_0071 + IT_0134) + 6*IT_0066*
      (0.166666666666667*IT_0039*IT_0107 + conj(IT_0048)*IT_0127 + IT_0130
      *IT_0131 + 0.166666666666667*conj(IT_0066)*IT_0132 + IT_0133 + IT_0134);
    return create_ccomplex_return(IT_0135);
}

