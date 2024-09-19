#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cpow(V_tb, 2);
    const ccomplex_t IT_0009 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0010 = cos(alpha);
    const ccomplex_t IT_0011 = IT_0001*IT_0003;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = 2*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0005 + 1./3
      *IT_0004*IT_0016) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)*conj(U_sb_01)
       + 1./2*IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_11)) + U_sb_10*(IT_0002
      *IT_0005*IT_0009*IT_0010*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0004
      *IT_0014*IT_0016 + (-3)*IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)*conj
      (U_sb_11));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(e_em*IT_0018 +
       1.4142135623731*A_b*IT_0003*IT_0008*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_sb_1, 2);
    const ccomplex_t IT_0023 = pow(m_sb_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0003*IT_0012;
    const ccomplex_t IT_0031 = IT_0001*IT_0010;
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = (-2)*U_sb_00*((IT_0000*IT_0002*IT_0005*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0016)*IT_0032)*conj
      (U_sb_01) + -1./2*IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_11)) + U_sb_10
      *(IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_01) + (-2)*(IT_0000*IT_0002
      *IT_0005*IT_0010*IT_0017 + -1./3*m_W*IT_0004*IT_0016*IT_0032)*conj(U_sb_11
      ));
    const ccomplex_t IT_0034 = (0 + _Complex_I*-1)*(e_em*IT_0033 + (
      -1.4142135623731)*A_b*IT_0008*IT_0010*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0029*IT_0035;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0026 + -IT_0038;
    const ccomplex_t IT_0040 = -s_14;
    const ccomplex_t IT_0041 = s_13 + IT_0040;
    const ccomplex_t IT_0042 = m_d*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0044 = cpow(IT_0015, -1);
    const ccomplex_t IT_0045 = IT_0004*IT_0044;
    const ccomplex_t IT_0046 = e_em*IT_0045;
    const ccomplex_t IT_0047 = IT_0005*IT_0015;
    const ccomplex_t IT_0048 = e_em*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0046 + 3*IT_0048);
    const ccomplex_t IT_0050 = (-0.166666666666667)*IT_0049;
    const ccomplex_t IT_0051 = IT_0043*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0015 + 1./3*IT_0004*IT_0044)*U_sb_00*conj(U_sb_01) + IT_0004*IT_0044
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0054 = (-0.666666666666667)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0044;
    const ccomplex_t IT_0058 = 0.333333333333333*IT_0057;
    const ccomplex_t IT_0059 = IT_0043*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = (-0.333333333333333)*IT_0053;
    const ccomplex_t IT_0064 = m_d*IT_0063;
    const ccomplex_t IT_0065 = m_d*IT_0054;
    const ccomplex_t IT_0066 = (IT_0052 + -IT_0060)*(IT_0064 + -IT_0065);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = pow(m_d, 2);
    const ccomplex_t IT_0069 = (-18)*IT_0068;
    const ccomplex_t IT_0070 = 18*s_34;
    const ccomplex_t IT_0071 = IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = 36*s_34;
    const ccomplex_t IT_0073 = m_d*s_14;
    const ccomplex_t IT_0074 = (-18)*IT_0066;
    const ccomplex_t IT_0075 = m_d*e_em*IT_0000*IT_0002*IT_0005*IT_0012;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0078 = IT_0005*IT_0077;
    const ccomplex_t IT_0079 = m_b*e_em*mu_h*IT_0000*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0080 = IT_0005*IT_0079;
    const ccomplex_t IT_0081 = A_b*IT_0008*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0082 = IT_0012*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = A_b*IT_0008*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0085 = IT_0012*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0078 + 0.5*IT_0080 + (-0.5)*IT_0083 +
       0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0076*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = pow(m_Z, 2);
    const ccomplex_t IT_0093 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + IT_0092 
      + -reg_prop, -1);
    const ccomplex_t IT_0094 = m_d*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0002*IT_0005*IT_0012;
    const ccomplex_t IT_0097 = IT_0077*IT_0096;
    const ccomplex_t IT_0098 = IT_0079*IT_0096;
    const ccomplex_t IT_0099 = IT_0001*IT_0081;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0001*IT_0084;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0097 + 0.5*IT_0098 + 0.5*IT_0100 + (
      -0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0095*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = -IT_0091 + -IT_0106;
    const ccomplex_t IT_0108 = (-18)*s_34 + IT_0069;
    const ccomplex_t IT_0109 = (-18)*conj(IT_0066);
    const ccomplex_t IT_0110 = (-18)*IT_0073;
    const ccomplex_t IT_0111 = m_d*s_13;
    const ccomplex_t IT_0112 = 18*IT_0111;
    const ccomplex_t IT_0113 = 18*IT_0068 + IT_0070;
    const ccomplex_t IT_0114 = s_13 + s_14;
    const ccomplex_t IT_0115 = m_d*IT_0114;
    const ccomplex_t IT_0116 = (-18)*IT_0115;
    const ccomplex_t IT_0117 = 18*IT_0115;
    const ccomplex_t IT_0118 = s_34*IT_0022;
    const ccomplex_t IT_0119 = (-18)*IT_0118;
    const ccomplex_t IT_0120 = s_13*s_14;
    const ccomplex_t IT_0121 = 36*IT_0120;
    const ccomplex_t IT_0122 = IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = IT_0022*IT_0068;
    const ccomplex_t IT_0124 = (-18)*IT_0123;
    const ccomplex_t IT_0125 = 18*conj(IT_0039)*IT_0042;
    const ccomplex_t IT_0126 = IT_0039*(IT_0042*(18*conj(IT_0056) + 18*conj
      (IT_0062)) + (conj(IT_0066) + conj(IT_0067))*IT_0071 + conj(IT_0039)*((-36
      )*IT_0068 + IT_0072)) + (conj(IT_0067)*IT_0068 + conj(IT_0056)*IT_0073)
      *IT_0074 + IT_0067*(conj(IT_0067)*IT_0070 + conj(IT_0039)*IT_0071 + conj
      (IT_0107)*IT_0108 + IT_0068*IT_0109 + conj(IT_0062)*IT_0110 + conj(IT_0056
      )*IT_0112) + IT_0066*(conj(IT_0066)*IT_0070 + conj(IT_0039)*IT_0071 + conj
      (IT_0062)*IT_0112 + conj(IT_0107)*IT_0113) + IT_0107*((36*IT_0068 +
       IT_0072)*conj(IT_0107) + conj(IT_0067)*IT_0108 + conj(IT_0066)*IT_0113 +
       conj(IT_0056)*IT_0116 + conj(IT_0062)*IT_0117) + IT_0062*(conj(IT_0067)
      *IT_0110 + 18*conj(IT_0066)*IT_0111 + conj(IT_0107)*IT_0117 + conj(IT_0062
      )*IT_0122 + conj(IT_0056)*IT_0124 + IT_0125) + IT_0056*(IT_0073*IT_0109 +
       conj(IT_0067)*IT_0112 + conj(IT_0107)*IT_0116 + conj(IT_0056)*IT_0122 +
       conj(IT_0062)*IT_0124 + IT_0125);
    return create_ccomplex_return(IT_0126);
}

