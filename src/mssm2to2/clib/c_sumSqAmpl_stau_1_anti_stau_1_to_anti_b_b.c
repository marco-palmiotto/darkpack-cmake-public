#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
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
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = mu_h*m_tau*IT_0000;
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = -IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0013*(IT_0005 + 
      -IT_0004*IT_0015) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_00) + -1./2*IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_10)) +
       U_stau_10*(IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_00) + 2*(m_W
      *IT_0004*IT_0013*IT_0015 + -IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_10));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0003*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_stau_1, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + (-2)*IT_0021 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0009;
    const ccomplex_t IT_0028 = IT_0003*IT_0010;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = 2*U_stau_00*((IT_0000*IT_0002*IT_0005*IT_0009
      *IT_0016 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0015)*IT_0029)*conj(U_stau_00
      ) + 1./2*IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_10)) + U_stau_10*
      (IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_00) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0009*IT_0016 + -m_W*IT_0004*IT_0015*IT_0029)*conj(U_stau_10));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0024 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_b, 2);
    const ccomplex_t IT_0039 = cpow(IT_0014, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0014 + -IT_0004*IT_0039)*U_stau_00*conj(U_stau_00) + IT_0004*IT_0039
      *U_stau_10*conj(U_stau_10));
    const ccomplex_t IT_0041 = (-2)*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*IT_0039;
    const ccomplex_t IT_0043 = e_em*IT_0042;
    const ccomplex_t IT_0044 = IT_0005*IT_0014;
    const ccomplex_t IT_0045 = e_em*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0043 + 3*IT_0045);
    const ccomplex_t IT_0047 = (-0.166666666666667)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0041*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0039;
    const ccomplex_t IT_0054 = 0.333333333333333*IT_0053;
    const ccomplex_t IT_0055 = IT_0048*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0041*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -s_14;
    const ccomplex_t IT_0060 = s_13 + IT_0059;
    const ccomplex_t IT_0061 = m_b*IT_0060;
    const ccomplex_t IT_0062 = -IT_0040;
    const ccomplex_t IT_0063 = m_b*IT_0062;
    const ccomplex_t IT_0064 = m_b*IT_0041;
    const ccomplex_t IT_0065 = (IT_0050 + -IT_0056)*(IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = -(IT_0050 + -IT_0056)*(IT_0063 + -IT_0064);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = -IT_0038;
    const ccomplex_t IT_0072 = s_34 + IT_0071;
    const ccomplex_t IT_0073 = IT_0021*IT_0038;
    const ccomplex_t IT_0074 = m_b*s_14;
    const ccomplex_t IT_0075 = (-6)*IT_0052;
    const ccomplex_t IT_0076 = m_b*e_em*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = e_em*mu_h*m_tau*IT_0000*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0079 = IT_0005*IT_0078;
    const ccomplex_t IT_0080 = e_em*mu_h*m_tau*IT_0000*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0081 = IT_0005*IT_0080;
    const ccomplex_t IT_0082 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0086 = IT_0010*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0079 + 0.5*IT_0081 + (-0.5)*IT_0084 +
       0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0077*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_12 + (-2)*IT_0021 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0096 = IT_0078*IT_0095;
    const ccomplex_t IT_0097 = IT_0080*IT_0095;
    const ccomplex_t IT_0098 = IT_0001*IT_0082;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = IT_0001*IT_0085;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0096 + 0.5*IT_0097 + 0.5*IT_0099 + (
      -0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0094*IT_0102;
    const ccomplex_t IT_0104 = pow(m_Z, 2);
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + (-2)*IT_0021 + IT_0104 + 
      -reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0092 + -IT_0107;
    const ccomplex_t IT_0109 = s_34 + IT_0038;
    const ccomplex_t IT_0110 = (-6)*conj(IT_0068);
    const ccomplex_t IT_0111 = s_13 + s_14;
    const ccomplex_t IT_0112 = m_b*IT_0111;
    const ccomplex_t IT_0113 = -conj(IT_0058);
    const ccomplex_t IT_0114 = cpow(s_12 + IT_0021 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0116 = (-0.333333333333333)*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0115;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = IT_0114*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = s_13*s_14;
    const ccomplex_t IT_0123 = s_34*IT_0021;
    const ccomplex_t IT_0124 = (-6)*IT_0068;
    const ccomplex_t IT_0125 = (-6)*IT_0123;
    const ccomplex_t IT_0126 = (-6)*IT_0073;
    const ccomplex_t IT_0127 = 12*IT_0122;
    const ccomplex_t IT_0128 = IT_0125 + IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = IT_0121*((-12)*IT_0073 + 24*IT_0122 + (-12)
      *IT_0123) + -IT_0061*((-12)*IT_0037 + 12*conj(IT_0037) + 6*conj(IT_0068) +
       (-6)*IT_0070 + 6*conj(IT_0070) + IT_0124) + (IT_0052 + -conj(IT_0052) +
       IT_0058 + IT_0113)*IT_0128;
    const ccomplex_t IT_0130 = m_b*s_13;
    const ccomplex_t IT_0131 = IT_0125 + IT_0127;
    const ccomplex_t IT_0132 = conj(IT_0037)*IT_0061;
    const ccomplex_t IT_0133 = (-6)*conj(IT_0052);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0135 = 6*s_34;
    const ccomplex_t IT_0136 = conj(IT_0037)*IT_0072;
    const ccomplex_t IT_0137 = conj(IT_0108)*IT_0109;
    const ccomplex_t IT_0138 = IT_0037*(conj(IT_0037)*(12*s_34 + (-12)*IT_0038
      ) + 6*(conj(IT_0052) + conj(IT_0058))*IT_0061 + 6*(conj(IT_0068) + conj
      (IT_0070))*IT_0072) + (conj(IT_0058)*IT_0073 + conj(IT_0070)*IT_0074)
      *IT_0075 + IT_0108*((12*s_34 + 12*IT_0038)*conj(IT_0108) + 6*IT_0109*(conj
      (IT_0070) + 0.166666666666667*IT_0110) + (-6)*IT_0112*(conj(IT_0052) +
       IT_0113)) + -IT_0121*IT_0129 + (-6)*IT_0052*(conj(IT_0108)*IT_0112 + 
      -conj(IT_0068)*IT_0130 + (-0.166666666666667)*conj(IT_0052)*IT_0131 + 
      -IT_0132) + 6*IT_0058*(conj(IT_0108)*IT_0112 + conj(IT_0070)*IT_0130 +
       0.166666666666667*conj(IT_0058)*IT_0131 + IT_0132 + 0.166666666666667
      *IT_0073*IT_0133 + IT_0110*IT_0134) + 6*IT_0068*(conj(IT_0052)*IT_0130 +
       0.166666666666667*conj(IT_0068)*IT_0135 + IT_0136) + IT_0124*(IT_0038
      *conj(IT_0070) + conj(IT_0058)*IT_0074 + IT_0137) + 6*IT_0070*
      (0.166666666666667*IT_0038*IT_0110 + conj(IT_0058)*IT_0130 + IT_0133
      *IT_0134 + 0.166666666666667*conj(IT_0070)*IT_0135 + IT_0136 + IT_0137);
    return create_ccomplex_return(IT_0138);
}

