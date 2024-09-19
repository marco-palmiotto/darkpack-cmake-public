#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_h_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_st_2, 2);
    const ccomplex_t IT_0001 = s_12 + IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = e_em*IT_0003*(IT_0004*IT_0005 + IT_0006*IT_0007
      )*IT_0009;
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0003*IT_0008 
      + (-3)*IT_0002*IT_0009)*U_st_01*conj(U_st_01) + IT_0003*IT_0008*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0013 = 0.666666666666667*IT_0012;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = s_13 + s_23;
    const ccomplex_t IT_0017 = -IT_0010;
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = 1.33333333333333*IT_0012;
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0020;
    const ccomplex_t IT_0024 = s_23*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0019 + (-0.5)*IT_0022 + 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0015 + IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = pow(m_h, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0030 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0032 = pow(m_W, -1);
    const ccomplex_t IT_0033 = m_t*e_em*mu_h*IT_0032*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = m_t*e_em*mu_h*IT_0032*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0036 = IT_0009*IT_0035;
    const ccomplex_t IT_0037 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0038 = IT_0004*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0041 = IT_0004*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0034 + 0.5*IT_0036 + (-0.5)*IT_0039 +
       0.5*IT_0042;
    const ccomplex_t IT_0044 = m_t*mu_h*IT_0032;
    const ccomplex_t IT_0045 = cpow(IT_0006, -1);
    const ccomplex_t IT_0046 = IT_0004*IT_0007;
    const ccomplex_t IT_0047 = IT_0005*IT_0006;
    const ccomplex_t IT_0048 = IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = cpow(IT_0002, -2);
    const ccomplex_t IT_0050 = pow(m_t, 2);
    const ccomplex_t IT_0051 = 2*U_st_01*(((-0.5)*m_W*IT_0048*(IT_0009 + -1./3
      *IT_0008*IT_0049) + IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*conj(U_st_01)
       + 1./2*IT_0007*IT_0009*IT_0044*IT_0045*conj(U_st_11)) + U_st_11*(IT_0007
      *IT_0009*IT_0044*IT_0045*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0008
      *IT_0048*IT_0049 + (-1.5)*IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*conj
      (U_st_11));
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(e_em*IT_0051 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0043*IT_0053;
    const ccomplex_t IT_0055 = IT_0031*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_13 + IT_0030 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0054*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0008*IT_0049;
    const ccomplex_t IT_0061 = IT_0009 + IT_0060;
    const ccomplex_t IT_0062 = cpow(IT_0006, 3);
    const ccomplex_t IT_0063 = cpow(IT_0006, 2);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_W*e_em*IT_0061*(IT_0005
      *IT_0062 + -IT_0004*(IT_0004*(IT_0005*IT_0006 + IT_0004*IT_0007) + 
      -IT_0007*IT_0063));
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0043*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = m_t*e_em*mu_h*IT_0032*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = m_t*e_em*mu_h*IT_0032*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0073 = IT_0009*IT_0072;
    const ccomplex_t IT_0074 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0075 = IT_0004*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0078 = IT_0004*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0071 + 0.5*IT_0073 + (-0.5)*IT_0076 +
       0.5*IT_0079;
    const ccomplex_t IT_0081 = 2*conj(U_st_00)*(((-0.5)*m_W*IT_0048*(IT_0009 +
       -1./3*IT_0008*IT_0049) + IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*U_st_01
       + 1./2*IT_0007*IT_0009*IT_0044*IT_0045*U_st_11) + conj(U_st_10)*(IT_0007
      *IT_0009*IT_0044*IT_0045*U_st_01 + (-1.33333333333333)*(m_W*IT_0008
      *IT_0048*IT_0049 + (-1.5)*IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*U_st_11);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(e_em*IT_0081 +
       1.4142135623731*A_t*IT_0005*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0030 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = pow(m_Z, 2);
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0088 + 
      -reg_prop, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*IT_0006*
      (IT_0005*IT_0006 + IT_0004*IT_0007)*IT_0061;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0004*IT_0009*IT_0045;
    const ccomplex_t IT_0093 = IT_0033*IT_0092;
    const ccomplex_t IT_0094 = IT_0035*IT_0092;
    const ccomplex_t IT_0095 = IT_0006*IT_0037;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = IT_0006*IT_0040;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0093 + (-0.5)*IT_0094 + (-0.5)*IT_0096 +
       0.5*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = IT_0089*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0104 = m_t*e_em*mu_h*IT_0032*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0105 = IT_0009*IT_0104;
    const ccomplex_t IT_0106 = m_t*e_em*mu_h*IT_0032*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0107 = IT_0009*IT_0106;
    const ccomplex_t IT_0108 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0109 = IT_0004*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0112 = IT_0004*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = 0.5*IT_0105 + (-0.5)*IT_0107 + 0.5*IT_0110 + (
      -0.5)*IT_0113;
    const ccomplex_t IT_0115 = 2*U_st_00*(((-0.5)*m_W*IT_0048*(IT_0009 + -1./3
      *IT_0008*IT_0049) + IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*conj(U_st_01)
       + 1./2*IT_0007*IT_0009*IT_0044*IT_0045*conj(U_st_11)) + U_st_10*(IT_0007
      *IT_0009*IT_0044*IT_0045*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0008
      *IT_0048*IT_0049 + (-1.5)*IT_0005*IT_0009*IT_0032*IT_0045*IT_0050)*conj
      (U_st_11));
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(e_em*IT_0115 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = IT_0114*IT_0117;
    const ccomplex_t IT_0119 = IT_0103*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (-2)*IT_0029 + IT_0056 + IT_0059 + -IT_0069 +
       IT_0087 + -IT_0102 + IT_0120;
    return create_ccomplex_return(3*IT_0121*conj(IT_0121));
}

