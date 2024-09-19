#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_h_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(alpha);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = pow(m_W, -2);
    const ccomplex_t IT_0006 = pow(m_tau, 2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0002*((
      -0.5)*(IT_0004 + -(1 + -IT_0005*IT_0006*IT_0008)*IT_0010)*U_stau_00*conj
      (U_stau_01) + (IT_0004 + (-0.5)*IT_0005*IT_0006*IT_0008*IT_0010)*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = pow(m_stau_1, 2);
    const ccomplex_t IT_0014 = pow(m_stau_2, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + -IT_0013 + -IT_0014 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = mu_h*m_tau*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0007, -1);
    const ccomplex_t IT_0019 = cpow(IT_0009, -1);
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = -IT_0002*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0007;
    const ccomplex_t IT_0023 = IT_0021 + IT_0022;
    const ccomplex_t IT_0024 = (-2)*U_stau_00*((IT_0001*IT_0006*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0023)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0006*IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0023)
      *conj(U_stau_11));
    const ccomplex_t IT_0025 = (0 + _Complex_I*-1)*(e_em*IT_0024 + (
      -1.4142135623731)*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0004*IT_0009;
    const ccomplex_t IT_0028 = IT_0019 + IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0002, 3);
    const ccomplex_t IT_0030 = cpow(IT_0002, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_W*e_em*IT_0028*(IT_0007
      *IT_0029 + 5*IT_0001*((-1)*IT_0001*(IT_0002*IT_0007 + 1./5*IT_0001*IT_0020
      ) + IT_0020*IT_0030));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0032;
    const ccomplex_t IT_0034 = IT_0015*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0013 + -IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0037 = cpow(IT_0001, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_W*e_em*IT_0028*(IT_0020
      *IT_0029 + (-5)*IT_0001*(IT_0001*IT_0002*IT_0020 + 1*IT_0007*(IT_0030 + 
      -1./5*IT_0037)));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0002*IT_0007;
    const ccomplex_t IT_0041 = IT_0001*IT_0020;
    const ccomplex_t IT_0042 = IT_0040 + IT_0041;
    const ccomplex_t IT_0043 = 2*U_stau_00*((IT_0002*IT_0006*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0042)*conj(U_stau_01) +
       1./2*IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) + U_stau_10*
      (IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + 2*(IT_0002*IT_0006
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0042)*conj(U_stau_11));
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(e_em*IT_0043 +
       1.4142135623731*A_tau*IT_0002*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0039*IT_0045;
    const ccomplex_t IT_0047 = IT_0036*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = pow(m_h, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0049 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0051 = 2*U_stau_00*((IT_0002*IT_0006*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0042)*conj(U_stau_00) +
       1./2*IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_10)) + U_stau_10*
      (IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_00) + 2*(IT_0002*IT_0006
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0042)*conj(U_stau_10));
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(e_em*IT_0051 +
       1.4142135623731*A_tau*IT_0002*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = (-2)*U_stau_00*((IT_0001*IT_0006*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0023)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0006*IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0023)
      *conj(U_stau_11));
    const ccomplex_t IT_0055 = (0 + _Complex_I*-1)*(e_em*IT_0054 + (
      -1.4142135623731)*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = IT_0050*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0049 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0061 = (-2)*U_stau_00*((IT_0001*IT_0006*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0023)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0006*IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0023)
      *conj(U_stau_11));
    const ccomplex_t IT_0062 = (0 + _Complex_I*-1)*(e_em*IT_0061 + (
      -1.4142135623731)*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = 2*U_stau_01*((IT_0002*IT_0006*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0042)*conj(U_stau_01) +
       1./2*IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) + U_stau_11*
      (IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + 2*(IT_0002*IT_0006
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0042)*conj(U_stau_11));
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(e_em*IT_0064 +
       1.4142135623731*A_tau*IT_0002*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0063*IT_0066;
    const ccomplex_t IT_0068 = IT_0060*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0014 + IT_0049 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0071 = (-2)*U_stau_00*((IT_0001*IT_0006*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0023)*conj
      (U_stau_00) + -1./2*IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_10)) +
       U_stau_10*(IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_00) + (-2)*
      (IT_0001*IT_0006*IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0023)
      *conj(U_stau_10));
    const ccomplex_t IT_0072 = (0 + _Complex_I*-1)*(e_em*IT_0071 + (
      -1.4142135623731)*A_tau*IT_0001*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = 2*U_stau_00*((IT_0002*IT_0006*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0042)*conj(U_stau_01) +
       1./2*IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) + U_stau_10*
      (IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + 2*(IT_0002*IT_0006
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0042)*conj(U_stau_11));
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(e_em*IT_0074 +
       1.4142135623731*A_tau*IT_0002*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0073*IT_0076;
    const ccomplex_t IT_0078 = IT_0070*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_13 + IT_0013 + IT_0049 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0081 = (-2)*U_stau_01*((IT_0001*IT_0006*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0023)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) +
       U_stau_11*(IT_0002*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0006*IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0023)
      *conj(U_stau_11));
    const ccomplex_t IT_0082 = (0 + _Complex_I*-1)*(e_em*IT_0081 + (
      -1.4142135623731)*A_tau*IT_0001*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = 2*U_stau_00*((IT_0002*IT_0006*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0019)*IT_0042)*conj(U_stau_01) +
       1./2*IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) + U_stau_10*
      (IT_0001*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + 2*(IT_0002*IT_0006
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0004*IT_0009*IT_0042)*conj(U_stau_11));
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(e_em*IT_0084 +
       1.4142135623731*A_tau*IT_0002*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0083*IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0012 + -IT_0035 + -IT_0048 + IT_0059 +
       IT_0069 + IT_0079 + IT_0089;
    return create_ccomplex_return(IT_0090*conj(IT_0090));
}

