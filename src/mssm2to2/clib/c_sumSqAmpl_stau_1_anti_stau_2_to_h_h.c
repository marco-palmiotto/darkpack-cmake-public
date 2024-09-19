#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_h_h(
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
    const ccomplex_t IT_0001 = sin(alpha);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -2);
    const ccomplex_t IT_0008 = pow(m_W, -2);
    const ccomplex_t IT_0009 = pow(m_tau, 2);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = (-0.5)*(IT_0005*(IT_0007 + -IT_0013) + IT_0003*
      (IT_0007 + -(1 + (-2)*IT_0008*IT_0009*IT_0011)*IT_0013))*U_stau_00*conj
      (U_stau_01) + ((IT_0003 + IT_0005)*IT_0007 + IT_0002*IT_0008*IT_0009
      *IT_0011*IT_0013)*U_stau_10*conj(U_stau_11);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0000*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_stau_1, 2);
    const ccomplex_t IT_0018 = pow(m_stau_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + -IT_0017 + -IT_0018 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0007*IT_0012;
    const ccomplex_t IT_0021 = cpow(IT_0012, -1);
    const ccomplex_t IT_0022 = IT_0020 + IT_0021;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0001, 3);
    const ccomplex_t IT_0025 = (0 + _Complex_I*-5)*m_W*e_em*IT_0022*(IT_0004*
      (1*(IT_0002 + -1./5*IT_0005)*IT_0010 + IT_0001*IT_0004*IT_0023) + -1./5
      *IT_0023*IT_0024);
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = mu_h*m_tau*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0010, -1);
    const ccomplex_t IT_0030 = -IT_0001*IT_0023;
    const ccomplex_t IT_0031 = IT_0004*IT_0010;
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = (-2)*U_stau_00*((IT_0004*IT_0009*IT_0021
      *IT_0027*IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0032)*conj
      (U_stau_01) + -1./2*IT_0001*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) +
       U_stau_10*(IT_0001*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + (-2)*
      (IT_0004*IT_0009*IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0032)
      *conj(U_stau_11));
    const ccomplex_t IT_0034 = (0 + _Complex_I*-1)*(e_em*IT_0033 + (
      -1.4142135623731)*A_tau*IT_0004*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = IT_0019*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0017 + -IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0004*
      (IT_0002*IT_0023 + -IT_0004*(IT_0001*IT_0010 + IT_0004*IT_0023)) + IT_0010
      *IT_0024);
    const ccomplex_t IT_0041 = 1.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0001*IT_0010;
    const ccomplex_t IT_0043 = IT_0004*IT_0023;
    const ccomplex_t IT_0044 = IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_01) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_11));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(e_em*IT_0045 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0041*IT_0047;
    const ccomplex_t IT_0049 = IT_0039*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = pow(m_h, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0018 + IT_0051 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0053 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_00) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_10)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_00) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_10));
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(e_em*IT_0053 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_01) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_11));
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(e_em*IT_0056 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0055*IT_0058;
    const ccomplex_t IT_0060 = IT_0052*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_13 + IT_0051 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0063 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_00) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_10)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_00) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_10));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(e_em*IT_0063 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_01) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_11));
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(e_em*IT_0066 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0065*IT_0068;
    const ccomplex_t IT_0070 = IT_0062*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0051 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0073 = 2*U_stau_01*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_01) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) + U_stau_11*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_11));
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(e_em*IT_0073 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = 2*U_stau_00*((IT_0001*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.5*m_W*(IT_0007*IT_0012 + -IT_0021)*IT_0044)*conj(U_stau_01) +
       1./2*IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_11)) + U_stau_10*
      (IT_0004*IT_0021*IT_0028*IT_0029*conj(U_stau_01) + 2*(IT_0001*IT_0009
      *IT_0021*IT_0027*IT_0029 + -m_W*IT_0007*IT_0012*IT_0044)*conj(U_stau_11));
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(e_em*IT_0076 +
       1.4142135623731*A_tau*IT_0001*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0075*IT_0078;
    const ccomplex_t IT_0080 = IT_0072*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0017 + IT_0051 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0083 = IT_0079*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0016 + -IT_0038 + -IT_0050 + IT_0061 +
       IT_0071 + IT_0081 + IT_0084;
    return create_ccomplex_return(IT_0085*conj(IT_0085));
}

