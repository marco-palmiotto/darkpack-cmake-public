#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_H0_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_H0_H0(
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
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_tau, 2);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (-0.5)*(IT_0004*(IT_0002 + -IT_0012) + -IT_0006
      *(IT_0002 + -(1 + (-2)*IT_0007*IT_0008*IT_0010)*IT_0012))*U_stau_00*conj
      (U_stau_00) + (IT_0002*IT_0004 + -IT_0006*(IT_0002 + -IT_0007*IT_0008
      *IT_0010*IT_0012))*U_stau_10*conj(U_stau_10);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0000*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = mu_h*m_tau*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0009, -1);
    const ccomplex_t IT_0019 = cpow(IT_0011, -1);
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = -IT_0003*IT_0020;
    const ccomplex_t IT_0022 = IT_0005*IT_0009;
    const ccomplex_t IT_0023 = IT_0021 + IT_0022;
    const ccomplex_t IT_0024 = (-2)*conj(U_stau_00)*((IT_0005*IT_0008*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0002*IT_0011 + -IT_0019)*IT_0023)*U_stau_01
       + -1./2*IT_0003*IT_0017*IT_0018*IT_0019*U_stau_11) + conj(U_stau_10)*
      (IT_0003*IT_0017*IT_0018*IT_0019*U_stau_01 + (-2)*(IT_0005*IT_0008*IT_0016
      *IT_0018*IT_0019 + -m_W*IT_0002*IT_0011*IT_0023)*U_stau_11);
    const ccomplex_t IT_0025 = (0 + _Complex_I*-1)*(e_em*IT_0024 + (
      -1.4142135623731)*A_tau*IT_0005*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = (-2)*U_stau_00*((IT_0005*IT_0008*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0002*IT_0011 + -IT_0019)*IT_0023)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0017*IT_0018*IT_0019*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0017*IT_0018*IT_0019*conj(U_stau_01) + (-2)*
      (IT_0005*IT_0008*IT_0016*IT_0018*IT_0019 + -m_W*IT_0002*IT_0011*IT_0023)
      *conj(U_stau_11));
    const ccomplex_t IT_0028 = (0 + _Complex_I*-1)*(e_em*IT_0027 + (
      -1.4142135623731)*A_tau*IT_0005*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = pow(m_H0, 2);
    const ccomplex_t IT_0032 = pow(m_stau_1, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0031 + IT_0032 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + (-2)*IT_0032 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0002*IT_0011;
    const ccomplex_t IT_0038 = IT_0019 + IT_0037;
    const ccomplex_t IT_0039 = cpow(IT_0003, 3);
    const ccomplex_t IT_0040 = (0 + _Complex_I*-1)*m_W*e_em*IT_0038*(IT_0005*(
      (IT_0004 + -IT_0006)*IT_0009 + IT_0003*IT_0005*IT_0020) + -IT_0020*IT_0039);
    const ccomplex_t IT_0041 = (-1.5)*IT_0040;
    const ccomplex_t IT_0042 = (-2)*U_stau_00*((IT_0005*IT_0008*IT_0016
      *IT_0018*IT_0019 + 0.5*m_W*(IT_0002*IT_0011 + -IT_0019)*IT_0023)*conj
      (U_stau_00) + -1./2*IT_0003*IT_0017*IT_0018*IT_0019*conj(U_stau_10)) +
       U_stau_10*(IT_0003*IT_0017*IT_0018*IT_0019*conj(U_stau_00) + (-2)*
      (IT_0005*IT_0008*IT_0016*IT_0018*IT_0019 + -m_W*IT_0002*IT_0011*IT_0023)
      *conj(U_stau_10));
    const ccomplex_t IT_0043 = (0 + _Complex_I*-1)*(e_em*IT_0042 + (
      -1.4142135623731)*A_tau*IT_0005*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = IT_0036*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0031 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0049 = cpow(IT_0044, 2);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*5)*m_W*e_em*IT_0038*(IT_0005*
      (IT_0004*IT_0020 + (-1)*IT_0005*(IT_0003*IT_0009 + 1./5*IT_0005*IT_0020)) 
      + 1./5*IT_0009*IT_0039);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0003*IT_0009;
    const ccomplex_t IT_0056 = IT_0005*IT_0020;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = 2*U_stau_00*((IT_0003*IT_0008*IT_0016*IT_0018
      *IT_0019 + 0.5*m_W*(IT_0002*IT_0011 + -IT_0019)*IT_0057)*conj(U_stau_00) +
       1./2*IT_0005*IT_0017*IT_0018*IT_0019*conj(U_stau_10)) + U_stau_10*
      (IT_0005*IT_0017*IT_0018*IT_0019*conj(U_stau_00) + 2*(IT_0003*IT_0008
      *IT_0016*IT_0018*IT_0019 + -m_W*IT_0002*IT_0011*IT_0057)*conj(U_stau_10));
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(e_em*IT_0058 +
       1.4142135623731*A_tau*IT_0003*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = IT_0052*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0031 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0065 = IT_0050*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0031 + IT_0032 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0067 = IT_0030*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0015 + IT_0035 + -IT_0047 + IT_0051 + 
      -IT_0063 + IT_0065 + IT_0068;
    return create_ccomplex_return(IT_0069*conj(IT_0069));
}

