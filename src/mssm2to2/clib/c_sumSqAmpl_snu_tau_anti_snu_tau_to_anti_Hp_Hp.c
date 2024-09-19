#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_tau, 2);
    const ccomplex_t IT_0011 = cpow(IT_0001, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + -IT_0008*(IT_0004 + -IT_0006*(1 + (-2)*IT_0009*IT_0010*IT_0011
      )));
    const ccomplex_t IT_0013 = 0.25*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cpow(IT_0001, -1);
    const ccomplex_t IT_0016 = cpow(IT_0005, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(A_tau*IT_0007*conj
      (U_stau_11) + (-1.4142135623731)*e_em*IT_0016*(m_W*IT_0001*IT_0007*conj
      (U_stau_01) + (-0.5)*m_tau*IT_0014*(m_tau*IT_0007*IT_0015*conj(U_stau_01) 
      + mu_h*conj(U_stau_11))));
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(A_tau*IT_0007*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0016*(m_W*IT_0001*IT_0007*U_stau_01 + (-0.5)
      *m_tau*IT_0014*(m_tau*IT_0007*IT_0015*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = pow(m_Hp, 2);
    const ccomplex_t IT_0021 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0020 + IT_0021 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0003, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*-1)*e_em*(IT_0003*IT_0016 + 
      -IT_0005*IT_0025);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*e_em*IT_0016*IT_0025;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = s_14 + s_24;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = s_12 + IT_0021;
    const ccomplex_t IT_0033 = 0.5*IT_0026;
    const ccomplex_t IT_0034 = -IT_0027;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = IT_0031 + IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0034;
    const ccomplex_t IT_0039 = s_24*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0033;
    const ccomplex_t IT_0041 = IT_0032*IT_0040;
    const ccomplex_t IT_0042 = -IT_0039 + (-2)*IT_0041;
    const ccomplex_t IT_0043 = IT_0037 + IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0004*IT_0005;
    const ccomplex_t IT_0048 = -IT_0016 + IT_0047;
    const ccomplex_t IT_0049 = cos(alpha);
    const ccomplex_t IT_0050 = cpow(IT_0007, 3);
    const ccomplex_t IT_0051 = 3*IT_0016 + IT_0047;
    const ccomplex_t IT_0052 = sin(alpha);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_W*e_em*(IT_0048*IT_0049
      *IT_0050 + IT_0001*(IT_0008*IT_0051*IT_0052 + -IT_0001*(IT_0007*IT_0049
      *IT_0051 + IT_0001*IT_0048*IT_0052)));
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0016 + IT_0047;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_W*e_em*(IT_0007*IT_0049 +
       IT_0001*IT_0052)*IT_0055;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0054*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0016;
    const ccomplex_t IT_0063 = IT_0047 + IT_0062;
    const ccomplex_t IT_0064 = -IT_0001;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_W*e_em*(IT_0050*IT_0052
      *IT_0063 + IT_0064*(IT_0008*IT_0049*IT_0051 + -(IT_0007*IT_0051*IT_0052 + 
      -IT_0001*IT_0049*IT_0063)*IT_0064));
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_W*e_em*(IT_0001*IT_0049 + 
      -IT_0007*IT_0052)*IT_0055;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + (-2)*IT_0021 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(A_tau*IT_0007*conj
      (U_stau_10) + (-1.4142135623731)*e_em*IT_0016*(m_W*IT_0001*IT_0007*conj
      (U_stau_00) + (-0.5)*m_tau*IT_0014*(m_tau*IT_0007*IT_0015*conj(U_stau_00) 
      + mu_h*conj(U_stau_10))));
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(A_tau*IT_0007*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0016*(m_W*IT_0001*IT_0007*U_stau_00 + (-0.5)
      *m_tau*IT_0014*(m_tau*IT_0007*IT_0015*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0020 + IT_0021 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0013 + IT_0024 + IT_0046 + -IT_0061 + 
      -IT_0072 + IT_0078;
    return create_ccomplex_return(IT_0079*conj(IT_0079));
}

