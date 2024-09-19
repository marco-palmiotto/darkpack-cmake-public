#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = mu_h*m_tau*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0002, -1);
    const ccomplex_t IT_0010 = -IT_0001*IT_0003;
    const ccomplex_t IT_0011 = IT_0000*IT_0002;
    const ccomplex_t IT_0012 = IT_0010 + IT_0011;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -2);
    const ccomplex_t IT_0015 = pow(m_tau, 2);
    const ccomplex_t IT_0016 = (-2)*U_stau_01*(((-0.5)*m_W*IT_0012*(IT_0005 + 
      -IT_0004*IT_0014) + IT_0000*IT_0005*IT_0007*IT_0009*IT_0015)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0005*IT_0008*IT_0009*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0005*IT_0008*IT_0009*conj(U_stau_01) + 2*(m_W
      *IT_0004*IT_0012*IT_0014 + -IT_0000*IT_0005*IT_0007*IT_0009*IT_0015)*conj
      (U_stau_11));
    const ccomplex_t IT_0017 = (0 + _Complex_I*-1)*(e_em*IT_0016 + (
      -1.4142135623731)*A_tau*IT_0000*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_stau_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + (-2)*IT_0019 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0002*IT_0003;
    const ccomplex_t IT_0027 = IT_0000*IT_0001;
    const ccomplex_t IT_0028 = IT_0026 + IT_0027;
    const ccomplex_t IT_0029 = 2*U_stau_01*((IT_0003*IT_0005*IT_0007*IT_0009
      *IT_0015 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0014)*IT_0028)*conj(U_stau_01
      ) + 1./2*IT_0000*IT_0005*IT_0008*IT_0009*conj(U_stau_11)) + U_stau_11*
      (IT_0000*IT_0005*IT_0008*IT_0009*conj(U_stau_01) + 2*(IT_0003*IT_0005
      *IT_0007*IT_0009*IT_0015 + -m_W*IT_0004*IT_0014*IT_0028)*conj(U_stau_11));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(e_em*IT_0029 +
       1.4142135623731*A_tau*IT_0003*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0025*IT_0034;
    const ccomplex_t IT_0036 = e_em*IT_0005;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = e_em*mu_h*m_tau*IT_0007*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0039 = IT_0005*IT_0038;
    const ccomplex_t IT_0040 = e_em*mu_h*m_tau*IT_0007*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0041 = IT_0005*IT_0040;
    const ccomplex_t IT_0042 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0043 = IT_0001*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0039 + 0.5*IT_0041 + (-0.5)*IT_0044 +
       0.5*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + (-2)*IT_0019 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0037*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *conj(U_stau_01);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(A_tau*IT_0001*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0005*(m_W*IT_0001*IT_0002*U_stau_01 + (-0.5)
      *m_tau*IT_0007*(m_tau*IT_0001*IT_0009*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0056 = pow(m_Hp, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_13 + IT_0019 + IT_0056 + reg_prop +
       -m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0054*IT_0059;
    const ccomplex_t IT_0061 = -IT_0023 + -IT_0035 + -IT_0052 + IT_0060;
    const ccomplex_t IT_0062 = pow(m_W, -2);
    const ccomplex_t IT_0063 = pow(s_24, 2);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = -IT_0023 + -IT_0035 + -IT_0052;
    const ccomplex_t IT_0066 = s_14*s_24*IT_0062;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = pow(s_14, 2);
    const ccomplex_t IT_0069 = IT_0062*IT_0068;
    return create_ccomplex_return(-IT_0061*(conj(IT_0061)*(IT_0019 + -IT_0064)
       + conj(IT_0065)*(s_12 + IT_0067)) + -IT_0065*(conj(IT_0061)*(s_12 +
       IT_0067) + conj(IT_0065)*(IT_0019 + -IT_0069)));
}

