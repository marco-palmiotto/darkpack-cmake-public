#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_anti_Hp_W(
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
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = e_em*mu_h*m_tau*IT_0004*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = e_em*mu_h*m_tau*IT_0004*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0012 +
       0.5*IT_0015;
    const ccomplex_t IT_0017 = pow(m_stau_1, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_00);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(A_tau*IT_0009*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0001*(m_W*IT_0009*IT_0024*U_stau_00 + (-0.5)
      *m_tau*IT_0004*(m_tau*IT_0009*IT_0025*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0027 = pow(m_Hp, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0017 + IT_0027 + reg_prop +
       -m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = cos(alpha);
    const ccomplex_t IT_0033 = sin(alpha);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0024
      *IT_0032 + IT_0009*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = mu_h*m_tau*IT_0004;
    const ccomplex_t IT_0037 = IT_0024*IT_0033;
    const ccomplex_t IT_0038 = IT_0009*IT_0032;
    const ccomplex_t IT_0039 = IT_0037 + IT_0038;
    const ccomplex_t IT_0040 = cos(theta_W);
    const ccomplex_t IT_0041 = cpow(IT_0040, -2);
    const ccomplex_t IT_0042 = pow(m_tau, 2);
    const ccomplex_t IT_0043 = 2*U_stau_00*(((-0.5)*m_W*IT_0039*(IT_0001 + 
      -IT_0000*IT_0041) + IT_0001*IT_0004*IT_0025*IT_0033*IT_0042)*conj
      (U_stau_00) + 1./2*IT_0001*IT_0025*IT_0032*IT_0036*conj(U_stau_10)) +
       U_stau_10*(IT_0001*IT_0025*IT_0032*IT_0036*conj(U_stau_00) + (-2)*(m_W
      *IT_0000*IT_0039*IT_0041 + -IT_0001*IT_0004*IT_0025*IT_0033*IT_0042)*conj
      (U_stau_10));
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(e_em*IT_0043 +
       1.4142135623731*A_tau*IT_0033*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0035*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0009
      *IT_0032 + -IT_0024*IT_0033);
    const ccomplex_t IT_0051 = -IT_0009*IT_0033;
    const ccomplex_t IT_0052 = IT_0024*IT_0032;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = (-2)*U_stau_00*((IT_0001*IT_0004*IT_0025
      *IT_0032*IT_0042 + (-0.5)*m_W*(IT_0001 + -IT_0000*IT_0041)*IT_0053)*conj
      (U_stau_00) + -1./2*IT_0001*IT_0025*IT_0033*IT_0036*conj(U_stau_10)) +
       U_stau_10*(IT_0001*IT_0025*IT_0033*IT_0036*conj(U_stau_00) + (-2)*
      (IT_0001*IT_0004*IT_0025*IT_0032*IT_0042 + -m_W*IT_0000*IT_0041*IT_0053)
      *conj(U_stau_10));
    const ccomplex_t IT_0055 = (0 + _Complex_I*-1)*(e_em*IT_0054 + (
      -1.4142135623731)*A_tau*IT_0032*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0050*IT_0059;
    const ccomplex_t IT_0061 = -IT_0021 + IT_0031 + -IT_0049 + -IT_0060;
    const ccomplex_t IT_0062 = pow(m_W, -2);
    const ccomplex_t IT_0063 = pow(s_24, 2);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (-2)*U_stau_00*((IT_0001*IT_0004*IT_0025
      *IT_0032*IT_0042 + (-0.5)*m_W*(IT_0001 + -IT_0000*IT_0041)*IT_0053)*conj
      (U_stau_00) + -1./2*IT_0001*IT_0025*IT_0033*IT_0036*conj(U_stau_10)) +
       U_stau_10*(IT_0001*IT_0025*IT_0033*IT_0036*conj(U_stau_00) + (-2)*
      (IT_0001*IT_0004*IT_0025*IT_0032*IT_0042 + -m_W*IT_0000*IT_0041*IT_0053)
      *conj(U_stau_10));
    const ccomplex_t IT_0066 = (0 + _Complex_I*-1)*(e_em*IT_0065 + (
      -1.4142135623731)*A_tau*IT_0032*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0057*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0050*IT_0069;
    const ccomplex_t IT_0071 = 2*U_stau_00*(((-0.5)*m_W*IT_0039*(IT_0001 + 
      -IT_0000*IT_0041) + IT_0001*IT_0004*IT_0025*IT_0033*IT_0042)*conj
      (U_stau_00) + 1./2*IT_0001*IT_0025*IT_0032*IT_0036*conj(U_stau_10)) +
       U_stau_10*(IT_0001*IT_0025*IT_0032*IT_0036*conj(U_stau_00) + (-2)*(m_W
      *IT_0000*IT_0039*IT_0041 + -IT_0001*IT_0004*IT_0025*IT_0033*IT_0042)*conj
      (U_stau_10));
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(e_em*IT_0071 +
       1.4142135623731*A_tau*IT_0033*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0046*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0035*IT_0075;
    const ccomplex_t IT_0077 = -IT_0021 + -IT_0070 + -IT_0076;
    const ccomplex_t IT_0078 = s_14*s_24*IT_0062;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = pow(s_14, 2);
    const ccomplex_t IT_0081 = IT_0062*IT_0080;
    return create_ccomplex_return(-IT_0061*(conj(IT_0061)*(IT_0017 + -IT_0064)
       + conj(IT_0077)*(s_12 + IT_0079)) + -IT_0077*(conj(IT_0061)*(s_12 +
       IT_0079) + conj(IT_0077)*(IT_0017 + -IT_0081)));
}

