#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_Hp_Hp(
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
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
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
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = IT_0007 + -IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -2);
    const ccomplex_t IT_0012 = pow(m_tau, 2);
    const ccomplex_t IT_0013 = cpow(IT_0006, -2);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0000*(IT_0005*IT_0010
      *U_stau_00*conj(U_stau_01) + (-2)*(IT_0002*IT_0010 + IT_0004*IT_0009
      *IT_0011*IT_0012*IT_0013)*U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0015 = 0.25*IT_0014;
    const ccomplex_t IT_0016 = pow(m_stau_1, 2);
    const ccomplex_t IT_0017 = pow(m_stau_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0002*IT_0003;
    const ccomplex_t IT_0020 = cpow(IT_0003, -1);
    const ccomplex_t IT_0021 = IT_0019 + -IT_0020;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = cpow(IT_0008, 3);
    const ccomplex_t IT_0024 = IT_0019 + 3*IT_0020;
    const ccomplex_t IT_0025 = sin(alpha);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_W*e_em*(IT_0021*IT_0022
      *IT_0023 + IT_0006*(IT_0009*IT_0024*IT_0025 + -IT_0006*(IT_0008*IT_0022
      *IT_0024 + IT_0006*IT_0021*IT_0025)));
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = mu_h*m_tau*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0006, -1);
    const ccomplex_t IT_0031 = IT_0006*IT_0025;
    const ccomplex_t IT_0032 = IT_0008*IT_0022;
    const ccomplex_t IT_0033 = IT_0031 + IT_0032;
    const ccomplex_t IT_0034 = 2*U_stau_00*((IT_0012*IT_0020*IT_0025*IT_0028
      *IT_0030 + 0.5*m_W*(IT_0002*IT_0003 + -IT_0020)*IT_0033)*conj(U_stau_01) +
       1./2*IT_0020*IT_0022*IT_0029*IT_0030*conj(U_stau_11)) + U_stau_10*
      (IT_0020*IT_0022*IT_0029*IT_0030*conj(U_stau_01) + 2*(IT_0012*IT_0020
      *IT_0025*IT_0028*IT_0030 + -m_W*IT_0002*IT_0003*IT_0033)*conj(U_stau_11));
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(e_em*IT_0034 +
       1.4142135623731*A_tau*IT_0025*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = IT_0018*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + -IT_0016 + -IT_0017 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0041 = -IT_0020;
    const ccomplex_t IT_0042 = IT_0019 + IT_0041;
    const ccomplex_t IT_0043 = -IT_0006;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_W*e_em*(IT_0023*IT_0025
      *IT_0042 + IT_0043*(IT_0009*IT_0022*IT_0024 + -(IT_0008*IT_0024*IT_0025 + 
      -IT_0006*IT_0022*IT_0042)*IT_0043));
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = -IT_0008*IT_0025;
    const ccomplex_t IT_0047 = IT_0006*IT_0022;
    const ccomplex_t IT_0048 = IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = (-2)*U_stau_00*((IT_0012*IT_0020*IT_0022
      *IT_0028*IT_0030 + 0.5*m_W*(IT_0002*IT_0003 + -IT_0020)*IT_0048)*conj
      (U_stau_01) + -1./2*IT_0020*IT_0025*IT_0029*IT_0030*conj(U_stau_11)) +
       U_stau_10*(IT_0020*IT_0025*IT_0029*IT_0030*conj(U_stau_01) + (-2)*
      (IT_0012*IT_0020*IT_0022*IT_0028*IT_0030 + -m_W*IT_0002*IT_0003*IT_0048)
      *conj(U_stau_11));
    const ccomplex_t IT_0050 = (0 + _Complex_I*-1)*(e_em*IT_0049 + (
      -1.4142135623731)*A_tau*IT_0022*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = IT_0040*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0056 = cpow(IT_0001, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*-1)*e_em*(IT_0001*IT_0020 + 
      -IT_0003*IT_0056);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0020 + -IT_0003*IT_0056)*U_stau_00*conj(U_stau_01) + IT_0003*IT_0056
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = s_14 + s_24;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0057;
    const ccomplex_t IT_0064 = (-2)*IT_0058;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = s_12 + IT_0017;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = IT_0062 + IT_0067;
    const ccomplex_t IT_0069 = IT_0057*IT_0064;
    const ccomplex_t IT_0070 = s_24*IT_0069;
    const ccomplex_t IT_0071 = IT_0059*IT_0063;
    const ccomplex_t IT_0072 = 0.5*IT_0016 + 0.5*IT_0017;
    const ccomplex_t IT_0073 = s_12 + IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = -IT_0070 + (-2)*IT_0074;
    const ccomplex_t IT_0076 = IT_0068 + IT_0075;
    const ccomplex_t IT_0077 = IT_0055*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(A_tau*IT_0008*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0020*(m_W*IT_0006*IT_0008*U_stau_00 + (-0.5)
      *m_tau*IT_0028*(m_tau*IT_0008*IT_0030*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(A_tau*IT_0008*conj
      (U_stau_11) + (-1.4142135623731)*e_em*IT_0020*(m_W*IT_0006*IT_0008*conj
      (U_stau_01) + (-0.5)*m_tau*IT_0028*(m_tau*IT_0008*IT_0030*conj(U_stau_01) 
      + mu_h*conj(U_stau_11))));
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = pow(m_Hp, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0016 + IT_0082 + reg_prop +
       -m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0015 + -IT_0039 + -IT_0054 + IT_0078 +
       IT_0085;
    return create_ccomplex_return(IT_0086*conj(IT_0086));
}

