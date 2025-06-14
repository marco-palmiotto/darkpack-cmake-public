#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_snu_tau_to_h_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_snu_tau_to_h_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = cpow(IT_0000, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(A_tau*IT_0002*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0005*(m_W*IT_0000*IT_0002*U_stau_01 + (-0.5)
      *m_tau*IT_0008*(m_tau*IT_0002*IT_0009*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0011 = pow(m_stau_2, 2);
    const ccomplex_t IT_0012 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *U_stau_01;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_W*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = pow(m_tau, 2);
    const ccomplex_t IT_0027 = cpow(IT_0000, 2);
    const ccomplex_t IT_0028 = cpow(IT_0002, 2);
    const ccomplex_t IT_0029 = A_tau*IT_0000*U_stau_11 + 0.707106781186548
      *e_em*IT_0005*((IT_0008*IT_0026 + -m_W*(IT_0027 + -IT_0028))*U_stau_01 + 
      -mu_h*m_tau*IT_0002*IT_0008*IT_0009*U_stau_11);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0032 
      + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = -IT_0016 + -IT_0023 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -2);
    const ccomplex_t IT_0039 = pow(s_14, 2);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = -IT_0017;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -2);
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = IT_0005 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_W*e_em*(IT_0001*IT_0002 +
       IT_0000*IT_0003)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = pow(m_h, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0048 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0041*IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = s_14*s_34*IT_0038;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = s_13 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *U_stau_00;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = mu_h*m_tau*IT_0008;
    const ccomplex_t IT_0060 = IT_0000*IT_0003;
    const ccomplex_t IT_0061 = IT_0001*IT_0002;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = 2*conj(U_stau_00)*((IT_0003*IT_0005*IT_0008
      *IT_0009*IT_0026 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0043)*IT_0062)
      *U_stau_01 + 1./2*IT_0001*IT_0005*IT_0009*IT_0059*U_stau_11) + conj
      (U_stau_10)*(IT_0001*IT_0005*IT_0009*IT_0059*U_stau_01 + 2*(IT_0003
      *IT_0005*IT_0008*IT_0009*IT_0026 + -m_W*IT_0004*IT_0043*IT_0062)*U_stau_11);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(e_em*IT_0063 +
       1.4142135623731*A_tau*IT_0003*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0011 + IT_0048 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0058*IT_0068;
    const ccomplex_t IT_0070 = 2*U_stau_01*((IT_0003*IT_0005*IT_0008*IT_0009
      *IT_0026 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0043)*IT_0062)*conj(U_stau_01
      ) + 1./2*IT_0001*IT_0005*IT_0009*IT_0059*conj(U_stau_11)) + U_stau_11*
      (IT_0001*IT_0005*IT_0009*IT_0059*conj(U_stau_01) + 2*(IT_0003*IT_0005
      *IT_0008*IT_0009*IT_0026 + -m_W*IT_0004*IT_0043*IT_0062)*conj(U_stau_11));
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(e_em*IT_0070 +
       1.4142135623731*A_tau*IT_0003*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0048 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0041*IT_0075;
    const ccomplex_t IT_0077 = -IT_0041;
    const ccomplex_t IT_0078 = IT_0018 + IT_0077;
    const ccomplex_t IT_0079 = IT_0022*IT_0078;
    const ccomplex_t IT_0080 = -IT_0016 + -IT_0036 + IT_0052 + IT_0069 +
       IT_0076 + -IT_0079;
    const ccomplex_t IT_0081 = s_14*s_24*IT_0038;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = s_12 + IT_0082;
    const ccomplex_t IT_0084 = -IT_0037;
    const ccomplex_t IT_0085 = pow(s_34, 2);
    const ccomplex_t IT_0086 = IT_0038*IT_0085;
    const ccomplex_t IT_0087 = s_24*s_34*IT_0038;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = s_23 + IT_0088;
    const ccomplex_t IT_0090 = -conj(IT_0037);
    const ccomplex_t IT_0091 = pow(s_24, 2);
    const ccomplex_t IT_0092 = IT_0038*IT_0091;
    const ccomplex_t IT_0093 = (conj(IT_0037)*(IT_0011 + -IT_0040) + conj
      (IT_0053)*IT_0056 + conj(IT_0080)*IT_0083)*IT_0084 + -IT_0053*(conj
      (IT_0053)*(IT_0048 + -IT_0086) + conj(IT_0080)*IT_0089 + -IT_0056*IT_0090)
       + -IT_0080*(conj(IT_0053)*IT_0089 + -IT_0083*IT_0090 + conj(IT_0080)*
      (IT_0012 + -IT_0092));
    return create_ccomplex_return(IT_0093);
}

