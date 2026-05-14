#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_snu_tau_to_H0_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_snu_tau_to_H0_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = cpow(IT_0002, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(A_tau*IT_0001*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0005*(m_W*IT_0001*IT_0002*U_stau_01 + (-0.5)
      *m_tau*IT_0007*(m_tau*IT_0001*IT_0008*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0010 = pow(m_stau_2, 2);
    const ccomplex_t IT_0011 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *U_stau_01;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0010 + -IT_0011 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0017*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_tau, 2);
    const ccomplex_t IT_0026 = cpow(IT_0002, 2);
    const ccomplex_t IT_0027 = cpow(IT_0001, 2);
    const ccomplex_t IT_0028 = A_tau*IT_0002*U_stau_11 + 0.707106781186548
      *e_em*IT_0005*((IT_0007*IT_0025 + -m_W*(IT_0026 + -IT_0027))*U_stau_01 + 
      -mu_h*m_tau*IT_0001*IT_0007*IT_0008*U_stau_11);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0024*IT_0034;
    const ccomplex_t IT_0036 = -IT_0015 + -IT_0022 + -IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, -2);
    const ccomplex_t IT_0038 = pow(s_14, 2);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = -IT_0016;
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -2);
    const ccomplex_t IT_0043 = IT_0004*IT_0042;
    const ccomplex_t IT_0044 = IT_0005 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0002 + 
      -IT_0001*IT_0003)*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = pow(m_H0, 2);
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0047 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0040*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = s_14*s_34*IT_0037;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = s_13 + IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *U_stau_00;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = mu_h*m_tau*IT_0007;
    const ccomplex_t IT_0059 = -IT_0001*IT_0003;
    const ccomplex_t IT_0060 = IT_0000*IT_0002;
    const ccomplex_t IT_0061 = IT_0059 + IT_0060;
    const ccomplex_t IT_0062 = (-2)*conj(U_stau_00)*((IT_0000*IT_0005*IT_0007
      *IT_0008*IT_0025 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0042)*IT_0061)
      *U_stau_01 + -1./2*IT_0003*IT_0005*IT_0008*IT_0058*U_stau_11) + conj
      (U_stau_10)*(IT_0003*IT_0005*IT_0008*IT_0058*U_stau_01 + (-2)*(IT_0000
      *IT_0005*IT_0007*IT_0008*IT_0025 + -m_W*IT_0004*IT_0042*IT_0061)*U_stau_11);
    const ccomplex_t IT_0063 = (0 + _Complex_I*-1)*(e_em*IT_0062 + (
      -1.4142135623731)*A_tau*IT_0000*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0010 + IT_0047 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0057*IT_0067;
    const ccomplex_t IT_0069 = (-2)*U_stau_01*((IT_0000*IT_0005*IT_0007
      *IT_0008*IT_0025 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0042)*IT_0061)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0005*IT_0008*IT_0058*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0005*IT_0008*IT_0058*conj(U_stau_01) + (-2)*
      (IT_0000*IT_0005*IT_0007*IT_0008*IT_0025 + -m_W*IT_0004*IT_0042*IT_0061)
      *conj(U_stau_11));
    const ccomplex_t IT_0070 = (0 + _Complex_I*-1)*(e_em*IT_0069 + (
      -1.4142135623731)*A_tau*IT_0000*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0047 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0040*IT_0074;
    const ccomplex_t IT_0076 = -IT_0040;
    const ccomplex_t IT_0077 = IT_0017 + IT_0076;
    const ccomplex_t IT_0078 = IT_0021*IT_0077;
    const ccomplex_t IT_0079 = -IT_0015 + -IT_0035 + IT_0051 + IT_0068 +
       IT_0075 + -IT_0078;
    const ccomplex_t IT_0080 = s_14*s_24*IT_0037;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = s_12 + IT_0081;
    const ccomplex_t IT_0083 = -IT_0036;
    const ccomplex_t IT_0084 = pow(s_34, 2);
    const ccomplex_t IT_0085 = IT_0037*IT_0084;
    const ccomplex_t IT_0086 = s_24*s_34*IT_0037;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = s_23 + IT_0087;
    const ccomplex_t IT_0089 = -conj(IT_0036);
    const ccomplex_t IT_0090 = pow(s_24, 2);
    const ccomplex_t IT_0091 = IT_0037*IT_0090;
    const ccomplex_t IT_0092 = (conj(IT_0036)*(IT_0010 + -IT_0039) + conj
      (IT_0052)*IT_0055 + conj(IT_0079)*IT_0082)*IT_0083 + -IT_0052*(conj
      (IT_0052)*(IT_0047 + -IT_0085) + conj(IT_0079)*IT_0088 + -IT_0055*IT_0089)
       + -IT_0079*(conj(IT_0052)*IT_0088 + -IT_0082*IT_0089 + conj(IT_0079)*
      (IT_0011 + -IT_0091));
    return create_ccomplex_return(IT_0092);
}

