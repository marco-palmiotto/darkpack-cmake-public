#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_stau_1_to_tau_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_stau_1_to_tau_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = N_d4*e_em*m_tau*IT_0010*U_stau_00;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = pow(m_N_4, 2);
    const ccomplex_t IT_0020 = pow(m_tau, 2);
    const ccomplex_t IT_0021 = pow(m_stau_1, 2);
    const ccomplex_t IT_0022 = cpow(s_12 + 0.5*IT_0019 + (-0.5)*IT_0020 + 0.5
      *IT_0021 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = 2*IT_0016;
    const ccomplex_t IT_0026 = cpow(s_13 + (-0.5)*IT_0019 + (-0.5)*IT_0020 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = IT_0024 + 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0016*IT_0028;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0037 = IT_0008*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_d4)*e_em*m_tau*IT_0010*U_stau_10;
    const ccomplex_t IT_0040 = IT_0009*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0035 + IT_0038 + 
      -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = IT_0022*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0026*IT_0043;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0025*IT_0048;
    const ccomplex_t IT_0050 = -IT_0046 + (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0016*IT_0048;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (-2)*s_13;
    const ccomplex_t IT_0054 = s_34*IT_0053;
    const ccomplex_t IT_0055 = IT_0020*IT_0053;
    const ccomplex_t IT_0056 = (-0.5)*IT_0029;
    const ccomplex_t IT_0057 = (-0.5)*IT_0046;
    const ccomplex_t IT_0058 = s_14*m_N_4*m_tau;
    const ccomplex_t IT_0059 = 0.5*IT_0024;
    const ccomplex_t IT_0060 = s_13*s_14;
    const ccomplex_t IT_0061 = s_34*IT_0019;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = s_14*IT_0053;
    const ccomplex_t IT_0065 = pow(s_13, 2);
    const ccomplex_t IT_0066 = (-2)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0049;
    const ccomplex_t IT_0068 = m_N_4*m_tau;
    const ccomplex_t IT_0069 = (-2)*s_14;
    const ccomplex_t IT_0070 = (-2)*s_34;
    const ccomplex_t IT_0071 = (-2)*IT_0019;
    const ccomplex_t IT_0072 = (-2)*IT_0020;
    const ccomplex_t IT_0073 = IT_0053*(conj(IT_0050)*IT_0056 + IT_0050*conj
      (IT_0056) + conj(IT_0030)*IT_0067 + IT_0030*conj(IT_0067)) + (conj(IT_0052
      )*IT_0056 + IT_0052*conj(IT_0056) + conj(IT_0032)*IT_0067 + IT_0032*conj
      (IT_0067))*IT_0069 + (conj(IT_0032)*IT_0050 + IT_0032*conj(IT_0050) + conj
      (IT_0030)*IT_0052 + IT_0030*conj(IT_0052))*IT_0070 + (conj(IT_0056)
      *IT_0067 + IT_0056*conj(IT_0067))*IT_0071 + (conj(IT_0030)*IT_0050 +
       IT_0030*conj(IT_0050))*IT_0072;
    const ccomplex_t IT_0074 = s_34*m_N_4*m_tau;
    const ccomplex_t IT_0075 = (-2)*conj(IT_0030);
    const ccomplex_t IT_0076 = s_14*IT_0020;
    const ccomplex_t IT_0077 = (-2)*conj(IT_0050);
    const ccomplex_t IT_0078 = s_14*s_34;
    const ccomplex_t IT_0079 = (-2)*IT_0050;
    const ccomplex_t IT_0080 = (-0.5)*IT_0076;
    const ccomplex_t IT_0081 = (-2)*IT_0030;
    const ccomplex_t IT_0082 = (-0.5)*IT_0074;
    const ccomplex_t IT_0083 = (conj(IT_0030)*IT_0032 + IT_0030*conj(IT_0032) 
      + conj(IT_0050)*IT_0052 + IT_0050*conj(IT_0052))*IT_0054 + (IT_0030*conj
      (IT_0030) + IT_0050*conj(IT_0050))*IT_0055 + (-2)*conj(IT_0056)*((-0.5)
      *IT_0019*IT_0053*IT_0056 + IT_0057*IT_0058 + 2*IT_0059*IT_0063 + (-0.5)
      *IT_0032*IT_0064 + (-0.5)*IT_0030*IT_0066) + IT_0056*(conj(IT_0032)
      *IT_0064 + conj(IT_0030)*IT_0066) + (conj(IT_0052)*IT_0064 + conj(IT_0050)
      *IT_0066)*IT_0067 + (-2)*(IT_0058*IT_0059 + 2*IT_0057*IT_0063 + (-0.5)
      *IT_0052*IT_0064 + (-0.5)*IT_0050*IT_0066 + (-0.5)*IT_0019*IT_0053*IT_0067
      )*conj(IT_0067) + IT_0068*IT_0073 + IT_0059*(IT_0075*IT_0076 + IT_0074
      *IT_0077) + IT_0057*(IT_0074*IT_0075 + IT_0076*IT_0077) + (-4)*conj
      (IT_0059)*(IT_0056*IT_0063 + 1./2*IT_0058*IT_0067 + (-2)*IT_0059*IT_0078 +
       1./2*IT_0080*IT_0081 + 1./2*IT_0079*IT_0082) + (-2)*conj(IT_0057)*
      (IT_0056*IT_0058 + 2*IT_0063*IT_0067 + (-4)*IT_0057*IT_0078 + IT_0079
      *IT_0080 + IT_0081*IT_0082);
    return create_ccomplex_return(IT_0083);
}

