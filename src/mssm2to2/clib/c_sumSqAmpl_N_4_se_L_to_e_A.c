#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_se_L_to_e_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_se_L_to_e_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_e, 2);
    const ccomplex_t IT_0001 = pow(m_N_4, 2);
    const ccomplex_t IT_0002 = pow(m_se_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 2*IT_0011;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (-2)*IT_0015 + -IT_0020;
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = conj(N_B4)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W4)*e_em;
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0031;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0017*IT_0031;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0016*IT_0036;
    const ccomplex_t IT_0038 = 2*IT_0034 + IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0000;
    const ccomplex_t IT_0040 = IT_0011*IT_0019;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0011*IT_0036;
    const ccomplex_t IT_0043 = (-2)*s_13;
    const ccomplex_t IT_0044 = s_34*IT_0043;
    const ccomplex_t IT_0045 = m_e*s_14*m_N_4;
    const ccomplex_t IT_0046 = -IT_0015;
    const ccomplex_t IT_0047 = s_13*s_14;
    const ccomplex_t IT_0048 = s_34*IT_0001;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = s_14*IT_0043;
    const ccomplex_t IT_0052 = IT_0001*IT_0043;
    const ccomplex_t IT_0053 = pow(s_13, 2);
    const ccomplex_t IT_0054 = (-2)*IT_0053;
    const ccomplex_t IT_0055 = -IT_0037;
    const ccomplex_t IT_0056 = m_e*m_N_4;
    const ccomplex_t IT_0057 = (-2)*s_14;
    const ccomplex_t IT_0058 = (-2)*s_34;
    const ccomplex_t IT_0059 = (-2)*IT_0001;
    const ccomplex_t IT_0060 = (conj(IT_0021)*IT_0038 + IT_0021*conj(IT_0038))
      *IT_0039 + IT_0043*(conj(IT_0020)*IT_0038 + IT_0020*conj(IT_0038) + conj
      (IT_0021)*IT_0055 + IT_0021*conj(IT_0055)) + (conj(IT_0020)*IT_0042 +
       IT_0020*conj(IT_0042) + conj(IT_0041)*IT_0055 + IT_0041*conj(IT_0055))
      *IT_0057 + (conj(IT_0038)*IT_0041 + IT_0038*conj(IT_0041) + conj(IT_0021)
      *IT_0042 + IT_0021*conj(IT_0042))*IT_0058 + (conj(IT_0020)*IT_0055 +
       IT_0020*conj(IT_0055))*IT_0059;
    const ccomplex_t IT_0061 = m_e*s_34*m_N_4;
    const ccomplex_t IT_0062 = (-2)*conj(IT_0021);
    const ccomplex_t IT_0063 = s_14*IT_0000;
    const ccomplex_t IT_0064 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0065 = s_14*s_34;
    const ccomplex_t IT_0066 = (-2)*IT_0038;
    const ccomplex_t IT_0067 = (-0.5)*IT_0063;
    const ccomplex_t IT_0068 = (-2)*IT_0021;
    const ccomplex_t IT_0069 = (-0.5)*IT_0061;
    const ccomplex_t IT_0070 = s_13*(IT_0021*conj(IT_0021) + IT_0038*conj
      (IT_0038))*IT_0039 + (conj(IT_0021)*IT_0041 + IT_0021*conj(IT_0041) + conj
      (IT_0038)*IT_0042 + IT_0038*conj(IT_0042))*IT_0044 + (-2)*conj(IT_0020)*
      (IT_0034*IT_0045 + 2*IT_0046*IT_0050 + (-0.5)*IT_0041*IT_0051 + (-0.5)
      *IT_0020*IT_0052 + (-0.5)*IT_0021*IT_0054) + IT_0020*(conj(IT_0041)
      *IT_0051 + conj(IT_0021)*IT_0054) + (conj(IT_0042)*IT_0051 + conj(IT_0038)
      *IT_0054)*IT_0055 + (-2)*(IT_0045*IT_0046 + 2*IT_0034*IT_0050 + (-0.5)
      *IT_0042*IT_0051 + (-0.5)*IT_0038*IT_0054 + (-0.5)*IT_0052*IT_0055)*conj
      (IT_0055) + IT_0056*IT_0060 + IT_0046*(IT_0062*IT_0063 + IT_0061*IT_0064) 
      + IT_0034*(IT_0061*IT_0062 + IT_0063*IT_0064) + (-4)*conj(IT_0046)*
      (IT_0020*IT_0050 + 1./2*IT_0045*IT_0055 + (-2)*IT_0046*IT_0065 + 1./2
      *IT_0067*IT_0068 + 1./2*IT_0066*IT_0069) + (-2)*conj(IT_0034)*(IT_0020
      *IT_0045 + 2*IT_0050*IT_0055 + (-4)*IT_0034*IT_0065 + IT_0066*IT_0067 +
       IT_0068*IT_0069);
    return create_ccomplex_return(IT_0070);
}

