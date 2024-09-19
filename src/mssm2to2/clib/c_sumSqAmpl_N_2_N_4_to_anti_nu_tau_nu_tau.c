#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_4_to_anti_nu_tau_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_4_to_anti_nu_tau_nu_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = N_B4*e_em;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = N_W4*e_em;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + -IT_0017);
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_N_4, 2);
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0021 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0005;
    const ccomplex_t IT_0027 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0000*IT_0006;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = e_em*IT_0026;
    const ccomplex_t IT_0040 = e_em*IT_0029;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0039 + IT_0040);
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = pow(m_N_2, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0021 + -IT_0043 + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0038*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = s_14*s_23;
    const ccomplex_t IT_0050 = s_13*s_24;
    const ccomplex_t IT_0051 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0052 = IT_0026*IT_0051;
    const ccomplex_t IT_0053 = IT_0029*IT_0051;
    const ccomplex_t IT_0054 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0055 = IT_0029*IT_0054;
    const ccomplex_t IT_0056 = IT_0026*IT_0054;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0052 + IT_0053 + 
      -IT_0055 + -IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0046*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = N_B2*e_em;
    const ccomplex_t IT_0063 = IT_0001*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = N_W2*e_em;
    const ccomplex_t IT_0066 = IT_0006*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B4)*e_em;
    const ccomplex_t IT_0071 = IT_0001*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W4)*e_em;
    const ccomplex_t IT_0074 = IT_0006*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + -IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0069*IT_0077;
    const ccomplex_t IT_0079 = cpow(s_13 + (-0.5)*IT_0043 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = 2*IT_0061 + IT_0082;
    const ccomplex_t IT_0084 = 4*IT_0082;
    const ccomplex_t IT_0085 = 0.125*IT_0084;
    const ccomplex_t IT_0086 = IT_0061 + IT_0085;
    const ccomplex_t IT_0087 = s_34*m_N_2*m_N_4;
    const ccomplex_t IT_0088 = 4*(IT_0025 + 2*IT_0048)*(conj(IT_0025) + 2*conj
      (IT_0048))*IT_0049 + 8*IT_0050*(conj(IT_0061)*IT_0083 + conj(IT_0082)
      *IT_0086) + 4*((IT_0025 + 2*IT_0048)*(conj(IT_0061) + 1./2*conj(IT_0082)) 
      + 1./2*conj(IT_0025)*IT_0083 + 2*conj(IT_0048)*IT_0086)*IT_0087;
    return create_ccomplex_return(IT_0088);
}

