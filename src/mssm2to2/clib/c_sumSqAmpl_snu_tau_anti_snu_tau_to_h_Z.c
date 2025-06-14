#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_snu_tau_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_snu_tau_to_h_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0003;
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0000, -2);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = IT_0003 + IT_0007;
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = sin(alpha);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_W*e_em*IT_0008*(IT_0009
      *IT_0010 + IT_0011*IT_0012);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0015 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0005*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = pow(m_Z, -2);
    const ccomplex_t IT_0022 = pow(s_34, 2);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = m_W*e_em;
    const ccomplex_t IT_0025 = cpow(IT_0000, -4);
    const ccomplex_t IT_0026 = IT_0002*IT_0009*IT_0010*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0011*IT_0012*IT_0025;
    const ccomplex_t IT_0028 = cpow(IT_0002, -3);
    const ccomplex_t IT_0029 = IT_0009*IT_0010*IT_0028;
    const ccomplex_t IT_0030 = IT_0011*IT_0012*IT_0028;
    const ccomplex_t IT_0031 = IT_0003*IT_0006*IT_0009*IT_0010;
    const ccomplex_t IT_0032 = IT_0003*IT_0006*IT_0011*IT_0012;
    const ccomplex_t IT_0033 = cpow(IT_0002, -2);
    const ccomplex_t IT_0034 = cpow(IT_0006 + IT_0033, -1);
    const ccomplex_t IT_0035 = (IT_0026 + -IT_0027 + IT_0029 + -IT_0030 + 2
      *IT_0031 + (-2)*IT_0032)*IT_0034;
    const ccomplex_t IT_0036 = IT_0024*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0038 + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0004;
    const ccomplex_t IT_0043 = -IT_0005;
    const ccomplex_t IT_0044 = IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0015 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0005*IT_0048;
    const ccomplex_t IT_0050 = IT_0019 + -IT_0045 + IT_0049;
    const ccomplex_t IT_0051 = s_24*s_34*IT_0021;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = s_23 + IT_0052;
    const ccomplex_t IT_0054 = s_14*s_34*IT_0021;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = s_13 + IT_0055;
    const ccomplex_t IT_0057 = IT_0041*IT_0042;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -conj(IT_0058);
    const ccomplex_t IT_0060 = pow(s_14, 2);
    const ccomplex_t IT_0061 = IT_0021*IT_0060;
    const ccomplex_t IT_0062 = s_14*s_24*IT_0021;
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = s_12 + IT_0063;
    const ccomplex_t IT_0065 = -conj(IT_0020);
    const ccomplex_t IT_0066 = pow(s_24, 2);
    const ccomplex_t IT_0067 = IT_0021*IT_0066;
    const ccomplex_t IT_0068 = -IT_0020*(conj(IT_0020)*(IT_0015 + -IT_0023) +
       conj(IT_0050)*IT_0053 + -IT_0056*IT_0059) + -IT_0058*(conj(IT_0058)*
      (IT_0038 + -IT_0061) + conj(IT_0050)*IT_0064 + -IT_0056*IT_0065) + IT_0050
      *(IT_0059*IT_0064 + IT_0053*IT_0065 + -conj(IT_0050)*(IT_0038 + -IT_0067));
    return create_ccomplex_return(IT_0068);
}

