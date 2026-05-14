#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_2_to_b_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_2_to_b_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_b = param->m_b;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0003 = pow(m_b, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = pow(m_sb_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0009 
      + -reg_prop, -1);
    const ccomplex_t IT_0011 = (-0.333333333333333)*IT_0000;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0008 + 2*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0016;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = -IT_0019 + (-2)*IT_0022;
    const ccomplex_t IT_0024 = (-8)*IT_0003;
    const ccomplex_t IT_0025 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0025;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = (-8)*s_13;
    const ccomplex_t IT_0030 = s_34*IT_0029;
    const ccomplex_t IT_0031 = -IT_0008;
    const ccomplex_t IT_0032 = -IT_0022;
    const ccomplex_t IT_0033 = m_b*m_sG*s_14;
    const ccomplex_t IT_0034 = s_34*IT_0004;
    const ccomplex_t IT_0035 = 8*IT_0034;
    const ccomplex_t IT_0036 = s_13*s_14;
    const ccomplex_t IT_0037 = (-16)*IT_0036;
    const ccomplex_t IT_0038 = IT_0035 + IT_0037;
    const ccomplex_t IT_0039 = (-8)*IT_0004;
    const ccomplex_t IT_0040 = s_13*IT_0039;
    const ccomplex_t IT_0041 = s_14*IT_0029;
    const ccomplex_t IT_0042 = pow(s_13, 2);
    const ccomplex_t IT_0043 = (-8)*IT_0042;
    const ccomplex_t IT_0044 = m_b*m_sG;
    const ccomplex_t IT_0045 = (-8)*s_14;
    const ccomplex_t IT_0046 = (-8)*s_34;
    const ccomplex_t IT_0047 = (conj(IT_0015)*IT_0023 + IT_0015*conj(IT_0023))
      *IT_0024 + IT_0029*(conj(IT_0015)*IT_0019 + IT_0015*conj(IT_0019) + conj
      (IT_0023)*IT_0031 + IT_0023*conj(IT_0031)) + (conj(IT_0019)*IT_0031 +
       IT_0019*conj(IT_0031))*IT_0039 + (conj(IT_0019)*IT_0026 + IT_0019*conj
      (IT_0026) + conj(IT_0028)*IT_0031 + IT_0028*conj(IT_0031))*IT_0045 + (conj
      (IT_0023)*IT_0026 + IT_0023*conj(IT_0026) + conj(IT_0015)*IT_0028 +
       IT_0015*conj(IT_0028))*IT_0046;
    const ccomplex_t IT_0048 = m_b*m_sG*s_34;
    const ccomplex_t IT_0049 = (-8)*conj(IT_0015);
    const ccomplex_t IT_0050 = s_14*IT_0003;
    const ccomplex_t IT_0051 = (-8)*conj(IT_0023);
    const ccomplex_t IT_0052 = s_14*s_34;
    const ccomplex_t IT_0053 = (-8)*IT_0023;
    const ccomplex_t IT_0054 = (-0.125)*IT_0050;
    const ccomplex_t IT_0055 = (-8)*IT_0015;
    const ccomplex_t IT_0056 = (-0.125)*IT_0048;
    const ccomplex_t IT_0057 = s_13*(IT_0015*conj(IT_0015) + IT_0023*conj
      (IT_0023))*IT_0024 + (conj(IT_0015)*IT_0026 + IT_0015*conj(IT_0026) + conj
      (IT_0023)*IT_0028 + IT_0023*conj(IT_0028))*IT_0030 + (-8)*conj(IT_0031)*
      (IT_0032*IT_0033 + (-0.125)*IT_0014*IT_0038 + (-0.125)*IT_0031*IT_0040 + (
      -0.125)*IT_0026*IT_0041 + (-0.125)*IT_0015*IT_0043) + IT_0031*(conj
      (IT_0026)*IT_0041 + conj(IT_0015)*IT_0043) + (-8)*conj(IT_0019)*(IT_0014
      *IT_0033 + (-0.125)*IT_0032*IT_0038 + (-0.125)*IT_0019*IT_0040 + (-0.125)
      *IT_0028*IT_0041 + (-0.125)*IT_0023*IT_0043) + IT_0019*(conj(IT_0028)
      *IT_0041 + conj(IT_0023)*IT_0043) + IT_0044*IT_0047 + IT_0014*(IT_0049
      *IT_0050 + IT_0048*IT_0051) + IT_0032*(IT_0048*IT_0049 + IT_0050*IT_0051) 
      + (-8)*conj(IT_0014)*(IT_0019*IT_0033 + (-0.125)*IT_0031*IT_0038 + (-4)
      *IT_0014*IT_0052 + IT_0054*IT_0055 + IT_0053*IT_0056) + (-8)*conj(IT_0032)
      *(IT_0031*IT_0033 + (-0.125)*IT_0019*IT_0038 + (-4)*IT_0032*IT_0052 +
       IT_0053*IT_0054 + IT_0055*IT_0056);
    return create_ccomplex_return(IT_0057);
}

