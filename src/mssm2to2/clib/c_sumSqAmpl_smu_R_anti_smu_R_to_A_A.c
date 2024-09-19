#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_A_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_A_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = tan(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cpow(IT_0001 + IT_0001*IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0009 = 2*IT_0008;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0009, 2);
    const ccomplex_t IT_0012 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_smu_R
      *Gamma_mr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (IT_0010 + -IT_0011)*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_smu_R
      *Gamma_mr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0008, 2);
    const ccomplex_t IT_0019 = (-0.5)*IT_0011 + (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0010 + IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = -IT_0018;
    const ccomplex_t IT_0023 = IT_0010 + IT_0022;
    const ccomplex_t IT_0024 = IT_0013*IT_0023;
    const ccomplex_t IT_0025 = IT_0021 + 0.5*IT_0024;
    const ccomplex_t IT_0026 = conj(IT_0015) + conj(IT_0025);
    const ccomplex_t IT_0027 = (IT_0010 + -IT_0011)*IT_0017;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0010*IT_0013;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = conj(IT_0028) + conj(IT_0030);
    const ccomplex_t IT_0032 = s_12*IT_0007;
    const ccomplex_t IT_0033 = (IT_0013 + IT_0017)*(IT_0010 + -IT_0018);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = pow(s_13, 2);
    const ccomplex_t IT_0036 = IT_0010*IT_0017;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = pow(s_23, 2);
    const ccomplex_t IT_0039 = 0.5*IT_0014;
    const ccomplex_t IT_0040 = s_13*s_23;
    const ccomplex_t IT_0041 = 0.5*IT_0036;
    const ccomplex_t IT_0042 = IT_0017*IT_0023;
    const ccomplex_t IT_0043 = IT_0013*IT_0018;
    const ccomplex_t IT_0044 = 0.5*IT_0042 + (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = conj(IT_0041) + conj(IT_0044);
    const ccomplex_t IT_0046 = s_13*IT_0007;
    const ccomplex_t IT_0047 = pow(m_smu_R, 2);
    const ccomplex_t IT_0048 = s_23*IT_0047;
    const ccomplex_t IT_0049 = IT_0039*IT_0048;
    const ccomplex_t IT_0050 = -s_23*IT_0007;
    const ccomplex_t IT_0051 = s_12*s_13;
    const ccomplex_t IT_0052 = conj(IT_0037)*IT_0051;
    const ccomplex_t IT_0053 = s_12*s_23;
    const ccomplex_t IT_0054 = s_13*IT_0047;
    const ccomplex_t IT_0055 = conj(IT_0039)*IT_0048;
    const ccomplex_t IT_0056 = IT_0007*IT_0047;
    const ccomplex_t IT_0057 = pow(m_smu_R, 4);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -conj(IT_0041);
    const ccomplex_t IT_0060 = s_12*IT_0047;
    const ccomplex_t IT_0061 = -conj(IT_0039)*IT_0060;
    const ccomplex_t IT_0062 = -conj(IT_0015);
    const ccomplex_t IT_0063 = -conj(IT_0044);
    const ccomplex_t IT_0064 = -conj(IT_0025);
    const ccomplex_t IT_0065 = -conj(IT_0034)*IT_0040;
    const ccomplex_t IT_0066 = conj(IT_0037)*IT_0060;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = pow(s_12, 2);
    const ccomplex_t IT_0069 = -IT_0007*IT_0047;
    const ccomplex_t IT_0070 = IT_0031*IT_0060;
    const ccomplex_t IT_0071 = -s_13*IT_0007;
    const ccomplex_t IT_0072 = conj(IT_0037)*IT_0054;
    const ccomplex_t IT_0073 = conj(IT_0039)*IT_0053;
    const ccomplex_t IT_0074 = (-4)*IT_0007*(IT_0007 + (-0.25)*s_23*IT_0026) +
       IT_0031*IT_0032 + IT_0034*(IT_0035*conj(IT_0037) + IT_0038*conj(IT_0039) 
      + IT_0031*IT_0040) + IT_0045*IT_0046 + IT_0026*IT_0049 + IT_0025*(conj
      (IT_0015)*IT_0038 + IT_0040*conj(IT_0041) + IT_0050 + IT_0052 + conj
      (IT_0030)*IT_0053 + conj(IT_0028)*IT_0054 + IT_0055) + IT_0015*(conj
      (IT_0025)*IT_0038 + IT_0040*conj(IT_0044) + IT_0050 + IT_0052 + conj
      (IT_0028)*IT_0053 + conj(IT_0030)*IT_0054 + IT_0055) + (conj(IT_0037) +
       conj(IT_0039))*IT_0056 + -IT_0030*(s_12*(IT_0007 + -s_12*conj(IT_0028)) +
       conj(IT_0030)*IT_0058 + IT_0051*IT_0059 + IT_0061 + IT_0054*IT_0062 +
       IT_0048*IT_0063 + IT_0053*IT_0064 + IT_0065 + IT_0067) + -IT_0028*(s_12*
      (IT_0007 + -s_12*conj(IT_0030)) + conj(IT_0028)*IT_0058 + IT_0048*IT_0059 
      + IT_0061 + IT_0053*IT_0062 + IT_0051*IT_0063 + IT_0054*IT_0064 + IT_0065 
      + IT_0067) + IT_0039*(conj(IT_0034)*IT_0038 + IT_0045*IT_0053 + conj
      (IT_0039)*IT_0057 + conj(IT_0037)*IT_0068 + IT_0069 + IT_0070) + IT_0037*
      (conj(IT_0034)*IT_0035 + IT_0026*IT_0051 + IT_0045*IT_0054 + conj(IT_0037)
      *IT_0057 + conj(IT_0039)*IT_0068 + IT_0069 + IT_0070) + IT_0044*(conj
      (IT_0015)*IT_0040 + IT_0035*conj(IT_0041) + conj(IT_0030)*IT_0048 + conj
      (IT_0028)*IT_0051 + IT_0071 + IT_0072 + IT_0073) + IT_0041*(conj(IT_0025)
      *IT_0040 + IT_0035*conj(IT_0044) + conj(IT_0028)*IT_0048 + conj(IT_0030)
      *IT_0051 + IT_0071 + IT_0072 + IT_0073);
    return create_ccomplex_return(IT_0074);
}

