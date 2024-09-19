#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_A_Z(
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
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -3);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = tan(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = cpow(1 + IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0003
      *IT_0006;
    const ccomplex_t IT_0008 = (-2)*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = pow(m_Z, -2);
    const ccomplex_t IT_0011 = pow(s_34, 2);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_smu_R
      *Gamma_mr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0016 = cpow(IT_0001, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (-2)*IT_0017;
    const ccomplex_t IT_0021 = IT_0015*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0019 + IT_0022;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_smu_R
      *Gamma_mr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019;
    const ccomplex_t IT_0028 = 2*IT_0015;
    const ccomplex_t IT_0029 = IT_0018*IT_0028;
    const ccomplex_t IT_0030 = IT_0027 + IT_0029;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0024 + (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0029;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = -s_13;
    const ccomplex_t IT_0036 = IT_0020*IT_0028;
    const ccomplex_t IT_0037 = IT_0014*(IT_0029 + -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = -s_23;
    const ccomplex_t IT_0040 = -IT_0008;
    const ccomplex_t IT_0041 = pow(m_smu_R, 2);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0032*(IT_0008 + conj(IT_0034)*IT_0035 + conj
      (IT_0038)*IT_0039) + conj(IT_0032)*IT_0040 + IT_0034*(conj(IT_0032)
      *IT_0035 + -s_12*conj(IT_0038) + conj(IT_0034)*IT_0042) + -IT_0038*(s_12
      *conj(IT_0034) + -conj(IT_0032)*IT_0039 + -conj(IT_0038)*IT_0042);
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0041 + IT_0046;
    const ccomplex_t IT_0048 = IT_0014*IT_0019;
    const ccomplex_t IT_0049 = 0.5*IT_0031 + (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0033;
    const ccomplex_t IT_0051 = s_13*IT_0050;
    const ccomplex_t IT_0052 = IT_0014*IT_0029;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = s_23*IT_0053;
    const ccomplex_t IT_0055 = IT_0051 + IT_0054;
    const ccomplex_t IT_0056 = s_13*conj(IT_0050);
    const ccomplex_t IT_0057 = s_23*conj(IT_0053);
    const ccomplex_t IT_0058 = IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = s_12*conj(IT_0050);
    const ccomplex_t IT_0060 = IT_0041*IT_0050;
    const ccomplex_t IT_0061 = s_12*IT_0050;
    const ccomplex_t IT_0062 = IT_0041*IT_0053;
    const ccomplex_t IT_0063 = s_14*s_24*IT_0010;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = s_12 + IT_0064;
    const ccomplex_t IT_0066 = IT_0026*(IT_0029 + -IT_0036);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = s_12*conj(IT_0053);
    const ccomplex_t IT_0069 = s_13*conj(IT_0049);
    const ccomplex_t IT_0070 = IT_0041*conj(IT_0050);
    const ccomplex_t IT_0071 = IT_0040 + IT_0068 + IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0037;
    const ccomplex_t IT_0073 = s_23*conj(IT_0049);
    const ccomplex_t IT_0074 = IT_0041*conj(IT_0053);
    const ccomplex_t IT_0075 = IT_0059 + IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = IT_0021 + IT_0029;
    const ccomplex_t IT_0077 = -IT_0019 + -IT_0036;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = IT_0026*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0024 + 0.5*IT_0079;
    const ccomplex_t IT_0081 = s_12*IT_0053;
    const ccomplex_t IT_0082 = s_13*IT_0049;
    const ccomplex_t IT_0083 = s_23*IT_0049;
    const ccomplex_t IT_0084 = IT_0040*IT_0053 + IT_0008*conj(IT_0053) +
       IT_0067*IT_0071 + IT_0072*IT_0075 + IT_0058*IT_0080 + IT_0055*conj
      (IT_0080) + conj(IT_0067)*(IT_0008 + IT_0060 + IT_0081 + IT_0082) + conj
      (IT_0072)*(IT_0061 + IT_0062 + IT_0083);
    const ccomplex_t IT_0085 = s_14*s_34*IT_0010;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = s_13 + IT_0086;
    const ccomplex_t IT_0088 = IT_0040*IT_0049 + IT_0008*conj(IT_0049) + conj
      (IT_0032)*IT_0055 + IT_0032*IT_0058 + IT_0034*IT_0071 + IT_0038*IT_0075 +
       conj(IT_0034)*(IT_0008 + IT_0060 + IT_0081 + IT_0082) + conj(IT_0038)*
      (IT_0061 + IT_0062 + IT_0083);
    const ccomplex_t IT_0089 = pow(s_24, 2);
    const ccomplex_t IT_0090 = IT_0010*IT_0089;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0041 + IT_0091;
    const ccomplex_t IT_0093 = s_13*IT_0067;
    const ccomplex_t IT_0094 = s_23*IT_0072;
    const ccomplex_t IT_0095 = IT_0093 + IT_0094;
    const ccomplex_t IT_0096 = s_13*conj(IT_0067);
    const ccomplex_t IT_0097 = s_23*conj(IT_0072);
    const ccomplex_t IT_0098 = IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = s_12*conj(IT_0067);
    const ccomplex_t IT_0100 = IT_0041*IT_0067;
    const ccomplex_t IT_0101 = s_12*IT_0067;
    const ccomplex_t IT_0102 = IT_0041*IT_0072;
    const ccomplex_t IT_0103 = s_24*s_34*IT_0010;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = s_23 + IT_0104;
    const ccomplex_t IT_0106 = IT_0040*IT_0080 + IT_0008*conj(IT_0080) +
       IT_0034*(IT_0041*conj(IT_0067) + s_12*conj(IT_0072) + s_13*conj(IT_0080))
       + conj(IT_0032)*IT_0095 + IT_0032*IT_0098 + IT_0038*(IT_0040 + IT_0041
      *conj(IT_0072) + s_23*conj(IT_0080) + IT_0099) + conj(IT_0034)*(s_12
      *IT_0072 + s_13*IT_0080 + IT_0100) + conj(IT_0038)*(IT_0008 + s_23*IT_0080
       + IT_0101 + IT_0102);
    const ccomplex_t IT_0107 = (-3)*IT_0009 + IT_0012*IT_0043 + IT_0047*
      (IT_0040*IT_0050 + conj(IT_0049)*IT_0055 + IT_0049*IT_0058 + IT_0053
      *IT_0059 + conj(IT_0050)*(IT_0008 + IT_0060) + conj(IT_0053)*(IT_0061 +
       IT_0062)) + IT_0065*IT_0084 + IT_0087*IT_0088 + IT_0092*(conj(IT_0080)
      *IT_0095 + IT_0080*IT_0098 + IT_0072*(IT_0040 + IT_0099) + conj(IT_0067)
      *IT_0100 + conj(IT_0072)*(IT_0008 + IT_0101 + IT_0102)) + IT_0105*IT_0106;
    return create_ccomplex_return(IT_0107);
}

