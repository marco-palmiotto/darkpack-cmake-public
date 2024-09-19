#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_d_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_d_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = s_13*s_14;
    const ccomplex_t IT_0003 = (-2)*IT_0002;
    const ccomplex_t IT_0004 = IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = pow(m_d, 2);
    const ccomplex_t IT_0006 = pow(m_su_L, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + -IT_0000 + IT_0005 + -IT_0006 
      + -reg_prop, -1);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0013 = (-0.333333333333333)*IT_0012;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (-1.33333333333333)*IT_0012;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0000 + IT_0005 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = 6*IT_0022;
    const ccomplex_t IT_0024 = IT_0011*IT_0012;
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (-2)*IT_0027;
    const ccomplex_t IT_0029 = (-6)*IT_0028;
    const ccomplex_t IT_0030 = IT_0023 + IT_0029;
    const ccomplex_t IT_0031 = 6*conj(IT_0022);
    const ccomplex_t IT_0032 = (-6)*conj(IT_0028);
    const ccomplex_t IT_0033 = IT_0031 + IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = cos(beta);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0009*IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0013*IT_0038;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0018*IT_0038;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0017*IT_0044;
    const ccomplex_t IT_0046 = 6*IT_0045;
    const ccomplex_t IT_0047 = IT_0012*IT_0038;
    const ccomplex_t IT_0048 = IT_0025*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = 2*IT_0049;
    const ccomplex_t IT_0051 = (-6)*IT_0050;
    const ccomplex_t IT_0052 = IT_0046 + IT_0051;
    const ccomplex_t IT_0053 = 6*conj(IT_0045);
    const ccomplex_t IT_0054 = (-6)*conj(IT_0050);
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = -IT_0049;
    const ccomplex_t IT_0057 = m_d*s_14*m_C_2;
    const ccomplex_t IT_0058 = -(conj(IT_0022) + -conj(IT_0028))*IT_0042 + -
      (IT_0022 + -IT_0028)*conj(IT_0042) + -conj(IT_0016)*(IT_0045 + -IT_0050) +
       -IT_0016*(conj(IT_0045) + -conj(IT_0050)) + 0.166666666666667*conj
      (IT_0027)*IT_0052 + 0.166666666666667*IT_0027*IT_0055 + 0.166666666666667
      *IT_0033*IT_0056 + 0.166666666666667*IT_0030*conj(IT_0056);
    const ccomplex_t IT_0059 = s_14*s_34;
    const ccomplex_t IT_0060 = 24*IT_0016;
    const ccomplex_t IT_0061 = (-0.666666666666667)*IT_0012;
    const ccomplex_t IT_0062 = -IT_0020*(IT_0017 + -IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = 6*s_13;
    const ccomplex_t IT_0065 = s_14*IT_0064;
    const ccomplex_t IT_0066 = (-6)*s_13;
    const ccomplex_t IT_0067 = s_14*IT_0066;
    const ccomplex_t IT_0068 = 2*IT_0016 + IT_0021;
    const ccomplex_t IT_0069 = s_34*IT_0066;
    const ccomplex_t IT_0070 = -IT_0044*(IT_0017 + -IT_0061);
    const ccomplex_t IT_0071 = (-2)*IT_0041 + -IT_0045;
    const ccomplex_t IT_0072 = IT_0000*IT_0064;
    const ccomplex_t IT_0073 = pow(s_13, 2);
    const ccomplex_t IT_0074 = 6*IT_0073;
    const ccomplex_t IT_0075 = IT_0000*IT_0066;
    const ccomplex_t IT_0076 = (-6)*IT_0005;
    const ccomplex_t IT_0077 = (-6)*IT_0073;
    const ccomplex_t IT_0078 = m_d*m_C_2;
    const ccomplex_t IT_0079 = 6*s_14;
    const ccomplex_t IT_0080 = (-6)*s_14;
    const ccomplex_t IT_0081 = (-6)*s_34;
    const ccomplex_t IT_0082 = 6*IT_0000;
    const ccomplex_t IT_0083 = (-6)*IT_0000;
    const ccomplex_t IT_0084 = conj(IT_0068)*(IT_0050*IT_0064 + IT_0045
      *IT_0066 + IT_0071*IT_0076 + IT_0070*IT_0081) + IT_0068*(conj(IT_0050)
      *IT_0064 + conj(IT_0045)*IT_0066 + conj(IT_0071)*IT_0076 + conj(IT_0070)
      *IT_0081) + conj(IT_0063)*(IT_0050*IT_0079 + IT_0045*IT_0080 + IT_0071
      *IT_0081) + IT_0063*(conj(IT_0050)*IT_0079 + conj(IT_0045)*IT_0080 + conj
      (IT_0071)*IT_0081) + conj(IT_0022)*(IT_0066*IT_0071 + IT_0070*IT_0080 +
       IT_0050*IT_0082 + IT_0045*IT_0083) + IT_0022*(IT_0066*conj(IT_0071) +
       conj(IT_0070)*IT_0080 + conj(IT_0050)*IT_0082 + conj(IT_0045)*IT_0083) +
       conj(IT_0028)*(IT_0064*IT_0071 + IT_0070*IT_0079 + IT_0045*IT_0082 +
       IT_0050*IT_0083) + IT_0028*(IT_0064*conj(IT_0071) + conj(IT_0070)*IT_0079
       + conj(IT_0045)*IT_0082 + conj(IT_0050)*IT_0083);
    const ccomplex_t IT_0085 = s_14*IT_0005;
    const ccomplex_t IT_0086 = (-6)*IT_0068;
    const ccomplex_t IT_0087 = (-6)*conj(IT_0068);
    const ccomplex_t IT_0088 = (-6)*IT_0071;
    const ccomplex_t IT_0089 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0090 = s_13*s_34;
    const ccomplex_t IT_0091 = (-2)*IT_0090;
    const ccomplex_t IT_0092 = IT_0085 + IT_0091;
    const ccomplex_t IT_0093 = m_d*s_34*m_C_2;
    const ccomplex_t IT_0094 = 6*(conj(IT_0056)*IT_0068 + IT_0056*conj(IT_0068
      ) + conj(IT_0027)*IT_0071 + IT_0027*conj(IT_0071) + 0.166666666666667*conj
      (IT_0042)*IT_0086 + 0.166666666666667*IT_0042*IT_0087 + 0.166666666666667
      *conj(IT_0016)*IT_0088 + 0.166666666666667*IT_0016*IT_0089)*IT_0093;
    const ccomplex_t IT_0095 = IT_0004*(conj(IT_0016)*IT_0030 + IT_0016
      *IT_0033 + conj(IT_0042)*IT_0052 + IT_0042*IT_0055) + IT_0001*(conj
      (IT_0027)*IT_0030 + IT_0027*IT_0033 + IT_0055*IT_0056 + IT_0052*conj
      (IT_0056)) + 6*IT_0057*IT_0058 + IT_0059*((IT_0042 + IT_0056)*(24*conj
      (IT_0042) + 24*conj(IT_0056)) + (conj(IT_0016) + conj(IT_0027))*(24
      *IT_0027 + IT_0060)) + conj(IT_0063)*(IT_0028*IT_0065 + IT_0022*IT_0067 +
       IT_0068*IT_0069) + conj(IT_0070)*(IT_0050*IT_0065 + IT_0045*IT_0067 +
       IT_0069*IT_0071) + conj(IT_0028)*(IT_0063*IT_0065 + IT_0022*IT_0072 +
       IT_0068*IT_0074 + IT_0028*IT_0075) + conj(IT_0050)*(IT_0065*IT_0070 +
       IT_0045*IT_0072 + IT_0071*IT_0074 + IT_0050*IT_0075) + conj(IT_0068)*
      (IT_0063*IT_0069 + IT_0028*IT_0074 + s_13*IT_0068*IT_0076 + IT_0022
      *IT_0077) + conj(IT_0071)*(IT_0069*IT_0070 + IT_0050*IT_0074 + s_13
      *IT_0071*IT_0076 + IT_0045*IT_0077) + conj(IT_0022)*(IT_0063*IT_0067 +
       IT_0028*IT_0072 + IT_0022*IT_0075 + IT_0068*IT_0077) + conj(IT_0045)*
      (IT_0067*IT_0070 + IT_0050*IT_0072 + IT_0045*IT_0075 + IT_0071*IT_0077) +
       IT_0078*IT_0084 + IT_0085*(conj(IT_0016)*IT_0086 + IT_0016*IT_0087 + conj
      (IT_0042)*IT_0088 + IT_0042*IT_0089) + (conj(IT_0027)*IT_0086 + IT_0027
      *IT_0087 + conj(IT_0056)*IT_0088 + IT_0056*IT_0089)*IT_0092 + IT_0094;
    return create_ccomplex_return(IT_0095);
}

