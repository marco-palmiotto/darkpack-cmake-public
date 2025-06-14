#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_e_to_anti_e_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_e_to_anti_e_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*(IT_0006 + -IT_0009);
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_e, 2);
    const ccomplex_t IT_0016 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0023 = IT_0005*IT_0022;
    const ccomplex_t IT_0024 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*(IT_0023 + -IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0019 + IT_0033;
    const ccomplex_t IT_0035 = m_C_1*IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, -1);
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0001*IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = V_u1*U_Wm1;
    const ccomplex_t IT_0042 = IT_0005*IT_0041;
    const ccomplex_t IT_0043 = U_d1*V_Wp1;
    const ccomplex_t IT_0044 = IT_0008*IT_0043;
    const ccomplex_t IT_0045 = IT_0001*(IT_0042 + -IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0040*IT_0048;
    const ccomplex_t IT_0050 = IT_0017*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0001*IT_0036*IT_0038;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = V_u1*U_Wm2;
    const ccomplex_t IT_0055 = IT_0005*IT_0054;
    const ccomplex_t IT_0056 = U_d2*V_Wp1;
    const ccomplex_t IT_0057 = IT_0008*IT_0056;
    const ccomplex_t IT_0058 = IT_0001*(IT_0055 + -IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = IT_0031*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0051 + IT_0064;
    const ccomplex_t IT_0066 = m_C_1*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = pow(m_C_1, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0015 + IT_0068 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0070 = cos(theta_W);
    const ccomplex_t IT_0071 = cpow(IT_0070, -2);
    const ccomplex_t IT_0072 = IT_0000*IT_0071;
    const ccomplex_t IT_0073 = IT_0001 + IT_0072;
    const ccomplex_t IT_0074 = sin(beta);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_W*e_em*IT_0073*(IT_0008
      *IT_0037 + IT_0005*IT_0074);
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0040*IT_0076;
    const ccomplex_t IT_0078 = IT_0069*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0013*IT_0040;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0082 = IT_0017*IT_0080*IT_0081;
    const ccomplex_t IT_0083 = IT_0029*IT_0053;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0085 = IT_0031*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_e*e_em*IT_0001*IT_0008
      *IT_0036*IT_0038;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0040*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + IT_0015 + -IT_0016 + -IT_0068 
      + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = m_e*IT_0091;
    const ccomplex_t IT_0093 = -IT_0079 + IT_0082 + IT_0085 + 2*IT_0092;
    const ccomplex_t IT_0094 = m_e*m_C_1;
    const ccomplex_t IT_0095 = IT_0003*IT_0076;
    const ccomplex_t IT_0096 = IT_0069*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0021*IT_0061;
    const ccomplex_t IT_0099 = IT_0031*IT_0084*IT_0098;
    const ccomplex_t IT_0100 = IT_0003*IT_0087;
    const ccomplex_t IT_0101 = IT_0089*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = m_e*IT_0102;
    const ccomplex_t IT_0104 = IT_0003*IT_0048;
    const ccomplex_t IT_0105 = IT_0017*IT_0081*IT_0104;
    const ccomplex_t IT_0106 = IT_0097 + -IT_0099 + (-2)*IT_0103 + -IT_0105;
    const ccomplex_t IT_0107 = pow(m_h, 2);
    const ccomplex_t IT_0108 = m_e*m_C_1*IT_0107;
    const ccomplex_t IT_0109 = conj(IT_0065) + conj(IT_0102);
    const ccomplex_t IT_0110 = -IT_0019 + -IT_0033;
    const ccomplex_t IT_0111 = -IT_0091;
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = IT_0065 + IT_0102;
    const ccomplex_t IT_0114 = conj(IT_0110) + conj(IT_0111);
    const ccomplex_t IT_0115 = m_e*s_14;
    const ccomplex_t IT_0116 = (-2)*IT_0115;
    const ccomplex_t IT_0117 = conj(IT_0093)*IT_0116;
    const ccomplex_t IT_0118 = s_34*m_C_1;
    const ccomplex_t IT_0119 = (-2)*IT_0118;
    const ccomplex_t IT_0120 = conj(IT_0106)*IT_0119;
    const ccomplex_t IT_0121 = IT_0117 + IT_0120;
    const ccomplex_t IT_0122 = (-2)*s_13;
    const ccomplex_t IT_0123 = s_13*IT_0107;
    const ccomplex_t IT_0124 = s_14*s_34;
    const ccomplex_t IT_0125 = (-2)*IT_0124;
    const ccomplex_t IT_0126 = IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = (-2)*IT_0102;
    const ccomplex_t IT_0128 = (-2)*conj(IT_0102);
    const ccomplex_t IT_0129 = (-2)*IT_0094;
    const ccomplex_t IT_0130 = (-2)*IT_0126;
    const ccomplex_t IT_0131 = IT_0093*IT_0116;
    const ccomplex_t IT_0132 = IT_0110*IT_0130;
    const ccomplex_t IT_0133 = 2*s_13*(IT_0035*conj(IT_0035) + IT_0067*conj
      (IT_0067) + IT_0093*conj(IT_0093)) + 2*IT_0094*(conj(IT_0035)*(IT_0067 + 
      -IT_0093) + IT_0035*(conj(IT_0067) + -conj(IT_0093)) + IT_0093*conj
      (IT_0106)) + 2*IT_0108*(IT_0109*IT_0112 + IT_0113*IT_0114) + conj(IT_0065)
      *IT_0093*IT_0119 + IT_0110*IT_0121 + (conj(IT_0067)*IT_0093 + IT_0067*conj
      (IT_0093) + IT_0035*conj(IT_0106))*IT_0122 + IT_0109*IT_0126*IT_0127 + 2
      *IT_0118*(IT_0067*IT_0109 + conj(IT_0035)*IT_0112 + conj(IT_0067)*IT_0113 
      + IT_0035*IT_0114 + 0.5*conj(IT_0093)*IT_0127 + 0.5*IT_0093*IT_0128) + 2
      *IT_0115*(IT_0035*IT_0109 + conj(IT_0067)*IT_0112 + conj(IT_0035)*IT_0113 
      + IT_0067*IT_0114 + 0.5*conj(IT_0106)*IT_0127 + 0.5*IT_0106*IT_0128) +
       IT_0067*conj(IT_0106)*IT_0129 + 2*IT_0106*(conj(IT_0093)*IT_0094 + s_13
      *conj(IT_0106) + 0.5*conj(IT_0065)*IT_0116 + 0.5*IT_0114*IT_0119 + 0.5
      *conj(IT_0035)*IT_0122 + 0.5*conj(IT_0067)*IT_0129) + IT_0065*(conj
      (IT_0106)*IT_0116 + conj(IT_0093)*IT_0119 + IT_0126*IT_0128 + conj(IT_0065
      )*IT_0130) + IT_0111*(IT_0121 + IT_0114*IT_0130) + IT_0114*(IT_0131 +
       IT_0132);
    return create_ccomplex_return(IT_0133);
}

