#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_R_to_anti_d_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_R_to_anti_d_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0000*IT_0005*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = pow(m_d, 2);
    const ccomplex_t IT_0014 = pow(m_C_1, 2);
    const ccomplex_t IT_0015 = pow(m_su_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + IT_0013 + -IT_0014 + -IT_0015 
      + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = m_d*s_14;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0000*IT_0005*IT_0009;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = IT_0005*(IT_0024 + -IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0021*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0013 + IT_0015 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0036 = IT_0022*IT_0035;
    const ccomplex_t IT_0037 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*(IT_0036 + -IT_0038);
    const ccomplex_t IT_0040 = 1.4142135623731*e_em*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0011*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0013 + IT_0015 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0034 + IT_0046;
    const ccomplex_t IT_0048 = m_C_1*IT_0047;
    const ccomplex_t IT_0049 = pow(m_h, 2);
    const ccomplex_t IT_0050 = m_d*m_C_1*IT_0049;
    const ccomplex_t IT_0051 = -IT_0034 + -IT_0046;
    const ccomplex_t IT_0052 = m_d*IT_0018;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0003*IT_0005*IT_0009;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = V_u1*U_Wm1;
    const ccomplex_t IT_0062 = IT_0022*IT_0061;
    const ccomplex_t IT_0063 = U_d1*V_Wp1;
    const ccomplex_t IT_0064 = IT_0003*IT_0063;
    const ccomplex_t IT_0065 = IT_0005*(IT_0062 + -IT_0064);
    const ccomplex_t IT_0066 = 1.4142135623731*e_em*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0011*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0071 = IT_0044*IT_0069*IT_0070;
    const ccomplex_t IT_0072 = cos(theta_W);
    const ccomplex_t IT_0073 = cpow(IT_0072, -2);
    const ccomplex_t IT_0074 = pow(m_u, 2);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*(IT_0001
      *IT_0003 + IT_0008*IT_0022)*IT_0073 + (-1.5)*IT_0000*IT_0005*IT_0009
      *IT_0022*IT_0074);
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = IT_0011*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = V_u1*U_Wm2;
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = U_d2*V_Wp1;
    const ccomplex_t IT_0084 = IT_0003*IT_0083;
    const ccomplex_t IT_0085 = IT_0005*(IT_0082 + -IT_0084);
    const ccomplex_t IT_0086 = 1.4142135623731*e_em*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0021*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0091 = IT_0032*IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (-2)*IT_0052 + IT_0060 + -IT_0071 + IT_0080 + 
      -IT_0091;
    const ccomplex_t IT_0093 = s_34*m_C_1;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0054*IT_0095;
    const ccomplex_t IT_0097 = IT_0058*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = s_13*IT_0049;
    const ccomplex_t IT_0101 = s_14*s_34;
    const ccomplex_t IT_0102 = (-2)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = (-6)*IT_0018;
    const ccomplex_t IT_0105 = (-6)*IT_0093;
    const ccomplex_t IT_0106 = (-6)*IT_0019;
    const ccomplex_t IT_0107 = m_d*m_C_1;
    const ccomplex_t IT_0108 = -conj(IT_0092);
    const ccomplex_t IT_0109 = conj(IT_0048) + IT_0108;
    const ccomplex_t IT_0110 = (-6)*conj(IT_0018);
    const ccomplex_t IT_0111 = (-0.166666666666667)*conj(IT_0051);
    const ccomplex_t IT_0112 = s_13*IT_0109;
    const ccomplex_t IT_0113 = 6*IT_0018*(IT_0019*conj(IT_0048) + IT_0050*conj
      (IT_0051)) + (IT_0019*conj(IT_0092) + IT_0093*conj(IT_0099) + conj(IT_0018
      )*IT_0103)*IT_0104 + 6*IT_0051*(conj(IT_0018)*IT_0050 + conj(IT_0048)
      *IT_0093 + -conj(IT_0051)*IT_0103 + 0.166666666666667*conj(IT_0092)
      *IT_0105 + 0.166666666666667*conj(IT_0099)*IT_0106) + 6*IT_0099*(s_13*conj
      (IT_0099) + -IT_0107*IT_0109 + 0.166666666666667*IT_0093*IT_0110 + 
      -IT_0106*IT_0111) + 6*IT_0048*(conj(IT_0018)*IT_0019 + conj(IT_0051)
      *IT_0093 + -conj(IT_0099)*IT_0107 + IT_0112) + 6*IT_0092*(conj(IT_0099)
      *IT_0107 + 0.166666666666667*IT_0019*IT_0110 + -IT_0105*IT_0111 + -IT_0112);
    return create_ccomplex_return(IT_0113);
}

