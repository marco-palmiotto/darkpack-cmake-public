#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_se_L_to_nu_e_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_se_L_to_nu_e_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = U_d1*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = V_u1*U_Wm1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*(IT_0005 + IT_0008);
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = pow(m_se_L, 2);
    const ccomplex_t IT_0015 = cpow(s_23 + (-0.5)*IT_0014 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0019 = U_d2*V_Wp1;
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = V_u1*U_Wm2;
    const ccomplex_t IT_0022 = IT_0006*IT_0021;
    const ccomplex_t IT_0023 = IT_0001*(IT_0020 + IT_0022);
    const ccomplex_t IT_0024 = 1.4142135623731*e_em*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0026;
    const ccomplex_t IT_0028 = cpow(s_23 + (-0.5)*IT_0014 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0017 + 0.5*IT_0030;
    const ccomplex_t IT_0032 = pow(m_H0, 2);
    const ccomplex_t IT_0033 = s_13*IT_0032;
    const ccomplex_t IT_0034 = s_14*s_34;
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0001
      *IT_0006*IT_0035*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0001*IT_0035*IT_0037;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = pow(m_C_1, 2);
    const ccomplex_t IT_0044 = cpow(s_13 + (-0.5)*IT_0043 + 0.5*m_se_R*(m_se_R
       + (0 + _Complex_I*-1)*Gamma_er) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = pow(m_e, 2);
    const ccomplex_t IT_0048 = sin(beta);
    const ccomplex_t IT_0049 = cos(theta_W);
    const ccomplex_t IT_0050 = cpow(IT_0049, -2);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0003
      *IT_0035*IT_0037*IT_0047 + (-0.5)*m_W*(IT_0003*IT_0036 + -IT_0006*IT_0048)
      *(IT_0001 + -IT_0000*IT_0050));
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0002*IT_0052;
    const ccomplex_t IT_0054 = cpow(s_13 + (-0.5)*IT_0043 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0060 = IT_0006*IT_0059;
    const ccomplex_t IT_0061 = IT_0001*(IT_0058 + IT_0060);
    const ccomplex_t IT_0062 = 1.4142135623731*e_em*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0002*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0067 = IT_0015*IT_0065*IT_0066;
    const ccomplex_t IT_0068 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0071 = IT_0006*IT_0070;
    const ccomplex_t IT_0072 = IT_0001*(IT_0069 + IT_0071);
    const ccomplex_t IT_0073 = 1.4142135623731*e_em*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0018*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0078 = IT_0028*IT_0076*IT_0077;
    const ccomplex_t IT_0079 = m_C_1*IT_0017;
    const ccomplex_t IT_0080 = m_C_1*IT_0030;
    const ccomplex_t IT_0081 = 0.5*IT_0046 + 0.5*IT_0056 + (-0.5)*IT_0067 + (
      -0.5)*IT_0078 + (-0.5)*IT_0079 + (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = s_34*m_C_1;
    const ccomplex_t IT_0083 = 2*IT_0082;
    return create_ccomplex_return(2*IT_0081*(s_13*conj(IT_0081) + 0.5*conj
      (IT_0031)*IT_0083) + (-2)*IT_0031*(conj(IT_0031)*(IT_0033 + (-2)*IT_0034) 
      + (-0.5)*conj(IT_0081)*IT_0083));
}

