#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_sb_1_to_d_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_sb_1_to_d_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_d, 2);
    const ccomplex_t IT_0001 = pow(m_su_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0007;
    const ccomplex_t IT_0011 = m_b*conj(U_d2)*V_tb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0013 + (-2)*IT_0015);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0009*IT_0017;
    const ccomplex_t IT_0019 = IT_0002*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = m_b*conj(U_d1)*V_tb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = IT_0023*IT_0029;
    const ccomplex_t IT_0031 = IT_0021*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0020 + IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0007;
    const ccomplex_t IT_0035 = sin(beta);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0007*IT_0036*U_sb_00;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0034*IT_0038;
    const ccomplex_t IT_0040 = IT_0002*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0007;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0007*IT_0036*U_sb_00;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = IT_0021*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0041 + IT_0047;
    const ccomplex_t IT_0049 = s_34*IT_0001;
    const ccomplex_t IT_0050 = (-18)*IT_0049;
    const ccomplex_t IT_0051 = s_13*s_14;
    const ccomplex_t IT_0052 = 36*IT_0051;
    const ccomplex_t IT_0053 = IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = m_d*IT_0047;
    const ccomplex_t IT_0055 = m_d*IT_0041;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0057 = IT_0009*IT_0038;
    const ccomplex_t IT_0058 = IT_0002*IT_0056*IT_0057;
    const ccomplex_t IT_0059 = IT_0054 + IT_0055 + -IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0061 = IT_0029*IT_0042;
    const ccomplex_t IT_0062 = IT_0021*IT_0060*IT_0061;
    const ccomplex_t IT_0063 = IT_0017*IT_0034;
    const ccomplex_t IT_0064 = IT_0002*IT_0056*IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = IT_0023*IT_0044;
    const ccomplex_t IT_0067 = IT_0021*IT_0060*IT_0066;
    const ccomplex_t IT_0068 = m_d*IT_0033;
    const ccomplex_t IT_0069 = 18*s_34;
    const ccomplex_t IT_0070 = (-18)*s_34;
    const ccomplex_t IT_0071 = m_d*m_t*IT_0001;
    const ccomplex_t IT_0072 = (-18)*IT_0033;
    const ccomplex_t IT_0073 = (-18)*conj(IT_0033);
    const ccomplex_t IT_0074 = m_t*s_13;
    const ccomplex_t IT_0075 = (-18)*IT_0065;
    const ccomplex_t IT_0076 = 18*IT_0068;
    const ccomplex_t IT_0077 = IT_0075 + IT_0076;
    const ccomplex_t IT_0078 = (-18)*conj(IT_0065);
    const ccomplex_t IT_0079 = 18*conj(IT_0068);
    const ccomplex_t IT_0080 = IT_0078 + IT_0079;
    const ccomplex_t IT_0081 = m_d*m_t;
    const ccomplex_t IT_0082 = (-18)*IT_0059;
    const ccomplex_t IT_0083 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0084 = m_d*s_14;
    const ccomplex_t IT_0085 = (IT_0033*conj(IT_0033) + IT_0048*conj(IT_0048))
      *IT_0053 + (IT_0059*conj(IT_0059) + IT_0065*conj(IT_0065) + IT_0067*conj
      (IT_0067) + IT_0068*conj(IT_0068))*IT_0069 + (conj(IT_0059)*IT_0067 +
       IT_0059*conj(IT_0067) + conj(IT_0065)*IT_0068 + IT_0065*conj(IT_0068))
      *IT_0070 + IT_0071*(conj(IT_0048)*IT_0072 + IT_0048*IT_0073) + IT_0074*(18
      *conj(IT_0033)*IT_0059 + 18*IT_0033*conj(IT_0059) + conj(IT_0067)*IT_0072 
      + IT_0067*IT_0073 + conj(IT_0048)*IT_0077 + IT_0048*IT_0080) + IT_0081*(18
      *conj(IT_0059)*IT_0065 + 18*IT_0059*conj(IT_0065) + conj(IT_0067)*IT_0077 
      + IT_0067*IT_0080 + conj(IT_0068)*IT_0082 + IT_0068*IT_0083) + (18*conj
      (IT_0033)*IT_0065 + 18*IT_0033*conj(IT_0065) + conj(IT_0068)*IT_0072 +
       IT_0068*IT_0073 + conj(IT_0048)*(18*IT_0067 + IT_0082) + IT_0048*(18*conj
      (IT_0067) + IT_0083))*IT_0084;
    return create_ccomplex_return(IT_0085);
}

