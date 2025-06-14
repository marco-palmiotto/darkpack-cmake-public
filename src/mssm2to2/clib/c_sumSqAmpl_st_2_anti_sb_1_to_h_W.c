#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_1_to_h_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_1_to_h_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0003 = (-0.5)*IT_0002;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*(IT_0004
      *IT_0005 + -IT_0006*IT_0007);
    const ccomplex_t IT_0009 = pow(m_sb_1, 2);
    const ccomplex_t IT_0010 = pow(m_st_2, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -2);
    const ccomplex_t IT_0017 = pow(s_14, 2);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = -IT_0002;
    const ccomplex_t IT_0020 = cpow(V_tb, 2);
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_b*mu_h*IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0006, -1);
    const ccomplex_t IT_0024 = IT_0006*IT_0007;
    const ccomplex_t IT_0025 = IT_0004*IT_0005;
    const ccomplex_t IT_0026 = IT_0024 + IT_0025;
    const ccomplex_t IT_0027 = cos(theta_W);
    const ccomplex_t IT_0028 = cpow(IT_0027, -2);
    const ccomplex_t IT_0029 = pow(m_b, 2);
    const ccomplex_t IT_0030 = 2*U_sb_00*(((-0.5)*m_W*IT_0026*(IT_0001 + 1./3
      *IT_0000*IT_0028) + IT_0001*IT_0007*IT_0021*IT_0023*IT_0029)*conj(U_sb_00)
       + 1./2*IT_0001*IT_0004*IT_0022*IT_0023*conj(U_sb_10)) + U_sb_10*(IT_0001
      *IT_0004*IT_0022*IT_0023*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0000
      *IT_0026*IT_0028 + (-3)*IT_0001*IT_0007*IT_0021*IT_0023*IT_0029)*conj
      (U_sb_10));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_b*IT_0007*IT_0020*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = pow(m_h, 2);
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0033 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0019*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0001*U_st_01;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0007*conj(U_sb_10);
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0009 + IT_0033 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0005 + -IT_0006*IT_0007);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0005, -1);
    const ccomplex_t IT_0049 = pow(m_t, 2);
    const ccomplex_t IT_0050 = cpow(IT_0006, 2);
    const ccomplex_t IT_0051 = cpow(IT_0005, 2);
    const ccomplex_t IT_0052 = m_b*mu_h*IT_0005*IT_0021*IT_0023*conj(U_sb_10)
      *U_st_01 + -conj(U_sb_00)*((IT_0021*(IT_0029 + -IT_0049) + -m_W*(IT_0050 +
       -IT_0051))*U_st_01 + m_t*mu_h*IT_0006*IT_0021*IT_0048*U_st_11);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0001*IT_0052 + A_b*IT_0006*conj(U_sb_10)*U_st_01 + -A_t*IT_0005
      *conj(U_sb_00)*U_st_11);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = pow(m_W, 2);
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + IT_0055 
      + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0006 + IT_0005*IT_0007);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = m_t*conj(U_sb_00)*(m_t*IT_0006*IT_0048*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0005*IT_0023*conj(U_sb_00)*U_st_01 + conj
      (U_sb_10)*(mu_h*U_st_01 + m_t*IT_0023*IT_0048*U_st_11));
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_01 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0021
      *IT_0062 + m_W*IT_0005*IT_0006*conj(U_sb_00)*U_st_01) + A_t*IT_0006*conj
      (U_sb_00)*U_st_11);
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = IT_0061*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = 2*conj(U_sb_00)*(((-0.5)*m_W*IT_0026*(IT_0001 +
       1./3*IT_0000*IT_0028) + IT_0001*IT_0007*IT_0021*IT_0023*IT_0029)*U_sb_01 
      + 1./2*IT_0001*IT_0004*IT_0022*IT_0023*U_sb_11) + conj(U_sb_10)*(IT_0001
      *IT_0004*IT_0022*IT_0023*U_sb_01 + (-0.666666666666667)*(m_W*IT_0000
      *IT_0026*IT_0028 + (-3)*IT_0001*IT_0007*IT_0021*IT_0023*IT_0029)*U_sb_11);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(e_em*IT_0070 +
       1.4142135623731*A_b*IT_0007*IT_0020*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0009 + IT_0033 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0069*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0001*U_st_01;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0007*conj(U_sb_10);
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0009 + IT_0033 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0078*IT_0083;
    const ccomplex_t IT_0085 = -IT_0037 + -IT_0045 + IT_0059 + IT_0067 + 
      -IT_0076 + -IT_0084;
    const ccomplex_t IT_0086 = s_14*s_34*IT_0016;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = s_13 + IT_0087;
    const ccomplex_t IT_0089 = -IT_0019;
    const ccomplex_t IT_0090 = IT_0003 + IT_0089;
    const ccomplex_t IT_0091 = IT_0013*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = m_t*mu_h*IT_0021;
    const ccomplex_t IT_0095 = conj(U_st_10)*(IT_0001*IT_0007*IT_0048*IT_0094
      *U_st_01 + (-1.33333333333333)*(m_W*IT_0000*IT_0026*IT_0028 + (-1.5)
      *IT_0001*IT_0004*IT_0021*IT_0048*IT_0049)*U_st_11) + 2*conj(U_st_00)*(((
      -0.5)*m_W*IT_0026*(IT_0001 + -1./3*IT_0000*IT_0028) + IT_0001*IT_0004
      *IT_0021*IT_0048*IT_0049)*U_st_01 + 1./2*IT_0001*IT_0007*IT_0048*IT_0094
      *U_st_11);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(e_em*IT_0095 +
       1.4142135623731*A_t*IT_0004*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_13 + IT_0010 + IT_0033 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0093*IT_0100;
    const ccomplex_t IT_0102 = U_st_11*(IT_0001*IT_0007*IT_0048*IT_0094*conj
      (U_st_01) + (-1.33333333333333)*(m_W*IT_0000*IT_0026*IT_0028 + (-1.5)
      *IT_0001*IT_0004*IT_0021*IT_0048*IT_0049)*conj(U_st_11)) + 2*U_st_01*(((
      -0.5)*m_W*IT_0026*(IT_0001 + -1./3*IT_0000*IT_0028) + IT_0001*IT_0004
      *IT_0021*IT_0048*IT_0049)*conj(U_st_01) + 1./2*IT_0001*IT_0007*IT_0048
      *IT_0094*conj(U_st_11));
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(e_em*IT_0102 +
       1.4142135623731*A_t*IT_0004*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_13 + IT_0033 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0019*IT_0107;
    const ccomplex_t IT_0109 = IT_0037 + IT_0045 + IT_0076 + IT_0084 + 
      -IT_0091 + IT_0101 + IT_0108;
    const ccomplex_t IT_0110 = s_14*s_24*IT_0016;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = s_12 + IT_0111;
    const ccomplex_t IT_0113 = (-3)*conj(IT_0015);
    const ccomplex_t IT_0114 = pow(s_34, 2);
    const ccomplex_t IT_0115 = IT_0016*IT_0114;
    const ccomplex_t IT_0116 = s_24*s_34*IT_0016;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = s_23 + IT_0117;
    const ccomplex_t IT_0119 = pow(s_24, 2);
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = (-3)*IT_0109;
    const ccomplex_t IT_0122 = (-3)*IT_0015*(conj(IT_0015)*(IT_0010 + -IT_0018
      ) + conj(IT_0085)*IT_0088 + conj(IT_0109)*IT_0112) + IT_0109*IT_0112
      *IT_0113 + (-3)*IT_0085*((-0.333333333333333)*IT_0088*IT_0113 + conj
      (IT_0085)*(IT_0033 + -IT_0115) + conj(IT_0109)*IT_0118) + (conj(IT_0085)
      *IT_0118 + conj(IT_0109)*(IT_0009 + -IT_0120))*IT_0121;
    return create_ccomplex_return(IT_0122);
}

