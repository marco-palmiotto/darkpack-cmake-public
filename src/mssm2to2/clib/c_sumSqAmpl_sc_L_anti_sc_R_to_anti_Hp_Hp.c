#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_R_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_R_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0000*IT_0006*(mu_h*U_sb_01 + m_b*IT_0002*IT_0004*U_sb_11);
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = pow(m_c, 2);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0006*(m_W*IT_0001*IT_0003*conj(U_sb_01) + (-0.5)*IT_0000*(IT_0001
      *IT_0004*IT_0009*conj(U_sb_01) + m_b*(m_b*IT_0002*IT_0003*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = pow(m_Hp, 2);
    const ccomplex_t IT_0014 = pow(m_sc_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_sc_L, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0014 + -IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = sin(alpha);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0004*IT_0006*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = cos(theta_W);
    const ccomplex_t IT_0024 = cpow(IT_0023, -2);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = -IT_0006 + IT_0025;
    const ccomplex_t IT_0027 = cos(alpha);
    const ccomplex_t IT_0028 = cpow(IT_0003, 3);
    const ccomplex_t IT_0029 = 3*IT_0006 + IT_0025;
    const ccomplex_t IT_0030 = cpow(IT_0003, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_W*e_em*(IT_0026*IT_0027
      *IT_0028 + -IT_0001*(IT_0001*(IT_0001*IT_0020*IT_0026 + IT_0003*IT_0027
      *IT_0029) + -IT_0020*IT_0029*IT_0030));
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0022*IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + -IT_0014 + -IT_0018 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0004*IT_0006*IT_0027;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = -IT_0006;
    const ccomplex_t IT_0040 = IT_0025 + IT_0039;
    const ccomplex_t IT_0041 = -IT_0001;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_W*e_em*(IT_0020*IT_0028
      *IT_0040 + IT_0041*(IT_0027*IT_0029*IT_0030 + -(IT_0003*IT_0020*IT_0029 + 
      -IT_0001*IT_0027*IT_0040)*IT_0041));
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*mu_h*IT_0000*IT_0006;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_c*m_d*V_cd
      *e_em*IT_0000*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0006;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*V_cs
      *e_em*IT_0000*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = pow(m_s, 2);
    const ccomplex_t IT_0064 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cs)
      *e_em*IT_0006*(IT_0001*(m_W*IT_0003 + (-0.5)*IT_0000*IT_0004*IT_0009) + 
      -1./2*IT_0000*IT_0002*IT_0003*IT_0063);
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0000*IT_0006;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = pow(m_d, 2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cd)
      *e_em*IT_0006*(IT_0001*(m_W*IT_0003 + (-0.5)*IT_0000*IT_0004*IT_0009) + 
      -1./2*IT_0000*IT_0002*IT_0003*IT_0072);
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*e_em
      *mu_h*IT_0000*IT_0006;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0006*(m_W*IT_0001*IT_0003*conj(U_sb_00) + (-0.5)*IT_0000*(IT_0001
      *IT_0004*IT_0009*conj(U_sb_00) + m_b*(m_b*IT_0002*IT_0003*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0000*IT_0006*(mu_h*U_sb_00 + m_b*IT_0002*IT_0004*U_sb_10);
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0081*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0017 + -IT_0035 + -IT_0046 + IT_0054 +
       IT_0062 + IT_0071 + IT_0080 + IT_0088;
    return create_ccomplex_return(3*IT_0089*conj(IT_0089));
}

