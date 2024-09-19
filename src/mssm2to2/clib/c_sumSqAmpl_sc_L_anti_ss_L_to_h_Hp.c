#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_h_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_h_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -2);
    const ccomplex_t IT_0004 = pow(m_c, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = cpow(IT_0001, -2);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *IT_0000*(IT_0001*(IT_0002 + -IT_0002*IT_0003*IT_0004*IT_0006) + -IT_0005*
      (1 + -IT_0003*IT_0007*IT_0008)*IT_0009)*IT_0011;
    const ccomplex_t IT_0013 = (-0.25)*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0010, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 +
       IT_0005*IT_0009)*IT_0014;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = pow(m_sc_L, 2);
    const ccomplex_t IT_0021 = pow(m_ss_L, 2);
    const ccomplex_t IT_0022 = 0.5*IT_0020 + 0.5*IT_0021;
    const ccomplex_t IT_0023 = s_12 + IT_0022;
    const ccomplex_t IT_0024 = IT_0019*IT_0023;
    const ccomplex_t IT_0025 = -IT_0017;
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0020 + 0.5
      *IT_0021;
    const ccomplex_t IT_0028 = s_12 + IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = s_12 + IT_0021;
    const ccomplex_t IT_0031 = -IT_0015;
    const ccomplex_t IT_0032 = IT_0018*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = IT_0025*IT_0031;
    const ccomplex_t IT_0035 = -s_23;
    const ccomplex_t IT_0036 = IT_0030 + IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = -IT_0029 + (-0.5)*IT_0033 + 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0024 + IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, -1);
    const ccomplex_t IT_0044 = cpow(IT_0001, -1);
    const ccomplex_t IT_0045 = cpow(IT_0005, -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0014*(IT_0005*IT_0007*IT_0043*IT_0044 + (-2)*IT_0001*(m_W*IT_0005
       + (-0.5)*IT_0004*IT_0043*IT_0045));
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = cos(theta_W);
    const ccomplex_t IT_0049 = cpow(IT_0048, -2);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0004
      *IT_0014*IT_0043*IT_0045 + (-0.5)*m_W*(IT_0002*IT_0005 + IT_0001*IT_0009)*
      (IT_0014 + -1./3*IT_0010*IT_0049));
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0047*IT_0051;
    const ccomplex_t IT_0053 = pow(m_h, 2);
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0053 + (0 + _Complex_I*1)
      *m_sc_L*Gamma_cl + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em*(IT_0007*IT_0009
      *IT_0014*IT_0043*IT_0044 + (-0.5)*m_W*(IT_0002*IT_0005 + IT_0001*IT_0009)*
      (IT_0014 + 1./3*IT_0010*IT_0049));
    const ccomplex_t IT_0058 = IT_0047*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0053 + (0 + _Complex_I*1)
      *m_ss_L*Gamma_sl + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0014*(m_W*IT_0001*IT_0005*conj(U_sb_00) + (-0.5)*IT_0043*(IT_0001
      *IT_0004*IT_0045*conj(U_sb_00) + m_b*(m_b*IT_0005*IT_0044*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_10;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0021 + IT_0053 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0010*IT_0049;
    const ccomplex_t IT_0071 = -IT_0014 + IT_0070;
    const ccomplex_t IT_0072 = cpow(IT_0005, 3);
    const ccomplex_t IT_0073 = 3*IT_0014 + IT_0070;
    const ccomplex_t IT_0074 = cpow(IT_0005, 2);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_W*e_em*(IT_0002*IT_0071
      *IT_0072 + -IT_0001*(IT_0001*(IT_0001*IT_0009*IT_0071 + IT_0002*IT_0005
      *IT_0073) + -IT_0009*IT_0073*IT_0074));
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0047*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0009
      *IT_0014*IT_0043*IT_0045;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0014*IT_0043;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0020 + IT_0053 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*mu_h*IT_0014*IT_0043;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0009;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0021 + IT_0053 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0014*IT_0043;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0009 + m_s*e_em*mu_h*IT_0002*IT_0014*IT_0043*IT_0044);
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0021 + IT_0053 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0014*(m_W*IT_0001*IT_0005*conj(U_sb_01) + (-0.5)*IT_0043*(IT_0001
      *IT_0004*IT_0045*conj(U_sb_01) + m_b*(m_b*IT_0005*IT_0044*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_11;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0106*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0021 + IT_0053 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = cpow(IT_0001, 2);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_W*e_em*(IT_0009*IT_0014
      *IT_0072 + (-2)*IT_0001*(IT_0005*(IT_0002*IT_0005 + IT_0001*IT_0009)*(0.5
      *IT_0014 + IT_0010*IT_0049) + -1./2*IT_0002*IT_0014*IT_0113));
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0014*((IT_0004 + -IT_0007)*IT_0043 + -m_W*(IT_0074 + -IT_0113));
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = pow(m_W, 2);
    const ccomplex_t IT_0120 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + IT_0119 
      + -reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0013 + (-2)*IT_0042 + IT_0056 + IT_0061 +
       IT_0069 + -IT_0080 + IT_0088 + IT_0096 + IT_0104 + IT_0112 + -IT_0122;
    return create_ccomplex_return(3*IT_0123*conj(IT_0123));
}

