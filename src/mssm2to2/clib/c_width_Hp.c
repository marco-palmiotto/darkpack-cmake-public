#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Hp, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = pow(m_e, 2);
    const ccomplex_t IT_0012 = pow(m_Hp, 2);
    const ccomplex_t IT_0013 = cpow(IT_0011 + -IT_0012, 2);
    const ccomplex_t IT_0014 = (-0.25)*IT_0001*(IT_0011 + IT_0012)*(IT_0011 + 
      -IT_0012) + 0.25*IT_0001*IT_0013;
    const ccomplex_t IT_0015 = cabs(0.25*IT_0001*IT_0013);
    const ccomplex_t IT_0016 = cpow(IT_0015, 0.5);
    const ccomplex_t IT_0017 = ( ((m_Hp
) >= (m_e
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0010*IT_0014*IT_0016
) : (0
));
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*-1)*e_em*(IT_0007*IT_0018 + 
      -IT_0006*IT_0019);
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0021, 2);
    const ccomplex_t IT_0023 = pow(m_Z, -2);
    const ccomplex_t IT_0024 = pow(m_Z, 2);
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 4);
    const ccomplex_t IT_0029 = cabs(-IT_0024 + 0.25*IT_0001*IT_0028);
    const ccomplex_t IT_0030 = cpow(IT_0029, 0.5);
    const ccomplex_t IT_0031 = ( ((m_Hp
) >= (m_Z + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0022*(IT_0012 + -IT_0027)*IT_0030
) : (0
));
    const ccomplex_t IT_0032 = cos(alpha);
    const ccomplex_t IT_0033 = sin(alpha);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*IT_0007*(IT_0003
      *IT_0032 + IT_0005*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0036, 2);
    const ccomplex_t IT_0038 = pow(m_W, -2);
    const ccomplex_t IT_0039 = pow(m_W, 2);
    const ccomplex_t IT_0040 = pow(m_h, 2);
    const ccomplex_t IT_0041 = 0.5*IT_0012 + 0.5*IT_0039 + (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0041, 2);
    const ccomplex_t IT_0043 = IT_0038*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0012 + IT_0039 + -IT_0040, 2);
    const ccomplex_t IT_0045 = cabs(-IT_0039 + 0.25*IT_0001*IT_0044);
    const ccomplex_t IT_0046 = cpow(IT_0045, 0.5);
    const ccomplex_t IT_0047 = ( ((m_Hp
) >= (m_W + m_h
)) ? (0.125*IT_0000
      *IT_0001*IT_0037*(IT_0012 + -IT_0043)*IT_0046
) : (0
));
    const ccomplex_t IT_0048 = cpow(IT_0005, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_d*V_td*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = m_d*m_t;
    const ccomplex_t IT_0054 = pow(m_d, 2);
    const ccomplex_t IT_0055 = pow(m_t, 2);
    const ccomplex_t IT_0056 = cpow(IT_0012 + IT_0054 + -IT_0055, 2);
    const ccomplex_t IT_0057 = -IT_0054 + 0.25*IT_0001*(IT_0012 + -IT_0054 +
       IT_0055)*(IT_0012 + IT_0054 + -IT_0055) + 0.25*IT_0001*IT_0056;
    const ccomplex_t IT_0058 = 6*IT_0057;
    const ccomplex_t IT_0059 = cabs(-IT_0054 + 0.25*IT_0001*IT_0056);
    const ccomplex_t IT_0060 = cpow(IT_0059, 0.5);
    const ccomplex_t IT_0061 = ( ((m_Hp
) >= (m_d + m_t
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0050)*(IT_0052*IT_0053 + (-0.166666666666667)*IT_0050
      *IT_0058) + conj(IT_0052)*(IT_0050*IT_0053 + (-0.166666666666667)*IT_0052
      *IT_0058))*IT_0060
) : (0
));
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_s*V_ts*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = m_s*m_t;
    const ccomplex_t IT_0067 = pow(m_s, 2);
    const ccomplex_t IT_0068 = cpow(IT_0012 + -IT_0055 + IT_0067, 2);
    const ccomplex_t IT_0069 = -IT_0067 + 0.25*IT_0001*(IT_0012 + -IT_0055 +
       IT_0067)*(IT_0012 + IT_0055 + -IT_0067) + 0.25*IT_0001*IT_0068;
    const ccomplex_t IT_0070 = 6*IT_0069;
    const ccomplex_t IT_0071 = cabs(-IT_0067 + 0.25*IT_0001*IT_0068);
    const ccomplex_t IT_0072 = cpow(IT_0071, 0.5);
    const ccomplex_t IT_0073 = ( ((m_Hp
) >= (m_s + m_t
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0063)*(IT_0065*IT_0066 + (-0.166666666666667)*IT_0063
      *IT_0070) + conj(IT_0065)*(IT_0063*IT_0066 + (-0.166666666666667)*IT_0065
      *IT_0070))*IT_0072
) : (0
));
    const ccomplex_t IT_0074 = e_em*IT_0007;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = cpow(IT_0076, 2);
    const ccomplex_t IT_0078 = pow(m_A0, 2);
    const ccomplex_t IT_0079 = 0.5*IT_0012 + 0.5*IT_0039 + (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0079, 2);
    const ccomplex_t IT_0081 = IT_0038*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0012 + IT_0039 + -IT_0078, 2);
    const ccomplex_t IT_0083 = cabs(-IT_0039 + 0.25*IT_0001*IT_0082);
    const ccomplex_t IT_0084 = cpow(IT_0083, 0.5);
    const ccomplex_t IT_0085 = ( ((m_Hp
) >= (m_W + m_A0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0077*(IT_0012 + -IT_0081)*IT_0084
) : (0
));
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*e_em*IT_0007*(IT_0005
      *IT_0032 + -IT_0003*IT_0033);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = cpow(IT_0087, 2);
    const ccomplex_t IT_0089 = pow(m_H0, 2);
    const ccomplex_t IT_0090 = 0.5*IT_0012 + 0.5*IT_0039 + (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = cpow(IT_0090, 2);
    const ccomplex_t IT_0092 = IT_0038*IT_0091;
    const ccomplex_t IT_0093 = cpow(IT_0012 + IT_0039 + -IT_0089, 2);
    const ccomplex_t IT_0094 = cabs(-IT_0039 + 0.25*IT_0001*IT_0093);
    const ccomplex_t IT_0095 = cpow(IT_0094, 0.5);
    const ccomplex_t IT_0096 = ( ((m_Hp
) >= (m_W + m_H0
)) ? (0.125*IT_0000
      *IT_0001*IT_0088*(IT_0012 + -IT_0092)*IT_0095
) : (0
));
    const ccomplex_t IT_0097 = cpow(IT_0018, -2);
    const ccomplex_t IT_0098 = IT_0006*IT_0097;
    const ccomplex_t IT_0099 = -IT_0007 + IT_0098;
    const ccomplex_t IT_0100 = cpow(IT_0005, 3);
    const ccomplex_t IT_0101 = 3*IT_0007 + IT_0098;
    const ccomplex_t IT_0102 = cpow(IT_0005, 2);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_W*e_em*(IT_0032*IT_0099
      *IT_0100 + -IT_0003*(IT_0003*(IT_0003*IT_0033*IT_0099 + IT_0005*IT_0032
      *IT_0101) + -IT_0033*IT_0101*IT_0102));
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0104, 2);
    const ccomplex_t IT_0106 = cpow((-2)*IT_0012 + IT_0040, 2);
    const ccomplex_t IT_0107 = cabs(-IT_0012 + 0.25*IT_0001*IT_0106);
    const ccomplex_t IT_0108 = cpow(IT_0107, 0.5);
    const ccomplex_t IT_0109 = ( ((m_Hp
) >= (m_h + m_Hp
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0105*IT_0108
) : (0
));
    const ccomplex_t IT_0110 = -IT_0007;
    const ccomplex_t IT_0111 = IT_0098 + IT_0110;
    const ccomplex_t IT_0112 = -IT_0003;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*m_W*e_em*(IT_0033*IT_0100
      *IT_0111 + IT_0112*(IT_0032*IT_0101*IT_0102 + -(IT_0005*IT_0033*IT_0101 + 
      -IT_0003*IT_0032*IT_0111)*IT_0112));
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = cpow(IT_0114, 2);
    const ccomplex_t IT_0116 = cpow((-2)*IT_0012 + IT_0089, 2);
    const ccomplex_t IT_0117 = cabs(-IT_0012 + 0.25*IT_0001*IT_0116);
    const ccomplex_t IT_0118 = cpow(IT_0117, 0.5);
    const ccomplex_t IT_0119 = ( ((m_Hp
) >= (m_H0 + m_Hp
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0115*IT_0118
) : (0
));
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = cpow(IT_0121, 2);
    const ccomplex_t IT_0123 = pow(m_c, 2);
    const ccomplex_t IT_0124 = pow(m_b, 2);
    const ccomplex_t IT_0125 = cpow(-IT_0012 + -IT_0123 + IT_0124, 2);
    const ccomplex_t IT_0126 = -IT_0123 + 0.25*IT_0001*(IT_0012 + -IT_0123 +
       IT_0124)*(IT_0012 + IT_0123 + -IT_0124) + 0.25*IT_0001*IT_0125;
    const ccomplex_t IT_0127 = (-6)*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = m_b*m_c;
    const ccomplex_t IT_0131 = cabs(-IT_0123 + 0.25*IT_0001*IT_0125);
    const ccomplex_t IT_0132 = cpow(IT_0131, 0.5);
    const ccomplex_t IT_0133 = ( ((m_Hp
) >= (m_b + m_c
)) ? (0.125*IT_0000
      *IT_0001*(IT_0122*IT_0127 + IT_0129*(IT_0127*IT_0129 + 12*IT_0121*IT_0130)
      )*IT_0132
) : (0
));
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*m_d*V_cd*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = m_c*m_d;
    const ccomplex_t IT_0139 = cpow(-IT_0012 + -IT_0054 + IT_0123, 2);
    const ccomplex_t IT_0140 = -IT_0054 + 0.25*IT_0001*(IT_0012 + -IT_0054 +
       IT_0123)*(IT_0012 + IT_0054 + -IT_0123) + 0.25*IT_0001*IT_0139;
    const ccomplex_t IT_0141 = 6*IT_0140;
    const ccomplex_t IT_0142 = cabs(-IT_0054 + 0.25*IT_0001*IT_0139);
    const ccomplex_t IT_0143 = cpow(IT_0142, 0.5);
    const ccomplex_t IT_0144 = ( ((m_Hp
) >= (m_c + m_d
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0135)*(IT_0137*IT_0138 + (-0.166666666666667)*IT_0135
      *IT_0141) + conj(IT_0137)*(IT_0135*IT_0138 + (-0.166666666666667)*IT_0137
      *IT_0141))*IT_0143
) : (0
));
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = m_c*m_s;
    const ccomplex_t IT_0150 = cpow(-IT_0012 + -IT_0067 + IT_0123, 2);
    const ccomplex_t IT_0151 = -IT_0067 + 0.25*IT_0001*(IT_0012 + -IT_0067 +
       IT_0123)*(IT_0012 + IT_0067 + -IT_0123) + 0.25*IT_0001*IT_0150;
    const ccomplex_t IT_0152 = 6*IT_0151;
    const ccomplex_t IT_0153 = cabs(-IT_0067 + 0.25*IT_0001*IT_0150);
    const ccomplex_t IT_0154 = cpow(IT_0153, 0.5);
    const ccomplex_t IT_0155 = ( ((m_Hp
) >= (m_c + m_s
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0146)*(IT_0148*IT_0149 + (-0.166666666666667)*IT_0146
      *IT_0152) + conj(IT_0148)*(IT_0146*IT_0149 + (-0.166666666666667)*IT_0148
      *IT_0152))*IT_0154
) : (0
));
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = m_b*m_t;
    const ccomplex_t IT_0161 = cpow(-IT_0012 + -IT_0055 + IT_0124, 2);
    const ccomplex_t IT_0162 = -IT_0055 + 0.25*IT_0001*(IT_0012 + -IT_0055 +
       IT_0124)*(IT_0012 + IT_0055 + -IT_0124) + 0.25*IT_0001*IT_0161;
    const ccomplex_t IT_0163 = 6*IT_0162;
    const ccomplex_t IT_0164 = cabs(-IT_0055 + 0.25*IT_0001*IT_0161);
    const ccomplex_t IT_0165 = cpow(IT_0164, 0.5);
    const ccomplex_t IT_0166 = ( ((m_Hp
) >= (m_b + m_t
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0157)*(IT_0159*IT_0160 + (-0.166666666666667)*IT_0157
      *IT_0163) + conj(IT_0159)*(IT_0157*IT_0160 + (-0.166666666666667)*IT_0159
      *IT_0163))*IT_0165
) : (0
));
    const ccomplex_t IT_0167 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007*IT_0167*V_ub_mod;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = cpow(IT_0169, 2);
    const ccomplex_t IT_0171 = pow(m_u, 2);
    const ccomplex_t IT_0172 = cpow(-IT_0012 + IT_0124 + -IT_0171, 2);
    const ccomplex_t IT_0173 = -IT_0171 + 0.25*IT_0001*(IT_0012 + -IT_0124 +
       IT_0171)*(IT_0012 + IT_0124 + -IT_0171) + 0.25*IT_0001*IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048*IT_0167*V_ub_mod;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = m_b*m_u;
    const ccomplex_t IT_0178 = cabs(-IT_0171 + 0.25*IT_0001*IT_0172);
    const ccomplex_t IT_0179 = cpow(IT_0178, 0.5);
    const ccomplex_t IT_0180 = ( ((m_Hp
) >= (m_b + m_u
)) ? (0.125*IT_0000
      *IT_0001*(IT_0170*IT_0174 + IT_0176*(IT_0174*IT_0176 + (-12)*IT_0169
      *IT_0177))*IT_0179
) : (0
));
    const ccomplex_t IT_0181 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = cpow(IT_0182, 2);
    const ccomplex_t IT_0184 = cpow(-IT_0012 + IT_0054 + -IT_0171, 2);
    const ccomplex_t IT_0185 = -IT_0171 + 0.25*IT_0001*(IT_0012 + -IT_0054 +
       IT_0171)*(IT_0012 + IT_0054 + -IT_0171) + 0.25*IT_0001*IT_0184;
    const ccomplex_t IT_0186 = (-6)*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = m_d*m_u;
    const ccomplex_t IT_0190 = cabs(-IT_0171 + 0.25*IT_0001*IT_0184);
    const ccomplex_t IT_0191 = cpow(IT_0190, 0.5);
    const ccomplex_t IT_0192 = ( ((m_Hp
) >= (m_d + m_u
)) ? (0.125*IT_0000
      *IT_0001*(IT_0183*IT_0186 + IT_0188*(IT_0186*IT_0188 + 12*IT_0182*IT_0189)
      )*IT_0191
) : (0
));
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0194 = 0.5*IT_0193;
    const ccomplex_t IT_0195 = cpow(IT_0194, 2);
    const ccomplex_t IT_0196 = cpow(-IT_0012 + IT_0067 + -IT_0171, 2);
    const ccomplex_t IT_0197 = -IT_0171 + 0.25*IT_0001*(IT_0012 + -IT_0067 +
       IT_0171)*(IT_0012 + IT_0067 + -IT_0171) + 0.25*IT_0001*IT_0196;
    const ccomplex_t IT_0198 = (-6)*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0002*IT_0003*IT_0007*IT_0048;
    const ccomplex_t IT_0200 = 0.5*IT_0199;
    const ccomplex_t IT_0201 = m_s*m_u;
    const ccomplex_t IT_0202 = cabs(-IT_0171 + 0.25*IT_0001*IT_0196);
    const ccomplex_t IT_0203 = cpow(IT_0202, 0.5);
    const ccomplex_t IT_0204 = ( ((m_Hp
) >= (m_s + m_u
)) ? (0.125*IT_0000
      *IT_0001*(IT_0195*IT_0198 + IT_0200*(IT_0198*IT_0200 + 12*IT_0194*IT_0201)
      )*IT_0203
) : (0
));
    const ccomplex_t IT_0205 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = cpow(IT_0206, 2);
    const ccomplex_t IT_0208 = pow(m_mu, 2);
    const ccomplex_t IT_0209 = cpow(IT_0012 + -IT_0208, 2);
    const ccomplex_t IT_0210 = 0.25*IT_0001*(IT_0012 + IT_0208)*(IT_0012 + 
      -IT_0208) + 0.25*IT_0001*IT_0209;
    const ccomplex_t IT_0211 = cabs(0.25*IT_0001*IT_0209);
    const ccomplex_t IT_0212 = cpow(IT_0211, 0.5);
    const ccomplex_t IT_0213 = ( ((m_Hp
) >= (m_mu
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0207*IT_0210*IT_0212
) : (0
));
    const ccomplex_t IT_0214 = IT_0005*IT_0007;
    const ccomplex_t IT_0215 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0216;
    const ccomplex_t IT_0218 = IT_0005*IT_0019;
    const ccomplex_t IT_0219 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0220 = IT_0218*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0223 = IT_0214*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = IT_0221 + IT_0224;
    const ccomplex_t IT_0226 = IT_0217 + IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = 0.5*IT_0227;
    const ccomplex_t IT_0229 = IT_0003*IT_0007;
    const ccomplex_t IT_0230 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0231 = IT_0229*IT_0230;
    const ccomplex_t IT_0232 = 2*IT_0231;
    const ccomplex_t IT_0233 = IT_0003*IT_0019;
    const ccomplex_t IT_0234 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0235 = IT_0233*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0238 = IT_0229*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = IT_0236 + IT_0239;
    const ccomplex_t IT_0241 = IT_0232 + IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = 0.5*IT_0242;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = m_C_1*m_N_1;
    const ccomplex_t IT_0246 = pow(m_N_1, 2);
    const ccomplex_t IT_0247 = pow(m_C_1, 2);
    const ccomplex_t IT_0248 = cpow(IT_0012 + IT_0246 + -IT_0247, 2);
    const ccomplex_t IT_0249 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0247)*(IT_0012 + IT_0246 + -IT_0247) + 0.25*IT_0001*IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0249;
    const ccomplex_t IT_0251 = cabs(-IT_0246 + 0.25*IT_0001*IT_0248);
    const ccomplex_t IT_0252 = cpow(IT_0251, 0.5);
    const ccomplex_t IT_0253 = ( ((m_Hp
) >= (m_C_1 + m_N_1
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0228)*(IT_0244*IT_0245 + (-0.5)*IT_0228*IT_0250)
       + conj(IT_0244)*(IT_0228*IT_0245 + (-0.5)*IT_0244*IT_0250))*IT_0252
) : 
      (0
));
    const ccomplex_t IT_0254 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0255 = IT_0214*IT_0254;
    const ccomplex_t IT_0256 = (-2)*IT_0255;
    const ccomplex_t IT_0257 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0258 = IT_0218*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0261 = IT_0214*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = IT_0259 + IT_0262;
    const ccomplex_t IT_0264 = IT_0256 + IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = 0.5*IT_0265;
    const ccomplex_t IT_0267 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0268 = IT_0229*IT_0267;
    const ccomplex_t IT_0269 = 2*IT_0268;
    const ccomplex_t IT_0270 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0271 = IT_0233*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0274 = IT_0229*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = IT_0272 + IT_0275;
    const ccomplex_t IT_0277 = IT_0269 + IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = 0.5*IT_0278;
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = m_C_2*m_N_1;
    const ccomplex_t IT_0282 = pow(m_C_2, 2);
    const ccomplex_t IT_0283 = cpow(IT_0012 + IT_0246 + -IT_0282, 2);
    const ccomplex_t IT_0284 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0282)*(IT_0012 + IT_0246 + -IT_0282) + 0.25*IT_0001*IT_0283;
    const ccomplex_t IT_0285 = 2*IT_0284;
    const ccomplex_t IT_0286 = cabs(-IT_0246 + 0.25*IT_0001*IT_0283);
    const ccomplex_t IT_0287 = cpow(IT_0286, 0.5);
    const ccomplex_t IT_0288 = ( ((m_Hp
) >= (m_C_2 + m_N_1
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0266)*(IT_0280*IT_0281 + (-0.5)*IT_0266*IT_0285)
       + conj(IT_0280)*(IT_0266*IT_0281 + (-0.5)*IT_0280*IT_0285))*IT_0287
) : 
      (0
));
    const ccomplex_t IT_0289 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0290 = IT_0214*IT_0289;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0293 = IT_0218*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0296 = IT_0214*IT_0295;
    const ccomplex_t IT_0297 = 1.4142135623731*IT_0296;
    const ccomplex_t IT_0298 = IT_0294 + IT_0297;
    const ccomplex_t IT_0299 = IT_0291 + IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = 0.5*IT_0300;
    const ccomplex_t IT_0302 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0303 = IT_0229*IT_0302;
    const ccomplex_t IT_0304 = 2*IT_0303;
    const ccomplex_t IT_0305 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0306 = IT_0233*IT_0305;
    const ccomplex_t IT_0307 = 1.4142135623731*IT_0306;
    const ccomplex_t IT_0308 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0309 = IT_0229*IT_0308;
    const ccomplex_t IT_0310 = 1.4142135623731*IT_0309;
    const ccomplex_t IT_0311 = IT_0307 + IT_0310;
    const ccomplex_t IT_0312 = IT_0304 + IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*IT_0312;
    const ccomplex_t IT_0314 = 0.5*IT_0313;
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = m_C_1*m_N_2;
    const ccomplex_t IT_0317 = pow(m_N_2, 2);
    const ccomplex_t IT_0318 = cpow(IT_0012 + -IT_0247 + IT_0317, 2);
    const ccomplex_t IT_0319 = -IT_0317 + 0.25*IT_0001*(IT_0012 + -IT_0247 +
       IT_0317)*(IT_0012 + IT_0247 + -IT_0317) + 0.25*IT_0001*IT_0318;
    const ccomplex_t IT_0320 = 2*IT_0319;
    const ccomplex_t IT_0321 = cabs(-IT_0317 + 0.25*IT_0001*IT_0318);
    const ccomplex_t IT_0322 = cpow(IT_0321, 0.5);
    const ccomplex_t IT_0323 = ( ((m_Hp
) >= (m_C_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0301)*(IT_0315*IT_0316 + (-0.5)*IT_0301*IT_0320)
       + conj(IT_0315)*(IT_0301*IT_0316 + (-0.5)*IT_0315*IT_0320))*IT_0322
) : 
      (0
));
    const ccomplex_t IT_0324 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0325 = IT_0214*IT_0324;
    const ccomplex_t IT_0326 = (-2)*IT_0325;
    const ccomplex_t IT_0327 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0328 = IT_0218*IT_0327;
    const ccomplex_t IT_0329 = 1.4142135623731*IT_0328;
    const ccomplex_t IT_0330 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0331 = IT_0214*IT_0330;
    const ccomplex_t IT_0332 = 1.4142135623731*IT_0331;
    const ccomplex_t IT_0333 = IT_0329 + IT_0332;
    const ccomplex_t IT_0334 = IT_0326 + IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = 0.5*IT_0335;
    const ccomplex_t IT_0337 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0338 = IT_0229*IT_0337;
    const ccomplex_t IT_0339 = 2*IT_0338;
    const ccomplex_t IT_0340 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0341 = IT_0233*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0344 = IT_0229*IT_0343;
    const ccomplex_t IT_0345 = 1.4142135623731*IT_0344;
    const ccomplex_t IT_0346 = IT_0342 + IT_0345;
    const ccomplex_t IT_0347 = IT_0339 + IT_0346;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0348;
    const ccomplex_t IT_0350 = -IT_0349;
    const ccomplex_t IT_0351 = m_C_2*m_N_2;
    const ccomplex_t IT_0352 = cpow(IT_0012 + -IT_0282 + IT_0317, 2);
    const ccomplex_t IT_0353 = -IT_0317 + 0.25*IT_0001*(IT_0012 + -IT_0282 +
       IT_0317)*(IT_0012 + IT_0282 + -IT_0317) + 0.25*IT_0001*IT_0352;
    const ccomplex_t IT_0354 = 2*IT_0353;
    const ccomplex_t IT_0355 = cabs(-IT_0317 + 0.25*IT_0001*IT_0352);
    const ccomplex_t IT_0356 = cpow(IT_0355, 0.5);
    const ccomplex_t IT_0357 = ( ((m_Hp
) >= (m_C_2 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0336)*(IT_0350*IT_0351 + (-0.5)*IT_0336*IT_0354)
       + conj(IT_0350)*(IT_0336*IT_0351 + (-0.5)*IT_0350*IT_0354))*IT_0356
) : 
      (0
));
    const ccomplex_t IT_0358 = N_d3*e_em*U_Wm1;
    const ccomplex_t IT_0359 = IT_0214*IT_0358;
    const ccomplex_t IT_0360 = (-2)*IT_0359;
    const ccomplex_t IT_0361 = N_B3*U_d1*e_em;
    const ccomplex_t IT_0362 = IT_0218*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = N_W3*U_d1*e_em;
    const ccomplex_t IT_0365 = IT_0214*IT_0364;
    const ccomplex_t IT_0366 = 1.4142135623731*IT_0365;
    const ccomplex_t IT_0367 = IT_0363 + IT_0366;
    const ccomplex_t IT_0368 = IT_0360 + IT_0367;
    const ccomplex_t IT_0369 = (0 + _Complex_I*1)*IT_0368;
    const ccomplex_t IT_0370 = 0.5*IT_0369;
    const ccomplex_t IT_0371 = conj(N_u3)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0372 = IT_0229*IT_0371;
    const ccomplex_t IT_0373 = 2*IT_0372;
    const ccomplex_t IT_0374 = conj(N_B3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0375 = IT_0233*IT_0374;
    const ccomplex_t IT_0376 = 1.4142135623731*IT_0375;
    const ccomplex_t IT_0377 = conj(N_W3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0378 = IT_0229*IT_0377;
    const ccomplex_t IT_0379 = 1.4142135623731*IT_0378;
    const ccomplex_t IT_0380 = IT_0376 + IT_0379;
    const ccomplex_t IT_0381 = IT_0373 + IT_0380;
    const ccomplex_t IT_0382 = (0 + _Complex_I*1)*IT_0381;
    const ccomplex_t IT_0383 = 0.5*IT_0382;
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = m_C_1*m_N_3;
    const ccomplex_t IT_0386 = pow(m_N_3, 2);
    const ccomplex_t IT_0387 = cpow(IT_0012 + -IT_0247 + IT_0386, 2);
    const ccomplex_t IT_0388 = -IT_0386 + 0.25*IT_0001*(IT_0012 + -IT_0247 +
       IT_0386)*(IT_0012 + IT_0247 + -IT_0386) + 0.25*IT_0001*IT_0387;
    const ccomplex_t IT_0389 = 2*IT_0388;
    const ccomplex_t IT_0390 = cabs(-IT_0386 + 0.25*IT_0001*IT_0387);
    const ccomplex_t IT_0391 = cpow(IT_0390, 0.5);
    const ccomplex_t IT_0392 = ( ((m_Hp
) >= (m_C_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0370)*(IT_0384*IT_0385 + (-0.5)*IT_0370*IT_0389)
       + conj(IT_0384)*(IT_0370*IT_0385 + (-0.5)*IT_0384*IT_0389))*IT_0391
) : 
      (0
));
    const ccomplex_t IT_0393 = conj(N_u3)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0394 = IT_0229*IT_0393;
    const ccomplex_t IT_0395 = 2*IT_0394;
    const ccomplex_t IT_0396 = conj(N_B3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0397 = IT_0233*IT_0396;
    const ccomplex_t IT_0398 = 1.4142135623731*IT_0397;
    const ccomplex_t IT_0399 = conj(N_W3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0400 = IT_0229*IT_0399;
    const ccomplex_t IT_0401 = 1.4142135623731*IT_0400;
    const ccomplex_t IT_0402 = IT_0398 + IT_0401;
    const ccomplex_t IT_0403 = IT_0395 + IT_0402;
    const ccomplex_t IT_0404 = (0 + _Complex_I*1)*IT_0403;
    const ccomplex_t IT_0405 = 0.5*IT_0404;
    const ccomplex_t IT_0406 = -IT_0405;
    const ccomplex_t IT_0407 = N_W3*U_d2*e_em;
    const ccomplex_t IT_0408 = IT_0214*IT_0407;
    const ccomplex_t IT_0409 = 1.4142135623731*IT_0408;
    const ccomplex_t IT_0410 = N_B3*U_d2*e_em;
    const ccomplex_t IT_0411 = IT_0218*IT_0410;
    const ccomplex_t IT_0412 = 1.4142135623731*IT_0411;
    const ccomplex_t IT_0413 = IT_0409 + IT_0412;
    const ccomplex_t IT_0414 = N_d3*e_em*U_Wm2;
    const ccomplex_t IT_0415 = IT_0214*IT_0414;
    const ccomplex_t IT_0416 = (-2)*IT_0415;
    const ccomplex_t IT_0417 = IT_0413 + IT_0416;
    const ccomplex_t IT_0418 = (0 + _Complex_I*1)*IT_0417;
    const ccomplex_t IT_0419 = 0.5*IT_0418;
    const ccomplex_t IT_0420 = m_C_2*m_N_3;
    const ccomplex_t IT_0421 = cpow(IT_0012 + -IT_0282 + IT_0386, 2);
    const ccomplex_t IT_0422 = -IT_0386 + 0.25*IT_0001*(IT_0012 + -IT_0282 +
       IT_0386)*(IT_0012 + IT_0282 + -IT_0386) + 0.25*IT_0001*IT_0421;
    const ccomplex_t IT_0423 = 2*IT_0422;
    const ccomplex_t IT_0424 = cabs(-IT_0386 + 0.25*IT_0001*IT_0421);
    const ccomplex_t IT_0425 = cpow(IT_0424, 0.5);
    const ccomplex_t IT_0426 = ( ((m_Hp
) >= (m_C_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0406)*(IT_0419*IT_0420 + (-0.5)*IT_0406*IT_0423)
       + conj(IT_0419)*(IT_0406*IT_0420 + (-0.5)*IT_0419*IT_0423))*IT_0425
) : 
      (0
));
    const ccomplex_t IT_0427 = N_d4*e_em*U_Wm1;
    const ccomplex_t IT_0428 = IT_0214*IT_0427;
    const ccomplex_t IT_0429 = (-2)*IT_0428;
    const ccomplex_t IT_0430 = N_B4*U_d1*e_em;
    const ccomplex_t IT_0431 = IT_0218*IT_0430;
    const ccomplex_t IT_0432 = 1.4142135623731*IT_0431;
    const ccomplex_t IT_0433 = N_W4*U_d1*e_em;
    const ccomplex_t IT_0434 = IT_0214*IT_0433;
    const ccomplex_t IT_0435 = 1.4142135623731*IT_0434;
    const ccomplex_t IT_0436 = IT_0432 + IT_0435;
    const ccomplex_t IT_0437 = IT_0429 + IT_0436;
    const ccomplex_t IT_0438 = (0 + _Complex_I*1)*IT_0437;
    const ccomplex_t IT_0439 = 0.5*IT_0438;
    const ccomplex_t IT_0440 = conj(N_u4)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0441 = IT_0229*IT_0440;
    const ccomplex_t IT_0442 = 2*IT_0441;
    const ccomplex_t IT_0443 = conj(N_B4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0444 = IT_0233*IT_0443;
    const ccomplex_t IT_0445 = 1.4142135623731*IT_0444;
    const ccomplex_t IT_0446 = conj(N_W4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0447 = IT_0229*IT_0446;
    const ccomplex_t IT_0448 = 1.4142135623731*IT_0447;
    const ccomplex_t IT_0449 = IT_0445 + IT_0448;
    const ccomplex_t IT_0450 = IT_0442 + IT_0449;
    const ccomplex_t IT_0451 = (0 + _Complex_I*1)*IT_0450;
    const ccomplex_t IT_0452 = 0.5*IT_0451;
    const ccomplex_t IT_0453 = -IT_0452;
    const ccomplex_t IT_0454 = m_C_1*m_N_4;
    const ccomplex_t IT_0455 = pow(m_N_4, 2);
    const ccomplex_t IT_0456 = cpow(IT_0012 + -IT_0247 + IT_0455, 2);
    const ccomplex_t IT_0457 = -IT_0455 + 0.25*IT_0001*(IT_0012 + -IT_0247 +
       IT_0455)*(IT_0012 + IT_0247 + -IT_0455) + 0.25*IT_0001*IT_0456;
    const ccomplex_t IT_0458 = 2*IT_0457;
    const ccomplex_t IT_0459 = cabs(-IT_0455 + 0.25*IT_0001*IT_0456);
    const ccomplex_t IT_0460 = cpow(IT_0459, 0.5);
    const ccomplex_t IT_0461 = ( ((m_Hp
) >= (m_C_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0439)*(IT_0453*IT_0454 + (-0.5)*IT_0439*IT_0458)
       + conj(IT_0453)*(IT_0439*IT_0454 + (-0.5)*IT_0453*IT_0458))*IT_0460
) : 
      (0
));
    const ccomplex_t IT_0462 = conj(N_B4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0463 = IT_0233*IT_0462;
    const ccomplex_t IT_0464 = 1.4142135623731*IT_0463;
    const ccomplex_t IT_0465 = conj(N_W4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0466 = IT_0229*IT_0465;
    const ccomplex_t IT_0467 = 1.4142135623731*IT_0466;
    const ccomplex_t IT_0468 = IT_0464 + IT_0467;
    const ccomplex_t IT_0469 = conj(N_u4)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0470 = IT_0229*IT_0469;
    const ccomplex_t IT_0471 = 2*IT_0470;
    const ccomplex_t IT_0472 = IT_0468 + IT_0471;
    const ccomplex_t IT_0473 = (0 + _Complex_I*1)*IT_0472;
    const ccomplex_t IT_0474 = 0.5*IT_0473;
    const ccomplex_t IT_0475 = -IT_0474;
    const ccomplex_t IT_0476 = N_d4*e_em*U_Wm2;
    const ccomplex_t IT_0477 = IT_0214*IT_0476;
    const ccomplex_t IT_0478 = (-2)*IT_0477;
    const ccomplex_t IT_0479 = N_B4*U_d2*e_em;
    const ccomplex_t IT_0480 = IT_0218*IT_0479;
    const ccomplex_t IT_0481 = 1.4142135623731*IT_0480;
    const ccomplex_t IT_0482 = N_W4*U_d2*e_em;
    const ccomplex_t IT_0483 = IT_0214*IT_0482;
    const ccomplex_t IT_0484 = 1.4142135623731*IT_0483;
    const ccomplex_t IT_0485 = IT_0481 + IT_0484;
    const ccomplex_t IT_0486 = IT_0478 + IT_0485;
    const ccomplex_t IT_0487 = (0 + _Complex_I*1)*IT_0486;
    const ccomplex_t IT_0488 = 0.5*IT_0487;
    const ccomplex_t IT_0489 = m_C_2*m_N_4;
    const ccomplex_t IT_0490 = cpow(IT_0012 + -IT_0282 + IT_0455, 2);
    const ccomplex_t IT_0491 = -IT_0455 + 0.25*IT_0001*(IT_0012 + -IT_0282 +
       IT_0455)*(IT_0012 + IT_0282 + -IT_0455) + 0.25*IT_0001*IT_0490;
    const ccomplex_t IT_0492 = 2*IT_0491;
    const ccomplex_t IT_0493 = cabs(-IT_0455 + 0.25*IT_0001*IT_0490);
    const ccomplex_t IT_0494 = cpow(IT_0493, 0.5);
    const ccomplex_t IT_0495 = ( ((m_Hp
) >= (m_C_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0475)*(IT_0488*IT_0489 + (-0.5)*IT_0475*IT_0492)
       + conj(IT_0488)*(IT_0475*IT_0489 + (-0.5)*IT_0488*IT_0492))*IT_0494
) : 
      (0
));
    const ccomplex_t IT_0496 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0497 = 0.5*IT_0496;
    const ccomplex_t IT_0498 = cpow(IT_0497, 2);
    const ccomplex_t IT_0499 = pow(m_tau, 2);
    const ccomplex_t IT_0500 = cpow(IT_0012 + -IT_0499, 2);
    const ccomplex_t IT_0501 = 0.25*IT_0001*(IT_0012 + IT_0499)*(IT_0012 + 
      -IT_0499) + 0.25*IT_0001*IT_0500;
    const ccomplex_t IT_0502 = cabs(0.25*IT_0001*IT_0500);
    const ccomplex_t IT_0503 = cpow(IT_0502, 0.5);
    const ccomplex_t IT_0504 = ( ((m_Hp
) >= (m_tau
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0498*IT_0501*IT_0503
) : (0
));
    const ccomplex_t IT_0505 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0007*(m_W*IT_0003*IT_0005*U_sb_00 + (-0.5)*IT_0002*(IT_0003*IT_0048
      *IT_0123*U_sb_00 + m_b*(m_b*IT_0004*IT_0005*U_sb_00 + mu_h*U_sb_10)));
    const ccomplex_t IT_0506 = -IT_0505;
    const ccomplex_t IT_0507 = pow(m_sc_L, 2);
    const ccomplex_t IT_0508 = pow(m_sb_1, 2);
    const ccomplex_t IT_0509 = cpow(IT_0012 + IT_0507 + -IT_0508, 2);
    const ccomplex_t IT_0510 = cabs(-IT_0507 + 0.25*IT_0001*IT_0509);
    const ccomplex_t IT_0511 = cpow(IT_0510, 0.5);
    const ccomplex_t IT_0512 = ( ((m_Hp
) >= (m_sb_1 + m_sc_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0506*conj(IT_0506)*IT_0511
) : (0
));
    const ccomplex_t IT_0513 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0007*(m_W*IT_0003*IT_0005*U_sb_01 + (-0.5)*IT_0002*(IT_0003*IT_0048
      *IT_0123*U_sb_01 + m_b*(m_b*IT_0004*IT_0005*U_sb_01 + mu_h*U_sb_11)));
    const ccomplex_t IT_0514 = -IT_0513;
    const ccomplex_t IT_0515 = pow(m_sb_2, 2);
    const ccomplex_t IT_0516 = cpow(IT_0012 + IT_0507 + -IT_0515, 2);
    const ccomplex_t IT_0517 = cabs(-IT_0507 + 0.25*IT_0001*IT_0516);
    const ccomplex_t IT_0518 = cpow(IT_0517, 0.5);
    const ccomplex_t IT_0519 = ( ((m_Hp
) >= (m_sb_2 + m_sc_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0514*conj(IT_0514)*IT_0518
) : (0
));
    const ccomplex_t IT_0520 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0002*IT_0007*(mu_h*U_sb_00 + m_b*IT_0004*IT_0048*U_sb_10);
    const ccomplex_t IT_0521 = 0.5*IT_0520;
    const ccomplex_t IT_0522 = pow(m_sc_R, 2);
    const ccomplex_t IT_0523 = cpow(IT_0012 + -IT_0508 + IT_0522, 2);
    const ccomplex_t IT_0524 = cabs(-IT_0522 + 0.25*IT_0001*IT_0523);
    const ccomplex_t IT_0525 = cpow(IT_0524, 0.5);
    const ccomplex_t IT_0526 = ( ((m_Hp
) >= (m_sb_1 + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0521*conj(IT_0521)*IT_0525
) : (0
));
    const ccomplex_t IT_0527 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0002*IT_0007*(mu_h*U_sb_01 + m_b*IT_0004*IT_0048*U_sb_11);
    const ccomplex_t IT_0528 = 0.5*IT_0527;
    const ccomplex_t IT_0529 = cpow(IT_0012 + -IT_0515 + IT_0522, 2);
    const ccomplex_t IT_0530 = cabs(-IT_0522 + 0.25*IT_0001*IT_0529);
    const ccomplex_t IT_0531 = cpow(IT_0530, 0.5);
    const ccomplex_t IT_0532 = ( ((m_Hp
) >= (m_sb_2 + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0528*conj(IT_0528)*IT_0531
) : (0
));
    const ccomplex_t IT_0533 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0007*(IT_0002*IT_0004*IT_0005*IT_0054 + (-2)*IT_0003*(m_W*IT_0005 + (
      -0.5)*IT_0002*IT_0048*IT_0123));
    const ccomplex_t IT_0534 = 0.5*IT_0533;
    const ccomplex_t IT_0535 = pow(m_sd_L, 2);
    const ccomplex_t IT_0536 = cpow(IT_0012 + IT_0507 + -IT_0535, 2);
    const ccomplex_t IT_0537 = cabs(-IT_0507 + 0.25*IT_0001*IT_0536);
    const ccomplex_t IT_0538 = cpow(IT_0537, 0.5);
    const ccomplex_t IT_0539 = ( ((m_Hp
) >= (m_sc_L + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0534*conj(IT_0534)*IT_0538
) : (0
));
    const ccomplex_t IT_0540 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0541 = 0.5*IT_0540;
    const ccomplex_t IT_0542 = cpow(IT_0012 + IT_0522 + -IT_0535, 2);
    const ccomplex_t IT_0543 = cabs(-IT_0522 + 0.25*IT_0001*IT_0542);
    const ccomplex_t IT_0544 = cpow(IT_0543, 0.5);
    const ccomplex_t IT_0545 = ( ((m_Hp
) >= (m_sc_R + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0541*conj(IT_0541)*IT_0544
) : (0
));
    const ccomplex_t IT_0546 = (0 + _Complex_I*1.4142135623731)*m_d*V_cd*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0547 = 0.5*IT_0546;
    const ccomplex_t IT_0548 = pow(m_sd_R, 2);
    const ccomplex_t IT_0549 = cpow(IT_0012 + IT_0507 + -IT_0548, 2);
    const ccomplex_t IT_0550 = cabs(-IT_0507 + 0.25*IT_0001*IT_0549);
    const ccomplex_t IT_0551 = cpow(IT_0550, 0.5);
    const ccomplex_t IT_0552 = ( ((m_Hp
) >= (m_sc_L + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0547*conj(IT_0547)*IT_0551
) : (0
));
    const ccomplex_t IT_0553 = (0 + _Complex_I*1.4142135623731)*m_c*m_d*V_cd
      *e_em*IT_0002*IT_0004*IT_0007*IT_0048;
    const ccomplex_t IT_0554 = 0.5*IT_0553;
    const ccomplex_t IT_0555 = cpow(IT_0012 + IT_0522 + -IT_0548, 2);
    const ccomplex_t IT_0556 = cabs(-IT_0522 + 0.25*IT_0001*IT_0555);
    const ccomplex_t IT_0557 = cpow(IT_0556, 0.5);
    const ccomplex_t IT_0558 = ( ((m_Hp
) >= (m_sc_R + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0554*conj(IT_0554)*IT_0557
) : (0
));
    const ccomplex_t IT_0559 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0007*(IT_0002*IT_0004*IT_0005*IT_0067 + (-2)*IT_0003*(m_W*IT_0005 + (
      -0.5)*IT_0002*IT_0048*IT_0123));
    const ccomplex_t IT_0560 = 0.5*IT_0559;
    const ccomplex_t IT_0561 = pow(m_ss_L, 2);
    const ccomplex_t IT_0562 = cpow(IT_0012 + IT_0507 + -IT_0561, 2);
    const ccomplex_t IT_0563 = cabs(-IT_0507 + 0.25*IT_0001*IT_0562);
    const ccomplex_t IT_0564 = cpow(IT_0563, 0.5);
    const ccomplex_t IT_0565 = ( ((m_Hp
) >= (m_sc_L + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0560*conj(IT_0560)*IT_0564
) : (0
));
    const ccomplex_t IT_0566 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0567 = 0.5*IT_0566;
    const ccomplex_t IT_0568 = cpow(IT_0012 + IT_0522 + -IT_0561, 2);
    const ccomplex_t IT_0569 = cabs(-IT_0522 + 0.25*IT_0001*IT_0568);
    const ccomplex_t IT_0570 = cpow(IT_0569, 0.5);
    const ccomplex_t IT_0571 = ( ((m_Hp
) >= (m_sc_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0567*conj(IT_0567)*IT_0570
) : (0
));
    const ccomplex_t IT_0572 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0573 = 0.5*IT_0572;
    const ccomplex_t IT_0574 = pow(m_ss_R, 2);
    const ccomplex_t IT_0575 = cpow(IT_0012 + IT_0507 + -IT_0574, 2);
    const ccomplex_t IT_0576 = cabs(-IT_0507 + 0.25*IT_0001*IT_0575);
    const ccomplex_t IT_0577 = cpow(IT_0576, 0.5);
    const ccomplex_t IT_0578 = ( ((m_Hp
) >= (m_sc_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0573*conj(IT_0573)*IT_0577
) : (0
));
    const ccomplex_t IT_0579 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*V_cs
      *e_em*IT_0002*IT_0004*IT_0007*IT_0048;
    const ccomplex_t IT_0580 = 0.5*IT_0579;
    const ccomplex_t IT_0581 = cpow(IT_0012 + IT_0522 + -IT_0574, 2);
    const ccomplex_t IT_0582 = cabs(-IT_0522 + 0.25*IT_0001*IT_0581);
    const ccomplex_t IT_0583 = cpow(IT_0582, 0.5);
    const ccomplex_t IT_0584 = ( ((m_Hp
) >= (m_sc_R + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0580*conj(IT_0580)*IT_0583
) : (0
));
    const ccomplex_t IT_0585 = m_t*U_sb_00*(m_t*IT_0003*IT_0048*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0004*IT_0005*U_sb_00*conj(U_st_00) +
       U_sb_10*(mu_h*conj(U_st_00) + m_t*IT_0004*IT_0048*conj(U_st_10)));
    const ccomplex_t IT_0586 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_10
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002*IT_0585 +
       m_W*IT_0003*IT_0005*U_sb_00*conj(U_st_00)) + A_t*IT_0003*U_sb_00*conj
      (U_st_10));
    const ccomplex_t IT_0587 = pow(m_st_1, 2);
    const ccomplex_t IT_0588 = cpow(IT_0012 + -IT_0508 + IT_0587, 2);
    const ccomplex_t IT_0589 = cabs(-IT_0587 + 0.25*IT_0001*IT_0588);
    const ccomplex_t IT_0590 = cpow(IT_0589, 0.5);
    const ccomplex_t IT_0591 = ( ((m_Hp
) >= (m_sb_1 + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0586*conj(IT_0586)*IT_0590
) : (0
));
    const ccomplex_t IT_0592 = m_t*U_sb_01*(m_t*IT_0003*IT_0048*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0004*IT_0005*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0004*IT_0048*conj(U_st_10)));
    const ccomplex_t IT_0593 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002*IT_0592 +
       m_W*IT_0003*IT_0005*U_sb_01*conj(U_st_00)) + A_t*IT_0003*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0594 = cpow(IT_0012 + -IT_0515 + IT_0587, 2);
    const ccomplex_t IT_0595 = cabs(-IT_0587 + 0.25*IT_0001*IT_0594);
    const ccomplex_t IT_0596 = cpow(IT_0595, 0.5);
    const ccomplex_t IT_0597 = ( ((m_Hp
) >= (m_sb_2 + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0593*conj(IT_0593)*IT_0596
) : (0
));
    const ccomplex_t IT_0598 = A_t*IT_0003*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0007*(m_W*IT_0003*IT_0005*conj(U_st_00) + (-0.5)*IT_0002*(IT_0004
      *IT_0005*IT_0054*conj(U_st_00) + m_t*(m_t*IT_0003*IT_0048*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0599 = (0 + _Complex_I*1)*V_td*IT_0598;
    const ccomplex_t IT_0600 = cpow(IT_0012 + IT_0535 + -IT_0587, 2);
    const ccomplex_t IT_0601 = cabs(-IT_0535 + 0.25*IT_0001*IT_0600);
    const ccomplex_t IT_0602 = cpow(IT_0601, 0.5);
    const ccomplex_t IT_0603 = ( ((m_Hp
) >= (m_sd_L + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0599*conj(IT_0599)*IT_0602
) : (0
));
    const ccomplex_t IT_0604 = (0 + _Complex_I*1)*V_td*(A_b*IT_0005*conj
      (U_st_00) + 0.707106781186548*m_d*e_em*IT_0002*IT_0007*(mu_h*conj(U_st_00)
       + m_t*IT_0004*IT_0048*conj(U_st_10)));
    const ccomplex_t IT_0605 = cpow(IT_0012 + IT_0548 + -IT_0587, 2);
    const ccomplex_t IT_0606 = cabs(-IT_0548 + 0.25*IT_0001*IT_0605);
    const ccomplex_t IT_0607 = cpow(IT_0606, 0.5);
    const ccomplex_t IT_0608 = ( ((m_Hp
) >= (m_sd_R + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0604*conj(IT_0604)*IT_0607
) : (0
));
    const ccomplex_t IT_0609 = A_t*IT_0003*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0007*(m_W*IT_0003*IT_0005*conj(U_st_00) + (-0.5)*IT_0002*(IT_0004
      *IT_0005*IT_0067*conj(U_st_00) + m_t*(m_t*IT_0003*IT_0048*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0610 = (0 + _Complex_I*1)*V_ts*IT_0609;
    const ccomplex_t IT_0611 = cpow(IT_0012 + IT_0561 + -IT_0587, 2);
    const ccomplex_t IT_0612 = cabs(-IT_0561 + 0.25*IT_0001*IT_0611);
    const ccomplex_t IT_0613 = cpow(IT_0612, 0.5);
    const ccomplex_t IT_0614 = ( ((m_Hp
) >= (m_ss_L + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0610*conj(IT_0610)*IT_0613
) : (0
));
    const ccomplex_t IT_0615 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0005*conj
      (U_st_00) + 0.707106781186548*m_s*e_em*IT_0002*IT_0007*(mu_h*conj(U_st_00)
       + m_t*IT_0004*IT_0048*conj(U_st_10)));
    const ccomplex_t IT_0616 = cpow(IT_0012 + IT_0574 + -IT_0587, 2);
    const ccomplex_t IT_0617 = cabs(-IT_0574 + 0.25*IT_0001*IT_0616);
    const ccomplex_t IT_0618 = cpow(IT_0617, 0.5);
    const ccomplex_t IT_0619 = ( ((m_Hp
) >= (m_ss_R + m_st_1
)) ? (0.375
      *IT_0000*IT_0001*IT_0615*conj(IT_0615)*IT_0618
) : (0
));
    const ccomplex_t IT_0620 = m_t*U_sb_00*(m_t*IT_0003*IT_0048*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0004*IT_0005*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0004*IT_0048*conj(U_st_11)));
    const ccomplex_t IT_0621 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002*IT_0620 +
       m_W*IT_0003*IT_0005*U_sb_00*conj(U_st_01)) + A_t*IT_0003*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0622 = pow(m_st_2, 2);
    const ccomplex_t IT_0623 = cpow(IT_0012 + -IT_0508 + IT_0622, 2);
    const ccomplex_t IT_0624 = cabs(-IT_0622 + 0.25*IT_0001*IT_0623);
    const ccomplex_t IT_0625 = cpow(IT_0624, 0.5);
    const ccomplex_t IT_0626 = ( ((m_Hp
) >= (m_sb_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0621*conj(IT_0621)*IT_0625
) : (0
));
    const ccomplex_t IT_0627 = m_t*U_sb_01*(m_t*IT_0003*IT_0048*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0004*IT_0005*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0004*IT_0048*conj(U_st_11)));
    const ccomplex_t IT_0628 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0007*((-0.5)*IT_0002*IT_0627 +
       m_W*IT_0003*IT_0005*U_sb_01*conj(U_st_01)) + A_t*IT_0003*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0629 = cpow(IT_0012 + -IT_0515 + IT_0622, 2);
    const ccomplex_t IT_0630 = cabs(-IT_0622 + 0.25*IT_0001*IT_0629);
    const ccomplex_t IT_0631 = cpow(IT_0630, 0.5);
    const ccomplex_t IT_0632 = ( ((m_Hp
) >= (m_sb_2 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0628*conj(IT_0628)*IT_0631
) : (0
));
    const ccomplex_t IT_0633 = A_t*IT_0003*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0007*(m_W*IT_0003*IT_0005*conj(U_st_01) + (-0.5)*IT_0002*(IT_0004
      *IT_0005*IT_0054*conj(U_st_01) + m_t*(m_t*IT_0003*IT_0048*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0634 = (0 + _Complex_I*1)*V_td*IT_0633;
    const ccomplex_t IT_0635 = cpow(IT_0012 + IT_0535 + -IT_0622, 2);
    const ccomplex_t IT_0636 = cabs(-IT_0535 + 0.25*IT_0001*IT_0635);
    const ccomplex_t IT_0637 = cpow(IT_0636, 0.5);
    const ccomplex_t IT_0638 = ( ((m_Hp
) >= (m_sd_L + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0634*conj(IT_0634)*IT_0637
) : (0
));
    const ccomplex_t IT_0639 = (0 + _Complex_I*1)*V_td*(A_b*IT_0005*conj
      (U_st_01) + 0.707106781186548*m_d*e_em*IT_0002*IT_0007*(mu_h*conj(U_st_01)
       + m_t*IT_0004*IT_0048*conj(U_st_11)));
    const ccomplex_t IT_0640 = cpow(IT_0012 + IT_0548 + -IT_0622, 2);
    const ccomplex_t IT_0641 = cabs(-IT_0548 + 0.25*IT_0001*IT_0640);
    const ccomplex_t IT_0642 = cpow(IT_0641, 0.5);
    const ccomplex_t IT_0643 = ( ((m_Hp
) >= (m_sd_R + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0639*conj(IT_0639)*IT_0642
) : (0
));
    const ccomplex_t IT_0644 = A_t*IT_0003*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0007*(m_W*IT_0003*IT_0005*conj(U_st_01) + (-0.5)*IT_0002*(IT_0004
      *IT_0005*IT_0067*conj(U_st_01) + m_t*(m_t*IT_0003*IT_0048*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0645 = (0 + _Complex_I*1)*V_ts*IT_0644;
    const ccomplex_t IT_0646 = cpow(IT_0012 + IT_0561 + -IT_0622, 2);
    const ccomplex_t IT_0647 = cabs(-IT_0561 + 0.25*IT_0001*IT_0646);
    const ccomplex_t IT_0648 = cpow(IT_0647, 0.5);
    const ccomplex_t IT_0649 = ( ((m_Hp
) >= (m_ss_L + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0645*conj(IT_0645)*IT_0648
) : (0
));
    const ccomplex_t IT_0650 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0005*conj
      (U_st_01) + 0.707106781186548*m_s*e_em*IT_0002*IT_0007*(mu_h*conj(U_st_01)
       + m_t*IT_0004*IT_0048*conj(U_st_11)));
    const ccomplex_t IT_0651 = cpow(IT_0012 + IT_0574 + -IT_0622, 2);
    const ccomplex_t IT_0652 = cabs(-IT_0574 + 0.25*IT_0001*IT_0651);
    const ccomplex_t IT_0653 = cpow(IT_0652, 0.5);
    const ccomplex_t IT_0654 = ( ((m_Hp
) >= (m_ss_R + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0650*conj(IT_0650)*IT_0653
) : (0
));
    const ccomplex_t IT_0655 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007
      *IT_0167*(m_W*IT_0003*IT_0005*U_sb_00 + (-0.5)*IT_0002*(IT_0003*IT_0048
      *IT_0171*U_sb_00 + m_b*(m_b*IT_0004*IT_0005*U_sb_00 + mu_h*U_sb_10)))
      *V_ub_mod;
    const ccomplex_t IT_0656 = -IT_0655;
    const ccomplex_t IT_0657 = pow(m_su_L, 2);
    const ccomplex_t IT_0658 = cpow(IT_0012 + -IT_0508 + IT_0657, 2);
    const ccomplex_t IT_0659 = cabs(-IT_0657 + 0.25*IT_0001*IT_0658);
    const ccomplex_t IT_0660 = cpow(IT_0659, 0.5);
    const ccomplex_t IT_0661 = ( ((m_Hp
) >= (m_sb_1 + m_su_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0656*conj(IT_0656)*IT_0660
) : (0
));
    const ccomplex_t IT_0662 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007
      *IT_0167*(m_W*IT_0003*IT_0005*U_sb_01 + (-0.5)*IT_0002*(IT_0003*IT_0048
      *IT_0171*U_sb_01 + m_b*(m_b*IT_0004*IT_0005*U_sb_01 + mu_h*U_sb_11)))
      *V_ub_mod;
    const ccomplex_t IT_0663 = -IT_0662;
    const ccomplex_t IT_0664 = cpow(IT_0012 + -IT_0515 + IT_0657, 2);
    const ccomplex_t IT_0665 = cabs(-IT_0657 + 0.25*IT_0001*IT_0664);
    const ccomplex_t IT_0666 = cpow(IT_0665, 0.5);
    const ccomplex_t IT_0667 = ( ((m_Hp
) >= (m_sb_2 + m_su_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0663*conj(IT_0663)*IT_0666
) : (0
));
    const ccomplex_t IT_0668 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0007*(IT_0002*IT_0004*IT_0005*IT_0054 + (-2)*IT_0003*(m_W*IT_0005 + (
      -0.5)*IT_0002*IT_0048*IT_0171));
    const ccomplex_t IT_0669 = 0.5*IT_0668;
    const ccomplex_t IT_0670 = cpow(IT_0669, 2);
    const ccomplex_t IT_0671 = cpow(IT_0012 + -IT_0535 + IT_0657, 2);
    const ccomplex_t IT_0672 = cabs(-IT_0657 + 0.25*IT_0001*IT_0671);
    const ccomplex_t IT_0673 = cpow(IT_0672, 0.5);
    const ccomplex_t IT_0674 = ( ((m_Hp
) >= (m_sd_L + m_su_L
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0670*IT_0673
) : (0
));
    const ccomplex_t IT_0675 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0676 = 0.5*IT_0675;
    const ccomplex_t IT_0677 = cpow(IT_0676, 2);
    const ccomplex_t IT_0678 = cpow(IT_0012 + -IT_0548 + IT_0657, 2);
    const ccomplex_t IT_0679 = cabs(-IT_0657 + 0.25*IT_0001*IT_0678);
    const ccomplex_t IT_0680 = cpow(IT_0679, 0.5);
    const ccomplex_t IT_0681 = ( ((m_Hp
) >= (m_sd_R + m_su_L
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0677*IT_0680
) : (0
));
    const ccomplex_t IT_0682 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0007*(IT_0002*IT_0004*IT_0005*IT_0067 + (-2)*IT_0003*(m_W*IT_0005 + (
      -0.5)*IT_0002*IT_0048*IT_0171));
    const ccomplex_t IT_0683 = 0.5*IT_0682;
    const ccomplex_t IT_0684 = cpow(IT_0683, 2);
    const ccomplex_t IT_0685 = cpow(IT_0012 + -IT_0561 + IT_0657, 2);
    const ccomplex_t IT_0686 = cabs(-IT_0657 + 0.25*IT_0001*IT_0685);
    const ccomplex_t IT_0687 = cpow(IT_0686, 0.5);
    const ccomplex_t IT_0688 = ( ((m_Hp
) >= (m_ss_L + m_su_L
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0684*IT_0687
) : (0
));
    const ccomplex_t IT_0689 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0690 = 0.5*IT_0689;
    const ccomplex_t IT_0691 = cpow(IT_0690, 2);
    const ccomplex_t IT_0692 = cpow(IT_0012 + -IT_0574 + IT_0657, 2);
    const ccomplex_t IT_0693 = cabs(-IT_0657 + 0.25*IT_0001*IT_0692);
    const ccomplex_t IT_0694 = cpow(IT_0693, 0.5);
    const ccomplex_t IT_0695 = ( ((m_Hp
) >= (m_ss_R + m_su_L
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0691*IT_0694
) : (0
));
    const ccomplex_t IT_0696 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0002*IT_0007*IT_0167*(mu_h*U_sb_00 + m_b*IT_0004*IT_0048*U_sb_10)
      *V_ub_mod;
    const ccomplex_t IT_0697 = 0.5*IT_0696;
    const ccomplex_t IT_0698 = pow(m_su_R, 2);
    const ccomplex_t IT_0699 = cpow(IT_0012 + -IT_0508 + IT_0698, 2);
    const ccomplex_t IT_0700 = cabs(-IT_0698 + 0.25*IT_0001*IT_0699);
    const ccomplex_t IT_0701 = cpow(IT_0700, 0.5);
    const ccomplex_t IT_0702 = ( ((m_Hp
) >= (m_sb_1 + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0697*conj(IT_0697)*IT_0701
) : (0
));
    const ccomplex_t IT_0703 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0002*IT_0007*IT_0167*(mu_h*U_sb_01 + m_b*IT_0004*IT_0048*U_sb_11)
      *V_ub_mod;
    const ccomplex_t IT_0704 = 0.5*IT_0703;
    const ccomplex_t IT_0705 = cpow(IT_0012 + -IT_0515 + IT_0698, 2);
    const ccomplex_t IT_0706 = cabs(-IT_0698 + 0.25*IT_0001*IT_0705);
    const ccomplex_t IT_0707 = cpow(IT_0706, 0.5);
    const ccomplex_t IT_0708 = ( ((m_Hp
) >= (m_sb_2 + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0704*conj(IT_0704)*IT_0707
) : (0
));
    const ccomplex_t IT_0709 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0710 = 0.5*IT_0709;
    const ccomplex_t IT_0711 = cpow(IT_0710, 2);
    const ccomplex_t IT_0712 = cpow(IT_0012 + -IT_0535 + IT_0698, 2);
    const ccomplex_t IT_0713 = cabs(-IT_0698 + 0.25*IT_0001*IT_0712);
    const ccomplex_t IT_0714 = cpow(IT_0713, 0.5);
    const ccomplex_t IT_0715 = ( ((m_Hp
) >= (m_sd_L + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0711*IT_0714
) : (0
));
    const ccomplex_t IT_0716 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0002*IT_0004*IT_0007*IT_0048;
    const ccomplex_t IT_0717 = 0.5*IT_0716;
    const ccomplex_t IT_0718 = cpow(IT_0717, 2);
    const ccomplex_t IT_0719 = cpow(IT_0012 + -IT_0548 + IT_0698, 2);
    const ccomplex_t IT_0720 = cabs(-IT_0698 + 0.25*IT_0001*IT_0719);
    const ccomplex_t IT_0721 = cpow(IT_0720, 0.5);
    const ccomplex_t IT_0722 = ( ((m_Hp
) >= (m_sd_R + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0718*IT_0721
) : (0
));
    const ccomplex_t IT_0723 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0724 = 0.5*IT_0723;
    const ccomplex_t IT_0725 = cpow(IT_0724, 2);
    const ccomplex_t IT_0726 = cpow(IT_0012 + -IT_0561 + IT_0698, 2);
    const ccomplex_t IT_0727 = cabs(-IT_0698 + 0.25*IT_0001*IT_0726);
    const ccomplex_t IT_0728 = cpow(IT_0727, 0.5);
    const ccomplex_t IT_0729 = ( ((m_Hp
) >= (m_ss_L + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0725*IT_0728
) : (0
));
    const ccomplex_t IT_0730 = (0 + _Complex_I*1.4142135623731)*m_s*m_u*V_us
      *e_em*IT_0002*IT_0004*IT_0007*IT_0048;
    const ccomplex_t IT_0731 = 0.5*IT_0730;
    const ccomplex_t IT_0732 = cpow(IT_0731, 2);
    const ccomplex_t IT_0733 = cpow(IT_0012 + -IT_0574 + IT_0698, 2);
    const ccomplex_t IT_0734 = cabs(-IT_0698 + 0.25*IT_0001*IT_0733);
    const ccomplex_t IT_0735 = cpow(IT_0734, 0.5);
    const ccomplex_t IT_0736 = ( ((m_Hp
) >= (m_ss_R + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0732*IT_0735
) : (0
));
    const ccomplex_t IT_0737 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *IT_0007*(m_W*IT_0003 + (-0.5)*IT_0002*IT_0004*IT_0011);
    const ccomplex_t IT_0738 = -IT_0737;
    const ccomplex_t IT_0739 = cpow(IT_0738, 2);
    const ccomplex_t IT_0740 = pow(m_se_L, 2);
    const ccomplex_t IT_0741 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0742 = cpow(IT_0012 + IT_0740 + -IT_0741, 2);
    const ccomplex_t IT_0743 = cabs(-IT_0740 + 0.25*IT_0001*IT_0742);
    const ccomplex_t IT_0744 = cpow(IT_0743, 0.5);
    const ccomplex_t IT_0745 = ( ((m_Hp
) >= (m_se_L + m_snu_mu
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0739*IT_0744
) : (0
));
    const ccomplex_t IT_0746 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0002*IT_0007;
    const ccomplex_t IT_0747 = 0.5*IT_0746;
    const ccomplex_t IT_0748 = cpow(IT_0747, 2);
    const ccomplex_t IT_0749 = pow(m_se_R, 2);
    const ccomplex_t IT_0750 = cpow(IT_0012 + -IT_0741 + IT_0749, 2);
    const ccomplex_t IT_0751 = cabs(-IT_0749 + 0.25*IT_0001*IT_0750);
    const ccomplex_t IT_0752 = cpow(IT_0751, 0.5);
    const ccomplex_t IT_0753 = ( ((m_Hp
) >= (m_se_R + m_snu_mu
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0748*IT_0752
) : (0
));
    const ccomplex_t IT_0754 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *IT_0007*(m_W*IT_0003 + (-0.5)*IT_0002*IT_0004*IT_0208);
    const ccomplex_t IT_0755 = -IT_0754;
    const ccomplex_t IT_0756 = cpow(IT_0755, 2);
    const ccomplex_t IT_0757 = pow(m_smu_L, 2);
    const ccomplex_t IT_0758 = cpow(IT_0012 + -IT_0741 + IT_0757, 2);
    const ccomplex_t IT_0759 = cabs(-IT_0757 + 0.25*IT_0001*IT_0758);
    const ccomplex_t IT_0760 = cpow(IT_0759, 0.5);
    const ccomplex_t IT_0761 = ( ((m_Hp
) >= (m_smu_L + m_snu_mu
)) ? ((-0.125
      )*IT_0000*IT_0001*IT_0756*IT_0760
) : (0
));
    const ccomplex_t IT_0762 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0002*IT_0007;
    const ccomplex_t IT_0763 = 0.5*IT_0762;
    const ccomplex_t IT_0764 = cpow(IT_0763, 2);
    const ccomplex_t IT_0765 = pow(m_smu_R, 2);
    const ccomplex_t IT_0766 = cpow(IT_0012 + -IT_0741 + IT_0765, 2);
    const ccomplex_t IT_0767 = cabs(-IT_0765 + 0.25*IT_0001*IT_0766);
    const ccomplex_t IT_0768 = cpow(IT_0767, 0.5);
    const ccomplex_t IT_0769 = ( ((m_Hp
) >= (m_smu_R + m_snu_mu
)) ? ((-0.125
      )*IT_0000*IT_0001*IT_0764*IT_0768
) : (0
));
    const ccomplex_t IT_0770 = (0 + _Complex_I*1)*(A_tau*IT_0005*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0007*(m_W*IT_0003*IT_0005*U_stau_00 + (-0.5)
      *m_tau*IT_0002*(m_tau*IT_0004*IT_0005*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0771 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0772 = pow(m_stau_1, 2);
    const ccomplex_t IT_0773 = cpow(IT_0012 + IT_0771 + -IT_0772, 2);
    const ccomplex_t IT_0774 = cabs(-IT_0771 + 0.25*IT_0001*IT_0773);
    const ccomplex_t IT_0775 = cpow(IT_0774, 0.5);
    const ccomplex_t IT_0776 = ( ((m_Hp
) >= (m_stau_1 + m_snu_tau
)) ? (0.125
      *IT_0000*IT_0001*IT_0770*conj(IT_0770)*IT_0775
) : (0
));
    const ccomplex_t IT_0777 = (0 + _Complex_I*1)*(A_tau*IT_0005*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0007*(m_W*IT_0003*IT_0005*U_stau_01 + (-0.5)
      *m_tau*IT_0002*(m_tau*IT_0004*IT_0005*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0778 = pow(m_stau_2, 2);
    const ccomplex_t IT_0779 = cpow(IT_0012 + IT_0771 + -IT_0778, 2);
    const ccomplex_t IT_0780 = cabs(-IT_0771 + 0.25*IT_0001*IT_0779);
    const ccomplex_t IT_0781 = cpow(IT_0780, 0.5);
    const ccomplex_t IT_0782 = ( ((m_Hp
) >= (m_stau_2 + m_snu_tau
)) ? (0.125
      *IT_0000*IT_0001*IT_0777*conj(IT_0777)*IT_0781
) : (0
));
    const ccomplex_t IT_0783 = IT_0017 + IT_0031 + IT_0047 + IT_0061 + IT_0073
       + IT_0085 + IT_0096 + IT_0109 + IT_0119 + IT_0133 + IT_0144 + IT_0155 +
       IT_0166 + IT_0180 + IT_0192 + IT_0204 + IT_0213 + IT_0253 + IT_0288 +
       IT_0323 + IT_0357 + IT_0392 + IT_0426 + IT_0461 + IT_0495 + IT_0504 +
       IT_0512 + IT_0519 + IT_0526 + IT_0532 + IT_0539 + IT_0545 + IT_0552 +
       IT_0558 + IT_0565 + IT_0571 + IT_0578 + IT_0584 + IT_0591 + IT_0597 +
       IT_0603 + IT_0608 + IT_0614 + IT_0619 + IT_0626 + IT_0632 + IT_0638 +
       IT_0643 + IT_0649 + IT_0654 + IT_0661 + IT_0667 + IT_0674 + IT_0681 +
       IT_0688 + IT_0695 + IT_0702 + IT_0708 + IT_0715 + IT_0722 + IT_0729 +
       IT_0736 + IT_0745 + IT_0753 + IT_0761 + IT_0769 + IT_0776 + IT_0782;
    return create_ccomplex_return(IT_0783);
}

