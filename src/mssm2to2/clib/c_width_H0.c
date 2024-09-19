#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_H0(
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
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t m_snu_tau = param->m_snu_tau;
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
    const ccomplex_t IT_0001 = pow(m_H0, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = IT_0005 + IT_0006;
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = cpow(IT_0009, 3);
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0011, 2);
    const ccomplex_t IT_0014 = cpow(IT_0009, 2);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0010 + -IT_0011*(IT_0008*IT_0009*IT_0011 + -IT_0012*(IT_0013 + 
      -IT_0014)));
    const ccomplex_t IT_0016 = (-1.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = pow(m_H0, 2);
    const ccomplex_t IT_0019 = cabs((-0.75)*IT_0018);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_H0
) >= (2*m_H0
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0017*IT_0020
) : (0
));
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_W*e_em*IT_0006*(IT_0008
      *IT_0009 + IT_0011*IT_0012);
    const ccomplex_t IT_0023 = cpow(IT_0022, 2);
    const ccomplex_t IT_0024 = pow(m_W, -4);
    const ccomplex_t IT_0025 = pow(m_W, 2);
    const ccomplex_t IT_0026 = 0.5*IT_0018 + -IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0026, 2);
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = cabs(0.25*IT_0018 + -IT_0025);
    const ccomplex_t IT_0030 = cpow(IT_0029, 0.5);
    const ccomplex_t IT_0031 = ( ((m_H0
) >= (2*m_W
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0023*(1 + 0.5*IT_0028)*IT_0030
) : (0
));
    const ccomplex_t IT_0032 = cpow(IT_0002, -4);
    const ccomplex_t IT_0033 = IT_0004*IT_0011*IT_0012*IT_0032;
    const ccomplex_t IT_0034 = IT_0004*IT_0008*IT_0009*IT_0032;
    const ccomplex_t IT_0035 = cpow(IT_0004, -3);
    const ccomplex_t IT_0036 = IT_0011*IT_0012*IT_0035;
    const ccomplex_t IT_0037 = IT_0008*IT_0009*IT_0035;
    const ccomplex_t IT_0038 = IT_0003*IT_0006*IT_0011*IT_0012;
    const ccomplex_t IT_0039 = IT_0003*IT_0006*IT_0008*IT_0009;
    const ccomplex_t IT_0040 = cpow(IT_0004, -2);
    const ccomplex_t IT_0041 = cpow(IT_0003 + IT_0040, -1);
    const ccomplex_t IT_0042 = (IT_0033 + IT_0034 + IT_0036 + IT_0037 + 2
      *IT_0038 + 2*IT_0039)*IT_0041;
    const ccomplex_t IT_0043 = m_W*e_em;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = cpow(IT_0045, 2);
    const ccomplex_t IT_0047 = pow(m_Z, -4);
    const ccomplex_t IT_0048 = pow(m_Z, 2);
    const ccomplex_t IT_0049 = 0.5*IT_0018 + -IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0049, 2);
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = cabs(0.25*IT_0018 + -IT_0048);
    const ccomplex_t IT_0053 = cpow(IT_0052, 0.5);
    const ccomplex_t IT_0054 = ( ((m_H0
) >= (2*m_Z
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0046*(1 + 0.5*IT_0051)*IT_0053
) : (0
));
    const ccomplex_t IT_0055 = pow(m_W, -1);
    const ccomplex_t IT_0056 = cpow(IT_0012, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_b*e_em*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0058, 2);
    const ccomplex_t IT_0060 = pow(m_b, 2);
    const ccomplex_t IT_0061 = 0.5*IT_0018 + -IT_0060;
    const ccomplex_t IT_0062 = cabs(0.25*IT_0018 + -IT_0060);
    const ccomplex_t IT_0063 = cpow(IT_0062, 0.5);
    const ccomplex_t IT_0064 = ( ((m_H0
) >= (2*m_b
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0059*((-12)*IT_0060 + 12*IT_0061)*IT_0063
) : (0
));
    const ccomplex_t IT_0065 = cpow(IT_0008, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_c*e_em*IT_0006*IT_0009
      *IT_0055*IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = cpow(IT_0067, 2);
    const ccomplex_t IT_0069 = pow(m_c, 2);
    const ccomplex_t IT_0070 = 0.5*IT_0018 + -IT_0069;
    const ccomplex_t IT_0071 = cabs(0.25*IT_0018 + -IT_0069);
    const ccomplex_t IT_0072 = cpow(IT_0071, 0.5);
    const ccomplex_t IT_0073 = ( ((m_H0
) >= (2*m_c
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0068*((-12)*IT_0069 + 12*IT_0070)*IT_0072
) : (0
));
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*m_d*e_em*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = cpow(IT_0075, 2);
    const ccomplex_t IT_0077 = pow(m_d, 2);
    const ccomplex_t IT_0078 = 0.5*IT_0018 + -IT_0077;
    const ccomplex_t IT_0079 = cabs(0.25*IT_0018 + -IT_0077);
    const ccomplex_t IT_0080 = cpow(IT_0079, 0.5);
    const ccomplex_t IT_0081 = ( ((m_H0
) >= (2*m_d
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0076*((-12)*IT_0077 + 12*IT_0078)*IT_0080
) : (0
));
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_e*e_em*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0083, 2);
    const ccomplex_t IT_0085 = pow(m_e, 2);
    const ccomplex_t IT_0086 = 0.5*IT_0018 + -IT_0085;
    const ccomplex_t IT_0087 = cabs(0.25*IT_0018 + -IT_0085);
    const ccomplex_t IT_0088 = cpow(IT_0087, 0.5);
    const ccomplex_t IT_0089 = ( ((m_H0
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *IT_0084*(IT_0085 + -IT_0086)*IT_0088
) : (0
));
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0010 + (-5)*IT_0011*(IT_0008*IT_0009*IT_0011 + (-0.2)*IT_0012*(IT_0013
       + (-5)*IT_0014)));
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0091, 2);
    const ccomplex_t IT_0093 = pow(m_h, 2);
    const ccomplex_t IT_0094 = cabs(0.25*IT_0018 + -IT_0093);
    const ccomplex_t IT_0095 = cpow(IT_0094, 0.5);
    const ccomplex_t IT_0096 = ( ((m_H0
) >= (2*m_h
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0092*IT_0095
) : (0
));
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_s*e_em*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = cpow(IT_0098, 2);
    const ccomplex_t IT_0100 = pow(m_s, 2);
    const ccomplex_t IT_0101 = 0.5*IT_0018 + -IT_0100;
    const ccomplex_t IT_0102 = cabs(0.25*IT_0018 + -IT_0100);
    const ccomplex_t IT_0103 = cpow(IT_0102, 0.5);
    const ccomplex_t IT_0104 = ( ((m_H0
) >= (2*m_s
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0099*((-12)*IT_0100 + 12*IT_0101)*IT_0103
) : (0
));
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0009
      *IT_0055*IT_0065;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0106, 2);
    const ccomplex_t IT_0108 = pow(m_t, 2);
    const ccomplex_t IT_0109 = 0.5*IT_0018 + -IT_0108;
    const ccomplex_t IT_0110 = cabs(0.25*IT_0018 + -IT_0108);
    const ccomplex_t IT_0111 = cpow(IT_0110, 0.5);
    const ccomplex_t IT_0112 = ( ((m_H0
) >= (2*m_t
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0107*((-12)*IT_0108 + 12*IT_0109)*IT_0111
) : (0
));
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*m_u*e_em*IT_0006*IT_0009
      *IT_0055*IT_0065;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = cpow(IT_0114, 2);
    const ccomplex_t IT_0116 = pow(m_u, 2);
    const ccomplex_t IT_0117 = 0.5*IT_0018 + -IT_0116;
    const ccomplex_t IT_0118 = cabs(0.25*IT_0018 + -IT_0116);
    const ccomplex_t IT_0119 = cpow(IT_0118, 0.5);
    const ccomplex_t IT_0120 = ( ((m_H0
) >= (2*m_u
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0115*((-12)*IT_0116 + 12*IT_0117)*IT_0119
) : (0
));
    const ccomplex_t IT_0121 = cpow(IT_0008, 3);
    const ccomplex_t IT_0122 = -IT_0012;
    const ccomplex_t IT_0123 = cpow(IT_0008, 2);
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0009
      *IT_0121 + -IT_0122*((IT_0008*IT_0009 + -IT_0011*IT_0012)*IT_0122 + 
      -IT_0011*IT_0123));
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = cpow(IT_0125, 2);
    const ccomplex_t IT_0127 = pow(m_A0, 2);
    const ccomplex_t IT_0128 = cabs(0.25*IT_0018 + -IT_0127);
    const ccomplex_t IT_0129 = cpow(IT_0128, 0.5);
    const ccomplex_t IT_0130 = ( ((m_H0
) >= (2*m_A0
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0126*IT_0129
) : (0
));
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0010
      *IT_0012 + 5*IT_0011*((-0.2)*IT_0011*(IT_0008*IT_0011 + 5*IT_0009*IT_0012)
       + IT_0008*IT_0014));
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = cpow(IT_0132, 2);
    const ccomplex_t IT_0134 = pow(m_h, 4);
    const ccomplex_t IT_0135 = cabs(-IT_0093 + 0.25*IT_0001*IT_0134);
    const ccomplex_t IT_0136 = cpow(IT_0135, 0.5);
    const ccomplex_t IT_0137 = ( ((m_H0
) >= (m_h + m_H0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0133*IT_0136
) : (0
));
    const ccomplex_t IT_0138 = cpow(IT_0002, -1);
    const ccomplex_t IT_0139 = e_em*IT_0006*(IT_0008*IT_0011 + -IT_0009
      *IT_0012)*IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0139, 2);
    const ccomplex_t IT_0141 = pow(m_Z, -2);
    const ccomplex_t IT_0142 = 0.5*IT_0018 + 0.5*IT_0048 + (-0.5)*IT_0127;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = cpow(IT_0018 + IT_0048 + -IT_0127, 2);
    const ccomplex_t IT_0146 = cabs(-IT_0048 + 0.25*IT_0001*IT_0145);
    const ccomplex_t IT_0147 = cpow(IT_0146, 0.5);
    const ccomplex_t IT_0148 = ( ((m_H0
) >= (m_Z + m_A0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0140*(IT_0018 + -IT_0144)*IT_0147
) : (0
));
    const ccomplex_t IT_0149 = -IT_0006;
    const ccomplex_t IT_0150 = IT_0005 + IT_0149;
    const ccomplex_t IT_0151 = IT_0005 + 3*IT_0006;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_W*e_em*(IT_0009*IT_0121
      *IT_0150 + IT_0122*(IT_0011*IT_0123*IT_0151 + IT_0122*(IT_0011*IT_0012
      *IT_0150 + -IT_0008*IT_0009*IT_0151)));
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = cpow(IT_0153, 2);
    const ccomplex_t IT_0155 = pow(m_Hp, 2);
    const ccomplex_t IT_0156 = cabs(0.25*IT_0018 + -IT_0155);
    const ccomplex_t IT_0157 = cpow(IT_0156, 0.5);
    const ccomplex_t IT_0158 = ( ((m_H0
) >= (2*m_Hp
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0154*IT_0157
) : (0
));
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*e_em*IT_0006*(IT_0008
      *IT_0011 + -IT_0009*IT_0012);
    const ccomplex_t IT_0160 = cpow(IT_0159, 2);
    const ccomplex_t IT_0161 = pow(m_W, -2);
    const ccomplex_t IT_0162 = 0.5*IT_0018 + 0.5*IT_0025 + (-0.5)*IT_0155;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = cpow(IT_0018 + IT_0025 + -IT_0155, 2);
    const ccomplex_t IT_0166 = cabs(-IT_0025 + 0.25*IT_0001*IT_0165);
    const ccomplex_t IT_0167 = cpow(IT_0166, 0.5);
    const ccomplex_t IT_0168 = ( ((m_H0
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0160*(IT_0018 + -IT_0164)*IT_0167
) : (0
));
    const ccomplex_t IT_0169 = -IT_0159;
    const ccomplex_t IT_0170 = cpow(IT_0169, 2);
    const ccomplex_t IT_0171 = ( ((m_H0
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*(IT_0018 + -IT_0164)*IT_0167*IT_0170
) : (0
));
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*e_em*m_mu*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = cpow(IT_0173, 2);
    const ccomplex_t IT_0175 = pow(m_mu, 2);
    const ccomplex_t IT_0176 = 0.5*IT_0018 + -IT_0175;
    const ccomplex_t IT_0177 = cabs(0.25*IT_0018 + -IT_0175);
    const ccomplex_t IT_0178 = cpow(IT_0177, 0.5);
    const ccomplex_t IT_0179 = ( ((m_H0
) >= (2*m_mu
)) ? (0.5*IT_0000*IT_0001
      *IT_0174*(IT_0175 + -IT_0176)*IT_0178
) : (0
));
    const ccomplex_t IT_0180 = U_d1*V_Wp1;
    const ccomplex_t IT_0181 = IT_0011*IT_0180;
    const ccomplex_t IT_0182 = V_u1*U_Wm1;
    const ccomplex_t IT_0183 = IT_0009*IT_0182;
    const ccomplex_t IT_0184 = IT_0006*(IT_0181 + IT_0183);
    const ccomplex_t IT_0185 = 1.4142135623731*e_em*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0190 = IT_0011*IT_0189;
    const ccomplex_t IT_0191 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0192 = IT_0009*IT_0191;
    const ccomplex_t IT_0193 = IT_0006*(IT_0190 + IT_0192);
    const ccomplex_t IT_0194 = 1.4142135623731*e_em*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = pow(m_C_1, 2);
    const ccomplex_t IT_0199 = 0.5*IT_0018 + -IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = cabs(0.25*IT_0018 + -IT_0198);
    const ccomplex_t IT_0202 = cpow(IT_0201, 0.5);
    const ccomplex_t IT_0203 = ( ((m_H0
) >= (2*m_C_1
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0188)*(IT_0197*IT_0198 + (-0.5)*IT_0188*IT_0200) + conj
      (IT_0197)*(IT_0188*IT_0198 + (-0.5)*IT_0197*IT_0200))*IT_0202
) : (0
));
    const ccomplex_t IT_0204 = U_d2*V_Wp2;
    const ccomplex_t IT_0205 = IT_0011*IT_0204;
    const ccomplex_t IT_0206 = V_u2*U_Wm2;
    const ccomplex_t IT_0207 = IT_0009*IT_0206;
    const ccomplex_t IT_0208 = IT_0006*(IT_0205 + IT_0207);
    const ccomplex_t IT_0209 = 1.4142135623731*e_em*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0214 = IT_0011*IT_0213;
    const ccomplex_t IT_0215 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0216 = IT_0009*IT_0215;
    const ccomplex_t IT_0217 = IT_0006*(IT_0214 + IT_0216);
    const ccomplex_t IT_0218 = 1.4142135623731*e_em*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = pow(m_C_2, 2);
    const ccomplex_t IT_0223 = 0.5*IT_0018 + -IT_0222;
    const ccomplex_t IT_0224 = 2*IT_0223;
    const ccomplex_t IT_0225 = cabs(0.25*IT_0018 + -IT_0222);
    const ccomplex_t IT_0226 = cpow(IT_0225, 0.5);
    const ccomplex_t IT_0227 = ( ((m_H0
) >= (2*m_C_2
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0212)*(IT_0221*IT_0222 + (-0.5)*IT_0212*IT_0224) + conj
      (IT_0221)*(IT_0212*IT_0222 + (-0.5)*IT_0221*IT_0224))*IT_0226
) : (0
));
    const ccomplex_t IT_0228 = U_d2*V_Wp1;
    const ccomplex_t IT_0229 = IT_0011*IT_0228;
    const ccomplex_t IT_0230 = V_u1*U_Wm2;
    const ccomplex_t IT_0231 = IT_0009*IT_0230;
    const ccomplex_t IT_0232 = IT_0006*(IT_0229 + IT_0231);
    const ccomplex_t IT_0233 = 1.4142135623731*e_em*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0238 = IT_0011*IT_0237;
    const ccomplex_t IT_0239 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0240 = IT_0009*IT_0239;
    const ccomplex_t IT_0241 = IT_0006*(IT_0238 + IT_0240);
    const ccomplex_t IT_0242 = 1.4142135623731*e_em*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = m_C_1*m_C_2;
    const ccomplex_t IT_0247 = cpow(IT_0018 + IT_0198 + -IT_0222, 2);
    const ccomplex_t IT_0248 = -IT_0198 + 0.25*IT_0001*(IT_0018 + -IT_0198 +
       IT_0222)*(IT_0018 + IT_0198 + -IT_0222) + 0.25*IT_0001*IT_0247;
    const ccomplex_t IT_0249 = 2*IT_0248;
    const ccomplex_t IT_0250 = cabs(-IT_0198 + 0.25*IT_0001*IT_0247);
    const ccomplex_t IT_0251 = cpow(IT_0250, 0.5);
    const ccomplex_t IT_0252 = ( ((m_H0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0236)*(IT_0245*IT_0246 + (-0.5)*IT_0236*IT_0249)
       + conj(IT_0245)*(IT_0236*IT_0246 + (-0.5)*IT_0245*IT_0249))*IT_0251
) : 
      (0
));
    const ccomplex_t IT_0253 = U_d1*V_Wp2;
    const ccomplex_t IT_0254 = IT_0011*IT_0253;
    const ccomplex_t IT_0255 = V_u2*U_Wm1;
    const ccomplex_t IT_0256 = IT_0009*IT_0255;
    const ccomplex_t IT_0257 = IT_0006*(IT_0254 + IT_0256);
    const ccomplex_t IT_0258 = 1.4142135623731*e_em*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0259;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0263 = IT_0011*IT_0262;
    const ccomplex_t IT_0264 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0265 = IT_0009*IT_0264;
    const ccomplex_t IT_0266 = IT_0006*(IT_0263 + IT_0265);
    const ccomplex_t IT_0267 = 1.4142135623731*e_em*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = ( ((m_H0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*IT_0251*(conj(IT_0261)*((-0.5)*IT_0249*IT_0261 + IT_0246
      *IT_0270) + (IT_0246*IT_0261 + (-0.5)*IT_0249*IT_0270)*conj(IT_0270))
) : 
      (0
));
    const ccomplex_t IT_0272 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0273 = IT_0011*IT_0138;
    const ccomplex_t IT_0274 = IT_0272*IT_0273;
    const ccomplex_t IT_0275 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0276 = IT_0009*IT_0138;
    const ccomplex_t IT_0277 = IT_0275*IT_0276;
    const ccomplex_t IT_0278 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0279 = IT_0006*IT_0011;
    const ccomplex_t IT_0280 = IT_0278*IT_0279;
    const ccomplex_t IT_0281 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0282 = IT_0006*IT_0009;
    const ccomplex_t IT_0283 = IT_0281*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*(IT_0274 + -IT_0277 + 
      -IT_0280 + IT_0283);
    const ccomplex_t IT_0285 = 0.5*IT_0284;
    const ccomplex_t IT_0286 = (-0.5)*IT_0284;
    const ccomplex_t IT_0287 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0288 = IT_0273*IT_0287;
    const ccomplex_t IT_0289 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0290 = IT_0276*IT_0289;
    const ccomplex_t IT_0291 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0292 = IT_0279*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0294 = IT_0282*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*(IT_0288 + -IT_0290 + 
      -IT_0292 + IT_0294);
    const ccomplex_t IT_0296 = 0.5*IT_0295;
    const ccomplex_t IT_0297 = (-0.5)*IT_0295;
    const ccomplex_t IT_0298 = pow(m_N_1, 2);
    const ccomplex_t IT_0299 = 0.5*IT_0018 + -IT_0298;
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = ((conj(IT_0285) + -conj(IT_0286))*(IT_0296 + 
      -IT_0297) + (IT_0285 + -IT_0286)*(conj(IT_0296) + -conj(IT_0297)))*IT_0298
       + -(IT_0285*conj(IT_0285) + IT_0286*conj(IT_0286) + IT_0296*conj(IT_0296)
       + IT_0297*conj(IT_0297))*IT_0299 + -(conj(IT_0285)*IT_0286 + IT_0285*conj
      (IT_0286) + conj(IT_0296)*IT_0297 + IT_0296*conj(IT_0297))*IT_0300;
    const ccomplex_t IT_0302 = cabs(0.25*IT_0018 + -IT_0298);
    const ccomplex_t IT_0303 = cpow(IT_0302, 0.5);
    const ccomplex_t IT_0304 = ( ((m_H0
) >= (2*m_N_1
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0301*IT_0303
) : (0
));
    const ccomplex_t IT_0305 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0306 = IT_0273*IT_0305;
    const ccomplex_t IT_0307 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0308 = IT_0276*IT_0307;
    const ccomplex_t IT_0309 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0310 = IT_0279*IT_0309;
    const ccomplex_t IT_0311 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0312 = IT_0282*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0306 + -IT_0308 + 
      -IT_0310 + IT_0312);
    const ccomplex_t IT_0314 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0315 = IT_0273*IT_0314;
    const ccomplex_t IT_0316 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0317 = IT_0276*IT_0316;
    const ccomplex_t IT_0318 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0319 = IT_0279*IT_0318;
    const ccomplex_t IT_0320 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0321 = IT_0282*IT_0320;
    const ccomplex_t IT_0322 = (0 + _Complex_I*1)*(IT_0315 + -IT_0317 + 
      -IT_0319 + IT_0321);
    const ccomplex_t IT_0323 = pow(m_N_2, 2);
    const ccomplex_t IT_0324 = 0.5*IT_0018 + -IT_0323;
    const ccomplex_t IT_0325 = cabs(0.25*IT_0018 + -IT_0323);
    const ccomplex_t IT_0326 = cpow(IT_0325, 0.5);
    const ccomplex_t IT_0327 = ( ((m_H0
) >= (2*m_N_2
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0313)*IT_0322 + IT_0313*conj(IT_0322))*IT_0323 + -
      (IT_0313*conj(IT_0313) + IT_0322*conj(IT_0322))*IT_0324)*IT_0326
) : (0
));
    const ccomplex_t IT_0328 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0329 = IT_0273*IT_0328;
    const ccomplex_t IT_0330 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0331 = IT_0282*IT_0330;
    const ccomplex_t IT_0332 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0333 = IT_0273*IT_0332;
    const ccomplex_t IT_0334 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0335 = IT_0282*IT_0334;
    const ccomplex_t IT_0336 = IT_0329 + IT_0331 + IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0338 = IT_0276*IT_0337;
    const ccomplex_t IT_0339 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0340 = IT_0279*IT_0339;
    const ccomplex_t IT_0341 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0342 = IT_0276*IT_0341;
    const ccomplex_t IT_0343 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0344 = IT_0279*IT_0343;
    const ccomplex_t IT_0345 = -IT_0338 + -IT_0340 + -IT_0342 + -IT_0344;
    const ccomplex_t IT_0346 = IT_0336 + IT_0345;
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0350 = IT_0273*IT_0349;
    const ccomplex_t IT_0351 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0352 = IT_0282*IT_0351;
    const ccomplex_t IT_0353 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0354 = IT_0273*IT_0353;
    const ccomplex_t IT_0355 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0356 = IT_0282*IT_0355;
    const ccomplex_t IT_0357 = IT_0350 + IT_0352 + IT_0354 + IT_0356;
    const ccomplex_t IT_0358 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0359 = IT_0276*IT_0358;
    const ccomplex_t IT_0360 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0361 = IT_0279*IT_0360;
    const ccomplex_t IT_0362 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0363 = IT_0276*IT_0362;
    const ccomplex_t IT_0364 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0365 = IT_0279*IT_0364;
    const ccomplex_t IT_0366 = -IT_0359 + -IT_0361 + -IT_0363 + -IT_0365;
    const ccomplex_t IT_0367 = IT_0357 + IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = 0.5*IT_0368;
    const ccomplex_t IT_0370 = m_N_1*m_N_2;
    const ccomplex_t IT_0371 = cpow(IT_0018 + IT_0298 + -IT_0323, 2);
    const ccomplex_t IT_0372 = -IT_0298 + 0.25*IT_0001*(IT_0018 + -IT_0298 +
       IT_0323)*(IT_0018 + IT_0298 + -IT_0323) + 0.25*IT_0001*IT_0371;
    const ccomplex_t IT_0373 = 2*IT_0372;
    const ccomplex_t IT_0374 = cabs(-IT_0298 + 0.25*IT_0001*IT_0371);
    const ccomplex_t IT_0375 = cpow(IT_0374, 0.5);
    const ccomplex_t IT_0376 = ( ((m_H0
) >= (m_N_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0348)*(IT_0369*IT_0370 + (-0.5)*IT_0348*IT_0373)
       + conj(IT_0369)*(IT_0348*IT_0370 + (-0.5)*IT_0369*IT_0373))*IT_0375
) : 
      (0
));
    const ccomplex_t IT_0377 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0378 = IT_0273*IT_0377;
    const ccomplex_t IT_0379 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0380 = IT_0276*IT_0379;
    const ccomplex_t IT_0381 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0382 = IT_0279*IT_0381;
    const ccomplex_t IT_0383 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0384 = IT_0282*IT_0383;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*(IT_0378 + -IT_0380 + 
      -IT_0382 + IT_0384);
    const ccomplex_t IT_0386 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0387 = IT_0273*IT_0386;
    const ccomplex_t IT_0388 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0389 = IT_0276*IT_0388;
    const ccomplex_t IT_0390 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0391 = IT_0279*IT_0390;
    const ccomplex_t IT_0392 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0393 = IT_0282*IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*(IT_0387 + -IT_0389 + 
      -IT_0391 + IT_0393);
    const ccomplex_t IT_0395 = pow(m_N_3, 2);
    const ccomplex_t IT_0396 = 0.5*IT_0018 + -IT_0395;
    const ccomplex_t IT_0397 = cabs(0.25*IT_0018 + -IT_0395);
    const ccomplex_t IT_0398 = cpow(IT_0397, 0.5);
    const ccomplex_t IT_0399 = ( ((m_H0
) >= (2*m_N_3
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0385)*IT_0394 + IT_0385*conj(IT_0394))*IT_0395 + -
      (IT_0385*conj(IT_0385) + IT_0394*conj(IT_0394))*IT_0396)*IT_0398
) : (0
));
    const ccomplex_t IT_0400 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0401 = IT_0273*IT_0400;
    const ccomplex_t IT_0402 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0403 = IT_0282*IT_0402;
    const ccomplex_t IT_0404 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0405 = IT_0273*IT_0404;
    const ccomplex_t IT_0406 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0407 = IT_0282*IT_0406;
    const ccomplex_t IT_0408 = IT_0401 + IT_0403 + IT_0405 + IT_0407;
    const ccomplex_t IT_0409 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0410 = IT_0276*IT_0409;
    const ccomplex_t IT_0411 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0412 = IT_0279*IT_0411;
    const ccomplex_t IT_0413 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0414 = IT_0276*IT_0413;
    const ccomplex_t IT_0415 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0416 = IT_0279*IT_0415;
    const ccomplex_t IT_0417 = -IT_0410 + -IT_0412 + -IT_0414 + -IT_0416;
    const ccomplex_t IT_0418 = IT_0408 + IT_0417;
    const ccomplex_t IT_0419 = (0 + _Complex_I*1)*IT_0418;
    const ccomplex_t IT_0420 = 0.5*IT_0419;
    const ccomplex_t IT_0421 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0422 = IT_0273*IT_0421;
    const ccomplex_t IT_0423 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0424 = IT_0282*IT_0423;
    const ccomplex_t IT_0425 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0426 = IT_0273*IT_0425;
    const ccomplex_t IT_0427 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0428 = IT_0282*IT_0427;
    const ccomplex_t IT_0429 = IT_0422 + IT_0424 + IT_0426 + IT_0428;
    const ccomplex_t IT_0430 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0431 = IT_0276*IT_0430;
    const ccomplex_t IT_0432 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0433 = IT_0279*IT_0432;
    const ccomplex_t IT_0434 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0435 = IT_0276*IT_0434;
    const ccomplex_t IT_0436 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0437 = IT_0279*IT_0436;
    const ccomplex_t IT_0438 = -IT_0431 + -IT_0433 + -IT_0435 + -IT_0437;
    const ccomplex_t IT_0439 = IT_0429 + IT_0438;
    const ccomplex_t IT_0440 = (0 + _Complex_I*1)*IT_0439;
    const ccomplex_t IT_0441 = 0.5*IT_0440;
    const ccomplex_t IT_0442 = m_N_1*m_N_3;
    const ccomplex_t IT_0443 = cpow(IT_0018 + IT_0298 + -IT_0395, 2);
    const ccomplex_t IT_0444 = -IT_0298 + 0.25*IT_0001*(IT_0018 + -IT_0298 +
       IT_0395)*(IT_0018 + IT_0298 + -IT_0395) + 0.25*IT_0001*IT_0443;
    const ccomplex_t IT_0445 = 2*IT_0444;
    const ccomplex_t IT_0446 = cabs(-IT_0298 + 0.25*IT_0001*IT_0443);
    const ccomplex_t IT_0447 = cpow(IT_0446, 0.5);
    const ccomplex_t IT_0448 = ( ((m_H0
) >= (m_N_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0420)*(IT_0441*IT_0442 + (-0.5)*IT_0420*IT_0445)
       + conj(IT_0441)*(IT_0420*IT_0442 + (-0.5)*IT_0441*IT_0445))*IT_0447
) : 
      (0
));
    const ccomplex_t IT_0449 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0450 = IT_0273*IT_0449;
    const ccomplex_t IT_0451 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0452 = IT_0282*IT_0451;
    const ccomplex_t IT_0453 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0454 = IT_0273*IT_0453;
    const ccomplex_t IT_0455 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0456 = IT_0282*IT_0455;
    const ccomplex_t IT_0457 = IT_0450 + IT_0452 + IT_0454 + IT_0456;
    const ccomplex_t IT_0458 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0459 = IT_0276*IT_0458;
    const ccomplex_t IT_0460 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0461 = IT_0279*IT_0460;
    const ccomplex_t IT_0462 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0463 = IT_0276*IT_0462;
    const ccomplex_t IT_0464 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0465 = IT_0279*IT_0464;
    const ccomplex_t IT_0466 = -IT_0459 + -IT_0461 + -IT_0463 + -IT_0465;
    const ccomplex_t IT_0467 = IT_0457 + IT_0466;
    const ccomplex_t IT_0468 = (0 + _Complex_I*1)*IT_0467;
    const ccomplex_t IT_0469 = 0.5*IT_0468;
    const ccomplex_t IT_0470 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0471 = IT_0273*IT_0470;
    const ccomplex_t IT_0472 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0473 = IT_0282*IT_0472;
    const ccomplex_t IT_0474 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0475 = IT_0273*IT_0474;
    const ccomplex_t IT_0476 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0477 = IT_0282*IT_0476;
    const ccomplex_t IT_0478 = IT_0471 + IT_0473 + IT_0475 + IT_0477;
    const ccomplex_t IT_0479 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0480 = IT_0276*IT_0479;
    const ccomplex_t IT_0481 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0482 = IT_0279*IT_0481;
    const ccomplex_t IT_0483 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0484 = IT_0276*IT_0483;
    const ccomplex_t IT_0485 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0486 = IT_0279*IT_0485;
    const ccomplex_t IT_0487 = -IT_0480 + -IT_0482 + -IT_0484 + -IT_0486;
    const ccomplex_t IT_0488 = IT_0478 + IT_0487;
    const ccomplex_t IT_0489 = (0 + _Complex_I*1)*IT_0488;
    const ccomplex_t IT_0490 = 0.5*IT_0489;
    const ccomplex_t IT_0491 = m_N_2*m_N_3;
    const ccomplex_t IT_0492 = cpow(IT_0018 + IT_0323 + -IT_0395, 2);
    const ccomplex_t IT_0493 = -IT_0323 + 0.25*IT_0001*(IT_0018 + -IT_0323 +
       IT_0395)*(IT_0018 + IT_0323 + -IT_0395) + 0.25*IT_0001*IT_0492;
    const ccomplex_t IT_0494 = 2*IT_0493;
    const ccomplex_t IT_0495 = cabs(-IT_0323 + 0.25*IT_0001*IT_0492);
    const ccomplex_t IT_0496 = cpow(IT_0495, 0.5);
    const ccomplex_t IT_0497 = ( ((m_H0
) >= (m_N_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0469)*(IT_0490*IT_0491 + (-0.5)*IT_0469*IT_0494)
       + conj(IT_0490)*(IT_0469*IT_0491 + (-0.5)*IT_0490*IT_0494))*IT_0496
) : 
      (0
));
    const ccomplex_t IT_0498 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0499 = IT_0273*IT_0498;
    const ccomplex_t IT_0500 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0501 = IT_0276*IT_0500;
    const ccomplex_t IT_0502 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0503 = IT_0279*IT_0502;
    const ccomplex_t IT_0504 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0505 = IT_0282*IT_0504;
    const ccomplex_t IT_0506 = (0 + _Complex_I*1)*(IT_0499 + -IT_0501 + 
      -IT_0503 + IT_0505);
    const ccomplex_t IT_0507 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0508 = IT_0273*IT_0507;
    const ccomplex_t IT_0509 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0510 = IT_0276*IT_0509;
    const ccomplex_t IT_0511 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0512 = IT_0279*IT_0511;
    const ccomplex_t IT_0513 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0514 = IT_0282*IT_0513;
    const ccomplex_t IT_0515 = (0 + _Complex_I*1)*(IT_0508 + -IT_0510 + 
      -IT_0512 + IT_0514);
    const ccomplex_t IT_0516 = pow(m_N_4, 2);
    const ccomplex_t IT_0517 = 0.5*IT_0018 + -IT_0516;
    const ccomplex_t IT_0518 = cabs(0.25*IT_0018 + -IT_0516);
    const ccomplex_t IT_0519 = cpow(IT_0518, 0.5);
    const ccomplex_t IT_0520 = ( ((m_H0
) >= (2*m_N_4
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0506)*IT_0515 + IT_0506*conj(IT_0515))*IT_0516 + -
      (IT_0506*conj(IT_0506) + IT_0515*conj(IT_0515))*IT_0517)*IT_0519
) : (0
));
    const ccomplex_t IT_0521 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0522 = IT_0273*IT_0521;
    const ccomplex_t IT_0523 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0524 = IT_0282*IT_0523;
    const ccomplex_t IT_0525 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0526 = IT_0273*IT_0525;
    const ccomplex_t IT_0527 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0528 = IT_0282*IT_0527;
    const ccomplex_t IT_0529 = IT_0522 + IT_0524 + IT_0526 + IT_0528;
    const ccomplex_t IT_0530 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0531 = IT_0276*IT_0530;
    const ccomplex_t IT_0532 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0533 = IT_0279*IT_0532;
    const ccomplex_t IT_0534 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0535 = IT_0276*IT_0534;
    const ccomplex_t IT_0536 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0537 = IT_0279*IT_0536;
    const ccomplex_t IT_0538 = -IT_0531 + -IT_0533 + -IT_0535 + -IT_0537;
    const ccomplex_t IT_0539 = IT_0529 + IT_0538;
    const ccomplex_t IT_0540 = (0 + _Complex_I*1)*IT_0539;
    const ccomplex_t IT_0541 = 0.5*IT_0540;
    const ccomplex_t IT_0542 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0543 = IT_0273*IT_0542;
    const ccomplex_t IT_0544 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0545 = IT_0282*IT_0544;
    const ccomplex_t IT_0546 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0547 = IT_0273*IT_0546;
    const ccomplex_t IT_0548 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0549 = IT_0282*IT_0548;
    const ccomplex_t IT_0550 = IT_0543 + IT_0545 + IT_0547 + IT_0549;
    const ccomplex_t IT_0551 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0552 = IT_0276*IT_0551;
    const ccomplex_t IT_0553 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0554 = IT_0279*IT_0553;
    const ccomplex_t IT_0555 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0556 = IT_0276*IT_0555;
    const ccomplex_t IT_0557 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0558 = IT_0279*IT_0557;
    const ccomplex_t IT_0559 = -IT_0552 + -IT_0554 + -IT_0556 + -IT_0558;
    const ccomplex_t IT_0560 = IT_0550 + IT_0559;
    const ccomplex_t IT_0561 = (0 + _Complex_I*1)*IT_0560;
    const ccomplex_t IT_0562 = 0.5*IT_0561;
    const ccomplex_t IT_0563 = m_N_1*m_N_4;
    const ccomplex_t IT_0564 = cpow(IT_0018 + IT_0298 + -IT_0516, 2);
    const ccomplex_t IT_0565 = -IT_0298 + 0.25*IT_0001*(IT_0018 + -IT_0298 +
       IT_0516)*(IT_0018 + IT_0298 + -IT_0516) + 0.25*IT_0001*IT_0564;
    const ccomplex_t IT_0566 = 2*IT_0565;
    const ccomplex_t IT_0567 = cabs(-IT_0298 + 0.25*IT_0001*IT_0564);
    const ccomplex_t IT_0568 = cpow(IT_0567, 0.5);
    const ccomplex_t IT_0569 = ( ((m_H0
) >= (m_N_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0541)*(IT_0562*IT_0563 + (-0.5)*IT_0541*IT_0566)
       + conj(IT_0562)*(IT_0541*IT_0563 + (-0.5)*IT_0562*IT_0566))*IT_0568
) : 
      (0
));
    const ccomplex_t IT_0570 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0571 = IT_0273*IT_0570;
    const ccomplex_t IT_0572 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0573 = IT_0282*IT_0572;
    const ccomplex_t IT_0574 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0575 = IT_0273*IT_0574;
    const ccomplex_t IT_0576 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0577 = IT_0282*IT_0576;
    const ccomplex_t IT_0578 = IT_0571 + IT_0573 + IT_0575 + IT_0577;
    const ccomplex_t IT_0579 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0580 = IT_0276*IT_0579;
    const ccomplex_t IT_0581 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0582 = IT_0279*IT_0581;
    const ccomplex_t IT_0583 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0584 = IT_0276*IT_0583;
    const ccomplex_t IT_0585 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0586 = IT_0279*IT_0585;
    const ccomplex_t IT_0587 = -IT_0580 + -IT_0582 + -IT_0584 + -IT_0586;
    const ccomplex_t IT_0588 = IT_0578 + IT_0587;
    const ccomplex_t IT_0589 = (0 + _Complex_I*1)*IT_0588;
    const ccomplex_t IT_0590 = 0.5*IT_0589;
    const ccomplex_t IT_0591 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0592 = IT_0273*IT_0591;
    const ccomplex_t IT_0593 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0594 = IT_0282*IT_0593;
    const ccomplex_t IT_0595 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0596 = IT_0273*IT_0595;
    const ccomplex_t IT_0597 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0598 = IT_0282*IT_0597;
    const ccomplex_t IT_0599 = IT_0592 + IT_0594 + IT_0596 + IT_0598;
    const ccomplex_t IT_0600 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0601 = IT_0276*IT_0600;
    const ccomplex_t IT_0602 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0603 = IT_0279*IT_0602;
    const ccomplex_t IT_0604 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0605 = IT_0276*IT_0604;
    const ccomplex_t IT_0606 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0607 = IT_0279*IT_0606;
    const ccomplex_t IT_0608 = -IT_0601 + -IT_0603 + -IT_0605 + -IT_0607;
    const ccomplex_t IT_0609 = IT_0599 + IT_0608;
    const ccomplex_t IT_0610 = (0 + _Complex_I*1)*IT_0609;
    const ccomplex_t IT_0611 = 0.5*IT_0610;
    const ccomplex_t IT_0612 = m_N_2*m_N_4;
    const ccomplex_t IT_0613 = cpow(IT_0018 + IT_0323 + -IT_0516, 2);
    const ccomplex_t IT_0614 = -IT_0323 + 0.25*IT_0001*(IT_0018 + -IT_0323 +
       IT_0516)*(IT_0018 + IT_0323 + -IT_0516) + 0.25*IT_0001*IT_0613;
    const ccomplex_t IT_0615 = 2*IT_0614;
    const ccomplex_t IT_0616 = cabs(-IT_0323 + 0.25*IT_0001*IT_0613);
    const ccomplex_t IT_0617 = cpow(IT_0616, 0.5);
    const ccomplex_t IT_0618 = ( ((m_H0
) >= (m_N_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0590)*(IT_0611*IT_0612 + (-0.5)*IT_0590*IT_0615)
       + conj(IT_0611)*(IT_0590*IT_0612 + (-0.5)*IT_0611*IT_0615))*IT_0617
) : 
      (0
));
    const ccomplex_t IT_0619 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0620 = IT_0273*IT_0619;
    const ccomplex_t IT_0621 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0622 = IT_0282*IT_0621;
    const ccomplex_t IT_0623 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0624 = IT_0273*IT_0623;
    const ccomplex_t IT_0625 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0626 = IT_0282*IT_0625;
    const ccomplex_t IT_0627 = IT_0620 + IT_0622 + IT_0624 + IT_0626;
    const ccomplex_t IT_0628 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0629 = IT_0276*IT_0628;
    const ccomplex_t IT_0630 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0631 = IT_0279*IT_0630;
    const ccomplex_t IT_0632 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0633 = IT_0276*IT_0632;
    const ccomplex_t IT_0634 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0635 = IT_0279*IT_0634;
    const ccomplex_t IT_0636 = -IT_0629 + -IT_0631 + -IT_0633 + -IT_0635;
    const ccomplex_t IT_0637 = IT_0627 + IT_0636;
    const ccomplex_t IT_0638 = (0 + _Complex_I*1)*IT_0637;
    const ccomplex_t IT_0639 = 0.5*IT_0638;
    const ccomplex_t IT_0640 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0641 = IT_0273*IT_0640;
    const ccomplex_t IT_0642 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0643 = IT_0282*IT_0642;
    const ccomplex_t IT_0644 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0645 = IT_0273*IT_0644;
    const ccomplex_t IT_0646 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0647 = IT_0282*IT_0646;
    const ccomplex_t IT_0648 = IT_0641 + IT_0643 + IT_0645 + IT_0647;
    const ccomplex_t IT_0649 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0650 = IT_0276*IT_0649;
    const ccomplex_t IT_0651 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0652 = IT_0279*IT_0651;
    const ccomplex_t IT_0653 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0654 = IT_0276*IT_0653;
    const ccomplex_t IT_0655 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0656 = IT_0279*IT_0655;
    const ccomplex_t IT_0657 = -IT_0650 + -IT_0652 + -IT_0654 + -IT_0656;
    const ccomplex_t IT_0658 = IT_0648 + IT_0657;
    const ccomplex_t IT_0659 = (0 + _Complex_I*1)*IT_0658;
    const ccomplex_t IT_0660 = 0.5*IT_0659;
    const ccomplex_t IT_0661 = m_N_3*m_N_4;
    const ccomplex_t IT_0662 = cpow(IT_0018 + IT_0395 + -IT_0516, 2);
    const ccomplex_t IT_0663 = -IT_0395 + 0.25*IT_0001*(IT_0018 + -IT_0395 +
       IT_0516)*(IT_0018 + IT_0395 + -IT_0516) + 0.25*IT_0001*IT_0662;
    const ccomplex_t IT_0664 = 2*IT_0663;
    const ccomplex_t IT_0665 = cabs(-IT_0395 + 0.25*IT_0001*IT_0662);
    const ccomplex_t IT_0666 = cpow(IT_0665, 0.5);
    const ccomplex_t IT_0667 = ( ((m_H0
) >= (m_N_3 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0639)*(IT_0660*IT_0661 + (-0.5)*IT_0639*IT_0664)
       + conj(IT_0660)*(IT_0639*IT_0661 + (-0.5)*IT_0660*IT_0664))*IT_0666
) : 
      (0
));
    const ccomplex_t IT_0668 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0011
      *IT_0055*IT_0056;
    const ccomplex_t IT_0669 = (-0.5)*IT_0668;
    const ccomplex_t IT_0670 = cpow(IT_0669, 2);
    const ccomplex_t IT_0671 = pow(m_tau, 2);
    const ccomplex_t IT_0672 = 0.5*IT_0018 + -IT_0671;
    const ccomplex_t IT_0673 = cabs(0.25*IT_0018 + -IT_0671);
    const ccomplex_t IT_0674 = cpow(IT_0673, 0.5);
    const ccomplex_t IT_0675 = ( ((m_H0
) >= (2*m_tau
)) ? (0.5*IT_0000
      *IT_0001*IT_0670*(IT_0671 + -IT_0672)*IT_0674
) : (0
));
    const ccomplex_t IT_0676 = cpow(V_tb, 2);
    const ccomplex_t IT_0677 = m_b*mu_h*IT_0055;
    const ccomplex_t IT_0678 = -IT_0008*IT_0009;
    const ccomplex_t IT_0679 = IT_0011*IT_0012;
    const ccomplex_t IT_0680 = IT_0678 + IT_0679;
    const ccomplex_t IT_0681 = (-2)*U_sb_00*((IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0060 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0680)
      *conj(U_sb_00) + -1./2*IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_10)) +
       U_sb_10*(IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_00) + (-2)*(IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0060 + -1./3*m_W*IT_0003*IT_0004*IT_0680)*conj
      (U_sb_10));
    const ccomplex_t IT_0682 = (0 + _Complex_I*-1)*(e_em*IT_0681 + (
      -1.4142135623731)*A_b*IT_0011*IT_0676*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0683 = (-0.5)*IT_0682;
    const ccomplex_t IT_0684 = pow(m_sb_1, 2);
    const ccomplex_t IT_0685 = cabs(0.25*IT_0018 + -IT_0684);
    const ccomplex_t IT_0686 = cpow(IT_0685, 0.5);
    const ccomplex_t IT_0687 = ( ((m_H0
) >= (2*m_sb_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0683*conj(IT_0683)*IT_0686
) : (0
));
    const ccomplex_t IT_0688 = (-2)*U_sb_01*((IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0060 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0680)
      *conj(U_sb_01) + -1./2*IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_11)) +
       U_sb_11*(IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_01) + (-2)*(IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0060 + -1./3*m_W*IT_0003*IT_0004*IT_0680)*conj
      (U_sb_11));
    const ccomplex_t IT_0689 = (0 + _Complex_I*-1)*(e_em*IT_0688 + (
      -1.4142135623731)*A_b*IT_0011*IT_0676*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0690 = (-0.5)*IT_0689;
    const ccomplex_t IT_0691 = pow(m_sb_2, 2);
    const ccomplex_t IT_0692 = cabs(0.25*IT_0018 + -IT_0691);
    const ccomplex_t IT_0693 = cpow(IT_0692, 0.5);
    const ccomplex_t IT_0694 = ( ((m_H0
) >= (2*m_sb_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0690*conj(IT_0690)*IT_0693
) : (0
));
    const ccomplex_t IT_0695 = (-2)*U_sb_00*((IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0060 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0680)
      *conj(U_sb_01) + -1./2*IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_11)) +
       U_sb_10*(IT_0006*IT_0009*IT_0056*IT_0677*conj(U_sb_01) + (-2)*(IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0060 + -1./3*m_W*IT_0003*IT_0004*IT_0680)*conj
      (U_sb_11));
    const ccomplex_t IT_0696 = (0 + _Complex_I*-1)*(e_em*IT_0695 + (
      -1.4142135623731)*A_b*IT_0011*IT_0676*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0697 = (-0.5)*IT_0696;
    const ccomplex_t IT_0698 = cpow(IT_0018 + IT_0684 + -IT_0691, 2);
    const ccomplex_t IT_0699 = cabs(-IT_0684 + 0.25*IT_0001*IT_0698);
    const ccomplex_t IT_0700 = cpow(IT_0699, 0.5);
    const ccomplex_t IT_0701 = ( ((m_H0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0697*conj(IT_0697)*IT_0700
) : (0
));
    const ccomplex_t IT_0702 = (-2)*conj(U_sb_00)*((IT_0006*IT_0011*IT_0055
      *IT_0056*IT_0060 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)
      *IT_0680)*U_sb_01 + -1./2*IT_0006*IT_0009*IT_0056*IT_0677*U_sb_11) + conj
      (U_sb_10)*(IT_0006*IT_0009*IT_0056*IT_0677*U_sb_01 + (-2)*(IT_0006*IT_0011
      *IT_0055*IT_0056*IT_0060 + -1./3*m_W*IT_0003*IT_0004*IT_0680)*U_sb_11);
    const ccomplex_t IT_0703 = (0 + _Complex_I*-1)*(e_em*IT_0702 + (
      -1.4142135623731)*A_b*IT_0011*IT_0676*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0704 = (-0.5)*IT_0703;
    const ccomplex_t IT_0705 = ( ((m_H0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0700*IT_0704*conj(IT_0704)
) : (0
));
    const ccomplex_t IT_0706 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0055*IT_0065*IT_0069);
    const ccomplex_t IT_0707 = -IT_0706;
    const ccomplex_t IT_0708 = cpow(IT_0707, 2);
    const ccomplex_t IT_0709 = pow(m_sc_L, 2);
    const ccomplex_t IT_0710 = cabs(0.25*IT_0018 + -IT_0709);
    const ccomplex_t IT_0711 = cpow(IT_0710, 0.5);
    const ccomplex_t IT_0712 = ( ((m_H0
) >= (2*m_sc_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0708*IT_0711
) : (0
));
    const ccomplex_t IT_0713 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0009*IT_0055
      *IT_0065*IT_0069);
    const ccomplex_t IT_0714 = (-0.666666666666667)*IT_0713;
    const ccomplex_t IT_0715 = cpow(IT_0714, 2);
    const ccomplex_t IT_0716 = pow(m_sc_R, 2);
    const ccomplex_t IT_0717 = cabs(0.25*IT_0018 + -IT_0716);
    const ccomplex_t IT_0718 = cpow(IT_0717, 0.5);
    const ccomplex_t IT_0719 = ( ((m_H0
) >= (2*m_sc_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0715*IT_0718
) : (0
));
    const ccomplex_t IT_0720 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0006
      *IT_0011*IT_0055*IT_0065;
    const ccomplex_t IT_0721 = 0.5*IT_0720;
    const ccomplex_t IT_0722 = cpow(IT_0721, 2);
    const ccomplex_t IT_0723 = cpow(IT_0018 + IT_0709 + -IT_0716, 2);
    const ccomplex_t IT_0724 = cabs(-IT_0709 + 0.25*IT_0001*IT_0723);
    const ccomplex_t IT_0725 = cpow(IT_0724, 0.5);
    const ccomplex_t IT_0726 = ( ((m_H0
) >= (m_sc_L + m_sc_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0722*IT_0725
) : (0
));
    const ccomplex_t IT_0727 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0055*IT_0056*IT_0077);
    const ccomplex_t IT_0728 = -IT_0727;
    const ccomplex_t IT_0729 = cpow(IT_0728, 2);
    const ccomplex_t IT_0730 = pow(m_sd_L, 2);
    const ccomplex_t IT_0731 = cabs(0.25*IT_0018 + -IT_0730);
    const ccomplex_t IT_0732 = cpow(IT_0731, 0.5);
    const ccomplex_t IT_0733 = ( ((m_H0
) >= (2*m_sd_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0729*IT_0732
) : (0
));
    const ccomplex_t IT_0734 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_10;
    const ccomplex_t IT_0735 = (-0.5)*IT_0734;
    const ccomplex_t IT_0736 = cpow(IT_0018 + -IT_0684 + IT_0730, 2);
    const ccomplex_t IT_0737 = cabs(-IT_0730 + 0.25*IT_0001*IT_0736);
    const ccomplex_t IT_0738 = cpow(IT_0737, 0.5);
    const ccomplex_t IT_0739 = ( ((m_H0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0735*conj(IT_0735)*IT_0738
) : (0
));
    const ccomplex_t IT_0740 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_10);
    const ccomplex_t IT_0741 = (-0.5)*IT_0740;
    const ccomplex_t IT_0742 = ( ((m_H0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0738*IT_0741*conj(IT_0741)
) : (0
));
    const ccomplex_t IT_0743 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0744 = (-0.5)*IT_0743;
    const ccomplex_t IT_0745 = cpow(IT_0018 + -IT_0691 + IT_0730, 2);
    const ccomplex_t IT_0746 = cabs(-IT_0730 + 0.25*IT_0001*IT_0745);
    const ccomplex_t IT_0747 = cpow(IT_0746, 0.5);
    const ccomplex_t IT_0748 = ( ((m_H0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0744*conj(IT_0744)*IT_0747
) : (0
));
    const ccomplex_t IT_0749 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_11;
    const ccomplex_t IT_0750 = (-0.5)*IT_0749;
    const ccomplex_t IT_0751 = ( ((m_H0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0747*IT_0750*conj(IT_0750)
) : (0
));
    const ccomplex_t IT_0752 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + 3*IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0077);
    const ccomplex_t IT_0753 = 0.333333333333333*IT_0752;
    const ccomplex_t IT_0754 = cpow(IT_0753, 2);
    const ccomplex_t IT_0755 = pow(m_sd_R, 2);
    const ccomplex_t IT_0756 = cabs(0.25*IT_0018 + -IT_0755);
    const ccomplex_t IT_0757 = cpow(IT_0756, 0.5);
    const ccomplex_t IT_0758 = ( ((m_H0
) >= (2*m_sd_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0754*IT_0757
) : (0
));
    const ccomplex_t IT_0759 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_00;
    const ccomplex_t IT_0760 = (-0.5)*IT_0759;
    const ccomplex_t IT_0761 = cpow(IT_0018 + -IT_0684 + IT_0755, 2);
    const ccomplex_t IT_0762 = cabs(-IT_0755 + 0.25*IT_0001*IT_0761);
    const ccomplex_t IT_0763 = cpow(IT_0762, 0.5);
    const ccomplex_t IT_0764 = ( ((m_H0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0760*conj(IT_0760)*IT_0763
) : (0
));
    const ccomplex_t IT_0765 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_00);
    const ccomplex_t IT_0766 = (-0.5)*IT_0765;
    const ccomplex_t IT_0767 = ( ((m_H0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0763*IT_0766*conj(IT_0766)
) : (0
));
    const ccomplex_t IT_0768 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0769 = (-0.5)*IT_0768;
    const ccomplex_t IT_0770 = cpow(IT_0018 + -IT_0691 + IT_0755, 2);
    const ccomplex_t IT_0771 = cabs(-IT_0755 + 0.25*IT_0001*IT_0770);
    const ccomplex_t IT_0772 = cpow(IT_0771, 0.5);
    const ccomplex_t IT_0773 = ( ((m_H0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0769*conj(IT_0769)*IT_0772
) : (0
));
    const ccomplex_t IT_0774 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_01;
    const ccomplex_t IT_0775 = (-0.5)*IT_0774;
    const ccomplex_t IT_0776 = ( ((m_H0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0772*IT_0775*conj(IT_0775)
) : (0
));
    const ccomplex_t IT_0777 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0011 + -m_d*e_em*mu_h*IT_0006*IT_0009*IT_0055*IT_0056);
    const ccomplex_t IT_0778 = (-0.5)*IT_0777;
    const ccomplex_t IT_0779 = cpow(IT_0778, 2);
    const ccomplex_t IT_0780 = cpow(IT_0018 + -IT_0730 + IT_0755, 2);
    const ccomplex_t IT_0781 = cabs(-IT_0755 + 0.25*IT_0001*IT_0780);
    const ccomplex_t IT_0782 = cpow(IT_0781, 0.5);
    const ccomplex_t IT_0783 = ( ((m_H0
) >= (m_sd_L + m_sd_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0779*IT_0782
) : (0
));
    const ccomplex_t IT_0784 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0085);
    const ccomplex_t IT_0785 = -IT_0784;
    const ccomplex_t IT_0786 = cpow(IT_0785, 2);
    const ccomplex_t IT_0787 = pow(m_se_L, 2);
    const ccomplex_t IT_0788 = cabs(0.25*IT_0018 + -IT_0787);
    const ccomplex_t IT_0789 = cpow(IT_0788, 0.5);
    const ccomplex_t IT_0790 = ( ((m_H0
) >= (2*m_se_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0786*IT_0789
) : (0
));
    const ccomplex_t IT_0791 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0085);
    const ccomplex_t IT_0792 = cpow(IT_0791, 2);
    const ccomplex_t IT_0793 = pow(m_se_R, 2);
    const ccomplex_t IT_0794 = cabs(0.25*IT_0018 + -IT_0793);
    const ccomplex_t IT_0795 = cpow(IT_0794, 0.5);
    const ccomplex_t IT_0796 = ( ((m_H0
) >= (2*m_se_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0792*IT_0795
) : (0
));
    const ccomplex_t IT_0797 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0006
      *IT_0009*IT_0055*IT_0056;
    const ccomplex_t IT_0798 = 0.5*IT_0797;
    const ccomplex_t IT_0799 = cpow(IT_0798, 2);
    const ccomplex_t IT_0800 = cpow(IT_0018 + -IT_0787 + IT_0793, 2);
    const ccomplex_t IT_0801 = cabs(-IT_0793 + 0.25*IT_0001*IT_0800);
    const ccomplex_t IT_0802 = cpow(IT_0801, 0.5);
    const ccomplex_t IT_0803 = ( ((m_H0
) >= (m_se_L + m_se_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0799*IT_0802
) : (0
));
    const ccomplex_t IT_0804 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0055*IT_0056*IT_0100);
    const ccomplex_t IT_0805 = -IT_0804;
    const ccomplex_t IT_0806 = cpow(IT_0805, 2);
    const ccomplex_t IT_0807 = pow(m_ss_L, 2);
    const ccomplex_t IT_0808 = cabs(0.25*IT_0018 + -IT_0807);
    const ccomplex_t IT_0809 = cpow(IT_0808, 0.5);
    const ccomplex_t IT_0810 = ( ((m_H0
) >= (2*m_ss_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0806*IT_0809
) : (0
));
    const ccomplex_t IT_0811 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_10);
    const ccomplex_t IT_0812 = (-0.5)*IT_0811;
    const ccomplex_t IT_0813 = cpow(IT_0018 + -IT_0684 + IT_0807, 2);
    const ccomplex_t IT_0814 = cabs(-IT_0807 + 0.25*IT_0001*IT_0813);
    const ccomplex_t IT_0815 = cpow(IT_0814, 0.5);
    const ccomplex_t IT_0816 = ( ((m_H0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0812*conj(IT_0812)*IT_0815
) : (0
));
    const ccomplex_t IT_0817 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_10;
    const ccomplex_t IT_0818 = (-0.5)*IT_0817;
    const ccomplex_t IT_0819 = ( ((m_H0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0815*IT_0818*conj(IT_0818)
) : (0
));
    const ccomplex_t IT_0820 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0821 = (-0.5)*IT_0820;
    const ccomplex_t IT_0822 = cpow(IT_0018 + -IT_0691 + IT_0807, 2);
    const ccomplex_t IT_0823 = cabs(-IT_0807 + 0.25*IT_0001*IT_0822);
    const ccomplex_t IT_0824 = cpow(IT_0823, 0.5);
    const ccomplex_t IT_0825 = ( ((m_H0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0821*conj(IT_0821)*IT_0824
) : (0
));
    const ccomplex_t IT_0826 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_11;
    const ccomplex_t IT_0827 = (-0.5)*IT_0826;
    const ccomplex_t IT_0828 = ( ((m_H0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0824*IT_0827*conj(IT_0827)
) : (0
));
    const ccomplex_t IT_0829 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0011;
    const ccomplex_t IT_0830 = (-0.5)*IT_0829;
    const ccomplex_t IT_0831 = cpow(IT_0018 + IT_0755 + -IT_0807, 2);
    const ccomplex_t IT_0832 = cabs(-IT_0755 + 0.25*IT_0001*IT_0831);
    const ccomplex_t IT_0833 = cpow(IT_0832, 0.5);
    const ccomplex_t IT_0834 = ( ((m_H0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0830*conj(IT_0830)*IT_0833
) : (0
));
    const ccomplex_t IT_0835 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0011;
    const ccomplex_t IT_0836 = (-0.5)*IT_0835;
    const ccomplex_t IT_0837 = ( ((m_H0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0833*IT_0836*conj(IT_0836)
) : (0
));
    const ccomplex_t IT_0838 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + 3*IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0100);
    const ccomplex_t IT_0839 = 0.333333333333333*IT_0838;
    const ccomplex_t IT_0840 = cpow(IT_0839, 2);
    const ccomplex_t IT_0841 = pow(m_ss_R, 2);
    const ccomplex_t IT_0842 = cabs(0.25*IT_0018 + -IT_0841);
    const ccomplex_t IT_0843 = cpow(IT_0842, 0.5);
    const ccomplex_t IT_0844 = ( ((m_H0
) >= (2*m_ss_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0840*IT_0843
) : (0
));
    const ccomplex_t IT_0845 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_00;
    const ccomplex_t IT_0846 = (-0.5)*IT_0845;
    const ccomplex_t IT_0847 = cpow(IT_0018 + -IT_0684 + IT_0841, 2);
    const ccomplex_t IT_0848 = cabs(-IT_0841 + 0.25*IT_0001*IT_0847);
    const ccomplex_t IT_0849 = cpow(IT_0848, 0.5);
    const ccomplex_t IT_0850 = ( ((m_H0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0846*conj(IT_0846)*IT_0849
) : (0
));
    const ccomplex_t IT_0851 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_00);
    const ccomplex_t IT_0852 = (-0.5)*IT_0851;
    const ccomplex_t IT_0853 = ( ((m_H0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0849*IT_0852*conj(IT_0852)
) : (0
));
    const ccomplex_t IT_0854 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0855 = (-0.5)*IT_0854;
    const ccomplex_t IT_0856 = cpow(IT_0018 + -IT_0691 + IT_0841, 2);
    const ccomplex_t IT_0857 = cabs(-IT_0841 + 0.25*IT_0001*IT_0856);
    const ccomplex_t IT_0858 = cpow(IT_0857, 0.5);
    const ccomplex_t IT_0859 = ( ((m_H0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0855*conj(IT_0855)*IT_0858
) : (0
));
    const ccomplex_t IT_0860 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_01;
    const ccomplex_t IT_0861 = (-0.5)*IT_0860;
    const ccomplex_t IT_0862 = ( ((m_H0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0858*IT_0861*conj(IT_0861)
) : (0
));
    const ccomplex_t IT_0863 = cpow(IT_0018 + -IT_0730 + IT_0841, 2);
    const ccomplex_t IT_0864 = cabs(-IT_0841 + 0.25*IT_0001*IT_0863);
    const ccomplex_t IT_0865 = cpow(IT_0864, 0.5);
    const ccomplex_t IT_0866 = ( ((m_H0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0830*conj(IT_0830)*IT_0865
) : (0
));
    const ccomplex_t IT_0867 = ( ((m_H0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0836*conj(IT_0836)*IT_0865
) : (0
));
    const ccomplex_t IT_0868 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0011 + -m_s*e_em*mu_h*IT_0006*IT_0009*IT_0055*IT_0056);
    const ccomplex_t IT_0869 = (-0.5)*IT_0868;
    const ccomplex_t IT_0870 = cpow(IT_0869, 2);
    const ccomplex_t IT_0871 = cpow(IT_0018 + -IT_0807 + IT_0841, 2);
    const ccomplex_t IT_0872 = cabs(-IT_0841 + 0.25*IT_0001*IT_0871);
    const ccomplex_t IT_0873 = cpow(IT_0872, 0.5);
    const ccomplex_t IT_0874 = ( ((m_H0
) >= (m_ss_L + m_ss_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0870*IT_0873
) : (0
));
    const ccomplex_t IT_0875 = m_t*mu_h*IT_0055;
    const ccomplex_t IT_0876 = -IT_0008*IT_0009 + IT_0679;
    const ccomplex_t IT_0877 = (-2)*U_st_00*((IT_0006*IT_0009*IT_0055*IT_0065
      *IT_0108 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0876)*conj(U_st_00) + -1./2*IT_0006*IT_0011*IT_0065*IT_0875*conj
      (U_st_10)) + U_st_10*(IT_0006*IT_0011*IT_0065*IT_0875*conj(U_st_00) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0055*IT_0065*IT_0108 + m_W
      *IT_0003*IT_0004*IT_0876)*conj(U_st_10));
    const ccomplex_t IT_0878 = (0 + _Complex_I*-1)*(e_em*IT_0877 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0879 = (-0.5)*IT_0878;
    const ccomplex_t IT_0880 = pow(m_st_1, 2);
    const ccomplex_t IT_0881 = cabs(0.25*IT_0018 + -IT_0880);
    const ccomplex_t IT_0882 = cpow(IT_0881, 0.5);
    const ccomplex_t IT_0883 = ( ((m_H0
) >= (2*m_st_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0879*conj(IT_0879)*IT_0882
) : (0
));
    const ccomplex_t IT_0884 = (-2)*U_st_01*((IT_0006*IT_0009*IT_0055*IT_0065
      *IT_0108 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0876)*conj(U_st_01) + -1./2*IT_0006*IT_0011*IT_0065*IT_0875*conj
      (U_st_11)) + U_st_11*(IT_0006*IT_0011*IT_0065*IT_0875*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0055*IT_0065*IT_0108 + m_W
      *IT_0003*IT_0004*IT_0876)*conj(U_st_11));
    const ccomplex_t IT_0885 = (0 + _Complex_I*-1)*(e_em*IT_0884 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0886 = (-0.5)*IT_0885;
    const ccomplex_t IT_0887 = pow(m_st_2, 2);
    const ccomplex_t IT_0888 = cabs(0.25*IT_0018 + -IT_0887);
    const ccomplex_t IT_0889 = cpow(IT_0888, 0.5);
    const ccomplex_t IT_0890 = ( ((m_H0
) >= (2*m_st_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0886*conj(IT_0886)*IT_0889
) : (0
));
    const ccomplex_t IT_0891 = (-2)*U_st_00*((IT_0006*IT_0009*IT_0055*IT_0065
      *IT_0108 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0876)*conj(U_st_01) + -1./2*IT_0006*IT_0011*IT_0065*IT_0875*conj
      (U_st_11)) + U_st_10*(IT_0006*IT_0011*IT_0065*IT_0875*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0055*IT_0065*IT_0108 + m_W
      *IT_0003*IT_0004*IT_0876)*conj(U_st_11));
    const ccomplex_t IT_0892 = (0 + _Complex_I*-1)*(e_em*IT_0891 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0893 = (-0.5)*IT_0892;
    const ccomplex_t IT_0894 = cpow(IT_0018 + IT_0880 + -IT_0887, 2);
    const ccomplex_t IT_0895 = cabs(-IT_0880 + 0.25*IT_0001*IT_0894);
    const ccomplex_t IT_0896 = cpow(IT_0895, 0.5);
    const ccomplex_t IT_0897 = ( ((m_H0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0893*conj(IT_0893)*IT_0896
) : (0
));
    const ccomplex_t IT_0898 = (-2)*conj(U_st_00)*((IT_0006*IT_0009*IT_0055
      *IT_0065*IT_0108 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)
      *IT_0006)*IT_0876)*U_st_01 + -1./2*IT_0006*IT_0011*IT_0065*IT_0875*U_st_11
      ) + conj(U_st_10)*(IT_0006*IT_0011*IT_0065*IT_0875*U_st_01 + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0055*IT_0065*IT_0108 + m_W
      *IT_0003*IT_0004*IT_0876)*U_st_11);
    const ccomplex_t IT_0899 = (0 + _Complex_I*-1)*(e_em*IT_0898 + (
      -1.4142135623731)*A_t*IT_0009*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0900 = (-0.5)*IT_0899;
    const ccomplex_t IT_0901 = ( ((m_H0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0896*IT_0900*conj(IT_0900)
) : (0
));
    const ccomplex_t IT_0902 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0055*IT_0065*IT_0116);
    const ccomplex_t IT_0903 = -IT_0902;
    const ccomplex_t IT_0904 = cpow(IT_0903, 2);
    const ccomplex_t IT_0905 = pow(m_su_L, 2);
    const ccomplex_t IT_0906 = cabs(0.25*IT_0018 + -IT_0905);
    const ccomplex_t IT_0907 = cpow(IT_0906, 0.5);
    const ccomplex_t IT_0908 = ( ((m_H0
) >= (2*m_su_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0904*IT_0907
) : (0
));
    const ccomplex_t IT_0909 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0009*IT_0055
      *IT_0065*IT_0116);
    const ccomplex_t IT_0910 = (-0.666666666666667)*IT_0909;
    const ccomplex_t IT_0911 = cpow(IT_0910, 2);
    const ccomplex_t IT_0912 = pow(m_su_R, 2);
    const ccomplex_t IT_0913 = cabs(0.25*IT_0018 + -IT_0912);
    const ccomplex_t IT_0914 = cpow(IT_0913, 0.5);
    const ccomplex_t IT_0915 = ( ((m_H0
) >= (2*m_su_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0911*IT_0914
) : (0
));
    const ccomplex_t IT_0916 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0006
      *IT_0011*IT_0055*IT_0065;
    const ccomplex_t IT_0917 = 0.5*IT_0916;
    const ccomplex_t IT_0918 = cpow(IT_0917, 2);
    const ccomplex_t IT_0919 = cpow(IT_0018 + IT_0905 + -IT_0912, 2);
    const ccomplex_t IT_0920 = cabs(-IT_0905 + 0.25*IT_0001*IT_0919);
    const ccomplex_t IT_0921 = cpow(IT_0920, 0.5);
    const ccomplex_t IT_0922 = ( ((m_H0
) >= (m_su_L + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0918*IT_0921
) : (0
));
    const ccomplex_t IT_0923 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0175);
    const ccomplex_t IT_0924 = -IT_0923;
    const ccomplex_t IT_0925 = cpow(IT_0924, 2);
    const ccomplex_t IT_0926 = pow(m_smu_L, 2);
    const ccomplex_t IT_0927 = cabs(0.25*IT_0018 + -IT_0926);
    const ccomplex_t IT_0928 = cpow(IT_0927, 0.5);
    const ccomplex_t IT_0929 = ( ((m_H0
) >= (2*m_smu_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0925*IT_0928
) : (0
));
    const ccomplex_t IT_0930 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006*IT_0011*IT_0055*IT_0056
      *IT_0175);
    const ccomplex_t IT_0931 = cpow(IT_0930, 2);
    const ccomplex_t IT_0932 = pow(m_smu_R, 2);
    const ccomplex_t IT_0933 = cabs(0.25*IT_0018 + -IT_0932);
    const ccomplex_t IT_0934 = cpow(IT_0933, 0.5);
    const ccomplex_t IT_0935 = ( ((m_H0
) >= (2*m_smu_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0931*IT_0934
) : (0
));
    const ccomplex_t IT_0936 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0006
      *IT_0009*IT_0055*IT_0056;
    const ccomplex_t IT_0937 = 0.5*IT_0936;
    const ccomplex_t IT_0938 = cpow(IT_0937, 2);
    const ccomplex_t IT_0939 = cpow(IT_0018 + -IT_0926 + IT_0932, 2);
    const ccomplex_t IT_0940 = cabs(-IT_0932 + 0.25*IT_0001*IT_0939);
    const ccomplex_t IT_0941 = cpow(IT_0940, 0.5);
    const ccomplex_t IT_0942 = ( ((m_H0
) >= (m_smu_L + m_smu_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0938*IT_0941
) : (0
));
    const ccomplex_t IT_0943 = (0 + _Complex_I*-1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0009 + -IT_0011*IT_0012);
    const ccomplex_t IT_0944 = (-0.5)*IT_0943;
    const ccomplex_t IT_0945 = cpow(IT_0944, 2);
    const ccomplex_t IT_0946 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0947 = cabs(0.25*IT_0018 + -IT_0946);
    const ccomplex_t IT_0948 = cpow(IT_0947, 0.5);
    const ccomplex_t IT_0949 = ( ((m_H0
) >= (2*m_snu_mu
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0945*IT_0948
) : (0
));
    const ccomplex_t IT_0950 = mu_h*m_tau*IT_0055;
    const ccomplex_t IT_0951 = U_stau_10*(IT_0006*IT_0009*IT_0056*IT_0950*conj
      (U_stau_00) + (-2)*(IT_0006*IT_0011*IT_0055*IT_0056*IT_0671 + -m_W*IT_0003
      *IT_0004*IT_0680)*conj(U_stau_10)) + (-2)*U_stau_00*((IT_0006*IT_0011
      *IT_0055*IT_0056*IT_0671 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0680)
      *conj(U_stau_00) + -1./2*IT_0006*IT_0009*IT_0056*IT_0950*conj(U_stau_10));
    const ccomplex_t IT_0952 = (0 + _Complex_I*-1)*(e_em*IT_0951 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0953 = (-0.5)*IT_0952;
    const ccomplex_t IT_0954 = pow(m_stau_1, 2);
    const ccomplex_t IT_0955 = cabs(0.25*IT_0018 + -IT_0954);
    const ccomplex_t IT_0956 = cpow(IT_0955, 0.5);
    const ccomplex_t IT_0957 = ( ((m_H0
) >= (2*m_stau_1
)) ? (0.125*IT_0000
      *IT_0001*IT_0953*conj(IT_0953)*IT_0956
) : (0
));
    const ccomplex_t IT_0958 = U_stau_11*(IT_0006*IT_0009*IT_0056*IT_0950*conj
      (U_stau_01) + (-2)*(IT_0006*IT_0011*IT_0055*IT_0056*IT_0671 + -m_W*IT_0003
      *IT_0004*IT_0680)*conj(U_stau_11)) + (-2)*U_stau_01*((IT_0006*IT_0011
      *IT_0055*IT_0056*IT_0671 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0680)
      *conj(U_stau_01) + -1./2*IT_0006*IT_0009*IT_0056*IT_0950*conj(U_stau_11));
    const ccomplex_t IT_0959 = (0 + _Complex_I*-1)*(e_em*IT_0958 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0960 = (-0.5)*IT_0959;
    const ccomplex_t IT_0961 = pow(m_stau_2, 2);
    const ccomplex_t IT_0962 = cabs(0.25*IT_0018 + -IT_0961);
    const ccomplex_t IT_0963 = cpow(IT_0962, 0.5);
    const ccomplex_t IT_0964 = ( ((m_H0
) >= (2*m_stau_2
)) ? (0.125*IT_0000
      *IT_0001*IT_0960*conj(IT_0960)*IT_0963
) : (0
));
    const ccomplex_t IT_0965 = U_stau_10*(IT_0006*IT_0009*IT_0056*IT_0950*conj
      (U_stau_01) + (-2)*(IT_0006*IT_0011*IT_0055*IT_0056*IT_0671 + -m_W*IT_0003
      *IT_0004*IT_0680)*conj(U_stau_11)) + (-2)*U_stau_00*((IT_0006*IT_0011
      *IT_0055*IT_0056*IT_0671 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0680)
      *conj(U_stau_01) + -1./2*IT_0006*IT_0009*IT_0056*IT_0950*conj(U_stau_11));
    const ccomplex_t IT_0966 = (0 + _Complex_I*-1)*(e_em*IT_0965 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0967 = (-0.5)*IT_0966;
    const ccomplex_t IT_0968 = cpow(IT_0018 + IT_0954 + -IT_0961, 2);
    const ccomplex_t IT_0969 = cabs(-IT_0954 + 0.25*IT_0001*IT_0968);
    const ccomplex_t IT_0970 = cpow(IT_0969, 0.5);
    const ccomplex_t IT_0971 = ( ((m_H0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0967*conj(IT_0967)*IT_0970
) : (0
));
    const ccomplex_t IT_0972 = conj(U_stau_10)*(IT_0006*IT_0009*IT_0056
      *IT_0950*U_stau_01 + (-2)*(IT_0006*IT_0011*IT_0055*IT_0056*IT_0671 + -m_W
      *IT_0003*IT_0004*IT_0680)*U_stau_11) + (-2)*conj(U_stau_00)*((IT_0006
      *IT_0011*IT_0055*IT_0056*IT_0671 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)
      *IT_0680)*U_stau_01 + -1./2*IT_0006*IT_0009*IT_0056*IT_0950*U_stau_11);
    const ccomplex_t IT_0973 = (0 + _Complex_I*-1)*(e_em*IT_0972 + (
      -1.4142135623731)*A_tau*IT_0011*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0974 = (-0.5)*IT_0973;
    const ccomplex_t IT_0975 = ( ((m_H0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0970*IT_0974*conj(IT_0974)
) : (0
));
    const ccomplex_t IT_0976 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0977 = cabs(0.25*IT_0018 + -IT_0976);
    const ccomplex_t IT_0978 = cpow(IT_0977, 0.5);
    const ccomplex_t IT_0979 = ( ((m_H0
) >= (2*m_snu_tau
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0945*IT_0978
) : (0
));
    const ccomplex_t IT_0980 = IT_0021 + IT_0031 + IT_0054 + IT_0064 + IT_0073
       + IT_0081 + IT_0089 + IT_0096 + IT_0104 + IT_0112 + IT_0120 + IT_0130 +
       IT_0137 + IT_0148 + IT_0158 + IT_0168 + IT_0171 + IT_0179 + IT_0203 +
       IT_0227 + IT_0252 + IT_0271 + IT_0304 + IT_0327 + IT_0376 + IT_0399 +
       IT_0448 + IT_0497 + IT_0520 + IT_0569 + IT_0618 + IT_0667 + IT_0675 +
       IT_0687 + IT_0694 + IT_0701 + IT_0705 + IT_0712 + IT_0719 + 2*IT_0726 +
       IT_0733 + IT_0739 + IT_0742 + IT_0748 + IT_0751 + IT_0758 + IT_0764 +
       IT_0767 + IT_0773 + IT_0776 + 2*IT_0783 + IT_0790 + IT_0796 + 2*IT_0803 +
       IT_0810 + IT_0816 + IT_0819 + IT_0825 + IT_0828 + IT_0834 + IT_0837 +
       IT_0844 + IT_0850 + IT_0853 + IT_0859 + IT_0862 + IT_0866 + IT_0867 + 2
      *IT_0874 + IT_0883 + IT_0890 + IT_0897 + IT_0901 + IT_0908 + IT_0915 + 2
      *IT_0922 + IT_0929 + IT_0935 + 2*IT_0942 + 2*IT_0949 + IT_0957 + IT_0964 +
       IT_0971 + IT_0975 + IT_0979;
    return create_ccomplex_return(IT_0980);
}

