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
    const ccomplex_t IT_0022 = pow(m_W, -1);
    const ccomplex_t IT_0023 = cpow(IT_0012, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_b*e_em*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = pow(m_b, 2);
    const ccomplex_t IT_0028 = 0.5*IT_0018 + -IT_0027;
    const ccomplex_t IT_0029 = cabs(0.25*IT_0018 + -IT_0027);
    const ccomplex_t IT_0030 = cpow(IT_0029, 0.5);
    const ccomplex_t IT_0031 = ( ((m_H0
) >= (2*m_b
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0026*((-12)*IT_0027 + 12*IT_0028)*IT_0030
) : (0
));
    const ccomplex_t IT_0032 = cpow(IT_0008, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_c*e_em*IT_0006*IT_0009
      *IT_0022*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = cpow(IT_0034, 2);
    const ccomplex_t IT_0036 = pow(m_c, 2);
    const ccomplex_t IT_0037 = 0.5*IT_0018 + -IT_0036;
    const ccomplex_t IT_0038 = cabs(0.25*IT_0018 + -IT_0036);
    const ccomplex_t IT_0039 = cpow(IT_0038, 0.5);
    const ccomplex_t IT_0040 = ( ((m_H0
) >= (2*m_c
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0035*((-12)*IT_0036 + 12*IT_0037)*IT_0039
) : (0
));
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_d*e_em*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = cpow(IT_0042, 2);
    const ccomplex_t IT_0044 = pow(m_d, 2);
    const ccomplex_t IT_0045 = 0.5*IT_0018 + -IT_0044;
    const ccomplex_t IT_0046 = cabs(0.25*IT_0018 + -IT_0044);
    const ccomplex_t IT_0047 = cpow(IT_0046, 0.5);
    const ccomplex_t IT_0048 = ( ((m_H0
) >= (2*m_d
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0043*((-12)*IT_0044 + 12*IT_0045)*IT_0047
) : (0
));
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_e*e_em*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = pow(m_e, 2);
    const ccomplex_t IT_0053 = 0.5*IT_0018 + -IT_0052;
    const ccomplex_t IT_0054 = cabs(0.25*IT_0018 + -IT_0052);
    const ccomplex_t IT_0055 = cpow(IT_0054, 0.5);
    const ccomplex_t IT_0056 = ( ((m_H0
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *IT_0051*(IT_0052 + -IT_0053)*IT_0055
) : (0
));
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0010 + (-5)*IT_0011*(IT_0008*IT_0009*IT_0011 + (-0.2)*IT_0012*(IT_0013
       + (-5)*IT_0014)));
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0058, 2);
    const ccomplex_t IT_0060 = pow(m_h, 2);
    const ccomplex_t IT_0061 = cabs(0.25*IT_0018 + -IT_0060);
    const ccomplex_t IT_0062 = cpow(IT_0061, 0.5);
    const ccomplex_t IT_0063 = ( ((m_H0
) >= (2*m_h
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0059*IT_0062
) : (0
));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_s*e_em*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0065, 2);
    const ccomplex_t IT_0067 = pow(m_s, 2);
    const ccomplex_t IT_0068 = 0.5*IT_0018 + -IT_0067;
    const ccomplex_t IT_0069 = cabs(0.25*IT_0018 + -IT_0067);
    const ccomplex_t IT_0070 = cpow(IT_0069, 0.5);
    const ccomplex_t IT_0071 = ( ((m_H0
) >= (2*m_s
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0066*((-12)*IT_0067 + 12*IT_0068)*IT_0070
) : (0
));
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0009
      *IT_0022*IT_0032;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = cpow(IT_0073, 2);
    const ccomplex_t IT_0075 = pow(m_t, 2);
    const ccomplex_t IT_0076 = 0.5*IT_0018 + -IT_0075;
    const ccomplex_t IT_0077 = cabs(0.25*IT_0018 + -IT_0075);
    const ccomplex_t IT_0078 = cpow(IT_0077, 0.5);
    const ccomplex_t IT_0079 = ( ((m_H0
) >= (2*m_t
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0074*((-12)*IT_0075 + 12*IT_0076)*IT_0078
) : (0
));
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*m_u*e_em*IT_0006*IT_0009
      *IT_0022*IT_0032;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0081, 2);
    const ccomplex_t IT_0083 = pow(m_u, 2);
    const ccomplex_t IT_0084 = 0.5*IT_0018 + -IT_0083;
    const ccomplex_t IT_0085 = cabs(0.25*IT_0018 + -IT_0083);
    const ccomplex_t IT_0086 = cpow(IT_0085, 0.5);
    const ccomplex_t IT_0087 = ( ((m_H0
) >= (2*m_u
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0082*((-12)*IT_0083 + 12*IT_0084)*IT_0086
) : (0
));
    const ccomplex_t IT_0088 = cpow(IT_0008, 3);
    const ccomplex_t IT_0089 = -IT_0012;
    const ccomplex_t IT_0090 = cpow(IT_0008, 2);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0009
      *IT_0088 + -IT_0089*((IT_0008*IT_0009 + -IT_0011*IT_0012)*IT_0089 + 
      -IT_0011*IT_0090));
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = cpow(IT_0092, 2);
    const ccomplex_t IT_0094 = pow(m_A0, 2);
    const ccomplex_t IT_0095 = cabs(0.25*IT_0018 + -IT_0094);
    const ccomplex_t IT_0096 = cpow(IT_0095, 0.5);
    const ccomplex_t IT_0097 = ( ((m_H0
) >= (2*m_A0
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0093*IT_0096
) : (0
));
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0010
      *IT_0012 + 5*IT_0011*((-0.2)*IT_0011*(IT_0008*IT_0011 + 5*IT_0009*IT_0012)
       + IT_0008*IT_0014));
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = cpow(IT_0099, 2);
    const ccomplex_t IT_0101 = pow(m_h, 4);
    const ccomplex_t IT_0102 = cabs(-IT_0060 + 0.25*IT_0001*IT_0101);
    const ccomplex_t IT_0103 = cpow(IT_0102, 0.5);
    const ccomplex_t IT_0104 = ( ((m_H0
) >= (m_h + m_H0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0100*IT_0103
) : (0
));
    const ccomplex_t IT_0105 = cpow(IT_0002, -1);
    const ccomplex_t IT_0106 = e_em*IT_0006*(IT_0008*IT_0011 + -IT_0009
      *IT_0012)*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0106, 2);
    const ccomplex_t IT_0108 = pow(m_Z, -2);
    const ccomplex_t IT_0109 = pow(m_Z, 2);
    const ccomplex_t IT_0110 = 0.5*IT_0018 + (-0.5)*IT_0094 + 0.5*IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0110, 2);
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = cpow(IT_0018 + -IT_0094 + IT_0109, 2);
    const ccomplex_t IT_0114 = cabs(-IT_0109 + 0.25*IT_0001*IT_0113);
    const ccomplex_t IT_0115 = cpow(IT_0114, 0.5);
    const ccomplex_t IT_0116 = ( ((m_H0
) >= (m_Z + m_A0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0107*(IT_0018 + -IT_0112)*IT_0115
) : (0
));
    const ccomplex_t IT_0117 = -IT_0006;
    const ccomplex_t IT_0118 = IT_0005 + IT_0117;
    const ccomplex_t IT_0119 = IT_0005 + 3*IT_0006;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_W*e_em*(IT_0009*IT_0088
      *IT_0118 + IT_0089*(IT_0011*IT_0090*IT_0119 + IT_0089*(IT_0011*IT_0012
      *IT_0118 + -IT_0008*IT_0009*IT_0119)));
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = cpow(IT_0121, 2);
    const ccomplex_t IT_0123 = pow(m_Hp, 2);
    const ccomplex_t IT_0124 = cabs(0.25*IT_0018 + -IT_0123);
    const ccomplex_t IT_0125 = cpow(IT_0124, 0.5);
    const ccomplex_t IT_0126 = ( ((m_H0
) >= (2*m_Hp
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0122*IT_0125
) : (0
));
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*e_em*IT_0006*(IT_0008
      *IT_0011 + -IT_0009*IT_0012);
    const ccomplex_t IT_0128 = cpow(IT_0127, 2);
    const ccomplex_t IT_0129 = pow(m_W, -2);
    const ccomplex_t IT_0130 = pow(m_W, 2);
    const ccomplex_t IT_0131 = 0.5*IT_0018 + (-0.5)*IT_0123 + 0.5*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0131, 2);
    const ccomplex_t IT_0133 = IT_0129*IT_0132;
    const ccomplex_t IT_0134 = cpow(IT_0018 + -IT_0123 + IT_0130, 2);
    const ccomplex_t IT_0135 = cabs(-IT_0130 + 0.25*IT_0001*IT_0134);
    const ccomplex_t IT_0136 = cpow(IT_0135, 0.5);
    const ccomplex_t IT_0137 = ( ((m_H0
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0128*(IT_0018 + -IT_0133)*IT_0136
) : (0
));
    const ccomplex_t IT_0138 = -IT_0127;
    const ccomplex_t IT_0139 = cpow(IT_0138, 2);
    const ccomplex_t IT_0140 = ( ((m_H0
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*(IT_0018 + -IT_0133)*IT_0136*IT_0139
) : (0
));
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*e_em*m_mu*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = pow(m_mu, 2);
    const ccomplex_t IT_0145 = 0.5*IT_0018 + -IT_0144;
    const ccomplex_t IT_0146 = cabs(0.25*IT_0018 + -IT_0144);
    const ccomplex_t IT_0147 = cpow(IT_0146, 0.5);
    const ccomplex_t IT_0148 = ( ((m_H0
) >= (2*m_mu
)) ? (0.5*IT_0000*IT_0001
      *IT_0143*(IT_0144 + -IT_0145)*IT_0147
) : (0
));
    const ccomplex_t IT_0149 = U_d1*V_Wp1;
    const ccomplex_t IT_0150 = IT_0011*IT_0149;
    const ccomplex_t IT_0151 = V_u1*U_Wm1;
    const ccomplex_t IT_0152 = IT_0009*IT_0151;
    const ccomplex_t IT_0153 = IT_0006*(IT_0150 + IT_0152);
    const ccomplex_t IT_0154 = 1.4142135623731*e_em*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0159 = IT_0011*IT_0158;
    const ccomplex_t IT_0160 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0161 = IT_0009*IT_0160;
    const ccomplex_t IT_0162 = IT_0006*(IT_0159 + IT_0161);
    const ccomplex_t IT_0163 = 1.4142135623731*e_em*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = pow(m_C_1, 2);
    const ccomplex_t IT_0168 = 0.5*IT_0018 + -IT_0167;
    const ccomplex_t IT_0169 = 2*IT_0168;
    const ccomplex_t IT_0170 = cabs(0.25*IT_0018 + -IT_0167);
    const ccomplex_t IT_0171 = cpow(IT_0170, 0.5);
    const ccomplex_t IT_0172 = ( ((m_H0
) >= (2*m_C_1
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0157)*(IT_0166*IT_0167 + (-0.5)*IT_0157*IT_0169) + conj
      (IT_0166)*(IT_0157*IT_0167 + (-0.5)*IT_0166*IT_0169))*IT_0171
) : (0
));
    const ccomplex_t IT_0173 = U_d2*V_Wp2;
    const ccomplex_t IT_0174 = IT_0011*IT_0173;
    const ccomplex_t IT_0175 = V_u2*U_Wm2;
    const ccomplex_t IT_0176 = IT_0009*IT_0175;
    const ccomplex_t IT_0177 = IT_0006*(IT_0174 + IT_0176);
    const ccomplex_t IT_0178 = 1.4142135623731*e_em*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = 0.5*IT_0179;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0183 = IT_0011*IT_0182;
    const ccomplex_t IT_0184 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0185 = IT_0009*IT_0184;
    const ccomplex_t IT_0186 = IT_0006*(IT_0183 + IT_0185);
    const ccomplex_t IT_0187 = 1.4142135623731*e_em*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = pow(m_C_2, 2);
    const ccomplex_t IT_0192 = 0.5*IT_0018 + -IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0192;
    const ccomplex_t IT_0194 = cabs(0.25*IT_0018 + -IT_0191);
    const ccomplex_t IT_0195 = cpow(IT_0194, 0.5);
    const ccomplex_t IT_0196 = ( ((m_H0
) >= (2*m_C_2
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0181)*(IT_0190*IT_0191 + (-0.5)*IT_0181*IT_0193) + conj
      (IT_0190)*(IT_0181*IT_0191 + (-0.5)*IT_0190*IT_0193))*IT_0195
) : (0
));
    const ccomplex_t IT_0197 = U_d2*V_Wp1;
    const ccomplex_t IT_0198 = IT_0011*IT_0197;
    const ccomplex_t IT_0199 = V_u1*U_Wm2;
    const ccomplex_t IT_0200 = IT_0009*IT_0199;
    const ccomplex_t IT_0201 = IT_0006*(IT_0198 + IT_0200);
    const ccomplex_t IT_0202 = 1.4142135623731*e_em*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0207 = IT_0011*IT_0206;
    const ccomplex_t IT_0208 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0209 = IT_0009*IT_0208;
    const ccomplex_t IT_0210 = IT_0006*(IT_0207 + IT_0209);
    const ccomplex_t IT_0211 = 1.4142135623731*e_em*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = m_C_1*m_C_2;
    const ccomplex_t IT_0216 = cpow(IT_0018 + IT_0167 + -IT_0191, 2);
    const ccomplex_t IT_0217 = -IT_0167 + 0.25*IT_0001*(IT_0018 + -IT_0167 +
       IT_0191)*(IT_0018 + IT_0167 + -IT_0191) + 0.25*IT_0001*IT_0216;
    const ccomplex_t IT_0218 = 2*IT_0217;
    const ccomplex_t IT_0219 = cabs(-IT_0167 + 0.25*IT_0001*IT_0216);
    const ccomplex_t IT_0220 = cpow(IT_0219, 0.5);
    const ccomplex_t IT_0221 = ( ((m_H0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0205)*(IT_0214*IT_0215 + (-0.5)*IT_0205*IT_0218)
       + conj(IT_0214)*(IT_0205*IT_0215 + (-0.5)*IT_0214*IT_0218))*IT_0220
) : 
      (0
));
    const ccomplex_t IT_0222 = U_d1*V_Wp2;
    const ccomplex_t IT_0223 = IT_0011*IT_0222;
    const ccomplex_t IT_0224 = V_u2*U_Wm1;
    const ccomplex_t IT_0225 = IT_0009*IT_0224;
    const ccomplex_t IT_0226 = IT_0006*(IT_0223 + IT_0225);
    const ccomplex_t IT_0227 = 1.4142135623731*e_em*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = 0.5*IT_0228;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0232 = IT_0011*IT_0231;
    const ccomplex_t IT_0233 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0234 = IT_0009*IT_0233;
    const ccomplex_t IT_0235 = IT_0006*(IT_0232 + IT_0234);
    const ccomplex_t IT_0236 = 1.4142135623731*e_em*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = ( ((m_H0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*IT_0220*(conj(IT_0230)*((-0.5)*IT_0218*IT_0230 + IT_0215
      *IT_0239) + (IT_0215*IT_0230 + (-0.5)*IT_0218*IT_0239)*conj(IT_0239))
) : 
      (0
));
    const ccomplex_t IT_0241 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0242 = IT_0011*IT_0105;
    const ccomplex_t IT_0243 = IT_0241*IT_0242;
    const ccomplex_t IT_0244 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0245 = IT_0009*IT_0105;
    const ccomplex_t IT_0246 = IT_0244*IT_0245;
    const ccomplex_t IT_0247 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0248 = IT_0006*IT_0011;
    const ccomplex_t IT_0249 = IT_0247*IT_0248;
    const ccomplex_t IT_0250 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0251 = IT_0006*IT_0009;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0243 + -IT_0246 + 
      -IT_0249 + IT_0252);
    const ccomplex_t IT_0254 = (-0.5)*IT_0253;
    const ccomplex_t IT_0255 = 0.5*IT_0253;
    const ccomplex_t IT_0256 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0257 = IT_0242*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0259 = IT_0245*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0261 = IT_0248*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0263 = IT_0251*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0257 + -IT_0259 + 
      -IT_0261 + IT_0263);
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = 0.5*IT_0264;
    const ccomplex_t IT_0267 = pow(m_N_1, 2);
    const ccomplex_t IT_0268 = 0.5*IT_0018 + -IT_0267;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = ((conj(IT_0254) + -conj(IT_0255))*(IT_0265 + 
      -IT_0266) + (IT_0254 + -IT_0255)*(conj(IT_0265) + -conj(IT_0266)))*IT_0267
       + -(IT_0254*conj(IT_0254) + IT_0255*conj(IT_0255) + IT_0265*conj(IT_0265)
       + IT_0266*conj(IT_0266))*IT_0268 + -(conj(IT_0254)*IT_0255 + IT_0254*conj
      (IT_0255) + conj(IT_0265)*IT_0266 + IT_0265*conj(IT_0266))*IT_0269;
    const ccomplex_t IT_0271 = cabs(0.25*IT_0018 + -IT_0267);
    const ccomplex_t IT_0272 = cpow(IT_0271, 0.5);
    const ccomplex_t IT_0273 = ( ((m_H0
) >= (2*m_N_1
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0270*IT_0272
) : (0
));
    const ccomplex_t IT_0274 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0275 = IT_0242*IT_0274;
    const ccomplex_t IT_0276 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0277 = IT_0245*IT_0276;
    const ccomplex_t IT_0278 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0279 = IT_0248*IT_0278;
    const ccomplex_t IT_0280 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0281 = IT_0251*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*(IT_0275 + -IT_0277 + 
      -IT_0279 + IT_0281);
    const ccomplex_t IT_0283 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0284 = IT_0242*IT_0283;
    const ccomplex_t IT_0285 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0286 = IT_0245*IT_0285;
    const ccomplex_t IT_0287 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0288 = IT_0248*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0290 = IT_0251*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*(IT_0284 + -IT_0286 + 
      -IT_0288 + IT_0290);
    const ccomplex_t IT_0292 = pow(m_N_2, 2);
    const ccomplex_t IT_0293 = 0.5*IT_0018 + -IT_0292;
    const ccomplex_t IT_0294 = cabs(0.25*IT_0018 + -IT_0292);
    const ccomplex_t IT_0295 = cpow(IT_0294, 0.5);
    const ccomplex_t IT_0296 = ( ((m_H0
) >= (2*m_N_2
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0282)*IT_0291 + IT_0282*conj(IT_0291))*IT_0292 + -
      (IT_0282*conj(IT_0282) + IT_0291*conj(IT_0291))*IT_0293)*IT_0295
) : (0
));
    const ccomplex_t IT_0297 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0298 = IT_0242*IT_0297;
    const ccomplex_t IT_0299 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0300 = IT_0251*IT_0299;
    const ccomplex_t IT_0301 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0302 = IT_0242*IT_0301;
    const ccomplex_t IT_0303 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0304 = IT_0251*IT_0303;
    const ccomplex_t IT_0305 = IT_0298 + IT_0300 + IT_0302 + IT_0304;
    const ccomplex_t IT_0306 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0307 = IT_0245*IT_0306;
    const ccomplex_t IT_0308 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0309 = IT_0248*IT_0308;
    const ccomplex_t IT_0310 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0311 = IT_0245*IT_0310;
    const ccomplex_t IT_0312 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0313 = IT_0248*IT_0312;
    const ccomplex_t IT_0314 = -IT_0307 + -IT_0309 + -IT_0311 + -IT_0313;
    const ccomplex_t IT_0315 = IT_0305 + IT_0314;
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*IT_0315;
    const ccomplex_t IT_0317 = 0.5*IT_0316;
    const ccomplex_t IT_0318 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0319 = IT_0242*IT_0318;
    const ccomplex_t IT_0320 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0321 = IT_0251*IT_0320;
    const ccomplex_t IT_0322 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0323 = IT_0242*IT_0322;
    const ccomplex_t IT_0324 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0325 = IT_0251*IT_0324;
    const ccomplex_t IT_0326 = IT_0319 + IT_0321 + IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0328 = IT_0245*IT_0327;
    const ccomplex_t IT_0329 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0330 = IT_0248*IT_0329;
    const ccomplex_t IT_0331 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0332 = IT_0245*IT_0331;
    const ccomplex_t IT_0333 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0334 = IT_0248*IT_0333;
    const ccomplex_t IT_0335 = -IT_0328 + -IT_0330 + -IT_0332 + -IT_0334;
    const ccomplex_t IT_0336 = IT_0326 + IT_0335;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*IT_0336;
    const ccomplex_t IT_0338 = 0.5*IT_0337;
    const ccomplex_t IT_0339 = m_N_1*m_N_2;
    const ccomplex_t IT_0340 = cpow(IT_0018 + IT_0267 + -IT_0292, 2);
    const ccomplex_t IT_0341 = -IT_0267 + 0.25*IT_0001*(IT_0018 + -IT_0267 +
       IT_0292)*(IT_0018 + IT_0267 + -IT_0292) + 0.25*IT_0001*IT_0340;
    const ccomplex_t IT_0342 = 2*IT_0341;
    const ccomplex_t IT_0343 = cabs(-IT_0267 + 0.25*IT_0001*IT_0340);
    const ccomplex_t IT_0344 = cpow(IT_0343, 0.5);
    const ccomplex_t IT_0345 = ( ((m_H0
) >= (m_N_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0317)*(IT_0338*IT_0339 + (-0.5)*IT_0317*IT_0342)
       + conj(IT_0338)*(IT_0317*IT_0339 + (-0.5)*IT_0338*IT_0342))*IT_0344
) : 
      (0
));
    const ccomplex_t IT_0346 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0347 = IT_0242*IT_0346;
    const ccomplex_t IT_0348 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0349 = IT_0245*IT_0348;
    const ccomplex_t IT_0350 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0351 = IT_0248*IT_0350;
    const ccomplex_t IT_0352 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0353 = IT_0251*IT_0352;
    const ccomplex_t IT_0354 = (0 + _Complex_I*1)*(IT_0347 + -IT_0349 + 
      -IT_0351 + IT_0353);
    const ccomplex_t IT_0355 = (-0.5)*IT_0354;
    const ccomplex_t IT_0356 = 0.5*IT_0354;
    const ccomplex_t IT_0357 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0358 = IT_0242*IT_0357;
    const ccomplex_t IT_0359 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0360 = IT_0245*IT_0359;
    const ccomplex_t IT_0361 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0362 = IT_0248*IT_0361;
    const ccomplex_t IT_0363 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0364 = IT_0251*IT_0363;
    const ccomplex_t IT_0365 = (0 + _Complex_I*1)*(IT_0358 + -IT_0360 + 
      -IT_0362 + IT_0364);
    const ccomplex_t IT_0366 = (-0.5)*IT_0365;
    const ccomplex_t IT_0367 = 0.5*IT_0365;
    const ccomplex_t IT_0368 = pow(m_N_3, 2);
    const ccomplex_t IT_0369 = 0.5*IT_0018 + -IT_0368;
    const ccomplex_t IT_0370 = -IT_0369;
    const ccomplex_t IT_0371 = ((conj(IT_0355) + -conj(IT_0356))*(IT_0366 + 
      -IT_0367) + (IT_0355 + -IT_0356)*(conj(IT_0366) + -conj(IT_0367)))*IT_0368
       + -(IT_0355*conj(IT_0355) + IT_0356*conj(IT_0356) + IT_0366*conj(IT_0366)
       + IT_0367*conj(IT_0367))*IT_0369 + -(conj(IT_0355)*IT_0356 + IT_0355*conj
      (IT_0356) + conj(IT_0366)*IT_0367 + IT_0366*conj(IT_0367))*IT_0370;
    const ccomplex_t IT_0372 = cabs(0.25*IT_0018 + -IT_0368);
    const ccomplex_t IT_0373 = cpow(IT_0372, 0.5);
    const ccomplex_t IT_0374 = ( ((m_H0
) >= (2*m_N_3
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0371*IT_0373
) : (0
));
    const ccomplex_t IT_0375 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0376 = IT_0242*IT_0375;
    const ccomplex_t IT_0377 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0378 = IT_0251*IT_0377;
    const ccomplex_t IT_0379 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0380 = IT_0242*IT_0379;
    const ccomplex_t IT_0381 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0382 = IT_0251*IT_0381;
    const ccomplex_t IT_0383 = IT_0376 + IT_0378 + IT_0380 + IT_0382;
    const ccomplex_t IT_0384 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0385 = IT_0245*IT_0384;
    const ccomplex_t IT_0386 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0387 = IT_0248*IT_0386;
    const ccomplex_t IT_0388 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0389 = IT_0245*IT_0388;
    const ccomplex_t IT_0390 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0391 = IT_0248*IT_0390;
    const ccomplex_t IT_0392 = -IT_0385 + -IT_0387 + -IT_0389 + -IT_0391;
    const ccomplex_t IT_0393 = IT_0383 + IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*IT_0393;
    const ccomplex_t IT_0395 = 0.5*IT_0394;
    const ccomplex_t IT_0396 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0397 = IT_0242*IT_0396;
    const ccomplex_t IT_0398 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0399 = IT_0251*IT_0398;
    const ccomplex_t IT_0400 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0401 = IT_0242*IT_0400;
    const ccomplex_t IT_0402 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0403 = IT_0251*IT_0402;
    const ccomplex_t IT_0404 = IT_0397 + IT_0399 + IT_0401 + IT_0403;
    const ccomplex_t IT_0405 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0406 = IT_0245*IT_0405;
    const ccomplex_t IT_0407 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0408 = IT_0248*IT_0407;
    const ccomplex_t IT_0409 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0410 = IT_0245*IT_0409;
    const ccomplex_t IT_0411 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0412 = IT_0248*IT_0411;
    const ccomplex_t IT_0413 = -IT_0406 + -IT_0408 + -IT_0410 + -IT_0412;
    const ccomplex_t IT_0414 = IT_0404 + IT_0413;
    const ccomplex_t IT_0415 = (0 + _Complex_I*1)*IT_0414;
    const ccomplex_t IT_0416 = 0.5*IT_0415;
    const ccomplex_t IT_0417 = m_N_1*m_N_3;
    const ccomplex_t IT_0418 = cpow(IT_0018 + IT_0267 + -IT_0368, 2);
    const ccomplex_t IT_0419 = -IT_0267 + 0.25*IT_0001*(IT_0018 + -IT_0267 +
       IT_0368)*(IT_0018 + IT_0267 + -IT_0368) + 0.25*IT_0001*IT_0418;
    const ccomplex_t IT_0420 = 2*IT_0419;
    const ccomplex_t IT_0421 = cabs(-IT_0267 + 0.25*IT_0001*IT_0418);
    const ccomplex_t IT_0422 = cpow(IT_0421, 0.5);
    const ccomplex_t IT_0423 = ( ((m_H0
) >= (m_N_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0395)*(IT_0416*IT_0417 + (-0.5)*IT_0395*IT_0420)
       + conj(IT_0416)*(IT_0395*IT_0417 + (-0.5)*IT_0416*IT_0420))*IT_0422
) : 
      (0
));
    const ccomplex_t IT_0424 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0425 = IT_0242*IT_0424;
    const ccomplex_t IT_0426 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0427 = IT_0251*IT_0426;
    const ccomplex_t IT_0428 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0429 = IT_0242*IT_0428;
    const ccomplex_t IT_0430 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0431 = IT_0251*IT_0430;
    const ccomplex_t IT_0432 = IT_0425 + IT_0427 + IT_0429 + IT_0431;
    const ccomplex_t IT_0433 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0434 = IT_0245*IT_0433;
    const ccomplex_t IT_0435 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0436 = IT_0248*IT_0435;
    const ccomplex_t IT_0437 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0438 = IT_0245*IT_0437;
    const ccomplex_t IT_0439 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0440 = IT_0248*IT_0439;
    const ccomplex_t IT_0441 = -IT_0434 + -IT_0436 + -IT_0438 + -IT_0440;
    const ccomplex_t IT_0442 = IT_0432 + IT_0441;
    const ccomplex_t IT_0443 = (0 + _Complex_I*1)*IT_0442;
    const ccomplex_t IT_0444 = 0.5*IT_0443;
    const ccomplex_t IT_0445 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0446 = IT_0242*IT_0445;
    const ccomplex_t IT_0447 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0448 = IT_0251*IT_0447;
    const ccomplex_t IT_0449 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0450 = IT_0242*IT_0449;
    const ccomplex_t IT_0451 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0452 = IT_0251*IT_0451;
    const ccomplex_t IT_0453 = IT_0446 + IT_0448 + IT_0450 + IT_0452;
    const ccomplex_t IT_0454 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0455 = IT_0245*IT_0454;
    const ccomplex_t IT_0456 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0457 = IT_0248*IT_0456;
    const ccomplex_t IT_0458 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0459 = IT_0245*IT_0458;
    const ccomplex_t IT_0460 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0461 = IT_0248*IT_0460;
    const ccomplex_t IT_0462 = -IT_0455 + -IT_0457 + -IT_0459 + -IT_0461;
    const ccomplex_t IT_0463 = IT_0453 + IT_0462;
    const ccomplex_t IT_0464 = (0 + _Complex_I*1)*IT_0463;
    const ccomplex_t IT_0465 = 0.5*IT_0464;
    const ccomplex_t IT_0466 = m_N_2*m_N_3;
    const ccomplex_t IT_0467 = cpow(IT_0018 + IT_0292 + -IT_0368, 2);
    const ccomplex_t IT_0468 = -IT_0292 + 0.25*IT_0001*(IT_0018 + -IT_0292 +
       IT_0368)*(IT_0018 + IT_0292 + -IT_0368) + 0.25*IT_0001*IT_0467;
    const ccomplex_t IT_0469 = 2*IT_0468;
    const ccomplex_t IT_0470 = cabs(-IT_0292 + 0.25*IT_0001*IT_0467);
    const ccomplex_t IT_0471 = cpow(IT_0470, 0.5);
    const ccomplex_t IT_0472 = ( ((m_H0
) >= (m_N_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0444)*(IT_0465*IT_0466 + (-0.5)*IT_0444*IT_0469)
       + conj(IT_0465)*(IT_0444*IT_0466 + (-0.5)*IT_0465*IT_0469))*IT_0471
) : 
      (0
));
    const ccomplex_t IT_0473 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0474 = IT_0242*IT_0473;
    const ccomplex_t IT_0475 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0476 = IT_0245*IT_0475;
    const ccomplex_t IT_0477 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0478 = IT_0248*IT_0477;
    const ccomplex_t IT_0479 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0480 = IT_0251*IT_0479;
    const ccomplex_t IT_0481 = (0 + _Complex_I*1)*(IT_0474 + -IT_0476 + 
      -IT_0478 + IT_0480);
    const ccomplex_t IT_0482 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0483 = IT_0242*IT_0482;
    const ccomplex_t IT_0484 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0485 = IT_0245*IT_0484;
    const ccomplex_t IT_0486 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0487 = IT_0248*IT_0486;
    const ccomplex_t IT_0488 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0489 = IT_0251*IT_0488;
    const ccomplex_t IT_0490 = (0 + _Complex_I*1)*(IT_0483 + -IT_0485 + 
      -IT_0487 + IT_0489);
    const ccomplex_t IT_0491 = pow(m_N_4, 2);
    const ccomplex_t IT_0492 = 0.5*IT_0018 + -IT_0491;
    const ccomplex_t IT_0493 = cabs(0.25*IT_0018 + -IT_0491);
    const ccomplex_t IT_0494 = cpow(IT_0493, 0.5);
    const ccomplex_t IT_0495 = ( ((m_H0
) >= (2*m_N_4
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0481)*IT_0490 + IT_0481*conj(IT_0490))*IT_0491 + -
      (IT_0481*conj(IT_0481) + IT_0490*conj(IT_0490))*IT_0492)*IT_0494
) : (0
));
    const ccomplex_t IT_0496 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0497 = IT_0242*IT_0496;
    const ccomplex_t IT_0498 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0499 = IT_0251*IT_0498;
    const ccomplex_t IT_0500 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0501 = IT_0242*IT_0500;
    const ccomplex_t IT_0502 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0503 = IT_0251*IT_0502;
    const ccomplex_t IT_0504 = IT_0497 + IT_0499 + IT_0501 + IT_0503;
    const ccomplex_t IT_0505 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0506 = IT_0245*IT_0505;
    const ccomplex_t IT_0507 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0508 = IT_0248*IT_0507;
    const ccomplex_t IT_0509 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0510 = IT_0245*IT_0509;
    const ccomplex_t IT_0511 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0512 = IT_0248*IT_0511;
    const ccomplex_t IT_0513 = -IT_0506 + -IT_0508 + -IT_0510 + -IT_0512;
    const ccomplex_t IT_0514 = IT_0504 + IT_0513;
    const ccomplex_t IT_0515 = (0 + _Complex_I*1)*IT_0514;
    const ccomplex_t IT_0516 = 0.5*IT_0515;
    const ccomplex_t IT_0517 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0518 = IT_0242*IT_0517;
    const ccomplex_t IT_0519 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0520 = IT_0251*IT_0519;
    const ccomplex_t IT_0521 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0522 = IT_0242*IT_0521;
    const ccomplex_t IT_0523 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0524 = IT_0251*IT_0523;
    const ccomplex_t IT_0525 = IT_0518 + IT_0520 + IT_0522 + IT_0524;
    const ccomplex_t IT_0526 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0527 = IT_0245*IT_0526;
    const ccomplex_t IT_0528 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0529 = IT_0248*IT_0528;
    const ccomplex_t IT_0530 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0531 = IT_0245*IT_0530;
    const ccomplex_t IT_0532 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0533 = IT_0248*IT_0532;
    const ccomplex_t IT_0534 = -IT_0527 + -IT_0529 + -IT_0531 + -IT_0533;
    const ccomplex_t IT_0535 = IT_0525 + IT_0534;
    const ccomplex_t IT_0536 = (0 + _Complex_I*1)*IT_0535;
    const ccomplex_t IT_0537 = 0.5*IT_0536;
    const ccomplex_t IT_0538 = m_N_1*m_N_4;
    const ccomplex_t IT_0539 = cpow(IT_0018 + IT_0267 + -IT_0491, 2);
    const ccomplex_t IT_0540 = -IT_0267 + 0.25*IT_0001*(IT_0018 + -IT_0267 +
       IT_0491)*(IT_0018 + IT_0267 + -IT_0491) + 0.25*IT_0001*IT_0539;
    const ccomplex_t IT_0541 = 2*IT_0540;
    const ccomplex_t IT_0542 = cabs(-IT_0267 + 0.25*IT_0001*IT_0539);
    const ccomplex_t IT_0543 = cpow(IT_0542, 0.5);
    const ccomplex_t IT_0544 = ( ((m_H0
) >= (m_N_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0516)*(IT_0537*IT_0538 + (-0.5)*IT_0516*IT_0541)
       + conj(IT_0537)*(IT_0516*IT_0538 + (-0.5)*IT_0537*IT_0541))*IT_0543
) : 
      (0
));
    const ccomplex_t IT_0545 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0546 = IT_0242*IT_0545;
    const ccomplex_t IT_0547 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0548 = IT_0251*IT_0547;
    const ccomplex_t IT_0549 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0550 = IT_0242*IT_0549;
    const ccomplex_t IT_0551 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0552 = IT_0251*IT_0551;
    const ccomplex_t IT_0553 = IT_0546 + IT_0548 + IT_0550 + IT_0552;
    const ccomplex_t IT_0554 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0555 = IT_0245*IT_0554;
    const ccomplex_t IT_0556 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0557 = IT_0248*IT_0556;
    const ccomplex_t IT_0558 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0559 = IT_0245*IT_0558;
    const ccomplex_t IT_0560 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0561 = IT_0248*IT_0560;
    const ccomplex_t IT_0562 = -IT_0555 + -IT_0557 + -IT_0559 + -IT_0561;
    const ccomplex_t IT_0563 = IT_0553 + IT_0562;
    const ccomplex_t IT_0564 = (0 + _Complex_I*1)*IT_0563;
    const ccomplex_t IT_0565 = 0.5*IT_0564;
    const ccomplex_t IT_0566 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0567 = IT_0242*IT_0566;
    const ccomplex_t IT_0568 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0569 = IT_0251*IT_0568;
    const ccomplex_t IT_0570 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0571 = IT_0242*IT_0570;
    const ccomplex_t IT_0572 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0573 = IT_0251*IT_0572;
    const ccomplex_t IT_0574 = IT_0567 + IT_0569 + IT_0571 + IT_0573;
    const ccomplex_t IT_0575 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0576 = IT_0245*IT_0575;
    const ccomplex_t IT_0577 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0578 = IT_0248*IT_0577;
    const ccomplex_t IT_0579 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0580 = IT_0245*IT_0579;
    const ccomplex_t IT_0581 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0582 = IT_0248*IT_0581;
    const ccomplex_t IT_0583 = -IT_0576 + -IT_0578 + -IT_0580 + -IT_0582;
    const ccomplex_t IT_0584 = IT_0574 + IT_0583;
    const ccomplex_t IT_0585 = (0 + _Complex_I*1)*IT_0584;
    const ccomplex_t IT_0586 = 0.5*IT_0585;
    const ccomplex_t IT_0587 = m_N_2*m_N_4;
    const ccomplex_t IT_0588 = cpow(IT_0018 + IT_0292 + -IT_0491, 2);
    const ccomplex_t IT_0589 = -IT_0292 + 0.25*IT_0001*(IT_0018 + -IT_0292 +
       IT_0491)*(IT_0018 + IT_0292 + -IT_0491) + 0.25*IT_0001*IT_0588;
    const ccomplex_t IT_0590 = 2*IT_0589;
    const ccomplex_t IT_0591 = cabs(-IT_0292 + 0.25*IT_0001*IT_0588);
    const ccomplex_t IT_0592 = cpow(IT_0591, 0.5);
    const ccomplex_t IT_0593 = ( ((m_H0
) >= (m_N_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0565)*(IT_0586*IT_0587 + (-0.5)*IT_0565*IT_0590)
       + conj(IT_0586)*(IT_0565*IT_0587 + (-0.5)*IT_0586*IT_0590))*IT_0592
) : 
      (0
));
    const ccomplex_t IT_0594 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0595 = IT_0242*IT_0594;
    const ccomplex_t IT_0596 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0597 = IT_0251*IT_0596;
    const ccomplex_t IT_0598 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0599 = IT_0242*IT_0598;
    const ccomplex_t IT_0600 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0601 = IT_0251*IT_0600;
    const ccomplex_t IT_0602 = IT_0595 + IT_0597 + IT_0599 + IT_0601;
    const ccomplex_t IT_0603 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0604 = IT_0245*IT_0603;
    const ccomplex_t IT_0605 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0606 = IT_0248*IT_0605;
    const ccomplex_t IT_0607 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0608 = IT_0245*IT_0607;
    const ccomplex_t IT_0609 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0610 = IT_0248*IT_0609;
    const ccomplex_t IT_0611 = -IT_0604 + -IT_0606 + -IT_0608 + -IT_0610;
    const ccomplex_t IT_0612 = IT_0602 + IT_0611;
    const ccomplex_t IT_0613 = (0 + _Complex_I*1)*IT_0612;
    const ccomplex_t IT_0614 = 0.5*IT_0613;
    const ccomplex_t IT_0615 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0616 = IT_0242*IT_0615;
    const ccomplex_t IT_0617 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0618 = IT_0251*IT_0617;
    const ccomplex_t IT_0619 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0620 = IT_0242*IT_0619;
    const ccomplex_t IT_0621 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0622 = IT_0251*IT_0621;
    const ccomplex_t IT_0623 = IT_0616 + IT_0618 + IT_0620 + IT_0622;
    const ccomplex_t IT_0624 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0625 = IT_0245*IT_0624;
    const ccomplex_t IT_0626 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0627 = IT_0248*IT_0626;
    const ccomplex_t IT_0628 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0629 = IT_0245*IT_0628;
    const ccomplex_t IT_0630 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0631 = IT_0248*IT_0630;
    const ccomplex_t IT_0632 = -IT_0625 + -IT_0627 + -IT_0629 + -IT_0631;
    const ccomplex_t IT_0633 = IT_0623 + IT_0632;
    const ccomplex_t IT_0634 = (0 + _Complex_I*1)*IT_0633;
    const ccomplex_t IT_0635 = 0.5*IT_0634;
    const ccomplex_t IT_0636 = m_N_3*m_N_4;
    const ccomplex_t IT_0637 = cpow(IT_0018 + IT_0368 + -IT_0491, 2);
    const ccomplex_t IT_0638 = -IT_0368 + 0.25*IT_0001*(IT_0018 + -IT_0368 +
       IT_0491)*(IT_0018 + IT_0368 + -IT_0491) + 0.25*IT_0001*IT_0637;
    const ccomplex_t IT_0639 = 2*IT_0638;
    const ccomplex_t IT_0640 = cabs(-IT_0368 + 0.25*IT_0001*IT_0637);
    const ccomplex_t IT_0641 = cpow(IT_0640, 0.5);
    const ccomplex_t IT_0642 = ( ((m_H0
) >= (m_N_3 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0614)*(IT_0635*IT_0636 + (-0.5)*IT_0614*IT_0639)
       + conj(IT_0635)*(IT_0614*IT_0636 + (-0.5)*IT_0635*IT_0639))*IT_0641
) : 
      (0
));
    const ccomplex_t IT_0643 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0011
      *IT_0022*IT_0023;
    const ccomplex_t IT_0644 = (-0.5)*IT_0643;
    const ccomplex_t IT_0645 = cpow(IT_0644, 2);
    const ccomplex_t IT_0646 = pow(m_tau, 2);
    const ccomplex_t IT_0647 = 0.5*IT_0018 + -IT_0646;
    const ccomplex_t IT_0648 = cabs(0.25*IT_0018 + -IT_0646);
    const ccomplex_t IT_0649 = cpow(IT_0648, 0.5);
    const ccomplex_t IT_0650 = ( ((m_H0
) >= (2*m_tau
)) ? (0.5*IT_0000
      *IT_0001*IT_0645*(IT_0646 + -IT_0647)*IT_0649
) : (0
));
    const ccomplex_t IT_0651 = cpow(V_tb, 2);
    const ccomplex_t IT_0652 = m_b*mu_h*IT_0022;
    const ccomplex_t IT_0653 = -IT_0008*IT_0009;
    const ccomplex_t IT_0654 = IT_0011*IT_0012;
    const ccomplex_t IT_0655 = IT_0653 + IT_0654;
    const ccomplex_t IT_0656 = (-2)*U_sb_00*((IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0027 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0655)
      *conj(U_sb_00) + -1./2*IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_10)) +
       U_sb_10*(IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_00) + (-2)*(IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0027 + -1./3*m_W*IT_0003*IT_0004*IT_0655)*conj
      (U_sb_10));
    const ccomplex_t IT_0657 = (0 + _Complex_I*-1)*(e_em*IT_0656 + (
      -1.4142135623731)*A_b*IT_0011*IT_0651*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0658 = (-0.5)*IT_0657;
    const ccomplex_t IT_0659 = pow(m_sb_1, 2);
    const ccomplex_t IT_0660 = cabs(0.25*IT_0018 + -IT_0659);
    const ccomplex_t IT_0661 = cpow(IT_0660, 0.5);
    const ccomplex_t IT_0662 = ( ((m_H0
) >= (2*m_sb_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0658*conj(IT_0658)*IT_0661
) : (0
));
    const ccomplex_t IT_0663 = (-2)*U_sb_01*((IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0027 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0655)
      *conj(U_sb_01) + -1./2*IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_11)) +
       U_sb_11*(IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_01) + (-2)*(IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0027 + -1./3*m_W*IT_0003*IT_0004*IT_0655)*conj
      (U_sb_11));
    const ccomplex_t IT_0664 = (0 + _Complex_I*-1)*(e_em*IT_0663 + (
      -1.4142135623731)*A_b*IT_0011*IT_0651*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0665 = (-0.5)*IT_0664;
    const ccomplex_t IT_0666 = pow(m_sb_2, 2);
    const ccomplex_t IT_0667 = cabs(0.25*IT_0018 + -IT_0666);
    const ccomplex_t IT_0668 = cpow(IT_0667, 0.5);
    const ccomplex_t IT_0669 = ( ((m_H0
) >= (2*m_sb_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0665*conj(IT_0665)*IT_0668
) : (0
));
    const ccomplex_t IT_0670 = (-2)*U_sb_00*((IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0027 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0655)
      *conj(U_sb_01) + -1./2*IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_11)) +
       U_sb_10*(IT_0006*IT_0009*IT_0023*IT_0652*conj(U_sb_01) + (-2)*(IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0027 + -1./3*m_W*IT_0003*IT_0004*IT_0655)*conj
      (U_sb_11));
    const ccomplex_t IT_0671 = (0 + _Complex_I*-1)*(e_em*IT_0670 + (
      -1.4142135623731)*A_b*IT_0011*IT_0651*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0672 = (-0.5)*IT_0671;
    const ccomplex_t IT_0673 = cpow(IT_0018 + IT_0659 + -IT_0666, 2);
    const ccomplex_t IT_0674 = cabs(-IT_0659 + 0.25*IT_0001*IT_0673);
    const ccomplex_t IT_0675 = cpow(IT_0674, 0.5);
    const ccomplex_t IT_0676 = ( ((m_H0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0672*conj(IT_0672)*IT_0675
) : (0
));
    const ccomplex_t IT_0677 = (-2)*conj(U_sb_00)*((IT_0006*IT_0011*IT_0022
      *IT_0023*IT_0027 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)
      *IT_0655)*U_sb_01 + -1./2*IT_0006*IT_0009*IT_0023*IT_0652*U_sb_11) + conj
      (U_sb_10)*(IT_0006*IT_0009*IT_0023*IT_0652*U_sb_01 + (-2)*(IT_0006*IT_0011
      *IT_0022*IT_0023*IT_0027 + -1./3*m_W*IT_0003*IT_0004*IT_0655)*U_sb_11);
    const ccomplex_t IT_0678 = (0 + _Complex_I*-1)*(e_em*IT_0677 + (
      -1.4142135623731)*A_b*IT_0011*IT_0651*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0679 = (-0.5)*IT_0678;
    const ccomplex_t IT_0680 = ( ((m_H0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0675*IT_0679*conj(IT_0679)
) : (0
));
    const ccomplex_t IT_0681 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0022*IT_0032*IT_0036);
    const ccomplex_t IT_0682 = -IT_0681;
    const ccomplex_t IT_0683 = cpow(IT_0682, 2);
    const ccomplex_t IT_0684 = pow(m_sc_L, 2);
    const ccomplex_t IT_0685 = cabs(0.25*IT_0018 + -IT_0684);
    const ccomplex_t IT_0686 = cpow(IT_0685, 0.5);
    const ccomplex_t IT_0687 = ( ((m_H0
) >= (2*m_sc_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0683*IT_0686
) : (0
));
    const ccomplex_t IT_0688 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0009*IT_0022
      *IT_0032*IT_0036);
    const ccomplex_t IT_0689 = (-0.666666666666667)*IT_0688;
    const ccomplex_t IT_0690 = cpow(IT_0689, 2);
    const ccomplex_t IT_0691 = pow(m_sc_R, 2);
    const ccomplex_t IT_0692 = cabs(0.25*IT_0018 + -IT_0691);
    const ccomplex_t IT_0693 = cpow(IT_0692, 0.5);
    const ccomplex_t IT_0694 = ( ((m_H0
) >= (2*m_sc_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0690*IT_0693
) : (0
));
    const ccomplex_t IT_0695 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0006
      *IT_0011*IT_0022*IT_0032;
    const ccomplex_t IT_0696 = 0.5*IT_0695;
    const ccomplex_t IT_0697 = cpow(IT_0696, 2);
    const ccomplex_t IT_0698 = cpow(IT_0018 + IT_0684 + -IT_0691, 2);
    const ccomplex_t IT_0699 = cabs(-IT_0684 + 0.25*IT_0001*IT_0698);
    const ccomplex_t IT_0700 = cpow(IT_0699, 0.5);
    const ccomplex_t IT_0701 = ( ((m_H0
) >= (m_sc_L + m_sc_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0697*IT_0700
) : (0
));
    const ccomplex_t IT_0702 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0022*IT_0023*IT_0044);
    const ccomplex_t IT_0703 = -IT_0702;
    const ccomplex_t IT_0704 = cpow(IT_0703, 2);
    const ccomplex_t IT_0705 = pow(m_sd_L, 2);
    const ccomplex_t IT_0706 = cabs(0.25*IT_0018 + -IT_0705);
    const ccomplex_t IT_0707 = cpow(IT_0706, 0.5);
    const ccomplex_t IT_0708 = ( ((m_H0
) >= (2*m_sd_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0704*IT_0707
) : (0
));
    const ccomplex_t IT_0709 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_10);
    const ccomplex_t IT_0710 = (-0.5)*IT_0709;
    const ccomplex_t IT_0711 = cpow(IT_0018 + -IT_0659 + IT_0705, 2);
    const ccomplex_t IT_0712 = cabs(-IT_0705 + 0.25*IT_0001*IT_0711);
    const ccomplex_t IT_0713 = cpow(IT_0712, 0.5);
    const ccomplex_t IT_0714 = ( ((m_H0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0710*conj(IT_0710)*IT_0713
) : (0
));
    const ccomplex_t IT_0715 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_10;
    const ccomplex_t IT_0716 = (-0.5)*IT_0715;
    const ccomplex_t IT_0717 = ( ((m_H0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0713*IT_0716*conj(IT_0716)
) : (0
));
    const ccomplex_t IT_0718 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0719 = (-0.5)*IT_0718;
    const ccomplex_t IT_0720 = cpow(IT_0018 + -IT_0666 + IT_0705, 2);
    const ccomplex_t IT_0721 = cabs(-IT_0705 + 0.25*IT_0001*IT_0720);
    const ccomplex_t IT_0722 = cpow(IT_0721, 0.5);
    const ccomplex_t IT_0723 = ( ((m_H0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0719*conj(IT_0719)*IT_0722
) : (0
));
    const ccomplex_t IT_0724 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_11;
    const ccomplex_t IT_0725 = (-0.5)*IT_0724;
    const ccomplex_t IT_0726 = ( ((m_H0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0722*IT_0725*conj(IT_0725)
) : (0
));
    const ccomplex_t IT_0727 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + 3*IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0044);
    const ccomplex_t IT_0728 = 0.333333333333333*IT_0727;
    const ccomplex_t IT_0729 = cpow(IT_0728, 2);
    const ccomplex_t IT_0730 = pow(m_sd_R, 2);
    const ccomplex_t IT_0731 = cabs(0.25*IT_0018 + -IT_0730);
    const ccomplex_t IT_0732 = cpow(IT_0731, 0.5);
    const ccomplex_t IT_0733 = ( ((m_H0
) >= (2*m_sd_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0729*IT_0732
) : (0
));
    const ccomplex_t IT_0734 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_00);
    const ccomplex_t IT_0735 = (-0.5)*IT_0734;
    const ccomplex_t IT_0736 = cpow(IT_0018 + -IT_0659 + IT_0730, 2);
    const ccomplex_t IT_0737 = cabs(-IT_0730 + 0.25*IT_0001*IT_0736);
    const ccomplex_t IT_0738 = cpow(IT_0737, 0.5);
    const ccomplex_t IT_0739 = ( ((m_H0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0735*conj(IT_0735)*IT_0738
) : (0
));
    const ccomplex_t IT_0740 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_00;
    const ccomplex_t IT_0741 = (-0.5)*IT_0740;
    const ccomplex_t IT_0742 = ( ((m_H0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0738*IT_0741*conj(IT_0741)
) : (0
));
    const ccomplex_t IT_0743 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0744 = (-0.5)*IT_0743;
    const ccomplex_t IT_0745 = cpow(IT_0018 + -IT_0666 + IT_0730, 2);
    const ccomplex_t IT_0746 = cabs(-IT_0730 + 0.25*IT_0001*IT_0745);
    const ccomplex_t IT_0747 = cpow(IT_0746, 0.5);
    const ccomplex_t IT_0748 = ( ((m_H0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0744*conj(IT_0744)*IT_0747
) : (0
));
    const ccomplex_t IT_0749 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_01;
    const ccomplex_t IT_0750 = (-0.5)*IT_0749;
    const ccomplex_t IT_0751 = ( ((m_H0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0747*IT_0750*conj(IT_0750)
) : (0
));
    const ccomplex_t IT_0752 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0011 + -m_d*e_em*mu_h*IT_0006*IT_0009*IT_0022*IT_0023);
    const ccomplex_t IT_0753 = (-0.5)*IT_0752;
    const ccomplex_t IT_0754 = cpow(IT_0753, 2);
    const ccomplex_t IT_0755 = cpow(IT_0018 + -IT_0705 + IT_0730, 2);
    const ccomplex_t IT_0756 = cabs(-IT_0730 + 0.25*IT_0001*IT_0755);
    const ccomplex_t IT_0757 = cpow(IT_0756, 0.5);
    const ccomplex_t IT_0758 = ( ((m_H0
) >= (m_sd_L + m_sd_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0754*IT_0757
) : (0
));
    const ccomplex_t IT_0759 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0052);
    const ccomplex_t IT_0760 = -IT_0759;
    const ccomplex_t IT_0761 = cpow(IT_0760, 2);
    const ccomplex_t IT_0762 = pow(m_se_L, 2);
    const ccomplex_t IT_0763 = cabs(0.25*IT_0018 + -IT_0762);
    const ccomplex_t IT_0764 = cpow(IT_0763, 0.5);
    const ccomplex_t IT_0765 = ( ((m_H0
) >= (2*m_se_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0761*IT_0764
) : (0
));
    const ccomplex_t IT_0766 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0052);
    const ccomplex_t IT_0767 = cpow(IT_0766, 2);
    const ccomplex_t IT_0768 = pow(m_se_R, 2);
    const ccomplex_t IT_0769 = cabs(0.25*IT_0018 + -IT_0768);
    const ccomplex_t IT_0770 = cpow(IT_0769, 0.5);
    const ccomplex_t IT_0771 = ( ((m_H0
) >= (2*m_se_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0767*IT_0770
) : (0
));
    const ccomplex_t IT_0772 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0006
      *IT_0009*IT_0022*IT_0023;
    const ccomplex_t IT_0773 = 0.5*IT_0772;
    const ccomplex_t IT_0774 = cpow(IT_0773, 2);
    const ccomplex_t IT_0775 = cpow(IT_0018 + -IT_0762 + IT_0768, 2);
    const ccomplex_t IT_0776 = cabs(-IT_0768 + 0.25*IT_0001*IT_0775);
    const ccomplex_t IT_0777 = cpow(IT_0776, 0.5);
    const ccomplex_t IT_0778 = ( ((m_H0
) >= (m_se_L + m_se_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0774*IT_0777
) : (0
));
    const ccomplex_t IT_0779 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0022*IT_0023*IT_0067);
    const ccomplex_t IT_0780 = -IT_0779;
    const ccomplex_t IT_0781 = cpow(IT_0780, 2);
    const ccomplex_t IT_0782 = pow(m_ss_L, 2);
    const ccomplex_t IT_0783 = cabs(0.25*IT_0018 + -IT_0782);
    const ccomplex_t IT_0784 = cpow(IT_0783, 0.5);
    const ccomplex_t IT_0785 = ( ((m_H0
) >= (2*m_ss_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0781*IT_0784
) : (0
));
    const ccomplex_t IT_0786 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_10);
    const ccomplex_t IT_0787 = (-0.5)*IT_0786;
    const ccomplex_t IT_0788 = cpow(IT_0018 + -IT_0659 + IT_0782, 2);
    const ccomplex_t IT_0789 = cabs(-IT_0782 + 0.25*IT_0001*IT_0788);
    const ccomplex_t IT_0790 = cpow(IT_0789, 0.5);
    const ccomplex_t IT_0791 = ( ((m_H0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0787*conj(IT_0787)*IT_0790
) : (0
));
    const ccomplex_t IT_0792 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_10;
    const ccomplex_t IT_0793 = (-0.5)*IT_0792;
    const ccomplex_t IT_0794 = ( ((m_H0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0790*IT_0793*conj(IT_0793)
) : (0
));
    const ccomplex_t IT_0795 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0796 = (-0.5)*IT_0795;
    const ccomplex_t IT_0797 = cpow(IT_0018 + -IT_0666 + IT_0782, 2);
    const ccomplex_t IT_0798 = cabs(-IT_0782 + 0.25*IT_0001*IT_0797);
    const ccomplex_t IT_0799 = cpow(IT_0798, 0.5);
    const ccomplex_t IT_0800 = ( ((m_H0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0796*conj(IT_0796)*IT_0799
) : (0
));
    const ccomplex_t IT_0801 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_11;
    const ccomplex_t IT_0802 = (-0.5)*IT_0801;
    const ccomplex_t IT_0803 = ( ((m_H0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0799*IT_0802*conj(IT_0802)
) : (0
));
    const ccomplex_t IT_0804 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0011;
    const ccomplex_t IT_0805 = (-0.5)*IT_0804;
    const ccomplex_t IT_0806 = cpow(IT_0018 + IT_0730 + -IT_0782, 2);
    const ccomplex_t IT_0807 = cabs(-IT_0730 + 0.25*IT_0001*IT_0806);
    const ccomplex_t IT_0808 = cpow(IT_0807, 0.5);
    const ccomplex_t IT_0809 = ( ((m_H0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0805*conj(IT_0805)*IT_0808
) : (0
));
    const ccomplex_t IT_0810 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0011;
    const ccomplex_t IT_0811 = (-0.5)*IT_0810;
    const ccomplex_t IT_0812 = ( ((m_H0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0808*IT_0811*conj(IT_0811)
) : (0
));
    const ccomplex_t IT_0813 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + 3*IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0067);
    const ccomplex_t IT_0814 = 0.333333333333333*IT_0813;
    const ccomplex_t IT_0815 = cpow(IT_0814, 2);
    const ccomplex_t IT_0816 = pow(m_ss_R, 2);
    const ccomplex_t IT_0817 = cabs(0.25*IT_0018 + -IT_0816);
    const ccomplex_t IT_0818 = cpow(IT_0817, 0.5);
    const ccomplex_t IT_0819 = ( ((m_H0
) >= (2*m_ss_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0815*IT_0818
) : (0
));
    const ccomplex_t IT_0820 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_00);
    const ccomplex_t IT_0821 = (-0.5)*IT_0820;
    const ccomplex_t IT_0822 = cpow(IT_0018 + -IT_0659 + IT_0816, 2);
    const ccomplex_t IT_0823 = cabs(-IT_0816 + 0.25*IT_0001*IT_0822);
    const ccomplex_t IT_0824 = cpow(IT_0823, 0.5);
    const ccomplex_t IT_0825 = ( ((m_H0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0821*conj(IT_0821)*IT_0824
) : (0
));
    const ccomplex_t IT_0826 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_00;
    const ccomplex_t IT_0827 = (-0.5)*IT_0826;
    const ccomplex_t IT_0828 = ( ((m_H0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0824*IT_0827*conj(IT_0827)
) : (0
));
    const ccomplex_t IT_0829 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0830 = (-0.5)*IT_0829;
    const ccomplex_t IT_0831 = cpow(IT_0018 + -IT_0666 + IT_0816, 2);
    const ccomplex_t IT_0832 = cabs(-IT_0816 + 0.25*IT_0001*IT_0831);
    const ccomplex_t IT_0833 = cpow(IT_0832, 0.5);
    const ccomplex_t IT_0834 = ( ((m_H0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0830*conj(IT_0830)*IT_0833
) : (0
));
    const ccomplex_t IT_0835 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_01;
    const ccomplex_t IT_0836 = (-0.5)*IT_0835;
    const ccomplex_t IT_0837 = ( ((m_H0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0833*IT_0836*conj(IT_0836)
) : (0
));
    const ccomplex_t IT_0838 = cpow(IT_0018 + -IT_0705 + IT_0816, 2);
    const ccomplex_t IT_0839 = cabs(-IT_0816 + 0.25*IT_0001*IT_0838);
    const ccomplex_t IT_0840 = cpow(IT_0839, 0.5);
    const ccomplex_t IT_0841 = ( ((m_H0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0805*conj(IT_0805)*IT_0840
) : (0
));
    const ccomplex_t IT_0842 = ( ((m_H0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0811*conj(IT_0811)*IT_0840
) : (0
));
    const ccomplex_t IT_0843 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0011 + -m_s*e_em*mu_h*IT_0006*IT_0009*IT_0022*IT_0023);
    const ccomplex_t IT_0844 = (-0.5)*IT_0843;
    const ccomplex_t IT_0845 = cpow(IT_0844, 2);
    const ccomplex_t IT_0846 = cpow(IT_0018 + -IT_0782 + IT_0816, 2);
    const ccomplex_t IT_0847 = cabs(-IT_0816 + 0.25*IT_0001*IT_0846);
    const ccomplex_t IT_0848 = cpow(IT_0847, 0.5);
    const ccomplex_t IT_0849 = ( ((m_H0
) >= (m_ss_L + m_ss_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0845*IT_0848
) : (0
));
    const ccomplex_t IT_0850 = m_t*mu_h*IT_0022;
    const ccomplex_t IT_0851 = -IT_0008*IT_0009 + IT_0654;
    const ccomplex_t IT_0852 = (-2)*U_st_00*((IT_0006*IT_0009*IT_0022*IT_0032
      *IT_0075 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0851)*conj(U_st_00) + -1./2*IT_0006*IT_0011*IT_0032*IT_0850*conj
      (U_st_10)) + U_st_10*(IT_0006*IT_0011*IT_0032*IT_0850*conj(U_st_00) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0022*IT_0032*IT_0075 + m_W
      *IT_0003*IT_0004*IT_0851)*conj(U_st_10));
    const ccomplex_t IT_0853 = (0 + _Complex_I*-1)*(e_em*IT_0852 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0854 = (-0.5)*IT_0853;
    const ccomplex_t IT_0855 = pow(m_st_1, 2);
    const ccomplex_t IT_0856 = cabs(0.25*IT_0018 + -IT_0855);
    const ccomplex_t IT_0857 = cpow(IT_0856, 0.5);
    const ccomplex_t IT_0858 = ( ((m_H0
) >= (2*m_st_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0854*conj(IT_0854)*IT_0857
) : (0
));
    const ccomplex_t IT_0859 = (-2)*U_st_01*((IT_0006*IT_0009*IT_0022*IT_0032
      *IT_0075 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0851)*conj(U_st_01) + -1./2*IT_0006*IT_0011*IT_0032*IT_0850*conj
      (U_st_11)) + U_st_11*(IT_0006*IT_0011*IT_0032*IT_0850*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0022*IT_0032*IT_0075 + m_W
      *IT_0003*IT_0004*IT_0851)*conj(U_st_11));
    const ccomplex_t IT_0860 = (0 + _Complex_I*-1)*(e_em*IT_0859 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0861 = (-0.5)*IT_0860;
    const ccomplex_t IT_0862 = pow(m_st_2, 2);
    const ccomplex_t IT_0863 = cabs(0.25*IT_0018 + -IT_0862);
    const ccomplex_t IT_0864 = cpow(IT_0863, 0.5);
    const ccomplex_t IT_0865 = ( ((m_H0
) >= (2*m_st_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0861*conj(IT_0861)*IT_0864
) : (0
));
    const ccomplex_t IT_0866 = (-2)*U_st_00*((IT_0006*IT_0009*IT_0022*IT_0032
      *IT_0075 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)*IT_0006)
      *IT_0851)*conj(U_st_01) + -1./2*IT_0006*IT_0011*IT_0032*IT_0850*conj
      (U_st_11)) + U_st_10*(IT_0006*IT_0011*IT_0032*IT_0850*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0022*IT_0032*IT_0075 + m_W
      *IT_0003*IT_0004*IT_0851)*conj(U_st_11));
    const ccomplex_t IT_0867 = (0 + _Complex_I*-1)*(e_em*IT_0866 + (
      -1.4142135623731)*A_t*IT_0009*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0868 = (-0.5)*IT_0867;
    const ccomplex_t IT_0869 = cpow(IT_0018 + IT_0855 + -IT_0862, 2);
    const ccomplex_t IT_0870 = cabs(-IT_0855 + 0.25*IT_0001*IT_0869);
    const ccomplex_t IT_0871 = cpow(IT_0870, 0.5);
    const ccomplex_t IT_0872 = ( ((m_H0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0868*conj(IT_0868)*IT_0871
) : (0
));
    const ccomplex_t IT_0873 = (-2)*conj(U_st_00)*((IT_0006*IT_0009*IT_0022
      *IT_0032*IT_0075 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + (-3)
      *IT_0006)*IT_0851)*U_st_01 + -1./2*IT_0006*IT_0011*IT_0032*IT_0850*U_st_11
      ) + conj(U_st_10)*(IT_0006*IT_0011*IT_0032*IT_0850*U_st_01 + (
      -1.33333333333333)*(1.5*IT_0006*IT_0009*IT_0022*IT_0032*IT_0075 + m_W
      *IT_0003*IT_0004*IT_0851)*U_st_11);
    const ccomplex_t IT_0874 = (0 + _Complex_I*-1)*(e_em*IT_0873 + (
      -1.4142135623731)*A_t*IT_0009*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0875 = (-0.5)*IT_0874;
    const ccomplex_t IT_0876 = ( ((m_H0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0871*IT_0875*conj(IT_0875)
) : (0
));
    const ccomplex_t IT_0877 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0022*IT_0032*IT_0083);
    const ccomplex_t IT_0878 = -IT_0877;
    const ccomplex_t IT_0879 = cpow(IT_0878, 2);
    const ccomplex_t IT_0880 = pow(m_su_L, 2);
    const ccomplex_t IT_0881 = cabs(0.25*IT_0018 + -IT_0880);
    const ccomplex_t IT_0882 = cpow(IT_0881, 0.5);
    const ccomplex_t IT_0883 = ( ((m_H0
) >= (2*m_su_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0879*IT_0882
) : (0
));
    const ccomplex_t IT_0884 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0009*IT_0022
      *IT_0032*IT_0083);
    const ccomplex_t IT_0885 = (-0.666666666666667)*IT_0884;
    const ccomplex_t IT_0886 = cpow(IT_0885, 2);
    const ccomplex_t IT_0887 = pow(m_su_R, 2);
    const ccomplex_t IT_0888 = cabs(0.25*IT_0018 + -IT_0887);
    const ccomplex_t IT_0889 = cpow(IT_0888, 0.5);
    const ccomplex_t IT_0890 = ( ((m_H0
) >= (2*m_su_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0886*IT_0889
) : (0
));
    const ccomplex_t IT_0891 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0006
      *IT_0011*IT_0022*IT_0032;
    const ccomplex_t IT_0892 = 0.5*IT_0891;
    const ccomplex_t IT_0893 = cpow(IT_0892, 2);
    const ccomplex_t IT_0894 = cpow(IT_0018 + IT_0880 + -IT_0887, 2);
    const ccomplex_t IT_0895 = cabs(-IT_0880 + 0.25*IT_0001*IT_0894);
    const ccomplex_t IT_0896 = cpow(IT_0895, 0.5);
    const ccomplex_t IT_0897 = ( ((m_H0
) >= (m_su_L + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0893*IT_0896
) : (0
));
    const ccomplex_t IT_0898 = (0 + _Complex_I*1)*e_em*((-0.5)*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0144);
    const ccomplex_t IT_0899 = -IT_0898;
    const ccomplex_t IT_0900 = cpow(IT_0899, 2);
    const ccomplex_t IT_0901 = pow(m_smu_L, 2);
    const ccomplex_t IT_0902 = cabs(0.25*IT_0018 + -IT_0901);
    const ccomplex_t IT_0903 = cpow(IT_0902, 0.5);
    const ccomplex_t IT_0904 = ( ((m_H0
) >= (2*m_smu_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0900*IT_0903
) : (0
));
    const ccomplex_t IT_0905 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + -IT_0011*IT_0012) + IT_0006*IT_0011*IT_0022*IT_0023
      *IT_0144);
    const ccomplex_t IT_0906 = cpow(IT_0905, 2);
    const ccomplex_t IT_0907 = pow(m_smu_R, 2);
    const ccomplex_t IT_0908 = cabs(0.25*IT_0018 + -IT_0907);
    const ccomplex_t IT_0909 = cpow(IT_0908, 0.5);
    const ccomplex_t IT_0910 = ( ((m_H0
) >= (2*m_smu_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0906*IT_0909
) : (0
));
    const ccomplex_t IT_0911 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0006
      *IT_0009*IT_0022*IT_0023;
    const ccomplex_t IT_0912 = 0.5*IT_0911;
    const ccomplex_t IT_0913 = cpow(IT_0912, 2);
    const ccomplex_t IT_0914 = cpow(IT_0018 + -IT_0901 + IT_0907, 2);
    const ccomplex_t IT_0915 = cabs(-IT_0907 + 0.25*IT_0001*IT_0914);
    const ccomplex_t IT_0916 = cpow(IT_0915, 0.5);
    const ccomplex_t IT_0917 = ( ((m_H0
) >= (m_smu_L + m_smu_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0913*IT_0916
) : (0
));
    const ccomplex_t IT_0918 = (0 + _Complex_I*-1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0009 + -IT_0011*IT_0012);
    const ccomplex_t IT_0919 = (-0.5)*IT_0918;
    const ccomplex_t IT_0920 = cpow(IT_0919, 2);
    const ccomplex_t IT_0921 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0922 = cabs(0.25*IT_0018 + -IT_0921);
    const ccomplex_t IT_0923 = cpow(IT_0922, 0.5);
    const ccomplex_t IT_0924 = ( ((m_H0
) >= (2*m_snu_mu
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0920*IT_0923
) : (0
));
    const ccomplex_t IT_0925 = mu_h*m_tau*IT_0022;
    const ccomplex_t IT_0926 = U_stau_10*(IT_0006*IT_0009*IT_0023*IT_0925*conj
      (U_stau_00) + (-2)*(IT_0006*IT_0011*IT_0022*IT_0023*IT_0646 + -m_W*IT_0003
      *IT_0004*IT_0655)*conj(U_stau_10)) + (-2)*U_stau_00*((IT_0006*IT_0011
      *IT_0022*IT_0023*IT_0646 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0655)
      *conj(U_stau_00) + -1./2*IT_0006*IT_0009*IT_0023*IT_0925*conj(U_stau_10));
    const ccomplex_t IT_0927 = (0 + _Complex_I*-1)*(e_em*IT_0926 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0928 = (-0.5)*IT_0927;
    const ccomplex_t IT_0929 = pow(m_stau_1, 2);
    const ccomplex_t IT_0930 = cabs(0.25*IT_0018 + -IT_0929);
    const ccomplex_t IT_0931 = cpow(IT_0930, 0.5);
    const ccomplex_t IT_0932 = ( ((m_H0
) >= (2*m_stau_1
)) ? (0.125*IT_0000
      *IT_0001*IT_0928*conj(IT_0928)*IT_0931
) : (0
));
    const ccomplex_t IT_0933 = U_stau_11*(IT_0006*IT_0009*IT_0023*IT_0925*conj
      (U_stau_01) + (-2)*(IT_0006*IT_0011*IT_0022*IT_0023*IT_0646 + -m_W*IT_0003
      *IT_0004*IT_0655)*conj(U_stau_11)) + (-2)*U_stau_01*((IT_0006*IT_0011
      *IT_0022*IT_0023*IT_0646 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0655)
      *conj(U_stau_01) + -1./2*IT_0006*IT_0009*IT_0023*IT_0925*conj(U_stau_11));
    const ccomplex_t IT_0934 = (0 + _Complex_I*-1)*(e_em*IT_0933 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0935 = (-0.5)*IT_0934;
    const ccomplex_t IT_0936 = pow(m_stau_2, 2);
    const ccomplex_t IT_0937 = cabs(0.25*IT_0018 + -IT_0936);
    const ccomplex_t IT_0938 = cpow(IT_0937, 0.5);
    const ccomplex_t IT_0939 = ( ((m_H0
) >= (2*m_stau_2
)) ? (0.125*IT_0000
      *IT_0001*IT_0935*conj(IT_0935)*IT_0938
) : (0
));
    const ccomplex_t IT_0940 = U_stau_10*(IT_0006*IT_0009*IT_0023*IT_0925*conj
      (U_stau_01) + (-2)*(IT_0006*IT_0011*IT_0022*IT_0023*IT_0646 + -m_W*IT_0003
      *IT_0004*IT_0655)*conj(U_stau_11)) + (-2)*U_stau_00*((IT_0006*IT_0011
      *IT_0022*IT_0023*IT_0646 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0655)
      *conj(U_stau_01) + -1./2*IT_0006*IT_0009*IT_0023*IT_0925*conj(U_stau_11));
    const ccomplex_t IT_0941 = (0 + _Complex_I*-1)*(e_em*IT_0940 + (
      -1.4142135623731)*A_tau*IT_0011*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0942 = (-0.5)*IT_0941;
    const ccomplex_t IT_0943 = cpow(IT_0018 + IT_0929 + -IT_0936, 2);
    const ccomplex_t IT_0944 = cabs(-IT_0929 + 0.25*IT_0001*IT_0943);
    const ccomplex_t IT_0945 = cpow(IT_0944, 0.5);
    const ccomplex_t IT_0946 = ( ((m_H0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0942*conj(IT_0942)*IT_0945
) : (0
));
    const ccomplex_t IT_0947 = conj(U_stau_10)*(IT_0006*IT_0009*IT_0023
      *IT_0925*U_stau_01 + (-2)*(IT_0006*IT_0011*IT_0022*IT_0023*IT_0646 + -m_W
      *IT_0003*IT_0004*IT_0655)*U_stau_11) + (-2)*conj(U_stau_00)*((IT_0006
      *IT_0011*IT_0022*IT_0023*IT_0646 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)
      *IT_0655)*U_stau_01 + -1./2*IT_0006*IT_0009*IT_0023*IT_0925*U_stau_11);
    const ccomplex_t IT_0948 = (0 + _Complex_I*-1)*(e_em*IT_0947 + (
      -1.4142135623731)*A_tau*IT_0011*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0949 = (-0.5)*IT_0948;
    const ccomplex_t IT_0950 = ( ((m_H0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0945*IT_0949*conj(IT_0949)
) : (0
));
    const ccomplex_t IT_0951 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0952 = cabs(0.25*IT_0018 + -IT_0951);
    const ccomplex_t IT_0953 = cpow(IT_0952, 0.5);
    const ccomplex_t IT_0954 = ( ((m_H0
) >= (2*m_snu_tau
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0920*IT_0953
) : (0
));
    const ccomplex_t IT_0955 = (0 + _Complex_I*1)*m_W*e_em*IT_0006*(IT_0008
      *IT_0009 + IT_0011*IT_0012);
    const ccomplex_t IT_0956 = cpow(IT_0955, 2);
    const ccomplex_t IT_0957 = pow(m_W, -4);
    const ccomplex_t IT_0958 = 0.5*IT_0018 + -IT_0130;
    const ccomplex_t IT_0959 = cpow(IT_0958, 2);
    const ccomplex_t IT_0960 = IT_0957*IT_0959;
    const ccomplex_t IT_0961 = cabs(0.25*IT_0018 + -IT_0130);
    const ccomplex_t IT_0962 = cpow(IT_0961, 0.5);
    const ccomplex_t IT_0963 = pow(m_W, 4);
    const ccomplex_t IT_0964 = pow(m_H0, -4);
    const ccomplex_t IT_0965 = clog(IT_0001*IT_0130);
    const ccomplex_t IT_0966 = cpow(-1 + 4*IT_0001*IT_0130, (-0.5));
    const ccomplex_t IT_0967 = pow(m_W, -3);
    const ccomplex_t IT_0968 = pow(m_H0, 3);
    const ccomplex_t IT_0969 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0130)*IT_0967
      *IT_0968);
    const ccomplex_t IT_0970 = IT_0018*IT_0129*(1 + -IT_0001*IT_0130)*(1 + (
      -6.5)*IT_0001*IT_0130 + 23.5*IT_0963*IT_0964) + 1.5*(1 + (-6)*IT_0001
      *IT_0130 + 4*IT_0963*IT_0964)*IT_0965 + (-3)*(1 + (-8)*IT_0001*IT_0130 +
       20*IT_0963*IT_0964)*IT_0966*IT_0969;
    const ccomplex_t IT_0971 = ( ((m_H0
) >= (2*m_W
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0956*(1 + 0.5*IT_0960)*IT_0962
) : (0.25*IT_0000*IT_0001
      *IT_0956*(1 + 0.5*IT_0960)*IT_0962*IT_0970
));
    const ccomplex_t IT_0972 = cpow(IT_0002, -4);
    const ccomplex_t IT_0973 = IT_0004*IT_0011*IT_0012*IT_0972;
    const ccomplex_t IT_0974 = IT_0004*IT_0008*IT_0009*IT_0972;
    const ccomplex_t IT_0975 = cpow(IT_0004, -3);
    const ccomplex_t IT_0976 = IT_0011*IT_0012*IT_0975;
    const ccomplex_t IT_0977 = IT_0008*IT_0009*IT_0975;
    const ccomplex_t IT_0978 = IT_0003*IT_0006*IT_0011*IT_0012;
    const ccomplex_t IT_0979 = IT_0003*IT_0006*IT_0008*IT_0009;
    const ccomplex_t IT_0980 = cpow(IT_0004, -2);
    const ccomplex_t IT_0981 = cpow(IT_0003 + IT_0980, -1);
    const ccomplex_t IT_0982 = (IT_0973 + IT_0974 + IT_0976 + IT_0977 + 2
      *IT_0978 + 2*IT_0979)*IT_0981;
    const ccomplex_t IT_0983 = m_W*e_em;
    const ccomplex_t IT_0984 = IT_0982*IT_0983;
    const ccomplex_t IT_0985 = (0 + _Complex_I*1)*IT_0984;
    const ccomplex_t IT_0986 = cpow(IT_0985, 2);
    const ccomplex_t IT_0987 = pow(m_Z, -4);
    const ccomplex_t IT_0988 = 0.5*IT_0018 + -IT_0109;
    const ccomplex_t IT_0989 = cpow(IT_0988, 2);
    const ccomplex_t IT_0990 = IT_0987*IT_0989;
    const ccomplex_t IT_0991 = cabs(0.25*IT_0018 + -IT_0109);
    const ccomplex_t IT_0992 = cpow(IT_0991, 0.5);
    const ccomplex_t IT_0993 = pow(m_Z, 4);
    const ccomplex_t IT_0994 = clog(IT_0001*IT_0109);
    const ccomplex_t IT_0995 = cpow(-1 + 4*IT_0001*IT_0109, (-0.5));
    const ccomplex_t IT_0996 = pow(m_Z, -3);
    const ccomplex_t IT_0997 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0109)*IT_0968
      *IT_0996);
    const ccomplex_t IT_0998 = IT_0018*IT_0108*(1 + -IT_0001*IT_0109)*(1 + (
      -6.5)*IT_0001*IT_0109 + 23.5*IT_0964*IT_0993) + 1.5*(1 + (-6)*IT_0001
      *IT_0109 + 4*IT_0964*IT_0993)*IT_0994 + (-3)*(1 + (-8)*IT_0001*IT_0109 +
       20*IT_0964*IT_0993)*IT_0995*IT_0997;
    const ccomplex_t IT_0999 = ( ((m_H0
) >= (2*m_Z
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0986*(1 + 0.5*IT_0990)*IT_0992
) : (0.125*IT_0000*IT_0001
      *IT_0986*(1 + 0.5*IT_0990)*IT_0992*IT_0998
));
    const ccomplex_t IT_1000 = IT_0021 + IT_0031 + IT_0040 + IT_0048 + IT_0056
       + IT_0063 + IT_0071 + IT_0079 + IT_0087 + IT_0097 + IT_0104 + IT_0116 +
       IT_0126 + IT_0137 + IT_0140 + IT_0148 + IT_0172 + IT_0196 + IT_0221 +
       IT_0240 + IT_0273 + IT_0296 + IT_0345 + IT_0374 + IT_0423 + IT_0472 +
       IT_0495 + IT_0544 + IT_0593 + IT_0642 + IT_0650 + IT_0662 + IT_0669 +
       IT_0676 + IT_0680 + IT_0687 + IT_0694 + 2*IT_0701 + IT_0708 + IT_0714 +
       IT_0717 + IT_0723 + IT_0726 + IT_0733 + IT_0739 + IT_0742 + IT_0748 +
       IT_0751 + 2*IT_0758 + IT_0765 + IT_0771 + 2*IT_0778 + IT_0785 + IT_0791 +
       IT_0794 + IT_0800 + IT_0803 + IT_0809 + IT_0812 + IT_0819 + IT_0825 +
       IT_0828 + IT_0834 + IT_0837 + IT_0841 + IT_0842 + 2*IT_0849 + IT_0858 +
       IT_0865 + IT_0872 + IT_0876 + IT_0883 + IT_0890 + 2*IT_0897 + IT_0904 +
       IT_0910 + 2*IT_0917 + 2*IT_0924 + IT_0932 + IT_0939 + IT_0946 + IT_0950 +
       IT_0954 + IT_0971 + IT_0999;
    return create_ccomplex_return(IT_1000);
}

