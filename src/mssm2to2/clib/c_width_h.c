#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_h(
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
    const ccomplex_t IT_0001 = pow(m_h, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = IT_0005 + IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = cpow(IT_0009, 3);
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = cpow(IT_0009, 2);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0010 + -IT_0011*(IT_0011*(IT_0008*IT_0009 + IT_0011*IT_0012) + 
      -IT_0012*IT_0013));
    const ccomplex_t IT_0015 = 1.5*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0015, 2);
    const ccomplex_t IT_0017 = pow(m_h, 2);
    const ccomplex_t IT_0018 = cabs((-0.75)*IT_0017);
    const ccomplex_t IT_0019 = cpow(IT_0018, 0.5);
    const ccomplex_t IT_0020 = ( ((m_h
) >= (2*m_h
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0016*IT_0019
) : (0
));
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = cpow(IT_0008, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_b*e_em*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = pow(m_b, 2);
    const ccomplex_t IT_0027 = 0.5*IT_0017 + -IT_0026;
    const ccomplex_t IT_0028 = cabs(0.25*IT_0017 + -IT_0026);
    const ccomplex_t IT_0029 = cpow(IT_0028, 0.5);
    const ccomplex_t IT_0030 = ( ((m_h
) >= (2*m_b
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0025*((-12)*IT_0026 + 12*IT_0027)*IT_0029
) : (0
));
    const ccomplex_t IT_0031 = cpow(IT_0012, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_c*e_em*IT_0006*IT_0011
      *IT_0021*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0033, 2);
    const ccomplex_t IT_0035 = pow(m_c, 2);
    const ccomplex_t IT_0036 = 0.5*IT_0017 + -IT_0035;
    const ccomplex_t IT_0037 = cabs(0.25*IT_0017 + -IT_0035);
    const ccomplex_t IT_0038 = cpow(IT_0037, 0.5);
    const ccomplex_t IT_0039 = ( ((m_h
) >= (2*m_c
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0034*((-12)*IT_0035 + 12*IT_0036)*IT_0038
) : (0
));
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_d*e_em*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0041, 2);
    const ccomplex_t IT_0043 = pow(m_d, 2);
    const ccomplex_t IT_0044 = 0.5*IT_0017 + -IT_0043;
    const ccomplex_t IT_0045 = cabs(0.25*IT_0017 + -IT_0043);
    const ccomplex_t IT_0046 = cpow(IT_0045, 0.5);
    const ccomplex_t IT_0047 = ( ((m_h
) >= (2*m_d
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0042*((-12)*IT_0043 + 12*IT_0044)*IT_0046
) : (0
));
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_e*e_em*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0049, 2);
    const ccomplex_t IT_0051 = pow(m_e, 2);
    const ccomplex_t IT_0052 = 0.5*IT_0017 + -IT_0051;
    const ccomplex_t IT_0053 = cabs(0.25*IT_0017 + -IT_0051);
    const ccomplex_t IT_0054 = cpow(IT_0053, 0.5);
    const ccomplex_t IT_0055 = ( ((m_h
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *IT_0050*(IT_0051 + -IT_0052)*IT_0054
) : (0
));
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_s*e_em*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = cpow(IT_0057, 2);
    const ccomplex_t IT_0059 = pow(m_s, 2);
    const ccomplex_t IT_0060 = 0.5*IT_0017 + -IT_0059;
    const ccomplex_t IT_0061 = cabs(0.25*IT_0017 + -IT_0059);
    const ccomplex_t IT_0062 = cpow(IT_0061, 0.5);
    const ccomplex_t IT_0063 = ( ((m_h
) >= (2*m_s
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0058*((-12)*IT_0059 + 12*IT_0060)*IT_0062
) : (0
));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0011
      *IT_0021*IT_0031;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0065, 2);
    const ccomplex_t IT_0067 = pow(m_t, 2);
    const ccomplex_t IT_0068 = 0.5*IT_0017 + -IT_0067;
    const ccomplex_t IT_0069 = cabs(0.25*IT_0017 + -IT_0067);
    const ccomplex_t IT_0070 = cpow(IT_0069, 0.5);
    const ccomplex_t IT_0071 = ( ((m_h
) >= (2*m_t
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0066*((-12)*IT_0067 + 12*IT_0068)*IT_0070
) : (0
));
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_u*e_em*IT_0006*IT_0011
      *IT_0021*IT_0031;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = cpow(IT_0073, 2);
    const ccomplex_t IT_0075 = pow(m_u, 2);
    const ccomplex_t IT_0076 = 0.5*IT_0017 + -IT_0075;
    const ccomplex_t IT_0077 = cabs(0.25*IT_0017 + -IT_0075);
    const ccomplex_t IT_0078 = cpow(IT_0077, 0.5);
    const ccomplex_t IT_0079 = ( ((m_h
) >= (2*m_u
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0074*((-12)*IT_0075 + 12*IT_0076)*IT_0078
) : (0
));
    const ccomplex_t IT_0080 = cpow(IT_0012, 3);
    const ccomplex_t IT_0081 = cpow(IT_0012, 2);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0011
      *IT_0080 + -IT_0008*(IT_0008*(IT_0008*IT_0009 + IT_0011*IT_0012) + 
      -IT_0009*IT_0081));
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0083, 2);
    const ccomplex_t IT_0085 = pow(m_A0, 2);
    const ccomplex_t IT_0086 = cabs(0.25*IT_0017 + -IT_0085);
    const ccomplex_t IT_0087 = cpow(IT_0086, 0.5);
    const ccomplex_t IT_0088 = ( ((m_h
) >= (2*m_A0
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0084*IT_0087
) : (0
));
    const ccomplex_t IT_0089 = cpow(IT_0002, -1);
    const ccomplex_t IT_0090 = e_em*IT_0006*(IT_0008*IT_0011 + IT_0009*IT_0012
      )*IT_0089;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0091, 2);
    const ccomplex_t IT_0093 = pow(m_Z, -2);
    const ccomplex_t IT_0094 = pow(m_Z, 2);
    const ccomplex_t IT_0095 = 0.5*IT_0017 + (-0.5)*IT_0085 + 0.5*IT_0094;
    const ccomplex_t IT_0096 = cpow(IT_0095, 2);
    const ccomplex_t IT_0097 = IT_0093*IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0017 + -IT_0085 + IT_0094, 2);
    const ccomplex_t IT_0099 = cabs(-IT_0094 + 0.25*IT_0001*IT_0098);
    const ccomplex_t IT_0100 = cpow(IT_0099, 0.5);
    const ccomplex_t IT_0101 = ( ((m_h
) >= (m_Z + m_A0
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0092*(IT_0017 + -IT_0097)*IT_0100
) : (0
));
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0010 + 5*IT_0011*((-1)*IT_0011*(IT_0008*IT_0009 + 1./5*IT_0011*IT_0012
      ) + IT_0012*IT_0013));
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0103, 2);
    const ccomplex_t IT_0105 = pow(m_H0, 2);
    const ccomplex_t IT_0106 = cabs(0.25*IT_0017 + -IT_0105);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_h
) >= (2*m_H0
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0104*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = cpow(IT_0011, 2);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0010
      *IT_0012 + (-5)*IT_0011*(IT_0009*IT_0011*IT_0012 + 1*IT_0008*(IT_0013 + 
      -1./5*IT_0109)));
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = cpow(IT_0111, 2);
    const ccomplex_t IT_0113 = cpow(IT_0017 + (-0.5)*IT_0105, 2);
    const ccomplex_t IT_0114 = cabs(-IT_0017 + IT_0001*IT_0113);
    const ccomplex_t IT_0115 = cpow(IT_0114, 0.5);
    const ccomplex_t IT_0116 = ( ((m_h
) >= (m_h + m_H0
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0112*IT_0115
) : (0
));
    const ccomplex_t IT_0117 = IT_0005 + -IT_0006;
    const ccomplex_t IT_0118 = IT_0005 + 3*IT_0006;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_W*e_em*(IT_0011*IT_0080
      *IT_0117 + IT_0008*(IT_0009*IT_0081*IT_0118 + -IT_0008*(IT_0008*IT_0009
      *IT_0117 + IT_0011*IT_0012*IT_0118)));
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = pow(m_Hp, 2);
    const ccomplex_t IT_0123 = cabs(0.25*IT_0017 + -IT_0122);
    const ccomplex_t IT_0124 = cpow(IT_0123, 0.5);
    const ccomplex_t IT_0125 = ( ((m_h
) >= (2*m_Hp
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0121*IT_0124
) : (0
));
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*e_em*IT_0006*(IT_0008
      *IT_0011 + IT_0009*IT_0012);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = cpow(IT_0127, 2);
    const ccomplex_t IT_0129 = pow(m_W, -2);
    const ccomplex_t IT_0130 = pow(m_W, 2);
    const ccomplex_t IT_0131 = 0.5*IT_0017 + (-0.5)*IT_0122 + 0.5*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0131, 2);
    const ccomplex_t IT_0133 = IT_0129*IT_0132;
    const ccomplex_t IT_0134 = cpow(IT_0017 + -IT_0122 + IT_0130, 2);
    const ccomplex_t IT_0135 = cabs(-IT_0130 + 0.25*IT_0001*IT_0134);
    const ccomplex_t IT_0136 = cpow(IT_0135, 0.5);
    const ccomplex_t IT_0137 = ( ((m_h
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0128*(IT_0017 + -IT_0133)*IT_0136
) : (0
));
    const ccomplex_t IT_0138 = -IT_0127;
    const ccomplex_t IT_0139 = cpow(IT_0138, 2);
    const ccomplex_t IT_0140 = ( ((m_h
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*(IT_0017 + -IT_0133)*IT_0136*IT_0139
) : (0
));
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*e_em*m_mu*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = pow(m_mu, 2);
    const ccomplex_t IT_0145 = 0.5*IT_0017 + -IT_0144;
    const ccomplex_t IT_0146 = cabs(0.25*IT_0017 + -IT_0144);
    const ccomplex_t IT_0147 = cpow(IT_0146, 0.5);
    const ccomplex_t IT_0148 = ( ((m_h
) >= (2*m_mu
)) ? (0.5*IT_0000*IT_0001
      *IT_0143*(IT_0144 + -IT_0145)*IT_0147
) : (0
));
    const ccomplex_t IT_0149 = V_u1*U_Wm1;
    const ccomplex_t IT_0150 = IT_0011*IT_0149;
    const ccomplex_t IT_0151 = U_d1*V_Wp1;
    const ccomplex_t IT_0152 = IT_0009*IT_0151;
    const ccomplex_t IT_0153 = IT_0006*(IT_0150 + -IT_0152);
    const ccomplex_t IT_0154 = 1.4142135623731*e_em*IT_0153;
    const ccomplex_t IT_0155 = 0.5*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0159 = IT_0011*IT_0158;
    const ccomplex_t IT_0160 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0161 = IT_0009*IT_0160;
    const ccomplex_t IT_0162 = IT_0006*(IT_0159 + -IT_0161);
    const ccomplex_t IT_0163 = 1.4142135623731*e_em*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = pow(m_C_1, 2);
    const ccomplex_t IT_0168 = 0.5*IT_0017 + -IT_0167;
    const ccomplex_t IT_0169 = 2*IT_0168;
    const ccomplex_t IT_0170 = cabs(0.25*IT_0017 + -IT_0167);
    const ccomplex_t IT_0171 = cpow(IT_0170, 0.5);
    const ccomplex_t IT_0172 = ( ((m_h
) >= (2*m_C_1
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0157)*(IT_0166*IT_0167 + (-0.5)*IT_0157*IT_0169) + conj
      (IT_0166)*(IT_0157*IT_0167 + (-0.5)*IT_0166*IT_0169))*IT_0171
) : (0
));
    const ccomplex_t IT_0173 = V_u2*U_Wm2;
    const ccomplex_t IT_0174 = IT_0011*IT_0173;
    const ccomplex_t IT_0175 = U_d2*V_Wp2;
    const ccomplex_t IT_0176 = IT_0009*IT_0175;
    const ccomplex_t IT_0177 = IT_0006*(IT_0174 + -IT_0176);
    const ccomplex_t IT_0178 = 1.4142135623731*e_em*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0183 = IT_0011*IT_0182;
    const ccomplex_t IT_0184 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0185 = IT_0009*IT_0184;
    const ccomplex_t IT_0186 = IT_0006*(IT_0183 + -IT_0185);
    const ccomplex_t IT_0187 = 1.4142135623731*e_em*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = pow(m_C_2, 2);
    const ccomplex_t IT_0192 = 0.5*IT_0017 + -IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0192;
    const ccomplex_t IT_0194 = cabs(0.25*IT_0017 + -IT_0191);
    const ccomplex_t IT_0195 = cpow(IT_0194, 0.5);
    const ccomplex_t IT_0196 = ( ((m_h
) >= (2*m_C_2
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0181)*(IT_0190*IT_0191 + (-0.5)*IT_0181*IT_0193) + conj
      (IT_0190)*(IT_0181*IT_0191 + (-0.5)*IT_0190*IT_0193))*IT_0195
) : (0
));
    const ccomplex_t IT_0197 = V_u1*U_Wm2;
    const ccomplex_t IT_0198 = IT_0011*IT_0197;
    const ccomplex_t IT_0199 = U_d2*V_Wp1;
    const ccomplex_t IT_0200 = IT_0009*IT_0199;
    const ccomplex_t IT_0201 = IT_0006*(IT_0198 + -IT_0200);
    const ccomplex_t IT_0202 = 1.4142135623731*e_em*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0207 = IT_0011*IT_0206;
    const ccomplex_t IT_0208 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0209 = IT_0009*IT_0208;
    const ccomplex_t IT_0210 = IT_0006*(IT_0207 + -IT_0209);
    const ccomplex_t IT_0211 = 1.4142135623731*e_em*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = m_C_1*m_C_2;
    const ccomplex_t IT_0216 = cpow(IT_0017 + IT_0167 + -IT_0191, 2);
    const ccomplex_t IT_0217 = -IT_0167 + 0.25*IT_0001*(IT_0017 + -IT_0167 +
       IT_0191)*(IT_0017 + IT_0167 + -IT_0191) + 0.25*IT_0001*IT_0216;
    const ccomplex_t IT_0218 = 2*IT_0217;
    const ccomplex_t IT_0219 = cabs(-IT_0167 + 0.25*IT_0001*IT_0216);
    const ccomplex_t IT_0220 = cpow(IT_0219, 0.5);
    const ccomplex_t IT_0221 = ( ((m_h
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0205)*(IT_0214*IT_0215 + (-0.5)*IT_0205*IT_0218)
       + conj(IT_0214)*(IT_0205*IT_0215 + (-0.5)*IT_0214*IT_0218))*IT_0220
) : 
      (0
));
    const ccomplex_t IT_0222 = V_u2*U_Wm1;
    const ccomplex_t IT_0223 = IT_0011*IT_0222;
    const ccomplex_t IT_0224 = U_d1*V_Wp2;
    const ccomplex_t IT_0225 = IT_0009*IT_0224;
    const ccomplex_t IT_0226 = IT_0006*(IT_0223 + -IT_0225);
    const ccomplex_t IT_0227 = 1.4142135623731*e_em*IT_0226;
    const ccomplex_t IT_0228 = 0.5*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0232 = IT_0011*IT_0231;
    const ccomplex_t IT_0233 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0234 = IT_0009*IT_0233;
    const ccomplex_t IT_0235 = IT_0006*(IT_0232 + -IT_0234);
    const ccomplex_t IT_0236 = 1.4142135623731*e_em*IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = ( ((m_h
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*IT_0220*(conj(IT_0230)*((-0.5)*IT_0218*IT_0230 + IT_0215
      *IT_0239) + (IT_0215*IT_0230 + (-0.5)*IT_0218*IT_0239)*conj(IT_0239))
) : 
      (0
));
    const ccomplex_t IT_0241 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0242 = IT_0011*IT_0089;
    const ccomplex_t IT_0243 = IT_0241*IT_0242;
    const ccomplex_t IT_0244 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0245 = IT_0009*IT_0089;
    const ccomplex_t IT_0246 = IT_0244*IT_0245;
    const ccomplex_t IT_0247 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0248 = IT_0006*IT_0011;
    const ccomplex_t IT_0249 = IT_0247*IT_0248;
    const ccomplex_t IT_0250 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0251 = IT_0006*IT_0009;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0243 + IT_0246 + 
      -IT_0249 + -IT_0252);
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0256 = IT_0242*IT_0255;
    const ccomplex_t IT_0257 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0258 = IT_0245*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0260 = IT_0248*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0262 = IT_0251*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*(IT_0256 + IT_0258 + 
      -IT_0260 + -IT_0262);
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = pow(m_N_1, 2);
    const ccomplex_t IT_0266 = 0.5*IT_0017 + -IT_0265;
    const ccomplex_t IT_0267 = cabs(0.25*IT_0017 + -IT_0265);
    const ccomplex_t IT_0268 = cpow(IT_0267, 0.5);
    const ccomplex_t IT_0269 = ( ((m_h
) >= (2*m_N_1
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0254)*IT_0264 + IT_0254*conj(IT_0264))*IT_0265 + -
      (IT_0254*conj(IT_0254) + IT_0264*conj(IT_0264))*IT_0266)*IT_0268
) : (0
));
    const ccomplex_t IT_0270 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0271 = IT_0242*IT_0270;
    const ccomplex_t IT_0272 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0273 = IT_0245*IT_0272;
    const ccomplex_t IT_0274 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0275 = IT_0248*IT_0274;
    const ccomplex_t IT_0276 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0277 = IT_0251*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*(IT_0271 + IT_0273 + 
      -IT_0275 + -IT_0277);
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0281 = IT_0242*IT_0280;
    const ccomplex_t IT_0282 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0283 = IT_0245*IT_0282;
    const ccomplex_t IT_0284 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0285 = IT_0248*IT_0284;
    const ccomplex_t IT_0286 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0287 = IT_0251*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*(IT_0281 + IT_0283 + 
      -IT_0285 + -IT_0287);
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = pow(m_N_2, 2);
    const ccomplex_t IT_0291 = 0.5*IT_0017 + -IT_0290;
    const ccomplex_t IT_0292 = cabs(0.25*IT_0017 + -IT_0290);
    const ccomplex_t IT_0293 = cpow(IT_0292, 0.5);
    const ccomplex_t IT_0294 = ( ((m_h
) >= (2*m_N_2
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0279)*IT_0289 + IT_0279*conj(IT_0289))*IT_0290 + -
      (IT_0279*conj(IT_0279) + IT_0289*conj(IT_0289))*IT_0291)*IT_0293
) : (0
));
    const ccomplex_t IT_0295 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0296 = IT_0242*IT_0295;
    const ccomplex_t IT_0297 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0298 = IT_0245*IT_0297;
    const ccomplex_t IT_0299 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0300 = IT_0242*IT_0299;
    const ccomplex_t IT_0301 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0302 = IT_0245*IT_0301;
    const ccomplex_t IT_0303 = IT_0296 + IT_0298 + IT_0300 + IT_0302;
    const ccomplex_t IT_0304 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0305 = IT_0248*IT_0304;
    const ccomplex_t IT_0306 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0307 = IT_0251*IT_0306;
    const ccomplex_t IT_0308 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0309 = IT_0248*IT_0308;
    const ccomplex_t IT_0310 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0311 = IT_0251*IT_0310;
    const ccomplex_t IT_0312 = -IT_0305 + -IT_0307 + -IT_0309 + -IT_0311;
    const ccomplex_t IT_0313 = IT_0303 + IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = (-0.5)*IT_0314;
    const ccomplex_t IT_0316 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0317 = IT_0242*IT_0316;
    const ccomplex_t IT_0318 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0319 = IT_0245*IT_0318;
    const ccomplex_t IT_0320 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0321 = IT_0242*IT_0320;
    const ccomplex_t IT_0322 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0323 = IT_0245*IT_0322;
    const ccomplex_t IT_0324 = IT_0317 + IT_0319 + IT_0321 + IT_0323;
    const ccomplex_t IT_0325 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0326 = IT_0248*IT_0325;
    const ccomplex_t IT_0327 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0328 = IT_0251*IT_0327;
    const ccomplex_t IT_0329 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0330 = IT_0248*IT_0329;
    const ccomplex_t IT_0331 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0332 = IT_0251*IT_0331;
    const ccomplex_t IT_0333 = -IT_0326 + -IT_0328 + -IT_0330 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0324 + IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = (-0.5)*IT_0335;
    const ccomplex_t IT_0337 = m_N_1*m_N_2;
    const ccomplex_t IT_0338 = cpow(IT_0017 + IT_0265 + -IT_0290, 2);
    const ccomplex_t IT_0339 = -IT_0265 + 0.25*IT_0001*(IT_0017 + -IT_0265 +
       IT_0290)*(IT_0017 + IT_0265 + -IT_0290) + 0.25*IT_0001*IT_0338;
    const ccomplex_t IT_0340 = 2*IT_0339;
    const ccomplex_t IT_0341 = cabs(-IT_0265 + 0.25*IT_0001*IT_0338);
    const ccomplex_t IT_0342 = cpow(IT_0341, 0.5);
    const ccomplex_t IT_0343 = ( ((m_h
) >= (m_N_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0315)*(IT_0336*IT_0337 + (-0.5)*IT_0315*IT_0340)
       + conj(IT_0336)*(IT_0315*IT_0337 + (-0.5)*IT_0336*IT_0340))*IT_0342
) : 
      (0
));
    const ccomplex_t IT_0344 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0345 = IT_0242*IT_0344;
    const ccomplex_t IT_0346 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0347 = IT_0245*IT_0346;
    const ccomplex_t IT_0348 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0349 = IT_0248*IT_0348;
    const ccomplex_t IT_0350 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0351 = IT_0251*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*(IT_0345 + IT_0347 + 
      -IT_0349 + -IT_0351);
    const ccomplex_t IT_0353 = -IT_0352;
    const ccomplex_t IT_0354 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0355 = IT_0242*IT_0354;
    const ccomplex_t IT_0356 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0357 = IT_0245*IT_0356;
    const ccomplex_t IT_0358 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0359 = IT_0248*IT_0358;
    const ccomplex_t IT_0360 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0361 = IT_0251*IT_0360;
    const ccomplex_t IT_0362 = (0 + _Complex_I*1)*(IT_0355 + IT_0357 + 
      -IT_0359 + -IT_0361);
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = pow(m_N_3, 2);
    const ccomplex_t IT_0365 = 0.5*IT_0017 + -IT_0364;
    const ccomplex_t IT_0366 = cabs(0.25*IT_0017 + -IT_0364);
    const ccomplex_t IT_0367 = cpow(IT_0366, 0.5);
    const ccomplex_t IT_0368 = ( ((m_h
) >= (2*m_N_3
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0353)*IT_0363 + IT_0353*conj(IT_0363))*IT_0364 + -
      (IT_0353*conj(IT_0353) + IT_0363*conj(IT_0363))*IT_0365)*IT_0367
) : (0
));
    const ccomplex_t IT_0369 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0370 = IT_0242*IT_0369;
    const ccomplex_t IT_0371 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0372 = IT_0245*IT_0371;
    const ccomplex_t IT_0373 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0374 = IT_0242*IT_0373;
    const ccomplex_t IT_0375 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0376 = IT_0245*IT_0375;
    const ccomplex_t IT_0377 = IT_0370 + IT_0372 + IT_0374 + IT_0376;
    const ccomplex_t IT_0378 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0379 = IT_0248*IT_0378;
    const ccomplex_t IT_0380 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0381 = IT_0251*IT_0380;
    const ccomplex_t IT_0382 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0383 = IT_0248*IT_0382;
    const ccomplex_t IT_0384 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0385 = IT_0251*IT_0384;
    const ccomplex_t IT_0386 = -IT_0379 + -IT_0381 + -IT_0383 + -IT_0385;
    const ccomplex_t IT_0387 = IT_0377 + IT_0386;
    const ccomplex_t IT_0388 = (0 + _Complex_I*1)*IT_0387;
    const ccomplex_t IT_0389 = (-0.5)*IT_0388;
    const ccomplex_t IT_0390 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0391 = IT_0242*IT_0390;
    const ccomplex_t IT_0392 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0393 = IT_0245*IT_0392;
    const ccomplex_t IT_0394 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0395 = IT_0242*IT_0394;
    const ccomplex_t IT_0396 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0397 = IT_0245*IT_0396;
    const ccomplex_t IT_0398 = IT_0391 + IT_0393 + IT_0395 + IT_0397;
    const ccomplex_t IT_0399 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0400 = IT_0248*IT_0399;
    const ccomplex_t IT_0401 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0402 = IT_0251*IT_0401;
    const ccomplex_t IT_0403 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0404 = IT_0248*IT_0403;
    const ccomplex_t IT_0405 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0406 = IT_0251*IT_0405;
    const ccomplex_t IT_0407 = -IT_0400 + -IT_0402 + -IT_0404 + -IT_0406;
    const ccomplex_t IT_0408 = IT_0398 + IT_0407;
    const ccomplex_t IT_0409 = (0 + _Complex_I*1)*IT_0408;
    const ccomplex_t IT_0410 = (-0.5)*IT_0409;
    const ccomplex_t IT_0411 = m_N_1*m_N_3;
    const ccomplex_t IT_0412 = cpow(IT_0017 + IT_0265 + -IT_0364, 2);
    const ccomplex_t IT_0413 = -IT_0265 + 0.25*IT_0001*(IT_0017 + -IT_0265 +
       IT_0364)*(IT_0017 + IT_0265 + -IT_0364) + 0.25*IT_0001*IT_0412;
    const ccomplex_t IT_0414 = 2*IT_0413;
    const ccomplex_t IT_0415 = cabs(-IT_0265 + 0.25*IT_0001*IT_0412);
    const ccomplex_t IT_0416 = cpow(IT_0415, 0.5);
    const ccomplex_t IT_0417 = ( ((m_h
) >= (m_N_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0389)*(IT_0410*IT_0411 + (-0.5)*IT_0389*IT_0414)
       + conj(IT_0410)*(IT_0389*IT_0411 + (-0.5)*IT_0410*IT_0414))*IT_0416
) : 
      (0
));
    const ccomplex_t IT_0418 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0419 = IT_0242*IT_0418;
    const ccomplex_t IT_0420 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0421 = IT_0245*IT_0420;
    const ccomplex_t IT_0422 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0423 = IT_0242*IT_0422;
    const ccomplex_t IT_0424 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0425 = IT_0245*IT_0424;
    const ccomplex_t IT_0426 = IT_0419 + IT_0421 + IT_0423 + IT_0425;
    const ccomplex_t IT_0427 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0428 = IT_0248*IT_0427;
    const ccomplex_t IT_0429 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0430 = IT_0251*IT_0429;
    const ccomplex_t IT_0431 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0432 = IT_0248*IT_0431;
    const ccomplex_t IT_0433 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0434 = IT_0251*IT_0433;
    const ccomplex_t IT_0435 = -IT_0428 + -IT_0430 + -IT_0432 + -IT_0434;
    const ccomplex_t IT_0436 = IT_0426 + IT_0435;
    const ccomplex_t IT_0437 = (0 + _Complex_I*1)*IT_0436;
    const ccomplex_t IT_0438 = (-0.5)*IT_0437;
    const ccomplex_t IT_0439 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0440 = IT_0242*IT_0439;
    const ccomplex_t IT_0441 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0442 = IT_0245*IT_0441;
    const ccomplex_t IT_0443 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0444 = IT_0242*IT_0443;
    const ccomplex_t IT_0445 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0446 = IT_0245*IT_0445;
    const ccomplex_t IT_0447 = IT_0440 + IT_0442 + IT_0444 + IT_0446;
    const ccomplex_t IT_0448 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0449 = IT_0248*IT_0448;
    const ccomplex_t IT_0450 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0451 = IT_0251*IT_0450;
    const ccomplex_t IT_0452 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0453 = IT_0248*IT_0452;
    const ccomplex_t IT_0454 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0455 = IT_0251*IT_0454;
    const ccomplex_t IT_0456 = -IT_0449 + -IT_0451 + -IT_0453 + -IT_0455;
    const ccomplex_t IT_0457 = IT_0447 + IT_0456;
    const ccomplex_t IT_0458 = (0 + _Complex_I*1)*IT_0457;
    const ccomplex_t IT_0459 = (-0.5)*IT_0458;
    const ccomplex_t IT_0460 = m_N_2*m_N_3;
    const ccomplex_t IT_0461 = cpow(IT_0017 + IT_0290 + -IT_0364, 2);
    const ccomplex_t IT_0462 = -IT_0290 + 0.25*IT_0001*(IT_0017 + -IT_0290 +
       IT_0364)*(IT_0017 + IT_0290 + -IT_0364) + 0.25*IT_0001*IT_0461;
    const ccomplex_t IT_0463 = 2*IT_0462;
    const ccomplex_t IT_0464 = cabs(-IT_0290 + 0.25*IT_0001*IT_0461);
    const ccomplex_t IT_0465 = cpow(IT_0464, 0.5);
    const ccomplex_t IT_0466 = ( ((m_h
) >= (m_N_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0438)*(IT_0459*IT_0460 + (-0.5)*IT_0438*IT_0463)
       + conj(IT_0459)*(IT_0438*IT_0460 + (-0.5)*IT_0459*IT_0463))*IT_0465
) : 
      (0
));
    const ccomplex_t IT_0467 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0468 = IT_0242*IT_0467;
    const ccomplex_t IT_0469 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0470 = IT_0245*IT_0469;
    const ccomplex_t IT_0471 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0472 = IT_0248*IT_0471;
    const ccomplex_t IT_0473 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0474 = IT_0251*IT_0473;
    const ccomplex_t IT_0475 = (0 + _Complex_I*1)*(IT_0468 + IT_0470 + 
      -IT_0472 + -IT_0474);
    const ccomplex_t IT_0476 = -IT_0475;
    const ccomplex_t IT_0477 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0478 = IT_0242*IT_0477;
    const ccomplex_t IT_0479 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0480 = IT_0245*IT_0479;
    const ccomplex_t IT_0481 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0482 = IT_0248*IT_0481;
    const ccomplex_t IT_0483 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0484 = IT_0251*IT_0483;
    const ccomplex_t IT_0485 = (0 + _Complex_I*1)*(IT_0478 + IT_0480 + 
      -IT_0482 + -IT_0484);
    const ccomplex_t IT_0486 = -IT_0485;
    const ccomplex_t IT_0487 = pow(m_N_4, 2);
    const ccomplex_t IT_0488 = 0.5*IT_0017 + -IT_0487;
    const ccomplex_t IT_0489 = cabs(0.25*IT_0017 + -IT_0487);
    const ccomplex_t IT_0490 = cpow(IT_0489, 0.5);
    const ccomplex_t IT_0491 = ( ((m_h
) >= (2*m_N_4
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0476)*IT_0486 + IT_0476*conj(IT_0486))*IT_0487 + -
      (IT_0476*conj(IT_0476) + IT_0486*conj(IT_0486))*IT_0488)*IT_0490
) : (0
));
    const ccomplex_t IT_0492 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0493 = IT_0242*IT_0492;
    const ccomplex_t IT_0494 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0495 = IT_0245*IT_0494;
    const ccomplex_t IT_0496 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0497 = IT_0242*IT_0496;
    const ccomplex_t IT_0498 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0499 = IT_0245*IT_0498;
    const ccomplex_t IT_0500 = IT_0493 + IT_0495 + IT_0497 + IT_0499;
    const ccomplex_t IT_0501 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0502 = IT_0248*IT_0501;
    const ccomplex_t IT_0503 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0504 = IT_0251*IT_0503;
    const ccomplex_t IT_0505 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0506 = IT_0248*IT_0505;
    const ccomplex_t IT_0507 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0508 = IT_0251*IT_0507;
    const ccomplex_t IT_0509 = -IT_0502 + -IT_0504 + -IT_0506 + -IT_0508;
    const ccomplex_t IT_0510 = IT_0500 + IT_0509;
    const ccomplex_t IT_0511 = (0 + _Complex_I*1)*IT_0510;
    const ccomplex_t IT_0512 = (-0.5)*IT_0511;
    const ccomplex_t IT_0513 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0514 = IT_0242*IT_0513;
    const ccomplex_t IT_0515 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0516 = IT_0245*IT_0515;
    const ccomplex_t IT_0517 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0518 = IT_0242*IT_0517;
    const ccomplex_t IT_0519 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0520 = IT_0245*IT_0519;
    const ccomplex_t IT_0521 = IT_0514 + IT_0516 + IT_0518 + IT_0520;
    const ccomplex_t IT_0522 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0523 = IT_0248*IT_0522;
    const ccomplex_t IT_0524 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0525 = IT_0251*IT_0524;
    const ccomplex_t IT_0526 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0527 = IT_0248*IT_0526;
    const ccomplex_t IT_0528 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0529 = IT_0251*IT_0528;
    const ccomplex_t IT_0530 = -IT_0523 + -IT_0525 + -IT_0527 + -IT_0529;
    const ccomplex_t IT_0531 = IT_0521 + IT_0530;
    const ccomplex_t IT_0532 = (0 + _Complex_I*1)*IT_0531;
    const ccomplex_t IT_0533 = (-0.5)*IT_0532;
    const ccomplex_t IT_0534 = m_N_1*m_N_4;
    const ccomplex_t IT_0535 = cpow(IT_0017 + IT_0265 + -IT_0487, 2);
    const ccomplex_t IT_0536 = -IT_0265 + 0.25*IT_0001*(IT_0017 + -IT_0265 +
       IT_0487)*(IT_0017 + IT_0265 + -IT_0487) + 0.25*IT_0001*IT_0535;
    const ccomplex_t IT_0537 = 2*IT_0536;
    const ccomplex_t IT_0538 = cabs(-IT_0265 + 0.25*IT_0001*IT_0535);
    const ccomplex_t IT_0539 = cpow(IT_0538, 0.5);
    const ccomplex_t IT_0540 = ( ((m_h
) >= (m_N_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0512)*(IT_0533*IT_0534 + (-0.5)*IT_0512*IT_0537)
       + conj(IT_0533)*(IT_0512*IT_0534 + (-0.5)*IT_0533*IT_0537))*IT_0539
) : 
      (0
));
    const ccomplex_t IT_0541 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0542 = IT_0242*IT_0541;
    const ccomplex_t IT_0543 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0544 = IT_0245*IT_0543;
    const ccomplex_t IT_0545 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0546 = IT_0242*IT_0545;
    const ccomplex_t IT_0547 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0548 = IT_0245*IT_0547;
    const ccomplex_t IT_0549 = IT_0542 + IT_0544 + IT_0546 + IT_0548;
    const ccomplex_t IT_0550 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0551 = IT_0248*IT_0550;
    const ccomplex_t IT_0552 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0553 = IT_0251*IT_0552;
    const ccomplex_t IT_0554 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0555 = IT_0248*IT_0554;
    const ccomplex_t IT_0556 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0557 = IT_0251*IT_0556;
    const ccomplex_t IT_0558 = -IT_0551 + -IT_0553 + -IT_0555 + -IT_0557;
    const ccomplex_t IT_0559 = IT_0549 + IT_0558;
    const ccomplex_t IT_0560 = (0 + _Complex_I*1)*IT_0559;
    const ccomplex_t IT_0561 = (-0.5)*IT_0560;
    const ccomplex_t IT_0562 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0563 = IT_0242*IT_0562;
    const ccomplex_t IT_0564 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0565 = IT_0245*IT_0564;
    const ccomplex_t IT_0566 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0567 = IT_0242*IT_0566;
    const ccomplex_t IT_0568 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0569 = IT_0245*IT_0568;
    const ccomplex_t IT_0570 = IT_0563 + IT_0565 + IT_0567 + IT_0569;
    const ccomplex_t IT_0571 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0572 = IT_0248*IT_0571;
    const ccomplex_t IT_0573 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0574 = IT_0251*IT_0573;
    const ccomplex_t IT_0575 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0576 = IT_0248*IT_0575;
    const ccomplex_t IT_0577 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0578 = IT_0251*IT_0577;
    const ccomplex_t IT_0579 = -IT_0572 + -IT_0574 + -IT_0576 + -IT_0578;
    const ccomplex_t IT_0580 = IT_0570 + IT_0579;
    const ccomplex_t IT_0581 = (0 + _Complex_I*1)*IT_0580;
    const ccomplex_t IT_0582 = (-0.5)*IT_0581;
    const ccomplex_t IT_0583 = m_N_2*m_N_4;
    const ccomplex_t IT_0584 = cpow(IT_0017 + IT_0290 + -IT_0487, 2);
    const ccomplex_t IT_0585 = -IT_0290 + 0.25*IT_0001*(IT_0017 + -IT_0290 +
       IT_0487)*(IT_0017 + IT_0290 + -IT_0487) + 0.25*IT_0001*IT_0584;
    const ccomplex_t IT_0586 = 2*IT_0585;
    const ccomplex_t IT_0587 = cabs(-IT_0290 + 0.25*IT_0001*IT_0584);
    const ccomplex_t IT_0588 = cpow(IT_0587, 0.5);
    const ccomplex_t IT_0589 = ( ((m_h
) >= (m_N_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0561)*(IT_0582*IT_0583 + (-0.5)*IT_0561*IT_0586)
       + conj(IT_0582)*(IT_0561*IT_0583 + (-0.5)*IT_0582*IT_0586))*IT_0588
) : 
      (0
));
    const ccomplex_t IT_0590 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0591 = IT_0242*IT_0590;
    const ccomplex_t IT_0592 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0593 = IT_0245*IT_0592;
    const ccomplex_t IT_0594 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0595 = IT_0242*IT_0594;
    const ccomplex_t IT_0596 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0597 = IT_0245*IT_0596;
    const ccomplex_t IT_0598 = IT_0591 + IT_0593 + IT_0595 + IT_0597;
    const ccomplex_t IT_0599 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0600 = IT_0248*IT_0599;
    const ccomplex_t IT_0601 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0602 = IT_0251*IT_0601;
    const ccomplex_t IT_0603 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0604 = IT_0248*IT_0603;
    const ccomplex_t IT_0605 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0606 = IT_0251*IT_0605;
    const ccomplex_t IT_0607 = -IT_0600 + -IT_0602 + -IT_0604 + -IT_0606;
    const ccomplex_t IT_0608 = IT_0598 + IT_0607;
    const ccomplex_t IT_0609 = (0 + _Complex_I*1)*IT_0608;
    const ccomplex_t IT_0610 = (-0.5)*IT_0609;
    const ccomplex_t IT_0611 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0612 = IT_0242*IT_0611;
    const ccomplex_t IT_0613 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0614 = IT_0245*IT_0613;
    const ccomplex_t IT_0615 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0616 = IT_0242*IT_0615;
    const ccomplex_t IT_0617 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0618 = IT_0245*IT_0617;
    const ccomplex_t IT_0619 = IT_0612 + IT_0614 + IT_0616 + IT_0618;
    const ccomplex_t IT_0620 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0621 = IT_0248*IT_0620;
    const ccomplex_t IT_0622 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0623 = IT_0251*IT_0622;
    const ccomplex_t IT_0624 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0625 = IT_0248*IT_0624;
    const ccomplex_t IT_0626 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0627 = IT_0251*IT_0626;
    const ccomplex_t IT_0628 = -IT_0621 + -IT_0623 + -IT_0625 + -IT_0627;
    const ccomplex_t IT_0629 = IT_0619 + IT_0628;
    const ccomplex_t IT_0630 = (0 + _Complex_I*1)*IT_0629;
    const ccomplex_t IT_0631 = (-0.5)*IT_0630;
    const ccomplex_t IT_0632 = m_N_3*m_N_4;
    const ccomplex_t IT_0633 = cpow(IT_0017 + IT_0364 + -IT_0487, 2);
    const ccomplex_t IT_0634 = -IT_0364 + 0.25*IT_0001*(IT_0017 + -IT_0364 +
       IT_0487)*(IT_0017 + IT_0364 + -IT_0487) + 0.25*IT_0001*IT_0633;
    const ccomplex_t IT_0635 = 2*IT_0634;
    const ccomplex_t IT_0636 = cabs(-IT_0364 + 0.25*IT_0001*IT_0633);
    const ccomplex_t IT_0637 = cpow(IT_0636, 0.5);
    const ccomplex_t IT_0638 = ( ((m_h
) >= (m_N_3 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0610)*(IT_0631*IT_0632 + (-0.5)*IT_0610*IT_0635)
       + conj(IT_0631)*(IT_0610*IT_0632 + (-0.5)*IT_0631*IT_0635))*IT_0637
) : 
      (0
));
    const ccomplex_t IT_0639 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0009
      *IT_0021*IT_0022;
    const ccomplex_t IT_0640 = 0.5*IT_0639;
    const ccomplex_t IT_0641 = cpow(IT_0640, 2);
    const ccomplex_t IT_0642 = pow(m_tau, 2);
    const ccomplex_t IT_0643 = 0.5*IT_0017 + -IT_0642;
    const ccomplex_t IT_0644 = cabs(0.25*IT_0017 + -IT_0642);
    const ccomplex_t IT_0645 = cpow(IT_0644, 0.5);
    const ccomplex_t IT_0646 = ( ((m_h
) >= (2*m_tau
)) ? (0.5*IT_0000*IT_0001
      *IT_0641*(IT_0642 + -IT_0643)*IT_0645
) : (0
));
    const ccomplex_t IT_0647 = cpow(V_tb, 2);
    const ccomplex_t IT_0648 = m_b*mu_h*IT_0021;
    const ccomplex_t IT_0649 = IT_0008*IT_0009;
    const ccomplex_t IT_0650 = IT_0011*IT_0012;
    const ccomplex_t IT_0651 = IT_0649 + IT_0650;
    const ccomplex_t IT_0652 = 2*U_sb_00*((IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0026 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0651)
      *conj(U_sb_00) + 1./2*IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_10)) +
       U_sb_10*(IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_00) + 2*(IT_0006
      *IT_0009*IT_0021*IT_0022*IT_0026 + -1./3*m_W*IT_0003*IT_0004*IT_0651)*conj
      (U_sb_10));
    const ccomplex_t IT_0653 = (0 + _Complex_I*1)*(e_em*IT_0652 +
       1.4142135623731*A_b*IT_0009*IT_0647*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0654 = 0.5*IT_0653;
    const ccomplex_t IT_0655 = pow(m_sb_1, 2);
    const ccomplex_t IT_0656 = cabs(0.25*IT_0017 + -IT_0655);
    const ccomplex_t IT_0657 = cpow(IT_0656, 0.5);
    const ccomplex_t IT_0658 = ( ((m_h
) >= (2*m_sb_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0654*conj(IT_0654)*IT_0657
) : (0
));
    const ccomplex_t IT_0659 = 2*U_sb_01*((IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0026 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0651)
      *conj(U_sb_01) + 1./2*IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_11)) +
       U_sb_11*(IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_01) + 2*(IT_0006
      *IT_0009*IT_0021*IT_0022*IT_0026 + -1./3*m_W*IT_0003*IT_0004*IT_0651)*conj
      (U_sb_11));
    const ccomplex_t IT_0660 = (0 + _Complex_I*1)*(e_em*IT_0659 +
       1.4142135623731*A_b*IT_0009*IT_0647*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0661 = 0.5*IT_0660;
    const ccomplex_t IT_0662 = pow(m_sb_2, 2);
    const ccomplex_t IT_0663 = cabs(0.25*IT_0017 + -IT_0662);
    const ccomplex_t IT_0664 = cpow(IT_0663, 0.5);
    const ccomplex_t IT_0665 = ( ((m_h
) >= (2*m_sb_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0661*conj(IT_0661)*IT_0664
) : (0
));
    const ccomplex_t IT_0666 = 2*U_sb_00*((IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0026 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)*IT_0651)
      *conj(U_sb_01) + 1./2*IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_11)) +
       U_sb_10*(IT_0006*IT_0011*IT_0022*IT_0648*conj(U_sb_01) + 2*(IT_0006
      *IT_0009*IT_0021*IT_0022*IT_0026 + -1./3*m_W*IT_0003*IT_0004*IT_0651)*conj
      (U_sb_11));
    const ccomplex_t IT_0667 = (0 + _Complex_I*1)*(e_em*IT_0666 +
       1.4142135623731*A_b*IT_0009*IT_0647*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0668 = 0.5*IT_0667;
    const ccomplex_t IT_0669 = cpow(IT_0017 + IT_0655 + -IT_0662, 2);
    const ccomplex_t IT_0670 = cabs(-IT_0655 + 0.25*IT_0001*IT_0669);
    const ccomplex_t IT_0671 = cpow(IT_0670, 0.5);
    const ccomplex_t IT_0672 = ( ((m_h
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0668*conj(IT_0668)*IT_0671
) : (0
));
    const ccomplex_t IT_0673 = 2*conj(U_sb_00)*((IT_0006*IT_0009*IT_0021
      *IT_0022*IT_0026 + (-0.166666666666667)*m_W*(IT_0003*IT_0004 + 3*IT_0006)
      *IT_0651)*U_sb_01 + 1./2*IT_0006*IT_0011*IT_0022*IT_0648*U_sb_11) + conj
      (U_sb_10)*(IT_0006*IT_0011*IT_0022*IT_0648*U_sb_01 + 2*(IT_0006*IT_0009
      *IT_0021*IT_0022*IT_0026 + -1./3*m_W*IT_0003*IT_0004*IT_0651)*U_sb_11);
    const ccomplex_t IT_0674 = (0 + _Complex_I*1)*(e_em*IT_0673 +
       1.4142135623731*A_b*IT_0009*IT_0647*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0675 = 0.5*IT_0674;
    const ccomplex_t IT_0676 = ( ((m_h
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0671*IT_0675*conj(IT_0675)
) : (0
));
    const ccomplex_t IT_0677 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0021*IT_0031*IT_0035);
    const ccomplex_t IT_0678 = -IT_0677;
    const ccomplex_t IT_0679 = cpow(IT_0678, 2);
    const ccomplex_t IT_0680 = pow(m_sc_L, 2);
    const ccomplex_t IT_0681 = cabs(0.25*IT_0017 + -IT_0680);
    const ccomplex_t IT_0682 = cpow(IT_0681, 0.5);
    const ccomplex_t IT_0683 = ( ((m_h
) >= (2*m_sc_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0679*IT_0682
) : (0
));
    const ccomplex_t IT_0684 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0011*IT_0021
      *IT_0031*IT_0035);
    const ccomplex_t IT_0685 = 0.666666666666667*IT_0684;
    const ccomplex_t IT_0686 = cpow(IT_0685, 2);
    const ccomplex_t IT_0687 = pow(m_sc_R, 2);
    const ccomplex_t IT_0688 = cabs(0.25*IT_0017 + -IT_0687);
    const ccomplex_t IT_0689 = cpow(IT_0688, 0.5);
    const ccomplex_t IT_0690 = ( ((m_h
) >= (2*m_sc_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0686*IT_0689
) : (0
));
    const ccomplex_t IT_0691 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0006
      *IT_0009*IT_0021*IT_0031;
    const ccomplex_t IT_0692 = (-0.5)*IT_0691;
    const ccomplex_t IT_0693 = cpow(IT_0692, 2);
    const ccomplex_t IT_0694 = cpow(IT_0017 + IT_0680 + -IT_0687, 2);
    const ccomplex_t IT_0695 = cabs(-IT_0680 + 0.25*IT_0001*IT_0694);
    const ccomplex_t IT_0696 = cpow(IT_0695, 0.5);
    const ccomplex_t IT_0697 = ( ((m_h
) >= (m_sc_L + m_sc_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0693*IT_0696
) : (0
));
    const ccomplex_t IT_0698 = (0 + _Complex_I*1)*e_em*((-0.166666666666667)
      *m_W*(IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0021*IT_0022*IT_0043);
    const ccomplex_t IT_0699 = cpow(IT_0698, 2);
    const ccomplex_t IT_0700 = pow(m_sd_L, 2);
    const ccomplex_t IT_0701 = cabs(0.25*IT_0017 + -IT_0700);
    const ccomplex_t IT_0702 = cpow(IT_0701, 0.5);
    const ccomplex_t IT_0703 = ( ((m_h
) >= (2*m_sd_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0699*IT_0702
) : (0
));
    const ccomplex_t IT_0704 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_10);
    const ccomplex_t IT_0705 = 0.5*IT_0704;
    const ccomplex_t IT_0706 = cpow(IT_0017 + -IT_0655 + IT_0700, 2);
    const ccomplex_t IT_0707 = cabs(-IT_0700 + 0.25*IT_0001*IT_0706);
    const ccomplex_t IT_0708 = cpow(IT_0707, 0.5);
    const ccomplex_t IT_0709 = ( ((m_h
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0705*conj(IT_0705)*IT_0708
) : (0
));
    const ccomplex_t IT_0710 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_10;
    const ccomplex_t IT_0711 = 0.5*IT_0710;
    const ccomplex_t IT_0712 = ( ((m_h
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0708*IT_0711*conj(IT_0711)
) : (0
));
    const ccomplex_t IT_0713 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_11);
    const ccomplex_t IT_0714 = 0.5*IT_0713;
    const ccomplex_t IT_0715 = cpow(IT_0017 + -IT_0662 + IT_0700, 2);
    const ccomplex_t IT_0716 = cabs(-IT_0700 + 0.25*IT_0001*IT_0715);
    const ccomplex_t IT_0717 = cpow(IT_0716, 0.5);
    const ccomplex_t IT_0718 = ( ((m_h
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0714*conj(IT_0714)*IT_0717
) : (0
));
    const ccomplex_t IT_0719 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_11;
    const ccomplex_t IT_0720 = 0.5*IT_0719;
    const ccomplex_t IT_0721 = ( ((m_h
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0717*IT_0720*conj(IT_0720)
) : (0
));
    const ccomplex_t IT_0722 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + (-3)*IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0043);
    const ccomplex_t IT_0723 = (-0.333333333333333)*IT_0722;
    const ccomplex_t IT_0724 = cpow(IT_0723, 2);
    const ccomplex_t IT_0725 = pow(m_sd_R, 2);
    const ccomplex_t IT_0726 = cabs(0.25*IT_0017 + -IT_0725);
    const ccomplex_t IT_0727 = cpow(IT_0726, 0.5);
    const ccomplex_t IT_0728 = ( ((m_h
) >= (2*m_sd_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0724*IT_0727
) : (0
));
    const ccomplex_t IT_0729 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_00);
    const ccomplex_t IT_0730 = 0.5*IT_0729;
    const ccomplex_t IT_0731 = cpow(IT_0017 + -IT_0655 + IT_0725, 2);
    const ccomplex_t IT_0732 = cabs(-IT_0725 + 0.25*IT_0001*IT_0731);
    const ccomplex_t IT_0733 = cpow(IT_0732, 0.5);
    const ccomplex_t IT_0734 = ( ((m_h
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0730*conj(IT_0730)*IT_0733
) : (0
));
    const ccomplex_t IT_0735 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_00;
    const ccomplex_t IT_0736 = 0.5*IT_0735;
    const ccomplex_t IT_0737 = ( ((m_h
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0733*IT_0736*conj(IT_0736)
) : (0
));
    const ccomplex_t IT_0738 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_01);
    const ccomplex_t IT_0739 = 0.5*IT_0738;
    const ccomplex_t IT_0740 = cpow(IT_0017 + -IT_0662 + IT_0725, 2);
    const ccomplex_t IT_0741 = cabs(-IT_0725 + 0.25*IT_0001*IT_0740);
    const ccomplex_t IT_0742 = cpow(IT_0741, 0.5);
    const ccomplex_t IT_0743 = ( ((m_h
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0739*conj(IT_0739)*IT_0742
) : (0
));
    const ccomplex_t IT_0744 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_01;
    const ccomplex_t IT_0745 = 0.5*IT_0744;
    const ccomplex_t IT_0746 = ( ((m_h
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0742*IT_0745*conj(IT_0745)
) : (0
));
    const ccomplex_t IT_0747 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0009 + m_d*e_em*mu_h*IT_0006*IT_0011*IT_0021*IT_0022);
    const ccomplex_t IT_0748 = 0.5*IT_0747;
    const ccomplex_t IT_0749 = cpow(IT_0748, 2);
    const ccomplex_t IT_0750 = cpow(IT_0017 + -IT_0700 + IT_0725, 2);
    const ccomplex_t IT_0751 = cabs(-IT_0725 + 0.25*IT_0001*IT_0750);
    const ccomplex_t IT_0752 = cpow(IT_0751, 0.5);
    const ccomplex_t IT_0753 = ( ((m_h
) >= (m_sd_L + m_sd_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0749*IT_0752
) : (0
));
    const ccomplex_t IT_0754 = (0 + _Complex_I*1)*e_em*(0.5*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) + IT_0006*IT_0009
      *IT_0021*IT_0022*IT_0051);
    const ccomplex_t IT_0755 = cpow(IT_0754, 2);
    const ccomplex_t IT_0756 = pow(m_se_L, 2);
    const ccomplex_t IT_0757 = cabs(0.25*IT_0017 + -IT_0756);
    const ccomplex_t IT_0758 = cpow(IT_0757, 0.5);
    const ccomplex_t IT_0759 = ( ((m_h
) >= (2*m_se_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0755*IT_0758
) : (0
));
    const ccomplex_t IT_0760 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + -IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0051);
    const ccomplex_t IT_0761 = -IT_0760;
    const ccomplex_t IT_0762 = cpow(IT_0761, 2);
    const ccomplex_t IT_0763 = pow(m_se_R, 2);
    const ccomplex_t IT_0764 = cabs(0.25*IT_0017 + -IT_0763);
    const ccomplex_t IT_0765 = cpow(IT_0764, 0.5);
    const ccomplex_t IT_0766 = ( ((m_h
) >= (2*m_se_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0762*IT_0765
) : (0
));
    const ccomplex_t IT_0767 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0006
      *IT_0011*IT_0021*IT_0022;
    const ccomplex_t IT_0768 = 0.5*IT_0767;
    const ccomplex_t IT_0769 = cpow(IT_0768, 2);
    const ccomplex_t IT_0770 = cpow(IT_0017 + -IT_0756 + IT_0763, 2);
    const ccomplex_t IT_0771 = cabs(-IT_0763 + 0.25*IT_0001*IT_0770);
    const ccomplex_t IT_0772 = cpow(IT_0771, 0.5);
    const ccomplex_t IT_0773 = ( ((m_h
) >= (m_se_L + m_se_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0769*IT_0772
) : (0
));
    const ccomplex_t IT_0774 = (0 + _Complex_I*1)*e_em*((-0.166666666666667)
      *m_W*(IT_0003*IT_0004 + 3*IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) +
       IT_0006*IT_0009*IT_0021*IT_0022*IT_0059);
    const ccomplex_t IT_0775 = cpow(IT_0774, 2);
    const ccomplex_t IT_0776 = pow(m_ss_L, 2);
    const ccomplex_t IT_0777 = cabs(0.25*IT_0017 + -IT_0776);
    const ccomplex_t IT_0778 = cpow(IT_0777, 0.5);
    const ccomplex_t IT_0779 = ( ((m_h
) >= (2*m_ss_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0775*IT_0778
) : (0
));
    const ccomplex_t IT_0780 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_10);
    const ccomplex_t IT_0781 = 0.5*IT_0780;
    const ccomplex_t IT_0782 = cpow(IT_0017 + -IT_0655 + IT_0776, 2);
    const ccomplex_t IT_0783 = cabs(-IT_0776 + 0.25*IT_0001*IT_0782);
    const ccomplex_t IT_0784 = cpow(IT_0783, 0.5);
    const ccomplex_t IT_0785 = ( ((m_h
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0781*conj(IT_0781)*IT_0784
) : (0
));
    const ccomplex_t IT_0786 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_10;
    const ccomplex_t IT_0787 = 0.5*IT_0786;
    const ccomplex_t IT_0788 = ( ((m_h
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0784*IT_0787*conj(IT_0787)
) : (0
));
    const ccomplex_t IT_0789 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_11);
    const ccomplex_t IT_0790 = 0.5*IT_0789;
    const ccomplex_t IT_0791 = cpow(IT_0017 + -IT_0662 + IT_0776, 2);
    const ccomplex_t IT_0792 = cabs(-IT_0776 + 0.25*IT_0001*IT_0791);
    const ccomplex_t IT_0793 = cpow(IT_0792, 0.5);
    const ccomplex_t IT_0794 = ( ((m_h
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0790*conj(IT_0790)*IT_0793
) : (0
));
    const ccomplex_t IT_0795 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_11;
    const ccomplex_t IT_0796 = 0.5*IT_0795;
    const ccomplex_t IT_0797 = ( ((m_h
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0793*IT_0796*conj(IT_0796)
) : (0
));
    const ccomplex_t IT_0798 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0009;
    const ccomplex_t IT_0799 = 0.5*IT_0798;
    const ccomplex_t IT_0800 = cpow(IT_0017 + IT_0725 + -IT_0776, 2);
    const ccomplex_t IT_0801 = cabs(-IT_0725 + 0.25*IT_0001*IT_0800);
    const ccomplex_t IT_0802 = cpow(IT_0801, 0.5);
    const ccomplex_t IT_0803 = ( ((m_h
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0799*conj(IT_0799)*IT_0802
) : (0
));
    const ccomplex_t IT_0804 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0009;
    const ccomplex_t IT_0805 = 0.5*IT_0804;
    const ccomplex_t IT_0806 = ( ((m_h
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0802*IT_0805*conj(IT_0805)
) : (0
));
    const ccomplex_t IT_0807 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + (-3)*IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0059);
    const ccomplex_t IT_0808 = (-0.333333333333333)*IT_0807;
    const ccomplex_t IT_0809 = cpow(IT_0808, 2);
    const ccomplex_t IT_0810 = pow(m_ss_R, 2);
    const ccomplex_t IT_0811 = cabs(0.25*IT_0017 + -IT_0810);
    const ccomplex_t IT_0812 = cpow(IT_0811, 0.5);
    const ccomplex_t IT_0813 = ( ((m_h
) >= (2*m_ss_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0809*IT_0812
) : (0
));
    const ccomplex_t IT_0814 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_00);
    const ccomplex_t IT_0815 = 0.5*IT_0814;
    const ccomplex_t IT_0816 = cpow(IT_0017 + -IT_0655 + IT_0810, 2);
    const ccomplex_t IT_0817 = cabs(-IT_0810 + 0.25*IT_0001*IT_0816);
    const ccomplex_t IT_0818 = cpow(IT_0817, 0.5);
    const ccomplex_t IT_0819 = ( ((m_h
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0815*conj(IT_0815)*IT_0818
) : (0
));
    const ccomplex_t IT_0820 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_00;
    const ccomplex_t IT_0821 = 0.5*IT_0820;
    const ccomplex_t IT_0822 = ( ((m_h
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0818*IT_0821*conj(IT_0821)
) : (0
));
    const ccomplex_t IT_0823 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_01);
    const ccomplex_t IT_0824 = 0.5*IT_0823;
    const ccomplex_t IT_0825 = cpow(IT_0017 + -IT_0662 + IT_0810, 2);
    const ccomplex_t IT_0826 = cabs(-IT_0810 + 0.25*IT_0001*IT_0825);
    const ccomplex_t IT_0827 = cpow(IT_0826, 0.5);
    const ccomplex_t IT_0828 = ( ((m_h
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0824*conj(IT_0824)*IT_0827
) : (0
));
    const ccomplex_t IT_0829 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0009*U_sb_01;
    const ccomplex_t IT_0830 = 0.5*IT_0829;
    const ccomplex_t IT_0831 = ( ((m_h
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0827*IT_0830*conj(IT_0830)
) : (0
));
    const ccomplex_t IT_0832 = cpow(IT_0017 + -IT_0700 + IT_0810, 2);
    const ccomplex_t IT_0833 = cabs(-IT_0810 + 0.25*IT_0001*IT_0832);
    const ccomplex_t IT_0834 = cpow(IT_0833, 0.5);
    const ccomplex_t IT_0835 = ( ((m_h
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0799*conj(IT_0799)*IT_0834
) : (0
));
    const ccomplex_t IT_0836 = ( ((m_h
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0805*conj(IT_0805)*IT_0834
) : (0
));
    const ccomplex_t IT_0837 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0009 + m_s*e_em*mu_h*IT_0006*IT_0011*IT_0021*IT_0022);
    const ccomplex_t IT_0838 = 0.5*IT_0837;
    const ccomplex_t IT_0839 = cpow(IT_0838, 2);
    const ccomplex_t IT_0840 = cpow(IT_0017 + -IT_0776 + IT_0810, 2);
    const ccomplex_t IT_0841 = cabs(-IT_0810 + 0.25*IT_0001*IT_0840);
    const ccomplex_t IT_0842 = cpow(IT_0841, 0.5);
    const ccomplex_t IT_0843 = ( ((m_h
) >= (m_ss_L + m_ss_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0839*IT_0842
) : (0
));
    const ccomplex_t IT_0844 = m_t*mu_h*IT_0021;
    const ccomplex_t IT_0845 = U_st_10*(IT_0006*IT_0009*IT_0031*IT_0844*conj
      (U_st_00) + (-1.33333333333333)*((-1.5)*IT_0006*IT_0011*IT_0021*IT_0031
      *IT_0067 + m_W*IT_0003*IT_0004*IT_0651)*conj(U_st_10)) + 2*U_st_00*(
      (IT_0006*IT_0011*IT_0021*IT_0031*IT_0067 + 0.166666666666667*m_W*(IT_0003
      *IT_0004 + (-3)*IT_0006)*IT_0651)*conj(U_st_00) + 1./2*IT_0006*IT_0009
      *IT_0031*IT_0844*conj(U_st_10));
    const ccomplex_t IT_0846 = (0 + _Complex_I*1)*(e_em*IT_0845 +
       1.4142135623731*A_t*IT_0011*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0847 = (-0.5)*IT_0846;
    const ccomplex_t IT_0848 = pow(m_st_1, 2);
    const ccomplex_t IT_0849 = cabs(0.25*IT_0017 + -IT_0848);
    const ccomplex_t IT_0850 = cpow(IT_0849, 0.5);
    const ccomplex_t IT_0851 = ( ((m_h
) >= (2*m_st_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0847*conj(IT_0847)*IT_0850
) : (0
));
    const ccomplex_t IT_0852 = U_st_11*(IT_0006*IT_0009*IT_0031*IT_0844*conj
      (U_st_01) + (-1.33333333333333)*((-1.5)*IT_0006*IT_0011*IT_0021*IT_0031
      *IT_0067 + m_W*IT_0003*IT_0004*IT_0651)*conj(U_st_11)) + 2*U_st_01*(
      (IT_0006*IT_0011*IT_0021*IT_0031*IT_0067 + 0.166666666666667*m_W*(IT_0003
      *IT_0004 + (-3)*IT_0006)*IT_0651)*conj(U_st_01) + 1./2*IT_0006*IT_0009
      *IT_0031*IT_0844*conj(U_st_11));
    const ccomplex_t IT_0853 = (0 + _Complex_I*1)*(e_em*IT_0852 +
       1.4142135623731*A_t*IT_0011*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0854 = (-0.5)*IT_0853;
    const ccomplex_t IT_0855 = pow(m_st_2, 2);
    const ccomplex_t IT_0856 = cabs(0.25*IT_0017 + -IT_0855);
    const ccomplex_t IT_0857 = cpow(IT_0856, 0.5);
    const ccomplex_t IT_0858 = ( ((m_h
) >= (2*m_st_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0854*conj(IT_0854)*IT_0857
) : (0
));
    const ccomplex_t IT_0859 = U_st_10*(IT_0006*IT_0009*IT_0031*IT_0844*conj
      (U_st_01) + (-1.33333333333333)*((-1.5)*IT_0006*IT_0011*IT_0021*IT_0031
      *IT_0067 + m_W*IT_0003*IT_0004*IT_0651)*conj(U_st_11)) + 2*U_st_00*(
      (IT_0006*IT_0011*IT_0021*IT_0031*IT_0067 + 0.166666666666667*m_W*(IT_0003
      *IT_0004 + (-3)*IT_0006)*IT_0651)*conj(U_st_01) + 1./2*IT_0006*IT_0009
      *IT_0031*IT_0844*conj(U_st_11));
    const ccomplex_t IT_0860 = (0 + _Complex_I*1)*(e_em*IT_0859 +
       1.4142135623731*A_t*IT_0011*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0861 = (-0.5)*IT_0860;
    const ccomplex_t IT_0862 = cpow(IT_0017 + IT_0848 + -IT_0855, 2);
    const ccomplex_t IT_0863 = cabs(-IT_0848 + 0.25*IT_0001*IT_0862);
    const ccomplex_t IT_0864 = cpow(IT_0863, 0.5);
    const ccomplex_t IT_0865 = ( ((m_h
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0861*conj(IT_0861)*IT_0864
) : (0
));
    const ccomplex_t IT_0866 = conj(U_st_10)*(IT_0006*IT_0009*IT_0031*IT_0844
      *U_st_01 + (-1.33333333333333)*((-1.5)*IT_0006*IT_0011*IT_0021*IT_0031
      *IT_0067 + m_W*IT_0003*IT_0004*IT_0651)*U_st_11) + 2*conj(U_st_00)*(
      (IT_0006*IT_0011*IT_0021*IT_0031*IT_0067 + 0.166666666666667*m_W*(IT_0003
      *IT_0004 + (-3)*IT_0006)*IT_0651)*U_st_01 + 1./2*IT_0006*IT_0009*IT_0031
      *IT_0844*U_st_11);
    const ccomplex_t IT_0867 = (0 + _Complex_I*1)*(e_em*IT_0866 +
       1.4142135623731*A_t*IT_0011*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0868 = (-0.5)*IT_0867;
    const ccomplex_t IT_0869 = ( ((m_h
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0864*IT_0868*conj(IT_0868)
) : (0
));
    const ccomplex_t IT_0870 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0003*IT_0004 + (-3)*IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) +
       IT_0006*IT_0011*IT_0021*IT_0031*IT_0075);
    const ccomplex_t IT_0871 = -IT_0870;
    const ccomplex_t IT_0872 = cpow(IT_0871, 2);
    const ccomplex_t IT_0873 = pow(m_su_L, 2);
    const ccomplex_t IT_0874 = cabs(0.25*IT_0017 + -IT_0873);
    const ccomplex_t IT_0875 = cpow(IT_0874, 0.5);
    const ccomplex_t IT_0876 = ( ((m_h
) >= (2*m_su_L
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0872*IT_0875
) : (0
));
    const ccomplex_t IT_0877 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + (-1.5)*IT_0006*IT_0011*IT_0021
      *IT_0031*IT_0075);
    const ccomplex_t IT_0878 = 0.666666666666667*IT_0877;
    const ccomplex_t IT_0879 = cpow(IT_0878, 2);
    const ccomplex_t IT_0880 = pow(m_su_R, 2);
    const ccomplex_t IT_0881 = cabs(0.25*IT_0017 + -IT_0880);
    const ccomplex_t IT_0882 = cpow(IT_0881, 0.5);
    const ccomplex_t IT_0883 = ( ((m_h
) >= (2*m_su_R
)) ? ((-0.375)*IT_0000
      *IT_0001*IT_0879*IT_0882
) : (0
));
    const ccomplex_t IT_0884 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0006
      *IT_0009*IT_0021*IT_0031;
    const ccomplex_t IT_0885 = (-0.5)*IT_0884;
    const ccomplex_t IT_0886 = cpow(IT_0885, 2);
    const ccomplex_t IT_0887 = cpow(IT_0017 + IT_0873 + -IT_0880, 2);
    const ccomplex_t IT_0888 = cabs(-IT_0873 + 0.25*IT_0001*IT_0887);
    const ccomplex_t IT_0889 = cpow(IT_0888, 0.5);
    const ccomplex_t IT_0890 = ( ((m_h
) >= (m_su_L + m_su_R
)) ? ((-0.375)
      *IT_0000*IT_0001*IT_0886*IT_0889
) : (0
));
    const ccomplex_t IT_0891 = (0 + _Complex_I*1)*e_em*(0.5*m_W*(IT_0003
      *IT_0004 + -IT_0006)*(IT_0008*IT_0009 + IT_0011*IT_0012) + IT_0006*IT_0009
      *IT_0021*IT_0022*IT_0144);
    const ccomplex_t IT_0892 = cpow(IT_0891, 2);
    const ccomplex_t IT_0893 = pow(m_smu_L, 2);
    const ccomplex_t IT_0894 = cabs(0.25*IT_0017 + -IT_0893);
    const ccomplex_t IT_0895 = cpow(IT_0894, 0.5);
    const ccomplex_t IT_0896 = ( ((m_h
) >= (2*m_smu_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0892*IT_0895
) : (0
));
    const ccomplex_t IT_0897 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*IT_0004*
      (IT_0008*IT_0009 + IT_0011*IT_0012) + -IT_0006*IT_0009*IT_0021*IT_0022
      *IT_0144);
    const ccomplex_t IT_0898 = -IT_0897;
    const ccomplex_t IT_0899 = cpow(IT_0898, 2);
    const ccomplex_t IT_0900 = pow(m_smu_R, 2);
    const ccomplex_t IT_0901 = cabs(0.25*IT_0017 + -IT_0900);
    const ccomplex_t IT_0902 = cpow(IT_0901, 0.5);
    const ccomplex_t IT_0903 = ( ((m_h
) >= (2*m_smu_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0899*IT_0902
) : (0
));
    const ccomplex_t IT_0904 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0006
      *IT_0011*IT_0021*IT_0022;
    const ccomplex_t IT_0905 = 0.5*IT_0904;
    const ccomplex_t IT_0906 = cpow(IT_0905, 2);
    const ccomplex_t IT_0907 = cpow(IT_0017 + -IT_0893 + IT_0900, 2);
    const ccomplex_t IT_0908 = cabs(-IT_0900 + 0.25*IT_0001*IT_0907);
    const ccomplex_t IT_0909 = cpow(IT_0908, 0.5);
    const ccomplex_t IT_0910 = ( ((m_h
) >= (m_smu_L + m_smu_R
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0906*IT_0909
) : (0
));
    const ccomplex_t IT_0911 = (0 + _Complex_I*1)*m_W*e_em*IT_0007*(IT_0008
      *IT_0009 + IT_0011*IT_0012);
    const ccomplex_t IT_0912 = 0.5*IT_0911;
    const ccomplex_t IT_0913 = cpow(IT_0912, 2);
    const ccomplex_t IT_0914 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0915 = cabs(0.25*IT_0017 + -IT_0914);
    const ccomplex_t IT_0916 = cpow(IT_0915, 0.5);
    const ccomplex_t IT_0917 = ( ((m_h
) >= (2*m_snu_mu
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0913*IT_0916
) : (0
));
    const ccomplex_t IT_0918 = mu_h*m_tau*IT_0021;
    const ccomplex_t IT_0919 = U_stau_10*(IT_0006*IT_0011*IT_0022*IT_0918*conj
      (U_stau_00) + 2*(IT_0006*IT_0009*IT_0021*IT_0022*IT_0642 + -m_W*IT_0003
      *IT_0004*IT_0651)*conj(U_stau_10)) + 2*U_stau_00*((IT_0006*IT_0009*IT_0021
      *IT_0022*IT_0642 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0651)*conj
      (U_stau_00) + 1./2*IT_0006*IT_0011*IT_0022*IT_0918*conj(U_stau_10));
    const ccomplex_t IT_0920 = (0 + _Complex_I*1)*(e_em*IT_0919 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0921 = 0.5*IT_0920;
    const ccomplex_t IT_0922 = pow(m_stau_1, 2);
    const ccomplex_t IT_0923 = cabs(0.25*IT_0017 + -IT_0922);
    const ccomplex_t IT_0924 = cpow(IT_0923, 0.5);
    const ccomplex_t IT_0925 = ( ((m_h
) >= (2*m_stau_1
)) ? (0.125*IT_0000
      *IT_0001*IT_0921*conj(IT_0921)*IT_0924
) : (0
));
    const ccomplex_t IT_0926 = U_stau_11*(IT_0006*IT_0011*IT_0022*IT_0918*conj
      (U_stau_01) + 2*(IT_0006*IT_0009*IT_0021*IT_0022*IT_0642 + -m_W*IT_0003
      *IT_0004*IT_0651)*conj(U_stau_11)) + 2*U_stau_01*((IT_0006*IT_0009*IT_0021
      *IT_0022*IT_0642 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0651)*conj
      (U_stau_01) + 1./2*IT_0006*IT_0011*IT_0022*IT_0918*conj(U_stau_11));
    const ccomplex_t IT_0927 = (0 + _Complex_I*1)*(e_em*IT_0926 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0928 = 0.5*IT_0927;
    const ccomplex_t IT_0929 = pow(m_stau_2, 2);
    const ccomplex_t IT_0930 = cabs(0.25*IT_0017 + -IT_0929);
    const ccomplex_t IT_0931 = cpow(IT_0930, 0.5);
    const ccomplex_t IT_0932 = ( ((m_h
) >= (2*m_stau_2
)) ? (0.125*IT_0000
      *IT_0001*IT_0928*conj(IT_0928)*IT_0931
) : (0
));
    const ccomplex_t IT_0933 = U_stau_10*(IT_0006*IT_0011*IT_0022*IT_0918*conj
      (U_stau_01) + 2*(IT_0006*IT_0009*IT_0021*IT_0022*IT_0642 + -m_W*IT_0003
      *IT_0004*IT_0651)*conj(U_stau_11)) + 2*U_stau_00*((IT_0006*IT_0009*IT_0021
      *IT_0022*IT_0642 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0651)*conj
      (U_stau_01) + 1./2*IT_0006*IT_0011*IT_0022*IT_0918*conj(U_stau_11));
    const ccomplex_t IT_0934 = (0 + _Complex_I*1)*(e_em*IT_0933 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0935 = 0.5*IT_0934;
    const ccomplex_t IT_0936 = cpow(IT_0017 + IT_0922 + -IT_0929, 2);
    const ccomplex_t IT_0937 = cabs(-IT_0922 + 0.25*IT_0001*IT_0936);
    const ccomplex_t IT_0938 = cpow(IT_0937, 0.5);
    const ccomplex_t IT_0939 = ( ((m_h
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0935*conj(IT_0935)*IT_0938
) : (0
));
    const ccomplex_t IT_0940 = conj(U_stau_10)*(IT_0006*IT_0011*IT_0022
      *IT_0918*U_stau_01 + 2*(IT_0006*IT_0009*IT_0021*IT_0022*IT_0642 + -m_W
      *IT_0003*IT_0004*IT_0651)*U_stau_11) + 2*conj(U_stau_00)*((IT_0006*IT_0009
      *IT_0021*IT_0022*IT_0642 + 0.5*m_W*(IT_0003*IT_0004 + -IT_0006)*IT_0651)
      *U_stau_01 + 1./2*IT_0006*IT_0011*IT_0022*IT_0918*U_stau_11);
    const ccomplex_t IT_0941 = (0 + _Complex_I*1)*(e_em*IT_0940 +
       1.4142135623731*A_tau*IT_0009*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0942 = 0.5*IT_0941;
    const ccomplex_t IT_0943 = ( ((m_h
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0938*IT_0942*conj(IT_0942)
) : (0
));
    const ccomplex_t IT_0944 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0945 = cabs(0.25*IT_0017 + -IT_0944);
    const ccomplex_t IT_0946 = cpow(IT_0945, 0.5);
    const ccomplex_t IT_0947 = ( ((m_h
) >= (2*m_snu_tau
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0913*IT_0946
) : (0
));
    const ccomplex_t IT_0948 = (0 + _Complex_I*-1)*m_W*e_em*IT_0006*(IT_0008
      *IT_0009 + -IT_0011*IT_0012);
    const ccomplex_t IT_0949 = cpow(IT_0948, 2);
    const ccomplex_t IT_0950 = pow(m_W, -4);
    const ccomplex_t IT_0951 = 0.5*IT_0017 + -IT_0130;
    const ccomplex_t IT_0952 = cpow(IT_0951, 2);
    const ccomplex_t IT_0953 = IT_0950*IT_0952;
    const ccomplex_t IT_0954 = cabs(0.25*IT_0017 + -IT_0130);
    const ccomplex_t IT_0955 = cpow(IT_0954, 0.5);
    const ccomplex_t IT_0956 = pow(m_W, 4);
    const ccomplex_t IT_0957 = pow(m_h, -4);
    const ccomplex_t IT_0958 = clog(IT_0001*IT_0130);
    const ccomplex_t IT_0959 = cpow(-1 + 4*IT_0001*IT_0130, (-0.5));
    const ccomplex_t IT_0960 = pow(m_W, -3);
    const ccomplex_t IT_0961 = pow(m_h, 3);
    const ccomplex_t IT_0962 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0130)*IT_0960
      *IT_0961);
    const ccomplex_t IT_0963 = IT_0017*IT_0129*(1 + -IT_0001*IT_0130)*(1 + (
      -6.5)*IT_0001*IT_0130 + 23.5*IT_0956*IT_0957) + 1.5*(1 + (-6)*IT_0001
      *IT_0130 + 4*IT_0956*IT_0957)*IT_0958 + (-3)*(1 + (-8)*IT_0001*IT_0130 +
       20*IT_0956*IT_0957)*IT_0959*IT_0962;
    const ccomplex_t IT_0964 = ( ((m_h
) >= (2*m_W
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0949*(1 + 0.5*IT_0953)*IT_0955
) : (0.25*IT_0000*IT_0001
      *IT_0949*(1 + 0.5*IT_0953)*IT_0955*IT_0963
));
    const ccomplex_t IT_0965 = m_W*e_em;
    const ccomplex_t IT_0966 = cpow(IT_0002, -4);
    const ccomplex_t IT_0967 = IT_0004*IT_0011*IT_0012*IT_0966;
    const ccomplex_t IT_0968 = IT_0004*IT_0008*IT_0009*IT_0966;
    const ccomplex_t IT_0969 = cpow(IT_0004, -3);
    const ccomplex_t IT_0970 = IT_0011*IT_0012*IT_0969;
    const ccomplex_t IT_0971 = IT_0008*IT_0009*IT_0969;
    const ccomplex_t IT_0972 = IT_0003*IT_0006*IT_0011*IT_0012;
    const ccomplex_t IT_0973 = IT_0003*IT_0006*IT_0008*IT_0009;
    const ccomplex_t IT_0974 = cpow(IT_0004, -2);
    const ccomplex_t IT_0975 = cpow(IT_0003 + IT_0974, -1);
    const ccomplex_t IT_0976 = (IT_0967 + -IT_0968 + IT_0970 + -IT_0971 + 2
      *IT_0972 + (-2)*IT_0973)*IT_0975;
    const ccomplex_t IT_0977 = IT_0965*IT_0976;
    const ccomplex_t IT_0978 = (0 + _Complex_I*1)*IT_0977;
    const ccomplex_t IT_0979 = cpow(IT_0978, 2);
    const ccomplex_t IT_0980 = pow(m_Z, -4);
    const ccomplex_t IT_0981 = 0.5*IT_0017 + -IT_0094;
    const ccomplex_t IT_0982 = cpow(IT_0981, 2);
    const ccomplex_t IT_0983 = IT_0980*IT_0982;
    const ccomplex_t IT_0984 = cabs(0.25*IT_0017 + -IT_0094);
    const ccomplex_t IT_0985 = cpow(IT_0984, 0.5);
    const ccomplex_t IT_0986 = pow(m_Z, 4);
    const ccomplex_t IT_0987 = clog(IT_0001*IT_0094);
    const ccomplex_t IT_0988 = cpow(-1 + 4*IT_0001*IT_0094, (-0.5));
    const ccomplex_t IT_0989 = pow(m_Z, -3);
    const ccomplex_t IT_0990 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0094)*IT_0961
      *IT_0989);
    const ccomplex_t IT_0991 = IT_0017*IT_0093*(1 + -IT_0001*IT_0094)*(1 + (
      -6.5)*IT_0001*IT_0094 + 23.5*IT_0957*IT_0986) + 1.5*(1 + (-6)*IT_0001
      *IT_0094 + 4*IT_0957*IT_0986)*IT_0987 + (-3)*(1 + (-8)*IT_0001*IT_0094 +
       20*IT_0957*IT_0986)*IT_0988*IT_0990;
    const ccomplex_t IT_0992 = ( ((m_h
) >= (2*m_Z
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0979*(1 + 0.5*IT_0983)*IT_0985
) : (0.125*IT_0000*IT_0001
      *IT_0979*(1 + 0.5*IT_0983)*IT_0985*IT_0991
));
    const ccomplex_t IT_0993 = IT_0020 + IT_0030 + IT_0039 + IT_0047 + IT_0055
       + IT_0063 + IT_0071 + IT_0079 + IT_0088 + IT_0101 + IT_0108 + IT_0116 +
       IT_0125 + IT_0137 + IT_0140 + IT_0148 + IT_0172 + IT_0196 + IT_0221 +
       IT_0240 + IT_0269 + IT_0294 + IT_0343 + IT_0368 + IT_0417 + IT_0466 +
       IT_0491 + IT_0540 + IT_0589 + IT_0638 + IT_0646 + IT_0658 + IT_0665 +
       IT_0672 + IT_0676 + IT_0683 + IT_0690 + 2*IT_0697 + IT_0703 + IT_0709 +
       IT_0712 + IT_0718 + IT_0721 + IT_0728 + IT_0734 + IT_0737 + IT_0743 +
       IT_0746 + 2*IT_0753 + IT_0759 + IT_0766 + 2*IT_0773 + IT_0779 + IT_0785 +
       IT_0788 + IT_0794 + IT_0797 + IT_0803 + IT_0806 + IT_0813 + IT_0819 +
       IT_0822 + IT_0828 + IT_0831 + IT_0835 + IT_0836 + 2*IT_0843 + IT_0851 +
       IT_0858 + IT_0865 + IT_0869 + IT_0876 + IT_0883 + 2*IT_0890 + IT_0896 +
       IT_0903 + 2*IT_0910 + 2*IT_0917 + IT_0925 + IT_0932 + IT_0939 + IT_0943 +
       IT_0947 + IT_0964 + IT_0992;
    return create_ccomplex_return(IT_0993);
}

