#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_A0(
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
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
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
    const ccomplex_t IT_0001 = pow(m_A0, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = m_b*e_em*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = pow(m_b, 2);
    const ccomplex_t IT_0012 = pow(m_A0, 2);
    const ccomplex_t IT_0013 = -IT_0011 + 0.5*IT_0012;
    const ccomplex_t IT_0014 = cabs(-IT_0011 + 0.25*IT_0012);
    const ccomplex_t IT_0015 = cpow(IT_0014, 0.5);
    const ccomplex_t IT_0016 = ( ((m_A0
) >= (2*m_b
)) ? (0.125*IT_0000
      *IT_0001*IT_0010*(12*IT_0011 + 12*IT_0013)*IT_0015
) : (0
));
    const ccomplex_t IT_0017 = cpow(IT_0005, -1);
    const ccomplex_t IT_0018 = m_c*e_em*IT_0002*IT_0003*IT_0007*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = pow(m_c, 2);
    const ccomplex_t IT_0022 = 0.5*IT_0012 + -IT_0021;
    const ccomplex_t IT_0023 = cabs(0.25*IT_0012 + -IT_0021);
    const ccomplex_t IT_0024 = cpow(IT_0023, 0.5);
    const ccomplex_t IT_0025 = ( ((m_A0
) >= (2*m_c
)) ? (0.125*IT_0000
      *IT_0001*IT_0020*(12*IT_0021 + 12*IT_0022)*IT_0024
) : (0
));
    const ccomplex_t IT_0026 = m_d*e_em*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0027, 2);
    const ccomplex_t IT_0029 = pow(m_d, 2);
    const ccomplex_t IT_0030 = 0.5*IT_0012 + -IT_0029;
    const ccomplex_t IT_0031 = cabs(0.25*IT_0012 + -IT_0029);
    const ccomplex_t IT_0032 = cpow(IT_0031, 0.5);
    const ccomplex_t IT_0033 = ( ((m_A0
) >= (2*m_d
)) ? (0.125*IT_0000
      *IT_0001*IT_0028*(12*IT_0029 + 12*IT_0030)*IT_0032
) : (0
));
    const ccomplex_t IT_0034 = m_e*e_em*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = pow(m_e, 2);
    const ccomplex_t IT_0038 = 0.5*IT_0012 + -IT_0037;
    const ccomplex_t IT_0039 = cabs(0.25*IT_0012 + -IT_0037);
    const ccomplex_t IT_0040 = cpow(IT_0039, 0.5);
    const ccomplex_t IT_0041 = ( ((m_A0
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *IT_0036*(IT_0037 + IT_0038)*IT_0040
) : (0
));
    const ccomplex_t IT_0042 = m_s*e_em*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0043, 2);
    const ccomplex_t IT_0045 = pow(m_s, 2);
    const ccomplex_t IT_0046 = 0.5*IT_0012 + -IT_0045;
    const ccomplex_t IT_0047 = cabs(0.25*IT_0012 + -IT_0045);
    const ccomplex_t IT_0048 = cpow(IT_0047, 0.5);
    const ccomplex_t IT_0049 = ( ((m_A0
) >= (2*m_s
)) ? (0.125*IT_0000
      *IT_0001*IT_0044*(12*IT_0045 + 12*IT_0046)*IT_0048
) : (0
));
    const ccomplex_t IT_0050 = m_t*e_em*IT_0002*IT_0003*IT_0007*IT_0017;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0051, 2);
    const ccomplex_t IT_0053 = pow(m_t, 2);
    const ccomplex_t IT_0054 = 0.5*IT_0012 + -IT_0053;
    const ccomplex_t IT_0055 = cabs(0.25*IT_0012 + -IT_0053);
    const ccomplex_t IT_0056 = cpow(IT_0055, 0.5);
    const ccomplex_t IT_0057 = ( ((m_A0
) >= (2*m_t
)) ? (0.125*IT_0000
      *IT_0001*IT_0052*(12*IT_0053 + 12*IT_0054)*IT_0056
) : (0
));
    const ccomplex_t IT_0058 = m_u*e_em*IT_0002*IT_0003*IT_0007*IT_0017;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0059, 2);
    const ccomplex_t IT_0061 = pow(m_u, 2);
    const ccomplex_t IT_0062 = 0.5*IT_0012 + -IT_0061;
    const ccomplex_t IT_0063 = cabs(0.25*IT_0012 + -IT_0061);
    const ccomplex_t IT_0064 = cpow(IT_0063, 0.5);
    const ccomplex_t IT_0065 = ( ((m_A0
) >= (2*m_u
)) ? (0.125*IT_0000
      *IT_0001*IT_0060*(12*IT_0061 + 12*IT_0062)*IT_0064
) : (0
));
    const ccomplex_t IT_0066 = cos(theta_W);
    const ccomplex_t IT_0067 = cpow(IT_0066, -1);
    const ccomplex_t IT_0068 = cos(alpha);
    const ccomplex_t IT_0069 = sin(alpha);
    const ccomplex_t IT_0070 = e_em*IT_0007*IT_0067*(IT_0003*IT_0068 + IT_0005
      *IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = cpow(IT_0072, 2);
    const ccomplex_t IT_0074 = pow(m_Z, -2);
    const ccomplex_t IT_0075 = pow(m_Z, 2);
    const ccomplex_t IT_0076 = pow(m_h, 2);
    const ccomplex_t IT_0077 = 0.5*IT_0012 + 0.5*IT_0075 + (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = cpow(IT_0077, 2);
    const ccomplex_t IT_0079 = IT_0074*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0012 + IT_0075 + -IT_0076, 2);
    const ccomplex_t IT_0081 = cabs(-IT_0075 + 0.25*IT_0001*IT_0080);
    const ccomplex_t IT_0082 = cpow(IT_0081, 0.5);
    const ccomplex_t IT_0083 = ( ((m_A0
) >= (m_Z + m_h
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0073*(IT_0012 + -IT_0079)*IT_0082
) : (0
));
    const ccomplex_t IT_0084 = cpow(IT_0066, -2);
    const ccomplex_t IT_0085 = IT_0006*IT_0084;
    const ccomplex_t IT_0086 = IT_0007 + IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0005, 3);
    const ccomplex_t IT_0088 = cpow(IT_0005, 2);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_W*e_em*IT_0086*(IT_0068
      *IT_0087 + -IT_0003*(IT_0003*(IT_0005*IT_0068 + IT_0003*IT_0069) + 
      -IT_0069*IT_0088));
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = cpow(IT_0090, 2);
    const ccomplex_t IT_0092 = cpow((-2)*IT_0012 + IT_0076, 2);
    const ccomplex_t IT_0093 = cabs(-IT_0012 + 0.25*IT_0001*IT_0092);
    const ccomplex_t IT_0094 = cpow(IT_0093, 0.5);
    const ccomplex_t IT_0095 = ( ((m_A0
) >= (m_h + m_A0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0091*IT_0094
) : (0
));
    const ccomplex_t IT_0096 = e_em*IT_0007*IT_0067*(IT_0005*IT_0068 + 
      -IT_0003*IT_0069);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0097, 2);
    const ccomplex_t IT_0099 = pow(m_H0, 2);
    const ccomplex_t IT_0100 = 0.5*IT_0012 + 0.5*IT_0075 + (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = cpow(IT_0100, 2);
    const ccomplex_t IT_0102 = IT_0074*IT_0101;
    const ccomplex_t IT_0103 = cpow(IT_0012 + IT_0075 + -IT_0099, 2);
    const ccomplex_t IT_0104 = cabs(-IT_0075 + 0.25*IT_0001*IT_0103);
    const ccomplex_t IT_0105 = cpow(IT_0104, 0.5);
    const ccomplex_t IT_0106 = ( ((m_A0
) >= (m_Z + m_H0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0098*(IT_0012 + -IT_0102)*IT_0105
) : (0
));
    const ccomplex_t IT_0107 = -IT_0003;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_W*e_em*IT_0086*(IT_0069
      *IT_0087 + IT_0107*(IT_0068*IT_0088 + (IT_0003*IT_0068 + -IT_0005*IT_0069)
      *IT_0107));
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0109, 2);
    const ccomplex_t IT_0111 = cpow(IT_0012 + (-0.5)*IT_0099, 2);
    const ccomplex_t IT_0112 = cabs(-IT_0012 + IT_0001*IT_0111);
    const ccomplex_t IT_0113 = cpow(IT_0112, 0.5);
    const ccomplex_t IT_0114 = ( ((m_A0
) >= (m_A0 + m_H0
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0110*IT_0113
) : (0
));
    const ccomplex_t IT_0115 = e_em*IT_0007;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = cpow(IT_0116, 2);
    const ccomplex_t IT_0118 = pow(m_W, -2);
    const ccomplex_t IT_0119 = pow(m_W, 2);
    const ccomplex_t IT_0120 = pow(m_Hp, 2);
    const ccomplex_t IT_0121 = 0.5*IT_0012 + 0.5*IT_0119 + (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = cpow(IT_0121, 2);
    const ccomplex_t IT_0123 = IT_0118*IT_0122;
    const ccomplex_t IT_0124 = cpow(IT_0012 + IT_0119 + -IT_0120, 2);
    const ccomplex_t IT_0125 = cabs(-IT_0119 + 0.25*IT_0001*IT_0124);
    const ccomplex_t IT_0126 = cpow(IT_0125, 0.5);
    const ccomplex_t IT_0127 = ( ((m_A0
) >= (m_W + m_Hp
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0117*(IT_0012 + -IT_0123)*IT_0126
) : (0
));
    const ccomplex_t IT_0128 = e_em*m_mu*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = cpow(IT_0129, 2);
    const ccomplex_t IT_0131 = pow(m_mu, 2);
    const ccomplex_t IT_0132 = 0.5*IT_0012 + -IT_0131;
    const ccomplex_t IT_0133 = cabs(0.25*IT_0012 + -IT_0131);
    const ccomplex_t IT_0134 = cpow(IT_0133, 0.5);
    const ccomplex_t IT_0135 = ( ((m_A0
) >= (2*m_mu
)) ? (0.5*IT_0000*IT_0001
      *IT_0130*(IT_0131 + IT_0132)*IT_0134
) : (0
));
    const ccomplex_t IT_0136 = V_u1*U_Wm1;
    const ccomplex_t IT_0137 = IT_0003*IT_0136;
    const ccomplex_t IT_0138 = U_d1*V_Wp1;
    const ccomplex_t IT_0139 = IT_0005*IT_0138;
    const ccomplex_t IT_0140 = IT_0007*(IT_0137 + IT_0139);
    const ccomplex_t IT_0141 = 1.4142135623731*e_em*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0145 = IT_0003*IT_0144;
    const ccomplex_t IT_0146 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0147 = IT_0005*IT_0146;
    const ccomplex_t IT_0148 = IT_0007*(IT_0145 + IT_0147);
    const ccomplex_t IT_0149 = 1.4142135623731*e_em*IT_0148;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = pow(m_C_1, 2);
    const ccomplex_t IT_0153 = 0.5*IT_0012 + -IT_0152;
    const ccomplex_t IT_0154 = 2*IT_0153;
    const ccomplex_t IT_0155 = cabs(0.25*IT_0012 + -IT_0152);
    const ccomplex_t IT_0156 = cpow(IT_0155, 0.5);
    const ccomplex_t IT_0157 = ( ((m_A0
) >= (2*m_C_1
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0143)*(IT_0151*IT_0152 + (-0.5)*IT_0143*IT_0154) + conj
      (IT_0151)*(IT_0143*IT_0152 + (-0.5)*IT_0151*IT_0154))*IT_0156
) : (0
));
    const ccomplex_t IT_0158 = V_u2*U_Wm2;
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = U_d2*V_Wp2;
    const ccomplex_t IT_0161 = IT_0005*IT_0160;
    const ccomplex_t IT_0162 = IT_0007*(IT_0159 + IT_0161);
    const ccomplex_t IT_0163 = 1.4142135623731*e_em*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0167 = IT_0003*IT_0166;
    const ccomplex_t IT_0168 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0169 = IT_0005*IT_0168;
    const ccomplex_t IT_0170 = IT_0007*(IT_0167 + IT_0169);
    const ccomplex_t IT_0171 = 1.4142135623731*e_em*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = pow(m_C_2, 2);
    const ccomplex_t IT_0175 = 0.5*IT_0012 + -IT_0174;
    const ccomplex_t IT_0176 = 2*IT_0175;
    const ccomplex_t IT_0177 = cabs(0.25*IT_0012 + -IT_0174);
    const ccomplex_t IT_0178 = cpow(IT_0177, 0.5);
    const ccomplex_t IT_0179 = ( ((m_A0
) >= (2*m_C_2
)) ? ((-0.25)*IT_0000
      *IT_0001*(conj(IT_0165)*(IT_0173*IT_0174 + (-0.5)*IT_0165*IT_0176) + conj
      (IT_0173)*(IT_0165*IT_0174 + (-0.5)*IT_0173*IT_0176))*IT_0178
) : (0
));
    const ccomplex_t IT_0180 = V_u1*U_Wm2;
    const ccomplex_t IT_0181 = IT_0003*IT_0180;
    const ccomplex_t IT_0182 = U_d2*V_Wp1;
    const ccomplex_t IT_0183 = IT_0005*IT_0182;
    const ccomplex_t IT_0184 = IT_0007*(IT_0181 + IT_0183);
    const ccomplex_t IT_0185 = 1.4142135623731*e_em*IT_0184;
    const ccomplex_t IT_0186 = 0.5*IT_0185;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0189 = IT_0003*IT_0188;
    const ccomplex_t IT_0190 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0191 = IT_0005*IT_0190;
    const ccomplex_t IT_0192 = IT_0007*(IT_0189 + IT_0191);
    const ccomplex_t IT_0193 = 1.4142135623731*e_em*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = m_C_1*m_C_2;
    const ccomplex_t IT_0197 = cpow(IT_0012 + IT_0152 + -IT_0174, 2);
    const ccomplex_t IT_0198 = -IT_0152 + 0.25*IT_0001*(IT_0012 + -IT_0152 +
       IT_0174)*(IT_0012 + IT_0152 + -IT_0174) + 0.25*IT_0001*IT_0197;
    const ccomplex_t IT_0199 = 2*IT_0198;
    const ccomplex_t IT_0200 = cabs(-IT_0152 + 0.25*IT_0001*IT_0197);
    const ccomplex_t IT_0201 = cpow(IT_0200, 0.5);
    const ccomplex_t IT_0202 = ( ((m_A0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0187)*(IT_0195*IT_0196 + (-0.5)*IT_0187*IT_0199)
       + conj(IT_0195)*(IT_0187*IT_0196 + (-0.5)*IT_0195*IT_0199))*IT_0201
) : 
      (0
));
    const ccomplex_t IT_0203 = V_u2*U_Wm1;
    const ccomplex_t IT_0204 = IT_0003*IT_0203;
    const ccomplex_t IT_0205 = U_d1*V_Wp2;
    const ccomplex_t IT_0206 = IT_0005*IT_0205;
    const ccomplex_t IT_0207 = IT_0007*(IT_0204 + IT_0206);
    const ccomplex_t IT_0208 = 1.4142135623731*e_em*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0212 = IT_0003*IT_0211;
    const ccomplex_t IT_0213 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0214 = IT_0005*IT_0213;
    const ccomplex_t IT_0215 = IT_0007*(IT_0212 + IT_0214);
    const ccomplex_t IT_0216 = 1.4142135623731*e_em*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = ( ((m_A0
) >= (m_C_1 + m_C_2
)) ? ((-0.25)
      *IT_0000*IT_0001*IT_0201*(conj(IT_0210)*((-0.5)*IT_0199*IT_0210 + IT_0196
      *IT_0218) + (IT_0196*IT_0210 + (-0.5)*IT_0199*IT_0218)*conj(IT_0218))
) : 
      (0
));
    const ccomplex_t IT_0220 = IT_0003*IT_0067;
    const ccomplex_t IT_0221 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0222 = IT_0220*IT_0221;
    const ccomplex_t IT_0223 = IT_0005*IT_0067;
    const ccomplex_t IT_0224 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0225 = IT_0223*IT_0224;
    const ccomplex_t IT_0226 = IT_0003*IT_0007;
    const ccomplex_t IT_0227 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = IT_0005*IT_0007;
    const ccomplex_t IT_0230 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0231 = IT_0229*IT_0230;
    const ccomplex_t IT_0232 = -IT_0222 + IT_0225 + IT_0228 + -IT_0231;
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0232;
    const ccomplex_t IT_0235 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0236 = IT_0220*IT_0235;
    const ccomplex_t IT_0237 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0238 = IT_0223*IT_0237;
    const ccomplex_t IT_0239 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0240 = IT_0226*IT_0239;
    const ccomplex_t IT_0241 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0242 = IT_0229*IT_0241;
    const ccomplex_t IT_0243 = IT_0236 + -IT_0238 + -IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = (-0.5)*IT_0243;
    const ccomplex_t IT_0246 = pow(m_N_1, 2);
    const ccomplex_t IT_0247 = 0.5*IT_0012 + -IT_0246;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = ((conj(IT_0233) + -conj(IT_0234))*(IT_0244 + 
      -IT_0245) + (IT_0233 + -IT_0234)*(conj(IT_0244) + -conj(IT_0245)))*IT_0246
       + -(IT_0233*conj(IT_0233) + IT_0234*conj(IT_0234) + IT_0244*conj(IT_0244)
       + IT_0245*conj(IT_0245))*IT_0247 + -(conj(IT_0233)*IT_0234 + IT_0233*conj
      (IT_0234) + conj(IT_0244)*IT_0245 + IT_0244*conj(IT_0245))*IT_0248;
    const ccomplex_t IT_0250 = cabs(0.25*IT_0012 + -IT_0246);
    const ccomplex_t IT_0251 = cpow(IT_0250, 0.5);
    const ccomplex_t IT_0252 = ( ((m_A0
) >= (2*m_N_1
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0249*IT_0251
) : (0
));
    const ccomplex_t IT_0253 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0254 = IT_0220*IT_0253;
    const ccomplex_t IT_0255 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0256 = IT_0223*IT_0255;
    const ccomplex_t IT_0257 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0258 = IT_0226*IT_0257;
    const ccomplex_t IT_0259 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0260 = IT_0229*IT_0259;
    const ccomplex_t IT_0261 = -IT_0254 + IT_0256 + IT_0258 + -IT_0260;
    const ccomplex_t IT_0262 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0263 = IT_0220*IT_0262;
    const ccomplex_t IT_0264 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0265 = IT_0223*IT_0264;
    const ccomplex_t IT_0266 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0267 = IT_0226*IT_0266;
    const ccomplex_t IT_0268 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0269 = IT_0229*IT_0268;
    const ccomplex_t IT_0270 = IT_0263 + -IT_0265 + -IT_0267 + IT_0269;
    const ccomplex_t IT_0271 = pow(m_N_2, 2);
    const ccomplex_t IT_0272 = 0.5*IT_0012 + -IT_0271;
    const ccomplex_t IT_0273 = cabs(0.25*IT_0012 + -IT_0271);
    const ccomplex_t IT_0274 = cpow(IT_0273, 0.5);
    const ccomplex_t IT_0275 = ( ((m_A0
) >= (2*m_N_2
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0261)*IT_0270 + IT_0261*conj(IT_0270))*IT_0271 + -
      (IT_0261*conj(IT_0261) + IT_0270*conj(IT_0270))*IT_0272)*IT_0274
) : (0
));
    const ccomplex_t IT_0276 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0277 = IT_0220*IT_0276;
    const ccomplex_t IT_0278 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0279 = IT_0223*IT_0278;
    const ccomplex_t IT_0280 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0281 = IT_0226*IT_0280;
    const ccomplex_t IT_0282 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0283 = IT_0229*IT_0282;
    const ccomplex_t IT_0284 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0285 = IT_0220*IT_0284;
    const ccomplex_t IT_0286 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0287 = IT_0223*IT_0286;
    const ccomplex_t IT_0288 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0289 = IT_0226*IT_0288;
    const ccomplex_t IT_0290 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0291 = IT_0229*IT_0290;
    const ccomplex_t IT_0292 = (-0.5)*IT_0277 + 0.5*IT_0279 + 0.5*IT_0281 + (
      -0.5)*IT_0283 + (-0.5)*IT_0285 + 0.5*IT_0287 + 0.5*IT_0289 + (-0.5)*IT_0291;
    const ccomplex_t IT_0293 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0294 = IT_0220*IT_0293;
    const ccomplex_t IT_0295 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0296 = IT_0223*IT_0295;
    const ccomplex_t IT_0297 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0298 = IT_0226*IT_0297;
    const ccomplex_t IT_0299 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0300 = IT_0229*IT_0299;
    const ccomplex_t IT_0301 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0302 = IT_0220*IT_0301;
    const ccomplex_t IT_0303 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0304 = IT_0223*IT_0303;
    const ccomplex_t IT_0305 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0306 = IT_0226*IT_0305;
    const ccomplex_t IT_0307 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0308 = IT_0229*IT_0307;
    const ccomplex_t IT_0309 = 0.5*IT_0294 + (-0.5)*IT_0296 + (-0.5)*IT_0298 +
       0.5*IT_0300 + 0.5*IT_0302 + (-0.5)*IT_0304 + (-0.5)*IT_0306 + 0.5*IT_0308;
    const ccomplex_t IT_0310 = m_N_1*m_N_2;
    const ccomplex_t IT_0311 = cpow(IT_0012 + IT_0246 + -IT_0271, 2);
    const ccomplex_t IT_0312 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0271)*(IT_0012 + IT_0246 + -IT_0271) + 0.25*IT_0001*IT_0311;
    const ccomplex_t IT_0313 = 2*IT_0312;
    const ccomplex_t IT_0314 = cabs(-IT_0246 + 0.25*IT_0001*IT_0311);
    const ccomplex_t IT_0315 = cpow(IT_0314, 0.5);
    const ccomplex_t IT_0316 = ( ((m_A0
) >= (m_N_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0292)*(IT_0309*IT_0310 + (-0.5)*IT_0292*IT_0313)
       + conj(IT_0309)*(IT_0292*IT_0310 + (-0.5)*IT_0309*IT_0313))*IT_0315
) : 
      (0
));
    const ccomplex_t IT_0317 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0318 = IT_0220*IT_0317;
    const ccomplex_t IT_0319 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0320 = IT_0223*IT_0319;
    const ccomplex_t IT_0321 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0322 = IT_0226*IT_0321;
    const ccomplex_t IT_0323 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0324 = IT_0229*IT_0323;
    const ccomplex_t IT_0325 = -IT_0318 + IT_0320 + IT_0322 + -IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = 0.5*IT_0325;
    const ccomplex_t IT_0328 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0329 = IT_0220*IT_0328;
    const ccomplex_t IT_0330 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0331 = IT_0223*IT_0330;
    const ccomplex_t IT_0332 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0333 = IT_0226*IT_0332;
    const ccomplex_t IT_0334 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0335 = IT_0229*IT_0334;
    const ccomplex_t IT_0336 = IT_0329 + -IT_0331 + -IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = (-0.5)*IT_0336;
    const ccomplex_t IT_0338 = 0.5*IT_0336;
    const ccomplex_t IT_0339 = pow(m_N_3, 2);
    const ccomplex_t IT_0340 = 0.5*IT_0012 + -IT_0339;
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = ((conj(IT_0326) + -conj(IT_0327))*(IT_0337 + 
      -IT_0338) + (IT_0326 + -IT_0327)*(conj(IT_0337) + -conj(IT_0338)))*IT_0339
       + -(IT_0326*conj(IT_0326) + IT_0327*conj(IT_0327) + IT_0337*conj(IT_0337)
       + IT_0338*conj(IT_0338))*IT_0340 + -(conj(IT_0326)*IT_0327 + IT_0326*conj
      (IT_0327) + conj(IT_0337)*IT_0338 + IT_0337*conj(IT_0338))*IT_0341;
    const ccomplex_t IT_0343 = cabs(0.25*IT_0012 + -IT_0339);
    const ccomplex_t IT_0344 = cpow(IT_0343, 0.5);
    const ccomplex_t IT_0345 = ( ((m_A0
) >= (2*m_N_3
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0342*IT_0344
) : (0
));
    const ccomplex_t IT_0346 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0347 = IT_0220*IT_0346;
    const ccomplex_t IT_0348 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0349 = IT_0223*IT_0348;
    const ccomplex_t IT_0350 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0351 = IT_0226*IT_0350;
    const ccomplex_t IT_0352 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0353 = IT_0229*IT_0352;
    const ccomplex_t IT_0354 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0355 = IT_0220*IT_0354;
    const ccomplex_t IT_0356 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0357 = IT_0223*IT_0356;
    const ccomplex_t IT_0358 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0359 = IT_0226*IT_0358;
    const ccomplex_t IT_0360 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0361 = IT_0229*IT_0360;
    const ccomplex_t IT_0362 = (-0.5)*IT_0347 + 0.5*IT_0349 + 0.5*IT_0351 + (
      -0.5)*IT_0353 + (-0.5)*IT_0355 + 0.5*IT_0357 + 0.5*IT_0359 + (-0.5)*IT_0361;
    const ccomplex_t IT_0363 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0364 = IT_0220*IT_0363;
    const ccomplex_t IT_0365 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0366 = IT_0223*IT_0365;
    const ccomplex_t IT_0367 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0368 = IT_0226*IT_0367;
    const ccomplex_t IT_0369 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0370 = IT_0229*IT_0369;
    const ccomplex_t IT_0371 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0372 = IT_0220*IT_0371;
    const ccomplex_t IT_0373 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0374 = IT_0223*IT_0373;
    const ccomplex_t IT_0375 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0376 = IT_0226*IT_0375;
    const ccomplex_t IT_0377 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0378 = IT_0229*IT_0377;
    const ccomplex_t IT_0379 = 0.5*IT_0364 + (-0.5)*IT_0366 + (-0.5)*IT_0368 +
       0.5*IT_0370 + 0.5*IT_0372 + (-0.5)*IT_0374 + (-0.5)*IT_0376 + 0.5*IT_0378;
    const ccomplex_t IT_0380 = m_N_1*m_N_3;
    const ccomplex_t IT_0381 = cpow(IT_0012 + IT_0246 + -IT_0339, 2);
    const ccomplex_t IT_0382 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0339)*(IT_0012 + IT_0246 + -IT_0339) + 0.25*IT_0001*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0382;
    const ccomplex_t IT_0384 = cabs(-IT_0246 + 0.25*IT_0001*IT_0381);
    const ccomplex_t IT_0385 = cpow(IT_0384, 0.5);
    const ccomplex_t IT_0386 = ( ((m_A0
) >= (m_N_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0362)*(IT_0379*IT_0380 + (-0.5)*IT_0362*IT_0383)
       + conj(IT_0379)*(IT_0362*IT_0380 + (-0.5)*IT_0379*IT_0383))*IT_0385
) : 
      (0
));
    const ccomplex_t IT_0387 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0388 = IT_0220*IT_0387;
    const ccomplex_t IT_0389 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0390 = IT_0223*IT_0389;
    const ccomplex_t IT_0391 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0392 = IT_0226*IT_0391;
    const ccomplex_t IT_0393 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0394 = IT_0229*IT_0393;
    const ccomplex_t IT_0395 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0396 = IT_0220*IT_0395;
    const ccomplex_t IT_0397 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0398 = IT_0223*IT_0397;
    const ccomplex_t IT_0399 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0400 = IT_0226*IT_0399;
    const ccomplex_t IT_0401 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0402 = IT_0229*IT_0401;
    const ccomplex_t IT_0403 = (-0.5)*IT_0388 + 0.5*IT_0390 + 0.5*IT_0392 + (
      -0.5)*IT_0394 + (-0.5)*IT_0396 + 0.5*IT_0398 + 0.5*IT_0400 + (-0.5)*IT_0402;
    const ccomplex_t IT_0404 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0405 = IT_0220*IT_0404;
    const ccomplex_t IT_0406 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0407 = IT_0223*IT_0406;
    const ccomplex_t IT_0408 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0409 = IT_0226*IT_0408;
    const ccomplex_t IT_0410 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0411 = IT_0229*IT_0410;
    const ccomplex_t IT_0412 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0413 = IT_0220*IT_0412;
    const ccomplex_t IT_0414 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0415 = IT_0223*IT_0414;
    const ccomplex_t IT_0416 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0417 = IT_0226*IT_0416;
    const ccomplex_t IT_0418 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0419 = IT_0229*IT_0418;
    const ccomplex_t IT_0420 = 0.5*IT_0405 + (-0.5)*IT_0407 + (-0.5)*IT_0409 +
       0.5*IT_0411 + 0.5*IT_0413 + (-0.5)*IT_0415 + (-0.5)*IT_0417 + 0.5*IT_0419;
    const ccomplex_t IT_0421 = m_N_2*m_N_3;
    const ccomplex_t IT_0422 = cpow(IT_0012 + IT_0271 + -IT_0339, 2);
    const ccomplex_t IT_0423 = -IT_0271 + 0.25*IT_0001*(IT_0012 + -IT_0271 +
       IT_0339)*(IT_0012 + IT_0271 + -IT_0339) + 0.25*IT_0001*IT_0422;
    const ccomplex_t IT_0424 = 2*IT_0423;
    const ccomplex_t IT_0425 = cabs(-IT_0271 + 0.25*IT_0001*IT_0422);
    const ccomplex_t IT_0426 = cpow(IT_0425, 0.5);
    const ccomplex_t IT_0427 = ( ((m_A0
) >= (m_N_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0403)*(IT_0420*IT_0421 + (-0.5)*IT_0403*IT_0424)
       + conj(IT_0420)*(IT_0403*IT_0421 + (-0.5)*IT_0420*IT_0424))*IT_0426
) : 
      (0
));
    const ccomplex_t IT_0428 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0429 = IT_0220*IT_0428;
    const ccomplex_t IT_0430 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0431 = IT_0223*IT_0430;
    const ccomplex_t IT_0432 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0433 = IT_0226*IT_0432;
    const ccomplex_t IT_0434 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0435 = IT_0229*IT_0434;
    const ccomplex_t IT_0436 = -IT_0429 + IT_0431 + IT_0433 + -IT_0435;
    const ccomplex_t IT_0437 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0438 = IT_0220*IT_0437;
    const ccomplex_t IT_0439 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0440 = IT_0223*IT_0439;
    const ccomplex_t IT_0441 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0442 = IT_0226*IT_0441;
    const ccomplex_t IT_0443 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0444 = IT_0229*IT_0443;
    const ccomplex_t IT_0445 = IT_0438 + -IT_0440 + -IT_0442 + IT_0444;
    const ccomplex_t IT_0446 = pow(m_N_4, 2);
    const ccomplex_t IT_0447 = 0.5*IT_0012 + -IT_0446;
    const ccomplex_t IT_0448 = cabs(0.25*IT_0012 + -IT_0446);
    const ccomplex_t IT_0449 = cpow(IT_0448, 0.5);
    const ccomplex_t IT_0450 = ( ((m_A0
) >= (2*m_N_4
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0436)*IT_0445 + IT_0436*conj(IT_0445))*IT_0446 + -
      (IT_0436*conj(IT_0436) + IT_0445*conj(IT_0445))*IT_0447)*IT_0449
) : (0
));
    const ccomplex_t IT_0451 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0452 = IT_0220*IT_0451;
    const ccomplex_t IT_0453 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0454 = IT_0223*IT_0453;
    const ccomplex_t IT_0455 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0456 = IT_0226*IT_0455;
    const ccomplex_t IT_0457 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0458 = IT_0229*IT_0457;
    const ccomplex_t IT_0459 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0460 = IT_0220*IT_0459;
    const ccomplex_t IT_0461 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0462 = IT_0223*IT_0461;
    const ccomplex_t IT_0463 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0464 = IT_0226*IT_0463;
    const ccomplex_t IT_0465 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0466 = IT_0229*IT_0465;
    const ccomplex_t IT_0467 = (-0.5)*IT_0452 + 0.5*IT_0454 + 0.5*IT_0456 + (
      -0.5)*IT_0458 + (-0.5)*IT_0460 + 0.5*IT_0462 + 0.5*IT_0464 + (-0.5)*IT_0466;
    const ccomplex_t IT_0468 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0469 = IT_0220*IT_0468;
    const ccomplex_t IT_0470 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0471 = IT_0223*IT_0470;
    const ccomplex_t IT_0472 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0473 = IT_0226*IT_0472;
    const ccomplex_t IT_0474 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0475 = IT_0229*IT_0474;
    const ccomplex_t IT_0476 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0477 = IT_0220*IT_0476;
    const ccomplex_t IT_0478 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0479 = IT_0223*IT_0478;
    const ccomplex_t IT_0480 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0481 = IT_0226*IT_0480;
    const ccomplex_t IT_0482 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0483 = IT_0229*IT_0482;
    const ccomplex_t IT_0484 = 0.5*IT_0469 + (-0.5)*IT_0471 + (-0.5)*IT_0473 +
       0.5*IT_0475 + 0.5*IT_0477 + (-0.5)*IT_0479 + (-0.5)*IT_0481 + 0.5*IT_0483;
    const ccomplex_t IT_0485 = m_N_1*m_N_4;
    const ccomplex_t IT_0486 = cpow(IT_0012 + IT_0246 + -IT_0446, 2);
    const ccomplex_t IT_0487 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0446)*(IT_0012 + IT_0246 + -IT_0446) + 0.25*IT_0001*IT_0486;
    const ccomplex_t IT_0488 = 2*IT_0487;
    const ccomplex_t IT_0489 = cabs(-IT_0246 + 0.25*IT_0001*IT_0486);
    const ccomplex_t IT_0490 = cpow(IT_0489, 0.5);
    const ccomplex_t IT_0491 = ( ((m_A0
) >= (m_N_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0467)*(IT_0484*IT_0485 + (-0.5)*IT_0467*IT_0488)
       + conj(IT_0484)*(IT_0467*IT_0485 + (-0.5)*IT_0484*IT_0488))*IT_0490
) : 
      (0
));
    const ccomplex_t IT_0492 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0493 = IT_0220*IT_0492;
    const ccomplex_t IT_0494 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0495 = IT_0223*IT_0494;
    const ccomplex_t IT_0496 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0497 = IT_0226*IT_0496;
    const ccomplex_t IT_0498 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0499 = IT_0229*IT_0498;
    const ccomplex_t IT_0500 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0501 = IT_0220*IT_0500;
    const ccomplex_t IT_0502 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0503 = IT_0223*IT_0502;
    const ccomplex_t IT_0504 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0505 = IT_0226*IT_0504;
    const ccomplex_t IT_0506 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0507 = IT_0229*IT_0506;
    const ccomplex_t IT_0508 = (-0.5)*IT_0493 + 0.5*IT_0495 + 0.5*IT_0497 + (
      -0.5)*IT_0499 + (-0.5)*IT_0501 + 0.5*IT_0503 + 0.5*IT_0505 + (-0.5)*IT_0507;
    const ccomplex_t IT_0509 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0510 = IT_0220*IT_0509;
    const ccomplex_t IT_0511 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0512 = IT_0223*IT_0511;
    const ccomplex_t IT_0513 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0514 = IT_0226*IT_0513;
    const ccomplex_t IT_0515 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0516 = IT_0229*IT_0515;
    const ccomplex_t IT_0517 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0518 = IT_0220*IT_0517;
    const ccomplex_t IT_0519 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0520 = IT_0223*IT_0519;
    const ccomplex_t IT_0521 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0522 = IT_0226*IT_0521;
    const ccomplex_t IT_0523 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0524 = IT_0229*IT_0523;
    const ccomplex_t IT_0525 = 0.5*IT_0510 + (-0.5)*IT_0512 + (-0.5)*IT_0514 +
       0.5*IT_0516 + 0.5*IT_0518 + (-0.5)*IT_0520 + (-0.5)*IT_0522 + 0.5*IT_0524;
    const ccomplex_t IT_0526 = m_N_2*m_N_4;
    const ccomplex_t IT_0527 = cpow(IT_0012 + IT_0271 + -IT_0446, 2);
    const ccomplex_t IT_0528 = -IT_0271 + 0.25*IT_0001*(IT_0012 + -IT_0271 +
       IT_0446)*(IT_0012 + IT_0271 + -IT_0446) + 0.25*IT_0001*IT_0527;
    const ccomplex_t IT_0529 = 2*IT_0528;
    const ccomplex_t IT_0530 = cabs(-IT_0271 + 0.25*IT_0001*IT_0527);
    const ccomplex_t IT_0531 = cpow(IT_0530, 0.5);
    const ccomplex_t IT_0532 = ( ((m_A0
) >= (m_N_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0508)*(IT_0525*IT_0526 + (-0.5)*IT_0508*IT_0529)
       + conj(IT_0525)*(IT_0508*IT_0526 + (-0.5)*IT_0525*IT_0529))*IT_0531
) : 
      (0
));
    const ccomplex_t IT_0533 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0534 = IT_0220*IT_0533;
    const ccomplex_t IT_0535 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0536 = IT_0223*IT_0535;
    const ccomplex_t IT_0537 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0538 = IT_0226*IT_0537;
    const ccomplex_t IT_0539 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0540 = IT_0229*IT_0539;
    const ccomplex_t IT_0541 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0542 = IT_0220*IT_0541;
    const ccomplex_t IT_0543 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0544 = IT_0223*IT_0543;
    const ccomplex_t IT_0545 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0546 = IT_0226*IT_0545;
    const ccomplex_t IT_0547 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0548 = IT_0229*IT_0547;
    const ccomplex_t IT_0549 = (-0.5)*IT_0534 + 0.5*IT_0536 + 0.5*IT_0538 + (
      -0.5)*IT_0540 + (-0.5)*IT_0542 + 0.5*IT_0544 + 0.5*IT_0546 + (-0.5)*IT_0548;
    const ccomplex_t IT_0550 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0551 = IT_0220*IT_0550;
    const ccomplex_t IT_0552 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0553 = IT_0223*IT_0552;
    const ccomplex_t IT_0554 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0555 = IT_0226*IT_0554;
    const ccomplex_t IT_0556 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0557 = IT_0229*IT_0556;
    const ccomplex_t IT_0558 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0559 = IT_0220*IT_0558;
    const ccomplex_t IT_0560 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0561 = IT_0223*IT_0560;
    const ccomplex_t IT_0562 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0563 = IT_0226*IT_0562;
    const ccomplex_t IT_0564 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0565 = IT_0229*IT_0564;
    const ccomplex_t IT_0566 = 0.5*IT_0551 + (-0.5)*IT_0553 + (-0.5)*IT_0555 +
       0.5*IT_0557 + 0.5*IT_0559 + (-0.5)*IT_0561 + (-0.5)*IT_0563 + 0.5*IT_0565;
    const ccomplex_t IT_0567 = m_N_3*m_N_4;
    const ccomplex_t IT_0568 = cpow(IT_0012 + IT_0339 + -IT_0446, 2);
    const ccomplex_t IT_0569 = -IT_0339 + 0.25*IT_0001*(IT_0012 + -IT_0339 +
       IT_0446)*(IT_0012 + IT_0339 + -IT_0446) + 0.25*IT_0001*IT_0568;
    const ccomplex_t IT_0570 = 2*IT_0569;
    const ccomplex_t IT_0571 = cabs(-IT_0339 + 0.25*IT_0001*IT_0568);
    const ccomplex_t IT_0572 = cpow(IT_0571, 0.5);
    const ccomplex_t IT_0573 = ( ((m_A0
) >= (m_N_3 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0549)*(IT_0566*IT_0567 + (-0.5)*IT_0549*IT_0570)
       + conj(IT_0566)*(IT_0549*IT_0567 + (-0.5)*IT_0566*IT_0570))*IT_0572
) : 
      (0
));
    const ccomplex_t IT_0574 = e_em*m_tau*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0575 = (-0.5)*IT_0574;
    const ccomplex_t IT_0576 = cpow(IT_0575, 2);
    const ccomplex_t IT_0577 = pow(m_tau, 2);
    const ccomplex_t IT_0578 = 0.5*IT_0012 + -IT_0577;
    const ccomplex_t IT_0579 = cabs(0.25*IT_0012 + -IT_0577);
    const ccomplex_t IT_0580 = cpow(IT_0579, 0.5);
    const ccomplex_t IT_0581 = ( ((m_A0
) >= (2*m_tau
)) ? (0.5*IT_0000
      *IT_0001*IT_0576*(IT_0577 + IT_0578)*IT_0580
) : (0
));
    const ccomplex_t IT_0582 = m_b*e_em*mu_h*IT_0002*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0583 = IT_0007*IT_0582;
    const ccomplex_t IT_0584 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0585 = IT_0007*IT_0584;
    const ccomplex_t IT_0586 = cpow(V_tb, 2);
    const ccomplex_t IT_0587 = A_b*IT_0586*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0588 = IT_0005*IT_0587;
    const ccomplex_t IT_0589 = 1.4142135623731*IT_0588;
    const ccomplex_t IT_0590 = A_b*IT_0586*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0591 = IT_0005*IT_0590;
    const ccomplex_t IT_0592 = 1.4142135623731*IT_0591;
    const ccomplex_t IT_0593 = (-0.5)*IT_0583 + 0.5*IT_0585 + (-0.5)*IT_0589 +
       0.5*IT_0592;
    const ccomplex_t IT_0594 = pow(m_sb_1, 2);
    const ccomplex_t IT_0595 = cabs(0.25*IT_0012 + -IT_0594);
    const ccomplex_t IT_0596 = cpow(IT_0595, 0.5);
    const ccomplex_t IT_0597 = ( ((m_A0
) >= (2*m_sb_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0593*conj(IT_0593)*IT_0596
) : (0
));
    const ccomplex_t IT_0598 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0599 = IT_0007*IT_0598;
    const ccomplex_t IT_0600 = m_b*e_em*mu_h*IT_0002*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0601 = IT_0007*IT_0600;
    const ccomplex_t IT_0602 = A_b*IT_0586*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0603 = IT_0005*IT_0602;
    const ccomplex_t IT_0604 = 1.4142135623731*IT_0603;
    const ccomplex_t IT_0605 = A_b*IT_0586*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0606 = IT_0005*IT_0605;
    const ccomplex_t IT_0607 = 1.4142135623731*IT_0606;
    const ccomplex_t IT_0608 = (-0.5)*IT_0599 + 0.5*IT_0601 + (-0.5)*IT_0604 +
       0.5*IT_0607;
    const ccomplex_t IT_0609 = pow(m_sb_2, 2);
    const ccomplex_t IT_0610 = cabs(0.25*IT_0012 + -IT_0609);
    const ccomplex_t IT_0611 = cpow(IT_0610, 0.5);
    const ccomplex_t IT_0612 = ( ((m_A0
) >= (2*m_sb_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0608*conj(IT_0608)*IT_0611
) : (0
));
    const ccomplex_t IT_0613 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0614 = IT_0007*IT_0613;
    const ccomplex_t IT_0615 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0616 = IT_0007*IT_0615;
    const ccomplex_t IT_0617 = A_b*IT_0586*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0618 = IT_0005*IT_0617;
    const ccomplex_t IT_0619 = 1.4142135623731*IT_0618;
    const ccomplex_t IT_0620 = A_b*IT_0586*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0621 = IT_0005*IT_0620;
    const ccomplex_t IT_0622 = 1.4142135623731*IT_0621;
    const ccomplex_t IT_0623 = (-0.5)*IT_0614 + 0.5*IT_0616 + (-0.5)*IT_0619 +
       0.5*IT_0622;
    const ccomplex_t IT_0624 = cpow(IT_0012 + IT_0594 + -IT_0609, 2);
    const ccomplex_t IT_0625 = cabs(-IT_0594 + 0.25*IT_0001*IT_0624);
    const ccomplex_t IT_0626 = cpow(IT_0625, 0.5);
    const ccomplex_t IT_0627 = ( ((m_A0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0623*conj(IT_0623)*IT_0626
) : (0
));
    const ccomplex_t IT_0628 = m_b*e_em*mu_h*IT_0002*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0629 = IT_0007*IT_0628;
    const ccomplex_t IT_0630 = m_b*e_em*mu_h*IT_0002*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0631 = IT_0007*IT_0630;
    const ccomplex_t IT_0632 = A_b*IT_0586*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0633 = IT_0005*IT_0632;
    const ccomplex_t IT_0634 = 1.4142135623731*IT_0633;
    const ccomplex_t IT_0635 = A_b*IT_0586*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0636 = IT_0005*IT_0635;
    const ccomplex_t IT_0637 = 1.4142135623731*IT_0636;
    const ccomplex_t IT_0638 = 0.5*IT_0629 + (-0.5)*IT_0631 + 0.5*IT_0634 + (
      -0.5)*IT_0637;
    const ccomplex_t IT_0639 = ( ((m_A0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0626*IT_0638*conj(IT_0638)
) : (0
));
    const ccomplex_t IT_0640 = m_c*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0641 = 0.5*IT_0640;
    const ccomplex_t IT_0642 = cpow(IT_0641, 2);
    const ccomplex_t IT_0643 = pow(m_sc_L, 2);
    const ccomplex_t IT_0644 = pow(m_sc_R, 2);
    const ccomplex_t IT_0645 = cpow(IT_0012 + IT_0643 + -IT_0644, 2);
    const ccomplex_t IT_0646 = cabs(-IT_0643 + 0.25*IT_0001*IT_0645);
    const ccomplex_t IT_0647 = cpow(IT_0646, 0.5);
    const ccomplex_t IT_0648 = ( ((m_A0
) >= (m_sc_L + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0642*IT_0647
) : (0
));
    const ccomplex_t IT_0649 = (-0.5)*IT_0640;
    const ccomplex_t IT_0650 = cpow(IT_0649, 2);
    const ccomplex_t IT_0651 = ( ((m_A0
) >= (m_sc_L + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0647*IT_0650
) : (0
));
    const ccomplex_t IT_0652 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_10);
    const ccomplex_t IT_0653 = 0.5*IT_0652;
    const ccomplex_t IT_0654 = pow(m_sd_L, 2);
    const ccomplex_t IT_0655 = cpow(IT_0012 + -IT_0594 + IT_0654, 2);
    const ccomplex_t IT_0656 = cabs(-IT_0654 + 0.25*IT_0001*IT_0655);
    const ccomplex_t IT_0657 = cpow(IT_0656, 0.5);
    const ccomplex_t IT_0658 = ( ((m_A0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0653*conj(IT_0653)*IT_0657
) : (0
));
    const ccomplex_t IT_0659 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_10;
    const ccomplex_t IT_0660 = (-0.5)*IT_0659;
    const ccomplex_t IT_0661 = ( ((m_A0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0657*IT_0660*conj(IT_0660)
) : (0
));
    const ccomplex_t IT_0662 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_11;
    const ccomplex_t IT_0663 = (-0.5)*IT_0662;
    const ccomplex_t IT_0664 = cpow(IT_0012 + -IT_0609 + IT_0654, 2);
    const ccomplex_t IT_0665 = cabs(-IT_0654 + 0.25*IT_0001*IT_0664);
    const ccomplex_t IT_0666 = cpow(IT_0665, 0.5);
    const ccomplex_t IT_0667 = ( ((m_A0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0663*conj(IT_0663)*IT_0666
) : (0
));
    const ccomplex_t IT_0668 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_11);
    const ccomplex_t IT_0669 = 0.5*IT_0668;
    const ccomplex_t IT_0670 = ( ((m_A0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0666*IT_0669*conj(IT_0669)
) : (0
));
    const ccomplex_t IT_0671 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_00);
    const ccomplex_t IT_0672 = (-0.5)*IT_0671;
    const ccomplex_t IT_0673 = pow(m_sd_R, 2);
    const ccomplex_t IT_0674 = cpow(IT_0012 + -IT_0594 + IT_0673, 2);
    const ccomplex_t IT_0675 = cabs(-IT_0673 + 0.25*IT_0001*IT_0674);
    const ccomplex_t IT_0676 = cpow(IT_0675, 0.5);
    const ccomplex_t IT_0677 = ( ((m_A0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0672*conj(IT_0672)*IT_0676
) : (0
));
    const ccomplex_t IT_0678 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_00;
    const ccomplex_t IT_0679 = 0.5*IT_0678;
    const ccomplex_t IT_0680 = ( ((m_A0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0676*IT_0679*conj(IT_0679)
) : (0
));
    const ccomplex_t IT_0681 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_01);
    const ccomplex_t IT_0682 = (-0.5)*IT_0681;
    const ccomplex_t IT_0683 = cpow(IT_0012 + -IT_0609 + IT_0673, 2);
    const ccomplex_t IT_0684 = cabs(-IT_0673 + 0.25*IT_0001*IT_0683);
    const ccomplex_t IT_0685 = cpow(IT_0684, 0.5);
    const ccomplex_t IT_0686 = ( ((m_A0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0682*conj(IT_0682)*IT_0685
) : (0
));
    const ccomplex_t IT_0687 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_01;
    const ccomplex_t IT_0688 = 0.5*IT_0687;
    const ccomplex_t IT_0689 = ( ((m_A0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0685*IT_0688*conj(IT_0688)
) : (0
));
    const ccomplex_t IT_0690 = m_d*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0691 = IT_0007*IT_0690;
    const ccomplex_t IT_0692 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0693 = IT_0005*IT_0692;
    const ccomplex_t IT_0694 = 1.4142135623731*IT_0693;
    const ccomplex_t IT_0695 = 0.5*IT_0691 + 0.5*IT_0694;
    const ccomplex_t IT_0696 = cpow(IT_0695, 2);
    const ccomplex_t IT_0697 = cpow(IT_0012 + -IT_0654 + IT_0673, 2);
    const ccomplex_t IT_0698 = cabs(-IT_0673 + 0.25*IT_0001*IT_0697);
    const ccomplex_t IT_0699 = cpow(IT_0698, 0.5);
    const ccomplex_t IT_0700 = ( ((m_A0
) >= (m_sd_L + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0696*IT_0699
) : (0
));
    const ccomplex_t IT_0701 = (-0.5)*IT_0691 + (-0.5)*IT_0694;
    const ccomplex_t IT_0702 = cpow(IT_0701, 2);
    const ccomplex_t IT_0703 = ( ((m_A0
) >= (m_sd_L + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0699*IT_0702
) : (0
));
    const ccomplex_t IT_0704 = m_e*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0705 = 0.5*IT_0704;
    const ccomplex_t IT_0706 = cpow(IT_0705, 2);
    const ccomplex_t IT_0707 = pow(m_se_R, 2);
    const ccomplex_t IT_0708 = pow(m_se_L, 2);
    const ccomplex_t IT_0709 = cpow(IT_0012 + IT_0707 + -IT_0708, 2);
    const ccomplex_t IT_0710 = cabs(-IT_0707 + 0.25*IT_0001*IT_0709);
    const ccomplex_t IT_0711 = cpow(IT_0710, 0.5);
    const ccomplex_t IT_0712 = ( ((m_A0
) >= (m_se_L + m_se_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0706*IT_0711
) : (0
));
    const ccomplex_t IT_0713 = (-0.5)*IT_0704;
    const ccomplex_t IT_0714 = cpow(IT_0713, 2);
    const ccomplex_t IT_0715 = ( ((m_A0
) >= (m_se_L + m_se_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0711*IT_0714
) : (0
));
    const ccomplex_t IT_0716 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_10);
    const ccomplex_t IT_0717 = 0.5*IT_0716;
    const ccomplex_t IT_0718 = pow(m_ss_L, 2);
    const ccomplex_t IT_0719 = cpow(IT_0012 + -IT_0594 + IT_0718, 2);
    const ccomplex_t IT_0720 = cabs(-IT_0718 + 0.25*IT_0001*IT_0719);
    const ccomplex_t IT_0721 = cpow(IT_0720, 0.5);
    const ccomplex_t IT_0722 = ( ((m_A0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0717*conj(IT_0717)*IT_0721
) : (0
));
    const ccomplex_t IT_0723 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_10;
    const ccomplex_t IT_0724 = (-0.5)*IT_0723;
    const ccomplex_t IT_0725 = ( ((m_A0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0721*IT_0724*conj(IT_0724)
) : (0
));
    const ccomplex_t IT_0726 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_11);
    const ccomplex_t IT_0727 = 0.5*IT_0726;
    const ccomplex_t IT_0728 = cpow(IT_0012 + -IT_0609 + IT_0718, 2);
    const ccomplex_t IT_0729 = cabs(-IT_0718 + 0.25*IT_0001*IT_0728);
    const ccomplex_t IT_0730 = cpow(IT_0729, 0.5);
    const ccomplex_t IT_0731 = ( ((m_A0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0727*conj(IT_0727)*IT_0730
) : (0
));
    const ccomplex_t IT_0732 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_11;
    const ccomplex_t IT_0733 = (-0.5)*IT_0732;
    const ccomplex_t IT_0734 = ( ((m_A0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0730*IT_0733*conj(IT_0733)
) : (0
));
    const ccomplex_t IT_0735 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0005;
    const ccomplex_t IT_0736 = 0.5*IT_0735;
    const ccomplex_t IT_0737 = cpow(IT_0012 + IT_0673 + -IT_0718, 2);
    const ccomplex_t IT_0738 = cabs(-IT_0673 + 0.25*IT_0001*IT_0737);
    const ccomplex_t IT_0739 = cpow(IT_0738, 0.5);
    const ccomplex_t IT_0740 = ( ((m_A0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0736*conj(IT_0736)*IT_0739
) : (0
));
    const ccomplex_t IT_0741 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0005;
    const ccomplex_t IT_0742 = (-0.5)*IT_0741;
    const ccomplex_t IT_0743 = ( ((m_A0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0739*IT_0742*conj(IT_0742)
) : (0
));
    const ccomplex_t IT_0744 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_00);
    const ccomplex_t IT_0745 = (-0.5)*IT_0744;
    const ccomplex_t IT_0746 = pow(m_ss_R, 2);
    const ccomplex_t IT_0747 = cpow(IT_0012 + -IT_0594 + IT_0746, 2);
    const ccomplex_t IT_0748 = cabs(-IT_0746 + 0.25*IT_0001*IT_0747);
    const ccomplex_t IT_0749 = cpow(IT_0748, 0.5);
    const ccomplex_t IT_0750 = ( ((m_A0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0745*conj(IT_0745)*IT_0749
) : (0
));
    const ccomplex_t IT_0751 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_00;
    const ccomplex_t IT_0752 = 0.5*IT_0751;
    const ccomplex_t IT_0753 = ( ((m_A0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0749*IT_0752*conj(IT_0752)
) : (0
));
    const ccomplex_t IT_0754 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_01;
    const ccomplex_t IT_0755 = 0.5*IT_0754;
    const ccomplex_t IT_0756 = cpow(IT_0012 + -IT_0609 + IT_0746, 2);
    const ccomplex_t IT_0757 = cabs(-IT_0746 + 0.25*IT_0001*IT_0756);
    const ccomplex_t IT_0758 = cpow(IT_0757, 0.5);
    const ccomplex_t IT_0759 = ( ((m_A0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0755*conj(IT_0755)*IT_0758
) : (0
));
    const ccomplex_t IT_0760 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_01);
    const ccomplex_t IT_0761 = (-0.5)*IT_0760;
    const ccomplex_t IT_0762 = ( ((m_A0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0758*IT_0761*conj(IT_0761)
) : (0
));
    const ccomplex_t IT_0763 = 0.5*IT_0741;
    const ccomplex_t IT_0764 = cpow(IT_0012 + -IT_0654 + IT_0746, 2);
    const ccomplex_t IT_0765 = cabs(-IT_0746 + 0.25*IT_0001*IT_0764);
    const ccomplex_t IT_0766 = cpow(IT_0765, 0.5);
    const ccomplex_t IT_0767 = ( ((m_A0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0763*conj(IT_0763)*IT_0766
) : (0
));
    const ccomplex_t IT_0768 = (-0.5)*IT_0735;
    const ccomplex_t IT_0769 = ( ((m_A0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0766*IT_0768*conj(IT_0768)
) : (0
));
    const ccomplex_t IT_0770 = m_s*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0771 = IT_0007*IT_0770;
    const ccomplex_t IT_0772 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0773 = IT_0005*IT_0772;
    const ccomplex_t IT_0774 = 1.4142135623731*IT_0773;
    const ccomplex_t IT_0775 = 0.5*IT_0771 + 0.5*IT_0774;
    const ccomplex_t IT_0776 = cpow(IT_0775, 2);
    const ccomplex_t IT_0777 = cpow(IT_0012 + -IT_0718 + IT_0746, 2);
    const ccomplex_t IT_0778 = cabs(-IT_0746 + 0.25*IT_0001*IT_0777);
    const ccomplex_t IT_0779 = cpow(IT_0778, 0.5);
    const ccomplex_t IT_0780 = ( ((m_A0
) >= (m_ss_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0776*IT_0779
) : (0
));
    const ccomplex_t IT_0781 = (-0.5)*IT_0771 + (-0.5)*IT_0774;
    const ccomplex_t IT_0782 = cpow(IT_0781, 2);
    const ccomplex_t IT_0783 = ( ((m_A0
) >= (m_ss_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0779*IT_0782
) : (0
));
    const ccomplex_t IT_0784 = m_t*e_em*mu_h*IT_0002*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0785 = IT_0007*IT_0784;
    const ccomplex_t IT_0786 = m_t*e_em*mu_h*IT_0002*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0787 = IT_0007*IT_0786;
    const ccomplex_t IT_0788 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0789 = IT_0003*IT_0788;
    const ccomplex_t IT_0790 = 1.4142135623731*IT_0789;
    const ccomplex_t IT_0791 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0792 = IT_0003*IT_0791;
    const ccomplex_t IT_0793 = 1.4142135623731*IT_0792;
    const ccomplex_t IT_0794 = (-0.5)*IT_0785 + 0.5*IT_0787 + (-0.5)*IT_0790 +
       0.5*IT_0793;
    const ccomplex_t IT_0795 = pow(m_st_1, 2);
    const ccomplex_t IT_0796 = cabs(0.25*IT_0012 + -IT_0795);
    const ccomplex_t IT_0797 = cpow(IT_0796, 0.5);
    const ccomplex_t IT_0798 = ( ((m_A0
) >= (2*m_st_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0794*conj(IT_0794)*IT_0797
) : (0
));
    const ccomplex_t IT_0799 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0800 = IT_0007*IT_0799;
    const ccomplex_t IT_0801 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0802 = IT_0007*IT_0801;
    const ccomplex_t IT_0803 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0804 = IT_0003*IT_0803;
    const ccomplex_t IT_0805 = 1.4142135623731*IT_0804;
    const ccomplex_t IT_0806 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0807 = IT_0003*IT_0806;
    const ccomplex_t IT_0808 = 1.4142135623731*IT_0807;
    const ccomplex_t IT_0809 = (-0.5)*IT_0800 + 0.5*IT_0802 + (-0.5)*IT_0805 +
       0.5*IT_0808;
    const ccomplex_t IT_0810 = pow(m_st_2, 2);
    const ccomplex_t IT_0811 = cabs(0.25*IT_0012 + -IT_0810);
    const ccomplex_t IT_0812 = cpow(IT_0811, 0.5);
    const ccomplex_t IT_0813 = ( ((m_A0
) >= (2*m_st_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0809*conj(IT_0809)*IT_0812
) : (0
));
    const ccomplex_t IT_0814 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0815 = IT_0007*IT_0814;
    const ccomplex_t IT_0816 = m_t*e_em*mu_h*IT_0002*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0817 = IT_0007*IT_0816;
    const ccomplex_t IT_0818 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0819 = IT_0003*IT_0818;
    const ccomplex_t IT_0820 = 1.4142135623731*IT_0819;
    const ccomplex_t IT_0821 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0822 = IT_0003*IT_0821;
    const ccomplex_t IT_0823 = 1.4142135623731*IT_0822;
    const ccomplex_t IT_0824 = 0.5*IT_0815 + (-0.5)*IT_0817 + 0.5*IT_0820 + (
      -0.5)*IT_0823;
    const ccomplex_t IT_0825 = cpow(IT_0012 + IT_0795 + -IT_0810, 2);
    const ccomplex_t IT_0826 = cabs(-IT_0795 + 0.25*IT_0001*IT_0825);
    const ccomplex_t IT_0827 = cpow(IT_0826, 0.5);
    const ccomplex_t IT_0828 = ( ((m_A0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0824*conj(IT_0824)*IT_0827
) : (0
));
    const ccomplex_t IT_0829 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0830 = IT_0007*IT_0829;
    const ccomplex_t IT_0831 = m_t*e_em*mu_h*IT_0002*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0832 = IT_0007*IT_0831;
    const ccomplex_t IT_0833 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0834 = IT_0003*IT_0833;
    const ccomplex_t IT_0835 = 1.4142135623731*IT_0834;
    const ccomplex_t IT_0836 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0837 = IT_0003*IT_0836;
    const ccomplex_t IT_0838 = 1.4142135623731*IT_0837;
    const ccomplex_t IT_0839 = (-0.5)*IT_0830 + 0.5*IT_0832 + (-0.5)*IT_0835 +
       0.5*IT_0838;
    const ccomplex_t IT_0840 = ( ((m_A0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0827*IT_0839*conj(IT_0839)
) : (0
));
    const ccomplex_t IT_0841 = m_u*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0842 = 0.5*IT_0841;
    const ccomplex_t IT_0843 = cpow(IT_0842, 2);
    const ccomplex_t IT_0844 = pow(m_su_L, 2);
    const ccomplex_t IT_0845 = pow(m_su_R, 2);
    const ccomplex_t IT_0846 = cpow(IT_0012 + IT_0844 + -IT_0845, 2);
    const ccomplex_t IT_0847 = cabs(-IT_0844 + 0.25*IT_0001*IT_0846);
    const ccomplex_t IT_0848 = cpow(IT_0847, 0.5);
    const ccomplex_t IT_0849 = ( ((m_A0
) >= (m_su_L + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0843*IT_0848
) : (0
));
    const ccomplex_t IT_0850 = (-0.5)*IT_0841;
    const ccomplex_t IT_0851 = cpow(IT_0850, 2);
    const ccomplex_t IT_0852 = ( ((m_A0
) >= (m_su_L + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0848*IT_0851
) : (0
));
    const ccomplex_t IT_0853 = e_em*m_mu*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0854 = 0.5*IT_0853;
    const ccomplex_t IT_0855 = cpow(IT_0854, 2);
    const ccomplex_t IT_0856 = pow(m_smu_R, 2);
    const ccomplex_t IT_0857 = pow(m_smu_L, 2);
    const ccomplex_t IT_0858 = cpow(IT_0012 + IT_0856 + -IT_0857, 2);
    const ccomplex_t IT_0859 = cabs(-IT_0856 + 0.25*IT_0001*IT_0858);
    const ccomplex_t IT_0860 = cpow(IT_0859, 0.5);
    const ccomplex_t IT_0861 = ( ((m_A0
) >= (m_smu_L + m_smu_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0855*IT_0860
) : (0
));
    const ccomplex_t IT_0862 = (-0.5)*IT_0853;
    const ccomplex_t IT_0863 = cpow(IT_0862, 2);
    const ccomplex_t IT_0864 = ( ((m_A0
) >= (m_smu_L + m_smu_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0860*IT_0863
) : (0
));
    const ccomplex_t IT_0865 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0866 = IT_0007*IT_0865;
    const ccomplex_t IT_0867 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0868 = IT_0007*IT_0867;
    const ccomplex_t IT_0869 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0870 = IT_0005*IT_0869;
    const ccomplex_t IT_0871 = 1.4142135623731*IT_0870;
    const ccomplex_t IT_0872 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0873 = IT_0005*IT_0872;
    const ccomplex_t IT_0874 = 1.4142135623731*IT_0873;
    const ccomplex_t IT_0875 = (-0.5)*IT_0866 + 0.5*IT_0868 + (-0.5)*IT_0871 +
       0.5*IT_0874;
    const ccomplex_t IT_0876 = pow(m_stau_1, 2);
    const ccomplex_t IT_0877 = cabs(0.25*IT_0012 + -IT_0876);
    const ccomplex_t IT_0878 = cpow(IT_0877, 0.5);
    const ccomplex_t IT_0879 = ( ((m_A0
) >= (2*m_stau_1
)) ? (0.125*IT_0000
      *IT_0001*IT_0875*conj(IT_0875)*IT_0878
) : (0
));
    const ccomplex_t IT_0880 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0881 = IT_0007*IT_0880;
    const ccomplex_t IT_0882 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0883 = IT_0007*IT_0882;
    const ccomplex_t IT_0884 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0885 = IT_0005*IT_0884;
    const ccomplex_t IT_0886 = 1.4142135623731*IT_0885;
    const ccomplex_t IT_0887 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0888 = IT_0005*IT_0887;
    const ccomplex_t IT_0889 = 1.4142135623731*IT_0888;
    const ccomplex_t IT_0890 = (-0.5)*IT_0881 + 0.5*IT_0883 + (-0.5)*IT_0886 +
       0.5*IT_0889;
    const ccomplex_t IT_0891 = pow(m_stau_2, 2);
    const ccomplex_t IT_0892 = cabs(0.25*IT_0012 + -IT_0891);
    const ccomplex_t IT_0893 = cpow(IT_0892, 0.5);
    const ccomplex_t IT_0894 = ( ((m_A0
) >= (2*m_stau_2
)) ? (0.125*IT_0000
      *IT_0001*IT_0890*conj(IT_0890)*IT_0893
) : (0
));
    const ccomplex_t IT_0895 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0896 = IT_0007*IT_0895;
    const ccomplex_t IT_0897 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0898 = IT_0007*IT_0897;
    const ccomplex_t IT_0899 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0900 = IT_0005*IT_0899;
    const ccomplex_t IT_0901 = 1.4142135623731*IT_0900;
    const ccomplex_t IT_0902 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0903 = IT_0005*IT_0902;
    const ccomplex_t IT_0904 = 1.4142135623731*IT_0903;
    const ccomplex_t IT_0905 = (-0.5)*IT_0896 + 0.5*IT_0898 + (-0.5)*IT_0901 +
       0.5*IT_0904;
    const ccomplex_t IT_0906 = cpow(IT_0012 + IT_0876 + -IT_0891, 2);
    const ccomplex_t IT_0907 = cabs(-IT_0876 + 0.25*IT_0001*IT_0906);
    const ccomplex_t IT_0908 = cpow(IT_0907, 0.5);
    const ccomplex_t IT_0909 = ( ((m_A0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0905*conj(IT_0905)*IT_0908
) : (0
));
    const ccomplex_t IT_0910 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0911 = IT_0007*IT_0910;
    const ccomplex_t IT_0912 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0913 = IT_0007*IT_0912;
    const ccomplex_t IT_0914 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0915 = IT_0005*IT_0914;
    const ccomplex_t IT_0916 = 1.4142135623731*IT_0915;
    const ccomplex_t IT_0917 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0918 = IT_0005*IT_0917;
    const ccomplex_t IT_0919 = 1.4142135623731*IT_0918;
    const ccomplex_t IT_0920 = 0.5*IT_0911 + (-0.5)*IT_0913 + 0.5*IT_0916 + (
      -0.5)*IT_0919;
    const ccomplex_t IT_0921 = ( ((m_A0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0908*IT_0920*conj(IT_0920)
) : (0
));
    const ccomplex_t IT_0922 = IT_0016 + IT_0025 + IT_0033 + IT_0041 + IT_0049
       + IT_0057 + IT_0065 + IT_0083 + IT_0095 + IT_0106 + IT_0114 + 2*IT_0127 +
       IT_0135 + IT_0157 + IT_0179 + IT_0202 + IT_0219 + IT_0252 + IT_0275 +
       IT_0316 + IT_0345 + IT_0386 + IT_0427 + IT_0450 + IT_0491 + IT_0532 +
       IT_0573 + IT_0581 + IT_0597 + IT_0612 + IT_0627 + IT_0639 + IT_0648 +
       IT_0651 + IT_0658 + IT_0661 + IT_0667 + IT_0670 + IT_0677 + IT_0680 +
       IT_0686 + IT_0689 + IT_0700 + IT_0703 + IT_0712 + IT_0715 + IT_0722 +
       IT_0725 + IT_0731 + IT_0734 + IT_0740 + IT_0743 + IT_0750 + IT_0753 +
       IT_0759 + IT_0762 + IT_0767 + IT_0769 + IT_0780 + IT_0783 + IT_0798 +
       IT_0813 + IT_0828 + IT_0840 + IT_0849 + IT_0852 + IT_0861 + IT_0864 +
       IT_0879 + IT_0894 + IT_0909 + IT_0921;
    return create_ccomplex_return(IT_0922);
}

