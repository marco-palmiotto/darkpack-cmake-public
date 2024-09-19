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
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0232;
    const ccomplex_t IT_0235 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0236 = IT_0220*IT_0235;
    const ccomplex_t IT_0237 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0238 = IT_0223*IT_0237;
    const ccomplex_t IT_0239 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0240 = IT_0226*IT_0239;
    const ccomplex_t IT_0241 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0242 = IT_0229*IT_0241;
    const ccomplex_t IT_0243 = IT_0236 + -IT_0238 + -IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0243;
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
    const ccomplex_t IT_0262 = 0.5*IT_0261;
    const ccomplex_t IT_0263 = (-0.5)*IT_0261;
    const ccomplex_t IT_0264 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0265 = IT_0220*IT_0264;
    const ccomplex_t IT_0266 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0267 = IT_0223*IT_0266;
    const ccomplex_t IT_0268 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0269 = IT_0226*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0271 = IT_0229*IT_0270;
    const ccomplex_t IT_0272 = IT_0265 + -IT_0267 + -IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = 0.5*IT_0272;
    const ccomplex_t IT_0274 = (-0.5)*IT_0272;
    const ccomplex_t IT_0275 = pow(m_N_2, 2);
    const ccomplex_t IT_0276 = 0.5*IT_0012 + -IT_0275;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = ((conj(IT_0262) + -conj(IT_0263))*(IT_0273 + 
      -IT_0274) + (IT_0262 + -IT_0263)*(conj(IT_0273) + -conj(IT_0274)))*IT_0275
       + -(IT_0262*conj(IT_0262) + IT_0263*conj(IT_0263) + IT_0273*conj(IT_0273)
       + IT_0274*conj(IT_0274))*IT_0276 + -(conj(IT_0262)*IT_0263 + IT_0262*conj
      (IT_0263) + conj(IT_0273)*IT_0274 + IT_0273*conj(IT_0274))*IT_0277;
    const ccomplex_t IT_0279 = cabs(0.25*IT_0012 + -IT_0275);
    const ccomplex_t IT_0280 = cpow(IT_0279, 0.5);
    const ccomplex_t IT_0281 = ( ((m_A0
) >= (2*m_N_2
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0278*IT_0280
) : (0
));
    const ccomplex_t IT_0282 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0283 = IT_0220*IT_0282;
    const ccomplex_t IT_0284 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0285 = IT_0223*IT_0284;
    const ccomplex_t IT_0286 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0287 = IT_0226*IT_0286;
    const ccomplex_t IT_0288 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0289 = IT_0229*IT_0288;
    const ccomplex_t IT_0290 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0291 = IT_0220*IT_0290;
    const ccomplex_t IT_0292 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0293 = IT_0223*IT_0292;
    const ccomplex_t IT_0294 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0295 = IT_0226*IT_0294;
    const ccomplex_t IT_0296 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0297 = IT_0229*IT_0296;
    const ccomplex_t IT_0298 = (-0.5)*IT_0283 + 0.5*IT_0285 + 0.5*IT_0287 + (
      -0.5)*IT_0289 + (-0.5)*IT_0291 + 0.5*IT_0293 + 0.5*IT_0295 + (-0.5)*IT_0297;
    const ccomplex_t IT_0299 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0300 = IT_0220*IT_0299;
    const ccomplex_t IT_0301 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0302 = IT_0223*IT_0301;
    const ccomplex_t IT_0303 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0304 = IT_0226*IT_0303;
    const ccomplex_t IT_0305 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0306 = IT_0229*IT_0305;
    const ccomplex_t IT_0307 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0308 = IT_0220*IT_0307;
    const ccomplex_t IT_0309 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0310 = IT_0223*IT_0309;
    const ccomplex_t IT_0311 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0312 = IT_0226*IT_0311;
    const ccomplex_t IT_0313 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0314 = IT_0229*IT_0313;
    const ccomplex_t IT_0315 = 0.5*IT_0300 + (-0.5)*IT_0302 + (-0.5)*IT_0304 +
       0.5*IT_0306 + 0.5*IT_0308 + (-0.5)*IT_0310 + (-0.5)*IT_0312 + 0.5*IT_0314;
    const ccomplex_t IT_0316 = m_N_1*m_N_2;
    const ccomplex_t IT_0317 = cpow(IT_0012 + IT_0246 + -IT_0275, 2);
    const ccomplex_t IT_0318 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0275)*(IT_0012 + IT_0246 + -IT_0275) + 0.25*IT_0001*IT_0317;
    const ccomplex_t IT_0319 = 2*IT_0318;
    const ccomplex_t IT_0320 = cabs(-IT_0246 + 0.25*IT_0001*IT_0317);
    const ccomplex_t IT_0321 = cpow(IT_0320, 0.5);
    const ccomplex_t IT_0322 = ( ((m_A0
) >= (m_N_1 + m_N_2
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0298)*(IT_0315*IT_0316 + (-0.5)*IT_0298*IT_0319)
       + conj(IT_0315)*(IT_0298*IT_0316 + (-0.5)*IT_0315*IT_0319))*IT_0321
) : 
      (0
));
    const ccomplex_t IT_0323 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0324 = IT_0220*IT_0323;
    const ccomplex_t IT_0325 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0326 = IT_0223*IT_0325;
    const ccomplex_t IT_0327 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0328 = IT_0226*IT_0327;
    const ccomplex_t IT_0329 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0330 = IT_0229*IT_0329;
    const ccomplex_t IT_0331 = -IT_0324 + IT_0326 + IT_0328 + -IT_0330;
    const ccomplex_t IT_0332 = 0.5*IT_0331;
    const ccomplex_t IT_0333 = (-0.5)*IT_0331;
    const ccomplex_t IT_0334 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0335 = IT_0220*IT_0334;
    const ccomplex_t IT_0336 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0337 = IT_0223*IT_0336;
    const ccomplex_t IT_0338 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0339 = IT_0226*IT_0338;
    const ccomplex_t IT_0340 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0341 = IT_0229*IT_0340;
    const ccomplex_t IT_0342 = IT_0335 + -IT_0337 + -IT_0339 + IT_0341;
    const ccomplex_t IT_0343 = 0.5*IT_0342;
    const ccomplex_t IT_0344 = (-0.5)*IT_0342;
    const ccomplex_t IT_0345 = pow(m_N_3, 2);
    const ccomplex_t IT_0346 = 0.5*IT_0012 + -IT_0345;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = ((conj(IT_0332) + -conj(IT_0333))*(IT_0343 + 
      -IT_0344) + (IT_0332 + -IT_0333)*(conj(IT_0343) + -conj(IT_0344)))*IT_0345
       + -(IT_0332*conj(IT_0332) + IT_0333*conj(IT_0333) + IT_0343*conj(IT_0343)
       + IT_0344*conj(IT_0344))*IT_0346 + -(conj(IT_0332)*IT_0333 + IT_0332*conj
      (IT_0333) + conj(IT_0343)*IT_0344 + IT_0343*conj(IT_0344))*IT_0347;
    const ccomplex_t IT_0349 = cabs(0.25*IT_0012 + -IT_0345);
    const ccomplex_t IT_0350 = cpow(IT_0349, 0.5);
    const ccomplex_t IT_0351 = ( ((m_A0
) >= (2*m_N_3
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0348*IT_0350
) : (0
));
    const ccomplex_t IT_0352 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0353 = IT_0220*IT_0352;
    const ccomplex_t IT_0354 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0355 = IT_0223*IT_0354;
    const ccomplex_t IT_0356 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0357 = IT_0226*IT_0356;
    const ccomplex_t IT_0358 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0359 = IT_0229*IT_0358;
    const ccomplex_t IT_0360 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0361 = IT_0220*IT_0360;
    const ccomplex_t IT_0362 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0363 = IT_0223*IT_0362;
    const ccomplex_t IT_0364 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0365 = IT_0226*IT_0364;
    const ccomplex_t IT_0366 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0367 = IT_0229*IT_0366;
    const ccomplex_t IT_0368 = (-0.5)*IT_0353 + 0.5*IT_0355 + 0.5*IT_0357 + (
      -0.5)*IT_0359 + (-0.5)*IT_0361 + 0.5*IT_0363 + 0.5*IT_0365 + (-0.5)*IT_0367;
    const ccomplex_t IT_0369 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0370 = IT_0220*IT_0369;
    const ccomplex_t IT_0371 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0372 = IT_0223*IT_0371;
    const ccomplex_t IT_0373 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0374 = IT_0226*IT_0373;
    const ccomplex_t IT_0375 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0376 = IT_0229*IT_0375;
    const ccomplex_t IT_0377 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0378 = IT_0220*IT_0377;
    const ccomplex_t IT_0379 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0380 = IT_0223*IT_0379;
    const ccomplex_t IT_0381 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0382 = IT_0226*IT_0381;
    const ccomplex_t IT_0383 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0384 = IT_0229*IT_0383;
    const ccomplex_t IT_0385 = 0.5*IT_0370 + (-0.5)*IT_0372 + (-0.5)*IT_0374 +
       0.5*IT_0376 + 0.5*IT_0378 + (-0.5)*IT_0380 + (-0.5)*IT_0382 + 0.5*IT_0384;
    const ccomplex_t IT_0386 = m_N_1*m_N_3;
    const ccomplex_t IT_0387 = cpow(IT_0012 + IT_0246 + -IT_0345, 2);
    const ccomplex_t IT_0388 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0345)*(IT_0012 + IT_0246 + -IT_0345) + 0.25*IT_0001*IT_0387;
    const ccomplex_t IT_0389 = 2*IT_0388;
    const ccomplex_t IT_0390 = cabs(-IT_0246 + 0.25*IT_0001*IT_0387);
    const ccomplex_t IT_0391 = cpow(IT_0390, 0.5);
    const ccomplex_t IT_0392 = ( ((m_A0
) >= (m_N_1 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0368)*(IT_0385*IT_0386 + (-0.5)*IT_0368*IT_0389)
       + conj(IT_0385)*(IT_0368*IT_0386 + (-0.5)*IT_0385*IT_0389))*IT_0391
) : 
      (0
));
    const ccomplex_t IT_0393 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0394 = IT_0220*IT_0393;
    const ccomplex_t IT_0395 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0396 = IT_0223*IT_0395;
    const ccomplex_t IT_0397 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0398 = IT_0226*IT_0397;
    const ccomplex_t IT_0399 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0400 = IT_0229*IT_0399;
    const ccomplex_t IT_0401 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0402 = IT_0220*IT_0401;
    const ccomplex_t IT_0403 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0404 = IT_0223*IT_0403;
    const ccomplex_t IT_0405 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0406 = IT_0226*IT_0405;
    const ccomplex_t IT_0407 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0408 = IT_0229*IT_0407;
    const ccomplex_t IT_0409 = (-0.5)*IT_0394 + 0.5*IT_0396 + 0.5*IT_0398 + (
      -0.5)*IT_0400 + (-0.5)*IT_0402 + 0.5*IT_0404 + 0.5*IT_0406 + (-0.5)*IT_0408;
    const ccomplex_t IT_0410 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0411 = IT_0220*IT_0410;
    const ccomplex_t IT_0412 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0413 = IT_0223*IT_0412;
    const ccomplex_t IT_0414 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0415 = IT_0226*IT_0414;
    const ccomplex_t IT_0416 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0417 = IT_0229*IT_0416;
    const ccomplex_t IT_0418 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0419 = IT_0220*IT_0418;
    const ccomplex_t IT_0420 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0421 = IT_0223*IT_0420;
    const ccomplex_t IT_0422 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0423 = IT_0226*IT_0422;
    const ccomplex_t IT_0424 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0425 = IT_0229*IT_0424;
    const ccomplex_t IT_0426 = 0.5*IT_0411 + (-0.5)*IT_0413 + (-0.5)*IT_0415 +
       0.5*IT_0417 + 0.5*IT_0419 + (-0.5)*IT_0421 + (-0.5)*IT_0423 + 0.5*IT_0425;
    const ccomplex_t IT_0427 = m_N_2*m_N_3;
    const ccomplex_t IT_0428 = cpow(IT_0012 + IT_0275 + -IT_0345, 2);
    const ccomplex_t IT_0429 = -IT_0275 + 0.25*IT_0001*(IT_0012 + -IT_0275 +
       IT_0345)*(IT_0012 + IT_0275 + -IT_0345) + 0.25*IT_0001*IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0429;
    const ccomplex_t IT_0431 = cabs(-IT_0275 + 0.25*IT_0001*IT_0428);
    const ccomplex_t IT_0432 = cpow(IT_0431, 0.5);
    const ccomplex_t IT_0433 = ( ((m_A0
) >= (m_N_2 + m_N_3
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0409)*(IT_0426*IT_0427 + (-0.5)*IT_0409*IT_0430)
       + conj(IT_0426)*(IT_0409*IT_0427 + (-0.5)*IT_0426*IT_0430))*IT_0432
) : 
      (0
));
    const ccomplex_t IT_0434 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0435 = IT_0220*IT_0434;
    const ccomplex_t IT_0436 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0437 = IT_0223*IT_0436;
    const ccomplex_t IT_0438 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0439 = IT_0226*IT_0438;
    const ccomplex_t IT_0440 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0441 = IT_0229*IT_0440;
    const ccomplex_t IT_0442 = -IT_0435 + IT_0437 + IT_0439 + -IT_0441;
    const ccomplex_t IT_0443 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0444 = IT_0220*IT_0443;
    const ccomplex_t IT_0445 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0446 = IT_0223*IT_0445;
    const ccomplex_t IT_0447 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0448 = IT_0226*IT_0447;
    const ccomplex_t IT_0449 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0450 = IT_0229*IT_0449;
    const ccomplex_t IT_0451 = IT_0444 + -IT_0446 + -IT_0448 + IT_0450;
    const ccomplex_t IT_0452 = pow(m_N_4, 2);
    const ccomplex_t IT_0453 = 0.5*IT_0012 + -IT_0452;
    const ccomplex_t IT_0454 = cabs(0.25*IT_0012 + -IT_0452);
    const ccomplex_t IT_0455 = cpow(IT_0454, 0.5);
    const ccomplex_t IT_0456 = ( ((m_A0
) >= (2*m_N_4
)) ? ((-0.125)*IT_0000
      *IT_0001*((conj(IT_0442)*IT_0451 + IT_0442*conj(IT_0451))*IT_0452 + -
      (IT_0442*conj(IT_0442) + IT_0451*conj(IT_0451))*IT_0453)*IT_0455
) : (0
));
    const ccomplex_t IT_0457 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0458 = IT_0220*IT_0457;
    const ccomplex_t IT_0459 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0460 = IT_0223*IT_0459;
    const ccomplex_t IT_0461 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0462 = IT_0226*IT_0461;
    const ccomplex_t IT_0463 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0464 = IT_0229*IT_0463;
    const ccomplex_t IT_0465 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0466 = IT_0220*IT_0465;
    const ccomplex_t IT_0467 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0468 = IT_0223*IT_0467;
    const ccomplex_t IT_0469 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0470 = IT_0226*IT_0469;
    const ccomplex_t IT_0471 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0472 = IT_0229*IT_0471;
    const ccomplex_t IT_0473 = (-0.5)*IT_0458 + 0.5*IT_0460 + 0.5*IT_0462 + (
      -0.5)*IT_0464 + (-0.5)*IT_0466 + 0.5*IT_0468 + 0.5*IT_0470 + (-0.5)*IT_0472;
    const ccomplex_t IT_0474 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0475 = IT_0220*IT_0474;
    const ccomplex_t IT_0476 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0477 = IT_0223*IT_0476;
    const ccomplex_t IT_0478 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0479 = IT_0226*IT_0478;
    const ccomplex_t IT_0480 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0481 = IT_0229*IT_0480;
    const ccomplex_t IT_0482 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0483 = IT_0220*IT_0482;
    const ccomplex_t IT_0484 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0485 = IT_0223*IT_0484;
    const ccomplex_t IT_0486 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0487 = IT_0226*IT_0486;
    const ccomplex_t IT_0488 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0489 = IT_0229*IT_0488;
    const ccomplex_t IT_0490 = 0.5*IT_0475 + (-0.5)*IT_0477 + (-0.5)*IT_0479 +
       0.5*IT_0481 + 0.5*IT_0483 + (-0.5)*IT_0485 + (-0.5)*IT_0487 + 0.5*IT_0489;
    const ccomplex_t IT_0491 = m_N_1*m_N_4;
    const ccomplex_t IT_0492 = cpow(IT_0012 + IT_0246 + -IT_0452, 2);
    const ccomplex_t IT_0493 = -IT_0246 + 0.25*IT_0001*(IT_0012 + -IT_0246 +
       IT_0452)*(IT_0012 + IT_0246 + -IT_0452) + 0.25*IT_0001*IT_0492;
    const ccomplex_t IT_0494 = 2*IT_0493;
    const ccomplex_t IT_0495 = cabs(-IT_0246 + 0.25*IT_0001*IT_0492);
    const ccomplex_t IT_0496 = cpow(IT_0495, 0.5);
    const ccomplex_t IT_0497 = ( ((m_A0
) >= (m_N_1 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0473)*(IT_0490*IT_0491 + (-0.5)*IT_0473*IT_0494)
       + conj(IT_0490)*(IT_0473*IT_0491 + (-0.5)*IT_0490*IT_0494))*IT_0496
) : 
      (0
));
    const ccomplex_t IT_0498 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0499 = IT_0220*IT_0498;
    const ccomplex_t IT_0500 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0501 = IT_0223*IT_0500;
    const ccomplex_t IT_0502 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0503 = IT_0226*IT_0502;
    const ccomplex_t IT_0504 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0505 = IT_0229*IT_0504;
    const ccomplex_t IT_0506 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0507 = IT_0220*IT_0506;
    const ccomplex_t IT_0508 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0509 = IT_0223*IT_0508;
    const ccomplex_t IT_0510 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0511 = IT_0226*IT_0510;
    const ccomplex_t IT_0512 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0513 = IT_0229*IT_0512;
    const ccomplex_t IT_0514 = (-0.5)*IT_0499 + 0.5*IT_0501 + 0.5*IT_0503 + (
      -0.5)*IT_0505 + (-0.5)*IT_0507 + 0.5*IT_0509 + 0.5*IT_0511 + (-0.5)*IT_0513;
    const ccomplex_t IT_0515 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0516 = IT_0220*IT_0515;
    const ccomplex_t IT_0517 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0518 = IT_0223*IT_0517;
    const ccomplex_t IT_0519 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0520 = IT_0226*IT_0519;
    const ccomplex_t IT_0521 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0522 = IT_0229*IT_0521;
    const ccomplex_t IT_0523 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0524 = IT_0220*IT_0523;
    const ccomplex_t IT_0525 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0526 = IT_0223*IT_0525;
    const ccomplex_t IT_0527 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0528 = IT_0226*IT_0527;
    const ccomplex_t IT_0529 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0530 = IT_0229*IT_0529;
    const ccomplex_t IT_0531 = 0.5*IT_0516 + (-0.5)*IT_0518 + (-0.5)*IT_0520 +
       0.5*IT_0522 + 0.5*IT_0524 + (-0.5)*IT_0526 + (-0.5)*IT_0528 + 0.5*IT_0530;
    const ccomplex_t IT_0532 = m_N_2*m_N_4;
    const ccomplex_t IT_0533 = cpow(IT_0012 + IT_0275 + -IT_0452, 2);
    const ccomplex_t IT_0534 = -IT_0275 + 0.25*IT_0001*(IT_0012 + -IT_0275 +
       IT_0452)*(IT_0012 + IT_0275 + -IT_0452) + 0.25*IT_0001*IT_0533;
    const ccomplex_t IT_0535 = 2*IT_0534;
    const ccomplex_t IT_0536 = cabs(-IT_0275 + 0.25*IT_0001*IT_0533);
    const ccomplex_t IT_0537 = cpow(IT_0536, 0.5);
    const ccomplex_t IT_0538 = ( ((m_A0
) >= (m_N_2 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0514)*(IT_0531*IT_0532 + (-0.5)*IT_0514*IT_0535)
       + conj(IT_0531)*(IT_0514*IT_0532 + (-0.5)*IT_0531*IT_0535))*IT_0537
) : 
      (0
));
    const ccomplex_t IT_0539 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0540 = IT_0220*IT_0539;
    const ccomplex_t IT_0541 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0542 = IT_0223*IT_0541;
    const ccomplex_t IT_0543 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0544 = IT_0226*IT_0543;
    const ccomplex_t IT_0545 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0546 = IT_0229*IT_0545;
    const ccomplex_t IT_0547 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0548 = IT_0220*IT_0547;
    const ccomplex_t IT_0549 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0550 = IT_0223*IT_0549;
    const ccomplex_t IT_0551 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0552 = IT_0226*IT_0551;
    const ccomplex_t IT_0553 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0554 = IT_0229*IT_0553;
    const ccomplex_t IT_0555 = (-0.5)*IT_0540 + 0.5*IT_0542 + 0.5*IT_0544 + (
      -0.5)*IT_0546 + (-0.5)*IT_0548 + 0.5*IT_0550 + 0.5*IT_0552 + (-0.5)*IT_0554;
    const ccomplex_t IT_0556 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0557 = IT_0220*IT_0556;
    const ccomplex_t IT_0558 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0559 = IT_0223*IT_0558;
    const ccomplex_t IT_0560 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0561 = IT_0226*IT_0560;
    const ccomplex_t IT_0562 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0563 = IT_0229*IT_0562;
    const ccomplex_t IT_0564 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0565 = IT_0220*IT_0564;
    const ccomplex_t IT_0566 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0567 = IT_0223*IT_0566;
    const ccomplex_t IT_0568 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0569 = IT_0226*IT_0568;
    const ccomplex_t IT_0570 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0571 = IT_0229*IT_0570;
    const ccomplex_t IT_0572 = 0.5*IT_0557 + (-0.5)*IT_0559 + (-0.5)*IT_0561 +
       0.5*IT_0563 + 0.5*IT_0565 + (-0.5)*IT_0567 + (-0.5)*IT_0569 + 0.5*IT_0571;
    const ccomplex_t IT_0573 = m_N_3*m_N_4;
    const ccomplex_t IT_0574 = cpow(IT_0012 + IT_0345 + -IT_0452, 2);
    const ccomplex_t IT_0575 = -IT_0345 + 0.25*IT_0001*(IT_0012 + -IT_0345 +
       IT_0452)*(IT_0012 + IT_0345 + -IT_0452) + 0.25*IT_0001*IT_0574;
    const ccomplex_t IT_0576 = 2*IT_0575;
    const ccomplex_t IT_0577 = cabs(-IT_0345 + 0.25*IT_0001*IT_0574);
    const ccomplex_t IT_0578 = cpow(IT_0577, 0.5);
    const ccomplex_t IT_0579 = ( ((m_A0
) >= (m_N_3 + m_N_4
)) ? ((-0.25)
      *IT_0000*IT_0001*(conj(IT_0555)*(IT_0572*IT_0573 + (-0.5)*IT_0555*IT_0576)
       + conj(IT_0572)*(IT_0555*IT_0573 + (-0.5)*IT_0572*IT_0576))*IT_0578
) : 
      (0
));
    const ccomplex_t IT_0580 = e_em*m_tau*IT_0002*IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0581 = (-0.5)*IT_0580;
    const ccomplex_t IT_0582 = cpow(IT_0581, 2);
    const ccomplex_t IT_0583 = pow(m_tau, 2);
    const ccomplex_t IT_0584 = 0.5*IT_0012 + -IT_0583;
    const ccomplex_t IT_0585 = cabs(0.25*IT_0012 + -IT_0583);
    const ccomplex_t IT_0586 = cpow(IT_0585, 0.5);
    const ccomplex_t IT_0587 = ( ((m_A0
) >= (2*m_tau
)) ? (0.5*IT_0000
      *IT_0001*IT_0582*(IT_0583 + IT_0584)*IT_0586
) : (0
));
    const ccomplex_t IT_0588 = m_b*e_em*mu_h*IT_0002*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0589 = IT_0007*IT_0588;
    const ccomplex_t IT_0590 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0591 = IT_0007*IT_0590;
    const ccomplex_t IT_0592 = cpow(V_tb, 2);
    const ccomplex_t IT_0593 = A_b*IT_0592*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0594 = IT_0005*IT_0593;
    const ccomplex_t IT_0595 = 1.4142135623731*IT_0594;
    const ccomplex_t IT_0596 = A_b*IT_0592*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0597 = IT_0005*IT_0596;
    const ccomplex_t IT_0598 = 1.4142135623731*IT_0597;
    const ccomplex_t IT_0599 = (-0.5)*IT_0589 + 0.5*IT_0591 + (-0.5)*IT_0595 +
       0.5*IT_0598;
    const ccomplex_t IT_0600 = pow(m_sb_1, 2);
    const ccomplex_t IT_0601 = cabs(0.25*IT_0012 + -IT_0600);
    const ccomplex_t IT_0602 = cpow(IT_0601, 0.5);
    const ccomplex_t IT_0603 = ( ((m_A0
) >= (2*m_sb_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0599*conj(IT_0599)*IT_0602
) : (0
));
    const ccomplex_t IT_0604 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0605 = IT_0007*IT_0604;
    const ccomplex_t IT_0606 = m_b*e_em*mu_h*IT_0002*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0607 = IT_0007*IT_0606;
    const ccomplex_t IT_0608 = A_b*IT_0592*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0609 = IT_0005*IT_0608;
    const ccomplex_t IT_0610 = 1.4142135623731*IT_0609;
    const ccomplex_t IT_0611 = A_b*IT_0592*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0612 = IT_0005*IT_0611;
    const ccomplex_t IT_0613 = 1.4142135623731*IT_0612;
    const ccomplex_t IT_0614 = (-0.5)*IT_0605 + 0.5*IT_0607 + (-0.5)*IT_0610 +
       0.5*IT_0613;
    const ccomplex_t IT_0615 = pow(m_sb_2, 2);
    const ccomplex_t IT_0616 = cabs(0.25*IT_0012 + -IT_0615);
    const ccomplex_t IT_0617 = cpow(IT_0616, 0.5);
    const ccomplex_t IT_0618 = ( ((m_A0
) >= (2*m_sb_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0614*conj(IT_0614)*IT_0617
) : (0
));
    const ccomplex_t IT_0619 = m_b*e_em*mu_h*IT_0002*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0620 = IT_0007*IT_0619;
    const ccomplex_t IT_0621 = m_b*e_em*mu_h*IT_0002*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0622 = IT_0007*IT_0621;
    const ccomplex_t IT_0623 = A_b*IT_0592*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0624 = IT_0005*IT_0623;
    const ccomplex_t IT_0625 = 1.4142135623731*IT_0624;
    const ccomplex_t IT_0626 = A_b*IT_0592*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0627 = IT_0005*IT_0626;
    const ccomplex_t IT_0628 = 1.4142135623731*IT_0627;
    const ccomplex_t IT_0629 = (-0.5)*IT_0620 + 0.5*IT_0622 + (-0.5)*IT_0625 +
       0.5*IT_0628;
    const ccomplex_t IT_0630 = cpow(IT_0012 + IT_0600 + -IT_0615, 2);
    const ccomplex_t IT_0631 = cabs(-IT_0600 + 0.25*IT_0001*IT_0630);
    const ccomplex_t IT_0632 = cpow(IT_0631, 0.5);
    const ccomplex_t IT_0633 = ( ((m_A0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0629*conj(IT_0629)*IT_0632
) : (0
));
    const ccomplex_t IT_0634 = m_b*e_em*mu_h*IT_0002*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0635 = IT_0007*IT_0634;
    const ccomplex_t IT_0636 = m_b*e_em*mu_h*IT_0002*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0637 = IT_0007*IT_0636;
    const ccomplex_t IT_0638 = A_b*IT_0592*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0639 = IT_0005*IT_0638;
    const ccomplex_t IT_0640 = 1.4142135623731*IT_0639;
    const ccomplex_t IT_0641 = A_b*IT_0592*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0642 = IT_0005*IT_0641;
    const ccomplex_t IT_0643 = 1.4142135623731*IT_0642;
    const ccomplex_t IT_0644 = 0.5*IT_0635 + (-0.5)*IT_0637 + 0.5*IT_0640 + (
      -0.5)*IT_0643;
    const ccomplex_t IT_0645 = ( ((m_A0
) >= (m_sb_1 + m_sb_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0632*IT_0644*conj(IT_0644)
) : (0
));
    const ccomplex_t IT_0646 = m_c*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0647 = 0.5*IT_0646;
    const ccomplex_t IT_0648 = cpow(IT_0647, 2);
    const ccomplex_t IT_0649 = pow(m_sc_L, 2);
    const ccomplex_t IT_0650 = pow(m_sc_R, 2);
    const ccomplex_t IT_0651 = cpow(IT_0012 + IT_0649 + -IT_0650, 2);
    const ccomplex_t IT_0652 = cabs(-IT_0649 + 0.25*IT_0001*IT_0651);
    const ccomplex_t IT_0653 = cpow(IT_0652, 0.5);
    const ccomplex_t IT_0654 = ( ((m_A0
) >= (m_sc_L + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0648*IT_0653
) : (0
));
    const ccomplex_t IT_0655 = (-0.5)*IT_0646;
    const ccomplex_t IT_0656 = cpow(IT_0655, 2);
    const ccomplex_t IT_0657 = ( ((m_A0
) >= (m_sc_L + m_sc_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0653*IT_0656
) : (0
));
    const ccomplex_t IT_0658 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_10;
    const ccomplex_t IT_0659 = (-0.5)*IT_0658;
    const ccomplex_t IT_0660 = pow(m_sd_L, 2);
    const ccomplex_t IT_0661 = cpow(IT_0012 + -IT_0600 + IT_0660, 2);
    const ccomplex_t IT_0662 = cabs(-IT_0660 + 0.25*IT_0001*IT_0661);
    const ccomplex_t IT_0663 = cpow(IT_0662, 0.5);
    const ccomplex_t IT_0664 = ( ((m_A0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0659*conj(IT_0659)*IT_0663
) : (0
));
    const ccomplex_t IT_0665 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_10);
    const ccomplex_t IT_0666 = 0.5*IT_0665;
    const ccomplex_t IT_0667 = ( ((m_A0
) >= (m_sb_1 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0663*IT_0666*conj(IT_0666)
) : (0
));
    const ccomplex_t IT_0668 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_11);
    const ccomplex_t IT_0669 = 0.5*IT_0668;
    const ccomplex_t IT_0670 = cpow(IT_0012 + -IT_0615 + IT_0660, 2);
    const ccomplex_t IT_0671 = cabs(-IT_0660 + 0.25*IT_0001*IT_0670);
    const ccomplex_t IT_0672 = cpow(IT_0671, 0.5);
    const ccomplex_t IT_0673 = ( ((m_A0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0669*conj(IT_0669)*IT_0672
) : (0
));
    const ccomplex_t IT_0674 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_11;
    const ccomplex_t IT_0675 = (-0.5)*IT_0674;
    const ccomplex_t IT_0676 = ( ((m_A0
) >= (m_sb_2 + m_sd_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0672*IT_0675*conj(IT_0675)
) : (0
));
    const ccomplex_t IT_0677 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_00;
    const ccomplex_t IT_0678 = 0.5*IT_0677;
    const ccomplex_t IT_0679 = pow(m_sd_R, 2);
    const ccomplex_t IT_0680 = cpow(IT_0012 + -IT_0600 + IT_0679, 2);
    const ccomplex_t IT_0681 = cabs(-IT_0679 + 0.25*IT_0001*IT_0680);
    const ccomplex_t IT_0682 = cpow(IT_0681, 0.5);
    const ccomplex_t IT_0683 = ( ((m_A0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0678*conj(IT_0678)*IT_0682
) : (0
));
    const ccomplex_t IT_0684 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_00);
    const ccomplex_t IT_0685 = (-0.5)*IT_0684;
    const ccomplex_t IT_0686 = ( ((m_A0
) >= (m_sb_1 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0682*IT_0685*conj(IT_0685)
) : (0
));
    const ccomplex_t IT_0687 = 1.4142135623731*A_b*V_tb*V_td*IT_0005*conj
      (U_sb_01);
    const ccomplex_t IT_0688 = (-0.5)*IT_0687;
    const ccomplex_t IT_0689 = cpow(IT_0012 + -IT_0615 + IT_0679, 2);
    const ccomplex_t IT_0690 = cabs(-IT_0679 + 0.25*IT_0001*IT_0689);
    const ccomplex_t IT_0691 = cpow(IT_0690, 0.5);
    const ccomplex_t IT_0692 = ( ((m_A0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0688*conj(IT_0688)*IT_0691
) : (0
));
    const ccomplex_t IT_0693 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_01;
    const ccomplex_t IT_0694 = 0.5*IT_0693;
    const ccomplex_t IT_0695 = ( ((m_A0
) >= (m_sb_2 + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0691*IT_0694*conj(IT_0694)
) : (0
));
    const ccomplex_t IT_0696 = m_d*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0697 = IT_0007*IT_0696;
    const ccomplex_t IT_0698 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0699 = IT_0005*IT_0698;
    const ccomplex_t IT_0700 = 1.4142135623731*IT_0699;
    const ccomplex_t IT_0701 = 0.5*IT_0697 + 0.5*IT_0700;
    const ccomplex_t IT_0702 = cpow(IT_0701, 2);
    const ccomplex_t IT_0703 = cpow(IT_0012 + -IT_0660 + IT_0679, 2);
    const ccomplex_t IT_0704 = cabs(-IT_0679 + 0.25*IT_0001*IT_0703);
    const ccomplex_t IT_0705 = cpow(IT_0704, 0.5);
    const ccomplex_t IT_0706 = ( ((m_A0
) >= (m_sd_L + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0702*IT_0705
) : (0
));
    const ccomplex_t IT_0707 = (-0.5)*IT_0697 + (-0.5)*IT_0700;
    const ccomplex_t IT_0708 = cpow(IT_0707, 2);
    const ccomplex_t IT_0709 = ( ((m_A0
) >= (m_sd_L + m_sd_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0705*IT_0708
) : (0
));
    const ccomplex_t IT_0710 = m_e*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0711 = 0.5*IT_0710;
    const ccomplex_t IT_0712 = cpow(IT_0711, 2);
    const ccomplex_t IT_0713 = pow(m_se_R, 2);
    const ccomplex_t IT_0714 = pow(m_se_L, 2);
    const ccomplex_t IT_0715 = cpow(IT_0012 + IT_0713 + -IT_0714, 2);
    const ccomplex_t IT_0716 = cabs(-IT_0713 + 0.25*IT_0001*IT_0715);
    const ccomplex_t IT_0717 = cpow(IT_0716, 0.5);
    const ccomplex_t IT_0718 = ( ((m_A0
) >= (m_se_L + m_se_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0712*IT_0717
) : (0
));
    const ccomplex_t IT_0719 = (-0.5)*IT_0710;
    const ccomplex_t IT_0720 = cpow(IT_0719, 2);
    const ccomplex_t IT_0721 = ( ((m_A0
) >= (m_se_L + m_se_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0717*IT_0720
) : (0
));
    const ccomplex_t IT_0722 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_10;
    const ccomplex_t IT_0723 = (-0.5)*IT_0722;
    const ccomplex_t IT_0724 = pow(m_ss_L, 2);
    const ccomplex_t IT_0725 = cpow(IT_0012 + -IT_0600 + IT_0724, 2);
    const ccomplex_t IT_0726 = cabs(-IT_0724 + 0.25*IT_0001*IT_0725);
    const ccomplex_t IT_0727 = cpow(IT_0726, 0.5);
    const ccomplex_t IT_0728 = ( ((m_A0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0723*conj(IT_0723)*IT_0727
) : (0
));
    const ccomplex_t IT_0729 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_10);
    const ccomplex_t IT_0730 = 0.5*IT_0729;
    const ccomplex_t IT_0731 = ( ((m_A0
) >= (m_sb_1 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0727*IT_0730*conj(IT_0730)
) : (0
));
    const ccomplex_t IT_0732 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_11);
    const ccomplex_t IT_0733 = 0.5*IT_0732;
    const ccomplex_t IT_0734 = cpow(IT_0012 + -IT_0615 + IT_0724, 2);
    const ccomplex_t IT_0735 = cabs(-IT_0724 + 0.25*IT_0001*IT_0734);
    const ccomplex_t IT_0736 = cpow(IT_0735, 0.5);
    const ccomplex_t IT_0737 = ( ((m_A0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0733*conj(IT_0733)*IT_0736
) : (0
));
    const ccomplex_t IT_0738 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_11;
    const ccomplex_t IT_0739 = (-0.5)*IT_0738;
    const ccomplex_t IT_0740 = ( ((m_A0
) >= (m_sb_2 + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0736*IT_0739*conj(IT_0739)
) : (0
));
    const ccomplex_t IT_0741 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0005;
    const ccomplex_t IT_0742 = 0.5*IT_0741;
    const ccomplex_t IT_0743 = cpow(IT_0012 + IT_0679 + -IT_0724, 2);
    const ccomplex_t IT_0744 = cabs(-IT_0679 + 0.25*IT_0001*IT_0743);
    const ccomplex_t IT_0745 = cpow(IT_0744, 0.5);
    const ccomplex_t IT_0746 = ( ((m_A0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0742*conj(IT_0742)*IT_0745
) : (0
));
    const ccomplex_t IT_0747 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0005;
    const ccomplex_t IT_0748 = (-0.5)*IT_0747;
    const ccomplex_t IT_0749 = ( ((m_A0
) >= (m_sd_R + m_ss_L
)) ? (0.375
      *IT_0000*IT_0001*IT_0745*IT_0748*conj(IT_0748)
) : (0
));
    const ccomplex_t IT_0750 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_00;
    const ccomplex_t IT_0751 = 0.5*IT_0750;
    const ccomplex_t IT_0752 = pow(m_ss_R, 2);
    const ccomplex_t IT_0753 = cpow(IT_0012 + -IT_0600 + IT_0752, 2);
    const ccomplex_t IT_0754 = cabs(-IT_0752 + 0.25*IT_0001*IT_0753);
    const ccomplex_t IT_0755 = cpow(IT_0754, 0.5);
    const ccomplex_t IT_0756 = ( ((m_A0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0751*conj(IT_0751)*IT_0755
) : (0
));
    const ccomplex_t IT_0757 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_00);
    const ccomplex_t IT_0758 = (-0.5)*IT_0757;
    const ccomplex_t IT_0759 = ( ((m_A0
) >= (m_sb_1 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0755*IT_0758*conj(IT_0758)
) : (0
));
    const ccomplex_t IT_0760 = 1.4142135623731*A_b*V_tb*V_ts*IT_0005*conj
      (U_sb_01);
    const ccomplex_t IT_0761 = (-0.5)*IT_0760;
    const ccomplex_t IT_0762 = cpow(IT_0012 + -IT_0615 + IT_0752, 2);
    const ccomplex_t IT_0763 = cabs(-IT_0752 + 0.25*IT_0001*IT_0762);
    const ccomplex_t IT_0764 = cpow(IT_0763, 0.5);
    const ccomplex_t IT_0765 = ( ((m_A0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0761*conj(IT_0761)*IT_0764
) : (0
));
    const ccomplex_t IT_0766 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0005
      *U_sb_01;
    const ccomplex_t IT_0767 = 0.5*IT_0766;
    const ccomplex_t IT_0768 = ( ((m_A0
) >= (m_sb_2 + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0764*IT_0767*conj(IT_0767)
) : (0
));
    const ccomplex_t IT_0769 = 0.5*IT_0747;
    const ccomplex_t IT_0770 = cpow(IT_0012 + -IT_0660 + IT_0752, 2);
    const ccomplex_t IT_0771 = cabs(-IT_0752 + 0.25*IT_0001*IT_0770);
    const ccomplex_t IT_0772 = cpow(IT_0771, 0.5);
    const ccomplex_t IT_0773 = ( ((m_A0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0769*conj(IT_0769)*IT_0772
) : (0
));
    const ccomplex_t IT_0774 = (-0.5)*IT_0741;
    const ccomplex_t IT_0775 = ( ((m_A0
) >= (m_sd_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0772*IT_0774*conj(IT_0774)
) : (0
));
    const ccomplex_t IT_0776 = m_s*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0777 = IT_0007*IT_0776;
    const ccomplex_t IT_0778 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0779 = IT_0005*IT_0778;
    const ccomplex_t IT_0780 = 1.4142135623731*IT_0779;
    const ccomplex_t IT_0781 = 0.5*IT_0777 + 0.5*IT_0780;
    const ccomplex_t IT_0782 = cpow(IT_0781, 2);
    const ccomplex_t IT_0783 = cpow(IT_0012 + -IT_0724 + IT_0752, 2);
    const ccomplex_t IT_0784 = cabs(-IT_0752 + 0.25*IT_0001*IT_0783);
    const ccomplex_t IT_0785 = cpow(IT_0784, 0.5);
    const ccomplex_t IT_0786 = ( ((m_A0
) >= (m_ss_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0782*IT_0785
) : (0
));
    const ccomplex_t IT_0787 = (-0.5)*IT_0777 + (-0.5)*IT_0780;
    const ccomplex_t IT_0788 = cpow(IT_0787, 2);
    const ccomplex_t IT_0789 = ( ((m_A0
) >= (m_ss_L + m_ss_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0785*IT_0788
) : (0
));
    const ccomplex_t IT_0790 = m_t*e_em*mu_h*IT_0002*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0791 = IT_0007*IT_0790;
    const ccomplex_t IT_0792 = m_t*e_em*mu_h*IT_0002*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0793 = IT_0007*IT_0792;
    const ccomplex_t IT_0794 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0795 = IT_0003*IT_0794;
    const ccomplex_t IT_0796 = 1.4142135623731*IT_0795;
    const ccomplex_t IT_0797 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0798 = IT_0003*IT_0797;
    const ccomplex_t IT_0799 = 1.4142135623731*IT_0798;
    const ccomplex_t IT_0800 = (-0.5)*IT_0791 + 0.5*IT_0793 + (-0.5)*IT_0796 +
       0.5*IT_0799;
    const ccomplex_t IT_0801 = pow(m_st_1, 2);
    const ccomplex_t IT_0802 = cabs(0.25*IT_0012 + -IT_0801);
    const ccomplex_t IT_0803 = cpow(IT_0802, 0.5);
    const ccomplex_t IT_0804 = ( ((m_A0
) >= (2*m_st_1
)) ? (0.375*IT_0000
      *IT_0001*IT_0800*conj(IT_0800)*IT_0803
) : (0
));
    const ccomplex_t IT_0805 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0806 = IT_0007*IT_0805;
    const ccomplex_t IT_0807 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0808 = IT_0007*IT_0807;
    const ccomplex_t IT_0809 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0810 = IT_0003*IT_0809;
    const ccomplex_t IT_0811 = 1.4142135623731*IT_0810;
    const ccomplex_t IT_0812 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0813 = IT_0003*IT_0812;
    const ccomplex_t IT_0814 = 1.4142135623731*IT_0813;
    const ccomplex_t IT_0815 = (-0.5)*IT_0806 + 0.5*IT_0808 + (-0.5)*IT_0811 +
       0.5*IT_0814;
    const ccomplex_t IT_0816 = pow(m_st_2, 2);
    const ccomplex_t IT_0817 = cabs(0.25*IT_0012 + -IT_0816);
    const ccomplex_t IT_0818 = cpow(IT_0817, 0.5);
    const ccomplex_t IT_0819 = ( ((m_A0
) >= (2*m_st_2
)) ? (0.375*IT_0000
      *IT_0001*IT_0815*conj(IT_0815)*IT_0818
) : (0
));
    const ccomplex_t IT_0820 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0821 = IT_0007*IT_0820;
    const ccomplex_t IT_0822 = m_t*e_em*mu_h*IT_0002*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0823 = IT_0007*IT_0822;
    const ccomplex_t IT_0824 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0825 = IT_0003*IT_0824;
    const ccomplex_t IT_0826 = 1.4142135623731*IT_0825;
    const ccomplex_t IT_0827 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0828 = IT_0003*IT_0827;
    const ccomplex_t IT_0829 = 1.4142135623731*IT_0828;
    const ccomplex_t IT_0830 = 0.5*IT_0821 + (-0.5)*IT_0823 + 0.5*IT_0826 + (
      -0.5)*IT_0829;
    const ccomplex_t IT_0831 = cpow(IT_0012 + IT_0801 + -IT_0816, 2);
    const ccomplex_t IT_0832 = cabs(-IT_0801 + 0.25*IT_0001*IT_0831);
    const ccomplex_t IT_0833 = cpow(IT_0832, 0.5);
    const ccomplex_t IT_0834 = ( ((m_A0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0830*conj(IT_0830)*IT_0833
) : (0
));
    const ccomplex_t IT_0835 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0836 = IT_0007*IT_0835;
    const ccomplex_t IT_0837 = m_t*e_em*mu_h*IT_0002*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0838 = IT_0007*IT_0837;
    const ccomplex_t IT_0839 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0840 = IT_0003*IT_0839;
    const ccomplex_t IT_0841 = 1.4142135623731*IT_0840;
    const ccomplex_t IT_0842 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0843 = IT_0003*IT_0842;
    const ccomplex_t IT_0844 = 1.4142135623731*IT_0843;
    const ccomplex_t IT_0845 = (-0.5)*IT_0836 + 0.5*IT_0838 + (-0.5)*IT_0841 +
       0.5*IT_0844;
    const ccomplex_t IT_0846 = ( ((m_A0
) >= (m_st_1 + m_st_2
)) ? (0.375
      *IT_0000*IT_0001*IT_0833*IT_0845*conj(IT_0845)
) : (0
));
    const ccomplex_t IT_0847 = m_u*e_em*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0848 = 0.5*IT_0847;
    const ccomplex_t IT_0849 = cpow(IT_0848, 2);
    const ccomplex_t IT_0850 = pow(m_su_L, 2);
    const ccomplex_t IT_0851 = pow(m_su_R, 2);
    const ccomplex_t IT_0852 = cpow(IT_0012 + IT_0850 + -IT_0851, 2);
    const ccomplex_t IT_0853 = cabs(-IT_0850 + 0.25*IT_0001*IT_0852);
    const ccomplex_t IT_0854 = cpow(IT_0853, 0.5);
    const ccomplex_t IT_0855 = ( ((m_A0
) >= (m_su_L + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0849*IT_0854
) : (0
));
    const ccomplex_t IT_0856 = (-0.5)*IT_0847;
    const ccomplex_t IT_0857 = cpow(IT_0856, 2);
    const ccomplex_t IT_0858 = ( ((m_A0
) >= (m_su_L + m_su_R
)) ? (0.375
      *IT_0000*IT_0001*IT_0854*IT_0857
) : (0
));
    const ccomplex_t IT_0859 = e_em*m_mu*mu_h*IT_0002*IT_0007;
    const ccomplex_t IT_0860 = 0.5*IT_0859;
    const ccomplex_t IT_0861 = cpow(IT_0860, 2);
    const ccomplex_t IT_0862 = pow(m_smu_R, 2);
    const ccomplex_t IT_0863 = pow(m_smu_L, 2);
    const ccomplex_t IT_0864 = cpow(IT_0012 + IT_0862 + -IT_0863, 2);
    const ccomplex_t IT_0865 = cabs(-IT_0862 + 0.25*IT_0001*IT_0864);
    const ccomplex_t IT_0866 = cpow(IT_0865, 0.5);
    const ccomplex_t IT_0867 = ( ((m_A0
) >= (m_smu_L + m_smu_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0861*IT_0866
) : (0
));
    const ccomplex_t IT_0868 = (-0.5)*IT_0859;
    const ccomplex_t IT_0869 = cpow(IT_0868, 2);
    const ccomplex_t IT_0870 = ( ((m_A0
) >= (m_smu_L + m_smu_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0866*IT_0869
) : (0
));
    const ccomplex_t IT_0871 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0872 = IT_0007*IT_0871;
    const ccomplex_t IT_0873 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0874 = IT_0007*IT_0873;
    const ccomplex_t IT_0875 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0876 = IT_0005*IT_0875;
    const ccomplex_t IT_0877 = 1.4142135623731*IT_0876;
    const ccomplex_t IT_0878 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0879 = IT_0005*IT_0878;
    const ccomplex_t IT_0880 = 1.4142135623731*IT_0879;
    const ccomplex_t IT_0881 = (-0.5)*IT_0872 + 0.5*IT_0874 + (-0.5)*IT_0877 +
       0.5*IT_0880;
    const ccomplex_t IT_0882 = pow(m_stau_1, 2);
    const ccomplex_t IT_0883 = cabs(0.25*IT_0012 + -IT_0882);
    const ccomplex_t IT_0884 = cpow(IT_0883, 0.5);
    const ccomplex_t IT_0885 = ( ((m_A0
) >= (2*m_stau_1
)) ? (0.125*IT_0000
      *IT_0001*IT_0881*conj(IT_0881)*IT_0884
) : (0
));
    const ccomplex_t IT_0886 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0887 = IT_0007*IT_0886;
    const ccomplex_t IT_0888 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0889 = IT_0007*IT_0888;
    const ccomplex_t IT_0890 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0891 = IT_0005*IT_0890;
    const ccomplex_t IT_0892 = 1.4142135623731*IT_0891;
    const ccomplex_t IT_0893 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0894 = IT_0005*IT_0893;
    const ccomplex_t IT_0895 = 1.4142135623731*IT_0894;
    const ccomplex_t IT_0896 = (-0.5)*IT_0887 + 0.5*IT_0889 + (-0.5)*IT_0892 +
       0.5*IT_0895;
    const ccomplex_t IT_0897 = pow(m_stau_2, 2);
    const ccomplex_t IT_0898 = cabs(0.25*IT_0012 + -IT_0897);
    const ccomplex_t IT_0899 = cpow(IT_0898, 0.5);
    const ccomplex_t IT_0900 = ( ((m_A0
) >= (2*m_stau_2
)) ? (0.125*IT_0000
      *IT_0001*IT_0896*conj(IT_0896)*IT_0899
) : (0
));
    const ccomplex_t IT_0901 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0902 = IT_0007*IT_0901;
    const ccomplex_t IT_0903 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0904 = IT_0007*IT_0903;
    const ccomplex_t IT_0905 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0906 = IT_0005*IT_0905;
    const ccomplex_t IT_0907 = 1.4142135623731*IT_0906;
    const ccomplex_t IT_0908 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0909 = IT_0005*IT_0908;
    const ccomplex_t IT_0910 = 1.4142135623731*IT_0909;
    const ccomplex_t IT_0911 = (-0.5)*IT_0902 + 0.5*IT_0904 + (-0.5)*IT_0907 +
       0.5*IT_0910;
    const ccomplex_t IT_0912 = cpow(IT_0012 + IT_0882 + -IT_0897, 2);
    const ccomplex_t IT_0913 = cabs(-IT_0882 + 0.25*IT_0001*IT_0912);
    const ccomplex_t IT_0914 = cpow(IT_0913, 0.5);
    const ccomplex_t IT_0915 = ( ((m_A0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0911*conj(IT_0911)*IT_0914
) : (0
));
    const ccomplex_t IT_0916 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0917 = IT_0007*IT_0916;
    const ccomplex_t IT_0918 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0919 = IT_0007*IT_0918;
    const ccomplex_t IT_0920 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0921 = IT_0005*IT_0920;
    const ccomplex_t IT_0922 = 1.4142135623731*IT_0921;
    const ccomplex_t IT_0923 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0924 = IT_0005*IT_0923;
    const ccomplex_t IT_0925 = 1.4142135623731*IT_0924;
    const ccomplex_t IT_0926 = 0.5*IT_0917 + (-0.5)*IT_0919 + 0.5*IT_0922 + (
      -0.5)*IT_0925;
    const ccomplex_t IT_0927 = ( ((m_A0
) >= (m_stau_1 + m_stau_2
)) ? (0.125
      *IT_0000*IT_0001*IT_0914*IT_0926*conj(IT_0926)
) : (0
));
    const ccomplex_t IT_0928 = IT_0016 + IT_0025 + IT_0033 + IT_0041 + IT_0049
       + IT_0057 + IT_0065 + IT_0083 + IT_0095 + IT_0106 + IT_0114 + 2*IT_0127 +
       IT_0135 + IT_0157 + IT_0179 + IT_0202 + IT_0219 + IT_0252 + IT_0281 +
       IT_0322 + IT_0351 + IT_0392 + IT_0433 + IT_0456 + IT_0497 + IT_0538 +
       IT_0579 + IT_0587 + IT_0603 + IT_0618 + IT_0633 + IT_0645 + IT_0654 +
       IT_0657 + IT_0664 + IT_0667 + IT_0673 + IT_0676 + IT_0683 + IT_0686 +
       IT_0692 + IT_0695 + IT_0706 + IT_0709 + IT_0718 + IT_0721 + IT_0728 +
       IT_0731 + IT_0737 + IT_0740 + IT_0746 + IT_0749 + IT_0756 + IT_0759 +
       IT_0765 + IT_0768 + IT_0773 + IT_0775 + IT_0786 + IT_0789 + IT_0804 +
       IT_0819 + IT_0834 + IT_0846 + IT_0855 + IT_0858 + IT_0867 + IT_0870 +
       IT_0885 + IT_0900 + IT_0915 + IT_0927;
    return create_ccomplex_return(IT_0928);
}

