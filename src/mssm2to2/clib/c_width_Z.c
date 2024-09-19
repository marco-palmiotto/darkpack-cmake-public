#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
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
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0004, -1);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0011, 2);
    const ccomplex_t IT_0013 = pow(m_Z, 2);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = IT_0013*IT_0014;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0015 + IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = cabs(0.25*IT_0013);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_Z
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0012*(IT_0014 + IT_0019)*IT_0021
) : (0
));
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0007;
    const ccomplex_t IT_0024 = (-2)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = pow(m_W, 2);
    const ccomplex_t IT_0028 = IT_0001*IT_0015;
    const ccomplex_t IT_0029 = pow(m_W, -4);
    const ccomplex_t IT_0030 = 0.5*IT_0013 + -IT_0027;
    const ccomplex_t IT_0031 = cpow(IT_0030, 2);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (IT_0027 + -IT_0028)*(1 + 0.5*IT_0032);
    const ccomplex_t IT_0034 = 2*IT_0023;
    const ccomplex_t IT_0035 = pow(m_W, -2);
    const ccomplex_t IT_0036 = IT_0001*IT_0015*IT_0035;
    const ccomplex_t IT_0037 = IT_0015*IT_0035;
    const ccomplex_t IT_0038 = (1 + 0.5*IT_0036)*(IT_0013 + -IT_0037);
    const ccomplex_t IT_0039 = cpow(IT_0014, 3);
    const ccomplex_t IT_0040 = IT_0001*IT_0035*IT_0039;
    const ccomplex_t IT_0041 = IT_0014*IT_0036;
    const ccomplex_t IT_0042 = IT_0014*IT_0032;
    const ccomplex_t IT_0043 = IT_0030*IT_0035;
    const ccomplex_t IT_0044 = IT_0014*IT_0043;
    const ccomplex_t IT_0045 = IT_0001*IT_0029*IT_0039;
    const ccomplex_t IT_0046 = IT_0030*IT_0045;
    const ccomplex_t IT_0047 = -IT_0023;
    const ccomplex_t IT_0048 = IT_0031*IT_0035;
    const ccomplex_t IT_0049 = (1 + 0.5*IT_0036)*(IT_0027 + -IT_0048);
    const ccomplex_t IT_0050 = (1 + 0.5*IT_0036)*(IT_0014 + -IT_0044);
    const ccomplex_t IT_0051 = IT_0014*IT_0029;
    const ccomplex_t IT_0052 = IT_0014*IT_0051;
    const ccomplex_t IT_0053 = IT_0001*IT_0014;
    const ccomplex_t IT_0054 = IT_0014*IT_0053;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0030 + IT_0055;
    const ccomplex_t IT_0057 = IT_0052*IT_0056;
    const ccomplex_t IT_0058 = -IT_0032 + -IT_0036;
    const ccomplex_t IT_0059 = IT_0001*IT_0029;
    const ccomplex_t IT_0060 = IT_0030*IT_0059;
    const ccomplex_t IT_0061 = IT_0014*IT_0060;
    const ccomplex_t IT_0062 = IT_0014*IT_0061;
    const ccomplex_t IT_0063 = IT_0058 + IT_0062;
    const ccomplex_t IT_0064 = 1 + IT_0063;
    const ccomplex_t IT_0065 = IT_0014*IT_0064;
    const ccomplex_t IT_0066 = (-2)*IT_0036;
    const ccomplex_t IT_0067 = IT_0062 + IT_0066;
    const ccomplex_t IT_0068 = 1 + IT_0067;
    const ccomplex_t IT_0069 = IT_0014*IT_0068;
    const ccomplex_t IT_0070 = IT_0001*IT_0035;
    const ccomplex_t IT_0071 = IT_0030*IT_0070;
    const ccomplex_t IT_0072 = IT_0014*IT_0071;
    const ccomplex_t IT_0073 = -IT_0044;
    const ccomplex_t IT_0074 = IT_0014 + IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0026*IT_0033 +
       0.666666666666667*IT_0034*(IT_0034*IT_0038 + -IT_0025*(IT_0014 + -IT_0040
       + -IT_0041 + -IT_0042 + IT_0044 + IT_0046)) + 0.666666666666667*IT_0047*
      (IT_0047*(IT_0038 + IT_0049 + 2*IT_0050) + IT_0034*(IT_0057 + -IT_0065) + 
      -IT_0025*(IT_0069 + -IT_0075));
    const ccomplex_t IT_0077 = cabs(0.25*IT_0013 + -IT_0027);
    const ccomplex_t IT_0078 = cpow(IT_0077, 0.5);
    const ccomplex_t IT_0079 = ( ((m_Z
) >= (2*m_W
)) ? (0.125*IT_0000*IT_0001
      *IT_0076*IT_0078
) : (0
));
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0081 = 0.333333333333333*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0081, 2);
    const ccomplex_t IT_0083 = pow(m_b, 2);
    const ccomplex_t IT_0084 = 0.5*IT_0013 + -IT_0083;
    const ccomplex_t IT_0085 = IT_0013*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0015 + IT_0086;
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = IT_0084 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0006 + 3*IT_0009);
    const ccomplex_t IT_0091 = (-0.166666666666667)*IT_0090;
    const ccomplex_t IT_0092 = cabs(0.25*IT_0013 + -IT_0083);
    const ccomplex_t IT_0093 = cpow(IT_0092, 0.5);
    const ccomplex_t IT_0094 = ( ((m_Z
) >= (2*m_b
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0089 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0083 + (-4)
      *IT_0089*IT_0091))*IT_0093
) : (0
));
    const ccomplex_t IT_0095 = (-0.666666666666667)*IT_0080;
    const ccomplex_t IT_0096 = cpow(IT_0095, 2);
    const ccomplex_t IT_0097 = pow(m_c, 2);
    const ccomplex_t IT_0098 = 0.5*IT_0013 + -IT_0097;
    const ccomplex_t IT_0099 = IT_0013*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0015 + IT_0100;
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = IT_0098 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0006 + (-3)*IT_0009);
    const ccomplex_t IT_0105 = (-0.166666666666667)*IT_0104;
    const ccomplex_t IT_0106 = cabs(0.25*IT_0013 + -IT_0097);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_Z
) >= (2*m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0096*IT_0103 + (-0.25)*IT_0105*((-12)*IT_0095*IT_0097 + (-4)
      *IT_0103*IT_0105))*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = pow(m_d, 2);
    const ccomplex_t IT_0110 = 0.5*IT_0013 + -IT_0109;
    const ccomplex_t IT_0111 = IT_0013*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0015 + IT_0112;
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = IT_0110 + IT_0114;
    const ccomplex_t IT_0116 = cabs(0.25*IT_0013 + -IT_0109);
    const ccomplex_t IT_0117 = cpow(IT_0116, 0.5);
    const ccomplex_t IT_0118 = ( ((m_Z
) >= (2*m_d
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0115 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0109 + (-4)
      *IT_0091*IT_0115))*IT_0117
) : (0
));
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009);
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = pow(m_e, 2);
    const ccomplex_t IT_0123 = 0.5*IT_0013 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0013*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0015 + IT_0125;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = IT_0123 + IT_0127;
    const ccomplex_t IT_0129 = cabs(0.25*IT_0013 + -IT_0122);
    const ccomplex_t IT_0130 = cpow(IT_0129, 0.5);
    const ccomplex_t IT_0131 = ( ((m_Z
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *((-1.33333333333333)*IT_0121*IT_0128 + (-4)*IT_0080*(IT_0120*IT_0122 +
       0.333333333333333*IT_0080*IT_0128))*IT_0130
) : (0
));
    const ccomplex_t IT_0132 = m_W*e_em;
    const ccomplex_t IT_0133 = cos(alpha);
    const ccomplex_t IT_0134 = cpow(IT_0002, -4);
    const ccomplex_t IT_0135 = sin(beta);
    const ccomplex_t IT_0136 = IT_0004*IT_0133*IT_0134*IT_0135;
    const ccomplex_t IT_0137 = cos(beta);
    const ccomplex_t IT_0138 = sin(alpha);
    const ccomplex_t IT_0139 = IT_0004*IT_0134*IT_0137*IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0004, -3);
    const ccomplex_t IT_0141 = IT_0133*IT_0135*IT_0140;
    const ccomplex_t IT_0142 = IT_0137*IT_0138*IT_0140;
    const ccomplex_t IT_0143 = cpow(IT_0002, -2);
    const ccomplex_t IT_0144 = IT_0007*IT_0133*IT_0135*IT_0143;
    const ccomplex_t IT_0145 = IT_0007*IT_0137*IT_0138*IT_0143;
    const ccomplex_t IT_0146 = cpow(IT_0004, -2);
    const ccomplex_t IT_0147 = cpow(IT_0143 + IT_0146, -1);
    const ccomplex_t IT_0148 = (IT_0136 + -IT_0139 + IT_0141 + -IT_0142 + 2
      *IT_0144 + (-2)*IT_0145)*IT_0147;
    const ccomplex_t IT_0149 = IT_0132*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = cpow(IT_0150, 2);
    const ccomplex_t IT_0152 = pow(m_Z, -4);
    const ccomplex_t IT_0153 = pow(m_h, 2);
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0013 + IT_0154;
    const ccomplex_t IT_0156 = cpow(IT_0155, 2);
    const ccomplex_t IT_0157 = IT_0152*IT_0156;
    const ccomplex_t IT_0158 = cpow(IT_0013 + (-0.5)*IT_0153, 2);
    const ccomplex_t IT_0159 = cabs(-IT_0013 + IT_0001*IT_0158);
    const ccomplex_t IT_0160 = cpow(IT_0159, 0.5);
    const ccomplex_t IT_0161 = ( ((m_Z
) >= (m_Z + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0151*(1 + 0.5*IT_0157)*IT_0160
) :
       (0
));
    const ccomplex_t IT_0162 = pow(m_s, 2);
    const ccomplex_t IT_0163 = 0.5*IT_0013 + -IT_0162;
    const ccomplex_t IT_0164 = IT_0013*IT_0163;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = IT_0015 + IT_0165;
    const ccomplex_t IT_0167 = IT_0001*IT_0166;
    const ccomplex_t IT_0168 = IT_0163 + IT_0167;
    const ccomplex_t IT_0169 = cabs(0.25*IT_0013 + -IT_0162);
    const ccomplex_t IT_0170 = cpow(IT_0169, 0.5);
    const ccomplex_t IT_0171 = ( ((m_Z
) >= (2*m_s
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0168 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0162 + (-4)
      *IT_0091*IT_0168))*IT_0170
) : (0
));
    const ccomplex_t IT_0172 = pow(m_t, 2);
    const ccomplex_t IT_0173 = 0.5*IT_0013 + -IT_0172;
    const ccomplex_t IT_0174 = IT_0013*IT_0173;
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = IT_0015 + IT_0175;
    const ccomplex_t IT_0177 = IT_0001*IT_0176;
    const ccomplex_t IT_0178 = IT_0173 + IT_0177;
    const ccomplex_t IT_0179 = cabs(0.25*IT_0013 + -IT_0172);
    const ccomplex_t IT_0180 = cpow(IT_0179, 0.5);
    const ccomplex_t IT_0181 = ( ((m_Z
) >= (2*m_t
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0096*IT_0178 + (-0.25)*IT_0105*((-12)*IT_0095*IT_0172 + (-4)
      *IT_0105*IT_0178))*IT_0180
) : (0
));
    const ccomplex_t IT_0182 = cpow(IT_0105, 2);
    const ccomplex_t IT_0183 = pow(m_u, 2);
    const ccomplex_t IT_0184 = 0.5*IT_0013 + -IT_0183;
    const ccomplex_t IT_0185 = IT_0013*IT_0184;
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0015 + IT_0186;
    const ccomplex_t IT_0188 = IT_0001*IT_0187;
    const ccomplex_t IT_0189 = IT_0184 + IT_0188;
    const ccomplex_t IT_0190 = cabs(0.25*IT_0013 + -IT_0183);
    const ccomplex_t IT_0191 = cpow(IT_0190, 0.5);
    const ccomplex_t IT_0192 = ( ((m_Z
) >= (2*m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0182*IT_0189 + (-0.25)*IT_0095*((-12)*IT_0105*IT_0183 + (-4)
      *IT_0095*IT_0189))*IT_0191
) : (0
));
    const ccomplex_t IT_0193 = e_em*IT_0003*IT_0007*(IT_0133*IT_0137 + IT_0135
      *IT_0138);
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = cpow(IT_0194, 2);
    const ccomplex_t IT_0196 = pow(m_A0, 2);
    const ccomplex_t IT_0197 = 0.5*IT_0013 + (-0.5)*IT_0153 + 0.5*IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0197, 2);
    const ccomplex_t IT_0199 = IT_0001*IT_0198;
    const ccomplex_t IT_0200 = cpow(IT_0013 + -IT_0153 + IT_0196, 2);
    const ccomplex_t IT_0201 = cabs(-IT_0196 + 0.25*IT_0001*IT_0200);
    const ccomplex_t IT_0202 = cpow(IT_0201, 0.5);
    const ccomplex_t IT_0203 = ( ((m_Z
) >= (m_h + m_A0
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0195*(IT_0196 + -IT_0199)*IT_0202

      ) : (0
));
    const ccomplex_t IT_0204 = e_em*IT_0003*IT_0007*(IT_0133*IT_0135 + 
      -IT_0137*IT_0138);
    const ccomplex_t IT_0205 = cpow(IT_0204, 2);
    const ccomplex_t IT_0206 = pow(m_H0, 2);
    const ccomplex_t IT_0207 = 0.5*IT_0013 + 0.5*IT_0196 + (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = cpow(IT_0207, 2);
    const ccomplex_t IT_0209 = IT_0001*IT_0208;
    const ccomplex_t IT_0210 = cpow(IT_0013 + IT_0196 + -IT_0206, 2);
    const ccomplex_t IT_0211 = cabs(-IT_0196 + 0.25*IT_0001*IT_0210);
    const ccomplex_t IT_0212 = cpow(IT_0211, 0.5);
    const ccomplex_t IT_0213 = ( ((m_Z
) >= (m_A0 + m_H0
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0205*(IT_0196 + -IT_0209)*IT_0212

      ) : (0
));
    const ccomplex_t IT_0214 = IT_0004*IT_0133*IT_0134*IT_0137;
    const ccomplex_t IT_0215 = IT_0004*IT_0134*IT_0135*IT_0138;
    const ccomplex_t IT_0216 = IT_0133*IT_0137*IT_0140;
    const ccomplex_t IT_0217 = IT_0135*IT_0138*IT_0140;
    const ccomplex_t IT_0218 = IT_0007*IT_0133*IT_0137*IT_0143;
    const ccomplex_t IT_0219 = IT_0007*IT_0135*IT_0138*IT_0143;
    const ccomplex_t IT_0220 = IT_0147*(IT_0214 + IT_0215 + IT_0216 + IT_0217 
      + 2*IT_0218 + 2*IT_0219);
    const ccomplex_t IT_0221 = IT_0132*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = cpow(IT_0222, 2);
    const ccomplex_t IT_0224 = (-0.5)*IT_0206;
    const ccomplex_t IT_0225 = IT_0013 + IT_0224;
    const ccomplex_t IT_0226 = cpow(IT_0225, 2);
    const ccomplex_t IT_0227 = IT_0152*IT_0226;
    const ccomplex_t IT_0228 = cpow(IT_0013 + (-0.5)*IT_0206, 2);
    const ccomplex_t IT_0229 = cabs(-IT_0013 + IT_0001*IT_0228);
    const ccomplex_t IT_0230 = cpow(IT_0229, 0.5);
    const ccomplex_t IT_0231 = ( ((m_Z
) >= (m_Z + m_H0
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0223*(1 + 0.5*IT_0227)*IT_0230
) :
       (0
));
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0004 + 
      -IT_0002*IT_0007);
    const ccomplex_t IT_0233 = cpow(IT_0232, 2);
    const ccomplex_t IT_0234 = pow(m_Hp, 2);
    const ccomplex_t IT_0235 = cabs(0.25*IT_0013 + -IT_0234);
    const ccomplex_t IT_0236 = cpow(IT_0235, 0.5);
    const ccomplex_t IT_0237 = ( ((m_Z
) >= (2*m_Hp
)) ? ((-0.0416666666666667
      )*IT_0000*IT_0001*IT_0233*(IT_0028 + -IT_0234)*IT_0236
) : (0
));
    const ccomplex_t IT_0238 = pow(m_mu, 2);
    const ccomplex_t IT_0239 = 0.5*IT_0013 + -IT_0238;
    const ccomplex_t IT_0240 = IT_0013*IT_0239;
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = IT_0015 + IT_0241;
    const ccomplex_t IT_0243 = IT_0001*IT_0242;
    const ccomplex_t IT_0244 = IT_0239 + IT_0243;
    const ccomplex_t IT_0245 = cabs(0.25*IT_0013 + -IT_0238);
    const ccomplex_t IT_0246 = cpow(IT_0245, 0.5);
    const ccomplex_t IT_0247 = ( ((m_Z
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0121*IT_0244 + (-4)*IT_0080*(IT_0120
      *IT_0238 + 0.333333333333333*IT_0080*IT_0244))*IT_0246
) : (0
));
    const ccomplex_t IT_0248 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0249 = IT_0005*IT_0248;
    const ccomplex_t IT_0250 = IT_0008*IT_0248;
    const ccomplex_t IT_0251 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0252 = IT_0008*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0249 + -IT_0250 + (-2)
      *IT_0252);
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = cpow(IT_0255, 2);
    const ccomplex_t IT_0257 = pow(m_C_1, 2);
    const ccomplex_t IT_0258 = 0.5*IT_0013 + -IT_0257;
    const ccomplex_t IT_0259 = IT_0013*IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = IT_0015 + IT_0260;
    const ccomplex_t IT_0262 = IT_0001*IT_0261;
    const ccomplex_t IT_0263 = IT_0258 + IT_0262;
    const ccomplex_t IT_0264 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0265 = IT_0005*IT_0264;
    const ccomplex_t IT_0266 = IT_0008*IT_0264;
    const ccomplex_t IT_0267 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0268 = IT_0008*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*(IT_0265 + -IT_0266 + (-2)
      *IT_0268);
    const ccomplex_t IT_0270 = (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = cabs(0.25*IT_0013 + -IT_0257);
    const ccomplex_t IT_0272 = cpow(IT_0271, 0.5);
    const ccomplex_t IT_0273 = ( ((m_Z
) >= (2*m_C_1
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0256*IT_0263 + (-4)*IT_0270*(IT_0255
      *IT_0257 + 0.333333333333333*IT_0263*IT_0270))*IT_0272
) : (0
));
    const ccomplex_t IT_0274 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0275 = IT_0005*IT_0274;
    const ccomplex_t IT_0276 = IT_0008*IT_0274;
    const ccomplex_t IT_0277 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0278 = IT_0008*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0275 + -IT_0276 + (-2)
      *IT_0278);
    const ccomplex_t IT_0280 = 0.5*IT_0279;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = cpow(IT_0281, 2);
    const ccomplex_t IT_0283 = pow(m_C_2, 2);
    const ccomplex_t IT_0284 = 0.5*IT_0013 + -IT_0283;
    const ccomplex_t IT_0285 = IT_0013*IT_0284;
    const ccomplex_t IT_0286 = (-0.5)*IT_0285;
    const ccomplex_t IT_0287 = IT_0015 + IT_0286;
    const ccomplex_t IT_0288 = IT_0001*IT_0287;
    const ccomplex_t IT_0289 = IT_0284 + IT_0288;
    const ccomplex_t IT_0290 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0291 = IT_0005*IT_0290;
    const ccomplex_t IT_0292 = IT_0008*IT_0290;
    const ccomplex_t IT_0293 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0294 = IT_0008*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*(IT_0291 + -IT_0292 + (-2)
      *IT_0294);
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = cabs(0.25*IT_0013 + -IT_0283);
    const ccomplex_t IT_0298 = cpow(IT_0297, 0.5);
    const ccomplex_t IT_0299 = ( ((m_Z
) >= (2*m_C_2
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0282*IT_0289 + (-4)*IT_0296*(IT_0281
      *IT_0283 + 0.333333333333333*IT_0289*IT_0296))*IT_0298
) : (0
));
    const ccomplex_t IT_0300 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0301 = IT_0005*IT_0300;
    const ccomplex_t IT_0302 = IT_0008*IT_0300;
    const ccomplex_t IT_0303 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0304 = IT_0008*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0301 + -IT_0302 + (-2)
      *IT_0304);
    const ccomplex_t IT_0306 = (-0.5)*IT_0305;
    const ccomplex_t IT_0307 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0308 = IT_0005*IT_0307;
    const ccomplex_t IT_0309 = IT_0008*IT_0307;
    const ccomplex_t IT_0310 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0311 = IT_0008*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*(IT_0308 + -IT_0309 + (-2)
      *IT_0311);
    const ccomplex_t IT_0313 = 0.5*IT_0312;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = m_C_1*m_C_2;
    const ccomplex_t IT_0316 = cpow(IT_0013 + IT_0257 + -IT_0283, 2);
    const ccomplex_t IT_0317 = -IT_0257 + 0.25*IT_0001*(IT_0013 + -IT_0257 +
       IT_0283)*(IT_0013 + IT_0257 + -IT_0283) + 0.25*IT_0001*IT_0316;
    const ccomplex_t IT_0318 = 0.5*IT_0013 + 0.5*IT_0257 + (-0.5)*IT_0283;
    const ccomplex_t IT_0319 = 0.5*IT_0013 + (-0.5)*IT_0257 + 0.5*IT_0283;
    const ccomplex_t IT_0320 = IT_0318*IT_0319;
    const ccomplex_t IT_0321 = IT_0013*IT_0317;
    const ccomplex_t IT_0322 = (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = IT_0001*IT_0323;
    const ccomplex_t IT_0325 = IT_0317 + IT_0324;
    const ccomplex_t IT_0326 = cabs(-IT_0257 + 0.25*IT_0001*IT_0316);
    const ccomplex_t IT_0327 = cpow(IT_0326, 0.5);
    const ccomplex_t IT_0328 = ( ((m_Z
) >= (m_C_1 + m_C_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0306)*(IT_0314*IT_0315 + 0.666666666666667*IT_0306
      *IT_0325) + conj(IT_0314)*(IT_0306*IT_0315 + 0.666666666666667*IT_0314
      *IT_0325))*IT_0327
) : (0
));
    const ccomplex_t IT_0329 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0330 = IT_0005*IT_0329;
    const ccomplex_t IT_0331 = IT_0008*IT_0329;
    const ccomplex_t IT_0332 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0333 = IT_0008*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*(IT_0330 + -IT_0331 + (-2)
      *IT_0333);
    const ccomplex_t IT_0335 = (-0.5)*IT_0334;
    const ccomplex_t IT_0336 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0337 = IT_0005*IT_0336;
    const ccomplex_t IT_0338 = IT_0008*IT_0336;
    const ccomplex_t IT_0339 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0340 = IT_0008*IT_0339;
    const ccomplex_t IT_0341 = (0 + _Complex_I*1)*(IT_0337 + -IT_0338 + (-2)
      *IT_0340);
    const ccomplex_t IT_0342 = 0.5*IT_0341;
    const ccomplex_t IT_0343 = -IT_0342;
    const ccomplex_t IT_0344 = ( ((m_Z
) >= (m_C_1 + m_C_2
)) ? (0.25*IT_0000
      *IT_0001*IT_0327*(conj(IT_0335)*(0.666666666666667*IT_0325*IT_0335 +
       IT_0315*IT_0343) + (IT_0315*IT_0335 + 0.666666666666667*IT_0325*IT_0343)
      *conj(IT_0343))
) : (0
));
    const ccomplex_t IT_0345 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0346 = IT_0005*IT_0345;
    const ccomplex_t IT_0347 = IT_0008*IT_0345;
    const ccomplex_t IT_0348 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0349 = IT_0005*IT_0348;
    const ccomplex_t IT_0350 = IT_0008*IT_0348;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*(IT_0346 + IT_0347 + 
      -IT_0349 + -IT_0350);
    const ccomplex_t IT_0352 = 0.5*IT_0351;
    const ccomplex_t IT_0353 = cpow(IT_0352, 2);
    const ccomplex_t IT_0354 = pow(m_N_1, 2);
    const ccomplex_t IT_0355 = 0.5*IT_0013 + -IT_0354;
    const ccomplex_t IT_0356 = IT_0013*IT_0355;
    const ccomplex_t IT_0357 = (-0.5)*IT_0356;
    const ccomplex_t IT_0358 = IT_0015 + IT_0357;
    const ccomplex_t IT_0359 = IT_0001*IT_0358;
    const ccomplex_t IT_0360 = IT_0355 + IT_0359;
    const ccomplex_t IT_0361 = (-0.5)*IT_0351;
    const ccomplex_t IT_0362 = cabs(0.25*IT_0013 + -IT_0354);
    const ccomplex_t IT_0363 = cpow(IT_0362, 0.5);
    const ccomplex_t IT_0364 = ( ((m_Z
) >= (2*m_N_1
)) ? (0.125*IT_0000
      *IT_0001*((-0.666666666666667)*IT_0353*IT_0360 + (-2)*IT_0361*(IT_0352
      *IT_0354 + 0.333333333333333*IT_0360*IT_0361))*IT_0363
) : (0
));
    const ccomplex_t IT_0365 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0366 = IT_0005*IT_0365;
    const ccomplex_t IT_0367 = IT_0008*IT_0365;
    const ccomplex_t IT_0368 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0369 = IT_0005*IT_0368;
    const ccomplex_t IT_0370 = IT_0008*IT_0368;
    const ccomplex_t IT_0371 = (0 + _Complex_I*1)*(IT_0366 + IT_0367 + 
      -IT_0369 + -IT_0370);
    const ccomplex_t IT_0372 = 0.5*IT_0371;
    const ccomplex_t IT_0373 = cpow(IT_0372, 2);
    const ccomplex_t IT_0374 = pow(m_N_2, 2);
    const ccomplex_t IT_0375 = 0.5*IT_0013 + -IT_0374;
    const ccomplex_t IT_0376 = IT_0013*IT_0375;
    const ccomplex_t IT_0377 = (-0.5)*IT_0376;
    const ccomplex_t IT_0378 = IT_0015 + IT_0377;
    const ccomplex_t IT_0379 = IT_0001*IT_0378;
    const ccomplex_t IT_0380 = IT_0375 + IT_0379;
    const ccomplex_t IT_0381 = (-0.5)*IT_0371;
    const ccomplex_t IT_0382 = cabs(0.25*IT_0013 + -IT_0374);
    const ccomplex_t IT_0383 = cpow(IT_0382, 0.5);
    const ccomplex_t IT_0384 = ( ((m_Z
) >= (2*m_N_2
)) ? (0.125*IT_0000
      *IT_0001*((-0.666666666666667)*IT_0373*IT_0380 + (-2)*IT_0381*(IT_0372
      *IT_0374 + 0.333333333333333*IT_0380*IT_0381))*IT_0383
) : (0
));
    const ccomplex_t IT_0385 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0386 = IT_0005*IT_0385;
    const ccomplex_t IT_0387 = IT_0008*IT_0385;
    const ccomplex_t IT_0388 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0389 = IT_0005*IT_0388;
    const ccomplex_t IT_0390 = IT_0008*IT_0388;
    const ccomplex_t IT_0391 = (0 + _Complex_I*1)*(IT_0386 + IT_0387 + 
      -IT_0389 + -IT_0390);
    const ccomplex_t IT_0392 = (-0.5)*IT_0391;
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0395 = IT_0005*IT_0394;
    const ccomplex_t IT_0396 = IT_0008*IT_0394;
    const ccomplex_t IT_0397 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0398 = IT_0005*IT_0397;
    const ccomplex_t IT_0399 = IT_0008*IT_0397;
    const ccomplex_t IT_0400 = (0 + _Complex_I*1)*(IT_0395 + IT_0396 + 
      -IT_0398 + -IT_0399);
    const ccomplex_t IT_0401 = (-0.5)*IT_0400;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = m_N_1*m_N_2;
    const ccomplex_t IT_0404 = cpow(IT_0013 + IT_0354 + -IT_0374, 2);
    const ccomplex_t IT_0405 = -IT_0354 + 0.25*IT_0001*(IT_0013 + -IT_0354 +
       IT_0374)*(IT_0013 + IT_0354 + -IT_0374) + 0.25*IT_0001*IT_0404;
    const ccomplex_t IT_0406 = 0.5*IT_0013 + 0.5*IT_0354 + (-0.5)*IT_0374;
    const ccomplex_t IT_0407 = 0.5*IT_0013 + (-0.5)*IT_0354 + 0.5*IT_0374;
    const ccomplex_t IT_0408 = IT_0406*IT_0407;
    const ccomplex_t IT_0409 = IT_0013*IT_0405;
    const ccomplex_t IT_0410 = (-0.5)*IT_0409;
    const ccomplex_t IT_0411 = IT_0408 + IT_0410;
    const ccomplex_t IT_0412 = IT_0001*IT_0411;
    const ccomplex_t IT_0413 = IT_0405 + IT_0412;
    const ccomplex_t IT_0414 = cabs(-IT_0354 + 0.25*IT_0001*IT_0404);
    const ccomplex_t IT_0415 = cpow(IT_0414, 0.5);
    const ccomplex_t IT_0416 = ( ((m_Z
) >= (m_N_1 + m_N_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0393)*(IT_0402*IT_0403 + 0.666666666666667*IT_0393
      *IT_0413) + conj(IT_0402)*(IT_0393*IT_0403 + 0.666666666666667*IT_0402
      *IT_0413))*IT_0415
) : (0
));
    const ccomplex_t IT_0417 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0418 = IT_0005*IT_0417;
    const ccomplex_t IT_0419 = IT_0008*IT_0417;
    const ccomplex_t IT_0420 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0421 = IT_0005*IT_0420;
    const ccomplex_t IT_0422 = IT_0008*IT_0420;
    const ccomplex_t IT_0423 = (0 + _Complex_I*1)*(IT_0418 + IT_0419 + 
      -IT_0421 + -IT_0422);
    const ccomplex_t IT_0424 = 0.5*IT_0423;
    const ccomplex_t IT_0425 = cpow(IT_0424, 2);
    const ccomplex_t IT_0426 = pow(m_N_3, 2);
    const ccomplex_t IT_0427 = 0.5*IT_0013 + -IT_0426;
    const ccomplex_t IT_0428 = IT_0013*IT_0427;
    const ccomplex_t IT_0429 = (-0.5)*IT_0428;
    const ccomplex_t IT_0430 = IT_0015 + IT_0429;
    const ccomplex_t IT_0431 = IT_0001*IT_0430;
    const ccomplex_t IT_0432 = IT_0427 + IT_0431;
    const ccomplex_t IT_0433 = (-0.5)*IT_0423;
    const ccomplex_t IT_0434 = cabs(0.25*IT_0013 + -IT_0426);
    const ccomplex_t IT_0435 = cpow(IT_0434, 0.5);
    const ccomplex_t IT_0436 = ( ((m_Z
) >= (2*m_N_3
)) ? (0.125*IT_0000
      *IT_0001*((-0.666666666666667)*IT_0425*IT_0432 + (-2)*IT_0433*(IT_0424
      *IT_0426 + 0.333333333333333*IT_0432*IT_0433))*IT_0435
) : (0
));
    const ccomplex_t IT_0437 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0438 = IT_0005*IT_0437;
    const ccomplex_t IT_0439 = IT_0008*IT_0437;
    const ccomplex_t IT_0440 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0441 = IT_0005*IT_0440;
    const ccomplex_t IT_0442 = IT_0008*IT_0440;
    const ccomplex_t IT_0443 = (0 + _Complex_I*1)*(IT_0438 + IT_0439 + 
      -IT_0441 + -IT_0442);
    const ccomplex_t IT_0444 = (-0.5)*IT_0443;
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0447 = IT_0005*IT_0446;
    const ccomplex_t IT_0448 = IT_0008*IT_0446;
    const ccomplex_t IT_0449 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0450 = IT_0005*IT_0449;
    const ccomplex_t IT_0451 = IT_0008*IT_0449;
    const ccomplex_t IT_0452 = (0 + _Complex_I*1)*(IT_0447 + IT_0448 + 
      -IT_0450 + -IT_0451);
    const ccomplex_t IT_0453 = (-0.5)*IT_0452;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = m_N_1*m_N_3;
    const ccomplex_t IT_0456 = cpow(IT_0013 + IT_0354 + -IT_0426, 2);
    const ccomplex_t IT_0457 = -IT_0354 + 0.25*IT_0001*(IT_0013 + -IT_0354 +
       IT_0426)*(IT_0013 + IT_0354 + -IT_0426) + 0.25*IT_0001*IT_0456;
    const ccomplex_t IT_0458 = 0.5*IT_0013 + 0.5*IT_0354 + (-0.5)*IT_0426;
    const ccomplex_t IT_0459 = 0.5*IT_0013 + (-0.5)*IT_0354 + 0.5*IT_0426;
    const ccomplex_t IT_0460 = IT_0458*IT_0459;
    const ccomplex_t IT_0461 = IT_0013*IT_0457;
    const ccomplex_t IT_0462 = (-0.5)*IT_0461;
    const ccomplex_t IT_0463 = IT_0460 + IT_0462;
    const ccomplex_t IT_0464 = IT_0001*IT_0463;
    const ccomplex_t IT_0465 = IT_0457 + IT_0464;
    const ccomplex_t IT_0466 = cabs(-IT_0354 + 0.25*IT_0001*IT_0456);
    const ccomplex_t IT_0467 = cpow(IT_0466, 0.5);
    const ccomplex_t IT_0468 = ( ((m_Z
) >= (m_N_1 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0445)*(IT_0454*IT_0455 + 0.666666666666667*IT_0445
      *IT_0465) + conj(IT_0454)*(IT_0445*IT_0455 + 0.666666666666667*IT_0454
      *IT_0465))*IT_0467
) : (0
));
    const ccomplex_t IT_0469 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0470 = IT_0005*IT_0469;
    const ccomplex_t IT_0471 = IT_0008*IT_0469;
    const ccomplex_t IT_0472 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0473 = IT_0005*IT_0472;
    const ccomplex_t IT_0474 = IT_0008*IT_0472;
    const ccomplex_t IT_0475 = (0 + _Complex_I*1)*(IT_0470 + IT_0471 + 
      -IT_0473 + -IT_0474);
    const ccomplex_t IT_0476 = (-0.5)*IT_0475;
    const ccomplex_t IT_0477 = -IT_0476;
    const ccomplex_t IT_0478 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0479 = IT_0005*IT_0478;
    const ccomplex_t IT_0480 = IT_0008*IT_0478;
    const ccomplex_t IT_0481 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0482 = IT_0005*IT_0481;
    const ccomplex_t IT_0483 = IT_0008*IT_0481;
    const ccomplex_t IT_0484 = (0 + _Complex_I*1)*(IT_0479 + IT_0480 + 
      -IT_0482 + -IT_0483);
    const ccomplex_t IT_0485 = (-0.5)*IT_0484;
    const ccomplex_t IT_0486 = -IT_0485;
    const ccomplex_t IT_0487 = m_N_2*m_N_3;
    const ccomplex_t IT_0488 = cpow(IT_0013 + IT_0374 + -IT_0426, 2);
    const ccomplex_t IT_0489 = -IT_0374 + 0.25*IT_0001*(IT_0013 + -IT_0374 +
       IT_0426)*(IT_0013 + IT_0374 + -IT_0426) + 0.25*IT_0001*IT_0488;
    const ccomplex_t IT_0490 = 0.5*IT_0013 + 0.5*IT_0374 + (-0.5)*IT_0426;
    const ccomplex_t IT_0491 = 0.5*IT_0013 + (-0.5)*IT_0374 + 0.5*IT_0426;
    const ccomplex_t IT_0492 = IT_0490*IT_0491;
    const ccomplex_t IT_0493 = IT_0013*IT_0489;
    const ccomplex_t IT_0494 = (-0.5)*IT_0493;
    const ccomplex_t IT_0495 = IT_0492 + IT_0494;
    const ccomplex_t IT_0496 = IT_0001*IT_0495;
    const ccomplex_t IT_0497 = IT_0489 + IT_0496;
    const ccomplex_t IT_0498 = cabs(-IT_0374 + 0.25*IT_0001*IT_0488);
    const ccomplex_t IT_0499 = cpow(IT_0498, 0.5);
    const ccomplex_t IT_0500 = ( ((m_Z
) >= (m_N_2 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0477)*(IT_0486*IT_0487 + 0.666666666666667*IT_0477
      *IT_0497) + conj(IT_0486)*(IT_0477*IT_0487 + 0.666666666666667*IT_0486
      *IT_0497))*IT_0499
) : (0
));
    const ccomplex_t IT_0501 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0502 = IT_0005*IT_0501;
    const ccomplex_t IT_0503 = IT_0008*IT_0501;
    const ccomplex_t IT_0504 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0505 = IT_0005*IT_0504;
    const ccomplex_t IT_0506 = IT_0008*IT_0504;
    const ccomplex_t IT_0507 = (0 + _Complex_I*1)*(IT_0502 + IT_0503 + 
      -IT_0505 + -IT_0506);
    const ccomplex_t IT_0508 = 0.5*IT_0507;
    const ccomplex_t IT_0509 = cpow(IT_0508, 2);
    const ccomplex_t IT_0510 = pow(m_N_4, 2);
    const ccomplex_t IT_0511 = 0.5*IT_0013 + -IT_0510;
    const ccomplex_t IT_0512 = IT_0013*IT_0511;
    const ccomplex_t IT_0513 = (-0.5)*IT_0512;
    const ccomplex_t IT_0514 = IT_0015 + IT_0513;
    const ccomplex_t IT_0515 = IT_0001*IT_0514;
    const ccomplex_t IT_0516 = IT_0511 + IT_0515;
    const ccomplex_t IT_0517 = (-0.5)*IT_0507;
    const ccomplex_t IT_0518 = cabs(0.25*IT_0013 + -IT_0510);
    const ccomplex_t IT_0519 = cpow(IT_0518, 0.5);
    const ccomplex_t IT_0520 = ( ((m_Z
) >= (2*m_N_4
)) ? (0.125*IT_0000
      *IT_0001*((-0.666666666666667)*IT_0509*IT_0516 + (-2)*IT_0517*(IT_0508
      *IT_0510 + 0.333333333333333*IT_0516*IT_0517))*IT_0519
) : (0
));
    const ccomplex_t IT_0521 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0522 = IT_0005*IT_0521;
    const ccomplex_t IT_0523 = IT_0008*IT_0521;
    const ccomplex_t IT_0524 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0525 = IT_0005*IT_0524;
    const ccomplex_t IT_0526 = IT_0008*IT_0524;
    const ccomplex_t IT_0527 = (0 + _Complex_I*1)*(IT_0522 + IT_0523 + 
      -IT_0525 + -IT_0526);
    const ccomplex_t IT_0528 = (-0.5)*IT_0527;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0531 = IT_0005*IT_0530;
    const ccomplex_t IT_0532 = IT_0008*IT_0530;
    const ccomplex_t IT_0533 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0534 = IT_0005*IT_0533;
    const ccomplex_t IT_0535 = IT_0008*IT_0533;
    const ccomplex_t IT_0536 = (0 + _Complex_I*1)*(IT_0531 + IT_0532 + 
      -IT_0534 + -IT_0535);
    const ccomplex_t IT_0537 = (-0.5)*IT_0536;
    const ccomplex_t IT_0538 = -IT_0537;
    const ccomplex_t IT_0539 = m_N_1*m_N_4;
    const ccomplex_t IT_0540 = cpow(IT_0013 + IT_0354 + -IT_0510, 2);
    const ccomplex_t IT_0541 = -IT_0354 + 0.25*IT_0001*(IT_0013 + -IT_0354 +
       IT_0510)*(IT_0013 + IT_0354 + -IT_0510) + 0.25*IT_0001*IT_0540;
    const ccomplex_t IT_0542 = 0.5*IT_0013 + 0.5*IT_0354 + (-0.5)*IT_0510;
    const ccomplex_t IT_0543 = 0.5*IT_0013 + (-0.5)*IT_0354 + 0.5*IT_0510;
    const ccomplex_t IT_0544 = IT_0542*IT_0543;
    const ccomplex_t IT_0545 = IT_0013*IT_0541;
    const ccomplex_t IT_0546 = (-0.5)*IT_0545;
    const ccomplex_t IT_0547 = IT_0544 + IT_0546;
    const ccomplex_t IT_0548 = IT_0001*IT_0547;
    const ccomplex_t IT_0549 = IT_0541 + IT_0548;
    const ccomplex_t IT_0550 = cabs(-IT_0354 + 0.25*IT_0001*IT_0540);
    const ccomplex_t IT_0551 = cpow(IT_0550, 0.5);
    const ccomplex_t IT_0552 = ( ((m_Z
) >= (m_N_1 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0529)*(IT_0538*IT_0539 + 0.666666666666667*IT_0529
      *IT_0549) + conj(IT_0538)*(IT_0529*IT_0539 + 0.666666666666667*IT_0538
      *IT_0549))*IT_0551
) : (0
));
    const ccomplex_t IT_0553 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0554 = IT_0005*IT_0553;
    const ccomplex_t IT_0555 = IT_0008*IT_0553;
    const ccomplex_t IT_0556 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0557 = IT_0005*IT_0556;
    const ccomplex_t IT_0558 = IT_0008*IT_0556;
    const ccomplex_t IT_0559 = (0 + _Complex_I*1)*(IT_0554 + IT_0555 + 
      -IT_0557 + -IT_0558);
    const ccomplex_t IT_0560 = (-0.5)*IT_0559;
    const ccomplex_t IT_0561 = -IT_0560;
    const ccomplex_t IT_0562 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0563 = IT_0005*IT_0562;
    const ccomplex_t IT_0564 = IT_0008*IT_0562;
    const ccomplex_t IT_0565 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0566 = IT_0005*IT_0565;
    const ccomplex_t IT_0567 = IT_0008*IT_0565;
    const ccomplex_t IT_0568 = (0 + _Complex_I*1)*(IT_0563 + IT_0564 + 
      -IT_0566 + -IT_0567);
    const ccomplex_t IT_0569 = (-0.5)*IT_0568;
    const ccomplex_t IT_0570 = -IT_0569;
    const ccomplex_t IT_0571 = m_N_2*m_N_4;
    const ccomplex_t IT_0572 = cpow(IT_0013 + IT_0374 + -IT_0510, 2);
    const ccomplex_t IT_0573 = -IT_0374 + 0.25*IT_0001*(IT_0013 + -IT_0374 +
       IT_0510)*(IT_0013 + IT_0374 + -IT_0510) + 0.25*IT_0001*IT_0572;
    const ccomplex_t IT_0574 = 0.5*IT_0013 + 0.5*IT_0374 + (-0.5)*IT_0510;
    const ccomplex_t IT_0575 = 0.5*IT_0013 + (-0.5)*IT_0374 + 0.5*IT_0510;
    const ccomplex_t IT_0576 = IT_0574*IT_0575;
    const ccomplex_t IT_0577 = IT_0013*IT_0573;
    const ccomplex_t IT_0578 = (-0.5)*IT_0577;
    const ccomplex_t IT_0579 = IT_0576 + IT_0578;
    const ccomplex_t IT_0580 = IT_0001*IT_0579;
    const ccomplex_t IT_0581 = IT_0573 + IT_0580;
    const ccomplex_t IT_0582 = cabs(-IT_0374 + 0.25*IT_0001*IT_0572);
    const ccomplex_t IT_0583 = cpow(IT_0582, 0.5);
    const ccomplex_t IT_0584 = ( ((m_Z
) >= (m_N_2 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0561)*(IT_0570*IT_0571 + 0.666666666666667*IT_0561
      *IT_0581) + conj(IT_0570)*(IT_0561*IT_0571 + 0.666666666666667*IT_0570
      *IT_0581))*IT_0583
) : (0
));
    const ccomplex_t IT_0585 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0586 = IT_0005*IT_0585;
    const ccomplex_t IT_0587 = IT_0008*IT_0585;
    const ccomplex_t IT_0588 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0589 = IT_0005*IT_0588;
    const ccomplex_t IT_0590 = IT_0008*IT_0588;
    const ccomplex_t IT_0591 = (0 + _Complex_I*1)*(IT_0586 + IT_0587 + 
      -IT_0589 + -IT_0590);
    const ccomplex_t IT_0592 = (-0.5)*IT_0591;
    const ccomplex_t IT_0593 = -IT_0592;
    const ccomplex_t IT_0594 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0595 = IT_0005*IT_0594;
    const ccomplex_t IT_0596 = IT_0008*IT_0594;
    const ccomplex_t IT_0597 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0598 = IT_0005*IT_0597;
    const ccomplex_t IT_0599 = IT_0008*IT_0597;
    const ccomplex_t IT_0600 = (0 + _Complex_I*1)*(IT_0595 + IT_0596 + 
      -IT_0598 + -IT_0599);
    const ccomplex_t IT_0601 = (-0.5)*IT_0600;
    const ccomplex_t IT_0602 = -IT_0601;
    const ccomplex_t IT_0603 = m_N_3*m_N_4;
    const ccomplex_t IT_0604 = cpow(IT_0013 + IT_0426 + -IT_0510, 2);
    const ccomplex_t IT_0605 = -IT_0426 + 0.25*IT_0001*(IT_0013 + -IT_0426 +
       IT_0510)*(IT_0013 + IT_0426 + -IT_0510) + 0.25*IT_0001*IT_0604;
    const ccomplex_t IT_0606 = 0.5*IT_0013 + 0.5*IT_0426 + (-0.5)*IT_0510;
    const ccomplex_t IT_0607 = 0.5*IT_0013 + (-0.5)*IT_0426 + 0.5*IT_0510;
    const ccomplex_t IT_0608 = IT_0606*IT_0607;
    const ccomplex_t IT_0609 = IT_0013*IT_0605;
    const ccomplex_t IT_0610 = (-0.5)*IT_0609;
    const ccomplex_t IT_0611 = IT_0608 + IT_0610;
    const ccomplex_t IT_0612 = IT_0001*IT_0611;
    const ccomplex_t IT_0613 = IT_0605 + IT_0612;
    const ccomplex_t IT_0614 = cabs(-IT_0426 + 0.25*IT_0001*IT_0604);
    const ccomplex_t IT_0615 = cpow(IT_0614, 0.5);
    const ccomplex_t IT_0616 = ( ((m_Z
) >= (m_N_3 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0593)*(IT_0602*IT_0603 + 0.666666666666667*IT_0593
      *IT_0613) + conj(IT_0602)*(IT_0593*IT_0603 + 0.666666666666667*IT_0602
      *IT_0613))*IT_0615
) : (0
));
    const ccomplex_t IT_0617 = pow(m_tau, 2);
    const ccomplex_t IT_0618 = 0.5*IT_0013 + -IT_0617;
    const ccomplex_t IT_0619 = IT_0013*IT_0618;
    const ccomplex_t IT_0620 = (-0.5)*IT_0619;
    const ccomplex_t IT_0621 = IT_0015 + IT_0620;
    const ccomplex_t IT_0622 = IT_0001*IT_0621;
    const ccomplex_t IT_0623 = IT_0618 + IT_0622;
    const ccomplex_t IT_0624 = cabs(0.25*IT_0013 + -IT_0617);
    const ccomplex_t IT_0625 = cpow(IT_0624, 0.5);
    const ccomplex_t IT_0626 = ( ((m_Z
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0121*IT_0623 + (-4)*IT_0080*(IT_0120
      *IT_0617 + 0.333333333333333*IT_0080*IT_0623))*IT_0625
) : (0
));
    const ccomplex_t IT_0627 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0003
      *IT_0004 + 3*IT_0002*IT_0007)*U_sb_00*conj(U_sb_00) + IT_0003*IT_0004
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0628 = (-0.666666666666667)*IT_0627;
    const ccomplex_t IT_0629 = cpow(IT_0628, 2);
    const ccomplex_t IT_0630 = pow(m_sb_1, 2);
    const ccomplex_t IT_0631 = cabs(0.25*IT_0013 + -IT_0630);
    const ccomplex_t IT_0632 = cpow(IT_0631, 0.5);
    const ccomplex_t IT_0633 = ( ((m_Z
) >= (2*m_sb_1
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0629*(IT_0028 + -IT_0630)*IT_0632
) : (0
));
    const ccomplex_t IT_0634 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0003
      *IT_0004 + 3*IT_0002*IT_0007)*U_sb_01*conj(U_sb_01) + IT_0003*IT_0004
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0635 = (-0.666666666666667)*IT_0634;
    const ccomplex_t IT_0636 = cpow(IT_0635, 2);
    const ccomplex_t IT_0637 = pow(m_sb_2, 2);
    const ccomplex_t IT_0638 = cabs(0.25*IT_0013 + -IT_0637);
    const ccomplex_t IT_0639 = cpow(IT_0638, 0.5);
    const ccomplex_t IT_0640 = ( ((m_Z
) >= (2*m_sb_2
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0636*(IT_0028 + -IT_0637)*IT_0639
) : (0
));
    const ccomplex_t IT_0641 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0003
      *IT_0004 + 3*IT_0002*IT_0007)*U_sb_00*conj(U_sb_01) + IT_0003*IT_0004
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0642 = (-0.666666666666667)*IT_0641;
    const ccomplex_t IT_0643 = 0.5*IT_0013 + 0.5*IT_0630 + (-0.5)*IT_0637;
    const ccomplex_t IT_0644 = cpow(IT_0643, 2);
    const ccomplex_t IT_0645 = IT_0001*IT_0644;
    const ccomplex_t IT_0646 = cpow(IT_0013 + IT_0630 + -IT_0637, 2);
    const ccomplex_t IT_0647 = cabs(-IT_0630 + 0.25*IT_0001*IT_0646);
    const ccomplex_t IT_0648 = cpow(IT_0647, 0.5);
    const ccomplex_t IT_0649 = ( ((m_Z
) >= (m_sb_1 + m_sb_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0642*conj(IT_0642)*(IT_0630 + -IT_0645)*IT_0648
) : 
      (0
));
    const ccomplex_t IT_0650 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0003
      *IT_0004 + 3*IT_0002*IT_0007)*conj(U_sb_00)*U_sb_01 + IT_0003*IT_0004*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0651 = (-0.666666666666667)*IT_0650;
    const ccomplex_t IT_0652 = -IT_0651;
    const ccomplex_t IT_0653 = ( ((m_Z
) >= (m_sb_1 + m_sb_2
)) ? ((-0.125)
      *IT_0000*IT_0001*(IT_0630 + -IT_0645)*IT_0648*IT_0652*conj(IT_0652)
) : 
      (0
));
    const ccomplex_t IT_0654 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0004 + (-3)
      *IT_0002*IT_0007);
    const ccomplex_t IT_0655 = 0.333333333333333*IT_0654;
    const ccomplex_t IT_0656 = cpow(IT_0655, 2);
    const ccomplex_t IT_0657 = pow(m_sc_L, 2);
    const ccomplex_t IT_0658 = cabs(0.25*IT_0013 + -IT_0657);
    const ccomplex_t IT_0659 = cpow(IT_0658, 0.5);
    const ccomplex_t IT_0660 = ( ((m_Z
) >= (2*m_sc_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0656*(IT_0028 + -IT_0657)*IT_0659
) : (0
));
    const ccomplex_t IT_0661 = 1.33333333333333*IT_0080;
    const ccomplex_t IT_0662 = cpow(IT_0661, 2);
    const ccomplex_t IT_0663 = pow(m_sc_R, 2);
    const ccomplex_t IT_0664 = cabs(0.25*IT_0013 + -IT_0663);
    const ccomplex_t IT_0665 = cpow(IT_0664, 0.5);
    const ccomplex_t IT_0666 = ( ((m_Z
) >= (2*m_sc_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0662*(IT_0028 + -IT_0663)*IT_0665
) : (0
));
    const ccomplex_t IT_0667 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0004 + 3
      *IT_0002*IT_0007);
    const ccomplex_t IT_0668 = 0.333333333333333*IT_0667;
    const ccomplex_t IT_0669 = cpow(IT_0668, 2);
    const ccomplex_t IT_0670 = pow(m_sd_L, 2);
    const ccomplex_t IT_0671 = cabs(0.25*IT_0013 + -IT_0670);
    const ccomplex_t IT_0672 = cpow(IT_0671, 0.5);
    const ccomplex_t IT_0673 = ( ((m_Z
) >= (2*m_sd_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0669*(IT_0028 + -IT_0670)*IT_0672
) : (0
));
    const ccomplex_t IT_0674 = pow(m_sd_R, 2);
    const ccomplex_t IT_0675 = cabs(0.25*IT_0013 + -IT_0674);
    const ccomplex_t IT_0676 = cpow(IT_0675, 0.5);
    const ccomplex_t IT_0677 = ( ((m_Z
) >= (2*m_sd_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0096*(IT_0028 + -IT_0674)*IT_0676
) : (0
));
    const ccomplex_t IT_0678 = -IT_0232;
    const ccomplex_t IT_0679 = cpow(IT_0678, 2);
    const ccomplex_t IT_0680 = pow(m_se_L, 2);
    const ccomplex_t IT_0681 = cabs(0.25*IT_0013 + -IT_0680);
    const ccomplex_t IT_0682 = cpow(IT_0681, 0.5);
    const ccomplex_t IT_0683 = ( ((m_Z
) >= (2*m_se_L
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0679*(IT_0028 + -IT_0680)*IT_0682

      ) : (0
));
    const ccomplex_t IT_0684 = (-2)*IT_0080;
    const ccomplex_t IT_0685 = cpow(IT_0684, 2);
    const ccomplex_t IT_0686 = pow(m_se_R, 2);
    const ccomplex_t IT_0687 = cabs(0.25*IT_0013 + -IT_0686);
    const ccomplex_t IT_0688 = cpow(IT_0687, 0.5);
    const ccomplex_t IT_0689 = ( ((m_Z
) >= (2*m_se_R
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0685*(IT_0028 + -IT_0686)*IT_0688

      ) : (0
));
    const ccomplex_t IT_0690 = pow(m_ss_L, 2);
    const ccomplex_t IT_0691 = cabs(0.25*IT_0013 + -IT_0690);
    const ccomplex_t IT_0692 = cpow(IT_0691, 0.5);
    const ccomplex_t IT_0693 = ( ((m_Z
) >= (2*m_ss_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0669*(IT_0028 + -IT_0690)*IT_0692
) : (0
));
    const ccomplex_t IT_0694 = pow(m_ss_R, 2);
    const ccomplex_t IT_0695 = cabs(0.25*IT_0013 + -IT_0694);
    const ccomplex_t IT_0696 = cpow(IT_0695, 0.5);
    const ccomplex_t IT_0697 = ( ((m_Z
) >= (2*m_ss_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0096*(IT_0028 + -IT_0694)*IT_0696
) : (0
));
    const ccomplex_t IT_0698 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0003*IT_0004 
      + (-3)*IT_0002*IT_0007)*U_st_00*conj(U_st_00) + IT_0003*IT_0004*U_st_10
      *conj(U_st_10));
    const ccomplex_t IT_0699 = 1.33333333333333*IT_0698;
    const ccomplex_t IT_0700 = cpow(IT_0699, 2);
    const ccomplex_t IT_0701 = pow(m_st_1, 2);
    const ccomplex_t IT_0702 = cabs(0.25*IT_0013 + -IT_0701);
    const ccomplex_t IT_0703 = cpow(IT_0702, 0.5);
    const ccomplex_t IT_0704 = ( ((m_Z
) >= (2*m_st_1
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0700*(IT_0028 + -IT_0701)*IT_0703
) : (0
));
    const ccomplex_t IT_0705 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0003*IT_0004 
      + (-3)*IT_0002*IT_0007)*U_st_01*conj(U_st_01) + IT_0003*IT_0004*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0706 = 1.33333333333333*IT_0705;
    const ccomplex_t IT_0707 = cpow(IT_0706, 2);
    const ccomplex_t IT_0708 = pow(m_st_2, 2);
    const ccomplex_t IT_0709 = cabs(0.25*IT_0013 + -IT_0708);
    const ccomplex_t IT_0710 = cpow(IT_0709, 0.5);
    const ccomplex_t IT_0711 = ( ((m_Z
) >= (2*m_st_2
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0707*(IT_0028 + -IT_0708)*IT_0710
) : (0
));
    const ccomplex_t IT_0712 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0003*IT_0004 
      + (-3)*IT_0002*IT_0007)*U_st_00*conj(U_st_01) + IT_0003*IT_0004*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0713 = 1.33333333333333*IT_0712;
    const ccomplex_t IT_0714 = 0.5*IT_0013 + 0.5*IT_0701 + (-0.5)*IT_0708;
    const ccomplex_t IT_0715 = cpow(IT_0714, 2);
    const ccomplex_t IT_0716 = IT_0001*IT_0715;
    const ccomplex_t IT_0717 = cpow(IT_0013 + IT_0701 + -IT_0708, 2);
    const ccomplex_t IT_0718 = cabs(-IT_0701 + 0.25*IT_0001*IT_0717);
    const ccomplex_t IT_0719 = cpow(IT_0718, 0.5);
    const ccomplex_t IT_0720 = ( ((m_Z
) >= (m_st_1 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0713*conj(IT_0713)*(IT_0701 + -IT_0716)*IT_0719
) : 
      (0
));
    const ccomplex_t IT_0721 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0003*IT_0004 
      + (-3)*IT_0002*IT_0007)*conj(U_st_00)*U_st_01 + IT_0003*IT_0004*conj
      (U_st_10)*U_st_11);
    const ccomplex_t IT_0722 = 1.33333333333333*IT_0721;
    const ccomplex_t IT_0723 = -IT_0722;
    const ccomplex_t IT_0724 = ( ((m_Z
) >= (m_st_1 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*(IT_0701 + -IT_0716)*IT_0719*IT_0723*conj(IT_0723)
) : 
      (0
));
    const ccomplex_t IT_0725 = pow(m_su_L, 2);
    const ccomplex_t IT_0726 = cabs(0.25*IT_0013 + -IT_0725);
    const ccomplex_t IT_0727 = cpow(IT_0726, 0.5);
    const ccomplex_t IT_0728 = ( ((m_Z
) >= (2*m_su_L
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0656*(IT_0028 + -IT_0725)*IT_0727
) : (0
));
    const ccomplex_t IT_0729 = pow(m_su_R, 2);
    const ccomplex_t IT_0730 = cabs(0.25*IT_0013 + -IT_0729);
    const ccomplex_t IT_0731 = cpow(IT_0730, 0.5);
    const ccomplex_t IT_0732 = ( ((m_Z
) >= (2*m_su_R
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0662*(IT_0028 + -IT_0729)*IT_0731
) : (0
));
    const ccomplex_t IT_0733 = pow(m_smu_L, 2);
    const ccomplex_t IT_0734 = cabs(0.25*IT_0013 + -IT_0733);
    const ccomplex_t IT_0735 = cpow(IT_0734, 0.5);
    const ccomplex_t IT_0736 = ( ((m_Z
) >= (2*m_smu_L
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0679*(IT_0028 + -IT_0733)*IT_0735

      ) : (0
));
    const ccomplex_t IT_0737 = pow(m_smu_R, 2);
    const ccomplex_t IT_0738 = cabs(0.25*IT_0013 + -IT_0737);
    const ccomplex_t IT_0739 = cpow(IT_0738, 0.5);
    const ccomplex_t IT_0740 = ( ((m_Z
) >= (2*m_smu_R
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0685*(IT_0028 + -IT_0737)*IT_0739

      ) : (0
));
    const ccomplex_t IT_0741 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0007;
    const ccomplex_t IT_0742 = -IT_0741;
    const ccomplex_t IT_0743 = cpow(IT_0742, 2);
    const ccomplex_t IT_0744 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0745 = cabs(0.25*IT_0013 + -IT_0744);
    const ccomplex_t IT_0746 = cpow(IT_0745, 0.5);
    const ccomplex_t IT_0747 = ( ((m_Z
) >= (2*m_snu_mu
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0743*(IT_0028 + -IT_0744)*IT_0746

      ) : (0
));
    const ccomplex_t IT_0748 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0003*IT_0004 +
       -IT_0002*IT_0007)*U_stau_00*conj(U_stau_00) + IT_0003*IT_0004*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0749 = (-2)*IT_0748;
    const ccomplex_t IT_0750 = cpow(IT_0749, 2);
    const ccomplex_t IT_0751 = pow(m_stau_1, 2);
    const ccomplex_t IT_0752 = cabs(0.25*IT_0013 + -IT_0751);
    const ccomplex_t IT_0753 = cpow(IT_0752, 0.5);
    const ccomplex_t IT_0754 = ( ((m_Z
) >= (2*m_stau_1
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0750*(IT_0028 + -IT_0751)*IT_0753

      ) : (0
));
    const ccomplex_t IT_0755 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0003*IT_0004 +
       -IT_0002*IT_0007)*U_stau_01*conj(U_stau_01) + IT_0003*IT_0004*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0756 = (-2)*IT_0755;
    const ccomplex_t IT_0757 = cpow(IT_0756, 2);
    const ccomplex_t IT_0758 = pow(m_stau_2, 2);
    const ccomplex_t IT_0759 = cabs(0.25*IT_0013 + -IT_0758);
    const ccomplex_t IT_0760 = cpow(IT_0759, 0.5);
    const ccomplex_t IT_0761 = ( ((m_Z
) >= (2*m_stau_2
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0757*(IT_0028 + -IT_0758)*IT_0760

      ) : (0
));
    const ccomplex_t IT_0762 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0003*IT_0004 +
       -IT_0002*IT_0007)*U_stau_00*conj(U_stau_01) + IT_0003*IT_0004*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0763 = (-2)*IT_0762;
    const ccomplex_t IT_0764 = 0.5*IT_0013 + 0.5*IT_0751 + (-0.5)*IT_0758;
    const ccomplex_t IT_0765 = cpow(IT_0764, 2);
    const ccomplex_t IT_0766 = IT_0001*IT_0765;
    const ccomplex_t IT_0767 = cpow(IT_0013 + IT_0751 + -IT_0758, 2);
    const ccomplex_t IT_0768 = cabs(-IT_0751 + 0.25*IT_0001*IT_0767);
    const ccomplex_t IT_0769 = cpow(IT_0768, 0.5);
    const ccomplex_t IT_0770 = ( ((m_Z
) >= (m_stau_1 + m_stau_2
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0763*conj(IT_0763)*(IT_0751 + 
      -IT_0766)*IT_0769
) : (0
));
    const ccomplex_t IT_0771 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0003*IT_0004 +
       -IT_0002*IT_0007)*conj(U_stau_00)*U_stau_01 + IT_0003*IT_0004*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0772 = (-2)*IT_0771;
    const ccomplex_t IT_0773 = -IT_0772;
    const ccomplex_t IT_0774 = ( ((m_Z
) >= (m_stau_1 + m_stau_2
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*(IT_0751 + -IT_0766)*IT_0769*IT_0773
      *conj(IT_0773)
) : (0
));
    const ccomplex_t IT_0775 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0776 = cabs(0.25*IT_0013 + -IT_0775);
    const ccomplex_t IT_0777 = cpow(IT_0776, 0.5);
    const ccomplex_t IT_0778 = ( ((m_Z
) >= (2*m_snu_tau
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0743*(IT_0028 + -IT_0775)*IT_0777

      ) : (0
));
    const ccomplex_t IT_0779 = 3*IT_0022 + IT_0079 + IT_0094 + IT_0108 +
       IT_0118 + IT_0131 + IT_0161 + IT_0171 + IT_0181 + IT_0192 + IT_0203 +
       IT_0213 + IT_0231 + IT_0237 + IT_0247 + IT_0273 + IT_0299 + IT_0328 +
       IT_0344 + IT_0364 + IT_0384 + IT_0416 + IT_0436 + IT_0468 + IT_0500 +
       IT_0520 + IT_0552 + IT_0584 + IT_0616 + IT_0626 + IT_0633 + IT_0640 +
       IT_0649 + IT_0653 + IT_0660 + IT_0666 + IT_0673 + IT_0677 + IT_0683 +
       IT_0689 + IT_0693 + IT_0697 + IT_0704 + IT_0711 + IT_0720 + IT_0724 +
       IT_0728 + IT_0732 + IT_0736 + IT_0740 + 2*IT_0747 + IT_0754 + IT_0761 +
       IT_0770 + IT_0774 + IT_0778;
    return create_ccomplex_return(IT_0779);
}

