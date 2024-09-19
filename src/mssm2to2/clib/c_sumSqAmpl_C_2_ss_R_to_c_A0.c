#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_ss_R_to_c_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_ss_R_to_c_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0000 = sin(beta);
    const ccomplex_t IT_0001 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0000;
    const ccomplex_t IT_0002 = (-0.5)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cpow(IT_0000, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0003*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = pow(m_c, 2);
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = m_s*e_em*mu_h*IT_0003;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0018 = IT_0000*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0016 + (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0003*IT_0004*IT_0006;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0003*IT_0004*IT_0006*U_sb_01;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*A_b*V_tb*V_ts*IT_0000*conj
      (U_sb_01);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0003*IT_0004*IT_0006*U_sb_00;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*A_b*V_tb*V_ts*IT_0000*conj
      (U_sb_00);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = IT_0035*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0014 + IT_0026 + IT_0034 + IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*V_cd*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0045 = IT_0002*IT_0044;
    const ccomplex_t IT_0046 = IT_0012*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0049 = IT_0020*IT_0048;
    const ccomplex_t IT_0050 = IT_0024*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cos(beta);
    const ccomplex_t IT_0053 = cpow(IT_0052, -1);
    const ccomplex_t IT_0054 = IT_0006*IT_0053;
    const ccomplex_t IT_0055 = m_b*conj(U_d2)*V_cb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = V_cb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0059 = IT_0006*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0057 + (-2)*IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0031*IT_0061;
    const ccomplex_t IT_0063 = IT_0027*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0066 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0067 = IT_0052*IT_0066;
    const ccomplex_t IT_0068 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0069 = IT_0000*IT_0068;
    const ccomplex_t IT_0070 = IT_0006*(IT_0067 + IT_0069);
    const ccomplex_t IT_0071 = 1.4142135623731*e_em*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0003*IT_0006*IT_0053;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = pow(m_ss_R, 2);
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0010 + IT_0076 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0065*IT_0075*IT_0077;
    const ccomplex_t IT_0079 = V_u2*U_Wm2;
    const ccomplex_t IT_0080 = IT_0052*IT_0079;
    const ccomplex_t IT_0081 = U_d2*V_Wp2;
    const ccomplex_t IT_0082 = IT_0000*IT_0081;
    const ccomplex_t IT_0083 = IT_0006*(IT_0080 + IT_0082);
    const ccomplex_t IT_0084 = 1.4142135623731*e_em*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = IT_0077*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = m_C_2*IT_0088;
    const ccomplex_t IT_0090 = m_b*conj(U_d2)*V_cb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0091 = IT_0054*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = V_cb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0094 = IT_0006*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0092 + (-2)*IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0039*IT_0096;
    const ccomplex_t IT_0098 = IT_0035*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0101 = cpow((-2)*s_23 + IT_0010 + IT_0076 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0102 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0103 = IT_0052*IT_0102;
    const ccomplex_t IT_0104 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0105 = IT_0000*IT_0104;
    const ccomplex_t IT_0106 = IT_0006*(IT_0103 + IT_0105);
    const ccomplex_t IT_0107 = 1.4142135623731*e_em*IT_0106;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0003*IT_0006*IT_0053;
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0101*IT_0111;
    const ccomplex_t IT_0113 = V_u2*U_Wm1;
    const ccomplex_t IT_0114 = IT_0052*IT_0113;
    const ccomplex_t IT_0115 = U_d1*V_Wp2;
    const ccomplex_t IT_0116 = IT_0000*IT_0115;
    const ccomplex_t IT_0117 = IT_0006*(IT_0114 + IT_0116);
    const ccomplex_t IT_0118 = 1.4142135623731*e_em*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0110*IT_0119;
    const ccomplex_t IT_0121 = IT_0101*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = m_C_2*IT_0122;
    const ccomplex_t IT_0124 = -IT_0047 + -IT_0051 + -IT_0064 + IT_0078 +
       IT_0089 + -IT_0099 + IT_0112 + IT_0123;
    const ccomplex_t IT_0125 = 6*s_13;
    const ccomplex_t IT_0126 = m_c*m_C_2;
    const ccomplex_t IT_0127 = 6*IT_0043;
    const ccomplex_t IT_0128 = -IT_0088 + -IT_0122;
    const ccomplex_t IT_0129 = m_c*s_14;
    const ccomplex_t IT_0130 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0076 
      + -reg_prop, -1);
    const ccomplex_t IT_0131 = m_c*e_em*IT_0003*IT_0004*IT_0006*IT_0052;
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0074*IT_0132;
    const ccomplex_t IT_0134 = IT_0130*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_34*m_C_2;
    const ccomplex_t IT_0138 = 6*conj(IT_0043);
    const ccomplex_t IT_0139 = 6*conj(IT_0124);
    const ccomplex_t IT_0140 = pow(m_A0, 2);
    const ccomplex_t IT_0141 = m_c*m_C_2*IT_0140;
    const ccomplex_t IT_0142 = s_13*IT_0140;
    const ccomplex_t IT_0143 = s_14*s_34;
    const ccomplex_t IT_0144 = (-2)*IT_0143;
    const ccomplex_t IT_0145 = IT_0142 + IT_0144;
    const ccomplex_t IT_0146 = 6*IT_0124;
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0137;
    const ccomplex_t IT_0148 = (IT_0043*conj(IT_0043) + IT_0124*conj(IT_0124))
      *IT_0125 + conj(IT_0124)*IT_0126*IT_0127 + (IT_0124*IT_0126 + IT_0128
      *IT_0129 + IT_0136*IT_0137)*IT_0138 + (IT_0129*IT_0136 + IT_0128*IT_0137)
      *IT_0139 + 6*conj(IT_0136)*(IT_0128*IT_0141 + -IT_0136*IT_0145 +
       0.166666666666667*IT_0129*IT_0146 + IT_0127*IT_0147) + 6*conj(IT_0128)*
      (0.166666666666667*IT_0127*IT_0129 + IT_0136*IT_0141 + -IT_0128*IT_0145 +
       IT_0146*IT_0147);
    return create_ccomplex_return(IT_0148);
}

