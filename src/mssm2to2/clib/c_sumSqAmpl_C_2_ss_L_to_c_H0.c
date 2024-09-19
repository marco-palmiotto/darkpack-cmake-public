#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_ss_L_to_c_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_ss_L_to_c_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_s, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0004*IT_0006 + (-0.5)*m_W*(IT_0002*IT_0004 + -IT_0007*IT_0008)
      *(IT_0006 + 1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0007, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = pow(m_c, 2);
    const ccomplex_t IT_0018 = pow(m_C_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = U_d1*V_Wp2;
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = V_u2*U_Wm1;
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = IT_0006*(IT_0026 + IT_0028);
    const ccomplex_t IT_0030 = 1.4142135623731*e_em*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0024*IT_0032;
    const ccomplex_t IT_0034 = pow(m_ss_L, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0017 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0022*IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0038 = U_d2*V_Wp2;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = V_u2*U_Wm2;
    const ccomplex_t IT_0041 = IT_0008*IT_0040;
    const ccomplex_t IT_0042 = IT_0006*(IT_0039 + IT_0041);
    const ccomplex_t IT_0043 = 1.4142135623731*e_em*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0015*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0017 + IT_0034 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0037*IT_0046*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + IT_0017 + -IT_0018 + -IT_0034 
      + -reg_prop, -1);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0006
      *IT_0008*IT_0013;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0015*IT_0051;
    const ccomplex_t IT_0053 = IT_0049*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = m_c*IT_0054;
    const ccomplex_t IT_0056 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0059 = IT_0008*IT_0058;
    const ccomplex_t IT_0060 = IT_0006*(IT_0057 + IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = IT_0024*IT_0063;
    const ccomplex_t IT_0065 = IT_0035*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = m_C_2*IT_0066;
    const ccomplex_t IT_0068 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = IT_0006*(IT_0069 + IT_0071);
    const ccomplex_t IT_0073 = 1.4142135623731*e_em*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0015*IT_0075;
    const ccomplex_t IT_0077 = IT_0047*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = m_C_2*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0000*IT_0006*IT_0013*U_sb_01;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0004*conj(U_sb_11);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = IT_0080*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0000*IT_0006*IT_0013*U_sb_00;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0004*conj(U_sb_10);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = IT_0088*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0021 + -IT_0036 + -IT_0048 + (-2)*IT_0055 + 
      -IT_0067 + -IT_0079 + IT_0087 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0004;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_cd*e_em*IT_0000*IT_0003*IT_0006;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0003*IT_0006;
    const ccomplex_t IT_0106 = m_b*conj(U_d2)*V_cb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = V_cb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0110 = IT_0006*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0108 + (-2)*IT_0110);
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0084*IT_0112;
    const ccomplex_t IT_0114 = IT_0080*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0117 = IT_0063*IT_0116;
    const ccomplex_t IT_0118 = IT_0022*IT_0035*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0120 = IT_0075*IT_0119;
    const ccomplex_t IT_0121 = IT_0037*IT_0047*IT_0120;
    const ccomplex_t IT_0122 = IT_0045*IT_0119;
    const ccomplex_t IT_0123 = IT_0047*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = m_C_2*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0004 + -m_s*e_em*mu_h*IT_0000*IT_0003*IT_0006*IT_0008);
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0000*IT_0003*IT_0006;
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = IT_0126*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0012*IT_0119;
    const ccomplex_t IT_0135 = IT_0019*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = m_b*conj(U_d2)*V_cb*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0138 = IT_0105*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = V_cb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0141 = IT_0006*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0139 + (-2)*IT_0141);
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0092*IT_0143;
    const ccomplex_t IT_0145 = IT_0088*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0051*IT_0119;
    const ccomplex_t IT_0148 = IT_0049*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = m_c*IT_0149;
    const ccomplex_t IT_0151 = IT_0032*IT_0116;
    const ccomplex_t IT_0152 = IT_0035*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = m_C_2*IT_0153;
    const ccomplex_t IT_0155 = -IT_0104 + -IT_0115 + IT_0118 + IT_0121 +
       IT_0125 + -IT_0133 + -IT_0136 + -IT_0146 + 2*IT_0150 + IT_0154;
    const ccomplex_t IT_0156 = 6*s_13;
    const ccomplex_t IT_0157 = m_c*m_C_2;
    const ccomplex_t IT_0158 = 6*IT_0096;
    const ccomplex_t IT_0159 = -IT_0054 + -IT_0124 + -IT_0153;
    const ccomplex_t IT_0160 = m_c*s_14;
    const ccomplex_t IT_0161 = IT_0066 + IT_0078 + IT_0149;
    const ccomplex_t IT_0162 = s_34*m_C_2;
    const ccomplex_t IT_0163 = 6*conj(IT_0096);
    const ccomplex_t IT_0164 = 6*conj(IT_0155);
    const ccomplex_t IT_0165 = pow(m_H0, 2);
    const ccomplex_t IT_0166 = m_c*m_C_2*IT_0165;
    const ccomplex_t IT_0167 = s_13*IT_0165;
    const ccomplex_t IT_0168 = s_14*s_34;
    const ccomplex_t IT_0169 = (-2)*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = 6*IT_0155;
    const ccomplex_t IT_0172 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0173 = (IT_0096*conj(IT_0096) + IT_0155*conj(IT_0155))
      *IT_0156 + conj(IT_0155)*IT_0157*IT_0158 + (IT_0155*IT_0157 + IT_0159
      *IT_0160 + IT_0161*IT_0162)*IT_0163 + (IT_0160*IT_0161 + IT_0159*IT_0162)
      *IT_0164 + 6*conj(IT_0161)*(IT_0159*IT_0166 + -IT_0161*IT_0170 +
       0.166666666666667*IT_0160*IT_0171 + IT_0158*IT_0172) + 6*conj(IT_0159)*
      (0.166666666666667*IT_0158*IT_0160 + IT_0161*IT_0166 + -IT_0159*IT_0170 +
       IT_0171*IT_0172);
    return create_ccomplex_return(IT_0173);
}

