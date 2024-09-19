#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sb_2_to_t_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sb_2_to_t_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0001
      *U_sb_11;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_t, 2);
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = V_u1*U_Wm1;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = U_d1*V_Wp1;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = IT_0004*(IT_0018 + IT_0020);
    const ccomplex_t IT_0022 = 1.4142135623731*e_em*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0000*IT_0002*IT_0004*U_sb_01;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = pow(m_sb_2, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0010 + IT_0027 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0015*IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0031 = V_u1*U_Wm2;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = U_d2*V_Wp1;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = IT_0004*(IT_0032 + IT_0034);
    const ccomplex_t IT_0036 = 1.4142135623731*e_em*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0000*IT_0002*IT_0004*U_sb_01;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0010 + IT_0027 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0030*IT_0040*IT_0041;
    const ccomplex_t IT_0043 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0044 = IT_0016*IT_0043;
    const ccomplex_t IT_0045 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = IT_0004*(IT_0044 + IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = IT_0025*IT_0049;
    const ccomplex_t IT_0051 = IT_0028*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = m_C_1*IT_0052;
    const ccomplex_t IT_0054 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0055 = IT_0016*IT_0054;
    const ccomplex_t IT_0056 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0057 = IT_0001*IT_0056;
    const ccomplex_t IT_0058 = IT_0004*(IT_0055 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0039*IT_0060;
    const ccomplex_t IT_0062 = IT_0041*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = m_C_1*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0066 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = m_b*e_em*mu_h*IT_0000*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = cpow(V_tb, 2);
    const ccomplex_t IT_0071 = A_b*IT_0070*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0072 = IT_0001*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = A_b*IT_0070*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0075 = IT_0001*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0067 + (-0.5)*IT_0069 + 0.5*IT_0073 + (
      -0.5)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0000*IT_0002*IT_0004*U_sb_00;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = IT_0065*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u1
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0001
      *U_sb_11;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0092 = IT_0004*IT_0091;
    const ccomplex_t IT_0093 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = A_b*IT_0070*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0096 = IT_0001*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = A_b*IT_0070*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0099 = IT_0001*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (-0.5)*IT_0092 + 0.5*IT_0094 + (-0.5)*IT_0097 +
       0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0025*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0014 + -IT_0029 + -IT_0042 + -IT_0053 + 
      -IT_0064 + IT_0082 + IT_0090 + IT_0105;
    const ccomplex_t IT_0107 = 6*s_13;
    const ccomplex_t IT_0108 = cpow((-2)*s_12 + -IT_0011 + -IT_0027 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0109 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0110 = IT_0004*IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0016, -1);
    const ccomplex_t IT_0112 = IT_0004*IT_0111;
    const ccomplex_t IT_0113 = m_b*conj(U_d1)*V_tb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0110 + (-0.5)*IT_0115);
    const ccomplex_t IT_0117 = m_t*e_em*IT_0000*IT_0002*IT_0004*IT_0016;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = IT_0108*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = IT_0052 + IT_0063 + -IT_0121;
    const ccomplex_t IT_0123 = m_t*s_14;
    const ccomplex_t IT_0124 = IT_0023*IT_0116;
    const ccomplex_t IT_0125 = IT_0028*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0128 = IT_0004*IT_0127;
    const ccomplex_t IT_0129 = m_b*conj(U_d2)*V_tb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0130 = IT_0112*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + (-0.5)*IT_0131);
    const ccomplex_t IT_0133 = IT_0037*IT_0132;
    const ccomplex_t IT_0134 = IT_0041*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0025*IT_0118;
    const ccomplex_t IT_0137 = IT_0108*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = -IT_0126 + -IT_0135 + -IT_0138;
    const ccomplex_t IT_0140 = s_34*m_C_1;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0000*IT_0004*IT_0111;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0001
      *U_sb_01;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = IT_0142*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0004*IT_0111;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0001
      *U_sb_01;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = IT_0150*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0060*IT_0132;
    const ccomplex_t IT_0158 = IT_0030*IT_0041*IT_0157;
    const ccomplex_t IT_0159 = m_C_1*IT_0126;
    const ccomplex_t IT_0160 = m_C_1*IT_0135;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0162 = IT_0086*IT_0161;
    const ccomplex_t IT_0163 = IT_0088*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0101*IT_0116;
    const ccomplex_t IT_0166 = IT_0103*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0169 = IT_0008*IT_0168;
    const ccomplex_t IT_0170 = IT_0012*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0173 = IT_0004*IT_0172;
    const ccomplex_t IT_0174 = m_b*conj(U_d1)*V_tb*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0175 = IT_0112*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + (-0.5)*IT_0176);
    const ccomplex_t IT_0178 = IT_0077*IT_0177;
    const ccomplex_t IT_0179 = IT_0065*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0049*IT_0116;
    const ccomplex_t IT_0182 = IT_0015*IT_0028*IT_0181;
    const ccomplex_t IT_0183 = -IT_0148 + -IT_0156 + IT_0158 + IT_0159 +
       IT_0160 + -IT_0164 + -IT_0167 + -IT_0171 + -IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = 6*conj(IT_0183);
    const ccomplex_t IT_0185 = m_t*m_C_1;
    const ccomplex_t IT_0186 = 6*IT_0106;
    const ccomplex_t IT_0187 = 6*conj(IT_0106);
    const ccomplex_t IT_0188 = pow(m_A0, 2);
    const ccomplex_t IT_0189 = m_t*m_C_1*IT_0188;
    const ccomplex_t IT_0190 = s_13*IT_0188;
    const ccomplex_t IT_0191 = s_14*s_34;
    const ccomplex_t IT_0192 = (-2)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = 6*IT_0183;
    const ccomplex_t IT_0195 = 0.166666666666667*IT_0140;
    const ccomplex_t IT_0196 = IT_0106*conj(IT_0106)*IT_0107 + (IT_0122
      *IT_0123 + IT_0139*IT_0140)*IT_0184 + conj(IT_0183)*(IT_0107*IT_0183 +
       IT_0185*IT_0186) + (IT_0123*IT_0139 + IT_0122*IT_0140 + IT_0183*IT_0185)
      *IT_0187 + 6*conj(IT_0122)*(IT_0139*IT_0189 + -IT_0122*IT_0193 +
       0.166666666666667*IT_0123*IT_0194 + IT_0186*IT_0195) + 6*conj(IT_0139)*
      (0.166666666666667*IT_0123*IT_0186 + IT_0122*IT_0189 + -IT_0139*IT_0193 +
       IT_0194*IT_0195);
    return create_ccomplex_return(IT_0196);
}

