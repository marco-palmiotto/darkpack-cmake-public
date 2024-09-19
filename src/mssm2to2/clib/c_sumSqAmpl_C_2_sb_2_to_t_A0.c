#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_t_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_t_A0(
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = V_u2*U_Wm1;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = U_d1*V_Wp2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0000, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = m_b*conj(U_d1)*V_tb*e_em*IT_0015*U_sb_11;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0012 + (-0.5)*IT_0018);
    const ccomplex_t IT_0020 = IT_0010*IT_0019;
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = pow(m_sb_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0027 = V_u2*U_Wm2;
    const ccomplex_t IT_0028 = IT_0000*IT_0027;
    const ccomplex_t IT_0029 = U_d2*V_Wp2;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = IT_0007*(IT_0028 + IT_0030);
    const ccomplex_t IT_0032 = 1.4142135623731*e_em*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = m_b*conj(U_d2)*V_tb*e_em*IT_0015*U_sb_11;
    const ccomplex_t IT_0035 = IT_0014*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0038 = IT_0007*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0036 + (-2)*IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0033*IT_0040;
    const ccomplex_t IT_0042 = IT_0026*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = pow(m_C_2, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + -IT_0022 + -IT_0044 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0046 = cpow(IT_0003, -1);
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0007*IT_0015*IT_0046*U_sb_01;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = m_t*e_em*IT_0000*IT_0007*IT_0015*IT_0046;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0025 + -IT_0043 + -IT_0053;
    const ccomplex_t IT_0055 = pow(m_A0, 2);
    const ccomplex_t IT_0056 = m_t*m_C_2*IT_0055;
    const ccomplex_t IT_0057 = IT_0040*IT_0050;
    const ccomplex_t IT_0058 = IT_0045*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0061 = IT_0000*IT_0060;
    const ccomplex_t IT_0062 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = IT_0007*(IT_0061 + IT_0063);
    const ccomplex_t IT_0065 = 1.4142135623731*e_em*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0007*IT_0015*IT_0046*U_sb_01;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0023*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0073 = IT_0000*IT_0072;
    const ccomplex_t IT_0074 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0075 = IT_0003*IT_0074;
    const ccomplex_t IT_0076 = IT_0007*(IT_0073 + IT_0075);
    const ccomplex_t IT_0077 = 1.4142135623731*e_em*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0048*IT_0078;
    const ccomplex_t IT_0080 = IT_0026*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0059 + IT_0071 + IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u2
      *e_em*IT_0007*IT_0015*IT_0046;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_11;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0083*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0092 = m_b*e_em*mu_h*IT_0015*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = m_b*e_em*mu_h*IT_0015*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0095 = IT_0007*IT_0094;
    const ccomplex_t IT_0096 = cpow(V_tb, 2);
    const ccomplex_t IT_0097 = A_b*IT_0096*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = A_b*IT_0096*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0101 = IT_0003*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0093 + 0.5*IT_0095 + (-0.5)*IT_0099 +
       0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0048*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0007*IT_0015*IT_0046;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0003
      *U_sb_11;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = IT_0107*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0116 = m_b*e_em*mu_h*IT_0015*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0117 = IT_0007*IT_0116;
    const ccomplex_t IT_0118 = m_b*e_em*mu_h*IT_0015*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0119 = IT_0007*IT_0118;
    const ccomplex_t IT_0120 = A_b*IT_0096*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0121 = IT_0003*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = A_b*IT_0096*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0124 = IT_0003*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = 0.5*IT_0117 + (-0.5)*IT_0119 + 0.5*IT_0122 + (
      -0.5)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0007*IT_0015*IT_0046*U_sb_00;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0126*IT_0128;
    const ccomplex_t IT_0130 = IT_0115*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0133 = IT_0033*IT_0048;
    const ccomplex_t IT_0134 = IT_0026*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0136 = IT_0010*IT_0068;
    const ccomplex_t IT_0137 = IT_0023*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = m_C_2*IT_0071;
    const ccomplex_t IT_0139 = m_C_2*IT_0081;
    const ccomplex_t IT_0140 = IT_0090 + IT_0106 + IT_0114 + IT_0131 + 
      -IT_0134 + -IT_0137 + -IT_0138 + -IT_0139;
    const ccomplex_t IT_0141 = m_t*s_14;
    const ccomplex_t IT_0142 = 6*IT_0141;
    const ccomplex_t IT_0143 = s_34*m_C_2;
    const ccomplex_t IT_0144 = 6*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_01;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0003
      *U_sb_01;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0154*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_13 + IT_0021 + IT_0044 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0040*IT_0103;
    const ccomplex_t IT_0162 = IT_0091*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(U_d2)*V_tb*e_em*IT_0015*U_sb_10;
    const ccomplex_t IT_0165 = IT_0014*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0168 = IT_0007*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0166 + (-2)*IT_0168);
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0126*IT_0170;
    const ccomplex_t IT_0172 = IT_0115*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0040*IT_0078;
    const ccomplex_t IT_0175 = IT_0026*IT_0132*IT_0174;
    const ccomplex_t IT_0176 = m_C_2*IT_0025;
    const ccomplex_t IT_0177 = m_C_2*IT_0043;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0179 = IT_0111*IT_0178;
    const ccomplex_t IT_0180 = IT_0107*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0183 = IT_0087*IT_0182;
    const ccomplex_t IT_0184 = IT_0083*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0019*IT_0066;
    const ccomplex_t IT_0187 = IT_0023*IT_0135*IT_0186;
    const ccomplex_t IT_0188 = -IT_0152 + -IT_0160 + -IT_0163 + -IT_0173 +
       IT_0175 + IT_0176 + IT_0177 + -IT_0181 + -IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = 0.166666666666667*conj(IT_0188);
    const ccomplex_t IT_0190 = s_13*IT_0055;
    const ccomplex_t IT_0191 = s_14*s_34;
    const ccomplex_t IT_0192 = (-2)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = (-6)*IT_0193;
    const ccomplex_t IT_0195 = 0.166666666666667*conj(IT_0054);
    const ccomplex_t IT_0196 = m_t*m_C_2;
    const ccomplex_t IT_0197 = 6*IT_0196;
    const ccomplex_t IT_0198 = 0.166666666666667*conj(IT_0082);
    const ccomplex_t IT_0199 = 6*IT_0054*(IT_0056*conj(IT_0082) +
       0.166666666666667*conj(IT_0140)*IT_0142 + IT_0144*IT_0189 + IT_0194
      *IT_0195) + 6*IT_0188*(s_13*conj(IT_0188) + IT_0144*IT_0195 +
       0.166666666666667*conj(IT_0140)*IT_0197 + IT_0142*IT_0198) + 6*IT_0140*
      (s_13*conj(IT_0140) + IT_0142*IT_0195 + IT_0189*IT_0197 + IT_0144*IT_0198)
       + 6*IT_0082*(conj(IT_0054)*IT_0056 + 0.166666666666667*conj(IT_0140)
      *IT_0144 + IT_0142*IT_0189 + IT_0194*IT_0198);
    return create_ccomplex_return(IT_0199);
}

