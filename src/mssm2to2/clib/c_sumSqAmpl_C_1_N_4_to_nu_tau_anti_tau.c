#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_4_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_4_to_nu_tau_anti_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + (-0.5)*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = pow(m_C_1, 2);
    const ccomplex_t IT_0014 = pow(m_N_4, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0013 + -IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cpow(s_23 + (-0.5)*IT_0014 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = conj(N_B4)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*e_em;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + -IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0002;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = IT_0018*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = 8*IT_0017;
    const ccomplex_t IT_0036 = s_23*m_C_1*m_tau;
    const ccomplex_t IT_0037 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0038 = IT_0002*IT_0037;
    const ccomplex_t IT_0039 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0040 = IT_0002*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 0.5*IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = IT_0015*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = pow(m_W, -1);
    const ccomplex_t IT_0048 = cos(beta);
    const ccomplex_t IT_0049 = cpow(IT_0048, -1);
    const ccomplex_t IT_0050 = sin(beta);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0002*IT_0047*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0002*IT_0048;
    const ccomplex_t IT_0054 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = 2*IT_0055;
    const ccomplex_t IT_0057 = IT_0020*IT_0048;
    const ccomplex_t IT_0058 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0060 + IT_0063;
    const ccomplex_t IT_0065 = IT_0056 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0052*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + -IT_0013 + -IT_0014 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = pow(m_W, 2);
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + -IT_0013 + -IT_0014 + IT_0072 
      + -reg_prop, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0002*IT_0047;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0002*IT_0050;
    const ccomplex_t IT_0077 = IT_0054*IT_0076;
    const ccomplex_t IT_0078 = 2*IT_0077;
    const ccomplex_t IT_0079 = IT_0020*IT_0050;
    const ccomplex_t IT_0080 = IT_0058*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0061*IT_0076;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = IT_0081 + IT_0083;
    const ccomplex_t IT_0085 = IT_0078 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0075*IT_0087;
    const ccomplex_t IT_0089 = IT_0073*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0071 + IT_0090;
    const ccomplex_t IT_0092 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0093 = IT_0053*IT_0092;
    const ccomplex_t IT_0094 = (-2)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0096 = IT_0057*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0099 = IT_0053*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0075*IT_0104;
    const ccomplex_t IT_0106 = IT_0073*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0076*IT_0092;
    const ccomplex_t IT_0109 = (-2)*IT_0108;
    const ccomplex_t IT_0110 = IT_0079*IT_0095;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = IT_0076*IT_0098;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = IT_0109 + IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0052*IT_0117;
    const ccomplex_t IT_0119 = IT_0069*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = -IT_0107 + -IT_0120;
    const ccomplex_t IT_0122 = conj(N_B4)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0123 = IT_0020*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = IT_0002*IT_0049;
    const ccomplex_t IT_0126 = conj(N_d4)*e_em*m_tau*IT_0047*conj(U_stau_01);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0124 + 0.5*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = conj(U_d1)*e_em*m_tau*IT_0047*U_stau_11;
    const ccomplex_t IT_0134 = IT_0125*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-0.5)*IT_0135);
    const ccomplex_t IT_0137 = IT_0130*IT_0136;
    const ccomplex_t IT_0138 = cpow(s_13 + (-0.5)*IT_0013 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0142 = IT_0020*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = conj(N_d4)*e_em*m_tau*IT_0047*conj(U_stau_00);
    const ccomplex_t IT_0145 = IT_0125*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + 0.5*IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0150 = IT_0002*IT_0149;
    const ccomplex_t IT_0151 = conj(U_d1)*e_em*m_tau*IT_0047*U_stau_10;
    const ccomplex_t IT_0152 = IT_0125*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + (-0.5)*IT_0153);
    const ccomplex_t IT_0155 = IT_0148*IT_0154;
    const ccomplex_t IT_0156 = cpow(s_13 + (-0.5)*IT_0013 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (-0.5)*IT_0140 + (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = 4*conj(IT_0017);
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0002*IT_0047*IT_0049;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0028*IT_0162;
    const ccomplex_t IT_0164 = IT_0018*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = 8*conj(IT_0017);
    const ccomplex_t IT_0168 = (-4)*IT_0036*(conj(IT_0046)*IT_0091 + IT_0046
      *conj(IT_0091) + conj(IT_0017)*IT_0121 + IT_0017*(conj(IT_0121) + -conj
      (IT_0159)) + (-0.25)*IT_0159*IT_0160 + (-0.25)*IT_0035*conj(IT_0166) + (
      -0.25)*IT_0166*IT_0167);
    const ccomplex_t IT_0169 = s_13*s_24;
    const ccomplex_t IT_0170 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0171 = IT_0020*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0174 = IT_0002*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = N_d4*e_em*m_tau*IT_0047*conj(U_stau_10);
    const ccomplex_t IT_0177 = IT_0125*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0172 + IT_0175 + 
      -IT_0178);
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = IT_0154*IT_0180;
    const ccomplex_t IT_0182 = IT_0156*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = N_B4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0185 = IT_0020*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = N_W4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0188 = IT_0002*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = N_d4*e_em*m_tau*IT_0047*conj(U_stau_11);
    const ccomplex_t IT_0191 = IT_0125*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0186 + IT_0189 + 
      -IT_0192);
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0136*IT_0194;
    const ccomplex_t IT_0196 = IT_0138*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0183 + 0.5*IT_0197;
    const ccomplex_t IT_0199 = 8*IT_0046;
    const ccomplex_t IT_0200 = s_34*m_C_1*m_N_4;
    const ccomplex_t IT_0201 = 2*IT_0200;
    const ccomplex_t IT_0202 = conj(IT_0121) + conj(IT_0166);
    const ccomplex_t IT_0203 = s_13*m_N_4*m_tau;
    const ccomplex_t IT_0204 = 2*IT_0203;
    const ccomplex_t IT_0205 = 4*IT_0203;
    const ccomplex_t IT_0206 = 4*IT_0169;
    const ccomplex_t IT_0207 = 8*conj(IT_0046);
    const ccomplex_t IT_0208 = (-2)*IT_0036;
    const ccomplex_t IT_0209 = 2*IT_0036;
    const ccomplex_t IT_0210 = 4*IT_0036;
    const ccomplex_t IT_0211 = 4*IT_0000;
    const ccomplex_t IT_0212 = s_12*s_34;
    const ccomplex_t IT_0213 = (-2)*IT_0000 + 2*IT_0169 + (-2)*IT_0212;
    const ccomplex_t IT_0214 = -IT_0000;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = -IT_0169;
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = (-2)*IT_0217;
    const ccomplex_t IT_0219 = IT_0169 + IT_0215;
    const ccomplex_t IT_0220 = 2*IT_0219;
    const ccomplex_t IT_0221 = 4*IT_0046;
    const ccomplex_t IT_0222 = 0.25*conj(IT_0034);
    const ccomplex_t IT_0223 = 4*IT_0200*((IT_0034 + 0.25*IT_0035)*conj
      (IT_0046) + (-0.5)*conj(IT_0091)*IT_0159 + (-0.5)*IT_0091*conj(IT_0159) +
       0.25*IT_0046*IT_0167 + 0.25*IT_0160*IT_0198 + IT_0017*conj(IT_0198) +
       IT_0221*IT_0222);
    const ccomplex_t IT_0224 = 0.25*conj(IT_0159);
    const ccomplex_t IT_0225 = (-4)*IT_0200;
    const ccomplex_t IT_0226 = 0.25*conj(IT_0198);
    const ccomplex_t IT_0227 = 0.25*conj(IT_0166);
    const ccomplex_t IT_0228 = IT_0000*(16*IT_0017*conj(IT_0017) + conj
      (IT_0034)*IT_0035) + IT_0168 + IT_0169*(16*IT_0046*conj(IT_0046) + conj
      (IT_0198)*IT_0199) + IT_0198*(conj(IT_0034)*IT_0201 + IT_0202*IT_0204 +
       conj(IT_0159)*IT_0205 + conj(IT_0198)*IT_0206 + IT_0169*IT_0207 + conj
      (IT_0091)*IT_0208) + IT_0034*(IT_0000*IT_0167 + conj(IT_0198)*IT_0201 +
       conj(IT_0091)*IT_0204 + conj(IT_0121)*IT_0208 + conj(IT_0159)*IT_0209 +
       conj(IT_0166)*IT_0210 + conj(IT_0034)*IT_0211) + IT_0166*(conj(IT_0091)
      *IT_0201 + conj(IT_0198)*IT_0204 + conj(IT_0034)*IT_0210 + conj(IT_0166)
      *IT_0211 + conj(IT_0121)*IT_0213 + conj(IT_0159)*IT_0218) + IT_0159*(conj
      (IT_0198)*IT_0205 + conj(IT_0159)*IT_0206 + conj(IT_0034)*IT_0209 + conj
      (IT_0166)*IT_0218 + conj(IT_0121)*IT_0220) + IT_0223 + 4*IT_0203*(IT_0017
      *conj(IT_0091) + 0.25*IT_0091*IT_0160 + conj(IT_0046)*(IT_0121 + IT_0166) 
      + 0.25*IT_0159*IT_0207 + 0.25*IT_0202*IT_0221 + IT_0199*IT_0224) + 4
      *IT_0091*(conj(IT_0091)*IT_0212 + IT_0204*IT_0222 + 0.25*conj(IT_0121)
      *IT_0225 + IT_0208*IT_0226 + IT_0201*IT_0227) + 4*IT_0121*(conj(IT_0121)
      *IT_0212 + IT_0208*IT_0222 + IT_0220*IT_0224 + 0.25*conj(IT_0091)*IT_0225 
      + IT_0204*IT_0226 + IT_0213*IT_0227);
    return create_ccomplex_return(IT_0228);
}

