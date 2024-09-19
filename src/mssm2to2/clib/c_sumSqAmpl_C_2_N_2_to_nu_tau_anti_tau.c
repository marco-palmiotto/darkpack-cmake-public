#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_nu_tau_anti_tau(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = conj(N_d2)*e_em*m_tau*IT_0010*conj(U_stau_00);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0017 = IT_0008*IT_0016;
    const ccomplex_t IT_0018 = conj(U_d2)*e_em*m_tau*IT_0010*U_stau_10;
    const ccomplex_t IT_0019 = IT_0009*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0017 + (-0.5)*IT_0020);
    const ccomplex_t IT_0022 = IT_0015*IT_0021;
    const ccomplex_t IT_0023 = pow(m_C_2, 2);
    const ccomplex_t IT_0024 = cpow(s_13 + (-0.5)*IT_0023 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0028 = IT_0001*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = conj(N_d2)*e_em*m_tau*IT_0010*conj(U_stau_01);
    const ccomplex_t IT_0031 = IT_0009*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0029 + 0.5*IT_0032);
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0036 = IT_0008*IT_0035;
    const ccomplex_t IT_0037 = conj(U_d2)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0038 = IT_0009*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-0.5)*IT_0039);
    const ccomplex_t IT_0041 = IT_0034*IT_0040;
    const ccomplex_t IT_0042 = cpow(s_13 + (-0.5)*IT_0023 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0026 + (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = pow(m_N_2, 2);
    const ccomplex_t IT_0047 = cpow(s_23 + (-0.5)*IT_0046 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0048 = conj(N_B2)*e_em;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W2)*e_em;
    const ccomplex_t IT_0052 = IT_0008*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + -IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = s_12*s_34;
    const ccomplex_t IT_0063 = s_14*s_23;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = s_13*s_24;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = IT_0065 + IT_0066;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0008*IT_0010;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0005*IT_0008;
    const ccomplex_t IT_0073 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (-2)*IT_0074;
    const ccomplex_t IT_0076 = IT_0001*IT_0005;
    const ccomplex_t IT_0077 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0081 = IT_0072*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = IT_0079 + IT_0082;
    const ccomplex_t IT_0084 = IT_0075 + IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0071*IT_0086;
    const ccomplex_t IT_0088 = pow(m_W, 2);
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + -IT_0023 + -IT_0046 + IT_0088 
      + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = sin(beta);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0006*IT_0008*IT_0010*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0008*IT_0092;
    const ccomplex_t IT_0096 = IT_0073*IT_0095;
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = IT_0001*IT_0092;
    const ccomplex_t IT_0099 = IT_0077*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0080*IT_0095;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = IT_0097 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0094*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_12 + -IT_0023 + -IT_0046 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = -IT_0091 + -IT_0110;
    const ccomplex_t IT_0112 = 2*conj(IT_0111);
    const ccomplex_t IT_0113 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0117 = IT_0008*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_d2*e_em*m_tau*IT_0010*conj(U_stau_10);
    const ccomplex_t IT_0120 = IT_0009*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0115 + IT_0118 + 
      -IT_0121);
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0021*IT_0123;
    const ccomplex_t IT_0125 = IT_0024*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0128 = IT_0001*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0131 = IT_0008*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = N_d2*e_em*m_tau*IT_0010*conj(U_stau_11);
    const ccomplex_t IT_0134 = IT_0009*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0129 + IT_0132 + 
      -IT_0135);
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0040*IT_0137;
    const ccomplex_t IT_0139 = IT_0042*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = 0.5*IT_0126 + 0.5*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0145 = IT_0008*IT_0144;
    const ccomplex_t IT_0146 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0147 = IT_0008*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + 0.5*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = IT_0143*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0023 + -IT_0046 + -reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = 4*IT_0141;
    const ccomplex_t IT_0156 = 0.125*IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = 4*IT_0045;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0055*IT_0160;
    const ccomplex_t IT_0162 = IT_0047*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0166 = IT_0008*IT_0165;
    const ccomplex_t IT_0167 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0168 = IT_0008*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + (-0.5)*IT_0169);
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = IT_0143*IT_0171;
    const ccomplex_t IT_0173 = IT_0152*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = 4*IT_0164;
    const ccomplex_t IT_0176 = 4*IT_0061;
    const ccomplex_t IT_0177 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0178 = IT_0072*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0181 = IT_0076*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0184 = IT_0072*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = IT_0182 + IT_0185;
    const ccomplex_t IT_0187 = IT_0179 + IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = IT_0094*IT_0189;
    const ccomplex_t IT_0191 = IT_0108*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0095*IT_0177;
    const ccomplex_t IT_0194 = 2*IT_0193;
    const ccomplex_t IT_0195 = IT_0098*IT_0180;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = IT_0095*IT_0183;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0194 + IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0201;
    const ccomplex_t IT_0203 = IT_0071*IT_0202;
    const ccomplex_t IT_0204 = IT_0089*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0192 + IT_0205;
    const ccomplex_t IT_0207 = 4*IT_0206;
    const ccomplex_t IT_0208 = IT_0062 + IT_0063 + IT_0067;
    const ccomplex_t IT_0209 = (-2)*conj(IT_0111);
    const ccomplex_t IT_0210 = (-2)*IT_0061;
    const ccomplex_t IT_0211 = s_13*m_N_2*m_tau;
    const ccomplex_t IT_0212 = 2*IT_0206;
    const ccomplex_t IT_0213 = 2*conj(IT_0206);
    const ccomplex_t IT_0214 = 4*conj(IT_0141);
    const ccomplex_t IT_0215 = 0.5*IT_0045;
    const ccomplex_t IT_0216 = 0.5*conj(IT_0045);
    const ccomplex_t IT_0217 = (conj(IT_0061) + 0.5*IT_0112)*IT_0141 + 
      (IT_0061 + IT_0111)*conj(IT_0141) + 4*(conj(IT_0045) + 1./2*conj(IT_0061) 
      + 1./2*conj(IT_0111))*IT_0154 + 4*conj(IT_0154)*(IT_0045 + 1./2*IT_0111 +
       0.125*IT_0176) + 2*IT_0174*conj(IT_0206) + 0.5*conj(IT_0174)*IT_0207 +
       0.5*conj(IT_0164)*IT_0212 + 0.5*IT_0164*IT_0213 + IT_0214*IT_0215 +
       IT_0155*IT_0216;
    const ccomplex_t IT_0218 = s_34*m_C_2*m_N_2;
    const ccomplex_t IT_0219 = 2*IT_0164;
    const ccomplex_t IT_0220 = 4*conj(IT_0164);
    const ccomplex_t IT_0221 = (-4)*IT_0206;
    const ccomplex_t IT_0222 = (-4)*conj(IT_0206);
    const ccomplex_t IT_0223 = (-2)*conj(IT_0206);
    const ccomplex_t IT_0224 = (-2)*IT_0206;
    const ccomplex_t IT_0225 = IT_0141*conj(IT_0164) + 4*IT_0157*conj(IT_0174)
       + 0.5*conj(IT_0154)*IT_0175 + 0.5*conj(IT_0061)*IT_0212 + 0.5*IT_0061
      *IT_0213 + 4*IT_0174*(conj(IT_0154) + 0.125*IT_0214) + 0.5*conj(IT_0141)
      *IT_0219 + 0.5*IT_0154*IT_0220 + 0.5*conj(IT_0111)*IT_0221 + 0.5*IT_0111
      *IT_0222 + IT_0215*IT_0223 + IT_0216*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*conj(IT_0061);
    const ccomplex_t IT_0227 = s_23*m_C_2*m_tau;
    const ccomplex_t IT_0228 = -(IT_0045 + -IT_0111)*conj(IT_0164) + (-2)*
      (conj(IT_0045) + 2*conj(IT_0061) + -conj(IT_0111))*IT_0174 + (-0.5)
      *IT_0164*IT_0209 + 2*conj(IT_0174)*(IT_0111 + (-0.25)*IT_0158 + IT_0210) +
       (-0.5)*conj(IT_0045)*IT_0219 + (-0.5)*IT_0061*IT_0220 + (-0.5)*conj
      (IT_0154)*IT_0221 + (-0.5)*IT_0154*IT_0222 + (-0.5)*IT_0141*IT_0223 + (
      -0.5)*conj(IT_0141)*IT_0224 + IT_0175*IT_0226;
    const ccomplex_t IT_0229 = (-2)*IT_0045*(conj(IT_0061)*IT_0068 + (-0.5)
      *IT_0069*IT_0112) + 8*IT_0066*((IT_0141 + 2*IT_0154)*conj(IT_0154) + conj
      (IT_0141)*IT_0157 + 0.125*conj(IT_0045)*IT_0158) + 8*IT_0063*((IT_0164 + 2
      *IT_0174)*conj(IT_0174) + conj(IT_0164)*(IT_0174 + 0.125*IT_0175) + 0.125
      *conj(IT_0061)*IT_0176) + IT_0062*conj(IT_0206)*IT_0207 + (conj(IT_0045)
      *IT_0068 + (-0.5)*IT_0208*IT_0209)*IT_0210 + 2*IT_0211*IT_0217 + 2*IT_0218
      *IT_0225 + 4*IT_0111*(IT_0062*conj(IT_0111) + IT_0069*IT_0216 + IT_0208
      *IT_0226) + (-2)*IT_0227*IT_0228;
    return create_ccomplex_return(IT_0229);
}

