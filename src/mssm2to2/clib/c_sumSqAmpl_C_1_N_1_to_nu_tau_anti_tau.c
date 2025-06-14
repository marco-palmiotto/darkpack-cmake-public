#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_1_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_1_to_nu_tau_anti_tau(
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
    const creal_t m_N_1 = param->m_N_1;
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
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + (-0.5)*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = pow(m_C_1, 2);
    const ccomplex_t IT_0013 = pow(m_N_1, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = s_14*s_23;
    const ccomplex_t IT_0018 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + 0.5*IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = s_34*m_C_1*m_N_1;
    const ccomplex_t IT_0029 = 8*IT_0016;
    const ccomplex_t IT_0030 = s_13*s_24;
    const ccomplex_t IT_0031 = 8*conj(IT_0016);
    const ccomplex_t IT_0032 = s_23*m_C_1*m_tau;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0001*IT_0033*IT_0035*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0001*IT_0034;
    const ccomplex_t IT_0040 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0041;
    const ccomplex_t IT_0043 = cos(theta_W);
    const ccomplex_t IT_0044 = cpow(IT_0043, -1);
    const ccomplex_t IT_0045 = IT_0034*IT_0044;
    const ccomplex_t IT_0046 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0050 = IT_0039*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = IT_0048 + IT_0051;
    const ccomplex_t IT_0053 = IT_0042 + IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0038*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = pow(m_W, 2);
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0060 
      + -reg_prop, -1);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0001*IT_0033;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0001*IT_0036;
    const ccomplex_t IT_0065 = IT_0040*IT_0064;
    const ccomplex_t IT_0066 = 2*IT_0065;
    const ccomplex_t IT_0067 = IT_0036*IT_0044;
    const ccomplex_t IT_0068 = IT_0046*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = IT_0049*IT_0064;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = IT_0066 + IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0063*IT_0075;
    const ccomplex_t IT_0077 = IT_0061*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0059 + IT_0078;
    const ccomplex_t IT_0080 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0081 = IT_0067*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0084 = IT_0064*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = IT_0082 + IT_0085;
    const ccomplex_t IT_0087 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0088 = IT_0064*IT_0087;
    const ccomplex_t IT_0089 = (-2)*IT_0088;
    const ccomplex_t IT_0090 = IT_0086 + IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0038*IT_0092;
    const ccomplex_t IT_0094 = IT_0057*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0039*IT_0087;
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = IT_0045*IT_0080;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = IT_0039*IT_0083;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = IT_0099 + IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = IT_0063*IT_0105;
    const ccomplex_t IT_0107 = IT_0061*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0095 + -IT_0108;
    const ccomplex_t IT_0110 = cpow(s_13 + (-0.5)*IT_0012 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0111 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0112 = IT_0044*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = IT_0001*IT_0035;
    const ccomplex_t IT_0115 = conj(N_d1)*e_em*m_tau*IT_0033*conj(U_stau_00);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0113 + 0.5*IT_0117);
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = conj(U_d1)*e_em*m_tau*IT_0033*U_stau_10;
    const ccomplex_t IT_0123 = IT_0114*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + (-0.5)*IT_0124);
    const ccomplex_t IT_0126 = IT_0119*IT_0125;
    const ccomplex_t IT_0127 = IT_0110*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = cpow(s_13 + (-0.5)*IT_0012 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0130 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0131 = IT_0044*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_d1)*e_em*m_tau*IT_0033*conj(U_stau_01);
    const ccomplex_t IT_0134 = IT_0114*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0139 = IT_0001*IT_0138;
    const ccomplex_t IT_0140 = conj(U_d1)*e_em*m_tau*IT_0033*U_stau_11;
    const ccomplex_t IT_0141 = IT_0114*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.5)*IT_0142);
    const ccomplex_t IT_0144 = IT_0137*IT_0143;
    const ccomplex_t IT_0145 = IT_0129*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = (-0.5)*IT_0128 + (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = 4*conj(IT_0016);
    const ccomplex_t IT_0149 = cpow(s_23 + (-0.5)*IT_0013 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0001*IT_0033*IT_0035;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B1)*e_em;
    const ccomplex_t IT_0153 = IT_0044*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W1)*e_em;
    const ccomplex_t IT_0156 = IT_0001*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + -IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0151*IT_0159;
    const ccomplex_t IT_0161 = IT_0149*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = (-4)*IT_0032*(conj(IT_0027)*IT_0079 + IT_0027
      *conj(IT_0079) + conj(IT_0016)*IT_0109 + IT_0016*(conj(IT_0109) + -conj
      (IT_0147)) + (-0.25)*IT_0147*IT_0148 + (-0.25)*IT_0031*IT_0163 + (-0.25)
      *IT_0029*conj(IT_0163));
    const ccomplex_t IT_0165 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0166 = IT_0044*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0169 = IT_0001*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = N_d1*e_em*m_tau*IT_0033*conj(U_stau_10);
    const ccomplex_t IT_0172 = IT_0114*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + IT_0170 + 
      -IT_0173);
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = IT_0125*IT_0175;
    const ccomplex_t IT_0177 = IT_0110*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0180 = IT_0044*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0183 = IT_0001*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = N_d1*e_em*m_tau*IT_0033*conj(U_stau_11);
    const ccomplex_t IT_0186 = IT_0114*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0181 + IT_0184 + 
      -IT_0187);
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0143*IT_0189;
    const ccomplex_t IT_0191 = IT_0129*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 0.5*IT_0178 + 0.5*IT_0192;
    const ccomplex_t IT_0194 = conj(IT_0109) + conj(IT_0163);
    const ccomplex_t IT_0195 = s_13*m_N_1*m_tau;
    const ccomplex_t IT_0196 = 2*IT_0195;
    const ccomplex_t IT_0197 = 4*IT_0195;
    const ccomplex_t IT_0198 = 8*conj(IT_0027);
    const ccomplex_t IT_0199 = (-2)*IT_0032;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0159*IT_0201;
    const ccomplex_t IT_0203 = IT_0149*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = 2*IT_0028;
    const ccomplex_t IT_0207 = 2*IT_0032;
    const ccomplex_t IT_0208 = 4*IT_0027;
    const ccomplex_t IT_0209 = 4*IT_0032;
    const ccomplex_t IT_0210 = 4*IT_0017;
    const ccomplex_t IT_0211 = 4*conj(IT_0027);
    const ccomplex_t IT_0212 = s_12*s_34;
    const ccomplex_t IT_0213 = (-2)*IT_0017 + 2*IT_0030 + (-2)*IT_0212;
    const ccomplex_t IT_0214 = -IT_0017;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = -IT_0030;
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = (-2)*IT_0217;
    const ccomplex_t IT_0219 = IT_0030 + IT_0215;
    const ccomplex_t IT_0220 = 2*IT_0219;
    const ccomplex_t IT_0221 = 4*IT_0030;
    const ccomplex_t IT_0222 = 8*IT_0027;
    const ccomplex_t IT_0223 = 0.25*IT_0079;
    const ccomplex_t IT_0224 = 0.25*IT_0147;
    const ccomplex_t IT_0225 = IT_0109 + IT_0163;
    const ccomplex_t IT_0226 = 0.25*IT_0225;
    const ccomplex_t IT_0227 = 0.25*conj(IT_0147);
    const ccomplex_t IT_0228 = (-4)*IT_0028;
    const ccomplex_t IT_0229 = 0.25*conj(IT_0163);
    const ccomplex_t IT_0230 = 16*IT_0016*conj(IT_0016)*IT_0017 + conj(IT_0027
      )*IT_0028*IT_0029 + IT_0027*(16*conj(IT_0027)*IT_0030 + IT_0028*IT_0031) +
       IT_0164 + IT_0193*(IT_0028*IT_0148 + IT_0194*IT_0196 + conj(IT_0147)
      *IT_0197 + IT_0030*IT_0198 + conj(IT_0079)*IT_0199) + conj(IT_0205)*
      (IT_0017*IT_0029 + IT_0079*IT_0196 + IT_0109*IT_0199 + IT_0193*IT_0206 +
       IT_0147*IT_0207 + IT_0028*IT_0208 + IT_0163*IT_0209 + IT_0205*IT_0210) +
       IT_0205*(IT_0017*IT_0031 + conj(IT_0079)*IT_0196 + conj(IT_0109)*IT_0199 
      + conj(IT_0147)*IT_0207 + conj(IT_0163)*IT_0209 + IT_0028*IT_0211) +
       IT_0163*(conj(IT_0079)*IT_0206 + conj(IT_0163)*IT_0210 + conj(IT_0109)
      *IT_0213 + conj(IT_0147)*IT_0218) + (-2)*IT_0147*(IT_0028*conj(IT_0079) + 
      (-0.5)*conj(IT_0163)*IT_0218 + (-0.5)*conj(IT_0109)*IT_0220 + (-0.5)*conj
      (IT_0147)*IT_0221) + 4*conj(IT_0193)*(IT_0016*IT_0028 + 0.25*IT_0205
      *IT_0206 + 0.25*IT_0193*IT_0221 + 0.25*IT_0030*IT_0222 + IT_0199*IT_0223 +
       IT_0197*IT_0224 + IT_0196*IT_0226) + 4*IT_0195*(IT_0016*conj(IT_0079) +
       0.25*IT_0194*IT_0208 + IT_0148*IT_0223 + IT_0198*IT_0224 + IT_0211
      *IT_0226 + IT_0222*IT_0227) + 4*IT_0079*((-0.5)*IT_0028*conj(IT_0147) +
       conj(IT_0079)*IT_0212 + 0.25*conj(IT_0109)*IT_0228 + IT_0206*IT_0229) + 4
      *IT_0109*(conj(IT_0109)*IT_0212 + IT_0220*IT_0227 + 0.25*conj(IT_0079)
      *IT_0228 + IT_0213*IT_0229);
    return create_ccomplex_return(IT_0230);
}

