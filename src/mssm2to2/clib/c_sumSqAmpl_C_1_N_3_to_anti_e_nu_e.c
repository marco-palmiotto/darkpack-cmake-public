#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_3_to_anti_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_3_to_anti_e_nu_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = pow(m_W, 2);
    const ccomplex_t IT_0002 = pow(m_C_1, 2);
    const ccomplex_t IT_0003 = pow(m_N_3, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_12 + IT_0001 + -IT_0002 + -IT_0003 
      + -reg_prop, -1);
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 2*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0022 = IT_0011*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0020 + IT_0023;
    const ccomplex_t IT_0025 = IT_0014 + IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0027;
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0005*IT_0007*IT_0010*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0007*IT_0031;
    const ccomplex_t IT_0036 = IT_0012*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0036;
    const ccomplex_t IT_0038 = IT_0016*IT_0031;
    const ccomplex_t IT_0039 = IT_0018*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = IT_0021*IT_0035;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = IT_0040 + IT_0042;
    const ccomplex_t IT_0044 = IT_0037 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0034*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + -IT_0002 + -IT_0003 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = conj(IT_0030) + conj(IT_0050);
    const ccomplex_t IT_0052 = 4*IT_0030;
    const ccomplex_t IT_0053 = 0.25*IT_0052;
    const ccomplex_t IT_0054 = IT_0050 + IT_0053;
    const ccomplex_t IT_0055 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0056 = IT_0011*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0059 = IT_0017*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = IT_0057 + IT_0060;
    const ccomplex_t IT_0062 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0063 = IT_0011*IT_0062;
    const ccomplex_t IT_0064 = (-2)*IT_0063;
    const ccomplex_t IT_0065 = IT_0061 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0034*IT_0067;
    const ccomplex_t IT_0069 = IT_0048*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0035*IT_0062;
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = IT_0038*IT_0058;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = IT_0035*IT_0055;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = IT_0075 + IT_0077;
    const ccomplex_t IT_0079 = IT_0073 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0009*IT_0081;
    const ccomplex_t IT_0083 = IT_0004*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0005*IT_0007*IT_0032;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0005*IT_0007*IT_0032;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = pow(m_e, 2);
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0003 + IT_0091 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0007;
    const ccomplex_t IT_0096 = N_B3*e_em;
    const ccomplex_t IT_0097 = IT_0016*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = N_W3*e_em;
    const ccomplex_t IT_0100 = IT_0007*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0095*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0003 + IT_0091 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0094 + -IT_0107;
    const ccomplex_t IT_0109 = m_e*s_14*m_N_3;
    const ccomplex_t IT_0110 = 2*IT_0109;
    const ccomplex_t IT_0111 = conj(IT_0108)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0115 = IT_0007*IT_0114;
    const ccomplex_t IT_0116 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0117 = IT_0007*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.5)*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0113*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0002 + -IT_0003 + -reg_prop, -1);
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = m_e*s_24*m_C_1;
    const ccomplex_t IT_0126 = (-4)*conj(IT_0124)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0016;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = IT_0087*IT_0128;
    const ccomplex_t IT_0130 = IT_0092*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = s_13*s_24;
    const ccomplex_t IT_0133 = s_14*s_23;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = -IT_0000;
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = (-2)*conj(IT_0131)*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0005*IT_0007*IT_0032;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0095*IT_0140;
    const ccomplex_t IT_0142 = IT_0105*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = (-2)*IT_0137*conj(IT_0143);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0007;
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = conj(N_B3)*e_em;
    const ccomplex_t IT_0148 = IT_0016*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W3)*e_em;
    const ccomplex_t IT_0151 = IT_0007*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + -IT_0152);
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0146*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_13 + IT_0002 + IT_0091 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = (-2)*IT_0125;
    const ccomplex_t IT_0161 = conj(IT_0159)*IT_0160;
    const ccomplex_t IT_0162 = IT_0087*IT_0154;
    const ccomplex_t IT_0163 = IT_0156*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0000 + IT_0135;
    const ccomplex_t IT_0166 = (-2)*conj(IT_0164)*IT_0165;
    const ccomplex_t IT_0167 = conj(IT_0071) + conj(IT_0085);
    const ccomplex_t IT_0168 = 4*IT_0000;
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = IT_0111 + IT_0126 + IT_0138 + IT_0144 + IT_0161
       + IT_0166 + IT_0169;
    const ccomplex_t IT_0171 = s_34*m_C_1*m_N_3;
    const ccomplex_t IT_0172 = 2*IT_0171;
    const ccomplex_t IT_0173 = IT_0132 + IT_0133 + IT_0136;
    const ccomplex_t IT_0174 = conj(IT_0131) + conj(IT_0143);
    const ccomplex_t IT_0175 = 2*IT_0164;
    const ccomplex_t IT_0176 = IT_0131 + IT_0143;
    const ccomplex_t IT_0177 = 2*conj(IT_0164);
    const ccomplex_t IT_0178 = (-4)*IT_0030;
    const ccomplex_t IT_0179 = (-4)*IT_0050;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = IT_0071 + IT_0085;
    const ccomplex_t IT_0182 = (-4)*conj(IT_0030);
    const ccomplex_t IT_0183 = (-4)*conj(IT_0050);
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (-2)*IT_0030;
    const ccomplex_t IT_0186 = (-2)*IT_0050;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = (-2)*conj(IT_0030);
    const ccomplex_t IT_0189 = (-2)*conj(IT_0050);
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = IT_0030 + IT_0050;
    const ccomplex_t IT_0192 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0193 = IT_0007*IT_0192;
    const ccomplex_t IT_0194 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0195 = IT_0007*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 0.5*IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0113*IT_0198;
    const ccomplex_t IT_0200 = IT_0122*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 4*IT_0201;
    const ccomplex_t IT_0203 = 4*conj(IT_0201);
    const ccomplex_t IT_0204 = 8*IT_0124;
    const ccomplex_t IT_0205 = 8*conj(IT_0124);
    const ccomplex_t IT_0206 = conj(IT_0108)*IT_0124 + IT_0108*conj(IT_0124) +
       0.25*IT_0051*IT_0175 + 0.25*IT_0167*IT_0180 + 0.25*IT_0181*IT_0184 + 0.25
      *IT_0174*IT_0187 + 0.25*IT_0176*IT_0190 + 0.25*IT_0177*IT_0191 + 0.25*conj
      (IT_0159)*IT_0202 + 0.25*IT_0159*IT_0203 + 0.25*conj(IT_0201)*IT_0204 +
       0.25*IT_0201*IT_0205;
    const ccomplex_t IT_0207 = 4*IT_0164;
    const ccomplex_t IT_0208 = (-2)*IT_0137*IT_0143;
    const ccomplex_t IT_0209 = (-2)*IT_0164*IT_0165;
    const ccomplex_t IT_0210 = 4*IT_0131;
    const ccomplex_t IT_0211 = 4*IT_0143;
    const ccomplex_t IT_0212 = 4*conj(IT_0131);
    const ccomplex_t IT_0213 = 0.25*IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0108;
    const ccomplex_t IT_0215 = IT_0054*conj(IT_0124) + IT_0051*(IT_0124 + 0.5
      *IT_0159) + 0.5*conj(IT_0159)*IT_0191 + (conj(IT_0164) + 2*IT_0174)
      *IT_0201 + 0.25*IT_0167*IT_0202 + 0.25*IT_0181*IT_0203 + 2*conj(IT_0201)*
      (IT_0176 + 0.125*IT_0207) + 0.25*conj(IT_0108)*(IT_0175 + IT_0210 +
       IT_0211) + 2*(conj(IT_0143) + 0.25*IT_0177 + IT_0213)*IT_0214;
    const ccomplex_t IT_0216 = 0.5*conj(IT_0108);
    const ccomplex_t IT_0217 = IT_0159*(conj(IT_0164) + 0.5*IT_0174) + 0.25
      *IT_0184*IT_0201 + 0.25*IT_0180*conj(IT_0201) + 0.5*conj(IT_0159)*(IT_0176
       + 0.5*IT_0207) + 2*conj(IT_0124)*(IT_0164 + 0.125*IT_0210 + 0.125*IT_0211
      ) + -IT_0124*(conj(IT_0071) + conj(IT_0085) + -conj(IT_0143) + (-2)*conj
      (IT_0164) + -IT_0213) + 0.5*IT_0190*IT_0214 + 0.5*IT_0187*IT_0216;
    const ccomplex_t IT_0218 = 4*IT_0000*IT_0051*IT_0054 + (IT_0071 + IT_0085)
      *IT_0170 + (conj(IT_0108)*IT_0159 + IT_0108*conj(IT_0159))*IT_0172 +
       IT_0173*(IT_0174*IT_0175 + IT_0176*IT_0177) + 4*IT_0171*IT_0206 + IT_0132
      *(16*IT_0124*conj(IT_0124) + conj(IT_0159)*IT_0204 + 4*IT_0159*(conj
      (IT_0159) + 0.25*IT_0205) + conj(IT_0164)*IT_0207) + (-2)*IT_0167*((-0.5)
      *IT_0108*IT_0110 + IT_0131*IT_0137 + (-0.5)*IT_0159*IT_0160 + (-0.5)
      *IT_0208 + (-0.5)*IT_0209) + 4*IT_0109*IT_0215 + 8*IT_0133*(IT_0201*(conj
      (IT_0108) + 2*conj(IT_0201)) + 1./8*IT_0174*IT_0210 + 1./8*conj(IT_0143)
      *IT_0211 + 1./8*IT_0143*IT_0212 + IT_0108*(conj(IT_0201) + IT_0216)) + 4
      *IT_0125*IT_0217;
    return create_ccomplex_return(IT_0218);
}

