#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_2_to_anti_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_2_to_anti_e_nu_e(
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
    const creal_t m_N_2 = param->m_N_2;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = pow(m_W, 2);
    const ccomplex_t IT_0002 = pow(m_C_1, 2);
    const ccomplex_t IT_0003 = pow(m_N_2, 2);
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
    const ccomplex_t IT_0012 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 2*IT_0013;
    const ccomplex_t IT_0015 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0010*IT_0019;
    const ccomplex_t IT_0021 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0017 + IT_0023;
    const ccomplex_t IT_0025 = IT_0014 + IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0027;
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + -IT_0002 + -IT_0003 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0005*IT_0007*IT_0010*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0032;
    const ccomplex_t IT_0037 = IT_0021*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0032;
    const ccomplex_t IT_0040 = IT_0015*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = IT_0038 + IT_0041;
    const ccomplex_t IT_0043 = IT_0012*IT_0039;
    const ccomplex_t IT_0044 = 2*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0035*IT_0047;
    const ccomplex_t IT_0049 = IT_0031*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0030 + IT_0050;
    const ccomplex_t IT_0052 = conj(IT_0030) + conj(IT_0050);
    const ccomplex_t IT_0053 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0054 = IT_0011*IT_0053;
    const ccomplex_t IT_0055 = (-2)*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0060 = IT_0011*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = IT_0058 + IT_0061;
    const ccomplex_t IT_0063 = IT_0055 + IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0035*IT_0065;
    const ccomplex_t IT_0067 = IT_0031*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0039*IT_0053;
    const ccomplex_t IT_0071 = (-2)*IT_0070;
    const ccomplex_t IT_0072 = IT_0036*IT_0056;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = IT_0039*IT_0059;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = IT_0073 + IT_0075;
    const ccomplex_t IT_0077 = IT_0071 + IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0009*IT_0079;
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = s_34*m_C_1*m_N_2;
    const ccomplex_t IT_0085 = (-4)*IT_0030*IT_0084;
    const ccomplex_t IT_0086 = m_e*s_24*m_C_1;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0090 = IT_0007*IT_0089;
    const ccomplex_t IT_0091 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0092 = IT_0007*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + (-0.5)*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0088*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0002 + -IT_0003 + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (-4)*IT_0086*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0005*IT_0007*IT_0033;
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0019;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0102*IT_0104;
    const ccomplex_t IT_0106 = pow(m_e, 2);
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0003 + IT_0106 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = s_13*s_24;
    const ccomplex_t IT_0111 = s_14*s_23;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0110 + IT_0112;
    const ccomplex_t IT_0114 = -IT_0000;
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = (-2)*IT_0109*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B2)*e_em;
    const ccomplex_t IT_0118 = IT_0019*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W2)*e_em;
    const ccomplex_t IT_0121 = IT_0007*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + -IT_0122);
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = IT_0102*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0002 + IT_0106 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0000 + IT_0113;
    const ccomplex_t IT_0130 = (-2)*IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0007;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0005*IT_0007*IT_0033;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = IT_0131*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0003 + IT_0106 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = (-2)*IT_0115*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0007;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0124*IT_0140;
    const ccomplex_t IT_0142 = IT_0126*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = (-2)*IT_0086;
    const ccomplex_t IT_0146 = IT_0144*IT_0145;
    const ccomplex_t IT_0147 = IT_0069 + IT_0083;
    const ccomplex_t IT_0148 = 4*IT_0000;
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0005*IT_0007*IT_0033;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0102*IT_0151;
    const ccomplex_t IT_0153 = IT_0107*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = N_B2*e_em;
    const ccomplex_t IT_0156 = IT_0019*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = N_W2*e_em;
    const ccomplex_t IT_0159 = IT_0007*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + IT_0160);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0131*IT_0162;
    const ccomplex_t IT_0164 = IT_0135*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0154 + -IT_0165;
    const ccomplex_t IT_0167 = m_e*s_14*m_N_2;
    const ccomplex_t IT_0168 = 2*IT_0167;
    const ccomplex_t IT_0169 = IT_0166*IT_0168;
    const ccomplex_t IT_0170 = IT_0085 + IT_0100 + IT_0116 + IT_0130 + IT_0138
       + IT_0146 + IT_0149 + IT_0169;
    const ccomplex_t IT_0171 = 2*IT_0084;
    const ccomplex_t IT_0172 = IT_0110 + IT_0111 + IT_0114;
    const ccomplex_t IT_0173 = conj(IT_0109) + conj(IT_0137);
    const ccomplex_t IT_0174 = 2*IT_0128;
    const ccomplex_t IT_0175 = IT_0109 + IT_0137;
    const ccomplex_t IT_0176 = 2*conj(IT_0128);
    const ccomplex_t IT_0177 = 4*IT_0128;
    const ccomplex_t IT_0178 = 8*IT_0099;
    const ccomplex_t IT_0179 = 8*conj(IT_0099);
    const ccomplex_t IT_0180 = (-4)*conj(IT_0030);
    const ccomplex_t IT_0181 = (-4)*conj(IT_0050);
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = (-2)*IT_0030;
    const ccomplex_t IT_0184 = (-2)*IT_0050;
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = (-2)*conj(IT_0030);
    const ccomplex_t IT_0187 = (-2)*conj(IT_0050);
    const ccomplex_t IT_0188 = IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0190 = IT_0007*IT_0189;
    const ccomplex_t IT_0191 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0192 = IT_0007*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0190 + 0.5*IT_0193);
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = IT_0088*IT_0195;
    const ccomplex_t IT_0197 = IT_0097*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = 4*IT_0198;
    const ccomplex_t IT_0200 = 4*conj(IT_0198);
    const ccomplex_t IT_0201 = conj(IT_0069) + conj(IT_0083);
    const ccomplex_t IT_0202 = (-4)*IT_0050;
    const ccomplex_t IT_0203 = conj(IT_0099)*IT_0166 + IT_0099*conj(IT_0166) +
       0.25*IT_0052*IT_0174 + 0.25*IT_0051*IT_0176 + 0.25*IT_0147*IT_0182 + 0.25
      *IT_0173*IT_0185 + 0.25*IT_0175*IT_0188 + 0.25*IT_0179*IT_0198 + 0.25
      *IT_0178*conj(IT_0198) + 0.25*conj(IT_0144)*IT_0199 + 0.25*IT_0144*IT_0200
       + 0.25*IT_0201*IT_0202;
    const ccomplex_t IT_0204 = conj(IT_0166)*IT_0168;
    const ccomplex_t IT_0205 = (-4)*IT_0086*conj(IT_0099);
    const ccomplex_t IT_0206 = (-2)*conj(IT_0109)*IT_0115;
    const ccomplex_t IT_0207 = conj(IT_0144)*IT_0145;
    const ccomplex_t IT_0208 = (-2)*conj(IT_0128)*IT_0129;
    const ccomplex_t IT_0209 = (-2)*IT_0115*conj(IT_0137);
    const ccomplex_t IT_0210 = 4*IT_0109;
    const ccomplex_t IT_0211 = 4*IT_0137;
    const ccomplex_t IT_0212 = 4*conj(IT_0109);
    const ccomplex_t IT_0213 = 0.5*IT_0166;
    const ccomplex_t IT_0214 = 0.5*IT_0052*(IT_0099 + 0.5*IT_0144) + 0.5
      *IT_0051*(conj(IT_0099) + 0.5*conj(IT_0144)) + (conj(IT_0109) + 0.5*conj
      (IT_0128) + conj(IT_0137))*IT_0198 + (IT_0175 + 0.125*IT_0177)*conj
      (IT_0198) + 0.125*IT_0147*IT_0200 + 0.125*IT_0199*IT_0201 + 0.125*conj
      (IT_0166)*(IT_0174 + IT_0210 + IT_0211) + (conj(IT_0137) + 0.25*IT_0176 +
       0.25*IT_0212)*IT_0213;
    const ccomplex_t IT_0215 = 0.5*conj(IT_0166);
    const ccomplex_t IT_0216 = (conj(IT_0109) + 2*conj(IT_0128) + conj(IT_0137
      ))*IT_0144 + conj(IT_0144)*(IT_0175 + 0.5*IT_0177) + 0.5*IT_0182*IT_0198 +
       (-2)*conj(IT_0198)*(IT_0030 + (-0.25)*IT_0202) + 4*conj(IT_0099)*(IT_0128
       + 0.125*IT_0210 + 0.125*IT_0211) + 4*IT_0099*(conj(IT_0128) + 0.5*conj
      (IT_0137) + 0.125*IT_0212) + IT_0188*IT_0213 + IT_0185*IT_0215;
    const ccomplex_t IT_0217 = 4*IT_0000*IT_0051*IT_0052 + (conj(IT_0069) +
       conj(IT_0083))*IT_0170 + (conj(IT_0144)*IT_0166 + IT_0144*conj(IT_0166))
      *IT_0171 + IT_0172*(IT_0173*IT_0174 + IT_0175*IT_0176) + IT_0110*(16
      *IT_0099*conj(IT_0099) + conj(IT_0128)*IT_0177 + conj(IT_0144)*IT_0178 + 4
      *IT_0144*(conj(IT_0144) + 0.25*IT_0179)) + 4*IT_0084*IT_0203 + IT_0147*
      (IT_0204 + IT_0205 + IT_0206 + IT_0207 + IT_0208 + IT_0209) + 8*IT_0167
      *IT_0214 + 8*IT_0111*(IT_0198*(conj(IT_0166) + 2*conj(IT_0198)) + 1./8
      *IT_0173*IT_0210 + 1./8*conj(IT_0137)*IT_0211 + 1./8*IT_0137*IT_0212 +
       IT_0166*(conj(IT_0198) + IT_0215)) + 2*IT_0086*IT_0216;
    return create_ccomplex_return(IT_0217);
}

