#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_4_to_anti_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_4_to_anti_e_nu_e(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_4 = param->m_N_4;
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = s_13*s_24;
    const ccomplex_t IT_0002 = s_14*s_23;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = IT_0000 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_e, 2);
    const ccomplex_t IT_0007 = pow(m_C_2, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = conj(N_B4)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = conj(N_W4)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0013 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0015*IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = sin(beta);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0015*IT_0021*IT_0023*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0015*IT_0029;
    const ccomplex_t IT_0033 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (-2)*IT_0034;
    const ccomplex_t IT_0036 = IT_0010*IT_0029;
    const ccomplex_t IT_0037 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0041 = IT_0032*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = IT_0039 + IT_0042;
    const ccomplex_t IT_0044 = IT_0035 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0031*IT_0046;
    const ccomplex_t IT_0048 = pow(m_N_4, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0007 + -IT_0048 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0052;
    const ccomplex_t IT_0054 = pow(m_W, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + -IT_0007 + -IT_0048 + IT_0054 
      + -reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0015*IT_0021;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0015*IT_0022;
    const ccomplex_t IT_0059 = IT_0033*IT_0058;
    const ccomplex_t IT_0060 = (-2)*IT_0059;
    const ccomplex_t IT_0061 = IT_0010*IT_0022;
    const ccomplex_t IT_0062 = IT_0037*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0040*IT_0058;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = IT_0060 + IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0057*IT_0069;
    const ccomplex_t IT_0071 = IT_0055*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = (-2)*IT_0073;
    const ccomplex_t IT_0075 = IT_0053 + IT_0074;
    const ccomplex_t IT_0076 = (-2)*conj(IT_0052);
    const ccomplex_t IT_0077 = (-2)*conj(IT_0073);
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = -IT_0000;
    const ccomplex_t IT_0080 = IT_0004 + IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0006 + IT_0048 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0010;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = IT_0025*IT_0083;
    const ccomplex_t IT_0085 = IT_0081*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0015;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0015*IT_0021*IT_0023;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0006 + IT_0048 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = conj(IT_0086) + conj(IT_0093);
    const ccomplex_t IT_0095 = IT_0086 + IT_0093;
    const ccomplex_t IT_0096 = IT_0001 + IT_0002 + IT_0079;
    const ccomplex_t IT_0097 = 2*IT_0028;
    const ccomplex_t IT_0098 = 2*conj(IT_0028);
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0015;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0020*IT_0100;
    const ccomplex_t IT_0102 = IT_0008*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0007 + -IT_0048 + -reg_prop, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0015;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0109 = IT_0015*IT_0108;
    const ccomplex_t IT_0110 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + (-0.5)*IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0107*IT_0114;
    const ccomplex_t IT_0116 = IT_0105*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 2*conj(IT_0117);
    const ccomplex_t IT_0119 = conj(IT_0104) + IT_0118;
    const ccomplex_t IT_0120 = 4*IT_0028;
    const ccomplex_t IT_0121 = conj(IT_0052) + conj(IT_0073);
    const ccomplex_t IT_0122 = 4*IT_0052;
    const ccomplex_t IT_0123 = 4*conj(IT_0052);
    const ccomplex_t IT_0124 = 4*IT_0073;
    const ccomplex_t IT_0125 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0126 = IT_0032*IT_0125;
    const ccomplex_t IT_0127 = 2*IT_0126;
    const ccomplex_t IT_0128 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0129 = IT_0036*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0132 = IT_0032*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = IT_0130 + IT_0133;
    const ccomplex_t IT_0135 = IT_0127 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0057*IT_0137;
    const ccomplex_t IT_0139 = IT_0055*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0058*IT_0125;
    const ccomplex_t IT_0142 = 2*IT_0141;
    const ccomplex_t IT_0143 = IT_0061*IT_0128;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0058*IT_0131;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = IT_0144 + IT_0146;
    const ccomplex_t IT_0148 = IT_0142 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0031*IT_0150;
    const ccomplex_t IT_0152 = IT_0049*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = conj(IT_0140) + conj(IT_0153);
    const ccomplex_t IT_0155 = 4*IT_0140;
    const ccomplex_t IT_0156 = 4*conj(IT_0140);
    const ccomplex_t IT_0157 = 4*IT_0153;
    const ccomplex_t IT_0158 = m_e*s_14*m_N_4;
    const ccomplex_t IT_0159 = 2*IT_0104;
    const ccomplex_t IT_0160 = IT_0140 + IT_0153;
    const ccomplex_t IT_0161 = 2*conj(IT_0104);
    const ccomplex_t IT_0162 = N_B4*e_em;
    const ccomplex_t IT_0163 = IT_0010*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = N_W4*e_em;
    const ccomplex_t IT_0166 = IT_0015*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0164 + IT_0167);
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0087*IT_0169;
    const ccomplex_t IT_0171 = IT_0091*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0015*IT_0021*IT_0023;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = IT_0025*IT_0174;
    const ccomplex_t IT_0176 = IT_0081*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = -IT_0172 + -IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = IT_0052 + IT_0073;
    const ccomplex_t IT_0181 = 2*conj(IT_0178);
    const ccomplex_t IT_0182 = 4*IT_0178;
    const ccomplex_t IT_0183 = 4*conj(IT_0178);
    const ccomplex_t IT_0184 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0185 = IT_0015*IT_0184;
    const ccomplex_t IT_0186 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0187 = IT_0015*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + 0.5*IT_0188);
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = IT_0107*IT_0190;
    const ccomplex_t IT_0192 = IT_0105*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = 4*IT_0117*(conj(IT_0153) + 0.25*IT_0156) + conj
      (IT_0117)*(IT_0155 + IT_0157) + IT_0154*IT_0159 + IT_0160*IT_0161 + (conj
      (IT_0028) + IT_0121)*IT_0179 + (IT_0028 + IT_0180)*IT_0181 + IT_0094
      *IT_0182 + IT_0095*IT_0183 + 4*(conj(IT_0028) + conj(IT_0073) + 2*conj
      (IT_0086) + 2*conj(IT_0093) + 0.25*IT_0123)*IT_0193 + 8*(IT_0095 + 0.125
      *IT_0120 + 0.125*IT_0122 + 0.125*IT_0124)*conj(IT_0193);
    const ccomplex_t IT_0195 = 4*IT_0086;
    const ccomplex_t IT_0196 = 4*conj(IT_0086);
    const ccomplex_t IT_0197 = 4*IT_0093;
    const ccomplex_t IT_0198 = m_e*s_24*m_C_2;
    const ccomplex_t IT_0199 = (-4)*conj(IT_0052);
    const ccomplex_t IT_0200 = (-4)*conj(IT_0073);
    const ccomplex_t IT_0201 = (-2)*conj(IT_0104);
    const ccomplex_t IT_0202 = (-2)*IT_0178;
    const ccomplex_t IT_0203 = (-2)*conj(IT_0178);
    const ccomplex_t IT_0204 = (-4)*IT_0052;
    const ccomplex_t IT_0205 = (-4)*IT_0073;
    const ccomplex_t IT_0206 = 0.25*IT_0095;
    const ccomplex_t IT_0207 = 0.25*IT_0094;
    const ccomplex_t IT_0208 = IT_0028*conj(IT_0104) + IT_0104*(conj(IT_0028) 
      + (-0.5)*IT_0121) + 2*IT_0117*(conj(IT_0028) + 0.5*conj(IT_0093) + 0.125
      *IT_0196 + 0.125*IT_0199 + 0.125*IT_0200) + 0.25*IT_0180*IT_0201 + 
      -IT_0154*(IT_0193 + (-0.25)*IT_0202) + -IT_0160*(conj(IT_0193) + (-0.25)
      *IT_0203) + IT_0118*(IT_0028 + 0.125*IT_0195 + 0.125*IT_0197 + 0.125
      *IT_0204 + 0.125*IT_0205) + IT_0161*IT_0206 + IT_0159*IT_0207;
    const ccomplex_t IT_0209 = s_34*m_C_2*m_N_4;
    const ccomplex_t IT_0210 = (-2)*conj(IT_0140);
    const ccomplex_t IT_0211 = (-2)*conj(IT_0153);
    const ccomplex_t IT_0212 = (-2)*IT_0140;
    const ccomplex_t IT_0213 = (-2)*IT_0153;
    const ccomplex_t IT_0214 = 0.25*IT_0161*IT_0178 + 0.25*IT_0159*conj
      (IT_0178) + 0.25*conj(IT_0117)*IT_0182 + 0.25*IT_0117*IT_0183 + IT_0119
      *IT_0193 + (IT_0104 + 2*IT_0117)*conj(IT_0193) + 0.25*IT_0160*(IT_0098 +
       IT_0199 + IT_0200) + 0.25*IT_0154*(IT_0097 + IT_0204 + IT_0205) + IT_0206
      *(IT_0210 + IT_0211) + IT_0207*(IT_0212 + IT_0213);
    const ccomplex_t IT_0215 = IT_0005*(conj(IT_0028)*IT_0075 + IT_0028
      *IT_0078) + IT_0080*(IT_0075*IT_0094 + IT_0078*IT_0095) + IT_0096*(IT_0094
      *IT_0097 + IT_0095*IT_0098) + 4*IT_0001*((IT_0104 + 2*IT_0117)*IT_0119 +
       0.25*conj(IT_0028)*IT_0120) + IT_0000*(IT_0121*IT_0122 + IT_0073*IT_0123 
      + conj(IT_0073)*IT_0124 + IT_0154*IT_0155 + IT_0153*IT_0156 + conj(IT_0153
      )*IT_0157) + IT_0158*IT_0194 + 8*IT_0002*(conj(IT_0178)*(0.125*IT_0182 +
       IT_0193) + (IT_0178 + 2*IT_0193)*conj(IT_0193) + 0.125*IT_0094*IT_0195 +
       0.125*IT_0093*IT_0196 + 0.125*conj(IT_0093)*IT_0197) + 4*IT_0198*IT_0208 
      + 4*IT_0209*IT_0214;
    return create_ccomplex_return(IT_0215);
}

