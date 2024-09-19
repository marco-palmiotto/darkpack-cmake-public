#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_R_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_R_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B1)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = pow(m_c, 2);
    const ccomplex_t IT_0019 = pow(m_sc_R, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B2)*e_em;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W2)*e_em;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 3*IT_0028);
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B3)*e_em;
    const ccomplex_t IT_0038 = IT_0001*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W3)*e_em;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 3*IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B4)*e_em;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W4)*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 3*IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0022 + -IT_0036 + -IT_0050 + -IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0001;
    const ccomplex_t IT_0067 = 0.666666666666667*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0015;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = IT_0020*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0001;
    const ccomplex_t IT_0073 = 0.666666666666667*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0031;
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = IT_0034*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0001;
    const ccomplex_t IT_0079 = 0.666666666666667*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0045;
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = IT_0048*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0001;
    const ccomplex_t IT_0085 = 0.666666666666667*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0059;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = IT_0062*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0071 + -IT_0077 + -IT_0083 + -IT_0089;
    const ccomplex_t IT_0091 = -s_14;
    const ccomplex_t IT_0092 = s_13 + IT_0091;
    const ccomplex_t IT_0093 = m_c*IT_0092;
    const ccomplex_t IT_0094 = IT_0016*IT_0068;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0096 = IT_0020*IT_0094*IT_0095;
    const ccomplex_t IT_0097 = IT_0032*IT_0074;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0099 = IT_0034*IT_0097*IT_0098;
    const ccomplex_t IT_0100 = IT_0046*IT_0080;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0102 = IT_0048*IT_0100*IT_0101;
    const ccomplex_t IT_0103 = IT_0060*IT_0086;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0105 = IT_0062*IT_0103*IT_0104;
    const ccomplex_t IT_0106 = m_c*IT_0022;
    const ccomplex_t IT_0107 = m_c*IT_0036;
    const ccomplex_t IT_0108 = m_c*IT_0050;
    const ccomplex_t IT_0109 = m_c*IT_0064;
    const ccomplex_t IT_0110 = -IT_0096 + -IT_0099 + -IT_0102 + -IT_0105 +
       IT_0106 + IT_0107 + IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = IT_0011*IT_0067;
    const ccomplex_t IT_0112 = IT_0020*IT_0095*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0114 = cpow(IT_0113, 2);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_sG*IT_0114;
    const ccomplex_t IT_0116 = pow(m_sG, 2);
    const ccomplex_t IT_0117 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -IT_0116 +
       reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = IT_0030*IT_0073;
    const ccomplex_t IT_0120 = IT_0034*IT_0098*IT_0119;
    const ccomplex_t IT_0121 = IT_0044*IT_0079;
    const ccomplex_t IT_0122 = IT_0048*IT_0101*IT_0121;
    const ccomplex_t IT_0123 = IT_0058*IT_0085;
    const ccomplex_t IT_0124 = IT_0062*IT_0104*IT_0123;
    const ccomplex_t IT_0125 = m_c*IT_0071;
    const ccomplex_t IT_0126 = m_c*IT_0077;
    const ccomplex_t IT_0127 = m_c*IT_0083;
    const ccomplex_t IT_0128 = m_c*IT_0089;
    const ccomplex_t IT_0129 = -IT_0112 + 0.166666666666667*IT_0118 + -IT_0120
       + -IT_0122 + -IT_0124 + IT_0125 + IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = -IT_0018;
    const ccomplex_t IT_0131 = s_34 + IT_0130;
    const ccomplex_t IT_0132 = sin(alpha);
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0006
      *IT_0012*IT_0014*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = cos(alpha);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_c*e_em*IT_0006*IT_0012
      *IT_0014*IT_0135;
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0134*IT_0137;
    const ccomplex_t IT_0139 = pow(m_sc_L, 2);
    const ccomplex_t IT_0140 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0019 + -IT_0139 + -reg_prop, -1);
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0006
      *IT_0012*IT_0014*IT_0135;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*m_c*e_em*IT_0006*IT_0012
      *IT_0014*IT_0132;
    const ccomplex_t IT_0146 = (-0.5)*IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_12 + -IT_0019 + -IT_0139 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = -IT_0142 + -IT_0150;
    const ccomplex_t IT_0152 = 6*IT_0151;
    const ccomplex_t IT_0153 = (-18)*IT_0110;
    const ccomplex_t IT_0154 = m_c*s_13;
    const ccomplex_t IT_0155 = (-0.5)*IT_0118;
    const ccomplex_t IT_0156 = 6*IT_0155;
    const ccomplex_t IT_0157 = m_c*s_14;
    const ccomplex_t IT_0158 = (-6)*IT_0155;
    const ccomplex_t IT_0159 = 36*s_34;
    const ccomplex_t IT_0160 = m_c*e_em*mu_h*IT_0006*IT_0012;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = cos(beta);
    const ccomplex_t IT_0163 = m_c*e_em*IT_0006*IT_0012*IT_0014*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = IT_0161*IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_12 + -IT_0019 + -IT_0139 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0167 = IT_0165*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = m_c*e_em*IT_0006*IT_0012;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = m_c*e_em*mu_h*IT_0006*IT_0012*IT_0014*IT_0162;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = pow(m_Z, 2);
    const ccomplex_t IT_0175 = cpow((-2)*s_12 + -IT_0019 + -IT_0139 + IT_0174 
      + -reg_prop, -1);
    const ccomplex_t IT_0176 = IT_0173*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = -IT_0168 + -IT_0177;
    const ccomplex_t IT_0179 = 18*IT_0018;
    const ccomplex_t IT_0180 = 18*s_34;
    const ccomplex_t IT_0181 = IT_0179 + IT_0180;
    const ccomplex_t IT_0182 = (-18)*IT_0018;
    const ccomplex_t IT_0183 = IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0185 = 0.166666666666667*IT_0184;
    const ccomplex_t IT_0186 = (-6)*IT_0155*(s_34*(IT_0129 + -conj(IT_0129) +
       3*IT_0155) + 0.166666666666667*(IT_0178 + -conj(IT_0178))*IT_0181 +
       0.166666666666667*(IT_0151 + -conj(IT_0151))*IT_0183 + -IT_0018*(IT_0110 
      + IT_0185));
    const ccomplex_t IT_0187 = s_13 + s_14;
    const ccomplex_t IT_0188 = m_c*IT_0187;
    const ccomplex_t IT_0189 = s_34 + IT_0018;
    const ccomplex_t IT_0190 = IT_0018*IT_0139;
    const ccomplex_t IT_0191 = s_34*IT_0139;
    const ccomplex_t IT_0192 = (-18)*IT_0191;
    const ccomplex_t IT_0193 = s_13*s_14;
    const ccomplex_t IT_0194 = 36*IT_0193;
    const ccomplex_t IT_0195 = IT_0192 + IT_0194;
    const ccomplex_t IT_0196 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0197 = 6*conj(IT_0151);
    const ccomplex_t IT_0198 = IT_0093*IT_0197;
    const ccomplex_t IT_0199 = (-18)*conj(IT_0110);
    const ccomplex_t IT_0200 = (-18)*conj(IT_0065);
    const ccomplex_t IT_0201 = IT_0131*IT_0197;
    const ccomplex_t IT_0202 = ((conj(IT_0065) + conj(IT_0090))*IT_0093 + 
      (conj(IT_0110) + conj(IT_0129))*IT_0131)*IT_0152 + IT_0018*conj(IT_0129)
      *IT_0153 + conj(IT_0090)*(IT_0154*IT_0156 + IT_0153*IT_0157) + conj
      (IT_0065)*IT_0157*IT_0158 + IT_0151*conj(IT_0151)*((-36)*IT_0018 + IT_0159
      ) + IT_0186 + (-6)*IT_0178*((-0.166666666666667)*(36*IT_0018 + IT_0159)
      *conj(IT_0178) + (conj(IT_0065) + -conj(IT_0090))*IT_0188 + -(conj(IT_0129
      ) + IT_0185)*IT_0189) + IT_0065*(18*conj(IT_0110)*IT_0154 + ((-18)*conj
      (IT_0129) + IT_0156)*IT_0157 + (-18)*conj(IT_0090)*IT_0190 + conj(IT_0065)
      *IT_0195 + IT_0188*IT_0196 + IT_0198) + 6*IT_0090*(0.166666666666667
      *IT_0154*(18*conj(IT_0129) + IT_0158) + conj(IT_0178)*IT_0188 +
       0.166666666666667*conj(IT_0090)*IT_0195 + 0.166666666666667*IT_0198 +
       0.166666666666667*IT_0157*IT_0199 + 0.166666666666667*IT_0190*IT_0200) +
       IT_0110*(18*conj(IT_0065)*IT_0154 + conj(IT_0110)*IT_0180 + IT_0189
      *IT_0196 + IT_0201) + IT_0129*(18*conj(IT_0090)*IT_0154 + conj(IT_0129)
      *IT_0180 + 6*conj(IT_0178)*IT_0189 + IT_0018*IT_0199 + IT_0157*IT_0200 +
       IT_0201);
    return create_ccomplex_return(IT_0202);
}

