#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_4_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_4_to_u_anti_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = s_12*s_34;
    const ccomplex_t IT_0002 = s_14*s_23;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = IT_0000 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = N_B4*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = N_W4*e_em;
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0017 + (-3)*IT_0020);
    const ccomplex_t IT_0022 = 0.166666666666667*IT_0021;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = pow(m_u, 2);
    const ccomplex_t IT_0025 = pow(m_C_2, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_N_4, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0025 + -IT_0029 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0006*IT_0008*IT_0010*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0010*IT_0031;
    const ccomplex_t IT_0035 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = 2*IT_0036;
    const ccomplex_t IT_0038 = IT_0014*IT_0031;
    const ccomplex_t IT_0039 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = IT_0041 + IT_0044;
    const ccomplex_t IT_0046 = IT_0037 + IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0033*IT_0048;
    const ccomplex_t IT_0050 = IT_0030*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = pow(m_W, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + -IT_0025 + -IT_0029 + IT_0052 
      + -reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0006*IT_0010;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0007*IT_0010;
    const ccomplex_t IT_0057 = IT_0035*IT_0056;
    const ccomplex_t IT_0058 = 2*IT_0057;
    const ccomplex_t IT_0059 = IT_0007*IT_0014;
    const ccomplex_t IT_0060 = IT_0039*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = IT_0042*IT_0056;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = IT_0058 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0055*IT_0067;
    const ccomplex_t IT_0069 = IT_0053*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0051 + IT_0070;
    const ccomplex_t IT_0072 = 6*IT_0071;
    const ccomplex_t IT_0073 = 6*conj(IT_0071);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0075 = cpow(IT_0031, -1);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0006*IT_0010*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0074*IT_0077;
    const ccomplex_t IT_0079 = IT_0026*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0006*IT_0010*IT_0075;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0014;
    const ccomplex_t IT_0085 = (-0.333333333333333)*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0081*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0080 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0006*IT_0010;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0093 = IT_0034*IT_0092;
    const ccomplex_t IT_0094 = (-2)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0096 = IT_0038*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0099 = IT_0034*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = IT_0053*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0056*IT_0092;
    const ccomplex_t IT_0109 = (-2)*IT_0108;
    const ccomplex_t IT_0110 = IT_0059*IT_0095;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = IT_0056*IT_0098;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = IT_0109 + IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0006*IT_0007*IT_0010*IT_0075;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = IT_0030*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0107 + -IT_0122;
    const ccomplex_t IT_0124 = 6*IT_0123;
    const ccomplex_t IT_0125 = 6*conj(IT_0123);
    const ccomplex_t IT_0126 = IT_0055*IT_0104;
    const ccomplex_t IT_0127 = IT_0053*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0033*IT_0117;
    const ccomplex_t IT_0130 = IT_0030*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = -IT_0128 + -IT_0131;
    const ccomplex_t IT_0133 = 12*IT_0132;
    const ccomplex_t IT_0134 = 12*IT_0071;
    const ccomplex_t IT_0135 = IT_0048*IT_0119;
    const ccomplex_t IT_0136 = IT_0030*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0067*IT_0091;
    const ccomplex_t IT_0139 = IT_0053*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0137 + IT_0140;
    const ccomplex_t IT_0142 = 12*IT_0141;
    const ccomplex_t IT_0143 = 12*IT_0123;
    const ccomplex_t IT_0144 = m_u*s_24*m_C_2;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0010;
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = conj(N_B4)*e_em;
    const ccomplex_t IT_0148 = IT_0014*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W4)*e_em;
    const ccomplex_t IT_0151 = IT_0010*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + 3*IT_0152);
    const ccomplex_t IT_0154 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0155 = IT_0146*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_23 + IT_0024 + IT_0029 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0012*IT_0160;
    const ccomplex_t IT_0162 = cpow((-2)*s_23 + IT_0024 + IT_0029 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0158 + -IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0025 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0167 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0168 = IT_0010*IT_0167;
    const ccomplex_t IT_0169 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0170 = IT_0010*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + (-0.5)*IT_0171);
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0010;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = IT_0173*IT_0175;
    const ccomplex_t IT_0177 = IT_0166*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0012*IT_0077;
    const ccomplex_t IT_0180 = IT_0026*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = IT_0083*IT_0154;
    const ccomplex_t IT_0184 = IT_0156*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0124;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 0.166666666666667*IT_0125;
    const ccomplex_t IT_0189 = conj(IT_0185) + IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0006*IT_0010*IT_0075;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = IT_0083*IT_0191;
    const ccomplex_t IT_0193 = IT_0081*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0022*IT_0074;
    const ccomplex_t IT_0196 = IT_0026*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = -IT_0194 + -IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0014;
    const ccomplex_t IT_0200 = 0.666666666666667*IT_0199;
    const ccomplex_t IT_0201 = IT_0012*IT_0200;
    const ccomplex_t IT_0202 = IT_0162*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0146*IT_0205;
    const ccomplex_t IT_0207 = IT_0156*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0203 + IT_0208;
    const ccomplex_t IT_0210 = 6*IT_0209;
    const ccomplex_t IT_0211 = IT_0072 + IT_0210;
    const ccomplex_t IT_0212 = 6*conj(IT_0209);
    const ccomplex_t IT_0213 = IT_0073 + IT_0212;
    const ccomplex_t IT_0214 = 6*conj(IT_0132);
    const ccomplex_t IT_0215 = IT_0083*IT_0205;
    const ccomplex_t IT_0216 = IT_0156*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = 6*IT_0141;
    const ccomplex_t IT_0220 = 6*conj(IT_0141);
    const ccomplex_t IT_0221 = 12*IT_0198;
    const ccomplex_t IT_0222 = 12*conj(IT_0198);
    const ccomplex_t IT_0223 = 12*IT_0182;
    const ccomplex_t IT_0224 = 12*conj(IT_0182);
    const ccomplex_t IT_0225 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0226 = IT_0010*IT_0225;
    const ccomplex_t IT_0227 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0228 = IT_0010*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*(IT_0226 + 0.5*IT_0229);
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = IT_0175*IT_0231;
    const ccomplex_t IT_0233 = IT_0166*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = 12*IT_0209;
    const ccomplex_t IT_0236 = 12*conj(IT_0071);
    const ccomplex_t IT_0237 = 12*conj(IT_0209);
    const ccomplex_t IT_0238 = IT_0132*conj(IT_0165) + 2*conj(IT_0132)*IT_0178
       + 0.166666666666667*IT_0133*conj(IT_0178) + conj(IT_0182)*IT_0187 +
       IT_0182*IT_0189 + 0.166666666666667*conj(IT_0198)*IT_0211 +
       0.166666666666667*IT_0198*IT_0213 + 0.166666666666667*IT_0165*IT_0214 +
       0.166666666666667*conj(IT_0218)*IT_0219 + 0.166666666666667*IT_0218
      *IT_0220 + 0.166666666666667*conj(IT_0028)*IT_0221 + 0.166666666666667
      *IT_0028*IT_0222 + 0.166666666666667*conj(IT_0089)*IT_0223 +
       0.166666666666667*IT_0089*IT_0224 + 4*conj(IT_0234)*(IT_0028 +
       0.0416666666666667*IT_0134 + 0.0416666666666667*IT_0235) + 4*IT_0234*
      (conj(IT_0028) + 0.0416666666666667*IT_0236 + 0.0416666666666667*IT_0237);
    const ccomplex_t IT_0239 = 24*IT_0234;
    const ccomplex_t IT_0240 = IT_0221 + IT_0239;
    const ccomplex_t IT_0241 = 12*IT_0089;
    const ccomplex_t IT_0242 = 12*IT_0028;
    const ccomplex_t IT_0243 = m_d*s_13*m_N_4;
    const ccomplex_t IT_0244 = 6*IT_0198;
    const ccomplex_t IT_0245 = 6*conj(IT_0198);
    const ccomplex_t IT_0246 = 6*IT_0218;
    const ccomplex_t IT_0247 = 6*conj(IT_0218);
    const ccomplex_t IT_0248 = 12*conj(IT_0141);
    const ccomplex_t IT_0249 = 12*IT_0185;
    const ccomplex_t IT_0250 = IT_0142*conj(IT_0178) + conj(IT_0182)*IT_0211 +
       IT_0182*IT_0213 + conj(IT_0165)*IT_0219 + IT_0165*IT_0220 + conj(IT_0089)
      *IT_0221 + IT_0089*IT_0222 + conj(IT_0028)*IT_0223 + IT_0028*IT_0224 + (24
      *conj(IT_0089) + 12*conj(IT_0123) + 12*conj(IT_0185))*IT_0234 + (conj
      (IT_0123) + conj(IT_0185))*IT_0244 + (IT_0123 + IT_0185)*IT_0245 + conj
      (IT_0132)*IT_0246 + IT_0132*IT_0247 + IT_0178*IT_0248 + conj(IT_0234)*(24
      *IT_0089 + IT_0143 + IT_0249);
    const ccomplex_t IT_0251 = 12*IT_0165;
    const ccomplex_t IT_0252 = 24*IT_0178;
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 12*IT_0218;
    const ccomplex_t IT_0255 = s_34*m_C_2*m_N_4;
    const ccomplex_t IT_0256 = 12*conj(IT_0165);
    const ccomplex_t IT_0257 = (-12)*IT_0071;
    const ccomplex_t IT_0258 = (-12)*conj(IT_0071);
    const ccomplex_t IT_0259 = (-12)*IT_0141;
    const ccomplex_t IT_0260 = (-12)*conj(IT_0141);
    const ccomplex_t IT_0261 = (-6)*IT_0132;
    const ccomplex_t IT_0262 = (-6)*conj(IT_0132);
    const ccomplex_t IT_0263 = (-6)*IT_0141;
    const ccomplex_t IT_0264 = (-6)*conj(IT_0141);
    const ccomplex_t IT_0265 = conj(IT_0185)*IT_0219 + IT_0185*IT_0220 + conj
      (IT_0178)*IT_0221 + IT_0178*IT_0222 + conj(IT_0165)*IT_0244 + IT_0165
      *IT_0245 + conj(IT_0182)*IT_0246 + IT_0182*IT_0247 + conj(IT_0234)*IT_0253
       + IT_0234*(24*conj(IT_0178) + IT_0256) + conj(IT_0132)*(IT_0210 + IT_0257
      ) + IT_0132*(IT_0212 + IT_0258) + conj(IT_0123)*IT_0259 + IT_0123*IT_0260 
      + conj(IT_0028)*IT_0261 + IT_0028*IT_0262 + conj(IT_0089)*IT_0263 +
       IT_0089*IT_0264;
    const ccomplex_t IT_0266 = -IT_0000;
    const ccomplex_t IT_0267 = IT_0001 + IT_0002 + IT_0266;
    const ccomplex_t IT_0268 = (-6)*IT_0071;
    const ccomplex_t IT_0269 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0270 = (-6)*IT_0123;
    const ccomplex_t IT_0271 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0272 = m_d*m_u*m_C_2*m_N_4;
    const ccomplex_t IT_0273 = IT_0235 + IT_0268;
    const ccomplex_t IT_0274 = IT_0237 + IT_0269;
    const ccomplex_t IT_0275 = 12*conj(IT_0218);
    const ccomplex_t IT_0276 = 12*conj(IT_0089);
    const ccomplex_t IT_0277 = (-6)*IT_0209;
    const ccomplex_t IT_0278 = (-6)*conj(IT_0209);
    const ccomplex_t IT_0279 = conj(IT_0132)*IT_0142 + conj(IT_0089)*IT_0211 +
       IT_0089*IT_0213 + 24*conj(IT_0178)*IT_0218 + 24*IT_0178*conj(IT_0218) +
       IT_0182*(IT_0222 + 24*conj(IT_0234)) + conj(IT_0182)*IT_0240 + 6*conj
      (IT_0028)*(IT_0187 + 0.166666666666667*IT_0241) + IT_0132*IT_0248 + conj
      (IT_0165)*IT_0254 + conj(IT_0185)*IT_0273 + IT_0185*IT_0274 + IT_0165
      *IT_0275 + 6*IT_0028*(IT_0189 + 0.166666666666667*IT_0276) + conj(IT_0123)
      *(IT_0134 + IT_0277) + IT_0123*(IT_0236 + IT_0278);
    const ccomplex_t IT_0280 = m_d*s_23*m_C_2;
    const ccomplex_t IT_0281 = 6*IT_0165;
    const ccomplex_t IT_0282 = 6*conj(IT_0165);
    const ccomplex_t IT_0283 = (-12)*conj(IT_0123);
    const ccomplex_t IT_0284 = conj(IT_0178)*((-12)*IT_0123 + 24*IT_0185 +
       IT_0241) + conj(IT_0028)*IT_0246 + IT_0028*IT_0247 + conj(IT_0185)
      *IT_0251 + IT_0185*IT_0256 + conj(IT_0234)*IT_0259 + IT_0234*IT_0260 +
       conj(IT_0182)*IT_0261 + IT_0182*IT_0262 + conj(IT_0198)*IT_0263 + IT_0198
      *IT_0264 + conj(IT_0165)*IT_0270 + IT_0165*IT_0271 + conj(IT_0218)*IT_0273
       + IT_0218*IT_0274 + conj(IT_0089)*IT_0281 + IT_0089*IT_0282 + IT_0178*(24
      *conj(IT_0185) + IT_0276 + IT_0283);
    const ccomplex_t IT_0285 = m_d*m_u*s_12;
    const ccomplex_t IT_0286 = (-12)*conj(IT_0132);
    const ccomplex_t IT_0287 = IT_0185*IT_0214 + conj(IT_0209)*IT_0219 +
       IT_0209*IT_0220 + conj(IT_0178)*IT_0223 + IT_0178*IT_0224 + conj(IT_0198)
      *IT_0246 + IT_0198*IT_0247 + conj(IT_0234)*IT_0254 + conj(IT_0141)*IT_0257
       + IT_0141*IT_0258 + conj(IT_0089)*IT_0261 + IT_0089*IT_0262 + conj
      (IT_0028)*IT_0263 + IT_0028*IT_0264 + IT_0234*IT_0275 + conj(IT_0182)
      *IT_0281 + IT_0182*IT_0282 + 6*IT_0132*(conj(IT_0185) + 0.166666666666667
      *IT_0283) + IT_0123*IT_0286;
    const ccomplex_t IT_0288 = m_u*s_14*m_N_4;
    const ccomplex_t IT_0289 = (-0.166666666666667)*IT_0028;
    const ccomplex_t IT_0290 = (-0.166666666666667)*conj(IT_0028);
    const ccomplex_t IT_0291 = conj(IT_0132)*IT_0198 + conj(IT_0123)*IT_0218 +
       IT_0123*conj(IT_0218) + IT_0132*(conj(IT_0198) + 2*conj(IT_0234)) + (
      -0.166666666666667)*conj(IT_0089)*IT_0246 + (-0.166666666666667)*IT_0089
      *IT_0247 + (-0.166666666666667)*conj(IT_0185)*IT_0254 + (-4)*conj(IT_0178)
      *(IT_0209 + 0.0416666666666667*IT_0242 + 0.0416666666666667*IT_0257) + (-4
      )*IT_0178*(0.5*conj(IT_0028) + conj(IT_0209) + 0.0416666666666667*IT_0258)
       + (-0.166666666666667)*conj(IT_0182)*IT_0263 + (-0.166666666666667)
      *IT_0182*IT_0264 + (-0.166666666666667)*conj(IT_0165)*IT_0273 + (
      -0.166666666666667)*IT_0165*IT_0274 + (-0.166666666666667)*IT_0185*IT_0275
       + (-0.166666666666667)*IT_0234*IT_0286 + IT_0282*IT_0289 + IT_0281*IT_0290;
    const ccomplex_t IT_0292 = IT_0004 + IT_0266;
    const ccomplex_t IT_0293 = IT_0005*(conj(IT_0028)*IT_0072 + IT_0028
      *IT_0073 + conj(IT_0089)*IT_0124 + IT_0089*IT_0125) + IT_0001*(conj
      (IT_0132)*IT_0133 + conj(IT_0071)*IT_0134 + conj(IT_0141)*IT_0142 + conj
      (IT_0123)*IT_0143) + 6*IT_0144*IT_0238 + IT_0000*(conj(IT_0182)*IT_0223 + 
      (24*IT_0198 + 48*IT_0234)*conj(IT_0234) + conj(IT_0198)*IT_0240 + conj
      (IT_0089)*IT_0241 + conj(IT_0028)*IT_0242) + IT_0243*IT_0250 + IT_0002*(
      (24*IT_0165 + 48*IT_0178)*conj(IT_0178) + conj(IT_0209)*IT_0235 + conj
      (IT_0185)*IT_0249 + conj(IT_0165)*IT_0253 + conj(IT_0218)*IT_0254) +
       IT_0255*IT_0265 + IT_0267*(conj(IT_0209)*IT_0268 + IT_0209*IT_0269 + conj
      (IT_0185)*IT_0270 + IT_0185*IT_0271) + IT_0272*IT_0279 + IT_0280*IT_0284 +
       IT_0285*IT_0287 + (-6)*IT_0288*IT_0291 + (-6)*(conj(IT_0089)*IT_0185 +
       IT_0089*conj(IT_0185) + IT_0278*IT_0289 + IT_0277*IT_0290)*IT_0292;
    return create_ccomplex_return(IT_0293);
}

