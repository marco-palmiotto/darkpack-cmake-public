#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0002;
    const ccomplex_t IT_0005 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0002, -1);
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009 + (-2)
      *IT_0011);
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = e_em*IT_0004;
    const ccomplex_t IT_0020 = e_em*IT_0008;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0019 + -IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_tau, 2);
    const ccomplex_t IT_0027 = s_12*IT_0026;
    const ccomplex_t IT_0028 = m_C_2*m_tau;
    const ccomplex_t IT_0029 = -s_24;
    const ccomplex_t IT_0030 = s_23 + IT_0029;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 4*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = cos(alpha);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*e_em*m_tau*IT_0007*IT_0033
      *IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = U_d2*V_Wp2;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = sin(alpha);
    const ccomplex_t IT_0042 = V_u2*U_Wm2;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0007*(IT_0040 + IT_0043);
    const ccomplex_t IT_0045 = 1.4142135623731*e_em*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0038*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*m_tau*IT_0007*IT_0033
      *IT_0035*IT_0041;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0036*IT_0042;
    const ccomplex_t IT_0055 = IT_0039*IT_0041;
    const ccomplex_t IT_0056 = IT_0007*(IT_0054 + -IT_0055);
    const ccomplex_t IT_0057 = 1.4142135623731*e_em*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0053*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = -IT_0051 + -IT_0063;
    const ccomplex_t IT_0065 = -s_14;
    const ccomplex_t IT_0066 = s_13 + IT_0065;
    const ccomplex_t IT_0067 = IT_0028*IT_0066;
    const ccomplex_t IT_0068 = (-4)*IT_0067;
    const ccomplex_t IT_0069 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0070 = IT_0036*IT_0069;
    const ccomplex_t IT_0071 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0072 = IT_0041*IT_0071;
    const ccomplex_t IT_0073 = IT_0007*(IT_0070 + IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0038*IT_0076;
    const ccomplex_t IT_0078 = IT_0049*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0036*IT_0071;
    const ccomplex_t IT_0081 = IT_0041*IT_0069;
    const ccomplex_t IT_0082 = IT_0007*(IT_0080 + -IT_0081);
    const ccomplex_t IT_0083 = 1.4142135623731*e_em*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0053*IT_0085;
    const ccomplex_t IT_0087 = IT_0061*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = -IT_0079 + -IT_0088;
    const ccomplex_t IT_0090 = IT_0015*IT_0026;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0007;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0007;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = cpow(s_13 + (-0.5)*IT_0015 + (-0.5)*IT_0026 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = s_13*s_24;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0007*IT_0033*IT_0035;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0007*IT_0033*IT_0035;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = IT_0095*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = s_23 + s_24;
    const ccomplex_t IT_0109 = IT_0028*IT_0108;
    const ccomplex_t IT_0110 = sin(beta);
    const ccomplex_t IT_0111 = e_em*m_tau*IT_0007*IT_0033*IT_0035*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0034*IT_0042;
    const ccomplex_t IT_0114 = IT_0039*IT_0110;
    const ccomplex_t IT_0115 = IT_0007*(IT_0113 + IT_0114);
    const ccomplex_t IT_0116 = 1.4142135623731*e_em*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0112*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_12 + (-2)*IT_0015 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = pow(m_Z, 2);
    const ccomplex_t IT_0123 = cpow((-2)*s_12 + (-2)*IT_0015 + IT_0122 + 
      -reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0034*IT_0039;
    const ccomplex_t IT_0125 = IT_0042*IT_0110;
    const ccomplex_t IT_0126 = IT_0007*(IT_0124 + -IT_0125);
    const ccomplex_t IT_0127 = 1.4142135623731*e_em*IT_0126;
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = e_em*m_tau*IT_0007*IT_0033;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = IT_0123*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0121 + -IT_0133;
    const ccomplex_t IT_0135 = s_13 + s_14;
    const ccomplex_t IT_0136 = IT_0028*IT_0135;
    const ccomplex_t IT_0137 = IT_0034*IT_0071;
    const ccomplex_t IT_0138 = IT_0069*IT_0110;
    const ccomplex_t IT_0139 = IT_0007*(IT_0137 + IT_0138);
    const ccomplex_t IT_0140 = 1.4142135623731*e_em*IT_0139;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0112*IT_0141;
    const ccomplex_t IT_0143 = IT_0119*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0034*IT_0069;
    const ccomplex_t IT_0146 = IT_0071*IT_0110;
    const ccomplex_t IT_0147 = IT_0007*(IT_0145 + -IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0130*IT_0149;
    const ccomplex_t IT_0151 = IT_0123*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = -IT_0144 + -IT_0152;
    const ccomplex_t IT_0154 = 8*IT_0018*(conj(IT_0025)*IT_0027 + 0.125
      *IT_0032*conj(IT_0064) + 0.125*IT_0068*conj(IT_0089) + -IT_0090*conj
      (IT_0098) + -IT_0099*conj(IT_0107) + 0.5*IT_0109*conj(IT_0134) + (-0.5)
      *IT_0136*conj(IT_0153));
    const ccomplex_t IT_0155 = s_14*s_23;
    const ccomplex_t IT_0156 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0157 = IT_0004*IT_0156;
    const ccomplex_t IT_0158 = IT_0008*IT_0156;
    const ccomplex_t IT_0159 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0160 = IT_0008*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + -IT_0158 + (-2)
      *IT_0160);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0003*IT_0162;
    const ccomplex_t IT_0164 = IT_0016*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = s_34*IT_0015;
    const ccomplex_t IT_0168 = IT_0022*IT_0162;
    const ccomplex_t IT_0169 = IT_0016*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = s_23*m_C_2*m_tau;
    const ccomplex_t IT_0173 = IT_0092*IT_0101;
    const ccomplex_t IT_0174 = IT_0095*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = s_14*m_C_2*m_tau;
    const ccomplex_t IT_0178 = IT_0093*IT_0103;
    const ccomplex_t IT_0179 = IT_0095*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0180;
    const ccomplex_t IT_0182 = 4*IT_0025;
    const ccomplex_t IT_0183 = s_24*m_C_2*m_tau;
    const ccomplex_t IT_0184 = s_13*m_C_2*m_tau;
    const ccomplex_t IT_0185 = 8*IT_0018;
    const ccomplex_t IT_0186 = 16*IT_0018;
    const ccomplex_t IT_0187 = (-4)*IT_0025;
    const ccomplex_t IT_0188 = s_12*s_34;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = -IT_0099;
    const ccomplex_t IT_0191 = IT_0027 + IT_0155 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 2*IT_0191;
    const ccomplex_t IT_0193 = (-2)*IT_0027 + 2*IT_0155 + 2*IT_0189 + 2*IT_0190;
    const ccomplex_t IT_0194 = s_34 + IT_0026;
    const ccomplex_t IT_0195 = IT_0015*IT_0194;
    const ccomplex_t IT_0196 = 2*IT_0195;
    const ccomplex_t IT_0197 = -IT_0026;
    const ccomplex_t IT_0198 = s_34 + IT_0197;
    const ccomplex_t IT_0199 = IT_0015*IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = 4*conj(IT_0025);
    const ccomplex_t IT_0202 = 4*IT_0090;
    const ccomplex_t IT_0203 = 4*IT_0172;
    const ccomplex_t IT_0204 = 4*IT_0099;
    const ccomplex_t IT_0205 = 8*conj(IT_0018);
    const ccomplex_t IT_0206 = 8*IT_0183;
    const ccomplex_t IT_0207 = (-4)*IT_0184;
    const ccomplex_t IT_0208 = (-4)*IT_0183;
    const ccomplex_t IT_0209 = conj(IT_0089)*IT_0192 + conj(IT_0153)*IT_0193 +
       conj(IT_0134)*IT_0196 + conj(IT_0064)*IT_0200 + IT_0177*IT_0201 + conj
      (IT_0176)*IT_0202 + conj(IT_0166)*IT_0203 + conj(IT_0181)*IT_0204 +
       IT_0184*IT_0205 + conj(IT_0171)*IT_0206 + conj(IT_0107)*IT_0207 + conj
      (IT_0098)*IT_0208;
    const ccomplex_t IT_0210 = s_12*IT_0194;
    const ccomplex_t IT_0211 = 2*IT_0136;
    const ccomplex_t IT_0212 = 2*IT_0109;
    const ccomplex_t IT_0213 = 4*IT_0109;
    const ccomplex_t IT_0214 = (-8)*IT_0195;
    const ccomplex_t IT_0215 = (-4)*IT_0109;
    const ccomplex_t IT_0216 = (-2)*IT_0195;
    const ccomplex_t IT_0217 = 0.125*conj(IT_0181)*IT_0193 + (-0.5)*IT_0136*
      (conj(IT_0018) + (-0.25)*IT_0201) + conj(IT_0153)*IT_0210 + 0.125*conj
      (IT_0107)*IT_0211 + 0.125*conj(IT_0098)*IT_0212 + 0.125*conj(IT_0166)
      *IT_0213 + 0.125*conj(IT_0134)*IT_0214 + 0.125*conj(IT_0171)*IT_0215 +
       0.125*conj(IT_0176)*IT_0216;
    const ccomplex_t IT_0218 = (-8)*conj(IT_0018);
    const ccomplex_t IT_0219 = (-8)*IT_0090;
    const ccomplex_t IT_0220 = (-4)*conj(IT_0025);
    const ccomplex_t IT_0221 = (-4)*IT_0167;
    const ccomplex_t IT_0222 = (-2)*IT_0109;
    const ccomplex_t IT_0223 = (-2)*IT_0067;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = 2*IT_0031;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = conj(IT_0098)*IT_0202 + conj(IT_0107)*IT_0204 +
       conj(IT_0181)*IT_0207 + conj(IT_0176)*IT_0208 + conj(IT_0153)*IT_0211 +
       IT_0099*IT_0218 + conj(IT_0171)*IT_0219 + IT_0027*IT_0220 + conj(IT_0166)
      *IT_0221 + conj(IT_0134)*IT_0222 + conj(IT_0089)*IT_0224 + conj(IT_0064)
      *IT_0226;
    const ccomplex_t IT_0228 = (-8)*IT_0099;
    const ccomplex_t IT_0229 = (-4)*IT_0027;
    const ccomplex_t IT_0230 = (-2)*IT_0136;
    const ccomplex_t IT_0231 = conj(IT_0107)*IT_0202 + conj(IT_0098)*IT_0204 +
       conj(IT_0176)*IT_0207 + conj(IT_0181)*IT_0208 + conj(IT_0153)*IT_0212 +
       IT_0090*IT_0218 + IT_0167*IT_0220 + conj(IT_0064)*IT_0224 + conj(IT_0089)
      *IT_0226 + conj(IT_0171)*IT_0228 + conj(IT_0166)*IT_0229 + conj(IT_0134)
      *IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0232*IT_0233;
    const ccomplex_t IT_0235 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0236 = IT_0234*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = -s_14 + -s_23;
    const ccomplex_t IT_0240 = s_13 + s_24;
    const ccomplex_t IT_0241 = IT_0239 + IT_0240;
    const ccomplex_t IT_0242 = IT_0028*IT_0241;
    const ccomplex_t IT_0243 = (-8)*IT_0242;
    const ccomplex_t IT_0244 = 8*IT_0167;
    const ccomplex_t IT_0245 = 16*IT_0090;
    const ccomplex_t IT_0246 = 8*IT_0027;
    const ccomplex_t IT_0247 = 16*IT_0155;
    const ccomplex_t IT_0248 = IT_0244 + IT_0245 + IT_0246 + IT_0247;
    const ccomplex_t IT_0249 = 2*IT_0099;
    const ccomplex_t IT_0250 = 2*IT_0090;
    const ccomplex_t IT_0251 = IT_0027 + IT_0167 + IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0167;
    const ccomplex_t IT_0253 = 0.5*IT_0027 + IT_0090 + IT_0099 + IT_0252;
    const ccomplex_t IT_0254 = 16*IT_0099;
    const ccomplex_t IT_0255 = IT_0244 + IT_0245 + IT_0246 + IT_0254;
    const ccomplex_t IT_0256 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0257 = IT_0240 + IT_0256;
    const ccomplex_t IT_0258 = IT_0028*IT_0257;
    const ccomplex_t IT_0259 = 8*IT_0258;
    const ccomplex_t IT_0260 = (-8)*IT_0258;
    const ccomplex_t IT_0261 = -IT_0243;
    const ccomplex_t IT_0262 = -(32*IT_0027 + 64*IT_0090 + 32*IT_0099 + 32
      *IT_0155 + 32*IT_0167)*IT_0238 + (conj(IT_0064) + conj(IT_0089))*IT_0243 +
       -(IT_0025 + -conj(IT_0025) + IT_0166 + -conj(IT_0166))*IT_0248 + 4*
      (IT_0107 + -conj(IT_0107))*IT_0251 + 8*(IT_0098 + -conj(IT_0098))*IT_0253 
      + -(IT_0018 + -conj(IT_0018) + IT_0171 + -conj(IT_0171))*IT_0255 + (conj
      (IT_0176) + conj(IT_0181))*IT_0259 + (IT_0176 + IT_0181)*IT_0260 + 
      (IT_0064 + IT_0089)*IT_0261;
    const ccomplex_t IT_0263 = 4*IT_0177;
    const ccomplex_t IT_0264 = (-4)*IT_0136;
    const ccomplex_t IT_0265 = IT_0032*conj(IT_0089);
    const ccomplex_t IT_0266 = conj(IT_0064)*IT_0068;
    const ccomplex_t IT_0267 = 16*conj(IT_0025)*IT_0090 + conj(IT_0181)
      *IT_0203 + IT_0167*IT_0205 + conj(IT_0153)*IT_0213 + conj(IT_0107)*IT_0221
       + conj(IT_0098)*IT_0229 + conj(IT_0171)*IT_0246 + conj(IT_0166)*IT_0247 +
       conj(IT_0176)*IT_0263 + conj(IT_0134)*IT_0264 + IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = 4*IT_0136;
    const ccomplex_t IT_0269 = 8*IT_0184;
    const ccomplex_t IT_0270 = 16*conj(IT_0018)*IT_0090 + 8*conj(IT_0025)
      *IT_0167 + conj(IT_0181)*IT_0206 + conj(IT_0153)*IT_0215 + conj(IT_0107)
      *IT_0219 + conj(IT_0098)*IT_0228 + conj(IT_0166)*IT_0246 + conj(IT_0171)
      *IT_0254 + IT_0265 + IT_0266 + conj(IT_0134)*IT_0268 + conj(IT_0176)
      *IT_0269;
    const ccomplex_t IT_0271 = conj(IT_0166) + conj(IT_0171);
    const ccomplex_t IT_0272 = (-8)*IT_0199;
    const ccomplex_t IT_0273 = s_12*IT_0198;
    const ccomplex_t IT_0274 = 8*IT_0273;
    const ccomplex_t IT_0275 = -IT_0027;
    const ccomplex_t IT_0276 = IT_0155 + IT_0189 + IT_0190 + IT_0275;
    const ccomplex_t IT_0277 = (-0.5)*conj(IT_0176);
    const ccomplex_t IT_0278 = 0.125*conj(IT_0181)*IT_0196 + conj(IT_0134)
      *IT_0210 + 0.125*conj(IT_0153)*IT_0214 + 0.5*IT_0109*(conj(IT_0018) + 0.25
      *IT_0220) + 0.125*conj(IT_0107)*IT_0222 + 0.125*conj(IT_0098)*IT_0230 +
       0.125*conj(IT_0166)*IT_0264 + 0.125*conj(IT_0171)*IT_0268 + 0.5*IT_0276
      *IT_0277;
    const ccomplex_t IT_0279 = (-0.5)*conj(IT_0064)*IT_0192 + (-0.5)*conj
      (IT_0089)*IT_0200 + (-0.5)*IT_0172*IT_0201 + (-0.5)*conj(IT_0181)*IT_0202 
      + (-0.5)*IT_0183*IT_0205 + (-0.5)*conj(IT_0098)*IT_0207 + (-0.5)*conj
      (IT_0107)*IT_0208 + (-0.5)*conj(IT_0153)*IT_0216 + (-0.5)*conj(IT_0166)
      *IT_0263 + (-0.5)*conj(IT_0171)*IT_0269 + conj(IT_0134)*IT_0276 + IT_0204
      *IT_0277;
    const ccomplex_t IT_0280 = IT_0154 + 8*IT_0025*(conj(IT_0018)*IT_0027 + 1.
      /8*IT_0032*conj(IT_0064) + 1./8*IT_0068*conj(IT_0089) + 2*conj(IT_0025)
      *IT_0155 + 2*IT_0090*conj(IT_0166) + IT_0167*conj(IT_0171)) + (IT_0136
      *conj(IT_0153) + IT_0172*conj(IT_0176) + IT_0177*conj(IT_0181))*IT_0182 + 
      (conj(IT_0166)*IT_0167 + conj(IT_0176)*IT_0183 + conj(IT_0181)*IT_0184)
      *IT_0185 + (conj(IT_0018)*IT_0099 + IT_0090*conj(IT_0171))*IT_0186 + 
      (IT_0027*conj(IT_0107) + IT_0109*conj(IT_0134) + conj(IT_0098)*IT_0167)
      *IT_0187 + IT_0181*IT_0209 + 8*IT_0153*IT_0217 + IT_0107*IT_0227 + IT_0098
      *IT_0231 + IT_0238*IT_0262 + IT_0166*IT_0267 + IT_0171*IT_0270 + IT_0064*(
      (conj(IT_0018) + conj(IT_0025))*IT_0032 + conj(IT_0176)*IT_0192 + conj
      (IT_0181)*IT_0200 + conj(IT_0098)*IT_0224 + conj(IT_0107)*IT_0226 +
       IT_0068*IT_0271 + conj(IT_0089)*IT_0272 + conj(IT_0064)*IT_0274) +
       IT_0089*((conj(IT_0018) + conj(IT_0025))*IT_0068 + conj(IT_0181)*IT_0192 
      + conj(IT_0176)*IT_0200 + conj(IT_0107)*IT_0224 + conj(IT_0098)*IT_0226 +
       IT_0032*IT_0271 + conj(IT_0064)*IT_0272 + conj(IT_0089)*IT_0274) + 8
      *IT_0134*IT_0278 + (-2)*IT_0176*IT_0279;
    return create_ccomplex_return(IT_0280);
}

