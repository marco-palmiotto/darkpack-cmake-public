#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_C_1, 2);
    const ccomplex_t IT_0001 = pow(m_tau, 2);
    const ccomplex_t IT_0002 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0003 = cpow(s_12 + 0.5*IT_0000 + (-0.5)*IT_0001 + 0.5
      *IT_0002 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = e_em*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0008;
    const ccomplex_t IT_0013 = e_em*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0011 + -IT_0013);
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0002 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0010*IT_0023;
    const ccomplex_t IT_0025 = IT_0012*IT_0023;
    const ccomplex_t IT_0026 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0025 + (-2)
      *IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0002 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = IT_0012*IT_0034;
    const ccomplex_t IT_0037 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0038 = IT_0012*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0036 + (-2)
      *IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (-0.5)*IT_0032 + (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0019 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0009;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = -IT_0018 + (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = s_13*s_34;
    const ccomplex_t IT_0054 = s_14*IT_0001;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0056;
    const ccomplex_t IT_0058 = conj(IT_0052)*IT_0057;
    const ccomplex_t IT_0059 = pow(m_W, -1);
    const ccomplex_t IT_0060 = cos(beta);
    const ccomplex_t IT_0061 = cpow(IT_0060, -1);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0005*IT_0059*IT_0061;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0009;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0070 = IT_0010*IT_0069;
    const ccomplex_t IT_0071 = IT_0012*IT_0069;
    const ccomplex_t IT_0072 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0073 = IT_0012*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + -IT_0071 + (-2)
      *IT_0073);
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0063*IT_0075;
    const ccomplex_t IT_0077 = IT_0033*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0005*IT_0059*IT_0061;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0082 = IT_0010*IT_0081;
    const ccomplex_t IT_0083 = IT_0012*IT_0081;
    const ccomplex_t IT_0084 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0085 = IT_0012*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + -IT_0083 + (-2)
      *IT_0085);
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0087;
    const ccomplex_t IT_0089 = IT_0020*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0078 + (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0078 + IT_0090;
    const ccomplex_t IT_0093 = s_13*s_14;
    const ccomplex_t IT_0094 = s_34*IT_0000;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0093 + IT_0095;
    const ccomplex_t IT_0097 = 2*IT_0096;
    const ccomplex_t IT_0098 = conj(IT_0092)*IT_0097;
    const ccomplex_t IT_0099 = IT_0048*IT_0063;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0047*IT_0100;
    const ccomplex_t IT_0102 = IT_0067 + 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0057*conj(IT_0102);
    const ccomplex_t IT_0104 = (-0.5)*IT_0101;
    const ccomplex_t IT_0105 = (-2)*IT_0096;
    const ccomplex_t IT_0106 = conj(IT_0104)*IT_0105;
    const ccomplex_t IT_0107 = IT_0015*IT_0063;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0109 = IT_0003*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = m_tau*IT_0067;
    const ccomplex_t IT_0111 = (-0.5)*IT_0109 + 0.5*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0113 = IT_0007*IT_0075;
    const ccomplex_t IT_0114 = IT_0033*IT_0112*IT_0113;
    const ccomplex_t IT_0115 = m_C_1*IT_0043;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0117 = IT_0022*IT_0087;
    const ccomplex_t IT_0118 = IT_0020*IT_0116*IT_0117;
    const ccomplex_t IT_0119 = m_C_1*IT_0032;
    const ccomplex_t IT_0120 = (-0.5)*IT_0114 + (-0.5)*IT_0115 + (-0.5)
      *IT_0118 + (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0111 + IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0051;
    const ccomplex_t IT_0123 = pow(m_Z, -2);
    const ccomplex_t IT_0124 = s_14*s_34*IT_0123;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = s_13 + IT_0125;
    const ccomplex_t IT_0127 = m_C_1*IT_0126;
    const ccomplex_t IT_0128 = 2*IT_0127;
    const ccomplex_t IT_0129 = conj(IT_0122)*IT_0128;
    const ccomplex_t IT_0130 = IT_0032 + IT_0043;
    const ccomplex_t IT_0131 = (-2)*IT_0127;
    const ccomplex_t IT_0132 = conj(IT_0130)*IT_0131;
    const ccomplex_t IT_0133 = pow(s_34, 2);
    const ccomplex_t IT_0134 = IT_0123*IT_0133;
    const ccomplex_t IT_0135 = -IT_0001;
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = m_C_1*IT_0136;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = conj(IT_0052)*IT_0138;
    const ccomplex_t IT_0140 = pow(s_14, 2);
    const ccomplex_t IT_0141 = IT_0123*IT_0140;
    const ccomplex_t IT_0142 = -IT_0000;
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = s_13*IT_0143;
    const ccomplex_t IT_0145 = 2*IT_0144;
    const ccomplex_t IT_0146 = s_13*IT_0126;
    const ccomplex_t IT_0147 = 2*IT_0146;
    const ccomplex_t IT_0148 = m_C_1*m_tau;
    const ccomplex_t IT_0149 = IT_0143*IT_0148;
    const ccomplex_t IT_0150 = 2*IT_0149;
    const ccomplex_t IT_0151 = IT_0126*IT_0148;
    const ccomplex_t IT_0152 = 2*IT_0151;
    const ccomplex_t IT_0153 = (-2)*IT_0144;
    const ccomplex_t IT_0154 = (-2)*IT_0149;
    const ccomplex_t IT_0155 = m_tau*IT_0143;
    const ccomplex_t IT_0156 = 2*IT_0155;
    const ccomplex_t IT_0157 = s_13*IT_0136;
    const ccomplex_t IT_0158 = 2*IT_0157;
    const ccomplex_t IT_0159 = IT_0136*IT_0148;
    const ccomplex_t IT_0160 = 2*IT_0159;
    const ccomplex_t IT_0161 = (-2)*IT_0146;
    const ccomplex_t IT_0162 = (-2)*IT_0151;
    const ccomplex_t IT_0163 = m_tau*IT_0126;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = (-2)*IT_0155;
    const ccomplex_t IT_0166 = IT_0007*IT_0064;
    const ccomplex_t IT_0167 = IT_0003*IT_0108*IT_0166;
    const ccomplex_t IT_0168 = m_tau*IT_0018;
    const ccomplex_t IT_0169 = 0.5*IT_0167 + (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0029*IT_0080;
    const ccomplex_t IT_0171 = IT_0020*IT_0116*IT_0170;
    const ccomplex_t IT_0172 = IT_0040*IT_0063;
    const ccomplex_t IT_0173 = IT_0033*IT_0112*IT_0172;
    const ccomplex_t IT_0174 = m_C_1*IT_0078;
    const ccomplex_t IT_0175 = m_C_1*IT_0090;
    const ccomplex_t IT_0176 = (-0.5)*IT_0171 + (-0.5)*IT_0173 + (-0.5)
      *IT_0174 + (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0169 + IT_0176;
    const ccomplex_t IT_0178 = conj(IT_0104)*IT_0128;
    const ccomplex_t IT_0179 = conj(IT_0130)*IT_0156;
    const ccomplex_t IT_0180 = conj(IT_0052)*IT_0164;
    const ccomplex_t IT_0181 = (-6)*IT_0148;
    const ccomplex_t IT_0182 = conj(IT_0111)*IT_0181;
    const ccomplex_t IT_0183 = conj(IT_0092)*IT_0131;
    const ccomplex_t IT_0184 = conj(IT_0102)*IT_0138;
    const ccomplex_t IT_0185 = conj(IT_0122)*IT_0165;
    const ccomplex_t IT_0186 = IT_0130*IT_0156;
    const ccomplex_t IT_0187 = IT_0052*IT_0164;
    const ccomplex_t IT_0188 = pow(m_Z, 2);
    const ccomplex_t IT_0189 = s_34*IT_0188;
    const ccomplex_t IT_0190 = IT_0123*IT_0189;
    const ccomplex_t IT_0191 = (-0.25)*IT_0190;
    const ccomplex_t IT_0192 = s_34 + IT_0191;
    const ccomplex_t IT_0193 = m_C_1*IT_0192;
    const ccomplex_t IT_0194 = (-4)*IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0194;
    const ccomplex_t IT_0196 = s_14*IT_0188;
    const ccomplex_t IT_0197 = IT_0123*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = s_14 + IT_0198;
    const ccomplex_t IT_0200 = m_tau*IT_0199;
    const ccomplex_t IT_0201 = 2*IT_0200;
    const ccomplex_t IT_0202 = 2*IT_0201;
    const ccomplex_t IT_0203 = s_13*IT_0188;
    const ccomplex_t IT_0204 = s_14*s_34;
    const ccomplex_t IT_0205 = (-2)*IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = IT_0123*IT_0206;
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = s_13 + IT_0208;
    const ccomplex_t IT_0210 = 4*IT_0209;
    const ccomplex_t IT_0211 = IT_0169*IT_0210;
    const ccomplex_t IT_0212 = IT_0092*IT_0131;
    const ccomplex_t IT_0213 = IT_0102*IT_0138;
    const ccomplex_t IT_0214 = IT_0122*IT_0165;
    const ccomplex_t IT_0215 = 0.5*IT_0190;
    const ccomplex_t IT_0216 = s_34 + IT_0215;
    const ccomplex_t IT_0217 = m_C_1*IT_0216;
    const ccomplex_t IT_0218 = 2*IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = (-0.25)*IT_0197;
    const ccomplex_t IT_0221 = s_14 + IT_0220;
    const ccomplex_t IT_0222 = m_tau*IT_0221;
    const ccomplex_t IT_0223 = (-4)*IT_0222;
    const ccomplex_t IT_0224 = (-2)*IT_0223;
    const ccomplex_t IT_0225 = (conj(IT_0169) + conj(IT_0176))*(IT_0104
      *IT_0128 + IT_0121*IT_0181 + IT_0186 + IT_0187 + IT_0091*IT_0195 + IT_0044
      *IT_0202 + IT_0176*IT_0210 + IT_0211 + IT_0212 + IT_0213 + IT_0214 +
       IT_0068*IT_0219 + IT_0019*IT_0224);
    const ccomplex_t IT_0226 = IT_0092*IT_0097;
    const ccomplex_t IT_0227 = IT_0057*IT_0102;
    const ccomplex_t IT_0228 = IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = IT_0148*IT_0188;
    const ccomplex_t IT_0230 = (-3)*IT_0229;
    const ccomplex_t IT_0231 = 2*IT_0230;
    const ccomplex_t IT_0232 = pow(m_Z, 4);
    const ccomplex_t IT_0233 = s_13*IT_0232;
    const ccomplex_t IT_0234 = IT_0123*IT_0233;
    const ccomplex_t IT_0235 = (-0.25)*IT_0234;
    const ccomplex_t IT_0236 = IT_0204 + IT_0235;
    const ccomplex_t IT_0237 = 8*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0234;
    const ccomplex_t IT_0239 = IT_0203 + IT_0205 + IT_0238;
    const ccomplex_t IT_0240 = (-4)*IT_0239;
    const ccomplex_t IT_0241 = 3*IT_0229;
    const ccomplex_t IT_0242 = (-2)*IT_0241;
    const ccomplex_t IT_0243 = IT_0052*IT_0057 + IT_0097*IT_0130;
    const ccomplex_t IT_0244 = IT_0105*IT_0122;
    const ccomplex_t IT_0245 = IT_0122*IT_0128;
    const ccomplex_t IT_0246 = IT_0068*IT_0224;
    const ccomplex_t IT_0247 = IT_0130*IT_0131;
    const ccomplex_t IT_0248 = IT_0052*IT_0138;
    const ccomplex_t IT_0249 = IT_0092*IT_0156;
    const ccomplex_t IT_0250 = IT_0102*IT_0164;
    const ccomplex_t IT_0251 = IT_0121*IT_0210;
    const ccomplex_t IT_0252 = IT_0091*IT_0202;
    const ccomplex_t IT_0253 = IT_0045*IT_0058 + (IT_0068 + IT_0091)*(IT_0098 
      + IT_0103 + IT_0106) + IT_0121*(IT_0129 + IT_0132 + IT_0139) + conj
      (IT_0130)*(IT_0045*IT_0097 + IT_0130*IT_0145 + IT_0052*IT_0147 + IT_0092
      *IT_0150 + IT_0102*IT_0152 + IT_0122*IT_0153 + IT_0104*IT_0154) + conj
      (IT_0092)*(IT_0092*IT_0145 + IT_0102*IT_0147 + IT_0130*IT_0150 + IT_0052
      *IT_0152 + IT_0104*IT_0153 + IT_0122*IT_0154 + IT_0121*IT_0156) + conj
      (IT_0122)*(IT_0045*IT_0105 + IT_0122*IT_0145 + IT_0104*IT_0150 + IT_0130
      *IT_0153 + IT_0092*IT_0154 + IT_0052*IT_0161 + IT_0102*IT_0162) + conj
      (IT_0052)*(IT_0130*IT_0147 + IT_0092*IT_0152 + IT_0052*IT_0158 + IT_0102
      *IT_0160 + IT_0122*IT_0161 + IT_0104*IT_0162) + conj(IT_0102)*(IT_0092
      *IT_0147 + IT_0130*IT_0152 + IT_0102*IT_0158 + IT_0052*IT_0160 + IT_0104
      *IT_0161 + IT_0122*IT_0162 + IT_0121*IT_0164) + conj(IT_0104)*(IT_0104
      *IT_0145 + IT_0122*IT_0150 + IT_0092*IT_0153 + IT_0130*IT_0154 + IT_0102
      *IT_0161 + IT_0052*IT_0162 + IT_0121*IT_0165) + IT_0177*(IT_0178 + IT_0179
       + IT_0180 + IT_0182 + IT_0183 + IT_0184 + IT_0185) + IT_0225 + conj
      (IT_0091)*(IT_0104*IT_0105 + IT_0177*IT_0195 + IT_0121*IT_0202 + IT_0228 +
       IT_0044*IT_0231 + IT_0068*IT_0237 + IT_0091*IT_0240 + IT_0019*IT_0242) +
       conj(IT_0068)*(IT_0104*IT_0105 + IT_0177*IT_0219 + IT_0121*IT_0224 +
       IT_0228 + IT_0019*IT_0231 + IT_0091*IT_0237 + IT_0068*IT_0240 + IT_0044
      *IT_0242) + conj(IT_0044)*(IT_0121*IT_0195 + IT_0177*IT_0202 + IT_0091
      *IT_0231 + IT_0019*IT_0237 + IT_0044*IT_0240 + IT_0068*IT_0242 + IT_0243 +
       IT_0244) + conj(IT_0019)*(IT_0121*IT_0219 + IT_0177*IT_0224 + IT_0068
      *IT_0231 + IT_0044*IT_0237 + IT_0019*IT_0240 + IT_0091*IT_0242 + IT_0243 +
       IT_0244) + conj(IT_0111)*(IT_0104*IT_0165 + IT_0044*IT_0195 + IT_0019
      *IT_0219 + IT_0245 + IT_0246 + IT_0247 + IT_0248 + IT_0249 + IT_0250 +
       IT_0251 + IT_0252) + conj(IT_0120)*(IT_0104*IT_0165 + IT_0177*IT_0181 +
       IT_0044*IT_0195 + IT_0019*IT_0219 + IT_0245 + IT_0246 + IT_0247 + IT_0248
       + IT_0249 + IT_0250 + IT_0251 + IT_0252);
    return create_ccomplex_return(IT_0253);
}

