#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_st_2
      *Gamma_t2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (-2)*IT_0002;
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0009;
    const ccomplex_t IT_0012 = 5.33333333333333*s_23;
    const ccomplex_t IT_0013 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_st_2
      *Gamma_t2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0018 = s_24*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0007;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = pow(s_23, 2);
    const ccomplex_t IT_0023 = 5.33333333333333*IT_0022;
    const ccomplex_t IT_0024 = s_23*IT_0017;
    const ccomplex_t IT_0025 = IT_0014*IT_0019;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = s_34*IT_0017;
    const ccomplex_t IT_0028 = -IT_0014*(IT_0008 + -IT_0019);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = 5.33333333333333*s_13;
    const ccomplex_t IT_0032 = s_24*IT_0031;
    const ccomplex_t IT_0033 = pow(m_st_2, 2);
    const ccomplex_t IT_0034 = 5.33333333333333*s_24;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0037 = s_24*IT_0036;
    const ccomplex_t IT_0038 = conj(IT_0010)*IT_0037;
    const ccomplex_t IT_0039 = s_34*IT_0036;
    const ccomplex_t IT_0040 = conj(IT_0021)*IT_0027;
    const ccomplex_t IT_0041 = 0.5*IT_0020;
    const ccomplex_t IT_0042 = cpow(s_12 + IT_0033 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0043 = -g_s;
    const ccomplex_t IT_0044 = IT_0007*IT_0043;
    const ccomplex_t IT_0045 = IT_0003*IT_0043;
    const ccomplex_t IT_0046 = IT_0042*(IT_0044 + -IT_0045);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0046;
    const ccomplex_t IT_0049 = (-0.125)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = IT_0012*IT_0033;
    const ccomplex_t IT_0052 = IT_0035 + IT_0051;
    const ccomplex_t IT_0053 = g_s*IT_0003;
    const ccomplex_t IT_0054 = IT_0042*(IT_0045 + IT_0053);
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = 0.666666666666667*s_12;
    const ccomplex_t IT_0057 = s_13*IT_0056;
    const ccomplex_t IT_0058 = s_14*IT_0056;
    const ccomplex_t IT_0059 = g_s*IT_0007;
    const ccomplex_t IT_0060 = IT_0042*(IT_0045 + IT_0053 + -IT_0059);
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = s_24*IT_0056;
    const ccomplex_t IT_0063 = 0.666666666666667*s_13;
    const ccomplex_t IT_0064 = IT_0033*IT_0063;
    const ccomplex_t IT_0065 = IT_0042*IT_0045;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = s_34*IT_0056;
    const ccomplex_t IT_0068 = 0.666666666666667*s_14;
    const ccomplex_t IT_0069 = s_24*IT_0068;
    const ccomplex_t IT_0070 = IT_0042*(IT_0044 + -IT_0045 + -IT_0053);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = s_23*IT_0056;
    const ccomplex_t IT_0074 = IT_0033*IT_0068;
    const ccomplex_t IT_0075 = pow(g_s, 2);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = s_14 + s_24;
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = g_s*IT_0078;
    const ccomplex_t IT_0080 = s_13 + s_23;
    const ccomplex_t IT_0081 = IT_0045*IT_0080;
    const ccomplex_t IT_0082 = IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = s_23*IT_0044;
    const ccomplex_t IT_0084 = g_s*s_24;
    const ccomplex_t IT_0085 = IT_0007*IT_0084;
    const ccomplex_t IT_0086 = -IT_0083 + -IT_0085;
    const ccomplex_t IT_0087 = IT_0082 + IT_0086;
    const ccomplex_t IT_0088 = IT_0042*IT_0087;
    const ccomplex_t IT_0089 = IT_0076 + (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = 0.666666666666667*s_23;
    const ccomplex_t IT_0091 = 0.666666666666667*s_24;
    const ccomplex_t IT_0092 = IT_0042*(IT_0044 + -IT_0045 + -IT_0053 +
       IT_0059);
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0033*IT_0090;
    const ccomplex_t IT_0095 = IT_0033*IT_0091;
    const ccomplex_t IT_0096 = 0.666666666666667*s_34;
    const ccomplex_t IT_0097 = IT_0033*IT_0096;
    const ccomplex_t IT_0098 = IT_0042*(IT_0053 + -IT_0059);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = 0.666666666666667*IT_0022;
    const ccomplex_t IT_0101 = 5.33333333333333*s_34;
    const ccomplex_t IT_0102 = IT_0033*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0054;
    const ccomplex_t IT_0104 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0105 = s_13*IT_0104;
    const ccomplex_t IT_0106 = s_14*IT_0104;
    const ccomplex_t IT_0107 = IT_0042*IT_0053;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = s_34*IT_0104;
    const ccomplex_t IT_0110 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0111 = s_23*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0060;
    const ccomplex_t IT_0113 = s_24*IT_0104;
    const ccomplex_t IT_0114 = IT_0033*IT_0110;
    const ccomplex_t IT_0115 = 0.5*IT_0065;
    const ccomplex_t IT_0116 = (-5.33333333333333)*s_14;
    const ccomplex_t IT_0117 = s_24*IT_0116;
    const ccomplex_t IT_0118 = -IT_0046;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = s_23*IT_0104;
    const ccomplex_t IT_0121 = IT_0033*IT_0116;
    const ccomplex_t IT_0122 = IT_0076 + 0.5*IT_0088;
    const ccomplex_t IT_0123 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0124 = (-5.33333333333333)*s_24;
    const ccomplex_t IT_0125 = 0.5*IT_0092;
    const ccomplex_t IT_0126 = IT_0033*IT_0123;
    const ccomplex_t IT_0127 = IT_0033*IT_0124;
    const ccomplex_t IT_0128 = (-5.33333333333333)*s_34;
    const ccomplex_t IT_0129 = IT_0033*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0098;
    const ccomplex_t IT_0131 = (-5.33333333333333)*IT_0022;
    const ccomplex_t IT_0132 = IT_0017*IT_0033;
    const ccomplex_t IT_0133 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0134 = s_24*IT_0133;
    const ccomplex_t IT_0135 = s_34*IT_0012*conj(IT_0021) + IT_0018*conj
      (IT_0026) + s_24*(10.6666666666667*s_23*conj(IT_0041) + (-5.33333333333333
      )*s_24*IT_0050) + conj(IT_0010)*IT_0052 + IT_0055*(IT_0057 + IT_0058) +
       IT_0061*(IT_0062 + IT_0064) + IT_0066*(IT_0067 + IT_0069) + IT_0072*
      (IT_0073 + IT_0074) + IT_0089*(IT_0090 + IT_0091) + IT_0093*(IT_0094 +
       IT_0095) + IT_0048*IT_0097 + IT_0099*(IT_0097 + IT_0100) + conj(IT_0011)*
      (IT_0023 + IT_0102) + IT_0103*(IT_0105 + IT_0106) + IT_0108*(IT_0109 +
       IT_0111) + IT_0112*(IT_0113 + IT_0114) + IT_0115*(IT_0109 + IT_0117) +
       IT_0119*(IT_0120 + IT_0121) + IT_0122*(IT_0123 + IT_0124) + IT_0125*
      (IT_0126 + IT_0127) + IT_0047*IT_0129 + IT_0130*(IT_0129 + IT_0131) + conj
      (IT_0016)*(IT_0037 + IT_0132) + conj(IT_0030)*IT_0134;
    const ccomplex_t IT_0136 = -IT_0014*(IT_0004 + -IT_0019);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = s_34*IT_0133;
    const ccomplex_t IT_0140 = s_23*IT_0063;
    const ccomplex_t IT_0141 = s_34*IT_0063;
    const ccomplex_t IT_0142 = 5.33333333333333*s_12;
    const ccomplex_t IT_0143 = s_23*IT_0142;
    const ccomplex_t IT_0144 = s_34*IT_0031;
    const ccomplex_t IT_0145 = s_34*IT_0110;
    const ccomplex_t IT_0146 = IT_0033*IT_0133;
    const ccomplex_t IT_0147 = (-0.666666666666667)*s_34;
    const ccomplex_t IT_0148 = IT_0033*IT_0147;
    const ccomplex_t IT_0149 = IT_0024*conj(IT_0030) + IT_0026*IT_0097 +
       IT_0010*IT_0120 + IT_0027*conj(IT_0138) + IT_0030*IT_0140 + IT_0138
      *IT_0141 + conj(IT_0010)*IT_0143 + conj(IT_0006)*IT_0144 + IT_0006*IT_0145
       + conj(IT_0016)*IT_0146 + conj(IT_0026)*IT_0148;
    const ccomplex_t IT_0150 = pow(s_24, 2);
    const ccomplex_t IT_0151 = 5.33333333333333*s_14;
    const ccomplex_t IT_0152 = IT_0033*IT_0151;
    const ccomplex_t IT_0153 = IT_0143 + IT_0152;
    const ccomplex_t IT_0154 = IT_0010 + IT_0125;
    const ccomplex_t IT_0155 = s_13*IT_0142;
    const ccomplex_t IT_0156 = s_14*IT_0142;
    const ccomplex_t IT_0157 = s_34*IT_0142;
    const ccomplex_t IT_0158 = s_23*IT_0031;
    const ccomplex_t IT_0159 = s_24*IT_0142;
    const ccomplex_t IT_0160 = IT_0031*IT_0033;
    const ccomplex_t IT_0161 = s_24*IT_0151;
    const ccomplex_t IT_0162 = s_13*IT_0036;
    const ccomplex_t IT_0163 = s_14*IT_0036;
    const ccomplex_t IT_0164 = (-0.666666666666667)*s_14;
    const ccomplex_t IT_0165 = s_24*IT_0164;
    const ccomplex_t IT_0166 = s_23*IT_0036;
    const ccomplex_t IT_0167 = IT_0033*IT_0164;
    const ccomplex_t IT_0168 = (-0.666666666666667)*s_24;
    const ccomplex_t IT_0169 = IT_0033*IT_0168;
    const ccomplex_t IT_0170 = (-0.666666666666667)*IT_0022;
    const ccomplex_t IT_0171 = 0.1875*IT_0018*IT_0026 + 0.1875*IT_0034*IT_0122
       + 0.1875*IT_0012*(s_34*IT_0021 + IT_0122) + 0.1875*IT_0102*(IT_0011 +
       IT_0047 + IT_0130) + 0.1875*IT_0030*IT_0134 + 0.1875*IT_0048*IT_0148 +
       IT_0050*IT_0150 + 0.1875*IT_0119*IT_0153 + 0.1875*IT_0052*IT_0154 +
       0.1875*IT_0103*(IT_0155 + IT_0156) + 0.1875*IT_0108*(IT_0157 + IT_0158) +
       0.1875*IT_0112*(IT_0159 + IT_0160) + 0.1875*IT_0115*(IT_0157 + IT_0161) +
       0.1875*IT_0055*(IT_0162 + IT_0163) + 0.1875*IT_0066*(IT_0039 + IT_0165) +
       0.1875*IT_0072*(IT_0166 + IT_0167) + 0.1875*IT_0089*(IT_0133 + IT_0168) +
       0.1875*IT_0093*(IT_0146 + IT_0169) + 0.1875*IT_0099*(IT_0148 + IT_0170);
    const ccomplex_t IT_0172 = s_24*IT_0090;
    const ccomplex_t IT_0173 = pow(m_st_2, 4);
    const ccomplex_t IT_0174 = 0.666666666666667*IT_0173;
    const ccomplex_t IT_0175 = 1.33333333333333*s_12;
    const ccomplex_t IT_0176 = s_23*IT_0175;
    const ccomplex_t IT_0177 = s_24*IT_0175;
    const ccomplex_t IT_0178 = IT_0112 + IT_0119;
    const ccomplex_t IT_0179 = IT_0033*IT_0175;
    const ccomplex_t IT_0180 = 1.33333333333333*s_23;
    const ccomplex_t IT_0181 = IT_0033*IT_0180;
    const ccomplex_t IT_0182 = 1.33333333333333*s_24;
    const ccomplex_t IT_0183 = IT_0033*IT_0182;
    const ccomplex_t IT_0184 = 1.33333333333333*IT_0033;
    const ccomplex_t IT_0185 = IT_0033*IT_0142;
    const ccomplex_t IT_0186 = conj(IT_0016)*IT_0185;
    const ccomplex_t IT_0187 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0188 = s_24*IT_0187;
    const ccomplex_t IT_0189 = IT_0033*IT_0187;
    const ccomplex_t IT_0190 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0191 = IT_0033*IT_0190;
    const ccomplex_t IT_0192 = (-10.6666666666667)*s_24;
    const ccomplex_t IT_0193 = IT_0033*IT_0192;
    const ccomplex_t IT_0194 = (-10.6666666666667)*IT_0033;
    const ccomplex_t IT_0195 = IT_0089*IT_0194;
    const ccomplex_t IT_0196 = IT_0033*IT_0104;
    const ccomplex_t IT_0197 = IT_0016*IT_0196;
    const ccomplex_t IT_0198 = s_24*IT_0123;
    const ccomplex_t IT_0199 = (-5.33333333333333)*IT_0173;
    const ccomplex_t IT_0200 = (-0.666666666666667)*IT_0173;
    const ccomplex_t IT_0201 = conj(IT_0030)*IT_0051 + IT_0011*IT_0094 +
       IT_0021*IT_0100 + IT_0026*IT_0113 + IT_0030*IT_0126 + conj(IT_0006)
      *IT_0134 + s_24*IT_0012*conj(IT_0138) + conj(IT_0011)*IT_0146 + conj
      (IT_0026)*IT_0159 + conj(IT_0021)*IT_0170 + IT_0006*IT_0172 + IT_0010
      *IT_0174 + IT_0108*IT_0176 + IT_0115*IT_0177 + IT_0178*IT_0179 + IT_0130
      *IT_0181 + IT_0047*IT_0183 + IT_0122*IT_0184 + IT_0186 + IT_0066*IT_0188 +
       IT_0072*IT_0189 + IT_0099*IT_0191 + IT_0048*IT_0193 + IT_0195 + IT_0197 +
       IT_0138*IT_0198 + IT_0093*IT_0199 + conj(IT_0010)*IT_0200;
    const ccomplex_t IT_0202 = s_24*IT_0063;
    const ccomplex_t IT_0203 = s_14*IT_0175;
    const ccomplex_t IT_0204 = IT_0089*IT_0175;
    const ccomplex_t IT_0205 = 1.33333333333333*s_14;
    const ccomplex_t IT_0206 = IT_0033*IT_0205;
    const ccomplex_t IT_0207 = conj(IT_0010)*IT_0185;
    const ccomplex_t IT_0208 = conj(IT_0021)*IT_0158;
    const ccomplex_t IT_0209 = s_14*IT_0187;
    const ccomplex_t IT_0210 = s_23*IT_0187;
    const ccomplex_t IT_0211 = IT_0122*IT_0187;
    const ccomplex_t IT_0212 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0213 = IT_0033*IT_0212;
    const ccomplex_t IT_0214 = IT_0010*IT_0196;
    const ccomplex_t IT_0215 = IT_0021*IT_0111;
    const ccomplex_t IT_0216 = s_24*IT_0110;
    const ccomplex_t IT_0217 = conj(IT_0006)*IT_0032 + IT_0030*IT_0064 +
       IT_0026*IT_0095 + IT_0011*IT_0120 + conj(IT_0030)*IT_0132 + IT_0018*conj
      (IT_0138) + conj(IT_0011)*IT_0143 + conj(IT_0026)*IT_0169 + IT_0016
      *IT_0174 + IT_0099*IT_0176 + IT_0108*IT_0191 + IT_0112*IT_0199 + conj
      (IT_0016)*IT_0200 + IT_0138*IT_0202 + IT_0066*IT_0203 + IT_0204 + IT_0048
      *IT_0206 + IT_0207 + IT_0208 + IT_0115*IT_0209 + IT_0130*IT_0210 + IT_0211
       + IT_0047*IT_0213 + IT_0214 + IT_0215 + IT_0006*IT_0216;
    const ccomplex_t IT_0218 = s_34*IT_0068;
    const ccomplex_t IT_0219 = 1.33333333333333*s_13;
    const ccomplex_t IT_0220 = s_24*IT_0219;
    const ccomplex_t IT_0221 = s_34*IT_0151;
    const ccomplex_t IT_0222 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0223 = s_24*IT_0222;
    const ccomplex_t IT_0224 = s_34*IT_0116;
    const ccomplex_t IT_0225 = s_34*IT_0164;
    const ccomplex_t IT_0226 = -conj(IT_0026);
    const ccomplex_t IT_0227 = IT_0026 + IT_0226;
    const ccomplex_t IT_0228 = 1.5*IT_0108;
    const ccomplex_t IT_0229 = conj(IT_0011)*IT_0032 + conj(IT_0030)*IT_0039 +
       IT_0016*IT_0058 + IT_0030*IT_0067 + IT_0010*IT_0113 + conj(IT_0021)
      *IT_0144 + IT_0021*IT_0145 + conj(IT_0010)*IT_0159 + conj(IT_0016)*IT_0163
       + IT_0011*IT_0216 + IT_0138*IT_0218 + IT_0099*IT_0220 + conj(IT_0006)
      *IT_0221 + IT_0130*IT_0223 + IT_0006*IT_0224 + conj(IT_0138)*IT_0225 +
       0.666666666666667*s_14*(s_14*IT_0227 + IT_0222*IT_0228);
    const ccomplex_t IT_0230 = (-0.1875)*IT_0010;
    const ccomplex_t IT_0231 = (-0.1875)*IT_0011;
    const ccomplex_t IT_0232 = (-0.1875)*conj(IT_0011);
    const ccomplex_t IT_0233 = (-0.1875)*IT_0006;
    const ccomplex_t IT_0234 = (-0.1875)*IT_0030;
    const ccomplex_t IT_0235 = (-0.1875)*conj(IT_0030);
    const ccomplex_t IT_0236 = (-0.1875)*IT_0138;
    const ccomplex_t IT_0237 = (-0.1875)*IT_0108;
    const ccomplex_t IT_0238 = (-0.1875)*IT_0130;
    const ccomplex_t IT_0239 = (-0.1875)*IT_0099;
    const ccomplex_t IT_0240 = IT_0038 + IT_0040 + IT_0016*IT_0106 + IT_0021
      *IT_0141 + conj(IT_0016)*IT_0156 + conj(IT_0138)*IT_0221 + conj(IT_0006)
      *IT_0225 + (-5.33333333333333)*IT_0062*IT_0230 + (-5.33333333333333)
      *IT_0202*IT_0231 + (-5.33333333333333)*IT_0018*IT_0232 + (
      -5.33333333333333)*IT_0218*IT_0233 + (-5.33333333333333)*IT_0109*IT_0234 +
       (-5.33333333333333)*IT_0157*IT_0235 + (-5.33333333333333)*IT_0224*IT_0236
       + (-5.33333333333333)*s_14*(s_14*IT_0227 + IT_0219*IT_0237) + (
      -5.33333333333333)*IT_0220*IT_0238 + (-5.33333333333333)*IT_0223*IT_0239;
    const ccomplex_t IT_0241 = pow(s_12, 2);
    const ccomplex_t IT_0242 = 0.666666666666667*IT_0241;
    const ccomplex_t IT_0243 = s_14*IT_0063;
    const ccomplex_t IT_0244 = s_13*IT_0175;
    const ccomplex_t IT_0245 = 1.33333333333333*IT_0241;
    const ccomplex_t IT_0246 = s_13*IT_0187;
    const ccomplex_t IT_0247 = (-10.6666666666667)*IT_0241;
    const ccomplex_t IT_0248 = s_14*IT_0110;
    const ccomplex_t IT_0249 = (-0.666666666666667)*IT_0241;
    const ccomplex_t IT_0250 = IT_0033*IT_0219;
    const ccomplex_t IT_0251 = pow(s_13, 2);
    const ccomplex_t IT_0252 = (IT_0021 + -conj(IT_0021))*IT_0251;
    const ccomplex_t IT_0253 = IT_0011*IT_0057 + IT_0030*IT_0105 + IT_0026
      *IT_0121 + s_14*(conj(IT_0006)*IT_0017 + IT_0031*conj(IT_0138)) + conj
      (IT_0026)*IT_0152 + conj(IT_0030)*IT_0155 + conj(IT_0011)*IT_0162 +
       IT_0178*IT_0179 + IT_0122*IT_0184 + IT_0186 + IT_0072*IT_0189 + IT_0195 +
       IT_0197 + IT_0055*IT_0199 + IT_0047*IT_0203 + IT_0115*IT_0206 + IT_0048
      *IT_0209 + IT_0066*IT_0213 + IT_0010*IT_0242 + IT_0006*IT_0243 + IT_0130
      *IT_0244 + IT_0125*IT_0245 + IT_0099*IT_0246 + IT_0093*IT_0247 + IT_0138
      *IT_0248 + conj(IT_0010)*IT_0249 + 0.666666666666667*IT_0228*IT_0250 +
       0.666666666666667*IT_0252;
    const ccomplex_t IT_0254 = IT_0033*IT_0056;
    const ccomplex_t IT_0255 = IT_0016*IT_0254;
    const ccomplex_t IT_0256 = IT_0089*IT_0184;
    const ccomplex_t IT_0257 = 1.33333333333333*IT_0173;
    const ccomplex_t IT_0258 = 5.33333333333333*IT_0241;
    const ccomplex_t IT_0259 = IT_0033*IT_0036;
    const ccomplex_t IT_0260 = conj(IT_0016)*IT_0259;
    const ccomplex_t IT_0261 = (-5.33333333333333)*IT_0241;
    const ccomplex_t IT_0262 = IT_0033*IT_0222;
    const ccomplex_t IT_0263 = IT_0026*IT_0074 + s_14*(conj(IT_0006)*IT_0031 +
       IT_0017*conj(IT_0138)) + conj(IT_0026)*IT_0167 + IT_0072*IT_0179 +
       IT_0178*IT_0189 + IT_0122*IT_0194 + IT_0103*IT_0199 + IT_0048*IT_0203 +
       IT_0066*IT_0206 + IT_0047*IT_0209 + IT_0115*IT_0213 + (-5.33333333333333)
      *IT_0105*IT_0231 + (-5.33333333333333)*IT_0155*IT_0232 + (
      -5.33333333333333)*IT_0057*IT_0234 + (-5.33333333333333)*IT_0162*IT_0235 +
       (-5.33333333333333)*IT_0236*IT_0243 + (-5.33333333333333)*IT_0239*IT_0244
       + IT_0093*IT_0245 + (-5.33333333333333)*IT_0238*IT_0246 + IT_0125*IT_0247
       + (-5.33333333333333)*IT_0233*IT_0248 + (-5.33333333333333)*IT_0252 +
       IT_0255 + IT_0256 + IT_0055*IT_0257 + conj(IT_0010)*IT_0258 + IT_0260 + (
      -5.33333333333333)*IT_0230*IT_0261 + (-5.33333333333333)*IT_0237*IT_0262;
    const ccomplex_t IT_0264 = IT_0021*IT_0140;
    const ccomplex_t IT_0265 = s_23*IT_0068;
    const ccomplex_t IT_0266 = IT_0122*IT_0175;
    const ccomplex_t IT_0267 = s_23*IT_0151;
    const ccomplex_t IT_0268 = IT_0089*IT_0187;
    const ccomplex_t IT_0269 = s_23*IT_0116;
    const ccomplex_t IT_0270 = conj(IT_0010)*IT_0259;
    const ccomplex_t IT_0271 = conj(IT_0021)*IT_0024;
    const ccomplex_t IT_0272 = s_23*IT_0164;
    const ccomplex_t IT_0273 = IT_0010*IT_0254;
    const ccomplex_t IT_0274 = IT_0011*IT_0064 + IT_0026*IT_0106 + IT_0030
      *IT_0120 + conj(IT_0011)*IT_0132 + conj(IT_0030)*IT_0143 + conj(IT_0026)
      *IT_0156 + IT_0047*IT_0177 + IT_0115*IT_0183 + IT_0048*IT_0188 + IT_0066
      *IT_0193 + IT_0072*IT_0199 + IT_0108*IT_0244 + IT_0112*IT_0245 + IT_0130
      *IT_0250 + conj(IT_0016)*IT_0258 + IT_0016*IT_0261 + IT_0099*IT_0262 +
       IT_0264 + IT_0006*IT_0265 + IT_0266 + conj(IT_0138)*IT_0267 + IT_0268 +
       IT_0138*IT_0269 + IT_0270 + IT_0271 + conj(IT_0006)*IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = s_34*IT_0090;
    const ccomplex_t IT_0276 = s_34*IT_0123;
    const ccomplex_t IT_0277 = conj(IT_0010)*IT_0146;
    const ccomplex_t IT_0278 = IT_0023*conj(IT_0030) + IT_0010*IT_0094 +
       IT_0026*IT_0109 + IT_0030*IT_0131 + s_34*IT_0012*conj(IT_0138) + conj
      (IT_0006)*IT_0139 + conj(IT_0016)*IT_0143 + conj(IT_0026)*IT_0157 +
       IT_0006*IT_0275 + IT_0138*IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = IT_0037*conj(IT_0041);
    const ccomplex_t IT_0280 = conj(IT_0006)*IT_0018;
    const ccomplex_t IT_0281 = conj(IT_0041)*IT_0132;
    const ccomplex_t IT_0282 = 5.33333333333333*IT_0173;
    const ccomplex_t IT_0283 = conj(IT_0016)*IT_0282;
    const ccomplex_t IT_0284 = conj(IT_0030)*IT_0160;
    const ccomplex_t IT_0285 = IT_0032*conj(IT_0138);
    const ccomplex_t IT_0286 = conj(IT_0026)*IT_0035 + IT_0011*IT_0073 +
       IT_0030*IT_0114 + IT_0026*IT_0127 + conj(IT_0011)*IT_0166 + IT_0130
      *IT_0176 + (IT_0103 + IT_0125)*IT_0179 + IT_0108*IT_0181 + (IT_0055 +
       IT_0093)*IT_0189 + (IT_0016 + IT_0061)*IT_0199 + IT_0006*IT_0202 +
       IT_0115*IT_0203 + IT_0047*IT_0206 + IT_0066*IT_0209 + IT_0099*IT_0210 +
       IT_0048*IT_0213 + IT_0138*IT_0216 + IT_0119*IT_0245 + IT_0072*IT_0247 +
       IT_0112*IT_0257 + IT_0264 + IT_0266 + IT_0268 + IT_0270 + IT_0271 +
       IT_0273 + IT_0279 + IT_0280 + IT_0281 + IT_0283 + IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = s_24*conj(IT_0006)*IT_0012;
    const ccomplex_t IT_0288 = conj(IT_0021)*IT_0023;
    const ccomplex_t IT_0289 = conj(IT_0010)*IT_0282;
    const ccomplex_t IT_0290 = conj(IT_0011)*IT_0051;
    const ccomplex_t IT_0291 = conj(IT_0030)*IT_0146;
    const ccomplex_t IT_0292 = conj(IT_0026)*IT_0037;
    const ccomplex_t IT_0293 = IT_0134*conj(IT_0138);
    const ccomplex_t IT_0294 = IT_0026*IT_0062 + IT_0030*IT_0094 + IT_0011
      *IT_0126 + IT_0021*IT_0131 + IT_0138*IT_0172 + IT_0066*IT_0177 + IT_0072
      *IT_0179 + IT_0099*IT_0181 + IT_0048*IT_0183 + IT_0115*IT_0188 + IT_0178
      *IT_0189 + IT_0130*IT_0191 + IT_0047*IT_0193 + IT_0122*IT_0194 + IT_0006
      *IT_0198 + IT_0154*IT_0199 + IT_0108*IT_0210 + IT_0255 + IT_0256 + IT_0093
      *IT_0257 + IT_0260 + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291 +
       IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = IT_0023*conj(IT_0041);
    const ccomplex_t IT_0296 = s_34*conj(IT_0006)*IT_0012;
    const ccomplex_t IT_0297 = conj(IT_0010)*IT_0051;
    const ccomplex_t IT_0298 = conj(IT_0016)*IT_0166;
    const ccomplex_t IT_0299 = conj(IT_0030)*IT_0170;
    const ccomplex_t IT_0300 = conj(IT_0026)*IT_0039;
    const ccomplex_t IT_0301 = conj(IT_0138)*IT_0139;
    const ccomplex_t IT_0302 = IT_0026*IT_0058;
    const ccomplex_t IT_0303 = IT_0030*IT_0073;
    const ccomplex_t IT_0304 = IT_0016*IT_0242;
    const ccomplex_t IT_0305 = IT_0011*IT_0114;
    const ccomplex_t IT_0306 = conj(IT_0026)*IT_0163;
    const ccomplex_t IT_0307 = conj(IT_0030)*IT_0166;
    const ccomplex_t IT_0308 = conj(IT_0016)*IT_0249;
    const ccomplex_t IT_0309 = IT_0115*IT_0193;
    const ccomplex_t IT_0310 = IT_0108*IT_0246;
    const ccomplex_t IT_0311 = IT_0130*IT_0262;
    const ccomplex_t IT_0312 = conj(IT_0011)*IT_0160 + IT_0048*IT_0177 +
       IT_0066*IT_0183 + IT_0047*IT_0188 + IT_0119*IT_0199 + IT_0204 + IT_0207 +
       IT_0208 + IT_0211 + IT_0214 + IT_0215 + IT_0112*IT_0247 + IT_0099*IT_0250
       + IT_0072*IT_0257 + IT_0138*IT_0265 + conj(IT_0006)*IT_0267 + IT_0006
      *IT_0269 + conj(IT_0138)*IT_0272 + IT_0302 + IT_0303 + IT_0304 + IT_0305 +
       IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = (-0.5)*IT_0107;
    const ccomplex_t IT_0314 = conj(IT_0016)*IT_0051 + IT_0010*IT_0073 + conj
      (IT_0026)*IT_0102 + IT_0030*IT_0111 + IT_0016*IT_0126 + IT_0026*IT_0129 +
       IT_0041*(IT_0067 + IT_0120 + IT_0121 + IT_0140) + conj(IT_0041)*(IT_0024 
      + IT_0039 + IT_0153) + conj(IT_0030)*IT_0158 + conj(IT_0011)*IT_0160 +
       conj(IT_0010)*IT_0166 + IT_0125*(IT_0176 + IT_0189) + IT_0204 + IT_0047*
      (IT_0188 + s_23*IT_0205) + IT_0207 + IT_0208 + IT_0093*(IT_0179 + IT_0210)
       + IT_0211 + IT_0048*(IT_0177 + s_23*IT_0212) + IT_0214 + IT_0215 + (s_14
      *IT_0115 + IT_0122)*IT_0219 + IT_0089*IT_0222 + IT_0066*(IT_0183 + s_14
      *IT_0222) + IT_0119*((-10.6666666666667)*IT_0173 + IT_0244) + IT_0061*
      (IT_0191 + IT_0245) + IT_0112*(IT_0181 + IT_0247) + IT_0099*IT_0250 +
       IT_0103*(IT_0189 + IT_0250) + IT_0072*(IT_0246 + IT_0257) + IT_0055*
      (IT_0179 + IT_0262) + IT_0138*(IT_0145 + IT_0265) + conj(IT_0006)*(IT_0027
       + IT_0267) + IT_0006*(IT_0141 + IT_0269) + conj(IT_0138)*(IT_0144 +
       IT_0272) + IT_0302 + IT_0303 + IT_0304 + IT_0305 + IT_0306 + IT_0307 +
       IT_0308 + IT_0309 + IT_0310 + IT_0311 + (IT_0199 + IT_0244)*IT_0313;
    const ccomplex_t IT_0315 = conj(IT_0011)*IT_0012;
    const ccomplex_t IT_0316 = IT_0130*IT_0190;
    const ccomplex_t IT_0317 = IT_0099*IT_0180;
    const ccomplex_t IT_0318 = conj(IT_0010)*IT_0035 + IT_0026*IT_0069 +
       IT_0016*IT_0074 + IT_0030*IT_0097 + IT_0010*IT_0127 + conj(IT_0030)
      *IT_0148 + conj(IT_0026)*IT_0165 + conj(IT_0016)*IT_0167 + s_34*(IT_0012
      *conj(IT_0021) + conj(IT_0006)*IT_0034 + IT_0006*IT_0124 + IT_0091*IT_0138
       + conj(IT_0138)*IT_0168) + IT_0011*IT_0198 + s_23*IT_0108*IT_0212 +
       IT_0021*IT_0276 + s_24*(IT_0315 + IT_0316 + IT_0317);
    const ccomplex_t IT_0319 = m_st_2*(IT_0010 + -conj(IT_0010));
    const ccomplex_t IT_0320 = conj(IT_0016)*IT_0036 + IT_0016*IT_0056 +
       IT_0026*IT_0068 + 5.33333333333333*IT_0089 + IT_0030*IT_0090 + conj
      (IT_0006)*IT_0101 + (-21.3333333333333)*IT_0122 + IT_0011*IT_0123 +
       IT_0006*IT_0128 + conj(IT_0030)*IT_0133 + IT_0096*IT_0138 + conj(IT_0138)
      *IT_0147 + conj(IT_0026)*IT_0164 + IT_0048*IT_0182 + IT_0047*IT_0192 +
       IT_0066*IT_0205 + IT_0115*IT_0212 + IT_0108*IT_0222 + IT_0315 + IT_0316 +
       IT_0317 + (-5.33333333333333)*m_st_2*IT_0319;
    const ccomplex_t IT_0321 = IT_0130*IT_0180;
    const ccomplex_t IT_0322 = IT_0099*IT_0190;
    const ccomplex_t IT_0323 = IT_0012*conj(IT_0030) + (-21.3333333333333)
      *IT_0089 + IT_0011*IT_0090 + IT_0006*IT_0096 + IT_0016*IT_0104 + IT_0026
      *IT_0116 + IT_0030*IT_0123 + conj(IT_0011)*IT_0133 + IT_0128*IT_0138 +
       IT_0101*conj(IT_0138) + conj(IT_0016)*IT_0142 + conj(IT_0006)*IT_0147 +
       conj(IT_0026)*IT_0151 + IT_0047*IT_0182 + IT_0048*IT_0192 + IT_0115
      *IT_0205 + IT_0066*IT_0212 + IT_0108*IT_0219 + 0.666666666666667*m_st_2
      *IT_0319 + IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = IT_0010*IT_0095 + conj(IT_0030)*IT_0102 +
       IT_0026*IT_0117 + IT_0016*IT_0121 + IT_0030*IT_0129 + conj(IT_0011)
      *IT_0134 + conj(IT_0021)*IT_0139 + conj(IT_0016)*IT_0152 + conj(IT_0026)
      *IT_0161 + s_34*(IT_0006*IT_0091 + IT_0124*IT_0138 + IT_0034*conj(IT_0138)
       + conj(IT_0006)*IT_0168) + conj(IT_0010)*IT_0169 + IT_0011*IT_0172 + s_23
      *IT_0108*IT_0205 + IT_0021*IT_0275 + s_24*(IT_0321 + IT_0322);
    const ccomplex_t IT_0325 = IT_0006*((s_24*conj(IT_0010) + s_34*conj
      (IT_0011))*IT_0012 + conj(IT_0016)*IT_0018) + IT_0021*(conj(IT_0010)
      *IT_0023 + conj(IT_0016)*IT_0024 + conj(IT_0026)*IT_0027) + IT_0026*(conj
      (IT_0030)*IT_0032 + conj(IT_0016)*IT_0035 + IT_0038 + conj(IT_0011)
      *IT_0039 + IT_0040) + IT_0041*IT_0135 + IT_0138*(conj(IT_0016)*IT_0032 +
       conj(IT_0010)*IT_0134 + conj(IT_0011)*IT_0139) + IT_0108*IT_0149 +
       5.33333333333333*conj(IT_0041)*IT_0171 + IT_0093*IT_0201 + IT_0112
      *IT_0217 + IT_0115*IT_0229 + IT_0066*IT_0240 + IT_0055*IT_0253 + IT_0103
      *IT_0263 + IT_0072*IT_0274 + IT_0030*(conj(IT_0026)*IT_0032 + conj(IT_0016
      )*IT_0160 + conj(IT_0011)*IT_0170 + IT_0277) + IT_0099*IT_0278 + IT_0016*
      (conj(IT_0026)*IT_0035 + IT_0094*IT_0108 + IT_0099*IT_0120 + IT_0073
      *IT_0130 + conj(IT_0011)*IT_0166 + IT_0270 + IT_0271 + IT_0279 + IT_0280 +
       IT_0281 + IT_0283 + IT_0284 + IT_0285) + IT_0061*IT_0286 + IT_0010*(conj
      (IT_0016)*IT_0259 + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291 +
       IT_0292 + IT_0293) + IT_0125*IT_0294 + IT_0011*(IT_0295 + IT_0296 +
       IT_0297 + IT_0298 + IT_0299 + IT_0300 + IT_0301) + IT_0130*(IT_0026
      *IT_0067 + IT_0030*IT_0100 + IT_0010*IT_0126 + IT_0138*IT_0275 + IT_0006
      *IT_0276 + IT_0295 + IT_0296 + IT_0297 + IT_0298 + IT_0299 + IT_0300 +
       IT_0301) + IT_0119*IT_0312 + IT_0313*IT_0314 + IT_0047*IT_0318 + IT_0122
      *IT_0320 + IT_0089*IT_0323 + IT_0048*IT_0324;
    return create_ccomplex_return(IT_0325);
}

