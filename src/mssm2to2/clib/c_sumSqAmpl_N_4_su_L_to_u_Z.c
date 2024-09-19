#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_su_L_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_su_L_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -1);
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0006*IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_u, 2);
    const ccomplex_t IT_0023 = pow(m_su_L, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0016*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0015;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0032 = IT_0007*IT_0031;
    const ccomplex_t IT_0033 = IT_0003*IT_0031;
    const ccomplex_t IT_0034 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = IT_0003*IT_0034;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0032 + IT_0033 + 
      -IT_0035 + -IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0006*IT_0017*IT_0019;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = IT_0007*IT_0048;
    const ccomplex_t IT_0051 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0052 = IT_0007*IT_0051;
    const ccomplex_t IT_0053 = IT_0003*IT_0051;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0049 + IT_0050 + 
      -IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0006*IT_0017*IT_0019;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = IT_0007*IT_0065;
    const ccomplex_t IT_0068 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = IT_0007*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0006*IT_0017*IT_0019;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0022 + IT_0023 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0073*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0030 + -IT_0047 + -IT_0064 + -IT_0079;
    const ccomplex_t IT_0081 = conj(N_B4)*e_em;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W4)*e_em;
    const ccomplex_t IT_0085 = IT_0006*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + 3*IT_0086);
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = IT_0024*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0029*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0016*IT_0090;
    const ccomplex_t IT_0094 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = IT_0007*IT_0094;
    const ccomplex_t IT_0097 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0098 = IT_0007*IT_0097;
    const ccomplex_t IT_0099 = IT_0003*IT_0097;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0095 + IT_0096 + 
      -IT_0098 + -IT_0099);
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = conj(N_B2)*e_em;
    const ccomplex_t IT_0104 = IT_0001*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W2)*e_em;
    const ccomplex_t IT_0107 = IT_0006*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 3*IT_0108);
    const ccomplex_t IT_0110 = 0.166666666666667*IT_0109;
    const ccomplex_t IT_0111 = IT_0061*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0102*IT_0112;
    const ccomplex_t IT_0114 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0115 = IT_0003*IT_0114;
    const ccomplex_t IT_0116 = IT_0007*IT_0114;
    const ccomplex_t IT_0117 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0118 = IT_0003*IT_0117;
    const ccomplex_t IT_0119 = IT_0007*IT_0117;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0115 + IT_0116 + 
      -IT_0118 + -IT_0119);
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B3)*e_em;
    const ccomplex_t IT_0123 = IT_0001*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W3)*e_em;
    const ccomplex_t IT_0126 = IT_0006*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + 3*IT_0127);
    const ccomplex_t IT_0129 = 0.166666666666667*IT_0128;
    const ccomplex_t IT_0130 = IT_0076*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0121*IT_0131;
    const ccomplex_t IT_0133 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0134 = IT_0003*IT_0133;
    const ccomplex_t IT_0135 = IT_0007*IT_0133;
    const ccomplex_t IT_0136 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0137 = IT_0003*IT_0136;
    const ccomplex_t IT_0138 = IT_0007*IT_0136;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0134 + IT_0135 + 
      -IT_0137 + -IT_0138);
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B1)*e_em;
    const ccomplex_t IT_0144 = IT_0001*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W1)*e_em;
    const ccomplex_t IT_0147 = IT_0006*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + 3*IT_0148);
    const ccomplex_t IT_0150 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0151 = IT_0044*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0142*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0093 + IT_0113 + IT_0132 + IT_0153;
    const ccomplex_t IT_0155 = pow(m_Z, -2);
    const ccomplex_t IT_0156 = pow(m_Z, 4);
    const ccomplex_t IT_0157 = s_13*IT_0156;
    const ccomplex_t IT_0158 = IT_0155*IT_0157;
    const ccomplex_t IT_0159 = (-6)*IT_0158;
    const ccomplex_t IT_0160 = s_14*s_34;
    const ccomplex_t IT_0161 = (-24)*IT_0160;
    const ccomplex_t IT_0162 = pow(m_Z, 2);
    const ccomplex_t IT_0163 = s_13*IT_0162;
    const ccomplex_t IT_0164 = 12*IT_0163;
    const ccomplex_t IT_0165 = IT_0159 + IT_0161 + IT_0164;
    const ccomplex_t IT_0166 = 6*IT_0158;
    const ccomplex_t IT_0167 = 24*IT_0160;
    const ccomplex_t IT_0168 = (-12)*IT_0163;
    const ccomplex_t IT_0169 = IT_0166 + IT_0167 + IT_0168;
    const ccomplex_t IT_0170 = (-2)*IT_0160;
    const ccomplex_t IT_0171 = IT_0163 + IT_0170;
    const ccomplex_t IT_0172 = IT_0155*IT_0171;
    const ccomplex_t IT_0173 = (-6)*IT_0172;
    const ccomplex_t IT_0174 = 12*s_13;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = pow(m_N_4, 2);
    const ccomplex_t IT_0177 = cpow((-2)*s_12 + IT_0022 + -IT_0023 + -IT_0176 
      + -reg_prop, -1);
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0179 = e_em*IT_0003;
    const ccomplex_t IT_0180 = e_em*IT_0007;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*(IT_0179 + (-3)*IT_0180);
    const ccomplex_t IT_0182 = (-0.166666666666667)*IT_0181;
    const ccomplex_t IT_0183 = IT_0021*IT_0182;
    const ccomplex_t IT_0184 = IT_0177*IT_0178*IT_0183;
    const ccomplex_t IT_0185 = m_N_4*IT_0016;
    const ccomplex_t IT_0186 = m_N_4*IT_0029;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = IT_0090*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0191 = (-0.666666666666667)*IT_0190;
    const ccomplex_t IT_0192 = IT_0021*IT_0191;
    const ccomplex_t IT_0193 = IT_0177*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = m_u*IT_0194;
    const ccomplex_t IT_0196 = -IT_0184 + (-0.5)*IT_0189 + IT_0195;
    const ccomplex_t IT_0197 = m_N_1*IT_0041;
    const ccomplex_t IT_0198 = m_N_4*IT_0142;
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = IT_0152*IT_0200;
    const ccomplex_t IT_0202 = m_N_2*IT_0058;
    const ccomplex_t IT_0203 = m_N_4*IT_0102;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0112*IT_0205;
    const ccomplex_t IT_0207 = m_N_3*IT_0073;
    const ccomplex_t IT_0208 = m_N_4*IT_0121;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = IT_0131*IT_0210;
    const ccomplex_t IT_0212 = (-0.5)*IT_0189 + IT_0201 + IT_0206 + IT_0211;
    const ccomplex_t IT_0213 = IT_0196 + IT_0212;
    const ccomplex_t IT_0214 = conj(IT_0196) + conj(IT_0212);
    const ccomplex_t IT_0215 = IT_0088*IT_0191;
    const ccomplex_t IT_0216 = IT_0177*IT_0178*IT_0215;
    const ccomplex_t IT_0217 = IT_0088*IT_0182;
    const ccomplex_t IT_0218 = IT_0177*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_u*IT_0219;
    const ccomplex_t IT_0221 = IT_0026*IT_0188;
    const ccomplex_t IT_0222 = IT_0216 + -IT_0220 + (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = m_N_1*IT_0142;
    const ccomplex_t IT_0224 = m_N_4*IT_0041;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = IT_0046*IT_0226;
    const ccomplex_t IT_0228 = m_N_2*IT_0102;
    const ccomplex_t IT_0229 = m_N_4*IT_0058;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0063*IT_0231;
    const ccomplex_t IT_0233 = m_N_3*IT_0121;
    const ccomplex_t IT_0234 = m_N_4*IT_0073;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = IT_0233 + IT_0235;
    const ccomplex_t IT_0237 = IT_0078*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0221 + -IT_0227 + -IT_0232 + -IT_0237;
    const ccomplex_t IT_0239 = IT_0222 + IT_0238;
    const ccomplex_t IT_0240 = conj(IT_0222) + conj(IT_0238);
    const ccomplex_t IT_0241 = m_u*m_N_4;
    const ccomplex_t IT_0242 = pow(s_14, 2);
    const ccomplex_t IT_0243 = IT_0155*IT_0242;
    const ccomplex_t IT_0244 = -IT_0176;
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0241*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + (-3)
      *IT_0000*IT_0006);
    const ccomplex_t IT_0248 = 0.333333333333333*IT_0247;
    const ccomplex_t IT_0249 = cpow((-2)*s_13 + IT_0022 + IT_0176 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0250 = IT_0088*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0248*IT_0251;
    const ccomplex_t IT_0253 = -IT_0091 + -IT_0252;
    const ccomplex_t IT_0254 = IT_0021*IT_0249;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0248*IT_0255;
    const ccomplex_t IT_0257 = IT_0027 + IT_0256;
    const ccomplex_t IT_0258 = 6*IT_0257;
    const ccomplex_t IT_0259 = IT_0030 + 2*IT_0047 + 2*IT_0064 + 2*IT_0079;
    const ccomplex_t IT_0260 = (-6)*IT_0259;
    const ccomplex_t IT_0261 = IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = 6*conj(IT_0257);
    const ccomplex_t IT_0263 = (-6)*conj(IT_0259);
    const ccomplex_t IT_0264 = IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = -IT_0093 + (-2)*IT_0113 + (-2)*IT_0132 + (-2)
      *IT_0153;
    const ccomplex_t IT_0266 = 6*IT_0259;
    const ccomplex_t IT_0267 = (-6)*IT_0257;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = 6*conj(IT_0259);
    const ccomplex_t IT_0270 = (-6)*conj(IT_0257);
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = m_u*IT_0245;
    const ccomplex_t IT_0273 = 6*IT_0253;
    const ccomplex_t IT_0274 = (-6)*IT_0265;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 6*conj(IT_0253);
    const ccomplex_t IT_0277 = (-6)*conj(IT_0265);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = s_14*s_34*IT_0155;
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = s_13 + IT_0280;
    const ccomplex_t IT_0282 = s_13*IT_0281;
    const ccomplex_t IT_0283 = (-2)*IT_0194 + -IT_0256;
    const ccomplex_t IT_0284 = 2*IT_0219 + IT_0252;
    const ccomplex_t IT_0285 = 6*IT_0265;
    const ccomplex_t IT_0286 = (-6)*IT_0253;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = 6*conj(IT_0265);
    const ccomplex_t IT_0289 = (-6)*conj(IT_0253);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = s_13*s_14;
    const ccomplex_t IT_0292 = s_34*IT_0176;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = IT_0291 + IT_0293;
    const ccomplex_t IT_0295 = -IT_0194;
    const ccomplex_t IT_0296 = conj(IT_0080) + conj(IT_0295);
    const ccomplex_t IT_0297 = IT_0080 + IT_0295;
    const ccomplex_t IT_0298 = conj(IT_0154) + conj(IT_0219);
    const ccomplex_t IT_0299 = IT_0154 + IT_0219;
    const ccomplex_t IT_0300 = s_34*IT_0162;
    const ccomplex_t IT_0301 = IT_0155*IT_0300;
    const ccomplex_t IT_0302 = (-0.25)*IT_0301;
    const ccomplex_t IT_0303 = s_34 + IT_0302;
    const ccomplex_t IT_0304 = m_N_4*IT_0303;
    const ccomplex_t IT_0305 = (-4)*IT_0304;
    const ccomplex_t IT_0306 = 6*IT_0196;
    const ccomplex_t IT_0307 = 6*IT_0212;
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = 6*conj(IT_0196);
    const ccomplex_t IT_0310 = 6*conj(IT_0212);
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = 6*IT_0222;
    const ccomplex_t IT_0313 = 6*IT_0238;
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = 6*conj(IT_0222);
    const ccomplex_t IT_0316 = 6*conj(IT_0238);
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = (-6)*IT_0196;
    const ccomplex_t IT_0319 = (-6)*IT_0212;
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = (-6)*conj(IT_0196);
    const ccomplex_t IT_0322 = (-6)*conj(IT_0212);
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = (-6)*IT_0222;
    const ccomplex_t IT_0325 = (-6)*IT_0238;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0328 = (-6)*conj(IT_0238);
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = s_14*IT_0162;
    const ccomplex_t IT_0331 = IT_0155*IT_0330;
    const ccomplex_t IT_0332 = 0.5*IT_0331;
    const ccomplex_t IT_0333 = s_14 + IT_0332;
    const ccomplex_t IT_0334 = m_u*IT_0333;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = IT_0162*IT_0241;
    const ccomplex_t IT_0337 = 3*IT_0336;
    const ccomplex_t IT_0338 = 6*IT_0295;
    const ccomplex_t IT_0339 = 6*conj(IT_0295);
    const ccomplex_t IT_0340 = (-18)*conj(IT_0196);
    const ccomplex_t IT_0341 = IT_0241*IT_0281;
    const ccomplex_t IT_0342 = 6*conj(IT_0283);
    const ccomplex_t IT_0343 = (-6)*IT_0283;
    const ccomplex_t IT_0344 = (-6)*conj(IT_0283);
    const ccomplex_t IT_0345 = pow(s_34, 2);
    const ccomplex_t IT_0346 = IT_0155*IT_0345;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = IT_0022 + IT_0347;
    const ccomplex_t IT_0349 = IT_0241*IT_0348;
    const ccomplex_t IT_0350 = s_13*IT_0245;
    const ccomplex_t IT_0351 = s_13*IT_0348;
    const ccomplex_t IT_0352 = (-6)*IT_0284;
    const ccomplex_t IT_0353 = m_u*IT_0281;
    const ccomplex_t IT_0354 = (-6)*conj(IT_0284);
    const ccomplex_t IT_0355 = s_14*IT_0022;
    const ccomplex_t IT_0356 = s_13*s_34;
    const ccomplex_t IT_0357 = -IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = conj(IT_0028)*IT_0283 + conj(IT_0092)*IT_0284 +
       IT_0092*conj(IT_0284) + 0.166666666666667*IT_0028*IT_0342 +
       0.166666666666667*IT_0296*IT_0343 + 0.166666666666667*IT_0297*IT_0344 +
       0.166666666666667*IT_0298*IT_0352 + 0.166666666666667*IT_0299*IT_0354;
    const ccomplex_t IT_0360 = m_N_4*IT_0281;
    const ccomplex_t IT_0361 = m_N_4*IT_0348;
    const ccomplex_t IT_0362 = (-3)*IT_0336;
    const ccomplex_t IT_0363 = 6*IT_0028;
    const ccomplex_t IT_0364 = (-0.166666666666667)*IT_0363;
    const ccomplex_t IT_0365 = IT_0080 + IT_0364;
    const ccomplex_t IT_0366 = 6*conj(IT_0028);
    const ccomplex_t IT_0367 = (-0.166666666666667)*IT_0366;
    const ccomplex_t IT_0368 = conj(IT_0080) + IT_0367;
    const ccomplex_t IT_0369 = (-0.166666666666667)*conj(IT_0219);
    const ccomplex_t IT_0370 = IT_0159 + IT_0167;
    const ccomplex_t IT_0371 = IT_0161 + IT_0166;
    const ccomplex_t IT_0372 = (-0.25)*IT_0331;
    const ccomplex_t IT_0373 = s_14 + IT_0372;
    const ccomplex_t IT_0374 = m_u*IT_0373;
    const ccomplex_t IT_0375 = (-4)*IT_0374;
    const ccomplex_t IT_0376 = (-0.166666666666667)*IT_0375;
    const ccomplex_t IT_0377 = 0.5*IT_0301;
    const ccomplex_t IT_0378 = s_34 + IT_0377;
    const ccomplex_t IT_0379 = m_N_4*IT_0378;
    const ccomplex_t IT_0380 = 2*IT_0379;
    const ccomplex_t IT_0381 = (-0.166666666666667)*IT_0380;
    const ccomplex_t IT_0382 = (conj(IT_0028)*IT_0080 + IT_0028*conj(IT_0080) 
      + conj(IT_0092)*IT_0154 + IT_0092*conj(IT_0154))*IT_0165 + (IT_0028*conj
      (IT_0028) + IT_0080*conj(IT_0080) + IT_0092*conj(IT_0092) + IT_0154*conj
      (IT_0154))*IT_0169 + IT_0175*(IT_0213*IT_0214 + IT_0239*IT_0240) + IT_0246
      *(conj(IT_0253)*IT_0261 + IT_0253*IT_0264 + conj(IT_0265)*IT_0268 +
       IT_0265*IT_0271) + IT_0272*(IT_0214*IT_0261 + IT_0213*IT_0264 + IT_0240
      *IT_0275 + IT_0239*IT_0278) + IT_0282*(IT_0271*IT_0283 + IT_0268*conj
      (IT_0283) + conj(IT_0284)*IT_0287 + IT_0284*IT_0290) + IT_0294*(conj
      (IT_0028)*IT_0261 + IT_0028*IT_0264 + conj(IT_0092)*IT_0275 + IT_0092
      *IT_0278 + IT_0268*IT_0296 + IT_0271*IT_0297 + IT_0287*IT_0298 + IT_0290
      *IT_0299) + IT_0305*(conj(IT_0092)*IT_0308 + IT_0092*IT_0311 + conj
      (IT_0028)*IT_0314 + IT_0028*IT_0317 + conj(IT_0154)*IT_0320 + IT_0154
      *IT_0323 + conj(IT_0080)*IT_0326 + IT_0080*IT_0329) + (conj(IT_0028)
      *IT_0308 + IT_0028*IT_0311 + conj(IT_0092)*IT_0314 + IT_0092*IT_0317 +
       conj(IT_0080)*IT_0320 + IT_0080*IT_0323 + conj(IT_0154)*IT_0326 + IT_0154
      *IT_0329)*IT_0335 + IT_0337*(conj(IT_0092)*IT_0338 + IT_0092*IT_0339) +
       IT_0241*(IT_0213*((-18)*conj(IT_0222) + (-18)*conj(IT_0238)) + IT_0239*((
      -18)*conj(IT_0212) + IT_0340)) + 6*IT_0341*(conj(IT_0265)*IT_0283 +
       0.166666666666667*IT_0271*IT_0284 + 0.166666666666667*IT_0268*conj
      (IT_0284) + 0.166666666666667*IT_0265*IT_0342 + 0.166666666666667*conj
      (IT_0253)*IT_0343 + 0.166666666666667*IT_0253*IT_0344) + (conj(IT_0284)
      *IT_0343 + IT_0284*IT_0344)*IT_0349 + (conj(IT_0257)*IT_0258 + conj
      (IT_0259)*IT_0268 + IT_0259*IT_0270 + conj(IT_0253)*IT_0273 + conj(IT_0265
      )*IT_0287 + IT_0265*IT_0289)*IT_0350 + IT_0351*(conj(IT_0283)*IT_0343 +
       conj(IT_0284)*IT_0352) + IT_0353*(IT_0214*IT_0343 + IT_0213*IT_0344 +
       IT_0240*IT_0352 + IT_0239*IT_0354) + 6*IT_0358*IT_0359 + (IT_0240*IT_0268
       + IT_0239*IT_0271 + conj(IT_0212)*IT_0287 + IT_0212*IT_0290 + conj
      (IT_0265)*IT_0306 + IT_0265*IT_0309 + conj(IT_0253)*IT_0318 + IT_0253
      *IT_0321)*IT_0360 + (conj(IT_0284)*IT_0318 + IT_0284*IT_0321 + IT_0240
      *IT_0343 + IT_0239*IT_0344 + conj(IT_0212)*IT_0352 + IT_0212*IT_0354)
      *IT_0361 + (-6)*IT_0362*((conj(IT_0028) + -conj(IT_0080))*IT_0154 + 
      (IT_0028 + -IT_0080)*conj(IT_0154) + (-0.166666666666667)*IT_0219*IT_0339 
      + conj(IT_0092)*IT_0365 + IT_0092*IT_0368 + IT_0338*IT_0369) + (-6)*conj
      (IT_0219)*(IT_0337*IT_0365 + -1./6*IT_0154*IT_0370 + -1./6*IT_0092*IT_0371
       + IT_0314*IT_0376 + IT_0308*IT_0381) + (-6)*IT_0219*(IT_0337*IT_0368 +
       IT_0169*IT_0369 + -1./6*conj(IT_0154)*IT_0370 + -1./6*conj(IT_0092)
      *IT_0371 + IT_0317*IT_0376 + IT_0311*IT_0381) + (-6)*conj(IT_0295)*((
      -0.166666666666667)*IT_0169*IT_0295 + IT_0154*IT_0337 + (
      -0.166666666666667)*IT_0080*IT_0370 + (-0.166666666666667)*IT_0028*IT_0371
       + IT_0308*IT_0376 + IT_0314*IT_0381) + (-6)*IT_0295*(conj(IT_0154)
      *IT_0337 + (-0.166666666666667)*conj(IT_0080)*IT_0370 + (
      -0.166666666666667)*conj(IT_0028)*IT_0371 + IT_0311*IT_0376 + IT_0317
      *IT_0381);
    return create_ccomplex_return(IT_0382);
}

