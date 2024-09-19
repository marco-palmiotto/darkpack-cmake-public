namespace scalar2to2
{
inline void update_kinematics(Param_t &input, const double sij[5][5])
{
  input.s_12 = sij[1][2];
  input.s_34 = sij[3][4];
}
inline cparam_s copy_to_c_struct(Param_t &input)
{
  cparam_s output;
#define CONVERT_COMPLEX(X) ( (static_cast<std::complex<double>>(X)).real() + _mty_I*(static_cast<std::complex<double>>(X)).imag() )
  output.pi = input.pi;
  output.g_d = input.g_d;
  output.g_l = input.g_l;
  output.g_u = input.g_u;
  output.V_cb = input.V_cb;
  output.V_tb = input.V_tb;
  output.V_ud = input.V_ud;
  output.V_us = input.V_us;
  output.e_em = input.e_em;
  output.s_12 = input.s_12;
  output.s_34 = input.s_34;
  output.g_chi = input.g_chi;
  output.theta_W = input.theta_W;
  output.reg_prop = input.reg_prop;
  output.Gamma_phi = input.Gamma_phi;
  output.g_s = input.g_s;
  output.m_W = input.m_W;
  output.m_Z = input.m_Z;
  output.m_b = input.m_b;
  output.m_c = input.m_c;
  output.m_d = input.m_d;
  output.m_e = input.m_e;
  output.m_h = input.m_h;
  output.m_s = input.m_s;
  output.m_t = input.m_t;
  output.m_u = input.m_u;
  output.m_mu = input.m_mu;
  output.m_chi = input.m_chi;
  output.m_phi = input.m_phi;
  output.m_tau = input.m_tau;
  output.V_cd = CONVERT_COMPLEX(input.V_cd);
  output.V_cs = CONVERT_COMPLEX(input.V_cs);
  output.V_td = CONVERT_COMPLEX(input.V_td);
  output.V_ts = CONVERT_COMPLEX(input.V_ts);
  output.V_ub = CONVERT_COMPLEX(input.V_ub);
#undef CONVERT_COMPLEX
  return output;
}
inline cparam_s update_c_kinematics(const Param_t &input, cparam_s output)
{
  output.s_12 = input.s_12;
  output.s_34 = input.s_34;
  return output;
}
} // End of namespace __SPEC_LIB_NAME__
