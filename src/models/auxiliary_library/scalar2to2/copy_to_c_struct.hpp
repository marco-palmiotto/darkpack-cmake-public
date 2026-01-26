namespace scalar2to2
{
inline void update_kinematics(Param_t &input, const double sij[5][5])
{
}
inline cparam_s copy_to_c_struct(Param_t &input)
{
  cparam_s output;
#define CONVERT_COMPLEX(X) ( (static_cast<std::complex<double>>(X)).real() + _mty_I*(static_cast<std::complex<double>>(X)).imag() )
#undef CONVERT_COMPLEX
  return output;
}
inline cparam_s update_c_kinematics(const Param_t &input, cparam_s output)
{
  return output;
}
} // End of namespace scalar2to2
