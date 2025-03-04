#ifndef LESHOUCESFROMMARTY_INCLUDED
#define LESHOUCESFROMMARTY_INCLUDED
// NOTE: DO NOT replace the include guard with pragma once!
//       The constructor Param_t(std::string) needs the macro to be defined

#ifndef __cplusplus
  #define __cplusplus
#endif

#include "RunningSM.hpp"
#include "config.hpp"
#include "correspondance.hpp"
#include "marty/lha/lha.h"
#include "params_new.hpp"


namespace u1f2to2::readmodule
{

  inline void computeThetaWandMw(const double alpha_em, const double Gf, const double Mz, double& thetaW, double& Mw)
  {
    thetaW = 0.5 * std::asin(std::sqrt(4 * M_PI * alpha_em / (std::sqrt(2) * Gf * Mz * Mz)));
    Mw = Mz * std::cos(thetaW);
    return;
  };

  void AssignSMparams(Param_t& param, const std::string namefile);
  void ReadBSMparams(Param_t& param, const std::string filename);
  void AssignBSMparams(Param_t& param);

  Param_t ReadLHA(const std::string name);

} // namespace u1f2to2::readmodule
#endif