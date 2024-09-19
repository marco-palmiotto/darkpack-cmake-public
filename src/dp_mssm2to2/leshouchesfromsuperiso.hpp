#pragma once

#include "input_manip.hpp"
#include "leshouches.h"
#include "params_new.hpp"

namespace mssm2to2::superisosupport
{
  int InitInterfaceStruct(const struct parameters * const param, struct param_t &input);

  int InitInterfaceStruct_Full(const struct parameters * const param, struct Param_t &input);

  struct parameters ReadLHA(struct Param_t &input, const char * name, int *err);
}
