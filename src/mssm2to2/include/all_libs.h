#ifndef ALL_LIBS_MSSM_INCLUDED
#define ALL_LIBS_MSSM_INCLUDED
#ifndef cplusplus__
#define cplusplus__
#endif

#include "mssm2to2.h"
#include "correspondance.h"
#include "advmath.h"
#include "params_new.h"
#include "compare.h"
#include "input_manip.h"
#include "leshouches.h"
#include "leshouchesfrommarty.h"
#include "leshouchesfromsuperiso.h"
#include "process.h"
#include "avgsvcalculator.h"
#include "RunningSM.h"
#include "boltzmann.h"

namespace mssm2to2{
std::string getSname(const Process2to2 &proc);
}
#endif











