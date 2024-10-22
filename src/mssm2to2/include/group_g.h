#ifndef CSL_LIB_mssm2to2_G_H_INCLUDED
#define CSL_LIB_mssm2to2_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_mssm2to2.h"

namespace mssm2to2 {


extern const std::array<Callable<complex_t, param_t>, 3456> f_G;


}
 // End of namespace mssm2to2

#endif
