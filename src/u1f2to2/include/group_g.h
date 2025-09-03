#ifndef CSL_LIB_u1f2to2_G_H_INCLUDED
#define CSL_LIB_u1f2to2_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_u1f2to2.h"

namespace u1f2to2 {


extern const std::array<Callable<complex_t, param_t>, 38> f_G;
extern const std::map<std::string, Callable<complex_t, param_t>> fmap_G;

}
 // End of namespace u1f2to2

#endif
