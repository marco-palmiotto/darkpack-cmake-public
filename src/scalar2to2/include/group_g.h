#ifndef CSL_LIB_scalar2to2_G_H_INCLUDED
#define CSL_LIB_scalar2to2_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_scalar2to2.h"

namespace scalar2to2 {


extern const std::array<Callable<complex_t, param_t>, 14> f_G;
extern const std::map<std::string, Callable<complex_t, param_t>> fmap_G;

}
 // End of namespace scalar2to2

#endif
