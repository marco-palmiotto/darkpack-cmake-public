#include "scalar2to2.h"
#include "config.hpp"
using namespace scalar2to2;

int main(int argn, const char ** args) {
    std::cout << "dp_scalar2to2 version " << __SPEC_LIB_NAME__VERSION_MAJOR << '.' << __SPEC_LIB_NAME__VERSION_MINOR << '\n';
    std::cout << "Running " << args[0]<< " with " << argn << " argument(s).\n";
    std::cout << "I=" << I << '\n';
    return 0;
}