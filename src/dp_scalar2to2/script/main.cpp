#include "scalar2to2.h"
#include "config.hpp"
using namespace scalar2to2;

int main(int argn, const char ** args) {
    std::cout << "dp_scalar2to2 version " << dp_scalar2to2_VERSION_MAJOR << '.' << dp_scalar2to2_VERSION_MINOR << '\n';
    std::cout << "Running " << args[0]<< " with " << argn << " argument(s).\n";
    std::cout << "I=" << I << '\n';
    return 0;
}