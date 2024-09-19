#include "MathLib/add.hpp"
#include "MathLib/subtract.hpp"
#include "MathLib/modulo.h"
#include "MathLib/multiply.h"
#include "FinalLib/entity.h"
#include <iostream>

int main()
{

  int a=2, b=3, c=5;

  std::cout << "Using add library: 2+3= " << MathLib::add(a,b) << '\n';

  std::cout << "Using subtract library: 2-3= " << MathLib::subtract(a,b) << '\n';

  std::cout << "Using add C library: 5 % 3= " << modulo(c,b) << '\n';
  
  std::cout << "Using multiply C library: 5 * 3= " << multiply(c,b) << '\n';

  spec::entity_t e(0,0,0);
  spec::display(e);
  spec::move(e, 1, -1, 0);
  spec::display(e);
  return 0;
   
}
