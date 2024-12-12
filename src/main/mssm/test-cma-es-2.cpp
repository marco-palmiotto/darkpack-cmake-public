#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp"
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "libcmaes/cmaes.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/esoptimizer.h>
#include <libcmaes/llogging.h>
using namespace scalar2to2;
using namespace advmath;

const std::string PATHPLOTS = std::string(OUTPATH);


using namespace libcmaes;


int main(int argc, char* argv[])
{
  FitFunc fsphere = [](const double* x, const int N)
  {
    double val = 0.0;
    for (int i = 0; i < N; i++)
      val += x[i] * x[i];
    return val;
  };

  int dim = 10;                      // problem dimensions.
  double sigma = 0.1;                // initial step-size, i.e. estimated initial parameter error.
  std::vector<double> x0(dim, 10.0); // initialize x0 as 10.0 in all 10 dimensions
  CMAParameters<> cmaparams(x0, sigma);
  cmaparams.set_algo(aCMAES); // select active CMA-ES as algorithm (default is CMA-ES).

  cmaparams.set_max_fevals(10000); // limits the number of function evaluations to 10000
  cmaparams.set_max_iter(100000);  // limits the number of iterations of the algorithms to 1000000
  cmaparams.set_ftarget(1e-8);     // stops the optimization whenever the objective function values gets below 1e-8

  CMASolutions cmasols = cmaes<>(fsphere, cmaparams);

  std::cout << "best solution: " << cmasols << std::endl;
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << std::endl;
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  std::cout << cmasols.run_status(); // the optimization status, failed if < 0

  return 0;
}