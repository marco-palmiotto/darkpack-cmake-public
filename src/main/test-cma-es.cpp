#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include <libcmaes/esoptimizer.h>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/llogging.h>

using namespace scalar2to2;
using namespace advmath;

const std::string PATHPLOTS=std::string(OUTPATH);


using namespace libcmaes;

FitFunc cigtab = [](const double *x, const int N)
{
  int i;
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)
    sum += x[i]*x[i];
  return sum;
};

int main(int argc, char *argv[])
{
#ifdef HAVE_GLOG
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr=1;
  google::SetLogDestination(google::INFO,"");
  //FLAGS_log_prefix=false;
#endif

  int dim = 10;
  std::vector<double> x0 = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double sigma = 0.2;
  int lambda = 10;
  CMAParameters<> cmaparams(x0,sigma,lambda);
  ESOptimizer<CMAStrategy<CovarianceUpdate>,CMAParameters<>> cmaes(cigtab,cmaparams);
  cmaes.optimize();
  double edm = cmaes.edm();
  std::cerr << "EDM=" << edm << " / EDM/fm=" << edm / cmaes.get_solutions().best_candidate().get_fvalue() << std::endl;
}