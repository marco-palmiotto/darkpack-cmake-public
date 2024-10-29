#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include "libcmaes/cmaes.h"
#include <libcmaes/esoptimizer.h>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/llogging.h>
#include <fstream>
#include <random>
using namespace scalar2to2;
using namespace advmath;
using namespace libcmaes;

const std::string PATHPLOTS=std::string(OUTPATH);
constexpr const double omega_h2_target=0.12;
constexpr const double omega_h2_err=0.00091;

std::string filename="";

static inline double sq(const double x)
{
  return x*x;
}

static inline double gaussian_shape(const double val, const double mean, const double stddev)
{
  return -std::exp( - sq(val-mean)/2./sq(stddev) );
}

double function_to_minimize(const double *parameters, const int npar=3)
{
  /* 
    BLOCK SCALARMODEL
      1  0   # g_d
      2  7.80000000E-1   # g_u    [0 -1]
      3  0   # g_l
      4  7.8000000E-1     # g_chi
      5  2.00000000E+03   # m_chi [100 GeV - 1 TeV]
      6  4.40000000E+03   # m_phi
  */
  // static unsigned int ncalls=0;
  static Param_t input(filename.c_str());
  input.g_d=0.;
  input.g_u=1.;
  input.refresh();
  static BoltzmannSolver boltz(input);

  // parameters in input will be
  input.g_chi = parameters[0];
  input.m_chi = parameters[1];
  input.m_phi = parameters[2];

  input.refresh();
  boltz.changeInput(input);
  const double omega=boltz.relic_density();
  return std::abs(omega - omega_h2_target);
}


int main(int argc, char *argv[])
{

  std::cout << "Running " << argv[0] << '\n';

  constexpr const double sigma=omega_h2_err;
  constexpr const int n_required_args=2;
  constexpr const double factor=2.2;
  constexpr const int npar=1;
  const double lbounds[1] ={1.0e-20},
               ubounds[1] ={1.0};

  constexpr const double init_seed_g_chi=0.3;
  constexpr const double init_m_chi=1.0e+3;
  constexpr const double end_m_chi=1.e+5;
  constexpr const int npoints=1000;
  const double stepsize=std::pow(end_m_chi/init_m_chi, 1./npoints);

  if(argc < n_required_args)
  {
    std::cout << "This function requires " << n_required_args << " arguments:\n"
              << " - The name of the input file\n"
              << " - The name of the output file\n";
    return 1;
  }
  filename=argv[1];
  
  std::cout << "The input file is " << filename << '\n';

  std::vector<double> x0;
  x0.reserve(npar);
  x0.emplace_back(init_seed_g_chi);

  if(x0.size() != npar)
  {
    std::cerr << "Initial parameters have wrong size: x0.size()=" << x0.size() << ", npar=" << npar << '\n';
    return 2;
  } 

  /* 
    BLOCK SCALARMODEL
    1  0   # g_d
    2  7.80000000E-1   # g_u    [0 -1]
    3  0   # g_l
    4  7.8000000E-1     # g_chi
    5  2.00000000E+03   # m_chi [100 GeV - 1 TeV]
    6  4.40000000E+03   # m_phi
  */ 
  Param_t input(filename.c_str());
  input.g_d=0.;
  input.g_u=1.;
  input.m_chi=init_m_chi;
  input.m_phi=factor*input.m_chi;
  input.refresh();
  BoltzmannSolver boltz(input);

  auto function_to_minimize = [&] (const double *parameters, const int npar=3)
  {
    // parameters in input will be
    input.g_chi = parameters[0];
    input.refresh();
    boltz.changeInput(input);
    const double omega=boltz.relic_density();
    return std::abs(omega - omega_h2_target);
  };

  GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,npar); // genotype / phenotype transform associated to bounds.

  // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator.           
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,sigma,-1,0,gp);     
  
  cmaparams.set_algo(aCMAES);
  FitFunc f = function_to_minimize;
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams);

  std::cout << "seed g_chi is\n";
  for(int i=0; i < npar ; i++) std::cout <<  x0[i] << '\n';
  std::cout << "\nbest solution: ";
  cmasols.print(std::cout,0,gp);
  std::cout << '\n';
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << '\n';
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";


  std::ofstream outfile{argv[2]};
  if(!outfile)
  {
    std::cerr << "Impossible to open " << argv[2] << " exiting\n";
    return 0;
  } 

  Candidate bcand = cmasols.best_candidate();
  double fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
  std::vector<double> x_stdv = bcand.get_x(); // vector of objective function parameters at minimum.
  const double* x_dptr = bcand.get_x_ptr(); // vector of objective function parameters at minimum, as C-style double array
  Eigen::VectorXd x_ev = bcand.get_x_dvec(); // vector of objective function parameters at minimum, as Eigen vector
  double edm = cmasols.edm(); // expected distance to the minimum.
  
  Eigen::VectorXd bestparameters = gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec());
  // Applying pheno to the best candidate output

  std::cout << gp.pheno(x_ev).transpose();
  std::cout << bestparameters.transpose();
  
  std::cout << "# m_chi   m_phi    g_chi   pull-Oh2\nRESULT=";
  std::cout << input.m_chi.get() << '\t' << input.m_phi.get() << '\t' << x_dptr[0] << '\t' << fmin << '\n';
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(std::log(lbounds[0]), std::log(ubounds[0]));

  while(input.m_chi < end_m_chi)
  {
    if( x_dptr[0] > ubounds[0] ||  x_dptr[0] < lbounds[0]) 
    {
      std::cerr << "Error: result is out of bound\n";
      return 1;
    }
    input.m_chi = input.m_chi*stepsize;
    x0[0]=(fmin < 1.0e-3) ? x_dptr[0] : std::exp(dis(gen));

    CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams1(x0,sigma,-1,0,gp); // -1 for automatically \
  decided lambda, 0 is for random	seeding	of the internal generator.                                              
    cmaparams1.set_algo(aCMAES);
    cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams1);

    bcand = cmasols.best_candidate();
    fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
    x_stdv = bcand.get_x(); // vector of objective function parameters at minimum.
    x_dptr = bcand.get_x_ptr(); // vector of objective function parameters at minimum, as C-style double array
    x_ev = bcand.get_x_dvec(); // vector of objective function parameters at minimum, as Eigen vector
    edm = cmasols.edm(); // expected distance to the minimum.

    std::cout << "RESULT=" <<input.m_chi.get() << '\t' << input.m_phi.get() << '\t' << x_dptr[0] << '\t' << fmin << '\n';
  }

  return cmasols.run_status();
}