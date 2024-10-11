#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"

// #define USE_BOLTZVECTOR

using namespace scalar2to2;
using namespace advmath;

const std::string PATHPLOTS=std::string(OUTPATH)+"plots/" ;

inline double sq(const double x)
{
    return x*x;
}

int main(int argc, char ** argv)
{
 
    if(argc < 2)
    {
      std::cout << "This function requires 1 argument:\n"
                << " 1 - The name of the input file\n";
      return 1;
    }
    std::cout << std::setprecision(5) << std::scientific; 
    
    std::cout << "Reading inputs from " << argv[1] << std::endl;
    
    struct Param_t input(argv[1]);
    
    input.Print(); // Printing out the input structure

    AvgSvCalculator avg(input);
    
    std::cout << "\nDefining boltzfixed\n";
    BoltzmannSolver boltzfixed(input, avg.getProcList());
    std::cout << "Class constructed\n";
    double benchmark_omega = boltzfixed.relic_density();
    std::cout <<  "x_fo=" << input.getLightestBSMmass()/boltzfixed.Tfo ;
    std::cout << "\nOmega h^2 = " << benchmark_omega << '\n';

    std::cout << "Re-calculating by increasing m_chi by 10%\n";
    input.g_chi=0.4;
    input.m_chi=input.m_chi+200;
    input.AssignMassesVector();// <=== the criminal line!
    boltzfixed.changeInput(input);
    benchmark_omega = boltzfixed.relic_density();
    std::cout << input.m_chi;
    std::cout <<  "x_fo=" << input.getLightestBSMmass()/boltzfixed.Tfo ;
    std::cout << "\nOmega h^2 = " << benchmark_omega << '\n';
    
    const std::string filename=OUTPATH+static_cast<std::string>("scans/5th_scalar_scan_lin_m_chi2TeVm_phi5TeV.out");
    
    std::ofstream tempfile(filename.c_str());
        
    if(!tempfile)
    {
      std::cout << "Impossible to open " << filename << '\n'; 
      exit(1);
    }
    
    tempfile << std::setprecision(15) << std::scientific; 
    
    
    #ifdef USE_BOLTZVECTOR
      std::vector<BoltzmannSolver> boltztemp;
      boltztemp.reserve(100000);
    #else
      std::cout << "\nSolving the Boltzmann equation using the same instance of the Boltzmann class\n"; 
    #endif
    
//     for( input.g_u=1. ; input.g_u > std::numerical_limits<double>::min() ; input.g_u*=0.9)
    
    const double ref_mchi=2.e+3;
    const double ref_mphi=5.0e+3;
    const double ref_gg = 0.78*0.78;
    
    const double halfrange_mchi=1.0e+3;
    
    const double step_gchi=1.0e-2;
    const double step_mchi=2.0e+1;
    
    input.g_u =std::sqrt(ref_gg);// g_u is fixed
      for( input.g_chi=0.1 ; input.g_chi < 0.9 ; input.g_chi=input.g_chi + step_gchi)
        for( input.m_chi=ref_mchi - halfrange_mchi ; input.m_chi < ref_mchi + halfrange_mchi ; input.m_chi=input.m_chi + step_mchi)
        {
          input.AssignMassesVector();
      
          #ifdef USE_BOLTZVECTOR
            boltztemp.emplace_back(input);
            const double relicdensity=boltztemp.back().relic_density(), Tfo = boltztemp.back().Tfo;
          #else
            unsigned int exit_code = boltzfixed.changeInput(input);
            if( exit_code != 0 ) 
            {
              std::cout << "Failed to change to:\n" 
                        << " - g_u = " << input.g_u << '\n'
                        << " - g_chi = " << input.g_chi << '\n'
                        << " - m_chi = " << input.m_chi << '\n';
              continue;
            }
            
            const double relicdensity=boltzfixed.relic_density();
            const double Tfo = boltzfixed.Tfo;
          #endif
          const double xfo = input.getLightestBSMmass()/Tfo;
          tempfile << static_cast<double>(input.g_chi) << '\t' 
                   << static_cast<double>(input.m_chi) << '\t' 
                   << relicdensity    << '\t' 
                   << xfo             << '\n';
        }
          
   

    return 0;
}
