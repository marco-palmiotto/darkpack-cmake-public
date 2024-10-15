#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"

using namespace scalar2to2;
using namespace advmath;

const std::string PATHPLOTS=std::string(OUTPATH)+"plots/" ; 

scalar2to2::RunningSM run;

inline double sq(const double x)
{
    return x*x;
}

double gamma_phi_pred_onlyup(const Param_t &input)
{
  const double v_h = std::pow(2., -0.25)*std::pow(input.Gfermi, -0.5);
  
  std::vector<int> upquarks = {corr::u, corr::c, corr::t};
  double quarkcontrib = 0.;
  
  for( int quark : upquarks )
  {
    double mf= input.masses_vector[quark];
    double Y = std::sqrt(2.)*mf/v_h;
    quarkcontrib+= 3.*sq(Y*input.g_u)*input.m_phi/16./M_PI*std::pow( 1. - 4.*sq(mf/input.m_phi), 3./2.);
  }
  
  return quarkcontrib;
}

// Computation of Gamma_phi
double gamma_phi_pred(const Param_t &input)
{
//   const double v_h = std::pow(2., -0.25)*std::pow(input.Gfermi, -0.5);
  
  const double quarkcontrib = gamma_phi_pred_onlyup(input);
  
//   input.Print();
  
  const double chicontrib= sq(input.g_chi)*input.m_phi/8./M_PI*std::pow( 1. - 4.*sq(input.m_chi/input.m_phi), 3./2.);


//   printf("%s = %.3e ", "input.g_chi", double(input.g_chi));
//   
//   std::cout << "chicontrib=" << chicontrib << std::endl;
  
  const double strongcontrib=0.;
  
  return quarkcontrib + chicontrib + strongcontrib;
}

// Computation of the sum of the squared amplitudes from a given formula (the one I found)
double sum_Squaredampl_pred( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double C=sq(mf*gf*input.g_chi*input.e_em/(input.m_W*std::sin(input.theta_W)));
  // Defining s
  const double s = Ecm*Ecm;
  
  
  double result =Ncolors*C;
  
  // Writing the numerator of the result
  result*=(s-4*sq(input.m_chi))*(s-4*sq(mf));
  
  // Writing the denominator of the result
  result/=( sq(s-sq(input.m_phi)) + sq(input.m_phi*input.Gamma_phi));
  
  return result;
}

// This function computes the contribution to dWeff/dcostheta
// for a single process (i.e. a single SM fermion),
// starting from the |M|^2 expression defined in DarkPACK 
double dW_dcos_pred( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 =sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);
  
  const double s = sq(Ecm);
  const double p12= 0.5*std::sqrt(s-4.*sq(input.m_chi));
  const double p34= 0.5*std::sqrt(s-4.*sq(mf));
  const double peff = p12;
  
  double result = p12*p34*M2;
  result/= 8.*M_PI*sq(corr::part_hel_dof[corr::chi])*peff*Ecm;
 
  
  return result;
}

// This function computes the contribution to dWeff/dcostheta without the global coefficient, 
// for a single process (i.e. a single SM fermion),
// starting from the |M|^2 expression defined in DarkPACK 
double dW_dcos_pred_nocoeff( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 =sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);
  
  const double s = sq(Ecm);
  const double p12= 0.5*std::sqrt(s-4.*sq(input.m_chi));
  const double p34= 0.5*std::sqrt(s-4.*sq(mf));
  
  return p12*p34*M2;
}

// This function computes the total Weff starting from the |M|^2 expression defined in DarkPACK 
double Weff_pred_fromM2( Param_t &input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);
  
  std::cout << "\nIn Weff_pred_fromM2, " << input.masses_vector[corr::c] << " " << input.m_c << " " <<
               input.m_c_m_c << "\n";
  
  double sum=0.;
  
  const std::array<int, 9> 
        part({corr::u ,corr::c,  corr::t,
              corr::d, corr::s,  corr::b,
              corr::e, corr::mu, corr::tau}),
        part_color({3,3,3,
                    3,3,3,
                    1,1,1});
  const std::array<double, 9> 
        part_gf({input.g_u,input.g_u,input.g_u,
                input.g_d,input.g_d,input.g_d,
                input.g_l,input.g_l,input.g_l});
         
  size_t counter=0;
  
  for(auto i = 0 ; i <9; i++)
    {
      Process2to2 singleproc({corr::chi, corr::chi, part[i],part[i]},{false, true, false,true}); 
      if(singleproc.checkExistance())
      {
        double contrib=dW_dcos_pred(input, Ecm, part_color[i], part_gf[i], input.masses_vector[part[i]]);
        sum+=contrib;
        std::cout << "For " << singleproc.getName() << " dW/dcos(" << Ecm << ")= " << contrib << '\n';
        counter++;
      }
    }
    
  std::cout << "In the sum there were " << counter << " processes.\n";
  
  return 2*sum;
}


// This function computes Weff using the formula derived by hand
double Weff_formula( Param_t &input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);
  
//   std::cout << input;
  
//   std::cout << "Defining prefactor\n";
  const double C = sq(input.g_chi.get()*input.e_em/input.m_W/std::sin(input.theta_W));
  
//   std::cout << "Defining s\n";
  const double s = Ecm*Ecm;
  
  double result =C/(8.*M_PI*Ecm*sq(input.getLightestBSMdof()))*(s-4.*sq(input.m_chi));
  // Writing the denominator of the result
  result/=( sq(s-sq(input.m_phi)) + sq(input.m_phi*input.Gamma_phi));
  
  double sum=0.;
  
    const std::array<int, 9> part({corr::u, corr::c, corr::t,corr::d, corr::s, corr::b,corr::e, corr::mu, corr::tau}),
                             part_color({3,3,3,3,3,3,1,1,1});
    const std::array<double, 9> part_gf({input.g_u,input.g_u, input.g_u,input.g_d,input.g_d,input.g_d,input.g_l,input.g_l,input.g_l});
//   std::cout << "calculating the sum over all fermions\n";
  for(auto i = 0 ; i <9; i++)
  {
    Process2to2 singleproc({corr::chi, corr::chi, part[i],part[i]},{false, true, false,true}); 
    if(singleproc.checkExistance())
    {
      sum+=sq(1.0*input.masses_vector[part[i]]*part_gf[i])*part_color[i]*std::pow(s-4.*sq(input.masses_vector[part[i]]), 3./2.);
    }
  }
    
//   std::cout << "Multiplying the sum\n";
  result*=sum;
  
  return result;
}


int main(int argc, char ** argv)
{
 
    if(argc < 2)
    {
      std::cout << "This function requires 1 argument\n";
      return 1;
    }
    std::cout << std::setprecision(5) << std::scientific; 
    
    std::cout << "Reading inputs from " << argv[1] << std::endl;
    
    struct Param_t input(argv[1]);
    
    input.Print(); // Printing out the input structure
    
    // run.init(input);
    // run.RunCharmMass(false);
    
    // std::cout << "v = " << 2*(input.m_W*std::sin(input.theta_W)/input.e_em);
    
    std::cout << "\nDefining BoltzmannSolver\n";
    BoltzmannSolver boltz(input);

    std::cout << "Oh2 =" << boltz.relic_density() << std::endl;
    
    
    return 0;
}
