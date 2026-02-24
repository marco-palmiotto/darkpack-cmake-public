#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include "dp_scalar2to2/correspondance.hpp"
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace scalar2to2;
using namespace advmath;

const std::string PATHPLOTS = std::string(OUTPATH) + "plots/"; //<- Path to the plots directory

RunningSM run; //<- Instance of the RunningSM class to handle running parameters

inline double sq(const double x) { return x * x; } //<- Function to compute the square of a number

/*
 * @brief Computes the decay width of phi into up-type quarks only.
 * @param input The parameters of the model
 */
double gamma_phi_pred_onlyup(const Param_t& input)
{
  const double v_h = std::pow(2., -0.25) * std::pow(input.Gfermi, -0.5); //<- Higgs vev in GeV

  const std::array<int, 3> upquarks = {corr::u, corr::c, corr::t}; //<- List of up-type quarks
  double quarkcontrib = 0.; //<- Variable to accumulate the contribution from up-type quarks

  for (const int quark : upquarks)
  {
    const double mf = input.masses_vector[quark]; //<- Mass of the current up-type quark
    const double Y = std::sqrt(2.) * mf / v_h;    //<- Yukawa coupling of the current up-type quark
    quarkcontrib +=
        3. * sq(Y * input.g_u) * input.m_phi / 16. / M_PI * std::pow(1. - 4. * sq(mf / input.m_phi), 3. / 2.);
  }

  return quarkcontrib;
}

// Computation of Gamma_phi
double gamma_phi_pred(const Param_t& input)
{
  //   const double v_h = std::pow(2., -0.25)*std::pow(input.Gfermi, -0.5);

  const double quarkcontrib = gamma_phi_pred_onlyup(input);

  //   input.Print();

  const double chicontrib =
      sq(input.g_chi) * input.m_phi / 8. / M_PI * std::pow(1. - 4. * sq(input.m_chi / input.m_phi), 3. / 2.);


  //   printf("%s = %.3e ", "input.g_chi", double(input.g_chi));
  //
  //   std::cout << "chicontrib=" << chicontrib << std::endl;

  const double strongcontrib = 0.; // This contributes at 1-loop, therefore we do not include it here

  return quarkcontrib + chicontrib + strongcontrib;
}

// Computation of the sum of the squared amplitudes from a given formula (the one I found)
double sum_Squaredampl_pred(Param_t& input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double C = sq(mf * gf * input.g_chi * input.e_em / (input.m_W * std::sin(input.theta_W)));
  // Defining s
  const double s = Ecm * Ecm;


  double result = Ncolors * C;

  // Writing the numerator of the result
  result *= (s - 4 * sq(input.m_chi)) * (s - 4 * sq(mf));

  // Writing the denominator of the result
  result /= (sq(s - sq(input.m_phi)) + sq(input.m_phi * input.Gamma_phi));

  return result;
}

// This function computes the contribution to dWeff/dcostheta
// for a single process (i.e. a single SM fermion),
// starting from the |M|^2 expression defined in DarkPACK
double dW_dcos_pred(Param_t& input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 = sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);

  const double s = sq(Ecm);
  const double p12 = 0.5 * std::sqrt(s - 4. * sq(input.m_chi));
  const double p34 = 0.5 * std::sqrt(s - 4. * sq(mf));
  const double peff = p12;

  double result = p12 * p34 * M2;
  result /= 8. * M_PI * sq(corr::part_hel_dof[corr::chi]) * peff * Ecm;


  return result;
}

// This function computes the contribution to dWeff/dcostheta without the global coefficient,
// for a single process (i.e. a single SM fermion),
// starting from the |M|^2 expression defined in DarkPACK
double dW_dcos_pred_nocoeff(Param_t& input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 = sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);

  const double s = sq(Ecm);
  const double p12 = 0.5 * std::sqrt(s - 4. * sq(input.m_chi));
  const double p34 = 0.5 * std::sqrt(s - 4. * sq(mf));

  return p12 * p34 * M2;
}

// This function computes the total Weff starting from the |M|^2 expression defined in DarkPACK
double Weff_pred_fromM2(Param_t& input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);

  // std::cout << "\nIn Weff_pred_fromM2, " << input.masses_vector[corr::c] << " " << input.m_c << " " << input.m_c_m_c
  //           << "\n";

  double sum = 0.;

  const std::array<int, 9> part({corr::u, corr::c, corr::t, corr::d, corr::s, corr::b, corr::e, corr::mu, corr::tau}),
      part_color({3, 3, 3, 3, 3, 3, 1, 1, 1});
  const std::array<double, 9> part_gf(
      {input.g_u, input.g_u, input.g_u, input.g_d, input.g_d, input.g_d, input.g_l, input.g_l, input.g_l});

  size_t counter = 0;

  for (auto i = 0; i < 9; i++)
  {
    Process2to2 singleproc({corr::chi, corr::chi, part[i], part[i]}, {false, true, false, true});
    if (singleproc.checkExistance())
    {
      double contrib = dW_dcos_pred(input, Ecm, part_color[i], part_gf[i], input.masses_vector[part[i]]);
      sum += contrib;
      // std::cout << "For " << singleproc.getName() << " dW/dcos(" << Ecm << ")= " << contrib << '\n';
      counter++;
    }
  }

  // std::cout << "In the sum there were " << counter << " processes.\n";

  return 2 * sum;
}


// This function computes Weff using the formula derived by hand
double Weff_formula(Param_t& input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);

  //   std::cout << input;

  //   std::cout << "Defining prefactor\n";
  const double C = sq(input.g_chi.get() * input.e_em / input.m_W / std::sin(input.theta_W));

  //   std::cout << "Defining s\n";
  const double s = Ecm * Ecm;

  double result = C / (8. * M_PI * Ecm * sq(input.getLightestBSMdof())) * (s - 4. * sq(input.m_chi));
  // Writing the denominator of the result
  result /= (sq(s - sq(input.m_phi)) + sq(input.m_phi * input.Gamma_phi));

  double sum = 0.;

  const std::array<int, 9> part({corr::u, corr::c, corr::t, corr::d, corr::s, corr::b, corr::e, corr::mu, corr::tau}),
      part_color({3, 3, 3, 3, 3, 3, 1, 1, 1});
  const std::array<double, 9> part_gf(
      {input.g_u, input.g_u, input.g_u, input.g_d, input.g_d, input.g_d, input.g_l, input.g_l, input.g_l});
  //   std::cout << "calculating the sum over all fermions\n";
  for (auto i = 0; i < 9; i++)
  {
    Process2to2 singleproc({corr::chi, corr::chi, part[i], part[i]}, {false, true, false, true});
    if (singleproc.checkExistance())
    {
      sum += sq(1.0 * input.masses_vector[part[i]] * part_gf[i]) * part_color[i] *
             std::pow(s - 4. * sq(input.masses_vector[part[i]]), 3. / 2.);
    }
  }

  //   std::cout << "Multiplying the sum\n";
  result *= sum;

  return result;
}


int main(int argc, char** argv)
{

  if (argc < 2)
  {
    std::cout << "This function requires 1 argument:\n"
              << " 1. The name of the input file with the parameters.\n";
    return 1;
  }
  std::cout << std::setprecision(5) << std::scientific;

  std::cout << "Reading inputs from " << argv[1] << std::endl;

  struct Param_t input(argv[1]);

  input.Print(); // Printing out the input structure

  // In this demo, we're comparing the results of the DarkPACK functions
  // with the results of the formulas derived by hand, or present in literature.
  // In order to compute the latters, we initalise the instance of the global
  // RunningSM class.
  run.init(input);
  run.RunCharmMass(false);

  std::cout << "The Higgs vev is:\n";
  std::cout << "v = " << 2 * (input.m_W * std::sin(input.theta_W) / input.e_em);

  std::cout << "\nDefining the instance of AvgSvCalculator\n";
  AvgSvCalculator avgsvcalc(input);
  avgsvcalc.setWeffcuts(false);
  std::cout << "AvgSvCalculator object instantiated\n";

  // Computing single 2 to 2 sum of M^2
  double sqrts = 10. * input.m_chi;
  constexpr double costheta = 0.5;
  constexpr const std::array<int, 9> part(
      {corr::u, corr::c, corr::t, corr::d, corr::s, corr::b, corr::e, corr::mu, corr::tau}),
      part_color({3, 3, 3, 3, 3, 3, 1, 1, 1});
  const std::array<double, 9> part_gf(
      {input.g_u, input.g_u, input.g_u, input.g_d, input.g_d, input.g_d, input.g_l, input.g_l, input.g_l});
  std::cout << std::setprecision(5) << std::scientific;
  constexpr int width_field = 18;
  for (auto i = 0; i < 9; i++)
  {
    Process2to2 singleproc({corr::chi, corr::chi, part[i], part[i]}, {false, true, false, true});
    if (!singleproc.checkExistance())
      continue;

    std::cout << '\n' << singleproc;
    run.HandleParamRunning(input, sqrts);
    Param_t input_m(input); //<- Copy of the input to pass to the Process2to2 instance

    // Optimisation for the running
    singleproc.setRunningData(&run);
    singleproc.setRunningExternal();

    const double marty = singleproc.getSumSquaredAmpl(input_m, sqrts, costheta); //<- Prediction from MARTY
    const double pred = sum_Squaredampl_pred(input, sqrts, part_color[i], part_gf[i],
                                             input.masses_vector[part[i]]); //<- Prediction from the formula

    const double martydweff = singleproc.get_g2_dweff_dcos(input_m, sqrts, 0.5) /
                              sq(corr::part_hel_dof[corr::chi]); //<- Prediction from DarkPACK
    const double predweff = dW_dcos_pred(input, sqrts, part_color[i], part_gf[i],
                                         input.masses_vector[part[i]]); //<- Prediction from the formula

    // Print the header row
    std::cout << std::setw(width_field) << ".getSumSquaredAmpl" << '\t' << std::setw(width_field)
              << "sum_Squaredampl_pred" << '\t' << std::setw(width_field) << "ratio" << '\t' << std::setw(width_field)
              << ".get_g2_dweff_dcos/g2" << '\t' << std::setw(width_field) << "dW_dcos_pred" << '\t'
              << std::setw(width_field) << "ratio" << '\n';
    std::cout << std::setw(width_field) << marty << '\t' << std::setw(width_field) << pred << '\t'
              << std::setw(width_field) << marty / pred << '\t' << std::setw(width_field) << martydweff << '\t'
              << std::setw(width_field) << predweff << '\t' << std::setw(width_field) << martydweff / predweff << '\n';
  }

  run.HandleParamRunning(input, sqrts);

  const double weff_from_darkpack =
      avgsvcalc.getWeff(sqrts); //<- Weff computed by DarkPACK, using the numerical integration
  std::cout << "\nWeff (DarkPACK) = " << weff_from_darkpack << std::endl; // Uses multithreading

  const double weff_formula_by_hand = Weff_formula(input, sqrts);
  std::cout << "Weff (full formula by hand) = " << weff_formula_by_hand << std::endl;
  const double weff_with_kinematics_by_hand = Weff_pred_fromM2(input, sqrts);
  std::cout << "Weff (from formula with M2 from MARTY) = " << weff_with_kinematics_by_hand << std::endl;

  std::cout << "Ratio Weff/formula  = " << weff_from_darkpack / weff_formula_by_hand << std::endl;

  std::cout << "Defining a low temperature, to compare computed and predicted values of <σv>\n";
  double T = 0.001 * input.m_chi;
  const double sigmav = avgsvcalc.getAverageSigmav(T);
  std::cout << "<σv>( " << T << " GeV) = " << sigmav << " GeV^-2\n\n";

  // Lambda function for computing the value of <σv> predicted by the literature
  auto sigmav_pred = [](const double Temp, const int Nc, const int part_l, const double gf, Param_t& input_l) mutable
  {
    run.HandleParamRunning(input_l, input_l.getLightestBSMmass());
    // Higgs vev in the SM
    const double v_h = std::pow(2., -0.25) * std::pow(input_l.Gfermi, -0.5);
    const double Y_f = std::sqrt(2.) * input_l.masses_vector[part_l] / v_h;
    // The temperature corresponds to the center of mass energy
    double numerator = Nc * 3. * sq(input_l.g_chi * gf * Y_f) * input_l.m_chi * Temp *
                       std::pow(1. - sq(input_l.masses_vector[part_l] / input_l.m_chi), 3. / 2.);
    double denom = 8 * M_PI * (sq(sq(input_l.m_phi) - 4. * sq(input_l.m_chi)) + sq(input_l.m_phi * input_l.Gamma_phi));

    return numerator / denom;
  };

  // Computing the value according to the literature
  double sigmavpred = 0.;

  for (auto count_part : {corr::e, corr::mu, corr::tau})
    sigmavpred += sigmav_pred(T, 1, count_part, input.g_l, input);

  // std::cout << "Computed sigmav for leptons\n";

  for (auto count_part : {corr::d, corr::s, corr::b})
    sigmavpred += sigmav_pred(T, 3, count_part, input.g_d, input);

  // std::cout << "Computed sigmav for d-type\n";

  for (auto count_part : {corr::u, corr::c, corr::t})
    sigmavpred += sigmav_pred(T, 3, count_part, input.g_u, input);

  // std::cout << "Computed sigmav for u-type\n";

  // Comparing the values
  std::cout << "The predicted sigmav is " << sigmavpred << " GeV^-2\n";
  std::cout << "The ratio is " << sigmavpred / sigmav << std::endl;

  // Computing the decay with of φ at the tree-level
  const double gamma_phi_up = width_phi(input).real();
  std::cout << "\nΓ_φ_up_computed = " << gamma_phi_up;
  std::cout << "\nΓ_φ_up_pred = " << gamma_phi_pred(input);


  // Let us plot W_eff
  // Creating plot files
  std::ofstream fout{PATHPLOTS + "Weff.out"};
  if (!fout)
  {
    std::cerr << "Impossible to open " << PATHPLOTS + "Weff.out" << std::endl;
    exit(1);
  }

  std::cout << "\n\nOpened " << PATHPLOTS << "Weff.out\n";

  fout << std::scientific << std::setprecision(5);

  std::vector<double> Ecmv;
  Ecmv.reserve(avgsvcalc.getWefftabsize());
  for (size_t i = 0; i < avgsvcalc.getWefftabsize(); i++)
  {

    double peff, g2_weff;
    avgsvcalc.get_g2_WeffTabElement(i, peff, Ecmv[i], g2_weff);

    const double Ecm = Ecmv[i];

    run.HandleParamRunning(input, Ecm);

    fout << std::scientific << Ecm << '\t' << 0.5 * Ecm * Ecm / avgsvcalc.getMassLBSM() / avgsvcalc.getMassLBSM() - 1.
         << '\t' << g2_weff / sq(corr::part_hel_dof[corr::chi]) << '\t' << Weff_formula(input, Ecm) << '\t'
         << avgsvcalc.getWeff(Ecm) << '\n'; // Here adding the top mass?
  }

  fout.close();

  // Creating output file for sigmav
  fout.open(PATHPLOTS + "scalar_sigmav_all.out");
  if (!fout)
  {
    std::cerr << "Impossible to open " << PATHPLOTS << "scalar_sigmav_all.out\n";
    exit(1);
  }

  std::cout << "\nOpened " << PATHPLOTS << "scalar_sigmav_all.out\n";


  for (T = 1.0e-3; T <= 100. * input.getLightestBSMmass(); T *= 1.1)
  {
    fout << std::scientific << T / input.getLightestBSMmass() << '\t' << avgsvcalc.getAverageSigmav_coan_hightemp(T)
         << '\t' << avgsvcalc.getAverageSigmav_coan_lowtemp(T) << '\t';

    sigmavpred = 0.;
    for (auto count_part : {corr::e, corr::mu, corr::tau})
      sigmavpred += sigmav_pred(T, 1, count_part, input.g_l, input);
    for (auto count_part : {corr::d, corr::s, corr::b})
      sigmavpred += sigmav_pred(T, 3, count_part, input.g_d, input);
    for (auto count_part : {corr::u, corr::c, corr::t})
      sigmavpred += sigmav_pred(T, 3, count_part, input.g_u, input);

    fout << sigmavpred << '\n';
  }

  fout.close();


  // Low temperature behaviour: let us print out the Taylor coefficients for the <σv> expansion
  avgsvcalc.print_coefficients();

  sigmavpred = 0.;
  T = 1.;
  for (auto count_part : {corr::e, corr::mu, corr::tau})
    sigmavpred += sigmav_pred(T, 1, count_part, input.g_l, input);
  for (auto count_part : {corr::d, corr::s, corr::b})
    sigmavpred += sigmav_pred(T, 3, count_part, input.g_d, input);
  for (auto count_part : {corr::u, corr::c, corr::t})
    sigmavpred += sigmav_pred(T, 3, count_part, input.g_u, input);


  std::cout << "Taylor_pred_1 =  " << sigmavpred * input.getLightestBSMmass() << '\n'
            << "Ratio ourTaylor/FormulaTaylor  = "
            << avgsvcalc.get_TaylorCoeffSVT(1) / (sigmavpred * input.getLightestBSMmass()) << '\n';


  // Let us now compute the relic density
  BoltzmannSolver boltz(avgsvcalc);

  const double relic_density = boltz.relic_density();
  std::cout << "The relic density is: " << relic_density << '\n';

  return 0;
}
