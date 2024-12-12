#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp"
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"

using namespace scalar2to2;
using namespace advmath;

std::string thisfolder = "/mnt/Data/projects/darkpack-cmake";

int main(int argc, char** argv)
{

  if (argc < 4)
  {
    std::cout << "This function requires 3 arguments:\n"
              << " - The name of the file with the input data\n"
              << " - The name of the file with the input to vary (i.e. m_chi, m_phi/m_chi, g_chi)\n"
              << " - The name of the file where to write the output\n";
    return 1;
  }
  std::cout << std::setprecision(5) << std::scientific;

  std::cout << "Reading inputs from " << argv[1] << std::endl;

  struct Param_t input(argv[1]);
  BoltzmannSolver boltz(input);
  // boltz.solver=2;

  std::ifstream filein{argv[2]};
  if (!filein)
  {
    std::cout << "Impossible to open " << argv[2] << '\n';
    return 1;
  }
  std::cout << "\nOpened " << argv[2];

  std::ofstream fileout{argv[3]};
  if (!fileout)
  {
    std::cout << "Impossible to open " << argv[3] << '\n';
    return 1;
  }
  std::cout << "\nOpened " << argv[3];

  std::vector<double> val_mchi, val_ratio, val_g;
  unsigned int nlines = 0;
  while (!filein.eof())
  {
    double temp;
    std::cout << "\nReading line " << nlines;
    filein >> temp;
    val_mchi.emplace_back(temp);
    filein >> temp;
    val_ratio.emplace_back(temp);
    filein >> temp;
    val_g.emplace_back(temp);
    nlines++;
  }

  for (size_t i = 0; i < val_mchi.size(); i++)
  {
    std::cout << "\nPoint " << i << '/' << nlines;
    std::cout << "\nm_chi=" << val_mchi[i] << " ratio=" << val_ratio[i] << " g=" << val_g[i];
  }

  fileout << "#\tm_chi\tm_phi/m_chi\tg_chi\tOh2\txfo\t<sv>filename\n";
  for (size_t i = 0; i < val_mchi.size(); i++)
  {
    std::cout << "\nComputing point " << i << '/' << nlines;
    std::cout << "\nm_chi=" << val_mchi[i] << " ratio=" << val_ratio[i] << " g=" << val_g[i];

    std::string filename = thisfolder + "/out/scans/sigmav/scalar_sv_" + std::to_string(i + 1) + ".out";
    std::ofstream sigmavfile(filename.c_str());
    if (!sigmavfile)
    {
      std::cerr << "Impossible to open " << filename << std::endl;
      return 1;
    }

    input.m_chi = val_mchi[i];
    input.m_phi = input.m_chi * val_ratio[i];
    input.g_chi = 0.3;
    input.refresh();
    boltz.changeInput(input);
    boltz.full_comput = 1;
    boltz.solver = 1;
    const real_t relicd = boltz.relic_density(), xfo = input.getLightestBSMmass() / boltz.Tfo;
    // Printing <sv> on file
    const real_t xmin = 1.0e-4;
    const real_t xmax = 1.0e+3;
    const int npoints_sv = 500;
    const real_t factor_sv = std::pow(xmax / xmin, 1. / npoints_sv);
    for (int j = 0; j <= npoints_sv; j++)
    {
      const real_t x = xmin * std::pow(factor_sv, j);
      sigmavfile << std::setprecision(5) << x << '\t';
      const real_t sigmav = boltz.getAverageSigmav_coan(input.getLightestBSMmass() / x);
      sigmavfile << sigmav << '\n';
    }
    std::cout << val_mchi[i] << '\t' << val_ratio[i] << '\t' << val_g[i] << '\t' << relicd << '\t' << xfo << '\n';
    fileout << val_mchi[i] << '\t' << val_ratio[i] << '\t' << val_g[i] << '\t' << relicd << '\t' << xfo << '\t'
            << filename << '\n';
  }


  return 0;
}
