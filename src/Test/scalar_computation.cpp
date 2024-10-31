#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include <fstream>
using namespace scalar2to2;
using namespace advmath;

int main(int argc, char ** argv)
{
 
    if(argc < 4)
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

    std::ifstream filein{argv[2]};
    if(!filein)
    {
        std::cout << "Impossible to open " << argv[2] << '\n';
        return 1;
    }

    std::ofstream fileout{argv[3]};
    if(!fileout)
    {
        std::cout << "Impossible to open " << argv[3] << '\n';
        return 1;
    }

    std::vector<double> val_mchi, val_ratio, val_g;
    unsigned int nlines=0;
    while(!filein.eof())
    {
        double temp;

        filein >> temp; val_mchi.emplace_back(temp);
        filein >> temp; val_ratio.emplace_back(temp);
        filein >> temp; val_g.emplace_back(temp);
        nlines++;
    }

    fileout << "#\tm_chi\tm_phi/m_chi\tg_chi\tOh2\txfo\n";
    for(size_t i=0; i < val_mchi.size(); i++)
    {
        std::cout << "\nComputing point " << i << '/' << nlines;
        input.m_chi=val_mchi[i];
        input.m_phi=input.m_chi*val_ratio[i];
        input.g_chi=val_g[i];
        input.refresh();
        boltz.changeInput(input);
        fileout << val_mchi[i] << '\t' << val_ratio[i] << '\t' << val_g[i] << '\t' <<
                   boltz.relic_density()<< '\t' << input.getLightestBSMmass()/boltz.Tfo << '\n';
    }
    
    
    return 0;
}
