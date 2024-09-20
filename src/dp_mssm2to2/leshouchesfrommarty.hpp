#ifndef MODEL_MSSM
#define MODEL_MSSM
#endif 

#ifndef LESHOUCESFROMMARTY_INCLUDED
#define LESHOUCESFROMMARTY_INCLUDED


#include "config.hpp"
#include "process.hpp"
#include "marty/lha/lha.h"

#include "input_manip.hpp"
#include "RunningSM.hpp"

namespace mssm2to2::readmodule
{

inline void computeThetaWandMw(const double alpha_em, const double Gf, const double Mz, double &thetaW, double &Mw)
{
    thetaW=0.5*std::asin(std::sqrt(4*M_PI*alpha_em/(std::sqrt(2)*Gf*Mz*Mz)));
    Mw = Mz*std::cos(thetaW);
    return;
}


void AssignSMparams(Param_t &param, const std::string &namefile);

void ReadBSMparams(Param_t &param, const std::string &namefile);

void ReadWidths(Param_t &param, const std::string &namefile);

Param_t ReadLHA(const std::string &name);

void createwidthsfile(const Param_t &param, const std::string &namefile);


}
#endif