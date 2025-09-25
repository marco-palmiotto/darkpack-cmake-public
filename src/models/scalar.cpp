#include "marty/models/sm.h"
#include "marty/sgl/sgl.h"
#include <filesystem>
#include <functional>
#include <locale>
#include <marty.h>

#include "marty/ufo/electriccharge.h"

#define LIBNAME         "scalar2to2" // Name of the output library
#define MAXNUMPROCESSES 5000         // Max number of processes to be stored in a Processes variable

#define CORRESPONDANCE // Defining this variable the creation of correspondance.hpp is enabled
#define SQUAREDAMP     // Defining this variable squared amplitudes are calculated and added to the
//                        library. Disabling this behaviour can be useful if the user just wants
//                        to generate the auxiliary files for the library

#define FULL_PARAM
#define DISPLAYSQUAREDAMPLITUDESEXPRESSIONS
// #define DEBUG

template <typename T> inline void printdebug(const T& s1)
{
#ifdef DEBUG
  std::cout << s1 << std::endl;
#endif
  return;
};

using namespace mty;
using namespace mty::experimental;
using namespace mty::sm_input;
using namespace csl;

using Process = std::vector<Insertion>;
using Processes = std::vector<Process>;

bool setShow = false;

#include "src/processName.cpp"

#include "src/generateCorrespondance.cpp"

#include "src/checkCPsymmetric.cpp"

#include "src/computeAndAddToLibFromList.cpp"

// Defining the scalar model as a class
#define BASE_MODEL SM_Model
class ScalarModel : public BASE_MODEL
{
  public:
  /*  It's better to define an inline constrcutor
   *  to tell it not to intialise the base model
   *  and handle here all the details of the initialisation
   */
  ScalarModel(bool init = true) : BASE_MODEL(false)
  {
    printdebug("Initialising the scalar model");
    initScalarModel();
  }

  private:
  void addDMcandidates();

  void getToLowEnergyLagrangian();

  void initScalarModel();
};

void ScalarModel::addDMcandidates()
{
  printdebug("Initialising Phi and DMchi");

  Particle Phi = scalarboson_s("phi;\\phi", *this);
  Phi->setSelfConjugate(true);

  Particle Chi = diracfermion_s("chi;\\chi", *this);

  Expr m_chi = constant_s("m_chi");
  Expr m_phi = constant_s("m_phi");

  Chi->setMass(m_chi);
  Phi->setMass(m_phi);

  addParticles({Phi, Chi});

#ifdef DEBUG
  std::cout << *this;
#endif

  printdebug("Defining constants to be put in the potential");

  Expr g_chi = constant_s("g_chi");
  Expr g_u = constant_s("g_u");
  Expr g_d = constant_s("g_d");
  Expr g_l = constant_s("g_l");

  // Getting expressions to write the potential

  Tensor Y_u = GetYukawa(*this, "Yu");
  Tensor Y_d = GetYukawa(*this, "Yd");
  Tensor Y_e = GetYukawa(*this, "Ye");

  Particle uL = getParticle("U_L");
  Particle uR = getParticle("U_R");
  Particle dL = getParticle("D_L");
  Particle dR = getParticle("D_R");
  Particle eL = getParticle("E_L");
  Particle eR = getParticle("E_R");

  // Defining indexes and their spaces
  const csl::Space* flavorSpace = GetSpace(Y_e);
  csl::Index f1 = GetIndex(flavorSpace);
  csl::Index f2 = GetIndex(flavorSpace);

  Index a = DiracIndex();
  Index A = GaugeIndex(*this, "C", uR); // Color index

  Expr VPhiChi = -g_chi * GetComplexConjugate(Chi(a)) * Chi(a) * Phi;
  addLagrangianTerm(VPhiChi);

  Expr VPhiMatter = -Phi / sqrt_s(2) *
                    (g_u * Y_u({f1, f2}) * GetComplexConjugate(uL({f1, A, a})) * uR({f2, A, a}) +
                     g_d * Y_d({f1, f2}) * GetComplexConjugate(dL({f1, A, a})) * dR({f2, A, a}) +
                     g_l * Y_e({f1, f2}) * GetComplexConjugate(eL({f1, a})) * eR({f2, a}));
  addLagrangianTerm(VPhiMatter, true);
  printdebug("Added SM-DM interaction terms");
};

void ScalarModel::getToLowEnergyLagrangian()
{
  printdebug("Calling gaugeSymmetryBreaking");
  gaugeSymmetryBreaking();
  //     std::cout << *this;

  printdebug("Calling HiggsVEVExpansion");
  HiggsVEVExpansion();
  //     std::cout << *this;

  printdebug("Calling diagonalizeMassMatrices");
  diagonalizeSMMassMatrices();
  //     std::cout << *this;

  printdebug("Calling addDMcandidates");
  // The fields have no EW charge, so they can be added
  // after the EW-SSB
  // Adding the potential terms before breaking the flavor
  // symmetry is convenient , since they're
  // symmetric under flavor
  // We do both things with the same method now
  addDMcandidates();
  // std::cout << *this;

  printdebug("Calling replaceLeptonYukawa");
  replaceLeptonYukawa();
  printdebug("Calling replaceUPYukawa");
  replaceUpYukawa();
  printdebug("Calling replaceDownYukawa");
  replaceDownYukawa();
  // std::cout << *this;
  printdebug("Calling flavorSymmetryBreaking");
  flavorSymmetryBreaking();
  // std::cout << *this;
  printdebug("Calling adjust");
  adjust();
}

void ScalarModel::initScalarModel()
{
BASE_MODEL:
  initContent();
  ScalarModel::getToLowEnergyLagrangian();

  refresh();
};
#undef BASE_MODEL


int main()
{
  std::cout << "Defining the model\n";
  ScalarModel scalar;
  std::cout << scalar;

  std::cout << "\nSuppressing verbose output during the calculation of the amplitudes\n";
  mty::option::verboseAmplitude = false;
  std::cout << "\nChoosing not to call gafed to show Feynman diagrams during the computation\n";
  std::cout << "of an amplitude. They will be printed in a directory\n";
  setShow = false;
  std::cout << "Creating the library\n";
  mty::Library lib(LIBNAME);

  std::cout << "\nListing all the Physical particles\n";
  std::vector<Particle> part_0 = scalar.getPhysicalParticles([&](Particle p) { return p->isPhysical(); });
  std::vector<Particle> part;
  for (size_t i = 0; i != part_0.size(); i++)
  {
    if (!(IsOfType<GhostBoson>(part_0[i]) || IsOfType<GoldstoneBoson>(part_0[i])))
      part.push_back(part_0[i]);
  }

  std::cout << "Re-fixing mass names:\n";
  for (size_t i = 0; i != part.size(); i++)
  {
    std::string tempstring;
    tempstring.clear();
    if (!part[i]->getMass()->getName().empty())
    {
      part[i]->getMass()->setName("m_" + part[i]->getName());
      tempstring = part[i]->getMass()->getName();
      std::cout << part[i]->getName() << " " << tempstring << std::endl;
    }
  }

  std::cout << "Generating Spectrum..." << std::endl;
  lib.generateSpectrum(scalar);
  std::cout << "Spectrum Generated" << std::endl;

  // The following replacement is necessary to have no symbol dependencies in the electric charge
  scalar.replace((1 + pow_s(tan_s(theta_W), 2)), 1 / pow_s(cos_s(theta_W), 2));

  std::cout << "Separating SM and BSM particles" << std::endl;
  std::vector<Particle> psm, pbsm;

  for (size_t i = 0; i != part.size(); i++)
  {
    if (part[i]->getName() == "h" || part[i]->getName() == "W" || part[i]->getName() == "Z" ||
        part[i]->getName() == "G" || part[i]->getName() == "A" || part[i]->getName() == "d" ||
        part[i]->getName() == "u" || part[i]->getName() == "s" || part[i]->getName() == "c" ||
        part[i]->getName() == "b" || part[i]->getName() == "t" || part[i]->getName() == "e" ||
        part[i]->getName() == "mu" || part[i]->getName() == "tau" || part[i]->getName() == "nu_e" ||
        part[i]->getName() == "nu_mu" || part[i]->getName() == "nu_tau" || part[i]->getName() == "phi")
    {
      psm.push_back(part[i]);
    }
    else
    {
      pbsm.push_back(part[i]);
    }
  }

  if (pbsm.size() > 0)
    generateCorrespondance(psm, pbsm);

  std::cout << "\nSetting widths symbolically\n";
  Expr Gamma_h = constant_s("Gamma_h");
  Expr Gamma_W = constant_s("Gamma_W");
  Expr Gamma_Z = constant_s("Gamma_Z");
  Expr Gamma_top = constant_s("Gamma_top");
  Expr Gamma_phi = constant_s("Gamma_phi");

  std::cout << "Constants have been defined\n";

  scalar.getParticle("h")->setWidth(Gamma_h);
  scalar.getParticle("W")->setWidth(Gamma_W);
  scalar.getParticle("Z")->setWidth(Gamma_Z);
  scalar.getParticle("t")->setWidth(Gamma_top);
  scalar.getParticle("phi")->setWidth(Gamma_phi);

  std::cout << "\nGetting electric charges: (Feynman rules are needed at this purpose)\n";
  std::vector<double> electric_charges;
  for (const auto& particle : part)
  {
    // Create the necessary data for the module to extract the charge
    ufo::ElectricChargeMessenger elec{
        particle.get(),                // Particle from which we need the charge
        scalar.getParticle("A").get(), // Vector boson carrying the charge
        scalar,                        // Model
        ufo::SMData{}                  // Model conventions
    };

    // The charge is read in the vertices
    ufo::Charge electronCharge = ufo::getElectricCharge(elec);
    // We convert it to a symbolic expression and display it
    csl::Expr symbCharge = csl::intfraction_s(electronCharge.num, electronCharge.denom);
    double float_charge = electronCharge.num * 1. / electronCharge.denom;
    std::cout << particle->getName() << " has charge: " << symbCharge << " = " << float_charge << std::endl;
    electric_charges.push_back(float_charge);
  }

  std::cout << "\nGetting degrees of freedom:\n";
  std::vector<int> dofs;
  for (const auto& particle : part)
  {
    std::cout << particle->getName() << " has dof: " << particle->getNDegreesOfFreedom() << std::endl;
  }

  std::vector<Process2to2ToCompute> list_of_processes;
  std::vector<Process2to2ToCompute> list_of_processes_1to2;

  /*
  Our list is formed by all the processes of the kind Chibar Chi -> Fbar F
  where F is a SM fermion, so we need to build the list of SM fermions.
  NOTE: we could already use psm, but all the other amplitudes will
        be null, and therefore not be present in the library.
  */
  std::cout << "Building the list of SM fermions\n";

  std::vector<Particle> smfermions;

  for (size_t i = 0; i != part.size(); i++)
  {
    if (part[i]->getName() == "d" || part[i]->getName() == "u" || part[i]->getName() == "s" ||
        part[i]->getName() == "c" || part[i]->getName() == "b" || part[i]->getName() == "t" ||
        part[i]->getName() == "e" || part[i]->getName() == "mu" || part[i]->getName() == "tau" ||
        part[i]->getName() == "nu_e" || part[i]->getName() == "nu_mu" || part[i]->getName() == "nu_tau")
    {
      smfermions.push_back(part[i]);
    }
  }

  for (auto field : smfermions)
  {
    Process2to2ToCompute proc_c;

    proc_c.process = {Incoming(AntiPart(scalar.getParticle("chi"))), Incoming(scalar.getParticle("chi")),
                      Outgoing(AntiPart(field)), Outgoing(field)};
    proc_c.order = mty::Order::TreeLevel;
    proc_c.leading_order = false;
    proc_c.Wgauge = mty::gauge::Type::Feynman;

    list_of_processes.push_back(proc_c);
  }

  for (auto field : part)
  {
    Process2to2ToCompute proc_c;

    proc_c.process = {Incoming(scalar.getParticle("h")), Outgoing(AntiPart(field)), Outgoing(field)};
    proc_c.order = mty::Order::TreeLevel;
    proc_c.leading_order = false;
    proc_c.Wgauge = mty::gauge::Type::Feynman;

    list_of_processes_1to2.push_back(proc_c);
  }

  computeAndAddToLibFromList(scalar, lib, list_of_processes, list_of_processes_1to2, "smBsm.hpp");

  // Computation of the widths at the three TreeLevel
  Expr width_h = scalar.computeWidth(Order::TreeLevel, "h");
  Expr width_W = scalar.computeWidth(Order::TreeLevel, "W");
  Expr width_Z = scalar.computeWidth(Order::TreeLevel, "Z");
  Expr width_t = scalar.computeWidth(Order::TreeLevel, "t");
  Expr width_phi = scalar.computeWidth(Order::TreeLevel, "phi");


  std::cout << "\n\n" << Evaluated(width_phi, eval::abbreviation) << std::endl << std::endl;

  /*
  We construct a vector of expressions to expand
  NOTE: this prevents hardcoding some quantities in the library
  */
  std::vector<Expr> expr_to_expand = {width_h, width_W, width_Z, width_t, width_phi};

  for (Expr expression : expr_to_expand)
  {
    auto temp = csl::Evaluated(expression, csl::eval::abbreviation);
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
  }

  lib.addFunction("width_h", width_h);
  lib.addFunction("width_W", width_W);
  lib.addFunction("width_Z", width_Z);
  lib.addFunction("width_t", width_t);
  lib.addFunction("width_phi", width_phi);

  lib.addDefaultParameter("g_s", false);

#ifdef SQUAREDAMP
  std::cout << "Writing the library on disk\n";
  lib.print();
#endif
  return 0;
}
