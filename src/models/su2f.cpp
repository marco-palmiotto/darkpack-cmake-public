#include <marty.h>
#include <marty/models/sm.h>
#include <marty/ufo/electriccharge.h>

using namespace mty;
using namespace mty::experimental;
using namespace mty::sm_input;
using namespace csl;

class SU2f_Model : public mty::SM_Model
{
  public:
  SU2f_Model(int model);

  void init();
  void initContent();
  void breakFlavor();
  void getToLowEnergyLagrangian();
  void initGauge();
  void initFermions();
  void initScalars();
  void initYukawas();
  void toFlavorCartan();
  void initFlavorCouplings();
  void expandFlavorVev();
  void adjust();
  // void addFlavorGaugeFixingTerms();

  private:
  int m_model;
};

SU2f_Model::SU2f_Model(int model) : SM_Model(false), m_model(model)
{
  this->init();
  this->refresh();
}

void SU2f_Model::init()
{
  this->initContent();
  this->breakFlavor();
  this->getToLowEnergyLagrangian();
}

void SU2f_Model::initContent()
{
  this->initGauge();
  this->initFermions();
  this->initScalars();
  this->initYukawas();
}

void SU2f_Model::breakFlavor()
{
  breakGaugeSymmetry("f");
  this->toFlavorCartan();
  this->initFlavorCouplings();
  this->expandFlavorVev();
}

void SU2f_Model::getToLowEnergyLagrangian()
{
  SM_Model::gaugeSymmetryBreaking();
  SM_Model::HiggsVEVExpansion();
  SM_Model::diagonalizeSMMassMatrices();
  SM_Model::replaceLeptonYukawa();
  SM_Model::replaceUpYukawa();
  SM_Model::replaceDownYukawa();
  breakFlavorSymmetry("su");
  SM_Model::flavorSymmetryBreaking();
  this->adjust();
}

void SU2f_Model::initGauge()
{
  addGaugedGroup(mty::group::Type::SU, "f", 2, csl::constant_s("g_f"));
  addFlavorGroup("su", 3);
  SM_Model::initGauge();
  renameParticle("A_f", "V");
}

void SU2f_Model::initFermions()
{
  SM_Model::initFermions();
  this->removeParticle("L");
  this->removeParticle("E_R");

  mty::Particle L = weylfermion_s("L", *this, mty::Chirality::Left);
  L->setGroupRep("L", 1);
  L->setGroupRep("Y", {-1, 2});
  L->setFundamentalFlavorRep("su");

  mty::Particle E = weylfermion_s("E_R", *this, mty::Chirality::Right);
  E->setGroupRep("Y", -1);
  E->setFundamentalFlavorRep("su");

  addParticles({L, E});

  csl::Expr m_cdm = csl::constant_s("m_cdm");
  mty::Particle Chi_dm_L = weylfermion_s("chi_dm_L;\\chi_{DM,L}", *this, mty::Chirality::Left);
  Chi_dm_L->setGroupRep("f", 1);
  Chi_dm_L->setMass(m_cdm);
  addParticle(Chi_dm_L);

  mty::Particle Chi_dm_R = weylfermion_s("chi_dm_R;\\chi_{DM,R}", *this, mty::Chirality::Right);
  Chi_dm_R->setGroupRep("f", 1);
  Chi_dm_R->setMass(m_cdm);
  addParticle(Chi_dm_R);
}

void SU2f_Model::initScalars()
{
  SM_Model::initHiggsPotential();

  mty::Particle phi = mty::scalarboson_s("Phi;\\Phi", *this);
  phi->setGroupRep("f", 1);
  addParticle(phi);

  csl::Index I = generateIndex("f", phi);

  csl::Expr phiSquared = csl::GetComplexConjugate(phi(I)) * phi(I);
  csl::Expr m_phi = csl::constant_s("m_phi;m_\\phi");
  csl::Expr v_phi = csl::constant_s("v_phi;v_\\phi");
  addScalarCoupling(m_phi);
  addScalarCoupling(v_phi);

  addLagrangianTerm(m_phi * m_phi / 2 * phiSquared);
  addLagrangianTerm(-csl::pow_s(m_phi * phiSquared / v_phi, 2) / 2);
}

void SU2f_Model::initYukawas()
{

  // Same as SM_Model

  auto* flavorSpace = mty::GetVectorSpace(*this, "SM_flavor");

  csl::Tensor Yu("Yu", {flavorSpace, flavorSpace});
  Yu->setComplexProperty(csl::ComplexProperty::Complex);
  csl::Tensor Yd("Yd", {flavorSpace, flavorSpace});
  Yd->setComplexProperty(csl::ComplexProperty::Complex);

  csl::Tensor eps = getVectorSpace("L", "Q")->getEpsilon();
  csl::Index I = flavorSpace->generateIndex();
  csl::Index J = flavorSpace->generateIndex();
  csl::Index a = generateIndex("C", "Q");
  csl::Index i = generateIndex("L", "Q");
  csl::Index j = generateIndex("L", "Q");
  csl::Index al = mty::DiracIndex();

  mty::Particle Q = getParticle("Q");
  mty::Particle U = getParticle("U_R");
  mty::Particle D = getParticle("D_R");
  mty::Particle H = getParticle("H");

  addLagrangianTerm(Yu({I, J}) * csl::GetComplexConjugate(H(i)) * eps({i, j}) * GetComplexConjugate(Q({I, a, j, al})) *
                        U({J, a, al}),
                    true);
  addLagrangianTerm(-Yd({I, J}) * H(i) * GetComplexConjugate(Q({I, a, i, al})) * D({J, a, al}), true);

  addTensorCoupling(Yu);
  addTensorCoupling(Yd);

  // Different from SM_Model (different flavor group)

  auto* gaugedFlavorSpace = mty::GetVectorSpace(*this, "su");

  csl::Tensor Ye("Ye", {gaugedFlavorSpace, gaugedFlavorSpace});
  Ye->setComplexProperty(csl::ComplexProperty::Complex);

  mty::Particle L = getParticle("L");
  mty::Particle E = getParticle("E_R");
  csl::Index II = gaugedFlavorSpace->generateIndex();
  csl::Index JJ = gaugedFlavorSpace->generateIndex();

  addLagrangianTerm(-Ye({II, JJ}) * H(i) * GetComplexConjugate(L({II, i, al})) * E({JJ, al}), true);

  addTensorCoupling(Ye);
}

void SU2f_Model::toFlavorCartan()
{
  mty::Particle V1 = GetParticle(*this, "V_1");
  mty::Particle V2 = GetParticle(*this, "V_2");
  mty::Particle V = V1->generateSimilar("V");
  V->setSelfConjugate(false);

  mty::Particle cV1 = getParticle("c_V_1");
  mty::Particle cV2 = getParticle("c_V_2");
  mty::Particle cp = V->getGhostBoson();
  cp->setName("c_Vp;c_p");
  mty::Particle cm = mty::ghostboson_s("c_Vm;c_m", V, true);
  V->setConjugatedGhostBoson(cm);

  csl::Index mu = mty::MinkowskiIndex();
  csl::Index nu = mty::MinkowskiIndex();
  csl::Expr V_p = V(+mu);
  csl::Expr V_m = csl::GetComplexConjugate(V(+mu));
  csl::Expr F_p = V({+mu, +nu});
  csl::Expr F_m = csl::GetComplexConjugate(V({+mu, +nu}));

  csl::Expr sq2 = csl::sqrt_s(2);
  replace(V1, (V_p + V_m) / sq2);
  replace(V2, CSL_I * (V_p - V_m) / sq2);
  replace(V1->getFieldStrength(), (F_p + F_m) / sq2);
  replace(V2->getFieldStrength(), CSL_I * (F_p - F_m) / sq2);
  replace(cV1, (static_cast<csl::Expr>(cp) + static_cast<csl::Expr>(cm)) / sq2);
  replace(cV2, CSL_I * (static_cast<csl::Expr>(cp) - static_cast<csl::Expr>(cm)) / sq2);
}

csl::Expr axisRotation(char axis, const csl::Expr& angle)
{
  csl::Expr rotMatrix{};
  switch (axis)
  {
  case 'x':
    rotMatrix = csl::matrix_s({{CSL_1, CSL_0, CSL_0},
                               {CSL_0, csl::cos_s(angle), -csl::sin_s(angle)},
                               {CSL_0, csl::sin_s(angle), csl::cos_s(angle)}});
    break;
  case 'y':
    rotMatrix = csl::matrix_s({{csl::cos_s(angle), CSL_0, -csl::sin_s(angle)},
                               {CSL_0, CSL_1, CSL_0},
                               {csl::sin_s(angle), CSL_0, csl::cos_s(angle)}});
    break;
  case 'z':
    rotMatrix = csl::matrix_s({{csl::cos_s(angle), -csl::sin_s(angle), CSL_0},
                               {csl::sin_s(angle), csl::cos_s(angle), CSL_0},
                               {CSL_0, CSL_0, CSL_1}});
    break;
  }
  return rotMatrix;
}

csl::Expr mtxProd(const csl::Expr& A, const csl::Expr& B)
{
  const int N = static_cast<int>(csl::Size(A));
  csl::Expr AB = csl::matrix_s(N, N);
  for (size_t i = 0; i != N; ++i)
  {
    for (size_t j = 0; j != N; ++j)
    {
      csl::Expr AB_ij = CSL_0;
      for (size_t k = 0; k != N; ++k)
      {
        AB_ij += A[i][k] * B[k][j];
      }
      AB[i][j] = AB_ij;
    }
  }
  return AB;
}

csl::Expr hermitian(const csl::Expr& A)
{
  const int N = static_cast<int>(csl::Size(A));
  csl::Expr A_T = csl::matrix_s(N, N);
  for (size_t i = 0; i != N; ++i)
  {
    for (size_t j = 0; j != N; ++j)
    {
      A_T[i][j] = csl::GetComplexConjugate(A[j][i]);
    }
  }
  return A_T;
}

csl::Expr mtxProd(const std::vector<csl::Expr>& factors)
{
  if (std::size(factors) == 1)
  {
    return factors[0];
  }
  else
  {
    return mtxProd(factors[0], mtxProd(std::vector<csl::Expr>(factors.begin() + 1, factors.end())));
  }
}

csl::Expr rotation(const csl::Expr& t_x, const csl::Expr& t_y, const csl::Expr& t_z)
{
  csl::Expr V_x = axisRotation('x', t_x);
  csl::Expr V_y = axisRotation('y', t_y);
  csl::Expr V_z = axisRotation('z', t_z);
  return mtxProd({V_z, V_x, V_y});
}

csl::Tensor Q(const std::string& name, const csl::Expr& T, const csl::Space* space, const csl::Expr& t_12,
              const csl::Expr& t_23, const csl::Expr& t_13)
{
  csl::Expr V_mtx = rotation(t_23, t_13, t_12);
  csl::Tensor Q = csl::tensor_s(name, {space, space}, mtxProd({hermitian(V_mtx), T, V_mtx}));
  return Q;
}

void SU2f_Model::initFlavorCouplings()
{
  auto flavorSpace = mty::GetVectorSpace(*this, "su");
  csl::Index I = flavorSpace->generateIndex();
  csl::Index J = flavorSpace->generateIndex();
  csl::Index i = generateIndex("L", "L");
  csl::Index mu = mty::MinkowskiIndex();
  csl::Index al = mty::DiracIndex();
  csl::Index be = mty::DiracIndex();
  csl::Tensor gamma = mty::DiracGamma();
  csl::Expr g = getScalarCoupling("g_f");

  mty::Particle V = getParticle("V");
  mty::Particle V_3 = getParticle("V_3");
  mty::Particle L = getParticle("L");
  mty::Particle E = getParticle("E_R");

  csl::Expr T_p, T_3;

  switch (this->m_model)
  {
  case 12:
    T_p = csl::matrix_s({{CSL_0, CSL_1, CSL_0}, {CSL_0, CSL_0, CSL_0}, {CSL_0, CSL_0, CSL_0}}) / csl::sqrt_s(2);

    T_3 = csl::matrix_s({{CSL_HALF, CSL_0, CSL_0}, {CSL_0, CSL_M_HALF, CSL_0}, {CSL_0, CSL_0, CSL_0}});
    break;
  case 13:
    T_p = csl::matrix_s({{CSL_0, CSL_0, CSL_1}, {CSL_0, CSL_0, CSL_0}, {CSL_0, CSL_0, CSL_0}}) / csl::sqrt_s(2);

    T_3 = csl::matrix_s({{CSL_HALF, CSL_0, CSL_0}, {CSL_0, CSL_0, CSL_0}, {CSL_0, CSL_0, CSL_M_HALF}});
    break;
  case 23:
    T_p = csl::matrix_s({{CSL_0, CSL_0, CSL_0}, {CSL_0, CSL_0, CSL_1}, {CSL_0, CSL_0, CSL_0}}) / csl::sqrt_s(2);

    T_3 = csl::matrix_s({{CSL_0, CSL_0, CSL_0}, {CSL_0, CSL_HALF, CSL_0}, {CSL_0, CSL_0, CSL_M_HALF}});
    break;
  }

  csl::Expr tl_12 = csl::constant_s("theta_12_L;\\theta_{12}^L");
  csl::Expr tr_12 = csl::constant_s("theta_12_R;\\theta_{12}^R");
  csl::Expr tl_13 = csl::constant_s("theta_13_L;\\theta_{13}^L");
  csl::Expr tr_13 = csl::constant_s("theta_13_R;\\theta_{13}^R");
  csl::Expr tl_23 = csl::constant_s("theta_23_L;\\theta_{23}^L");
  csl::Expr tr_23 = csl::constant_s("theta_23_R;\\theta_{23}^R");
  csl::Tensor Q_L = Q("Q_L", T_p, flavorSpace, tl_12, tl_23, tl_13);
  csl::Tensor Q_3_L = Q("Q_3_L", T_3, flavorSpace, tl_12, tl_23, tl_13);
  csl::Tensor Q_R = Q("Q_R", T_p, flavorSpace, tr_12, tr_23, tr_13);
  csl::Tensor Q_3_R = Q("Q_3_R", T_3, flavorSpace, tr_12, tr_23, tr_13);

  // Coupling between Vpm and L
  addLagrangianTerm(
      -g * V(+mu) * csl::GetComplexConjugate(L({I, i, al})) * Q_L({I, J}) * gamma({mu, al, be}) * L({J, i, be}), true);

  // Coupling between V3 and L
  addLagrangianTerm(-g * V_3(+mu) * csl::GetComplexConjugate(L({I, i, al})) * Q_3_L({I, J}) * gamma({mu, al, be}) *
                    L({J, i, be}));

  // Coupling between Vpm and E
  addLagrangianTerm(-g * V(+mu) * csl::GetComplexConjugate(E({I, al})) * Q_R({I, J}) * gamma({mu, al, be}) * E({J, be}),
                    true);

  // Coupling between V3 and E
  addLagrangianTerm(-g * V_3(+mu) * csl::GetComplexConjugate(E({I, al})) * Q_3_R({I, J}) * gamma({mu, al, be}) *
                    E({J, be}));

  addTensorCoupling(Q_L);
  addTensorCoupling(Q_3_L);
  addTensorCoupling(Q_R);
  addTensorCoupling(Q_3_R);
}

void SU2f_Model::expandFlavorVev()
{
  mty::Particle phi_1 = getParticle("Phi_1");
  mty::Particle phi_2 = getParticle("Phi_2");
  csl::Expr v = getScalarCoupling("v_phi");

  mty::Particle phi = mty::scalarboson_s("phi;\\phi", *this); // Flavon
  mty::Particle chi_0 = mty::scalarboson_s("chi_0;\\chi^0", *this);
  mty::Particle chi_p = mty::scalarboson_s("chi_p;\\chi^p", *this);
  phi->setSelfConjugate(true);
  chi_0->setSelfConjugate(true);

  replace(phi_1, chi_p());
  replace(phi_2, (v + phi() + CSL_I * chi_0()) / csl::sqrt_s(2));
}

void SU2f_Model::adjust()
{

  csl::Expr M_V = csl::constant_s("M_V"); // V mass
  csl::Expr g_f = getScalarCoupling("g_f");
  csl::Expr v_phi = getScalarCoupling("v_phi");

  getParticle("V")->setMass(M_V);
  replace(getParticle("V")->getMass(), M_V);

  getParticle("V_3")->setMass(M_V);
  replace(getParticle("V_3")->getMass(), M_V);
  replace(v_phi, 2 * M_V / g_f);

  promoteToGoldstone("chi_0", "V_3");
  promoteToGoldstone("chi_p", "V");

  diracFermionEmbedding("chi_dm_L_1", "chi_dm_R_1");
  diracFermionEmbedding("chi_dm_L_2", "chi_dm_R_2");

  SM_Model::adjust();
}

#include "marty/ufo/electriccharge.h"

#define LIBNAME         "su2f2to2" // Name of the output library
#define MAXNUMPROCESSES 5000       // Max number of processes to be stored in a Processes variable

#define CORRESPONDANCE // Defining this variable the creation of correspondance.h is enabled
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


using Process = std::vector<Insertion>;
using Processes = std::vector<Process>;


bool setShow = false;

#include "src/processName.cpp"

#include "src/generateCorrespondance.cpp"

#include "src/checkCPsymmetric.cpp"

#include "src/computeAndAddToLibFromList.cpp"


int main()
{
  SU2f_Model su2f_model(23);

  std::cout << su2f_model;


  std::cout << "\nSuppressing verbose output during the calculation of the amplitudes\n";
  mty::option::verboseAmplitude = false;
  std::cout << "\nChoosing not to call gafed to show Feynman diagrams during the computation\n";
  std::cout << "of an amplitude. They will be printed in a directory\n";
  setShow = false;
  std::cout << "Creating the library\n";
  mty::Library lib(LIBNAME);

  std::cout << "Adding the flag for using posix threads in the C++ files in the library\n";
  lib.addLibrary("-lpthread");

  std::cout << "\nListing all the Physical particles\n";
  std::vector<Particle> part_0 = su2f_model.getPhysicalParticles([&](Particle p) { return p->isPhysical(); });
  std::vector<Particle> part;
  std::cout << "Printing Physical particles before filtering\n";
  for (size_t i = 0; i != part_0.size(); i++)
  {
    std::cout << part_0[i]->getName() << '\n';
    if (!(IsOfType<GhostBoson>(part_0[i]) || IsOfType<GoldstoneBoson>(part_0[i])))
      part.push_back(part_0[i]);
  }

  std::cout << "Printing Physical particles after filtering\n";
  for (auto elem : part)
    std::cout << elem->getName() << '\n';

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
  lib.generateSpectrum(su2f_model);
  std::cout << "Spectrum Generated" << std::endl;

  // The following replacement is necessary to have no symbol dependencies in the electric charge
  su2f_model.replace((1 + pow_s(tan_s(theta_W), 2)), 1 / pow_s(cos_s(theta_W), 2));

  std::cout << "Separating SM and BSM particles" << std::endl;
  std::vector<Particle> psm;
  std::vector<Particle> pbsm({su2f_model.getParticle("chi_dm_1"), su2f_model.getParticle("chi_dm_2")});
  for (size_t i = 0; i != part.size(); i++)
  {
    bool istoadd = true;
    for (size_t j = 0; j != pbsm.size() && istoadd; j++)
    {
      istoadd = part[i]->getName() != pbsm[j]->getName();
    }
    if (istoadd)
      psm.push_back(part[i]);
  }

  if (pbsm.size() > 0)
    generateCorrespondance(psm, pbsm);

  std::cout << "\nSetting widths symbolically\n";
  Expr Gamma_h = constant_s("Gamma_h");
  Expr Gamma_W = constant_s("Gamma_W");
  Expr Gamma_Z = constant_s("Gamma_Z");
  Expr Gamma_top = constant_s("Gamma_top");
  Expr Gamma_V = constant_s("Gamma_V");
  Expr Gamma_V3 = constant_s("Gamma_V_3");

  std::cout << "Constants have been defined\n";

  su2f_model.getParticle("h")->setWidth(Gamma_h);
  su2f_model.getParticle("W")->setWidth(Gamma_W);
  su2f_model.getParticle("Z")->setWidth(Gamma_Z);
  su2f_model.getParticle("t")->setWidth(Gamma_top);
  su2f_model.getParticle("V")->setWidth(Gamma_V);
  su2f_model.getParticle("V_3")->setWidth(Gamma_V3);


  std::cout << "\nGetting electric charges: (Feynman rules are needed at this purpose)\n";
  std::vector<double> electric_charges;
  for (const auto& particle : part)
  {
    // Create the necessary data for the module to extract the charge
    ufo::ElectricChargeMessenger elec{
        particle.get(),                    // Particle from which we need the charge
        su2f_model.getParticle("A").get(), // Vector boson carrying the charge
        su2f_model,                        // Model
        ufo::SMData{}                      // Model conventions
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


  // Our list is formed by all the processes of the kind Chibar_{1,2,3} Chi_{1,2,3} -> Fbar F
  // where F is a SM lepton

  std::cout << "Building the list of SM leptons\n";

  std::array<Particle, 2> bsm_fields({su2f_model.getParticle("chi_dm_1"), su2f_model.getParticle("chi_dm_2")});

  for (size_t i1 = 0; i1 < bsm_fields.size(); i1++)
  {
    for (mty::Insertion state1 :
         std::array<mty::Insertion, 2>({Incoming(bsm_fields[i1]), Incoming(AntiPart(bsm_fields[i1]))}))
    {
      for (size_t i2 = i1; i2 < bsm_fields.size(); i2++)
      {
        for (mty::Insertion state2 :
             std::array<mty::Insertion, 2>({Incoming(bsm_fields[i2]), Incoming(AntiPart(bsm_fields[i2]))}))
        {
          for (size_t j1 = 0; j1 < psm.size(); j1++)
          {
            for (mty::Insertion state3 :
                 std::array<mty::Insertion, 2>({Outgoing(psm[j1]), Outgoing(AntiPart(psm[j1]))}))
            {
              for (size_t j2 = j1; j2 < psm.size(); j2++)
              {
                for (mty::Insertion state4 :
                     std::array<mty::Insertion, 2>({Outgoing(psm[j2]), Outgoing(AntiPart(psm[j2]))}))
                {
                  Process2to2ToCompute proc_c;

                  proc_c.process = {(state1), (state2), (state3), (state4)};
                  proc_c.order = mty::Order::TreeLevel;
                  proc_c.leading_order = false;
                  proc_c.Wgauge = mty::gauge::Type::Feynman;

                  list_of_processes.push_back(proc_c);
                }
              }
            }
          }
        }
      }
    }
  }


  computeAndAddToLibFromList(su2f_model, lib, list_of_processes);


  // Computation of the widths at the three TreeLevel
  Expr width_h = su2f_model.computeWidth(Order::TreeLevel, "h");
  Expr width_W = su2f_model.computeWidth(Order::TreeLevel, "W");
  Expr width_Z = su2f_model.computeWidth(Order::TreeLevel, "Z");
  Expr width_t = su2f_model.computeWidth(Order::TreeLevel, "t");
  Expr width_V = su2f_model.computeWidth(Order::TreeLevel, "V");
  Expr width_V3 = su2f_model.computeWidth(Order::TreeLevel, "V_3");


  std::vector<Expr> expr_to_expand = {width_h, width_W, width_Z, width_t, width_V, width_V3};

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
  lib.addFunction("width_V", width_V);
  lib.addFunction("width_V_3", width_V3);
  lib.addDefaultParameter("g_s", false);

#ifdef SQUAREDAMP
  std::cout << "Writing the library on disk\n";
  lib.print();
#endif

  return 0;
}
