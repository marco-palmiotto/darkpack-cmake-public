#pragma once

#include "functions.cpp"
#include <marty.h>
#include <marty/models/sm.h>

class SU2f_Model : public virtual mty::SM_Model
{
  public:
  SU2f_Model(const int model, const bool refresh = true);

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

  protected:
  int m_model; // Variable to determine which families are coupled
};

SU2f_Model::SU2f_Model(const int model, const bool refresh) : SM_Model(false), m_model(model)
{
  this->init();
  this->adjust();
  if (refresh)
  {
    this->refresh();
  }
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