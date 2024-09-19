/*--------------------------------------------------------------------*/
/*                             AlterBBN                               */
/*--------------------------------------------------------------------*/

// #define pi    		3.1415926535897932385
// #define zeta3 		1.2020569031595942855
// #define Gn    		6.67428e-8  /* in cm^3.g^-1.s^-2 */
// #define K_to_eV     8.617330637338339e-05 /* conversion factor T(10**9 K) * K_to_GeV = T(GeV) or T(K) * K_to_eV = T(eV) */
// #define m_e    		510.9989461e-6 /* electron mass in GeV */
// #define M_u         1.660538921e-24 /* unit atomic mass M_u in g */
// #define g_to_GeV    5.60958884538932e+23 /* conversion factor M(g) * g_to_GeV = M(GeV) */
// #define kg_to_GeV   5.60958884538932e+26 /* conversion factor M(kg) * kg_to_GeV = M(GeV) */
// #define sigma_SB    0.16449340668482282 /* Stefan-Boltzmann constant = pi^2/60 */
// #define m_to_GeV    5.067730582705779e+15 /* conversion factor L(m) * m_to_GeV = L(GeV^-1) */
// #define cm_to_GeV   5.067730582705779e+13 /* conversion factor L(cm) * cm_to_GeV = L(GeV^-1) */
// #define s_to_GeV    1.5192674078711377e+24 /* conversion factor t(s) * s_to_GeV = t(GeV^-1) */
// #define G           6.708609142443796e-39 /* Gn*pow(m_to_GeV,3.)*pow(g_to_GeV,-1.)*pow(s_to_GeV,-2.) Newton constant in GeV^-2 */
// #define Mplanck     1.2209102930946623e+19 /* in GeV, more precise definition than before */
// #define DMpn        0.0012934 /* mass difference between neutron and proton in GeV */
// #define zeta        1.6103162253325862 /* 3*k_B/(2*c^2*Mu) in GeV^-1 */
// #define k_B			8.617330e-5		/* Boltzmann's constant in GeV/GK */
// #define alphaem  		0.007297353 /* fine-structure constant */

#include "params_new.hpp"
using namespace __SPEC_LIB_NAME__::mathematicalconstants;
using namespace __SPEC_LIB_NAME__::pdg2016Value;
