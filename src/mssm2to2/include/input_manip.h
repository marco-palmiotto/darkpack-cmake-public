#ifndef INPUT_MANIP
#define INPUT_MANIP

#ifndef cplusplus__
#define cplusplus__
#endif

#include "mssm2to2.h"
#include "params_new.h"
#include "correspondance.h"

namespace mssm2to2::input_manip
{
  // The following is a MSSM specific function
  // It adapts the values of the trilinear couplings 
  // from the lha convention to the convention that works in MARTY
  
  void FixTrilinears(mssm2to2::param_t &input);
  
  // The following function replaces the CKM matrix
  // with the identity matrix

  void CKMIdentity(mssm2to2::param_t &input);
  
  // The following functions can be used to check if the CKM 
  // matrix in the input is unitary, and to replace it with the 
  // PDG values if needed

  bool CheckUnitarity(mssm2to2::param_t &input);

  bool ImposeUnitarity(mssm2to2::param_t &input);
}
#endif
