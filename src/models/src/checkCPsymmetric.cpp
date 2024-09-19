// This function takes as input a model m and a process p of the kind 
// A, B -> C, D in such model.
// It returns true if the process anti A, anti B -> anti C, anti D 
// is possible, otherise it returns false
bool checkCPsymmetric(mty::Model &m, Process &p)
{
  bool cptrue=false;
  bool exactSameIncoming=true, exactSameOutgoing=true;
  bool sameIncoming=true, sameOutgoing=true;
  Process cpsymmetric;
  
  for(auto j : {0,1})
  {
    if(!(p[j].getField()->isSelfConjugate()))
    {
      exactSameIncoming=false;
    }
  }
  
  for(auto j : {2,3})
  {
    if(!(p[j].getField()->isSelfConjugate()))
    {
      exactSameOutgoing=false;
    }
  }
  
  cpsymmetric = AntiPart(p);
  
  // If the state is the exact one, then we have the same external state
  // Otherwise, we have to check permutations
  if(!exactSameIncoming && !(cpsymmetric[0] == p[1] && cpsymmetric[1] == p[0])) 
  {
    sameIncoming = false;
  }
  
  if(!exactSameOutgoing && !(cpsymmetric[2] == p[3] && cpsymmetric[3] == p[2])) 
  {
    sameOutgoing = false;
  }
  
  if( !sameIncoming || !sameOutgoing)
  {
    auto cpampl= m.computeAmplitude( Order::TreeLevel, cpsymmetric) ;
    cptrue = !cpampl.empty();
  }
  return cptrue;
}
