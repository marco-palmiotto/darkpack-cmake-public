#include "leshouches_patch.hpp"

// static inline double creal(Complex_C_t &z)
// {
//   return (static_cast<std::complex<double>>(z)).real();
// }

// static inline double cimag(Complex_C_t &z)
// {
//   return (static_cast<std::complex<double>>(z)).imag();
// }

void print_parameters_siso(struct parameters *p, FILE *fptr) 
{
  if(fptr==NULL) fptr=stdout;
  fprintf(fptr, "SM %d\n", p->SM);
  fprintf(fptr, "model %d\n", p->model);
  fprintf(fptr, "generator %d\n", p->generator);
  fprintf(fptr, "Q %.4e\n", p->Q);
  fprintf(fptr, "m0 %.4e\n", p->m0);
  fprintf(fptr, "m12 %.4e\n", p->m12);
  fprintf(fptr, "tan_beta %.4e\n", p->tan_beta);
  fprintf(fptr, "sign_mu %.4e\n", p->sign_mu);
  fprintf(fptr, "A0 %.4e\n", p->A0);
  fprintf(fptr, "Lambda %.4e\n", p->Lambda);
  fprintf(fptr, "Mmess %.4e\n", p->Mmess);
  fprintf(fptr, "N5 %.4e\n", p->N5);
  fprintf(fptr, "cgrav %.4e\n", p->cgrav);
  fprintf(fptr, "m32 %.4e\n", p->m32);
  fprintf(fptr, "mass_Z %.4e\n", p->mass_Z);
  fprintf(fptr, "mass_W %.4e\n", p->mass_W);
  fprintf(fptr, "mass_b %.4e\n", p->mass_b);
  fprintf(fptr, "mass_top_pole %.4e\n", p->mass_top_pole);
  fprintf(fptr, "mass_tau_pole %.4e\n", p->mass_tau_pole);
  fprintf(fptr, "inv_alpha_em %.4e\n", p->inv_alpha_em);
  fprintf(fptr, "alphas_MZ %.4e\n", p->alphas_MZ);
  fprintf(fptr, "Gfermi %.4e\n", p->Gfermi);
  fprintf(fptr, "GAUGE_Q %.4e\n", p->GAUGE_Q);
  fprintf(fptr, "charg_Umix:\n");
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      fprintf(fptr, "%.4e ", p->charg_Umix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fprintf(fptr, "charg_Vmix:\n");
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      fprintf(fptr, "%.4e ", p->charg_Vmix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fprintf(fptr, "stop_mix:\n");
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      fprintf(fptr, "%.4e ", p->stop_mix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fprintf(fptr, "sbot_mix:\n");
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      fprintf(fptr, "%.4e ", p->sbot_mix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fprintf(fptr, "stau_mix:\n");
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 3; j++) 
    {
      fprintf(fptr, "%.4e ", p->stau_mix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fprintf(fptr, "neut_mix:\n");
  for (int i = 0; i < 6; i++) 
  {
    for (int j = 0; j < 6; j++) 
    {
      fprintf(fptr, "%.4e ", p->neut_mix[i][j]);
    }
    fprintf(fptr, "\n");
  }
  
  fprintf(fptr, "mass_neut[0] %.4e\n", p->mass_neut[0]);
  fprintf(fptr, "mass_neut[1] %.4e\n", p->mass_neut[1]);
  fprintf(fptr, "mass_neut[2] %.4e\n", p->mass_neut[2]);
  fprintf(fptr, "mass_neut[3] %.4e\n", p->mass_neut[3]);
  fprintf(fptr, "mass_neut[4] %.4e\n", p->mass_neut[4]);
  fprintf(fptr, "mass_neut[5] %.4e\n", p->mass_neut[5]);
  fprintf(fptr, "alpha %.4e\n", p->alpha);
  fprintf(fptr, "Min %.4e\n", p->Min);
  fprintf(fptr, "M1_Min %.4e\n", p->M1_Min);
  fprintf(fptr, "M2_Min %.4e\n", p->M2_Min);
  fprintf(fptr, "M3_Min %.4e\n", p->M3_Min);
  fprintf(fptr, "At_Min %.4e\n", p->At_Min);
  fprintf(fptr, "Ab_Min %.4e\n", p->Ab_Min);
  fprintf(fptr, "Atau_Min %.4e\n", p->Atau_Min);
  fprintf(fptr, "M2H1_Min %.4e\n", p->M2H1_Min);
  fprintf(fptr, "M2H2_Min %.4e\n", p->M2H2_Min);
  fprintf(fptr, "mu_Min %.4e\n", p->mu_Min);
  fprintf(fptr, "M2A_Min %.4e\n", p->M2A_Min);
  fprintf(fptr, "tb_Min %.4e\n", p->tb_Min);
  fprintf(fptr, "mA_Min %.4e\n", p->mA_Min);
  fprintf(fptr, "MeL_Min %.4e\n", p->MeL_Min);
  fprintf(fptr, "MmuL_Min %.4e\n", p->MmuL_Min);
  fprintf(fptr, "MtauL_Min %.4e\n", p->MtauL_Min);
  fprintf(fptr, "MeR_Min %.4e\n", p->MeR_Min);
  fprintf(fptr, "MmuR_Min %.4e\n", p->MmuR_Min);
  fprintf(fptr, "MtauR_Min %.4e\n", p->MtauR_Min);
  fprintf(fptr, "MqL1_Min %.4e\n", p->MqL1_Min);
  fprintf(fptr, "MqL2_Min %.4e\n", p->MqL2_Min);
  fprintf(fptr, "MqL3_Min %.4e\n", p->MqL3_Min);
  fprintf(fptr, "MuR_Min %.4e\n", p->MuR_Min);
  fprintf(fptr, "McR_Min %.4e\n", p->McR_Min);
  fprintf(fptr, "MtR_Min %.4e\n", p->MtR_Min);
  fprintf(fptr, "MdR_Min %.4e\n", p->MdR_Min);
  fprintf(fptr, "MsR_Min %.4e\n", p->MsR_Min);
  fprintf(fptr, "MbR_Min %.4e\n", p->MbR_Min);
  fprintf(fptr, "N51 %.4e\n", p->N51);
  fprintf(fptr, "N52 %.4e\n", p->N52);
  fprintf(fptr, "N53 %.4e\n", p->N53);
  fprintf(fptr, "M2H1_Q %.4e\n", p->M2H1_Q);
  fprintf(fptr, "M2H2_Q %.4e\n", p->M2H2_Q);
  fprintf(fptr, "mass_d %.4e\n", p->mass_d);
  
  fprintf(fptr, "mass_u %.4e\n",p->mass_u);
  fprintf(fptr, "mass_s %.4e\n",p->mass_s);
  fprintf(fptr, "mass_c %.4e\n",p->mass_c);
  fprintf(fptr, "mass_t %.4e\n",p->mass_t);
  fprintf(fptr, "mass_e %.4e\n",p->mass_e);
  fprintf(fptr, "mass_nue %.4e\n",p->mass_nue);
  fprintf(fptr, "mass_mu %.4e\n",p->mass_mu);
  fprintf(fptr, "mass_num %.4e\n",p->mass_num);
  fprintf(fptr, "mass_tau %.4e\n",p->mass_tau);
  fprintf(fptr, "mass_nut %.4e\n",p->mass_nut);
  fprintf(fptr, "mass_gluon %.4e\n",p->mass_gluon);
  fprintf(fptr, "mass_photon %.4e\n",p->mass_photon);
  fprintf(fptr, "mass_Z0 %.4e\n",p->mass_Z0);
  fprintf(fptr, "mass_h0 %.4e\n",p->mass_h0);
  fprintf(fptr, "mass_H0 %.4e\n",p->mass_H0);
  fprintf(fptr, "mass_A0 %.4e\n",p->mass_A0);
  fprintf(fptr, "mass_H %.4e\n",p->mass_H);
  fprintf(fptr, "mass_dnl %.4e\n",p->mass_dnl);
  fprintf(fptr, "mass_upl %.4e\n",p->mass_upl);
  fprintf(fptr, "mass_stl %.4e\n",p->mass_stl);
  fprintf(fptr, "mass_chl %.4e\n",p->mass_chl);
  fprintf(fptr, "mass_b1 %.4e\n",p->mass_b1);
  fprintf(fptr, "mass_t1 %.4e\n",p->mass_t1);
  fprintf(fptr, "mass_el %.4e\n",p->mass_el);
  fprintf(fptr, "mass_nuel %.4e\n",p->mass_nuel);
  fprintf(fptr, "mass_mul %.4e\n",p->mass_mul);
  fprintf(fptr, "mass_numl %.4e\n",p->mass_numl);
  fprintf(fptr, "mass_tau1 %.4e\n",p->mass_tau1);
  fprintf(fptr, "mass_nutl %.4e\n",p->mass_nutl);
  fprintf(fptr, "mass_gluino %.4e\n",p->mass_gluino);
  fprintf(fptr, "mass_cha1 %.4e\n",p->mass_cha1);
  fprintf(fptr, "mass_cha2 %.4e\n",p->mass_cha2);
  fprintf(fptr, "mass_dnr %.4e\n",p->mass_dnr);
  fprintf(fptr, "mass_upr %.4e\n",p->mass_upr);
  fprintf(fptr, "mass_str %.4e\n",p->mass_str);
  fprintf(fptr, "mass_chr %.4e\n",p->mass_chr);
  fprintf(fptr, "mass_b2 %.4e\n",p->mass_b2);
  fprintf(fptr, "mass_t2 %.4e\n",p->mass_t2);
  fprintf(fptr, "mass_er %.4e\n",p->mass_er);
  fprintf(fptr, "mass_mur %.4e\n",p->mass_mur);
  fprintf(fptr, "mass_tau2 %.4e\n",p->mass_tau2);
  fprintf(fptr, "mass_nuer %.4e\n",p->mass_nuer);
  fprintf(fptr, "mass_numr %.4e\n",p->mass_numr);
  fprintf(fptr, "mass_nutr %.4e\n",p->mass_nutr);
  fprintf(fptr, "mass_graviton %.4e\n",p->mass_graviton);
  fprintf(fptr, "mass_gravitino %.4e\n",p->mass_gravitino);
  
  fprintf(fptr, "gp %.4e\n", p->gp);
  fprintf(fptr, "g2 %.4e\n", p->g2);
  fprintf(fptr, "gp_Q %.4e\n", p->gp_Q);
  fprintf(fptr, "g2_Q %.4e\n", p->g2_Q);
  fprintf(fptr, "g3_Q %.4e\n", p->g3_Q);
  fprintf(fptr, "YU_Q %.4e\n", p->YU_Q);
  fprintf(fptr, "yut[4] %.4e %.4e %.4e %.4e\n", p->yut[0], p->yut[1], p->yut[2], p->yut[3]);
  fprintf(fptr, "YD_Q %.4e\n", p->YD_Q);
  fprintf(fptr, "yub[4] %.4e %.4e %.4e %.4e\n", p->yub[0], p->yub[1], p->yub[2], p->yub[3]);
  fprintf(fptr, "YE_Q %.4e\n", p->YE_Q);
  fprintf(fptr, "yutau[4] %.4e %.4e %.4e %.4e\n", p->yutau[0], p->yutau[1], p->yutau[2], p->yutau[3]);
  fprintf(fptr, "HMIX_Q %.4e\n", p->HMIX_Q);
  fprintf(fptr, "mu_Q %.4e\n", p->mu_Q);
  fprintf(fptr, "tanb_GUT %.4e\n", p->tanb_GUT);
  fprintf(fptr, "Higgs_VEV %.4e\n", p->Higgs_VEV);
  fprintf(fptr, "mA2_Q %.4e\n", p->mA2_Q);
  fprintf(fptr, "MSOFT_Q %.4e\n", p->MSOFT_Q);
  fprintf(fptr, "M1_Q %.4e\n", p->M1_Q);
  fprintf(fptr, "M2_Q %.4e\n", p->M2_Q);
  fprintf(fptr, "M3_Q %.4e\n", p->M3_Q);
  fprintf(fptr, "MeL_Q %.4e\n", p->MeL_Q);
  fprintf(fptr, "MmuL_Q %.4e\n", p->MmuL_Q);
  fprintf(fptr, "MtauL_Q %.4e\n", p->MtauL_Q);
  fprintf(fptr, "MeR_Q %.4e\n", p->MeR_Q);
  fprintf(fptr, "MmuR_Q %.4e\n", p->MmuR_Q);
  fprintf(fptr, "MtauR_Q %.4e\n", p->MtauR_Q);
  fprintf(fptr, "MqL1_Q %.4e\n", p->MqL1_Q);
  fprintf(fptr, "MqL2_Q %.4e\n", p->MqL2_Q);
  fprintf(fptr, "MqL3_Q %.4e\n", p->MqL3_Q);
  fprintf(fptr, "MuR_Q %.4e\n", p->MuR_Q);
  fprintf(fptr, "McR_Q %.4e\n", p->McR_Q);
  fprintf(fptr, "MtR_Q %.4e\n", p->MtR_Q);
  fprintf(fptr, "MdR_Q %.4e\n", p->MdR_Q);
  fprintf(fptr, "MsR_Q %.4e\n", p->MsR_Q);
  fprintf(fptr, "MbR_Q %.4e\n", p->MbR_Q);
  
  fprintf(fptr, "AU_Q %.4e\n", p->AU_Q);
  fprintf(fptr, "A_u %.4e\n", p->A_u);
  fprintf(fptr, "A_c %.4e\n", p->A_c);
  fprintf(fptr, "A_t %.4e\n", p->A_t);
  fprintf(fptr, "AD_Q %.4e\n", p->AD_Q);
  fprintf(fptr, "A_d %.4e\n", p->A_d);
  fprintf(fptr, "A_s %.4e\n", p->A_s);
  fprintf(fptr, "A_b %.4e\n", p->A_b);
  fprintf(fptr, "AE_Q %.4e\n", p->AE_Q);
  fprintf(fptr, "A_e %.4e\n", p->A_e);
  fprintf(fptr, "A_mu %.4e\n", p->A_mu);
  fprintf(fptr, "A_tau %.4e\n", p->A_tau);
  fprintf(fptr, "NMSSM %d\n", p->NMSSM);
  fprintf(fptr, "RV %d\n", p->RV);
  fprintf(fptr, "CPV %d\n", p->CPV);
  fprintf(fptr, "FV %d\n", p->FV);
  fprintf(fptr, "CKM_lambda %.4e\n", p->CKM_lambda);
  fprintf(fptr, "CKM_A %.4e\n", p->CKM_A);
  fprintf(fptr, "CKM_rhobar %.4e\n", p->CKM_rhobar);
  fprintf(fptr, "CKM_etabar %.4e\n", p->CKM_etabar);
  fprintf(fptr, "PMNS_theta12 %.4e\n", p->PMNS_theta12);
  fprintf(fptr, "PMNS_theta23 %.4e\n", p->PMNS_theta23);
  fprintf(fptr, "PMNS_theta13 %.4e\n", p->PMNS_theta13);
  fprintf(fptr, "PMNS_delta13 %.4e\n", p->PMNS_delta13);
  fprintf(fptr, "PMNS_alpha1 %.4e\n", p->PMNS_alpha1);
  fprintf(fptr, "PMNS_alpha2 %.4e\n", p->PMNS_alpha2);
  
  
  fprintf(fptr, "NMSSM parameters:\n");
  
  fprintf(fptr, "lambdaNMSSM_Min %.4e\n", p->lambdaNMSSM_Min);
  fprintf(fptr, "kappaNMSSM_Min %.4e\n", p->kappaNMSSM_Min);
  fprintf(fptr, "AlambdaNMSSM_Min %.4e\n", p->AlambdaNMSSM_Min);
  fprintf(fptr, "AkappaNMSSM_Min %.4e\n", p->AkappaNMSSM_Min);
  fprintf(fptr, "lambdaSNMSSM_Min %.4e\n", p->lambdaSNMSSM_Min);
  fprintf(fptr, "xiFNMSSM_Min %.4e\n", p->xiFNMSSM_Min);
  fprintf(fptr, "xiSNMSSM_Min %.4e\n", p->xiSNMSSM_Min);
  fprintf(fptr, "mupNMSSM_Min %.4e\n", p->mupNMSSM_Min);
  fprintf(fptr, "mSp2NMSSM_Min %.4e\n", p->mSp2NMSSM_Min);
  fprintf(fptr, "mS2NMSSM_Min %.4e\n", p->mS2NMSSM_Min);
  fprintf(fptr, "mass_H03 %.4e\n", p->mass_H03);
  fprintf(fptr, "mass_A02 %.4e\n", p->mass_A02);
  
  fprintf(fptr, "NMSSMRUN_Q %.4e\n", p->NMSSMRUN_Q);
  fprintf(fptr, "lambdaNMSSM %.4e\n", p->lambdaNMSSM);
  fprintf(fptr, "kappaNMSSM %.4e\n", p->kappaNMSSM);
  fprintf(fptr, "AlambdaNMSSM %.4e\n", p->AlambdaNMSSM);
  fprintf(fptr, "AkappaNMSSM %.4e\n", p->AkappaNMSSM);
  fprintf(fptr, "lambdaSNMSSM %.4e\n", p->lambdaSNMSSM);
  fprintf(fptr, "xiFNMSSM %.4e\n", p->xiFNMSSM);
  fprintf(fptr, "xiSNMSSM %.4e\n", p->xiSNMSSM);
  fprintf(fptr, "mupNMSSM %.4e\n", p->mupNMSSM);
  fprintf(fptr, "mSp2NMSSM %.4e\n", p->mSp2NMSSM);
  fprintf(fptr, "mS2NMSSM %.4e\n", p->mS2NMSSM);
  fprintf(fptr, "\n");

  fprintf(fptr, "CKM matrix:\n");
  for (int i=0; i<=3; i++){
    for (int j=0; j<=3; j++){
      fprintf(fptr, "CKM(%d,%d) %.4e\n", i, j, p->CKM[i][j]);
    }
  }
//   fprintf(fptr, "CKM(1,4) %.4e\n", p->CKM[1][4]);
//   fprintf(fptr, "CKM(2,4) %.4e\n", p->CKM[2][4]);
//   fprintf(fptr, "CKM(3,4) %.4e\n", p->CKM[3][4]);
//   fprintf(fptr, "CKM(4,4) %.4e\n", p->CKM[4][4]);
  fprintf(fptr, "\n");

  fprintf(fptr, "IMCKM matrix:\n");
  for (int i=0; i<=3; i++){
    for (int j=0; j<=3; j++){
      fprintf(fptr, "IMCKM(%d,%d) %.4e\n", i, j, p->IMCKM[i][j]);
    }
  }
//   fprintf(fptr, "IMCKM(1,4) %.4e\n", p->IMCKM[1][4]);
//   fprintf(fptr, "IMCKM(2,4) %.4e\n", p->IMCKM[2][4]);
//   fprintf(fptr, "IMCKM(3,4) %.4e\n", p->IMCKM[3][4]);
//   fprintf(fptr, "IMCKM(4,4) %.4e\n", p->IMCKM[4][4]);
  fprintf(fptr, "\n");

  fprintf(fptr, "Higgs mixing matrices:\n");
  for (int i=0; i<=3; i++){
    for (int j=0; j<=3; j++){
      fprintf(fptr, "H0_mix(%d,%d) %.4e\n", i, j, p->H0_mix[i][j]);
    }
  }
//   fprintf(fptr, "H0_mix(1,4) %.4e\n", p->H0_mix[1][4]);
//   fprintf(fptr, "H0_mix(2,4) %.4e\n", p->H0_mix[2][4]);
//   fprintf(fptr, "H0_mix(3,4) %.4e\n", p->H0_mix[3][4]);
//   fprintf(fptr, "H0_mix(4,4) %.4e\n", p->H0_mix[4][4]);
  
  fprintf(fptr, "A0_mix %.4e\n", p->A0_mix[0][0]);
  fprintf(fptr, "sU_mix %.4e\n", p->sU_mix[0][0]);
  fprintf(fptr, "sD_mix %.4e\n", p->sD_mix[0][0]);
  fprintf(fptr, "sE_mix %.4e\n", p->sE_mix[0][0]);
  fprintf(fptr, "sNU_mix %.4e\n", p->sNU_mix[0][0]);
  fprintf(fptr, "sCKM_msq2 %.4e\n", p->sCKM_msq2[0][0]);
  fprintf(fptr, "sCKM_msl2 %.4e\n", p->sCKM_msl2[0][0]);
  fprintf(fptr, "sCKM_msd2 %.4e\n", p->sCKM_msd2[0][0]);
  fprintf(fptr, "sCKM_msu2 %.4e\n", p->sCKM_msu2[0][0]);
  fprintf(fptr, "sCKM_mse2 %.4e\n", p->sCKM_mse2[0][0]);
  fprintf(fptr, "PMNS_U %.4e\n", p->PMNS_U[0][0]);
  fprintf(fptr, "TU %.4e\n", p->TU[0][0]);
  fprintf(fptr, "TD %.4e\n", p->TD[0][0]);
  fprintf(fptr, "TE %.4e\n", p->TE[0][0]);
  fprintf(fptr, "mass_c_pole %.4e\n", p->mass_c_pole);
  fprintf(fptr, "mass_b_1S %.4e\n", p->mass_b_1S);
  fprintf(fptr, "mass_b_pole %.4e\n", p->mass_b_pole);
  fprintf(fptr, "mtmt %.4e\n", p->mtmt);
  fprintf(fptr, "scheme_c_mass %d\n", p->scheme_c_mass);
  fprintf(fptr, "Lambda3 %.4e\n", p->Lambda3);
  fprintf(fptr, "Lambda4 %.4e\n", p->Lambda4);
  fprintf(fptr, "Lambda5 %.4e\n", p->Lambda5);
  fprintf(fptr, "Lambda6 %.4e\n", p->Lambda6);
  fprintf(fptr, "alphasMZ_Lambda3 %.4e\n", p->alphasMZ_Lambda3);
  fprintf(fptr, "alphasMZ_Lambda4 %.4e\n", p->alphasMZ_Lambda4);
  fprintf(fptr, "alphasMZ_Lambda5 %.4e\n", p->alphasMZ_Lambda5);
  fprintf(fptr, "alphasMZ_Lambda6 %.4e\n", p->alphasMZ_Lambda6);
  fprintf(fptr, "Vud %.4e + %.4ei\n", creal(p->Vud), cimag(p->Vud));
  fprintf(fptr, "Vus %.4e + %.4ei\n", creal(p->Vus), cimag(p->Vus));
  fprintf(fptr, "Vub %.4e + %.4ei\n", creal(p->Vub), cimag(p->Vub));
  fprintf(fptr, "Vcd %.4e + %.4ei\n", creal(p->Vcd), cimag(p->Vcd));
  fprintf(fptr, "Vcs %.4e + %.4ei\n", creal(p->Vcs), cimag(p->Vcs));
  fprintf(fptr, "Vcb %.4e + %.4ei\n", creal(p->Vcb), cimag(p->Vcb));
  
  fprintf(fptr, "Vtd %.4e\n", creal(p->Vtd));
  fprintf(fptr, "Vts %.4e\n", creal(p->Vts));
  fprintf(fptr, "Vtb %.4e\n", creal(p->Vtb));

  fprintf(fptr, "widthcalc %d\n", p->widthcalc);
  fprintf(fptr, "width_h0 %.4e\n", p->width_h0);
  fprintf(fptr, "width_H0 %.4e\n", p->width_H0);
  fprintf(fptr, "width_A0 %.4e\n", p->width_A0);
  fprintf(fptr, "width_H %.4e\n", p->width_H);
  fprintf(fptr, "width_Z %.4e\n", p->width_Z);
  fprintf(fptr, "width_W %.4e\n", p->width_W);
  fprintf(fptr, "width_top %.4e\n", p->width_top);
  fprintf(fptr, "width_H03 %.4e\n", p->width_H03);
  fprintf(fptr, "width_A02 %.4e\n", p->width_A02);
  fprintf(fptr, "width_gluino %.4e\n", p->width_gluino);
  fprintf(fptr, "width_t1 %.4e\n", p->width_t1);
  fprintf(fptr, "width_t2 %.4e\n", p->width_t2);
  fprintf(fptr, "width_b1 %.4e\n", p->width_b1);
  fprintf(fptr, "width_b2 %.4e\n", p->width_b2);
  fprintf(fptr, "width_ul %.4e\n", p->width_ul);
  fprintf(fptr, "width_ur %.4e\n", p->width_ur);
  fprintf(fptr, "width_dl %.4e\n", p->width_dl);
  fprintf(fptr, "width_dr %.4e\n", p->width_dr);
  fprintf(fptr, "width_cl %.4e\n", p->width_cl);
  fprintf(fptr, "width_cr %.4e\n", p->width_cr);
  fprintf(fptr, "width_sl %.4e\n", p->width_sl);
  fprintf(fptr, "width_sr %.4e\n", p->width_sr);
  fprintf(fptr, "width_el %.4e\n", p->width_el);
  fprintf(fptr, "width_er %.4e\n", p->width_er);
  fprintf(fptr, "width_ml %.4e\n", p->width_ml);
  fprintf(fptr, "width_mr %.4e\n", p->width_mr);
  fprintf(fptr, "width_tau1 %.4e\n", p->width_tau1);
  fprintf(fptr, "width_tau2 %.4e\n", p->width_tau2);
  fprintf(fptr, "width_gravitino %.4e\n", p->width_gravitino);
  fprintf(fptr, "width_nuel %.4e\n", p->width_nuel);
  fprintf(fptr, "width_numl %.4e\n", p->width_numl);
  fprintf(fptr, "width_nutaul %.4e\n", p->width_nutaul);
  fprintf(fptr, "width_c1 %.4e\n", p->width_c1);
  fprintf(fptr, "width_c2 %.4e\n", p->width_c2);
  fprintf(fptr, "width_o1 %.4e\n", p->width_o1);
  fprintf(fptr, "width_o2 %.4e\n", p->width_o2);
  fprintf(fptr, "width_o3 %.4e\n", p->width_o3);
  fprintf(fptr, "width_o4 %.4e\n", p->width_o4);
  fprintf(fptr, "width_o5 %.4e\n", p->width_o5);
}

