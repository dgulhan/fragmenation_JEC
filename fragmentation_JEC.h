#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;
struct fragmentation_JEC
{
 private:
  static const int ncent=4;
  static const int nstepmax=10;
  static const double lower_pt_cut=15;
  static const double higher_pt_cut=700;
  static const int ntrkmax=14;
  int radius;
  int nstep;
  int cent_max[ncent];
  int cent_min[ncent];
  double PF_pt_cut;
  double PF_eta_cut;
  int jetType;
  bool do_PbPb;
  bool do_pp_tracking;
  bool do_residual_correction;
  TString algo_corr; 
  TF1 *correction_matrix[ntrkmax][ncent];
  TF1 *residual_correction_function[ncent][nstepmax];
  TH1D *hist_eff_pt[ncent];
  TH1D *hist_eff_eta[ncent];
  TFile *correction_file;
  TFile *residual_correction_file[nstepmax];
  TFile *eff_file; 

  public:
  void reset()
  { 
   for(int iNpf = 0; iNpf < ntrkmax ;iNpf++){
    for(int icent=0;icent<ncent;icent++){
     correction_matrix[icent]=NULL;
    }
   }
   correction_file=NULL;
   nstep=1;
   PF_eta_cut=2.4;
   cent_min[0]=0;
   cent_max[0]=cent_min[1]=20;
   cent_max[1]=cent_min[2]=60;
   cent_max[2]=cent_min[3]=100;
   cent_max[3]=200;
  }
  
  
  
  fragmentation_JEC(int radius=3, bool do_PbPb=1, bool do_pp_tracking=0, bool do_residual_correction=1,
  int nstep=1, double PF_pt_cut=2, int type=0)
  {
   reset();
   jetType=type;
   cout <<"JetType = "<<jetType<<endl;
   if(do_PbPb==1){
    do_pp_tracking=0;
   }
   this->do_PbPb=do_PbPb;
   this->radius=radius;
   this->PF_pt_cut=PF_pt_cut;
   this->do_pp_tracking=do_pp_tracking;
   this->do_residual_correction=do_residual_correction;
   // if(PF_pt_cut==3) ntrkmax=26;
   // else if(PF_pt_cut==2) ntrkmax=31;
   // else if(PF_pt_cut==1) ntrkmax=41;
   this->nstep=nstep;
   if (jetType==0) {
    algo_corr=Form("akVs%dCalo",radius);
   } else if (jetType==1) {
    algo_corr=Form("akVs%dPF",radius);
   } else if (jetType==2) {
    algo_corr=Form("akPu%dCalo",radius);
   }    
  }

  
  void set_correction()
  {
   cout<<"setting correction"<<endl;
   if(do_PbPb){
    if(jetType==0) {
   if(do_PbPb){
    correction_file = new TFile(Form("corrections_2015_06_09/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    for(int icent=0;icent<ncent;icent++){
     for(int iNpf=0; iNpf<ntrkmax;iNpf++){	  
	  correction_matrix[iNpf][icent]=(TF1*)correction_file->Get(Form("fit_total_NPF_%d_cent_%d_%d",iNpf,(int)(cent_min[icent]*0.5),(int)(cent_max[icent]*0.5)));
	 }
    } 
     if(do_residual_correction){
      for(int istep=0;istep<nstep;istep++){
       residual_correction_file[istep] = new TFile(Form("corrections_2014_12_12_PbPb/residualcorr%d_%s.root",istep,algo_corr.Data()));
       for(int icent=0;icent<ncent;icent++){
        residual_correction_function[icent][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",icent));
       }
      }
     } 
    } 
	// else if (jetType==1) {
     // correction_file = new TFile(Form("corrections_2015_02_09_PbPb_PF/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     // for(int icent=0;icent<ncent;icent++){
 	 // correction_matrix[icent]=(TH2D*)correction_file->Get(Form("pNtrk_pt%d",icent));
     // } 
     
     // if(do_residual_correction){
      // for(int istep=0;istep<nstep;istep++){
       // residual_correction_file[istep] = new TFile(Form("corrections_2015_02_09_PbPb_PF/residualcorr%d_%s.root",istep,algo_corr.Data()));
       // for(int icent=0;icent<ncent;icent++){
        // residual_correction_function[icent][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",icent));
       // }
      // }
     // } 
    // } else if (jetType==2) {
     // correction_file = new TFile(Form("corrections_2015_02_09_PbPb_Pu/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     // for(int icent=0;icent<ncent;icent++){
 	 // correction_matrix[icent]=(TH2D*)correction_file->Get(Form("pNtrk_pt%d",icent));
     // } 
     
     if(do_residual_correction){
      for(int istep=0;istep<nstep;istep++){
       residual_correction_file[istep] = new TFile(Form("corrections_2015_02_09_PbPb_Pu/residualcorr%d_%s.root",istep,algo_corr.Data()));
       for(int icent=0;icent<ncent;icent++){
        residual_correction_function[icent][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",icent));
       }
      }
     } 
        
    }	
   }
   // else{
  	// algo_corr=Form("ak%dCalo",radius);
    // if(do_pp_tracking){
     // correction_file = new TFile(Form("corrections_2014_12_20_pp/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     // correction_matrix[0][0]=(TF1*)correction_file->Get("pNtrk_pt");    
     
     // if(do_residual_correction){
      // for(int istep=0;istep<nstep;istep++){
       // residual_correction_file[istep] = new TFile(Form("corrections_2014_12_20_pp/residualcorr%d_%s.root",istep,algo_corr.Data()));
       // residual_correction_function[0][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",0));
	  // }
     // }
    // }else{// correction for all R values are not available for HI tracking for the moment     
     // correction_file = new TFile(Form("corrections_2015_02_02_pp_HI_tracking/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     // correction_matrix[0][0]=(TF1*)correction_file->Get("pNtrk_pt");
     // if(do_residual_correction){
      // for(int istep=0;istep<nstep;istep++){
       // residual_correction_file[istep] = new TFile(Form("corrections_2015_02_02_pp_HI_tracking/residualcorr%d_%s.root",istep,algo_corr.Data()));
       // residual_correction_function[0][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",0));
	    // }
     // }
    // }
   // }
  }
  
  
  bool passes_PF_selection(double PF_pt, double PF_eta, double PF_phi, int PF_id, double jet_eta, double jet_phi)
  { 
   if(PF_pt<PF_pt_cut || fabs(PF_eta)>PF_eta_cut || PF_id!=1) return false;
   
   double r=sqrt(pow(jet_eta-PF_eta,2)+pow(acos(cos(jet_phi-PF_phi)),2));
   
   if(r<((double)(radius)*0.1)) return true;
   else return false;
  }
  
  double get_corrected_pt(double jetpt, int ntrk, int cent=0)
  {
   //correction for fragmentation dependent JEC as a function of number of charged particle flow candidates and reconstructed jet pt
   double correction=1;
   
   int cent_bin=0;
   if(do_PbPb){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   if(jetpt<lower_pt_cut) return jetpt; 
   
   double jetpt_for_correction=jetpt;
   if(jetpt<25) jetpt_for_correction=25;
   if(jetpt>=700) jetpt_for_correction=699;
   
   if(ntrk>14) ntrk = 14;
   correction=correction_matrix[ntrk][cent_bin]->Eval(jetpt_for_correction);
   
   return (1/correction)*jetpt;
  }
    
  double get_residual_corrected_pt(double corrected_jetpt, int cent=0)
  {
   // residual correction to correct for the effects of jet resolution in fragmentation jec with a simple centrality binned fit function
   double residual_correction=1;
   
   int cent_bin=0;
   if(do_PbPb){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   if(corrected_jetpt<lower_pt_cut) return corrected_jetpt;
   
   double jetpt_for_correction=corrected_jetpt;
   if(corrected_jetpt<25) jetpt_for_correction=25;
   if(corrected_jetpt>=700) jetpt_for_correction=699;

   for(int istep=0;istep<nstep;istep++){
    residual_correction=residual_correction*residual_correction_function[cent_bin][istep]->Eval((1/(residual_correction))*jetpt_for_correction);
   }
   return (1/(residual_correction))*corrected_jetpt;
  }
  
  
};
