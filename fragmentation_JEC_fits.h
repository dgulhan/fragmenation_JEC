#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;
class fragmentation_JEC
{
 private:
  static const int ncent=4;
  static const int nstepmax=10;
  static const double lower_pt_cut=5;
  static const double second_jtpt_cut = 25;
  static const double higher_pt_cut=700;
  int radius;
  int ntrkmax;
  int nstep;
  int cent_max[ncent];
  int cent_min[ncent];
  double PF_pt_cut;
  double PF_eta_cut;
  bool do_PbPb;
  bool do_pp_tracking;
  bool do_residual_correction;
  TString algo_corr; 
  TF1 *correction_matrix[15][ncent];
  TFile *correction_file;
  TFile *residual_correction_file[nstepmax];
  TFile *eff_file_pt; 
  TFile *eff_file_eta; 
  TF1 * residual_correction_function[ncent][nstepmax];
  public:
  void reset()
  { 
   for(int iNpf = 0; iNpf < 15 ;iNpf++){
    for(int icent=0;icent<ncent;icent++){
     correction_matrix[iNpf][icent]=NULL;
    }
   }
   correction_file=NULL;
   ntrkmax=15;
   nstep=1;
   PF_eta_cut=2.4;
   cent_min[0]=0;
   cent_max[0]=cent_min[1]=20;
   cent_max[1]=cent_min[2]=60;
   cent_max[2]=cent_min[3]=100;
   cent_max[3]=200;
  }
  
  fragmentation_JEC(int radius=3, bool do_PbPb=1, bool do_pp_tracking=0, bool do_residual_correction=1, int nstep=1, double PF_pt_cut=2)
  {
   reset();
   if(do_PbPb==1){
    do_pp_tracking=0;
   }
   this->do_PbPb=do_PbPb;
   this->radius=radius;
   this->PF_pt_cut=PF_pt_cut;
   this->do_pp_tracking=do_pp_tracking;
   this->do_residual_correction=do_residual_correction;
   if(PF_pt_cut==3) ntrkmax=14;
   else if(PF_pt_cut==2) ntrkmax=14;
   else if(PF_pt_cut==1) ntrkmax=14;
   this->nstep=nstep;
   algo_corr=Form("akVs%dCalo",radius);
   // algo_corr=Form("akPu%dCalo",radius);
  }
   
  void set_correction()
  {
   //std::cout<<"setting correction"<<endl;
   if(do_PbPb){
    correction_file = new TFile(Form("corrections_2015_06_08_PbPb/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    // correction_file = new TFile(Form("corrections_2015_02_09_PbPb/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    // correction_file = new TFile(Form("corrections_2015_05_01_PbPb/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    for(int icent=0;icent<ncent;icent++){
     for(int iNpf=0; iNpf<15;iNpf++){	  
	  correction_matrix[iNpf][icent]=(TF1*)correction_file->Get(Form("fit_total_NPF_%d_cent_%d_%d",iNpf,(int)(cent_min[icent]*0.5),(int)(cent_max[icent]*0.5)));
	 }
    } 
    
    if(do_residual_correction){
     for(int istep=0;istep<nstep;istep++){
      // residual_correction_file[istep] = new TFile(Form("corrections_2014_12_12_PbPb/residualcorr%d_%s.root",istep,algo_corr.Data()));
      // residual_correction_file[istep] = new TFile(Form("corrections_2015_06_08_PbPb/residualcorr%d_%s.root",istep,algo_corr.Data()));
      residual_correction_file[istep] = new TFile(Form("corrections_2015_06_14_PbPb/residualcorr%d_%s.root",istep,algo_corr.Data()));
      for(int icent=0;icent<ncent;icent++){
       residual_correction_function[icent][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",icent));
      }
     }
    }
	
   }else{
  	algo_corr=Form("ak%dCalo",radius);
    if(do_pp_tracking){
     correction_file = new TFile(Form("corrections_2015_06_11_pp/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     for(int iNpf=0; iNpf<15;iNpf++){	  
      correction_matrix[iNpf][0]=(TF1*)correction_file->Get(Form("fit_total_NPF_%d",iNpf));    
     }
     if(do_residual_correction){
      for(int istep=0;istep<nstep;istep++){
       // residual_correction_file[istep] = new TFile(Form("corrections_2015_06_11_pp/residualcorr%d_%s.root",istep,algo_corr.Data()));
       residual_correction_file[istep] = new TFile(Form("corrections_2015_07_20_pp/residualcorr%d_%s.root",istep,algo_corr.Data()));
       residual_correction_function[0][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",0));
	  }
     }
    }else{//! correction for all R values are not available for HI tracking for the moment
     correction_file = new TFile(Form("corrections_2015_02_02_pp_HI_tracking/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     correction_matrix[0][0]=(TF1*)correction_file->Get("pNtrk_pt");
     if(do_residual_correction){
      for(int istep=0;istep<nstep;istep++){
       residual_correction_file[istep] = new TFile(Form("corrections_2015_02_02_pp_HI_tracking/residualcorr%d_%s.root",istep,algo_corr.Data()));
       residual_correction_function[0][istep] = (TF1*)residual_correction_file[istep]->Get(Form("fit%d",0));
	    }
     }
    }
   }
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
   //correction for fragmentation dependent JEC as a function of efficiency corrected number of charged particle flow candidates and reconstructed jet pt
   double correction=1;
   
   int cent_bin=0;
   if(do_PbPb){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   if(jetpt<lower_pt_cut) return jetpt; 
   
   double jetpt_for_correction=jetpt;
   if(jetpt<second_jtpt_cut) jetpt_for_correction=second_jtpt_cut;
   if(jetpt>=700) jetpt_for_correction=699;
   
   if(ntrk>14) ntrk = 14;
   //std::cout<<jetpt<<" "<<ntrk<<" " <<cent<<endl;
   correction=correction_matrix[ntrk][cent_bin]->Eval(jetpt_for_correction);
   //std::cout<<jetpt<<" "<<ntrk<<endl;
   
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
      
   for(int istep=0;istep<nstep;istep++){
    double jetpt_for_correction=corrected_jetpt;
    if(jetpt_for_correction<second_jtpt_cut) jetpt_for_correction=second_jtpt_cut;
    if(jetpt_for_correction>=700) jetpt_for_correction=699;
    residual_correction=residual_correction_function[cent_bin][istep]->Eval(jetpt_for_correction);
    jetpt_for_correction=(1/(residual_correction))*jetpt_for_correction;
    if(jetpt_for_correction<second_jtpt_cut) jetpt_for_correction=second_jtpt_cut;
    if(jetpt_for_correction>=700) jetpt_for_correction=699;
    residual_correction=residual_correction_function[cent_bin][istep]->Eval(jetpt_for_correction);
    corrected_jetpt=(1/(residual_correction))*corrected_jetpt;
   }
   
   return corrected_jetpt;
  }
  
  double get_residual_corrected_pt_single_step(double corrected_jetpt, int cent=0)
  {
   // residual correction to correct for the effects of jet resolution in fragmentation jec with a simple centrality binned fit function
   double residual_correction=1;
   
   int cent_bin=0;
   if(do_PbPb){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   double jetpt_for_correction=corrected_jetpt;
   if(corrected_jetpt<lower_pt_cut) return corrected_jetpt;

   if(corrected_jetpt<second_jtpt_cut) jetpt_for_correction=second_jtpt_cut;
   if(corrected_jetpt>700) jetpt_for_correction=700;

   residual_correction=residual_correction*residual_correction_function[cent_bin][0]->Eval(jetpt_for_correction);
	
   return (1/(residual_correction))*corrected_jetpt;
  }
  
};