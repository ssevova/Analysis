//================================================================================================
//
// Performs preselection to produce bacon bits for training and testing top MVA ID
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                          // access to gROOT, entry point to ROOT system
#include <TSystem.h>                        // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                          // file handle class
#include <TTree.h>                          // class to access ntuples
#include <TH1D.h>                           // 1D histogram class
#include <TLorentzVector.h>                 // 4-vector class
#include <TVector2.h>                       // 2-vector class
#include <vector>                           // STL vector class
#include <iostream>                         // standard I/O
#include <iomanip>                          // functions to format standard I/O
#include <fstream>                          // functions for file I/O
#include <string>                           // C++ string class
#include <cmath>                            // C++ math library
#include <cassert>
#include <utility>                          //std::pair

#include "/tthome/ssevova/CMSSW_5_3_14_patch2/src/DMSAna/Utils/interface/QGSyst.h"  // class to handle smearing correction for MC q/g discriminator

#include <TMVA/Reader.h>

#include "CPlot.hh"
#include "KStyle.hh"
#include "CSample.hh"

#endif

using namespace std;


//=== FUNCTION DECLARATIONS ======================================================================================
template <typename Iterator>
inline bool next_combination(const Iterator first, Iterator k, const Iterator last)
{
  if((first == last) || (first == k) || (last == k)) return false;

  Iterator itr1 = first;
  Iterator itr2 = last;
  ++itr1;

  if(last == itr1) return false;
  
  itr1 = last;
  --itr1;
  itr1 = k;
  --itr2;
  
  while (first != itr1)
    {
      if(*--itr1 < *itr2)
	{
	  Iterator j = k;
	  while (!(*itr1 < *j)) ++j;
	  std::iter_swap(itr1, j);
	  ++itr1;
	  ++j;
	  itr2 = k;
	  std::rotate(itr1,j,last);
	  while(last != j)
	    {
	      ++j;
	      ++itr2;
	    }
	  std::rotate(k, itr2, last);
	  return true;
	}
    }
  std::rotate(first,k,last);
  return false;
}
//------------------------------------------------------------------------------------------------

double deltaPhi(const double phi1, const double phi2) {
  double result = phi1 - phi2;
  if     (result >  TMath::Pi()) { result = result - 2*TMath::Pi(); }
  else if(result < -TMath::Pi()) { result = result + 2*TMath::Pi(); }
  return result;
}

//--------------------------------------------------------------------------------------------------
float maximum( float a, float b, float c )
  {
    float max = ( a < b ) ? b : a;
    return ( ( max < c ) ? c : max );
  }
//--------------------------------------------------------------------------------------------------

//=== MAIN MACRO =================================================================================================

void mvaSelect()
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  
  //
  // Preselection cuts
  //
  const unsigned int NJETS_CUT   = 4;
  const unsigned int NBJETS_CUT  = 2;
  const double       WMASSLOW    = 60;
  const double       WMASSHIGH   = 110;
  const double       TOPMASSLOW  = 100;
  const double       TOPMASSHIGH = 300;

  //
  // input/output file
  //  
  const string infilename("Summer12_TTJets_SemiLeptMGDecays_ttbarbits.root");
  const string outfilename("trainingbits.root");

  const string weightfile_lowpt("/tthome/ssevova/CMSSW_5_3_14_patch2/src/DMSAna/TTbar/topmva/wmva_weights/wtraining_lowpt_cen.root_BDTG.weights.xml");
  const string weightfile_highpt("/tthome/ssevova/CMSSW_5_3_14_patch2/src/DMSAna/TTbar/topmva/wmva_weights/wtraining_highpt_cen.root_BDTG.weights.xml");

  const string puWeightFilename("/tthome/ssevova/CMSSW_5_3_14_patch2/src/DMSAna/Utils/data/PUWeights_2012.root");
  const string qgidCorrFilename("/tthome/ssevova/CMSSW_5_3_14_patch2/src/DMSAna/Utils/data/SystDatabase.txt");  
  const string muonEffFilename("");
  const string eleEffFilename("");
  
  // plot output directory
  const string outputDir("TopPlots");

  gSystem->mkdir(outputDir.c_str(), true);
  CPlot::sOutDir = outputDir;
  const string format("png");
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================
  
  //
  // Declare variables to read in ntuple
  //

  unsigned int    metfilter;                                                 // MET filter bits
  unsigned int    npv, npu;                                                  // number of PV / PU
  unsigned int    njets, nbjets;                                             // jet multiplicity
  float           rho;                                                       // event energy density
  float           scale1fb;                                                  // cross section scale factor per 1/fb
  float           ttWeight;                                                  // ttbar event weight (needed for MadGraph ttbar)
  int             lepId;                                                     // lepton PDG ID
  TLorentzVector *lepton=0;                                                  //lep from leptonic decay
  float           pfmet;
  int             jet1flavGen,  jet2flavGen,  jet3flavGen,  jet4flavGen;
  float           jet1qgid,     jet2qgid,     jet3qgid,     jet4qgid;      // jet q/g discriminant
  float           jet1csv,      jet2csv,      jet3csv,      jet4csv;       // jet CSV b-tagger
  float           jet1flavAlgo, jet2flavAlgo, jet3flavAlgo, jet4flavAlgo;  // jet flavor, algorithmic definition (MC only)
  float           jet1q,        jet2q,        jet3q,        jet4q;         // jet charge (kappa=1)
  TLorentzVector *jet1=0,      *jet2=0,      *jet3=0,      *jet4=0;        // jet 4-vector
  TVector2       *jet1pull=0,  *jet2pull=0,  *jet3pull=0,  *jet4pull=0;    // jet pull vector
 

  //
  // Set up output file
  //
  unsigned int isSig, b_mis, w_mis, wb_mis;
  float wmva, pttop, mtop, detaj1b, detaj2b, dphij1b, dphij2b, drj1b, drj2b;
  float bjcsv, j1csv, j2csv;
  float weight;
  TLorentzVector *bjet=0;
  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  TTree *outTree = new TTree("Events","Events");  

  outTree->Branch("isSig",  &isSig,  "isSig/i");
  outTree->Branch("b_mis",  &b_mis,  "b_mis/i");
  outTree->Branch("w_mis",  &w_mis,  "w_mis/i");
  outTree->Branch("wb_mis", &wb_mis, "wb_mis/i");
  outTree->Branch("wmva",   &wmva,   "wmva/F");
  outTree->Branch("pttop",  &pttop,  "pttop/F");
  outTree->Branch("mtop",   &mtop,   "mtop/F");
  outTree->Branch("detaj1b",&detaj1b,"detaj1b/F");
  outTree->Branch("detaj2b",&detaj2b,"detaj2b/F");
  outTree->Branch("dphij1b",&dphij1b,"dphij1b/F");
  outTree->Branch("dphij2b",&dphij2b,"dphij2b/F");
  
  outTree->Branch("drj1b",  &drj1b,  "drj1b/F");
  outTree->Branch("drj2b",  &drj2b,  "drj2b/F");
  outTree->Branch("bjcsv",  &bjcsv,  "bjcsv/F");
  outTree->Branch("j1csv",  &j1csv,  "j1csv/F");
  outTree->Branch("j2csv",  &j2csv,  "j2csv/F");
  outTree->Branch("weight", &weight, "weight/F");
  outTree->Branch("jet1", "TLorentzVector", &jet1);
  outTree->Branch("jet2", "TLorentzVector", &jet2);
  outTree->Branch("jet3", "TLorentzVector", &jet3);
  outTree->Branch("jet4", "TLorentzVector", &jet4);
  
  outTree->Branch("bjet", "TLorentzVector", &bjet);
  // 
  // Declare histograms
  //
  vector<TH1D*> hTopMassCombov;
  TH1D* hTopMass;
  char hname[100];
  hTopMass = new TH1D("hTopMass","",50,100,300);
  for(unsigned int i=0; i<4; ++i){
    sprintf(hname,"hTopMassCombo_%i",i); hTopMassCombov.push_back(new TH1D(hname,"",50,100,300)); hTopMassCombov[i]->Sumw2();
  }
  
  // set up pile-up reweighting
  TFile puWeightFile(puWeightFilename.c_str());
  TH1D *hPUWeights = (TH1D*)puWeightFile.Get("pileup");
  hPUWeights->SetDirectory(0);
  puWeightFile.Close();

  // set up MC smearing corrections for q/g discriminator
  QGSyst qgsyst;
  qgsyst.ReadDatabaseDoubleMin(qgidCorrFilename);

   // set up reader/book mva for WMVA
  int dummy; 
  float ptmjj, qjj, pull1ang, pull2ang, qgid1, qgid2, mdrop;
  TMVA::Reader *reader_lowpt, *reader_highpt;
  
   reader_lowpt = new TMVA::Reader("");
   reader_lowpt->AddSpectator("isSig",  &dummy);
   reader_lowpt->AddVariable("ptmjj",    &ptmjj);      
   reader_lowpt->AddVariable("qjj",      &qjj);      
   reader_lowpt->AddVariable("pull1ang", &pull1ang);      
   reader_lowpt->AddVariable("pull2ang", &pull2ang);      
   reader_lowpt->AddVariable("qgid1",    &qgid1);
   reader_lowpt->AddVariable("qgid2",    &qgid2);
   reader_lowpt->AddVariable("mdrop",    &mdrop);
   reader_lowpt->BookMVA("BDTG", weightfile_lowpt.c_str());

   reader_highpt = new TMVA::Reader("");
   reader_highpt->AddSpectator("isSig",  &dummy);
   reader_highpt->AddVariable("ptmjj",    &ptmjj);
   reader_highpt->AddVariable("qjj",      &qjj);
   reader_highpt->AddVariable("pull1ang", &pull1ang);
   reader_highpt->AddVariable("pull2ang", &pull2ang);
   reader_highpt->AddVariable("qgid1",    &qgid1);
   reader_highpt->AddVariable("qgid2",    &qgid2);
   reader_highpt->AddVariable("mdrop",    &mdrop);
   reader_highpt->BookMVA("BDTG", weightfile_highpt.c_str());
   

  TFile *infile=0;
  TTree *intree=0;
  
  cout << " ==> Processing " << infilename << "..." << endl;
  infile = new TFile(infilename.c_str()); assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
  
  intree->SetBranchAddress("metfilter",    &metfilter);
  intree->SetBranchAddress("npv",          &npv);
  intree->SetBranchAddress("npu",          &npu);
  intree->SetBranchAddress("njets",        &njets);
  intree->SetBranchAddress("nbjets",       &nbjets);
  // intree->SetBranchAddress("rho",          &rho);
  intree->SetBranchAddress("scale1fb",     &scale1fb);
  intree->SetBranchAddress("ttWeight",     &ttWeight);
  intree->SetBranchAddress("lepId",        &lepId);
  intree->SetBranchAddress("lepton",          &lepton);

  intree->SetBranchAddress("jet1flavGen",  &jet1flavGen);
  intree->SetBranchAddress("jet1qgid",     &jet1qgid);
  intree->SetBranchAddress("jet1csv",      &jet1csv);
  intree->SetBranchAddress("jet1flavAlgo", &jet1flavAlgo);
  intree->SetBranchAddress("jet1q",        &jet1q);
  intree->SetBranchAddress("jet1",         &jet1);
  intree->SetBranchAddress("jet1pull",     &jet1pull);
  
  intree->SetBranchAddress("jet2flavGen",  &jet2flavGen);
  intree->SetBranchAddress("jet2qgid",     &jet2qgid);
  intree->SetBranchAddress("jet2csv",      &jet2csv);
  intree->SetBranchAddress("jet2flavAlgo", &jet2flavAlgo);
  intree->SetBranchAddress("jet2q",        &jet2q);
  intree->SetBranchAddress("jet2",         &jet2);
  intree->SetBranchAddress("jet2pull",     &jet2pull);

  intree->SetBranchAddress("jet3flavGen",  &jet3flavGen);
  intree->SetBranchAddress("jet3qgid",     &jet3qgid);
  intree->SetBranchAddress("jet3csv",      &jet3csv);
  intree->SetBranchAddress("jet3flavAlgo", &jet3flavAlgo);
  intree->SetBranchAddress("jet3q",        &jet3q);
  intree->SetBranchAddress("jet3",         &jet3);
  intree->SetBranchAddress("jet3pull",     &jet3pull);

  intree->SetBranchAddress("jet4flavGen",  &jet4flavGen);
  intree->SetBranchAddress("jet4qgid",     &jet4qgid);
  intree->SetBranchAddress("jet4csv",      &jet4csv);
  intree->SetBranchAddress("jet4flavAlgo", &jet4flavAlgo);
  intree->SetBranchAddress("jet4q",        &jet4q);
  intree->SetBranchAddress("jet4",         &jet4);
  intree->SetBranchAddress("jet4pull",     &jet4pull);

  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(metfilter!=0)        continue;
    if(abs(lepId)!=13)      continue;
    if(njets  < NJETS_CUT)  continue;
            
    // QGL correction 
    int iflav1=0, iflav2=0, iflav3=0, iflav4=0;
    string sflav1="", sflav2="", sflav3="", sflav4="";
    if     (jet1flavAlgo==0)	   { iflav1 = 1; }		      // unmatched
    else if(abs(jet1flavAlgo) < 4) { iflav1 = 5; sflav1 = "quark"; }  // u,d,s
    else if(abs(jet1flavAlgo)== 4) { iflav1 = 4; sflav1 = "quark"; }  // c
    else if(abs(jet1flavAlgo)== 5) { iflav1 = 3; sflav1 = "quark"; }  // b
    else if(abs(jet1flavAlgo)==21) { iflav1 = 2; sflav1 = "gluon"; }  // gluon
    
    if     (jet2flavAlgo==0)	   { iflav2 = 1; }		      // unmatched
    else if(abs(jet2flavAlgo) < 4) { iflav2 = 5; sflav2 = "quark"; }  // u,d,s
    else if(abs(jet2flavAlgo)== 4) { iflav2 = 4; sflav2 = "quark"; }  // c
    else if(abs(jet2flavAlgo)== 5) { iflav2 = 3; sflav2 = "quark"; }  // b
    else if(abs(jet2flavAlgo)==21) { iflav2 = 2; sflav2 = "gluon"; }  // gluon

    if     (jet3flavAlgo==0)	   { iflav3 = 1; }		      // unmatched
    else if(abs(jet3flavAlgo) < 4) { iflav3 = 5; sflav3 = "quark"; }  // u,d,s
    else if(abs(jet3flavAlgo)== 4) { iflav3 = 4; sflav3 = "quark"; }  // c
    else if(abs(jet3flavAlgo)== 5) { iflav3 = 3; sflav3 = "quark"; }  // b
    else if(abs(jet3flavAlgo)==21) { iflav3 = 2; sflav3 = "gluon"; }  // gluon

    if     (jet4flavAlgo==0)	   { iflav4 = 1; }		      // unmatched
    else if(abs(jet4flavAlgo) < 4) { iflav4 = 5; sflav4 = "quark"; }  // u,d,s
    else if(abs(jet4flavAlgo)== 4) { iflav4 = 4; sflav4 = "quark"; }  // c
    else if(abs(jet4flavAlgo)== 5) { iflav4 = 3; sflav4 = "quark"; }  // b
    else if(abs(jet4flavAlgo)==21) { iflav4 = 2; sflav4 = "gluon"; }  // gluon

    if(iflav1!=1) jet1qgid = qgsyst.Smear(jet1->Pt(), jet1->Eta(), rho, jet1qgid, sflav1);
    if(iflav2!=1) jet2qgid = qgsyst.Smear(jet2->Pt(), jet2->Eta(), rho, jet2qgid, sflav2);	      
    if(iflav3!=1) jet3qgid = qgsyst.Smear(jet3->Pt(), jet3->Eta(), rho, jet3qgid, sflav3);	  
    if(iflav4!=1) jet4qgid = qgsyst.Smear(jet4->Pt(), jet4->Eta(), rho, jet4qgid, sflav4);	  

    double wgt = 1;
    wgt *= ttWeight;
    wgt *= hPUWeights->GetBinContent(hPUWeights->FindBin(npu));
    
    if( jet1csv ==-1 || jet2csv == -1 || jet3csv == -1 || jet4csv == -1 ) continue;
 
    vector <float> vjetcsv;
    vector < pair <TLorentzVector,float> > vjetinfo;
    vector < pair <float,float> > vjetGen;
    vector < pair <float,TVector2 > > vjetpull;    
    vjetcsv.push_back(jet1csv);
    vjetcsv.push_back(jet2csv);
    vjetcsv.push_back(jet3csv);
    vjetcsv.push_back(jet4csv);

    vjetinfo.push_back(make_pair(*jet1,jet1csv));
    vjetinfo.push_back(make_pair(*jet2,jet2csv));
    vjetinfo.push_back(make_pair(*jet3,jet3csv));
    vjetinfo.push_back(make_pair(*jet4,jet4csv));
    
    vjetpull.push_back(make_pair(jet1csv,*jet1pull));
    vjetpull.push_back(make_pair(jet2csv,*jet2pull));
    vjetpull.push_back(make_pair(jet3csv,*jet3pull));
    vjetpull.push_back(make_pair(jet4csv,*jet4pull));

    vector <float> vjetarr[4];
    vjetarr[0].push_back(jet1csv);     vjetarr[0].push_back(jet1flavGen);     vjetarr[0].push_back(jet1q);     vjetarr[0].push_back(jet1qgid);
    vjetarr[1].push_back(jet2csv);     vjetarr[1].push_back(jet2flavGen);     vjetarr[1].push_back(jet2q);     vjetarr[1].push_back(jet2qgid);
    vjetarr[2].push_back(jet3csv);     vjetarr[2].push_back(jet3flavGen);     vjetarr[2].push_back(jet3q);     vjetarr[2].push_back(jet3qgid);
    vjetarr[3].push_back(jet4csv);     vjetarr[3].push_back(jet4flavGen);     vjetarr[3].push_back(jet4q);     vjetarr[3].push_back(jet4qgid);

    // //---------debug----------
    // vector <int> MyVec;
    // MyVec.push_back(2);     MyVec.push_back(4);     MyVec.push_back(6);     MyVec.push_back(8); 
    // cout << "Debug: Test next_combinations" << endl;
    // do {
    //   for(unsigned int i=0; i < 3; i++){
    // 	cout << "My Vec " << i << ": " << MyVec[i] << endl; 
    //   }
    // }while(next_combination(MyVec.begin(), MyVec.begin()+3, MyVec.end()));
    // cout << "Done debug test" << endl;
    // //----------------------
    
    // cout << "start combinations" << endl;
    // consider all possible 3-jet combinations
    float max_csv = 0;
    TLorentzVector vbjet;
    TLorentzVector j1, j2;
    float bjetcsv, J1csv, J2csv, bjetflav, j1flav, j2flav, j1q, j2q, bjetq, j1qgid, j2qgid, bjetqgid;
    TVector2 j1pull, j2pull, bjetpull;
    
    do{
      max_csv = maximum(vjetcsv[0],vjetcsv[1],vjetcsv[2]);
      sort(vjetcsv.begin(), vjetcsv.begin()+3);
      
      for(unsigned int j=0; j < vjetcsv.size(); ++j)
       	{
      	  if (vjetcsv[0] == vjetinfo[j].second && vjetcsv[0] == vjetpull[j].first && vjetcsv[0] == vjetarr[j][0]) {
      	    j1 = vjetinfo[j].first;
	    J1csv = vjetcsv[0];
      	    j1flav = vjetarr[j][1];
      	    j1q = vjetarr[j][2];
	    j1qgid = vjetarr[j][3];
      	    j1pull = vjetpull[j].second;
	  }
      	  else if (vjetcsv[1] == vjetinfo[j].second && vjetcsv[1] == vjetpull[j].first &&  vjetcsv[1] == vjetarr[j][0]) {
      	    j2 = vjetinfo[j].first;
	    J2csv = vjetcsv[1];
      	    j2flav = vjetarr[j][1];
      	    j2q = vjetarr[j][2];
	    j2qgid = vjetarr[j][3];
      	    j2pull = vjetpull[j].second;
	  }
       	  else if (max_csv == vjetinfo[j].second && max_csv == vjetarr[j][0]) {
      	    vbjet = vjetinfo[j].first;
	    bjetcsv = max_csv;
      	    bjetflav = vjetarr[j][1];
      	    bjetq = vjetarr[j][2];
	    bjetqgid = vjetarr[j][3];
      	    bjetpull = vjetpull[j].second;
	  }
      	}
          
      //
      // extra cuts
      //
      if(fabs(j1.Eta())>2.4) continue;
      if(fabs(j2.Eta())>2.4) continue;
      
      // NOTE: assumes sample is semi-leptonic ttbar, so that only one top decays hadronically...
      
      bool wmatch = (j1flav!=0 && j2flav!=0               // jets have matches 
		     && j1flav!=j2flav                      // jets not matched to the same quark
		     && abs(j1flav)<5 && abs(j2flav)<5);    // jets not matched to b-quarks
      
      
      // change coordinates of jet axes for jet pull angle computation
      TVector2 jet2_in_jet1coord(j2.Rapidity() - j1.Rapidity(), deltaPhi(j2.Phi(),j1.Phi()));
      TVector2 jet1_in_jet2coord(j1.Rapidity() - j2.Rapidity(), deltaPhi(j1.Phi(),j2.Phi()));
      
      TLorentzVector dijet = (j1 + j2);
      if(dijet.M()<WMASSLOW || dijet.M()>WMASSHIGH) continue;
      
      // input variables for W MVA ID
      ptmjj    = dijet.Pt()/dijet.M();
      qjj      = j1q + j2q;
      pull1ang = j1pull.DeltaPhi(jet2_in_jet1coord);//jet1pull->DeltaPhi(jet2_in_jet1coord);
      pull2ang = j2pull.DeltaPhi(jet1_in_jet2coord);//jet2pull->DeltaPhi(jet1_in_jet2coord);
      qgid1    = j1qgid;//jet1qgid;
      qgid2    = j2qgid;//jet2qgid;
      mdrop    = TMath::Max(j1.M(), j2.M())/dijet.M()*(j1.DeltaR(j2));
      
      //
      // Fill output tree
      // 
      wmva   = (dijet.Pt()<160) ? reader_lowpt->EvaluateMVA("BDTG") : reader_highpt->EvaluateMVA("BDTG");
      weight = wgt;
      
      TLorentzVector vtop = dijet + vbjet;      
      
      if(vtop.M()>=TOPMASSLOW && vtop.M()<=TOPMASSHIGH) {
	bool bmatch = false;
	if( ((j1flav== 1 || j1flav== 3) && bjetflav==-5) ||
	    ((j1flav==-1 || j1flav==-3) && bjetflav== 5) ||
	    ((j1flav== 2 || j1flav== 4) && bjetflav== 5) ||
	    ((j1flav==-2 || j1flav==-4) && bjetflav==-5) 
	    ) {
	  bmatch = true;
	}

      if(wmatch && !bmatch){ hTopMassCombov[0]->Fill(vtop.M()); }
      else if(!wmatch && bmatch) { hTopMassCombov[1]->Fill(vtop.M()); }
      else if(!wmatch && !bmatch){ hTopMassCombov[2]->Fill(vtop.M()); }
      else if(wmatch && bmatch){ hTopMassCombov[3]->Fill(vtop.M()); }
      
      isSig   = wmatch && bmatch;
      b_mis   = wmatch && !bmatch;
      w_mis   = !wmatch && bmatch;
      wb_mis  = !wmatch && !bmatch;
      pttop   = vtop.Pt();
      mtop    = vtop.M();
      detaj1b = fabs(j1.Eta() - vbjet.Eta());
      detaj2b = fabs(j2.Eta() - vbjet.Eta());
      dphij1b = fabs(j1.DeltaPhi(vbjet));
      dphij2b = fabs(j2.DeltaPhi(vbjet));
      drj1b   = j1.DeltaR(vbjet);
      drj2b   = j2.DeltaR(vbjet);
      bjet    = &vbjet;
      bjcsv   = bjetcsv;
      j1csv   = J1csv;
      j2csv   = J2csv;
      outTree->Fill();
    }
    } while(next_combination(vjetcsv.begin(),vjetcsv.begin() + 3,vjetcsv.end()));      
    
  }
  
  delete infile;
  infile=0;
  intree=0; 

  //
  //Make plots
  //
  
  gStyle->SetTitleOffset(1.600,"Y");
  gStyle->SetPalette(1);
  TCanvas *c = MakeCanvas("c","c",800,600);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetGridx(1);
  c->SetGridy(1);
  
  float norm_TopM;
  norm_TopM = (1.0/hTopMass->Integral()); hTopMass->Scale(norm_TopM);

  vector<float> norm_TopMCombo;
  for(unsigned int k=0; k<4; ++k){
    norm_TopMCombo.push_back(1.0/hTopMassCombov[k]->Integral());
    hTopMassCombov[k]->Scale(norm_TopMCombo[k]);
  }

  char ylabel[100];
  sprintf(ylabel, "Fraction / %.2f",hTopMass->GetBinWidth(1));
  CPlot plotTopMass("TopMass","Top Mass","Mass [GeV]",ylabel);
  plotTopMass.AddHist1D(hTopMass,"","hist",kBlue,1,0);
  plotTopMass.Draw(c,true,format.c_str());
  
  sprintf(ylabel,"Fraction / %.2f",hTopMassCombov[0]->GetBinWidth(1));
  CPlot plotTopMassCombo("TopMassCombinations","","Mass [GeV]",ylabel);
  plotTopMassCombo.AddHist1D(hTopMassCombov[0],"B mismatch","hist",kGreen+2,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[1],"W mismatch","hist",kBlue,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[2],"B & W mismatch","hist",kOrange+7,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[3],"Signal Top","hist",kRed,1,0);
  plotTopMassCombo.TransLegend(-0.05,-0.02);
  plotTopMassCombo.Draw(c,true,format.c_str());

  outFile->Write();
  outFile->Close();
}
