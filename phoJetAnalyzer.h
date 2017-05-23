//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 19 15:38:26 2017 by ROOT version 6.06/01
// from TTree eventTree/event tree for analysis
// found on file: testOutput_data.root
//////////////////////////////////////////////////////////

#ifndef phoJetAnalyzer_h
#define phoJetAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <TLorentzVector.h>
#include <TPie.h>
using namespace std;
// Header file for the classes stored in the TTree if any.

class phoJetAnalyzer {
public :
   TH1F *weights_, *MC_distr_, *Data_distr_, *den, *hist;
   std::vector<float> MCpileup;
   std::vector<float> datapileup;

   TString pileFile;

// ***********************Input Tree****************************//

   TFile *fileName;
   TTree *tree;
   TTree *treeG;
   TTree *treeT;

   double event_weight;
   double neg_weight;
   bool isOverlap;
   int isHLT165;
   int isPrescaled;
   double TrigPhoton_pt;

   double Pho_pt;
   double Pho_eta;
   double Pho_phi;

   int NJets;
   int Jet_isLoose;
   double Jet_pt;
   double Jet_eta;
   double Jet_phi;

   bool Jet_isLoosePho;
   bool Jet_isMediumPho;
   bool Jet_isTightPho;

   double Jet_CHF;
   double Jet_CEF;
   double Jet_NHF;
   double Jet_NEF;
   int Jet_NCH;
   int Jet_NNH;

   int Jet_NConst;
   int Jet_NCharged;
   double Jet_area;
   double Jet_EtaWidth;
   double Jet_PhiWidth;
   double Jet_EtaPhiSpread;
   double Jet_ptDist;
   int Jet_n60;
   int Jet_n70;
   int Jet_n80;
   int Jet_n90;
   
   int Jet_NChargedHad;
   int Jet_NNeutralHad;
   int Jet_NPhoton;
   int Jet_NChargedHad10;
   int Jet_NNeutralHad10;
   int Jet_NPhoton10;
   int Jet_NChargedHad20;
   int Jet_NNeutralHad20;
   int Jet_NPhoton20;
   int Jet_NChargedHad50;
   int Jet_NNeutralHad50;
   int Jet_NPhoton50;
   int Jet_NChargedHad100;
   int Jet_NNeutralHad100;
   int Jet_NPhoton100;

   double PhoJet_dPhi;
   double PhoJet_dR;

   double JetGenPho_dR;
   double JetDarkPho_dR;
   double GenDarkPho_pt;
   double GenDarkPho_eta;
   double GenDarkPho_phi;
   double GenPhoton_pt;
   double GenPhoton_eta;
   double GenPhoton_phi;
   int GenDarkPho_isPho;
   int GenDarkPho_isJet;
   int GenDarkPho_isMet;
  
   double VertR; 
   double VertZ; 

// ***********************Uncertainty Tree****************************//
   double Jet_JEC_u;
   double Jet_JEC_d;
   double Jet_Res_u;
   double Jet_Res_d;
   double Pho_Sta_u;
   double Pho_Sta_d;
   double Pho_Sys_u;
   double Pho_Sys_d;
   double Pho_Gan_u;
   double Pho_Gan_d;
   double Pho_Res_u;
   double Pho_Res_d;

// ***********************Input Tree****************************//

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Int_t           nVtx;
   vector<float>   *vtxX;
   vector<float>   *vtxY;
   vector<float>   *vtxZ;
   vector<int>     *vtxNtrks;
   vector<bool>    *vtx_isFake;
   vector<int>     *vtx_ndof;
   vector<float>   *vtx_rho;
   vector<bool>    *isGoodVtx;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTPho;
   ULong64_t       HLTPhoIsPrescaled;
   vector<int>     *HLTPhoPrescale;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoPx;
   vector<float>   *phoPy;
   vector<float>   *phoPz;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoHoverE;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIEtaMapFull5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoIDMVA;
   vector<unsigned short> *phoIDbit;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<unsigned int> *phoFiredSingleTrgs;
   vector<unsigned int> *phoFiredDoubleTrgs;
   vector<unsigned int> *phoFiredL1Trgs;
   vector<float>   *phoScale_stat_up;
   vector<float>   *phoScale_stat_dn;
   vector<float>   *phoScale_syst_up;
   vector<float>   *phoScale_syst_dn;
   vector<float>   *phoScale_gain_up;
   vector<float>   *phoScale_gain_dn;
   vector<float>   *phoResol_rho_up;
   vector<float>   *phoResol_rho_dn;
   vector<float>   *phoResol_phi_up;
   vector<float>   *phoResol_phi_dn;
   Int_t           nJet;
   Int_t           jetTotalPt;
   Int_t           jetHT;
   Int_t           nGoodJets;
   vector<float>   *jetPt;
   vector<float>   *jetPx;
   vector<float>   *jetPy;
   vector<float>   *jetPz;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetMass;
   vector<float>   *jetMaxDistance;
   vector<float>   *jetPhiPhiMoment;
   vector<float>   *jetConstituentEtaPhiSpread;
   vector<float>   *jetConstituentPtDistribution;
   vector<float>   *jetPileup;
   vector<int>     *jetPartonID;
   vector<int>     *jetHadFlvr;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetLeadTrackEta;
   vector<float>   *jetLeadTrackPhi;
   vector<int>     *jetLepTrackPID;
   vector<float>   *jetLepTrackPt;
   vector<float>   *jetLepTrackEta;
   vector<float>   *jetLepTrackPhi;
   vector<float>   *jetChargedEmF;
   vector<float>   *jetNeutralEmF;
   vector<float>   *jetChargedHadF;
   vector<float>   *jetNeutralHadF;
   vector<int>     *jetPhotonEnF;
   vector<int>     *jetElectronEnF;
   vector<float>   *jetMuonEnF;
   vector<float>   *jetChargedMuonEnF;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConst;
   vector<int>     *jetNConstituents;
   vector<int>     *jetNConst60;
   vector<int>     *jetNConst70;
   vector<int>     *jetNConst80;
   vector<int>     *jetNConst90;
   vector<int>     *jetNConst92;
   vector<int>     *jetNConst94;
   vector<int>     *jetNConst96;
   vector<int>     *jetNConst98;
   vector<int>     *jetNCharged;
   vector<int>     *jetNNeutral;
   vector<int>     *jetNChargedHad;
   vector<int>     *jetNNeutralHad;
   vector<int>     *jetNPhoton;
   vector<int>     *jetNElectron;
   vector<int>     *jetNMuon;
   vector<float>   *jetetaWidth;
   vector<float>   *jetphiWidth;
   vector<float>   *jetPFCand12PtSum;
   vector<float>   *jetPFCandsPtSum;
   vector<float>   *jetPFCand12Ratio;
   vector<vector<double> > *jetConstEt;
   vector<vector<double> > *jetConstPt;
   vector<vector<int> > *jetConstPdgId;
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtxNtrks;
   vector<float>   *jetVtx3DVal;
   vector<float>   *jetVtx3DSig;
   vector<float>   *jetCSV2BJetTags;
   vector<float>   *jetJetProbabilityBJetTags;
   vector<float>   *jetpfCombinedMVAV2BJetTags;
   vector<bool>    *jetPFLooseId;
   vector<int>     *jetID;
   vector<float>   *jetPUID;
   vector<int>     *jetPUFullID;
   vector<float>   *jetJECUnc;
   vector<float>   *jetP4Smear;
   vector<float>   *jetP4SmearUp;
   vector<float>   *jetP4SmearDo;
   vector<float>   *jetGenJetEn;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenEn;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   Int_t           nEle;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleR9Full5x5;
   vector<float>   *eleEn;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaAtCalo;
   vector<float>   *eledEtaseedAtVtx;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *eleIDMVA;
   vector<unsigned short> *eleIDbit;
   vector<float>   *eleScale_stat_up;
   vector<float>   *eleScale_stat_dn;
   vector<float>   *eleScale_syst_up;
   vector<float>   *eleScale_syst_dn;
   vector<float>   *eleScale_gain_up;
   vector<float>   *eleScale_gain_dn;
   vector<float>   *eleResol_rho_up;
   vector<float>   *eleResol_rho_dn;
   vector<float>   *eleResol_phi_up;
   vector<float>   *eleResol_phi_dn;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<unsigned short> *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           metFilters;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         pdfWeight;
   vector<float>   *pdfSystWeight;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcStatus;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcIndex;
   Double_t        vertr;
   Double_t        vertz;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_vtxNtrks;   //!
   TBranch        *b_vtx_isFake;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_rho;   //!
   TBranch        *b_isGoodVtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTPhoPrescale;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIEtaMapFull5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoScale_stat_up;   //!
   TBranch        *b_phoScale_stat_dn;   //!
   TBranch        *b_phoScale_syst_up;   //!
   TBranch        *b_phoScale_syst_dn;   //!
   TBranch        *b_phoScale_gain_up;   //!
   TBranch        *b_phoScale_gain_dn;   //!
   TBranch        *b_phoResol_rho_up;   //!
   TBranch        *b_phoResol_rho_dn;   //!
   TBranch        *b_phoResol_phi_up;   //!
   TBranch        *b_phoResol_phi_dn;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTotalPt;   //!
   TBranch        *b_jetHT;   //!
   TBranch        *b_nGoodJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetPx;   //!
   TBranch        *b_jetPy;   //!
   TBranch        *b_jetPz;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetMaxDistance;   //!
   TBranch        *b_jetPhiPhiMoment;   //!
   TBranch        *b_jetConstituentEtaPhiSpread;   //!
   TBranch        *b_jetConstituentPtDistribution;   //!
   TBranch        *b_jetPileup;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetHadFlvr;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetLeadTrackEta;   //!
   TBranch        *b_jetLeadTrackPhi;   //!
   TBranch        *b_jetLepTrackPID;   //!
   TBranch        *b_jetLepTrackPt;   //!
   TBranch        *b_jetLepTrackEta;   //!
   TBranch        *b_jetLepTrackPhi;   //!
   TBranch        *b_jetChargedEmF;   //!
   TBranch        *b_jetNeutralEmF;   //!
   TBranch        *b_jetChargedHadF;   //!
   TBranch        *b_jetNeutralHadF;   //!
   TBranch        *b_jetPhotonEnF;   //!
   TBranch        *b_jetElectronEnF;   //!
   TBranch        *b_jetMuonEnF;   //!
   TBranch        *b_jetChargedMuonEnF;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConst;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetNConst60;   //!
   TBranch        *b_jetNConst70;   //!
   TBranch        *b_jetNConst80;   //!
   TBranch        *b_jetNConst90;   //!
   TBranch        *b_jetNConst92;   //!
   TBranch        *b_jetNConst94;   //!
   TBranch        *b_jetNConst96;   //!
   TBranch        *b_jetNConst98;   //!
   TBranch        *b_jetNCharged;   //!
   TBranch        *b_jetNNeutral;   //!
   TBranch        *b_jetNChargedHad;   //!
   TBranch        *b_jetNNeutralHad;   //!
   TBranch        *b_jetNPhoton;   //!
   TBranch        *b_jetNElectron;   //!
   TBranch        *b_jetNMuon;   //!
   TBranch        *b_jetetaWidth;   //!
   TBranch        *b_jetphiWidth;   //!
   TBranch        *b_jetPFCand12PtSum;   //!
   TBranch        *b_jetPFCandsPtSum;   //!
   TBranch        *b_jetPFCand12Ratio;   //!
   TBranch        *b_jetConstEt;   //!
   TBranch        *b_jetConstPt;   //!
   TBranch        *b_jetConstPdgId;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtxNtrks;   //!
   TBranch        *b_jetVtx3DVal;   //!
   TBranch        *b_jetVtx3DSig;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetpfCombinedMVAV2BJetTags;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetPUFullID;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetP4Smear;   //!
   TBranch        *b_jetP4SmearUp;   //!
   TBranch        *b_jetP4SmearDo;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eledEtaseedAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_eleIDMVA;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleScale_stat_up;   //!
   TBranch        *b_eleScale_stat_dn;   //!
   TBranch        *b_eleScale_syst_up;   //!
   TBranch        *b_eleScale_syst_dn;   //!
   TBranch        *b_eleScale_gain_up;   //!
   TBranch        *b_eleScale_gain_dn;   //!
   TBranch        *b_eleResol_rho_up;   //!
   TBranch        *b_eleResol_rho_dn;   //!
   TBranch        *b_eleResol_phi_up;   //!
   TBranch        *b_eleResol_phi_dn;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_pdfWeight;   //!
   TBranch        *b_pdfSystWeight;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_vertr;   //!
   TBranch        *b_vertz;   //!

   phoJetAnalyzer(const char* file1, const char* runner);
   virtual ~phoJetAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int isMonteCarlo);//,const char*pileupfile);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     BookHistos(const char* file2);
   bool to_bool(int atta);

   void SetSystematics();
   double passAnalysisCuts(int jetsystem, int phosystem);
   int NConstituents(int jet_index, int pdg, double pt);
   bool isPhoJetOverlap(int pho);
   bool OverlapWithElectron(double eta, double phi);
   bool OverlapWithMuon(double eta, double phi);
   double DeltaR(double eta1, double eta2, double phi1,double phi2);
   double DeltaPhi(double phi1,double phi2);
   bool MediumPhotonIdDecision(int &pho_index);
   void JetIsPhoton(int jet_index,int &jetpho_index);
   bool IsLoosePhoton(int &pho_index,double pt);
   bool IsMediumPhoton(int &pho_index,double pt);
   bool IsTightPhoton(int &pho_index,double pt);
   bool IsPixelPhoton(int &pho_index,double pt);
   bool VertexDecision(int &foundVertex);
   int JetDecision(int &pho_index, int &jet_index);
   void DarkPhotonType(double eta, double phi, int &isPho, int &isJet, int &isMet);
   double  EAElectroncharged(double eta);
   double  EAPFWorstElectroncharged(double eta);
   double  EAElectronneutral(double eta);
   double  EAElectronphoton(double eta);
   void LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr);
   double puweight(float npv);
};

#endif

#ifdef phoJetAnalyzer_cxx
phoJetAnalyzer::phoJetAnalyzer(const char* file1, const char* chrunner){
   gErrorIgnoreLevel = kError;
   pileFile = file1;
   TChain *chain = new TChain("phoJetNtuplizer/eventTree");
   int runner = atoi(chrunner);
   int fileNumber = 0;
   int runNumber =0;

   TString xrootder = "root://cmsxrootd.fnal.gov/";
   TString fileOut = file1;
   TString srunner = chrunner;
   TString infile = fileOut+".txt";
   ifstream in(infile);
   string fname;
   while (std::getline(in,fname)){
        if(fileNumber>-1 && runNumber==runner){
                TString FullPathInputFile = (xrootder+fname);
                std::cout<<FullPathInputFile<<std::endl;
                chain->Add(FullPathInputFile);
                std::cout<<chain->GetEntries()<<std::endl;
        }
        if(fileNumber%20==0) runNumber++;
        fileNumber++;
   }
   cout<<runNumber<<endl;
   Init(chain);
   cout<<"ADDED files:"<<endl;
   TString file2 = (fileOut+"_"+srunner+".root");
   BookHistos(file2);
}

phoJetAnalyzer::~phoJetAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   tree->Write();
   treeT->Write();
   treeG->Write();
   fileName->Close();
}

Int_t phoJetAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t phoJetAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void phoJetAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vtxX = 0;
   vtxY = 0;
   vtxZ = 0;
   vtxNtrks = 0;
   vtx_isFake = 0;
   vtx_ndof = 0;
   vtx_rho = 0;
   isGoodVtx = 0;
   HLTPhoPrescale = 0;
   phoE = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoR9Full5x5 = 0;
   phoHoverE = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIEtaMapFull5x5 = 0;
   phoPFChIso = 0;
   phoPFChWorstIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoIDMVA = 0;
   phoIDbit = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoScale_stat_up = 0;
   phoScale_stat_dn = 0;
   phoScale_syst_up = 0;
   phoScale_syst_dn = 0;
   phoScale_gain_up = 0;
   phoScale_gain_dn = 0;
   phoResol_rho_up = 0;
   phoResol_rho_dn = 0;
   phoResol_phi_up = 0;
   phoResol_phi_dn = 0;
   jetPt = 0;
   jetPx = 0;
   jetPy = 0;
   jetPz = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetMass = 0;
   jetMaxDistance = 0;
   jetPhiPhiMoment = 0;
   jetConstituentEtaPhiSpread = 0;
   jetConstituentPtDistribution = 0;
   jetPileup = 0;
   jetPartonID = 0;
   jetHadFlvr = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetChargedEmF = 0;
   jetNeutralEmF = 0;
   jetChargedHadF = 0;
   jetNeutralHadF = 0;
   jetPhotonEnF = 0;
   jetElectronEnF = 0;
   jetMuonEnF = 0;
   jetChargedMuonEnF = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConst = 0;
   jetNConstituents = 0;
   jetNConst60 = 0;
   jetNConst70 = 0;
   jetNConst80 = 0;
   jetNConst90 = 0;
   jetNConst92 = 0;
   jetNConst94 = 0;
   jetNConst96 = 0;
   jetNConst98 = 0;
   jetNCharged = 0;
   jetNNeutral = 0;
   jetNChargedHad = 0;
   jetNNeutralHad = 0;
   jetNPhoton = 0;
   jetNElectron = 0;
   jetNMuon = 0;
   jetetaWidth = 0;
   jetphiWidth = 0;
   jetPFCand12PtSum = 0;
   jetPFCandsPtSum = 0;
   jetPFCand12Ratio = 0;
   jetConstEt = 0;
   jetConstPt = 0;
   jetConstPdgId = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtxNtrks = 0;
   jetVtx3DVal = 0;
   jetVtx3DSig = 0;
   jetCSV2BJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetpfCombinedMVAV2BJetTags = 0;
   jetPFLooseId = 0;
   jetID = 0;
   jetPUID = 0;
   jetPUFullID = 0;
   jetJECUnc = 0;
   jetP4Smear = 0;
   jetP4SmearUp = 0;
   jetP4SmearDo = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleR9Full5x5 = 0;
   eleEn = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCRawEn = 0;
   eleSCEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eledEtaseedAtVtx = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   eleIDMVA = 0;
   eleIDbit = 0;
   eleScale_stat_up = 0;
   eleScale_stat_dn = 0;
   eleScale_syst_up = 0;
   eleScale_syst_dn = 0;
   eleScale_gain_up = 0;
   eleScale_gain_dn = 0;
   eleResol_rho_up = 0;
   eleResol_rho_dn = 0;
   eleResol_phi_up = 0;
   eleResol_phi_dn = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   pdf = 0;
   pdfSystWeight = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcStatus = 0;
   mcStatusFlag = 0;
   mcIndex = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("vtxNtrks", &vtxNtrks, &b_vtxNtrks);
   fChain->SetBranchAddress("vtx_isFake", &vtx_isFake, &b_vtx_isFake);
   fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_rho", &vtx_rho, &b_vtx_rho);
   fChain->SetBranchAddress("isGoodVtx", &isGoodVtx, &b_isGoodVtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTPhoPrescale", &HLTPhoPrescale, &b_HLTPhoPrescale);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaMapFull5x5", &phoSigmaIEtaIEtaMapFull5x5, &b_phoSigmaIEtaIEtaMapFull5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoScale_stat_up", &phoScale_stat_up, &b_phoScale_stat_up);
   fChain->SetBranchAddress("phoScale_stat_dn", &phoScale_stat_dn, &b_phoScale_stat_dn);
   fChain->SetBranchAddress("phoScale_syst_up", &phoScale_syst_up, &b_phoScale_syst_up);
   fChain->SetBranchAddress("phoScale_syst_dn", &phoScale_syst_dn, &b_phoScale_syst_dn);
   fChain->SetBranchAddress("phoScale_gain_up", &phoScale_gain_up, &b_phoScale_gain_up);
   fChain->SetBranchAddress("phoScale_gain_dn", &phoScale_gain_dn, &b_phoScale_gain_dn);
   fChain->SetBranchAddress("phoResol_rho_up", &phoResol_rho_up, &b_phoResol_rho_up);
   fChain->SetBranchAddress("phoResol_rho_dn", &phoResol_rho_dn, &b_phoResol_rho_dn);
   fChain->SetBranchAddress("phoResol_phi_up", &phoResol_phi_up, &b_phoResol_phi_up);
   fChain->SetBranchAddress("phoResol_phi_dn", &phoResol_phi_dn, &b_phoResol_phi_dn);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTotalPt", &jetTotalPt, &b_jetTotalPt);
   fChain->SetBranchAddress("jetHT", &jetHT, &b_jetHT);
   fChain->SetBranchAddress("nGoodJets", &nGoodJets, &b_nGoodJets);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetPx", &jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", &jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", &jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetMaxDistance", &jetMaxDistance, &b_jetMaxDistance);
   fChain->SetBranchAddress("jetPhiPhiMoment", &jetPhiPhiMoment, &b_jetPhiPhiMoment);
   fChain->SetBranchAddress("jetConstituentEtaPhiSpread", &jetConstituentEtaPhiSpread, &b_jetConstituentEtaPhiSpread);
   fChain->SetBranchAddress("jetConstituentPtDistribution", &jetConstituentPtDistribution, &b_jetConstituentPtDistribution);
   fChain->SetBranchAddress("jetPileup", &jetPileup, &b_jetPileup);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetHadFlvr", &jetHadFlvr, &b_jetHadFlvr);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetChargedEmF", &jetChargedEmF, &b_jetChargedEmF);
   fChain->SetBranchAddress("jetNeutralEmF", &jetNeutralEmF, &b_jetNeutralEmF);
   fChain->SetBranchAddress("jetChargedHadF", &jetChargedHadF, &b_jetChargedHadF);
   fChain->SetBranchAddress("jetNeutralHadF", &jetNeutralHadF, &b_jetNeutralHadF);
   fChain->SetBranchAddress("jetPhotonEnF", &jetPhotonEnF, &b_jetPhotonEnF);
   fChain->SetBranchAddress("jetElectronEnF", &jetElectronEnF, &b_jetElectronEnF);
   fChain->SetBranchAddress("jetMuonEnF", &jetMuonEnF, &b_jetMuonEnF);
   fChain->SetBranchAddress("jetChargedMuonEnF", &jetChargedMuonEnF, &b_jetChargedMuonEnF);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConst", &jetNConst, &b_jetNConst);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetNConst60", &jetNConst60, &b_jetNConst60);
   fChain->SetBranchAddress("jetNConst70", &jetNConst70, &b_jetNConst70);
   fChain->SetBranchAddress("jetNConst80", &jetNConst80, &b_jetNConst80);
   fChain->SetBranchAddress("jetNConst90", &jetNConst90, &b_jetNConst90);
   fChain->SetBranchAddress("jetNConst92", &jetNConst92, &b_jetNConst92);
   fChain->SetBranchAddress("jetNConst94", &jetNConst94, &b_jetNConst94);
   fChain->SetBranchAddress("jetNConst96", &jetNConst96, &b_jetNConst96);
   fChain->SetBranchAddress("jetNConst98", &jetNConst98, &b_jetNConst98);
   fChain->SetBranchAddress("jetNCharged", &jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetNNeutral", &jetNNeutral, &b_jetNNeutral);
   fChain->SetBranchAddress("jetNChargedHad", &jetNChargedHad, &b_jetNChargedHad);
   fChain->SetBranchAddress("jetNNeutralHad", &jetNNeutralHad, &b_jetNNeutralHad);
   fChain->SetBranchAddress("jetNPhoton", &jetNPhoton, &b_jetNPhoton);
   fChain->SetBranchAddress("jetNElectron", &jetNElectron, &b_jetNElectron);
   fChain->SetBranchAddress("jetNMuon", &jetNMuon, &b_jetNMuon);
   fChain->SetBranchAddress("jetetaWidth", &jetetaWidth, &b_jetetaWidth);
   fChain->SetBranchAddress("jetphiWidth", &jetphiWidth, &b_jetphiWidth);
   fChain->SetBranchAddress("jetPFCand12PtSum", &jetPFCand12PtSum, &b_jetPFCand12PtSum);
   fChain->SetBranchAddress("jetPFCandsPtSum", &jetPFCandsPtSum, &b_jetPFCandsPtSum);
   fChain->SetBranchAddress("jetPFCand12Ratio", &jetPFCand12Ratio, &b_jetPFCand12Ratio);
   fChain->SetBranchAddress("jetConstEt", &jetConstEt, &b_jetConstEt);
   fChain->SetBranchAddress("jetConstPt", &jetConstPt, &b_jetConstPt);
   fChain->SetBranchAddress("jetConstPdgId", &jetConstPdgId, &b_jetConstPdgId);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
   fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
   fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags, &b_jetpfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetPUFullID", &jetPUFullID, &b_jetPUFullID);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetP4Smear", &jetP4Smear, &b_jetP4Smear);
   fChain->SetBranchAddress("jetP4SmearUp", &jetP4SmearUp, &b_jetP4SmearUp);
   fChain->SetBranchAddress("jetP4SmearDo", &jetP4SmearDo, &b_jetP4SmearDo);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("eleIDMVA", &eleIDMVA, &b_eleIDMVA);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleScale_stat_up", &eleScale_stat_up, &b_eleScale_stat_up);
   fChain->SetBranchAddress("eleScale_stat_dn", &eleScale_stat_dn, &b_eleScale_stat_dn);
   fChain->SetBranchAddress("eleScale_syst_up", &eleScale_syst_up, &b_eleScale_syst_up);
   fChain->SetBranchAddress("eleScale_syst_dn", &eleScale_syst_dn, &b_eleScale_syst_dn);
   fChain->SetBranchAddress("eleScale_gain_up", &eleScale_gain_up, &b_eleScale_gain_up);
   fChain->SetBranchAddress("eleScale_gain_dn", &eleScale_gain_dn, &b_eleScale_gain_dn);
   fChain->SetBranchAddress("eleResol_rho_up", &eleResol_rho_up, &b_eleResol_rho_up);
   fChain->SetBranchAddress("eleResol_rho_dn", &eleResol_rho_dn, &b_eleResol_rho_dn);
   fChain->SetBranchAddress("eleResol_phi_up", &eleResol_phi_up, &b_eleResol_phi_up);
   fChain->SetBranchAddress("eleResol_phi_dn", &eleResol_phi_dn, &b_eleResol_phi_dn);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("pdfWeight", &pdfWeight, &b_pdfWeight);
   fChain->SetBranchAddress("pdfSystWeight", &pdfSystWeight, &b_pdfSystWeight);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("vertr", &vertr, &b_vertr);
   fChain->SetBranchAddress("vertz", &vertz, &b_vertz);
   Notify();
}

Bool_t phoJetAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void phoJetAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t phoJetAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef phoJetAnalyzer_cxx
