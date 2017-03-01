#define darkPhotonAnalyzer_cxx
#include "darkPhotonAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include "TLorentzVector.h"

using namespace std;

int main(int argc, const char* argv[]){
  const char* pileupfile = argv[1];
  int isMonteCarlo = atoi(argv[2]);
  darkPhotonAnalyzer t(argv[1],argv[3]);
  t.Loop(isMonteCarlo,pileupfile);
  return 0;
}

void darkPhotonAnalyzer::Loop(int isMonteCarlo, const char* pileupfile){
   pho_ptcut = 175.0;
   jet_ptcut = 170.0;
   track_ptcut = 20.0;
   mu_ptcut = 10.0;
   ele_ptcut = 10.0;
   bool isMC = to_bool(isMonteCarlo);
   TString pileFile = pileupfile;

   if(isMC){
      //TFile *f2 = TFile::Open("/uscms_data/d2/sushil/ForDusty/hist_DataPU.root");
      TFile *f2 = TFile::Open("/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/anaRoots/dataPileup.root");
      TH1F* dataPU = (TH1F*)f2->Get("hist");
      TString dir = "/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/forPileup/";
      TFile *fmc = TFile::Open(dir+pileupfile+".root");
      TH1F* mcPU = (TH1F*)fmc->Get("histo");
      for(int i=0; i<60; i++){
              MCpileup.push_back(mcPU->GetBinContent(i+1));
              datapileup.push_back(dataPU->GetBinContent(i+1));
      }
      LumiReWeighting(MCpileup,datapileup);
   }

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   GenDarkPho_isPho=0;
   GenDarkPho_isJet=0;
   GenDarkPho_isMet=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int pho_index = -1, jet_index = -1;
      int isPho = 0, isJet = 0, isMet = 0;
      int jetmultiplicity = 0;
      bool phoIdDecision = MediumPhotonIdDecision(pho_index);
      if(pileFile.Contains("df")){
	      TLorentzVector dp,hp;
	      dp.SetPxPyPzE(gen_DarkPho_px,gen_DarkPho_py,gen_DarkPho_pz,gen_DarkPho_E);
	      hp.SetPxPyPzE(gen_HardPho_px,gen_HardPho_py,gen_HardPho_pz,gen_HardPho_E);
	      GenDarkPho_eta = dp.Eta();
	      GenDarkPho_phi = dp.Phi();
	      GenDarkPho_pt = dp.Pt();
	      GenPhoton_eta = hp.Eta();
	      GenPhoton_pt = hp.Pt();
	      DarkPhotonType(dp.Eta(),dp.Phi(),isPho, isJet, isMet);
	      GenDarkPho_isPho+= isPho;
	      GenDarkPho_isJet+= isJet;
	      GenDarkPho_isMet+= isMet;
	      treeG->Fill();
      }
      neg_weight = 1.0;
      event_weight = 1.0;
      isOverlap = false;
      Jet_isLoosePho = false;
      Jet_isMediumPho= false;
      Jet_isTightPho = false;
      Jet_isPixelPho = false;
      if(isMC){
	      event_weight = puweight(puTrue->at(0));
	      if(pileFile == "dipho"){
		      fabs(genWeight) > 0.0 ? neg_weight *= genWeight/fabs(genWeight) : neg_weight = 0;
	      }
      	      hist->Fill(puTrue->at(0),(event_weight*neg_weight));
      }
      else{
	      if(HLTPhoton165_HE10 < 1) continue;
//	      if(HLTPhoton120 < 1 || HLTPhoton90 < 1 || HLTPhoton75 < 1) continue;
      }
      JetPho_dR=999.0;
      if(pho_index>=0) jetmultiplicity = JetDecision(pho_index,jet_index);
      if(phoIdDecision!=true) continue;
      isOverlap = isPhoJetOverlap(pho_index,pileFile);
      if(jet_index<0) continue;
      JetIsPhoton(jet_index,Jet_isLoosePho,Jet_isMediumPho,Jet_isTightPho,Jet_isPixelPho);

      if(pileFile == "dipho"){
      	for(int i=0;i<mcPID->size();i++){
//           if(mcPID->at(i)==22 && mcStatus->at(i)==1 && mcEt->at(i) > 100) cout<<"DeltaR: "<<DeltaR(mcEta->at(i),jetEta->at(jet_index),mcPhi->at(i),jetPhi->at(jet_index))<<endl;
           if(mcPID->at(i)==22 && mcStatus->at(i)==1 && mcEt->at(i) > 100 && JetPho_dR>DeltaR(mcEta->at(i),jetEta->at(jet_index),mcPhi->at(i),jetPhi->at(jet_index))){
              JetPho_dR = DeltaR(mcEta->at(i),jetEta->at(jet_index),mcPhi->at(i),jetPhi->at(jet_index));
           }
        }
      }

      Pho_pt  = phopt->at(pho_index);
      Pho_eta = phoeta->at(pho_index);
      Jet_isLoose = jetPFLooseId->at(jet_index);
      Jet_pt  = jetPt->at(jet_index);
      Jet_eta = jetEta->at(jet_index);
      Jet_NConst = jetNConstituents->at(jet_index);
      Jet_CHF = jetCHF->at(jet_index);
      Jet_CEF = jetCEF->at(jet_index);
      Jet_NCH = jetNCH->at(jet_index);
      Jet_NHF = jetNHF->at(jet_index);
      Jet_NEF = jetNEF->at(jet_index);
      Jet_NNH = jetNConstituents->at(jet_index)-jetNCH->at(jet_index);
      Jet_EtaWidth = jetEtaWidth->at(jet_index);
      Jet_PhiWidth = jetPhiWidth->at(jet_index);
      Jet_EtaWidthInECal = jetEtaWidthInECal->at(jet_index);
      Jet_PhiWidthInECal = jetPhiWidthInECal->at(jet_index);
      Jet_EtaWidthInHCal = jetEtaWidthInHCal->at(jet_index);
      Jet_PhiWidthInHCal = jetPhiWidthInHCal->at(jet_index);
      Jet_area = jetJetArea->at(jet_index);
      Jet_ptDist = jetConstituentPtDistribution->at(jet_index);
      Jet_n60 = jetN60->at(jet_index);
      Jet_n90 = jetN90->at(jet_index);
      Jet_PEF = jetPhotonEnergyFraction->at(jet_index);
      PhoJet_dPhi = TMath::Abs(DeltaPhi(phophi->at(pho_index),jetPhi->at(jet_index)));
      PhoJet_dR = DeltaR(phoeta->at(pho_index),jetEta->at(jet_index),phophi->at(pho_index),jetPhi->at(jet_index));
      NJets = jetmultiplicity;
      if(pileFile.Contains("df")) GenJetDarkPho_dR = DeltaR(GenDarkPho_eta,jetEta->at(jet_index),GenDarkPho_phi,jetPhi->at(jet_index));
//      JetDarkPho_dR = DeltaR(
//      if(NJets < 2 && Jet_NConst<=12 && Jet_CHF<=0.1){
//      	cout<<"run: "<<run<<"\tevent: "<<event<<"\tlumis: "<<lumis<<endl;
//      }
      tree->Fill();
   }
   treeN->Fill();
}

void darkPhotonAnalyzer::BookHistos(const char* file){
   Float_t PtBins[7]={145., 160., 190., 250., 400., 700.0,1000.0};
   Float_t MetBins[8]={130., 150., 170., 190., 250., 400., 700.0,1000.0};
   fileName = new TFile(file, "RECREATE");
   fileName->cd();
   tree = new TTree("darkPho","darkPho");
   tree->Branch("event_weight",&event_weight,"event_weight/D"); 
   tree->Branch("neg_weight",&neg_weight,"neg_weight/D");
   tree->Branch("isOverlap",&isOverlap,"isOverlap/B");
   tree->Branch("Jet_isLoosePho",&Jet_isLoosePho,"Jet_isLoosePho/B");
   tree->Branch("Jet_isMediumPho",&Jet_isMediumPho,"Jet_isMediumPho/B");
   tree->Branch("Jet_isTightPho",&Jet_isTightPho,"Jet_isTightPho/B");
   tree->Branch("Jet_isPixelPho",&Jet_isPixelPho,"Jet_isPixelPho/B");
   tree->Branch("Pho_pt",&Pho_pt,"Pho_pt/D"); 
   tree->Branch("Pho_eta",&Pho_eta,"Pho_eta/D");
   tree->Branch("Jet_isLoose",&Jet_isLoose,"Jet_isLoose/I");
   tree->Branch("Jet_pt",&Jet_pt,"Jet_pt/D");
   tree->Branch("Jet_eta",&Jet_eta,"Jet_eta/D");
   tree->Branch("Jet_CHF",&Jet_CHF,"Jet_CHF/D");
   tree->Branch("Jet_CEF",&Jet_CEF,"Jet_CEF/D");
   tree->Branch("Jet_NHF",&Jet_NHF,"Jet_NHF/D");
   tree->Branch("Jet_NEF",&Jet_NEF,"Jet_NEF/D");
   tree->Branch("Jet_NConst",&Jet_NConst,"Jet_NConst/I");
   tree->Branch("Jet_NCH",&Jet_NCH,"Jet_NCH/I");
   tree->Branch("Jet_NNH",&Jet_NNH,"Jet_NNH/I");
   tree->Branch("Jet_EtaWidth",&Jet_EtaWidth,"Jet_EtaWidth/D");
   tree->Branch("Jet_PhiWidth",&Jet_PhiWidth,"Jet_PhiWidth/D");
   tree->Branch("Jet_EtaWidthInECal",&Jet_EtaWidthInECal,"Jet_EtaWidthInECal/D");
   tree->Branch("Jet_PhiWidthInECal",&Jet_PhiWidthInECal,"Jet_PhiWidthInECal/D");
   tree->Branch("Jet_EtaWidthInHCal",&Jet_EtaWidthInHCal,"Jet_EtaWidthInHCal/D");
   tree->Branch("Jet_PhiWidthInHCal",&Jet_PhiWidthInHCal,"Jet_PhiWidthInHCal/D");
   tree->Branch("Jet_area",&Jet_area,"Jet_area/D");
   tree->Branch("Jet_ptDist",&Jet_ptDist,"Jet_ptDist/D");
   tree->Branch("Jet_n60",&Jet_n60,"Jet_n60/I");
   tree->Branch("Jet_n90",&Jet_n90,"Jet_n90/I");
   tree->Branch("Jet_PEF",&Jet_PEF,"Jet_PEF/D");
   tree->Branch("PhoJet_dPhi",&PhoJet_dPhi,"PhoJet_dPhi/D");
   tree->Branch("PhoJet_dR",&PhoJet_dR,"PhoJet_dR/D");
   tree->Branch("NJets",&NJets,"NJets/I");
   tree->Branch("GenDarkPho_pt",&GenDarkPho_pt,"GenDarkPho_pt/D");
   tree->Branch("GenDarkPho_eta",&GenDarkPho_eta,"GenDarkPho_eta/D");
   tree->Branch("GenPhoton_pt",&GenPhoton_pt,"GenPhoton_pt/D");
   tree->Branch("GenPhoton_eta",&GenPhoton_eta,"GenPhoton_eta/D");
   tree->Branch("JetPho_dR",&JetPho_dR,"JetPho_dR/D");
   tree->Branch("GenJetDarkPho_dR",&GenJetDarkPho_dR,"GenJetDarkPho_dR/D");
    
   treeN = new TTree("darkTypes","darkTypes");
   treeN->Branch("GenDarkPho_isPho",&GenDarkPho_isPho,"GenDarkPho_isPho/I");
   treeN->Branch("GenDarkPho_isJet",&GenDarkPho_isJet,"GenDarkPho_isJet/I");
   treeN->Branch("GenDarkPho_isMet",&GenDarkPho_isMet,"GenDarkPho_isMet/I");

   treeG = new TTree("darkGens","darkGens");
   treeG->Branch("GenDarkPho_pt",&GenDarkPho_pt,"GenDarkPho_pt/D");
   treeG->Branch("GenPhoton_pt",&GenPhoton_pt,"GenPhoton_pt/D");
   treeG->Branch("GenDarkPho_eta",&GenDarkPho_eta,"GenDarkPho_eta/D");
   treeG->Branch("GenPhoton_eta",&GenPhoton_eta,"GenPhoton_eta/D");


   hist = new TH1F("hist","Number of Primary Vertices",60,0.0,60.0);
}

bool darkPhotonAnalyzer::MediumPhotonIdDecision(int &pho_index){
   bool mediumPhotonID=false;
   for(int i=0; i< phopt->size()  ;i++) {
       if(IsMediumPhoton(i) && phopt->at(i)> pho_ptcut && fabs(phoeta->at(i))<1.4442){
	  mediumPhotonID = true;
	  pho_index = i;
	  break;
       }
   }
   return mediumPhotonID;  
}

bool darkPhotonAnalyzer::IsLoosePhoton(int &in){
   bool outVal = false;
   bool HoE = false;
   bool sieie = false;
   bool charged = false;
   bool neutral = false;
   bool photon = false;
   bool pixel = false;
//   bool ResSpikeCut     = (phoSigmaIEtaIEtaFull5x5->at(in)>0.001);
   if(phoeta->at(in)<1.4442){
	   HoE = (phohOverE->at(in)  < 0.0597);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01031);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 1.295 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 10.910+0.0148*phopt->at(in)+0.000017*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 3.630+0.0047*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   else{
	   HoE = (phohOverE->at(in)  < 0.0481);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03013);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 1.011 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 5.931+0.0163*phopt->at(in)+0.000014*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 6.641+0.0034*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   //cout<<"HoE: "<<HoE<<"\tSigmaIetaIeta: "<<sieie<<"\tcharged: "<<charged<<"\tneutral: "<<neutral<<"\tphoton: "<<photon<<"\tpixel: "<<pixel<<endl;
   if(HoE && sieie && charged && neutral && photon && pixel){
       outVal = true;
   }
   return outVal;
}

bool darkPhotonAnalyzer::IsMediumPhoton(int &in){
   bool outVal = false;
   bool HoE = false;
   bool sieie = false;
   bool charged = false;
   bool neutral = false;
   bool photon = false;
   bool pixel = false;
//   bool ResSpikeCut     = (phoSigmaIEtaIEtaFull5x5->at(in)>0.001);
   if(phoeta->at(in)<1.4442){
	   HoE = (phohOverE->at(in)  < 0.0396);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01022);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.441 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 2.725+0.0148*phopt->at(in)+0.000017*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 2.571+0.0047*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   else{
	   HoE = (phohOverE->at(in)  < 0.0209);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03001);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.442 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 1.715+0.0163*phopt->at(in)+0.000014*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 3.863+0.0034*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   //cout<<"HoE: "<<HoE<<"\tSigmaIetaIeta: "<<sieie<<"\tcharged: "<<charged<<"\tneutral: "<<neutral<<"\tphoton: "<<photon<<"\tpixel: "<<pixel<<endl;
   if(HoE && sieie && charged && neutral && photon && pixel){
       outVal = true;
   }
   return outVal;
}


bool darkPhotonAnalyzer::IsTightPhoton(int &in){
   bool outVal = false;
   bool HoE = false;
   bool sieie = false;
   bool charged = false;
   bool neutral = false;
   bool photon = false;
   bool pixel = false;
//   bool ResSpikeCut     = (phoSigmaIEtaIEtaFull5x5->at(in)>0.001);
   if(phoeta->at(in)<1.4442){
	   HoE = (phohOverE->at(in)  < 0.0269);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.00994);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.202 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 0.264+0.0148*phopt->at(in)+0.000017*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 2.362+0.0047*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   else{
	   HoE = (phohOverE->at(in)  < 0.0213);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03000);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.034 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 0.586+0.0163*phopt->at(in)+0.000014*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 2.617+0.0034*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 0);
   }
   //cout<<"HoE: "<<HoE<<"\tSigmaIetaIeta: "<<sieie<<"\tcharged: "<<charged<<"\tneutral: "<<neutral<<"\tphoton: "<<photon<<"\tpixel: "<<pixel<<endl;
   if(HoE && sieie && charged && neutral && photon && pixel){
       outVal = true;
   }
   return outVal;
}

bool darkPhotonAnalyzer::IsPixelPhoton(int &in){
   bool outVal = false;
   bool HoE = false;
   bool sieie = false;
   bool charged = false;
   bool neutral = false;
   bool photon = false;
   bool pixel = false;
//   bool ResSpikeCut     = (phoSigmaIEtaIEtaFull5x5->at(in)>0.001);
   if(phoeta->at(in)<1.4442){
	   HoE = (phohOverE->at(in)  < 0.0396);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01022);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.441 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 2.725+0.0148*phopt->at(in)+0.000017*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 2.571+0.0047*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 1);
   }
   else{
	   HoE = (phohOverE->at(in)  < 0.0209);
	   sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03001);
	   charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.442 );
	   neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 1.715+0.0163*phopt->at(in)+0.000014*phopt->at(in)*phopt->at(in) );
	   photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0)  < 3.863+0.0034*phopt->at(in) );
	   pixel =   (phohasPixelSeed->at(in)  == 1);
   }
   //cout<<"HoE: "<<HoE<<"\tSigmaIetaIeta: "<<sieie<<"\tcharged: "<<charged<<"\tneutral: "<<neutral<<"\tphoton: "<<photon<<"\tpixel: "<<pixel<<endl;
   if(HoE && sieie && charged && neutral && photon && pixel){
       outVal = true;
   }
   return outVal;
}

/*
bool darkPhotonAnalyzer::IsMediumPhoton(int &in){
   bool outVal = false;
//   bool ResSpikeCut     = (phoSigmaIEtaIEtaFull5x5->at(in)>0.001);
   bool HoE = (phohOverE->at(in)  < 0.0396);
   bool sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01022);
   bool charged =   (TMath::Max(((phoisoChargedHadrons->at(in))-rho*EAElectroncharged(phoeta->at(in))),0.0) < 0.441 );
   bool neutral =   (TMath::Max(((phoisoNeutralHadrons->at(in))-rho*EAElectronneutral(phoeta->at(in))),0.0) < 2.725+0.0148*phopt->at(in)+0.000017*phopt->at(in)*phopt->at(in) );
   bool photon =   (TMath::Max(((phoisoPhotons->at(in))       - rho*EAElectronphoton(phoeta->at(in))),0.0) < 2.571+0.0047*phopt->at(in) );
   bool pixel =   (phohasPixelSeed->at(in)  == 0);
   //cout<<"HoE: "<<HoE<<"\tSigmaIetaIeta: "<<sieie<<"\tcharged: "<<charged<<"\tneutral: "<<neutral<<"\tphoton: "<<photon<<"\tpixel: "<<pixel<<endl;
   if(HoE && sieie && charged && neutral && photon && pixel){
       outVal = true;
   }
   return outVal;
}
*/
int darkPhotonAnalyzer::JetDecision(int &pho_index, int &jet_index){
   int nJets = 0;
   for(int j=0;j<jetPt->size();j++){
//       if(jetPFLooseId->at(j) == 0) continue;	   
       if(jetPt->at(j) >30.0 && fabs(jetEta->at(j))<2.4) nJets++;
       if(OverlapWithMuon(jetEta->at(j),jetPhi->at(j))) continue;
       if(OverlapWithElectron(jetEta->at(j),jetPhi->at(j))) continue;
       if(DeltaR(jetEta->at(j),phoeta->at(pho_index),jetPhi->at(j),phophi->at(pho_index)) < 0.5) continue;
       if(jetPt->at(j) >30.0 && fabs(jetEta->at(j))<2.4){
//          nJets++;
	  if(jetPt->at(j)>jet_ptcut && jet_index<0){
		  jet_index = j;
	  }
       }
   }
   return nJets;
}

void darkPhotonAnalyzer::DarkPhotonType(double darketa, double darkphi, int &isPho, int &isJet, int &isMet){
   double dR;
   for(int p=0;p<phopt->size();p++){
       dR =DeltaR(darketa,phoeta->at(p),darkphi,phophi->at(p)); 
       if(dR<0.5){
          isPho=1;
	  break;
       }
   }
   if(isPho<1){
       for(int j=0;j<jetPt->size();j++){
	  dR = DeltaR(darketa,jetEta->at(j),darkphi,jetPhi->at(j));
          if(dR<0.5){
            isJet=1;
	    break;
          }
       }
       if(isJet<1) isMet=1;
   }
}

bool darkPhotonAnalyzer::OverlapWithElectron(double eta, double phi){
  bool overlap = false;
  for(int i=0;i<elePt->size();++i){
    if(ElectronPassVetoID->at(i) == true && elePt->at(i)>ele_ptcut){
      float dRJetEle = DeltaR(eleEta->at(i),eta,elePhi->at(i),phi);
      if(dRJetEle<0.5){
        overlap = true;
        break;
      }
    }
  }
  return overlap;
}

bool darkPhotonAnalyzer::OverlapWithMuon(double eta, double phi){
  bool overlap = false;
  for(int k=0;k<muPt->size();++k){
    if(muIsLooseID->at(k) == true && muPt->at(k) >mu_ptcut){
      float dRJetMu = DeltaR(muEta->at(k),eta,muPhi->at(k),phi);
      if(dRJetMu<0.5) {
        overlap = true;
        break;
      }
    }
  }
  return overlap;
}

double darkPhotonAnalyzer::DeltaR(double eta1, double eta2, double phi1, double phi2){
   double dPhi = fabs(phi1-phi2);
   Double_t pi =3.141592654;
   Double_t twopi =2.0*pi;

   if(dPhi<0) dPhi=-dPhi;
   if(dPhi>=(2*pi-dPhi))dPhi= 2.0*pi-dPhi;

   double dEta = fabs(eta1-eta2);

   double DR =1.0;
   DR= pow((dPhi*dPhi + dEta*dEta),0.5);
   return DR;
}

double darkPhotonAnalyzer::DeltaPhi(double phi1, double phi2){
  double result = -999.;
  result=(phi1-phi2);
  if(result > M_PI) result -= 2*M_PI;
  if(result <= -M_PI) result += 2*M_PI;
  return result;
}

bool darkPhotonAnalyzer::to_bool(int s){
        bool rets = (s == 1);
        return rets;
}

void darkPhotonAnalyzer::JetIsPhoton(int jet_index,bool &isLoose,bool &isMedium,bool &isTight,bool &isPixel){
	isLoose = false;
	isMedium = false;
	isTight = false;
	isPixel = false;
	for(int g=0;g<phopt->size();g++){
//		cout<<"is a medium ID Photon: "<<IsMediumPhoton(g)<<"\tPhoton has energy: "<<phopt->at(g)<<"\tDeltaR to Jet: "<<DeltaR(phoeta->at(g),jetEta->at(jet_index),phophi->at(g),jetPhi->at(jet_index))<<endl;
		if( phopt->at(g)>140
		   && DeltaR(phoeta->at(g),jetEta->at(jet_index),phophi->at(g),jetPhi->at(jet_index))<0.1){
			if(IsLoosePhoton(g))  isLoose = true;
			if(IsMediumPhoton(g)) isMedium = true;
			if(IsTightPhoton(g))  isTight = true;
			if(IsPixelPhoton(g))  isPixel = true;
		}
	}
}

bool darkPhotonAnalyzer::isPhoJetOverlap(int pho, TString run){
	bool retval = false;
	int genPho=-1;
	if(!run.Contains("qcd")) return retval;
	for(int g=0;g<mcPID->size();g++){
		if(mcPID->at(g) !=22 
		   || DeltaR(mcEta->at(g),phoeta->at(pho),mcPhi->at(g),phophi->at(pho)) >= 0.1 
		   || fabs(mcPt->at(g)-phopt->at(pho))/mcPt->at(g) >= 0.1 
		   || mcStatus->at(g) != 1) continue;
		genPho=g;
		break;
	}
	if(genPho<0) return retval;
	for(int p=0;p<mcPID->size();p++){
		if((fabs(mcPID->at(p))<9 || fabs(mcPID->at(p)) == 21)  
		   && DeltaR(mcEta->at(p),mcEta->at(genPho),mcPhi->at(p),mcPhi->at(genPho)) > 0.05) retval = true;
	}
	return retval;
}

void darkPhotonAnalyzer::LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr){
   if( MC_distr.size() != Lumi_distr.size() ){

     std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
     return;

   }

   Int_t NBins = MC_distr.size();

     MC_distr_ = new TH1F("MC_distr","MC dist",NBins,0.0, 60.0);
     Data_distr_ = new TH1F("Data_distr","Data dist",NBins,0.0, 60.0);

     weights_ = new TH1F("luminumer","luminumer",NBins,0.0,60.0);
     den = new TH1F("lumidenom","lumidenom",NBins,0.0,60.0);

     for(int ibin = 1; ibin<NBins+1; ++ibin ) {
       weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
       Data_distr_->SetBinContent(ibin, Lumi_distr[ibin-1]);
       den->SetBinContent(ibin,MC_distr[ibin-1]);
       MC_distr_->SetBinContent(ibin,MC_distr[ibin-1]);
     }

     float deltaH = weights_->Integral();
     if(fabs(1.0 - deltaH) > 0.02 ) {
       weights_->Scale( 1.0/ deltaH );
       Data_distr_->Scale( 1.0/ deltaH );
     }
     float deltaMC = den->Integral();
     if(fabs(1.0 - deltaMC) > 0.02 ) {
       den->Scale(1.0/ deltaMC );
       MC_distr_->Scale(1.0/ deltaMC );
     }
     weights_->Divide( den );
}                           
                           
                           
double darkPhotonAnalyzer::puweight(Float_t npv)
{                          
    Int_t bin = weights_->GetXaxis()->FindBin( npv );
    return weights_->GetBinContent( bin );
 }                       

double darkPhotonAnalyzer::EAPFWorstElectroncharged(double eta){
  float EffectiveArea=0.;
  if (fabs(eta) < 1.0 )   EffectiveArea = 0.075;
  if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.062;

  return EffectiveArea;
}

double darkPhotonAnalyzer::EAElectroncharged(double eta){
   float EffectiveArea=0.;
   if (fabs(eta) < 1.0 )   EffectiveArea = 0.0360;
   if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0377;

   return EffectiveArea;
}

double darkPhotonAnalyzer::EAElectronneutral(double eta){
   float EffectiveArea=0.;
   if (fabs(eta) < 1.0 )   EffectiveArea = 0.0597;
   if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0807;

   return EffectiveArea;
}

double darkPhotonAnalyzer::EAElectronphoton(double eta){
   float EffectiveArea=0.;
   if (fabs(eta) < 1.0 )   EffectiveArea = 0.1210;
   if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1107;

   return EffectiveArea;
}
