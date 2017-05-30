#define phoJetAnalyzer_cxx
#include "phoJetAnalyzer.h"
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

int main(int argc, const char* argv[]){
	phoJetAnalyzer t(argv[1],argv[2],argv[3]);
	t.Loop();
	return 0;
}

void phoJetAnalyzer::Loop(){
	//Set up pileup, scalefactor, and kfactor histograms
	if(isMC){
		TFile *f2 = TFile::Open("/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/Weights/forPileup/Data.root");
		TH1F* dataPU = (TH1F*)f2->Get("histo");
		TString dir = "/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/Weights/forPileup/";
		TFile *fmc = TFile::Open(dir+pileFile+".root");
		TH1F* mcPU = (TH1F*)fmc->Get("histo");
		for(int i=0; i<60; i++){
			MCpileup.push_back(mcPU->GetBinContent(i+1));
			datapileup.push_back(dataPU->GetBinContent(i+1));
		}
		LumiReWeighting(MCpileup,datapileup);
	}
	Sfac = TFile::Open("/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/Weights/ScaleFactors/egammaEffi.txt_EGM2D.root");
	Kfac = TFile::Open("/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/Weights/KFactors/KFactor.root");
	scaleHist= (TH2F*)Sfac->Get("EGamma_SF2D");
	kfacHist = (TH1F*)Kfac->Get("pt_hist");
	GenDarkPho_isPho=0;
	GenDarkPho_isJet=0;
	GenDarkPho_isMet=0;

	//For systematics
	Jet_JEC_u=0;
	Jet_JEC_d=0;
//	Jet_Res_u=0;
//	Jet_Res_d=0;
	Pho_Scale_u=0;
	Pho_Scale_d=0;
	Pho_IDScale_u=0;
	Pho_IDScale_d=0;
	Sig_kfac_u=0;
	Sig_kfac_d=0;
//	Pho_Gan_u=0;
//	Pho_Gan_d=0;
//	Pho_Res_u=0;
//	Pho_Res_d=0;
	vector<int> GenPhos;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		int pho_index=-1, jet_index=-1, jetpho_index=-1;
		int isPho=0,isJet=0,isMet=0;
		int jetmultiplicity=0;

		//For weights
		neg_weight = 1.0;
		event_weight = 1.0;
		isOverlap = 0;
		IDScaleFactor = 1.0;
		kfactor = 1.0;

      		JetGenPho_dR=999.0;
		isHLT165=0;
		isPrescaled=0;
		double jetPtSmear = 0.0;
		GenPhos.clear();
		TLorentzVector dp;
		NVerts = nVtx;
		if(isMC){
			for(int i=0;i<mcPID->size();i++){
				if(mcPID->at(i) == 1072000) dp.SetPtEtaPhiE(mcPt->at(i),mcEta->at(i),mcPhi->at(i),mcE->at(i));
				if(mcPID->at(i) == 22 && mcStatus->at(i) == 1 && mcPt->at(i) > 100){
					GenPhos.push_back(i);
				}
			}
			if(pileFile.Contains("df")){
				GenDarkPho_eta = dp.Eta();
				GenDarkPho_phi = dp.Phi();
				GenDarkPho_pt = dp.Pt();
				if(GenPhos.size()>0){
					GenPhoton_eta = mcEta->at(GenPhos.at(0));
					GenPhoton_pt = mcPt->at(GenPhos.at(0));
					GenPhoton_phi = mcPhi->at(GenPhos.at(0));
				}
				DarkPhotonType(dp.Eta(),dp.Phi(),isPho, isJet, isMet);
				GenDarkPho_isPho+= isPho;
				GenDarkPho_isJet+= isJet;
				GenDarkPho_isMet+= isMet;
				VertR = vertr;
				VertZ = fabs(vertz);
				kfactor = GetKFactor(GenDarkPho_pt,0);
			}
			jetPtSmear = 1.0;
			event_weight = puweight(puTrue->at(0));
//			if(pileFile=="DiPhoton"){
//				fabs(genWeight) > 0.0 ? neg_weight *= genWeight/fabs(genWeight) : neg_weight = 0;
//			}
			NVerts = puTrue->at(0);
		}

		//For Trigger turn on.
		if(HLTPho>>7 &1) isHLT165 = 1;
		if(HLTPho>>0 &1 
		   || HLTPho>>1 &1 
		   || HLTPho>>2 &1 
		   || HLTPho>>3 &1 
		   || HLTPho>>4 &1 
		   || HLTPho>>5 &1 
		   || HLTPho>>6 &1) isPrescaled = 1;
		(phoEt->size()<1) ? TrigPhoton_pt = 0.0 : TrigPhoton_pt = phoEt->at(0);
		treeT->Fill();
//		if(isMC) SetSystematics();
//		if(!isMC && HLTPho==false) continue;
		
		//Good Vertex requirement		
		if(!isGoodVtx) continue;
		//Photon cuts
		bool phoIdDecision = MediumPhotonIdDecision(pho_index);
	        if(phoIdDecision!=true) continue;
		//Jet cuts
            	jetmultiplicity = JetDecision(pho_index,jet_index,jetPtSmear);
		if(jet_index<0) continue;
		//find generated photon closest to jet
		for(int i=0;i<GenPhos.size();i++){
			double delR = DeltaR(mcEta->at(GenPhos.at(i)),jetEta->at(jet_index),mcPhi->at(GenPhos.at(i)),jetPhi->at(jet_index));
			if(JetGenPho_dR>delR) JetGenPho_dR=delR;
		}
		//find closest photon to jet
		JetIsPhoton(jet_index,jetpho_index);
      		if(jetpho_index>=0){
                	Jet_isLoosePho = IsLoosePhoton(jetpho_index,phoEt->at(jetpho_index));
                	Jet_isMediumPho = IsMediumPhoton(jetpho_index,phoEt->at(jetpho_index));
                	Jet_isTightPho = IsTightPhoton(jetpho_index,phoEt->at(jetpho_index));
		}

		//Set scalings
		if(isMC) IDScaleFactor = GetIDScaleFactor(pho_index,0);
		if(pileFile.Contains("QCD")) isOverlap = isPhoJetOverlap(pho_index);

		//Set vairables to store
		if(pileFile.Contains("df")) JetDarkPho_dR = DeltaR(GenDarkPho_eta,jetEta->at(jet_index),GenDarkPho_phi,jetPhi->at(jet_index));
		Pho_pt  = phoEt->at(pho_index);
		Pho_eta = phoEta->at(pho_index);
		Pho_phi = phoPhi->at(pho_index);
	
		NJets = jetmultiplicity;
		Jet_isLoose = jetPFLooseId->at(jet_index);
		Jet_pt = jetPtSmear;
		Jet_eta = jetEta->at(jet_index);
		Jet_phi = jetPhi->at(jet_index);
		Jet_CHF = jetChargedHadF->at(jet_index);
		Jet_CEF = jetChargedEmF->at(jet_index);
		Jet_NHF = jetNeutralHadF->at(jet_index);
		Jet_NEF = jetNeutralEmF->at(jet_index);
		Jet_NCH = jetNChargedHad->at(jet_index);
		Jet_NNH = jetNNeutralHad->at(jet_index);

		Jet_NConst = jetNConstituents->at(jet_index);
		Jet_NCharged = jetNConstituents->at(jet_index);
		Jet_area = jetArea->at(jet_index);
		Jet_EtaWidth = jetetaWidth->at(jet_index);
		Jet_PhiWidth = jetphiWidth->at(jet_index);
		Jet_EtaPhiSpread = jetConstituentEtaPhiSpread->at(jet_index);
		Jet_ptDist = jetConstituentPtDistribution->at(jet_index);
		Jet_n60 = jetNConst60->at(jet_index);
		Jet_n70 = jetNConst70->at(jet_index);
		Jet_n80 = jetNConst80->at(jet_index);
		Jet_n90 = jetNConst90->at(jet_index);

		Jet_NChargedHad    = jetNChargedHad->at(jet_index);
		Jet_NNeutralHad    = jetNNeutralHad->at(jet_index);
		Jet_NPhoton        = jetNPhoton->at(jet_index);
		Jet_NChargedHad0   = NConstituents(jet_index,211,0.);
		Jet_NNeutralHad0   = NConstituents(jet_index,130,0.);
		Jet_NPhoton0       = NConstituents(jet_index,22,0.);
		Jet_NChargedHad1   = NConstituents(jet_index,211,1.);
		Jet_NNeutralHad1   = NConstituents(jet_index,130,1.);
		Jet_NPhoton1       = NConstituents(jet_index,22,1.);
		Jet_NChargedHad2   = NConstituents(jet_index,211,2.);
		Jet_NNeutralHad2   = NConstituents(jet_index,130,2.);
		Jet_NPhoton2       = NConstituents(jet_index,22,2.);
		Jet_NChargedHad5   = NConstituents(jet_index,211,5.);
		Jet_NNeutralHad5   = NConstituents(jet_index,130,5.);
		Jet_NPhoton5       = NConstituents(jet_index,22,5.);
		Jet_NChargedHad10  = NConstituents(jet_index,211,10.);
		Jet_NNeutralHad10  = NConstituents(jet_index,130,10.);
		Jet_NPhoton10      = NConstituents(jet_index,22,10.);
		Jet_NChargedHad20  = NConstituents(jet_index,211,20.);
		Jet_NNeutralHad20  = NConstituents(jet_index,130,20.);
		Jet_NPhoton20      = NConstituents(jet_index,22,20.);
		Jet_NChargedHad50  = NConstituents(jet_index,211,50.);
		Jet_NNeutralHad50  = NConstituents(jet_index,130,50.);
		Jet_NPhoton50      = NConstituents(jet_index,22,50.);
		Jet_NChargedHad100 = NConstituents(jet_index,211,100.);
		Jet_NNeutralHad100 = NConstituents(jet_index,130,100.);
		Jet_NPhoton100     = NConstituents(jet_index,22,100.);

		Jet_NConstituents0   = NConstituents(jet_index,-1,0.);
		Jet_NConstituents1   = NConstituents(jet_index,-1,1.);
		Jet_NConstituents2   = NConstituents(jet_index,-1,2.);
		Jet_NConstituents5   = NConstituents(jet_index,-1,5.);
		Jet_NConstituents10  = NConstituents(jet_index,-1,10.);
		Jet_NConstituents20  = NConstituents(jet_index,-1,20.);
		Jet_NConstituents50  = NConstituents(jet_index,-1,50.);
		Jet_NConstituents100 = NConstituents(jet_index,-1,100.);
		Jet_NConstituents120 = NConstituents(jet_index,-1,120.);
		Jet_NConstituents140 = NConstituents(jet_index,-1,140.);
		Jet_NConstituents150 = NConstituents(jet_index,-1,150.);
		Jet_NConstituents160 = NConstituents(jet_index,-1,160.);
		Jet_NConstituents170 = NConstituents(jet_index,-1,170.);

		Jet_NConst70 = NConstituents(jet_index,-1,0.7*jetPt->at(jet_index));
		Jet_NConst80 = NConstituents(jet_index,-1,0.8*jetPt->at(jet_index));
		Jet_NConst82 = NConstituents(jet_index,-1,0.82*jetPt->at(jet_index));
		Jet_NConst85 = NConstituents(jet_index,-1,0.85*jetPt->at(jet_index));
		Jet_NConst87 = NConstituents(jet_index,-1,0.87*jetPt->at(jet_index));
		Jet_NConst90 = NConstituents(jet_index,-1,0.9*jetPt->at(jet_index));
		Jet_NConst95 = NConstituents(jet_index,-1,0.95*jetPt->at(jet_index));

		PhoJet_dPhi = fabs(DeltaPhi(jetPhi->at(jet_index),phoPhi->at(pho_index)));
		PhoJet_dR =  DeltaR(jetEta->at(jet_index),phoEta->at(pho_index),jetPhi->at(jet_index),phoPhi->at(pho_index));

		tree->Fill();
	}
	treeG->Fill();
}

void phoJetAnalyzer::BookHistos(const char* file){
	fileName = new TFile(file, "RECREATE");
	tree = new TTree("darkPho","darkPho");
	tree->Branch("event_weight",&event_weight,"event_weight/D");
	tree->Branch("neg_weight",&neg_weight,"neg_weight/D");
	tree->Branch("IDScaleFactor",&IDScaleFactor,"IDScaleFactor/D");
	tree->Branch("kfactor",&kfactor,"kfactor/D");
	tree->Branch("isOverlap",&isOverlap,"isOverlap/I");
	tree->Branch("Pho_pt",&Pho_pt,"Pho_pt/D");
	tree->Branch("Pho_phi",&Pho_phi,"Pho_phi/D");
	tree->Branch("Pho_eta",&Pho_eta,"Pho_eta/D");
	tree->Branch("Jet_isLoose",&Jet_isLoose,"Jet_isLoose/I");
	tree->Branch("Jet_isLoosePho",&Jet_isLoosePho,"Jet_isLoosePho/B");
	tree->Branch("Jet_isMediumPho",&Jet_isMediumPho,"Jet_isMediumPho/B");
	tree->Branch("Jet_isTightPho",&Jet_isTightPho,"Jet_isTightPho/B");
	tree->Branch("Jet_pt",&Jet_pt,"Jet_pt/D");
	tree->Branch("Jet_phi",&Jet_phi,"Jet_phi/D");
	tree->Branch("Jet_eta",&Jet_eta,"Jet_eta/D");
	tree->Branch("Jet_CHF",&Jet_CHF,"Jet_CHF/D");
	tree->Branch("Jet_CEF",&Jet_CEF,"Jet_CEF/D");
	tree->Branch("Jet_NHF",&Jet_NHF,"Jet_NHF/D");
	tree->Branch("Jet_NEF",&Jet_NEF,"Jet_NEF/D");
	tree->Branch("Jet_NCharged",&Jet_NCharged,"Jet_NCharged/I");
	tree->Branch("Jet_NConst",&Jet_NConst,"Jet_NConst/I");
	tree->Branch("Jet_NCH",&Jet_NCH,"Jet_NCH/I");
	tree->Branch("Jet_NNH",&Jet_NNH,"Jet_NNH/I");
	tree->Branch("Jet_EtaWidth",&Jet_EtaWidth,"Jet_EtaWidth/D");
	tree->Branch("Jet_PhiWidth",&Jet_PhiWidth,"Jet_PhiWidth/D");
	tree->Branch("Jet_area",&Jet_area,"Jet_area/D");
	tree->Branch("Jet_ptDist",&Jet_ptDist,"Jet_ptDist/D");
	tree->Branch("Jet_n60",&Jet_n60,"Jet_n60/I");
	tree->Branch("Jet_n70",&Jet_n70,"Jet_n70/I");
	tree->Branch("Jet_n80",&Jet_n80,"Jet_n80/I");
	tree->Branch("Jet_n90",&Jet_n90,"Jet_n90/I");
	tree->Branch("PhoJet_dPhi",&PhoJet_dPhi,"PhoJet_dPhi/D");
	tree->Branch("PhoJet_dR",&PhoJet_dR,"PhoJet_dR/D");
	tree->Branch("JetGenPho_dR",&JetGenPho_dR,"JetGenPho_dR/D");
	tree->Branch("NJets",&NJets,"NJets/I");

	tree->Branch("Jet_NConstituents0",&Jet_NConstituents0,"Jet_NConstituents0/I");
	tree->Branch("Jet_NConstituents1",&Jet_NConstituents1,"Jet_NConstituents1/I");
	tree->Branch("Jet_NConstituents2",&Jet_NConstituents2,"Jet_NConstituents2/I");
	tree->Branch("Jet_NConstituents5",&Jet_NConstituents5,"Jet_NConstituents5/I");
	tree->Branch("Jet_NConstituents10",&Jet_NConstituents10,"Jet_NConstituents10/I");
	tree->Branch("Jet_NConstituents20",&Jet_NConstituents20,"Jet_NConstituents20/I");
	tree->Branch("Jet_NConstituents50",&Jet_NConstituents50,"Jet_NConstituents50/I");
	tree->Branch("Jet_NConstituents100",&Jet_NConstituents100,"Jet_NConstituents100/I");
	tree->Branch("Jet_NConstituents120",&Jet_NConstituents120,"Jet_NConstituents120/I");
	tree->Branch("Jet_NConstituents140",&Jet_NConstituents140,"Jet_NConstituents140/I");
	tree->Branch("Jet_NConstituents150",&Jet_NConstituents150,"Jet_NConstituents150/I");
	tree->Branch("Jet_NConstituents160",&Jet_NConstituents160,"Jet_NConstituents160/I");
	tree->Branch("Jet_NConstituents170",&Jet_NConstituents170,"Jet_NConstituents170/I");

	tree->Branch("Jet_NConst70",&Jet_NConst70,"Jet_NConst70/I");
	tree->Branch("Jet_NConst80",&Jet_NConst80,"Jet_NConst80/I");
	tree->Branch("Jet_NConst82",&Jet_NConst82,"Jet_NConst82/I");
	tree->Branch("Jet_NConst85",&Jet_NConst85,"Jet_NConst85/I");
	tree->Branch("Jet_NConst87",&Jet_NConst87,"Jet_NConst87/I");
	tree->Branch("Jet_NConst90",&Jet_NConst90,"Jet_NConst90/I");
	tree->Branch("Jet_NConst95",&Jet_NConst95,"Jet_NConst95/I");

	tree->Branch("Jet_NChargedHad",&Jet_NChargedHad,"Jet_NChargedHad/I");
	tree->Branch("Jet_NNeutralHad",&Jet_NNeutralHad,"Jet_NNeutralHad/I");
	tree->Branch("Jet_NPhoton",&Jet_NPhoton,"Jet_NPhoton/I");
	tree->Branch("Jet_NChargedHad1",&Jet_NChargedHad1,"Jet_NChargedHad1/I");
	tree->Branch("Jet_NNeutralHad1",&Jet_NNeutralHad1,"Jet_NNeutralHad1/I");
	tree->Branch("Jet_NPhoton1",&Jet_NPhoton1,"Jet_NPhoton1/I");
	tree->Branch("Jet_NChargedHad2",&Jet_NChargedHad2,"Jet_NChargedHad2/I");
	tree->Branch("Jet_NNeutralHad2",&Jet_NNeutralHad2,"Jet_NNeutralHad2/I");
	tree->Branch("Jet_NPhoton2",&Jet_NPhoton2,"Jet_NPhoton2/I");
	tree->Branch("Jet_NChargedHad5",&Jet_NChargedHad5,"Jet_NChargedHad5/I");
	tree->Branch("Jet_NNeutralHad5",&Jet_NNeutralHad5,"Jet_NNeutralHad5/I");
	tree->Branch("Jet_NPhoton5",&Jet_NPhoton5,"Jet_NPhoton5/I");
	tree->Branch("Jet_NChargedHad10",&Jet_NChargedHad10,"Jet_NChargedHad10/I");
	tree->Branch("Jet_NNeutralHad10",&Jet_NNeutralHad10,"Jet_NNeutralHad10/I");
	tree->Branch("Jet_NPhoton10",&Jet_NPhoton10,"Jet_NPhoton10/I");
	tree->Branch("Jet_NChargedHad20",&Jet_NChargedHad20,"Jet_NChargedHad20/I");
	tree->Branch("Jet_NNeutralHad20",&Jet_NNeutralHad20,"Jet_NNeutralHad20/I");
	tree->Branch("Jet_NPhoton20",&Jet_NPhoton20,"Jet_NPhoton20/I");
	tree->Branch("Jet_NChargedHad50",&Jet_NChargedHad50,"Jet_NChargedHad50/I");
	tree->Branch("Jet_NNeutralHad50",&Jet_NNeutralHad50,"Jet_NNeutralHad50/I");
	tree->Branch("Jet_NPhoton50",&Jet_NPhoton50,"Jet_NPhoton50/I");
	tree->Branch("Jet_NChargedHad100",&Jet_NChargedHad100,"Jet_NChargedHad100/I");
	tree->Branch("Jet_NNeutralHad100",&Jet_NNeutralHad100,"Jet_NNeutralHad100/I");
	tree->Branch("Jet_NPhoton100",&Jet_NPhoton100,"Jet_NPhoton100/I");

	tree->Branch("GenDarkPho_pt",&GenDarkPho_pt,"GenDarkPho_pt/D");
	tree->Branch("GenDarkPho_eta",&GenDarkPho_eta,"GenDarkPho_eta/D");
	tree->Branch("GenPhoton_pt",&GenPhoton_pt,"GenPhoton_pt/D");
	tree->Branch("GenPhoton_eta",&GenPhoton_eta,"GenPhoton_eta/D");
	tree->Branch("JetDarkPho_dR",&JetDarkPho_dR,"JetDarkPho_dR/D");
	tree->Branch("isHLT165",&isHLT165,"isHLT165/I");

	treeG = new TTree("darkTypes","darkTypes");
	treeG->Branch("GenDarkPho_isPho",&GenDarkPho_isPho,"GenDarkPho_isPho/I");
	treeG->Branch("GenDarkPho_isJet",&GenDarkPho_isJet,"GenDarkPho_isJet/I");
	treeG->Branch("GenDarkPho_isMet",&GenDarkPho_isMet,"GenDarkPho_isMet/I");

	treeG->Branch("Jet_JEC_u",&Jet_JEC_u,"Jet_JEC_u/D"); 
	treeG->Branch("Jet_JEC_d",&Jet_JEC_d,"Jet_JEC_d/D"); 
//	treeG->Branch("Jet_Res_u",&Jet_Res_u,"Jet_Res_u/D");
//	treeG->Branch("Jet_Res_d",&Jet_Res_d,"Jet_Res_d/D");
	treeG->Branch("Pho_Scale_u",&Pho_Scale_u,"Pho_Scale_u/D"); 
	treeG->Branch("Pho_Scale_d",&Pho_Scale_d,"Pho_Scale_d/D"); 
	treeG->Branch("Pho_IDScale_u",&Pho_IDScale_u,"Pho_IDScale_u/D"); 
	treeG->Branch("Pho_IDScale_d",&Pho_IDScale_d,"Pho_IDScale_d/D"); 
	treeG->Branch("Sig_kfac_u",&Sig_kfac_u,"Sig_kfac_u/D"); 
	treeG->Branch("Sig_kfac_d",&Sig_kfac_d,"Sig_kfac_d/D"); 
//	treeG->Branch("Pho_Sys_u",&Pho_Sys_u,"Pho_Sys_u/D"); 
//	treeG->Branch("Pho_Sys_d",&Pho_Sys_d,"Pho_Sys_d/D"); 
//	treeG->Branch("Pho_Gan_u",&Pho_Gan_u,"Pho_Gan_u/D"); 
//	treeG->Branch("Pho_Gan_d",&Pho_Gan_d,"Pho_gan_d/D"); 
//	treeG->Branch("Pho_Res_u",&Pho_Res_u,"Pho_Res_u/D"); 
//	treeG->Branch("Pho_Res_d",&Pho_Res_d,"Pho_Res_d/D"); 

	treeT = new TTree("Trigger","Trigger");
	treeT->Branch("NVerts",&NVerts,"NVerts/I");
	treeT->Branch("event_weight",&event_weight,"event_weight/D");
	treeT->Branch("GenDarkPho_pt",&GenDarkPho_pt,"GenDarkPho_pt/D");
	treeT->Branch("GenPhoton_pt",&GenPhoton_pt,"GenPhoton_pt/D");
	treeT->Branch("GenDarkPho_eta",&GenDarkPho_eta,"GenDarkPho_eta/D");
	treeT->Branch("GenPhoton_eta",&GenPhoton_eta,"GenPhoton_eta/D");
	treeT->Branch("isHLT165",&isHLT165,"isHLT165/I");
	treeT->Branch("isPrescaled",&isPrescaled,"isPrescaled/I");
	treeT->Branch("TrigPhoton_pt",&TrigPhoton_pt,"TrigPhoton_pt/D");
	treeT->Branch("VertR",&VertR,"VertR/D");
	treeT->Branch("VertZ",&VertZ,"VertZ/D");
}

//Find Medium Id photon with the kinematic requirements
bool phoJetAnalyzer::MediumPhotonIdDecision(int &pho_index){
	bool mediumPhotonID=false;
	for(int i=0; i< phoEt->size()  ;i++) {
		if(IsMediumPhoton(i,phoEt->at(i)) && phoEt->at(i)>=175 && fabs(phoEta->at(i))<1.4442){
			mediumPhotonID = true;
			pho_index = i;
			break;
		}
	}
	return mediumPhotonID;
}

//Loose photon 2016 ID
bool phoJetAnalyzer::IsLoosePhoton(int &in,double pt){
	bool outVal = false;
	bool HoE = false;
	bool sieie = false;
	bool charged = false;
	bool neutral = false;
	bool photon = false;
	bool pixel = false;
	if(phoEta->at(in)<1.4442){
		HoE = (phoHoverE->at(in)  < 0.0597);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01031);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 1.295 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 10.910+0.0148*pt+0.000017*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 3.630+0.0047*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	else{
		HoE = (phoHoverE->at(in)  < 0.0481);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03013);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 1.011 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 5.931+0.0163*pt+0.000014*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 6.641+0.0034*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	if(HoE && sieie && charged && neutral && photon && pixel) outVal = true;
	return outVal;
}

//Medium photon 2016 ID
bool phoJetAnalyzer::IsMediumPhoton(int &in,double pt){
	bool outVal = false;
	bool HoE = false;
	bool sieie = false;
	bool charged = false;
	bool neutral = false;
	bool photon = false;
	bool pixel = false;
	if(phoEta->at(in)<1.4442){
		HoE = (phoHoverE->at(in)  < 0.0396);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.01022);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 0.441 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 2.725+0.0148*pt+0.000017*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 2.571+0.0047*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	else{
		HoE = (phoHoverE->at(in)  < 0.0209);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03001);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 0.442 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 1.715+0.0163*pt+0.000014*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 3.863+0.0034*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	if(HoE && sieie && charged && neutral && photon && pixel) outVal = true;
	return outVal;
}

//Tight photon 216 ID
bool phoJetAnalyzer::IsTightPhoton(int &in,double pt){
	bool outVal = false;
	bool HoE = false;
	bool sieie = false;
	bool charged = false;
	bool neutral = false;
	bool photon = false;
	bool pixel = false;
	if(phoEta->at(in)<1.4442){
		HoE = (phoHoverE->at(in)  < 0.0269);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.00994);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 0.202 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 0.264+0.0148*pt+0.000017*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 2.362+0.0047*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	else{
		HoE = (phoHoverE->at(in)  < 0.0213);
		sieie = (phoSigmaIEtaIEtaFull5x5->at(in)  < 0.03000);
		charged =   (TMath::Max(((phoPFChIso->at(in))-rho*EAElectroncharged(phoEta->at(in))),0.0) < 0.034 );
		neutral =   (TMath::Max(((phoPFNeuIso->at(in))-rho*EAElectronneutral(phoEta->at(in))),0.0) < 0.586+0.0163*pt+0.000014*(pt*pt));
		photon =   (TMath::Max(((phoPFPhoIso->at(in))       - rho*EAElectronphoton(phoEta->at(in))),0.0)  < 2.617+0.0034*pt );
		pixel =   (phohasPixelSeed->at(in)  == 0);
	}
	if(HoE && sieie && charged && neutral && photon && pixel) outVal = true;
	return outVal;
}

//Find jet with correct kinematics
int phoJetAnalyzer::JetDecision(int &pho_index, int &jet_index, double &jetPtSmear){
	int nJets = 0;
	double tempPt=0;
	for(int j=0;j<jetPt->size();j++){
		newJet.SetPxPyPzE(jetPx->at(j),
				  jetPy->at(j),
				  jetPz->at(j),
				  jetEn->at(j));
		tempPt = newJet.Pt();
		
		if(isMC) tempPt = (newJet*jetP4Smear->at(j)).Pt();
		if(tempPt >30.0 && fabs(jetEta->at(j))<2.4) nJets++;
		if(OverlapWithMuon(jetEta->at(j),jetPhi->at(j))) continue;
		if(OverlapWithElectron(jetEta->at(j),jetPhi->at(j))) continue;
		if(DeltaR(jetEta->at(j),phoEta->at(pho_index),jetPhi->at(j),phoPhi->at(pho_index)) < 0.5) continue;
//		if(DeltaPhi(jetPhi->at(j),phoPhi->at(pho_index))<2.0) continue;
		if(tempPt >=170.0 && fabs(jetEta->at(j))<2.4 && jet_index<0){
//			if((jetPt->at(j))>=170.0 && jet_index<0){
			jet_index = j;
			jetPtSmear = tempPt;
//			}
		}
	}
	return nJets;
}

//Find the type of PF object that the Dark Photon is reconstructed
void phoJetAnalyzer::DarkPhotonType(double darketa, double darkphi, int &isPho, int &isJet, int &isMet){
	double dR;
	for(int p=0;p<phoEt->size();p++){
		dR =DeltaR(darketa,phoEta->at(p),darkphi,phoPhi->at(p));
		if(dR<0.5){
			isPho=1;
			break;
		}
	}
	if(isPho<1){
	for(int j=0;j<jetPt->size();j++){
		dR = DeltaR(darketa,jetEta->at(j),darkphi,jetPhi->at(j));
		if(dR<0.1){
			isJet=1;
	  		break;
	  	}
	}
	if(isJet<1) isMet=1;
	}
}

//Check if an electron is nearby
bool phoJetAnalyzer::OverlapWithElectron(double eta, double phi){
	bool overlap = false;
	for(int i=0;i<elePt->size();++i){
		if(eleIDbit->at(i)>>1 == 1 && elePt->at(i)>20.){
			float dRJetEle = DeltaR(eleEta->at(i),eta,elePhi->at(i),phi);
			if(dRJetEle<0.1){
				overlap = true;
				break;
			}
		}
	}
	return overlap;
}

//Check if there is a muon nearby
bool phoJetAnalyzer::OverlapWithMuon(double eta, double phi){
	bool overlap = false;
	for(int k=0;k<muPt->size();++k){
		if(muIDbit->at(k)>>1 == 1 && muPt->at(k) >20.){
			float dRJetMu = DeltaR(muEta->at(k),eta,muPhi->at(k),phi);
			if(dRJetMu<0.1) {
				overlap = true;
				break;
			}
		}
	}
	return overlap;
}

//Separation calculator
double phoJetAnalyzer::DeltaR(double eta1, double eta2, double phi1, double phi2){
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

//Delta phi separation
double phoJetAnalyzer::DeltaPhi(double phi1, double phi2){
	double result = -999.;
	result=(phi1-phi2);
	if(result > M_PI) result -= 2*M_PI;
	if(result <= -M_PI) result += 2*M_PI;
	return result;
}

//Convert intiger to boolean
bool phoJetAnalyzer::to_bool(int s){
        bool rets = (s == 1);
        return rets;
}

//Find a nearby photon
void phoJetAnalyzer::JetIsPhoton(int jet_index,int &jetphoIndex){
        jetphoIndex = -1;
        for(int g=0;g<phoEt->size();g++){
                if(phoEt->at(g)>140
                   && DeltaR(phoEta->at(g),jetEta->at(jet_index),phoPhi->at(g),jetPhi->at(jet_index))<0.05){
                        jetphoIndex = g;
                }
        }
}

//Calculate the number of constituents of a certain type above a pt threshold
int phoJetAnalyzer::NConstituents(int jet_index, int pdgid, double ptcut){
	int count = 0;
	for(int i=0;i<jetConstPdgId->at(jet_index).size();i++){
		if((abs(jetConstPdgId->at(jet_index)[i]) == pdgid || pdgid==-1) && jetConstEt->at(jet_index)[i]>=ptcut) count++;
	}
	return count;
}

//For QCD overlap removal
int phoJetAnalyzer::isPhoJetOverlap(int pho){
        bool retval = 0;
        int genPho=-1;
        for(int g=0;g<mcPID->size();g++){
                if(mcPID->at(g) !=22
                   || DeltaR(mcEta->at(g),phoEta->at(pho),mcPhi->at(g),phoPhi->at(pho)) >= 0.1
                   || fabs(mcPt->at(g)-phoEt->at(pho))/mcPt->at(g) >= 0.1
                   || mcStatus->at(g) != 1) continue;
                genPho=g;
                break;
        }
        if(genPho<0) return 0;
        for(int p=0;p<mcPID->size();p++){
                if((fabs(mcPID->at(p))<9 || fabs(mcPID->at(p)) == 21)
                   && DeltaR(mcEta->at(p),mcEta->at(genPho),mcPhi->at(p),mcPhi->at(genPho)) > 0.05) retval = 1;
        }
        return retval;
}

//For Luminosity calculator
void phoJetAnalyzer::LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr){
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

double phoJetAnalyzer::puweight(Float_t npv){
	Int_t bin = weights_->GetXaxis()->FindBin( npv );
	return weights_->GetBinContent( bin );
}

//Effective Areas
double phoJetAnalyzer::EAPFWorstElectroncharged(double eta){
	float EffectiveArea=0.;
	if (fabs(eta) < 1.0 )   EffectiveArea = 0.075;
	if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.062;
	return EffectiveArea;
}

double phoJetAnalyzer::EAElectroncharged(double eta){
	float EffectiveArea=0.;
	if (fabs(eta) < 1.0 )   EffectiveArea = 0.0360;
	if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0377;
	return EffectiveArea;
}

double phoJetAnalyzer::EAElectronneutral(double eta){
	float EffectiveArea=0.;
	if (fabs(eta) < 1.0 )   EffectiveArea = 0.0597;
	if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0807;
	return EffectiveArea;
}

double phoJetAnalyzer::EAElectronphoton(double eta){
	float EffectiveArea=0.;
	if (fabs(eta) < 1.0 )   EffectiveArea = 0.1210;
	if (fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1107;
	return EffectiveArea;
}

//For systematics calculation
double phoJetAnalyzer::passAnalysisCuts(int jetsystem, int phosystem, int idadj, int kfactadj){
	double rets = 0.0;
	bool isOverlap=false;
	double kfactorScale = 1.0;
	for(int i=0; i< phoEt->size()  ;i++) {
		double phoScale [] = {1.0
				      ,TMath::Sqrt(phoScale_stat_up->at(i)*phoScale_stat_up->at(i)
				      			+phoScale_syst_up->at(i)*phoScale_syst_up->at(i)
							+phoScale_gain_up->at(i)*phoScale_gain_up->at(i))
				      ,TMath::Sqrt(phoScale_stat_dn->at(i)*phoScale_stat_dn->at(i)
				      			+phoScale_syst_dn->at(i)*phoScale_syst_dn->at(i)
							+phoScale_gain_dn->at(i)*phoScale_gain_dn->at(i))};
/*				      ,phoScale_stat_dn->at(i)
				      ,phoScale_syst_up->at(i)
				      ,phoScale_syst_dn->at(i)
				      ,phoScale_gain_up->at(i)
				      ,phoScale_gain_dn->at(i)
				      ,phoResol_rho_up->at(i)
				      ,phoResol_rho_dn->at(i)};
*/
		if(IsMediumPhoton(i,phoEt->at(i)*phoScale[phosystem])
		   && (phoEt->at(i)*phoScale[phosystem])>=175.0
		   && fabs(phoEta->at(i))<1.4442){
			for(int j=0;j<jetPt->size();j++){
				double jetScale [] = {1.0
						      ,1+jetJECUnc->at(j)
						      ,1-jetJECUnc->at(j)};
//						      ,jetP4SmearUp->at(j)
//						      ,jetP4SmearDo->at(j)};
				if(!OverlapWithMuon(jetEta->at(j),jetPhi->at(j)) 
				   && !OverlapWithElectron(jetEta->at(j),jetPhi->at(j)) 
				   && DeltaR(jetEta->at(j),phoEta->at(i),jetPhi->at(j),phoPhi->at(i)) > 0.5
				   && (jetPt->at(j)*jetScale[jetsystem])>=170.0 && fabs(jetEta->at(j))<2.4){
					if(jetNConstituents->at(j)<5){
						if(pileFile.Contains("QCD")) isOverlap = isPhoJetOverlap(i);
						if(pileFile.Contains("df")) kfactorScale = GetKFactor(GenDarkPho_pt,kfactadj);
						rets = GetIDScaleFactor(i,idadj)*event_weight*(!isOverlap)*kfactorScale;
						if(rets>0) return rets;
					}
				}
			}
		}
	}
	return 0.0;
}

//For calculating systematics
void phoJetAnalyzer::SetSystematics(){
	Jet_JEC_u += passAnalysisCuts(1,0,0,0);
	Jet_JEC_d += passAnalysisCuts(2,0,0,0);
	Pho_Scale_u += passAnalysisCuts(0,1,0,0);
	Pho_Scale_d += passAnalysisCuts(0,2,0,0);
	Pho_IDScale_u += passAnalysisCuts(0,0,1,0);
	Pho_IDScale_d += passAnalysisCuts(0,0,-1,0);
//	Sig_kfac_u += passAnalysisCuts(0,0,0,1);
//	Sig_kfac_d += passAnalysisCuts(0,0,0,-1);
}

//Calculate scale factors
double phoJetAnalyzer::GetIDScaleFactor(int pho_index,double adj){
	int x_bin = scaleHist->GetXaxis()->FindBin(phoEta->at(pho_index));
	double Scale = scaleHist->GetBinContent(x_bin,5) + adj*scaleHist->GetBinError(x_bin,5);
	return Scale;
}
	
double phoJetAnalyzer::GetKFactor(double pt,double adj){
	Int_t bin = kfacHist->GetXaxis()->FindBin(pt);
        return kfacHist->GetBinContent(bin) + adj*kfacHist->GetBinError(bin);
}



