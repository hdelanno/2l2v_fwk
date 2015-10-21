// -*- C++ -*-
//
// Analysis for the difference in kinematics for 8 and 13 TeV at gen lvl
// 

#include <iostream>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/TrigUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

// constants

const int ELECTRON_PDGID = 11;
const int MUON_PDGID = 13;
const int PHOTON_PDGID = 22;

using namespace trigUtils;

SmartSelectionMonitor
initHistograms(){
  SmartSelectionMonitor mon;

  // jet 
  mon.addHistogram(new TH1F("njet", ";Number of Jets;Events", 20, 0, 20) );
  mon.addHistogram(new TH1F("jetpt", ";Jet pT [GeV];Events", 40, 0, 400) ); 
  mon.addHistogram(new TH1F("jeteta", ";Jet pseudo-rapidity;Events", 50, 0, 5) );
  mon.addHistogram(new TH1F("jetrawen", ";Jet raw energy;Events", 100, 0, 1000) );
  mon.addHistogram(new TH1F("jetnhf", ";Jet neutral hadron energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetnef", ";Jet electromagnetic energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetcef", ";Jet charged electromagnetic energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetchf", ";Jet charged hadron energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetnch", ";Jet charged multiplicity;Events", 100, 0, 100) );
  mon.addHistogram(new TH1F("jetnconst", ";Jet number of constitutes;Events", 100, 0, 100) );
  mon.addHistogram(new TH1F("jetpudsct", ";Jet pileup ID discriminant;Events", 100, -1, 1) );
	//genParticles
  mon.addHistogram(new TH1F("ngenParticles", ";Number of genParticles;Events", 200, 0, 200) );

  //jet no nu
  mon.addHistogram(new TH1F("njet_noNu", ";Number of Jets;Events", 20, 0, 20) );
  mon.addHistogram(new TH1F("jetpt_noNU", ";Jet pT [GeV];Events", 40, 0, 400) );
  mon.addHistogram(new TH1F("jeteta_noNU", ";Jet pseudo-rapidity;Events", 50, 0, 5) );

	//leptons
  mon.addHistogram(new TH1F("nleptons", ";Number of Leptons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("leptonpt", ";Lepton pT [GeV];Events", 40, 0, 400) );

	//photons
  mon.addHistogram(new TH1F("nphotons", ";Number of Photons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("photonpt", ";Photon pT [GeV];Events", 40, 0, 400) );

	//neutrinos
  mon.addHistogram(new TH1F("nneutrinos", ";Number of Neutrinos;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("neutrinospt", ";Neutrinos pT [GeV];Events", 40, 0, 400) );

	//Counter
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("jetcat", ";;Events", 3,0,3) );
  h->GetXaxis()->SetBinLabel(1,"= 0 jet");
  h->GetXaxis()->SetBinLabel(2,"#geq 1 jet");
  h->GetXaxis()->SetBinLabel(3,"VBF");

  TH1F *h2=(TH1F*) mon.addHistogram( new TH1F ("VBFcriteria", ";;Events", 8,0,8) );
  h2->GetXaxis()->SetBinLabel(1,"VBF");
  h2->GetXaxis()->SetBinLabel(2,"non-VBF");
  h2->GetXaxis()->SetBinLabel(3,"njet = 0");
  h2->GetXaxis()->SetBinLabel(4,"njet = 1");
  h2->GetXaxis()->SetBinLabel(5,"|#Delta #eta| < 4.0");
  h2->GetXaxis()->SetBinLabel(6,"M_{jj} < 500 GeV");
  h2->GetXaxis()->SetBinLabel(7,"# central jets > 0");
  h2->GetXaxis()->SetBinLabel(8,"# central bosons != 1");

  TH1F *h3=(TH1F*) mon.addHistogram( new TH1F ("VBFcriteriaLoose", ";;Events", 8,0,8) );
  h3->GetXaxis()->SetBinLabel(1,"VBF");
  h3->GetXaxis()->SetBinLabel(2,"non-VBF");
  h3->GetXaxis()->SetBinLabel(3,"njet = 0");
  h3->GetXaxis()->SetBinLabel(4,"njet = 1");
  h3->GetXaxis()->SetBinLabel(5,"|#Delta #eta| < 3.5");
  h3->GetXaxis()->SetBinLabel(6,"M_{jj} < 200 GeV");
  h3->GetXaxis()->SetBinLabel(7,"# central jets > 0");
  h3->GetXaxis()->SetBinLabel(8,"# central bosons != 1");

  //TProfile between pt of reco and gen jets
  mon.addHistogram(new TProfile("profilePt_pt_gen_vs_reco_above30_percent", ";p_{T}^{genJet};#frac{p_{T}^{genJet} - p_{T}^{recoJet}}{p_{T}^{genJet}}", 100, 0, 100, -1, 1 )) ;
  mon.addHistogram(new TProfile("profileEta_pt_gen_vs_reco_above30_percent", ";#eta^{genJet};#frac{p_{T}^{genJet} - p_{T}^{recoJet}}{p_{T}^{genJet}}", 20, -5, 5, -1, 1 )) ;

  //TH1 in eta bins for jet reso
  mon.addHistogram(new TH1F("jetReso_eta0_2p5", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}};Events", 100, -1, 1)) ;
  mon.addHistogram(new TH1F("jetReso_eta2p5_2p75", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}};Events", 100, -1, 1)) ;
  mon.addHistogram(new TH1F("jetReso_eta2p75_3", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}};Events", 100, -1, 1)) ;
  mon.addHistogram(new TH1F("jetReso_eta3_inf", ";#frac{p_{T}^{recoJet} - p_{T}^{rgenJet}}{p_{T}^{genJet}};Events", 100, -1, 1)) ;

	//Study of dEta
	mon.addHistogram(new TH1F("dEta_Study", ";Pseudo-rapidity span;Events", 20, 0, 10));


	//OBVIOUSLY it has to be at the very end of this function!!!
  return mon; 
}


int main(int argc, char* argv[])
{
  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
	// load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  bool debug = runProcess.getParameter<bool>("debug");
  int maxEvents = runProcess.getParameter<int>("maxevents");
  bool isMC = runProcess.getParameter<bool>("isMC");
  bool photonTriggerStudy = runProcess.getParameter<bool>("triggerstudy");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  edm::ParameterSet pujetidparas = runProcess.getParameter<edm::ParameterSet>("pujetidparas"); 
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
  // TString dtag=runProcess.getParameter<std::string>("dtag");

  // if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; } //FIXME
  // TString outdir=runProcess.getParameter<std::string>("outdir");
  // TString outUrl( outdir );
  // gSystem->Exec("mkdir -p " + outUrl);

  TString outUrl = runProcess.getParameter<std::string>("outfile");

  // std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  // std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");

  // initiating histograms
  SmartSelectionMonitor mon = initHistograms();
  // Loop over files for reliable performance on CRAB
  for(unsigned int iFile=0; iFile<urls.size(); ++iFile){
    // open input file
		
		TFile* inFile = TFile::Open(urls[iFile].c_str());
    if (!inFile) continue; 
  	
		// get ready for the event loop
  	// fwlite::ChainEvent ev(urls);
    fwlite::Event ev(inFile); 
  	size_t totalEntries(0);
  	if (maxEvents == -1) totalEntries = ev.size();
  	else totalEntries = maxEvents;

  	//MC normalization (to 1/pb)
  	double xsecWeight = xsec/totalEntries;

  	if(!isMC) xsecWeight=1.0;
  	if (debug) {
    	printf("DEBUG: xsec= %f\n", xsec);
    	printf("DEBUG: xsecWeight = %f\n", xsecWeight);
  	}


  	// make sure that histogram internally produced in 
  	// lumireweighting are not destroyed when closing the file
  	gROOT->cd();  

  	higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning

  	// ----------------------------------------------------------------------------------------
  	// event loop  
  	// loop on all the events
  	printf("\nTotal entries to process: %zu \n", totalEntries);
  	printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  	printf("Scanning the ntuple :");
  	int treeStep(totalEntries/50);

  	TString tag("all");
  	double weight(1.0); // for testing now. FIXME 
  	//double weight(xsecWeight);

  	for( size_t iev=0; iev<totalEntries; iev++){
    	if(iev%treeStep==0){printf(".");fflush(stdout);}
    	// load the event content from the EDM file
    	ev.to(iev);


    	// fwlite::Handle<edm::TriggerResults> triggerBits;
    	// //    fwlite::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    	// fwlite::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    	// triggerBits.getByLabel(ev,"TriggerResults::HLT");
    	// //  triggerObjects.getByLabel(ev,"selectedPatTrigger");
    	// triggerPrescales.getByLabel(ev,"patTrigger");

    	// const edm::TriggerNames &names = ev.triggerNames(*triggerBits);
    	// std::cout << "\n === TRIGGER PATHS === " << std::endl;
    	// for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    	//   std::cout << "Trigger " << names.triggerName(i) << 
    	// 	", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    	// 	": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
    	// 		<< std::endl;
    	// }
			reco::GenJetCollection jets; //here we analyze only GenJet :)
			fwlite::Handle< reco::GenJetCollection > jetsHandle;
			jetsHandle.getByLabel(ev, "ak4GenJets");
			if(!jetsHandle.isValid()) jetsHandle.getByLabel(ev, "slimmedGenJets"); //small fix for this analysis since we don't have the same collections
			if(jetsHandle.isValid()){ jets = *jetsHandle;}

			reco::GenJetCollection genJets_noNu; //here we analyze only GenJet :)
			fwlite::Handle< reco::GenJetCollection > noNu_Handle;
			noNu_Handle.getByLabel(ev, "ak4GenJetsNoNu");
			//There is no equivalent for miniAOD.
			if(noNu_Handle.isValid()){ genJets_noNu = *noNu_Handle;}

			reco::GenParticleCollection genParticles;
			fwlite::Handle< reco::GenParticleCollection > genParticlesHandle;
			genParticlesHandle.getByLabel(ev, "genParticles");
			if(!genParticlesHandle.isValid()) genParticlesHandle.getByLabel(ev, "prunedGenParticles"); //small fix for this analysis since we don't have the same collections
			if(genParticlesHandle.isValid()){ genParticles = *genParticlesHandle;}

			pat::JetCollection recoJets; 
			fwlite::Handle< pat::JetCollection > recoJetsHandle;
			recoJetsHandle.getByLabel(ev, "slimmedJets");
			if(recoJetsHandle.isValid()){ recoJets = *recoJetsHandle;}



			//create leptons and photons and neutrinos
			reco::GenParticleCollection genLeptons;
			reco::GenParticleCollection genPhotons;
			reco::GenParticleCollection genNeutrinos;
			for (unsigned int i =0; i < genParticles.size(); i++){
				reco::GenParticle genParticle = genParticles[i];
				
				if(genParticle.status() == 2) continue; //drop the shower part of the history
				
				if(genParticle.pt() < 5) continue; //ask a minimal pt to remove garbage
				//let's be even more strict. I look at jet of 20GeV, so let's ask for no particles with a pT below 15
				if(genParticle.pt() <15) continue;

				//leptons
				if(genParticle.status()==1 && (fabs(genParticle.pdgId())==11 || fabs(genParticle.pdgId())==13 || fabs(genParticle.pdgId())==15)){
					genLeptons.push_back(genParticle);
				}
				//photons
				else if(genParticle.status()==1 && genParticle.pdgId()==22){
				  genPhotons.push_back(genParticle);
				}
				//neutrinos
				else if(genParticle.status()==1 && (fabs(genParticle.pdgId())==12 || fabs(genParticle.pdgId())==14 || fabs(genParticle.pdgId())==16)){
					genNeutrinos.push_back(genParticle);
				}


			}

 			std::sort(genLeptons.begin(), genLeptons.end(), utils::sort_CandidatesByPt);
 			std::sort(genPhotons.begin(), genPhotons.end(), utils::sort_CandidatesByPt);
 			std::sort(genNeutrinos.begin(), genNeutrinos.end(), utils::sort_CandidatesByPt);


			//FIXME No clean on photons
			genPhotons.clear();

			tag = "inclusive";

			mon.fillHisto("nleptons", tag, genLeptons.size(), weight);
			for(unsigned int i =0; i < genLeptons.size(); i++){
				TString number = "_"; number +=i;
				mon.fillHisto("leptonpt", tag, genLeptons[i].pt(), weight);
			}
			
			mon.fillHisto("nphotons", tag, genPhotons.size(), weight);
			for(unsigned int i =0; i < genPhotons.size(); i++){
				TString number = "_"; number +=i;
				mon.fillHisto("photonpt", tag, genPhotons[i].pt(), weight);
			}
			
			mon.fillHisto("nneutrinos", tag, genNeutrinos.size(), weight);
			for(unsigned int i =0; i < genNeutrinos.size(); i++){
				TString number = "_"; number +=i;
				mon.fillHisto("neutrinospt", tag, genNeutrinos[i].pt(), weight);
			}




			//jets characteristics			
			mon.fillHisto("njet", tag, jets.size(), weight);
			mon.fillHisto("ngenParticles", tag, genParticles.size(), weight);
			mon.fillHisto("njet_noNu", tag, genJets_noNu.size(), weight);

			for(unsigned int i = 0; i < jets.size(); i++){
				reco::GenJet jet = jets[i];
				TString number = "_" ; number += i;
				mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
				mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
			}

			for (unsigned int i = 0; i < genJets_noNu.size(); i++){
				reco::GenJet jet = genJets_noNu[i];
        TString number = "_" ; number += i;
				mon.fillHisto("jetpt_noNu", tag+number, jet.pt(), weight);
				mon.fillHisto("jeteta_noNu", tag+number, jet.eta(), weight);
			}
			//mon.fillHisto("jetrawen", "all", jets.size(), weight);
			//mon.fillHisto("jetnhf", "all", jets.size(), weight);
			//mon.fillHisto("jetnef", "all", jets.size(), weight);
			//mon.fillHisto("jetcef", "all", jets.size(), weight);
			//mon.fillHisto("jetchf", "all", jets.size(), weight);
			//mon.fillHisto("jetnch", "all", jets.size(), weight);
			//mon.fillHisto("jetnconst", "all", jets.size(), weight);
			//mon.fillHisto("jetpudsct", "all", jets.size(), weight);

		  //Define the bosons from the leptons
		  LorentzVector boson(0,0,0,0);
			if (genLeptons.size() < 2) continue;
			if (genLeptons[0].pt() < 20 || genLeptons[1].pt() < 20) continue; //sort of dimuon trigger... sort of...
			boson = genLeptons[0].p4() + genLeptons[1].p4();
			if ( fabs(boson.mass() - 91) > 15 || boson.pt() < 55) continue; //boosted Z boson






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Firt case : pt of 30GeV for both //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			reco::GenJetCollection selGenJets;
			for(unsigned int i = 0; i < jets.size();i++){
				double pt = jets[i].pt();
				double eta = jets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(jets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(jets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(jets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<30) continue;

				selGenJets.push_back(jets[i]);

			}
 			std::sort(selGenJets.begin(), selGenJets.end(), utils::sort_CandidatesByPt);
			tag = "selection";
      mon.fillHisto("njet", tag, selGenJets.size(), weight);

      for(unsigned int i = 0; i < selGenJets.size(); i++){
        reco::GenJet jet = selGenJets[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF(false);
      if(selGenJets.size()>=2){
				LorentzVector VBFSyst = selGenJets[0].p4() + selGenJets[1].p4();
				double j1eta=selGenJets[0].eta();
				double j2eta=selGenJets[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selGenJets.size(); ijet++){
	  			float jpt=selGenJets[ijet].pt();
	  			float jeta=selGenJets[ijet].eta();
	  			if(jpt<30)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) && (NCentralBoson==1) );
				//isVBF=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) ); //This is like the real one
				//isVBF=( (dEta>3.5) && (VBFSyst.M()>200) && (NCentralJet==0) ); //This is a looser version

				if(dEta<4.0) mon.fillHisto("VBFcriteria", "counter", 4, weight);
				if(VBFSyst.M()<500) mon.fillHisto("VBFcriteria", "counter", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteria", "counter", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteria", "counter", 7, weight);

      }
			else{
				if(selGenJets.size()==0) mon.fillHisto("VBFcriteria", "counter", 2, weight);		
				else mon.fillHisto("VBFcriteria", "counter", 3, weight);		
			}


			if(isVBF){
				mon.fillHisto("VBFcriteria", "counter", 0, weight);		
				mon.fillHisto("jetcat", "counter", 2, weight);

				tag = "vbf-like";
      	mon.fillHisto("njet", tag, selGenJets.size(), weight);

      	for(unsigned int i = 0; i < selGenJets.size(); i++){
        	reco::GenJet jet = selGenJets[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteria", "counter", 1, weight);		
				//other jets category
				if(selGenJets.size()==0) mon.fillHisto("jetcat", "counter", 0, weight);
				else mon.fillHisto("jetcat", "counter", 1, weight);


			}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Second case : pt of 20GeV for both ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			reco::GenJetCollection selGenJets_20;
			for(unsigned int i = 0; i < jets.size();i++){
				double pt = jets[i].pt();
				double eta = jets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(jets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 15) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(jets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(jets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<20) continue;

				selGenJets_20.push_back(jets[i]);

			}
 			std::sort(selGenJets_20.begin(), selGenJets_20.end(), utils::sort_CandidatesByPt);
			tag = "selection20";
      mon.fillHisto("njet", tag, selGenJets_20.size(), weight);

      for(unsigned int i = 0; i < selGenJets_20.size(); i++){
        reco::GenJet jet = selGenJets_20[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF_20(false);
      if(selGenJets_20.size()>=2){
				LorentzVector VBFSyst = selGenJets_20[0].p4() + selGenJets_20[1].p4();
				double j1eta=selGenJets_20[0].eta();
				double j2eta=selGenJets_20[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selGenJets_20.size(); ijet++){
	  			float jpt=selGenJets_20[ijet].pt();
	  			float jeta=selGenJets_20[ijet].eta();
	  			if(jpt<20)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF_20=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) && (NCentralBoson==1) );
				//isVBF_20=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) ); //This is like the real one
				//isVBF_20=( (dEta>3.5) && (VBFSyst.M()>200) && (NCentralJet==0) ); //This is a looser version

				if(dEta<4.0) mon.fillHisto("VBFcriteria", "counter_20", 4, weight);
				if(VBFSyst.M()<500) mon.fillHisto("VBFcriteria", "counter_20", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteria", "counter_20", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteria", "counter_20", 7, weight);

      }
			else{
				if(selGenJets_20.size()==0) mon.fillHisto("VBFcriteria", "counter_20", 2, weight);		
				else mon.fillHisto("VBFcriteria", "counter_20", 3, weight);		
			}


			if(isVBF_20){
				mon.fillHisto("VBFcriteria", "counter_20", 0, weight);		
				mon.fillHisto("jetcat", "counter_20", 2, weight);

				tag = "vbf-like20";
      	mon.fillHisto("njet", tag, selGenJets_20.size(), weight);

      	for(unsigned int i = 0; i < selGenJets_20.size(); i++){
        	reco::GenJet jet = selGenJets_20[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteria", "counter_20", 1, weight);		
				//other jets category
				if(selGenJets_20.size()==0) mon.fillHisto("jetcat", "counter_20", 0, weight);
				else mon.fillHisto("jetcat", "counter_20", 1, weight);


			}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Third case : pt of 20GeV for both + Loose VBF//////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			reco::GenJetCollection selGenJets_20Loose;
			for(unsigned int i = 0; i < jets.size();i++){
				double pt = jets[i].pt();
				double eta = jets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(jets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 15) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(jets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(jets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<20) continue;

				selGenJets_20Loose.push_back(jets[i]);

			}
 			std::sort(selGenJets_20Loose.begin(), selGenJets_20Loose.end(), utils::sort_CandidatesByPt);
			tag = "selection20";
      mon.fillHisto("njet", tag, selGenJets_20Loose.size(), weight);

      for(unsigned int i = 0; i < selGenJets_20Loose.size(); i++){
        reco::GenJet jet = selGenJets_20Loose[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF_20Loose(false);
      if(selGenJets_20Loose.size()>=2){
				LorentzVector VBFSyst = selGenJets_20Loose[0].p4() + selGenJets_20Loose[1].p4();
				double j1eta=selGenJets_20Loose[0].eta();
				double j2eta=selGenJets_20Loose[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selGenJets_20Loose.size(); ijet++){
	  			float jpt=selGenJets_20Loose[ijet].pt();
	  			float jeta=selGenJets_20Loose[ijet].eta();
	  			if(jpt<20)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF_20Loose=( (dEta>3.5) && (VBFSyst.M()>200) && (NCentralJet==0) && (NCentralBoson==1) );
				//isVBF_20Loose=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) ); //This is like the real one
				//isVBF_20Loose=( (dEta>3.5) && (VBFSyst.M()>200) && (NCentralJet==0) ); //This is a looser version

				if(dEta<3.5) mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 4, weight);
				if(VBFSyst.M()<200) mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 7, weight);

      }
			else{
				if(selGenJets_20Loose.size()==0) mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 2, weight);		
				else mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 3, weight);		
			}


			if(isVBF_20Loose){
				mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 0, weight);		
				mon.fillHisto("jetcat", "counter_20Loose", 2, weight);

				tag = "vbf-like20";
      	mon.fillHisto("njet", tag, selGenJets_20Loose.size(), weight);

      	for(unsigned int i = 0; i < selGenJets_20Loose.size(); i++){
        	reco::GenJet jet = selGenJets_20Loose[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteriaLoose", "counter_20Loose", 1, weight);		
				//other jets category
				if(selGenJets_20Loose.size()==0) mon.fillHisto("jetcat", "counter_20Loose", 0, weight);
				else mon.fillHisto("jetcat", "counter_20Loose", 1, weight);


			}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////// Fourth case : RECO LEVEL ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







			pat::JetCollection selRecoJets;
			for(unsigned int i = 0; i < recoJets.size();i++){
				double pt = recoJets[i].pt();
				double eta = recoJets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(recoJets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(recoJets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(recoJets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<30) continue;

				selRecoJets.push_back(recoJets[i]);

			}
 			std::sort(selRecoJets.begin(), selRecoJets.end(), utils::sort_CandidatesByPt);
			tag = "selection_reco";
      mon.fillHisto("njet", tag, selRecoJets.size(), weight);

      for(unsigned int i = 0; i < selRecoJets.size(); i++){
        pat::Jet jet = selRecoJets[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF_reco(false);
      if(selRecoJets.size()>=2){
				LorentzVector VBFSyst = selRecoJets[0].p4() + selRecoJets[1].p4();
				double j1eta=selRecoJets[0].eta();
				double j2eta=selRecoJets[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selRecoJets.size(); ijet++){
	  			float jpt=selRecoJets[ijet].pt();
	  			float jeta=selRecoJets[ijet].eta();
	  			if(jpt<30)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF_reco=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) && (NCentralBoson==1) );
				//isVBF_reco=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) ); //This is like the real one
				//isVBF_reco=( (dEta>3.5) && (VBFSyst.M()>200) && (NCentralJet==0) ); //This is a looser version

				if(dEta<4.0) mon.fillHisto("VBFcriteria", "counter_reco", 4, weight);
				if(VBFSyst.M()<500) mon.fillHisto("VBFcriteria", "counter_reco", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteria", "counter_reco", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteria", "counter_reco", 7, weight);

      }
			else{
				if(selRecoJets.size()==0) mon.fillHisto("VBFcriteria", "counter_reco", 2, weight);		
				else mon.fillHisto("VBFcriteria", "counter_reco", 3, weight);		
			}


			if(isVBF_reco){
				mon.fillHisto("VBFcriteria", "counter_reco", 0, weight);		
				mon.fillHisto("jetcat", "counter_reco", 2, weight);

				tag = "vbf-like_reco";
      	mon.fillHisto("njet", tag, selRecoJets.size(), weight);

      	for(unsigned int i = 0; i < selRecoJets.size(); i++){
        	pat::Jet jet = selRecoJets[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteria", "counter_reco", 1, weight);		
				//other recoJets category
				if(selRecoJets.size()==0) mon.fillHisto("jetcat", "counter_reco", 0, weight);
				else mon.fillHisto("jetcat", "counter_reco", 1, weight);


			}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////// Fifth case : TESTING VBF CUTS (LOIC) ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////        At gen and at reco lvl        ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//Parameters of the selection
			double leadPt = 45;
			double subleadPt = 35;
			double deltaEta = 3.25;
			double VBFMassCut = 600;
			double jetEntreLeadingPt = 30;
			
			reco::GenJetCollection selGenJets_custom;
			for(unsigned int i = 0; i < jets.size();i++){
				double pt = jets[i].pt();
				double eta = jets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(jets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(jets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(jets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<20) continue;

				selGenJets_custom.push_back(jets[i]);

			}

 			std::sort(selGenJets_custom.begin(), selGenJets_custom.end(), utils::sort_CandidatesByPt);
			if( selGenJets_custom.size() > 1){
				if( selGenJets_custom[0].pt() < leadPt || selGenJets_custom[1].pt() < subleadPt) selGenJets_custom.clear();
			}

			tag = "selection_custom";
      mon.fillHisto("njet", tag, selGenJets_custom.size(), weight);

      for(unsigned int i = 0; i < selGenJets_custom.size(); i++){
        reco::GenJet jet = selGenJets_custom[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF_custom(false);
      if(selGenJets_custom.size()>=2){
				LorentzVector VBFSyst = selGenJets_custom[0].p4() + selGenJets_custom[1].p4();
				double j1eta=selGenJets_custom[0].eta();
				double j2eta=selGenJets_custom[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selGenJets_custom.size(); ijet++){
	  			float jpt=selGenJets_custom[ijet].pt();
	  			float jeta=selGenJets_custom[ijet].eta();
	  			if(jpt<jetEntreLeadingPt)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF_custom=( (dEta>deltaEta) && (VBFSyst.M()>VBFMassCut) && (NCentralJet==0) && (NCentralBoson==1) );

				if(dEta<deltaEta) mon.fillHisto("VBFcriteria", "counter_custom", 4, weight);
				if(VBFSyst.M()<VBFMassCut) mon.fillHisto("VBFcriteria", "counter_custom", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteria", "counter_custom", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteria", "counter_custom", 7, weight);

				//Look at dEta
				mon.fillHisto("dEta_Study", tag, dEta, weight);
      }
			else{
				if(selGenJets_custom.size()==0) mon.fillHisto("VBFcriteria", "counter_custom", 2, weight);		
				else mon.fillHisto("VBFcriteria", "counter_custom", 3, weight);		
			}


			if(isVBF_custom){
				mon.fillHisto("VBFcriteria", "counter_custom", 0, weight);		
				mon.fillHisto("jetcat", "counter_custom", 2, weight);

				tag = "vbf-like_custom";
      	mon.fillHisto("njet", tag, selGenJets_custom.size(), weight);

      	for(unsigned int i = 0; i < selGenJets_custom.size(); i++){
        	reco::GenJet jet = selGenJets_custom[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteria", "counter_custom", 1, weight);		
				//other jets category
				if(selGenJets_custom.size()==0) mon.fillHisto("jetcat", "counter_custom", 0, weight);
				else mon.fillHisto("jetcat", "counter_custom", 1, weight);


			}


			////////
			//  RECO
			////////


			pat::JetCollection selRecoJets_custom;
			for(unsigned int i = 0; i < recoJets.size();i++){
				double pt = recoJets[i].pt();
				double eta = recoJets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(recoJets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(recoJets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(recoJets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<20) continue;

				selRecoJets_custom.push_back(recoJets[i]);

			}
 			std::sort(selRecoJets_custom.begin(), selRecoJets_custom.end(), utils::sort_CandidatesByPt);
			if(selRecoJets_custom.size() >1){
				if( selRecoJets_custom[0].pt() < leadPt || selRecoJets_custom[1].pt() < subleadPt) selRecoJets_custom.clear();
			}

			tag = "selection_reco_custom";
      mon.fillHisto("njet", tag, selRecoJets_custom.size(), weight);

      for(unsigned int i = 0; i < selRecoJets_custom.size(); i++){
        pat::Jet jet = selRecoJets_custom[i];
        TString number = "_" ; number += i;
        mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      }



      //VBF tag
      bool isVBF_reco_custom(false);
      if(selRecoJets_custom.size()>=2){
				LorentzVector VBFSyst = selRecoJets_custom[0].p4() + selRecoJets_custom[1].p4();
				double j1eta=selRecoJets_custom[0].eta();
				double j2eta=selRecoJets_custom[1].eta();
				double dEta = fabs(j1eta-j2eta);

				int NCentralJet(0), NCentralBoson(0);
				double MaxEta, MinEta;
				if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
				else            { MinEta=j2eta; MaxEta=j1eta;}
				for(size_t ijet=2; ijet<selRecoJets_custom.size(); ijet++){
	  			float jpt=selRecoJets_custom[ijet].pt();
	  			float jeta=selRecoJets_custom[ijet].eta();
	  			if(jpt<jetEntreLeadingPt)continue; 
	  			if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
				}

				if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
				isVBF_reco_custom=( (dEta>deltaEta) && (VBFSyst.M()>VBFMassCut) && (NCentralJet==0) && (NCentralBoson==1) );

				if(dEta<deltaEta) mon.fillHisto("VBFcriteria", "counter_reco_custom", 4, weight);
				if(VBFSyst.M()<VBFMassCut) mon.fillHisto("VBFcriteria", "counter_reco_custom", 5, weight);
				if(NCentralJet > 0) mon.fillHisto("VBFcriteria", "counter_reco_custom", 6, weight);
				if(NCentralBoson != 1) mon.fillHisto("VBFcriteria", "counter_reco_custom", 7, weight);

				//Look at dEta
				mon.fillHisto("dEta_Study", tag, dEta, weight);

      }
			else{
				if(selRecoJets_custom.size()==0) mon.fillHisto("VBFcriteria", "counter_reco_custom", 2, weight);		
				else mon.fillHisto("VBFcriteria", "counter_reco_custom", 3, weight);		
			}


			if(isVBF_reco_custom){
				mon.fillHisto("VBFcriteria", "counter_reco_custom", 0, weight);		
				mon.fillHisto("jetcat", "counter_reco_custom", 2, weight);

				tag = "vbf-like_reco_custom";
      	mon.fillHisto("njet", tag, selRecoJets_custom.size(), weight);

      	for(unsigned int i = 0; i < selRecoJets_custom.size(); i++){
        	pat::Jet jet = selRecoJets_custom[i];
        	TString number = "_" ; number += i;
        	mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
        	mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
      	}

			}
			else{
				mon.fillHisto("VBFcriteria", "counter_reco_custom", 1, weight);		
				//other recoJets category
				if(selRecoJets_custom.size()==0) mon.fillHisto("jetcat", "counter_reco_custom", 0, weight);
				else mon.fillHisto("jetcat", "counter_reco_custom", 1, weight);


			}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Energy resolution of jets ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      tag = "selection_jets";
      for(unsigned int i =0; i <selRecoJets.size(); i++){
        pat::Jet jet = selRecoJets[i];
        double pt = jet.pt();
        double eta = jet.eta();

        const reco::GenJet* associatedGenJet = jet.genJet();

        if(!associatedGenJet || associatedGenJet->pt() < 20 || deltaR(jet, *associatedGenJet) > 0.25) continue;
        //let's check if this associatedGenJet is compatible with one of my selGenJet
        double minDR(9999.);
        reco::GenJet matchedGenJet;
        for(unsigned int j =0; j < selGenJets.size(); j++){
          if ( deltaR(*associatedGenJet, selGenJets[j]) < minDR){
            minDR = deltaR(*associatedGenJet, selGenJets[j]);
            matchedGenJet = selGenJets[j];
          }
        }
        if( minDR > 0.25 || deltaR(matchedGenJet, jet) > 0.25) continue;

        //NB: Vu ma selection, on est d'office dans le cas ou le pt est superieur a 30 GeV
        mon.fillProfile("profilePt_pt_gen_vs_reco_above30_percent", tag, associatedGenJet->pt(), (associatedGenJet->pt() - pt)/associatedGenJet->pt(), weight);
        mon.fillProfile("profileEta_pt_gen_vs_reco_above30_percent", tag, associatedGenJet->eta(), (associatedGenJet->pt() - pt)/associatedGenJet->pt(), weight);
        if(fabs(associatedGenJet->eta())<2.5) mon.fillHisto("jetReso_eta0_2p5", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
        else if(fabs(associatedGenJet->eta())<2.75) mon.fillHisto("jetReso_eta2p5_2p75", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
        else if(fabs(associatedGenJet->eta())<3) mon.fillHisto("jetReso_eta2p75_3", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
        else mon.fillHisto("jetReso_eta3_inf", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);


      }




		} // end event loop
  } // end files loop 
  printf(" done.\n"); 

  //save control plots to file
  // outUrl += "/";
  // outUrl += outFileUrl + ".root";
  // printf("Results saved in %s\n", outUrl.Data());
  // TFile *ofile=TFile::Open(outUrl, "recreate");

  printf("Results saved in %s\n", outUrl.Data());
  TFile *ofile=TFile::Open(outUrl, "recreate");

  mon.Write();
  ofile->Close();
}  







