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
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"

// constants

const int ELECTRON_PDGID = 11;
const int MUON_PDGID = 13;
const int PHOTON_PDGID = 22;

using namespace trigUtils;

SmartSelectionMonitor
initHistograms(){
  SmartSelectionMonitor mon;

  //counter of events for normalization
  mon.addHistogram(new TH1F("nevent", ";;Number of events", 1, 0, 1) );

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

  // puppi jets
  mon.addHistogram(new TH1F("njetpuppi", ";Number of PUPPI Jets;Events", 20, 0, 20) );
  mon.addHistogram(new TH1F("jetPuppipt", ";PUPPI Jet pT [GeV];Events", 40, 0, 400) ); 
  mon.addHistogram(new TH1F("jetPuppieta", ";PUPPI Jet pseudo-rapidity;Events", 50, 0, 5) );

	// genJets
  mon.addHistogram(new TH1F("ngenjet", ";Number of genJets;Events", 20, 0, 20) );
  mon.addHistogram(new TH1F("genJetpt", ";genJet pT [GeV];Events", 40, 0, 400) ); 
  mon.addHistogram(new TH1F("genJeteta", ";genJet pseudo-rapidity;Events", 50, 0, 5) );

	//genParticles
  mon.addHistogram(new TH1F("ngenParticles", ";Number of genParticles;Events", 200, 0, 200) );

	//leptons
  mon.addHistogram(new TH1F("nleptons", ";Number of Leptons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("leptonpt", ";Lepton pT [GeV];Events", 40, 0, 400) );

	//photons
  mon.addHistogram(new TH1F("nphotons", ";Number of Photons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("photonpt", ";Photon pT [GeV];Events", 40, 0, 400) );

	//neutrinos
  mon.addHistogram(new TH1F("nneutrinos", ";Number of Neutrinos;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("neutrinospt", ";Neutrinos pT [GeV];Events", 40, 0, 400) );

  //bosons
  mon.addHistogram(new TH1F("bosonpt", ";boson pT [GeV];Events", 40, 0, 400) );

  //TProfile between pt of reco and gen jets
  mon.addHistogram(new TProfile("profilePt_pt_gen_vs_reco_above30_percent", ";p_{T}^{genJet};#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, 0, 100, -4, 4 )) ;
  mon.addHistogram(new TProfile("profileEta_pt_gen_vs_reco_above30_percent", ";#eta^{genJet};#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 20, -5, 5, -4, 4 )) ;
  mon.addHistogram(new TProfile("profilePt_pt_jet_vs_pt_Z_above30_percent", ";#eta^{selJet 0};#frac{p_{T}^{selJet 0}}{p_{T}^{Z}}", 20, -5, 5, -4, 4 )) ;
  mon.addHistogram(new TProfile("profilePt_pt_jetPuppi_vs_pt_Z_above30_percent", ";#eta^{selJet 0};#frac{p_{T}^{selJet 0}}{p_{T}^{Z}}", 20, -5, 5, -4, 4 )) ;

  //Correlation plots
  mon.addHistogram(new TH2F("corr_j0_Z_pt", ";p_{T}^{Z};p_{T}^{j0}", 40, 0, 400, 40, 0, 400));

  //Compa entre reco et gen jets
  mon.addHistogram(new TH1F("jetpt", ";Jet pT [GeV];Events", 40, 0, 400) ); 

	//TH1 in eta bins for jet reso
	mon.addHistogram(new TH1F("jetReso_eta0_2p5", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, -1, 1)) ;
	mon.addHistogram(new TH1F("jetReso_eta2p5_2p75", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, -1, 1)) ;
	mon.addHistogram(new TH1F("jetReso_eta2p75_3", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, -1, 1)) ;
	mon.addHistogram(new TH1F("jetReso_eta3_inf", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, -1, 1)) ;

	mon.addHistogram(new TH1F("jetReso_eta3_3p5", ";#frac{p_{T}^{recoJet} - p_{T}^{genJet}}{p_{T}^{genJet}}", 100, -1, 1)) ;
	


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
    gStyle->SetOptStat(111111);
		gROOT->ForceStyle();


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
			mon.fillHisto("nevent", "totalEvent", 0, weight);


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
			pat::JetCollection jets; 
			fwlite::Handle< pat::JetCollection > jetsHandle;
			jetsHandle.getByLabel(ev, "slimmedJets");
			if(jetsHandle.isValid()){ jets = *jetsHandle;}

      pat::JetCollection jetsPuppi;
			fwlite::Handle< pat::JetCollection > jetsPuppiHandle;
      jetsPuppiHandle.getByLabel(ev, "slimmedJetsPuppi");
			if(jetsPuppiHandle.isValid()){ jetsPuppi = *jetsPuppiHandle;}

			reco::GenJetCollection genJets; 
			fwlite::Handle< reco::GenJetCollection > genJetsHandle;
			genJetsHandle.getByLabel(ev, "slimmedGenJets");
			if(genJetsHandle.isValid()){ genJets = *genJetsHandle;}

			reco::GenParticleCollection genParticles;
			fwlite::Handle< reco::GenParticleCollection > genParticlesHandle;
			genParticlesHandle.getByLabel(ev, "prunedGenParticles");
			if(genParticlesHandle.isValid()){ genParticles = *genParticlesHandle;}


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
			mon.fillHisto("ngenjet", tag, genJets.size(), weight);
			mon.fillHisto("njetpuppi", tag, jetsPuppi.size(), weight);
			mon.fillHisto("ngenParticles", tag, genParticles.size(), weight);

			for(unsigned int i = 0; i < jets.size(); i++){
				pat::Jet jet = jets[i];
				TString number = "_" ; number += i;
				mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
				mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
			}

			for(unsigned int i = 0; i < genJets.size(); i++){
				reco::GenJet jet = genJets[i];
				TString number = "_" ; number += i;
				mon.fillHisto("genJetpt", tag+number, jet.pt(), weight);
				mon.fillHisto("genJeteta", tag+number, jet.eta(), weight);
			}

			for(unsigned int i = 0; i < jetsPuppi.size(); i++){
				pat::Jet jet = jetsPuppi[i];
				TString number = "_" ; number += i;
				mon.fillHisto("jetPuppipt", tag+number, jet.pt(), weight);
				mon.fillHisto("jetPuppieta", tag+number, jet.eta(), weight);
			}

		  //Define the boson from the leptons
		  LorentzVector boson(0,0,0,0);
			if (genLeptons.size() < 2) continue;
			if (genLeptons[0].pt() < 20 || genLeptons[1].pt() < 20) continue; //sort of dimuon trigger... sort of...
			boson = genLeptons[0].p4() + genLeptons[1].p4();
			//if ( fabs(boson.mass() - 91) > 15 || boson.pt() < 55) continue; //boosted Z boson
			if ( fabs(boson.mass() - 91) > 15 /*|| boson.pt() < 55*/) continue; //boosted Z boson
			mon.fillHisto("bosonpt", tag, boson.pt(), weight);

			//Clean the different jet collections
			reco::GenJetCollection selGenJets;
			for(unsigned int i = 0; i < genJets.size();i++){
				double pt = genJets[i].pt();
				double eta = genJets[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(genJets[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(genJets[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(genJets[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<30) continue;

				selGenJets.push_back(genJets[i]);

			}
 			std::sort(selGenJets.begin(), selGenJets.end(), utils::sort_CandidatesByPt);


			pat::JetCollection selJets;
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

				selJets.push_back(jets[i]);

			}
 			std::sort(selJets.begin(), selJets.end(), utils::sort_CandidatesByPt);


			pat::JetCollection selJetsPuppi;
			for(unsigned int i = 0; i < jetsPuppi.size();i++){
				double pt = jetsPuppi[i].pt();
				double eta = jetsPuppi[i].eta();

				//cross-clean with leptons and photons and neutrinos
				double minDRlj(9999.), minDRlg(9999.), minDRlnu(9999.);
    		for(size_t ilep=0; ilep<genLeptons.size(); ilep++){
					if(genLeptons[ilep].pt() < 18) continue;	
      		minDRlj = TMath::Min( minDRlj, deltaR(jetsPuppi[i],genLeptons[ilep]) );
    		}
				for(size_t ipho=0; ipho<genPhotons.size(); ipho++){
					if(genPhotons[ipho].pt() < 25) continue;	
      		minDRlg = TMath::Min( minDRlg, deltaR(jetsPuppi[i],genPhotons[ipho]) );
				}
				for(size_t ineu=0; ineu< genNeutrinos.size(); ineu++){
					if(genNeutrinos[ineu].pt() < 18) continue;	
      		minDRlnu = TMath::Min( minDRlnu, deltaR(jetsPuppi[i],genNeutrinos[ineu]) );
				}
    		if(minDRlj<0.4 || minDRlg<0.4 || minDRlnu<0.4) continue;


				if(pt<30) continue;

				selJetsPuppi.push_back(jetsPuppi[i]);

			}
 			std::sort(selJetsPuppi.begin(), selJetsPuppi.end(), utils::sort_CandidatesByPt);


			tag = "selection";
			//jets characteristics after selection	
			mon.fillHisto("njet", tag, selJets.size(), weight);
			mon.fillHisto("ngenjet", tag, selGenJets.size(), weight);
			mon.fillHisto("njetpuppi", tag, selJetsPuppi.size(), weight);

			for(unsigned int i = 0; i < selJets.size(); i++){
				pat::Jet jet = selJets[i];
				TString number = "_" ; number += i;
				mon.fillHisto("jetpt", tag+number, jet.pt(), weight);
				mon.fillHisto("jeteta", tag+number, jet.eta(), weight);
			}

			for(unsigned int i = 0; i < selGenJets.size(); i++){
				reco::GenJet jet = selGenJets[i];
				TString number = "_" ; number += i;
				mon.fillHisto("genJetpt", tag+number, jet.pt(), weight);
				mon.fillHisto("genJeteta", tag+number, jet.eta(), weight);
			}

			for(unsigned int i = 0; i < selJetsPuppi.size(); i++){
				pat::Jet jet = selJetsPuppi[i];
				TString number = "_" ; number += i;
				mon.fillHisto("jetPuppipt", tag+number, jet.pt(), weight);
				mon.fillHisto("jetPuppieta", tag+number, jet.eta(), weight);
			}



			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////// Resolution study ////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////

			tag = "selection_jets";
			for(unsigned int i =0; i <selJets.size(); i++){
				pat::Jet jet = selJets[i];
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
				mon.fillProfile("profilePt_pt_gen_vs_reco_above30_percent", tag, associatedGenJet->pt(), (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				mon.fillProfile("profileEta_pt_gen_vs_reco_above30_percent", tag, associatedGenJet->eta(), (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				if(fabs(associatedGenJet->eta())<2.5) mon.fillHisto("jetReso_eta0_2p5", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				else if(fabs(associatedGenJet->eta())<2.75) mon.fillHisto("jetReso_eta2p5_2p75", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				else if(fabs(associatedGenJet->eta())<3) mon.fillHisto("jetReso_eta2p75_3", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				else mon.fillHisto("jetReso_eta3_inf", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);
				if( associatedGenJet->eta() < 3.5 && associatedGenJet->eta() >3.) mon.fillHisto("jetReso_eta3_3p5", tag, (pt - associatedGenJet->pt())/associatedGenJet->pt(), weight);

			}
			if(selJets.size() >0){
 			 	mon.fillProfile("profilePt_pt_jet_vs_pt_Z_above30_percent", tag, selJets[0].eta(), selJets[0].pt()/boson.pt(), weight);
 			 	mon.fillHisto("corr_j0_Z_pt", tag, selJets[0].pt(), boson.pt(), weight);
			}

			tag = "selection_jetsPuppi";
			for(unsigned int i =0; i <selJetsPuppi.size(); i++){
				pat::Jet jet = selJetsPuppi[i];
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

			}
			if(selJetsPuppi.size() >0){
				mon.fillProfile("profilePt_pt_jetPuppi_vs_pt_Z_above30_percent", tag, selJetsPuppi[0].eta(), selJetsPuppi[0].pt()/boson.pt(), weight);
 			 	mon.fillHisto("corr_j0_Z_pt", tag, selJetsPuppi[0].pt(), boson.pt(), weight);
			}

			//Additional cross-check




		} // end event loop


  	TGraphAsymmErrors *diff_reco_gen_jet_0 = new TGraphAsymmErrors(1); //sera mon graphe d'efficacite
  	TGraphAsymmErrors *diff_reco_gen_jet_1 = new TGraphAsymmErrors(1); //sera mon graphe d'efficacite
  	diff_reco_gen_jet_0->BayesDivide(mon.getHisto("jetpt", "selection_0"), mon.getHisto("genJetpt", "selection_0"));
  	diff_reco_gen_jet_1->BayesDivide(mon.getHisto("jetpt", "selection_1"), mon.getHisto("genJetpt", "selection_1"));
  	diff_reco_gen_jet_0->Write();
  	diff_reco_gen_jet_1->Write();



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







