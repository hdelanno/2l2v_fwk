#include <TH1.h>
#include <TFile.h>
#include <TLine.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TMarker.h>

#include <iostream>

using namespace std ;

TGraph* MakeROCCurve(TH1F* hSigEff, TH1F* hBkgEff){
  // We need to make a TGraph, so we need two arrays of values: (x,y) points, and the total number of points
  // Each TH1 contains the efficiency in a bin for a cut (so 10 bins = 10 different cuts)
  const Int_t nBins = hSigEff->GetNbinsX()+1 ;
  Float_t sigEff[nBins] ;
  Float_t bkgEff[nBins] ;

  // Fill in the (x,y) points from the histograms
  sigEff[0] = 0.0 ;
  bkgEff[0] = 0.0 ;
  for(Int_t bin=1 ; bin<=nBins ; ++bin){
    sigEff[bin] = hSigEff->GetBinContent(bin) ;
    //bkgEff[bin] = hBkgEff->GetBinContent(bin) ;
    bkgEff[bin] = 1. - hBkgEff->GetBinContent(bin) ;
  }

  // Make the TGraph.  We usually put sigEff on the y axis
  //TGraph* ROC_graph = new TGraph(nBins,bkgEff,sigEff) ;
  TGraph* ROC_graph = new TGraph(nBins,sigEff,bkgEff) ;
  return ROC_graph ;
}

TMarker* MakePoint(TH1F* hSigEff, TH1F* hBkgEff, float sigEffTarget){
  float sigEff = -1 ;
  float bkgEff = -1 ;
  for(Int_t bin=hSigEff->GetNbinsX()+1 ; bin>=1 ; --bin){
  	if(hSigEff->GetBinContent(bin)>sigEffTarget){
      sigEff = hSigEff->GetBinContent(bin) ;
      bkgEff = 1-hBkgEff->GetBinContent(bin) ;
      break ;
    }
  }
  TMarker* marker = new TMarker(sigEff,bkgEff,20) ;
  return marker ;
}

void ROC_and_TH1_eff(){

	gROOT->SetBatch();      //COMMANDE MAGIQUE QUI PERMET DE NE PAS AFFICHER LES CANVAS :D

	// Open the file and get the histograms
  // hSig and hBkg contain information about the signal efficiency of a given cut
  // I took them from my pixel match study at HLT
  TFile* file_in = new TFile("plotter.root", "READ") ;
  TFile* file_out = new TFile("ROC.root", "RECREATE") ;
	// Create a TH1 eff by bin from the MVA
  TH1F* hSig[4][3];
  TH1F* hBkg[4][3];
  TGraph* ROC_curve[4][3];
  TCanvas* canvas[4][3];
  TCanvas* cMVA[4][3];

	TString name;
	TString name_tmp;
  for (int eta = 0; eta < 4; eta++){
    //Loop over categories in eta (0-2.5 ; 2.5-2.75 ; 2.75-3 ; 3-4.7 )
  	name = "MVA_eta";
    name += eta;
    name += "_pt";
    name_tmp = name;
    for (int pt = 0; pt < 3; pt++){
      //Loop over categories in pt ( 15-20 ; 20-30 ; 30-inf )
      name = name_tmp;
      name += pt; 
  		file_in->cd();
  		TH1F* hSig_tmp = (TH1F*) file_in->Get("Z#rightarrow ll/truejetsid_"+name) ;
  		TH1F* hBkg_tmp = (TH1F*) file_in->Get("Z#rightarrow ll/pujetsid_"+name) ;
			file_out->cd();

  		cMVA[eta][pt] = new TCanvas("cMVA"+name, name, 100, 100, 800, 600) ;
  		cMVA[eta][pt]->cd();
  		hSig_tmp->Draw();
  		hSig_tmp->SetLineColor(kRed);
  		hSig_tmp->SetMarkerColor(kRed);
  		hBkg_tmp->Draw("same");
  		hBkg_tmp->SetLineColor(kGreen);

  		hSig[eta][pt] = new TH1F("hSig"+name, "MVA efficiency for true jets ("+name+");Position of the cut;Efficiency", 200, 0, 200);
      hBkg[eta][pt] = new TH1F("hBkg"+name, "MVA efficiency for PU jets ("+name+");Position of the cut;Efficiency", 200, 0, 200);

			for (int i = 0 ; i <= 200 ; i++){
				hSig[eta][pt]->SetBinContent(i, hSig_tmp->Integral(i, 201)/hSig_tmp->Integral()); 
				hBkg[eta][pt]->SetBinContent(i, hBkg_tmp->Integral(i, 201)/hBkg_tmp->Integral()); 
			}

  		// Create a ROC curve from the histogram
  		ROC_curve[eta][pt] = MakeROCCurve(hSig[eta][pt], hBkg[eta][pt]) ;

  		// We need to draw the TGraph before we can change its axes, style etc
  		canvas[eta][pt] = new TCanvas("canvas"+name, name, 100, 100, 800, 600) ;
  		ROC_curve[eta][pt]->Draw("AP") ;
  		ROC_curve[eta][pt]->GetXaxis()->SetTitle("#epsilon_{signal}") ;
  		ROC_curve[eta][pt]->GetYaxis()->SetTitle("1 - #epsilon_{background}"    ) ;
  		ROC_curve[eta][pt]->SetMarkerStyle(20) ;
  		ROC_curve[eta][pt]->SetMarkerColor(kRed) ;
  		ROC_curve[eta][pt]->SetLineColor(kRed) ;
  		ROC_curve[eta][pt]->SetTitle("ROC curve for "+name) ;

  		// Make some markers
  		const Int_t nMarkers = 3 ;
  		float effTargets[nMarkers] = { 0.99 , 0.95 , 0.90 } ;
  		TMarker* markers[nMarkers] ;
  		TLatex*  markerLabels[nMarkers] ;
  		TLine* lineX[nMarkers] ;
  		TLine* lineY[nMarkers] ;
  		for(Int_t i=0 ; i<nMarkers ; ++i){
    		float eff = effTargets[i] ;
    		markers[i] = MakePoint(hSig[eta][pt], hBkg[eta][pt], eff) ;
    		markers[i]->SetMarkerStyle(34) ;
    		markers[i]->SetMarkerColor(kBlue) ;
    		float x = markers[i]->GetX() ;
    		float y = markers[i]->GetY() ;
    		markerLabels[i] = new TLatex(x, y-0.05, Form("#varepsilon = %.0f%%", 100*eff)) ;
    		markerLabels[i]->SetTextSize(0.04) ;
    		markerLabels[i]->SetTextAngle(-45) ;

    		lineX[i] = new TLine(0,y,x,y) ;
    		lineY[i] = new TLine(x,0,x,y) ;
    		lineX[i]->SetLineStyle(kDashed) ;
    		lineY[i]->SetLineStyle(kDashed) ;
    		lineX[i]->SetLineWidth(2) ;
    		lineY[i]->SetLineWidth(2) ;
    		lineX[i]->SetLineColor(kBlue) ;
    		lineY[i]->SetLineColor(kBlue) ;
  		}

  		// Save it to canvas and to file
  		ROC_curve[eta][pt]->Draw("ACP") ;
  		for(Int_t i=0 ; i<nMarkers ; ++i){
    		markers[i]->Draw() ;
    		markerLabels[i]->Draw() ;
    		lineX[i]->Draw() ;
    		lineY[i]->Draw() ;
  		}

  		// Add the grid
  		canvas[eta][pt]->SetGridx() ;
  		canvas[eta][pt]->SetGridy() ;

  		//canvas[eta][pt]->Print("ROC/ROC"+name+".eps") ;
  		//canvas[eta][pt]->Print("RPC/ROC"+name+".png") ;

  		ROC_curve[eta][pt]->Write() ;
      canvas[eta][pt]->Write();
      cMVA[eta][pt]->Write();

		}
  }
  file_out->Close() ;
	exit(0);
}



