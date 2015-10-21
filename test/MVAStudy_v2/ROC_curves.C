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
    bkgEff[bin] = hBkgEff->GetBinContent(bin) ;
  }

  // Make the TGraph.  We usually put sigEff on the y axis
  TGraph* ROC_graph = new TGraph(nBins,bkgEff,sigEff) ;
  return ROC_graph ;
}

TMarker* MakePoint(TH1F* hSigEff, TH1F* hBkgEff, float sigEffTarget){
  float sigEff = -1 ;
  float bkgEff = -1 ;
  for(Int_t bin=1 ; bin<=hSigEff->GetNbinsX()+1 ; ++bin){
    if(hSigEff->GetBinContent(bin)>sigEffTarget){
      sigEff = hSigEff->GetBinContent(bin) ;
      bkgEff = hBkgEff->GetBinContent(bin) ;
      break ;
    }
  }
  TMarker* marker = new TMarker(bkgEff,sigEff,20) ;
  return marker ;
}

void ROC_curves(){
  // Open the file and get the histograms
  // hSig and hBkg contain information about the signal efficiency of a given cut
  // I took them from my pixel match study at HLT
  TFile* file_in = new TFile("eff_13xb25_trigger_17_8_phiB.root", "READ") ;
  TH1F* hSig = (TH1F*) file_in->Get("h_sig") ;
  TH1F* hBkg = (TH1F*) file_in->Get("h_bkg") ;

  // Create a ROC curve from the histogram
  TGraph* ROC_curve = MakeROCCurve(hSig, hBkg) ;

  // We need to draw the TGraph before we can change its axes, style etc
  TCanvas* canvas = new TCanvas("canvas", "", 100, 100, 800, 600) ;
  ROC_curve->Draw("AP") ;
  ROC_curve->GetXaxis()->SetTitle("background efficiency") ;
  ROC_curve->GetYaxis()->SetTitle("signal efficiency"    ) ;
  ROC_curve->SetMarkerStyle(20) ;
  ROC_curve->SetMarkerColor(kRed) ;
  ROC_curve->SetLineColor(kRed) ;
  ROC_curve->SetTitle("ROC curve for #Delta#phi_{1}") ;

  // Make some markers
  const Int_t nMarkers = 3 ;
  float effTargets[nMarkers] = { 0.99 , 0.95 , 0.90 } ;
  TMarker* markers[nMarkers] ;
  TLatex*  markerLabels[nMarkers] ;
  TLine* lineX[nMarkers] ;
  TLine* lineY[nMarkers] ;
  for(Int_t i=0 ; i<nMarkers ; ++i){
    float eff = effTargets[i] ;
    markers[i] = MakePoint(hSig, hBkg, eff) ;
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
  ROC_curve->Draw("ACP") ;
  for(Int_t i=0 ; i<nMarkers ; ++i){
    markers[i]->Draw() ;
    markerLabels[i]->Draw() ;
    lineX[i]->Draw() ;
    lineY[i]->Draw() ;
  }

  // Add the grid
  canvas->SetGridx() ;
  canvas->SetGridy() ;

  canvas->Print("ROC.eps") ;
  canvas->Print("ROC.png") ;

  TFile* file_out = new TFile("ROC.root", "RECREATE") ;
  ROC_curve->Write() ;
  file_out->Close() ;
}



