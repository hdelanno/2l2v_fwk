//To execute : .X MakeBeauty.C

#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMath.h>
#include <math.h>
#include <vector>
#include <TFile.h>

void MakeBeauty(){

  gROOT->SetBatch();      //COMMANDE MAGIQUE QUI PERMET DE NE PAS AFFICHER LES CANVAS :D


	//Open the files and recreate the output file
  TString path = "/localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/CompaCinematique/";

  TFile *f_plotter = TFile::Open(path+"plotter.root");

	TFile *f_beauty = new TFile(path+"beauty.root","RECREATE");


	//Def of the histo of interest
	/// njet ///
	TH1F *njet_8;
	TH1F *njet_13;

	/// pt ///
	TH1F *pt0_8, *pt1_8;
	TH1F *pt0_13, *pt1_13;

	/// eta ///
	TH1F *eta0_8, *eta1_8;
	TH1F *eta0_13, *eta1_13;

	//Legend
//	TLegend *leg;

	/// My canvas for the output ///
	TCanvas *c_njet = new TCanvas("c_njet");
	TCanvas *c_pt0 = new TCanvas("c_pt1");
	TCanvas *c_pt1 = new TCanvas("c_pt2");
	TCanvas *c_eta0 = new TCanvas("c_eta0");
	TCanvas *c_eta1 = new TCanvas("c_eta1");

	//Take the histos of interest
	f_plotter->cd("qqH(600)#rightarrow ZZ (8TeV)");
	njet_8 = (TH1F*) gDirectory->Get("vbf-like_njet");
	pt0_8 = (TH1F*) gDirectory->Get("vbf-like_0_jetpt");
	pt1_8 = (TH1F*) gDirectory->Get("vbf-like_1_jetpt");
	eta0_8 = (TH1F*) gDirectory->Get("vbf-like_0_jeteta");
	eta1_8 = (TH1F*) gDirectory->Get("vbf-like_1_jeteta");

	f_plotter->cd("qqH(600)#rightarrow ZZ (13TeV)");
	njet_13 = (TH1F*) gDirectory->Get("vbf-like_njet");
	pt0_13 = (TH1F*) gDirectory->Get("vbf-like_0_jetpt");
	pt1_13 = (TH1F*) gDirectory->Get("vbf-like_1_jetpt");
	eta0_13 = (TH1F*) gDirectory->Get("vbf-like_0_jeteta");
	eta1_13 = (TH1F*) gDirectory->Get("vbf-like_1_jeteta");

	/*
	//Rebinning (have to be done first)
  ggH400_pt->Rebin(2);
	ggH600_pt->Rebin(2);
	VBF400_pt->Rebin(2);
	VBF600_pt->Rebin(2);
	ggH400_pt_sub->Rebin(2);
	ggH600_pt_sub->Rebin(2);
	VBF400_pt_sub->Rebin(2);
	VBF600_pt_sub->Rebin(2);
  ggH400_eta->Rebin(4);
	ggH600_eta->Rebin(4);
	VBF400_eta->Rebin(4);
	VBF600_eta->Rebin(4);
	ggH400_eta_sub->Rebin(4);
	ggH600_eta_sub->Rebin(4);
	VBF400_eta_sub->Rebin(4);
	VBF600_eta_sub->Rebin(4);
	*/	
	/*
	//Choose the right style
	/// ggH400 -> black ///
	ggH400_pt->SetLineColor(kBlack);
	ggH400_pt_sub->SetLineColor(kBlack);
	ggH400_eta->SetLineColor(kBlack);
	ggH400_eta_sub->SetLineColor(kBlack);
	/// ggH600 -> red ///
	ggH600_pt->SetLineColor(kRed);
	ggH600_pt_sub->SetLineColor(kRed);
	ggH600_eta->SetLineColor(kRed);
	ggH600_eta_sub->SetLineColor(kRed);
	/// VBF400 -> green ///
  VBF400_pt->SetLineColor(kGreen);
  VBF400_pt_sub->SetLineColor(kGreen);
  VBF400_eta->SetLineColor(kGreen);
  VBF400_eta_sub->SetLineColor(kGreen);
	/// VBF600 -> blue ///
  VBF600_pt->SetLineColor(kBlue);
  VBF600_pt_sub->SetLineColor(kBlue);
  VBF600_eta->SetLineColor(kBlue);
  VBF600_eta_sub->SetLineColor(kBlue);
	*/
/*	
	/// first plot has to : ///
	//start at 0
	njet_8->SetMinimum(0);
	pt0_8->SetMinimum(0);
	pt1_8->SetMinimum(0);
	eta0_8->SetMinimum(0);
	eta1_8->SetMinimum(0);

	//No stats
  njet_8->SetStats(kFALSE);
	pt0_8->SetStats(kFALSE);
	pt1_8->SetStats(kFALSE);
	eta0_8->SetStats(kFALSE);
	eta1_8->SetStats(kFALSE);
	*/
	/*
	//No title
	njet_8->SetTitle();
	pt0_8->SetTitle();
	pt1_8->SetTitle();
	eta0_8->SetTitle();
	eta1_8->SetTitle();
	*/
/*
	//Have the good range in X :
  //VBF400_pt->GetXaxis()->SetRangeUser(0, 400);
  VBF400_pt->SetAxisRange(0, 400);
  ggH400_pt->SetAxisRange(0, 400);
  VBF400_pt_sub->SetAxisRange(0, 250);
  ggH400_pt_sub->SetAxisRange(0, 250);
  VBF400_eta->SetAxisRange(-6, 6);
  ggH400_eta->SetAxisRange(-6, 6);
  VBF400_eta_sub->SetAxisRange(-6, 6);
  ggH400_eta_sub->SetAxisRange(-6, 6);

*/
	//Make the legends
	leg = new TLegend(0.6,0.7,0.89,0.89);
	leg->AddEntry(pt0_8, "8TeV : qqH(600)#rightarrow ZZ #rightarrow 2l2#nu", "l");
	leg->AddEntry(pt0_13, "13TeV : qqH(600)#rightarrow ZZ #rightarrow 2l2#nu", "l");

	//Let's fill the canvas
	c_njet->cd();
	njet_8->DrawNormalized();
	njet_13->DrawNormalized("same");
	leg->Draw();

	c_pt0->cd();
	pt0_8->DrawNormalized();
	pt0_13->DrawNormalized("same");
	leg->Draw();

	c_pt1->cd();
	pt1_8->DrawNormalized();
	pt1_13->DrawNormalized("same");
	leg->Draw();

	c_eta0->cd();
	eta0_8->DrawNormalized();
	eta0_13->DrawNormalized("same");
	leg->Draw();

	c_eta1->cd();
	eta1_8->DrawNormalized();
	eta1_13->DrawNormalized("same");
	leg->Draw();

 	//Write the canvas
	f_beauty->cd();
	c_njet->Write();
	c_pt0->Write();
	c_pt1->Write();
	c_eta0->Write();
	c_eta1->Write();
	
     c_njet->SaveAs(path+"beauty/njet.pdf");
     c_pt0->SaveAs(path+"beauty/pt0.pdf");
     c_pt1->SaveAs(path+"beauty/pt1.pdf");
     c_eta0->SaveAs(path+"beauty/eta0.pdf");
     c_eta1->SaveAs(path+"beauty/eta1.pdf");
	//exit(0);

}
