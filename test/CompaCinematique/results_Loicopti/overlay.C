//To execute : .X overlay.C

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

void overlay(){

  gROOT->SetBatch();      //COMMANDE MAGIQUE QUI PERMET DE NE PAS AFFICHER LES CANVAS :D


	//Open the files and recreate the output file
  TString path = "/localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/CompaCinematique/results_Loicopti/";


  TFile *f_plotter1 = TFile::Open(path+"set1.root");
  TFile *f_plotter2 = TFile::Open(path+"set2.root");
  TFile *f_plotter3 = TFile::Open(path+"set3.root");
	
	TFile *f_overlay = new TFile(path+"combination.root","RECREATE");


	//Def of the histo of interest
	//dEta
	TH1F *dEta_1, *dEta_2, *dEta_3;
	TH1F *dEta_Reco1, *dEta_Reco2, *dEta_Reco3;

	//jetcat_VBFcriteria
	TH1F *VBFcriteria_1, *VBFcriteria_2, *VBFcriteria_3;
	TH1F *VBFcriteria_Reco1, *VBFcriteria_Reco2, *VBFcriteria_Reco3;

	//Legend
	//TLegend *leg;

	/// My canvas for the output ///
	TCanvas *c_dEta = new TCanvas("c_dEta");
	TCanvas *c_dEta_gen = new TCanvas("c_dEta_gen");
	TCanvas *c_dEta_reco = new TCanvas("c_dEta_reco");

	TCanvas *c_VBFcriteria = new TCanvas("c_VBFcriteria");
	TCanvas *c_VBFcriteria_gen = new TCanvas("c_VBFcriteria_gen");
	TCanvas *c_VBFcriteria_reco = new TCanvas("c_VBFcriteria_reco");
	
	
	//Take the histos of interest
	f_plotter1->cd();
	dEta_1 = (TH1F*) gDirectory->Get("selection_custom_dEta_Study");
  VBFcriteria_1 = (TH1F*) gDirectory->Get("counter_custom_VBFcriteria");
	dEta_Reco1 = (TH1F*) gDirectory->Get("selection_reco_custom_dEta_Study");
  VBFcriteria_Reco1 = (TH1F*) gDirectory->Get("counter_reco_custom_VBFcriteria");

	f_plotter2->cd();
	dEta_2 = (TH1F*) gDirectory->Get("selection_custom_dEta_Study");
  VBFcriteria_2 = (TH1F*) gDirectory->Get("counter_custom_VBFcriteria");
	dEta_Reco2 = (TH1F*) gDirectory->Get("selection_reco_custom_dEta_Study");
  VBFcriteria_Reco2 = (TH1F*) gDirectory->Get("counter_reco_custom_VBFcriteria");

	

	f_plotter3->cd();
	dEta_3 = (TH1F*) gDirectory->Get("selection_custom_dEta_Study");
  VBFcriteria_3 = (TH1F*) gDirectory->Get("counter_custom_VBFcriteria");
	dEta_Reco3 = (TH1F*) gDirectory->Get("selection_reco_custom_dEta_Study");
  VBFcriteria_Reco3 = (TH1F*) gDirectory->Get("counter_reco_custom_VBFcriteria");



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
	
	//Choose the right style
	/// Set1 -> black ///
	dEta_1->SetLineColor(kBlack);
	dEta_Reco1->SetLineColor(kBlack);
	VBFcriteria_1->SetLineColor(kBlack);
	VBFcriteria_Reco1->SetLineColor(kBlack);
	/// Set2 -> red ///
	dEta_2->SetLineColor(kRed);
	dEta_Reco2->SetLineColor(kRed);
	VBFcriteria_2->SetLineColor(kRed);
	VBFcriteria_Reco2->SetLineColor(kRed);
	/// Set3 -> green ///
 	dEta_3->SetLineColor(kGreen);
	dEta_Reco3->SetLineColor(kGreen);
	VBFcriteria_3->SetLineColor(kGreen);
	VBFcriteria_Reco3->SetLineColor(kGreen);

	//Change reco in dash
	dEta_Reco1->SetLineStyle(7);
	dEta_Reco2->SetLineStyle(7);
	dEta_Reco3->SetLineStyle(7);
	VBFcriteria_Reco1->SetLineStyle(7);
	VBFcriteria_Reco2->SetLineStyle(7);
	VBFcriteria_Reco3->SetLineStyle(7);

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
	leg->AddEntry(dEta_1, "Set 1: jet1pT>45, jet2pT>30, dEta>3.25, Mjj>600", "l");
	leg->AddEntry(dEta_2, "Set 2: jet1pT>25, jet2pT>20, dEta>4.00, Mjj>700", "l");
	leg->AddEntry(dEta_3, "Set 3: jet1pT>45, jet2pT>20, dEta>4.25, Mjj>850", "l");

	//Make the legends
	leg_all = new TLegend(0.6,0.7,0.89,0.89);
	leg_all->AddEntry(dEta_1, "Set 1: Gen, jet1pT>45, jet2pT>30, dEta>3.25, Mjj>600", "l");
	leg_all->AddEntry(dEta_2, "Set 2: Gen, jet1pT>25, jet2pT>20, dEta>4.00, Mjj>700", "l");
	leg_all->AddEntry(dEta_3, "Set 3: Gen, jet1pT>45, jet2pT>20, dEta>4.25, Mjj>850", "l");
	leg_all->AddEntry(dEta_Reco1, "Set 1: Rec, jet1pT>45, jet2pT>30, dEta>3.25, Mjj>600", "l");
	leg_all->AddEntry(dEta_Reco2, "Set 2: Rec, jet1pT>25, jet2pT>20, dEta>4.00, Mjj>700", "l");
	leg_all->AddEntry(dEta_Reco3, "Set 3: Rec, jet1pT>45, jet2pT>20, dEta>4.25, Mjj>850", "l");


	//Let's fill the canvas
	
	/// dEta ///
	c_dEta_gen->cd();
	dEta_1->DrawNormalized();
	dEta_2->DrawNormalized("same");
	dEta_3->DrawNormalized("same");
	leg->Draw();

	c_dEta_reco->cd();
	dEta_Reco1->DrawNormalized();
	dEta_Reco2->DrawNormalized("same");
	dEta_Reco3->DrawNormalized("same");
	leg->Draw();

	c_dEta->cd();
	dEta_1->DrawNormalized();
	dEta_2->DrawNormalized("same");
	dEta_3->DrawNormalized("same");
  dEta_Reco1->DrawNormalized("same");
  dEta_Reco2->DrawNormalized("same");
  dEta_Reco3->DrawNormalized("same");
	leg_all->Draw();

	
	/*
	/// dEta ///
	c_dEta_gen->cd();
	dEta_1->Draw();
	dEta_2->Draw("same");
	dEta_3->Draw("same");
	leg->Draw();

	c_dEta_reco->cd();
	dEta_Reco1->Draw();
	dEta_Reco2->Draw("same");
	dEta_Reco3->Draw("same");
	leg->Draw();

	c_dEta->cd();
	dEta_1->Draw();
	dEta_2->Draw("same");
	dEta_3->Draw("same");
  dEta_Reco1->Draw("same");
  dEta_Reco2->Draw("same");
  dEta_Reco3->Draw("same");
	leg_all->Draw();
	*/

	/// VBFcriteria ///
	c_VBFcriteria_gen->cd();
	VBFcriteria_1->Draw();
	VBFcriteria_2->Draw("same");
	VBFcriteria_3->Draw("same");
	leg->Draw();

	c_VBFcriteria_reco->cd();
	VBFcriteria_Reco1->Draw();
	VBFcriteria_Reco2->Draw("same");
	VBFcriteria_Reco3->Draw("same");
	leg->Draw();

	c_VBFcriteria->cd();
	VBFcriteria_1->Draw();
	VBFcriteria_2->Draw("same");
	VBFcriteria_3->Draw("same");
  VBFcriteria_Reco1->Draw("same");
  VBFcriteria_Reco2->Draw("same");
  VBFcriteria_Reco3->Draw("same");
	leg_all->Draw();


	
	//Write the canvas
	f_overlay->cd();
	c_dEta->Write();
	c_dEta_gen->Write();
	c_dEta_reco->Write();
	c_VBFcriteria->Write();
	c_VBFcriteria_gen->Write();
	c_VBFcriteria_reco->Write();
	
     //c_njet->SaveAs(path+"beauty/njet.pdf");
     //c_pt0->SaveAs(path+"beauty/pt0.pdf");
     //c_pt1->SaveAs(path+"beauty/pt1.pdf");
     //c_eta0->SaveAs(path+"beauty/eta0.pdf");
     //c_eta1->SaveAs(path+"beauty/eta1.pdf");
	//exit(0);

}
