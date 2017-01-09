#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TGraphSmooth.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <regex>

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/RootUtils.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

using namespace std;

TObjArray toSave;
std::map<TString,TString> systForClosure;

void closureTest(TFile *gfile,string &distr,string &ch, string &cat, bool purePhoton);
std::vector<TH1F *> getRatioOnly(TFile *llF,TFile *gF,TString distr,TString ch, TString cat, bool purePhoton);
//void runFinalClosure();
//void runFinalHZZClosure();

bool purePhoton=false;
bool smoothFakesHisto=false;

string dilCh="ll";                                                                                                                           

string mode="MC";     

//
//int main()
int main(int argc,const char* argv[])
{
  systForClosure.clear();

  TFile *gInF;

	std::vector<string> dilCh={"ll", "ee", "mumu"};
	//std::vector<string> dilCh={"mumu"};

  //open the files with the input plots
  string gDataFile="plotter.root";
  //  string zDataFile="plotter.root";

  string outDir = "photonZclosure/";  
cout<< __LINE__ <<endl;
  
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")!=string::npos){
      printf("--help   --> print this helping text\n");
      printf("--inFile  --> path of the plotter .root file with inputs\n");
      printf("--outDir  --> path of the output directory to save the plots\n");
      printf("--mode  --> use MC or DATA to run the closure test in MC or Data ");
      return 0;
    }
    if(arg.find("--inFile" )!=string::npos && i+1<argc){ gDataFile= argv[i+1];  gDataFile= argv[i+1]; i++; printf("input file = %s\n", gDataFile.c_str()); }
    if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir = argv[i+1];  i++; printf("outDir = %s\n", outDir.c_str());  }
    if(arg.find("--mode" )!=string::npos && i+1<argc) { mode = argv[i+1];  i++; printf("mode = %s\n", mode.c_str());  }  
  }
  
  //  TFile *llInF=TFile::Open(gDataFile.c_str());
  gInF=TFile::Open(gDataFile.c_str());

  //std::vector<string> distr = {"met","mt","axialmet","mindphijmet","balance", "hugoTest_pt0_30smaller50_qt", "hugoTest_pt0_30smaller40_qt", "hugoTest_dphi_j_boson", "hugoTest_MET_50_80_zmass", "hugoTest_MET_0_80_zmass", "hugoTest_MET_0_30_zmass", "hugoTest_MET_30_50_zmass", "hugoTest_MET_80_120_zmass", "hugoTest_MET_80_100_zmass", "hugoTest_dphijboson_3_inf_met", "hugoTest_dphijboson_3_3p5_met", "hugoTest_dphijboson_2_inf_met", "hugoTest__pt0_70smaller100_met", "hugoTest__pt0_50smaller90_met", "hugoTest__pt0_30smaller50_met"};
 // std::vector<string> distr = {"met","mt","axialmet","mindphijmet","balance", "hugoTest_dphi_j_boson", "hugoTest_dphijboson_3_inf_met", "hugoTest_dphijboson_3_3p5_met", "hugoTest_dphijboson_2_inf_met", "hugoTest_dphijboson_2p5_3_met", "hugoTest_dphijboson_2_2p5_met", "hugoTest_dphijboson_1p5_2_met", "hugoTest_dphijboson_0_2_met", "hugoTest_dphijboson_0_1_met", "hugoTest_dphijboson_0_0p5_met", "hugoTest_dphijboson_1_1p5_met", "hugoTest_dphijboson_0p5_1_met", "hugoTest__pt0_90smaller200_met", "hugoTest__pt0_100smaller150_met", "hugoTest__pt0_300smaller500_met", "hugoTest__pt0_200smallerInf_met", "hugoTest__pt0_70smaller100_met", "hugoTest__pt0_50smaller90_met", "hugoTest__pt0_50smaller70_met", "hugoTest__pt0_200smaller300_met", "hugoTest__pt0_150smaller200_met", "hugoTest__pt0_30smaller50_met", "hugoTest__pt0_30smaller40_met", "hugoTest__pt0_40smaller50_met", "hugoTest__pt0_500smallerInf_met", "hugoTest__nbreJets_above_30__1_met", "hugoTest__nbreJets_above_30__2_met", "hugoTest__nbreJets_above_30__3_met", "hugoTest__nbreJets_above_30__4_met", "hugoTest__nbreJets_above_30__5_met", "hugoTest__nbreJets_above_30__6_met", "hugoTest__nbreJets_above_30__7_met", "hugoTest__nbreJets_above_15__1_met", "hugoTest__nbreJets_above_15__2_met", "hugoTest__nbreJets_above_15__3_met", "hugoTest__nbreJets_above_15__4_met", "hugoTest__nbreJets_above_15__5_met", "hugoTest__nbreJets_above_15__6_met", "hugoTest__nbreJets_above_15__7_met"};
  //std::vector<string> distr = {"met","mt","axialmet","mindphijmet","balance"};
  std::vector<string> distr = {"met"};
  std::vector<string> cat = {"eq0jets","geq1jets","vbf", ""};
  //std::vector<string> cat = {"geq1jets"};
cout<< __LINE__ <<endl;

for( unsigned int lep=0; lep < dilCh.size(); lep++){
  for(unsigned int icat=0; icat<cat.size(); icat++)
    {
      for(unsigned int ich=0; ich<distr.size(); ich++) 
	{
	/*	if(cat[icat] == "eq0jets" && (distr[ich].find("hugoTest__pt0_") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above") != std::string::npos) ) continue;
	  if(cat[icat] == "vbf" && (distr[ich].find("hugoTest__nbreJets_above_30__1") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_30__5") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_30__6") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_30__7") != std::string::npos ) )  continue ;
	  if(cat[icat] == "vbf" && (distr[ich].find("hugoTest__nbreJets_above_15__1") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_15__5") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_15__6") != std::string::npos || distr[ich].find("hugoTest__nbreJets_above_15__7") != std::string::npos ) )  continue ;
	  */
	  std::cout << "Running for distr : " <<  distr[ich] << " , event flavor : " << dilCh[lep] << 
	    " and category : " << cat[icat] << std::endl;
cout<< __LINE__ <<endl;

closureTest(gInF,distr[ich],dilCh[lep],cat[icat],purePhoton);
	  
	}
   	} }

  system( (string("mkdir -p ") + outDir).c_str());
  system( (string(" mv ") + string("*closure*.* ") + outDir).c_str());
  system( (string(" mv ")  + string("*syst*.* ") + outDir ).c_str());
  system( (string(" mv ")  + outDir +string("/*.root . ")).c_str());
  
  //close all opened files
  //  llInF->Close();
  gInF->Close();
  toSave.Clear();
cout<< __LINE__ <<endl;

  return 0;
}

//
std::vector<TH1F *> getRatioOnly(TFile *llF,TFile *gF,TString distr,TString ch, TString cat, bool purePhoton)
{
  std::vector<TH1F *> toReturn;

  bool rebin(distr.Contains("jetdeta") || distr.Contains("spt"));

  //
  //GET HISTOS FROM FILES
  //
  TString mcdy("Ztoll_M-50");
  TH1 *hdy = 0;
  if(ch=="ll")
    {
      hdy=(TH1*) llF->Get(mcdy+"/ee"+cat+"_"+distr);
      if(hdy)
	{
	  hdy= (TH1 *)hdy->Clone();
	  hdy->Add( (TH1 *) llF->Get(mcdy+"/mumu"+cat+"_"+distr) );
	}
    }
  else
    {
      hdy=(TH1 *) llF->Get(mcdy+"/"+ch+cat+"_"+distr);
      if(hdy) hdy=(TH1 *)hdy->Clone();
    }
  if(hdy==0) return toReturn;
  hdy->SetDirectory(0);
  if(rebin)hdy->Rebin();


  std::vector<TString> mcg;
  mcg.push_back("#gamma+jets");
  if(!purePhoton) mcg.push_back("Multijets");
  TH1 *hg=0;
  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      if(ch=="ll")
	{
	  if(hg) {
	    hg->Add( (TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+cat+"_"+distr);
	  }
	  hg->Add( (TH1 *)gF->Get(mcg[ig]+"/ee"+cat+"_"+distr) );
	  
	}
      else
	{
	  if(hg){
	    hg->Add( (TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+ch+cat+"_"+distr);
	  }
	}
    }
  hg->SetDirectory(0);
  if(rebin)hg->Rebin();
  
  //compute the ratio
  TH1F *hratio=(TH1F *)hdy->Clone(distr+"_ratio");
  hratio->Divide(hg);
  toReturn.push_back((TH1F *)hdy);
  toReturn.push_back((TH1F *)hg);
  toReturn.push_back((TH1F *)hratio);
  return toReturn;
}


//
void closureTest(TFile *gF,string &distr,string &ch,string &cat, bool purePhoton)
{
cout<< __LINE__ <<endl;
	//  bool rebin(distr.Contains("jetdeta") || distr.Contains("spt") || distr.Contains("qgmva"));
  bool rebin=false;
  //
  //GET HISTOS FROM FILES
  //
  string mcdy;
  //if (mode=="MC") { mcdy="Z#rightarrow ee-#mu#mu_filt1113"; }
  if (mode=="MC") { mcdy="Z#rightarrow ee-#mu#mu_filt1113_dataZpt_reweighted"; }
  else if (mode=="DATA") { mcdy = "data"; }
cout<< __LINE__ <<endl;
  
  TH1D *hdy = NULL;
  if(ch=="ll")
    {
      hdy=(TH1D *) gF->Get( (mcdy+"/ee"+cat+"_"+distr).c_str() );
cout<< __LINE__ <<endl;
      hdy=(TH1D *) hdy->Clone( ("mcdy_"+cat+"_"+distr).c_str() );
cout<< __LINE__ <<endl;
      hdy->Add((TH1D *) gF->Get( (mcdy+"/mumu"+cat+"_"+distr).c_str()) );
    }
  else
    {
      hdy=(TH1D *) gF->Get( (mcdy+"/"+ch+cat+"_"+distr).c_str() );
      hdy=(TH1D *) hdy->Clone( ("mcdy_"+ch+cat+"_"+distr).c_str() ); 
    }
  if(hdy==0) return;
  hdy->SetDirectory(0);
  if(rebin) hdy->Rebin();

  // bool useFakes=true;
  
  std::vector<string> mcg;
cout<< __LINE__ <<endl;

  //  if (mode=="MC") { mcg.push_back("#gamma+jets_mc_reweighted"); } // incl. DY
  if (mode=="MC") { mcg.push_back("#gamma+jets_reweighted"); }  // DY in HTbins
  if (mode=="MC") { mcg.push_back("#gamma+jets_reweighted"); }  // DY in HTbins
  //else if (mode=="DATA") { mcg.push_back("Instr. MET"); }
  else if (mode=="DATA") { mcg.push_back("#gamma data_reweighted"); } //Here I don't play with Genuine MET so gamma reweighted is fine (Hugo)
  
cout<< __LINE__ <<endl;
  purePhoton=true;
	//purePhoton=false; 

  if(!purePhoton) { // Add EWK
    // mcg.push_back("W#rightarrow l#nu");
    // mcg.push_back("Top");
    // mcg.push_back("Z#gamma #rightarrow ll#gamma");
    // mcg.push_back("W#gamma #rightarrow l#nu#gamma");
    // mcg.push_back("Z#rightarrow #nu#nu");
    //mcg.push_back("QCD, HT>100"); //Should not be commented right?
    //mcg.push_back("QCD_EMEnr"); //Should not be commented right?
  }
 
  TH1D *hg=NULL, *hpureg=NULL;
  TH1D *hn=NULL;

  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      if(ch=="ll")
	{// gamma selection
	  if(hg) { // include impurity of photon events
	    hg->Add( (TH1D *)gF->Get( (mcg[ig]+"/gamma"+cat+"_"+distr).c_str()) );
	    hn->Add( (TH1D *)gF->Get( (mcg[ig]+"/gamma"+cat+"_"+distr).c_str()) );
	  }
	  else{ // 1st gamma plot
	    hg=(TH1D *)gF->Get( (mcg[ig]+"/gamma"+cat+"_"+distr).c_str() );
	    std::cout << "test2 " << ("mcg_"+cat+"_"+distr).c_str()  << std::endl;

	    hg=(TH1D *)hg->Clone( (cat+"mcg_"+cat+"_"+distr).c_str() );
	    hn=(TH1D *)hg->Clone( (cat+"mcg_"+cat+"_"+distr).c_str() ); hn->Reset();
	  }
	  // hg->Add( (TH1D *)gF->Get(mcg[ig]+"/ee"+cat+"_"+distr) );
	
	}
      else
	{
cout<< __LINE__ <<endl;
	  if(hg){
	    hg->Add( (TH1D *) gF->Get( (mcg[ig]+"/"+ch+cat+"_"+distr).c_str()) );
	    //	    hn->Add( (TH1D *)gF->Get( (mcg[ig]+"/gamma"+ch+cat+"_"+distr).c_str()) ); 
	    hn->Add( (TH1D *)gF->Get( (mcg[ig]+"/"+ch+cat+"_"+distr).c_str()) ); //Add by Hugo 
	  }
	  else{
	    std::cout << "test3: " << ("mcg_"+ch+cat+"_"+distr).c_str()  << std::endl; 
cout<< __LINE__ <<endl;
	    hg=(TH1D *) gF->Get( (mcg[ig]+"/"+ch+cat+"_"+distr).c_str() );
cout<< __LINE__ <<endl;
cout << "mcg_"+ch+cat+"_"+distr << endl;
	    hg=(TH1D *) hg->Clone(("mcg_"+ch+cat+"_"+distr).c_str());
cout<< __LINE__ <<endl;
	    hn=(TH1D *) hg->Clone(("mcg_"+ch+cat+"_"+distr).c_str()); hn->Reset(); //Add by Hugo
	  }
	}
cout<< __LINE__ <<endl;
      
      if(ig==0){
cout<< __LINE__ <<endl;
	TString pureName(hg->GetName());
	pureName.ReplaceAll("mcg","mcpureg");
	hpureg=(TH1D *)hg->Clone(pureName);
      }
    }
cout<< __LINE__ <<endl;
  if(hg==0 || hpureg==0) return;
cout<< __LINE__ <<endl;
  hg->SetDirectory(0);
cout<< __LINE__ <<endl;
  hpureg->SetDirectory(0);
cout<< __LINE__ <<endl;
  hn->SetDirectory(0);
cout<< __LINE__ <<endl;
  if(rebin) {
    hg->Rebin();
    hpureg->Rebin();
  }
  //if(distr=="qt") { hg->Rebin(4); hdy->Rebin(4); hpureg->Rebin(4); }
  //  else if(!distr.Contains("cjv")) { hg->Rebin(2); hdy->Rebin(2); hpureg->Rebin(2); }

cout<< __LINE__ <<endl;
//find limits
  float dyscale(hdy->Integral());
  // hdy->Scale(1./hdy->Integral());
  Double_t xmin(hdy->GetXaxis()->GetXmin());
  Double_t xmax(hdy->GetXaxis()->GetXmax());
  if( (distr=="met") && (distr!="axial")) {xmin=0;   xmax=600;}
  if( (distr=="mt") )  {xmin=120; xmax=1500;}
  float ymin(3e-5),ymax(hdy->GetMaximum()*4.0);

  //draw
  float gscale(hg->Integral());
  hg->Scale(dyscale/gscale);
  
  hpureg->Scale(dyscale/gscale);
  // hg->Scale(1./gscale);
  // hpureg->Scale(1./gscale);
  float nscale(hn->Integral());
  hn->Scale(dyscale/(gscale-nscale));


  // Estimate systematic
  // Print number of events for DY sample in MET>125:
  //Int_t mbinlast = hdy->GetXaxis()->FindBin(125.);
  Int_t mbinlast = hdy->FindLastBinAbove();
	cout << "LastBin" << endl;
	cout << mbinlast << endl;
	cout << hdy->GetXaxis()->GetBinCenter(mbinlast) << endl;


  Int_t mbin;
  Double_t ndy,ng;
  Double_t errdy,errg;

  TH1D *hsyst=0;
  TH1D *hevt1=0, *hevt2=0;
  hevt1=(TH1D*)hdy->Clone("hevt1"); hevt1->Reset();
  hevt2=(TH1D*)hdy->Clone("hevt2"); hevt2->Reset();
  hsyst=(TH1D*)hdy->Clone("hsyst"); hsyst->Reset();

  Double_t syst, esyst;
  Double_t diff, nmax, enmax;
  
cout<< __LINE__ <<endl;
  for (size_t i=1; i<=mbinlast; i++) {
  
    ndy=hdy->IntegralAndError(i,-1,errdy);
    ng=hg->IntegralAndError(i,-1,errg);

    hevt1->SetBinContent(i,ndy); hevt1->SetBinError(i,errdy);
    hevt2->SetBinContent(i,ng); hevt2->SetBinError(i,errg);
    //std::cout << "Nevents DY for MET>125 (bin = "<< mbin << ") = " << ndy << " +/-" << errdy << std::endl;
    //std::cout << "Nevents gamma for MET>125 (bin = "<< mbin << ") = " << ng << " +/-" << errg << std::endl;
    diff=(ndy-ng);
    //nmax=max(ndy,ng);
    nmax=ng;
    //enmax=max(errdy,errg);
    enmax=errg;
    
    //syst=fabs(diff)/nmax;
    syst=diff/nmax;
    //esyst=((errdy*errdy + errg*errg)/(nmax*nmax)) + ((enmax*enmax)*(diff*diff)/(nmax*nmax*nmax*nmax));
    esyst=((errdy*errdy)/(nmax*nmax)) + ((enmax*enmax)*(ndy*ndy)/(nmax*nmax*nmax*nmax));
    // esyst=sqrt(esyst);
    
    hsyst->SetBinContent(i,syst);
    hsyst->SetBinError(i,sqrt(esyst));
    //std::cout << "Systematic assigned = " << syst << std::endl;
  }

cout<< __LINE__ <<endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  // SHOW SYSTEMATIC
  TCanvas *cs = new TCanvas( ("cs"+cat+"_"+distr,"cs"+cat+"_"+distr).c_str() );
  cs->SetWindowSize(750,900);
  cs->cd();
  
  //distributions
  TPad *st1 = new TPad("sp1","sp1",0,0.3,1.0,1.0);
  st1->Draw();
  st1->cd();
  st1->SetTopMargin(0.08);
  //st1->SetBottomMargin(0);
  st1->SetRightMargin(0.05);

  st1->SetGridx(); st1->SetGridy();
  
  //hevt1->Draw("E"); hevt2->Draw("ESAME");

  hsyst->GetYaxis()->SetLabelSize(0.06);
  hsyst->GetYaxis()->SetTitleSize(0.05);      
  hsyst->GetYaxis()->SetTitleOffset(1.0);
  hsyst->GetYaxis()->SetLabelSize(0.05);
  hsyst->GetXaxis()->SetRangeUser(0.,200.); 
  hsyst->GetXaxis()->SetLabelSize(0.06);
  hsyst->GetXaxis()->SetTitleSize(0.05);      
  hsyst->GetXaxis()->SetTitleOffset(1.0);
  hsyst->GetXaxis()->SetLabelSize(0.005);
  // hsyst->GetYaxis()->SetRangeUser(ymin,ymax);
  hsyst->GetYaxis()->SetTitle("#delta N / N");
  hsyst->SetLineColor(1);
  hsyst->SetMarkerColor(2);
  hsyst->SetMarkerStyle(1);
  //hsyst->SetFillStyle(3001);
  // hsyst->SetFillColor(kGreen-10);
cout<< __LINE__ <<endl;
  
  hsyst->Draw("EP");
  hsyst->SetMarkerColor(9);
  hsyst->SetMarkerStyle(20);
  hsyst->SetMarkerSize(0.6);
  hsyst->GetYaxis()->SetRangeUser(-2.,2.);

  hsyst->GetXaxis()->SetTitle("MET threshold [GeV]");//hdy->GetXaxis()->GetTitle());
  //hsyst->GetXaxis()->SetTitle("MET [GeV]");//hdy->GetXaxis()->GetTitle());
  hsyst->GetXaxis()->SetLabelSize(0.12);
  hsyst->GetXaxis()->SetTitleSize(0.14);
  hsyst->GetXaxis()->SetTitleOffset(0.8);
  
  cs->Modified();
  cs->Update();

  cs->SaveAs( (ch+cat+"_"+distr+"_syst.png").c_str());
  cs->SaveAs( (ch+cat+"_"+distr+"_syst.C").c_str() );
  cs->SaveAs( (ch+cat+"_"+distr+"_syst.pdf").c_str() );

  //return;
  //
  //SHOW CLOSURE TEST
  //
  
 
cout<< __LINE__ <<endl;

  TCanvas *c = new TCanvas( ("c"+cat+"_"+distr,"c"+cat+"_"+distr).c_str() );
  c->SetWindowSize(500,600);
  c->cd();
  
  //distributions
  TPad *t1 = new TPad("p1","p1",0,0.3,1.0,1.0);
  t1->Draw();
  t1->cd();
  t1->SetTopMargin(0.08);
  t1->SetBottomMargin(0);
  t1->SetRightMargin(0.05);

  
  
  TH1D *hfakes=(TH1D *)hn->Clone("fakes");
  hfakes->SetDirectory(0);
  if(smoothFakesHisto) hfakes->Smooth();


  
  hg->SetTitle("#gamma + jets");
  hg->GetYaxis()->SetLabelSize(0.06);
  hg->GetYaxis()->SetTitleSize(0.05);      
  hg->GetYaxis()->SetTitleOffset(1.0);
  hg->GetYaxis()->SetLabelSize(0.05);
  hg->GetXaxis()->SetRangeUser(xmin,xmax); 
  if(distr=="qt")
    {
      ymin=3e-6;
      t1->SetLogy();
      t1->SetLogx();
      hg->GetXaxis()->SetRangeUser(50,1e3);
    }
  hg->GetYaxis()->SetRangeUser(ymin,ymax);
  hg->GetYaxis()->SetTitle("Events (a.u.)");
  hg->SetLineColor(1);
  hg->SetMarkerColor(1);
  hg->SetMarkerStyle(1);
  hg->SetFillStyle(3001);
  hg->SetFillColor(kGreen-10);
  hg->Draw("ehist");


  if(!purePhoton){
    //    hfakes->Add(hpureg,-1);
    hfakes->SetTitle("Fakes");
    hfakes->SetLineColor(1);
    hfakes->SetMarkerColor(1);
    hfakes->SetMarkerStyle(1);
    hfakes->SetFillStyle(3001);
    hfakes->SetFillColor(kGray);
    hfakes->Draw("histsame");   
  }

cout<< __LINE__ <<endl;

  hdy->SetTitle("DY #rightarrow ll");
  hdy->Draw("esame");
  hdy->SetMarkerColor(9);
  hdy->SetMarkerStyle(20);
  hdy->SetMarkerSize(0.8);
  
  bool setLogY(false);
  TString mjjCat("M_{jj}>1000");
  if(cat=="mjjq016") mjjCat="M_{jj}<250";
  if(cat=="mjjq033") mjjCat="250<M_{jj}<350";
  if(cat=="mjjq049") mjjCat="350<M_{jj}<450";
  if(cat=="mjjq066") mjjCat="450<M_{jj}<550";
  if(cat=="mjjq083") mjjCat="550<M_{jj}<750";
  if(cat=="mjjq092") mjjCat="750<M_{jj}<1000";
  if(cat=="mjjgt092") mjjCat="M_{jj}>750";
  if(cat=="highmjj") mjjCat="M_{jj}>1250";
  if(cat=="highhardpt") mjjCat="Hard p_{T}>50";
  if(cat=="lowhardpt") mjjCat="Hard p_{T}<50";
  if(cat=="") { mjjCat="inclusive"; setLogY=true; }
  if(cat=="eq0jets")  {mjjCat="=0 jets"; setLogY=true;}
  if(cat=="geq1jets") {mjjCat="#geq1 jets"; setLogY=true;}
  if(cat=="vbf")      {mjjCat="VBF"; setLogY=true;}
  if(setLogY)t1->SetLogy();

  TPaveText *pave = new TPaveText(0.7,0.85,0.95,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(32);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->SetTextColor(kBlue);
  pave->AddText("["+mjjCat+"]");
  pave->Draw();

  TLegend *leg=new TLegend(0.6,0.6,0.9,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  //leg->SetNColumns(3);
  leg->SetNColumns(2);
  if(!setLogY)
    {
      leg->AddEntry(hdy,hdy->GetTitle(),"P");
      leg->AddEntry(hg,hg->GetTitle(),"F");
      if(!purePhoton) {
	leg->AddEntry(hfakes,hfakes->GetTitle(),"F");
	leg->SetNColumns(3);
      }
    }
  else
    {
      leg->AddEntry(hdy,"Z","P");
      leg->AddEntry(hg,"#gamma","F");
    }
  leg->Draw("same");


cout<< __LINE__ <<endl;
  pave = new TPaveText(0.94,0.4,1.0,0.83,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(21);
  pave->SetTextFont(42);
  pave->SetTextColor(kGray+2);
  pave->SetTextSize(0.05);
  char buf[1000];  
  sprintf(buf,"#chi^{2}/ndof : %3.2f , K-S prob : %3.2f",hdy->Chi2Test(hg,"WWCHI2/NDF"),hdy->KolmogorovTest(hg,"") );
  pave->AddText(buf)->SetTextAngle(-90);
  pave->Draw();

  pave = new TPaveText(0.09,0.95,0.35,0.98,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(12);
  pave->SetTextSize(0.05);
  pave->AddText("CMS simulation, #sqrt{s}=13 TeV");
  pave->Draw();

  //closure
  c->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->SetRightMargin(0.05);
  t2->SetGridx(); t2->SetGridy();
  t2->Draw();
  t2->cd();

  leg = new TLegend(0.6,0.52,0.9,0.92,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(3001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.11);
  leg->SetTextAlign(12);
cout<< __LINE__ <<endl;


  //mc stats (of the photon sample here...)
  TH1D *denRelUncH=(TH1D *) hg->Clone("mcrelunc");
  for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
    {
      if(denRelUncH->GetBinContent(xbin)==0) continue;
      Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
      denRelUncH->SetBinContent(xbin,1);
      denRelUncH->SetBinError(xbin,err);
    }
  TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
  denRelUnc->SetLineColor(1);
  denRelUnc->SetFillStyle(3001);
  denRelUnc->SetFillColor(kGray);
  denRelUnc->SetMarkerColor(1);
  denRelUnc->SetMarkerStyle(1);
  denRelUncH->Reset("ICE");       
  denRelUncH->Draw();
  denRelUnc->Draw("3");
  leg->AddEntry(denRelUnc,"stat unc.","f");
  denRelUncH->GetYaxis()->SetRangeUser(-0.5,3.);
   //  denRelUncH->GetYaxis()->SetRangeUser(-0.2,1.74);
  denRelUncH->GetXaxis()->SetTitle(hdy->GetXaxis()->GetTitle());
  denRelUncH->GetXaxis()->SetLabelSize(0.12);
  denRelUncH->GetXaxis()->SetTitleSize(0.14);
  denRelUncH->GetXaxis()->SetTitleOffset(0.8);
  denRelUncH->GetYaxis()->SetLabelSize(0.12);
  denRelUncH->GetYaxis()->SetNdivisions(5);
  denRelUncH->GetYaxis()->SetTitleSize(0.12);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  denRelUncH->GetYaxis()->SetTitle("Ratio");
  denRelUncH->GetYaxis()->SetTitleOffset(0.4);
  denRelUncH->SetMarkerStyle(20);
  denRelUncH->SetMarkerColor(1);
  denRelUncH->SetLineColor(1);
  if(distr=="qt") { t2->SetLogx(); denRelUncH->GetXaxis()->SetRangeUser(50,1000); }

cout<< __LINE__ <<endl;
   
  //now the actual ratio 
  TH1D *gr = (TH1D *) hdy->Clone("closuregr");
  gr->SetDirectory(0);
  gr->Divide(hg);
 
  //smooth curve and symmetrize it
  TGraphErrors *uncGr=new TGraphErrors(gr);   
  TGraph *uncGrUp=new TGraph(gr); uncGrUp->SetName( (ch+cat+"_"+distr+"_unc").c_str());  uncGrUp->SetLineWidth(2);
  TGraph *uncGrDown=new TGraph;  uncGrDown->SetLineWidth(2);
  for(int ip=0; ip<uncGrUp->GetN(); ip++)
    {
      Double_t x,y;
       uncGr->GetPoint(ip,x,y);
      //  uncGrUp->GetPoint(ip,x,y);
      // uncGr->SetPoint(ip,x,1+TMath::Abs(y-1));
      // uncGrUp->SetPoint(ip,x,1+TMath::Abs(y-1));
      // uncGrDown->SetPoint(ip,x,1-TMath::Abs(y-1));
    }
  uncGr->SetMarkerColor(1); //uncGr->SetMarkerSize(0.3);
  uncGr->Draw("P");
 
  uncGrUp->SetLineWidth(1);
  uncGrDown->SetLineWidth(1);
  if(setLogY)   leg->AddEntry(uncGrUp,"Z/#gamma","l");
  else          leg->AddEntry(uncGrUp,"Z /#gamma","l");
  // uncGrUp->Draw("l");
  // uncGrDown->Draw("l");
  toSave.Add(uncGrUp);

  //return;
  
  leg->Draw();
  leg->SetNColumns(4);

  c->Modified();
  c->Update();

cout<< __LINE__ <<endl;
  c->SaveAs( (ch+cat+"_"+distr+"_closure.png").c_str());
  c->SaveAs( (ch+cat+"_"+distr+"_closure.C").c_str() );
  c->SaveAs( (ch+cat+"_"+distr+"_closure.pdf").c_str() );
}
