//macro to superimpose histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"

void overlay()
{
  //Author: Hugo 

  std::vector<string> file = {"data", "Z#rightarrow ee-#mu#mu_filt1113"};
	std::vector<string> cat   = {"mumu", "ee", "ll"};
  std::vector<string> bin   = {"", "eq0jets","geq1jets","vbf"};
  std::vector<string> var   = {"_MET_0_30","_MET_30_50","_MET_50_80","_MET_80_120"};

  TFile* dFile = new TFile("plotter.root", "READ");
  if(!dFile || dFile->IsZombie() || !dFile->IsOpen() || dFile->TestBit(TFile::kRecovered) ){
    printf("File plotter.root coudn't be opened\n");
  }


	for (unsigned int ifile = 0; ifile < file.size(); ifile++){
  	for (unsigned int icat = 0; icat < cat.size(); icat++){
			for (unsigned int ibin = 0; ibin < bin.size(); ibin++){

				std::vector<TH1D* > v_hist;
				for (unsigned int ivar = 0; ivar <var.size(); ivar++){
					v_hist.push_back((TH1D*) dFile->Get((file[ifile]+"/"+cat[icat]+bin[ibin]+"_hugoTest"+var[ivar]+"_zmass").c_str()));
				}
			}
		}
	}




}
