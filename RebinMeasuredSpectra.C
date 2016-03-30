#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

//Number of original pT bins
const int NBINS = 9;

//Number of bins in rebinned spectra
const int NBINS_RB = 12;

//Number of sectors
const int NSECT = 4;

//Number of events in each file
int nevents_truth = 447000;
int nevents_reco  = 0;
int nevents_data  = 0;

//Raw pT yield
TH1F *h_DeltaNDeltapT_truth;
TH1F *h_DeltaNDeltapT_reco;
TH1F *h_DeltaNDeltapT_data;

//Rebinned raw pT yield
TH1F *h_DeltaNDeltapT_truth_rebinned;
TH1F *h_DeltaNDeltapT_reco_rebinned;
TH1F *h_DeltaNDeltapT_data_rebinned;

//Acceptance and Efficiency Correction
TH1F *h_correction;

string zvtxcut  = "TMath::Abs(vtx[2]) < 5";
string chisqcut = "chisq/ndf < 3";
string dcacut   = "TMath::Abs(dca) < 0.15";
string dca2dcut = "TMath::Abs(dca2d) < 0.05";
string momcut   = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.2";
string nhitscut = "nhits[0]+nhits[1]+nhits[2]+nhits[3] == 4";
string eta1cut = "TMath::Abs(eta) < 0.35";
string eta2cut = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 0.35";

//Sector cuts
string etcut    = "TMath::ATan2(mom[1],mom[0]) < -2";
string ebcut    = "TMath::ATan2(mom[1],mom[0]) > 2";
string wtcut    = "TMath::ATan2(mom[1],mom[0]) > 0 && TMath::ATan2(mom[1],mom[0]) < 2";
string wbcut    = "TMath::ATan2(mom[1],mom[0]) < 0 && TMath::ATan2(mom[1],mom[0]) > -2";

string sectors[NSECT] = {"ET", "EB", "WT", "WB"};
string sectorCuts[NSECT] = {etcut, ebcut, wtcut, wbcut};

//---------------------------------------
// Functions
//---------------------------------------

void readFiles()
{
	//Read AMPT truth information
	TFile *f_true = new TFile("Data/ampt_pp_true.root");
	TTree *ntp_svxseg_true = (TTree*) f_true->Get("ntp_svxseg_true");
	ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_truth(%i,0.2,2)", NBINS), (momcut + "&&" + eta1cut).c_str(), "goff");
	h_DeltaNDeltapT_truth = (TH1F*) gDirectory->FindObject("h_DeltaNDeltapT_truth");

	//Read reconstructed AMPT information
	TFile *f_reco    = new TFile("Data/423844_ampt_smeared_105.root");
	TTree *ntp_svxseg_reco  = (TTree*) f_reco->Get("ntp_svxseg");
	TTree *ntp_event_reco = (TTree*) f_reco->Get("ntp_event");
	ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco(%i,0.2,2)", NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut).c_str(), "goff");
	h_DeltaNDeltapT_reco = (TH1F*) gDirectory->FindObject("h_DeltaNDeltapT_reco");
	nevents_reco = ntp_event_reco->GetEntries();

	//Read actual data
	TFile *f_data  = new TFile("Data/423844_data_105.root");
	TTree *ntp_svxseg_data  = (TTree*) f_data->Get("ntp_svxseg");
	TTree *ntp_event_data   = (TTree*) f_data->Get("ntp_event");
	ntp_svxseg_data->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_data(%i,0.2,2)", NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut).c_str() , "goff");
	h_DeltaNDeltapT_data = (TH1F*) gDirectory->FindObject("h_DeltaNDeltapT_data");
	nevents_data = ntp_event_data->GetEntries();

	//Set errors on these histograms since they just contain counts
	for(int i=1; i<=NBINS; i++)
	{
		h_DeltaNDeltapT_truth->SetBinError(i, TMath::Sqrt(h_DeltaNDeltapT_truth->GetBinContent(i)));
		h_DeltaNDeltapT_reco->SetBinError(i, TMath::Sqrt(h_DeltaNDeltapT_reco->GetBinContent(i)));
		h_DeltaNDeltapT_data->SetBinError(i, TMath::Sqrt(h_DeltaNDeltapT_data->GetBinContent(i)));
	}
}

void rebinHistograms()
{
	double newBins[NBINS_RB+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2.0};

	h_DeltaNDeltapT_truth_rebinned = (TH1F*) h_DeltaNDeltapT_truth->Rebin(NBINS_RB, "h_DeltaNDeltapT_truth_rebinned", newBins);
	h_DeltaNDeltapT_reco_rebinned = (TH1F*) h_DeltaNDeltapT_reco->Rebin(NBINS_RB, "h_DeltaNDeltapT_reco_rebinned", newBins);
	h_DeltaNDeltapT_data_rebinned = (TH1F*) h_DeltaNDeltapT_data->Rebin(NBINS_RB, "h_DeltaNDeltapT_data_rebinned", newBins);
}

void divideHistos()
{
	h_correction = (TH1F*) h_DeltaNDeltapT_reco->Clone("h_correction");
	h_correction->Divide(h_DeltaNDeltapT_truth);
}

void plot()
{
	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1","Truth",600,600);
	h_DeltaNDeltapT_truth->SetLineColor(kBlack);
	h_DeltaNDeltapT_truth->Draw();

	TCanvas *c2 = new TCanvas("c2","Reco",600,600);
	h_DeltaNDeltapT_reco->SetLineColor(kBlue);
	h_DeltaNDeltapT_reco->Draw();

	TCanvas *c3 = new TCanvas("c3","Data",600,600);
	h_DeltaNDeltapT_data->SetLineColor(kRed);
	h_DeltaNDeltapT_data->Draw();

	TCanvas *c4 = new TCanvas("c4","Acceptance and Efficiency Correction");
	h_correction->Draw();
}

void RebinMeasuredSpectra()
{
	readFiles();
	//rebinHistograms();
	divideHistos();
	plot();
}