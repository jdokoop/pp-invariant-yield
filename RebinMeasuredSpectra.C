#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

//Number of original pT bins
const int NBINS = 18;

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

string zvtxcut  = "TMath::Abs(vtx[2]) < 1";
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
}

void RebinMeasuredSpectra()
{
	readFiles();
}