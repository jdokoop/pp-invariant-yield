//-----------------------------------------------
// Macro to compute the BBC z-vertex resolution
// by fitting the distribution to a Gaussian.
// This macro compares the resolution of data
// and simulated p+p events.
//
// J.O. Koop
// April 22 2016
//-----------------------------------------------

#include <iostream>

using namespace std;

//---------------------------------
// Variables
//---------------------------------

//NTuples with event information
TTree *ntpEventData;
TTree *ntpEventSims;

//Distribution of BBC z-vertices in data and simulations
TH1F *hBBCzData;
TH1F *hBBCzSims;

//---------------------------------
// Functions
//---------------------------------

void plotResolutions()
{
	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	c1->SetLogy();
	hBBCzData->SetTitle("DATA");
	hBBCzData->GetXaxis()->SetTitle("BBC z-vertex [cm]");
	hBBCzData->GetXaxis()->SetTitleFont(62);
	hBBCzData->GetXaxis()->SetLabelFont(62);
	hBBCzData->GetYaxis()->SetTitle("Counts");
	hBBCzData->GetYaxis()->SetTitleFont(62);
	hBBCzData->GetYaxis()->SetLabelFont(62);
	hBBCzData->Draw();
	TF1 *f = (TF1*) hBBCzData->GetFunction("fit");
	f->Draw("same");
	TLatex *tl1 = new TLatex(0.4, 0.2, Form("Resol. = %g cm", f->GetParameter(2)));
	tl1->SetNDC(kTRUE);
	tl1->Draw("same");

	TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
	hBBCzSims->SetTitle("SIMULATIONS");
	hBBCzSims->GetXaxis()->SetTitle("BBC z-vertex [cm]");
	hBBCzSims->GetXaxis()->SetTitleFont(62);
	hBBCzSims->GetXaxis()->SetLabelFont(62);
	hBBCzSims->GetYaxis()->SetTitle("Counts");
	hBBCzSims->GetYaxis()->SetTitleFont(62);
	hBBCzSims->GetYaxis()->SetLabelFont(62);
	hBBCzSims->Draw();
}

void fitGaussian(TH1F *& h)
{
	float parHeight = h->GetMaximum();
	float parMean   = h->GetBinCenter(h->GetMaximumBin());
	float parWidth  = h->GetRMS();

	TF1 *fit = new TF1("fit", "gaus", -1.5 * parWidth, 1.5 * parWidth);

	fit->SetParameter(0, parHeight);
	fit->SetParameter(1, parMean);
	fit->SetParameter(2, parWidth);

	h->Fit(fit, "Q0R");
}

void ComputeBBCVtxResolution()
{
	//Read files and extract BBC z-vertex resolutions
	TFile *fData = new TFile("Data/423844_data_0_1.root");
	TFile *fSims = new TFile("Data/423844_reco_1_0_1_0.root");

	ntpEventData = (TTree*) fData->Get("ntp_event");
	ntpEventSims = (TTree*) fSims->Get("ntp_event");

	ntpEventData->Draw("vtx_bbc[2]>>hBBCzData(200,-100,100)", "", "goff");
	hBBCzData = (TH1F*) gDirectory->FindObject("hBBCzData");

	ntpEventSims->Draw("vtx_bbc[2]>>hBBCzSims(50,-25,25)", "", "goff");
	hBBCzSims = (TH1F*) gDirectory->FindObject("hBBCzSims");

	fitGaussian(hBBCzData);

	plotResolutions();
}