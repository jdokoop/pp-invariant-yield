#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

const int NBINS = 18;

TGraph *g_resolution;

TH1F *h_dNdEta_smeared;
TH1F *h_dNdEta_true;
TH1F *h_ratio;

//---------------------------------------
// Functions
//---------------------------------------

void fitResolution()
{
	//Take pT resolution for 4-hit tracks from Theo
	//Fit the points to a polynomial

	float pT[4] = {0.25, 0.76, 1.25, 1.75};
	float res[4] = {0.109, 0.110, 0.113, 0.119};

	g_resolution = new TGraph(4, pT, res);

	TF1 *fit = new TF1("fit","pol3", 0.25, 1.8);
	fit->SetParameter(1, 0.110045);
	fit->SetParameter(2, -0.00601476);
	fit->SetParameter(3, 0.00840499);
	fit->SetParameter(4, -0.00120621);

	g_resolution->Fit("fit","Q0R");
	g_resolution->SetMarkerStyle(20);
}

void smearTruth()
{
	TFile *fTruth = new TFile("Data/ampt_pp_true.root");
	TTree *ntp_svxseg_true = (TTree*) fTruth->Get("ntp_svxseg_true");

	//Extract true dNdpT
	h_dNdEta_true = new TH1F("h_dNdEta_true","h_dNdEta_true",NBINS,0.2,2.0);
	h_dNdEta_smeared = new TH1F("h_dNdEta_smeared","h_dNdEta_smeared",NBINS,0.2,2.0);

	h_dNdEta_true->Sumw2();
	h_dNdEta_smeared->Sumw2();

	//Initialize random number generator
	gRandom = new TRandom3();

	float mom[3];
	float eta;
	ntp_svxseg_true->SetBranchAddress("mom", &mom);
	ntp_svxseg_true->SetBranchAddress("eta", &eta);

	for(int i=1; i<= ntp_svxseg_true->GetEntries(); i++)
	{
		ntp_svxseg_true->GetEntry(i);

		float pT = TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]);
		float delta = gRandom->Gaus(0, g_resolution->Eval(pT));
		float pT_smeared = pT + delta;

		h_dNdEta_smeared->Fill(pT_smeared);
		h_dNdEta_true->Fill(pT);
	}

	h_ratio = (TH1F*) h_dNdEta_smeared->Clone("h_ratio");
	h_ratio->Divide(h_dNdEta_true);
}

void plot()
{
	gStyle->SetOptStat(0);
	TCanvas *cRatio = new TCanvas("cRatio","Ratio",600,600);
	h_ratio->GetXaxis()->SetTitle("");
	h_ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_ratio->GetYaxis()->SetTitle("Smeared/ Truth");
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("CP");

	TCanvas *cResolution = new TCanvas("cResolution", "Resolution", 600, 600);
	g_resolution->SetTitle("");
	g_resolution->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_resolution->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");
	TF1 *polFit = (TF1*) g_resolution->GetFunction("fit");
	g_resolution->Draw("AP");
	polFit->Draw("same");
}

void InvestigateResolutionEffect()
{
	fitResolution();
	smearTruth();
	plot();
}