#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

//Number of sectors
// --> Inclusive phi, ET, EB, WT, WB
const int NSECT = 5;

//Invariant yield
TH1F *h_dNdpT_truth[NSECT];
TH1F *h_dNdpT_reco[NSECT];
TH1F *h_dNdpT_data[NSECT];

//Acceptance and Efficiency Correction
TH1F *h_correction[NSECT];

//Label for each azimuthal sector
string sectorLabel[NSECT] = {"INCLUSIVE", "ET", "EB", "WT", "WB"};

//Colors for plotting yields
int yieldColor[NSECT] = {kBlack, kOrange + 1, kBlue, kRed, kGreen + 3};

//Tsallis fits to spectra in each sector
TF1 *tsallisSectorFit[5];

//Ratio of Tsallis fits to spectrum in each sector
TGraphErrors *tsallisSectorRatio[5];

//---------------------------------------
// Functions
//---------------------------------------

void readFiles()
{
	TFile *fin = new TFile("WorkingFiles/normalizedSpectra.root");

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i] = (TH1F*) fin->Get(Form("h_dNdpT_truth_%s", sectorLabel[i].c_str()));
		h_dNdpT_reco[i] = (TH1F*) fin->Get(Form("h_dNdpT_reco_%s", sectorLabel[i].c_str()));
		h_dNdpT_data[i] = (TH1F*) fin->Get(Form("h_dNdpT_data_%s", sectorLabel[i].c_str()));

		h_correction[i] = (TH1F*) fin->Get(Form("h_correction_%s", sectorLabel[i].c_str()));
	}
}

void applyAcceptanceCorrection()
{
	//Apply correction and recompute errors
	for (int i = 0; i < NSECT; i++)
	{
		for (int j = 1; j <= h_dNdpT_data[i]->GetNbinsX(); j++)
		{
			float dataContent = h_dNdpT_data[i]->GetBinContent(j);
			float correction = h_correction[i]->GetBinContent(j);

			float error = (dataContent / correction) * (h_dNdpT_data[i]->GetBinError(j) / dataContent + h_correction[i]->GetBinError(j) / correction);

			h_dNdpT_data[i]->SetBinError(j, error);
			h_dNdpT_data[i]->SetBinContent(j, dataContent / correction);
		}
	}
}

void fitSpectrum()
{
	//Fit spectra with Tsallis form
	float m = 0.15;
	TF1 *tsallisFit = new TF1("tsallisFit", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.140*([1] - 1))*([1]*[2] + 0.140)) * pow(([1]*[2] + TMath::Sqrt(0.140*0.140 + x*x))/([1]*[2]+0.140),-1*[1])", 0.3, 1.9);
	//tsallisFit->SetParameter(0, 1.28);
	//tsallisFit->SetParameter(1, 9.67);
	//tsallisFit->SetParameter(2, 121.077);

	tsallisFit->SetParameter(0, 0.78);
	tsallisFit->SetParameter(1, 38.79);
	tsallisFit->SetParameter(2, 0.145);

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_data[i]->Fit("tsallisFit", "Q0R");
		tsallisSectorFit[i] = (TF1*) h_dNdpT_data[i]->GetFunction("tsallisFit");
		tsallisSectorFit[i]->SetLineColor(yieldColor[i]);
	}
}

void divideFits()
{
	double x[9];
	double y[9];
	double ex[9] = {0};
	double ey[9] = {0};

	for (int j = 0; j < NSECT; j++)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] = h_dNdpT_data[j]->GetBinCenter(i + 1);

			y[i] = h_dNdpT_data[j]->GetBinContent(i + 1) / tsallisSectorFit[j]->Eval(x[i]);
			ey[i] = h_dNdpT_data[j]->GetBinError(i + 1) / tsallisSectorFit[j]->Eval(x[i]);
		}

		tsallisSectorRatio[j] = new TGraphErrors(9, x, y, ex, ey);
	}
}

void plotTruthReco()
{
	TCanvas *c0 = new TCanvas("c0", "AMPT True and Reco Spectra", 800, 800);
	c0->Divide(2, 2);

	for (int i = 1; i < NSECT; i++)
	{
		c0->cd(i);
		gPad->SetLogy();

		h_dNdpT_truth[i]->SetLineColor(kBlack);
		h_dNdpT_reco[i]->SetLineColor(kBlue);

		h_dNdpT_truth[i]->SetMarkerColor(kBlack);
		h_dNdpT_reco[i]->SetMarkerColor(kBlue);

		h_dNdpT_truth[i]->SetMarkerStyle(20);
		h_dNdpT_reco[i]->SetMarkerStyle(20);

		h_dNdpT_truth[i]->SetTitle(sectorLabel[i].c_str());

		h_dNdpT_truth[i]->GetXaxis()->SetTitle("p_{T}");
		h_dNdpT_truth[i]->GetXaxis()->SetTitleSize(0.05);
		h_dNdpT_truth[i]->GetXaxis()->SetLabelSize(0.05);
		h_dNdpT_truth[i]->GetXaxis()->SetTitleFont(62);
		h_dNdpT_truth[i]->GetXaxis()->SetLabelFont(62);
		h_dNdpT_truth[i]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");
		h_dNdpT_truth[i]->GetYaxis()->SetTitleSize(0.05);
		h_dNdpT_truth[i]->GetYaxis()->SetLabelSize(0.05);
		h_dNdpT_truth[i]->GetYaxis()->SetTitleFont(62);
		h_dNdpT_truth[i]->GetYaxis()->SetLabelFont(62);
		h_dNdpT_truth[i]->GetYaxis()->SetTitleOffset(1.3);

		gPad->SetLeftMargin(0.25);

		h_dNdpT_truth[i]->Draw("P");
		h_dNdpT_reco[i]->Draw("P,same");
	}
}

void plotInvariantYield()
{
	TCanvas *c1 = new TCanvas("c1", "Invariant Yield for all Sectors", 600, 600);

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_data[i]->SetLineColor(yieldColor[i]);
		h_dNdpT_data[i]->SetMarkerColor(yieldColor[i]);
		h_dNdpT_data[i]->SetMarkerStyle(20);
	}

	h_dNdpT_data[0]->Draw();
	h_dNdpT_data[1]->Draw("PE,same");
	h_dNdpT_data[2]->Draw("PE,same");
	h_dNdpT_data[3]->Draw("PE,same");
	h_dNdpT_data[4]->Draw("PE,same");

	tsallisSectorFit[0]->Draw("same");
	tsallisSectorFit[1]->Draw("same");
	tsallisSectorFit[2]->Draw("same");
	tsallisSectorFit[3]->Draw("same");
	tsallisSectorFit[4]->Draw("same");
}

void ComputeInvariantYield()
{
	gStyle->SetOptStat(0);
	readFiles();
	applyAcceptanceCorrection();
	fitSpectrum();
	divideFits();
	plotTruthReco();
	plotInvariantYield();
}