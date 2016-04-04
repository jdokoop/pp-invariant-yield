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
TGraphErrors *g_dNdpT_data[NSECT];

//Acceptance and Efficiency Correction
TH1F *h_correction[NSECT];
TGraphErrors *g_correction[NSECT];

//Azimuthal track distribution
TH1F *h_phi_data[NSECT];
TH1F *h_phi_reco[NSECT];

//Chisd distributions
TH1F *h_chisqndf_data[NSECT];

//Label for each azimuthal sector
string sectorLabel[NSECT] = {"INCLUSIVE", "ET", "EB", "WT", "WB"};

//Colors for plotting yields
int yieldColor[NSECT] = {kBlack, kOrange + 1, kBlue, kRed, kGreen + 3};

//Index for offsetting points in x
int offsetX[NSECT] = {1, 4, 0, 2, 3};

//Tsallis fits to spectra in each sector
TF1 *tsallisSectorFit[5];

//Ratio of Tsallis fits to spectrum in each sector
TGraphErrors *tsallisSectorRatio[5];

//Published charged hadron yield from PPG030
TGraphErrors *gHadrons;

//---------------------------------------
// Functions
//---------------------------------------

void readFiles()
{
	TFile *fin = new TFile("WorkingFiles/normalizedSpectra_phi_chisq.root");

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i] = (TH1F*) fin->Get(Form("h_dNdpT_truth_%s", sectorLabel[i].c_str()));
		h_dNdpT_reco[i] = (TH1F*) fin->Get(Form("h_dNdpT_reco_%s", sectorLabel[i].c_str()));
		h_dNdpT_data[i] = (TH1F*) fin->Get(Form("h_dNdpT_data_%s", sectorLabel[i].c_str()));

		h_correction[i] = (TH1F*) fin->Get(Form("h_correction_%s", sectorLabel[i].c_str()));

		h_phi_data[i] = (TH1F*) fin->Get(Form("h_phi_data_%s", sectorLabel[i].c_str()));
		h_phi_reco[i] = (TH1F*) fin->Get(Form("h_phi_reco_%s", sectorLabel[i].c_str()));

		h_chisqndf_data[i] = (TH1F*) fin->Get(Form("h_chisqndf_data_%s", sectorLabel[i].c_str()));
	}

	//Get hadron spectrum from PPG030
	TFile *f_ppg = new TFile("Data/PPG030_hadron_spectrum.root");
	gHadrons = (TGraphErrors*) f_ppg->Get("Graph");
	gHadrons->SetMarkerColor(kRed);
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

void applyTriggerCorrection()
{
	//From PPG030, the correction is 1/(0.8 pm 0.02) = 1.25
	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_data[i]->Scale(0.55 / 0.79);
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

	float epsilon = 0.01;

	for (int j = 0; j < NSECT; j++)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] = h_dNdpT_data[j]->GetBinCenter(i + 1) + offsetX[j] * epsilon;

			y[i] = h_dNdpT_data[j]->GetBinContent(i + 1) / tsallisSectorFit[0]->Eval(x[i]);
			ey[i] = h_dNdpT_data[j]->GetBinError(i + 1) / tsallisSectorFit[0]->Eval(x[i]);
		}

		tsallisSectorRatio[j] = new TGraphErrors(9, x, y, ex, ey);
		tsallisSectorRatio[j]->SetLineColor(yieldColor[j]);
		tsallisSectorRatio[j]->SetMarkerColor(yieldColor[j]);
		tsallisSectorRatio[j]->SetMarkerStyle(20);
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
	TCanvas *c1 = new TCanvas("c1", "Invariant Yield for all Sectors", 800, 600);
	c1->Divide(1, 2);

	c1->cd(1);
	gPad->SetLogy();

	for (int i = 0; i < NSECT; i++)
	{
		g_dNdpT_data[i]->SetLineColor(yieldColor[i]);
		g_dNdpT_data[i]->SetMarkerColor(yieldColor[i]);
		g_dNdpT_data[i]->SetMarkerStyle(20);
	}

	g_dNdpT_data[0]->SetTitle("");

	g_dNdpT_data[0]->GetXaxis()->SetRangeUser(0.2, 2.0);
	g_dNdpT_data[0]->GetXaxis()->SetLabelFont(62);
	g_dNdpT_data[0]->GetXaxis()->SetTitleFont(62);
	g_dNdpT_data[0]->GetXaxis()->SetLabelSize(62);
	g_dNdpT_data[0]->GetXaxis()->SetLabelSize(0.07);
	g_dNdpT_data[0]->GetXaxis()->SetTitleSize(0.07);
	g_dNdpT_data[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	g_dNdpT_data[0]->GetYaxis()->SetRangeUser(1e-4, 10.0);
	g_dNdpT_data[0]->GetYaxis()->SetLabelFont(62);
	g_dNdpT_data[0]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[0]->GetYaxis()->SetLabelSize(62);
	g_dNdpT_data[0]->GetYaxis()->SetLabelSize(0.07);
	g_dNdpT_data[0]->GetYaxis()->SetTitleSize(0.07);
	g_dNdpT_data[0]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");

	g_dNdpT_data[0]->Draw("AP");
	g_dNdpT_data[1]->Draw("P,same");
	g_dNdpT_data[2]->Draw("P,same");
	g_dNdpT_data[3]->Draw("P,same");
	g_dNdpT_data[4]->Draw("P,same");

	tsallisSectorFit[0]->Draw("same");
	tsallisSectorFit[1]->Draw("same");
	tsallisSectorFit[2]->Draw("same");
	tsallisSectorFit[3]->Draw("same");
	tsallisSectorFit[4]->Draw("same");

	c1->cd(2);
	tsallisSectorRatio[0]->SetTitle();
	tsallisSectorRatio[0]->GetXaxis()->SetRangeUser(0.2, 2.0);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(0.07);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleSize(0.07);
	tsallisSectorRatio[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	tsallisSectorRatio[0]->GetYaxis()->SetRangeUser(0.0, 10.0);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(0.07);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleSize(0.07);
	tsallisSectorRatio[0]->GetYaxis()->SetTitle("Points / Inclusive Fit");

	tsallisSectorRatio[0]->SetLineWidth(2);
	tsallisSectorRatio[1]->SetLineWidth(2);
	tsallisSectorRatio[2]->SetLineWidth(2);
	tsallisSectorRatio[3]->SetLineWidth(2);
	tsallisSectorRatio[4]->SetLineWidth(2);

	tsallisSectorRatio[0]->Draw("ALP");
	tsallisSectorRatio[1]->Draw("LP,same");
	tsallisSectorRatio[2]->Draw("LP,same");
	tsallisSectorRatio[3]->Draw("LP,same");
	tsallisSectorRatio[4]->Draw("LP,same");

	TCanvas *c2 = new TCanvas("c2", "Invariant Yield for all Sectors", 600, 600);
	c2->SetLogy();

	g_dNdpT_data[0]->Draw("AP");
	g_dNdpT_data[1]->Draw("P,same");
	g_dNdpT_data[2]->Draw("P,same");
	g_dNdpT_data[3]->Draw("P,same");
	g_dNdpT_data[4]->Draw("P,same");

	tsallisSectorFit[0]->Draw("same");
	tsallisSectorFit[1]->Draw("same");
	tsallisSectorFit[2]->Draw("same");
	tsallisSectorFit[3]->Draw("same");
	tsallisSectorFit[4]->Draw("same");

	gHadrons->Draw("P,same");

	TCanvas *c3 = new TCanvas("c3", "Invariant Yield for Inclusive Azimuth", 600, 600);
	g_dNdpT_data[0]->Draw("AP");
	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
}

void convertCorrectionToGraph()
{
	float x[9];
	float y[9];
	float ex[9];
	float ey[9];

	float epsilon = 0.01;

	for (int i = 0; i < NSECT; i++)
	{
		for (int j = 1; j <= 9; j++)
		{
			x[j - 1]  = h_correction[i]->GetBinCenter(j) + offsetX[i] * epsilon;
			y[j - 1]  = h_correction[i]->GetBinContent(j);
			ex[j - 1] = 0;
			ey[j - 1] = h_correction[i]->GetBinError(j);
		}

		g_correction[i] = new TGraphErrors(9, x, y, ex, ey);
	}
}

void convertYieldToGraph()
{
	float x[9];
	float y[9];
	float ex[9];
	float ey[9];

	float epsilon = 0.01;

	for (int i = 0; i < NSECT; i++)
	{
		for (int j = 1; j <= 9; j++)
		{
			x[j - 1]  = h_dNdpT_data[i]->GetBinCenter(j) + offsetX[i] * epsilon;
			y[j - 1]  = h_dNdpT_data[i]->GetBinContent(j);
			ex[j - 1] = 0;
			ey[j - 1] = h_dNdpT_data[i]->GetBinError(j);
		}

		g_dNdpT_data[i] = new TGraphErrors(9, x, y, ex, ey);
	}
}

void plotEfficiencyCorrection()
{
	TCanvas *cEff = new TCanvas("cEff", "cEff", 600, 600);

	for (int i = 0; i < NSECT; i++)
	{
		g_correction[i]->SetLineColor(yieldColor[i]);
		g_correction[i]->SetMarkerColor(yieldColor[i]);
		g_correction[i]->SetMarkerStyle(20);
	}

	g_correction[0]->SetTitle();
	g_correction[0]->GetXaxis()->SetRangeUser(0.2, 2.0);
	g_correction[0]->GetXaxis()->SetLabelFont(62);
	g_correction[0]->GetXaxis()->SetTitleFont(62);
	g_correction[0]->GetXaxis()->SetLabelSize(62);
	g_correction[0]->GetXaxis()->SetLabelSize(0.04);
	g_correction[0]->GetXaxis()->SetTitleSize(0.04);
	g_correction[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	g_correction[0]->GetYaxis()->SetLabelFont(62);
	g_correction[0]->GetYaxis()->SetTitleFont(62);
	g_correction[0]->GetYaxis()->SetLabelSize(62);
	g_correction[0]->GetYaxis()->SetLabelSize(0.04);
	g_correction[0]->GetYaxis()->SetTitleSize(0.04);
	g_correction[0]->GetYaxis()->SetTitle("#varepsilon_{acc}");

	g_correction[0]->Draw("AP");
	g_correction[1]->Draw("P,same");
	g_correction[2]->Draw("P,same");
	g_correction[3]->Draw("P,same");
	g_correction[4]->Draw("P,same");
}

void plotPhi()
{
	TCanvas *cPhiInclusive = new TCanvas("cPhiInclusive", "cPhiInclusive", 600, 600);

	h_phi_data[0]->GetXaxis()->SetLabelFont(62);
	h_phi_data[0]->GetXaxis()->SetTitleFont(62);
	h_phi_data[0]->GetXaxis()->SetLabelSize(62);
	h_phi_data[0]->GetXaxis()->SetLabelSize(0.04);
	h_phi_data[0]->GetXaxis()->SetTitleSize(0.04);
	h_phi_data[0]->GetXaxis()->SetTitle("Phi [rad]");

	h_phi_data[0]->GetYaxis()->SetRangeUser(0, 0.05);
	h_phi_data[0]->GetYaxis()->SetLabelFont(62);
	h_phi_data[0]->GetYaxis()->SetTitleFont(62);
	h_phi_data[0]->GetYaxis()->SetLabelSize(62);
	h_phi_data[0]->GetYaxis()->SetLabelSize(0.04);
	h_phi_data[0]->GetYaxis()->SetTitleSize(0.04);
	h_phi_data[0]->GetYaxis()->SetTitle("A.U.");

	h_phi_data[0]->Draw();
	h_phi_reco[0]->SetLineColor(kRed);
	h_phi_reco[0]->Draw("same");
}

void plotChisq()
{
	TCanvas *cChisq = new TCanvas("cChisq", "Chisq/ndf for Each Sector", 800, 800);
	cChisq->Divide(2, 2);

	for (int i = 1; i < NSECT; i++)
	{
		cChisq->cd(i);
		gPad->SetLogy();
		h_chisqndf_data[i]->Rebin(4);
		h_chisqndf_data[i]->SetLineColor(yieldColor[i]);
		h_chisqndf_data[i]->SetLineWidth(2);
		h_chisqndf_data[i]->SetTitle(sectorLabel[i].c_str());
		h_chisqndf_data[i]->GetXaxis()->SetLabelFont(62);
		h_chisqndf_data[i]->GetXaxis()->SetTitleFont(62);
		h_chisqndf_data[i]->GetXaxis()->SetLabelSize(62);
		h_chisqndf_data[i]->GetXaxis()->SetLabelSize(0.04);
		h_chisqndf_data[i]->GetXaxis()->SetTitleSize(0.04);
		h_chisqndf_data[i]->GetXaxis()->SetTitle("#chi^{2}/ndf");

		h_chisqndf_data[i]->GetYaxis()->SetRangeUser(1e-4, 10);
		h_chisqndf_data[i]->GetYaxis()->SetLabelFont(62);
		h_chisqndf_data[i]->GetYaxis()->SetTitleFont(62);
		h_chisqndf_data[i]->GetYaxis()->SetLabelSize(62);
		h_chisqndf_data[i]->GetYaxis()->SetLabelSize(0.04);
		h_chisqndf_data[i]->GetYaxis()->SetTitleSize(0.04);
		h_chisqndf_data[i]->GetYaxis()->SetTitle("A.U.");
		h_chisqndf_data[i]->Draw();
	}

	TCanvas *cChisqAll = new TCanvas("cChisqAll","cChisqAll",600,600);
	h_chisqndf_data[1]->Draw();
	h_chisqndf_data[2]->Draw("same");
	h_chisqndf_data[3]->Draw("same");
	h_chisqndf_data[4]->Draw("same");

}

void ComputeInvariantYield()
{
	gStyle->SetErrorX(0.0001);
	gStyle->SetOptStat(0);

	//Read data files
	readFiles();

	//Apply corrections to final invariant yield
	applyAcceptanceCorrection();
	applyTriggerCorrection();

	//Fit resulting spectra and find deviations from fit
	fitSpectrum();
	divideFits();

	//For plotting, plot TGraphs instead of histograms
	convertCorrectionToGraph();
	convertYieldToGraph();
	plotTruthReco();
	plotInvariantYield();
	plotEfficiencyCorrection();

	//Plot other track variables
	plotPhi();
	plotChisq();
}