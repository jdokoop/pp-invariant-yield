#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

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

//Label for each azimuthal sector
string sectorLabel[NSECT] = {"INCLUSIVE", "ET", "EB", "WT", "WB"};

//Colors for plotting yields
int yieldColor[NSECT] = {kBlack, kOrange + 1, kBlue, kRed, kGreen + 3};

//Index for offsetting points in x
int offsetX[NSECT] = {1, 4, 0, 2, 3};

//Tsallis fits to spectra in each sector
TF1 *tsallisSectorFit[5];

//Tsallis fit to PPG030 spectrum
TF1 *tsallisPublishedFit;

//Ratio of Tsallis fits to spectrum in each sector
TGraphErrors *tsallisSectorRatio[5];
TGraphErrors *tsallisSectorRatioPublished[5];

//Published charged hadron yield from PPG030
TGraphErrors *gHadrons;

TGraphErrors *gHadronsRatio;

//---------------------------------------
// Functions
//---------------------------------------

void readFiles()
{
	TFile *fin = new TFile("WorkingFiles/normSpectra_temp.root");

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i] = (TH1F*) fin->Get(Form("h_dNdpT_truth_%s", sectorLabel[i].c_str()));
		h_dNdpT_reco[i] = (TH1F*) fin->Get(Form("h_dNdpT_reco_%s", sectorLabel[i].c_str()));
		h_dNdpT_data[i] = (TH1F*) fin->Get(Form("h_dNdpT_data_%s", sectorLabel[i].c_str()));

		h_correction[i] = (TH1F*) fin->Get(Form("h_correction_%s", sectorLabel[i].c_str()));
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

	//tsallisFit->SetParameter(0, 2.61119e-01);
	//tsallisFit->SetParameter(1, 2.09755e+01);
	//tsallisFit->SetParameter(2, 1.68027e-01);

	tsallisFit->SetParameter(0, 0.24);
	tsallisFit->SetParameter(1, 17.62);
	tsallisFit->SetParameter(2, 0.18576);

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_data[i]->Fit("tsallisFit", "Q0R");
		tsallisSectorFit[i] = (TF1*) h_dNdpT_data[i]->GetFunction("tsallisFit");
		tsallisSectorFit[i]->SetLineColor(yieldColor[i]);
	}

	TF1 *tsallisFit2 = new TF1("tsallisFit2", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.140*([1] - 1))*([1]*[2] + 0.140)) * pow(([1]*[2] + TMath::Sqrt(0.140*0.140 + x*x))/([1]*[2]+0.140),-1*[1])", 0.3, 1.9);

	//tsallisFit2->SetParameter(0, 0.26);
	//tsallisFit2->SetParameter(1, 15.65);
	//tsallisFit2->SetParameter(2, 0.17);

	tsallisFit2->SetParameter(0, 0.25);
	tsallisFit2->SetParameter(1, 15.02);
	tsallisFit2->SetParameter(2, 0.17);

	//gHadrons->Fit("tsallisFit2", "0R");
	//tsallisPublishedFit = (TF1*) gHadrons->GetFunction("tsallisFit2");
	tsallisPublishedFit = (TF1*) tsallisFit2->Clone("tsallisPublishedFit");
}

void divideFits()
{
	//Divide yield in each sector to the inclusive fit
	double x[9];
	double y[9];
	double ex[9] = {0};
	double ey[9] = {0};

	double xpub[9];
	double ypub[9];
	double expub[9] = {0};
	double eypub[9] = {0};

	float epsilon = 0.01;

	for (int j = 0; j < NSECT; j++)
	{
		for (int i = 0; i < 9; i++)
		{
			x[i] = h_dNdpT_data[j]->GetBinCenter(i + 1) + offsetX[j] * epsilon;

			y[i] = (h_dNdpT_data[j]->GetBinContent(i + 1) / tsallisSectorFit[0]->Eval(x[i])) / 1.09;
			ey[i] = h_dNdpT_data[j]->GetBinError(i + 1) / tsallisSectorFit[0]->Eval(x[i]);

			xpub[i] = h_dNdpT_data[j]->GetBinCenter(i + 1) + offsetX[j] * epsilon;

			ypub[i] = (h_dNdpT_data[j]->GetBinContent(i + 1) / tsallisPublishedFit->Eval(xpub[i])) / 1.09;
			eypub[i] = h_dNdpT_data[j]->GetBinError(i + 1) / tsallisPublishedFit->Eval(xpub[i]);
		}

		tsallisSectorRatio[j] = new TGraphErrors(9, x, y, ex, ey);
		tsallisSectorRatio[j]->SetLineColor(yieldColor[j]);
		tsallisSectorRatio[j]->SetMarkerColor(yieldColor[j]);
		tsallisSectorRatio[j]->SetMarkerStyle(20);

		tsallisSectorRatioPublished[j] = new TGraphErrors(9, xpub, ypub, expub, eypub);
		tsallisSectorRatioPublished[j]->SetLineColor(yieldColor[j]);
		tsallisSectorRatioPublished[j]->SetMarkerColor(yieldColor[j]);
		tsallisSectorRatioPublished[j]->SetMarkerStyle(20);
	}

	//Divide published spectrum by inclusive fit
	gHadronsRatio = (TGraphErrors*) gHadrons->Clone();

	for (int i = 0; i < gHadronsRatio->GetN(); i++)
	{
		Double_t x_val;
		Double_t y_val;
		Double_t ratio;

		gHadronsRatio->GetPoint(i, x_val, y_val);

		ratio = (y_val / tsallisPublishedFit->Eval(x_val));

		gHadronsRatio->SetPoint(i, x_val, ratio);
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

		TLegend *tlegTruthRecoYield = new TLegend(0.55, 0.7, 0.85, 0.8);
		tlegTruthRecoYield->SetLineColor(kWhite);
		tlegTruthRecoYield->AddEntry(h_dNdpT_truth[i], "AMPT Truth", "P");
		tlegTruthRecoYield->AddEntry(h_dNdpT_reco[i], "AMPT Reco", "P");
		if (i == 1) tlegTruthRecoYield->Draw("same");
	}
}

void plotInvariantYield()
{
	TCanvas *c1 = new TCanvas("c1", "Invariant Yield for all Sectors", 500, 900);
	c1->Divide(1, 2, 0, 0);

	c1->cd(1);
	gPad->SetPad(.005, .3, .9, .95 );
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();
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
	g_dNdpT_data[0]->GetXaxis()->SetLabelSize(0.04);
	g_dNdpT_data[0]->GetXaxis()->SetTitleSize(0.04);
	g_dNdpT_data[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	g_dNdpT_data[0]->GetYaxis()->SetRangeUser(1e-4 + 0.0001, 10.0);
	g_dNdpT_data[0]->GetYaxis()->SetLabelFont(62);
	g_dNdpT_data[0]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[0]->GetYaxis()->SetLabelSize(62);
	g_dNdpT_data[0]->GetYaxis()->SetLabelSize(0.04);
	g_dNdpT_data[0]->GetYaxis()->SetTitleSize(0.04);
	g_dNdpT_data[0]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");

	g_dNdpT_data[0]->Draw("AP");
	g_dNdpT_data[1]->Draw("P,same");
	g_dNdpT_data[2]->Draw("P,same");
	g_dNdpT_data[3]->Draw("P,same");
	g_dNdpT_data[4]->Draw("P,same");

	tsallisSectorFit[0]->Draw("same");

	TLegend *tlegSectorYield = new TLegend(0.15, 0.1, 0.5, 0.4);
	tlegSectorYield->SetLineColor(kWhite);

	for (int i = 0; i < NSECT; i++)
	{
		tlegSectorYield->AddEntry(g_dNdpT_data[i], sectorLabel[i].c_str(), "P");
	}
	tlegSectorYield->AddEntry(tsallisSectorFit[0], "FIT TO INCLUSIVE", "L");

	tlegSectorYield->Draw("same");

	gPad->SetTopMargin(0.05);

	c1->cd(2);
	gPad->SetPad(.005, .008, .9, .3);
	gPad->SetTickx();
	gPad->SetTicky();

	tsallisSectorRatio[0]->SetTitle();
	tsallisSectorRatio[0]->GetXaxis()->SetRangeUser(0.2, 2.0);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(0.07);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleSize(0.07);
	tsallisSectorRatio[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	tsallisSectorRatio[0]->GetYaxis()->SetRangeUser(0.0, 2.44);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(0.07);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleSize(0.07);
	tsallisSectorRatio[0]->GetYaxis()->SetRangeUser(0.7, 1.4);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleOffset(0.7);
	tsallisSectorRatio[0]->GetYaxis()->SetTitle("Points / Inclusive Fit");

	tsallisSectorRatio[0]->SetLineWidth(2);
	tsallisSectorRatio[1]->SetLineWidth(2);
	tsallisSectorRatio[2]->SetLineWidth(2);
	tsallisSectorRatio[3]->SetLineWidth(2);
	tsallisSectorRatio[4]->SetLineWidth(2);

	tsallisSectorRatio[0]->Draw("AP");
	tsallisSectorRatio[1]->Draw("P,same");
	tsallisSectorRatio[2]->Draw("P,same");
	tsallisSectorRatio[3]->Draw("P,same");
	tsallisSectorRatio[4]->Draw("P,same");

	TLine *tlUnity = new TLine(0.2, 1.0, 2.0, 1.0);
	tlUnity->SetLineStyle(7);
	tlUnity->Draw("same");

	gPad->SetTopMargin(0.01);
	gPad->SetBottomMargin(0.2);

	TCanvas *c3 = new TCanvas("c3", "Invariant Yield for Inclusive Azimuth", 500, 600);
	c3->Divide(1, 2, 0, 0);

	c3->cd(1);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetLogy();
	gPad->SetTopMargin(0.05);
	gPad->SetPad(.005, .3, .9, .92 );
	g_dNdpT_data[0]->Draw("AP");
	g_dNdpT_data[0]->GetXaxis()->SetRangeUser(0, 2);
	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
	tsallisPublishedFit->Draw("same");
	tsallisPublishedFit->SetLineWidth(1);
	TLegend *tlegYield = new TLegend(0.2, 0.1, 0.55, 0.4);
	tlegYield->SetLineColor(kWhite);
	tlegYield->AddEntry(g_dNdpT_data[0],"PHI INCLUSIVE","P");
	tlegYield->AddEntry(gHadrons,"PPG030","P");
	tlegYield->AddEntry(tsallisPublishedFit,"FIT TO PPG030","L");
	tlegYield->Draw("same");

	c3->cd(2);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetBottomMargin(0.3);
	tsallisSectorRatioPublished[0]->SetTitle("");
	tsallisSectorRatioPublished[0]->Draw("AP");
	tsallisSectorRatioPublished[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	tsallisSectorRatioPublished[0]->GetXaxis()->SetRangeUser(0, 2);
	tsallisSectorRatioPublished[0]->GetXaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[0]->GetXaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[0]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[0]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetTitle("Points/ Fit to Published");
	tsallisSectorRatioPublished[0]->GetYaxis()->SetRangeUser(0.4, 1.4);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetTitleOffset(0.5);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[0]->GetYaxis()->SetTitleFont(62);
	gHadronsRatio->SetMarkerStyle(28);
	gHadronsRatio->SetMarkerColor(kMagenta);
	gHadronsRatio->Draw("P,same");
	tlUnity->Draw("same");

	TLine *tl40pct = new TLine(0.2, 0.6, 2.0, 0.6);
	tl40pct->SetLineStyle(7);
	tl40pct->Draw("same");
}

void plotInvariantYieldSectors()
{
	//-------------------
	// EAST TOP
	//-------------------
	TCanvas *cET = new TCanvas("cET", "Invariant Yield for ET", 500, 600);
	cET->Divide(1, 2, 0, 0);

	cET->cd(1);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetLogy();
	gPad->SetTopMargin(0.05);
	gPad->SetPad(.005, .3, .9, .92 );
	g_dNdpT_data[1]->SetTitle("");
	g_dNdpT_data[1]->Draw("AP");
	g_dNdpT_data[1]->GetXaxis()->SetRangeUser(0, 2);
	g_dNdpT_data[1]->GetXaxis()->SetTitleFont(62);
	g_dNdpT_data[1]->GetXaxis()->SetLabelFont(62);
	g_dNdpT_data[1]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");
	g_dNdpT_data[1]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[1]->GetYaxis()->SetLabelFont(62);
	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
	tsallisPublishedFit->Draw("same");
	tsallisPublishedFit->SetLineWidth(1);
	TLegend *tlegYield = new TLegend(0.2, 0.1, 0.55, 0.4);
	tlegYield->SetLineColor(kWhite);
	tlegYield->AddEntry(g_dNdpT_data[1],"ET","P");
	tlegYield->AddEntry(gHadrons,"PPG030","P");
	tlegYield->AddEntry(tsallisPublishedFit,"FIT TO PPG030","L");
	tlegYield->Draw("same");

	cET->cd(2);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetBottomMargin(0.3);
	tsallisSectorRatioPublished[1]->SetTitle("");
	tsallisSectorRatioPublished[1]->Draw("AP");
	tsallisSectorRatioPublished[1]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	tsallisSectorRatioPublished[1]->GetXaxis()->SetRangeUser(0, 2);
	tsallisSectorRatioPublished[1]->GetXaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[1]->GetXaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[1]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[1]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetTitle("Points/ Fit to Published");
	tsallisSectorRatioPublished[1]->GetYaxis()->SetRangeUser(0.4, 1.4);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetTitleOffset(0.5);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[1]->GetYaxis()->SetTitleFont(62);
	gHadronsRatio->SetMarkerStyle(28);
	gHadronsRatio->SetMarkerColor(kMagenta);
	gHadronsRatio->Draw("P,same");

	TLine *tlUnity = new TLine(0.2, 1.0, 2.0, 1.0);
	tlUnity->SetLineStyle(7);
	tlUnity->Draw("same");
	tlUnity->Draw("same");

	TLine *tl40pct = new TLine(0.2, 0.6, 2.0, 0.6);
	tl40pct->SetLineStyle(7);
	tl40pct->Draw("same");

	//-------------------
	// EAST BOTTOM
	//-------------------
	TCanvas *cEB = new TCanvas("cEB", "Invariant Yield for EB", 500, 600);
	cEB->Divide(1, 2, 0, 0);

	cEB->cd(1);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetLogy();
	gPad->SetTopMargin(0.05);
	gPad->SetPad(.005, .3, .9, .92 );
	g_dNdpT_data[2]->SetTitle("");
	g_dNdpT_data[2]->Draw("AP");
	g_dNdpT_data[2]->GetXaxis()->SetRangeUser(0, 2);
	g_dNdpT_data[2]->GetXaxis()->SetTitleFont(62);
	g_dNdpT_data[2]->GetXaxis()->SetLabelFont(62);
	g_dNdpT_data[2]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");
	g_dNdpT_data[2]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[2]->GetYaxis()->SetLabelFont(62);
	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
	tsallisPublishedFit->Draw("same");
	tsallisPublishedFit->SetLineWidth(1);
	TLegend *tlegYieldEB = new TLegend(0.2, 0.1, 0.55, 0.4);
	tlegYieldEB->SetLineColor(kWhite);
	tlegYieldEB->AddEntry(g_dNdpT_data[2],"EB","P");
	tlegYieldEB->AddEntry(gHadrons,"PPG030","P");
	tlegYieldEB->AddEntry(tsallisPublishedFit,"FIT TO PPG030","L");
	tlegYieldEB->Draw("same");

	cEB->cd(2);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetBottomMargin(0.3);
	tsallisSectorRatioPublished[2]->SetTitle("");
	tsallisSectorRatioPublished[2]->Draw("AP");
	tsallisSectorRatioPublished[2]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	tsallisSectorRatioPublished[2]->GetXaxis()->SetRangeUser(0, 2);
	tsallisSectorRatioPublished[2]->GetXaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[2]->GetXaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[2]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[2]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetTitle("Points/ Fit to Published");
	tsallisSectorRatioPublished[2]->GetYaxis()->SetRangeUser(0.4, 1.4);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetTitleOffset(0.5);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[2]->GetYaxis()->SetTitleFont(62);
	gHadronsRatio->SetMarkerStyle(28);
	gHadronsRatio->SetMarkerColor(kMagenta);
	gHadronsRatio->Draw("P,same");
	tlUnity->Draw("same");
	tl40pct->Draw("same");

	//-------------------
	// WEST TOP
	//-------------------
	TCanvas *cWT = new TCanvas("cWT", "Invariant Yield for WT", 500, 600);
	cWT->Divide(1, 2, 0, 0);

	cWT->cd(1);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetLogy();
	gPad->SetTopMargin(0.05);
	gPad->SetPad(.005, .3, .9, .92 );
	g_dNdpT_data[3]->SetTitle("");
	g_dNdpT_data[3]->Draw("AP");
	g_dNdpT_data[3]->GetXaxis()->SetRangeUser(0, 2);
	g_dNdpT_data[3]->GetXaxis()->SetTitleFont(62);
	g_dNdpT_data[3]->GetXaxis()->SetLabelFont(62);
	g_dNdpT_data[3]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");
	g_dNdpT_data[3]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[3]->GetYaxis()->SetLabelFont(62);	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
	tsallisPublishedFit->Draw("same");
	tsallisPublishedFit->SetLineWidth(1);
	TLegend *tlegYieldWT = new TLegend(0.2, 0.1, 0.55, 0.4);
	tlegYieldWT->SetLineColor(kWhite);
	tlegYieldWT->AddEntry(g_dNdpT_data[3],"WT","P");
	tlegYieldWT->AddEntry(gHadrons,"PPG030","P");
	tlegYieldWT->AddEntry(tsallisPublishedFit,"FIT TO PPG030","L");
	tlegYieldWT->Draw("same");

	cWT->cd(2);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetBottomMargin(0.3);
	tsallisSectorRatioPublished[3]->SetTitle("");
	tsallisSectorRatioPublished[3]->Draw("AP");
	tsallisSectorRatioPublished[3]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	tsallisSectorRatioPublished[3]->GetXaxis()->SetRangeUser(0, 2);
	tsallisSectorRatioPublished[3]->GetXaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[3]->GetXaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[3]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[3]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetTitle("Points/ Fit to Published");
	tsallisSectorRatioPublished[3]->GetYaxis()->SetRangeUser(0.4, 1.4);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetTitleOffset(0.5);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[3]->GetYaxis()->SetTitleFont(62);
	gHadronsRatio->SetMarkerStyle(28);
	gHadronsRatio->SetMarkerColor(kMagenta);
	gHadronsRatio->Draw("P,same");
	tlUnity->Draw("same");
	tl40pct->Draw("same");

	//-------------------
	// WEST BOTTOM
	//-------------------
	TCanvas *cWB = new TCanvas("cWB", "Invariant Yield for WB", 500, 600);
	cWB->Divide(1, 2, 0, 0);

	cWB->cd(1);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetLogy();
	gPad->SetTopMargin(0.05);
	gPad->SetPad(.005, .3, .9, .92 );
	g_dNdpT_data[4]->SetTitle("");
	g_dNdpT_data[4]->Draw("AP");
	g_dNdpT_data[4]->GetXaxis()->SetRangeUser(0, 2);
	g_dNdpT_data[4]->GetXaxis()->SetTitleFont(62);
	g_dNdpT_data[4]->GetXaxis()->SetLabelFont(62);
	g_dNdpT_data[4]->GetYaxis()->SetTitle("1/N_{evt} 1/2#pi p_{T} dN/dp_{T}d#eta");
	g_dNdpT_data[4]->GetYaxis()->SetTitleFont(62);
	g_dNdpT_data[4]->GetYaxis()->SetLabelFont(62);	gHadrons->SetMarkerStyle(28);
	gHadrons->SetMarkerColor(kViolet);
	gHadrons->Draw("P,same");
	tsallisPublishedFit->Draw("same");
	tsallisPublishedFit->SetLineWidth(1);
	TLegend *tlegYieldWB = new TLegend(0.2, 0.1, 0.55, 0.4);
	tlegYieldWB->SetLineColor(kWhite);
	tlegYieldWB->AddEntry(g_dNdpT_data[4],"WB","P");
	tlegYieldWB->AddEntry(gHadrons,"PPG030","P");
	tlegYieldWB->AddEntry(tsallisPublishedFit,"FIT TO PPG030","L");
	tlegYieldWB->Draw("same");

	cWB->cd(2);
	gPad->SetTicky();
	gPad->SetTickx();
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetBottomMargin(0.3);
	tsallisSectorRatioPublished[4]->SetTitle("");
	tsallisSectorRatioPublished[4]->Draw("AP");
	tsallisSectorRatioPublished[4]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	tsallisSectorRatioPublished[4]->GetXaxis()->SetRangeUser(0, 2);
	tsallisSectorRatioPublished[4]->GetXaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[4]->GetXaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[4]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[4]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetTitle("Points/ Fit to Published");
	tsallisSectorRatioPublished[4]->GetYaxis()->SetRangeUser(0.4, 1.4);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetLabelSize(0.08);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetTitleSize(0.08);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetTitleOffset(0.5);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatioPublished[4]->GetYaxis()->SetTitleFont(62);
	gHadronsRatio->SetMarkerStyle(28);
	gHadronsRatio->SetMarkerColor(kMagenta);
	gHadronsRatio->Draw("P,same");
	tlUnity->Draw("same");
	tl40pct->Draw("same");
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

	gPad->SetLeftMargin(0.15);

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
	g_correction[0]->GetYaxis()->SetTitle("Eff #times Acc");
	g_correction[0]->GetYaxis()->SetRangeUser(0, 0.5);
	g_correction[0]->GetYaxis()->SetTitleOffset(1.6);


	TLegend *tlegSectorAccEff = new TLegend(0.2, 0.5, 0.55, 0.8);
	tlegSectorAccEff->SetLineColor(kWhite);

	for (int i = 0; i < NSECT; i++)
	{
		tlegSectorAccEff->AddEntry(g_correction[i], sectorLabel[i].c_str(), "P");
	}

	g_correction[0]->Draw("AP");
	g_correction[1]->Draw("P,same");
	g_correction[2]->Draw("P,same");
	g_correction[3]->Draw("P,same");
	g_correction[4]->Draw("P,same");
	tlegSectorAccEff->Draw("same");
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
	plotInvariantYieldSectors();
	plotEfficiencyCorrection();
}