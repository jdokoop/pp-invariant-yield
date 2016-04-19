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
TH2F *h_eta_phi_data;
TH2F *h_eta_phi_reco;

//Eta distributions for odd phi regions
TH1F *h_eta_regions_data[4];
TH1F *h_eta_regions_reco[4];

//Chisd distributions
TH1F *h_chisqndf_data[NSECT];

//Label for each azimuthal sector
string sectorLabel[NSECT] = {"INCLUSIVE", "ET", "EB", "WT", "WB"};

//Label for selected phi regions (2 good, 2 bad)
string phiLabel[4] = {"2.4 < #phi < 3.1", "-0.2 < #phi < 0.2", "0.65 < #phi < 0.85", "-0.75 < #phi < -0.6"};

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
	TFile *fin = new TFile("WorkingFiles/normalizedSpectra_1.root");

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i] = (TH1F*) fin->Get(Form("h_dNdpT_truth_%s", sectorLabel[i].c_str()));
		h_dNdpT_reco[i] = (TH1F*) fin->Get(Form("h_dNdpT_reco_%s", sectorLabel[i].c_str()));
		h_dNdpT_data[i] = (TH1F*) fin->Get(Form("h_dNdpT_data_%s", sectorLabel[i].c_str()));

		h_correction[i] = (TH1F*) fin->Get(Form("h_correction_%s", sectorLabel[i].c_str()));

		h_chisqndf_data[i] = (TH1F*) fin->Get(Form("h_chisqndf_data_%s", sectorLabel[i].c_str()));
	}

	h_eta_phi_data = (TH2F*) fin->Get("h_eta_phi_data");
	h_eta_phi_reco = (TH2F*) fin->Get("h_eta_phi_sims");

	h_eta_phi_data->Scale(1.0/h_eta_phi_data->Integral());
	h_eta_phi_reco->Scale(1.0/h_eta_phi_reco->Integral());

	h_phi_data[0] = (TH1F*) h_eta_phi_data->ProjectionY("h_phi_data_inclusive");
	h_phi_reco[0] = (TH1F*) h_eta_phi_reco->ProjectionY("h_phi_reco_inclusive");

	//Create 1D eta histograms from 2D eta v. phi
	h_eta_regions_data[0] = (TH1F*) h_eta_phi_data->ProjectionX("h_eta_data_phi1",h_eta_phi_data->GetYaxis()->FindBin(2.4),h_eta_phi_data->GetYaxis()->FindBin(3.1));
	h_eta_regions_data[1] = (TH1F*) h_eta_phi_data->ProjectionX("h_eta_data_phi2",h_eta_phi_data->GetYaxis()->FindBin(-0.2),h_eta_phi_data->GetYaxis()->FindBin(0.2));
	h_eta_regions_data[2] = (TH1F*) h_eta_phi_data->ProjectionX("h_eta_data_phi3",h_eta_phi_data->GetYaxis()->FindBin(0.65),h_eta_phi_data->GetYaxis()->FindBin(0.85));
	h_eta_regions_data[3] = (TH1F*) h_eta_phi_data->ProjectionX("h_eta_data_phi4",h_eta_phi_data->GetYaxis()->FindBin(-0.75),h_eta_phi_data->GetYaxis()->FindBin(-0.6));

	h_eta_regions_reco[0] = (TH1F*) h_eta_phi_reco->ProjectionX("h_eta_reco_phi1",h_eta_phi_reco->GetYaxis()->FindBin(2.4),h_eta_phi_reco->GetYaxis()->FindBin(3.1));
	h_eta_regions_reco[1] = (TH1F*) h_eta_phi_reco->ProjectionX("h_eta_reco_phi2",h_eta_phi_reco->GetYaxis()->FindBin(-0.2),h_eta_phi_reco->GetYaxis()->FindBin(0.2));
	h_eta_regions_reco[2] = (TH1F*) h_eta_phi_reco->ProjectionX("h_eta_reco_phi3",h_eta_phi_reco->GetYaxis()->FindBin(0.65),h_eta_phi_reco->GetYaxis()->FindBin(0.85));
	h_eta_regions_reco[3] = (TH1F*) h_eta_phi_reco->ProjectionX("h_eta_reco_phi4",h_eta_phi_reco->GetYaxis()->FindBin(-0.75),h_eta_phi_reco->GetYaxis()->FindBin(-0.6));


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
	TCanvas *c1 = new TCanvas("c1", "Invariant Yield for all Sectors", 500, 900);
	c1->Divide(1, 2, 0, 0);

	c1->cd(1);
	gPad->SetPad(.005, .3, .9, .92 );
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
	//tsallisSectorFit[1]->Draw("same");
	//tsallisSectorFit[2]->Draw("same");
	//tsallisSectorFit[3]->Draw("same");
	//tsallisSectorFit[4]->Draw("same");

	c1->cd(2);
	gPad->SetPad(.005, .005, .9, .3);
	gPad->SetTickx();
	gPad->SetTicky();

	tsallisSectorRatio[0]->SetTitle();
	tsallisSectorRatio[0]->GetXaxis()->SetRangeUser(0.2, 2.0);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetXaxis()->SetLabelSize(0.04);
	tsallisSectorRatio[0]->GetXaxis()->SetTitleSize(0.04);
	tsallisSectorRatio[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	tsallisSectorRatio[0]->GetYaxis()->SetRangeUser(0.0, 2.44);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleFont(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(62);
	tsallisSectorRatio[0]->GetYaxis()->SetLabelSize(0.04);
	tsallisSectorRatio[0]->GetYaxis()->SetTitleSize(0.04);
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
	c3->SetLogy();
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
	//Vertical lines at the location of anomalous ladders
	TLine *tlcB0[3];
	tlcB0[0] = new TLine(2.21, 0, 2.21, 0.015);
	tlcB0[1] = new TLine(0.47, 0, 0.47, 0.015);
	tlcB0[2] = new TLine(0.00, 0, 0.00, 0.015);

	TBox *tbB0[3];
	tbB0[0] = new TBox(2.21-0.57/2, 0.015, 2.21+0.57/2, 0.015);
	tbB0[1] = new TBox(0.47-0.57/2, 0.015, 0.47+0.57/2, 0.015);
	tbB0[2] = new TBox(0.0-0.57/2, 0.015, 0.0+0.57/2, 0.015);

	TLine *tlcB1[8];
	tlcB1[0] = new TLine(2.78, 0, 2.78, 0.015);
	tlcB1[1] = new TLine(2.55, 0, 2.55, 0.015);
	tlcB1[2] = new TLine(3.73, 0, 3.73, 0.015);
	tlcB1[3] = new TLine(0.12, 0, 0.12, 0.015);
	tlcB1[4] = new TLine(0.59, 0, 0.59, 0.015);
	tlcB1[5] = new TLine(3.26, 0, 3.26, 0.015);
	tlcB1[6] = new TLine(3.50, 0, 3.50, 0.015);
	tlcB1[7] = new TLine(3.73, 0, 3.73, 0.015);
	
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

	h_phi_data[0]->Scale(1.0/h_phi_data[0]->Integral());
	h_phi_reco[0]->Scale(1.0/h_phi_reco[0]->Integral());

	h_phi_data[0]->SetTitle("");
	h_phi_data[0]->Draw();
	h_phi_reco[0]->SetLineColor(kRed);
	h_phi_reco[0]->Draw("same");

	tlcB0[0]->SetLineWidth(4);
	tlcB0[0]->Draw("same");

	tlcB0[1]->SetLineWidth(4);
	tlcB0[1]->Draw("same");

	tlcB0[2]->SetLineWidth(4);
	tlcB0[2]->Draw("same");

	//tbB0[0]->SetFillStyle(0);
    //tbB0[0]->SetLineColor(2);
    //tbB0[0]->SetLineWidth(2);

	//tbB0[0]->Draw("same");
	//tbB0[1]->Draw("l,same");
	//tbB0[2]->Draw("l,same");

	for(int i=0; i<8; i++)
	{
		tlcB1[i]->SetLineWidth(4);
		tlcB1[i]->SetLineColor(kGreen+3);
		tlcB1[i]->Draw("same");
	}
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

	TCanvas *cChisqAll = new TCanvas("cChisqAll", "cChisqAll", 600, 600);
	cChisqAll->SetLogy();
	h_chisqndf_data[1]->Draw();
	h_chisqndf_data[2]->Draw("same");
	h_chisqndf_data[3]->Draw("same");
	h_chisqndf_data[4]->Draw("same");

	TLine *tlChisq = new TLine(3.0, 1e-4, 3.0, 10);
	tlChisq->SetLineStyle(7);
	tlChisq->Draw("same");
}

void plotEta()
{
	TCanvas *cEta = new TCanvas("cEta", "Eta for Each Sector", 800, 800);
	cEta->Divide(2, 2);

	for (int i = 0; i < 4; i++)
	{
		cEta->cd(i + 1);

		h_eta_regions_data[i]->Rebin(2);
		h_eta_regions_reco[i]->Rebin(2);

		h_eta_regions_data[i]->SetLineWidth(2);
		h_eta_regions_data[i]->SetTitle(phiLabel[i].c_str());
		h_eta_regions_data[i]->GetXaxis()->SetLabelFont(62);
		h_eta_regions_data[i]->GetXaxis()->SetTitleFont(62);
		h_eta_regions_data[i]->GetXaxis()->SetLabelSize(62);
		h_eta_regions_data[i]->GetXaxis()->SetLabelSize(0.04);
		h_eta_regions_data[i]->GetXaxis()->SetTitleSize(0.04);
		h_eta_regions_data[i]->GetXaxis()->SetTitle("#eta");

		h_eta_regions_data[i]->GetYaxis()->SetLabelFont(62);
		h_eta_regions_data[i]->GetYaxis()->SetTitleFont(62);
		h_eta_regions_data[i]->GetYaxis()->SetLabelSize(62);
		h_eta_regions_data[i]->GetYaxis()->SetLabelSize(0.04);
		h_eta_regions_data[i]->GetYaxis()->SetTitleSize(0.04);
		h_eta_regions_data[i]->GetYaxis()->SetTitle("A.U.");

		h_eta_regions_reco[i]->SetLineColor(kBlue);
		h_eta_regions_data[i]->SetLineColor(kRed);

		h_eta_regions_data[i]->Rebin(2);
		h_eta_regions_reco[i]->Rebin(2);

		//h_eta_regions_data[i]->Scale(1.0/h_eta_regions_data[i]->Integral());
		//h_eta_regions_reco[i]->Scale(1.0/h_eta_regions_reco[i]->Integral());

		float areaData = h_eta_regions_data[i]->Integral(h_eta_regions_data[i]->FindBin(-0.4), h_eta_regions_data[i]->FindBin(0.4));
		float areaReco = h_eta_regions_reco[i]->Integral(h_eta_regions_reco[i]->FindBin(-0.4), h_eta_regions_reco[i]->FindBin(0.4));

		cout << areaReco/areaData << endl;

		h_eta_regions_data[i]->GetYaxis()->SetRangeUser(0,0.11);
		h_eta_regions_reco[i]->GetYaxis()->SetRangeUser(0,0.11);

		h_eta_regions_data[i]->GetXaxis()->SetRangeUser(-0.4, 0.4);
		h_eta_regions_reco[i]->GetXaxis()->SetRangeUser(-0.4, 0.4);

		h_eta_regions_data[i]->Draw();
		h_eta_regions_reco[i]->Draw("same");
	}
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
	plotEta();
	
}