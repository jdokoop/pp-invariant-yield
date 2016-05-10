#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace std;

//-----------------------------
// Variables
//-----------------------------

//Flags
bool normIntegral = true;

//Data histograms
TH1F *h_dPhi_data;
TH1F *h_dPhi_data_lowpT;
TH1F *h_dPhi_data_highpT;

TH1F *h_dNdEta_WB_data;
TH1F *h_dNdEta_WT_data;
TH1F *h_dNdEta_EB_data;
TH1F *h_dNdEta_ET_data;

TH1F *h_dNdEta_lowpT_data;
TH1F *h_dNdEta_midpT_data;
TH1F *h_dNdEta_hipT_data;

TH1F *h_dNdEta_data;

TH1F *h_chisqndf_data;
TH1F *h_dca_data;
TH1F *h_dca2d_data;
TH1F *h_pT_data;

TH2F *h_clusters_B0_data;
TH2F *h_clusters_B1_data;
TH2F *h_clusters_B2_data;

TH1F *h_phi_clusters_B0_data;
TH1F *h_phi_clusters_B1_data;
TH1F *h_phi_clusters_B2_data;
TH1F *h_phi_clusters_B3_data;

TH1F *h_zed_clusters_B0_data;
TH1F *h_zed_clusters_B1_data;
TH1F *h_zed_clusters_B2_data;
TH1F *h_zed_clusters_B3_data;

TH2F *h_eta_phi_data;

//Simulation histograms
TH1F *h_dPhi_sims;
TH1F *h_dPhi_sims_lowpT;
TH1F *h_dPhi_sims_highpT;

TH1F *h_dNdEta_WB_sims;
TH1F *h_dNdEta_WT_sims;
TH1F *h_dNdEta_EB_sims;
TH1F *h_dNdEta_ET_sims;

TH1F *h_dNdEta_lowpT_sims;
TH1F *h_dNdEta_midpT_sims;
TH1F *h_dNdEta_hipT_sims;

TH1F *h_dNdEta_sims;

TH1F *h_chisqndf_sims;
TH1F *h_dca_sims;
TH1F *h_dca2d_sims;
TH1F *h_pT_sims;

TH2F *h_clusters_B0_sims;
TH2F *h_clusters_B1_sims;
TH2F *h_clusters_B2_sims;

TH1F *h_phi_clusters_B0_sims;
TH1F *h_phi_clusters_B1_sims;
TH1F *h_phi_clusters_B2_sims;
TH1F *h_phi_clusters_B3_sims;

TH1F *h_zed_clusters_B0_sims;
TH1F *h_zed_clusters_B1_sims;
TH1F *h_zed_clusters_B2_sims;
TH1F *h_zed_clusters_B3_sims;

TH2F *h_eta_phi_sims;

//Clock noise histograms
TH1F *h_phi_clusters_B0_clock;
TH1F *h_phi_clusters_B1_clock;
TH1F *h_phi_clusters_B2_clock;
TH1F *h_phi_clusters_B3_clock;

TH1F *h_zed_clusters_B0_clock;
TH1F *h_zed_clusters_B1_clock;
TH1F *h_zed_clusters_B2_clock;
TH1F *h_zed_clusters_B3_clock;

//Ratio histograms
TH1F *h_pT_ratio;

//AMPT histograms
TTree *ntp_svxseg_true;
TH1F *h_pT_ampt;

//Azimuthal region limits
const float WB_LO_PHI = -1 * TMath::Pi() / 2.0;
const float WB_HI_PHI = 0;

const float WT_LO_PHI = 0;
const float WT_HI_PHI = TMath::Pi() / 2.0;

const float ET_LO_PHI = TMath::Pi() / 2.0;
const float ET_HI_PHI = TMath::Pi();

const float EB_LO_PHI = TMath::Pi();
const float EB_HI_PHI = 1.5 * TMath::Pi();

//-----------------------------
// Functions
//-----------------------------

void addClockNoiseToClusters()
{
	h_phi_clusters_B0_sims->Add(h_phi_clusters_B0_clock);
	h_phi_clusters_B1_sims->Add(h_phi_clusters_B1_clock);
	h_phi_clusters_B2_sims->Add(h_phi_clusters_B2_clock);
	h_phi_clusters_B3_sims->Add(h_phi_clusters_B3_clock);

	h_zed_clusters_B0_sims->Add(h_zed_clusters_B0_clock);
	h_zed_clusters_B1_sims->Add(h_zed_clusters_B1_clock);
	h_zed_clusters_B2_sims->Add(h_zed_clusters_B2_clock);
	h_zed_clusters_B3_sims->Add(h_zed_clusters_B3_clock);
}

void normalizeEtaPhi()
{
	if (normIntegral)
	{
		//Normalize eta and phi to unity simultaneously using only the 'good' part of ET (2.4 < phi < pi)
		float normalizationSims = h_eta_phi_sims->Integral(1, h_eta_phi_sims->GetNbinsX(), h_eta_phi_sims->GetYaxis()->FindBin(2.4), h_eta_phi_sims->GetYaxis()->FindBin(TMath::Pi()));
		float normalizationData = h_eta_phi_data->Integral(1, h_eta_phi_data->GetNbinsX(), h_eta_phi_data->GetYaxis()->FindBin(2.4), h_eta_phi_data->GetYaxis()->FindBin(TMath::Pi()));

		//Normalize to unity using all of the acceptance
		//float normalizationSims = h_eta_phi_sims->Integral();
		//float normalizationData = h_eta_phi_data->Integral();

		h_eta_phi_sims->Scale(1.0 / normalizationSims);
		h_eta_phi_data->Scale(1.0 / normalizationData);
	}

	//Project phi for inclusive eta
	TH1F *hPhiSims = (TH1F*) h_eta_phi_sims->ProjectionY("hPhiSims");
	TH1F *hPhiData = (TH1F*) h_eta_phi_data->ProjectionY("hPhiData");
	hPhiData->SetLineColor(kRed);

	hPhiSims->GetYaxis()->SetTitleFont(62);
	hPhiSims->GetYaxis()->SetLabelFont(62);
	hPhiSims->GetXaxis()->SetTitleFont(62);
	hPhiSims->GetXaxis()->SetLabelFont(62);

	hPhiSims->SetTitle("");

	if (normIntegral)
	{
		hPhiSims->GetYaxis()->SetTitle("AU");
		hPhiSims->GetYaxis()->SetRangeUser(0, 0.06);
	}
	else
	{
		hPhiSims->GetYaxis()->SetTitle("Tracks / Event");
		hPhiSims->GetYaxis()->SetRangeUser(0, 0.008);
	}

	TCanvas *cPhi = new TCanvas("cPhi", "cPhi", 1200, 600);
	hPhiSims->SetLineWidth(2);
	hPhiData->SetLineWidth(2);
	hPhiSims->Draw();
	hPhiData->Draw("same");

	//Draw lines at edges between azimuthal regions
	float lineHeight = 0.008;
	if (normIntegral)
	{
		lineHeight = 0.06;
	}

	TLine *tl1_azimuth = new TLine(0, 0, 0, lineHeight);
	TLine *tl2_azimuth = new TLine(TMath::Pi() / 2, 0, TMath::Pi() / 2, lineHeight);
	TLine *tl3_azimuth = new TLine(TMath::Pi(), 0, TMath::Pi(), lineHeight);

	tl1_azimuth->SetLineWidth(5);
	tl2_azimuth->SetLineWidth(5);
	tl3_azimuth->SetLineWidth(5);

	tl1_azimuth->Draw("same");
	tl2_azimuth->Draw("same");
	tl3_azimuth->Draw("same");

	//Draw labels to indicate azimuthal regions
	TLatex *tl1_wb = new TLatex(0.15, 0.8, "WB");
	tl1_wb->SetNDC(kTRUE);
	tl1_wb->Draw("same");

	TLatex *tl1_wt = new TLatex(0.4, 0.8, "WT");
	tl1_wt->SetNDC(kTRUE);
	tl1_wt->Draw("same");

	TLatex *tl1_et = new TLatex(0.57, 0.8, "ET");
	tl1_et->SetNDC(kTRUE);
	tl1_et->Draw("same");

	TLatex *tl1_eb = new TLatex(0.75, 0.8, "EB");
	tl1_eb->SetNDC(kTRUE);
	tl1_eb->Draw("same");

	TLegend *tlegPhi = new TLegend(0.12, 0.6, 0.22, 0.7);
	tlegPhi->SetLineColor(kWhite);
	tlegPhi->AddEntry(h_dPhi_data, "DATA", "L");
	tlegPhi->AddEntry(h_dPhi_sims, "AMPT RECO", "L");
	tlegPhi->Draw("same");

	//Draw labels to indicate the ratio of data/ampt area in each sector
	float ratio_wb = (hPhiSims->Integral(hPhiSims->FindBin(WB_LO_PHI), hPhiSims->FindBin(WB_HI_PHI)))/(hPhiData->Integral(hPhiData->FindBin(WB_LO_PHI), hPhiData->FindBin(WB_HI_PHI)));
	TLatex *tl1_wb_ratio = new TLatex(0.15, 0.75, Form("%g",ratio_wb));
	tl1_wb_ratio->SetNDC(kTRUE);
	tl1_wb_ratio->SetTextSize(0.03);
	tl1_wb_ratio->SetTextColor(kGreen + 3);
	tl1_wb_ratio->Draw("same");

	float ratio_wt = (hPhiSims->Integral(hPhiSims->FindBin(WT_LO_PHI), hPhiSims->FindBin(WT_HI_PHI)))/(hPhiData->Integral(hPhiData->FindBin(WT_LO_PHI), hPhiData->FindBin(WT_HI_PHI)));
	TLatex *tl1_wt_ratio = new TLatex(0.4, 0.75, Form("%g",ratio_wt));
	tl1_wt_ratio->SetNDC(kTRUE);
	tl1_wt_ratio->SetTextSize(0.03);
	tl1_wt_ratio->SetTextColor(kGreen + 3);
	tl1_wt_ratio->Draw("same");

	float ratio_et = (hPhiSims->Integral(hPhiSims->FindBin(ET_LO_PHI), hPhiSims->FindBin(ET_HI_PHI)))/(hPhiData->Integral(hPhiData->FindBin(ET_LO_PHI), hPhiData->FindBin(ET_HI_PHI)));
	TLatex *tl1_et_ratio = new TLatex(0.57, 0.75, Form("All ET: %g",ratio_et));
	tl1_et_ratio->SetNDC(kTRUE);
	tl1_et_ratio->SetTextSize(0.03);
	tl1_et_ratio->SetTextColor(kGreen + 3);
	tl1_et_ratio->Draw("same");

	float ratio_et_good = (hPhiSims->Integral(hPhiSims->FindBin(2.4), hPhiSims->FindBin(ET_HI_PHI)))/(hPhiData->Integral(hPhiData->FindBin(2.4), hPhiData->FindBin(ET_HI_PHI)));
	TLatex *tl1_et_ratio_good = new TLatex(0.57, 0.7, Form("Good ET: %g",ratio_et_good));
	tl1_et_ratio_good->SetNDC(kTRUE);
	tl1_et_ratio_good->SetTextSize(0.03);
	tl1_et_ratio_good->SetTextColor(kGreen + 3);
	//tl1_et_ratio_good->Draw("same");

	float ratio_eb = (hPhiSims->Integral(hPhiSims->FindBin(EB_LO_PHI), hPhiSims->FindBin(EB_HI_PHI)))/(hPhiData->Integral(hPhiData->FindBin(EB_LO_PHI), hPhiData->FindBin(EB_HI_PHI)));
	TLatex *tl1_eb_ratio = new TLatex(0.75, 0.75, Form("%g",ratio_eb));
	tl1_eb_ratio->SetNDC(kTRUE);
	tl1_eb_ratio->SetTextSize(0.03);
	tl1_eb_ratio->SetTextColor(kGreen + 3);
	tl1_eb_ratio->Draw("same");

	//Project eta for each of the azimuthal regions of interest
	TH1F *hEta_ETData = (TH1F*) h_eta_phi_data->ProjectionX("hEta_ETData", h_eta_phi_data->GetYaxis()->FindBin(TMath::Pi() / 2), h_eta_phi_data->GetYaxis()->FindBin(TMath::Pi()));
	TH1F *hEta_EBData = (TH1F*) h_eta_phi_data->ProjectionX("hEta_EBData", h_eta_phi_data->GetYaxis()->FindBin(TMath::Pi()), h_eta_phi_data->GetNbinsY());
	TH1F *hEta_WTData = (TH1F*) h_eta_phi_data->ProjectionX("hEta_WTData", h_eta_phi_data->GetYaxis()->FindBin(0.0), h_eta_phi_data->GetYaxis()->FindBin(TMath::Pi() / 2));
	TH1F *hEta_WBData = (TH1F*) h_eta_phi_data->ProjectionX("hEta_WBData", 1, h_eta_phi_data->GetYaxis()->FindBin(0.0));

	TH1F *hEta_ETSims = (TH1F*) h_eta_phi_sims->ProjectionX("hEta_ETSims", h_eta_phi_sims->GetYaxis()->FindBin(TMath::Pi() / 2), h_eta_phi_sims->GetYaxis()->FindBin(TMath::Pi()));
	TH1F *hEta_EBSims = (TH1F*) h_eta_phi_sims->ProjectionX("hEta_EBSims", h_eta_phi_sims->GetYaxis()->FindBin(TMath::Pi()), h_eta_phi_sims->GetNbinsY());
	TH1F *hEta_WTSims = (TH1F*) h_eta_phi_sims->ProjectionX("hEta_WTSims", h_eta_phi_sims->GetYaxis()->FindBin(0.0), h_eta_phi_sims->GetYaxis()->FindBin(TMath::Pi() / 2));
	TH1F *hEta_WBSims = (TH1F*) h_eta_phi_sims->ProjectionX("hEta_WBSims", 1, h_eta_phi_sims->GetYaxis()->FindBin(0.0));

	hEta_ETData->Rebin(4);
	hEta_WTData->Rebin(4);
	hEta_EBData->Rebin(4);
	hEta_WBData->Rebin(4);

	hEta_ETSims->Rebin(4);
	hEta_WTSims->Rebin(4);
	hEta_EBSims->Rebin(4);
	hEta_WBSims->Rebin(4);

	hEta_ETData->SetLineWidth(2);
	hEta_WTData->SetLineWidth(2);
	hEta_EBData->SetLineWidth(2);
	hEta_WBData->SetLineWidth(2);

	hEta_ETSims->SetLineWidth(2);
	hEta_WTSims->SetLineWidth(2);
	hEta_EBSims->SetLineWidth(2);
	hEta_WBSims->SetLineWidth(2);

	TCanvas *cEta_Sectors = new TCanvas("cEta_Sectors", "cEta_Sectors", 1000, 1200);
	cEta_Sectors->Divide(2, 2);

	cEta_Sectors->cd(1);
	gPad->SetLeftMargin(0.2);
	hEta_WTData->SetLineColor(kRed);
	hEta_WTData->SetTitle("WT");
	hEta_WTData->GetXaxis()->SetLabelFont(62);
	hEta_WTData->GetXaxis()->SetTitleFont(62);
	hEta_WTData->GetXaxis()->SetTitleSize(0.05);
	hEta_WTData->GetXaxis()->SetLabelSize(0.05);
	hEta_WTData->GetYaxis()->SetLabelFont(62);
	hEta_WTData->GetYaxis()->SetTitleFont(62);
	hEta_WTData->GetYaxis()->SetTitleSize(0.05);
	hEta_WTData->GetYaxis()->SetLabelSize(0.05);
	hEta_WTData->GetYaxis()->SetTitle("Tracks / Event");
	if (normIntegral) {hEta_WTData->GetYaxis()->SetTitle("AU");}

	hEta_WTData->Draw();
	hEta_WTSims->Draw("same");
	TLegend *tlegEtaSector = new TLegend(0.65, 0.75, 0.85, 0.85);
	tlegEtaSector->SetLineColor(kWhite);
	tlegEtaSector->AddEntry(hEta_WTData, "DATA", "L");
	tlegEtaSector->AddEntry(hEta_WTSims, "AMPT RECO", "L");
	tlegEtaSector->Draw("same");

	cEta_Sectors->cd(2);
	gPad->SetLeftMargin(0.2);
	hEta_ETData->SetLineColor(kRed);
	hEta_ETData->SetTitle("ET");
	hEta_ETData->GetXaxis()->SetLabelFont(62);
	hEta_ETData->GetXaxis()->SetTitleFont(62);
	hEta_ETData->GetXaxis()->SetTitleSize(0.05);
	hEta_ETData->GetXaxis()->SetLabelSize(0.05);
	hEta_ETData->GetYaxis()->SetLabelFont(62);
	hEta_ETData->GetYaxis()->SetTitleFont(62);
	hEta_ETData->GetYaxis()->SetTitleSize(0.05);
	hEta_ETData->GetYaxis()->SetLabelSize(0.05);
	hEta_ETData->GetYaxis()->SetTitle("Tracks / Event");
	if (normIntegral) {hEta_ETData->GetYaxis()->SetTitle("AU");}

	hEta_ETData->Draw();
	hEta_ETSims->Draw("same");

	cEta_Sectors->cd(3);
	gPad->SetLeftMargin(0.2);
	hEta_WBData->SetLineColor(kRed);
	hEta_WBData->SetTitle("WB");
	hEta_WBData->GetXaxis()->SetLabelFont(62);
	hEta_WBData->GetXaxis()->SetTitleFont(62);
	hEta_WBData->GetXaxis()->SetTitleSize(0.05);
	hEta_WBData->GetXaxis()->SetLabelSize(0.05);
	hEta_WBData->GetYaxis()->SetLabelFont(62);
	hEta_WBData->GetYaxis()->SetTitleFont(62);
	hEta_WBData->GetYaxis()->SetTitleSize(0.05);
	hEta_WBData->GetYaxis()->SetLabelSize(0.05);
	hEta_WBData->GetYaxis()->SetTitle("Tracks / Event");
	if (normIntegral) {hEta_WBData->GetYaxis()->SetTitle("AU");}

	hEta_WBData->Draw();
	hEta_WBSims->Draw("same");

	cEta_Sectors->cd(4);
	gPad->SetLeftMargin(0.2);
	hEta_EBData->SetLineColor(kRed);
	hEta_EBData->SetTitle("EB");
	hEta_EBData->GetXaxis()->SetLabelFont(62);
	hEta_EBData->GetXaxis()->SetTitleFont(62);
	hEta_EBData->GetXaxis()->SetTitleSize(0.05);
	hEta_EBData->GetXaxis()->SetLabelSize(0.05);
	hEta_EBData->GetYaxis()->SetLabelFont(62);
	hEta_EBData->GetYaxis()->SetTitleFont(62);
	hEta_EBData->GetYaxis()->SetTitleSize(0.05);
	hEta_EBData->GetYaxis()->SetLabelSize(0.05);
	hEta_EBData->GetYaxis()->SetTitle("Tracks / Event");
	if (normIntegral) {hEta_EBData->GetYaxis()->SetTitle("AU");}

	hEta_EBData->Draw();
	hEta_EBSims->Draw("same");
}

void plotPhi()
{
	h_dPhi_sims_highpT->Rebin(4);
	h_dPhi_data_highpT->Rebin(4);

	//Normalize to unit integral if indicated
	if (normIntegral)
	{
		h_dPhi_sims->Scale(1.0 / h_dPhi_sims->Integral());
		h_dPhi_data->Scale(1.0 / h_dPhi_data->Integral());

		h_dPhi_sims_lowpT->Scale(1.0 / h_dPhi_sims_lowpT->Integral());
		h_dPhi_data_lowpT->Scale(1.0 / h_dPhi_data_lowpT->Integral());

		h_dPhi_sims_highpT->Scale(1.0 / h_dPhi_sims_highpT->Integral());
		h_dPhi_data_highpT->Scale(1.0 / h_dPhi_data_highpT->Integral());
	}

	TCanvas *cPhi = new TCanvas("cPhi", "cPhi", 1200, 600);
	h_dPhi_data->Draw();
	h_dPhi_sims->Draw("same");

	if (normIntegral)
	{
		h_dPhi_data->GetYaxis()->SetTitle("AU");
	}
	else
	{
		h_dPhi_data->GetYaxis()->SetTitleOffset(1.5);
		h_dPhi_data->GetYaxis()->SetTitle("Tracks/Event");
	}

	//Draw lines at edges between azimuthal regions
	TLine *tl1_azimuth = new TLine(0, 0, 0, 0.0135);
	TLine *tl2_azimuth = new TLine(TMath::Pi() / 2, 0, TMath::Pi() / 2, 0.0135);
	TLine *tl3_azimuth = new TLine(TMath::Pi(), 0, TMath::Pi(), 0.0135);

	tl1_azimuth->SetLineWidth(5);
	tl2_azimuth->SetLineWidth(5);
	tl3_azimuth->SetLineWidth(5);

	tl1_azimuth->Draw("same");
	tl2_azimuth->Draw("same");
	tl3_azimuth->Draw("same");

	//Draw labels to indicate azimuthal regions
	TLatex *tl1_wb = new TLatex(0.15, 0.8, "WB");
	tl1_wb->SetNDC(kTRUE);
	tl1_wb->Draw("same");

	TLatex *tl1_wt = new TLatex(0.4, 0.8, "WT");
	tl1_wt->SetNDC(kTRUE);
	tl1_wt->Draw("same");

	TLatex *tl1_et = new TLatex(0.57, 0.8, "ET");
	tl1_et->SetNDC(kTRUE);
	tl1_et->Draw("same");

	TLatex *tl1_eb = new TLatex(0.75, 0.8, "EB");
	tl1_eb->SetNDC(kTRUE);
	tl1_eb->Draw("same");

	TLegend *tlegPhi = new TLegend(0.12, 0.6, 0.22, 0.7);
	tlegPhi->SetLineColor(kWhite);
	tlegPhi->AddEntry(h_dPhi_data, "DATA", "L");
	tlegPhi->AddEntry(h_dPhi_sims, "AMPT RECO", "L");
	tlegPhi->Draw("same");

	//Compare ratios in selected good and bad regions
	// -->Good
	float areaRecoGood1 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(2.4), h_dPhi_sims->FindBin(3.1));
	float areaRecoGood2 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(-0.2), h_dPhi_sims->FindBin(0.2));
	float areaRecoGood3 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(0.65), h_dPhi_sims->FindBin(0.85));
	float areaRecoGood4 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(-0.75), h_dPhi_sims->FindBin(-0.6));

	float areaDataGood1 = h_dPhi_data->Integral(h_dPhi_data->FindBin(2.4), h_dPhi_data->FindBin(3.1));
	float areaDataGood2 = h_dPhi_data->Integral(h_dPhi_data->FindBin(-0.2), h_dPhi_data->FindBin(0.2));
	float areaDataGood3 = h_dPhi_data->Integral(h_dPhi_data->FindBin(0.65), h_dPhi_data->FindBin(0.85));
	float areaDataGood4 = h_dPhi_data->Integral(h_dPhi_data->FindBin(-0.75), h_dPhi_data->FindBin(-0.6));

}

void plotEta()
{
	//Normalize to unit integral to compare shapes
	if (normIntegral)
	{
		h_dNdEta_lowpT_data->Scale(1.0 / h_dNdEta_lowpT_data->Integral());
		h_dNdEta_lowpT_sims->Scale(1.0 / h_dNdEta_lowpT_sims->Integral());

		h_dNdEta_midpT_data->Scale(1.0 / h_dNdEta_midpT_data->Integral());
		h_dNdEta_midpT_sims->Scale(1.0 / h_dNdEta_midpT_sims->Integral());

		h_dNdEta_hipT_data->Scale(1.0 / h_dNdEta_hipT_data->Integral());
		h_dNdEta_hipT_sims->Scale(1.0 / h_dNdEta_hipT_sims->Integral());

		h_dNdEta_WT_data->Scale(1.0 / h_dNdEta_WT_data->Integral());
		h_dNdEta_WT_sims->Scale(1.0 / h_dNdEta_WT_sims->Integral());

		h_dNdEta_ET_data->Scale(1.0 / h_dNdEta_ET_data->Integral());
		h_dNdEta_ET_sims->Scale(1.0 / h_dNdEta_ET_sims->Integral());

		h_dNdEta_WB_data->Scale(1.0 / h_dNdEta_WB_data->Integral());
		h_dNdEta_WB_sims->Scale(1.0 / h_dNdEta_WB_sims->Integral());

		h_dNdEta_EB_data->Scale(1.0 / h_dNdEta_EB_data->Integral());
		h_dNdEta_EB_sims->Scale(1.0 / h_dNdEta_EB_sims->Integral());
	}

	TCanvas *cEta_pT = new TCanvas("cEta_pT", "cEta_pT", 1800, 400);
	cEta_pT->Divide(3, 1);

	cEta_pT->cd(1);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_lowpT_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_lowpT_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_lowpT_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_lowpT_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_lowpT_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_lowpT_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_lowpT_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_lowpT_data->GetYaxis()->SetLabelSize(0.05);
	h_dNdEta_lowpT_data->GetYaxis()->SetTitleOffset(1.85);
	if (normIntegral) {h_dNdEta_lowpT_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_lowpT_data->Draw();
	h_dNdEta_lowpT_sims->Draw("same");
	TLegend *tlegEtapT = new TLegend(0.65, 0.75, 0.85, 0.85);
	tlegEtapT->SetLineColor(kWhite);
	tlegEtapT->AddEntry(h_dNdEta_lowpT_data, "DATA", "L");
	tlegEtapT->AddEntry(h_dNdEta_lowpT_sims, "AMPT RECO", "L");
	tlegEtapT->Draw("same");


	cEta_pT->cd(2);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_midpT_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_midpT_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_midpT_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_midpT_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_midpT_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_midpT_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_midpT_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_midpT_data->GetYaxis()->SetLabelSize(0.05);
	h_dNdEta_midpT_data->GetYaxis()->SetTitleOffset(1.85);
	if (normIntegral) {h_dNdEta_midpT_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_midpT_data->Draw();
	h_dNdEta_midpT_sims->Draw("same");

	cEta_pT->cd(3);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_hipT_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_hipT_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_hipT_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_hipT_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_hipT_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_hipT_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_hipT_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_hipT_data->GetYaxis()->SetLabelSize(0.05);
	h_dNdEta_hipT_data->GetYaxis()->SetTitleOffset(1.85);
	if (normIntegral) {h_dNdEta_hipT_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_hipT_data->Draw();
	h_dNdEta_hipT_sims->Draw("same");

	TCanvas *cEta_Sectors = new TCanvas("cEta_Sectors", "cEta_Sectors", 1000, 1200);
	cEta_Sectors->Divide(2, 2);

	cEta_Sectors->cd(1);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_WT_data->SetTitle("WT");
	h_dNdEta_WT_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_WT_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_WT_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_WT_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_WT_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_WT_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_WT_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_WT_data->GetYaxis()->SetLabelSize(0.05);
	if (normIntegral) {h_dNdEta_WT_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_WT_data->Draw();
	h_dNdEta_WT_sims->Draw("same");
	TLegend *tlegEtaSector = new TLegend(0.65, 0.75, 0.85, 0.85);
	tlegEtaSector->SetLineColor(kWhite);
	tlegEtaSector->AddEntry(h_dNdEta_WT_data, "DATA", "L");
	tlegEtaSector->AddEntry(h_dNdEta_WT_sims, "AMPT RECO", "L");
	tlegEtaSector->Draw("same");

	cEta_Sectors->cd(2);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_ET_data->SetTitle("ET");
	h_dNdEta_ET_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_ET_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_ET_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_ET_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_ET_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_ET_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_ET_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_ET_data->GetYaxis()->SetLabelSize(0.05);
	if (normIntegral) {h_dNdEta_ET_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_ET_data->Draw();
	h_dNdEta_ET_sims->Draw("same");

	cEta_Sectors->cd(3);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_WB_data->SetTitle("WB");
	h_dNdEta_WB_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_WB_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_WB_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_WB_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_WB_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_WB_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_WB_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_WB_data->GetYaxis()->SetLabelSize(0.05);
	if (normIntegral) {h_dNdEta_WB_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_WB_data->Draw();
	h_dNdEta_WB_sims->Draw("same");

	cEta_Sectors->cd(4);
	gPad->SetLeftMargin(0.2);
	h_dNdEta_EB_data->SetTitle("EB");
	h_dNdEta_EB_data->GetXaxis()->SetLabelFont(62);
	h_dNdEta_EB_data->GetXaxis()->SetTitleFont(62);
	h_dNdEta_EB_data->GetXaxis()->SetTitleSize(0.05);
	h_dNdEta_EB_data->GetXaxis()->SetLabelSize(0.05);
	h_dNdEta_EB_data->GetYaxis()->SetLabelFont(62);
	h_dNdEta_EB_data->GetYaxis()->SetTitleFont(62);
	h_dNdEta_EB_data->GetYaxis()->SetTitleSize(0.05);
	h_dNdEta_EB_data->GetYaxis()->SetLabelSize(0.05);
	if (normIntegral) {h_dNdEta_EB_data->GetYaxis()->SetTitle("AU");}
	h_dNdEta_EB_data->Draw();
	h_dNdEta_EB_sims->Draw("same");
}

void plotClusters()
{
	//2D plot of spatial cluster distribution
	TCanvas *cClustersB0 = new TCanvas("cClustersB0", "cClustersB0", 1000, 500);
	cClustersB0->Divide(2, 1);

	cClustersB0->cd(1);
	h_clusters_B0_data->SetTitle("B0 - DATA");
	h_clusters_B0_data->Draw("COLZ");
	cClustersB0->cd(2);
	h_clusters_B0_sims->SetTitle("B0 - AMPT");
	h_clusters_B0_sims->Draw("COLZ");

	TCanvas *cClustersB1 = new TCanvas("cClustersB1", "cClustersB1", 1000, 500);
	cClustersB1->Divide(2, 1);

	cClustersB1->cd(1);
	h_clusters_B1_data->SetTitle("B1 - DATA");
	h_clusters_B1_data->Draw("COLZ");
	cClustersB1->cd(2);
	h_clusters_B1_sims->SetTitle("B1 - AMPT");
	h_clusters_B1_sims->Draw("COLZ");

	TCanvas *cClustersB2 = new TCanvas("cClustersB2", "cClustersB2", 1000, 500);
	cClustersB2->Divide(2, 1);

	cClustersB2->cd(1);
	h_clusters_B2_data->SetTitle("B2 - DATA");
	h_clusters_B2_data->Draw("COLZ");
	cClustersB2->cd(2);
	h_clusters_B2_sims->SetTitle("B2 - AMPT");
	h_clusters_B2_sims->Draw("COLZ");

	//Azimuthal distribution of clusters
	TCanvas *cClustersPhi = new TCanvas("cClustersPhi", "cClustersPhi", 700, 700);
	cClustersPhi->Divide(2, 2);

	cClustersPhi->cd(1);
	gPad->SetLeftMargin(0.2);
	h_phi_clusters_B0_data->GetXaxis()->SetLabelFont(62);
	h_phi_clusters_B0_data->GetXaxis()->SetTitleFont(62);
	h_phi_clusters_B0_data->GetXaxis()->SetTitleSize(0.05);
	h_phi_clusters_B0_data->GetXaxis()->SetLabelSize(0.05);
	h_phi_clusters_B0_data->GetYaxis()->SetLabelFont(62);
	h_phi_clusters_B0_data->GetYaxis()->SetTitleFont(62);
	h_phi_clusters_B0_data->GetYaxis()->SetTitleSize(0.05);
	h_phi_clusters_B0_data->GetYaxis()->SetLabelSize(0.05);
	h_phi_clusters_B0_data->GetYaxis()->SetTitleOffset(1.6);
	h_phi_clusters_B0_data->Draw();
	h_phi_clusters_B0_sims->Draw("same");
	TLegend *tlegClustersPhi = new TLegend(0.65, 0.75, 0.85, 0.85);
	tlegClustersPhi->SetLineColor(kWhite);
	tlegClustersPhi->AddEntry(h_phi_clusters_B0_data, "DATA", "L");
	tlegClustersPhi->AddEntry(h_phi_clusters_B0_sims, "AMPT RECO", "L");
	tlegClustersPhi->Draw("same");

	cClustersPhi->cd(2);
	gPad->SetLeftMargin(0.2);
	h_phi_clusters_B1_data->GetXaxis()->SetLabelFont(62);
	h_phi_clusters_B1_data->GetXaxis()->SetTitleFont(62);
	h_phi_clusters_B1_data->GetXaxis()->SetTitleSize(0.05);
	h_phi_clusters_B1_data->GetXaxis()->SetLabelSize(0.05);
	h_phi_clusters_B1_data->GetYaxis()->SetLabelFont(62);
	h_phi_clusters_B1_data->GetYaxis()->SetTitleFont(62);
	h_phi_clusters_B1_data->GetYaxis()->SetTitleSize(0.05);
	h_phi_clusters_B1_data->GetYaxis()->SetLabelSize(0.05);
	h_phi_clusters_B1_data->GetYaxis()->SetTitleOffset(1.6);
	h_phi_clusters_B1_data->Draw();
	h_phi_clusters_B1_sims->Draw("same");

	cClustersPhi->cd(3);
	gPad->SetLeftMargin(0.2);
	h_phi_clusters_B2_data->GetXaxis()->SetLabelFont(62);
	h_phi_clusters_B2_data->GetXaxis()->SetTitleFont(62);
	h_phi_clusters_B2_data->GetXaxis()->SetTitleSize(0.05);
	h_phi_clusters_B2_data->GetXaxis()->SetLabelSize(0.05);
	h_phi_clusters_B2_data->GetYaxis()->SetLabelFont(62);
	h_phi_clusters_B2_data->GetYaxis()->SetTitleFont(62);
	h_phi_clusters_B2_data->GetYaxis()->SetTitleSize(0.05);
	h_phi_clusters_B2_data->GetYaxis()->SetLabelSize(0.05);
	h_phi_clusters_B2_data->GetYaxis()->SetTitleOffset(1.6);
	h_phi_clusters_B2_data->Draw();
	h_phi_clusters_B2_sims->Draw("same");

	cout << "B2 Ratio = " << h_zed_clusters_B2_data->Integral() / h_zed_clusters_B2_sims->Integral() << endl;

	cClustersPhi->cd(4);
	gPad->SetLeftMargin(0.2);
	h_phi_clusters_B3_data->GetXaxis()->SetLabelFont(62);
	h_phi_clusters_B3_data->GetXaxis()->SetTitleFont(62);
	h_phi_clusters_B3_data->GetXaxis()->SetTitleSize(0.05);
	h_phi_clusters_B3_data->GetXaxis()->SetLabelSize(0.05);
	h_phi_clusters_B3_data->GetYaxis()->SetLabelFont(62);
	h_phi_clusters_B3_data->GetYaxis()->SetTitleFont(62);
	h_phi_clusters_B3_data->GetYaxis()->SetTitleSize(0.05);
	h_phi_clusters_B3_data->GetYaxis()->SetLabelSize(0.05);
	h_phi_clusters_B3_data->GetYaxis()->SetTitleOffset(1.6);
	h_phi_clusters_B3_data->Draw();
	h_phi_clusters_B3_sims->Draw("same");

	cout << "B3 Ratio = " << h_zed_clusters_B3_data->Integral() / h_zed_clusters_B3_sims->Integral() << endl;

	//Longitudinal distribution of clusters
	TCanvas *cClustersZed = new TCanvas("cClustersEta", "cClustersEta", 700, 700);
	cClustersZed->Divide(2, 2);

	cClustersZed->cd(1);
	gPad->SetLeftMargin(0.2);
	h_zed_clusters_B0_data->GetXaxis()->SetLabelFont(62);
	h_zed_clusters_B0_data->GetXaxis()->SetTitleFont(62);
	h_zed_clusters_B0_data->GetXaxis()->SetTitleSize(0.05);
	h_zed_clusters_B0_data->GetXaxis()->SetLabelSize(0.05);
	h_zed_clusters_B0_data->GetYaxis()->SetLabelFont(62);
	h_zed_clusters_B0_data->GetYaxis()->SetTitleFont(62);
	h_zed_clusters_B0_data->GetYaxis()->SetTitleSize(0.05);
	h_zed_clusters_B0_data->GetYaxis()->SetLabelSize(0.05);
	h_zed_clusters_B0_data->GetYaxis()->SetTitleOffset(1.6);
	h_zed_clusters_B0_data->Draw();
	h_zed_clusters_B0_sims->Draw("same");
	TLegend *tlegClustersZed = new TLegend(0.65, 0.75, 0.85, 0.85);
	tlegClustersZed->SetLineColor(kWhite);
	tlegClustersZed->AddEntry(h_zed_clusters_B0_data, "DATA", "L");
	tlegClustersZed->AddEntry(h_zed_clusters_B0_sims, "AMPT RECO", "L");
	tlegClustersZed->Draw("same");

	cClustersZed->cd(2);
	gPad->SetLeftMargin(0.2);
	h_zed_clusters_B1_data->GetXaxis()->SetLabelFont(62);
	h_zed_clusters_B1_data->GetXaxis()->SetTitleFont(62);
	h_zed_clusters_B1_data->GetXaxis()->SetTitleSize(0.05);
	h_zed_clusters_B1_data->GetXaxis()->SetLabelSize(0.05);
	h_zed_clusters_B1_data->GetYaxis()->SetLabelFont(62);
	h_zed_clusters_B1_data->GetYaxis()->SetTitleFont(62);
	h_zed_clusters_B1_data->GetYaxis()->SetTitleSize(0.05);
	h_zed_clusters_B1_data->GetYaxis()->SetLabelSize(0.05);
	h_zed_clusters_B1_data->GetYaxis()->SetTitleOffset(1.6);
	h_zed_clusters_B1_data->Draw();
	h_zed_clusters_B1_sims->Draw("same");

	cClustersZed->cd(3);
	gPad->SetLeftMargin(0.2);
	h_zed_clusters_B2_data->GetXaxis()->SetLabelFont(62);
	h_zed_clusters_B2_data->GetXaxis()->SetTitleFont(62);
	h_zed_clusters_B2_data->GetXaxis()->SetTitleSize(0.05);
	h_zed_clusters_B2_data->GetXaxis()->SetLabelSize(0.05);
	h_zed_clusters_B2_data->GetYaxis()->SetLabelFont(62);
	h_zed_clusters_B2_data->GetYaxis()->SetTitleFont(62);
	h_zed_clusters_B2_data->GetYaxis()->SetTitleSize(0.05);
	h_zed_clusters_B2_data->GetYaxis()->SetLabelSize(0.05);
	h_zed_clusters_B2_data->GetYaxis()->SetTitleOffset(1.6);
	h_zed_clusters_B2_data->Draw();
	h_zed_clusters_B2_sims->Draw("same");

	cClustersZed->cd(4);
	gPad->SetLeftMargin(0.2);
	h_zed_clusters_B3_data->GetXaxis()->SetLabelFont(62);
	h_zed_clusters_B3_data->GetXaxis()->SetTitleFont(62);
	h_zed_clusters_B3_data->GetXaxis()->SetTitleSize(0.05);
	h_zed_clusters_B3_data->GetXaxis()->SetLabelSize(0.05);
	h_zed_clusters_B3_data->GetYaxis()->SetLabelFont(62);
	h_zed_clusters_B3_data->GetYaxis()->SetTitleFont(62);
	h_zed_clusters_B3_data->GetYaxis()->SetTitleSize(0.05);
	h_zed_clusters_B3_data->GetYaxis()->SetLabelSize(0.05);
	h_zed_clusters_B3_data->GetYaxis()->SetTitleOffset(1.6);
	h_zed_clusters_B3_data->Draw();
	h_zed_clusters_B3_sims->Draw("same");
}

void plotTrackVariables()
{
	TLatex *tl_4hit = new TLatex(0.6, 0.5, "4-Hit Tracks");
	tl_4hit->SetNDC(kTRUE);
	TLatex *tl_eta = new TLatex(0.6, 0.45, "|#eta| < 1.0");
	tl_eta->SetNDC(kTRUE);
	TLatex *tl_pT = new TLatex(0.6, 0.4, "p_{T} > 0.2 GeV/c");
	tl_pT->SetNDC(kTRUE);
	TLatex *tl_dca = new TLatex(0.6, 0.35, "|DCA| < 0.15");
	tl_dca->SetNDC(kTRUE);
	TLatex *tl_dca2d = new TLatex(0.6, 0.3, "|DCA2D| < 0.05");
	tl_dca2d->SetNDC(kTRUE);

	TCanvas *cChisq = new TCanvas("cChisq", "cChisq", 600, 600);
	cChisq->SetLogy();
	h_chisqndf_data->Scale(1.0 / h_chisqndf_data->GetMaximum());
	h_chisqndf_sims->Scale(1.0 / h_chisqndf_sims->GetMaximum());
	h_chisqndf_data->GetXaxis()->SetTitle("#chi^{2}/ndf");
	h_chisqndf_data->GetYaxis()->SetTitle("AU");
	h_chisqndf_data->SetTitle("");
	h_chisqndf_data->Draw();
	h_chisqndf_sims->Draw("same");
	TLegend *tlegChisq = new TLegend(0.55, 0.7, 0.85, 0.8);
	tlegChisq->SetLineColor(kWhite);
	tlegChisq->AddEntry(h_chisqndf_data, "Data", "L");
	tlegChisq->AddEntry(h_chisqndf_sims, "AMPT Reco", "L");
	tlegChisq->Draw("same");
	TLine *lChisqCut = new TLine(3, 0, 3, 1);
	lChisqCut->SetLineStyle(7);
	lChisqCut->Draw("same");
	tl_pT->Draw("same");
	tl_4hit->Draw("same");
	tl_dca->Draw("same");
	tl_dca2d->Draw("same");
	tl_eta->Draw("same");

	TCanvas *cDCA = new TCanvas("cDCA", "cDCA", 600, 600);
	cDCA->SetLogy();
	h_dca_data->Scale(1.0 / h_dca_data->GetMaximum());
	h_dca_sims->Scale(1.0 / h_dca_sims->GetMaximum());
	h_dca_data->GetXaxis()->SetTitle("DCA");
	h_dca_data->GetYaxis()->SetTitle("AU");
	h_dca_data->GetXaxis()->SetRangeUser(-1, 1);
	h_dca_data->SetTitle("");
	h_dca_data->Draw();
	h_dca_sims->Draw("same");
	TLegend *tlegDCA = new TLegend(0.55, 0.7, 0.85, 0.8);
	tlegDCA->SetLineColor(kWhite);
	tlegDCA->AddEntry(h_dca_data, "Data", "L");
	tlegDCA->AddEntry(h_dca_sims, "AMPT Reco", "L");
	tlegDCA->Draw("same");
	TLine *lDCACut1 = new TLine(0.15, 0, 0.15, 1);
	lDCACut1->SetLineStyle(7);
	TLine *lDCACut2 = new TLine(-0.15, 0, -0.15, 1);
	lDCACut2->SetLineStyle(7);
	lDCACut1->Draw("same");
	lDCACut2->Draw("same");
	tl_eta->Draw("same");
	tl_4hit->Draw("same");

	TCanvas *cDCA2D = new TCanvas("cDCA2D", "cDCA2D", 600, 600);
	cDCA2D->SetLogy();
	h_dca2d_data->Scale(1.0 / h_dca2d_data->GetMaximum());
	h_dca2d_sims->Scale(1.0 / h_dca2d_sims->GetMaximum());
	h_dca2d_data->GetXaxis()->SetTitle("DCA2D");
	h_dca2d_data->GetYaxis()->SetTitle("AU");
	h_dca2d_data->GetXaxis()->SetRangeUser(-0.6, 0.6);
	h_dca2d_data->Draw();
	h_dca2d_data->SetTitle("");
	h_dca2d_sims->Draw("same");
	TLegend *tlegDCA2D = new TLegend(0.55, 0.7, 0.85, 0.8);
	tlegDCA2D->SetLineColor(kWhite);
	tlegDCA2D->AddEntry(h_dca2d_data, "Data", "L");
	tlegDCA2D->AddEntry(h_dca2d_sims, "AMPT Reco", "L");
	tlegDCA2D->Draw("same");
	TLine *lDCA2DCut1 = new TLine(0.05, 0, 0.05, 1);
	lDCA2DCut1->SetLineStyle(7);
	TLine *lDCA2DCut2 = new TLine(-0.05, 0, -0.05, 1);
	lDCA2DCut2->SetLineStyle(7);
	lDCA2DCut1->Draw("same");
	lDCA2DCut2->Draw("same");
	tl_eta->Draw("same");
	tl_4hit->Draw("same");
}

void CompareTrackVariables()
{
	gStyle->SetOptStat(0);

	TFile *f_data = new TFile("WorkingFiles/trackvars_data_def_temp.root");

	h_dPhi_data          = (TH1F*) f_data->Get("h_dPhi");
	h_dPhi_data_lowpT    = (TH1F*) f_data->Get("h_dPhi_lowpT");
	h_dPhi_data_highpT   = (TH1F*) f_data->Get("h_dPhi_highpT");
	h_dNdEta_lowpT_data  = (TH1F*) f_data->Get("htmp6");
	h_dNdEta_midpT_data  = (TH1F*) f_data->Get("htmp7");
	h_dNdEta_hipT_data   = (TH1F*) f_data->Get("htmp8");
	h_dNdEta_WB_data     = (TH1F*) f_data->Get("htmp4");
	h_dNdEta_WT_data     = (TH1F*) f_data->Get("htmp3");
	h_dNdEta_EB_data     = (TH1F*) f_data->Get("htmp2");
	h_dNdEta_ET_data     = (TH1F*) f_data->Get("htmp1");
	h_dNdEta_data        = (TH1F*) f_data->Get("htmp5");
	h_chisqndf_data      = (TH1F*) f_data->Get("htmpchisq");
	h_dca_data           = (TH1F*) f_data->Get("htmpdca");
	h_dca2d_data         = (TH1F*) f_data->Get("htmpdca2d");
	h_pT_data            = (TH1F*) f_data->Get("htmppt");
	h_clusters_B0_data   = (TH2F*) f_data->Get("h_clusters_B0");
	h_clusters_B1_data   = (TH2F*) f_data->Get("h_clusters_B1");
	h_clusters_B2_data   = (TH2F*) f_data->Get("h_clusters_B2");
	h_phi_clusters_B0_data    = (TH1F*) f_data->Get("h_phi_clusters_B0");
	h_phi_clusters_B1_data    = (TH1F*) f_data->Get("h_phi_clusters_B1");
	h_phi_clusters_B2_data    = (TH1F*) f_data->Get("h_phi_clusters_B2");
	h_phi_clusters_B3_data    = (TH1F*) f_data->Get("h_phi_clusters_B3");
	h_zed_clusters_B0_data    = (TH1F*) f_data->Get("h_zed_clusters_B0");
	h_zed_clusters_B1_data    = (TH1F*) f_data->Get("h_zed_clusters_B1");
	h_zed_clusters_B2_data    = (TH1F*) f_data->Get("h_zed_clusters_B2");
	h_zed_clusters_B3_data    = (TH1F*) f_data->Get("h_zed_clusters_B3");
	h_eta_phi_data            = (TH2F*) f_data->Get("h_eta_phi");

	TFile *f_sims = new TFile("WorkingFiles/trackvars_sims_def_temp.root");

	h_dPhi_sims          = (TH1F*) f_sims->Get("h_dPhi");
	h_dPhi_sims_lowpT    = (TH1F*) f_sims->Get("h_dPhi_lowpT");
	h_dPhi_sims_highpT   = (TH1F*) f_sims->Get("h_dPhi_highpT");
	h_dNdEta_WB_sims     = (TH1F*) f_sims->Get("htmp4");
	h_dNdEta_WT_sims     = (TH1F*) f_sims->Get("htmp3");
	h_dNdEta_EB_sims     = (TH1F*) f_sims->Get("htmp2");
	h_dNdEta_ET_sims     = (TH1F*) f_sims->Get("htmp1");
	h_dNdEta_lowpT_sims  = (TH1F*) f_sims->Get("htmp6");
	h_dNdEta_midpT_sims  = (TH1F*) f_sims->Get("htmp7");
	h_dNdEta_hipT_sims   = (TH1F*) f_sims->Get("htmp8");
	h_dNdEta_sims        = (TH1F*) f_sims->Get("htmp5");
	h_chisqndf_sims      = (TH1F*) f_sims->Get("htmpchisq");
	h_dca_sims           = (TH1F*) f_sims->Get("htmpdca");
	h_dca2d_sims         = (TH1F*) f_sims->Get("htmpdca2d");
	h_pT_sims            = (TH1F*) f_sims->Get("htmppt");
	h_clusters_B0_sims   = (TH2F*) f_sims->Get("h_clusters_B0");
	h_clusters_B1_sims   = (TH2F*) f_sims->Get("h_clusters_B1");
	h_clusters_B2_sims   = (TH2F*) f_sims->Get("h_clusters_B2");
	h_phi_clusters_B0_sims    = (TH1F*) f_sims->Get("h_phi_clusters_B0");
	h_phi_clusters_B1_sims    = (TH1F*) f_sims->Get("h_phi_clusters_B1");
	h_phi_clusters_B2_sims    = (TH1F*) f_sims->Get("h_phi_clusters_B2");
	h_phi_clusters_B3_sims    = (TH1F*) f_sims->Get("h_phi_clusters_B3");
	h_zed_clusters_B0_sims    = (TH1F*) f_sims->Get("h_zed_clusters_B0");
	h_zed_clusters_B1_sims    = (TH1F*) f_sims->Get("h_zed_clusters_B1");
	h_zed_clusters_B2_sims    = (TH1F*) f_sims->Get("h_zed_clusters_B2");
	h_zed_clusters_B3_sims    = (TH1F*) f_sims->Get("h_zed_clusters_B3");
	h_eta_phi_sims            = (TH2F*) f_sims->Get("h_eta_phi");

	TFile *f_clock = new TFile("WorkingFiles/trackvars_clock_def.root");
	h_phi_clusters_B0_clock    = (TH1F*) f_clock->Get("h_phi_clusters_B0");
	h_phi_clusters_B1_clock    = (TH1F*) f_clock->Get("h_phi_clusters_B1");
	h_phi_clusters_B2_clock    = (TH1F*) f_clock->Get("h_phi_clusters_B2");
	h_phi_clusters_B3_clock    = (TH1F*) f_clock->Get("h_phi_clusters_B3");
	h_zed_clusters_B0_clock    = (TH1F*) f_clock->Get("h_zed_clusters_B0");
	h_zed_clusters_B1_clock    = (TH1F*) f_clock->Get("h_zed_clusters_B1");
	h_zed_clusters_B2_clock    = (TH1F*) f_clock->Get("h_zed_clusters_B2");
	h_zed_clusters_B3_clock    = (TH1F*) f_clock->Get("h_zed_clusters_B3");

	//AMPT data
	TFile *f_ampt = new TFile("Data/ampt_pp_true_1.root");
	ntp_svxseg_true = (TTree*) f_ampt->Get("ntp_svxseg_true");
	ntp_svxseg_true->Draw("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>hampt(400,0,10)", "", "goff");
	h_pT_ampt = (TH1F*) gDirectory->FindObject("hampt");
	h_pT_ampt->Scale(1.0 / h_pT_ampt->Integral());
	h_pT_ampt->SetLineColor(kBlack);
	h_pT_ampt->SetLineWidth(2);

	//Format data histograms
	h_dPhi_data->SetLineColor(kRed);
	h_dNdEta_lowpT_data->SetLineColor(kRed);
	h_dNdEta_midpT_data->SetLineColor(kRed);
	h_dNdEta_hipT_data->SetLineColor(kRed);
	h_chisqndf_data->SetLineColor(kRed);
	h_dca_data->SetLineColor(kRed);
	h_dca2d_data->SetLineColor(kRed);
	h_pT_data->SetLineColor(kRed);
	h_phi_clusters_B0_data->SetLineColor(kRed);
	h_phi_clusters_B1_data->SetLineColor(kRed);
	h_phi_clusters_B2_data->SetLineColor(kRed);
	h_phi_clusters_B3_data->SetLineColor(kRed);
	h_zed_clusters_B0_data->SetLineColor(kRed);
	h_zed_clusters_B1_data->SetLineColor(kRed);
	h_zed_clusters_B2_data->SetLineColor(kRed);
	h_zed_clusters_B3_data->SetLineColor(kRed);
	h_dNdEta_WB_data->SetLineColor(kRed);
	h_dNdEta_WT_data->SetLineColor(kRed);
	h_dNdEta_EB_data->SetLineColor(kRed);
	h_dNdEta_ET_data->SetLineColor(kRed);
	h_dNdEta_WB_sims->SetLineColor(kBlue);
	h_dNdEta_WT_sims->SetLineColor(kBlue);
	h_dNdEta_EB_sims->SetLineColor(kBlue);
	h_dNdEta_ET_sims->SetLineColor(kBlue);

	h_dPhi_data->SetLineWidth(2);
	h_dNdEta_lowpT_data->SetLineWidth(2);
	h_dNdEta_midpT_data->SetLineWidth(2);
	h_dNdEta_hipT_data->SetLineWidth(2);
	h_chisqndf_data->SetLineWidth(2);
	h_dca_data->SetLineWidth(2);
	h_dca2d_data->SetLineWidth(2);
	h_pT_data->SetLineWidth(2);

	h_dPhi_sims->SetLineWidth(2);
	h_dNdEta_lowpT_sims->SetLineWidth(2);
	h_dNdEta_midpT_sims->SetLineWidth(2);
	h_dNdEta_hipT_sims->SetLineWidth(2);
	h_chisqndf_sims->SetLineWidth(2);
	h_dca_sims->SetLineWidth(2);
	h_dca2d_sims->SetLineWidth(2);
	h_pT_sims->SetLineWidth(2);

	normalizeEtaPhi(); //Simultaneous normalization of eta and phi distributions
	addClockNoiseToClusters();
	plotClusters();
	plotTrackVariables();
}