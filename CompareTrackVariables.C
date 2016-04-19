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

//Ratio histograms
TH1F *h_pT_ratio;

//AMPT histograms
TTree *ntp_svxseg_true;
TH1F *h_pT_ampt;

//-----------------------------
// Functions
//-----------------------------
void plotPhi()
{
	h_dPhi_sims_highpT->Rebin(4);
	h_dPhi_data_highpT->Rebin(4);

	//Normalize to unit integral if indicated
	if (normIntegral)
	{
		cout << "PHI SIMS NORM = " << h_dPhi_sims->Integral() << endl;
		cout << "PHI DATA NORM = " << h_dPhi_data->Integral() << endl << endl;

		h_dPhi_sims->Scale(1.0 / h_dPhi_sims->Integral());
		h_dPhi_data->Scale(1.0 / h_dPhi_data->Integral());

		h_dPhi_sims_lowpT->Scale(1.0 / h_dPhi_sims_lowpT->Integral());
		h_dPhi_data_lowpT->Scale(1.0 / h_dPhi_data_lowpT->Integral());

		h_dPhi_sims_highpT->Scale(1.0 / h_dPhi_sims_highpT->Integral());
		h_dPhi_data_highpT->Scale(1.0 / h_dPhi_data_highpT->Integral());
	}

	TCanvas *cPhi = new TCanvas("cPhi", "cPhi", 600, 600);
	h_dPhi_data->Draw();
	h_dPhi_sims->Draw("same");

	TCanvas *cPhiMom = new TCanvas("cPhiMom","cPhiMom",600,800);
	cPhiMom->Divide(1,2);

	cPhiMom->cd(1);
	h_dPhi_sims->SetTitle("p_{T} > 0.2");
	h_dPhi_sims->GetXaxis()->SetTitle("Phi [rad]");
	h_dPhi_sims->SetLineColor(kBlue);
	h_dPhi_data->SetLineColor(kRed);

	h_dPhi_sims->Draw();
	h_dPhi_data->Draw("same");

	cPhiMom->cd(2);
	h_dPhi_sims_highpT->SetTitle("p_{T} > 0.75 GeV/c");
	h_dPhi_sims_highpT->GetXaxis()->SetTitle("Phi [rad]");
	h_dPhi_sims_highpT->SetLineColor(kBlue);
	h_dPhi_data_highpT->SetLineColor(kRed);

	h_dPhi_sims_highpT->Draw();
	h_dPhi_data_highpT->Draw("same");	

	//Compare ratios in selected good and bad regions
	// -->Good
	float areaRecoGood1 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(2.4),h_dPhi_sims->FindBin(3.1));
	float areaRecoGood2 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(-0.2),h_dPhi_sims->FindBin(0.2));
	float areaRecoGood3 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(0.65),h_dPhi_sims->FindBin(0.85));
	float areaRecoGood4 = h_dPhi_sims->Integral(h_dPhi_sims->FindBin(-0.75),h_dPhi_sims->FindBin(-0.6));

	float areaDataGood1 = h_dPhi_data->Integral(h_dPhi_data->FindBin(2.4),h_dPhi_data->FindBin(3.1));
	float areaDataGood2 = h_dPhi_data->Integral(h_dPhi_data->FindBin(-0.2),h_dPhi_data->FindBin(0.2));
	float areaDataGood3 = h_dPhi_data->Integral(h_dPhi_data->FindBin(0.65),h_dPhi_data->FindBin(0.85));
	float areaDataGood4 = h_dPhi_data->Integral(h_dPhi_data->FindBin(-0.75),h_dPhi_data->FindBin(-0.6));

	cout << "*** " << areaRecoGood1/areaDataGood1 << endl;
	cout << "*** " << areaRecoGood2/areaDataGood2 << endl;
	cout << "*** " << areaRecoGood3/areaDataGood3 << endl;
	cout << "*** " << areaRecoGood4/areaDataGood4 << endl;

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
	h_dNdEta_lowpT_data->Draw();
	h_dNdEta_lowpT_sims->Draw("same");

	cout << "-----------------------------" << endl;
	cout << "AMPT/DATA Eta Low pT Ratio = " << (float) h_dNdEta_lowpT_sims->Integral() / h_dNdEta_lowpT_data->Integral() << endl;
	cout << "-----------------------------" << endl;

	cEta_pT->cd(2);
	h_dNdEta_midpT_data->Draw();
	h_dNdEta_midpT_sims->Draw("same");

	cout << "-----------------------------" << endl;
	cout << "AMPT/DATA Eta Mid pT Ratio = " << (float) h_dNdEta_midpT_sims->Integral() / h_dNdEta_midpT_data->Integral() << endl;
	cout << "-----------------------------" << endl;

	cEta_pT->cd(3);
	h_dNdEta_hipT_data->Draw();
	h_dNdEta_hipT_sims->Draw("same");

	cout << "-----------------------------" << endl;
	cout << "AMPT/DATA Eta High pT Ratio = " << (float) h_dNdEta_hipT_sims->Integral() / h_dNdEta_hipT_data->Integral() << endl;
	cout << "-----------------------------" << endl;

	TCanvas *cEta_Sectors = new TCanvas("cEta_Sectors", "cEta_Sectors", 1000, 1200);
	cEta_Sectors->Divide(2, 2);

	cEta_Sectors->cd(1);
	h_dNdEta_WT_data->SetTitle("WT");
	h_dNdEta_WT_data->Draw();
	h_dNdEta_WT_sims->Draw("same");

	cEta_Sectors->cd(2);
	h_dNdEta_ET_data->SetTitle("ET");
	h_dNdEta_ET_data->Draw();
	h_dNdEta_ET_sims->Draw("same");

	cEta_Sectors->cd(3);
	h_dNdEta_WB_data->SetTitle("WB");
	h_dNdEta_WB_data->Draw();
	h_dNdEta_WB_sims->Draw("same");

	cEta_Sectors->cd(4);
	h_dNdEta_EB_data->SetTitle("EB");
	h_dNdEta_EB_data->Draw();
	h_dNdEta_EB_sims->Draw("same");
}

void plotClusters()
{
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

	TCanvas *cClustersPhi = new TCanvas("cClustersPhi", "cClustersPhi", 1300, 300);
	cClustersPhi->Divide(4, 1);

	cClustersPhi->cd(1);
	h_phi_clusters_B0_data->Draw();
	h_phi_clusters_B0_sims->Draw("same");

	cClustersPhi->cd(2);
	h_phi_clusters_B1_data->Draw();
	h_phi_clusters_B1_sims->Draw("same");

	cClustersPhi->cd(3);
	h_phi_clusters_B2_data->Draw();
	h_phi_clusters_B2_sims->Draw("same");

	cClustersPhi->cd(4);
	h_phi_clusters_B3_data->Draw();
	h_phi_clusters_B3_sims->Draw("same");

	TCanvas *cClustersZed = new TCanvas("cClustersEta", "cClustersEta", 1300, 300);
	cClustersZed->Divide(4, 1);

	cClustersZed->cd(1);
	h_zed_clusters_B0_data->Draw();
	h_zed_clusters_B0_sims->Draw("same");

	cClustersZed->cd(2);
	h_zed_clusters_B1_data->Draw();
	h_zed_clusters_B1_sims->Draw("same");

	cClustersZed->cd(3);
	h_zed_clusters_B2_data->Draw();
	h_zed_clusters_B2_sims->Draw("same");

	cClustersZed->cd(4);
	h_zed_clusters_B3_data->Draw();
	h_zed_clusters_B3_sims->Draw("same");
}

void plotTrackVariables()
{
	TCanvas *cChisq = new TCanvas("cChisq", "cChisq", 600, 600);
	h_chisqndf_data->Scale(1.0 / h_chisqndf_data->GetMaximum());
	h_chisqndf_sims->Scale(1.0 / h_chisqndf_sims->GetMaximum());
	h_chisqndf_data->Draw();
	h_chisqndf_sims->Draw("same");

	TCanvas *cDCA = new TCanvas("cDCA", "cDCA", 600, 600);
	h_dca_data->Scale(1.0 / h_dca_data->GetMaximum());
	h_dca_sims->Scale(1.0 / h_dca_sims->GetMaximum());
	h_dca_data->Draw();
	h_dca_sims->Draw("same");

	TCanvas *cDCA2D = new TCanvas("cDCA2D", "cDCA2D", 600, 600);
	h_dca2d_data->Scale(1.0 / h_dca2d_data->GetMaximum());
	h_dca2d_sims->Scale(1.0 / h_dca2d_sims->GetMaximum());
	h_dca2d_data->Draw();
	h_dca2d_sims->Draw("same");

	TCanvas *cpT = new TCanvas("cpT", "cpT", 600, 600);
	cpT->Divide(1, 2);

	cpT->cd(1);
	h_pT_data->Scale(1.0 / h_pT_data->GetMaximum());
	h_pT_sims->Scale(1.0 / h_pT_sims->GetMaximum());
	h_pT_ampt->Scale(1.0 / h_pT_ampt->GetMaximum());

	h_pT_data->SetTitle("");
	h_pT_data->Draw();
	h_pT_sims->Draw("same");
	h_pT_ampt->Draw("same");

	//Take ratio between data and sims pT
	h_pT_ratio = (TH1F*) h_pT_sims->Clone("h_pT_ratio");
	h_pT_ratio->Divide(h_pT_data);
	cpT->cd(2);
	h_pT_ratio->SetTitle("");
	h_pT_ratio->SetLineColor(kOrange - 3);
	h_pT_ratio->Draw("same");
}

void CompareTrackVariables()
{
	gStyle->SetOptStat(0);

	TFile *f_data = new TFile("WorkingFiles/dataTrackVariables_1_1_1.root");

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

	TFile *f_sims = new TFile("WorkingFiles/simsTrackVariables_newmat_1.root");

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

	//AMPT data
	TFile *f_ampt = new TFile("Data/ampt_pp_true_phifixed.root");
	ntp_svxseg_true = (TTree*) f_ampt->Get("ntp_svxseg_true");
	ntp_svxseg_true->Draw("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>hampt(400,0,10)");
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

	//Plot things
	plotPhi();
	plotEta();
	plotClusters();
	plotTrackVariables();
}