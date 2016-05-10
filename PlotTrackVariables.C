//------------------------------------------
// Macro to plot track variables such as
// pseudorapidity, DCA(2D), pT, phi, etc.
//
// Input file: result from running the
// Run15pp200_Issues_Check module on a DST,
// either from data or PISA simulation.
//
// Will save a ROOT file with the plots,
// which can be read by PlotTrackComparison.C
// to compare data and simulations.
//------------------------------------------

#include <iostream>

using namespace std;

//------------------------------------------
// Variables
//------------------------------------------

//Flags to indicate data under consideration
// 1 = clock trigger data
// 2 = mb data
// 3 = mb sims
const int DATA = 3;

//NTuples to be read from file, with event, cluster and track information
TTree *ntp_event;
TTree *ntp_svxseg;
TTree *ntp_cluster;
TTree *ntp_svxseg_true;

//Pseudorapidity distributions in four azimuthal sectors and in 3 pT categories
TH1F *h_dNdEta_WT;
TH1F *h_dNdEta_WB;
TH1F *h_dNdEta_ET;
TH1F *h_dNdEta_EB;
TH1F *h_dNdEta;
TH1F *h_dNdEta_lowpT;  // pT < 0.5 GeV
TH1F *h_dNdEta_midpT;  // 0.5 < pT < 1.0 GeV
TH1F *h_dNdEta_highpT; // pT > 1.0 GeV

//Azimuthal track distribution
TH1F *h_dPhi;
TH1F *h_dPhi_lowpT;
TH1F *h_dPhi_highpT;

//Azimuthal and longitudinal cluster distributions in each layer
TH1F *h_phi_clusters_B0;
TH1F *h_phi_clusters_B1;
TH1F *h_phi_clusters_B2;
TH1F *h_phi_clusters_B3;

TH1F *h_zed_clusters_B0;
TH1F *h_zed_clusters_B1;
TH1F *h_zed_clusters_B2;
TH1F *h_zed_clusters_B3;

TH2F *h_clusters_B0;
TH2F *h_clusters_B1;
TH2F *h_clusters_B2;

//2D eta-phi histogram
TH2F *h_eta_phi;

//Histograms for plotting variables
TH1F *h_chisqndf;
TH1F *h_dca;
TH1F *h_dca2d;
TH1F *h_pT;

int nevents;

const float PI = TMath::Pi();

string zvtxcut  = "TMath::Abs(vtx_bbc[2]) < 1";
string chisqcut = "chisq/ndf < 3";
string dcacut   = "TMath::Abs(dca) < 0.15";
string dca2dcut = "TMath::Abs(dca2d) < 0.05";
string momcut   = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.2";
string nhitscut = "nhits[0] >= 1 && nhits[1] >= 1 && nhits[2] >= 1 && nhits[3] >= 1";
string ebcut    = "TMath::ATan2(mom[1],mom[0]) < -1*TMath::Pi()/2";
string etcut    = "TMath::ATan2(mom[1],mom[0]) > TMath::Pi()/2";
string wtcut    = "TMath::ATan2(mom[1],mom[0]) > 0 && TMath::ATan2(mom[1],mom[0]) < TMath::Pi()/2";
string wbcut    = "TMath::ATan2(mom[1],mom[0]) < 0 && TMath::ATan2(mom[1],mom[0]) > -1*TMath::Pi()/2";
string lopTcut  = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.2 && TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) < 0.5";
string midpTcut = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) < 1.0 && TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.5";
string hipTcut  = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 1.0";
string etacutWide   = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 1.0";
string etacutNarrow   = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 1.0";
string bbccut   = "pmtbbcn > 0 && pmtbbcs > 0";

//------------------------------------------
// Functions
//------------------------------------------
int determineNumEvents()
{
	if (DATA == 1)
	{
		return ntp_event->GetEntries();
	}
	else
	{
		return ntp_event->GetEntries((zvtxcut + "&&" + bbccut).c_str());
	}
}

void plotVariables()
{
	//Chisq
	ntp_svxseg->Draw("chisq/ndf>>htmpchisq(100,0,10)", (nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + momcut + "&&" + etacutNarrow + "&&" + bbccut).c_str(), "goff");
	h_chisqndf = (TH1F*) gDirectory->FindObject("htmpchisq");
	h_chisqndf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_chisqndf->Scale(1.0 / nevents);

	//DCA
	ntp_svxseg->Draw("dca>>htmpdca(400,-3,3)", (nhitscut + "&&" + etacutNarrow + "&&" + bbccut + "&&" + momcut).c_str(), "goff");
	h_dca = (TH1F*) gDirectory->FindObject("htmpdca");
	h_dca->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_dca->Scale(1.0 / nevents);

	//DCA2D
	ntp_svxseg->Draw("dca2d>>htmpdca2d(400,-3,3)", (nhitscut + "&&" + etacutNarrow + "&&" + bbccut + "&&" + momcut).c_str(), "goff");
	h_dca2d = (TH1F*) gDirectory->FindObject("htmpdca2d");
	h_dca2d->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_dca2d->Scale(1.0 / nevents);

	//pT
	ntp_svxseg->Draw("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>htmppt(400,0,10)", (nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + etacutNarrow + "&&" + bbccut).c_str(), "goff");
	h_pT = (TH1F*) gDirectory->FindObject("htmppt");
	h_pT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_pT->Scale(1.0 / nevents);

	gStyle->SetOptStat(0);

	TCanvas *c_chisq = new TCanvas("c_chisq", "c_chisq", 600, 600);
	c_chisq->SetLogy();
	h_chisqndf->SetTitle("Chisq/NDF");
	h_chisqndf->Draw();

	TCanvas *c_dca = new TCanvas("c_dca", "c_dca", 600, 600);
	c_dca->SetLogy();
	h_dca->SetTitle("DCA");
	h_dca->Draw();

	TCanvas *c_dca2d = new TCanvas("c_dca2d", "c_dca2d", 600, 600);
	c_dca2d->SetLogy();
	h_dca2d->SetTitle("DCA2D");
	h_dca2d->Draw();

	TCanvas *c_pT = new TCanvas("c_pT", "c_pT", 600, 600);
	c_pT->SetLogy();
	h_pT->SetTitle("pT");
	h_pT->Draw();
}

void getEtaPhiDistribution()
{
	h_eta_phi = new TH2F("h_eta_phi", "h_eta_phi;#eta;#phi", 200, -1.5, 1.5, 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

	float mom[3];
	float vtx_bbc[3];
	float chisq;
	float dca;
	float dca2d;
	int pmtbbcs;
	int pmtbbcn;
	int nhits[4];
	float ndf;

	ntp_svxseg->SetBranchAddress("mom", &mom);
	ntp_svxseg->SetBranchAddress("vtx_bbc", &vtx_bbc);
	ntp_svxseg->SetBranchAddress("nhits", &nhits);
	ntp_svxseg->SetBranchAddress("dca2d", &dca2d);
	ntp_svxseg->SetBranchAddress("dca", &dca);
	ntp_svxseg->SetBranchAddress("pmtbbcs", &pmtbbcs);
	ntp_svxseg->SetBranchAddress("pmtbbcs", &pmtbbcs);
	ntp_svxseg->SetBranchAddress("chisq", &chisq);
	ntp_svxseg->SetBranchAddress("ndf", &ndf);

	for (int i = 0; i < ntp_svxseg->GetEntries(); i++)
	{
		ntp_svxseg->GetEntry(i);

		float phi = TMath::ATan2(mom[1], mom[0]);
		float eta = TMath::ATanH(mom[2] / TMath::Sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]));

		if (!(nhits[0] >= 1 && nhits[1] >= 1 && nhits[2] >= 1 && nhits[3] >= 1) || TMath::Abs(vtx_bbc[2]) > 1 || TMath::Abs(dca) > 0.15 || TMath::Abs(dca2d) > 0.05 || chisq / ndf > 3 || TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) < 0.2 || TMath::Abs(TMath::ATanH(mom[2] / TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) > 1.0 || pmtbbcn <= 0 || pmtbbcs <= 0)		{
			continue;
		}

		if (phi > 1.5 * PI)
		{
			phi = phi - 2 * PI;
		}
		else if (phi < -0.5 * PI)
		{
			phi = phi + 2 * PI;
		}

		h_eta_phi->Fill(eta, phi);
	}

	h_eta_phi->Scale(1.0 / nevents);
}

void getAzimuthalDistribution()
{
	h_dPhi = new TH1F("h_dPhi", "h_dPhi", 400, -0.5 * PI, 1.5 * PI);
	h_dPhi_lowpT = new TH1F("h_dPhi_lowpT", "h_dPhi_lowpT", 400, -0.5 * PI, 1.5 * PI);
	h_dPhi_highpT = new TH1F("h_dPhi_highpT", "h_dPhi_highpT", 400, -0.5 * PI, 1.5 * PI);

	float mom[3];
	float vtx_bbc[3];
	int nhits[4];
	float chisq;
	float ndf;
	float dca;
	float dca2d;
	int pmtbbcn;
	int pmtbbcs;
	ntp_svxseg->SetBranchAddress("mom", &mom);
	ntp_svxseg->SetBranchAddress("vtx_bbc", &vtx_bbc);
	ntp_svxseg->SetBranchAddress("nhits", &nhits);
	ntp_svxseg->SetBranchAddress("chisq", &chisq);
	ntp_svxseg->SetBranchAddress("ndf", &ndf);
	ntp_svxseg->SetBranchAddress("dca", &dca);
	ntp_svxseg->SetBranchAddress("dca2d", &dca2d);
	ntp_svxseg->SetBranchAddress("pmtbbcs", &pmtbbcs);
	ntp_svxseg->SetBranchAddress("pmtbbcn", &pmtbbcn);

	for (int i = 0; i < ntp_svxseg->GetEntries(); i++)
	{
		ntp_svxseg->GetEntry(i);

		float phi = TMath::ATan2(mom[1], mom[0]);

		if (!(nhits[0] >= 1 && nhits[1] >= 1 && nhits[2] >= 1 && nhits[3] >= 1) || TMath::Abs(vtx_bbc[2]) > 1 || TMath::Abs(dca) > 0.15 || TMath::Abs(dca2d) > 0.05 || chisq / ndf > 3 || TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) < 0.2 || TMath::Abs(TMath::ATanH(mom[2] / TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) > 1.0 || pmtbbcn <= 0 || pmtbbcs <= 0)
		{
			continue;
		}

		if (phi > 1.5 * PI)
		{
			phi = phi - 2 * PI;
		}
		else if (phi < -0.5 * PI)
		{
			phi = phi + 2 * PI;
		}

		h_dPhi->Fill(phi);

		if (TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) < 0.75)
		{
			h_dPhi_lowpT->Fill(phi);
		}
		else if (TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.75)
		{
			h_dPhi_highpT->Fill(phi);
		}
	}

	TCanvas *cAzimuth = new TCanvas("cAzimuth", "cAzimuth", 600, 600);
	gStyle->SetOptStat(0);
	h_dPhi->SetTitle("");
	h_dPhi->GetXaxis()->SetTitle("Phi [rad]");
	h_dPhi->GetXaxis()->SetTitleFont(62);
	h_dPhi->GetXaxis()->SetLabelFont(62);
	h_dPhi->GetYaxis()->SetTitleFont(62);
	h_dPhi->GetYaxis()->SetLabelFont(62);
	h_dPhi->SetLineWidth(2);
	h_dPhi->Scale(1.0 / nevents);
	h_dPhi->Draw();
}

void getEtaDistribution()
{
	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp1(200,-1.5,1.5)", (zvtxcut + "&&" + etcut + "&&" + nhitscut + "&&" + dca2dcut + "&&" + dcacut  + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_ET = (TH1F*) gDirectory->FindObject("htmp1");
	h_dNdEta_ET->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp2(200,-1.5,1.5)", (zvtxcut + "&&" + ebcut + "&&" + nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_EB = (TH1F*) gDirectory->FindObject("htmp2");
	h_dNdEta_EB->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp3(200,-1.5,1.5)", (zvtxcut + "&&" + wtcut + "&&" + nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_WT = (TH1F*) gDirectory->FindObject("htmp3");
	h_dNdEta_WT->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp4(200,-1.5,1.5)", (zvtxcut + "&&" + wbcut + "&&" + nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_WB = (TH1F*) gDirectory->FindObject("htmp4");
	h_dNdEta_WB->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp5(200,-1.5,1.5)", (zvtxcut + "&&" + nhitscut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta = (TH1F*) gDirectory->FindObject("htmp5");
	h_dNdEta->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp6(200,-1.5,1.5)", (zvtxcut + "&&" + nhitscut + "&&" + lopTcut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_lowpT = (TH1F*) gDirectory->FindObject("htmp6");
	h_dNdEta_lowpT->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp7(200,-1.5,1.5)", (zvtxcut + "&&" + nhitscut + "&&" + midpTcut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_midpT = (TH1F*) gDirectory->FindObject("htmp7");
	h_dNdEta_midpT->Scale(1.0 / nevents);

	ntp_svxseg->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>htmp8(200,-1.5,1.5)", (zvtxcut + "&&" + nhitscut + "&&" + hipTcut + "&&" + dca2dcut + "&&" + dcacut + "&&" + chisqcut + "&&" + momcut + "&&" + etacutWide + "&&" + bbccut).c_str(), "goff");
	h_dNdEta_highpT = (TH1F*) gDirectory->FindObject("htmp8");
	h_dNdEta_highpT->Scale(1.0 / nevents);

	TCanvas *cEta = new TCanvas("cEta", "cEta", 600, 600);
	gStyle->SetOptStat(0);
	h_dNdEta_ET->SetTitle("Run 423844");
	h_dNdEta_ET->GetXaxis()->SetTitle("#eta");
	h_dNdEta_ET->GetXaxis()->SetTitleFont(62);
	h_dNdEta_ET->GetXaxis()->SetLabelFont(62);
	h_dNdEta_ET->GetYaxis()->SetTitle("Tracks / Event");
	h_dNdEta_ET->GetYaxis()->SetTitleOffset(1.3);
	h_dNdEta_ET->GetYaxis()->SetTitleFont(62);
	h_dNdEta_ET->GetYaxis()->SetLabelFont(62);

	h_dNdEta_ET->SetLineColor(kRed);
	h_dNdEta_EB->SetLineColor(kBlue);
	h_dNdEta_WT->SetLineColor(kBlack);
	h_dNdEta_WB->SetLineColor(kGreen + 3);

	h_dNdEta_ET->SetLineWidth(2);
	h_dNdEta_EB->SetLineWidth(2);
	h_dNdEta_WT->SetLineWidth(2);
	h_dNdEta_WB->SetLineWidth(2);

	h_dNdEta_ET->Rebin(2);
	h_dNdEta_WT->Rebin(2);
	h_dNdEta_EB->Rebin(2);
	h_dNdEta_WB->Rebin(2);

	h_dNdEta_ET->Draw();
	h_dNdEta_EB->Draw("same");
	h_dNdEta_WT->Draw("same");
	h_dNdEta_WB->Draw("same");

	TLegend *leg1 = new TLegend(0.4, 0.8, 0.5, 0.9);
	leg1->AddEntry(h_dNdEta_ET, "ET", "L");
	leg1->AddEntry(h_dNdEta_WT, "WT", "L");
	leg1->AddEntry(h_dNdEta_EB, "EB", "L");
	leg1->AddEntry(h_dNdEta_WB, "WB", "L");
	leg1->SetLineColor(kWhite);
	leg1->Draw("same");

	TCanvas *cEtaAll = new TCanvas("cEtaAll", "cEtaAll", 600, 600);
	h_dNdEta->SetTitle("Run 423844");
	h_dNdEta->GetXaxis()->SetTitle("#eta");
	h_dNdEta->GetXaxis()->SetTitleFont(62);
	h_dNdEta->GetXaxis()->SetLabelFont(62);
	h_dNdEta->GetYaxis()->SetTitle("Tracks / Event");
	h_dNdEta->GetYaxis()->SetTitleOffset(1.3);
	h_dNdEta->GetYaxis()->SetTitleFont(62);
	h_dNdEta->GetYaxis()->SetLabelFont(62);

	h_dNdEta->Rebin(2);
	h_dNdEta->SetLineWidth(2);

	h_dNdEta->Draw();

	TCanvas *cEta_pT = new TCanvas("cEta_pT", "cEta_pT", 1800, 400);
	cEta_pT->Divide(3, 1);

	cEta_pT->cd(1);
	h_dNdEta_lowpT->SetTitle("0.2 < p_{T} < 0.5 GeV");
	h_dNdEta_lowpT->GetXaxis()->SetTitle("#eta");
	h_dNdEta_lowpT->GetXaxis()->SetTitleFont(62);
	h_dNdEta_lowpT->GetXaxis()->SetLabelFont(62);
	h_dNdEta_lowpT->GetYaxis()->SetTitle("Tracks / Event");
	h_dNdEta_lowpT->GetYaxis()->SetTitleOffset(1.3);
	h_dNdEta_lowpT->GetYaxis()->SetTitleFont(62);
	h_dNdEta_lowpT->GetYaxis()->SetLabelFont(62);

	h_dNdEta_lowpT->Rebin(2);
	h_dNdEta_lowpT->SetLineWidth(2);

	h_dNdEta_lowpT->Draw();

	cEta_pT->cd(2);
	h_dNdEta_midpT->SetTitle("0.5 < p_{T} < 1.0 GeV");
	h_dNdEta_midpT->GetXaxis()->SetTitle("#eta");
	h_dNdEta_midpT->GetXaxis()->SetTitleFont(62);
	h_dNdEta_midpT->GetXaxis()->SetLabelFont(62);
	h_dNdEta_midpT->GetYaxis()->SetTitle("Tracks / Event");
	h_dNdEta_midpT->GetYaxis()->SetTitleOffset(1.3);
	h_dNdEta_midpT->GetYaxis()->SetTitleFont(62);
	h_dNdEta_midpT->GetYaxis()->SetLabelFont(62);

	h_dNdEta_midpT->Rebin(2);
	h_dNdEta_midpT->SetLineWidth(2);

	h_dNdEta_midpT->Draw();

	cEta_pT->cd(3);
	h_dNdEta_highpT->SetTitle("p_{T} > 1.0 GeV");
	h_dNdEta_highpT->GetXaxis()->SetTitle("#eta");
	h_dNdEta_highpT->GetXaxis()->SetTitleFont(62);
	h_dNdEta_highpT->GetXaxis()->SetLabelFont(62);
	h_dNdEta_highpT->GetYaxis()->SetTitle("Tracks / Event");
	h_dNdEta_highpT->GetYaxis()->SetTitleOffset(1.3);
	h_dNdEta_highpT->GetYaxis()->SetTitleFont(62);
	h_dNdEta_highpT->GetYaxis()->SetLabelFont(62);

	h_dNdEta_highpT->Rebin(2);
	h_dNdEta_highpT->SetLineWidth(2);

	h_dNdEta_highpT->Draw();
}

void plotClusters()
{
	h_clusters_B0 = new TH2F("h_clusters_B0", "h_clusters_B0", 100, -0.5 * PI, 1.5 * PI, 100, -20, 20);
	h_clusters_B1 = new TH2F("h_clusters_B1", "h_clusters_B1", 100, -0.5 * PI, 1.5 * PI, 100, -20, 20);
	h_clusters_B2 = new TH2F("h_clusters_B2", "h_clusters_B2", 100, -0.5 * PI, 1.5 * PI, 100, -20, 20);

	h_phi_clusters_B0 = new TH1F("h_phi_clusters_B0", "h_phi_clusters_B0", 100, -0.5 * PI, 1.5 * PI);
	h_phi_clusters_B1 = new TH1F("h_phi_clusters_B1", "h_phi_clusters_B1", 100, -0.5 * PI, 1.5 * PI);
	h_phi_clusters_B2 = new TH1F("h_phi_clusters_B2", "h_phi_clusters_B2", 100, -0.5 * PI, 1.5 * PI);
	h_phi_clusters_B3 = new TH1F("h_phi_clusters_B3", "h_phi_clusters_B3", 100, -0.5 * PI, 1.5 * PI);

	h_zed_clusters_B0 = new TH1F("h_zed_clusters_B0", "h_zed_clusters_B0", 100, -20, 20);
	h_zed_clusters_B1 = new TH1F("h_zed_clusters_B1", "h_zed_clusters_B1", 100, -20, 20);
	h_zed_clusters_B2 = new TH1F("h_zed_clusters_B2", "h_zed_clusters_B2", 100, -20, 20);
	h_zed_clusters_B3 = new TH1F("h_zed_clusters_B3", "h_zed_clusters_B3", 100, -20, 20);

	float cx;
	float cy;
	float cz;
	float vtx[3];
	int layer;

	ntp_cluster->SetBranchAddress("cx", &cx);
	ntp_cluster->SetBranchAddress("cy", &cy);
	ntp_cluster->SetBranchAddress("cz", &cz);
	ntp_cluster->SetBranchAddress("vtx", &vtx);
	ntp_cluster->SetBranchAddress("layer", &layer);

	for (int i = 0; i < ntp_cluster->GetEntries(); i++)
	{
		ntp_cluster->GetEntry(i);

		if (TMath::Abs(vtx[2]) > 1 && (DATA == 2 || DATA == 3))
		{
			continue;
		}

		float phi = TMath::ATan2(cy, cx);

		if (phi > 1.5 * PI)
		{
			phi = phi - 2 * PI;
		}
		else if (phi < -0.5 * PI)
		{
			phi = phi + 2 * PI;
		}


		if (layer == 0)
		{
			h_clusters_B0->Fill(phi, cz);
			h_phi_clusters_B0->Fill(phi);
			h_zed_clusters_B0->Fill(cz);
		}
		else if (layer == 1)
		{
			h_clusters_B1->Fill(phi, cz);
			h_phi_clusters_B1->Fill(phi);
			h_zed_clusters_B1->Fill(cz);
		}
		else if (layer == 2)
		{
			h_clusters_B2->Fill(phi, cz);
			h_phi_clusters_B2->Fill(phi);
			h_zed_clusters_B2->Fill(cz);
		}
		else if (layer == 3)
		{
			h_phi_clusters_B3->Fill(phi);
			h_zed_clusters_B3->Fill(cz);
		}
	}

	h_clusters_B0->Scale(1.0 / nevents);
	h_clusters_B1->Scale(1.0 / nevents);
	h_clusters_B2->Scale(1.0 / nevents);

	h_phi_clusters_B0->Scale(1.0 / nevents);
	h_phi_clusters_B1->Scale(1.0 / nevents);
	h_phi_clusters_B2->Scale(1.0 / nevents);
	h_phi_clusters_B3->Scale(1.0 / nevents);

	h_zed_clusters_B0->Scale(1.0 / nevents);
	h_zed_clusters_B1->Scale(1.0 / nevents);
	h_zed_clusters_B2->Scale(1.0 / nevents);
	h_zed_clusters_B3->Scale(1.0 / nevents);

	TCanvas *cClustersB0 = new TCanvas("cClustersB0", "cClustersB0", 600, 600);
	h_clusters_B0->SetTitle("B0 Clusters");
	h_clusters_B0->GetXaxis()->SetTitle("#phi [rad]");
	h_clusters_B0->GetYaxis()->SetTitle("z [cm]");
	h_clusters_B0->Draw("COLZ");

	TCanvas *cClustersB1 = new TCanvas("cClustersB1", "cClustersB1", 600, 600);
	h_clusters_B1->SetTitle("B1 Clusters");
	h_clusters_B1->GetXaxis()->SetTitle("#phi [rad]");
	h_clusters_B1->GetYaxis()->SetTitle("z [cm]");
	h_clusters_B1->Draw("COLZ");

	TCanvas *cClustersB2 = new TCanvas("cClustersB2", "cClustersB2", 600, 600);
	h_clusters_B2->SetTitle("B2 Clusters");
	h_clusters_B2->GetXaxis()->SetTitle("#phi [rad]");
	h_clusters_B2->GetYaxis()->SetTitle("z [cm]");
	h_clusters_B2->Draw("COLZ");

	TCanvas *cClustersPhi = new TCanvas("cClustersPhi", "cClustersPhi", 1300, 500);
	cClustersPhi->Divide(4, 1);

	cClustersPhi->cd(1);
	h_phi_clusters_B0->SetTitle("B0");
	h_phi_clusters_B0->GetXaxis()->SetTitle("#phi [rad]");
	h_phi_clusters_B0->GetYaxis()->SetTitle("Clusters / Event");
	h_phi_clusters_B0->SetLineWidth(2);
	h_phi_clusters_B0->Draw();

	cClustersPhi->cd(2);
	h_phi_clusters_B1->SetTitle("B1");
	h_phi_clusters_B1->GetXaxis()->SetTitle("#phi [rad]");
	h_phi_clusters_B1->GetYaxis()->SetTitle("Clusters / Event");
	h_phi_clusters_B1->SetLineWidth(2);
	h_phi_clusters_B1->Draw();

	cClustersPhi->cd(3);
	h_phi_clusters_B2->SetTitle("B2");
	h_phi_clusters_B2->GetXaxis()->SetTitle("#phi [rad]");
	h_phi_clusters_B2->GetYaxis()->SetTitle("Clusters / Event");
	h_phi_clusters_B2->SetLineWidth(2);
	h_phi_clusters_B2->Draw();

	cClustersPhi->cd(4);
	h_phi_clusters_B3->SetTitle("B3");
	h_phi_clusters_B3->GetXaxis()->SetTitle("#phi [rad]");
	h_phi_clusters_B3->GetYaxis()->SetTitle("Clusters / Event");
	h_phi_clusters_B3->SetLineWidth(2);
	h_phi_clusters_B3->Draw();

	TCanvas *cClustersZed = new TCanvas("cClustersEta", "cClustersEta", 1300, 500);
	cClustersZed->Divide(4, 1);

	cClustersZed->cd(1);
	h_zed_clusters_B0->SetTitle("B0");
	h_zed_clusters_B0->GetXaxis()->SetTitle("z[cm]");
	h_zed_clusters_B0->GetYaxis()->SetTitle("Clusters / Event");
	h_zed_clusters_B0->SetLineWidth(2);
	h_zed_clusters_B0->Draw();

	cClustersZed->cd(2);
	h_zed_clusters_B1->SetTitle("B1");
	h_zed_clusters_B1->GetXaxis()->SetTitle("z[cm]");
	h_zed_clusters_B1->GetYaxis()->SetTitle("Clusters / Event");
	h_zed_clusters_B1->SetLineWidth(2);
	h_zed_clusters_B1->Draw();

	cClustersZed->cd(3);
	h_zed_clusters_B2->SetTitle("B2");
	h_zed_clusters_B2->GetXaxis()->SetTitle("z[cm]");
	h_zed_clusters_B2->GetYaxis()->SetTitle("Clusters / Event");
	h_zed_clusters_B2->SetLineWidth(2);
	h_zed_clusters_B2->Draw();

	cClustersZed->cd(4);
	h_zed_clusters_B3->SetTitle("B3");
	h_zed_clusters_B3->GetXaxis()->SetTitle("z[cm]");
	h_zed_clusters_B3->GetYaxis()->SetTitle("Clusters / Event");
	h_zed_clusters_B3->SetLineWidth(2);
	h_zed_clusters_B3->Draw();
}

void PlotTrackVariables()
{
	//Read file
	TFile *fin;
	if (DATA == 1)
	{
		fin  = new TFile("Data/423844_reco_clock.root");
	}
	else if (DATA == 2)
	{
		fin  = new TFile("Data/423844_data_0_1_050216.root");
	}
	else if (DATA == 3)
	{
		fin  = new TFile("Data/423844_reco_1_0_1_0_050416.root");
	}

	ntp_svxseg  = (TTree*) fin->Get("ntp_svxseg");
	ntp_event   = (TTree*) fin->Get("ntp_event");
	ntp_cluster = (TTree*) fin->Get("ntp_cluster");

	//Get number of events
	nevents = determineNumEvents();

	plotVariables();
	getAzimuthalDistribution();
	getEtaDistribution();
	getEtaPhiDistribution();
	plotClusters();

	//Write histograms to file
	TFile *fout;
	if (DATA == 1)
	{
		fout = new TFile("WorkingFiles/trackvars_clock_def_temp.root", "RECREATE");
	}
	else if (DATA == 2)
	{
		fout = new TFile("WorkingFiles/trackvars_data_def_temp.root", "RECREATE");
	}
	else if (DATA == 3)
	{
		fout = new TFile("WorkingFiles/trackvars_sims_def_temp.root", "RECREATE");
	}

	h_dNdEta_WT->Write();
	h_dNdEta_WB->Write();
	h_dNdEta_ET->Write();
	h_dNdEta_EB->Write();
	h_dNdEta->Write();
	h_dNdEta_lowpT->Write();
	h_dNdEta_midpT->Write();
	h_dNdEta_highpT->Write();
	h_dPhi->Write();
	h_dPhi_lowpT->Write();
	h_dPhi_highpT->Write();
	h_chisqndf->Write();
	h_dca->Write();
	h_dca2d->Write();
	h_pT->Write();
	h_clusters_B0->Write();
	h_clusters_B1->Write();
	h_clusters_B2->Write();
	h_phi_clusters_B0->Write();
	h_phi_clusters_B1->Write();
	h_phi_clusters_B2->Write();
	h_phi_clusters_B3->Write();
	h_zed_clusters_B0->Write();
	h_zed_clusters_B1->Write();
	h_zed_clusters_B2->Write();
	h_zed_clusters_B3->Write();
	h_eta_phi->Write();
}
