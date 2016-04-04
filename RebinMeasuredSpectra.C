//------------------------------------------------------------------------------
// This macro takes the data, PISA output and true AMPT SvxSegment information.
// It then produces the segment pT distribution with quality cuts applied.
// The distributions are normalized to compute the invariant yield, and the
// true and reconstructed AMPT information is used to determine the acceptance
// and efficiency correction to apply to data
//------------------------------------------------------------------------------

#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

//Number of original pT bins
const int NBINS = 9;

//Number of bins in rebinned spectra
const int NBINS_RB = 12;

//Number of sectors
// --> Inclusive phi, ET, EB, WT, WB
const int NSECT = 5;

//Number of events in each file
int nevents_truth = 447000;
int nevents_reco  = 0;
int nevents_data  = 0;

//Raw pT yield
TH1F *h_DeltaNDeltapT_truth[NSECT];
TH1F *h_DeltaNDeltapT_reco[NSECT];
TH1F *h_DeltaNDeltapT_data[NSECT];

//Rebinned raw pT yield
TH1F *h_DeltaNDeltapT_truth_rebinned[NSECT];
TH1F *h_DeltaNDeltapT_reco_rebinned[NSECT];
TH1F *h_DeltaNDeltapT_data_rebinned[NSECT];

//Invariant yield
TH1F *h_dNdpT_truth[NSECT];
TH1F *h_dNdpT_reco[NSECT];
TH1F *h_dNdpT_data[NSECT];

//Acceptance and Efficiency Correction
TH1F *h_correction[NSECT];

//Phi distribution
TH1F *h_phi_reco[NSECT];
TH1F *h_phi_data[NSECT];

//Chisq distribution
TH1F *h_chisqndf_reco[NSECT];
TH1F *h_chisqndf_data[NSECT];

//Pseudorapidity distributions for selected azimuthal regions
TH1F *h_eta_regions_data[4];
TH1F *h_eta_regions_reco[4];

string zvtxcut  = "TMath::Abs(vtx[2]) < 5";
string chisqcut = "chisq/ndf < 3";
string dcacut   = "TMath::Abs(dca) < 0.15";
string dca2dcut = "TMath::Abs(dca2d) < 0.05";
string momcut   = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.2";
string nhitscut = "nhits[0]+nhits[1]+nhits[2]+nhits[3] == 4";
string eta1cut = "TMath::Abs(eta) < 0.35";
string eta2cut = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 0.35";

//Sector cuts
string inclusivephi = "TMath::ATan2(mom[1],mom[0]) < TMath::Pi()";
string ebcut    = "TMath::ATan2(mom[1],mom[0]) < -1*TMath::Pi()/2";
string etcut    = "TMath::ATan2(mom[1],mom[0]) > TMath::Pi()/2";
string wtcut    = "TMath::ATan2(mom[1],mom[0]) > 0 && TMath::ATan2(mom[1],mom[0]) < TMath::Pi()/2";
string wbcut    = "TMath::ATan2(mom[1],mom[0]) < 0 && TMath::ATan2(mom[1],mom[0]) > -1*TMath::Pi()/2";

string sectorLabel[NSECT] = {"INCLUSIVE", "ET", "EB", "WT", "WB"};
string sectorCut[NSECT] = {inclusivephi, etcut, ebcut, wtcut, wbcut};

//Cuts for selected regions with phi discrepancy
string phiregion1cut = "TMath::ATan2(mom[1],mom[0]) > -2.6 && TMath::ATan2(mom[1],mom[0]) < -2.0";
string phiregion2cut = "TMath::ATan2(mom[1],mom[0]) > -1.0 && TMath::ATan2(mom[1],mom[0]) < -0.6";
string phiregion3cut = "TMath::ATan2(mom[1],mom[0]) > -0.2 && TMath::ATan2(mom[1],mom[0]) < 0.2";
string phiregion4cut = "TMath::ATan2(mom[1],mom[0]) > 0.6 && TMath::ATan2(mom[1],mom[0]) < 1.0";

//---------------------------------------
// Functions
//---------------------------------------

void readFiles()
{
	//Read in files
	TFile *f_true = new TFile("Data/ampt_pp_true.root");
	TFile *f_reco    = new TFile("Data/423844_ampt_smeared_105.root");
	TFile *f_data  = new TFile("Data/423844_data_105.root");

	//Extract relevant NTuples
	TTree *ntp_svxseg_true = (TTree*) f_true->Get("ntp_svxseg_true");
	TTree *ntp_svxseg_reco  = (TTree*) f_reco->Get("ntp_svxseg");
	TTree *ntp_event_reco = (TTree*) f_reco->Get("ntp_event");
	TTree *ntp_svxseg_data  = (TTree*) f_data->Get("ntp_svxseg");
	TTree *ntp_event_data   = (TTree*) f_data->Get("ntp_event");

	//Get number of events from NTuples for each case
	nevents_reco = ntp_event_reco->GetEntries();
	nevents_data = ntp_event_data->GetEntries();

	//Extract raw pT distributions for each case
	for (int i = 0; i < NSECT; i++)
	{
		//Spectra
		ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_truth_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), Form("%s && %s && %s", momcut.c_str(), eta1cut.c_str(), sectorCut[i].c_str()), "goff");
		h_DeltaNDeltapT_truth[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_truth_%s", sectorLabel[i].c_str()));

		ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_DeltaNDeltapT_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_reco_%s", sectorLabel[i].c_str()));

		ntp_svxseg_data->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_data_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_DeltaNDeltapT_data[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_data_%s", sectorLabel[i].c_str()));

		//Azimuthal angle distributions
		ntp_svxseg_reco->Draw(Form("TMath::ATan2(mom[1],mom[0])>>h_phi_reco_%s(200,-3.5,3.5)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_phi_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_phi_reco_%s", sectorLabel[i].c_str()));

		ntp_svxseg_data->Draw(Form("TMath::ATan2(mom[1],mom[0])>>h_phi_data_%s(200,-3.5,3.5)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_phi_data[i] = (TH1F*) gDirectory->FindObject(Form("h_phi_data_%s", sectorLabel[i].c_str()));

		//Chi square distributions
		ntp_svxseg_reco->Draw(Form("chisq/ndf>>h_chisqndf_reco_%s(100,0,10)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_chisqndf_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_chisqndf_reco_%s", sectorLabel[i].c_str()));

		ntp_svxseg_data->Draw(Form("chisq/ndf>>h_chisqndf_data_%s(100,0,10)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_chisqndf_data[i] = (TH1F*) gDirectory->FindObject(Form("h_chisqndf_data_%s", sectorLabel[i].c_str()));
	}

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_1(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion1cut).c_str(), "goff");
	h_eta_regions_data[0] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_1");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_2(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion2cut).c_str(), "goff");
	h_eta_regions_data[1] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_2");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_3(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion3cut).c_str(), "goff");
	h_eta_regions_data[2] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_3");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_4(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion4cut).c_str(), "goff");
	h_eta_regions_data[3] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_4");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_1(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion1cut).c_str(), "goff");
	h_eta_regions_reco[0] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_1");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_2(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion2cut).c_str(), "goff");
	h_eta_regions_reco[1] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_2");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_3(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion3cut).c_str(), "goff");
	h_eta_regions_reco[2] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_3");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_4(200,-0.4,0.4)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion4cut).c_str(), "goff");
	h_eta_regions_reco[3] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_4");

	//Set errors on these histograms since they just contain counts
	for (int i = 0; i < NSECT; i++)
	{
		for (int j = 1; j <= NBINS; j++)
		{
			h_DeltaNDeltapT_truth[i]->SetBinError(j, TMath::Sqrt(h_DeltaNDeltapT_truth[i]->GetBinContent(j)));
			h_DeltaNDeltapT_reco[i]->SetBinError(j, TMath::Sqrt(h_DeltaNDeltapT_reco[i]->GetBinContent(j)));
			h_DeltaNDeltapT_data[i]->SetBinError(j, TMath::Sqrt(h_DeltaNDeltapT_data[i]->GetBinContent(j)));
		}
	}
}

void rebinHistograms()
{
	double newBins[NBINS_RB + 1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2.0};

	for (int i = 0; i < NSECT; i++)
	{
		h_DeltaNDeltapT_truth_rebinned[i] = (TH1F*) h_DeltaNDeltapT_truth[i]->Rebin(NBINS_RB, Form("h_DeltaNDeltapT_truth_rebinned_%s", sectorLabel[i].c_str()), newBins);
		h_DeltaNDeltapT_reco_rebinned[i] = (TH1F*) h_DeltaNDeltapT_reco[i]->Rebin(NBINS_RB, Form("h_DeltaNDeltapT_reco_rebinned_%s", sectorLabel[i].c_str()), newBins);
		h_DeltaNDeltapT_data_rebinned[i] = (TH1F*) h_DeltaNDeltapT_data[i]->Rebin(NBINS_RB, Form("h_DeltaNDeltapT_data_rebinned_%s", sectorLabel[i].c_str()), newBins);
	}
}

void normalizeHistograms()
{
	//Normalize by the number of events
	//Divide by the bin width, the acceptance and bin pT
	//Since the Sumw2 structure is already in place, errors will be recomputed correctly
	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i] = (TH1F*) h_DeltaNDeltapT_truth[i]->Clone(Form("h_dNdpT_truth_%s", sectorLabel[i].c_str()));
		h_dNdpT_reco[i] = (TH1F*) h_DeltaNDeltapT_reco[i]->Clone(Form("h_dNdpT_reco_%s", sectorLabel[i].c_str()));
		h_dNdpT_data[i] = (TH1F*) h_DeltaNDeltapT_data[i]->Clone(Form("h_dNdpT_data_%s", sectorLabel[i].c_str()));

		h_dNdpT_truth[i]->Scale(1.0 / nevents_truth);
		h_dNdpT_reco[i]->Scale(1.0 / nevents_reco);
		h_dNdpT_data[i]->Scale(1.0 / nevents_data);

		h_dNdpT_truth[i]->Scale(1.0 / 0.7);
		h_dNdpT_reco[i]->Scale(1.0 / 0.7);
		h_dNdpT_data[i]->Scale(1.0 / 0.7);

		h_dNdpT_truth[i]->Scale(1.0 / 2 * TMath::Pi());
		h_dNdpT_reco[i]->Scale(1.0 / 2 * TMath::Pi());
		h_dNdpT_data[i]->Scale(1.0 / 2 * TMath::Pi());

		if (i > 0)
		{
			h_dNdpT_truth[i]->Scale(4.0);
			h_dNdpT_reco[i]->Scale(4.0);
			h_dNdpT_data[i]->Scale(4.0);
		}

		for (int j = 1; j <= NBINS; j++)
		{
			h_dNdpT_truth[i]->SetBinContent(j, h_dNdpT_truth[i]->GetBinContent(j) / h_dNdpT_truth[i]->GetBinCenter(j));
			h_dNdpT_reco[i]->SetBinContent(j, h_dNdpT_reco[i]->GetBinContent(j) / h_dNdpT_reco[i]->GetBinCenter(j));
			h_dNdpT_data[i]->SetBinContent(j, h_dNdpT_data[i]->GetBinContent(j) / h_dNdpT_data[i]->GetBinCenter(j));
		}

		//Normalize phi distributions to unit integral
		h_phi_data[i]->Scale(1.0 / h_phi_data[i]->Integral());
		h_phi_reco[i]->Scale(1.0 / h_phi_reco[i]->Integral());

		//Normalize chisq distributions to unit at the maximum point
		h_chisqndf_data[i]->Scale(1.0 / h_chisqndf_data[i]->GetMaximum());
		h_chisqndf_reco[i]->Scale(1.0 / h_chisqndf_reco[i]->GetMaximum());
	}

	for (int i = 0; i < 4; i++)
	{
		h_eta_regions_reco[i]->Scale(1.0 / h_eta_regions_reco[i]->Integral());
		h_eta_regions_data[i]->Scale(1.0 / h_eta_regions_data[i]->Integral());
	}
}

void computeCorrection()
{
	for (int i = 0; i < NSECT; i++)
	{
		h_correction[i] = (TH1F*) h_dNdpT_reco[i]->Clone(Form("h_correction_%s", sectorLabel[i].c_str()));
		h_correction[i]->Divide(h_dNdpT_truth[i]);
	}
}

void plot()
{
	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "Truth", 600, 600);
	h_dNdpT_truth[0]->SetLineColor(kBlack);
	h_dNdpT_truth[1]->Draw();

	TCanvas *c2 = new TCanvas("c2", "Reco", 600, 600);
	h_dNdpT_reco[1]->SetLineColor(kBlue);
	h_dNdpT_reco[1]->Draw();

	TCanvas *c3 = new TCanvas("c3", "Data", 600, 600);
	h_dNdpT_data[1]->SetLineColor(kRed);
	h_dNdpT_data[1]->Draw();

	TCanvas *c4 = new TCanvas("c4", "Acceptance and Efficiency Correction");
	h_correction[1]->Draw();
}

void writeToFile()
{
	TFile *fout = new TFile("WorkingFiles/normalizedSpectra_phi_chisq.root", "RECREATE");

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i]->Write();
		h_dNdpT_reco[i]->Write();
		h_dNdpT_data[i]->Write();

		h_correction[i]->Write();

		h_phi_reco[i]->Write();
		h_phi_data[i]->Write();

		h_chisqndf_reco[i]->Write();
		h_chisqndf_data[i]->Write();
	}

	for(int i=0; i<4; i++)
	{	
		h_eta_regions_reco[i]->Write();
		h_eta_regions_data[i]->Write();
	}
}

void RebinMeasuredSpectra()
{
	readFiles();
	//rebinHistograms();
	normalizeHistograms();
	computeCorrection();
	plot();
	writeToFile();
}