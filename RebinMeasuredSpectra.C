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

//Number of sectors
// --> Inclusive phi, ET, EB, WT, WB
const int NSECT = 5;

//Number of events in each file
int nevents_truth = 0;
int nevents_reco  = 0;
int nevents_data  = 0;

//General flags
bool reqPrecise = !true; //Require a precise vertex in reconstructed AMPT
bool reqSeed    = true;  //Require a seed and no precise vertex in reconstructed AMPT

//Raw pT yield
TH1F *h_DeltaNDeltapT_truth[NSECT];
TH1F *h_DeltaNDeltapT_reco[NSECT];
TH1F *h_DeltaNDeltapT_data[NSECT];

//Invariant yield
TH1F *h_dNdpT_truth[NSECT];
TH1F *h_dNdpT_reco[NSECT];
TH1F *h_dNdpT_data[NSECT];

//Acceptance and Efficiency Correction
TH1F *h_correction[NSECT];

//General track cuts
string precvtxcut   = "vtx_prec[2] == vtx_prec[2]";
string seedvtxcut   = "vtx_seed[2] == vtx_seed[2] && vtx_prec[2] != vtx_prec[2]";
string zvtxcut      = "TMath::Abs(vtx_bbc[2]) < 10";
string zvtxcuttruth = "TMath::Abs(vtx_bbc[2]) < 10";
string chisqcut     = "chisq/ndf < 3";
string dcacut       = "TMath::Abs(dca) < 0.15";
string dca2dcut     = "TMath::Abs(dca2d) < 0.05";
string momcut       = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0";
string momtruthcut  = "TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1]) > 0";
string nhitscut     = "nhits[0] > 0 && nhits[1]> 0 && nhits[2] > 0 && nhits[3] > 0";
string eta1cut      = "TMath::Abs(eta) < 0.35";
string eta2cut      = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 0.35";
string etatruthcut  = "TMath::Abs(TMath::ATanH(mom_truth[2]/TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1] + mom_truth[2]*mom_truth[2]))) < 0.35";
string bbccut       = "pmtbbcs > 0 && pmtbbcn > 0";

//Sector cuts
string inclusivephi = "TMath::ATan2(mom[1],mom[0]) < TMath::Pi()";
string ebcut        = "TMath::ATan2(mom[1],mom[0]) < -1*TMath::Pi()/2";
string etcut        = "TMath::ATan2(mom[1],mom[0]) > TMath::Pi()/2";
string wtcut        = "TMath::ATan2(mom[1],mom[0]) > 0 && TMath::ATan2(mom[1],mom[0]) < TMath::Pi()/2";
string wbcut        = "TMath::ATan2(mom[1],mom[0]) > -TMath::Pi()/2 && TMath::ATan2(mom[1],mom[0]) < 0";

string inclusivephitruth = "TMath::ATan2(mom_truth[1],mom_truth[0]) < TMath::Pi()";
string ebcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) < -1*TMath::Pi()/2";
string etcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) > TMath::Pi()/2";
string wtcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) > 0 && TMath::ATan2(mom_truth[1],mom_truth[0]) < TMath::Pi()/2";
string wbcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) > -TMath::Pi()/2 && TMath::ATan2(mom_truth[1],mom_truth[0]) < 0";


string sectorLabel[NSECT]    = {"INCLUSIVE", "ET", "EB", "WT", "WB"};
string sectorCut[NSECT]      = {inclusivephi, etcut, ebcut, wtcut, wbcut};
string sectorCutTruth[NSECT] = {inclusivephitruth, etcuttruth, ebcuttruth, wtcuttruth, wbcuttruth};

//TTrees to read in from file
TTree *ntp_svxseg_true;
TTree *ntp_svxseg_reco;
TTree *ntp_event_reco;
TTree *ntp_svxseg_data;
TTree *ntp_event_data;

//Files to read
TFile *f_reco;
TFile *f_data;

//---------------------------------------
// Functions
//---------------------------------------

void getTruthInformation()
{
	//The truth information is contained in the npt_svxseg_ampt tree
	ntp_svxseg_true = (TTree*) f_reco->Get("ntp_svxseg_ampt");

	for (int i = 0; i < NSECT; i++)
	{
		if (reqPrecise)
		{
			ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1])>>h_DeltaNDeltapT_truth_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), Form("%s && %s && %s && %s && %s && %s", momtruthcut.c_str(), etatruthcut.c_str(), sectorCutTruth[i].c_str(), bbccut.c_str(), zvtxcut.c_str(), precvtxcut.c_str()), "goff");
		}
		else if (reqSeed)
		{
			ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1])>>h_DeltaNDeltapT_truth_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), Form("%s && %s && %s && %s && %s && %s", momtruthcut.c_str(), etatruthcut.c_str(), sectorCutTruth[i].c_str(), bbccut.c_str(), zvtxcut.c_str(), seedvtxcut.c_str()), "goff");
		}
		else
		{
			ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1])>>h_DeltaNDeltapT_truth_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), Form("%s && %s && %s && %s && %s", momtruthcut.c_str(), etatruthcut.c_str(), sectorCutTruth[i].c_str(), bbccut.c_str(), zvtxcut.c_str()), "goff");
		}

		h_DeltaNDeltapT_truth[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_truth_%s", sectorLabel[i].c_str()));
	}

	//Count events that fire the bbc trigger and have a narrow vertex
	float vtx_bbc[3];
	float vtx_prec[3];
	float vtx_seed[3];
	int eventno;
	int pmtbbcs;
	int pmtbbcn;

	ntp_svxseg_true->SetBranchAddress("vtx_prec", &vtx_prec);
	ntp_svxseg_true->SetBranchAddress("vtx_bbc", &vtx_bbc);
	ntp_svxseg_true->SetBranchAddress("vtx_seed", &vtx_seed);
	ntp_svxseg_true->SetBranchAddress("eventno", &eventno);
	ntp_svxseg_true->SetBranchAddress("pmtbbcs", &pmtbbcs);
	ntp_svxseg_true->SetBranchAddress("pmtbbcn", &pmtbbcn);

	int lastEvent = -9999;
	for (int i = 0; i < ntp_svxseg_true->GetEntries(); i++)
	{
		ntp_svxseg_true->GetEntry(i);

		if (eventno != lastEvent && TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcs > 0 && pmtbbcn > 0)
		{
			if (reqPrecise && vtx_prec[2] != vtx_prec[2])
			{
				continue;
			}

			if (reqSeed && vtx_seed[2] != vtx_seed[2])
			{
				continue;
			}

			nevents_truth++;
			lastEvent = eventno;
		}
	}
}

void readFiles()
{
	//Read in files
	f_reco    = new TFile("Data/423844_reco_1_0_1_0_050416.root");
	f_data  = new TFile("Data/423844_data_0_1_050216.root");

	//Extract relevant NTuples
	ntp_svxseg_reco = (TTree*) f_reco->Get("ntp_svxseg");
	ntp_event_reco  = (TTree*) f_reco->Get("ntp_event");
	ntp_svxseg_data = (TTree*) f_data->Get("ntp_svxseg");
	ntp_event_data  = (TTree*) f_data->Get("ntp_event");

	//Get number of events from NTuples for each case

	nevents_data = ntp_event_data->GetEntries(Form("%s && %s", bbccut.c_str(), zvtxcut.c_str()));

	if (reqPrecise)
	{
		nevents_reco = ntp_event_reco->GetEntries(Form("%s && %s && %s", bbccut.c_str(), zvtxcut.c_str(), precvtxcut.c_str()));
	}
	else if (reqSeed)
	{
		nevents_reco = ntp_event_reco->GetEntries(Form("%s && %s && %s", bbccut.c_str(), zvtxcut.c_str(), seedvtxcut.c_str()));
	}
	else
	{
		nevents_reco = ntp_event_reco->GetEntries(Form("%s && %s", bbccut.c_str(), zvtxcut.c_str()));
	}

	//Extract the spectra from truth AMPT events that fire the BBC trigger
	getTruthInformation();

	//Extract raw pT distributions for each case
	for (int i = 0; i < NSECT; i++)
	{
		ntp_svxseg_data->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_data_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i]).c_str(), "goff");
		h_DeltaNDeltapT_data[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_data_%s", sectorLabel[i].c_str()));

		if (reqPrecise)
		{
			ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i] + "&&" + precvtxcut).c_str(), "goff");
			h_DeltaNDeltapT_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_reco_%s", sectorLabel[i].c_str()));
		}
		else if (reqSeed)
		{
			ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i] + "&&" + seedvtxcut).c_str(), "goff");
			h_DeltaNDeltapT_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_reco_%s", sectorLabel[i].c_str()));
		}
		else
		{
			ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i]).c_str(), "goff");
			h_DeltaNDeltapT_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_reco_%s", sectorLabel[i].c_str()));
		}
	}

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

void normalizeHistograms()
{
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

		h_dNdpT_truth[i]->Scale(1.0 / (2 * TMath::Pi()));
		h_dNdpT_reco[i]->Scale(1.0 / (2 * TMath::Pi()));
		h_dNdpT_data[i]->Scale(1.0 / (2 * TMath::Pi()));

		for (int j = 1; j <= NBINS; j++)
		{
			h_dNdpT_truth[i]->SetBinContent(j, h_dNdpT_truth[i]->GetBinContent(j) / h_dNdpT_truth[i]->GetBinWidth(j));
			h_dNdpT_reco[i]->SetBinContent(j, h_dNdpT_reco[i]->GetBinContent(j) / h_dNdpT_reco[i]->GetBinWidth(j));
			h_dNdpT_data[i]->SetBinContent(j, h_dNdpT_data[i]->GetBinContent(j) / h_dNdpT_data[i]->GetBinWidth(j));
		}

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
	TFile *fout;
	if (reqPrecise)
	{
		fout = new TFile("WorkingFiles/normSpectra_temp_prec.root", "RECREATE");
	}
	else if (reqSeed)
	{
		fout = new TFile("WorkingFiles/normSpectra_temp_seed.root", "RECREATE");
	}
	else
	{
		fout = new TFile("WorkingFiles/normSpectra_temp.root", "RECREATE");
	}

	for (int i = 0; i < NSECT; i++)
	{
		h_dNdpT_truth[i]->Write();
		h_dNdpT_reco[i]->Write();
		h_dNdpT_data[i]->Write();
		h_correction[i]->Write();
	}
}

void RebinMeasuredSpectra()
{
	readFiles();
	normalizeHistograms();
	computeCorrection();
	plot();
	writeToFile();
}