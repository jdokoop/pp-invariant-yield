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

//Phi distribution
TH1F *h_phi_reco[NSECT];
TH1F *h_phi_data[NSECT];

//Chisq distribution
TH1F *h_chisqndf_reco[NSECT];
TH1F *h_chisqndf_data[NSECT];

//Pseudorapidity distributions for |zprec| < 1
TH1F *h_eta_preccut_data;
TH1F *h_eta_preccut_reco;

//Pseudorapidity distributions for selected azimuthal regions
TH1F *h_eta_regions_data[4];
TH1F *h_eta_regions_reco[4];

//Eta vs phi distribution
TH2F *h_eta_phi_data;
TH2F *h_eta_phi_sims;

//General track cuts
string zvtxpreccut  = "TMath::Abs(vtx_prec[2]) < 1";
string zvtxcut      = "TMath::Abs(vtx_bbc[2]) < 10";
string zvtxcuttruth = "TMath::Abs(vtx_bbc[2]) < 10";
string chisqcut     = "chisq/ndf < 3";
string dcacut       = "TMath::Abs(dca) < 0.15";
string dca2dcut     = "TMath::Abs(dca2d) < 0.05";
//string momcut       = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0.2";
//string momtruthcut  = "TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1]) > 0.2";
string momcut       = "TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]) > 0";
string momtruthcut  = "TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1]) > 0";
string nhitscut     = "nhits[0] > 0 && nhits[1]> 0 && nhits[2] > 0 && nhits[3] > 0";
string eta1cut      = "TMath::Abs(eta) < 0.35";
string eta2cut      = "TMath::Abs(TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))) < 0.35";
string etatruthcut  = "TMath::Abs(TMath::ATanH(mom_truth[2]/TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1] + mom_truth[2]*mom_truth[2]))) < 0.35";
string bbccut       = "pmtbbcs > 0 && pmtbbcn > 0";

//Track cut to remove anomaly in ET (2 < phi(wrapped around = not wrapped around) < 2.5)
string anomalycutet = "TMath::ATan2(mom[1],mom[0]) < 2 && TMath::ATan2(mom[1],mom[0]) > 2.25";

//Sector cuts
string inclusivephi = "TMath::ATan2(mom[1],mom[0]) < TMath::Pi()";
string ebcut        = "TMath::ATan2(mom[1],mom[0]) < -1*TMath::Pi()/2";
string etcut        = "TMath::ATan2(mom[1],mom[0]) > TMath::Pi()/2";
string wtcut        = "TMath::ATan2(mom[1],mom[0]) > 0 && TMath::ATan2(mom[1],mom[0]) < TMath::Pi()/2";
string wbcut        = "TMath::ATan2(mom[1],mom[0]) < 0 && TMath::ATan2(mom[1],mom[0]) > -1*TMath::Pi()/2";

string inclusivephitruth = "TMath::ATan2(mom_truth[1],mom_truth[0]) < TMath::Pi()";
string ebcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) < -1*TMath::Pi()/2";
string etcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) > TMath::Pi()/2";
string wtcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) > 0 && TMath::ATan2(mom_truth[1],mom_truth[0]) < TMath::Pi()/2";
string wbcuttruth        = "TMath::ATan2(mom_truth[1],mom_truth[0]) < 0 && TMath::ATan2(mom_truth[1],mom_truth[0]) > -1*TMath::Pi()/2";

string sectorLabel[NSECT]    = {"INCLUSIVE", "ET", "EB", "WT", "WB"};
string sectorCut[NSECT]      = {inclusivephi, etcut, ebcut, wtcut, wbcut};
string sectorCutTruth[NSECT] = {inclusivephitruth, etcuttruth, ebcuttruth, wtcuttruth, wbcuttruth};

//Cuts for selected regions with phi discrepancy
string phiregion1cut = "TMath::ATan2(mom[1],mom[0]) > 2.4 && TMath::ATan2(mom[1],mom[0]) < 3.1";
string phiregion2cut = "TMath::ATan2(mom[1],mom[0]) > -0.2 && TMath::ATan2(mom[1],mom[0]) < 0.2";
string phiregion3cut = "TMath::ATan2(mom[1],mom[0]) > 0.65 && TMath::ATan2(mom[1],mom[0]) < 0.85";
string phiregion4cut = "TMath::ATan2(mom[1],mom[0]) > -0.75 && TMath::ATan2(mom[1],mom[0]) < -0.6";

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
	//Select only tracks for events with BBCS + BBCN > 2 and bbcz_vtx < 1 cm
	ntp_svxseg_true = (TTree*) f_reco->Get("ntp_svxseg_ampt");

	for (int i = 0; i < NSECT; i++)
	{
		ntp_svxseg_true->Draw(Form("TMath::Sqrt(mom_truth[0]*mom_truth[0] + mom_truth[1]*mom_truth[1])>>h_DeltaNDeltapT_truth_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), Form("%s && %s && %s && %s && %s", momtruthcut.c_str(), etatruthcut.c_str(), sectorCutTruth[i].c_str(), bbccut.c_str(), zvtxcut.c_str()), "goff");
		h_DeltaNDeltapT_truth[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_truth_%s", sectorLabel[i].c_str()));
	}

	//Count events that fire the bbc trigger and have a narrow vertex
	float vtx[3];
	int eventno;
	int pmtbbcs;
	int pmtbbcn;
	ntp_svxseg_true->SetBranchAddress("vtx_bbc", &vtx);
	ntp_svxseg_true->SetBranchAddress("eventno", &eventno);
	ntp_svxseg_true->SetBranchAddress("pmtbbcs", &pmtbbcs);
	ntp_svxseg_true->SetBranchAddress("pmtbbcn", &pmtbbcn);

	int lastEvent = -9999;
	for (int i = 0; i < ntp_svxseg_true->GetEntries(); i++)
	{
		ntp_svxseg_true->GetEntry(i);

		if (eventno != lastEvent && TMath::Abs(vtx[2]) < 10 && pmtbbcs > 0 && pmtbbcn > 0)
		{
			nevents_truth++;
			lastEvent = eventno;
		}
	}

	cout << "NEVENTS TRUTH = " << nevents_truth << endl;
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
	nevents_reco = ntp_event_reco->GetEntries(Form("%s && %s", bbccut.c_str(), zvtxcut.c_str()));
	nevents_data = ntp_event_data->GetEntries(Form("%s && %s", bbccut.c_str(), zvtxcut.c_str()));

	cout << "NEVENTS RECO = " << nevents_reco << endl;
	cout << "NEVENTS DATA = " << nevents_data << endl << endl;

	//Extract the spectra from truth AMPT events that fire the BBC trigger
	getTruthInformation();

	//Extract raw pT distributions for each case
	for (int i = 0; i < NSECT; i++)
	{
		ntp_svxseg_reco->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_reco_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i]).c_str(), "goff");
		h_DeltaNDeltapT_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_reco_%s", sectorLabel[i].c_str()));

		ntp_svxseg_data->Draw(Form("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h_DeltaNDeltapT_data_%s(%i,0.2,2)", sectorLabel[i].c_str(), NBINS), (eta2cut + "&&" + zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut + "&&" + sectorCut[i]).c_str(), "goff");
		h_DeltaNDeltapT_data[i] = (TH1F*) gDirectory->FindObject(Form("h_DeltaNDeltapT_data_%s", sectorLabel[i].c_str()));

		//Chi square distributions
		ntp_svxseg_reco->Draw(Form("chisq/ndf>>h_chisqndf_reco_%s(100,0,10)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i] + "&&" + bbccut).c_str(), "goff");
		h_chisqndf_reco[i] = (TH1F*) gDirectory->FindObject(Form("h_chisqndf_reco_%s", sectorLabel[i].c_str()));

		ntp_svxseg_data->Draw(Form("chisq/ndf>>h_chisqndf_data_%s(100,0,10)", sectorLabel[i].c_str()), (eta2cut + "&&" + zvtxcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + sectorCut[i]).c_str(), "goff");
		h_chisqndf_data[i] = (TH1F*) gDirectory->FindObject(Form("h_chisqndf_data_%s", sectorLabel[i].c_str()));
	}

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_1(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion1cut).c_str(), "goff");
	h_eta_regions_data[0] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_1");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_2(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion2cut).c_str(), "goff");
	h_eta_regions_data[1] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_2");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_3(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion3cut).c_str(), "goff");
	h_eta_regions_data[2] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_3");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_data_4(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion4cut).c_str(), "goff");
	h_eta_regions_data[3] = (TH1F*) gDirectory->FindObject("h_eta_regions_data_4");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_1(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion1cut + "&&" + bbccut).c_str(), "goff");
	h_eta_regions_reco[0] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_1");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_2(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion2cut + "&&" + bbccut).c_str(), "goff");
	h_eta_regions_reco[1] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_2");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_3(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion3cut + "&&" + bbccut).c_str(), "goff");
	h_eta_regions_reco[2] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_3");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_regions_reco_4(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + phiregion4cut + "&&" + bbccut).c_str(), "goff");
	h_eta_regions_reco[3] = (TH1F*) gDirectory->FindObject("h_eta_regions_reco_4");

	ntp_svxseg_reco->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_preccut_reco(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut).c_str(), "goff");
	h_eta_preccut_reco = (TH1F*) gDirectory->FindObject("h_eta_preccut_reco");

	ntp_svxseg_data->Draw("TMath::ATanH(mom[2]/TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]))>>h_eta_preccut_data(200,-1.5,1.5)", (zvtxcut + "&&" + chisqcut + "&&" + dcacut + "&&" + dca2dcut + "&&" + nhitscut + "&&" + momcut + "&&" + bbccut).c_str(), "goff");
	h_eta_preccut_data = (TH1F*) gDirectory->FindObject("h_eta_preccut_data");

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

	cout << "Raw Entries Reco  = " << h_DeltaNDeltapT_reco[0]->Integral() << endl;
	cout << "Raw Entries Data  = " << h_DeltaNDeltapT_data[0]->Integral() << endl;
	cout << "Raw Entries Truth = " << h_DeltaNDeltapT_truth[0]->Integral() << endl << endl;
}

void extractAzimuthalDistribution()
{
	h_phi_data[0] = new TH1F("h_phi_data_INCLUSIVE", "h_phi_data_INCLUSIVE", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_data[1] = new TH1F("h_phi_data_ET", "h_phi_data_ET", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_data[2] = new TH1F("h_phi_data_EB", "h_phi_data_EB", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_data[3] = new TH1F("h_phi_data_WT", "h_phi_data_WT", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_data[4] = new TH1F("h_phi_data_WB", "h_phi_data_WB", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

	h_phi_reco[0] = new TH1F("h_phi_reco_INCLUSIVE", "h_phi_reco_INCLUSIVE", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_reco[1] = new TH1F("h_phi_reco_ET", "h_phi_reco_ET", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_reco[2] = new TH1F("h_phi_reco_EB", "h_phi_reco_EB", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_reco[3] = new TH1F("h_phi_reco_WT", "h_phi_reco_WT", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_phi_reco[4] = new TH1F("h_phi_reco_WB", "h_phi_reco_WB", 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

	h_eta_phi_sims = new TH2F("h_eta_phi_sims", "h_eta_phi_sims;#eta;#phi", 200, -1.5, 1.5, 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	h_eta_phi_data = new TH2F("h_eta_phi_data", "h_eta_phi_data;#eta;#phi", 200, -1.5, 1.5, 400, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

	float mom_reco[3];
	float vtx_reco[3];
	int nhits_reco[4];
	float chisq_reco;
	float ndf_reco;
	float dca_reco;
	float dca2d_reco;
	int pmtbbcn_reco;
	int pmtbbcs_reco;

	float mom_data[3];
	float vtx_data[3];
	int nhits_data[4];
	float chisq_data;
	float ndf_data;
	float dca_data;
	float dca2d_data;

	ntp_svxseg_reco->SetBranchAddress("mom", &mom_reco);
	ntp_svxseg_reco->SetBranchAddress("vtx_bbc", &vtx_reco);
	ntp_svxseg_reco->SetBranchAddress("nhits", &nhits_reco);
	ntp_svxseg_reco->SetBranchAddress("chisq", &chisq_reco);
	ntp_svxseg_reco->SetBranchAddress("ndf", &ndf_reco);
	ntp_svxseg_reco->SetBranchAddress("dca", &dca_reco);
	ntp_svxseg_reco->SetBranchAddress("dca2d", &dca2d_reco);
	ntp_svxseg_reco->SetBranchAddress("pmtbbcs", &pmtbbcs_reco);
	ntp_svxseg_reco->SetBranchAddress("pmtbbcn", &pmtbbcn_reco);

	for (int i = 0; i < ntp_svxseg_reco->GetEntries(); i++)
	{
		ntp_svxseg_reco->GetEntry(i);

		float phi = TMath::ATan2(mom_reco[1], mom_reco[0]);
		float eta = TMath::ATanH(mom_reco[2] / TMath::Sqrt(mom_reco[0] * mom_reco[0] + mom_reco[1] * mom_reco[1] + mom_reco[2] * mom_reco[2]));

		if ((nhits_reco[0] + nhits_reco[1] + nhits_reco[2] + nhits_reco[3] != 4) || TMath::Abs(vtx_reco[2]) > 10 || TMath::Abs(dca_reco) > 0.15 || TMath::Abs(dca2d_reco) > 0.05 || chisq_reco / ndf_reco > 3 || TMath::Sqrt(mom_reco[0]*mom_reco[0] + mom_reco[1]*mom_reco[1]) < 0.2 || TMath::Abs(TMath::ATanH(mom_reco[2] / TMath::Sqrt(mom_reco[0]*mom_reco[0] + mom_reco[1]*mom_reco[1] + mom_reco[2]*mom_reco[2]))) > 0.35 || pmtbbcn_reco == 0 || pmtbbcs_reco == 0)
		{
			continue;
		}

		if (phi > 1.5 * TMath::Pi())
		{
			phi = phi - 2 * TMath::Pi();
		}
		else if (phi < -0.5 * TMath::Pi())
		{
			phi = phi + 2 * TMath::Pi();
		}

		h_phi_reco[0]->Fill(phi);
		h_eta_phi_sims->Fill(eta, phi);
	}

	ntp_svxseg_data->SetBranchAddress("mom", &mom_data);
	ntp_svxseg_data->SetBranchAddress("vtx_bbc", &vtx_data);
	ntp_svxseg_data->SetBranchAddress("nhits", &nhits_data);
	ntp_svxseg_data->SetBranchAddress("chisq", &chisq_data);
	ntp_svxseg_data->SetBranchAddress("ndf", &ndf_data);
	ntp_svxseg_data->SetBranchAddress("dca", &dca_data);
	ntp_svxseg_data->SetBranchAddress("dca2d", &dca2d_data);

	for (int i = 0; i < ntp_svxseg_data->GetEntries(); i++)
	{
		ntp_svxseg_data->GetEntry(i);

		float phi = TMath::ATan2(mom_data[1], mom_data[0]);
		float eta = TMath::ATanH(mom_data[2] / TMath::Sqrt(mom_data[0] * mom_data[0] + mom_data[1] * mom_data[1] + mom_data[2] * mom_data[2]));

		if ((nhits_data[0] + nhits_data[1] + nhits_data[2] + nhits_data[3] != 4) || TMath::Abs(vtx_data[2]) > 1 || TMath::Abs(dca_data) > 0.15 || TMath::Abs(dca2d_data) > 0.05 || chisq_data / ndf_data > 3 || TMath::Sqrt(mom_data[0]*mom_data[0] + mom_data[1]*mom_data[1]) < 0.2 || TMath::Abs(TMath::ATanH(mom_data[2] / TMath::Sqrt(mom_data[0]*mom_data[0] + mom_data[1]*mom_data[1] + mom_data[2]*mom_data[2]))) > 0.35)
		{
			continue;
		}

		if (phi > 1.5 * TMath::Pi())
		{
			phi = phi - 2 * TMath::Pi();
		}
		else if (phi < -0.5 * TMath::Pi())
		{
			phi = phi + 2 * TMath::Pi();
		}

		h_phi_data[0]->Fill(phi);
		h_eta_phi_data->Fill(eta, phi);
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

		if (i == 0)
		{
			cout << "Ntracks_reco after event norm = " << h_dNdpT_reco[i]->Integral() << endl;
			cout << "Ntracks_data after event norm = " << h_dNdpT_data[i]->Integral() << endl;
			cout << "Ntracks_truth after event norm = " << h_dNdpT_truth[i]->Integral() << endl << endl;
		}

		h_dNdpT_truth[i]->Scale(1.0 / 0.7);
		h_dNdpT_reco[i]->Scale(1.0 / 0.7);
		h_dNdpT_data[i]->Scale(1.0 / 0.7);

		if (i == 0)
		{
			cout << "Ntracks_reco after Deta norm = " << h_dNdpT_reco[i]->Integral() << endl;
			cout << "Ntracks_data after Deta norm = " << h_dNdpT_data[i]->Integral() << endl;
			cout << "Ntracks_truth after Deta norm = " << h_dNdpT_truth[i]->Integral() << endl << endl;
		}

		h_dNdpT_truth[i]->Scale(1.0 / (2 * TMath::Pi()));
		h_dNdpT_reco[i]->Scale(1.0 / (2 * TMath::Pi()));
		h_dNdpT_data[i]->Scale(1.0 / (2 * TMath::Pi()));

		if (i == 0)
		{
			cout << "Ntracks_reco after 2pi norm = " << h_dNdpT_reco[i]->Integral() << endl;
			cout << "Ntracks_data after 2pi norm = " << h_dNdpT_data[i]->Integral() << endl;
			cout << "Ntracks_truth after 2pi norm = " << h_dNdpT_truth[i]->Integral() << endl << endl;
		}

		for (int j = 1; j <= NBINS; j++)
		{
			h_dNdpT_truth[i]->SetBinContent(j, h_dNdpT_truth[i]->GetBinContent(j) / h_dNdpT_truth[i]->GetBinWidth(j));
			h_dNdpT_reco[i]->SetBinContent(j, h_dNdpT_reco[i]->GetBinContent(j) / h_dNdpT_reco[i]->GetBinWidth(j));
			h_dNdpT_data[i]->SetBinContent(j, h_dNdpT_data[i]->GetBinContent(j) / h_dNdpT_data[i]->GetBinWidth(j));
		}

		if (i == 0)
		{
			cout << "Ntracks_reco after norm = " << h_dNdpT_reco[i]->Integral() << endl;
			cout << "Ntracks_data after norm = " << h_dNdpT_data[i]->Integral() << endl;
			cout << "Ntracks_truth after norm = " << h_dNdpT_truth[i]->Integral() << endl << endl;
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

		//Normalize phi distributions to unit integral
		h_phi_data[i]->Scale(1.0 / h_phi_data[i]->Integral());
		h_phi_reco[i]->Scale(1.0 / h_phi_reco[i]->Integral());

		//Normalize chisq distributions to unit at the maximum point
		h_chisqndf_data[i]->Scale(1.0 / h_chisqndf_data[i]->GetMaximum());
		h_chisqndf_reco[i]->Scale(1.0 / h_chisqndf_reco[i]->GetMaximum());
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
	TFile *fout = new TFile("WorkingFiles/normSpectra_1_0_1_0_v_0_1_050216.root", "RECREATE");

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

	for (int i = 0; i < 4; i++)
	{
		h_eta_regions_reco[i]->Write();
		h_eta_regions_data[i]->Write();
	}

	h_eta_phi_sims->Write();
	h_eta_phi_data->Write();

	h_eta_preccut_reco->Write();
	h_eta_preccut_data->Write();
}

void RebinMeasuredSpectra()
{
	readFiles();
	extractAzimuthalDistribution();
	normalizeHistograms();
	computeCorrection();
	plot();
	writeToFile();
}