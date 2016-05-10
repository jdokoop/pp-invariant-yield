//-----------------------------------
// Macro to plot the charged hadron
// spectra published in PPG30
//-----------------------------------

#include <iostream>
#include <vector>

using namespace std;

//-----------------------------------
// Variables
//-----------------------------------

TGraphErrors *g_pi;

TGraphErrors *g_p;

TGraphErrors *g_k;

TGraphErrors *g_hadrons;

float x_pi[18];
float y_pi[18];
float ex_pi[18];
float ey_pi[18];

float x_p[17];
float y_p[17];
float ex_p[17];
float ey_p[17];

float x_k[13];
float y_k[13];
float ex_k[13];
float ey_k[13];

int totalParticles = 0;

vector<TBox*> pion_systematics;

//-----------------------------------
// Functions
//-----------------------------------

void createPionSpectrum()
{
	float x_plus[18];
	float y_plus[18];
	float ex_plus[18];
	float ey_plus[18];

	float x_minus[18];
	float y_minus[18];
	float ex_minus[18];
	float ey_minus[18];

	//Pi plus
	x_plus[0] = 0.35; y_plus[0] = 25.6287; ex_plus[0] = 0; ey_plus[0] = 0.246673;
	x_plus[1] = 0.45; y_plus[1] = 14.6531; ex_plus[1] = 0; ey_plus[1] = 0.14773;
	x_plus[2] = 0.55; y_plus[2] = 8.62473; ex_plus[2] = 0; ey_plus[2] = 0.0943008;
	x_plus[3] = 0.65; y_plus[3] = 4.80858; ex_plus[3] = 0; ey_plus[3] = 0.0579757;
	x_plus[4] = 0.75; y_plus[4] = 2.85599; ex_plus[4] = 0; ey_plus[4] = 0.03862;
	x_plus[5] = 0.85; y_plus[5] = 1.71022; ex_plus[5] = 0; ey_plus[5] = 0.0258761;
	x_plus[6] = 0.95; y_plus[6] = 1.08056; ex_plus[6] = 0; ey_plus[6] = 0.0183583;
	x_plus[7] = 1.05; y_plus[7] = 0.68761; ex_plus[7] = 0; ey_plus[7] = 0.013222;
	x_plus[8] = 1.15; y_plus[8] = 0.431937; ex_plus[8] = 0; ey_plus[8] = 0.00938769;
	x_plus[9] = 1.25; y_plus[9] = 0.308209; ex_plus[9] = 0; ey_plus[9] = 0.00732418;
	x_plus[10] = 1.35; y_plus[10] = 0.221919; ex_plus[10] = 0; ey_plus[10] = 0.00582217;
	x_plus[11] = 1.45; y_plus[11] = 0.147961; ex_plus[11] = 0; ey_plus[11] = 0.00437927;
	x_plus[12] = 1.6; y_plus[12] = 0.0825923; ex_plus[12] = 0; ey_plus[12] = 0.00205954;
	x_plus[13] = 1.8; y_plus[13] = 0.0410802; ex_plus[13] = 0; ey_plus[13] = 0.00131989;
	x_plus[14] = 2; y_plus[14] = 0.0211058; ex_plus[14] = 0; ey_plus[14] = 0.000876262;
	x_plus[15] = 2.2; y_plus[15] = 0.0116705; ex_plus[15] = 0; ey_plus[15] = 0.000638873;
	x_plus[16] = 2.4; y_plus[16] = 0.00787318; ex_plus[16] = 0; ey_plus[16] = 0.000557839;
	x_plus[17] = 2.6; y_plus[17] = 0.00447392; ex_plus[17] = 0; ey_plus[17] = 0.000458999;


	//Pi minus
	x_minus[0] = 0.35; y_minus[0] = 25.8453; ex_minus[0] = 0; ey_minus[0] = 0.21199;
	x_minus[1] = 0.45; y_minus[1] = 14.5864; ex_minus[1] = 0; ey_minus[1] = 0.127527;
	x_minus[2] = 0.55; y_minus[2] = 8.27867; ex_minus[2] = 0; ey_minus[2] = 0.0789266;
	x_minus[3] = 0.65; y_minus[3] = 4.72545; ex_minus[3] = 0; ey_minus[3] = 0.0498388;
	x_minus[4] = 0.75; y_minus[4] = 2.756; ex_minus[4] = 0; ey_minus[4] = 0.032014;
	x_minus[5] = 0.85; y_minus[5] = 1.66385; ex_minus[5] = 0; ey_minus[5] = 0.0215601;
	x_minus[6] = 0.95; y_minus[6] = 1.05331; ex_minus[6] = 0; ey_minus[6] = 0.0152551;
	x_minus[7] = 1.05; y_minus[7] = 0.648492; ex_minus[7] = 0; ey_minus[7] = 0.0105207;
	x_minus[8] = 1.15; y_minus[8] = 0.431804; ex_minus[8] = 0; ey_minus[8] = 0.0079342;
	x_minus[9] = 1.25; y_minus[9] = 0.286134; ex_minus[9] = 0; ey_minus[9] = 0.00585618;
	x_minus[10] = 1.35; y_minus[10] = 0.199094; ex_minus[10] = 0; ey_minus[10] = 0.00462401;
	x_minus[11] = 1.45; y_minus[11] = 0.133006; ex_minus[11] = 0; ey_minus[11] = 0.00350976;
	x_minus[12] = 1.6; y_minus[12] = 0.0795068; ex_minus[12] = 0; ey_minus[12] = 0.00178784;
	x_minus[13] = 1.8; y_minus[13] = 0.0422071; ex_minus[13] = 0; ey_minus[13] = 0.00119538;
	x_minus[14] = 2; y_minus[14] = 0.0208344; ex_minus[14] = 0; ey_minus[14] = 0.000781061;
	x_minus[15] = 2.2; y_minus[15] = 0.0131348; ex_minus[15] = 0; ey_minus[15] = 0.000625054;
	x_minus[16] = 2.4; y_minus[16] = 0.0070146; ex_minus[16] = 0; ey_minus[16] = 0.000479789;
	x_minus[17] = 2.6; y_minus[17] = 0.00425537; ex_minus[17] = 0; ey_minus[17] = 0.000408;

	//Systematics
	float esys_pi[18] = {0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.156016, 0.153756, 0.153756, 0.164913, 0.180853, 0.200503, 0.222883};

	//Divide by 42 to get (1/2pi)(1/pT)dN/dpTdy
	for (int i = 0; i < 18; i++)
	{
		x_pi[i] = x_minus[i];
		y_pi[i] = (y_minus[i] + y_plus[i]) / 42.0;

		ex_pi[i] = ex_minus[i];
		ey_pi[i] = 0;//TMath::Sqrt(ey_minus[i] * ey_minus[i] + ey_plus[i] * ey_plus[i]);

		totalParticles += (y_pi[i] * 2 * TMath::Pi() * x_pi[i]);

		TBox *box = new TBox(x_pi[i] - 0.025, y_pi[i] - esys_pi[i]*y_pi[i], x_pi[i] + 0.025, y_pi[i] + esys_pi[i]*y_pi[i]);
		box->SetLineColor(kGreen + 3);
		box->SetFillColorAlpha(kWhite, 0.0);
		pion_systematics.push_back(box);
	}

	g_pi = new TGraphErrors(18, x_pi, y_pi, ex_pi, ey_pi);
	g_pi->SetName("g_pi");
}

void createProtonSpectrum()
{
	float x_plus[17];
	float y_plus[17];
	float ex_plus[17];
	float ey_plus[17];

	float x_minus[17];
	float y_minus[17];
	float ex_minus[17];
	float ey_minus[17];

	//Protons
	x_plus[0] = 0.65; y_plus[0] = 0.467596; ex_plus[0] = 0; ey_plus[0] = 0.0124957;
	x_plus[1] = 0.75; y_plus[1] = 0.361103; ex_plus[1] = 0; ey_plus[1] = 0.0100168;
	x_plus[2] = 0.85; y_plus[2] = 0.250201; ex_plus[2] = 0; ey_plus[2] = 0.00749516;
	x_plus[3] = 0.95; y_plus[3] = 0.181417; ex_plus[3] = 0; ey_plus[3] = 0.00593336;
	x_plus[4] = 1.05; y_plus[4] = 0.126166; ex_plus[4] = 0; ey_plus[4] = 0.00445335;
	x_plus[5] = 1.15; y_plus[5] = 0.0902262; ex_plus[5] = 0; ey_plus[5] = 0.00351448;
	x_plus[6] = 1.25; y_plus[6] = 0.0651584; ex_plus[6] = 0; ey_plus[6] = 0.00275751;
	x_plus[7] = 1.35; y_plus[7] = 0.0448966; ex_plus[7] = 0; ey_plus[7] = 0.00214243;
	x_plus[8] = 1.45; y_plus[8] = 0.0328418; ex_plus[8] = 0; ey_plus[8] = 0.00172412;
	x_plus[9] = 1.6; y_plus[9] = 0.0203516; ex_plus[9] = 0; ey_plus[9] = 0.000881656;
	x_plus[10] = 1.8; y_plus[10] = 0.0111527; ex_plus[10] = 0; ey_plus[10] = 0.000593141;
	x_plus[11] = 2; y_plus[11] = 0.00477154; ex_plus[11] = 0; ey_plus[11] = 0.000353859;
	x_plus[12] = 2.2; y_plus[12] = 0.00309521; ex_plus[12] = 0; ey_plus[12] = 0.000265977;
	x_plus[13] = 2.4; y_plus[13] = 0.0018804; ex_plus[13] = 0; ey_plus[13] = 0.0001962;
	x_plus[14] = 2.6; y_plus[14] = 0.00122164; ex_plus[14] = 0; ey_plus[14] = 0.000150502;
	x_plus[15] = 2.9; y_plus[15] = 0.000420402; ex_plus[15] = 0; ey_plus[15] = 5.69246e-05;
	x_plus[16] = 3.4; y_plus[16] = 8.66954e-05; ex_plus[16] = 0; ey_plus[16] = 1.85193e-05;

	//Antiprotons
	x_minus[0] = 0.65; y_minus[0] = 0.306682; ex_minus[0] = 0; ey_minus[0] = 0.00840118;
	x_minus[1] = 0.75; y_minus[1] = 0.241129; ex_minus[1] = 0; ey_minus[1] = 0.00669015;
	x_minus[2] = 0.85; y_minus[2] = 0.175801; ex_minus[2] = 0; ey_minus[2] = 0.00512706;
	x_minus[3] = 0.95; y_minus[3] = 0.124961; ex_minus[3] = 0; ey_minus[3] = 0.00389379;
	x_minus[4] = 1.05; y_minus[4] = 0.0902856; ex_minus[4] = 0; ey_minus[4] = 0.00301377;
	x_minus[5] = 1.15; y_minus[5] = 0.0675941; ex_minus[5] = 0; ey_minus[5] = 0.00250251;
	x_minus[6] = 1.25; y_minus[6] = 0.0451574; ex_minus[6] = 0; ey_minus[6] = 0.00186547;
	x_minus[7] = 1.35; y_minus[7] = 0.0329953; ex_minus[7] = 0; ey_minus[7] = 0.00152812;
	x_minus[8] = 1.45; y_minus[8] = 0.0219945; ex_minus[8] = 0; ey_minus[8] = 0.00118843;
	x_minus[9] = 1.6; y_minus[9] = 0.0151914; ex_minus[9] = 0; ey_minus[9] = 0.000644838;
	x_minus[10] = 1.8; y_minus[10] = 0.00740469; ex_minus[10] = 0; ey_minus[10] = 0.000417221;
	x_minus[11] = 2; y_minus[11] = 0.0037345; ex_minus[11] = 0; ey_minus[11] = 0.000274817;
	x_minus[12] = 2.2; y_minus[12] = 0.00251093; ex_minus[12] = 0; ey_minus[12] = 0.000213118;
	x_minus[13] = 2.4; y_minus[13] = 0.00116624; ex_minus[13] = 0; ey_minus[13] = 0.000138554;
	x_minus[14] = 2.6; y_minus[14] = 0.000707305; ex_minus[14] = 0; ey_minus[14] = 0.000103754;
	x_minus[15] = 2.9; y_minus[15] = 0.000333785; ex_minus[15] = 0; ey_minus[15] = 4.60119e-05;
	x_minus[16] = 3.4; y_minus[16] = 6.25373e-05; ex_minus[16] = 0; ey_minus[16] = 1.43608e-05;

	//Divide by 42 mb
	for (int i = 0; i < 17; i++)
	{
		x_p[i] = x_minus[i];
		y_p[i] = (y_minus[i] + y_plus[i]) / 42.0;

		ex_p[i] = ex_minus[i];
		ey_p[i] = TMath::Sqrt(ey_minus[i] * ey_minus[i] + ey_plus[i] * ey_plus[i]);

		totalParticles += (y_p[i] * 2 * TMath::Pi() * x_p[i]);
	}

	g_p = new TGraphErrors(17, x_p, y_p, ex_p, ey_p);
	g_p->SetName("g_p");
}

void createKaonSpectrum()
{
	float x_plus[13];
	float y_plus[13];
	float ex_plus[13];
	float ey_plus[13];

	float x_minus[13];
	float y_minus[13];
	float ex_minus[13];
	float ey_minus[13];

	//K+
	x_plus[0] = 0.45; y_plus[0] = 1.95293; ex_plus[0] = 0; ey_plus[0] = 0.0724524;
	x_plus[1] = 0.55; y_plus[1] = 1.40404; ex_plus[1] = 0; ey_plus[1] = 0.0466304;
	x_plus[2] = 0.65; y_plus[2] = 0.925078; ex_plus[2] = 0; ey_plus[2] = 0.0296139;
	x_plus[3] = 0.75; y_plus[3] = 0.624339; ex_plus[3] = 0; ey_plus[3] = 0.020472;
	x_plus[4] = 0.85; y_plus[4] = 0.403909; ex_plus[4] = 0; ey_plus[4] = 0.0142644;
	x_plus[5] = 0.95; y_plus[5] = 0.281567; ex_plus[5] = 0; ey_plus[5] = 0.0107128;
	x_plus[6] = 1.05; y_plus[6] = 0.199669; ex_plus[6] = 0; ey_plus[6] = 0.00826777;
	x_plus[7] = 1.15; y_plus[7] = 0.137853; ex_plus[7] = 0; ey_plus[7] = 0.00620796;
	x_plus[8] = 1.25; y_plus[8] = 0.101317; ex_plus[8] = 0; ey_plus[8] = 0.00483626;
	x_plus[9] = 1.35; y_plus[9] = 0.0694105; ex_plus[9] = 0; ey_plus[9] = 0.00362463;
	x_plus[10] = 1.45; y_plus[10] = 0.0526523; ex_plus[10] = 0; ey_plus[10] = 0.00301769;
	x_plus[11] = 1.6; y_plus[11] = 0.0324873; ex_plus[11] = 0; ey_plus[11] = 0.00151623;
	x_plus[12] = 1.8; y_plus[12] = 0.0179392; ex_plus[12] = 0; ey_plus[12] = 0.00100964;

	//K-
	x_minus[0] = 0.45; y_minus[0] = 1.66952; ex_minus[0] = 0; ey_minus[0] = 0.0568774;
	x_minus[1] = 0.55; y_minus[1] = 1.25238; ex_minus[1] = 0; ey_minus[1] = 0.0374863;
	x_minus[2] = 0.65; y_minus[2] = 0.828785; ex_minus[2] = 0; ey_minus[2] = 0.0243554;
	x_minus[3] = 0.75; y_minus[3] = 0.604526; ex_minus[3] = 0; ey_minus[3] = 0.0177572;
	x_minus[4] = 0.85; y_minus[4] = 0.406696; ex_minus[4] = 0; ey_minus[4] = 0.0125949;
	x_minus[5] = 0.95; y_minus[5] = 0.261138; ex_minus[5] = 0; ey_minus[5] = 0.00878815;
	x_minus[6] = 1.05; y_minus[6] = 0.189842; ex_minus[6] = 0; ey_minus[6] = 0.00685834;
	x_minus[7] = 1.15; y_minus[7] = 0.118611; ex_minus[7] = 0; ey_minus[7] = 0.00483604;
	x_minus[8] = 1.25; y_minus[8] = 0.0977945; ex_minus[8] = 0; ey_minus[8] = 0.00407271;
	x_minus[9] = 1.35; y_minus[9] = 0.0673309; ex_minus[9] = 0; ey_minus[9] = 0.00315524;
	x_minus[10] = 1.45; y_minus[10] = 0.0439384; ex_minus[10] = 0; ey_minus[10] = 0.00239946;
	x_minus[11] = 1.6; y_minus[11] = 0.0283855; ex_minus[11] = 0; ey_minus[11] = 0.00125486;
	x_minus[12] = 1.8; y_minus[12] = 0.0181311; ex_minus[12] = 0; ey_minus[12] = 0.000911162;

	//Divide by 42 mb
	for (int i = 0; i < 13; i++)
	{
		x_k[i] = x_minus[i];
		y_k[i] = (y_minus[i] + y_plus[i]) / 42.0;

		ex_k[i] = ex_minus[i];
		ey_k[i] = TMath::Sqrt(ey_minus[i] * ey_minus[i] + ey_plus[i] * ey_plus[i]);

		totalParticles += (y_k[i] * 2 * TMath::Pi() * x_k[i]);
	}

	g_k = new TGraphErrors(13, x_k, y_k, ex_k, ey_k);
	g_k->SetName("g_k");
}

void createHadronSpectrum()
{
	float x[11];
	float y[11];
	float ex[11];
	float ey[11];

	for (int i = 0; i < 12; i++)
	{
		x[i]  = x_pi[i + 3];
		y[i]  = y_pi[i + 3] + y_p[i] + y_k[i + 2];
		ex[i] = ex_pi[i];
		ey[i] = TMath::Sqrt(ex_pi[i] * ex_pi[i] + ex_p[i] * ex_p[i] + ex_k[i] * ex_k[i]);
	}

	g_hadrons = new TGraphErrors(11, x, y, ex, ey);
}

void plotHadronComparison()
{
	TH1F *h_ampt_pt;
	TFile *f_ampt = new TFile("Data/ampt_pp_true.root");
	TTree *ntp_ampt = (TTree*) f_ampt->Get("ntp_svxseg_true");

	ntp_ampt->Draw("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h(100,0,2)", "TMath::Abs(y) < 0.5", "goff");
	h_ampt_pt = (TH1F*) gDirectory->FindObject("h");

	h_ampt_pt->Scale(1.0 / 2000000);            //Number of events
	h_ampt_pt->Scale(1.0 / 0.02);              //Bin width
	h_ampt_pt->Scale(1.0 / (2 * TMath::Pi())); //2Pi factor

	for (int i = 1; i <= h_ampt_pt->GetNbinsX(); i++)
	{
		float binCenter = h_ampt_pt->GetBinCenter(i);

		float binContent = h_ampt_pt->GetBinContent(i);
		binContent = binContent / binCenter;

		h_ampt_pt->SetBinContent(i, binContent);
	}

	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c", "c", 600, 600);

	g_hadrons->SetLineColor(kOrange + 3);
	g_hadrons->SetLineWidth(4);
	g_hadrons->SetMarkerStyle(20);
	h_ampt_pt->Draw();
	g_hadrons->Draw("CP,same");
	h_ampt_pt->SetMarkerStyle(20);
	h_ampt_pt->SetMarkerColor(kBlue);
	h_ampt_pt->SetMarkerSize(0.6);
	h_ampt_pt->SetTitle("");
	h_ampt_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_ampt_pt->GetXaxis()->SetLabelFont(62);
	h_ampt_pt->GetXaxis()->SetTitleFont(62);
	h_ampt_pt->GetYaxis()->SetTitle("(1/2#pi p_{T})(1/N_{ev})dN/dp_{T}dy");
	h_ampt_pt->GetYaxis()->SetLabelFont(62);
	h_ampt_pt->GetYaxis()->SetTitleFont(62);
}

void plotPionComparison()
{
	TH1F *h_ampt_pt;
	TFile *f_ampt = new TFile("Data/ampt_pp_true_1.root");
	TTree *ntp_ampt = (TTree*) f_ampt->Get("ntp_svxseg_true");

	ntp_ampt->Draw("TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1])>>h(100,0,2)", "TMath::Abs(y) < 0.5 && species == 211", "goff");
	h_ampt_pt = (TH1F*) gDirectory->FindObject("h");

	h_ampt_pt->Scale(1.0 / 2000000);            //Number of events
	h_ampt_pt->Scale(1.0 / 0.02);              //Bin width
	h_ampt_pt->Scale(1.0 / (2 * TMath::Pi())); //2Pi factor

	for (int i = 1; i <= h_ampt_pt->GetNbinsX(); i++)
	{
		float binCenter = h_ampt_pt->GetBinCenter(i);

		float binContent = h_ampt_pt->GetBinContent(i);
		binContent = binContent / binCenter;

		h_ampt_pt->SetBinContent(i, binContent);
	}

	gStyle->SetOptStat(0);
	TCanvas *cPions = new TCanvas("cPions", "cPions", 600, 600);

	h_ampt_pt->SetMarkerStyle(20);
	h_ampt_pt->SetMarkerColor(kBlue);
	h_ampt_pt->SetMarkerSize(0.6);
	h_ampt_pt->SetTitle("");
	h_ampt_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_ampt_pt->GetXaxis()->SetLabelFont(62);
	h_ampt_pt->GetXaxis()->SetTitleFont(62);
	h_ampt_pt->GetYaxis()->SetTitle("(1/2#pi p_{T})(1/N_{ev})dN/dp_{T}dy");
	h_ampt_pt->GetYaxis()->SetLabelFont(62);
	h_ampt_pt->GetYaxis()->SetTitleFont(62);

	g_pi->SetLineWidth(2);
	g_pi->SetMarkerStyle(20);
	g_pi->SetLineColor(kSpring - 6);
	g_pi->SetMarkerColor(kSpring - 6);

	//Fit with Tsallis functional form
	float m = 0.140;
	TF1 *myfit = new TF1("myfit", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.140*([1] - 1))*([1]*[2] + 0.140)) * pow(([1]*[2] + TMath::Sqrt(0.140*0.140 + x*x))/([1]*[2]+0.140),-1*[1])", 0, 1.9);
	myfit->SetParameter(0, 1.28);
	myfit->SetParameter(1, 9.67);
	myfit->SetParameter(2, 121.077);
	//g_pi->Fit("myfit","R");

	h_ampt_pt->Draw("P");
	g_pi->Draw("CP,same");

	//Plot systematics
	for (int i = 0; i < 18; i++)
	{
		TBox *b = (TBox*) pion_systematics[i];
		b->Draw("same,l");
	}

	TLegend *tleg = new TLegend(0.55, 0.7, 0.85, 0.8);
	tleg->AddEntry(g_pi, "PPG030", "P");
	tleg->AddEntry(h_ampt_pt, "AMPT", "P");
	tleg->SetLineColor(kWhite);
	tleg->Draw("same");

	TLatex *tl1_1 = new TLatex(0.15, 0.4, "p+p at 200 GeV");
	tl1_1->SetNDC(kTRUE);
	tl1_1->Draw("same");

	TLatex *tl1_2 = new TLatex(0.15, 0.3, "Charged Pions");
	tl1_2->SetNDC(kTRUE);
	tl1_2->Draw("same");

}

void PlotPublishedSpectrum()
{
	//Take the data points from PPG30 and create a TGraph for every hadron species
	createPionSpectrum();
	createProtonSpectrum();
	createKaonSpectrum();

	//Combine all hadron species into a single charged hadron spectrum
	createHadronSpectrum();

	//Compare the published data with the measurements from Run15 data
	plotPionComparison();
	//plotHadronComparison();
}