#include <iostream>

using namespace std;

//---------------------------------------
// Variables
//---------------------------------------

TGraph *g_resolution;

//---------------------------------------
// Functions
//---------------------------------------

void fitResolution()
{
	//Take pT resolution for 4-hit tracks from Theo
	//Fit the points to a polynomial

	float pT[5] = {0.25, 0.76, 1.25, 1.75, 2.23};
	float res[5] = {0.109, 0.110, 0.113, 0.119, 0.125};

	g_resolution = new TGraph(5, pT, res);

	g_resolution->Draw("ACP");
}

void InvestigateResolutionEffect()
{
	fitResolution();
}