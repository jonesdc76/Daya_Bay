#include "SNF.C"
#include "TGraph.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TCanvas.h"


int SNFcontribution(int det, int week)
{
  SNF snf;
  int nCore = snf.kReactor;
  snf.Configure("test.config");
  snf.LoadReactorPowerTree("WeeklyReactorData.root", "tree");
  snf.ParameterizeSNFvsT();
  for(int icore=0;icore<nCore;++icore)
    snf.FindRefuelTimes(icore);

  //Add refuels prior to P15A data at expected rate for each reactor

  //Add 7 refuels (~10 yr)for Daya Bay 1,2 and Ling Ao 1,2
  for(int i=0;i<7;++i){
    snf.InsertRefuelTime(0);
    snf.InsertRefuelTime(1);
    snf.InsertRefuelTime(2);
    snf.InsertRefuelTime(3);
  }
  //Add 5 refuels (~5 yr)for Ling Ao 3,4
  for(int i=0;i<7;++i){
    snf.InsertRefuelTime(4);
    snf.InsertRefuelTime(5);
  }
  int n = int(snf.vEnergy.size());
  cout<<n<<endl;
  double *energy = snf.vEnergy.data();
  double size = energy[1] - energy[0];
  double *relativeFlux = new double[n];
  for(int i=0;i<n;++i){
    relativeFlux[i] = snf.RelSpectrumAtDet(det, energy[i], week);
  }
  TGraph *gr = new TGraph(n, energy, relativeFlux);
  gr->Draw("ap");
  gr->SetTitle(Form("Relative SNF Antineutrino Flux at Detector %i for Week %i", det, week));
  gROOT->SetStyle("Plain");
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gPad->Update();
  gr->GetXaxis()->SetTitle("Antineutrino Energy (MeV)");
  gr->GetYaxis()->SetTitle("SNF Flux/Reactor Flux (%)");
  gPad->Update();
  return 0;
}
