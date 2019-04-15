#include <iostream>
#include <fstream>
#include <sstream>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "SNF.h"
SNF::SNF()
{
  RName[0] = "DB1";
  RName[1] = "DB2";
  RName[2] = "LA1";
  RName[3] = "LA2";
  RName[4] = "LA3";
  RName[5] = "LA4";
  RProperName[0] = "Daya Bay 1";
  RProperName[1] = "Daya Bay 2";
  RProperName[2] = "Ling Ao 1";
  RProperName[3] = "Ling Ao 2";
  RProperName[4] = "Ling Ao 3";
  RProperName[5] = "Ling Ao 4";
  vRefuelTimes.resize(kReactor);
  vRefuelFreq.push_back(75);//every 75 weeks for Daya Bay 1
  vRefuelFreq.push_back(75);//every 75 weeks for Daya Bay 2
  vRefuelFreq.push_back(75);//every 75 weeks for Ling Ao 1
  vRefuelFreq.push_back(75);//every 75 weeks for Ling Ao 2
  vRefuelFreq.push_back(52);//every 52 weeks for Ling Ao 3
  vRefuelFreq.push_back(52);//every 52 weeks for Ling Ao 4
}

SNF::~SNF()
{
  if(!RPTree)
    return;
  delete RPTree->GetCurrentFile();
}

double SNF::Alpha(int core){
  if(RName[core].Contains("DB") || RName[core].Contains("LA1")
     || RName[core].Contains("LA2"))
    return 1.0;
  else
    return 1.0;
}

double SNF::Baseline(int det, int core)//distance from reactor to core
{
  return vBaseline[det][core];
}

double SNF::Beta(int core, double E, double days)//fractional contribution to neutrino flux spectrum from a given reactor as a function of days from beginning of P15A
{
  double frac = 0;
  if(vRefuelTimes[core].size()==0)
    FindRefuelTimes(core);
  int week = int(days/double(kDaysInWeek));
  double PofT = ReactorPower(core, week);
  if(PofT < MinOnPower())PofT = 1.0;
  double Pmax = RatedOutputPower(core);
  double alpha = Alpha(core);
  double del = 2.0/(1.0+Delta(core, week));
  double re1_spec = RelSpectrumSim(E, days);
  int n = 0, nv = 0; 
  for(int i=0;i<int(vRefuelTimes[core].size());++i){
    double t =  vRefuelTimes[core][i];
    if(t*kDaysInWeek<days){
      //      cout<<days-t*kDaysInWeek<<" days since refuel.\n";
      frac += alpha*del*Pmax/PofT*RelSpectrumSim(E, days-t*kDaysInWeek);
      if(t>0)++n;
      else ++nv;
    }
  }
  //  cout<<"Summed over "<<n<<" known refueling periods and "
  //  cout<<nv<<" virtual ones.\n";
  return frac;
}
 
int SNF::Configure(const char *filename)//Read parameters from config file
{
  ifstream file(filename);
  if(!(file.good()&&file.is_open())){
    cout<<"Config file not found exiting.\n";
    exit(0);
  }

  string line;
  TString st[kDetector] = {"EH1-AD1","EH1-AD2", "EH2-AD1","EH2-AD2"};
  vBaseline.resize(kDetector);
  for(int idet=0;idet<kDetector;++idet)
    vBaseline[idet].resize(kReactor);

  //Read in baseline distances
  for(int idet=0;idet<kDetector;++idet){
    while(file.good()){
      getline(file, line);
      if(line.find(st[idet].Data())<100){
	char temp[255];
	double distance;
		cout<<st[idet]<<endl;
	for(int icore = 0;icore<kReactor;++icore){
	  file>>temp>>distance;
	  int i=0;
	  while(RName[i].CompareTo(temp)!=0){
	    ++i;
	    if(i==kReactor)break;
	  }
	  if(i>=kReactor){
	    cout<<"Issue finding reactor index. Exiting.\n";
	    exit(0);
	  }
	  vBaseline[idet][i]=distance;
	}
	break;
      }
      file.peek();
      if(file.eof()){
	cout<<"Error. End of config file reached.\n";
	exit(0);
      }
    }
  }
  file.seekg(0);
  //Get rated power output of reactors
  char ln[33] = "Reactor rated output power in GW";
  vReactorPower.resize(kReactor);
  bool found = 0;
  while(file.good()){
    getline(file, line);
    if(line.find("Reactor rated output power in GW")<100){
      found = 1;
      break;
    }
  }
  if(!found){
    cout<<"Reactor power data not found. Exiting.\n";
    exit(0);
  }
  int n = 0;
  char temp[255];
  for(int i=0;i<kReactor;++i){
    double pow = 0;
    file>>temp>>pow;
    for(int j=0;j<kReactor;++j){
      if(RName[j].CompareTo(temp)==0){
	//	cout<<RName[j].Data()<<" reactor power "<<pow<<endl;
	vReactorPower[j] = pow;
	++n;
	break;
      }
    }
  }
  if(n<kReactor){
    cout<<"Not all reactor power data found. Exiting.\n";
    exit(0);
  }
  file.seekg(0);

  while(file.good()){
    getline(file, line);
    if(line.find("#day\\\\MeV")<100)break;
    file.peek();
    if(file.eof()){
      cout<<"Error. End of configuration file reached.\n";
      exit(0);
    }
  }
  getline(file, line);
  stringstream ss(line);

  //Get the energy bin low edges
  while(ss.getline(temp,10,' ')){
    if(temp[0] != '|' && temp[0] != 0){
      vEnergy.push_back(atof(temp));
    }
  }
  //Skip the line of hyphens
  while(getline(file,line)){
    if(line.find_first_of("-----")<100){break;}
    if(file.eof()){
      cout<<"Error. End of configuration file reached.\n";
      exit(0);
    }
  }

  //Get time and spectra from table
  int size = vSpectrum.size();
  while(getline(file,line)){
    if(line.find("-----")<10)break;
    stringstream s(line);
    s.getline(temp,10,' ');
    vTime.push_back(atof(temp));

    vSpectrum.resize(size+1);
    while(s.getline(temp,10,' ')){
      if(temp[0] != '|' && temp[0] != 0){
	vSpectrum[size].push_back(atof(temp));
	//cout<<"j "<<vSpectrum[size].back()<<endl;
      }
    }
    file.peek();
    if(file.eof())break;
    size = vSpectrum.size();
  }
  file.close();
  fNTimeBins = vTime.size();
  cout<<"time bins: "<<fNTimeBins<<endl;
  fNEnergyBins = vEnergy.size();
  return 0;
}

double SNF::Delta(int core, int week)//scale factor to account for the
//greater SNF contribution when the reactor is off for refueling etc.
{
  if(ReactorPower(core, week) > 0.1)return 1.0;
  else return 0.0;
}

int SNF::FindRefuelTimes(int core)//Find refueling times for each reactor
{
  double pow, u235frac;
  RPTree->ResetBranchAddresses();
  RPTree->SetBranchAddress(Form("%s_Power", RName[core].Data()), &pow);
  RPTree->SetBranchAddress(Form("fracU235_%s", RName[core].Data()), &u235frac);
  for(int i=0;i<RPTree->GetEntries();++i){
    RPTree->GetEntry(i);
    double week = i;
    if(pow < 0.1){
      RPTree->GetEntry(i-2);
      double endfrac = u235frac;
      RPTree->GetEntry(i);
      while(pow < 0.5){
	++i;
	if(i>=RPTree->GetEntries())
	  break;
	RPTree->GetEntry(i);
      }
      //      cout<<"Week "<<i<<" "<<u235frac-endfrac<<endl;
      if(u235frac-endfrac > 0.1){
	vRefuelTimes[core].push_back(week);
	cout<<"Refuel core "<<RName[core].Data()<<" at week "<<week<<endl;
      }
    }
  }
  RPTree->ResetBranchAddresses();
  return (int)vRefuelTimes[core].size();
}
TString SNF::GetReactorName(int core)
{
  return RProperName[core];
}

void SNF::InsertRefuelTime(int core)//Use this function to add virtual refuels at a
//set frequency for the years prior to P15A
{
  vRefuelTimes[core].resize(vRefuelTimes[core].size()+1);
  for(int i=(int)vRefuelTimes[core].size()-1;i>0;--i)
    vRefuelTimes[core][i] = vRefuelTimes[core][i-1];
  vRefuelTimes[core][0] = vRefuelTimes[core][1] - vRefuelFreq[core];
}

int SNF::LoadReactorPowerTree(const char *filename, const char *treename)
{
  TFile *file = new TFile(filename);
  RPTree = (TTree*)file->Get(treename);
  // file->Close();
  return (int)RPTree->GetEntries();
}

int SNF::ParameterizeSNFvsT()//parameterize decay versus time for each E-bin
{
  TVirtualFitter::SetPrecision(1e-12);
  //Quadruple exponential decay forces fit through all data points by design.
   TString fnc = "[0]*exp(-[4]*x)+[1]*exp(-[5]*x)+[2]*exp(-[6]*x)"
     "+[3]*exp(-[7]*x)";
  double *x = new double[fNTimeBins];
  double *xe = new double[fNTimeBins];
  double *ye = new double[fNTimeBins];
  for(int it = 0;it<fNTimeBins;++it){
    x[it] = vTime[it];
    xe[it] = 0;
    ye[it] = 0.0001;
  }  
  double *y = new double[fNTimeBins];
  TF1 *f;
  for(int ien = 0;ien<fNEnergyBins;++ien){  
    f = new TF1(Form("f%i",ien), fnc.Data(), 0, 2500);
    for(int it = 0;it<fNTimeBins;++it){
      y[it] = vSpectrum[it][ien];
      xe[it] = 0;
      ye[it] = 0.001*y[it];
    }  
    TGraphErrors gr = TGraphErrors(fNTimeBins, x, y, xe, ye);
    f->SetParameters(0,0,0,0,0,0,0,0);
    f->SetNpx(10000);
    TRandom r(0);
    TFitResultPtr ftres = gr.Fit(f,"SM");
    int n = 0;
    while(int(gr.Fit(f,"S"))==4 && n < 1000){
      for(int i=0; i<f->GetNpar();++i){
	f->SetParameter(i,f->GetParameter(i)*(1+r.Gaus(0,0.1)));
      }
      ++n;
    }
    if(n==1000){
      cout<<"YIKES. Fit did not converge.\n";
      exit(0);
    }
    vSNFofT.push_back(f);
  }
  return (int)vSNFofT.size();
}

double SNF::RatedOutputPower(int core)//reactor rated output power in GW
{
  return vReactorPower[core];
}

double SNF::ReactorPower(int core, int week)
//Reactor power for week since 12/24/2011
{
  RPTree->GetEntry(week);
  return RPTree->GetLeaf(Form("%s_Power", RName[core].Data()))->GetValue();
}

double SNF::RelSpectrumSim(double E, double T)//SNF relative to full power spectrum
{
  //energy bin width
  double E_bw = (vEnergy[fNEnergyBins-1]-vEnergy[0])/double(fNEnergyBins-1);

  //assume 0 for out of range energies
  if(E<vEnergy[0] || E>(vEnergy[fNEnergyBins-1] + E_bw))
    return 0;

  //Use discrete values for energy as in http://arxiv.org/pdf/1512.07353v1.pdf
  int idx_E = (int)((E-vEnergy[0])/E_bw); 

  //Use parameterization of time decay.
  return vSNFofT[idx_E]->Eval(T);

}

double SNF::RelSpectrumAtDet(int det, double E, int week)//fractional SNF flux at E (in MeV) and week since beginning of P15A dataset at specified detector location (not convoluted with IBD X-section)
{
  double frac = 0, norm = 0;
  for(int icore=0;icore<kReactor;++icore){
    //Average over days in week. Relevant only just after refuel when change
    //is rapid. Perhaps at this time even finer time bins should be used.
    for(int t=0;t<1;++t){
      //    for(int t=0;t<kDaysInWeek;++t){
      //double day = t + 0.5 + week * kDaysInWeek;
      double day = t + 0.000001 + week * kDaysInWeek;
      double RPow = ReactorPower(icore, week);
      if(RPow < MinOnPower())RPow = 1.0;
      double weight = RPow/pow(Baseline(det, icore),2);
      frac += weight * Beta(icore, E, day);
      norm += weight;
    }
  }
  //  frac /= (norm * kDaysInWeek);
  frac /= norm;
  return frac;
}
