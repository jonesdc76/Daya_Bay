#include <vector>
#include <map>
#include "TROOT.h"
#include "TTree.h"
#include "TF1.h"

class SNF
{
 public:
  SNF();
  ~SNF();
  static const int kDetector = 4;
  static const int kReactor = 6;
  static const int kDaysInWeek = 7;
  static const double kMinOnPower;
  static double MinOnPower(){return 0.01;}
  vector<double>vReactorPower;
  vector<double>vRefuelFreq;
  vector<double>vEnergy;
  vector<double>vTime;
  vector<vector<double> >vRefuelTimes;//[core][week]
  vector<vector<double> >vBaseline;//[detector][core]
  vector<vector<double> >vSpectrum;//[time][energy]
  double Alpha(int);//scale factor accounts for SNF with irradiation != 45GWd/t
  double Baseline(int, int);//distance from reactor to core
  double Beta(int, double, double);//fractional contribution to neutrino flux spectrum from a given reactor as a function of energy (in GeV) and time (in days) 
  double Delta(int, int);//scale factor to account for reactor off SNF
  int Configure(const char *);//Read parameters from config file
  int FindRefuelTimes(int);//Find refueling times for each reactor
  TString GetReactorName(int);//Useful for labeling graphs
  void InsertRefuelTime(int);//function to add refuel times before start of P15A
  int LoadReactorPowerTree(const char *, const char *);//Root tree of reactor power
  int ParameterizeSNFvsT();//parameterize decay versus time for each energy bin
  double RatedOutputPower(int);//reactor rated output power in GW
  double ReactorPower(int, int);//Reactor power as function of core# and week#
  double RelSpectrumSim(double, double);//Simulated fractional SNF spectrum of E (in MeV) and time (in days) 
  double RelSpectrumAtDet(int, double, int);//fractional SNF flux at E (in MeV) and t (in days) at specified detector location (not convoluted with IBD X-section)
private:
  TString RName[kReactor];
  TString RProperName[kReactor];
  TTree *RPTree;
  int fNEnergyBins;
  int fNTimeBins;
  vector<TF1*>vSNFofT;
};
