{
  SNF snf;
  int fNEnergyBins=9;
  snf.Configure("test.config");
  int fNTimeBins = 8;
  TString fnc="[0]*exp(-[4]*x)+[1]*exp(-[5]*x)+[2]*exp(-[6]*x)+[3]*exp(-[7]*x)";

  TF1 *f = new TF1("f", fnc.Data(), 0, 2500);
  f->SetParameters(0,0,0,0,0,0,0,0,0);
  //f->SetParameters(0.15,0.54,0.64,0,1.28,0,0.1,0.1);
  f->SetNpx(10000);
  double *x = new double[fNTimeBins];
  double *xe = new double[fNTimeBins];
  double *ye = new double[fNTimeBins];
  for(int it = 0;it<fNTimeBins;++it){
    x[it] = snf.vTime[it];
    xe[it] = 0;
  }  
  double *y = new double[fNTimeBins];
  double err = 0.0001;
  en = 8;
  for(int it = 0;it<fNTimeBins;++it){
    y[it] = snf.vSpectrum[it][en];
    //ye[it] = (it>5? 0.37*err:err);
    ye[it] = err;
    if(it>5)ye[it]=err*0.4;
    //if(it>4)ye[it]*=0.5;
  }
 
  TCanvas *c = new TCanvas("c","c",0,0,1400,600);
  c->Divide(2,1);
  c->cd(1)->SetLogx();
  TGraphErrors *gr = new TGraphErrors(fNTimeBins, x, y, xe, ye);
  gr->Draw("ap");
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gr->SetLineColor(kRed);
  TVirtualFitter::SetPrecision(1e-12);

  gr->Fit(f,"MR");
  gr->Fit(f,"MR");
  gr->Fit(f,"MR");
  c->cd(2);
  double *yr = new double[fNTimeBins];
  for(int i=0;i<fNTimeBins;++i)
    yr[i] = (y[i] - f->Eval(x[i]))/y[i]*100;
  TGraphErrors *grres = new TGraphErrors(fNTimeBins, x, yr, xe, ye);
  grres->Draw("ap");
  grres->SetMarkerColor(kBlue);
  grres->SetMarkerStyle(20);
  grres->SetLineColor(kBlue);

}  
