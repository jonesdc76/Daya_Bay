{
  c->cd(1);
  gr->Fit(f,"MR");
  c->cd(2);
  for(int i=0;i<fNTimeBins;++i)
    yr[i] = (y[i] - f->Eval(x[i]))/y[i]*100;
  grres = new TGraphErrors(fNTimeBins, x, yr, xe, ye);
  grres->Draw("ap");
  grres->SetMarkerColor(kBlue);
  grres->SetMarkerStyle(20);
  grres->SetLineColor(kBlue);
  cout<<Form("%13.7f",snf.vEnergy[en])<<endl;
  cout<<"-------------"<<endl;
  for(int i=0; i<8;++i)
    cout<<Form("%13.7f",f->GetParameter(i))<<endl;
}
