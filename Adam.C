{
  gStyle->SetOptStat("e");
  gStyle->SetPadLeftMargin(0.11);
  //first load the tree
  TFile::Open("EH1Po215.root");
  cout<<Po215->GetEntries()<<" entries"<<endl;
  Po215->SetAlias("R","sqrt(py*py+px*px)");
  TCanvas *c = new TCanvas("c","c",0,0,1800,1000);
  c->Divide(3,2);


  c->cd(1);
  //draw events prompt versus delayed
  gPad->SetRightMargin(0.12);
  Po215->Draw("ep:ed>>Ep_vs_Ed(150,0.55,2.05,150,0.59,1.11)","","colz");
  TH2D *hv = (TH2D*)gDirectory->Get("Ep_vs_Ed");
  hv->SetTitle("Number of Events vs E_{prompt} and E_{delayed}");
  gPad->Update();
  hv->GetXaxis()->SetTitle("Delayed Energy (MeV)");
  hv->GetYaxis()->SetTitle("Prompt Energy (MeV)");
  hv->GetXaxis()->SetTitleSize(0.04);
  hv->GetYaxis()->SetTitleSize(0.04);


  c->cd(2);
  gPad->SetRightMargin(0.12);
  //now try to reject Po-215 decays
  Po215->Draw("ep:ed>>hx(150,0.55,2.05,150,0.59,1.11)","pow(ed-0.9,2)+pow(ep-0.8,2)>0.008","colz");
  TH2D *hx = (TH2D*)gDirectory->Get("hx");
  hx->SetTitle("Actinium Event Selection Removed");
  gPad->Update();
  hx->GetXaxis()->SetTitle("Delayed Energy (MeV)");
  hx->GetYaxis()->SetTitle("Prompt Energy (MeV)");
  hx->GetXaxis()->SetTitleSize(0.04);
  hx->GetYaxis()->SetTitleSize(0.04);


  c->cd(3);
  //now try to select only Po-215 decays
  //TCut cut("pow(ed-0.9,2)+pow(ep-0.8,2)<0.004");
  double center_ep = 0.83, center_ed = 0.934;
  TCut cut(Form("pow(ed-%f,2)+pow(ep-%f,2)<0.004", center_ed, center_ep));
  Po215->Draw("ep:ed>>Actinium_Cut(120,0.6,2,120)", cut, "colz");
  TH2D *hi = (TH2D*)gDirectory->Get("Actinium_Cut");
  hi->SetTitle("Actinium Event Selection");
  gPad->Update();
  hi->GetXaxis()->SetTitle("Delayed Energy (MeV)");
  hi->GetYaxis()->SetTitle("Prompt Energy (MeV)");
  hi->GetXaxis()->SetTitleSize(0.04);
  hi->GetYaxis()->SetTitleSize(0.04);
  hi->GetYaxis()->SetTitleOffset(1.1);


  c->cd(4)->SetLogy();
  //check to see if half-life is right
  Po215->Draw("TimeInt/1e6-0.5>>he", cut, "colz");
  TH1D *he = (TH1D*)gDirectory->Get("he");
  he->SetName("Po-215 Decay");
  he->SetTitle("Time Between Successive Events");
  he->GetXaxis()->SetTitle("Time Interval (ms)");
  he->GetYaxis()->SetTitle("Number of Events");
  he->GetXaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetTitleSize(0.04);
  he->Draw();
  gPad->Update();
  TF1 *f = new TF1("f","[0]*exp(log(0.5)*x/[1])",0, 0.495);
  f->SetParName(0,"N_{0}");
  f->SetParName(1,"T_{1/2} (ms)");
  f->SetParName(2,"Constant");
  f->SetParameters(h1->GetBinContent(1),2);
  //make sure function fit shows up on histogram
  gStyle->SetOptFit(1111);
  //Make sure not to include the last bin in the fit
  f->SetParameters(1e4,2);
  he->Fit(f,"r");
  //he->GetYaxis()->SetRangeUser(8500,1.2e4);
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)(he->GetListOfFunctions()->FindObject("stats"));
  ps->SetX1NDC(0.55);
  ps->SetY1NDC(0.65);
  ps->Draw();
 

  c->cd(5);
  //show where actinium events show up
  Po215->Draw("py:px>>hs(100,-1600,1600,100,-1600,1600)", cut, "colz");
  TH2D *hs = (TH2D*)gDirectory->Get("hs");
  gPad->Update();
  hs->SetTitle("XY-position of ^{215}Po Decays in Detector");
  hs->GetXaxis()->SetTitle("X-position (m)");
  hs->GetYaxis()->SetTitle("Y-position (m)");
  hs->GetYaxis()->SetTitleOffset(1.2);
  hs->GetXaxis()->SetTitleSize(0.04);
  hs->GetYaxis()->SetTitleSize(0.04);
  gPad->Update();

  c->cd(6);
  //Is there any r vs z dependence?
  Po215->Draw("pz:R>>hr(100,-100,1600,100,-1600,1600)", cut, "colz");
  TH2D *hr = (TH2D*)gDirectory->Get("hr");
  gPad->Update();
  hr->SetTitle("Radial vs Z-position of ^{215}Po Decays in Detector");
  hr->GetXaxis()->SetTitle("Radial position (m)");
  hr->GetYaxis()->SetTitle("Z-position (m)");
  hr->GetYaxis()->SetTitleOffset(1.2);
  hr->GetXaxis()->SetTitleSize(0.04);
  hr->GetYaxis()->SetTitleSize(0.04);
  gPad->Update();


  TCanvas *c2 = new TCanvas("c2","c2",0,0,1200,1000);
  c2->Divide(2,2);


  c2->cd(1);
  //Is there any r vs energy dependence?
  Po215->Draw("ep:R>>hp", cut,"prof");
  TProfile *hp = (TProfile*)gDirectory->Get("hp");
  hp->SetLineColor(kBlue);
  hp->SetMarkerColor(kBlue);
  hp->Draw();
  gPad->Update();
  hp->SetTitle("E_{prompt} vs Radial Position of ^{215}Po Decays in Detector");
  hp->GetXaxis()->SetTitle("Radial position (m)");
  hp->GetYaxis()->SetTitle("E_{prompt} (MeV)");
  hp->GetXaxis()->SetTitleSize(0.04);
  hp->GetYaxis()->SetTitleSize(0.04);
  hp->GetYaxis()->SetTitleOffset(1.35);
  // hp->GetYaxis()->SetRangeUser(0.792,0.81);
  hp->GetYaxis()->SetRangeUser(center_ep-0.01, center_ep+0.005);
  hp->Draw();
  gPad->Update();


  c2->cd(2);
  Po215->Draw("ed:R>>hd", cut,"prof");
  TProfile *hd = (TProfile*)gDirectory->Get("hd");
  hd->SetMarkerColor(kRed);
  hd->SetLineColor(kRed);
  hd->Draw();
  gPad->Update();
  //  hd->GetYaxis()->SetRangeUser(0.89,0.915);
  hd->GetYaxis()->SetRangeUser(center_ed-0.01, center_ed+0.005);
  hd->SetTitle("Profile Plot of E_{delayed} vs Radial Position of ^{215}Po Decays");
  hd->GetXaxis()->SetTitle("Radial position (m)");
  hd->GetYaxis()->SetTitle("E_{delayed} (MeV)");
  hd->GetXaxis()->SetTitleSize(0.04);
  hd->GetYaxis()->SetTitleSize(0.04);
  hd->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();


  c2->cd(3);
  Po215->Draw("ep:pz>>hpz", cut, "prof");
  TProfile *hpz = (TProfile*)gDirectory->Get("hpz");
  hpz->SetMarkerColor(kBlue);
  hpz->SetLineColor(kBlue);
  hpz->Draw();
  gPad->Update();
  //  hpz->GetYaxis()->SetRangeUser(0.794,0.804);
  hpz->GetYaxis()->SetRangeUser(center_ep-0.01, center_ep+0.005);
  hpz->SetTitle("Profile Plot of E_{prompt} vs Z-position of ^{215}Po Decays");
  hpz->GetXaxis()->SetTitle("Z-position (m)");
  hpz->GetYaxis()->SetTitle("E_{prompt} (MeV)");
  hpz->GetXaxis()->SetTitleSize(0.04);
  hpz->GetYaxis()->SetTitleSize(0.04);
  hpz->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();


  c2->cd(4);
  Po215->Draw("ed:pz>>hdz", cut,"prof");
  TProfile *hdz = (TProfile*)gDirectory->Get("hdz");
  hdz->SetMarkerColor(kRed);
  hdz->SetLineColor(kRed);
  hdz->Draw();
  gPad->Update();
  hdz->GetYaxis()->SetRangeUser(center_ed-0.01, center_ed+0.005);
  hdz->SetTitle("Profile Plot of E_{delayed} vs Z-position of ^{215}Po Decays");
  hdz->GetXaxis()->SetTitle("Z-position (m)");
  hdz->GetYaxis()->SetTitle("E_{delayed} (MeV)");
  hdz->GetXaxis()->SetTitleSize(0.04);
  hdz->GetYaxis()->SetTitleSize(0.04);
  hdz->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();
  c2->ForceUpdate();

}
