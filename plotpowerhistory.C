{
  TFile f("WeeklyReactorData.root");
  TTree *tree = (TTree*)f.Get("tree");
  tree->SetLineColor(kRed);
  TCanvas *c = new TCanvas("c","c",0,0,800,1200);
  c->Divide(2,3);
  const int n=6;
  TString str[n] = {"DB1","DB2","LA1","LA2","LA3","LA4"};
  TString name[n] = {"Daya Bay Reactor 1","Daya Bay Reactor 2","Ling Ao Reactor 1","Ling Ao Reactor 2","Ling Ao Reactor 3","Ling Ao Reactor 4"};
  TGraph *h[n];
  gStyle->SetOptStat(0);
  for(int i=0;i<n;++i){
    c->cd(i+1);
    tree->Draw(Form("%s_Power:weeks_since_Dec_24_2011",str[i].Data()),"","l");
    gPad->Update();
    h[i] = (TGraph*)gPad->GetPrimitive("Graph");
    h[i]->Draw("l");
    h[i]->SetTitle(Form("%s Fractional Power",name[i].Data()));
    h[i]->GetYaxis()->SetTitle("Fractional Power");
    gPad->Update();
    tree->Draw(Form("fracU235_%s:weeks_since_Dec_24_2011",str[i].Data()),"","same");
  }
  c->SaveAs("powerhistory.pdf");
}
