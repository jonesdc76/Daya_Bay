{
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",0,0,1800,600);
  c->Divide(3,1);
  TChain *ch = new TChain("tree");
  ch->Add("WeeklyReactorData.root");
  c->cd(1)->SetGrid();
  ch->Draw("fracU235_DB1:Entry$-100>>hDB1(75,0,75,1000,0,0.94)","Entry$>=100&&Entry$<168");
  TH2D *hDB1 = (TH2D*)gDirectory->Get("hDB1");
  hDB1->Draw();
  hDB1->SetTitle("Fuel Fractions Daya Bay Reactors");
  hDB1->GetYaxis()->SetTitle("Fuel Fractions");
  hDB1->GetYaxis()->SetTitleOffset(1.3);
  hDB1->GetXaxis()->SetTitle("Weeks Since Refuel");
  hDB1->SetMarkerStyle(kMultiply);
  hDB1->SetMarkerColor(kRed);
  ch->Draw("fracU238_DB1:Entry$-100>>hDB1_1(75,0,75,1000,0,0.94)","Entry$>=100&&Entry$<168", "goff");
  TH2D *hDB1_1 = (TH2D*)gDirectory->Get("hDB1_1");
  hDB1_1->SetMarkerStyle(kMultiply);
  hDB1_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_DB1:Entry$-100>>hDB1_2(75,0,75,1000,0,0.94)","Entry$>=100&&Entry$<168","goff");
  TH2D *hDB1_2 = (TH2D*)gDirectory->Get("hDB1_2");
  hDB1_2->SetMarkerStyle(kMultiply);
  hDB1_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_DB1:Entry$-100>>hDB1_3(75,0,75,1000,0,0.94)","Entry$>=100&&Entry$<168","goff");
  TH2D *hDB1_3 = (TH2D*)gDirectory->Get("hDB1_3");
  hDB1_3->SetMarkerStyle(kMultiply);
  hDB1_3->SetMarkerColor(kBlack);
  hDB1_1->Draw("same");
  hDB1_2->Draw("same");
  hDB1_3->Draw("same");

  gPad->Update();


  ch->Draw("fracU235_DB2:Entry$-75>>hDB2(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<143","goff");
  TH2D *hDB2 = (TH2D*)gDirectory->Get("hDB2");
  // hDB2->Draw();
  // hDB2->SetTitle("Fuel Fractions Daya Bay Reactor 2");
  // hDB2->GetYaxis()->SetTitle("Fuel Fractions");
  // hDB2->GetXaxis()->SetTitle("Weeks Since Refuel");
  hDB2->SetMarkerStyle(kPlus);
  hDB2->SetMarkerColor(kRed);
  ch->Draw("fracU238_DB2:Entry$-75>>hDB2_1(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<142", "goff");
  TH2D *hDB2_1 = (TH2D*)gDirectory->Get("hDB2_1");
  hDB2_1->SetMarkerStyle(kPlus);
  hDB2_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_DB2:Entry$-75>>hDB2_2(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<142","goff");
  TH2D *hDB2_2 = (TH2D*)gDirectory->Get("hDB2_2");
  hDB2_2->SetMarkerStyle(kPlus);
  hDB2_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_DB2:Entry$-75>>hDB2_3(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<142","goff");
  TH2D *hDB2_3 = (TH2D*)gDirectory->Get("hDB2_3");
  hDB2_3->SetMarkerStyle(kPlus);
  hDB2_3->SetMarkerColor(kBlack);
  hDB2->Draw("same");
  hDB2_1->Draw("same");
  hDB2_2->Draw("same");
  hDB2_3->Draw("same");

  gPad->Update();
  TPaveText *pt[8];
  double size = 0.06;
  TString name[8] = {"DB1  ^{235}U: #times", "DB1  ^{238}U: #times","DB1 ^{239}Pu: #times","DB1 ^{241}Pu: #times","DB2  ^{235}U: +", "DB2  ^{238}U: +","DB2 ^{239}Pu: +","DB2 ^{241}Pu: +"};
  int col[4] = {kRed, kMagenta, kBlue, kBlack};
  for(int i=0;i<8;++i){
    pt[i]= new TPaveText(0.83,0.83-i*size,0.99,0.89-i*size,"ndc");
    pt[i]->SetFillColor(0);
    pt[i]->SetShadowColor(0);
    pt[i]->SetBorderSize(0);
    pt[i]->SetTextColor(col[i%4]);
    pt[i]->AddText(name[i].Data());
    pt[i]->Draw();
  }


  c->cd(2)->SetGrid();
  ch->Draw("fracU235_LA1:Entry$-65>>hLA1(75,0,75,1000,0,0.94)","Entry$>=65&&Entry$<125");
  TH2D *hLA1 = (TH2D*)gDirectory->Get("hLA1");
  hLA1->Draw();
  hLA1->SetTitle("Fuel Fractions Ling Ao 1,2 Reactors");
  hLA1->GetYaxis()->SetTitle("Fuel Fractions");
  hLA1->GetYaxis()->SetTitleOffset(1.3);
  hLA1->GetXaxis()->SetTitle("Weeks Since Refuel");
  hLA1->SetMarkerStyle(kMultiply);
  hLA1->SetMarkerColor(kRed);
  ch->Draw("fracU238_LA1:Entry$-65>>hLA1_1(75,0,75,1000,0,0.94)","Entry$>=65&&Entry$<125", "goff");
  TH2D *hLA1_1 = (TH2D*)gDirectory->Get("hLA1_1");
  hLA1_1->SetMarkerStyle(kMultiply);
  hLA1_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_LA1:Entry$-65>>hLA1_2(75,0,75,1000,0,0.94)","Entry$>=65&&Entry$<125","goff");
  TH2D *hLA1_2 = (TH2D*)gDirectory->Get("hLA1_2");
  hLA1_2->SetMarkerStyle(kMultiply);
  hLA1_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_LA1:Entry$-65>>hLA1_3(75,0,75,1000,0,0.94)","Entry$>=65&&Entry$<125","goff");
  TH2D *hLA1_3 = (TH2D*)gDirectory->Get("hLA1_3");
  hLA1_3->SetMarkerStyle(kMultiply);
  hLA1_3->SetMarkerColor(kBlack);
  hLA1_1->Draw("same");
  hLA1_2->Draw("same");
  hLA1_3->Draw("same");

  gPad->Update();

  //  c->cd(2)->SetGrid();
  ch->Draw("fracU235_LA2:Entry$-108>>hLA2(75,0,75,1000,0,0.94)","Entry$>=108&&Entry$<172","goff");
  TH2D *hLA2 = (TH2D*)gDirectory->Get("hLA2");
  // hLA2->Draw();
  // hLA2->SetTitle("Fuel Fractions Daya Bay Reactor 2");
  // hLA2->GetYaxis()->SetTitle("Fuel Fractions");
  // hLA2->GetXaxis()->SetTitle("Weeks Since Refuel");
  hLA2->SetMarkerStyle(kPlus);
  hLA2->SetMarkerColor(kRed);
  ch->Draw("fracU238_LA2:Entry$-108>>hLA2_1(75,0,75,1000,0,0.94)","Entry$>=108&&Entry$<172", "goff");
  TH2D *hLA2_1 = (TH2D*)gDirectory->Get("hLA2_1");
  hLA2_1->SetMarkerStyle(kPlus);
  hLA2_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_LA2:Entry$-108>>hLA2_2(75,0,75,1000,0,0.94)","Entry$>=108&&Entry$<172","goff");
  TH2D *hLA2_2 = (TH2D*)gDirectory->Get("hLA2_2");
  hLA2_2->SetMarkerStyle(kPlus);
  hLA2_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_LA2:Entry$-108>>hLA2_3(75,0,75,1000,0,0.94)","Entry$>=108&&Entry$<172","goff");
  TH2D *hLA2_3 = (TH2D*)gDirectory->Get("hLA2_3");
  hLA2_3->SetMarkerStyle(kPlus);
  hLA2_3->SetMarkerColor(kBlack);
  hLA2->Draw("same");
  hLA2_1->Draw("same");
  hLA2_2->Draw("same");
  hLA2_3->Draw("same");

  gPad->Update();
  TPaveText *ptLA[8];
  TString nameLA[8] = {"LA1  ^{235}U: #times", "LA1  ^{238}U: #times","LA1 ^{239}Pu: #times","LA1 ^{241}Pu: #times","LA2  ^{235}U: +", "LA2  ^{238}U: +","LA2 ^{239}Pu: +","LA2 ^{241}Pu: +"};
  for(int i=0;i<8;++i){
    ptLA[i]= new TPaveText(0.83,0.83-i*size,0.99,0.89-i*size,"ndc");
    ptLA[i]->SetFillColor(0);
    ptLA[i]->SetShadowColor(0);
    ptLA[i]->SetBorderSize(0);
    ptLA[i]->SetTextColor(col[i%4]);
    ptLA[i]->AddText(nameLA[i].Data());
    ptLA[i]->Draw();
  }

  c->cd(3)->SetGrid();
  ch->Draw("fracU235_LA3:Entry$-75>>hLA3(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<115");
  TH2D *hLA3 = (TH2D*)gDirectory->Get("hLA3");
  hLA3->Draw();
  hLA3->SetTitle("Fuel Fractions Ling Ao 3,4 Reactors");
  hLA3->GetYaxis()->SetTitle("Fuel Fractions");
  hLA3->GetYaxis()->SetTitleOffset(1.3);
  hLA3->GetXaxis()->SetTitle("Weeks Since Refuel");
  hLA3->SetMarkerStyle(kMultiply);
  hLA3->SetMarkerColor(kRed);
  ch->Draw("fracU238_LA3:Entry$-75>>hLA3_1(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<115", "goff");
  TH2D *hLA3_1 = (TH2D*)gDirectory->Get("hLA3_1");
  hLA3_1->SetMarkerStyle(kMultiply);
  hLA3_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_LA3:Entry$-75>>hLA3_2(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<115","goff");
  TH2D *hLA3_2 = (TH2D*)gDirectory->Get("hLA3_2");
  hLA3_2->SetMarkerStyle(kMultiply);
  hLA3_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_LA3:Entry$-75>>hLA3_3(75,0,75,1000,0,0.94)","Entry$>=75&&Entry$<115","goff");
  TH2D *hLA3_3 = (TH2D*)gDirectory->Get("hLA3_3");
  hLA3_3->SetMarkerStyle(kMultiply);
  hLA3_3->SetMarkerColor(kBlack);
  hLA3_1->Draw("same");
  hLA3_2->Draw("same");
  hLA3_3->Draw("same");

  gPad->Update();

  //  c->cd(2)->SetGrid();
  ch->Draw("fracU235_LA4:Entry$-113>>hLA4(75,0,75,1000,0,0.94)","Entry$>=113&&Entry$<163","goff");
  TH2D *hLA4 = (TH2D*)gDirectory->Get("hLA4");
  // hLA4->Draw();
  // hLA4->SetTitle("Fuel Fractions Daya Bay Reactor 2");
  // hLA4->GetYaxis()->SetTitle("Fuel Fractions");
  // hLA4->GetXaxis()->SetTitle("Weeks Since Refuel");
  hLA4->SetMarkerStyle(kPlus);
  hLA4->SetMarkerColor(kRed);
  ch->Draw("fracU238_LA4:Entry$-113>>hLA4_1(75,0,75,1000,0,0.94)","Entry$>=113&&Entry$<163", "goff");
  TH2D *hLA4_1 = (TH2D*)gDirectory->Get("hLA4_1");
  hLA4_1->SetMarkerStyle(kPlus);
  hLA4_1->SetMarkerColor(kBlue);
  ch->Draw("fracPu239_LA4:Entry$-113>>hLA4_2(75,0,75,1000,0,0.94)","Entry$>=113&&Entry$<163","goff");
  TH2D *hLA4_2 = (TH2D*)gDirectory->Get("hLA4_2");
  hLA4_2->SetMarkerStyle(kPlus);
  hLA4_2->SetMarkerColor(kMagenta);
  ch->Draw("fracPu241_LA4:Entry$-113>>hLA4_3(75,0,75,1000,0,0.94)","Entry$>=113&&Entry$<163","goff");
  TH2D *hLA4_3 = (TH2D*)gDirectory->Get("hLA4_3");
  hLA4_3->SetMarkerStyle(kPlus);
  hLA4_3->SetMarkerColor(kBlack);
  hLA4->Draw("same");
  hLA4_1->Draw("same");
  hLA4_2->Draw("same");
  hLA4_3->Draw("same");

  gPad->Update();
  TPaveText *ptLA3[8];
  TString nameLA3[8] = {"LA3  ^{235}U: #times", "LA3   ^{238}U: #times","LA3 ^{239}Pu: #times","LA3 ^{241}Pu: #times","LA4 ^{235}U: +", "LA4 ^{238}U: +","LA4 ^{239}Pu: +","LA4 ^{241}Pu: +"};
  for(int i=0;i<8;++i){
    ptLA3[i]= new TPaveText(0.83,0.83-i*size,0.99,0.89-i*size,"ndc");
    ptLA3[i]->SetFillColor(0);
    ptLA3[i]->SetShadowColor(0);
    ptLA3[i]->SetBorderSize(0);
    ptLA3[i]->SetTextColor(col[i%4]);
    ptLA3[i]->AddText(nameLA3[i].Data());
    ptLA3[i]->Draw();
  }

}
