{
  gROOT->Reset();
  gStyle->SetOptStat(1111);
  
  //Get the histograms
  TFile f("ExampleDSTAnalysis.root");
  TH1F *h1 = (TH1F*)f.Get("fFTVNNbb");
  TH1F *h2 = (TH1F*)f.Get("fFTVNN4all");
  TH1F *h3 = (TH1F*)f.Get("fFTVNNbb");
  TCanvas *c1 = new TCanvas("c1","B-Tagging Purity Vs. Efficiency",200,10,700,500);

   c1->SetGrid();

 //GRAPH
  int bins = h1->GetSize();
  bins -= 2;
  Double_t fy[bins], fx[bins];
  float entries1=0.0, integral1=0.0, fint1=0.0, integral2=0.0, integral3=0.0, fint2=0.0;
  //efficiency
  for (int i=0; i<bins; i++){
    entries1 = h1->GetEntries();
    integral1 = h1->Integral(i,bins);
    fint1 = integral1/entries1;
    fx[i] = fint1;
  }
  //purity
  for (int i=0; i<bins; i++){
    //number of events that passed the cut
    integral2 = h2->Integral(i,bins);
    //number of events that passed the cut and were true b-jets
    integral3 = h3->Integral(i,bins);
    fint2 = integral3/integral2;
    fy[i] = fint2;
  }

  for (int i=0; i<bins; i++){
    cout << "fx["<<i<<"] = " << fx[i] << "    fy["<<i<<"] = " << fy[i] << " bin num " << i << " of total number of bins " << bins << endl;
}
   const Int_t n = 100;
   auto gr = new TGraph(n,fx,fy);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("B-Tagging Purity Vs. Efficiency");
   gr->GetXaxis()->SetTitle("Efficiency");
   gr->GetYaxis()->SetTitle("Purity");
   gr->Draw("ACP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
//   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();

  gPad->Update();
//  gPad.Print("~/Desktop/b-tagging_purity_vs_eff_WWsig_10000.png");

}

