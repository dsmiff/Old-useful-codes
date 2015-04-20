{
  gROOT->Reset();
  gStyle->SetOptStat(1111);
  // Get Mbb histogram from file
  TFile f("FILENAME");
  TH1F *Mbb = (TH1F*)f.Get("HISTNAME");
  Mbb->Draw();

  // Histo parameters
  int N = Mbb->GetEntries();
  int bins = Mbb->GetSize();
  bins -= 2;
  cout << "Number of bins: " << bins << endl;

  TCanvas c("c1", "CANVASNAME", 10, 10, 600, 600); 
  c1->cd();

  // Fs(x) = (1/N)*\sum_{j=1}^{m} n_j Gauss(x, t_j, h_j)
  // h_j = (4/3/N)^{1/5} \Delta x \sqrt{N/n_j}
  // n_j - number of events in bin j
  // t_j - center of bin j
  // h_j - smoothing parameter (sigma)
  Double_t fy[bins], fx[bins];
  double fs, hj, nj;
  double Dx = Mbb->GetBinWidth(1);
  double smf = 1.0; // smearing factor (added by Dan)
  for (int i=0; i<bins; i++){
    fx[i] = Mbb->GetBinCenter(i+1);
  }
  for (int i=0; i<bins; i++){
    fs = 0.0;
    for(int j=0; j<bins; j++){
       nj = Mbb->GetBinContent(j+1);
       hj = TMath::Power((4./3./N), 0.2)*Dx;
       if(smf<1) // this factor redistributes the weights
         hj *= TMath::Sqrt(N/nj);
       fs += nj*TMath::Gaus(fx[i], fx[j], smf*hj);
    }
    fy[i] = fs/N;
    cout << "SKE x=" << fx[i] << " Fs(x)=" << fy[i] << endl;
  }
  TGraph *gr = new TGraph(bins,fx,fy);

  // Scale fit for visual comparison
  if(smf==1)
    Mbb->Scale(1./Mbb->Integral());
  else	// gr is not normalised anymore
    Mbb->Scale(gr->Integral()/Dx/Mbb->Integral());
  Mbb->SetTitle("HISTTITLE");
  Mbb->SetXTitle("Mass");
  Mbb->Draw();

  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->Draw("CSAME");

  TLegend *legend = new TLegend(0.6,0.5,0.8,0.7);
  legend->AddEntry(Mbb, "Data");
  legend->AddEntry(gr, "SKE");
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->Draw();
  
  gPad->Update();
  gPad.Print("FILEOUTPUT");
}
