#include "TH1.h";
#include <valarray>;
#include <cmath>;
#include "TVector3.h"
#include "TLorentzVector.h"

void Analysis()
{
 // Load shared library
gSystem->Load("lib/libExRootAnalysis.so");
gSystem->Load("libPhysics");
using namespace std;

//use vectors as (E,px,py,pz,pt)

double positron[4], eneutrino[4], bquark[4], wboson[4], tquark[4], weight;
double antimuon[4], muneutrino[4];
Int_t numberOfParticles;
bool printout = 0;
double pi = 3.14159265358979312;
cout << "Macro to make histograms of Madgraph generated events" << endl;

// Create chain of root trees
TChain chain("LHEF");
chain.Add("results/ht_events.root");

// Create object of class ExRootTreeReader
ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
Long64_t numberOfEntries = treeReader->GetEntries();

// Open output file
output=new TFile("results/ht.root","RECREATE");

// Get pointers to branches used in this analysis
TClonesArray *branchEvent = treeReader->UseBranch("Event");
TClonesArray *branchParticle = treeReader->UseBranch("Particle");
Long64_t numberOfEntries = treeReader->GetEntries();
  // Book histograms
// Define necessary histograms



TH1 *h1 = new TH1F("h1", "Elep", 50,0.0,500);
TH1 *h2 = new TH1F("h2", "Ptlep", 50,0.0,300);
TH1 *h3 = new TH1F("h3", "phi", 50,0.0,2*pi);
TH1 *h4 = new TH1F("h4","theta",50,0.0,pi);
TH1 *h5 = new TH1F("h5","phi1",2,0.0,1.0);
TH1 *h6 = new TH1F("h6","theta1",2,0.0,1.0);

// Loop over all events

for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
  //for(Int_t entry = 0; entry < 1000; ++entry) {
  if( printout ) cout << "**************************************************" << endl;
  if( printout ) cout << "Analyzing Event " << entry << endl;
  treeReader->ReadEntry(entry);
  TRootLHEFEvent *eventInfo = (TRootLHEFEvent*) branchEvent->At(0);
  weight = eventInfo->Weight;

  numberOfParticles = eventInfo->Nparticles;
  if( printout ) cout << "Number of particles in Event = " << numberOfParticles << endl;
  for( Int_t pcount=0 ; pcount < numberOfParticles ; pcount++){
    TRootLHEFParticle *particleMG = (TRootLHEFParticle*) branchParticle->At(pcount);

    //cout << "Particle id = " << particleMG->PID << endl;
    if( particleMG->Status == 1 ) {
     //final state particle
    if( particleMG->PID == -11 ) {
     if( printout) cout << "found a positron" << endl;
        positron[0]=particleMG->E;
        positron[1]=particleMG->Px;
        positron[2]=particleMG->Py;
        positron[3]=particleMG->Pz;
       }
       elseif( particleMG->PID == 12 ) {
        if( printout)  cout << "found a electron-neutrino" << endl;
        eneutrino[0]=particleMG->E;
        eneutrino[1]=particleMG->Px;
        eneutrino[2]=particleMG->Py;
        eneutrino[3]=particleMG->Pz;
       }
￼
 }

 }




 }
elseif( particleMG->PID == 5 ) {
 if( printout ) cout << "found a b-quark" << endl;
  bquark[0]=particleMG->E;
 bquark[1]=particleMG->Px;
 bquark[2]=particleMG->Py;
 bquark[3]=particleMG->Pz;
elseif( particleMG->PID == -13 ) {
 if( printout ) cout << "found an antimuon" << endl;
  positron[0]=particleMG->E;
 positron[1]=particleMG->Px;
 positron[2]=particleMG->Py;
 positron[3]=particleMG->Pz;
elseif( particleMG->PID == 14 ) {
 if( printout ) cout << "found a muon neutrino" << endl;
 eneutrino[0]=particleMG->E;
 eneutrino[1]=particleMG->Px;
 eneutrino[2]=particleMG->Py;
 eneutrino[3]=particleMG->Pz;
elseif( particleMG->PID == -37 ) {
  if( printout ) cout << "found a charged higgs" << endl;
}
      }
      if( particleMG->Status == 2 ) {
         if( particleMG->PID == 6 ){
           if( printout ) cout << "Found a top quark " << endl;
           tquark[0]=particleMG->E;
           tquark[1]=particleMG->Px;
           tquark[2]=particleMG->Py;
           tquark[3]=particleMG->Pz;
} }
    }  // looping over particles in event
      // Now calculate various observables
  
  // Define dot product, modulus, azimuthal angle and polar angle
    double Elep;
    Elep=abs(positron[0]);
    double Ptlep;
    Ptlep=sqrt(pow(positron[1],2)+pow(positron[2],2));
    double Pttop;
    Pttop=sqrt(pow(tquark[1],2)+pow(tquark[2],2));
    double Ptlep3;
    Ptlep3=sqrt(pow(positron[1],2)+pow(positron[2],2)+pow(positron[3],2));
    double Pttop3;
    Pttop3=sqrt(pow(tquark[1],2)+pow(tquark[2],2)+pow(tquark[3],2));
    double modPP;
    modPP=Ptlep*Pttop;
    double modPP3;
    modPP3=Ptlep3*Pttop3;
    double dotPP;
    dotPP=(tquark[1]*positron[1])+(tquark[2]*positron[2]);
    double dotPP3;
    dotPP3=(tquark[1]*positron[1])+(tquark[2]*positron[2])+(tquark[3]*positron[3]);
    double phi;
    phi=acos(dotPP/modPP);
    double theta;
    theta=acos(dotPP3/modPP3);

// Define angle between lepton and y-axis for diff orientations
    double alpha;
if (positron[1]>0 && positron[2]>0){
! alpha=atan(positron[1]/positron[2]);
}
else if (positron[1]>0 && positron[2]<0){
! alpha=pi/2-atan(positron[2]/positron[1]);
}
else if (positron[1]<0 && positron[2]<0){
! alpha=atan(positron[1]/positron[2])+pi;
}
else if (positron[1]<0 && positron[2]>0){
! alpha=3*pi/2-atan(positron[2]/positron[1]); }
// Rotate x component of top quark by alpha when lepton is orientated to y-axis
double Pxtopnew;
  Pxtopnew=tquark[1]*cos(alpha)+tquark[2]*sin(alpha);

// Find angle between lepton and top quark
if (Pxtopnew<0){
￼! phi=2*pi-phi; }
double phi1;
if (phi<pi/2 || phi>3*pi/2){
    phi1=0.25;}
else phi1=0.75;
double theta1;
if (theta>pi/4){
! theta1=0.75;} else theta1=0.25;

  // Fill histograms
    h1->Fill(Elep,weight);
    h2->Fill(Ptlep,weight);
    h3->Fill(phi,weight);
    h4->Fill(theta,weight);
    h5->Fill(phi1,weight);
    h6->Fill(theta1,weight);

// Define asymmetry paramters and conditions
double sigma1a;
sigma1a=h5->GetBinContent(1);
double sigma2a;
sigma2a=h5->GetBinContent(2);
double Aa;
Aa=(sigma1a-sigma2a)/(sigma1a+sigma2a);
double sigma1p;
sigma1p=h6->GetBinContent(1);
double sigma2p;
sigma2p=h6->GetBinContent(2);
double Ap;
Ap=(sigma1p-sigma2p)/(sigma1p+sigma2p);
  } // looping over events
cout << "Aa" << endl;
cout << Aa << endl;
cout << "Ap" << endl;
cout << Ap << endl;
Double_t scale;
 if( 1 ) {
TCanvas *c1 = new TCanvas("c1","Test",1);
c1->SetFillColor(0);
c1->SetGrid();
c1->GetFrame()->SetBorderSize(20);
c1->SetFrameFillColor(0);
c1->Divide(3,2);
c1->cd(1);
h1->GetXaxis()->SetTitle("E_{lep}");
h1->GetYaxis()->SetTitle("1/#sigma d#sigma/d#E_{lep}");
scale = 1/h1->Integral("width");
h1->Scale(scale);
h1->Draw();
h2->GetXaxis()->SetTitle("Pt_{lep}");
h2->GetYaxis()->SetTitle("1/#sigma d#sigma/d#Pt_{lep}");
scale = 1/h2->Integral("width");
h2->Scale(scale);
h2->Draw();
h3->GetXaxis()->SetTitle("phi");
h3->GetYaxis()->SetTitle("1/#sigma d#sigma/d#phi");
scale = 1/h3->Integral("width");
h3->Scale(scale);
h3->Draw();
h4->GetXaxis()->SetTitle("theta");
h4->GetYaxis()->SetTitle("1/#sigma d#sigma/d#theta");
scale = 1/h4->Integral("width");
h4->Scale(scale);
h4->Draw();
h5->GetXaxis()->SetTitle("phi1");
h5->GetYaxis()->SetTitle("1/#sigma d#sigma/d#phi1");
scale = 1/h5->Integral("width");
h5->Scale(scale);
h5->Draw();
h6->GetXaxis()->SetTitle("theta1");
h6->GetYaxis()->SetTitle("1/#sigma d#sigma/d#theta1");
scale = 1/h6->Integral("width");
h6->Scale(scale);
h6->Draw();
 c1->Print("plots.eps");
 output->Write();
 output->Close();
} }
