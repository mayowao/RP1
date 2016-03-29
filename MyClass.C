#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <strstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <list>
#include <map>
#include <vector>
#include <iterator>

using namespace std;

//Create struct
   struct evaluatejet {
     int number;
     float phip;
     float etap;
     float ptp;
     float massp; 
   }particles; 

//Add lists
    list<evaluatejet>jet;

//variables being used
   float phii,phij,etai,etaj,pti,ptj,massi,massj;
   Double_t R,dij,diB,djB,mind, combp, combpt, thetai, thetaj;
   Double_t combeta,combtheta,combpx,combpy,combpz,ppi,ppj;

//Creating the structure with starting set of particles from the original vector 
   for (Long64_t k=0;k<nentries;k++) {
     jet[k].number = k;
     jet[k].phip = phi->GetEntry(k);
     jet[k].etap = eta->GetEntry(k);
     jet[k].ptp = pt->GetEntry(k);
     jet[k].massp = mass->GetEntry(k);
   }



void MyClass::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MyClass.C
//      Root > MyClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   TBranch *jetp = tree->Branch("pt",&pt,"pt/F");
   //tree->SetBranchAddress("jetptp",&pt);
   //tree->SetBranchAddress("jetetap",&eta);
   //tree->SetBranchAddress("jetphip",&eta);
   //tree->SetBranchAddress("jetmassp",&eta);
   TH1F *jetphip = new TH1F("jetphip","phi",-4,0.,4);
   TH1F *jeteta = new TH1F("jeteta","eta", -4,0.,4);
   TH1F *jetptp = new TH1F("jetptp","pt",0,0.,500);
   TH1F *jetp = new TH1F("jetp","p",0,0,500);
   TH1F *jetmassp = new TH1F("jetmassp","mass",0,0.,500);

   //Need new tree so that values can be deleted/altered and clone the data into these friends??

   int stoploop =nentries;
   //   while (jetempty==false){//only evaluate whe our original strcut is nonempy??
   for (Long64_t i=0;i<stoploop;i++){
     for (Long64_t j=i+1; j<stoploop;j++){//start j=i accounts for all particles??
             phii= jet[i].phip;
	     phij= jet[j].phip;
	     etai= jet[i].etap;
	     etaj= jet[j].etap;
	     pti= jet[i].ptp;
	     ptj= jet[j].ptp;
	     massi=jet[i].massp;
	     massj=jet[j].massp;

	     R=TMath::Sqrt((phii-phij)^2+(etai-etaj)^2);
	    

	     dij=TMath::(R^2)*min(pti^(-2),ptj^(-2));
	     diB=TMath::(pti)^(-2);
	     djB=TMath::(ptj)^(-2);
             mind=TMath::Min(dij,diB,djB);

      //if the min is i with beam:
	     if  (mind==diB){
	     jetp->Fill(i));
	     jetphip->Fill(jet[i].phip);
	     jetptp->Fill(jet[i].ptp);
             jetetap->Fill(jet[i].etap);
	     jetmassp->Fill(jet[i].massp);

     	     //Get rid of particle from list
             
     }

      //if the min is j with beam
     if  (mind==djB){
	     jetp->Fill(j));
             jetphip->Fill(jet[j].phip);
             jetptp->Fill(jet[j].ptp);
             jetetap->Fill(jet[j].etap);
	     jetmassp->Fill(jet[j].massp);    

              //Get rid of particle from struct
	     for (Long64_t k=j;k<stoploop;k++){
	       jet[k].phip=jet[k+1].phip;
	       jet[k].etap=jet[k+1].etap;
	       jet[k].ptp=jet[k+1].ptp;
	       jet[k].massp=jet[k+1].massp;
	     }
	     stoploop = stoploop-1;
   
	     // jet[j].phip=[];
	     // jet[j].etap=[];
	     // jet[j].ptp=[];
   }

      //if the min is distance between particles
	  if  (mind==dij){
	    //make the ith particle the combined particle
	    thetai=TMath::2*(ATan(E(-etai)));
	    thetaj=TMath::2*(ATan(E(-etaj)));
	    ppi=TMath::pti/Sin(thetai);
	    ppj=TMath::ptj/Sin(thetaj);
	    combp=TMath::Sqrt((pti*Cos(phii)+ptj*Cos(phij))^2+(pti*Sin(phii)+ptj*Sin(phij))^2+(ppi*Cos(thetai)+ppj*Cos(thetaj))^2);
	    combpx=TMath::ppi*Sin(thetai)*Cos(phii)+ptj*Sin(thetaj)*Cos(phij);
	    combpy=TMath::ppi*Sin(thetai)*Cos(phii)+ptj*Sin(thetaj)*Sin(phij);

	    combpt=TMath::Sqrt(combpx^2+combpy^2);
	    combtheta=TMath::ASin(combpt/combp);
	    combeta=TMath::-LogE(Tan(combtheta/2);
	    jet[i].ptp=combpt;
	    jet[i].massp=TMath::massi+massj;
	    jet[i].etap=combeta;
	    jet[i].phip=TMath::ASin(combpx/combpt);

	 
	  }//for j 

}//fori

//create histograms

}//fornbytes 
}//void
