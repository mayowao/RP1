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

using std::vector;

//Create struct
   struct evaluatejet {
     int number;
     float phip;
     float etap;
     float ptp;
     float massp; 
   };
   
   evaluatejet particle; 

//Add lists
    vector<evaluatejet> jet; //jet is a vector of structs
    

//variables being used
   float phii,phij,etai,etaj,pti,ptj,massi,massj;
   Double_t R,dij,diB,djB,mind, combp, combpt, thetai, thetaj;
   Double_t combeta,combtheta,combpx,combpy,combpz,ppi,ppj;




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

   TH1F *jetphip = new TH1F("jetphip","phip",-4,0.,4);
   TH1F *jetetap = new TH1F("jeteta","etap", -4,0.,4);
   TH1F *jetptp = new TH1F("jetptp","ptp",0,0.,500);
   TH1F *jetp = new TH1F("jetp","p",0,0,500);
   TH1F *jetmassp = new TH1F("jetmassp","massp",0,0.,500);

   Long64_t nentries = fChain->GetEntriesFast();

//Loops through the EVENTS. And we want to find number of jets of a certain type in each event
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


   //IDEA?:Need new tree so that values can be deleted/altered and clone the data into these friends??

//Creating the list of the starting set of particles from the original vector 
int numentries = pt->size();

   for (Long64_t k=0;k<numentries;k++) {   
     particle.number = k;
     particle.phip = (*phi)[k];   
     particle.etap = (*eta)[k];
     particle.ptp = (*pt)[k];
     particle.massp = (*mass)[k];
     jet.push_back(particle);
     std::cout<<"We added"<<jet.number<<;
   }
   
/*

   int stoploop =numentries;
   
   while (!jet.empty()){//only evaluate while our list is nonempy??  
   for (Long64_t i=0;i<stoploop;i++){
     for (Long64_t j=i+1; j<stoploop;j++){//start j=i accounts for all particles??
             phii= jet[i].phip;   //ERROR: need to change to iterators??
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
                jet.erase(i); 
             }

      //if the min is j with beam
             if  (mind==djB){
	        jetp->Fill(j));
                jetphip->Fill(jet[j].phip);
                jetptp->Fill(jet[j].ptp);
                jetetap->Fill(jet[j].etap);
	        jetmassp->Fill(jet[j].massp);    

              //Get rid of particle from list
	        jet.erase(i);
    	    
             }

      //if the min is distance between particles
	     if  (mind==dij){   //ERROR: from TMath, and Double_t(line 163) sol'n need c++ math???
	    //make the ith particle the combined particle
	       thetai=TMath::2*(ATan(E(-etai)));   //ERROR: expected unqualified id before num constnt, expected ;??
	       thetaj=TMath::2*(ATan(E(-etaj)));
	       ppi=TMath::pti/Sin(thetai);   //ERROR : pti not a member of TMath
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

	 //delete the jth particle
	       jet.erase(j);
	     }
	  }//for j 
       }//fori
}//for while
//create histograms
jetp->Draw();
jetphip->Draw();
jetetap->Draw();
jetmassp->Draw();
*/

std::cout<<"Hi";
}//fornbytes 
}//void

