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
#include <stdio.h>      /* printf */
#include <math.h>

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
   float R,dij,diB,djB,mind, combp, combpt, thetai, thetaj;
   float combeta,combtheta,combpx,combpy,combpz,ppi,ppj;




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
     std::cout<<"We added"<<jet[k].phip;
   }
   
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
	     
	     R=sqrt(pow((phii-phij),2)+pow((etai-etaj),2));
	    
	     dij=(pow(R,2))*min(pow(pti,-2),pow(ptj,-2));
	     diB=pow(pti,-2);
	     djB=pow(ptj,-2);
             mind=min(min(dij,diB),djB);

             //if the min is i with beam:
	     if  (mind==diB){
	        jetp->Fill(i);
	        jetphip->Fill(jet[i].phip);
	        jetptp->Fill(jet[i].ptp);
                jetetap->Fill(jet[i].etap);
	        jetmassp->Fill(jet[i].massp);
 
     	     //Get rid of particle from list
                jet.erase(jet.begin()+i); 
             }

      //if the min is j with beam
             if  (mind==djB){
	        jetp->Fill(j);
                jetphip->Fill(jet[j].phip);
                jetptp->Fill(jet[j].ptp);
                jetetap->Fill(jet[j].etap);
	        jetmassp->Fill(jet[j].massp);    

              //Get rid of particle from list
	        jet.erase(jet.begin()+j);
    	    
             }

      //if the min is distance between particles
	     if  (mind==dij){   //ERROR: from TMath, and Double_t(line 163) sol'n need c++ math???
	    //make the ith particle the combined particle
	      //calculate values for the ith entries
	       thetai=2*(atan(exp(-etai)));   
	       thetaj=2*(atan(exp(-etaj)));
	       ppi=pti/sin(thetai);   
	       ppj=ptj/sin(thetaj);
	       combp=sqrt(pow((pti*cos(phii)+ptj*cos(phij)),2)+pow((pti*sin(phii)+ptj*sin(phij)),2)+pow((ppi*cos(thetai)+ppj*cos(thetaj)),2))*(dij/R);
	       combpx=ppi*sin(thetai)*cos(phii)+ptj*sin(thetaj)*cos(phij);
	       combpy=ppi*sin(thetai)*cos(phii)+ptj*sin(thetaj)*sin(phij);

	       combpt=sqrt(pow(combpx,2)+pow(combpy,2));
	       combtheta=asin(combpt/combp);
	       combeta=-log(tan(combtheta/2));
	       
	       //std::cout<<"thetai="<<thetai;
	       //std::cout<<"thetaj="<<thetaj;

	       //add to the ith spot
	       particle.number = i;
               particle.phip = asin(combpx/combpt);   
               particle.etap = combeta;
               particle.ptp = combpt;
               particle.massp = massi+massj;
               //std::cout<<combpt;
               std::cout<<" i is "<<particle.number<<" ";
               //jet.insert(jet.begin()+i, evaluatejet& particle);
               jet.push_back(particle);
               jet.erase(jet.begin()+i);
	       
	       std::cout<<"We changed "<<jet[i].number;
	       
	       //delete the jth particle
	       jet.erase(jet.begin()+j);
	     } //if combined particle
	    /* */ 
	     //std::cout<<"Hello";
	  }//for j 
       }//fori
}//for while
//create histograms
/*jetp->Draw();
jetphip->Draw();
jetetap->Draw();
jetmassp->Draw();
*/


std::cout<<"Hi";
}//fornbytes 
}//void

