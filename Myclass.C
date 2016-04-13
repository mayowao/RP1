#define Myclass_cxx
#include "Myclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <string>
#include <strstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <list>
#include <vector>
#include <iterator>

//Denotes namespace
using std::vector;
using std::cout;
using std::list;

//Create index variables for vector                                                                  
unsigned int it, jt; 
int minindex_i, minindex_j, minindex_iB;

//Create Instance Variables                                                                          
float px, py, pz, newphi, neweta, newmass, newtheta, totalpm, newpT, Dib, Dij, Rij, Rrr, deltaphi, deltaeta, Dijmin, Dibmin;

//Set Speed of light for Energy Calc                                                    
const float csq = pow(299792458, 2);

//Set Pi for Phi Wrap                                                                   
const double PI = 4 * atan(1);
//cout::<<"PI is "<< PI;

//Function for Phi Wrap Soln                                                                         
float PhiWrap( float val )
{
  float output;
  if( val > PI)
    output = ((2 * PI) - val);
  else
    output = val;
  return output;
}

//Function to find distance between particles                                           
float Distance_P (float aphi, float bphi, float ceta, float deta, float ept, float fpt)
{
  deltaphi = PhiWrap (aphi - bphi) ;
  deltaeta = ceta - deta;
  Rij = hypot( deltaphi, deltaeta);
  Dij = std::min( pow( ept, -2), pow( fpt, -2));
  Dij = Dij * pow( ( Rij / Rrr), 2);
  return Dij;
}

//Function for Beam Distance                                                            
float Distance_J (float gpt)
{
  Dib = pow( gpt, -2);
  return Dib;
}

//Declare a struct (AllParticles) to hold arrays which will hold all of the values. The object particle of type P can acess all of these values                                           
struct AllParticles
{
  float pt;
  float phi;
  float eta;
  float mass;
  float energy;
  float widestphi;
  int wphifirst;
  int wphisecond;
};

//Create a struct object that can access the member variables
AllParticles item;

//Create Vectors of Struct type
vector <AllParticles> Particles;
vector <AllParticles> Jets;
vector <AllParticles> HighEnergyJets;
vector <AllParticles> OriginalParticles;

//List fo Find Highest pT Jet
list <float> SortingHEJets;

//Very Important Method
void Myclass::Loop()
{

//   In a ROOT session, you can do:
//      Root > .L Myclass.C
//      Root > Myclass t
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

  //Create Canvas
  TCanvas *c1 = new TCanvas("c1", "Jet Clustering", 200, 10, 1400, 700);
  c1 -> SetFillColor(42);
  c1 -> Divide(2, 2);
  TPad *p1 = (TPad *)(c1 -> cd(1));
  p1 -> SetLogy();
  
  //Set up Histograms
  
   //pT of jets 
  TH1F *histo3 = new TH1F("histo3", "pT Spectrum Jets over 50 GeV", 100, 0, 2500);
  histo3 -> SetMarkerStyle(4);
 
  //Number of jets in an event
  TH1F* histo2 = new TH1F("histo2", "Number of Jets", 100, 0, 55);
  histo2 -> SetMarkerStyle(4);
  
  //pT of High Energy Jets
  TH1F *histo1 = new TH1F("histo1", "pT Spectrum  Jets", 100, -200, 1300);
  histo1 -> SetMarkerStyle(4);

  //pT of highest Energy Jet
  TH1F *histo4 = new TH1F("histo4", "pT of Highest Energy Jet per Event", 100, 0, 3000);
  histo4 -> SetMarkerStyle(4);

  //Eta plot
  TH1F *histo5 = new TH1F("histo5", "Eta Plot", 100, -10, 10);
  histo5 -> SetMarkerStyle(4);

  //Phi plot
  TH1F *histo6 = new TH1F("histo5", "Phi Plot", 100, -3.5, 3.5);
  histo6 -> SetMarkerStyle(4);

  //Set Cut for pT HEJets
  float cutpT = 50;

  //Fail-Safe
  if (fChain == 0) 
    return;

  //This value can be changed easily and it denotes the size of the cones          
  Rrr = .4;

  //Get the numner of events
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //EVENT LOOP
  for (Long64_t jentry=0; jentry < nentries; jentry++) 
    {
      //Load the event in to memory
      Long64_t ientry = LoadTree(jentry);

      //Print the event number
      std::cout<< "jentry = " << jentry << std::endl;
      
      //Fail-safe
      if (ientry < 0) 
	break;
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;
      
      // if (Cut(ientry) < 0) continue;
      
      //Get the number of elements in one of the vectors
      int numentry = pt->size();
      
      //Fill the array "particles" with the values from the vectors                             
      //Check pointer to array syntax for accesing element                                      
      for (int wh = 0; wh < numentry - 1; wh++)
        {
          item.pt = (*pt)[wh];
          item.phi = (*phi)[wh];
          item.eta = (*eta)[wh];
          item.mass = (*mass)[wh];
          float thetas = 2*atan(exp(-(*eta)[wh]));
          float calcpx = ((*pt)[wh])*cos((*phi)[wh]);
          float calcpy = ((*pt)[wh])*sin((*phi)[wh]);
          float calcpz = ((*pt)[wh])*cos(thetas);
          float calcp = sqrt(pow(calcpx,2)+pow(calcpy,2)+pow(calcpz,2));
          item.energy = sqrt(pow(calcp,2)+pow((*mass)[wh],2));
          item.widestphi = 0;
          item.wphifirst = wh;
          item.wphisecond = wh;
          Particles.push_back( item );
          OriginalParticles.push_back( item );
        }

      
      //This is to make sure that the first value checked is assigned to the smallest           
      Dijmin = 100000;
      Dibmin = 100000;

      //Ensure that the loop keeps going until every particle is in a beam
      while ( !Particles.empty() )     
	{
	  //Loop over the first particle in the pair
	  for( it = 0; it < Particles.size(); it++)
	    { 
	      //Loop over the second particle in the pair
	      for (jt = it + 1; jt < Particles.size(); jt++)
		{
		  //Calculate the distance between the particles
		  Dij = Distance_P( Particles[it].phi, Particles[jt].phi, Particles[it].eta, Particles[jt].eta, Particles[it].pt, Particles[jt].pt );
		  //Determine if this is the smallest Dij so far
		  if (Dij < Dijmin)                 
		    {
		      Dijmin = Dij;
		      minindex_i = it;//Save indeces for later
		      minindex_j = jt;
		    }
		  
		}//Exit inner (jt) for loop

	      //Calculate beam distance
	      Dib = Distance_J( Particles[it].pt );
	      //Figure out if the Dib is the smallest so far
	      if (Dib < Dibmin)
		{
		  Dibmin = Dib;
		  minindex_iB = it; //Save the particle so it can be removed later
		}//if (Dib < Dibmin)

	    }// exit outer (it) for loop

      //If the smallest is a beam, add to beam list and remove from particle list          
      if (Dibmin < Dijmin)
	{
	  Jets.push_back( Particles[minindex_iB] );
	  Particles.erase( Particles.begin() + minindex_iB );
	 
	  //Cut for HEJets
	  if ( Jets.back().pt > cutpT )
	    {
	      HighEnergyJets.push_back( Jets.back() );
	      SortingHEJets.push_back( Jets.back().pt );
	      //Book a histogram for the pT of high energy jets                                                      
	      histo3 -> Fill( HighEnergyJets.back().pt );
	    }
	  //Book a histogram for the pT of each jet
	  histo2 -> Fill( Jets.back().pt );
          histo5 -> Fill( Jets.back().eta );
          histo6 -> Fill( Jets.back().phi );

	}

      //If the smallest is not a beam, add momenta, remove other two particles, and add to vector
      else
	{
	  //ADDING MOMENTUM                                                                                                     
	  //Adds the components of momentum in the x, y, and z directions     
	  px = ( ( (Particles[minindex_i].pt * cos((Particles[minindex_i].phi) ) ) + (Particles[minindex_j].pt * cos(Particles[minindex_j].phi) ) ) );
	  py = ( ( (Particles[minindex_i].pt * sin((Particles[minindex_i].phi) ) ) + (Particles[minindex_j].pt * sin( Particles[minindex_j].phi) ) ) );
	  pz = ( ( (Particles[minindex_i].pt * sinh(Particles[minindex_i].eta) ) ) + (Particles[minindex_j].pt * sinh(Particles[minindex_j].eta) ) );

	  //Calculate new values for new paricle         
	  newpT = hypot(px, py);
	  totalpm = sqrt((px * px) + (py * py) + (pz * pz));
	  newtheta = atan2( newpT , pz ); //use atan2
	  neweta = -log( tan( newtheta/ 2));
	  newphi = atan2( py, px );
          
	  float newcalcpx = (newpT)*cos(newphi);
          float newcalcpy = (newpT)*sin(newphi);
          float newcalcpz = (newpT)*cos(newtheta);
          float newcalcp = sqrt(pow(newcalcpx,2)+pow(newcalcpy,2)+pow(newcalcpz,2));
          float newenergy = Particles[minindex_i].energy + Particles[minindex_j].energy;
	  //Add energies not mass                        
	  newmass = sqrt(pow(newenergy,2)-pow(newcalcp,2)); // Use energy to calculate new mass
	  
          //Firgure out the new widestphi
          //int wphifindex = Particles[
          //float phibetweenf = max(OriginalParticles[lastphii1].phi-OriginalParticles[lastphij1].phi,OriginalParticles[l;
	  float phibetween = 0;
          float newwidestphi = max(Particles[minindex_i].widestphi,max(Particles[minindex_j].widestphi,phibetween));

	  //Build New Struct with all new values                                          
	  item.pt = newpT;
	  item.eta = neweta;
	  item.phi = newphi;
          item.mass = newmass;
	  item.energy = newenergy;
	  item.widestphi = newwidestphi;
	  
	  //Remove the Particles
	  //Check to see which one is bigger before removing so as so not mess up indeces
	  if ( minindex_j < minindex_i )
	    {
	      Particles.erase( Particles.begin() + (minindex_i) );
	      Particles.erase( Particles.begin() + (minindex_j) );
	    }
	  else
	    {
	      Particles.erase( Particles.begin() + (minindex_j) );
	      Particles.erase( Particles.begin() + (minindex_i) );
	    }
	  //Add new Particle
	  Particles.push_back( item );
	  
	}//Exit ELSE statement for when the Dijmin is the smallest
      
      //Reset values for next iteration
      Dijmin = 100000;
      Dibmin = 100000;
      minindex_j = 0;
      minindex_i = 0;
      minindex_iB = 0;
	}//while loop
      
      //Fill Histogram with Highest pT Jet from each event
      SortingHEJets.sort();
      histo4 -> Fill( SortingHEJets.back() );
      
      //Fill histogram with the number of Jets per event
      histo1 -> Fill( Jets.size() );
      //histo5 -> Fill ( Jets.eta() );

      //Clear Jets vectors for next event
      SortingHEJets.clear();
      HighEnergyJets.clear();
      Jets.clear();
    }//Exit event (jentry) for loop
  
  //Histogram stuff
  histo3 -> Draw("");

  //c1 -> cd(2);
  TPad *p2 = (TPad *)(c1 -> cd(2));
  p2 -> SetLogy();
  histo2 -> Draw("");
 
  c1 -> cd(3);
  histo6 -> Draw("");
  c1 -> cd(4);
  histo5 -> Draw("");
  c1 -> SaveAs("prettypic.gif");
}//Void Loop()
