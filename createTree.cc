#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <math.h>       
#include <cmath>        
#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "TVector3.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TBrowser.h"
#include "TApplication.h"
#include <vector>
#include "TRandom.h"
#include <sstream>
#include <algorithm>
using namespace std;

vector<double> Energy;
vector<double> dEdXBins;

double findWET(double, double);

int main(int argc, char *argv[]) {

  int detectorVersion = 1; //*** 1 for ideal trackers, 2 for Bergen DTC
  
  //Variable definitions			
  TVector3 TPSsource(0,-6600,0); //*** Position of scanning magnets [mm] and where proton angles are calculated from. Every position coming from GATE is in mm
  float D_t = 50.0; //*** mm between trackers
  float TPSx0, TPSy0, TPSz0, TPSpx0, TPSpy0, TPSpz0, 
        x11,y11,z11, x12,y12,z12, p1x,p1y,p1z, 
        x21,y21,z21, x22,y22,z22, p2x,p2y,p2z, 
        Einit,Estop;
  float WEPLideal, WEPLreal;      
  float tempX, tempY, tempZ;
  int IDenter11, IDenter12, IDexit21, IDexit22, IDevent; 
  vector<float> SpotX,SpotZ;
  float s_pos;
  float scatter_in, scatter_out;
  float matThick1,matThick2,matBudget1,matBudget2, X01,X02; 
  
  //-----------------Plan description with beam spot locations for TPS----------------------
  int NbProtonsperSpot;
  std::string linep;
  std::ifstream TPSGen("../gate/PlanDescriptionToGate.txt"); //*** TPS-plan file
  double datap[3];
  int j=0;
  while(getline(TPSGen, linep)) {
	if(j>=40){
      stringstream ss(linep);
      for(int i=0;i<3;i++) ss >> datap[i];
      SpotX.push_back(datap[0]);
      SpotZ.push_back(datap[1]);
    }
    if (j==40) NbProtonsperSpot = datap[2]; //beam spot positions and weights start at line 40!
    j++;
  }  
  cout<<NbProtonsperSpot<<endl;
    
  //----------------------------
  //For handling the proton WEPL. These .dat files must be prepared beforehand!
  //----------------------------
  //For getting the WEPL from the MC calculated energy difference
  std::string line;
  std::ifstream SPWater("Water_Proton.dat"); //file containing St.Pow. from underlying simulation (geant4 is made with 78eV!)
  double data[2];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<2;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }
  //For getting the systematic WEPL error of the detector
  vector<float> WEPLsys;
  vector<float> WEPLsysErr;
  std::string line1;
  std::ifstream errWEPL("weplErr.dat"); //file containing systematic WEPL error from DTC characterization
  double data1[2];
  while(getline(errWEPL, line1)) {
    stringstream ss(line1);
    for(int i=0;i<2;i++) ss >> data1[i];
    WEPLsys.push_back(data1[0]);
    WEPLsysErr.push_back(data1[1]);
  }  
  //For getting the sigma of the detector estimated WEPL, range straggling (noise)
  vector<float> WEPLp;
  vector<float> sigmap;
  std::string line2;
  std::ifstream sigmaWEPL("weplSigma.dat"); //file containing range straggling from DTC characterization
  double data2[2];
  while(getline(sigmaWEPL, line2)) {
    stringstream ss(line2);
    for(int i=0;i<2;i++) ss >> data2[i];
    WEPLp.push_back(data2[0]);
    sigmap.push_back(data2[1]);
  }   
  
//----------------------------
//Prepare and combine the root files from simulations/phasespace detectors
//---------------------------- 
  char enter1f[50] = "../output/entering11Ideal_"; //***
  char enter2f[50] = "../output/entering12Ideal_"; //***
  char exit1f[50] = "../output/exiting21Ideal_"; //***
  char exit2f[50] = "../output/exiting22Ideal_"; //***
  
  char outf[50] = "../outputFiles/yourFileName_"; //***
  
  strcat(enter1f,argv[1]);
  strcat(enter2f,argv[1]);
  strcat(exit1f,argv[1]);
  strcat(exit2f,argv[1]);
  
  strcat(outf,argv[1]);
  char rootf[10] = ".root";

  strcat(enter1f,rootf);
  strcat(enter2f,rootf); 
  strcat(exit1f,rootf);
  strcat(exit2f,rootf);
  
  strcat(outf,rootf);

  TFile *entering11 = new TFile(enter1f,"update"); 
  TFile *entering12 = new TFile(enter2f,"update"); 
  TFile *exiting21 = new TFile(exit1f,"update"); 
  TFile *exiting22 = new TFile(exit2f,"update"); 
  
  TFile* fileOut = new TFile(outf,"RECREATE"); 

  TTree *t11 = (TTree*)entering11->Get("PhaseSpace"); 
  TTree *t12 = (TTree*)entering12->Get("PhaseSpace"); 
  TTree *t21 = (TTree*)exiting21->Get("PhaseSpace");
  TTree *t22 = (TTree*)exiting22->Get("PhaseSpace"); 
  
  TTree *tOut = new TTree("OutTree", "Collection of tracker variables"); //***SET MANUALLY*** name of tree in the output root file 
  
//----------------------Read data from Monte Carlo made root files--------------------------
  //Front tracker data  
  t11->SetBranchAddress("X",&x11); 
  t11->SetBranchAddress("Y",&y11);
  t11->SetBranchAddress("Z",&z11);
  t11->SetBranchAddress("EventID",&IDenter11);
  
  t12->SetBranchAddress("X",&x12);
  t12->SetBranchAddress("Y",&y12);
  t12->SetBranchAddress("Z",&z12);
  t12->SetBranchAddress("EventID",&IDenter12);
  
  t12->SetBranchAddress("Ekine", &Einit);

  //Rear tracker data:
  t21->SetBranchAddress("Ekine", &Estop);
  
  //Tracker21
  t21->SetBranchAddress("X",&x21);
  t21->SetBranchAddress("Y",&y21);
  t21->SetBranchAddress("Z",&z21);
  t21->SetBranchAddress("EventID",&IDexit21);

  //Tracker22
  t22->SetBranchAddress("X",&x22);
  t22->SetBranchAddress("Y",&y22);
  t22->SetBranchAddress("Z",&z22);
  t22->SetBranchAddress("EventID",&IDexit22);

//----------------------Compile new root file from Monte Carlo root files-----------------
  //Beam data:
  tOut->Branch("Einit",&Einit);  
  tOut->Branch("Estop",&Estop);
  tOut->Branch("scatter_in",&scatter_in);
  tOut->Branch("scatter_out",&scatter_out);
  tOut->Branch("eventID",&IDexit22);  
  tOut->Branch("WEPLideal",&WEPLideal);
  tOut->Branch("WEPLreal",&WEPLreal);
    
  //TPS data:    
  tOut->Branch("TPSx0",&TPSx0);
  tOut->Branch("TPSy0",&TPSy0);  
  tOut->Branch("TPSz0",&TPSz0);
  
  tOut->Branch("TPSpx0",&TPSpx0);
  tOut->Branch("TPSpy0",&TPSpy0);
  tOut->Branch("TPSpz0",&TPSpz0);
  
  //Front tracker data
  tOut->Branch("x11",&x11);
  tOut->Branch("y11",&y11);
  tOut->Branch("z11",&z11);
  
  tOut->Branch("x12",&x12);
  tOut->Branch("y12",&y12);
  tOut->Branch("z12",&z12);
  
  tOut->Branch("p1x",&p1x);
  tOut->Branch("p1y",&p1y);
  tOut->Branch("p1z",&p1z);
  
  //Rear tracker data
  tOut->Branch("x21",&x21); 
  tOut->Branch("y21",&y21);
  tOut->Branch("z21",&z21);
  
  tOut->Branch("x22",&x22);
  tOut->Branch("y22",&y22);
  tOut->Branch("z22",&z22);

  tOut->Branch("p2x",&p2x);
  tOut->Branch("p2y",&p2y);
  tOut->Branch("p2z",&p2z);  
  
  int Nprotons = t22->GetEntriesFast(); //Needs to be t22(exiting protons) otherwise the following for loop tries to find protons that has been lost 
  int enterID11 = 0; int enterID12 = 0; 
  int exitID21 = 0; int exitID22 = 0;
  int ProcessedSpots = 0;
  
  cout<<Nprotons<<endl;
    
  for(int i=0;i<Nprotons-100;i++){ //Loop over each and every proton and get their information
	//GetEntries for initializing!!
	  t11->GetEntry(enterID11); 
	  t12->GetEntry(enterID12);	  
    t21->GetEntry(exitID21);
    t22->GetEntry(exitID22); 

    //We align everything to the protons that actually made it to the detector. Meaning highest eventID, which is "IDexit22"     
	  while (IDexit22!=IDenter11 || IDexit22!=IDenter12 || IDexit22!=IDexit21){
      while(IDexit22<IDenter11 || IDexit22<IDenter12 || IDexit22<IDexit21){
	    exitID22++;
	    t22->GetEntry(exitID22); 
	    }			
	    while(IDexit22>IDenter11){
	      enterID11++;
	      t11->GetEntry(enterID11);
        }
        while(IDexit22>IDenter12){
		    enterID12++;
	      t12->GetEntry(enterID12);
	    }
	    while(IDexit22>IDexit21){
		    exitID21++;
		    t21->GetEntry(exitID21);
	    }
    }   
  //The proton is identified	  
	if(IDexit22==IDenter11 && IDexit22==IDenter12 && IDexit22==IDexit21) {

	  t11->GetEntry(enterID11);
	  enterID11++;
	  t12->GetEntry(enterID12);
	  enterID12++;
	  t21->GetEntry(exitID21);
	  exitID21++;
	  t22->GetEntry(exitID22);
	  exitID22++;
	  
	  //Make sure we get the correct TPS spots from the plan description file
	  if(IDexit22 >= NbProtonsperSpot*ProcessedSpots){   
		  if(IDexit22 < NbProtonsperSpot*(ProcessedSpots+1)){
			  TPSx0 = SpotX.at(ProcessedSpots);
			  TPSz0 = SpotZ.at(ProcessedSpots);
			  TPSy0 = 0.0; //mm ***SET MANUALLY*** The longitidunal position of the isocenter
			  TVector3 TPSspot(TPSx0,TPSy0,TPSz0);
		  	TVector3 dir = TPSspot - TPSsource;
		  	dir.SetMag(1); //Normalizes the dir vector
			
			  TPSpx0 = dir.x();
			  TPSpy0 = dir.y();
			  TPSpz0 = dir.z();
			
		   	//Project to tracker position, at -100 mm 
			  TPSx0 = TPSx0 - ((TMath::Tan(TPSpx0))*100);
			  TPSz0 = TPSz0 - ((TMath::Tan(TPSpz0))*100);
			  TPSy0 = -100.0;
		  } else {ProcessedSpots++;}
    }
      
    WEPLideal=findWET(Einit, Estop);  //Finding the ideal WEPL
    if(WEPLideal>250.0)continue; //only have data up to 250 mm WEPL
      
    if(detectorVersion==1){ //ideal detector
		  WEPLreal = WEPLideal;   
	    scatter_in = 0.0;
	    scatter_out = 0.0;
	  
	    p1x = atan((x12-x11)/D_t);
      p1y = 1.0; //Forward facing
      p1z = atan((z12-z11)/D_t);  
    
      p2x = atan((x22-x21)/D_t);
      p2y = 1.0; //Forward facing
      p2z = atan((z22-z21)/D_t);
	  } 
      
    if(detectorVersion==2){ // Bergen detector
      //Bergen properties
	    s_pos =0.005; //mm
		  //Prepare for scattering calculation using Highland   
		  matThick1 = 0.02; //cm of Carbon fiber
      X01 = 25.723; //cm Carbon
		  matBudget1 = matThick1/X01;
		  matThick2 = 0.002; //cm of Silicon
		  X02 = 9.365; //cm Silicon
		  matBudget2 = matThick2/X02;
		  //Calculate Sigma_scatter_in and out      
		  scatter_in = 0.0; //no front tracker
	
	    //Calculate Sigma_scatter_out using Highland  
	    double E = Estop;   
	    double tau = E/938.3;
	    double pv = ((tau+2)/(tau+1))*E;
	    double scat_out1 = (14.1/(pv)*sqrt(matBudget1)*(1+((1/9)*log10(matBudget1))));
	    double scat_out2 = (14.1/(pv)*sqrt(matBudget2)*(1+((1/9)*log10(matBudget2))));
	    scatter_out = sqrt((scat_out1*scat_out1)+(scat_out2*scat_out2));

      //model new WEPL!
      int it_weplErr = lower_bound(WEPLsys.begin(), WEPLsys.end(), WEPLideal) - WEPLsys.begin();
      float weplDectErr = WEPLsysErr.at(it_weplErr);
      
      float WEPLdect = WEPLideal+weplDectErr;

      int it_weplSigma = lower_bound(WEPLp.begin(), WEPLp.end(), WEPLdect) - WEPLp.begin();
      float funnet_sigma = sigmap.at(it_weplSigma); 
      
      WEPLreal = gRandom->Gaus(WEPLdect, funnet_sigma); //Finding the modelled WEPL
      	
	    //create the direction angles and uncertain positions using Gaussian based on position resolution from ALPIDE      
	    double newXenter1 = gRandom->Gaus(x11, s_pos);
	    double newZenter1 = gRandom->Gaus(z11, s_pos);      
	    double newXenter2 = gRandom->Gaus(x12, s_pos);    
	    double newZenter2 = gRandom->Gaus(z12, s_pos);     
	    double newXexit1 = gRandom->Gaus(x21, s_pos);
	    double newZexit1 = gRandom->Gaus(z21, s_pos);     
	    double newXexit2 = gRandom->Gaus(x22, s_pos);   
	    double newZexit2 = gRandom->Gaus(z22, s_pos);
	
	    p1x = atan((newXenter2-newXenter1)/D_t);
	    p1y = 1.0; //Forward facing
      p1z = atan((newZenter2-newZenter1)/D_t);         
		
	    p2x = atan((newXexit2-newXexit1)/D_t);
      p2y = 1.0; //Forward facing
	    p2z = atan((newZexit2-newZexit1)/D_t);  
				
	    x11 = newXenter1;
	    x12 = newXenter2;
	    z11 = newZenter1;
	    z12 = newZenter2;
	    x21 = newXexit1;
	    x22 = newXexit2;
	    z21 = newZexit1;
	    z22 = newZexit2;     
    }
	          
    tOut->Fill(); //fill the tree with the data from latest GetEntries!
    if (i%100000==0) cout << i << endl;
    }

  }
  cout<<"done"<<endl;
  
  entering11->Close();
  entering12->Close();
  exiting21->Close();
  exiting22->Close();
  
  fileOut->cd(); //set which file to write to, by default it is the last open file
  tOut->Write("", TObject::kOverwrite);
  fileOut->Close();
  
  return 0; 
}

////////////////////////////////////////////
// Extract WET
////////////////////////////////////////////
double findWET(double Einit,double Estop){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
  int it_Estop = lower_bound(Energy.begin(), Energy.end(), Estop)-Energy.begin();
  double step = Energy.at(1)-Energy.at(0);
  double WET = 0 ;
  for(int i=it_Estop;i<it_Einit;i++){
    WET += 1./dEdXBins[i];
  }
  return WET*step; // to get WET in cm, divide this by 10;
}
