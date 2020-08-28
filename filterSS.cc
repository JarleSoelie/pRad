#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
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
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>
#include <random>
using namespace std;

//All of these are in order of cm and so should all input to the ComputeMLP function!! (Alternatively, one can change the polyonomial factors into orders of mm)
//230MeV
#define azero   5.77059e-6
#define aone    2.74001e-7
#define atwo   -2.49026e-8
#define athree  4.63381e-9
#define afour  (-2.65153e-10)
#define afive   6.22291e-12

#define X_0 36.1

double d_source = 660.0; //*** cm Distance to scanning magnets/source
double s_angle = 0.003; //*** rad This should be the uncertainty on the TPS angle 

double s_pos_out = 0.0005; //*** cm Position uncertainty in tracking layers. 0.0 = ideal detector, 0.0005 = Bergen detector

double d_entry = 0.0; //*** cm OBS! If not using manual hull, change accordingly (=0.0 if manual hull is used or found)
double d_exit = 0.0; //*** cm
double d_T = 5.0; //*** cm Distance between trackers 

int Nsteps = 512; //Number of steps to take in MLP

struct Proton{
  
  Float_t TPSx0,TPSy0,TPSz0;
  Float_t TPSpx0,TPSpy0,TPSpz0;
  
  Float_t x12_h,y12_h,z12_h;
  Float_t p1x_h,p1y_h,p1z_h; 
   
  Float_t x21,y21,z21, x21_h,y21_h,z21_h; 
  Float_t p2x,p2y,p2z, p2x_h,p2y_h,p2z_h;
  
  Float_t scatter_out, scatter_out_h;
  
  Float_t weplReal, weplReal_h;
  
};

struct inPos{
  double X_opt, Z_opt, theta_X_opt, theta_Z_opt;
};

struct inPos ComputeOptIn(Proton *, float);
double Sigmat1(double);
double Sigmay1(double);
double Sigmaty1(double);
double Sigmat2(double,double);
double Sigmay2(double,double);
double Sigmaty2(double,double);

int main(int argc, char *argv[]){
  // Load Proton data
  Proton Point;
  
  float d_phantom = 150.0; //*** mm between inner tracker and phantom edge. This is used to reduce the hull calculation time at later steps
  
  //For defining the reconstruction area for filtering, 0.75x1.25 mm pixels in this case, but change to whatever you want
  int NbinsX = 272; //*** 
  int NbinsZ = 128; //*** 
  float  Xmin = -102; //*** mm 
  float  Xmax = 102; //*** mm 
  float  Zmin = -80; //*** mm 
  float  Zmax = 80; //*** mm 
  
  //Prepare for hull (Head-phantom)
  TFile* phantomFile = new TFile("HeadPhantom.root","update");
  TH3S* hull = new TH3S("hull", "Hounsfield Units", 1024,-90,102, 1024,-90,102, 128,-166.25,-6.25); //mm
  phantomFile->GetObject("hu",hull); 
  
  double x,y,z, xRot,yRot,zRot, xOrient,yOrient,zOrient;
    
  //Prepare in file from createTree	
  char rootf[10] = ".root";
  char inf[50] = "../outputFiles/yourFileName_";//***
  strcat(inf,argv[1]);
  strcat(inf,rootf);   
  TFile* phaseFile = new TFile(inf,"update"); //Associate any new objects to the proper working file
 
  TTree* t = (TTree*)phaseFile->Get("OutTree"); //name from create Tree!   
  TTree *tt = new TTree("filteredTree", "Updated collection of tracker variables"); //***SET MANUALLY*** name of the new filtered tree
  
//--------Head phantom hull map with rotation------------
  float headAngle = (float)atoi(argv[1]);
  double radConv = 3.14159265/180.0;
  TH3S* HUMap = new TH3S("HullMap", "Hounsfield Units", 1088,-102,102, 160,-100,100, 1600,-150,150); //Head phantom   
  for( int binzid = 1; binzid<=(HUnits->GetZaxis()->GetNbins()); binzid++){
	  for( int binyid = 1; binyid<=(HUnits->GetYaxis()->GetNbins()); binyid++){
      for( int binxid = 1; binxid<=(HUnits->GetXaxis()->GetNbins()); binxid++ ){ 
	   
	      int HU = (int)HUnits->GetBinContent(binxid,binyid,binzid); 	  
	      x = (HUnits->GetXaxis()->GetBinCenter(binxid))-6.0; //translate midpoint to origin (=simulation origin)
	      y = (HUnits->GetYaxis()->GetBinCenter(binyid))-6.0;
        z = (HUnits->GetZaxis()->GetBinCenter(binzid))+86.25; //Move it down one Zbin! It is off by one pixel in z direction 86.25mm -> 85.0mm
        
        xRot = (x*cos(double(headAngle)*radConv) - (y*sin(double(headAngle)*radConv))); //Rotate around z-axis
        yRot = (x*sin(double(headAngle)*radConv) + (y*cos(double(headAngle)*radConv))); 
	
        xOrient = xRot; //Translate phantom placement to simulation coordinates
        yOrient = yRot;
        zOrient = z*(-1.); //The phantom is upside-down in the HU-file
 
        HUMap->Fill(xOrient,zOrient,yOrient,HU);
      }
    }
  }
  phantomFile->Close(); //This deletes any objects associated with the phantom file from memory  
//--------------------------------------  
  
  t->SetBranchAddress("TPSx0",&Point.TPSx0);
  t->SetBranchAddress("TPSy0",&Point.TPSy0);
  t->SetBranchAddress("TPSz0",&Point.TPSz0);
  
  t->SetBranchAddress("TPSpx0",&Point.TPSpx0);
  t->SetBranchAddress("TPSpy0",&Point.TPSpy0);
  t->SetBranchAddress("TPSpz0",&Point.TPSpz0);

  t->SetBranchAddress("x21",&Point.x21);
  t->SetBranchAddress("y21",&Point.y21);
  t->SetBranchAddress("z21",&Point.z21);
  
  t->SetBranchAddress("p2x",&Point.p2x);
  t->SetBranchAddress("p2y",&Point.p2y);
  t->SetBranchAddress("p2z",&Point.p2z); 
  
  t->SetBranchAddress("WEPLreal",&Point.weplReal);
  t->SetBranchAddress("scatter_out",&Point.scatter_out);
  
  //To be updated after hull
  tt->Branch("x12_h",&Point.TPSx0);
  tt->Branch("y12_h",&Point.TPSy0);
  tt->Branch("z12_h",&Point.TPSz0);
  
  tt->Branch("p1x_h",&Point.TPSpx0);
  tt->Branch("p1y_h",&Point.TPSpy0);
  tt->Branch("p1z_h",&Point.TPSpz0);
  
  tt->Branch("x21_h",&Point.x21);
  tt->Branch("y21_h",&Point.y21);
  tt->Branch("z21_h",&Point.z21);
  
  tt->Branch("p2x_h",&Point.p2x);
  tt->Branch("p2y_h",&Point.p2y);
  tt->Branch("p2z_h",&Point.p2z); 
  
  tt->Branch("weplReal_h",&Point.weplReal);
  tt->Branch("scatter_out_h",&Point.scatter_out);
  
  int NEntries = t->GetEntries(); 
 
  cout << "Begin filter initialization" << endl;   //We filter on information from protons exiting phantom! (exit-enter)

  TH3F* angularDistX = new TH3F("angularDistX", "Angular distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-0.2,0.2);
  TH3F* angularDistZ = new TH3F("angularDistZ", "Angular distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-0.2,0.2);
  TH3F* WEPLDist = new TH3F("WEPLDist", "WEPL distribution", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 125,-50,200); //only wepl values up to 200 mm is expected

  float distAngleX, distAngleZ, WEPLvalue;
  float meanAngleX, sigAngleX, meanAngleZ, sigAngleZ, meanWEPL, sigWEPL;
  float s_pos;
  
  cout<<NEntries<<endl;
  for(int i=0;i<NEntries;i++){//Loop over all protons to fill distributions
    t->GetEntry(i);
    if(i%100000 == 0) cout<<i<<endl;
    
    Point.x21 = Point.x21 - (Point.p2x*d_phantom);
    Point.z21 = Point.z21 - (Point.p2z*d_phantom);
    Point.y21 = Point.y21 - d_phantom;
   
    s_pos = 0.3;//cm beam sigma at surface of phantom (known from simulations)
    
		float step = 0.1; // mm Step in that you project the vectors onto the hull. The finer the step, the better the hull algorithm, but the more time consuming
		TVector3 X2_hull = TVector3(Point.x21, Point.z21, Point.y21); //End-point (for that particular proton)
		TVector3 dirX2 = TVector3(Point.p2x, Point.p2z, Point.p2y); //final direction
		int counter1 = 0;  //Counter is a sefety measure to prevent infinite loops 
		int HUbin = HUMap->FindBin(X2_hull.x(), X2_hull.y(), X2_hull.z()); // assuming both HU and measured positions are in the same coordinates
		int HU = HUMap->GetBinContent(HUbin);
		if(HU == 0){HU=(-998);}
		while(HU==(-998) && counter1 < 2000){ // project as long as the HU value is air. Counter is a sefety measure to prevent infinite loops 
	      X2_hull = X2_hull - step*dirX2; // Project in straight line along initial direction NOTE THE MINUS!
	      HUbin = HUMap->FindBin(X2_hull.x(), X2_hull.y(), X2_hull.z());
	      HU = HUMap->GetBinContent(HUbin);
	      if(HU == 0 || HU==2310){HU=(-998);}
	      counter1++;
		}
		 
	  if(X2_hull.z()>(-50.0)){ //hit the hull from the back
      Point.x21 = X2_hull.x();
      Point.z21 = X2_hull.y();
      Point.y21 = X2_hull.z();
    } 
	
	  struct inPos posIn = ComputeOptIn(&Point, s_pos); //optimized starting point at hull depth if found, otherwise uses the projected tracker position
    TVector3 X0_hull = TVector3(posIn.X_opt, posIn.Z_opt, Point.TPSy0);
    TVector3 dirX0 = TVector3(posIn.theta_X_opt, posIn.theta_Z_opt, Point.TPSpy0);
	  int counter2 = 0; //resetting counter
	  HUbin = HUMap->FindBin(X0_hull.x(), X0_hull.y(), X0_hull.z()); 
	  HU = HUMap->GetBinContent(HUbin);
	  if(HU == 0)HU=(-998);
		while(HU==(-998) && counter2 < 2000){ // project as long as the HU value in the protons way is smaller than lung material. Counter is a sefety measure to prevent infinite loops 
	      X0_hull = X0_hull + step*dirX0; // Project in straight line along initial direction
	      HUbin = HUMap->FindBin(X0_hull.x(), X0_hull.y(), X0_hull.z());
		    HU = HUMap->GetBinContent(HUbin);
	      if(HU == 0 || HU==2310){HU=(-998);}
	      counter2++;
		}
	
 	  if(X0_hull.z()<50.0 && X2_hull.z()>(-50.0)){//Both sides have hit the pahantom	
		
	    Point.TPSx0 = X0_hull.x(); 
      Point.TPSz0 = X0_hull.y(); 
      Point.TPSy0 = X0_hull.z(); 

	    Point.TPSpx0 = posIn.theta_X_opt; 
      Point.TPSpz0 = posIn.theta_Z_opt; 
                
      Point.x21 = X2_hull.x();
      Point.z21 = X2_hull.y();
      Point.y21 = X2_hull.z();
       
      tt->Fill();
       
      distAngleX = Point.p2x - Point.TPSpx0;
      angularDistX->Fill(Point.x21, Point.z21, distAngleX);    
      distAngleZ = Point.p2z - Point.TPSpz0;
      angularDistZ->Fill(Point.x21, Point.z21, distAngleZ);
    }
    
    else if(X0_hull.z()>50.0 && X2_hull.z()<(-50.0)){//Completely missed the phantom from both sides
		
	    Point.TPSx0 = posIn.X_opt; 
      Point.TPSz0 = posIn.Z_opt;
       
      Point.TPSpx0 = posIn.theta_X_opt; 
      Point.TPSpz0 = posIn.theta_Z_opt;
       
      tt->Fill();
       
      distAngleX = Point.p2x - Point.TPSpx0;
      angularDistX->Fill(Point.x21, Point.z21, distAngleX);    
      distAngleZ = Point.p2z - Point.TPSpz0;
      angularDistZ->Fill(Point.x21, Point.z21, distAngleZ);
	  }
  }
  
  //Filtering on angle all done and we have the best possible positions 
  //Now we have a new file where positions are on the hull and we have prepared for the filtering
  tt->Write("", TObject::kOverwrite);
  
  tt->SetBranchAddress("x12_h",&Point.x12_h);
  tt->SetBranchAddress("y12_h",&Point.y12_h);
  tt->SetBranchAddress("z12_h",&Point.z12_h);
  
  tt->SetBranchAddress("p1x_h",&Point.p1x_h);
  tt->SetBranchAddress("p1y_h",&Point.p1y_h);
  tt->SetBranchAddress("p1z_h",&Point.p1z_h);
  
  tt->SetBranchAddress("x21_h",&Point.x21_h);
  tt->SetBranchAddress("y21_h",&Point.y21_h);
  tt->SetBranchAddress("z21_h",&Point.z21_h);
  
  tt->SetBranchAddress("p2x_h",&Point.p2x_h);
  tt->SetBranchAddress("p2y_h",&Point.p2y_h);
  tt->SetBranchAddress("p2z_h",&Point.p2z_h); 
  
  tt->SetBranchAddress("weplReal_h",&Point.weplReal_h);
  tt->SetBranchAddress("scatter_out_h",&Point.scatter_out_h);
    
  NEntries = tt->GetEntries(); //update entries
  
  cout<<NEntries<<" new tree"<<endl;
 
  TH2F* angleMeanX = new TH2F("angleMeanX", "Angular distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* angleMeanZ = new TH2F("angleMeanZ", "Angular distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* WEPLMean = new TH2F("WEPLMean", "WEPL distribution", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
	
  TH2F* angleSigX = new TH2F("angleSigX", "Angular sigma distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* angleSigZ = new TH2F("angleSigZ", "Angular sigma distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* WEPLSig = new TH2F("WEPLSig", "WEPL sigma distribution", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
	
  for( int binxid = 1; binxid<=(WEPLDist->GetXaxis()->GetNbins()); binxid++ ){  if(binxid%100 == 0) cout<<binxid<<endl; //For visualizing progress
    for( int binzid = 1; binzid<=(WEPLDist->GetYaxis()->GetNbins()); binzid++){
		
      float x = WEPLDist->GetXaxis()->GetBinCenter(binxid); //The "local" bin's center is requested so that the pixel from the 3D histogram can be translated to 2D
	    float z = WEPLDist->GetYaxis()->GetBinCenter(binzid);
		  
	    //Begin filling up the 2D histograms where each bin is a double image pixel and filled/weighted by the mean or mean-error
		
	    meanAngleX = angularDistX->ProjectionZ("pz1",binxid,binxid,binzid,binzid)->GetMean();
	    angleMeanX->Fill(x,z,meanAngleX);
		
	    sigAngleX = angularDistX->ProjectionZ("pz2",binxid,binxid,binzid,binzid)->GetStdDev();
	    angleSigX->Fill(x,z,sigAngleX);
		
	    meanAngleZ = angularDistZ->ProjectionZ("pz3",binxid,binxid,binzid,binzid)->GetMean();
	    angleMeanZ->Fill(x,z,meanAngleZ);
		
	    sigAngleZ = angularDistZ->ProjectionZ("pz4",binxid,binxid,binzid,binzid)->GetStdDev();
	    angleSigZ->Fill(x,z,sigAngleZ);		
	  }
  }
 
  angleMeanX->Write("MeanAx",TObject::kOverwrite);
  angleSigX->Write("SigAx",TObject::kOverwrite);
  angleMeanZ->Write("MeanAz",TObject::kOverwrite);
  angleSigZ->Write("SigAz",TObject::kOverwrite);

  //Filter the wepl!
  TH1D* WEPLdistf;
  for(int i=0;i<NEntries;i++){//Loop over all protons in the new file/tree to fill distributions
    tt->GetEntry(i); 
            
    int binGlobal = angleMeanX->FindBin(Point.x21_h, Point.z21_h); 
    float angleX = Point.p2x_h - Point.p1x_h;
    float angleZ = Point.p2z_h - Point.p1z_h;
 
    if( angleX >= (angleMeanX->GetBinContent(binGlobal) - 2.5*angleSigX->GetBinContent(binGlobal)) && angleX <= (angleMeanX->GetBinContent(binGlobal) + 2.5*angleSigX->GetBinContent(binGlobal)) && 
        angleZ >= (angleMeanZ->GetBinContent(binGlobal) - 2.5*angleSigZ->GetBinContent(binGlobal)) && angleZ <= (angleMeanZ->GetBinContent(binGlobal) + 2.5*angleSigZ->GetBinContent(binGlobal))
      ) 
      {
       WEPLDist->Fill(Point.x21_h, Point.z21_h, Point.weplReal_h); 
       if(i%100000 == 0) cout<<i<<endl;
    }
  }
  
  for( int binxid = 1; binxid<=(WEPLDist->GetXaxis()->GetNbins()); binxid++ ){  if(binxid%100 == 0) cout<<binxid<<endl; //For visualizing progress
	  for( int binzid = 1; binzid<=(WEPLDist->GetYaxis()->GetNbins()); binzid++){
		
		  double x = WEPLDist->GetXaxis()->GetBinCenter(binxid); //The "local" bin's center is requested so that the pixel from the 3D histogram can be translated to 2D
		  double z = WEPLDist->GetYaxis()->GetBinCenter(binzid);
		  
		  //Finding the maximum WEPL peak inside a bin/pixel	  
		  TH1D* WEPLdistPix = WEPLDist->ProjectionZ("pz15",binxid,binxid,binzid,binzid);
		  int binmax = WEPLdistPix->GetMaximumBin(); 
		  int bin = binmax;
	    float mode = WEPLdistPix->GetXaxis()->GetBinCenter(binmax);
	    float maxValue = WEPLdistPix->GetBinContent(binmax); 
	    float binValue = WEPLdistPix->GetBinContent(bin);
	    int nbin = 0; //prevent infinite loop
	      
	    while(binValue>=(maxValue/2.0) && nbin<11){ 
			  bin=bin+1;
		    binValue = WEPLdistPix->GetBinContent(bin);
		    nbin++;
		  }
		  int topBin = bin; 
	    float modeTop = WEPLdistPix->GetXaxis()->GetBinCenter(topBin); 
	    nbin=0;
	    bin=binmax;
	    binValue = WEPLdistPix->GetBinContent(bin);
	    while(binValue>=(maxValue/2.0) && nbin<11){
			  bin=bin-1;
		    binValue = WEPLdistPix->GetBinContent(bin);
		    nbin++;
		  }
		  int botBin = bin; 
	    float modeBot = WEPLdistPix->GetXaxis()->GetBinCenter(botBin);
	      
	    float fwhm = modeTop - modeBot;
	      
	    WEPLdistf = new TH1D("weplDistf", "distribution of wepl", 50, (mode-(fwhm)), (mode+(fwhm)));
	    for( int binpixid = 1; binpixid<=(WEPLdistPix->GetXaxis()->GetNbins()); binpixid++ ){ 
			  float WEPLvalue = WEPLdistPix->GetBinContent(binpixid);
			  float xpos = WEPLdistPix->GetXaxis()->GetBinCenter(binpixid);
			  WEPLdistf->Fill(xpos,WEPLvalue);
		  }
	      
	    float mean = WEPLdistf->GetMean();
		  WEPLMean->Fill(x,z,mean);
		  float sig = WEPLdistf->GetStdDev();
	   	WEPLSig->Fill(x,z,sig);
	   
	   	delete WEPLdistf;
    }
  }
   
  WEPLMean->Write("MeanWEPL",TObject::kOverwrite);
  WEPLSig->Write("SigWEPL",TObject::kOverwrite); 
  
  phaseFile->Close();
  
  cout << "Filter initialization all done" << endl;
  return 0; 
}

////////////////////////////////////////////
// Compute optimized entrance point
////////////////////////////////////////////
struct inPos ComputeOptIn(Proton *Point, float s_pos){
	
	//---------------------------------------------------------------
    //transform all mm into cm 
    //---------------------------------------------------------------
    Point->TPSx0 /=10.; Point->TPSy0 /=10.; Point->TPSz0 /=10.;
    Point->x21	/=10.; Point->y21 /=10.; Point->z21 /=10.; 
    //------------------------------------------------------------  
	
  TVector3 m0(Point->TPSx0, Point->TPSy0, Point->TPSz0);
  TVector3 p0(Point->TPSpx0, Point->TPSpy0, Point->TPSpz0);
  
  TVector3 m1(Point->x21, Point->y21, Point->z21);
  TVector3 p1(Point->p2x, Point->p2y, Point->p2z);
   
        //Initialization of the MLP
    float X_mlp,Z_mlp, theta_X_mlp,theta_Z_mlp, X_mlp_sigma,Z_mlp_sigma;
    
    //Iterator for the matrix operations (needed later)
    int a;
    
    //Matrices    
    double T_beam[4];
      T_beam[0] = 1;
      T_beam[1] = 0;
      T_beam[2] = 1/d_source;
      T_beam[3] = 1;  
    double T_beam_transpose[4];
      T_beam_transpose[0] = 1;
      T_beam_transpose[1] = 1/d_source;
      T_beam_transpose[2] = 0;
      T_beam_transpose[3] = 1;
    
    double T_out[4];
      T_out[0] = 1; //eq. 26 in Krah
      T_out[1] = 0; 
      T_out[2] = -1/(d_T);
      T_out[3] = 1/(d_T);
    double T_out_transpose[4];
      T_out_transpose[0] = T_out[0];
      T_out_transpose[1] = T_out[2];
      T_out_transpose[2] = T_out[1];
      T_out_transpose[3] = T_out[3];
      
    double beam_uncert[4]; //sigma spot-size
      beam_uncert[0] = s_pos*s_pos;
      beam_uncert[1] = 0;
      beam_uncert[2] = 0;
      beam_uncert[3] = 0;
    
    double beam_1[4]; //T_beam*beam_uncert1
      beam_1[0] = ((T_beam[0] * beam_uncert[0]) + (T_beam[1]*beam_uncert[2]));
      beam_1[1] = ((T_beam[0] * beam_uncert[1]) + (T_beam[1]*beam_uncert[3]));
      beam_1[2] = ((T_beam[2] * beam_uncert[0]) + (T_beam[3]*beam_uncert[2]));
      beam_1[3] = ((T_beam[2] * beam_uncert[1]) + (T_beam[3]*beam_uncert[3]));  

    double Sigma_in[4]={0}; //Single-sided: eq.29 in Krah, now eq.28
      Sigma_in[0] = ( (beam_1[0] * T_beam_transpose[0]) + (beam_1[1] * T_beam_transpose[2]) );
      Sigma_in[1] = ( (beam_1[0] * T_beam_transpose[1]) + (beam_1[1] * T_beam_transpose[3]) );
      Sigma_in[2] = ( (beam_1[2] * T_beam_transpose[0]) + (beam_1[3] * T_beam_transpose[2]) );
      Sigma_in[3] = ( ((beam_1[2] * T_beam_transpose[1]) + (beam_1[3] * T_beam_transpose[3]) )) + pow(s_angle, 2); 
    
    double Sigma_out[4];
      Sigma_out[0] = pow(s_pos_out, 2) * ( (T_out[0] * T_out_transpose[0]) + (T_out[1]*T_out_transpose[2]) );
      Sigma_out[1] = pow(s_pos_out, 2) * ( (T_out[0] * T_out_transpose[1]) + (T_out[1]*T_out_transpose[3]) );
      Sigma_out[2] = pow(s_pos_out, 2) * ( (T_out[2] * T_out_transpose[0]) + (T_out[3]*T_out_transpose[2]) );
      Sigma_out[3] = (pow(s_pos_out, 2) * ( (T_out[2] * T_out_transpose[1]) + (T_out[3]*T_out_transpose[3]) )) + pow(Point->scatter_out, 2);  
    
    double S_in[4]; // eq. 14 in Krah
      S_in[0] = 1;
      S_in[1] = d_entry;
      S_in[2] = 0;
      S_in[3] = 1;
    double S_in_transpose[4];
      S_in_transpose[0] = S_in[0];
      S_in_transpose[1] = S_in[2];
      S_in_transpose[2] = S_in[1];
      S_in_transpose[3] = S_in[3];
      
    double S_out_inverse[4]; //inverse of a two by two matrix ((a b),(c d)): 1/determinant ((d,-b),(-c,a))
      S_out_inverse[0] = 1;
      S_out_inverse[1] = -d_exit;
      S_out_inverse[2] = 0;
      S_out_inverse[3] = 1;
    double S_out_inverse_transpose[4];
      S_out_inverse_transpose[0] = 1;
      S_out_inverse_transpose[1] = 0;
      S_out_inverse_transpose[2] = -d_exit;
      S_out_inverse_transpose[3] = 1;
        
      //Can calculate parts of the C1 and C2 terms for later use  
    double SS_in[4]; //S_in*Sigma_in
      SS_in[0] = (S_in[0] * Sigma_in[0]) + (S_in[1] * Sigma_in[2]);
      SS_in[1] = (S_in[0] * Sigma_in[1]) + (S_in[1] * Sigma_in[3]);
      SS_in[2] = (S_in[2] * Sigma_in[0]) + (S_in[3] * Sigma_in[2]);
      SS_in[3] = (S_in[2] * Sigma_in[1]) + (S_in[3] * Sigma_in[3]);
    double SSS_in[4]; //S_in*Sigma_in*S_in_transpose (to be multiplied with R_0 and R_0_transpose later)
      SSS_in[0] = (SS_in[0] * S_in_transpose[0]) + (SS_in[1] * S_in_transpose[2]);
      SSS_in[1] = (SS_in[0] * S_in_transpose[1]) + (SS_in[1] * S_in_transpose[3]);
      SSS_in[2] = (SS_in[2] * S_in_transpose[0]) + (SS_in[3] * S_in_transpose[2]);
      SSS_in[3] = (SS_in[2] * S_in_transpose[1]) + (SS_in[3] * S_in_transpose[3]);

    double SS_out[4]; //S_out_inverse*Sigma_out
      SS_out[0] = (S_out_inverse[0] * Sigma_out[0]) + (S_out_inverse[1] * Sigma_out[2]);
      SS_out[1] = (S_out_inverse[0] * Sigma_out[1]) + (S_out_inverse[1] * Sigma_out[3]);
      SS_out[2] = (S_out_inverse[2] * Sigma_out[0]) + (S_out_inverse[3] * Sigma_out[2]);
      SS_out[3] = (S_out_inverse[2] * Sigma_out[1]) + (S_out_inverse[3] * Sigma_out[3]);
    double SSS_out[4]; //S_out_inverse*Sigma_out*S_out_inverse_transpose (to be multiplied with R_1_inverse and R_1_inverse_transpose later)
      SSS_out[0] = (SS_out[0] * S_out_inverse_transpose[0]) + (SS_out[1] * S_out_inverse_transpose[2]);
      SSS_out[1] = (SS_out[0] * S_out_inverse_transpose[1]) + (SS_out[1] * S_out_inverse_transpose[3]);
      SSS_out[2] = (SS_out[2] * S_out_inverse_transpose[0]) + (SS_out[3] * S_out_inverse_transpose[2]);
      SSS_out[3] = (SS_out[2] * S_out_inverse_transpose[1]) + (SS_out[3] * S_out_inverse_transpose[3]);    
    
    double R_0[4]={0};
    double R_0_transpose[4]={0};
    double RSSS_0[4]={0};
    double RS_0[4]={0};
    double R_1_inverse[4]={0};
    double R_1_inverse_transpose[4]={0};
    double RSSS_1[4]={0};
    double RS_1[4]={0};
    double RS_2[4]={0};
    
    double Sigma_1[4]={0};
    double Sigma_2[4]={0};
    
    double y_0[2]={0};
    double y_2[2]={0};
    
    double C1_1[4]={0};
    double C1[4]={0};
    
    double C2_1[4]={0};
    double C2_2[4]={0};
    double C2[4]={0};
    
    double C12[4]={0};
    double C12_inverse[4]={0};
    
    double first_first[4]={0};
    double second_first[4]={0};
    double first_second[2]={0};
    double second_second[2]={0};
    double first[2] = {0};
    double second[2]={0};
    
    //Parameter initialization
    double sy1, sy2, st1, st2, sty1, sty2;
    double determinant_1, determinant_2, determinant_C12;
    
    //Initialize the MLP iterators
    double step_length = (m1.y()-m0.y())/Nsteps;
    double posy=m0.y()+step_length;
                    
        //Transvection matrices, eq.8 in Krah  
        R_0[0] = 1;
        R_0[1] = posy-m0.y();
        R_0[2] = 0;
        R_0[3] = 1;
        R_0_transpose[0] = 1;
        R_0_transpose[1] = 0;
        R_0_transpose[2] = posy-m0.y();
        R_0_transpose[3] = 1;

        R_1_inverse[0] = 1;
        R_1_inverse[1] = -(m1.y()-posy);
        R_1_inverse[2] = 0;
        R_1_inverse[3] = 1;
        R_1_inverse_transpose[0] = 1; 
        R_1_inverse_transpose[1] = 0;
        R_1_inverse_transpose[2] = -(m1.y()-posy);
        R_1_inverse_transpose[3] = 1;
        
        //need these for later:R_0*S_in, R_0*SSS_in, and R_1_inverse*SSS_out
        RS_0[0]= (R_0[0]*S_in[0]) + (R_0[1]*S_in[2]);
        RS_0[1]= (R_0[0]*S_in[1]) + (R_0[1]*S_in[3]);
        RS_0[2]= (R_0[2]*S_in[0]) + (R_0[3]*S_in[2]);
        RS_0[3]= (R_0[2]*S_in[1]) + (R_0[3]*S_in[3]);
        
        RSSS_0[0]= (R_0[0]*SSS_in[0]) + (R_0[1]*SSS_in[2]);
        RSSS_0[1]= (R_0[0]*SSS_in[1]) + (R_0[1]*SSS_in[3]);
        RSSS_0[2]= (R_0[2]*SSS_in[0]) + (R_0[3]*SSS_in[2]);
        RSSS_0[3]= (R_0[2]*SSS_in[1]) + (R_0[3]*SSS_in[3]);
        
        RS_1[0]= (R_1_inverse[0]*S_out_inverse[0]) + (R_1_inverse[1]*S_out_inverse[2]);
        RS_1[1]= (R_1_inverse[0]*S_out_inverse[1]) + (R_1_inverse[1]*S_out_inverse[3]);
        RS_1[2]= (R_1_inverse[2]*S_out_inverse[0]) + (R_1_inverse[3]*S_out_inverse[2]);
        RS_1[3]= (R_1_inverse[2]*S_out_inverse[1]) + (R_1_inverse[3]*S_out_inverse[3]);
        
        RSSS_1[0]= (R_1_inverse[0]*SSS_out[0]) + (R_1_inverse[1]*SSS_out[2]);
        RSSS_1[1]= (R_1_inverse[0]*SSS_out[1]) + (R_1_inverse[1]*SSS_out[3]);
        RSSS_1[2]= (R_1_inverse[2]*SSS_out[0]) + (R_1_inverse[3]*SSS_out[2]);
        RSSS_1[3]= (R_1_inverse[2]*SSS_out[1]) + (R_1_inverse[3]*SSS_out[3]);
                            
        // First do everything not depending on the direction of interest (either X_mlp or Z_mlp)
        //scattering sigma matrices
        sy1 = Sigmay1(posy-m0.y());
        st1 = Sigmat1(posy-m0.y());
        sty1 = Sigmaty1(posy-m0.y());
            
        sy2 = Sigmay2(m1.y()-m0.y(),posy-m0.y());
        sty2 = Sigmaty2(m1.y()-m0.y(),posy-m0.y());
        st2 = Sigmat2(m1.y()-m0.y(),posy-m0.y());

        // Scattering sigma matrices
        Sigma_1[0] = sy1;
        Sigma_1[1] = sty1;
        Sigma_1[2] = sty1;
        Sigma_1[3] = st1;
        
        Sigma_2[0] = sy2;
        Sigma_2[1] = sty2;
        Sigma_2[2] = sty2;
        Sigma_2[3] = st2;
        
        // Calculate the pre factors C1 and C2 as in Krah et al. (2018): C2*(C1+C2)^-1 R0*S0*Y0 + C1*(C1+C2)^1 R1^-1*S1^-1*Y2
        // First start with the C1 = ((R_0*S_in*Sigma_in*S_in_transpose)*(R_0_transpose))+Sigma_1
        C1_1[0]=(RSSS_0[0]*R_0_transpose[0])+(RSSS_0[1]*R_0_transpose[2]);
        C1_1[1]=(RSSS_0[0]*R_0_transpose[1])+(RSSS_0[1]*R_0_transpose[3]);
        C1_1[2]=(RSSS_0[2]*R_0_transpose[0])+(RSSS_0[3]*R_0_transpose[2]);
        C1_1[3]=(RSSS_0[2]*R_0_transpose[1])+(RSSS_0[3]*R_0_transpose[3]);

        for (a=0;a<4;a++){
            C1[a] = C1_1[a]+Sigma_1[a];
        }
        
        //Now calculate C2 = (R_1_inverse*S_out_inverse*Sigma_out*S_out_inverse_transpose*R_1_inverse_transpose) + (R_1_inverse*Sigma_2*R_1_inverse_transpose)        
        C2_1[0] = (RSSS_1[0] * R_1_inverse_transpose[0]) + (RSSS_1[1] * R_1_inverse_transpose[2]);
        C2_1[1] = (RSSS_1[0] * R_1_inverse_transpose[1]) + (RSSS_1[1] * R_1_inverse_transpose[3]);
        C2_1[2] = (RSSS_1[2] * R_1_inverse_transpose[0]) + (RSSS_1[3] * R_1_inverse_transpose[2]);
        C2_1[3] = (RSSS_1[2] * R_1_inverse_transpose[1]) + (RSSS_1[3] * R_1_inverse_transpose[3]);

        RS_2[0] = (R_1_inverse[0] * Sigma_2[0]) + (R_1_inverse[1] * Sigma_2[2]);
        RS_2[1] = (R_1_inverse[0] * Sigma_2[1]) + (R_1_inverse[1] * Sigma_2[3]);
        RS_2[2] = (R_1_inverse[2] * Sigma_2[0]) + (R_1_inverse[3] * Sigma_2[2]);
        RS_2[3] = (R_1_inverse[2] * Sigma_2[1]) + (R_1_inverse[3] * Sigma_2[3]);
        
        C2_2[0] = (RS_2[0] * R_1_inverse_transpose[0]) + (RS_2[1] * R_1_inverse_transpose[2]);
        C2_2[1] = (RS_2[0] * R_1_inverse_transpose[1]) + (RS_2[1] * R_1_inverse_transpose[3]);
        C2_2[2] = (RS_2[2] * R_1_inverse_transpose[0]) + (RS_2[3] * R_1_inverse_transpose[2]);
        C2_2[3] = (RS_2[2] * R_1_inverse_transpose[1]) + (RS_2[3] * R_1_inverse_transpose[3]);
        
        C2[0] = C2_1[0]+C2_2[0];
        C2[1] = C2_1[1]+C2_2[1];
        C2[2] = C2_1[2]+C2_2[2];
        C2[3] = C2_1[3]+C2_2[3];
        
        //Add the second factor to the first to yield C1+C2
        for(a=0;a<4;a++){
            C12[a]=C1[a]+C2[a];
        }

        //invert so to get the prefactor (C1+C2)^-1
        determinant_C12=(C12[0]*C12[3])-(C12[1]*C12[2]);
        C12_inverse[0]=C12[3]/determinant_C12;
        C12_inverse[1]=-C12[1]/determinant_C12;
        C12_inverse[2]=-C12[2]/determinant_C12;
        C12_inverse[3]=C12[0]/determinant_C12;

        //Multiply C2 to yield the first prefactor C2*(C1+C2)^-1
        first_first[0]=(C2[0]*C12_inverse[0])+(C2[1]*C12_inverse[2]);
        first_first[1]=(C2[0]*C12_inverse[1])+(C2[1]*C12_inverse[3]);
        first_first[2]=(C2[2]*C12_inverse[0])+(C2[3]*C12_inverse[2]);
        first_first[3]=(C2[2]*C12_inverse[1])+(C2[3]*C12_inverse[3]);

        //Same with C1 to yield the second prefactor C1*(C1+C2)^-1
        second_first[0]=(C1[0]*C12_inverse[0])+(C1[1]*C12_inverse[2]);
        second_first[1]=(C1[0]*C12_inverse[1])+(C1[1]*C12_inverse[3]);
        second_first[2]=(C1[2]*C12_inverse[0])+(C1[3]*C12_inverse[2]);
        second_first[3]=(C1[2]*C12_inverse[1])+(C1[3]*C12_inverse[3]);

        //Now the second part of each factor (start with X drection)
        y_0[0] = m0.x();
        y_0[1] = tan(p0.x());
        
        y_2[0] = m1.x();
        y_2[1] = tan(p1.x());

        // Start with R_0*S_in*Y0
        first_second[0] = (RS_0[0]*y_0[0])+(RS_0[1]*y_0[1]);
        first_second[1] = (RS_0[2]*y_0[0])+(RS_0[3]*y_0[1]);

        // Now R1_inverse*S_out_inverse*Y2
        second_second[0]=(RS_1[0]*y_2[0])+(RS_1[1]*y_2[1]);
        second_second[1]=(RS_1[2]*y_2[0])+(RS_1[3]*y_2[1]);
        
        // Put Everything together: (C2*(C1+C2)^-1)*(R_0*S_in*Y0)
        first[0]=(first_first[0]*first_second[0])+(first_first[1]*first_second[1]);
        first[1]=(first_first[2]*first_second[0])+(first_first[3]*first_second[1]);
        //+(C1*(C1+C2)^-1)*(R_1_inverse*S_out_inverse*Y2)
        second[0]=(second_first[0]*second_second[0])+(second_first[1]*second_second[1]);
        second[1]=(second_first[2]*second_second[0])+(second_first[3]*second_second[1]);
        
      X_mlp=(first[0]+second[0]); //=C2*(C1+C2)^-1 R0*S0*Y0 + C1*(C1+C2)^1 R1^-1*S1^-1*Y2
      theta_X_mlp=(first[1]+second[1]);
        
        // Now do the other direction
        y_0[0] = m0.z();
        y_0[1] = tan(p0.z());
        
        y_2[0] = m1.z();
        y_2[1] = tan(p1.z());

        // Again with respect to the other direction
        first_second[0] = (RS_0[0]*y_0[0])+(RS_0[1]*y_0[1]);
        first_second[1] = (RS_0[2]*y_0[0])+(RS_0[3]*y_0[1]);

        second_second[0]=(RS_1[0]*y_2[0])+(RS_1[1]*y_2[1]);
        second_second[1]=(RS_1[2]*y_2[0])+(RS_1[3]*y_2[1]);
        
        // Put Everything together again
        first[0]=(first_first[0]*first_second[0])+(first_first[1]*first_second[1]);
        first[1]=(first_first[2]*first_second[0])+(first_first[3]*first_second[1]);

        second[0]=(second_first[0]*second_second[0])+(second_first[1]*second_second[1]);
        second[1]=(second_first[2]*second_second[0])+(second_first[3]*second_second[1]);
        
      Z_mlp=(first[0]+second[0]);
      theta_Z_mlp=(first[1]+second[1]);
      
   	//---------------------------------------------------------------
    //transform all cm back into mm again 
    //---------------------------------------------------------------
    Point->TPSx0 *=10.; Point->TPSy0 *=10.; Point->TPSz0 *=10.;
    Point->x21	*=10.; Point->y21 *=10.; Point->z21 *=10.; 
    //------------------------------------------------------------   

   struct inPos posIn = {X_mlp*10.0, Z_mlp*10.0, theta_X_mlp, theta_Z_mlp};
   return posIn;
}	

double Sigmat1(double position)
{
	double p = position;
	double sigt1 = (azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigt1/X_0);
}

double Sigmay1(double position)
{
	double p = position;
	double sigy1 = (azero*p*p*p/3)+(aone*p*p*p*p/12)+(atwo*p*p*p*p*p/30)+(athree*p*p*p*p*p*p/60)+(afour*p*p*p*p*p*p*p/105)+(afive*p*p*p*p*p*p*p*p/168);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigy1/X_0);
}

double Sigmaty1(double position)
{
	double p = position;
	double sigty1 = (azero*p*p/2)+(aone*p*p*p/6)+(atwo*p*p*p*p/12)+(athree*p*p*p*p*p/20)+(afour*p*p*p*p*p*p/30)+(afive*p*p*p*p*p*p*p/42);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigty1/X_0);
}

double Sigmat2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigt2 = ((azero*s)+(aone*s*s/2)+(atwo*s*s*s/3)+(athree*s*s*s*s/4)+(afour*s*s*s*s*s/5)+(afive*s*s*s*s*s*s/6))-((azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigt2/X_0);
}

double Sigmay2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigy2 = (azero*s*s*s/3)+(aone*s*s*s*s/12)+(atwo*s*s*s*s*s/30)+(athree*s*s*s*s*s*s/60)+(afour*s*s*s*s*s*s*s/105)+(afive*s*s*s*s*s*s*s*s/168)-((azero*s*s*p)+(((aone*s*s/2)-(azero*s))*p*p)+(((atwo*s*s/3)-(2*aone*s/3)+(azero/3))*p*p*p)+(((athree*s*s/4)-(atwo*s/2)+(aone/4))*p*p*p*p)+(((afour*s*s/5)-(2*athree*s/5)+(atwo/5))*p*p*p*p*p)+(((afive*s*s/6)-(afour*s/3)+(athree/6))*p*p*p*p*p*p)+(((afour/7)-(2*afive*s/7))*p*p*p*p*p*p*p)+(afive*p*p*p*p*p*p*p*p/8));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigy2/X_0);
}

double Sigmaty2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigty2 = ((azero*s*s/2)+(aone*s*s*s/6)+(atwo*s*s*s*s/12)+(athree*s*s*s*s*s/20)+(afour*s*s*s*s*s*s/30)+(afive*s*s*s*s*s*s*s/42))-((azero*s*p)+(((aone*s)-azero)*p*p/2)+(((atwo*s)-aone)*p*p*p/3)+(((athree*s)-atwo)*p*p*p*p/4)+(((afour*s)-athree)*p*p*p*p*p/5)+(((afive*s)-afour)*p*p*p*p*p*p/6)-(afive*p*p*p*p*p*p*p/7));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigty2/X_0);
}

