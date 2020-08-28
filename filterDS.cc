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

struct Proton{
 
  Float_t x12,y12,z12, x12_h,y12_h,z12_h;
  Float_t p1x,p1y,p1z, p1x_h,p1y_h,p1z_h; 
   
  Float_t x21,y21,z21, x21_h,y21_h,z21_h; 
  Float_t p2x,p2y,p2z, p2x_h,p2y_h,p2z_h;
  
  Float_t scatter_in, scatter_in_h;
  Float_t scatter_out, scatter_out_h;
  
  Float_t weplReal, weplReal_h;
  
};

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
  
  t->SetBranchAddress("x12",&Point.x12);
  t->SetBranchAddress("y12",&Point.y12);
  t->SetBranchAddress("z12",&Point.z12);
  
  t->SetBranchAddress("p1x",&Point.p1x);
  t->SetBranchAddress("p1y",&Point.p1y);
  t->SetBranchAddress("p1z",&Point.p1z);

  t->SetBranchAddress("x21",&Point.x21);
  t->SetBranchAddress("y21",&Point.y21);
  t->SetBranchAddress("z21",&Point.z21);
  
  t->SetBranchAddress("p2x",&Point.p2x);
  t->SetBranchAddress("p2y",&Point.p2y);
  t->SetBranchAddress("p2z",&Point.p2z); 
  
  t->SetBranchAddress("WEPLreal",&Point.weplReal);
  t->SetBranchAddress("scatter_in",&Point.scatter_in);
  t->SetBranchAddress("scatter_out",&Point.scatter_out);
  
  //Update after hull
  tt->Branch("x12_h",&Point.x12);
  tt->Branch("y12_h",&Point.y12);
  tt->Branch("z12_h",&Point.z12);
  
  tt->Branch("p1x_h",&Point.p1x);
  tt->Branch("p1y_h",&Point.p1y);
  tt->Branch("p1z_h",&Point.p1z);
  
  tt->Branch("x21_h",&Point.x21);
  tt->Branch("y21_h",&Point.y21);
  tt->Branch("z21_h",&Point.z21);
  
  tt->Branch("p2x_h",&Point.p2x);
  tt->Branch("p2y_h",&Point.p2y);
  tt->Branch("p2z_h",&Point.p2z); 
  
  tt->Branch("weplReal_h",&Point.weplReal);
  tt->Branch("scatter_in_h",&Point.scatter_in);
  tt->Branch("scatter_out_h",&Point.scatter_out);
  
  int NEntries = t->GetEntries(); 
 
  cout << "Begin filter initialization" << endl;   //We filter on information from protons exiting phantom! (exit-enter) Two pixels are put together for increased accuracy

  TH3F* angularDistX = new TH3F("angularDistX", "Angular distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-0.2,0.2);
  TH3F* angularDistZ = new TH3F("angularDistZ", "Angular distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-0.2,0.2);
  //TH3F* lateralDistX = new TH3F("lateralDistX", "Lateral distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-10.0,10.0);
  //TH3F* lateralDistZ = new TH3F("lateralDistZ", "Lateral distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 200,-10.0,10.0);
  TH3F* WEPLDist = new TH3F("WEPLDist", "WEPL distribution", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax, 125,-50,200);

  float distAngleX, distAngleZ, distDeltaX, distDeltaZ, WEPLvalue;
  float meanAngleX, sigAngleX, meanAngleZ, sigAngleZ, meanLatX, sigLatX, meanLatZ, sigLatZ, meanWEPL, sigWEPL;
  float s_pos;
  
  cout<<NEntries<<endl;
  for(int i=0;i<NEntries;i++){//Loop over all protons to fill distributions
    t->GetEntry(i);
    if(i%100000 == 0) cout<<i<<endl;
   
    //Hull algorithm placed here:  		
    Point.x21 = Point.x21 - (tan(Point.p2x)*d_phantom);
    Point.z21 = Point.z21 - (tan(Point.p2z)*d_phantom);
    Point.y21 = Point.y21 - d_phantom;
    
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

    Point.x12 = Point.x12 + (tan(Point.p1x)*d_phantom);
    Point.z12 = Point.z12 + (tan(Point.p1z)*d_phantom);
    Point.y12 = Point.y12 + d_phantom;

		TVector3 X0_hull = TVector3(Point.x12, Point.z12, Point.y12); //End-point (for that particular proton)
		TVector3 dirX0 = TVector3(Point.p1x, Point.p1z, Point.p1y); //final direction
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
			
		  Point.x12 = X0_hull.x(); 
	    Point.z12 = X0_hull.y(); 
	    Point.y12 = X0_hull.z(); 
	                
	    Point.x21 = X2_hull.x();
	    Point.z21 = X2_hull.y();
	    Point.y21 = X2_hull.z();
	       
	    tt->Fill();   
	    
	    distAngleX = Point.p2x - Point.p1x;
	    angularDistX->Fill(Point.x21, Point.z21, distAngleX);    
	    distAngleZ = Point.p2z - Point.p1z;
	    angularDistZ->Fill(Point.x21, Point.z21, distAngleZ);
	    
	    //distDeltaX = Point.x21 - Point.x12;
      //lateralDistX->Fill(Point.x21, Point.z21, distDeltaX);    
      //distDeltaZ = Point.z21 - Point.z12;
      //lateralDistZ->Fill(Point.x21, Point.z21, distDeltaZ); 
	  }
	    
	  else if(X0_hull.z()>50.0 && X2_hull.z()<(-50.0)){//Completely missed the phantom from both sides
	       
	    tt->Fill();
	    
	    distAngleX = Point.p2x - Point.p1x;
	    angularDistX->Fill(Point.x21, Point.z21, distAngleX);    
	    distAngleZ = Point.p2z - Point.p1z;
	    angularDistZ->Fill(Point.x21, Point.z21, distAngleZ);
	    
	    //distDeltaX = Point.x21 - Point.x12;
      //lateralDistX->Fill(Point.x21, Point.z21, distDeltaX);    
      //distDeltaZ = Point.z21 - Point.z12;
      //lateralDistZ->Fill(Point.x21, Point.z21, distDeltaZ); 
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
  tt->SetBranchAddress("scatter_in_h",&Point.scatter_in_h);
  tt->SetBranchAddress("scatter_out_h",&Point.scatter_out_h);
    
  NEntries = tt->GetEntries(); //update entries
  
  cout<<NEntries<<endl;
    
  TH2F* angleMeanX = new TH2F("angleMeanX", "Angular distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* angleMeanZ = new TH2F("angleMeanZ", "Angular distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  //TH2F* latMeanX = new TH2F("latMeanX", "Lateral distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  //TH2F* latMeanZ = new TH2F("latMeanZ", "Lateral distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* WEPLMean = new TH2F("WEPLMean", "WEPL distribution", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
	
  TH2F* angleSigX = new TH2F("angleSigX", "Angular sigma distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  TH2F* angleSigZ = new TH2F("angleSigZ", "Angular sigma distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  //TH2F* latSigX = new TH2F("latSigX", "Lateral sigma distribution in X direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
  //TH2F* latSigZ = new TH2F("latSigZ", "Lateral sigma distribution in Z direction", NbinsX, Xmin, Xmax, NbinsZ, Zmin, Zmax);
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
		  
		  //meanLatX = lateralDistX->ProjectionZ("pz5",binxid,binxid,binzid,binzid)->GetMean();
		  //latMeanX->Fill(x,z,meanLatZ);
		  
		  //sigLatX = lateralDistX->ProjectionZ("pz6",binxid,binxid,binzid,binzid)->GetStdDev();
		  //latSigX->Fill(x,z,sigLatX);
		  
		  //meanLatZ = lateralDistZ->ProjectionZ("pz7",binxid,binxid,binzid,binzid)->GetMean();
		  //latMeanZ->Fill(x,z,meanLatZ);
		  
		  //sigLatZ = lateralDistZ->ProjectionZ("pz8",binxid,binxid,binzid,binzid)->GetStdDev();
		  //latSigZ->Fill(x,z,sigLatZ); 
	  }
  }
  
  angleMeanX->Write("MeanAx_D",TObject::kOverwrite);
  angleSigX->Write("SigAx_D",TObject::kOverwrite);
  angleMeanZ->Write("MeanAz_D",TObject::kOverwrite);
  angleSigZ->Write("SigAz_D",TObject::kOverwrite);
  //latMeanX->Write("MeanLatx_D",TObject::kOverwrite);
  //latSigX->Write("SigLatx_D",TObject::kOverwrite);
  //latMeanZ->Write("MeanLatz_D",TObject::kOverwrite);
  //latSigZ->Write("SigLatz_D",TObject::kOverwrite);

  //Filter the wepl!
  TH1D* WEPLdistf;
  for(int i=0;i<NEntries;i++){//Loop over all protons in the new file/tree to fill distributions
    tt->GetEntry(i); 
            
    int binGlobal = angleMeanX->FindBin(Point.x21_h, Point.z21_h); 
    float angleX = Point.p2x_h - Point.p1x_h;
    float angleZ = Point.p2z_h - Point.p1z_h;
    float deltaX = Point.x21_h - Point.x12_h; 
    float deltaZ = Point.z21_h - Point.z12_h;
 
    if( angleX >= (angleMeanX->GetBinContent(binGlobal) - 2.5*angleSigX->GetBinContent(binGlobal)) && angleX <= (angleMeanX->GetBinContent(binGlobal) + 2.5*angleSigX->GetBinContent(binGlobal)) && 
        angleZ >= (angleMeanZ->GetBinContent(binGlobal) - 2.5*angleSigZ->GetBinContent(binGlobal)) && angleZ <= (angleMeanZ->GetBinContent(binGlobal) + 2.5*angleSigZ->GetBinContent(binGlobal)) 
      //  deltaX >= (latMeanX->GetBinContent(binGlobal) - 2.5*latSigX->GetBinContent(binGlobal)) && deltaX <= (latMeanX->GetBinContent(binGlobal) + 2.5*latSigX->GetBinContent(binGlobal)) && 
      //  deltaZ >= (latMeanZ->GetBinContent(binGlobal) - 2.5*latSigZ->GetBinContent(binGlobal)) && deltaZ <= (latMeanZ->GetBinContent(binGlobal) + 2.5*latSigZ->GetBinContent(binGlobal)) 
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
  WEPLMean->Write("MeanWEPL_D",TObject::kOverwrite);
  WEPLSig->Write("SigWEPL_D",TObject::kOverwrite); 
  
  phaseFile->Close(); 
  
  cout << "Filter initialization all done" << endl;
  return 0; 
}
