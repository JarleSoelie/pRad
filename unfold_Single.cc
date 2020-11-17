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

double d_source = 660.0; //*** cm (6.6m distance to scanning magnets/source)
double s_pos = 0.3; //*** cm beam spot size at isocenter (from simulations)
double s_angle = 0.003; //*** rad This should be the uncertainty on the TPS angle  

double s_pos_out = 0.0005; //*** cm Position unc. in tracking layers, 0 = ideal detector, 0.0005 = Bergen detector

double d_entry = 0; //*** cm OBS If not using manual hull, change accordingly (=0.0 if manual hull is used or found)
double d_exit = 0; //*** cm
double d_T = 5.0; //*** cm Distance between trackers 

int Nsteps = 512;

//For defining the reconstruction area, in this case 0.75x1.25 mm pixel size, but change to whatever you need:
int NbinsX = 272; //***
int NbinsZ = 128; //***
float  Xmin = -102; //*** mm 
float  Xmax = 102; //*** mm 
float  Zmin = -80; //*** mm 
float  Zmax = 80; //*** mm 

#define IMAGENAME "FullPath_Single" //*** name of image/radiography

struct Proton{
  Float_t Einit;
  
  Float_t TPSx0,TPSy0,TPSz0;
  Float_t TPSpx0,TPSpy0,TPSpz0;
  
  Float_t x11,y11,z11;
  Float_t x12,y12,z12;
  Float_t p1x,p1y,p1z;
  
  Float_t x21,y21,z21;
  Float_t x22,y22,z22;
  Float_t p2x,p2y,p2z;
  
  Float_t scatter_out;
  Float_t weplReal;
};

void ComputeMLP(Proton *, TProfile2D *, float);
double Sigmat1(double);
double Sigmay1(double);
double Sigmaty1(double);
double Sigmat2(double,double);
double Sigmay2(double,double);
double Sigmaty2(double,double);

int main(int argc, char *argv[]){
  // Load Proton data
  Proton Point;
  
  //Prepare in file from createTree	
  char rootf[10] = ".root";
  char inf[50] = "../outputFiles/yourFileName_";//***
  strcat(inf,argv[1]);
  strcat(inf,rootf);   
  TFile* phaseFile = new TFile(inf,"update"); //Associate any new objects to the proper working file
  
  TTree* tt = (TTree*)phaseFile->Get("filteredTree"); 

  TProfile2D* project2D = new TProfile2D("Title of image",IMAGENAME,NbinsX,Xmin,Xmax, NbinsZ,Zmin,Zmax, "s");
 
  project2D->GetXaxis()->SetTitle("Lateral position [mm]");
  project2D->GetYaxis()->SetTitle("Vertical position [mm]");  
  project2D->GetZaxis()->SetTitle("WEPL [mm]");
   
  tt->SetBranchAddress("x12_h",&Point.x12);
  tt->SetBranchAddress("y12_h",&Point.y12);
  tt->SetBranchAddress("z12_h",&Point.z12);
  
  tt->SetBranchAddress("p1x_h",&Point.p1x);
  tt->SetBranchAddress("p1y_h",&Point.p1y);
  tt->SetBranchAddress("p1z_h",&Point.p1z);
  
  tt->SetBranchAddress("x21_h",&Point.x21);
  tt->SetBranchAddress("y21_h",&Point.y21);
  tt->SetBranchAddress("z21_h",&Point.z21);
  
  tt->SetBranchAddress("p2x_h",&Point.p2x);
  tt->SetBranchAddress("p2y_h",&Point.p2y);
  tt->SetBranchAddress("p2z_h",&Point.p2z); 
  
  tt->SetBranchAddress("weplReal_h",&Point.weplReal);
  tt->SetBranchAddress("scatter_out_h",&Point.scatter_out);
  
  int NEntries = tt->GetEntries(); 
 
  //Find the filter maps from filter.cc
  TH2F* angleMeanX = ((TH2F *)(gROOT->FindObject("MeanAx"))); 
  TH2F* angleMeanZ = ((TH2F *)(gROOT->FindObject("MeanAz")));
  TH2F* WEPLMean = ((TH2F *)(gROOT->FindObject("MeanWEPL")));
	
  TH2F* angleSigX = ((TH2F *)(gROOT->FindObject("SigAx")));
  TH2F* angleSigZ = ((TH2F *)(gROOT->FindObject("SigAz")));
  TH2F* WEPLSig = ((TH2F *)(gROOT->FindObject("SigWEPL")));
  TH2F* WEPLfwhm = ((TH2F *)(gROOT->FindObject("FWHMWEPL")));

  if(angleMeanX && angleSigX){cout<<"Angle in x filter found"<<endl;}
  if(angleMeanZ && angleSigZ){cout<<"Angle in z filter found"<<endl;}
  if(WEPLMean && WEPLSig){cout<<"WEPL filter found"<<endl;}

  cout << "Beginning reconstructing/unfolding" << endl;
 
  int NPerformed = 0;
  int binGlobal;
  double angleX, angleZ, WEPLn;
  float s_pos;
 
  for(int i=0;i<NEntries;i++){//Loop over all protons again
	if(i%100000 == 0) cout<<i<<" "<<NPerformed<<endl;
    tt->GetEntry(i);  
     
    binGlobal = angleMeanX->FindBin(Point.x21, Point.z21); //Exact same bin for all the other histograms and binGlobal = angleMeanX->GetBin((angleMeanX->GetXaxis()->FindBin(Point.z21)), (angleMeanX->GetYaxis()->FindBin(Point.x21))); 
    angleX = Point.p2x - Point.p1x;
    angleZ = Point.p2z - Point.p1z;
    WEPLn = Point.weplReal;
    //I prefer using the fwhm for the wepl filter instead of the 2.5 sigma   
    if( WEPLn >= (WEPLMean->GetBinContent(binGlobal) - 1*WEPLfwhm->GetBinContent(binGlobal)) && WEPLn <= (WEPLMean->GetBinContent(binGlobal) + 1*WEPLfwhm->GetBinContent(binGlobal)) &&
	//WEPLn >= (WEPLMean->GetBinContent(binGlobal) - 2.5*WEPLSig->GetBinContent(binGlobal)) && WEPLn <= (WEPLMean->GetBinContent(binGlobal) + 2.5*WEPLSig->GetBinContent(binGlobal)) &&
        angleX >= (angleMeanX->GetBinContent(binGlobal) - 2.5*angleSigX->GetBinContent(binGlobal)) && angleX <= (angleMeanX->GetBinContent(binGlobal) + 2.5*angleSigX->GetBinContent(binGlobal)) && 
        angleZ >= (angleMeanZ->GetBinContent(binGlobal) - 2.5*angleSigZ->GetBinContent(binGlobal)) && angleZ <= (angleMeanZ->GetBinContent(binGlobal) + 2.5*angleSigZ->GetBinContent(binGlobal))
      ){  
	      ComputeMLP(&Point, project2D, s_pos); 
		    NPerformed++;
    }   
  }  
   
  project2D->Write(IMAGENAME,TObject::kOverwrite);

  phaseFile->Close();
  cout << "Radiography ready" << endl;  
  return 0; 
}

void ComputeMLP(Proton *Point, TProfile2D* project2D, float s_pos){	
	  //---------------------------------------------------------------
    //transform all mm into cm 
    //---------------------------------------------------------------
    Point->x12/=10.; Point->y12/=10.; Point->z12/=10.;
    Point->x21/=10.; Point->y21/=10.; Point->z21/=10.; 
    //------------------------------------------------------------  	

    TVector3 m0(Point->x12,Point->y12,Point->z12);
    TVector3 p0(Point->p1x, Point->p1y, Point->p1z);
    
    TVector3 m1(Point->x21, Point->y21, Point->z21);
    TVector3 p1(Point->p2x, Point->p2y, Point->p2z);

        //Initialization of the MLP
    float X_mlp,Z_mlp, theta_X_mlp,theta_Z_mlp, X_mlp_sigma,Z_mlp_sigma;
    TVector3 p;
    TVector3 p_old = m0;
    
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
      T_out[0] = 1; //Eq. 26 in Krah
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
    
    double S_in[4]; //Was not included before, eq. 14 in Krah
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

    //-------------Initialize the image reconstruction variables!-------
    std::map<pair<int,int>,double> Lengthmap; // Image grid that will contain all path information (length spent in each image column)
    std::pair<std::map<std::pair<int,int>,double>::iterator,bool> ret;
    int binx,biny,binz;
    double TotL = 0;
    //------------------------------------
    //Parameter initialization
    double sy1, sy2, st1, st2, sty1, sty2;
    double determinant_1, determinant_2, determinant_C12;
    
    //Initialize the MLP iterators
    double step_length = (m1.y()-m0.y())/Nsteps;
    double posy=m0.y()+step_length; //skip the first track value since not defined here
    
    //Step through until reach the exit depth
    while(posy<m1.y()){  //y is in beam direction
                    
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
      X_mlp_sigma = second_first[2] * C2[1] + second_first[3] * C2[3];
        
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
      Z_mlp_sigma = second_first[2] * C2[1] + second_first[3] * C2[3];
   
    // This is where the reconstruction is happening where the proton length spent inside each channel is accounted for and weights the proton wepl accordingly
    p = TVector3(X_mlp, posy, Z_mlp);
      
    binx = project2D->GetXaxis()->FindBin(p.x()*10.0); //find the image bin we are in in the two traverse directions z and x
    binz = project2D->GetYaxis()->FindBin(p.z()*10.0);  
    
    double L = TVector3(p-p_old).Mag(); // This is the length (magnitude) of the proton path inside the step  (/10 -> mm-cm  conversion)         
    pair<int,int> bin2dID = pair<int,int>(binx,binz); //find the key to the bin positions the proton is in (bin2dID) 
    
    ret = Lengthmap.insert(pair<pair<int,int>,double>(bin2dID,L)); //checks using insert whether the bin2dID key (column) has changed (proton has crossed it) and only connects L to the bin pos if it has
    if ( ret.second==false ) Lengthmap[bin2dID] += L; //if (false=bin2dID is the same), column corresponding to bin2dID has not been crossed by the proton investigated, L is added to the length traversed so far
        
    TotL+=L; //total length, sum up   
    p_old = p;
    posy+=step_length;
  }

  std::map<std::pair<int,int>,double>::iterator it;
 
  //This is where the image is filled up and weighted!
  for(it = Lengthmap.begin(); it != Lengthmap.end(); it++) { //Goes through the Lengthmap
    int BinX = it->first.first;
    int BinZ = it->first.second;
    double L = it->second; 
    double x = project2D->GetXaxis()->GetBinCenter(BinX);
    double z = project2D->GetYaxis()->GetBinCenter(BinZ);
    project2D->Fill(x,z,Point->weplReal, pow((L/TotL),2)); 
  }     
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
