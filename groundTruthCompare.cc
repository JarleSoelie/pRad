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
#include "TH2F.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>
using namespace std;

int main(){

//same dimensions as your reconstructed pRad and ground truth	
int recXmin = -102;
int recXmax = 102;
int recZmin = -80;
int recZmax = 80;
int NbinsXrec = 272;
int NbinsZrec = 128;	
  
  TFile* groundTruthFile = new TFile("HeadPhantom.root","update");
  TH2D *groundTruth = (TH2D*)gROOT->FindObject("groundTruth90deg");
  
  
  TFile* imgFile = new TFile("../outputFiles/yourFileNameWithImage", "update");
  TH2D *img = (TH2D*)gROOT->FindObject("FullPath_Single"); //name of reconstructed pRad
      
 
  TH2D* weplDiff = new TH2D("WEPLErrorDistribution", "WEPL errors", NbinsXrec,recXmin,recXmax, NbinsZrec,recZmin,recZmax);
  
  int nx = weplDiff->GetXaxis()->GetNbins();
  int ny = weplDiff->GetYaxis()->GetNbins();
  for (int i=1;i<=nx;i++) {
    for (int j=1;j<=ny;j++) {
      double c1 =groundTruth->GetBinContent(i,j);
      double c2 =img->GetBinContent(i,j);
      weplDiff->SetBinContent(i,j,(c2-c1));
    }  
  }

  weplDiff->Write("",TObject::kOverwrite);
  groundTruthFile->Close();
  imgFile->Close();
  
  return 0; 
}
