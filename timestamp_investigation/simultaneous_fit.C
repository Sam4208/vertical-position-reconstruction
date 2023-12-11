#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include<TApplication.h>
#include <THStack.h>
#include <stdio.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
#include "TImageDump.h"
#include "TImage.h"
#include "TMath.h"
#include "TPoint.h"
#include "TColor.h"
#include "TVirtualPad.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TText.h"
#include "RStipples.h"
#include "TList.h"
//#include <bits/stdc++.h> 
#include <iostream>
#include <chrono>
#include <ctime>   
#include "Math/ProbFunc.h"
#include <TGraph.h>
#include <random>
//#include "/sps/nemo/scratch/spratt/analysis/analysis_library.h"

using namespace std;








// combined fit of two histogram with separate functions

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "TH1.h"
#include "Math/WrappedTF1.h"
#include "HFitInterface.h"



string name_folder(string run_num){
  
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
    std::string folder = directory+run_num;

    return folder;
}


struct GlobalChi2 { 
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,  ROOT::Math::IMultiGenFunction & f2, ROOT::Math::IMultiGenFunction & f3,ROOT::Math::IMultiGenFunction & f4) : 
      fChi2_1(f1), fChi2_2(f2), fChi2_3(f3), fChi2_4(f4) {}

   // parameter vector is first background (in common 1 and 2) and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[4]; 
      p1[0] = par[0]; // mean same for all
      p1[1] = par[1]; // width same for all
      p1[2] = par[2]; // height different
      p1[2] = par[3]; // shift different

      double p2[4]; 
      p2[0] = par[0]; // mean same for all
      p2[1] = par[1]; // width same for all
      p2[2] = par[2]; // height different 
      p2[3] = par[4]; // shift different

      double p3[4]; 
      p2[0] = par[0]; // mean same for all
      p2[1] = par[1]; // width same for all
      p2[2] = par[2]; // height different 
      p2[3] = par[5]; // shift different

      double p4[4]; 
      p2[0] = par[0]; // mean same for all
      p2[1] = par[1]; // width same for all
      p2[2] = par[2]; // height different 
      p2[3] = par[6]; // shift different


      return fChi2_1(p1) + fChi2_2(p2) + fChi2_3(p3) + fChi2_4(p4);
   } 

   const  ROOT::Math::IMultiGenFunction & fChi2_1;
   const  ROOT::Math::IMultiGenFunction & fChi2_2;
   const  ROOT::Math::IMultiGenFunction & fChi2_3;
   const  ROOT::Math::IMultiGenFunction & fChi2_4;
};





double gaussian_fit_function(double* x, double* par){
  double X=x[0]*0.0001;
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double shift=par[3];
  double a = (((shift+1)/X) - m ) / s;
  double y = h* std::exp(-0.5f * a * a);
  return y;

}



int main2() {
  Double_t angle = TMath::Pi()/4.;
  Double_t cosphi = TMath::Cos(angle);
  Double_t sinphi = TMath::Sin(angle);
   
  TH2F *h2 = new TH2F("h2","h2",100,0,0.8,100,0,0.8);
  TRandom r;     
  for (Int_t i=0;i<1000;i++) {
    Double_t x = r.Uniform(0,1);
    Double_t y = r.Gaus(0,0.02);
    Double_t u = cosphi*x -sinphi*y;
    Double_t v = 2*sinphi*x +cosphi*y;
    h2->Fill(u,v);
  }



  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,1);
  f1->SetParameters(0.,1.);
  f1->SetLineColor(kRed);
  h2->Fit(f1);

   

  string run_num = "1051";
  string folder = name_folder(run_num);
  std::string title_root_file = "/simultaneous_fits.root";
  TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
  rootfile->cd();



  TCanvas* simultaneous_fit = new TCanvas("simultaneous_fit");
  simultaneous_fit ->cd();
  h2->Draw();
  f1->Draw("same");
  simultaneous_fit->Write();
  

  rootfile->Close();

  
  

  return 1;
}




int other_main(){

string run_num = "1051";
string folder = name_folder(run_num);
string general_inves = (folder+"/general_inves.root").c_str();

TFile *general_inves_file= new TFile(general_inves.c_str(), "READ");
TH2D* z_vs_prop_time_hist_inner_good_cells_data = (TH2D*)general_inves_file->Get("z_vs_prop_time_hist_inner_good_cells");
TH1D** y_projections = new TH1D*[z_vs_prop_time_hist_inner_good_cells_data->GetNbinsX()];

ROOT::Fit::BinData data; 

for (int i=0;i<z_vs_prop_time_hist_inner_good_cells_data->GetNbinsX();i++){
y_projections[i] = z_vs_prop_time_hist_inner_good_cells_data->ProjectionY(("velocity_1D_hist_"+to_string(i+1)).c_str(),i,i+1);

ROOT::Fit::FillData(data, y_projections[i]); 
}

  cout << "data size is " << data.Size() << endl;

  TF1 * f1 = new TF1("f1",gaussian_fit_function,3000,7000,4);
  f1->SetParameters(pow(5000,-1),0.1*(pow(4000,-1) - pow(6000,-1)),10000,0.1);

  ROOT::Math::WrappedTF1 wf(*f1);

  ROOT::Fit::Fitter fitter;
  fitter.SetFunction(wf);

  fitter.Config().ParSettings(0).SetLimits(2.053*pow(10,-4) ,  2.05400001*pow(10,-4));
  fitter.Config().ParSettings(1).SetLimits(9.92*pow(10,-6)  ,  9.920000001*pow(10,-6));
  fitter.Config().ParSettings(2).SetLimits(0,1000000);
  fitter.Config().ParSettings(3).SetLimits(0,0.01);

  fitter.Fit(data);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  double mean       = result.Parameter(0);
  double sd         = result.Parameter(1);
  double height     = result.Parameter(2);
  double shift      = result.Parameter(3);

  cout<<mean<<"      "<<sd<<"    "<<height<<"     "<<shift<<endl; 

  
  std::string title_root_file = "/simultaneous_fits.root";
  TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
  rootfile->cd();

  TCanvas* simultaneous_fit = new TCanvas("simultaneous_fit");
  simultaneous_fit ->cd();

  y_projections[10]->Draw();
  f1->Draw("Same");
  simultaneous_fit->Update();
  simultaneous_fit->Write();

  rootfile->Close();
  
return 1;
}



double gaussian_fit_function_new(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double phi=par[3];
  double a = ((((1+phi)/(X)) - (m))/s);
  double y = h * std::exp(-0.5f * a * a);
  return y;

}


void fit_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,1/6,1/3);//mean
    fitFunction->SetParLimits(1,0.01, 10);//sd
    fitFunction->SetParLimits(2,0, 1000000);//hight

    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);

    cout<<endl;
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;

}

int main3() { 



string run_num = "1051";
string folder = name_folder(run_num);
string general_inves = (folder+"/general_inves.root").c_str();

TFile *general_inves_file= new TFile(general_inves.c_str(), "READ");
TH2D* z_vs_prop_time_hist_inner_good_cells_data = (TH2D*)general_inves_file->Get("z_vs_prop_time_hist_inner_good_cells");
TH1D** y_projections = new TH1D*[z_vs_prop_time_hist_inner_good_cells_data->GetNbinsX()];
Double_t factor = 1000;

TH1D* shift_vs_z = new TH1D("shift_vs_z","shift_vs_z",100,0,1);

for (int i=0;i<z_vs_prop_time_hist_inner_good_cells_data->GetNbinsX();i++){ 
y_projections[i] = z_vs_prop_time_hist_inner_good_cells_data->ProjectionY(("velocity_1D_hist_"+to_string(i+1)).c_str(),i,i+1);
y_projections[i]->Scale(factor/(y_projections[i]->GetBinContent(y_projections[i]->GetMaximumBin())), "width");

}


  TH1D * h1 = y_projections[20];
  TH1D * h2 = y_projections[25];
  TH1D * h3 = y_projections[30];
  TH1D * h4 = y_projections[35];

  TH1D * h_test = y_projections[50];

  TF1* simple_fit_function = new TF1("simple_fit_function", gaussian_fit_function_new,3000,7000,3);
  fit_gaussian_function(h_test,simple_fit_function);
  double standard_mean = simple_fit_function->GetParameter(0);
  double standard_sd = simple_fit_function->GetParameter(1);
  double standard_height = simple_fit_function->GetParameter(2);


  TF1 * f1 = new TF1("f1",gaussian_fit_function,3000,7000,4);
  //f1->SetParameters(pow(5000,-1),0.1*(pow(4000,-1) - pow(6000,-1)),10000,0.001);

  TF1 * f2 = new TF1("f2",gaussian_fit_function,3000,7000,4);
  //f2->SetParameters(pow(5000,-1),0.1*(pow(4000,-1) - pow(6000,-1)),10000,0.001);

  TF1 * f3 = new TF1("f3",gaussian_fit_function,3000,7000,4);
  //f3->SetParameters(pow(5000,-1),0.1*(pow(4000,-1) - pow(6000,-1)),10000,0.001);

  TF1 * f4 = new TF1("f4",gaussian_fit_function,3000,7000,4);
  //f3->SetParameters(pow(5000,-1),0.1*(pow(4000,-1) - pow(6000,-1)),10000,0.001);

  ROOT::Math::WrappedMultiTF1 wf1(*f1,1);
  ROOT::Math::WrappedMultiTF1 wf2(*f2,1);
  ROOT::Math::WrappedMultiTF1 wf3(*f3,1);
  ROOT::Math::WrappedMultiTF1 wf4(*f4,1);

  ROOT::Fit::DataOptions opt; 

  ROOT::Fit::DataRange range1; 
  range1.SetRange(3000,7000);
  ROOT::Fit::BinData data1(opt,range1); 
  ROOT::Fit::FillData(data1, h1);

  ROOT::Fit::DataRange range2; 
  range2.SetRange(3000,7000);
  ROOT::Fit::BinData data2(opt,range2); 
  ROOT::Fit::FillData(data2, h2);

  ROOT::Fit::DataRange range3; 
  range3.SetRange(3000,7000);
  ROOT::Fit::BinData data3(opt,range3); 
  ROOT::Fit::FillData(data3, h3);

  ROOT::Fit::DataRange range4; 
  range4.SetRange(3000,7000);
  ROOT::Fit::BinData data4(opt,range4); 
  ROOT::Fit::FillData(data4, h4);

  ROOT::Fit::Chi2Function chi2_1(data1, wf1);
  ROOT::Fit::Chi2Function chi2_2(data2, wf2);
  ROOT::Fit::Chi2Function chi2_3(data3, wf3);
  ROOT::Fit::Chi2Function chi2_4(data4, wf4);

  GlobalChi2 globalChi2(chi2_1, chi2_2, chi2_3, chi2_4);

  ROOT::Fit::Fitter fitter;

  //double mean_ll = 2.053*pow(10,-4);      strict mean bounds
  //double mean_ul =  2.05400001*pow(10,-4);

  double mean_ll = standard_mean - standard_mean*0.3;
  double mean_ul =  standard_mean + standard_mean*0.3;
  double sd_ll = standard_sd - standard_sd*0.3;
  double sd_ul = standard_sd + standard_sd*0.3;
  double h_ll= standard_height - standard_height*0.3;
  double h_ul= standard_height + standard_height*0.3;
  double shift_ll = 0 ;
  double shift_ul = 0.05;

  const int Npar = 7; 
  double par0[Npar] = {standard_mean,standard_sd,standard_height,0,0,0,0};


  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(7,par0);

  //fitter.Config().ParSettings(0).SetLimits(0,100000);//mean
  //fitter.Config().ParSettings(1).SetLimits(0,10000);//sd
  //fitter.Config().ParSettings(2).SetLimits(0,10000);//height1

  fitter.Config().ParSettings(0).SetLimits(mean_ll,mean_ul);//mean
  fitter.Config().ParSettings(1).SetLimits(sd_ll,sd_ul);//sd
  fitter.Config().ParSettings(2).SetLimits(h_ll,h_ul);//height1
  fitter.Config().ParSettings(3).SetLimits(shift_ll,shift_ul);//shift1
  fitter.Config().ParSettings(4).SetLimits(shift_ll,shift_ul);//shift2
  fitter.Config().ParSettings(5).SetLimits(shift_ll,shift_ul);//shift3
  fitter.Config().ParSettings(6).SetLimits(shift_ll,shift_ul);//shift3
  //fitter.Config().ParSettings(5).Fix();

  fitter.FitFCN(7,globalChi2,par0,data1.Size()+data2.Size()+data3.Size());
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);


  double mean        = result.Parameter(0);
  double sd          = result.Parameter(1);
  double height      = result.Parameter(2);
  double shift1      = result.Parameter(3);
  double shift2      = result.Parameter(4);
  double shift3      = result.Parameter(5);
  double shift4      = result.Parameter(6);



  cout<<"Shared Mean:  "<<mean<<endl; 
  cout<<"Shared SD:    "<<sd<<endl; 
  cout<<"height:       "<<height<<endl; 
  cout<<"shift1:       "<<shift1<<endl; 
  cout<<"shift2:       "<<shift2<<endl; 
  cout<<"shift3:       "<<shift3<<endl; 
  cout<<"shift3:       "<<shift4<<endl; 

  //shift_vs_z->

  std::string title_root_file = "/simultaneous_fits.root";
  TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
  rootfile->cd();


  TCanvas * c1 = new TCanvas("Simultaneous fit");
  c1->Divide(1,5);
  c1->cd(1);
  gStyle->SetOptFit(1111);

  f1->SetParameters(result.Parameter(0), result.Parameter(1),result.Parameter(2),result.Parameter(3));
  h1->Draw(); 
  f1->Draw("Same");
  c1->Update();


  c1->cd(2);
  f2->SetParameters(result.Parameter(0), result.Parameter(1),result.Parameter(2),result.Parameter(4));
  h2->Draw(); 
  f2->Draw("Same");


  c1->cd(3);
  f3->SetParameters(result.Parameter(0), result.Parameter(1),result.Parameter(2),result.Parameter(5));
  h3->Draw(); 
  f3->Draw("Same");

  c1->cd(4);
  f4->SetParameters(result.Parameter(0), result.Parameter(1),result.Parameter(2),result.Parameter(6));
  h4->Draw(); 
  f4->Draw("Same");

  c1->cd(5);
  simple_fit_function->Draw();
  h_test->Draw();

  c1->Update();
  c1->Write();

  f1->Write();
  f2->Write();
  f3->Write();

  rootfile->Close();
  
  return 1;



}





void run_fit_function(TH2D* hist2d, TF2* fit_function){

    //fit_function->SetParLimits(0,2.17*pow(10,-4),2.4300001*pow(10,-4));//mean
    //fit_function->SetParLimits(1,7*pow(10,-6), 1.3*pow(10,-5));//sd
    //fit_function->SetParLimits(2,0, 1000000);//hight
    //fit_function->SetParLimits(3,0, 1);//phi


    fit_function->SetParLimits(0,0,100);//mean
    fit_function->SetParLimits(1,0,100);//sd
    fit_function->SetParLimits(2,0,1000000);//hight
    fit_function->SetParLimits(3,0, 0.00001);//phi

  hist2d->Fit(fit_function); // "R" for range-based fit

  double standard_mean = fit_function->GetParameter(0);
  double standard_sd = fit_function->GetParameter(1);
  double standard_height = fit_function->GetParameter(2);
  double shift_parameter = fit_function->GetParameter(3);

  double mean_ll = standard_mean - standard_mean*0.3;
  double mean_ul =  standard_mean + standard_mean*0.3;
  double sd_ll = standard_sd - standard_sd*0.3;
  double sd_ul = standard_sd + standard_sd*0.3;
  double h_ll= standard_height - standard_height*0.3;
  double h_ul= standard_height + standard_height*0.3;
  double shift_ll = 0 ;
  double shift_ul = 0.05;


}





double fit_function_equation(double* x, double* par){

  double z=x[0];
  double tp=x[1]*0.001;

  double mean =par[0];
  double sd = par[1];
  double height = par[2];
  double shift = par[3];

  double distance_term = (pow((z+1)*0.5,2)+pow((1-z)*0.5,2));

  double a = (((1+(shift*distance_term))/tp) - mean ) / sd;
  double y = height* std::exp(-0.5f * a * a);
  return y;


}




void combine_runs(int run_num, TH2D** prop_vs_z_2D_hist_array,TH1D** prop_vs_z_1D_hist_output_array,TH1D** prop_vs_z_1D_hist_counter_array){

cout<<"Combining Run "<<run_num<<endl;

string folder = name_folder(to_string(run_num).c_str());
string general_inves = (folder+"/general_inves.root").c_str();
string prop_time_vs_z_individual_cells = (folder+"/prop_time_vs_z_individual_cells.root").c_str();

TFile *prop_time_vs_z_individual_cells_file= new TFile(prop_time_vs_z_individual_cells.c_str(), "READ");
int no_cells=2034;
TH1D** average_hists = new TH1D*[no_cells];
TH1D** counter_hists = new TH1D*[no_cells];
TH2D** heatmap_hists = new TH2D*[no_cells];
for (int i=0;i<no_cells;i++){
average_hists[i] =  (TH1D*)prop_time_vs_z_individual_cells_file->Get(("output"+to_string(i)).c_str());
counter_hists[i] =  (TH1D*)prop_time_vs_z_individual_cells_file->Get(("counter"+to_string(i)).c_str());
heatmap_hists[i] =  (TH2D*)prop_time_vs_z_individual_cells_file->Get(("full_2d_heatmap_prop_vs_z"+to_string(i)).c_str());

prop_vs_z_2D_hist_array[i]->Add(heatmap_hists[i]);
prop_vs_z_1D_hist_output_array[i]->Add(counter_hists[i]);
}
}




TH1D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound){
    
      int no_cells=2034;
      TH1D** h = new TH1D*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TH1D(name.c_str(),name.c_str(), bins, lower_bound, upper_bound);
                        }

      return h;
}




TH2D** make_array_of_hist_2d(string name_of_hist, double bins,double lower_bound,double upper_bound,double bins2,double lower_bound2,double upper_bound2){
      
      int no_cells=2034;
      TH2D** h = new TH2D*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TH2D(name.c_str(),name.c_str(), bins, lower_bound, upper_bound ,bins2, lower_bound2, upper_bound2);
                        }

      return h;
}






int main(){




string run_num = "1051";
string folder = name_folder(run_num);

TH2D** heatmap = make_array_of_hist_2d("prop_vs_z_2D_hist_array",100,-1,1,1000,3000,7000);
TH1D** output =  make_array_of_hist("prop_vs_z_1D_hist_output_array",100,-1,1);
TH1D** counter = make_array_of_hist("prop_vs_z_1D_hist_counter_array",100,-1,1);

combine_runs(1051,heatmap,output,counter);
//combine_runs(1054,heatmap,output,counter);
//combine_runs(1058,heatmap,output,counter);
//combine_runs(1062,heatmap,output,counter);
//combine_runs(1066,heatmap,output,counter);

TH2D* z_vs_prop_time_hist_inner_good_cells_data = heatmap[100];
TH1D* hist_1d = heatmap[100]->ProjectionY(("velocity_1D_hist_"+to_string(100)).c_str(),20,21);
TF1*  fit_1d = new TF1(("fitFunction_gauss_with_shift"+to_string(20)).c_str(), gaussian_fit_function_new,4000, 5000,4);
z_vs_prop_time_hist_inner_good_cells_data->Rebin2D(1,10);




TF2* fit_function = new TF2("fit_function",fit_function_equation,-0.8,0.8,4000,5000,4);
run_fit_function(z_vs_prop_time_hist_inner_good_cells_data,fit_function);

fit_1d->SetParameters(fit_function->GetParameter(0), fit_function->GetParameter(1), fit_function->GetParameter(2), fit_function->GetParameter(3)*(pow((1-0.2),2)+pow((1+0.2),2)));


  std::string title_root_file = "/simultaneous_fits.root";
  TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
  rootfile->cd();
  

  TCanvas * c1 = new TCanvas("Simultaneous fit");
  c1->cd();
  hist_1d->Draw();
  fit_1d->Draw("Same");
  c1->Write();
  fit_1d->Write();

  z_vs_prop_time_hist_inner_good_cells_data->Write();
  fit_function->Write();



  rootfile->Close();
  
  return 1;



}

























#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TRandom3.h>

// Define the Gaussian function
Double_t gaussian(Double_t* x, Double_t* par) {
    // par[0]: amplitude, par[1]: mean, par[2]: width
    return par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / (2 * par[2] * par[2]));
}

// Function to set parameters for simultaneous fit
void setParameters(Double_t* params, Int_t i, Double_t sharedAmplitude, Double_t sharedWidth) {
    params[0] = sharedAmplitude;  // Shared amplitude
    params[1] = i;                // Individual mean
    params[2] = sharedWidth;      // Shared width
}


int main_extra() {



  string run_num = "1051";
string folder = name_folder(run_num);

    const Int_t numHistograms = 10;
    const Int_t nBins = 100;

    // Create example histograms
    TH1D* histograms[numHistograms];
    TF1* globalFitFunction = new TF1("globalFitFunction", gaussian, -5, 15, 3 * numHistograms);  // One fit function for all histograms

    for (Int_t i = 0; i < numHistograms; ++i) {
        histograms[i] = new TH1D(Form("hist%d", i), Form("Histogram %d", i), nBins, -5, 15);

TRandom3 randomGenerator;
         for (Int_t j = 0; j < 100; ++j) {
        Double_t value = randomGenerator.Gaus(10, 2);
        histograms[i]->Fill(value);
    }
        // Fill histograms with your data...
    }

    // Set initial parameters for the entire fit
    Double_t params[3 * numHistograms];
    for (Int_t i = 0; i < numHistograms; ++i) {
        setParameters(params + 3 * i, i, 1.0, 1.0);  // Initial values for shared amplitude and width
    }
    globalFitFunction->SetParameters(params);

    // Fit all histograms simultaneously with the same fit function
    for (Int_t i = 0; i < numHistograms; ++i) {
        // Set individual parameters for each histogram
        setParameters(params + 3 * i, i, globalFitFunction->GetParameter(3 * i), globalFitFunction->GetParameter(3 * i + 2));
    }
    globalFitFunction->SetParameters(params);

    // Fit all histograms simultaneously
    for (Int_t i = 0; i < numHistograms; ++i) {
        histograms[i]->Fit(globalFitFunction, "Q"); // "Q" option for quiet mode
    }

    // Access fit results and parameters
    std::cout << "Shared Parameters:" << std::endl;
    for (Int_t i = 0; i < numHistograms; ++i) {
        std::cout << "Histogram " << i << " - Amp: " << globalFitFunction->GetParameter(3 * i) << ", Wid: " << globalFitFunction->GetParameter(3 * i + 2) << std::endl;
    }




  std::string title_root_file = "/simultaneous_fits.root";
  TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
  rootfile->cd();
    // Optionally, you can create a canvas to visualize the fits
    TCanvas* canvas = new TCanvas("canvas", "Fits", 800, 600);
    canvas->Divide(2, 5);

    for (Int_t i = 0; i < numHistograms; ++i) {
        canvas->cd(i + 1);
        histograms[i]->Draw();
        globalFitFunction->Draw("same");
    }

    canvas->Update();
    canvas->Write();

     rootfile->Close();


    return 0;
}
