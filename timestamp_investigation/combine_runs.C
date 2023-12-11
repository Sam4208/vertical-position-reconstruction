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


//////////////////////////////////////

double calc_eps(double d,double vel,double acc){

double L=1;


double term1 = ((2*d)-L)*pow(vel,2);
double term2 = ((2*acc)*(pow(d,2)-pow((L-d),2)));
double term3 = (L*pow(vel,2));
double term4 = ((2*acc)*(pow(d,2)+pow((L-d),2)));

double eps = (term1+term2)/(term3+term4);

//cout<<term1<<"    "<<term2<<"    "<<term3<<"    "<<term4<<endl;

return eps;

}



double calc_tp(double d,double vel,double acc){

double L=1;

double tp = (L/vel) + 2*(acc/pow(vel,3))*(pow(d,2)+pow((L-d),2));

return tp;

}




bool check_cell_on_edge(int i){

bool output = false;

int side_;
int column_;
int layer_;

int cell_num;

for (int side=0;side<2;side++){
    for (int column=0;column<113;column++){
        for (int layer=0;layer<9;layer++){
        cell_num =  side*9*113 + column*9 + layer;

//cout<<cell_num<<endl;

if (cell_num == i){
side_=side;
column_=column;
layer_=layer;
}


}}}

if (layer_ == 0 || layer_== 8 ){
output=true;
}
if (column_==0 || column_==112 ){
output=true;
}

return output;

}




bool is_cell_new(double side, double column, double layer){

      std::vector<double> bad_cell_column={9,106,99,53,82,102,99,26};
      std::vector<double> bad_cell_layer={7,2,7,1,1,2,8,3};
      std::vector<double> bad_cell_side={1,0,0,0,1,0,0,0};
      bool output=false;
   
      for (int i=0; i< bad_cell_layer.size(); i++){
         
      if (side==bad_cell_side.at(i)  &&  layer==bad_cell_layer.at(i)   &&  column  == bad_cell_column.at(i)){

         output=true;}}
       
      return output;

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


bool is_cell_other(double side, double column, double layer){

      std::vector<double> bad_cell_column={56,56,56,56};
      std::vector<double> bad_cell_layer={4,5,4,5};
      std::vector<double> bad_cell_side={0,0,1,1};
      bool output=false;
   
      for (int i=0; i< bad_cell_layer.size(); i++){
         
      if (side==bad_cell_side.at(i)  &&  layer==bad_cell_layer.at(i)   &&  column  == bad_cell_column.at(i)){

         output=true;}}
       
      return output;

}


bool is_cell_noisy(double side, double column, double layer){

      std::vector<double> bad_cell_column={1,2,4};
      std::vector<double> bad_cell_layer={1,1,3};
      std::vector<double> bad_cell_side={0,0,0};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}



bool is_cell_HV_trip(double side, double column, double layer){

      std::vector<double> bad_cell_column={3,9,21,84,32,70,75,73,80,100};
      std::vector<double> bad_cell_layer={0,0,0,0,0,0,4,2,2,2};
      std::vector<double> bad_cell_side={0,0,0,0,1,0,0,1,1,1};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}
         

         
         }


      return output;
}



bool is_cell_internal_cabling(double side, double column, double layer){

      std::vector<double> bad_cell_column={86,87,63,57,56,47,79,84,91,99,107,110};
      std::vector<double> bad_cell_layer={8,0,8,0,3,8,6,6,2,8,8,0};
      std::vector<double> bad_cell_side={0,0,0,0,0,1,1,1,1,1,1,1};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}



bool is_cell_broken_anode(double side, double column, double layer){

      std::vector<double> bad_cell_column={11,77,89,101};
      std::vector<double> bad_cell_layer={0,8,1,3};
      std::vector<double> bad_cell_side={0,0,0,0};
      bool output=false;

      for (int i=0; i< bad_cell_layer.size(); i++){
      if (side==bad_cell_side[i]  &&  layer==bad_cell_layer[i]   &&  column  == bad_cell_column[i]  ){
         output=true;}}

      return output;
}




bool is_cell_bad(double side, double column, double layer){
      bool output = false;

      bool broken_anode = is_cell_broken_anode(side, column, layer);
      bool internal_cab = is_cell_internal_cabling(side, column, layer);
      bool HV_trip      = is_cell_HV_trip(side, column, layer);
      bool noisy        = is_cell_noisy(side, column, layer);
      bool other        = is_cell_other(side, column, layer);
      bool new_         = is_cell_new(side, column, layer);

      if (broken_anode == true  ||  internal_cab  == true    ||    HV_trip == true     ||   noisy == true    ||      other == true ||  new_ == true){
         output = true;
      }

      
      return output;
}




bool is_cell_bad(int i){
      bool output = false;
      int side_;
      int column_;
      int layer_;
      int cell_num;

      for (int side=0;side<2;side++){
         for (int column=0;column<113;column++){
            for (int layer=0;layer<9;layer++){
            cell_num =  side*9*113 + column*9 + layer;


      if (cell_num == i){
      side_=side;
      column_=column;
      layer_=layer;
      }

}}}
      int side=side_;
      int column=column_;
      int layer=layer_;

      bool broken_anode = is_cell_broken_anode(side, column, layer);
      bool internal_cab = is_cell_internal_cabling(side, column, layer);
      bool HV_trip      = is_cell_HV_trip(side, column, layer);
      bool noisy        = is_cell_noisy(side, column, layer);
      bool other        = is_cell_other(side, column, layer);
      bool new_         = is_cell_new(side, column, layer);

      if (broken_anode == true  ||  internal_cab  == true    ||    HV_trip == true     ||   noisy == true    ||      other == true ||  new_ == true){
         output = true;
      }

      return output;
}




bool is_next_to_bad_cell(int i){

 bool output = false;
      int side_;
      int column_;
      int layer_;
      int cell_num;

      for (int side=0;side<2;side++){
         for (int column=0;column<113;column++){
            for (int layer=0;layer<9;layer++){
            cell_num =  side*9*113 + column*9 + layer;


      if (cell_num == i){
      side_=side;
      column_=column;
      layer_=layer;
      }

}}}
      int side=side_;
      int column=column_;
      int layer=layer_;


         if (is_cell_bad(side, column+1, layer) == true ||  is_cell_bad(side, column-1, layer)  ==true  ||   is_cell_bad(side, column, layer+1) ==true   ||      is_cell_bad(side, column, layer-1)==true ){

               output = true;

         }
   return output;

}




bool is_next_to_bad_cell(double side, double column, double layer){

   bool output = false;
         if (is_cell_bad(side, column+1, layer) == true ||  is_cell_bad(side, column-1, layer)  ==true  ||   is_cell_bad(side, column, layer+1) ==true   ||      is_cell_bad(side, column, layer-1)==true ){

               output = true;

         }
   return output;

}





string name_folder(string run_num){
  
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
    std::string folder = directory+run_num;

    return folder;
}






///////////////////////////////////







void average_over_hist(TH1D* h_full,TH1D* h_counter,TH1D* h_output){
   
   double number_bins = h_full->GetNbinsX();

	for(int i=0;i<number_bins;i++){

	double total = h_full->GetBinContent(i);
	double no = h_counter->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	h_output->SetBinContent(i,set_value);
	}}}





double get_tp_given_eps(Double_t *x,Double_t *par) {
    double L = 1;
    double eps = x[0]-0.5;

    double vel = par[0];
    double acc = par[1];

    double best_d;


    for (int i=0;i<=10000;i++){
    double d = i*0.0001;
    double eps_test_value = calc_eps(d,vel,acc);
    if ((eps_test_value-eps)>0){
      best_d = d;
     // cout<<d<<endl;
      break;
    }
    }


    double tp = calc_tp(best_d,vel,acc);

    return tp;
   }




void fit_iteration_eps_vs_z(TGraph *graph1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,0, 0.1);
    fitFunction->SetParLimits(1,0, pow(10,-8));

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;

}



void fit_gaussian_function_with_shift(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,2.053*pow(10,-4),2.05400001*pow(10,-4));//mean
    fitFunction->SetParLimits(1,9.92*pow(10,-6), 9.920000001*pow(10,-6));//sd
    fitFunction->SetParLimits(2,0, 1000000);//hight
    fitFunction->SetParLimits(3,-10, 10);//phi

    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    double par_3 = fitFunction->GetParameter(3);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;
    std::cout << "shift : " << par_3 << std::endl;

}

void fit_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,1*pow(10,-4),3*pow(10,-4));//mean
    fitFunction->SetParLimits(1,0, 1);//sd
    fitFunction->SetParLimits(2,0, 1000000);//hight


    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;

}



void fit_simple_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,2.4000,6000);//mean
    fitFunction->SetParLimits(1,50, 1000);//sd
    fitFunction->SetParLimits(2,0, 1000000);//hight

    hist1->Fit(fitFunction,"Q","0"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "mean: " << par_0 << std::endl;
    std::cout << "sd: " << par_1 << std::endl;
    std::cout << "height_constant: " << par_2 << std::endl;

}



void fit_double_gaussian_function(TH1D *hist1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,4000,6000);
    fitFunction->SetParLimits(1,10, 1000);
    fitFunction->SetParLimits(2,10, 10000000);

    fitFunction->SetParLimits(3,4000,6000);
    fitFunction->SetParLimits(4,10, 1000);
    fitFunction->SetParLimits(5,10, 10000000);

    hist1->Fit(fitFunction,"Q"); // "R" for range-based fit

}



double calc_average_of_hist(TH1D *hist,int bin_lower, int bin_higher){
  double sum=0;
  double no_bins =  hist->GetNbinsX();
  for (int i=bin_lower;i<bin_higher;i++){
  sum = sum + hist->GetBinContent(i);
  }
  double average = sum/(bin_higher-bin_lower);
  
return average;
}



void set_hist_limits(TH1D* hist, double lower, double upper){

hist->SetMaximum(8000);
hist->SetMinimum(4000);


}


void add_to_hist(TH1D* hist_average, TH1D* hist_counter, TH1D* hist_full_output,TH1D* hist_full_counter){

  int no_bins = hist_average->GetNbinsX();
  double average;
  double count;
  double input;
  for (int i=0;i<no_bins;i++){
       average = hist_average->GetBinContent(i);
       count = hist_counter->GetBinContent(i);
       input = average*count;
       hist_full_output->AddBinContent(i,input);
       hist_full_counter->AddBinContent(i,count);
  }}



TH1D** generate_full_distribution(TH1D** average_hists,TH1D** counter_hists){
int no_cells = 2034;
TH1D** output_hist_array = new TH1D*[4];

TH1D* bad = new TH1D("bad","bad", 100,0,1);
TH1D* bad_counter = new TH1D("bad_counter","bad_counter", 100,0,1);
TH1D* bad_output = new TH1D("bad_output","bad_output", 100,0,1);

TH1D* edge = new TH1D("edge","edge", 100,0,1);
TH1D* edge_counter = new TH1D("edge_counter","edge_counter", 100,0,1);
TH1D* edge_output = new TH1D("edge_output","edge_output", 100,0,1);

TH1D* inner = new TH1D("inner","inner", 100,0,1);
TH1D* inner_counter = new TH1D("inner_counter","inner_counter", 100,0,1);
TH1D* inner_output = new TH1D("inner_output","inner_output", 100,0,1);

TH1D* near = new TH1D("near","near", 100,0,1);
TH1D* near_counter = new TH1D("near_counter","near_counter", 100,0,1);
TH1D* near_output = new TH1D("near_output","near_output", 100,0,1);

    for (int i=0; i<no_cells;i++){
        
        if (is_cell_bad(i)==true){
         add_to_hist(average_hists[i],counter_hists[i],bad,bad_counter);
        }

        if (is_next_to_bad_cell(i)==true && is_cell_bad(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],near,near_counter);
        }

        if (check_cell_on_edge(i)==false  && is_cell_bad(i)==false && is_next_to_bad_cell(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],inner,inner_counter);
        }

        if (check_cell_on_edge(i)==true  && is_cell_bad(i)==false && is_next_to_bad_cell(i)==false){
         add_to_hist(average_hists[i],counter_hists[i],edge,edge_counter);
        }}

       // set_hist_limits(inner_output,4000,8000);
       // set_hist_limits(edge_output,4000,8000);
       // set_hist_limits(near_output,4000,8000);
       // set_hist_limits(bad_output,4000,8000);

        average_over_hist(inner,inner_counter,inner_output);
        average_over_hist(edge,edge_counter,edge_output);
        average_over_hist(near,near_counter,near_output);
        average_over_hist(bad,bad_counter,bad_counter);

        output_hist_array[0]=inner_output;
        output_hist_array[1]=edge_output;
        output_hist_array[2]=near_output;
        output_hist_array[3]=bad_output;
        return output_hist_array;
        }




double pick_from_norm_distribution(double v_central, double width){
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::normal_distribution<double> distribution(v_central,width);
  double vel = distribution(generator);
  return vel;
}



double pick_random_height(){
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> distribution(150,850);
  double height = distribution(generator);
  height = height*0.001;
  return height;
}



void fill_model_2d_hist(TH2D* model_prop_vs_z_hist,double acc_dis_width){

double v_central = 2.17*pow(10,-4);
double v_width = 1*pow(10,-5);
double height = pick_random_height();
double vel    = pick_from_norm_distribution(v_central,v_width);
double k    = pick_from_norm_distribution(5.10710,5.10710*acc_dis_width);
double acc = k*pow(10,-6)*vel;
double eps = calc_eps(height,vel,acc);
double tp =calc_tp(height,vel,acc);
model_prop_vs_z_hist->Fill(eps,tp,1);

//cout<<"  height: "<<height<<"  vel: "<<vel<<"  eps: "<<eps<<"  tp: "<<tp<<endl;

}


void fill_model_1d_hist(TH1D* model_prop_1D_hist,double height,double acc_dis_width,double acc_factor){

double v_central = 2.17*pow(10,-4);
double v_width = 1*pow(10,-5);
height = height*0.01;
double vel    = pick_from_norm_distribution(v_central,v_width);
double k_centre_value = 5.10710*acc_factor;
double k    = pick_from_norm_distribution(k_centre_value,k_centre_value*acc_dis_width);
double acc = k*pow(10,-6)*vel;
double eps = calc_eps(height,vel,acc);
double tp =calc_tp(height,vel,acc);
model_prop_1D_hist->Fill(tp,1);

}




double gaussian_fit_function_for_plot(double X, double m, double s, double h, double shift){
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((pow((X-shift),-1)) - pow(m,-1)) *s;
  double y = h*(inv_sqrt_2pi *s) * std::exp(-0.5f * a * a);
  return y;
}

double gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((pow((X),-1)) - (m)) / s;
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}




double gaussian_fit_function_new(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double phi=par[3];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((((1+phi)/(X)) - (m))/s);
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}



double simple_gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((X) - (m)) / s;
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}



double gaussian_fit_function_with_shift(double* x, double* par){
  double X=x[0];
  double m=par[0];
  double s=par[1];
  double h=par[2];
  double shift=par[3];
  double inv_sqrt_2pi = 0.3989422804014327;
  double a = ((pow((X-shift),-1)) - (m)) / s;
  double y = h*(inv_sqrt_2pi) * std::exp(-0.5f * a * a);
  return y;

}


double double_gaussian_fit_function(double* x, double* par){
  double X=x[0];
  double m1=par[0];
  double s1=par[1];
  double h1=par[2];

  double m2=par[3];
  double s2=par[4];


  double h2=par[5];

  double inv_sqrt_2pi = 0.3989422804014327;
  double a1 = (X - m1) / s1;
  double a2 = (X - m2) / s2;
  double y1 = h1*(inv_sqrt_2pi / s1) * std::exp(-0.5f * a1 * a1);
  double y2 = h2*(inv_sqrt_2pi / s2) * std::exp(-0.5f * a2 * a2);

  double y_output = y1+y2;

  return y_output;

}




void combine_runs(std::string run_num, TH2D** prop_vs_z_2D_hist_array,TH1D** prop_vs_z_1D_hist_output_array,TH1D** prop_vs_z_1D_hist_counter_array){


string run_num = to_string(run_num);
string folder = name_folder(run_num);
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

prop_vs_z_2D_hist_array->Add(heatmap_hists[i]);
prop_vs_z_1D_hist_output_array->Add(counter_hists[i]);

}}



/////////////////////////////////////Main Function//////////////////////////////
int main() {


TH2D** heatmap = make_array_of_hist_2d["prop_vs_z_2D_hist_array",100,0,1,500,3000,7000];
TH1D** output =  make_array_of_hist["prop_vs_z_1D_hist_output_array",100,0,1];
TH1D** counter = make_array_of_hist["prop_vs_z_1D_hist_counter_array",100,0,1];

combine_runs(1051,heatmap,output,counter);
combine_runs(1066,heatmap,output,counter);

///////////////////////////Load File one////////////////////////////////////




  return 0;
}















