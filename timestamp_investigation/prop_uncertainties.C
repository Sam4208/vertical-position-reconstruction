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
#include "/sps/nemo/scratch/spratt/sndisplay/sndisplay.cc"
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
#include "/sps/nemo/scratch/spratt/analysis/analysis_library.h"




using namespace std;



double improved_z_calc(double R5,double R6){

   double k = 1.243*pow(10,-5);
   double term1 = R5-R6;
   double term2 = R5+R6;
   double term3 = pow(R5,2)+pow(R6,2);
   double z = (term1/(term2+((k*0.5)*(term3))))  +  0.5*(term1)/((1/k)+((term3)/(term2)));
   return z;
}



//////////////////////////////////////



double calc_tp(double d,double vel,double acc){

double L=1;

double tp = (L/vel) + 2*(acc/pow(vel,3))*(pow(d,2)+pow((L-d),2));

return tp;

}




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






void save_snd(std::string run_num,sndisplay::tracker *snd_prop, string name, int lower_limit, int upper_limit){
    string folder =  "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
    string title ="/snd_";
    string file_type =".png";
    snd_prop->setrange(lower_limit,upper_limit);
    snd_prop->draw();//saving all sndisplay plits
    snd_prop->update(false);
    snd_prop->canvas->SaveAs((folder+run_num+title+name+file_type).c_str());
}






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






int what_is_crate_number(int column){
double output=-1;

if (column<38){

   output=1;
}
if (column>37 && column<75){
   output=2;
}
if (column>74){
   output=3;
}

return output;

}




TH1D** create_3_hists_for_averaging_function(string name_of_hist, double no_bins,double lower_bound,double upper_bound){

TH1D** h = new TH1D*[3];

h[0] = new TH1D(TString::Format(name_of_hist.c_str()),TString::Format((name_of_hist).c_str()), no_bins, lower_bound, upper_bound);
h[1] = new TH1D(TString::Format((name_of_hist+"_counter").c_str()),TString::Format((name_of_hist+"_counter").c_str()), no_bins, lower_bound, upper_bound);
h[2]= new TH1D(TString::Format((name_of_hist+"_output").c_str()),TString::Format((name_of_hist+"_output").c_str()), no_bins, lower_bound, upper_bound);


return h;

}



int cell_num_calc(double side, double column, double layer){

    int cell_num =  side*9*113 + column*9 + layer;
    return cell_num;

}




void make_dir(string run_num){
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";//code to create folder that will contain sndisplay plots for each run
    string folder = run_num;
    if(mkdir((directory+folder+"/").c_str(), 0777) == -1);
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










void fill_averaging_histograms(TH1D** h,double x,double y){

h[0]->Fill(x,y);
h[1]->Fill(x,1);

}




void prop_time_asymmetry_plot(TH1D *prop_vs_h_output_hist,TH1D *output_hist, TH1D *prop_vs_h_count_hist){

    double number_bins = prop_vs_h_output_hist->GetNbinsX();
    for (int i=0;i<number_bins/2;i++){

     double no_counts_pos = prop_vs_h_count_hist->GetBinContent(i);
     double no_counts_neg = prop_vs_h_count_hist->GetBinContent(number_bins-i);
       
     double pos_h_prop_average = prop_vs_h_output_hist->GetBinContent(i);
     double neg_h_prop_average = prop_vs_h_output_hist->GetBinContent(number_bins-i);
     double diff = pos_h_prop_average-neg_h_prop_average;

     double uncertainty = sqrt(((pow(pos_h_prop_average,2))/no_counts_pos)+((pow(neg_h_prop_average,2))/no_counts_neg));
     output_hist->SetBinContent(i,diff);
     output_hist->SetBinError(i,uncertainty);

    }

}



void average_over_hist(TH1D** h){
   
   double number_bins = h[0]->GetNbinsX();

	for(int i=0;i<number_bins+1;i++){

	double total = h[0]->GetBinContent(i);
	double no = h[1]->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	//cout<<set_value<<endl;
	h[2]->SetBinContent(i,set_value);
	}}}



void average_over_hist(TH1D *hist_sum,TH1D *hist_count,TH1D *hist_out,int no_points){

	for(int i=0;i<no_points+1;i++){
	double total = hist_sum->GetBinContent(i);
	double no = hist_count->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	//cout<<set_value<<endl;
	hist_out->SetBinContent(i,set_value);
	}}}




////////////////////
//Formatting////
///////////////////

void PaintBin (TH1 *h, Int_t bin, Int_t color) 
{
   printf("%d %d %d\n", bin, color, h->GetBinContent(bin));
   TBox *b = new TBox(h->GetBinLowEdge(bin),
                      0,
                      h->GetBinWidth(bin)+h->GetBinLowEdge(bin),
                      h->GetBinContent(bin));
   b->SetFillColor(color);
   b->Draw("same");

}


////////////////////
//Checks////
///////////////////


double findAverage_vec(vector<double>* vec) {
    double sum = 0.0;
    int size = vec->size();
double output = 0;
if (size>0){

for (int i = 0; i < size; i++) {
        sum += vec->at(i);
        //cout<<vec->at(i)<<endl;
    }

 output = sum / vec->size();
}

    return output;
}

bool check_ts_valid(double R){
               bool output = true;
               if (R<0 || R>100000){
                  output=false;
               }
               return output;
               }

bool check_TS_paired(double R1,double R3){
                       bool output = false;
                        if (check_ts_valid(R1)==true || check_ts_valid(R3)==true){
                           if (abs(R1-R3)<200){
                              output=true;
                           }

                        }
                        return output;
                        }

int calc_number(bool r1, bool r2,bool r3, bool r4,bool r5, bool r6){
int output = 1;

  if (r1 == false){
       output = 2*output;
  }
   if (r2 == false){
       output = 3*output;
  }
   if (r3 == false){
       output = 4*output;
  }
   if (r4 == false){
       output = 5*output;
  }
   if (r5 == false){
       output = 6*output;
  }
   if (r6 == false){
       output = 7*output;
  }


  return output;

}

int pair_timestamps(double R5, double R6, double R_check){

               int pair;


               if (R_check<R5+300 && R_check>R5-300){
               pair = 1;

               }


               if (R_check<R6+300 && R_check>R6-300){
                  pair = 2;

               }


               if (R_check<R5+300 && R_check>R5-300 && R_check<R6+300 && R_check>R6-300){
               pair = 3;
               }


               return pair;

               }




int check_scenario(bool r1, bool r2,bool r3, bool r4,bool r5, bool r6){

  int output = 11;

  if (r1=true && r2==true && r3==true && r4==true && r5==true && r6==true){
   output=0;
  }

  if (r1=true && r2==false && r3==true &&  r4==true &&  r5==true && r6==true){//R2
   output=1;
  }
 
  if (r1=true && r2==true && r3==true &&  r4==true && r5==false && r6==true){//R3
   output=2;
  }

  if (r1=true && r2==true && r3==true &&  r4==true &&  r5==true && r6==false){//R6
   output=4;
  }

  if (r1=true && r2==false && r3==true &&  r4==false &&  r5==true && r6==true){//R2 R4
   output=3;
  }

  if (r1=true && r2==false && r3==true &&  r4==true && r5==false && r6==true){//R2 R5
   output=5;
  }

  if (r1=true && r2==false && r3==true && r4==true && r5==true && r6==false){//R2 R6
   output=6;
  }

  if (r1=true && r2==true && r3==true && r4==true && r5==false && r6==false){// R5 R6
   output=7;
  }

  if (r1=true && r2==false && r3==true && r4==false && r5==false && r6==true){//R2 R4 R5
   output=8;
  }

  if (r1=true && r2==false && r3==true && r4==false &&  r5==true && r6==false){//R2 R4 R6
   output=9;
  }

  if (r1=true && r2==false && r3==true && r4==true &&  r5==false && r6==false){//R2 R5 R6
   output=10;
  }

  return output;
}


int locate_deadtime(double dead_time, double max_dead_time, int num_points){


double vec_number=-10;

double interval = max_dead_time/num_points;

//cout<<shifted_z<<endl;
for (int i=0;i<num_points;i++){
       if (i*interval<dead_time && dead_time<(i+1)*interval){
       //  cout<<i*interval<<"    "<<dead_time<<"   "<<(i+1)*interval<<"    "<<i<<endl;
              vec_number = i;}
}

return vec_number;
}



void prop_time_vs_dead_time_function_modified(TH1D *dead_time_h1_min_h2_hist , TH1D *dead_time_h1_min_h2_counter_hist, TH2D *dead_time_h1_vs_h2_hist, TH2D *dead_time_h1_vs_h2_counter_hist, std::vector<double>** vector_prop_time_vs_dead_time_half_triggered , std::vector<double>** vector_prop_time_vs_dead_time_modified, std::vector<double>** vector_prop_time_vs_dead_time, double R0_previous, double R5_previous, double R6_previous,double R0_sec, double R5, double R6, int num_points, double max_dead_time,double &counter_all, double &counter_select){

counter_all=counter_all+1;
double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;
double dead_time = R0_sec - R0_previous;

double z_previous = (R5_previous-R6_previous)/(R5_previous+R6_previous);
double prop_time_previous = R5_previous+R6_previous;

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true ){
   if (check_ts_valid(R5_previous)==true && check_ts_valid(R6_previous)==false){
   

     int vec_number = locate_deadtime(dead_time, max_dead_time,num_points);

      if (vec_number != -10 && vec_number <num_points  && vec_number >-1 ){
      vector_prop_time_vs_dead_time_half_triggered[vec_number]->push_back(prop_time);
      
}

   }

   if (check_ts_valid(R5_previous)==false && check_ts_valid(R6_previous)==true){
   
     int vec_number = locate_deadtime(dead_time, max_dead_time,num_points);

      if (vec_number != -10 && vec_number <num_points  && vec_number >-1 ){
      vector_prop_time_vs_dead_time_half_triggered[vec_number]->push_back(prop_time);

   }

}

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true && check_ts_valid(R5_previous)==true && check_ts_valid(R6_previous)==true){

if (prop_time<0.5*pow(10,6)){
dead_time_h1_vs_h2_hist->Fill(z,z_previous,prop_time);
dead_time_h1_vs_h2_counter_hist->Fill(z,z_previous);
dead_time_h1_min_h2_hist->Fill(z-z_previous,prop_time);
dead_time_h1_min_h2_counter_hist->Fill(z-z_previous,1);
}

counter_select=counter_select+1;

int vec_number = locate_deadtime(dead_time, max_dead_time,num_points);

if (vec_number != -10 && vec_number <num_points  && vec_number >-1 ){
vector_prop_time_vs_dead_time_modified[vec_number]->push_back(prop_time-prop_time_previous);
vector_prop_time_vs_dead_time[vec_number]->push_back(prop_time);

}

}

}}









void prop_time_vs_dead_time_function(std::vector<double>** vector_prop_time_vs_dead_time,  double previous_ts_hit,double R0_sec, double R5, double R6, int num_points, double max_dead_time){

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;
double dead_time = R0_sec - previous_ts_hit;
//double max_dead_time = 0.01;
//cout<<dead_time<<"       "<<prop_time<<endl;
int vec_number = locate_deadtime(dead_time, max_dead_time,num_points);
//cout<<vec_number<<endl;
if (vec_number != -10 && vec_number <num_points  && vec_number >-1 ){
vector_prop_time_vs_dead_time[vec_number]->push_back(prop_time);
//cout<<prop_time<<"   prop_time      dead_time"<<dead_time<<endl;
}}}



void R5_vs_R6_function(std::vector<double>** vector, double R5, double R6, int num_points, double max_R5){

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

    if (R5>100 && R5<5500 && R6>100 && R6<5500){

int vec_number = locate_deadtime(R6, max_R5,num_points);
if (vec_number != -10 && vec_number <num_points  && vec_number >-1 ){
vector[vec_number]->push_back(R5);
}}}}


void FindMinMax(TGraph *graph, double &minX, double &maxX, double &minY, double &maxY) {
    int nPoints = graph->GetN();
    double *xValues = graph->GetX();
    double *yValues = graph->GetY();

    if (nPoints > 0) {
        minX = maxX = xValues[0];
        minY = maxY = yValues[0];

        for (int i = 1; i < nPoints; ++i) {
            if (xValues[i] < minX)
                minX = xValues[i];
            if (xValues[i] > maxX)
                maxX = xValues[i];
            if (yValues[i] < minY)
                minY = yValues[i];
            if (yValues[i] > maxY)
                maxY = yValues[i];
        }}}


void graph_function_2(TH1D *hist, std::vector<double>** vector_prop_time_vs_z,int num_points, string title, string x_axis, string y_axis, double max_value){

hist->SetTitle(title.c_str());
hist->SetName(title.c_str());
TAxis *yAxis = hist->GetYaxis();
TAxis *xAxis = hist->GetXaxis();
yAxis->SetLabelSize(0.02);       // Set the font size of the X-axis labels (in pixels)
yAxis->SetTitleSize(0.04); 
xAxis->SetTitle(x_axis.c_str());
yAxis->SetTitle(y_axis.c_str());

for (int vec_num=0;vec_num<num_points;vec_num++){

size_t no_events_at_this_z = vector_prop_time_vs_z[vec_num]->size();
//cout<< vec_num<<"        " <<findAverage_vec(vector_prop_time_vs_z[vec_num])<<endl;
if (std::isnan(no_events_at_this_z) == false){

   double num_points_double = num_points;
   double vec_num_double = vec_num;
   double x_axis_value = (0.5*max_value/num_points)+(vec_num_double/num_points_double)*max_value;
  // cout<<x_axis_value<<"       h"<<endl;
cout<<vec_num_double<<"      "<<no_events_at_this_z<<endl;
hist->SetBinContent(vec_num_double+1,no_events_at_this_z); 

}}

hist->SetMarkerStyle(20);  // Circle marker
hist->SetMarkerSize(0.5);
}







void graph_function(TGraph *graph, std::vector<double>** vector_prop_time_vs_z,int num_points, string title, string x_axis, string y_axis, double max_value, double factor){

max_value = max_value*factor;
graph->SetTitle(title.c_str());
graph->SetName(title.c_str());
TAxis *yAxis = graph->GetYaxis();
TAxis *xAxis = graph->GetXaxis();
yAxis->SetLabelSize(0.02);       // Set the font size of the X-axis labels (in pixels)
yAxis->SetTitleSize(0.04); 
xAxis->SetTitle(x_axis.c_str());
yAxis->SetTitle(y_axis.c_str());


for (int vec_num=0;vec_num<num_points;vec_num++){

double avg_plasma_prop_time = findAverage_vec(vector_prop_time_vs_z[vec_num]);
//cout<< vec_num<<"        " <<findAverage_vec(vector_prop_time_vs_z[vec_num])<<endl;
if (std::isnan(avg_plasma_prop_time) == false){

   double num_points_double = num_points;
   double vec_num_double = vec_num;
   double x_axis_value = (0.5*max_value/num_points)+(vec_num_double/num_points_double)*max_value;
  // cout<<x_axis_value<<"       h"<<endl;

graph->SetPoint(vec_num, x_axis_value, avg_plasma_prop_time*factor); 
}}

double minX, maxX, minY, maxY;
FindMinMax(graph, minX, maxX, minY, maxY);
//yAxis->SetRangeUser(minY-abs(0.3*minY),maxY+abs(0.3*minY));
//xAxis->SetRangeUser(0,max_value);
yAxis->SetRangeUser(4000,6000);
xAxis->SetRangeUser(0,max_value);
graph->SetMarkerStyle(20);  // Circle marker
graph->SetMarkerSize(0.5);
//yAxis->SetRange(0, 7000);
}






void fit_plasma_prop_time_vs_eps(TGraph *graph1,TF1 *fitFunction){

    fitFunction->SetParLimits(0, 0.3,1.9);
    fitFunction->SetParLimits(1, 0,10);

    //fitFunction->SetParLimits(2, 0.01,10);

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit
    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
  //  double par_2 = fitFunction->GetParameter(2);
  
    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;
  //  std::cout << "acc: " << par_2 << std::endl;

}


void fit_quadratic(TGraph *graph1,TF1 *fitFunction){


    // [0] is the slope, [1] is the intercept

    // Fit the graph with the function
    fitFunction->SetParLimits(0, 0,10000);
    fitFunction->SetParLimits(1, 0,10000);
    fitFunction->SetParLimits(2, 0,10000);
    graph1->Fit(fitFunction, "R"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    double par_2 = fitFunction->GetParameter(2);

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "a: " << par_0 << std::endl;
    std::cout << "b: " << par_1 << std::endl;
    std::cout << "c: " << par_2 << std::endl;

}


void fit_iteration_eps_vs_z(TGraph *graph1,TF1 *fitFunction){

    fitFunction->SetParLimits(0,0, 0.1);
    fitFunction->SetParLimits(1,0, 0.01);

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit

    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
    cout<<"Iteration Fit"<<endl;

    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;

}



void fit_distance_function_to_t_graph(TGraph *graph1,TF1 *fitFunction){


    // [0] is the slope, [1] is the intercept

    // Fit the graph with the function
    fitFunction->SetParLimits(0, 0,10);
    fitFunction->SetParLimits(1, 0,1);
    //fitFunction->SetParLimits(2, 0.01,10);

    graph1->Fit(fitFunction, "R"); // "R" for range-based fit
    double par_0 = fitFunction->GetParameter(0);
    double par_1 = fitFunction->GetParameter(1);
  //  double par_2 = fitFunction->GetParameter(2);
  
    std::cout << "Fit Parameters:" << std::endl;
    std::cout << "vel: " << par_0 << std::endl;
    std::cout << "acc: " << par_1 << std::endl;
  //  std::cout << "acc: " << par_2 << std::endl;

}


void graph_function_flip_x_y_axis(TGraph *graph1, TGraph *graph2){


    for (int i = 0; i < graph1->GetN(); ++i) {
        double x, y;
        graph1->GetPoint(i, x, y);
        // Apply your manipulation (e.g., multiply y-values by 2)
       // y= 3-(3*(y*(1/pow(4,-5))));
        graph2->SetPoint(i, y/1000, x);
        graph2->SetMarkerStyle(20);  // Circle marker
        graph2->SetMarkerSize(0.5);
}}


void graph_function_split_graph_and_flip_x_y(TGraph *graph1, TGraph *graph2, TGraph *graph3){

    int pointIndex;



    for (int i = graph1->GetN()/10; i < graph1->GetN()/2; ++i) {
        double x, y;
        graph1->GetPoint(i, x, y);

        pointIndex = graph2->GetN(); // Get the current number of point
        graph2->SetPoint(pointIndex, y/1000, x-0.5);
        }

        graph2->SetMarkerStyle(20);  // Circle marker
        graph2->SetMarkerSize(0.5);

    for (int i = graph1->GetN()/2 ; i < graph1->GetN()-graph1->GetN()/10; ++i) {
        double x, y;
        graph1->GetPoint(i, x, y);
        
        pointIndex = graph3->GetN(); 
        graph3->SetPoint(pointIndex, y/1000, x-0.5);
        }

        graph3->SetMarkerStyle(20);  // Circle marker
        graph3->SetMarkerSize(0.5);
        }


void graph_function_plasma_speed_plot(TGraph *graph1, TGraph *graph2){


    for (int i = 0; i < graph1->GetN(); ++i) {
        double x, y;
        graph1->GetPoint(i, x, y);

        graph2->SetPoint(i, x*100000, y*100000
        );
        graph2->SetMarkerStyle(20);  // Circle marker
        graph2
        ->SetMarkerSize(0.5);
      }}





int locate_z(double z, double max_z, int num_points){
double vec_number=-10;
double interval = max_z/num_points;
double shifted_z = (z+1)*0.5;
//cout<<shifted_z<<endl;
for (int i=0;i<num_points;i++){
       if (i*interval<shifted_z && shifted_z<(i+1)*interval){
       //  cout<<i*interval<<"    "<<dead_time<<"   "<<(i+1)*interval<<"    "<<i<<endl;
              vec_number = i;}
}
return vec_number;
}



void R5_vs_R6_hist_function(TH2D *hist, double R5, double R6){

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
hist->Fill(R5,R6);
}}




void prop_vs_h_function(TH1D *prop_vs_h_hist,TH1D *prop_vs_h_counter_hist,double R5,double R6){

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;
prop_vs_h_hist->Fill(z,prop_time);
prop_vs_h_counter_hist->Fill(z,1);
}}


void prop_time_vs_z_function(std::vector<double>** vector_prop_time_vs_dead_time, double previous_ts_hit,double R0_sec, double R5, double R6, int num_points){

//if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;
double dead_time = R0_sec - previous_ts_hit;
double max_z = 1;

int vec_number = locate_z(z, max_z ,num_points);
if (prop_time<pow(10,8) && prop_time>1){
if (vec_number != -10 && vec_number <num_points  && vec_number >-1){
   
vector_prop_time_vs_dead_time[vec_number]->push_back(prop_time);
}}}



void prop_time_vs_z_function_modified(std::vector<double>** vector_prop_time_vs_R5,std::vector<double>** vector_prop_time_vs_R6, double previous_ts_hit,double R0_sec, double R5, double R6, int num_points){
if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;
double dead_time = R0_sec - previous_ts_hit;
double max_z = 1;
int vec_number = locate_z(z, max_z ,num_points);
//cout<<vec_number<<endl;
//if (vec_number==num_points){cout<<vec_number;}
if (vec_number != -10 && vec_number <num_points  && vec_number >-1){
vector_prop_time_vs_R5[vec_number]->push_back(R5);
vector_prop_time_vs_R6[vec_number]->push_back(R6);
}}}



//note used anymore
void mean_prop_time_z_calc_vectors(double R5, double R6,std::vector<double>* vector_0, std::vector<double>* vector_1,std::vector<double>* vector_2,std::vector<double>* vector_3,std::vector<double>* vector_4,std::vector<double>* vector_5,std::vector<double>* vector_6,std::vector<double>* vector_7,std::vector<double>* vector_8,std::vector<double>* vector_9){
if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;

if (z>-1 && z<-0.8){
   vector_0->push_back(prop_time);
}
if (z>-0.8 && z<-0.6){
   vector_1->push_back(prop_time);
}
if (z>-0.6 && z<-0.4){
   vector_2->push_back(prop_time);
}
if (z>-0.4 && z<0.2){
   vector_3->push_back(prop_time);
}
if (z>-0.2 && z<0){
   vector_4->push_back(prop_time);
}
if (z>0 && z<0.2){
   vector_5->push_back(prop_time);
}
if (z>0.2 && z<0.4){
   vector_6->push_back(prop_time);
}
if (z>0.4 && z<0.6){
   vector_7->push_back(prop_time);
}
if (z>0.6 && z<0.8){
   vector_8->push_back(prop_time);
}
if (z>0.8 && z<1){
   vector_9->push_back(prop_time);
}}}
//not used anymore
void z_vs_prop_time_plot_function(double R5, double R6, TH2D* hist){

if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

double z = (R5-R6)/(R5+R6);
double prop_time = R5+R6;

hist->Fill(z,prop_time);

}}


void categorise_peaks(TH1D *hist,TH1D *hist_lower, TH1D *hist_mid, TH1D *hist_higher, int title, int bins, int lower_bounds, int upper_bounds){

    double centre_point = (hist->GetMaximumBin()*(upper_bounds-lower_bounds)/bins)+lower_bounds;
    if (centre_point>4500){
    hist_higher -> Add(hist);
    }
    if (centre_point<4500 && centre_point>4000){
    hist_mid -> Add(hist);
    }
    if (centre_point<4000){
    hist_lower -> Add(hist);
    }
    //snd_prop->update_canvas();
    //tree_new->Fill();

}

void prop_time_vs_time_hit_function(double R0, double R5, double R6, TH2D* hist,int cell_num){
if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
if (cell_num==248){
   hist->Fill(R5+R6,R0/pow(10,9),1);
   //cout<<R0/pow(10,9)<<endl;
}}

}

bool print_out_timestamp_scenario_sub_function(double R, string name){
    bool output = true;
if (check_ts_valid(R)==false){
//cout<<name.c_str()<<"      ";
output=false;
}
return output;
}

void print_out_timestamp_scenario(double R1, double R2, double R3, double R4, double R5, double R6, TH1D* scenario_counter_hist, TH1D* scenario_basic_hist, TH1D* R3_minus_R1_hist_case_1){
   
                  bool r1 = print_out_timestamp_scenario_sub_function(R1,"R1");
                  bool r2 = print_out_timestamp_scenario_sub_function(R2,"R2");
                  bool r3 = print_out_timestamp_scenario_sub_function(R3,"R3");
                  bool r4 = print_out_timestamp_scenario_sub_function(R4,"R4");
                  bool r5 = print_out_timestamp_scenario_sub_function(R5,"R5");
                  bool r6 = print_out_timestamp_scenario_sub_function(R6,"R6");

                  int number = calc_number(r1,r2,r3,r4,r5,r6);
                  // cout<<number<<endl;
                  scenario_counter_hist->Fill(number,1);
                  int scenario = check_scenario(r1,r2,r3,r4,r5,r6);
                  if (scenario == 7){

                     R3_minus_R1_hist_case_1->Fill(R3-R1,1);
                  }
                  scenario_basic_hist->Fill(scenario,1);

                  check_TS_paired(R1,R3);
                  check_TS_paired(R2,R4);
               }


THStack* stack_histograms_6(TH1D *hist1, TH1D *hist2, TH1D *hist3,  TH1D *hist4){



               THStack* hs = new THStack("hs","Stacked 1D histograms");
                  
                  hist1->SetFillColor(kRed);
                  hist1->SetMarkerStyle(21);
                  hist1->SetMarkerColor(kRed);
                  hs->Add(hist1);
                  
                  hist2->SetFillColor(kBlue);
                  hist2->SetMarkerStyle(21);
                  hist2->SetMarkerColor(kBlue);
                  hs->Add(hist2);

                  hist3->SetFillColor(kGreen);
                  hist3->SetMarkerStyle(21);
                  hist3->SetMarkerColor(kGreen);
                  hs->Add(hist3);

                  hist4->SetFillColor(kYellow);
                  hist4->SetMarkerStyle(21);
                  hist4->SetMarkerColor(kYellow);
                  hs->Add(hist4);

                  return hs;
  
}

void stack_histograms_4(TCanvas* c1, TH1D* R1, TH1D* R2, TH1D* R3, TH1D* R4){
   
   R1->SetTitle("R-R1");
   
   R2->SetTitle("R-R2");
   
   R3->SetTitle("R-R3");
  
   R4->SetTitle("R-R4");
   
   gStyle->SetOptStat(0);

   c1->Divide(2,2);
   c1->cd(1); R1->Draw();
   c1->cd(2); R2->Draw(); 
   c1->cd(3); R3->Draw(); 
   c1->cd(4); R4->Draw(); 
  
}

void stack_histograms_3(TCanvas* c10, TH1D* R1_hist, TH1D* R1_hist_b, TH1D* R1_hist_c, TH1D* R1_hist_d, TH1D* R2_hist, TH1D* R2_hist_b, TH1D* R2_hist_c ,TH1D* R2_hist_d){

   THStack* R1_st = stack_histograms_6(R1_hist_c,R1_hist_b,R1_hist,R1_hist_d);
   R1_st->SetTitle("R3-R1");
   THStack* R2_st = stack_histograms_6(R2_hist_c,R2_hist_b,R2_hist,R2_hist_d);
   R2_st->SetTitle("R4-R2");

   gStyle->SetOptStat(0);

   c10->Divide(2,1);
   c10->cd(1); R1_st->Draw();
   c10->cd(2); R2_st->Draw(); 

   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(R1_hist_c,"R5 R6 val");
   legend->AddEntry(R1_hist_b,"R6 inv");
   legend->AddEntry(R1_hist,"R5 inv");
   legend->AddEntry(R1_hist_d,"R5 R6 inv");
   legend->Draw();


}

void stack_histograms_2(TCanvas* c12, TH1D* R1_hist, TH1D* R1_hist_b, TH1D* R1_hist_c, TH1D* R1_hist_d, TH1D* R2_hist, TH1D* R2_hist_b, TH1D* R2_hist_c ,TH1D* R2_hist_d , TH1D* R3_hist, TH1D* R3_hist_b, TH1D* R3_hist_c,TH1D* R3_hist_d , TH1D* R4_hist, TH1D* R4_hist_b, TH1D* R4_hist_c, TH1D* R4_hist_d){



   THStack* R1_st = stack_histograms_6(R1_hist_c,R1_hist_b,R1_hist,R1_hist_d);
   R1_st->SetTitle("R1");
   THStack* R2_st = stack_histograms_6(R2_hist_c,R2_hist_b,R2_hist,R2_hist_d);
   R2_st->SetTitle("R2");
   THStack* R3_st = stack_histograms_6(R3_hist_c,R3_hist_b,R3_hist,R3_hist_d);
   R3_st->SetTitle("R3");
   THStack* R4_st = stack_histograms_6(R4_hist_c,R4_hist_b,R4_hist,R3_hist_d);
   R4_st->SetTitle("R4");


   
   gStyle->SetOptStat(0);

   c12->Divide(2,2);
   c12->cd(1); R1_st->Draw();
   c12->cd(2); R2_st->Draw(); 
   c12->cd(3); R3_st->Draw(); 
   c12->cd(4); R4_st->Draw(); 

   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(R1_hist_c,"R5 R6 val");
   legend->AddEntry(R1_hist_b,"R6 inv");
   legend->AddEntry(R1_hist,"R5 inv");
   legend->AddEntry(R1_hist_d,"R5 R6 inv");
   legend->Draw();



}

Double_t prop_vs_eps_fit_upper(Double_t *x,Double_t *par) {

   Double_t prop_time = x[0];   
   Double_t vel = par[0];
   Double_t acc = par[1];

   Double_t L=3;



   Double_t factor = 4*pow(L,2) - 8*((((L/vel) - prop_time)*(pow(vel,3)/(2*acc)))+ pow(L,2));

   Double_t d_true = (L/2) + (sqrt(factor)/4);


   Double_t eps = ((((2*d_true-L)*pow(vel,2))+(((2*acc)*(pow(d_true,2) - pow((L-d_true),2)))))*(L/2))/(prop_time*pow(vel,3));

   return eps;
   }


Double_t prop_vs_eps_fit_lower(Double_t *x,Double_t *par) {

   Double_t prop_time = x[0];   
   Double_t vel = par[0];
   Double_t acc = par[1];

   Double_t L=3;



   Double_t factor = 4*pow(L,2) - 8*((((L/vel) - prop_time)*(pow(vel,3)/(2*acc)))+ pow(L,2));

   Double_t d_true = (L/2) - (sqrt(factor)/4
   );


   Double_t eps = (((2*d_true-L)*pow(vel,2))+(((2*acc)*(pow(d_true,2) - pow((L-d_true),2)))*(L/2)))/(prop_time*pow(vel,3));

   return eps;
   }


Double_t quadratic_fitting_function(Double_t *x,Double_t *par) {
    Double_t L = 3;
    Double_t X = x[0]-0.5;

    Double_t vel = par[0];
    Double_t acc = par[1];
 

    Double_t fitval_2 = ((L/vel)+(2*(acc)/pow(vel,3))*pow(L,2)) +  (((-4*acc*L)/pow(vel,3))*X)   +   ((4*acc)/pow(vel,3))*X*X;

    Double_t fitval = par[0]+(par[1]*(X-0.5))+(par[2]*X*X);

    return fitval;
   }



Double_t skewed_gaussian(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];

    Double_t fitval = 2*par[0]*TMath::Exp(-0.5*arg*arg)*ROOT::Math::normal_cdf(par[3]*x[0],1,0);
    return fitval;
   }



Double_t double_gaussian(Double_t *x,Double_t *par) {
    Double_t arg1 = 0;
    Double_t arg2 = 0;

    if (par[2]!=0) 
    arg1 = (x[0] - par[1])/par[2];
    arg2 = (x[0] - par[4])/par[2];

    Double_t fitval = par[0]*TMath::Exp(-0.5*arg1*arg1)+par[3]*TMath::Exp(-0.5*arg2*arg2);

    return fitval;
   }


   Double_t linearly_decelerating_plasma(Double_t *x,Double_t *par) {



    Double_t fitval =   (3/par[0]) +   x[0]   - (par[1]/par[0])*x[0]*x[0]    -     (1/(8*pow(par[0],3)))*pow((    -3  +    (x[0]*par[0])    -  par[1]*x[0]*x[0]  ) ,2 )                 ;

   
    return fitval;
   }


void perform_fit_1(TH1D *hist, int title){

   TF1 *func = new TF1(TString::Format("fit_%d", title),skewed_gaussian,3000,6000,4);
      
      func->SetParameters(hist->GetMaximum(),(hist->GetMaximumBin()*(6000-3000)/300)+3000,hist->GetRMS(),1);
      func->SetParLimits(3,0.01,10);
      func->SetParLimits(2,0,4000);
      func->SetParNames ("Constant","Mean_value","Sigma","Skew");
       hist->Fit(TString::Format("fit_%d", title));
       //cout<<title<<endl;   
     
}

void perform_fit_2(TH1D *hist, int title, int bins, int lower_bounds, int upper_bounds){

    TF1 *func = new TF1(TString::Format("fit_%d", title),double_gaussian,3000,6000,5);

    double centre_point = (hist->GetMaximumBin()*(upper_bounds-lower_bounds)/bins)+lower_bounds;

    func->SetParameters(hist->GetMaximum(),centre_point,hist->GetRMS(),hist->GetMaximum()/10,centre_point+10);

    func->SetParLimits(0,hist->GetMaximum()-20,hist->GetMaximum()+20);//height1 bounds
    func->SetParLimits(1,centre_point-100,centre_point+100);//centre1 bounds
    func->SetParLimits(2,0,4000);//sd bounds
    func->SetParLimits(3,0,hist->GetMaximum()+20);//height2 bounds
    func->SetParLimits(4,centre_point-500,centre_point+500);//centre2 bounds
    
    func->SetParNames ("Height","Centre Point","sigma","Height2","Centre Point2");
    hist->Fit(TString::Format("fit_%d", title),"q");
    //cout<<title<<endl;  
    
}



















int main() {

string run_num = "1051";
//std::string file_1 ="/sps/nemo/scratch/chauveau/commissioning/sam/snemo_run-1544_track-udd.root";
//std::string file_1 = "/sps/nemo/scratch/golivier/aussois_tuto_data/bi207_run/snemo_run-840_udd.root";//Chosen run file 
 
//std::string file_1 = "/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/snemo_run-1066_udd.root";

std::string file_1 ="/sps/nemo/scratch/golivier/data/UDD/snemo_run-1051_udd.root";
//std::string file_1 ="/sps/nemo/scratch/golivier/data/UDD/snemo_run-974_udd.root";
//std::string file_1 ="/sps/nemo/scratch/spratt/analysis/UDD_Data/snemo_run-974_udd.root";

//std::string file_1 ="/sps/nemo/scratch/spratt/analysis/UDD_Data/snemo_run-1018_udd.root";
//std::string file_1 ="/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/snemo_run-1116_udd.root";

make_dir(run_num);


TFile *file= new TFile(file_1.c_str(), "READ");


  std::vector<std::vector<short>> *wave = new std::vector<std::vector<short>>;//initiating vectors that will be filled with data from the event file
  std::vector<int> *calo_wall   = new std::vector<int>;
  std::vector<int> *calo_side   = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row    = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;
  std::vector<int> *calo_ampl   = new std::vector<int>;
  std::vector<int> *calo_ts     = new std::vector<int>;
  std::vector<int> *tracker_side   = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_row    = new std::vector<int>;
  std::vector<int> *tracker_layer  = new std::vector<int>;
  std::vector<std::vector<int64_t> >*tracker_botcathtime  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_topcathtime  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_R0  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_R1  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_R2  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_R3  = new std::vector<std::vector<int64_t> >;
  std::vector<std::vector<int64_t> > *tracker_R4  = new std::vector<std::vector<int64_t> >;

  int ncalo_main_wall = 520;

  int eventnumber = 0;
  int calo_nohits = 0;
  int tracker_nohits = 0;

  TTree* tree = (TTree*)file->Get("SimData");//getting event data from the input file
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("header.eventnumber",1);
  tree->SetBranchAddress("header.eventnumber", &eventnumber);
  
  //tree->SetBranchStatus("digicalo.nohits",1);
  //tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
  //tree->SetBranchStatus("digicalo.waveform",1);
  //tree->SetBranchAddress("digicalo.waveform", &wave);
  //tree->SetBranchStatus("digicalo.wall",1);
  //tree->SetBranchAddress("digicalo.wall", &calo_wall);
  //tree->SetBranchStatus("digicalo.side",1);
  //tree->SetBranchAddress("digicalo.side", &calo_side);
  //tree->SetBranchStatus("digicalo.column",1);
  //tree->SetBranchAddress("digicalo.column", &calo_column);
  //tree->SetBranchStatus("digicalo.row",1);
  //tree->SetBranchAddress("digicalo.row", &calo_row);
  //tree->SetBranchStatus("digicalo.charge",1);
  //tree->SetBranchAddress("digicalo.charge", &calo_charge);
  //tree->SetBranchStatus("digicalo.peakamplitude",1);
  //tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
  //tree->SetBranchStatus("digicalo.timestamp",1);
  //tree->SetBranchAddress("digicalo.timestamp", &calo_ts);


  tree->SetBranchStatus("digitracker.nohits",1);
  tree->SetBranchAddress("digitracker.nohits", &tracker_nohits);
  tree->SetBranchStatus("digitracker.side",1);
  tree->SetBranchAddress("digitracker.side", &tracker_side);
  tree->SetBranchStatus("digitracker.layer",1);
  tree->SetBranchAddress("digitracker.layer", &tracker_layer);
  tree->SetBranchStatus("digitracker.column",1);
  tree->SetBranchAddress("digitracker.column", &tracker_column);
  tree->SetBranchStatus("digitracker.anodetimestampR0",1);
  tree->SetBranchAddress("digitracker.anodetimestampR0", &tracker_R0);
  tree->SetBranchStatus("digitracker.anodetimestampR1",1);
  tree->SetBranchAddress("digitracker.anodetimestampR1", &tracker_R1);
  tree->SetBranchStatus("digitracker.anodetimestampR2",1);
  tree->SetBranchAddress("digitracker.anodetimestampR2", &tracker_R2);
  tree->SetBranchStatus("digitracker.anodetimestampR3",1);
  tree->SetBranchAddress("digitracker.anodetimestampR3", &tracker_R3);
  tree->SetBranchStatus("digitracker.anodetimestampR4",1);
  tree->SetBranchAddress("digitracker.anodetimestampR4", &tracker_R4);
  tree->SetBranchStatus("digitracker.bottomcathodetimestamp",1);
  tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &tracker_botcathtime);
  tree->SetBranchStatus("digitracker.topcathodetimestamp",1);
  tree->SetBranchAddress("digitracker.topcathodetimestamp", &tracker_topcathtime);

  //gStyle->SetOptStat(0);








     //////////////////////////
     ////    SNDISPLAY    ////
     //////////////////////////

  sndisplay::tracker *snd_prop= new sndisplay::tracker ("tracker_prop_times",true);//initiating sndisplay tracker canvas
  sndisplay::tracker *snd_check= new sndisplay::tracker ("check",true);//initiating sndisplay tracker canvas
  sndisplay::tracker *bad_cells= new sndisplay::tracker ("check_bad",true);//initiating sndisplay tracker canvas
  sndisplay::tracker *array_check= new sndisplay::tracker ("array_check",true);//initiating sndisplay tracker canvas
  sndisplay::tracker *activity_snd= new sndisplay::tracker ("activity_snd",true);//initiating sndisplay tracker canvas

//////BAD CELL CHECK////

for (int side=0; side<2;side++){
for (int layer=0;layer<9;layer++){
for (int column=0;column<113;column++){

activity_snd->fill(cell_num_calc(side,column,layer),1);
if (is_cell_bad(side, column, layer) == true){
bad_cells->setcontent(cell_num_calc(side,column,layer),1);
}


}}}





  TCanvas *c100 = new TCanvas("c100","hists with different scales",600,400);
  
  double get_entries = (tree->GetEntries() - 1);
  int no_cells = 2035;

  double mean = 4214.52;
  double counter = 0;
  double prop_time_sum = 0;
  double prop_time=0;
  double sd_sum = 0;
  int cell_num;
  int bins=500;
  int lower_bounds = 3000;
  int upper_bounds = 6000;
  double centre_point;




  std::vector<double> *v0   = new std::vector<double>;
  std::vector<double> *v1   = new std::vector<double>;
  std::vector<double> *v2   = new std::vector<double>;
  std::vector<double> *v3   = new std::vector<double>;
  std::vector<double> *v4   = new std::vector<double>;
  std::vector<double> *v5   = new std::vector<double>;
  std::vector<double> *v6   = new std::vector<double>;
  std::vector<double> *v7   = new std::vector<double>;
  std::vector<double> *v8   = new std::vector<double>;
  std::vector<double> *v9   = new std::vector<double>;

TMultiGraph* multiGraph = new TMultiGraph();

int num_points_prop_vs_z_cell_2000 = 10;//30
TGraph* mean_prop_time_vs_z_cell_2000_graph = new TGraph(num_points_prop_vs_z_cell_2000);

int num_points_R5 = 100;//30
TGraph* R5_vs_R6_graph = new TGraph(num_points_R5);
TGraph* plasma_speed_graph = new TGraph(num_points_R5);




int num_points_prop_vs_z_cell_1000 = 100;//30
TGraph* mean_prop_time_vs_z_cell_1000_graph = new TGraph(num_points_prop_vs_z_cell_1000);

int num_points_prop_vs_z_cell_100 = 100;//30

TGraph* mean_prop_time_vs_z_cell_100_graph = new TGraph(num_points_prop_vs_z_cell_100);

int num_points_prop_vs_z = 100;//300
TGraph* mean_prop_time_vs_z_graph = new TGraph(num_points_prop_vs_z);
TGraph* mean_prop_time_vs_z_inner_graph = new TGraph(num_points_prop_vs_z);
TGraph* mean_prop_time_vs_z_outer_graph = new TGraph(num_points_prop_vs_z);
TGraph* mean_prop_time_vs_z_near_bad_cell_graph = new TGraph(num_points_prop_vs_z);
TGraph* mean_prop_time_vs_z_recent_hit_graph = new TGraph(num_points_prop_vs_z);

TGraph* mean_prop_time_vs_z_graph_lower_z = new TGraph();
TGraph* mean_prop_time_vs_z_graph_upper_z = new TGraph();

TGraph* mean_prop_time_vs_R5_graph = new TGraph(num_points_prop_vs_z);
TGraph* mean_prop_time_vs_R6_graph = new TGraph(num_points_prop_vs_z);

int num_points = 100;//120
TGraph* prop_time_vs_dead_time_graph = new TGraph(num_points);

int num_points_dead_time_modified = 100;//120
TGraph* prop_time_vs_dead_time_modified_graph = new TGraph(num_points_dead_time_modified);

TGraph* prop_time_vs_dead_time_half_triggered_graph = new TGraph(num_points);

int num_points_dead_time_neighbour = 100;//300
TGraph* prop_time_vs_dead_time_neighbour_graph = new TGraph(num_points_dead_time_neighbour);

TH1D* hist_z_no = new TH1D("hist_z_no","hist_z_no",num_points_prop_vs_z, 0, 1);



     //////////////////////////
     ////    HISTOGRAMS    ////
     //////////////////////////

  

  TH1D* R1_hist = new TH1D("R1","R1", 100, -1000, 7000);
  TH1D* R2_hist = new TH1D("R2","R2", 100, -1000, 7000);
  TH1D* R3_hist = new TH1D("R3","R3", 100, -1000, 7000);
  TH1D* R4_hist = new TH1D("R4","R4", 100, -1000, 7000);
  TH1D* R5_hist = new TH1D("R5","R5", 100, -1000, 7000);
  TH1D* R6_hist = new TH1D("R6","R6", 100, -1000, 7000);

  TH1D* R1_hist_b = new TH1D("R1b","R1b", 100, -1000, 7000);
  TH1D* R2_hist_b = new TH1D("R2b","R2b", 100, -1000, 7000);
  TH1D* R3_hist_b = new TH1D("R3b","R3b", 100, -1000, 7000);
  TH1D* R4_hist_b = new TH1D("R4b","R4b", 100, -1000, 7000);
  TH1D* R5_hist_b = new TH1D("R5b","R5b", 100, -1000, 7000);
  TH1D* R6_hist_b = new TH1D("R6b","R6b", 100, -1000, 7000);

  TH1D* R1_hist_c = new TH1D("R1c","R1c", 100, -1000, 7000);
  TH1D* R2_hist_c = new TH1D("R2c","R2c", 100, -1000, 7000);
  TH1D* R3_hist_c = new TH1D("R3c","R3c", 100, -1000, 7000);
  TH1D* R4_hist_c = new TH1D("R4c","R4c", 100, -1000, 7000);
  TH1D* R5_hist_c = new TH1D("R5c","R5c", 100, -1000, 7000);
  TH1D* R6_hist_c = new TH1D("R6c","R6c", 100, -1000, 7000);

  TH1D* R1_hist_d = new TH1D("R1d","R1d", 100, -1000, 7000);
  TH1D* R2_hist_d = new TH1D("R2d","R2d", 100, -1000, 7000);
  TH1D* R3_hist_d = new TH1D("R3d","R3d", 100, -1000, 7000);
  TH1D* R4_hist_d = new TH1D("R4d","R4d", 100, -1000, 7000);
  TH1D* R5_hist_d = new TH1D("R5d","R5d", 100, -1000, 7000);
  TH1D* R6_hist_d = new TH1D("R6d","R6d", 100, -1000, 7000);
  gStyle->SetOptStat(0);
  TH1D* scenario_counter_hist = new TH1D("scenerio_counter_hist","scenerio_counter_hist", 6000, 0, 6000);
  TH1D* scenario_basic_hist = new TH1D("scenerio_basic_hist","scenerio_basic_hist", 12, 0, 12);
  
  TH1D* R3_minus_R1_hist_case_1 = new TH1D("R3_minus_R1_hist_case_1","R3_minus_R1_hist_case_1", 100, -1000, 7000);
  TH1D* R6_minus_R1_hist_case_1 = new TH1D("R6_minus_R1_hist_case_1","R6_minus_R1_hist_case_1", 100, -1000, 7000);
  TH1D* R5_minus_R1_hist_case_1 = new TH1D("R5_minus_R1_hist_case_1","R5_minus_R1_hist_case_1", 100, -1000, 7000);

//why not working - misundertanding about how getentry-> works???
   //int entry_final = get_entries-100;
   //tree->GetEntry(entry_final);
   //int final_R0 = tracker_R0->at(1).at(0);
   //cout<<final_R0<<endl;
   //cout<<round(1.1*(final_R0/pow(10,8)))<<endl;
   //int hist_bound = (1.1*(final_R0/pow(10,9)));

  TH2D* prop_time_vs_time_of_hit = new TH2D("prop_time_vs_time_of_hit","prop_time_vs_time_of_hit", 50, 3000, 7000, 250, 0, 10000);

  TH2D* z_vs_prop_time_hist = new TH2D("z_vs_prop_time_hist","z_vs_prop_time_hist", 100, -1, 1, 500, 3000, 7000);
  TH2D* z_vs_prop_time_hist_inner_good_cells = new TH2D("z_vs_prop_time_hist_inner_good_cells","z_vs_prop_time_hist_inner_good_cells", 100, -1, 1, 500, 3000, 7000);
  TH2D* z_vs_t12_t22_hist_inner_good_cells = new TH2D("z_vs_t12_t22_hist_inner_good_cells","z_vs_t12_t22_hist_inner_good_cells", 200, 8, 20, 500, 4, 6);
 

     /////////////////////////////////////////////////////////////////////////////////////////////////
     ////    subtraction histograms - a,b,c,d refers to different timestamp validity scenarios    ////
     /////////////////////////////////////////////////////////////////////////////////////////////////



  TH1D* R1_R3_hist = new TH1D("R1_R3","R1_R3", 300, -5500, -10);
  TH1D* R2_R4_hist = new TH1D("R2_R4","R2_R4", 300,  -5500, -10);

  TH1D* R1_R3_hist_b = new TH1D("R1_R3_b","R1_R3_b", 300,  -5500, -10);
  TH1D* R2_R4_hist_b = new TH1D("R2_R4_b","R2_R4_b", 300,  -5500, -10);

  TH1D* R1_R3_hist_c = new TH1D("R1_R3_c","R1_R3_c", 300,  -5500, -10);
  TH1D* R2_R4_hist_c = new TH1D("R2_R4_c","R2_R4_c", 300, -5500, -10);

  TH1D* R1_R3_hist_d = new TH1D("R1_R3_d","R1_R3_d", 300,  -5500, -10);
  TH1D* R2_R4_hist_d = new TH1D("R2_R4_d","R2_R4_d", 300, -5500, -10);


     /////////////////////////////////////////////////////////////////////////////////////////////////
     ////    subtraction histograms - z refers to zoom - zoomed on one location of histogram    ////
     /////////////////////////////////////////////////////////////////////////////////////////////////




  TH1D* R1_R3_histz = new TH1D("R1_R3z","R1_R3z", 300, -100, 200);
  TH1D* R2_R4_histz = new TH1D("R2_R4z","R2_R4z", 300,  -100, 200);

  TH1D* R1_R3_hist_bz = new TH1D("R1_R3_bz","R1_R3_bz", 300,  -100, 200);
  TH1D* R2_R4_hist_bz = new TH1D("R2_R4_bz","R2_R4_bz", 300,  -100, 200);

  TH1D* R1_R3_hist_cz = new TH1D("R1_R3_cz","R1_R3_cz", 300,  -100, 200);
  TH1D* R2_R4_hist_cz = new TH1D("R2_R4_cz","R2_R4_cz", 300, -100, 200);

  TH1D* R1_R3_hist_dz = new TH1D("R1_R3_dz","R1_R3_dz", 300,  -100, 200);
  TH1D* R2_R4_hist_dz = new TH1D("R2_R4_dz","R2_R4_dz", 300, -100, 200);


     /////////////////////////////////////////////////////////////////////////////////////////////////
     ////    subtraction histograms - prior to timestamp pairing                                   ////
     /////////////////////////////////////////////////////////////////////////////////////////////////




  TH1D* R5_R1_hist = new TH1D("R5_R1","R5_R1", 100, -50, 50);
  TH1D* R5_R2_hist = new TH1D("R5_R2","R5_R2", 100,  -50, 50);
  TH1D* R5_R3_hist = new TH1D("R5_R3","R5_R3", 100,  -150, 0);
  TH1D* R5_R4_hist = new TH1D("R5_R4","R5_R4", 100,  -150, 0);

  TH1D* R6_R1_hist = new TH1D("R6_R1","R6_R1", 100, -50, 50);
  TH1D* R6_R2_hist = new TH1D("R6_R2","R6_R2", 100, -50, 50);
  TH1D* R6_R3_hist = new TH1D("R6_R3","R6_R3", 100, -150, 0);
  TH1D* R6_R4_hist = new TH1D("R6_R4","R6_R4", 100,  -150, 0);


  TH1D* R_hist_diff_R1 = new TH1D("R_hist_diff_R1","R_hist_diff_R1", 100,  -150, 150);
  TH1D* R_hist_diff_R2 = new TH1D("R_hist_diff_R2","R_hist_diff_R2", 100,  -150, 150);
  TH1D* R_hist_diff_R3 = new TH1D("R_hist_diff_R3","R_hist_diff_R3", 100,  -150, 150);
  TH1D* R_hist_diff_R4 = new TH1D("R_hist_diff_R4","R_hist_diff_R4", 100,  -150, 150);









  TH1D* all_cells = new TH1D("all_cells","all_cells", bins, lower_bounds, upper_bounds);
  TH1D* hist_lower = new TH1D("hist_lower","hist_lower", bins, lower_bounds, upper_bounds);
  TH1D* hist_mid = new TH1D("hist_mid","hist_mid", bins, lower_bounds, upper_bounds);
  TH1D* hist_higher = new TH1D("hist_higher","hist_higher", bins, lower_bounds, upper_bounds);

  TH1D* trigger_difference_hist = new TH1D("trigger_difference_hist","trigger_difference_hist", 1000, 0, 2.5);
  TH1D* trigger_difference_hist_zoom = new TH1D("trigger_difference_hist_zoom","trigger_difference_hist_zoom", 50, 0, 0.1);

  TH1D* central_points = new TH1D("central_points","central_points", 100, lower_bounds, upper_bounds);
  TH1D* central_points_edge_cells = new TH1D("central_points_edge_cells","central_points_edge_cells", 100, lower_bounds, upper_bounds);

  TH2D* prop_time_vs_dead_time = new TH2D("prop_time_vs_dead_time","prop_time_vs_dead_time", 10000, 0, 100000000 , 100, 0, 10000);
  TH2D* prop_time_vs_dead_time_replaced_TS = new TH2D("prop_time_vs_dead_time_replaced_TS","prop_time_vs_dead_time_replaced_TS", 500, 0, 3, 500, 0, 6000);
  TH2D* fringes = new TH2D("fringes","fringes", 500, 1000, 5000, 500, 1000, 5000);
  TH2D* fringes_selection = new TH2D("fringes_selection","fringes_selection", 500, 1000, 5000, 500, 1000, 5000);

  TH2D* R5_vs_R6 = new TH2D("R5_vs_R6","R5_vs_R6", 100, 0, 6000, 100, 0, 6000);


  
  TH2D* dead_time_h1_vs_h2_output_hist = new TH2D("dead_time_h1_vs_h2_output_hist","dead_time_h1_vs_h2_output_hist", 100, -1, 1, 100, -1, 1);
  TH2D* dead_time_h1_vs_h2_hist = new TH2D("dead_time_h1_vs_h2_hist","dead_time_h1_vs_h2_hist", 100, -1, 1, 100, -1, 1);
  TH2D* dead_time_h1_vs_h2_counter_hist = new TH2D("dead_time_h1_vs_h2_counter_hist","dead_time_h1_vs_h2_counter_hist", 100, -1, 1, 100, -1, 1);

  TH1D* dead_time_h1_min_h2_output_hist = new TH1D("dead_time_h1_min_h2_output_hist","dead_time_h1_min_h2_output_hist", 100, -1, 1);
  TH1D* dead_time_h1_min_h2_hist = new TH1D("dead_time_h1_min_h2_hist","dead_time_h1_min_h2_hist",  100, -1, 1);
  TH1D* dead_time_h1_min_h2_counter_hist = new TH1D("dead_time_h1_min_h2_counter_hist","dead_time_min_vs_h2_counter_hist",  100, -1, 1);

  TH1D* prop_vs_h_hist = new TH1D("prop_vs_h_hist","prop_vs_h_hist", 100, -1, 1);
  TH1D* prop_vs_h_output_hist = new TH1D("prop_vs_h_output_hist","prop_vs_h_output_hist",  100, -1, 1);
  TH1D* prop_vs_h_counter_hist = new TH1D("prop_vs_h_counter_hist","prop_vs_h_counter_hist",  100, -1, 1);

  TH1D* test_histogram = new TH1D("test_histogram","test_histogram", 100, 0, 10000000);

  TH1D** crate_1_histograms = create_3_hists_for_averaging_function("crate_1",100,-1,1);
  TH1D** crate_2_histograms = create_3_hists_for_averaging_function("crate_2",100,-1,1);
  TH1D** crate_3_histograms = create_3_hists_for_averaging_function("crate_3",100,-1,1);

  TH2D* z_vs_z_improved = new TH2D("z_vs_z_improved","z_vs_z_improved", 1000, -1, 1, 1000, 1, 1);
  TH2D* z_minus_z_improved_2d = new TH2D("z_minus_z_improved_2d","z_minus_z_improved_2d", 100, -1, 1, 100, -0.01, 0.01);
  TH1D* z_minus_z_improved_1d = new TH1D("z_minus_z_improved_1d","z_minus_z_improved_1d", 100, -0.01, 0.01);
///////////////////////////////


////vectors/array generation///
///////////////////////////////

  TH1D** h = new TH1D*[no_cells];
  for (int i=0;i<no_cells;i++) {
      h[i] = new TH1D(TString::Format("h0_%d", i),TString::Format("tracker cell %d", i), bins, lower_bounds, upper_bounds);
                    }

    std::vector<double>** vector_prop_time_vs_dead_time_neighbour = new std::vector<double>*[num_points_dead_time_neighbour];
  for (int i=0;i<num_points_dead_time_neighbour;i++) {
      vector_prop_time_vs_dead_time_neighbour[i] = new std::vector<double>;
                    }

    std::vector<double>** vector_prop_time_vs_dead_time = new std::vector<double>*[num_points];
  for (int i=0;i<num_points;i++) {
      vector_prop_time_vs_dead_time[i] = new std::vector<double>;
  }

    std::vector<double>** vector_R5_vs_R6 = new std::vector<double>*[num_points_R5];
  for (int i=0;i<num_points_R5;i++) {
      vector_R5_vs_R6[i] = new std::vector<double>;
  }

     std::vector<double>** vector_prop_time_vs_dead_time_modified = new std::vector<double>*[num_points];
  for (int i=0;i<num_points_dead_time_modified;i++) {
      vector_prop_time_vs_dead_time_modified[i] = new std::vector<double>;
  }


     std::vector<double>** vector_prop_time_vs_dead_time_half_triggered = new std::vector<double>*[num_points];
  for (int i=0;i<num_points;i++) {
      vector_prop_time_vs_dead_time_half_triggered[i] = new std::vector<double>;
                    }






    std::vector<double>** vector_prop_time_vs_z_outer = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_z_outer[i] = new std::vector<double>;
                    }

     std::vector<double>** vector_prop_time_vs_z_inner = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_z_inner[i] = new std::vector<double>;
                    }


     std::vector<double>** vector_prop_time_vs_z = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_z[i] = new std::vector<double>;
                    }


   std::vector<double>** vector_prop_time_vs_z_near_bad_cell = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_z_near_bad_cell[i] = new std::vector<double>;
                    }

    std::vector<double>** vector_prop_time_vs_z_recent_hit = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_z_recent_hit[i] = new std::vector<double>;
                    }




    std::vector<double>** vector_prop_time_vs_R5 = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_R5[i] = new std::vector<double>;
                    }

     std::vector<double>** vector_prop_time_vs_R6 = new std::vector<double>*[num_points_prop_vs_z];
  for (int i=0;i<num_points_prop_vs_z;i++){
      vector_prop_time_vs_R6[i] = new std::vector<double>;
                    }

    std::vector<double>** vector_prop_time_vs_z_cell_100 = new std::vector<double>*[num_points_prop_vs_z_cell_100];
  for (int i=0;i<num_points_prop_vs_z_cell_100;i++){
      vector_prop_time_vs_z_cell_100[i] = new std::vector<double>;
                    }

    std::vector<double>** vector_prop_time_vs_z_cell_1000 = new std::vector<double>*[num_points_prop_vs_z_cell_1000];
  for (int i=0;i<num_points_prop_vs_z_cell_1000;i++){
      vector_prop_time_vs_z_cell_1000[i] = new std::vector<double>;
                    }

    std::vector<double>** vector_prop_time_vs_z_cell_2000 = new std::vector<double>*[num_points_prop_vs_z_cell_2000];
  for (int i=0;i<num_points_prop_vs_z_cell_2000;i++){
      vector_prop_time_vs_z_cell_2000[i] = new std::vector<double>;
                    }


/////////////////////////////////
///prop_time vector initiation///
/////////////////////////////////


//  TH1D** hist_prop_time_vs_ep = new TH1D*[no_cells];
//  for (int i=0;i<no_cells;i++) {
//      hist_prop_time_vs_ep[i] = new TH1D(TString::Format("h0_%d", i),TString::Format("tracker cell %d", i), bins, lower_bounds, upper_bounds);
//                    }





  string title_root = "/cell_props.root";
  string title_root_2nd_file = "/general_inves.root";
  string folder = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
  //string folder ="/sps/nemo/scratch/spratt/analysis/timestamp_investigation/";
  TFile *newfile = new TFile((folder+run_num+title_root).c_str(), "RECREATE");//saving all histograms
  newfile->cd();


    float       R5_;
    float       R6_;
    

   TTree *tree_new = new TTree("T","CERN 1988 staff data");
   tree_new->Branch("R5_",&R5_,"R5_");
   tree_new->Branch("R6_",&R6_,"R6_");

  int counter_2 =0;

  double R0_previous=0;
  double count_less_than=0;
  double count_total=0;

  double previous_TS[2034] = {0};
  double previous_TS_for_function[2034] = {0};
  double previous_TS_for_function_neighbour[2034] = {0};
  double previous_TS_for_function_modified[2034][3]={0};
  double time_previous_hit;
  double max_dead_time_prop_time_vs_dead_time = 1000000;
  double max_dead_time_prop_time_vs_dead_time_modified = 1000000;
  double max_dead_time_prop_time_vs_dead_time_neighbour = 0.1;
  double counter_all = 1;
  double counter_select = 1;




cout<<get_entries<<endl;



   //int entry_final = get_entries-10;
   //tree->GetEntry(entry_final);
   //int final_R0 = tracker_R0->at(0).at(0);
   //cout<<final_R0<<endl;
   //cout<<round(1.1*(final_R0/pow(10,8)))<<endl;
   //int hist_bound = round(1.1*(final_R0/pow(10,8)));

  

     /////////////////////////////////////////
     ////    LOOPING OVER TRACKER HITS    ////
     /////////////////////////////////////////





   double R0=0;
   double R0_previous__;
   double bot_time ;
   double top_time ;
   double R1;
   double R2;
   double R3;
   double R4;
   double R5;
   double R6;
   double norm_factor_tracker =1;


  for (int i = 0; i < round((get_entries)); i++) {  //looping over events
    //cout<<" -------"<<endl;

    tree->GetEntry(i); // for (int i =0; i< 50;i++){
    for (int j = 0; j < tracker_nohits; j++) {    //looping over tracker hits for each event 
        
        R0_previous__= R0;;

        R0 = tracker_R0->at(j).at(0)*norm_factor_tracker;// assigning new labels for anode and cathode timestamps and 'zeroing' them wrt R0 - should be changed to the calo hit that inititated gathering the data 
        bot_time = tracker_botcathtime->at(j).at(0)*norm_factor_tracker-R0;
        top_time = tracker_topcathtime->at(j).at(0)*norm_factor_tracker-R0;
        R1 = tracker_R1->at(j).at(0)*norm_factor_tracker-R0;
        R2 = tracker_R2->at(j).at(0)*norm_factor_tracker-R0;
        R3 = tracker_R3->at(j).at(0)*norm_factor_tracker-R0;

        R4 = tracker_R4->at(j).at(0)*norm_factor_tracker-R0;
        
        R5 = bot_time;
        R6 = top_time;

        double z=(R5-R6)/(R5+R6);
        double z_improved = improved_z_calc(R5,R6);
        double prop_time = R5+R6;

        double side = tracker_side->at(j);
        double column = tracker_column->at(j);
        double layer = tracker_layer->at(j);
        
        double R0_sec = R0/pow(10,8);
        cell_num = cell_num_calc(side, column, layer);
      

     //////////////////////////////////////////////////
     ////        PROPORTION INVESTIGATION          ////
     //////////////////////////////////////////////////

print_out_timestamp_scenario(R1,R2,R3,R4,R5,R6,scenario_counter_hist,scenario_basic_hist,R3_minus_R1_hist_case_1);

prop_time_vs_time_hit_function(R0,R5,R6,prop_time_vs_time_of_hit,cell_num);

z_vs_prop_time_plot_function(R5,R6,z_vs_prop_time_hist);

mean_prop_time_z_calc_vectors(R5,R6,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9);

prop_time_vs_z_function(vector_prop_time_vs_z,0,R0_sec,R5,R6,num_points_prop_vs_z);




if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){

//cout<<"z:  "<<z<<"z_improved:  "<<z_improved<<endl;
z_vs_z_improved->Fill(z,z_improved,1);
z_minus_z_improved_2d->Fill(z,z-z_improved,1);
z_minus_z_improved_1d->Fill(z-z_improved,1);
}




double previous_ts_hit_modified_R0 = previous_TS_for_function_modified[cell_num][0];
double previous_ts_hit_modified_R5 = previous_TS_for_function_modified[cell_num][1];
double previous_ts_hit_modified_R6 = previous_TS_for_function_modified[cell_num][2];

prop_time_vs_dead_time_function_modified(dead_time_h1_min_h2_hist,dead_time_h1_min_h2_counter_hist,dead_time_h1_vs_h2_hist,dead_time_h1_vs_h2_counter_hist,vector_prop_time_vs_dead_time_half_triggered,vector_prop_time_vs_dead_time_modified,vector_prop_time_vs_dead_time,previous_ts_hit_modified_R0,previous_ts_hit_modified_R5,previous_ts_hit_modified_R6,R0,R5,R6,num_points_dead_time_modified,max_dead_time_prop_time_vs_dead_time_modified,counter_all,counter_select);


prop_time = R5+R6;
double dead_time_ = R0 - previous_ts_hit_modified_R0;
if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
prop_time_vs_dead_time->Fill(dead_time_,prop_time,1);}

previous_TS_for_function_modified[cell_num][0] = R0;
previous_TS_for_function_modified[cell_num][1] = R5; 
previous_TS_for_function_modified[cell_num][2] = R6;

double dead_time = R0 - previous_ts_hit_modified_R0;

///////////think of best way to do this!! 
double previous_ts_hit_neighbour = previous_TS_for_function_neighbour[cell_num];

prop_time_vs_dead_time_function(vector_prop_time_vs_dead_time_neighbour,previous_ts_hit_neighbour,R0_sec,R5,R6,num_points_dead_time_neighbour,max_dead_time_prop_time_vs_dead_time_neighbour);

previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column, layer+1)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column, layer-1)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer-1)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer+1)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer-1)] = R0_sec;
previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer+1)] = R0_sec;







if (dead_time>0.5*pow(10,6)){

   if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
      if (what_is_crate_number(column)==1){
fill_averaging_histograms(crate_1_histograms,z,R5+R6);}

if (what_is_crate_number(column)==2){
fill_averaging_histograms(crate_2_histograms,z,R5+R6);}

if (what_is_crate_number(column)==3){
fill_averaging_histograms(crate_3_histograms,z,R5+R6);}


}






 if (what_is_crate_number(column)==2){

if (is_next_to_bad_cell(side,column,layer)==false && is_cell_bad(side,column,layer)==false){

if (check_cell_on_edge(cell_num)==false){
prop_time_vs_z_function(vector_prop_time_vs_z_inner,0,R0_sec,R5,R6,num_points_prop_vs_z);
prop_vs_h_function(prop_vs_h_hist,prop_vs_h_counter_hist,R5,R6);
z_vs_prop_time_plot_function(R5,R6,z_vs_prop_time_hist_inner_good_cells);  
   if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){ 
      if(z<0.6 && z>0.5){ 
z_vs_t12_t22_hist_inner_good_cells->Fill((pow((R5*0.001),2)+pow((R6*0.001),2)),prop_time*0.001);}}
}


if (check_cell_on_edge(cell_num)==true){
prop_time_vs_z_function(vector_prop_time_vs_z_outer,0,R0_sec,R5,R6,num_points_prop_vs_z);
}

prop_time_vs_z_function_modified(vector_prop_time_vs_R5,vector_prop_time_vs_R6,0,R0_sec,R5,R6,num_points_prop_vs_z);

if (cell_num ==594){

prop_time_vs_z_function(vector_prop_time_vs_z_cell_100,0,R0_sec,R5,R6,num_points_prop_vs_z_cell_100);
}

if (cell_num ==500){
prop_time_vs_z_function(vector_prop_time_vs_z_cell_1000,0,R0_sec,R5,R6,num_points_prop_vs_z_cell_1000);
}

if (cell_num ==600){
prop_time_vs_z_function(vector_prop_time_vs_z_cell_2000,0,R0_sec,R5,R6,num_points_prop_vs_z_cell_2000);
}

R5_vs_R6_function(vector_R5_vs_R6,R5,R6,num_points_R5,7000);

R5_vs_R6_hist_function(R5_vs_R6,R5,R6); 

}


if (is_next_to_bad_cell(side,column,layer)==true && is_cell_bad(side,column,layer)==false){

prop_time_vs_z_function(vector_prop_time_vs_z_near_bad_cell,0,R0_sec,R5,R6,num_points_prop_vs_z);

}

}


if (dead_time<0.5*pow(10,6)){
prop_time_vs_z_function(vector_prop_time_vs_z_recent_hit,0,R0_sec,R5,R6,num_points_prop_vs_z);
}}



if (R5<0 && R6>0){
R1_hist->Fill(R1,1);
R2_hist->Fill(R2,1);
R3_hist->Fill(R3,1);

R4_hist->Fill(R4,1);

R1_R3_hist->Fill(R3-R1,1);
R2_R4_hist->Fill(R4-R2,1);

R1_R3_histz->Fill(R3-R1,1);
R2_R4_histz->Fill(R4-R2,1);
}

if (R5>0 && R6<0){
R1_hist_b->Fill(R1,1);
R2_hist_b->Fill(R2,1);
R3_hist_b->Fill(R3,1);
R4_hist_b->Fill(R4,1);

R1_R3_hist_b->Fill(R3-R1,1);
R2_R4_hist_b->Fill(R4-R2,1);

R1_R3_hist_bz->Fill(R3-R1,1);
R2_R4_hist_bz->Fill(R4-R2,1);
}

if (R5>0 && R6>0){
R1_hist_c->Fill(R1,1);
R2_hist_c->Fill(R2,1);
R3_hist_c->Fill(R3,1);
R4_hist_c->Fill(R4,1);

R5_R1_hist->Fill(R5-R1,1);
R5_R2_hist->Fill(R5-R2,1);
R5_R3_hist->Fill(R5-R3,1);
R5_R4_hist->Fill(R5-R4,1);

R6_R1_hist->Fill(R6-R1,1);
R6_R2_hist->Fill(R6-R2,1);
R6_R3_hist->Fill(R6-R3,1);
R6_R4_hist->Fill(R6-R4,1);



int pair = pair_timestamps(R5,R6,R1);
if (pair == 1){
R_hist_diff_R1->Fill(R5-R1);
R_hist_diff_R3->Fill(R5-R3);
}
if (pair == 2){
R_hist_diff_R1->Fill(R6-R1);
R_hist_diff_R3->Fill(R6-R3);
}

int pair_2 = pair_timestamps(R5,R6,R2);
if (pair_2 == 1){
R_hist_diff_R2->Fill(R5-R2);
R_hist_diff_R4->Fill(R5-R4);
}
if (pair_2 == 2){
R_hist_diff_R2->Fill(R6-R2);
R_hist_diff_R4->Fill(R6-R4);
}

}



if (R5<0 && R6<0){
R1_hist_d->Fill(R1,1);
R2_hist_d->Fill(R2,1);
R3_hist_d->Fill(R3,1);
R4_hist_d->Fill(R4,1);

R1_R3_hist_d->Fill(R3-R1,1);
R2_R4_hist_d->Fill(R4-R2,1);

R1_R3_hist_dz->Fill(R3-R1,1);
R2_R4_hist_dz->Fill(R4-R2,1);
}




R5_hist->Fill(R5,1);
R6_hist->Fill(R6,1);






     //////////////////////////////////////////////////
     ////   DEADTIME AND PROPTIME INVESTIGATION    ////
     //////////////////////////////////////////////////


if (R5>0 && R6>0 && R5<100000 && R6<100000){

  
cell_num = cell_num_calc(side, column, layer);
array_check->fill(cell_num,1);

prop_time = R5+R6;

counter +=1;
all_cells->Fill(prop_time,1);
h[cell_num]->Fill(prop_time,1); 
tree_new->Fill();

time_previous_hit = previous_TS[cell_num];

 // cout << R0_sec <<" seconds" <<endl;
 trigger_difference_hist->Fill(R0_sec-time_previous_hit,1);
 trigger_difference_hist_zoom->Fill(R0_sec-time_previous_hit,1);
 

if (R5>R1 && R1>0 && R1<100000 && abs(R1-R5)>200){
 prop_time_vs_dead_time_replaced_TS->Fill(R0_sec-time_previous_hit,R5+R1);
 fringes->Fill(R1,R5,1);
 if (R0_sec-time_previous_hit<0.05 ){

    fringes_selection->Fill(R1,R5,1);

 }
}

previous_TS[cell_num] = R0_sec;


}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }}
  //cout<<"---------"<<endl;
  //cout<<counter_select/counter_all<<endl;
  //cout<<counter_all<<"     "<<counter_select<<endl;

       /////////////////////////////////////////
      ////    END OF TRACKER HIT LOOP     ////
     /////////////////////////////////////////


























       /////////////////////////////////////////
      ////    2ND TRACKER HIT LOOP     ////
     /////////////////////////////////////////

for (int i = 0; i < round((get_entries)); i++) {  //looping over events
    std::vector<int> *cell_nums    = new std::vector<int>;
    tree->GetEntry(i); // for (int i =0; i< 50;i++){
    
    
    for (int j = 0; j < tracker_nohits; j++) {    //lo

     cell_nums->push_back(cell_num_calc(tracker_side->at(j),tracker_column->at(j),tracker_layer->at(j)));

    }
    
     for (int u=0;u<cell_nums->size();u++){

      
        for (int k=0;k<cell_nums->size();k++){
             if (cell_nums->at(k) == cell_nums->at(u) && k!=u){
               for (int h=0;h<cell_nums->size();h++){

       // cout << cell_nums->at(h) << " "<<endl;
             }

        }

    }
    
    
    }

}


















    /////////////////////////////////////////
    ////      FITTING FOR EACH CELL      ////
    /////////////////////////////////////////

    for (int i=0;i<no_cells;i++) {
    if (h[i]->GetRMS()>0){
                    
                    

      //perform_fit_2(h[i],i,bins, lower_bounds, upper_bounds);
      categorise_peaks(h[i],hist_lower,hist_mid,hist_higher,i,bins, lower_bounds, upper_bounds);
      
      centre_point = (h[i]->GetMaximumBin()*(upper_bounds-lower_bounds)/bins)+lower_bounds;
      if (centre_point>3500 && centre_point<5500){
      snd_prop->setcontent(i,centre_point);
      }
      central_points->Fill(centre_point,1);

    bool cell_on_edge = check_cell_on_edge(i); 
      if (cell_on_edge ==true){
        snd_check->fill(i,1);
        central_points_edge_cells->Fill(centre_point,1);
      }

      h[i]->GetXaxis()->SetTitle("Propagation Time");
      h[i]->Write();            
    
        }
    }



       


TGraph* test_plot = new TGraph(10000);





for (int i=0;i<102;i++){
   double vel = 0.005;
   double acc=0.00001
   

   ;
   double d = i*0.01;

   
   double eps = calc_eps(d,vel,acc);



   double tp = calc_tp(d,vel,acc);

   test_plot->SetPoint(i, eps, tp);
   //cout<<eps<<"         "<<tp<<"    "<<d<<endl;

}










       
    /////////////////////////////////////////
    ////      PLOT SETUP      ////
    /////////////////////////////////////////


         TCanvas *c12=new TCanvas();

         all_cells->SetLineColor(kRed);
         hist_lower->SetLineColor(kBlue);
         hist_mid->SetLineColor(kGreen);
         hist_higher->SetLineColor(kOrange);

         all_cells->Draw();
         hist_lower->Draw("same");
         hist_higher->Draw("same");
         hist_mid->Draw("same");
         c12->Write();



   
   all_cells->GetXaxis()->SetTitle("Propagation Time");
   all_cells->Write(); 

   central_points->GetXaxis()->SetTitle("Propagation Time Centre Point of Each Cell");
   central_points->Write(); 

   central_points_edge_cells->GetXaxis()->SetTitle("Propagation Time Centre Point of Each Cell");
   central_points_edge_cells->Write(); 

   hist_lower->Write();
   hist_mid->Write();
   hist_higher->Write(); 

   trigger_difference_hist->Write(); 
   trigger_difference_hist_zoom->Write();
   prop_time_vs_dead_time_replaced_TS->Write();
   fringes->Write();
   fringes_selection->Write();

   tree_new->Write();
   newfile->Close();

  
  

     // TCanvas *c12=new TCanvas();
  // hs->Draw();

  TCanvas* stacked_ts_canvas = new TCanvas("stacked_ts_canvas","stacked_ts_canvas",600,500);  
  stack_histograms_2(stacked_ts_canvas,R1_hist,R1_hist_b,R1_hist_c,R1_hist_d,R2_hist,R2_hist_b,R2_hist_c,R2_hist_d,R3_hist,R3_hist_b,R3_hist_c,R3_hist_d,R4_hist,R4_hist_b,R4_hist_c,R4_hist_d);

  TCanvas* stacked_ts_difference_canvas_1 = new TCanvas("stacked_ts_difference_canvas_1","stacked_ts_difference_canvas_1",600,500);  
  stack_histograms_3(stacked_ts_difference_canvas_1,R1_R3_hist,R1_R3_hist_b,R1_R3_hist_c,R1_R3_hist_d,R2_R4_hist,R2_R4_hist_b,R2_R4_hist_c,R2_R4_hist_d);

  TCanvas* stacked_ts_difference_canvas_2 = new TCanvas("stacked_ts_difference_canvas_2","stacked_ts_difference_canvas_2",600,500);  
  stack_histograms_3(stacked_ts_difference_canvas_2,R1_R3_histz,R1_R3_hist_bz,R1_R3_hist_cz,R1_R3_hist_dz,R2_R4_histz,R2_R4_hist_bz,R2_R4_hist_cz,R2_R4_hist_dz);

  TCanvas* R5_ts_difference_canvas = new TCanvas("R5_ts_difference_canvas","R5_ts_difference_canvas",600,500); 
  stack_histograms_4(R5_ts_difference_canvas,R5_R1_hist,R5_R2_hist,R5_R3_hist,R5_R4_hist);

  TCanvas* R6_ts_difference_canvas = new TCanvas("R6_ts_difference_canvas","R6_ts_difference_canvas",600,500); 
  stack_histograms_4(R6_ts_difference_canvas,R6_R1_hist,R6_R2_hist,R6_R3_hist,R6_R4_hist);

  TCanvas* paired_ts_difference_canvas = new TCanvas("paired_ts_difference_canvas","paired_ts_difference_canvas",600,500); 
  stack_histograms_4(paired_ts_difference_canvas,R_hist_diff_R1,R_hist_diff_R2,R_hist_diff_R3,R_hist_diff_R4);
 
  


   //Double_t sum = scenario_basic_hist->Integral();

   // Scale the histogram by dividing each bin content by the total sum
   //scenario_basic_hist->Scale(1.0 / sum);

   // Set the y-axis range to accommodate the scaled histogram
   //Double_t maxY = scenario_basic_hist->GetMaximum();
   //scenario_basic_hist->GetYaxis()->SetRangeUser(0, maxY * 1.1); // Adjust the factor (1.1) as needed






  TFile *newnewfile = new TFile((folder+run_num+title_root_2nd_file).c_str(), "RECREATE");//saving all histograms
  newnewfile->cd();

   TCanvas *c15=new TCanvas();
   scenario_basic_hist->Draw();
   PaintBin(scenario_basic_hist, 1, 38);

   PaintBin(scenario_basic_hist, 5, 30);
   PaintBin(scenario_basic_hist, 6, 30);
   PaintBin(scenario_basic_hist, 7, 30);
   PaintBin(scenario_basic_hist, 8, 30);
   PaintBin(scenario_basic_hist, 9, 30);
   PaintBin(scenario_basic_hist, 10, 30);
   PaintBin(scenario_basic_hist, 11, 30);

   PaintBin(scenario_basic_hist, 2, 25);
   PaintBin(scenario_basic_hist, 3, 25);
   PaintBin(scenario_basic_hist, 4, 25);


   PaintBin(scenario_basic_hist, 12, 35);

  ///////////////////////////////////
  //////proptime vs z tgraphs////////
  ///////////////////////////////////

  graph_function(mean_prop_time_vs_z_cell_100_graph,vector_prop_time_vs_z_cell_100,num_points_prop_vs_z_cell_100,"mean_prop_time_vs_z_cell_100_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_cell_1000_graph,vector_prop_time_vs_z_cell_1000,num_points_prop_vs_z_cell_1000,"mean_prop_time_vs_z_cell_1000_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_cell_2000_graph,vector_prop_time_vs_z_cell_2000,num_points_prop_vs_z_cell_2000,"mean_prop_time_vs_z_cell_2000_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_outer_graph,vector_prop_time_vs_z_outer,num_points_prop_vs_z,"mean_prop_time_vs_z_outer_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_inner_graph,vector_prop_time_vs_z_inner,num_points_prop_vs_z,"mean_prop_time_vs_z_inner_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_graph,vector_prop_time_vs_z,num_points_prop_vs_z,"mean_prop_time_vs_z_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_near_bad_cell_graph,vector_prop_time_vs_z_near_bad_cell,num_points_prop_vs_z,"mean_prop_time_vs_z_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_z_recent_hit_graph,vector_prop_time_vs_z_recent_hit,num_points_prop_vs_z,"mean_prop_time_vs_z_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  
   //////////////////////////
  //////Other tgraphs////////
  //////////////////////////

  graph_function(mean_prop_time_vs_R5_graph,vector_prop_time_vs_R5,num_points_prop_vs_z,"R5_vs_z_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(mean_prop_time_vs_R6_graph,vector_prop_time_vs_R6,num_points_prop_vs_z,"R6_vs_z_graph","z position (proportion)","total prop time (10-4 sec)",1,1);
  graph_function(prop_time_vs_dead_time_graph,vector_prop_time_vs_dead_time,num_points,"prop_time_vs_dead_time_graph","deadtime of cell (sec)","total prop time (sec)",max_dead_time_prop_time_vs_dead_time,pow(10,-8));
  graph_function(prop_time_vs_dead_time_neighbour_graph,vector_prop_time_vs_dead_time_neighbour,num_points_dead_time_neighbour,"prop_time_vs_dead_time_neighbour_graph","deadtime of latest neighbour cell (10-4 sec)","total prop time (10-4 sec)",max_dead_time_prop_time_vs_dead_time_neighbour,1);
  graph_function(prop_time_vs_dead_time_modified_graph,vector_prop_time_vs_dead_time_modified,num_points,"prop_time_vs_dead_time_modified_graph","deadtime (sec)","total prop time (sec)",max_dead_time_prop_time_vs_dead_time_modified,pow(10,-8));
  graph_function(prop_time_vs_dead_time_half_triggered_graph,vector_prop_time_vs_dead_time_half_triggered,num_points,"prop_time_vs_dead_time_modified_graph","deadtime (sec)","total prop time (sec)",max_dead_time_prop_time_vs_dead_time_modified,pow(10,-8));
  graph_function(R5_vs_R6_graph,vector_R5_vs_R6,num_points_R5,"R5 vs R6","R5","R6",5500,pow(10,-8));
     
  graph_function_plasma_speed_plot(R5_vs_R6_graph, plasma_speed_graph);

  graph_function_2(hist_z_no,vector_prop_time_vs_z,num_points_prop_vs_z,"no events for different epsilon values","z value (proportion)","no events",1);

  graph_function_split_graph_and_flip_x_y(mean_prop_time_vs_z_graph,mean_prop_time_vs_z_graph_lower_z,mean_prop_time_vs_z_graph_upper_z);

  TF1 *fitFunction = new TF1("fitFunction", "(3/[0]) -   x   + ([1]/[0])*x*x    +     [1]*(1/(8*pow([0],3)))*pow((    -3  +    (x*[0])    -  [1]*x*x  ) ,2 ) "  , 1,3);
  fit_distance_function_to_t_graph(plasma_speed_graph,fitFunction);

  //TF1 *fitFunction_quadratic = new TF1("fitFunction_quadratic_name", quadratic_fitting_function ,0.1,0.9,3);
  //fit_quadratic(mean_prop_time_vs_z_graph,fitFunction_quadratic);

  TF1 *fitFunction_iteration_eps_vs_z = new TF1("fitFunction_iteration_eps_vs_z", get_tp_given_eps ,0.15,0.85,2);
  fit_iteration_eps_vs_z(mean_prop_time_vs_z_graph,fitFunction_iteration_eps_vs_z);

  TF1 *fitFunction_iteration_eps_vs_z_inner_cells = new TF1("fitFunction_iteration_eps_vs_z_inner_cells", get_tp_given_eps ,0.15,0.85,2);
  fit_iteration_eps_vs_z(mean_prop_time_vs_z_inner_graph,fitFunction_iteration_eps_vs_z_inner_cells);

  

   for (int i=0;i<100;i++){
      for (int j=0;j<100;j++){
   double no_content = dead_time_h1_vs_h2_counter_hist->GetBinContent(i,j);
   double total_content =dead_time_h1_vs_h2_hist->GetBinContent(i,j);
   if (no_content>0){
      double set_value =total_content/no_content;
      //cout<<set_value<<"   "<<total_content<<"      "<<no_content<<endl;
   dead_time_h1_vs_h2_output_hist->SetBinContent(i,j,set_value);
   }}}

 for (int i=0;i<101;i++){
   double no_content = dead_time_h1_min_h2_counter_hist->GetBinContent(i);
   double total_content =dead_time_h1_min_h2_hist->GetBinContent(i);
   if (no_content>0){
      double set_value =total_content/no_content;
      //cout<<set_value<<"   "<<total_content<<"      "<<no_content<<endl;
   dead_time_h1_min_h2_output_hist->SetBinContent(i,set_value);
   }}


//generalise this fitting function code in future so fits saves and does everything in one big function! this is definitly possible !

  TF1 *fitFunction_prop_vs_eps_lower = new TF1("prop_vs_eps_fit_lower_z",prop_vs_eps_fit_lower,4,5,2);
  fit_plasma_prop_time_vs_eps(mean_prop_time_vs_z_graph_lower_z,fitFunction_prop_vs_eps_lower);


  TF1 *fitFunction_prop_vs_eps_upper = new TF1("prop_vs_eps_fit_upper_z",prop_vs_eps_fit_upper,4,5,2);
  fit_plasma_prop_time_vs_eps(mean_prop_time_vs_z_graph_upper_z,fitFunction_prop_vs_eps_upper);



   z_vs_z_improved->Write();
   z_minus_z_improved_1d->Write();
   z_minus_z_improved_2d->Write();

   z_vs_prop_time_hist_inner_good_cells->Write();
   z_vs_t12_t22_hist_inner_good_cells->Write();
   mean_prop_time_vs_z_inner_graph->Write();
   test_plot->Write();

   dead_time_h1_min_h2_output_hist->Write();
   dead_time_h1_vs_h2_output_hist->Write();
   dead_time_h1_vs_h2_counter_hist->Write();

   mean_prop_time_vs_z_graph_lower_z->Write();
   mean_prop_time_vs_z_graph_upper_z->Write();
   TCanvas* plasma_vel_canvas = new TCanvas("plasma_vel","plasma_vel");  
   fitFunction->SetLineColor(kBlue);
   plasma_speed_graph->Draw();
   fitFunction->Draw("Same");
   double xMin=0;
   double xMax=6;
   plasma_vel_canvas->Draw();
   plasma_vel_canvas->Write();
   fitFunction->Write();



   average_over_hist(prop_vs_h_hist,prop_vs_h_counter_hist,prop_vs_h_output_hist,100);
   average_over_hist(crate_1_histograms);
   average_over_hist(crate_2_histograms);
   average_over_hist(crate_3_histograms);

   TCanvas* crate_histograms_canvas = new TCanvas("crate_histograms","crate_histograms");  
   crate_1_histograms[2]->Draw();
   crate_2_histograms[2]->Draw("Same");
   crate_3_histograms[2]->Draw("Same");
   crate_histograms_canvas->Draw();
   crate_histograms_canvas->Write();

   crate_1_histograms[2]->Write();
   crate_2_histograms[2]->Write();
   crate_3_histograms[2]->Write();




   TH1D* prop_time_assymetry = new TH1D("prop_time_assymetry","prop_time_assymetry", prop_vs_h_output_hist->GetNbinsX()/2, 0, 1);
   prop_time_asymmetry_plot(prop_vs_h_output_hist,prop_time_assymetry,prop_vs_h_counter_hist);

   prop_time_assymetry->Write();
   prop_vs_h_output_hist->Write();

   TCanvas* plasma_vs_eps_canvas_lower = new TCanvas("plasma_vs_eps_canvas_lower","plasma_vs_eps_canvas_lower");  
   fitFunction_prop_vs_eps_lower->SetLineColor(kBlue);
   mean_prop_time_vs_z_graph_lower_z->Draw();
   fitFunction_prop_vs_eps_lower->Draw("Same");
   plasma_vs_eps_canvas_lower->Draw();
   plasma_vs_eps_canvas_lower->Write();
   fitFunction_prop_vs_eps_lower->Write();


   TCanvas* plasma_vs_eps_canvas_upper = new TCanvas("plasma_vs_eps_canvas_upper","plasma_vs_eps_canvas_upper");  
   fitFunction_prop_vs_eps_upper->SetLineColor(kBlue);
   mean_prop_time_vs_z_graph_upper_z->Draw();
   fitFunction_prop_vs_eps_upper->Draw("Same");
   plasma_vs_eps_canvas_upper->Draw();
   plasma_vs_eps_canvas_upper->Write();
   fitFunction_prop_vs_eps_upper->Write();


   TCanvas* fit_quadratic_canvas = new TCanvas("fit_quadratic_canvas","fit_quadratic_canvas");  
   fitFunction_iteration_eps_vs_z->SetLineColor(kBlue);
   mean_prop_time_vs_z_graph->Draw();
   fitFunction_iteration_eps_vs_z->Draw("Same");
   fit_quadratic_canvas->Draw();
   fit_quadratic_canvas->Write();
   fitFunction_iteration_eps_vs_z->Write();



   TCanvas* prop_time_inner_vs_outer_cells = new TCanvas("prop_time_inner_vs_outer","prop_time_inner_vs_outer");  
   mean_prop_time_vs_z_outer_graph->SetLineColor(kBlue);
   mean_prop_time_vs_z_outer_graph->SetMarkerColor(kBlue);
   mean_prop_time_vs_z_inner_graph->SetLineColor(kGreen);
   mean_prop_time_vs_z_inner_graph->SetMarkerColor(kGreen);
   mean_prop_time_vs_z_near_bad_cell_graph->SetLineColor(kRed);
   mean_prop_time_vs_z_near_bad_cell_graph->SetMarkerColor(kRed);
   mean_prop_time_vs_z_recent_hit_graph->SetMarkerColor(kOrange);
   mean_prop_time_vs_z_recent_hit_graph->SetLineColor(kOrange);


   mean_prop_time_vs_z_recent_hit_graph->Draw();
   mean_prop_time_vs_z_outer_graph->Draw("Same");

   mean_prop_time_vs_z_inner_graph->Draw("Same");
   mean_prop_time_vs_z_near_bad_cell_graph->Draw("Same");
   mean_prop_time_vs_z_recent_hit_graph->Draw("Same");

   auto legend = new TLegend(0.1,0.1);
   legend->AddEntry(mean_prop_time_vs_z_outer_graph,"outer cells");
   legend->AddEntry(mean_prop_time_vs_z_inner_graph,"inner cells");
   legend->AddEntry(mean_prop_time_vs_z_near_bad_cell_graph,"near bad cells");
   legend->AddEntry(mean_prop_time_vs_z_recent_hit_graph,"recent hit");
   legend->Draw("Same");
   prop_time_inner_vs_outer_cells->Update();
   prop_time_inner_vs_outer_cells->Draw();
   prop_time_inner_vs_outer_cells->Write();




  R5_vs_R6->Write();
  //plasma_speed_graph->Write();

  R5_vs_R6_graph->Write();
  prop_time_vs_dead_time->Write();
  prop_time_vs_dead_time_modified_graph->Write();
  prop_time_vs_dead_time_half_triggered_graph->Write();
  stacked_ts_canvas->Write();
  stacked_ts_difference_canvas_1->Write();
  stacked_ts_difference_canvas_2->Write();
  R5_ts_difference_canvas->Write();
  R6_ts_difference_canvas->Write();

  paired_ts_difference_canvas->Write();
  scenario_counter_hist->Write();
  scenario_basic_hist->Write();
  R3_minus_R1_hist_case_1->Write();
  c15->Write();
  prop_time_vs_time_of_hit->Write();
  z_vs_prop_time_hist->Write();
  mean_prop_time_vs_z_graph->Write();
  mean_prop_time_vs_R5_graph->Write();
  mean_prop_time_vs_R6_graph->Write();
  prop_time_vs_dead_time_graph->Write();
  prop_time_vs_dead_time_neighbour_graph->Write();
  hist_z_no->Write();

  mean_prop_time_vs_z_cell_100_graph->Write();
  mean_prop_time_vs_z_cell_1000_graph->Write();
  mean_prop_time_vs_z_cell_2000_graph->Write();

  

   //TCanvas* c1 = new TCanvas("c1","c1",600,500);
   //gStyle->SetOptStat(0);

  newnewfile->Close();



   save_snd(run_num,snd_prop,"prop",3500,5500);
   save_snd(run_num,snd_check,"check",0,3);   
   save_snd(run_num,bad_cells,"bad_cells",0,2); 
   save_snd(run_num,array_check,"array_check",0,10); 
   save_snd(run_num,activity_snd,"activity_snd",0,1000); 
   
   

  return 0;
}












