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
#include "/sps/nemo/scratch/spratt/analysis/analysis_library.h"
using namespace std;











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



bool check_ts_valid(double R){
               bool output = true;
               if (R<0 || R>100000){
                  output=false;
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
  
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/"+run_num+"/individual_cell_data/";

    return directory;
}






///////////////////////////////////














void average_over_hist(TH1D* h_sum,TH1D* h_counter,TH1D* h_output){
   
   double number_bins = h_sum->GetNbinsX();

	for(int i=0;i<number_bins+1;i++){

	double total = h_sum->GetBinContent(i);
	double no = h_counter->GetBinContent(i);
	if (no>0){
	double set_value = total/no;
	//cout<<set_value<<endl;
	h_output->SetBinContent(i,set_value);
	}}}


void fill_averaging_histograms(TH1D* h_sum,TH1D* h_count,double x,double y){

      h_sum->Fill(x,y);
      h_count->Fill(x,1);

}


void average_all_cell_hists(TH1D** h_sum,TH1D** h_counter,TH1D** h_output){

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
average_over_hist(h_sum[i],h_counter[i],h_output[i]);
}
}


void add_cell_description(TH1D* hist,TCanvas* canvas,int i){

string edge_cell="NO";
string bad_cell="NO";
string near_bad_cell="NO";
string crate=to_string(what_is_crate_number(i));

if (check_cell_on_edge(i)==true){
edge_cell = "YES";
}
if (is_cell_bad(i)==true){
bad_cell = "YES";
}
if (is_next_to_bad_cell(i)==true){
near_bad_cell = "YES";

}

TPaveText *t = new TPaveText(0.9,0.9,0.9,0.9);
   t->AddText(("Edge Cell: "+edge_cell).c_str());
   t->AddText(("Bad Cell: "+bad_cell).c_str());
   t->AddText(("Near Bad Cell: "+near_bad_cell).c_str());
   t->AddText(("Crate: "+crate).c_str());
  


canvas->cd();
hist->Draw();

t->Draw("Same");

}


void make_dir(string run_num){

    string folder = name_folder(run_num);
    if(mkdir((folder).c_str(), 0777) == -1);
}

void make_root_direc(string name,TFile *file){

  TDirectory *direc = file->mkdir(name.c_str());
    direc->cd();

}


TFile* make_root_file(string name,string run_num){

  string folder = name_folder(run_num);
  TFile *rootfile = new TFile((folder+name+".root").c_str(), "RECREATE");//saving all histograms
  rootfile->cd();

  return rootfile;

}


void save_all_cell_hists(TH1D** output_hists,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
output_hists[i]->Write();
}
root_file->Close();
}


void save_all_cell_hists(TH2D** output_hists,string name,string run_num){

TFile* root_file =  make_root_file(name,run_num);

int num_cells = 2034;
for (int i=0;i<num_cells;i++){
output_hists[i]->Write();
}
root_file->Close();
}  







int main() {


//std::string file_1 ="/sps/nemo/scratch/chauveau/commissioning/sam/snemo_run-1544_track-udd.root";
//
string run_num = "1051";
cout<<run_num<<endl;

//std::string file_1 = "/sps/nemo/scratch/golivier/aussois_tuto_data/bi207_run/snemo_run-840_udd.root";//Chosen run file 
//std::string file_1 ="/sps/nemo/scratch/golivier/data/UDD/snemo_run-840_udd.root";
//std::string file_1 ="/sps/nemo/scratch/golivier/data/UDD/snemo_run-974_udd.root";
//std::string file_1 ="/sps/nemo/scratch/spratt/analysis/UDD_Data/snemo_run-974_udd.root";
//std::string file_1 ="/sps/nemo/scratch/spratt/analysis/UDD_Data/snemo_run-1018_udd.root";
std::string file_1 ="/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/snemo_run-1051_udd.root";


make_dir(run_num);
TFile *file= new TFile(file_1.c_str(), "READ");

  std::vector<std::vector<short>> *wave = new std::vector<std::vector<short>>;//initiating vectors that will be filled with data from the event file
  std::vector<int> *calo_wall   = new std::vector<int>;
  std::vector<int> *calo_side   = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row    = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;
  std::vector<int> *calo_ampl   = new std::vector<int>;
  std::vector<long> *calo_ts     = new std::vector<long>;
  std::vector<int> *tracker_side   = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_row    = new std::vector<int>;
  std::vector<int> *tracker_layer  = new std::vector<int>;
  std::vector<std::vector<long> >*tracker_botcathtime  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_topcathtime  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_R0  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_R1  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_R2  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_R3  = new std::vector<std::vector<long> >;
  std::vector<std::vector<long> > *tracker_R4  = new std::vector<std::vector<long> >;

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
  tree->SetBranchStatus("digicalo.timestamp",1);
  tree->SetBranchAddress("digicalo.timestamp", &calo_ts);

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

   TH1D** sum_hists = make_array_of_hist("sum",100,-1,1);   
   TH1D** counter_hists = make_array_of_hist("counter",100,-1,1);
   TH1D** output_hists = make_array_of_hist("output",100,-1,1);

   TH2D** dead_time_vs_prop_time = make_array_of_hist_2d("dead_time_vs_prop_time",10000, 0, 100000000 , 100, 3000, 7000);
   TH2D** dead_time_vs_prop_time_neighbour = make_array_of_hist_2d("dead_time_vs_prop_time_neighbour",1000, 0, 10000 , 100, 3000, 7000);
   TH2D** full_2d_heatmap_prop_vs_z = make_array_of_hist_2d("full_2d_heatmap_prop_vs_z",100,-1,1,1000,3000,7000);
   TH1D** prop_time_over_run = make_array_of_hist_2d("prop_time_over_run",1000,3000,7000,100,start_R0,end_R0);


   double get_entries = (tree->GetEntries() - 1);
   double R0;
   double bot_time ;
   double top_time ;

   double R1;
   double R2;
   double R3;

   double R4;
   double R5;
   double R6;
   double RC;
   double RC_zero;
   double t_drift;
   double t_drift_0;
   double t_drift_1;

   double previous_TS_for_function_modified[2034][3] = {0};
   double previous_TS_for_function_neighbour[2034]   = {0};
   double time_previous_hit;

   double side;

   double column;
   double layer;

   int cell_num;

   double norm_factor_calo = 1;//6.25;
   double norm_factor_tracker =1;//12.5;

   cout<<run_num<<endl;

   cout<<get_entries<<endl;

   double end_R0=0;
   double end_R0_test;

   for (int i=round(get_entries*0.97);i<get_entries;i++){
     // cout<<i<<endl;
   tree->GetEntry(i); 
    for(int j=0;j<tracker_nohits;j++){
   end_R0_test = tracker_R0->at(j).at(0)*norm_factor_tracker;
   //cout<<end_R0_test<<endl;
   if (end_R0_test>end_R0){
           end_R0=end_R0_test;
   }}}

    cout<<end_R0<<endl;
   
   tree->GetEntry(0); 
   int start_R0 = (tracker_R0->at(0).at(0)*norm_factor_tracker);
   cout<<start_R0<<endl;

   
  for (int i = 0; i < round((get_entries)); i++) {  //looping over events
    tree->GetEntry(i); 

    for (int j = 0; j < tracker_nohits; j++) {    //looping over tracker hits for each event 
        
        R0 = tracker_R0->at(j).at(0)*norm_factor_tracker;// assigning new labels for anode and cathode timestamps and 'zeroing' them wrt R0 - should be changed to the calo hit that inititated gathering the data 
        bot_time = tracker_botcathtime->at(j).at(0)*norm_factor_tracker;
        top_time = tracker_topcathtime->at(j).at(0)*norm_factor_tracker;
        R1 = tracker_R1->at(j).at(0)*norm_factor_tracker;
        R2 = tracker_R2->at(j).at(0)*norm_factor_tracker;
        R3 = tracker_R3->at(j).at(0)*norm_factor_tracker;
        R4 = tracker_R4->at(j).at(0)*norm_factor_tracker;
        R5 = bot_time;
        R6 = top_time;
        
     //   cout<<R0<<endl;
        double calo_time_correction_term = round((R0-RC)/(2.6843*pow(10,10)));
        //RC = RC + calo_time_correction_term*(2.6843*pow(10,10));

        R1=R1-R0;
        R2=R2-R0;
        R3=R3-R0;
        R4=R4-R0;
        R5=R5-R0;
        R6=R6-R0;

        side = tracker_side->at(j);
        column = tracker_column->at(j);
        layer = tracker_layer->at(j);
        
        cell_num = cell_num_calc(side, column, layer);
        double z = (R5-R6)/(R5+R6);
        double prop_time = R5+R6;
        if (check_ts_valid(R5)==true && check_ts_valid(R6)==true){
                  
                  double previous_ts_hit_modified_R0 = previous_TS_for_function_modified[cell_num][0];
                  double previous_ts_hit_modified_R5 = previous_TS_for_function_modified[cell_num][1];
                  double previous_ts_hit_modified_R6 = previous_TS_for_function_modified[cell_num][2];
                  double previous_ts_hit_neighbour_R0 = previous_TS_for_function_neighbour[cell_num];

                  double dead_time = R0 - previous_ts_hit_modified_R0;
                  double dead_time_neighbour = R0 - previous_ts_hit_neighbour_R0;
               
                  dead_time_vs_prop_time[cell_num]->Fill(dead_time,prop_time,1);
               if (dead_time>0.7*pow(10,6)*norm_factor_tracker){
                  dead_time_vs_prop_time_neighbour[cell_num]->Fill(dead_time_neighbour,prop_time,1);}

                  previous_TS_for_function_modified[cell_num][0] = R0;
                  previous_TS_for_function_modified[cell_num][1] = R5; 
                  previous_TS_for_function_modified[cell_num][2] = R6;

                  previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column, layer+1)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column, layer-1)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer-1)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column+1, layer+1)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer-1)] = R0;
                  previous_TS_for_function_neighbour[cell_num_calc(side, column-1, layer+1)] = R0;












            if(dead_time>0.7*pow(10,6)*norm_factor_tracker && dead_time_neighbour>0.7*pow(10,6)*norm_factor_tracker){
                  fill_averaging_histograms(sum_hists[cell_num],counter_hists[cell_num],z,prop_time);
                  full_2d_heatmap_prop_vs_z[cell_num]->Fill(z,prop_time,1);
                  prop_time_over_run[cell_num]->Fill(prop_time,R0,1);
                         }
                  









    }}}

 
average_all_cell_hists(sum_hists,counter_hists,output_hists);

save_all_cell_hists(full_2d_heatmap_prop_vs_z,"prop_vs_z_2d",run_num);  
save_all_cell_hists(prop_time_over_run,"prop_time_over_run",run_num);
save_all_cell_hists(output_hists,"prop_time_vz_z_average",run_num);
save_all_cell_hists(counter_hists,"prop_time_vs_z_counter",run_num);
save_all_cell_hists(dead_time_vs_prop_time,"dead_time_vs_prop_time",run_num);
save_all_cell_hists(dead_time_vs_prop_time_neighbour,"dead_time_vs_prop_time_neighbour",run_num);




  return 0;
}















