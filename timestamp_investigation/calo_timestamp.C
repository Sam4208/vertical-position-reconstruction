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
#include "/sps/nemo/scratch/spratt/sndisplay/sndisplay.cc"
using namespace std;

void calo_hit_relative_times_function(TH2D *hist,double calo_nohit , std::vector<int> *calo_ts){
  
}




int cell_num_calc(double side, double column, double layer){

    int cell_num =  side*9*113 + column*9 + layer;
    return cell_num;

}

string name_folder(string run_num){
  
    string directory = "/sps/nemo/scratch/spratt/analysis/timestamp_investigation/run_output_analysis_files/";
    std::string folder = directory+run_num;

    return folder;
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




TH2D** make_array_of_hist(string name_of_hist, double bins,double lower_bound,double upper_bound,double bins2,double lower_bound2,double upper_bound2){
    
      int no_cells=2034;
      TH2D** h = new TH2D*[no_cells];
      for (int i=0;i<no_cells;i++){
            string cell_num = to_string(i);
            string name = (name_of_hist+cell_num);
            h[i] = new TH2D(name.c_str(),name.c_str(), bins, lower_bound, upper_bound ,bins2, lower_bound2, upper_bound2);
                        }

      return h;
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
  TFile *rootfile = new TFile((folder+"/"+name+".root").c_str(), "RECREATE");//saving all histograms
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





int calo_track_corresponder(int calo_column, int track_layer){
  if (calo_column == 0 && track_layer >= 0 && track_layer <= 6) return 1;
  if (calo_column == 1 && track_layer >= 2 && track_layer <= 11) return 1;
  if (calo_column == 2 && track_layer >= 8 && track_layer <= 17) return 1;
  if (calo_column == 3 && track_layer >= 14 && track_layer <= 24) return 1;
  if (calo_column == 4 && track_layer >= 19 && track_layer <= 29) return 1;
  if (calo_column == 5 && track_layer >= 25 && track_layer <= 34) return 1;
  if (calo_column == 6 && track_layer >= 31 && track_layer <= 40) return 1;
  if (calo_column == 7 && track_layer >= 37 && track_layer <= 46) return 1;
  if (calo_column == 8 && track_layer >= 43 && track_layer <= 52) return 1;
  if (calo_column == 9 && track_layer >= 49 && track_layer <= 58) return 1;
  if (calo_column == 10 && track_layer >= 55 && track_layer <= 64) return 1;
  if (calo_column == 11 && track_layer >= 61 && track_layer <= 70) return 1;
  if (calo_column == 12 && track_layer >= 66 && track_layer <= 75) return 1;
  if (calo_column == 13 && track_layer >= 72 && track_layer <= 81) return 1;
  if (calo_column == 14 && track_layer >= 78 && track_layer <= 87) return 1;
  if (calo_column == 15 && track_layer >= 84 && track_layer <= 93) return 1;
  if (calo_column == 16 && track_layer >= 90 && track_layer <= 99) return 1;
  if (calo_column == 17 && track_layer >= 96 && track_layer <= 105) return 1;
  if (calo_column == 18 && track_layer >= 101 && track_layer <= 110) return 1;
  if (calo_column == 19 && track_layer >= 107 && track_layer <= 112) return 1;
  else return 0;
}













int main() {

string run_num = "1051";
//std::string file_1 = "/sps/nemo/scratch/golivier/aussois_tuto_data/bi207_run/snemo_run-840_udd.root";//Chosen run file 
std::string file_1 = "/sps/nemo/snemo/snemo_data/reco_data/UDD/delta-tdc-10us-v3/snemo_run-1051_udd.root";

make_dir(run_num);

TFile *file= new TFile(file_1.c_str(), "READ");

  std::vector<std::vector<short>> *wave = new std::vector<std::vector<short>>;//initiating vectors that will be filled with data from the event file
  std::vector<int> *calo_wall   = new std::vector<int>;
  std::vector<int> *calo_side   = new std::vector<int>;
  std::vector<int> *calo_column = new std::vector<int>;
  std::vector<int> *calo_row    = new std::vector<int>;
  std::vector<int> *calo_charge = new std::vector<int>;

  std::vector<int> *calo_ampl   = new std::vector<int>;
  std::vector<int> *tracker_side   = new std::vector<int>;
  std::vector<int> *tracker_column = new std::vector<int>;
  std::vector<int> *tracker_row    = new std::vector<int>;
  std::vector<int> *tracker_layer  = new std::vector<int>;
  std::vector<long> *calo_ts     = new std::vector<long>;
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
  tree->SetBranchStatus("digicalo.nohits",1);
  tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
  tree->SetBranchStatus("digicalo.waveform",1);
  tree->SetBranchAddress("digicalo.waveform", &wave);
  tree->SetBranchStatus("digicalo.wall",1);
  tree->SetBranchAddress("digicalo.wall", &calo_wall);
  tree->SetBranchStatus("digicalo.side",1);
  tree->SetBranchAddress("digicalo.side", &calo_side);
  tree->SetBranchStatus("digicalo.column",1);
  tree->SetBranchAddress("digicalo.column", &calo_column);
  tree->SetBranchStatus("digicalo.row",1);
  tree->SetBranchAddress("digicalo.row", &calo_row);
  tree->SetBranchStatus("digicalo.charge",1);
  tree->SetBranchAddress("digicalo.charge", &calo_charge);
  tree->SetBranchStatus("digicalo.peakamplitude",1);
  tree->SetBranchAddress("digicalo.peakamplitude", &calo_ampl);
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

   bool plot_events=false;

   double get_entries = (tree->GetEntries() - 1);
   double norm_factor_calo = 6.25;
   double norm_factor_tracker =12.5;

   int no_events_loop = 100;
   TH2D** time_stamps_plot = make_array_of_hist("time_stamps_plot",2000,-20000,100000,15,0,15);
   TH2D** drift_time_vs_prop_time_per_cell = make_array_of_hist("drift_time_vs_prop_time_per_cell",1000,-1000,5000,1000,30000,80000);

   TH2D* drift_time_vs_prop_time = new TH2D("drift_time_vs_prop_time","drift_time_vs_prop_time",1000,-1000,5000,1000,30000,80000);

   TH1D* R0_vs_RC = new TH1D("R0_vs_RC","R0_vs_RC",1000,-1000,5000);

   TH2D** full_2d_heatmap_prop_vs_z = make_array_of_hist("full_2d_heatmap_prop_vs_z",100,-1,1,1000,30000,80000);

   TH2D** drift_time_vs_z_array = make_array_of_hist("drift_time_vs_z_array",100,-1,1,1000,-1000,5000);


   double R0;
   double bot_time;
   double top_time;
   double R1;
   double R2;
   double R3;
   double R4;
   double R5;
   double R6;
   double RC;
   double RC_ground;
   double RC_zero;
   double t_drift;
   double t_drift_0;
   double t_drift_1;

   double side;
   double column;
   double layer;

   int cell_num;



   double prop_time;


string folder = name_folder(run_num);

make_dir((run_num+"/event_displays").c_str());

sndisplay::demonstrator *sndemonstrator = new sndisplay::demonstrator ("demonstrator_test");
TCanvas* full_canvas = new TCanvas("full_canvas");



int no_events_to_run = round((get_entries));

  for (int i = 0; i < no_events_to_run; i++) {  //looping over events
    tree->GetEntry(i); 
    






    for (int k = 0; k < calo_nohits; k++) {      //k : loop on calo hit number
      int om_num = calo_side->at(k)*260 + calo_column->at(k)*13 + calo_row->at(k);

      if (-calo_ampl->at(k) > 200 && om_num < 520) {      // condition to cut small charge and keep only MW OM

        int layer_tab[9]{0};
        int layer_sum = 0;
        int last_layer = 0;


        for (int j = 0; j < tracker_nohits; j++) {        //j  loop on tracker hit number

          if (layer_tab[tracker_layer->at(j)] == 0 && tracker_side->at(j) == calo_side->at(k)) {
            layer_tab[tracker_layer->at(j)] = 1;             // loop to count the number of layer with a trace on the same side than the calo hit
            if (tracker_layer->at(j) == 8) {
              last_layer = j;                               // to compare the layer with calo column
            }
          }
        }

        for (size_t layer_num = 0; layer_num < 9; layer_num++) {
          layer_sum += layer_tab[layer_num];
        }

        if(layer_sum > 7){                                  // at least 7 layer between the source and the calorimeter
          if (layer_tab[0] == 1 || layer_tab[1] == 1) {                       //at least a track hit near the source foil

            if(calo_track_corresponder(calo_column->at(k), tracker_column->at(last_layer)) == 1){     // condition between calo column and tracker layer
             



            
            sndemonstrator->setomcontent(om_num,  1.0); // M:0.1.10 => 260*0 + 13*1 + 10 = 23
    
      

           


             for (int j = 0; j < tracker_nohits; j++) {  

            RC = calo_ts->at(k)*norm_factor_calo;
            RC_zero = calo_ts->at(0)*norm_factor_calo;
            RC_ground = calo_ts->at(0)*norm_factor_calo;

            side = tracker_side->at(j);
            column = tracker_column->at(j);
            layer = tracker_layer->at(j);
            
            cell_num = cell_num_calc(side, column, layer);


              R0 = tracker_R0->at(j).at(0)*norm_factor_tracker;// assigning new labels for anode and cathode timestamps and 'zeroing' them wrt R0 - should be changed to the calo hit that inititated gathering the data 
              bot_time = tracker_botcathtime->at(j).at(0)*norm_factor_tracker;
              top_time = tracker_topcathtime->at(j).at(0)*norm_factor_tracker;
              R1 = tracker_R1->at(j).at(0)*norm_factor_tracker;
              R2 = tracker_R2->at(j).at(0)*norm_factor_tracker;
              R3 = tracker_R3->at(j).at(0)*norm_factor_tracker;
              R4 = tracker_R4->at(j).at(0)*norm_factor_tracker;
              R5 = bot_time;
              R6 = top_time;
              prop_time = R5-R0+R6-R0;
              

       
            if(R5-R0>0 && R6-R0>0 && R5-R0<1000000 && R6-R0<1000000){
              double z=(R5-R6)/(R5+R6-(2*R0));
            //  cout<<R0-RC<<endl;
            //  cout<<"proptime"<<R5-R0+R6-R0<<endl;


           if(i<no_events_loop){
              time_stamps_plot[i]->Fill(R0-RC_ground,0);
              time_stamps_plot[i]->Fill(R1-RC_ground,1);
              time_stamps_plot[i]->Fill(R2-RC_ground,2);
              time_stamps_plot[i]->Fill(R3-RC_ground,3);
              time_stamps_plot[i]->Fill(R4-RC_ground,4);
              time_stamps_plot[i]->Fill(R5-RC_ground,5);
              time_stamps_plot[i]->Fill(R6-RC_ground,6);
              time_stamps_plot[i]->Fill(RC-RC_ground,8+k);//if k is 0 goes to 8n if k is 3 goes to 11 ... 
                                                          }



              sndemonstrator->setggcontent(tracker_side->at(j), tracker_column->at(j), tracker_layer->at(j), 1);
              


              drift_time_vs_prop_time->Fill(R0-RC,prop_time,1);
              if (z<0.1 && z>-0.1){
              drift_time_vs_prop_time_per_cell[cell_num]->Fill(R0-RC,prop_time,1);}
              R0_vs_RC->Fill(R0-RC);
              drift_time_vs_z_array[cell_num]->Fill(z,R0-RC,1);
         //   cout<<"hello"<<endl;
            if(R0-RC<200){
              full_2d_heatmap_prop_vs_z[cell_num]->Fill(z,prop_time,1);}
             


            
            }
            
           
            }
            if (plot_events==true){
            if (i<100){
            sndemonstrator->settitle(("Full detector view event "+to_string(eventnumber)+"calo hit "+to_string(k)).c_str());
            sndemonstrator->draw_top();
            sndemonstrator->canvas->SaveAs((folder+"/event_displays/display_event_"+to_string(eventnumber)+"_calo_hit_"+to_string(k)+".png").c_str());
            }}
            sndemonstrator->reset();



            


            }
          }
        }
        layer_sum = 0;
      }
    }

           // full_canvas->cd();
           // time_stamps_plot[i]->Draw();
           // full_canvas->Update();
           // full_canvas->SaveAs((folder+"/event_order_event_"+to_string(eventnumber)+".png").c_str());
           // full_canvas->Clear();

  }







cout<<"hello"<<endl;
save_all_cell_hists(time_stamps_plot,"event_displays/time_stamps_by_event",run_num);
save_all_cell_hists(drift_time_vs_prop_time_per_cell,"individual_cell_data/drift_time_vs_prop_time_per_cell",run_num);
save_all_cell_hists(full_2d_heatmap_prop_vs_z,"individual_cell_data/prop_vs_z_2d_with_drift_cut",run_num);  
save_all_cell_hists(drift_time_vs_z_array,"individual_cell_data/drift_time_vs_z",run_num);  

string title_root_file = "/calo_info.root";
TFile *rootfile = new TFile((folder+title_root_file).c_str(), "RECREATE");//saving all histograms
rootfile->cd();

drift_time_vs_prop_time->Write();
R0_vs_RC->Write();

rootfile->Close();


  return 0;
}















