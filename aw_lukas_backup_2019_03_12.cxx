/**************************************
Initial code by Peter.Drexler@exp2.physik.uni-giessen.de and Markus.Moritz@exp2.physik.uni-giessen.de
Written by Lukas.Nies@physik.uni-giessen.de
modified for single board / self triggered readout 
**************************************/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <string.h>
#include "aw_lib.h"
#include "tagger_lib.h"
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "Math/Interpolator.h"
#include "Math/Polynomial.h"
#include "TStyle.h"
#include "TList.h"
#include "TExec.h"
#include "TText.h"
#include "TGraphPainter.h"

using namespace std;

// global constants
int TRACELEN = 250;
int TRACELEN_EFF;
int MATRIX_J = 1;
int MATRIX_K = 2;
int VERBOSE=1;
int CHANNELS = 8;
int CHANNELS_EFF;
int SAMPLE_t = 10; // Sample frequency
double SAMPLE_t_eff; // Effective sample frequency
Int_t N_E_WINDOW = 12;
Int_t E_WINDOW_LEFT = 0; 
Int_t E_WINDOW_RIGHT = 6000; 
Int_t E_WINDOW_LENGTH = (int) (E_WINDOW_RIGHT-E_WINDOW_LEFT)/N_E_WINDOW; 
int BASELINE_CUT = 70;
int ENERGY_WINDOW_MAX = 100;  
Int_t N_INTPOL_SAMPLES = 10; // must be even
Int_t NB = 2;
Int_t ENERGY_NORM = 2500; // Cosmis peak will be normed to this channel
// Multi sampling calibration
int MULTIS = 1;
bool MULTIS_CALIB_MODE = false;
// global settings
int THRESHOLD_MULTIPLICY= 1;
bool GLITCH_FILTER = false;
int GLITCH_FILTER_RANGE = 0;
double CFD_fraction = 0.5;
int L=5;     // Length of moving average intervals, centered around current value (mus be odd)
int DELAY=5; // DELAY for the CFD
int M=5; // Window for MWD
int TAU=1; // Impact value for MA part of MWD
// global counters
unsigned int NOE=0;

ReadSystem_class DETECTOR;
// Taggerfile_class tagger;	

struct hist_struct
{
  TF1 *fit;
  TH1D *hist;
  vector<Double_t> fit_params;
};

struct calib_struct
{
  vector<double> multis; // calibration for inter sampling mode
  vector<double> energy; // calibration between different xtals
};

struct multis_norm_struct
{
  hist_struct h_hist;
  Double_t ratio = 0;
  Double_t ratio_err = 0;
};

struct time_struct
{
  vector<hist_struct> h_timing;
  vector<double> timing; // Is N_E_WINDOWS long when built
};

struct baseline_struct
{
  // Container for storing the baseline traces
  vector<Double_t> trace;
  // Variables for storing the baseline statistics
  double mean = 0;
  double std = 0;
  // Store the thresholds 
  double TH = 0;
  hist_struct h_mean;
  hist_struct h_std;
};

struct CFD_struct
{
  vector<Double_t> trace;
  double max = 0;
  int Xzero = 0; 
  int Xzero_int = 0; 
  double int_x0 = 0;
  double int_b = 0;
  double int_m = 0;
};

struct signal_struct
{
  vector<Double_t> trace;
  baseline_struct base;
  CFD_struct CFD;
  double energy = 0; 
  int is_signal = 1;
  hist_struct h_energy;
};

// Construct containers for storing all wave forms + histograms + various informations 
vector<signal_struct> RAW;
vector<signal_struct> RAW_CALIB;
vector<signal_struct> MA;
vector<signal_struct> MWD;
vector<signal_struct> TMAX;
// Build struct to save all time related histograms
vector<vector<time_struct> > MA_TIMES;
vector<vector<time_struct> > MWD_TIMES;
vector<vector<time_struct> > TMAX_TIMES;
// Struct to save all multi sampling renormalisation histograms and parameters
// Dimensions: [eff. channels][MULTIS][MULTIS]
vector<vector<vector<multis_norm_struct> > > MULTIS_NORM;
// Build root file
TFile *hfile;
// Build the calibration struct
calib_struct CALIB;

// Definition of functions
void extraction();
void multis_calib();
void plot_waves(vector<signal_struct> &array, char const *name, char const *modus, double sample_freq);
void plot_interpol(vector<vector<double> > &x, vector<vector<double> > &y);
void plot_time_energy(const vector<vector<Double_t> > &params);
void plot_norm_hists();
double awg_MA(const vector<double> &array, int ch=1, int down=0, int up=0);
double randit(int ini=0);
double array_mean(int start, int end, const vector<double> &array);
double array_std(int start, int end, double mean, const vector<double> &array);
void interpolate(vector<signal_struct> &signal, vector<vector<time_struct> > &times);
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func = "gaus", double range = 200);
void print_usage();
void build_structure();
void print_final_statistics();
void print_stat_multis_calib();
void initialize_histograms(int cahnnels);
void init_intersamp_hist(int channels);
void init_signal(vector<signal_struct> &signal, int channels);
void init_times(vector<vector<time_struct> > &array, int channels);
void reset_signal(vector<signal_struct> &signal);
void reset_times(vector<vector<time_struct> > &array);
void init_multis_norm(vector<vector<vector<multis_norm_struct> > > &array, int channels);
void fill_hists();
void init_hists(int channels);
bool read_config(char const *file);
bool linreg(vector<double> &x, vector<double> &y, double *m, double *b);
Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
Int_t largest_1Dbin(TH1D *hist, Int_t lower, Int_t upper);

  
int main(int argc, char *argv[])
{
  ///////////////////////////////////////////////////////////////
  // Read all existing command line argunments and file names
  ///////////////////////////////////////////////////////////////
  char inputfile[200];
  char outputfile[200];
  char configfile[200];
  unsigned int no_of_events=0, do_no_of_events=0;
  int realign_first_NOE=0;
  // Parse the input file
  if(argc<=1){ 
    printf("No datafile set!!!\n");
    print_usage();
    return(-1);
  }
  strcpy(inputfile, argv[1]);
  // Parse the output file
  if(argc<=2){ 
    printf("No outputfile declared.\n");
    sprintf(outputfile,"%s.root",argv[1]);
    printf("outputfile set to: %s\n",outputfile);
  }
  else strcpy(outputfile, argv[2]);
  // Parse the config file
  bool conf_exists = false;
  if(argc<=3){ 
    printf("No config file set, using standard or command line parameters.\n");
  }
  else {
    conf_exists = true;
    strcpy(configfile, argv[3]);    
    printf("Config file read in: %s\n",configfile);
  }
  // Now parse all command line options
  for(int n=0; n<argc; n++){
    if(strstr(argv[n],"-n")!=NULL){  // set stop after # of counts
      n++;  
      printf("Max count enabled!\n");
      if(n<argc){
        do_no_of_events=1;
        no_of_events=atoi(argv[n]);
      }
      else{
        printf("Missing max. number of events!\n");
        return(-1);
      }
    }
    // Reading of Digital threshold
    if(strstr(argv[n],"-t")!=NULL){  
      n++;
      if(n<argc){
        THRESHOLD_MULTIPLICY=atoi(argv[n]);
        printf("Software threshold multiplicy set to: %i !\n", THRESHOLD_MULTIPLICY);
      }  
      else{
        printf("Missing treshold multiplicy!\n");
        return(-1);
      }
    }
    // Reading of glitch filter settings
    if(strstr(argv[n],"-g")!=NULL){  // set stop after # of counts
      n++;
      if(n<argc){
        GLITCH_FILTER_RANGE=atoi(argv[n]);
        GLITCH_FILTER = true;
        printf("Glitch filter range set to: %i !\n", GLITCH_FILTER_RANGE);
      }  
      else{
        printf("Missing glitch filter parameter!\n");
        return(-1);
      }
    }
    if(strstr(argv[n],"-L")!=NULL){  // set stop after # of counts
      n++;
      if(n<argc){
        L=atoi(argv[n]);
        printf("Moving average range set to: %i !\n", L);
      }  
      // Check if L is odd
      if (L%2 == 0){
        printf("SETTINGS ERROR: MA window must be odd!\n");
        return(-1);
      }

    }
    if(strstr(argv[n],"-d")!=NULL){  // set stop after # of counts
      n++;
      if(n<argc){
        DELAY=atoi(argv[n]);
        printf("DELAY for CFD set to: %i !\n", DELAY);
      }  
    }
    if(strstr(argv[n],"-f")!=NULL){  // set stop after # of counts
      n++;
      if(n<argc){
        CFD_fraction=atof(argv[n]);
        if ( CFD_fraction < 1. && CFD_fraction > 0. ){
          printf("Fraction for CFD set to: %f !\n", CFD_fraction);
        }
        else{
          printf("SETTINGS ERROR: CFD fraction must be set between 0 and 1!\n");
          return(-1);
        }
      }  
    }
    if(strstr(argv[n],"-M")!=NULL){  // set stop after # of counts
      n++;
      if(n<argc){
        M=atoi(argv[n]);
        printf("Window length for MWD set to: %i !\n", M);
      }  
      // Check if M is odd
      if (M%2 == 0){
        printf("SETTINGS ERROR: MWD window must be odd!\n");
        return(-1);
      }
    }
    if(strstr(argv[n],"-I")!=NULL){  // set the number of inter samples
      // Check if MULTIS is 1, 2, or 4
      n++;
      if (atoi(argv[n]) == 1 ||  atoi(argv[n]) == 2 || atoi(argv[n]) == 4){
        if(n<argc){
          MULTIS=atoi(argv[n]);
          printf("Window length for MWD set to: %i !\n", M);
        }  
      }
      else{
        printf("SETTINGS ERROR: Intersampling number not 1,2, or 4!\n");
        return(-1);
      }
    }
    if(strstr(argv[n],"-r")!=NULL){  // realign_first_NOE
      realign_first_NOE=1;
    }
    if(strstr(argv[n],"-m")!=NULL){  // multiple files with startfile
      realign_first_NOE=42;
    }
    if(strstr(argv[n],"-v")!=NULL){  // turn verbose output on
      VERBOSE=1;
    }
    if(strstr(argv[n],"-e")!=NULL){  // turn verbose output on
      MULTIS_CALIB_MODE = true;
      printf("\n\n++++++++++++++++++++++++++++++++\n");
      printf("+++ ENERGY CALIBRATION MODE  +++\n");
      printf("++++++++++++++++++++++++++++++++\n");
    }
    if(strstr(argv[n],"-h")!=NULL){  // print help page
      print_usage();
      return(-1);
    }
  } 
  // Check if MA interval length is larger than sample of baseline cut 
  if (L > BASELINE_CUT){
    printf("SETTINGS ERROR: Moving average interval larger than baseline! Not allowed!\n");
    return(-1);
  }
  // Check if MWD interval length is larger than sample of baseline cut 
  if (M > BASELINE_CUT){
    printf("SETTINGS ERROR: MWD moving average interval larger than baseline! Not allowed!\n");
    return(-1);
  }

  //////////////////////////////////////////////////////////
  // Initialize detector readout
  //////////////////////////////////////////////////////////

  // Read the config file if file is given
  if(conf_exists){ 
    bool conf_healthy = read_config(configfile);
    if (conf_healthy == false){
      printf("ERROR: Config file error!\n");
      return(-1);
    }
  }
  
  // Set up detector read out
  DETECTOR.init_readout(inputfile, realign_first_NOE);
  //TRACELEN=DETECTOR.get_maxevents(1);  

  randit(1);

  //////////////////////////////////////////////////////////
  // Construct TFile for storing all plot/histograms and define folder structure
  //////////////////////////////////////////////////////////

  //TFile hfile(outputfile,"RECREATE","NTEC analysis");
  hfile = new TFile(outputfile,"RECREATE","NTEC analysis");
  hfile->SetCompressionLevel(1);
  // Create folder structure
  build_structure();

  ////////////////////////////////////////////////////////
  // BUILD HISTOGRAMS
  ////////////////////////////////////////////////////////

  // Reset/initialize signal container
  // If in normal extraction mode, use CHANNELS_EFF
  if (MULTIS_CALIB_MODE == false){
    // Initialize the signal containers
    init_signal(RAW, CHANNELS); // Has to be CHANNELS 
    init_signal(RAW_CALIB, CHANNELS_EFF);
    init_signal(MA, CHANNELS_EFF);
    init_signal(MWD, CHANNELS_EFF);
    init_signal(TMAX, CHANNELS_EFF);
    // Initialize the times container
    init_times(MA_TIMES, CHANNELS_EFF);
    init_times(MWD_TIMES, CHANNELS_EFF);
    init_times(TMAX_TIMES, CHANNELS_EFF);
    // Initialize histograms, done in extra function
    init_hists(CHANNELS_EFF);
  }
  else{
    init_signal(RAW, CHANNELS);
    // Initialize histograms, done in extra function
    init_multis_norm(MULTIS_NORM, CHANNELS_EFF);
    init_hists(CHANNELS);
  }

  /////////////////////////////////////////////////////
  // BEGIN SIGNAL READOUT
  /////////////////////////////////////////////////////

  // Initialize the loop condition
  int m=1;
  // Reset NOE number of events counter
  NOE=0;
  do{
    // Keep reading as long there are unread events
    m=DETECTOR.read_one_event(NOE);
    // Print some verbose information
    if (VERBOSE==true){
      if(NOE%1000==0) printf("Analysing event: %i\n",NOE);
    }
    // Increase Global number of events counter
    NOE++;
    // If there is an event, do the extraction
    if(m==1){ // && t!=-13){  // not eof for either of these files
      // MULTISAMPLING CALIBRATION MODE
      // Depending on the program mode, do either calibration or final extraction
      //
      if (MULTIS_CALIB_MODE == true){
        // Extract calibration information
        multis_calib();
        // Fill histograms with result, respecting the correct program mode
        fill_hists();
      }
      // FEATURE EXTRACTION MODE
      // Depending on the program mode, do either calibration or final extraction
      //
      else{
        // Extract energy and timing information
        extraction();
        // Fill histograms with result, respecting the correct program mode
        fill_hists();
        // Every 1000 events plot a random event
        if (NOE%1000==0){
          // Plot every 1000st raw signal
          hfile->cd("WAVE_FORMS/RAW"); 
          plot_waves(RAW, "Signal_RAW", "TRACE", SAMPLE_t);  
        }
      }
    }
    if(do_no_of_events==1){
      if(NOE>=no_of_events) m=0;
    }
  }while(m==1);

  /////////////////////////////////////////////////////
  // END SIGNAL READOUT
  /////////////////////////////////////////////////////

  // Print final statistics
  if (MULTIS_CALIB_MODE==false){
    print_final_statistics();
  }
  else{
    print_stat_multis_calib();
    plot_norm_hists();
  }



  // Write and save and delete the root file element in program
  printf("Writing root file...\n"); 
  hfile->Write();   
  // If put stuff after deleting hfile, software might crash. Beware!
  delete hfile;
  printf("\nRoot file written! Program ends here.\n");
  // If put stuff after deleting hfile, software might crash. Beware!
  
  /////////////////////////////////////////////////////
  // END OF PROGRAM
  /////////////////////////////////////////////////////
}

/////////////////////////////////////////////////////
// ENERGY AND TIMING EXTRACITON
/////////////////////////////////////////////////////
void extraction(){
  /////////////////////////////////////////////////////
  // INITIALIZE SIGNAL CONTAINERS
  /////////////////////////////////////////////////////
  reset_signal(RAW); // raw has to be of lentgh CHANNELS 
  reset_signal(RAW_CALIB);
  reset_signal(MA);
  reset_signal(MWD);
  reset_signal(TMAX);
  //
  reset_times(MA_TIMES);
  reset_times(MWD_TIMES);
  reset_times(TMAX_TIMES);
  // Create dummy container for fetching the signal
  unsigned int dumm_cont[CHANNELS][TRACELEN];
  // Fetch the signal
  for(int board=1; board<=1; board++){
    for(int channel=0; channel<CHANNELS; channel++){ 
      DETECTOR.get_trace(board, channel, dumm_cont[channel]);
    }
  }
  /////////////////////////////////////////////////////
  // FILL THE RAW CONTAINER WITH DUMMY CONTENT
  /////////////////////////////////////////////////////
  for(int i=0; i<CHANNELS; i++){
    for(int n=0; n<TRACELEN; n++){
      RAW[i].trace.push_back((double) dumm_cont[i][n]);
    }
    // Calculate the baseline information for the RAW and RAW_Calib
    RAW[i].base.mean = array_mean(0,BASELINE_CUT,RAW[i].trace);
    // Already sutract the baseline from the signals
    for (int n = 0; n<TRACELEN; n++){
      RAW[i].trace[n] -= RAW[i].base.mean;
    }
    // Recalculate the baseline information for the RAW and RAW_Calib
    RAW[i].base.mean = array_mean(0,BASELINE_CUT,RAW[i].trace);
    // Now calculate the std and TH for the baselines
    RAW[i].base.std = array_std(0, BASELINE_CUT,RAW[i].base.mean,RAW[i].trace);
    RAW[i].base.TH = THRESHOLD_MULTIPLICY * RAW[i].base.std; 
  }
  /////////////////////////////////////////////////////
  // FILL THE CALIBRATED RAW CONTAINER
  /////////////////////////////////////////////////////
  // If multi sampling is used, sort samples correctly 
  // Already apply multi sampling and inter channel energy calibration
  for(int i=0; i<CHANNELS_EFF; i++){
    for(int n=0; n<TRACELEN; n++){
      if (MULTIS == 1){
        RAW_CALIB[i].trace.push_back((double) RAW[i].trace[n] * CALIB.multis[i] * CALIB.energy[i]);
      }
      else if (MULTIS == 2){
        RAW_CALIB[i].trace.push_back((double) RAW[2*i+1].trace[n] * CALIB.multis[2*i+1] * CALIB.energy[i]);       
        RAW_CALIB[i].trace.push_back((double) RAW[2*i].trace[n] * CALIB.multis[2*i] * CALIB.energy[i]);
      }
      else if (MULTIS == 4){
        RAW_CALIB[i].trace.push_back((double) RAW[4*i+3].trace[n] * CALIB.multis[4*i+3] * CALIB.energy[i]);
        RAW_CALIB[i].trace.push_back((double) RAW[4*i+2].trace[n] * CALIB.multis[4*i+2] * CALIB.energy[i]);
        RAW_CALIB[i].trace.push_back((double) RAW[4*i+1].trace[n] * CALIB.multis[4*i+1] * CALIB.energy[i]);
        RAW_CALIB[i].trace.push_back((double) RAW[4*i].trace[n] * CALIB.multis[4*i] * CALIB.energy[i]);
      }
    }
    // Calculate the RAW_CALIB baseline statistics
    RAW_CALIB[i].base.mean = array_mean(0,BASELINE_CUT*MULTIS,RAW_CALIB[i].trace);
    RAW_CALIB[i].base.std = array_std(0, BASELINE_CUT*MULTIS,RAW_CALIB[i].base.mean,RAW_CALIB[i].trace);
    RAW_CALIB[i].base.TH = THRESHOLD_MULTIPLICY * RAW_CALIB[i].base.std; 
  }
  /////////////////////////////////////////////////////
  // SAMPLE SIGNAL
  /////////////////////////////////////////////////////
  // Loop across all channels
	for(int i=0; i<CHANNELS_EFF; i++){
    /////////////////////////////////////////////////////
    // BASELINE SECTION
    /////////////////////////////////////////////////////
		for(int n=0; n<BASELINE_CUT*MULTIS; n++){
      // MA FILTER
      // Apply moving average filter to baseline, take care that baseline does NOT filter also parts of the rising edge of signal!
      if( n - ((L-1)/2) < 1 ){
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace,i, 0, n + (L-1)/2));
      } 
      if( n - ((L-1)/2) >= 1 && n + ((L-1)/2) <= TRACELEN_EFF ){ 
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace,i, n - ((L-1)/2), n + ((L-1)/2)));
      }
      if( n + ((L-1)/2) > TRACELEN_EFF ){ 
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace, i, n - ((L-1)/2), TRACELEN_EFF));
      }
      MA[i].base.trace.push_back(MA[i].trace[n]); 
      // MWD FILTER
      // Apply Moving Window Deconvolution Filter, take care of boundaries
      if (n == 0 || n == 1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n]);
      }
      if (n-M <= 0 && n>1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[0] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, 0, n-1));        
      }
      if (n-M > 0){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[n-M] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, n-M, n-1));
      }
      MWD[i].base.trace.push_back(MWD[i].trace[n]);
      // TMAX Filter 
      // To be implemented
      TMAX[i].trace.push_back(0);
      TMAX[i].base.trace.push_back(0);
		}
    // Calculate statistics and software threshold for signals
    MA[i].base.mean = array_mean(0,BASELINE_CUT*MULTIS-(L-1)/2,MA[i].base.trace);    
    MA[i].base.std = array_std(0,BASELINE_CUT*MULTIS-(L-1)/2,MA[i].base.mean,MA[i].base.trace);
    MWD[i].base.mean = array_mean(M,BASELINE_CUT*MULTIS,MWD[i].base.trace);    
    MWD[i].base.std = array_std(M,BASELINE_CUT*MULTIS,MWD[i].base.mean,MWD[i].base.trace);
    TMAX[i].base.mean = array_mean(0,BASELINE_CUT*MULTIS,TMAX[i].base.trace);    
    TMAX[i].base.std = array_std(0,BASELINE_CUT*MULTIS,TMAX[i].base.mean,TMAX[i].base.trace);
    // Calculate software th based on multiplicy of baseline RMS
    MA[i].base.TH = THRESHOLD_MULTIPLICY * MA[i].base.std;
    MWD[i].base.TH = THRESHOLD_MULTIPLICY * MWD[i].base.std; 
    TMAX[i].base.TH = THRESHOLD_MULTIPLICY * TMAX[i].base.std; 
    // Now substract baseline froms samples
    for (int n = 0; n<BASELINE_CUT*MULTIS; n++ ){
      MA[i].trace[n] -= MA[i].base.mean;
      MWD[i].trace[n] -= MWD[i].base.mean;
      TMAX[i].trace[n] -= TMAX[i].base.mean;
    }
    /////////////////////////////////////////////////////
    // SIGNAL REGION
    /////////////////////////////////////////////////////
    for(int n=BASELINE_CUT*MULTIS; n<ENERGY_WINDOW_MAX*MULTIS; n++){
      // Start  with the RAW extraction
      if(RAW[i].trace[n]>RAW[i].energy){RAW[i].energy=RAW[i].trace[n];}
      // Check if new value is higher than previous max value and above TH
      if(RAW_CALIB[i].trace[n]>RAW_CALIB[i].energy){
        // If glitch filter is not activated or trace is almost at the and of range:
        if (GLITCH_FILTER == false){// || n+GLITCH_FILTER_RANGE >= TRACELEN){
          RAW_CALIB[i].energy=RAW_CALIB[i].trace[n];
        }
        // If glitch filter is activated
        else{
          int is_glitch_array[GLITCH_FILTER_RANGE];
          int is_glitch = 1;
          // Check for the next few samples if all are above threshold
          for(int k=n; k<n+GLITCH_FILTER_RANGE; k++){
            // If k-th sample is greater than TH, then save 0 in array
            // Test with RAW signal because RAW_CALIB is not yet calculated
            if (RAW_CALIB[i].trace[k]>RAW_CALIB[i].base.TH){
              is_glitch_array[k-n] = 1;
            }
            // If not, then set it save 1 in array
            else{
              is_glitch_array[k-n] = 0;
            }
          }
          // Check of samples vary around baseline TH 
          for(int k=0; k<GLITCH_FILTER_RANGE; k++){
            is_glitch *= is_glitch_array[k]; 
          }
          // If not a glitch, save energy
          if (is_glitch==1){
            RAW_CALIB[i].energy = RAW_CALIB[i].trace[n];
          }
        }
      }
      // MA FILTER
      // Continue with calculating the moving average (MA)
      if( n - ((L-1)/2) >= 1 && n + ((L-1)/2) <= TRACELEN_EFF ){ 
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace, i, n - ((L-1)/2), n + ((L-1)/2)));
      }
      if (MA[i].trace[n] > MA[i].energy && MA[i].trace[n]>MA[i].base.TH){
        MA[i].energy=MA[i].trace[n];
      }
      // MWD FILTER
      // Apply Moving Window Deconvolution Filter, take care of boundaries
      if (n == 0 || n == 1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n]);
      }
      if (n-M <= 0 && n>1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[0] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, 0, n-1));        
      }
      if (n-M > 0){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[n-M] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, n-M, n-1));
      }
      // Look for maximum
      if (MWD[i].trace[n] > MWD[i].energy && MWD[i].trace[n]>MWD[i].base.TH){
        MWD[i].energy=MWD[i].trace[n];
      }
      // TMAX FILTER
      // To be implemented
      TMAX[i].trace.push_back(0);
      TMAX[i].base.trace.push_back(0);
    }
    /////////////////////////////////////////////////////
    // REST OF SIGNAL
    /////////////////////////////////////////////////////
    for(int n=ENERGY_WINDOW_MAX*MULTIS; n<TRACELEN_EFF; n++){
      // MA FILTER
      // Continue with calculating the moving average (MA)
      if( n - ((L-1)/2) >= 1 && n + ((L-1)/2) <= TRACELEN_EFF ){ 
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace, i, n - ((L-1)/2), n + ((L-1)/2)));
      } 
       if( n >= TRACELEN_EFF - ((L-1)/2) ){ 
        MA[i].trace.push_back(awg_MA(RAW_CALIB[i].trace, i, n - ((L-1)/2), n));
      }
      // MWD FILTER
      // Apply Moving Window Deconvolution Filter, take care of boundaries
      if (n == 0 || n == 1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n]);
      }
      if (n-M <= 0 && n>1){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[0] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, 0, n-1));        
      }
      if (n-M > 0){
        MWD[i].trace.push_back(RAW_CALIB[i].trace[n] - RAW_CALIB[i].trace[n-M] + (1/TAU)*awg_MA(RAW_CALIB[i].trace, i, n-M, n-1));
      }
      // TMAX Filter
      // To be implemented
      TMAX[i].trace.push_back(0);
      TMAX[i].base.trace.push_back(0);
    }
    /////////////////////////////////////////////////////
    // EXTRACTION OF ENERGY FEATURES 
    /////////////////////////////////////////////////////
    if ( RAW_CALIB[i].energy < RAW_CALIB[i].base.TH) {
      RAW_CALIB[i].energy = 0;
      RAW_CALIB[i].is_signal = 0;
      RAW[i].is_signal = 0;
    }
    if ( MA[i].energy < MA[i].base.TH) {
      MA[i].energy = 0;
      MA[i].is_signal = 0;
    }
    if ( MWD[i].energy < MWD[i].base.TH) {
      MWD[i].energy = 0;
      MWD[i].is_signal = 0;
    }
    if ( TMAX[i].energy < TMAX[i].base.TH) {
      TMAX[i].energy = 0;
      TMAX[i].is_signal = 0;
    }
    /////////////////////////////////////////////////////
    // APPLYING CONSTANT FRACTION DISCRIMINATOR
    /////////////////////////////////////////////////////
    for (int n = 0; n<TRACELEN_EFF; n++){
      // Calculate the CFD signal of the RAW signals 
      if (n-DELAY>0){
        RAW[i].CFD.trace.push_back(RAW[i].trace[n-DELAY] - CFD_fraction * RAW[i].trace[n]);
        RAW_CALIB[i].CFD.trace.push_back(RAW_CALIB[i].trace[n-DELAY] - CFD_fraction * RAW_CALIB[i].trace[n]);
        MA[i].CFD.trace.push_back(MA[i].trace[n-DELAY] - CFD_fraction * MA[i].trace[n]);
        MWD[i].CFD.trace.push_back(MWD[i].trace[n-DELAY] - CFD_fraction * MWD[i].trace[n]);
        TMAX[i].CFD.trace.push_back(TMAX[i].trace[n-DELAY] - CFD_fraction * TMAX[i].trace[n]);
      }
      else{
        RAW[i].CFD.trace.push_back(0);
        RAW_CALIB[i].CFD.trace.push_back(0);
        MA[i].CFD.trace.push_back(0);
        MWD[i].CFD.trace.push_back(0);
        TMAX[i].CFD.trace.push_back(0);
      }
      // Normalize CFD 
      // RAW[i].CFD.trace[n] *= (double)(1000/RAW[i].energy);
      // RAW_CALIB[i].CFD.trace[n] *= (double)(1000/RAW_CALIB[i].energy);
      // MA[i].CFD.trace[n] *= (double)(1000/MA[i].energy);
      // MWD[i].CFD.trace[n] *= (double)(1000/MWD[i].energy);
      // TMAX[i].CFD.trace[n] *= (double)(1000/TMAX[i].energy);

      // Look for the sample, where the CFD signal intersects the x-axis for later interpolation analysis
      // All samples after baseline cut should be negative at first and then turn positive once zero crossing point is passed
      if (n > BASELINE_CUT*MULTIS && n<ENERGY_WINDOW_MAX*MULTIS){
        if (RAW[i].CFD.trace[n-1] < 0 && RAW[i].CFD.trace[n] > 0 && RAW[i].CFD.Xzero==0){RAW[i].CFD.Xzero=n;}
        if (RAW_CALIB[i].CFD.trace[n-1] < 0 && RAW_CALIB[i].CFD.trace[n] > 0 && RAW_CALIB[i].CFD.Xzero==0){RAW_CALIB[i].CFD.Xzero=n;}
        if (MA[i].CFD.trace[n-1] < 0 && MA[i].CFD.trace[n] > 0 && MA[i].CFD.Xzero==0){MA[i].CFD.Xzero=n;}
        if (MWD[i].CFD.trace[n-1] < 0 && MWD[i].CFD.trace[n] > 0 && MWD[i].CFD.Xzero==0){MWD[i].CFD.Xzero=n;}
        if (TMAX[i].CFD.trace[n-1] < 0 && TMAX[i].CFD.trace[n] > 0 && TMAX[i].CFD.Xzero==0){TMAX[i].CFD.Xzero=n;}          
      }
    }
  }
  /////////////////////////////////////////////////////
  // TEST FOR COINCIDENCE
  /////////////////////////////////////////////////////
  // Only take events which have been seen at least in two crystals. 
  // Check coincidence only for RAW_CALIB, all other Filters must also be coincident then
  // Count the number of coincident channels 
  int n_coincidences = 0;
  for(int i=0; i<(int)RAW_CALIB.size(); i++){
    for(int j=i+1; j<(int)RAW_CALIB.size(); j++){
      if (RAW_CALIB[i].is_signal == 1 && RAW_CALIB[j].is_signal == 1){
        n_coincidences++;
      }
    }
  }
  // Throw away non-coinc events
  if (n_coincidences == 0) {
    for(int i=0; i<CHANNELS_EFF; i++){
      RAW_CALIB[i].energy = 0;
      MA[i].energy = 0;
      MWD[i].energy = 0;
      TMAX[i].energy = 0;
    }
  }
  // If coinc, event is kept and every k events, a waveform is saved
  // Also, interpolate around the zero crossing point for extracting timing information
  else{
    // print coincident waveforms
    if(NOE%1000==0){
      hfile->cd("WAVE_FORMS/RAW_CALIB");
      plot_waves(RAW_CALIB, "SIGNAL_RAW_CALIB", "TRACE", SAMPLE_t_eff);
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/MA");
      plot_waves(MA, "SIGNAL_MA", "TRACE", SAMPLE_t_eff);
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/MWD");
      plot_waves(MWD, "SIGNAL_MWD", "TRACE", SAMPLE_t_eff);
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/TMAX");
      plot_waves(TMAX, "SIGNAL_TMAX", "TRACE", SAMPLE_t_eff);
      // Set printing folder to CFD 
      hfile->cd("WAVE_FORMS/CFD/RAW_CALIB");
      plot_waves(RAW_CALIB, "CFD_RAW_CALIB", "CFD", SAMPLE_t_eff);
      hfile->cd("WAVE_FORMS/CFD/MA");
      plot_waves(MA, "CFD_MA", "CFD", SAMPLE_t_eff);
    }
    // start the interpolation by looping over all channels
    interpolate(MA, MA_TIMES);
    interpolate(MWD, MWD_TIMES);
    interpolate(TMAX, TMAX_TIMES);
  }
} 

/////////////////////////////////////////////////////
// MULTI SAMPLING CALIBRATION EXTRACTION
/////////////////////////////////////////////////////
// Extract information from the raw data for the inter sampling calibration
void multis_calib(){
  /////////////////////////////////////////////////////
  // INITIALIZE SIGNAL CONTAINERS AND FETCH EVENTS
  /////////////////////////////////////////////////////
  // only initialize the RAW container
  reset_signal(RAW);
  // Create dummy container for fetching the signal
  unsigned int dumm_cont[CHANNELS][TRACELEN];
  // Initialization of the boards
  for(int board=1; board<=1; board++){
    for(int channel=0; channel<CHANNELS; channel++){ 
      DETECTOR.get_trace(board, channel, dumm_cont[channel]);
    }
  }
  for(int i=0; i<CHANNELS; i++){
    for(int n=0; n<TRACELEN; n++){
      RAW[i].trace.push_back((double) dumm_cont[i][n]);
    }
  }
  /////////////////////////////////////////////////////
  // SAMPLE SIGNAL
  /////////////////////////////////////////////////////
  // Variable declaration for feature extraction
  int is_signal[CHANNELS];
  int is_coinc = 1;
  for(int i=0; i<CHANNELS; i++){
    is_signal[i] = 1;
    /////////////////////////////////////////////////////
    // BASELINE SECTION
    /////////////////////////////////////////////////////
    for(int n=0; n<BASELINE_CUT; n++){
      // Save the raw baselines
      RAW[i].base.trace.push_back(RAW[i].trace[n]);
    }
    // Calculate statistics and software threshold for signals
    RAW[i].base.mean = array_mean(0, BASELINE_CUT,RAW[i].base.trace);
    RAW[i].base.std = array_std(0, BASELINE_CUT,RAW[i].base.mean,RAW[i].base.trace);
    // Calculate software th based on multiplicy of baseline RMS
    RAW[i].base.TH = THRESHOLD_MULTIPLICY * RAW[i].base.std;
    // Subtract the baseline from the data 
    for(int n=0; n<BASELINE_CUT; n++){
      RAW[i].trace[n] -= RAW[i].base.mean;
    }
    /////////////////////////////////////////////////////
    // SIGNAL REGION
    /////////////////////////////////////////////////////
    for(int n=BASELINE_CUT; n<ENERGY_WINDOW_MAX; n++){
      // Subtract the baseline
      RAW[i].trace[n] -= RAW[i].base.mean;
      // Check for energy
      if(RAW[i].trace[n]>RAW[i].energy){RAW[i].energy=RAW[i].trace[n];}
      // Check if new value is higher than previous max value and above TH
      if(RAW[i].trace[n]>RAW[i].energy){
        // If glitch filter is not activated or trace is almost at the and of range:
        if (GLITCH_FILTER == false){// || n+GLITCH_FILTER_RANGE >= TRACELEN){
          RAW[i].energy=RAW[i].trace[n];
        }
        // If glitch filter is activated
        else{
          int is_glitch_array[GLITCH_FILTER_RANGE];
          int is_glitch = 1;
          // Check for the next few samples if all are above threshold
          for(int k=n; k<n+GLITCH_FILTER_RANGE; k++){
            // If k-th sample is greater than TH, then save 0 in array
            // Test with RAW signal because RAW_CALIB is not yet calculated
            if (RAW[i].trace[k]>RAW[i].base.TH){
              is_glitch_array[k-n] = 1;
            }
            // If not, then set it save 1 in array
            else{
              is_glitch_array[k-n] = 0;
            }
          }
          // Check of samples vary around baseline TH 
          for(int k=0; k<GLITCH_FILTER_RANGE; k++){
            is_glitch *= is_glitch_array[k]; 
          }
          // If not a glitch, save energy
          if (is_glitch==1){
            RAW[i].energy = RAW[i].trace[n];
          }
        }
      }
    }
    /////////////////////////////////////////////////////
    // REST OF SIGNAL
    /////////////////////////////////////////////////////
    for(int n=ENERGY_WINDOW_MAX; n<TRACELEN; n++){
      RAW[i].trace[n] -= RAW[i].base.mean;
    }
    /////////////////////////////////////////////////////
    // EXTRACTION OF FEATURES 
    /////////////////////////////////////////////////////
    if (RAW[i].energy < RAW[i].base.TH){ is_signal[i] = 0; }
  }
  /////////////////////////////////////////////////////
  // TEST FOR COINCIDENCE
  /////////////////////////////////////////////////////
  for(int i=0; i<CHANNELS; i++){
    is_coinc *= is_signal[i];
  }
  // Throw away non-coinc events
  if (is_coinc == 0) {
    for(int i=0; i<CHANNELS; i++){
      RAW[i].energy = 0;
    }
  }
  else{
    // print coincident waveforms
    if(NOE%1000==0){
      hfile->cd("WAVE_FORMS/RAW");
      plot_waves(RAW, "SIGNAL_RAW", "TRACE", SAMPLE_t);
    }
    // For Coincident events, also calculate the 
    for(int i=0; i<CHANNELS_EFF; i++){
      for(int j=0; j<MULTIS; j++){
        for(int k=0; k<MULTIS; k++){
          MULTIS_NORM[i][j][k].ratio = RAW[i*MULTIS+j].energy / RAW[i*MULTIS+k].energy;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////
// INTERPOLATE zero crossing for given trace
/////////////////////////////////////////////////////
void interpolate(vector<signal_struct> &signal, vector<vector<time_struct> > &times)
{
  // For time extraction, compare all channels with each other and extract the time information and store 
  //  it in the correct field of the TIMING matrix
  // Outer channel loop: Compare ith channel with jth channel
  int warning_counter = 0;
  for(int i=0; i<(int)signal.size(); i++){
    for (int j=i+1; j<(int)signal.size(); j++){
      // Check if channels saw signal, if not do:
      if ( signal[i].is_signal == 0 || signal[j].is_signal == 0){
        for (Int_t k = 0; k < N_E_WINDOW; k++){
          times[i][j-i-1].timing[k] =  0.0;
        }
        continue; 
      }
      // Reset the Xzero marker for both xtals
      signal[i].CFD.Xzero_int=0; signal[j].CFD.Xzero_int=0; 
      // Check for the left interpolation boundary
      int int_smallest;
      int int_largest;
      if (signal[i].CFD.Xzero<signal[j].CFD.Xzero) {
        int_smallest = signal[i].CFD.Xzero;
        int_largest = signal[j].CFD.Xzero;
      }
      else { 
        int_smallest = signal[j].CFD.Xzero; 
        int_largest = signal[i].CFD.Xzero; 
      }
      // Now set the interpolation interval 
      Int_t int_left = int_smallest - N_INTPOL_SAMPLES/2;
      Int_t int_right = int_largest + N_INTPOL_SAMPLES/2;
      Int_t int_range = int_right - int_left;
      // initialize the new containers for storing the interpolated signal
      Double_t xi[2][int_range]; 
      Double_t yi[2][int_range]; 
      vector<vector<Double_t> > x;
      vector<vector<Double_t> > y;
      // Initilaize the interpolator
      ROOT::Math::Interpolator inter0(int_range, ROOT::Math::Interpolation::kPOLYNOMIAL);
      ROOT::Math::Interpolator inter1(int_range, ROOT::Math::Interpolation::kPOLYNOMIAL);
      // Fill arrays with data
      for ( Int_t k = int_left; k < int_right; k++)
      {
        xi[0][k-int_left]  = (Double_t)k; 
        yi[0][k-int_left]  = (Double_t)signal[i].CFD.trace[k];
        xi[1][k-int_left]  = (Double_t)k; 
        yi[1][k-int_left]  = (Double_t)signal[j].CFD.trace[k];
      }
      // Set the Data
      inter0.SetData(int_range, xi[0], yi[0]);
      inter1.SetData(int_range, xi[1], yi[1]);
      // printf("%i %i\n", int_left, int_right);
      // Be careful with the range switching from one grid to the other
      vector<Double_t> tempx;
      vector<Double_t> tempy1;
      vector<Double_t> tempy2;
      for ( Int_t k = 0; k < (Int_t)(NB * int_range - NB + 1); k++ )
      {
        Double_t x_value = (Double_t) int_left + (Double_t) k/NB;
        tempx.push_back(x_value);
        // printf("%f\n", tempx[k]);
        tempy1.push_back( (Double_t) inter0.Eval(x_value));
        tempy2.push_back( (Double_t) inter1.Eval(x_value));
        // Already look for the zero crossing point
        if (tempy1[k] > 0 && signal[i].CFD.Xzero_int==0){signal[i].CFD.Xzero_int=k;}
        if (tempy2[k] > 0 && signal[j].CFD.Xzero_int==0){signal[j].CFD.Xzero_int=k;}
      }
      x.push_back(tempx); x.push_back(tempx); 
      y.push_back(tempy1); y.push_back(tempy2);
      // Old method:
      // Calculate zero crossing for first trace
      // signal[i].CFD.int_m = (y[0][signal[i].CFD.Xzero_int]-y[0][signal[i].CFD.Xzero_int-1])/(x[0][signal[i].CFD.Xzero_int]-x[0][signal[i].CFD.Xzero_int-1]);
      // signal[i].CFD.int_b = y[0][signal[i].CFD.Xzero_int] - signal[i].CFD.int_m*x[0][signal[i].CFD.Xzero_int];
      // signal[i].CFD.int_x0 = - (signal[i].CFD.int_b/signal[i].CFD.int_m);
      // // Calculate zero crossing for second trace
      // signal[j].CFD.int_m = (y[1][signal[j].CFD.Xzero_int]-y[1][signal[j].CFD.Xzero_int-1])/(x[0][signal[j].CFD.Xzero_int]-x[0][signal[j].CFD.Xzero_int-1]);
      // signal[j].CFD.int_b = y[1][signal[j].CFD.Xzero_int] - signal[j].CFD.int_m*x[0][signal[j].CFD.Xzero_int];
      // signal[j].CFD.int_x0 = - (signal[j].CFD.int_b/signal[j].CFD.int_m);

      // Remove some interpolated samples close to the boundaries where the interpolation fails
      for (int k = 0; k < NB; k++){
        x[0].erase(x[0].begin()); x[0].erase(x[0].end()-1);
        x[1].erase(x[1].begin()); x[1].erase(x[1].end()-1);
        y[0].erase(y[0].begin()); y[0].erase(y[0].end()-1);
        y[1].erase(y[1].begin()); y[1].erase(y[1].end()-1);
      }
      // Use fast linear regression fit to find time value
      if (!linreg(x[0],y[0],&signal[i].CFD.int_m,&signal[i].CFD.int_b)) {
        warning_counter++;
        if ( warning_counter < 50){
          printf("WARNING (interpolate): Linreg: Singular matrix, can't solve problem\n");
        }
        else{
          printf("WARNING (interpolate): Linreg: Further warnings suppressed.\n");
        }
      }
      if (!linreg(x[1],y[1],&signal[j].CFD.int_m,&signal[j].CFD.int_b)) {
        warning_counter++;
        if ( warning_counter < 50){
          printf("WARNING (interpolate): Linreg: Singular matrix, can't solve problem\n");
        }
        else{
          printf("WARNING (interpolate): Linreg: Further warnings suppressed.\n");
        }
      }
      signal[i].CFD.int_x0 = - (signal[i].CFD.int_b/signal[i].CFD.int_m);
      signal[j].CFD.int_x0 = - (signal[j].CFD.int_b/signal[j].CFD.int_m);
      // Print evert 1000th event
      if(NOE%1000==0){
        hfile->cd("WAVE_FORMS/CFD/MA/INTERPOL");
        plot_interpol(x, y);
      }
      // Now extract the timing from the interpolated zero_crossings by "fitting" a linear curve a la f(x)=mx+b
      // Scan through all energy windows and extract timing
      for (Int_t k = 0; k < N_E_WINDOW; k++){
        if (signal[i].energy >= (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH) && 
            signal[i].energy < (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH+E_WINDOW_LENGTH) && 
            signal[j].energy >= (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH) && 
            signal[j].energy < (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH+E_WINDOW_LENGTH)){
          // Calculate and norm time difference
          double time_diff = (signal[i].CFD.int_x0 - signal[j].CFD.int_x0) * SAMPLE_t_eff;
          times[i][j-i-1].timing[k] = time_diff;
          // times[i][j-i-1].timing.push_back( time_diff );
        }
        // else{ times[i][j-i-1].timing.push_back( 0.0 ); }
      }
    }
  }
}

void print_final_statistics(){
  Double_t total_counts=0;
  for (Int_t i = 0; i<CHANNELS_EFF; i++){
    total_counts += RAW[i].h_energy.hist->Integral();
  }
  printf("\n\n\n++++++++     FINAL STATISTICS     ++++++++\n");
  printf("+\n");
  printf("+ TOTAL COUNTS: %i\n", (Int_t) total_counts);
  for (Int_t i = 0; i<CHANNELS_EFF; i++){
    printf("-    Channel %i: %i\n", i, (Int_t) RAW[i].h_energy.hist->Integral());
  }
  printf("+\n");
  Int_t total_coincidents=MA[0].h_energy.hist->Integral();
  Double_t total_efficiency = total_coincidents/(total_counts/CHANNELS_EFF);
  printf("+ TOTAL COINCIDENCES: %i\n", (Int_t)total_coincidents);
  printf("-    Coincidence efficiency %i/%i: %3.1f%%\n", 
        (Int_t)total_coincidents, 
        (Int_t)(total_counts/CHANNELS_EFF), 
        total_efficiency*100);
  printf("-\n");

  if (total_coincidents==0){
    printf("\n\nCoincidence error, check code\n\n");
  }
  // Print Energy Statistics for RAW_CALIB
  for (Int_t i = 0; i < CHANNELS_EFF; i++){
    // fit the energy histograms
    RAW_CALIB[i].h_energy.fit_params = fit_hist(RAW_CALIB[i].h_energy.hist, RAW_CALIB[i].h_energy.fit, "langaus", 200);
    MA[i].h_energy.fit_params = fit_hist(MA[i].h_energy.hist, MA[i].h_energy.fit, "langaus", 200);
    MWD[i].h_energy.fit_params = fit_hist(MWD[i].h_energy.hist, MWD[i].h_energy.fit, "langaus", 200);
  }
  // Print Energy extraction 
  printf("+ COSMIC ENERGY DISTRIBUTIONS FOR MA\n");
  for (int i = 0; i < CHANNELS_EFF; i++){
    printf("-    Channel %i:\n", i);
    printf("-       Pos : %4.1f\n", MA[i].h_energy.fit_params[12]);
    printf("-       FWHM: %4.1f\n", MA[i].h_energy.fit_params[13]);
    printf("-       Calib_param: %4.3f\n", ENERGY_NORM/MA[i].h_energy.fit_params[12]);
    printf("-\n");
  }
  printf("-\n");

  // Start extracting information for all windows
  // Do the fitting for the histograms
  vector<vector<Double_t> > params;
  for (Int_t k = 0; k < N_E_WINDOW; k++){
    params.push_back(fit_hist(MA_TIMES[0][0].h_timing[k].hist, MA_TIMES[0][0].h_timing[k].fit, "gaus", 200));
  }
  hfile->cd("TIMING/MA/ENERGY_WINDOWS");
  plot_time_energy(params);
  // Don't forget to delete the arrays created with new
  Int_t e_matching_efficiency = 0;
  printf("+ TIMING\n");
  for (Int_t i = 0; i < N_E_WINDOW; i++){
    printf("-   TIMING ENERGY WINDOW %i-%i : %3.3f+-%3.3f (%i ENTRIES)\n", 
      (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*i, 
      (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*i +E_WINDOW_LENGTH, 
      params[i][4], 
      params[i][5], 
      (int) MA_TIMES[0][0].h_timing[i].hist->Integral());
    e_matching_efficiency += (Int_t) MA_TIMES[0][0].h_timing[i].hist->Integral();
  }
  printf("-\n");
  printf("-   Timing/Energy matching efficiency: %i/%i: %3.1f%%\n", 
        e_matching_efficiency, 
        total_coincidents, 
        (Double_t)e_matching_efficiency/total_coincidents * 100);
  printf("+\n");
  printf("+++++++++++++++++++++++++++++++++++++++++\n");
} 

void print_stat_multis_calib(){
  Double_t total_counts=0;
  for (int i = 0; i<CHANNELS; i++){
    total_counts += RAW[i].h_energy.hist->Integral();
  }
  printf("\n\n\n++++++++     FINAL STATISTICS     ++++++++\n");
  printf("++++  MULTISAMPLING NORM MODE     ++++++++\n");
  printf("+\n");
  printf("+ TOTAL COUNTS (COINC): %i\n", (Int_t) total_counts);
  for (Int_t i = 0; i<CHANNELS; i++){
    printf("-    Channel %i: %i\n", i, (Int_t) RAW[i].h_energy.hist->Integral());
  }
  printf("+\n");
  printf("+ CALIBRATION RATIOS\n");
  // Do the fitting for the calibration matrices
  for (int i = 0; i<CHANNELS_EFF; i++){
    printf("-   Effective Channel %i\n", i);
    vector<vector<vector<double> > > params;
    for (int j = 0; j<MULTIS; j++){
      vector<vector<double> > dim2;
      for (int k = 0; k<MULTIS; k++){
        dim2.push_back(fit_hist(MULTIS_NORM[i][j][k].h_hist.hist, MULTIS_NORM[i][j][k].h_hist.fit, "gaus", 2));
      }
      params.push_back(dim2);
    }
    // Now print the results
    for (int j = 0; j<MULTIS; j++){
      for (int k = 0; k<MULTIS; k++){
        printf("%f ", params[j][k][2]);
      }
      printf("\n");
    }
    printf("\n");
  }

  printf("+\n");
  printf("+++++++++++++++++++++++++++++++++++++++++\n");
}

void plot_waves(vector<signal_struct> &array, char const *name, char const *modus, double sampling_freq) {
  TCanvas *c_waves = new TCanvas("c1","Wave_forms",200,10,500,300);
  c_waves->SetGrid();
  TMultiGraph *mg_waves = new TMultiGraph();
  mg_waves->SetTitle("Signal example; Time [1ns]; ADC channel [arb. unit]");
  TGraph *tg_waves[array.size()];
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader("ADC Digitization"); // option "C" allows to center the header
  for(Int_t i=0; i<(int)array.size(); i++){   //Channel loop 
    if (strcmp(modus, "TRACE") == 0){
      Double_t wave_y[array[i].trace.size()];
      Double_t wave_x[array[i].trace.size()];
      for(Int_t n = 0; n< (Int_t) array[i].trace.size(); n++){
        wave_y[n] = array[i].trace[n]; 
        wave_x[n] = n * sampling_freq; // Calibrate to the sampling rate
      }
      tg_waves[i] = new TGraph((Int_t) array[i].trace.size(),wave_x,wave_y);
    }
    else if (strcmp(modus, "CFD") == 0){
      Double_t wave_y[array[i].CFD.trace.size()];
      Double_t wave_x[array[i].CFD.trace.size()];
      for(Int_t n = 0; n< (Int_t) array[i].CFD.trace.size(); n++){
        wave_y[n] = array[i].CFD.trace[n];
        wave_x[n] = n * sampling_freq; // Calibrate to the sampling rate
      }
      tg_waves[i] = new TGraph((Int_t) array[i].CFD.trace.size(),wave_x,wave_y);
    } 
    else{ printf("ERROR (plot_waves): Plot option invalid\n");}
    tg_waves[i]->SetLineColor(i+1);
    tg_waves[i]->SetMarkerColor(i+1);
    tg_waves[i]->SetLineWidth(2);
    tg_waves[i]->SetTitle("Channel i");
    //tg_waves[i-1]->Draw("AL*");
    char c_number[12];
    char text[12] = "Channel_";
    sprintf(c_number,"%i",i);
    strcat(text,c_number);
    legend->AddEntry(tg_waves[i],text,"f");
    mg_waves->Add(tg_waves[i]);
  }
  mg_waves->Draw("AL*");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_waves->Write(name);
  delete c_waves;
}

void plot_time_energy(const vector<vector<Double_t> > &params) {
  TCanvas *c_waves = new TCanvas("c1","Wave_forms",200,10,500,300);
  c_waves->SetGrid();
  TGraphErrors *tg_waves;
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader(""); // option "C" allows to center the header
  Double_t wave_x[N_E_WINDOW];
  for(Int_t n = 0; n<N_E_WINDOW; n++){
    wave_x[n] = 5 + n*5;
  }
  Double_t wave_y[N_E_WINDOW];
  for(Int_t n = 0; n<N_E_WINDOW; n++){
    wave_y[n] = params[n][4];
  }
  Double_t wave_yerr[N_E_WINDOW];
  for(Int_t n = 0; n<N_E_WINDOW; n++){
    wave_yerr[n] = params[n][5];
  }
  tg_waves = new TGraphErrors(N_E_WINDOW,wave_x,wave_y,0,wave_yerr);
  tg_waves->GetXaxis()->SetTitle("Energy [MeV]");
  tg_waves->GetYaxis()->SetTitle("Time resolution [ns]");
  tg_waves->SetLineColor(1);
  tg_waves->SetMarkerColor(2);
  tg_waves->SetLineWidth(2);
  tg_waves->SetTitle("Time resolution for different energy ranges");
  

  // tg_waves->SetHighlight();
  // TExec *ex = new TExec("ex", "HighlightHisto()");
  // tg_waves->GetListOfFunctions()->Add(ex);

  // TPad *ph = new TPad("ph", "ph", 0.3, 0.4, 1.0, 1.0);
  // ph->SetFillColor(kBlue-10);
  // ph->Draw();
  // ph->cd();
  // TText *info = new TText(0.5, 0.5, "please move the mouse over the graph");
  // info->SetTextAlign(22);
  // info->Draw();
  // ch->cd();
  

  legend->AddEntry(tg_waves,"Ch1/Ch2","f");
  tg_waves->Draw("AL*");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_waves->Write("Signal_RAW");
  delete c_waves;
}

// void HighlightHisto()
// {
//    TVirtualPad *ph = (TVirtualPad *)gPad->FindObject("ph");
//    if (!ph) return;
//    Int_t ih = g->GetHighlightPoint();
//    if (ih == -1) return;

//    TVirtualPad *savepad = gPad;
//    ph->cd();
//    l->At(ih)->Draw();
//    savepad->cd();
// }

void plot_interpol(vector<vector<double> > &x, vector<vector<double> > &y){
  TCanvas *c_waves = new TCanvas("c1","Wave_interpolation",200,10,500,300);
  c_waves->SetGrid();
  TMultiGraph *mg_waves = new TMultiGraph();
  mg_waves->SetTitle("Interpolation exaple; Time [ns]; interpol. ADC channel [arb. unit]");
  TGraph *tg_waves[CHANNELS_EFF];
  TGraph *tg_fits[CHANNELS_EFF];
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  int warning_counter = 0;
  legend->SetHeader("Interpolated ADC values"); // option "C" allows to center the header
  for(Int_t i=0; i<(Int_t)x.size(); i++){   //Channel loop
    // Fit the passed waves with a linear regression
    double m = 0; double b = 0;
    if (!linreg(x[i],y[i],&m,&b)) {
      warning_counter++;
      if (warning_counter < 50){
        printf("WARNING (plot_interpol): Linreg: Singular matrix, can't solve problem\n");
      }
      else{
        printf("WARNING (plot_interpol): Linreg: Further warnigs supressed.\n");
      }
    }
    // Transform vector to array, for handing it to incompetent root TGraph
    Double_t wave_x[x[i].size()];
    Double_t wave_y[y[i].size()];
    Double_t fit_y[y[i].size()];
    for(Int_t n = 0; n<(Int_t)x[i].size(); n++){
      wave_x[n] = (Double_t) x[i][n] * SAMPLE_t_eff; // Calibrate to the sampling rate
      wave_y[n] = (Double_t) y[i][n];
      fit_y[n] = (Double_t) m * x[i][n] + b; // linear fit for wave
    }
    // Build TGraph
    tg_waves[i] = new TGraph(x[i].size(),wave_x,wave_y);
    tg_waves[i]->SetLineColor(i+1);
    tg_waves[i]->SetMarkerColor(i+1);
    tg_waves[i]->SetLineWidth(2);
    tg_waves[i]->SetTitle("Channel i");
    tg_waves[i]->SetMarkerStyle(kOpenSquare); // Asterisk
    tg_fits[i] = new TGraph(x[i].size(),wave_x,fit_y);
    tg_fits[i]->SetLineColor(i+1);
    tg_fits[i]->SetMarkerColor(i+1);
    tg_fits[i]->SetLineWidth(2);
    tg_fits[i]->SetTitle("Channel i");
    tg_fits[i]->SetMarkerStyle(kDot);
    //tg_waves[i-1]->Draw("AL*");
    char waves_number[25];
    char fit_number[25];
    char waves_text[25] = "Channel_";
    char fit_text[25] = "Channel_";
    sprintf(waves_number,"%i",i);
    sprintf(fit_number,"%i_fit",i);
    strcat(waves_text,waves_number);
    strcat(fit_text,fit_number);
    legend->AddEntry(tg_waves[i],waves_text,"f");
    legend->AddEntry(tg_fits[i],fit_text,"f");
    mg_waves->Add(tg_waves[i]);
    mg_waves->Add(tg_fits[i]);
  }
  mg_waves->Draw("AL");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_waves->Write("Signal_Interpolated"); 
  delete c_waves;
}

void plot_norm_hists(){
  hfile->cd("ENERGY/RAW/Intersampling_calibration");
  // THStack *hs[CHANNELS_EFF];
  for (int i=0; i<CHANNELS_EFF; i++){
    // THStack hs[i] = new THStack("hs","");
    TCanvas *cs = new TCanvas("cs","cs",10,10,700,900);
    int counter = 1;
    cs->Divide(MULTIS,MULTIS);
    for (int j=0; j<MULTIS; j++){
      for (int k=0; k<MULTIS; k++){
        cs->cd(counter); MULTIS_NORM[i][j][k].h_hist.hist->Draw();
        counter++;
      }
    }
    char name[100];
    sprintf(name, "EFF_CHANNEL%i", i);
    cs->Write(name);
    delete cs;
  }
}

double awg_MA(const vector<double> &array, int ch, int down, int up){ 
  double n_samples=0.0;   
  double summe=0.0;
  for(int k=down; k<up; k++){
    n_samples++;
  	summe+=array[k];
  }
  if (n_samples != 0){
    return summe/n_samples; 
  }     			
  else{
    return summe;
  }		
}          

double randit(int ini){
  if(ini==1) srand(time(NULL));
  return (((rand()%100) -50.) /100.);
}

// Calculates the mean of all elements between start and end of an array
double array_mean(int start, int end, const vector<double> &array){
  double sum = 0;
  for (int i = start; i<end; i++){
    sum += array[i];
  }
  return (double) sum/(end-start);
}

// Calculates the std of all elements between start and end of an array
double array_std(int start, int end, double mean, const vector<double> &array){
  double sum = 0;
  for (int i = start; i<end; i++){
    sum += (double) pow(array[i]-mean, 2.);
  }
  return (double) TMath::Sqrt(sum/(end-start));
}

// Function for initializing the global signal contaiers
void init_signal(vector<signal_struct> &signal, int channels){
  // create new items in vector
  for(int i=0; i<channels; i++){
    signal.push_back(signal_struct());
  }
}

void init_times(vector<vector<time_struct> > &array, int channels){
  // For a crystal matrix of MATRIX_J x MATRIX_K xtals, initialize the time structure
  // Fill row j of matrix
  for (int i = 0; i<channels; i++){
    // Fill column k of matrix
    vector<time_struct> tempk;
    for (int j = i+1; j<channels; j++){
      tempk.push_back( time_struct() );
      // Create pointers to N_E_WINDOW histograms in field jk of crystal matrix
      for (int k = 0; k<N_E_WINDOW; k++){
        tempk[j-i-1].h_timing.push_back( hist_struct() );
        tempk[j-i-1].timing.push_back( 0.0 );
      }
    }
    array.push_back(tempk);
  }
}

// Function for initializing the global multisampling normalization contaiers
void init_multis_norm(vector<vector<vector<multis_norm_struct> > > &array, int channels){
  // create new items in vector
  // Fill each channel with 2D vectors
  for(int i=0; i<channels; i++){
    // Create MULTIS columns for MULTIS rows
    vector<vector<multis_norm_struct> > temp2;
    for (int j=0; j<MULTIS; j++){
      // Create MULTIS rows
      vector<multis_norm_struct> temp1;
      for (int k=0; k<MULTIS; k++){
        temp1.push_back(multis_norm_struct());
      }
      temp2.push_back(temp1);
    }
    array.push_back(temp2);
  }
}

///////////////////////////////////////////////////////////////
// INITIALIZE THE HISTOGRAMS
///////////////////////////////////////////////////////////////
void init_hists(int channels){
  // Some variables
  char name[100];
  int BINS=20000;
  int RANGE=20000;
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(1111);
  // If in multi sampling calibration mode, initialize the following
  if (MULTIS_CALIB_MODE == true){
    hfile->cd("ENERGY/RAW/");  
    for (Int_t i=0; i<channels; i++){
      hfile->cd("ENERGY/RAW");
      sprintf(name,"ENERGY_RAW%02d",i);
      RAW[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      RAW[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
    }
    hfile->cd("ENERGY/RAW/Intersampling_calibration");
    for (Int_t i=0; i<CHANNELS_EFF; i++){
      for (Int_t j=0; j<MULTIS; j++){ 
        for (Int_t k=0; k<MULTIS; k++){ 
          sprintf(name,"MULTIS_NORM_HIST_CH%02d_R%i/%i",i,j,k);
          MULTIS_NORM[i][j][k].h_hist.hist=new TH1D(name,"",500,0.0,2.0);
          MULTIS_NORM[i][j][k].h_hist.hist->GetXaxis()->SetTitle("Ratio");
          MULTIS_NORM[i][j][k].h_hist.hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }  
  }
  // Or initilize the histograms for the normal program execution
  else{
    for (Int_t i=0; i<channels; i++){
      hfile->cd("ENERGY/RAW");
      sprintf(name,"ENERGY_RAW%02d",i);
      RAW[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      RAW[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"BASELINE_RAW_MEAN%02d",i);
      RAW[i].base.h_mean.hist=new TH1D(name,"",BINS+20000,0,RANGE+20000);
      RAW[i].base.h_mean.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].base.h_mean.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"BASELINE_RAW_STD%02d",i);
      RAW[i].base.h_std.hist=new TH1D(name,"",BINS-0,0,RANGE);
      RAW[i].base.h_std.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].base.h_std.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/RAW_CALIB");
      // Calibrated raw extraction
      sprintf(name,"ENERGY_RAW_CALIB%02d",i);
      RAW_CALIB[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      RAW_CALIB[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW_CALIB[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      // Enter the MA folder
      hfile->cd("ENERGY/MA");
      sprintf(name,"ENERGY_MA%02d",i);
      gStyle->SetOptFit(1112);
      MA[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      MA[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      MA[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      // Enter the MWD folder
      hfile->cd("ENERGY/MWD");
      sprintf(name,"ENERGY_MWD%02d",i);
      MWD[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      MWD[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      MWD[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      // Enter the TMAX folder
      hfile->cd("ENERGY/TMAX");
      sprintf(name,"ENERGY_TMAX%02d",i);
      TMAX[i].h_energy.hist=new TH1D(name,"",BINS-0,0,RANGE);
      TMAX[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      TMAX[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
    }  
    // Initialize the histograms for the time extraction. 
    // There are Channels_eff*Channels_eff channels to compare
    for (Int_t i=0; i<(int)MA_TIMES.size(); i++){  
      for (Int_t j=i+1; j<(int)MA_TIMES.size(); j++){ 
        // Go to the correct subfolder
        hfile->cd("TIMING/MA/ENERGY_WINDOWS");
        // Now in field ij fill N_E_WINDOW histograms
        for (Int_t k = 0; k<N_E_WINDOW; k++){
          sprintf(name,"MA_TIMING_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          MA_TIMES[i][j-i-1].h_timing[k].hist=new TH1D(name,"",1000,-100,100);
          MA_TIMES[i][j-i-1].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          MA_TIMES[i][j-i-1].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }    
    for (Int_t i=0; i<(int)MWD_TIMES.size(); i++){  
      for (Int_t j=i+1; j<(int)MWD_TIMES.size(); j++){ 
        // Go to the correct subfolder
        hfile->cd("TIMING/MWD/ENERGY_WINDOWS");
        // Now in field ij fill N_E_WINDOW histograms
        for (Int_t k = 0; k<N_E_WINDOW; k++){
          sprintf(name,"MWD_TIMING_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          MWD_TIMES[i][j-i-1].h_timing[k].hist=new TH1D(name,"",1000,-100,100);
          MWD_TIMES[i][j-i-1].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          MWD_TIMES[i][j-i-1].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }
    for (Int_t i=0; i<(int)TMAX_TIMES.size(); i++){  
      for (Int_t j=i+1; j<(int)TMAX_TIMES.size(); j++){ 
        // Go to the correct subfolder
        hfile->cd("TIMING/TMAX/ENERGY_WINDOWS");
        // Now in field ij fill N_E_WINDOW histograms
        for (Int_t k = 0; k<N_E_WINDOW; k++){
          sprintf(name,"TMAX_TIMING_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          TMAX_TIMES[i][j-i-1].h_timing[k].hist=new TH1D(name,"",1000,-100,100);
          TMAX_TIMES[i][j-i-1].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          TMAX_TIMES[i][j-i-1].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }
  }
}

void fill_hists(){
  // Fill histograms respecting the running program mode
  if (MULTIS_CALIB_MODE == true){
    // Fill the calibration histograms
    for (int i = 0; i<CHANNELS; i++){
      if(RAW[i].energy>0) RAW[i].h_energy.hist->Fill(RAW[i].energy); 
    }
    // Fill the calibration histograms
    for (int i = 0; i<CHANNELS_EFF; i++){
      for(int j=0; j< MULTIS; j++){
        for(int k=0; k< MULTIS; k++){
          if (MULTIS_NORM[i][j][k].ratio>0.0){ 
            MULTIS_NORM[i][j][k].h_hist.hist->Fill((Double_t)MULTIS_NORM[i][j][k].ratio, 1);
          }
        }
      }
    }
  }
  else{
    // Fill histograms
    for(int i=0; i<CHANNELS_EFF; i++){
      if(RAW[i].energy>0) RAW[i].h_energy.hist->Fill(RAW[i].energy);
      if(RAW[i].base.mean>0) RAW[i].base.h_mean.hist->Fill(RAW[i].base.mean);
      if(RAW[i].base.std>0) RAW[i].base.h_std.hist->Fill(RAW[i].base.std);
      if(RAW_CALIB[i].energy>0) RAW_CALIB[i].h_energy.hist->Fill(RAW_CALIB[i].energy);
      if(MA[i].energy>0) MA[i].h_energy.hist->Fill(MA[i].energy);
      if(MWD[i].energy>0) MWD[i].h_energy.hist->Fill(MWD[i].energy);
      if(TMAX[i].energy>0) TMAX[i].h_energy.hist->Fill(TMAX[i].energy);
    }
    for(int i=0; i<(int)MA_TIMES.size(); i++){
      for(int j=i+1; j<(int)MA_TIMES.size(); j++){
        for(int k=0; k<N_E_WINDOW; k++){
          // printf("%d %d %d %f\n", i, j, k, MA_TIMES[i][j-i-1].timing[k]);
          if(MA_TIMES[i][j-i-1].timing[k] != 0.0) MA_TIMES[i][j-i-1].h_timing[k].hist->Fill(MA_TIMES[i][j-i-1].timing[k]);
        }
      }
    } 
    for(int i=0; i<(int)MWD_TIMES.size(); i++){
      for(int j=i+1; j<(int)MWD_TIMES.size(); j++){
        for(int k=0; k<N_E_WINDOW; k++){
          // printf("%d %d %d %f\n", i, j, k, MA_TIMES[i][j-i-1].timing[k]);
          if(MWD_TIMES[i][j-i-1].timing[k] != 0.0) MWD_TIMES[i][j-i-1].h_timing[k].hist->Fill(MWD_TIMES[i][j-i-1].timing[k]);
        }
      }
    } 
    for(int i=0; i<(int)TMAX_TIMES.size(); i++){
      for(int j=i+1; j<(int)TMAX_TIMES.size(); j++){
        for(int k=0; k<N_E_WINDOW; k++){
          // printf("%d %d %d %f\n", i, j, k, MA_TIMES[i][j-i-1].timing[k]);
          if(TMAX_TIMES[i][j-i-1].timing[k] != 0.0) TMAX_TIMES[i][j-i-1].h_timing[k].hist->Fill(TMAX_TIMES[i][j-i-1].timing[k]);
        }
      }
    } 
  }
}

// Function for resetting all values in the singal container after one event has been processed
void reset_signal(vector<signal_struct> &signal){
  // create new items in vector
  for(int i=0; i<(int)signal.size(); i++){
    // Reset the energy 
    signal[i].energy = 0;
    // Reset signal variable
    signal[i].is_signal = 1;
    // Reset the baseline statistics arrays
    signal[i].base.mean = 0;
    signal[i].base.std = 0; 
    signal[i].base.TH = 0; 
    // Reset CFD struct
    signal[i].CFD.max = 0;
    signal[i].CFD.Xzero = 0;
    signal[i].CFD.Xzero_int = 0;
    signal[i].CFD.int_x0 = 0;
    signal[i].CFD.int_b = 0; 
    signal[i].CFD.int_m = 0; 
    // Initialize the traces
    signal[i].base.trace.clear(); 
    signal[i].trace.clear(); 
    signal[i].CFD.trace.clear();  
  }
}

void reset_times(vector<vector<time_struct> > &array){
  for (int i = 0; i < (int)array.size(); i++){
    for (int j = i+1; j < (int)array.size(); j++){
      for ( int k = 0; k < (int)array[i][j-i-1].timing.size(); k++){
        array[i][j-i-1].timing[k] = 0.0;
      }
    }
  }
}

// From $ROOTSYS/tutorials/fit/langaus.C
Double_t langaufun(Double_t *x, Double_t *par) {
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  // MP shift correctionj
  mpc = par[1] - mpshift * par[0];
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
  // Convolution integral of Landau and Gaussian by sum
  //   Gaus(Double_t x, Double_t mean, Double_t sigma, Bool_t norm), norm==true is normalizing the gaus
  //   Landau(Double_t x, Double_t mpv, Double_t sigma, Bool_t norm), norm==true is normalizing the landau
  for(i=0.; i<np; i++) {
  // for(i=1.0; i<=np/2; i++) {
     // xx = xlow + (i-.5) * step;
     xx = xlow + (i) * step;
     fland = TMath::Landau(xx,mpc,par[0],true);//+ (par[4]/(pow(x[0], par[5])));
     sum += fland * TMath::Gaus(x[0],xx,par[3], true);

     // xx = xupp - (i-.5) * step;
     // fland = TMath::Landau(xx,mpc,par[0],true);
     // sum += fland * TMath::Gaus(x[0],xx,par[3], true);
  }
  return (par[2] * step * sum  / par[3]);
}

// From $ROOTSYS/tutorials/fit/langaus.C
TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF){
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf
  Int_t i;
  Char_t FunName[100];
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  // TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma", "A", "m");
  // ffit->SetParNames("Width","MP","Area","GSigma");
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  his->Fit(FunName,"RB");   // fit within specified range, use ParLimits, "Q" for quite output
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  return (ffit);              // return fit function
}

// From $ROOTSYS/tutorials/fit/langaus.C
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  // Seaches for the location (x value) at the maximum of the
  // Landau-Gaussian convolute and its full width at half-maximum.
  //
  // The search is probably not very efficient, but it's a first try.
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;
  // Search for maximum
  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    if (l < lold){
      step = -step/10;
    }
    p += step;
  }
  if (i == MAXCALLS)
    return (-1);
  maxx = x;
  fy = l/2;
  // Search for right x location of fy
  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
       step = -step/10;
    p += step;
  }
  if (i == MAXCALLS)
    return (-2);
  fxr = x;
  // Search for left x location of fy
  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
       step = -step/10;
    p += step;
  }
  if (i == MAXCALLS)
    return (-3);
  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

// Look for the bins in a histogram with the most entries between lower and upper
Int_t largest_1Dbin(TH1D *hist, Int_t lower = 0, Int_t upper = 100000){
  // Get number of cells plus overflow and underfolw
  Int_t ncells = hist->GetSize();
  // Check if hist has entries different from under and overflow
  if (ncells < 3){
    printf("WARNING (largest_1Dbin): Histogram empty, check under/overflow!\n");
    return(0);
  } 
  // Check if lower and upper are included in the cells
  if (upper > ncells) {
    upper = ncells -1;
  }
  if (lower > ncells) {
    lower = 1;
  }

  // Now loop over the cells and look for highest entry
  Int_t largest = 0;
  Int_t largest_bin = 0;
  for (Int_t i = lower; i<upper; i++){
    if (largest < hist->GetBinContent(i)){ 
      largest_bin = i; 
      largest = hist->GetBinContent(i); 
    }
  }
  // Return index of largest bin
  return(largest_bin);
}

// Do various fits for a 1D histogram 
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func, double range){
  vector<Double_t> params;
  // Set the fit area around the mean value of the gaussian
  //hist->GetXaxis()->SetRange(hist->GetMean()-range,hist->GetMean()+range);
  // Check for which fit function is chosen
  int n;
  // If gauss
  // Field 0,1: amp, amp_err
  // Field 2,3: mean, mean_err
  // Field 4,5: std, std_err 
  if (strcmp(func, "gaus")==0){n = 3;}
  else if (strcmp(func, "langaus")==0){n = 4;}
  else{ printf("ERROR (fit_hist): Wrong fit function name passed! Currently available: gaus, langaus.\n"); }
  // Check if the histograms were filled
  if ( (int) hist->Integral() < 1 ){
    if (VERBOSE > 0){
      printf("WARNING (fit_hist): Histogram empty!\n");     
    }
    for (int i = 0; i<n; i++){
      params.push_back(0);
      params.push_back(0);
    } 
    return(params);   
  }
  // Now start the fit routines
  // Gaus fit
  if (strcmp(func, "gaus")==0){
    n = 3;
    hist->Fit("gaus", "Q");
    fit = hist->GetFunction("gaus");
    for (Int_t i = 0; i<n; i++){
      params.push_back(fit->GetParameter(i));
      params.push_back(fit->GetParError(i));
    }
    return(params);  
  }
  // Langaus (landau convoluted with gaus) fit
  if (strcmp(func, "langaus")==0){
    n = 6;
    // Rebin a temp histogram before fitting
    Int_t nrebin = 20;
    // Setting fit range and start values
    hist->Rebin(nrebin);
    Double_t largest_bin = largest_1Dbin(hist, 200/nrebin)*nrebin; // heaviest bin between lower and upper bound
    printf("%f\n", largest_bin);
    Double_t fr[2]; // fit boundaries
    Double_t sv[6], pllo[6], plhi[6], fp[6], fpe[6]; 
    fr[0]=0.75*largest_bin; // Lower fit boundary
    fr[1]=3.0*largest_bin; // Upper fit boundary
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //par[4]=A from A/(x^(m))
    //par[5]=m from A/(x^(m))
    //ADDED LATER: par[6]= Maximum of convoluted function
    //ADDED LATER: par[7]= FWHM of convoluted function
    pllo[0]=10.     ; pllo[1]=largest_bin-1000.; pllo[2]=1.0        ; pllo[3]=0.1    ; pllo[4]=-1000.0 ; pllo[5]=0.0001;  // Lower parameter limits
    plhi[0]=1000.   ; plhi[1]=largest_bin+1000.; plhi[2]=10000000000.0; plhi[3]=10000.0; plhi[4]=100000.0; plhi[5]=5.0; // Upper parameter limits
    sv[0]  =500.    ; sv[1]  =largest_bin      ; sv[2]  =5000000.0    ; sv[3]  =1000.0 ; sv[4]  =100.0   ; sv[5]=0.05;// Start values
    Double_t chisqr; // Chi squared
    Int_t ndf; // # degrees of freedom
    fit = langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf); // Fit
    Double_t SNRPeak, SNRFWHM; 
    langaupro(fp,SNRPeak,SNRFWHM); // Search for peak and FWHM
    // Save the parameters in the return container
    for (int i = 0; i<n; i++){ 
      params.push_back(fp[i]);
      params.push_back(fpe[i]);
    }
    // Add the parameters from the langaupro search
    params.push_back(SNRPeak); params.push_back(SNRFWHM);
    return(params);
  }
  else {return(params);}
}

// Fast linear regression fit 
// Source: http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
bool linreg(vector<double> &x, vector<double> &y, double *m, double *b){
  double r; // Correlation coefficient
  int n = (int)x.size();
  double sumx = 0.0;                      /* sum of x     */
  double sumx2 = 0.0;                     /* sum of x**2  */
  double sumxy = 0.0;                     /* sum of x * y */
  double sumy = 0.0;                      /* sum of y     */
  double sumy2 = 0.0;                     /* sum of y**2  */
  for (int i=0; i < n; i++){ 
    sumx  += x[i];       
    sumx2 += pow(x[i], 2.);  
    sumxy += x[i] * y[i];
    sumy  += y[i];      
    sumy2 += pow(y[i], 2.); 
  } 
  double denom = (n * sumx2 - pow(sumx, 2.));
  if (denom == 0) {
    // singular matrix. can't solve the problem.
    *m = 0;
    *b = 0;
    return false;
  }
  *m = (n * sumxy  -  sumx * sumy) / denom;
  *b = (sumy * sumx2  -  sumx * sumxy) / denom;
  r = (sumxy - sumx * sumy / n) /    /* compute correlation coeff */
        sqrt((sumx2 - pow(sumx, 2.)/n) *
        (sumy2 - pow(sumy, 2.)/n));
  return true; 
}


// Prints usage informaion for user
void print_usage(){
  printf("Use the following: aw_lukas datfile.dat rootfile.root (optional: configfile.ini) [-command line options]\n");
  printf("Command line options:\n"); 
  printf("  - n (int)      : Max number of events to be read in.\n"); 
  printf("  - r (int)      : Realign first data file (?).\n"); 
  printf("  - m            : Read in multiple files, starting with first given input file.\n"); 
  printf("  - t (int)      : Set software threshold to t times the base line std.\n"); 
  printf("  - g (int)      : Turns on glitch filter and sets its range.\n"); 
  printf("  - L (odd int)  : Range of moving average interval.\n"); 
  printf("  - d (int)      : Number of sample delay for CFD.\n"); 
  printf("  - v            : Turns on verbose mode. Because more ouput is always better!\n"); 
  printf("  - M (odd int)  : Range of moving window deconvolution window.\n"); 
  printf("  - I (1,2, or 4): Number of inter samples.\n"); 
  printf("  - f (double,0<f<1): Fraction for CFD algorithm.\n"); 
  printf("  - e            : Runs the program in inter sampling calibration mode. Needs option I activated.\n"); 
  printf("  - h            : Shows this help page.\n"); 
  printf("\n\n");
}

// Builds folder structure
void build_structure(){
  hfile->mkdir("ENERGY");
  hfile->mkdir("ENERGY/RAW");
  hfile->mkdir("ENERGY/RAW/Intersampling_calibration");
  hfile->mkdir("ENERGY/RAW_CALIB");
  hfile->mkdir("ENERGY/MA");
  hfile->mkdir("ENERGY/MWD");
  hfile->mkdir("ENERGY/TMAX");
  // Timing hirachy
  hfile->mkdir("TIMING");
  hfile->mkdir("TIMING/MA");
  hfile->mkdir("TIMING/MA/INTERPOL");
  hfile->mkdir("TIMING/MA/ENERGY_WINDOWS");
  hfile->mkdir("TIMING/MWD");
  hfile->mkdir("TIMING/MWD/INTERPOL");
  hfile->mkdir("TIMING/MWD/ENERGY_WINDOWS");
  hfile->mkdir("TIMING/TMAX");
  hfile->mkdir("TIMING/TMAX/INTERPOL");
  hfile->mkdir("TIMING/TMAX/ENERGY_WINDOWS");
  // Wave forms
  hfile->mkdir("WAVE_FORMS");
  hfile->mkdir("WAVE_FORMS/RAW");
  hfile->mkdir("WAVE_FORMS/RAW_CALIB");
  hfile->mkdir("WAVE_FORMS/MA");
  hfile->mkdir("WAVE_FORMS/MWD");
  hfile->mkdir("WAVE_FORMS/TMAX");
  hfile->mkdir("WAVE_FORMS/CFD");
  hfile->mkdir("WAVE_FORMS/CFD/RAW");
  hfile->mkdir("WAVE_FORMS/CFD/RAW_CALIB");
  hfile->mkdir("WAVE_FORMS/CFD/MA");
  hfile->mkdir("WAVE_FORMS/CFD/MA/INTERPOL");
  hfile->mkdir("WAVE_FORMS/CFD/MWD");
  hfile->mkdir("WAVE_FORMS/CFD/MWD/INTERPOL");
  hfile->mkdir("WAVE_FORMS/CFD/TMAX");
  hfile->mkdir("WAVE_FORMS/CFD/TMAX/INTERPOL");
}

bool read_config(const char *file){
  // set status for return function
  ifstream file_in(file);
  // Check if file exists
  if (!file_in.is_open()){
    printf("ERROR (read_config): No config file found under this name!\n");
    return(false);
  }
  // If exists, loop over all lines
  string line;
  while(getline(file_in, line) )
  {
    istringstream is_line(line);
    string key;
    if( getline(is_line, key, '=') )
    {
      string value;
      if( getline(is_line, value, '#') )
      {
        //Now test for all possible parameter cases 
        if(strcmp(key.c_str(), "CHANNELS") == 0) CHANNELS = stoi(value);
        if(strcmp(key.c_str(), "TRACELEN") == 0) TRACELEN = stoi(value);
        if(strcmp(key.c_str(), "SAMPLE_t") == 0) SAMPLE_t = stoi(value);
        
        if(strcmp(key.c_str(), "VERBOSE") == 0) VERBOSE = stoi(value);

        if(strcmp(key.c_str(), "N_E_WINDOW") == 0) N_E_WINDOW = stoi(value);
        if(strcmp(key.c_str(), "E_WINDOW_LEFT") == 0) E_WINDOW_LEFT = stoi(value);
        if(strcmp(key.c_str(), "E_WINDOW_RIGHT") == 0) E_WINDOW_RIGHT = stoi(value);
        if(strcmp(key.c_str(), "BASELINE_CUT") == 0) BASELINE_CUT = stoi(value);
        if(strcmp(key.c_str(), "ENERGY_WINDOW_MAX") == 0) ENERGY_WINDOW_MAX = stoi(value);
        if(strcmp(key.c_str(), "ENERGY_NORM") == 0) ENERGY_NORM = stoi(value);
        if(strcmp(key.c_str(), "N_INTPOL_SAMPLES") == 0) N_INTPOL_SAMPLES = stoi(value);
        if(strcmp(key.c_str(), "NB") == 0) NB = stoi(value);
        if(strcmp(key.c_str(), "MULTIS") == 0) MULTIS = stoi(value);
        if(strcmp(key.c_str(), "THRESHOLD_MULTIPLICY") == 0) THRESHOLD_MULTIPLICY = stoi(value);
        if(strcmp(key.c_str(), "GLITCH_FILTER") == 0) GLITCH_FILTER = stoi(value);
        if(strcmp(key.c_str(), "GLITCH_FILTER_RANGE") == 0) GLITCH_FILTER_RANGE = stoi(value);
        // Feature extraction algorithm parameters
        if(strcmp(key.c_str(), "L") == 0) L = stoi(value);
        if(strcmp(key.c_str(), "DELAY") == 0) DELAY = stoi(value);
        if(strcmp(key.c_str(), "M") == 0) M = stoi(value); 
        if(strcmp(key.c_str(), "TAU") == 0) TAU = stoi(value);
        if(strcmp(key.c_str(), "CFD_fraction") == 0) CFD_fraction = stod(value);
        // Calibration parameters
        // Check for the MULTIS_CALIB's
        string subkey = key.substr(0, (int)key.size()-2);
        if(strcmp(subkey.c_str(), "MULTIS_CALIB") == 0){
          // get the calib channel
          CALIB.multis.push_back( stod(value) );
        } 
        if(strcmp(subkey.c_str(), "ENERGY_CALIB") == 0){
          // get the calib channel
          CALIB.energy.push_back( stod(value) );
        } 
      }
    }
  }
  // Calculate effective channel number, effective trace length, and effective sample frequency
  TRACELEN_EFF = (int) TRACELEN*MULTIS;
  CHANNELS_EFF = (int) CHANNELS/MULTIS;
  SAMPLE_t_eff = (double) SAMPLE_t/MULTIS;
  if (VERBOSE > 0){
    printf("\n");
    printf("+++ CONFIG FILE PARAMETERS SET ++\n");
    printf("Number of channels: %i\n", CHANNELS);
    printf("Lentgh of moving average interval: %i\n", L);
    printf("Lentgh of moving window interval: %i\n", M);
    printf("Fraction of CFD algorithm: %3.2f\n", CFD_fraction);
    printf("Length of trace: %i\n", TRACELEN);
    printf("Effective length of trace: %i\n", TRACELEN_EFF);
    printf("Effective number of channels: %i\n", CHANNELS_EFF);
    printf("+++++++++++++++++++++++++++\n\n");
  }
  // Some tests for the read configuration parameters
  // Test if the dimension of the calibration parameters fit the number of channels, etc.
  if ( (int)CALIB.energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every channel has its energy calibration factor!\n");
    return(false);
  }
  if ( (int)CALIB.multis.size() < CHANNELS ){
    printf("ERROR (read_config): Not every channel has its energy multi sampling factor!\n");
    return(false);
  }
  if ( MULTIS > CHANNELS ){
    printf("ERROR (read_config): Not enough channels for given value of MULTIS!\n");
    return(false);
  }
  if ( N_INTPOL_SAMPLES % 2 != 0 ){
    printf("ERROR (read_config): N_INTPOL_SAMPLES must be even!\n");
    return(false);
  }


  return(true);
}

