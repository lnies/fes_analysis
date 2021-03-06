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
#include <cstddef>
#include <chrono> // measuring high-res execution time 
#include "Math/Interpolator.h"
#include "Math/Polynomial.h"
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
#include "TStyle.h"
#include "TList.h"
#include "TExec.h"
#include "TText.h"
#include "TGraphPainter.h"
#include "TSpectrum.h"  // For spectrum and peak analysis
#include "TVirtualFitter.h" // fitting
#include "TPaveStats.h"
// RooFit Framework
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "TLatex.h"
#include "RooNovosibirsk.h"

using namespace RooFit ;

using namespace std;
using namespace chrono;
using namespace RooFit;


int doof = 0;


// global constants
int TRACELEN = 250;
int MATRIX_J = 1;
int MATRIX_K = 2;
char VERBOSE[200];
int CHANNELS = 8;
int TAG_CHANNELS = 16;
int TAGGER_WINDOW = 100;
int CENTRAL = 4;
int BOARDS = 1;
int CHANNELS_EFF;
int SAMPLE_t = 10; // Sample frequency in ns
char MODE[100]; 
double SAMPLE_t_eff; // Effective sample frequency
Int_t N_E_WINDOW = 12;
Int_t E_WINDOW_LEFT = 0; 
Int_t E_WINDOW_RIGHT = 6000; 
Int_t E_WINDOW_LENGTH = (int) (E_WINDOW_RIGHT-E_WINDOW_LEFT)/N_E_WINDOW; 
int COINC_LEVEL = 1; // number of coincidences needed per event
int BASELINE_CUT = 70;
int ENERGY_WINDOW_MAX = 100;  
Int_t N_INTPOL_SAMPLES = 10; // must be even
Int_t NB = 2;
Int_t ENERGY_NORM = 2500; // Cosmis peak will be normed to this channel
int NB_ACT_CHANNELS = 0; // Number of active channels
double GENERAL_SCALING = 1.; // General histogram scaling for energy calibration
int EXTRACT_PROTO=1; // Extract proto trace and save it to file
double LOWER_RATIO = 0.;
double UPPER_RATIO = 100000.;
// Multi sampling calibration
int MULTIS = 1;
bool MULTIS_CALIB_MODE = false;
// global settings
double THRESHOLD_MULTIPLICY= 1;
bool GLITCH_FILTER = false;
bool SATURATION_FILTER = false;
int GLITCH_FILTER_RANGE = 0;
double CFD_fraction = 0.5;
int L=5;     // Length of moving average intervals, centered around current value (mus be odd)
int DELAY=5; // DELAY for the CFD
int M=5; // Window for MWD
double TAU=1.; // Impact value for MA part of MWD
int nTabs = 20; // Number of tabs for the FIR filter
vector<double> FIR_COEF; // Array for the FIR filter coefficients
double LIN_COMP = 0.; // Fit Parameter for linearizing the detector energy sum
// global counters
unsigned int NOE=0;

ReadSystem_class DETECTOR;
// Taggerfile_class tagger;	

struct mapping_struct
{
  int multis = 1; // multiplication for the multisampling (multisampling)
  int polarity = +1; // Polarity of the input signal
  // if multisampling is true, save start and end h_channel
  int h_channel_start = 0;
  int h_channel_stop = 0;
  int s_channel = 0; // software channel address
  int board_nb = 0; // board number
};

struct hist_struct_TH1D
{
  TF1 *fit;
  TH1D *hist;
  vector<Double_t> params;
};

struct hist_struct_TH2D
{
  TF1 *fit;
  TH2D *hist;
  vector<Double_t> params;
};

struct calib_struct
{
  vector<double> multis; // calibration for inter sampling mode
  vector<double> RAW_energy; // calibration between different xtals
  vector<double> MA_energy; // calibration between different xtals
  vector<double> MWD_energy; // calibration between different xtals
  vector<double> TMAX_energy; // calibration between different xtals
  vector<double> NMO_energy; // calibration between different xtals
  hist_struct_TH1D h_RAW_energy; // hist the calibration parameters
  hist_struct_TH1D h_MA_energy;
  hist_struct_TH1D h_MWD_energy;
  hist_struct_TH1D h_TMAX_energy;
  hist_struct_TH1D h_NMO_energy;
};

struct multis_norm_struct
{
  hist_struct_TH1D h_hist;
  Double_t ratio = 0.;
  Double_t ratio_err = 0.;
};

struct time_struct
{
  vector<hist_struct_TH1D> h_timing;
  vector<double> timing; // Is N_E_WINDOWS long when built
};

struct tagger_energy
{
  hist_struct_TH1D h_energy; // General energy tagged histogram
  hist_struct_TH1D h_energy_m; // Histogram w/o multiples
  hist_struct_TH1D h_energy_mt; // Histogram w/o multiples and w/ timing cut
  double energy = 0.0; 
  double energy_m = 0.0; 
  double energy_mt = 0.0; 
  // 
  hist_struct_TH1D h_integral; // General integral tagged histogram
  hist_struct_TH1D h_integral_m; // Histogram w/o multiples
  hist_struct_TH1D h_integral_mt; // Histogram w/o multiples and w/ timing cut
  double integral = 0.0; 
  double integral_m = 0.0; 
  double integral_mt = 0.0; 

};

struct baseline_struct
{
  // Container for storing the baseline traces
  vector<Double_t> trace;
  // Variables for storing the baseline statistics
  double mean = 0.;
  double std = 0.;
  // Store the thresholds 
  double TH = 0.;
  double TH_multiplicity = 0.;
  hist_struct_TH1D h_mean;
  hist_struct_TH1D h_std;
  hist_struct_TH1D h_samples;
};

struct CFD_struct
{
  vector<Double_t> trace; // CFD trace
  vector<double> x_interpol; //x- interpolated section of the signal zero crossing
  vector<double> y_interpol; //y- interpolated section of the signal zero crossing
  double max = 0.;
  int Xzero = 0; 
  int Xzero_int = 0; 
  double int_x0 = 0.; // the important time information
  double int_b = 0.;
  double int_m = 0.;
};

struct signal_struct
{
  // Trace is saved 
  vector<Double_t> trace;
  // weighted sum of all traces, not to be reset
  hist_struct_TH2D TH2D_proto_trace;
  vector<Double_t> proto_trace;
  vector<Double_t> proto_trace_fit;
  vector<Double_t> proto_trace_maxbin;
  vector<hist_struct_TH1D> TH1D_proto_trace;
  // Baseline is copied and saved separatly here
  baseline_struct base;
  // CFD conversion of trace
  CFD_struct CFD;
  // Energy of the signal by pulsehight
  double energy = 0.; 
  int energy_n = 0; // sample number of peak 
  // Energy of the signal
  double integral = 0.; 
  double ratio = 0.; // ratio of integral/pulseheight
  // Marker for signal or non-signal
  int is_signal = 0;
  // General energy histogram
  hist_struct_TH1D h_energy;
  // General integral histogram
  hist_struct_TH1D h_integral;
  hist_struct_TH1D h_ratio; // ratio of integral ober pulseheight
  // Tagged energy histograms and energies
  vector<tagger_energy> tagged;
  // Timing content
  vector<time_struct> time; 
  // Some information about the harware characteritics of the channel
  bool is_valid = false; // if this channel is not used for feature extraction, then config file will set this to false
  int hardware_addr = 0;
  int clock_speed = 100; // in MHz
  double sample_t = 10.; // Sampling time in ns
  int multis = 1; // Multisampling, number of hardware channels per software channel
  int tracelen = 250; // number of samples in trace
  int polarity = 1; // Polarity of input signal
  bool is_raw = false; // Flag for the RAW container
};

struct tagger_struct
{
  double time[16]; // Time information of tagger (arbitrary information)
  int counts[16]; // Counts per tagger energy
  int multiples_per_count[16]; // Multiples per counts 
  int multiples_per_channel[16]; // Multiples per tagger channel
  vector<hist_struct_TH1D> t_hist; // Histogram for tagger timing distribution
  hist_struct_TH2D h_tagger_vs_energy; // Tagger time vs energy in ECAL
  int cut[16]; // Mean tagger time for cutting times off the mean time
  double energy[16] = {41.3,44.8,51.6,69.14,79.3,99.8,149.37,200.,249.8,350.0,450.2,550.37,650.3,699.7,724.4,725.5}; // Energies corresponding to each tagger channel
};


// Construct containers for storing all wave forms + histograms + various informations 
vector<signal_struct> RAW;
vector<signal_struct> RAW_CALIB;
vector<signal_struct> MA;
vector<signal_struct> MWD;
vector<signal_struct> TMAX;
vector<signal_struct> NMO; // NMO: Nelder-Mead Optimization Algorithm
//
vector<signal_struct> ECAL25; // Sum of all channels, only one element in vector
//
// Proto trace
vector<signal_struct> PROTO;
// Struct to save all multi sampling renormalisation histograms and parameters
// Dimensions: [eff. channels][MULTIS][MULTIS]
vector<vector<vector<multis_norm_struct> > > MULTIS_NORM;
// Build root file
TFile *hfile;
// output stream for writing the proto trace
ofstream *proto_out;
// Build the calibration struct
calib_struct CALIB;
// Build the Mapping struct
vector<vector<mapping_struct> > MAPPING;
// Initialize a tagger
tagger_struct TAGGER;
// Container for costum multiplicities
vector<double> TH_MULTIPLICITY;

// Definition of functions
void extraction();
void multis_calib();
void plot_waves(vector<signal_struct> &array, char const *name, char const *modus);
void plot_waves_compare(char const *name, char const *modus);
void plot_interpol(vector<double> &x, vector<double> &y);
void plot_time_energy(time_struct &array);
void plot_energy_hist(vector<signal_struct> &array, char const *path, const char *mode);
void plot_tagger_hist(vector<signal_struct> &array, char const *path, const char *mode);
void plot_timing_hist(vector<signal_struct> &signal, char const *path);
void plot_multis_hist();
void plot_TH2D_hist(TH2D *hist, char const *path, const char *name);
void plot_TH2D_hist_graph(TH2D *hist, vector<double> trace, char const *path, const char *name);
void plot_energy_vs_tagged(signal_struct &signal, char const *path, const char *name);
void plot_energy_resolution(signal_struct &signal, char const *path, const char *name);
void plot_trace(vector<double> &trace, char const *name, char const *modus);
double randit(int ini=0);
double array_mean(const vector<double> &array, int start, int end);
double array_std(const vector<double> &array, int start, int end, double mean);
int array_largest(vector<double> &array, int lower, int upper);
int array_zero_xing(vector<double> &array, int lower, int upper);
double array_compare(vector<double> &array1, vector<double> &array2, vector<double> &weigths, vector<double> &par, int start, int end, int debug);
vector<double> array_adjust(vector<double> &array, vector<double> &x, int debug);
vector<double> array_sum(vector<double> &array1, vector<double> &array2, double factor);
vector<double> array_simulate_proto();
vector<double> array_smooth(vector<double> &array, int s, int L);
void interpolate(vector<signal_struct> &signal);
void time_compare(vector<signal_struct> &signal);
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
vector<Double_t> fit_graph_err(TGraphErrors *graph, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
void print_usage();
void build_structure();
void print_final_statistics();
void print_energy_statistics(vector<signal_struct> &array, const char *name);
void print_energy_calib();
void print_timing_statistics(time_struct &array, Int_t total_coincidents, const char *name);
void print_stat_multis_calib();
void init_intersamp_hist(int channels);
void init_signal(vector<signal_struct> &signal, int channels, bool is_raw = false);
void init_times(vector<vector<time_struct> > &array, int channels);
void reset_signal(vector<signal_struct> &signal);
void reset_times(vector<vector<time_struct> > &array);
void init_multis_norm(vector<vector<vector<multis_norm_struct> > > &array, int channels);
void fill_hists();
void init_hists(int channels);
bool read_file(string file);
bool read_config(char const *file);
bool linreg(vector<double> &x, vector<double> &y, double *m, double *b);
// Double_t fpeaks(Double_t *x, Double_t *par);
Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, bool silent);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
Int_t largest_1Dbin(TH1D *hist, Int_t lower, Int_t upper);
vector<double> FIR_filter(vector<double> &trace, double calib);
vector<double> MA_filter(vector<double> &trace, double calib);
vector<double> MWD_filter(vector<double> &trace, double calib);
vector<double> CFD(vector<double> &trace, double multiplier);
void print_detector_config();
bool is_in_string(char const *character, char const *letter);
bool is_glitch(vector<double> &trace, double TH, int n);
bool is_saturation(signal_struct &signal, int n);
bool baseline_weird(signal_struct &signal);
int is_valid_max(signal_struct &signal, int n);
double signal_integral(signal_struct &signal, int debug);
double polnx(double x, vector<double> &par);
double log3x(double x, vector<double> &par);
double ExpDecay1(double x, vector<double> &par);
double ExpGro1(double x, vector<double> &par);
Double_t SIPMpixel(Double_t *x, Double_t *par);
Double_t SIPMlinearize(Double_t x, Double_t A);
Double_t resolution(Double_t *x, Double_t *par);

#include "simplex.h"

  
int main(int argc, char *argv[])
{
  // Timer for program execution time measurement
  high_resolution_clock::time_point t_total_begin = high_resolution_clock::now();
  ///////////////////////////////////////////////////////////////
  // Read all existing command line argunments and file names
  ///////////////////////////////////////////////////////////////
  char inputfile[200];
  char outputfile[200];
  char outputfile_proto_trace[200];
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
        printf("Software threshold multiplicy set to: %3.2f !\n", THRESHOLD_MULTIPLICY);
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
    if(strstr(argv[n],"-v")!=NULL){  // turn verbose output on
      n++;
      if(n<argc){
        strcpy(VERBOSE, argv[n]);
        printf("Verbose options set to: %s !\n", VERBOSE);
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

  // Standard detector analysis mode is cosmics mode, is changed by config file
  sprintf(MODE, "COSMICS"); // standard setting

  // Read the config file if file is given
  if(conf_exists){ 
    bool conf_healthy = read_config(configfile);
    if (conf_healthy == false){
      printf("ERROR: Config file error!\n");
      return(-1);
    }
  }

  // Print mapping if verbose flag is set
  if (is_in_string(VERBOSE, "c")) print_detector_config();

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

  // Proto signal is to be extracted
  if (EXTRACT_PROTO == 1){
    sprintf(outputfile_proto_trace, "%s_proto_trace.dat", outputfile);
    proto_out = new ofstream(outputfile_proto_trace);
    if (proto_out->is_open()){
      printf("NOTICE(main): Opened proto trace output file %s.\n", outputfile_proto_trace);
    }
    else printf("WARNING(statistics): Unable to write proto trace file.\n");
  }

  // array_simulate_proto();


  ////////////////////////////////////////////////////////
  // BUILD HISTOGRAMS / LOAD PROTO TRACES FROM FILE
  ////////////////////////////////////////////////////////

  // Reset/initialize signal container
  // If in normal extraction mode, use CHANNELS_EFF
  if (MULTIS_CALIB_MODE == false){
    // Initialize the signal containers
    init_signal(RAW, CHANNELS, "true"); // Has to be CHANNELS 
    init_signal(RAW_CALIB, CHANNELS_EFF);
    init_signal(MA, CHANNELS_EFF);
    init_signal(MWD, CHANNELS_EFF);
    init_signal(TMAX, CHANNELS_EFF);
    init_signal(NMO, CHANNELS_EFF);
    // Init the Calorimeter sum (only one element per filter type in the vector)
    init_signal(ECAL25, 5);
    // Initialize histograms, done in extra functionS
    init_hists(CHANNELS_EFF);
  }
  else{
    init_signal(RAW, CHANNELS);
    // Initialize histograms, done in extra function  
    // Initialize the MULTIS container
    init_multis_norm(MULTIS_NORM, CHANNELS);
    init_hists(CHANNELS);
  }
  // Save the Proto Trace 
  if (EXTRACT_PROTO == 0){
    for (int i = 0; i<(int)RAW_CALIB.size(); i++){
      RAW_CALIB[i].proto_trace_fit = PROTO[i].proto_trace_fit;
    }
  }
  // Save costum threshold multiplicies
  if (TH_MULTIPLICITY.size() != RAW_CALIB.size() && THRESHOLD_MULTIPLICY == -1){
    printf("ERROR: Check your costum threshold multiplicies from the ini file!\n");
    return(-1);
  }
  if (THRESHOLD_MULTIPLICY == -1){
    for (int i = 0; i<(int)RAW_CALIB.size(); i++){
      RAW_CALIB[i].base.TH_multiplicity = TH_MULTIPLICITY[i];
    }
  }
  else{
    for (int i = 0; i<(int)RAW_CALIB.size(); i++){
      RAW_CALIB[i].base.TH_multiplicity = THRESHOLD_MULTIPLICY;
    }
  }


  printf("INITIALIZATION COMPLETED! BEGIN READOUT!\n");



  /////////////////////////////////////////////////////
  // BEGIN SIGNAL READOUT
  /////////////////////////////////////////////////////

  // Initialize the loop condition
  int m=1;
  // Reset NOE number of events counter
  // loop timer
  high_resolution_clock::time_point t_loop_begin = high_resolution_clock::now();
  NOE=0;
  do{
    // Keep reading as long there are unread events
    m=DETECTOR.read_one_event(NOE);
    // Print some verbose information
    if (is_in_string(VERBOSE, "p")){ // "p" for progress report
      if (no_of_events != 0) {
        if(NOE%((int)no_of_events/10)==0){
		      high_resolution_clock::time_point t_loop_end = high_resolution_clock::now();
		      // duration<double, milli> dur = ( t_loop_end - t_loop_begin );
		      auto duration = duration_cast<milliseconds>( t_loop_end - t_loop_begin ).count();
          cout << "Analysing event: " << NOE << " (" << duration/1000 << "s per cycle)" << endl;
        } 
      }
      else { 
        if(NOE%1000 ==0) {
		      high_resolution_clock::time_point t_loop_end = high_resolution_clock::now();
		      auto duration = duration_cast<milliseconds>( t_loop_end - t_loop_begin ).count();
        	cout << "Analysing event: " << NOE << " (" << duration/1000 << "s per cycle)" << endl;
        }
      }
      // reset timer 
      t_loop_begin = high_resolution_clock::now();
    }
    // Increase Global number of events counter
    NOE++;
    // If there is an event, do the extraction
    if(m==1){ // && t!=-13){  // not eof for either of these files
      // MULTISAMPLING CALIBRATION MODE
      // Depending on the program mode, do either calibration or final extraction
      //
      if(strcmp(MODE, "MULTIS") == 0 || MULTIS_CALIB_MODE == true) {
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
      }
      // Every 1000 events plot a random event
      if (NOE%1000==0){
        // Plot every 100th raw signal
        hfile->cd("WAVE_FORMS/RAW"); 
        plot_waves(RAW, "Signal_RAW", "TRACE");  
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

  if(strcmp(MODE, "MULTIS") == 0 || MULTIS_CALIB_MODE == true) {
    print_stat_multis_calib();
    plot_multis_hist();
  }
  else{
    print_final_statistics();
  }

  /////////////////////////////////////////////////////
  // END PHYSICS PROGRAM
  /////////////////////////////////////////////////////

  // End programm timer
  high_resolution_clock::time_point t_total_end = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>( t_total_end - t_total_begin ).count();
  if (NB_ACT_CHANNELS == 0) NB_ACT_CHANNELS++; // Avoid deviding by zero
  cout << endl << "Program exectuion time: " << (duration/1000) 
      << "s (" << duration/NB_ACT_CHANNELS/1000 << "s per active channel)" << endl << endl;

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
  reset_signal(NMO);
  reset_signal(ECAL25);
  /////////////////////////////////////////////////////
  // FETCH RAW CONTENT FROM DATA
  /////////////////////////////////////////////////////
  // Create dummy container for fetching the signal
  int entry = 0;
  unsigned int dumm_cont[CHANNELS][TRACELEN];
  // Initialization of the boards (has to start with 1, since 0 is usually the IOL board)
  //    and fetching the ADC content
  for(int board=1; board<=BOARDS; board++){
    for(int channel=0; channel<8; channel++){
      // Convert the board channel to overall channel 
      entry = (board-1)*8 + channel;
      // Fetch the channels
      DETECTOR.get_trace(board, channel, dumm_cont[entry]);
      // Save the entries in the RAW container
      for(int n=0; n<TRACELEN; n++){
        double dumm = (double) dumm_cont[entry][n]; // Has no function but is needed to work (dunno why)
        RAW[entry].trace.push_back((double) dumm_cont[entry][n]);
        dumm++; // Has no function but is needed to work (dunno why) 
      }
      // Calculate the baseline information for the RAW and RAW_Calib
      RAW[entry].base.mean = array_mean(RAW[entry].trace, 0,BASELINE_CUT);
      // Already subtract the baseline from the signals
      for (int n = 0; n<TRACELEN; n++){
        RAW[entry].trace[n] -= RAW[entry].base.mean;
      }
      // Recalculate the baseline information for the RAW and RAW_Calib
      RAW[entry].base.mean = array_mean(RAW[entry].trace, 0,BASELINE_CUT);
      // Now calculate the std and TH for the baselines
      RAW[entry].base.std = array_std(RAW[entry].trace,0, BASELINE_CUT,RAW[entry].base.mean);
      RAW[entry].base.TH = THRESHOLD_MULTIPLICY * RAW[entry].base.std; 
    }
  }    
  // According to polarity and baseline, transform negative signals into positive signals
  for (int a = 0; a<(int)MAPPING.size(); a++){

    for (int b = 0; b<(int)MAPPING[a].size(); b++){

      // Mapped hardware channel in RAW container
      int h = MAPPING[a][b].board_nb * 8 + MAPPING[a][b].h_channel_start;
      // If it's a solo channel
      if (MAPPING[a][b].multis == 1 && MAPPING[a][b].polarity == -1){
        for (int n = 0; n < (int)RAW[h].trace.size(); n++){
          RAW[h].trace[n] = 2 * RAW[h].base.mean + (RAW[h].trace[n]*(-1));
        }
      } 
      // If not:

      if (MAPPING[a][b].multis == 2 && MAPPING[a][b].polarity == -1){
        for (int n = 0; n < (int)RAW[h].trace.size(); n++){
          RAW[h].trace[n] = 2 * RAW[h].base.mean + (RAW[h].trace[n]*(-1));
          RAW[h+1].trace[n] = 2 * RAW[h+1].base.mean + (RAW[h+1].trace[n]*(-1));
        }
      } 
      if (MAPPING[a][b].multis == 4 && MAPPING[a][b].polarity == -1){
        for (int n = 0; n < (int)RAW[h].trace.size(); n++){
          RAW[h].trace[n] = 2 * RAW[h].base.mean + (RAW[h].trace[n]*(-1));
          RAW[h+1].trace[n] = 2 * RAW[h+1].base.mean + (RAW[h+1].trace[n]*(-1));
          RAW[h+2].trace[n] = 2 * RAW[h+2].base.mean + (RAW[h+2].trace[n]*(-1));
          RAW[h+3].trace[n] = 2 * RAW[h+3].base.mean + (RAW[h+3].trace[n]*(-1));
        }
      } 
    }
  }
  // i is the software channel, including all multisampling channels. mapping from hardware channels has to be done
  for (int a = 0; a<(int)MAPPING.size(); a++){
    for (int b = 0; b<(int)MAPPING[a].size(); b++){
      // Chrystal channel (Channel of matrix)
      int ch = b*(int)MAPPING.size()+a; 
      // Mapped hardware channel in RAW container
      int h = MAPPING[a][b].board_nb * 8 + MAPPING[a][b].h_channel_start;
      // Check if hardware channel is valid channel, otherwise leave that channel out and set the flag
      if (MAPPING[a][b].board_nb != 99) {
        RAW_CALIB[ch].is_valid = true;      
        MA[ch].is_valid = true;      
        MWD[ch].is_valid = true;      
        TMAX[ch].is_valid = true;  
        NMO[ch].is_valid = true;  
        // Aditionally for the RAW channels
        RAW[h].is_valid = true;
        if (MAPPING[a][b].multis == 2){
          RAW[h+1].is_valid = true;
        }
        if (MAPPING[a][b].multis == 4){
          RAW[h+1].is_valid = true;
          RAW[h+2].is_valid = true;
          RAW[h+3].is_valid = true;
        }
      }
      else continue;
      // Now fill the software channels according to the mapping and multisampling
      for(int n=0; n<TRACELEN; n++){
        if (MAPPING[a][b].multis == 1){
          RAW_CALIB[ch].trace.push_back( ((double) RAW[h].trace[n] )* CALIB.multis[h] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          // Update the channel info
          RAW_CALIB[ch].hardware_addr = h;
          RAW_CALIB[ch].polarity = MAPPING[a][b].polarity;
        }
        else if (MAPPING[a][b].multis == 2){
          // Save the samples correctly 
          RAW_CALIB[ch].trace.push_back((double) RAW[h+1].trace[n] * CALIB.multis[h+1] * CALIB.RAW_energy[ch] * GENERAL_SCALING);       
          RAW_CALIB[ch].trace.push_back((double) RAW[h].trace[n] * CALIB.multis[h] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          // Update the channel information
          RAW_CALIB[ch].multis = 2;
          RAW_CALIB[ch].clock_speed = 200; // in MHz
          RAW_CALIB[ch].sample_t = 5.; // in MHz
          RAW_CALIB[ch].tracelen = TRACELEN * 2;
          RAW_CALIB[ch].hardware_addr = h;
          RAW_CALIB[ch].polarity = MAPPING[a][b].polarity;
        }
        else if (MAPPING[a][b].multis == 4){
          // Save the samples correctly
          RAW_CALIB[ch].trace.push_back((double) RAW[h+3].trace[n] * CALIB.multis[h+3] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          RAW_CALIB[ch].trace.push_back((double) RAW[h+2].trace[n] * CALIB.multis[h+2] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          RAW_CALIB[ch].trace.push_back((double) RAW[h+1].trace[n] * CALIB.multis[h+1] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          RAW_CALIB[ch].trace.push_back((double) RAW[h].trace[n] * CALIB.multis[h] * CALIB.RAW_energy[ch] * GENERAL_SCALING);
          // Update the channel information
          RAW_CALIB[ch].multis = 4;
          RAW_CALIB[ch].clock_speed = 400; // in MHz
          RAW_CALIB[ch].sample_t = 2.5; // in MHz
          RAW_CALIB[ch].tracelen = TRACELEN * 4;
          RAW_CALIB[ch].hardware_addr = h;
          RAW_CALIB[ch].polarity = MAPPING[a][b].polarity;
        }
      }
    }
  }

  for(int board=1; board<=BOARDS; board++){
    for(int channel=0; channel<8; channel++){
      // Convert the board channel to overall channel 
      entry = (board-1)*8 + channel;
      RAW[entry].CFD.trace.clear();
      RAW[entry].CFD.trace = CFD(RAW[entry].trace, 1);
    }
  }

  for (int i = 0; i < (int)RAW_CALIB.size(); i++){
	  /////////////////////////////////////////////////////
	  // APPLICATION OF THE FILTER TO RAW_CALIB 
	  /////////////////////////////////////////////////////
	  MA[i].trace.clear();
	  MWD[i].trace.clear();
	  TMAX[i].trace.clear();
	  MA[i].trace = MA_filter(RAW_CALIB[i].trace, CALIB.MA_energy[i]);
	  MWD[i].trace = MWD_filter(RAW_CALIB[i].trace, CALIB.MWD_energy[i]);
	  TMAX[i].trace = FIR_filter(RAW_CALIB[i].trace, CALIB.TMAX_energy[i]);
    for (int n = 0; n < (int)TMAX[i].trace.size(); n++){
      // printf("%3.1f %3.1f %3.1f\n", RAW_CALIB[i].trace[n], TMAX[i].trace[n], CALIB.TMAX_energy[i] );
    }

	  /////////////////////////////////////////////////////
	  // APPLICATION OF THE CONSTANT FRACTION DISCRIMINATOR 
	  /////////////////////////////////////////////////////
	  RAW_CALIB[i].CFD.trace.clear();
	  MA[i].CFD.trace.clear();
	  MWD[i].CFD.trace.clear();
	  TMAX[i].CFD.trace.clear();
	  RAW_CALIB[i].CFD.trace = CFD(RAW_CALIB[i].trace, 1);
	  MA[i].CFD.trace = CFD(MA[i].trace, 1);
	  MWD[i].CFD.trace = CFD(MWD[i].trace, 1);
	  TMAX[i].CFD.trace = CFD(TMAX[i].trace, 1);
	}

  /////////////////////////////////////////////////////
  // CALCULATE BASELINE STATISTICS 
  /////////////////////////////////////////////////////
  // Calculate the baseline statistics
  for (int i = 0; i < (int)RAW_CALIB.size(); i++){
    // Check if channel is valid
    if (RAW_CALIB[i].is_valid == false) continue;
    // otherwise do the calculation
    RAW_CALIB[i].base.mean = array_mean(RAW_CALIB[i].trace, 0,BASELINE_CUT*RAW_CALIB[i].multis);
    RAW_CALIB[i].base.std = array_std(RAW_CALIB[i].trace, 0, BASELINE_CUT*RAW_CALIB[i].multis,RAW_CALIB[i].base.mean);
    RAW_CALIB[i].base.TH = RAW_CALIB[i].base.TH_multiplicity * RAW_CALIB[i].base.std; 
    //
    MA[i].base.mean = array_mean(MA[i].trace, 0,BASELINE_CUT*MA[i].multis);
    MA[i].base.std = array_std(MA[i].trace, 0, BASELINE_CUT*MA[i].multis,MA[i].base.mean);
    MA[i].base.TH = RAW_CALIB[i].base.TH_multiplicity * MA[i].base.std; 
    //
    MWD[i].base.mean = array_mean(MWD[i].trace, 0,BASELINE_CUT*MWD[i].multis);
    MWD[i].base.std = array_std(MWD[i].trace, 0, BASELINE_CUT*MWD[i].multis,MWD[i].base.mean);
    MWD[i].base.TH = RAW_CALIB[i].base.TH_multiplicity * MWD[i].base.std; 
    //
    TMAX[i].base.mean = array_mean(TMAX[i].trace, 0,BASELINE_CUT*TMAX[i].multis);
    TMAX[i].base.std = array_std(TMAX[i].trace, 0, BASELINE_CUT*TMAX[i].multis,TMAX[i].base.mean);
    TMAX[i].base.TH = RAW_CALIB[i].base.TH_multiplicity * TMAX[i].base.std; 
  }
  /////////////////////////////////////////////////////
  // EXTRACT FEATURES FROM TRACES
  /////////////////////////////////////////////////////
  // First the RAW traces
  for(int i=0; i<(int)RAW.size(); i++){
    // Extract the maximum
  	RAW[i].energy = RAW[i].trace[array_largest(RAW[i].trace, BASELINE_CUT, ENERGY_WINDOW_MAX)];
  	// Look for CFD Zero Crossing
  	RAW[i].CFD.Xzero = array_zero_xing(RAW[i].CFD.trace, BASELINE_CUT, ENERGY_WINDOW_MAX);
    // Already fill the samples into a histogram for baseline noise estimation
    for(int n=0; n < BASELINE_CUT; n++){
      RAW[i].base.h_samples.hist->Fill(RAW[i].trace[n]);   
    }
  }

	int plot = 0;

  //
  // Then the Filtered traces
  // 
  // In beam mode first check if the central crystal is healthy
  int central_healty = 0;
  // Left and Right borders for maximum detection
  int left = BASELINE_CUT * RAW_CALIB[0].multis;
  int right = ENERGY_WINDOW_MAX * RAW_CALIB[0].multis;
  // Already extract features of RAW_CALIB
  if (strcmp(MODE, "BEAM")==0){
    RAW_CALIB[CENTRAL].energy_n = array_largest(RAW_CALIB[CENTRAL].trace, left, right); // sample number of largest sample
    RAW_CALIB[CENTRAL].energy = RAW_CALIB[CENTRAL].trace[RAW_CALIB[CENTRAL].energy_n];
    RAW_CALIB[CENTRAL].integral = signal_integral(RAW_CALIB[CENTRAL], 0);
    RAW_CALIB[CENTRAL].ratio = RAW_CALIB[CENTRAL].integral / RAW_CALIB[CENTRAL].energy;
    if ( is_valid_max( RAW_CALIB[CENTRAL], RAW_CALIB[CENTRAL].energy_n ) == 0 ){
      if ( RAW_CALIB[CENTRAL].ratio < UPPER_RATIO && RAW_CALIB[CENTRAL].ratio > LOWER_RATIO){
        RAW_CALIB[CENTRAL].is_signal = 1;
        central_healty = 1;
      }
      else central_healty = 0;
    }
    else central_healty = 0;

  }
  else central_healty = 1;
  // Now loop through all the channels
	for(int i=0; i<(int)RAW_CALIB.size(); i++){
    // plot = 1;
    // Only extract features from valid channels
    if (RAW_CALIB[i].is_valid == false) continue;
		// Left and Right borders for maximum detection
		int left = BASELINE_CUT * RAW_CALIB[i].multis;
		int right = ENERGY_WINDOW_MAX * RAW_CALIB[i].multis;
		// Already extract features of RAW_CALIB
    RAW_CALIB[i].energy_n = array_largest(RAW_CALIB[i].trace, RAW_CALIB[i].multis*BASELINE_CUT, RAW_CALIB[i].multis*ENERGY_WINDOW_MAX); // sample number of largest sample
    RAW_CALIB[i].energy = RAW_CALIB[i].trace[RAW_CALIB[i].energy_n];
    RAW_CALIB[i].integral = signal_integral(RAW_CALIB[i], 0);
    RAW_CALIB[i].ratio = RAW_CALIB[i].integral / RAW_CALIB[i].energy;
    // printf("%3.1f %3.1f %d\n", RAW_CALIB[i].energy, RAW_CALIB[i].base.TH, RAW_CALIB[i].energy_n);

    // Search for maximum and check if valid
    if ( is_valid_max( RAW_CALIB[i], RAW_CALIB[i].energy_n ) == 0 && central_healty == 1 ){ // The maximum is valid
      // Signal is good
      RAW_CALIB[i].is_signal = 1;
      // Save the sum of all signals per channel
      RAW_CALIB[i].proto_trace = array_sum(RAW_CALIB[i].proto_trace, RAW_CALIB[i].trace, 1);
      // Fill the proto trace hist
      for (int n = 0; n < (int)RAW_CALIB[i].proto_trace.size(); n++){
        RAW_CALIB[i].TH2D_proto_trace.hist->Fill(n, RAW_CALIB[i].trace[n]);
      }
      //
    	// Only look for the filter traces if raw trace is not a glitch
      MA[i].energy_n = array_largest(MA[i].trace, left, right);
      if ( is_valid_max( MA[i], MA[i].energy_n ) == 0 ){
        MA[i].energy = MA[i].trace[MA[i].energy_n];
        MA[i].integral = signal_integral(MA[i], 0);
        MA[i].ratio = MA[i].integral / MA[i].energy;
        MA[i].is_signal = 1;
      }
      else{ MA[i].energy = 0; MA[i].integral = 0; MA[i].ratio = 0; MA[i].is_signal = 0;}
      //
      MWD[i].energy_n = array_largest(MWD[i].trace, left, right);
      if ( is_valid_max( MWD[i], MWD[i].energy_n ) == 0 ){
        MWD[i].energy = MWD[i].trace[MWD[i].energy_n];
        MWD[i].integral = signal_integral(MWD[i], 0);
        MWD[i].ratio = MWD[i].integral / MWD[i].energy;
        MWD[i].is_signal = 1;
      }
      else{ MWD[i].energy = 0; MWD[i].integral = 0; MWD[i].ratio = 0; MWD[i].is_signal = 0;}
      //
      TMAX[i].energy_n = array_largest(TMAX[i].trace, left, right);
      if ( is_valid_max( TMAX[i], TMAX[i].energy_n ) == 0 ){
        TMAX[i].energy = TMAX[i].trace[TMAX[i].energy_n];
        TMAX[i].integral = signal_integral(TMAX[i], 0);
        TMAX[i].ratio = TMAX[i].integral / TMAX[i].energy;
        TMAX[i].is_signal = 1;
      }
      else{ TMAX[i].energy = 0; TMAX[i].integral = 0; TMAX[i].ratio = 0; TMAX[i].is_signal = 0;}
      // For the sake of computational efficiency, only do the Nelder-Mead optimization for signal events

      // If the proto traces are already extracted, do the optimization
      if (EXTRACT_PROTO==0){
        // Now do optimization for prototraces and RAW_CALIB


        vector<double> init; init.push_back(1.0); init.push_back(1.0);

        double a0[]={+1, -1}, a1[]={0.5, -2}, a2[]={ 1.5, 1.0};
        vector<vector<double> > simplex;
        simplex.push_back( vector<double>(a0,a0+2) );
        simplex.push_back( vector<double>(a1,a1+2) );
        simplex.push_back( vector<double>(a2,a2+2) );

        vector<double> weights;

        for (int n = 0; n < 100; n++){
          weights.push_back(1);
        }
        for (int n = 100; n < 120; n++){
          weights.push_back(1);
        }
        for (int n = 120; n < (int)RAW_CALIB[i].trace.size(); n++){
          weights.push_back(1);
        }

        vector<double> opt;
        // printf("\n%d %d\n", (int)RAW_CALIB[i].trace.size(), (int)RAW_CALIB[i].proto_trace_fit.size());
        int debug = 0;
        // if (i == 6) debug = 2;
        opt = BT::Simplex(array_compare, 
                      RAW_CALIB[i].trace,
                      RAW_CALIB[i].proto_trace_fit,
                      weights,
                      50, 250,
                      init,
                      (double) 1e-7,
                      simplex,
                      1E5,
                      debug
        );
        // printf("++++++++++\n %3.3f %3.3f\n++++++++++", opt[0], opt[1]);
        // Save optimized wave form in NMO trace
        
        NMO[i].trace = array_adjust(RAW_CALIB[i].proto_trace_fit, opt, 0);
        NMO[i].is_signal = 1;
        NMO[i].is_valid = true;

        int len1 = (int)NMO[i].trace.size();
        int len2 = (int)RAW_CALIB[i].trace.size();
        double frac = (double)len1/(double)len2;


        NMO[i].energy_n = array_largest(NMO[i].trace, (int)left*frac, (int)right*frac); // sample number of largest sample
        NMO[i].energy = NMO[i].trace[NMO[i].energy_n];
        NMO[i].integral = signal_integral(NMO[i], 0);
        NMO[i].ratio = NMO[i].integral / NMO[i].energy;

        // if (i == CENTRAL ) plot = 1;

        // plot = 1;

      }
      else{
        // For the NMO, fill the non-optimized traces with 0
        for (int n = 0; n < (int)RAW_CALIB[0].trace.size(); n++){
          NMO[i].trace.push_back(0.0);
        }
      }
    }
    else { 
    	// If not valid, set all energies and integrals to 0
    	RAW_CALIB[i].energy = 0; RAW_CALIB[i].integral = 0;RAW_CALIB[i].ratio = 0; RAW_CALIB[i].is_signal = 0;
    	MA[i].energy = 0; MA[i].integral = 0; MA[i].ratio = 0; MA[i].is_signal = 0;
    	MWD[i].energy = 0; MWD[i].integral = 0; MWD[i].ratio = 0; MWD[i].is_signal = 0;
    	TMAX[i].energy = 0; TMAX[i].integral = 0; TMAX[i].ratio = 0; TMAX[i].is_signal = 0;
      NMO[i].energy = 0; NMO[i].integral = 0; NMO[i].ratio = 0; NMO[i].is_signal = 0;
      // For the NMO, fill the non-optimized traces with 0
      for (int n = 0; n < (int)RAW_CALIB[0].trace.size(); n++){
        NMO[i].trace.push_back(0.0);
      }
    }

    // Calculate sample_t of NMO which can be different than 
    double nmo_tracelen = (double)NMO[i].trace.size();
    NMO[i].sample_t = (double)( SAMPLE_t * (TRACELEN/nmo_tracelen) ); // NMO can have a higher sampling rate, calculate it based on the ADC sampling
    // printf("%3.2f\n", NMO[i].sample_t);

    // Calculate the CFD traces
    RAW_CALIB[i].CFD.trace.clear();
    RAW_CALIB[i].CFD.trace = CFD(RAW_CALIB[i].trace, 1);
    MA[i].CFD.trace.clear();
    MA[i].CFD.trace = CFD(MA[i].trace, 1);
    MWD[i].CFD.trace.clear();
    MWD[i].CFD.trace = CFD(MWD[i].trace, 1);
    TMAX[i].CFD.trace.clear();
    TMAX[i].CFD.trace = CFD(TMAX[i].trace, 1);
    // Calculate the multiplication adjustment to the higher res NMO trace
    double multiplier = (double)(NMO[i].trace.size() / RAW_CALIB[i].trace.size());
    //
    NMO[i].CFD.trace.clear();
    NMO[i].CFD.trace = CFD(NMO[i].trace, multiplier);
  

    // if ( RAW_CALIB[i].ratio == 0 ){
    // 	if ( i == 4 ) plot = 1;
    // }
    // if ( is_valid_max( RAW_CALIB[i], array_largest(RAW_CALIB[i].trace, left, right) ) == 4 ){
    // 	if ( i == 4 ) plot = 2;
    // }

  }



  /////////////////////////////////////////////////////
  // TEST FOR COINCIDENCE
  /////////////////////////////////////////////////////
  // Check for coincidences according to the analysis mode
  // + Column cut: only look for coincident events in stacked crystals

  // Column cut / cosmics mode:
  bool is_coinc = false;
  //
  if(strcmp(MODE, "COSMICS") == 0) {
    // Array for saving column information
    int column = 0; 
    // Check for each column if there was any event
    for (int b = 0; b < (int)MAPPING[0].size(); b++){
      column = 0;
      for (int a = 0; a < (int)MAPPING.size(); a++){
        // Chrystal channel (Channel of matrix)
        int ch = b*(int)MAPPING.size()+a; 
        // Check if there was a signal in channel ch
        if (RAW_CALIB[ch].is_signal == 1 && RAW_CALIB[ch].is_valid == 1){
          column++;
        }
      }
      // printf("%d %d\n", column, (int)MAPPING.size() - 1);
      // Now check how many signals per column
      if (column < COINC_LEVEL ){//(int)MAPPING[0].size()-1){
        // Throw away events in this column 
        for (int a = 0; a < (int)MAPPING.size(); a++){
          // Chrystal channel (Channel of matrix)
          int ch = b*(int)MAPPING.size()+a; 
          RAW_CALIB[ch].energy = 0;
          MA[ch].energy = 0;
          MWD[ch].energy = 0;
          TMAX[ch].energy = 0;
          NMO[ch].energy = 0;
          //
          RAW_CALIB[ch].integral = 0;
          MA[ch].integral = 0;
          MWD[ch].integral = 0;
          TMAX[ch].integral = 0;
          NMO[ch].integral = 0;
        }
      }
      else is_coinc = true;
      // If it's only one crystal per column, set the coincidencec setting still to true
      if ( (int)MAPPING.size() == 1){
        is_coinc = true;
      }
      // If there are coincident events, do timing extraction
      else{
        // start the interpolation by looping over all channels
        // interpolate is searching for the zero crossing point of the CFD signals
        // interpolate(RAW_CALIB);
        // interpolate(MA);
        // interpolate(MWD);
        // interpolate(TMAX);
        // // And now compare the timing
        // time_compare(RAW_CALIB);
        // time_compare(MA);
        // time_compare(MWD);
        // time_compare(TMAX);
      }
    }
  } 
  if (plot == true){
    printf("PLOTTED\n");
    hfile->cd("JUNK/");
    plot_waves_compare("PLOT1", "TRACE");
  }
  if (plot == 2){
    hfile->cd("JUNK/");
    plot_waves_compare("PLOT2", "TRACE");
  }

  //
  //
  if(strcmp(MODE, "BEAM") == 0) {
    is_coinc = false;
    // If xtal in beam is not valid
    // Read out the tagger statistics
    int tags_per_event = 0;
    for(int n=0; n<TAG_CHANNELS; n++){
      TAGGER.time[n]=DETECTOR.get_value(BOARDS+1, n, 1);
      // Check for multiple taggs or non tags
      if (TAGGER.time[n] != 0.0 ){
        tags_per_event++;
        TAGGER.counts[n]++;
      } 
    }
    // Enter number of tags per event into the counter
    TAGGER.multiples_per_count[tags_per_event]++;
    // Check in witch channel the multi counts appeared
    if (tags_per_event > 1){
      for(int n=0; n<TAG_CHANNELS; n++){
        if (TAGGER.time[n] != 0.0 ){
          TAGGER.multiples_per_channel[n]++;
        } 
      }
    }
    // Enter the tagged energy in the right histogram and do the timing
    for (int k = 0; k < N_E_WINDOW; k ++){
      // If the right tagger channel is found
      for (int i = 0; i < (int)RAW_CALIB.size(); i++){
        // Check for the central crystal
        if ( i==CENTRAL && RAW_CALIB[i].is_signal == 1 ) is_coinc = true; 
        // Fill all events for the right tagger
        if ( RAW_CALIB[i].is_signal == true && // If signal is valid
             TAGGER.time[k] != 0 // and if tagger is set for energy k
              ){ 
          RAW_CALIB[i].tagged[k].energy = RAW_CALIB[i].energy;
          MA[i].tagged[k].energy = MA[i].energy;
          MWD[i].tagged[k].energy = MWD[i].energy;
          TMAX[i].tagged[k].energy = TMAX[i].energy;
          NMO[i].tagged[k].energy = NMO[i].energy;
          //
          RAW_CALIB[i].tagged[k].integral = RAW_CALIB[i].integral;
          MA[i].tagged[k].integral = MA[i].integral;
          MWD[i].tagged[k].integral = MWD[i].integral;
          TMAX[i].tagged[k].integral = TMAX[i].integral;
          NMO[i].tagged[k].integral = NMO[i].integral;
          //
          ECAL25[0].tagged[k].energy += RAW_CALIB[i].energy;
          ECAL25[1].tagged[k].energy += MA[i].energy;
          ECAL25[2].tagged[k].energy += MWD[i].energy;
          ECAL25[3].tagged[k].energy += TMAX[i].energy;
          ECAL25[4].tagged[k].energy += NMO[i].energy;
          //
          ECAL25[0].tagged[k].integral += RAW_CALIB[i].integral;
          ECAL25[1].tagged[k].integral += MA[i].integral;
          ECAL25[2].tagged[k].integral += MWD[i].integral;
          ECAL25[3].tagged[k].integral += TMAX[i].integral;
          ECAL25[4].tagged[k].integral += NMO[i].integral;
        }
        else{
          RAW_CALIB[i].tagged[k].energy = 0.0;
          MA[i].tagged[k].energy = 0.0;
          MWD[i].tagged[k].energy = 0.0;
          TMAX[i].tagged[k].energy = 0.0;
          NMO[i].tagged[k].energy = 0.0;
          //
          RAW_CALIB[i].tagged[k].integral = 0.0;
          MA[i].tagged[k].integral = 0.0;
          MWD[i].tagged[k].integral = 0.0;
          TMAX[i].tagged[k].integral = 0.0;
          NMO[i].tagged[k].integral = 0.0;
        }
        // Fill energies only if no multiples are detected
        if ( RAW_CALIB[i].is_signal == true && // If signal is valid
             TAGGER.time[k] != 0 && // and if tagger is set for energy k
             tags_per_event == 1 // and if only one tagged energy per tag
              ){    
          RAW_CALIB[i].tagged[k].energy_m = RAW_CALIB[i].energy;
          MA[i].tagged[k].energy_m = MA[i].energy;
          MWD[i].tagged[k].energy_m = MWD[i].energy;
          TMAX[i].tagged[k].energy_m = TMAX[i].energy;
          NMO[i].tagged[k].energy_m = NMO[i].energy;
          //
          RAW_CALIB[i].tagged[k].integral_m = RAW_CALIB[i].integral;
          MA[i].tagged[k].integral_m = MA[i].integral;
          MWD[i].tagged[k].integral_m = MWD[i].integral;
          TMAX[i].tagged[k].integral_m = TMAX[i].integral;
          NMO[i].tagged[k].integral_m = NMO[i].integral;
          //
          ECAL25[0].tagged[k].energy_m += RAW_CALIB[i].tagged[k].energy_m;
          ECAL25[1].tagged[k].energy_m += MA[i].tagged[k].energy_m;
          ECAL25[2].tagged[k].energy_m += MWD[i].tagged[k].energy_m;
          ECAL25[3].tagged[k].energy_m += TMAX[i].tagged[k].energy_m;
          ECAL25[4].tagged[k].energy_m += NMO[i].tagged[k].energy_m;
          //
          ECAL25[0].tagged[k].integral_m += RAW_CALIB[i].tagged[k].integral_m;
          ECAL25[1].tagged[k].integral_m += MA[i].tagged[k].integral_m;
          ECAL25[2].tagged[k].integral_m += MWD[i].tagged[k].integral_m;
          ECAL25[3].tagged[k].integral_m += TMAX[i].tagged[k].integral_m;
          ECAL25[4].tagged[k].integral_m += NMO[i].tagged[k].integral_m;

        }
        else{
          RAW_CALIB[i].tagged[k].energy_m = 0.0;
          MA[i].tagged[k].energy_m = 0.0;
          MWD[i].tagged[k].energy_m = 0.0;
          TMAX[i].tagged[k].energy_m = 0.0;
          NMO[i].tagged[k].energy_m = 0.0;
          //
          RAW_CALIB[i].tagged[k].integral_m = 0.0;
          MA[i].tagged[k].integral_m = 0.0;
          MWD[i].tagged[k].integral_m = 0.0;
          TMAX[i].tagged[k].integral_m = 0.0;
          NMO[i].tagged[k].integral_m = 0.0;
        }
        // Fill energies for constrained timing window and without multiples
        if ( RAW_CALIB[i].is_signal == true && // If signal is valid
             TAGGER.time[k] != 0 && // and if tagger is set for energy k
             tags_per_event == 1 && // and if only one tagged energy per tag
             TAGGER.time[k] >= (TAGGER.cut[k] - TAGGER_WINDOW) &&
             TAGGER.time[k] <= (TAGGER.cut[k] + TAGGER_WINDOW) // and if tagger time is in the timing window cut
              ){ 
          RAW_CALIB[i].tagged[k].energy_mt = RAW_CALIB[i].energy;
          MA[i].tagged[k].energy_mt = MA[i].energy;
          MWD[i].tagged[k].energy_mt = MWD[i].energy;
          TMAX[i].tagged[k].energy_mt = TMAX[i].energy;
          NMO[i].tagged[k].energy_mt = NMO[i].energy;
          //
          RAW_CALIB[i].tagged[k].integral_mt = RAW_CALIB[i].integral;
          MA[i].tagged[k].integral_mt = MA[i].integral;
          MWD[i].tagged[k].integral_mt = MWD[i].integral;
          TMAX[i].tagged[k].integral_mt = TMAX[i].integral;
          NMO[i].tagged[k].integral_mt = NMO[i].integral;
          //
          ECAL25[0].tagged[k].energy_mt += RAW_CALIB[i].tagged[k].energy_mt;
          ECAL25[1].tagged[k].energy_mt += MA[i].tagged[k].energy_mt;
          ECAL25[2].tagged[k].energy_mt += MWD[i].tagged[k].energy_mt;
          ECAL25[3].tagged[k].energy_mt += TMAX[i].tagged[k].energy_mt;
          ECAL25[4].tagged[k].energy_mt += NMO[i].tagged[k].energy_mt;
          //
          ECAL25[0].tagged[k].integral_mt += RAW_CALIB[i].tagged[k].integral_mt;
          ECAL25[1].tagged[k].integral_mt += MA[i].tagged[k].integral_mt;
          ECAL25[2].tagged[k].integral_mt += MWD[i].tagged[k].integral_mt;
          ECAL25[3].tagged[k].integral_mt += TMAX[i].tagged[k].integral_mt;
          ECAL25[4].tagged[k].integral_mt += NMO[i].tagged[k].integral_mt;

          // if  (k == 0 && 
          //     (NOE % 10) == 0 &&
          //     i == 0 
          //   ){
          //   doof++;
          //   hfile->cd("JUNK");
          //   plot_waves_compare("k0", "TRACE");
          //   printf("CANDIDATE: %d\n", doof);
          //   for (int t = 0; t < N_E_WINDOW; t ++){ 
          //     printf("%3.1f ", TAGGER.time[t]);
          //   }
          //   printf("\n\n");
          // }
        }
        else{
          RAW_CALIB[i].tagged[k].energy_mt = 0.0;
          MA[i].tagged[k].energy_mt = 0.0;
          MWD[i].tagged[k].energy_mt = 0.0;
          TMAX[i].tagged[k].energy_mt = 0.0;
          NMO[i].tagged[k].energy_mt = 0.0;
          //
          RAW_CALIB[i].tagged[k].integral_mt = 0.0;
          MA[i].tagged[k].integral_mt = 0.0;
          MWD[i].tagged[k].integral_mt = 0.0;
          TMAX[i].tagged[k].integral_mt = 0.0;
          NMO[i].tagged[k].integral_mt = 0.0;
        }
      }
    }
    // After all single channels were looked at, now look at the detector sum
    //
    // First for the general untagged energies
    // If central crystal has an event
    if (RAW_CALIB[CENTRAL].is_signal == 1){
      // look for the whole matrix (or submatrix)
      for (int i = 0; i < (int)RAW_CALIB.size(); i++){
        // If signal, sum up
        if (RAW_CALIB[i].is_signal == 1){
          //
          ECAL25[0].energy += RAW_CALIB[i].energy; 
          ECAL25[0].integral += RAW_CALIB[i].integral; 
          //
          for (int k = 0; k < N_E_WINDOW; k++){

          }
        }
        if (MA[i].is_signal == 1){
          ECAL25[1].energy += MA[i].energy; 
          ECAL25[1].integral += MA[i].integral; 
        }
        if (MWD[i].is_signal == 1){
          ECAL25[2].energy += MWD[i].energy; 
          ECAL25[2].integral += MWD[i].integral; 
        }
        // printf("%d %d %d\n", i, (int)RAW_CALIB.size(),(int)TMAX.size());
        if (TMAX[i].is_signal == 1){
          ECAL25[3].energy += TMAX[i].energy; 
          ECAL25[3].integral += TMAX[i].integral; 
        }
        if (NMO[i].is_signal == 1){
          ECAL25[4].energy += NMO[i].energy; 
          ECAL25[4].integral += NMO[i].integral; 
        }
      }
    }
    // If central crystal does not have a valid event, don't sum up
    else{
      for (int i = 0; i < (int)ECAL25.size(); i++){ 
        //
        ECAL25[i].energy = 0; 
        ECAL25[i].integral = 0; 
        //
        for (int k = 0; k < N_E_WINDOW; k++){
          ECAL25[i].tagged[k].energy = 0;
          ECAL25[i].tagged[k].energy_m = 0;
          ECAL25[i].tagged[k].energy_mt = 0;
          ECAL25[i].tagged[k].integral = 0;
          ECAL25[i].tagged[k].integral_m = 0;
          ECAL25[i].tagged[k].integral_mt = 0;
        }
      }
    }
    // Now enter the the sum energy of the ECAL into the tagger time vs ecal energy 2D histogram
    // Only if tagger cut and multiplicity is obeyed
    tags_per_event = 0;
    int tagger_channel = 0;
    for (int k = 0; k < N_E_WINDOW; k++){
      if ( 
        TAGGER.time[k] != 0 && // and if tagger is set for energy k
        TAGGER.time[k] >= (TAGGER.cut[k] - TAGGER_WINDOW) &&
        TAGGER.time[k] <= (TAGGER.cut[k] + TAGGER_WINDOW) // and if tagger time is in the timing window cut
      ){ 
        // Check if only one tagger channel fired
        tagger_channel = TAGGER.time[k];
        tags_per_event++;
      }
    }
    // Only fill if only one tag per event
    if ( tags_per_event == 1){
      if (ECAL25[4].energy != 0) TAGGER.h_tagger_vs_energy.hist->Fill( ECAL25[4].energy, tagger_channel);
    }

  } 
  // Now print wave forms
  // print coincident waveforms
  if( (strcmp(MODE, "COSMICS") == 0 && is_coinc == true) || // either Cosmics mode or
      (strcmp(MODE, "BEAM") == 0 && is_coinc == true) || // Beam mode or
      (strcmp(MODE, "PULSER") == 0) ) { // Pulser mode
    if(NOE%1000==0){
      hfile->cd("WAVE_FORMS/");
      plot_waves_compare("FILTER_COMPARE", "TRACE");
      hfile->cd("WAVE_FORMS/RAW_CALIB");
      plot_waves(RAW_CALIB, "SIGNAL_RAW_CALIB", "TRACE");
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/MA");
      plot_waves(MA, "SIGNAL_MA", "TRACE");
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/MWD");
      plot_waves(MWD, "SIGNAL_MWD", "TRACE");
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/TMAX");
      plot_waves(TMAX, "SIGNAL_TMAX", "TRACE");
      // Set printing folder to MA coincident waves
      hfile->cd("WAVE_FORMS/NMO");
      plot_waves(NMO, "SIGNAL_NMO", "TRACE");
      // Set printing folder to CFD 
      hfile->cd("WAVE_FORMS/CFD/RAW_CALIB");
      plot_waves(RAW_CALIB, "CFD_RAW_CALIB", "CFD");
      hfile->cd("WAVE_FORMS/CFD/MA");
      plot_waves(MA, "CFD_MA", "CFD");
      hfile->cd("WAVE_FORMS/CFD/MWD");
      plot_waves(MWD, "CFD_MWD", "CFD");
      hfile->cd("WAVE_FORMS/CFD/TMAX");
      plot_waves(TMAX, "CFD_TMAX", "CFD");
      hfile->cd("WAVE_FORMS/CFD/NMO");
      plot_waves(NMO, "CFD_NMO", "CFD");

      // if( (strcmp(MODE, "COSMICS") == 0 && is_coinc == true ) && ((int)MAPPING.size() > 1) ){
      //   // The interpolation
      //   hfile->cd("WAVE_FORMS/CFD/MA/INTERPOL");
      //   plot_interpol(MA[0].CFD.x_interpol, MA[0].CFD.y_interpol);
      // }
    }
  }
} 

/////////////////////////////////////////////////////
// MULTI SAMPLING CALIBRATION EXTRACTION
/////////////////////////////////////////////////////
// Extract information from the raw data for the inter sampling calibration
void multis_calib(){
  /////////////////////////////////////////////////////
  // FETCH RAW CONTENT FROM DATA
  /////////////////////////////////////////////////////
  reset_signal(RAW);
  // Create dummy container for fetching the signal
  int entry = 0;
  unsigned int dumm_cont[CHANNELS][TRACELEN];
  // Initialization of the boards (has to start with 1, since 0 is usually the IOL board)
  //    and fetching the ADC content
  for(int board=1; board<=BOARDS; board++){
    for(int channel=0; channel<8; channel++){
      // Convert the board channel to overall channel 
      entry = (board-1)*8 + channel;
      // Fetch the channels
      DETECTOR.get_trace(board, channel, dumm_cont[entry]);
      // Save the entries in the RAW container
      for(int n=0; n<TRACELEN; n++){
        RAW[entry].trace.push_back((double) dumm_cont[entry][n] * CALIB.multis[entry]);
      }
      // Calculate the baseline information for the RAW and RAW_Calib
      RAW[entry].base.mean = array_mean(RAW[entry].trace,0,BASELINE_CUT);
      // Already subtract the baseline from the signals
      for (int n = 0; n<TRACELEN; n++){
        RAW[entry].trace[n] -= RAW[entry].base.mean;
      }
      // Recalculate the baseline information for the RAW and RAW_Calib
      // RAW[entry].base.mean = array_mean(0,BASELINE_CUT,RAW[entry].trace);
      // Now calculate the std and TH for the baselines
      RAW[entry].base.std = array_std(RAW[entry].trace, 0, BASELINE_CUT,RAW[entry].base.mean);
      RAW[entry].base.TH = THRESHOLD_MULTIPLICY * RAW[entry].base.std; 
    }
  }
  // Number of multisampled effective channels
  int multis_nb = 0;
  // i is the software channel, including all multisampling channels. mapping from hardware channels has to be done
  for (int a = 0; a<(int)MAPPING.size(); a++){
    for (int b = 0; b<(int)MAPPING[a].size(); b++){
      // Mapped hardware channel in RAW container
      int h = MAPPING[a][b].board_nb * 8 + MAPPING[a][b].h_channel_start;
      // Now fill the software channels according to the mapping and multisampling
      if (MAPPING[a][b].multis == 1){
        // 
        RAW[h].multis = 1;
      }
      else if (MAPPING[a][b].multis == 2){
        //
        RAW[h].multis = 2;
        RAW[h].is_valid = true;
        RAW[h+1].is_valid = true;
        multis_nb++;
      }
      else if (MAPPING[a][b].multis == 4){
        //
        RAW[h].multis = 4;
        RAW[h].is_valid = true;
        RAW[h+1].is_valid = true;
        RAW[h+2].is_valid = true;
        RAW[h+3].is_valid = true;
        multis_nb++;
      }
    }
  } 

  /////////////////////////////////////////////////////
  // SAMPLE ONLY THE MULTISAMPLED SIGNALS
  /////////////////////////////////////////////////////
  int is_signal[CHANNELS];
  // Variable declaration for feature extraction
  for(int i=0; i<CHANNELS; i++){
    // only look at valid channels
    if (RAW[i].is_valid == false) {
      continue;
    }
    is_signal[i] = 1;    
    /////////////////////////////////////////////////////
    // SIGNAL REGION
    /////////////////////////////////////////////////////
    for(int n=BASELINE_CUT; n<ENERGY_WINDOW_MAX; n++){
      // Subtract the baseline
      // RAW[i].trace[n] -= RAW[i].base.mean;
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
    else { is_signal[i]++; }
  }
  /////////////////////////////////////////////////////
  // TEST FOR COINCIDENCE ONLY IN THE MULTISAMPLED CHANNELS
  /////////////////////////////////////////////////////
  // print waveforms
  if(NOE%1000==0){
    hfile->cd("WAVE_FORMS/RAW");
    plot_waves(RAW, "SIGNAL_RAW", "TRACE");
  }
  // For all events, do the calculation 
  //    An event in a multisampled channel must be in all partaking channels
  for(int i=0; i<CHANNELS; i++){
    // Check if channel takes part in multisampling
    if (RAW[i].multis == 1) continue;
    // if not do calculation
    for(int j=0; j<RAW[i].multis; j++){
      for(int k=0; k<RAW[i].multis; k++){
        if (RAW[i+k].energy == 0.0) MULTIS_NORM[i][j][k].ratio = 0.0;
        else {
          MULTIS_NORM[i][j][k].ratio = RAW[i+j].energy / RAW[i+k].energy;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////
// INTERPOLATE zero crossing for given trace
/////////////////////////////////////////////////////
void interpolate(vector<signal_struct> &signal)
{
  // For time extraction, compare all channels with each other and extract the time information and store 
  //  it in the correct field
  // Outer channel loop: Interpolate each channel trace and calculate the zero crossing point
  int warning_counter = 0;
  for(int i=0; i<(int)signal.size(); i++){
    // If channel is not valid or there's no event, save zero crossing of 0.0
    if (signal[i].is_valid == false || signal[i].is_signal == 0 ){
      for (int j=0; j<(int)signal.size(); j++){
        for (Int_t k = 0; k < N_E_WINDOW; k++){
          signal[i].time[j].timing[k] =  0.0; // When filling the histograms, 0.0 is not entered in the histogram 
        }
      }
      continue;
    } 
    //
    // Else, go on
    // Reset the Xzero marker for both xtals
    signal[i].CFD.Xzero_int=0;
    // Now set the interpolation interval 
    Int_t int_left = signal[i].CFD.Xzero - N_INTPOL_SAMPLES/2;
    Int_t int_right = signal[i].CFD.Xzero + N_INTPOL_SAMPLES/2;
    // Check if interval is out of range
    if (int_left < 0) int_left = 0;
    if (int_right > (int)signal[i].CFD.trace.size()) int_right = (int)signal[i].CFD.trace.size();
    Int_t int_range = int_right - int_left;
    // initialize the new containers for storing the interpolated signal
    Double_t xi[int_range]; 
    Double_t yi[int_range]; 
    signal[i].CFD.x_interpol.clear();
    signal[i].CFD.y_interpol.clear();
    // Initilaize the interpolator
    ROOT::Math::Interpolator inter(int_range, ROOT::Math::Interpolation::kPOLYNOMIAL);
    // Fill arrays with data
    for ( Int_t k = int_left; k < int_right; k++)
    {
      xi[k-int_left]  = (Double_t) k; 
      yi[k-int_left]  = (Double_t) signal[i].CFD.trace[k];
    }
    // Set the Data
    inter.SetData(int_range, xi, yi);
    // printf("%i %i\n", int_left, int_right);
    // Be careful with the range switching from one grid to the other
    for ( Int_t k = 0; k < (Int_t)(NB * int_range - NB + 1); k++ )
    {
      Double_t x_value = (Double_t) int_left + (Double_t) k/NB;
      signal[i].CFD.x_interpol.push_back(x_value);
      // printf("%f\n", tempx[k]);
      signal[i].CFD.y_interpol.push_back( (Double_t) inter.Eval(x_value));
      // Already look for the zero crossing point
      if (signal[i].CFD.y_interpol[k] > 0 && signal[i].CFD.Xzero_int==0){signal[i].CFD.Xzero_int=k;}
    }
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
      signal[i].CFD.x_interpol.erase(signal[i].CFD.x_interpol.begin()); signal[i].CFD.x_interpol.erase(signal[i].CFD.x_interpol.end()-1);
      signal[i].CFD.y_interpol.erase(signal[i].CFD.y_interpol.begin()); signal[i].CFD.y_interpol.erase(signal[i].CFD.y_interpol.end()-1);
    }
    // Use fast linear regression fit to find time value
    if (!linreg(signal[i].CFD.x_interpol,signal[i].CFD.y_interpol,&signal[i].CFD.int_m,&signal[i].CFD.int_b)) {
      warning_counter++;
      if ( warning_counter < 50){
        printf("WARNING (interpolate): Linreg: Singular matrix, can't solve problem\n");
      }
      else{
        printf("WARNING (interpolate): Linreg: Further warnings suppressed.\n");
      }
    }
    // Save the timing information of the zero crossing point
    signal[i].CFD.int_x0 = - (signal[i].CFD.int_b/signal[i].CFD.int_m);
  }
}

// Function for comparing a signal with all other signal and saving the time difference of occurance
void time_compare(vector<signal_struct> &signal){
  // Now compare all the channel
  for (int i = 0; i < (int)signal.size(); i++){
    for (int j = 0; j < (int)signal.size(); j++){
      // Scan through all energy windows and extract timing
      for (Int_t k = 0; k < N_E_WINDOW; k++){
        // Check if channel is valid/is a signal
        if (signal[i].is_valid == false || signal[j].is_valid == false || signal[i].is_signal == 0 || signal[i].is_signal == 0){
          signal[i].time[j].timing[k] = 0.0;
          continue;
        }
        // If all good, proceed
        if (signal[i].energy >= (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH) && 
            signal[i].energy < (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH+E_WINDOW_LENGTH) && 
            signal[j].energy >= (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH) && 
            signal[j].energy < (double) (E_WINDOW_LEFT+k*E_WINDOW_LENGTH+E_WINDOW_LENGTH)){
          // Calculate and norm time difference
          double time_diff = (signal[i].CFD.int_x0 - signal[j].CFD.int_x0) * signal[i].sample_t;
          signal[i].time[j].timing[k] = time_diff;
        }
      }
    }
  }
}



void print_final_statistics(){
  Double_t total_counts=0;
  // 
  //  GENERAL STATISTICS
  //
  for (Int_t i = 0; i<CHANNELS_EFF; i++){
    total_counts += RAW_CALIB[i].h_energy.hist->Integral();
  }
  printf("\n\n\n++++++++     FINAL STATISTICS     ++++++++\n");
  printf("+\n");
  printf("+ TOTAL COUNTS: %i\n", (Int_t) total_counts);
  // Print counts per channel in square form
  for (int a = 0; a<(int)MAPPING.size(); a++){
    printf("-    ");
    for (int b = 0; b<(int)MAPPING[a].size(); b++){
      int ch = b*(int)MAPPING.size()+a; 
      printf("%6d ", (Int_t) RAW_CALIB[ch].h_energy.hist->Integral());
    }
    printf("\n");
  }
  printf("+\n");

  //
  //  COSMICS MODE STATISTICS  
  //
  if(strcmp(MODE, "COSMICS") == 0) {
    Int_t total_coincidents=RAW_CALIB[0].h_energy.hist->Integral();
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
    // Do the langaus fits for the energy histograms
    for (Int_t i = 0; i < (int)RAW_CALIB.size(); i++){
      if (RAW_CALIB[i].is_valid == false) continue;
      // fit the energy histograms
      RAW_CALIB[i].h_energy.params = fit_hist(RAW_CALIB[i].h_energy.hist, RAW_CALIB[i].h_energy.fit, "langaus");
      MA[i].h_energy.params = fit_hist(MA[i].h_energy.hist, MA[i].h_energy.fit, "langaus");
      MWD[i].h_energy.params = fit_hist(MWD[i].h_energy.hist, MWD[i].h_energy.fit, "langaus");
      TMAX[i].h_energy.params = fit_hist(TMAX[i].h_energy.hist, TMAX[i].h_energy.fit, "langaus");
      NMO[i].h_energy.params = fit_hist(NMO[i].h_energy.hist, NMO[i].h_energy.fit, "langaus");
      // Fill the calibration parameters into the histogram
      CALIB.h_RAW_energy.hist->Fill((ENERGY_NORM/RAW_CALIB[i].h_energy.params[8])*CALIB.RAW_energy[i]);
      // Some runtime information
      if (is_in_string(VERBOSE,"p")){
        printf("+ Energy Fit for channel %d done.\n", i);
      }
    }
    // Print Energy extraction 
    if (is_in_string(VERBOSE,"e")){
      // Statistics for the fit
      print_energy_statistics(RAW_CALIB, "RAW_CALIB");
      print_energy_statistics(MA, "MA");
      print_energy_statistics(MWD, "MWD");
      print_energy_statistics(TMAX, "TMAX");
      print_energy_statistics(NMO, "NMO");
      // Resulting calibration parameters
      print_energy_calib();
    }

    //////// TIMING (if there is more than one channel)

    // if (CHANNELS_EFF > 1){
    //   // Start extracting information for all windows
    //   // Do the fitting for the histograms and plotting
    //   plot_timing_hist(RAW_CALIB, "TIMING/RAW_CALIB");
    //   if (is_in_string(VERBOSE,"p")) printf("+ Fitting time hists finished (RAW_CALIB)\n");
    //   plot_timing_hist(MA, "TIMING/MA");
    //   if (is_in_string(VERBOSE,"p")) printf("+ Fitting time hists finished (MA)\n");
    //   plot_timing_hist(MWD, "TIMING/MWD");
    //   if (is_in_string(VERBOSE,"p")) printf("+ Fitting time hists finished (MWD)\n");
    //   plot_timing_hist(TMAX, "TIMING/TMAX");
    //   if (is_in_string(VERBOSE,"p")) printf("+ Fitting time hists finished (TMAX)\n");

    //   // Plot the time vs energy graph
    //   hfile->cd("TIMING/RAW_CALIB");
    //   plot_time_energy(RAW_CALIB[0].time[1]);
    //   hfile->cd("TIMING/MA");
    //   plot_time_energy(MA[0].time[1]); 
    //   hfile->cd("TIMING/MWD");
    //   plot_time_energy(MWD[0].time[1]); 
    //   hfile->cd("TIMING/TMAX");
    //   plot_time_energy(TMAX[0].time[1]); 


    //   // Print out the timing for a filter type
    //   if (is_in_string(VERBOSE,"t")){
    //     print_timing_statistics(RAW_CALIB[0].time[1], total_coincidents, "RAW_CALIB");
    //     print_timing_statistics(MA[0].time[1], total_coincidents, "MA");
    //     print_timing_statistics(MWD[0].time[1], total_coincidents, "MWD");
    //     print_timing_statistics(TMAX[0].time[1], total_coincidents, "TMAX");
    //   }
    // }
  }
  //
  //  PULSER MODE STATISTICS  
  //
  if(strcmp(MODE, "PULSER") == 0){
    for (int i = 0; i < (int)RAW_CALIB.size(); i++){
      fit_hist(RAW_CALIB[i].h_energy.hist, RAW_CALIB[i].h_energy.fit, "multigaus");  
      fit_hist(MA[i].h_energy.hist, MA[i].h_energy.fit, "multigaus");  
      fit_hist(MWD[i].h_energy.hist, MWD[i].h_energy.fit, "multigaus");  
      fit_hist(TMAX[i].h_energy.hist, TMAX[i].h_energy.fit, "multigaus");  
      fit_hist(NMO[i].h_energy.hist, NMO[i].h_energy.fit, "multigaus");  
    }
  }

  //
  //  BEAM MODE STATISTICS
  //
  if(strcmp(MODE, "BEAM") == 0) {
    // Tagger energy statistics if verbose is set
    int sum = 0;
    if (is_in_string(VERBOSE, "g")){
      printf("+++ Tagger Energy statistics +++\n");
      printf("+Tagger #: total counts (multiple couts) (counts w/o multiples and w/ timing cut)\n");
      for (int k = 0; k < TAG_CHANNELS; k++){
        printf("+ Tagger %2d: %d (%d) (%d)\n", k, TAGGER.counts[k], TAGGER.multiples_per_channel[k], (Int_t) RAW_CALIB[CENTRAL].tagged[k].h_energy_mt.hist->Integral());
        sum += TAGGER.counts[k];
      }
      printf("+\n+ Sum of all tagged channels: %d (includes multiples)\n", sum);
      printf("+\n");
      // Tagger multiplicity statistics
      printf("+++ Multiple tagger counts +++\n");
      for (int i = 0; i < TAG_CHANNELS; i++){
        if (TAGGER.multiples_per_count[i] != 0) printf("+ %d events per tag: %2d\n", i, TAGGER.multiples_per_count[i]);
      }
      printf("+\n");
      for (int k = 0; k < TAG_CHANNELS; k++){
        printf("+ TAGGER ENERGY %2d:\n", k);
        for (int a = 0; a<(int)MAPPING.size(); a++){
          printf("-    ");
          for (int b = 0; b<(int)MAPPING[a].size(); b++){
            int ch = b*(int)MAPPING.size()+a; 
            printf("%6d ", (Int_t) RAW_CALIB[ch].tagged[k].h_energy.hist->Integral());
          }
          printf("\n");
        }
      }
      printf("+\n");
      // Tagger timing histogram: with with most entries for cutting decision
      for (int k = 0; k < TAG_CHANNELS; k++){
        printf("TAGGER_CUT%02d=%i#\n", k, largest_1Dbin(TAGGER.t_hist[k].hist,-5000,5000));  
      }
    }
    // Print the split screen tagger energy histograms
    plot_tagger_hist(RAW_CALIB, "ENERGY/TAGGER/PULSE_HIGHT/RAW_CALIB", "pulseheight");
    plot_tagger_hist(MA, "ENERGY/TAGGER/PULSE_HIGHT/MA", "pulseheight");
    plot_tagger_hist(MWD, "ENERGY/TAGGER/PULSE_HIGHT/MWD", "pulseheight");
    plot_tagger_hist(TMAX, "ENERGY/TAGGER/PULSE_HIGHT/TMAX", "pulseheight");
    plot_tagger_hist(NMO, "ENERGY/TAGGER/PULSE_HIGHT/NMO", "pulseheight");
    // Print the split screen tagger energy histograms
    plot_tagger_hist(RAW_CALIB, "ENERGY/TAGGER/INTEGRAL/RAW_CALIB", "integral");
    plot_tagger_hist(MA, "ENERGY/TAGGER/INTEGRAL/MA", "integral");
    plot_tagger_hist(MWD, "ENERGY/TAGGER/INTEGRAL/MWD", "integral");
    plot_tagger_hist(TMAX, "ENERGY/TAGGER/INTEGRAL/TMAX", "integral");
    plot_tagger_hist(NMO, "ENERGY/TAGGER/INTEGRAL/NMO", "integral");

    // Extract TRACELEN 1D hists out of the 2D hist
    char name[100];
    if (EXTRACT_PROTO == 1){
      for (int i = 0; i < (int)RAW_CALIB.size(); i++){
        sprintf(name,"WAVE_FORMS/RAW_CALIB/PROTO_TRACE_%02d", i);
        hfile->cd(name);
        double largest = 0;
        int min = 0;
        int max = 0;
        int nbins = 0;
        for (int n = 1; n <= (int)RAW_CALIB[i].trace.size(); n++){
          sprintf(name, "SLICE%2d", n);
          RAW_CALIB[i].TH1D_proto_trace[n].hist = RAW_CALIB[i].TH2D_proto_trace.hist->ProjectionY(name, n,n);
          RAW_CALIB[i].TH1D_proto_trace[n].hist->Rebin(10);
          // Search for the largest bin
          largest = largest_1Dbin( RAW_CALIB[i].TH1D_proto_trace[n].hist, -5000, 10000);
          // Recalculate bin to normal grid
          min = RAW_CALIB[i].TH1D_proto_trace[n].hist->GetXaxis()->GetXmin();
          max = RAW_CALIB[i].TH1D_proto_trace[n].hist->GetXaxis()->GetXmax();
          nbins = RAW_CALIB[i].TH1D_proto_trace[n].hist->GetSize()-2;
          largest = min + ((max+abs(min))/nbins * largest);
          RAW_CALIB[i].proto_trace_maxbin.push_back(largest);
          // Fit gauss aroundthe largest value
          if ( 110 < n && n < 130) RAW_CALIB[i].TH1D_proto_trace[n].params = fit_hist(RAW_CALIB[i].TH1D_proto_trace[n].hist, RAW_CALIB[i].TH1D_proto_trace[n].fit, "gaus", largest-300,largest+300,0);
          else RAW_CALIB[i].TH1D_proto_trace[n].params = fit_hist(RAW_CALIB[i].TH1D_proto_trace[n].hist, RAW_CALIB[i].TH1D_proto_trace[n].fit, "gaus", largest-300,largest+100,0);
          RAW_CALIB[i].proto_trace_fit.push_back(RAW_CALIB[i].TH1D_proto_trace[n].params[2]);
        }
        plot_trace(RAW_CALIB[i].proto_trace, "PROTO_TRACE", "AL*");
        plot_trace(RAW_CALIB[i].proto_trace_maxbin, "PROTO_TRACE_MAXBIN", "AL*");
        plot_trace(RAW_CALIB[i].proto_trace_fit, "PROTO_TRACE_FIT", "AL*");

        sprintf(name,"WAVE_FORMS/RAW_CALIB/PROTO_TRACE_%02d", i);
        plot_TH2D_hist_graph(RAW_CALIB[i].TH2D_proto_trace.hist, RAW_CALIB[i].proto_trace_maxbin, name, "PROTO_TRACE_2D_MAXBIN");
        plot_TH2D_hist_graph(RAW_CALIB[i].TH2D_proto_trace.hist, RAW_CALIB[i].proto_trace_fit, name, "PROTO_TRACE_2D_FIT");
        // adjust trace for baseline offset
        double mean = array_mean(RAW_CALIB[i].proto_trace_fit, 0, BASELINE_CUT);
        for (int n = 0; n < (int)RAW_CALIB[i].proto_trace_fit.size(); n++ ){
          RAW_CALIB[i].proto_trace_fit[n] -= mean;  
        }
        plot_trace(RAW_CALIB[i].proto_trace_fit, "PROTO_TRACE_FIT_ADJUSTED", "AL*");
      }
      // Extract proto trace to text file for later reading 
      if (proto_out->is_open()){
        for (int n = 0; n < (int)RAW_CALIB[0].proto_trace_fit.size(); n++){
          for (int i = 0; i < (int)RAW_CALIB.size(); i++){
            *proto_out << RAW_CALIB[i].proto_trace_fit[n] << ";"; 
          }
          *proto_out << endl;
        }
        proto_out->close();
      }
      else printf("WARNING(statistics): Unable to write proto trace file.\n");

    }
    if (EXTRACT_PROTO == 0){
      for (int i = 0; i < (int)RAW_CALIB.size(); i++){
        RAW_CALIB[i].proto_trace_fit = PROTO[i].proto_trace_fit;
        plot_trace(PROTO[i].proto_trace_fit, "PROTO_TRACE_FIT", "AL*");
        sprintf(name,"WAVE_FORMS/RAW_CALIB/PROTO_TRACE_%02d", i);
        plot_TH2D_hist_graph(RAW_CALIB[i].TH2D_proto_trace.hist, PROTO[i].proto_trace_fit,name, "PROTO_TRACE_2D_FIT");
      }
    }
    /// Fit for the energy sum
    /// 
/*    for (int i = 0; i < (int)ECAL25.size(); i++){
      for (int k = 0; k < N_E_WINDOW; k++){
        if (k<3) ECAL25[i].tagged[k].h_energy_mt.hist->Rebin(50);
        if (k>=3) ECAL25[i].tagged[k].h_energy_mt.hist->Rebin(50);

        if (k<3) ECAL25[i].tagged[k].h_energy_mt.params = fit_hist(ECAL25[i].tagged[k].h_energy_mt.hist, ECAL25[i].tagged[k].h_energy_mt.fit, "gaus", 0, 5000,  0); 
        if (k==3) ECAL25[i].tagged[k].h_energy_mt.params = fit_hist(ECAL25[i].tagged[k].h_energy_mt.hist, ECAL25[i].tagged[k].h_energy_mt.fit, "gaus", 1000, 8000,  0); 
        if (k==4) ECAL25[i].tagged[k].h_energy_mt.params = fit_hist(ECAL25[i].tagged[k].h_energy_mt.hist, ECAL25[i].tagged[k].h_energy_mt.fit, "gaus", 1000, 8000,  0); 
        if (k>4) ECAL25[i].tagged[k].h_energy_mt.params = fit_hist(ECAL25[i].tagged[k].h_energy_mt.hist, ECAL25[i].tagged[k].h_energy_mt.fit, "gaus", 0, 80000,  0); 
      }
    }
    // Plot the tagger time versus ECAL energy 
    plot_TH2D_hist(TAGGER.h_tagger_vs_energy.hist, "JUNK", "Tagger_time_vs_energy");
    // Plot ECAL energy vs tagger energy
    plot_energy_vs_tagged(ECAL25[3], "JUNK", "energy_vs_tagged");
    // Plot ECAL energy resolution
    plot_energy_resolution(ECAL25[3], "JUNK", "energy_resolution");
*/    // Plot bar chart for the tagger counting frequency
    vector<double> tagger_frequencies;
    for (int k = 0; k<16; k++){
      tagger_frequencies.push_back( TAGGER.counts[k] );
      printf("%d\n", TAGGER.counts[k]);
    }
    hfile->cd("ENERGY/TAGGER/");
    plot_trace(tagger_frequencies, "Tagger Frequency", "AB");
  }

  //

  //
  //  Print the split screen histograms  
  //
  // plot_energy_hist(RAW, "ENERGY/PULSE_HIGHT/RAW", "pulseheight");
  plot_energy_hist(RAW_CALIB, "ENERGY/PULSE_HIGHT/RAW_CALIB", "pulseheight");
  plot_energy_hist(MA, "ENERGY/PULSE_HIGHT/MA", "pulseheight");
  plot_energy_hist(MWD, "ENERGY/PULSE_HIGHT/MWD", "pulseheight");
  plot_energy_hist(TMAX, "ENERGY/PULSE_HIGHT/TMAX", "pulseheight");
  plot_energy_hist(NMO, "ENERGY/PULSE_HIGHT/NMO", "pulseheight");
  //
  // plot_energy_hist(RAW, "ENERGY/INTEGRAL/RAW", "integral");
  plot_energy_hist(RAW_CALIB, "ENERGY/INTEGRAL/RAW_CALIB", "integral");
  plot_energy_hist(MA, "ENERGY/INTEGRAL/MA", "integral");
  plot_energy_hist(MWD, "ENERGY/INTEGRAL/MWD", "integral");
  plot_energy_hist(TMAX, "ENERGY/INTEGRAL/TMAX", "integral");
  plot_energy_hist(NMO, "ENERGY/INTEGRAL/NMO", "integral");
  
  printf("+ END OF STATISTICS +\n");

} 

void print_energy_statistics(vector<signal_struct> &array, const char *name){
  // Print Energy extraction 
  printf("+ COSMIC ENERGY DISTRIBUTIONS FOR %s\n", name);
  for (int i = 0; i < (int)array.size(); i++){
    if (array[i].is_valid == false) continue;
    // Check if histograms/parameters are emtpty are emtpy 
    printf("-    Channel %i:\n", i);
    if (array[i].h_energy.hist->Integral() == 0){
      printf("-       Histogram empty!\n");
      continue;
    } 
    else{
      printf("-       Pos : %4.1f\n", array[i].h_energy.params[8]);
      printf("-       FWHM: %4.1f\n", array[i].h_energy.params[9]);
      printf("-       Calib_param: %4.3f\n", ENERGY_NORM/array[i].h_energy.params[8]);
      printf("-\n"); 
    }
  }
  printf("-\n");
}

void print_energy_calib(){
  for (int i = 0; i < (int)RAW_CALIB.size(); i++){
    if (RAW_CALIB[i].is_valid == false){
      printf("ENERGY_CALIB0%d=1.000,1.000,1.000,1.000 (invalid channel)\n", i);
    }
    else {
      if (i == 0) printf("+ Absolut calibration parameters\n");
      printf("ENERGY_CALIB0%d=%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n", 
      i,
      (ENERGY_NORM/RAW_CALIB[i].h_energy.params[8])*CALIB.RAW_energy[i],
      (ENERGY_NORM/MA[i].h_energy.params[8])*(CALIB.MA_energy[i]*CALIB.RAW_energy[i]),
      (ENERGY_NORM/MWD[i].h_energy.params[8])*CALIB.MWD_energy[i]*CALIB.RAW_energy[i],
      (ENERGY_NORM/TMAX[i].h_energy.params[8]*CALIB.TMAX_energy[i]*CALIB.RAW_energy[i]),
      (ENERGY_NORM/NMO[i].h_energy.params[8]*CALIB.NMO_energy[i]*CALIB.RAW_energy[i])
      );
    }
  }
  printf("+\n");
  for (int i = 0; i < (int)RAW_CALIB.size(); i++){
    if (RAW_CALIB[i].is_valid == false){
      printf("ENERGY_CALIB0%d=1.000,1.000,1.000,1.000 (invalid channel)\n", i);
    }
    else{
      if (i == 0) printf("+ Relative calibration parameters\n");
      printf("ENERGY_CALIB0%d=%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n", 
      i,
      (ENERGY_NORM/RAW_CALIB[i].h_energy.params[8]),
      (ENERGY_NORM/MA[i].h_energy.params[8]),
      (ENERGY_NORM/MWD[i].h_energy.params[8]),
      (ENERGY_NORM/TMAX[i].h_energy.params[8]),
      (ENERGY_NORM/NMO[i].h_energy.params[8])
      );
    }
  }
  printf("+\n");
  for (int i = 0; i < (int)RAW_CALIB.size(); i++){
    if (RAW_CALIB[i].is_valid == false){
      printf("ENERGY_CALIB0%d=1.000,1.000,1.000,1.000 (invalid channel)\n", i);
    }
    else{
      if (i == 0) printf("+ Multiplication calibration factors\n");
      printf("ENERGY_CALIB0%d=%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n", 
      i,
      CALIB.RAW_energy[i],
      CALIB.MA_energy[i],
      CALIB.MWD_energy[i],
      CALIB.TMAX_energy[i],
      CALIB.NMO_energy[i]
      );
    }
  }
  printf("+\n");
}



// Function for printing timing fit information from a time_struct element
void print_timing_statistics(time_struct &array, Int_t total_coincidents, const char *name){
  // Print out the timing for a filter type
  Int_t e_matching_efficiency = 0;
  printf("+ TIMING FOR %s FILTER\n", name);
  for (Int_t k = 0; k < N_E_WINDOW; k++){
    printf("-   TIMING ENERGY WINDOW %i-%i : %3.3f+-%3.3f (%i ENTRIES)\n", 
      (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k, 
      (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
      array.h_timing[k].params[4], 
      array.h_timing[k].params[5], 
      (int) array.h_timing[k].hist->Integral());
    e_matching_efficiency += (Int_t) array.h_timing[k].hist->Integral();
  }
  printf("-\n");
  printf("-   Timing/Energy matching efficiency: %i/%i: %3.1f%%\n", 
        e_matching_efficiency, 
        total_coincidents, 
        (Double_t)e_matching_efficiency/total_coincidents * 100);
  printf("+\n");
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
  for (int i = 0; i<CHANNELS; i++){
    // Only print the multi sampled channels
    if (RAW[i].multis == 1) continue;
    // Else, continue
    // Calculate the real channel from the hardware channel
    //  
    printf("-   Effective Channel %i\n", i);
    vector<vector<vector<double> > > params;
    for (int j = 0; j<RAW[i].multis; j++){
      vector<vector<double> > dim2;
      for (int k = 0; k<RAW[i].multis; k++){
        dim2.push_back(fit_hist(MULTIS_NORM[i][j][k].h_hist.hist, MULTIS_NORM[i][j][k].h_hist.fit, "gaus", 0.,3.));
      }
      params.push_back(dim2);
    }
    // Now print the results
    for (int j = 0; j<RAW[i].multis; j++){
      for (int k = 0; k<RAW[i].multis; k++){
        printf("%f ", params[j][k][2]);
      }
      printf("\n");
    }
    printf("\n");
  }

  printf("+\n");
  printf("+++++++++++++++++++++++++++++++++++++++++\n");
}

void plot_waves(vector<signal_struct> &array, char const *name, char const *modus) {
  // Canvas for the combined waveforms
  TCanvas *c_combined = new TCanvas("c_combined","Wave_forms",200,10,1280,1024);
  // Canvas for split wavedforms
  TCanvas *c_split = new TCanvas("c_split","Wave_forms",10,10,700,900);
  // Check if array is a RAW trace
  if (array[0].is_raw) c_split->Divide(BOARDS, 8);
  else { c_split->Divide((int)MAPPING[0].size(), (int)MAPPING.size()); }
  // Set grid
  c_combined->SetGrid();
  TMultiGraph *mg_combined = new TMultiGraph();
  mg_combined->SetTitle("Signal example; Time [1ns]; ADC channel [arb. unit]");
  TGraph *tg_combined[array.size()];
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader("ADC Digitization"); // option "C" allows to center the header
  for(Int_t i=0; i<(int)array.size(); i++){   //Channel loop 
    // Only paint vaild channels / non-empty channels
    if (array[i].is_valid == false || (int)array[i].trace.size() == 0) continue;
    // If valid, paint
    if (strcmp(modus, "TRACE") == 0){
      Double_t wave_y[array[i].trace.size()];
      Double_t wave_x[array[i].trace.size()];
      for(Int_t n = 0; n< (Int_t) array[i].trace.size(); n++){
        wave_y[n] = array[i].trace[n]; 
        wave_x[n] = n * array[i].sample_t; // Calibrate to the sampling rate
      }
      tg_combined[i] = new TGraph((Int_t) array[i].trace.size(),wave_x,wave_y);
    }
    else if (strcmp(modus, "CFD") == 0){
      Double_t wave_y[array[i].CFD.trace.size()];
      Double_t wave_x[array[i].CFD.trace.size()];
      for(Int_t n = 0; n< (Int_t) array[i].CFD.trace.size(); n++){
        wave_y[n] = array[i].CFD.trace[n];
        wave_x[n] = n * array[i].sample_t; // Calibrate to the sampling rate
      }
      tg_combined[i] = new TGraph((Int_t) array[i].CFD.trace.size(),wave_x,wave_y);
    } 
    else{ printf("ERROR (plot_waves): Plot option invalid\n");}
    // tg_combined[i]->SetLineColor(i+1);
    tg_combined[i]->SetLineColor(1);
    // tg_combined[i]->SetMarkerColor(i+1);
    tg_combined[i]->SetMarkerColor(1);
    tg_combined[i]->SetMarkerSize(0.35);
    tg_combined[i]->SetLineWidth(2);
    //tg_waves[i-1]->Draw("AL*");
    char c_number[12];
    char text[12] = "Channel_";
    sprintf(c_number,"%i",i);
    strcat(text,c_number);
    tg_combined[i]->SetTitle(text);
    legend->AddEntry(tg_combined[i],text,"f");
    mg_combined->Add(tg_combined[i]);
    // Now draw the single wave into the split screen view
    // Check if RAW container
    int ch;
    if (array[i].is_raw) {
      ch = (i / 8) + (i%8)*BOARDS;
    } 
    else {ch = (i / (int)MAPPING.size()) + (i%(int)MAPPING.size())*(int)MAPPING[0].size(); }
    c_split->cd(ch+1);
    gPad->SetGrid();
    gPad->SetBorderSize(5);
    // For the splitscreen, only draw black waves
    tg_combined[i]->Draw("AL*");
    // Reset to the normal color scheme
  }
  c_combined->cd();
  mg_combined->Draw("AL*");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_combined->Write(name);
  char name_split[100];
  sprintf(name_split, "%s_split", name);
  c_split->Write(name_split);
  // delete legend;
  // delete[] tg_combined;
  // delete mg_combined;
  delete c_combined;
  delete c_split;
}

void plot_waves_compare(char const *name, char const *modus) {
  // Canvas for signals split into channels
  TCanvas *c_split = new TCanvas("c_split","Wave_forms",10,10,1280,1024);
  c_split->Divide((int)MAPPING[0].size(), (int)MAPPING.size());
  // Set grid
  c_split->SetGrid();
  // Now loop through all channel
  for(Int_t i=0; i < (int)RAW_CALIB.size(); i++){   //Channel loop 
    int tracelen = (int)RAW_CALIB[i].trace.size();
    // Only paint vaild channels / non-empty channels
    if (RAW_CALIB[i].is_valid == false || tracelen == 0) continue;
      // Create a multigraph canvas
    TMultiGraph *mg_combined = new TMultiGraph();
    // Set its title
    mg_combined->SetTitle("Signal example; Time [1ns]; ADC channel [arb. unit]");
    // Create legend
    TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetHeader("ADC Digitization"); // option "C" allows to center the header
    // Create different graphs
    TGraph *tg_combined[5]; // For the 4 Filter types 
    // If valid, paint
    if (strcmp(modus, "TRACE") == 0){
      // Container for passing to the Root TGraph because vectors don't work
      Double_t wave_y[4][tracelen];
      Double_t wave_x[4][tracelen];
      Double_t wave_y_long[(int)NMO[i].trace.size()];
      Double_t wave_x_long[(int)NMO[i].trace.size()];
      //
      for(Int_t n = 0; n < tracelen; n++){
        wave_y[0][n] = RAW_CALIB[i].trace[n]; 
        wave_x[0][n] = n * RAW_CALIB[i].sample_t; // Calibrate to the sampling rate
        wave_y[1][n] = MA[i].trace[n]; 
        wave_x[1][n] = n * MA[i].sample_t; // Calibrate to the sampling rate
        wave_y[2][n] = MWD[i].trace[n]; 
        wave_x[2][n] = n * MWD[i].sample_t; // Calibrate to the sampling rate
        wave_y[3][n] = TMAX[i].trace[n]; 
        wave_x[3][n] = n * TMAX[i].sample_t; // Calibrate to the sampling rate
      }
      for(Int_t n = 0; n < (int)NMO[i].trace.size(); n++){
        wave_y_long[n] = NMO[i].trace[n]; 
        wave_x_long[n] = n * NMO[i].sample_t; // Calibrate to the sampling rate
        // printf("%3.2f %3.2f %3.2f\n", wave_x_long[n], wave_y_long[n], NMO[i].sample_t );
      }

      for (int j = 0; j < 4; j++){
        if ( j == 1 || j == 2) continue; // dont print MA, MWD
        tg_combined[j] = new TGraph(tracelen,wave_x[j],wave_y[j]);
      }
      tg_combined[4] = new TGraph((int)NMO[i].trace.size(),wave_x_long,wave_y_long);
    }
    else if (strcmp(modus, "CFD") == 0){
      Double_t wave_y[4][tracelen];
      Double_t wave_x[4][tracelen];
      Double_t wave_y_long[(int)NMO[i].trace.size()];
      Double_t wave_x_long[(int)NMO[i].trace.size()];
      //
      for(Int_t n = 0; n< (Int_t) RAW_CALIB[i].CFD.trace.size(); n++){
        wave_y[0][n] = RAW_CALIB[i].CFD.trace[n];
        wave_x[0][n] = n * RAW_CALIB[i].sample_t; // Calibrate to the sampling rate
        wave_y[1][n] = MA[i].CFD.trace[n];
        wave_x[1][n] = n * MA[i].sample_t; // Calibrate to the sampling rate
        wave_y[2][n] = MWD[i].CFD.trace[n];
        wave_x[2][n] = n * MWD[i].sample_t; // Calibrate to the sampling rate
        wave_y[3][n] = TMAX[i].CFD.trace[n];
        wave_x[3][n] = n * TMAX[i].sample_t; // Calibrate to the sampling rate
      }
      for(Int_t n = 0; n < (int)NMO[i].trace.size(); n++){
        wave_y_long[n] = NMO[i].CFD.trace[n];
        wave_x_long[n] = n * NMO[i].sample_t; // Calibrate to the sampling rate
      }

      for (int j = 0; j < 4; j++){
        if ( j == 1 || j == 2) continue; // dont print MA, MWD
        tg_combined[j] = new TGraph(tracelen,wave_x[j],wave_y[j]);
      }
      tg_combined[4] = new TGraph((int)NMO[i].trace.size(),wave_x_long,wave_y_long);
    } 
    else{ printf("ERROR (plot_waves): Plot option invalid\n");}
    int color = 0;
    for (int j = 0; j < 5; j++){
    	if (j == 1 || j == 2) continue;
      if (j==0) color = 1;
      if (j==3) color = 4;
      if (j==4) color = 2;
      tg_combined[j]->SetLineColor(color);
      tg_combined[j]->SetMarkerColor(color);
      tg_combined[j]->SetMarkerSize(0.35);
      tg_combined[j]->SetLineWidth(2);
      char text[100];
      sprintf(text,"SIGNAL_%i",j);
      tg_combined[j]->SetTitle(text);
      legend->AddEntry(tg_combined[j],text,"f");
      mg_combined->Add(tg_combined[j]);
    }
    int ch = (i / (int)MAPPING.size()) + (i%(int)MAPPING.size())*(int)MAPPING[0].size();
    c_split->cd(ch+1);
    gPad->SetGrid();
    // For the splitscreen, only draw black waves
    mg_combined->Draw("AL*");
    legend->Draw();
    // delete legend;
    // delete[] tg_combined;
    // delete mg_combined;
  }
  // gPad->Modified();
  // gPad->Update();
  char name_split[100];
  sprintf(name_split, "%s_split", name);
  c_split->Write(name_split);
  delete c_split;
}


void plot_time_energy(time_struct &array) {
  TCanvas *c_waves = new TCanvas("c1","Wave_forms",200,10,500,300);
  c_waves->SetGrid();
  TGraphErrors *tg_waves;
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader(""); // option "C" allows to center the header
  // Calculate the energy steps
  double e_step = (E_WINDOW_RIGHT - E_WINDOW_LEFT)/N_E_WINDOW/100; // channels calculated to 100ch=1MeV
  Double_t wave_x[N_E_WINDOW];
  for(Int_t k = 0; k<N_E_WINDOW; k++){
    wave_x[k] = e_step + k*e_step;
  }
  Double_t wave_y[N_E_WINDOW];
  for(Int_t k = 0; k<N_E_WINDOW; k++){
    // Check if timing array is empty
    if ((int)array.h_timing[k].params.size() < 1){
      printf("WARNING (plot_time_energy): Timing array is empty. Maybe there is only one channel in the readout?\n");
      return;
    }
    wave_y[k] = array.h_timing[k].params[4];
  }
  Double_t wave_yerr[N_E_WINDOW];
  for(Int_t k = 0; k<N_E_WINDOW; k++){
    wave_yerr[k] = array.h_timing[k].params[5];
  }
  tg_waves = new TGraphErrors(N_E_WINDOW,wave_x,wave_y,0,wave_yerr);
  tg_waves->GetXaxis()->SetTitle("Energy [MeV]");
  tg_waves->GetYaxis()->SetTitle("Time resolution [ns]");
  tg_waves->SetLineColor(1);
  tg_waves->SetMarkerColor(2);
  tg_waves->SetLineWidth(2);
  tg_waves->SetTitle("Time resolution for different energy ranges");  
  legend->AddEntry(tg_waves,"Ch1/Ch2","f");
  tg_waves->Draw("AL*");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_waves->Write("Signal_RAW");
  delete c_waves;
}

void plot_interpol(vector<double> &x, vector<double> &y){
  TCanvas *c_waves = new TCanvas("c1","Wave_interpolation",200,10,500,300);
  c_waves->SetGrid();
  TMultiGraph *mg_waves = new TMultiGraph();
  mg_waves->SetTitle("Interpolation exaple; Time [ns]; interpol. ADC channel [arb. unit]");
  TGraph *tg_waves;
  TGraph *tg_fits;
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  int warning_counter = 0;
  legend->SetHeader("Interpolated ADC values"); // option "C" allows to center the header

  // Fit the passed waves with a linear regression
  double m = 0; double b = 0;
  if (!linreg(x,y,&m,&b)) {
    warning_counter++;
    if (warning_counter < 50){
      printf("WARNING (plot_interpol): Linreg: Singular matrix, can't solve problem\n");
    }
    else{
      printf("WARNING (plot_interpol): Linreg: Further warnigs supressed.\n");
    }
  }
  // Transform vector to array, for handing it to incompetent root TGraph
  Double_t wave_x[x.size()];
  Double_t wave_y[y.size()];
  Double_t fit_y[y.size()];
  for(Int_t n = 0; n<(Int_t)x.size(); n++){
    wave_x[n] = (Double_t) x[n] * SAMPLE_t_eff; // Calibrate to the sampling rate
    wave_y[n] = (Double_t) y[n];
    fit_y[n] = (Double_t) m * x[n] + b; // linear fit for wave
  }
  // Build TGraph
  tg_waves = new TGraph(x.size(),wave_x,wave_y);
  tg_waves->SetLineColor(1);
  tg_waves->SetMarkerColor(1);
  tg_waves->SetLineWidth(2);
  tg_waves->SetTitle("Channel i");
  tg_waves->SetMarkerStyle(kOpenSquare); // Asterisk
  tg_fits = new TGraph(x.size(),wave_x,fit_y);
  tg_fits->SetLineColor(1);
  tg_fits->SetMarkerColor(1);
  tg_fits->SetLineWidth(2);
  tg_fits->SetTitle("Channel i");
  tg_fits->SetMarkerStyle(kDot);
  //tg_waves[i-1]->Draw("AL*");
  char waves_number[25];
  char fit_number[25];
  char waves_text[25] = "Channel_";
  char fit_text[25] = "Channel_";
  sprintf(waves_number,"1");
  sprintf(fit_number,"1_fit");
  strcat(waves_text,waves_number);
  strcat(fit_text,fit_number);
  legend->AddEntry(tg_waves,waves_text,"f");
  legend->AddEntry(tg_fits,fit_text,"f");
  mg_waves->Add(tg_waves);
  mg_waves->Add(tg_fits);

  mg_waves->Draw("AL");
  legend->Draw();
  gPad->Modified();
  gPad->Update();
  c_waves->Write("Signal_Interpolated"); 
  delete c_waves;
}

void plot_multis_hist(){
  hfile->cd("ENERGY/PULSE_HIGHT/RAW/Intersampling_calibration");
  // THStack *hs[CHANNELS_EFF];
  for (int i=0; i<CHANNELS; i++){
    // Only plot valid channels
    if (RAW[i].multis == 1) continue;
    // THStack hs[i] = new THStack("hs","");
    TCanvas *cs = new TCanvas("cs","cs",10,10,700,900);
    int counter = 1;
    cs->Divide(4,4);
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        cs->cd(counter); 
        MULTIS_NORM[i][j][k].h_hist.hist->Draw();
        counter++;
      }
    }
    char name[100];
    sprintf(name, "EFF_CHANNEL%i", i);
    cs->Write(name);
    delete cs;
  }
}

void plot_energy_hist(vector<signal_struct> &array, char const *path, const char *mode){
  // GENERAL Energy Histograms
  int ch;
  hfile->cd(path);
  TCanvas *cs = new TCanvas("cs","cs",10,10,1280,1024);
  // Check if array is a RAW trace
  if (array[0].is_raw) cs->Divide(BOARDS, 8);
  else { cs->Divide((int)MAPPING[0].size(), (int)MAPPING.size()); }
  // THStack *hs[CHANNELS_EFF];
  for (int i=0; i<(int)array.size(); i++){
    // Only plot valid channels
    if (array[i].is_valid == false) continue;
    // Check if RAW container
    if (array[i].is_raw) {
      ch = (i / 8) + (i%8)*BOARDS;
    }
    else {ch = (i / (int)MAPPING.size()) + (i%(int)MAPPING.size())*(int)MAPPING[0].size(); }
    // Choose the right field for the canvas
    cs->cd(ch+1);
    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptFit(1111);
    gStyle->SetLabelSize(0.05,"X");
    if (strcmp(mode,"pulseheight")==0){
      array[i].h_energy.hist->Rebin(1);
      array[i].h_energy.hist->SetFillColor(kBlue-7);
      array[i].h_energy.hist->SetFillStyle(4050);
      array[i].h_energy.hist->GetXaxis()->SetNdivisions(5);
      array[i].h_energy.hist->Draw();
    }
    if (strcmp(mode,"integral")==0){
      array[i].h_integral.hist->Rebin(1);
      array[i].h_integral.hist->SetFillColor(kBlue-7);
      array[i].h_integral.hist->SetFillStyle(4050);
      array[i].h_integral.hist->GetXaxis()->SetNdivisions(5);
      array[i].h_integral.hist->Draw();
    }
  }
  if (strcmp(mode,"pulseheight")==0){ cs->Write("ENERGY_SPLIT_SCREEN"); }
  if (strcmp(mode,"integral")==0){ cs->Write("INTEGRAL_SPLIT_SCREEN"); }
  delete cs;
}

void plot_TH2D_hist(TH2D *hist, char const *path, const char *name){
  // GENERAL Energy Histograms
  hfile->cd(path);
  TCanvas *cs = new TCanvas("cs","cs",10,10,1280,1024);
  // Check if array is a RAW trace
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptFit(1111);
  gStyle->SetLabelSize(0.05,"X");
  hist->Draw("COLZ");

  cs->Write(name); 
  delete cs;
}

void plot_TH2D_hist_graph(TH2D *hist, vector<double> trace, char const *path, const char *name){
  // GENERAL Energy Histograms
  hfile->cd(path);
  TCanvas *cs = new TCanvas("cs","cs",10,10,1280,1024);
  // Check if array is a RAW trace
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptFit(1111);
  gStyle->SetLabelSize(0.05,"X");
  hist->Draw("COLZ");
  // Now the graph part
  TGraph *graph;
  Double_t wave_y[trace.size()];
  Double_t wave_x[trace.size()];
  for(Int_t n = 0; n< (Int_t)trace.size(); n++){
    wave_y[n] = trace[n]; 
    wave_x[n] = n; // Calibrate to the sampling rate
  }

  graph = new TGraph((Int_t)trace.size(),wave_x,wave_y);
  
  graph->SetLineColor(kRed);
  // graph->SetMarkerColor(i+1);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerSize(0.35);
  graph->SetLineWidth(2);

  graph->Draw("SAME");


  cs->Write(name); 
  delete cs;
}


void plot_energy_vs_tagged(signal_struct &signal, char const *path, const char *name){
  // GENERAL Energy Histograms
  hfile->cd(path);
  TCanvas *cs = new TCanvas("cs","cs",10,10,1280,1024);
  // Check if array is a RAW trace
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptFit(1111);
  gStyle->SetLabelSize(0.05,"X");
  // Now the graph part
  TGraphErrors *graph;
  TGraph *graph2;
  Double_t wave_y[(int)signal.tagged.size()];
  Double_t wave_y_err[(int)signal.tagged.size()];
  Double_t wave_x[(int)signal.tagged.size()];

  for(Int_t n = 0; n < (int)signal.tagged.size(); n++){
    wave_y[n] = signal.tagged[n].h_energy_mt.params[2] / 100; //  Norm to MeV 
    wave_y_err[n] = signal.tagged[n].h_energy_mt.params[4] /100 ; // Norm to MeV
    wave_x[n] = TAGGER.energy[n]; // Tagger energies
  }

  graph = new TGraphErrors((int)signal.tagged.size(),wave_x,wave_y,0,wave_y_err);

  graph->RemovePoint(5);
  
  graph->SetLineColor(kBlack);
  // graph->SetMarkerColor(i+1);
  graph->SetMarkerColor(kBlack);
  graph->SetMarkerStyle(24);
  graph->SetMarkerSize(2);
  graph->SetLineWidth(2);
  graph->GetXaxis()->SetTitle("Tagger Energy [MeV]");
  graph->GetYaxis()->SetTitle("Detector Energy [MeV]");

  graph->Draw("ap");

  Double_t wave_y2[800];
  Double_t wave_x2[800];

  for(Int_t n = 0; n < 800; n++){
    wave_y2[n] = n;
    wave_x2[n] = n;
  }

  graph2 = new TGraph(800,wave_x2,wave_y2);

  graph2->SetMarkerColor(kGray+2);
  graph2->SetMarkerStyle(31);
  graph2->SetMarkerSize(0.35);
  graph2->SetLineWidth(2);
  graph2->SetLineStyle(9);

  graph2->Draw("SAME");

  vector<Double_t> fit_params;

  fit_params = fit_graph_err(graph, "SIPMpixel", 0, 500, 0);


  // printf("+++++++++ FIT +++++++++++\n");
  // for (int i = 0; i < (int)fit_params.size(); i++){
  //   printf("%3.3f\n", fit_params[i]);
  // }



  cs->Write(name); 
  delete cs;
}

// Energy resolution
void plot_energy_resolution(signal_struct &signal, char const *path, const char *name){
  // GENERAL Energy Histograms
  hfile->cd(path);
  TCanvas *cs = new TCanvas("cs","cs",10,10,1280,1024);
  // Check if array is a RAW trace
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptFit(1111);
  gStyle->SetLabelSize(0.05,"X");
  // Now the graph part
  TGraphErrors *graph;
  Double_t wave_y[(int)signal.tagged.size()];
  Double_t wave_y_err[(int)signal.tagged.size()];
  Double_t wave_x[(int)signal.tagged.size()];

  for(Int_t n = 0; n < (int)signal.tagged.size(); n++){
    // From gaus fit:
    // Field 0,1: amp, amp_err
    // Field 2,3: mean, mean_err
    // Field 4,5: std, std_err 
    wave_y[n] = signal.tagged[n].h_energy_mt.params[4] / (TAGGER.energy[n]*100) *100; //  Norm to MeV and to percent
    wave_y_err[n] = signal.tagged[n].h_energy_mt.params[5] / (TAGGER.energy[n]*100 ) *100 ; // Norm to MeV and to percent
    wave_x[n] = TAGGER.energy[n]; // Tagger energies
  }

  graph = new TGraphErrors((int)signal.tagged.size(),wave_x,wave_y,0,wave_y_err);

  graph->RemovePoint(5); // Tagger energy does not exist
  graph->SetName("data_graph");
  
  graph->SetLineColor(kBlack);
  graph->SetTitle("Data");
  // graph->SetMarkerColor(i+1);
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(23);
  graph->SetMarkerSize(4);
  graph->SetLineWidth(2);
  graph->GetXaxis()->SetTitle("Tagger Energy [MeV]");
  graph->GetYaxis()->SetTitle("\\frac{\\sigma}{E} \\left[\\%\\right]");

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);



  graph->Draw("ap");

  vector<Double_t> fit_params;

  fit_params = fit_graph_err(graph, "resolution", 10, 300, 0);

  graph->GetFunction("resolution")->SetLineColor(kRed);
  graph->GetFunction("resolution")->SetLineWidth(4);

  cs->BuildLegend(0.15, 0.7, 0.4, 0.9,"");

  // TPaveStats st = ((TPaveStats)(graph->GetListOfFunctions()->FindObject("resolution")));
  // if (st) {
  //   st->SetTextColor(graph->GetLineColor());
  //   st->SetX1NDC(0.64); st->SetX2NDC(0.99);
  //   st->SetY1NDC(0.4); st->SetY2NDC(0.6);
  // }

  cs->Modified(); cs->Update(); // make sure it’s really (re)drawn

  cs->Write(name); 
  delete cs;
}


void plot_tagger_hist(vector<signal_struct> &array, char const *path, const char *mode){
  // GENERAL Energy Histograms
  int ch;
  hfile->cd(path);
  for (int k = 0; k < TAG_CHANNELS; k++){
    TCanvas *cs = new TCanvas("cs","cs",10,10,700,900);
    // Check if array is a RAW trace
    if (array[0].is_raw) cs->Divide(BOARDS, 8);
    else { cs->Divide((int)MAPPING[0].size(), (int)MAPPING.size()); }
    // THStack *hs[CHANNELS_EFF];
    for (int i=0; i<(int)array.size(); i++){
      // Only plot valid channels
      if (array[i].is_valid == false) continue;
      // Check if RAW container
      if (array[i].is_raw) {
        ch = (i / 8) + (i%8)*BOARDS;
      }
      else {ch = (i / (int)MAPPING.size()) + (i%(int)MAPPING.size())*(int)MAPPING[0].size(); }
      // Choose the right field for the canvas
      cs->cd(ch+1);

      gPad->SetGridx();
      gPad->SetGridy();
      if (strcmp(mode,"pulseheight")==0){
        array[i].tagged[k].h_energy.hist->Rebin(20);
        array[i].tagged[k].h_energy.hist->SetLineColor(1);
        // array[i].tagged[k].h_energy.hist->SetFillColorAlpha(5, 0.35);
        // array[i].tagged[k].h_energy.hist->SetFillStyle(4050);
        array[i].tagged[k].h_energy_m.hist->Rebin(20);
        array[i].tagged[k].h_energy_m.hist->SetLineColor(2);
        // array[i].tagged[k].h_energy_m.hist->SetFillColorAlpha(8, 0.35);
        // array[i].tagged[k].h_energy_m.hist->SetFillStyle(4050);
        array[i].tagged[k].h_energy_mt.hist->Rebin(20);
        array[i].tagged[k].h_energy_mt.hist->SetLineColor(4);
        // array[i].tagged[k].h_energy_mt.hist->SetFillColorAlpha(10, 0.35);
        // array[i].tagged[k].h_energy_mt.hist->SetFillStyle(4050);
        array[i].tagged[k].h_energy.hist->Draw();
        array[i].tagged[k].h_energy_m.hist->Draw("same");
        array[i].tagged[k].h_energy_mt.hist->Draw("same");
      }
      if (strcmp(mode,"integral")==0){
        array[i].tagged[k].h_integral.hist->Rebin(20);
        array[i].tagged[k].h_integral.hist->SetLineColor(1);
        // array[i].tagged[k].h_integral.hist->SetFillColorAlpha(5, 0.35);
        // array[i].tagged[k].h_integral.hist->SetFillStyle(4050);
        array[i].tagged[k].h_integral_m.hist->Rebin(20);
        array[i].tagged[k].h_integral_m.hist->SetLineColor(2);
        // array[i].tagged[k].h_integral_m.hist->SetFillColorAlpha(8, 0.35);
        // array[i].tagged[k].h_integral_m.hist->SetFillStyle(4050);
        array[i].tagged[k].h_integral_mt.hist->Rebin(20);
        array[i].tagged[k].h_integral_mt.hist->SetLineColor(4);
        // array[i].tagged[k].h_integral_mt.hist->SetFillColorAlpha(10, 0.35);
        // array[i].tagged[k].h_integral_mt.hist->SetFillStyle(4050);
        array[i].tagged[k].h_integral.hist->Draw();
        array[i].tagged[k].h_integral_m.hist->Draw("same");
        array[i].tagged[k].h_integral_mt.hist->Draw("same");
      }

    }
    char name[200];
    if (strcmp(mode,"integral")==0){ sprintf(name, "TAGGER_%i_INTEGRAL", k); }
    if (strcmp(mode,"pulseheight")==0){ sprintf(name, "TAGGER_%i_ENERGY", k); }
    cs->Write(name);
    delete cs;
  }
}





void plot_timing_hist(vector<signal_struct> &signal, char const *path){
  // Do the gaus fitting of the histograms
  for (int i=0; i<(int)signal.size(); i++){
    for (int j=0; j<(int)signal.size(); j++){
      for (Int_t k = 0; k < N_E_WINDOW; k++){
        // Only fit valid channels
        if (signal[i].is_valid == false || signal[j].is_valid == false) continue;
        signal[i].time[j].h_timing[k].params = fit_hist(signal[i].time[j].h_timing[k].hist, signal[i].time[j].h_timing[k].fit, "gaus", -5.,5.);
      }
    }
  }
  //
  for (int i=0; i<(int)signal.size(); i++){
    for (Int_t k = 0; k < N_E_WINDOW; k++){
      // For each channel and energy, create one split screen
      int ch;
      char name[100];
      sprintf(name,"%s/CHANNEL_%i", path, i);
      hfile->cd(name);
      TCanvas *cs = new TCanvas("cs","cs",10,10,700,900);
      cs->Divide((int)MAPPING[0].size(), (int)MAPPING.size());
      for (int j=0; j<(int)signal.size(); j++){
        // Only plot valid channels
        if (signal[j].is_valid == false) continue;
        // Check if RAW container
        ch = (j / (int)MAPPING.size()) + (j%(int)MAPPING.size())*(int)MAPPING[0].size(); 
        // Choose the right field for the canvas
        cs->cd(ch+1);
        gPad->SetGridx();
        gPad->SetGridy();
        signal[i].time[j].h_timing[k].hist->Draw();
      }
      sprintf(name,"CHANNEL_%i_ENERGY_%i", i, (int)(E_WINDOW_RIGHT-E_WINDOW_LEFT)/N_E_WINDOW*k);
      cs->Write(name);
      delete cs;
    }
  }
}

void plot_trace(vector<double> &trace, char const *name, char const *modus) {
  // Canvas for the combined waveforms
  TCanvas *canvas = new TCanvas("canvas","Wave_forms",200,10,1280,1024);
  canvas->SetGrid();
  canvas->SetTitle("Signal example; Time [1ns]; ADC channel [arb. unit]");

  TGraph *graph;

  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader("ADC Digitization"); // option "C" allows to center the header

  Double_t wave_y[trace.size()];
  Double_t wave_x[trace.size()];
  for(Int_t n = 0; n< (Int_t)trace.size(); n++){
    wave_y[n] = trace[n]; 
    wave_x[n] = n; // Calibrate to the sampling rate
  }

  graph = new TGraph((Int_t)trace.size(),wave_x,wave_y);
  
  graph->SetLineColor(1);
  // graph->SetMarkerColor(i+1);
  graph->SetMarkerColor(1);
  graph->SetMarkerSize(0.35);
  graph->SetLineWidth(2);

  if (is_in_string(modus, "B")) graph->SetFillColor(38);

  graph->Draw(modus);

  canvas->Write(name);
  // delete legend;
  // delete[] tg_combined;
  // delete mg_combined;
  delete canvas;
}


double randit(int ini){
  if(ini==1) srand(time(NULL));
  return (((rand()%100) -50.) /100.);
}

// Calculates the mean of all elements between start and end of an array
double array_mean(const vector<double> &array, int start, int end){
  double sum = 0;
  for (int i = start; i<end; i++){
    sum += array[i];
  }
  return (double) sum/(end-start);
}

// Calculates the std of all elements between start and end of an array
double array_std(const vector<double> &array, int start, int end, double mean){
  double sum = 0;
  for (int i = start; i<end; i++){
    sum += (double) pow(array[i]-mean, 2.);
  }
  return (double) TMath::Sqrt(sum/(end-start));
}

// Return index of the largest value in array
int array_largest(vector<double> &array, int lower, int upper){
	//
	int index = 0;
	int max = 0;
	// Check boundaries
	if (lower < 0) lower = 0;
	if (upper > (int)array.size()) upper = (int)array.size(); 
	// Loop over all indices
	for (int n = lower; n < upper; n++){
		//
		if ( array[n] > max ) {
			max = array[n];
			index = n;
		}
	}
	// 
	return(index);

}// Return an array of the sum of to arrays. each field is summed individually and can be automatically divided by a factor
vector<double> array_sum(vector<double> &array1, vector<double> &array2, double factor = 1){
  //
  vector<double> returned; 
  // Check boundaries if arrays have the same length
  if ( array1.size() != array2.size() ){
    printf("WARNING(array_sum): Arrays do not have the same lentgh! Nul array of lentgh 0 returned!\n");
    return returned;
  }
  // If ok, do the sum
  for (int i = 0; i < (int)array1.size(); i++){
    returned.push_back( (array1[i]+array2[i]) * factor );
  }
  // 
  return(returned);
}

// Looks for the first zero crossing between lower and upper of an array, returns index right to the xing
int array_zero_xing(vector<double> &array, int lower, int upper){
	// Check boundaries
	if (lower < 0) lower = 0;
	if (upper > (int)array.size()) upper = (int)array.size(); 
	// Loop over all indices
	for (int n = lower; n < upper; n++){
		//
    if ( array[n-1] < 0 && array[n] > 0 ){
    	return(n);
    }
	}
	// if no xing is found, return 0
	return(0);
}

// Returns the sum of the weighted differences of each element of array1 and array2 between index start and index end
double array_compare(vector<double> &array1, vector<double> &array2, vector<double> &weigths, vector<double> &par, int start, int end, int debug){
  // Some debug output
  debug=0;
  if (debug == 1){
    printf("DEBUG(array_compare): %d %d %d\n", (int)array1.size(), (int)array2.size(), (int)weigths.size());
  }
  // Conversion factor to compare two non equal arrays (who have the same time frame)
  double conversion = 1.;
  double size1 = (double)array1.size();
  double size2 = (double)array2.size();
  // Check all arrays have the same size
  if ( array1.size() != array2.size() ){
    conversion = size1 / size2 ;
    if (debug == 1) printf("WARNING(array_compare): array1 (%d) and array2 (%d) of different length! Using conversion factor %3.2f to compare arrays.\n", (int)array1.size(), (int)array2.size(), conversion);
  }
  // Check if the start and stop values are valid 
  if ( start < 0 ) start = 0;
  if ( start > (int)array1.size() ) start = 0;
  if ( end > (int)array1.size() ) end = (int)array1.size();
  // if no weigths are passed, initialize the weigths as 1 for every sample
  if ( weigths.size() == 0 ){
    for (int n = 0; n < (int)array1.size(); n++) weigths.push_back(1.0);
  }
  if ( array1.size() != weigths.size() ){
    printf("WARNING(array_compare): array1 and weigths of different length! returned value of 1000000.\n");
    return(1000000.0);
  }
  // Check position adjustment parameter is an "integer"
  if ( floor(par[1]) != ceil(par[1]) ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed position parameter par[1] is not an integer. Rounded to %d.\n", (int)round(par[1]));
    par[1] = (double)round(par[1]);
  }

  // 
  double sum = 0.0;
  int n2 = 0; 
  // No everything should be of the correct size
  for (int n1 = start; n1 < end; n1++){
    // Using array1 as reference frame, calculate equivalent sample value for array2
    n2 = (int)( (int)(round(n1/conversion)) - (int)par[1] );
    // 
    if (n2 > (int)array2.size()) printf("WARNING(array_compare): Conversion went wrong, calculated sample n2=%d is larger than array2 size %d!\n", n2, (int)array2.size());
    // Do the sum
    sum += weigths[n1] * abs( array1[n1]- par[0] * array2[n2] );
  }
  return(sum/(end-start));
}

/* Adjust vector according to adjustment parameters vector<double> x:
    - x[0]: A is multiplication factor for amplitude adjustment
    - x[1]: p is position adjustment, should be integer. Shifts array by p positions towards higher or lower indices
*/
vector<double> array_adjust(vector<double> &array, vector<double> &x, int debug){
  // Return array
  vector<double> returned ((int)array.size(), 0.0);
  debug=0;
  // Parameters:
  double A = x[0]; // Amplitude
  double m = x[1]; // Time shift
  // for (int i = 0; i < (int)array.size(); i++) returned.push_back(0.0);
  vector<double> dummy ((int)array.size(), 0.0);
  // for (int i = 0; i < (int)array.size(); i++) dummy.push_back(0.0);
  // Check if passed parameters are healthy
  if ( (int)array.size() == 0 ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed array is empty. Return empty array.\n"); 
    returned.clear();
    return(returned);
  }
  if ( (int)x.size() != 2 ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed parameter array is of wrong dimension (has to be N=2). Return empty array.\n"); 
    returned.clear();
    return(returned);
  }
  // Check position adjustment parameter is an "integer"
  if ( floor(m) != ceil(m) ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed position parameter x[1] is not an integer. Rounded to %d.\n", (int)round(x[1]));
    m = (double)round(x[1]);
  }
  // Check amplitude adjustment parameter is an greater 0 (not 0 or negative)
  if ( A <= 0 ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed amplitude parameter x[0] is not greter 0. Return x[0] = 1.\n"); 
    A = 0.0;
  }
  // Adjust amplitude 
  for (int i = 0; i < (int)array.size(); i++){
    dummy[i] = array[i] * A;
  }
  // If shift adjustment is zero then dont shift and return array here
  if ( m == 0 ){
    if (debug == 1 ) if (debug == 1) printf("DEBUG(array_adjust): m = 0\n");
    return(dummy);
  }  
  // Shift array by m indizes, the sign is important: - is shift to lower indices, + is shift to higher indices
  // If shift is negative:
  if ( m < 0 ) {
    for (int i = 0; i < (int)array.size()-(int)(abs(m)); i++){
      returned[i] = dummy[i+abs(m)]; // The rightmost m values of the array are 0
    } 
  }
  // If shift is positive:
  if ( m > 0 ) {
    for (int i = 0; i < (int)array.size()-(int)(abs(m)); i++){
      returned[i+abs(m)] = dummy[i]; // The leftmost m values of the array are 0
    } 
  }
  // Test if the adjusted array has the same dimension as the input array
  if ( (int)returned.size() != (int)array.size() ){
    if (debug == 1 ) printf("WARNING(array_adjust): After processing, the adjusted arrray lost initial dimension (bug in code).\n"); 
    returned.clear();
    return(returned);
  }
  // if everything is fine, return valid array
  else{
    return(returned);
  }
}

// Function for initializing the global signal contaiers
void init_signal(vector<signal_struct> &signal, int channels, bool is_raw){
  // create new items in vector
  for(int i=0; i<channels; i++){
    // for each channel, add a channel item
    signal.push_back( signal_struct() );
    // Initialize the proto trace with 0s
    for ( int n = 0; n < TRACELEN; n++) {
      signal[i].proto_trace.push_back(0.0); 
      signal[i].TH1D_proto_trace.push_back( hist_struct_TH1D() );
    }
    // If RAW container, set the RAW flag
    signal[i].is_raw = is_raw;
    // for each channel, add relative time items
    for (int j=0; j<channels; j++){
      signal[i].time.push_back( time_struct() );
      // for each time item, add different energy time items
      for (int k=0; k<N_E_WINDOW; k++){
        signal[i].time[j].h_timing.push_back( hist_struct_TH1D() );
        signal[i].time[j].timing.push_back( 0.0 );
      }
    }
    // In beam mode, also add tagger energy histograms
    if (strcmp(MODE, "BEAM") == 0){
      for (Int_t k=0; k<N_E_WINDOW; k++){
        signal[i].tagged.push_back(tagger_energy());
      }
    }
  }
}

// Function for initializing the global multisampling normalization contaiers
void init_multis_norm(vector<vector<vector<multis_norm_struct> > > &array, int channels){
  // create new items in vector
  // Fill each channel with 2D vectors
  for(int i=0; i<channels; i++){
    // Create MULTIS columns for MULTIS rows
    vector<vector<multis_norm_struct> > temp2;
    for (int j=0; j<4; j++){
      // Create MULTIS rows
      vector<multis_norm_struct> temp1;
      for (int k=0; k<4; k++){
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
  int bins_pulsheight=40000;
  int range_pulsheight=80000*GENERAL_SCALING;
  int bins_integral=10000;
  int range_integral=1000000*GENERAL_SCALING;
  int bins_ratio=4000;
  int range_ratio=4000;
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(1111);
  // If in multi sampling calibration mode, initialize the following
  if (MULTIS_CALIB_MODE == true){
    hfile->cd("ENERGY/RAW/PULSE_HIGHT/");  
    for (Int_t i=0; i<channels; i++){
      hfile->cd("ENERGY/RAW/PULSE_HIGHT/");
      sprintf(name,"ENERGY_RAW%02d",i);
      RAW[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,(int)(range_pulsheight/GENERAL_SCALING));
      RAW[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
    }
    hfile->cd("ENERGY/RAW/Intersampling_calibration");
    for (Int_t i=0; i<CHANNELS; i++){
      for (Int_t j=0; j<4; j++){ 
        for (Int_t k=0; k<4; k++){ 
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
    // The RAW containers always need full channel range
    for (Int_t i=0; i<CHANNELS; i++){
      hfile->cd("ENERGY/PULSE_HIGHT/RAW");
      sprintf(name,"ENERGY_RAW%02d",i);
      RAW[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,(int)(range_pulsheight/GENERAL_SCALING));
      RAW[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"BASELINE_RAW_MEAN%02d",i);
      RAW[i].base.h_mean.hist=new TH1D(name,"",bins_pulsheight+20000,0,(int)(range_pulsheight/GENERAL_SCALING)+20000*GENERAL_SCALING);
      RAW[i].base.h_mean.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].base.h_mean.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"BASELINE_RAW_STD%02d",i);
      RAW[i].base.h_std.hist=new TH1D(name,"",bins_pulsheight,0,(int)(range_pulsheight/GENERAL_SCALING));
      RAW[i].base.h_std.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].base.h_std.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"BASELINE_SAMPLE_DISTRIBUTION%02d",i);
      RAW[i].base.h_samples.hist=new TH1D(name,"",20000,-10000,10000);
      RAW[i].base.h_samples.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW[i].base.h_samples.hist->GetYaxis()->SetTitle("Counts");
    }
    for (Int_t i=0; i<channels; i++){
      //
      // Calibrated raw extraction
      //
      // General energy
      hfile->cd("ENERGY/PULSE_HIGHT/RAW_CALIB");
      sprintf(name,"ENERGY_RAW_CALIB%02d",i);
      RAW_CALIB[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
      RAW_CALIB[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      RAW_CALIB[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/INTEGRAL/RAW_CALIB");
      sprintf(name,"INTEGRAL_RAW_CALIB%02d",i);
      RAW_CALIB[i].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
      RAW_CALIB[i].h_integral.hist->GetXaxis()->SetTitle("Pulse_integral / a.u.");
      RAW_CALIB[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"RATIO_RAW_CALIB%02d",i);
      RAW_CALIB[i].h_ratio.hist=new TH1D(name,"",bins_ratio,0,range_ratio);
      RAW_CALIB[i].h_ratio.hist->GetXaxis()->SetTitle("Pulse_ratio / a.u.");
      RAW_CALIB[i].h_ratio.hist->GetYaxis()->SetTitle("Counts");
      // Tagged energy
      if (strcmp(MODE, "BEAM") == 0){
        hfile->cd("ENERGY/TAGGER/PULSE_HIGHT/RAW_CALIB");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"ENERGY_RAW_CALIB_TAGGER%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          RAW_CALIB[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"ENERGY_RAW_CALIB_TAGGER_M%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          RAW_CALIB[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"ENERGY_RAW_CALIB_TAGGER_MT%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          RAW_CALIB[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        //
        hfile->cd("ENERGY/TAGGER/INTEGRAL/RAW_CALIB");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"INTEGRAL_RAW_CALIB_TAGGER%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          RAW_CALIB[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"INTEGRAL_RAW_CALIB_TAGGER_M%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          RAW_CALIB[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"INTEGRAL_RAW_CALIB_TAGGER_MT%02d_%02d",i, k);
          RAW_CALIB[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          RAW_CALIB[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          RAW_CALIB[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // 2D Hists for the proto trace
      sprintf(name,"WAVE_FORMS/RAW_CALIB/PROTO_TRACE_%02d", i);
      hfile->cd(name);
      sprintf(name,"PROTO_TRACE%02d",i);
      RAW_CALIB[i].TH2D_proto_trace.hist=new TH2D(name,"",
                      300,0,300, // x-dimension
                      1000,-5000,10000); // y-dimension

      // Timing
      hfile->cd("TIMING/RAW_CALIB/HISTOGRAMS");
      for (Int_t j=0; j<(int)RAW_CALIB[i].time.size(); j++){ 
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"RAW_CALIB_TIMINGg_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          RAW_CALIB[i].time[j].h_timing[k].hist = new TH1D(name,"",1000,-100,100);
          RAW_CALIB[i].time[j].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          RAW_CALIB[i].time[j].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
      //
      // Enter the MA folder
      //
      hfile->cd("ENERGY/PULSE_HIGHT/MA");
      sprintf(name,"ENERGY_MA%02d",i);
      gStyle->SetOptFit(1112);
      MA[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
      MA[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      MA[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/INTEGRAL/MA");
      sprintf(name,"INTEGRAL_MA%02d",i);
      MA[i].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
      MA[i].h_integral.hist->GetXaxis()->SetTitle("Pulse_integral / a.u.");
      MA[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"RATIO_MA%02d",i);
      MA[i].h_ratio.hist=new TH1D(name,"",bins_ratio,0,range_ratio);
      MA[i].h_ratio.hist->GetXaxis()->SetTitle("Pulse_ratio / a.u.");
      MA[i].h_ratio.hist->GetYaxis()->SetTitle("Counts");
      // Tagged energy
      if (strcmp(MODE, "BEAM") == 0){
        hfile->cd("ENERGY/TAGGER/PULSE_HIGHT/MA");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MA_TAGGER%02d_%02d",i, k);
          MA[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MA[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_TAGGER_M%02d_%02d",i, k);
          MA[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MA[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_TAGGER_MT%02d_%02d",i, k);
          MA[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MA[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        hfile->cd("ENERGY/TAGGER/INTEGRAL/MA");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MA_TAGGER%02d_%02d",i, k);
          MA[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MA[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_TAGGER_M%02d_%02d",i, k);
          MA[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MA[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_TAGGER_MT%02d_%02d",i, k);
          MA[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MA[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MA[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // Timing
      hfile->cd("TIMING/MA/HISTOGRAMS");
      for (Int_t j=0; j<(int)MA[i].time.size(); j++){ 
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MA_TIMINGg_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          MA[i].time[j].h_timing[k].hist = new TH1D(name,"",1000,-100,100);
          MA[i].time[j].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          MA[i].time[j].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // Enter the MWD folder
      hfile->cd("ENERGY/PULSE_HIGHT/MWD");
      sprintf(name,"ENERGY_MWD%02d",i);
      MWD[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
      MWD[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      MWD[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/INTEGRAL/MWD");
      sprintf(name,"INTEGRAL_MWD%02d",i);
      MWD[i].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
      MWD[i].h_integral.hist->GetXaxis()->SetTitle("Pulse_integral / a.u.");
      MWD[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"RATIO_MWD%02d",i);
      MWD[i].h_ratio.hist=new TH1D(name,"",bins_ratio,0,range_ratio);
      MWD[i].h_ratio.hist->GetXaxis()->SetTitle("Pulse_ratio / a.u.");
      MWD[i].h_ratio.hist->GetYaxis()->SetTitle("Counts");
      // Tagged energy
      if (strcmp(MODE, "BEAM") == 0){
        hfile->cd("ENERGY/TAGGER/PULSE_HIGHT/MWD");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MWD_TAGGER%02d_%02d",i, k);
          MWD[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MWD[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_TAGGER_M%02d_%02d",i, k);
          MWD[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MWD[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_TAGGER_MT%02d_%02d",i, k);
          MWD[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          MWD[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        hfile->cd("ENERGY/TAGGER/INTEGRAL/MWD");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MWD_TAGGER%02d_%02d",i, k);
          MWD[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MWD[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_TAGGER_M%02d_%02d",i, k);
          MWD[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MWD[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_TAGGER_MT%02d_%02d",i, k);
          MWD[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          MWD[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          MWD[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // Timing
      hfile->cd("TIMING/MWD/HISTOGRAMS");
      for (Int_t j=0; j<(int)MWD[i].time.size(); j++){ 
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"MWD_TIMINGg_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          MWD[i].time[j].h_timing[k].hist = new TH1D(name,"",1000,-100,100);
          MWD[i].time[j].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          MWD[i].time[j].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
      //
      // Enter the TMAX folder
      //
      hfile->cd("ENERGY/PULSE_HIGHT/TMAX");
      sprintf(name,"ENERGY_TMAX%02d",i);
      TMAX[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
      TMAX[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      TMAX[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/INTEGRAL/TMAX");
      sprintf(name,"INTEGRAL_TMAX%02d",i);
      TMAX[i].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
      TMAX[i].h_integral.hist->GetXaxis()->SetTitle("Pulse_integral / a.u.");
      TMAX[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"RATIO_TMAX%02d",i);
      TMAX[i].h_ratio.hist=new TH1D(name,"",bins_ratio,0,range_ratio);
      TMAX[i].h_ratio.hist->GetXaxis()->SetTitle("Pulse_ratio / a.u.");
      TMAX[i].h_ratio.hist->GetYaxis()->SetTitle("Counts");
      // Tagged energy
      if (strcmp(MODE, "BEAM") == 0){
        hfile->cd("ENERGY/TAGGER/PULSE_HIGHT/TMAX");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"TMAX_TAGGER%02d_%02d",i, k);
          TMAX[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          TMAX[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_TAGGER_M%02d_%02d",i, k);
          TMAX[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          TMAX[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_TAGGER_MT%02d_%02d",i, k);
          TMAX[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          TMAX[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        hfile->cd("ENERGY/TAGGER/INTEGRAL/TMAX");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"TMAX_TAGGER%02d_%02d",i, k);
          TMAX[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          TMAX[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_TAGGER_M%02d_%02d",i, k);
          TMAX[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          TMAX[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_TAGGER_MT%02d_%02d",i, k);
          TMAX[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          TMAX[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          TMAX[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // Timing
      hfile->cd("TIMING/TMAX/HISTOGRAMS");
      for (Int_t j=0; j<(int)TMAX[i].time.size(); j++){ 
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"TMAX_TIMINGg_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          TMAX[i].time[j].h_timing[k].hist = new TH1D(name,"",1000,-100,100);
          TMAX[i].time[j].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          TMAX[i].time[j].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
      //
      // Enter the NMO folder
      //
      hfile->cd("ENERGY/PULSE_HIGHT/NMO");
      sprintf(name,"ENERGY_NMO%02d",i);
      NMO[i].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
      NMO[i].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
      NMO[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
      hfile->cd("ENERGY/INTEGRAL/NMO");
      sprintf(name,"INTEGRAL_NMO%02d",i);
      NMO[i].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
      NMO[i].h_integral.hist->GetXaxis()->SetTitle("Pulse_integral / a.u.");
      NMO[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      sprintf(name,"RATIO_NMO%02d",i);
      NMO[i].h_ratio.hist=new TH1D(name,"",bins_ratio,0,range_ratio);
      NMO[i].h_ratio.hist->GetXaxis()->SetTitle("Pulse_ratio / a.u.");
      NMO[i].h_ratio.hist->GetYaxis()->SetTitle("Counts");
      // Tagged energy
      if (strcmp(MODE, "BEAM") == 0){
        hfile->cd("ENERGY/TAGGER/PULSE_HIGHT/NMO");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"NMO_TAGGER%02d_%02d",i, k);
          NMO[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          NMO[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_TAGGER_M%02d_%02d",i, k);
          NMO[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          NMO[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_TAGGER_MT%02d_%02d",i, k);
          NMO[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          NMO[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        hfile->cd("ENERGY/TAGGER/INTEGRAL/NMO");
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"NMO_TAGGER%02d_%02d",i, k);
          NMO[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          NMO[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_TAGGER_M%02d_%02d",i, k);
          NMO[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          NMO[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_TAGGER_MT%02d_%02d",i, k);
          NMO[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          NMO[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Pulse_height / a.u.");
          NMO[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
      // Timing
      hfile->cd("TIMING/NMO/HISTOGRAMS");
      for (Int_t j=0; j<(int)NMO[i].time.size(); j++){ 
        for (Int_t k=0; k<N_E_WINDOW; k++){ 
          sprintf(name,"NMO_TIMINGg_ENERGY_WINDOW %i-%i XTAL%i/%i",
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k,
            (int) E_WINDOW_LEFT + E_WINDOW_LENGTH*k +E_WINDOW_LENGTH, 
            i, 
            j);
          NMO[i].time[j].h_timing[k].hist = new TH1D(name,"",1000,-100,100);
          NMO[i].time[j].h_timing[k].hist->GetXaxis()->SetTitle("Time difference / ns");
          NMO[i].time[j].h_timing[k].hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }
  }
  // Additional stuff
  if (strcmp(MODE, "BEAM") == 0){
    // also initialize the TAGGER and hists counter (sorry for doing it here...)
    for (int k = 0; k < TAG_CHANNELS; k++){
      TAGGER.t_hist.push_back(hist_struct_TH1D());
      hfile->cd("ENERGY/TAGGER/TIME");
      sprintf(name,"TAGGER_TIME_HIST%02d",k);
      TAGGER.t_hist[k].hist = new TH1D(name,"",10000,-5000,5000);
      TAGGER.t_hist[k].hist->GetXaxis()->SetTitle("Time difference / arb. units");
      TAGGER.t_hist[k].hist->GetYaxis()->SetTitle("Counts");
      // Also initialize couters
      TAGGER.multiples_per_count[k] = 0;
      TAGGER.multiples_per_channel[k] = 0;
    }
    // Tagger time vs full energy in ECAL
    hfile->cd("ENERGY/TAGGER/TIME");
    TAGGER.h_tagger_vs_energy.hist = new TH2D("Tagger_time_vs_energy","Tagger_time_vs_energy",
                      bins_pulsheight,0,range_pulsheight, // x-dimension
                      1000,500,1500); // y-dimension
    TAGGER.h_tagger_vs_energy.hist->GetXaxis()->SetTitle("Energy in ECAL25 [MeV]");
    TAGGER.h_tagger_vs_energy.hist->GetYaxis()->SetTitle("Tagger time in channels");


    //
    // Initialize the detector sum stuff
    //
    for (int i = 0; i < (int)ECAL25.size(); i++){
      // For the RAW_CALIB filter
      if (i == 0){
        // Initialize the energy histograms
        hfile->cd("ENERGY/SUM/PULSE_HIGHT/RAW_CALIB");
        ECAL25[i].h_energy.hist=new TH1D("ENERGY_SUM_RAW_CALIB","",bins_pulsheight,0,range_pulsheight);
        ECAL25[i].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
        hfile->cd("ENERGY/SUM/INTEGRAL/RAW_CALIB");
        ECAL25[i].h_integral.hist=new TH1D("INTEGRAL_SUM_RAW_CALIB","",bins_integral,0,range_integral);
        ECAL25[i].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      }
      if (i == 1){
        // Initialize the energy histograms
        hfile->cd("ENERGY/SUM/PULSE_HIGHT/MA");
        ECAL25[i].h_energy.hist=new TH1D("ENERGY_SUM_MA","",bins_pulsheight,0,range_pulsheight);
        ECAL25[i].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
        hfile->cd("ENERGY/SUM/INTEGRAL/MA");
        ECAL25[i].h_integral.hist=new TH1D("INTEGRAL_SUM_MA","",bins_integral,0,range_integral);
        ECAL25[i].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      }
      if (i == 2){
        // Initialize the energy histograms
        hfile->cd("ENERGY/SUM/PULSE_HIGHT/MWD");
        ECAL25[i].h_energy.hist=new TH1D("ENERGY_SUM_MWD","",bins_pulsheight,0,range_pulsheight);
        ECAL25[i].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
        hfile->cd("ENERGY/SUM/INTEGRAL/MWD");
        ECAL25[i].h_integral.hist=new TH1D("INTEGRAL_SUM_MWD","",bins_integral,0,range_integral);
        ECAL25[i].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      }
      if (i == 3){
        // Initialize the energy histograms
        hfile->cd("ENERGY/SUM/PULSE_HIGHT/TMAX");
        ECAL25[i].h_energy.hist=new TH1D("ENERGY_SUM_TMAX","",bins_pulsheight,0,range_pulsheight);
        ECAL25[i].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
        hfile->cd("ENERGY/SUM/INTEGRAL/TMAX");
        ECAL25[i].h_integral.hist=new TH1D("INTEGRAL_SUM_TMAX","",bins_integral,0,range_integral);
        ECAL25[i].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      }
      if (i == 4){
        // Initialize the energy histograms
        hfile->cd("ENERGY/SUM/PULSE_HIGHT/NMO");
        ECAL25[i].h_energy.hist=new TH1D("ENERGY_SUM_NMO","",bins_pulsheight,0,range_pulsheight);
        ECAL25[i].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_energy.hist->GetYaxis()->SetTitle("Counts");
        hfile->cd("ENERGY/SUM/INTEGRAL/NMO");
        ECAL25[i].h_integral.hist=new TH1D("INTEGRAL_SUM_NMO","",bins_integral,0,range_integral);
        ECAL25[i].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
        ECAL25[i].h_integral.hist->GetYaxis()->SetTitle("Counts");
      }
    }
    // The tagger histograms
    for (Int_t k=0; k<N_E_WINDOW; k++){ 
      for (int i = 0; i < (int)ECAL25.size(); i++){
        if (i==0){
          hfile->cd("ENERGY/SUM/TAGGER/PULSE_HIGHT/RAW_CALIB");
          sprintf(name,"RAW_CALIB_PH_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"RAW_CALIB_PH_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"RAW_CALIB_PH_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
          hfile->cd("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");       
          sprintf(name,"RAW_CALIB_INT_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"RAW_CALIB_INT_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"RAW_CALIB_INT_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        if (i==1){
          hfile->cd("ENERGY/SUM/TAGGER/PULSE_HIGHT/MA");
          sprintf(name,"MA_PH_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_PH_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_PH_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
          hfile->cd("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");
          sprintf(name,"MA_INT_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_INT_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MA_INT_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        if (i==2){
          hfile->cd("ENERGY/SUM/TAGGER/PULSE_HIGHT/MWD");
          sprintf(name,"MWD_PH_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_PH_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_PH_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
          hfile->cd("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");
          sprintf(name,"MWD_INT_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_INT_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"MWD_INT_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        if (i==3){
          hfile->cd("ENERGY/SUM/TAGGER/PULSE_HIGHT/TMAX");
          sprintf(name,"TMAX_PH_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_PH_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_PH_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
          hfile->cd("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");
          sprintf(name,"TMAX_INT_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_INT_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"TMAX_INT_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
        if (i==4){
          hfile->cd("ENERGY/SUM/TAGGER/PULSE_HIGHT/NMO");
          sprintf(name,"NMO_PH_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_energy.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_PH_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_energy_m.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_PH_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt.hist->GetYaxis()->SetTitle("Counts");
          hfile->cd("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");
          sprintf(name,"NMO_INT_SUM_TAGGER_%02d", k);
          ECAL25[i].tagged[k].h_integral.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_INT_SUM_TAGGER_M_%02d", k);
          ECAL25[i].tagged[k].h_integral_m.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_m.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_m.hist->GetYaxis()->SetTitle("Counts");
          sprintf(name,"NMO_INT_SUM_TAGGER_MT_%02d", k);
          ECAL25[i].tagged[k].h_integral_mt.hist=new TH1D(name,"",bins_integral,0,range_integral);
          ECAL25[i].tagged[k].h_integral_mt.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_integral_mt.hist->GetYaxis()->SetTitle("Counts");
        }
      }
    }
  }
  if (strcmp(MODE, "COSMICS") == 0){
    hfile->cd("ENERGY/PULSE_HIGHT/RAW_CALIB");
    sprintf(name,"RAW_CALIB_PARAMETERS");
    CALIB.h_RAW_energy.hist = new TH1D(name,"",100,0,10);
    CALIB.h_RAW_energy.hist->GetXaxis()->SetTitle("Calibration factor / arb. units");
    CALIB.h_RAW_energy.hist->GetYaxis()->SetTitle("Counts");
    CALIB.h_RAW_energy.hist->SetFillColor(kBlue-7);
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
    for (int i = 0; i<CHANNELS; i++){
      // only fill valid channels
      if (RAW[i].is_valid == false) continue;
      for(int j=0; j< 4; j++){
        for(int k=0; k< 4; k++){
          if (MULTIS_NORM[i][j][k].ratio>0.0){ 
            MULTIS_NORM[i][j][k].h_hist.hist->Fill((Double_t)MULTIS_NORM[i][j][k].ratio, 1);
          }
        }
      }
    }
  }
  else{
    // Always fill RAW  all histograms
    for(int i=0; i<(int)RAW.size(); i++){
      if(RAW[i].energy>0) RAW[i].h_energy.hist->Fill(RAW[i].energy);
      if(RAW[i].base.mean>0) RAW[i].base.h_mean.hist->Fill(RAW[i].base.mean);
      if(RAW[i].base.std>0) RAW[i].base.h_std.hist->Fill(RAW[i].base.std);
    }
    // Fill histograms
    for(int i=0; i<(int)RAW_CALIB.size(); i++){
      // Fill the energy histograms
      if(RAW_CALIB[i].energy>0) RAW_CALIB[i].h_energy.hist->Fill(RAW_CALIB[i].energy);
      if(MA[i].energy>0) MA[i].h_energy.hist->Fill(MA[i].energy);
      if(MWD[i].energy>0) MWD[i].h_energy.hist->Fill(MWD[i].energy);
      if(TMAX[i].energy>0) TMAX[i].h_energy.hist->Fill(TMAX[i].energy);
      if(NMO[i].energy>0) NMO[i].h_energy.hist->Fill(NMO[i].energy);
      // Fill the integral and ratio histograms
      if(RAW_CALIB[i].integral>0) RAW_CALIB[i].h_integral.hist->Fill(RAW_CALIB[i].integral);
      if(MA[i].integral>0) MA[i].h_integral.hist->Fill(MA[i].integral);
      if(MWD[i].integral>0) MWD[i].h_integral.hist->Fill(MWD[i].integral);
      if(TMAX[i].integral>0) TMAX[i].h_integral.hist->Fill(TMAX[i].integral);
      if(NMO[i].integral>0) NMO[i].h_integral.hist->Fill(NMO[i].integral);
      //
      if(RAW_CALIB[i].ratio>0) RAW_CALIB[i].h_ratio.hist->Fill(RAW_CALIB[i].ratio);
      if(MA[i].ratio>0) MA[i].h_ratio.hist->Fill(MA[i].ratio);
      if(MWD[i].ratio>0) MWD[i].h_ratio.hist->Fill(MWD[i].ratio);
      if(TMAX[i].ratio>0) TMAX[i].h_ratio.hist->Fill(TMAX[i].ratio);
      if(NMO[i].ratio>0) NMO[i].h_ratio.hist->Fill(NMO[i].ratio);
      // Fill the timing histograms
      for(int j=0; j<(int)RAW_CALIB.size(); j++){
        // Each energy bin
        for(int k=0; k<N_E_WINDOW; k++){
          if (RAW_CALIB[i].time[j].timing[k] != 0.0) RAW_CALIB[i].time[j].h_timing[k].hist->Fill(RAW_CALIB[i].time[j].timing[k]);
          if (MA[i].time[j].timing[k] != 0.0) MA[i].time[j].h_timing[k].hist->Fill(MA[i].time[j].timing[k]);
          if (MWD[i].time[j].timing[k] != 0.0) MWD[i].time[j].h_timing[k].hist->Fill(MWD[i].time[j].timing[k]);
          if (TMAX[i].time[j].timing[k] != 0.0) TMAX[i].time[j].h_timing[k].hist->Fill(TMAX[i].time[j].timing[k]);
          if (NMO[i].time[j].timing[k] != 0.0) NMO[i].time[j].h_timing[k].hist->Fill(NMO[i].time[j].timing[k]);
        }
      }
      // In beam mode also fill the tagger energy histograms
      if (strcmp(MODE, "BEAM") == 0){
        for(int k=0; k<N_E_WINDOW; k++){    
          // General hists
          if (RAW_CALIB[i].tagged[k].energy != 0.0) RAW_CALIB[i].tagged[k].h_energy.hist->Fill(RAW_CALIB[i].tagged[k].energy);
          if (MA[i].tagged[k].energy != 0.0) MA[i].tagged[k].h_energy.hist->Fill(MA[i].tagged[k].energy);
          if (MWD[i].tagged[k].energy != 0.0) MWD[i].tagged[k].h_energy.hist->Fill(MWD[i].tagged[k].energy);
          if (TMAX[i].tagged[k].energy != 0.0) TMAX[i].tagged[k].h_energy.hist->Fill(TMAX[i].tagged[k].energy);
          if (NMO[i].tagged[k].energy != 0.0) NMO[i].tagged[k].h_energy.hist->Fill(NMO[i].tagged[k].energy);
          //
          if (RAW_CALIB[i].tagged[k].integral != 0.0) RAW_CALIB[i].tagged[k].h_integral.hist->Fill(RAW_CALIB[i].tagged[k].integral);
          if (MA[i].tagged[k].integral != 0.0) MA[i].tagged[k].h_integral.hist->Fill(MA[i].tagged[k].integral);
          if (MWD[i].tagged[k].integral != 0.0) MWD[i].tagged[k].h_integral.hist->Fill(MWD[i].tagged[k].integral);
          if (TMAX[i].tagged[k].integral != 0.0) TMAX[i].tagged[k].h_integral.hist->Fill(TMAX[i].tagged[k].integral);
          if (NMO[i].tagged[k].integral != 0.0) NMO[i].tagged[k].h_integral.hist->Fill(NMO[i].tagged[k].integral);
          // W/o multiples
          if (RAW_CALIB[i].tagged[k].energy_m != 0.0) RAW_CALIB[i].tagged[k].h_energy_m.hist->Fill(RAW_CALIB[i].tagged[k].energy_m);
          if (MA[i].tagged[k].energy_m != 0.0) MA[i].tagged[k].h_energy_m.hist->Fill(MA[i].tagged[k].energy_m);
          if (MWD[i].tagged[k].energy_m != 0.0) MWD[i].tagged[k].h_energy_m.hist->Fill(MWD[i].tagged[k].energy_m);
          if (TMAX[i].tagged[k].energy_m != 0.0) TMAX[i].tagged[k].h_energy_m.hist->Fill(TMAX[i].tagged[k].energy_m);
          if (NMO[i].tagged[k].energy_m != 0.0) NMO[i].tagged[k].h_energy_m.hist->Fill(NMO[i].tagged[k].energy_m);
          //
          if (RAW_CALIB[i].tagged[k].integral_m != 0.0) RAW_CALIB[i].tagged[k].h_integral_m.hist->Fill(RAW_CALIB[i].tagged[k].integral_m);
          if (MA[i].tagged[k].integral_m != 0.0) MA[i].tagged[k].h_integral_m.hist->Fill(MA[i].tagged[k].integral_m);
          if (MWD[i].tagged[k].integral_m != 0.0) MWD[i].tagged[k].h_integral_m.hist->Fill(MWD[i].tagged[k].integral_m);
          if (TMAX[i].tagged[k].integral_m != 0.0) TMAX[i].tagged[k].h_integral_m.hist->Fill(TMAX[i].tagged[k].integral_m);
          if (NMO[i].tagged[k].integral_m != 0.0) NMO[i].tagged[k].h_integral_m.hist->Fill(NMO[i].tagged[k].integral_m);
          // w/o multiples and w/ timing window
          if (RAW_CALIB[i].tagged[k].energy_mt != 0.0) RAW_CALIB[i].tagged[k].h_energy_mt.hist->Fill(RAW_CALIB[i].tagged[k].energy_mt);
          if (MA[i].tagged[k].energy_mt != 0.0) MA[i].tagged[k].h_energy_mt.hist->Fill(MA[i].tagged[k].energy_mt);
          if (MWD[i].tagged[k].energy_mt != 0.0) MWD[i].tagged[k].h_energy_mt.hist->Fill(MWD[i].tagged[k].energy_mt);
          if (TMAX[i].tagged[k].energy_mt != 0.0) TMAX[i].tagged[k].h_energy_mt.hist->Fill(TMAX[i].tagged[k].energy_mt);
          if (NMO[i].tagged[k].energy_mt != 0.0) NMO[i].tagged[k].h_energy_mt.hist->Fill(NMO[i].tagged[k].energy_mt);
          //
          if (RAW_CALIB[i].tagged[k].integral_mt != 0.0) RAW_CALIB[i].tagged[k].h_integral_mt.hist->Fill(RAW_CALIB[i].tagged[k].integral_mt);
          if (MA[i].tagged[k].integral_mt != 0.0) MA[i].tagged[k].h_integral_mt.hist->Fill(MA[i].tagged[k].integral_mt);
          if (MWD[i].tagged[k].integral_mt != 0.0) MWD[i].tagged[k].h_integral_mt.hist->Fill(MWD[i].tagged[k].integral_mt);
          if (TMAX[i].tagged[k].integral_mt != 0.0) TMAX[i].tagged[k].h_integral_mt.hist->Fill(TMAX[i].tagged[k].integral_mt);
          if (NMO[i].tagged[k].integral_mt != 0.0) NMO[i].tagged[k].h_integral_mt.hist->Fill(NMO[i].tagged[k].integral_mt);
        }
      }
    }
    // Save the tagger timing in the histograms
    if (strcmp(MODE, "BEAM") == 0){
      for(int k=0; k<TAG_CHANNELS; k++){
        if (TAGGER.time[k] != 0.0) TAGGER.t_hist[k].hist->Fill(TAGGER.time[k]);
      }
    }
    // Save the ECAL sum in the histogram
    if (strcmp(MODE, "BEAM") == 0){
      for(int i=0; i<(int)ECAL25.size(); i++){
        if (ECAL25[i].energy > 0) ECAL25[i].h_energy.hist->Fill(ECAL25[i].energy);
        if (ECAL25[i].integral > 0) ECAL25[i].h_integral.hist->Fill(ECAL25[i].integral);
      }
      // The tagged energy sum
      for(int k=0; k<N_E_WINDOW; k++){
        for(int i=0; i<(int)ECAL25.size(); i++){
          if (ECAL25[i].tagged[k].energy > 0) ECAL25[i].tagged[k].h_energy.hist->Fill( ECAL25[i].tagged[k].energy );
          if (ECAL25[i].tagged[k].energy_m > 0) ECAL25[i].tagged[k].h_energy_m.hist->Fill( ECAL25[i].tagged[k].energy_m );
          if ( i == 3 && LIN_COMP != 0.0) {
            if (ECAL25[i].tagged[k].energy_mt > 0) ECAL25[i].tagged[k].h_energy_mt.hist->Fill( SIPMlinearize(ECAL25[i].tagged[k].energy_mt, LIN_COMP) );
          }
          else {
            if (ECAL25[i].tagged[k].energy_mt > 0) ECAL25[i].tagged[k].h_energy_mt.hist->Fill( ECAL25[i].tagged[k].energy_mt );
          }
          
          //
          if (ECAL25[i].tagged[k].integral > 0) ECAL25[i].tagged[k].h_integral.hist->Fill( ECAL25[i].tagged[k].integral );
          if (ECAL25[i].tagged[k].integral_m > 0) ECAL25[i].tagged[k].h_integral_m.hist->Fill( ECAL25[i].tagged[k].integral_m );
          if (ECAL25[i].tagged[k].integral_mt > 0) ECAL25[i].tagged[k].h_integral_mt.hist->Fill( ECAL25[i].tagged[k].integral_mt );
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
    signal[i].energy_n = 0;
    signal[i].integral = 0;
    signal[i].ratio = 0;
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
    // Reset the timing 
    for (int j = 0; j < (int)signal[i].time.size(); j++){
      for ( int k = 0; k < N_E_WINDOW; k++){
        signal[i].time[j].timing[k] = 0.0;
      }
    }
    // In beam mode reset the tagged energies
    if (strcmp(MODE, "BEAM") == 0){
      for(int k=0; k<N_E_WINDOW; k++){    
        signal[i].tagged[k].energy = 0.0;
        signal[i].tagged[k].energy_m = 0.0;
        signal[i].tagged[k].energy_mt = 0.0;
        signal[i].tagged[k].integral = 0.0;
        signal[i].tagged[k].integral_m = 0.0;
        signal[i].tagged[k].integral_mt = 0.0;
      }
    }
  }
}

// Superimposes f gaussian distributions
// Double_t fpeaks(Double_t *x, Double_t *par) {
//    Double_t result = par[0] + par[1]*x[0];
//    for (Int_t p=0;p<npeaks;p++) {
//       Double_t norm  = par[3*p+2];
//       Double_t mean  = par[3*p+3];
//       Double_t sigma = par[3*p+4];
//       result += norm*TMath::Gaus(x[0],mean,sigma);
//    }
//    return result;
// }

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
  // Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
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
TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, bool silent){
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
  ffit->SetParNames("Width","MP","Area","GSigma");//, "A", "m");
  // ffit->SetParNames("Width","MP","Area","GSigma");
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  // Do the fit (if suppress output, use silent=true aka "Q")
  if (silent) his->Fit(FunName,"RBQ");
  else his->Fit(FunName,"RB");
  // fit within specified range, use ParLimits, "Q" for quite output
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
  Int_t MAXCALLS = 50000;
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
  // Int_t xmin = hist->GetXaxis()->GetXmin(); // Xmin of histogram
  // Int_t xmax = hist->GetXaxis()->GetXmax(); // Xmin of histogram
  for (Int_t i = 0; i < (upper+abs(lower)); i++){
    if (largest < hist->GetBinContent(i)){ 
      largest_bin = i ; // Norm it to the hist range
      largest = hist->GetBinContent(i);
    }
  }
  // Return index of largest bin
  // printf("%d %d %d %d %d\n", xmin, xmax, ncells, lower, largest_bin);
  return(largest_bin);
}

// Do various fits for a 1D histogram 
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func, Double_t lower, Double_t upper, int verbose){
  vector<Double_t> params;
  // Set the fit area around the mean value of the gaussian
  //hist->GetXaxis()->SetRange(hist->GetMean()-range,hist->GetMean()+range);
  // Check for which fit function is chosen
  int n;
  if (strcmp(func, "gaus")==0){n = 3;}
  else if (strcmp(func, "multigaus")==0){n = 3;}
  else if (strcmp(func, "langaus")==0){n = 6;}
  else{ printf("ERROR (fit_hist): Wrong fit function name passed! Currently available: gaus, multigaus, langaus.\n"); }
  // Check if the histograms were filled
  if ( (int) hist->Integral() < 2 ){
    if (is_in_string(VERBOSE, "h")){
      printf("WARNING (fit_hist: %s): Histogram empty!\n", func);     
    }
    for (int i = 0; i<n+2; i++){
      params.push_back(1.);
      params.push_back(1.);
    } 
    return(params);   
  }
  // Now start the fit routines
  // Gaus fit
  if (strcmp(func, "gaus")==0){
    n = 3;
    // Field 0,1: amp, amp_err
    // Field 2,3: mean, mean_err
    // Field 4,5: std, std_err 
    hist->SetAxisRange(lower, upper);
    // Check if subrange of hist is empty
    if ( (int) hist->Integral() < 2 ){
      if (is_in_string(VERBOSE, "h")){
        printf("WARNING (fit_hist: %s, range %d to %d): Histogram empty!\n", func, (int)lower, (int)upper);     
      }
      for (int i = 0; i<n+2; i++){
        params.push_back(1.);
        params.push_back(1.);
      } 
      return(params);   
    }
    // If not, fit the sub range
    hist->Fit("gaus", "Q", "SAME", lower, upper);
    fit = hist->GetFunction("gaus");
    for (Int_t i = 0; i<n; i++){
      params.push_back(fit->GetParameter(i));
      params.push_back(fit->GetParError(i));
    }
    return(params);  
  }
  // Multigaus: uses peaksensing to find peaks in a histgrams and then fits them with a superposition of gaussians
  if (strcmp(func, "multigaus")==0){
    //Use TSpectrum to find the peak candidates
    TSpectrum *s = new TSpectrum(5);
    Int_t upper_range = hist->GetSize() - 2;
    hist->Rebin(20);
    hist->SetAxisRange(200, upper_range);
    Int_t nfound = s->Search(hist,3,"",0.50);
    // printf("Found %d candidate peaks to fit\n",nfound);
    // nfound = 3;
    //Loop on all found peaks. Eliminate peaks at the background level
    Double_t *xpeaks = s->GetPositionX();
    TF1 *g[nfound]; 
    char t_name[1000]="gaus(0)";
    for (int p = 0; p < nfound; p++){
      char g_name[100];
      sprintf(g_name, "m%d", p);
      if (p > 0) sprintf(t_name, "%s+gaus(%d)", t_name, p*3);
      g[p] = new TF1(g_name,"gaus",xpeaks[p]-200,xpeaks[p]+200);
      // printf("%s %s\n", g_name, t_name);
    }
    // The total is the sum of the three, each has 3 parameters
    fit = new TF1("mstotal",t_name,xpeaks[0]-500,xpeaks[nfound-1]+500);
    Double_t par[nfound*3];
    for (int p = 0; p < nfound; p++){
      char g_name[100];
      sprintf(g_name, "m%d", p);
      if (p == 0) hist->Fit(g_name,"0R");
      if (p > 0) hist->Fit(g_name,"0R+");
      g[p]->GetParameters(&par[p*3]);
    }
    fit->SetParameters(par);
    hist->Fit(fit,"R+");
    fit->GetParameters(par);
    for (int i = 0; i < nfound*3; i++){
      params.push_back(fit->GetParameter(i));
      params.push_back(fit->GetParError(i));
    }
    return(params);
  }
  // Langaus (landau convoluted with gaus) fit
  if (strcmp(func, "langaus")==0){
    n = 4;
    // Rebin a temp histogram before fitting
    Int_t nrebin = 100;
    int errors = 0;
    // Setting fit range and start values
    hist->Rebin(nrebin);
    Double_t largest_bin = largest_1Dbin(hist, 2000*GENERAL_SCALING/nrebin)*nrebin*GENERAL_SCALING; // heaviest bin between lower and upper bound
    printf("%3.1f\n", largest_bin);
    Double_t fr[2]; // fit boundaries
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4]; 
    fr[0]=0.75*largest_bin; // Lower fit boundary
    fr[1]=2.0*largest_bin; // Upper fit boundary
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //par[4]=A from A/(x^(m))
    //par[5]=m from A/(x^(m))
    //ADDED LATER: par[6]= Maximum of convoluted function
    //ADDED LATER: par[7]= FWHM of convoluted function
    pllo[0]=1.      ; pllo[1]=0.               ; pllo[2]=1.0          ; pllo[3]=1.     ;// pllo[4]=-10.0 ; pllo[5]=0.0001;  // Lower parameter limits
    plhi[0]=20000.   ; plhi[1]=largest_bin+20000.; plhi[2]=10000000000.0; plhi[3]=10000.0;// plhi[4]=100000.0; plhi[5]=5.0; // Upper parameter limits
    sv[0]  =100.    ; sv[1]  =largest_bin      ; sv[2]  =5000000.0    ; sv[3]  =1000.0 ;// sv[4]  =100.0   ; sv[5]=0.05;// Start values
    Double_t chisqr; // Chi squared
    Int_t ndf; // # degrees of freedom
    bool silent = true;
    if (is_in_string(VERBOSE, "f")) silent = false; 
    fit = langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf, silent); // true/false for quite mode
    Double_t SNRPeak, SNRFWHM; 
    errors = langaupro(fp,SNRPeak,SNRFWHM); // Search for peak and FWHM
    if (errors != 0) printf("WARNING(fit_hist): ERROR(langaupro): %d\n", errors);
    // Save the parameters in the return container
    for (int i = 0; i<n; i++){ 
      params.push_back(fp[i]);
      params.push_back(fpe[i]);
    }
    // Add the parameters from the langaupro search
    params.push_back(SNRPeak); params.push_back(SNRFWHM);
    return(params);
  }
  // Novosibirsk function
  // if (strcmp(func, "novosibirsk")==0){
  //   // Declare observable x
  //   int lower = hist->GetXaxis()->GetXmin();
  //   int upper = hist->GetXaxis()->GetXmax();
  //   RooRealVar x("x","x",lower,upper);
  //   // --- Import data ---
  //   RooDataHist data("data","data",x,Import(*hist)) ;
  //   // --- Parameters ---
  //   RooRealVar mean("sigmean","mean",5.28,5.20,5.30);
  //   RooRealVar width("sigwidth","width",0.0027,0.001,1.);
  //   RooRealVar tail("sigwidth","tail",0.0027,0.001,1.);

  //   return(params);
  // }


  // 
  else {return(params);}
}

vector<Double_t> fit_graph_err(TGraphErrors *graph, char const *func, Double_t lower, Double_t upper, int verbose){
  //
  //
  vector<Double_t> params;
  int n = 0;
  TF1 *fit;
  //
  if (strcmp(func, "SIPMpixel")==0){
    n = 1;
    Double_t par[1];
    par[0] = 1000.;
    fit = new TF1( "SIPMpixel_fit", SIPMpixel,  lower, upper, 1);
    // printf("%3.3f\n", par[0]);
    fit->SetParameters(par);
    fit->SetParNames("N_tot");
    graph->Fit("SIPMpixel_fit", "R", "SAME", lower, upper);

    for (int i = 0; i < n; i ++){
      params.push_back(fit->GetParameter(i));
      params.push_back(fit->GetParError(i));
    }
  }
  //
  //
  if (strcmp(func, "resolution")==0){
    // par[0]: A is stochastic shower development parameter
    // par[1]: B is parameter for electric noice  
    // par[2]: C is free parameter
    // par[3]: D is third degree parameter
    // par[4]: x_0 is pole shift parameter
    n = 5;
    //
    double par_start[n+1], par_low[n+1], par_up[n+1];
    par_low[0] = 0.001; par_start[0] = 430; par_up[0] = 1000;
    par_low[1] = 0.001; par_start[1] = 2200; par_up[1] = 50000;
    par_low[2] = 1; par_start[2] = 15; par_up[2] = 50;
    par_low[3] = 0.001; par_start[3] = 1000; par_up[3] = 50000;
    par_low[4] = 0.01; par_start[4] = 1; par_up[4] = 25;
    //
    fit = new TF1( "Resolution_fit", resolution,  lower, upper, n);
    // printf("%3.3f\n", par[0]);
    fit->SetParameters(par_start);
    fit->SetParNames("A", "B", "C", "D", "x_{0}");
    //
    for (int i=0; i<n; i++) {
      fit->SetParLimits(i, par_low[i], par_up[i]);
    }

    //
    graph->Fit("Resolution_fit", "R", "SAME", lower, upper);

    for (int i = 0; i < n; i ++){
      params.push_back(fit->GetParameter(i));
      params.push_back(fit->GetParError(i));
    }



  }
  return params;
}

// Fast linear regression fit 
// Source: http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
bool linreg(vector<double> &x, vector<double> &y, double *m, double *b){
  // double r; // Correlation coefficient
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
  // r = (sumxy - sumx * sumy / n) /    /* compute correlation coeff */
  //       sqrt((sumx2 - pow(sumx, 2.)/n) *
  //       (sumy2 - pow(sumy, 2.)/n));
  return true; 
}

// Finite Impulse Response Filter
vector<double> FIR_filter(vector<double> &trace, double calib){
  // Return array
  vector<double> filtered;
  // Number of tabs
  nTabs = (int)FIR_COEF.size();
  // total sum
  double sum = 0.;
  // For each sample of the trace
  for (int n = 0; n < (int)trace.size(); n++){
    // reset the sum
    sum = 0.;
    // For each tab
    for (int k = 0; k < nTabs; k++){
      // Check the border condition
      if ( (n-k) < 0 ) continue;
      sum += trace[n - k] * FIR_COEF[k];
    }
    // Push back the convolution sum multiplied with the energy calibration
    sum *= calib;
    filtered.push_back(sum);
  }
  return filtered;
}

// Moving average Filter
vector<double> MA_filter(vector<double> &trace, double calib){
  // Return array
  vector<double> filtered;
  double value = 0;
  // For each sample of the trac do:
	for(int n=0; n < (int)trace.size(); n++){
		value = 0;
		// Take care of boundary conditions
		// Left boundary
	  if( n - ((L-1)/2) < 1 ){
	    value = array_mean(trace, 0, n + (L-1)/2);
	  } 
	  // Center values
	  if( n - ((L-1)/2) >= 1 && n + ((L-1)/2) <= (int)trace.size() ){ 
	    value = array_mean(trace, n - ((L-1)/2), n + ((L-1)/2));
	  }
	  // Right boundary
	  if( n + ((L-1)/2) > (int)trace.size() ){ 
	    value = array_mean(trace, n - ((L-1)/2), (int)trace.size());
	  }
	  // Calibration
	  value *= calib;
	  // Save value in array
	  filtered.push_back(value);
	}
	return filtered;
}

// Moving average Filter
vector<double> MWD_filter(vector<double> &trace, double calib){
  // Return array
  vector<double> filtered;
  double value = 0;
  // For each sample of the trac do:
	for(int n=0; n < (int)trace.size(); n++){
		value = 0;
    // Apply Moving Window Deconvolution Filter, take care of boundaries
    if (n == 0 || n == 1){
      value = trace[n];
    }
    if (n-M <= 0 && n>1){
      value = trace[n] - trace[0] + (1/TAU)*array_mean(trace, 0, n-1);        
    }
    if (n-M > 0){
      value = trace[n] - trace[n-M] + (1/TAU)*array_mean(trace, n-M, n-1);
    }
	  // Calibration
	  value *= calib;
	  // Save value in array
	  filtered.push_back(value);
	}
	return filtered;
}

// Digital Constant Fraction Discriminator 
vector<double> CFD(vector<double> &trace, double multiplier){
  // Return array
  vector<double> converted;
  // Adjust DELAY to the tracelength
  int delay = (int)((multiplier)*DELAY);
  // For each sample of the trac do:
	for(int n=0; n < (int)trace.size(); n++){
    // Calculate the Constant Fraction Trace
		if (n-delay>0) converted.push_back(trace[n-delay] - CFD_fraction * trace[n]);
    else {converted.push_back(0);}	  // Calibration
	}
	return converted;
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
  printf("  - v (str)      : Turns on verbose mode. Because more ouput is always better!\n"); 
  printf("      + p: runtime progress reports\n");
  printf("      + c: config file parameters\n");
  printf("      + f: fit routine output\n");
  printf("      + h: histogram operation output\n");
  printf("      + e: energy calibration/extraction ouput (only cosmics mode)\n");
  printf("      + t: timing statistics\n");
  printf("      + g: tagger statistics (only beam mode)\n");
  printf("  - M (odd int)  : Range of moving window deconvolution window.\n"); 
  printf("  - I (1,2, or 4): Number of inter samples.\n"); 
  printf("  - f (double,0<f<1): Fraction for CFD algorithm.\n"); 
  printf("  - e            : Runs the program in inter sampling calibration mode. Needs option I activated.\n"); 
  printf("  - h            : Shows this help page.\n"); 
  printf("\n\n");
}

// Builds folder structure
void build_structure(){
  // Energy
  //
  hfile->mkdir("ENERGY");
  hfile->mkdir("ENERGY/PULSE_HIGHT");
  hfile->mkdir("ENERGY/PULSE_HIGHT/RAW");
  hfile->mkdir("ENERGY/PULSE_HIGHT/RAW/Intersampling_calibration");
  hfile->mkdir("ENERGY/PULSE_HIGHT/RAW_CALIB");
  hfile->mkdir("ENERGY/PULSE_HIGHT/MA");
  hfile->mkdir("ENERGY/PULSE_HIGHT/MWD");
  hfile->mkdir("ENERGY/PULSE_HIGHT/TMAX");
  hfile->mkdir("ENERGY/PULSE_HIGHT/NMO");
  hfile->mkdir("ENERGY/INTEGRAL");
  hfile->mkdir("ENERGY/INTEGRAL/RAW_CALIB");
  hfile->mkdir("ENERGY/INTEGRAL/MA");
  hfile->mkdir("ENERGY/INTEGRAL/MWD");
  hfile->mkdir("ENERGY/INTEGRAL/TMAX");
  hfile->mkdir("ENERGY/INTEGRAL/NMO");
  // In case a tagger is used in beam mode
  if (strcmp(MODE, "BEAM") == 0){
    hfile->mkdir("ENERGY/TAGGER");
    //
    hfile->mkdir("ENERGY/TAGGER/TIME");
    //
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT");
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT/RAW_CALIB");
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT/MA");
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT/MWD");
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT/TMAX");
    hfile->mkdir("ENERGY/TAGGER/PULSE_HIGHT/NMO");
    //
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL");
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL/RAW_CALIB");
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL/MA");
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL/MWD");
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL/TMAX");
    hfile->mkdir("ENERGY/TAGGER/INTEGRAL/NMO");
    // For the energy sum
    hfile->mkdir("ENERGY/SUM");
    //
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/RAW");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/RAW/Intersampling_calibration");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/RAW_CALIB");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/MA");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/MWD");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/TMAX");
    hfile->mkdir("ENERGY/SUM/PULSE_HIGHT/NMO");
    hfile->mkdir("ENERGY/SUM/INTEGRAL");
    hfile->mkdir("ENERGY/SUM/INTEGRAL/RAW_CALIB");
    hfile->mkdir("ENERGY/SUM/INTEGRAL/MA");
    hfile->mkdir("ENERGY/SUM/INTEGRAL/MWD");
    hfile->mkdir("ENERGY/SUM/INTEGRAL/TMAX");
    hfile->mkdir("ENERGY/SUM/INTEGRAL/NMO");
    //
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT");
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT/RAW_CALIB");
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT/MA");
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT/MWD");
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT/TMAX");
    hfile->mkdir("ENERGY/SUM/TAGGER/PULSE_HIGHT/NMO");
    //
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL");
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL/RAW_CALIB");
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL/MA");
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL/MWD");
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL/TMAX");
    hfile->mkdir("ENERGY/SUM/TAGGER/INTEGRAL/NMO");


  }
  // Timing hirachy
  // Timing hirachy
  //
  hfile->mkdir("TIMING");
  hfile->mkdir("TIMING/RAW_CALIB");
  hfile->mkdir("TIMING/RAW_CALIB/HISTOGRAMS");
  hfile->mkdir("TIMING/MA");
  hfile->mkdir("TIMING/MA/HISTOGRAMS");
  hfile->mkdir("TIMING/MWD");
  hfile->mkdir("TIMING/MWD/HISTOGRAMS");
  hfile->mkdir("TIMING/TMAX");
  hfile->mkdir("TIMING/TMAX/HISTOGRAMS");
  hfile->mkdir("TIMING/NMO");
  hfile->mkdir("TIMING/NMO/HISTOGRAMS");
  // Create subfolders for each channel
  char name[100];
  for (int i = 0; i < CHANNELS_EFF; i++){
    sprintf(name,"TIMING/RAW_CALIB/CHANNEL_%i", i);
    hfile->mkdir(name);
    sprintf(name,"TIMING/MA/CHANNEL_%i", i);
    hfile->mkdir(name);
    sprintf(name,"TIMING/MWD/CHANNEL_%i", i);
    hfile->mkdir(name);
    sprintf(name,"TIMING/TMAX/CHANNEL_%i", i);
    hfile->mkdir(name);
    sprintf(name,"TIMING/NMO/CHANNEL_%i", i);
    hfile->mkdir(name);
  }

  // Wave forms
  //
  hfile->mkdir("WAVE_FORMS");
  hfile->mkdir("WAVE_FORMS/RAW");
  hfile->mkdir("WAVE_FORMS/RAW_CALIB");
  for (int i = 0; i < CHANNELS_EFF; i++){
    sprintf(name,"WAVE_FORMS/RAW_CALIB/PROTO_TRACE_%02d", i);
    hfile->mkdir(name);
  }
  hfile->mkdir("WAVE_FORMS/MA");
  hfile->mkdir("WAVE_FORMS/MWD");
  hfile->mkdir("WAVE_FORMS/TMAX");
  hfile->mkdir("WAVE_FORMS/NMO");
  hfile->mkdir("WAVE_FORMS/CFD");
  hfile->mkdir("WAVE_FORMS/CFD/RAW");
  hfile->mkdir("WAVE_FORMS/CFD/RAW_CALIB");
  hfile->mkdir("WAVE_FORMS/CFD/MA");
  hfile->mkdir("WAVE_FORMS/CFD/MA/INTERPOL");
  hfile->mkdir("WAVE_FORMS/CFD/MWD");
  hfile->mkdir("WAVE_FORMS/CFD/MWD/INTERPOL");
  hfile->mkdir("WAVE_FORMS/CFD/TMAX");
  hfile->mkdir("WAVE_FORMS/CFD/TMAX/INTERPOL");
  hfile->mkdir("WAVE_FORMS/CFD/NMO");
  hfile->mkdir("WAVE_FORMS/CFD/NMO/INTERPOL");

  // Junk
  hfile->mkdir("JUNK");
}

bool read_file(const char *file, const char *mode){
  // set status for return function
  ifstream file_in(file);
  // Check if file exists
  if (!file_in.is_open()){
    printf("ERROR (read_file): File %s does not exist!\n", file);
    return(false);
  }
  else {
    if (strcmp( mode, "FIR" )==0) printf("read_file: FIR coefficients read in under: %s\n", file);
    if (strcmp( mode, "PROTO" )==0) printf("read_file: PROTO trace read in under: %s\n", file);
  }
  // If yes, fill the parameters in the vector
  string line;
  int i = 0;
  int n = 0;
  //
  while(getline(file_in, line) )
  {
    i = 0;
    if (strcmp( mode, "FIR" )==0) FIR_COEF.push_back( stof(line) );
    // Proto traces from file (multiple columns)
    if (strcmp( mode, "PROTO" )==0) {
      istringstream is_line(line);
      string value;
      // cout << line << endl;
      while( getline(is_line, value, ';') )
      {
        // cout << value << endl;
        // printf("%d %d\n", i, n);
        // For the first line, initialize the proto container
        if ( n == 0) PROTO.push_back( signal_struct() ); 
        // Push back the values of each cell
        PROTO[i].proto_trace_fit.push_back( stof(value) );
        i++;
      }
      n++;
    }
  }
  return(true);
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
        istringstream value_line(value);
        //Now test for all possible parameter cases 
        if(strcmp(key.c_str(), "CHANNELS") == 0) CHANNELS = stoi(value);
        if(strcmp(key.c_str(), "BOARDS") == 0) BOARDS = stoi(value);
        if(strcmp(key.c_str(), "TAG_CHANNELS") == 0) TAG_CHANNELS = stoi(value);
        if(strcmp(key.c_str(), "TAGGER_WINDOW") == 0) TAGGER_WINDOW = stoi(value);
        if(strcmp(key.c_str(), "CENTRAL") == 0) CENTRAL = stoi(value);
        if(strcmp(key.c_str(), "TRACELEN") == 0) TRACELEN = stoi(value);
        if(strcmp(key.c_str(), "SAMPLE_t") == 0) SAMPLE_t = stoi(value);
        if(strcmp(key.c_str(), "VERBOSE") == 0){ strcpy(VERBOSE, value.c_str());} 
        if(strcmp(key.c_str(), "MODE") == 0){ strcpy(MODE, value.c_str());} 
        if(strcmp(key.c_str(), "N_E_WINDOW") == 0) N_E_WINDOW = stoi(value);
        if(strcmp(key.c_str(), "E_WINDOW_LEFT") == 0) E_WINDOW_LEFT = stoi(value);
        if(strcmp(key.c_str(), "E_WINDOW_RIGHT") == 0) E_WINDOW_RIGHT = stoi(value);
        if(strcmp(key.c_str(), "BASELINE_CUT") == 0) BASELINE_CUT = stoi(value);
        if(strcmp(key.c_str(), "ENERGY_WINDOW_MAX") == 0) ENERGY_WINDOW_MAX = stoi(value);
        if(strcmp(key.c_str(), "ENERGY_NORM") == 0) ENERGY_NORM = stoi(value);
        if(strcmp(key.c_str(), "N_INTPOL_SAMPLES") == 0) N_INTPOL_SAMPLES = stoi(value);
        if(strcmp(key.c_str(), "NB") == 0) NB = stoi(value);
        if(strcmp(key.c_str(), "COINC_LEVEL") == 0) COINC_LEVEL = stoi(value);
        if(strcmp(key.c_str(), "EXTRACT_PROTO") == 0) EXTRACT_PROTO = stoi(value);
        if(strcmp(key.c_str(), "LIN_COMP") == 0) LIN_COMP = stod(value);
        // if(strcmp(key.c_str(), "MULTIS") == 0) MULTIS = stoi(value);
        if(strcmp(key.c_str(), "THRESHOLD_MULTIPLICY") == 0) THRESHOLD_MULTIPLICY = stod(value);
        if(strcmp(key.c_str(), "SATURATION_FILTER") == 0) SATURATION_FILTER = stoi(value);
        if(strcmp(key.c_str(), "GLITCH_FILTER") == 0) GLITCH_FILTER = stoi(value);
        if(strcmp(key.c_str(), "GLITCH_FILTER_RANGE") == 0) GLITCH_FILTER_RANGE = stoi(value);
        // Feature extraction algorithm parameters
        if(strcmp(key.c_str(), "L") == 0) L = stoi(value);
        if(strcmp(key.c_str(), "DELAY") == 0) DELAY = stoi(value);
        if(strcmp(key.c_str(), "M") == 0) M = stoi(value); 
        if(strcmp(key.c_str(), "TAU") == 0) TAU = stod(value);
        if(strcmp(key.c_str(), "CFD_fraction") == 0) CFD_fraction = stod(value);
        if(strcmp(key.c_str(), "GENERAL_SCALING") == 0) GENERAL_SCALING = stod(value);
        if(strcmp(key.c_str(), "LOWER_RATIO") == 0) LOWER_RATIO = stod(value);
        if(strcmp(key.c_str(), "UPPER_RATIO") == 0) UPPER_RATIO = stod(value);
        // Read in FIR filter coefficients if chosen
        if(strcmp(key.c_str(), "FIR_FILE") == 0) {
          if (!read_file(value.c_str(), "FIR")) {
            printf("ERROR (read_config): FIR filter coef. file error.\n");
            return(false);
          }
        }
        // Read in PROTO trace if chosen
        if(strcmp(key.c_str(), "PROTO_FILE") == 0 && EXTRACT_PROTO == 0) {
          if (!read_file(value.c_str(), "PROTO")) {
            printf("ERROR (read_config): PROTO trace file error.\n");
            return(false);
          }
        }
        // Calibration parameters
        // Check for the MULTIS_CALIB's
        string subkey = key.substr(0, (int)key.size()-2);
        string subcounter;
        if ((int)key.size()-2 > 2) subcounter = key.substr((int)key.size()-2, (int)key.size()-1); // If variable read in is only one letter, key - 2 would throw error
        if(strcmp(subkey.c_str(), "MULTIS_CALIB") == 0){
          // get the calib channel
          CALIB.multis.push_back( stod(value) );
        } 
        if(strcmp(subkey.c_str(), "ENERGY_CALIB") == 0){
          // get the calib channel
          int value_counter = 0;
          while( getline(value_line, value, ',') ){
            // First value is RAW, then MA, MWD, TMAX
            if (value_counter == 0) CALIB.RAW_energy.push_back( stod(value) );
            if (value_counter == 1) CALIB.MA_energy.push_back( stod(value) );
            if (value_counter == 2) CALIB.MWD_energy.push_back( stod(value) );
            if (value_counter == 3) CALIB.TMAX_energy.push_back( stod(value) );
            if (value_counter == 4) CALIB.NMO_energy.push_back( stod(value) ); 
            value_counter++;
          }
        }
        // Check for Tagger time cutting parameters
        if(strcmp(subkey.c_str(), "TAGGER_CUT") == 0){
          for (int k = 0; k < TAG_CHANNELS; k++){
            if ( stoi(subcounter) == k ){
              TAGGER.cut[k] = stoi(value);
            }
          }
        } 
        // If individual channel TH is given, read in the thresh hold multiplicity
        if(strcmp(subkey.c_str(), "CHANNEL_TH") == 0 && THRESHOLD_MULTIPLICY == -1 ){
          TH_MULTIPLICITY.push_back(stod(value));
        } 
        // Read the channel mapping (mapping from hardware to software channel)
        subkey = key.substr(0, (int)key.size()-2);
        if(strcmp(subkey.c_str(), "MAPPING") == 0){
          // New row of crystals, create temporary vector for pushing
          vector<mapping_struct> temp;
          // Now read whole row of crystals
          while( getline(value_line, value, ',') ){
            // Create new element in temp
            temp.push_back(mapping_struct());
            // The hardware-software mapping depends on the multichannel mode
            // If channel is sampled once, then the "X-Y" is X==Y, if not, X!=Y (see config file)
            string delim = "/";
            auto start = 0U;
            auto end = value.find(delim);
            int counter = 0;
            while (end != string::npos)
            {
              // String 0 is the polarity
              if (counter == 0) {
                temp.back().polarity = stoi(value.substr(start, end - start));
                // Count number of effective hardware channels 
                CHANNELS_EFF++;
              }
              // String 1 is the board number 
              if (counter == 1) {
                temp.back().board_nb = stoi(value.substr(start, end - start));
                // Count number of valid channels
                if (stoi(value.substr(start, end - start)) != 99) NB_ACT_CHANNELS++;
              }
              if (counter == 2){
                temp.back().h_channel_start = stoi(value.substr(start, end - start));
              }
              counter++;            
              start = end + delim.length();
              end = value.find(delim, start);
            }
            if (counter == 3){
              temp.back().h_channel_stop = stoi(value.substr(start, end - start));
            }
            // Test if start channel and end channel are the same
            if (temp.back().h_channel_start != temp.back().h_channel_stop) {
              temp.back().multis = temp.back().h_channel_stop - temp.back().h_channel_start + 1;
            }
          }
          // Now push the row of crystals into the mapping field
          MAPPING.push_back(temp);
        }       
      }
    }
  }
  // Calculate effective channel number, effective trace length, and effective sample frequency
  SAMPLE_t_eff = (double) SAMPLE_t;
  if (is_in_string(VERBOSE, "c")){
    printf("\n");
    printf("+++ CONFIG FILE PARAMETERS SET ++\n");
    printf("Analysis mode: %s\n", MODE);
    printf("Number of total channels: %i\n", CHANNELS);
    printf("Number of ADC boards: %i\n", BOARDS);
    printf("Lentgh of moving average interval: %i\n", L);
    printf("Lentgh of moving window interval: %i\n", M);
    printf("Fraction of CFD algorithm: %3.2f\n", CFD_fraction);
    printf("Length of trace: %i\n", TRACELEN);
    printf("Effective number of channels: %i\n", CHANNELS_EFF);
    printf("+++++++++++++++++++++++++++\n\n");
  }
  // Since MA, MWD, and TMAX are calculated from a calibrated RAW trace, the filter calibration vlues
  //    have to be adjusted to the initial calibration to avoid false multiplication of two factors
  for (int i = 0; i<(int)CALIB.RAW_energy.size(); i++){
    CALIB.MA_energy[i] /= CALIB.RAW_energy[i];
    CALIB.MWD_energy[i] /= CALIB.RAW_energy[i];
    CALIB.TMAX_energy[i] /= CALIB.RAW_energy[i];
    CALIB.NMO_energy[i] /= CALIB.RAW_energy[i];
  }
  // Some tests for the read configuration parameters
  // Test if the dimension of the calibration parameters fit the number of channels, etc.
  if ( (int)CALIB.RAW_energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every RAW filter channel has its energy calibration factor!\n");
    return(false);
  }
  if ( (int)CALIB.MA_energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every MA filter channel has its energy calibration factor!\n");
    return(false);
  }
  if ( (int)CALIB.MWD_energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every MWD filter channel has its energy calibration factor!\n");
    return(false);
  }
  if ( (int)CALIB.TMAX_energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every TMAX filter channel has its energy calibration factor!\n");
    return(false);
  }
  if ( (int)CALIB.NMO_energy.size() < CHANNELS ){
    printf("ERROR (read_config): Not every NMO filter channel has its energy calibration factor!\n");
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
  if (strcmp(MODE, "COSMICS") != 0 &&
      strcmp(MODE, "PULSER") != 0 && 
      strcmp(MODE, "MULTIS") != 0 &&
      strcmp(MODE, "BEAM") != 0) {
    printf("ERROR (read_config): Wrong MODE setting! Available: COSMICS, MULTIS, PULSER, BEAM. And don't forget the white space at the end of the line!\n");
    return(false);
  }


  return(true);
}

// Prints the hardware to software mapping of the ADC
void print_detector_config(){
  printf("++ HARDWARE TO SOFTWARE MAPPING ++\n");
  printf("+ HARDWARE CHANNEL +\n");
  for (int i = 0; i<(int)MAPPING.size(); i++){
    for (int j = 0; j<(int)MAPPING[i].size(); j++){
      printf("%i/%i/%i/%i ", MAPPING[i][j].polarity, MAPPING[i][j].board_nb, MAPPING[i][j].h_channel_start, MAPPING[i][j].h_channel_stop);
    }
    printf("\n");
  }
  printf("\n");
  printf("+ SOFTWARE CHANNEL +\n");
  for (int i = 0; i<(int)MAPPING.size(); i++){
    for (int j = 0; j<(int)MAPPING[i].size(); j++){
      printf("%2d ", (int)MAPPING.size()*j + i);
    }
    printf("\n");
  }
  printf("\n");
}

// Searches string if character is present
bool is_in_string(char const *character, char const *letter){
  // convert char to str
  string str(character);
  // look for character
  size_t found = str.find_first_of(letter);
  // printf("%d\n", (int)found);
  // cout << character << " " << letter << " " << found << endl;
  if ((int)found == -1) return false;
  else return true;
}

// Checks if datapoint (assumed maximum) is glitch
bool is_glitch(vector<double> &trace, double TH, int n){
  int is_glitch_array[GLITCH_FILTER_RANGE];
  int is_glitch = 1;
  // Check for the next few samples if all are above threshold
  for(int k=n; k<n+GLITCH_FILTER_RANGE; k++){
    // If k-th sample is greater than TH, then save 0 in array
    // Test with RAW signal because RAW_CALIB is not yet calculated
    if ( abs(trace[k]) > TH){
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
  // If not a glitch, return false
  if (is_glitch==1){
    return false;
  }
  // If it's a glitch return true
  else{return(true);}
}

// Checks if datapoint (assumed maximum) is saturation
	// Look at k% how wide the pulse is and decide
bool is_saturation(signal_struct &signal, int n){
	double max = signal.trace[n];
	double fraction = 0.95; // fraction of the pulse hight
	int width = 7; // width at this fraction
	int rising_width = 7; // width of signal rising edge
	int width_left = 0;  
	int width_right = 0;
	int rising_left = 0;
	int rising_right = n;
	int rising_start = 80; // powerful, but use carefully
	double sat_threshold = 10000;
	// Search for the rising edge crossing the kth % of the maximum
	for (int i = 0; i < ENERGY_WINDOW_MAX; i++){
		if ( signal.trace[i] > fraction * max) {
			width_left = i-1;
			break;
		}
	}
	// Search for the falling edge crossing the kth % of the maximum
	for (int i = width_left + 1; i < ENERGY_WINDOW_MAX; i++){
		if ( signal.trace[i] < fraction * max) {
			width_right = i;
			break;
		}
	}
	// Check how fast the rising edge is by going left starting at max
	for (int i = 0; i < n; i++){
		if ( signal.trace[n-i] < signal.base.TH ){
			rising_left = n-i;
			break;
		}
	}
	// now check if the width is smaller than the typical PWO shape width 
	// printf("%d\n", width_right-width_left);
	if ( (width_right - width_left < width && // Width at k% smaller than normal signal
         signal.trace[n] > sat_threshold && // Max larger than certain fixed TH
         rising_right - rising_left < rising_width) || // Rising edge smaller than normal signal
			 (rising_left < rising_start) // Signal starts before a certain sample number
			){
		// printf("%3.2f %d %d\n", trace[n], width_left, width_right);
		return(true);
	} 
	else{
		return(false);		
	}
}

// Check if baseline is weird (signal starts befor baseline cut, baseline is tilted,...)
bool baseline_weird(signal_struct &signal){
  //
  int rising_left = 0;
  // double mean;
  // Check where the rising edge starts
  for (int i = 0; i < signal.energy_n; i++){
    if ( signal.trace[signal.energy_n - i] < signal.base.TH ){
      rising_left = signal.energy_n - i;
      break;
    }
  }
  // If start is before baseline cut, signal is weird
  if ( rising_left < BASELINE_CUT ) return(true);
  // //
  // // Check if baseline is by calculating the actual mean of all samples before baseline starts
  // mean = array_mean(signal.trace, 0, rising_left);
  // // if the mean value is larger than a fraction of the initial TH calculation then baseline is weird
  // if ( abs(mean) > abs(0.20 * signal.base.TH) ) return(true);
  // //
  // //if nothing mathes the criteria, baseline is ok
  return(false);
}


// Test, if the sample is a valid sample and above TH.
// 	Includes glitch and saturation filter
// Return parameter code: 
		// 0: Is valid max 
		// 1: Glitch detected 
		// 2: Saturation detected 
		// 3: Trace is not not above TH
    // 4: Baseline is wierd
int is_valid_max(signal_struct &signal, int n){
	// Test if the current sample is larger than the max energy
  if(signal.trace[n] > signal.base.TH){
    // If glitch filter and saturation filer are not active
    if (GLITCH_FILTER == false && SATURATION_FILTER == false){
      return(0);
    }
    // FILTER
    else{
    	// Glitch filter
    	//
    	if (GLITCH_FILTER == true){
        if (is_glitch(signal.trace, signal.base.TH, n)){
        	return(1);
        }
    	}
    	// Saturation filter
    	//
    	if (SATURATION_FILTER == true){
        if (is_saturation(signal, n)){
        	return(2);
        }
    	}
      // Saturation filter
      //
      if (baseline_weird(signal)){
        return(4);
      }
    	return(0);
    }
  }
  return(3);
}

// returns integral of a signal (but not over/under-shoot)
double signal_integral(signal_struct &signal, int debug = 0){
	//
	int left = 0;
	int right = 0;
	double integral = 0;
  double fraction = 0.1;
	// Find start of signal by looking for first sample above threshold left of the maximum
	for (int n = signal.energy_n; n > 0; n--){
		if (signal.trace[n] < fraction * signal.energy) {
			left = n;
			break;
		}
	}
  if (debug == 1) printf("%d\n", left);
	// If no sample above TH is found, return 0
	if (left == 0) return(0.0);
	// If a sample is found, look for the next sample below baseline (due to undershoot)
	for (int n = left; n < (int)signal.trace.size(); n++){
		// sum up the integral
		integral += signal.trace[n];
		// look for the zero crossing point after the maximum sample
		if (signal.trace[n] < 0 && n > signal.energy_n) {
			right = n;
			break;
		}
	}
	// if valid signal, return integral
  if (debug == 1) printf("%d %3.1f %d\n", signal.energy_n, signal.energy, right - left);
	if ( right - left > 0 ) return(integral);
	else { return(0.0); }
}

// Polynomial of nth degree, dimension given by par.size()
double polnx(double x, vector<double> &par){
  // Dimension
  int n = (int)par.size();
  // Return value
  double f = 0;
  // Execute sum of polynomials
  for (int i = 0; i < n; i++){
    f += par[i] * pow(x, i);
  }
  // Return f
  return(f);
}

// Logarithm with 3 parameters
double log3x(double x, vector<double> &par){
  // Dimension
  if ((int)par.size() != 3){
    printf("WARNING(log3x): Wrong amount of parameters passed (size: %d, must be 3!\n", (int)par.size());
    return(0);
  } 
  // Return value
  double f = par[0] - par[1] * log(x+par[2]);
  // Return f
  return(f);
}

// Exponential Growth
double ExpGro1(double x, vector<double> &par){
  // Dimension
  if ((int)par.size() != 3){
    printf("WARNING(ExpGro1): Wrong amount of parameters passed (size: %d, must be 3!\n", (int)par.size());
    return(0);
  } 
  // Return value
  double f = par[0] * exp( x / par[1] ) + par[2];
  // Return f
  return(f);
}

// Exponential Decay with undershoot compensation
double ExpDecay1(double x, vector<double> &par){
  // Dimension
  if ((int)par.size() != 5){
    printf("WARNING(ExpGro1): Wrong amount of parameters passed (size: %d, must be 55!\n", (int)par.size());
    return(0);
  } 
  // Return value
  double f = par[0] * exp( - x / par[1] ) + par[2] + par[3] * x + par[4] * pow(x, 2);
  // Return f
  return(f);
}

// Fit function for pixel saturation
Double_t SIPMpixel(Double_t *x, Double_t *par){
  // Parameter carriers
  // par[0]: corresponds to total amount of pixels, N_total
  // x: N of pixels fired 
  double value = par[0] * ( 1 - exp( - x[0] / par[0] ) );
  return value;
}

// Fit function for pixel saturation
Double_t SIPMlinearize(Double_t x, Double_t A){
  // Parameter carriers
  // par[0]: corresponds to total amount of pixels, N_total
  // x: N of pixels fired 
  double value = - A * log( 1 - ( x / A ) );
  return value;
}

// Energy  and time resolution fit
Double_t resolution(Double_t *x, Double_t *par){
  // Parameter carriers
  // par[0]: A is stochastic shower development parameter
  // par[1]: B is parameter for electric noice  
  // par[2]: C is free parameter
  // par[3]: D is third degree parameter
  // par[4]: x_0 is pole shift parameter
  // x: Energy
  double value = sqrt( (pow(par[0], 2)/(x[0]-par[4])) + (pow(par[1], 2)/pow((x[0]-par[4]),2)) + (pow(par[3], 3)/pow((x[0]-par[4]),3)) + pow(par[2], 2) );
  return value;
}

// Smooth array around passed value with smoothing interval given
vector<double> array_smooth(vector<double> &array, int s, int L){
  // copy array
  vector<double> returned = array;
  // do the smoothing around sample s
  for (int n = s-L; n < (s+L); n++){
    // adjust smoothing interval around current sample n to the distance to s
    // smoothing interval length
    int length = L;//- (int)( 0.5 *abs(s-n))+ L;
    int counter = 0;
    double sum = 0;
    for (int k = n-length ; k < n+length; k++){
      sum += array[k];
      counter++;
    }
    sum /= counter;
    returned[n] = sum;
  }
  return(returned);
}



vector<double> array_simulate_proto(){
  //
  //  Channel 13
  //
  int multiplier = 100;
  int junc1 = 108 * multiplier;
  int junc2 = 112 * multiplier;
  int junc3 = 127 * multiplier;
  // 
  vector<double> par1;
  par1.push_back(-10.9343);
  par1.push_back(0.22214);
  // 
  double value;
  //
  vector<double> ch12;
  for (int x = 0; x < junc1; x++){
    value = (double)(x) / (double)(multiplier);
    ch12.push_back(polnx(value, par1));
  }
  vector<double> par2;
  par2.push_back(-1308.47);
  par2.push_back(-1998.47);
  par2.push_back(-106.078);
  for (int x = junc1; x < junc2; x++){
    value = (double)(x) / (double)(multiplier);
    ch12.push_back(log3x(value, par2));
  }
  vector<double> par3;
  par3.push_back(-953116.94);
  par3.push_back(23054.8573);
  par3.push_back(-184.806);
  par3.push_back(0.49213);
  for (int x = junc2; x < junc3; x++){
    value = (double)(x) / (double)(multiplier);
    ch12.push_back(polnx(value, par3));
  }
  vector<double> par4;
  par4.push_back(31867.61);
  par4.push_back(-426.061);
  par4.push_back(1.755);
  par4.push_back(-0.001);
  par4.push_back(-8.7766e-6);
  par4.push_back(1.4731e-8);
  for (int x = junc3; x < 300*multiplier; x++){
    value = (double)(x) / (double)(multiplier);
    ch12.push_back(polnx(value, par4));
  }
  //
  // Smooth out the overlaps
  ch12 = array_smooth(ch12, junc1, multiplier/1);
  ch12 = array_smooth(ch12, junc2, multiplier/1);
  ch12 = array_smooth(ch12, junc3, multiplier/1);


  //
  //  CH07
  //
  // multiplier = 100;
  junc1 = 110 * multiplier;
  // int junc2 = 112 * multiplier;
  // int junc3 = 127 * multiplier;
  // 
  par1.clear();
  par1.push_back(2.1488e-112);
  par1.push_back(0.41476);
  par1.push_back(-4.13246);
  // 
  value = 0.0;
  //
  vector<double> ch07;
  for (int x = 0; x < junc1; x++){
    value = (double)(x) / (double)(multiplier);
    ch07.push_back(ExpGro1(value, par1));
  }
  par2.clear();
  par2.push_back(543803.86499);
  par2.push_back(330.92119);
  par2.push_back(-520005.51566);
  par2.push_back(1330.21254);
  par2.push_back(-1.10282);
  for (int x = junc1; x < 300*multiplier; x++){
    value = (double)(x) / (double)(multiplier);
    ch07.push_back(ExpDecay1(value, par2));
  }
  //
  // Smooth out the overlaps
  ch07 = array_smooth(ch07, junc1, multiplier/10);


  if (proto_out->is_open()){
    for (int n = 0; n < (int)ch07.size(); n++){
      *proto_out << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch07[n] << ";" 
                 << ch07[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n] << ";" 
                 << ch12[n]; 
      *proto_out << endl;
    }
    proto_out->close();
  }

  //
  hfile->cd("JUNK");
  plot_trace(ch07, "ch07", "AL*");
  plot_trace(ch12, "ch12", "AL*");
  return ch12;
}