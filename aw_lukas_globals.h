// Sub-header file for aw_lukas
//
// Globals for data readout


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
int ZERO_XING_CUT=70;
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

