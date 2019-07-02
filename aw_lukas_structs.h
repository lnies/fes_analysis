// Sub-header file for aw_lukas
//
// Structs for data readout

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
  hist_struct_TH1D h_energy_mt_ring; // Histogram w/o multiples and w/ timing cut w/o the central crystal
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
