// Sub-header file for aw_lukas
//
// Function definitions for data readout

// Plot functions
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

// Main "Physics" functions
void extraction();
void multis_calib();
void interpolate(vector<signal_struct> &signal);
void print_final_statistics();
void print_energy_statistics(vector<signal_struct> &array, const char *name);
void print_energy_calib();
void print_timing_statistics(time_struct &array, Int_t total_coincidents, const char *name);
void print_stat_multis_calib();
void print_detector_config();

// "Janitor" functions
double randit(int ini=0);
void print_usage();
void build_structure();
bool read_file(string file);
bool read_config(char const *file);

// Math and fit functions
double polnx(double x, vector<double> &par);
double log3x(double x, vector<double> &par);
double ExpDecay1(double x, vector<double> &par);
double ExpGro1(double x, vector<double> &par);
Double_t SIPMpixel(Double_t *x, Double_t *par);
Double_t SIPMlinearize(Double_t x, Double_t A);
Double_t resolution(Double_t *x, Double_t *par);
Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, bool silent);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
bool linreg(vector<double> &x, vector<double> &y, double *m, double *b);


// Array operations
double array_mean(const vector<double> &array, int start, int end);
double array_std(const vector<double> &array, int start, int end, double mean);
int array_largest(vector<double> &array, int lower, int upper);
int array_zero_xing(vector<double> &array, int lower, int upper);
double array_compare(vector<double> &array1, vector<double> &array2, vector<double> &weigths, vector<double> &par, int start, int end, int debug);
vector<double> array_adjust(vector<double> &array, vector<double> &x, int debug);
vector<double> array_sum(vector<double> &array1, vector<double> &array2, double factor);
vector<double> array_simulate_proto();
vector<double> array_smooth(vector<double> &array, int s, int L);
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
vector<Double_t> fit_graph_err(TGraphErrors *graph, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
Int_t largest_1Dbin(TH1D *hist, Int_t lower, Int_t upper);
vector<double> FIR_filter(vector<double> &trace, double calib);
vector<double> MA_filter(vector<double> &trace, double calib);
vector<double> MWD_filter(vector<double> &trace, double calib);
vector<double> CFD(vector<double> &trace, double multiplier);

// Program routines
void init_intersamp_hist(int channels);
void init_signal(vector<signal_struct> &signal, int channels, bool is_raw = false);
void init_times(vector<vector<time_struct> > &array, int channels);
void reset_signal(vector<signal_struct> &signal);
void reset_times(vector<vector<time_struct> > &array);
void init_multis_norm(vector<vector<vector<multis_norm_struct> > > &array, int channels);
void fill_hists();
void init_hists(int channels);

// Signal_struct routines
bool is_in_string(char const *character, char const *letter);
bool is_glitch(vector<double> &trace, double TH, int n);
bool is_saturation(signal_struct &signal, int n);
bool baseline_weird(signal_struct &signal);
int is_valid_max(signal_struct &signal, int n);
double signal_integral(signal_struct &signal, int debug);
void time_compare(vector<signal_struct> &signal);
