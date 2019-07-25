// Sub-header file for aw_lukas
//
// Function definitions for data readout

// Plot functions
void plot_waves(vector<signal_struct> &array, char const *name, char const *modus);
void plot_waves_compare(char const *name, char const *modus);
void plot_interpol(vector<double> &x, vector<double> &y, double calib);
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
int array_zero_xing(vector<double> &array, int lower, int upper, int direction = 1);
double array_compare(vector<double> &array1, vector<double> &array2, vector<double> &weigths, vector<double> &par, int start, int end, int debug);
vector<double> array_adjust(vector<double> &array, vector<double> &x, int debug);
vector<double> array_sum(vector<double> &array1, vector<double> &array2, double factor);
vector<double> array_simulate_proto();
vector<double> array_smooth(vector<double> &array, int s, int L);
vector<Double_t> fit_hist(TH1D *hist, TF1 *fit, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
vector<Double_t> fit_graph_err(TGraphErrors *graph, char const *func, Double_t lower = 0, Double_t upper = 1, int verbose = 0);
Int_t largest_1Dbin(TH1D *hist, Int_t lower, Int_t upper);
vector<double> FIR_filter(vector<double> &trace, double calib);
vector<double> MA_filter(vector<double> &trace, double calib, int window);
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




/// BODY OF FUNCTIONS HERE
// Some functions still in main, have to be migrated later when done (or not)




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
        if (j == 2 || j == 3 || j == 4) continue; // dont print MA, MWD
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
        if (j == 2 || j == 3 || j == 4) continue; // dont print MA, MWD
        tg_combined[j] = new TGraph(tracelen,wave_x[j],wave_y[j]);
      }
      tg_combined[4] = new TGraph((int)NMO[i].trace.size(),wave_x_long,wave_y_long);
    } 
    else{ printf("ERROR (plot_waves): Plot option invalid\n");}
    for (int j = 0; j < 5; j++){
		if (j == 2 || j == 3 || j == 4) continue;
      if (j==0) {tg_combined[j]->SetLineColor(kBlack); tg_combined[j]->SetMarkerColor(kBlack);}
      if (j==1) {tg_combined[j]->SetLineColor(kGreen); tg_combined[j]->SetMarkerColor(kGreen);}
      if (j==2) {tg_combined[j]->SetLineColor(kOrange); tg_combined[j]->SetMarkerColor(kOrange);}
      if (j==3) {tg_combined[j]->SetLineColor(kBlue); tg_combined[j]->SetMarkerColor(kBlue);}
      if (j==4) {tg_combined[j]->SetLineColor(kRed); tg_combined[j]->SetMarkerColor(kRed);}      
      tg_combined[j]->SetMarkerSize(0.35);
      tg_combined[j]->SetLineWidth(2);
      char text[100];
      if (j==0) sprintf(text,"RAW_CALIB");
      if (j==1) sprintf(text,"MA");
      if (j==2) sprintf(text,"MWD");
      if (j==3) sprintf(text,"TMAX");
      if (j==4) sprintf(text,"NMO");
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
  TCanvas *c_waves = new TCanvas("c1","Wave_forms",10,10,1024,1024);
	gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptFit(1111);
  gStyle->SetLabelSize(0.05,"X");
  c_waves->SetGrid();
  TGraphErrors *tg_waves;
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetHeader(""); // option "C" allows to center the header
  // Calculate the energy steps
  double e_step = (E_WINDOW_RIGHT - E_WINDOW_LEFT)/N_E_WINDOW/100; // channels calculated to 100ch=1MeV
  Double_t wave_x[N_E_WINDOW];
  for(Int_t k = 0; k<N_E_WINDOW; k++){
    wave_x[k] = E_WINDOW_LEFT/100 + k*e_step;
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

  // tg_waves->RemovePoint(0); 
  // tg_waves->RemovePoint(0); 
  // tg_waves->RemovePoint(0); 

  tg_waves->SetName("data_graph");
  
  tg_waves->SetLineColor(kBlack);
  // tg_waves->SetMarkerColor(i+1);
  tg_waves->SetMarkerColor(kBlack);
  tg_waves->SetMarkerStyle(20);
  tg_waves->SetMarkerSize(3);
  tg_waves->SetLineWidth(2);

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);

  c_waves->SetTickx();
  c_waves->SetTicky();


  tg_waves->GetXaxis()->SetTitle("Energy [MeV]");
  tg_waves->GetYaxis()->SetTitle("Time resolution #sigma [ns]");
  tg_waves->SetTitle("Time resolution for different energy ranges");  

  tg_waves->Draw("ap");

  vector<Double_t> fit_params;

  fit_params = fit_graph_err(tg_waves, "resolution", 20, 30, 0);



  gPad->Modified();
  gPad->Update();
  c_waves->Write("Time_Resolution");

  delete c_waves;
}

void plot_interpol(vector<double> &x, vector<double> &y, double calib){
  TCanvas *c_waves = new TCanvas("c1","Wave_interpolation",200,10,500,300);
  TMarker *xing; 
  c_waves->SetGrid();
  TMultiGraph *mg_waves = new TMultiGraph();
  mg_waves->SetTitle("Interpolation exaple; Time [ns]; interpol. ADC channel [arb. unit]");
  TGraph *tg_waves;
  TGraph *tg_fits;
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
  int warning_counter = 0;
  legend->SetHeader("Interpolated ADC values"); // option "C" allows to center the header
  // Calibrate the x-values
  for (int n = 0; n < (int)x.size(); n++){
  	x[n] *= calib;
  }

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
  vector<double> fit_y_vector;
  for(Int_t n = 0; n<(Int_t)x.size(); n++){
    wave_x[n] = (Double_t) x[n];
    wave_y[n] = (Double_t) y[n];
    fit_y[n] = (Double_t) m * x[n] + b; // linear fit for wave
    fit_y_vector.push_back( (Double_t) m * x[n] + b ); // linear fit for wave
  }
  // Build TGraph
  tg_waves = new TGraph(x.size(),wave_x,wave_y);
  tg_waves->SetLineColor(1);
  tg_waves->SetMarkerColor(1);
  tg_waves->SetLineWidth(2);
  tg_waves->SetTitle("Channel i");
  tg_waves->SetMarkerStyle(kOpenSquare); // Asterisk
  tg_fits = new TGraph(x.size(),wave_x,fit_y);
  tg_fits->SetLineColor(kRed);
  tg_fits->SetMarkerColor(kRed);
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

  double xcrossing = -b/m;

  xing = new TMarker(xcrossing, 0, 1);
  xing->SetMarkerStyle(5);
  xing->SetMarkerColor(kRed);
  xing->SetMarkerSize(3);

  mg_waves->Draw("AL");
  xing->Draw();
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
      // array[i].h_energy.hist->SetFillColor(38);
      array[i].h_energy.hist->SetFillStyle(1001);
      array[i].h_integral.hist->SetLineWidth(3);
      array[i].h_integral.hist->SetLineColor(kBlack);
      array[i].h_energy.hist->GetXaxis()->SetNdivisions(20);
      array[i].h_energy.hist->GetXaxis()->SetRangeUser(0,20000);
      array[i].h_energy.hist->Draw();
    }
    if (strcmp(mode,"integral")==0){
      array[i].h_integral.hist->Rebin(1);
      array[i].h_integral.hist->SetFillColor(38);
      array[i].h_integral.hist->SetFillStyle(1001);
      array[i].h_integral.hist->SetLineWidth(3);
      array[i].h_integral.hist->SetLineColor(kBlack);      
      array[i].h_integral.hist->GetXaxis()->SetNdivisions(20);
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
  graph->GetYaxis()->SetTitle("\\frac{\\sigma}{E} \\left[%\\right]");

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);



  graph->Draw("ap");

  vector<Double_t> fit_params;

  fit_params = fit_graph_err(graph, "resolution", 10, 300, 0);

  // graph->GetFunction("resolution")->SetLineColor(kRed);
  // graph->GetFunction("resolution")->SetLineWidth(4);

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
int array_zero_xing(vector<double> &array, int lower, int upper, int direction){
	// Check boundaries
	if (lower < 0) lower = 0;
	if (upper > (int)array.size()) upper = (int)array.size(); 
	// For Forward direction:
	if ( direction > 0 ){
		// Loop over all indices
		for (int n = lower; n < upper; n++){
			//
		    if ( array[n-1] < 0 && array[n] > 0 ){
		    	return(n);
		    }
		}
	}
	// For Revers direction:
	if ( direction < 0 ){
		// Loop over all indices
		for (int n = upper; n > lower; n--){
			//
		    if ( array[n-1] < 0 && array[n] > 0 ){
		    	return(n);
		    }
		}
	}
	// if no xing is found, return 0
	return(0);
}

// Returns the sum of the weighted differences of each element of array1 and array2 between index start and index end
double array_compare(vector<double> &array1, vector<double> &array2, vector<double> &weigths, vector<double> &par, int start, int end, int debug){
  // Parameter
  double A = par[0]; // Multiplication
  double m = par[1]; // Shift
  // Some debug output
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
  if ( floor(m) != ceil(m) ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed position parameter par[1] is not an integer. Rounded to %d.\n", (int)round(m));
    m = (double)round(m);
  }

  // 
  double sum = 0.0;
  int n2 = 0; 
  // No everything should be of the correct size
  for (int n1 = start; n1 < end; n1++){
    // Using array1 as reference frame, calculate equivalent sample value for array2
    n2 = (int)( (int)(round(n1/conversion)) - (int)m );
    // 
    if (n2 > (int)array2.size()){ // If shift is out of range, give penalty back
      if (debug == 1) printf("WARNING(array_compare): Conversion went wrong, calculated sample n2=%d is larger than array2 size %d! Array value is %3.3f!\n", n2, (int)array2.size(), array2[n2]);
      sum += 100000;
    }
    // Do the sum
    else {
      sum += weigths[n1] * abs( array1[n1]- A * array2[n2] );
    }
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
  // Parameters:
  double A = x[0]; // Amplitude
  double m = x[1]; // Time shift
  // x[2] is the error value from array_compare, not being used here
  // for (int i = 0; i < (int)array.size(); i++) returned.push_back(0.0);
  vector<double> dummy ((int)array.size(), 0.0);
  // for (int i = 0; i < (int)array.size(); i++) dummy.push_back(0.0);
  // Check if passed parameters are healthy
  if ( (int)array.size() == 0 ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed array is empty. Return empty array.\n"); 
    returned.clear();
    return(returned);
  }
  if ( (int)x.size() > 3 ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed parameter array is of wrong dimension (has to be N=2 (or 3) ). Return empty array.\n"); 
    returned.clear();
    return(returned);
  }
  // Check position adjustment parameter is an "integer"
  if ( floor(m) != ceil(m) ){
    if (debug == 1 ) printf("WARNING(array_adjust): Passed position parameter x[1] is not an integer. Rounded to %d.\n", (int)round(x[1]));
    m = (double)round(m);
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
  	// Get mapping 
  	for (int a = 0; a<(int)MAPPING.size(); a++){
	    for (int b = 0; b<(int)MAPPING[a].size(); b++){
	      // Chrystal channel (Channel of matrix)
	      int ch = b*(int)MAPPING.size()+a; 
	      // if ch == i check the multiplicity
	      if (i == ch) signal[i].multis = MAPPING[a][b].multis;
	    }
	  }
    // Initialize the proto trace with 0s
    for ( int n = 0; n < TRACELEN*signal[i].multis; n++) {
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
    // In cosmics mode, initialize all energy fit parameters with 0
    if (strcmp(MODE, "COSMICS") == 0){
      for (int t = 0; t < 10; t++)
        signal[i].h_energy.params.push_back(1.0);
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
                      TRACELEN*RAW_CALIB[i].multis,0,TRACELEN*RAW_CALIB[i].multis, // x-dimension
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
          sprintf(name,"TMAX_PH_SUM_TAGGER_MT_RING_%02d", k);
          ECAL25[i].tagged[k].h_energy_mt_ring.hist=new TH1D(name,"",bins_pulsheight,0,range_pulsheight);
          ECAL25[i].tagged[k].h_energy_mt_ring.hist->GetXaxis()->SetTitle("Energy / MeV (NOT CALIB YET)");
          ECAL25[i].tagged[k].h_energy_mt_ring.hist->GetYaxis()->SetTitle("Counts");
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
  for (Int_t i = abs(lower); i < (upper+abs(lower)); i++){
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
  else if (strcmp(func, "langaus_roofit")==0){n = 6;}
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
    Double_t divide = hist->GetXaxis()->GetXmax() / hist->GetNbinsX();
    Double_t largest_bin = largest_1Dbin(hist, 2000*GENERAL_SCALING/divide, 80000*GENERAL_SCALING/nrebin)*divide*GENERAL_SCALING; // heaviest bin between lower and upper bound
    // printf("%3.1f\n", largest_bin);
    Double_t fr[2]; // fit boundaries
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4]; 
    fr[0]=0.65*largest_bin; // Lower fit boundary
    fr[1]=3.0*largest_bin; // Upper fit boundary
    printf("%3.3f %3.3f %3.3f %3.3f\n", largest_bin, fr[0], fr[1], divide);

    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //par[4]=A from A/(x^(m))
    //par[5]=m from A/(x^(m))
    //ADDED LATER: par[6]= Maximum of convoluted function
    //ADDED LATER: par[7]= FWHM of convoluted function
    pllo[0]=0.01      ; pllo[1]=0.               ; pllo[2]=1.0          ; pllo[3]=1.     ;// pllo[4]=-10.0 ; pllo[5]=0.0001;  // Lower parameter limits
    plhi[0]=20000.   ; plhi[1]=largest_bin+20000.; plhi[2]=1e12; plhi[3]=10000.0;// plhi[4]=100000.0; plhi[5]=5.0; // Upper parameter limits
    sv[0]  =200.    ; sv[1]  =largest_bin      ; sv[2]  =5000000.0    ; sv[3]  =1000.0 ;// sv[4]  =100.0   ; sv[5]=0.05;// Start values
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

  // Stable langaus with RooFit 
  if (strcmp(func, "langaus_roofit")==0){
    //
    // Look for the maximum value
    hfile->cd("JUNK");
    int nrebin = 1;
    Double_t divide = hist->GetXaxis()->GetXmax() / hist->GetNbinsX();  
    Double_t maxBin = largest_1Dbin(hist, 2000*GENERAL_SCALING/divide, 80000*GENERAL_SCALING/nrebin)*divide*GENERAL_SCALING; // heaviest bin between lower and upper bound
    printf("+++++++++++++++++++ MAXBIN: %3.3f\n", maxBin);
    Double_t minX = 0.75 * maxBin;
    Double_t maxX = 3.00 * maxBin;
    Double_t leftX = 0.75 * maxBin;
    Double_t rightX = 1.25 * maxBin;
    //assuming your histogram is the variable "hist".
    Double_t sigma = (rightX-leftX)/2.35;
    // Construct observable
    // RooRealVar t("t","t",minX,maxX);
    RooRealVar t("t","t",hist->GetXaxis()->GetBinLowEdge(1),hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));
    // Define fit range
    t.setRange("ROI_1",minX,maxX);
    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg","mg",0) ;
    RooRealVar sg("sg","sg",sigma,0.1*sigma,5.*sigma) ;
    RooGaussian gauss("gauss","gauss",t,mg,sg) ;

    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml","mean landau",maxBin,maxBin-sigma,maxBin+sigma) ;
    RooRealVar sl("sl","sigma landau",sigma,0.1*sigma,5.*sigma) ;
    // RooRealVar sl("sl","sigma landau",0.04,0.,0.2) ;
    RooLandau landau("lx","lx",t,ml,sl) ;

    // C o n s t r u c t   c o n v o l u t i o n   p d f 
    // ---------------------------------------

    // Set #bins to be used for FFT sampling
    t.setBins(5000,"cache") ; 

    // Construct landau (x) gauss
    RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;

    // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
    // ----------------------------------------------------------------------

    RooDataHist* data = new RooDataHist("dh","dh",t,RooFit::Import(*hist)) ;

    // Fit gxlx to data
    lxg.fitTo(*data,RooFit::Range("ROI_1"));
    
    // Plot data, landau pdf, landau (X) gauss pdf
    TCanvas *canvas;
    canvas = new TCanvas("Langaus","Langaus",800,800);
    RooPlot* frame = t.frame(RooFit::Title((TString)"FitProjection")) ;
    // data->plotOn(frame) ;
    RooPlot* fitLine=lxg.plotOn(frame) ;

    lxg.paramOn(frame);
    data->statOn(frame);
    landau.plotOn(frame,RooFit::LineStyle(kDashed)) ;
    // gauss.plotOn(frame,RooFit::LineStyle(kDashed)) ;

    TF1* flxg = lxg.asTF (RooArgList(t)) ;//NOT NORMALIZED
      
    Double_t modX=flxg->GetMaximumX();
    
    RooCurve* fitCurve = (RooCurve*)fitLine->findObject("lxg_Norm[t]");
    // RooCurve* fitCurve = (RooCurve*)fitLine->findObject("lxg_Norm[t]_Range[fit_nll_lxg_dh]_NormRange[fit_nll_lxg_dh]");
    Double_t maxpico=fitCurve?fitCurve->getYAxisMax():0; 
    Double_t maxpicoU=flxg->Eval(modX);//unnormalized
    Double_t extrems[2];
    
    extrems[0]=flxg->GetX(maxpicoU/2.0,minX,modX);
    extrems[1]=flxg->GetX(maxpicoU/2.0,modX,maxX);

    printf("%3.3f %3.3f %3.3f %3.3f\n", maxpico, maxpicoU, extrems[0], extrems[1]);
    
    // Draw frame on canvas
    frame->Draw();

    TH1F *clone = (TH1F*)(hist->Clone("clone"));
    Double_t norm = clone->GetEntries();
    clone->Scale(1/clone->GetBinContent(maxBin/100));
    printf("\n\n\n+++++ BINCONTENT: %3.3f\n\n\n\n\n", clone->GetBinContent(maxBin));

    clone->Draw("SAME");

    canvas->Write("RooFit_Langaus");
  
    TLine* l;
    l = new TLine(extrems[0],maxpico/2.0,extrems[1],maxpico/2.0);
    l->SetLineWidth(1);
    l->SetLineColor(1);
    l->Draw("same");

    l = new TLine(extrems[0],hist->GetMinimum(),extrems[1],hist->GetMinimum());
    l->SetLineWidth(1);
    l->SetLineColor(1);
    l->Draw("same");
    
    
    l = new TLine(modX,hist->GetMinimum(),modX,maxpico);
    l->SetLineWidth(0.5);
    l->SetLineColor(1);
    l->Draw("same");
    l=NULL;

    delete canvas;



    return(params);


  }

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
    par_low[0] = 0.00001; par_start[0] = 430; par_up[0] = 1000;
    par_low[1] = 0.00001; par_start[1] = 2200; par_up[1] = 50000;
    par_low[2] = 0.00001; par_start[2] = 15; par_up[2] = 50;
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
vector<double> MA_filter(vector<double> &trace, double calib, int window){
  // Return array
  vector<double> filtered;
  // Check if window is odd
  if ((window % 2) == 0) window += 1;
  double value = 0;
  // For each sample of the trac do:
	for(int n=0; n < (int)trace.size(); n++){
		value = 0;
		// Take care of boundary conditions
		// Left boundary
	  if( n - ((window-1)/2) < 1 ){
	    value = array_mean(trace, 0, n + (window-1)/2);
	  } 
	  // Center values
	  if( n - ((window-1)/2) >= 1 && n + ((window-1)/2) <= (int)trace.size() ){ 
	    value = array_mean(trace, n - ((window-1)/2), n + ((window-1)/2));
	  }
	  // Right boundary
	  if( n + ((window-1)/2) > (int)trace.size() ){ 
	    value = array_mean(trace, n - ((window-1)/2), (int)trace.size());
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
  hfile->mkdir("WAVE_FORMS/CFD/RAW_CALIB/INTERPOL");
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
        if(strcmp(key.c_str(), "ZERO_XING_CUT") == 0) ZERO_XING_CUT = stoi(value);
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
    if ( trace[k] > TH){ // Either test if above TH or above 0, depending on noise level
    // if ( trace[k] > 0){
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
	int rising_start = 100; // powerful, but use carefully
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
    // 4: Baseline is weird
    // 5: Maximum not in energy range
int is_valid_max(signal_struct &signal, int n){
	// Test if the current sample is larger than the max energy
  if(signal.trace[n] > signal.base.TH){
    // maximum must be in the specified energy window
    if ( ENERGY_WINDOW_MAX*signal.multis < signal.energy_n || BASELINE_CUT*signal.multis > signal.energy_n ){
      return(5); 
    }
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