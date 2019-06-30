// Main header file for aw_lukas
//
// Include all other headers


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
// aw_lukas specific header
#include "simplex.h"
#include "aw_lukas_structs.h"
#include "aw_lukas_globals.h"
#include "aw_lukas_functions.h"
