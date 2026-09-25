#ifndef CRADLE_CONFIG_CONTAINER_H
#define CRADLE_CONFIG_CONTAINER_H

#include "CLI11.hpp"

#include <string>
#include <complex>

namespace CRADLE {

struct NuclearOptions {
  std::string Name = "";
  int Charge = 0;
  int Nucleons = 0;
  double Energy = 0.0;
};

struct General {
  int Verbosity = 1;
  int Verbosity_file = 1;
  int Loop = 0;
  int Threads = 8;
  std::string Output = "output.txt";
  int Seed = 0;
};

struct CouplingConstants {
  std::complex<double> CS = std::complex<double>(0.0,0.0);
  std::complex<double> CSP = std::complex<double>(0.0,0.0);
  std::complex<double> CV = std::complex<double>(1.0,0.0);
  std::complex<double> CVP = std::complex<double>(1.0,0.0);
  std::complex<double> CT = std::complex<double>(0.0,0.0);
  std::complex<double> CTP = std::complex<double>(0.0,0.0);
  std::complex<double> CA = std::complex<double>(1.2754,0.0);
  std::complex<double> CAP = std::complex<double>(1.2754,0.0);
  double a = std::nan("");
  double b = std::nan("");
};


struct FormFactors {
    float fb = 4.0; // by default 4.0
    float fc1 = 1.0; // by default 1.0
    float fd = 1.0 ; // by default 1.0
};


struct Cuts {
  double Distance = 1.E10;
  double Lifetime = 1.E40;
  double Energy = 1.E10;
};

struct BetaDecay {
  std::string Default = "Auto";
  std::string FermiFunction = "";
  std::string aCustom = "";
  double PolarisationX = 0;
  double PolarisationY = 0;
  double PolarisationZ = 0;
  double PolarisationMag = 0;
  double Alignment = 0;
  bool RadiativeCorrection = true ;
  double OmegaValue = 0.001 ;
  bool ElectronCapture = true ;
};

struct EnvOptions {
  std::string AMEdata="../Nuclear_Databases/AMEdata_.txt"; 
  std::string Gammadata="../GammaData"; 
  std::string Radiationdata="../RadiationData"; 
  std::string MixingRatiodata="../MixingRatio/MixingRatio_.txt";
};

struct ConfigOptions{
  NuclearOptions nuclearOptions;
  General general;
  CouplingConstants couplingConstants;
  FormFactors formfactors;
  Cuts cuts;
  BetaDecay betaDecay;
  EnvOptions envOptions;
};

ConfigOptions ParseOptions(std::string, int argc = 0, const char** argv = nullptr);

void SetCmdOptions(CLI::App&, NuclearOptions&);
void SetGeneralOptions(CLI::App&, General&);
void SetCouplingConstants(CLI::App&, CouplingConstants&);
void SetFormFactors(CLI::App&, FormFactors&);
void SetCuts(CLI::App&, Cuts&);
void SetBetaDecayOptions(CLI::App&, BetaDecay&);
void SetEnvironmentOptions(CLI::App&, EnvOptions&);

}

#endif
