#include "CRADLE/DecayManager.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/DecayChannel.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/SpectrumGenerator.hh"
#include "CRADLE/RadiativeCorrections.hh"
#include "CRADLE/Polarisation.hh"
#include "CRADLE/ECshell.hh"

// ROOT headers for .root output format
#include "TFile.h"
#include "TTree.h"

#define BOOST_TIMER_ENABLE_DEPRECATED
#include <boost/progress.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <future>
#include <complex>
#include <chrono>
#include <iomanip>
// #include <iostream>

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}

namespace CRADLE {

using std::map;
using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::endl;

DecayManager::~DecayManager() {
  // cout << "Destroying decaymanager" << endl;
  /*for (vector<Particle*>::iterator it = particleStack.begin();
       it != particleStack.end(); ++it) {
    delete *it;
  }*/
  for (map<const string, Particle*>::iterator it = registeredParticles.begin();
       it != registeredParticles.end(); ++it) {
    for (vector<DecayChannel*>::iterator it2 =
             (it->second)->GetDecayChannels().begin();
         it2 != (it->second)->GetDecayChannels().end(); ++it2) {
      delete *it2;
    }
    delete it->second;
  }
  registeredParticles.clear();

  for (map<const string, vector<vector<double> >*>::iterator it =
           registeredDistributions.begin();
       it != registeredDistributions.end(); ++it) {
    delete it->second;
  }
  registeredDistributions.clear();
}

void DecayManager::RegisterDecayMode(const string name, DecayMode& dm) {
  registeredDecayModes.insert(pair<string, DecayMode&>(name, dm));
  if (configOptions.general.Verbosity > 0)
    cout << "Registered DecayMode " << name << endl;
}

DecayMode& DecayManager::GetDecayMode(const string name) {
  if (registeredDecayModes.count(name) == 0) {
    throw std::invalid_argument("DecayMode " + name + " not registered. Aborting.");
  }
  return registeredDecayModes.at(name);
}

void DecayManager::RegisterParticle(Particle* p) {
  registeredParticles.insert(pair<string, Particle*>(p->GetRawName(), p));
  // cout << "Registered particle " << p->GetRawName() << endl;
}

Particle* DecayManager::GetNewParticle(const string name, int Z, int A) {
  if (registeredParticles.count(name) == 0) {
    GenerateNucleus(name, Z, A);
  }
  return new Particle(*(registeredParticles.at(name)));
}

void DecayManager::RegisterDistribution(const string name,
                                        vector<vector<double> >* dist) {
  registeredDistributions.insert(
      pair<string, vector<vector<double> >*>(name, dist));
  // cout << "Registered distribution " << name << endl;
}

vector<vector<double> >* DecayManager::GetDistribution(const string name) {
  if (registeredDistributions.count(name) == 0) {
    throw std::invalid_argument("Distribution not registered.");
  }
  return registeredDistributions.at(name);
}

 /////// ajout de SL 12/05/2023//////////////////////////////////////
void DecayManager::RegisterBetaType(const string name,
                                        const string nameType) {
  registeredBetaType.insert(
      pair<string, string>(name, nameType));
  // cout << "Registered BetaType " << name << endl;
}

string DecayManager::GetBetaType(const string name) {
  if (registeredBetaType.count(name) == 0) {
    throw std::invalid_argument("Transition not registered.");
  }
  return registeredBetaType.at(name);
}


void DecayManager::RegisterParameterMC(const string name, std::vector<double> W) {
  registeredParameterMC.insert(pair<string, std::vector<double> > (name, W)) ;  
}

std::vector<double> DecayManager::GetParameterMC(const string name) {
  if (registeredParameterMC.count(name)==0) {
    throw std::invalid_argument("not registered.");
  }
  return registeredParameterMC.at(name) ;
}
/////////////////////////////////////////////////////////////////////

void DecayManager::RegisterBasicParticles() {
  RegisterParticle(new Particle("e-", utilities::EMASSC2, -1, 0, 0.5, 0.));
  RegisterParticle(new Particle("e+", utilities::EMASSC2, 1, 0, 0.5, 0.));
  RegisterParticle(new Particle("p", utilities::PMASSC2, 1, 0, 0.5, 0.));
  RegisterParticle(new Particle("n", utilities::NMASSC2, 0, 1, 0.5, 0.));
  RegisterParticle(new Particle("alpha", utilities::ALPHAMASSC2, 2, 2, 0., 0.));
  RegisterParticle(new Particle("enu", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("enubar", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("gamma", 0., 0, 0, 0., 0.));

  RegisterParticle(new Particle("e-RC", utilities::EMASSC2, -1, 0, 0.5, 0.));
  RegisterParticle(new Particle("e+RC", utilities::EMASSC2, 1, 0, 0.5, 0.));
  RegisterParticle(new Particle("enuRC", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("enubarRC", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("gammaBR", 0., 0, 0, 0., 0.));

  RegisterParticle(new Particle("enuEC", 0., 0, 0, 0.5, 0.));

  // For Xavier simulations
  RegisterParticle(new Particle("e-90Y", utilities::EMASSC2, -1, 0, 0.5, 0.));
  RegisterParticle(new Particle("enubar90Y", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("gammaBR90Y", 0., 0, 0, 0., 0.));

  RegisterParticle(new Particle("e-90Zr", utilities::EMASSC2, -1, 0, 0.5, 0.));
  RegisterParticle(new Particle("enubar90Zr", 0., 0, 0, 0.5, 0.));
  RegisterParticle(new Particle("gammaBR90Zr", 0., 0, 0, 0., 0.));
}

void DecayManager::RegisterBasicDecayModes() {
  RegisterDecayMode("BetaMinus", BetaMinus::GetInstance());
  RegisterDecayMode("BetaPlus", BetaPlus::GetInstance());
  // RegisterDecayMode("ConversionElectron", ConversionElectron::GetInstance());
  RegisterDecayMode("Proton", Proton::GetInstance());
  RegisterDecayMode("Alpha", Alpha::GetInstance());
  RegisterDecayMode("Gamma", Gamma::GetInstance());
  RegisterDecayMode("IT", Gamma::GetInstance());
  RegisterDecayMode("BetaMinusRadiative", BetaMinusRadiative::GetInstance());
  RegisterDecayMode("BetaPlusRadiative", BetaPlusRadiative::GetInstance());
  RegisterDecayMode("BetaMinusVirtualSoft", BetaMinusVirtualSoft::GetInstance());
  RegisterDecayMode("BetaPlusVirtualSoft", BetaPlusVirtualSoft::GetInstance());
  RegisterDecayMode("BetaMinusPolarised", BetaMinusPolarised::GetInstance());
  RegisterDecayMode("BetaPlusPolarised", BetaPlusPolarised::GetInstance());
  RegisterDecayMode("KshellEC", ShellEC::GetInstance());
  RegisterDecayMode("LshellEC", ShellEC::GetInstance());
  RegisterDecayMode("MshellEC", ShellEC::GetInstance());
}

void DecayManager::RegisterSpectrumGenerator(const string decayMode, SpectrumGenerator& sg) {
  try {
    DecayMode& dm = GetDecayMode(decayMode);
    dm.SetSpectrumGenerator(&sg);
    if (configOptions.general.Verbosity > 0)
      cout << "Registered " << decayMode << " Spectrum Generator " << typeid(sg).name() << endl;
  }
  catch (const std::invalid_argument &e) {
    cout << "Cannot register" <<  typeid(sg).name() << "spectrum generator. Decay mode "
    << decayMode << " not registered." << endl;
  }
}

void DecayManager::RegisterBasicSpectrumGenerators() {
  RegisterSpectrumGenerator("Proton", DeltaSpectrumGenerator::GetInstance());
  RegisterSpectrumGenerator("Alpha", DeltaSpectrumGenerator::GetInstance());
  RegisterSpectrumGenerator("Gamma", DeltaSpectrumGenerator::GetInstance());
  RegisterSpectrumGenerator("IT", DeltaSpectrumGenerator::GetInstance());
  RegisterSpectrumGenerator("BetaPlus", SimpleBetaDecay::GetInstance());
  RegisterSpectrumGenerator("BetaMinus", SimpleBetaDecay::GetInstance());
  RegisterSpectrumGenerator("BetaPlusPolarised", SimpleBetaDecay::GetInstance()); //remove once proper generator implemented
  RegisterSpectrumGenerator("BetaMinusPolarised", SimpleBetaDecay::GetInstance()); //remove once proper generator implemented
  RegisterSpectrumGenerator("ShellEC", DeltaSpectrumGenerator::GetInstance());
}

void DecayManager::ListRegisteredParticles() {
  cout << "--------------------------------------------------------\n";
  cout << " List of registered particles\n";
  cout << "--------------------------------------------------------\n\n";
  for (map<const string, Particle*>::iterator it = registeredParticles.begin();
       it != registeredParticles.end(); ++it) {
    cout << it->second->ListInformation();
    cout << "\n";
  }
  cout << "--------------------------------------------------------\n\n"
       << endl;
}

bool DecayManager::GenerateNucleus(string name, int Z, int A) {
  std::ostringstream filename;
  filename << configOptions.envOptions.Radiationdata;
  filename << "/z" << Z << ".a" << A;
  std::ifstream radDataFile((filename.str()).c_str());

  //cout << "Generating nucleus " << name << endl;

  string line;
  double excitationEnergy = 0.;
  double lifetime;
  double atomicMass = utilities::GetAMEMass(configOptions.envOptions.AMEdata, Z, A);

  if (atomicMass == 0) {
    atomicMass = utilities::GetApproximateMass(Z, A);
  }

  Particle* p = new Particle(name, atomicMass, Z, (A - Z), 0., 0);

  //cout << filename.str() << endl;

  while (getline(radDataFile, line)) {
    //cout<<line<<endl;

    if (!line.compare(0, 1, "#")) {
      // Comment line
      continue;
    } else if (!line.compare(0, 1, "P")) {
      // Parent line
      std::istringstream iss(line);
      string p;
      string flag;

      iss >> p >> excitationEnergy >> flag >> lifetime;
      continue;

    }
    //cout << "Lifetime: " << lifetime << endl;
    string mode;
    double daughterExcitationEnergy = 0.;
    double intensity = 0;
    double Q = 0.;
    string modifier;
    string flag;

    std::istringstream iss(line);
    iss >> mode >> daughterExcitationEnergy >> flag >> intensity >> Q >> modifier;
    //cout << "Mode : " << mode << endl;
    //cout << "Daughter Energy :" << daughterExcitationEnergy << endl;
    //cout << "Flag :" << flag <<endl;
    //cout << "Intensity : " << intensity << endl;
    //cout << "Q : " << Q << endl;
    //cout << "Modifier : " << modifier << endl;
    //cout << "\n" <<endl;

    //bool PolarisedNuclei = false;
    //if ((configOptions.betaDecay.PolarisationX != 0 || configOptions.betaDecay.PolarisationY != 0 || configOptions.betaDecay.PolarisationZ != 0) && configOptions.betaDecay.PolarisationMag != 0)
      //PolarisedNuclei = true;

    //NEW CONDITION TO HAVE POLARISATION AXIS WITH ALIGNMENT WO POLARISATION MAG
    bool PolarisedNuclei = false;

    bool hasPolarisationAxis =
    (configOptions.betaDecay.PolarisationX != 0 ||
    configOptions.betaDecay.PolarisationY != 0 ||
    configOptions.betaDecay.PolarisationZ != 0);

    bool hasPolarisation =
    (configOptions.betaDecay.PolarisationMag != 0);

    bool hasAlignment =
    (configOptions.betaDecay.Alignment != 0);

    if (hasPolarisationAxis && (hasPolarisation || hasAlignment)) {
      PolarisedNuclei = true;
    }


    if (Q > 0. ) {
      //std::cout << "Q > 0" << "\n" ;
      /*cout << "Adding DecayChannel " << mode << " Excitation Energy " <<
      excitationEnergy << " to " << daughterExcitationEnergy << endl;*/
      if ((mode.find("shellEC") != string::npos) && (configOptions.betaDecay.ElectronCapture == true)) {
        cout << "Mode : " << mode << endl;
        std::cout << "Intensity : " << intensity << "\n" ;
        
        int Z = p->GetCharge() ;
        double E_exc = 0. ;
        int NbShells = ecshell::fNumberOfShells[Z] ;
        int IndexShell = 0 ; 
        for (int i=0 ; i < Z-1 ; i++) {
          IndexShell += ecshell::fNumberOfShells[i+1] ;
        }
        //std::cout << "Z : " << Z << "\n";
        //std::cout << "Nb shell : " << NbShells << "\n";
        //std::cout << "index shell : " << IndexShell << "\n";
        //std::cout << "index shell : " << IndexShell << "\n";

        if (mode == "KshellEC") {
          E_exc = ecshell::fBindingEnergies[IndexShell + 1] * std::pow(10, -3) ; 
          //excitationEnergy = E_exc ;
          //excitationEnergy = 0 ;
          std::cout << "E_exc : " << E_exc << "\n" ;
          DecayChannel* dcK = new DecayChannel("ShellEC", &GetDecayMode(mode), Q - E_exc, intensity, lifetime, excitationEnergy, daughterExcitationEnergy);
          p->AddDecayChannel(dcK) ;
          std::cout << "\n" ;

        } else if (mode == "LshellEC") {
          E_exc = ecshell::fBindingEnergies[IndexShell + 2] * std::pow(10, -3) ;
          //excitationEnergy = E_exc ;
          //excitationEnergy = 0 ;
          std::cout << "E_exc : " << E_exc << "\n" ;
          double PL1 = ecshell::ProbabilityL1(Z) * 1;
          std::cout << "PL1 : " << PL1 << "\n" ;
          DecayChannel* dcL1 = new DecayChannel("ShellEC", &GetDecayMode(mode), Q - E_exc, intensity*PL1, lifetime, excitationEnergy, daughterExcitationEnergy);
          p->AddDecayChannel(dcL1) ;

          if (NbShells > 2) {
            E_exc = ecshell::fBindingEnergies[IndexShell + 3] * std::pow(10, -3) ;
            //excitationEnergy = E_exc ;
            //excitationEnergy = 0 ;
            std::cout << "E_exc : " << E_exc << "\n" ;
            double PL2 = 1 - PL1 ;
            std::cout << "PL2 : " << PL2 << "\n" ;
            DecayChannel* dcL2 = new DecayChannel("ShellEC", &GetDecayMode(mode), Q- E_exc, intensity*PL2, lifetime, excitationEnergy, daughterExcitationEnergy);
            p->AddDecayChannel(dcL2) ;
          }
          std::cout << "\n" ;

        } else if (mode == "MshellEC") {
          E_exc = ecshell::fBindingEnergies[IndexShell + 5] * std::pow(10, -3) ;
          //excitationEnergy = E_exc ; 
          std::cout << "E_exc : " << E_exc << "\n" ;
          double PM1 = ecshell::ProbabilityM1(Z) * 1 ;
          std::cout << "PM1 : " << PM1 << "\n" ;
          DecayChannel* dcM1 = new DecayChannel("ShellEC", &GetDecayMode(mode), Q- E_exc, intensity*PM1, lifetime, excitationEnergy, daughterExcitationEnergy);
          p->AddDecayChannel(dcM1) ;

          if (NbShells > 4) {
            E_exc = ecshell::fBindingEnergies[IndexShell + 6] * std::pow(10, -3) ;
            //excitationEnergy = E_exc ;
            std::cout << "E_exc : " << E_exc << "\n" ;
            double PM2 = 1 - PM1 ;
            std::cout << "PM2 : " << PM2 << "\n" ;
            DecayChannel* dcM2 = new DecayChannel("ShellEC", &GetDecayMode(mode), Q - E_exc, intensity*PM2, lifetime, excitationEnergy, daughterExcitationEnergy);
            p->AddDecayChannel(dcM2) ;
          }
          std::cout << "\n" ;

        }
        
      }
      else if (configOptions.betaDecay.RadiativeCorrection == false && (PolarisedNuclei == false) && ((mode.find("BetaPlus")!= std::string::npos)||(mode.find("BetaMinus")!= std::string::npos))) {
        std::cout << "RC false" << "\n" ;
        std::cout << "Intensity : " << intensity << "\n" ;
        std::cout << "\n" ;
        DecayChannel* dc = new DecayChannel(mode, &GetDecayMode(mode), Q, intensity, lifetime, excitationEnergy, daughterExcitationEnergy);
        p->AddDecayChannel(dc) ;
      }
      else if (configOptions.betaDecay.RadiativeCorrection == false && (PolarisedNuclei == true) && ((mode.find("BetaPlus")!= std::string::npos)||(mode.find("BetaMinus")!= std::string::npos) )) {
        std::cout << "RC false" << "\n" ;
        std::cout << "Intensity : " << intensity << "\n" ;

	std::complex<double> CS = configOptions.couplingConstants.CS;
        std::complex<double> CSP = configOptions.couplingConstants.CSP;
        std::complex<double> CV = configOptions.couplingConstants.CV;
        std::complex<double> CVP = configOptions.couplingConstants.CVP;
        std::complex<double> CA = configOptions.couplingConstants.CA;
        std::complex<double> CAP = configOptions.couplingConstants.CAP;
        std::complex<double> CT = configOptions.couplingConstants.CT;
        std::complex<double> CTP = configOptions.couplingConstants.CTP;

	double mf = 0.;
	double mgt = 0.;

	int Z = p->GetCharge() ;
        int A = p->GetCharge() + p->GetNeutrons(); 
        int betaType = Z/std::abs(Z);
  
        if (mode == "BetaPlus") {
          betaType = -1 ;
        }
	
	double j_i = utilities::GetJpi(A,Z,excitationEnergy);
        j_i = std::abs(j_i); //polarity not needed, only absolute value of J
        
        int Z_d = Z + betaType;
        double j_f = utilities::GetJpi(A,Z_d,daughterExcitationEnergy);
        j_f = std::abs(j_f);


	if (configOptions.betaDecay.Default == "Fermi") {
	  mf = 1.;
	}
	else if (configOptions.betaDecay.Default == "Gamow-Teller") {
	  mgt = 1.;
	}
	else if (configOptions.betaDecay.Default == "Mixed") {
	  mf = 1. ;
	  mgt = 1. ;
	}
	else if (configOptions.betaDecay.Default == "Auto") {
	  std::string Type = utilities::FindBetaType(A, Z, Z_d, excitationEnergy, daughterExcitationEnergy);
	  if (Type == "Fermi") {
	    mf = 1.;
	  } else if (Type == "Gamow-Teller") {
	    mgt = 1.;
	  } else { //Mixed
	    mgt = std::stod(Type.substr(5))/1.2754;
	    mf = 1.;
	  }
	}

	if (configOptions.general.Verbosity > 0){
	  std::cout << "Z : " << Z << " A : " << A << " J_in : " << j_i << " Decay Type: " << betaType << " Q : " << Q << std::endl;
	  std::cout << "Z : " << Z_d << " A : " << A << " J_f : " << j_f << " Level Energy: " << daughterExcitationEnergy << std::endl;
	  std::cout << "M_GT: " <<  mgt << ", M_F: " << mf << std::endl;
        }

	double xi = utilities::CalculateXiBetaDecay(CS,CSP,CT,CTP,CV,CVP,CA,CAP,mf,mgt);

	if (xi != 0){
	  DecayChannel* dc = new DecayChannel(mode+"Polarised", &GetDecayMode(mode+"Polarised"), Q, intensity, lifetime, excitationEnergy, daughterExcitationEnergy);
	  p->AddDecayChannel(dc);
	  if (configOptions.general.Verbosity > 0){
	    //std::cout << CT << std::endl;
	    double a_conf = configOptions.couplingConstants.a;
	    double b_conf = configOptions.couplingConstants.b;
	    double E = Q/3 + utilities::EMASSC2;
	    double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, E, Z_d, betaType);
	    double b = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, Z_d, betaType);
	    double c = 0;
	    double A = 0;
	    double B = 0;
	    double D = 0;
	    if (j_i > 0){
	      A = polarisation::CalculateBetaAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, betaType, Z_d, E);
	      B = polarisation::CalculateNeutrinoAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, betaType, Z_d, E);
	      D = polarisation::CalculateDTripleCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, betaType, Z_d, E);
	      if (j_i > 0.5){
		c = polarisation::CalculateAlignmentCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, betaType, Z_d, E);
	      }
	    }
	    std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", A: " << A << ", B: " << B << ", D: " << D << std::endl;
	  }
	} else {
	  std::cout << "The current choice of coupling constants makes this decay impossible" << std::endl;
	}
      } 
      else if (configOptions.betaDecay.RadiativeCorrection == true && ((mode.find("BetaPlus")!= std::string::npos)||(mode.find("BetaMinus")!= std::string::npos) )) {
        std::cout << "RC true" << "\n" ;
        double mass_i = p->GetMass() ;
        double mass_f = mass_i - Q ;
        DecayManager& dm = DecayManager::GetInstance();
        double mf = 0.;
        double mgt = 0.;
        if (dm.configOptions.betaDecay.Default == "Fermi") {
          mf = 1.;
        } 
        else {
          mgt = 1. ;
        }
        double CS = configOptions.couplingConstants.CS.real();
        double CSP = configOptions.couplingConstants.CSP.real();
        double CV = configOptions.couplingConstants.CV.real();
        double CVP = configOptions.couplingConstants.CVP.real();
        double CA = configOptions.couplingConstants.CA.real();
        double CAP = configOptions.couplingConstants.CAP.real();
        double CT = configOptions.couplingConstants.CT.real();
        double CTP = configOptions.couplingConstants.CTP.real();
        double a_conf = configOptions.couplingConstants.a;
        double b_conf = configOptions.couplingConstants.b;

        double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
        
        double Cs = configOptions.betaDecay.OmegaValue ;
        int Z = p->GetCharge() ;
        int A = p->GetCharge() + p->GetNeutrons() ; 
        double R = utilities::ApproximateRadius(A) ;
        int betaType = Z/std::abs(Z);
        if (mode == "BetaPlus") {
          betaType = -1 ;
          Q = Q - 2*utilities::EMASSC2 ;
        }

        double Jpi_init = utilities::GetJpi(A, Z, 0);
        double Jpi_final = utilities::GetJpi(A, Z + betaType, daughterExcitationEnergy);
        int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));
        std::cout << "Jpi init : " << Jpi_init << "\t Jpi final : " << Jpi_final << "\t Labs : " << Labs << std::endl;

        bool advanced = false ;
        if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
            advanced = true;
          }
        
        double ph = radiativecorrections::PH(Cs, mf, mgt, a, mass_i, std::abs(Z), A, Q, Labs, advanced, betaType) ;

        if (std::isnan(ph)) {
          std::cout << "Intensity : " << intensity << "\n";
        } else {

        //DecayChannel* dc = new DecayChannel(mode, &GetDecayMode(mode), Q, intensity, lifetime, excitationEnergy, daughterExcitationEnergy);
        //p->AddDecayChannel(dc) ;

        std::cout << mode+"Radiative" << "\n";
        std::cout << "Intensity : " << intensity << "\n";
        std::cout << "PH (en %) : " << ph*100 << "\n" ;
        std::cout << "Q : " << Q/utilities::EMASSC2 + 1 << "\n" ;
        std::cout << "\n" ;
        
        DecayChannel* dc1 = new DecayChannel(mode+"VirtualSoft", &GetDecayMode(mode+"VirtualSoft"), Q, intensity*(1.-ph), lifetime, excitationEnergy, daughterExcitationEnergy);
        p->AddDecayChannel(dc1) ;

        DecayChannel* dc2 = new DecayChannel(mode+"Radiative", &GetDecayMode(mode+"Radiative"), Q, intensity*ph, lifetime, excitationEnergy, daughterExcitationEnergy);
        p->AddDecayChannel(dc2) ;
        }
      } 
      else if (mode == "Proton") {
        std::cout << mode << "\n";
        std::cout << "Intensity : " << intensity << "\n";
        std::cout << "\n" ;
        
        DecayChannel* dcp = new DecayChannel("Proton", &GetDecayMode(mode), Q, intensity, lifetime, excitationEnergy, daughterExcitationEnergy);
        p->AddDecayChannel(dcp);
      }
    } 
  }

  std::ostringstream gammaFileSS;
  gammaFileSS << configOptions.envOptions.Gammadata;
  gammaFileSS << "z" << Z << ".a" << A;
  std::ifstream gammaDataFile(gammaFileSS.str().c_str());
  if (gammaDataFile.is_open()) {
    //std::cout << "here 1" << "\n";
    while (getline(gammaDataFile, line)) {
      //std::cout << "here 2" << "\n";
      int levelNr;
      double initEnergy, E;
      double intensity;
      double convIntensity;
      double kCoeff, lCoeff1, lCoeff2, lCoeff3, mCoeff1, mCoeff2, mCoeff3, mCoeff4, mCoeff5;
      double lifetime;
      string angMom;
      string polarity;
      string flag;

      int nGammas;

      std::istringstream iss(line);
      iss >> levelNr >> flag >> initEnergy >> lifetime >> angMom >> nGammas;
      for (int i = 0; i < nGammas; ++i) {
        getline(gammaDataFile, line);
        int daughterLevelNr;
        int multipolarity;
        double multipolarityMixing;
	std::istringstream issLevel(line);
        issLevel >> daughterLevelNr >> E >> intensity >> multipolarity >> multipolarityMixing >> convIntensity >> kCoeff >> lCoeff1 >> lCoeff2 >> lCoeff3 >> mCoeff1 >> mCoeff2 >> mCoeff3 >> mCoeff4 >> mCoeff5;

        // cout << "Adding gamma decay level " << initEnergy << " " << E << endl;
	if ((initEnergy - E) >= 0)  {
	  //std::cout << "here 3" << "\n";
	  //std::cout << "energy E : " << E << "\n";
	  //std::cout << "init E : " << initEnergy << "\n";
	  //std::cout << "diff : " << initEnergy - E << "\n" ;
          DecayChannel* dcGamma = new DecayChannel("Gamma", &GetDecayMode("Gamma"), E, intensity / (1. + convIntensity), lifetime, initEnergy, initEnergy - E);
          p->AddDecayChannel(dcGamma);
	} else {
	  std::cerr << "WARNING: Attempted to add gamma branch to a final state with negative excitation energy. Please check you are using the correct version of PhotonEvaporation.\nCurrent filename: " << gammaFileSS.str() << std::endl;
	}
      }
      // iss >> initEnergy >> Q >> intensity >> polarity >> lifetime >> angMom >>
      //    convIntensity >> kCoeff >> lCoeff1 >> lCoeff2 >> lCoeff3 >> mCoeff1 >>
      //    mCoeff2 >> mCoeff3 >> mCoeff4 >> mCoeff5;
      // TODO Implement conversion electrons
      // p->AddDecayChannel(new DecayChannel("ConversionElectron", Q))
    }
  }
  RegisterParticle(p);
  return true;
}

  bool DecayManager::Initialise(std::string configFilename, int argc, const char** argv) {
    ConfigOptions configOptions = ParseOptions(configFilename, argc, argv);
    return Initialise(configOptions);
  }

  bool DecayManager::Initialise(ConfigOptions _configOptions) {
    //cout << "Initialising..." << endl;
    configOptions = _configOptions;
    initStateName = configOptions.nuclearOptions.Name;
    initExcitationEn = configOptions.nuclearOptions.Energy;
    outputName = configOptions.general.Output;
    NRTHREADS = configOptions.general.Threads;
    generator.seed(configOptions.general.Seed);
    
    if (initStateName != "" && configOptions.nuclearOptions.Nucleons > 0) {
      struct stat infoRD;
      struct stat infoG;
      int i = stat(
          configOptions.envOptions.Radiationdata.c_str(),
          &infoRD);
      int j = stat(
          configOptions.envOptions.Gammadata.c_str(),
          &infoG);

      if (i == 0 && j == 0 && S_ISDIR(infoRD.st_mode) && S_ISDIR(infoG.st_mode)) {
        RegisterBasicParticles();
	RegisterBasicDecayModes();
	RegisterBasicSpectrumGenerators();
        return GenerateNucleus(initStateName, configOptions.nuclearOptions.Charge, configOptions.nuclearOptions.Nucleons);
      } else {
              std::cerr << "ERROR: Data files not found. Set Radiationdata and "
                "Gammadata to their correct folders." << std::endl;
        return false;
      }
    } else {
      std::cerr << "ERROR: Initial nucleus is not defined." << std::endl;
      return false;
    }
  }


// ---------------------------------------------------------------------------
// GenerateEvent_ROOT: runs one decay chain and returns a vector of
// ParticleData structs ready to be filled into a ROOT TTree.
// verbosity == 0 : only final-state particles (e-, e+, gamma, alpha, p, n)
// verbosity != 0 : all particles in the decay chain
// ---------------------------------------------------------------------------
std::vector<ParticleData> DecayManager::GenerateEvent_ROOT(int eventNr, int verbosity) {
  std::vector<ParticleData> vec;

  // Template for a zeroed-out ParticleData entry
  ParticleData pData_blank;
  pData_blank.event             = 0;
  pData_blank.time              = 0;
  pData_blank.name              = "";
  pData_blank.code              = 0;               
  pData_blank.excitation_energy = 0;
  pData_blank.kinetic_energy    = 0;
  pData_blank.total_energy      = 0;
  pData_blank.px                = 0;
  pData_blank.py                = 0;
  pData_blank.pz                = 0;

  double time    = 0.;
  int totEvents  = 0;

  std::vector<Particle*> particleStack;
  Particle* initState = GetNewParticle(initStateName);
  initState->SetExcitationEnergy(initExcitationEn);
  particleStack.push_back(initState);


  while (!particleStack.empty()) {
    Particle* p = particleStack.back();
    vector<Particle*> finalStates;
    double decayTime = p->GetDecayTime();

    // Determine whether this particle should be recorded
    bool isDetectable = (p->GetRawName() == "p"     ||
                         p->GetRawName() == "e+"    ||
                         p->GetRawName() == "e-"    ||
                         p->GetRawName() == "enu"    ||
                         p->GetRawName() == "enubar"    ||
                         p->GetRawName() == "gamma" ||
                         p->GetRawName() == "gammaBR" ||
                         p->GetRawName() == "alpha" ||
                         p->GetRawName() == "2He"   ||
                         p->GetRawName() == "n");

    if (verbosity == 0) {
      // Record only detectable final-state particles
      if (isDetectable) {
        vec.push_back(pData_blank);
        double mom = std::sqrt(p->GetMomentum()[1]*p->GetMomentum()[1] +
                               p->GetMomentum()[2]*p->GetMomentum()[2] +
                               p->GetMomentum()[3]*p->GetMomentum()[3]);
        vec[totEvents].event             = eventNr;
        vec[totEvents].time              = roundf(time * 10000) / 10000.;
        vec[totEvents].name              = p->GetRawName();
        vec[totEvents].code              = screening::NametoPDG(p->GetRawName()); 
        vec[totEvents].excitation_energy = p->GetExcitationEnergy();
        vec[totEvents].kinetic_energy    = p->GetKinEnergy();
        vec[totEvents].total_energy    = p->GetKinEnergy()+p->GetMass();
        vec[totEvents].px = (mom > 0.) ? p->GetMomentum()[1] / mom : 0.;
        vec[totEvents].py = (mom > 0.) ? p->GetMomentum()[2] / mom : 0.;
        vec[totEvents].pz = (mom > 0.) ? p->GetMomentum()[3] / mom : 0.;
        ++totEvents;
      }
    } else {
      // Verbose mode: record every particle in the chain
      vec.push_back(pData_blank);
      double mom = std::sqrt(p->GetMomentum()[1]*p->GetMomentum()[1] +
                             p->GetMomentum()[2]*p->GetMomentum()[2] +
                             p->GetMomentum()[3]*p->GetMomentum()[3]);
      vec[totEvents].event             = eventNr;
      vec[totEvents].time              = roundf(time * 10000) / 10000.;
      vec[totEvents].name              = p->GetRawName();
      vec[totEvents].code              = screening::NametoPDG(p->GetRawName()); 
      vec[totEvents].excitation_energy = p->GetExcitationEnergy();
      vec[totEvents].kinetic_energy    = p->GetKinEnergy();
      vec[totEvents].total_energy    = p->GetKinEnergy()+p->GetMass();
      vec[totEvents].px = (mom > 0.) ? p->GetMomentum()[1] / mom : 0.;
      vec[totEvents].py = (mom > 0.) ? p->GetMomentum()[2] / mom : 0.;
      vec[totEvents].pz = (mom > 0.) ? p->GetMomentum()[3] / mom : 0.;
      ++totEvents;
    }

    if ((time + decayTime) <= configOptions.cuts.Lifetime) {
      try {
        time += decayTime;
        finalStates = p->Decay();
      } catch (const std::invalid_argument& e) {
        std::cerr << "Decay Mode for particle " << p->GetName()
                  << " not found. Aborting." << endl;
      }
    }
    delete particleStack.back();
    particleStack.pop_back();
    if (!finalStates.empty()) {
      particleStack.insert(particleStack.end(), finalStates.begin(),
                           finalStates.end());
    }
  }
  return vec;
};


// ---------------------------------------------------------------------------
// GenerateEvent_TXT: runs one decay chain and returns formatted text lines
// for appending to a plain-text output file.
// verbosity == 0 : only detectable final-state particles are written
// verbosity != 0 : full sub-event structure (all intermediate states)
// ---------------------------------------------------------------------------
/*std::string DecayManager::GenerateEvent_TXT(int eventNr, int verbosity) {
  double time         = 0.;
  double checkTime    = 0.;
  int subEventNr      = 0;
  int totSubEvents    = 0;
  int totEvents       = 0;
  std::ostringstream eventData;
  std::ostringstream subHeader;
  std::ostringstream subEventData;
  std::vector<Particle*> particleStack;
  Particle* ini = GetNewParticle(initStateName);
  ini->SetExcitationEnergy(initExcitationEn);
  particleStack.push_back(ini);

  if (verbosity == 0) {
    // Simple mode: write only detectable particles
    while (!particleStack.empty()) {
      Particle* p = particleStack.back();
      vector<Particle*> finalStates;
      double decayTime = p->GetDecayTime();

      bool isDetectable = (p->GetRawName() == "p"     ||
                           p->GetRawName() == "e+"    ||
                           p->GetRawName() == "e-"    ||
                           p->GetRawName() == "gamma" ||
                           p->GetRawName() == "alpha" ||
                           p->GetRawName() == "2He"   ||
                           p->GetRawName() == "n");

      if (decayTime > configOptions.cuts.Lifetime && isDetectable) {
        ++totSubEvents;
        subEventData << eventNr << "\t\t"
                     << std::fixed << std::setprecision(4)
                     << roundf(time * 10000) / 10000.
                     << "\t" << p->GetInfoForFile() << "\n";
      }

      if ((time + decayTime) <= configOptions.cuts.Lifetime) {
        try {
          time += decayTime;
          finalStates = p->Decay();
        } catch (const std::invalid_argument& e) {
          cout << "Decay Mode for particle " << p->GetName()
               << " not found. Aborting." << endl;
          return "";
        }
      }
      delete particleStack.back();
      particleStack.pop_back();
      if (!finalStates.empty()) {
        particleStack.insert(particleStack.end(), finalStates.begin(),
                             finalStates.end());
      }
    }
  } else {
    // Verbose mode: write full sub-event structure
    while (!particleStack.empty()) {
      Particle* p = particleStack.back();
      vector<Particle*> finalStates;
      double decayTime = p->GetDecayTime();

      ++totSubEvents;
      subEventData << eventNr << "\t\t"
                   << std::fixed << std::setprecision(4)
                   << roundf(time * 10000) / 10000.
                   << "\t" << p->GetInfoForFile() << "\n";

      if ((time + decayTime) <= configOptions.cuts.Lifetime) {
        try {
          time += decayTime;
          finalStates = p->Decay();
          subHeader << eventNr << std::setw(8) << subEventNr
                    << "\t\t" << totSubEvents << "\n" << subEventData.str();
          subEventData.str(std::string());
          totEvents    += totSubEvents;
          totSubEvents  = 0;
          ++subEventNr;
        } catch (const std::invalid_argument& e) {
          cout << "Decay Mode for particle " << p->GetName()
               << " not found. Aborting." << endl;
          return "";
        }
      }
      delete particleStack.back();
      particleStack.pop_back();
      if (!finalStates.empty()) {
        particleStack.insert(particleStack.end(), finalStates.begin(),
                             finalStates.end());
      }
    }
  }

  totEvents += totSubEvents;
  subHeader << eventNr << std::setw(8) << subEventNr
            << "\t\t" << totSubEvents << "\n" << subEventData.str();
  eventData << eventNr << "\t\t" << totEvents << "\n" << subHeader.str();

  return eventData.str();
};*/


std::string DecayManager::GenerateEvent_TXT(int eventNr) {
  double time = 0.;
  std::ostringstream eventDataSS;
  std::vector<Particle*> particleStack;
  Particle* ini = GetNewParticle(initStateName);
  ini->SetExcitationEnergy(initExcitationEn);
  particleStack.push_back(ini);
  while (!particleStack.empty()) {
    Particle* p = particleStack.back();
    //cout << "Decaying particle " << p->GetName() << endl;
    vector<Particle*> finalStates;
    double decayTime = p->GetDecayTime();


    if ((time + decayTime) <= configOptions.cuts.Lifetime) {

      try {
        finalStates = p->Decay();
        time += decayTime;
        //cout << "Decay finished" << endl;
      } catch (const std::invalid_argument& e) {
        cout << "Decay Mode for particle " << p->GetName() << " not found. Aborting." << endl;
        return "";
      }
    } else {
      //cout << "Particle " << p->GetName() << " is stable" << endl;
      eventDataSS << eventNr << "\t" << time << "\t" << p->GetInfoForFile() << "\n";
    }
    delete particleStack.back();
    particleStack.pop_back();
    if (!finalStates.empty()) {
      particleStack.insert(particleStack.end(), finalStates.begin(),
                           finalStates.end());
    }
  }
  return eventDataSS.str();
}


bool DecayManager::MainLoop() {
  int nrParticles = configOptions.general.Loop;
  int verbosity   = configOptions.general.Verbosity_file;

  if (nrParticles < 1) {
    std::cerr << "ERROR: Incorrect number of events (" << nrParticles << ")" << std::endl;
    return true;
  }
  cout << "Starting Main Loop (" << nrParticles << " events)" << endl;

  std::ios::sync_with_stdio(false);
  boost::progress_display show_progress(nrParticles);
  boost::progress_timer t;

  // ----- ROOT output -----
  if (outputName.find("root") != std::string::npos) {
    TFile* outputFile = new TFile(outputName.c_str(), "RECREATE");
    TTree* tree = new TTree("ParticleTree", "Tree for Particle Data");

    ParticleData pData;
    tree->Branch("event",              &pData.event);
    tree->Branch("time",               &pData.time);
    tree->Branch("name",               &pData.name);
    tree->Branch("code",               &pData.code);
    tree->Branch("kinetic_energy",     &pData.kinetic_energy);
    tree->Branch("excitation_energy",  &pData.excitation_energy);
    tree->Branch("total_energy",       &pData.total_energy);
    tree->Branch("px",                 &pData.px);
    tree->Branch("py",                 &pData.py);
    tree->Branch("pz",                 &pData.pz);

    for (int i = 0; i < nrParticles; i += NRTHREADS) {
      int threads = std::min(NRTHREADS, nrParticles - i);
      std::future<std::vector<ParticleData>> f[threads];
      for (int t = 0; t < threads; t++) {
        f[t] = std::async(std::launch::async,
                          &DecayManager::GenerateEvent_ROOT, this, i + t, verbosity);
      }
      for (int t = 0; t < threads; t++) {
        for (const auto& particle : f[t].get()) {
          pData = particle;
          tree->Fill();
        }
        ++show_progress;
      }
    }
    outputFile->Write();
    outputFile->Close();
    delete outputFile;

  // ----- Plain-text output -----
  } else if (outputName.find("txt") != std::string::npos) {
    std::ofstream fileStream;
    fileStream.open(outputName.c_str());
    for (int i = 0; i < nrParticles; i += NRTHREADS) {
      int threads = std::min(NRTHREADS, nrParticles - i);
      std::future<std::string> f[threads];
      for (int t = 0; t < threads; t++) {
        f[t] = std::async(std::launch::async,
                          &DecayManager::GenerateEvent_TXT, this, i + t);
      }
      for (int t = 0; t < threads; t++) {
        fileStream << f[t].get();
        ++show_progress;
      }
    }
    fileStream.flush();
    fileStream.close();

  // ----- Unknown extension -----
  } else {
    std::cerr << "ERROR: Output file extension not recognised. "
                 "Please use .txt or .root for the output filename." << std::endl;
    return true;
  }

  std::cout << "Done! Time taken: ";
  return true;
};


}//End of CRADLE namespace
