#include "CRADLE/DecayMode.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/SpectrumGenerator.hh"
#include "CRADLE/RadiativeCorrections.hh"
#include "CRADLE/Polarisation.hh"

#include <string>
#include <sstream>
#include <complex>
#include <chrono>
#include "exprtk.hpp"

namespace CRADLE {

void DecayMode::FourBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, Particle* finalState4) //, ublas::vector<double>& dir1, ublas::vector<double>& dir2, ublas::vector<double>& dirg, double EnergyElectron, double EnergyNeutrino, double EnergyBrPhoton) 
{
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);
  ublas::vector<double> momentumg (4);

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  double massg = finalState4->GetMass();

  ublas::vector<double> p1 = finalState1->Get3Momentum();
  ublas::vector<double> p2 = finalState2->Get3Momentum();
  ublas::vector<double> pg = finalState4->Get3Momentum();
  ublas::vector<double> p3 = -(p1+p2+pg) ;

  double p1Norm = utilities::GetNorm(p1);
  double p2Norm = utilities::GetNorm(p2);
  double p3Norm = utilities::GetNorm(p3);
  double pgNorm = utilities::GetNorm(pg);

  momentum1(0) = std::sqrt(mass1*mass1 + std::pow(p1Norm, 2)) ;
  momentum1(1) = p1[0] ;
  momentum1(2) = p1[1] ;
  momentum1(3) = p1[2] ;

  momentum2(0) = std::sqrt(mass2*mass2 + std::pow(p2Norm, 2)) ;
  momentum2(1) = p2[0] ;
  momentum2(2) = p2[1] ;
  momentum2(3) = p2[2] ;

  momentumg(0) = std::sqrt(massg*massg + std::pow(pgNorm, 2)) ;
  momentumg(1) = pg[0] ;
  momentumg(2) = pg[1] ;
  momentumg(3) = pg[2] ;

  momentum3(0) = std::sqrt(mass3*mass3 + std::pow(p3Norm, 2)) ;
  momentum3(1) = p3[0] ;
  momentum3(2) = p3[1] ;
  momentum3(3) = p3[2] ;

  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));
  finalState4->SetMomentum(utilities::LorentzBoost(velocity, momentumg));
}

void DecayMode::ThreeBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, ublas::vector<double>& dir2, double Q) {
  //Perform decay in CoM frame
  ublas::vector<double> momentum1 = finalState1->GetMomentum();
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);

  ublas::vector<double> p2 (3);
  double p2Norm = 0.;

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  ublas::vector<double> p1 = finalState1->Get3Momentum();

  double a = mass2*mass2;
  double b = mass3*mass3;
  double c = utilities::GetNorm(p1);
  double d = Q + mass1 + mass2 + mass3 - momentum1(0);
  double e = inner_prod(p1, dir2)/c;

  double first = 1./2./(c*c*e*e-d*d);
  double second = a*a*d*d-2*a*b*d*d+4.*a*c*c*d*d*e*e-2.*a*c*c*d*d-2.*a*d*d*d*d+b*b*d*d+2*b*c*c*d*d-2.*b*d*d*d*d+c*c*c*c*d*d-2.*c*c*d*d*d*d+d*d*d*d*d*d;
  double third = a*c*e-b*c*e-c*c*c*e+c*d*d*e;

  p2Norm = first*(-std::sqrt(second)+third);
  p2 = p2Norm*dir2;

  ublas::vector<double> p3 = -(p1+p2);
  double p3Norm = utilities::GetNorm(p3);

  momentum2(0) = std::sqrt(a+p2Norm*p2Norm);
  momentum2(1) = p2(0);
  momentum2(2) = p2(1);
  momentum2(3) = p2(2);

  momentum3(0) = std::sqrt(b+p3Norm*p3Norm);
  momentum3(1) = p3(0);
  momentum3(2) = p3(1);
  momentum3(3) = p3(2);

  //std::cout << "\t" << inner_prod(p1, p2)/p2Norm/c << std::endl;

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));

}

void DecayMode::TwoBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, double Q) {
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);

  ublas::vector<double> dir = utilities::RandomDirection();

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  
  double M = Q + mass1 + mass2;

  double p = 1./(2.*M)*std::sqrt((M*M-std::pow(mass1-mass2, 2.))*(M*M-std::pow(mass1+mass2, 2.)));
  
  double p1 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass1*mass1 - mass2*mass2, 2) - 4*M*M*mass1*mass1 );
  double p2 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass2*mass2 - mass1*mass1, 2) - 4*M*M*mass2*mass2 );

  double energy1 = std::sqrt(mass1*mass1+p1*p1);
  double energy2 = std::sqrt(mass2*mass2+p2*p2);

  
  momentum1(0) = energy1;
  momentum1(1) = p1*dir[0];
  momentum1(2) = p1*dir[1];
  momentum1(3) = p1*dir[2];

  momentum2(0) = energy2;
  momentum2(1) = -p2*dir[0];
  momentum2(2) = -p2*dir[1];
  momentum2(3) = -p2*dir[2];

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  
}

std::vector<Particle*> BetaMinus::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  // std::cout << "In BetaMinus Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
    
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enubar");

  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;

  oss.str("");
  oss.clear();
  oss << "BetaMinus:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double mf = 0.;
  double mgt = 0.;
  
  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      double rho = std::stod(Type.substr(5));
      mgt = utilities::CalculateMatrixElementRatio(rho);
      mf = 1.;
    }
  } 

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout << "fierz " << fierz << " a " << a << std::endl;
  
  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));
  
  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, Q, Labs); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
  }

  //std::cout << "Found distribution" << std::endl;

  std::vector<double> p;

  double elEnergy = utilities::RandomFromDistribution(*dist)+utilities::EMASSC2;
  double elMomentum = std::sqrt(elEnergy*elEnergy-std::pow(utilities::EMASSC2, 2.));
  ublas::vector<double> enuDir = utilities::RandomDirection();

  p.push_back(1.);
  p.push_back(a*elMomentum/elEnergy);
  ublas::vector<double> eDir = utilities::GetParticleDirection(enuDir, p);

  elFourMomentum(0) = elEnergy;
  elFourMomentum(1) = elMomentum*eDir[0];
  elFourMomentum(2) = elMomentum*eDir[1];
  elFourMomentum(3) = elMomentum*eDir[2];

  e->SetMomentum(elFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, e, enu, recoil, enuDir, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(e);
  finalStates.push_back(enu);

  return finalStates;
}

std::vector<Particle*> BetaPlus::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;

  //std::cout << "In BetaPlus Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enu");

  oss.str("");
  oss.clear();

  oss << "BetaPlus:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> enubarDir = utilities::RandomDirection();

  ublas::vector<double> posFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;
  double rho = std::nan("");

  /////////ajout de SL 12/05/2023//////////////

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      rho = std::stod(Type.substr(5));
      mgt = utilities::CalculateMatrixElementRatio(rho);
      mf = 1.;
    }
  }

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  // ==========================================
  // Custom a expression with ExprTk
  // ==========================================

  double a_exprtk = std::nan("");
  std::string aCustomExpression = dm.configOptions.betaDecay.aCustom;

  if(!aCustomExpression.empty()){

    double rM = mgt / mf;

    exprtk::symbol_table<double> symbolTable;

    symbolTable.add_variable("rho", rho);
    symbolTable.add_variable("rM", rM);
    symbolTable.add_variable("CS",  CS);
    symbolTable.add_variable("CSP", CSP);
    symbolTable.add_variable("CV",  CV);
    symbolTable.add_variable("CVP", CVP);
    symbolTable.add_variable("CT",  CT);
    symbolTable.add_variable("CTP", CTP);
    symbolTable.add_variable("CA",  CA);
    symbolTable.add_variable("CAP", CAP);

    exprtk::expression<double> expression;
    expression.register_symbol_table(symbolTable);

    exprtk::parser<double> parser;
    if (parser.compile(aCustomExpression, expression)) {
      a_exprtk = expression.value();
    }
    else {
      std::cerr << "ERROR: Invalid aCustom expression: "
                << aCustomExpression << std::endl;

      std::cerr << "ExprTk parser error: "
                << parser.error() << std::endl;

      throw std::runtime_error(
        "Failed to parse aCustom expression.");
    }

  }

  // ==========================================
  // Standard CRADLE calculation
  // ==========================================

  double a = utilities::CalculateBetaNeutrinoAsymmetry(
      CS, CSP, CT, CTP, CV, CVP, CA, CAP,
      mf, mgt, a_conf, b_conf);

  // ==========================================
  // Override a if aCustom was provided
  // ==========================================

  if (!std::isnan(a_exprtk)) {
    a = a_exprtk;
  }

  double fierz = utilities::CalculateFierz(
      CS, CSP, CT, CTP, CV, CVP, CA, CAP,
      mf, mgt, a_conf, b_conf);
  
  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, E0, Labs); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    int i=0;
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }


  // ublas::vector<std::vector<double> >* dist;
  double posEnergy = utilities::RandomFromDistribution(*dist) + utilities::EMASSC2;
  double posMomentum = std::sqrt(posEnergy*posEnergy-std::pow(utilities::EMASSC2, 2.));

  std::vector<double> p;
  p.push_back(1.);
  p.push_back(a*posMomentum/posEnergy);
  ublas::vector<double> posDir = utilities::GetParticleDirection(enubarDir, p);
  posFourMomentum(0) = posEnergy;
  posFourMomentum(1) = posMomentum*posDir[0];
  posFourMomentum(2) = posMomentum*posDir[1];
  posFourMomentum(3) = posMomentum*posDir[2];

  pos->SetMomentum(posFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, pos, enubar, recoil, enubarDir, E0);


  finalStates.push_back(recoil);
  finalStates.push_back(pos);
  finalStates.push_back(enubar);

  return finalStates;
}

////////////////////////////////////////////////////////////
std::vector<Particle*> BetaMinusRadiative::Decay(Particle* initState, double Q, double daughterExEn) {
  
  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];

  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"", initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  
  std::string parentName = oss.str();
  
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-"); //+parentName);
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enubar"); //+parentName);
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gammaBR"); //+parentName);

  oss.str("");
  oss.clear();
  oss << "BetaMinusRadiative:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  ublas::vector<double> eFourMomentum (4);
  ublas::vector<double> enubarFourMomentum (4);
  ublas::vector<double> gammaFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  else {
    mgt = 1.0;
    mf = 1.0;
    } 

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  //std::cout << "fierz " << fierz << " a " << a << std::endl;
 
  int nH = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ;
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;
  
  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  bool advanced = false ;
  if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advanced = true;
    }

  auto& gen = radiativecorrections::get_thread_local_generator();
  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double WHmax = radiativecorrections::WH_max(1000000, Cs, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType);
    W = {WHmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_wHmax(0.0, W[0]); 

  double DELTA = Q/utilities::EMASSC2 + 1. ;
  int n_verif = nH ;
  
  while (n_verif == nH) {
    double wH_NR = distribution_wHmax(gen) ;

    double U[8];
    for (int j = 0; j < 8; j++) U[j] = distribution(gen);
        
    double E2 = 1. + (DELTA - 1.) * U[0] ;
    double E10 = radiativecorrections::E10_variable(E2, Q); 
    double omega = radiativecorrections::OMEGA_variable(E10, Cs) ;
    double K = omega * exp(-U[1] * log(Cs)) ;
    double E1 = E10 - K;

    double BETA = radiativecorrections::BETA_variable(E2) ;
    double N = radiativecorrections::N_variable(E2) ; 

    double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2]))/BETA ;
    double COS_NEUTRINO = 2. * U[3] - 1. ;
    double COS_ELECTRON = 2. * U[4] - 1. ;

    double PHI_GAMMA = 2. * utilities::PI * U[5] ;
    double PHI_NEUTRINO = 2. * utilities::PI * U[6] ;
    double PHI_ELECTRON = 2. * utilities::PI * U[7] ;

    double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2))) ;
    double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2))) ;
    double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2))) ;

    double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON} ; 
    double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.} ;
    double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON} ;

    double n_PERPENDICULAIRE_GAMMA[3] ; 
    double n_PERPENDICULAIRE_NEUTRINO[3] ;  
    double n_GAMMA[3] ;
    double n_NEUTRINO[3] = {SIN_NEUTRINO*cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};
    for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA) ;
        n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA ;
    }

    double N1_N2 = n_NEUTRINO[0]*n_ELECTRON[0] + n_NEUTRINO[1]*n_ELECTRON[1] + n_NEUTRINO[2]*n_ELECTRON[2] ;
    double N1_K = n_NEUTRINO[0]*n_GAMMA[0] + n_NEUTRINO[1]*n_GAMMA[1] + n_NEUTRINO[2]*n_GAMMA[2] ;

    double val = (radiativecorrections::WH(E2, K, COS_GAMMA, N1_K, N1_N2, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType)
                  /(radiativecorrections::g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(utilities::PI, 8) 
                  * std::pow((mass_i/utilities::EMASSC2), 2))) ;

    if (wH_NR <= val) {
      nH += 1 ;

      ublas::vector<double> velocity = -initState->GetVelocity();
      double eMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enubarMomentum = E1 * utilities::EMASSC2 ;
      double gammaMomentum = K * utilities::EMASSC2 ;
      
      eFourMomentum(0) = E2 * utilities::EMASSC2;
      eFourMomentum(1) = eMomentum*n_ELECTRON[0];
      eFourMomentum(2) = eMomentum*n_ELECTRON[1];
      eFourMomentum(3) = eMomentum*n_ELECTRON[2];

      enubarFourMomentum(0) = enubarMomentum ; 
      enubarFourMomentum(1) = enubarMomentum*n_NEUTRINO[0];
      enubarFourMomentum(2) = enubarMomentum*n_NEUTRINO[1];
      enubarFourMomentum(3) = enubarMomentum*n_NEUTRINO[2];

      gammaFourMomentum(0) = gammaMomentum ;
      gammaFourMomentum(1) = gammaMomentum*n_GAMMA[0];
      gammaFourMomentum(2) = gammaMomentum*n_GAMMA[1];
      gammaFourMomentum(3) = gammaMomentum*n_GAMMA[2];

      e->SetMomentum(eFourMomentum);
      enubar->SetMomentum(enubarFourMomentum);
      gamma->SetMomentum(gammaFourMomentum);
      FourBodyDecay(velocity, e, enubar, recoil, gamma);

      finalStates.push_back(e);
      finalStates.push_back(enubar);
      finalStates.push_back(recoil);
      finalStates.push_back(gamma) ;
    } 
  }
  return finalStates;
}

std::vector<Particle*> BetaMinusVirtualSoft::Decay(Particle* initState, double Q, double daughterExEn) {

  std::vector<Particle*> finalStates;
  DecayManager& dm = DecayManager::GetInstance();

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];

  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  
  std::string parentName = oss.str();
  //std::cout << "Parent name : " << parentName << std::endl;
  //std::cout << "e-"+parentName << std::endl;
  
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-"); //+parentName);
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enubar"); //+parentName);

    
  oss.str("");
  oss.clear();
  oss << "BetaMinusVirtualSoft:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  ublas::vector<double> eFourMomentum (4);
  ublas::vector<double> enubarFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  
  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  //std::cout << "a : " << a << "\n";
  //std::cout << "fierz " << fierz << " a " << a << std::endl;
 
  int n0VS = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ;
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  //std::cout << "mf : " << mf << " mgt : " << mgt << std::endl;
  //std::cout << "Q : " << Q << " A : " << A << " Z : " << Z << " betatype : " << betaType << " labs : " << Labs << std::endl;
  bool advanced = false ;
  if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advanced = true;
    }

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double W0VSmax = radiativecorrections::W0VS_max(Cs, mf, mgt, a, mass_i, Z, A, Q, Labs,advanced, betaType);
    W = {W0VSmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_w0VSmax(0.0, W[0]); 

  double DELTA = Q /utilities::EMASSC2 + 1. ;
  int n_verif = n0VS ;
  auto& gen = radiativecorrections::get_thread_local_generator();
  while (n_verif == n0VS) {
    double w0VS_NR = distribution_w0VSmax(gen) ;

    double U[5];
    for (int j = 0; j < 5; j++) U[j] = distribution(gen);

    double E2 = 1. + (DELTA - 1.) * U[0] ;
    double COS_NEUTRINO = 2. * U[1] - 1. ;
    
    if (w0VS_NR <= radiativecorrections::W0VS(E2, COS_NEUTRINO, Cs, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType)) {
      n0VS += 1 ;
      double E10 = radiativecorrections::E10_variable(E2, Q) ;
      double BETA = radiativecorrections::BETA_variable(E2);
      
      double COS_ELECTRON = 2. * U[2] - 1. ;
      double PHI_NEUTRINO = 2. * utilities::PI * U[3] ;
      double PHI_ELECTRON = 2. * utilities::PI * U[4] ;
      double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
      double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;

      ublas::vector<double> n_ELECTRON (3); 
      n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
      double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
      double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

      double n_PERPENDICULAIRE_NEUTRINO[3] ; 
      ublas::vector<double> n_NEUTRINO (3);
      for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_NEUTRINO[j] = n_ELECTRON_PRIME[j] * cos(PHI_NEUTRINO) + n_ELECTRON_SECOND[j] * sin(PHI_NEUTRINO) ;
        n_NEUTRINO[j] = n_ELECTRON[j] * COS_NEUTRINO + n_PERPENDICULAIRE_NEUTRINO[j] * SIN_NEUTRINO ;
      }

      ublas::vector<double> velocity = -initState->GetVelocity();
      double eMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enubarMomentum = E10*utilities::EMASSC2 ;
      
      eFourMomentum(0) = E2*utilities::EMASSC2;
      eFourMomentum(1) = eMomentum*n_ELECTRON[0];
      eFourMomentum(2) = eMomentum*n_ELECTRON[1];
      eFourMomentum(3) = eMomentum*n_ELECTRON[2];

      e->SetMomentum(eFourMomentum);
      ThreeBodyDecay(velocity, e, enubar, recoil, n_NEUTRINO, Q);

      finalStates.push_back(e);
      finalStates.push_back(enubar);
      finalStates.push_back(recoil);
    } 
  }
  return finalStates;
}

std::vector<Particle*> BetaMinusPolarised::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  // std::cout << "In BetaMinus Polarised Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
    
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enubar");

  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;

  oss.str("");
  oss.clear();
  oss << "BetaMinusPolarised:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  std::complex<double> CS = dm.configOptions.couplingConstants.CS;
  std::complex<double> CSP = dm.configOptions.couplingConstants.CSP;
  std::complex<double> CV = dm.configOptions.couplingConstants.CV;
  std::complex<double> CVP = dm.configOptions.couplingConstants.CVP;
  std::complex<double> CA = dm.configOptions.couplingConstants.CA;
  std::complex<double> CAP = dm.configOptions.couplingConstants.CAP;
  std::complex<double> CT = dm.configOptions.couplingConstants.CT;
  std::complex<double> CTP = dm.configOptions.couplingConstants.CTP;
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double mf = 0.;
  double mgt = 0.;
  
  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = std::stod(Type.substr(5))/1.2754; //data is experimental mixing ratio
      mf = 1.;
    }
  } 
  //mgt = 0. ;
  //mf = 1. ;
  
  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  //std::cout << "fierz " << fierz << " a " << a << std::endl;
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, recoil->GetCharge(), +1); //beta_minus 
   
  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, Q, Labs); //// changement de Q en E0 par SL 10/05/2023
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }

  double elEnergy = utilities::RandomFromDistribution(*dist)+utilities::EMASSC2;
  
  double elMomentum = std::sqrt(elEnergy*elEnergy-std::pow(utilities::EMASSC2, 2.));
  
  double align = dm.configOptions.betaDecay.Alignment;
  ublas::vector<double> polDir(3);
  polDir(0) = dm.configOptions.betaDecay.PolarisationX;
  polDir(1) = dm.configOptions.betaDecay.PolarisationY;
  polDir(2) = dm.configOptions.betaDecay.PolarisationZ;
  polDir = utilities::NormaliseVector(polDir);
  double polMag = dm.configOptions.betaDecay.PolarisationMag;

  std::vector<double> angular_mom;
  double j_i, j_f;
  try {
    angular_mom = DecayManager::GetInstance().GetParameterMC(oss.str());
    j_i = std::abs(angular_mom[0]);
    j_f = std::abs(angular_mom[1]);
  } catch (const std::invalid_argument& e){
    double j_i = utilities::GetJpi(initState->GetCharge()+initState->GetNeutrons(),initState->GetCharge(),initState->GetExcitationEnergy());
    double j_f = utilities::GetJpi(recoil->GetCharge()+recoil->GetNeutrons(),recoil->GetCharge(),recoil->GetExcitationEnergy());
    angular_mom = {j_i,j_f};
    DecayManager::GetInstance().RegisterParameterMC(oss.str(),angular_mom);
    j_i = std::abs(j_i);
    j_f = std::abs(j_f);
  }
  
  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, elEnergy, recoil->GetCharge(), +1);  
  double c = 0;
  double A = 0;
  double B = 0;
  double D = 0;

  if (j_i > 0){
    A = polarisation::CalculateBetaAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, +1, recoil->GetCharge(), elEnergy);
    B = polarisation::CalculateNeutrinoAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, +1, recoil->GetCharge(), elEnergy);
    D = polarisation::CalculateDTripleCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, +1, recoil->GetCharge(), elEnergy);
    if (j_i > 0.5){
      c = polarisation::CalculateAlignmentCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, +1, recoil->GetCharge(), elEnergy);
    }
  }
  c *= align*(-1);
  A *= polMag;
  B *= polMag;
  D *= polMag;
  
  double F_max = polarisation::AnalyticalMaximumAngCorrFactor(a, fierz, c, A, B, D, elEnergy);
  //std::cout << "a " << a <<  ", c " << c << ", A " << A << ", B " << B << ", D " << D << std::endl;
  //std::cout << "F_max:" << F_max << std::endl;
  //sampling
  double F = F_max;
  double F_point = 0;
  
 
  std::uniform_real_distribution<double> distribution_F(0.0, F_max);
  while (F > F_point){
    double F = distribution_F(dm.generator);
    ublas::vector<double> eDir = utilities::RandomDirection();
    ublas::vector<double> enuDir = utilities::RandomDirection();
    double F_point = polarisation::CalculateAngularCorrelationFactor(a, fierz, c, A, B, D, elEnergy, eDir, enuDir, polDir);
    //std::cout << "F " << F << ", F_point " << F_point << std::endl;
    if (F < F_point){
      elFourMomentum(0) = elEnergy;
      elFourMomentum(1) = elMomentum*eDir[0];
      elFourMomentum(2) = elMomentum*eDir[1];
      elFourMomentum(3) = elMomentum*eDir[2];

      e->SetMomentum(elFourMomentum);

      ublas::vector<double> velocity = -initState->GetVelocity();
      ThreeBodyDecay(velocity, e, enu, recoil, enuDir, Q);

      finalStates.push_back(recoil);
      finalStates.push_back(e);
      finalStates.push_back(enu);
      break;
    }
  }
  return finalStates;
}


////////////////////////////////////////////////////////////////

std::vector<Particle*> BetaPlusRadiative::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enu");
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gammaBR");

  oss.str("");
  oss.clear();

  oss << "BetaPlusRadiative:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame

  double mf = 0.;
  double mgt = 0.;

  ublas::vector<double> posFourMomentum (4);
  ublas::vector<double> enuFourMomentum (4);
  ublas::vector<double> gammaFourMomentum (4);

  /////////ajout de SL 12/05/2023//////////////   

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 

  ////////////////////////////////////////////////
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  //std::cout << "MF : " << mf << "\n";
  //std::cout << "MGT : " << mgt << "\n";
  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;

  int nH = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ; 
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  bool advanced = false ;
  if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advanced = true;
    }

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double WHmax = radiativecorrections::WH_max(1000000, Cs, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType);
    W = {WHmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_wHmax(0.0, W[0]); 
  
  auto& gen = radiativecorrections::get_thread_local_generator();
  double DELTA = Q /utilities::EMASSC2 + 1. ;
  int n_verif = nH ;
  while (n_verif == nH) {
    double wH_NR = distribution_wHmax(gen) ;

    double U[8];
    for (int j = 0; j < 8; j++) U[j] = distribution(gen);
    
    double E2 = 1. + (DELTA - 1.) * U[0] ;
    double E10 = radiativecorrections::E10_variable(E2, Q); 
    double omega = radiativecorrections::OMEGA_variable(E10, Cs) ;
    double K = omega * exp(-U[1] * log(Cs)) ;
    double E1 = E10 - K;

    double BETA = radiativecorrections::BETA_variable(E2) ;
    double N = radiativecorrections::N_variable(E2) ; 

    double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2]))/BETA ;
    double COS_NEUTRINO = 2. * U[3] - 1. ;
    double COS_ELECTRON = 2. * U[4] - 1. ;

    double PHI_GAMMA = 2. * utilities::PI * U[5] ;
    double PHI_NEUTRINO = 2. * utilities::PI * U[6] ;
    double PHI_ELECTRON = 2. * utilities::PI * U[7] ;

    double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2))) ;
    double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2))) ;
    double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2))) ;

    double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON} ; 
    double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.} ;
    double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON} ;

    double n_PERPENDICULAIRE_GAMMA[3] ; 
    double n_PERPENDICULAIRE_NEUTRINO[3] ;  
    double n_GAMMA[3] ;
    double n_NEUTRINO[3] = {SIN_NEUTRINO*cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};
    for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA) ;
        n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA ;
    }

    double N1_N2 = n_NEUTRINO[0]*n_ELECTRON[0] + n_NEUTRINO[1]*n_ELECTRON[1] + n_NEUTRINO[2]*n_ELECTRON[2] ;
    double N1_K = n_NEUTRINO[0]*n_GAMMA[0] + n_NEUTRINO[1]*n_GAMMA[1] + n_NEUTRINO[2]*n_GAMMA[2] ;

    double val = (radiativecorrections::WH(E2, K, COS_GAMMA, N1_K, N1_N2, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType)
                  /(radiativecorrections::g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(utilities::PI, 8) 
                  * std::pow((mass_i/utilities::EMASSC2), 2))) ;
    if (wH_NR <= val) {
      nH += 1 ;

      ublas::vector<double> velocity = -initState->GetVelocity();
      double posMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enuMomentum = E1*utilities::EMASSC2 ;
      double gammaMomentum = K*utilities::EMASSC2 ;
      
      posFourMomentum(0) = E2*utilities::EMASSC2;
      posFourMomentum(1) = posMomentum*n_ELECTRON[0];
      posFourMomentum(2) = posMomentum*n_ELECTRON[1];
      posFourMomentum(3) = posMomentum*n_ELECTRON[2];

      enuFourMomentum(0) = enuMomentum ; 
      enuFourMomentum(1) = enuMomentum*n_NEUTRINO[0];
      enuFourMomentum(2) = enuMomentum*n_NEUTRINO[1];
      enuFourMomentum(3) = enuMomentum*n_NEUTRINO[2];

      gammaFourMomentum(0) = gammaMomentum ;
      gammaFourMomentum(1) = gammaMomentum*n_GAMMA[0];
      gammaFourMomentum(2) = gammaMomentum*n_GAMMA[1];
      gammaFourMomentum(3) = gammaMomentum*n_GAMMA[2];

      pos->SetMomentum(posFourMomentum);
      enu->SetMomentum(enuFourMomentum);
      gamma->SetMomentum(gammaFourMomentum);
      FourBodyDecay(velocity, pos, enu, recoil, gamma);

      finalStates.push_back(pos);
      finalStates.push_back(enu);
      finalStates.push_back(recoil);
      finalStates.push_back(gamma) ;
    } 
  }
  return finalStates;
}


std::vector<Particle*> BetaPlusVirtualSoft::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enu");

  oss.str("");
  oss.clear();

  oss << "BetaPlusVirtualSoft:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame

  double mf = 0.;
  double mgt = 0.;

  ublas::vector<double> posFourMomentum (4);
  ublas::vector<double> enuFourMomentum (4);

  /////////ajout de SL 12/05/2023//////////////   

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 

  ////////////////////////////////////////////////
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;
  //std::cout << "a : " << a << "\n";
  
  int n0VS = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ; 
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));
  
  bool advanced = false ;
  if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advanced = true;
    }

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double W0VSmax = radiativecorrections::W0VS_max(Cs, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType);
    W = {W0VSmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_w0VSmax(0.0, W[0]); 
  
  auto& gen = radiativecorrections::get_thread_local_generator();
  double DELTA = Q /utilities::EMASSC2 + 1. ;
  int n_verif = n0VS ;
  while (n_verif == n0VS) {
    double w0VS_NR = distribution_w0VSmax(gen) ;

    double U[5];
    for (int j = 0; j < 5; j++) U[j] = distribution(gen);
    double E2 = 1. + (DELTA - 1 ) * U[0] ;
    double COS_NEUTRINO = 2. * U[1] - 1. ;    

    if (w0VS_NR < radiativecorrections::W0VS(E2, COS_NEUTRINO, Cs, mf, mgt, a, mass_i, Z, A, Q, Labs, advanced, betaType)) {
      n0VS += 1 ;
      double E10 = radiativecorrections::E10_variable(E2, Q);
      double BETA = radiativecorrections::BETA_variable(E2) ;
      double COS_ELECTRON = 2. * U[2] - 1. ;

      double PHI_NEUTRINO = 2. * utilities::PI * U[3] ;
      double PHI_ELECTRON = 2. * utilities::PI * U[4] ;
      double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
      double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;

      ublas::vector<double> n_ELECTRON (3); 
      n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
      double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
      double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

      double n_PERPENDICULAIRE_NEUTRINO[3] ; 
      ublas::vector<double> n_NEUTRINO (3);
      for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_NEUTRINO[j] = n_ELECTRON_PRIME[j] * cos(PHI_NEUTRINO) + n_ELECTRON_SECOND[j] * sin(PHI_NEUTRINO) ;
        n_NEUTRINO[j] = n_ELECTRON[j] * COS_NEUTRINO + n_PERPENDICULAIRE_NEUTRINO[j] * SIN_NEUTRINO ;
      }

      ublas::vector<double> velocity = -initState->GetVelocity();
      double posMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enuMomentum = E10*utilities::EMASSC2 ;

      posFourMomentum(0) = E2*utilities::EMASSC2;
      posFourMomentum(1) = posMomentum*n_ELECTRON[0];
      posFourMomentum(2) = posMomentum*n_ELECTRON[1];
      posFourMomentum(3) = posMomentum*n_ELECTRON[2];

      enuFourMomentum(0) = enuMomentum ; 
      enuFourMomentum(1) = enuMomentum*n_NEUTRINO[0];
      enuFourMomentum(2) = enuMomentum*n_NEUTRINO[1];
      enuFourMomentum(3) = enuMomentum*n_NEUTRINO[2];

      pos->SetMomentum(posFourMomentum);
      enu->SetMomentum(enuFourMomentum);
      ThreeBodyDecay(velocity, pos, enu, recoil, n_NEUTRINO, Q);

      finalStates.push_back(pos);
      finalStates.push_back(enu);
      finalStates.push_back(recoil);
    } 
  }
  return finalStates;
}


std::vector<Particle*> BetaPlusPolarised::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;
  
  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;

  //std::cout << "In BetaPlusPolarised Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enu");

  oss.str("");
  oss.clear();

  oss << "BetaPlusPolarised:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> enubarDir = utilities::RandomDirection();

  ublas::vector<double> posFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  /////////ajout de SL 12/05/2023//////////////

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = std::stod(Type.substr(5))/1.2754;
      mf = 1.;
    }
  }
  //std::cout << "mf : " << mf << "\n";
  //std::cout << "mgt : " << mgt << "\n";
  //mf = 1.;
  //mgt = 0.;
  ////////////////////////////////////////////////
  std::complex<double> CS = dm.configOptions.couplingConstants.CS;
  std::complex<double> CSP = dm.configOptions.couplingConstants.CSP;
  std::complex<double> CV = dm.configOptions.couplingConstants.CV;
  std::complex<double> CVP = dm.configOptions.couplingConstants.CVP;
  std::complex<double> CA = dm.configOptions.couplingConstants.CA;
  std::complex<double> CAP = dm.configOptions.couplingConstants.CAP;
  std::complex<double> CT = dm.configOptions.couplingConstants.CT;
  std::complex<double> CTP = dm.configOptions.couplingConstants.CTP;
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, recoil->GetCharge(), -1);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;
  //std::cout << "a : " << a << "\n";
  
  double Jpi_init = utilities::GetJpi(initState->GetNeutrons() + initState->GetCharge(), initState->GetCharge(), initState->GetExcitationEnergy());
  double Jpi_final = utilities::GetJpi(recoil->GetNeutrons() + recoil->GetCharge(), recoil->GetCharge(), recoil->GetExcitationEnergy());
  int Labs = std::abs(std::abs(Jpi_final) - std::abs(Jpi_init));

  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, E0, Labs); //// changement de Q en E0 par SL 10/05/2023
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }
  double posEnergy = utilities::RandomFromDistribution(*dist)+utilities::EMASSC2;
  double posMomentum = std::sqrt(posEnergy*posEnergy-std::pow(utilities::EMASSC2, 2.));
  
  double align = dm.configOptions.betaDecay.Alignment;
  ublas::vector<double> polDir(3);
  polDir(0) = dm.configOptions.betaDecay.PolarisationX;
  polDir(1) = dm.configOptions.betaDecay.PolarisationY;
  polDir(2) = dm.configOptions.betaDecay.PolarisationZ;
  polDir = utilities::NormaliseVector(polDir);
  double polMag = dm.configOptions.betaDecay.PolarisationMag;

  std::vector<double> angular_mom;
  double j_i, j_f;
  try {
    angular_mom = DecayManager::GetInstance().GetParameterMC(oss.str());
    j_i = std::abs(angular_mom[0]);
    j_f = std::abs(angular_mom[1]);
  } catch (const std::invalid_argument& e){
    double j_i = utilities::GetJpi(initState->GetCharge()+initState->GetNeutrons(),initState->GetCharge(),initState->GetExcitationEnergy());
    double j_f = utilities::GetJpi(recoil->GetCharge()+recoil->GetNeutrons(),recoil->GetCharge(),recoil->GetExcitationEnergy());
    angular_mom = {j_i,j_f};
    DecayManager::GetInstance().RegisterParameterMC(oss.str(),angular_mom);
    j_i = std::abs(j_i);
    j_f = std::abs(j_f);
  }

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf, posEnergy, recoil->GetCharge(), -1);
  //TESTING PURPOSES ONLY
  std::cout << "ANGULAR DEBUG: "
          << "mf = " << mf
          << ", mgt = " << mgt
          << ", posEnergy = " << posEnergy
          << ", a = " << a
          << std::endl;
  double c = 0;
  double A = 0;
  double B = 0;
  double D = 0;

  if (j_i > 0){
    A = polarisation::CalculateBetaAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, -1, recoil->GetCharge(), posEnergy);
    B = polarisation::CalculateNeutrinoAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, -1, recoil->GetCharge(), posEnergy);
    D = polarisation::CalculateDTripleCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, -1, recoil->GetCharge(), posEnergy);
    if (j_i > 0.5){
      c = polarisation::CalculateAlignmentCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, j_i, j_f, -1, recoil->GetCharge(), posEnergy);
    }
  }
  
  c *= align*(-1);
  A *= polMag;
  B *= polMag;
  D *= polMag;
  double F_max = polarisation::AnalyticalMaximumAngCorrFactor(a, fierz, c, A, B, D, posEnergy);

  //sampling
  double F = F_max;
  double F_point = 0;


  std::uniform_real_distribution<double> distribution_F(0.0, F_max);
  while (F > F_point){
    double F = distribution_F(dm.generator);
    ublas::vector<double> posDir = utilities::RandomDirection();
    ublas::vector<double> enubarDir = utilities::RandomDirection();
    double F_point = polarisation::CalculateAngularCorrelationFactor(a, fierz, c, A, B, D, posEnergy, posDir, enubarDir, polDir);
    if (F < F_point){
      posFourMomentum(0) = posEnergy;
      posFourMomentum(1) = posMomentum*posDir[0];
      posFourMomentum(2) = posMomentum*posDir[1];
      posFourMomentum(3) = posMomentum*posDir[2];

      pos->SetMomentum(posFourMomentum);

      ublas::vector<double> velocity = -initState->GetVelocity();
      ThreeBodyDecay(velocity, pos, enubar, recoil, enubarDir, Q);

      finalStates.push_back(recoil);
      finalStates.push_back(pos);
      finalStates.push_back(enubar);
      break;
    }
  }

  return finalStates;
}

/////////////////////////////////////////////////////

std::vector<Particle*> ShellEC::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates ;
  ublas::vector<double> velocity = -initState->GetVelocity();
  
  std::ostringstream oss ;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"EC", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enuEC");  
  recoil->SetExcitationEnergy(daughterExEn);
  
  //double mass_f = recoil -> GetMass() ;
  //std::cout << "mass f ec : " << mass_f << "\n";
  //double Excitation = decayChannel::GetParentExcitationEnergy() ;
  //std::cout << "excitation : " << excitation << "\n";
  
  TwoBodyDecay(velocity, enu, recoil, Q);
  
  finalStates.push_back(recoil) ;
  finalStates.push_back(enu) ;
  
  return finalStates;
}

/////////////////////////////////////////////////////

std::vector<Particle*> ConversionElectron::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  ublas::vector<double> velocity = -initState->GetVelocity();
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(initState->GetRawName());
  Particle* e = DecayManager::GetInstance().GetNewParticle("e+");
  recoil->SetExcitationEnergy(daughterExEn);

  TwoBodyDecay(velocity, recoil, e, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(e);

  return finalStates;
}

std::vector<Particle*> Proton::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() - 1 << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons()-1);
  Particle* p = DecayManager::GetInstance().GetNewParticle("p");
  recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, recoil, p, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(p);

  return finalStates;
}

std::vector<Particle*> Alpha::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() - 4 << utilities::atoms[initState->GetCharge()-3];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-2, initState->GetCharge()+initState->GetNeutrons()-4);
  Particle* alpha = DecayManager::GetInstance().GetNewParticle("alpha");
  recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, recoil, alpha, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(alpha);

  return finalStates;
}

std::vector<Particle*> Gamma::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  ublas::vector<double> velocity = -initState->GetVelocity();
  //std::cout << velocity[0] << "\n" ;
  //std::cout << velocity[1] << "\n" ;
  //std::cout << velocity[2] << "\n" ;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(initState->GetRawName());
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gamma");
  recoil->SetExcitationEnergy(daughterExEn);
  
  /*ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);

  ublas::vector<double> dir = utilities::RandomDirection();

  double mass1 = recoil->GetMass();
  double mass2 = gamma->GetMass();
  
  double M = Q + mass1 + mass2;

  double p = 1./(2.*M)*std::sqrt((M*M-std::pow(mass1-mass2, 2.))*(M*M-std::pow(mass1+mass2, 2.)));
  
  double p1 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass1*mass1 - mass2*mass2, 2) - 4*M*M*mass1*mass1 );
  double p2 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass2*mass2 - mass1*mass1, 2) - 4*M*M*mass2*mass2 );

  double energy1 = std::sqrt(mass1*mass1+p1*p1);
  double energy2 = std::sqrt(mass2*mass2+p2*p2);
  
  momentum1(0) = energy1;
  momentum1(1) = p1*dir[0];
  momentum1(2) = p1*dir[1];
  momentum1(3) = p1*dir[2];
  
  momentum2(0) = energy2;
  momentum2(1) = -p2*dir[0];
  momentum2(2) = -p2*dir[1];
  momentum2(3) = -p2*dir[2];
  
  recoil->SetMomentum(momentum1) ;
  gamma->SetMomentum(utilities::LorentzBoost(velocity, momentum2)) ;*/
  //double mass_i = initState -> GetMass() ;
  //double mass_f = recoil -> GetMass() ;
  //std::cout << "mass i : " << mass_i << "\n";
  //std::cout << "mass f : " << mass_f << "\n" ;
  //std::cout << "diff mass : " << mass_i - mass_f - Q << "\n";
  //Q = daughterExEn ;
  //std::cout << "Q : " << Q << "\n";
  //std::cout << "daughter exc : " << daughterExEn << "\n";

  TwoBodyDecay(velocity, recoil, gamma, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(gamma);

  return finalStates;
}

DecayMode::DecayMode() { }

DecayMode::~DecayMode() { }

void DecayMode::SetSpectrumGenerator(SpectrumGenerator* sg) {
  spectrumGen = sg;
}

BetaMinus::BetaMinus() { }

BetaPlus::BetaPlus() { }

BetaMinusRadiative::BetaMinusRadiative() { }

BetaPlusRadiative::BetaPlusRadiative() { }

BetaMinusVirtualSoft::BetaMinusVirtualSoft() { }

BetaPlusVirtualSoft::BetaPlusVirtualSoft() { }

BetaMinusPolarised::BetaMinusPolarised() { }

BetaPlusPolarised::BetaPlusPolarised() { }

ShellEC::ShellEC () { } 

ConversionElectron::ConversionElectron() { }

Proton::Proton () { }

Alpha::Alpha () { }

Gamma::Gamma () { }

}//End of CRADLE namespace
