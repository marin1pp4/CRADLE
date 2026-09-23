#ifndef RADIATIVECORRECTIONS
#define RADIATIVECORRECTIONS

#include <math.h>
#include <thread>    
#include "gsl/gsl_sf_dilog.h"

#include "CRADLE/Utilities.hh"


namespace CRADLE {

namespace radiativecorrections {
  using namespace boost::numeric::ublas;

  const double PI = 3.14159265359;
  const double EMASSC2 = 510.9989461;//keV
  const double PMASSC2 = 938272.046;//keV
  const double FINESTRUCTURE = 0.0072973525664;
  const double e = std::sqrt( (4 * PI * FINESTRUCTURE));
  const double LAMBDA = -1.2754;
  const double Vud = 0.97435 ;
  const double GF = 1.166378 * std::pow(10, -5) ;
  const double FERMICONSTANT =  1. ; //GF ; //GF * Vud ;

//////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::mt19937& get_thread_local_generator() {
        thread_local std::mt19937 gen(
            std::random_device{}() ^
            static_cast<uint32_t>(
                std::hash<std::thread::id>{}(std::this_thread::get_id())
            )
        );
        return gen;
};


inline double Spence_function(double x) {
    return -gsl_sf_dilog(x);
};


inline double DELTA_variable(double MIMASSC2, double MFMASSC2, std::string mode) {
    double E_transition = (MIMASSC2 - MFMASSC2) ;
    double Er_max = (( std::pow( (E_transition - EMASSC2)/EMASSC2, 2 ) 
                       + 2. * (E_transition - EMASSC2)/EMASSC2 )/(2 * MFMASSC2 / EMASSC2)) ;

    double d ;
    if (mode == "BetaMinus") {
        d = (E_transition - Er_max*EMASSC2)/EMASSC2 + 1.;
    } else if (mode == "BetaPlus") {
        d = (E_transition - Er_max*EMASSC2)/EMASSC2 - 1.;
    }
    //std::cout << "delta : " << d << "\n";
    return d ;
};


inline double BETA_variable(double E2) {
    return std::sqrt(1. - 1./std::pow(E2, 2)) ;
};


inline double E10_variable(double E2, double Q) {
    double DELTA = Q / EMASSC2 + 1. ;
    std::cout << "Delta : " << DELTA << std::endl;
    std::cout << "Q : " << Q << std::endl;
    return DELTA - E2;
};


inline double E1_variable(double E2, double K, double Q) {
    double DELTA = Q / EMASSC2 + 1. ;
    std::cout << "Delta : " << DELTA << std::endl;
    std::cout << "Q : " << Q << std::endl;
    return DELTA - K - E2 ;
};


inline double XHI(double MF, double MGT) {
    double x = std::pow(MF, 2) + std::pow(LAMBDA, 2) * std::pow(MGT, 2) ;
    return x ;
};


inline double XHIa(double MF, double MGT) {
    double xa = (std::pow(MF, 2) - std::pow(LAMBDA, 2) * std::pow(MGT, 2)/3.) ;
    return xa ;
};


inline double XHIa_XHI(double MF, double MGT) { // gives the electron-neutrino angular correlation a
    double a = XHIa(MF, MGT) / XHI(MF, MGT) ;
    return a ;
};


inline double p2_k_variable(double E2, double K, double COS_GAMMA) {
    double BETA = BETA_variable(E2) ;
    return E2 * K - BETA * E2 * K * COS_GAMMA ;
};


inline double P2_variable(double E2, double K, double COS_GAMMA) {
    double p2_k = p2_k_variable(E2, K, COS_GAMMA);
    return 1. / K / K + 1. / p2_k / p2_k - (2. * E2)/(K * p2_k);
};


inline double p1_p2_variable(double E2, double E1, double N1_N2) {
    double BETA = BETA_variable(E2) ;
    return BETA * E1 * E2 * N1_N2 ;
};


inline double p1_k_variable(double E1, double K, double N1_K) {
    return E1 * K *  N1_K;
};


inline double H0(double E2,double K, double COS_GAMMA, double Q) {
    double E1 = E1_variable(E2, K, Q);
    double p2_k = p2_k_variable(E2, K, COS_GAMMA);
    double P2 = P2_variable(E2, K, COS_GAMMA);
    return E1 * (-(E2 + K) * P2 + K / p2_k) ;
};


inline double H1(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double Q) {
    double E1 = E1_variable(E2, K, Q);
    double p2_k = p2_k_variable(E2, K, COS_GAMMA);
    double P2 = P2_variable(E2, K, COS_GAMMA);
    double p1_p2 = p1_p2_variable(E2, E1, N1_N2);
    double p1_k = p1_k_variable(E1, K, N1_K) ; 
    return p1_p2 * (-P2 + 1. / p2_k) + p1_k * ((E2 + K)/(K * p2_k) - 1. / p2_k / p2_k) ;
};


inline double MBR(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double MF, double MGT, double a, 
                  double MIMASSC2, int Z, int A, double Q, int Labs, bool advanced, 
                  int betaType) {
    return (16. * std::pow(FERMICONSTANT, 2) * std::pow((MIMASSC2/EMASSC2) , 2) 
            * std::pow(e, 2) * (1. * H0(E2, K, COS_GAMMA, Q) 
            + a * H1(E2, K, COS_GAMMA, N1_K, N1_N2, Q)) 
            * utilities::GetBetaCorrections(Z, A, Q, E2, betaType, Labs, advanced)
            );
}; 


inline double M0(double E2, double COS, double MF, double MGT, double a, double MIMASSC2, int Z, 
                 int A, double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2);
    
    return (16. * std::pow(FERMICONSTANT, 2) * std::pow( MIMASSC2/EMASSC2 , 2) 
            * E10 * E2 * (1. + a * BETA * COS) 
            * utilities::GetBetaCorrections(Z, A, Q, E2, betaType, Labs, advanced)
            );
};


inline double N_variable(double E2) {
    double BETA = BETA_variable(E2) ;
    return 0.5 * log( (1. + BETA)/(1. - BETA) ) ;
};


inline double Mtilde(double E2, double MF, double MGT, double MIMASSC2, int Z, int A, 
                     double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2) ;
    double N = N_variable(E2) ;
    return (- (FINESTRUCTURE/PI) * (16. * std::pow(FERMICONSTANT, 2) 
            * ((1 - std::pow(BETA, 2))/BETA) * N * std::pow((MIMASSC2/EMASSC2), 2)
            * E10 * E2) 
            * utilities::GetBetaCorrections(Z, A, Q, E2, betaType, Labs, advanced)
            );
};


inline double OMEGA_variable(double E10, double Cs){
    return E10 * Cs ;
};


inline double zVS_variable(double E2, double E10, double Cs) {
    double BETA = BETA_variable(E2) ;
    double N = N_variable(E2) ;
    double omega = OMEGA_variable(E10, Cs);
    return ((FINESTRUCTURE/PI) * (1.5 * log(PMASSC2/EMASSC2) + 2.*((N/BETA) - 1.) * log( (2.*omega)) 
            + 2. * (N/BETA) * (1. - N) + (2./BETA) * (Spence_function(2. * BETA/(1. + BETA))) - (3./8.)));
};


inline double MVS(double E2, double COS, double Cs,double MF, double MGT, double a, double MIMASSC2, 
                  int Z, int A, double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA = BETA_variable(E2) ;
    double N = N_variable(E2) ;
    double omega = OMEGA_variable(E10, Cs) ;
    double zVS = zVS_variable(E2, E10, Cs) ;
    return (zVS * M0(E2, COS, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) 
            + Mtilde(E2, MF, MGT, MIMASSC2, Z, A, Q, Labs, advanced, betaType)) ;
};


inline double g_weight(double E2, double K, double COS_GAMMA) {
    double BETA = BETA_variable(E2);
    double N = N_variable(E2) ;
    double p2_k = p2_k_variable(E2, K, COS_GAMMA) ;
    return (BETA * E2)/(2. * N * p2_k) ; 
};


inline double WH(double E2, double K, double COS_GAMMA, double N1_K, double N1_N2, double MF, double MGT, double a, 
                 double MIMASSC2, int Z, int A, double Q, int Labs, bool advanced, int betaType) {
    double E1 = E1_variable(E2, K, Q); 
    double BETA = BETA_variable(E2);
    return (K * BETA * E1 * E2 
            * MBR(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ) ;
}; 


inline double Vg_variable(double Cs, double Q){
    double DELTA = Q / EMASSC2 + 1. ;
    return -32. * std::pow(PI, 3) * (DELTA - 1.) * log(Cs);
};


inline double rho_H(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, int Z, 
                    int A, double Q, int Labs, bool advanced, int betaType) {
    double somme_rhoH = 0. ;
    double Vg = Vg_variable(Cs, Q) ;
    int nout = 0 ;

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    
    double DELTA = Q / EMASSC2 + 1. ;
    auto& gen = get_thread_local_generator();
    //for (int i = 0 ; i < n ; i++) {
    int i = 0 ;
    while (i < n) {
        double U[8];
        for (int j = 0; j < 8; j++) U[j] = distribution(gen);
        
        double E2 = 1. + (DELTA - 1.) * U[0] ;
        double E10 = E10_variable(E2, Q); 
        double omega = OMEGA_variable(E10, Cs) ;
        double K = omega * exp(-U[1] * log(Cs)) ;

        double BETA = BETA_variable(E2) ;
        double N = N_variable(E2) ; 

        double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2]))/BETA ;
        double COS_NEUTRINO = 2. * U[3] - 1. ;
        double COS_ELECTRON = 2. * U[4] - 1. ;

        double PHI_GAMMA = 2. * PI * U[5] ;
        double PHI_NEUTRINO = 2. * PI * U[6] ;
        double PHI_ELECTRON = 2. * PI * U[7] ;

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

        double wh = (WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType)
                    /(g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2/EMASSC2), 2))) ;

        if (std::isnan(wh) || std::isinf(wh)) {
            //nout += 1 ;
        } else {
            i += 1 ;
            somme_rhoH += wh ;  
        } 
    }
    //std::cout << "nout : " << nout << "\n";
    double RHOH = (Vg * somme_rhoH)/(n-nout) ;
    return RHOH ;
};


inline double WH_max(int n, double Cs, double MF, double MGT, double a, double MIMASSC2, int Z, 
                     int A, double Q, int Labs, bool advanced, int betaType) {
    double max = 0. ;

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    
    double DELTA = Q / EMASSC2 + 1. ;
    auto& gen = get_thread_local_generator();
    for (int i = 0 ; i < n ; i++) {
        double U[8];
        for (int j = 0; j < 8; j++) U[j] = distribution(gen);
        
        double E2 = 1. + (DELTA - 1.) * U[0] ;
        double E10 = E10_variable(E2, Q); 
        double omega = OMEGA_variable(E10, Cs) ;
        double K = omega * exp(-U[1] * log(Cs)) ;

        double BETA = BETA_variable(E2) ;
        double N = N_variable(E2) ; 

        double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2]))/BETA ;
        double COS_NEUTRINO = 2. * U[3] - 1. ;
        double COS_ELECTRON = 2. * U[4] - 1. ;

        double PHI_GAMMA = 2. * PI * U[5] ;
        double PHI_NEUTRINO = 2. * PI * U[6] ;
        double PHI_ELECTRON = 2. * PI * U[7] ;

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

        double wh = (WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType)
                    /(g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2/EMASSC2), 2))) ;
        if (wh > max ) {
            if (std::isnan(wh) || std::isinf(wh)) {

            } else {
                max = wh ;
            }
        }
    }
    return max ;
}; 


// Equation 5.18
inline double w0(double E2, double MF, double MGT, int Z, int A, 
                 double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2) ;
    double weight0 = (std::pow(FERMICONSTANT, 2) * BETA * E10 * E10 * E2 * E2)/(2. * std::pow(PI, 3)) ;

    return (weight0 
            * utilities::GetBetaCorrections(Z, A, Q, E2, betaType, Labs, advanced)
            );
};


// Calcul du max de w0 pour Neumann Rejection
inline double w0_max(double MF, double MGT, int Z, int A, double Q, 
                     int Labs, bool advanced, int betaType) {
    double w0max = 0 ;
    double n = 1000. ;
    
    double DELTA = Q / EMASSC2 + 1. ;
    double intervalle_E2 = (Q - 1.) / (n) ;
    std::vector<double> tableau_E2(n); 
    for (int i=0 ; i<=n ; i++) {
        tableau_E2[i] = 1. + i * intervalle_E2 ;
        
        double val = w0(tableau_E2[i], MF, MGT, Z, A, Q, Labs, advanced, betaType);
        if (w0max < val) {
            w0max = val ;
        }
    }
    std::cout << "w0max : " << w0max << "\n";
    return w0max ;
};


inline double rho0(int n, double MF, double MGT, int Z, int A, double Q, 
                   int Labs, bool advanced, int betaType) {
    double DELTA = Q / EMASSC2 + 1. ;
    double h = (DELTA - 1.) / (n) ; 
    double result = 0 ;
    for (int i=0 ; i<=n ; i++) {
        if (!std::isnan(w0(1+h*i, MF, MGT, Z, A, Q, Labs, advanced, betaType))) {
            
            double dx1 = 1. + i*h ;
            double dx2 = dx1 + h ;
            result += (h)*w0(0.5*(dx1+dx2), MF, MGT, Z, A, Q, Labs, advanced, betaType);
        }
    }
    return result ;
};


// Equation 5.20
inline double wVS(double E2, double Cs, double MF, double MGT, int Z, int A, 
                  double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2) ;
    double N = N_variable(E2) ; 
    double omega = OMEGA_variable(E10, Cs) ;
    double zVS = zVS_variable(E2, E10, Cs);
    double weightVS = (w0(E2, MF, MGT, Z, A, Q, Labs, advanced, betaType) 
                       * (zVS - (FINESTRUCTURE * N/PI) * (1. - std::pow(BETA, 2))/BETA)) ;
    
    return weightVS ;
};


// Calcul du max de wVS pour Neumann Rejection
inline double wVS_max(double Cs, double MF, double MGT, int Z, int A, 
                      double Q, int Labs, bool advanced, int betaType) {
    double wVSmax = 0 ;
    double n = 1000. ;
    double DELTA = Q / EMASSC2 + 1. ;
    double intervalle_E2 = (DELTA - 1.) / (n) ;
    std::vector<double> tableau_E2(n); 
    for (int i=0 ; i<=n ; i++) {
        tableau_E2[i] = 1. + i * intervalle_E2 ;
        
        double val = wVS(tableau_E2[i], Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType);
        if (wVSmax < val) {
            wVSmax = val ;
        }
    }
    return wVSmax ; 
};


inline double rhoVS(int n, double Cs, double MF, double MGT, int Z, int A, 
                    double Q, int Labs, bool advanced, int betaType) {
    double DELTA = Q / EMASSC2 + 1. ;
    double h = (DELTA - 1) / (n) ; 
    double result = 0 ;
    for (int i=0 ; i<=n ; i++) {
        if(!std::isnan(wVS(1 + h*i, Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType))) {

            double dx1 = 1. + i*h ;
            double dx2 = dx1 + h ;
            result += (h)*wVS(0.5*(dx1+dx2), Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType);

        }
    }
    return result ;
};


// Equation 5.18 + 5.20
inline double w0VS(double E2, double Cs, double MF, double MGT, int Z, int A, 
                   double Q, int Labs, bool advanced, int betaType) {
    return (w0(E2, MF, MGT, Z, A, Q, Labs, advanced, betaType) 
            + wVS(E2, Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType));
};


// Calcul du max de w0VS pour Neumann Rejection (spectre)
inline double w0VS_max(double Cs, double MF, double MGT, int Z, int A, 
                       double Q, int Labs, bool advanced, int betaType) {
    double w0VSmax = 0 ;
    double n = 1000. ;
    double DELTA = Q / EMASSC2 + 1. ;
    double intervalle_E2 = (DELTA - 1.) / (n) ;
    std::vector<double> tableau_E2(n); 
    for (int i=0 ; i<=n ; i++) {
        tableau_E2[i] = 1. + i * intervalle_E2 ;
        
        double val = w0VS(tableau_E2[i], Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType);
        if (w0VSmax < val) {
            w0VSmax = val ;
        }
    }
    return w0VSmax ;
};


inline double W0(double E2, double COS, double MF, double MGT, double a, double MIMASSC2, int Z, 
                 int A, double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2) ;
    return BETA * E10 * E2 * M0(E2, COS, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType);
};


inline double W0_max(double MF, double MGT, double a, double MIMASSC2, int Z, int A, 
                     double Q, int Labs, bool advanced, int betaType) {
    double W0max = 0 ;
    double n = 1000. ;
    double DELTA = Q / EMASSC2 + 1. ;
    double intervalle_E2 = (DELTA - 1.) / (n) ;
    double intervalle_cos = 2. / (n) ;
    std::vector<double> tableau_E2(n); 
    std::vector<double> tableau_cos(n); 
    for (int i=0 ; i<=n ; i++) {
        tableau_E2[i] = 1. + i * intervalle_E2 ;
        tableau_cos[i] = 1. - i * intervalle_cos ;
        std::cout << tableau_E2[i] << "\n";
    }
    for (int i=0 ; i<=n ; i++) {
        for (int j=0 ; j<=n ; j++) {
            double val = W0(tableau_E2[i], tableau_cos[j], MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType);
            if (W0max < val) {
                W0max = val ;
            }
        }
    }
    return W0max ;
};


inline double W0VS(double E2, double C, double Cs, double MF, double MGT, double a, double MIMASSC2, 
                   int Z, int A, double Q, int Labs, bool advanced, int betaType) {
    double E10 = E10_variable(E2, Q) ;
    double BETA = BETA_variable(E2) ;
    return BETA * E10 * E2 * (M0(E2, C, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) 
                              + MVS(E2, C, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ); 
};


inline double W0VS_max(double Cs, double MF, double MGT, double a, double MIMASSC2, int Z, int A, 
                       double Q, int Labs, bool advanced, int betaType) {
    double W0VSmax = 0 ;
    double n = 1000. ;
    double DELTA = Q / EMASSC2 + 1. ;
    double intervalle_E2 = (DELTA - 1.) / (n) ;
    double intervalle_cos = 2. / (n) ;
    std::vector<double> tableau_E2(n+1); // n -> n+1 CORRECTION
    std::vector<double> tableau_cos(n+1); // n -> n+1 CORRECTION
    for (int i=0 ; i<=n ; i++) {
        tableau_E2[i] = 1. + i * intervalle_E2 ;
        tableau_cos[i] = 1. - i * intervalle_cos ;
    }
    for (int i=0 ; i<=n ; i++) {
        for (int j=0 ; j<=n ; j++) {
            double val = W0VS(tableau_E2[i], tableau_cos[j], Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType);
            if (W0VSmax < val) {
                W0VSmax = val ;
            }
        }
    }
    return W0VSmax ;
};


// Efficiency of the “soft” Monte Carlo method (0VS) by Neumann rejection
inline double Efficiency_0VS(int ns, double Cs, double MF, double MGT, double a, double MIMASSC2,
                             int Z, int A, double Q, int Labs, bool advanced, int betaType) {
    double somme_W0VS = 0. ;
    int nout = 0 ;

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double DELTA = Q / EMASSC2 + 1. ;
    auto& gen = get_thread_local_generator();
    //for (int i = 0 ; i < ns ; i++) {
    int i = 0 ;
    while (i < ns) {
        double U_E2 = distribution(gen) ;
        double U_COS = distribution(gen) ;

        double E2 = 1. + (DELTA - 1.) * U_E2 ;
        double COS = 1. - 2. * U_COS ;

        double val = W0VS(E2, COS, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;

        if (std::isnan(val) || std::isinf(val)) {
            //nout += 1 ;
        } else {
            i += 1 ;
            somme_W0VS += val ;
        }
    }

    double mean_W0VS = somme_W0VS / (ns - nout) ;
    double W0VSmax = W0VS_max(Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;

    double E_0VS = 100. * mean_W0VS / W0VSmax ;
    return E_0VS ;
};


// Efficiency of the “hard” Monte Carlo method (H) by Neumann rejection
inline double Efficiency_H(int nH, double Cs, double MF, double MGT, double a, double MIMASSC2,
                           int Z, int A, double Q, int Labs, bool advanced, int betaType) {
    double somme_w = 0. ;
    int nout = 0 ;

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double DELTA = Q / EMASSC2 + 1. ;
    auto& gen = get_thread_local_generator();
    //for (int i = 0 ; i < nH ; i++) {
    int i = 0 ;
    while (i < nH) {
        double U[8];
        for (int j = 0; j < 8; j++) U[j] = distribution(gen);

        double E2 = 1. + (DELTA - 1.) * U[0] ;
        double E10 = E10_variable(E2, Q);
        double omega = OMEGA_variable(E10, Cs) ;
        double K = omega * exp(-U[1] * log(Cs)) ;

        double BETA = BETA_variable(E2) ;
        double N = N_variable(E2) ;

        double COS_GAMMA = (1. - (1. + BETA) * exp(-2. * N * U[2]))/BETA ;
        double COS_NEUTRINO = 2. * U[3] - 1. ;
        double COS_ELECTRON = 2. * U[4] - 1. ;

        double PHI_GAMMA = 2. * PI * U[5] ;
        double PHI_NEUTRINO = 2. * PI * U[6] ;
        double PHI_ELECTRON = 2. * PI * U[7] ;

        double SIN_GAMMA = std::sqrt((1. - std::pow(COS_GAMMA, 2))) ;
        double SIN_NEUTRINO = std::sqrt((1. - std::pow(COS_NEUTRINO, 2))) ;
        double SIN_ELECTRON = std::sqrt((1. - std::pow(COS_ELECTRON, 2))) ;

        double n_ELECTRON[3] = {SIN_ELECTRON * cos(PHI_ELECTRON), SIN_ELECTRON * sin(PHI_ELECTRON), COS_ELECTRON} ;
        double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0.} ;
        double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON} ;

        double n_PERPENDICULAIRE_GAMMA[3] ;
        double n_GAMMA[3] ;
        double n_NEUTRINO[3] = {SIN_NEUTRINO*cos(PHI_NEUTRINO), SIN_NEUTRINO * sin(PHI_NEUTRINO), COS_NEUTRINO};
        for (int j = 0; j < 3; j++) {
            n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA) ;
            n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA ;
        }

        double N1_N2 = n_NEUTRINO[0]*n_ELECTRON[0] + n_NEUTRINO[1]*n_ELECTRON[1] + n_NEUTRINO[2]*n_ELECTRON[2] ;
        double N1_K = n_NEUTRINO[0]*n_GAMMA[0] + n_NEUTRINO[1]*n_GAMMA[1] + n_NEUTRINO[2]*n_GAMMA[2] ;

        double w = (WH(E2, K, COS_GAMMA, N1_K, N1_N2, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType)
                   /(g_weight(E2, K, COS_GAMMA) * std::pow(2, 13) * std::pow(PI, 8) * std::pow((MIMASSC2/EMASSC2), 2))) ;

        if (std::isnan(w) || std::isinf(w)) {
            //nout += 1 ;
        } else {
            i += 1 ;
            somme_w += w ;
        }
    }

    double mean_w = somme_w / (nH - nout) ;
    double w_max = WH_max(nH, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;

    double E_H = 100. * mean_w / w_max ;
    return E_H ;
};


inline double PH(double Cs, double MF, double MGT, double a, double MIMASSC2, int Z, int A, 
                 double Q, int Labs, bool advanced, int betaType) {

    double RHOH = rho_H(1000000, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;
    double RHO0 = rho0(5000, MF, MGT, Z, A, Q, Labs, advanced, betaType) ;
    double RHOVS = rhoVS(5000, Cs, MF, MGT, Z, A, Q, Labs, advanced, betaType) ;
    double RHO0VS = RHO0 + RHOVS ;

    double E_0VS = Efficiency_0VS(100000, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;
    double E_H = Efficiency_H(100000, Cs, MF, MGT, a, MIMASSC2, Z, A, Q, Labs, advanced, betaType) ;
    std::cout << "Efficacite 0VS (%) : " << E_0VS << "\n";
    std::cout << "Efficacite H (%) : " << E_H << "\n";
    
    //std::cout << "rhoh : " << RHOH << "\n";
    //std::cout << "rho0 : " << RHO0 << "\n";
    //std::cout << "rhoVS : " << RHOVS << "\n";
    //std::cout << "rho0VS : " << RHO0VS << "\n";
    //std::cout << "pH : " << RHOH/(RHO0VS + RHOH) << "\n";
    //std::cout << "mf : " << MF_2 << "\n";
    //std::cout << "mgt : " << MGT_2 << "\n";  
    //std::cout << "wh max : " << WH_max(1000000, Cs, a, MIMASSC2, MFMASSC2, Z, betaType, mode) << "\n";
    //std::cout << "w0vs max : " << W0VS_max(MF_2, MGT_2, MIMASSC2, MFMASSC2) << "\n";"
    std::cout << "r_rho :" << 100*(RHOVS+RHOH)/RHO0 << "\n";

    return RHOH/(RHO0VS + RHOH) ; 
};

   
}//End of radiativecorrections namespace
}//End of CRADLE namespace
#endif
