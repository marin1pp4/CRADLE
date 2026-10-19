#ifndef POLARISATION
#define POLARISATION

#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <complex>
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_result.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_sf_dilog.h"

#include "CRADLE/Utilities.hh"
/*#include "CRADLE/Screening.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/ConfigParser.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/Particle.hh"*/


namespace CRADLE {

namespace polarisation {
  using namespace boost::numeric::ublas;

  const double PI = 3.14159265359;
  const double C = 299792458;//m/s
  const double EMASSC2 = 510.9989461;//keV
  const double PMASSC2 = 938272.046;//keV
  const double NMASSC2 = 939565.4133;//keV
  const double UMASSC2 = 931494.10242;//keV
  const double ALPHAMASSC2 = 3727379.4066;//keV
  const double FINESTRUCTUREMASSC2 = 3727379.508;//keV
  const double FINESTRUCTURE = 0.0072973525664;
  const double E = 2.718281828459045;
  const double HBAR = 6.58211889e-16;//ev*s
  const double NATURALLENGTH = HBAR*C/EMASSC2/1000.;//m
  const double a_CORR = -1.0;
  const double EULER_MASCHERONI_CONSTANT = 0.577215664901532;  /**< the Euler-Mascheroni constant */

  //Polarisation Variables (from Nuclear Physics 4 (1957) 206-212; by J. D. Jackson)

  inline double SmallLambdaJiJfFactor(double j_i, double j_f){
    if (j_f == (j_i - 1)){
      return 1;
    }
    else if (j_f == j_i){
      return 1/(j_i+1);
    }
    else if (j_f == (j_i + 1)){
      return -j_i/j_f;
    }
    else
      return 0;
  }

  inline double BigLambdaJiJfFactor(double j_i, double j_f){
    if (j_f == (j_i - 1)){
      return 1;
    }
    else if (j_f == j_i){
      return -(2*j_i-1)/(j_i+1);
    }
    else if (j_f == (j_i + 1)){
      return j_i*(2*j_i-1)/j_f/(2*j_f+1);
    }
    else
      return 0;
  }
  
  inline double CalculateBetaAssymetry(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1); //alpha*Z*m_e/p_e
    double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.))*EMASSC2/energy; //gamma*m_e/E
    double A = mgt*mgt*SmallLambdaJiJfFactor(j_i,j_f)*(2.*betaType*(ct*conj(ctp)-ca*conj(cap)).real()+2.*coulombCorr*(ct*conj(cap)+ctp*conj(ca)).imag()); 
    if (j_i == j_f){ //Mixed decay factors
      A += mf*mgt*std::sqrt(j_i/(j_i+1))*(2.*(cs*conj(ctp)+csp*conj(ct)-cv*conj(cap)-cvp*conj(ca)).real()+
                                          2*betaType*coulombCorr*(cs*conj(cap)+csp*conj(ca)-cv*conj(ctp)-cvp*conj(ct)).imag()); 
      }
    A /= utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    return A;
  }

  inline double CalculateNeutrinoAssymetry(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.))*EMASSC2/energy; //gamma*m_e/E
    double B = 2.*mgt*mgt*SmallLambdaJiJfFactor(j_i,j_f)*(gamma_ratio*(ct*conj(cap)+ctp*conj(ca)).real()+betaType*(ct*conj(ctp)+ca*conj(cap)).real());
    if (j_i == j_f){ //Mixed decay factors
      B -= mf*mgt*std::sqrt(j_i/(j_i+1))*(2.*(cs*conj(ctp)+csp*conj(ct)+cv*conj(cap)+cvp*conj(ca)).real()+
                                              2.*betaType*gamma_ratio*(cs*conj(cap)+csp*conj(ca)+cv*conj(ctp)+cvp*conj(ct)).real()); 
    }
    B /= utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    return B;
  }

  inline double CalculateDTripleCorrelation(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    //Only mixed decay factors
    if (j_i == j_f){ 
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1);
    double D = mf*mgt*std::sqrt(j_i/(j_i+1))*(-2.*coulombCorr*betaType*(cs*conj(ca)+csp*conj(cap)-cv*conj(ct)-cvp*conj(ctp)).real()+
                                              2.*(cs*conj(ct)+csp*conj(ctp)-cv*conj(ca)-cvp*conj(cap)).imag());
    return D/utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    }
    else
      return 0;
  }

  inline double CalculateAlignmentCorrelation(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1);
    double c = mgt*mgt*BigLambdaJiJfFactor(j_i,j_f)*(norm(ct)+norm(ctp)-norm(ca)-norm(cap)+2*betaType*coulombCorr*(ct*conj(ca)+ctp*conj(cap)).imag());
    return c/utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
  }
  

  //Maximum Functions for c, A, B, D

  //A: iterate over the energy due to non-zero coulomb Corrections, dependnt on energy
  inline double* MaximumA(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
    double* max_A_point = new double[2]; // (energy of the maximum, value of the maximum for the abs)
    max_A_point[1] = -1.0; //since we look at the abs, first one will trigger; 
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < E0+0.1){
      double en = kin_en + EMASSC2;
      double A = CalculateBetaAssymetry(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(A) > max_A_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous B: " << max_B << " ; new B: " << B << std::endl;
        max_A_point[0] = en;
        max_A_point[1] = std::abs(A);
      }
      kin_en += 1;
    }
    return max_A_point;
  }

  //B: iterate over the energy
  inline double* MaximumB(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
    double* max_B_point = new double[2]; // (energy of the maximum, value of the maximum for the abs)
    max_B_point[1] = -1.0; //since we look at the abs
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < E0+0.1){
      double en = kin_en + EMASSC2;
      double B = CalculateNeutrinoAssymetry(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(B) > max_B_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous B: " << max_B << " ; new B: " << B << std::endl;
        max_B_point[0] = en;
        max_B_point[1] = std::abs(B);
      }
      kin_en += 1;
    }
    return max_B_point;
  }

  //D: iterate over the energy
  inline double* MaximumD(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
     double* max_D_point = new double[2];  // (energy of the maximum, value of the maximum)
    max_D_point[1] = -1.0;
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < (E0 + 0.1)){
      double en = kin_en + EMASSC2;
      double D = CalculateDTripleCorrelation(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(D) > max_D_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous D: " << max_D << " ; new D: " << D << std::endl;
        max_D_point[0] = en;
        max_D_point[1] = std::abs(D);
      }
      kin_en += 1;
    }
    return max_D_point;
  }
  
  inline double CalculateAngularCorrelationFactor(double a, double A, double B, double D, double E, double cosTheta_e, double cosTheta_enu, double phi){
    /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. No c term, kept for compatibility purposes*/
    ublas::vector<double> elDir (3);
    ublas::vector<double> enuDir (3);
    
    //electron in XZ plane, using axial symmetry in Z direction (direction of J)
    double sinTheta_e = std::sqrt(1-cosTheta_e*cosTheta_e);
    elDir(0) = sinTheta_e;
    elDir(1) = 0; 
    elDir(2) = cosTheta_e;

    double sinTheta_enu = std::sqrt(1-cosTheta_enu*cosTheta_enu);
    enuDir(0) = sinTheta_enu*std::cos(phi);
    enuDir(1) = sinTheta_enu*std::sin(phi);
    enuDir(2) = cosTheta_enu;
    
    double beta_e = std::sqrt(1-EMASSC2*EMASSC2/E/E); // p_e/E_e; p_nu/E_nu = 1
    
    double angCorrFactor = 1.;
    angCorrFactor += a*beta_e*(elDir(0)*enuDir(0)+elDir(1)*enuDir(1)+elDir(2)*enuDir(2));
    angCorrFactor += A*beta_e*elDir(2);
    angCorrFactor += B*enuDir(2);
    angCorrFactor += D*beta_e*(elDir(0)*enuDir(1)-elDir(1)*enuDir(0));
    
    
    return angCorrFactor;
  }

  inline double CalculateAngularCorrelationFactor(double a, double b, double c, double A, double B, double D, double E, double cosTheta_e, double cosTheta_enu, double phi){
    /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. Includes b and c term, and assumes perfect alignment. Assumes J/|J| is a unit vector in the z component. Used in tests*/
    vector<double> elDir (3);
    vector<double> enuDir (3);
    
    //electron in XZ plane, using axial symmetry in Z direction (direction of J)
    double sinTheta_e = std::sqrt(1-cosTheta_e*cosTheta_e);
    elDir(0) = sinTheta_e;
    elDir(1) = 0; 
    elDir(2) = cosTheta_e;

    double sinTheta_enu = std::sqrt(1-cosTheta_enu*cosTheta_enu);
    enuDir(0) = sinTheta_enu*std::cos(phi);
    enuDir(1) = sinTheta_enu*std::sin(phi);
    enuDir(2) = cosTheta_enu;

    double beta_e = std::sqrt(1-EMASSC2*EMASSC2/E/E); // p_e/E_e; p_nu/E_nu = 1
    
    double angCorrFactor = 1.;
    angCorrFactor += b*EMASSC2/E;
    angCorrFactor += a*beta_e*(elDir(0)*enuDir(0)+elDir(1)*enuDir(1)+elDir(2)*enuDir(2));
    angCorrFactor += A*beta_e*elDir(2);
    angCorrFactor += B*enuDir(2);
    angCorrFactor += D*beta_e*(elDir(0)*enuDir(1)-elDir(1)*enuDir(0));
    angCorrFactor -= c*beta_e*(elDir(0)*enuDir(0)+elDir(1)*enuDir(1)-2*elDir(2)*enuDir(2))/3;
    

    return angCorrFactor;
  }

  inline double CalculateAngularCorrelationFactor(double a, double b, double c, double A, double B, double D, double E, vector<double> elDir, vector<double> enuDir, vector<double> polDir){
    /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above. Includes b and c term, noting that c, A, B and D already account for the polarisation and alignment dependent factors. Here starting from vectors themselves, so j isn't bound to z axis. polDir is unit vector*/

    double beta_e = std::sqrt(1-EMASSC2*EMASSC2/E/E); // p_e/E_e; p_nu/E_nu = 1
    
    double angCorrFactor = 1.;
    angCorrFactor += b*EMASSC2/E;
    angCorrFactor += a*beta_e*inner_prod(elDir,enuDir);
    angCorrFactor += A*beta_e*inner_prod(elDir,polDir);
    angCorrFactor += B*inner_prod(enuDir,polDir);
    angCorrFactor += D*beta_e*inner_prod(polDir,utilities::CrossProduct(elDir,enuDir));
    angCorrFactor += c*beta_e*(inner_prod(elDir,enuDir)/3-inner_prod(elDir,polDir)*inner_prod(enuDir,polDir));

    return angCorrFactor;
  }

  inline double MaximumAngCorrFactor(double a, double A, double B, double D, double E){
    /*Search of the maximum using grid search (Obsolete, use the Analytical method for performance, kept for compatibility in tests)*/
    double maxAngCorrFactor =  std::numeric_limits<double>::lowest(); //not sure if value is nonpositive
    for (int z_e = -200; z_e <= 200; z_e++){
      for (int z_enu = -200; z_enu <= 200; z_enu++){
        for (int phi = 0; phi < 720; phi++){
          double cosTheta_e = ((double) z_e)/200;
          double cosTheta_enu = ((double) z_enu)/200;
          double radPhi = phi*PI/360;
          double angCorrFactor = CalculateAngularCorrelationFactor(a, A, B, D, E, cosTheta_e, cosTheta_enu, radPhi);
          if (angCorrFactor > maxAngCorrFactor)
	          maxAngCorrFactor = angCorrFactor;
        }
      }
      if (z_e%10 == 0) std::cout << "cos(theta_e) = "<< ((double) z_e)/200 << std::endl;
    }
    return maxAngCorrFactor;
  }

  inline double MaximumF(double a, double A, double B, double K, double znu){
    return std::sqrt(std::pow(a*znu+A,2)+K*K*(1-znu*znu))+B*znu;
  }

  inline double AnalyticalMaximumAngCorrFactor(double a, double b, double c, double A, double B, double D, double E){                     
    /*Search of the maximum value analitically*/
    double beta = std::sqrt(1-EMASSC2*EMASSC2/E/E);
    //scale a, c, A and D by beta. Note that c is multiplied already by the alignment, and A, B and D are multiplied by J as well
    a *= beta;
    c *= beta;
    A *= beta;
    D *= beta;

    double K = std::sqrt(D*D+(a+c/3)*(a+c/3));
    double a_st = a-2.*c/3;

    double A_m = (a_st*a_st-K*K)*(a_st*a_st-K*K-B*B);
    double B_m = 2*a_st*A*(a_st*a_st-K*K-B*B);
    double C_m = a_st*a_st*A*A-A*A*B*B-B*B*K*K;
    double znu_m, znu_m2; //candidates for maximum
    vector<double> F_cand(4, 0);

    //computing the values at the extrema of the interval
    F_cand(0) = MaximumF(a_st, A, B, K, 1);
    F_cand(1) = MaximumF(a_st, A, B, K, -1);

    if (A_m == 0) {
      znu_m = -C_m/B_m;
      if ((znu_m > -1) && (znu_m < 1)) {
        F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
            }
          } else if (B_m == 0 && C_m == 0){
            F_cand(2) = MaximumF(a_st, A, B, K, 0); //cover some edge cases like only D non-zero that should never happen in reality
          } else {
            double det = B_m*B_m - 4*A_m*C_m;
            if ((det < 0) && (det > -1e-16)) det = 0; //avoid spurious cases
            if (det >= 0){
        znu_m = (-B_m+std::sqrt(det))/2/A_m;
        if ((znu_m > -1) && (znu_m < 1)) {
          F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
        } 
        znu_m2 = (-B_m-std::sqrt(det))/2/A_m;
        if ((znu_m2 > -1) && (znu_m2 < 1)) {
          F_cand(3) = MaximumF(a_st, A, B, K, znu_m2);
        }
      }
    }
    // std::cout << F_cand << std::endl;
    double F_max = *std::max_element(F_cand.begin(),F_cand.end());
    F_max += 1 + b*EMASSC2/E; //adding the constant terms
    return F_max;
  }
  
  inline vector<double> MaximumAngCorrFactorPos(double a, double b, double c, double A, double B, double D, double E){
    /*Search of the position of the maximum analitically*/
    double beta = std::sqrt(1-EMASSC2*EMASSC2/E/E);
    //scale a, c, A and D by beta. Note that c is multiplied already by the alignment, and A, B and D are multiplied by J as well
    a *= beta;
    c *= beta;
    A *= beta;
    D *= beta;
      
    vector<double> max_pos(3);

    double K = std::sqrt(D*D+(a+c/3)*(a+c/3));
    double a_st = a-2.*c/3;

    max_pos(2) = std::atan2(D,a+c/3); //phi

    //F_max = max(sqrt((a_st*u+A)**2+K**2*(1-u**2))+B*u)

    double A_m = (a_st*a_st-K*K)*(a_st*a_st-K*K-B*B);
    double B_m = 2*a_st*A*(a_st*a_st-K*K-B*B);
    double C_m = a_st*a_st*A*A-A*A*B*B-B*B*K*K;
    double znu_m, znu_m2; //candidates for maximum
    vector<double> F_cand(4, 0);
    double F_max;
    //computing the values at the extrema of the interval
    F_cand(0) = MaximumF(a_st, A, B, K, 1);
    F_cand(1) = MaximumF(a_st, A, B, K, -1);

    if (A_m == 0 && B_m != 0) {
      znu_m = -C_m/B_m;
      if ((znu_m > -1) && (znu_m < 1)) {
        F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
            }
          } else if (B_m == 0 && C_m == 0){
            F_cand(2) = MaximumF(a_st, A, B, K, 0); //cover some edge cases like only D non-zero that should never happen in reality
            znu_m = 0;
          } else{
            double det = B_m*B_m - 4*A_m*C_m;
            if ((det < 0) && (det > -1e-16)) det = 0; //avoid spurious cases
            if (det >= 0){
        znu_m = (-B_m+std::sqrt(det))/2/A_m;
        if ((znu_m > -1) && (znu_m < 1)) {
          F_cand(2) = MaximumF(a_st, A, B, K, znu_m);
        } 
        znu_m2 = (-B_m-std::sqrt(det))/2/A_m;
        if ((znu_m2 > -1) && (znu_m2 < 1)) {
          F_cand(3) = MaximumF(a_st, A, B, K, znu_m2);
        }
      }
    }

    std::cout << znu_m << ", " <<  znu_m2 << std::endl;
    
    int indexMax = std::distance(F_cand.begin(),std::max_element(F_cand.begin(),F_cand.end()));
    switch (indexMax){
    case 0:
      max_pos(1) = 1; //znu
      max_pos(0) = std::copysign(1.,A + a_st); //ze
      break;
    case 1:
      max_pos(1) = -1;
      max_pos(0) = std::copysign(1.,A - a_st);
      break;
    case 2:
      max_pos(1) = znu_m;
      max_pos(0) = (a_st*znu_m+A)/(MaximumF(a_st,A,0,K,znu_m));
      break;
    case 3:
      max_pos(1) = znu_m2;
      max_pos(0) = (a_st*znu_m+A)/(MaximumF(a_st,A,0,K,znu_m2));
      break;
    }
    return max_pos;
  }

  inline double AnalyticalMinimumAngCorrFactor(double a, double b, double c, double A, double B, double D, double E, double u,double v){

    // Electron/positron velocity
    const double beta = std::sqrt(
      1.0 - EMASSC2 * EMASSC2 / (E * E)
    );

    // Geometric factor
    const double s = std::sqrt(
      (1.0 - u * u) * (1.0 - v * v)
    );

    // Phi-independent contribution
    const double W0 =
    1.0
    + b * EMASSC2 / E
    + beta * (a - 2.0 * c / 3.0) * u * v
    + A * beta * u
    + B * v;

    // Coefficients of cos(phi) and sin(phi)
    const double C =
    beta * s * (a + c / 3.0);

    const double S =
    D * beta * s;

    // Exact minimum with respect to phi
    return W0 - std::sqrt(C * C + S * S);
  }
  
}//closing polarisation namespace
}//closing CRADLE namespace
#endif
