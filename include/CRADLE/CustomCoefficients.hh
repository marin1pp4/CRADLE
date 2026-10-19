#ifndef CRADLE_CUSTOM_COEFFICIENTS_HH
#define CRADLE_CUSTOM_COEFFICIENTS_HH

#include <complex>
#include <limits>
#include <string>

namespace CRADLE {
namespace customcoefficients {

// Physical context used to evaluate custom beta-decay coefficients.
struct CustomCoefficientContext {

    double Ee = 0.;
    double rho = std::numeric_limits<double>::quiet_NaN();
    double rM = std::numeric_limits<double>::quiet_NaN();

    double mf = 0.;
    double mgt = 0.;
    double ji = 0.;
    double jf = 0.;
    double Z = 0.;
    double betaType = 0.;
    double xi = 0.;
    double gamma = 1.;
    double gammaRatio = 0.;
    double coulombCorr = 0.;
    double bigLambda = 0.;
    double smallLambda = 0.;

    std::complex<double> CS{0., 0.};
    std::complex<double> CSP{0., 0.};
    std::complex<double> CV{0., 0.};
    std::complex<double> CVP{0., 0.};
    std::complex<double> CT{0., 0.};
    std::complex<double> CTP{0., 0.};
    std::complex<double> CA{0., 0.};
    std::complex<double> CAP{0., 0.};
};

// Read a coefficient expression from a JSON file.
std::string GetCustomCoefficientExpression(
    const std::string& filename,
    const std::string& transition,
    const std::string& coefficient
);

// Evaluate a coefficient expression using ExprTk.
double EvaluateCustomCoefficient(
    const std::string& expressionString,
    const CustomCoefficientContext& ctx
);

struct BetaTransitionParameters {

    // Bypass coefficient configuration
    double a_conf = std::numeric_limits<double>::quiet_NaN();
    double b_conf = std::numeric_limits<double>::quiet_NaN();

    // Nuclear polarization and alignment
    double polMag = 0.0;
    double alignment = 0.0;

};

struct CustomCoefficientExpressions {

    std::string a;
    std::string b;
    std::string c;
    std::string A;
    std::string B;
    std::string D;

};

struct CorrelationCoefficients {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double A = 0.0;
    double B = 0.0;
    double D = 0.0;
};

// Calculate standard beta-decay correlation coefficients,
// apply custom expressions when provided, and include
// polarization and alignment factors.
//
// Returns effective coefficients ready for angular sampling.
CorrelationCoefficients EvaluateCorrelationCoefficients(
    double Ee,
    const CustomCoefficientContext& ctx,
    const BetaTransitionParameters& parameters,
    const CustomCoefficientExpressions& expressions
);

} // namespace customcoefficients
} // namespace CRADLE

#endif
