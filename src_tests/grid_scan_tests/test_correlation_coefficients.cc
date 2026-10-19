#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "CRADLE/CustomCoefficients.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/Polarisation.hh"

int main() {

    namespace cc = CRADLE::customcoefficients;

    // ============================================
    // 1. Physical parameters: 11C -> 11B
    // ============================================

    const double Ee = 1000.0;  // keV

    cc::CustomCoefficientContext ctx;

    ctx.Ee = Ee;
    ctx.rho = -0.75442;

    ctx.mf = 1.0;
    ctx.mgt = CRADLE::utilities::CalculateMatrixElementRatio(
        ctx.rho
    );

    ctx.rM = ctx.mgt / ctx.mf;

    ctx.ji = 1.5;
    ctx.jf = 1.5;

    ctx.Z = 5.0;
    ctx.betaType = -1.0;

    ctx.CV  = {1.0, 0.0};
    ctx.CVP = {1.0, 0.0};

    ctx.CA  = {1.2754, 0.0};
    ctx.CAP = {1.2754, 0.0};

    ctx.CT  = {0.1, 0.0};
    ctx.CTP = {0.1, 0.0};

    ctx.CS  = {0.0, 0.0};
    ctx.CSP = {0.0, 0.0};

    ctx.xi = CRADLE::utilities::CalculateXiBetaDecay(
        ctx.CS, ctx.CSP,
        ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP,
        ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt
    );

    ctx.gamma = std::sqrt(
        1.0 - std::pow(
            CRADLE::utilities::FINESTRUCTURE * ctx.Z, 2
        )
    );

    ctx.bigLambda =
        CRADLE::polarisation::BigLambdaJiJfFactor(
            ctx.ji, ctx.jf
        );

    ctx.smallLambda =
        CRADLE::polarisation::SmallLambdaJiJfFactor(
            ctx.ji, ctx.jf
        );

    // ============================================
    // 2. Transition parameters
    // ============================================

    cc::BetaTransitionParameters parameters;

    parameters.a_conf = std::numeric_limits<double>::quiet_NaN();
    parameters.b_conf = std::numeric_limits<double>::quiet_NaN();

    parameters.polMag = 1.0;
    parameters.alignment = 0.0;

    // No custom expressions for this first test
    cc::CustomCoefficientExpressions expressions;

    // ============================================
    // 3. Evaluate using the new function
    // ============================================

    const auto coeff = cc::EvaluateCorrelationCoefficients(
        Ee,
        ctx,
        parameters,
        expressions
    );

    // ============================================================
    // 4. Reference coefficients using the original CRADLE++ code
    // ============================================================

    const int betaType = -1; // Beta-plus decay

    const double a_ref =
    CRADLE::utilities::CalculateBetaNeutrinoAsymmetry(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        parameters.a_conf, parameters.b_conf,
        Ee, static_cast<int>(ctx.Z), betaType
    );

    const double b_ref =
    CRADLE::utilities::CalculateFierz(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        parameters.a_conf, parameters.b_conf,
        static_cast<int>(ctx.Z), betaType
    );

    const double A_ref =
    CRADLE::polarisation::CalculateBetaAssymetry(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        ctx.ji, ctx.jf,
        betaType, static_cast<int>(ctx.Z), Ee
    ) * parameters.polMag;

    const double B_ref =
    CRADLE::polarisation::CalculateNeutrinoAssymetry(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        ctx.ji, ctx.jf,
        betaType, static_cast<int>(ctx.Z), Ee
    ) * parameters.polMag;

    const double D_ref =
    CRADLE::polarisation::CalculateDTripleCorrelation(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        ctx.ji, ctx.jf,
        betaType, static_cast<int>(ctx.Z), Ee
    ) * parameters.polMag;

    const double c_ref =
    CRADLE::polarisation::CalculateAlignmentCorrelation(
        ctx.CS, ctx.CSP, ctx.CT, ctx.CTP,
        ctx.CV, ctx.CVP, ctx.CA, ctx.CAP,
        ctx.mf, ctx.mgt,
        ctx.ji, ctx.jf,
        betaType, static_cast<int>(ctx.Z), Ee
    ) * (-parameters.alignment);

    const double tolerance = 1e-12;

    bool passed =
    std::abs(coeff.a - a_ref) < tolerance &&
    std::abs(coeff.b - b_ref) < tolerance &&
    std::abs(coeff.c - c_ref) < tolerance &&
    std::abs(coeff.A - A_ref) < tolerance &&
    std::abs(coeff.B - B_ref) < tolerance &&
    std::abs(coeff.D - D_ref) < tolerance;

    std::cout << "\nRegression test: "
    << (passed ? "PASS" : "FAIL")
    << '\n';


    // ============================================================
    // 6. Test custom coefficient expressions
    // ============================================================

    cc::CustomCoefficientExpressions customExpressions;

    customExpressions.a = "0.2 + 0.1*Ee/2000";
    customExpressions.b = "-0.05";
    customExpressions.c = "0.1";
    customExpressions.A = "0.4";
    customExpressions.B = "-0.2";
    customExpressions.D = "0.1 + 0.2*Ee/2000";

    cc::BetaTransitionParameters customParameters = parameters;

    customParameters.polMag = 0.8;
    customParameters.alignment = 0.5;

    const auto customCoeff = cc::EvaluateCorrelationCoefficients(
        Ee,
        ctx,
        customParameters,
        customExpressions
    );

    std::cout << "\n===== Custom Coefficient Test =====\n";

    std::cout << "a = " << customCoeff.a << '\n';
    std::cout << "b = " << customCoeff.b << '\n';
    std::cout << "c = " << customCoeff.c << '\n';
    std::cout << "A = " << customCoeff.A << '\n';
    std::cout << "B = " << customCoeff.B << '\n';
    std::cout << "D = " << customCoeff.D << '\n';

    const double customTolerance = 1e-12;

    const bool customPassed =
    std::isfinite(customCoeff.a) &&
    std::isfinite(customCoeff.b) &&
    std::isfinite(customCoeff.c) &&
    std::isfinite(customCoeff.A) &&
    std::isfinite(customCoeff.B) &&
    std::isfinite(customCoeff.D) &&
    std::abs(customCoeff.a - 0.25) < customTolerance &&
    std::abs(customCoeff.b + 0.05) < customTolerance &&
    std::abs(customCoeff.c + 0.05) < customTolerance &&
    std::abs(customCoeff.A - 0.32) < customTolerance &&
    std::abs(customCoeff.B + 0.16) < customTolerance &&
    std::abs(customCoeff.D - 0.16) < customTolerance;

    std::cout << "Custom test: "
    << (customPassed ? "PASS" : "FAIL")
    << '\n';


    return (passed && customPassed) ? 0 : 1;
}
