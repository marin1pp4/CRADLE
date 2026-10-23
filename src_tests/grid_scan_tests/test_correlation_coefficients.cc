#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <complex>
#include <vector>

#include "CRADLE/CustomCoefficients.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/Polarisation.hh"
#include "CRADLE/CustomCoefficientValidator.hh"

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

    // ============================================================
    // Test: LegacyReal coefficient calculation
    // ============================================================

    cc::BetaTransitionParameters legacyParameters = parameters;

    legacyParameters.calculationMode =
    cc::CoefficientCalculationMode::LegacyReal;

    legacyParameters.polMag = 0.0;
    legacyParameters.alignment = 0.0;

    cc::CustomCoefficientExpressions noExpressions;

    const auto legacyCoeff =
    cc::EvaluateCorrelationCoefficients(
        Ee,
        ctx,
        legacyParameters,
        noExpressions
    );

    // Reference: original real-valued functions
    const double a_legacy_ref =
    CRADLE::utilities::CalculateBetaNeutrinoAsymmetry(
        ctx.CS.real(), ctx.CSP.real(),
                                                      ctx.CT.real(), ctx.CTP.real(),
                                                      ctx.CV.real(), ctx.CVP.real(),
                                                      ctx.CA.real(), ctx.CAP.real(),
                                                      ctx.mf, ctx.mgt,
                                                      legacyParameters.a_conf,
                                                      legacyParameters.b_conf
    );

    const double b_legacy_ref =
    CRADLE::utilities::CalculateFierz(
        ctx.CS.real(), ctx.CSP.real(),
                                      ctx.CT.real(), ctx.CTP.real(),
                                      ctx.CV.real(), ctx.CVP.real(),
                                      ctx.CA.real(), ctx.CAP.real(),
                                      ctx.mf, ctx.mgt,
                                      legacyParameters.a_conf,
                                      legacyParameters.b_conf
    );

    const bool legacyPassed =
    std::isfinite(legacyCoeff.a) &&
    std::isfinite(legacyCoeff.b) &&
    std::abs(legacyCoeff.a - a_legacy_ref) < 1e-12 &&
    std::abs(legacyCoeff.b - b_legacy_ref) < 1e-12;

    std::cout << "\n===== LegacyReal Test =====\n";
    std::cout << "a (new) = " << legacyCoeff.a << '\n';
    std::cout << "a (ref) = " << a_legacy_ref << '\n';
    std::cout << "b (new) = " << legacyCoeff.b << '\n';
    std::cout << "b (ref) = " << b_legacy_ref << '\n';

    std::cout << "LegacyReal test: "
    << (legacyPassed ? "PASS" : "FAIL")
    << '\n';


    // ============================================================
    // Test: Integrated angular-correlation validation
    // ============================================================

    const double E_min = 600.0;   // keV
    const double E_max = 2000.0;  // keV

    const int N_E = 50;
    const int N_u = 100;
    const int N_v = 100;

    // Use the same custom expressions and polarization parameters
    // already tested above.
    const auto scanResult =
    CRADLE::polarisation::ValidateCorrelationCoefficients(
        E_min,
        E_max,
        N_E,
        N_u,
        N_v,
        ctx,
        customParameters,
        customExpressions
    );

    std::cout << "\n===== Integrated Grid Scan Test =====\n";

    std::cout << "Minimum W = "
    << scanResult.minimumValue << '\n';

    std::cout << "Ee = "
    << scanResult.energy << " keV\n";

    std::cout << "u = " << scanResult.u << '\n';
    std::cout << "v = " << scanResult.v << '\n';

    std::cout << "Status: "
    << (scanResult.negativeFound
    ? "NEGATIVE REGION FOUND"
    : "NO NEGATIVE POINTS FOUND")
    << '\n';

    const bool scanPassed =
    std::isfinite(scanResult.minimumValue) &&
    !scanResult.negativeFound;

    std::cout << "Integrated scan test: "
    << (scanPassed ? "PASS" : "FAIL")
    << '\n';

    // ============================================================
    // Test: Energy dependence of the custom D coefficient
    // ============================================================

    std::cout << "\n===== Energy Dependence of D =====\n";

    // Physical context for 11C -> 11B beta+ decay.
    CRADLE::customcoefficients::CustomCoefficientContext dCtx;

    dCtx.mf = 1.0;

    dCtx.mgt = CRADLE::utilities::CalculateMatrixElementRatio(-0.75442);

    dCtx.rho = -0.75442;
    dCtx.rM = dCtx.mgt / dCtx.mf;

    dCtx.ji = 1.5;
    dCtx.jf = 1.5;

    dCtx.Z = 5.0;
    dCtx.betaType = -1.0;

    // Weak couplings from validation_D.txt.
    dCtx.CS  = std::complex<double>(0.0, 0.0);
    dCtx.CSP = std::complex<double>(0.0, 0.0);

    dCtx.CV  = std::complex<double>(1.0, 0.0);
    dCtx.CVP = std::complex<double>(1.0, 0.0);

    dCtx.CT  = std::complex<double>(0.0, 0.1);
    dCtx.CTP = std::complex<double>(0.0, 0.1);

    dCtx.CA  = std::complex<double>(1.2754, 0.0);
    dCtx.CAP = std::complex<double>(1.2754, 0.0);

    // Nuclear and Coulomb factors.
    dCtx.xi = CRADLE::utilities::CalculateXiBetaDecay(
        dCtx.CS, dCtx.CSP,
        dCtx.CT, dCtx.CTP,
        dCtx.CV, dCtx.CVP,
        dCtx.CA, dCtx.CAP,
        dCtx.mf, dCtx.mgt
    );

    dCtx.gamma = std::sqrt(
        1.0 - std::pow(
            CRADLE::utilities::FINESTRUCTURE * dCtx.Z, 2.0
        )
    );

    dCtx.bigLambda =
    CRADLE::polarisation::BigLambdaJiJfFactor(
        dCtx.ji, dCtx.jf
    );

    dCtx.smallLambda =
    CRADLE::polarisation::SmallLambdaJiJfFactor(
        dCtx.ji, dCtx.jf
    );

    CRADLE::customcoefficients::CustomCoefficientExpressions dExpressions;

    const std::string jsonFile =
    "../config/validation/coefficients_D.json";

    dExpressions.a =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "a"
    );

    dExpressions.b =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "b"
    );

    dExpressions.c =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "c"
    );

    dExpressions.A =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "A"
    );

    dExpressions.B =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "B"
    );

    dExpressions.D =
    CRADLE::customcoefficients::GetCustomCoefficientExpression(
        jsonFile, "11C", "D"
    );

    std::cout << "\n===== JSON EXPRESSION CHECK =====\n";

    std::cout
    << "D expression length = "
    << dExpressions.D.size() << '\n';

    std::cout
    << "D expression empty  = "
    << std::boolalpha
    << dExpressions.D.empty() << '\n';

    std::cout
    << "CT  = " << dCtx.CT << '\n'
    << "CTP = " << dCtx.CTP << '\n'
    << "CV  = " << dCtx.CV << '\n'
    << "CA  = " << dCtx.CA << '\n';


    CRADLE::customcoefficients::BetaTransitionParameters dParameters;

    dParameters.a_conf = std::numeric_limits<double>::quiet_NaN();
    dParameters.b_conf = std::numeric_limits<double>::quiet_NaN();

    dParameters.polMag = 1.0;
    dParameters.alignment = 0.0;

    dParameters.calculationMode =
    CRADLE::customcoefficients::CoefficientCalculationMode::Complex;

    const std::vector<double> energies = {
        600.0,
        800.0,
        1000.0,
        1200.0,
        1400.0
    };

    std::cout << std::fixed << std::setprecision(8);

    std::cout
    << "Ee (keV)\t"
    << "D(Ee)\t\t"
    << "beta\t\t"
    << "W_min(Ee)\n";

    for (const double Ee : energies) {

        const auto coefficients =
        CRADLE::customcoefficients::EvaluateCorrelationCoefficients(
            Ee,
            dCtx,
            dParameters,
            dExpressions
        );

        std::cout << std::scientific << std::setprecision(15);

        const double beta = std::sqrt(
            1.0 - std::pow(
                CRADLE::utilities::EMASSC2 / Ee, 2.0
            )
        );

        const double Wmin =
        1.0 - beta * std::abs(coefficients.D);

        std::cout
        << Ee << '\t'
        << coefficients.D << '\t'
        << beta << '\t'
        << Wmin << '\n';

        std::cout
        << "rho         = " << dCtx.rho << '\n'
        << "rM          = " << dCtx.rM << '\n'
        << "mf          = " << dCtx.mf << '\n'
        << "mgt         = " << dCtx.mgt << '\n'
        << "xi          = " << dCtx.xi << '\n'
        << "ji          = " << dCtx.ji << '\n'
        << "jf          = " << dCtx.jf << '\n'
        << "bigLambda   = " << dCtx.bigLambda << '\n'
        << "smallLambda = " << dCtx.smallLambda << '\n';
    }


    // ============================================================
    // Diagnostic: direct ExprTk evaluation of D
    // ============================================================

    const double E_test = 1000.0;

    dCtx.Ee = E_test;

    dCtx.gammaRatio =
    dCtx.gamma * CRADLE::utilities::EMASSC2 / E_test;

    dCtx.coulombCorr =
    CRADLE::utilities::FINESTRUCTURE * dCtx.Z /
    std::sqrt(
        E_test * E_test /
        (CRADLE::utilities::EMASSC2 *
        CRADLE::utilities::EMASSC2) - 1.0
    );

    const double D_direct =
    CRADLE::customcoefficients::EvaluateCustomCoefficient(
        dExpressions.D,
        dCtx
    );

    const auto dCoeff =
    CRADLE::customcoefficients::EvaluateCorrelationCoefficients(
        E_test,
        dCtx,
        dParameters,
        dExpressions
    );

    std::cout << std::scientific << std::setprecision(15);

    std::cout
    << "\n===== DIRECT D EVALUATION =====\n"
    << "D (ExprTk directly) = " << D_direct << '\n'
    << "D (central evaluator)= " << dCoeff.D << '\n'
    << "===============================\n";



    return (passed && customPassed && legacyPassed && scanPassed) ? 0 : 1;
}
