#include "CRADLE/CustomCoefficientValidator.hh"
#include "CRADLE/Polarisation.hh"
#include "CRADLE/Utilities.hh"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <chrono>

namespace CRADLE {
namespace polarisation {

GridScanResult GridScanMinimum(
    double E_min,
    double E_max,
    int N_E,
    int N_u,
    int N_v,
    const std::function<CRADLE::customcoefficients::CorrelationCoefficients(double)>& evaluateCoefficients)
{
    // Validate grid parameters.
    if (N_E < 2 || N_u < 2 || N_v < 2 ||
        !std::isfinite(E_min) ||
        !std::isfinite(E_max) ||
        E_min <= EMASSC2 ||
        E_max <= E_min ||
        !evaluateCoefficients) {

        throw std::invalid_argument(
            "Invalid grid scan parameters"
        );
    }

    GridScanResult result;

    // Loop over electron energies.
    for (int i = 0; i < N_E; ++i) {

        const double Ee = E_min + i * (E_max - E_min) / (N_E - 1);

        // Evaluate correlation coefficients once per energy.
        const customcoefficients::CorrelationCoefficients coeff = evaluateCoefficients(Ee);

        // Check the evaluated coefficients.
        if (!std::isfinite(coeff.a) ||
            !std::isfinite(coeff.b) ||
            !std::isfinite(coeff.c) ||
            !std::isfinite(coeff.A) ||
            !std::isfinite(coeff.B) ||
            !std::isfinite(coeff.D)) {

            throw std::runtime_error("Non-finite correlation coefficient");
        }

        // Loop over electron polar angle.
        for (int j = 0; j < N_u; ++j) {

            const double u = -1.0 + 2.0 * j / (N_u - 1);

            // Loop over neutrino polar angle.
            for (int k = 0; k < N_v; ++k) {

                const double v = -1.0 + 2.0 * k / (N_v - 1);

                // Exact minimization over azimuthal angle phi.
                const double W = AnalyticalMinimumAngCorrFactor(coeff.a,coeff.b,coeff.c,coeff.A,coeff.B,coeff.D,Ee,u,v);

                // Check numerical validity.
                if (!std::isfinite(W)) {

                    throw std::runtime_error(
                        "Non-finite angular correlation factor"
                    );
                }

                // Update global minimum.
                if (W < result.minimumValue) {
                    result.minimumValue = W;
                    result.energy = Ee;
                    result.u = u;
                    result.v = v;
                }
            }
        }
    }

    // Check whether a negative region was detected.
    result.negativeFound = (result.minimumValue < 0.0);

    return result;
}

GridScanResult ValidateCorrelationCoefficients(
    double E_min,
    double E_max,
    int N_E,
    int N_u,
    int N_v,
    const customcoefficients::CustomCoefficientContext& ctx,
    const customcoefficients::BetaTransitionParameters& parameters,
    const customcoefficients::CustomCoefficientExpressions& expressions)
{
    auto evaluateCoefficients = [&](double Ee) {
        return customcoefficients::EvaluateCorrelationCoefficients(Ee, ctx, parameters, expressions);
    };

    return GridScanMinimum(E_min, E_max, N_E, N_u, N_v, evaluateCoefficients);
}

void ValidateTransitionOnce(
    const std::string& transitionKey,
    const std::vector<std::vector<double>>& spectrum,
    const customcoefficients::CustomCoefficientContext& ctx,
    const customcoefficients::BetaTransitionParameters& parameters,
    const customcoefficients::CustomCoefficientExpressions& expressions,
    int N_E,
    int N_u,
    int N_v)
{
    static std::mutex validationMutex;
    static std::set<std::string> validatedConfigurations;

    // Protect the complete check-validation-register sequence.
    std::lock_guard<std::mutex> lock(validationMutex);

    // ============================================================
    // 1. Determine the energy interval
    // ============================================================

    double Tmin = std::numeric_limits<double>::infinity();
    double Tmax = -std::numeric_limits<double>::infinity();

    for (const auto& point : spectrum) {

        if (point.empty()) {
            throw std::runtime_error("GridScan: invalid spectrum entry.");
        }

        const double T = point[0];

        if (!std::isfinite(T)) {
            throw std::runtime_error("GridScan: non-finite spectrum energy.");
        }

        if (T > 0.0) {
            Tmin = std::min(Tmin, T);
        }

        Tmax = std::max(Tmax, T);
    }

    if (!std::isfinite(Tmin) || !std::isfinite(Tmax) || Tmin > Tmax) {
        throw std::runtime_error("GridScan: invalid energy interval.");
    }

    const double E_min = Tmin + utilities::EMASSC2;
    const double E_max = Tmax + utilities::EMASSC2;

    // ============================================================
    // 2. Construct a configuration key
    // ============================================================

    std::ostringstream keyStream;
    keyStream << std::setprecision(17);

    keyStream << transitionKey << "|" << E_min << "|" << E_max << "|" << N_E << "|" << N_u
              << "|" << N_v << "|" << ctx.rho << "|" << ctx.rM << "|" << ctx.mf << "|" << ctx.mgt
              << "|" << ctx.ji << "|" << ctx.jf << "|" << ctx.Z << "|" << ctx.betaType << "|" << ctx.CS
              << "|" << ctx.CSP << "|" << ctx.CV << "|" << ctx.CVP << "|" << ctx.CT<< "|" << ctx.CTP
              << "|" << ctx.CA << "|" << ctx.CAP << "|" << parameters.a_conf << "|" << parameters.b_conf
              << "|" << parameters.polMag << "|" << parameters.alignment
              << "|" << static_cast<int>(parameters.calculationMode) << "|" << expressions.a
              << "|" << expressions.b << "|" << expressions.c << "|" << expressions.A
              << "|" << expressions.B << "|" << expressions.D;

    const std::string validationKey = keyStream.str();

    // ============================================================
    // 3. Skip previously validated configurations
    // ============================================================

    if (validatedConfigurations.count(validationKey) != 0) {
        return;
    }

    // ============================================================
    // 4. Run the positivity grid scan
    // ============================================================

    std::ostringstream report;

    report << std::fixed << std::setprecision(10);

    report
    << "\n========================================================\n"
    << "          ANGULAR DISTRIBUTION GRID SCAN\n"
    << "========================================================\n\n"

    << "Transition      : " << transitionKey << "\n\n"

    << "Energy range    : ["
    << E_min << ", " << E_max << "] keV\n\n"

    << "Grid resolution :\n"
    << "    Energy (NE) : " << N_E << '\n'
    << "    u      (NU) : " << N_u << '\n'
    << "    v      (NV) : " << N_v << "\n\n"

    << "Grid points     : "
    << static_cast<long long>(N_E) * N_u * N_v
    << "\n\n";

    // To measure the execution time
    const auto scanStart = std::chrono::steady_clock::now();

    // Execute the grid scan.
    const auto result = ValidateCorrelationCoefficients(E_min, E_max, N_E, N_u,N_v, ctx, parameters, expressions);

    const auto scanEnd = std::chrono::steady_clock::now();
    const double scanTime = std::chrono::duration<double>(scanEnd - scanStart).count();

    // ============================================================
    // 5. Print validation results
    // ============================================================

    const bool finiteMinimum = std::isfinite(result.minimumValue);
    const bool validationPassed =
    finiteMinimum &&
    !result.negativeFound &&
    result.minimumValue >= 0.0;

    report
    << "--------------------------------------------------------\n"
    << "MINIMUM ANGULAR DISTRIBUTION\n"
    << "--------------------------------------------------------\n\n"

    << "W_min           : " << result.minimumValue << '\n'
    << "Electron energy : " << result.energy << " keV\n"
    << "u               : " << result.u << '\n'
    << "v               : " << result.v << "\n\n"

    << "--------------------------------------------------------\n"

    << "Validation      : "
    << (validationPassed ? "PASS" : "FAIL")
    << '\n'

    << (validationPassed
    ? "No negative values found in the scanned grid.\n"
    : "Negative or non-finite angular distribution detected.\n")

    << "Grid scan time  : "
    << scanTime << " s\n"

    << "========================================================\n";

    std::cout << report.str();

    // Reject invalid distributions.
    if (!finiteMinimum) {
        throw std::runtime_error(
            "GridScan: non-finite angular distribution."
        );
    }

    if (!validationPassed) {
        throw std::runtime_error(
            "GridScan failed: negative angular distribution."
        );
    }

    // Register only successful validations.
    validatedConfigurations.insert(validationKey);

}

} // namespace polarisation
} // namespace CRADLE
