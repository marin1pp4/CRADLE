#include "CRADLE/CustomCoefficientValidator.hh"
#include "CRADLE/Polarisation.hh"

#include <cmath>
#include <limits>
#include <stdexcept>

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

            throw std::runtime_error(
                "Non-finite correlation coefficient"
            );
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

} // namespace polarisation
} // namespace CRADLE
