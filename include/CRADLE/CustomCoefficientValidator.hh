#ifndef CRADLE_CUSTOM_COEFFICIENT_VALIDATOR_HH
#define CRADLE_CUSTOM_COEFFICIENT_VALIDATOR_HH

#include <functional>
#include <limits>
#include <string>
#include <vector>
#include "CRADLE/CustomCoefficients.hh"

namespace CRADLE{
    namespace polarisation {

        struct GridScanResult {
            double minimumValue = std::numeric_limits<double>::infinity();
            double energy = 0.0;
            double u = 0.0;
            double v = 0.0;
            bool negativeFound = false;
        };

        GridScanResult GridScanMinimum(
            double E_min,
            double E_max,
            int N_E,
            int N_u,
            int N_v,
            const std::function<
            customcoefficients::CorrelationCoefficients(double)
            >& evaluateCoefficients
        );

        // Validate the angular distribution using the same
        // coefficient evaluator employed by the event generator.
        GridScanResult ValidateCorrelationCoefficients(
            double E_min,
            double E_max,
            int N_E,
            int N_u,
            int N_v,
            const customcoefficients::CustomCoefficientContext& ctx,
            const customcoefficients::BetaTransitionParameters& parameters,
            const customcoefficients::CustomCoefficientExpressions& expressions
        );

        // Validate the angular distribution once per transition/configuration.
        // A negative distribution raises std::runtime_error.
        void ValidateTransitionOnce(
            const std::string& transitionKey,
            const std::vector<std::vector<double>>& spectrum,
            const customcoefficients::CustomCoefficientContext& ctx,
            const customcoefficients::BetaTransitionParameters& parameters,
            const customcoefficients::CustomCoefficientExpressions& expressions,
            int N_E,
            int N_u,
            int N_v
        );
    }
}

#endif
