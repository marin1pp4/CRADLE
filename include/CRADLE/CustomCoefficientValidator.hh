#ifndef CRADLE_CUSTOM_COEFFICIENT_VALIDATOR_HH
#define CRADLE_CUSTOM_COEFFICIENT_VALIDATOR_HH

#include <functional>
#include <limits>
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
    }
}

#endif
