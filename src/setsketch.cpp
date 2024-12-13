#include <setsketch.h>

namespace sketch {
namespace setsketch {
namespace detail {

/**
 * Calculates optimal parameters for set sketch compression
 * 
 * @param maxreg Maximum register value
 * @param minreg Minimum register value 
 * @param q Number of quantiles/registers
 * @return Pair containing optimal b and a parameters:
 *         b = base for exponential spacing between registers
 *         a = maxreg/b = coefficient for register values
 */
std::pair<long double, long double> optimal_parameters(const long double maxreg, const long double minreg, const long double q) noexcept {
    // Calculate base b as q-th root of ratio between max and min register values
    const long double b = std::exp(std::log(maxreg / minreg) / q);
    // Return optimal b and a=maxreg/b parameters
    return {b, maxreg / b};
}
}}}
