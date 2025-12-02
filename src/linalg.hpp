#ifndef NAIAD_LINALG_HPP
#define NAIAD_LINALG_HPP

#include <vector>

namespace naiad
{

// NOTE: this makes some careful copies!
std::vector<double> trid(
    const std::vector<double> & sub,
    std::vector<double> dia,
    const std::vector<double> & sup,
    std::vector<double> b);

} // namespace naiad

#endif
