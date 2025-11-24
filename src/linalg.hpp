#ifndef NAIAD_LINALG_HPP
#define NAIAD_LINALG_HPP

#include <vector>

namespace naiad
{

std::vector<double> trid(const std::vector<double> & sub, const std::vector<double> & dia, const std::vector<double> & sup, const std::vector<double> & b);

} // namespace naiad

#endif
