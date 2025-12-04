#include "linalg.hpp"

#include <iostream>

namespace naiad
{

std::vector<double> trid(
    const std::vector<double> & sub,
    std::vector<double> dia,
    const std::vector<double> & sup,
    std::vector<double> b)
{
  for (std::size_t i = 1; i < dia.size(); ++i)
  {
    const double w{sub[i-1] / dia[i-1]};
    dia[i] -= w * sup[i-1];
    b[i] -= w * b[i-1];
  }

  std::vector<double> x;
  x.resize(b.size());
  x.back() = b.back() / dia.back();
  for (int64_t i = dia.size() - 3; i >= 0; --i)
    x[i] = (b[i] - sup[i] * x[i+1]) / dia[i];

  return x;
}

} // namespace naiad
