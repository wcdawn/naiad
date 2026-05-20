#include "result.hpp"

#include "xslibrary.hpp"

namespace naiad
{

void Result::calculate_power(const std::vector<const XSMaterial *> & mat_map)
{
  power.resize(mat_map.size());
  for (std::size_t i = 0; i < mat_map.size(); ++i)
  {
    if (mat_map[i]->isfis)
    {
      for (std::size_t g = 0; g < phi.size(); ++g)
      {
        // NOTE: using nusf reaction rate for now.
        // This is representative of power if kappa/nu is constant.
        if (!mat_map[i]->sigma_f.empty())
          power[i] += mat_map[i]->sigma_f[g] * phi[g][i * (pnorder + 1) + 0];
        else
          power[i] += mat_map[i]->nusf[g] * phi[g][i * (pnorder + 1) + 0];
      }
    }
  }
}

} // namespace naiad
