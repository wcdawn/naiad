#include "diffusion.hpp"

#include "exception_handler.hpp"
#include "linalg.hpp"

namespace naiad
{

std::vector<Tridiagonal_matrix> Diffusion_solver::build_matrix() const
{

  const auto nx{geo.dx.size()};

  std::vector<Tridiagonal_matrix> A;
  A.resize(xslib.ngroup());
  for (auto & Amat : A)
  {
    Amat.sub.resize(nx-1);
    Amat.dia.resize(nx);
    Amat.sup.resize(nx-1);
  }

  // BC at x=0, i=0 (mirror)
  const auto mthis{geo.mat_map[0]};
  const auto mnext{geo.mat_map[1]};
  for (int g = 0; g < xslib.ngroup(); ++g)
  {
    const double dnext{2 
      * (xslib(mthis).diffusion[g]/geo.dx[0] * xslib(mnext).diffusion[g]/geo.dx[1])
      / (xslib(mthis).diffusion[g]/geo.dx[0] + xslib(mnext).diffusion[g]/geo.dx[1])};
    A[g].dia[0] = dnext + (xslib(mthis).sigma_t[g] - xslib(mthis).scatter[0](g,g)) * geo.dx[0];
  }

  return {};
}

Result Diffusion_solver::solve() const
{

  const std::vector<Tridiagonal_matrix> trimat{build_matrix()};

  return {};
}

} // namespace naiad
