#include "diffusion.hpp"

#include "exception_handler.hpp"
#include "linalg.hpp"

namespace naiad
{

Tridiagonal_multigroup_matrix Diffusion_solver::build_matrix() const
{

  const auto nx{geo.dx.size()};

  Tridiagonal_multigroup_matrix A{};
  A.sub.resize(xslib.ngroup());
  A.dia.resize(xslib.ngroup());
  A.sup.resize(xslib.ngroup());
  for (auto & sub: A.sub)
    sub.resize(nx-1);
  for (auto & dia : A.dia)
    dia.resize(nx);
  for (auto & sup : A.sup)
    sup.resize(nx-1);

  // BC at x=0, i=0
  const auto mthis{geo.mat_map[0]};
  const auto mnext{geo.mat_map[1]};
  for (int g = 0; g < xslib.ngroup(); ++g)
  {
    const double dnext{2 
      * (xslib(mthis).diffusion[g]/geo.dx[0] * xslib(mnext).diffusion[g]/geo.dx[1])
      / (xslib(mthis).diffusion[g]/geo.dx[0] + xslib(mnext).diffusion[g]/geo.dx[1])};
    A.dia[g][0] = dnext + (xslib(mthis).sigma_t[g] - xslib(mthis).scatter[0](g,g)) * geo.dx[0];
  }

  return {};
}

Result Diffusion_solver::solve() const
{

  const Tridiagonal_multigroup_matrix trimat{build_matrix()};

  return {};
}

} // namespace naiad
