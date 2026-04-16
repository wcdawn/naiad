#include "transport.hpp"

#include "exception_handler.hpp"

namespace naiad
{

Transport_solver::Transport_solver(
  const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
  const XSLibrary & xslib_, const Tolerance & tol_, const Quadrature1d * const quad_)
  : geo{geo_}, bc_right{bc_right_}, xslib{xslib_}, tol{tol_}, quad{quad_}
{
  if (bc_left_ != Boundary_condition::mirror)
  {
    exception.fatal(
      std::string{"Left boundary condition must be mirror for transport solver."}
      + " specified=" + enum2str(bc_left_));
  }
}

Result Transport_solver::solve() const
{
  // TODO this will have to grow another loop for moments
  // with anisotropic scattering
  std::vector<std::vector<double>> flux;
  flux.resize(xslib.ngroup());
  for (auto & f : flux)
  {
    f.resize(geo.dx.size());
    for (auto & x : f)
      x = 1.0;
  }
  double keff{1.0};

  return {std::vector<std::vector<double>>{},1.0};
}

} // namespace naiad
