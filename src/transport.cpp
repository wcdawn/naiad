#include "transport.hpp"

#include "exception_handler.hpp"

namespace naiad
{

Transport_solver::Transport_solver(
  const Geometry & geo_, const Spatial_method & spatial_method,
  const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
  const XSLibrary & xslib_, const Tolerance & tol_, const Quadrature1d * const quad_)
  : geo{geo_}, bc_left{bc_left_}, bc_right{bc_right_}, xslib{xslib_}, tol{tol_}, quad{quad_}
{
  switch (spatial_method)
  {
    case (Spatial_method::diamond_difference):
      sweeper = std::make_unique<Diamond_difference_sweeper>(geo, bc_left, bc_right, xslib, quad);
      break;
    default:
      exception.fatal(std::string{"Unimplemented spatial transport method. requested="} + enum2str(spatial_method));
  }
}

std::vector<std::vector<double>> Transport_solver::build_fsource(
  const std::vector<std::vector<double>> & flux) const
{
  std::vector<std::vector<double>> fsource;
  fsource.resize(xslib.ngroup());
  for (auto & f : fsource)
    f.resize(geo.dx.size());
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    if (xsthis.isfis)
    {
      double xsum{0.0};
      for (int g = 0; g < xslib.ngroup(); ++g)
        xsum += xsthis.nusf[g] * flux[g][i];
      xsum *= geo.dx[i];
      for (int g = 0; g < xslib.ngroup(); ++g)
        fsource[g][i] = xsthis.chi[g] * xsum;
    }
  }
  return fsource;
}

// TODO this does not support anisotropic scattering
std::vector<std::vector<double>> Transport_solver::build_upscatter(
  const std::vector<std::vector<double>> & flux) const
{
  std::vector<std::vector<double>> upscatter;
  upscatter.resize(xslib.ngroup());
  for (auto & u : upscatter)
    u.resize(geo.dx.size());
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    for (int g = 0; g < xslib.ngroup(); ++g)
    {
      for (int gprime = g+1; gprime < xslib.ngroup(); ++gprime)
        upscatter[g][i] += xsthis.scatter[0](gprime,g) * flux[gprime][i];
      upscatter[g][i] *= geo.dx[i];
    }
  }
  return upscatter;
}

// TODO this does not support anisotropic scattering
std::vector<double> Transport_solver::build_downscatter(
  const std::vector<std::vector<double>> & flux,
  const int g) const
{
  std::vector<double> downscatter;
  downscatter.resize(geo.dx.size());
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis {xslib(geo.mat_map[i])};
    for (int gprime = 0; gprime < g; ++gprime)
      downscatter[i] += xsthis.scatter[0](gprime,g) * flux[gprime][i];
    downscatter[i] *= geo.dx[i];
  }
  return downscatter;
}

double Transport_solver::fission_summation(
  const std::vector<std::vector<double>> & flux) const
{
  double fsum{0.0};
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    if (xsthis.isfis)
    {
      for (int g = 0; g < xslib.ngroup(); ++g)
        fsum += xsthis.nusf[g] * flux[g][i] * geo.dx[i];
    }
  }
  return fsum;
}

double Transport_solver::convergence_phi_scat(
  const std::vector<double> & fluxg,
  const std::vector<double> & fluxg_old)
{
  double xdif{0.0};
  double xmax{0.0};
  for (std::size_t i = 0; i < fluxg.size(); ++i)
  {
    xdif = std::max(xdif, std::abs(fluxg[i] - fluxg_old[i]));
    xmax = std::max(xmax, fluxg[i]);
  }
  return xdif / xmax;
}

double Transport_solver::convergence_phi(
  const std::vector<std::vector<double>> & flux,
  const std::vector<std::vector<double>> & flux_old)
{
  double xdif{0.0};
  double xmax{0.0};
  for (std::size_t g = 0; g < flux.size(); ++g)
  {
    for (std::size_t i = 0; i < flux[g].size(); ++i)
    {
      xdif = std::max(xdif, std::abs(flux[g][i] - flux_old[g][i]));
      xmax = std::max(xmax, flux[g][i]);
    }
  }
  return xdif / xmax;
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
  double fsum{1.0};

  naiad::out << "=== TRANSPORT POWER ITERATION ===" << std::endl;

  for (int iter = 0; iter < tol.max_iter_phi; ++iter)
  {
    const auto flux_old{flux};
    const double k_old{keff};
    const double fsum_old{fsum};

    const auto fsource{build_fsource(flux)};
    const auto upscatter{build_upscatter(flux)};

    for (int g = 0; g < xslib.ngroup(); ++g)
    {
      const auto downscatter{build_downscatter(flux, g)};

      std::vector<double> qmost;
      qmost.resize(geo.dx.size());
      for (std::size_t i = 0; i < geo.dx.size(); ++i)
        qmost[i] = fsource[g][i]/keff + upscatter[g][i] + downscatter[i];

      for (int inner = 0; inner < tol.max_iter_scatter; ++inner)
      {
        const auto fluxg_old{flux[g]};
        flux[g] = sweeper->sweep(flux[g], qmost, g);
        const double dphi{convergence_phi_scat(flux[g], fluxg_old)};
        naiad::out << "   ... scatter " << inner << " dphi=" << std::format("{:7.1e}", dphi) << std::endl;
        if (dphi < tol.scatter)
          break;
      }
    }

    fsum = fission_summation(flux);
    if (iter > 0)
      keff *= fsum / fsum_old;

    const double delta_k{std::abs(keff - k_old)};
    const double delta_phi{convergence_phi(flux, flux_old)};

    naiad::out << "it=" << std::format("{:4d}", iter)
               << " dk=" << std::format("{:7.1e}", delta_k)
               << " dphi=" << std::format("{:7.1e}", delta_phi)
               << " keff=" << std::format("{:8.6f}", keff)
               << std::endl;

    if ((delta_k < tol.k) && (delta_phi < tol.phi))
    {
      naiad::out << "CONVERGENCE!" << std::endl;
      break;
    }

    if (iter == (tol.max_iter_phi-1))
      exception.warning("failed to converge");
  }

  naiad::out << std::endl;

  return {flux, keff};
}

Transport_sweeper::Transport_sweeper(const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
    const XSLibrary & xslib_, const Quadrature1d * const quad_)
  : geo{geo_}, bc_left{bc_left_}, bc_right{bc_right_}, xslib{xslib_}, quad{quad_}
{
  psi_left.resize(xslib.ngroup());
  for (auto & psi : psi_left)
    psi.resize(quad->get_npoints());
  psi_right.resize(xslib.ngroup());
  for (auto & psi : psi_right)
    psi.resize(quad->get_npoints());
}

std::vector<double> Diamond_difference_sweeper::sweep(const std::vector<double> & fluxg_in, const std::vector<double> & qmost, const int g)
{
  const std::vector<double> fluxg_old{fluxg_in};
  std::vector<double> fluxg;
  fluxg.resize(geo.dx.size());

  for (std::size_t j = 0; j < quad->get_npoints(); ++j)
  {
    // TODO negative flux fixup
    const auto qp{quad->get_points()[j]};
    if (qp.x > 0.0)
    {
      double psi_edge;
      switch (bc_left)
      {
        case (Boundary_condition::vacuum):
        {
          psi_edge = 0.0;
          break;
        }
        case (Boundary_condition::mirror):
        {
          const std::size_t jmirror{quad->get_npoints()-1ul-j};
          psi_edge = psi_left[g][jmirror];
          break;
        }
        default:
          exception.fatal(std::string{"Unknown bc_left= "} + enum2str(bc_left));
      }
      for (std::size_t i = 0; i < geo.dx.size(); ++i)
      {
        const auto & xsthis{xslib(geo.mat_map[i])};
        const double q{0.5*(qmost[i] + fluxg_old[i]*xsthis.scatter[0](g,g)*geo.dx[i])};
        const double psi_center{psi_edge / (1.0 + 0.5*xsthis.sigma_t[g]*geo.dx[i]/qp.x) 
          + q/(xsthis.sigma_t[g]*geo.dx[i] + 2.0*qp.x)};
        fluxg[i] += qp.w * psi_center;
        psi_edge = 2.0 * psi_center - psi_edge;
        if (i == geo.dx.size()-1ul)
          psi_right[g][j] = psi_edge;
      }
    }
    else
    {
      double psi_edge;
      switch (bc_right)
      {
        case (Boundary_condition::vacuum):
        {
          psi_edge = 0.0;
          break;
        }
        case (Boundary_condition::mirror):
        {
          const std::size_t jmirror{quad->get_npoints()-1ul-j};
          psi_edge = psi_right[g][jmirror];
          break;
        }
        default:
          exception.fatal(std::string{"Unknown bc_right= "} + enum2str(bc_right));
      }
      for (long i = static_cast<long>(geo.dx.size())-1l; i >= 0; --i)
      {
        const auto & xsthis{xslib(geo.mat_map[i])};
        const double q{0.5*(qmost[i] + fluxg_old[i]*xsthis.scatter[0](g,g)*geo.dx[i])};
        const double psi_center{psi_edge / (1.0 - 0.5*xsthis.sigma_t[g]*geo.dx[i]/qp.x)
          + q/(xsthis.sigma_t[g]*geo.dx[i] - 2.0*qp.x)};
        fluxg[i] += qp.w * psi_center;
        psi_edge = 2.0 * psi_center - psi_edge;
        // make sure to store the left boundary for the opposite directions
        if (i == 0)
          psi_left[g][j] = psi_edge;
      }
    }
  }
  return fluxg;
}

} // namespace naiad
