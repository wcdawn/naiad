#include "transport.hpp"

#include <omp.h>

#include <cmath>

#include "exception_handler.hpp"

namespace naiad
{

Transport_solver::Transport_solver(const Geometry & geo_, const Spatial_method & spatial_method,
                                   const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
                                   const XSLibrary & xslib_, const Tolerance & tol_, const Quadrature1d * const quad_,
                                   const int pnorder_)
  : geo{geo_}, bc_left{bc_left_}, bc_right{bc_right_}, xslib{xslib_}, tol{tol_}, quad{quad_}, pnorder{pnorder_}
{
  switch (spatial_method)
  {
    case (Spatial_method::diamond_difference):
      sweeper = std::make_unique<Diamond_difference_sweeper>(geo, bc_left, bc_right, xslib, quad, pnorder);
      break;
    default:
      exception.fatal(std::string{"Unimplemented spatial transport method. requested="} + enum2str(spatial_method));
  }
}

std::vector<std::vector<double>> Transport_solver::build_fsource(const std::vector<std::vector<double>> & flux,
                                                                 const double keff) const
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
        xsum += xsthis.nusf[g] * flux[g][i * (pnorder + 1) + 0]; // always scalar flux
      xsum *= geo.dx[i];
      for (int g = 0; g < xslib.ngroup(); ++g)
        fsource[g][i] = xsthis.chi[g] / keff * xsum;
    }
  }
  return fsource;
}

// TODO this needs work for anisotropic scattering
std::vector<std::vector<double>> Transport_solver::build_upscatter(const std::vector<std::vector<double>> & flux) const
{
  std::vector<std::vector<double>> upscatter;
  upscatter.resize(xslib.ngroup());
  for (auto & u : upscatter)
    u.resize(geo.dx.size() * (pnorder + 1));
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    for (int g = 0; g < xslib.ngroup(); ++g)
    {
      // only considering isotropic contribution (for now)
      for (int gprime = g + 1; gprime < xslib.ngroup(); ++gprime)
        upscatter[g][i * (pnorder + 1) + 0] += xsthis.scatter[0](gprime, g) * flux[gprime][i * (pnorder + 1) + 0];
      upscatter[g][i * (pnorder + 1) + 0] *= geo.dx[i];
    }
  }
  return upscatter;
}

// TODO this needs work for anisotropic scattering
std::vector<double> Transport_solver::build_downscatter(const std::vector<std::vector<double>> & flux,
                                                        const int g) const
{
  std::vector<double> downscatter;
  downscatter.resize(geo.dx.size() * (pnorder + 1));
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    // only considering isotropic contribution for now
    for (int gprime = 0; gprime < g; ++gprime)
      downscatter[i * (pnorder + 1) + 0] += xsthis.scatter[0](gprime, g) * flux[gprime][i * (pnorder + 1) + 0];
    downscatter[i * (pnorder + 1) + 0] *= geo.dx[i];
  }
  return downscatter;
}

double Transport_solver::fission_summation(const std::vector<std::vector<double>> & flux) const
{
  double fsum{0.0};
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    if (xsthis.isfis)
    {
      for (int g = 0; g < xslib.ngroup(); ++g)
        fsum += xsthis.nusf[g] * flux[g][i * (pnorder + 1) + 0] * geo.dx[i];
    }
  }
  return fsum;
}

// I'm only going to look for relative change in the scalar flux.
// The flux moments will be allowed to "float."
double Transport_solver::convergence_phi_scat(const std::vector<double> & fluxg,
                                              const std::vector<double> & fluxg_old) const
{
  double xdif{0.0};
  double xmax{0.0};
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const double fg{fluxg[i * (pnorder + 1) + 0]};
    xdif = std::max(xdif, std::abs(fg - fluxg_old[i * (pnorder + 1) + 0]));
    xmax = std::max(xmax, fg);
  }
  return xdif / xmax;
}

// I'm only going to look for relative change in the scalar flux.
// The flux moments will be allowed to "float."
// In the future, maybe only look at convergence of fission reaction rate.
double Transport_solver::convergence_phi(const std::vector<std::vector<double>> & flux,
                                         const std::vector<std::vector<double>> & flux_old) const
{
  double xdif{0.0};
  double xmax{0.0};
  for (std::size_t g = 0; g < flux.size(); ++g)
  {
    for (std::size_t i = 0; i < geo.dx.size(); ++i)
    {
      const double fg{flux[g][i * (pnorder + 1) + 0]};
      xdif = std::max(xdif, std::abs(fg - flux_old[g][i * (pnorder + 1) + 0]));
      xmax = std::max(xmax, fg);
    }
  }
  return xdif / xmax;
}

Result Transport_solver::solve() const
{
  std::vector<std::vector<double>> flux;
  flux.resize(xslib.ngroup());
  for (auto & f : flux)
  {
    f.resize(geo.dx.size() * (pnorder + 1));
    // initialize scalar flux (all groups) to unity
    // all higher moments initialized to zero
    for (std::size_t i = 0; i < geo.dx.size(); ++i)
      f[i * (pnorder + 1) + 0] = 1.0;
  }

  double keff{1.0};
  double fsum{1.0};

  naiad::out << "=== TRANSPORT POWER ITERATION ===" << std::endl;
  naiad::out << std::endl;
  naiad::out << "  Iter.    Max. Scat.    Max. Scat.    Delta    Delta     keff  " << std::endl;
  naiad::out << "             Iter.          dphi       keff      phi            " << std::endl;
  naiad::out << " -------  ------------  ------------  -------  -------  --------" << std::endl;

  for (int iter = 0; iter < tol.max_iter_phi; ++iter)
  {
    const auto flux_old{flux};
    const double k_old{keff};
    const double fsum_old{fsum};

    const auto fsource{build_fsource(flux, keff)};
    const auto upscatter{build_upscatter(flux)};

    int max_scat_iter{0};
    double max_scat_dphi{0.0};

    for (int g = 0; g < xslib.ngroup(); ++g)
    {
      const auto downscatter{build_downscatter(flux, g)};

      // TODO this should be down-sized to std::min(xslib.nmoment()+1, pnorder+1)
      // That way, we would only store the non-zero sources.
      // Right now, we would be storing lots of zero sources for pnorder > nmoment.
      std::vector<double> qmost;
      qmost.resize(geo.dx.size() * (pnorder + 1));
      for (std::size_t i = 0; i < geo.dx.size(); ++i)
        for (int ell = 0; ell < pnorder + 1; ++ell)
          qmost[i * (pnorder + 1) + ell] = (ell == 0 ? fsource[g][i] : 0.0) + upscatter[g][i * (pnorder + 1) + ell]
                                           + downscatter[i * (pnorder + 1) + ell];

      for (int inner = 0; inner < tol.max_iter_scatter; ++inner)
      {
        const auto fluxg_old{flux[g]};
        flux[g] = sweeper->sweep(flux[g], qmost, g);
        const double dphi{convergence_phi_scat(flux[g], fluxg_old)};
        if (dphi < tol.scatter)
        {
          max_scat_iter = std::max(max_scat_iter, inner);
          max_scat_dphi = std::max(max_scat_dphi, dphi);
          break;
        }
        else if (inner == tol.max_iter_scatter - 1)
        {
          max_scat_iter = tol.max_iter_scatter; // can't get bigger than this
          max_scat_dphi = std::max(max_scat_dphi, dphi);
        }
      }
    }

    fsum = fission_summation(flux);
    if (iter > 0)
      keff *= fsum / fsum_old;

    const double delta_k{std::abs(keff - k_old)};
    const double delta_phi{convergence_phi(flux, flux_old)};

    naiad::out << std::format(" {:4d} ", iter) << "       " << std::format("{:5d}", max_scat_iter) << "        "
               << std::format("{:7.1e}", max_scat_dphi) << "     " << std::format("{:7.1e}", delta_k) << "  "
               << std::format("{:7.1e}", delta_phi) << "  " << std::format("{:8.6f}", keff) << std::endl;

    if ((delta_k < tol.k) && (delta_phi < tol.phi))
    {
      naiad::out << std::endl << "CONVERGENCE!" << std::endl;
      break;
    }

    if (iter == (tol.max_iter_phi - 1))
      exception.warning("failed to converge");
  }

  naiad::out << std::endl;

  return {flux, keff, pnorder};
}

Transport_sweeper::Transport_sweeper(const Geometry & geo_, const Boundary_condition & bc_left_,
                                     const Boundary_condition & bc_right_, const XSLibrary & xslib_,
                                     const Quadrature1d * const quad_, const int pnorder_)
  : geo{geo_}, bc_left{bc_left_}, bc_right{bc_right_}, xslib{xslib_}, quad{quad_}, pnorder{pnorder_}
{
  psi_left.resize(xslib.ngroup());
  for (auto & psi : psi_left)
    psi.resize(quad->get_npoints());
  psi_right.resize(xslib.ngroup());
  for (auto & psi : psi_right)
    psi.resize(quad->get_npoints());
}

std::vector<double> Diamond_difference_sweeper::sweep(const std::vector<double> & fluxg_in,
                                                      const std::vector<double> & qmost, const int g)
{
  const std::vector<double> fluxg_old{fluxg_in};

  int nthread;
#pragma omp parallel default(none) shared(nthread)
  {
    nthread = omp_get_num_threads();
  }

  // NOTE: I have stored a copy for each thread.
  // This could be performed with locking instead.
  std::vector<std::vector<double>> parfluxg; // [nthread][nx]
  parfluxg.resize(nthread);
  for (auto & fluxg : parfluxg)
    fluxg.resize(geo.dx.size() * (pnorder + 1));

#pragma omp parallel for default(none) shared(quad, bc_left, bc_right, exception, psi_left, psi_right) \
    shared(fluxg_old, g, qmost) shared(parfluxg) shared(std::cout)
  for (std::size_t j = 0; j < quad->get_npoints(); ++j)
  {
    const int myid{omp_get_thread_num()};
    std::vector<double> & fluxg{parfluxg[myid]};
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
          const std::size_t jmirror{quad->get_npoints() - 1ul - j};
          psi_edge = psi_left[g][jmirror];
          break;
        }
        default:
          exception.fatal(std::string{"Unknown bc_left= "} + enum2str(bc_left));
          psi_edge = 0.0;
      }
      for (std::size_t i = 0; i < geo.dx.size(); ++i)
      {
        const auto & xsthis{xslib(geo.mat_map[i])};
        const double q{
            0.5
            * (qmost[i * (pnorder + 1) + 0] + fluxg_old[i * (pnorder + 1) + 0] * xsthis.scatter[0](g, g) * geo.dx[i])};
        const double psi_center{psi_edge / (1.0 + 0.5 * xsthis.sigma_t[g] * geo.dx[i] / qp.x)
                                + q / (xsthis.sigma_t[g] * geo.dx[i] + 2.0 * qp.x)};
        for (int n = 0; n < pnorder + 1; ++n)
          fluxg[i * (pnorder + 1) + n] += std::pow(qp.x, n) * qp.w * psi_center;
        psi_edge = 2.0 * psi_center - psi_edge;
        if (i == geo.dx.size() - 1ul)
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
          const std::size_t jmirror{quad->get_npoints() - 1ul - j};
          psi_edge = psi_right[g][jmirror];
          break;
        }
        default:
          exception.fatal(std::string{"Unknown bc_right= "} + enum2str(bc_right));
          psi_edge = 0.0;
      }
      for (long i = static_cast<long>(geo.dx.size()) - 1l; i >= 0; --i)
      {
        const auto & xsthis{xslib(geo.mat_map[i])};
        const double q{
            0.5
            * (qmost[i * (pnorder + 1) + 0] + fluxg_old[i * (pnorder + 1) + 0] * xsthis.scatter[0](g, g) * geo.dx[i])};
        const double psi_center{psi_edge / (1.0 - 0.5 * xsthis.sigma_t[g] * geo.dx[i] / qp.x)
                                + q / (xsthis.sigma_t[g] * geo.dx[i] - 2.0 * qp.x)};
        for (int n = 0; n < pnorder + 1; ++n)
          fluxg[i * (pnorder + 1) + n] += std::pow(qp.x, n) * qp.w * psi_center;
        psi_edge = 2.0 * psi_center - psi_edge;
        // make sure to store the left boundary for the opposite directions
        if (i == 0)
          psi_left[g][j] = psi_edge;
      }
    }
  }

  // parallel reduction
  std::vector<double> fluxg;
  fluxg.resize(geo.dx.size() * (pnorder + 1));
  for (int n = 0; n < nthread; ++n)
    for (std::size_t i = 0; i < geo.dx.size(); ++i)
      for (int ell = 0; ell < pnorder + 1; ++ell)
        fluxg[i * (pnorder + 1) + ell] += parfluxg[n][i * (pnorder + 1) + ell];

  return fluxg;
}

} // namespace naiad
