#include "diffusion.hpp"

#include "exception_handler.hpp"
#include "linalg.hpp"

namespace naiad
{

Diffusion_solver::Diffusion_solver(
  const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
  const XSLibrary & xslib_, const Tolerance & tol_)
  : geo{geo_}, bc_right{bc_right_}, xslib{xslib_}, tol{tol_}
{
  if (bc_left_ != Boundary_condition::mirror)
  {
    exception.fatal(
      std::string{"Left boundary condition must be mirror for diffusion solver."}
      + " specified=" + enum2str(bc_left_));
  }
}

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
  {
    const auto & xsthis{xslib(geo.mat_map[0])};
    const auto & xsnext{xslib(geo.mat_map[1])};
    for (int g = 0; g < xslib.ngroup(); ++g)
    {
      const double dnext{2.0
        * (xsthis.diffusion[g]/geo.dx[0] * xsnext.diffusion[g]/geo.dx[1])
        / (xsthis.diffusion[g]/geo.dx[0] + xsnext.diffusion[g]/geo.dx[1])};
      A[g].dia[0] = dnext + 
        (xsthis.sigma_t[g] - xsthis.scatter[0](g,g)) * geo.dx[0];
      A[g].sup[0] = -dnext;
    }
  }

  for (int g = 0; g < xslib.ngroup(); ++g)
  {
    for (std::size_t i = 1; i < nx-1; ++i)
    {
      const auto & xsprev{xslib(geo.mat_map[i-1])};
      const auto & xsthis{xslib(geo.mat_map[i])};
      const auto xsnext{xslib(geo.mat_map[i+1])};

      const double dprev{2.0
        * (xsthis.diffusion[g] / geo.dx[i] * xsprev.diffusion[g] / geo.dx[i-1])
        / (xsthis.diffusion[g] / geo.dx[i] + xsprev.diffusion[g] / geo.dx[i-1])};
      const double dnext{2.0
        * (xsthis.diffusion[g] / geo.dx[i] * xsnext.diffusion[g] / geo.dx[i+1])
        / (xsthis.diffusion[g] / geo.dx[i] + xsnext.diffusion[g] / geo.dx[i+1])};

      A[g].sub[i-1] = -dprev;
      A[g].dia[i] = dprev + dnext
        + (xsthis.sigma_t[g] - xsthis.scatter[0](g,g)) * geo.dx[i];
      A[g].sup[i] = -dnext;
    }
  }

  // BC at x=L, i=nx-1
  {
    const auto & xsprev{xslib(geo.mat_map[nx-2])};
    const auto & xsthis{xslib(geo.mat_map[nx-1])};
    switch (bc_right)
    {
      case (Boundary_condition::mirror):
      {
        for (int g = 0; g < xslib.ngroup(); ++g)
        {
          const double dprev{2.0
            * (xsthis.diffusion[g] / geo.dx[nx-1] * xsprev.diffusion[g] / geo.dx[nx-2])
            / (xsthis.diffusion[g] / geo.dx[nx-1] + xsprev.diffusion[g] / geo.dx[nx-2])};
          A[g].sub[nx-2] = -dprev;
          A[g].dia[nx-1] = dprev
            + (xsthis.sigma_t[g] - xsthis.scatter[0](g,g)) * geo.dx[nx-1];
        }
        break;
      }
      case (Boundary_condition::zero):
      {
        for (int g = 0; g < xslib.ngroup(); ++g)
        {
          const double dprev{2.0
            * (xsthis.diffusion[g] / geo.dx[nx-1] * xsprev.diffusion[g] / geo.dx[nx-2])
            / (xsthis.diffusion[g] / geo.dx[nx-1] + xsprev.diffusion[g] / geo.dx[nx-2])};
          A[g].sub[nx-2] = -dprev;
          A[g].dia[nx-1] = dprev
            + (xsthis.sigma_t[g] - xsthis.scatter[0](g,g)) * geo.dx[nx-1];
          A[g].dia[nx-1] += 2.0 * xsthis.diffusion[g] / geo.dx[nx-1];
        }
        break;
      }
      default:
        exception.fatal(
          std::string{"Unsupported right boundary condition in diffusion solver."}
          + " specified=" + enum2str(bc_right));
    }
  }

  return A;
}

std::vector<std::vector<double>> Diffusion_solver::build_fsource(
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
};

std::vector<std::vector<double>> Diffusion_solver::build_upscatter(
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
      for (int gprime = g+1; gprime < xslib.ngroup(); ++gprime)
        upscatter[g][i] = xsthis.scatter[0](gprime,g) * flux[gprime][i]; // TODO transpose?
  }

  return upscatter;
}

std::vector<double> Diffusion_solver::build_downscatter(const std::vector<std::vector<double>> & flux, const int g) const
{
  std::vector<double> downscatter;
  downscatter.resize(geo.dx.size());
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    for (int gprime = 0; gprime < g; ++gprime)
      downscatter[i] += xsthis.scatter[0](gprime,g) * flux[i][gprime]; // TODO transpose?
  }
  return downscatter;
}

double Diffusion_solver::fission_summation(const std::vector<std::vector<double>> & flux) const
{
  double xsum{0.0};
  for (std::size_t i = 0; i < geo.dx.size(); ++i)
  {
    const auto & xsthis{xslib(geo.mat_map[i])};
    if (xsthis.isfis)
    {
      for (int g = 0; g < xslib.ngroup(); ++g)
        xsum += xsthis.nusf[g] * flux[g][i];
    }
  }
  return xsum;
}

double Diffusion_solver::convergence_phi(
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
  return xdif/xmax;
}

Result Diffusion_solver::solve() const
{
  const std::vector<Tridiagonal_matrix> trimat{build_matrix()};

  constexpr bool matrix_dump{false};
  if (matrix_dump)
  {
    std::ofstream ofs{"Amat.dat"};
    ofs << std::format("{:.16e} , {:.16e} , {:.16e}\n", 0.0, trimat[0].dia[0], trimat[0].sup[0]);
    for (std::size_t i = 1; i < trimat[0].dia.size()-1; ++i)
      ofs << std::format("{:.16e} , {:.16e} , {:.16e}\n",
          trimat[0].sub[i-1], trimat[0].dia[i], trimat[0].sup[i]);
    ofs << std::format("{:.16e} , {:.16e} , {:.16e}\n", trimat[0].sub.back(), trimat[0].dia.back(), 0.0);
  }

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

  naiad::out << "=== DIFFUSION POWER ITERATION ===" << std::endl;

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

      std::vector<double> q;
      q.resize(geo.dx.size());
      for (std::size_t i = 0; i < geo.dx.size(); ++i)
        q[i] = fsource[g][i]/keff + upscatter[g][i]  + downscatter[g];

      flux[g] = trid(trimat[g].sub, trimat[g].dia, trimat[g].sup, q);
    }

    fsum = fission_summation(flux);
    if (iter > 0)
      keff *= fsum / fsum_old;

    const double delta_k{std::abs(keff - k_old)};
    const double delta_phi{convergence_phi(flux, flux_old)};

    naiad::out << "it=" << std::format("{:4d}", iter)
               << " dk=" << std::format("{:8.1e}", delta_k)
               << " dphi=" << std::format("{:8.1e}", delta_phi)
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

} // namespace naiad
