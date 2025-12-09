#include "analysis.hpp"

#include "exception_handler.hpp"

namespace naiad
{

std::string enum2str(const Analysis_reference & analysis_ref)
{
  switch (analysis_ref)
  {
    case (Analysis_reference::none):
      return "none";
    case (Analysis_reference::onegroup):
      return "onegroup";
    case (Analysis_reference::twogroup):
      return "twogroup";
    case (Analysis_reference::tworegion):
      return "tworegion";
    case (Analysis_reference::critical):
      return "critical";
    default:
      exception.fatal(std::format("Unknown Analysis_reference: {:d}", static_cast<int>(analysis_ref)));
  }
  return "unknown";
}

Analysis_reference str2enum_analysis_reference(const std::string & str)
{
  if (str == "none")
    return Analysis_reference::none;
  if (str == "onegroup")
    return Analysis_reference::onegroup;
  if (str == "twogroup")
    return Analysis_reference::twogroup;
  if (str == "tworegion")
    return Analysis_reference::tworegion;
  if (str == "critical")
    return Analysis_reference::critical;
  exception.fatal(std::string{"Unknown Analysis_reference in conversion : "} + str);
  return Analysis_reference::onegroup;
}

std::vector<double> Analysis::error_linf() const
{
  constexpr bool dump{true};
  const std::vector<std::vector<double>> flux_exact{exact_flux()};

  // First, extrapolate to estimate phi(0.0) in the first group.
  // Then, normalize the flux for comparison such that phi(0.0) = 1.0 for the first group.
  std::vector<std::vector<double>> flux_analysis{res.phi};
  const double phi_zero{-(flux_analysis[0][1] - flux_analysis[0][0])
    * geo.dx[0] / (geo.dx[0] + geo.dx[1]) + flux_analysis[0][0]};
  for (auto & f : flux_analysis)
    for (double & x : f)
      x /= phi_zero;

  const std::vector<double> xcenter{geo.xcenter()};

  std::vector<double> linf_group;
  linf_group.resize(flux_exact.size());
  for (std::size_t g = 0; g < flux_exact.size(); ++g)
  {
    std::unique_ptr<std::ofstream> ofs;
    if (dump)
      ofs = std::make_unique<std::ofstream>(std::format("analytic_g{:d}.csv", g));
    if (ofs)
      *ofs << "x [cm] , exact , calculated\n";
    for (std::size_t i = 0; i < flux_exact[g].size(); ++i)
    {
      if (ofs)
        *ofs << std::format("{:.16e} , {:.16e} , {:.16e}\n", xcenter[i], flux_exact[g][i], flux_analysis[g][i]);
      linf_group[g] = std::max(linf_group[g], std::abs(flux_exact[g][i] - flux_analysis[g][i]));
    }
  }
  return linf_group;
}

double Analysis::error_keff() const
{
  return std::abs(exact_keff() - res.keff); // absolute units
}

void Analysis::summary(std::ostream & os) const
{
  os << "=== ANALYSIS SUMMARY ===" << std::endl;

  os << "keff_exact = " << std::format("{:.16f}", exact_keff()) << std::endl;
  os << "keff_calc  = " << std::format("{:.16f}", res.keff) << std::endl;
  os << "keff_diff  = " << std::format("{:.4f}", error_keff()*1e5) << " [pcm]" << std::endl;

  const std::vector<double> linf{error_linf()};

  for (std::size_t g = 0; g < linf.size(); ++g)
    os << std::format("linf_g{:d}", g) << " = " << std::format("{:.4e}", linf[g]) << std::endl;
  os << std::endl;
}

std::vector<std::vector<double>> Analysis_critical::exact_flux() const
{
  return {}; // return zero-length (no groups)
}

double Analysis_critical::exact_keff() const
{
  return 1.0;
}

std::vector<std::vector<double>> Analysis_onegroup::exact_flux() const
{
  const std::vector<double> xcenter{geo.xcenter()};
  const double Lx{2.0 * geo.xright().back()};
  constexpr double phi0{1.0};

  std::vector<std::vector<double>> flux;
  flux.resize(1);
  for (auto & f : flux)
    f.resize(xcenter.size());

  for (std::size_t i = 0; i < xcenter.size(); ++i)
    flux[0][i] = phi0 * std::cos(M_PI * xcenter[i] / Lx);

  return flux;
}

double Analysis_onegroup::exact_keff() const 
{
  const double Lx{2.0 * geo.xright().back()};
  const double Bsq{std::pow(M_PI / Lx, 2)};
  const auto & xs{xslib("FUEL")};
  const double rem{xs.sigma_t[0] - xs.scatter[0](0,0)};
  return xs.nusf[0] / (xs.diffusion[0] * Bsq + rem);
}

std::vector<std::vector<double>> Analysis_twogroup::exact_flux() const
{
  const double Lx{2.0 * geo.xright().back()};
  const double Bsq{std::pow(M_PI / Lx, 2)};

  std::vector<std::vector<double>> flux;
  flux.resize(2);
  for (auto & f : flux)
    f.resize(geo.dx.size());

  constexpr double phi0{1.0};
  const auto & xs{xslib("FUEL")};
  const double rem2{xs.sigma_t[1] - xs.scatter[0](1,1)};
  const double ratio{xs.scatter[0](0,1) / (xs.diffusion[1] * Bsq + rem2)};
  const std::array<double,2> factor{phi0, phi0*ratio};

  const std::vector<double> xcenter{geo.xcenter()};
  for (std::size_t g = 0; g < flux.size(); ++g)
    for (std::size_t i = 0; i < flux[g].size(); ++i)
      flux[g][i] = factor[g] * std::cos(M_PI * xcenter[i] / Lx);
  return flux;
}

double Analysis_twogroup::exact_keff() const
{
  const double Lx{2.0 * geo.xright().back()};
  const double Bsq{std::pow(M_PI / Lx, 2)};
  const auto & xs{xslib("FUEL")};
  const double rem1{xs.sigma_t[0] - xs.scatter[0](0,0)};
  const double rem2{xs.sigma_t[1] - xs.scatter[0](1,1)};
  const double ratio{xs.scatter[0](0,1) / (xs.diffusion[1] * Bsq + rem2)};
  return (xs.nusf[0] + xs.nusf[1]*ratio) / (xs.diffusion[0] * Bsq + rem1);
}

std::vector<std::vector<double>> Analysis_tworegion::exact_flux() const
{
  constexpr double BF{0.01716859736590981012};
  constexpr double LF{80.0};
  constexpr double LR{100.0};
  const auto & xs_refl{xslib("REFL")};
  const double xs_refl_rem{xs_refl.sigma_t[0] - xs_refl.scatter[0](0,0)};
  const auto KR{std::sqrt(xs_refl_rem / xs_refl.diffusion[0])};
  constexpr double phi0{1.0};
  const std::vector<double> xcenter{geo.xcenter()};
  std::vector<std::vector<double>> flux;
  flux.resize(1);
  flux[0].resize(xcenter.size());
  for (std::size_t i = 0; i < xcenter.size(); ++i)
  {
    if (xcenter[i] <= LF)
      flux[0][i] = phi0 * std::cos(BF * xcenter[i]);
    else
      flux[0][i] = phi0 * std::cos(BF * LF)
        * (std::cosh(KR * (xcenter[i] - LF)) - std::sinh(KR * (xcenter[i] - LF)) / std::tanh(KR*(LR - LF)));
  }
  return flux;
}

double Analysis_tworegion::exact_keff() const
{
  constexpr double BF{0.01716859736590981012};
  const auto & xs_fuel{xslib("FUEL")};
  const double xs_fuel_rem{xs_fuel.sigma_t[0] - xs_fuel.scatter[0](0,0)};
  return xs_fuel.nusf[0] / (xs_fuel.diffusion[0] * std::pow(BF, 2) + xs_fuel_rem);
}

} // namespace naiad
