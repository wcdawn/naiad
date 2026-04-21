#include "writer.hpp"
#include "output.hpp"

#include <fstream>

namespace naiad
{

void Writer::write_flux(const std::string & fname) const
{
  const std::vector<double> xcenter{geo.xcenter()};
  std::ofstream ofs{fname};

  ofs << "x [cm]";
  for (std::size_t g = 0; g < res.phi.size(); ++g)
    ofs << " , flux_g" << g;
  ofs << std::endl;

  for (std::size_t i = 0; i < xcenter.size(); ++i)
  {
    ofs << std::format("{:.16e}", xcenter[i]);
    for (std::size_t g = 0; g < res.phi.size(); ++g)
      ofs << " , " << std::format("{:.16e}", res.phi[g][i]);
    ofs << std::endl;
  }
}

void Writer::write_phi(const std::string & fname) const
{
  const std::vector<double> xcenter{geo.xcenter()};
  std::ofstream ofs{fname};

  constexpr int pnorder{0};

  ofs << "x [cm]";
  for (int n = 0; n < pnorder+1; ++n)
    for (std::size_t g = 0; g < res.phi.size(); ++g)
      ofs << " , phi_n" << n << "_g" << g;
  ofs << std::endl;

  for (std::size_t i = 0; i < xcenter.size(); ++i)
  {
    ofs << std::format("{:.16e}", xcenter[i]);
    for (int n = 0; n < pnorder+1; ++n)
      for (std::size_t g = 0; g < res.phi.size(); ++g)
        ofs << " , " << std::format("{:.16e}", res.phi[g][i*(pnorder+1) + n]);
    ofs << std::endl;
  }
}

} // namespace naiad
