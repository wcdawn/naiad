#include "writer.hpp"

#include <fstream>

namespace naiad
{

void Writer::write_flux(const std::string & fname) const
{
  const std::vector<double> xcenter{geo.xcenter()};
  std::ofstream ofs{fname};

  ofs << "x [cm]";
  for (std::size_t g = 0; g < res.phi.size(); ++g)
    ofs << " , flux_g" << std::format("{:d}", g);
  ofs << std::endl;

  for (std::size_t i = 0; i < xcenter.size(); ++i)
  {
    ofs << std::format("{:.16e}", xcenter[i]);
    for (std::size_t g = 0; g < res.phi.size(); ++g)
      ofs << " , " << std::format("{:.16e}", res.phi[g][i]);
    ofs << std::endl;
  }
}

} // namespace naiad
