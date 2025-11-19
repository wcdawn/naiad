#include "geometry.hpp"

namespace naiad
{

void Geometry::refine()
{
  const std::vector<double> dx_old{dx};
  const std::vector<int> mat_map_old{mat_map};
  dx.resize(std::size_t{2}*dx_old.size());
  mat_map.resize(std::size_t{2}*mat_map_old.size());
  for (std::size_t i = 0; i < dx_old.size(); ++i)
  {
    dx[2*i+0] = 0.5*dx_old[i];
    dx[2*i+1] = 0.5*dx_old[i];
    mat_map[2*i+0] = mat_map_old[i];
    mat_map[2*i+1] = mat_map_old[i];
  }
}

} // namespace naiad
