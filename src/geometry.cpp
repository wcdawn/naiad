#include "geometry.hpp"
#include "exception_handler.hpp"

#include <iostream>
#include <format>
#include <algorithm>

namespace naiad
{

Geometry::Geometry(const std::vector<double> & dx_, const std::vector<int> & mat_map_)
  : dx{dx_}, mat_map{mat_map_}
{
  if (dx.size() != mat_map.size())
  {
    exception.fatal(std::string{"Inconsistent sizes of dx and mat_map.\n"}
              + "dx.size()= " + std::format("{:d}", dx.size()) + "\n"
              + "mat_map.size()= " + std::format("{:d}", mat_map.size()));
  }
}

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

std::vector<double> Geometry::xleft(double xinit) const
{
  std::vector<double> xl;
  xl.resize(dx.size());
  xl[0] = xinit;
  for (std::size_t i = 1; i < xl.size(); ++i)
    xl[i] = xl[i-1] + dx[i-1];
  return xl;
}

std::vector<double> Geometry::xright(double xinit) const
{
  std::vector<double> xr;
  xr.resize(dx.size());
  xr[0] = xinit + dx[0];
  for (std::size_t i = 1; i < xr.size(); ++i)
    xr[i] = xr[i-1] + dx[i];
  return xr;
}

std::vector<double> Geometry::xcenter(double xinit) const
{
  std::vector<double> xc;
  xc.resize(dx.size());
  xc[0] = xinit + 0.5*dx[0];
  for (std::size_t i = 1; i < xc.size(); ++i)
    xc[i] = xc[i-1] + 0.5*(dx[i-1]+dx[i]);
  return xc;
}

void Geometry::summary(std::ostream & os) const
{
  os << "=== GEOMETRY SUMMARY ===" << std::endl;
  os << "nx= " << dx.size() << std::endl;
  if (dx.size() <= std::size_t{10})
  {
    os << "dx=";
    for (const double & x : dx)
      os << std::format(" {:.2e}", x);
    os << std::endl;
  }
  else
  {
    os << "dx omitted for space" << std::endl;
  }
  os << "maximum dx= " << std::format("{:.2e}", *std::max_element(dx.begin(), dx.end())) << std::endl;
  os << "minimum dx= " << std::format("{:.2e}", *std::min_element(dx.begin(), dx.end())) << std::endl;
  if (mat_map.size() <= std::size_t{10})
  {
    os << "mat_map=";
    for (const auto & x : mat_map)
      os << " " << x;
    os << std::endl;
  }
  else
  {
    os << "mat_map omitted for space" << std::endl;
  }
  os << std::endl;
}

} // namespace naiad
