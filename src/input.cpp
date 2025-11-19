#include "input.hpp"

#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include <iostream>

namespace naiad
{

std::string get_card(std::ifstream & ifs)
{
  std::string card;
  ifs >> card;
  return card;
}

std::string slurp(std::ifstream & ifs)
{
  std::ostringstream sstr;
  sstr << ifs.rdbuf();
  return sstr.str();
}

std::string slurp(const std::string & filename)
{
  std::ifstream ifs{filename};
  return slurp(ifs);
}

Input::Input(const std::string & filename_)
  :filename{filename_}, echo_str{slurp(filename_)}
{

  std::ifstream ifs{filename};

  std::string fname_xslib;

  std::vector<double> dx;
  std::vector<int> mat_map;
  int refine;

  while (ifs.good())
  {

    const std::string card{get_card(ifs)};
    if (card[0] == '#')
    {
      for (char x{'a'}; x != '\n'; ifs.get(x)) {}
      continue;
    }

    if (card == "")
      break;

    if (card == "xslib")
    {
      ifs >> fname_xslib;
    }

    else if (card == "nx")
    {
      std::size_t nx;
      ifs >> nx;
      dx.resize(nx);
      mat_map.resize(nx);
    }
    else if (card == "dx")
    {
      for (double & x : dx)
        ifs >> x;
    }
    else if (card == "mat_map")
    {
      for (int & x : mat_map)
        ifs >> x;
    }
    else if (card == "refine")
    {
      ifs >> refine;
    }
    
    else if (card == "snorder")
    {
      ifs >> snorder;
    }
    else if (card == "pnorder")
    {
      ifs >> pnorder;
    }

    else if (card == "spatial_method")
    {
      const std::string spat{get_card(ifs)};
      spatial_method = str2enum_spatial_method(spat);
    }

    else if (card == "bc_left")
    {
      const std::string bc{get_card(ifs)};
      bc_left = str2enum_boundary_condition(bc);
    }
    else if (card == "bc_right")
    {
      const std::string bc{get_card(ifs)};
      bc_right = str2enum_boundary_condition(bc);
    }
    else
    {
      std::cout << "UNKNOWN INPUT CARD: " << card << std::endl;
    }

  }

  geo = Geometry{dx, mat_map};
  xs = XSLibrary{fname_xslib};

}

void Input::check() const
{
  if (snorder % 2 != 0)
  {
    std::cerr << "ERROR: SN order must be even! snorder=" << snorder << std::endl;
    std::abort();
  }
}

void Input::echo(std::ostream & os) const
{
  os << "=== INPUT ECHO (BEGIN) ===" << std::endl;
  os << echo_str << std::endl;
  os << "=== INPUT ECHO (END) ===" << std::endl;
  os << std::endl;
}

void Input::summary(std::ostream & os) const
{
  os << "=== INPUT SUMMARY ===" << std::endl;
  os << "input filename: " << filename << std::endl;
  os << "xslib filename: " << xs.filename() << std::endl;
  os << "snorder= " << snorder << std::endl;
  os << "pnorder= " << pnorder << std::endl;
  os << "spatial method: " << enum2str(spatial_method) << std::endl;
  os << "BC left: " << enum2str(bc_left) << std::endl;
  os << "BC right: " << enum2str(bc_right) << std::endl;
  os << std::endl;
}

Spatial_method str2enum_spatial_method(const std::string & str)
{
  if (str == "diamond_difference")
    return Spatial_method::diamond_difference;
  if (str == "linear_discontinuous")
    return Spatial_method::linear_discontinuous;
  if (str == "step_characteristic")
    return Spatial_method::step_characteristic;
  if (str == "linear_characteristic")
    return Spatial_method::linear_characteristic;
  if (str == "quadratic_characteristic")
    return Spatial_method::quadratic_characteristic;
  // TODO exception handling
  std::cerr << "Failure to identify spatial method: " << str << std::endl;
  std::abort();
  return Spatial_method::diamond_difference;
}

std::string enum2str(const Spatial_method spatial_method)
{
  switch (spatial_method)
  {
    case (Spatial_method::diamond_difference):
      return "diamond_difference";
    case (Spatial_method::linear_discontinuous):
      return "linear_discontinuous";
    case (Spatial_method::step_characteristic):
      return "step_characteristic";
    case (Spatial_method::linear_characteristic):
      return "linear_characteristic";
    case (Spatial_method::quadratic_characteristic):
      return "quadratic_characteristic";
    default:
      // TODO exception handling
      std::cerr << "Unable to determine spatial_method string: " << static_cast<int>(spatial_method) << std::endl;
      std::abort();
      return "unknown";
  }
}

Boundary_condition str2enum_boundary_condition(const std::string & str)
{
  if (str == "vacuum")
    return Boundary_condition::vacuum;
  if (str == "mirror")
    return Boundary_condition::mirror;
  if (str == "zero")
    return Boundary_condition::zero;
  std::cerr << "Unable to identify boundary condition: " << str << std::endl;
  std::abort();
  return Boundary_condition::vacuum;
}

std::string enum2str(const Boundary_condition bc)
{
  switch (bc)
  {
    case (Boundary_condition::vacuum):
      return "vacuum";
    case (Boundary_condition::mirror):
      return "mirror";
    case (Boundary_condition::zero):
      return "zero";
    default:
      std::cerr << "Unable to identify boundary condition name: " << static_cast<int>(bc) << std::endl;
      std::abort();
      return "unknown";
  }
}

} // namespace naiad
