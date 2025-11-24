#ifndef NAIAD_INPUT_HPP
#define NAIAD_INPUT_HPP

#include <string>
#include <fstream>

#include "geometry.hpp"
#include "xslibrary.hpp"

namespace naiad
{

std::string get_card(std::ifstream & ifs);
std::string slurp(std::ifstream & ifs);
std::string slurp(const std::string & filename);

class Tolerance
{
  public:
    double phi{1e-7};
    double k{1e-6};
    double scatter{1e-7};

    int max_iter_phi{1'000};
    int max_iter_scatter{100};
};

enum class Spatial_method
{
  diamond_difference,
  linear_discontinuous,
  step_characteristic,
  linear_characteristic,
  quadratic_characteristic,
};

Spatial_method str2enum_spatial_method(const std::string & str);
std::string enum2str(const Spatial_method spatial_method);

enum class Boundary_condition
{
  vacuum,
  mirror,
  zero,
};

Boundary_condition str2enum_boundary_condition(const std::string & str);
std::string enum2str(const Boundary_condition bc);

class Input
{
  public:
    Input(const std::string & filename_);

    void check() const;
    void echo(std::ostream & os) const;
    void summary(std::ostream & os) const;

    const auto & xslibrary() const { return xs; }
    const auto & geometry() const { return geo; }
    const auto & tolerance() const { return tol; }

    int refine{0};

    int pnorder{0};
    int snorder{0}; // snorder==0 will be used to access the diffusion solver

    Spatial_method spatial_method{Spatial_method::diamond_difference};

    Boundary_condition bc_left{Boundary_condition::mirror};
    Boundary_condition bc_right{Boundary_condition::mirror};

  private:
    const std::string filename;
    const std::string echo_str;
    Geometry geo;
    XSLibrary xs;
    Tolerance tol;

};

} // namespace naiad

#endif
