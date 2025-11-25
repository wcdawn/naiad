#ifndef NAIAD_DIFFUSION_HPP
#define NAIAD_DIFFUSION_HPP

#include "input.hpp"
#include "geometry.hpp"
#include "xslibrary.hpp"
#include "result.hpp"

#include <vector>

namespace naiad
{

class Tridiagonal_matrix
{
  public:
    std::vector<double> sub;
    std::vector<double> dia;
    std::vector<double> sup;
};

class Diffusion_solver
{
  public:
    // TODO need to somehow pass around the boundary condition
    // and need to check if they are reasonable values
    // i.e., no vaccum, zero and mirror, must be mirror at left
    Diffusion_solver(
        const Geometry & geo_, const XSLibrary & xslib_, const Tolerance & tol_)
      : geo{geo_}, xslib{xslib_}, tol{tol_}
    {}

    Result solve() const;

  private:
    const Geometry & geo;
    const XSLibrary & xslib;
    const Tolerance & tol;

    // sub, dia, sup
    std::vector<Tridiagonal_matrix> build_matrix() const;

};

} // namespace naiad

#endif
