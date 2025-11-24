#ifndef NAIAD_DIFFUSION_HPP
#define NAIAD_DIFFUSION_HPP

#include "input.hpp"
#include "geometry.hpp"
#include "xslibrary.hpp"
#include "result.hpp"

namespace naiad
{

class Diffusion_solver
{
  public:
    Diffusion_solver(
        const Geometry & geo_, const XSLibrary & xslib_, const Tolerance & tol_)
      : geo{geo_}, xslib{xslib_}, tol{tol_}
    {}

    Result solve() const;

  private:
    const Geometry & geo;
    const XSLibrary & xslib;
    const Tolerance & tol;

};

} // namespace naiad

#endif
