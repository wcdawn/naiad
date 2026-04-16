#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "input.hpp"
#include "geometry.hpp"
#include "xslibrary.hpp"
#include "result.hpp"
#include "quadrature1d.hpp"

#include <vector>

namespace naiad
{

class Transport_solver
{
  public:
    Transport_solver(const Geometry & geo_, const Boundary_condition & bc_left, const Boundary_condition & bc_right,
      const XSLibrary & xslib_, const Tolerance & tol_, const Quadrature1d * const quad_);
    Result solve() const;
  private:
    const Geometry & geo;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Tolerance & tol;
    const Quadrature1d * const quad;
};

} // namespace naiad

#endif
