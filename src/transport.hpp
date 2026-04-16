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
    const Boundary_condition bc_left;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Tolerance & tol;
    const Quadrature1d * const quad;

    std::vector<std::vector<double>> build_fsource(const std::vector<std::vector<double>> & flux) const;
    std::vector<std::vector<double>> build_upscatter(const std::vector<std::vector<double>> & flux) const;
    std::vector<double> build_downscatter(const std::vector<std::vector<double>> & flux, const int g) const;

    double fission_summation(const std::vector<std::vector<double>> & flux) const;

    static double convergence_phi_scat(const std::vector<double> & flux, const std::vector<double> & flux_old);
    static double convergence_phi(const std::vector<std::vector<double>> & flux, const std::vector<std::vector<double>> & flux_old);
};

} // namespace naiad

#endif
