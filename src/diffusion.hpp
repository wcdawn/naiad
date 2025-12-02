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
    Diffusion_solver(
        const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right_,
        const XSLibrary & xslib_, const Tolerance & tol_);

    Result solve() const;

  private:
    const Geometry & geo;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Tolerance & tol;

    // sub, dia, sup
    std::vector<Tridiagonal_matrix> build_matrix() const;

    std::vector<std::vector<double>> build_fsource(const std::vector<std::vector<double>> & flux) const;
    std::vector<std::vector<double>> build_upscatter(const std::vector<std::vector<double>> & flux) const;
    std::vector<double> build_downscatter(const std::vector<std::vector<double>> & flux, const int g) const;

    double fission_summation(const std::vector<std::vector<double>> & flux) const;

    static double convergence_phi(const std::vector<std::vector<double>> & flux, const std::vector<std::vector<double>> & flux_old);

};

} // namespace naiad

#endif
