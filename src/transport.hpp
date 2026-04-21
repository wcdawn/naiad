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

class Transport_sweeper // abstract base class
{
  public:
    Transport_sweeper(const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right,
        const XSLibrary & xslib_, const Quadrature1d * const quad_, const int pnorder_);
    virtual std::vector<double> sweep(const std::vector<double> & flux, const std::vector<double> & qmost, int g) = 0; // may update psi_left and psi_right
    virtual ~Transport_sweeper(){}
  protected:
    const Geometry & geo;
    const Boundary_condition bc_left;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Quadrature1d * const quad;
    const int pnorder;

    // [group][quadrature_jidx]
    std::vector<std::vector<double>> psi_left;
    std::vector<std::vector<double>> psi_right;
};

class Diamond_difference_sweeper : public Transport_sweeper
{
  public:
    using Transport_sweeper::Transport_sweeper;
    std::vector<double> sweep(const std::vector<double> & flux, const std::vector<double> & qmost, int g);
};

class Transport_solver
{
  public:
    Transport_solver(const Geometry & geo_, const Spatial_method & spatial_method,
      const Boundary_condition & bc_left_, const Boundary_condition & bc_right,
      const XSLibrary & xslib_, const Tolerance & tol_, const Quadrature1d * const quad_, const int pnorder_);
    Result solve() const;
  private:
    const Geometry & geo;
    const Boundary_condition bc_left;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Tolerance & tol;
    const Quadrature1d * const quad;
    const int pnorder;

    std::unique_ptr<Transport_sweeper> sweeper;

    std::vector<std::vector<double>> build_fsource(const std::vector<std::vector<double>> & flux) const;
    std::vector<std::vector<double>> build_upscatter(const std::vector<std::vector<double>> & flux) const;
    std::vector<double> build_downscatter(const std::vector<std::vector<double>> & flux, const int g) const;

    double fission_summation(const std::vector<std::vector<double>> & flux) const;

    static double convergence_phi_scat(const std::vector<double> & flux, const std::vector<double> & flux_old);
    static double convergence_phi(const std::vector<std::vector<double>> & flux, const std::vector<std::vector<double>> & flux_old);
};

} // namespace naiad

#endif
