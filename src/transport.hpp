#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include <vector>

#include "geometry.hpp"
#include "input.hpp"
#include "quadrature1d.hpp"
#include "result.hpp"
#include "xslibrary.hpp"

namespace naiad
{

class Transport_sweeper // abstract base class
{
  public:

    Transport_sweeper(const Geometry & geo_, const Boundary_condition & bc_left_, const Boundary_condition & bc_right,
                      const XSLibrary & xslib_, const Quadrature1d * const quad_, const int pnorder_,
                      const int src_order_);
    virtual std::vector<double> sweep(const std::vector<double> & flux, const std::vector<double> & qmost, int g)
        = 0; // may update psi_left and psi_right
    virtual ~Transport_sweeper() {}

  protected:

    const Geometry & geo;
    const Boundary_condition bc_left;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Quadrature1d * const quad;
    const int pnorder;
    const int src_order;

    double get_psi_left(const std::size_t j, const int g) const;
    double get_psi_right(const std::size_t j, const int g) const;

    void set_psi_left(const std::size_t j, const int g, const double psi);
    void set_psi_right(const std::size_t j, const int g, const double psi);

    // TODO can these be private?
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

class Step_characteristic_sweeper : public Transport_sweeper
{
  public:

    using Transport_sweeper::Transport_sweeper;
    std::vector<double> sweep(const std::vector<double> & flux, const std::vector<double> & qmost, int g);
};

class Transport_solver
{
  public:

    Transport_solver(const Geometry & geo_, const Spatial_method & spatial_method, const Boundary_condition & bc_left_,
                     const Boundary_condition & bc_right, const XSLibrary & xslib_, const Tolerance & tol_,
                     const Quadrature1d * const quad_, const int pnorder_);
    Result solve() const;

  private:

    const Geometry & geo;
    const Boundary_condition bc_left;
    const Boundary_condition bc_right;
    const XSLibrary & xslib;
    const Tolerance & tol;
    const Quadrature1d * const quad;
    const int pnorder;

    // maximum order of scattering source
    const int src_order;

    std::unique_ptr<Transport_sweeper> sweeper;

    std::vector<std::vector<double>> build_fsource(const std::vector<std::vector<double>> & flux,
                                                   const double keff) const;
    std::vector<std::vector<double>> build_upscatter(const std::vector<std::vector<double>> & flux) const;
    std::vector<double> build_downscatter(const std::vector<std::vector<double>> & flux, const int g) const;

    double fission_summation(const std::vector<std::vector<double>> & flux) const;

    double convergence_phi_scat(const std::vector<double> & flux, const std::vector<double> & flux_old) const;
    double convergence_phi(const std::vector<std::vector<double>> & flux,
                           const std::vector<std::vector<double>> & flux_old) const;
};

} // namespace naiad

#endif
