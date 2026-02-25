#ifndef QUADRATURE1D_HPP
#define QUADRATURE1D_HPP

#include <vector>
#include <iostream>
#include <functional>

namespace naiad
{

class Quadrature_point
{
  public:
    double x;
    double w;
};

class Quadrature1d
{
  public:
    std::vector<Quadrature_point> get_points() const { return points; }
    auto get_npoints() const { return points.size(); }
    virtual double integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const = 0;
    virtual ~Quadrature1d() {}
  protected:
    std::vector<Quadrature_point> points;
};

} // namespace naiad

#endif
