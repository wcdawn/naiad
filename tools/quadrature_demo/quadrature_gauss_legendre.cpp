#include "quadrature_gauss_legendre.hpp"

const std::vector<std::vector<Quadrature_point>> Quadrature_gauss_legendre::quad = std::vector<std::vector<Quadrature_point>>{
#include "../gauss_legendre.txt"
};

double Quadrature_gauss_legendre::integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const
{
  double xsum{0.0};
  for (const auto & qp : points)
  {
    const double xi{0.5*(xhi-xlo)*qp.x + 0.5*(xhi+xlo)}; // coordinate transform
    xsum += qp.w * f(xi);
  }
  return (xhi-xlo)*0.5*xsum;
}
