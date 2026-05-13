#include "quadrature_uniform.hpp"

namespace naiad
{

double Quadrature_uniform::integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const
{
  double xsum{0.0};
  for (const auto & qp : points)
  {
    const double xi{0.5 * (xhi - xlo) * qp.x + 0.5 * (xhi + xlo)}; // coordinate transform
    xsum += qp.w * f(xi);
  }
  return (xhi - xlo) * 0.5 * xsum;
}

void Quadrature_uniform::populate(int order)
{
  points.reserve(order);
  const double dx{2.0 / static_cast<double>(order)};
  for (int i = 0; i < order; ++i)
    points.emplace_back(Quadrature_point{.x = (-1.0 + dx * (i + 0.5)), .w = dx});
}

} // namespace naiad
