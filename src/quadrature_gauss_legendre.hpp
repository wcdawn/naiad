#ifndef QUADRATURE_GAUSS_LEGENDRE
#define QUADRATURE_GAUSS_LEGENDRE

#include <vector>
#include <iostream>
#include <functional>

#include "quadrature1d.hpp"

namespace naiad
{

class Quadrature_gauss_legendre : public Quadrature1d
{
  public:
    Quadrature_gauss_legendre(int order)
    {
      populate(order);
    }
    double integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const override;
  private:
    static const std::vector<std::vector<Quadrature_point>> quad;
    void populate(int order);
};

} // namespace naiad

#endif
