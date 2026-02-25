#ifndef QUADRATURE_GAUSS_LEGENDRE
#define QUADRATURE_GAUSS_LEGENDRE

#include <vector>
#include <iostream>
#include <functional>

#include "quadrature1d.hpp"

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
    void populate(int order)
    {
      try
      {
        points = quad.at(order-1);
      }
      catch (...)
      {
        std::cerr << "Unacceptable Gauss-Legendre order: " << std::to_string(order) << std::endl;
        std::abort();
      }
    }
};

#endif
