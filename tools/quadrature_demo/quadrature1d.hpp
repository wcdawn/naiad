#ifndef QUADRATURE1D_HPP
#define QUADRATURE1D_HPP

#include <vector>
#include <iostream>
#include <functional>

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
    virtual double integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const = 0;
    virtual ~Quadrature1d() {}
  protected:
    std::vector<Quadrature_point> points;
};

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
