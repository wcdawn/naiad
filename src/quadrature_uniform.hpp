#ifndef QUADRATURE_UNIFORM
#define QUADRATURE_UNIFORM

#include <functional>
#include <vector>

#include "quadrature1d.hpp"

namespace naiad
{

class Quadrature_uniform : public Quadrature1d
{
  public:

    Quadrature_uniform(int order) { populate(order); }
    double integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const override;

  private:

    static const std::vector<std::vector<Quadrature_point>> quad;
    void populate(int order);
};

} // namespace naiad

#endif
