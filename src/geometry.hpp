#ifndef NAIAD_GEOMETRY_HPP
#define NAIAD_GEOMETRY_HPP

#include <vector>

namespace naiad
{

class Geometry
{
  public:
    Geometry(const std::vector<double> & dx_, const std::vector<int> & mat_map_);

    Geometry()
      : dx{}, mat_map{}
    {}

    void refine();

    std::vector<double> xleft(double xinit = 0.0) const;
    std::vector<double> xright(double xinit = 0.0) const;
    std::vector<double> xcenter(double xinit = 0.0) const;

  private:
    std::vector<double> dx;
    std::vector<int> mat_map;
};

} // namespace naiad

#endif
