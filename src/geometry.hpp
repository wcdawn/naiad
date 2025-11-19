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

  private:
    std::vector<double> dx;
    std::vector<int> mat_map;
};

} // namespace naiad

#endif
