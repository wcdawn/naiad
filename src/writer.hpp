#ifndef NAIAD_WRITER_HPP
#define NAIAD_WRITER_HPP

#include "geometry.hpp"
#include "result.hpp"

#include <string>

namespace naiad
{

class Writer
{
  public:
    Writer(const Geometry & geo_, const Result & res_)
      : geo{geo_}, res{res_}
    {}

    void write_flux(const std::string & fname) const;

  private:
    const Geometry & geo;
    const Result & res;


};

} // namespace naiad

#endif
